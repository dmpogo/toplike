MODULE Topology_map_mod
  USE Topology_types
  USE FITSTOOLS
  USE HEAD_FITS
  USE PIX_TOOLS
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: make_fake_map, Write_map, ReadExpData
  
CONTAINS

  SUBROUTINE make_fake_map(ampl,map)
    ! Works on eigenvalue decomposition of CTpp in CTpp_eval and CTpp_evec
    USE beams
    USE ALM_TOOLS
    USE RNGMOD
    real(DP), intent(in)                               :: ampl
    real(DP), intent(out), allocatable, dimension(:,:) :: map

    type(planck_rng)        :: rng_handle
    real(DP),    allocatable, dimension(:,:)   :: random_numbers
    complex(DP), allocatable, dimension(:,:,:) :: alm
    integer(I4B)            :: i, iter_order=5
 
    if (allocated(map)) deallocate(map)
    allocate(map(0:npix_fits-1,1:1))
! if we add noise, generate random full sky noise map, else set it to zero
    call rand_init(rng_handle,iseed)
    if (add_map_noise) then
       do_Gsmooth=.false.
       do_mask=.false.
       call ReadExpData(expdata_format)
       if (npix_cut /= npix_fits) stop 'Must be full sky. Check noise masking'

       do i = 0, npix_fits-1
          map(i,1)=rand_gauss(rng_handle)*sqrt(map_npp(i,i))
       enddo
       write(0,*)'random noise map has been generated'
    else
       map=0.0_dp
    endif

! Generate random set of CTpp eigenvalues (sqrt)
    allocate(random_numbers(0:n_evalues-1,1))
    do i = 0, n_evalues - 1
       random_numbers(i,1) = rand_gauss(rng_handle)*sqrt(ampl*CTpp_eval(i))
    enddo

! Combine correlated and noise map on full sky
    call dgemm('N','N',npix_fits,1,n_evalues,1.0_dp,CTpp_evec,npix_fits,random_numbers,n_evalues,1.0_dp,map,npix_fits)
    deallocate(random_numbers)
    write(0,*)'random signal map has been generated and added to noise'

! Smooth full-sky map if needed
    if (beam_fwhm > 0.0_dp) then
       call collect_beams(Wl,lmax,G_fwhm=beam_fwhm,reset=.true.)
       allocate( alm(1:1,0:lmax,0:lmax) )
       call map2alm_iterative(nside,lmax,lmax,iter_order,map,alm,(/0.0_dp, 0.0_dp/),w8ring)
       call alter_alm(nside,lmax,lmax,beam_fwhm,alm,window=Wl)
       call alm2map(nside,lmax,lmax,alm,map(:,1))
       write(0,*)'random map has been smoothed with Gaussian beam'
       deallocate( alm )
    endif

    return
  END SUBROUTINE make_fake_map

  SUBROUTINE Write_map(heal_map)
    REAL(DP), INTENT(IN), DIMENSION(:,:)   :: heal_map
    !-----------------------------------------------------------------------
    !                        generates header
    !-----------------------------------------------------------------------
    INTEGER, PARAMETER :: nlheader = 80
    CHARACTER(LEN=80), DIMENSION(1:nlheader) :: header
    INTEGER :: iostatus, npix_loc, nside_loc
    REAL    :: nullval
    LOGICAL :: anynull,filefound

    !heal_map(:,1) = unpack(map_signal,map_mask,-1.6375d30)
    npix_loc =size(heal_map,1)
    nside_loc=npix2nside(npix_loc)

    CALL write_minimal_header(header,'MAP',nside=nside_loc,ordering='RING',creator='Topology_make_map',coordsys='G',randseed=iseed,units='mK',nlmax=lmax)
!    CALL add_card(header,'COMMENT','-----------------------------------------------')
    INQUIRE(file=TRIM(ADJUSTL(fake_file)),exist=filefound)
    IF(filefound) THEN
       open(26,file=TRIM(ADJUSTL(fake_file)),form='unformatted')
       close(26,status='delete')
       write(0,*) 'Overwriting existing fits file'
    ELSE
       open(26,file=TRIM(ADJUSTL(fake_file)),form='unformatted',status='new',iostat=iostatus)
       IF ( iostatus > 0 ) THEN
          write(0,*) 'Unable to open output ',TRIM(ADJUSTL(fake_file))
          stop
       ELSE
          close(26,status='delete')
       ENDIF
    ENDIF
    
    write(0,*)'Writing bintab'
    write(0,*)size(heal_map,1),size(heal_map,2),npix_loc,nlheader
    CALL write_bintab(heal_map, npix_loc, 1, header, nlheader,TRIM(ADJUSTL(fake_file)))
    write(0,*)'Done'

    RETURN
  END SUBROUTINE Write_map

  SUBROUTINE ReadExpData(format_choice)
    character*(*), intent(in)    :: format_choice

    if (index(ADJUSTL(format_choice),'WMAP')) then
       call ReadWMAP_map()
    else if (index(ADJUSTL(format_choice),'PLANCK')) then
       call ReadPlanck_map()
    else
       stop 'Unknown input format'
    endif
  END SUBROUTINE ReadExpData

  SUBROUTINE ReadWMAP_map()
    USE ctpplm_transforms
    USE beams
    USE ALM_TOOLS
    !Global map_signal,map_npp,diag_noise,map_signal_file,map_mask_file,nside
    
    INTEGER :: i,j,ordering,lcount,nmaps,iter_order=5
    logical :: convert_from_nested=.false.
    complex(DP), DIMENSION(:,:,:), ALLOCATABLE :: alm
    complex(DP), DIMENSION(:,:),   ALLOCATABLE :: clm
    REAL(DP), DIMENSION(:),   ALLOCATABLE :: diag_noise
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: wmap_noise, wmap_signal
    REAL(SP), DIMENSION(:,:), ALLOCATABLE :: wmap_data, wmap_mask 
    
! Input must be full-sky in globally accessible files
    npix_fits=getsize_fits(map_signal_file,nmaps=nmaps,ordering=ordering)
    if ( nside2npix(nside) /= npix_fits ) then
       stop 'Mismatch between map size and expected nside'
    endif

! Allocate arrays and input necessary data. 
    ALLOCATE(wmap_data(0:npix_fits-1,1:nmaps))
    WRITE(0,'(a,a)') '> ', TRIM(ADJUSTL(map_signal_file))
    CALL input_map(TRIM(ADJUSTL(map_signal_file)),wmap_data,npix_fits,nmaps)

! Check ordering, CTpp is probably given in RING (to do: auto-synchronized)
    if ( ordering == 0 ) then
       write(0,*)'Ordering of the input map is unknown, assumed RING'
    else if ( ordering == 2 ) then
       write(0,*)'Input converted from NESTED to RING'
       call convert_nest2ring(nside,wmap_data)
       convert_from_nested=.true.
    else
       write(0,*)'Input is in RING pixelization'
    endif

    ALLOCATE(wmap_noise(0:npix_fits-1,1:1))
    IF( add_noise .and. nmaps > 1 ) THEN
       ! Assumes diagonal noise, 
       ! and wmap_data to contain inverse noise in the last map
       where( wmap_data(:,nmaps) > 0.0 ) wmap_noise(:,1) = 1.d0/wmap_data(:,nmaps)
    ELSE
       wmap_noise = 0.0_dp
    ENDIF

    ALLOCATE(map_mask(0:npix_fits-1))
    IF(do_mask) THEN
       ALLOCATE(wmap_mask(0:npix_fits-1,1:1))
       WRITE(0,'(a,a)') '> ', TRIM(ADJUSTL(map_mask_file))
       CALL input_map(TRIM(ADJUSTL(map_mask_file)),wmap_mask,npix_fits,1)
       map_mask = ( wmap_mask(:,1) /= 0 ) 
       DEALLOCATE(wmap_mask)
    ELSE
       map_mask = .true.
    ENDIF

! Count unmasked pixels 
    npix_cut = count(map_mask)
    write(0,'(a,i7,a,i7)') 'Found ',npix_cut,' unmasked pixels from ',npix_fits
    if (npix_cut == 0) STOP 'All pixels are masked'

! data is in, now process it and store in global arrays

!Smooth and pack the signal.

    ALLOCATE(map_signal(0:npix_cut-1))
    IF (do_Gsmooth) THEN
       ALLOCATE( alm(1:1,0:lmax,0:lmax) )
       ALLOCATE( wmap_signal(0:npix_fits-1,1:1) )
       wmap_signal(:,1:1) = wmap_data(:,1:1) ! unfortunate copy to go to DP
       call map2alm_iterative(nside,lmax,lmax,iter_order,wmap_signal,alm,(/0.0_dp, 0.0_dp/),w8ring)
       call alter_alm(nside,lmax,lmax,beam_fwhm,alm,window=Wl)
       call alm2map(nside,lmax,lmax,alm,wmap_signal(:,1))
       map_signal = pack(wmap_signal(:,1),map_mask)
       DEALLOCATE(alm, wmap_signal)
    ELSE
       map_signal = pack(wmap_data(:,1),map_mask)
    ENDIF
    DEALLOCATE(wmap_data)

!Smooth noise if needed and store it in a cut sky matrix map_npp
    ALLOCATE(map_npp(0:npix_cut-1,0:npix_cut-1))
    IF (add_noise.and.do_Gsmooth) THEN
       call getclm(clm, lcount, wmap_noise, npix_fits, lmax, w8_file=w8_file)
       call smooth_clm(clm, lcount, Wl(:,1), lmax)
       call getcpp(map_npp, npix_cut, clm, lcount, nside, mask=map_mask)
    ELSE
       map_npp=0.0_dp
    ENDIF

! diagonal noise needs always to be defined, although may be zero
    ALLOCATE(diag_noise(0:npix_cut-1))
    if ( epsil > 0.0 ) then
        diag_noise = epsil
    else 
        diag_noise = 0.0d0
    endif

    ! If noise was smoothed, it is in map_npp, otherwise add sigma^2/hit_counts
    if (add_noise.and.(.not.do_Gsmooth) ) then
       diag_noise = diag_noise + pack(wmap_noise(:,1),map_mask)
    endif
    DEALLOCATE(wmap_noise)

    FORALL(i=0:npix_cut-1) map_npp(i,i) = map_npp(i,i) + diag_noise(i)

    RETURN
  END SUBROUTINE ReadWMAP_map

  SUBROUTINE ReadPlanck_map()
    USE ctpplm_transforms
    USE beams
    USE ALM_TOOLS
    !Global map_signal,map_npp,diag_noise,map_signal_file,map_mask_file,nside
    
    INTEGER :: i,j,ordering,lcount,npix,nmaps,iter_order=5
    logical :: convert_from_nested=.false.
    complex(DP), DIMENSION(:,:,:), ALLOCATABLE :: alm
    complex(DP), DIMENSION(:,:),   ALLOCATABLE :: clm
    REAL(DP), DIMENSION(:),   ALLOCATABLE :: diag_noise, WORK
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: exp_noise,exp_data,exp_mask,bpp
    REAL(DP), DIMENSION(0:3)              :: mpoles
    
! Input must be full-sky in globally accessible files
    npix=getsize_fits(map_signal_file,nmaps=nmaps,ordering=ordering)
    npix_fits=nside2npix(nside)
    if (npix > npix_fits) then
       write(0,*)'Input map will be coarsened from nside=',npix2nside(npix)
    elseif (npix < npix_fits) then
       stop 'Analysis on map coarser than CTpp is not implemented'
    endif

! Allocate arrays and input necessary data. 
    ALLOCATE(exp_data(0:npix-1,1:1))
    WRITE(0,'(a,a)') '> ', TRIM(ADJUSTL(map_signal_file))
    CALL input_map(TRIM(ADJUSTL(map_signal_file)),exp_data,npix,1)
    exp_data=exp_data*expdata_scale   ! convert to mK, from whatever input is


    IF( add_noise) THEN
       ALLOCATE(exp_noise(0:npix-1,1:1))
       ! Assumes diagonal noise, 
       ! and wmap_data to contain noise per pixel in same units as map
       WRITE(0,'(a,a)') '> ', TRIM(ADJUSTL(map_noise_file))
       CALL input_map(TRIM(ADJUSTL(map_noise_file)),exp_noise,npix,1)
       exp_noise=exp_noise*expdata_scale  ! convert to mK
       exp_noise=exp_noise**2             ! make the variance
    ELSE
       ALLOCATE(exp_noise(0:npix_fits-1,1:1))
       exp_noise = 0.0_dp
    ENDIF

    IF(do_mask) THEN
       ALLOCATE(exp_mask(0:npix-1,1:1))
       WRITE(0,'(a,a)') '> ', TRIM(ADJUSTL(map_mask_file))
       CALL input_map(TRIM(ADJUSTL(map_mask_file)),exp_mask,npix,1)
    ELSE
       ALLOCATE(exp_mask(0:npix_fits-1,1:1))
       exp_mask = 1
    ENDIF

! Check ordering, CTpp is probably given in RING (to do: auto-synchronized)
    if ( ordering == 0 ) then
       write(0,*)'Ordering of the input map is unknown, assumed RING'
    else if ( ordering == 2 ) then
       write(0,*)'Input converted from NESTED to RING'
       convert_from_nested=.true.
    else
       write(0,*)'Input is in RING pixelization'
    endif

    if (npix > npix_fits) then 
       CALL REBIN_MAP(exp_data,exp_noise,exp_mask,npix,npix_fits,convert_from_nested)
    endif

    if ( convert_from_nested ) then
       call convert_nest2ring(nside,exp_data)
       call convert_nest2ring(nside,exp_mask)
       call convert_nest2ring(nside,exp_noise)
    endif

    ALLOCATE(map_mask(0:npix_fits-1))
    ! This threshold interplace with one in REBIN, if mask is not scaled to 0-1
    ! It was introduced to work agains smoothing of masked map issues
    ! to mask edge pixels - that's why a bit over 50% is a good criterium
    map_mask = ( exp_mask(:,1) >= 0.6_dp )
    DEALLOCATE(exp_mask)

! Count unmasked pixels 
    npix_cut = count(map_mask)
    write(0,'(a,i7,a,i7)') 'Found ',npix_cut,' unmasked pixels from ',npix_fits
    if (npix_cut == 0) STOP 'All pixels are masked'

! data is in, now process it and store in global arrays

! Smooth and pack the signal.
    if (.false.) then       ! Hack: hard choice of real space smoothing
    if (do_Gsmooth) then
       ALLOCATE ( bpp(0:npix_fits-1,0:npix_fits-1) )
       call create_beam(bpp,Wl(:,1),lmax,nside=nside)
       call weight_beam(bpp,1)

       ! Cut sky weighting (seems like after pixel weightening ?)
       do j=0,npix_fits-1
          bpp(:,j) = bpp(:,j)/sum(bpp(:,j),mask=map_mask)
       enddo
       bpp = bpp*(FOURPI*npix_cut/npix_fits)

       do j=0,npix_fits-1
          bpp(:,j) = pack(bpp(:,j),map_mask)
       enddo
       do i=0,npix_fits-1
          bpp(i,:) = pack(bpp(i,:),map_mask)
       enddo

    endif

    ALLOCATE(map_signal(0:npix_cut-1))
    IF (do_Gsmooth) THEN
       allocate(WORK(0:npix_cut-1))
       WORK = pack(exp_data(:,1),map_mask)
       call DGEMV('T',npix_cut,npix_cut,1.d0,bpp,npix_fits,WORK,1,0.d0,map_signal,1)
    ELSE
       map_signal = pack(exp_data(:,1),map_mask)
    ENDIF

!Smooth noise if needed and store it in a cut sky matrix map_npp
    ALLOCATE(map_npp(0:npix_cut-1,0:npix_cut-1))
    IF (add_noise.and.do_Gsmooth) THEN
       WORK = pack(exp_noise(:,1),map_mask)
       do i=0,npix_cut-1
          bpp(i,:)=sqrt(WORK(i))*bpp(i,:)
       enddo
       call DGEMM('T','N',npix_cut,npix_cut,npix_cut,1.d0,bpp,npix_fits,bpp,npix_fits,0.d0,map_npp,npix_cut)
       deallocate(WORK)
    ELSE
       map_npp=0.0_dp
    ENDIF

    else      ! Hack:   to choose lm smoothing

    ALLOCATE(map_signal(0:npix_cut-1))
    IF (do_Gsmooth) THEN
       ALLOCATE( alm(1:1,0:lmax,0:lmax) )
       call map2alm_iterative(nside,lmax,lmax,iter_order,exp_data,alm,(/0.0_dp, 0.0_dp/),w8ring)
       write(0,*) alm(1:1,0:1,0:1)
       ! alm(1:1,0:1,0:1) = cmplx(0.0_dp,0.0_dp)
       call alter_alm(nside,lmax,lmax,beam_fwhm,alm,window=Wl)
       call alm2map(nside,lmax,lmax,alm,exp_data(:,1))
       map_signal = pack(exp_data(:,1),map_mask)
       DEALLOCATE(alm)
    ELSE
       map_signal = pack(exp_data(:,1),map_mask)
    ENDIF

!Smooth noise if needed and store it in a cut sky matrix map_npp
   ALLOCATE(map_npp(0:npix_cut-1,0:npix_cut-1))
    IF (add_noise.and.do_Gsmooth) THEN
       call getclm(clm, lcount, exp_noise, npix_fits, lmax, w8_file=w8_file)
       call smooth_clm(clm, lcount, Wl(:,1), lmax)
       call getcpp(map_npp, npix_cut, clm, lcount, nside, mask=map_mask)
    ELSE
       map_npp=0.0_dp
    ENDIF

    endif         ! End the Hack

! Subtract monopole and dipole
    allocate(WORK(0:npix_fits-1))
    WORK = unpack(map_signal,map_mask,-1.6375d30)
    call remove_dipole(nside,WORK,1,2,mpoles,(/0._dp,0._dp/))
    map_signal = pack(WORK,map_mask)
    deallocate(WORK)

!    fake_file='ppp_0.1_0.1_rs'
!    exp_data(:,1) = unpack(map_signal,map_mask,-1.6375d30)
!    call Write_map(exp_data)
    DEALLOCATE(exp_data)

! diagonal noise needs always to be defined, although may be zero
    ALLOCATE(diag_noise(0:npix_cut-1))
    if ( epsil > 0.0 ) then
        diag_noise = epsil
    else 
        diag_noise = 0.0d0
    endif

    ! If noise was smoothed, it is in map_npp, otherwise add sigma^2/hit_counts
    if (add_noise.and.(.not.do_Gsmooth) ) then
       diag_noise = diag_noise + pack(exp_noise(:,1),map_mask)
    endif
    DEALLOCATE(exp_noise)

    FORALL(i=0:npix_cut-1) map_npp(i,i) = map_npp(i,i) + diag_noise(i)

    RETURN
  END SUBROUTINE ReadPlanck_map

  SUBROUTINE REBIN_MAP(exp_data,exp_noise,exp_mask,npixin,npixout,ordering)
    USE udgrade_nr
    REAL(DP), INTENT(INOUT), DIMENSION(:,:), ALLOCATABLE :: exp_noise,exp_data,exp_mask
    INTEGER(I4B), INTENT(IN)          :: npixin,npixout
    LOGICAL,      INTENT(IN)          :: ordering ! NESTED=TRUE, RING=FALSE 

    real(DP), dimension(:,:),allocatable :: work
    real(DP), parameter                  :: good_fraction=0.1_dp
    real(DP), parameter                  :: BAD_PIXEL=-1.6375d30
    integer(I4B)                         :: nsidein,nsideout
    logical,  dimension(:,:),allocatable :: mask

    nsidein =npix2nside(npixin)
    nsideout=npix2nside(npixout)
    allocate(work(0:npixin-1,1:1))

    ! Output is the masked, rebinned, noise weighted data, noise and mask. 
    ! Rebin sets valid data to any big pixel that had at least one good highres
    ! Rebined mask is equal the fraction of good highres pixel in a lowres.
    if (do_mask) then
       ! We need to prepare the missing values in data and noise
       allocate( mask(0:npixin-1,1:1) )
       mask = ( exp_mask <  0.5_dp )      ! Input should be 0 or 1
       where( mask ) exp_data=BAD_PIXEL   ! mark masked pixels as bad

       ! Now we degrade the mask itself
       work = exp_mask
       call REBIN_FULL_SKY_WORK_TO(exp_mask)
    endif

    if (add_noise) then
       ! Average data wth inverse noise weights, first step
       if ( do_mask ) where( mask ) exp_noise=1.0_dp
       work = exp_data/exp_noise
       call REBIN_FULL_SKY_WORK_TO(exp_data)

       ! Noise is coadded in inverse
       work=1.0_dp/exp_noise
       if ( do_mask ) where( mask ) work=BAD_PIXEL ! mark bad pixels
       call REBIN_FULL_SKY_WORK_TO(exp_noise)
       ! reweight the data with inverse averaged noise
       exp_data=exp_data/exp_noise
       ! The noise is extensive quantity, undo averaging with final inverse
       exp_noise=(real(npixout,DP)/real(npixin,DP))/exp_noise
       if ( do_mask ) where( mask ) exp_noise=exp_noise/exp_mask
    else
       work = exp_data
       call REBIN_FULL_SKY_WORK_TO(exp_data)
    endif

    deallocate(mask,work)

    if ( do_mask ) then
       where ( exp_mask < good_fraction )
          exp_data = 0.0_dp ! we zero rather than BAD_PIXEL for smoothing
          exp_noise= 0.0_dp
          exp_mask = 0.0_dp
       elsewhere
!          exp_mask = 1.0_dp
       end where
    endif

    if ( DEBUG .and. STORE_REBINNED_MAPS ) then
       write(0,*)'Rebinned maps stored in default localtion'
       fake_file='rebinned_cmb.fits'
       call Write_map(exp_data)
       if ( add_noise ) then
          fake_file='rebinned_noise.fits'
          call Write_map(exp_noise)
       endif
       if ( do_mask ) then
          fake_file='rebinned_mask.fits'
          call Write_map(exp_mask)
          write(0,*)count(exp_mask > good_fraction),' good pixels'
       endif
       if ( STOP_AFTER_STORING_REBINNED_MAPS ) stop
    endif

    return

    CONTAINS
        SUBROUTINE REBIN_FULL_SKY_WORK_TO(map_out)
        real(DP), intent(inout), allocatable, dimension(:,:) :: map_out

           deallocate(map_out)
           allocate(map_out(0:npixout-1,1:1))
           if (ordering) then
              call udgrade_nest(work,nsidein,map_out,nsideout)
           else
              call udgrade_ring(work,nsidein,map_out,nsideout)
           endif
           return
        END SUBROUTINE REBIN_FULL_SKY_WORK_TO
  END SUBROUTINE REBIN_MAP

SUBROUTINE healpix_euler(alpha,beta,gamma,amat)
! see Varshalovich, Moskalev & Kershonskii
  IMPLICIT NONE
  REAL, INTENT(in) :: alpha,beta,gamma
  REAL :: s1,s2,s3,c1,c2,c3
  REAL, DIMENSION(1:3,1:3) :: m1,m2,m3,mat
  REAL, DIMENSION(1:3,1:3), INTENT(out) :: amat

  c1 = COS(alpha)
  s1 = SIN(alpha)
  c2 = COS(beta)
  s2 = SIN(beta)
  c3 = COS(gamma)
  s3 = SIN(gamma)
  
!!$  amat(1,1) = c1*c2*c3-s1*s3
!!$  amat(1,2) = -s1*c2*c3-c1*s3
!!$  amat(1,3) = s2*c3
!!$
!!$  amat(2,1) = c1*c2*s3+s1*c3
!!$  amat(2,2) = -s1*c2*s3+c1*c3
!!$  amat(2,3) = s2*s3
!!$  
!!$  amat(3,1) = c1*s2
!!$  amat(3,2) = s1*s2
!!$  amat(3,3) = c2


  amat(1, 1) = c1*c2*c3-s1*s3
  amat(1, 2) = -c1*c2*s3-s1*c3
  amat(1, 3) = c1*s2

  amat(2, 1) = s1*c2*c3+c1*s3
  amat(2, 2) = -s1*c2*s3+c1*c3
  amat(2, 3) = s1*s2

  amat(3, 1) = -s2*c3
  amat(3, 2) = s2*s3
  amat(3, 3) = c2

END SUBROUTINE  healpix_euler

END MODULE Topology_map_mod
