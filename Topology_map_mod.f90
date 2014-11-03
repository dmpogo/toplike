MODULE Topology_map_mod
  USE Topology_types
  USE FITSTOOLS
  USE HEAD_FITS
  USE PIX_TOOLS
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: make_fake_map, Write_map, ReadExpData, Read_CovNoise
  
CONTAINS

  SUBROUTINE make_fake_map(map)
    ! Works on eigenvalue decomposition of CTpp in CTpp_eval and CTpp_evec
    USE beams
    USE ALM_TOOLS
    USE RNGMOD
    real(DP), intent(out), allocatable, dimension(:,:) :: map

    type(planck_rng)        :: rng_handle
    real(DP),    allocatable, dimension(:)     :: WORK
    real(DP),    allocatable, dimension(:)     :: random_numbers
    complex(DP), allocatable, dimension(:,:,:) :: alm
    integer(I4B)            :: i, p, np, iter_order=5
 
    write(0,*)'Generating random map'
    if (allocated(map)) deallocate(map)
    allocate(map(0:npix_fits-1,1:npol))

! Generate random set of CTpp eigenvalues (sqrt)
    allocate(random_numbers(0:n_evalues-1))
    call rand_init(rng_handle,iseed)
    do i = 0, n_evalues - 1
       random_numbers(i) = rand_gauss(rng_handle)*sqrt(CTpp_eval(i))
    enddo
    
! Combine correlated and noise map on full sky
    allocate( WORK(0:ntot-1) )
    call dgemv('N',ntot,n_evalues,1.0_dp,CTpp_evec,ntot,random_numbers,1,0.0_dp,WORK,1)
    map=reshape( WORK, (/ npix_fits, npol /) )
    deallocate(WORK,random_numbers)
    write(0,*)'random map has been generated'

! Smooth full-sky map if needed
    if ( do_smooth_data ) then
       allocate( alm(1:npol,0:lmax,0:lmax) )
       call map2alm_iterative(nside,lmax,lmax,iter_order,map,alm,(/0.0_dp, 0.0_dp/),w8ring)
       call alter_alm(nside,lmax,lmax,0.0_dp,alm,window=Wl)
       call alm2map(nside,lmax,lmax,alm,map)
       write(0,*)'random map has been smoothed with window in Wl'
       deallocate( alm )
    endif

    return
  END SUBROUTINE make_fake_map

  SUBROUTINE Write_map(heal_map,map_file)
    REAL(DP), INTENT(IN), DIMENSION(:,:)   :: heal_map
    CHARACTER(*), INTENT(IN)               :: map_file
    !-----------------------------------------------------------------------
    !                        generates header
    !-----------------------------------------------------------------------
    INTEGER, PARAMETER :: nlheader = 80
    CHARACTER(LEN=80), DIMENSION(1:nlheader) :: header
    INTEGER :: iostatus, npix_loc, nmap_loc, nside_loc
    REAL    :: nullval
    LOGICAL :: anynull,filefound, polarization=.true.

    !heal_map(:,1) = unpack(map_signal,map_mask,-1.6375d30)
    npix_loc =size(heal_map,1)
    nmap_loc =size(heal_map,2)
    nside_loc=npix2nside(npix_loc)
    if ( npol == 1 ) polarization=.false.

    CALL write_minimal_header(header,'MAP',nside=nside_loc,ordering='RING',creator='Topology_make_map',coordsys='G',randseed=iseed,units='mK',nlmax=lmax,polar=polarization,fwhm_degree=map_beam_fwhm)

    INQUIRE(file=TRIM(ADJUSTL(map_file)),exist=filefound)
    IF(filefound) THEN
       open(26,file=TRIM(ADJUSTL(map_file)),form='unformatted')
       close(26,status='delete')
       write(0,*) 'Overwriting existing fits file'
    ELSE
       open(26,file=TRIM(ADJUSTL(map_file)),form='unformatted',status='new',iostat=iostatus)
       IF ( iostatus > 0 ) THEN
          write(0,*) 'Unable to open output ',TRIM(ADJUSTL(map_file))
          stop
       ELSE
          close(26,status='delete')
       ENDIF
    ENDIF
    
    write(0,*)'Writing bintab'
    write(0,*)size(heal_map,1),size(heal_map,2),npix_loc,nmap_loc,nlheader
    CALL write_bintab(heal_map, npix_loc, nmap_loc, header, nlheader,TRIM(ADJUSTL(map_file)))
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

  SUBROUTINE ReadWMAP_map()     ! Not setup for polarization
    USE ctpplm_transforms
    USE beams
    USE ALM_TOOLS
    !Global map_signal,map_npp,diag_noise,map_signal_file,map_mask_file,nside
    
    INTEGER :: i,j,ordering,nmaps,lcount,iter_order=5
    logical :: convert_from_nested=.false.
    logical,     DIMENSION(:,:),   ALLOCATABLE :: bool_mask
    complex(DP), DIMENSION(:,:,:), ALLOCATABLE :: alm
    complex(DP), DIMENSION(:,:,:), ALLOCATABLE :: clm
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
    wmap_noise = 0.0_dp
    IF( add_noise_diag .and. nmaps > 1 ) THEN
       ! Assumes diagonal noise, 
       ! and wmap_data to contain inverse noise in the last map
       where( wmap_data(:,nmaps) > 0.0 ) wmap_noise(:,1) = 1.d0/wmap_data(:,nmaps)
    ENDIF

    ALLOCATE(map_mask(0:npix_fits-1))
    IF(do_mask) THEN
       ALLOCATE(wmap_mask(0:npix_fits-1,1:1))
       WRITE(0,'(a,a)') '> ', TRIM(ADJUSTL(map_mask_file))
       CALL input_map(TRIM(ADJUSTL(map_mask_file)),wmap_mask,npix_fits,1)
       map_mask = ( wmap_mask(:,1) > 0.5 ) 
!       DEALLOCATE(wmap_mask)
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
    IF (do_smooth_data) THEN
       ALLOCATE( alm(1:1,0:lmax,0:lmax) )
       ALLOCATE( wmap_signal(0:npix_fits-1,1:1) )
       wmap_signal(:,1:1) = wmap_data(:,1:1) ! unfortunate copy to go to DP
       call map2alm_iterative(nside,lmax,lmax,iter_order,wmap_signal,alm,(/0.0_dp, 0.0_dp/),w8ring)
       call alter_alm(nside,lmax,lmax,0.0_dp,alm,window=Wl)
       call alm2map(nside,lmax,lmax,alm,wmap_signal(:,1))
       map_signal = pack(wmap_signal(:,1),map_mask)
       DEALLOCATE(alm, wmap_signal)
    ELSE
       map_signal = pack(wmap_data(:,1),map_mask)
    ENDIF
    DEALLOCATE(wmap_data)

!Smooth noise if needed and store it in a cut sky matrix map_npp
    ALLOCATE(map_npp(0:npix_cut-1,0:npix_cut-1))
    IF ( add_noise_diag.and.do_smooth_noise ) THEN
       call getclm(clm, lcount, wmap_noise, npix_fits, 1, lmax, w8_file=w8_file)
       call smooth_clm(clm, lcount, 1, Wl, lmax)
       call getcpp(map_npp, clm, lcount, 1, nside, mask=( wmap_mask > 0.5 ))
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
    if ( add_noise_diag.and.(.not.do_smooth_noise) ) then
       diag_noise = diag_noise + pack(wmap_noise(:,1),map_mask)
    endif
    DEALLOCATE(wmap_noise)

    FORALL(i=0:npix_cut-1) map_npp(i,i) = map_npp(i,i) + diag_noise(i)

    RETURN
  END SUBROUTINE ReadWMAP_map

  SUBROUTINE ReadPlanck_map()
    USE healpix_extras
    USE ALM_TOOLS
    !Global map_signal,map_npp,diag_noise,map_signal_file,map_mask_file,nside
    
    INTEGER :: i,j,ordering,ifpol,lcount,npix,nmaps,nnoise,nmasks,iter_order=5
    logical :: convert_from_nested=.false.
    logical,     DIMENSION(:,:),     ALLOCATABLE :: bool_mask
    complex(DP), DIMENSION(:,:,:),   ALLOCATABLE :: alm
    complex(DP), DIMENSION(:,:,:),   ALLOCATABLE :: clm
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: exp_noise,exp_data,exp_mask,bpp
    
! Global npol=1,2 or 3 controls data used (I, QU, or IQU) 
! Input must be full-sky in globally accessible files 
! Signal - single file with either 1 (T-only), 2 (QU-only) or 3 (IQU) maps.
! Other formats must be converted into this
    npix=getsize_fits(map_signal_file,nmaps=nmaps,ordering=ordering,polarisation=ifpol)
    if (nmaps < npol) then
       write(0,*)'nmaps=',nmaps,' < npol=',npol
       stop 'Not enough maps to do what requested'
    endif
    if (ifpol == 1 .and. nmaps < 2) stop 'For polarization need at least 2 maps in the input'
    if (nmaps == 2) stop '2 map (polarization only) format is not yet supported'
    if (ifpol == 0 .and. npol > 1) then 
       stop 'Polarization is explicitly unset, check data files'
    endif

    npix_fits=nside2npix(nside)
    if (npix > npix_fits) then
       write(0,*)'Input map will be coarsened from nside=',npix2nside(npix)
       if (ifpol == 1) stop 'However coarsening is not yet supported for polarization'
    elseif (npix < npix_fits) then
       stop 'Analysis on map coarser than CTpp is not implemented'
    endif

! Allocate arrays and input necessary data. 
    ALLOCATE(exp_data(0:npix-1,1:nmaps))
    WRITE(0,'(a,a)') '> ', TRIM(ADJUSTL(map_signal_file))
    CALL input_map(TRIM(ADJUSTL(map_signal_file)),exp_data,npix,nmaps)
    exp_data=exp_data*expdata_scale   ! convert to mK, from whatever input is
                                      ! assumes same scale for I and QU

    allocate(exp_noise(0:npix-1,1:nmaps))
    exp_noise=1.0_dp
    IF ( add_noise_diag ) THEN
       ! Assumes exp_noise to contain noise per pixel in same units as map
       write(0,'(a,a)') '> ', TRIM(ADJUSTL(map_noise_file))
       if (getsize_fits(map_noise_file,nmaps=nnoise) /= npix) &
          stop 'sizes of the signal and noise files are inconsistent'
       call input_map(TRIM(ADJUSTL(map_noise_file)),exp_noise,npix,nnoise)
       if ( nnoise == 2 .and. nmaps == 3 ) exp_mask(:,3) = exp_mask(:,2)
       exp_noise=exp_noise*expdata_scale  ! convert to mK
       exp_noise=exp_noise**2             ! make the variance
    ELSE IF ( add_noise_cov ) THEN
       ! If noise is given as covariance, it should be in the  final
       ! working pixelization. No infrastructure noise covariance coarsening.
       write(0,*)'Covariance matrix will be read in separately'
    ENDIF

    ALLOCATE(exp_mask(0:npix-1,1:nmaps))
    IF ( do_mask ) THEN
       write(0,'(a,a)') '> ', TRIM(ADJUSTL(map_mask_file))
       if (getsize_fits(map_mask_file,nmaps=nmasks) /= npix) &
          stop 'sizes of the signal and mask files are inconsistent'
       call input_map(TRIM(ADJUSTL(map_mask_file)),exp_mask,npix,nmasks)
       if ( nmasks < nmaps ) then  ! copy the mask from the last present
          forall(i=nmasks+1:nmaps) exp_mask(:,i) = exp_mask(:,nmasks)
       endif
    ELSE
       exp_mask = 1
    ENDIF

! if npol < nmaps, use mask to short out unused data
    if ( npol == 1 .and. nmaps > 1) then  ! I only
       exp_mask(:,2:nmaps) = 0.0_dp
       npol = 3 
    else if ( npol == 2 ) then ! QU only can only be done by masking 
       exp_mask(:,1) = 0.0_dp
       npol = 3 
    endif

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
       if ( add_noise_diag ) call convert_nest2ring(nside,exp_noise)
    endif

! Convert coarse grained mask to  boolean
    ALLOCATE(bool_mask(0:npix_fits-1,1:nmaps))
    bool_mask = .false.
    ! This threshold interplace with one in REBIN, if mask is not scaled to 0-1
    ! It was introduced to work agains smoothing of masked map issues
    ! to mask edge pixels - that's why a bit over 50% is a good criterium
    ! Allow for different criterium for I and QU
    if (nmaps /= 2) then   ! I
       bool_mask(:,1) = ( exp_mask(:,1) >= 0.6_dp )  
    endif

    if (nmaps /= 1) then   ! QU
       bool_mask(:,nmaps-1:nmaps) = ( exp_mask(:,nmaps-1:nmaps) >= 0.6_dp )
    endif
    DEALLOCATE(exp_mask)

! Data in, set ntot and store it in the linear arrays for analysis
    ntot = npix_fits*nmaps

    allocate(map_mask(0:ntot-1))
    map_mask = reshape(bool_mask, (/ ntot /))

! Count unmasked pixels 
    npix_cut = count(bool_mask)
    write(0,'(a,i7,a,i7)') 'Found ',npix_cut,' unmasked pixels from ',ntot
    if (npix_cut == 0) STOP 'All pixels are masked'
  
! Smooth (in lm space) and pack signal

    IF (do_smooth_data) THEN
       ALLOCATE( alm(1:nmaps,0:lmax,0:lmax) )
       call map2alm_iterative(nside,lmax,lmax,iter_order,exp_data,alm,(/0.0_dp, 0.0_dp/),w8ring)
       write(0,*) 'Monopole and Dipole ',alm(1:nmaps,0:1,0:1)
       alm(1:nmaps,0:1,0:1) = cmplx(0.0_dp,0.0_dp)
       call alter_alm(nside,lmax,lmax,0.0_dp,alm,window=Wl)
       call alm2map(nside,lmax,lmax,alm,exp_data)
       DEALLOCATE(alm)
    ENDIF

    ALLOCATE(diag_noise(0:npix_cut-1))
    IF ( add_noise_diag ) THEN
       diag_noise=PackMasked(exp_noise,bool_mask)
       DEALLOCATE(exp_noise)
    ELSE
      diag_noise = 0.0_dp
    ENDIF

    ALLOCATE(map_signal(0:npix_cut-1))
    map_signal=PackMasked(exp_data,bool_mask)
    DEALLOCATE(exp_data,bool_mask)

    RETURN
  END SUBROUTINE ReadPlanck_map

  SUBROUTINE Read_CovNoise()
    USE healpix_extras
    USE ctpplm_transforms
    USE beams
    USE ALM_TOOLS
    USE ct_io
    !Global  map_npp,diag_noise,,map_mask,nside,npix_fits,map_noise_file
    ! Assumes map_mask to be defined, everything is in RING
    
    INTEGER                                      :: i,j,lcount
    complex(DP), DIMENSION(:,:,:),   ALLOCATABLE :: clm

! Allocate arrays and input necessary data. 

    allocate( map_npp(0:ntot-1,0:ntot-1) )
    IF ( add_noise_cov ) THEN
       ! If noise is given as covariance, it should be in the  final
       ! working pixelization. No infrastructure noise covariance coarsening.
       call ReadCTpp(map_noise_file,map_npp,npix_fits,npol)
       map_npp=map_npp*(expdata_scale**2)
    ELSE
       map_npp=0.0_dp
    ENDIF

!Smooth noise if needed and store it in a cut sky matrix map_npp

    IF (  add_noise_cov ) THEN
       IF ( do_smooth_noise ) THEN
          ! cov matrix smoothed and masked on reconstruction
          call getclm(clm, lcount, map_npp, npix_fits, npol, lmax, w8_file=w8_file)
          call smooth_clm(clm, lcount, npol, Wl, lmax)
          call getcpp(map_npp, clm, lcount, npol, nside)
          deallocate( clm )
          write(0,*) 'Noise covariance has been smoothed'
       ENDIF

       IF ( do_mask ) THEN
          i = 0
          do j=0,ntot-1
             if ( map_mask(j) ) then
                map_npp(0:npix_cut-1,i) = pack( map_npp(:,j), map_mask )
                i = i+1
             endif
          enddo
       ENDIF
    ENDIF

! if diagonal noise was defined, add it unsmoothed
    IF ( add_noise_diag ) THEN
       FORALL(i=0:npix_cut-1) map_npp(i,i) = map_npp(i,i) + diag_noise(i)
    ENDIF

! if we have regularizer, add it
    IF ( epsil > 0.0_dp ) THEN
       FORALL(i=0:npix_cut-1) map_npp(i,i) = map_npp(i,i) + epsil
    ENDIF

    RETURN
  END SUBROUTINE Read_CovNoise

  SUBROUTINE REBIN_MAP(exp_data,exp_noise,exp_mask,npixin,npixout,ordering)
    USE udgrade_nr
    REAL(DP), INTENT(INOUT), DIMENSION(:,:), ALLOCATABLE :: exp_noise,exp_data,exp_mask
    INTEGER(I4B), INTENT(IN)          :: npixin,npixout
    LOGICAL,      INTENT(IN)          :: ordering ! NESTED=TRUE, RING=FALSE 

    real(DP), dimension(:,:),allocatable :: work
    real(DP), parameter                  :: good_fraction=0.1_dp
    real(DP), parameter                  :: BAD_PIXEL=-1.6375d30
    integer(I4B)                         :: nsidein,nsideout,nmapsin
    logical,  dimension(:,:),allocatable :: mask

    nsidein =npix2nside(npixin)
    nsideout=npix2nside(npixout)
    nmapsin =size(exp_data,2)
    allocate(work(0:npixin-1,1:nmapsin))

    ! Output is the masked, rebinned, noise weighted data, noise and mask. 
    ! Rebin sets valid data to any big pixel that had at least one good highres
    ! Rebined mask is equal the fraction of good highres pixel in a lowres.
    if (do_mask) then
       ! We need to prepare the missing values in data and noise
       allocate( mask(0:npixin-1,1:nmapsin) )
       mask = ( exp_mask <  0.5_dp )      ! Input should be 0 or 1
       where( mask ) exp_data=BAD_PIXEL   ! mark masked pixels as bad

       ! Now we degrade the mask itself
       work = exp_mask
       call REBIN_FULL_SKY_WORK_TO(exp_mask)
    endif

    if ( add_noise_diag ) then
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
       call Write_map(exp_data,'rebinned_cmb_fits')
       if ( add_noise_diag ) then
          call Write_map(exp_noise,'rebinned_noise.fits')
       endif
       if ( do_mask ) then
          call Write_map(exp_mask,'rebinned_mask_fits')
          write(0,*)count(exp_mask > good_fraction),' good pixels'
       endif
       if ( STOP_AFTER_STORING_REBINNED_MAPS ) stop
    endif

    return

    CONTAINS
        SUBROUTINE REBIN_FULL_SKY_WORK_TO(map_out)
        real(DP), intent(inout), allocatable, dimension(:,:) :: map_out

           deallocate(map_out)
           allocate(map_out(0:npixout-1,1:nmapsin))
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
