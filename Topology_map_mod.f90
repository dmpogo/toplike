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
    USE beams
    USE ALM_TOOLS
    USE RAN_TOOLS, ONLY : randgauss_boxmuller
    real(DP), intent(in)                               :: ampl
    real(DP), intent(out), allocatable, dimension(:,:) :: map

    real(DP),    allocatable, dimension(:,:)   :: random_numbers
    complex(DP), allocatable, dimension(:,:,:) :: alm
    integer(I4B)            :: i, iter_order=5
 
    if (allocated(map)) deallocate(map)
    allocate(map(0:npix_fits-1,1:1))
! if we add noise, generate random full sky noise map, else set it to zero
    if (add_map_noise) then
       do_Gsmooth=.false.
       do_mask=.false.
       CALL ReadExpData(expdata_format)
       if (npix_cut /= npix_fits) stop 'Must be full sky. Check noise masking'

       allocate(random_numbers(0:npix_fits-1,1))
       do i = 0, npix_fits-1
          random_numbers(i,1) = randgauss_boxmuller(iseed)
          map(i,1)=random_numbers(i,1)*sqrt(map_npp(i,i))
       enddo
       deallocate(random_numbers)
    else
       map=0.0_dp
    endif

! Generate random set of CTpp eigenvalues (sqrt)
    allocate(random_numbers(0:n_evalues-1,1))
    do i = 0, n_evalues - 1
       random_numbers(i,1) = randgauss_boxmuller(iseed)
       random_numbers(i,1) = random_numbers(i,1)*sqrt(ampl*CTpp_eval(i))
    enddo

! Combine correlated and noise map on full sky
    call dgemm('N','N',npix_fits,1,n_evalues,1.0_dp,CTpp_evec,npix_fits,random_numbers,n_evalues,1.0_dp,map,npix_fits)
    deallocate(random_numbers)

! Smooth full-sky map if needed
    if (beam_fwhm > 0.0_dp) then
       call collect_beams(Wl,lmax,G_fwhm=beam_fwhm,reset=.true.)
       allocate( alm(1:1,0:lmax,0:lmax) )
       call map2alm_iterative(nside,lmax,lmax,iter_order,map,alm,(/0.0_dp, 0.0_dp/),w8ring)
       call alter_alm(nside,lmax,lmax,beam_fwhm,alm,window=Wl)
       call alm2map(nside,lmax,lmax,alm,map(:,1))
       deallocate( alm )
    endif

    write(0,*)'random map has been generated'
    return
  END SUBROUTINE make_fake_map

  SUBROUTINE Write_map(heal_map)
    REAL(DP), INTENT(IN), DIMENSION(:,:)   :: heal_map
    !-----------------------------------------------------------------------
    !                        generates header
    !-----------------------------------------------------------------------
    INTEGER, PARAMETER :: nlheader = 30
    CHARACTER(LEN=80), DIMENSION(1:nlheader) :: header
    INTEGER :: iostatus
    REAL    :: nullval
    LOGICAL :: anynull,filefound

    !heal_map(:,1) = unpack(map_signal,map_mask,0.d0)

    header = ''

    CALL add_card(header,'COMMENT','-----------------------------------------------')
    CALL add_card(header,'COMMENT','     Sky Map Pixelisation Specific Keywords    ')
    CALL add_card(header,'COMMENT','-----------------------------------------------')
    CALL add_card(header,'PIXTYPE','HEALPIX','HEALPIX Pixelisation')
    CALL add_card(header,'ORDERING','RING',  'Pixel ordering scheme, either RING or NESTED')
    CALL add_card(header,'NSIDE'   ,nside,   'Resolution parameter for HEALPIX')
    CALL add_card(header,'FIRSTPIX',0,'First pixel # (0 based)')
    CALL add_card(header,'LASTPIX',npix_fits-1,'Last pixel # (0 based)')
    CALL add_card(header) ! blank line
    CALL add_card(header,'COMMENT','-----------------------------------------------')
    CALL add_card(header,'COMMENT','     Planck Simulation Specific Keywords      ')
    CALL add_card(header,'COMMENT','-----------------------------------------------')
    CALL add_card(header,'EXTNAME','''SMOOTHED DATA''')
    CALL add_card(header,'CREATOR','Carlo and Dima',        'Software creating the FITS file')
    CALL add_card(header,'VERSION','Ugly Hack',     'Version of the simulation software')
    CALL add_card(header)
    CALL add_card(header,'COMMENT','-----------------------------------------------')
    CALL add_card(header,'COMMENT','     Data Description Specific Keywords       ')
    CALL add_card(header,'COMMENT','-----------------------------------------------')
    CALL add_card(header,'INDXSCHM','IMPLICIT',' Indexing : IMPLICIT or EXPLICIT')
    CALL add_card(header,'GRAIN', 0, ' Grain of pixel indexing')
    CALL add_card(header,'COMMENT','GRAIN=0 : no indexing of pixel data                         (IMPLICIT)')
    CALL add_card(header,'COMMENT','GRAIN=1 : 1 pixel index -> 1 pixel data                     (EXPLICIT)')
    CALL add_card(header,'COMMENT','GRAIN>1 : 1 pixel index -> data of GRAIN consecutive pixels (EXPLICIT)')
    CALL add_card(header) ! blank line
!    fake_file='./fake.fits'
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
    write(0,*)size(heal_map,1),size(heal_map,2),npix_fits,nmaps, nlheader
    CALL write_bintab(heal_map, npix_fits, nmaps, header, nlheader,TRIM(ADJUSTL(fake_file)))
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
    
    INTEGER :: i,j,ordering,lcount,iter_order=5
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
    
    INTEGER :: i,j,ordering,lcount,iter_order=5
    logical :: convert_from_nested=.false.
    complex(DP), DIMENSION(:,:,:), ALLOCATABLE :: alm
    complex(DP), DIMENSION(:,:),   ALLOCATABLE :: clm
    REAL(DP), DIMENSION(:),   ALLOCATABLE :: diag_noise
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: exp_noise,exp_data,exp_mask
    
! Input must be full-sky in globally accessible files
    npix_fits=getsize_fits(map_signal_file,nmaps=nmaps,ordering=ordering)
    if ( nside2npix(nside) /= npix_fits ) then
       stop 'Mismatch between map size and expected nside'
    endif

! Allocate arrays and input necessary data. 
    ALLOCATE(exp_data(0:npix_fits-1,1:1))
    WRITE(0,'(a,a)') '> ', TRIM(ADJUSTL(map_signal_file))
    CALL input_map(TRIM(ADJUSTL(map_signal_file)),exp_data,npix_fits,1)
    exp_data=exp_data*1.0d3      ! convert to mK, assuming input is in K

! Check ordering, CTpp is probably given in RING (to do: auto-synchronized)
    if ( ordering == 0 ) then
       write(0,*)'Ordering of the input map is unknown, assumed RING'
    else if ( ordering == 2 ) then
       write(0,*)'Input converted from NESTED to RING'
       call convert_nest2ring(nside,exp_data)
       convert_from_nested=.true.
    else
       write(0,*)'Input is in RING pixelization'
    endif

    ALLOCATE(exp_noise(0:npix_fits-1,1:1))
    IF( add_noise) THEN
       ! Assumes diagonal noise, 
       ! and wmap_data to contain noise per pixel in same units as map
       WRITE(0,'(a,a)') '> ', TRIM(ADJUSTL(map_noise_file))
       CALL input_map(TRIM(ADJUSTL(map_noise_file)),exp_noise,npix_fits,1)
       exp_noise=exp_noise*1.0d3   ! convert to mK, it is in K (check !)
       exp_noise=exp_noise**2      ! make the variance
       if (convert_from_nested) call convert_nest2ring(nside,exp_noise)
    ELSE
       exp_noise = 0.0_dp
    ENDIF

    ALLOCATE(map_mask(0:npix_fits-1))
    IF(do_mask) THEN
       ALLOCATE(exp_mask(0:npix_fits-1,1:1))
       WRITE(0,'(a,a)') '> ', TRIM(ADJUSTL(map_mask_file))
       CALL input_map(TRIM(ADJUSTL(map_mask_file)),exp_mask,npix_fits,1)
       map_mask = ( exp_mask(:,1) /= 0 ) 
       DEALLOCATE(exp_mask)
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
       call map2alm_iterative(nside,lmax,lmax,iter_order,exp_data,alm,(/0.0_dp, 0.0_dp/),w8ring)
       call alter_alm(nside,lmax,lmax,beam_fwhm,alm,window=Wl)
       call alm2map(nside,lmax,lmax,alm,exp_data(:,1))
       map_signal = pack(exp_data(:,1),map_mask)
       DEALLOCATE(alm)
    ELSE
       map_signal = pack(exp_data(:,1),map_mask)
    ENDIF
    DEALLOCATE(exp_data)

!Smooth noise if needed and store it in a cut sky matrix map_npp
    ALLOCATE(map_npp(0:npix_cut-1,0:npix_cut-1))
    IF (add_noise.and.do_Gsmooth) THEN
       call getclm(clm, lcount, exp_noise, npix_fits, lmax, w8_file=w8_file)
       call smooth_clm(clm, lcount, Wl(:,1), lmax)
       call getcpp(map_npp, npix_cut, clm, lcount, nside, mask=map_mask)
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
       diag_noise = diag_noise + pack(exp_noise(:,1),map_mask)
    endif
    DEALLOCATE(exp_noise)

    FORALL(i=0:npix_cut-1) map_npp(i,i) = map_npp(i,i) + diag_noise(i)

    RETURN
  END SUBROUTINE ReadPlanck_map

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
