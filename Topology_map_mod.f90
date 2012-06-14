MODULE Topology_map_mod
  USE TOPOLOGY_TYPES
  USE FITSTOOLS
  USE HEAD_FITS
  USE PIX_TOOLS
  USE RAN_TOOLS, ONLY : randgauss_boxmuller
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: make_fake_map, WriteWmap_map, ReadWmap_map, Read_w8ring
  
CONTAINS

  SUBROUTINE make_fake_map(ampl)
    REAL(DP), INTENT(IN)                 ::ampl

    REAL(DP), DIMENSION(0:npix_cut*10-1)  :: ework
    REAL(DP), DIMENSION(0:npix_cut-1)     :: evals
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: matb
    REAL(DP), DIMENSION(:), ALLOCATABLE   :: map_cut
    INTEGER :: INFO,i,iring,j


    REAL(DP), ALLOCATABLE, DIMENSION(:) :: WORKNEL,D
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: U,VT
!    GOTO 696

    IF(iseed == 0) THEN
       CALL SYSTEM_CLOCK(count = iseed)
       IF (MOD(iseed, 2) .EQ. 0) iseed = iseed + 1
    ENDIF

    ALLOCATE(matb(0:npix_cut-1,0:npix_cut-1))

! calculate Hermitean square root of correlation matrix. Spoils
    CNTpp=CTpp*ampl
    IF(add_noise.and.do_smooth) THEN
      CNTpp=CNTpp+map_npp
    ENDIF
    FORALL(i=0:npix_cut-1)  CNTpp(i,i) = CNTpp(i,i) + diag_noise(i)
!    IF(SVD) THEN
!       WRITE(0,*)"Doing SVD MAP"
!       ALLOCATE(U(0:npix_cut-1,0:npix_cut-1))
!       ALLOCATE(VT(0:npix_cut-1,0:npix_cut-1))
!       ALLOCATE(WORKNEL(0:5*npix_cut))
!       INFO = 0
!       DO i = 0, npix_cut-1
!          DO j = i, npix_cut-1
!             CNTpp(i,j) = CNTpp(j,i)
!          ENDDO
!       ENDDO
!!    Do general SVD 
!       CALL DGESVD('A','A',npix_cut,npix_cut,CNTpp,npix_cut,evals,&
!                  & U,npix_cut,VT,npix_cut,WORKNEL,5*npix_cut,INFO)
!       IF(INFO/=0) THEN
!          write(0,*) "DGESVD info=", INFO
!          STOP 'Error SVD DGESVD'
!       ENDIF
!       DO i = 0, mode_number
!          IF(abs(evals(i)) == 0.0) THEN
!             VT(i,:) = 0.0d0
!          ELSE
!             VT(i,:) = VT(i,:)*SQRT(evals(i))
!          ENDIF
!       ENDDO
!       IF(mode_number<npix_cut-1) THEN
!          DO i = mode_number+1, npix_cut-1
!             VT(i,:) = 0.0d0
!          ENDDO
!       ENDIF
!
!       CALL DGEMM('N','N',npix_cut,npix_cut,npix_cut,1.0d0,U,npix_cut,&
!                 & VT,npix_cut,0.0d0,matb,npix_cut)
!       DEALLOCATE(U)
!       DEALLOCATE(VT)
!       DEALLOCATE(WORKNEL)
!    ELSE !Cholesky method
       CALL dsyev('V', 'L', npix_cut, CNTpp, npix_cut, evals, ework, &
                 & 10*npix_cut,INFO)
       IF(INFO == 0) THEN
          DO i=0,npix_cut-1
             IF (evals(i) < 0.d0) THEN
                WRITE(0,'(a,i7,1x,1pe12.4)') 'Warning negative eigenvalue ', i, evals(i)
                evals(i) = 0.d0
             ENDIF
             CNTpp(:,i) = CNTpp(:,i) * evals(i)**0.25d0
          ENDDO
       ELSE
          WRITE(0,*)'INFO=',INFO
          STOP 'Failed on DSYEV'
       ENDIF

       CALL dsyrk('L', 'N', npix_cut, npix_cut, 1.0d0, CNTpp, npix_cut,&
            &0.0d0, matb, npix_cut)
!    ENDIF

    ALLOCATE(map_cut(0:npix_cut-1))
    
    DO i = 0, npix_cut - 1
       map_cut(i) = DBLE(randgauss_boxmuller(iseed))
    ENDDO

!   Generate random realization for correlated gaussian map
    
    CALL dsymv('L',npix_cut,1.0d0,matb,npix_cut,map_cut,1,0.0d0,map_signal,1)
    
!   Uncomment to print real map with cut
!696 CONTINUE
!   ALLOCATE(map_cut2(0:npix_cut-1))
!   ALLOCATE(heal_map(0:npix_fits-1,1))
!   map_cut2=wmap_signal 
!   map_cut2=0.0d0 
!   map_cut2(600)=wmap_noise(600)

    DEALLOCATE(map_cut,matb)

    return
  END SUBROUTINE make_fake_map

  SUBROUTINE WriteWMAP_map()
    !-----------------------------------------------------------------------
    !                        generates header
    !-----------------------------------------------------------------------
    REAL(SP), DIMENSION(:,:), ALLOCATABLE    :: heal_map
    INTEGER, PARAMETER :: nlheader = 30
    CHARACTER(LEN=80), DIMENSION(1:nlheader) :: header
    INTEGER :: iostatus
    REAL    :: nullval
    LOGICAL :: anynull,filefound

    ALLOCATE(heal_map(0:npix_fits-1,1))
    heal_map(:,1) = unpack(map_signal,map_mask,0.d0)

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
          close(26)
       ENDIF
    ENDIF
    
    write(0,*)'Writing bintab'
    write(0,*)size(heal_map,1),size(heal_map,2),npix_fits,nmaps, nlheader
    CALL write_bintab(heal_map, npix_fits, nmaps, header, nlheader,TRIM(ADJUSTL(fake_file)))
    write(0,*)'Done'

    DEALLOCATE(heal_map)

    RETURN
  END SUBROUTINE WriteWMAP_map

  SUBROUTINE Read_w8ring()

    real(DP)    :: nullval
    logical     :: found,anynull
    
    allocate(w8ring(1:2*nside,1:1))

    write(0,*)w8_file
    INQUIRE(file=TRIM(w8_file),exist=found)
    if (found == .TRUE.) then
        call read_bintab(w8_file,w8ring,2*nside,1,nullval,anynull)
        write(0,*)w8ring
        w8ring=1.d0+w8ring
    else
        w8ring=1.d0
        write(0,*), 'Ring weights file is not found, using uniform weights'
    endif

  END SUBROUTINE Read_w8ring  


  SUBROUTINE ReadWMAP_map()
    !Global map_signal, map_signal_file, map_noise_file, map_mask_file, nside
    
    INTEGER :: i,j
    REAL(DP), DIMENSION(:,:), ALLOCATABLE   :: wmap_beam
    REAL(DP), DIMENSION(:),   ALLOCATABLE   :: wmap_signal
    REAL(SP), DIMENSION(:,:), ALLOCATABLE   :: wmap_data, wmap_noise, wmap_mask 
    
! Allocate arrays and input necessary data. 
! Input must be full-sky in ring format in globally accessible files
    ALLOCATE(wmap_data(0:npix_fits-1,1:nmaps))
    ALLOCATE(wmap_noise(0:npix_fits-1,1:nmaps))
    ALLOCATE(wmap_mask(0:npix_fits-1,1:nmaps))
    ALLOCATE(wmap_beam(0:npix_fits-1,0:npix_fits-1))

    WRITE(0,'(a,a)') '> ', TRIM(ADJUSTL(map_signal_file))
    WRITE(0,'(a,a)') '> ', TRIM(ADJUSTL(map_noise_file))
    CALL input_map(TRIM(ADJUSTL(map_signal_file)),wmap_data,npix_fits,nmaps)

! Assumes diagonal noise, and wmap_noise to contain hit counts
    IF(add_noise) THEN
       WRITE(0,'(a,a)') '> ', TRIM(ADJUSTL(map_noise_file))
       CALL input_map(TRIM(ADJUSTL(map_noise_file)),wmap_noise,npix_fits,nmaps)
       WHERE( wmap_noise > 0.0 ) wmap_noise = 1./wmap_noise
    ELSE
       wmap_noise(:,1) = 0.0
    ENDIF

    IF(do_mask) THEN
       WRITE(0,'(a,a)') '> ', TRIM(ADJUSTL(map_mask_file))
       CALL input_map(TRIM(ADJUSTL(map_mask_file)),wmap_mask,npix_fits,nmaps)
    ELSE
       wmap_mask(:,1) = 1.0
    ENDIF

    IF(do_smooth) THEN
       WRITE(0,'(a,a)') '> ', TRIM(ADJUSTL(beam_file))
       open(100,file=TRIM(ADJUSTL(beam_file)), status='unknown', form='unformatted')
       read(100)wmap_beam
    ENDIF

! data is in, now process it and store in global arrays

! Count unmasked pixels and store the mask
    ALLOCATE(map_mask(0:npix_fits-1))
    map_mask = ( wmap_mask(:,1) /= 0 ) 
    npix_cut = count(map_mask)

    write(0,'(a,i7,a,i7)') 'Found ',npix_cut,' unmasked pixels from ',npix_fits
    if (npix_cut == 0) STOP 'All pixels are masked'

!   Pack wmap_beam into (npix_cut,npix_fits) size
    FORALL(i=0:npix_fits-1) wmap_beam(0:npix_cut-1,i)=pack(wmap_beam(:,i),map_mask)

!Smooth and pack the signal.
    ALLOCATE(map_signal(0:npix_cut-1))
    IF (do_smooth) THEN
       ALLOCATE(wmap_signal(0:npix_fits-1))
       wmap_signal(:) = wmap_data(:,1)
       CALL DGEMM('N','N',npix_cut,1,npix_fits,1.0d0,wmap_beam,npix_fits,wmap_signal,npix_fits,0.d0,map_signal,npix_cut)
       DEALLOCATE(wmap_signal)
    ELSE
       map_signal = pack(wmap_data(:,1),map_mask)
    ENDIF


!Smooth noise if needed and store it
    IF (add_noise.and.do_smooth) THEN
       ALLOCATE(map_npp(0:npix_cut-1,0:npix_cut-1))
       ! map_npp(cut,cut) = map_beam(cut,fits) x transpose(map_beam)(fits,cut)
       FORALL(i=0:npix_fits-1) wmap_beam(0:npix_cut-1,i)=wmap_beam(0:npix_cut-1,i)*sqrt(wmap_noise(i,1))
       call DSYRK('L','N',npix_cut,npix_fits,1.0d0,wmap_beam,npix_fits,0.0d0,map_npp,npix_cut)
    ENDIF

! diagonal needs always to be defined, although may be zero
    ALLOCATE(diag_noise(0:npix_cut-1))
    if ( epsil > 0.0 ) then
        diag_noise = epsil
    else 
        diag_noise = 0.0d0
    endif

    ! If noise was smoothed, it is in map_npp, otherwise add 1/hit_counts
    if (add_noise.and.(.not.do_smooth) ) then
       diag_noise = diag_noise + pack(wmap_noise(:,1),map_mask)
    endif

    ! write(0,'(2(e10.5,1x))') (wmap_noise(i),wmap_npp(i,i),i=0,npix_cut-1)
    ! Deallocate read-in arrays
    ! DEALLOCATE(wmap_data,wmap_noise,wmap_mask,wmap_beam,wmap_signal)
    RETURN
  END SUBROUTINE ReadWMAP_map

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
