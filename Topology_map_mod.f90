MODULE Topology_map_mod
  USE TOPOLOGY_TYPES
  USE FITSTOOLS
  USE HEAD_FITS
  USE PIX_TOOLS
  USE RAN_TOOLS, ONLY : randgauss_boxmuller
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: make_fake_map, ReadWmap_map, Read_w8ring
  
CONTAINS

  SUBROUTINE make_fake_map(ampl)
    IMPLICIT NONE

    REAL(DP), INTENT(IN)                 ::ampl

    REAL(DP), DIMENSION(0:npix_cut*10-1) :: ework
    REAL(DP), DIMENSION(0:npix_cut-1)    :: evals
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: mata,matb
    REAL,     DIMENSION(:,:), ALLOCATABLE :: heal_map
    REAL(DP), DIMENSION(:), ALLOCATABLE :: map_cut, map_cut2
    INTEGER :: INFO,i,iring,j


    INTEGER, PARAMETER :: nlheader = 30
    CHARACTER(LEN=80), DIMENSION(1:nlheader) :: header
    REAL    :: nullval
    LOGICAL :: anynull,filefound

    REAL(DP), ALLOCATABLE, DIMENSION(:) :: WORKNEL,D
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: U,VT
!    GOTO 696
    IF(iseed == 0) THEN
       CALL SYSTEM_CLOCK(count = iseed)
       IF (MOD(iseed, 2) .EQ. 0) iseed = iseed + 1
    ENDIF

    ALLOCATE(mata(0:npix_cut-1,0:npix_cut-1))
    ALLOCATE(matb(0:npix_cut-1,0:npix_cut-1))

    !calculate Hermitean square root of correlation matrix. Spoils
    mata=CTpp*ampl
!    IF(SVD) THEN
!       WRITE(0,*)"Doing SVD MAP"
!       ALLOCATE(U(0:npix_cut-1,0:npix_cut-1))
!       ALLOCATE(VT(0:npix_cut-1,0:npix_cut-1))
!       ALLOCATE(WORKNEL(0:5*npix_cut))
!       INFO = 0
!       DO i = 0, npix_cut-1
!          DO j = i, npix_cut-1
!             mata(i,j) = mata(j,i)
!          ENDDO
!       ENDDO
!!    Do general SVD 
!       CALL DGESVD('A','A',npix_cut,npix_cut,mata,npix_cut,evals,&
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
!    ELSE !Cholesky methode
       CALL dsyev('V', 'L', npix_cut, mata, npix_cut, evals, ework, &
                 & 10*npix_cut,INFO)
       IF(INFO == 0) THEN
          DO i=0,npix_cut-1
             IF (evals(i) < 0.d0) THEN
                WRITE(0,'(a,i7,1x,1pe12.4)') 'Warning negative eigenvalue ', i, evals(i)
                evals(i) = 0.d0
             ENDIF
             mata(:,i) = mata(:,i) * evals(i)**0.25d0
          ENDDO
       ELSE
          WRITE(0,*)'INFO=',INFO
          STOP 'Failed on DSYEV'
       ENDIF

       CALL dsyrk('L', 'N', npix_cut, npix_cut, 1.0d0, mata, npix_cut,&
            &0.0d0, matb, npix_cut)
!    ENDIF

    DEALLOCATE(mata)

    ALLOCATE(map_cut(0:npix_cut-1))
    ALLOCATE(map_cut2(0:npix_cut-1))
    ALLOCATE(heal_map(0:npix_fits-1,1))
    
    DO i = 0, npix_cut - 1
       map_cut(i) = DBLE(randgauss_boxmuller(iseed))
    ENDDO

!   Generate random realization from purly gaussian map
    
    CALL dsymv('L',npix_cut,1.0d0,matb,npix_cut,map_cut,1,0.0d0,map_cut2,1)

    IF(add_map_noise) THEN !If add_noise=.FALSE => add_map_noise=.FALSE.
       DO i = 0, npix_cut - 1
          map_cut2(i) = map_cut2(i) + DBLE(randgauss_boxmuller(iseed)*SQRT(wmap_noise(i)))
       ENDDO
    ENDIF
    
    !Uncomment to print real map with cut
!696 CONTINUE
!    ALLOCATE(map_cut2(0:npix_cut-1))
!    ALLOCATE(heal_map(0:npix_fits-1,1))
!    map_cut2=wmap_signal 
    heal_map(:,1) = 0.d0
    DO i=0,npix_cut-1
       CALL vec2pix_ring(nside, DBLE(wmap_qhat(:,i)), iring)
       heal_map(iring,1) = map_cut2(i)
    ENDDO

    !-----------------------------------------------------------------------
    !                        generates header
    !-----------------------------------------------------------------------
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
    INQUIRE(file=TRIM(ADJUSTL(fake_file)),exist=filefound)
    IF(filefound) THEN
       OPEN(26,file=TRIM(ADJUSTL(fake_file)),form='unformatted')
       CLOSE(26,status='delete')
       WRITE(0,*) 'Overwriting existing fits file'
    ENDIF

    write(0,*)'Writing bintab'
    write(0,*)size(heal_map,1),size(heal_map,2),npix_fits,nmaps, nlheader
    CALL write_bintab(heal_map, npix_fits, nmaps, header, nlheader,TRIM(ADJUSTL(fake_file)))
    write(0,*)'Done'

    DEALLOCATE(map_cut,matb,heal_map)

    RETURN
  END SUBROUTINE make_fake_map

  SUBROUTINE Read_w8ring()
    IMPLICIT NONE
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
    !wmap_signal_file, wmap_noise_file, wmap_mask_file, nside
    IMPLICIT NONE
    
    INTEGER :: i,j
    DOUBLE PRECISION :: vec(0:2)
    REAL(DP), DIMENSION(:,:), ALLOCATABLE   :: map_beam,wmap_npp_temp
    REAL(DP), DIMENSION(:), ALLOCATABLE   :: signal,map_signal
    REAL(SP), DIMENSION(:,:), ALLOCATABLE   :: map, map_noise, map_mask 
    
    !maps must be in ring format
    ! Allocate arrays
    ALLOCATE(map(0:npix_fits-1,1:nmaps))
    ALLOCATE(map_signal(0:npix_fits-1))
    ALLOCATE(map_noise(0:npix_fits-1,1:nmaps))
    ALLOCATE(map_mask(0:npix_fits-1,1:nmaps))
    ALLOCATE(map_beam(0:npix_fits-1,0:npix_fits-1))

    WRITE(0,'(a,a)') '> ', TRIM(ADJUSTL(wmap_signal_file))
    WRITE(0,'(a,a)') '> ', TRIM(ADJUSTL(wmap_noise_file))
    CALL input_map(TRIM(ADJUSTL(wmap_signal_file)),map,npix_fits,nmaps)

    IF(add_noise) THEN
       WRITE(0,'(a,a)') '> ', TRIM(ADJUSTL(wmap_noise_file))
       CALL input_map(TRIM(ADJUSTL(wmap_noise_file)),map_noise,npix_fits,nmaps)
    ELSE
       map_noise(:,1) = 0.d0
    ENDIF

    IF(do_mask) THEN
       WRITE(0,'(a,a)') '> ', TRIM(ADJUSTL(wmap_mask_file))
       CALL input_map(TRIM(ADJUSTL(wmap_mask_file)),map_mask,npix_fits,nmaps)
    ELSE
       map_mask(:,1) = 1.d0
    ENDIF

    IF(do_smooth) THEN
       WRITE(0,'(a,a)') '> ', TRIM(ADJUSTL(beam_file))
       open(100,file=TRIM(ADJUSTL(beam_file)), status='unknown', form='unformatted')
       read(100)map_beam
    ENDIF
    


    ALLOCATE(wmap_mask(0:npix_fits-1))

    npix_cut = 0
    DO i=0,npix_fits-1
       IF(map_mask(i,1) /=0.0) then
         npix_cut = npix_cut + 1
         wmap_mask(i)=.TRUE.
       ELSE
         wmap_mask(i)=.FALSE.
       ENDIF
    ENDDO

    WRITE(0,'(a,i7,a,i7)') 'Found ', npix_cut, ' unmasked pixels from ', npix_fits
    if (npix_cut == 0) STOP 'All pixels are masked'
     
    ! Allocate output arrays
    ALLOCATE(wmap_signal(0:npix_cut-1))
    ALLOCATE(wmap_noise(0:npix_cut-1))
    ALLOCATE(diag_noise(0:npix_cut-1))
    ALLOCATE(wmap_qhat(0:2,0:npix_cut-1))
    ALLOCATE(wmap_npp_temp(0:npix_fits-1,0:npix_fits-1))
    ALLOCATE(wmap_npp(0:npix_cut-1,0:npix_cut-1))
!Smooth signal add a check or fix so it won't double smooth (works!!!!)
    !IF (do_smooth) THEN
    !   ALLOCATE(signal(0:npix_fits-1))
    !   signal(:) = map(:,1)
    !   CALL DGEMM('N','N',npix_fits,1,npix_fits,1.0d0,map_beam,npix_fits,signal,npix_fits,0,map_signal,npix_fits)
    !   DEALLOCATE(signal)
    !ELSE
       map_signal(:) = map(:,1)
    !ENDIF

!Smooth noise
    IF (add_noise.and.do_smooth) THEN
       DO i=0,npix_fits-1
          map_beam(:,i)=map_beam(:,i)/sqrt(map_noise(i,1))
       ENDDO
    ! wmap_npp = map_beam x transpose(map_beam)
       call DSYRK('L', 'N', npix_fits, npix_fits, 1.0d0, map_beam, npix_fits, 0.0d0, wmap_npp_temp, npix_fits)
    ENDIF

    npix_cut = 0
    IF (add_noise) THEN
       DO i=0,npix_fits-1
          IF(map_mask(i,1) /= 0.0) THEN
             IF (do_smooth) THEN
                wmap_npp_temp(:,npix_cut) = wmap_npp_temp(:,i)
             ENDIF
             wmap_signal(npix_cut) = map_signal(i)
             !the maps actually store the coadding weight factor which
             !is 1/noise variance of each pixel
             wmap_noise(npix_cut) = 1.d0/map_noise(i,1)
             CALL pix2vec_ring(nside, i, vec)
             wmap_qhat(:,npix_cut) = REAL(vec(:))
             npix_cut = npix_cut + 1
          ENDIF
       ENDDO
       IF (do_smooth) THEN
          npix_cut = 0
          DO i=0,npix_fits-1
             IF (map_mask(i,1) /= 0.0) THEN
                wmap_npp_temp(npix_cut,:) = wmap_npp_temp(i,:)
                npix_cut = npix_cut + 1
             ENDIF
          ENDDO
       ENDIF
    ELSE
       DO i=0,npix_fits-1
          IF(map_mask(i,1) /=0.0) THEN
             wmap_signal(npix_cut) = map_signal(i)
             CALL pix2vec_ring(nside, i, vec)
             wmap_qhat(:,npix_cut) = REAL(vec(:))
             npix_cut = npix_cut + 1
          ENDIF
       ENDDO
    ENDIF

    IF ((do_smooth.and.add_noise).or.(epsil /= 0.0)) THEN
       DO i=0, npix_cut-1
          DO j=0, npix_cut-1
             wmap_npp(i,j)=wmap_npp_temp(i,j)
          ENDDO
       ENDDO

       diag_noise(:) = epsil
    ELSE
       diag_noise(:) = wmap_noise(:)
    ENDIF

    
    ! write(0,'(2(e10.5,1x))') (wmap_noise(i),wmap_npp(i,i),i=0,npix_cut-1)
    ! Deallocate read-in arrays
    DEALLOCATE(map,map_noise,map_mask,map_beam,wmap_npp_temp,map_signal)
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
