MODULE Topology_map_mod_nel
  USE TOPOLOGY_TYPES
  USE Topology_map_mod
  USE FITSTOOLS
  USE HEAD_FITS
  USE PIX_TOOLS
  USE RAN_TOOLS, ONLY : randgauss_boxmuller
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: make_fake_mode_map, ReadWmap_map, Read_w8ring
  
CONTAINS

  SUBROUTINE make_fake_mode_map(ampl)
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
    CHARACTER(LEN=80) :: modefile
    REAL    :: nullval
    LOGICAL :: anynull,filefound

    REAL(DP), ALLOCATABLE, DIMENSION(:) :: WORKNEL,D
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: U,VT
    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:) :: LOSTMAPS
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: LOSTMAPSREND
    INTEGER :: nk,nk1,nk2
    
    IF(iseed == 0) THEN
       CALL SYSTEM_CLOCK(count = iseed)
       IF (MOD(iseed, 2) .EQ. 0) iseed = iseed + 1
    ENDIF

    ALLOCATE(mata(0:npix_cut-1,0:npix_cut-1))
    ALLOCATE(matb(0:npix_cut-1,0:npix_cut-1))

    !calculate Hermitean square root of correlation matrix. Spoils

    mata=CTpp*exp(ampl)
    
   IF(SVD) THEN
       WRITE(0,*)"Doing SVD MODE MAP"
       ALLOCATE(U(0:npix_cut-1,0:npix_cut-1))
       ALLOCATE(VT(0:npix_cut-1,0:npix_cut-1))
       ALLOCATE(WORKNEL(0:5*npix_cut))
       ALLOCATE(LOSTMAPS(0:npix_cut-1,0:npix_cut-1,1:10))
       INFO = 0
       DO i = 0, npix_cut-1
          DO j = i, npix_cut-1
             mata(i,j) = mata(j,i)
          ENDDO
       ENDDO
!    Do general SVD 
       CALL DGESVD('A','A',npix_cut,npix_cut,mata,npix_cut,evals,&
                  & U,npix_cut,VT,npix_cut,WORKNEL,5*npix_cut,INFO)
       IF(INFO/=0) THEN
          write(0,*) "DGESVD info=", INFO
          STOP 'Error SVD DGESVD'
       ENDIF
       !Visualization of lost modes
       nk=1
       !DO i = 0, npix_cut-1
       !   WRITE(0,*)i,evals(i)
       !ENDDO
       DO i = 0, npix_cut-1
          IF (nk>10) GOTO 9000
          IF(i>mode_number) THEN
!          IF(evals(i)<1.0E-8) THEN
           WRITE(0,*)'evals=',evals(i)
           IF(nk>0) THEN
             WRITE(0,*)'!!!!!!!!!!!!!!!!!!!mod =',nk
             DO nk1 = 0, npix_cut-1
                DO nk2 = 0, npix_cut-1
                       LOSTMAPS(nk1,nk2,nk) = U(nk1,i)*VT(i,nk2)
                ENDDO
             ENDDO
           ENDIF
           nk=nk+1
          ENDIF
       ENDDO
9000   CONTINUE
       DEALLOCATE(mata)

       ALLOCATE(map_cut(0:npix_cut-1))
       ALLOCATE(map_cut2(0:npix_cut-1))
       ALLOCATE(heal_map(0:npix_fits-1,1))
!       nk=nk-5
       ALLOCATE(LOSTMAPSREND(0:npix_fits-1,1:nk-1)) 
       DO nk2 = 1, nk-1
          map_cut2(:)=LOSTMAPS(:,1000,nk2)
          heal_map(:,1) = 0.d0
          DO i=0,npix_cut-1
             CALL vec2pix_ring(nside, DBLE(wmap_qhat(:,i)), iring)
             LOSTMAPSREND(iring,nk2) = map_cut2(i)
          ENDDO
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

    modefile="./modes_maps/modes1_10_"//TRIM(ADJUSTL(infile(25:100)))//"modenum.fits"
    INQUIRE(file=modefile,exist=filefound)
    IF(filefound) THEN
       OPEN(26,file=modefile,form='unformatted')
       CLOSE(26,status='delete')
       WRITE(0,*) 'Overwriting existing fits file'
    ENDIF
    write(0,*)'Writing bintab'

    write(0,*)size(LOSTMAPSREND,1),size(LOSTMAPSREND,2),npix_fits,nk-1, nlheader
    CALL write_bintab(LOSTMAPSREND, npix_fits, nk-1, header, nlheader,modefile)
    !CALL write_bintab(heal_map, npix_fits, nmaps, header, nlheader,TRIM(ADJUSTL(fake_file)))
    !CALL output_map(heal_map(:,1),header,TRIM(ADJUSTL(fake_file)))
    write(0,*)'Done'

    DEALLOCATE(map_cut,matb,heal_map)

  ENDIF
    RETURN
  END SUBROUTINE make_fake_mode_map

END MODULE Topology_map_mod_nel
