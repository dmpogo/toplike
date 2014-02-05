PROGRAM merge_maps
  USE HEALPIX_TYPES
  USE HEAD_FITS
  USE FITSTOOLS
  USE PIX_TOOLS

  IMPLICIT NONE

  INTEGER, PARAMETER :: nlheader = 80
  CHARACTER(LEN=80), DIMENSION(1:nlheader) :: header

  CHARACTER(LEN=255) :: infile1, infile2, outfile, inputline
  CHARACTER(LEN=10)  :: sufx
  REAL(DP), allocatable, dimension(:,:)    :: map1, map2, mapout
  INTEGER(I4B)       :: nfiles=1, narg, nside, nmapsout, iostatus
  INTEGER(I4B)       :: npix1, nmaps1, ordering1, ifpol1
  INTEGER(I4B)       :: npix2, nmaps2, ordering2, ifpol2
  LOGICAL            :: filefound

    narg = iargc()
    if ( narg < 1 ) then
       write(0,*)'Usage: merge_maps infile1 [infile2 [manual]]'
       stop
    endif

    CALL getarg(1,infile1)
    WRITE(0,'(a14,a)') 'Input1  file : ', TRIM(infile1)

    if ( narg >= 2 ) then
       CALL getarg(2,infile2)
    else
       infile2=infile1
    endif
    WRITE(0,'(a14,a)') 'Input2  file : ', TRIM(infile2)

    outfile = TRIM(ADJUSTL(infile1))//'_IQU.fits'
    if ( narg <= 2 ) then
       infile1 = TRIM(ADJUSTL(infile1))//'_I.fits'
       INQUIRE(file=TRIM(ADJUSTL(infile2))//'_QU.fits',exist=filefound)
       if ( filefound ) then
          infile2 = TRIM(ADJUSTL(infile2))//'_QU.fits'
       else
          INQUIRE(file=TRIM(ADJUSTL(infile2))//'_P.fits',exist=filefound)
          if ( filefound ) then
             infile2 = TRIM(ADJUSTL(infile2))//'_P.fits'
          else
             write(0,*) 'Polarization file not found ',infile2
          endif
       endif
    endif

    npix1=getsize_fits(infile1,nmaps=nmaps1,ordering=ordering1) 
    npix2=getsize_fits(infile2,nmaps=nmaps2,ordering=ordering2) 

    write(0,*) npix1,npix2,nmaps1,nmaps2,ifpol1,ifpol2

    if (npix1 /= npix2) stop 'Input maps are of different size'

    allocate( map1(0:npix1-1,nmaps1), map2(0:npix2-1,nmaps2) )
    call input_map(TRIM(ADJUSTL(infile1)),map1,npix1,nmaps1)
    call input_map(TRIM(ADJUSTL(infile2)),map2,npix2,nmaps2)

    nside=npix2nside(npix1)

! Check ordering, convert to RING if needed
    if ( ordering1 == 0 ) then
      write(0,*)'Ordering of the input map 1 is unknown, assumed RING'
    else if ( ordering1 == 2 ) then
      write(0,*)'Input 1 converted from NESTED to RING'
      call convert_nest2ring(nside,map1)
    else
       write(0,*)'Input 1 is in RING pixelization'
    endif

    if ( ordering2 == 0 ) then
      write(0,*)'Ordering of the input map 2 is unknown, assumed RING'
    else if ( ordering2 == 2 ) then
      write(0,*)'Input 2 converted from NESTED to RING'
      call convert_nest2ring(nside,map2)
    else
       write(0,*)'Input 2 is in RING pixelization'
    endif

    nmapsout = nmaps1+nmaps2
    if ( nmapsout > 3) stop 'More than three maps are unsupported'
    allocate( mapout(0: npix1-1, nmapsout) )

    if (nmaps1 <= nmaps2) then
       mapout(:,1:nmaps1) = map1
       mapout(:,nmaps1+1:nmapsout) = map2
    else
       mapout(:,1:nmaps2) = map2
       mapout(:,nmaps2+1:nmapsout) = map1
    endif

    deallocate(map1, map2)

    CALL write_minimal_header(header,'MAP',nside=nside,ordering='RING',creator='merge_maps',coordsys='G',units='mK',polar=.true.)

    INQUIRE(file=TRIM(ADJUSTL(outfile)),exist=filefound)
    IF(filefound) THEN
       open(26,file=TRIM(ADJUSTL(outfile)),form='unformatted')
       close(26,status='delete')
       write(0,*) 'Overwriting existing fits file'
    ELSE
       open(26,file=TRIM(ADJUSTL(outfile)),form='unformatted',status='new',iostat=iostatus)
       IF ( iostatus > 0 ) THEN
          write(0,*) 'Unable to open output ',TRIM(ADJUSTL(outfile))
          stop
       ELSE
          close(26,status='delete')
       ENDIF
    ENDIF
    
    write(0,*)'Writing bintab'
    write(0,*)size(mapout,1),size(mapout,2),npix1,nmapsout,nlheader
    CALL write_bintab(mapout, npix1, nmapsout, header, nlheader,TRIM(ADJUSTL(outfile)))
    write(0,*)'Done'

END PROGRAM merge_maps
