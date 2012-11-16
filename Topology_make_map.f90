PROGRAM Topology_make_map
  !Program to make random map realization from CTpp, and optional amplitude
  !Also optional is mask, rotation of CTpp and smoothing
  !
  USE Topology_types
  USE Topology_map_mod
  USE ctpp_eigen_mod
  USE healpix_extras, ONLY : Read_w8ring, ring2pixw8
  USE beams,          ONLY : collect_beams, smooth_ctpp_lm
  USE lm_rotate,      ONLY : getcplm, rotate_ctpp
  USE PIX_TOOLS
  IMPLICIT NONE

  LOGICAL :: found
 
  real(DP)     :: ampl, ang(3)
  real(DP),    allocatable, dimension(:,:)   :: map

!------------------------------------------------------------------------
!  Input Parameters
!------------------------------------------------------------------------
! Read in files (even if not being used) 

! data and noise file
  read(*,'(a)') map_signal_file
! Map modification files
  read(*,'(a)') map_noise_file
  read(*,'(a)') beam_file
! Ring Weights file
  read(*,'(a)') w8_file
! CTpp file
  read(*,'(a)') infile
! fake map output file
  read(*,'(a)') fake_file

  read(*,*) add_map_noise
  read(*,*) iseed

  read(*,*) beam_fwhm
  read(*,*) lmax

  read(*,*) ampl
  read(*,*) do_rotate
  read(*,*) ang(1),ang(2),ang(3)

!======================================================================

  IF (beam_fwhm == 0.0) THEN
     do_Gsmooth=.FALSE.
  ELSE
     do_Gsmooth=.TRUE.
  ENDIF

  INQUIRE(file=TRIM(map_signal_file),exist=found)
  WRITE(0,*) 'Signal file', TRIM(map_signal_file)
  IF (.NOT.found) THEN
     WRITE(0,*) 'Can not find file', TRIM(map_signal_file)
     STOP "No signal file"
  ENDIF

  INQUIRE(file=TRIM(infile),exist=found)
  IF (found) THEN
     WRITE(0,*) 'CTpp:', TRIM(infile)
  ELSE
     WRITE(0,*) 'Can not find file', TRIM(infile)
     STOP "No CTpp file"
  ENDIF

  INQUIRE(file=TRIM(beam_file),exist=do_expsmooth)
  WRITE(0,*) 'CTpp smoothed by:'
  WRITE(0,*) '     - Pixel window:' 
  IF (do_expsmooth .and. do_Gsmooth) THEN
     WRITE(0,*) '     - Gaussian beam (arcmin):', beam_fwhm 
     WRITE(0,*) '     - experimantal beam:', TRIM(beam_file)
  ELSEIF (do_expsmooth .and. .not.do_Gsmooth) THEN
     WRITE(0,*) '     - experimantal beam:', TRIM(beam_file)
  ELSEIF (.not.do_expsmooth .and. do_Gsmooth) THEN
     WRITE(0,*) '     - Gaussian beam (arcmin):', beam_fwhm 
  ENDIF

  INQUIRE(file=TRIM(map_noise_file),exist=found)
  IF (found) THEN
     WRITE(0,*) 'Using a noise file', TRIM(map_noise_file)
  ELSE
     WRITE(0,*) 'No separate noise file'  ! Some formats (WMAP) do not require 
                                 ! separate noise file, so add_map_noise
                                 ! will decide if any of the noise is used
  ENDIF

  IF (add_map_noise) THEN
     WRITE(0,*) 'Adding noise to the map'
  ELSE
     WRITE(0,*) 'Not adding noise'
  ENDIF

!-------------------------------------------------------------------
! Read full sky CTpp in
!
  open(102,file=TRIM(infile),status='old',form='unformatted')
  read(102) npix_fits
  nside = npix2nside(npix_fits)
  write(0,'(2(a6,I6))')'nside=',nside,' npix=',npix_fits
  if ( nside == -1 ) stop 'Size of Ctpp array does not match any nside'

!-------------------------------------------------------------------
! Allocate global working array for full sky manipulations
  allocate(FullSkyWorkSpace(0:npix_fits-1,0:npix_fits-1))
  CTpp_full => FullSkyWorkSpace
  read(102)CTpp_full
  close(102)

! Set experimental beam and pixel window to smooth CTpp
  allocate ( Wl(0:lmax,1:1) )
  if (do_expsmooth) then
     CALL collect_beams(Wl,lmax,beamfile=beam_file,nside=nside,reset=.true.)
  else
     CALL collect_beams(Wl,lmax,nside=nside,reset=.true.)
  endif
  CALL smooth_ctpp_lm(CTpp_full,lmax,window=Wl)
  
!-------------------------------------------------------------------
!     CTpp_eval  - (with CTpp_evec) - full sky theoretical normalized
!                  pixel-pixel correlations in eigenvector decomposition

  call Read_w8ring(nside,w8ring,w8_file)
  call ring2pixw8(w8ring,w8pix)
  allocate(CTpp_eval(0:npix_fits-1))

! Decompose CTpp_full into eigenfuctions stored in CTpp_evec and CTpp_eval
! CTpp_full is destroyed and disassociated in favour of CTpp_evec
  write(0,*)'Eigenvalue decomposition started'
  CALL DECOMPOSE_AND_SAVE_EIGENVALUES()
  CALL SORT_AND_LIMIT_EIGENVALUES()
  write(0,*)'Eigenvalue decomposition completed'

!-------------------------------------------------------------------
! Optionally rotate full sky eigevectors into revised CTpp_evec
  if (do_rotate) then
     ! Decompose CTpp_evec into multipoles, stored in CTpp_cplm
     CTpp_evec => FullSkyWorkSpace
     allocate(CTpp_cplm(0:lmax*(lmax+2),0:n_evalues-1))
     call getcplm(CTpp_cplm,CTpp_evec,nside,n_evalues,lmax,w8ring)
     call rotate_ctpp(CTpp_evec,CTpp_cplm,nside,n_evalues,lmax,ang(1),ang(2),ang(3),.TRUE.)
     deallocate(CTpp_cplm)
     write(0,'((a16),3(1x,e11.4))')'CTpp rotated to',ang
  endif

! make map.
  call make_fake_map(ampl,map)
  call Write_map(map)

END PROGRAM Topology_make_map
