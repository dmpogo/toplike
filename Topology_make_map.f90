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
  USE ct_io
  USE PIX_TOOLS
  IMPLICIT NONE

  LOGICAL :: found

  real(DP)     :: ampl, ang(3), CTpp_norm
  real(DP),    allocatable, dimension(:,:)   :: map

!------------------------------------------------------------------------
!  Input Parameters
!------------------------------------------------------------------------
! Read in files (even if not being used) 

! data and noise file
  read(*,'(a)') expdata_format
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

  read(*,*) add_noise_diag
  read(*,*) add_noise_cov
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
                                 ! separate noise file, so add_noise
                                 ! will decide if any of the noise is used
  ENDIF

  IF ( add_noise_diag .or. add_noise_cov ) THEN
     add_noise = .true.
     WRITE(0,*) 'Adding noise to the map, diag=',add_noise_diag,' cov=',add_noise_cov
  ELSE
     add_noise = .false.
     WRITE(0,*) 'Not adding noise'
  ENDIF

  IF ( add_noise_diag .and. add_noise_cov ) THEN
     WRITE(0,*) 'Warning: both diag and covariance noise are defined, use diag'
     add_noise_cov = .false.
  ENDIF

!-------------------------------------------------------------------
! Read full sky CTpp in while allocating global working array
!
  call ReadCTpp(infile,FullSkyWorkSpace,npix_fits,npol,overwrite=.true.)
  nside = npix2nside(npix_fits)
  write(0,'(2(a6,I6))')'nside=',nside,' npix=',npix_fits
  if ( nside == -1 ) stop 'Size of Ctpp array does not match any nside'
  ntot=npix_fits*npol
  CTpp_full => FullSkyWorkSpace
  CTpp_full =  ampl*CTpp_full

! Set experimental beam and pixel window to smooth CTpp
  call Read_w8ring(nside,w8ring,w8_file)
  call ring2pixw8(w8ring,w8pix)
  allocate ( Wl(0:lmax,1:npol) )
  if (do_expsmooth) then
     CALL collect_beams(Wl,lmax,beamfile=beam_file,nside=nside,reset=.true.)
  else
     CALL collect_beams(Wl,lmax,nside=nside,reset=.true.)
  endif
  CALL smooth_ctpp_lm(CTpp_full,npix_fits,npol,lmax,window=Wl)
  
! Get map_npp, ensure it is full sky
  if ( add_noise ) then
     do_Gsmooth=.false.
     do_mask = .false.
     CALL ReadExpData(TRIM(ADJUSTL(expdata_format)))
     if ( npix_cut /= ntot ) then
        stop 'Ensure that noise is not masked'
     endif
     write(0,*)'Read the noise in'
     CTpp_full = CTpp_full+map_npp
  endif

!-------------------------------------------------------------------
!     CTpp_eval  - (with CTpp_evec) - full sky theoretical normalized
!                  pixel-pixel correlations in eigenvector decomposition

  allocate(CTpp_eval(0:ntot-1))

! Decompose CTpp_full into eigenfuctions stored in CTpp_evec and CTpp_eval
! CTpp_full is destroyed and disassociated in favour of CTpp_evec
  write(0,*)'Eigenvalue decomposition started'
  CALL DECOMPOSE_AND_SAVE_EIGENVALUES()
  CALL SORT_AND_LIMIT_EIGENVALUES()
  if ( ampl == 1.0_dp ) then
     CALL NORMALIZE_EIGENVALUES(CTpp_norm)
     write(0,*)'Eigenvalue decomposition completed, scaled by ',CTpp_norm
  endif

!  fake_file='ttt_testmap_000'
!  call make_fake_map(map)
!  call Write_map(map)

!-------------------------------------------------------------------
! Optionally rotate full sky eigevectors into revised CTpp_evec
  if (do_rotate) then
     ! Decompose CTpp_evec into multipoles, stored in CTpp_cplm
     CTpp_evec => FullSkyWorkSpace
     allocate(CTpp_cplm(1:npol,0:lmax*(lmax+2),0:n_evalues-1))
     call getcplm(CTpp_cplm,CTpp_evec,nside,n_evalues,npol,lmax,w8ring)
     call rotate_ctpp(CTpp_evec,CTpp_cplm,nside,n_evalues,npol,lmax,ang(1),ang(2),ang(3),.TRUE.)
     deallocate(CTpp_cplm)
     write(0,'((a16),3(1x,e11.4))')'CTpp rotated to',ang
  endif

! make map.
!  fake_file='ttt_testmap_turn'
  call make_fake_map(map)
  call Write_map(map)

END PROGRAM Topology_make_map
