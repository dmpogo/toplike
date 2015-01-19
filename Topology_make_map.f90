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

  real(DP)     :: ampl, ang(3), CTpp_norm, noise_scale
  real(DP),    allocatable, dimension(:,:)   :: map
  character(255)                             :: map_output_file

!------------------------------------------------------------------------
!  Input Parameters
!------------------------------------------------------------------------
! Read in files (even if not being used) 

! Map modification files
  read(*,'(a)') map_noise_file
! Ring Weights file
  read(*,'(a)') w8_file
! CTpp file
  read(*,'(a)') infile
  read(*,'(a)') ctpp_beam_file
  read(*,*)     ctpp_beam_fwhm
! fake map output file
  read(*,'(a)') map_output_file

  read(*,'(a)') map_beam_file
  read(*,*)     map_beam_fwhm

  read(*,*) add_noise_cov
  read(*,*) noise_scale
  read(*,*) iseed

  read(*,*) lmax

  read(*,*) ampl
  read(*,*) do_rotate
  read(*,*) ang(1),ang(2),ang(3)

!======================================================================

  INQUIRE(file=TRIM(infile),exist=found)
  IF (found) THEN
     WRITE(0,*) 'CTpp:', TRIM(infile)
  ELSE
     WRITE(0,*) 'Can not find file', TRIM(infile)
     STOP "No CTpp file"
  ENDIF

  INQUIRE(file=TRIM(ctpp_beam_file),exist=do_smooth_ctpp)
  WRITE(0,*) 'CTpp smoothed by:'
  WRITE(0,*) '     - Pixel window:' 
  WRITE(0,*) 'Before adding noise, CTpp is also smoothed with:'
  IF (do_smooth_ctpp) THEN
     WRITE(0,*) '     - experimental beam:', TRIM(ctpp_beam_file)
  ENDIF
  IF (ctpp_beam_fwhm > 0.0_dp) THEN
     do_smooth_ctpp=.true.
     WRITE(0,*) '     - Gaussian beam (arcmin):', ctpp_beam_fwhm 
  ENDIF

  INQUIRE(file=TRIM(map_beam_file),exist=do_smooth_data)
  WRITE(0,*) 'Output map smoothed by:'
  IF (do_smooth_data) THEN
     WRITE(0,*) '     - experimental beam:', TRIM(map_beam_file)
  ENDIF
  IF (map_beam_fwhm > 0.0_dp) THEN
     do_smooth_data=.true.
     WRITE(0,*) '     - Gaussian beam (arcmin):', map_beam_fwhm 
  ENDIF

  INQUIRE(file=TRIM(map_noise_file),exist=found)
  IF (found) THEN
     WRITE(0,*) 'Using a noise file', TRIM(map_noise_file)
  ELSE
     WRITE(0,*) 'No separate noise file'  ! Some formats (WMAP) do not require 
                                 ! separate noise file, so add_noise
                                 ! will decide if any of the noise is used
  ENDIF

  IF ( add_noise_cov ) THEN
     add_noise = .true.
     WRITE(0,*) 'Adding noise to the map'
  ELSE
     add_noise = .false.
     WRITE(0,*) 'Not adding noise'
  ENDIF

!-------------------------------------------------------------------
! Read full sky CTpp in while allocating global working array
!
  write(0,*)'reading in CTpp ==========================================='

  call ReadCTpp(infile,FullSkyWorkSpace,npix_fits,npol,overwrite=.true.)
  nside = npix2nside(npix_fits)
  write(0,'(2(a6,I6))')'nside=',nside,' npix=',npix_fits
  if ( nside == -1 ) stop 'Size of Ctpp array does not match any nside'
  ntot=npix_fits*npol
  CTpp_full => FullSkyWorkSpace

! Set experimental beam and pixel window to smooth CTpp
  write(0,*)'smoothing and normalizing CTpp ============================'
  call Read_w8ring(nside,w8ring,w8_file)
  call ring2pixw8(w8ring,w8pix)
  allocate ( Wl(0:lmax,1:npol) )
  CALL collect_beams(Wl,lmax,nside=nside,reset=.true.)
  if (do_smooth_ctpp) then
     CALL collect_beams(Wl,lmax,beamfile=ctpp_beam_file,G_fwhm=ctpp_beam_fwhm,reset=.false.)
  endif
  CALL smooth_ctpp_lm(CTpp_full,npix_fits,npol,lmax,window=Wl)
  CALL normalize_ctpp(CTpp_norm)

  CTpp_full =  ampl*CTpp_full


!-------------------------------------------------------------------
! Read in noise convariance. Defines map_npp

   if ( add_noise ) then
      ! Set up full sky parameters
      do_smooth_noise = .false.
      add_noise_diag = .false.
      do_mask = .false.
      epsil = -1.0_dp
      npix_cut = ntot

      CALL READ_CovNoise()
      map_npp=map_npp*(noise_scale**2)
      CTpp_full = CTpp_full+map_npp
      write(0,*)'Unsmoothed noise scaled by ',noise_scale,'added to CTpp'
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

!  call make_fake_map(map)
!  call Write_map(map,'ttt_testmap_600')

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
  if (do_smooth_data) then
     CALL collect_beams(Wl,lmax,beamfile=map_beam_file,G_fwhm=map_beam_fwhm,reset=.true.)
  endif
  call make_fake_map(map)
  call Write_map(map,map_output_file)

CONTAINS

  subroutine normalize_ctpp(ctppnorm_out)
  real(DP), intent(out), optional  :: ctppnorm_out

  integer(I4B)            :: i, l, lmcount
  real(DP)                :: flatnorm, ctppnorm

  ! Normalize to power over l=2,lnorm being equal to that with flat curlCl
  ! Cl = 2pi/(l*(l+1)) curlCl,   4pi/Npix sum_p Cpp = \sum_l (2l+1) Cl

  flatnorm=0.0_dp
  do l=2,lnorm
     flatnorm = flatnorm + (l+0.5_dp)/((l+1.0_dp)*l)*Wl(l,1)**2
  enddo
  flatnorm = flatnorm*curlCl_in_mK

  do i=1,npix_fits
     ctppnorm = ctppnorm + CTpp_full(i,i)
  enddo

  ctppnorm=flatnorm*npix_fits/ctppnorm
  CTpp_full = CTpp_full*ctppnorm

  write(0,*)'Normalized over l=',lnorm,'to curlCl(mK)=',curlCl_in_mK
  write(0,*)'Normalization factors, flatnorm=',flatnorm,'ctppnorm=',ctppnorm

  if (present(ctppnorm_out)) ctppnorm_out=ctppnorm
  end subroutine normalize_ctpp


END PROGRAM Topology_make_map
