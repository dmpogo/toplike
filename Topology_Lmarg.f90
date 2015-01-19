PROGRAM Topology_Lmarg
  !Program to calculate the likelihood of a topology model
  !wrt the WMAP data with various cuts.
  !Also makes realizations of the topology
  !
  !Preliminary version          CC 14Jun04 CITA
  !Further development          DP June-Aug 2004 UofA
  !
  USE Topology_types
  USE Topology_map_mod
  USE Topology_Lmarg_mod
  USE healpix_extras, ONLY : Read_w8ring, ring2pixw8
  USE beams,          ONLY : collect_beams, smooth_ctpp_lm
  USE lm_rotate,      ONLY : getcplm
  USE ct_io
  USE PIX_TOOLS
  IMPLICIT NONE

  LOGICAL :: found, do_nice_out_file, random_angles
 
  type(LikelihoodData) :: LnL_max
  real(DP)             :: ampl_var, ampl_curv, CTpp_norm
  real(DP), allocatable, dimension(:,:) :: pixels
  integer(I4B)                          :: iter, niter, seed_size
  CHARACTER(LEN=255) :: nice_out_file, inputang

!------------------------------------------------------------------------
!  Input Parameters for likelihood run
!------------------------------------------------------------------------
! Read in files (even if not being used) 
  read(*,'(a)') nice_out_file

! data and noise file
  read(*,'(a)') expdata_format
  read(*,*)     expdata_scale
  read(*,'(a)') map_signal_file
  read(*,'(a)') map_noise_file
! Map modification files
  read(*,'(a)') map_mask_file
  read(*,'(a)') map_beam_file
! Ring Weights file
  read(*,'(a)') w8_file
! Fiducial CTpp
  read(*,'(a)') fidfile
! CTpp file
  read(*,'(a)') infile
! Theory smoothing
  read(*,'(a)') ctpp_beam_file

  read(*,*) do_nice_out_file
! Read in parameters
  read(*,*) nside
  read(*,*) npol
  read(*,*) lmax
! Next four parameters are strictly for printout
  read(*,*) nsh
  read(*,*) OmegaL
  read(*,*) H0
  read(*,*) Ok

  read(*,*) map_beam_fwhm
  read(*,*) ctpp_beam_fwhm

  read(*,*) do_rotate
  read(*,'(a)') inputang
  read(*,*) find_best_angles
  read(*,*) niter
  read(*,*) iseed

  read(*,*) add_noise_diag
  read(*,*) add_noise_cov
  read(*,*) epsil

!======================================================================

  call random_seed(seed_size)
  call random_seed(put=iseed+37*(/ (iter-1, iter=1,seed_size) /))

  INQUIRE(file=TRIM(map_signal_file),exist=found)
  WRITE(0,*) 'Signal:'
  WRITE(0,'(5x,255a)') TRIM(ADJUSTL(map_signal_file))
  IF (.NOT.found) THEN
     WRITE(0,*) 'Can not find file', TRIM(map_signal_file)
     STOP "No signal file"
  ENDIF

  INQUIRE(file=TRIM(map_beam_file),exist=do_smooth_data)
  IF (do_smooth_data) THEN
     WRITE(0,*) 'Signal map smoothed by:',TRIM(map_beam_file)
  ENDIF
  IF (map_beam_fwhm > 0.0_dp) THEN
     do_smooth_data = .TRUE.
     WRITE(0,*) 'Signal map smoothed by Gaussian beam:', map_beam_fwhm
  ENDIF



  INQUIRE(file=TRIM(fidfile),exist=found)
  IF (found) THEN
     WRITE(0,*) 'fiducial CTpp:'
     WRITE(0,'(5x,255a)') TRIM(ADJUSTL(fidfile))
  ELSE
     WRITE(0,*) 'Can not find file', TRIM(fidfile)
     STOP "No fiducial CTpp file"
  ENDIF

  INQUIRE(file=TRIM(infile),exist=found)
  IF (found) THEN
     WRITE(0,*) 'CTpp:'
     WRITE(0,'(5x,255a)') TRIM(ADJUSTL(infile))
  ELSE
     WRITE(0,*) 'Can not find file', TRIM(infile)
     STOP "No CTpp file"
  ENDIF

  INQUIRE(file=TRIM(ctpp_beam_file),exist=do_smooth_ctpp)
  WRITE(0,*) 'CTpp smoothed by:'
  WRITE(0,*) '     - Pixel window:', nside 
  IF (do_smooth_ctpp) THEN
     WRITE(0,*) '     - beam file:', TRIM(ctpp_beam_file)
  ENDIF
  IF (ctpp_beam_fwhm > 0.0_dp) THEN
     do_smooth_ctpp = .TRUE.
     WRITE(0,*) '     - Gaussian beam:', ctpp_beam_fwhm
  ENDIF

  INQUIRE(file=TRIM(map_mask_file),exist=found)
  IF (found) THEN
     do_mask = .TRUE.
     WRITE(0,*) 'Using a mask ', TRIM(map_mask_file)
  ELSE
     do_mask = .FALSE.
     WRITE(0,*) 'Not using a mask'
  ENDIF

  INQUIRE(file=TRIM(map_noise_file),exist=found)
  IF (found) THEN
     WRITE(0,*) 'Using a noise file ', TRIM(map_noise_file)
  ELSE
     WRITE(0,*) 'No noise file'  ! Some formats (WMAP) do not require 
                                 ! separate noise file, so add_noise
                                 ! will decide if any of the noise is used
     WRITE(0,*) TRIM(map_noise_file)
     IF ( add_noise_cov ) THEN
        WRITE(0,*) 'Noise via covariance cant be used'
        add_noise_cov = .false.
     ENDIF
  ENDIF

  IF ( add_noise_diag .or. add_noise_cov ) THEN
     add_noise = .true.
     do_smooth_noise = do_smooth_data
     WRITE(0,*) 'Adding noise to the map, diag=',add_noise_diag,' cov=',add_noise_cov
     WRITE(0,*) 'Smoothing noise ',do_smooth_noise
  ELSE
     add_noise = .false.
     WRITE(0,*) 'Not adding noise'
  ENDIF

  IF ( add_noise_diag .and. add_noise_cov ) THEN
     WRITE(0,*) 'Warning: both diag and covariance noise are defined, use diag'
     add_noise_cov = .false.
  ENDIF

  IF (epsil == 0.0) THEN
     WRITE(0,*) 'Not using regularization option'
  ELSE
     WRITE(0,*) 'Using regularization option, epsil =', epsil
  ENDIF

  IF (do_rotate) THEN
     if ( index(inputang,'random') /= 0 ) then
        random_angles=.true.
        LnL_max%ang=(/0.0_dp,0.0_dp,0.0_dp/)
        WRITE(0,*)'Rotating to random angles'
     else
        random_angles=.false.
        read(inputang,*) LnL_max%ang
        WRITE(0,'(1x,"Rotating to:",3(1x,d12.4))') LnL_max%ang
     endif 

     if (find_best_angles) then
        WRITE(0,*)'find_best_angles :', find_best_angles
     endif
     WRITE(0,*)'lmax :', lmax
  ELSE
     WRITE(0,*)'No rotation'
  ENDIF

  IF (do_nice_out_file) THEN
     WRITE(0,*) 'Using nice out file'
  ELSE
     WRITE(0,*) 'Not using nice out file'
  ENDIF

  IF(do_nice_out_file) THEN
    OPEN(103,file=TRIM(nice_out_file),status='Unknown')
    WRITE(103,'(1Xa,a)') 'CTpp file:', TRIM(infile)
    WRITE(103,'(1Xa,a)') 'Signal file   :', TRIM(map_signal_file)

    WRITE(103,'(1Xa,I4)')'Nside  :', nside
    WRITE(103,'(1Xa,I4)')'nsh    :', nsh
    WRITE(103,'(1Xa,F9.4)')'OmegaL :', OmegaL
    WRITE(103,'(1Xa,F9.4)')'H0     :', H0

    WRITE(103,'(1Xa)') 'CTpp smoothed by:'
    WRITE(103,'(1Xa,I4)') '     - Pixel window:', nside 
    IF (do_smooth_ctpp ) THEN
       WRITE(103,'(1Xa,a)') '    - beam:', TRIM(ctpp_beam_file)
    ENDIF

    IF (do_smooth_data) THEN
       WRITE(103,'(1Xa)') 'Signal map smoothed by:'
       WRITE(103,'(1Xa,a)') '    - beam:', TRIM(map_beam_file)
    ENDIF
  
    IF (do_mask) THEN
       WRITE(103,'(1Xa)') 'Using a mask'
    ELSE
       WRITE(103,'(1Xa)') 'Not using a mask'
    ENDIF
    IF (add_noise) THEN
      IF ( do_smooth_noise ) THEN
         WRITE(103,'(1Xa)') 'Noise smoothed by:'
         WRITE(103,'(1Xa,a)') '    - beam:', TRIM(map_beam_file)
      ENDIF
    ELSE
       WRITE(103,'(1Xa)') 'Not using noise'
    ENDIF
   
    IF (epsil == 0.0) THEN
       WRITE(103,'(1Xa)') 'Not using regularization option'
    ELSE
       WRITE(103,'(1Xa,E9.2E1)') 'Using regularization option, epsil =', epsil
    ENDIF

    IF(do_rotate) THEN

       if ( index(inputang,'random') /= 0 ) then
          WRITE(103,'(1Xa)')'Rotating to random angles'
       else
          WRITE(103,'(1Xa, 2XE11.5E2, 2XE11.5E2, 2XE11.5E2)')'Rotating to :',LnL_max%ang
       endif 

       IF(find_best_angles) THEN
          WRITE(103,'(1Xa,1XL1)')'find_best_angles :', find_best_angles
       ENDIF
       WRITE(103,'(1Xa,I4)')'lmax :', lmax
    ENDIF
  ENDIF

!-------------------------------------------------------------------
! Read and sets: 
!       mask,npix_cut,map_signal(npix_cut),diag_npp(npix_cut)
!       signal is optionally smoothed 

  write(0,*)'Reading the data in ===================================='
  CALL Read_w8ring(nside,w8ring,w8_file)
  CALL ring2pixw8(w8ring,w8pix)
  allocate ( Wl(0:lmax,1:3) )
  CALL collect_beams(Wl,lmax,reset=.true.)

  if (do_smooth_data) then
     CALL collect_beams(Wl,lmax,beamfile=map_beam_file,G_fwhm=map_beam_fwhm,reset=.false.)
  endif
  CALL ReadExpData(TRIM(ADJUSTL(expdata_format)))
 
!-------------------------------------------------------------------
! Read in and smooth noise convariance. Defines map_npp

  write(0,*)'Reading the noise in ===================================='
  if (do_smooth_noise) then
     CALL collect_beams(Wl,lmax,beamfile=map_beam_file,G_fwhm=map_beam_fwhm,reset=.true.)
  endif
  CALL READ_CovNoise()

!-------------------------------------------------------------------
! Read in fiducial model and set up cut-sky mode basis

  write(0,*)'Reading the fiducial model in ============================'
  CALL ReadCTpp(fidfile,FullSkyWorkSpace,npix_fits,npol,overwrite=.true.)
  write(0,'(2(a6,I6))')'nside=',nside,' npix=',npix_fits
  if ( nside /= npix2nside(npix_fits) ) then
     stop 'Size of fiducial Ctpp array does not match requested NSIDE'
  endif
  ntot=npix_fits*npol
  CTpp_fid => FullSkyWorkSpace

  write(0,*)'Smoothing the fiducial model and setting basis modes ====='
  CALL collect_beams(Wl,lmax,nside=nside,reset=.true.)
  if (do_smooth_ctpp) then
     CALL collect_beams(Wl,lmax,beamfile=ctpp_beam_file,G_fwhm=ctpp_beam_fwhm,reset=.false.)
  endif
  CALL smooth_ctpp_lm(CTpp_fid,npix_fits,npol,lmax,window=Wl)

  CALL SET_BASIS_MODES()
 
!-------------------------------------------------------------------
! Expand the data in the basis. Deallocates map_npp and creates CNpp

  write(0,*)'Projecting data and noise onto basis modes================='

  CALL PROJECT_VECTOR_ONTO_BASIS_MODES(map_signal)
  CALL PROJECT_NOISE_ONTO_BASIS_MODES()

!-------------------------------------------------------------------
! Read full sky CTpp in
!
  write(0,*)'reading in CTpp ==========================================='
 
  CALL ReadCTpp(infile,FullSkyWorkSpace,npix_fits,npol,overwrite=.false.)
  CTpp_full => FullSkyWorkSpace

  write(0,*)'smoothing and normalizing CTpp ============================'
  CALL collect_beams(Wl,lmax,nside=nside,reset=.true.)
  if (do_smooth_ctpp) then
     CALL collect_beams(Wl,lmax,beamfile=ctpp_beam_file,G_fwhm=ctpp_beam_fwhm,reset=.false.)
  endif
  CALL smooth_ctpp_lm(CTpp_full,npix_fits,npol,lmax,window=Wl)

  CALL normalize_ctpp(CTpp_norm)

!-------------------------------------------------------------------
! Allocate main data blocks
!     CTpp_evec  - (with CTpp_eval) - full sky theoretical normalized
!                  pixel-pixel correlations in eigenvector decomposition.
!                  eigenvectors storage CTpp_evec points to FullSkyWorkSpace
!     CTpp_cplm  - lm decomposition of CTpp_evec truncated to largest evalues
!     CTpp       - Cut sky pixel-pixel correlation for given rotation, norm ampl
!     CNTpp      - at the end of calculations = (ampl_best*Ctpp+N)^{-1}

  allocate(CTpp_eval(0:ntot-1))
  allocate(CTpp(0:nmode_cut-1,0:nmode_cut-1))
  allocate(CNTpp(0:nmode_cut-1,0:nmode_cut-1))

  write(0,*)'Decomposing CTpp in eigenmodes and Ylm  =============='
! Decompose CTpp_full into eigenfuctions stored in CTpp_evec and CTpp_eval
! CTpp_full is destroyed and disassociated in favour of CTpp_evec
  CALL DECOMPOSE_AND_SAVE_EIGENVALUES()
  CALL SORT_AND_LIMIT_EIGENVALUES()
!  CALL NORMALIZE_EIGENVALUES(CTpp_norm)
! Decompose CTpp_evec into multipoles, stored in CTpp_cplm
  allocate(CTpp_cplm(1:npol,0:lmax*(lmax+2),0:n_evalues-1))
  CALL GETCPLM(CTpp_cplm,CTpp_evec,nside,n_evalues,npol,lmax,w8ring)

  write(0,*) 'Preparation of CTpp completed ======================='
!-------------------------------------------------------------------
! Main calls to determine best fit parameters
  do iter=1,niter

  if (do_rotate) then
     write(0,*) 'Starting angles ', TRIM(inputang)
     if (find_best_angles) then
        CALL FIND_BEST_ANGLES_AND_AMPLITUDE(LnL_max,random_angles)
     else
        CALL ROTATE_AND_FIND_BEST_AMPLITUDE(LnL_max,random_angles)
     endif
  else
     CALL FIND_BEST_AMPLITUDE(LnL_max)
  endif

!-------------------------------------------------------------------
! At the result of LnLikelihood call, CNTpp contains (Abest*CTpp+Npp)^-1
! Use it together with CTpp to calculate Fisher matrix at the bestfit amplitude

  ampl_var =LmaxFisher()
  ampl_curv=LmaxCurv()
  write(0,*) 'LmaxFisher',ampl_var
  write(0,*) 'LmaxCurv_part',ampl_curv
  ampl_curv=ampl_curv-ampl_var
  ampl_var =1.d0/sqrt(ampl_var)
  ampl_curv=1.d0/sqrt(ampl_curv)
 
! Debug output, will be written at the end of all convergence steps
  WRITE(0,'(a, 1pd15.7)') '-Ln(L) max    : ', LnL_max%LnL
  WRITE(0,'(a, 1pd15.7)') '-Ln(L) marg(F): ', LnL_max%LnL-log(ampl_var)
  WRITE(0,'(a, 1pd15.7)') '-Ln(L) marg(C): ', LnL_max%LnL-log(ampl_curv)
  WRITE(0,'(a, 1pd15.7)') ' Ampl  best   : ', exp(LnL_max%ampl)
  WRITE(0,'(a, 1pd15.7)') ' Ampl  var(F) : ', ampl_var
  WRITE(0,'(a, 1pd15.7)') ' Ampl  var(C) : ', ampl_curv
  WRITE(0,'(a, 1pd15.7)') ' CTpp  norm   : ', CTpp_norm
  WRITE(0,'(a, 3(1x,d12.4))') ' Angles best  :', LnL_max%ang
 
! Final one line answer to the standard output
  WRITE(*,'(f9.4,4(1x,d15.7),3(1x,d12.4),5(1x,d15.7))')                        &
        Ok,                                                  &
        LnL_max,LnL_max%LnL-log(ampl_var),LnL_max%LnL-log(ampl_curv),&
        ampl_var,ampl_curv,CTpp_norm

  enddo

! Archive for storage in the nice commented file
  if (do_nice_out_file) then
     WRITE(103,'(1Xa,I)') 'npix_cut      : ', npix_cut
     WRITE(103,'(1Xa,I)') 'nmode_cut     : ', nmode_cut
     WRITE(103,'(a, 1pd15.7)') '-Ln(L) max    : ', LnL_max%LnL
     WRITE(103,'(a, 1pd15.7)') '-Ln(L) marg(F): ', LnL_max%LnL-log(ampl_var)
     WRITE(103,'(a, 1pd15.7)') '-Ln(L) marg(C): ', LnL_max%LnL-log(ampl_curv)
     WRITE(103,'(a, 1pd15.7)') ' Ampl  best   : ', exp(LnL_max%ampl)
     WRITE(103,'(a, 1pd15.7)') ' Ampl  var(F) : ', ampl_var
     WRITE(103,'(a, 1pd15.7)') ' Ampl  var(C) : ', ampl_curv
     WRITE(103,'(a, I)')'Normalization range l=2,',lnorm
     WRITE(103,'(a, 1pd15.7)')'  curlCl(mK)   :',curlCl_in_mK
     WRITE(103,'(a, 1pd15.7)') ' CTpp  norm   : ', CTpp_norm
     close(103)
  endif

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

END PROGRAM Topology_Lmarg
