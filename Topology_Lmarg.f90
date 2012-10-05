PROGRAM Topology_Lmarg
  !Program to calculate th likelihood of a topology model
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
  USE lm_rotate, ONLY : getcplm
  USE PIX_TOOLS
  IMPLICIT NONE

  LOGICAL :: found, do_nice_out_file
 
  real(DP) :: ampl_best, ampl_var, ampl_curv, LnL_max, alpha, beta, gamma
  !Dreal(DP) :: amp, lnamp !NELSON LOOP
  real(DP), allocatable, dimension(:,:) :: pixels
  CHARACTER(LEN=120) :: nice_out_file

!------------------------------------------------------------------------
!  Input Parameters for likelihood run
!------------------------------------------------------------------------
! Read in files (even if not being used) 
  read(*,'(a)') nice_out_file

! data and noise file
  read(*,'(a)') map_signal_file
! Map modification files
  read(*,'(a)') map_mask_file
  read(*,'(a)') beam_file
! Ring Weights file
  read(*,'(a)') w8_file
! CTpp file
  read(*,'(a)') infile
! Makefake map output file
  read(*,'(a)') fake_file

  read(*,*) do_nice_out_file
! Read in parameters
  read(*,*) nside
  read(*,*) nsh
  read(*,*) OmegaL
  read(*,*) H0
  read(*,*) Ok

  read(*,*) make_map
  read(*,*) add_map_noise
  read(*,*) iseed
  read(*,*) make_map_only

  read(*,*) beam_fwhm
  
  read(*,*) lmax

  read(*,*) do_rotate
  read(*,*) alpha,beta,gamma
  read(*,*) find_best_angles

  read(*,*) add_noise
  read(*,*) epsil

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
  IF (do_Gsmooth) THEN
     WRITE(0,*) 'Signal map smoothed by:'
     WRITE(0,*) '     - Gaussian beam (arcmin):', beam_fwhm 
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
  WRITE(0,*) '     - Pixle window:', nside 
  IF (do_expsmooth .and. do_Gsmooth) THEN
     WRITE(0,*) '     - Gaussian beam (arcmin):', beam_fwhm 
     WRITE(0,*) '     - experimantal beam:', TRIM(beam_file)
  ELSEIF (do_expsmooth .and. .not.do_Gsmooth) THEN
     WRITE(0,*) '     - experimantal beam:', TRIM(beam_file)
  ELSEIF (.not.do_expsmooth .and. do_Gsmooth) THEN
     WRITE(0,*) '     - Gaussian beam (arcmin):', beam_fwhm 
  ENDIF

  INQUIRE(file=TRIM(map_mask_file),exist=found)
  IF (found) THEN
     do_mask = .TRUE.
     WRITE(0,*) 'Using a mask', TRIM(map_mask_file)
  ELSE
     do_mask = .FALSE.
     WRITE(0,*) 'Not using a mask'
  ENDIF

  IF (add_noise) THEN
     IF (do_Gsmooth .or. do_expsmooth) THEN
        WRITE(0,*) 'Using full smoothed noise matrix'
     ELSE
        WRITE(0,*) 'Using only diagonal noise'
     ENDIF
  ELSE
     add_map_noise = .FALSE. !Cannot add noise if no noise file
     WRITE(0,*) 'Not using noise'
  ENDIF

  IF (epsil == 0.0) THEN
     WRITE(0,*) 'Not using regularization option'
  ELSE
     WRITE(0,*) 'Using regularization option, epsil =', epsil
  ENDIF

  IF (make_map_only) THEN
     WRITE(0,*)'Making map only'
     GO TO 990
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
    WRITE(103,'(1Xa,I4)') '     - Pixle window:', nside 
    IF (do_expsmooth .and. do_Gsmooth) THEN
       WRITE(103,'(1Xa,F9.4)') '     - Gaussian beam (arcmin):', beam_fwhm 
       WRITE(103,'(1Xa,a)') '     - experimantal beam:', TRIM(beam_file)
    ELSEIF (do_expsmooth .and. .not.do_Gsmooth) THEN
       WRITE(103,'(1Xa,a)') '     - experimantal beam:', TRIM(beam_file)
    ELSEIF (.not.do_expsmooth .and. do_Gsmooth) THEN
       WRITE(103,'(1Xa,F9.4)') '     - Gaussian beam (arcmin):', beam_fwhm 
    ENDIF
    IF (do_Gsmooth) THEN
       WRITE(103,'(1Xa)') 'Signal map smoothed by:'
       WRITE(103,'(1Xa,F9.4)') '     - Gaussian beam (arcmin):', beam_fwhm 
    ENDIF
  
    IF (do_mask) THEN
       WRITE(103,'(1Xa)') 'Using a mask'
    ELSE
       WRITE(103,'(1Xa)') 'Not using a mask'
    ENDIF
    IF (add_noise) THEN
      IF (do_Gsmooth .or. do_expsmooth) THEN
         WRITE(103,'(1Xa)') 'Using full noise matrix'
      ELSE
         WRITE(103,'(1Xa)') 'Using only diagonal of noise matrix'
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
       WRITE(103,'(1Xa,1XL1)')'do_rotate :', do_rotate
       IF(find_best_angles) THEN
          WRITE(103,'(1Xa,1XL1)')'find_best_angles :', find_best_angles
       ELSE
          WRITE(103,'(1Xa, 2XE11.5E2, 2XE11.5E2, 2XE11.5E2)')'Rotating to :',alpha,beta,gamma
       ENDIF
       WRITE(103,'(1Xa,I4)')'lmax :', lmax
    ENDIF
    IF(make_map) THEN
       WRITE(103,'(1Xa,1Xa)')'Making map    :', TRIM(fake_file)
       WRITE(103,'(1Xa,1XL1)')'add_map_noise :', add_map_noise
       WRITE(103,'(1Xa,1XL1)')'smooth_map    :', do_Gsmooth
    ENDIF
  ENDIF
990 CONTINUE
  !IF (SVD) THEN
 !    WRITE(0,*) 'Using singular value decomposition method'
 ! ELSE
 !    WRITE(0,*) 'Using Cholesky decomposition method'
 ! ENDIF
!-------------------------------------------------------------------
! Read full sky CTpp in, put it temporarily into CTpp_evec
!
  open(102,file=TRIM(infile),status='old',form='unformatted')
  read(102) npix_fits
  write(0,*)'npix=',npix_fits
  if ( nside /= npix2nside(npix_fits) ) then
     write(0,*)'Size of Ctpp array does not match requested NSIDE',npix_fits,nside
     stop
  endif
  allocate(CTpp_evec(0:npix_fits-1,0:npix_fits-1))
  read(102)CTpp_evec
  close(102)

!-------------------------------------------------------------------
! Read data:  signal map map_signal, noise  map_npp and mask
!             
! signal and noise are smoothed which must coincide with smoothing
! of CTpp. Calling shell script should check for that.
!

  CALL Read_w8ring(nside,w8ring,w8_file)
  CALL ring2pixw8(w8ring,w8pix)
  allocate ( Wl(0:lmax,1:1) )
  if (do_Gsmooth) then
     CALL collect_beams(Wl,lmax,G_fwhm=beam_fwhm,reset=.true.)
  else
     CALL collect_beams(Wl,lmax,reset=.true.)
  endif
  CALL ReadWMAP_map()
  write(0,*)'Read the data in'

  !GOTO 991 ! Go straight to make map
!-------------------------------------------------------------------
! Allocate main data blocks
!     CTpp_evec  - (with CTpp_eval) - full sky theoretical normalized
!                  pixel-pixel correlations in eigenvector decomposition
!     CTpp_cplm  - lm decomposition of CTpp_evec truncated to largest evalues
!     CTpp       - Cut sky pixel-pixel correlation for given rotation, norm ampl
!     CNTpp      - at the end of calculations = (ampl_best*Ctpp+N)^{-1}

  allocate(CTpp_eval(0:npix_fits-1))
  allocate(CTpp(0:npix_cut-1,0:npix_cut-1))
  allocate(CNTpp(0:npix_cut-1,0:npix_cut-1))

! Add experimental beam and pixel window to preset Gaussian and smooth CTpp
  if (do_expsmooth) then
     CALL collect_beams(Wl,lmax,beamfile=beam_file,nside=nside,reset=.false.)
  else
     CALL collect_beams(Wl,lmax,nside=nside,reset=.false.)
  endif

  CALL smooth_ctpp_lm(CTpp_evec,lmax,window=Wl)
! Decompose CTpp(_evec) into eigenfuctions stored in CTpp_evec and CTpp_eval
  CALL DECOMPOSE_AND_SAVE_EIGENVALUES()
  CALL SORT_AND_LIMIT_EIGENVALUES()
  CALL NORMALIZE_EIGENVALUES(CTpp_eval)

!-------------------------------------------------------------------
! Main calls to determine best fit parameters
  IF (make_map_only) THEN 
     call RECONSTRUCT_FROM_EIGENVALUES()
     ampl_best= 1.0d0
     GO TO 991
  ENDIF

  if (do_rotate) then
     ! Decompose CTpp_evec into multipoles, stored in CTpp_cplm
     allocate(CTpp_cplm(0:lmax*(lmax+2),0:n_evalues-1))
     call getcplm(CTpp_cplm,CTpp_evec,nside,n_evalues,lmax,w8ring)

     if (find_best_angles) then
        CALL FIND_BEST_ANGLES_AND_AMPLITUDE(ampl_best,alpha,beta,gamma,LnL_max)
     else
        CALL ROTATE_AND_FIND_BEST_AMPLITUDE(ampl_best,alpha,beta,gamma,LnL_max)
     endif
  else
     CALL FIND_BEST_AMPLITUDE(ampl_best,LnL_max)
  endif
  ampl_best=exp(ampl_best)

!-------------------------------------------------------------------
! At the result of LnLikelihood call, CNTpp contains (CTpp+Npp)^-1
! Use it together with CTpp to calculate Fisher matrix at the bestfit amplitude

  ampl_var =LmaxFisher()
  ampl_curv=LmaxCurv()
  write(0,*) 'LmaxFisher',ampl_var
  write(0,*) 'LmaxCurv_part',ampl_curv
  ampl_curv=ampl_curv-ampl_var
  ampl_var =1.d0/sqrt(ampl_var)
  ampl_curv=1.d0/sqrt(ampl_curv)
 
!  write(0,*) ampl_curv

! Debug output, will be written at the end of all convergence steps
  WRITE(0,'(a, 1pd15.7)') '-Ln(L) max    : ', LnL_max
  WRITE(0,'(a, 1pd15.7)') '-Ln(L) marg(F): ', LnL_max-log(ampl_var)
  WRITE(0,'(a, 1pd15.7)') '-Ln(L) marg(C): ', LnL_max-log(ampl_curv)
  WRITE(0,'(a, 1pd15.7)') ' Ampl  best   : ', ampl_best
  WRITE(0,'(a, 1pd15.7)') ' Ampl  var(F) : ', ampl_var
  WRITE(0,'(a, 1pd15.7)') ' Ampl  var(C) : ', ampl_curv
  WRITE(0,'(a, 3(1x,d12.4))') ' Angles best  :', alpha,beta,gamma
 
! Final one line answer to the standard output
  WRITE(*,'(f7.4,6(1x,d15.7),3(1x,d12.4))')                        &
        Ok,                                                  &
        LnL_max,LnL_max-log(ampl_var),LnL_max-log(ampl_curv),&
        ampl_best,ampl_var,ampl_curv,                        &
        alpha,beta,gamma

! Archive for storage in the nice commented file
  if (do_nice_out_file) then
     WRITE(103,'(1Xa,I)') 'npix_cut      : ', npix_cut
     WRITE(103,'(a, 1pd15.7)') '-Ln(L) max    : ', LnL_max
     WRITE(103,'(a, 1pd15.7)') '-Ln(L) marg(F): ', LnL_max-log(ampl_var)
     WRITE(103,'(a, 1pd15.7)') '-Ln(L) marg(C): ', LnL_max-log(ampl_curv)
     WRITE(103,'(a, 1pd15.7)') ' Ampl  best   : ', ampl_best
     WRITE(103,'(a, 1pd15.7)') ' Ampl  var(F) : ', ampl_var
     WRITE(103,'(a, 1pd15.7)') ' Ampl  var(C) : ', ampl_curv
  endif

! Optional output of cut-sky realization from CTpp
!  Call make_fake_mode_map(ampl_best)

991 CONTINUE
  IF (make_map) then
     CALL make_fake_map(ampl_best)
     CALL WriteWMAP_map()
  endif
  close(103)

END PROGRAM Topology_Lmarg
