PROGRAM Topology_Lmarg
  !Program to calculate th likelihood of a topology model
  !wrt the WMAP data with various cuts.
  !Also makes realizations of the topology
  !
  !Preliminary version          Carlo 14Jun04 CITA
  !Further development          Dmitri June-Aug 2004 UofA
  !
  USE nrtype
  USE Topology_types
  USE Topology_map_mod
  USE Topology_Lmarg_mod
  USE lm_rotate, ONLY : getcplm
  USE PIX_TOOLS
  IMPLICIT NONE

  INTEGER :: iargc
  CHARACTER(len=100) :: parfile
  LOGICAL :: found
 
  EXTERNAL iargc
 
  real(DP), allocatable, dimension(:,:) :: pixels   ! For CTpp
  real(DP) :: ampl_best, ampl_var, ampl_curv, LnL_max, alpha, beta, gamma
  integer :: iter
  CHARACTER(LEN=80) :: modefile

!  character(len=100) :: infile

  IF (iargc() == 0) STOP 'Usage: topmarg <input file>'

! --------------------------------------------------------------------  
! Processing configuration file
  CALL getarg(1,parfile)
  WRITE(0,'(a13,a)') 'Input file : ', TRIM(parfile)
  found = .FALSE.
  INQUIRE(file=TRIM(parfile),exist=found)
  IF (found) THEN
     CALL make_namelist(parfile, 'toplike', './tmp.ini')
     OPEN(2,file='./tmp.ini',STATUS='OLD')
     READ(2,NML=toplike)
     CLOSE(2,STATUS='DELETE')
  ELSE
     WRITE(0,*) 'input file not found : ', TRIM(parfile)
     STOP
  ENDIF

  INQUIRE(file=TRIM(wmap_mask_file),exist=found)
  IF (found) THEN
     do_mask = .TRUE.
     WRITE(0,*) 'Using a mask'
  ELSE
     do_mask = .FALSE.
     WRITE(0,*) 'Not using a mask'
  ENDIF

  INQUIRE(file=TRIM(wmap_noise_file),exist=found)
  IF (found) THEN
     add_noise = .TRUE.
     WRITE(0,*) 'Using the noise'
  ELSE
     add_noise = .FALSE.
     add_map_noise = .FALSE. !Cannot add noise if no noise file
     WRITE(0,*) 'Not using the noise'
  ENDIF

  INQUIRE(file=TRIM(beam_file),exist=found)
  IF (found) THEN
     do_smooth = .TRUE.
  ELSE
     do_smooth = .FALSE.
     WRITE(0,*) 'Not smoothing'
  ENDIF

  IF (epsil == 0.0) THEN
     WRITE(0,*) 'Not using regularization option'
  ELSE
     WRITE(0,*) 'Using regularization option'
  ENDIF

  IF (SVD) THEN
     WRITE(0,*) 'Using singular value decomposition method'
  ELSE
     WRITE(0,*) 'Using Cholesky decomposition method'
  ENDIF
!-------------------------------------------------------------------
  write(0,'(''File with CTpp - '',$)')
  read(*,'(a)') infile
  open(102,file=TRIM(infile),status='old',form='unformatted')
  read(102) npix_fits
  write(0,*)'npix=',npix_fits
  if ( nside /= npix2nside(npix_fits) ) then
     write(0,*)'Size of Ctpp array does not match requested NSIDE',npix_fits,nside
     stop
  endif
  allocate(pixels(3,npix_fits))
  allocate(CTpp_evec(0:npix_fits-1,0:npix_fits-1))
  read(102)pixels
  read(102)CTpp_evec
  close(102)
  deallocate(pixels)
!-------------------------------------------------------------------
  write(0,'(''Map Output (make fake map)- '',$)')
  read(*,'(a)') fake_file
!-------------------------------------------------------------------

  nmaps     = 1

  CALL ReadWMAP_map()


  write(0,*)'Read the data in'
!-------------------------------------------------------------------
! Allocate main data blocks
!     CTpp_evec  - (with CTpp_eval) - full sky theoretical normalized
!                  pixel-pixel correlations in eigenvector decomposition
!     CTpp_cplm  - lm decomposition of CTpp_evec
!     CTpp       - Cut sky pixel-pixel correlation for given rotation, norm ampl
!     CNTpp      - at the end of calculations = (ampl_best*Ctpp+N)^{-1}

  allocate(CTpp_eval(0:npix_fits-1))
  allocate(CTpp(0:npix_cut-1,0:npix_cut-1))
  allocate(CNTpp(0:npix_cut-1,0:npix_cut-1))
! Decompose CTpp(_evec) into eigenfuctions stored in CTpp_evec and CTpp_eval
  CALL DECOMPOSE_AND_SAVE_EIGENVALUES()
  CALL NORMALIZE_EIGENVALUES()

! Determines the number of modes of Ctpp and keeps the number constant
!-------------------------------------------------------------------
! Main calls to determine best fit parameters
  if (do_rotate) then
     ! Decompose CTpp_evec into multipoles, stored in CTpp_cplm
     allocate(CTpp_cplm(0:lmax*(lmax+2),0:npix_fits-1))
     !allocate(CTpp_clm(0:lmax*(lmax+2),0:lmax*(lmax+2)))
     CALL Read_w8ring()
     CALL GETCPLM(CTpp_cplm,CTpp_evec,nside,lmax,w8ring)
     if (SVD) then
        CALL MODECOUNTER()
     endif

     if (find_best_angles) then
        CALL FIND_BEST_ANGLES_AND_AMPLITUDE(ampl_best,alpha,beta,gamma,LnL_max)
     else
        write(0,'(''alpha,beta,gamma - '',$)')
        read(*,*) alpha,beta,gamma
        CALL ROTATE_AND_FIND_BEST_AMPLITUDE(ampl_best,alpha,beta,gamma,LnL_max)
     endif
  else
     if (SVD) then
        CALL MODECOUNTER()
     endif
     CALL FIND_BEST_AMPLITUDE(ampl_best,LnL_max)
  endif
  ampl_best=exp(ampl_best)

!-------------------------------------------------------------------
! At the result of LnLikelihood call, CNTpp contains (CTpp+Npp)^-1
! Use it together with CTpp to calculate Fisher matrix at the bestfit amplitude

  ampl_var =LmaxFisher()
  ampl_curv=LmaxCurv()-ampl_var
  ampl_var =1.d0/sqrt(ampl_var)
  ampl_curv=1.d0/sqrt(ampl_curv)
 
  write(0,*) ampl_curv

  if (nice_output) then
     WRITE(0,'(a, 1pd15.7)') '-Ln(L) max    : ', LnL_max
     WRITE(0,'(a, 1pd15.7)') '-Ln(L) marg(F): ', LnL_max-log(ampl_var)
     WRITE(0,'(a, 1pd15.7)') '-Ln(L) marg(C): ', LnL_max-log(ampl_curv)
     WRITE(0,'(a, 1pd15.7)') ' Ampl  best   : ', ampl_best
     WRITE(0,'(a, 1pd15.7)') ' Ampl  var(F) : ', ampl_var
     WRITE(0,'(a, 1pd15.7)') ' Ampl  var(C) : ', ampl_curv
     WRITE(0,'(a, 3(1x,d12.4))') ' Angles best  :', alpha,beta,gamma
  else
     WRITE(0,'(6(1x,d15.7),3(1x,d12.4))')                             &
                 LnL_max,LnL_max-log(ampl_var),LnL_max-log(ampl_curv),&
                 ampl_best,ampl_var,ampl_curv,                        &
                 alpha,beta,gamma
  endif

  !optional output of cut-sky realization from CTpp
  Call make_fake_mode_map(ampl_best)
  IF (make_map) then
     CALL make_fake_map(ampl_best)
  endif

END PROGRAM Topology_Lmarg
