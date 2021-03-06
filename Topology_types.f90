MODULE TOPOLOGY_TYPES
  USE HEALPIX_TYPES

  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER :: code = 'TOPMARG'
  CHARACTER(LEN=*), PARAMETER :: version = 'Winter2014'
  
! ====================================================
! Output
  TYPE LikelihoodData
       REAL(DP)  :: LnL    ! minus log Likelhood LnL = 0.5*( Det + chisq )
       REAL(DP)  :: Det    ! log of the determinant
       REAL(DP)  :: chisq  ! chi square part
       REAL(DP)  :: ampl   ! amplitude
       REAL(DP)  :: ang(3) ! orientation Euler angles
  END TYPE LikelihoodData
     
! ====================================================
! Debugging output
  LOGICAL  :: DEBUG=.false.
  LOGICAL  :: STORE_REBINNED_MAPS=.true., STOP_AFTER_STORING_REBINNED_MAPS=.false.

! ====================================================
! Control

! Eigenmode basis type:   0 - S/N, 1 - C Bweight, negative - scaled P

  INTEGER(I1B), PARAMETER :: BASIS_TYPE = 1

! Eigenmode cut strategy, for cut sky basis and full sky CTpp_full
  INTEGER(I1B), PARAMETER :: S_NONE=0, S_LCUT=1, S_NCUT=2, S_CONDITIONING=3
  TYPE EVALUE_CUT
       INTEGER(I1B)    :: STRATEGY=S_NONE
       INTEGER(I4B)    :: LMAX=40
       INTEGER(I4B)    :: NMAX=1200
       REAL(DP)        :: CONDITION_NUMBER=1.d-4
       LOGICAL         :: SET_BAD_MODES_HIGH=.false.
       REAL(DP)        :: BAD_MODES_NOISE=1.d4
  END TYPE EVALUE_CUT

  TYPE(EVALUE_CUT) :: evalue_cut_fsky=EVALUE_CUT(S_NONE,0,0,0.0,.false.,0.0)
  TYPE(EVALUE_CUT) :: evalue_cut_csky=EVALUE_CUT(S_NCUT,0,1085,0.0,.false.,0.0)
!  TYPE(EVALUE_CUT) :: evalue_cut_csky=EVALUE_CUT(S_NCUT,0,837,0.0,.false.,0.0)
!  TYPE(EVALUE_CUT) :: evalue_cut_csky=EVALUE_CUT(S_NCUT,0,1200,0.0,.false.,0.0)

! ====================================================


  REAL(DP),    PARAMETER :: Top_bad_value = 100000.d0

! ====================================================
! CTpp normalization
  integer(I4B)  :: lnorm=64
  real(DP)      :: curlCl_in_mK = 1.d-2  ! 10x real to have ampl~-2
!  real(DP)      :: curlCl_in_mK = 1.d-3  !  temporary for maps

! ====================================================
! Data
  CHARACTER(LEN=6)                           :: expdata_format='WMAP'
  REAL(DP)                                   :: expdata_scale=1.0_dp
  REAL(DP),    DIMENSION(:),   ALLOCATABLE   :: map_signal, diag_noise
  REAL(DP),    DIMENSION(:,:), ALLOCATABLE   :: map_npp
  LOGICAL,     DIMENSION(:),   ALLOCATABLE   :: map_mask

! ====================================================
! Global arrays to hold different matrices

  ! Intermittent data that operates in full sky (npix_fits, npix_fits) space
  REAL(DP),    DIMENSION(:,:), ALLOCATABLE, TARGET  :: FullSkyWorkSpace
  REAL(DP),    DIMENSION(:,:), POINTER  :: CTpp_full => NULL()
  REAL(DP),    DIMENSION(:,:), POINTER  :: CTpp_evec => NULL()
  REAL(DP),    DIMENSION(:,:), POINTER  :: CTpp_fid  => NULL()

  ! Persistent are full-sky CTpp_eval and CTpp_cplm
  ! Intermittent are nmode size CTpp, CNTpp
  REAL(DP),    DIMENSION(:,:),   ALLOCATABLE   :: CTpp, CNTpp, CNpp
  REAL(DP),    DIMENSION(:),     ALLOCATABLE   :: CTpp_eval
  COMPLEX(DP), DIMENSION(:,:,:), ALLOCATABLE   :: CTpp_cplm

  ! Auxilliary arrays
  REAL(DP),    DIMENSION(:,:), ALLOCATABLE   :: VM

  REAL(DP),    DIMENSION(:,:), ALLOCATABLE   :: Wl, w8ring, w8pix
  REAL(DP),    DIMENSION(1:3,1:3)            :: euler_mat

! ====================================================
! Files
  CHARACTER(len=255)  :: map_signal_file, map_noise_file, map_mask_file
  CHARACTER(len=255)  :: w8_file, map_beam_file, ctpp_beam_file
  CHARACTER(len=255)  :: infile, fidfile

! ====================================================
! Global numerical variables
  REAL(DP) :: epsil, Ok, OmegaL, H0, map_beam_fwhm, ctpp_beam_fwhm, logdetCTpp
  INTEGER  :: nside, npix_fits, npol, ntot
  INTEGER  :: npix_cut, npix_cut_I, npix_cut_QU
  INTEGER  :: lmax, iseed, nsh
  INTEGER  :: n_evalues, nmode_cut

! ====================================================
! Job Control variables
  LOGICAL  :: do_mask, do_rotate, find_best_angles
  LOGICAL  :: add_noise, add_noise_diag, add_noise_cov
  LOGICAL  :: do_smooth_data, do_smooth_ctpp, do_smooth_noise

END MODULE TOPOLOGY_TYPES
