MODULE TOPOLOGY_TYPES
  USE HEALPIX_TYPES

  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER :: code = 'TOPLIKE'
  CHARACTER(LEN=*), PARAMETER :: version = 'Summer2012'
  
! Control
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
  TYPE(EVALUE_CUT) :: evalue_cut_csky=EVALUE_CUT(S_NCUT,0,1200,0.0,.false.,0.0)

  REAL(DP),    PARAMETER :: Top_bad_value = 100000.d0

! CTpp normalization
  integer(I4B)  :: lnorm=10
  real(DP)      :: curlCl_in_mK = 1.d-2  ! 10x real to have ampl~-2

! Data
  CHARACTER(len=6), PARAMETER                :: expdata_format='PLANCK'
  REAL(DP),    DIMENSION(:),   ALLOCATABLE   :: map_signal
  REAL(DP),    DIMENSION(:,:), ALLOCATABLE   :: map_npp
  LOGICAL,     DIMENSION(:),   ALLOCATABLE   :: map_mask

! Global arrays to hold different matrices

  ! Intermittent data that operates in full sky (npix_fits, npix_fits) space
  REAL(DP),    DIMENSION(:,:), ALLOCATABLE, TARGET  :: FullSkyWorkSpace
  REAL(DP),    DIMENSION(:,:), POINTER  :: CTpp_full => NULL()
  REAL(DP),    DIMENSION(:,:), POINTER  :: CTpp_evec => NULL()
  REAL(DP),    DIMENSION(:,:), POINTER  :: CTpp_fid  => NULL()

  ! Persistent are full-sky CTpp_eval and CTpp_cplm
  ! Intermittent are nmode size CTpp, CNTpp
  REAL(DP),    DIMENSION(:,:), ALLOCATABLE   :: CTpp, CNTpp
  REAL(DP),    DIMENSION(:),   ALLOCATABLE   :: CTpp_eval
  COMPLEX(DP), DIMENSION(:,:), ALLOCATABLE   :: CTpp_cplm

  ! Auxilliary arrays
  REAL(DP),    DIMENSION(:,:), ALLOCATABLE   :: VM

  REAL(DP),    DIMENSION(:,:), ALLOCATABLE   :: Wl, w8ring, w8pix
  REAL(DP),    DIMENSION(1:3,1:3)            :: euler_mat

! Files
  CHARACTER(len=255)  :: map_signal_file, map_noise_file, map_mask_file
  CHARACTER(len=255)  :: w8_file, beam_file
  CHARACTER(len=255)  :: infile, fidfile
  CHARACTER(len=255)  :: fake_file

! Global numerical variables
  REAL(DP) :: epsil, Ok, OmegaL, H0, beam_fwhm, logdetCTpp
  INTEGER  :: nside, npix_cut, npix_fits, lmax, iseed, nsh
  INTEGER  :: n_evalues, nmode_cut

! Control variables
  LOGICAL  :: do_mask, do_rotate, find_best_angles, add_noise
  LOGICAL  :: do_Gsmooth, do_expsmooth
  LOGICAL  :: make_map, make_map_only, add_map_noise

END MODULE TOPOLOGY_TYPES
