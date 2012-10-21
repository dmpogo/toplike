MODULE TOPOLOGY_TYPES
  USE HEALPIX_TYPES

  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER :: code = 'TOPLIKE'
  CHARACTER(LEN=*), PARAMETER :: version = 'Summer2012'
  
  TYPE TOP_EVALUE_CUT
       INTEGER(I4B)    :: NONE=0, LCUT=1, CONDITIONING=2
       INTEGER(I4B)    :: STRATEGY=1
       INTEGER(I4B)    :: LMAX=35
       REAL(DP)        :: CONDITION_NUMBER=1.d-4
  END TYPE TOP_EVALUE_CUT

  TYPE(TOP_EVALUE_CUT)   :: evalue_cut

  REAL(DP),    PARAMETER :: Top_bad_value = 100000.d0

! Data
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

  REAL(DP),    DIMENSION(:,:), ALLOCATABLE   :: Wl, w8ring, w8pix
  REAL(DP),    DIMENSION(1:3,1:3)            :: euler_mat

  CHARACTER(len=255)  :: map_signal_file, map_noise_file, map_mask_file
  CHARACTER(len=255)  :: w8_file, beam_file
  CHARACTER(len=255)  :: infile
  CHARACTER(len=255)  :: fake_file

  REAL(DP) :: epsil, Ok, OmegaL, H0, beam_fwhm
  !REAL(DP) :: epsil, condition_num
  INTEGER  :: nside, npix_cut, npix_fits, lmax, nmaps, iseed, nsh
  INTEGER  :: n_evalues

  LOGICAL  :: do_mask, do_rotate, find_best_angles, add_noise
  LOGICAL  :: do_Gsmooth, do_expsmooth
  LOGICAL  :: make_map, make_map_only, add_map_noise
  !LOGICAL  :: nice_output, SVD, First_time

END MODULE TOPOLOGY_TYPES
