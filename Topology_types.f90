MODULE TOPOLOGY_TYPES
  USE nrtype

  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER :: code = 'TOPLIKE'
  CHARACTER(LEN=*), PARAMETER :: version = 'Summer2012'
  
  REAL(DP),    PARAMETER :: Top_bad_value = 100000.d0
  REAL(DP),    PARAMETER :: Top_Evalue_precision = 1.d-4

! Data
  REAL(DP),    DIMENSION(:),   ALLOCATABLE   :: map_signal
  REAL(DP),    DIMENSION(:,:), ALLOCATABLE   :: map_npp
  LOGICAL,     DIMENSION(:),   ALLOCATABLE   :: map_mask

! Global arrays to hold different matrices
! Persistent are CTpp_eval and CTpp_cplm
! Variable are CTpp, CNTpp and CTpp_evec
  REAL(DP),    DIMENSION(:,:), ALLOCATABLE   :: CTpp, CNTpp
  REAL(DP),    DIMENSION(:,:), ALLOCATABLE   :: CTpp_evec
  REAL(DP),    DIMENSION(:),   ALLOCATABLE   :: CTpp_eval
  COMPLEX(DP), DIMENSION(:,:), ALLOCATABLE   :: CTpp_cplm

  REAL(DP),    DIMENSION(:,:), ALLOCATABLE   :: w8ring, w8pix
  REAL(DP),    DIMENSION(1:3,1:3)            :: euler_mat

  CHARACTER(len=255)  :: map_signal_file, map_noise_file, map_mask_file
  CHARACTER(len=255)  :: w8_file, beam_file
  CHARACTER(len=255)  :: infile
  CHARACTER(len=255)  :: fake_file

  REAL(DP) :: epsil, Ok, OmegaL, H0
  !REAL(DP) :: epsil, condition_num
  INTEGER  :: nside, npix_cut, npix_fits, lmax, nmaps, iseed, nsh
  INTEGER  :: n_evalues

  LOGICAL  :: do_mask, do_rotate, find_best_angles, add_noise, do_smooth
  LOGICAL  :: make_map, make_map_only, add_map_noise
  !LOGICAL  :: nice_output, SVD, First_time

END MODULE TOPOLOGY_TYPES
