MODULE TOPOLOGY_TYPES
  USE nrtype

  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER :: code = 'TOPLIKE'
  CHARACTER(LEN=*), PARAMETER :: version = 'Summer2010'
  
  REAL(DP),    PARAMETER :: Top_bad_value = 100000.d0

  REAL,        DIMENSION(:,:), ALLOCATABLE   :: map, map_noise, map_mask
  REAL(DP),    DIMENSION(:),   ALLOCATABLE   :: wmap_signal, wmap_noise, diag_noise
  REAL(DP),    DIMENSION(:,:), ALLOCATABLE   :: wmap_npp
  REAL(DP),    DIMENSION(:,:), ALLOCATABLE   :: wmap_qhat
  LOGICAL,     DIMENSION(:),   ALLOCATABLE   :: wmap_mask

  REAL(DP),    DIMENSION(:,:), ALLOCATABLE   :: CTpp, CNTpp
  REAL(DP),    DIMENSION(:,:), ALLOCATABLE   :: CTpp_evec
  REAL(DP),    DIMENSION(:),   ALLOCATABLE   :: CTpp_eval
  COMPLEX(DP), DIMENSION(:,:), ALLOCATABLE   :: CTpp_cplm
!  COMPLEX(DP), DIMENSION(:,:), ALLOCATABLE   :: CTpp_clm

  REAL(DP),    DIMENSION(:,:), ALLOCATABLE   :: w8ring
  REAL(DP),    DIMENSION(1:3,1:3)            :: euler_mat

  CHARACTER(len=80) :: wmap_signal_file, wmap_noise_file, wmap_mask_file
  CHARACTER(len=80) :: w8_file, beam_file
  CHARACTER(len=100) :: infile
  CHARACTER(len=100) :: fake_file

  REAL(DP) :: epsil
  !REAL(DP) :: epsil, condition_num
  INTEGER  :: nside, npix_cut, npix_fits, lmax, nmaps, iseed
  !INTEGER  :: nside, npix_cut, npix_fits, lmax, nmaps, iseed, mode_number

  LOGICAL  :: do_mask, make_map, add_noise, do_rotate, add_map_noise
  LOGICAL  :: do_smooth, find_best_angles, nice_output, First_time
  LOGICAL  :: make_map_only
  !LOGICAL  :: do_smooth, find_best_angles, nice_output, SVD, First_time

  NAMELIST /toplike/  wmap_signal_file, wmap_noise_file,  wmap_mask_file, &
       nside, lmax, w8_file, beam_file, do_rotate, find_best_angles, &
       nice_output, make_map, make_map_only, add_map_noise, do_smooth, &
       iseed, epsil
       !iseed, SVD, epsil, mode_number

  
END MODULE TOPOLOGY_TYPES
