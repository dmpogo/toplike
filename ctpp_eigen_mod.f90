MODULE ctpp_eigen_mod
   ! Ctpp manipulations routines
   USE nrtype
   USE Topology_types
   USE nml_mod
   USE lm_rotate , only : rotate_ctpp
   IMPLICIT NONE

CONTAINS

   subroutine DECOMPOSE_AND_SAVE_EIGENVALUES()
   ! Full sky CTpp stored in CTpp_evec(npix_fits,npix_fits) -> 
   !        eigenfunctions in CTpp_evec + eigenvalues in CTpp_eval

   integer :: ival,INFO
   real(DP), allocatable, dimension(:) :: WORK
      allocate(WORK(0:3*npix_fits-1))
      call DSYEV('V','L',npix_fits,CTpp_evec,npix_fits,CTpp_eval,WORK,3*npix_fits,INFO)
      deallocate(WORK)

      do ival=0,npix_fits-1
         if (CTpp_eval(ival) < 0.0_dp ) CTpp_eval(ival)=0.0_dp
      enddo
      return
   end subroutine DECOMPOSE_AND_SAVE_EIGENVALUES


   subroutine RECONSTRUCT_FROM_EIGENVALUES(CTpp_evec_work)
   ! Eigenvalues in CTpp_eval + eigenfunctions in CTpp_evec_work -> 
   !        cut sky CTpp(npix_cut,npix_cut)

   real(DP), intent(inout), dimension(0:npix_fits-1,0:npix_fits-1) :: CTpp_evec_work

   integer :: ipix,jpix,ne,ic,jc

   ! Each masked eigenvector is packed into first npix_cut elements of a column
   ! The rest of a column is garbage.

      do ne=0,npix_fits-1
         CTpp_evec_work(0:npix_cut-1,ne)=pack(CTpp_evec_work(:,ne),wmap_mask)
         CTpp_evec_work(0:npix_cut-1,ne)=sqrt(CTpp_eval(ne))*CTpp_evec_work(0:npix_cut-1,ne)
      enddo

   ! Now npix_fits rows and npix_cut columns contain valid set of 
   ! normalized eigenvectors

      call DGEMM('N','T',npix_cut,npix_cut,npix_fits,1.0d0,CTpp_evec_work,npix_fits,CTpp_evec_work,npix_fits,0.0d0,CTpp,npix_cut)

      return
   end subroutine RECONSTRUCT_FROM_EIGENVALUES

   subroutine NORMALIZE_EIGENVALUES(eval)  ! Normalizes CTpp \sum(eigenvalues)=4 Pi
   real(DP) :: norm
   real(DP), DIMENSION(:) :: eval
      norm = 16.0d0*atan(1.0d0)/SUM(eval)
      eval = eval*norm

      return
   end subroutine NORMALIZE_EIGENVALUES

!   subroutine MODECOUNTER()
!   integer :: i,j,INFO
!   real(DP), dimension(:), allocatable :: D,WORK
!   real(DP), dimension(:,:), allocatable :: VT,U
!   real(DP), dimension(:,:), allocatable :: Ctpp_evec_0
!       IF (do_rotate) THEN
!          allocate(CTpp_evec_0(0:npix_fits-1,0:npix_fits-1))
!          call rotate_ctpp(CTpp_evec_0,CTpp_cplm,nside,lmax,0.0d0,0.0d0,0.0d0,.TRUE.)
!          call RECONSTRUCT_FROM_EIGENVALUES(Ctpp_evec_0)
!          deallocate(CTpp_evec_0)
!       ENDIF
!
!       IF(SVD) THEN
!          ALLOCATE(D(0:npix_cut-1))
!          ALLOCATE(WORK(0:5*npix_cut))
!          ALLOCATE(VT(0:npix_cut-1,0:npix_cut-1))
!          ALLOCATE(U(0:npix_cut-1,0:npix_cut-1))
!          INFO = 0
!          DO i = 0, npix_cut-1
!             DO j = i, npix_cut-1
!                CTpp(i,j) = CTpp(j,i)
!             ENDDO
!          ENDDO
!          INFO = 0
!    Do general SVD 
!          CALL DGESVD('A','A',npix_cut,npix_cut,CTpp,npix_cut,D,U,npix_cut,VT,npix_cut,WORK,5*npix_cut,INFO)
!          IF(INFO/=0) THEN
!             write(0,*) "DGESVD info=", INFO
!             STOP 'Error SVD DGESVD'
!          ENDIF
!    Determine the number of usable modes
!          DO i = 0, npix_cut-1
!             if((D(i)/D(0))<condition_num) then
!                mode_number = i-1
!                goto 9001
!             endif
!          ENDDO
!          mode_number = npix_cut-1
!9001      continue
!          DEALLOCATE(D)
!          DEALLOCATE(WORK)
!          DEALLOCATE(VT)
!          DEALLOCATE(U)
!       ELSE
!NELSON ADD CHOLESKY mode number
!       ENDIF
!       write(0,*) "Number of modes =",mode_number
!       write(0,*) "Condition number =",condition_num
!       return
!   end subroutine MODECOUNTER
END MODULE ctpp_eigen_mod
