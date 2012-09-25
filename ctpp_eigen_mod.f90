MODULE ctpp_eigen_mod
   ! Ctpp manipulations routines
   USE Topology_types
   USE nml_mod
   USE lm_rotate , only : rotate_ctpp
   IMPLICIT NONE

CONTAINS

   subroutine DECOMPOSE_AND_SAVE_EIGENVALUES()
   USE healpix_extras
   ! Full sky CTpp stored in CTpp_evec(npix_fits,npix_fits) -> 
   !        eigenfunctions in CTpp_evec + eigenvalues in CTpp_eval

   integer(I4B)                          :: INFO,ipix
   real(DP), allocatable, dimension(:)   :: WORK
   real(DP), allocatable, dimension(:,:) :: Bweights

      allocate(WORK(0:3*npix_fits-1))
      if ( w8ring(1,1) == 1.d0 ) then  ! Weights are trivial
         call DSYEV('V','L',npix_fits,CTpp_evec,npix_fits,CTpp_eval,WORK,3*npix_fits,INFO)
      else    ! General case, define eigenvectors orthonormal wrt weights

         ! Set diagonal matrix of weights Bweights
         allocate( Bweights(0:npix_fits-1, 0:npix_fits-1) )
         Bweights = 0.d0
         forall(ipix=0:npix_fits-1) Bweights(ipix,ipix)=w8pix(ipix,1)

         ! Solve generalized problem CTpp*B*evec = eval*evec
         call DSYGV(2,'V','L',npix_fits,CTpp_evec,npix_fits,Bweights,npix_fits,CTpp_eval,WORK,3*npix_fits,INFO)
         deallocate(Bweights)
      endif   
      if ( INFO /= 0 ) stop 'CTpp eigenvalue decomposition failed, terminating'

      return
   end subroutine DECOMPOSE_AND_SAVE_EIGENVALUES

   subroutine SORT_AND_LIMIT_EIGENVALUES()
      ! Reverse order of eigenvalues, and do not forget eigenvectors
      ! Each column is an eigenvector, columns correspond to different evalues
   integer(I4B)                          :: INFO,ipix,iev,inverse_ev
   real(DP), allocatable, dimension(:)   :: WORK
   real(DP)                              :: evalue_min, eval_temp

      allocate(WORK(0:npix_fits-1))
      do iev = 0, npix_fits/2-1
         inverse_ev = npix_fits-iev-1
         WORK(0:npix_fits-1)     = CTpp_evec(:,iev)
         CTpp_evec(:,iev)        = CTpp_evec(:,inverse_ev)
         CTpp_evec(:,inverse_ev) = WORK(0:npix_fits-1)
         eval_temp               = CTpp_eval(iev)
         CTpp_eval(iev)          = CTpp_eval(inverse_ev)
         CTpp_eval(inverse_ev)   = eval_temp
      enddo
      deallocate(WORK)

! Test output =======
!      do ipix=0,npix_fits-1
!         write(0,*) CTpp_eval(ipix)
!      enddo
! ====================

      where(CTpp_eval < 0.0_dp) CTpp_eval = 0.0_dp
      if ( evalue_cut%STRATEGY == evalue_cut%CONDITIONING ) then 
          evalue_min   = evalue_cut%condition_number*maxval(CTpp_eval)
          n_evalues    = count(CTpp_eval >= evalue_min)
! Test case
      else if ( evalue_cut%STRATEGY == evalue_cut%LCUT ) then
          n_evalues = (evalue_cut%lmax + 1)**2 - 4
          evalue_min = CTpp_eval(n_evalues-1)
      endif
! ====================
      write(0,*)evalue_min, n_evalues
      return
   end subroutine SORT_AND_LIMIT_EIGENVALUES


   subroutine RECONSTRUCT_FROM_EIGENVALUES()
   ! Eigenvalues in CTpp_eval + eigenfunctions in CTpp_evec -> 
   !        cut sky CTpp(npix_cut,npix_cut)
   ! Note: only significant eigenvalues (up to n_evalues) are used
   ! Note: CTpp_evec is corrupted on return

   integer :: ipix,jpix,ne,ic,jc

   ! Each masked eigenvector is packed into first npix_cut elements of a column
   ! The rest of a column is garbage.

      do ne=0,n_evalues-1
         CTpp_evec(0:npix_cut-1,ne)=pack(CTpp_evec(:,ne),map_mask)
         CTpp_evec(0:npix_cut-1,ne)=sqrt(CTpp_eval(ne))*CTpp_evec(0:npix_cut-1,ne)
      enddo

   ! Now n_evalues colums contain valid set of 
   ! normalized eigenvectors that are restricted to npix_cut length

      call DGEMM('N','T',npix_cut,npix_cut,n_evalues,1.0d0,CTpp_evec,npix_fits,CTpp_evec,npix_fits,0.0d0,CTpp,npix_cut)

      return
   end subroutine RECONSTRUCT_FROM_EIGENVALUES

   subroutine NORMALIZE_EIGENVALUES(eval)  ! Normalizes CTpp \sum(eigenvalues)=4 Pi
   real(DP) :: norm
   real(DP), DIMENSION(0:) :: eval
!      norm = 16.0d0*atan(1.0d0)/SUM(eval)
      norm = 1.d0/eval(0)
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
!          call rotate_ctpp(CTpp_evec_0,CTpp_cplm,nside,n_evalues,lmax,0.0d0,0.0d0,0.0d0,.TRUE.)
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
