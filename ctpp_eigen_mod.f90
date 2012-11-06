MODULE ctpp_eigen_mod
   ! Ctpp manipulations routines
   USE Topology_types
   IMPLICIT NONE

CONTAINS

   subroutine DECOMPOSE_AND_SAVE_EIGENVALUES()
   USE healpix_extras
   ! Full sky CTpp stored in CTpp_evec(npix_fits,npix_fits) -> 
   !        eigenfunctions in CTpp_evec + eigenvalues in CTpp_eval

   integer(I4B)                          :: INFO,ipix
   real(DP), allocatable, dimension(:)   :: WORK
   real(DP), allocatable, dimension(:,:) :: Bweights

      if ( .not.associated(CTpp_full,FullSkyWorkSpace) ) then
         stop 'Full sky CTpp matrix has not been set in FullSkyWorkSpace'
      endif
         
      allocate(WORK(0:3*npix_fits-1))
      if ( w8ring(1,1) == 1.d0 ) then  ! Weights are trivial
         call DSYEV('V','L',npix_fits,CTpp_full,npix_fits,CTpp_eval,WORK,3*npix_fits,INFO)
      else    ! General case, define eigenvectors orthonormal wrt weights

         ! Set diagonal matrix of weights Bweights
         allocate( Bweights(0:npix_fits-1, 0:npix_fits-1) )
         Bweights = 0.d0
         forall(ipix=0:npix_fits-1) Bweights(ipix,ipix)=w8pix(ipix,1)

         ! Solve generalized problem CTpp*B*evec = eval*evec
         call DSYGV(2,'V','L',npix_fits,CTpp_full,npix_fits,Bweights,npix_fits,CTpp_eval,WORK,3*npix_fits,INFO)
         deallocate(Bweights)
      endif   
      if ( INFO /= 0 ) stop 'CTpp eigenvalue decomposition failed, terminating'

! Full sky work space now contains eigenvectors
      CTpp_evec => FullSkyWorkSpace
      CTpp_full => NULL()

      return
   end subroutine DECOMPOSE_AND_SAVE_EIGENVALUES

   subroutine SORT_AND_LIMIT_EIGENVALUES()
      ! Reverse order of eigenvalues, and do not forget eigenvectors
      ! Each column is an eigenvector, columns correspond to different evalues
   integer(I4B)                          :: INFO,ipix,iev,inverse_ev
   real(DP), allocatable, dimension(:)   :: WORK
   real(DP)                              :: evalue_min, eval_temp

      if ( .not.associated(CTpp_evec,FullSkyWorkSpace) ) then
         stop 'CTpp_evec has not been set in FullSkyWorkSpace'
      endif

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
      call SET_N_EVALUES(CTpp_eval,evalue_cut_fsky,n_evalues,evalue_min)
! ====================
      write(0,*)evalue_min, n_evalues
      return
   end subroutine SORT_AND_LIMIT_EIGENVALUES


   subroutine RECONSTRUCT_FROM_EIGENVALUES()
   USE basis_modes
   ! Eigenvalues in CTpp_eval + eigenfunctions in CTpp_evec -> 
   !        cut sky CTpp(nmode_cut,nmode_cut) 
   ! Note: only significant eigenvalues (up to n_evalues) are used
   ! Note: CTpp_evec is corrupted on return

   integer :: ipix,jpix,ne,ic,jc

      if ( .not.associated(CTpp_evec,FullSkyWorkSpace) ) then
         stop 'CTpp_evec has not been set in FullSkyWorkSpace'
      endif
   ! Each masked eigenvector is packed into first npix_cut elements of a column
   ! The rest of a column is garbage.

      do ne=0,n_evalues-1
         CTpp_evec(0:npix_cut-1,ne)=pack(CTpp_evec(:,ne),map_mask)
         call project_vector_onto_basis_modes(CTpp_evec(0:npix_cut-1,ne))
         CTpp_evec(0:nmode_cut-1,ne)=sqrt(CTpp_eval(ne))*CTpp_evec(0:nmode_cut-1,ne)
      enddo

   ! Now n_evalues colums contain valid set of 
   ! normalized eigenvectors that are restricted to nmode_cut length

      call DGEMM('N','T',nmode_cut,nmode_cut,n_evalues,1.0d0,CTpp_evec,npix_fits,CTpp_evec,npix_fits,0.0d0,CTpp,nmode_cut)

   ! CTpp_evec data has been destroyed
      CTpp_evec => NULL()

      return
   end subroutine RECONSTRUCT_FROM_EIGENVALUES

   subroutine NORMALIZE_EIGENVALUES(eval)  ! Normalizes CTpp \sum(eigenvalues)=4 Pi
   real(DP), intent(inout), dimension(0:) :: eval

   integer(I4B)            :: l, lmcount
   real(DP)                :: flatnorm

   ! Normalize to power over l=2,lnorm being equal to that with flat curlCl
   ! Note: eigenvalues are eval ~ Cl * 4pi/npix
      flatnorm=0.0_dp
      do l=2,lnorm
         flatnorm = flatnorm + (l+0.5_dp)/((l+1.0_dp)*l)
      enddo
      flatnorm = flatnorm*curlCl_in_mK

      ! assuming eigenmodes are lm's from l=2
      lmcount=(lnorm+3)*(lnorm-1)
      write(0,*)'norm',lmcount,flatnorm,SUM(eval(0:lmcount-1)),eval(0)
      eval = eval/SUM(eval(0:lmcount-1))*flatnorm*npix_fits

      write(0,*)'Normalized over l=',lnorm,'to curlCl(mK)=',curlCl_in_mK

      return
   end subroutine NORMALIZE_EIGENVALUES

   subroutine INVERT_CTPP_SVD()
   ! Invert cut-sky CTpp using SVD decomposition

   real(DP), allocatable, dimension(:)   :: D, WORK
   real(DP), allocatable, dimension(:,:) :: U,VT
   integer(I4B)                          :: INFO, i, n_evalues_csky

      allocate(D(0:npix_cut-1))
      allocate(WORK(0:5*npix_cut))
      allocate(VT(0:npix_cut-1,0:npix_cut-1))
      allocate(U(0:npix_cut-1,0:npix_cut-1))

      call DGESVD('A','A',npix_cut,npix_cut,CTpp,npix_cut,D,U,npix_cut,VT,npix_cut,WORK,5*npix_cut,INFO)
      if(INFO/=0) then
         write(0,*) "DGESVD info=", INFO
         stop 'Error SVD DGESVD'
      endif

      call SET_N_EVALUES(D,evalue_cut_csky,n_evalues_csky)
      write(0,*) n_evalues_csky

      logdetCTpp = 0.0d0
      do i = 0, n_evalues_csky-1
         VT(i,:) = VT(i,:)/D(i)
         logdetCTpp=logdetCTpp+LOG(D(i))
!         write(0,*) i,D(i)
      enddo
      VT(n_evalues_csky:npix_cut-1,:) = 0.0d0

      call DGEMM('T','T',npix_cut,npix_cut,n_evalues_csky,1.0d0,VT,npix_cut,&
                   & U,npix_cut,0.0d0,CTpp,npix_cut)
      deallocate(WORK,VT,U)

      return
   end subroutine INVERT_CTPP_SVD

   subroutine SET_N_EVALUES(D,evalue_cut_in,n_evalues,evalue_min_out)
   REAL(DP),         intent(inout), dimension(0:) :: D
   TYPE(EVALUE_CUT), intent(in)                   :: evalue_cut_in
   INTEGER(I4B),     intent(out)                  :: n_evalues
   REAL(DP),         intent(out), optional        :: evalue_min_out

   REAL(DP)          :: evalue_min
   INTEGER(I4B)      :: npix

      npix=size(D,1)
      if ( D(0) < D(npix-1) ) then
         stop 'Eigenvalue array must be sorted in decreasing order'
      endif

      if ( evalue_cut_in%STRATEGY == S_NONE ) then
         n_evalues = npix
         evalue_min=0.0_dp
      else
         if ( evalue_cut_in%STRATEGY == S_CONDITIONING ) then 
            evalue_min   = evalue_cut_in%condition_number*D(0)
            n_evalues    = count(D >= evalue_min)
         else if ( evalue_cut_in%STRATEGY == S_LCUT ) then
            n_evalues = (evalue_cut_in%lmax + 1)**2 - 4
            evalue_min = D(n_evalues-1)
         else if ( evalue_cut_in%STRATEGY == S_NCUT ) then
            n_evalues = evalue_cut_in%nmax
            evalue_min = D(n_evalues-1)
         endif
         if ( evalue_cut_in%SET_BAD_MODES_HIGH ) then
            D(n_evalues:npix-1)=evalue_cut_in%BAD_MODES_NOISE
            n_evalues = npix
         endif
      endif

      if (n_evalues > npix) n_evalues=npix
      if (present(evalue_min_out)) evalue_min_out=evalue_min
        
      return
   end subroutine SET_N_EVALUES

END MODULE ctpp_eigen_mod
