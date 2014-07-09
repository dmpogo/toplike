MODULE ctpp_eigen_mod
   ! Ctpp manipulations routines
   USE Topology_types
   IMPLICIT NONE

CONTAINS

   subroutine DECOMPOSE_AND_SAVE_EIGENVALUES()
   USE healpix_extras
   ! Full sky CTpp stored in CTpp_evec(ntot,ntot) -> 
   !        eigenfunctions in CTpp_evec + eigenvalues in CTpp_eval

   integer(I4B)                          :: INFO,ip
   real(DP), allocatable, dimension(:)   :: WORK
   real(DP), allocatable, dimension(:,:) :: Bweights

      if ( .not.associated(CTpp_full,FullSkyWorkSpace) ) then
         stop 'Full sky CTpp matrix has not been set in FullSkyWorkSpace'
      endif
         
      allocate(WORK(0:34*ntot-1))
      if ( w8ring(1,1) == 1.d0 ) then 
         ! Trivial weights
         call DSYEV('V','L',ntot,CTpp_full,ntot,CTpp_eval,WORK,34*ntot,INFO)
      else   
         ! General case,  C B x = lambda x eigenproblem
         write(0,*)'Using ring weights'
         allocate( Bweights(0:ntot-1, 0:ntot-1) )
         Bweights=0.0_dp
         forall(ip=0:ntot-1)
            Bweights(ip,ip)=w8pix(mod(ip,npix_fits),ip/npix_fits+1)
         end forall
         call DSYGV(2,'V','L',ntot,CTpp_full,ntot,Bweights,ntot,CTpp_eval,WORK,34*ntot,INFO)
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

      allocate(WORK(0:ntot-1))
      do iev = 0, ntot/2-1
         inverse_ev = ntot-iev-1
         WORK(0:ntot-1)          = CTpp_evec(:,iev)
         CTpp_evec(:,iev)        = CTpp_evec(:,inverse_ev)
         CTpp_evec(:,inverse_ev) = WORK(0:ntot-1)
         eval_temp               = CTpp_eval(iev)
         CTpp_eval(iev)          = CTpp_eval(inverse_ev)
         CTpp_eval(inverse_ev)   = eval_temp
      enddo
      deallocate(WORK)

! Test output =======
!      do ipix=0,ntot-1
!         write(0,*) CTpp_eval(ipix)
!      enddo
! ====================

      where(CTpp_eval < 0.0_dp) CTpp_eval = 0.0_dp
      call SET_N_EVALUES(CTpp_eval,evalue_cut_fsky,n_evalues,evalue_min)

      write(0,*)evalue_min, n_evalues
      return
   end subroutine SORT_AND_LIMIT_EIGENVALUES


   subroutine RECONSTRUCT_FROM_EIGENVALUES()
   ! Eigenvalues in CTpp_eval + eigenfunctions in CTpp_evec -> 
   !        cut sky CTpp(nmode_cut,nmode_cut) 
   ! Note: only significant eigenvalues (up to n_evalues) are used
   ! Note: CTpp_evec is corrupted on return

   integer :: ne

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

      call DGEMM('N','T',nmode_cut,nmode_cut,n_evalues,1.0d0,CTpp_evec,ntot,CTpp_evec,ntot,0.0d0,CTpp,nmode_cut)

   ! CTpp_evec data has been destroyed
      CTpp_evec => NULL()

      return
   end subroutine RECONSTRUCT_FROM_EIGENVALUES

   subroutine NORMALIZE_EIGENVALUES(ctppnorm_out)
   real(DP), intent(out), optional  :: ctppnorm_out

   integer(I4B)            :: l, lmcount
   real(DP)                :: flatnorm, ctppnorm

   ! Normalize to power over l=2,lnorm being equal to that with flat curlCl
   ! Note: eigenvalues are eval ~ Cl * 4pi/npix
      flatnorm=0.0_dp
      do l=2,lnorm
         flatnorm = flatnorm + (l+0.5_dp)/((l+1.0_dp)*l)*Wl(l,1)**2
      enddo
      flatnorm = flatnorm*curlCl_in_mK

      ! assuming eigenmodes are lm's from l=2
      lmcount=(lnorm+3)*(lnorm-1)
      ctppnorm=flatnorm*npix_fits/SUM(CTpp_eval(0:lmcount-1))
      write(0,*)'norm ',flatnorm, SUM(CTpp_eval(0:lmcount-1)), ctppnorm
      CTpp_eval = CTpp_eval*ctppnorm

      write(0,*)'Normalized over l=',lnorm,'to curlCl(mK)=',curlCl_in_mK

      if (present(ctppnorm_out)) ctppnorm_out=ctppnorm
      return
   end subroutine NORMALIZE_EIGENVALUES

   subroutine INVERT_CTPP_SVD()
   ! Invert cut-sky CTpp using SVD decomposition

   real(DP), allocatable, dimension(:)   :: D, WORK
   real(DP), allocatable, dimension(:,:) :: U,VT
   integer(I4B)                          :: INFO, i, n_evalues_csky

      allocate(D(0:npix_cut-1))
      allocate(WORK(0:34*npix_cut))
      allocate(VT(0:npix_cut-1,0:npix_cut-1))
      allocate(U(0:npix_cut-1,0:npix_cut-1))

      call DGESVD('A','A',npix_cut,npix_cut,CTpp,npix_cut,D,U,npix_cut,VT,npix_cut,WORK,34*npix_cut,INFO)
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
   ! Returns number of n_evalues, not index in the array
   REAL(DP),         intent(inout), dimension(0:) :: D
   TYPE(EVALUE_CUT), intent(in)                   :: evalue_cut_in
   INTEGER(I4B),     intent(out)                  :: n_evalues
   REAL(DP),         intent(out), optional        :: evalue_min_out

   REAL(DP)          :: evalue_min
   INTEGER(I4B)      :: npix
   LOGICAL           :: descending=.false.

      npix=size(D,1)
      if ( D(0) < D(npix-1) ) descending=.true.

      if ( evalue_cut_in%STRATEGY == S_NONE ) then
         n_evalues = count(D > 0.0_dp)
         evalue_min= minval(D, mask = D > 0.0_dp)
      else
         if ( evalue_cut_in%STRATEGY == S_CONDITIONING ) then 
            evalue_min   = evalue_cut_in%condition_number*maxval(D)
            n_evalues    = count(D >= evalue_min)
         else if ( evalue_cut_in%STRATEGY == S_LCUT ) then
            n_evalues = (evalue_cut_in%lmax + 1)**2 - 4
            if ( descending ) then
               evalue_min = D(n_evalues-1)
            else
               evalue_min = D(npix-n_evalues)
            endif
         else if ( evalue_cut_in%STRATEGY == S_NCUT ) then
            n_evalues = evalue_cut_in%nmax
            if ( descending ) then
               evalue_min = D(n_evalues-1)
            else
               evalue_min = D(npix-n_evalues)
            endif
         endif
         if ( evalue_cut_in%SET_BAD_MODES_HIGH ) then
            where(D < evalue_min) D=evalue_cut_in%BAD_MODES_NOISE
            n_evalues = npix
         endif
      endif

      if (n_evalues > npix) n_evalues=npix
      if (present(evalue_min_out)) evalue_min_out=evalue_min
        
      return
   end subroutine SET_N_EVALUES


   subroutine SET_BASIS_MODES()
   ! Basis modes are either
   !     1) eigenmodes of fiducial CTpp_fid on a cut sky
   ! or  2) signal to noise eigenmode for N^{-1} CTpp_fid on a cut sky

   real(DP), dimension(0:npix_cut-1)     :: eval, Bweights_diag
   real(DP), dimension(:,:), allocatable :: Bweights
   real(DP), dimension(34*npix_cut)      :: WORK
   real(DP)                              :: scaleI,scaleP
   integer(I4B)                          :: INFO,ip,ic,iev

      if ( .not.associated(CTpp_fid,FullSkyWorkSpace) ) then
         stop 'CTpp_fid has not been set in FullSkyWorkSpace'
      endif

      ! Possibly scale polarization part to have same power as temperature
      if ( BASIS_TYPE < 0 ) then
         write(0,*)'Scaling polarization to temperature power for mode basis'
         scaleI = 0.0_dp
         scaleP = 0.0_dp
         do ip=0,npix_fits-1
            scaleI=scaleI+CTpp_fid(ip,ip)
            scaleP=scaleP+CTpp_fid(npix_fits+ip,npix_fits+ip)+CTpp_fid(2*npix_fits+ip,2*npix_fits+ip)
         enddo
         scaleP=sqrt(scaleI/scaleP)
         write(0,*)'Scale factor=',scaleP
         do ip=npix_fits,ntot-1
            CTpp_fid(:,ip) = CTpp_fid(:,ip)*scaleP
            CTpp_fid(ip,:) = CTpp_fid(ip,:)*scaleP
         enddo
      endif

      ! Compactify CTpp_fid and Bweights to cut sky
      ic=0
      do ip=0,ntot-1
         if (map_mask(ip)) then
            CTpp_fid(0:npix_cut-1,ic)=pack(CTpp_fid(:,ip),map_mask)
            Bweights_diag(ic)=w8pix(mod(ip,npix_fits),ip/npix_fits+1)
            ic=ic+1
         endif
      enddo

      if ( BASIS_TYPE == 0 ) then
         ! Solve generalized problem CTpp*evec = eval*Npp*evec
         ! Will fail if map_npp is not strictily positive definite
         write(0,*)'signal to noise basis'
         allocate ( Bweights(0:npix_cut-1, 0:npix_cut-1) )
         Bweights = map_npp
         
         ! Instead using DSYGV driver, make explicit steps, omitting
         ! the last redefinition of the eigenvector x=L^T{-1} y
         ! returned y is normalized as y^T y = 1, while x^T N x = 1
         call DPOTRF(    'L',npix_cut,Bweights,npix_cut,INFO )
         if ( INFO /= 0 ) stop 'map_npp is not strictly positive definite'
         call DSYGST(  1,'L',npix_cut,CTpp_fid,ntot,Bweights,npix_cut,INFO )
         call DSYEV ('V','L',npix_cut,CTpp_fid,ntot,eval,WORK,34*npix_cut,INFO )

!         call DSYGV(1,'V','L',npix_cut,CTpp_fid,ntot,Bweights,npix_cut,eval,WORK,3*npix_cut,INFO)
         ! Normalize eigenvectors to x B x = I
         do iev=0,npix_cut-1
            do ic = 0,npix_cut-1
               CTpp_fid(ic,iev) = CTpp_fid(ic,iev)/sqrt(Bweights_diag(ic))
            enddo
         enddo
      else if ( ABS(BASIS_TYPE) == 1 ) then
         ! Solve generalized problem CTpp*B*evec = eval*evec
         write(0,*)'weighted eigenvalue basis'
         allocate ( Bweights(0:npix_cut-1, 0:npix_cut-1) )
         forall (ic=0:npix_cut-1) Bweights(ic,ic)=Bweights_diag(ic)
         call DSYGV(2,'V','L',npix_cut,CTpp_fid,ntot,Bweights,npix_cut,eval,WORK,34*npix_cut,INFO)
         deallocate(Bweights)
      else
         ! Direct, unweighted eigenvector decomposition
         write(0,*)'unweighted eigenvalue basis'
         call DSYEV ('V','L',npix_cut,CTpp_fid,ntot,eval,WORK,3*npix_cut,INFO )
      endif

      ! Set nmode_cut count of eigenvalues left after the cut
      call SET_N_EVALUES(eval,evalue_cut_csky,nmode_cut)
      write(0,*) 'selected number of modes=',nmode_cut

      ! Reorder the sorting to decreasing and store eigenvectors
      if (allocated(VM)) deallocate(VM)
      allocate( VM(0:npix_cut-1,0:nmode_cut-1) )
      forall(iev=0:nmode_cut-1) VM(:,iev) = CTpp_fid(:,npix_cut-1-iev)
      
      CTpp_fid => NULL()
      return
   end subroutine SET_BASIS_MODES

   subroutine PROJECT_MATRIX_ONTO_BASIS_MODES(mat)
   ! Projects matrix onto a set of globally saved basis modes
   real(DP), intent(inout), dimension(:,:) :: mat
   
   real(DP), allocatable, dimension(:,:)   :: mat_temp

      allocate(mat_temp(npix_cut,nmode_cut))
      call DSYMM('L','L',npix_cut,nmode_cut,1.0_dp,mat,npix_cut,VM,npix_cut,0.0_dp,mat_temp,npix_cut)
      call DGEMM('T','N',nmode_cut,nmode_cut,npix_cut,1.0_dp,VM,npix_cut,mat_temp,npix_cut,0.0_dp,mat,npix_cut)
      deallocate(mat_temp)

      return
   end subroutine PROJECT_MATRIX_ONTO_BASIS_MODES

   subroutine PROJECT_VECTOR_ONTO_BASIS_MODES(vec)
   ! Projects vector onto a set of globally saved basis modes
   real(DP), intent(inout), dimension(0:) :: vec
   
   real(DP), allocatable, dimension(:)   :: vec_temp

      allocate(vec_temp(nmode_cut))
      call DGEMV('T',npix_cut,nmode_cut,1.0_dp,VM,npix_cut,vec,1,0.0_dp,vec_temp,1)
      vec(0:nmode_cut-1)=vec_temp
      deallocate(vec_temp)

      return
   end subroutine PROJECT_VECTOR_ONTO_BASIS_MODES

END MODULE ctpp_eigen_mod
