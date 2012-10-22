MODULE basis_modes
   USE Topology_types
   IMPLICIT NONE

CONTAINS

   subroutine SET_BASIS_MODES()
   ! Basis modes are eigenmodes of fiducial CTpp_fid on a cut sky

   real(DP),allocatable,dimension(:,:) :: Bweights
   real(DP),   dimension(0:npix_cut-1) :: eval
   real(DP),   dimension(3*npix_cut)   :: WORK
   integer(I4B)                        :: INFO,ip,ic

      if ( .not.associated(CTpp_fid,FullSkyWorkSpace) ) then
         stop 'CTpp_fid has not been set in FullSkyWorkSpace'
      endif

      allocate( Bweights(0:npix_cut-1, 0:npix_cut-1) )
      Bweights = 0.0_dp

      ic=0
      do ip=0,npix_fits-1
         if (map_mask(ip)) then
            CTpp_fid(0:npix_cut-1,ic)=pack(CTpp_fid(:,ip),map_mask)
            Bweights(ic,ic)=w8pix(ip,1)
            ic=ic+1
         endif
      enddo
   
      ! Solve generalized problem CTpp*B*evec = eval*evec
      call DSYGV(2,'V','L',npix_cut,CTpp_fid,npix_fits,Bweights,npix_cut,eval,WORK,3*npix_cut,INFO)
      deallocate(Bweights)

      !short term fix, consider using more sophisticated cut
      nmode_cut=evalue_cut_csky%NMAX

      if (allocated(VM)) deallocate(VM)
      allocate( VM(0:npix_cut-1,0:nmode_cut-1) )
      forall(ic=0:nmode_cut-1) VM(:,ic) = CTpp_fid(:,npix_cut-1-ic)
      
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

END MODULE basis_modes
