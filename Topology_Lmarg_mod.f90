MODULE Topology_Lmarg_mod
  !
  !Program to calculate th likelihood of a topology model
  !wrt the CMB data with various cuts.
  !
  USE nrtype
  USE Topology_types
  USE ctpp_eigen_mod
  USE Topology_map_mod_nel
  USE nr_minimization
  IMPLICIT NONE

  real(DP) :: ampl                  ! amplitude, available for local routines

  PRIVATE ampl

  CONTAINS

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! High level wrappers - they use only fixed global CTpp_evec/eval/lm 
! Can be safely called from outside module without prelim global set up
! All of them set CTpp and CNTpp=(CTpp+N)^-1 as byproduct
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

     SUBROUTINE  FIND_BEST_ANGLES_AND_AMPLITUDE(ampl_best,alpha,beta,gamma,LnL_max)
     real(DP), intent(out) :: ampl_best,alpha,beta,gamma,LnL_max

     integer(I4B)      :: iter,i
     real(DP)          :: tolerance=1.d-5
     real(DP)          :: y(4),p(1:4,1:3)      ! For amoeba,contains angles

        ampl=-1.0d0
        write(0,'(''Set simplex of angles - '',$)')
        CALL SET_amoeba_simplex(p,y)
        write(0,*)'Done'

        write(0,*)'Going into Amoeba',size(p,1),size(p,2),size(y)
        CALL amoeba(p,y,tolerance,LnLrotated_at_best_amplitude,iter)
        write(0,*)'Out of Amoeba',iter

        write(0,*)'#####################################################'
        do i=1,4
           write(0,'(a9,3(f7.5,1x),a9,f10.4)')' Angles: ',p(i,:),' LnL is: ',y(i)
        enddo
        write(0,*)'#####################################################'

        alpha=p(1,1)
        beta=p(1,2)
        gamma=p(1,3)

        ! The best model is in p(1,*),y(1) but the last calculation wasn't
        ! the best, so Ctpp and ampl is not correct. 
        ! Repeat the run with the best angles

        LnL_max=LnLrotated_at_best_amplitude(p(1,:))
        ampl_best=ampl
        return  
     END SUBROUTINE  FIND_BEST_ANGLES_AND_AMPLITUDE


     SUBROUTINE ROTATE_AND_FIND_BEST_AMPLITUDE(ampl_best,alpha,beta,gamma,LnL_max)
     real(DP), intent(in)    ::  alpha,beta,gamma
     real(DP), intent(out)   ::  ampl_best,LnL_max

     real(DP), dimension(3)  ::  p

     ampl=-1.0d0
     p(1)=alpha
     p(2)=beta
     p(3)=gamma

     LnL_max=LnLrotated_at_best_amplitude(p)
     ampl_best=ampl

     END SUBROUTINE ROTATE_AND_FIND_BEST_AMPLITUDE

     SUBROUTINE FIND_BEST_AMPLITUDE(ampl_best,LnL_max) 
     real(DP), intent(out)   :: ampl_best,LnL_max

     real(DP), allocatable,dimension(:,:) :: CTpp_evec_temp

! From global CTpp_evec and Ctpp_eval  reconstruct  Ctpp
! Reconstruction corrupts CTpp_evec, so work on copy
       allocate(CTpp_evec_temp(0:npix_fits-1,0:npix_fits-1))
       CTpp_evec_temp=CTpp_evec
       call RECONSTRUCT_FROM_EIGENVALUES(CTpp_evec_temp)
       deallocate(CTpp_evec_temp)

! Find best amplitude and store (Abest*C+N)^-1 in CNTpp. 
       ampl_best=-1.0d0                !Ininial guess for the amplitude
       LnL_max=LnL_bestampl(ampl_best) !Found ampl_best

       return
     END SUBROUTINE FIND_BEST_AMPLITUDE

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Low level drivers. Depend on the current state of intermittent global
! variables. Care must be taken to ensure they are correct before
! the routines are called.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

     FUNCTION LnLrotated_at_best_amplitude(ang)   ! driver for amoeba
     use nrtype
     use Topology_types
     use Topology_map_mod
     use LM_ROTATE

! 'ampl' is private global, invisible from outside, so calling this
! routine directly means no handle on initial guess for the amplitude.
! This is just to be able to pass amplitude from amoeba and as a bonus
! to have next amoeba step to use previous step amplitude as initial
! guess for speed up.

     real(DP), intent(in), dimension(:) :: ang  ! Vector of (alpha,beta,gamma)
     real(DP)                           ::  LnLrotated_at_best_amplitude

     real(DP), allocatable,dimension(:,:) :: CTpp_evec_rotated

! Rotate CTpp_evec (in CTpp_cplm form) into temporary CTpp_evec_rotated
       allocate(CTpp_evec_rotated(0:npix_fits-1,0:npix_fits-1))
       call rotate_ctpp(CTpp_evec_rotated,CTpp_cplm,nside,lmax,ang(1),ang(2),ang(3),.TRUE.)
!       CTpp_evec_rotated=Ctpp_evec
       write(0,'(a,$)')'I have rotated '
  
! From rotated CTpp_evec_rotated and global Ctpp_eval reconstruct cut-sky Ctpp
! CTpp_evec_rotated is corrupted on output
       call RECONSTRUCT_FROM_EIGENVALUES(CTpp_evec_rotated)
       deallocate(CTpp_evec_rotated)
       write(0,*)'and reconstructed'

! Find best amplitude and store (Abest*C+N)^-1 in CNTpp. 
! We use private global ampl to set initial guess and store the best-fit result 
       LnLrotated_at_best_amplitude=LnL_bestampl(ampl)
       write(0,'(a5,f8.4,a9,3(f7.5,1x),a9,f10.4)')'Ampl ',ampl,' Angles: ',ang,' LnL is: ',LnLrotated_at_best_amplitude
       return
     END FUNCTION LnLrotated_at_best_amplitude
  
     FUNCTION LnL_bestampl(ampl_best) 
     use Topology_map_mod
     ! Finds best amplitude and Log likelihood for globally given CTpp.
     ! as well as (CTpp+N)^-1 return as CNTpp

     real(DP), intent(inout) :: ampl_best
     real(DP)                :: LnL_bestampl
                              ! CTpp is global input
                              ! CNTpp is global output

     real(DP), parameter  :: err=1.d-3
     real(DP)             :: ax,bx,cx,fa,fb,fc,relerr

       ax=ampl_best       !Input value of ampl_best is used as initial guess
       bx=ax-0.1d0
       CALL mnbrak(ax,bx,cx,fa,fb,fc,LnLikelihood)
       write(0,*)'Bracketing:',bx,ax,cx,fb,fa,fc

       relerr=abs(err/bx)
       LnL_bestampl=brent(fb,ax,bx,cx,LnLikelihood,relerr,ampl_best)

       return
     END FUNCTION LnL_bestampl


     FUNCTION LnLikelihood(ampl_in)
     !Find -Ln(Likelihood) for global CTpp and set (CTpp+N)^-1

     REAL(DP), intent(in) :: ampl_in
                           ! CTpp is global input
                           ! CNTpp is global output

     INTEGER :: INFO, i, j 
     REAL(DP), DIMENSION(npix_cut) :: vec
     REAL(DP) :: LnLikelihood
     REAL(DP) :: trace, LnL_exp, ampl_0noisebest, LnL_0noisebest, DDOT
     REAL(DP), allocatable, dimension(:) :: eigen,WORK,D
     REAL(DP), allocatable, dimension(:,:) :: U,VT

!   scale CTpp and add noise.
!   Caution - only 'L' triangualr part in wmap_npp and thus CNTpp is valid
       CNTpp=CTpp*exp(ampl_in)

       IF (add_noise.and.do_smooth) THEN
          CNTpp=CNTpp+wmap_npp
       ENDIF
!   add diagonal noise. Perhaps set to 0 or set as regularizer or the noise weights of the pixle
       DO i = 0, npix_cut-1
          CNTpp(i,i) = CNTpp(i,i) + diag_noise(i)
       ENDDO

       ALLOCATE(D(0:npix_cut-1))
!       IF(SVD) THEN
!          ALLOCATE(WORK(0:5*npix_cut))
!          ALLOCATE(VT(0:npix_cut-1,0:npix_cut-1))
!          ALLOCATE(U(0:npix_cut-1,0:npix_cut-1))
!          INFO = 0
!          DO i = 0, npix_cut-1
!             DO j = i, npix_cut-1
!                CNTpp(i,j) = CNTpp(j,i)
!             ENDDO
!          ENDDO
!          INFO = 0
!!    Do general SVD 
!          CALL DGESVD('A','A',npix_cut,npix_cut,CNTpp,npix_cut,D,U,npix_cut,VT,npix_cut,WORK,5*npix_cut,INFO)
!          IF(INFO/=0) THEN
!             write(0,*) "DGESVD info=", INFO
!             STOP 'Error SVD DGESVD'
!          ENDIF
!
!          trace = 0.0d0
!          DO i = 0, mode_number
!             IF(abs(D(i)) == 0) THEN
!                WRITE(0,*) "Condition number not sufficient"
!                VT(i,:) = 0.0d0
!             ELSE
!                VT(i,:) = VT(i,:)/D(i)
!                trace=trace+0.5d0*LOG(D(i))
!             ENDIF
!          ENDDO
!          IF(mode_number<npix_cut-1) THEN
!             DO i=mode_number+1, npix_cut-1
!                VT(i,:) = 0.0d0
!             ENDDO
!          ENDIF
!!    Complete the invertion
!          CALL DGEMM('T','T',npix_cut,npix_cut,npix_cut,1.0d0,VT,npix_cut,&
!                    & U,npix_cut,0.0d0,CNTpp,npix_cut)
!          DEALLOCATE(WORK)
!          DEALLOCATE(VT)
!          DEALLOCATE(U)
!       ELSE
!   Cholesky decomposition of CNTpp
          INFO = 0
          CALL DPOTRF( 'L', npix_cut, CNTpp, npix_cut, INFO )
          IF (INFO/=0) THEN
             WRITE(0,'(a,i0)') 'NBad matrix ',INFO
             allocate(eigen(npix_cut),WORK(3*npix_cut))
             call DSYEV('V', 'L',npix_cut,CNTpp,npix_cut,eigen,WORK,3*npix_cut,INFO)
             write(0,*) eigen
             LnLikelihood = Top_bad_value
             deallocate(eigen,WORK)
             STOP
          ENDIF
!   Determinant calculation is now in Cholesky form, trace=log[sqrt(det(Ctpp))]
          trace = 0.0d0
          DO i = 0, npix_cut-1
             trace = trace + LOG(CNTpp(i,i))
          ENDDO

!   Finish inverting CNTpp
          INFO = 0
          CALL DPOTRI( 'L', npix_cut, CNTpp, npix_cut, INFO )
          IF(INFO/=0) STOP 'Error on DPOTRI'
!       ENDIF

!   Compute exponential part of -LnL
       call DSYMV('L',npix_cut,1.0d0,CNTpp,npix_cut,wmap_signal,1,0.0d0,vec,1)
       LnL_exp = DDOT(npix_cut,wmap_signal,1,vec,1)

!   trace has already been divided by 2
       LnLikelihood = 0.5d0*LnL_exp + trace
       write(0,*) 'Full noise:',ampl_in,LnLikelihood,trace,LnL_exp
!       write(0,*) 'D(mode_number+1):',D(mode_number+1)
       DEALLOCATE(D)

!   Analytic best amplitude and likelihood for zero noise - for information only
!       ampl_0noisebest= LnL_exp/npix_cut
!       LnL_0noisebest = 0.5d0*npix_cut*LOG(ampl_0noisebest)+trace
!       write(0,*) 'Zero noise:',ampl_0noisebest,LnL_0noisebest

       RETURN
     END FUNCTION LnLikelihood

     FUNCTION LmaxFisher() 
                       !Get 'Fisher matrix' diag element for amplitude at maxL
                       ! F(A) = 1/2 \sum_ij (Abest*C+N)^{-1}_ij C_ij
                       !Implicit input: (C+N)^-1 (at ampl=best) in CNtpp,
                       !                 C       (at ampl=1) in Ctpp
                       !                 npix_cut
     INTEGER :: i,j
     REAL(DP)    :: LmaxFisher
     REAL(DP), allocatable,dimension(:,:) :: CTN_1CT

       allocate(CTN_1CT(0:npix_cut-1,0:npix_cut-1))

       call DSYMM('L','L',npix_cut,npix_cut,1.0d0,CNTpp,npix_cut,CTpp,npix_cut,0.0d0,CTN_1CT,npix_cut)
    
       LmaxFisher=0.0d0
       do j=0,npix_cut-1
          do i=j+1,npix_cut-1
             LmaxFisher=LmaxFisher+2.0d0*CTN_1CT(i,j)**2
          enddo
          LmaxFisher=LmaxFisher+CTN_1CT(j,j)**2
       enddo

       LmaxFisher=0.5d0*LmaxFisher

       deallocate(CTN_1CT)
       RETURN
     END FUNCTION LmaxFisher

     FUNCTION LmaxCurv() 
                             !Get data projected part of curvature matrix
                             !Implicit input: (C+N)^-1 (ampl=best) in CNtpp,
                             !                 C       (ampl=1) in Ctpp
                             !                 wmap_signal
                             !                 npix_cut
     REAL(DP)    :: LmaxCurv,DDOT
     REAL(DP), allocatable,dimension(:) :: vec1,vec2

       allocate(vec1(0:npix_cut-1),vec2(0:npix_cut-1))

       call DSYMV('L',npix_cut,1.0d0,CNTpp,npix_cut,wmap_signal,1,0.0d0,vec1,1)
       call DSYMV('L',npix_cut,1.0d0,CTpp,npix_cut,vec1,1,0.0d0,vec2,1)
       call DSYMV('L',npix_cut,1.0d0,CNTpp,npix_cut,vec2,1,0.0d0,vec1,1)
       LmaxCurv=DDOT(npix_cut,vec1,1,vec2,1)
    
       deallocate(vec1,vec2)
       RETURN
     END FUNCTION LmaxCurv

     SUBROUTINE SET_amoeba_simplex(p,y)
     real(DP), intent(out), dimension(:,:) :: p
     real(DP), intent(out), dimension(:)   :: y

     integer :: i

        p(:,1) = (/ 0.0_dp, 0.0_dp, 0.0_dp, 0.1_dp /)    ! alphas
        p(:,2) = (/ 0.0_dp, 0.0_dp, 0.1_dp, 0.0_dp /)    ! betas
        p(:,3) = (/ 0.0_dp, 0.1_dp, 0.0_dp, 0.0_dp /)    ! gammas

        do i=1,4
           y(i)=LnLrotated_at_best_amplitude(p(i,:))
        enddo

        return
     END SUBROUTINE SET_amoeba_simplex

END MODULE Topology_Lmarg_mod
