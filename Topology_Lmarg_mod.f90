MODULE Topology_Lmarg_mod
  !
  !Program to calculate th likelihood of a topology model
  !wrt the CMB data with various cuts.
  !
  USE Topology_types
  USE ctpp_eigen_mod
  USE nr_minimization
  IMPLICIT NONE

  real(DP) :: ampl                  ! amplitude, available for local routines

  PRIVATE ampl

  CONTAINS

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! High level wrappers - they use only static global CTpp_eval/lm 
! Can be safely called from outside of the module without prelim set up
! All of them set CTpp and CNTpp=(CTpp+N)^-1 as byproduct
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

     SUBROUTINE  FIND_BEST_ANGLES_AND_AMPLITUDE(ampl_best,ang,LnL_max,ifrandom)
     real(DP), intent(inout)        :: ampl_best,ang(3)
     real(DP), intent(out)          :: LnL_max
     logical,  intent(in), optional :: ifrandom
! Works on global CTpp_cplm and CTpp_eval, 
! produces at best angles cut sky CTpp at unit amplitude
! and (Abest*CTpp + N)^-1 in CNTpp at best amplitude

     integer(I4B)      :: iter,i
     real(DP)          :: tolerance=1.d-5
     real(DP)          :: y(4),p(1:4,1:3)      ! For amoeba,contains angles
     logical           :: ifsuccess

        ampl=-1.0d0
        write(0,'(''Set simplex of angles - '',$)')
        CALL SET_amoeba_simplex(ang,p,y,ifrandom)
        write(0,*)'Done'

        write(0,*)'Going into Amoeba'
        CALL amoeba(p,y,tolerance,LnLrotated_at_best_amplitude,iter)
        write(0,*)'Out of Amoeba',iter

        write(0,*)'#####################################################'
        do i=1,4
          write(0,'(a9,3(f8.5,1x),a9,f12.4)')' Angles: ',p(i,:),' LnL is: ',y(i)
        enddo
        write(0,*)'#####################################################'

        ! The best model is in p(1,*),y(1) but the last calculation wasn't
        ! the best, so Ctpp and ampl is not correct. 
        ! Repeat the run with the best angles

        LnL_max=LnLrotated_at_best_amplitude(p(1,:))
        ampl_best=ampl
        call projectedS3_to_angles(p(1,1:3),ang,ifsuccess)

        return  
     END SUBROUTINE  FIND_BEST_ANGLES_AND_AMPLITUDE


     SUBROUTINE ROTATE_AND_FIND_BEST_AMPLITUDE(ampl_best,ang,LnL_max,ifrandom)
     real(DP), intent(inout)        :: ang(3)
     real(DP), intent(out)          :: ampl_best,LnL_max
     logical,  intent(in), optional :: ifrandom

! Works on global CTpp_cplm and CTpp_eval, 
! produces after rotation cutsky CTpp at unit amplitude
! and (Abest*CTpp + N)^-1 in CNTpp at best amplitude
     real(DP)  :: st_ang2=3.0_dp-2.0_dp*sqrt(2.0_dp)
     real(DP)  :: u(3),u2,q(0:3),spsi,stheta
     logical   :: ifsuccess

       ampl=-1.0d0
       if ( present(ifrandom) .and. ifrandom ) then
          do
            call random_number(u)

            !u = 2.0_dp*(u - 0.5_dp)
            !u2 = dot_product(u,u)
            !if (u2 > st_ang2) cycle

            q(0:1)=2.0_dp*u(1:2)-1.0_dp
            q(2)=cos(u(3)*PI)           
            q(3)=sin(u(3)*PI)           
            spsi=sqrt(1.0_dp-q(0)**2)
            stheta=sqrt(1.0_dp-q(1)**2)
            q(2:3)=q(2:3)*stheta
            q(1:3)=q(1:3)*spsi
            call quartenion_to_projectedS3(q,u,ifsuccess)

            if ( ifsuccess ) then 
               call projectedS3_to_angles(u,ang,ifsuccess)
               exit
            endif
          enddo
       else
          call angles_to_projectedS3(ang,u,ifsuccess)
       endif

       LnL_max=LnLrotated_at_best_amplitude(u)
       ampl_best=ampl

     END SUBROUTINE ROTATE_AND_FIND_BEST_AMPLITUDE

     SUBROUTINE FIND_BEST_AMPLITUDE(ampl_best,LnL_max) 
     real(DP), intent(out)   :: ampl_best,LnL_max

     real(DP), dimension(3)  :: ang=(/0.0_dp,0.0_dp,0.0_dp/)

       call ROTATE_AND_FIND_BEST_AMPLITUDE(ampl_best,ang,LnL_max)

     END SUBROUTINE FIND_BEST_AMPLITUDE

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Low level drivers. Depend on the current state of intermittent global
! variables. Care must be taken to ensure they are correct before
! the routines are called.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

     FUNCTION LnLrotated_at_best_amplitude(u)   ! driver for amoeba
     use Topology_types
     use Topology_map_mod
     use LM_ROTATE

! 'ampl' is private global, invisible from outside, so calling this
! routine directly means no handle on initial guess for the amplitude.
! This is just to be able to pass amplitude from amoeba and as a bonus
! to have next amoeba step to use previous step amplitude as initial
! guess for speed up.

     real(DP), intent(in), dimension(:) :: u ! Vector on S3 projection
     real(DP)                           :: LnLrotated_at_best_amplitude

     real(DP) :: ang(3),  DDOT
     logical  :: ifsuccess

! Convert projected coordinates into Euler angles and ensure the valid range
       write(0,*)'Rotated to u and angles'
       write(0,'(a,3(1x,e14.7))') 'u',u
       call projectedS3_to_angles(u,ang,ifsuccess)
       if ( .not.ifsuccess ) then
          LnLrotated_at_best_amplitude = Top_bad_value
          write(0,*)'Failed S3 projection, LnL set to bad value'
          return
       endif
          
! Rotate full sky eigevectors in CTpp_cplm form into current CTpp_evec
       CTpp_evec => FullSkyWorkSpace
       call rotate_ctpp(CTpp_evec,CTpp_cplm,nside,n_evalues,npol,lmax,ang(1),ang(2),ang(3),.TRUE.)
       write(0,'(a,3(1x,f10.7))') 'a',ang

! From rotated CTpp_evec and global Ctpp_eval reconstruct cut-sky CTpp 
! then project CTpp onto mode basis
! CTpp_evec is deassociated on output
       call RECONSTRUCT_FROM_EIGENVALUES()
       write(0,*)'and reconstructed'

! Find best amplitude and store (Abest*C+N)^-1 in CNTpp. 
! We use private global ampl to set initial guess and store the best-fit result 
       LnLrotated_at_best_amplitude=LnL_bestampl(ampl)
       write(0,'(a5,f8.4,a9,3(f8.5,1x),a9,f12.4)')'Ampl ',ampl,' Angles: ',ang,' LnL is: ',LnLrotated_at_best_amplitude
       return
     END FUNCTION LnLrotated_at_best_amplitude
  
     FUNCTION LnL_bestampl(ampl_best) 
     use Topology_map_mod
     ! Finds best amplitude and Log likelihood for globally given CTpp.
     ! as well defines (Abest*CTpp+N)^-1 in CNTpp

     real(DP), intent(inout) :: ampl_best
     real(DP)                :: LnL_bestampl
                              ! CTpp is global input
                              ! CNTpp is global output

     real(DP), parameter  :: err=1.d-3
     real(DP)             :: ax,bx,cx,fa,fb,fc,relerr

       ax=ampl_best       !Input value of ampl_best is used as initial guess
       bx=ax-0.5d0
       write(0,*)'Start bracketing Min'
       CALL mnbrak(ax,bx,cx,fa,fb,fc,LnLikelihood)
       write(0,*)'Bracketing:',bx,ax,cx,fb,fa,fc

       relerr=abs(err/bx)
       LnL_bestampl=brent(fb,ax,bx,cx,LnLikelihood,relerr,ampl_best)
       write(0,*) 'Best amplitude', ampl_best
       return
     END FUNCTION LnL_bestampl


     FUNCTION LnLikelihood(ampl_in)
     !Find -Ln(Likelihood) for global CTpp and set (ampl CTpp + N)^-1
     USE ct_io
     REAL(DP), intent(in) :: ampl_in
                           ! CTpp, map_npp and map_signal is global input
                           ! All are projected on a chosen basis
                           ! CNTpp is global output in the same basis

     INTEGER :: INFO, i, j 
     REAL(DP), DIMENSION(nmode_cut) :: vec
     REAL(DP) :: LnLikelihood
     REAL(DP) :: trace, LnL_exp, ampl_0noisebest, LnL_0noisebest, DDOT
     REAL(DP), allocatable, dimension(:) :: eigen,WORK
     
!   scale CTpp and add noise.
       CNTpp=CTpp*exp(ampl_in)
       CNTpp=CNTpp+CNpp

       INFO = 0
       CALL DPOTRF( 'L', nmode_cut, CNTpp, nmode_cut, INFO )
       IF (INFO/=0) THEN
          WRITE(0,'(a,i0)') 'NBad matrix ',INFO
          WRITE(0,*) "ampl_in", ampl_in
          allocate(eigen(nmode_cut),WORK(34*nmode_cut))
       ! Resets the corrupted CNTpp matrix to get the real eigen values
          CNTpp=CTpp*exp(ampl_in)
          CNTpp=CNTpp+CNpp
          call WriteCTpp('cntpp',CNTpp,nmode_cut,1)
          call WriteCTpp('cnpp',CNpp,nmode_cut,1)
          call WriteCTpp('ctpp',CTpp,nmode_cut,1)
          call DSYEV('V', 'L',nmode_cut,CNTpp,nmode_cut,eigen,WORK,34*nmode_cut,INFO)
          write(0,*) eigen
          LnLikelihood = Top_bad_value
          deallocate(eigen,WORK)
          STOP
       ENDIF
!   Determinant calculation is now in Cholesky form, trace=log[sqrt(det(Ctpp))]
       trace = 0.0d0
       DO i = 0, nmode_cut-1
          trace = trace + LOG(CNTpp(i,i))
       ENDDO

!   Finish inverting CNTpp
       INFO = 0
       CALL DPOTRI( 'L', nmode_cut, CNTpp, nmode_cut, INFO )
       IF(INFO/=0) STOP 'Error on DPOTRI'

!   Compute exponential part of -LnL
       call DSYMV('L',nmode_cut,1.0d0,CNTpp,nmode_cut,map_signal,1,0.0d0,vec,1)
       LnL_exp = DDOT(nmode_cut,map_signal,1,vec,1)

!   trace has already been divided by 2
       LnLikelihood = 0.5d0*LnL_exp + trace
       write(0,*) 'Full noise:',ampl_in,LnLikelihood,trace,LnL_exp

!   Analytic best amplitude and likelihood for zero noise - for information only
!       ampl_0noisebest= LnL_exp/nmode_cut
!       LnL_0noisebest = 0.5d0*nmode_cut*LOG(ampl_0noisebest)+trace
!       write(0,*) 'Zero noise:',ampl_0noisebest,LnL_0noisebest

       RETURN
     END FUNCTION LnLikelihood

     FUNCTION LmaxFisher() 
                       !Get 'Fisher matrix' diag element for amplitude at maxL
                       ! F(A) = 1/2 \sum_ij (Abest*C+N)^{-1}_ij C_ij
                       !Implicit input: (C+N)^-1 (at ampl=best) in CNtpp,
                       !                 C       (at ampl=1) in Ctpp
                       !                 nmode_cut
     INTEGER :: i,j
     REAL(DP)    :: LmaxFisher, trace
     REAL(DP), allocatable,dimension(:,:) :: CTN_1CT

       allocate(CTN_1CT(0:nmode_cut-1,0:nmode_cut-1))

       call DSYMM('L','L',nmode_cut,nmode_cut,1.0d0,CNTpp,nmode_cut,CTpp,nmode_cut,0.0d0,CTN_1CT,nmode_cut)
    
       LmaxFisher=0.0d0
       trace=0.d0
       do j=0,nmode_cut-1
          do i=j+1,nmode_cut-1
             LmaxFisher=LmaxFisher+2.0d0*CTN_1CT(i,j)**2
          enddo
          LmaxFisher=LmaxFisher+CTN_1CT(j,j)**2
          trace = trace + CTN_1CT(j,j)
       enddo

       write(0,*)'Gradient balance: Trace', trace

       LmaxFisher=0.5d0*LmaxFisher

       deallocate(CTN_1CT)
       RETURN
     END FUNCTION LmaxFisher

     FUNCTION LmaxCurv() 
                             !Get data projected part of curvature matrix
                             !Implicit input: (C+N)^-1 (ampl=best) in CNtpp,
                             !                 C       (ampl=1) in Ctpp
                             !                 map_signal
                             !                 nmode_cut
     REAL(DP)    :: LmaxCurv,DDOT
     REAL(DP), allocatable,dimension(:) :: vec1,vec2

       allocate(vec1(0:nmode_cut-1),vec2(0:nmode_cut-1))

       call DSYMV('L',nmode_cut,1.0d0,CNTpp,nmode_cut,map_signal,1,0.0d0,vec1,1)
       write(0,*)'exp logL', DDOT(nmode_cut,vec1,1,map_signal,1)
       call DSYMV('L',nmode_cut,1.0d0,CTpp,nmode_cut,vec1,1,0.0d0,vec2,1)
       write(0,*)'Gradient balance: exponential part:', DDOT(nmode_cut,vec1,1,vec2,1)
       call DSYMV('L',nmode_cut,1.0d0,CNTpp,nmode_cut,vec2,1,0.0d0,vec1,1)
       LmaxCurv=DDOT(nmode_cut,vec1,1,vec2,1)
    
       deallocate(vec1,vec2)
       RETURN
     END FUNCTION LmaxCurv

     SUBROUTINE SET_amoeba_simplex(ang,p,y,ifrandom)
     real(DP), intent(in)                  :: ang(3)
     real(DP), intent(out), dimension(:,:) :: p
     real(DP), intent(out), dimension(:)   :: y
     logical,  intent(in),  optional       :: ifrandom

     real(DP) :: st_ang=sqrt(2.0_dp)-1.0_dp, sq3=sqrt(3.0_dp), st_ang2
     real(DP) :: u(3), u2
     integer  :: i 
     logical  :: ifsuccess

         ! Amoeba start in terms of Euler angles
!        p(:,1) = (/ 0.0_dp, 0.0_dp, 0.0_dp, HALFPI /)    ! alphas
!        p(:,2) = (/ 0.0_dp, 0.0_dp, st_ang, st_ang /)    ! betas
!        p(:,3) = (/ 0.0_dp, st_ang, 0.0_dp, -HALFPI/)    ! gammas

        ! S3 coordinates of amoeba are rotated by p/6 in beta and gamma
!        p(:,1) = (/ 0.0_dp, 0.25_dp*sq3, -0.5_dp,     0.75_dp /)
!        p(:,2) = (/ 0.0_dp, 0.25_dp,      0.5_dp*sq3, 0.25_dp*sq3 /) 
!        p(:,3) = (/ 0.0_dp, 0.5_dp*sq3,   0.0_dp,    -0.5_dp /) 

        ! S3 coordinates of amoeba aligned with coordinate axes
!        p(:,1) = (/ 0.0_dp, 0.0_dp,  0.0_dp, 1.0_dp /)
!        p(:,2) = (/ 0.0_dp, 0.0_dp,  1.0_dp, 0.0_dp /) 
!        p(:,3) = (/ 0.0_dp, 1.0_dp,  0.0_dp, 0.0_dp /) 
!        p = p*st_ang

        ! Random projected S3 simplex
        st_ang2=st_ang**2
        if ( present(ifrandom) .and. ifrandom ) then
           i=0
        else
           ! except the first vertex
           call angles_to_projectedS3(ang,u,ifsuccess)
           if (.not.ifsuccess) then
              write(0,*)'Unable to initialize amoeba with supplied angles, abort'
              stop
           endif
           i=1
           p(i,:)=u
        endif
        do
          call random_number(u)
          u = 2.0_dp*(u - 0.5_dp)
          u2 = dot_product(u,u)
          if (u2 > st_ang2) cycle
          i=i+1
          p(i,:)=u
          if (i == 4) exit
        enddo

        do i=1,4
           y(i)=LnLrotated_at_best_amplitude(p(i,:))
        enddo

        write(0,*) 'Starting simplex:'
        do i=1,4
           write(0,'(5x,3(d12.4,1x),d15.7)') p(i,:),y(i)
        enddo

        return
     END SUBROUTINE SET_amoeba_simplex

     SUBROUTINE projectedS3_to_angles(u,ang,ifsuccess)
     real(DP), intent(in),  dimension(3)  :: u   ! Projected S3 cartezian
     real(DP), intent(out), dimension(3)  :: ang ! Euler angles
     logical,  intent(out)                :: ifsuccess

     real(DP)  :: q(0:3)

     ! First step - from S3 projected onto equatorial plane to quartenion
     call projectedS3_to_quartenion(u,q,ifsuccess)
     if ( .not.ifsuccess ) return

     ! Second step - from quartenion to Euler angles,
     call quartenion_to_angles(q,ang)

     return
     END SUBROUTINE projectedS3_to_angles

     SUBROUTINE angles_to_projectedS3(ang,u,ifsuccess)
     real(DP), intent(in),  dimension(3)  :: ang
     real(DP), intent(out), dimension(3)  :: u
     logical,  intent(out)                :: ifsuccess

     real(DP)  :: q(0:3)

     ! First step - from Euler angles to quartenion
     call angles_to_quartenion(ang,q)
     
     ! Second step - from quartenion to projected u on equatorial plane
     call quartenion_to_projectedS3(q,u,ifsuccess)

     return
     END SUBROUTINE angles_to_projectedS3

     SUBROUTINE quartenion_to_projectedS3(q,u,ifsuccess)
     real(DP), intent(in),  dimension(0:3) :: q
     real(DP), intent(out), dimension(3)   :: u
     logical,  intent(out)                 :: ifsuccess

     if ( q(0) < 0.0_dp ) then
        write(0,*)'quartenion does not correspond to northen semisphere of S3'
        ifsuccess = .false.
     endif
 
     u(1:3) = q(1:3)/(1.0_dp+q(0))

     ifsuccess = .true.
     return
     END SUBROUTINE quartenion_to_projectedS3

     SUBROUTINE projectedS3_to_quartenion(u,q,ifsuccess)
     real(DP), intent(in),  dimension(3)   :: u   ! Projected S3 cartezian
     real(DP), intent(out), dimension(0:3) :: q   ! quartenion
     logical,  intent(out)                 :: ifsuccess

     real(DP)  :: u2, DDOT

     u2 = DDOT(3,u,1,u,1)
     if (u2 > 1.0_dp) then
        ifsuccess = .false.
        return
     endif
     q(0)   = (1.0_dp-u2)/(1.0_dp+u2)
     q(1:3) = (2.0_dp/(1.0_dp+u2))*u(1:3)

     ifsuccess = .true.
     return

     END SUBROUTINE projectedS3_to_quartenion

     SUBROUTINE quartenion_to_angles(q,ang)
     real(DP), intent(in),  dimension(0:3) :: q   ! quartenion
     real(DP), intent(out), dimension(3)   :: ang ! Euler angles

     real(DP)  :: q03, q12, ppf, pmf

     ! from quartenion to Euler angles,
     ! ang=(phi,theta,psi) in Healpix ZYZ (fix axis) convention

     q03 = sqrt(q(0)**2 + q(3)**2)
     q12 = sqrt(q(1)**2 + q(2)**2)
     ang(2) = 2.0_dp*acos(q03)

     if (q03 > q12) then
        ppf = atan2(q(3),q(0)) 
        if ( q12 < 1.0d-6 ) then
           pmf = 0.0_dp
        else
           pmf = atan2(q(1),q(2))
        endif
     else
        pmf = atan2(q(1),q(2))
        if ( q03 < 1.0d-6) then
           ppf = 0.0_dp
        else
           ppf = atan2(q(3),q(0))
        endif
     endif

     ang(1) = ppf - pmf
     ang(3) = ppf + pmf

     return
     END SUBROUTINE quartenion_to_angles

     SUBROUTINE angles_to_quartenion(ang,q)
     real(DP), intent(in),  dimension(3)   :: ang
     real(DP), intent(out), dimension(0:3) :: q

     real(DP)  :: ppf,pmf,cost2,sint2,cosppf,sinppf,cospmf,sinpmf

     ppf = 0.5_dp*(ang(1)+ang(3))
     pmf = 0.5_dp*(ang(3)-ang(1))
     cost2  = cos(ang(2)/2.0_dp)
     sint2  = sin(ang(2)/2.0_dp)
     cosppf = cos(ppf)
     sinppf = sin(ppf)
     cospmf = cos(pmf)
     sinpmf = sin(pmf)
     
     q(0) = cost2*cosppf
     q(1) = sint2*sinpmf
     q(2) = sint2*cospmf
     q(3) = cost2*sinppf

     ! Choose quartenion in Nothern hemisphere of S3
     if ( q(0) < 0.0_dp ) q = -q
     return

     END SUBROUTINE angles_to_quartenion

END MODULE Topology_Lmarg_mod
