MODULE LM_ROTATE
  USE ALM_TOOLS
  USE PIX_TOOLS
  USE HEALPIX_TYPES

  REAL(DP), DIMENSION(:, :), ALLOCATABLE :: dlmm
  REAL(DP), DIMENSION(:, :), ALLOCATABLE :: plm
  INTEGER                                :: n_plm=0

  PRIVATE dlmm,plm,n_plm
contains  
    
  !-----------------------------------------------------------------------------
  ! ROTATE_CTPP - takes in the matrix of eigenvectors of the covariance matrix
  ! and rotates them in using D^l_mm' rotations. Eigenvectors are calculated 
  ! before calling this subroutine with e.g. dsyev
  !
  ! a,b,g  : euler angles for the rotation
  ! CUT_MD : decides if the monopole and dipole terms are set to  zero
  ! kernel : a function to smooth the alms with e.g. B(l) must be of length 
  !          0-->lmax. Set them all to 1.d0 if you don't want any smoothing
  !
  ! WARNING : dimesnion of ctpp must be for full sky. If you want to cut the 
  ! sky do it after calling this subroutine. This is so that we can use the 
  ! HEALPIX rshts.
  !
  ! Carlo      24Jun04     CITA
  ! Carlo      02Jul04     IAP   fixed small bug 
  !                              checked everything works nicely
  ! Carlo      04Jul04     IoA   added monopole and dipole cut 
  !                              added kernel smoothing
  ! Carlo      13Oct04     CITA  added precomputed plm support
  !                              must use topology_types now
  ! Carlo      26Oct04     CITA  Changed to single precision dlmm
  ! Dmitri     19Jun05     CITA  
  !				 1) Make dlmm local to the module
  !				    New flag redo_dlmm to rotate_ctpp() 
  !				    says whether to compute new dlmm's.
  !				    Should be .TRUE. when 'b' changes
  ! Dmitri     17Nov05     UofA  Conversion to HealPix 2.01
  ! Dmitri     25Jun07     IAP   
  !                              1) Speed up rotate_ctpp() by factor 5
  !                                 by cleaning inner loop
  !                              2) Corrected parallelization
  ! Dmitri     Jan03       UofA  include polarization
  !-----------------------------------------------------------------
  
  SUBROUTINE rotate_ctpp(ctpp, cplm, nside, n_evalues, npol, lmax, a, b, g, redo_dlmm)
    IMPLICIT NONE

    REAL(DP),    intent(out), DIMENSION(0:, 0:)     :: ctpp
    COMPLEX(DP), intent(in),  DIMENSION(1:, 0:, 0:) :: cplm
    INTEGER,     intent(in) :: nside,npol,lmax,n_evalues
    REAL(DP),    intent(in) :: a, b, g
    LOGICAL,     intent(in) :: redo_dlmm
    
    INTEGER                                         :: npix, l, m, mp, p, ip, indl
    REAL(DP),    DIMENSION(0:12*nside**2-1, 1:npol) :: map
    COMPLEX(DP), DIMENSION(1:npol, 0:lmax, 0:lmax)  :: alm
    COMPLEX(DP), DIMENSION(0: lmax)                 :: ea,eg
    REAL(DP)                                        :: sgn

    integer :: OMP_GET_THREAD_NUM

    IF (redo_dlmm.or.(.not.allocated(dlmm)).or.(size(dlmm,1)/=lmax*(lmax+2)+1)) THEN
       CALL get_dlmm(lmax,b)
       ! redo_dlmm=.TRUE.        To return diagnostics, intent should be inout
    ENDIF

!    write(0,*)'Writing out dlmm'
!    do m=0,2
!       do mp=-2,2
!          write(0,*) m,mp,dlmm(2+m,2+mp)
!       enddo
!    enddo

    ! Rotate the lm eigenvectors and simultaneously back transform in place 
    ! into ctpp()

    do m = 0,lmax
       ea(m) = exp(DCMPLX(0.d0, -DBLE(m)*a))
       eg(m) = exp(DCMPLX(0.d0, -DBLE(m)*g))
    enddo

    npix=nside2npix(nside)

!$OMP PARALLEL DO DEFAULT(SHARED),SCHEDULE(STATIC), PRIVATE(p,map,alm,indl,l,m,mp,sgn)

    DO p = 0,n_evalues-1

!$OMP CRITICAL
!test       write(0,*)'---After rotation-------', p
       alm = DCMPLX(0.d0, 0.d0)
       DO l = 0, lmax
          indl = l**2 + l
          DO m = 0, l
             alm(:, l, m) = ea(m) * dlmm(indl + m, indl) * cplm(:, indl, p)  ! mp=0
             sgn = -1.d0
             DO mp = 1, l
                alm(:, l, m) = alm(:, l, m) + ea(m) * &
                     (dlmm(indl+m,indl+mp)*cplm(:,indl+mp,p)*eg(mp) + &
                     sgn*dlmm(indl+m,indl-mp)*CONJG(cplm(:,indl+mp,p)*eg(mp)))
                sgn = -sgn
             ENDDO
          ENDDO
!test       if (l <= lmax  ) then
!test          write(0,*) l, alm(1,l,0:l)
!test       endif
       ENDDO
       IF(n_plm .eq. nside*(lmax+1)*(lmax+2) ) THEN
          CALL alm2map(nside, lmax, lmax, alm, map, plm(0:n_plm-1,1))
       ELSE
          CALL alm2map(nside, lmax, lmax, alm, map)
       ENDIF

       DO ip=1,npol
          ctpp((ip-1)*npix:ip*npix-1,p) = map(:,ip)
       ENDDO
!$OMP END CRITICAL

!test          write(0,*)'----------------------'
    ENDDO
!$OMP END PARALLEL DO
    !DEALLOCATE(dlmm)

!    stop
    RETURN
  END SUBROUTINE rotate_ctpp

  SUBROUTINE getcplm(cplm, ctpp, nside, n_evalues, npol, lmax, w8ring, kernel, cut_md)
    ! The lm-transform of the eigenvector matrix.
    IMPLICIT NONE
    
    COMPLEX(DP), intent(out), allocatable, DIMENSION(:,:,:) :: cplm
    REAL(DP), intent(in), DIMENSION(0:, 0:)       :: ctpp
    INTEGER,  intent(in)                          :: nside, npol, lmax, n_evalues
    REAL(DP), intent(in), DIMENSION(:,:)          :: w8ring
    REAL(DP), intent(in), optional, DIMENSION(:)  :: kernel
    LOGICAL,  intent(in), optional                :: cut_md

    COMPLEX(DP), DIMENSION(1:npol, 0:lmax, 0:lmax)   :: alm
    REAL(DP),    DIMENSION(0:12*nside**2 - 1,1:npol) :: map
    REAL(DP),    DIMENSION(2)                        :: zbounds=(/0.0_dp, 0.0_dp/)
    INTEGER(I4B)                                     :: npix, l, m, p, ip, indl, lmin
    INTEGER(I4B), PARAMETER                          :: iter_order=5
    
    ! Cut out the monopole and dipole if requested.
    IF (present(cut_md).and.cut_md) THEN
      lmin = 2
    ELSE
      lmin = 0
    END IF

    allocate( cplm(1:npol,0:lmax*(lmax+2),0:n_evalues-1) )
    cplm = DCMPLX(0.d0,0.d0)
   
    ! Precompute plm's   (should check if it is faster)
    ! n_plm=nside*(lmax+1)*(lmax+2)
    ! if (.not.allocated(plm)) then
    !   allocate(plm(0:n_plm-1,1))
    ! endif

    !call plm_gen(nside,lmax,lmax,plm)

    npix=nside2npix(nside)

!$OMP PARALLEL DO DEFAULT(SHARED),SCHEDULE(STATIC), PRIVATE(p,map,alm,indl,l,m)
    DO p = 0, n_evalues - 1

!$OMP CRITICAL       
       do ip=1,npol
          map(:,ip) = ctpp((ip-1)*npix:ip*npix-1,p)
       enddo
       IF (n_plm .ne. 0) THEN
          CALL map2alm_iterative(nside,lmax,lmax,iter_order,map,alm,zbounds,w8ring,plm(0:n_plm-1,1:1))
       ELSE
          CALL map2alm_iterative(nside,lmax,lmax,iter_order,map,alm,zbounds,w8ring)
       ENDIF

       if (present(kernel) ) then             ! smooth
          DO l = 0, lmax
             alm(:,l,0:l) = kernel(l) * alm(:, l, 0:l)
          ENDDO
       endif

       DO l = lmin, lmax
          indl = l**2 + l
          DO m = 0, l
             cplm(:, indl+m, p) =  alm(:, l, m)
          ENDDO
       ENDDO
!$OMP END CRITICAL       

    ENDDO
!$OMP END PARALLEL DO

    RETURN
  END SUBROUTINE getcplm

  SUBROUTINE smooth_ctpp(ctpp, nside, n_evalues, lmax, w8ring, kernel, cut_md)
    ! Smoothing of the eigenvectors. This routine needs a thought
    IMPLICIT NONE
    
    REAL(DP), DIMENSION(0:12*nside**2-1, 0:n_evalues-1), intent(inout) :: ctpp
    INTEGER,  intent(in) :: nside, lmax, n_evalues
    REAL(DP), intent(in), DIMENSION(1: 2 * nside, 1: 1)  :: w8ring
    REAL(DP), intent(in), DIMENSION(0: lmax), optional   :: kernel
    LOGICAL,  intent(in), optional :: cut_md

    COMPLEX(DP), DIMENSION(1: 1, 0: lmax, 0: lmax) :: alm
    REAL(DP),    DIMENSION(0: 12 * nside**2 - 1)   :: map
    REAL(DP),    DIMENSION(2)                   :: zbounds=(/0.0_dp, 0.0_dp/)
    INTEGER :: l, m, p
    
    ! The lm-transform of the eigenvector matrix.
    DO p = 0, n_evalues - 1
       map(:) = ctpp(:,p)
       IF(n_plm .ne. 0) THEN
          CALL map2alm(nside, lmax, lmax, map, alm, zbounds, w8ring, plm(0:n_plm-1,1))
       ELSE
          CALL map2alm(nside, lmax, lmax, map, alm, zbounds, w8ring)
       ENDIF

       if(present(cut_md).and.CUT_MD) THEN   ! cut monopole and dipole
          alm(1,0,0) = CMPLX(0.d0,0.d0)
          alm(1,1,0) = CMPLX(0.d0,0.d0)
          alm(1,1,1) = CMPLX(0.d0,0.d0)
       ENDIF

       if (present(kernel) ) then             ! smooth
          DO l = 0, lmax
             alm(1,l,0:l) = kernel(l) * alm(1, l, 0:l)
          ENDDO
       endif

       IF(n_plm .ne. 0) THEN
          CALL alm2map(nside, lmax, lmax, alm, map, plm(0:n_plm-1,1))
       ELSE
          CALL alm2map(nside, lmax, lmax, alm, map)
       ENDIF
       ctpp(:,p) = map(:)
    ENDDO

    RETURN
  END SUBROUTINE smooth_ctpp
    
  SUBROUTINE add2inverse(mat,sigma,n)
    !=================================================================
    ! Special case of Sherman-Morrison formula for a diagonal
    ! addition to the matrix e.g.
    !   A^-1 ---> (A + diag(u))^-1
    ! here sigma(i) is assumed to be the square root of u
    ! and mat is the pre computed inverse of A
    !=================================================================
    IMPLICIT NONE
    
    INTEGER :: n
    DOUBLE PRECISION, DIMENSION(1:n,1:n) :: mat,v
    DOUBLE PRECISION, DIMENSION(1:n) :: sigma
    
    DOUBLE PRECISION :: lambda
    INTEGER :: i,k

    lambda = 0.d0
    DO k=1,n
       lambda = lambda + mat(k,k)*sigma(k)**2
       v(k,:) = sigma(k)*mat(k,:)
    ENDDO
    
    DO k=1,n
       DO i=1,n
          mat(i,:) = mat(i,:) - v(k,i)*v(k,:)/(1.d0 + lambda)
       ENDDO
    ENDDO

    RETURN
  END SUBROUTINE add2inverse

  SUBROUTINE get_dlmm(lmax,b)
    !=================================================================
    !  Subroutine to calculate the d^l_mm'(beta) for the computation of
    !  the rotation matrix a_lm' = sum_m D^l_m'm (alpha,beta,gamma) a_lm
    !
    !  with D^l_mm' = exp(-i.m.alpha) d^l_mm'(beta) exp(-i.m'.gamma)
    !  i.e. you have to sandwich the output of this routine byt the
    !  exponential to get the actual rotation matrix
    !
    !  Uses recursive relations in m and l so best strategy is to run 
    !  all ls upto lmax
    !
    !  Refs: 
    !  Blanco, Florez & Bermejo (ECCC3 1997)
    !  Varshalovich, Moskalev & Khersonskii
    !
    !  Carlo      21Jun04     CITA
    !  Carlo      18Oct04     CITA       explicitly solved l*cos(B)-m
    !                                    conflict
    !  Dmitri     19Jun05     CITA       dlmm is global for the module
    !=================================================================
    IMPLICIT NONE
    INTEGER,          INTENT(in) :: lmax
    DOUBLE PRECISION, INTENT(in) :: b

    INTEGER :: l,m,mp, indl,indlm1,indlm2
    DOUBLE PRECISION :: cb,sb,cb2,sb2,tb2
    DOUBLE PRECISION :: dl,dm,dmp

    IF (allocated(dlmm)) THEN
       deallocate(dlmm)            ! start fresh
    ENDIF
    allocate(dlmm(0:lmax*(lmax+2),0:lmax*(lmax+2)))

    !some trig functions
    cb = dcos(dble(b))
    sb = dsin(dble(b))
    cb2 = dcos(dble(b)/2.d0)
    sb2 = dsin(dble(b)/2.d0)
    tb2 = dtan(dble(b)/2.d0)

    dlmm(:,:) = 0.d0

    !start with the smallest l
    dlmm(0,0) =  1.d0
    dlmm(2+0,2+0) = cb
    dlmm(2+1,2-1) = sb2**2
    dlmm(2+1,2+0) = -1.d0/dsqrt(2.d0) * sb
    dlmm(2+0,2-1) = dlmm(2+1,2+0)
    dlmm(2+1,2+1) = cb2**2
    dlmm(2-1,2-1) = dlmm(2+1,2+1)

    DO l = 2,lmax
       dl = DBLE(l)
       indl = l**2+l
       indlm1 = l**2-l
       indlm2 = l**2-3*l+2
       DO m= 0,l-2
          dm = DBLE(m)
          DO mp = -m, m
             dmp = DBLE(mp)
             dlmm(indl+m,indl+mp) = dl*(2.d0*dl-1.d0)/SQRT((dl**2-dm**2)*(dl**2-dmp**2))*&
                  ((dlmm(2+0,2+0)-dm*dmp/(dl*(dl-1.d0)))*dlmm(indlm1+m,indlm1+mp) - &
                  SQRT(((dl-1.d0)**2-dm**2)*((dl-1.d0)**2-dmp**2))/(dl-1.d0)/       &
                  (2.d0*dl-1.d0)*dlmm(indlm2+m,indlm2+mp))
             dlmm(indl-mp,indl-m) = dlmm(indl+m,indl+mp)
          ENDDO
       ENDDO
       dlmm(indl+l,indl+l) = dlmm(3,3)*dlmm(indlm1+l-1,indlm1+l-1)
       dlmm(indl-l,indl-l) = dlmm(indl+l,indl+l) 
       dlmm(indl+l-1,indl+l-1) = (dl*dlmm(2,2)-dl+1.d0)*dlmm(indlm1+l-1,indlm1+l-1)
       dlmm(indl-l+1,indl-l+1) = dlmm(indl+l-1,indl+l-1)

       DO mp = l,-l+1,-1
          dmp = DBLE(mp)
          dlmm(indl+l,indl+mp-1) = -dsqrt((dl+dmp)/(dl-dmp+1.d0))*tb2*dlmm(indl+l,indl+mp)
          dlmm(indl-mp+1,indl-l) = dlmm(indl+l,indl+mp-1)
       ENDDO
       DO mp=l-1,2-l,-1
          dmp = DBLE(mp)
          if(dlmm(indl+l-1,indl+mp).ne.0.d0) then
             dlmm(indl+l-1,indl+mp-1) = -(dl*cb-dmp+1.d0)/(dl*cb-dmp)*dsqrt((dl+dmp)/(dl-dmp+1.d0))*tb2*dlmm(indl+l-1,indl+mp)
             dlmm(indl-mp+1,indl-l+1) = dlmm(indl+l-1,indl+mp-1)
          else
             dlmm(indl+l-1,indl+mp-1) = 0.d0
             dlmm(indl-mp+1,indl-l+1) = dlmm(indl+l-1,indl+mp-1)
          endif
       ENDDO
       !Final symmetrization
       DO m = -l,l
          DO mp = m+1,l
             dlmm(indl+m,indl+mp) = (-1.d0)**(m-mp)*dlmm(indl+mp,indl+m)
          ENDDO
       ENDDO
    ENDDO

    !symmetrize l=1 block
    l=1
    indl = 2
    DO m = -l,l
       DO mp = m+1,l
          dlmm(indl+m,indl+mp) = (-1.d0)**(m-mp)*dlmm(indl+mp,indl+m)
       ENDDO
    ENDDO

    RETURN
  END SUBROUTINE get_dlmm

END MODULE LM_ROTATE
