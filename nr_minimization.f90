MODULE nr_minimization
   USE nrtype

   INTERFACE amoeba
       MODULE PROCEDURE amoeba_r, amoeba_d
   END INTERFACE
   INTERFACE brent
       MODULE PROCEDURE brent_r, brent_d
   END INTERFACE
   INTERFACE mnbrak
       MODULE PROCEDURE mnbrak_r, mnbrak_d
   END INTERFACE

CONTAINS
	SUBROUTINE amoeba_r(p,y,ftol,func,iter)
	USE nrtype; USE nrutil, ONLY : assert_eq,imaxloc,iminloc,nrerror,swap
	IMPLICIT NONE
	INTEGER(I4B), INTENT(OUT) :: iter
	REAL(SP), INTENT(IN) :: ftol
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: p
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: ITMAX=5000
	REAL(SP), PARAMETER :: TINY=1.0e-10
	INTEGER(I4B) :: ihi,ndim
	REAL(SP), DIMENSION(size(p,2)) :: psum
	call amoeba_private
	CONTAINS
!BL
	SUBROUTINE amoeba_private
	IMPLICIT NONE
	INTEGER(I4B) :: i,ilo,inhi
	REAL(SP) :: rtol,ysave,ytry,ytmp
	ndim=assert_eq(size(p,2),size(p,1)-1,size(y)-1,'amoeba')
	iter=0
	psum(:)=sum(p(:,:),dim=1)
	do
	  ilo=iminloc(y(:))
	  ihi=imaxloc(y(:))
	  ytmp=y(ihi)
	  y(ihi)=y(ilo)
	  inhi=imaxloc(y(:))
	  y(ihi)=ytmp
	  rtol=2.0_sp*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo))+TINY)
	  if (rtol < ftol) then
	     call swap(y(1),y(ilo))
	     call swap(p(1,:),p(ilo,:))
	     RETURN
	  end if
	  if (iter >= ITMAX) call nrerror('ITMAX exceeded in amoeba')
	  ytry=amotry(-1.0_sp)
	  iter=iter+1
	  if (ytry <= y(ilo)) then
	     ytry=amotry(2.0_sp)
	     iter=iter+1
	  else if (ytry >= y(inhi)) then
	     ysave=y(ihi)
	     ytry=amotry(0.5_sp)
	     iter=iter+1
	     if (ytry >= ysave) then
		p(:,:)=0.5_sp*(p(:,:)+spread(p(ilo,:),1,size(p,1)))
		do i=1,ndim+1
		   if (i /= ilo) y(i)=func(p(i,:))
		end do
		iter=iter+ndim
		psum(:)=sum(p(:,:),dim=1)
	     end if
	  end if
	end do
	END SUBROUTINE amoeba_private
!BL
	FUNCTION amotry(fac)
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: fac
	REAL(SP) :: amotry
	REAL(SP) :: fac1,fac2,ytry
	REAL(SP), DIMENSION(size(p,2)) :: ptry
	fac1=(1.0_sp-fac)/ndim
	fac2=fac1-fac
	ptry(:)=psum(:)*fac1-p(ihi,:)*fac2
	ytry=func(ptry)
	if (ytry < y(ihi)) then
		y(ihi)=ytry
		psum(:)=psum(:)-p(ihi,:)+ptry(:)
		p(ihi,:)=ptry(:)
	end if
	amotry=ytry
	END FUNCTION amotry
	END SUBROUTINE amoeba_r

	SUBROUTINE amoeba_d(p,y,ftol,func,iter)
	USE nrtype; USE nrutil, ONLY : assert_eq,imaxloc,iminloc,nrerror,swap
	IMPLICIT NONE
	INTEGER(I4B), INTENT(OUT) :: iter
	REAL(DP), INTENT(IN) :: ftol
	REAL(DP), DIMENSION(:), INTENT(INOUT) :: y
	REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: p
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), DIMENSION(:), INTENT(IN) :: x
		REAL(DP) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: ITMAX=5000
	REAL(DP), PARAMETER :: TINY=1.0d-10
	INTEGER(I4B) :: ihi,ndim
	REAL(DP), DIMENSION(size(p,2)) :: psum
	call amoeba_private
	CONTAINS
!BL
	SUBROUTINE amoeba_private
	IMPLICIT NONE
	INTEGER(I4B) :: i,ilo,inhi
	REAL(DP) :: rtol,ysave,ytry,ytmp
	ndim=assert_eq(size(p,2),size(p,1)-1,size(y)-1,'amoeba')
	iter=0
	psum(:)=sum(p(:,:),dim=1)
	do
	  ilo=iminloc(y(:))
	  ihi=imaxloc(y(:))
	  ytmp=y(ihi)
	  y(ihi)=y(ilo)
	  inhi=imaxloc(y(:))
	  y(ihi)=ytmp
	  rtol=2.0_dp*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo))+TINY)
	  if (rtol < ftol) then
	     call swap(y(1),y(ilo))
	     call swap(p(1,:),p(ilo,:))
	     RETURN
	  end if
	  if (iter >= ITMAX) call nrerror('ITMAX exceeded in amoeba')
	  ytry=amotry(-1.0_dp)
	  iter=iter+1
	  if (ytry <= y(ilo)) then
	     ytry=amotry(2.0_dp)
	     iter=iter+1
	  else if (ytry >= y(inhi)) then
	     ysave=y(ihi)
	     ytry=amotry(0.5_dp)
	     iter=iter+1
	     if (ytry >= ysave) then
		p(:,:)=0.5_dp*(p(:,:)+spread(p(ilo,:),1,size(p,1)))
		do i=1,ndim+1
		   if (i /= ilo) y(i)=func(p(i,:))
		end do
		iter=iter+ndim
		psum(:)=sum(p(:,:),dim=1)
	     end if
	  end if
	end do
	END SUBROUTINE amoeba_private
!BL
	FUNCTION amotry(fac)
	IMPLICIT NONE
	REAL(DP), INTENT(IN) :: fac
	REAL(DP) :: amotry
	REAL(DP) :: fac1,fac2,ytry
	REAL(DP), DIMENSION(size(p,2)) :: ptry
	fac1=(1.0_dp-fac)/ndim
	fac2=fac1-fac
	ptry(:)=psum(:)*fac1-p(ihi,:)*fac2
	ytry=func(ptry)
	if (ytry < y(ihi)) then
		y(ihi)=ytry
		psum(:)=psum(:)-p(ihi,:)+ptry(:)
		p(ihi,:)=ptry(:)
	end if
	amotry=ytry
	END FUNCTION amotry
	END SUBROUTINE amoeba_d

	FUNCTION brent_r(ax,bx,cx,func,tol,xmin)
	USE nrtype; USE nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: ax,bx,cx,tol
	REAL(SP), INTENT(OUT) :: xmin
	REAL(SP) :: brent_r
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: ITMAX=100
	REAL(SP), PARAMETER :: CGOLD=0.3819660_sp,ZEPS=1.0e-3_sp*epsilon(ax)
	INTEGER(I4B) :: iter
	REAL(SP) :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
	a=min(ax,cx)
	b=max(ax,cx)
	v=bx
	w=v
	x=v
	e=0.0
	fx=func(x)
	fv=fx
	fw=fx
	do iter=1,ITMAX
		xm=0.5_sp*(a+b)
		tol1=tol*abs(x)+ZEPS
		tol2=2.0_sp*tol1
		if (abs(x-xm) <= (tol2-0.5_sp*(b-a))) then
			xmin=x
			brent_r=fx
			RETURN
		end if
		if (abs(e) > tol1) then
			r=(x-w)*(fx-fv)
			q=(x-v)*(fx-fw)
			p=(x-v)*q-(x-w)*r
			q=2.0_sp*(q-r)
			if (q > 0.0) p=-p
			q=abs(q)
			etemp=e
			e=d
			if (abs(p) >= abs(0.5_sp*q*etemp) .or. &
				p <= q*(a-x) .or. p >= q*(b-x)) then
				e=merge(a-x,b-x, x >= xm )
				d=CGOLD*e
			else
				d=p/q
				u=x+d
				if (u-a < tol2 .or. b-u < tol2) d=sign(tol1,xm-x)
			end if
		else
			e=merge(a-x,b-x, x >= xm )
			d=CGOLD*e
		end if
		u=merge(x+d,x+sign(tol1,d), abs(d) >= tol1 )
		fu=func(u)
		if (fu <= fx) then
			if (u >= x) then
				a=x
			else
				b=x
			end if
			call shft(v,w,x,u)
			call shft(fv,fw,fx,fu)
		else
			if (u < x) then
				a=u
			else
				b=u
			end if
			if (fu <= fw .or. w == x) then
				v=w
				fv=fw
				w=u
				fw=fu
			else if (fu <= fv .or. v == x .or. v == w) then
				v=u
				fv=fu
			end if
		end if
	end do
	call nrerror('brent: exceed maximum iterations')
	CONTAINS
!BL
	SUBROUTINE shft(a,b,c,d)
	REAL(SP), INTENT(OUT) :: a
	REAL(SP), INTENT(INOUT) :: b,c
	REAL(SP), INTENT(IN) :: d
	a=b
	b=c
	c=d
	END SUBROUTINE shft
	END FUNCTION brent_r

	FUNCTION brent_d(fbx,ax,bx,cx,func,tol,xmin)
	USE nrtype; USE nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(DP), INTENT(IN) :: ax,bx,cx,tol,fbx
	REAL(DP), INTENT(OUT) :: xmin
	REAL(DP) :: brent_d
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), INTENT(IN) :: x
		REAL(DP) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: ITMAX=100
	REAL(DP), PARAMETER :: CGOLD=0.3819660_dp,ZEPS=1.0e-3_dp*epsilon(ax)
	INTEGER(I4B) :: iter
	REAL(DP) :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
	a=min(ax,cx)
	b=max(ax,cx)
	v=bx
	w=v
	x=v
	e=0.0
	fx=fbx
	fv=fx
	fw=fx
	do iter=1,ITMAX
		xm=0.5_dp*(a+b)
		tol1=tol*abs(x)+ZEPS
		tol2=2.0_dp*tol1
		if (abs(x-xm) <= (tol2-0.5_dp*(b-a))) then
			xmin=x
			brent_d=fx
			RETURN
		end if
		if (abs(e) > tol1) then
			r=(x-w)*(fx-fv)
			q=(x-v)*(fx-fw)
			p=(x-v)*q-(x-w)*r
			q=2.0_dp*(q-r)
			if (q > 0.0) p=-p
			q=abs(q)
			etemp=e
			e=d
			if (abs(p) >= abs(0.5_dp*q*etemp) .or. &
				p <= q*(a-x) .or. p >= q*(b-x)) then
				e=merge(a-x,b-x, x >= xm )
				d=CGOLD*e
			else
				d=p/q
				u=x+d
				if (u-a < tol2 .or. b-u < tol2) d=sign(tol1,xm-x)
			end if
		else
			e=merge(a-x,b-x, x >= xm )
			d=CGOLD*e
		end if
		u=merge(x+d,x+sign(tol1,d), abs(d) >= tol1 )
		fu=func(u)
		if (fu <= fx) then
			if (u >= x) then
				a=x
			else
				b=x
			end if
			call shft(v,w,x,u)
			call shft(fv,fw,fx,fu)
		else
			if (u < x) then
				a=u
			else
				b=u
			end if
			if (fu <= fw .or. w == x) then
				v=w
				fv=fw
				w=u
				fw=fu
			else if (fu <= fv .or. v == x .or. v == w) then
				v=u
				fv=fu
			end if
		end if
	end do
	call nrerror('brent: exceed maximum iterations')
	CONTAINS
!BL
	SUBROUTINE shft(a,b,c,d)
	REAL(DP), INTENT(OUT) :: a
	REAL(DP), INTENT(INOUT) :: b,c
	REAL(DP), INTENT(IN) :: d
	a=b
	b=c
	c=d
	END SUBROUTINE shft
	END FUNCTION brent_d

	SUBROUTINE mnbrak_r(ax,bx,cx,fa,fb,fc,func)
	USE nrtype; USE nrutil, ONLY : swap
	IMPLICIT NONE
	REAL(SP), INTENT(INOUT) :: ax,bx
	REAL(SP), INTENT(OUT) :: cx,fa,fb,fc
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: func
		END FUNCTION func
	END INTERFACE
	REAL(SP), PARAMETER :: GOLD=1.618034_sp,GLIMIT=100.0_sp,TINY=1.0e-20_sp
	REAL(SP) :: fu,q,r,u,ulim
	fa=func(ax)
	fb=func(bx)
	if (fb > fa) then
		call swap(ax,bx)
		call swap(fa,fb)
	end if
	cx=bx+GOLD*(bx-ax)
	fc=func(cx)
	do
		if (fb < fc) RETURN
		r=(bx-ax)*(fb-fc)
		q=(bx-cx)*(fb-fa)
		u=bx-((bx-cx)*q-(bx-ax)*r)/(2.0_sp*sign(max(abs(q-r),TINY),q-r))
		ulim=bx+GLIMIT*(cx-bx)
		if ((bx-u)*(u-cx) > 0.0) then
			fu=func(u)
			if (fu < fc) then
				ax=bx
				fa=fb
				bx=u
				fb=fu
				RETURN
			else if (fu > fb) then
				cx=u
				fc=fu
				RETURN
			end if
			u=cx+GOLD*(cx-bx)
			fu=func(u)
		else if ((cx-u)*(u-ulim) > 0.0) then
			fu=func(u)
			if (fu < fc) then
				bx=cx
				cx=u
				u=cx+GOLD*(cx-bx)
				call shft(fb,fc,fu,func(u))
			end if
		else if ((u-ulim)*(ulim-cx) >= 0.0) then
			u=ulim
			fu=func(u)
		else
			u=cx+GOLD*(cx-bx)
			fu=func(u)
		end if
		call shft(ax,bx,cx,u)
		call shft(fa,fb,fc,fu)
	end do
	CONTAINS
!BL
	SUBROUTINE shft(a,b,c,d)
	REAL(SP), INTENT(OUT) :: a
	REAL(SP), INTENT(INOUT) :: b,c
	REAL(SP), INTENT(IN) :: d
	a=b
	b=c
	c=d
	END SUBROUTINE shft
	END SUBROUTINE mnbrak_r

	SUBROUTINE mnbrak_d(ax,bx,cx,fa,fb,fc,func)
	USE nrtype; USE nrutil, ONLY : swap
	IMPLICIT NONE
	REAL(DP), INTENT(INOUT) :: ax,bx
	REAL(DP), INTENT(OUT) :: cx,fa,fb,fc
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), INTENT(IN) :: x
		REAL(DP) :: func
		END FUNCTION func
	END INTERFACE
	REAL(DP), PARAMETER :: GOLD=1.618034_dp,GLIMIT=100.0_dp,TINY=1.0e-20_dp
	REAL(DP) :: fu,q,r,u,ulim
	fa=func(ax)
	fb=func(bx)
	if (fb > fa) then
		call swap(ax,bx)
		call swap(fa,fb)
	end if
	cx=bx+GOLD*(bx-ax)
	fc=func(cx)
	do
		if (fb < fc) RETURN
		r=(bx-ax)*(fb-fc)
		q=(bx-cx)*(fb-fa)
		u=bx-((bx-cx)*q-(bx-ax)*r)/(2.0_dp*sign(max(abs(q-r),TINY),q-r))
		ulim=bx+GLIMIT*(cx-bx)
		if ((bx-u)*(u-cx) > 0.0) then
			fu=func(u)
			if (fu < fc) then
				ax=bx
				fa=fb
				bx=u
				fb=fu
				RETURN
			else if (fu > fb) then
				cx=u
				fc=fu
				RETURN
			end if
			u=cx+GOLD*(cx-bx)
			fu=func(u)
		else if ((cx-u)*(u-ulim) > 0.0) then
			fu=func(u)
			if (fu < fc) then
				bx=cx
				cx=u
				u=cx+GOLD*(cx-bx)
				call shft(fb,fc,fu,func(u))
			end if
		else if ((u-ulim)*(ulim-cx) >= 0.0) then
			u=ulim
			fu=func(u)
		else
			u=cx+GOLD*(cx-bx)
			fu=func(u)
		end if
		call shft(ax,bx,cx,u)
		call shft(fa,fb,fc,fu)
	end do
	CONTAINS
!BL
	SUBROUTINE shft(a,b,c,d)
	REAL(DP), INTENT(OUT) :: a
	REAL(DP), INTENT(INOUT) :: b,c
	REAL(DP), INTENT(IN) :: d
	a=b
	b=c
	c=d
	END SUBROUTINE shft
	END SUBROUTINE mnbrak_d

END MODULE nr_minimization
