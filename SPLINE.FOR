!Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e., yi = f(xi), with
!x1 < x2 < . . . < xN , and given values yp1 and ypn for the first derivative of the interpolating function 
!at points 1 and n, respectively, this routine returns an array y2(1:n) of
!length n which contains the second derivatives of the interpolating function at the tabulated
!points xi. If yp1 and/or ypn are equal to 1 กั 1030 or larger, the routine is signaled to set
!the corresponding boundary condition for a natural spline, with zero second derivative on
!that boundary.
!Parameter: NMAX is the largest anticipated value of n.


      SUBROUTINE SPLINE(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      DOUBLE PRECISION yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=2500)
      INTEGER i,k
      DOUBLE PRECISION p,qn,sig,un,u(NMAX)
      if (yp1.gt..99d30) then
        y2(1)=0.d0
        u(1)=0.d0
      else
        y2(1)=-0.5d0
        u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.d0
        y2(i)=(sig-1.d0)/p
        u(i)=(6.d0*((y(i+1)-y(i))/(x(i+ 
     *1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*  
     *u(i-1))/p
       end do
      if (ypn.gt..99d30) then
        qn=0.d0
        un=0.d0
      else
        qn=0.5d0
        un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
      do k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
      end do
      return
      END
