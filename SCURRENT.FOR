       SUBROUTINE SCURRENT(L,TL,TR,DT,x1,x4up,x4dn,SEE,SEE1,SEE2,KQ)
!      USE IMSL
	IMPLICIT NONE     
	logical alive
	integer i,j,L,n
      integer error,count 
	real*8 EG,HG,EG1,HG1,SG,FERMIR,FERMIL,DT, yp1,ypn,xup,xdn,AA,BB,U
	real*8 K1ULdn,K1URdn,K2ULdn,K2URdn
	real*8 K1ULup,K1URup,K2ULup,K2URup,KQ
      parameter(n=300)  
      real*8 x1(L),x2(L),x3(L),x4up(L),x4dn(L),x5(L),x6(L),x7(L)
	real*8 y2up(L),y2dn(L),CC,DD,CC1,DD1, SG1

  !   common /XX/x1(L), x2(L),x3(L),x4(L),x5(L),x6(L),x7(L)

      Real*8   y(n)
      Real*8   w(n) 
	real*8 VL,VR,KT1,KT2,VBIAS
      common /VL/VL/VR/VR/KT1/KT1/KT2/KT2/VBIAS/VBIAS
      

      real*8 SEE,SEE1,SEE2,T,TL,TR
       

 !!    NOW WE DEFINE KT1 AND KT2 (K=BOLTZMANN EV/K, 1 RYDBERG=13.6058 EV)
 !!    INITIALIZE ALL VARIABLES
 !!    L= ENERGY CUT POINTS, N= GAUSSLEG INTEGRATION POINTS
  

      KT1=8.6170d0*1E-5*(TL)/13.60580d0
      KT2=8.6170d0*1E-5*(TR)/13.60580d0
      CC=0.0d0
      DD=0.0d0
      CC1=0.0D0
      DD1=0.0D0
      SEE=0.0d0 
      SEE1=0.0D0
     
      
      CALL GaussLeg(x1(1),x1(L),n,y,w)  
        
        CALL SPLINE(x1,x4up,L,yp1,ypn,y2up)
        CALL SPLINE(x1,x4dn,L,yp1,ypn,y2dn)

        do j=1,n
      call SPLINT(x1,x4up,y2up,L,y(j),xup)
      call SPLINT(x1,x4dn,y2dn,L,y(j),xdn)
 
      CC=CC+w(j)*xup*
     *(1.0d0/TL*FERMIL(y(j))*(1.0d0-FERMIL(y(j)))/(8.6170d0*1E-5*(TL))
     **(-VL+y(j))
     *+1.0d0/TR*FERMIR(y(j))*(1.0d0-FERMIR(y(j)))/(8.6170d0*1E-5*(TR))
     **(-VR+y(j)))

      CC1=CC1+w(j)*xdn*
     *(1.0d0/TL*FERMIL(y(j))*(1.0d0-FERMIL(y(j)))/(8.6170d0*1E-5*(TL))
     **(-VL+y(j))
     *+1.0d0/TR*FERMIR(y(j))*(1.0d0-FERMIR(y(j)))/(8.6170d0*1E-5*(TR))
     **(-VR+y(j)))

      DD=DD+w(j)*xup*(FERMIR(y(j))*(1.0d0-FERMIR(y(j)))
     */(8.6170d0*1E-5*(TR))
     *+FERMIL(y(j))*(1.0d0-FERMIL(y(j)))/(8.6170d0*1E-5*(TL)))
      
      DD1=DD1+w(j)*xdn*(FERMIR(y(j))*(1.0d0-FERMIR(y(j)))
     */(8.6170d0*1E-5*(TR))
     *+FERMIL(y(j))*(1.0d0-FERMIL(y(j)))/(8.6170d0*1E-5*(TL)))	
      
      end do
       
      SEE=-(CC)/(DD)*13.60580d0*1E6
      SEE1=-(CC1)/(DD1)*13.60580D0*1E6
      SEE2=-(CC-CC1)/(DD+DD1)*13.60580D0*1E6
       
      K1ULup=0.0d0 
	K1URup=0.0d0
	K2ULup=0.0d0 
	K2URup=0.0d0
      K1ULdn=0.0d0 
	K1URdn=0.0d0
	K2ULdn=0.0d0 
	K2URdn=0.0d0
	 
      do i=1,L-1

	call GAULEG(x1(i),x1(i+1),y,w,n) 

      call  SPLINE(x1,x4up,L,yp1,ypn,y2up)
      call  SPLINE(x1,x4dn,L,yp1,ypn,y2dn)
!______________________________________________________________________
!         CONVERT EVERYTHING TO EV            
       
      do j=1,n
      call SPLINT(x1,x4up,y2up,L,y(j),xup)
      call SPLINT(x1,x4dn,y2dn,L,y(j),xdn)
       
      K1ULup=K1ULup+w(j)*xup*(FERMIL(y(j))*(1.0d0-FERMIL(y(j))))
     **(13.60580d0)
     */(8.6170d0*1E-5*(TL))*(-VL+y(j))/151.0d0
	K1URup=K1URup+w(j)*xup*(FERMIR(y(j))*(1.0d0-FERMIR(y(j))))
     **(13.60580d0)
     */(8.6170d0*1E-5*(TR))*(-VR+y(j))/151.0d0

      K1ULdn=K1ULdn+w(j)*xdn*(FERMIL(y(j))*(1.0d0-FERMIL(y(j))))
     **(13.60580d0)
     */(8.6170d0*1E-5*(TL))*(-VL+y(j))/151.0d0
	K1URdn=K1URdn+w(j)*xdn*(FERMIR(y(j))*(1.0d0-FERMIR(y(j))))
     **(13.60580d0)
     */(8.6170d0*1E-5*(TR))*(-VR+y(j))/151.0d0

	K2ULup=K2ULup+w(j)*xup*(FERMIL(y(j))*(1.0d0-FERMIL(y(j))))
     */(8.6170d0*1E-5*(TL))*(-VL+y(j))**2*(13.60580d0)**2/151.0d0
	K2URup=K2URup+w(j)*xup*(FERMIR(y(j))*(1.0d0-FERMIR(y(j))))
     */(8.6170d0*1E-5*(TR))*(-VR+y(j))**2*(13.60580d0)**2/151.0d0
      
      K2ULdn=K2ULdn+w(j)*xdn*(FERMIL(y(j))*(1.0d0-FERMIL(y(j))))
     */(8.6170d0*1E-5*(TL))*(-VL+y(j))**2*(13.60580d0)**2/151.0d0
	K2URdn=K2URdn+w(j)*xdn*(FERMIR(y(j))*(1.0d0-FERMIR(y(j))))
     */(8.6170d0*1E-5*(TR))*(-VR+y(j))**2*(13.60580d0)**2/151.0d0

	END DO 
            
	END DO  
  
        
       KQ=2.0d0*(((K1ULup+K1URup-K1ULdn-K1URdn))*SEE*1E-6
     *	+K2ULup/TL+K2URup/TR-K2ULdn/TL-K2URdn/TR)

      END

 
