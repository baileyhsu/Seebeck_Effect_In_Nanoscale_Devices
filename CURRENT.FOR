      SUBROUTINE CURRENT(L,TL,TR,DT,x1,x4,CC,DD,EG,HG,HH,SEE,KQ)
!     USE IMSL
	implicit none
      
	logical alive
	integer i,j,L,n
      integer error,count 
	real*8 EG,HG,EG1,HG1,SG,FERMIR,FERMIL,DT, yp1,ypn,x,AA,BB,U,HH
      real*8 K1UL,K1UR,K2UL,K2UR,KQ,ETT,EGG,EG11
      parameter(n=200) 
      
  
      real*8 x1(L),x2(L),x3(L),x4(L) ,x5(L),x6(L),x7(L),y2(L),CC,DD
     *,SG1
      

  !   common /XX/x1(L), x2(L),x3(L),x4(L),x5(L),x6(L),x7(L)

      Real*8   y(n)
      Real*8   w(n) 

	real*8 VL,VR,KT1,KT2,VBIAS
      common /VL/VL/VR/VR/KT1/KT1/KT2/KT2/VBIAS/VBIAS
      

      real*8 SEE,T,TL,TR

  

    !   U=VL
 !      do while (DT.LT.10.0d0)

 !!    NOW WE DEFINE KT1 AND KT2 (K=BOLTZMANN EV/K, 1 RYDBERG=13.6058 EV)
 !!    INITIALIZE ALL VARIABLES
 !!    L= ENERGY CUT POINTS, N= GAUSSLEG INTEGRATION POINTS

        
       KT1=8.6170d0*1E-5*(TL)/13.60580d0
	 KT2=8.6170d0*1E-5*(TR)/13.60580d0

	EG=0.0d0 
	HG=0.0d0 
	CC=0.0d0
	DD=0.0d0
	HH=0.0d0 
      SEE=0.0d0 
	ETT=0.0d0
	EG11=0.0d0  

 !    FIRST WE GO OVER ALL ENERGY POINTS WHICH IS L LONG AND FIND THE 2ND DERIVATIVE Y2 USING SPLINE

	do i=1,L-1
	
 !    HERE WE CALL GAUSSLEG INTEGRATION WHERE BETWEEN TWO ENERGY POINTS, WE SELECT "n" POINTS AND GAUSSLEG
 !    RETURNS THE ABSCISSA Y(N) and Weight(W(N))

	CALL GAULEG(x1(i),x1(i+1),y,w,n) 
 
 !    HERE WE CALL SPLINE TO DO INTERPOLATION, X1, X4, L CORRESPOND TO THE ORIGINAL DATA
 !    WHERE YP1 AND YPN AND Y2 
 
 !    SPINLESS
      CALL SPLINE(x1,x4,L,yp1,ypn,y2)
       
 
 !    Y(J) IS THE ABCISSA, THE 200 ENERGY CUT POINTS, AND THE X CORRESPONDS TO 
 !    INTERPOLATED CURRENT (TRANSMISSION)
             
      do j=1,n
      call SPLINT(x1,x4,y2,L,y(j),x)
       
	AA=(FERMIR(y(j))-FERMIL(y(j)))
      EG=EG+x*w(j)*AA
	EG11=EG11+x*w(j)*AA*(VR-VL)
	HG=HG+x*w(j)*AA*(-VL+y(j))

!     THE CC AND DD TERMS ARE USED TO CALCULATE THE SEEBECK TERM
!     CC ACCOUNTS FOR THE K1/T WHERE K1=- INT( (E-EFL)*F(1-F)/(KB*T)*X (SUMMING LEFT AND RIGHT)
!     DD ACCOUNTS FOR THE K0= -INT ( F(1-F)/(KB*T) )*X (SUMMING LEFT AND RIGHT)

      CC=CC+w(j)*x*
     *(1.0d0/TL*FERMIL(y(j))*(1.0d0-FERMIL(y(j)))/(8.6170d0*1E-5*(TL))
     **(-VL+y(j))
     *+1.0d0/TR*FERMIR(y(j))*(1.0d0-FERMIR(y(j)))/(8.6170d0*1E-5*(TR))
     **(-VR+y(j)))
	DD=DD+w(j)*x*(FERMIR(y(j))*(1.0d0-FERMIR(y(j)))
     */(8.6170d0*1E-5*(TR))
     *+FERMIL(y(j))*(1.0d0-FERMIL(y(j)))/(8.6170d0*1E-5*(TL)))


	HH=HH+0.50d0*w(j)*x*(FERMIL(y(j))*(1.0d0-FERMIL(y(j)))
     */(8.6170d0*1E-5*(TL))+FERMIR(y(j))*(1.0d0-FERMIR(y(j)))
     */(8.6170d0*1E-5*(TR)))
     */151.0d0
	
	
      
	end do 
	end do  
      SEE=-(CC)/(DD)*13.60580d0*1E6
	ETT=EG11/HG 	

      EG=EG/151.0d0 
       
      K1UL=0.0d0 
	K1UR=0.0d0
	K2UL=0.0d0 
	K2UR=0.0d0

	 
      do i=1,L-1

	call GAULEG(x1(i),x1(i+1),y,w,n) 
      call  SPLINE(x1,x4,L,yp1,ypn,y2)
      
!______________________________________________________________________
!         CONVERT EVERYTHING TO EV            
       
      do j=1,n
      call SPLINT(x1,x4,y2,L,y(j),x)
       
      K1UL=K1UL+w(j)*x*(FERMIL(y(j))*(1.0d0-FERMIL(y(j))))
     **(13.60580d0)
     */(8.6170d0*1E-5*(TL))*(-VL+y(j))/151.0d0
	K1UR=K1UR+w(j)*x*(FERMIR(y(j))*(1.0d0-FERMIR(y(j))))
     **(13.60580d0)
     */(8.6170d0*1E-5*(TR))*(-VR+y(j))/151.0d0

	K2UL=K2UL+w(j)*x*(FERMIL(y(j))*(1.0d0-FERMIL(y(j))))
     */(8.6170d0*1E-5*(TL))*(-VL+y(j))**2*(13.60580d0)**2/151.0d0
	K2UR=K2UR+w(j)*x*(FERMIR(y(j))*(1.0d0-FERMIR(y(j))))
     */(8.6170d0*1E-5*(TR))*(-VR+y(j))**2*(13.60580d0)**2/151.0d0

      
	end do 
      
      

       
	end do  
   

       
       KQ=2.0d0*(((K1UL+K1UR))*SEE*1E-6+K2UL/TL+K2UR/TR)



	END 
