      SUBROUTINE CURRENT_EP(L,TL,TR,DT,fr,x1,x4,CC,DD,EG,HG,HH,SEE,KQ,
     *CLR,CRR,Tw,SEE1)

!      USE IMSL
	implicit none
      
	logical alive
	integer i,j,k
      integer error,count,NMODE 
	real*8 EG,HG,EG1,HG1,SG,FERMIR,FERMIL,DT, yp1,ypn,x,AA,BB,U,HH,Tw
      real*8 K1UL,K1UR,K2UL,K2UR,KQ,ETT,EGG,EG11,Temp,AKB,BMAX,TEMPL
	real*8 K0LL,K0LR
	real*8 TEMPR,AVG1,fre,CLRc,CRRc
      integer L
 !     parameter(L=116) 
      integer n
      parameter(n=20, NMODE=3, AKB=8.617d-5)   
      real*8 x1(L),x2(L),x3(L),x4(L) ,x5(L),x6(L),x7(L),y2(L),CC,DD,CC1
      real*8 DD1,CC1F,DD1F
	real*8 fr(NMODE),CLR(NMODE),CRR(NMODE),AVG(NMODE)
     *,SG1
  !    common /XX/x1(L), x2(L),x3(L),x4(L),x5(L),x6(L),x7(L)
      Real*8   y(n)
      Real*8   w(n) 
	real*8 VL,VR,KT1,KT2,VBIAS
      common /VL/VL/VR/VR/KT1/KT1/KT2/KT2/VBIAS/VBIAS
      real*8 SEE,T,TL,TR,SEE1

        
       KT1=8.6170d0*1E-5*(TL)/13.6058
	 KT2=8.6170d0*1E-5*(TR)/13.6058      
  !    write(*,*) 	-0.2544352939d0-8.1670d0*1E-5*(400.0d0)*5.0d0/13.60580
!	  pause
	EG=0.0d0 
	HG=0.0d0 

	HH=0.0d0 
      SEE=0.0d0 
	ETT=0.0d0
	EG11=0.0d0
      CLRc=0.0d0
	CRRc=0.0d0
      CC1F=0.d0
      DD1F=0.d0
      SEE1=0.d0 
      

       do k=1,NMODE
        IF(Tw.ge.1.d-2)THEN                          
        TEMP=fr(k)/(AKB*Tw)
        BMAX=MIN(TEMP,400.d0)
        AVG(k)=1.d0/(DEXP(BMAX)-1.d0) 
        ELSE
        AVG(k)=0.d0
        END IF      
        
	  AVG1=AVG(k)
	  fre=fr(k)
        CLRc=CLR(k)
	  CRRc=CRR(k)
       
        CC1=0.d0
	  DD1=0.d0
	  CC=0.d0
	  DD=0.d0

       
      do i=1,L-1
	call GAULEG(x1(i),x1(i+1),y,w,n) 
      call SPLINE(x1,x4,L,yp1,ypn,y2)

      
      do j=1,n
      call SPLINT(x1,x4,y2,L,y(j),x)

      TEMPL=8.6170d0*1E-5*TL/13.6058 
	TEMPR=8.6170d0*1E-5*TR/13.6058 
      
      
	AA=(FERMIR(y(j))-FERMIL(y(j)))

      EG=EG+x*w(j)*AA
	EG11=EG11+x*w(j)*AA*(VR-VL)
	HG=HG+x*w(j)*AA*(-VL+y(j))
      
	CC=CC+w(j)*x*
     *(1.0d0/TL*FERMIL(y(j))*(1.0d0-FERMIL(y(j)))/TEMPL
     **(-VL+y(j))
     *+1.0d0/TR*FERMIR(y(j))*(1.0d0-FERMIR(y(j)))/TEMPR
     **(-VR+y(j)))


      DD=DD+w(j)*x*(FERMIR(y(j))*(1.0d0-FERMIR(y(j)))
     */TEMPR
     *+FERMIL(y(j))*(1.0d0-FERMIL(y(j)))/TEMPL)


      CC1=CC1+w(j)*x*
     *(CLRc*AVG1/TR*FERMIR(y(j))*(1.0d0-FERMIR(y(j)))
     **(y(j)-VR)*(1.0d0-FERMIL(y(j)+fre))/TEMPR	
     *+CLRc*(1.0d0+AVG1)/TR*FERMIR(y(j))*(1.0d0-FERMIR
     *(y(j)))*(y(j)-VR)*(1.0d0-FERMIL(y(j)-fre))/TEMPR
     *-CRRc*AVG1/TR*FERMIR(y(j))*(1.0d0-FERMIR(y(j)))
     **(-VR+y(j))*(1.0d0-FERMIR(y(j)+fre))/TEMPR
     *-CRRc*(1.0d0+AVG1)/TR*FERMIR(y(j))*(1.0d0-FERMIR(y(j)))
     **(y(j)-VR)*(1.0d0-FERMIR(y(j)-fre))/TEMPR
     *+CLRc*AVG1/TL*FERMIL(y(j))*(1.0d0-FERMIL(y(j)))/TEMPL*
     *(y(j)-VL)*(1.0d0-FERMIR(y(j)+fre))
     *+CLRc*(1.0d0+AVG1)/TL*FERMIL(y(j))*(1.0d0-FERMIL(y(j)))/TEMPL*
     *(y(j)-VL)*(1.0d0-FERMIR(y(j)-fre))
     *-CRRc*AVG1/TR*FERMIL(y(j))*(1.0d0-FERMIL(y(j)))
     */TEMPL*(-VL+y(j))*(1.0d0-FERMIL(y(j)+fre))
     *-CRRc*(1.0d0+AVG1)/TR*FERMIL(y(j))*(1.0d0-FERMIL(y(j)))
     */TEMPR*(y(j)-VL)*(1.0d0-FERMIL(y(j)-fre))
     *+CLRc*AVG1*(y(j)-VR+fre)/TR*FERMIR(y(j)+fre)*(1.0d0
     *-FERMIR(y(j)+fre))/TEMPR*FERMIL(y(j))
     *+CLRc*AVG1*(y(j)-VL+fre)/TL*FERMIL(y(j)+fre)*(1.0d0-FERMIL(y(j)
     *+fre))/TEMPL*FERMIR(y(j))   
     *+CLRc*(1.0d0+AVG1)*(y(j)-VR-fre)/TR*FERMIR(y(j)-fre)*
     *(1.0d0-FERMIR(y(j)-fre))/TEMPR*FERMIL(y(j))
     *+CLRc*(1.0d0+AVG1)*(y(j)-VL-fre)/TL*FERMIL(y(j)-fre)*
     *(1.d0-FERMIL(y(j)-fre))/TEMPL*FERMIR(y(j))
     *+CRRc*AVG1*(y(j)-VR+fre)/TR*FERMIR(y(j)+fre)*
     *(1.0d0-FERMIR(y(j)+fre))/TEMPR*FERMIR(y(j))
     *+CRRc*AVG1*(y(j)-VL+fre)/TL*FERMIL(y(j)+fre)*
     *(1.d0-FERMIL(y(j)+fre))/TEMPL*FERMIL(y(j))
     *+CRRc*(1.0d0+AVG1)*(y(j)-VR-fre)/TR*FERMIR(y(j)-fre)*
     *(1.0d0-FERMIR(y(j)-fre))/TEMPR*FERMIR(y(j))
     *+CRRc*(1.0d0+AVG1)*(y(j)-VL-fre)/TL*FERMIL(y(j)-fre)*
     *(1.0d0-FERMIL(y(j)-fre))*FERMIL(y(j))/TEMPL)

     
    
     	DD1=DD1+w(j)*x*
     *(CLRc*AVG1*FERMIR(y(j))*(1.0d0-FERMIR(y(j)))*
     *(1.d0-FERMIL(y(j)+fre))/TEMPR
     *+CLRc*(1+AVG1)*FERMIR(y(j))*(1.0d0-FERMIR(y(j)))*
     *(1-FERMIL(y(j)-fre))/TEMPR
     *-CRRc*AVG1*FERMIR(y(j))*(1.0d0-FERMIR(y(j)))
     **(1-FERMIR(y(j)+fre))/TEMPR
     *-CRRc*(1+AVG1)*FERMIR(y(j))*(1.0d0-FERMIR(y(j)))
     **(1-FERMIR(y(j)-fre))/TEMPR
     *+CLRc*AVG1*FERMIL(y(j))*(1.0d0-FERMIL(y(j)))*
     *(1.0d0-FERMIR(y(j)+fre))/TEMPL
     *+CLRc*(1+AVG1)*FERMIL(y(j))*(1.0d0-FERMIL(y(j)))*
     *(1.0d0-FERMIR(y(j)-fre))/TEMPL
     *-CRRc*AVG1*FERMIL(y(j))*(1.0d0-FERMIL(y(j)))
     **(1-FERMIL(y(j)+fre))/TEMPL
     *-CRRc*(1+AVG1)*FERMIL(y(j))*(1.0d0-FERMIL(y(j)))
     **(1-FERMIL(y(j)-fre))/TEMPR
     *+CLRc*AVG1*FERMIR(y(j)+fre)*(1.d0-FERMIR(y(j)+fre)) 
     **FERMIL(y(j))/TEMPR
     *+CLRc*AVG1*FERMIL(y(j)+fre)*(1.d0-FERMIL(y(j)+fre)) 
     **FERMIR(y(j))/TEMPL   
     *+CLRc*(1.d0+AVG1)*FERMIR(y(j)-fre)*(1.d0-FERMIR(y(j)
     *-fre))*FERMIL(y(j))/TEMPR
     *+CLRc*(1.d0+AVG1)*FERMIL(y(j)-fre)*(1.d0-FERMIL(y(j)
     *-fre))*FERMIR(y(j))/TEMPL
     *+CRRc*AVG1*FERMIR(y(j)+fre)*(1.d0-FERMIR(y(j)+fre)) 
     **FERMIR(y(j))/TEMPR
     *+CRRc*AVG1*FERMIL(y(j)+fre)*(1.d0-FERMIL(y(j)+fre)) 
     **FERMIL(y(j))/TEMPL
     *+CRRc*(1.d0+AVG1)*FERMIR(y(j)-fre)*(1.d0-FERMIR(y(j)
     *-fre))*FERMIR(y(j))/TEMPR
     *+CRRc*(1+AVG1)*FERMIL(y(j)-fre)*(1.d0-FERMIL(y(j)
     *-fre))*FERMIL(y(j))/TEMPL)
       
   
	HH=HH+0.50d0*w(j)*x*(FERMIL(y(j))*(1.0d0-FERMIL(y(j)))
     */(8.6170d0*1E-5*(TL))+FERMIR(y(j))*(1.0d0-FERMIR(y(j)))
     */(8.6170d0*1E-5*(TR)))
     */151.0d0
	
	
      
	
	end do 
	
	end do  
   !  SUMMING ALL MODE CONTRIBUTION   
  
      CC1F=CC1F+CC1
	DD1F=DD1F+DD1 
	 

       
	end do 
	
		
   ! 	PRINT*, CC, CC1F, DD, DD1F  
          
      SEE=-(CC)/(DD)*13.60580d0*1E6
	SEE1=-(CC-CC1F)/(DD-DD1F)*13.6058d0*1E6
	ETT=EG11/HG 	

      EG=EG/151.0d0 
      
	K0LL=0.0d0
	K0LR=0.0d0	 
      K1UL=0.0d0 
	K1UR=0.0d0
	K2UL=0.0d0 
	K2UR=0.0d0
	 
      do i=1,L-1

	call GAULEG(x1(i),x1(i+1),y,w,n) 
      call  SPLINE(x1,x4,L,yp1,ypn,y2)
       
   
       
      do j=1,n
      call SPLINT(x1,x4,y2,L,y(j),x)

      K0LL=K0LL+w(j)*x*(FERMIL(y(j))*(1.0d0-FERMIL(y(j))))
     */(8.6170d0*1E-5*(TL))/151.0d0
	
	K0LR=K0LR+w(j)*x*(FERMIR(y(j))*(1.0d0-FERMIR(y(j))))
     */(8.6170d0*1E-5*(TR))/151.0d0

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
          
       KQ=1.0d0/2.0d0*(((K1UL+K1UR)**2)/(K0LL+K0LR)+K2UL/TL+K2UR/TR)


	END 
