      SUBROUTINE thermcd(NCH2,TL,TR,scaleM,KQQ)
       

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      REAL*8 hbar,AKB,PI,Te,Tw,Ef_L,Ef_R,YP1,YPN,AVGI
      REAL*8 ABSERR,RELERR,AXSI1,AXSI2,Twmax,Twmin
      REAL*8 temp1,temp2,temp3,temp4,WTOT
      REAL*8 VT,VL,RHO,Ymdl,dcross,dlong,CR,G1,G2,G3
      COMPLEX*16 ctemp1,ctemp2,ctemp3,ctemp4,CAJ
      LOGICAL alive,scaleH
      REAL*8, EXTERNAL:: F_R,F_L,FCR ,FCNPOWER2,FCNPOWERnod
      REAL*8  vTw(2000),vPTOE(2000)  !store Power dissipated to electrode
   
      COMMON/SURF/VL,VT,CR
      COMMON/TVSPOWER/vTw,vPTOE
	  REAL*8 TRR1,Cth,KK,CC,KB,TL,TR,KQQ

!Introduce scaling factor beta=0.78(A)^-1 from Mark Reed PRB 68 035416 (2003)
       betaM=0.78d0*0.529d0
  !     ddd1=9.4162d0
   !    ddd2=9.4162d0+DFLOAT(NCH2-1)*(((23.23d0-2.74d0)/0.529d0)/16.d0)
  !        scaleM=DEXP(-1.d0*betaM*(ddd2-ddd1))
!We scale I->scalef*I, Heating from WTOT->scalef*WTOT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!rescale the conductance to 1.1G0
!	  scaleI= 1.d0
!     WRITE(*,*)'Ef_L=',Ef_L,Ef_R,NCH2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      AKB=(8.6178D-5)/(13.6D0) !BOLTZMANN CONST.
	hbar=1.d0
      PI=ACOS(-1.D0)  !PI=3.1415926
!	write(*,*) sqrt(4.0d0*21.40d0/PI) 
!	pause

!!!!!!!!!!!!INPUT for Geller's paper: (heat dissipate to ETECTRODES)
       VT=1.2d5   ! FOR AU
       VL=3.2d5   ! FOR AU 

 !	  VT=3.1d5    !FOR AL 
 !       VL=6.35d5   !FOR AL
       
        

      RHO=19.3d0 !for AU 
  !  RHO=2.7d0  !For AL

        
       Ymdl=5.260d12 !from the total energy caculation 
 !     Ymdl=2.30d12 
 !     dcross=6.828d0 
 !     dcross=5.2198990d0 ! FOR C
	 dcross=8.00d0   
        
  !      dlong=21.59960d0+(NCH2-4)*2.830d0
	 !9.4162d0  !+DFLOAT(NCH2-1)*(((10.42d0-2.74d0)/0.529d0)/6.d0)
  !     scaleM=exp(-0.950d0*(NCH2-4))  
  !     scaleM=DEXP(-1.d0*betaM*(NCH2-4)*2.830d0)
!	 dlong=10.80d0  !FOR 2AL 
!	 dlong=16.5680d0  !FOR 3AL 
!	 dlong=23.940d0 ! FOR 4AL 
       dlong=16.40d0 !FOR NO 
!	 dlong=17.20d0 !FOR NH2
!Twrite(*,*)'dlong=',dlong*0.529
!Get Parameters use by Geller's paper:
!FOR  SURFACE PHONON
           DHBAR=1.05459d-27
           ABSERR=1.D-16
           RELERR=1.D-16
           MAXFCT=1000
           ITAPE=1005
           AXSI1=0.D0
           AXSI2=1.D0

       
   


      CALL ZEROIN (FCR, ABSERR, RELERR, MAXFCT, ITAPE, AXSI1, AXSI2,
     * FB,NUMFCT, IERR) 
           CR=AXSI2*VT
      CALL SURFPH(G1,G2,G3)
         ! TRR1=0.0d0 
	   ! DO WHILE(TRR1.LT.50)
!	    TE=TRR1
!	    Twin=(TE+TRR1)/2.0d0
          Cth=(G1*PI/CR**3+G2/VL**3+G3/(2.d0*VT**3))/(4.d0*PI*PI*RHO)
     */DHBAR
          KB=1.38065D-16
      !   Cth=1.3D8 
!		write(*,*) Cth 
!		pause 
        !KK=PI*dcross**2/(4.0d0*dlong)*0.529d-8*Ymdl   !For C  !!!
	KK=21717.05263  !2Benzene
!	   KK=PI*dcross**2/(4.0d0*dlong)*0.529d-8*Ymdl !FOR C 
  !      write(*,*) PI*dcross**2/(4.0d0)
!	  pause 
!	KK=1.930d0*1.602190d0*1E-12/(0.529**2)*1E16 !FOR 2 AL 
   !   KK=1.760d0*1.602190d0*1E-12/(0.529**2)*1E16 !FOR 3 AL 
   !   KK=1.590d0*1.602190d0*1E-12/(0.529**2)*1E16 !FOR 4 AL 
	!	write(*,*) KK 
!		pause
		!KK=1.0D4
       ! CC=8.0d0*PI**5*KK**2*Cth**2*KB**4*1.D-7/15.0d0/DHBAR 
!	  PRINT*, KK, Cth
!        PRINT*, KB, DHBAR
!	PRINT*, PI
!       KQQ=8.0d0*PI**5*KK**2*KB**4*Cth**2/15.0d0/DHBAR
!     **(TL**3+TR**3)/2.0d0*1.0E-7  !J S-1 K-1 
 




       end 

      FUNCTION FCR(X)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       IMPLICIT INTEGER (I-N)
       COMMON/SURF/VL,VT,CR
        PI=ACOS(-1.D0)
           V=VT/VL
           FCR=X**6-8.D0*X**4+8.D0*(3.D0-2.D0*V*V)*X*X-16.D0*(1.D0-V*V)
       RETURN
       END


      INTEGER FUNCTION MACHPD (X)
      DOUBLEPRECISION X
      MACHPD = 0
      IF (1.0D0.LT.X) MACHPD = 1
      RETURN
      END FUNCTION MACHPD