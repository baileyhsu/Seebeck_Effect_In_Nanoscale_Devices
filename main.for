      PROGRAM Seebeck  
      !USE IMSL
	IMPLICIT NONE
      integer L, NMODE
	INTEGER CASES
      PARAMETER(L=254, NMODE=3)
	INTEGER i,j,NCH2
	real*8 EG,HG,EG1,HG1,SG,DT,yp1,ypn,x,AA,BB,U,DI,GG
      real*8 HH,EE,XX,x1(L),x2(L),x3(L),x4(L),x5(L),x6(L),x7(L)
	real*8 x4up(L),x4dn(L)
      real*8 CC,DD,SG1,DE,SEE,SEE1,DIANDAO,ddd,KQ,KEQ,KEQ1,KQQ
	real*8 CLR(NMODE),CRR(NMODE),fr(NMODE)
	real*8 Tw
	real*8 VL,VR,KT1,KT2,VBIAS,T,HH1,diandao1 
      common /VL/VL/VR/VR/KT1/KT1/KT2/KT2/VBIAS/VBIAS     
	real*8 TL,TR,ZT,scaleM,Kphonon,ET,ETT
      
   
   !   THE VOLTAGE UNIT HERE IS RYDBERG, AND SETTING THE VR TO BE THE FERMI LEVEL
   !   AND THE VL IS SHIFTED DOWNWARDS BY A BIAS VOLTAGE
    
       VBIAS=0.010d0
	 VR=-0.2536999997d0
	 VL=VR-VBIAS/13.60580d0 
   !	 Write(*,*) VL,VR
      WRITE(*,*) 'BITTE WAEHLEN COMPUTATION MODE'
      WRITE(*,*) 'SPINLESS=(1), SPIN=(2), EP=(3)'
      READ(*,*) CASES
    
      SELECT CASE(CASES)
!_________________________________________________________________
	CASE(1)
   !  SPINLESS CASE
      WRITE(*,*) 'NOW SPINLESS CALCULATION' 
	open (600,file='abscurrent.s1',status='old')
	do i=1,L
	read (600,*)  x1(i), x2(i),x3(i),x4(i),x5(i),x6(i),x7(i)
      end do  
       
!_____________________________________________________________	 
	CASE(2)    
   !  SPIN CASE
      WRITE(*,*) 'NOW SPIN CALCULATION'
	open (6000,file='abscurrent_UP_benzene.ss1',status='old')
	do i=1,L
	read (6000,*)  x1(i), x2(i),x3(i),x4up(i), x5(i),x6(i),x7(i)
      end do  


      open (6001,file='abscurrent_DN_benzene.ss1',status='old')
	do i=1,L
	read (6001,*)  x1(i), x2(i),x3(i),x4dn(i), x5(i),x6(i),x7(i)
      end do  


!________________________________________________________________
      CASE(3)
!     WITH ELECTRON-PHONON, SPINLESS 
      WRITE(*,*) 'WITH ELECTRON-PHONON, SPINLESS'

      open (600,file='abscurrent.s1',status='old')
	do i=1,L
	read (600,*)  x1(i), x2(i),x3(i),x4(i),x5(i),x6(i),x7(i)
      end do  

      open (601,file='couple.dat', status='old')
	open (602,file='w.input',status='old')
	do i=1,NMODE
      read (601,*) CLR(i), CRR(i)
	read (602,*) fr(i)
	fr(i)=fr(i)/(8.0665d3)/0.9612d0/13.6058 
	END DO
!____________________________________________________________________
      CASE DEFAULT
	WRITE(*,*) 'TRY AGAIN' 
      STOP
      END SELECT


      do i=1,L
      open (UNIT=7,file='1.dat')
      write(7,1300)  x1(i)*13.60580d0-VL*13.60580d0,x4(i)  
 1300 format (f16.10,2x,f16.10,2x,f16.10,2x,f16.10,2x,f16.10) 
	end do 
 

      SELECT CASE(CASES) 
 !_____________________________________________________________________
       
      CASE(1)
      TR=1.0D0
      DT=0.0d0 
      do while (TR.LT.3.0D0) 
      TL=TR 

      CALL CURRENT(L,TL,TR,DT,x1,x4,CC,DD,EG,HG,HH,SEE,KEQ)
      
      open (UNIT=10,file='SEE.dat')   
      write(10,2400) TR, SEE 
 2400	Format(f20.10,2x,f20.15)

	TR=TR+1.0    
	end do 
      
      write(*,*) 'Task done!' 
	STOP
 
!__________________________________________________________________
      
      CASE(2)
	TR=0.0D0

      WRITE(*,5401) 'TR','SEEBECK','KEL'
5401  FORMAT(A4,20x,A9,20x,A4)

      do while (TR.LE.300) 

      TL=TR 

      CALL SCURRENT(L,TL,TR,DT,x1,x4up,x4dn,CC,DD,SEE,KQ)

      open (UNIT=10,file='SEE.dat')
      WRITE(10,2401) TR, SEE

	WRITE(*,2402) TR, SEE, KQ
	 
 2401	Format(f20.10,2x,f20.15)
 2402 Format(f20.10,2x,f20.15,2x,e8.2)
 
	TR=TR+1.0D0    
	end do 
      
      write(*,*) 'Task done!' 
	STOP
!_________________________________________________________________	 
      CASE(3)

      TL=4.0d0
	TR=TL
       
	VBIAS=0.d0
      DO WHILE (VBIAS.LT.20.d0)
	VR=-0.2536999997d0 
	VL=VR-(VBIAS/13.60580d3) 
	
      CALL CURRENT_EP(L,TL,TR,DT,fr,x1,x4,CC,DD,EG,HG,HH,SEE,KEQ,CLR,CRR
     *,Tw,SEE1)
      CALL thermcd(NCH2,TL,TR,scaleM,KQQ)
      ZT=SEE**2*1.D-12*HH*(TL+TR)/2.0d0/(KEQ) !,KQ

      OPEN (UNIT=9,file='ZT.dat')   
      WRITE(9,1402) TR,  ZT
 1402 Format(f20.10,2x,f20.15)
	WRITE(*,*) VBIAS, SEE, SEE1
	VBIAS=VBIAS+1.0d0 
	END DO

	WRITE(*, *) 'TASK DONE'
	STOP
!____________________________________________________________________-
      CASE DEFAULT
	STOP

      END SELECT  
!_____________________________________________________________________________


  
	END PROGRAM 
