!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE SURFPH(GG1,GG2,GG3)
!For calculating the surface phonon spectral function by Geller's paper
!PRB 64,155320(2001)
!11/18/2002
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      REAL*8 SC3(2048),WGTES(2048),TC(2048),WGTET(2048)
      COMMON/SURF/VL,VT,CR
!INPUT:VL,VT,RHO
        N=2048
!parameter: HBAR
        PI=ACOS(-1.D0)
!
!    RAYLEIGH BRANCH
!TEST
!       WRITE(*,*)'Xsi=',B
!       WRITE(*,*)'FB=',FB
!       WRITE(*,*)'NUMFCT=',NUMFCT
!       WRITE(*,*)'ERR=',IERR
!ENDTEST
           GAMA=DSQRT(1.D0-(CR/VL)**2)
           ETA=DSQRT(1.D0-(CR/VT)**2)
           ALFARSQ=2.D0*(GAMA**3)*(ETA**2)*(1.D0-(2.D0/(1.D0+ETA**2)))**2   &
           /((GAMA-ETA)*(GAMA-ETA+2.D0*GAMA*ETA**2))
!TEST
        GG1=ALFARSQ
!ENDTEST
           SPRAY=ALFARSQ/CR**3
! PLUS-MINUS BRANCH
                X1=0.D0
                X2=1.D0/VL
        CALL GaussLeg(X1,X2,N,SC3,WGTES)
                SPPLUS=0.D0
                SPMINUS=0.D0
        DO J=1,N
                F1Z=1.D0/SC3(J)
                CALL FZPNSQ(F1Z,FPSQ,FNSQ)
                SPPLUS=SPPLUS+(FPSQ/(SC3(J)**2))*WGTES(J)
                SPMINUS=SPMINUS+(FNSQ/(SC3(J)**2))*WGTES(J)
        END DO
!TEST
        GG2=(SPPLUS+SPMINUS)*VL**3*PI
!ENDTEST
! ZERO BRANCH
                 X3=VT
                 X4=VL
         CALL GaussLeg(X3,X4,N,TC,WGTET)
                 SP0=0.D0
        DO J=1,N
                FZ=TC(J)
                CALL FZ0SQ(FZ,BBB)
                SP0=SP0+WGTET(J)*BBB
        END DO
!
        GG3=SP0*2.d0*VT**3*PI
         END