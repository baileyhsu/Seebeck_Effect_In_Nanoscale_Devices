!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE FZPNSQ(C,FPSQ,FNSQ)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       IMPLICIT INTEGER (I-N)
       COMPLEX*16 I,TEMP_P,TEMP_N
       COMMON/SURF/VL,VT,CR
        PI=ACOS(-1.D0)
           I=(0.D0,1.D0)
           ALFA=DSQRT((C/VL)**2-1.D0)
           BETA=DSQRT((C/VT)**2-1.D0)
           A=((BETA**2-1.D0)**2-4.D0*ALFA*BETA)&
             /((BETA**2-1.D0)**2+4.D0*ALFA*BETA)
           B=(4.D0*DSQRT(ALFA*BETA)*(BETA**2-1.D0))   &
           /((BETA**2-1.D0)**2+4.D0*ALFA*BETA)
           TEMP_P=DSQRT(ALFA)*(1.D0+A+I*B)+I*(1.D0-A-I*B)/DSQRT(BETA)
           TEMP_N=-DSQRT(ALFA)*(1.D0+A-I*B)+I*(1.D0-A+I*B)/DSQRT(BETA)
           FPSQ=DREAL(TEMP_P*DCONJG(TEMP_P))/(4.D0*PI*C)/C**3
           FNSQ=DREAL(TEMP_N*DCONJG(TEMP_N))/(4.D0*PI*C)/C**3
           RETURN
           END
