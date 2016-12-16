!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE FZ0SQ(C,RESULT)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       IMPLICIT INTEGER (I-N)
       COMPLEX*16 I,TEMP_P,TEMP_N,D,E
       COMMON/SURF/VL,VT,CR
       PI=ACOS(-1.D0)
       I=(0.D0,1.D0)
       GAMA=DSQRT(1.D0-(CR/VL)**2)
       BETA=DSQRT((C/VT)**2-1.D0)
       D=((4.D0*BETA*(BETA**2-1.D0)**3)-16.d0*I*GAMA*BETA**2*(BETA**2-1.D0)) &
       /((BETA**2-1.D0)**4+16.D0*GAMA**2*BETA**2)
       E=((BETA**2-1.D0)**4-16*GAMA**2*BETA**2-8.D0*I*GAMA*BETA*(BETA**2-1.D0)**2) &
        /((BETA**2-1.D0)**4+16.D0*GAMA**2*BETA**2)
       TEMP_P=-GAMA*D+I*(1.D0-E)
       F0SQ=DREAL(TEMP_P*DCONJG(TEMP_P))/(2.D0*PI*C*BETA)/C**3
       RESULT=F0SQ
       RETURN
       END
