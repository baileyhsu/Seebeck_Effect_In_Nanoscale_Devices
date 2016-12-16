   SUBROUTINE ZEROIN (FCT, ABSERR, RELERR, MAXFCT, ITAPE, A, B, FB,  &
      NUMFCT, IERR)
!
!*****************************************************************
!                                                                *
!  The SUBROUTINE ZEROIN finds a zero of odd order of a          *
!  continuous function FCT in the interval [A,B] provided that   *
!  FCT(A) and FCT(B) have different signs, i.e.,                 *
!  FCT(A)*FCT(B) < 0.                                            *
!  The Zeroin method combines the bisection and the secant method*
!  with inverse quadratic interpolation.                         *
!                                                                *
!                                                                *
!  INPUT PARAMETERS:                                             *
!  =================                                             *
!  FCT    : Function, whose zero one wants to find. It has the   *
!           form   DOUBLE PRECISION FUNCTION FCT(X)  and must be *
!           declared as EXTERNAL in the calling program, or as   *
!           INTRINSIC, if it is described by standard FORTRAN-77 *
!           functions.                                           *
!  ABSERR : ) error bounds with ABSERR >= 0 and RELERR >= 0.     *
!  RELERR : ) Their sum must be positive. The break-off criterion*
!             used is as follows:                                *
!             ABS(XM) <= 0.5*(ABSERR+ABS(B)*RELERR), where       *
!             XM denotes half the interval length, XM = (B-A)/2. *
!             If RELERR=0, then we only test the absolute error; *
!             if ABSERR=0, we only test the relative arror.      *
!             The input values for ABSERR and RELERR are only    *
!             used by the subroutine if each exceeds four times  *
!             the machine constant, or if one is zero then the   *
!             other must exceed that bound. Otherwise one or both*
!             of them are internally adjusted to four times the  *
!             machine constant.                                  *
!  MAXFCT : Maximal number of function evaluations allowed       *
!  ITAPE  : > 0, Number for a data set that will absorb inter-   *
!                mediate results                                 *
!  A, B   : endpoints of the interval that contains a zero of FCT*
!                                                                *
!                                                                *
!  OUTPUT PARAMETERS:                                            *
!  ==================                                            *
!  ABSERR : ) actually used error bounds                         *
!  RELERR : )                                                    *
!  B      : approximate zero                                     *
!  FB     : functionmal value at the approximate zero B          *
!  NUMFCT : number of actual functional evaluations performed    *
!  IERR   : error parameter:                                     *
!           =-2, ABSERR or RELERR is negative or both are zero   *
!                or  MAXFCT < 1                                  *
!           =-1, FCT(A)*FCT(B) < 0.0 is not true                 *
!           = 0, A or B are numerical zeros of FCT               *
!           = 1, B is a zero with FCT(B)=0.0                     *
!           = 2, the desired accuracy has been achieved:         *
!                ABS(XM) <= error bound                          *
!           = 3, maximal number of function evaluations has been *
!                reached without meeting the break-off criterion *
!                                                                *
!----------------------------------------------------------------*
!                                                                *
!  Required subroutines: MACHPD, BI                              *
!                                                                *
!*****************************************************************
!                                                                *
!  Authors     : Siegmar Neuner, Gisela Engeln-M\201llges           *
!  Date        : 06.01.1992                                      *
!  Source      : FORTRAN 77                                      *
!                                                                *
!*****************************************************************
!                
      DOUBLEPRECISION FMACHP, FA, FB, FC, A, B, C, D, E, ABSERR, RELERR,&
      TOL1, EPS, XM, R, Q, S, P, FCT, HELP
      INTEGER MAXFCT, ITAPE, NUMFCT, IERR, MACHPD

!                
!  Compute four times the machine constant FMACHP
!                
      FMACHP = 1.0D0
   10 FMACHP = 0.5D0 * FMACHP
      IF (MACHPD (1.0D0 + FMACHP) .EQ.1) GOTO 10
      FMACHP = 2.0D0 * FMACHP
      EPS = 4.0D0 * FMACHP
!                
!  Compute FCT at A, B
!                
      FA = FCT (A)
      FB = FCT (B)
      NUMFCT = 2
!                
!  Check whether  FCT(A)*FCT(B) < 0.0D0
!                
      HELP = FA * FB
      IF (HELP.GT.0.0D0) THEN
         IERR = - 1
         RETURN
      ELSEIF (HELP.EQ.0.0D0) THEN
         IERR = 0
         RETURN
      ENDIF
!                
!  Check input eror parameters
!                
      IF (ABSERR.LT.0.0D0.OR.RELERR.LT.0.0D0.OR.ABSERR +                &
      RELERR.LE.0.0D0.OR.MAXFCT.LT.1) THEN
         IERR = - 2
         RETURN
      ENDIF
      IF (RELERR.EQ.0.0D0) THEN
         IF (ABSERR.LT.EPS) ABSERR = EPS
      ELSEIF (ABSERR.EQ.0.0D0) THEN
         IF (RELERR.LT.EPS) RELERR = EPS
      ELSE
         IF (ABSERR.LT.EPS) ABSERR = EPS
         IF (RELERR.LT.EPS) RELERR = EPS
      ENDIF
!                
!  No zero between  B and C, set C so that there is a
!  zero between  B and C
!                
   20 C = A
      FC = FA
      D = B - A
      E = D
!                
!  If  FC is the smaller sized function value,
!  swap the interval ends
!                
   30 IF (DABS (FC) .LT.DABS (FB) ) THEN
         A = B
         B = C
         C = A
         FA = FB
         FB = FC
         FC = FA
      ENDIF
!                
!  TOL1 is an auxiliary variable used for the mixed error test
!                
      TOL1 = 0.5D0 * (ABSERR + RELERR * DABS (B) )
!                
!  take half the interval length XM
!                
      XM = 0.5D0 * (C - B)
      IF (ITAPE.GT.0) THEN
         WRITE (ITAPE, 900) A, B, C
         WRITE (ITAPE, 910) FA, FB, FC
      ENDIF
      R = 0.0D0
!                
!  If XM is less than TOL1, we have achieved a sufficiently
!  good approximate zero
!                
      IF (DABS (XM) .LE.TOL1) THEN
         IERR = 2
         RETURN
!                
!  Check whether the value  FB  of the best approximation to
!  a zero already vanishes
!                
      ELSEIF (FB.EQ.0.0D0) THEN
         IERR = 1
         RETURN
      ENDIF
      IF (DABS (E) .LT.TOL1) THEN
         CALL BI (XM, ITAPE, D, E)
      ELSE
         IF (DABS (FA) .LE.DABS (FB) ) THEN
            CALL BI (XM, ITAPE, D, E)
         ELSE
!                
!  If  A and C are different, then together with B one can use these thr
!  points for an inverse quadratse interpolation
!                
            IF (A.NE.C) THEN
               Q = FA / FC
               R = FB / FC
               S = FB / FA
               P = S * (2.0D0 * XM * Q * (Q - R) - (B - A) * (R - 1.0D0)&
               ) 
               Q = (Q - 1.0D0) * (R - 1.0D0) * (S - 1.0D0)
            ELSE
!                
!  Here we use the secant method or interpolate linearly
!                
               S = FB / FA
               P = 2.0D0 * XM * S
               Q = 1.0D0 - S
            ENDIF
!                
!  The sign of P/Q is reversed for the following division
!                
            IF (P.GT.0.0D0) THEN
               Q = - Q
            ELSE
               P = DABS (P)
            ENDIF
            IF ( (2.0D0 * P) .GE. (3.0D0 * XM * Q - DABS (TOL1 * Q) ) ) &
            THEN 
               CALL BI (XM, ITAPE, D, E)
            ELSE
               IF (P.GE.DABS (0.5D0 * E * Q) ) THEN
                  CALL BI (XM, ITAPE, D, E)
               ELSE
!                
!  For either interpolation we compute the quotient P/Q
!  which shall be added to B
!                
                  E = D
                  D = P / Q
                  IF (ITAPE.GT.0) THEN
                     IF (R.EQ.0.0D0) THEN
                        WRITE (ITAPE, * ) 'Secant method'
                     ELSE
                        WRITE (ITAPE, * ) 'Inverse quadratic ',         &
                        'interpolation'
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF
!                
!  The previous best approximate zero B is stored in A
!  and the function value FB in FA
!                
      A = B
      FA = FB
!                
!  If D exceeds TOL1, it is added to B
!                
      IF (DABS (D) .GT.TOL1) THEN
         IF (ITAPE.GT.0) WRITE (ITAPE, 920) D
         B = B + D
      ELSE
!                
!  The desired accuracy has been achieved.
!  The best approximate zero B is improved by adding the error bound.
!                
         IF (ITAPE.GT.0) THEN
            IF (XM.LT.0.0D0) THEN
               WRITE (ITAPE, 930) TOL1
            ELSE
               WRITE (ITAPE, 940) TOL1
            ENDIF
         ENDIF
         B = B + DSIGN (TOL1, XM)
      ENDIF
!                
!  Compute the new value FB at B
!                
      FB = FCT (B)
!                
!  The iterationcounter is upped by 1
!                
      NUMFCT = NUMFCT + 1
      IF (ITAPE.GT.0) WRITE (ITAPE, 950) B, FB, NUMFCT
      IF (NUMFCT.GT.MAXFCT) THEN
         IERR = 3
         RETURN
      ENDIF
!                
!  If the signs of the function at  B and C are opposite, then
!  there is a zero between B and C
!                
      IF ( (FB * (FC / DABS (FC) ) ) .GT.0.0D0) GOTO 20
      GOTO 30
!                
  900 FORMAT(1X,'A = ',D20.14,'  B = ',D20.14,'  C = ',D20.14)
  910 FORMAT(1X,'FA= ',D20.14,'  FB= ',D20.14,'  FC= ',D20.14)
  920 FORMAT(1X,'distance to the new B:  D= ',D20.14)
  930 FORMAT(1X,'error bound is subtracted:  D= -',D20.14)
  940 FORMAT(1X,'error bound is added:  D= ',D20.14)
  950 FORMAT(1X,'B = ',D20.14,'  FB= ',D20.14,/,                        &
     &       1X,'number of function evaluations: ',I4,//)
      
END SUBROUTINE ZEROIN