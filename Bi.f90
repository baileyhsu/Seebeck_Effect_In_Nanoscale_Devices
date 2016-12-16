     SUBROUTINE BI (XM, ITAPE, D, E)
!                
!*****************************************************************
!                                                                *
!  The SUBROUTINE BI is a bisection algorithm for the SUBROUTINE *
!  ZEROIN.                                                       *
!                                                                *
!*****************************************************************
!                
      DOUBLEPRECISION XM, D, E
      D = XM
      E = D
      IF (ITAPE.GT.0) WRITE (ITAPE, * ) 'Bisection'
      RETURN
      END SUBROUTINE BI  