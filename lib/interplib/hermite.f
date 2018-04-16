C+----------------------------------------------------------------------
C
      SUBROUTINE HERMITE (Y1, Y2, D1, D2, NEVAL, XEVAL, YEVAL)
C
C     Hermite evaluates the Hermite cubic on the unit interval defined by
C     two points (0,Y1), (1,Y2), and the corresponding derivatives.
C
C     History:
C
C     08/22/97  D.A.Saunders  Implemented to illustrate problems with
C                             Hermite-based TFI.
C
C-----------------------------------------------------------------------
    
      IMPLICIT NONE

C     Arguments:

      REAL       Y1, Y2, D1, D2  ! Y and dY/dX at X = 0 and X = 1
      INTEGER    NEVAL           ! # evaluations
      REAL       XEVAL(NEVAL)    ! Target abscissas in [0, 1]
      REAL       YEVAL(NEVAL)    ! Output evaluations

C     Local constants:

      REAL       ONE, TWO, THREE
      PARAMETER (ONE = 1., TWO = 2., THREE = 3.)

C     Local variables:

      INTEGER    I
      REAL       CY1, CY2, CD1, CD2, XNORM

C     Execution:

      DO I = 1, NEVAL
         XNORM = XEVAL(I)
         CY1 = (ONE + TWO*XNORM) * (ONE - XNORM)**2
         CY2 = (THREE - TWO*XNORM) * XNORM**2
         CD1 = (ONE - XNORM)**2 * XNORM
         CD2 = (XNORM - ONE) * XNORM**2
         YEVAL(I) = CY1 * Y1  + CY2 * Y2 + CD1 * D1 + CD2 * D2
      END DO

      END
