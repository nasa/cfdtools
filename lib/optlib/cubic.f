************************************************************************
*
      SUBROUTINE CUBIC (A, B, C, ROOT, IER)
*
*     Calculate the roots of the cubic equation x^3 + ax^2 + bx + c = 0
*     for real coefficients a, b, c.
*
*     The method is that from Numerical Recipes:
*
*     Letting  Q = (a^2 - 3b) / 9  and  R = (2a^3 - 9ab + 27c) / 54,
*     then if  R^2 < Q^3 there are three real roots.
*
*     Note that this excludes Q = 0, so degenerate cases such as
*     (X - 1)^3 = 0 are evidently not handled.
*
*     Other cases are not treated here yet.
*
*     02/13/01  David Saunders  Initial implementation (real roots only).
*               ELORET/NASA Ames
*
************************************************************************

      IMPLICIT NONE

*     Arguments:

      REAL,    INTENT (IN) ::
     >   A, B, C             ! Cubic coefficients

      REAL,    INTENT (OUT) ::
     >   ROOT(3)             ! Three real roots (undefined if IER /= 0)

      INTEGER, INTENT (OUT) ::
     >   IER                 ! 0 means three real roots were found;
                             ! 1 means roots are complex (or possibly
                             ! not distinct) and are not returned

*     Local constants:

      REAL, PARAMETER ::
     >   THIRD = 1./3.,
     >   TWOPI = 6.28318530717958647692528

*     Local variables:

      REAL
     >   Q, R, ROOTQ, THETA

*     Execution:

      Q = (A * A  - 3.* B) / 9.
      R = (A * (2.* A * A - 9.* B) + 27.* C) / 54.

      IF (Q > 0. .AND. R ** 2 < Q ** 3) THEN
         ROOTQ = SQRT (Q)
         THETA = ACOS (R / (Q * ROOTQ))
         ROOTQ = -2.* ROOTQ

         ROOT(1) = ROOTQ * COS ( THETA          * THIRD) - (A * THIRD)
         ROOT(2) = ROOTQ * COS ((THETA + TWOPI) * THIRD) - (A * THIRD)
         ROOT(3) = ROOTQ * COS ((THETA - TWOPI) * THIRD) - (A * THIRD)

         IER = 0
      ELSE
         IER = 1
      END IF

      END SUBROUTINE CUBIC
