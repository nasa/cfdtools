C+----------------------------------------------------------------------
C
      SUBROUTINE VINOKUR (I1, I2, D1, D2, X, LUNOUT, IER)
C
C     VINOKUR imposes the indicated 2-sided Vinokur distribution on the
C     interior points of the grid line between I1 and I2. The end points
C     are assumed to be in place, and may be in either order: this
C     routine simplifies use of the underlying HTDIS4 for the case where
C     X(I2) is less than X(I1). D1 and D2 should be positive, full scale.
C     See HTDIS4 for further details, including LUNOUT & IER arguments.
C
C     07/07/97  DAS  2-sided variant of IC1 from SYN87.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER  I1, I2, LUNOUT, IER
      REAL     D1, D2, X(1:I2)

C     Local constants.

      REAL     ONE, ZERO
      PARAMETER (ONE = 1., ZERO = 0.)

C     Local variables.

      INTEGER  I
      REAL     RANGE, XI1, XI2

C     Execution.

      XI1 = X(I1)
      XI2 = X(I2)
      RANGE = XI2 - XI1

C     Normalized Vinokur distribution in ascending order:

      CALL HTDIS4 (.TRUE., ZERO, ONE, ABS (D1/RANGE), ABS (D2/RANGE),
     >             I2 - I1 + 1, X(I1), -LUNOUT, IER)

C     Denormalize:

      DO I = I1, I2
         X(I) = XI1 + X(I) * RANGE
      END DO

      X(I2) = XI2  ! Exactly

      RETURN
      END
