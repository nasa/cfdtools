C+------------------------------------------------------------------------------
C
      SUBROUTINE PARAMXY (IMIN, IMAX, JMIN, JMAX, I1, I2, J1, J2,
     >                    X, Y, S)
C
C  ONE-LINER: Relative arc-lengths for all lines of a 2-space grid
C
C  DESCRIPTION:
C
C        PARAMXY parameterizes a regular 2-space grid by setting up
C     the normalized arc-length increments in both index directions.
C
C  ENVIRONMENT:  FORTRAN 77 with minor extensions
C
C     12/13/94  DAS  Initial code.
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER IMIN, IMAX, JMIN, JMAX      !  I  Grid array dimensions.
      INTEGER I1, I2, J1, J2              !  I  Active area is (I1:I2, J1:J2).
      REAL    X (IMIN:IMAX, JMIN:JMAX),   !  I  Grid coordinates.
     >        Y (IMIN:IMAX, JMIN:JMAX)
      REAL    S (IMIN:IMAX, JMIN:JMAX, 2) !  O  Relative arc-length increments
                                          !     for lines in the I, J directions
                                          !     with S(I1,J,1) = S(I,J1,2) = 0.
                                          !     and  S(I2,J,1) = S(I,J2,2) = 1.
C-------------------------------------------------------------------------------

C     Local constants.

      REAL    ONE, ZERO
      PARAMETER (ONE = 1.E+0, ZERO = 0.E+0)

C     Local variables.

      INTEGER I, J

C     Local functions.

      REAL    DELI, DELJ

      DELI (I, J) = SQRT ((X (I, J) - X (I - 1, J)) ** 2 +
     >                    (Y (I, J) - Y (I - 1, J)) ** 2)

      DELJ (I, J) = SQRT ((X (I, J) - X (I, J - 1)) ** 2 +
     >                    (Y (I, J) - Y (I, J - 1)) ** 2)

C     Execution.
C     ----------

C     Zero the two low-end edges.

      DO J = J1, J2
         S (I1, J, 1) = ZERO
      END DO

      DO I = I1, I2
         S (I, J1, 2) = ZERO
      END DO

C     Set up the low-end edge lines because they are missed by the
C     following loops over most of the interior:

      DO I = I1 + 1, I2
         S (I, J1, 1) =  S (I - 1, J1, 1) + DELI (I, J1)
      END DO

      DO J = J1 + 1, J2
         S (I1, J, 2) =  S (I1, J - 1, 2) + DELJ (I1, J)
      END DO

C     Traverse the grid just once for all lines except those within
C     the low-end edges.

      DO J = J1 + 1, J2
         DO I = I1 + 1, I2
            S (I, J, 1) =  S (I - 1, J, 1) + DELI (I, J)
            S (I, J, 2) =  S (I, J - 1, 2) + DELJ (I, J)
         END DO
      END DO

C     Normalizing requires another pass through the grid:

      DO J = J1, J2
         DO I = I1, I2
            S (I, J, 1) = S (I, J, 1) / S (I2, J, 1)
            S (I, J, 2) = S (I, J, 2) / S (I, J2, 2)
         END DO
      END DO

C     Finally, precise 1s for the two high-end edges:

      DO J = J1, J2
         S (I2, J, 1) = ONE
      END DO

      DO I = I1, I2
         S (I, J2, 2) = ONE
      END DO

      RETURN
      END
