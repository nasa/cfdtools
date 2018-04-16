C+------------------------------------------------------------------------------
C
      SUBROUTINE PARAM2D (IDIM, JDIM, I1, I2, J1, J2, X, Y, Z, U, V)
C
C ONE-LINER: PARAMeterization in 2-space of an XYZ surface
C            -----               -
C PURPOSE:
C
C        PARAM2D parameterizes the given (sub)mesh on a surface using
C     the usual chord-length approximation to arc-length between grid
C     points.  The arc length of each row and column is generally
C     normalized to 1, but this version can suppress normalization.
C     U(I1,J1) should be input as 0. or -999. respectively.
C
C        This version also handles degenerate lines by inserting uniform
C     u or v for the normalized case.
C
C        For the unnormalized case, this version checks each edge and
C     replaces arc lengths of a collapsed edge with a scaled form of
C     the adjacent arc lengths (assuming no interior line is collapsed).
C
C        It is hoped that tying the dimensions of U and V to those of
C     X and Y (to keep the argument list modest) is not unwise, though
C     it could be if small submeshes are commonly involved.
C
C HISTORY:
C
C     11/19/93  DAS  Initial implementation.
C     12/03/93   "   Normalized each row and column.
C     03/20/95   "   Suppress normalization kludge: U(I1,J1) = -999.
C     02/06/98   "   Handled degenerate lines by inserting uniform u or v.
C     09/19/07   "   The unnormalized case also handles collapsed edge now.
C
C AUTHOR: David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER :: IDIM, JDIM     !(I) Grid dimensions in the calling program.

      INTEGER :: I1, I2, J1, J2 !(I) Define the submesh to be parameterized.

      REAL :: X (IDIM, JDIM),   !(I) The surface grid coordinates.
     >        Y (IDIM, JDIM),
     >        Z (IDIM, JDIM)

      REAL :: U (IDIM, JDIM),   !(O) The chord-length-based parameterization:
     >        V (IDIM, JDIM)    !    U (I, J) = SUM (K=2:I) dS (K) / S (J) where
                                !    dS (K) ** 2 = (X (K, J) - X (K-1, J) ** 2 +
                                !                  (Y (K, J) - Y (K-1, J) ** 2 +
                                !                  (Z (K, J) - Z (K-1, J) ** 2
                                !    and S (J) is the total length of row J;
                                !    similarly for V (I, J).
                                !(I) KLUDGE: If U(I1,J1) = -999. on input,
                                !    the normalization is suppressed.
C-------------------------------------------------------------------------------

C     Local constants.

      REAL, PARAMETER :: EPS = 1.E-6, FLAG = -999., ONE = 1., ZERO = 0.

C     Local variables.

      INTEGER :: I, INC, J
      REAL    :: RLENGTH
      LOGICAL :: NORM

C     Execution.

      NORM = U (I1, J1) /= FLAG

      DO J = J1, J2

         U (I1, J) = ZERO

         DO I = I1 + 1, I2
            U (I, J) = U (I - 1, J)  +  SQRT (
     >         (X (I, J) - X (I - 1, J)) ** 2 +
     >         (Y (I, J) - Y (I - 1, J)) ** 2 +
     >         (Z (I, J) - Z (I - 1, J)) ** 2)
         END DO

         IF (NORM) THEN

            IF (U (I2, J) > EPS) THEN
               RLENGTH = ONE / U (I2, J)
               DO I = I1 + 1, I2 - 1
                  U (I, J) = U (I, J) * RLENGTH
               END DO
            ELSE
               RLENGTH = ONE / REAL (I2 - I1)
               DO I = I1 + 1, I2 - 1
                  U (I, J) = REAL (I - I1) * RLENGTH
               END DO
            END IF
            U (I2, J) = ONE

         END IF

      END DO

      IF (.NOT. NORM) THEN  ! Check for a collapsed edge

         INC = 1
         DO J = J1, J2, J2 - J1
            IF (U (I2, J) < EPS) THEN  ! Use scaled-down adjacent values
               DO I = I1 + 1, I2
                  U (I, J) = EPS * U (I, J + INC)
               END DO
            END IF
            INC = -INC
         END DO

      END IF

      DO I = I1, I2

         V (I, J1) = ZERO

         DO J = J1 + 1, J2
            V (I, J) = V (I, J - 1)  +  SQRT (
     >         (X (I, J) - X (I, J - 1)) ** 2 +
     >         (Y (I, J) - Y (I, J - 1)) ** 2 +
     >         (Z (I, J) - Z (I, J - 1)) ** 2)
         END DO

         IF (NORM) THEN

            IF (V (I, J2) > EPS) THEN
               RLENGTH = ONE / V (I, J2)
               DO J = J1 + 1, J2 - 1
                  V (I, J) = V (I, J) * RLENGTH
               END DO
            ELSE
               RLENGTH = ONE / REAL (J2 - J1)
               DO J = J1 + 1, J2 - 1
                  V (I, J) = REAL (J - J1) * RLENGTH
               END DO
            END IF
            V (I, J2) = ONE

         END IF

      END DO

      IF (.NOT. NORM) THEN  ! Check for a collapsed edge

         INC = 1
         DO I = I1, I2, I2 - I1
            IF (V (I, J2) < EPS) THEN
               DO J = J1 + 1, J2
                  V (I, J) = EPS * V (I + INC, J)
               END DO
            END IF
            INC = -INC
         END DO

      END IF

      END SUBROUTINE PARAM2D
