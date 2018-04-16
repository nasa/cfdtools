C+------------------------------------------------------------------------------
C
      SUBROUTINE ARCS2D (IMIN, IMAX, JMIN, JMAX, I1, I2, J1, J2, X,
     >                   U, V)
C
C     One-liner: Normalized arc-lengths for a surface subgrid stored as triples
C
C     ARCS2D calculates normalized arc-lengths for the indicated portion of a
C     regular 3-space surface grid stored as (x,y,z) triples. (PARAM2D separates
C     x, y, and z; ARCS3D doesn't quite work either.)
C
C     Degenerate lines are handled with uniform distributions.
C     Note that U and V are compatible with RIPPLE2D searching.
C
C     11/02/99  DAS  Surface-only variant of ARCS3D; handled degeneracies.
C
C  AUTHOR:  David Saunders, Raytheon/NASA Ames, Moffett Field, CA.
C
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   IMIN, IMAX, JMIN, JMAX,       ! Grid array dimensions
     >   I1, I2, J1, J2                ! Active area

      REAL, INTENT (IN) ::
     >   X (1:3, IMIN:IMAX, JMIN:JMAX) ! Grid coordinates as triples

      REAL, INTENT (OUT) ::
     >   U (IMIN:IMAX, JMIN:JMAX),     ! Normalized arc-lengths
     >   V (IMIN:IMAX, JMIN:JMAX)      ! for the I & J lines:
                                       ! U(I1,J) = 0. = V(I,J1)
                                       ! U(I2,J) = 1. = V(I,J2)
C-------------------------------------------------------------------------------

C     Local constants:

      REAL, PARAMETER ::
     >   EPS = 1.E-6, ONE = 1.E+0, ZERO = 0.E+0

C     Local variables:

      INTEGER
     >   I, J

      REAL
     >   RTOTAL

C     Local functions:

      REAL
     >   DELI, DELJ

      DELI (I, J) = SQRT ((X (1, I, J) - X (1, I - 1, J)) ** 2 +
     >                    (X (2, I, J) - X (2, I - 1, J)) ** 2 +
     >                    (X (3, I, J) - X (3, I - 1, J)) ** 2)

      DELJ (I, J) = SQRT ((X (1, I, J) - X (1, I, J - 1)) ** 2 +
     >                    (X (2, I, J) - X (2, I, J - 1)) ** 2 +
     >                    (X (3, I, J) - X (3, I, J - 1)) ** 2)

C     Execution:

C     Zero the two low-end edges:

      U (I1, J1:J2) = ZERO
      V (I1:I2, J1) = ZERO

C     Set up the low-end edge lines because they are missed by the
C     following loops over most of the interior:

      DO I = I1 + 1, I2
         U (I, J1) =  U (I - 1, J1) + DELI (I, J1)
      END DO

      DO J = J1 + 1, J2
         V (I1, J) =  V (I1, J - 1) + DELJ (I1, J)
      END DO

C     Traverse the grid just once for the unnormalized arc lengths:

      DO J = J1 + 1, J2
         DO I = I1 + 1, I2
            U (I, J) =  U (I - 1, J) + DELI (I, J)
            V (I, J) =  V (I, J - 1) + DELJ (I, J)
         END DO
      END DO

C     Normalize in each direction, handling possible degeneracies:

      DO J = J1, J2
         RTOTAL = U (I2, J)

         IF (RTOTAL > EPS) THEN
            RTOTAL = ONE / RTOTAL
            DO I = I1 + 1, I2 - 1
               U (I, J) = U (I, J) * RTOTAL
            END DO
         ELSE
            RTOTAL = ONE / REAL (I2 - I1)
            DO I = I1 + 1, I2 - 1
               U (I, J) = REAL (I - I1) * RTOTAL
            END DO
         END IF

         U (I2, J) = ONE
      END DO

      DO I = I1, I2
         RTOTAL = V (I, J2)

         IF (RTOTAL > EPS) THEN
            RTOTAL = ONE / RTOTAL
            DO J = J1 + 1, J2 - 1
               V (I, J) = V (I, J) * RTOTAL
            END DO
         ELSE
            RTOTAL = ONE / REAL (J2 - J1)
            DO J = J1 + 1, J2 - 1
               V (I, J) = REAL (J - J1) * RTOTAL
            END DO
         END IF

         V (I, J2) = ONE
      END DO

      END SUBROUTINE ARCS2D
