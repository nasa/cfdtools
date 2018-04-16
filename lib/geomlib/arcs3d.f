C+------------------------------------------------------------------------------
C
      SUBROUTINE ARCS3D (IMIN, IMAX, JMIN, JMAX, KMIN, KMAX,
     >                   I1, I2, J1, J2, K1, K2, X, U, V, W)
C
C     One-liner: Normalized arc-lengths for a 3-space subgrid stored as triples
C
C     ARCS3D calculates normalized arc-lengths for the indicated portion of a
C     regular 3-space grid stored as (x,y,z) triples.  If a single plane is
C     specified (e.g., I1 = I2), the degenerate direction is not assigned.
C     Collapsed lines are handled with uniform distributions.
C
C     Note that U, V, & W are separated for compatibility with RIPPLE2D & -3D.
C
C     12/94 - 08/96  DAS  Original PARAMXYZ assumed x(i,j,k), y(i,j,k), etc.
C     10/12/99        "   ARCS3D adapted from PARAMXYZ for X(1:3,i,j,k) storage.
C     11/02/99        "   Belatedly realized S(i,j,k,1:3) is not tenable.
C
C  AUTHOR:  David Saunders, Raytheon/NASA Ames, Moffett Field, CA.
C
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   IMIN, IMAX, JMIN, JMAX, KMIN, KMAX,      ! Grid array dimensions
     >   I1, I2, J1, J2, K1, K2                   ! Active volume; one pair
                                                  ! may be equal
      REAL, INTENT (IN) ::
     >   X (1:3, IMIN:IMAX, JMIN:JMAX, KMIN:KMAX) ! Grid coordinates as triples

      REAL, INTENT (OUT), DIMENSION (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX) ::
     >   U, V, W                                  ! Normalized arc-lengths
                                                  ! for the I,J,K lines:
                                                  !   U(I1,J,K) = 0.,
                                                  !   V(I,J1,K) = 0.,
                                                  !   W(I,J,K1) = 0.,
                                                  !   U(I2,J,K) = 1., etc.
C-------------------------------------------------------------------------------

C     Local constants:

      REAL, PARAMETER ::
     >   EPS = 1.E-6, ONE = 1.E+0, ZERO = 0.E+0

C     Local variables:

      INTEGER
     >   I, J, K

      REAL
     >   RTOTAL

C     Local functions:

      REAL
     >   DELI, DELJ, DELK

      DELI (I, J, K) = SQRT ((X (1, I, J, K) - X (1, I - 1, J, K)) ** 2
     >                     + (X (2, I, J, K) - X (2, I - 1, J, K)) ** 2
     >                     + (X (3, I, J, K) - X (3, I - 1, J, K)) ** 2)

      DELJ (I, J, K) = SQRT ((X (1, I, J, K) - X (1, I, J - 1, K)) ** 2
     >                     + (X (2, I, J, K) - X (2, I, J - 1, K)) ** 2
     >                     + (X (3, I, J, K) - X (3, I, J - 1, K)) ** 2)

      DELK (I, J, K) = SQRT ((X (1, I, J, K) - X (1, I, J, K - 1)) ** 2
     >                     + (X (2, I, J, K) - X (2, I, J, K - 1)) ** 2
     >                     + (X (3, I, J, K) - X (3, I, J, K - 1)) ** 2)

C     Execution:

C     Zero the three low-end faces (or edges if one plane is specified):

      DO K = K1, K2
         U (I1, J1:J2, K)  = ZERO
         V (I1:I2, J1, K)  = ZERO
      END DO

      W (I1:I2, J1:J2, K1) = ZERO

C     Set up the low-end edge lines because they are missed by the
C     following loops over most of the low-end faces:

      DO I = I1 + 1, I2
         U (I, J1, K1) =  U (I - 1, J1, K1) + DELI (I, J1, K1)
      END DO

      DO J = J1 + 1, J2
         V (I1, J, K1) =  V (I1, J - 1, K1) + DELJ (I1, J, K1)
      END DO

      DO K = K1 + 1, K2
         W (I1, J1, K) =  W (I1, J1, K - 1) + DELK (I1, J1, K)
      END DO

C     Set up the rest of the low-end face lines because they are
C     missed by the the main loop over most of the volume:

      DO K = K1 + 1, K2
         DO J = J1 + 1, J2
            V (I1, J, K) =  V (I1, J - 1, K) + DELJ (I1, J, K)
            W (I1, J, K) =  W (I1, J, K - 1) + DELK (I1, J, K)
         END DO

         DO I = I1 + 1, I2
            U (I, J1, K) =  U (I - 1, J1, K) + DELI (I, J1, K)
            W (I, J1, K) =  W (I, J1, K - 1) + DELK (I, J1, K)
         END DO
      END DO

      DO J = J1 + 1, J2
         DO I = I1 + 1, I2
            U (I, J, K1) =  U (I - 1, J, K1) + DELI (I, J, K1)
            V (I, J, K1) =  V (I, J - 1, K1) + DELJ (I, J, K1)
         END DO
      END DO

C     Traverse the block just once for all lines except those within
C     the low-end faces:

      DO K = K1 + 1, K2
         DO J = J1 + 1, J2
            DO I = I1 + 1, I2
               U (I, J, K) =  U (I - 1, J, K) + DELI (I, J, K)
               V (I, J, K) =  V (I, J - 1, K) + DELJ (I, J, K)
               W (I, J, K) =  W (I, J, K - 1) + DELK (I, J, K)
            END DO
         END DO
      END DO

C     Normalizing requires another pass through the volume.
C     Degenerate cases must be avoided, so separate the 3 directions:

      IF (I1 /= I2) THEN

         DO K = K1, K2
            DO J = J1, J2
               RTOTAL = U (I2, J, K)

               IF (RTOTAL < EPS) THEN
                  DO I = I1 + 1, I2
                     U (I, J, K) = REAL (I - I1)
                  END DO
                  RTOTAL = U (I2, J, K)
               END IF

               RTOTAL = ONE / RTOTAL
               DO I = I1 + 1, I2 - 1
                  U (I, J, K) = U (I, J, K) * RTOTAL
               END DO
               U (I2, J, K) = ONE
            END DO
         END DO

      END IF

      IF (J1 /= J2) THEN

         DO K = K1, K2
            DO I = I1, I2
               RTOTAL = V (I, J2, K)

               IF (RTOTAL < EPS) THEN
                  DO J = J1 + 1, J2
                     V (I, J, K) = REAL (J - J1)
                  END DO
                  RTOTAL = V (I, J2, K)
               END IF

               RTOTAL = ONE / RTOTAL
               DO J = J1 + 1, J2 - 1
                  V (I, J, K) = V (I, J, K) * RTOTAL
               END DO
               V (I, J2, K) = ONE
            END DO
         END DO

      END IF

      IF (K1 /= K2) THEN

         DO J = J1, J2
            DO I = I1, I2
               RTOTAL = W (I, J, K2)

               IF (RTOTAL < EPS) THEN
                  DO K = K1 + 1, K2
                     W (I, J, K) = REAL (K - K1)
                  END DO
                  RTOTAL = W (I, J, K2)
               END IF

               RTOTAL = ONE / RTOTAL
               DO K = K1 + 1, K2 - 1
                  W (I, J, K) = W (I, J, K) * RTOTAL
               END DO
               W (I, J, K2) = ONE
            END DO
         END DO

      END IF

      END SUBROUTINE ARCS3D
