C+------------------------------------------------------------------------------
C
      SUBROUTINE PARAMXYZ (IMIN, IMAX, JMIN, JMAX, KMIN, KMAX,
     >                     I1, I2, J1, J2, K1, K2, X, Y, Z, S)
C
C  ONE-LINER: Relative arc-lengths for all lines of a 3-space subgrid
C
C  DESCRIPTION:
C
C        PARAMXYZ parameterizes the indicated portion of a regular
C     3-space grid by setting up the normalized arc-length increments
C     in all three index directions.  If a single plane is specified
C     (e.g. I1 = I2), the degenerate direction is assumed to be unused.
C     (All 1s are returned since no attempt is made to suppress them.)
C
C  ENVIRONMENT:  FORTRAN 77 with minor extensions
C
C     12/06/94  DAS/JJR  Separated set-up of original arc-lengths from
C                        subsequent grid perturbations; normalized them.
C                        Including all the face lines proved tedious!
C     12/16/94   "   "   Allowed for sub-volumes (or planes).
C     08/27/96    DAS    Degenerate directions are now handled properly.
C
C  AUTHOR:  David Saunders/James Reuther, NASA Ames, Mt. View, CA.
C
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER IMIN, IMAX, JMIN, JMAX, KMIN, KMAX     ! I Grid array dimensions.
      INTEGER I1, I2, J1, J2, K1, K2                 ! I Indicate active volume.
      REAL    X (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX),   ! I Grid coordinates.
     >        Y (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX),
     >        Z (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX)
      REAL    S (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX, 3) ! O Relative arc-lengths
                                                     !   for the I,J,K lines:
                                                     !   S(I1,J,K,1) = 0.,
                                                     !   S(I,J1,K,2) = 0.,
                                                     !   S(I,J,K1,3) = 0.,
                                                     !   S(I2,J,K,1) = 1., etc.
C-------------------------------------------------------------------------------

C     Local constants.

      REAL    EPS, ONE, ZERO
      PARAMETER (EPS = 1.E-30, ONE = 1.E+0, ZERO = 0.E+0)

C     Local variables.

      INTEGER I, J, K
      REAL    RTOTAL

C     Local functions.

      REAL    DELI, DELJ, DELK

      DELI (I, J, K) = SQRT ((X (I, J, K) - X (I - 1, J, K)) ** 2 +
     >                       (Y (I, J, K) - Y (I - 1, J, K)) ** 2 +
     >                       (Z (I, J, K) - Z (I - 1, J, K)) ** 2)

      DELJ (I, J, K) = SQRT ((X (I, J, K) - X (I, J - 1, K)) ** 2 +
     >                       (Y (I, J, K) - Y (I, J - 1, K)) ** 2 +
     >                       (Z (I, J, K) - Z (I, J - 1, K)) ** 2)

      DELK (I, J, K) = SQRT ((X (I, J, K) - X (I, J, K - 1)) ** 2 +
     >                       (Y (I, J, K) - Y (I, J, K - 1)) ** 2 +
     >                       (Z (I, J, K) - Z (I, J, K - 1)) ** 2)

C     Execution.
C     ----------

C     Zero the three low-end faces (or edges if one plane is specified).

      DO K = K1, K2
         DO J = J1, J2
            S (I1, J, K, 1) = ZERO
         END DO

         DO I = I1, I2
            S (I, J1, K, 2) = ZERO
         END DO
      END DO

      DO J = J1, J2
         DO I = I1, I2
            S (I, J, K1, 3) = ZERO
         END DO
      END DO

C     Set up the low-end edge lines because they are missed by the
C     following loops over most of the low-end faces:

      DO I = I1 + 1, I2
         S (I, J1, K1, 1) =  S (I - 1, J1, K1, 1) + DELI (I, J1, K1)
      END DO

      DO J = J1 + 1, J2
         S (I1, J, K1, 2) =  S (I1, J - 1, K1, 2) + DELJ (I1, J, K1)
      END DO

      DO K = K1 + 1, K2
         S (I1, J1, K, 3) =  S (I1, J1, K - 1, 3) + DELK (I1, J1, K)
      END DO

C     Set up the rest of the low-end face lines because they are
C     missed by the the main loop over most of the volume.

      DO K = K1 + 1, K2
         DO J = J1 + 1, J2
            S (I1, J, K, 2) =  S (I1, J - 1, K, 2) + DELJ (I1, J, K)
            S (I1, J, K, 3) =  S (I1, J, K - 1, 3) + DELK (I1, J, K)
         END DO

         DO I = I1 + 1, I2
            S (I, J1, K, 1) =  S (I - 1, J1, K, 1) + DELI (I, J1, K)
            S (I, J1, K, 3) =  S (I, J1, K - 1, 3) + DELK (I, J1, K)
         END DO
      END DO

      DO J = J1 + 1, J2
         DO I = I1 + 1, I2
            S (I, J, K1, 1) =  S (I - 1, J, K1, 1) + DELI (I, J, K1)
            S (I, J, K1, 2) =  S (I, J - 1, K1, 2) + DELJ (I, J, K1)
         END DO
      END DO

C     Traverse the block just once for all lines except those within
C     the low-end faces.

      DO K = K1 + 1, K2
         DO J = J1 + 1, J2
            DO I = I1 + 1, I2
               S (I, J, K, 1) =  S (I - 1, J, K, 1) + DELI (I, J, K)
               S (I, J, K, 2) =  S (I, J - 1, K, 2) + DELJ (I, J, K)
               S (I, J, K, 3) =  S (I, J, K - 1, 3) + DELK (I, J, K)
            END DO
         END DO
      END DO

C     Normalizing requires another pass through the volume.
C     Degenerate cases must be avoided, so separate the 3 directions:

      DO K = K1, K2
         DO J = J1, J2
            RTOTAL = ONE / MAX (S (I2, J, K, 1), EPS)
            DO I = I1 + 1, I2 - 1  ! No trips if I1 = I2
               S (I, J, K, 1) = S (I, J, K, 1) * RTOTAL
            END DO
            S (I2, J, K, 1) = ONE
         END DO
      END DO

      DO K = K1, K2
         DO I = I1, I2
            RTOTAL = ONE / MAX (S (I, J2, K, 2), EPS)
            DO J = J1 + 1, J2 - 1
               S (I, J, K, 2) = S (I, J, K, 2) * RTOTAL
            END DO
            S (I, J2, K, 2) = ONE
         END DO
      END DO

      DO J = J1, J2
         DO I = I1, I2
            RTOTAL = ONE / MAX (S (I, J, K2, 3), EPS)
            DO K = K1 + 1, K2 - 1
               S (I, J, K, 3) = S (I, J, K, 3) * RTOTAL
            END DO
            S (I, J, K2, 3) = ONE
         END DO
      END DO

      RETURN
      END
