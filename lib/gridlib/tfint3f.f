C+------------------------------------------------------------------------------
C
      SUBROUTINE TFINT3F (IMIN, IMAX, JMIN, JMAX, KMIN, KMAX,
     >                    I1, I2, J1, J2, K1, K2, BF, DIJ, DIK, DJK,
     >                    FUN1, FUN2, FUN3)
C
C  ONE-LINER: 3-space interpolation from faces to interior (3-function case)
C
C  METHOD:
C
C     TFINT3F is the 3-function version of TFINT3D.
C
C     The 3-stage method is based on transfinite interpolation as described
C     in the "projectors" section of "Numerical Grid Generation" by Thompson,
C     Warsi, and Mastin (1985), pp. 315-326:
C
C        (1) Interpolate the function in the I direction (from I face values).
C        (2) Interpolate the deviation of this result from the J faces, in the
C            J direction.
C        (3) Interpolate the combined deviation from the K faces, in the
C            K direction (and combine with the results of the first two stages).
C
C     The 3 stages are merged for efficiency by storing just face values.
C
C  HISTORY:
C
C     9/11/96  SJE  Initial implementation.
C
C  AUTHOR:  Stephen J. Edwards, NASA Ames, Mtn. View, CA.
C
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER IMIN, IMAX, JMIN, JMAX, KMIN, KMAX  ! I Grid array dimensions

      INTEGER I1, I2, J1, J2, K1, K2              ! I Define active volume

      REAL    BF(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX,3) ! I Blending functions in
                                                  !   [0,1] for the I,J,K lines

      REAL    DIJ(2,I1:I2,K1:K2,3),               ! S J face storage for I sweep
     >        DIK(2,I1:I2,J1:J2,3),               ! S K   "   "   "   "  I  "
     >        DJK(2,I1:I2,J1:J2,3)                ! S K   "   "   "   "  J  "

      REAL    FUN1(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX),! I/O Functions being
     >        FUN2(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX),!     interpolated into
     >        FUN3(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) !     the interior

C-------------------------------------------------------------------------------

C     Local constants:

      REAL, PARAMETER :: ONE = 1.

C     Local variables:

      INTEGER I, J, K

      REAL    P1, P2, P3, Q1, Q2, Q3

C     Execution:

C     Calculate the deviations from the I sweep on (just) the J faces:

      DO K = K1, K2
         DO I = I1, I2
            P1 = BF(I,J1,K,1)
            Q1 = ONE - P1

            DIJ(1,I,K,1) = FUN1(I,J1,K) -
     >                     Q1 * FUN1(I1,J1,K) - P1 * FUN1(I2,J1,K)
            DIJ(1,I,K,2) = FUN2(I,J1,K) -
     >                     Q1 * FUN2(I1,J1,K) - P1 * FUN2(I2,J1,K)
            DIJ(1,I,K,3) = FUN3(I,J1,K) -
     >                     Q1 * FUN3(I1,J1,K) - P1 * FUN3(I2,J1,K)

            P1 = BF(I,J2,K,1)
            Q1 = ONE - P1

            DIJ(2,I,K,1) = FUN1(I,J2,K) -
     >                     Q1 * FUN1(I1,J2,K) - P1 * FUN1(I2,J2,K)
            DIJ(2,I,K,2) = FUN2(I,J2,K) -
     >                     Q1 * FUN2(I1,J2,K) - P1 * FUN2(I2,J2,K)
            DIJ(2,I,K,3) = FUN3(I,J2,K) -
     >                     Q1 * FUN3(I1,J2,K) - P1 * FUN3(I2,J2,K)
          END DO
      END DO

C     Deviations from I sweep values for the K faces:

      DO J = J1, J2
         DO I = I1, I2
            P1 = BF(I,J,K1,1)
            Q1 = ONE - P1

            DIK(1,I,J,1) = FUN1(I,J,K1) -
     >                     Q1 * FUN1(I1,J,K1) - P1 * FUN1(I2,J,K1)
            DIK(1,I,J,2) = FUN2(I,J,K1) -
     >                     Q1 * FUN2(I1,J,K1) - P1 * FUN2(I2,J,K1)
            DIK(1,I,J,3) = FUN3(I,J,K1) -
     >                     Q1 * FUN3(I1,J,K1) - P1 * FUN3(I2,J,K1)

            P1 = BF(I,J,K2,1)
            Q1 = ONE - P1

            DIK(2,I,J,1) = FUN1(I,J,K2) -
     >                     Q1 * FUN1(I1,J,K2) - P1 * FUN1(I2,J,K2)
            DIK(2,I,J,2) = FUN2(I,J,K2) -
     >                     Q1 * FUN2(I1,J,K2) - P1 * FUN2(I2,J,K2)
            DIK(2,I,J,3) = FUN3(I,J,K2) -
     >                     Q1 * FUN3(I1,J,K2) - P1 * FUN3(I2,J,K2)
         END DO
      END DO

C     Set up the combined I & J sweep deviations on the K faces for the K pass.
C     (We could avoid DJK altogether by overwriting DIK here, but enough
C     clarity has been lost already.)

      DO J = J1, J2
         DO I = I1, I2
            P2 = BF(I,J,K1,2)
            Q2 = ONE - P2

            DJK(1,I,J,1) = DIK(1,I,J,1) - 
     >                     Q2 * DIJ(1,I,K1,1) - P2 * DIJ(2,I,K1,1)
            DJK(1,I,J,2) = DIK(1,I,J,2) -
     >                     Q2 * DIJ(1,I,K1,2) - P2 * DIJ(2,I,K1,2)
            DJK(1,I,J,3) = DIK(1,I,J,3) -
     >                     Q2 * DIJ(1,I,K1,3) - P2 * DIJ(2,I,K1,3)

            P2 = BF(I,J,K2,2)
            Q2 = ONE - P2

            DJK(2,I,J,1) = DIK(2,I,J,1) -
     >                     Q2 * DIJ(1,I,K2,1) - P2 * DIJ(2,I,K2,1)
            DJK(2,I,J,2) = DIK(2,I,J,2) -
     >                     Q2 * DIJ(1,I,K2,2) - P2 * DIJ(2,I,K2,2)
            DJK(2,I,J,3) = DIK(2,I,J,3) -
     >                     Q2 * DIJ(1,I,K2,3) - P2 * DIJ(2,I,K2,3)
          END DO
      END DO

C     A single pass through the volume completes all 3 sweeps from the
C     appropriate face quantities:

      DO K = K1 + 1, K2 - 1
         DO J = J1 + 1, J2 - 1
            DO I = I1 + 1, I2 - 1
               P1 = BF(I,J,K,1)
               Q1 = ONE - P1
               P2 = BF(I,J,K,2)
               Q2 = ONE - P2
               P3 = BF(I,J,K,3)
               Q3 = ONE - P3

               FUN1(I,J,K) = Q1 * FUN1(I1,J,K) + P1 * FUN1(I2,J,K) +
     >                       Q2 * DIJ(1,I,K,1) + P2 * DIJ(2,I,K,1) +
     >                       Q3 * DJK(1,I,J,1) + P3 * DJK(2,I,J,1)

               FUN2(I,J,K) = Q1 * FUN2(I1,J,K) + P1 * FUN2(I2,J,K) +
     >                       Q2 * DIJ(1,I,K,2) + P2 * DIJ(2,I,K,2) +
     >                       Q3 * DJK(1,I,J,2) + P3 * DJK(2,I,J,2)

               FUN3(I,J,K) = Q1 * FUN3(I1,J,K) + P1 * FUN3(I2,J,K) +
     >                       Q2 * DIJ(1,I,K,3) + P2 * DIJ(2,I,K,3) +
     >                       Q3 * DJK(1,I,J,3) + P3 * DJK(2,I,J,3)
            END DO
         END DO
      END DO

      END SUBROUTINE TFINT3F
