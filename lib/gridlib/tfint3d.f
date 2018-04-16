C+------------------------------------------------------------------------------
C
      SUBROUTINE TFINT3D (IMIN, IMAX, JMIN, JMAX, KMIN, KMAX,
     >                    I1, I2, J1, J2, K1, K2, S, DIJ, DIK, DJK, F)
C
C  ONE-LINER: 3-space interpolation from faces to interior (1-function case)
C
C  METHOD:
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
C  AUTHOR:  Stephen J. Edwards, NASA Ames, Mt. View, CA.
C
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER IMIN, IMAX, JMIN, JMAX, KMIN, KMAX  ! I Grid array dimensions

      INTEGER I1, I2, J1, J2, K1, K2              ! I Define active volume

      REAL    S(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX,3)  ! I Relative arc-lengths for
                                                  !   the I,J,K lines in [0,1]

      REAL    DIJ(2,IMIN:IMAX,KMIN:KMAX),         ! S J face storage for I sweep
     >        DIK(2,IMIN:IMAX,JMIN:JMAX),         ! S K   "   "   "   "  I  "
     >        DJK(2,IMIN:IMAX,JMIN:JMAX)          ! S K   "   "   "   "  J  "

      REAL    F(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)    ! I/O Function being interpol-
                                                  !     ated into the interior
C-------------------------------------------------------------------------------

C     Local constants.

      REAL    ONE
      PARAMETER (ONE = 1.E+0)

C     Local variables.

      INTEGER I, J, K
      REAL    P1, P2, P3

C     Execution:

C     Calculate the deviations from the I sweep on (just) the J faces:

      DO K = K1, K2
         DO I = I1, I2
            P1 = S(I,J1,K,1)

            DIJ(1,I,K) = F(I,J1,K) - (ONE-P1)*F(I1,J1,K) - P1*F(I2,J1,K)

            P1 = S(I,J2,K,1)

            DIJ(2,I,K) = F(I,J2,K) - (ONE-P1)*F(I1,J2,K) - P1*F(I2,J2,K)
          END DO
      END DO

C     Deviations from I sweep values for the K faces:

      DO J = J1, J2
         DO I = I1, I2
            P1 = S(I,J,K1,1)

            DIK(1,I,J) = F(I,J,K1) - (ONE-P1)*F(I1,J,K1) - P1*F(I2,J,K1)

            P1 = S(I,J,K2,1)

            DIK(2,I,J) = F(I,J,K2) - (ONE-P1)*F(I1,J,K2) - P1*F(I2,J,K2)
         END DO
      END DO

C     Interpolate the combined J & K sweep deviations on (just) the K faces.
C     (We could avoid DJK altogether by overwriting DIK here, but enough
C     clarity has been lost already.)

      DO J = J1, J2
         DO I = I1, I2
            P2 = S(I,J,K1,2)

            DJK(1,I,J) = DIK(1,I,J) -
     >                   (ONE - P2) * DIJ(1,I,K1) - P2  * DIJ(2,I,K1)

            P2 = S(I,J,K2,2)

            DJK(2,I,J) = DIK(2,I,J) -
     >                   (ONE - P2) * DIJ(1,I,K2) - P2  * DIJ(2,I,K2)
          END DO
      END DO

C     A single pass through the volume completes all 3 sweeps from the
C     appropriate face quantities:

      DO K = K1 + 1, K2 - 1
         DO J = J1 + 1, J2 - 1
            DO I = I1 + 1, I2 - 1
               P1 = S(I,J,K,1)
               P2 = S(I,J,K,2)
               P3 = S(I,J,K,3)

               F(I,J,K) = (ONE - P1) * F(I1,J,K)  + P1 * F(I2,J,K)  +
     >                    (ONE - P2) * DIJ(1,I,K) + P2 * DIJ(2,I,K) +
     >                    (ONE - P3) * DJK(1,I,J) + P3 * DJK(2,I,J)
            END DO
         END DO
      END DO

      RETURN
      END
