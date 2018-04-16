C+------------------------------------------------------------------------------
C
      SUBROUTINE CHANGEN (I1, I2, X1, Y1, Z1, IA, IB, X2, Y2, Z2,
     >                    METHOD)
C
C        CHANGEN interpolates X1, Y1, Z1 along a curve between points I1 and I2
C     in order to redistribute them as a different number of points IA : IB with
C     the end points matching.
C
C     12/22/97  DAS  3-space grid line utility wrapped around PLSCRV3D.
C     01/11/08   "   Added option to make the output distribution uniform.
C                    See expanded use of METHOD.
C
C     Author:  David Saunders, Sterling Software/NASA Ames, Moffett Field, CA
C
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER    I1, I2     ! I   First & last data points
      REAL       X1(I2),    ! I   Data points
     >           Y1(I2), Z1(I2)
      INTEGER    IA, IB     ! I   First & last redistributed points
      REAL       X2(IB),    !   O Redistributed coordinates
     >           Y2(IB), Z2(IB)
      CHARACTER  METHOD*1   ! I   Type of fit to be used by PLSCRV3D, normally
                            !     M (monotonic), B ("loose"), or L (linear);
                            !     lower case m, b, or l indicates that the
                            !     output line is (essentially) uniform
C     Procedures:

      EXTERNAL   CHORDS3D,  ! Arc-length utility
     >           PLSCRV3D,  ! 3-space local cubic spline utility
     >           UPCASE     ! Character string utility ensuring uppercase

C-------------------------------------------------------------------------------

C     Local constants:

      REAL,    PARAMETER :: DERIVS = -999., ! Suppresses them
     >                      ONE = 1.
      LOGICAL, PARAMETER :: NORM = .FALSE., ! No need to normalize
     >                      CLOSED = .FALSE.
C     Local variables:

      INTEGER    I, IEVAL, IP, NDATA
      REAL       T(I2), DT, P, R, RI1, RIP, TEVAL, TOTAL
      LOGICAL    NEW, UNIFORM
      CHARACTER  MUPPER*1

C     Execution:

      MUPPER = METHOD

      CALL UPCASE (MUPPER)

      UNIFORM = MUPPER /= METHOD

      NDATA = I2 - I1 + 1

      CALL CHORDS3D (NDATA, X1(I1), Y1(I1), Z1(I1), NORM, TOTAL, T(I1))

      IEVAL  = I1
      NEW    = .TRUE.

      X2(IA) = X1(I1)
      Y2(IA) = Y1(I1)
      Z2(IA) = Z1(I1)

      IF (UNIFORM) THEN

         DT = TOTAL / REAL (IB - IA)

         DO I = IA + 1, IB - 1

            TEVAL = DT * REAL (I - IA)

            CALL PLSCRV3D (NDATA, X1(I1), Y1(I1), Z1(I1), T(I1),
     >                     MUPPER, NEW, CLOSED, TEVAL, IEVAL,
     >                     X2(I), Y2(I), Z2(I), DERIVS)
            NEW = .FALSE.

         END DO

      ELSE  ! Preserve original form of spacing as much as possible

         R   = REAL (I2 - I1) / REAL (IB - IA)
         RI1 = REAL (I1)

         DO I = IA + 1, IB - 1

            RIP   = RI1 + R * REAL (I - IA)
            IP    = INT (RIP)
            P     = RIP - REAL (IP)
            TEVAL = (ONE - P) * T(IP) + P * T(IP+1)

            CALL PLSCRV3D (NDATA, X1(I1), Y1(I1), Z1(I1), T(I1),
     >                     MUPPER, NEW, CLOSED, TEVAL, IEVAL,
     >                     X2(I), Y2(I), Z2(I), DERIVS)
            NEW = .FALSE.

         END DO

      END IF

      X2(IB) = X1(I2)
      Y2(IB) = Y1(I2)
      Z2(IB) = Z1(I2)

      END
