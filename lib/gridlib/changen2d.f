C+------------------------------------------------------------------------------
C
      SUBROUTINE CHANGEN2D (I1, I2, X1, Y1, IA, IB, X2, Y2, METHOD)
C
C        CHANGEN2D interpolates X1 & Y1 along a curve between points I1 & I2
C     in order to redistribute them as a different number of points IA : IB
C     with the end points matching.
C
C     12/22/97  DAS  3-space grid line utility CHANGEN wrapped around PLSCRV3D.
C     10/23/02   "   2-space variant: CHANGEN2D.
C
C     Author:  David Saunders, ELORET/NASA Ames Research Ctr., Moffett Field, CA
C
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER    I1, I2     ! I   First & last data points
      REAL       X1(I2),    ! I   Data points
     >           Y1(I2)
      INTEGER    IA, IB     ! I   First & last redistributed points
      REAL       X2(IB),    !   O Redistributed coordinates
     >           Y2(IB)
      CHARACTER  METHOD*1   ! I   Type of fits to be used by LCSFIT, q.v.

C     Procedures:

      EXTERNAL   CHORDS2D,  ! Arc-length utility
     >           LCSFIT     ! 2-space local cubic spline utility

C-------------------------------------------------------------------------------

C     Local constants:

      REAL,    PARAMETER :: ONE = 1.
      LOGICAL, PARAMETER :: NEW = .TRUE., NORM = .FALSE.

C     Local variables:

      INTEGER  I, IP, NDATA, NINTERP
      REAL     T(I1:I2), TEVAL(IA:IB), DERIV(IA:IB),
     >         P, R, RI1, RIP, TOTAL

C     Execution:

      NDATA = I2 - I1 + 1

      CALL CHORDS2D (NDATA, X1(I1), Y1(I1), NORM, TOTAL, T)

      R   = REAL (I2 - I1) / REAL (IB - IA)
      RI1 = REAL (I1)

      DO I = IA, IB - 1
         RIP      = RI1 + R * REAL (I - IA)
         IP       = INT (RIP)
         P        = RIP - REAL (IP)
         TEVAL(I) = (ONE - P) * T(IP) + P * T(IP+1)
      END DO

      NINTERP = IB - IA

      CALL LCSFIT (NDATA, T(I1), X1(I1), NEW, METHOD, NINTERP, TEVAL,
     >             X2(IA), DERIV)

      X2(IB) = X1(I2)

      CALL LCSFIT (NDATA, T(I1), Y1(I1), NEW, METHOD, NINTERP, TEVAL,
     >             Y2(IA), DERIV)

      Y2(IB) = Y1(I2)

      END SUBROUTINE CHANGEN2D
