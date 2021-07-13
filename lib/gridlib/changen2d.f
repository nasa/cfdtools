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
C     07/10/21   "   The original algebraic method is deficient in that, for
C                    instance, doubling the number of cells produces mid-points
C                    of the given cells, which means that the cell growth rate
C                    is stair-stepped.  Fix this simply by invoking the
C                    corrected CHANGEN1D variant, which really does preserve
C                    the relative distribution of the input points.
C
C     Author:  David Saunders, ELORET/NASA Ames Research Ctr., Moffett Field, CA
C              Later with AMA, Inc. at NASA ARC.
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

      CALL CHORDS2D (NDATA, X1(I1), Y1(I1), NORM, TOTAL, T(I1))

      CALL CHANGEN1D (I1, I2, T(I1), IA, IB, TEVAL(IA))

      NINTERP = IB - IA

      CALL LCSFIT (NDATA, T(I1), X1(I1), NEW, METHOD, NINTERP,
     >             TEVAL(IA), X2(IA), DERIV)

      X2(IB) = X1(I2)

      CALL LCSFIT (NDATA, T(I1), Y1(I1), NEW, METHOD, NINTERP,
     >             TEVAL(IA), Y2(IA), DERIV)

      Y2(IB) = Y1(I2)

      END SUBROUTINE CHANGEN2D
