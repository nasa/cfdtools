C+----------------------------------------------------------------------
C
      SUBROUTINE BILINEAR (MODE, IDIM, JDIM, I1, I2, J1, J2, X, Y, F,
     >                     XEVAL, YEVAL, IEVAL, JEVAL, EPS, P, Q,
     >                     FEVAL, FXEVAL, FYEVAL, IER)
C
C        BILINEAR is the bilinear analogue of the bicubic LCSFIT2D, q.v.
C     It performs bilinear interpolation within a regular grid for a
C     discretized function of X and Y at the given point (XEVAL, YEVAL),
C     with the option to estimate partial derivatives.  The formulation
C     for the derivatives allows for arbitrary cell shapes.
C
C        Efficient application to more than one function at a time is
C     provided for.  See LCSFIT2D for further usage details.
C
C     11/22/93  DAS  Initial LCSFIT2D in conjunction with PLBICUBE.
C     07/20/99   "   BILINEAR adapted from LCSFIT2D and BILINT.
C     08/06/99   "   Returned (p,q) for consistency with LCSFIT2D, etc.
C
C     AUTHOR: David Saunders, Raytheon/NASA Ames, Moffett Field, CA.
C
C ----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   MODE,            ! MODE = 0 means return a function value only;
                          !      = 1 means return first partial derivatives too
     >   IDIM,            ! Max. # points provided for in the data arrays.
     >   JDIM,            ! (JDIM is superfluous but retained for symmetry.)
     >   I1,              ! First & last points in each grid direction
     >   I2,              ! eligible for the search and interpolation:
     >   J1,              !    1 <= I1 < I1 + 1 < I2 <= IDIM, etc.
     >   J2               ! Be careful not to include singular grid cells.

      REAL, INTENT (IN), DIMENSION (IDIM, JDIM) ::
     >   X,               ! X (I1:I2, J1:J2) are the meaningful X coordinates
     >   Y,               ! at the data points; likewise for the Y coordinates.
     >   F                ! F may be the third coordinate of a surface or
                          ! it may be some other function of X and Y
      REAL, INTENT (IN) ::
     >   XEVAL,           ! Coordinates of the target point
     >   YEVAL

      INTEGER, INTENT (INOUT) ::
     >   IEVAL,           ! Indices of the cell containing the target point;
     >   JEVAL            ! on input:  the values from a previous call, or
                          !            I1, J1 if no better estimate is known;
                          ! on output: the "lower left" indices of the cell
                          !            containing the target;
                          ! I1 <= IEVAL < I2 on input & output; sim. for JEVAL
      REAL, INTENT (IN) ::
     >   EPS              ! Tolerance used by the search utility - see BILINT

      REAL, INTENT (OUT) ::
     >   P, Q,            ! Fractional cell coords. (p,q) <-> (XEVAL, YEVAL)
     >   FEVAL,           ! Interpolated function value    "      "     "
     >   FXEVAL,          ! Partial df/dx and df/dy (if MODE = 1) "     "
     >   FYEVAL

      INTEGER, INTENT (OUT) ::
     >   IER              ! IER = 0 if the interpolation was successful;
                          !     = 1 if the target point was out of range,
                          !         but results may still be adequate
C     Procedures:

      EXTERNAL
     >   LUSOLVE,         ! For 2 x 2 linear system
     >   RIPPLE2D         ! 2D grid search utility

C-----------------------------------------------------------------------

C     Local constants:

      REAL, PARAMETER ::
     >   ONE = 1.

C     Local variables:

      INTEGER
     >   I, J

      REAL
     >   AJ (2, 2), DP (2), FP, FQ, FPQ, PM1, QM1
 
C     Execution:

C     Identify the cell containing the target point:

      CALL RIPPLE2D (IDIM, JDIM, I1, I2, J1, J2, X, Y, XEVAL, YEVAL,
     >               IEVAL, JEVAL, EPS, P, Q, IER)

C*****IF (IER /= 0) GO TO 99 ! Search failed, but results may be usable

      PM1 = ONE - P
      QM1 = ONE - Q
      I   = IEVAL
      J   = JEVAL

C     Bilinear function interpolation:

      FEVAL = (ONE - Q) * (PM1 * F (I, J)     + P * F (I + 1, J)) +
     >               Q  * (PM1 * F (I, J + 1) + P * F (I + 1, J + 1))
 
C     Interpolated derivatives as well?

      IF (MODE /= 0) THEN ! Adapt the essence of BILINT

         FP  = F (I + 1, J) - F (I, J)
         FQ  = F (I, J + 1) - F (I, J)
         FPQ = F (I + 1, J + 1) - F (I, J + 1) - FP

         DP (1) = FP + Q * FPQ    ! df/dp (RHS vector)
         DP (2) = FQ + P * FPQ    ! df/dq

         FP  = X (I + 1, J) - X (I, J)
         FQ  = X (I, J + 1) - X (I, J)
         FPQ = X (I + 1, J + 1) - X (I, J + 1) - FP

         AJ (1, 1) = FP + Q * FPQ ! dx/dp
         AJ (2, 1) = FQ + P * FPQ ! dx/dq

         FP  = Y (I + 1, J) - Y (I, J)
         FQ  = Y (I, J + 1) - Y (I, J)
         FPQ = Y (I + 1, J + 1) - Y (I, J + 1) - FP

         AJ (1, 2) = FP + Q * FPQ ! dy/dp
         AJ (2, 2) = FQ + P * FPQ ! dy/dq

         CALL LUSOLVE (2, 2, AJ, DP, IER) ! Solve J dfxy = dfpq

         FXEVAL = DP (1)          ! df/dx
         FYEVAL = DP (2)          ! df/dy

      END IF

      END SUBROUTINE BILINEAR
