C+----------------------------------------------------------------------
C
      SUBROUTINE QHSFIT (N, X, Y, ENDMODE, B, C, D, IER)
C
C  ACRONYM:  Quasi-Hermite cubic Spline FIT
C            -     -             -      ---
C  DESCRIPTION:
C
C     QHSFIT computes the coefficients B(I), C(I), D(I), I=1,2,...,N
C     for the quasi-Hermite cubic spline described by H. Akima in the
C     Journal of the ACM (Oct. 1970, Vol. 17, No. 4).  This is the
C     algorithm of the IMSL library's IQHSCU subroutine.
C
C     As for any cubic spline, the definition in terms of B, C, D, is
C
C     S(X) = Y(I) + B(I)*(X-X(I)) + C(I)*(X-X(I))**2 + D(I)*(X-X(I))**3
C
C     for X(I) <= X <= X(I+1).  The spline interpolates the data points
C     (X(I), Y(I)) (also called the knots).  Akima's method is a five-
C     point "local" method providing (only) first derivative continuity.
C
C     An Nth set of coefficients is computed, corresponding to X beyond
C     X(N).  (This was needed by one version of CSEVAL and is retained
C     because some applications expect it.)
C
C     The search in CSEVAL now permits decreasing abscissas; the algebra
C     in QHSFIT applies without change in this case.
C
C  METHOD:
C
C     The key idea is to define the spline gradient at the middle point,
C     3, in terms of the slopes m1, m2, m3, m4 of the line segments
C     around point 3 as
C                          |m4 - m3|m2 + |m2 - m1|m3
C                   s3  =  -------------------------
C                            |m4 - m3| + |m2 - m1|
C
C     and to use the average of m2 amd m3 if the denominator goes to 0.
C
C     Such slopes along with the ys at the ends of a given interval
C     define the (Hermite) cubic for that interval.
C
C     If point 3 is at the right end of the dataset, the Akima method
C     introduces two new points 4 and 5 such that
C
C                   x5 - x3  =  x4 - x2  =  x3 - x1
C
C     which allows y3 and y4 to be determined conveniently from the
C     quadratic through points 1, 2, 3 considered as
C
C                   y = y2 + c1 (x - x2) + c2 (x - x2) ** 2
C
C     (meaning all expressions are in terms of the dxs and dys).
C
C     However, this can produce poor results if either of the two end
C     slopes is wild.  The monotonic option should be used for such
C     datasets.
C
C  ARGUMENTS:
C
C     NAME  DIM TYPE I/O/S DESCRIPTION
C     N       -   I    I   Number of data points, or knots.  N >= 2.
C     X       N   R    I   Abscissas of the knots, strictly increasing
C                          or strictly decreasing.
C     Y       N   R    I   Ordinates of the knots.  Y(N)=Y(1) for the
C                          periodic case.
C     ENDMODE -  C*1   I   'P' means a Parabola through the first/last
C                              3 data points controls end behavior (as
C                              in Akima's paper, normally acceptable);
C                          'M' means Monotonicity is forced in the first
C                              two and last two intervals (using the
C                              technique of Brodlie, appropriate where
C                              wild gradients cause 'P' to misbehave);
C                          'C' means CYCLIC (periodic) end conditions
C                              are imposed.
C     B,C,D   N   R    O   Spline coefficients (see above and NOTES).
C     IER     -   I    O   =0: No errors were detected.
C                          =1: Too few data points; N < 2.
C                          =2: Cyclic mode: Y(N) is not equal to Y(1).
C  NOTES:
C     (1)  Y(I) = S(X(I))
C          B(I) = S'(X(I))
C          C(I) = S''(X(I))/2
C          D(I) = S'''(X(I))/6
C
C     (2)  To evaluate the spline, use subroutine CSEVAL.
C
C  PROCEDURES:
C
C     AKIMAG   Akima's formula for the spline gradient at each knot was
C              modularized because of all the end-point cases.
C     BRODLIE  Monotone formula from LCSFIT and PLSFIT, used with
C              THREEPT for gradients at the first and last two data
C              points if ENDMODE = 'M'.
C     HCOEFS   Modular form of the Hermite cubic coefficient formulas.
C     BUTLAND  Forward/backward 3-point formula for first derivative
C              with monotonic adjustment, as for BRODLIE.
C
C     AKIMAG and HCOEFS are in the same module as QHSFIT.
C
C  ENVIRONMENT:  VAX/VMS; FORTRAN 77
C
C  HISTORY:
C   06/06/91  D.A.Saunders   Reimplementation of IQHSCU from Akima's
C                            paper, to bypass the IMSL library.  The
C                            handling of N=2 and 3 and the periodic
C                            case (needed for parametric interpolation
C                            of closed curves) is additional.
C  06/15/91    "      "      Introduced a monotonic end-point option to
C                            overcome the weakness that is possible with
C                            the parabolic extrapolation used by Akima
C                            to add two phantom points off each end.
C  06/20/91    "      "      THREEPT module was renamed BUTLAND.
C
C  AUTHOR:  David Saunders, NASA Ames/Sterling Software, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C ... Arguments:

      INTEGER
     >   N, IER
      REAL
     >   X (N), Y (N), B (N), C (N), D (N)
      CHARACTER
     >   ENDMODE * 1

C ... Local constants:

      REAL
     >   ZERO, ONE, TWO, THREE
      PARAMETER
     >  (ZERO=0.E+0, ONE=1.E+0, TWO=2.E+0, THREE=3.E+0)

C ... Local variables:

      INTEGER
     >   I, I1
      REAL
     >   DXM, DX0, DX1, DX2, DYM, DY0, DY1, DY2, SLM, SL0, SL1, SL2,
     >   DEL (-1 : 1), H (-1 : 1)
      LOGICAL
     >   CYCLIC, MONOTONE

C ... Procedures:

      REAL
     >   AKIMAG, BRODLIE, BUTLAND
      EXTERNAL
     >   AKIMAG, BRODLIE, HCOEFS, BUTLAND

C ... Execution:

C ... Initial error trapping:

      IER = 0
      IF (N .LT. 2) THEN
         IER = 1
         GO TO 99
      END IF

      MONOTONE = ENDMODE .EQ. 'M'
      CYCLIC   = ENDMODE .EQ. 'C'

      IF (CYCLIC) THEN
         IF (Y (N) .NE. Y (1)) THEN
            IER = 2
            GO TO 99
         END IF
      END IF

C ... Initialization.
C     We normally have to introduce points 0 and -1 in that order.

      DX1 = X (2) - X (1)
      DY1 = Y (2) - Y (1)
      SL1 = DY1 / DX1

      IF (N .EQ. 2) THEN
         B (1) = SL1
         C (1) = ZERO
         D (1) = ZERO
         GO TO 50
      END IF

      DX2 = X (3) - X (2)
      SL2 = (Y (3) - Y (2)) / DX2

      IF (CYCLIC) THEN
         DX0 = X (N) - X (N - 1)
         DY0 = Y (N) - Y (N - 1)
         DXM = X (N - 1) - X (N - 2)
         DYM = Y (N - 1) - Y (N - 2)
      ELSE
         DX0 = DX2
         DY0 = DX2 * (TWO * SL1 - SL2)
         DXM = DX1
         DYM = DX1 * (TWO * SL0 - SL1)
      END IF
      SLM = DYM / DXM
      SL0 = DY0 / DX0

      IF (MONOTONE) THEN

C        First derivatives at first 2 points use first 3 data points:

         H (0) = DX1
         H (1) = DX2
         DEL (0) = SL1
         DEL (1) = SL2

         B (1) = BUTLAND (0, H, DEL)
         B (2) = BRODLIE (1, H, DEL)

C        Remaining Hermite coefficients for the first interval:

         CALL HCOEFS (DX1, SL1, B (1), B (2), C (1), D (1))

         I1 = 3
         IF (N .GT. 3) THEN    ! Need to be sure X (4) exists

C           Set up for the normal Akima method starting at the 3rd point:

            DX1 = DX2
            DX2 = X (4) - X (3)
            SL1 = SL2
            SL2 = (Y (4) - Y (3)) / DX2
         END IF

      ELSE                     ! Normal Akima method

C        Since N >= 3, we can now apply the Akima formula at point I = 1:

         B (1) = AKIMAG (SLM, SL0, SL1, SL2)
         I1 = 2

      END IF

C     Main loop (Akima's 5-point formula):

      DO 20, I = I1, N - 2    ! If N = 3, I = N - 1 normally covers I = 2.
         DX0 = DX1
         DX1 = DX2
         DX2 = X (I + 2) - X (I + 1)
         SLM = SL0
         SL0 = SL1
         SL1 = SL2
         SL2 = (Y (I + 2) - Y (I + 1)) / DX2

         B (I) = AKIMAG (SLM, SL0, SL1, SL2)

C        Hermite cubic coefs. for interval between X (I - 1) and X (I).
C        (Do it here to save redoing the bulk of DX, DY, and slope.)

         CALL HCOEFS (DX0, SL0, B (I-1), B (I), C (I-1), D (I-1))

   20 CONTINUE

C     Introduce an (N + 1)th point (not used for monotonic case):

      DX0 = DX1
      DX1 = DX2
      SLM = SL0
      SL0 = SL1
      SL1 = SL2

      IF (CYCLIC) THEN
         DX2 = X (2) - X (1)
         DY2 = Y (2) - Y (1)
      ELSE
         DX2 = DX0
         DY2 = DX0 * (TWO * SL1 - SL0)
      END IF
      SL2 = DY2 / DX2

      IF (MONOTONE) THEN
         H(-1) = DX0
         H (0) = DX1
         DEL(-1) = SL0
         DEL (0) = SL1
         B (N-1) = BRODLIE (0, H, DEL)
         B (N)   = BUTLAND (1, H, DEL)

         CALL HCOEFS (DX0, SL0, B (N-2), B (N-1), C (N-2), D (N-2))
         CALL HCOEFS (DX1, SL1, B (N-1), B (N),   C (N-1), D (N-1))

      ELSE                     ! Normal Akima method

         B (N-1) = AKIMAG (SLM, SL0, SL1, SL2)

         CALL HCOEFS (DX0, SL0, B (N-2), B (N-1), C (N-2), D (N-2))

C        Finally, introduce an (N + 2)th point:

         DX0 = DX1
         SLM = SL0
         SL0 = SL1
         SL1 = SL2

         IF (CYCLIC) THEN
            DX2 = X (3) - X (2)
            DY2 = Y (3) - Y (2)
         ELSE
            DX2 = DX0
            DY2 = DX0 * (TWO * SL1 - SL0)
         END IF
         SL2 = DY2 / DX2

         B (N) = AKIMAG (SLM, SL0, SL1, SL2)

         CALL HCOEFS (DX0, SL0, B (N-1), B (N), C (N-1), D (N-1))

      END IF

   50 CONTINUE

C ... Set coefficients for X > X (N) in case they're ever needed:

      IF (CYCLIC) THEN

         B (N) = B (1)
         C (N) = C (1)
         D (N) = D (1)

      ELSE

C ...    Derivs. at X (N) are the right-end derivs. for X in X (N-1):X (N)

         B (N) = B (N-1) + C (N-1) * DX0 * TWO + D (N-1) * DX0**2 *THREE
         C (N) = C (N-1) + D (N-1) * DX0 * THREE
         D (N) = D (N-1)

      END IF

   99 RETURN
      END
C+----------------------------------------------------------------------
C
      REAL FUNCTION AKIMAG (SL1, SL2, SL3, SL4)
C
C  PURPOSE:
C
C     AKIMAG returns a value to be used as the gradient of a spline at
C     the mid-point of the indicated 5 points, according to the method
C     of Akima.  See QHSFIT for further details.
C
C  ARGUMENTS:
C     NAME  DIM TYPE I/O/S DESCRIPTION
C     SL1,   -    R    I   Slopes of the four line segments surrounding
C     SL2,                 the point in question
C     SL3,
C     SL4
C
C  NOTES:
C     The formula breaks down when both SL1 = SL2 and SL3 = SL4.  (The
C     average of SL2 and SL3 is returned in this case.)  In the presence
C     of round-off, one normally avoids testing against zero, but for
C     extreme cases (very large slopes), the formula appears to degrade
C     quite favorably, even if it SHOULD take the singular branch and
C     average.  But maybe there's room for refinement here.  It seems
C     that the spline could vary discontinuously with the data no matter
C     how clever one is, but with sensible scaling of the data prior
C     to fitting the spline, this should rarely be a problem.
C
C-----------------------------------------------------------------------

C     Arguments:

      REAL
     >   SL1, SL2, SL3, SL4

C     Local constants:

      REAL
     >   HALF, ZERO
      PARAMETER
     >  (HALF = 0.5E+0, ZERO = 0.E+0)

C     Local variables:

      REAL
     >   DENOM, S

C     Execution:

      DENOM = ABS (SL4 - SL3) + ABS (SL2 - SL1)

      IF (DENOM .GT. ZERO) THEN   ! See notes above.
         S = (ABS (SL4 - SL3) * SL2 + ABS (SL2 - SL1) * SL3) / DENOM
      ELSE
         S = (SL2 + SL3) * HALF
      END IF

      AKIMAG = S

      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE HCOEFS (DX, DYDX, SL1, SL2, C, D)
C
C  PURPOSE:
C     HCOEFS modularizes the formulas for the coefficients of the
C     quadratic and cubic terms in the Hermite cubic defined by two
C     function values and two corresponding slopes.  It is needed in
C     this form by the quasi-Hermite spline routine QHSFIT, q.v.
C
C  ARGUMENTS:
C     NAME  DIM TYPE I/O/S DESCRIPTION
C     DX     -    R    I   X length of interval between the two points
C     DYDX   -    R    I   dY/dX for the straight line joining the pts.
C     SL1,   -    R    I   Slopes of the cubic at the two points
C     SL2
C     C,     -    R    O   Quadratic and cubic coefficients
C     D
C
C  NOTES:
C     HCOEFS is in the same source modules as QHSFIT and AKIMAG.
C      
C-----------------------------------------------------------------------

C     Arguments:

      REAL
     >   DX, DYDX, SL1, SL2, C, D

C     Local constants:

      REAL
     >   ONE, TWO, THREE
      PARAMETER
     >  (ONE = 1.E+0, TWO = 2.E+0, THREE = 3.E+0)

C     Local variables:

      REAL
     >   H

C     Execution:

      H = ONE / DX
      C = (THREE * DYDX - TWO * SL1 - SL2) * H
      D = ( -TWO * DYDX +       SL1 + SL2) * H * H

      RETURN
      END
