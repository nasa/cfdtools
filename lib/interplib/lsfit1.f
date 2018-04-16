C+----------------------------------------------------------------------
C
      SUBROUTINE LSFIT1 (N, X, F, DEGREE, NEARPTS, METHOD, CYCLIC,
     >                   G, B, C, D, IER)
C
C  ACRONYM:  Hybrid Local Least Squares Spline FIT, 1 dimensional data
C                   -     -     -       -      ---  -
C
C  DESCRIPTION:
C
C        For given data points X(I), F(I), I=1,2,...,N, LSFIT1 computes
C     smoothed function values G(I) and (if requested) the coefficients
C     B(I), C(I), and D(I) which define a non-interpolating cubic spline as
C
C     S(X) = G(I) + B(I)*(X-X(I)) + C(I)*(X-X(I))**2 + D(I)*(X-X(I))**3
C
C     for  X(I) <= X <= X(I+1).  (The analysis applies if the X(I) are
C     monotonically DEcreasing as well.)
C
C        Intended for use on noisy data F(I), this hybrid algorithm uses
C     a local least squares polynomial to generate smoothed function values
C     G(I) at each X(I) from the suitably weighted nearby data points (see
C     algorithm outline below).  The second step, involving fitting of a
C     spline to the smoothed values, offers a choice of two similar splines
C     or may be skipped altogether:  METHOD='G' returns just the G(I) values.
C
C        METHOD='H' means the smoothed values (X(I), G(I)) and corresponding
C     first derivatives from the local polynomial fits are used to define a
C     piecewise cubic Hermite spline;  METHOD='M' means just the smoothed
C     function values are used to determine a corresponding monotone spline
C     via MSFIT.  CSEVAL is appropriate for evaluating either of these splines.
C     
C        This class of spline ensures continuity of the first derivative
C     at the knots (X(I), G(I)) but not of the second derivative.  It
C     should be insensitive to the amount of noise in the data.
C
C        A periodic (cyclic) option is provided in order to close curves
C     smoothly as part of parametric fits to data where the X(I) are not
C     necessarily monotonic.
C
C
C  ARGUMENTS:
C     NAME   DIM TYPE I/O/S DESCRIPTION
C     N           I     I   Number of data points;  N >= NEARPTS >= 2
C     X       N   R     I   Abscissas, monotonic increasing or decreasing
C     F       N   R     I   Ordinates (presumably noisy)
C     DEGREE      I     I   Degree of polynomials for local weighted least
C                           squares fits;  DEGREE = 1 (more smoothing)
C                           through 4 (less smoothing);  2 or 3 are normally
C                           the best choices; actual degree used may be
C                           less than the input DEGREE if N is small.
C     NEARPTS     I     I   Number of points in the neighborhood of each
C                           data point to use for each weighted local fit;
C                           normally DEGREE + 2 <= M <= DEGREE + 4 where
C                           M = NEARPTS, unless N is small, when M may be
C                           adjusted so as not to exceed N.  NEARPTS <= 10.
C     METHOD      C*1   I   Controls second step of hybrid method:
C                           'G' suppresses this step, returning just G values;
C                           'H' gives the Hermite spline defined by the
C                               locally smoothed G and G' values;
C                           'M' gives the monotone spline defined by G.
C     CYCLIC  -   L     I   .TRUE. means apply periodic end conditions;
C                           ignored if METHOD='G'
C     G       N   R     O   Smoothed function values corresponding to X(I)
C     B,C,D   N   R     O   Spline coefficients; not used if METHOD='G';
C                           see DESCRIPTION and NOTES
C     IER     -   I     O   0: No errors
C                           1: N < 2                 (too few data points)
C                           2: DEGREE < 1            (should we allow 0 here?)
C                           3  NEARPTS < DEGREE + 2  (no smoothing is feasible)
C
C
C  PROCEDURES:
C
C     MSFIT   Monotone spline fit used for METHOD='M' option
C     PNWFIT  Weighted least squares polynomial fit
C     PNDVAL  Evaluation of polynomial and its derivatives
C
C
C  REFERENCE:
C
C     Foley, Thomas A.  "Scattered data interpolation and approximation
C     with error bounds, Computer Aided Geometric Design 3" (1986) p.163
C
C
C  NOTES:
C
C     (1) The spline, coefficients, and derivatives are related in the usual
C         way as follows:
C
C         G(I) = S (X(I))      (smooth)
C         B(I) = S' (X(I))     (continuous)
C         C(I) = S'' (X(I))/2  (not necessarily continuous)
C         D(I) = S''' (X(I))/6 ( "       "          "     )
C
C     (2) Algorithm outline:
C
C         >  For I = 1,...,N
C            >  Identify m points nearest to X(I)  (m = NEARPTS).
C            >  For each point k of these m points, calculate a weight
C                  w(k) = (1 - s(k)/s)  +  c
C               where s(k) = | X - X(I) | for kth neighbor X,  s = max s(k),
C               and  c ~ 1 / m serves to ensure that all m points contribute
C               non-trivially.
C            >  Solve the weighted least squares polynomial problem.
C            >  Evaluate the polynomial and its first derivative at X(I).
C
C         >  In each interval, the smoothed function value and derivative
C            at the two end-points define a cubic.  Form the coefficients
C            of the piecewise cubic Hermite interpolant to the smoothed data
C            in the usual way, as spline coefficients.
C
C         >  Alternatively, fit a monotone spline to the smoothed function
C            values (ignoring the analytic derivatives from the local poly-
C            nomial fits).  (Or skip the second step altogether.)
C
C     (3) Implementation details:
C
C         >  Since NEARPTS ~ 10 at most, the work-space needed for each
C            local fit is considered small enough to be declared locally.
C
C         >  Degenerate cases are handled as far as possible by using
C            reduced values of DEGREE and NEARPTS if necessary.
C
C
C  ENVIRONMENT:   VAX/VMS FORTRAN 77.  Extensions used:
C
C      IMPLICIT NONE; 8-character names; trailing ! comments.
C
C  HISTORY:
C
C     11/18/87  D. Saunders  Initial implementation of 1-D analog of 2-D
C                            algorithm described in ref. 1 above, with a
C                            view to an eventual 2-D implementation.
C     02/17/89    "    "     Provided option to suppress the second step.
C     05/08/89    "    "     Small values of N handled (belatedly).
C     04/04/90    "    "     Deliberate unused reference to X (N+1) gave
C                            trouble with bounds-checking - put in a MIN.
C     12/04/90    "    "     FLOAT --> REAL (prompted by DOUBLE version).
C     01/28/04    "    "     Increase MAXM from 10 to 100.
C     01/29/04    "    "     WORK dimension = MAXM*(MAXDEG+3), not +4.
C                            See also LSFIT2 for just the 'G' option.
C
C  AUTHOR:  David Saunders, NASA Ames/Sterling Software, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   DEGREE, IER, N, NEARPTS
      REAL
     >   B (N), C (N), D (N), F (N), G (N), X (N)
      CHARACTER
     >   METHOD * 1
      LOGICAL
     >   CYCLIC

C     Local constants:

      INTEGER
     >   LENWORK, MAXDEG, MAXM
      REAL
     >   ONE, TWO, THREE
      PARAMETER
     >  (MAXDEG  = 4,            ! Max. degree of local polynomials
     >   MAXM    = 100,          ! Max. no. of neighbors handled
     >   LENWORK = MAXM * (MAXDEG + 3),
     >   ONE     = 1.E+0,
     >   TWO     = 2.E+0,
     >   THREE   = 3.E+0)

C     Local variables:

      LOGICAL
     >   HERMITE
      INTEGER
     >   DEG, I, J, K, K1, K2, M
      REAL
     >   DFDX, H,
     >   PNCOEFS (0:MAXDEG),    ! Polynomial coefs.
     >   RMSDEV, RSMAX, WMIN, XI,
     >   WEIGHT (MAXM),
     >   WORK ((MAXDEG+3)*MAXM) ! Work-space for solving each
                                ! local overdetermined system

C     Procedures:

      EXTERNAL
     >   MSFIT, PNDVAL, PNWFIT

C     Execution:

      IER = 0
      IF (N .LT. 2)                IER = 1
      IF (DEGREE .LT. 1)           IER = 2  ! 0 works but is rather extreme.
      IF (NEARPTS .LT. DEGREE + 2) IER = 3  ! No smoothing is feasible.
      IF (IER .NE. 0) GO TO 999

C     Allow for small N:

      DEG = MIN (DEGREE, N - 1)
      M = MIN (N, NEARPTS)

      HERMITE = METHOD .EQ. 'H'           ! Only one needing G' evaluated
      WMIN = ONE / MAX (REAL (M), THREE)  ! Should this constant be tuned a bit?


C     Perform local least squares smoothing in neighborhood of each data pt.:

      K1 = 1                    ! Define the nearest M points to X (I).
      K2 = M                    ! These are they for X (1).
      XI = X (1)

      DO 400, I = 1, N

         RSMAX = ONE / MAX (ABS (X (K1) - XI), ABS (X (K2) - XI))
         J = 0

         DO 200, K = K1, K2
            J = J + 1
            WEIGHT (J) = (ONE - RSMAX * ABS (X (K) - XI)) + WMIN
  200    CONTINUE

C        Solve the weighted least squares polynomial problem:

         CALL PNWFIT (M, X (K1), F (K1), WEIGHT, DEG, LENWORK, WORK,
     >                PNCOEFS, RMSDEV, IER)

C        Evaluate the polynomial and its derivative(s) at XI:

         CALL PNDVAL (DEG, PNCOEFS, 1, XI, G (I), H, DFDX)
                                                ! H, DFDX are temporaries
         IF (HERMITE) B (I) = H

C        Adjust for the next polynomial fit:

         XI = X (MIN (I+1, N))          ! Won't be used if I is already N.
  300    CONTINUE
            IF (K2 .LT. N) THEN
               IF (ABS (X (K2 + 1) - XI) .LT. ABS (X (K1) - XI)) THEN
                  K1 = K1 + 1
                  K2 = K2 + 1
                  GO TO 300
               END IF
            END IF

  400 CONTINUE


      IF (METHOD .EQ. 'G') GO TO 999    ! Just smoothed values G(I) reqd.


C     Second part of hybrid method:  spline the smoothed values


C     <The following needs to be beefed up to wrap-around properly>


      IF (CYCLIC) THEN                     ! Average the end pt. values:
         G (1) = (G (1) + G (N)) / TWO
         G (N) = G (1)
         IF (HERMITE) THEN
            B (1) = (B (1) + B (N)) / TWO
            B (N) = B (1)
         END IF
      END IF


      IF (HERMITE) THEN

C        In each interval, the smoothed function value and 1st derivative
C        at the two end-points define a cubic.  Form the quadratic and cubic
C        spline coefficients for the piecewise cubic Hermite interpolant:

         DO 500, I = 1, N-1
            H =  ONE / (X (I+1) - X (I))
            DFDX = H * (G (I+1) - G (I))
            C (I) = (THREE * DFDX - TWO * B (I) - B (I+1)) * H
            D (I) = ( -TWO * DFDX +       B (I) + B (I+1)) * H * H
  500    CONTINUE

      ELSE

C        Fit a monotone spline through the smoothed function values:

         CALL MSFIT (N, X, G, CYCLIC, B, C, D, IER)

      END IF


  999 RETURN   ! End of LSFIT1
      END
