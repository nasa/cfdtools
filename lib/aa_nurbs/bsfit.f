C+------------------------------------------------------------------------------
C
      SUBROUTINE BSFIT (NPTS, T, Y, NROWY, NDIM, NDEG, IWT, WHT, ISTART,
     >                  IEND, MXITER, NCPTS, NWORK, WORK, MAXC, C, NC,
     >                  LUNOUT, IER)
C
C     One-liner: B-Spline curve FIT (n-space, with iterative refinement option)
C     ---------  - -            ---
C
C     Purpose:
C     --------
C
C        BSFIT fits a (non-rational) B-spline curve to data points in n-space.
C     It has the option to improve the standard linear least squares fit of
C     the control point coordinates by optimizing the knots as well.  This is
C     done by adjusting the parametric variable (t) values at which the curve
C     corresponds to the data points, deriving knots from those values,
C     recalculating the control points, and repeating until some iteration limit
C     (or goodness-of-fit tolerance) has been met.
C
C        The curve's end points may be constrained to match the end data points.
C
C        This version has a restart capability for two purposes:
C     (a) to continue the iterative refinement (since a few thousand iterations
C         may be desirable); and
C     (b) to obtain a cheap, good fit to a new dataset by reusing the knots and
C         t values from the good fit to some similar dataset.
C
C        Both of these situations have arisen in wing-section applications.
C     In particular, if all sections have the same knots, then conventional
C     lofting techniques can be applied to the control points of the defining
C     sections - a significant simplification compared with general NURBS
C     surface techniques.
C
C     Method:
C     -------
C
C        The initial fit uses the standard global method of the DT_NURBS
C     routine, DTLSA.  The from-scratch choice of knots and t values (left
C     up to the application by DT_NURBS) is heuristic - cumulative chord
C     lengths for the "t"s, and some scheme for deriving knots - see LSKNOTS.
C
C        If the number of control points specified is the same as the number
C     of data points, the fit is interpolatory, with no room for improvement.
C     Otherwise, the limitations caused by the use of cumulative chord length
C     parameterization (that is, by minimizing the (sum of the squares of
C     the) distances between the data points and the points on the spline
C     corresponding to the chord-length-based parameter values of the data
C     points) may be overcome by treating these spline-evaluation points as
C     (nonlinear) variables, and optimizing them to improve the fit.
C
C                                                                  3N      2
C        In 3-space, the function being minimized here is  F(t) = SUM f (t)
C                                                                 i=1  i
C
C     where the first N functions f(t) are the x-deviations of the spline
C     from the data points, the second N are the y-deviations, and so on.
C     The N values of t at which to match the spline to the data points
C     as well as possible are the N nonlinear variables.
C
C        This initial implementation makes use of the Gauss-Newton method,
C     which is well suited to sum-of-squares functions where the residuals
C     at the minimum are small.  (It approximates the Hessian matrix of
C     second derivatives of the full Newton method with just the part
C     involving first derivatives, namely 2J'J, where J is the Jacobian
C     matrix of the f(t) functions being squared and summed.)
C
C        It turns out that J in this case is composed of diagonal blocks
C     (one block for each dimension x, y, z, ...), with diagonal element i
C     of block j being the derivative of the spline's jth coordinate
C     evaluated at the current value of the parameter value associated
C     with data point i.
C
C        Singularity in J'J cannot be a problem, since not all partial
C     derivatives of x, y(, z) with respect to t can be zero at the same
C     time for a sensible curve.  [Would a cusp give trouble?]
C
C        The iteration is identical to that for the Newton solution of a
C     square system of nonlinear equations.  The main difference is that the
C     system is not square, and the norm of g = 2J'f is being driven to zero
C     rather than the norm of f.  The overdetermined system, Jp ~ -f, is
C     solved in the Normal Equations form J'Jp = -J'p for the search direction
C     because this generalizes nicely to any number of dimensions.
C
C        The iterative steps are safeguarded by halving each computed step
C     p until || f || is reduced from its value for the previous iteration.
C     (This is preferable to using || g || because the latter doesn't
C     necessarily decrease every step even for Newton directions. Thus the
C     update step is t <-- t + alpha p, where alpha <= 1.  Some day, a
C     more clever line search might be incorporated (but that is nontrivial
C     in itself).
C
C        Updating the control points every m iterations might offer a
C     savings for some optimal choice of m, which would surely differ from
C     application to application.  Experience indicates that this strategy
C     is more inclined to lead to a spurious local minimum for the nonlinear
C     problem - such as a solution with the deviations in x, say, as close
C     to zero as possible but the deviations in y still poor.  It appears
C     better to be as conservative as possible, even though perturbing the
C     problem by recomputing the knots and control points at every step
C     hardly helps the behavior of the Gauss-Newton iteration.
C
C        Tests for convergence are subject to the choice of tolerance, which
C     is dependent upon the data scaling.  The tolerance chosen (1.E-7 for the
C     2-norm of the gradient at a minimum with 64-bit precision) is used
C     relative to the (presumably) arc-length-based parameter value associated
C     with the last data point.
C
C     Environment:     FORTRAN 77, with...
C     ------------        > IMPLICIT NONE
C                         > Trailing ! comments
C                         > Names up to 8 characters
C
C     History:
C     --------
C
C     04/17/92  D.A.Saunders  Initial design and implementation, with good
C                             B-spline representations of airfoils in mind.
C                             (We need to keep the number of control points
C                             down for design-by-optimization applications.)
C
C     04/26/92     "     "    Gave up on efforts to accelerate convergence
C                             by doubling or incrementing alpha each step.
C                             Typical minima seem extremely flat.  But even
C                             unconverged solutions are normally much better
C                             than the initial standard linear solution.
C
C     03/05/93     "     "    Incorporated a restart capability for the two
C                             situations described above.
C
C     01/23/94     "     "    Normalized the knots in keeping with convention.
C                             (CHORDS2D/3D needed new TOTAL argument besides.)
C
C     02/01/94     "     "    Dispensed with DOUBLE PRECISION.  Compile with
C                             /REAL_LENGTH=64 or -r8 switches as necessary.
C
C     Author:  David Saunders, Sterling Software/NASA Ames Research Center,
C     -------                                    Mt. View, CA.
C
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER NPTS             ! (I)   Number of data points, >= NKNOTS + K.

      REAL    T (NPTS)
                               ! (S/O) Work-space for the values of the
                               !       parametric variable corresponding
                               !       to the data points.

      INTEGER NROWY            ! (I)   Row dimension of Y, >= NPTS.

      INTEGER NDIM             ! (I)   Number of coordinate dimensions, >= 1.

      REAL    Y (NROWY, NDIM)  ! (I)   Data point coordinates X, Y, ... in
                               !       columns corresponding to T (*).

      INTEGER NDEG             ! (I)   Degree of the spline, >= 0. K = NDEG + 1.

      INTEGER IWT              ! (I)   Weight option:
                               !       0 means 1. is used for all weights and
                               !         WHT (*) is ignored;
                               !       1 means WHT (*) values are used.

      REAL    WHT (NPTS)       ! (I)   Non-negative values used for weighting
                               !       the least squares fits (not to be
                               !       confused with the weights of rational
                               !       B-splines).  It IWT = 1, at least
                               !       NKNOTS + K of the WHT (*) values must
                               !       be positive.

      INTEGER ISTART           ! (I)   (Re)start option:
                               !       0 means use chord-length-based T values
                               !         and knots derived from them;
                               !       1 means start with the T values input
                               !         in T (*) and with the knots input in
                               !         the standard portion of C (*).
                               !       If MXITER = 0, these are unchanged.

      INTEGER IEND             ! (I)   End point option:
                               !       0 means the end-points of the curve are
                               !         not constrained - they are calculated
                               !         along with the rest of the control
                               !         points (a.k.a as spline coefficients).
                               !       1 means the end-points of the curve are
                               !         constrained as the end data points.

      INTEGER MXITER           ! (I)   Iterative refinement option:
                               !       0 means no iterative refinement is
                               !         requested - just use the basic global
                               !         least squares method (which may in
                               !         fact be interpolation if NCPTS = NPTS).
                               !      >0 means perform this many refinement
                               !         iterations (or fewer if convergence
                               !         to a minimum is achieved).  If NCPTS
                               !         = NPTS, though, MXITER is ignored.
                               !         Large numbers of iterations are
                               !         possible (in the thousands).  See
                               !         the ISTART argument.

      INTEGER NCPTS            ! (I)   Number of control points specified.
                               !       NDEG + 1 <= NCPTS <= NPTS.

      INTEGER NWORK            ! (I)   Length of WORK (*) provided. <TBD>

      REAL    WORK (NWORK)     ! (S)   Work-space needed by DTLSA and for
                               !       the iterative refinement.

      INTEGER MAXC             ! (I)   Storage allocated to spline vector C (*).

      REAL    C (MAXC)         ! (I/O) Computed spline vector in DT_NURBS form.
                               !       Input with the appropriate header and
                               !       knots if ISTART = 1, else output only.

      INTEGER NC               ! (O)   Length of spline vector C (*):
                               !             NC = 5 + (NDIM + 1) * NCOEF + K
                               !       where NCOEF = NKNOTS + K
                               !       where NKNOTS = number of INTERIOR knots
                               !       and   K = NDEG + 1

      INTEGER LUNOUT           ! (I)   Logical unit for displaying iterations,
                               !       which are suppressed if LUNOUT < 0.
                               !       |LUNOUT| is used for error messages.

      INTEGER IER              ! (O)   Success/error code as returned by DTLSA
                               !       with the following additional values
                               !       possible from the refinement iteration:
                               !        0 means no problem was encountered;
                               !       12 means the iteration limit was reached
                               !          without convergence; the result will
                               !          normally be usable still;
                               !       13 means the step-halving failed to
                               !          reduce || f ||;
                               !       14 means the step-halving failed to
                               !          keep the Ts monotonic;
C***** This shouldn't happen.  !       15 means that J'J was singular.

C     Procedures.

      EXTERNAL  CHORDS2D       ! Cumulative chord lengths in 2-space.
      EXTERNAL  CHORDS3D       ! Cumulative chord lengths in 3-space.
      EXTERNAL  DTLSA          ! Weighted global least sqrs, B-spline curve fit.
      EXTERNAL  DTSPDR         ! B-spline derivatives of order >= 0.
      EXTERNAL  LSKNOTS        ! Heuristic derivation of knots interleaving
                               ! parameter values associated with data points.

C-------------------------------------------------------------------------------

C     Local constants.

      REAL
     >   HALF, ONE, TOL, ZERO
      LOGICAL
     >   NORMLIZ
      PARAMETER
     >  (HALF    = 0.5,
     >   NORMLIZ = .TRUE.,      ! Knots
     >   ONE     = 1.0,
     >   TOL     = 0.0000001,   ! Tolerance for convergence tests, used
     >   ZERO    = 0.0)         ! relative to total arc length (|| g ||)
                                ! or absolutely (halving of step-length).
C     Local variables.

      INTEGER
     >   I, IFAIL, J, ITER, LUNERR, M, NM1, NKNOTS
      REAL
     >   ALPHA, FI, GI, FNORM, FNORM0, GNORM, TOLER, TTOTAL

C     TEMPORARY LOCAL STORAGE TILL THE METHOD IS KNOWN TO BE VALID.
C     -------------------------------------------------------------

      INTEGER
     >   MXPTS, MXDIM
      PARAMETER
     >  (MXPTS = 256, MXDIM = 2)

      REAL
     >   ALHS (MXPTS), BRHS (MXPTS), DERIVS (MXDIM, 2), P (MXPTS),
     >   TLAST (MXPTS)

      LOGICAL
     >   MONO, RESTART

C     Execution.

      LUNERR = ABS (LUNOUT)
      IF (NPTS .GT. MXPTS) GO TO 800

      RESTART = ISTART .NE. 0

      IF (.NOT. RESTART) THEN

C        Compute the chord lengths associated with the data points:

         IF (NDIM .EQ. 2) THEN
            CALL CHORDS2D (NPTS, Y (1, 1), Y (1, 2), NORMLIZ, TTOTAL, T)

         ELSE IF (NDIM .EQ. 3) THEN
            CALL CHORDS3D (NPTS, Y (1, 1), Y (1, 2), Y (1, 3), NORMLIZ,
     >                     TTOTAL, T)
         ELSE
            WRITE (LUNERR, 1010) 'MXDIM is too small.'
            STOP ' '
         END IF
      END IF

      ITER = 0
      ALPHA = ZERO            ! For printout purposes only
      FNORM0 = 1.E+20         ! I.e., big to avoid iteration 0 test
      TOLER = TOL * TTOTAL    ! Scale tolerance by total chord length
      P (1) = ZERO
      P (NPTS) = ZERO
      M = (NPTS +1) / 2
      NM1 = NPTS - 1
      MONO = .TRUE.

  300 CONTINUE

         IF (RESTART) THEN

            RESTART = .FALSE.          ! Update the knots once we're iterating
            NKNOTS = NCPTS - NDEG - 1  ! # interior knots

         ELSE

C           Derive interior knots from the data-point-related parameter values.
C           Put them in their eventual place (even though DTLSA makes a
C           redundant copy into C (*)):

            CALL LSKNOTS (NDEG, NPTS, NCPTS, T, NKNOTS, C (7 + NDEG),
     >                    IER)
         END IF

C        Compute the best control points for the current knots and Ts
C        (global linear least squares method):

         CALL DTLSA (NPTS, T, Y, NROWY, NDIM, NDEG, IWT, WHT,
     >               IEND, C (7 + NDEG), NKNOTS, WORK, NWORK,
     >               MAXC, C, NC, IFAIL, IER)
         IF (IER .NE. 0) GO TO 999               ! DTLSA gives diagnostics

         IF (MXITER .EQ. 0 .OR. NCPTS .EQ. NPTS) GO TO 999


C        Perform a safeguarded Newton-type iteration to improve the fit by
C        optimizing the values of t at which the spline matches the data pts.

C        Evaluate first derivatives w.r.t. t of the f(I) at the current T (*).
C        These are the diagonals of the blocks of the Jacobian J.

         FNORM = ZERO
         GNORM = ZERO

         DO 310, I = 1, NPTS

C           DTSPDR evaluates the spline as well as its 1st deriv. at one U:

            CALL DTSPDR (T (I), 1, C, WORK, NWORK, DERIVS, MXDIM, IER)
            IF (IER .NE. 0) THEN
               WRITE (LUNERR, 1010) 'Bad return from DTSPDR.'
               GO TO 999
            END IF

C           Find the norm of the gradient vector ~J'f, which should converge
C           to ~0.  Set up the (diagonal) system J'J p = -J'f in the process.

            ALHS (I) = ZERO
            BRHS (I) = ZERO
            DO 305, J = 1, NDIM
               FI = DERIVS (J, 1) - Y (I, J)
               GI = DERIVS (J, 2)
               ALHS (I) = ALHS (I) + GI * GI
               BRHS (I) = BRHS (I) - GI * FI
               FNORM = FNORM + FI ** 2
  305       CONTINUE

            GNORM = GNORM + BRHS (I) ** 2

  310    CONTINUE

         GNORM = SQRT (GNORM)
         FNORM = SQRT (FNORM)

C        Halve the step until || f || is reduced (except first time through).

         IF (FNORM .GT. FNORM0) THEN
            IF (ALPHA .GT. TOL) THEN   ! Presume TOL > 2 * machine epsilon
               ALPHA = HALF * ALPHA
               DO 320, I = 2, NM1
                  T (I) = TLAST (I) + ALPHA * P (I)
  320          CONTINUE

C              Monotonicity has already been checked below.

               GO TO 300    ! Update the spline as well.
            END IF
            GO TO 830
         END IF

         IF (LUNOUT .GT. 0) THEN
            WRITE (LUNOUT,
     >         '(A,I5,A,1P,E10.3,A,E9.2,A,E9.2,A,2E12.5,A,2E10.2)')
     >         ' Itn.', ITER, ' ||f||', FNORM, ' ||g||', GNORM, ' step',
     >         ALPHA, ' T2&M', T (2), T (M), ' P2&M', P (2), P (M)
         END IF

         IF (GNORM .LT. TOLER) GO TO 900  ! Success


         ITER = ITER + 1
         IF (ITER .GT. MXITER) THEN  ! Too many iterations needed, but
            ITER = MXITER            ! solution should be OK anyway.
            GO TO 810
         END IF

         FNORM0 = FNORM

C        Gauss-Newton step:  Solve  Jp ~ -f using the Normal Equations.

         DO 400, I = 2, NM1
C*****      IF (ABS (ALHS (I)) .LT. TOL) GO TO 850  ! Not all derivs. can be 0.
            P (I) = BRHS (I) / ALHS (I)
            TLAST (I) = T (I)
            T (I) = T (I) + P (I)
  400    CONTINUE

C        We can't use this T (*) vector if it's not monotonic.

         ALPHA = ONE
  500    CONTINUE

            MONO = .FALSE.
            DO 510, I = 2, NPTS
               IF (T (I) .LT. T (I - 1) + TOLER) GO TO 520
  510       CONTINUE
            MONO = .TRUE.
            GO TO 600

  520       ALPHA = HALF * ALPHA
            IF (ALPHA .LT. TOL) GO TO 830

            DO 530, I = 2, NM1
               T (I) = TLAST (I) + ALPHA * P (I)
  530       CONTINUE

         GO TO 500

  600    CONTINUE
         
      GO TO 300     ! Do another iteration, including new knots and ctl. pts.


C     Error handling.

  800 WRITE (LUNERR, 1010) 'Too many data points for local storage.'
      STOP 'Aborting.'
      
  810 IER = 12
      WRITE (LUNERR, 1010) 'Iteration limit reached.'
      GO TO 900

  830 IF (MONO) THEN
         IER = 13
         WRITE (LUNERR, 1010) 'Step-halving failed to reduce || f ||.'
      ELSE
         IER = 14
         WRITE (LUNERR, 1010) 'Step-halving failed to keep T monotonic.'
         IF (LUNOUT .GT. 0) WRITE (LUNOUT, '(/, A, /, 1P, (5D15.8))')
     >      ' Bad T values:', (T (I), I = 1, NPTS)
      END IF
      GO TO 900

C***  850 IER = 15
C***      WRITE (LUNERR, 1010) 'Singular system. ...'
C***      GO TO 999
      

  900 IF (LUNOUT .LT. 0) THEN  ! Show at least the final iteration.
         WRITE (LUNERR,
     >      '(/,A,/,A,I4,A,1P,E10.3,A,E9.2,A,E9.2,A,2E12.5,A,2E10.2)')
     >      ' BSFIT:', ' It.', ITER, ' ||f||', FNORM, ' ||g||', GNORM,
     >      ' step', ALPHA, ' T2&M', T (2), T (M), ' P2&M', P (2), P (M)
      END IF

  999 RETURN

C     Formats.

 1010 FORMAT (/, ' BSFIT: ', A)

      END
