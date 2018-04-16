C+----------------------------------------------------------------------
C
      SUBROUTINE INTSEC3 (C, CALCT2, N2, X2, Y2, T2, I2,
     >                    XINT, YINT, TINT, LUNOUT, IER)
C
C ONE-LINER: Intersection of B-spline and local spline curves in 2-space
C
C PURPOSE:
C
C        INTSEC3 determines where a B-spline curve and a curve defined by
C     discrete points in 2-space meet.  The first (B-spline) curve is in
C     DT_NURBS format (polynomial or rational); the second curve is treated
C     by "local" spline techniques, meaning its piecewise cubic spline
C     coefficients are calculable as needed from local data (by four-point
C     formulas).  Extrapolation of this second curve is permitted.
C
C        INTSEC3 was prompted by a need to adjust a 2D mesh following a
C     perturbation of an associated surface such as an airfoil represented
C     by a B-spline curve for design optimization purposes.  Use of B-spline
C     techniques for each mesh line as well would be considerably more
C     cumbersome than the (parametric) local spline techniques of INTSEC2
C     (q.v.), from which INTSEC3 was derived.
C
C METHOD:
C
C        The curves are represented parametrically as
C
C           x1 = x1 (t1), y1 = y1 (t1) and x2 = x2 (t2), y2 = y2 (t2).
C
C        The problem then becomes that of solving a nonlinear system of
C     two equations in two variables:
C
C           f1 (t1, t2) = x1 (t1) - x2 (t2) = 0    and
C           f2 (t1, t2) = y1 (t1) - y2 (t2) = 0
C
C     which is solved by a safeguarded Newton iteration.  Note that the
C     partial derivatives of f1, f2, are functions of the derivatives of
C     x1, y1, x2, y2 which are readily determined by the spline utilities.
C     Safeguarding involves halving each pure Newton step until || f ||
C     is reduced.
C
C        The starting guess is determined from the input values of TINT and
C     I2 - see these arguments below.
C
C        Parallel lines may cause failure, but this should be the only way
C     in which the maximum iteration count could be reached.  Tests for
C     convergence, however, are subject to the choice of tolerance, which is
C     dependent upon the data scaling.  The modest tolerance chosen (1.E-6)
C     is used relative to the arc length of the discrete curve, which is
C     computed for other reasons anyway.  Since the returned point of inter-
C     section is guaranteed to be on the B-spline curve, anything more
C     precise or consistent across calls should not be necessary.
C
C ENVIRONMENT:
C
C     VAX/VMS; FORTRAN 77, with...
C     > IMPLICIT NONE
C     > Trailing ! comments
C     > Names up to 8 characters
C
C HISTORY:
C
C     06/28/93  D.A.Saunders  Initial adaptation of INTSEC2 (dated 10/91),
C                             which treats two discrete-point curves.
C     07/31/93      "         Substituted LUSOLVE for DECOMP & SOLVE.
C
C AUTHOR: David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C ----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      REAL    C (*)            ! (I) B-spline curve (polynomial or rational)
                               !     in DT_NURBS format.

      LOGICAL CALCT2           ! (I) The calling program may well need the
                               !     cumulative chord lengths T2 (*).  Use
                               !     .FALSE. to suppress the calculations here.

      INTEGER N2               ! (I) Number of points in the discrete curve.
                               !     N2 >= 2.

      REAL    X2 (N2), Y2 (N2) ! (I) Coordinates defining the discrete curve.

      REAL    T2 (N2)          ! (S) Workspace for cumulative chord lengths.

      INTEGER I2               ! (I) Estimate of the index near the point of
                               !     intersection on the discrete curve.
                               !     If 0, the middle point is used.
                               ! (O) Actual index nearest to the point of
                               !     intersection in INTERVAL's sense (that
                               !     is, between 1 and N2 - 1 inclusive).

      REAL    XINT, YINT       ! (O) Estimated point of intersection,
                               !     guaranteed to be on the B-spline curve.

      REAL    TINT             ! (I) Estimate of the B-spline curve's dependent
                               !     variable at the point of intersection,
                               !     typically 0. on the first call and the
                               !     output from the previous call thereafter.
                               ! (O) The value corresponding to (XINT, YINT).
      INTEGER LUNOUT           ! (I) Logical unit for displaying iterations,
                               !     which are suppressed if LUNOUT < 0.
                               !     |LUNOUT| is used for error messages.

      INTEGER IER              ! (O) 0 means no problem was encountered;
                               !    -1 means B-spline trouble (further diag-
                               !       nostics provided by DT_NURBS);
                               !     1 means the safeguarded iteration failed
                               !       to satisfy TOL but did satisfy 10*TOL;
                               !     2 means the iteration limit was reached
                               !       without (close) convergence;
                               !     3 means the step-halving failed - somehow
                               !       the iteration diverged;
                               !     4 means the linearized system was singular:
                               !       the lines were probably parallel.
C     Procedures.

      REAL      CHORD
      EXTERNAL  CHORD          ! Chord-length utility.
      EXTERNAL  DTSPDR         ! Derivatives 0, 1, ... for a B-spline curve.
      EXTERNAL  INTERVAL       ! Search utility, also used by LCSFIT.
      EXTERNAL  LCSFIT         ! Local spline utility.  "Loose" fits are used.
      EXTERNAL  LUSOLVE        ! Square system solver (LU decomposition).

C-----------------------------------------------------------------------

C     Local constants.

      INTEGER   ITMAX, MAXDIM, MAXORDER, NT1
      REAL      HALF, ONE, TOL, ZERO
      LOGICAL   NEW
      CHARACTER METHOD * 1
      PARAMETER
     >  (ITMAX  = 30,          ! Outer iteration limit.  (Halving the step
     >   HALF   = 0.5E+0,      !   is the inner iteration, using TOL.)
     >   MAXDIM = 3,           ! 2 for polynomial 2-space curve; 3 for rational.
     >   MAXORDER = 6,         ! Max. order (degree + 1) of the B-spline curve.
     >   NT1    = MAXORDER * (MAXORDER + 4),   ! Enough for 1st derivatives.
     >   METHOD = 'B',         ! "Bessel" = "loose" spline fits.
     >   NEW    = .TRUE.,      ! Always alternating between X1, Y1, etc.
     >   ONE    = 1.E+0,
     >   TOL    = 1.E-6,       ! Tolerance for convergence tests, used
     >   ZERO   = 0.E+0)       ! relative to average arc length (|| f ||)
                               ! or absolutely (halving of step-length).
C     Local variables.

      INTEGER   I, ITER, LUNERR
      REAL      A (2, 2), BSDERIVS (3, 2), DT (2), F (2), T (2),
     >          T1 (NT1), TLAST (2), X (2), XP (2), Y (2), YP (2)
      REAL      ALPHA, FNORM, FNORM0, TOLER

C     Execution.

      LUNERR = ABS (LUNOUT)
      IF (I2 .EQ. 0) I2 = (N2 + 1) / 2

      IF (CALCT2) THEN

C        Set up the cumulative chord lengths for the discrete curve.

         T2 (2) = ZERO
         DO 220, I = 2, N2
            T2 (I) = CHORD (X2, Y2, I - 1, I) + T2 (I - 1)
  220    CONTINUE

      ELSE   ! T2 (*) was input by the calling program
      END IF

C     Perform a safeguarded Newton iteration to solve the 2x2 nonlinear system.

      ITER = 0
      T (1) = TINT
      T (2) = T2 (I2)      
      ALPHA = ZERO            ! For printout purposes only
      FNORM0 = 1.E+20         ! I.e., big to avoid iteration 0 test
      TOLER = TOL * T2 (N2)   ! Scale tolerance by arc length

  300 CONTINUE

C        Evaluate x1, dx1/dt, etc., at the current t = (t(1), t(2)).

         CALL DTSPDR (T (1), 1, C, T1, NT1, BSDERIVS, MAXDIM, IER)
         IF (IER .NE. 0) GO TO 800

         X (1)  = BSDERIVS (1, 1)
         Y (1)  = BSDERIVS (2, 1)
         XP (1) = BSDERIVS (1, 2)
         YP (1) = BSDERIVS (2, 2)

         CALL LCSFIT (N2, T2, X2, NEW, METHOD, 1, T (2), X (2), XP (2))
         CALL LCSFIT (N2, T2, Y2, NEW, METHOD, 1, T (2), Y (2), YP (2))

C        Find the norm of the residual vector f, which should converge to ~0.
C        Set up the RHS vector for the system J dT = f in the process.

         F (1) = X (1) - X (2)
         FNORM = ABS (F (1))
         F (2) = Y (1) - Y (2)
         FNORM = MAX (FNORM, ABS (F (2)))

C        Halve the step until || f || is reduced (except first time through).

         IF (FNORM .GT. FNORM0) THEN
            IF (ALPHA .GT. TOL) THEN   ! Presume TOL > 2 * machine epsilon
               ALPHA = HALF * ALPHA
               T (1) = TLAST (1) - ALPHA * DT (1)
               T (2) = TLAST (2) - ALPHA * DT (2)
               GO TO 300
            END IF
            GO TO 830
         END IF

         IF (LUNOUT .GT. 0) THEN
            WRITE (LUNOUT, '(A, I3, A, 1P, E9.2, A, E9.2, A, 2E14.6)')
     >         ' INTSEC3:', ITER, '  ||f||:', FNORM, '  step:', ALPHA,
     >         '  t:', T
         END IF

         IF (FNORM .LT. TOLER) GO TO 900  ! Success


         ITER = ITER + 1
         IF (ITER .EQ. ITMAX) GO TO 810

         FNORM0 = FNORM

C        Set up the LHS matrix for the next iteration.

         A (1, 1) = XP (1)
         A (2, 1) = YP (1)
         A (1, 2) = -XP (2)
         A (2, 2) = -YP (2)

C        Solve  J dt = f;  dt overwrites f.

         CALL LUSOLVE (2, 2, A, F, IER)
         IF (IER .NE. 0) GO TO 840

         DO 500, I = 1, 2         
            DT (I) = F (I)
            TLAST (I) = T (I)
            T (I) = T (I) - DT (I)
  500    CONTINUE

         ALPHA = ONE

      GO TO 300                           ! Do another iteration


C     Error handling.

  800 WRITE (LUNERR, 1010) 'B-spline curve trouble.'
      IER = -1
      GO TO 999

  810 IER = 1
      WRITE (LUNERR, 1010) 'Iteration limit reached.'
      IF (FNORM .LT. 10. * TOLER) GO TO 900

      IER = 2
      GO TO 999

  830 IER = 3
      WRITE (LUNERR, 1010) 'Step-halving iteration failed.'
      GO TO 999

  840 IER = 4
      WRITE (LUNERR, 1010) 'Singular system. Parallel lines?'
      GO TO 999
      

  900 CONTINUE
C     Wrap up, having converged or almost converged.

      XINT = X (1)
      YINT = Y (1)

      CALL INTERVAL (N2, T2, T (2), ONE, I2)

  999 RETURN

C     Formats.

 1010 FORMAT ('0INTSEC3: ', A)

      END
