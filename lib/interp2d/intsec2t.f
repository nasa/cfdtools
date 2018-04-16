C---+----------------------------------------------------------------------
C
      SUBROUTINE INTSEC2T (N1, X1, Y1, T1, I1, CALCT1,
     >                     N2, X2, Y2, T2, I2, CALCT2, TOL,
     >                     XINT, YINT, TINT1, TINT2, LUNOUT, IER)
C
C ONE-LINER: INTerSECtion of two curves in 2-space, allowing extrapolation
C
C PURPOSE:
C
C        INTSEC2T is a variant of INTSEC2 that returns the intersection Ts.
C     INTSEC should have done so from the start.
C
C        INTSEC2 determines where two curves in 2-space meet.  The curves
C     are defined as discrete points, and extrapolation of either curve or
C     both is permitted.  It was prompted by a need to sharpen the nose
C     of a rounded airfoil.  The earlier INTSEC does not allow extrapolation,
C     and is less precise in its use of linear interpolation.
C
C        This version allows the calling program to perform the arc-length
C     computations if that is more efficient (as it would be for multiple
C     calls where one of the curves is not changing).
C
C METHOD:
C
C        Each curve is represented parametrically by "local" splines
C
C           x1 = x1 (t1), y1 = y1 (t1) and x2 = x2 (t2), y2 = y2 (t2).
C
C     (Local spline coefficients may be calculated as needed from local data,
C     thus avoiding storage problems.  They are normally quite adequate for
C     interpolations.)
C
C        The problem then becomes that of solving a nonlinear system of
C     two equations in two variables:
C
C           f1 (t1, t2) = x1 (t1) - x2 (t2) = 0    and
C           f2 (t1, t2) = y1 (t1) - y2 (t2) = 0
C
C     which is solved by a safeguarded Newton iteration.  Note that the
C     partial derivatives of f1, f2, are functions of the derivatives of
C     x1, y1, x2, y2 which are readily determined by the spline utility.
C     Safeguarding involves halving each pure Newton step until || f ||
C     is reduced.
C
C        The starting guess is determined from the input values of I1 and
C     I2 (pointers to the data point on each curve as near as is known to
C     the point of intersection).  If 0 is entered for I1 or I2, then the
C     corresponding middle point is used for a starting guess.
C
C        The parametric technique avoids difficulties with vertical lines
C     which would affect INTSEC.  Parallel lines may still cause failure,
C     and this should be the only situation that could cause the maximum
C     iteration count to be reached.  Tests for convergence, however, are
C     subject to the choice of tolerance, which is dependent upon the data
C     scaling.  The input TOL is used relative to the average of the total
C     arc lengths of the two curves.
C
C HISTORY:
C     10/08/91  D.A.Saunders  Initial INTSEC2, for sharpening airfoils.
C     07/31/93       "        Substituted LUSOLVE for DECOMP & SOLVE, and
C                             added the CALCT1, CALCT2 arguments.
C     12/08/93       "        Substituted CHORDS2D for CHORD; made TOL an
C                             argument.
C     02/03/98       "        Changed FNORM .GT. FNORM0 test to .GE.
C     10/25/02       "        Added missing TINT1 & -2 output arguments;
C                             renamed as INTSEC2T.
C
C AUTHOR: David Saunders, ELORET/NASA Ames Research Ctr., Moffett Field, CA
C
C -------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER N1, N2           ! (I) Number of points in each curve, >= 2.

      REAL    X1 (N1), Y1 (N1) ! (I) Coordinates defining the two curves.
      REAL    X2 (N2), Y2 (N2)

      REAL    T1 (N1), T2 (N2) ! (I/S) Cumulative chord lengths, either input
                               !       or computed here - see CALCT1, CALCT2.

      INTEGER I1, I2           ! (I) Estimates of indices near the point of
                               !     intersection.  If 0, middle point is used.
                               ! (O) Actual indices nearest to the point of
                               !     intersection in INTERVAL's sense (that
                               !     is, between 1 and N - 1 inclusive).

      LOGICAL CALCT1, CALCT2   ! (I) .FALSE. means the calling program has
                               !     supplied T1 (*) or T2 (*).

      REAL    TOL              ! (I) Tolerance on || f || used relative to
                               !     the average curve length.  Try 1.E-6 or
                               !     1.E-14 for 32- or 64-bit arithmetic.

      REAL    XINT,  YINT,     ! (O) Estimated point of intersection ...
     >        TINT1, TINT2     ! (O) ... and corresponding Ts

      INTEGER LUNOUT           ! (I) Logical unit for displaying iterations,
                               !     which are suppressed if LUNOUT < 0.
                               !     |LUNOUT| is used for error messages.

      INTEGER IER              ! (O) 0 means no problem was encountered;
                               !     1 means the safeguarded iteration failed
                               !       to satisfy TOL but did satisfy 10*TOL;
                               !     2 means the iteration limit was reached
                               !       without (close) convergence;
                               !     3 means the step-halving failed - somehow
                               !       the iteration diverged;
                               !     4 means the linearized system was singular:
                               !       the lines were probably parallel.
C     Procedures.

      EXTERNAL  CHORDS2D       ! Chord-length utility.
      EXTERNAL  INTERVAL       ! Search utility, also used by LCSFIT.
      EXTERNAL  LCSFIT         ! Local spline utility.  "Loose" fits are used.
      EXTERNAL  LUSOLVE        ! Square system solver (LU decomposition).

C--------------------------------------------------------------------------

C     Local constants:

      INTEGER,   PARAMETER ::
     >   ITMAX    = 12         ! Outer iteration limit.  (Halving the step
                               ! is the inner iteration.)
      REAL,      PARAMETER ::
     >   HALF     = 0.5E+0,
     >   ONE      = 1.E+0,
     >   STEPMIN  = 1.E-3,
     >   ZERO     = 0.E+0

      LOGICAL,   PARAMETER ::
     >   NEW      = .TRUE.,    ! Always alternating between X1, Y1, etc.
     >   NORM     = .FALSE.    ! No need to normalize arc lengths

      CHARACTER, PARAMETER ::
     >   METHOD*1 = 'B'        ! "Bessel" = "loose" spline fits.

C     Local variables:

      INTEGER   I, ITER, LUNERR
      REAL      A (2, 2), DT (2), F (2), T (2), TLAST (2),
     >          X (2), XP (2), Y (2), YP (2)
      REAL      ALPHA, FNORM, FNORM0, TOLER, TTOTAL1, TTOTAL2

C     Execution:

      LUNERR = ABS (LUNOUT)
      IF (I1 == 0) I1 = (N1 + 1) / 2
      IF (I2 == 0) I2 = (N2 + 1) / 2

C     Set up the cumulative chord lengths unless they are supplied
C     by the higher level for efficiency reasons:

      IF (CALCT1) THEN
         CALL CHORDS2D (N1, X1, Y1, NORM, TTOTAL1, T1)
      ELSE
         TTOTAL1 = T1 (N1)
      END IF

      IF (CALCT2) THEN
         CALL CHORDS2D (N2, X2, Y2, NORM, TTOTAL2, T2)
      ELSE
         TTOTAL2 = T2 (N2)
      END IF

C     Perform a safeguarded Newton iteration to solve the 2x2 nonlinear system.

      IER = 0
      ITER = 0
      T (1) = T1 (I1)
      T (2) = T2 (I2)      
      ALPHA = ZERO            ! For printout purposes only
      FNORM0 = 1.E+20         ! I.e., big to avoid iteration 0 test
      TOLER = TOL * HALF * (TTOTAL1 + TTOTAL2)  ! Scale tolerance by avg. length

  300 CONTINUE

C        Evaluate x1, dx1/dt, etc., at the current t = (t(1), t(2)).

         CALL LCSFIT (N1, T1, X1, NEW, METHOD, 1, T (1), X (1), XP (1))
         CALL LCSFIT (N1, T1, Y1, NEW, METHOD, 1, T (1), Y (1), YP (1))
         CALL LCSFIT (N2, T2, X2, NEW, METHOD, 1, T (2), X (2), XP (2))
         CALL LCSFIT (N2, T2, Y2, NEW, METHOD, 1, T (2), Y (2), YP (2))

C        Find the norm of the residual vector f, which should converge to ~0.
C        Set up the RHS vector for the system J dT = f in the process.

         F (1) = X (1) - X (2)
         FNORM = ABS (F (1))
         F (2) = Y (1) - Y (2)
         FNORM = MAX (FNORM, ABS (F (2)))

C        Halve the step until || f || is reduced (except first time through).

         IF (FNORM >= FNORM0) THEN
            IF (ALPHA > STEPMIN) THEN
               ALPHA = HALF * ALPHA
               T (1) = TLAST (1) - ALPHA * DT (1)
               T (2) = TLAST (2) - ALPHA * DT (2)
               GO TO 300
            END IF
            GO TO 830
         END IF

         IF (LUNOUT > 0) THEN
            WRITE (LUNOUT, '(A, I3, A, 1P, E13.6, A, E9.2, A, 2E14.6)')
     >         ' INTSEC2T:', ITER, '  ||f||:', FNORM, '  step:', ALPHA,
     >         '  t:', T
         END IF

         IF (FNORM < TOLER) GO TO 900  ! Success


         ITER = ITER + 1;  IF (ITER == ITMAX) GO TO 810

         FNORM0 = FNORM

C        Set up the LHS matrix for the next iteration.

         A (1, 1) =  XP (1)
         A (2, 1) =  YP (1)
         A (1, 2) = -XP (2)
         A (2, 2) = -YP (2)

C        Solve  J dt = f; dt overwrites f.

         CALL LUSOLVE (2, 2, A, F, IER)
         IF (IER /= 0) GO TO 840

         DO I = 1, 2         
            DT (I) = F (I)
            TLAST (I) = T (I)
            T (I) = T (I) - DT (I)
         END DO

         ALPHA = ONE

      GO TO 300                           ! Do another iteration


C     Error handling:

  810 IER = 1
      WRITE (LUNERR, 1010) 'Iteration limit reached.'
      IF (FNORM < 10. * TOLER) GO TO 900

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

      XINT  = HALF * (X (1) + X (2))
      YINT  = HALF * (Y (1) + Y (2))

      TINT1 = T (1);   CALL INTERVAL (N1, T1, T (1), ONE, I1)
      TINT2 = T (2);   CALL INTERVAL (N2, T2, T (2), ONE, I2)

  999 RETURN

C     Formats:

 1010 FORMAT (/, ' INTSEC2T: ', A)

      END SUBROUTINE INTSEC2T
