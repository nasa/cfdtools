C+----------------------------------------------------------------------
C
      SUBROUTINE PLXCUT (N, X, Y, T, IX, CALCT, XCUT, TOL, TCUT, LUNOUT,
     >                   IER)
C
C ONE-LINER: Parametric local spline: Find T for give X
C
C DESCRIPTION:
C
C        PLXCUT estimates the value of the parametric variable T where
C     a curve has abscissa XCUT.  It is analogous to the earlier PLBICUT
C     surface utility.  (INTSEC2 could have done the job with a 2-point
C     line, but it fails to return the intersection T values.)
C
C        A safeguarded Newton iteration is used to solve the nonlinear
C     equation  F (T) = X (T) - XCUT = 0.
C
C ENVIRONMENT:  FORTRAN 77, with minor extensions
C
C HISTORY:
C     03/24/95  D.A.Saunders  Initial adaptation of INTSEC2.
C
C AUTHOR: David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C ----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER N            ! (I) Number of points on curve, >= 2.

      REAL    X (N), Y (N) ! (I) Curve coordinates.

      REAL    T (N)        ! (I/S) Cumulative chord lengths, either input
                           !       or computed here - see CALCT.

      INTEGER IX           ! (I) Estimate of index near the point of
                           !     intersection.  If 0, middle point is used.
                           ! (O) Actual index nearest to the point of
                           !     intersection in INTERVAL's sense (that
                           !     is, between 1 and N - 1 inclusive).

      LOGICAL CALCT        ! (I) .FALSE. means the calling program has
                           !     supplied T (*).

      REAL    XCUT         ! (I) X for which T is desired.

      REAL    TOL          ! (I) Tolerance on || f || used relative to
                           !     the curve length.  Try 1.E-6 or 1.E-14
                           !     for 32- or 64-bit arithmetic.

      REAL    TCUT         ! (O) Estimated T corresponding to X = XCUT.

      INTEGER LUNOUT       ! (I) Logical unit for displaying iterations,
                           !     which are suppressed if LUNOUT < 0.
                           !     |LUNOUT| is used for error messages.

      INTEGER IER          ! (O) 0 means no problem was encountered;
                           !     1 means the safeguarded iteration failed
                           !       to satisfy TOL but did satisfy 10*TOL;
                           !     2 means the iteration limit was reached
                           !       without (close) convergence;
                           !     3 means the step-halving failed - somehow
                           !       the iteration diverged;
                           !     4 means the curve was probably parallel
                           !       to X = XCUT.
C     Procedures.

      EXTERNAL  CHORDS2D   ! Chord-length utility.
      EXTERNAL  INTERVAL   ! Search utility, also used by LCSFIT.
      EXTERNAL  LCSFIT     ! Local spline utility.  "Loose" fits are used.

C-----------------------------------------------------------------------

C     Local constants.

      INTEGER   ITMAX
      REAL      HALF, ONE, ZERO
      CHARACTER METHOD * 1
      PARAMETER
     >  (ITMAX  = 15,          ! Outer iteration limit.  (Halving the step
     >   HALF   = 0.5E+0,      !   is the inner iteration, using TOL.)
     >   METHOD = 'B',         ! "Bessel" = "loose" spline fits.
     >   ONE    = 1.E+0,
     >   ZERO   = 0.E+0)

C     Local variables.

      INTEGER   ITER, LUNERR
      REAL      DT, FI, TI, TLAST, XI, XP
      REAL      ALPHA, FNORM, FNORM0, TOLER, TTOTAL
      LOGICAL   NEW

C     Execution.

      LUNERR = ABS (LUNOUT)
      IF (IX .EQ. 0) IX = (N + 1) / 2

C     Set up the cumulative chord lengths unless they are supplied
C     by the higher level for efficiency reasons:

      IF (CALCT) THEN
         CALL CHORDS2D (N, X, Y, .FALSE., TTOTAL, T)
      ELSE
         TTOTAL = T (N)
      END IF

C     Perform a safeguarded Newton iteration:

      IER = 0
      ITER = 0
      TI = T (IX)
      NEW = .TRUE.
      ALPHA = ZERO            ! For printout purposes only
      FNORM0 = 1.E+20         ! I.e., big to avoid iteration 0 test
      TOLER = TOL * TTOTAL

  300 CONTINUE

C        Evaluate x, dx/dt at the current t.

         CALL LCSFIT (N, T, X, NEW, METHOD, 1, TI, XI, XP)

         NEW = .FALSE.
         FI = XI - XCUT
         FNORM = ABS (FI)

C        Halve the step until | f | is reduced (except first time through).

         IF (FNORM .GT. FNORM0) THEN
            IF (ALPHA .GT. TOL) THEN   ! Presume TOL > 2 * machine epsilon
               ALPHA = HALF * ALPHA
               TI = TLAST - ALPHA * DT
               GO TO 300
            END IF
            GO TO 830
         END IF

         IF (LUNOUT .GT. 0) THEN
            WRITE (LUNOUT, '(A, I3, A, 1P, E9.2, A, E9.2, A, E14.6)')
     >         ' PLXCUT:', ITER, '  |f|:', FNORM, '  step:', ALPHA,
     >         '  t:', TI
         END IF

         IF (FNORM .LT. TOLER) GO TO 900  ! Success


         ITER = ITER + 1
         IF (ITER .EQ. ITMAX) GO TO 810

         FNORM0 = FNORM

C        Newton step:

         IF (XP .EQ. ZERO) GO TO 840

         DT = FI / XP
         TLAST = TI
         TI = TI - DT

         ALPHA = ONE

      GO TO 300                           ! Do another iteration


C     Error handling.

  810 IER = 1
      WRITE (LUNERR, 1010) 'Iteration limit reached.'
      IF (FNORM .LT. 10. * TOLER) GO TO 900

      IER = 2
      GO TO 999

  830 IER = 3
      WRITE (LUNERR, 1010) 'Step-halving iteration failed.'
      GO TO 999

  840 IER = 4
      WRITE (LUNERR, 1010) 'Curve parallel to cut?'
      GO TO 999
      

  900 CONTINUE
C     Wrap up, having converged or almost converged.

      TCUT = TI
      CALL INTERVAL (N, T, TI, ONE, IX)

  999 RETURN

C     Formats.

 1010 FORMAT (/, ' PLXCUT: ', A)

      END
