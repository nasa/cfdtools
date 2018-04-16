C+------------------------------------------------------------------------------
C
      SUBROUTINE INTSEC5 (IDIM, JDIM, I1, I2, J1, J2, XS, YS, ZS,
     >                    US, VS, IS, JS, EPS, K1, K2, XC, YC, ZC,
     >                    TC, KC, TTOTAL, METHOD, CLOSED,
     >                    TINT, UINT, VINT, XINT, YINT, ZINT,
     >                    LUNOUT, IER)
C
C ONE-LINER: INTerSECtion of a curve and a surface in 3-space (bilinear version)
C
C PURPOSE:
C
C        INTSEC5 determines where a curve in 3-space meets a surface.  The
C     curve and surface are defined as discrete points, and extrapolation of
C     the curve (but not the surface) is permitted.
C
C        This is an adaptation of INTSEC4, using parametric bilinear
C     interpolation in place of bicubic, with analytic derivatives from
C     the formulation of BILINT for variation with respect to p and q.
C
C        This version also treats the curve as piecewise linear.  See INTSEC4
C     for further details.
C
C HISTORY:
C     09/23/96  DAS  Adaptation of INTSEC4.
C     09/24/97   "   Replaced PLSCURVE with PLSCRV3D, specifying 'L'
C                    (appropriate for wing surface grids) rather than
C                    adding a METHOD argument.
C     11/28/97   "   Generalized WBINTR forced a METHOD argument.
C     02/06/98   "   Convergence test now uses relative errors in each of x,y,z.
C     08/06/99   "   PBILINT now has analytic derivatives, allowing removal of
C                    the finite differencing originally implemented here.
C
C AUTHOR: David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   IDIM, JDIM,           ! Max. # pts. provided for in the surface arrays
     >   I1, I2, J1, J2        ! Grid index range eligible for searching;
                               ! 1 <= I1 < I1 + 1 < I2 <= IDIM, etc.

      REAL, INTENT (IN), DIMENSION (IDIM, JDIM) ::
     >   XS, YS, ZS            ! Surface grid coordinates

      REAL, INTENT (INOUT), DIMENSION (IDIM, JDIM) ::
     >   US, VS                ! Parametric variables at the surface grid pts.:
                               ! calculated internally if requested via
                               ! US (I1, J1) < 0 on input, else input from a
                               ! prior call to INTSEC5 or PARAM2D (or PLBICUBE)

      INTEGER, INTENT (INOUT) ::
     >   IS, JS                ! Indices of the surface cell containing the
                               ! intersection point:
                               ! on input:  the values from a previous
                               !            call, or I1, J1, if no better
                               !            estimate is known;
                               ! on output: the "lower left" indices of
                               !            the cell (unless IER > 0);
                               ! I1 <= IS < I2 input & output; likewise for JS

      REAL, INTENT (IN) ::
     >   EPS                   ! Tolerance for convergence test here and in
                               ! PBILINT/RIPPLE2D, relative to 1.;
                               ! try 5.*machine eps for 32-bit arithmetic,
                               ! or 1.E-7 should be OK for 64-bit "

      INTEGER, INTENT (IN) ::
     >   K1, K2                ! Index range for the curve arrays: 1 <= K1 < K2

      REAL, INTENT (IN), DIMENSION (K2) ::
     >   XC, YC, ZC            ! Curve coordinates

      REAL, INTENT (INOUT) ::
     >   TC (K2)               ! Parameterization of the curve (normalized
                               ! cumulative chord lengths) - see TTOTAL.
                               ! TC (K1) = 0 and TC (K2) = 1.  TC is an
                               ! argument to avoid likely recalculation.

      INTEGER, INTENT (INOUT) ::
     >   KC                    ! Index of curve interval containing the
                               ! intersection point.  Use K1 if no better
                               ! estimate is known on input.  On output,
                               ! K1 <= KC < K2 by INTERVAL's convention.

      REAL, INTENT (INOUT) ::
     >   TTOTAL                ! Total cumulative chord length of the
                               ! curve (K1 : K2); accompanies TC (*)
                               ! because the normalization loses it.
                               ! Input TTOTAL = 0. if the calling program
                               ! is not supplying TC (*) and TTOTAL, in
                               ! which case INTSEC5 calculates them.

      CHARACTER, INTENT (IN) ::
     >   METHOD * 1            ! Curve fit method; use 'L', 'M', or 'B'
                               ! for piecewise linear, monotonic (tight)
                               ! cubic, or "Bessel" (smooth) cubic resp.

      LOGICAL, INTENT  (IN) ::
     >   CLOSED                ! .TRUE. means the curve is closed, with
                               ! input end pts. matching - see PLSCRV3D.

      REAL, INTENT (INOUT) ::
     >   TINT, UINT, VINT      ! Estimated point of intersection in the
                               ! parameter space.  On input, enter 0. if
                               ! unknown, in which case TC (KC) and
                               ! US/VS (IS, JS) are used as starting
                               ! guesses, else enter the values from a
                               ! previous call.  All are in [0, 1].

      REAL, INTENT (OUT) ::
     >   XINT, YINT, ZINT      ! Estimated point of intersection in real space

      INTEGER, INTENT (IN) ::
     >   LUNOUT                ! Logical unit for displaying iterations,
                               ! which are suppressed if LUNOUT < 0.
                               ! |LUNOUT| is used for error messages.

      INTEGER, INTENT (OUT) ::
     >   IER                   ! 0 means no problem was encountered;
                               ! 1 means trouble in the surface interpolation
                               !   but results may still be usable;
                               ! 2 means the safeguarded iteration failed to
                               !   satisfy EPS but did satisfy 10*EPS;
                               ! 3 means the iteration limit was reached without
                               !   (close) convergence;
                               ! 4 means the step-halving failed - somehow the
                               !   iteration diverged;
                               ! 5 means the linearized system was singular:
                               !   the surface and curve could be parallel.

C     Procedures:

      EXTERNAL
     >   CHORDS3D,             ! Cumulative chord lengths for a line
     >   PARAM2D,              ! Cumulative chord lengths for a surface mesh
     >   PBILINT,              ! Parametric bilinear interpolation
     >   PLSCRV3D,             ! Local cubic spline utility for a 3-space curve
     >   LUSOLVE               ! Square system solver (LU decomposition)

C-------------------------------------------------------------------------------

C     Local constants:

      INTEGER, PARAMETER ::
     >   ITMAX   = 10      ! Outer iteration limit.  (Halving the step
                           ! is the inner iteration.)
      REAL, PARAMETER ::
     >   HALF    = 0.5,
     >   ONE     = 1.0,
     >   STEPMIN = 0.001,  ! Limit on repeated step halvings
     >   ZERO    = 0.0

C     Local variables:

      INTEGER
     >   I, ITER, LUNERR, NPTS

      REAL
     >   A (3, 3), DT (3), F (3), T (3), TLAST (3),
     >   ABSF, ALPHA, FNORM, FNORM0, P, Q, TOLER, XINTC, YINTC, ZINTC,
     >   XINTS, YINTS, ZINTS, XSCALE, YSCALE, ZSCALE

      LOGICAL
     >   NEW


C     Initialization:
C     ---------------

      LUNERR = ABS (LUNOUT)
      NPTS = K2 - K1 + 1
      NEW = .TRUE.     ! For the curve
      A (1, 3) = ZERO  ! -999. would suppress curve derivatives

C     Parameterize the surface?
C     -------------------------

      IF (US (I1, J1) < ZERO) THEN ! Must be a first call

         CALL PARAM2D (IDIM, JDIM, I1, I2, J1, J2, XS, YS, ZS, US, VS)

      END IF


C     Parameterize the curve?
C     -----------------------

      IF (TTOTAL == ZERO) THEN ! Specify normalized chord lengths

         CALL CHORDS3D (NPTS, XC (K1), YC (K1), ZC (K1), .TRUE., TTOTAL,
     >                  TC (K1))
      END IF


C     Safeguarded Newton iteration to solve the 3x3 nonlinear system:
C     ---------------------------------------------------------------

      ITER   = 0
      ALPHA  = ZERO      ! Step length initialized for print purposes only
      FNORM0 = 1.E+20    ! I.e., big to avoid iteration 0 test
      TOLER  = 3.* EPS   ! Used for || f || test; avoids divide by 3

      T (1) = UINT       ! Starting guesses for u, v, t
      T (2) = VINT
      T (3) = TINT
      IF (T (1) == ZERO) T (1) = US (IS, JS)
      IF (T (2) == ZERO) T (2) = VS (IS, JS)      
      IF (T (3) == ZERO) T (3) = TC (KC)

  300 CONTINUE

C        Evaluate the surface and its u/v derivatives at the current (u,v):

         CALL PBILINT (1, IDIM, JDIM, I1, I2, J1, J2, XS, YS, ZS,
     >                 US, VS, T (1), T (2), IS, JS, EPS, P, Q,
     >                 XINTS, YINTS, ZINTS, A (1, 1), A (1, 2), IER)

         IF (IER /= 0) WRITE (LUNERR, 1000) T

C        Evaluate the curve coordinates and derivatives at the current t:

         CALL PLSCRV3D (NPTS, XC (K1), YC (K1), ZC (K1), TC (K1),
     >                  METHOD, NEW, CLOSED, T (3), KC,
     >                  XINTC, YINTC, ZINTC, A (1, 3))

         IF (ITER == 0) THEN  ! Get the scaling from the initial curve point:
            NEW = .FALSE.
            XSCALE = MIN (ONE, ONE / MAX (EPS, ABS (XINTC)))
            YSCALE = MIN (ONE, ONE / MAX (EPS, ABS (YINTC)))
            ZSCALE = MIN (ONE, ONE / MAX (EPS, ABS (ZINTC)))
         END IF

C        Find the norm of the residual vector f, which should converge to ~0.
C        Set up the RHS vector for the system J dT = f in the process.

         F (1) = XINTS - XINTC
         F (2) = YINTS - YINTC
         F (3) = ZINTS - ZINTC
         FNORM = ABS (F (1)) * XSCALE + ABS (F (2)) * YSCALE +
     >           ABS (F (3)) * ZSCALE

C        Halve the step until || f || is reduced (except first time through).
C        See comment below about constraining u and v but not t.

         IF (FNORM >= FNORM0) THEN
            IF (ALPHA > STEPMIN) THEN
               ALPHA = HALF * ALPHA
               T (1) = MAX (ZERO, MIN (ONE, TLAST (1) - ALPHA * DT (1)))
               T (2) = MAX (ZERO, MIN (ONE, TLAST (2) - ALPHA * DT (2)))
               T (3) = TLAST (3) - ALPHA * DT (3)
               GO TO 300
            END IF
            GO TO 830
         END IF

         IF (LUNOUT > 0) THEN
            WRITE (LUNOUT, '(A, I3, A, 1P, E14.6, A, E9.2, A, 3E14.6)')
     >         ' INTSEC5:', ITER, '  ||f||:', FNORM, '  step:', ALPHA,
     >         '  t:', T
         END IF

         IF (FNORM < TOLER) GO TO 900 ! Success


         ITER = ITER + 1
         IF (ITER == ITMAX) GO TO 810

         FNORM0 = FNORM

C        Fix up the LHS matrix for the next iteration.  The first
C        two columns are already the X/Y/Z vectors of surface derivatives:

         A (1, 3) = -A (1, 3)
         A (2, 3) = -A (2, 3)
         A (3, 3) = -A (3, 3)

C        Solve  J dt = f; dt overwrites f.

         CALL LUSOLVE (3, 3, A, F, IER)

         IF (IER /= 0) GO TO 840

         DO I = 1, 3
            DT (I) = F (I)
            TLAST (I) = T (I)
            T (I) = T (I) - DT (I)
         END DO

C        Extrapolation by the curve is permitted:

         T (1) = MAX (ZERO, MIN (ONE, T (1)))
         T (2) = MAX (ZERO, MIN (ONE, T (2)))
         ALPHA = ONE

      GO TO 300                           ! Do another iteration


C     Error handling:
C     ---------------

  810 IER = 2
      WRITE (LUNERR, 1010) 'Iteration limit reached.'
      IF (FNORM < 10.* TOLER) GO TO 900

      IER = 3
      GO TO 999

  830 IER = 4
      WRITE (LUNERR, 1010) 'Step-halving iteration failed.'
      GO TO 999

  840 IER = 5
      WRITE (LUNERR, 1010) 'Singular system. Parallel?'
      GO TO 999
      
  900 CONTINUE
C     Wrap up, having converged or almost converged.
C     Trying to average the surface & curve results somehow is too awkward.

      UINT = T (1)
      VINT = T (2)
      TINT = T (3)
      XINT = XINTS
      YINT = YINTS
      ZINT = ZINTS

  999 RETURN

C     Formats:

 1000 FORMAT (/, ' INTSEC5:  PBILINT target out of range.  Proceeding',
     >        ' with u, v, t =', 3F10.6)
 1010 FORMAT (/, ' INTSEC5: ', A)

      END SUBROUTINE INTSEC5
