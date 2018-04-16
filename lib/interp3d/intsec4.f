C+----------------------------------------------------------------------
C
      SUBROUTINE INTSEC4 (IDIM, JDIM, I1, I2, J1, J2, XS, YS, ZS,
     >                    US, VS, IS, JS, ISMAT, JSMAT, SMATRIX, EPS,
     >                    K1, K2, XC, YC, ZC, TC, KC, TTOTAL, METHOD,
     >                    CLOSED, TINT, UINT, VINT, XINT, YINT, ZINT,
     >                    LUNOUT, IER)
C
C ONE-LINER: INTerSECtion of a curve and a surface in 3-space
C
C PURPOSE:
C
C        INTSEC4 determines where a curve in 3-space meets a surface.  The
C     curve and surface are defined as discrete points, and extrapolation of
C     the curve (but not the surface) is permitted.
C
C        The calling routine may perform the arc-length-based parameterization
C     of the surface if that is more efficient (as it would be for multiple
C     calls where the surface is not changing).
C
C METHOD:
C
C        The curve and surface are represented parametrically by "local" splines
C
C           xc = xc (tc),     yc = yc (tc),     zc = zc (tc)        and
C           xs = xs (us, vs), ys = ys (us, vs), zs = zs (us, vs)
C
C     (Local spline coefficients may be calculated as needed from local data,
C     thus avoiding storage problems.  They are normally quite adequate for
C     interpolations.)
C
C        The problem then becomes that of solving a nonlinear system of
C     three equations in three variables t1 = u, t2 = v, t3 = t:
C
C           f1 (t1, t2, t3) = xs (u, v) - xc (t) = 0    and
C           f2 (t1, t2, t3) = ys (u, v) - yc (t) = 0    and
C           f3 (t1, t2, t3) = zs (u, v) - zc (t) = 0    and
C
C     which is solved by a safeguarded Newton iteration.  Note that the
C     partial derivatives of f1/2/3 are functions of the derivatives of
C     xs, ys, zs, and xc, yc which are readily determined by the spline
C     utilities.  Safeguarding involves halving each pure Newton step
C     until || f || is reduced.
C
C        The starting guess is determined from the input values of the
C     parametric variables TINT, UINT, VINT unless they are zero, in
C     which case grid parameter values corresponding to the pointers
C     IS, JS and KC are used.
C
C        The curve is treated as "new" for each call to INTSEC4.  The surface
C     may be treated as new or not new according to US (I1, J1) < 0 or not.
C
C        The standard convergence test on || f || has been adjusted to include
C     the errors in each of X, Y, and Z relative to the magnitudes of the
C     initial curve evaluation.  (Originally, the curve length was used to
C     scale EPS for a tolerance, but the present approach is closer to the
C     ideal of transforming all data to the unit cube, intersecting, then
C     transforming back, which is too cumbersome.)
C
C        This version retries the Newton iteration using finite difference
C     approximations for the surface partial derivatives if the analytic
C     surface derivatives don't work.  (They may be poor if the surface
C     dataset is not smooth.  Smoothness is highly desirable, but finite
C     differencing can sometimes help the iteration avoid getting stuck in
C     a cell making very small steps.  The curve is assumed to be smooth,
C     since this is more easily checked, so its derivatives remain analytic.)
C
C        NOTE:  Since local splines are constructed from finite difference
C     approximations to the derivatives at the data points, the use of
C     finite differences for derivatives at arbitrary target points may
C     be a little confusing.  The hope is that the function values at the
C     target points behave better than the spline derivatives at those
C     points if the underlying data is not as smooth as it should be.
C
C HISTORY:
C     12/05/93  DAS  Initial adaptation of INTSEC2.
C     12/09/93   "   Added TINT, UINT, VINT arguments because
C                    MODGRID4 needs TINT.
C     03/27/94   "   Added the finite-difference retry iteration.
C     01/24/96   "   Use PLBICUBE result even if IER = 1.
C     09/24/97   "   Replaced PLSCURVE with PLSCRV3D, specifying
C                    'M' rather than adding a METHOD argument.
C     11/28/97   "   Generalized WBINTR forced a METHOD argument.
C     02/06/98   "   Convergence test now uses relative errors in x,y,z.
C     08/08/99   "   PLBICUBE now returns (p,q).
C
C AUTHOR: David Saunders, Sterling Software/NASA Ames, Moffett Field, CA.
C
C ----------------------------------------------------------------------

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
                               ! prior call to INTSEC4 or PARAM2D (or PLBICUBE)
      INTEGER, INTENT (INOUT) ::
     >   IS, JS                ! Indices of the surface cell containing the
                               ! intersection point:
                               ! on input:  the values from a previous
                               !            call, or I1, J1, if no better
                               !            estimate is known;
                               ! on output: the "lower left" indices of
                               !            the cell (unless IER > 0);
                               ! I1 <= IS < I2 input & output; likewise for JS

      INTEGER, INTENT (INOUT) ::
     >   ISMAT, JSMAT          ! Target cell indices on a previous call,
                               ! employed for efficiency.  Initialize them
                               ! as I1 - 1, J1 - 1 if US, VS are set up
                               ! externally by a call to PARAM2D; no need
                               ! to initialize them otherwise.  DO NOT
                               ! CHANGE THEIR VALUES BETWEEN CALLS for a
                               ! given surface.
      REAL, INTENT (INOUT) ::
     >   SMATRIX (4, 4, 3)     ! Cell-specific intermediate quantities,
                               ! passed in and out for efficiency
      REAL, INTENT (IN) ::
     >   EPS                   ! Tolerance for convergence test here and in
                               ! PLBICUBE/RIPPLE2D, relative to 1.;
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
                               ! which case INTSEC4 calculates them.
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
                               ! 1 means trouble in the surface interpola-
                               !   tion but results may still be usable;
                               ! 2 means the safeguarded iteration failed
                               !   to satisfy EPS but did satisfy 10*EPS;
                               ! 3 means the iteration limit was reached
                               !   without (close) convergence;
                               ! 4 means the step-halving failed - somehow
                               !   the iteration diverged;
                               ! 5 means the linearized system was singular:
                               !   the surface and curve could be parallel.
C     Procedures:

      EXTERNAL
     >   CHORDS3D,             ! Cumulative chord length utility
     >   PARAM2D,              ! Surface analogue of CHORDS3D
     >   PLBICUBE,             ! Local bicubic spline utility for a surface
     >   PLSCRV3D,             ! Local cubic spline utility for 3-space curve
     >   LUSOLVE               ! Square system solver (LU decomposition)

C-----------------------------------------------------------------------

C     Local constants:

      INTEGER, PARAMETER ::
     >   ITMAX   = 10      ! Outer iteration limit.  (Halving the step
                           ! is the inner iteration.)
      REAL, PARAMETER ::
     >   HALF    = 0.5,
     >   DUV     = 0.0001, ! Differencing interval for approx. u & v derivs.
     >   ONE     = 1.0,
     >   STEPMIN = 0.001,
     >   ZERO    = 0.0

C     Local variables:

      INTEGER
     >   I, J, ILEFT, IRIGHT, IS0, ITER, JLEFT, JRIGHT, JS0, LUNERR,
     >   NPTS

      REAL
     >   A (3, 3), DT (3), F (3), T (3), TLAST (3), UNUSED (3),
     >   UVEVAL (-1:1), XEVAL (-1:1), YEVAL (-1:1), ZEVAL (-1:1),
     >   ABSF, ALPHA, DELUV, FNORM, FNORM0, P, Q, TOLER,
     >   XINTC, YINTC, ZINTC, XINTS, YINTS, ZINTS,
     >   XSCALE, YSCALE, ZSCALE

      LOGICAL
     >   ANALYTIC, NEW


C     Initialization:
C     ---------------

      LUNERR = ABS (LUNOUT)
      NPTS = K2 - K1 + 1
      NEW = .TRUE.      ! For the curve
      A (1, 3) = ZERO   ! -999. would suppress curve derivatives
      ANALYTIC = .TRUE. ! For the surface derivatives
      IS0 = IS
      JS0 = JS

C     Parameterize the surface?
C     -------------------------

      IF (US (I1, J1) < ZERO) THEN   ! Must be a first call

         CALL PARAM2D (IDIM, JDIM, I1, I2, J1, J2, XS, YS, ZS, US, VS)

         ISMAT = I1 - 1  ! So initial SMATRIX won't be used by PLBICUBE
         JSMAT = J1 - 1
      END IF


C     Parameterize the curve?
C     -----------------------

      IF (TTOTAL == ZERO) THEN  ! Specify normalized chord lengths

         CALL CHORDS3D (NPTS, XC (K1), YC (K1), ZC (K1), .TRUE., TTOTAL,
     >                  TC (K1))
      END IF


C     Safeguarded Newton iteration to solve the 3x3 nonlinear system:
C     ---------------------------------------------------------------

  200 CONTINUE
      ITER   = 0
      ALPHA  = ZERO        ! Step length initialized for print purposes only
      FNORM0 = 1.E+20      ! I.e., big to avoid iteration 0 test
      TOLER  = 3.* EPS     ! Used for || f || test; avoids divide by 3

      T (1) = UINT         ! Starting guesses for u, v, t
      T (2) = VINT
      T (3) = TINT
      IF (T (1) == ZERO) T (1) = US (IS, JS)
      IF (T (2) == ZERO) T (2) = VS (IS, JS)      
      IF (T (3) == ZERO) T (3) = TC (KC)

  300 CONTINUE

         IF (ANALYTIC) THEN

C           Evaluate the surface coordinates & derivatives at this (u, v):

            CALL PLBICUBE (1, IDIM, JDIM, I1, I2, J1, J2, XS, YS, ZS,
     >                     US, VS, T (1), T (2), IS, JS, ISMAT, JSMAT,
     >                     SMATRIX, EPS, P, Q, XINTS, YINTS, ZINTS,
     >                     A (1, 1), A (1, 2), IER)

            IF (IER /= 0) WRITE (LUNERR, 1000) T

         ELSE

C           Evaluate the surface, and approximate its derivatives (u first):

            ILEFT = -1                   ! Normally central differences
            IRIGHT = 1
            IF (T (1) == ONE) IRIGHT = 0 ! May have to be one-sided
            IF (T (1) == ZERO) ILEFT = 0

            DO I = ILEFT, IRIGHT

               UVEVAL (I) = T (1) + REAL (I) * DUV

               CALL PLBICUBE (0, IDIM, JDIM, I1, I2, J1, J2, XS, YS, ZS,
     >                        US, VS, UVEVAL (I), T (2), IS, JS, ISMAT,
     >                        JSMAT, SMATRIX, EPS, P, Q, XEVAL (I),
     >                        YEVAL (I), ZEVAL (I), UNUSED, UNUSED, IER)

               IF (IER /= 0) WRITE (LUNERR, 1000) T

            END DO

C           X, Y, Z at current (U, V):

            XINTS = XEVAL (0)
            YINTS = YEVAL (0)
            ZINTS = ZEVAL (0)

C           Approximate U derivatives of X, Y, Z at (U, V):

            DELUV = ONE / (UVEVAL (IRIGHT) - UVEVAL (ILEFT))
            A (1, 1) = (XEVAL (IRIGHT) - XEVAL (ILEFT)) * DELUV
            A (2, 1) = (YEVAL (IRIGHT) - YEVAL (ILEFT)) * DELUV
            A (3, 1) = (ZEVAL (IRIGHT) - ZEVAL (ILEFT)) * DELUV

C           Repeat for V derivatives:

            JLEFT = -1
            JRIGHT = 1
            IF (T (2) == ONE) JRIGHT = 0
            IF (T (2) == ZERO) JLEFT = 0

            DO J = JLEFT, JRIGHT

               UVEVAL (J) = T (2) + REAL (J) * DUV

               CALL PLBICUBE (0, IDIM, JDIM, I1, I2, J1, J2, XS, YS, ZS,
     >                        US, VS, T (1), UVEVAL (J), IS, JS, ISMAT,
     >                        JSMAT, SMATRIX, EPS, P, Q, XEVAL (J),
     >                        YEVAL (J), ZEVAL (J), UNUSED, UNUSED, IER)

               IF (IER /= 0) WRITE (LUNERR, 1000) T

            END DO

C           Approximate V derivatives of X, Y, Z at (U, V):

            DELUV = ONE / (UVEVAL (JRIGHT) - UVEVAL (JLEFT))
            A (1, 2) = (XEVAL (JRIGHT) - XEVAL (JLEFT)) * DELUV
            A (2, 2) = (YEVAL (JRIGHT) - YEVAL (JLEFT)) * DELUV
            A (3, 2) = (ZEVAL (JRIGHT) - ZEVAL (JLEFT)) * DELUV

         END IF

C        Evaluate the curve coordinates and derivatives at the current t:

         CALL PLSCRV3D (NPTS, XC (K1), YC (K1), ZC (K1), TC (K1),
     >                  METHOD, NEW, CLOSED, T (3), KC,
     >                  XINTC, YINTC, ZINTC, A (1, 3))

         IF (ITER == 0) THEN ! Get the scaling from the initial curve point
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
     >         ' INTSEC4:', ITER, '  ||f||:', FNORM, '  step:', ALPHA,
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

C        The surface utilities fail for an out of range u or v,
C        but extrapolation by the curve is permitted:

         T (1) = MAX (ZERO, MIN (ONE, T (1)))
         T (2) = MAX (ZERO, MIN (ONE, T (2)))
         ALPHA = ONE

      GO TO 300                           ! Do another iteration


C     Error handling:
C     ---------------

  810 IER = 2
      WRITE (LUNERR, 1010) 'Iteration limit reached.'
      IF (ANALYTIC) GO TO 890
      IF (FNORM < 10.* TOLER) GO TO 900

      IER = 3
      GO TO 999

  830 IER = 4
      WRITE (LUNERR, 1010) 'Step-halving iteration failed.'
      IF (ANALYTIC) GO TO 890
      GO TO 999

  840 IER = 5
      WRITE (LUNERR, 1010) 'Singular system. Parallel?'
      IF (ANALYTIC) GO TO 890
      GO TO 999
      
  890 CONTINUE
C     Set up for retrying with finite difference derivatives:

      ANALYTIC = .FALSE.
      WRITE (LUNERR, 1020) IS, JS
      IS = IS0
      JS = JS0
      GO TO 200


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

 1000 FORMAT (/, ' INTSEC4:  PLBICUBE target out of range.  Proceeding',
     >        ' with u, v, t =', 3F10.6)
 1010 FORMAT (/, ' INTSEC4: ', A)
 1020 FORMAT (' Retrying with approximate derivatives. i, j: ', 2I6)

      END SUBROUTINE INTSEC4
