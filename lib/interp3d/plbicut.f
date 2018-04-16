C+----------------------------------------------------------------------
C
      SUBROUTINE PLBICUT (KASE, IDIM, JDIM, I1, I2, J1, J2, XS, YS, ZS,
     >                    US, VS, IS, JS, ISMAT, JSMAT, SMATRIX, EPS,
     >                    DRANGE, XYZCUT, UVRANGE, UCUT, VCUT, P, Q,
     >                    LUNOUT, IER)
C
C ONE-LINER: Parametric local bicubic surface cut (Z for given X, Y, etc.)
C
C PURPOSE:
C
C        For a parametric surface x = x(u,v), y = y(u,v), z = z(u,v),
C     PLBICUT determines the point (u,v) and corresponding z which matches
C     a given x and y.  This case and the other two analogous cases are
C     all treated, according to argument KASE.
C
C METHOD:
C
C        For KASE = 1, the problem requires solving the following two
C     nonlinear equations in variables t1 = u, t2 = v:
C
C           f1 (t1, t2) = xs (u, v) - x = 0    and
C           f2 (t1, t2) = ys (u, v) - y = 0
C
C     then (if successful) evaluating zs (u, v) at the solution, which
C     is obtained by a safeguarded Newton iteration.  Note that the
C     partial derivatives of f1/2 are functions of the derivatives of
C     xs and ys which are provided by the spline utilities. Safeguarding
C     involves halving each pure Newton step until || f || is reduced.
C
C        The starting guess is determined from the input values of the
C     arguments UCUT and VCUT, which probably correspond to the pointers
C     IS and JS unless better values (as from a previous call) are known.
C     See INTSEC4 and/or PLBICUBE for further details.
C
C HISTORY:
C
C     01/15/95  DAS  Initial adaptation of INTSEC4.
C     03/09/95   "   Switched from DECOMP/SOLVE to LUSOLVE
C                    and abandoned returning of COND(J).
C     03/17/95   "   Version for use with unnormalized (u,v)s required
C                    inputting the (u,v) data range and abandoning the
C                    option to call PARAM2D.
C     01/24/95   "   Use PLBICUBE results even if IER = 1.
C     02/06/98   "   FNORM is now average of both "f"s;
C                    EPS input matches that for INTSEC4/5.
C     08/06/99   "   The (p,q) fractional coordinates are now returned.
C
C AUTHOR: David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C ----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   KASE                  ! (I)   KASE = 1 means find ZCUT given
                               !                XCUT and YCUT;
                               !       KASE = 2 means find XCUT given
                               !                YCUT and ZCUT;
                               !       KASE = 3 means find YCUT given
                               !                ZCUT and XCUT.
      INTEGER, INTENT (IN) ::
     >   IDIM, JDIM,           ! Max. # pts. provided for in the surface arrays
     >   I1, I2, J1, J2        ! Grid index range eligible for searching;
                               ! 1 <= I1 < I1 + 1 < I2 <= IDIM, etc.

      REAL, INTENT (IN), DIMENSION (IDIM, JDIM) ::
     >   XS, YS, ZS,           ! Surface grid coordinates ...
     >   US, VS                ! ... and parametric variables

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
                               ! as I1 - 1, J1 - 1 on the first call.
                               ! DO NOT CHANGE THEIR VALUES BETWEEN CALLS
                               ! for a given surface.
      REAL, INTENT (INOUT) ::
     >   SMATRIX (4, 4, 3)     ! Cell-specific intermediate quantities,
                               ! passed in and out for efficiency.
      REAL, INTENT (IN) ::
     >   EPS                   ! Tolerance for convergence test here and in
                               ! PLBICUBE/RIPPLE2D, relative to 1.;
                               ! try 5.*machine eps for 32-bit arithmetic,
                               ! or 1.E-7 should be OK for 64-bit "
      REAL, INTENT (IN) ::
     >   DRANGE                ! A measure of the XYZ data range, applied
                               ! to EPS for the convergence test on ||f||
      REAL, INTENT (INOUT) ::
     >   XYZCUT (3)            ! Target point in XYZ space.  XYZCUT(1:2)
                               ! are input, XYZCUT(3) output if KASE = 1;
                               ! other permutations for the other KASEs
      REAL, INTENT (IN) ::
     >   UVRANGE (4)           ! UVRANGE(1:2) = minimum & maximum u;
                               ! UVRANGE(3:4) =    "    "    "    v;
                               ! originally 0 and 1, but now need not be
      REAL, INTENT (INOUT) ::
     >   UCUT, VCUT            ! Estimated target point in the parameter
                               ! space [umin,umax] x [vmin,vmax], usually
                               ! input as US/VS (IS, JS) on a first call,
                               ! then input as previous call outputs. On
                               ! output, the computed u and v.  See IER.
      REAL, INTENT (OUT) ::
     >   P, Q                  ! Fractional cell coords. (p,q) <-> (UCUT,VCUT)

      INTEGER, INTENT (IN) ::
     >   LUNOUT                ! Logical unit for displaying iterations,
                               ! which are suppressed if LUNOUT < 0
      INTEGER, INTENT (OUT) ::
     >   IER                   ! 0 means no problem was encountered;
                               ! 1 means trouble in the surface interpolation
                               !   but results may be usable anyway;
                               ! 2 means the safeguarded iteration failed
                               !   to satisfy EPS but did satisfy 10*EPS;
                               ! 3 means the iteration limit was reached
                               !   without (close) convergence, but the
                               !   solution may still be useful;
                               ! 4 means the step-halving failed - somehow
                               !   the iteration diverged - or the first
                               !   call to PLBICUBE had (u,v) out of range;
                               ! 5 means the linearized system was singular:
                               !   the target X/Y/Z is off the surface.
C     Procedures:

      EXTERNAL
     >   PLBICUBE,             ! Parametric bicubic surface interpation
     >   LUSOLVE               ! Solution of Ax = b via LU decomposition

C-----------------------------------------------------------------------

C     Local constants:

      INTEGER, PARAMETER ::
     >   ITMAX   = 10      ! Outer iteration limit.  (Halving the step
                           ! is the inner iteration.)
      REAL, PARAMETER ::
     >   HALF    = 0.5,
     >   ONE     = 1.0,
     >   STEPMIN = 0.001,
     >   ZERO    = 0.E+0,
     >   BIG     = 1.E+32

C     Local variables:

      INTEGER
     >   ITER, K1, K2, K3

      REAL
     >   A (2, 2), DT (2), XYZU (3), XYZV (3), F (2), T (2), TLAST (2),
     >   XYZS (3), ABSF, ALPHA, FNORM, FNORM0, TOLER

C     Execution:
C     ----------

      K1 = KASE             ! Eqn. 1 involves X if KASE = 1, etc.
      K2 = K1 + 1           ! Eqn. 2 involves Y if KASE = 1, etc.
      IF (K2 == 4) K2 = 1
      K3 = K2 + 1
      IF (K3 == 4) K3 = 1

C     Safeguarded Newton iteration to solve the 2x2 nonlinear system:
C     ---------------------------------------------------------------

      ITER   = 0
      ALPHA  = ZERO         ! Step length initialized for print purposes only
      FNORM0 = BIG          ! To avoid iteration 0 test
      TOLER  = DRANGE * EPS ! Used for || f || test

      T (1) = UCUT          ! Starting guesses for u, v, t
      T (2) = VCUT
C *** IF (T (1) < UVRANGE (1)) GO TO 800   ! These are only rough safeguards
C *** IF (T (1) > UVRANGE (2)) GO TO 800   ! if (u,v) are not normalized to
C *** IF (T (2) < UVRANGE (3)) GO TO 800   ! the [0,1] interval. Experience
C *** IF (T (2) > UVRANGE (4)) GO TO 800   ! indicates normalization is not
                                           ! essential & may be undesirable.
  300 CONTINUE

C        Evaluate the surface coordinates & derivatives at this (u, v):

         CALL PLBICUBE (1, IDIM, JDIM, I1, I2, J1, J2, XS, YS, ZS,
     >                  US, VS, T (1), T (2), IS, JS, ISMAT, JSMAT,
     >                  SMATRIX, EPS, P, Q, XYZS (1), XYZS (2),
     >                  XYZS (3), XYZU, XYZV, IER)

C ***    IF (IER == 1) THEN ... ! Proceed anyway - probably near a boundary

C        Find the residual vector norm, which should converge to ~0.
C        Set up the RHS vector for the system J dT = f in the process.

         F (1) = XYZS (K1) - XYZCUT (K1)
         F (2) = XYZS (K2) - XYZCUT (K2)
         FNORM = (ABS (F (1)) + ABS (F (2))) * HALF

C        Halve the step until || f || is reduced (except first time through).
C        Constrain u and v to [u/vmin, u/vmax].

         IF (FNORM >= FNORM0) THEN
            IF (ALPHA > STEPMIN) THEN
               ALPHA = HALF * ALPHA
               T (1) = MAX (UVRANGE (1),
     >                 MIN (UVRANGE (2), TLAST (1) - ALPHA * DT(1)))
               T (2) = MAX (UVRANGE (3),
     >                 MIN (UVRANGE (4), TLAST (2) - ALPHA * DT(2)))
               GO TO 300
            END IF
            GO TO 830
         END IF

         IF (LUNOUT > 0) THEN
            WRITE (LUNOUT, '(A, I3, A, 1P, E14.6, A, E10.3, A, 2E14.6)')
     >         ' PLBICUT:', ITER, '  ||f||:', FNORM, '  step:', ALPHA,
     >         '  u,v:', T
         END IF

         IF (FNORM < TOLER) GO TO 900 ! Success


         ITER = ITER + 1
         IF (ITER == ITMAX) GO TO 810

         FNORM0 = FNORM

C        Set up the LHS matrix for the next iteration.  The u, v derivatives
C        of functions f1, f2 are those of the surface:

         A (1, 1) = XYZU (K1)
         A (2, 1) = XYZU (K2)
         A (1, 2) = XYZV (K1)
         A (2, 2) = XYZV (K2)

C        Solve  J dt = f; dt overwrites f.

         CALL LUSOLVE (2, 2, A, F, IER)

         IF (IER /= 0) GO TO 840

         DT (1) = F (1)
         TLAST (1) = T (1)
         T (1) = MAX (UVRANGE (1), MIN (UVRANGE (2), T (1) - DT (1)))

         DT (2) = F (2)
         TLAST (2) = T (2)
         T (2) = MAX (UVRANGE (3), MIN (UVRANGE (4), T (2) - DT (2)))
         ALPHA = ONE

      GO TO 300                           ! Do another iteration


C     Error handling:
C     ---------------

C*800 IER = 1       ! PLBICUBE target out of range.
C**** GO TO 999

  810 IER = 2       ! Iteration limit reached.
      IF (FNORM < 10.* TOLER) GO TO 900

      IER = 3
      GO TO 900

  830 IER = 4       ! Step-halving iteration failed
      GO TO 900

  840 IER = 5       ! Singular system; point off the surface?
      GO TO 999
      

  900 CONTINUE
C     Wrap up, having converged or almost converged.

      UCUT = T (1)
      VCUT = T (2)
      XYZCUT (K3) = XYZS (K3)

  999 RETURN

      END SUBROUTINE PLBICUT
