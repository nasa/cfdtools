C+----------------------------------------------------------------------
C
      SUBROUTINE HTDIS4 (DSINPUT, XA, XB, D1, D2, N, X, LUNOUT, IER)
C
C  ACRONYM:  Hyperbolic Tangent-type DIStribution, 2-sided, 4th version
C            -          -            ---                    -
C  PURPOSE
C  -------
C
C     HTDIS4 is the version of HTDIS2 intended for 64-bit systems.
C     It differs only in its internal procedure, HTDIS3, which need
C     not use DOUBLE PRECISION the way it does for HTDIS2.
C
C     The remaining documentation is that of HTDIS2.
C
C     HTDIS2 generates abscissas on the interval [XA, XB] using the
C     method of Marcel Vinokur to achieve asymmetric bunching.  The
C     input D1 and D2 may be either the desired initial and final
C     INCREMENTS ("deltas") if DSINPUT is .TRUE., or the desired
C     end-point SLOPES of the (normalized) stretching function if
C     DSINPUT is .FALSE.
C
C  METHOD
C  ------
C
C     The mathematics appear in the paper
C
C        "On One-Dimensional Stretching Functions Functions for Finite-
C        Difference Calculations" by Marcel Vinokur,
C        Journal of Computational Physics, Vol. 50, No. 2, May 1983.
C
C     The paper develops criteria for stretching functions that optimize
C     in some sense the truncation errors inherent in finite-difference
C     approximations to derivatives.  The simplest "universal function"
C     (of which a scaled portion is employed) satisfying the criteria is
C     w = tan z, where  z = x + iy is complex.
C
C     The analysis uses the quantity B = SQRT (S0 * S1), where S0 and S1
C     are dimensionless end-point SLOPES of the stretching function.
C     I.e.,  S0 = dXI / dT  at  T = 0,  and  S1 = dXI / dT  at  T = 1,
C     where XI and T are normalized variables and XI = XI (T) is the
C     stretching function.
C
C     For the cases of usual interest, B > 1 and z is purely imaginary,
C     leading to relationships involving sinh and tanh - hence the name
C     HTDIS2.  However, for B < 1, z is real and the analogous solution
C     involves sine and tangent.  Furthermore, for B ~ 1, both formula-
C     tions break down, so a third relation is used involving B - 1.
C
C     In this implementation, a Newton iteration is used to solve
C
C        SINH (DEL Y) / (DEL Y) = B   if  B > 1, or
C
C        SIN (DEL X) / (DEL X)  = B   if  B < 1.
C
C     Then for each XI = (I - 1) / (N - 1) an ANTISYMMETRIC function is
C     calculated:
C
C        U = 0.5 + 0.5 TANH [DEL Y * (XI - 0.5)] / TANH (0.5 * DEL Y) or
C
C        U = 0.5 + 0.5 TAN  [DEL X * (XI - 0.5)] / TAN  (0.5 * DEL X)
C
C     from which the unsymmetric general case is derived as
C
C        T = U / [A + (1 - A) * U]    where  A = SQRT (S0 / S1).
C
C     For B ~ 1. (actually | B - 1. | < 1.E-5), the approximation is
C
C        U = XI * [1. + 2. * (B - 1.) * (XI - 0.5) * (1. - XI)]
C
C     with no nonlinear equation to solve.  B = 1 if the distribution
C     is uniform, and this is handled as a special case if D1 = D2.
C
C  CLARIFICATION OF WHAT D1 AND D2 REALLY MEAN
C  -------------------------------------------
C
C     Confusion has arisen from the tendency to specify first and last
C     INTERVALS (D1 and D2), while the method calls for specifying SLOPES
C     at the stretching function's end-points.  In fact, one is related
C     to the reciprocal of the other.  We have, for the end-point slope,
C
C        Slope ~ [XI(2) - XI(1)] / [T(2) - T(1)]
C              = [(2 - 1) / (N - 1)] / dT(1)
C              = [1 / (N - 1)] * [1 / dT(1)]
C
C     Vinokur's use of the end slopes was intended to enable blending
C     of two distributions end-on-end smoothly.  However, there may be
C     situations where specifying and achieving precise increments is
C     preferred (for instance if it is desired to prepare a table
C     showing the behavior of the distribution in terms of largest and
C     smallest increments for varying N or for fixed N and varying D1, D2).
C
C     Therefore, this version has an option to perform an outer iteration
C     (via the secant method) in order to achieve specified increments
C     more precisely.  DSINPUT = .TRUE. invokes the outer iteration.
C
C     [This outer iteration really has to solve a nonlinear equation
C     in TWO variables - the D1, D2 arguments to lower-level routine
C     HTDIS3.  But independent parallel iterations with the secant
C     (non-derivative) method - one for adjusting D1 and one for
C     adjusting D2 - appear to be reliable in this peculiar circumstance.
C     Extreme cases can cause these iterations to converge before
C     the normal tolerances are met, but the results should be close
C     anyway.  IER = 3 in these cases.]
C
C  ARGUMENTS
C  ---------
C
C   NAME   DIM   TYPE I/O/S DESCRIPTION
C
C   DSINPUT -     L     I   .TRUE. means D1 and D2 are the desired first
C                                  and last increments;
C                           .FALSE. means they are the desired end-point
C                                  slopes of the normalized stretching
C                                  function.
C   XA      -     R     I   Desired X(1); passing X(1) here is safe.
C   XB      -     R     I   Desired X(N); passing X(N) here is safe.
C                           Note that XA = 0. and XB = 1. can provide one
C                           normalized distribution for multiple reuse,
C                           in some cases.  XA > XB is required.
C   D1,     -     R     I   See the DSINPUT description.  D1 and D2 should
C   D2                      be positive in either mode.  If they represent
C                           slopes, they refer to the end-slopes of the
C                           curve formed by plotting (I - 1) / (N - 1)
C                           (vertical axis) vs. (X (I) - XA) / (XB - XA)
C                           (horizontal axis).  This curve is independent
C                           of N for given slopes.  A larger slope means
C                           a higher density of points at the corresponding
C                           end, in normalized space.
C   N       -     I     I   Number of points; N > 3.
C   X       N     R     O   Reqd. distribn., with X(1) = XA, X(N) = XB.
C   LUNOUT  -     I     I   LUNOUT > 0 shows the iteration history;
C                           LUNOUT < 0 suppresses it (but any error
C                                      message goes to |LUNOUT|).
C   IER     -     I     O   IER = 0 means no problems;
C                                 1 means bad combination of D1 and D2,
C                                   such as D1 + D2 > XB - XA;
C                                 2 means the inner Newton iteration
C                                   failed - must have been bad inputs.
C                                   X (*) is not usable in this case;
C                                 3 means the outer iteration by the
C                                   secant method did not meet the
C                                   normal tolerances, presumably
C                                   because the case was extreme,
C                                   but things stopped changing so
C                                   it did the best it could.
C                                   X (*) should still be usable;
C                                 4 means the outer iteration did not
C                                   converge after 20 passes.  X (*)
C                                   is not usable in this case.
C  ERROR HANDLING:
C
C     See LUNOUT and IER.  Failure handling is left to the calling program,
C     which can proceed or reprompt for parameters and try again.
C
C  INTERNAL PROCEDURE:  HTDIS3 (the original HTDIS2 with additional arguments
C                       for performing an efficient outer iteration).
C
C  ENVIRONMENT:  FORTRAN 90
C
C  HISTORY:
C  c. 1980   M.Vinokur   Original analysis.
C  04/01/89  J.E.Melton  Initial implementation of HTDIS2 (bisection
C                        used for the inner iteration).
C  05/01/89  R.G.Langhi  Introduced error handling; patterned after
C                        one-sided routine HTDIS.
C  08/10/89  D.Saunders  Internal arithmetic is now in DOUBLE PRECISION
C                        to help application to Navier-Stokes grids.
C                        Results prove to be very similar, with the
C                        bisection just less likely to fail.  It seems
C                        great precision in DEL does not help:  X(2) and
C                        X(N-1) are still strangely imprecise except for
C                        the uniform case (where DEL in SINGLE and DOUBLE
C                        are quite different but give similar X(I)s).
C  11/30/90  D.Saunders  Safeguarded the B <= 1 cases, which have no soln.
C  04/23/91   "    "     Introduced an outer iteration for precise D1, D2:
C                        pushed the main algorithm down a level with added
C                        arguments for more efficient solution of the
C                        nonlinear equation.  6 steps and a fixed starting
C                        guess for DEL seem to suffice for likely cases.
C  04/28/91   "    "     Introduced a Newton inner iteration.  (Bisection
C                        hardly benefits from good starting guesses.)
C  06/04/91   "    "     Laborious explanation of the misconception that
C                        had arisen from the apparently imprecise results.
C                        Retained the more-precise option by setting the
C                        number of outer iterations to 1 or 6 according
C                        to the signs of the input D1, D2.
C  08/20/91   "    "     Incorporated the B < 1 and B ~ 1 cases.
C  08/28/91   "    "     The small-correction outer iteration was found to
C                        diverge on extreme cases with small N.  Replaced
C                        it with the secant method (up to 20 iterations
C                        for each of D1 and D2 in parallel).  IER = 3 or
C                        4 are new possibilities.
C  09/10/91   "    "     Introduced the DSINPUT argument to clarify the
C                        "deltas" or "slopes" options.
C  03/03/95   "    "     Encountered a near-singular case (D1 * D2 ~ Du**2)
C                        during iteration, which improperly terminated.
C                        EPS = 1.E-5 instead of 1.E-3 in HTDIS3 helps, and
C                        the outer iteration may now continue for CASE = 3.
C                        Improved the second estimates derived from the
C                        first call to HTDIS3.
C  03/03/95   "    "     HTDIS4 variation for 64-bit systems.
C  08/07/99   "    "     Fortran 90 upgrade: HTDIS3 is internal now, using
C                        just 2 arguments.  The other 10 become global.
C  01/23/05   "    "     One case out of nearly 28,000 hit the 30-iteration
C                        limit in HTDIS3.  It needed 36 iterations, so the
C                        limit is now 50.  Also, lowering EPS from 1.E-5 to
C                        1.E-6 gave identical results without encountering
C                        the case 3 situation that led to 36 iterations.
C
C  AUTHORS: Marcel Vinokur, NASA/Ames Research Center, Moffett Field, CA
C           and John Melton, Ron Langhi, David Saunders
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   IER, LUNOUT, N
      REAL
     >   D1, D2, XA, XB, X (N)
      LOGICAL
     >   DSINPUT

C     Local constants:

      INTEGER, PARAMETER ::
     >   MXITER = 20      ! Max. number of secant method iterations

      REAL, PARAMETER ::
     >   DLOW   = 1.E-2,  ! Smallest reduction in trial D1, D2 inputs to HTDIS3.
     >   ONE    = 1.E+0,
     >   THREE  = 3.E+0,
     >   TOL    = 1.E-4,  ! Relative tolerance on D1, D2
     >   ZERO   = 0.E+0

C     Local variables:

      INTEGER
     >   I

      REAL
     >   DEL0, DELSOLN, DMIN1, DMIN2, DX1 (3), DX2 (3), F1 (2), F2 (2),
     >   R, TOL1, TOL2, X1TARG, X2TARG

      LOGICAL
     >   CLOSE, TRIAL

C     Execution:
C     ----------

      TRIAL   = DSINPUT        ! Forces iteration to match Ds precisely
      DX1 (1) = D1
      DX2 (1) = D2

      IF (TRIAL) THEN
         TOL1   = DX1 (1) * TOL
         DMIN1  = DX1 (1) * DLOW
         X1TARG = XA + DX1 (1)      ! Target X (2)

         TOL2   = DX2 (1) * TOL
         DMIN2  = DX2 (1) * DLOW
         X2TARG = XB - DX2 (1)      ! Target X (N-1)

         IF (X1TARG > X2TARG) THEN  ! Bad D1 and D2
            IER = 1
            GO TO 99
         END IF
      END IF

      DEL0 = THREE  ! B > 1 gives DEL ~ 0.1 to 10 (increasingly non-uniform);
                    ! B < 1 gives DEL in the range [0+, PI-].
      CLOSE = .FALSE.     ! Used to detect "converged but above tolerance"

C     The first call to the basic Vinokur method is all we need if we are
C     matching slopes, else we need it to initialize the outer iteration:

C*****CALL HTDIS3 (DSINPUT, XA, XB, DX1 (1), DX2 (1), N, X, TRIAL, DEL0,
C    >             DELSOLN, LUNOUT, IER)

      CALL HTDIS3 (DX1 (1), DX2 (1)) ! Internal procedure now

      IF (IER /= 0) GO TO 99

      IF (.NOT. DSINPUT) GO TO 99    ! Slopes are matched without iteration.

      IF (DELSOLN == ZERO) GO TO 99  ! DEL = 0. signals uniform distribution
                                     ! which is specially handled.


C     Otherwise apply two secant iterations in parallel - one for D1 and
C     one for D2 - to solve what is really a nonlinear equation in TWO
C     variables.  Finish the first function evaluation:

      F1 (1) = X (2) - X1TARG      
      F2 (1) = X2TARG - X (N - 1)

C     Now for the second evaluation.  Attempt a better estimate.

      DX1 (2) = DX1 (1) ** 2 / (X (2) - XA)
      DX2 (2) = DX2 (1) ** 2 / (XB - X (N - 1))
      DEL0 = DELSOLN

      CALL HTDIS3 (DX1 (2), DX2 (2))

      IF (IER /= 0) GO TO 99

      F1 (2) = X (2) - X1TARG      
      F2 (2) = X2TARG - X (N - 1)

C     Here is the rest of the two secant iterations in parallel:

      DO I = 1, MXITER

         DEL0 = DELSOLN

         R = F1 (1) - F1 (2)
         IF (R /= ZERO) THEN
            R = F1 (1) / R
            DX1 (3) = R * DX1 (2) + (ONE - R) * DX1 (1)
            DX1 (3) = MAX (DX1 (3), DMIN1)
         ELSE
            DX1 (3) = DX1 (2)
         END IF

         R = F2 (1) - F2 (2)
         IF (R /= ZERO) THEN
            R = F2 (1) / R
            DX2 (3) = R * DX2 (2) + (ONE - R) * DX2 (1)
            DX2 (3) = MAX (DX2 (3), DMIN2)
         ELSE
            DX2 (3) = DX2 (2)
         END IF

         CALL HTDIS3 (DX1 (3), DX2 (3))

         IF (IER /= 0) GO TO 99

         IF (.NOT. TRIAL) THEN 
            IF (CLOSE) IER = 3
            GO TO 99
         END IF

C        Normal convergence test:
C        Reset TRIAL to get the full distribution on a final pass.

         F1 (1) = F1 (2)
         F1 (2) = X (2) - X1TARG      
         F2 (1) = F2 (2)
         F2 (2) = X2TARG - X (N - 1)

         IF (ABS (F1 (2)) < TOL1 .AND. ABS (F2 (2)) < TOL2) THEN
            TRIAL = .FALSE.
         ELSE IF (DX1 (3) == DX1 (2) .AND. DX2 (3) == DX2 (2)) THEN
            CLOSE = .TRUE.
            TRIAL = .FALSE.
         END IF

         DX1 (1) = DX1 (2)
         DX1 (2) = DX1 (3)
         DX2 (1) = DX2 (2)
         DX2 (2) = DX2 (3)

      END DO

      IER = 4  ! No convergence - X (*) is not usable.

   99 RETURN   ! From HTDIS4


C     Internal procedure for HTDIS4:
C     ------------------------------

      CONTAINS

!        -----------------------------------------------------------------------
!
         SUBROUTINE HTDIS3 (D1, D2)  ! The essence of the Vinokur method
!
!        -----------------------------------------------------------------------
!
!        ACRONYM:  Hyperbolic Tangent-type DIStribution, 2-sided (3rd version)
!                  -          -            ---                    -
!
!        HTDIS3 is the essence of the original HTDIS2.  It allows the specified
!        initial and final increments to be obtained with minimal work during
!        the intermediate iterations.
!
!        ADDITIONAL VARIABLES (See HTDIS4 for the others):
!
!        VAR      TYPE I/O/S DESCRIPTION
!
!        TRIAL     L     I   TRIAL = .TRUE. suppresses calculation of the
!                            interior X (*), since only X (2) and X (N-1)
!                            are being made use of at the higher level.
!        DEL0      R     I   Initial guess for the DELSOLN which solves
!                            the relevant nonlinear equation.  If B > 1,
!                            DEL is somewhere around 1. to 10. (larger
!                            for more nonuniform distributions); 1.E-10
!                            and 25. are used as bounds if necessary.
!                            If B < 1, the bounds are 1.E-10 and PI - eps.
!        DELSOLN   R     O   This solution is returned for reuse in
!                            a possible next call.
!
!        -----------------------------------------------------------------------

!        Arguments:

         REAL
     >      D1, D2

!        Local constants:

         INTEGER, PARAMETER ::
     >      MXITER = 50

         REAL, PARAMETER ::
     >      EPS = 1.E-6, HALF = 0.5E+0, ONE = 1.E+0,
     >      PIMINUS = 3.14159E+0, TOL = 1.E-12, TWO = 2.E+0,
     >      ZERO = 0.E+0

!        Local variables:

         INTEGER
     >      CASE, I, INC, ITER, LUNMSG

         REAL
     >      A, B, DEL, DELHI, DELLO, DELINV, DUNIFM, EXPDEL, FACTOR,
     >      FPRIME, FRACT, FUNDEL, RANGE, RFNM1, SLOPE0, SLOPE1, U, XI


!        Execution:
!        ----------

         IER = 0

!        Normalize the input control parameters:

         RANGE  = XB - XA
         RFNM1  = ONE / REAL (N - 1)
         DUNIFM = RANGE * RFNM1
         SLOPE0 = D1
         SLOPE1 = D2

         IF (DSINPUT) THEN
            SLOPE0 = DUNIFM / SLOPE0
            SLOPE1 = DUNIFM / SLOPE1
         END IF

         A = SQRT (SLOPE0 / SLOPE1)
         B = SQRT (SLOPE0 * SLOPE1)

!        B distinguishes the cases.

         IF (B == ONE) THEN

            IF (SLOPE0 == SLOPE1) THEN  ! D1 = D2 = uniform DX
               DEL = ZERO               ! Signal no need to iterate
               DO I = 1, N - 1
                  X (I) = XA + REAL (I-1) * DUNIFM
               END DO
               GO TO 40
            ELSE
               CASE = 3                 ! Singular case, but not uniform
               DEL = -ONE               ! Something other than zero
            END IF

         ELSE IF (B > ONE + EPS) THEN   ! Need to solve SINH (DEL) / DEL = B

            CASE = 1
            DELHI = 25.E+0
            DELLO = 1.E-10

         ELSE IF (B < ONE - EPS) THEN   ! Need to solve SIN (DEL) / DEL = B

            CASE = 2
            DELHI = PIMINUS
            DELLO = EPS

         ELSE                           ! Simpler approximation

            CASE = 3
            DEL = -ONE

         END IF


         IF (CASE /= 3) THEN

!           Newton iteration for the two main cases:

            IF (LUNOUT > 0) WRITE (LUNOUT, 1001) CASE

            DEL = DEL0

            DO ITER = 1, MXITER

               IF (CASE == 1) THEN  ! B > 1 + EPS

!                 Solve   SINH (DEL) / DEL = B  in the form
!
!                 F (DEL) = SINH (DEL) / (DEL * B) - 1 = 0  (better scaled?)
!                 The exponential form is used since EXP (DEL) and its
!                 reciprocal should be more efficient than SINH and COSH.

                  EXPDEL = EXP (DEL)
                  DELINV = ONE / DEL
                  FRACT  = HALF * DELINV / B
                  FUNDEL = FRACT * (EXPDEL - ONE / EXPDEL) - ONE
                  FPRIME = FRACT * ((ONE - DELINV) * EXPDEL +
     >                              (ONE + DELINV) / EXPDEL)

               ELSE  ! CASE = 2; B < 1 - EPS
         
!                 Solve   SIN (DEL) / DEL = B  in the form
!
!                 F (DEL) = SIN (DEL) / (DEL * B) - 1 = 0  (better scaled?)

                  FUNDEL = SIN (DEL) / (DEL * B) - ONE
                  FPRIME = (COS (DEL) / B - FUNDEL - ONE) / DEL

               END IF

               IF (LUNOUT > 0)
     >            WRITE (LUNOUT, 1002) ITER, DEL, FUNDEL, FPRIME, TOL

               IF (ABS (FUNDEL) < TOL) GO TO 20

               DEL = MIN (MAX (DEL - FUNDEL / FPRIME, DELLO), DELHI)

            END DO

            LUNMSG = ABS (LUNOUT)
            WRITE (LUNMSG, 1003) MXITER, DEL, FUNDEL, TOL
            IER = 2
            GO TO 99

   20       CONTINUE

         END IF

!        Calculate the antisymmetric function, from which the desired
!        stretching function is derived by another simple transformation.
!        The three cases are becoming cumbersome...

         IF (CASE == 1) THEN

            FACTOR = ONE / TANH (HALF * DEL)

         ELSE IF (CASE == 2) THEN

            FACTOR = ONE / TAN (HALF * DEL)

         ELSE

            FACTOR = TWO * (B - ONE)

         END IF

         X (1) = XA

         IF (TRIAL) THEN       ! Only the 1st & last increments are looked at
            INC = MAX (1, N - 3)
         ELSE
            INC = 1
         END IF

         DO I = 2, N - 1, INC

            XI = REAL (I-1) * RFNM1 - HALF  ! Subtracting 0.5 helps a bit

            IF (CASE == 1) THEN

               U = HALF * (ONE + TANH (DEL * XI) * FACTOR)

            ELSE IF (CASE == 2) THEN

               U = HALF * (ONE + TAN (DEL * XI) * FACTOR)

            ELSE

               U = (XI + HALF) * (ONE + FACTOR * XI * (HALF - XI))

            END IF

            FRACT = U / (A + (ONE - A) * U)
            X (I) = XA + FRACT * RANGE

         END DO

   40    X (N) = XB
         DELSOLN = DEL

         IF (LUNOUT > 0)
     >      WRITE (LUNOUT, 1004) B, X (2) - X (1), X (N) - X (N - 1)

   99    RETURN

 1001    FORMAT (/, ' ITER (Case = ', I1, ')    DEL', 7X, 'F(DEL)', 11X,
     >           'F''', 10X, 'TOL')
 1002    FORMAT (1X, I2, 1P, E20.12, 3E13.5)
 1003    FORMAT (/' *** HTDIS4:  Newton iteration failed ***',
     >           /' Number of iterations:', I4,
     >           /' DEL, F(DEL), TOL:', 1P, 3E14.6,
     >           /' No distribution computed.'/)
 1004    FORMAT (/, ' B:', 1P, E14.6, '  dX(1):', E14.6, '  dX(N-1):',
     >           E14.6)

         END SUBROUTINE HTDIS3

      END SUBROUTINE HTDIS4
