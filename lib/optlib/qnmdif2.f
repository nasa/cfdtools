C+----------------------------------------------------------------------
C
C                          CENDIF2 + QNMDIF2
C
C     UNCONSTRAINED OPTIMIZATION PACKAGE FOR EXPENSIVE FUNCTIONS
C         
C     This optimization package is intended for applications involving
C     CPU-intensive objective functions (such as those requiring CFD
C     flow solutions).
C
C     This version allows for use of two ways of estimating gradients:
C     the original way (by finite differencing) and an alternative way
C     (as by means of an "adjoint" solver).  The effect within QNMDIF2
C     of the alternative method is merely to by-pass any switching from
C     forward to central differences, because these would be identical.
C     (Use must still be made, in the calling program, of the fact that
C     QNMDIF2 performs finite differencing no matter where the gradient
C     information comes from.)
C
C     This package consists of two sub-packages, CENDIF2 and QNMDIF2.

C     CENDIF2 applies only to the original finite-differencing approach.
C     As a by-product of its estimation of good differencing intervals,
C     CENDIF2 determines the initial gradients and second partial
C     derivatives (diagonal elements only) of the objective function
C     with respect to the design variables.  CENDIF2 is adapted from
C     the original CENDIF package implemented by Robert Kennelly.
C
C     QNMDIF2 is the main optimization routine: a quasi-Newton algorithm
C     adapted from the QNMDIF of Gill, Murray, and Picken which itself
C     was adapted by Kennelly to handle objective functions calculable
C     to less than full machine precision in conjunction with the
C     implementation of CENDIF.
C
C     The term quasi-Newton means the Hessian matrix is not calculated
C     directly.  Rather, (factors of) a Hessian-like matrix are updated
C     by a rank two procedure (BFGS) at each optimization step.  Use of
C     the Cholesky factors allows safeguarding to ensure descent search
C     directions (unless the function is particularly poorly behaved).
C
C     The algorithm is unconstrained, meaning any desired constraints
C     must be imposed through quadratic penalties added to the objective
C     function.  Options in the original code include central differencing
C     or one-sided differencing as well as a random search both to check for
C     unimodality and to recover from difficult areas in the design space.
C
C     The two-part package has the following components:
C
C        CENDIF2:     Initial gradient and diagonal Hessian
C                     package (slightly adapted from CENDIF).
C           FDSTEP2:  Calculates optimal step sizes.
C           DIFFER2:  Calculates finite difference gradient.
C
C        QNMDIF2:     Main optimization routine (adapted from QNMDIF).
C           APROXG:   Calculates finite difference gradient.
C           C1MDCH:   Forms the Cholesky factored update to the
C                     Hessian.
C           LINSCH:   Performs a line search.
C           OUTPUT2:  Outputs optimization information.
C           UPCHOL:   Updates the Cholesky factored Hessian.
C
C-----------------------------------------------------------------------


C+----------------------------------------------------------------------
C
      SUBROUTINE CENDIF2 (N, X, FX, EPSOBJ, H, GRAD, DIAG, NFUNCT,
     >                    NALLOW, ZETA, LUNLOG, FUN, URH, NFDEL, IERROR)
C
C
C     Description:
C
C           CENDIF2 calculates central difference approximations to the
C        gradient and Hessian of an optimization objective function.  If
C        requested, improved finite difference step sizes are calculated.
C        CENDIF/FDSTEP/DIFFER were developed for use with QNMDIF (a
C        quasi-Newton optimizer), but the step sizes may be used by any
C        forward-difference method.
C
C
C     Arguments:
C
C        Name   Dimension   Type  I/O/S   Description
C
C        N                   I    I       Number of optimization variables.
C        X       N           R    I       Optimization variables.
C        FX                  R    I       Initial function value.
C        EPSOBJ              R    I       Estimate of absolute error in
C                                         objective.
C        H       N           R    I/O     Initial step sizes used for
C                                         differencing if >= 0.  If H(I)
C                                         is negative, indicates that an
C                                         improved value is to be estimated
C                                         by FDSTEP and returned.
C        GRAD    N           R      O     Central difference approximation to
C                                         gradient.
C        DIAG    N           R      O     Central difference approximation
C                                         to diagonal elements of Hessian.
C        NFUNCT              I      O     Number of function evaluations
C                                         (increment).
C        NALLOW              I    I       Max. # fn. evaluations allowed.
C        ZETA                R     I      The factor by which the the step
C                                         size is reduced in this subroutine.
C        LUNLOG              I     I      Logical unit printable output.
C        FUN                              Objective function routine:
C
C                                         CALL FUN (N, X, F, NCALL, FAIL).
C
C                                  NOTE:  NCALL is set internally by CENDIF2
C                                         to I = 1 : N, indicating the current
C                                         element of the gradient.
C
C        URH     N           R        S   Copy of input Hs.
C        NFDEL   N           I        S   Function count for each variable.
C        IERROR  N           I        S   Error flag for each variable.
C
C     Procedures:
C
C        DIFFER2   Computes first and second derivative approximations
C                  by central differences.
C        FDSTEP2   Produces estimates of optimal one sided finite
C                  difference step sizes.
C
C
C     Environment:  FORTRAN 77 with minor extensions.
C
C
C     Author: Robert Kennelly, Sterling Software/NASA Ames
C
C
C     History:
C
C        16 July 1982    RAK    Complete revision of earlier version by
C                               Trosin - now use FDSTEP and DIFFER.
C         2 Sep. 1982    RAK    Reordered parameters, redefined NFUNCT.
C        19 Jan. 1983    RAK    Updated NFDEL following DIFFER.
C        10 Mar. 1983    RAK    Repaired call to DIFFER (... DIAG(I) ...),
C                               added IF-THEN-ELSE blocks.
C         9 Aug. 1986    RAK    Added LUNLOG input parameter (had been
C                               in COMMON); reordered parameters again
C                               Use PARAMETER for local array dimensions,
C                               and test input N for consistency.  Increase
C                               MAXN from 30 to 40.  Add IMPLICIT NONE.
C           ?      J.J.Reuther  Added the NALLOW and ZETA arguments and
C                               common /OPT/.
C        Aug. 1994  D.Saunders  Effort to recover generality.
C        Feb. 1995   "    "     FUN now has NCALL as an argument;
C                               CENDIF2 & QNMDIF2 don't need NCALL argument.
C        Apr. 1995   "    "     All local arrays are now provided as
C                               arguments to avoid "MAXVAR" include file.
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   N, IERROR (N), NFDEL (N), LUNLOG, NFUNCT, NALLOW
      REAL
     >   DIAG (N), FX, EPSOBJ, GRAD (N), H (N), URH (N), X (N), ZETA
      EXTERNAL
     >   FUN

C     Local constants.

      REAL
     >   ZERO, HALF
      PARAMETER
     >  (ZERO = 0.0E+0,
     >   HALF = 0.5E+0)

C     Local variables.

      INTEGER
     >   I
      REAL
     >   CFPBAC, CFPFOR, CFPP, ERRBND, FPBAC, FPFOR, HI, URDIAG
      LOGICAL 
     >   FAIL

C     Procedures.

      EXTERNAL
     >   DIFFER2, FDSTEP2

      SAVE

C     Execution:


      WRITE (LUNLOG, 1000)

C     Loop over components of H and compute derivatives or improved
C     step sizes (as well as first and second derivatives).

      NFUNCT = 0
      DO 10 I = 1, N
         HI = H (I)
         URH (I) = HI
         IF (HI .GT. ZERO) THEN

C           Use the step size passed in for differencing.

            CALL DIFFER2 (I, N, X, EPSOBJ, HI, FX, FPFOR, CFPFOR, FPBAC,
     >                    CFPBAC, DIAG (I), CFPP, FAIL, LUNLOG, FUN)

            NFUNCT = NFUNCT + 2
            IERROR (I) = 0
            NFDEL (I) = 2
            GRAD (I) = HALF * (FPFOR + FPBAC)

         ELSE

C           Generate H from scratch, with GRAD and DIAG as by-products.

            NFDEL (I) = NFUNCT

            CALL FDSTEP2 (I, N, X, EPSOBJ, FX, GRAD (I), DIAG (I),
     >                    H (I), ERRBND, IERROR (I), NFUNCT, NALLOW,
     >                    LUNLOG, ZETA, FUN)

            NFDEL (I) = NFUNCT - NFDEL (I)

         END IF

         WRITE (LUNLOG, 1040) NFUNCT
         IF (NFUNCT .GE. NALLOW) THEN
            WRITE (LUNLOG, '(/, A)')
     >         ' CENDIF2: Max. # function evaluations reached.'
            STOP
         END IF

C        Safeguard Hessian approx. by bounding diagonals away from zero
C        to preserve positive definiteness.

         IF (DIAG (I) .LT. EPSOBJ) THEN
            URDIAG = DIAG (I)
            DIAG (I) = MAX (ABS (DIAG (I)), EPSOBJ)
            WRITE (LUNLOG, 1010) I, URDIAG, DIAG (I)
         END IF
   10 CONTINUE

C     Recap results.

      WRITE (LUNLOG, 1020)
      WRITE (LUNLOG, 1030) (I, X (I), GRAD (I), DIAG (I), URH (I),
     >   H (I), NFDEL (I), IERROR (I), I = 1, N)
      WRITE (LUNLOG, 1040) NFUNCT
      WRITE (LUNLOG, 1050)

      RETURN

C     Formats.
C     --------

 1000 FORMAT (/, ' CENDIF2: Begin gradient calculation.')
 1010 FORMAT (/, ' CENDIF2: WARNING - DIAG(', I3, ') changed from ',
     >   E12.6, ' to ', E12.6)
 1020 FORMAT ('1CENDIF2:', //, '     I    V', 20X, 'Gradient', 13X,
     >   'Diag. of Hessian     Old H', 16X, 'New H', 14X,
     >   '# Fns.  IERROR')
 1030 FORMAT (4X, I2, 3X, E18.12, 3X, E18.12, 3X, E18.12, 3X, E18.12,
     >   3X, E18.12, 4X, I2, 6X, I1)
 1040 FORMAT (/, ' CENDIF2: Total # function evaluations so far:', I6)
 1050 FORMAT ('1')

      END


C+----------------------------------------------------------------------
C
      SUBROUTINE FDSTEP2 (I, N, X, EPSOBJ, FX, FP, FPP, HOPT, ERRBND,
     >                    IERROR, NFUNCT, NALLOW, LUNLOG, ZETA, FUN)
C
C
C     Description:
C
C           This routine is intended for automatic estimation of optimal
C        finite difference intervals for numerical optimization of functions
C        whose derivatives cannot be calculated analytically.  Each call to
C        FDSTEP2 will handle one component of the step size vector.
C
C
C     Arguments:
C
C        Name    Dimension   Type  I/O/S  Description
C        I                    I    I      Index of the component of X to be
C                                         varied.
C        N                    I    I      Dimension of X.
C        X          N         R    I      Vector of optimization variables.
C        EPSOBJ               R    I      Smallest significant absolute change
C                                         in the objective function FUN.
C        FX                   R    I      Initial value of objective function.
C        FP                   R      O    Best central difference approximation
C                                         to first partial derivative of
C                                         objective with respect to the I-th
C                                         component of X.
C        FPP                  R      O    Best central difference approximation
C                                         to second partial derivative of
C                                         objective with respect to the I-th
C                                         component of X.
C        HOPT                 R      O    Estimate of optimal finite difference
C                                         step size for I-th component of X.
C        ERRBND               R      O    Estimated bound on (absolute) error
C                                         in first derivatives calculated
C                                         using one sided differencing with
C                                         step size HOPT.
C        IERROR               I      O    Error status upon termination. Zero
C                                         means OK; for other values see
C                                         error handling section below.
C        NFUNCT               I      O    Number of calls to FUN (increment).
C        NALLOW               I    I      Max. # fn. evals. allowed.
C        LUNLOG               I    I      Logical unit for output.
C        ZETA                 R    I      The factor by which the step size
C                                         is reduced.
C        FUN                              Routine for calculating objective
C                                         function (EXTERNAL).
C
C     Procedures:
C
C        DIFFER2    Utility for calculating derivatives by finite differences,
C                   along with cancellation error estimates.
C        FUN        User-supplied objective function subroutine.
C
C
C     Notes:
C
C        (1)  The goal is to find the smallest step size H(KPHI) which yields
C             acceptable relative cancellation errors in the forward and
C             backward difference first derivatives and the central difference
C             second derivative of FUN.  The second derivative computed with
C             H(KPHI) is then used to obtain HOPT, an interval which should
C             result in good forward difference approximations to the first
C             derivative.  The "optimality" of the interval chosen is based
C             on several assumptions - consult reference (1) for details.
C
C        (2)  The algorithm also returns estimates FP, the first derivative,
C             and FPP, the second derivative, which are both computed by
C             central differences.
C
C        (3)  The step size and derivatives calculated each iteration are
C             saved in arrays for possible later use to save function
C             evaluations.  These are relevant to the step associated with
C             optimization variable I only, and should not be confused with
C             the step size and gradient arrays in the calling routine.
C
C        (4)  Constant TINY, set below, is a machine dependent quantity
C             intended to provide some protection against division by zero.
C             The arithmetic expressions involved could still overflow for
C             badly scaled problems in which EPSOBJ is much larger than one.
C
C        (5)  FDSTEP and DIFFER are based on Algorithm FD in Reference (1).
C             The parenthesized labels in the code (e.g. FD1) refer to steps
C             of the algorithm as originally published.
C
C
C     Bibliography:
C
C        (1)  Gill, P., Murray, W., and Wright, M.  Practical Optimization,
C                pp. 341-344.  London: Academic Press, 1981.
C
C        (2)  Gill, P., Murray, W., Saunders, M., and Wright, M.  "A Procedure
C                for Computing Forward Difference Intervals for Numerical
C                Optimization."  Tech. Rep. SOL 81-25, Dept. of Operations
C                Research, Stanford University, Dec. 1981.
C
C
C     Environment:  FORTRAN 77
C
C
C     Author: Robert Kennelly, NASA Ames Research Center
C
C
C     History:
C
C        16 June 1982   RAK   Original coding ("literal" transcription).
C        23 June 1982   RAK   Revised and extensively restructured.
C        13 July 1982   RAK   Section FD3 revised (small CFPP now OK).
C         2 Sep. 1982   RAK   Redefined function counter NFUNCT.
C        24 Sep. 1982   RAK   Some minor changes adapted from Ref. (2).
C        11 Aug. 1986   RAK   Add IMPLICIT NONE.
C            ?       J.J.Reuther  Added ZETA argument and /OPT/.
C        17 Aug. 1994   DAS   Some clean-up.  No need for an INCLUDE.
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   I, IERROR, LUNLOG, N, NFUNCT, NALLOW
      REAL
     >   EPSOBJ, ERRBND, FP, FPP, FX, HOPT, X (N), ZETA
      EXTERNAL
     >   FUN

C     Local constants.

      INTEGER
     >   KBOUND
      PARAMETER
     >  (KBOUND = 6)

      REAL
     >   ZERO, ONE, TWO, HALF, BNDFP, BNDLO, BNDUP, ETA, OMEGA, RHO,
     >   TINY, BNTOP
      PARAMETER
     >  (ZERO  = 0.0E+0,
     >   ONE   = 1.0E+0,
     >   TWO   = 2.0E+0,
     >   HALF  = ONE / TWO,
     >   BNDFP = 1.0E-1,
     >   BNDLO = 1.0E-3,
     >   BNDUP = 1.0E-1,
     >   BNTOP = 1.0, 
     >   ETA   = ONE,
     >   OMEGA = ONE,
     >   TINY  = 1.0E-20)

C     Local variables.

      INTEGER
     >   KOUNT, KPHI, KSAVE
      REAL
     >   CFPP, CFPBAC, CFPFOR, H (0:KBOUND), HBAR, FPFOR (0:KBOUND),
     >   FPBAC (0:KBOUND), FPPCEN (0:KBOUND)
      LOGICAL
     >   DONE, FAIL, LFAIL

C     Procedures.

      EXTERNAL
     >   DIFFER2

C     Storage.

      SAVE

C     Statement functions.

      REAL
     >   A, B
      LOGICAL
     >   INSIDE, MAXLE, RLARGE

      INSIDE (A)   = (A .GE. BNDLO) .AND. (A .LE. BNDUP)
      MAXLE (A, B) = MAX (A, B) .LE. BNDFP
      RLARGE (A)   = A .GT. BNTOP

C     Execution.
C     ----------

C     Initialization (FD1).

      KOUNT = 0
      IERROR = 0
      DONE = .FALSE.
      RHO = ZETA

      HBAR = TWO * (ETA + ABS (X (I))) *
     >   SQRT (EPSOBJ / (OMEGA + ABS (FX)))
      H (KOUNT) = HBAR * 10.0

      CALL DIFFER2 (I, N, X, EPSOBJ, H (KOUNT), FX, FPFOR (KOUNT),
     >              CFPFOR, FPBAC (KOUNT), CFPBAC, FPPCEN (KOUNT), CFPP,
     >              FAIL, LUNLOG, FUN)

      NFUNCT = NFUNCT + 2
      IF (NFUNCT .GE. NALLOW) THEN
         WRITE (LUNLOG, '(/, A)') ' FDSTEP2:  Fn. eval. limit reached.'
         STOP
      END IF

C     Accept, increase, or decrease H?  (FD2).

      IF ((MAXLE (CFPFOR, CFPBAC) .AND. .NOT.RLARGE (CFPP)) .OR. FAIL)
     >   THEN
         IF (INSIDE (CFPP) .AND. .NOT. FAIL) THEN

C           Accept H and fall through to estimate of optimal H.

            KPHI = KOUNT
         ELSE

C           Decrease H either to reduce truncation error (FD4) or simply
CJJR        to allow the flow solver to converge.
C           The cancellation errors in the first derivatives are OK, while
C           that for the second derivative is smaller than necessary.
C           Try to reduce H without letting the errors get too big.
C
CJJR              NOTE: If the program has entered this IF THEN block  
CJJR              due to the failure of the flow solver then it must be 
CJJR              noted that the primary assumption is that the initial 
CJJR              solution is a highly convergent case.  Therefore the 
CJJR              smaller the deviation from the initial solution the 
CJJR              more likely that convergence will be obtained, hence 
CJJR              the program will always try to decrease the step size 
CJJR              to secure convergence.  This enforced flow solver 
CJJR              convergence criterion overrides the possible fact that 
CJJR              a particular design variable may be suffering from 
CJJR              cancellation error due to too small a step size, and 
CJJR              thus may be ill-scaled for QNMDIF.

            IF (FAIL) THEN
               WRITE (LUNLOG, '(/, A, A, I6, /, A, E20.10)')
     >            ' FDSTEP2: Function failure for initial H.',
     >            ' Design variable: ', I,
     >            ' Attempting to decrease H from ', H (KOUNT)
            END IF

   10       CONTINUE
               KOUNT = KOUNT + 1
               H (KOUNT) = H (KOUNT - 1) / RHO
               LFAIL = FAIL

               CALL DIFFER2 (I, N, X, EPSOBJ, H (KOUNT), FX,
     >                       FPFOR (KOUNT), CFPFOR, FPBAC (KOUNT),
     >                       CFPBAC, FPPCEN (KOUNT), CFPP, FAIL, LUNLOG,
     >                       FUN)

               IF (FAIL) THEN
                  WRITE (LUNLOG, '(/, A, I6, A, E20.10)')
     >               ' FDSTEP2: Function failure for variable ', I,
     >               '  for current reduced H', H (KOUNT)
               END IF

               NFUNCT = NFUNCT + 2
               IF (NFUNCT .GE. NALLOW) THEN
                  WRITE (LUNLOG, '(/, A)')
     >               ' FDSTEP2:  Fn. eval. limit reached.'
                  STOP
               END IF

               IF (.NOT. FAIL) THEN
                  IF ((.NOT. MAXLE (CFPFOR, CFPBAC) .OR. RLARGE (CFPP))
     >               .AND. .NOT. LFAIL) THEN

C                    We've gone too far - back up one iteration and quit.

                     DONE = .TRUE.
                     KPHI = KOUNT - 1
                  ELSE IF (INSIDE (CFPP)) THEN

C                    The current step size H is acceptable.

                     DONE = .TRUE.
                     KPHI = KOUNT
                  ELSE

C                    Second derivative cancellation error is still smaller
C                    than necessary, or the objective fn. has still failed to
CJJR                 converge, but check iteration counter before continuing.
C
                     IF (KOUNT .GE. KBOUND) THEN

C                       The iteration limit has been reached.  The error
C                       flag is set as a warning that the step size may not
C                       be as small as possible, or objective fn. convergence
CJJR                    may not have been reached.  Note: this can also
C                       indicate trouble - a slope discontinuity can also
C                       give this error, so the objective function should be
C                       double-checked in this case.

                        IERROR = 1
                        KPHI = KOUNT
                     END IF
                  END IF
               ELSE
                  IF (.NOT. LFAIL) THEN
                     WRITE (LUNLOG, '(/, A, I6, A, /, A)')
     >             ' FDSTEP2: Variable ', I, '  is poorly chosen.',
     >             ' Program will use previous step size.'
                     DONE = .TRUE.
                     KPHI = KOUNT - 1
                  ELSE
                     WRITE (LUNLOG, '(A, A)')
     >                  ' FDSTEP2: After iteration check,',
     >                  ' step size will be further reduced.'
                     IF (KOUNT .GE. KBOUND) THEN
                        IERROR = 1
                        KPHI = KOUNT
                     END IF
                  END IF
               END IF
               IF (.NOT.DONE .AND. IERROR .EQ. 0) GO TO 10

         END IF
      ELSE

C        Increase H to reduce cancellation error bounds in first or second
C        derivatives (FD3).

         KSAVE = -1
   20    CONTINUE

C           Keep track of index of smallest H (if any) with sufficiently
C           small cancellation errors in the one sided first derivatives.

            IF (MAXLE (CFPFOR, CFPBAC) .AND. KSAVE .LT. 0) KSAVE = KOUNT
            KOUNT = KOUNT + 1
            H (KOUNT) = H (KOUNT - 1) * RHO

            CALL DIFFER2 (I, N, X, EPSOBJ, H (KOUNT), FX, FPFOR (KOUNT),
     >                    CFPFOR, FPBAC (KOUNT), CFPBAC, FPPCEN (KOUNT),
     >                    CFPP, FAIL, LUNLOG, FUN)

            NFUNCT = NFUNCT + 2
            IF (NFUNCT .GE. NALLOW) THEN
               WRITE (LUNLOG, '(/, A)')
     >            ' FDSTEP2:  Fn. eval. limit reached.'
               STOP
            END IF

            IF (FAIL) THEN
              WRITE (LUNLOG, '(/, A, I6, /, A, E20.10)')
     >           ' FDSTEP2: Function failure for design variable ', I,
     >           ' after increasing step size to ', H (KOUNT),
     >           ' Reverting to previous step size.'

              DONE = .TRUE.
              KPHI = KOUNT - 1

            ELSE IF (MAXLE(CFPFOR,CFPBAC) .AND. .NOT.RLARGE (CFPP)) THEN

C              Current H is acceptable.

               DONE = .TRUE.
               KPHI = KOUNT
            ELSE IF (KOUNT .GE. KBOUND) THEN

C              No satisfactory interval has been found (FD6).

               IF (.NOT. MAXLE (CFPFOR, CFPBAC)) THEN

C                 Cancellation errors in first derivatives are still too
C                 big - the function is nearly constant.  Use the step
C                 appropriate for the "well scaled" case.

                  IERROR = 2
                  HOPT = HBAR
                  FP = ZERO
                  FPP = ZERO
               ELSE

C                 Second derivative cancellation error is excessive - the
C                 first derivative appears to be nearly constant.  Use the
C                 smallest H with acceptable cancellation error in first
C                 derivatives.

                  IERROR = 3
                  HOPT = H (KSAVE)
                  FP = FPFOR (KSAVE)
                  FPP = ZERO
               END IF
            END IF
            IF (.NOT.DONE .AND. IERROR .EQ. 0) GO TO 20

      END IF
      IF (IERROR .LE. 1) THEN

C        Set output values of derivatives and estimate of optimal H (FD5).

         FP = HALF * (FPFOR (KPHI) + FPBAC (KPHI))
         FPP = FPPCEN (KPHI)
         HOPT = MAX (TWO * SQRT (EPSOBJ / MAX (ABS (FPP), TINY)), TINY)
      END IF

C     Calculate estimate of error bound on one sided derivatives based
C     on finite difference interval HOPT.  If interval should be OK but
C     bound on relative error is large, set warning flag.

      ERRBND = HALF * HOPT * ABS (FPP) + TWO * EPSOBJ / HOPT
      IF ((IERROR .EQ. 0) .AND. (ERRBND .GT. HALF * ABS (FP)))
     >   IERROR = 4

      RETURN
      END


C+----------------------------------------------------------------------
C
      SUBROUTINE DIFFER2 (I, N, X, EPSOBJ, H, FX, FPFOR, CFPFOR, FPBAC,
     >                    CFPBAC, FPPCEN, CFPP, FAIL, LUNLOG, FUN)
C
C
C     Description:
C
C           Utility for calculating derivatives and cancellation error
C        estimates by finite differences.  Called repeatedly by subroutine
C        FDSTEP during estimation of optimal step sizes to be passed to a
C        nonlinear optimization routine.  See CENDIF/FDSTEP headers for
C        more details.
C
C
C     Arguments:
C
C        Name    Dimension   Type  I/O/S  Description
C        I                    I    I      Index of the component of X to be
C                                         varied.
C        N                    I    I      Dimension of X in calling routine.
C        X          N         R    I      Vector of optimization variables.
C        EPSOBJ               R    I      Smallest significant absolute change
C                                         in the objective function FUN.
C        H                    R    I      Finite difference step size for I-th
C                                         component of X.  Assumed > 0.
C        FX                   R    I      Initial value of objective function.
C        FPFOR                R      O    Forward difference approximation
C                                         to first partial derivative of
C                                         objective with respect to the I-th
C                                         component of X.
C        CFPFOR               R      O    Relative cancellation error in FPFOR.
C        FPBAC                R      O    Backward difference approximation
C                                         to first partial derivative of
C                                         objective with respect to the I-th
C                                         component of X.
C        CFPBAC               R      O    Relative cancellation error in FPBAC.
C        FPPCEN               R      O    Central difference approximation
C                                         to second partial derivative of
C                                         objective with respect to the I-th
C                                         component of X.
C        CFPP                 R      O    Relative cancellation error in
C                                         FPPCEN.
C        FAIL                 L      O    .TRUE. means objective fn. failed.
C        LUNLOG               I      I    Logical unit for output
C        FUN                              Routine for calculating objective
C                                         function (EXTERNAL).
C
C
C     Notes:
C
C        (1)  Constant TINY, set below, is a machine dependent quantity
C             intended to provide some protection against division by zero.
C             The arithmetic expressions involved could still overflow for
C             badly scaled problems in which EPSOBJ is much larger than one.
C
C
C     Environment:  FORTRAN 77
C
C
C     Author:  Robert Kennelly, NASA Ames
C
C
C     History:
C
C        23 June 1982    RAK    Original coding.
C         2 Sep. 1982    RAK    Deleted local function counter.
C        11 Aug. 1986    RAK    Add IMPLICIT NONE.
C            ?     J.J.Reuther  Significant changes as for CENDIF2, FDSTEP2
C           Aug. 1994    DAS    No need for a common block; some clean-up.
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   I, N, LUNLOG
      REAL
     >   CFPBAC, CFPFOR, CFPP, EPSOBJ, FPBAC, FPFOR, FPPCEN, FX, H,
     >   X (N)
      LOGICAL
     >   FAIL
      EXTERNAL
     >   FUN

C     Local constants.

      REAL
     >   TWO, FOUR, TINY
      PARAMETER
     >  (TWO  = 2.0E+0,
     >   FOUR = 4.0E+0,
     >   TINY = 1.E-20)

C     Local variables.

      INTEGER
     >   MFAIL
      REAL
     >   FXMH, FXPH, XINIT

C     Storage.

      SAVE


C     Execution.
C     ----------

      XINIT = X (I)

C     Compute first derivative by forward difference...

      X (I) = XINIT + H

      CALL FUN (N, X, FXPH, I, FAIL)

      IF (FAIL) THEN
         WRITE (LUNLOG, '(/, A)')
     >      ' DIFFER2: Function failure - continuing.'
         FXPH = FX
         MFAIL = 1
      END IF
      FPFOR = (FXPH - FX) / H
      CFPFOR = TWO * EPSOBJ / MAX (H * ABS (FPFOR), TINY)

C     ... and by backward difference.

      X (I) = XINIT - H

      CALL FUN (N, X, FXMH, I, FAIL)

      IF (FAIL) THEN
         WRITE (LUNLOG, '(/, A)')
     >      ' DIFFER2: Function failure; continuing.'
         FXMH = FX
      END IF
      FPBAC = (FX - FXMH) / H
      CFPBAC = TWO * EPSOBJ / MAX (H * ABS (FPBAC), TINY)

C     Calculate second derivative by central differences.

      FPPCEN = (FXPH - TWO * FX + FXMH) / MAX (H ** 2, TINY)
      CFPP = FOUR * EPSOBJ / MAX (H ** 2 * ABS (FPPCEN), TINY)

      WRITE (LUNLOG, '(/, A)')
     >   ' DIFFER2: Gradient is calculated by central differences.'
      WRITE (LUNLOG, '(1X, A, E17.8)')
     >   'Both current step-sizes   :', H,
     >   'Current forward  gradient :', FPFOR,
     >   'Current backward gradient :', FPBAC,
     >   'Current central 2nd deriv.:', FPPCEN

C     Restore X (I) to original value.

      X (I) = XINIT
      IF (MFAIL .EQ. 1) FAIL = .TRUE.

      RETURN
      END


C+----------------------------------------------------------------------
C
      SUBROUTINE QNMDIF2 (N, NLDIM, NFTOTL, KITER, NTYPE, LUNLOG, X, F,
     >                    FM, G, H, EL, D, ETA, TOL, STEPMX, EPSMCH,
     >                    EPSOBJ, UNITL, LOCAL, CONV, CONTIN, RESCUE,
     >                    PRINT, FUN, USER, OLDG, P, FFWD, SH, Y, Z,
     >                    PWORK, LFAIL, GTP, GNM, ZETA, NALLOW, IGRAD,
     >                    ISTART)

C     This is James Reuther's version of QNMDIF, q.v., as overhauled
C     for reusability by David Saunders.
C
C     Arguments retained from the original QNMDIF are described below
C     with minor streamlining and updating as necessary.
C
C      SUBROUTINE QNMDIF (N, NLDIM, NFTOTL, NITER, NTYPE, LUNLOG,
C     >   X, F, FM, G, H, EL, D, ETA, TOL, STEPMX, EPSMCH, EPSOBJ,
C     >   UNITL, LOCAL, CONV, CONTIN, RESCUE, PRINT, FUN, USER)
C
C     N - The number of independent variables, an integer constant.
CDAS     N >= 1 is not limited by any local storage.  For very large N,
CDAS     NLDIM (next) may be the limiting factor.
C
C     NLDIM - The value N * (N-1)/ 2, which is the number of elements
C        in the matrix EL, stored as a vector.  Note, however, that
C        NLDIM must always be at least 1, so that if QNMDIF is called
C        with N = 1, NLDIM should be set to 1.
C
C     NFTOTL - An integer, which is incremented by QNMDIF with the total
C        number of calls to the objective function executed during the
CDAS     location of the minimum.  It should be initialized before calling
CDAS     QNMDIF (possibly allowing for function calls from CENDIF).
C
C     NITER - An integer, which is returned by QNMDIF with the number of
CDAS     iterations (design steps) executed.  It should be input with
CDAS     the desired number of steps, possibly small for expensive functions.
C
CRAK  NTYPE - An integer flag set by QNMDIF to indicate whether the
CRAK     gradient evaluated before returning used forward (NTYPE = 0)
CRAK     or central differences (NTYPE = 1).  On input, must be set to
CRAK     correspond to the user-supplied gradient if UNITL = .FALSE.
CRAK     or to the desired initial gradient type if starting from
CRAK     scratch (UNITL = .TRUE.).
C
CDAS  LUNLOG - Logical unit number for printed output.
C
C     X - The array of independent variables, a vector of length N.
C        On entry, it should be set to the initial estimate of the
C        location of the minimum.  On output it will contain the best
C        estimate found of the minimum.
C
CRAK  F - A real scalar.  On entry, it must be set to the initial
CRAK     function value.  On exit it will be the lowest value found
CRAK     of the objective function.
CRAK
CRAK  FM - A real scalar.  On entry, it should be set to a lower
CRAK     bound on the value of the function at the minimum, or to
CRAK     a large negative value if no estimate is available.
CRAK
CRAK  G - A real array of length N, containing the gradient, which
CRAK     must be initialized if UNITL = .FALSE., and otherwise not.
CRAK     On exit, equals the last gradient evaluated, which may be
CRAK     useful.
C
C     H - A difference array, of length N, where H (I) is the interval
C        used to compute the partial derivative of F with respect to
C        X (I).  Great care should be used in choosing the H array,
CDAS     which is the reason for the existence of the CENDIF package.
C
C     EL - A vector of length NLDIM (at least 1), used to hold the
C        lower triangular factor of the approximate Hessian.  If the
C        user wishes to provide an initial estimate of the Hessian, EL
C        should be filled with the elements of the lower triangle of the
C        Cholesky factorization EL * D * ELT, where ELT is EL tranpose,
C        stored by rows.  Otherwise, set UNITL to .TRUE., and QNMDIF will
C        use the identity matrix as the initial approximation (i.e., the
C        first iteration will be steepest descent).  On exit, EL contains
C        the lower triangle of the factors of the final approximate Hessian.
C
C     D - A vector of length N, used to hold the diagonal matrix of the
C        factors of the approximate Hessian.  If the user wishes to
C        provide an initial estimate of the Hessian, D should be set in
C        the same way as EL (described above).  Otherwise, when UNITL is
C        set to .TRUE., the identity is used for the first iteration, and
C        the user need not set EL and D. On exit, D contains the diagonal
C        factor of the final approximate Hessian.
C
C     ETA - A scalar constant, set by the user to a number strictly
C        between 0 and 1.  It controls the accuracy of the search for the
C        minimum of the function along the current direction of search,
C        in the following way:  if the exact minimum along the line were
C        found at each step, the projected gradient (gradient of the
C        function along the line) would be zero.  However, it has been
C        found that finding the exact minimum at each step is not
C        computationally efficient, in that it requires many extra
C        function evaluations, and does not in general assist the speed
C        of convergence dramatically.  The parameter ETA is multiplied by
C        the projected gradient at the start of each minimization along a
C        line.  When the projected gradient at a point is less than ETA
C        times the initial value, the linear search is terminated.  Thus,
C        a small value of ETA means that the projected gradient must be
C        very small, i.e., close to the minimum along the line.  A value
C        of ETA close to 1 means that the projected gradient need
C        decrease only slightly from the starting point, so that the
C        resulting point is only a rough approximation to the minimum
C        along the line.  A guide to its value can be found in the NPL
C        report mentioned above, but the user may wish to try different
C        values on experimental runs to find the best value for his problem.
C        A value of 0.2E+0 is generally acceptable for QNMDIF.
CDAS     For very expensive functions being finite-differenced, 0.1E+0
CDAS     would normally be more appropriate in view of the cost of the
CDAS     the gradient estimate needed to set up the line search.
C
C     TOL - A real positive number, which specifies the accuracy to
C        which the minimum is to be found.  QNMDIF terminates when the
C        norm of the approximate gradient is less than TOL ** (1/ 3) and
C        the distance moved during the last iteration is less than a
C        number close to TOL.  TOL should be chosen as approximately
C        the Euclidean distance which can be allowed between the point
C        returned by QNMDIF and the actual minimum.  Its value should
C        be chosen with care, keeping in mind that the error in the
C        approximate gradient is of order MAXH (I) ** 2 at best (when
C        central differences are used), where MAXH is the maximum
C        component of H in absolute value, that the function may be
C        known only to a given number of figures, and that the
C        parameters may be known only to a given number of figures
C        (for example, when the problems arises from a physical
C        situation, or a problem involving measured data).  A lower
C        bound for an acceptable value of TOL for most problems
C        is SQRT (EPSOBJ), but a larger value may be required.
C
C     STEPMX - The maximum allowed step from any point to the minimum
C        along a line from that point.  This argument should be used
C        when the user is confident that the X vector should not
C        move more than the distance STEPMX from any particular
C        location.  If the linear search fails to find a minimum closer
C        to the current point than STEPMX, then the distance moved is
C        STEPMX.  It should be set to an upper bound on the distance
C        between the starting point and the minimum.  If the user does
c        not have a reasonable guess for this value, it can simply be
C        set to a very large number (say 1.0E+10), and then will have
C        essentially no effect on the linear search.
C
C     EPSMCH -  The machine precision of the computer to be used, i.e.,
C        the smallest value such that, in floating point arithmetic,
C        1 + EPSMCH > 1.  On the IBM 360, EPSMCH is 16.0 ** (-13).
CDAS     An available EPSLN utility could be used by the calling program.
C
CRAK  EPSOBJ - Must be input by the calling program as the smallest
CRAK     significant (absolute) change in the objective function.  QNMDIF
CRAK     resets EPSOBJ to EPSMCH if the input value is less than EPSMCH.
C
C     UNITL - A logical parameter, set by the user.  In general, the
C        user will not know an approximation to the Hessian matrix
C        at the starting point, and in this event UNITL should be set
C        to .TRUE., so that the identity matrix is used as the starting
C        approximation to the Hessian.  If, however, the user does have
C        an idea of the approximate Hessian, UNITL should be set to
C        .FALSE., and the EL and D arrays should be set to the Cholesky
C        factors of the initial approximation (see comments under EL
C        and D, above).
C
CRAK     MODIFICATION NOTE:  UNITL = .FALSE now implies that all of the
CRAK     information needed to begin an iteration (F, X, G, EL, and D)
CRAK     has been supplied (not merely EL and D).
C
C     LOCAL - A logical parameter, to be set to .TRUE. or .FALSE. by the
C        user.  Gill and Murray have developed a sophisticated (?) "local
C        search" which follows the normal minimization, and which checks
C        through a random search whether a point can be found lower than
C        the estimated minimum.  The purpose of this search is to avoid
C        convergence to a saddle point, where the gradient is zero, but
C        which is not a minimum.  This search requires a fair number of
C        extra function evaluations, but will essentially guarantee
C        convergence to a true minimum.  If the user feels that there is
C        no danger of convergence to a saddle point, and in particular,
C        if the function is extremely expensive to evaluate, then LOCAL
C        should be set to .FALSE., and the local search will be skipped.
C
C     CONV - A logical flag set by QNMDIF to .TRUE. if there has been
C        satisfactory convergence to a minimum, and to .FALSE. if not.
C        If LOCAL is .TRUE., then CONV = .TRUE. means that the local
C        search failed to find a significantly lower point.  If LOCAL
C        was set to .FALSE., then CONV = .TRUE. means that the gradient
C        at the final point was essentially zero, and that the steps to
C        the final point had converged.  Note that if LOCAL is set to
C        .FALSE., CONV = .TRUE. can happen if there is convergence to
C        a saddle point, since there is no way to check the second order
C        conditions to distinguish a minimum from a saddle point.
C
CRAK  CONTIN - A logical input flag which must be .TRUE. if QNMDIF is
CRAK     to update the gradient and Hessian before returning due to
CRAK     iteration count.  Efficiency is not compromised if the final
CRAK     values are saved and later used by the calling routine for
CRAK     restarting.
CRAK
CRAK  RESCUE - A logical input flag which may be set .TRUE. if the user
CRAK     wishes to begin the run with the local search procedure, which
CRAK     counts as two iterations.  Normally, use RESCUE = .FALSE. and set
CRAK     LOCAL = .TRUE. if some assurance against premature convergence
CRAK     is desired.  The RESCUE option may be useful if the objective
CRAK     function routine has blundered into a bad region.
C
C     PRINT - A logical flag, set by the user to control printout.
C        for PRINT = .FALSE., most printing is suppressed.
C
CRAK  USER - A subroutine provided by the user, which must be declared
CRAK     as EXTERNAL in the calling subprogram.  Its call must be of
CRAK     the form
CRAK
CDAS     CALL USER (IENTRY, N, NLDIM, X, F, G, H, EL, D, OLDG, P, GTP,
CDAS                GNORM, NFTOTL, NALLOW, NCALL, NITER, NTYPE, CONTIN,
CDAS                PRINT, LUNLOG)
CRAK
CRAK     where the variables are the current values of those described
CRAK     above.  QNMDIF calls USER at the end of each successful line
CRAK     search with IENTRY = 1, at which point X, F, H, NFTOTL, and
CRAK     NITER are current, but G, EL, and D have not yet been updated.
CRAK     This permits certain housekeeping chores to be performed before
CRAK     the calculations for the next iteration are begun.  QNMDIF sets
CRAK     IENTRY = 2 and calls USER again following calculation of gradient
CRAK     and Hessian updates, as well as after the local search if
CRAK     requested.  Disk files for saving restart information should be
CRAK     written when IENTRY = 2, with an initial REWIND command.
CRAK     Different applications may dictate other choices for parameter
CRAK     list and calling points.  Note that either or both of these
CRAK     options may be simply turned off by providing subroutine USER
CRAK     with RETURN statements for the unneeded cases.
C
C
C     Arguments in addition to those of QNMDIF, intended to enable
C     making partial steps in a run and restarting from them, and
C     to eliminate local storage:
C
C     Argument Dim Type  I/O/S  Description
C
C     OLDG      N    R       S  Workspace used by QNMDIF2
C     P         N    R       S  "    "    "    "  QNMDIF2
C     FFWD      N    R       S  "    "    "    "  APROXG
C     SH        N    R       S  Copy of original Hs; may be used by APROXG
C     Y         N    R       S  Workspace used by QNMDIF2
C     Z         N    R       S  Tentative new X used by LINSCH
C     PWORK     N    R       S  Workspace used by UPCHOL
C     LFAIL     N    L       S  "    "    "    "  APROXG
C     GTP            R       S  Are these really needed?
C     GNM            R       S
C     ZETA           R   I      The factor by which the the step
C                               size is reduced in APROXG.
C     NALLOW         I   I      Max. # fn. evaluations allowed.
C
C     IGRAD          I   I      Indicates the application's method of
C                               obtaining gradient information:
C                               IGRAD = 0 means normal finite differencing;
C                               IGRAD = 1 means an alternative means (still
C                                         employing QNMDIF2's finite diff-
C                                         erencing, but inappropriate for
C                                         its option to switch from forward
C                                         to central differences).
C
C     ISTART         I   I      Controls some (re)starting choices,
C                               mostly within the application:
C
C                               ISTART = 0 with UNITL = T means start
C                                          with the identity matrix for
C                                          the Hessian approximation.
C                               ISTART = 1 means the same as ISTART = 0
C                                          except the application reads a
C                                          partially optimized geometry.
C                               ISTART = 2 means in addition that the
C                                          application reads previously
C                                          saved gradient & Hessian info.
C                               ISTART = 3 means in addition that stored
C                                          Hs are also read and used.
C     
C     FUN                       Objective function routine:
C
C                               CALL FUN (N, X, F, NCALL, FAIL)
C
C                               NCALL tells the application what situation
C                               the function evaluation is needed for:
C
C                               NCALL = -3: Initial fn. evaluation (before
C                                           the call to QNMDIF2);
C                               NCALL = -2: "Central" fn. evaluation - a
C                                           repeat of the best following a
C                                           successful line search; may be
C                                           invoked from USER for IENTRY=1;
C                               NCALL = -1: Local search evaluation;
C                               NCALL =  0: Line search evaluation;
C                               NCALL >  0: Gradient element number.
C
C     History:
C
C        Date           Initials   Description
C            ?   1972              Initial coding in ALGOL at the NPL
C                                     (Gill, Murray, and Pitfield)
C            ?                     FORTRAN translation, Stanford Univ.
C                                     (Margaret Wright)
C        20 Sep. 1982   RAK        Version 1.1, modified as noted for
C                                     NASA release by Robert Kennelly.
C        13 Aug. 1986   RAK        Version 1.2 - mostly cosmetic mods.
C                                     Included CENDIF/FDSTEP/DIFFER for
C                                     stepsize selection (though they must
C                                     be called separately).  Added PARAMETER
C                                     statements for local array dimensions,
C                                     with test that N <= MAXN on entry.
C                                     Increase MAXN from 30 to 40.  Some
C                                     protection added in LINSCH.
C        Pre-Aug. '94 J.Reuther Various mods. for recovering from a
C                               bad function evaluation.
C        Aug. 1994  D.Saunders  Effort to recover the generality, etc.
C                               How come the argument description is gone????
C        Feb. 1995   JJR/DAS    IGRAD argument introduced to indicate
C                               use of an alternative gradient scheme.
C        Apr. 1995     DAS      All local arrays are now provided as
C                               arguments to avoid "MAXVAR" include file.
C        Oct. 1995      "       APROXG printing streamlined with large
C                               numbers of variables in mind.  Restored
C                               QNMDIF argument descriptions.
C        Dec. 1995      "       ALPMIN = Hmax, not Hmax ** 2./3., to
C                               accept otherwise "unsuccessful" searches.
C                               Avoid repeating a steepest descent search.
C                               Print line-search convergence criteria
C                               if N >= 10.
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   N, NLDIM, NFTOTL, KITER, NTYPE, LUNLOG, NALLOW, NCALL, IGRAD,
     >   ISTART
      REAL
     >   D (N), EL (NLDIM), FFWD (N), G (N), H (N), OLDG (N), P (N),
     >   PWORK (N), SH (N), X (N), Y (N), Z (N),
     >   F, FM, ETA, TOL, STEPMX, EPSMCH, EPSOBJ, GTP, GNM, ZETA
      LOGICAL
     >   CONTIN, CONV, LFAIL (N), LOCAL, PRINT, RESCUE, UNITL
      EXTERNAL
     >   FUN, USER

C     Local constants.

      REAL
     >   ONE, TINY, TWO3RDS, ZERO
      PARAMETER
     >  (ONE  = 1.E+0, TINY = 1.E-30, TWO3RDS = 2.E+0 / 3.E+0,
     >   ZERO = 0.E+0)

C     Local variables.

      INTEGER
     >   I, IENTRY, IFAIL, IN, IR, IS, IT, ITYPE, K, KGRADS, KN, NITMAX
      REAL
     >   ALPHA, ALPMIN, BL, BOUNDK, DELTAF, DMAX, EPS10, FNEW, GNMSQ,
     >   GTPFIN, HMAX, OLDF, PNORM, RTEPS, SFTBND, STEPM1, SUM, T, U, V
      LOGICAL
     >   ALREADY, COUNT, DONE, FAIL, FINAL, STEEPEST, SUCCES

C     Storage.

      SAVE

C     Execution.


      NITMAX = KITER
      KITER  = 0
      KGRADS = 0
      FNEW   = F
      OLDF   = F
      STEPM1 = STEPMX
      CONV   = .TRUE.
      FINAL  = .FALSE.
      COUNT  = .FALSE.
      BOUNDK = 0.01E+0 / (SQRT (REAL (N)) * EPSMCH)

      DO 5 I = 1, N
         SH(I) = H(I)
    5 CONTINUE

      IF (EPSOBJ .LT. EPSMCH) EPSOBJ = EPSMCH
      RTEPS = SQRT (EPSOBJ)

      ALREADY = .FALSE.
      IF (ISTART .EQ. 0 .AND. UNITL) THEN
         DO 10 I = 1, N
            G(I) = ZERO
            D(I) = ONE
   10    CONTINUE
         DO 20 I = 1, NLDIM
            EL(I) = ZERO
   20    CONTINUE
         IF (.NOT. RESCUE) THEN
            CALL APROXG (N, NTYPE, X, G, H, FFWD, FUN, FNEW, NFTOTL,
     >                   NALLOW, PRINT, LUNLOG, ZETA, SH, LFAIL, TINY,
     >                   EPSOBJ)
            KGRADS = KGRADS + 1
            ALREADY = .TRUE.    ! Already steepest descent - can't switch below
         END IF
      END IF

      HMAX = ZERO
      DO 40 I = 1, N
         U = ABS (H(I))
         HMAX = MAX (U, HMAX)
         P(I) = ZERO
   40 CONTINUE

      ALPMIN = HMAX !!! ** TWO3RDS !!! Too many "unsuccessful" line searches

      IF (RESCUE) GO TO 290

C     --------------------
C     Main iteration loop.
C     --------------------

   50 CONTINUE

      OLDF = FNEW
      DO 60 I = 1, N
         OLDG (I) = G (I)
         Y (I) = X (I)
   60 CONTINUE

      P (1) = -OLDG (1)
      IF (N .GT. 1) THEN
         IR = 1
         DO 80 I = 2, N
            SUM = -OLDG (I)
            IT  = I - 1
            DO 70 K = 1, IT
               SUM = SUM - P (K) * EL (IR)
               IR  = IR + 1
   70       CONTINUE
            P (I) = SUM
   80    CONTINUE
      END IF

      P (N) = P (N) / D (N)
      IF (N .GT. 1) THEN
         DO 110 I = 2, N
            IN = N + 1 - I
            IR = IR - 1
            IS = IR
            IT = IN + 1
            SUM = P (IN) / D (IN)
            DO 100 K = IT, N
               KN = N + IT - K
               SUM = SUM - P (KN) * EL (IS)
               IS = IS + 2 - KN
  100       CONTINUE
            P (IN) = SUM
  110    CONTINUE
      END IF

      IFAIL = 0
      STEPMX = STEPM1

      STEEPEST = .FALSE.
      IF (ALREADY) THEN
         STEEPEST = .TRUE.
         ALREADY  = .FALSE.
      END IF

  125 CONTINUE

      SUM = ZERO
      GTP = ZERO
      GNM = ZERO

      DO 130 I = 1, N
         SUM = SUM + P (I) ** 2
         GTP = GTP + G (I) * P (I)
         GNM = GNM + G (I) ** 2
  130 CONTINUE

      PNORM = SQRT (SUM) + TINY
      GNM   = SQRT (GNM) + TINY

      IF (ABS (GTP) .LT. TINY) GTP = SIGN (TINY, GTP)

      DELTAF = MAX (ABS (FNEW - FM), EPSOBJ)
      ALPHA = MIN (-2.0E+0 * DELTAF / GTP, ONE)

      T = EPSOBJ / PNORM
      SFTBND = ZERO
      IF (NTYPE .EQ. 0) SFTBND = ALPMIN / PNORM

      CALL LINSCH (N, FUN, RTEPS, T, ETA, SFTBND, STEPMX / PNORM, P, X,
     >             Z, FNEW, ALPHA, GTP, NFTOTL, NALLOW, SUCCES, PRINT,
     >             LUNLOG, FAIL)

      IF (SUCCES .OR. NTYPE.EQ.1 .OR. IFAIL.EQ.3) KITER = KITER + 1

      COUNT = KITER .GE. NITMAX

      IF (SUCCES) THEN
         STEPMX = STEPM1
         GO TO 150
      END IF

      WRITE (LUNLOG, 1000)
 1000 FORMAT (/, ' QNMDIF2: Unsuccessful line search.')

      IF (FAIL) THEN
         IF (IFAIL .LT. 3) THEN
            WRITE (LUNLOG, '(A, A)')
     >         ' Line search failure was probably due to function',
     >         ' failure.  Attempting to reduce STEPMX.'
            WRITE (LUNLOG, '(A, E20.10)') ' Previous STEPMX: ', STEPMX
            STEPMX = STEPMX * 0.5
            IFAIL = IFAIL + 1
            WRITE (LUNLOG, '(A, E20.10)') ' New STEPMX     : ', STEPMX
            GO TO 125
         ELSE
            WRITE (LUNLOG, '(A, A)')
     >         ' Three reductions in STEPMX failed to fix the function',
     >         ' problems.  STEPMX is being reset.', ' '
            STEPMX = STEPM1
         END IF
      END IF

      IF (.NOT. STEEPEST) THEN
         STEEPEST = .TRUE.
         WRITE (LUNLOG, '(/, A)')
     >      ' QNMDIF2: Resorting to steepest descent.'
         DO I = 1, N
            P (I) = -G (I)
         END DO
         GO TO 125
      ELSE
         WRITE (LUNLOG, '(/, A, A)') ' QNMDIF2: No new minumum found',
     >      ' after switching to steepest descent.'
      END IF

      IF (NTYPE .EQ. 1 .OR. IGRAD .EQ. 1 .OR. FAIL) GO TO 290

      NTYPE = 1
      IF (ISTART .LT. 2 .AND. .NOT. UNITL) THEN
         ITYPE = 1
      ELSE
         ITYPE = 2
      END IF

      CALL APROXG (N, ITYPE, X, G, H, FFWD, FUN, FNEW, NFTOTL, NALLOW,
     >             PRINT, LUNLOG, ZETA, SH, LFAIL, TINY, EPSOBJ)
      GO TO 50

  150 CONTINUE

      IF (PRINT) THEN
         WRITE (LUNLOG, 1010)
 1010    FORMAT ('1')

         CALL OUTPUT2 (N, NLDIM, Y, G, H, P, OLDF, ALPHA, D, EL,
     >                 KITER, NFTOTL, NALLOW, FINAL, LUNLOG)
         WRITE (LUNLOG, '(A)')
      END IF

      IF (IGRAD .EQ. 1 .AND. NTYPE .GE. 1) THEN
         WRITE (LUNLOG, '(A)') ' QNMDIF2: Mismatched IGRAD and NTYPE.'
         WRITE (LUNLOG, '(A, I6)') ' IGRAD: ', IGRAD, ' NTYPE: ', NTYPE
         STOP
      END IF

      IENTRY = 1
      CALL USER (IENTRY, N, NLDIM, X, FNEW, G, H, EL, D, OLDG, P, GTP,
     >           GNM, NFTOTL, NALLOW, NCALL, KITER, NTYPE, CONTIN,
     >           PRINT, LUNLOG)

      IF (COUNT .AND. .NOT. CONTIN) GO TO 290


C     -------------------
C     Calculate gradient.
C     -------------------

C     Check whether to switch back to forward differences following
C     use of central differences for a sufficiently big step:

      IF (NTYPE .EQ. 1 .AND. ALPHA .GT. ALPMIN / PNORM) NTYPE = 0

      CALL APROXG (N, NTYPE, X, G, H, FFWD, FUN, FNEW, NFTOTL, NALLOW,
     >             PRINT, LUNLOG, ZETA, SH, LFAIL, TINY, EPSOBJ)

      GNMSQ = ZERO
      DO 170 I = 1, N
         GNMSQ = GNMSQ + G (I) ** 2
  170 CONTINUE

C     If we are using forward differences and the norm of the gradient
C     is small, recalculate the gradient using central differences:

      IF (GNMSQ .LE. HMAX .AND. NTYPE .EQ. 0 .AND. IGRAD .EQ. 0) THEN
         NTYPE = 1
         CALL APROXG (N, 2, X, G, H, FFWD, FUN, FNEW, NFTOTL, NALLOW,
     >                PRINT, LUNLOG, ZETA, SH, LFAIL, TINY, EPSOBJ)

         GNMSQ = ZERO
         DO 180 I = 1, N
            GNMSQ = GNMSQ + G (I) ** 2
  180    CONTINUE
      END IF

C     ---------------
C     Update Hessian.
C     ---------------

      GTPFIN = ZERO
      DO 200 I = 1, N
         GTPFIN = GTPFIN + G (I) * P (I)
  200 CONTINUE

      V = ALPHA * (GTPFIN - GTP)
      IF (V .LE. ZERO) THEN
         IF (PRINT) WRITE (LUNLOG, 1020) KITER, GTP, GTPFIN
 1020    FORMAT (/, ' QNMDIF2: Warning - on iteration no.', I6,
     >      ', |GTP| did not decrease.', /,
     >      10X, 'GTP (initial):', E25.15, /,
     >      10X, 'GTP (final):  ', E25.15, /,
     >      10X, 'Updates to the Hessian will be bypassed.')
         GO TO 270
      END IF

      V = SQRT (V)
      DO 220 I = 1, N
         P (I) = (G (I) - OLDG (I)) / V
  220 CONTINUE

      CALL C1MDCH (N, NLDIM, D, EL, ONE, P, IFAIL)

      IF (IFAIL .NE. 0) GO TO 290

      V = SQRT (ABS (GTP))
      DO 230 I = 1, N
         OLDG (I) = -OLDG (I) / V
  230 CONTINUE

      CALL UPCHOL (N, NLDIM, EPSMCH, D, EL, OLDG, PWORK)

      DMAX = D (1)
      DO 240 I = 2, N
         DMAX = MAX (D (I), DMAX)
  240 CONTINUE
      BL = DMAX / BOUNDK
      DO 260 I = 1, N
         IF (D (I) .LT. BL)  D (I) = BL
  260 CONTINUE

  270 CONTINUE

      IENTRY = 2
      CALL USER (IENTRY, N, NLDIM, X, FNEW, G, H, EL, D, OLDG, P, GTP,
     >           GNM, NFTOTL, NALLOW, NCALL, KITER, NTYPE, CONTIN,
     >           PRINT, LUNLOG)

C     Check convergence.
C     ------------------

      SUM = ZERO
      DO 280 I = 1, N
         SUM = SUM + X (I) ** 2
  280 CONTINUE

      IF (GNMSQ .LT. TOL ** TWO3RDS .AND.
     >   ALPHA * PNORM .LT. SQRT (EPSMCH * SUM) + TOL) GO TO 300
      IF (OLDF .GT. FNEW .AND. .NOT. COUNT) GO TO 50

  290 CONTINUE

      CONV = .FALSE.
  300 CONTINUE

      IF (.NOT. (LOCAL .OR. RESCUE) .OR. COUNT) GO TO 500

C     -------------
C     Local search.
C     -------------

      IF (PRINT) WRITE (LUNLOG, 1030)
 1030 FORMAT (//, ' QNMDIF2: Begin local search procedure.')
      U = ALPMIN
      EPS10 = 10.E+0 * EPSOBJ * ABS (FNEW)

  310 CONTINUE
      DO 320 I = 1, N
         Y (I) = X (I) + U
  320 CONTINUE

      NCALL = -1
      CALL FUN (N, Y, V, NCALL, FAIL)

      IF (FAIL) THEN
         WRITE (LUNLOG, '(A)')
     >      ' QNMDIF2: Function failure during the local search.',
     >      ' The program will attempt to continue if RESCUE = .TRUE.',
     >      ' but it may fail. Either way, the program is exiting the',
     >      ' local search.'
        V = FNEW
        GO TO 480
      END IF
      NFTOTL = NFTOTL + 1
      WRITE (LUNLOG, '(A, I6)')
     >   ' QNMDIF2: Total # function evaluations so far:', NFTOTL
      IF (NFTOTL .GE. NALLOW) THEN
         WRITE (LUNLOG, '(/, A)') ' Max. # reached - aborting.'
         STOP
      END IF

      IF (ABS (V - OLDF) .LT. EPS10) THEN
         U = 5.0E+0 * U
         GO TO 310
      END IF

      P (1) = U
      IF (N .GT. 1) THEN
         DO 340 I = 2, N
            P (I) = -P (I - 1)
  340    CONTINUE
         IF (MOD (N, 2) .EQ. 1)  P (N) = ZERO
      END IF

      U = SQRT (ALPMIN)
      EPS10 = EPS10 * ABS (V)

      DO 360 I = 1, N
         OLDG (I) = Y (I) + U * P (I)
  360 CONTINUE

      NCALL = -1
      CALL FUN (N, OLDG, OLDF, NCALL, FAIL)

      IF (FAIL) THEN
         WRITE (LUNLOG, '(A)')
     >      ' QNMDIF2: Function failure during the local search.',
     >      ' The program will attempt to continue if RESCUE = .TRUE.',
     >      ' but it may fail. Either way, the program is exiting the',
     >      ' local search.'
        V = FNEW
        GO TO 480
      END IF

      NFTOTL = NFTOTL + 1
      WRITE (LUNLOG, '(A, I6)')
     >   ' QNMDIF2: Total # function evaluations so far:', NFTOTL
      IF (NFTOTL .GE. NALLOW) THEN
         WRITE (LUNLOG, '(/, A)') ' Max. # reached - aborting.'
         STOP
      END IF

      GTP = OLDF - V

      IF (ABS (GTP) .LT. EPS10) THEN
         U = 5.0E+0 * U
         GO TO 310
      END IF

      GTP = - ABS (GTP / U) - ONE
      ALPHA = U
      IF (OLDF .GT. V) THEN
         V = OLDF
         DO 380 I = 1, N
            Y (I) = OLDG (I)
            P (I) = -P (I)
  380    CONTINUE
      END IF

      SUM = ZERO
      DO 400 I = 1, N
         SUM = SUM + P (I) ** 2
  400 CONTINUE

      PNORM = SQRT (SUM) + TINY
      SFTBND = ZERO
      T = RTEPS / PNORM

      CALL LINSCH (N, FUN, RTEPS, T, EPSMCH, SFTBND, STEPMX / PNORM, P,
     >             Y, Z, V, ALPHA, GTP, NFTOTL, NALLOW, SUCCES, PRINT,
     >             LUNLOG, FAIL)

      DO 410 I = 1, N
         P (I) = X (I) - Y(I)
         Y (I) = X (I)
  410 CONTINUE

      IF (V .GE. FNEW) THEN
         SUM = ZERO
         DO 420 I = 1, N
            SUM = SUM + P (I) ** 2
  420    CONTINUE
         PNORM = SQRT (SUM) + TINY
         U = HMAX / PNORM
         DO 430 I = 1, N
            OLDG (I) = X (I) + U * P (I)
  430    CONTINUE

         NCALL = -1
         CALL FUN (N, OLDG, OLDF, NCALL, FAIL)

         IF (FAIL) THEN
            WRITE (LUNLOG, '(A)')
     >      ' QNMDIF2: Function failure during the local search.',
     >      ' The program will attempt to continue if RESCUE = .TRUE.',
     >      ' but it may fail. Either way, the program is exiting the',
     >      ' local search.'
            V = FNEW
            GO TO 480
         END IF

         NFTOTL = NFTOTL + 1
         WRITE (LUNLOG, '(A, I6)')
     >      ' QNMDIF2: Total # function evaluations so far:', NFTOTL
         IF (NFTOTL .GE. NALLOW) THEN
            WRITE (LUNLOG, '(/, A)') ' Max. # reached - aborting.'
            STOP
         END IF

         GTP = - ABS ((OLDF - FNEW) / U)
      ELSE
         GTP = (V - FNEW) * 4.0E+0 / 3.0E+0
      END IF

      ALPHA = ONE
      IF (V. LT. FNEW .OR. OLDF .GT. FNEW) THEN
         DO 460 I = 1, N
            P (I) = -P (I)
  460    CONTINUE
      END IF

      V = FNEW
      SFTBND = ZERO
      T = RTEPS / PNORM

      CALL LINSCH (N, FUN, RTEPS, T, EPSMCH, SFTBND, STEPMX / PNORM, P,
     >             Y, Z, V, ALPHA, GTP, NFTOTL, NALLOW, SUCCES, PRINT,
     >             LUNLOG, FAIL)

  480 CONTINUE
      KITER = KITER + 2
      COUNT = KITER .GE. NITMAX

      IENTRY = 1
      CALL USER (IENTRY, N, NLDIM, Y, V, G, H, EL, D, OLDG, P, GTP, GNM,
     >           NFTOTL, NALLOW, NCALL, KITER, NTYPE, CONTIN, PRINT,
     >           LUNLOG)
      IENTRY = 2
      CALL USER (IENTRY, N, NLDIM, Y, V, G, H, EL, D, OLDG, P, GTP, GNM,
     >           NFTOTL, NALLOW, NCALL, KITER, NTYPE, CONTIN, PRINT,
     >           LUNLOG)

      IF (PRINT) WRITE (LUNLOG, 1040)
 1040 FORMAT (/, ' QNMDIF2: End local search.')

      IF (V .GE. FNEW) GO TO 500

      FNEW = V
      SUM = ZERO
      DO 490 I = 1, N
         X (I) = Y (I)
         SUM = SUM + Y (I) ** 2
  490 CONTINUE
      DONE = (ALPHA * PNORM .LT. SQRT (EPSMCH * SUM) + TOL) .AND. CONV
      IF ((DONE .OR. COUNT) .AND. .NOT.CONTIN) GO TO 500

      NTYPE = 1
      IF (RESCUE) NTYPE = 0

      CALL APROXG (N, NTYPE, X, G, H, FFWD, FUN, FNEW, NFTOTL, NALLOW,
     >             PRINT, LUNLOG, ZETA, SH, LFAIL, TINY, EPSOBJ)

      IENTRY = 2
      CALL USER (IENTRY, N, NLDIM, X, V, G, H, EL, D, OLDG, P, GTP, GNM,
     >           NFTOTL, NALLOW, NCALL, KITER, NTYPE, CONTIN, PRINT,
     >           LUNLOG)
      IF (DONE .OR. COUNT) GO TO 500
      RESCUE = .FALSE.
      CONV = .TRUE.
      GO TO 50

C     ------------
C     Termination.
C     ------------

  500 CONTINUE
      F = FNEW

      WRITE (LUNLOG, 1010)
      IF (COUNT) WRITE (LUNLOG, 1050)
 1050 FORMAT (/,' QNMDIF2: Maximum number of iterations reached.')
      IF (CONV) WRITE (LUNLOG, 1060)
 1060 FORMAT (/,' QNMDIF2: The convergence test has been satisfied.')
      IF (.NOT. CONV) WRITE (LUNLOG, 1070)
 1070 FORMAT (/,' QNMDIF2 has not satisfied the convergence test.')

      FINAL = .TRUE.
      CALL OUTPUT2 (N, NLDIM, X, G, H, P, FNEW, ALPHA, D, EL, KITER,
     >              NFTOTL, NALLOW, FINAL, LUNLOG)

      RETURN
      END


C+----------------------------------------------------------------------
C
      SUBROUTINE APROXG (N, NTYPE, X, G, H, FFWD, FUN, FNEW, NFTOTL,
     >                   NALLOW, PRINT, LUNLOG, ZETA, SH, LFAIL, TINY,
     >                   EPSOBJ)
C
C        APROXG is called by subroutine QNMDIF2.  It assigns a finite-
C     difference approximation to the gradient of the function F (X)
C     the array G (1:N). The intervals for differencing F along each
C     of the coordinate directions are given in array H (I).  Argument
C     NTYPE determines which type of approximation is given.
C
C        NTYPE=0: Forward differences are used, and the N function
C           values are stored in FFWD (1:N).
C
C        NTYPE=1: Central differences requiring the full 2N function
C           evaluations are computed.
C
C        NTYPE=2: Central differences are evaluated using the forward
C           points previously stored in array FFWD.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   N, NTYPE, NFTOTL, NALLOW, LUNLOG
      REAL
     >   FNEW, ZETA, TINY, EPSOBJ,
     >   G (N), H (N), FFWD (N), X (N), SH (N)
      LOGICAL
     >   PRINT, LFAIL (N)
      EXTERNAL
     >   FUN

C     Local variables.

      INTEGER
     >   I, IFAIL, NCALL
      REAL
     >   BNDUP, CFPCD, CFPFD, FUNX, HI, XI
      LOGICAL
     >   FAIL

C     Storage.

      SAVE

C     Execution.
C     ----------

      IF (PRINT) THEN
         IF (NTYPE .EQ. 0) THEN
            WRITE (LUNLOG, 1000) 'forward', NTYPE
         ELSE 
            WRITE (LUNLOG, 1000) 'central', NTYPE
         END IF
      END IF

 1000 FORMAT (/, ' APROXG: Calculate gradient by ', A7, ' differences',
     >   ' (NTYPE = ', I1, ')')

      BNDUP = .25
      IF (NTYPE .EQ. 2) GO TO 40

C     Compute forward points.

      DO 20 I = 1, N
         IFAIL = 0
         LFAIL (I) = .FALSE.
         XI = X (I)

   15    CONTINUE

         IFAIL = IFAIL + 1
         X (I) = XI + H (I)
         NCALL = I

         CALL FUN (N, X, FUNX, NCALL, FAIL)

         NFTOTL = NFTOTL + 1
         IF (NFTOTL .GE. NALLOW) THEN
            WRITE (LUNLOG, '(/, A)')
     >         ' APROXG: Max. # function evaluations reached.'
            STOP
         END IF

         IF (FAIL) THEN
            IF (IFAIL .GE. 6) THEN
               WRITE (LUNLOG, '(/, A, I5)')
     >            ' APROXG: Too many reductions in H for variable ', I,
     >            ' Set gradient to zero and H to original value.'
               FUNX = FNEW
               H (I) = SH (I)
               LFAIL (I) = .FALSE.
            ELSE
               WRITE (LUNLOG, '(/, A, A, I5)')
     >            ' APROXG: Function failure in forward diff. ',
     >            'for design variable', I,
     >            ' Attempting to reduce step size.'
               WRITE (LUNLOG, '(A, E25.15)') ' Old step size: ', H (I)
               H (I) = H (I) / ZETA
               WRITE (LUNLOG, '(A, E25.15)') ' New step size: ', H (I)
               LFAIL (I) = FAIL
               GO TO 15
            END IF
         END IF
         FFWD (I) = FUNX
         X (I) = XI
   20 CONTINUE

      IF (NTYPE .GT. 0) GO TO 40


C     Compute forward difference gradient.

      DO 30 I = 1, N
         G(I) = (FFWD (I) - FNEW) / H (I)

         IF (LFAIL (I)) THEN
            CFPFD = 2.0 * EPSOBJ / MAX (H (I) * ABS (G (I)), TINY)
            IF (CFPFD .GT. BNDUP) THEN
               WRITE (LUNLOG, '(/, A, I5)')
     >           ' APROXG: Cancellation error too high for variable', I,
     >           ' Set gradient to zero and H(I) to original value.'
               H (I) = SH (I)
               G (I) = 0.0
               FFWD (I) = FNEW
            END IF
         END IF 
   30 CONTINUE

      WRITE (LUNLOG, '(/, A)') ' APROXG: Forward difference gradient:',
     >   '     I            G (I)            H (I)'
      WRITE (LUNLOG, '(I6, 2E17.8)') (I, G (I), H (I), I = 1, N)
      GO TO 90


   40 CONTINUE

C     Compute backward points and central difference gradient.

      WRITE (LUNLOG, '(/, A)') ' APROXG: Central difference gradient:',
     >   '     I            G (I)           HF (I)           HB (I)'

      DO 50 I = 1, N
         IFAIL = 0
         LFAIL (I) = .FALSE.
         HI = H (I)
         XI = X (I)

   45    CONTINUE

         IFAIL = IFAIL + 1
         X (I) = XI - H (I)
         NCALL = I

         CALL FUN (N, X, FUNX, NCALL, FAIL)

         NFTOTL = NFTOTL + 1
         IF (NFTOTL .GE. NALLOW) THEN
            WRITE (LUNLOG, '(/, A)')
     >         ' APROXG: Max. # function evaluations reached.'
            STOP
         END IF
 
        IF (FAIL) THEN
            IF (IFAIL .GE. 6) THEN
               WRITE (LUNLOG, '(/, A, I5)')
     >            ' APROXG: Too many reductions in H for variable ', I,
     >            ' Set gradient to zero and H to original value.'
               FUNX = FFWD (I)
               H (I) = SH (I)
               LFAIL (I) = .FALSE.
            ELSE
               WRITE (LUNLOG, '(/, A, A, I5)')
     >            ' APROXG: Function failure in backward diff. ',
     >            'for design variable', I,
     >            ' Attempting to reduce step size.'
               WRITE (LUNLOG, '(A, E25.15)') ' Old step size: ', H (I)
               H (I) = H (I) / ZETA
               WRITE (LUNLOG, '(A, E25.15)') ' New step size: ', H (I)
               LFAIL (I) = FAIL
               GO TO 45
            END IF
         END IF

         G (I) = (FFWD (I) - FUNX) / (HI + H (I))
         X (I) = XI

         IF (LFAIL (I)) THEN
            CFPCD = EPSOBJ / MAX (H (I) * ABS (G (I)), TINY)
            IF (CFPCD .GT. BNDUP) THEN
               WRITE (LUNLOG, '(/, A, I5)')
     >           ' APROXG: Cancellation error too high for variable', I,
     >           ' Set gradient to zero and H(I) to original value.'
               H (I) = SH (I)
               G (I) = 0.
            END IF
         END IF 

         WRITE (LUNLOG, '(I6, 3E17.8)') I, G (I), HI, H (I)
   50 CONTINUE

   90 WRITE (LUNLOG, '(/, A, I6)')
     >   ' APROXG: Total no. of function evaluations so far:', NFTOTL

      RETURN
      END


C+----------------------------------------------------------------------
C
      SUBROUTINE C1MDCH (N, NLDIM, D, EL, ALPHA, Z, IFAIL)
C
C
C        The subroutine C1MDCH updates the Cholesky factorization of the
C     matrix EL * D * ELT  + ALPHA * Z * ZT, where EL is a unit lower
C     triangular matrix, D is a diagonal matrix with positive elements,
C     ELT is the transpose of EL, ALPHA is a positive scalar, Z is a
C     vector, and ZT is the transpose of Z, so that ALPHA * Z * ZT is
C     a positive rank-one change.
C
C        The algorithm used is described in Gill, Golub, Murray, and
C     Saunders, "Methods for Modifying Matrix factorizations", Stanford
C     Computer Science Report CS-72-322, Stanford University, Dec. 1972.
C
C        N is the dimension of the matrices and vectors, NLDIM is
C     N * (N-1)/ 2, the length of the vector EL in which the strict
C     lower triangle of the matrix EL is stored by rows, D is a
C     vector of length N containing the elements of the matrix D,
C     EL is a vector of length NLDIM containing the elements of EL,
C     ALPHA is the positive scalar multiple of the rank-one change,
C     Z is the vector involved, and IFAIL is an integer flag set
C     to zero if there are no errors, and to 1 if overflow occurred.
C
C        Both EL and D are overwritten with the modified factors of the
C     altered matrix, and the values of ALPHA and Z are not retained.
C
C        The subroutine sets IFAIL to zero if there are no problems.
C     IFAIL is set to 1 if any element of the diagonal formed is zero
C     or negative, or when the ratio of the new diagonal element to
C     the old is too large and would cause overflow.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   N, NLDIM, IFAIL
      REAL
     >   D (N), EL (NLDIM), ALPHA, Z (N)

C     Local constants.

C     RMAX should be set to the largest positive floating point value such
C     that +RMAX and -RMAX can both be represented in the computer.  For
C     convenience, the smallest useful value (VAX) has been used - it is
C     unlikely to matter, but could be increased on another system.

      REAL
     >   RMAX
      PARAMETER
     >  (RMAX = 1.E+38)

C     Local variables.

      INTEGER
     >   I, IB, J, K
      REAL
     >   A, BETA, DB, DI, GAMMA, P1, T, W

C     Storage.

      SAVE

C     Execution.

      A = ALPHA
      K = 0
      IFAIL = 1
      DO 70 I = 1, N
         P1 = Z (I)
         DI = D (I)
         T = A * P1
         D (I) = DI + T * P1
         DB = D (I)

C        Exit with IFAIL = 1 if any of the new diagonal elements is
C        non-positive, or if the ratio of diagonals will overflow.

         IF (DB .LT. 1.0E+0) THEN
            IF (DB .LT. 0.0E+0 .OR. DI .GT. RMAX * DB)  RETURN
         END IF

         GAMMA = DI / DB
         BETA = T / DB
         K = K + I
         J = K
         A = A * GAMMA
         IF (I .LT. N) THEN
            IF (GAMMA .LE. 0.25E+0) THEN

C              If ratio of diagonals is less than 4.0, proceed with
C              normal update.

               DO 20 IB = I + 1, N
                  T = EL (J)
                  EL (J) = T * GAMMA + BETA * Z (IB)
                  Z (IB) = Z (IB) - P1 * T
                  J = J + IB - 1
   20          CONTINUE

            ELSE

C              Use alternative update if ratio of diagonals exceeds 4.0.

               DO 40 IB = I + 1, N
                  T = EL (J)
                  Z (IB) = Z (IB) - P1 * T
                  W = Z (IB)
                  EL (J) = BETA * W + T
                  J = J + IB - 1
   40          CONTINUE
            END IF
         END IF
   70 CONTINUE
      IFAIL = 0

      RETURN
      END


C+----------------------------------------------------------------------
C
      SUBROUTINE LINSCH (N, FUN, EPS, T, ETA, SFTBND, STEPMX, P, X, Z,
     >                   F, ALPHA, GTP, NFTOTL, NALLOW, SUCCES, PRINT,
     >                   LUNLOG, FAIL)
C
C        LINSCH is called by QNMDIF2 to execute a linear search along a
C     given direction P, starting from a given point X, to locate an
C     approximation ALPHA to the point at which the objective function
C     attains its minimum value along the direction P.  The method used
C     is that of successive quadratic interpolation with safeguards.
C
C     History:
C
C        12 Aug. 1986  R.Kennelly  Print and test T on entry to monitor
C                                  conflicts with SFTBND and STEPMX.
C                                  Check D and GTP before main loop.
C        Pre-Aug. 1994 J.Reuther   Various mods. for recovering from a
C                                  bad function evaluation.
C        August 1994   D.Saunders  Tried to recover the generality, etc.
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   N, NFTOTL, NALLOW, LUNLOG
      REAL
     >   EPS, T, ETA, SFTBND, STEPMX, P (N), X (N), Z (N), F, ALPHA, GTP
      LOGICAL
     >   SUCCES, PRINT, FAIL
      EXTERNAL
     >   FUN

C     Local variables.

      INTEGER
     >   I, IFAILZ, MAXLIN, NCALL, NLIN, NPCNT
      REAL
     >   A, A1, B, B1, D, D1, D2, E, FA, FACT, FMIN, FOLD, FU, FV, FW,
     >   GTEST1, GTEST2, Q, R, S, SCXBD, T2, TOL, U, V, W, XM, XMIN

C     Storage.

      SAVE

C     Execution.
C     ----------

      WRITE (LUNLOG, 1000) SFTBND, STEPMX, T
 1000 FORMAT (/, ' LINSCH: Begin search.  Min. step = ', E15.8, /
     >   24X, 'Max. step = ', E15.8, /
     >   24X, 'Tolerance = ', E15.8)

      IF (T .GT. STEPMX) THEN
         WRITE (LUNLOG, 1010)
 1010    FORMAT (/, ' LINSCH: Max. step must exceed tolerance!')
         SUCCES = .FALSE.
         GO TO 330
      END IF
      FAIL = .FALSE.
         

C     Initialization.
C     ---------------

C     NLIN is the number of function evaluations during the linear
C     search.  It is used locally as a flag to compute the parabolic
C     step the first time using GTP, the directional derivative along
C     the search vector.

CRAK  MAXLIN is a (fairly loose) upper limit on the number of function
CRAK  evaluations permitted before declaring the search a failure.

      NLIN = 0
      MAXLIN = 20
      V = 0.0E+0
      W = 0.0E+0
      XMIN = 0.0E+0
      FA = F
      FV = F
      FW = F
      FMIN = F
      FOLD = F
      D = ALPHA
      TOL = T
      T2 = TOL + TOL

C     All points in the linear search are scaled (shifted) so that
C     the currently lowest point (XMIN) is the origin.  A and B define
C     the interval of uncertainty, W is the last value of XMIN, and V
C     is the highest of the three points through which a parabola may
C     be fitted.

      A = 0.0E+0
      B = STEPMX + EPS * ABS (STEPMX) + T
      E = B
      SCXBD = STEPMX

C     NPCNT is a local count of the number of non-parabolic interpolation
C     steps used in LINSCH.  It may be useful for diagnostic purposes if
C     there are problems locating the minimum of a particular function.
C     GTEST1 and GTEST2 are used to test convergence.

      NPCNT = 0
      B1 = B
      A1 = A
      GTEST1 = -1.E-4 * GTP
      GTEST2 = -ETA * GTP
      ALPHA = 0.0E+0

C     D is the estimate of the next point at which the function is to be
C     evaluated, with respect to the current origin.  Make sure that it,
C     and GTP, are sensible before proceeding.

      IF (D .LE. 0.0E+0 .OR. GTP .GE. 0.0E+0) THEN
         WRITE (LUNLOG, 1020) D, GTP
 1020    FORMAT (/, ' LINSCH: Bad initial data: D = ', E25.15,
     >      ', GTP = ', E25.15)
         SUCCES = .FALSE.
         GO TO 330
      END IF


C     Begin main loop.
C     ----------------

   10 CONTINUE
      IF (D .GE. SCXBD) THEN

C        D exceeds the shifted bound on the minimum, so adjust D and
C        the bound.

         D = SCXBD
         SCXBD = (SCXBD - TOL) / (1.0E+0 + EPS)
      END IF

C     U is the point at which the function will actually be evaluated,
C     scaled with respect to the shifted origin.  Thus XMIN, the actual
C     step from the original starting point, must be added to obtain
C     the true step.  The estimate D will be used as the value for U
C     only if it is sufficiently far from zero (the current minimum),
C     since the function is not allowed to be evaluated at points closer
C     together than TOL.

      U = D
      IF (ABS (D) .LT. TOL) U = SIGN (TOL, D)
      R = XMIN + U
      DO 30 I = 1, N
         Z (I) = X (I) + R * P (I)
   30 CONTINUE


C     Evaluate function at new point.
C     -------------------------------

      IF (PRINT) WRITE (LUNLOG, 1030) R
 1030 FORMAT (/, ' LINSCH: Step =', E22.15)

      NCALL = 0
      CALL FUN (N, Z, FU, NCALL, FAIL)

      IF (FAIL) THEN
         WRITE (LUNLOG, '(/, A)')
     >      ' LINSCH: Function failure during line search.'
         WRITE (LUNLOG, '(A, E25.15)')
     >      ' Previous best function: ', FOLD,
     >      ' Current  best solution: ', FMIN,
     >      ' Exiting line search.'
         GO TO 280
      END IF

      NLIN = NLIN + 1

C     Update A, B, V, W, and XMIN.  Check whether new function value is
C     lower than previous minimum.

      IF (FU .GT. FMIN) GO TO 60

C     The following code treats the case where FU is the new lowest
C     point, so that the new point, U, becomes the next origin, and the
C     other points are shifted accordingly.

C     Shift left or right endpoint depending on which side of origin
C     the new point is.

      IF (U .GE. 0.0E+0) THEN
         A = 0.0E+0
         FA = FMIN
      ELSE
         B = 0.0E+0
      END IF


C     Shift all points with respect to new minimum.
C     ---------------------------------------------

      V = W
      FV = FW
      W = 0.0E+0
      FW = FMIN
      FMIN = FU
      XMIN = XMIN + U
      A = A - U
      B = B - U
      V = V - U
      W = W - U
      SCXBD = SCXBD - U
      TOL = EPS * ABS (XMIN) + T
      T2 = 2.0E+0 * TOL
      GO TO 110

C     In the following code, the new function value was greater than
C     the previous lowest.  Check the relationship of the new point
C     to the other values (next lowest, etc.).  The origin remains
C     unchanged, but other points may be interchanged.

C     Shift either the left or right endpoint of the interval of
C     uncertainty, depending on whether the new point, U, was to the
C     left or right of the origin.

   60 IF (U .LE. 0.0E+0) THEN
         A = U
         FA = FU
      ELSE
         B = U
      END IF
      IF (FU .GT. FW .AND. W .NE. 0.0E+0) GO TO 90

C     FU is less than or equal to previous second best point, or W = 0,
C     i.e., this is the first time through this section of code.

      V = W
      FV = FW
      W = U
      FW = FU
      GO TO 100

   90 IF (FU .GT. FV .AND. V .NE. 0.0E+0 .AND. V .NE. W) GO TO 100

C     FU .LE. FW, or V = 0, or V = W (first or second time through).

      V = U
      FV = FU
  100 U = 0.0E+0

C     Compute midpoint of interval of uncertainty.

  110 XM = 0.5E+0 * (A + B)


C     Check stopping criteria.
C     ------------------------

C     Stop if interval of uncertainty is sufficiently small, or
C     if the best point is less than the required bound, or
C     if function value has decreased sufficiently, or
C     if MAXLIN function evaluations have already been performed.

      IF (N .GE. 99) THEN  ! Changed from 10 on 08/22/11
         WRITE (LUNLOG, '(/,A,/,(2(10X,A,E20.10)))')
     >      ' Line search stopping critera:',
     >      'SFTBND:', SFTBND, 'STEPMX:', STEPMX,
     >      'A:     ', A,      'B:     ', B,
     >      'XM:    ', XM,     'XMIN:  ', XMIN,
     >      'FA:    ', FA,     'FMIN:  ', FMIN,
     >      'FOLD:  ', FOLD,   'GTEST2:', GTEST2,
     >      'T2:    ', T2,     'TOL:   ', TOL
      END IF

      IF ((ABS (XM) .LE. (T2 - 0.5E+0 * (B - A))) .OR.
     >   ((XMIN + B) .LE. SFTBND) .OR.
     >   ((FA - FMIN).LE.(ABS (A) * GTEST2) .AND. FMIN.LT.FOLD) .OR.
     >   (NLIN .GE. MAXLIN)) GO TO 280

C     If stopping criteria are not met, continue.

      S = 0.0E+0
      Q = 0.0E+0
      R = 0.0E+0
      IF (ABS (E) .LE. TOL) GO TO 170

C     Otherwise, fit parabola through XMIN, V, and W.

      IF (NLIN .EQ. 1) THEN

C        This is the first parabolic fit and the (known)
C        approximate gradient at the initial point can be used.

         Q = 2.0E+0 * (FW - FMIN - W * GTP)
         IF (XMIN .NE. 0.0E+0) THEN
            S = (2.0E+0 * (FMIN - FW) + W * GTP) * W
         ELSE
            S = GTP * W * W
         END IF

      ELSE  ! NLIN > 1, so the fit uses function values only.

         R = W * (FV - FMIN)
         Q = V * (FW - FMIN)
         S = R * W - Q * V
         Q = 2.0E+0 * (Q - R)
      END IF

      IF (Q .GT. 0.0E+0) THEN
         S = -S
      ELSE
         Q = -Q
      END IF

      R = E
      IF (D .NE. B1 .OR. B .LE. SCXBD) E = D

C     Construct an artificial bound on the estimated step length.

  170 A1 = A
      B1 = B
      IF (XMIN .EQ. 0.0E+0) THEN
         D = XM
         GO TO 230
      END IF

      IF (B .GT. SCXBD) THEN

C        Expand interval by 4 if minimum is still not bracketed.

         D = -4.0E+0 * A
         GO TO 230
      END IF

C     B is less than or equal to SCXBD.

      D1 = A

C     Determine interval of length D2 in which to set artificial bound.

      D2 = B
      IF (ABS (D2) .GT. TOL .AND.
     >  ((W .LE. 0.0E+0) .OR. (ABS (D1) .LE. TOL))) GO TO 200

C     ABS (D2) is less than or equal to TOL, or W is positive and
C     ABS (D1) is greater than TOL.  In either case, interchange D1, D2.

      U = D1
      D1 = D2
      D2 = U

C     Use parabolic interpolation only if new point lies in (A1,B1).

  200 U = -D1 / D2
      IF (U .GE. 1.0E+0) THEN  ! Extrapolation
         FACT = 5.0E+0 * (0.1E+0 + 1.0E+0 / U) / 11.0E+0
      ELSE                     ! Interpolation
         FACT = 0.5E+0 * SQRT (U)
      END IF

      D = D2 * FACT

  230 IF (D .GT. 0.0E+0) THEN
         B1 = D
      ELSE
         A1 = D
      END IF

      IF (ABS (S) .GE. ABS (0.5E+0 * Q * R) .OR. (S .LE. (Q * A1))
     >  .OR. (S .GE. (Q * B1))) GO TO 260

C     A parabolic interpolation step.

      D = S / Q

C     F must not be evaluated too close to A or B.

      IF ((D - A) .GE. T2 .AND. (B - D) .GE. T2) GO TO 10

C     Otherwise, set D to plus or minus TOL.

      D = SIGN (TOL, XM)
      GO TO 10

C     A non-interpolation step.

  260 NPCNT = NPCNT + 1
      IF (XM .GT. 0.0E+0) THEN
         E = B
      ELSE
         E = A
      END IF
      GO TO 10


C     Check safeguards.
C     -----------------

  280 CONTINUE
      SUCCES = .FALSE.
      IFAILZ = 1

C     Check that new point satisfies safeguard conditions.  If the
C     function did not decrease, or step to the minimum was less than
C     SFTBND, LINSCH has failed to locate an acceptable minimum.

  290 CONTINUE
         IF (XMIN .EQ. 0.0) GO TO 330
         IF (XMIN + B .LE. SFTBND) GO TO 330
         IF ((FOLD - FMIN .GT. GTEST1 * XMIN) .AND. (.NOT. FAIL .OR.
     >       IFAILZ .EQ. 1)) GO TO 310
         IF (NLIN .GE. MAXLIN) GO TO 330

C        A sufficiently lower point was not found - try halving step.

         IFAILZ = 0
         XMIN = XMIN * 0.5E+0
         IF (XMIN .LE. T) GO TO 330
         DO 300 I = 1, N
            Z (I) = X (I) + XMIN * P (I)
  300    CONTINUE
         IF (PRINT) WRITE (LUNLOG, 1030) XMIN

         NCALL = 0
         CALL FUN (N, Z, FMIN, NCALL, FAIL)
         NLIN = NLIN + 1
         GO TO 290

C     A sufficiently lower point was found - set output values.

  310 CONTINUE
      SUCCES = .TRUE.
      FAIL  = .FALSE.
      ALPHA = XMIN
      IF (SCXBD .LE. 0.0E+0) ALPHA = STEPMX
      DO 320 I = 1, N
         X (I) = X (I) + ALPHA * P (I)
  320 CONTINUE
      F = FMIN


C     Termination.
C     ------------

  330 CONTINUE
      NFTOTL = NFTOTL + NLIN
      WRITE (LUNLOG, '(/, A, I6)')
     >   ' LINSCH: Total # fn. evaluations after line search: ', NFTOTL
      IF (NFTOTL .GE. NALLOW) THEN
         WRITE (LUNLOG, '(/, A)') ' Max. # reached - aborting.'
         STOP
      END IF

      RETURN
      END


C+----------------------------------------------------------------------
C
      SUBROUTINE OUTPUT2 (N, NLDIM, X, G, H, P, FNEW, ALPHA, D, EL,
     >                    KITER, NFTOTL, NALLOW, FINAL, LUNLOG)
C
C     This routine prints details of an iteration by QNMDIF2.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   N, NLDIM, KITER, NFTOTL, NALLOW, LUNLOG
      REAL
     >   D (N), EL (NLDIM), G (N), H (N), P (N), X (N),
     >   FNEW, ALPHA
      LOGICAL
     >   FINAL

C     Local variables.

      INTEGER
     >   I, J, K1, K2
      REAL
     >   BOUNDK, DMAX, DMIN, GNORM

C     Storage.

      SAVE

C     Execution.

      IF (.NOT. FINAL) THEN
         WRITE (LUNLOG, 1000) KITER, FNEW
 1000    FORMAT (/, ' AT ITERATION', I5, ' THE FUNCTION VALUE IS',
     >      E25.15)
         WRITE (LUNLOG, 1010)
 1010    FORMAT (/, '            CURRENT SOLUTION         GRADIENT',
     >        '                 SEARCH DIRECTION         STEP SIZE')

         WRITE (LUNLOG, 1020) (I, X (I), G (I), P (I), H (I), I = 1, N)
 1020    FORMAT (1X, I5, 4E25.15)
         WRITE (LUNLOG, 1030) ALPHA, NFTOTL
 1030    FORMAT (/, ' LINEAR SEARCH STEP', E25.15, 5X,
     >      'NUMBER OF FUNCTION EVALUATIONS AT END OF ITERATION', I9)

      ELSE

C        Final iteration gets special treatment.

         WRITE (LUNLOG, 1040) KITER, FNEW
 1040    FORMAT (//, ' FINAL ITERATION', I6, '     FUNCTION VALUE',
     >      E25.15)
         WRITE (LUNLOG, 1050)
 1050    FORMAT (/, '            FINAL SOLUTION           GRADIENT')
         WRITE (LUNLOG, 1055) (I, X (I), G (I), I = 1, N)
 1055    FORMAT (1X, I5, 2E25.15)
         WRITE (LUNLOG, 1060) NFTOTL
 1060    FORMAT (/, ' NUMBER OF FUNCTION EVALUATIONS', I9)

      END IF

      WRITE (LUNLOG, 1070)
 1070 FORMAT (/, ' CHOLESKY FACTORS OF APPROXIMATE HESSIAN',
     >        //, '    ELEMENTS OF DIAGONAL MATRIX')
      WRITE (LUNLOG, 1080) (D (I), I = 1, N)
 1080 FORMAT (1X, 5E25.15)

      IF (N .GT. 1 .AND. N .LE. 20) THEN
         WRITE (LUNLOG, 1090)
 1090    FORMAT (/, '    LOWER TRIANGULAR FACTOR')
         K2 = 0
         DO 20 I = 1, N - 1
            K1 = K2 + 1
            K2 = K2 + I
            WRITE (LUNLOG, 1080) (EL (J), J = K1, K2)
   20    CONTINUE
      END IF

CRAK  Compute and print estimate of condition number of Hessian and
CRAK  norm of gradient (all iterations).

      DMAX = D (1)
      DMIN = D (1)
      GNORM = 0.E+0
      DO 70 I = 1, N
         DMAX = MAX (D (I), DMAX)
         DMIN = MIN (D (I), DMIN)
         GNORM = GNORM + G (I) ** 2
   70 CONTINUE
      GNORM = SQRT (GNORM)
      BOUNDK = DMAX / DMIN
      WRITE (LUNLOG, 1100) BOUNDK, GNORM
 1100 FORMAT (/, ' LOWER BOUND ON CONDITION NUMBER OF HESSIAN', E25.15,
     >        //, ' NORM OF GRADIENT', E25.15)

      RETURN
      END


C+----------------------------------------------------------------------
C
      SUBROUTINE UPCHOL (N, NLDIM, EPSMCH, D, EL, Z, P)
C
C
C        The subroutine UPCHOL forms the updated Cholesky factorization
C     of the matrix EL * D * ELT - Z * ZT, where EL is a unit lower
C     triangular matrix, D is diagonal matrix with positive elements,
C     ELT is the transpose of EL, Z is a vector, and ZT is its transpose
C     (so that -Z * ZT is a negative rank-one correction).
C
C        The algorithm used is described in Gill, Golub, Murray, and
C     Saunders, "Methods for Modifying Matrix Factorizations", Stanford
C     Computer Science Report CS-72-322, Stanford University, Dec. 1972.
C
C        N is the dimension of the matrices and vectors, NLDIM is
C     N * (N-1)/ 2, the length of the vector EL in which the strict
C     lower triangle of the matrix EL is stored by rows, EPSMCH is
C     "machine epsilon", used to ensure that the resulting matrix
C     is sufficiently positive definite, D is a vector of length N
C     containing the elements of the matrix D, EL is a vector of length
C     NLDIM containing the elements of EL, and Z is the vector to be
C     used in the update.
C
C       The vectors EL and D are overwritten with the updated Cholesky
C    factors, and the elements of Z are overwritten.  UPCHOL modifies
C    the Cholesky factors of a matrix that could be indefinite, but it
C    ensures that the new matrix is positive definite.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   N, NLDIM
      REAL
     >   EPSMCH, D (N), EL (NLDIM), Z (N), P (N)

C     Local variables.

      INTEGER
     >   I, IB, IQ, J, JJ, K
      REAL
     >   BETA, DJ, G, GAMMA, PJ, T

C     Storage.

      SAVE

C     Execution.

C     Solve EL * P = Z.

      J = 1
      P (1) = Z (1)
      GAMMA = P (1) ** 2 / D (1)

      DO 20 I = 2, N
         K = I - 1
         T = Z (I)
         DO 10 IB = 1, K
            T = T - P (IB) * EL (J)
            J = J + 1
   10    CONTINUE
         P (I) = T
         GAMMA = GAMMA + T * T / D (I)
   20 CONTINUE

C     If 1.0 - GAMMA < EPSMCH, then the modified matrix is not sufficiently
C     positive definite.  GAMMA is replaced by a quantity which ensures
C     that the modified matrix is positive definite regardless of subsequent 
C     rounding error.

      GAMMA = 1.0E+0 - GAMMA
      IF (GAMMA .LE. EPSMCH) THEN
         IF (-GAMMA .LE. EPSMCH) THEN
            GAMMA = EPSMCH
         ELSE
            GAMMA = -GAMMA
         END IF
      END IF

      K = J + N + N

C     Solve D * ELTRANSPOSE * Z = P.

      DO 90 JJ = 1, N
         J = N + 1 - JJ
         PJ = P (J)
         DJ = D (J)
         T = PJ / DJ
         Z (J) = PJ
         BETA = -T / GAMMA
         G = GAMMA + PJ * T
         D (J) = DJ * GAMMA / G
         GAMMA = G
         K = K - J - 1
         IQ = K

         DO 70 IB = J + 1, N
            T = EL (IQ)
            EL (IQ) = T + BETA * Z (IB)
            Z (IB) = Z (IB) + PJ * T
            IQ = IQ + IB - 1
   70    CONTINUE
   90 CONTINUE

      RETURN
      END
