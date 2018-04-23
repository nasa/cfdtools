C+----------------------------------------------------------------------
C
C                          CENDIF2 + QNMDIF2
C
C     (QNMDIF2 has been removed from this file for Traj_opt purposes.)
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


C+------------------------------------------------------------------------------
C
      SUBROUTINE LCSAREAS (NDATA, X, Y, METHOD, OFFSET, AREAS)
C
C  One-liner:  Local cubic spline integration of discrete data (partial sums)
C
C  Description:
C
C        LCSAREAS calculates the definite integrals up to each knot for the
C     local cubic spline defined by discrete data points in 2-space, with an
C     optional constant added to all partial integrals.  Local cubic splines
C     use 4-pt. methods to determine spline coefficients on the fly as needed.
C
C        This variation of LCSQUAD returns all partial integrals on the
C     intervals [X(1), X(I)], I = 1 : NDATA, whereas LCSQUAD returns a single
C     definite integral on the (arbitrary) interval [XA, XB].  The "CYCLIC"
C     option is also omitted here for added efficiency.
C
C  Arguments:
C
C     Name    Dimension  I/O/S  Description
C
C     NDATA              I      Length of X, Y input data arrays.
C
C     X,      (NDATA)    I      Input data coordinates.  The Xs
C     Y                         must be distinct and monotonic,
C                               either ascending or descending.
C                               (No check here.) 
C
C     METHOD  C*1        I      (Uppercase) Type of fit to be used:
C                               'M' means Monotonic piecewise cubics;
C                               'B' means non-monotonic "Bessel"-type
C                                   piecewise cubics (looser fit);
C                               'L' means piecewise Linear fit
C
C     OFFSET             I      Value assigned to AREAS(1);
C                               added to all partial integrals
C
C     AREAS   (NDATA)      O    Integrals of Y spline up to each knot.
C                               AREAS(1) = OFFSET;
C                               AREAS(I) = OFFSET + integral from X(1) to X(I).
C
C  Significant local variables:
C
C     H, DEL         Delta X and forward difference derivative arrays.
C     B, C, D        Coefficients of cubic on the bracketing interval.
C
C  Procedures:
C
C     BESSEL    First derivative (central 3-point formula).
C     BRODLIE   First derivative (central), adjusted for monotonicity.
C     BUTLAND   First derivative (non-central), adjusted for monotonicity.
C     THREEPT   First derivative (non-central 3-point formula).
C
C  Environment:  Fortran 90
C
C  History:
C
C     06/10/92  D.A.Saunders  LCSQUAD  adapted from LCSFIT.
C     06/12/02   "    "    "  LCSAREAS adapted from LCSQUAD and CSQUAD.
C                             Take advantage of marching through the
C                             data intervals by reusing the various
C                             2- and 3-point slope calculations,
C                             at the expense of repeated code.
C
C  Author:  David Saunders, ELORET/NASA Ames, Moffett Field, CA.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER,   INTENT (IN)  :: NDATA
      REAL,      INTENT (IN)  :: X (NDATA), Y (NDATA)
      CHARACTER, INTENT (IN)  :: METHOD * 1
      REAL,      INTENT (IN)  :: OFFSET
      REAL,      INTENT (OUT) :: AREAS (NDATA)

      REAL, PARAMETER ::
     &   ONE    = 1.0E+0,
     &   TWO    = 2.0E+0,
     &   THREE  = 3.0E+0,
     &   HALF   = ONE / TWO,
     &   THIRD  = ONE / THREE,
     &   FOURTH = ONE / (TWO + TWO)

C     Local variables:

      LOGICAL
     &   MONO
      INTEGER
     &   L
      REAL
     &   B (0:1), C, DELY (-1:1), D, DX, H (-1:1)

C     Procedures:

      REAL, EXTERNAL :: BESSEL, BRODLIE, BUTLAND, THREEPT

C     Execution:
C     ----------

C     Accumulate integrals for the cubics in successive intervals.

      AREAS (1) = OFFSET

      IF (METHOD .EQ. 'L' .OR. NDATA == 2) THEN ! Trapezoidal rule

         DO L = 2, NDATA
            AREAS (L) = (Y (L) + Y (L-1)) * HALF *
     &                  (X (L) - X (L-1)) + AREAS (L-1)
         END DO

      ELSE ! NDATA >= 3

         MONO = METHOD .EQ. 'M'

C        Left-most interval:

         H (0)    =  X (2) - X (1)
         DELY (0) = (Y (2) - Y (1)) / H (0) ! 2-pt. slopes
         H (1)    =  X (3) - X (2)
         DELY (1) = (Y (3) - Y (2)) / H (1)

         IF (MONO) THEN                     ! 3-pt. slopes
            B (0) = BUTLAND (0, H, DELY)
            B (1) = BRODLIE (1, H, DELY)
         ELSE
            B (0) = THREEPT (0, H, DELY)
            B (1) = BESSEL  (1, H, DELY)
         END IF

         DX = H (0)
         C  = (THREE * DELY (0) - TWO * B (0) - B (1)) / DX
         D  = ( -TWO * DELY (0) + B (0) + B (1)) / DX ** 2
         AREAS (2) = DX * (Y (1) + DX * (HALF * B (0) +
     &               DX * (THIRD * C + DX * FOURTH * D))) + AREAS (1)

C        Interior intervals:

         DO L = 3, NDATA - 1

            H (-1)    =  H (0)  ! Avoid recalculating earlier differences
            H ( 0)    =  H (1)
            H ( 1)    =  X (L+1) - X (L)
            DELY (-1) = DELY (0)
            DELY ( 0) = DELY (1)
            DELY ( 1) = (Y (L+1) - Y (L)) / H (1)

            IF (MONO) THEN
               B (0) = BRODLIE (0, H, DELY)
               B (1) = BRODLIE (1, H, DELY)
            ELSE
               B (0) = BESSEL  (0, H, DELY)
               B (1) = BESSEL  (1, H, DELY)
            END IF

            DX = H (0)
            C  = (THREE * DELY (0) - TWO * B (0) - B (1)) / DX
            D  = ( -TWO * DELY (0) + B (0) + B (1)) / DX ** 2
            AREAS (L) = DX * (Y (L-1) + DX * (HALF * B (0) +
     &                  DX * (THIRD * C + DX * FOURTH * D))) +
     &                  AREAS (L-1)
         END DO

C        Final interval:

         H (-1)    =  H (0)
         H ( 0)    =  H (1)
         DELY (-1) = DELY (0)
         DELY ( 0) = DELY (1)

         IF (MONO) THEN
            B (0) = BRODLIE (0, H, DELY)
            B (1) = BUTLAND (1, H, DELY)
         ELSE
            B (0) = BESSEL  (0, H, DELY)
            B (1) = THREEPT (1, H, DELY)
         END IF

         DX = H (0)
         C  = (THREE * DELY (0) - TWO * B (0) - B (1)) / DX
         D  = ( -TWO * DELY (0) + B (0) + B (1)) / DX ** 2
         L  = NDATA
         AREAS (L) = DX * (Y (L-1) + DX * (HALF * B (0) +
     &               DX * (THIRD * C + DX * FOURTH * D))) + AREAS (L-1)

      END IF

      END SUBROUTINE LCSAREAS


C+----------------------------------------------------------------------
C
      SUBROUTINE LCSFIT (NDATA, X, Y, NEW, METHOD, NEVAL, XEVAL, YEVAL,
     &   YPEVAL)
C
C     Two-liner:  Storage-efficient local cubic spline fit (2-space)
C     ----------  (monotonic and piecewise linear options too)
C
C     Description and usage:
C     ----------------------
C
C        LCSFIT is the non-parametric analog of PLSFIT (parametric).
C     It is intended for spline applications which do not require the
C     spline coefficients as output.  It is efficient for repeated
C     calls with the same data, so repeated use with NEVAL = 1 may be
C     preferable to storing vectors of results.
C
C        LCSFIT offers monotonic spline and piecewise linear options
C     also.  And it returns an interpolated first derivative along
C     with the function value.  (The second derivative is omitted
C     because Y" is not guaranteed to be continuous by local methods.)
C
C        See PLSFIT for more details on local methods.  As with most
C     numerical methods, scaling of the data to the unit interval (and
C     unscaling of the result) is recommended to avoid unnecessary
C     effects attributable to the data units.  Utilities GETSCALE and
C     USESCALE from the present authors are appropriate.  The data
C     abscissas should be distinct and either ascending or descending.
C     PROTECT is available to check this.  Extrapolation is permitted
C     (mainly in case of round-off; it is normally inadvisable).
C
C        The CSFIT/CSEVAL or CSDVAL pair are probably preferable if
C     efficiency is not an issue, since CSFIT gives Y" continuity.
C
C     Arguments:
C     ----------
C
C     Name    Type/Dimension  I/O/S  Description
C     NDATA   I               I      Length of X, Y input data arrays.
C
C     X,      R (NDATA)       I      Input data coordinates.  The Xs
C     Y                              must be distinct and monotonic,
C                                    either ascending or descending.
C                                    (No check here.) 
C
C     NEW     L               I      If control flag NEW is .TRUE., the
C                                    search for a bracket starts from
C                                    scratch, otherwise locally-saved
C                                    search and fit information will be
C                                    assumed to be correct. If calling
C                                    LCSFIT from within a loop, set
C                                    NEW = .FALSE. after the first call.
C
C     METHOD   C*1            I      (Uppercase) Type of fit to be used:
C                                    'M' means Monotonic piecewise cubics;
C                                    'B' means non-monotonic "Bessel"-type
C                                        piecewise cubics (looser fit);
C                                    'L' means piecewise Linear fit;
C                                    'C' means Cyclic (periodic) end
C                                        conditions: loose fit assumed.
C
C     NEVAL   I               I      Number of interpolations requested.
C                                    NEVAL >= 1.  One call per result
C                                    (NEVAL = 1) may save storage, and is
C                                    not too inefficient as long as NEW
C                                    is set to .FALSE. after the first.
C
C     XEVAL   R (NEVAL)       I      Abscissa(s) to interpolate to.  These
C                                    are normally in the data range, but
C                                    extrapolation - probably due to
C                                    round-off - is not prevented.
C
C     YEVAL   R (NEVAL)       O      Interpolated function value(s).
C
C     YPEVAL  R (NEVAL)       O      Interpolated 1st derivative value(s).
C                                    Pass the same storage as for YEVAL
C                                    if no derivatives are required.
C
C     Significant local variables:
C     ----------------------------
C
C     MEMORY         Indicates that coefficients are correct for the
C                    current point.
C
C     H, DEL         Delta X and forward difference derivative arrays.
C
C     B, C, D        Coefficients of cubic on the bracketing interval.
C
C     Procedures:
C     -----------
C
C     INTERVAL  1-D "interpolation" search.
C     BESSEL    First derivative (central 3-point formula).
C     BRODLIE   First derivative (central), adjusted for monotonicity.
C     BUTLAND   First derivative (non-central), adjusted for monotonicity.
C     THREEPT   First derivative (non-central 3-point formula).
C
C     Environment:  FORTRAN 90
C     ------------
C
C     Error handling:  None
C     ---------------
C
C     Notes:
C     ------
C
C     (1)  Since many of the calculations must be repeated at both ends
C          of an interval, the various finite difference quantities used
C          are stored as arrays. The following "map" of a typical interior
C          interval and its neighbors should help in understanding the
C          notation.  The local array indices are all numbered relative
C          to the left-hand end of the interval which brackets the point
C          to be evaluated.
C
C                                  LEFT       RIGHT
C
C          Point         -1          0         +1          +2
C
C          Data           x -------- x -------- x --------- x
C
C          Interval      -1          0         +1
C
C
C     Author: Robert Kennelly, Sterling Software/NASA Ames  (PLSFIT)
C     -------
C
C     History:
C     --------
C
C     27 Feb. 1987  R.A.Kennelly  Initial implementation of PLSFIT.
C     23 Aug. 1989  D.A.Saunders  LCSFIT adapted as non-parametric form,
C                                 for embedding in other utilities where
C                                 minimizing work-space is desirable.
C     20 June 1991    "    "      THREEPT (monotonic) renamed BUTLAND;
C                                 THREEPT (pure 3-pt. formula) now used
C                                 for nonmonotonic end-point handling;
C                                 METHOD='C' case belatedly added, as
C                                 needed by PLSINTRP for closed curves.
C     23 July 1991    "    "      The tests for being in the same interval
C                                 as before were not allowing for the
C                                 descending-Xs case.
C     06 May  1998    "    "      Minor Fortran 90 updates.
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER,   INTENT (IN)  :: NDATA, NEVAL
      REAL,      INTENT (IN)  :: X (NDATA), Y (NDATA), XEVAL (NEVAL)
      REAL,      INTENT (OUT) :: YEVAL (NEVAL), YPEVAL (NEVAL)
      LOGICAL,   INTENT (IN)  :: NEW
      CHARACTER, INTENT (IN)  :: METHOD * 1

C     Local constants:

      REAL, PARAMETER :: ZERO = 0., ONE = 1., TWO = 2., THREE = 3.

C     Local variables:

      LOGICAL
     &   CYCLIC, LINEAR, MEMORY, MONO
      INTEGER
     &   IEVAL, J, K, LEFT, RIGHT
      REAL
     &   ARROW, B (0:1), C, DELY (-1:1), D, DX, H (-1:1), XBYARROW, XE

C     Procedures:

      REAL, EXTERNAL ::
     &   BESSEL, BRODLIE, BUTLAND, THREEPT

C     Storage:

      SAVE
     &   ARROW, B, C, D, LEFT, RIGHT

C     Execution:
C     ----------

      MONO   = METHOD == 'M'
      CYCLIC = METHOD == 'C'
      LINEAR = METHOD == 'L'

      IF (CYCLIC) THEN
         IF (Y (NDATA) /= Y (1)) STOP 'LCSFIT: End points must match.'
      END IF

C     Initialize search or avoid it if possible:

      IF (NEW) THEN
         MEMORY = .FALSE.
         ARROW  = SIGN (ONE, X (2) - X (1))
         LEFT   = 1
      END IF

      IEVAL = 1
      XE = XEVAL (1)
      XBYARROW = XE * ARROW

      IF (.NOT. NEW) THEN
      
C        We can save a lot of time when LCSFIT is being called from within
C        a loop by setting MEMORY if possible. The out-of-range checking
C        relies on the fact that RIGHT = LEFT + 1 after the last return.
C        Cater to the more likely case of XE in the previous, interior
C        interval.

         MEMORY = XBYARROW >= X (LEFT)  * ARROW .AND.
     &            XBYARROW <  X (RIGHT) * ARROW

         IF (.NOT. MEMORY) THEN
            MEMORY =
     &         LEFT  == 1     .AND. XBYARROW <  X (RIGHT) * ARROW
     &         .OR.
     &         RIGHT == NDATA .AND. XBYARROW >= X (LEFT)  * ARROW
         END IF

      END IF

      IF (MEMORY) GO TO 70 ! Skip the bulk of the computation

C     Loop over evaluation points requiring a new search:
C     ---------------------------------------------------

   10 CONTINUE

C        Interpolation search for bracketing interval:
C        ---------------------------------------------

         CALL INTERVAL (NDATA, X, XE, ARROW, LEFT)

         RIGHT = LEFT + 1

C         -------------------------------------------
C        |                                           |
C        |   1 <= LEFT < RIGHT = LEFT + 1 <= NDATA   |
C        |                                           |
C         -------------------------------------------

C        Compute derivatives by finite-differences:
C        ------------------------------------------

         IF (NDATA > 2 .AND. .NOT. LINEAR) THEN

C           Interval and derivative approximations:
C           ---------------------------------------

C           The following duplicates more code than PLSFIT's approach,
C           but eliminates some indirection - no need to wrap-around here.
C           Handle the end conditions first to minimize testing LEFT, RIGHT.

            IF (LEFT == 1) THEN

               H (0) = X (2) - X (1)
               DELY (0) = (Y (2) - Y (1)) / H (0)
               H (1) = X (3) - X (2)
               DELY (1) = (Y (3) - Y (2)) / H (1)

               IF (CYCLIC) THEN ! Loose fit assumed
                  H (-1) = X (NDATA) - X (NDATA - 1)
                  DELY (-1) = (Y (NDATA) - Y (NDATA - 1)) / H (-1)
                  B (0) = BESSEL (0, H, DELY)
                  B (1) = BESSEL (1, H, DELY)
               ELSE
                  IF (MONO) THEN
                     B (0) = BUTLAND (0, H, DELY)
                     B (1) = BRODLIE (1, H, DELY)
                  ELSE
                     B (0) = THREEPT (0, H, DELY)
                     B (1) = BESSEL  (1, H, DELY)
                  END IF
               END IF

            ELSE IF (RIGHT == NDATA) THEN

               H(-1) = X (LEFT) - X (LEFT-1)
               DELY(-1) = (Y (LEFT) - Y (LEFT-1)) / H (-1)
               H (0) = X (RIGHT) - X (LEFT)
               DELY (0) = (Y (RIGHT) - Y (LEFT))  / H (0)

               IF (CYCLIC) THEN
                  H (1) = X (2) - X (1)
                  DELY (1) = (Y (2) - Y (1)) / H (1)
                  B (0) = BESSEL (0, H, DELY)
                  B (1) = BESSEL (1, H, DELY)
               ELSE

                  IF (MONO) THEN
                     B (0) = BRODLIE (0, H, DELY)
                     B (1) = BUTLAND (1, H, DELY)
                  ELSE
                     B (0) = BESSEL  (0, H, DELY)
                     B (1) = THREEPT (1, H, DELY)
                  END IF
               END IF

            ELSE

               K = LEFT
               DO J = -1, +1
                  H (J)    =  X (K) - X (K-1)
                  DELY (J) = (Y (K) - Y (K-1)) / H (J)
                  K = K + 1
               END DO

C              Select interpolation scheme:
C              ----------------------------

C              Compute (possibly adjusted) first derivatives at both
C              left- and right-hand endpoints of the interval.

               IF (MONO) THEN

C                 Monotone - use Brodlie modification of Butland's
C                 formula to adjust the derivatives at the knots.

                  B (0) = BRODLIE (0, H, DELY)
                  B (1) = BRODLIE (1, H, DELY)

               ELSE ! METHOD = 'B'

C                 Bessel - use central difference formula at the knots.

                  B (0) = BESSEL (0, H, DELY)
                  B (1) = BESSEL (1, H, DELY)

               END IF

            END IF

C           Compute the remaining cubic coefficients.

            C = (THREE * DELY (0) - TWO * B (0) - B (1)) / H (0)
            D = ( -TWO * DELY (0) + B (0) + B (1)) / H (0) ** 2

         ELSE ! NDATA = 2 .OR. METHOD = 'L'

C           Degenerate case (linear).
C           -------------------------

            B (0) = (Y (RIGHT) - Y (LEFT)) / (X (RIGHT) - X (LEFT))
            C     = ZERO
            D     = ZERO

         END IF

C        Evaluate the cubic (derivative first in case only YEVAL is reqd.):
C        ------------------------------------------------------------------

   70    CONTINUE ! Start of same-interval loop inside new-interval loop

            DX = XE - X (LEFT)
            YPEVAL (IEVAL) = B (0) + DX * (TWO * C + DX * THREE * D)
            YEVAL  (IEVAL) = Y (LEFT) + DX * (B (0) + DX * (C + DX * D))

C           The next evaluation point may be in the same interval:
C           ------------------------------------------------------

            IF (IEVAL < NEVAL) THEN ! Skips this if NEVAL = 1

               IEVAL = IEVAL + 1
               XE = XEVAL (IEVAL)
               XBYARROW = XE * ARROW
               IF (XBYARROW >= X (LEFT)  * ARROW  .AND.
     &             XBYARROW <  X (RIGHT) * ARROW) GO TO 70
            
               GO TO 10 ! Else much more work required.

            END IF

C     Termination:
C     ------------

      RETURN

      END SUBROUTINE LCSFIT


C+----------------------------------------------------------------------
C
      SUBROUTINE LCSQUAD (NDATA, X, Y, XA, XB, METHOD, AREA)
C
C  One-liner:  Storage-efficient integration of discrete data
C
C  Description:
C
C        LCSQUAD estimates the integral on [XA, XB] defined by a discrete
C     distribution, using the storage-efficient spline interpolation of
C     LCSFIT.
C
C        Apart from giving results which are smoother with respect to
C     data perturbations than (say) the trapezoidal rule, the monotonic
C     option available for the underlying curve fit guards against possibly
C     wild excursions which may affect other polynomial-based methods.
C        
C  Arguments:
C
C     Name    Type/Dimension  I/O/S  Description
C
C     NDATA   I               I      Length of X, Y input data arrays.
C
C     X,      R (NDATA)       I      Input data coordinates.  The Xs
C     Y                              must be distinct and monotonic,
C                                    either ascending or descending.
C                                    (No check here.) 
C
C     XA,     R               I      Limits of integration, in either
C     XB                             order, possibly different from the
C                                    data abscissas.  The curve fit will
C                                    extrapolate if necessary, so beware.
C
C     METHOD  C*1             I      (Uppercase) Type of fit to be used:
C                                    'M' means Monotonic piecewise cubics;
C                                    'B' means non-monotonic "Bessel"-type
C                                        piecewise cubics (looser fit);
C                                    'L' means piecewise Linear fit;
C                                    'C' means Cyclic (periodic) end
C                                        conditions: loose fit assumed.
C
C     AREA    R                 O    Estimate of the integral of Y between
C                                    XA and XB. 
C
C  Significant local variables:
C
C     H, DEL         Delta X and forward difference derivative arrays.
C     B, C, D        Coefficients of cubic on the bracketing interval.
C
C  Procedures:
C
C     INTERVAL  1-D "interpolation" search utility.
C     BESSEL    First derivative (central 3-point formula).
C     BRODLIE   First derivative (central), adjusted for monotonicity.
C     BUTLAND   First derivative (non-central), adjusted for monotonicity.
C     THREEPT   First derivative (non-central 3-point formula).
C
C  Environment:  Fortran 90
C
C  History:
C
C     06/06/92  D.A.Saunders  Initial version using LCSFIT & QUANC8RC.
C     06/10/92    "    "      Adapted LCSFIT to do it more directly.
C     04/03/97    "    "      J. Reuther found pitched wing sections
C                             gave the wrong sign for ARROW. Now use
C                             X(NDATA) - X(1) instead of X(2) - X(1)
C                             to help wing fuel volume calculations.
C     05/06/98    "    "      Minor Fortran 90 updates.
C
C     Author: David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER,   INTENT (IN)  :: NDATA
      REAL,      INTENT (IN)  :: X (NDATA), Y (NDATA), XA, XB
      CHARACTER, INTENT (IN)  :: METHOD * 1
      REAL,      INTENT (OUT) :: AREA

      REAL, PARAMETER ::
     &   ZERO   = 0.0E+0,
     &   ONE    = 1.0E+0,
     &   TWO    = 2.0E+0,
     &   THREE  = 3.0E+0,
     &   HALF   = ONE / TWO,
     &   THIRD  = ONE / THREE,
     &   FOURTH = ONE / (TWO + TWO)

C     Local variables:

      LOGICAL
     &   CYCLIC, LINEAR, MONO
      INTEGER
     &   J, K, LEFT, LEFTA, LEFTB, RIGHT
      REAL
     &   ARROW, B (0:1), C, DELY (-1:1), D, H (-1:1), SUM, XLEFT, XR,
     &   XRIGHT

C     Procedures:

      REAL, EXTERNAL :: BESSEL, BRODLIE, BUTLAND, THREEPT

C     Statement function:

C     (This is the integral of the cubic for the interval
C     [X (LEFT), X (LEFT+1)] from X (LEFT) to X (LEFT) + DX.)

      REAL
     &   QUAD, DX
      QUAD (DX) = DX * (Y (LEFT) + DX * (HALF * B (0) +
     &            DX * (THIRD * C + DX * FOURTH * D)))

C     Execution:
C     ----------

      MONO   = METHOD .EQ. 'M'
      CYCLIC = METHOD .EQ. 'C'
      LINEAR = METHOD .EQ. 'L'

      IF (CYCLIC) THEN
         IF (Y (NDATA) /= Y (1)) STOP 'LCSQUAD: End pts. must match.'
      END IF

      ARROW = SIGN (ONE, X (NDATA) - X (1))

C     Locate the intervals containing XA and XB:

      LEFTA = 1
      CALL INTERVAL (NDATA, X, XA, ARROW, LEFTA)

      LEFTB = NDATA - 1
      CALL INTERVAL (NDATA, X, XB, ARROW, LEFTB)

C     Ensure that we traverse the X intervals in ascending order
C     (else we get the partial-interval handling messed up):

      IF (LEFTA <= LEFTB) THEN ! Reuse ARROW for sign of integral
         ARROW  = ONE
         XLEFT  = XA
         XRIGHT = XB
      ELSE
         ARROW  = -ONE
         LEFT   = LEFTA
         LEFTA  = LEFTB
         LEFTB  = LEFT
         XLEFT  = XB
         XRIGHT = XA
      END IF

C     Accumulate integrals for the cubics in successive intervals.
C     (Most of this loop is lifted from LCSFIT.)

      SUM = ZERO

      DO LEFT = LEFTA, LEFTB

         RIGHT = LEFT + 1

C        Compute derivatives by finite-differences:
C        ------------------------------------------

         IF (NDATA > 2 .AND. .NOT. LINEAR) THEN

            IF (LEFT == 1) THEN

               H (0) = X (2) - X (1)
               DELY (0) = (Y (2) - Y (1)) / H (0)
               H (1) = X (3) - X (2)
               DELY (1) = (Y (3) - Y (2)) / H (1)

               IF (CYCLIC) THEN ! Loose fit assumed
                  H (-1) = X (NDATA) - X (NDATA - 1)
                  DELY (-1) = (Y (NDATA) - Y (NDATA - 1)) / H (-1)
                  B (0) = BESSEL (0, H, DELY)
                  B (1) = BESSEL (1, H, DELY)
               ELSE
                  IF (MONO) THEN
                     B (0) = BUTLAND (0, H, DELY)
                     B (1) = BRODLIE (1, H, DELY)
                  ELSE
                     B (0) = THREEPT (0, H, DELY)
                     B (1) = BESSEL  (1, H, DELY)
                  END IF
               END IF

            ELSE IF (RIGHT == NDATA) THEN

               H(-1) = X (LEFT) - X (LEFT-1)
               DELY(-1) = (Y (LEFT) - Y (LEFT-1)) / H (-1)
               H (0) = X (RIGHT) - X (LEFT)
               DELY (0) = (Y (RIGHT) - Y (LEFT))  / H (0)

               IF (CYCLIC) THEN
                  H (1) = X (2) - X (1)
                  DELY (1) = (Y (2) - Y (1)) / H (1)
                  B (0) = BESSEL (0, H, DELY)
                  B (1) = BESSEL (1, H, DELY)
               ELSE

                  IF (MONO) THEN
                     B (0) = BRODLIE (0, H, DELY)
                     B (1) = BUTLAND (1, H, DELY)
                  ELSE
                     B (0) = BESSEL  (0, H, DELY)
                     B (1) = THREEPT (1, H, DELY)
                  END IF
               END IF

            ELSE

               K = LEFT
               DO J = -1, +1
                  H (J)    =  X (K) - X (K-1)
                  DELY (J) = (Y (K) - Y (K-1)) / H (J)
                  K = K + 1
               END DO

C              Select interpolation scheme:
C              ----------------------------

C              Compute (possibly adjusted) first derivatives at both
C              left- and right-hand endpoints of the interval:

               IF (MONO) THEN

C                 Monotone - use Brodlie modification of Butland's
C                 formula to adjust the derivatives at the knots.

                  B (0) = BRODLIE (0, H, DELY)
                  B (1) = BRODLIE (1, H, DELY)

               ELSE ! METHOD = 'B'

C                 Bessel - use central difference formula at the knots.

                  B (0) = BESSEL (0, H, DELY)
                  B (1) = BESSEL (1, H, DELY)

               END IF

            END IF

C           Compute the remaining cubic coefficients.

            C = (THREE * DELY (0) - TWO * B (0) - B (1)) / H (0)
            D = ( -TWO * DELY (0) + B (0) + B (1)) / H (0) ** 2

         ELSE ! NDATA = 2 .OR. METHOD = 'L'

C           Degenerate case (linear):
C           -------------------------

            B (0) = (Y (RIGHT) - Y (LEFT)) / (X (RIGHT) - X (LEFT))
            C     = ZERO
            D     = ZERO

         END IF

C        Evaluate the integral of the cubic for this interval:
C        -----------------------------------------------------

C        (Make end adjustments in case XA or XB is not at an X (I).)

         IF (LEFT == LEFTA) SUM = SUM - QUAD (XLEFT - X (LEFT))

         IF (LEFT /= LEFTB) THEN
            XR = X (RIGHT)  ! For complete interior intervals
         ELSE
            XR = XRIGHT     ! No matter where XB is
         END IF

         SUM = SUM + QUAD (XR - X (LEFT))

      END DO
   
      AREA = SUM * ARROW

      RETURN

      END SUBROUTINE LCSQUAD


C+----------------------------------------------------------------------
C
      FUNCTION BESSEL (J, H, DEL)
C
C     One-liner: First derivative using central 3-point formula
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        Computes a first derivative approximation using the central
C     3-point formula.  The data must be in the form of arrays containing
C     finite difference interval lengths and 2-point forward difference
C     derivatives.  BESSEL is intended to be used by PLSFIT for determin-
C     ing end conditions on an interval for (non-monotonic) interpolation
C     by piecewise cubics.  See the PLSFIT header for more details.
C
C     Arguments:
C     ----------
C
C     Name    Type/Dimension  I/O/S  Description
C     J       I               I      Indicates at which end of the
C                                    interval the derivative is to be
C                                    estimated. J = 0 means left-hand
C                                    side, J = 1 means right.
C
C     H       R (-1:1)        I      Array of interval lengths. The 0th
C                                    element is the length of the interval
C                                    on which the cubic is to be deter-
C                                    mined.
C
C     DEL     R (-1:1)        I      Array of derivative estimates. The
C                                    0th element is the forward difference
C                                    derivative over the interval on which
C                                    the cubic is to be determined.
C                                     
C     BESSEL  R                 O    The function value is the adjusted
C                                    derivative.
C
C     Notes:
C     ------
C
C     (1)  IMPLICIT NONE is non-standard.
C
C     Author:  Robert Kennelly, Sterling Federal Systems/NASA-Ames
C     -------
C
C     Development history:
C     --------------------
C
C     18 Feb. 1987    RAK    Initial design and coding.
C
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     &   ONE
      PARAMETER
     &  (ONE = 1.0E+0)

C     Arguments.

      INTEGER
     &   J
      REAL
     &   H (-1:1), DEL (-1:1), BESSEL

C     Local variables.

      REAL
     &   WEIGHT

C     Execution.
C     ----------

C     Estimate first derivative on left (J = 0) or right side (J = 1) of
C     an interval.

      WEIGHT = H (J) / (H (J) + H (J - 1))
      BESSEL = WEIGHT * DEL (J - 1) + (ONE - WEIGHT) * DEL (J)

C     Termination.
C     ------------

      RETURN
      END


C+----------------------------------------------------------------------
C
      FUNCTION BRODLIE (J, H, DEL)
C
C     One-liner: First derivative, adjusted for monotonicity
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        BRODLIE is intended to be used by PLSFIT for determining end
C     conditions on an interval for monotonic interpolation by piecewise
C     cubics. The data must be in the form of arrays containing finite
C     difference interval lengths and 2-point forward difference deriva-
C     tives. See the PLSFIT header for more details.
C
C        The method is due to Brodlie, Butland, Carlson, and Fritsch,
C     as referenced in the PLSFIT header.
C
C     Arguments:
C     ----------
C
C     Name    Type/Dimension  I/O/S  Description
C     J       I               I      Indicates at which end of the
C                                    interval the derivative is to be
C                                    estimated. J = 0 means left-hand
C                                    side, J = 1 means right.
C
C     H       R (-1:1)        I      Array of interval lengths. The 0th
C                                    element is the length of the interval
C                                    on which the cubic is to be deter-
C                                    mined.
C
C     DEL     R (-1:1)        I      Array of derivative estimates. The
C                                    0th element is the forward difference
C                                    derivative over the interval on which
C                                    the cubic is to be determined.
C
C     BRODLIE R                 O    The function value is the adjusted
C                                    derivative.
C
C     Notes:
C     ------
C
C     (1)  IMPLICIT NONE and 8-character symbolic names are non-standard.
C
C     Author:  Robert Kennelly, Sterling Federal Systems/NASA-Ames
C     -------
C
C     Development history:
C     --------------------
C
C     18 Feb. 1987    RAK    Initial design and coding.
C     04 Dec. 2002    DAS    SIGN work-around for Intel compiler bug.
C
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     &   ZERO, ONE, THIRD
      PARAMETER
     &  (ZERO   = 0.0E+0,
     &   ONE    = 1.0E+0,
     &   THIRD  = ONE / 3.0E+0)

C     Arguments.

      INTEGER
     &   J
      REAL
     &   BRODLIE, H (-1:1), DEL (-1:1)

C     Local variables.

      REAL
     &   ALPHA, PRODUCT

C     Execution.
C     ----------

C     Compare the algebraic signs of the two DEL's.  Have to test that
C     at least one is positive to avoid a zero denominator (this fancy
C     test permits one term to be zero, but the answer below is zero
C     anyway in these cases).  The trick is to work around the SIGN
C     function, which returns positive even if its 2nd argument is zero.

C**** NO:  SIGN misbehaves on Intel systems when the 2nd argument is zero.

      PRODUCT = DEL (J - 1) * DEL (J)

      IF (PRODUCT == ZERO) THEN

         BRODLIE = ZERO

      ELSE IF (SIGN (ONE, -DEL (J - 1)) .NE. SIGN (ONE, DEL (J))) THEN

C        Form "weighted harmonic mean" of the two finite-difference
C        derivative approximations.  Note that we try to avoid overflow
C        by not multiplying them together directly.

         ALPHA   = THIRD * (ONE + H (J) / (H (J - 1) + H (J)))
         BRODLIE = PRODUCT / (ALPHA * DEL (J) +
     &      (ONE - ALPHA) * DEL (J - 1))
      ELSE

C        The signs differ, so make this point a local extremum.

         BRODLIE = ZERO
      END IF

C     Termination.
C     ------------

      RETURN
      END


C+----------------------------------------------------------------------
C
      FUNCTION BUTLAND (J, H, DEL)
C
C     One-liner: First derivative, non-central 3-point formula, adjusted
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        Computes a first derivative approximation for PLSFIT over an
C     interval at a data boundary, using a modified forward or backward
C     3-point formula.  The data must be in the form of arrays containing
C     finite difference interval lengths and 2-point forward difference
C     derivatives, and the differencing direction is controlled by a flag.
C     See PLSFIT for more details, or THREEPT for the pure 3-pt. formula.
C
C        The "shape preserving adjustments" are from PCHIP, a monotone
C     piecewise cubic interpolation package by F. N. Fritsch.
C
C     Arguments:
C     ----------
C
C     Name    Type/Dimension  I/O/S  Description
C     J       I               I      Indicates at which end of the
C                                    interval the derivative is to be
C                                    estimated. J = 0 means left-hand
C                                    side, J = 1 means right.
C
C     H       R (-1:1)        I      Array of interval lengths. The 0th
C                                    element is the length of the interval
C                                    on which the cubic is to be deter-
C                                    mined.
C
C     DEL     R (-1:1)        I      Array of derivative estimates. The
C                                    0th element is the forward difference
C                                    derivative over the interval on which
C                                    the cubic is to be determined.
C
C     BUTLAND R                 O    The function value is the adjusted
C                                    derivative.
C
C     Environment:  Fortran 90
C     ------------
C
C     Author:  Robert Kennelly, Sterling Federal Systems/NASA-Ames
C     -------
C
C     History:
C     --------
C
C     18 Feb. 1987    RAK    Initial design and coding, as THREEPT.
C     20 June 1991    DAS    Monotonic form renamed BUTLAND; THREEPT
C                            is now the pure 3-point formula.
C     04 Dec. 2002     "     SIGN work-arounds for Intel compiler bug.
C
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     &   J
      REAL
     &   H (-1:1), DEL (-1:1), BUTLAND

C     Local constants.

      REAL, PARAMETER ::
     &   ZERO  = 0.0E+0,
     &   ONE   = 1.0E+0,
     &   THREE = 3.0E+0

C     Local variables.

      INTEGER
     &   STEP
      REAL
     &   DMAX, WEIGHT
      LOGICAL
     &   CONSTRAIN

C     Execution.
C     ----------

C     Estimate first derivative on a left-hand boundary using a 3-point
C     forward difference (STEP = +1), or with a backward difference for
C     the right-hand boundary (STEP = -1).

      STEP = 1 - J - J   ! J here is consistent with related modules.

C     In {H, DEL} form, the derivative looks like a weighted average.

C***  Avoid zero as the second argument of SIGN: Intel compiler misbehaves.

      IF (DEL (0) == ZERO) THEN

         BUTLAND = ZERO

      ELSE ! BUTLAND below cannot be zero

         WEIGHT  = -H (0) / (H (0) + H (STEP))
         BUTLAND = WEIGHT * DEL (STEP) + (ONE - WEIGHT) * DEL (0)

C        Shape-preserving adjustments.  Note that we try to avoid overflow
C        by not multiplying quantities directly.

         IF (SIGN (ONE, BUTLAND) /= SIGN (ONE, DEL (0))) THEN

C           Defer to the estimate closest to the boundary.

            BUTLAND = ZERO

C******  ELSE IF (SIGN (ONE, DEL (0)) .NE. SIGN (ONE, DEL (STEP))) THEN
         ELSE

            IF (DEL (STEP) == ZERO) THEN
               CONSTRAIN = DEL (0) < ZERO
            ELSE
               CONSTRAIN = SIGN (ONE, DEL (0)) /= SIGN (ONE, DEL (STEP))
            END IF

            IF (CONSTRAIN) THEN

C              If the monotonicity switches, may need to bound the estimate.

               DMAX = THREE * DEL (0)
               IF (ABS (BUTLAND) > ABS (DMAX)) BUTLAND = DMAX
            END IF

         END IF

      END IF

C     Termination.
C     ------------

      END


C+----------------------------------------------------------------------
C
      FUNCTION THREEPT (J, H, DEL)
C
C     One-liner: First derivative, non-central 3-point formula
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        Computes a first derivative approximation for PLSFIT over an
C     interval at a data boundary, using a forward or backward 3-point
C     formula.  The data must be in the form of arrays containing finite
C     difference interval lengths and 2-point forward difference deriva-
C     tives, and the differencing direction is controlled by a flag. See
C     PLSFIT for more details.
C
C        See module BUTLAND for a version with "shape-preserving"
C     adjustments.
C
C     Arguments:
C     ----------
C
C     Name    Type/Dimension  I/O/S  Description
C     J       I               I      Indicates at which end of the
C                                    interval the derivative is to be
C                                    estimated. J = 0 means left-hand
C                                    side, J = 1 means right. 
C
C     H       R (-1:1)        I      Array of interval lengths. The 0th
C                                    element is the length of the interval
C                                    on which the cubic is to be deter-
C                                    mined.
C
C     DEL     R (-1:1)        I      Array of derivative estimates. The
C                                    0th element is the forward difference
C                                    derivative over the interval on which
C                                    the cubic is to be determined.
C                                     
C     THREEPT R                 O    The function value is the derivative.
C
C     Environment:  VAX/VMS; FORTRAN 77
C     ------------
C
C     IMPLICIT NONE and 8-character symbolic names are non-standard.
C
C     Author:  Robert Kennelly, Sterling Federal Systems/NASA-Ames
C     -------
C
C     History:
C     --------
C
C     18 Feb. 1987    RAK    Initial design and coding.
C     06 June 1991    DAS    Original THREEPT renamed BUTLAND; THREEPT
C                            now gives unmodified 1-sided 3-pt. results.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     &   J
      REAL
     &   H (-1:1), DEL (-1:1), THREEPT

C     Local constants.

      REAL
     &   ONE
      PARAMETER
     &  (ONE = 1.0E+0)

C     Local variables.

      INTEGER
     &   STEP
      REAL
     &   WEIGHT

C     Execution.
C     ----------

C     Estimate first derivative on a left-hand boundary using a 3-point
C     forward difference (STEP = +1), or with a backward difference for
C     the right-hand boundary (STEP = -1).

      STEP = 1 - J - J   ! J here is consistent with related modules.

C     In {H, DEL} form, the derivative looks like a weighted average.

      WEIGHT  = -H (0) / (H (0) + H (STEP))
      THREEPT = WEIGHT * DEL (STEP) + (ONE - WEIGHT) * DEL (0)

C     Termination.
C     ------------

      RETURN
      END


C+----------------------------------------------------------------------
C
      SUBROUTINE INTERVAL (NX, X, XFIND, ARROW, LEFT)
C
C     One-liner: Interpolation search for interval containing a point.
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        Written primarily for interval-based interpolations such as
C     piecewise linear or cubic spline, INTERVAL performs a search to
C     locate the best interval for evaluating the interpolant at a
C     given point. The normal case returns the "left-hand" endpoint of
C     the interval bracketing the point, but for the out-of-range cases
C     below or above the range of the knots, the interval to be used is
C     the first or last. The array of knots must be monotonic, either
C     increasing or decreasing. Diagrammatically, LEFT is returned as
C     shown below for the normal case (no extrapolation):
C
C          X (1)  ...   X (LEFT)   X (LEFT+1)   ...      X (NX)
C                                ^
C                              XFIND
C
C     And for extrapolation:
C
C                     X (LEFT = 1)  ...   X (NX)
C             ^
C           XFIND
C
C     or,
C                X (1)  ...   X (LEFT = NX-1)    X (NX)
C                                                           ^
C                                                         XFIND
C
C     If the point to be bracketed (XFIND) matches one of the knots, the
C     index of that knot is returned as LEFT, i.e., the condition for a
C     bracket of an interior point is:
C
C        X (LEFT) <= XFIND < X (LEFT+1)  if  ARROW = +1.0,  or
C        X (LEFT) >= XFIND > X (LEFT+1)  if  ARROW = -1.0.
C
C        This is a low-level routine with minimal error checking. The
C     calling program is assumed to have verified the following:
C
C     (1)  NX >= 2
C     (2)  X strictly monotonic
C     (3)  ARROW = +1.0 or -1.0
C
C     Subroutine PROTECT is available from the author for easily checking
C     conditions (2) and (3). LEFT is verified on input, but efficiency in
C     loops will benefit from passing the best estimate available, usually
C     just the result of the last call.
C
C        INTERVAL was originally written for use with CSEVAL and TABLE1.
C     The interpolation search was adapted from ideas in Sedgewick's book
C     referenced below.
C
C     Arguments:
C     ----------
C
C     Name  Dimension  Type  I/O/S  Description
C     NX                I    I      Number of points in array X; must
C                                   be >= 2 (no check performed).
C
C     X        NX       R    I      Array of points defining the set
C                                   of intervals to be examined. Only
C                                   the first NX-1 points are required.
C
C     XFIND             R    I      The point for which a bracketing
C                                   interval is sought.
C
C     ARROW             R    I      Monotonicity indicator for input
C                                   array X:
C                                     -1.0  strictly decreasing
C                                      0.0  NOT ALLOWED!
C                                     +1.0  strictly increasing
C                                   Supplied by the calling routine for
C                                   reasons of speed (not checked).
C
C     LEFT              I    I/O    Input: guessed index of left-hand
C                                   endpoint of the interval containing
C                                   the specified point.
C
C                                   Output: index of the largest array
C                                   value <= specified point (if ARROW=+1.0).
C                                   Special case for data out of range:
C                                   return left endpoint of closest interval.
C                                   Thus, LEFT = 1 for XFIND < X (2), and
C                                   LEFT = NX-1 for XFIND >= X (NX-1).
C                                   (If ARROW=-1.0, reverse the inequalities.)
C
C     Environment:  Digital VAX-11/780, VMS FORTRAN
C     ------------  Apple Macintosh, Absoft MacFORTRAN/020 v2.3
C
C     Notes:
C     ------
C
C     (1)  IMPLICIT NONE and eight character symbols are not (yet) standard.
C
C     (2)  In speed-critical applications, it might be a good idea to build
C          this algorithm in-line since it is typically called many times
C          from within a loop. Another potential speed-up is removal of the
C          ARROW multiplies, which restricts the method to increasing data.
C          So far, the simplicity of separating out the messy search details
C          and the generality of bi-directional searching have outweighed
C          the modest speed penalty incurred.
C
C     Bibliography:
C     -------------
C
C     (1) Sedgewick, R.  Algorithms.  Reading: Addison-Wesley, 1983.
C            (Chap. 14)
C
C     Author:  Robert Kennelly and David Saunders, Sterling Federal Systems
C     -------
C
C     Development history:
C     --------------------
C
C     20 Oct. 1987    RAK    Interpolation search adapted (with mods.
C                            for bidirectional search and some minor
C                            repair) from CSEVAL (RAK) and TABLE1 (DAS).
C     08 Aug. 1988    DAS    Clarified descriptions of bracketing, where
C                            the inequalities depend upon ARROW.
C
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     &   ONE
      PARAMETER
     &  (ONE = 1.0E+0)

C     Arguments.

      INTEGER
     &   LEFT, NX
      REAL
     &   ARROW, X (NX), XFIND

C     Local variables.

      INTEGER
     &   LENGTH, NXLESS1, RIGHT, TRIAL
      REAL
     &   XBYARROW

C     Execution.
C     ----------

      XBYARROW = XFIND * ARROW

C     Simplify things by disposing of two important special cases so that
C     X (LEFT) and X (RIGHT) can really bracket XFIND. As a by-product,
C     this also takes care of the NX = 2, 3 cases.

      NXLESS1 = NX - 1

      IF (XBYARROW .GE. X (NXLESS1) * ARROW) THEN
         LEFT = NXLESS1
         GO TO 990
      ELSE IF (XBYARROW .LT. X (2) * ARROW) THEN
         LEFT = 1
         GO TO 990
      END IF

C      ---------------------------------
C     |                                 |
C     |   X (2) <= XFIND < X (NX - 1)   |
C     |            - or -               |
C     |   X (2) > XFIND >= X (NX - 1)   |
C     |                                 |
C     |   NX > 3                        |
C     |                                 |
C      ---------------------------------

C     Adjust the pointers. We hope that the calling routine has provided
C     a reasonable guess (since it's probably working on an ordered array
C     of points to evaluate), but check anyway.

      LEFT = MIN (MAX (2, LEFT), NX - 2)

      IF (XBYARROW .GE. X (LEFT) * ARROW) THEN
         IF (XBYARROW .LT. X (LEFT + 1) * ARROW) THEN

C           XFIND is in the original guessed-at interval.

            GO TO 990
         ELSE

C           We'll look farther to the right. Evidently LEFT was < NX - 2.

            RIGHT = NXLESS1
            LEFT  = LEFT + 1
         END IF
      ELSE

C        Look to the left of the guess. Evidently LEFT was > 2.

         RIGHT = LEFT
         LEFT  = 2
      END IF

C      ----------------------------------
C     |                                  |
C     |   2 <= LEFT < RIGHT <= NX - 1    |
C     |                                  |
C      ----------------------------------

C     The interval length must decrease each time through - terminate
C     when the correct interval is found or when the interval length
C     cannot be decreased.

   10 CONTINUE
         LENGTH = RIGHT - LEFT
         IF (LENGTH .GT. 1) THEN

C           The trial value is a "linear" estimate of the left-hand endpoint
C           of the interval bracketing the target XFIND, with protection
C           against round-off (which can affect convergence).

            TRIAL = MIN (RIGHT - 1, LEFT + MAX (0, INT (REAL (LENGTH) *
     &         (XFIND - X (LEFT)) / (X (RIGHT) - X (LEFT)))))

C            ------------------------------------------
C           |                                          |
C           |   2 <= LEFT <= TRIAL < RIGHT <= NX - 1   |
C           |                                          |
C            ------------------------------------------

C           Adjust pointers. Increase LEFT or decrease RIGHT until done.

            IF (XBYARROW .GE. X (TRIAL + 1) * ARROW) THEN
               LEFT  = TRIAL + 1
            ELSE IF (XBYARROW .LT. X (TRIAL) * ARROW) THEN
               RIGHT = TRIAL
            ELSE

C              We're done: XFIND is in the interval [X (TRIAL), X (TRIAL+1)).

               LEFT  = TRIAL
               GO TO 990
            END IF
            GO TO 10

         END IF

C     Termination.
C     ------------

  990 CONTINUE
      RETURN
      END


C+-----------------------------------------------------------------------------+
C
      subroutine HSORTRI ( rlist, ilist, Length )
C
C Title: Sort Real list and corresponding integer list
C
C Purpose:
C   Sort a list (rlist) of real values, and rearrange a second integer list
C   so it matches the first.  If ilist is initialized to 1, 2, 3, etc. before
C   this module is called, consecutive elements of ilist will identify
C   monotonically increasing elements of rlist on return.  This is useful
C   for sorting one or more vectors with regard to another.
C
C Method:
C   Use the HeapSort algorithm of J.W.J. Williams to sort rlist, and change
C   ilist whenever rlist is changed.  Since HeapSort is a sort-in-place
C   algorithm, no scratch storage is needed, and the original order of the
C   arrays is lost.
C
C References:
C   Press, et al. "Numerical Recipes", chapter 8.
C   Try Knuth, "Searching and Sorting" for a full discussion.
C 
C Error Handling: 
C   Length <=0 will probably cause a Fortran array dimensions error, which 
C   cannot be caught here, so the calling application should make sure
C   this doesn't happen.
C
C Notes: This is a simple modification of the HSORTRR subroutine to define
C   a sort permutation index array to be applied to any number of associated
C   vectors that are to be sorted in terms of one vector.
C
C Standards violations: (any or all of the following)
C   >lower case in code
C   >long variable names
C   >underbar "_" used in variable names
C   >"!" used for end-of-line comments
C
C Coding conventions:
C   UPPERCASE names are external references (both code and data)
C   Capitalized names are constants (parameter, read-only arguments, etc.)
C   lowercase names are local variables, read-write arguments, Fortran keywords
C   ">"s prefix normal comments, ">" for high level comments, ">>" for details
C   { } is used for assertions (statements assumed to be true when they appear)
C   C#  comments for noteworthy code (modifications, machine dependent
C         code, Fortran extensions, debugging code, etc)
C   Labels below 999 are for normal code (999 is always for STOP or RETURN)
C   Labels from 1000 to 1999 are format statements
C   Labels from 2000 on      are for error handling
C   Syntax for code modifications:
C#mod0.0new# -- comment about modification 
C     ...new code
C#mod0.0old# 
C#    ...replaced code
C#end0.0#
C
C Environment:  FORTRAN 77 + extensions
C
C Author:  David Serafini, Sterling Software, Palo Alto, CA
C 
C Modification history:
C   ver 0.0  Jul 87  D.B.Serafini  Initial design and coding
C            Jul 88  P.J.Trosin    Built "ilist=integer" version
C            Jul 91  D.B.Serafini  Added check for Length=1
C            Dec 91  D.A.Saunders  Use on a Macintosh by Roy Hampton
C                                  encountered failure.  Changed a
C                                  composite IF test into two IFs.
C            Oct 99  D.L.Hermstad  Macintosh/Absoft version honors case
C                                  of variables.  Therefore, changed
C                                  argument length to Length.  (Prior
C                                  version should have failed to compile
C                                  "real    rlist(Length)"!)
C+---------------------------------------------------------------------+

C     >declare arguments
      integer Length            !length of rlist, ilist
      real    rlist(Length)     !list of real values to sort, used as sort key
      integer ilist(Length)     !list of (integer) values, not used for key

C-----------------------------------------------------------------------

C     >declare local variables
      integer lu ,ld            !pointers to the left and right boundaries of
     &       ,ru ,rd            ! the heap during sift-up and sift-down passes
      real    rtemp             !temporaries used to swap elements 
      integer itemp             !in the lists


C     >begin execution

C     >handle special case of 1 element list
      if( Length .EQ. 1 ) goto 999

C     >initialize
      lu = ( Length / 2 ) +1
      ru = Length

C     >loop until whole array is sorted, i.e. both down & up sifts are complete
 10   continue
C         >if lu>1, still sifting up
          if( lu .GT. 1 )then
C             >>sift up phase: move the left boundary down, and set the temp
              lu = lu -1
              rtemp = rlist(lu)
              itemp = ilist(lu)

          else
C             >>sift down: move the min in the heap to the low end of the
C               sorted array
              rtemp = rlist(ru)
              itemp = ilist(ru)
              rlist(ru) = rlist(1)
              ilist(ru) = ilist(1)
              ru = ru -1
C             >>see if we're done
              if( ru .EQ. 1 )then
                  rlist(1) = rtemp
                  ilist(1) = itemp
                  goto 999   !return
              endif

          endif

C         >>set up for sift-down
          ld = lu
          rd = lu + lu

C         >>loop while right end of sorted region is before end of heap
 20       if( rd .LE. ru )then

C             >>if not at the end of the heap, compare this element to the
C               next higher element, and expand the sorted region by 1 if its
C               less
C****         if( rd .LT. ru  .AND.
C****&            rlist(rd) .LT. rlist(rd+1) ) rd = rd +1
              if( rd .LT. ru )then
                  if( rlist(rd) .LT. rlist(rd+1) ) rd = rd +1
              endif

C             >>sift down until the right element is sorted
              if( rlist(rd) .GE. rtemp )then
                  rlist(ld) = rlist(rd)
                  ilist(ld) = ilist(rd)
                  ld = rd
                  rd = rd + rd

              else
                  rd = ru +1    !indicates end of sift-down

              endif

C             >>endloop 20
              goto 20

          endif      !{rd .GT. ru}

C         >>save the temp values
          rlist(ld) = rtemp
          ilist(ld) = itemp

C     >>endloop 10 - keep looping until the smallest value has sifted to the top
      goto 10


C     >terminate execution
 999  return

      end
