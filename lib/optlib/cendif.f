*DECK CENDIF
C+----------------------------------------------------------------------
C
      SUBROUTINE CENDIF (N, X, FX, EPSOBJ, H, GRAD, DIAG, NFUNCT,
     >   LWRITE, FUNCT)
C
C
C     Description:
C
C           CENDIF calculates central difference approximations to the
C        gradient and Hessian of an optimization objective function.  If
C        requested, improved finite difference stepsizes are calculated.
C        CENDIF/FDSTEP/DIFFER were developed for use with QNMDIF (a
C        quasi-Newton optimizer), but the stepsizes may be used by any
C        forward-difference method.
C
C
C     Parameters:
C
C        Name   Dimension   Type  I/O/S   Description
C        N                   I    I       Number of optimization variables.
C        X       N           R    I       Optimization variables.
C        FX                  R    I       Initial function value.
C        EPSOBJ              R    I       Estimate of absolute error in
C                                         objective.
C        H       N           R    I/O     Initial stepsizes used for
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
C        LWRITE              I      I     Logical unit number for output
C                                         diagnostics, provided >= 0.
C                                         Use a negative unit number to
C                                         suppress all output.
C        FUNCT                            Objective function subroutine
C                                         (EXTERNAL).
C
C     External References:
C
C        Name     Description
C        DIFFER   Computes first and second derivative approximations
C                 by central differences.
C        FDSTEP   Produces estimates of optimal one sided finite
C                 difference stepsizes.
C
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77).
C
C
C     Author:  Robert Kennelly, Sterling Software.
C
C
C     Development History:
C
C        16 July 1982    RAK    Complete revision of earlier version by
C                               Trosin - now use FDSTEP and DIFFER.
C         2 Sep. 1982    RAK    Reordered parameters, redefined NFUNCT.
C        19 Jan. 1983    RAK    Updated NFDEL following DIFFER.
C        10 Mar. 1983    RAK    Repaired call to DIFFER (... DIAG(I) ...),
C                               added IF-THEN-ELSE blocks.
C         9 Aug. 1986    RAK    Added LWRITE input parameter (had been
C                               in COMMON); reordered parameters again
C                               (sorry!).  If LWRITE is negative,
C                               suppress all output.  Use PARAMETER
C                               for local array dimensions, and test
C                               input N for consistency.  Increase
C                               MAXN from 30 to 40.  Add IMPLICIT NONE.
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      INTEGER
     >   MAXN
      PARAMETER
     >   (MAXN = 40)

      REAL
     >   ZERO, HALF
      PARAMETER
     >   (ZERO = 0.0E+0,
     >    HALF = 0.5E+0)

C     Arguments.

      INTEGER
     >   LWRITE, N, NFUNCT
      REAL
     >   DIAG (N), FX, EPSOBJ, GRAD (N), H (N), X (N)

C     Local variables.

      INTEGER
     >   I, IERROR (MAXN), NFDEL (MAXN)
      REAL
     >   CFPBAC, CFPFOR, CFPP, ERRBND, FPBAC, FPFOR, HI, URDIAG,
     >   URH (MAXN)

C     Procedures.

      EXTERNAL
     >   DIFFER, FDSTEP, FUNCT


C     Execution.
C     ----------

      IF (LWRITE .GE. 0) WRITE (LWRITE, 1000)
      IF (N .GT. MAXN) STOP 'CENDIF:  ERROR - too many variables!'

C     Loop over components of H and compute derivatives or improved
C     stepsizes (as well as first and second derivatives).

      NFUNCT = 0
      DO 10 I = 1, N
         HI = H (I)
         URH (I) = HI
         IF (HI .GT. ZERO) THEN

C           Use the stepsize passed in for differencing.

            CALL DIFFER (I, N, X, EPSOBJ, HI, FX, FPFOR, CFPFOR,
     >         FPBAC, CFPBAC, DIAG (I), CFPP, FUNCT)
            IERROR (I) = 0
            NFDEL (I) = 2
            GRAD (I) = HALF * (FPFOR + FPBAC)
         ELSE

C           Generate H from scratch, with GRAD and DIAG as by-products.

            CALL FDSTEP (I, N, X, EPSOBJ, FX, GRAD (I), DIAG (I), H (I),
     >         ERRBND, IERROR (I), NFDEL (I), FUNCT)
         END IF
         NFUNCT = NFUNCT + NFDEL (I)

C        Safeguard Hessian approx. by bounding diagonals away from zero
C        (to preserve positive definitness).

         IF (DIAG (I) .LT. EPSOBJ) THEN
            URDIAG = DIAG (I)
            DIAG (I) = MAX (ABS (DIAG (I)), EPSOBJ)
            IF (LWRITE .GE. 0) WRITE (LWRITE, 1010) I, URDIAG, DIAG (I)
         END IF
   10 CONTINUE

C     Recap results.

      IF (LWRITE .GE. 0) THEN
         WRITE (LWRITE, 1020)
         WRITE (LWRITE, 1030) (I, X (I), GRAD (I), DIAG (I), URH (I),
     >      H (I), NFDEL (I), IERROR (I), I = 1, N)
         WRITE (LWRITE, 1040)
      END IF


C     Termination.
C     ------------

      RETURN


C     Formats.
C     --------

 1000 FORMAT (37H0CENDIF:  Begin gradient calculation.)
 1010 FORMAT (25H0CENDIF:  WARNING - DIAG(, I3, 15H) changed from ,
     >   E12.6, 4H to , E12.6)
 1020 FORMAT (8H1CENDIF:/ 1H0, 4X, 1HI, 4X, 1HV, 20X, 8HGradient,
     >   13X, 15HDiag of Hessian, 6X, 5HOld H, 16X, 5HNew H, 14X,
     >   6HNo Fun, 2X, 6HIERROR)
 1030 FORMAT (1H , 3X, I2, 3X, E18.12, 3X, E18.12, 3X, E18.12, 3X,
     >   E18.12, 3X, E18.12, 4X, I2, 6X, I1)
 1040 FORMAT (1H1)

      END
*DECK FDSTEP
C+----------------------------------------------------------------------
C
      SUBROUTINE FDSTEP (I, N, X, EPSOBJ, FX, FP, FPP, HOPT, ERRBND,
     >   IERROR, NFUNCT, FUNCT)
C
C
C     Description:
C
C           This routine is intended for automatic estimation of optimal
C        finite difference intervals for numerical optimization of functions
C        whose derivatives cannot be calculated analytically.  Each call to
C        FDSTEP will handle one component of the stepsize vector.
C
C
C     Parameters:
C
C        Name    Dimension   Type  I/O/S  Description
C        I                    I    I      Index of the component of X to be
C                                         varied.
C        N                    I    I      Dimension of X.
C        X          N         R    I      Vector of optimization variables.
C        EPSOBJ               R    I      Smallest significant absolute change
C                                         in the objective function FUNCT.
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
C                                         stepsize for I-th component of X.
C        ERRBND               R      O    Estimated bound on (absolute) error
C                                         in first derivatives calculated
C                                         using one sided differencing with
C                                         stepsize HOPT.
C        IERROR               I      O    Error status upon termination. Zero
C                                         means OK; for other values see
C                                         error handling section below.
C        NFUNCT               I      O    Number of calls to FUNCT (increment).
C        FUNCT                     I      Routine for calculating objective
C                                         function (EXTERNAL).
C
C
C     External References:
C
C        DIFFER     Utility for calculating derivatives by finite differences,
C                   along with cancellation error estimates.
C        FUNCT      User-supplied objective function subroutine.
C
C
C     Notes:
C
C        (1)  The goal is to find the smallest stepsize H(KPHI) which yields
C             acceptable relative cancellation errors in the forward and
C             backward difference first derivatives and the central difference
C             second derivative of FUNCT.  The second derivative computed with
C             H(KPHI) is then used to obtain HOPT, an interval which should
C             result in good forward difference approximations to the first
C             derivative.  The "optimality" of the interval chosen is based
C             on several assumptions - consult reference (1) for details.
C
C        (2)  The algorithm also returns estimates FP, the first derivative,
C             and FPP, the second derivative, which are both computed by
C             central differences.
C
C        (3)  The stepsize and derivatives calculated each iteration are
C             saved in arrays for possible later use to save function
C             evaluations.  These are relevant to the step associated with
C             optimization variable I only, and should not be confused with
C             the stepsize and gradient arrays in the calling routine.
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
C                pp. 341-344.  London:  Academic Press, 1981.
C
C        (2)  Gill, P., Murray, W., Saunders, M., and Wright, M.  "A Procedure
C                for Computing Forward Difference Intervals for Numerical
C                Optimization."  Tech. Rep. SOL 81-25, Dept. of Operations
C                Research, Stanford University, Dec. 1981.
C
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77).
C
C
C     Author:  Robert Kennelly, Informatics General Corporation.
C
C
C     Development History:
C
C        16 June 1982   RAK   Original coding ("literal" transcription).
C        23 June 1982   RAK   Revised and extensively restructured.
C        13 July 1982   RAK   Section FD3 revised (small CFPP now OK).
C         2 Sep. 1982   RAK   Redefined function counter NFUNCT.
C        24 Sep. 1982   RAK   Some minor changes adapted from Ref. (2).
C        11 Aug. 1986   RAK   Add IMPLICIT NONE.
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      INTEGER
     >   KBOUND
      PARAMETER
     >   (KBOUND = 6)

      REAL
     >   ZERO, ONE, TWO, HALF, BNDFP, BNDLO, BNDUP, ETA, OMEGA, RHO,
     >   TINY
      PARAMETER
     >   (ZERO  = 0.0E+0,
     >    ONE   = 1.0E+0,
     >    TWO   = 2.0E+0,
     >    HALF  = ONE / TWO,
     >    BNDFP = 1.0E-1,
     >    BNDLO = 1.0E-3,
     >    BNDUP = 1.0E-1, 
     >    ETA   = ONE,
     >    OMEGA = ONE,
     >    RHO   = 1.0E+1,
     >    TINY  = 1.0E-32)

C     Arguments.

      INTEGER
     >   I, IERROR, N, NFUNCT
      REAL
     >   EPSOBJ, ERRBND, FP, FPP, FX, HOPT, X (N)

C     Local variables.

      LOGICAL
     >   DONE
      INTEGER
     >   KOUNT, KPHI, KSAVE
      REAL
     >   CFPP, CFPBAC, CFPFOR, H (0:KBOUND), HBAR, FPFOR (0:KBOUND),
     >   FPBAC (0:KBOUND), FPPCEN (0:KBOUND)

C     Procedures.

      EXTERNAL
     >   DIFFER, FUNCT

C     Statement functions.

      LOGICAL
     >   INSIDE, LARGE, MAXLE, SMALL
      REAL
     >   A, B

      INSIDE (A)   = (A .GE. BNDLO) .AND. (A .LE. BNDUP)
      LARGE (A)    = A .GT. BNDUP
      MAXLE (A, B) = MAX (A, B) .LE. BNDFP
      SMALL (A)    = A .LT. BNDLO


C     Execution.
C     ----------

C     Initialization (FD1).

      KOUNT = 0
      NFUNCT = 0
      IERROR = 0
      DONE = .FALSE.

      HBAR = TWO * (ETA + ABS (X (I))) *
     >   SQRT (EPSOBJ / (OMEGA + ABS (FX)))
      H (KOUNT) = HBAR * RHO

      CALL DIFFER (I, N, X, EPSOBJ, H (KOUNT), FX, FPFOR (KOUNT),
     >   CFPFOR, FPBAC (KOUNT), CFPBAC, FPPCEN (KOUNT), CFPP, FUNCT)
      NFUNCT = NFUNCT + 2

C     Accept, increase, or decrease H?  (FD2).

      IF (MAXLE (CFPFOR, CFPBAC) .AND. .NOT.LARGE (CFPP)) THEN
         IF (INSIDE (CFPP)) THEN

C           Accept H and fall through to estimate of optimal H.

            KPHI = KOUNT
         ELSE

C           Decrease H to reduce truncation error (FD4).

C           The cancellation errors in the first derivatives are OK, while
C           that for the second derivative is smaller than necessary.
C           Try to reduce H without letting the errors get too big.

   10       CONTINUE
               KOUNT = KOUNT + 1
               H (KOUNT) = H (KOUNT - 1) / RHO
               CALL DIFFER (I, N, X, EPSOBJ, H (KOUNT), FX,
     >            FPFOR (KOUNT), CFPFOR, FPBAC (KOUNT), CFPBAC,
     >            FPPCEN (KOUNT), CFPP, FUNCT)
               NFUNCT = NFUNCT + 2
               IF (.NOT.MAXLE (CFPFOR, CFPBAC) .OR. LARGE (CFPP)) THEN

C                 We've gone too far - back up one iteration and quit.

                  DONE = .TRUE.
                  KPHI = KOUNT - 1
               ELSE IF (INSIDE (CFPP)) THEN

C                 The current stepsize H is acceptable.

                  DONE = .TRUE.
                  KPHI = KOUNT
               ELSE

C                 Second derivative cancellation error is still smaller
C                 than necessary, but check iteration counter before
C                 continuing.

                  IF (KOUNT .GE. KBOUND) THEN

C                    The iteration limit has been reached.  The error
C                    flag is set as a warning that the stepsize may not
C                    be as small as possible.  Note:  this can also
C                    indicate trouble - a slope discontinuity can also
C                    give this error, so the objective function should be
C                    double checked in this case.

                     IERROR = 1
                     KPHI = KOUNT
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

            IF (MAXLE (CFPFOR, CFPBAC) .AND. KSAVE .LT. 0)
     >         KSAVE = KOUNT
            KOUNT = KOUNT + 1
            H (KOUNT) = H (KOUNT - 1) * RHO
            CALL DIFFER (I, N, X, EPSOBJ, H (KOUNT), FX, FPFOR (KOUNT),
     >         CFPFOR, FPBAC (KOUNT), CFPBAC, FPPCEN (KOUNT), CFPP,
     >         FUNCT)
            NFUNCT = NFUNCT + 2
            IF (MAXLE (CFPFOR, CFPBAC) .AND. .NOT.LARGE (CFPP)) THEN

C              Current H is acceptable.

               DONE = .TRUE.
               KPHI = KOUNT
            ELSE IF (KOUNT .GE. KBOUND) THEN

C              No satisfactory interval has been found (FD6).

               IF (.NOT.MAXLE (CFPFOR, CFPBAC)) THEN

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


C     Termination.
C     ------------

      RETURN
      END
*DECK DIFFER
C+----------------------------------------------------------------------
C
      SUBROUTINE DIFFER (I, N, X, EPSOBJ, H, FX, FPFOR, CFPFOR,
     >   FPBAC, CFPBAC, FPPCEN, CFPP, FUNCT)
C
C
C     Description:
C
C           Utility for calculating derivatives and cancellation error
C        estimates by finite differences.  Called repeatedly by subroutine
C        FDSTEP during estimation of optimal stepsizes to be passed to a
C        nonlinear optimization routine.  See CENDIF/FDSTEP headers for
C        more details.
C
C
C     Parameters:
C
C        Name    Dimension   Type  I/O/S  Description
C        I                    I    I      Index of the component of X to be
C                                         varied.
C        N                    I    I      Dimension of X in calling routine.
C        X          N         R    I      Vector of optimization variables.
C        EPSOBJ               R    I      Smallest significant absolute change
C                                         in the objective function FUNCT.
C        H                    R    I      Finite difference stepsize for I-th
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
C        FUNCT                     I      Routine for calculating objective
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
C     Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77).
C
C
C     Author:  Robert Kennelly, Informatics General Corporation.
C
C
C     Development History:
C
C        23 June 1982    RAK    Original coding.
C         2 Sep. 1982    RAK    Deleted local function counter.
C        11 Aug. 1986    RAK    Add IMPLICIT NONE.
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     >   TWO, FOUR, TINY
      PARAMETER
     >   (TWO  = 2.0E+0,
     >    FOUR = 4.0E+0,
     >    TINY = 1.E-32)

C     Arguments.

      INTEGER
     >   I, N
      REAL
     >   CFPBAC, CFPFOR, CFPP, EPSOBJ, FPBAC, FPFOR, FPPCEN, FX, H,
     >   X (N)

C     Local variables.

      REAL
     >   FXMH, FXPH, XINIT

C     Procedures.

      EXTERNAL
     >   FUNCT


C     Execution.
C     ----------

      XINIT = X (I)

C     Compute first derivative by forward difference...

      X (I) = XINIT + H
      CALL FUNCT (N, X, FXPH)
      FPFOR = (FXPH - FX) / H
      CFPFOR = TWO * EPSOBJ / MAX (H * ABS (FPFOR), TINY)

C     ... and by backward difference.

      X (I) = XINIT - H
      CALL FUNCT (N, X, FXMH)
      FPBAC = (FX - FXMH) / H
      CFPBAC = TWO * EPSOBJ / MAX (H * ABS (FPBAC), TINY)

C     Calculate second derivative by central differences.

      FPPCEN = (FXPH - TWO * FX + FXMH) / MAX (H ** 2, TINY)
      CFPP = FOUR * EPSOBJ / MAX (H ** 2 * ABS (FPPCEN), TINY)

C     Restore X (I) to original value.

      X (I) = XINIT


C     Termination.
C     ------------

      RETURN
      END
