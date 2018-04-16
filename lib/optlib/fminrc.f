C+------------------------------------------------------------------------------
C
      SUBROUTINE FMINRC (AX, BX, XMIN, FMIN, TOL, NUMFUN, CALLER,
     >                   LUNOUT, ISTAT)
C
C     Oneliner:
C
C     Function MINimization (1-D; no derivatives; Reverse Communication)
C     -        ---                                -       -
C
C     Description:
C
C        FMINRC estimates the point where a function of 1 variable has a
C     minimum value on the interval (AX, BX).  The corresponding function
C     value is also returned upon convergence.  Function derivatives are
C     not required.
C
C        This is an adaptation of FMIN77, which is itself a restructured,
C     subroutine version of Brent's FMIN function routine, referenced below.
C     FMINRC avoids the problems of passing as arguments the names of modules
C     with assumed calling sequences (such as FMIN77's FUN (X)) at the expense
C     of having to return to the calling program for every function evaluation. 
C     Such reverse communication is usually preferable to the use of common
C     blocks forced on nontrivial function routines by the likes of FMIN.
C
C
C     Arguments:
C
C     Name   Type   I/O/S   Description
C
C     AX      R     I       Left endpoint of initial interval.
C
C     BX      R     I       Right endpoint of initial interval.
C
C     XMIN    R       O     XMIN is output with the point at which the
C                           next function evaluation is required, until
C                           termination, when XMIN is the best estimate
C                           of the location of the minimum.
C
C     FMIN    R     I/O     After the initial call, FMIN should be input
C                           with the function value corresponding to XMIN.
C                           Upon termination, FMIN is output with the lowest
C                           function value found, corresponding to the
C                           final value of XMIN.
C
C     TOL     R     I       Desired length of the interval of uncertainty
C                           of the final result.  TOL >= 0.
C
C     NUMFUN  I     I/O     At the start of a minimization, NUMFUN should
C                           be input with the maximum number of function
C                           evaluations to be permitted during this
C                           minimization.  On output, NUMFUN contains
C                           the number of function evaluations so far.
C
C     CALLER  C*(*) I       Name of the calling program, printed if LUNOUT > 0
C                           and also used to warn of the (rare) possibility of
C                           nested minimizations, which are not supported
C                           by this version of the FMIN algorithm.  No more
C                           than 12 characters are significant.
C
C     LUNOUT  I     I       Logical unit number for output of the iteration
C                           history.  Use LUNOUT < 0 if you don't want the
C                           evaluations printed.  Any diagnostics go to
C                           the unit defined by ABS (LUNOUT).
C
C     ISTAT   I     I/O     Status flag used as follows:
C                           ISTAT = +2 on input means this call initiates
C                                      a new minimization.
C                           ISTAT = +1 on return means the calling program
C                                      should evaluate the function at the
C                                      current value of XMIN, and call FMINRC
C                                      again with this value in argument FMIN
C                                      and ISTAT = +1 still.
C                           ISTAT =  0 on output means a minimum has been
C                                      found to the specified tolerance.
C                           ISTAT = -1 on output means the limit on the
C                                      number of function evaluations
C                                      has been reached but the convergence
C                                      criteria have not been satisfied.
C                                      XMIN, FMIN are the best values found.
C                           ISTAT = -2 means AX > BX.  This is checked only
C                                      at the start of a minimization.
C                           ISTAT = -3 means the input value of NUMFUN on
C                                      the first call is either non-positive
C                                      or absurdly large.
C                           ISTAT = -4 means the calling program's name
C                                      has changed in the middle of a
C                                      minimization.  Such nested use
C                                      of FMINRC is not supported because
C                                      local variables are reused until a
C                                      minimization is complete.
C
C     FMINRC Usage (1):  FORTRAN 77 style
C
C        ISTAT = 2          ! Initialize the minimization
C        AX = ...
C        BX = ...
C        NUMFUN = ...       ! Max. no. of fn. evals. allowed
C                           <Define TOL and SUBNAME as well.>
C     10 CONTINUE
C
C           CALL FMINRC (AX, BX, XMIN, FMIN, TOL, NUMFUN, SUBNAME,
C       >                LUNOUT, ISTAT)
C
C           IF (ISTAT .LT. -1) THEN      ! Fatal error
C
C              <Handle it>
C
C           ELSE IF (ISTAT .LT. 0) THEN  ! Iteration limit reached
C
C              <Unconverged, but XMIN, FMIN may still be usable>
C
C           ELSE IF (ISTAT .GT. 0) THEN  ! Evaluate the function
C
C              CALL fun (XMIN, FMIN)     ! Or whatever
C              GO TO 10
C
C           ELSE ! ISTAT = 0 (success).  Ensure that FMIN is at XMIN.
C
C              CALL fun (XMIN, FMIN)     ! Or whatever
C
C           END IF
C
C 
C     FMINRC Usage (2):  Fortran 90 style
C
C           This example illustrates least-squares fitting of a mathematical
C           model with one nonlinear coefficient B and two linear coefs. A & C,
C           such as  y = Ae^(Bx) + C.
C
C           call estimateB ()
C
C           B1 = scale1 * B  ! Search interval for optimized B
C           B2 = scale2 * B
C
C        !  A minimum can be found to within sqrt (epsilon), but avoid the sqrt:
C
C           if (epsilon (tol) < 1.e-10) then
C              tol = 1.e-8
C           else
C              tol = 1.e-4
C           end if
C
C           lunerr = abs (lunout)
C           numfun = nfmax      ! Limit; FMINRC typically takes 12-14 iterations
C           istat  = 2          ! Initialize the minimization
C           ier    = 0          ! Upon exit, normally
C
C           do while (istat > 0)
C
C              call fminrc (B1, B2, B, ssq, tol, numfun, caller, lunout, istat)
C
C              if (istat < -1) then
C
C                 write (lunerr, '(/, 2a)') caller, ': FMINRC fatal error'
C                 ier = -1
C
C              else if (istat < 0) then ! Iteration limit; may be usable
C
C                 write (lunerr, '(/, 2a)') caller, ': Iteration limit.'
C                 istat = 0
C
C              else  ! istat >= 0       ! Evaluate the objective function at B
C
C                 call computeAC ()     ! Linear least squares A, C for this B
C                 call sumofsquares ()  ! Sum of squared devs. being minimized
C
C                 if (istat == 0) exit  ! Success
C              end if
C
C           end do
C
C     Note that A and C and ssq are (indeed) updated when istat = 0 for B above
C     to ensure that the BEST coefficients are obtained, not those from the LAST
C     iterate, which typically gives a slightly (negligibly) higher ssq.
C
C     Notes on the FMIN algorithm:
C
C        The method used is a combination of Golden Section search and
C     successive parabolic interpolation.  Convergence is never much slower
C     than that for a Fibonacci search.  If the function has a continuous
C     second derivative which is positive at the minimum (which is not at
C     AX or BX), then convergence is superlinear, and usually of the order
C     of about 1.324.
C
C        The function is never evaluated at two points closer together
C     than RTEPS * ABS (FMIN) + (TOL/3), where RTEPS is approximately the
C     square root of the relative machine precision.  If the function is
C     unimodal and the computed values are always unimodal when separated
C     by at least RTEPS * ABS (X) + (TOL/3), then the FMIN algorithm
C     approximates the abcissa of the global minimum of the function on the
C     interval (AX, BX) with an error less than 3 * RTEPS * ABS (FMIN) + TOL.
C     If the function is not unimodal, then FMIN may approximate a local, but
C     perhaps non-global, minimum to the same accuracy.
C
C        FMIN77 is a modified version of the ALGOL 60 procedure LOCALMIN
C     in Richard Brent, Algorithms for Minimization without Derivatives,
C     Prentice-Hall, Inc. (1973).  The FORTRAN 77 translation started from
C     the FORTRAN 66 version, FMIN, in Forsythe, Malcolm, and Moler but
C     also made use of the original ALGOL.
C
C        FMINRC is a minimal modification which takes advantage of the SAVE
C     feature of FORTRAN 77 to ensure that local variables are preserved
C     between related calls.
C
C
C     Further implementation notes:
C
C        (1)  A more definitive adaptation permitting nested minimizations
C             would require the local variables to be transferred to and
C             from an additional work-space array argument as in ZERORC.
C
C     History:
C
C     12/31/83  R.A.Kennelly  FMIN restructured and (partially) cleaned up
C               NASA Ames     as FMIN77.  Notes on unanswered questions are
C                             sprinkled throughout the code.
C
C     04/10/92  D.A.Saunders  FMINRC adapted from FMIN77 to avoid the FUN
C         Sterling Software/  argument problems.  Added the argument for
C               NASA Ames     counting/limiting the no. of function evals.
C                             RTEPS is now calculated during initialization.
C
C     02/28/93  DAS           Ensuring that FMIN is at XMIN was not shown
C                             in the code example above.
C
C     11/12/10  D.A.Saunders  Added the Fortran 90-style usage example above;
C               ERC, Inc/ARC  made use of the intrinsic for machine epsilon.
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      REAL
     >   AX, BX, XMIN, FMIN, TOL
      INTEGER
     >   NUMFUN, LUNOUT, ISTAT
      CHARACTER
     >   CALLER * (*)

C     Local constants:

C     RATIO is ~(3 - SQRT (5)) / 2, the squared inverse of the Golden Ratio.

      REAL, PARAMETER ::
     >   ZERO = 0.0E+0, ONE = 1.0E+0, TWO = 2.0E+0, THREE = 3.0E+0,
     >   HALF = ONE / TWO, THIRD = ONE / THREE,
     >   RATIO = .381966011250105E+0, SMALL = 1.E+0 / 2048.E+0

C     Local variables:

      REAL
     >   A, B, DELTA, E, F, FU, FV, FW, FX, MIDDLE, P, Q, R, RTEPS,
     >   TOL1, TOL2, U, V, W, X
      INTEGER
     >   LENNAM, LUNERR, MAXFUN
      LOGICAL
     >   GOLDEN, TEST
      CHARACTER
     >   NAME * 12

C     System functions:

      INTRINSIC
     >   ABS, SIGN, SQRT

      SAVE         ! Vital for repeated calls during one minimization


C     Execution:
C     ----------

CRAK
CRAK  How about checking all inputs, e.g. TOL, etc.?
CRAK

C     Notation:  At the start of a cycle, we have
C
C        A, B   Left and right endpoints of current interval of uncertainty
C        X      Point at which the function is smallest, or the most recent
C               point if there was a tie
C        U      Last point at which the function was evaluated
C        W      Second best point thus far
C        V      Third best point thus far
C        DELTA  Proposed offset from latest point (U = X + DELTA)
C        E      A measure of the previous step size (precisely old DELTA when
C               a parabolic step was taken)


      IF (ISTAT .GT. 1) THEN

C        Initialize the iteration.
C        First, compute the square root of the relative machine precision:

CC       RTEPS = SMALL
CC  5    CONTINUE
CC          RTEPS = RTEPS / TWO
CC          TOL1 = ONE + RTEPS
CC          IF (TOL1 .GT. ONE)
CC   >   GO TO 5
CC
CC       RTEPS = SQRT (RTEPS)

         RTEPS = SQRT (EPSILON (RTEPS))  ! Differs by a factor of 2 though!

C        Initialize a new minimization.
C        ------------------------------
               
         LUNERR = ABS (LUNOUT)
         LENNAM = MIN (LEN (CALLER), LEN (NAME))
         NAME = CALLER
         IF (LUNOUT .GT. 0) WRITE (LUNOUT, 1010) NAME

         IF (AX .GT. BX) THEN
            ISTAT = -2
            WRITE (LUNERR, 1030) AX, BX

         ELSE IF (NUMFUN .LE. 1 .OR. NUMFUN .GT. 200) THEN
            ISTAT = -3
            WRITE (LUNERR, 1040) NUMFUN

         ELSE
            ISTAT = 1
            MAXFUN = NUMFUN
            NUMFUN = 0
            A  = AX
            B  = BX
            E  = ZERO
            X  = A + RATIO * (B - A)
            W  = X
            V  = X
            XMIN = X

C           Return early for the first function evaluation.
         END IF

         GO TO 99


      ELSE IF (NUMFUN .EQ. 0) THEN

C        Second call to FMINRC after the first function evaluation.

         NUMFUN = 1
         FX = FMIN
         FW = FX
         FV = FX
         IF (LUNOUT .GT. 0) WRITE (LUNOUT, 1020) NUMFUN, X, FX

         TEST = .TRUE.
      END IF


      IF (CALLER (1 : LENNAM) .NE. NAME (1 : LENNAM)) THEN  ! Fatal error
         ISTAT = -4
         WRITE (LUNERR, 1050) CALLER
         GO TO 99
      END IF


C     Start of a regular iteration.
C     -----------------------------

   10 CONTINUE

      IF (TEST) THEN   ! "TEST" was introduced to allow jumping out for a
C                      ! function value then continuing on the next entry.

C        Test the two stopping criteria.
C        -------------------------------

C        The distance from X to the midpoint of the interval of uncertainty
C        plus half the width of the interval must be less than twice the
C        tolerance.

         MIDDLE = HALF * (A + B)
         TOL1 = RTEPS * ABS (X) + THIRD * TOL
         TOL2 = TWO * TOL1

CRAK     Note Brent doesn't use THIRD above (added in FM&M) ?
CRAK     Equivalent to:  X-A .LE. TOL2  when X > MIDDLE   ...and...
CRAK                     B-X .LE. TOL2  when X < MIDDLE
CRAK     This isn't the same as the error estimate in the header. (It's
CRAK     actually a little tighter.)

         IF (ABS (X - MIDDLE) .LE. (TOL2 - HALF * (B - A))) THEN ! Success
            ISTAT = 0
            GO TO 20
         ELSE IF (NUMFUN .EQ. MAXFUN) THEN ! Too many function evaluations
            ISTAT = -1
            GO TO 20
         END IF


C        Compute a trial step.
C        ---------------------

C        Use parabolic interpolation if possible.

         GOLDEN = (NUMFUN .LT. 3 .OR. ABS (E) .LE. TOL1)
         IF (.NOT.GOLDEN) THEN


C           Fit a parabola.
C           ---------------

            R = (X - W) * (FX - FV)
            Q = (X - V) * (FX - FW)
            P = (X - V) * Q - (X - W) * R
            Q = TWO * (Q - R)

            IF (Q .GT. ZERO) THEN
               P = -P
            ELSE
               Q = -Q
            END IF

C           Is the parabola acceptable ?  Require that the proposed step be:
C
C           (a) Less than half of the second-to-last one (measured
C               roughly by E in case the last step was by Golden Section),
C               *** I don't like this one !!!  No good reason for it that
C               I can see, and not discussed in Brent except a mention that
C               parabolic steps should be small near a minimum, p.74  *** 
C
C           (b) Smaller than the steps to the left or right endpoints of
C               the current interval of uncertainty.  *** Note that (b)
C               also tosses out the parabolic choice on the second step,
C               as must happen since the fit is ill-defined. ***  NUMFUN
C               makes this redundant, though. ***

            IF (ABS (P) .LT. ABS (HALF * Q * E) .AND.
     >         P .GT. Q * (A - X) .AND.
     >         P .LT. Q * (B - X)) THEN


C              Take a parabolic interpolation step.
C              ------------------------------------

C              Here, E is saved as the step taken last iteration.

CRAK           This is suspect, since when parab. follows golden, we end up
CRAK           using E in the first test, then saving corresponding DELTA as
CRAK           next E, thus lagging SAME THING in for TWO iterations !?

               E = DELTA
               DELTA = P / Q
               U = X + DELTA

C              The function must not be evaluated too close to the endpoints.

               IF ((U - A) .LT. TOL2 .OR. (B - U) .LT. TOL2)
     >            DELTA = SIGN (TOL1, MIDDLE - X)

            ELSE

C              Proposed parabolic point failed.

               GOLDEN = .TRUE.
            END IF
         END IF

         IF (GOLDEN) THEN


C           Use a Golden-Section step.
C           --------------------------

C           Here, E is the step from X to the more distant endpoint.

            IF (X .GE. MIDDLE) THEN
               E = A - X
            ELSE
               E = B - X
            END IF
            DELTA = RATIO * E
         END IF

C        The function must not be evaluated too close to X.

         IF (ABS (DELTA) .LT. TOL1) THEN
            U = X + SIGN (TOL1, DELTA)
         ELSE
            U = X + DELTA
         END IF


         XMIN = U
         TEST = .FALSE.

C        Return to caller to evaluate the function at the chosen point.
C        --------------------------------------------------------------


      ELSE

C        Pick up where we left off, with the new function value in hand.
C        ---------------------------------------------------------------

         FU = FMIN
         NUMFUN = NUMFUN + 1
         IF (LUNOUT .GT. 0) WRITE (LUNOUT, 1020) NUMFUN, U, FU

         IF (FU .LE. FX) THEN

C           The new value is an improvement.

            IF (U .LT. X) THEN
               B = X
            ELSE
               A = X
            END IF

            V  = W
            FV = FW
            W  = X
            FW = FX
            X  = U
            FX = FU

         ELSE

C           The new value is not as good as previous best.

            IF (U .LT. X) THEN
               A = U
            ELSE
               B = U
            END IF

            IF (FU .LE. FW .OR. W .EQ. X) THEN
               V  = W
               FV = FW
               W  = U
               FW = FU
            ELSE IF (FU .LE .FV .OR. V .EQ. X .OR. V .EQ. W) THEN
               V  = U
               FV = FU
            ELSE

C              What is left over here ???  Looks like:
C
C                 FU > FX, FW, FV  &  W >< X  &  V >< X, W
C
C              This DOES happen !  Why is it not discussed ?

            END IF

         END IF

         TEST = .TRUE.          ! Go test for convergence.
         GO TO 10

      END IF

      GO TO 99


C     Terminate.
C     ----------

   20 XMIN = X
      FMIN = FX


   99 RETURN

C     Formats:

 1010 FORMAT (/, ' FMINRC:  Iteration history for calls by ', A, //
     >        ' Iter', 23X, 'X', 21X, 'F(X)')
 1020 FORMAT (1X, I3, 1P, 2E25.15)
 1030 FORMAT (/, ' FMINRC:  AX > BX is illegal: ', 1P, 2E24.16)
 1040 FORMAT (/, ' FMINRC:  Suspicious input for (max.) NUMFUN: ', I10)
 1050 FORMAT (/, ' FMINRC:  Cannot start a new minimization before the',
     >       ' previous one is done.', /, ' Offending module: ', A)

      END SUBROUTINE FMINRC
