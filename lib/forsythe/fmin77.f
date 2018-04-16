C+----------------------------------------------------------------------
C
      SUBROUTINE FMIN77 (LUNOUT, AX, BX, TOL, OPTX, OPTF, FUN)
C
C        FMIN77 estimates a point OPTX where a function of 1 variable,
C     supplied as FUN, attains a minimum on the interval (AX, BX).  The
C     corresponding function value is returned as OPTF.  This is a
C     restructured, subroutine version of Brent's FMIN function routine,
C     referenced below.
C
C     Input..
C
C     LUNOUT Logical unit number for output of iteration history.
C            Use LUNOUT < 0 if you don't want the evaluations printed.
C     AX     Left endpoint of initial interval.
C     BX     Right endpoint of initial interval.
C     TOL    Desired length of the interval of uncertainty of the final
C            result (.GE. 0.0E+0).
C     FUN    FUNCTION subprogram which evaluates FUN (X) for any X
C            in the interval (AX, BX).
C
C     Output..
C
C     OPTX   Approximate point where FUN attains a minimum.
C     OPTF   Corresponding function value, FUN (OPTX).
C
C        The method used is a combination of Golden Section search and
C     successive parabolic interpolation.  Convergence is never much slower
C     than that for a Fibonacci search.  If FUN has a continuous second
C     derivative which is positive at the minimum (which is not at AX or
C     BX), then convergence is superlinear, and usually of the order of
C     about 1.324.
C
C        The function FUN is never evaluated at two points closer together
C     than RTEPS * ABS (FMIN) + (TOL/3), where RTEPS is approximately the
C     square root of the relative machine precision.  If FUN is a unimodal
C     function and the computed values of FUN are always unimodal when
C     separated by at least RTEPS * ABS (X) + (TOL/3), then FMIN
C     approximates the abcissa of the global minimum of FUN on the interval
C     (AX, BX) with an error less than 3 * RTEPS * ABS (FMIN) + TOL.  If
C     FUN is not unimodal, then FMIN may approximate a local, but perhaps
C     non-global, minimum to the same accuracy.
C
C        This is a modified version of the ALGOL 60 procedure LOCALMIN given
C     in Richard Brent, Algorithms for Minimization without Derivatives,
C     Prentice-Hall, Inc. (1973).  The FORTRAN 77 translation started from
C     the FORTRAN 66 version, FMIN, in Forsythe, Malcolm, and Moler but
C     also made use of the original ALGOL.
C
C     Notes:
C
C        (1)  IMPLICIT NONE is non-standard.
C
C        (2)  RTEPS is machine- and precision-dependent.  This version is
C             set up for a VAX using SINGLE precision.
C
CRAK  May want to provide for max. number of cycles, etc.
CRAK  What about handling termination where the test is, rather than
CRAK  at the end ?
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT NONE

C     RTEPS is approximately the square root of the relative machine
C     precision, and RATIO is (3 - SQRT (5)) / 2, the squared inverse
C     of the Golden Ratio.

      REAL
     >   ZERO, ONE, TWO, THREE, HALF, THIRD, RTEPS, RATIO
      PARAMETER
     >   (ZERO = 0.0E+0, ONE = 1.0E+0, TWO = 2.0E+0, THREE = 3.0E+0,
     >   HALF = ONE / TWO, THIRD = ONE / THREE,
     >   RTEPS = 2.44E-4, RATIO = .381966011250105E+0)

      LOGICAL
     >   GOLDEN
      INTEGER
     >   LUNOUT, NUMFUN
      REAL
     >   AX, BX, FUN, TOL, A, B, DELTA, E, MIDDLE, OPTF, OPTX,
     >   P, Q, R, TOL1, TOL2, U, V, W, X, FU, FV, FW, FX

      INTRINSIC
     >   ABS, SIGN
      EXTERNAL
     >   FUN


C     Initialization.
C     ---------------

      IF (LUNOUT .GT. 0) WRITE (LUNOUT, 1010)
CRAK
CRAK  How about checking inputs, e.g. length of initial interval, TOL, etc.
CRAK  Must we have B > A ?
CRAK

C     Notation:  At the start of a cycle, we have
C
C        A, B   Left and right endpoints of current interval of uncertainty
C        X      Point at which FUN is smallest, or the most recent point
C               if there was a tie
C        U      Last point at which FUN was evaluated
C        W      Second best point thus far
C        V      Third best point thus far
C        DELTA  Proposed offset from latest point (U = X + DELTA)
C        E      A measure of the previous step size (precisely old DELTA when
C               a parabolic step was taken)


      A  = AX
      B  = BX
      E  = ZERO

      X  = A + RATIO * (B - A)
      W  = X
      V  = X

      FX = FUN (X)
      NUMFUN = 1
      FW = FX
      FV = FX
      IF (LUNOUT .GT. 0) WRITE (LUNOUT, 1020) NUMFUN, X, FX

C*****write (LUNOUT, 3456)
C3456 format (1X, t10, 'a', t20, 'b', t25, 'delta', t40, 'e', t50,
C****>   'x', t62, 'w', t73, 'v', t84, 'fx', t95, 'fw', t106, 'fv'//)


   10 CONTINUE


C*****   write (LUNOUT, 2345) numfun, a, b, delta, e, x, w, v, fx,fw,fv
C2345    format (1x, i2, 10(2x, e9.3))


         MIDDLE = HALF * (A + B)
CRAK
CRAK     Note Brent doesn't use THIRD below (added in FM&M) ?
CRAK
         TOL1 = RTEPS * ABS (X) + THIRD * TOL
         TOL2 = TWO * TOL1


C        Check stopping criterion.
C        -------------------------

C        The distance from X to the midpoint of the interval of uncertainty
C        plus half the width of the interval must be less than twice the
C        tolerance.

CRAK
CRAK     Equivalent to:  X-A .LE. TOL2  when X > MIDDLE   ...and...
CRAK                     B-X .LE. TOL2  when X < MIDDLE
CRAK     This isn't the same as the error estimate in the header. (It's
CRAK     actually a little tighter.)
CRAK

         IF (ABS (X - MIDDLE) .LE. (TOL2 - HALF * (B - A))) GO TO 20

C        Compute a trial step, by parabolic interpolation if possible.

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

C              FUN  must not be evaluated too close to the endpoints.

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

C        FUN must not be evaluated too close to X.

         IF (ABS (DELTA) .LT. TOL1) THEN
            U = X + SIGN (TOL1, DELTA)
         ELSE
            U = X + DELTA
         END IF


C        Evaluate the function at the chosen point.
C        ------------------------------------------

         FU = FUN (U)
         NUMFUN = NUMFUN + 1
         IF (LUNOUT .GT. 0) WRITE (LUNOUT, 1020) NUMFUN, U, FU

C*****   if (golden) then
C*****      write (LUNOUT, *) 'GOLDEN step, NUMFUN = ', NUMFUN
C*****   else
C*****      write (LUNOUT, *) 'PARABOLIC step, NUMFUN = ', NUMFUN
C*****   end if

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

C*****         WRITE (LUNOUT, *) 'CATCH-ALL:  U, FU = ', U, FU
            END IF

         END IF
         GO TO 10


C     Terminate.
C     ----------

   20 CONTINUE
      OPTX = X
      OPTF = FX

      RETURN

C     FORMAT Statements.
C     ------------------

 1010 FORMAT ('0Iteration history from FMIN77:'//
     >        ' Iter', 15X, 'X', 24X, 'F(X)')
 1020 FORMAT (1X, I3, 1P, 2E14.6)
      END
