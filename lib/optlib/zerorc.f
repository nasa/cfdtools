C+----------------------------------------------------------------------
C
      SUBROUTINE ZERORC (AX, BX, X, F, TOL, NUMFUN, CALLER, LUNOUT,
     >                   HOLD, ISTAT)
C
C     Oneliner:
C
C     ZERO of a function (1-D; no derivatives; Reverse Communication)
C     ----                                     -       -
C
C     Description:
C
C        ZERORC estimates a zero of the function F(X) in the interval
C     (AX, BX).  Function derivatives are not required.
C
C        ZERORC is an adaptation of ZEROIN (see below), intended to
C     avoid the problems of passing as arguments the names of modules
C     with assumed calling sequences (such as ZEROIN's F(X)) at the
C     expense of having to return to the calling program for every
C     function evaluation.  Such reverse communication is usually
C     preferable to the use of common blocks forced on nontrivial
C     function routines by the likes of ZEROIN.
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
C     X       R       O     X is output with the point at which the
C                           next function evaluation is required, until
C                           termination, when X is the best estimate
C                           of the location of the zero in (AX, BX).
C
C     F       R     I/O     After the initial call, F should be input
C                           with the function value corresponding to X.
C                           Upon termination, F is output with the
C                           function value corresponding to X, which
C                           should be essentially zero.
C
C     TOL     R     I       Desired length of the interval of uncertainty
C                           of the final result.  TOL >= 0.
C
C     NUMFUN  I     I/O     At the start of an iteration, NUMFUN should
C                           be input with the maximum number of function
C                           evaluations to be permitted.  On output, it
C                           contains the number of function evaluations so far.
C
C     CALLER  C*(*) I       Name of the calling program, printed if LUNOUT > 0.
C                           No more than 8 characters are significant.
C
C     LUNOUT  I     I       Logical unit number for output of the iteration
C                           history.  Use LUNOUT < 0 if you don't want the
C                           evaluations printed.  Any diagnostics go to
C                           the unit defined by ABS (LUNOUT).
C
C     HOLD (13)  R  I/O     Holds copies of the local real variables reused
C                           during the calls to ZERORC in an iteration.
C                           This permits finding a zero of a function which
C                           in turn involves finding a zero of another function.
C                           The contents should be manipulated only by ZERORC.
C
C     ISTAT   I     I/O     Status flag used as follows:
C                           ISTAT = +2 on input means this call initiates
C                                      a new iteration.
C                           ISTAT = +1 on return means the calling program
C                                      should evaluate the function at the
C                                      current value of X, and call ZERORC
C                                      again with this value in argument F
C                                      and ISTAT = +1 still.
C                           ISTAT =  0 on output means a zero has been
C                                      found to the specified tolerance.
C                           ISTAT = -1 on output means the limit on the
C                                      number of function evaluations
C                                      has been reached but the convergence
C                                      criteria have not been satisfied.
C                                      X, F are the best values found.
C                           ISTAT = -2 means AX >= BX.  This is checked only
C                                      at the start of an iteration.
C                           ISTAT = -3 means the input value of NUMFUN on
C                                      the first call is either non-positive
C                                      or absurdly large.
C                           ISTAT = -4 means the function has the same sign
C                                      at AX and BX.
C
C     Usage:
C
C        ISTAT = 2          ! Initialize the iteration
C        AX = ...
C        BX = ...
C        NUMFUN = ...       ! Max. no. of fn. evals. allowed
C                           <Define TOL and SUBNAME as well.>
C     10 CONTINUE
C
C           CALL ZERORC (AX, BX, X, F, TOL, NUMFUN, SUBNAME,
C       >                LUNOUT, HOLD, ISTAT)
C
C           IF (ISTAT .LT. -1) THEN      ! Fatal error
C
C              <Handle it>
C
C           ELSE IF (ISTAT .LT. 0) THEN  ! Iteration limit reached
C
C              <Unconverged, but X may still be usable>
C
C           ELSE IF (ISTAT .GT. 0) THEN  ! Evaluate the function
C
C              CALL fun (X, F)           ! Or whatever
C              GO TO 10
C
C           ! ELSE ISTAT = 0 - a zero has been found.
C
C           END IF
C
C 
C     Notes on the ZEROIN algorithm:
C
C       ZEROIN determines a zero X in the given interval (AX, BX) to within
C     a tolerance 4 * MACHEPS * ABS(X) + TOL, where MACHEPS is the relative
C     machine precision.
C
C       ZEROIN is discussed by Forsythe, Malcolm, and Moler in Computer
C     Methods for Mathematical Computations, Prentice-Hall, Inc. (1977).
C     It is a slightly modified translation of the ALGOL 60 procedure
C     ZERO given in Richard Brent, Algorithms for Minimization Without
C     Derivatives, Prentice-Hall, Inc. (1973).
C
C     History:
C
C     c. 1973  R. Brent            ALGOL 60 procedure ZERO.
C     c. 1977  M. Malcolm, et al.  ZEROIN version (FORTRAN 66).
C
C     04/11/92 D. Saunders         ZERORC version, to avoid the need for
C     Sterling Software/NASA Ames  common blocks for nontrivial functions.
C
C     10/29/10 D. Saunders         Matt MacLean reported that the machine
C              ERC, Inc./NASA ARC  epsilon iteration produced zero (!) with
C                                  the Sun/Oracle compiler switch for faster
C                                  but less accurate floating point arithmetic.
C                                  Therefore, use Fortran 90's intrinsic now.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   ISTAT, LUNOUT, NUMFUN
      REAL
     >   AX, BX, F, TOL, X, HOLD (13)
      CHARACTER
     >   CALLER * (*)

C     Local constants:

      REAL, PARAMETER ::
     >   HALF = 0.5E+0, ONE = 1.E+0, TWO = 2.E+0, THREE = 3.E+0,
     >   ZERO = 0.0E+0

C     Local variables:

      INTEGER
     >   LUNERR, MAXFUN
      REAL
     >   A, B, C, D, E, EPS, FA, FB, FC, TOL1, XM, P, Q, R, S
      LOGICAL
     >   BEGIN

C     System functions:

      INTRINSIC
     >   ABS, EPSILON, SIGN

      SAVE
     >   LUNERR, MAXFUN, EPS   ! Assume these are the same for any
                               ! iterations-within-iterations

C     Execution:

      IF (ISTAT .GT. 1) THEN   ! Initialize the iteration

         EPS = EPSILON (EPS)

C        Check input arguments:

         LUNERR = ABS (LUNOUT)
         IF (LUNOUT .GT. 0) WRITE (LUNOUT, 1010) CALLER

         IF (AX .GE. BX) THEN
            ISTAT = -2
            WRITE (LUNERR, 1030) AX, BX, CALLER

         ELSE IF (NUMFUN .LE. 1 .OR. NUMFUN .GT. 200) THEN
            ISTAT = -3
            WRITE (LUNERR, 1040) NUMFUN, CALLER

         ELSE
            ISTAT = 1
            MAXFUN = NUMFUN
            NUMFUN = 0
            A = AX
            B = BX
            X = A
            HOLD (1) = A
            HOLD (2) = B

         END IF

         GO TO 99  ! Return early for the first function evaluation

      ELSE IF (NUMFUN .EQ. 0) THEN

C        Second call to ZERORC.  We have F(AX) but need F(BX).

         NUMFUN = 1
         FA = F
         HOLD (3) = FA
         IF (LUNOUT .GT. 0) WRITE (LUNOUT, 1020) NUMFUN, X, FA, CALLER

         X = BX
         GO TO 99

      ELSE IF (NUMFUN .EQ. 1) THEN
 
C        We now have F at both end points.
C        These function values must have opposite signs, or at least
C        one must be nonzero.

         FB = F
!!!      HOLD (4) = FB   ! Redundant
         FA = HOLD (3)

         IF (FA .EQ. ZERO .AND. FB .EQ. ZERO) THEN
            ISTAT = -4
         ELSE IF (FA * FB .GT. ZERO) THEN
            ISTAT = -4
         END IF

         IF (ISTAT .NE. 1) THEN
            WRITE (LUNERR, 1050) FA, FB, CALLER
            GO TO 99
         END IF

      END IF

      A  = HOLD (1)
      B  = HOLD (2)
      FA = HOLD (3)

      NUMFUN = NUMFUN + 1
      BEGIN  = NUMFUN .EQ. 2

      IF (.NOT. BEGIN) THEN   ! Restore local variables from previous call
!!!      FB = HOLD (4)        ! Redundant
         FC = HOLD (5)
         C  = HOLD (6)
         D  = HOLD (7)
         E  = HOLD (8)
         XM = HOLD (9)
         P  = HOLD (10)
         Q  = HOLD (11)
         R  = HOLD (12)
         S  = HOLD (13)

C        The next test had been at the end, after FB = F(B):

         FB = F
         BEGIN = FB * (FC /ABS (FC)) .GT. ZERO
      END IF

      IF (LUNOUT .GT. 0) WRITE (LUNOUT, 1020) NUMFUN, X, FB, CALLER

      IF (BEGIN) THEN  ! Begin step

         C  = A
         FC = FA
         D  = B - A
         E  = D

      END IF

      IF (ABS (FC) .LT. ABS (FB)) THEN
         A  = B
         B  = C
         C  = A
         FA = FB
         FB = FC
         FC = FA
      END IF

C     Convergence test:

      TOL1 = TWO * EPS * ABS (B) + HALF * TOL
      XM = HALF * (C - B)
      IF (ABS (XM) .LE. TOL1) GO TO 90
      IF (FB .EQ. ZERO) GO TO 90

C     Is bisection necessary?

      IF (ABS (E) .LT. TOL1) GO TO 70
      IF (ABS (FA) .LE. ABS (FB)) GO TO 70

C     Is quadratic interpolation possible?

      IF (A .EQ. C) THEN  ! No - do linear interpolation

         S = FB / FA
         P = TWO * XM * S
         Q = ONE - S

      ELSE                ! Inverse quadratic interpolation

         Q = FA / FC
         R = FB / FC
         S = FB / FA
         P = S * (TWO * XM * Q * (Q - R) - (B - A) * (R - ONE))
         Q = (Q - ONE) * (R - ONE) * (S - ONE)

      END IF

C     Adjust signs:

      IF (P .GT. ZERO) Q = -Q
      P = ABS (P)

C     Is interpolation acceptable?

      IF ((P + P) .GE. (THREE * XM * Q - ABS (TOL1 * Q))) GO TO 70
      IF (P .GE. ABS (HALF * E * Q)) GO TO 70

      E = D
      D = P / Q
      GO TO 80

C     Bisection:

   70 D = XM
      E = D

C     Complete step:

   80 A = B
      FA = FB
      IF (ABS (D) .GT. TOL1) B = B + D
      IF (ABS (D) .LE. TOL1) B = B + SIGN (TOL1, XM)

C*****FB = F(B)

C     Return for the next function value.

      X = B
      HOLD (1) = A
      HOLD (2) = B
      HOLD (3) = FA
      HOLD (4) = FB
      HOLD (5) = FC
      HOLD (6) = C
      HOLD (7) = D
      HOLD (8) = E
      HOLD (9) = XM
      HOLD (10) = P
      HOLD (11) = Q
      HOLD (12) = R
      HOLD (13) = S

      GO TO 99


C     Done:

   90 X = B
      F = FB
      ISTAT = 0


   99 RETURN


C     Formats:

 1010 FORMAT (/, ' ZERORC:  Iteration history for calls by ', A, //
     >        ' Iter', 15X, 'X', 24X, 'F(X)')
 1020 FORMAT (1X, I3, 1P, 2E14.6, 3X, A)
 1030 FORMAT (/, ' ZERORC:  AX >= BX is illegal: ', 1P, 2E15.6,
     >        '    Caller:  ', A)
 1040 FORMAT (/, ' ZERORC:  Suspicious input for (max.) NUMFUN: ', I10,
     >        '    Caller:  ', A)
 1050 FORMAT (/, ' ZERORC:  End-point function values have same sign: ',
     >        /, 1P, 2E15.6, '   Caller: ', A)

      END SUBROUTINE ZERORC
