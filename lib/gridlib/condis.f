C+----------------------------------------------------------------------
C
      SUBROUTINE CONDIS (N, M, X1, XM, XN, X)
C
C  PURPOSE:
C     CONDIS (CONstrained DIStribution) generates an exponential-type
C     distribution of points X (1:N) in the range [X1, XN] such that
C     interior point X (M) equals XM exactly, for specified M and XM.
C     The spacing X (I+1) - X (I) may steadily increase or steadily
C     decrease with I - see the method description below.
C
C  METHOD:
C     The mesh points are defined by the formula (for I = 1:N),
C
C     X (I) - X1     EXP (ALPHA * T (I)) - 1.                 I - 1
C     ----------  =  ------------------------  where  T (I) = ----- .
C     XN  -  X1         EXP (ALPHA) - 1.                      N - 1
C
C     To constrain X (M), the I = M equation may be solved for ALPHA.
C     A standard Newton iteration is used, but there are three cases
C     to deal with for this nonlinear equation, according to how the
C     left-hand side above compares with T (M).  Without loss of
C     generality, consider the [0, 1] interval for [X1, XN]:
C
C     (1) If XM < T (M) then ALPHA > 0; use ALPHA = 1.0 as starting guess.
C     (2) If XM > T (M) then ALPHA < 0; use ALPHA = -1.0  "    "    ".
C     (3) If XM = T (M) then the distribution must be uniform.  (ALPHA=0
C         is a singular case, and must be avoided - here we use the same
C         tolerance that is used for the convergence test.)
C
C  ARGUMENTS:
C     ARG  TYPE  I/O/S  DIM   DESCRIPTION
C     N      I     I     -    Required number of points. N > 2.
C     M      I     I     -    Interior point to be constrained. 1 < M < N.
C     X1     R     I     -    Defines X (1).
C     XM     R     I     -     "   "  X (M);  X1 < XM < XN.
C     XN     R     I     -     "   "  X (N).
C     X      R     O     N    Desired distribution.
C
C  ERROR HANDLING:
C    Failure to converge is almost certain to mean bad input data (though
C    ALPHA = 1.0 or -1.0 for the starting guesses have not been PROVED to
C    guarantee convergence).  Actually, |ALPHA| < 10. for reasonable dis-
C    tributions, and testing on cases where ALPHA ~ +/-10 suggests that
C    this starting guess is adequate.
C
C    To avoid logical unit number problems, just STOP 'CONDIS: .....' if
C    the iteration count limit (20) is reached.
C
C  ENVIRONMENT:  VAX/VMS FORTRAN 77
C
C  DEVELOPMENT HISTORY:
C    DATE    DESCRIPTION
C  05/13/86  Jeff Cordova (Sterling Software):
C            Original implementation as CLUSTER:  ALPHA > 0 case only.
C  10/09/87  David Saunders (Sterling):
C            Refined for inclusion in library as CONDIS:  bunching at
C            either end permitted; iteration safeguarded.
C  10/28/87  Phillip Snyder (Sterling):
C            Testing suggests iteration will converge for reasonable
C            distributions (|ALPHA| < 10 or so);  convergence test needed
C            |ALPHA| included; |ALPHA| also constrained with TOOBIG ~ 80. to
C            prevent overflow for extreme cases (which will diverge anyway).
C
C  CONCERNS:
C     (1) Choice of subscript M prevents direct control (via ALPHA) of
C         the amount of stretching.  Reasonable choices are expected.
C     (2) Should ALPHA be input as a potentially better starting guess?
C         (Fairly safe to say convergence won't be a problem normally.)
C     (3) Should ALPHA be returned as a measure of the stretching?
C         (Ratio of first and last interval sizes serves similar purpose.)
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   M, N
      REAL
     >   X (N), X1, XM, XN

C     Local constants:

      INTEGER
     >   MXITER
      REAL
     >   ONE, TOL, TOOBIG
      PARAMETER
     >  (MXITER = 20,
     >   ONE = 1.E+0,
     >   TOL = 1.E-6,
     >   TOOBIG = 80.E+0)   ! EXP (87...) ~ 1.E+37 ~ VAX limit

C     Local variables:

      INTEGER
     >   I, ITER
      REAL
     >   A, B, C, D, ALPHA, DALPHA, DX, RHS, RNM1, TERM, TM
      LOGICAL
     >   BUNCHX1

C     System functions:

      INTRINSIC
     >   ABS, EXP, MAX, MIN

C     Execution:

      RHS  = (XM - X1) / (XN - X1)
      RNM1 = ONE / (N - 1)
      TM   = (M - 1) * RNM1
      ITER = 0

      IF (ABS (TM - RHS) .LT. TOL) THEN   ! Singular case - set uniform X (*):

         DX = (XN - X1) * RNM1
         DO 20, I = 2, N - 1
            X (I) = X1 + (I - 1) * DX
   20    CONTINUE

      ELSE                                ! Need to solve nonlinear equation
                                          ! for ALPHA.  Set up starting guess:
         BUNCHX1 = TM .GE. RHS

         IF (BUNCHX1) THEN                ! Bunching towards X1:
            ALPHA = ONE
         ELSE                             ! Bunching towards XN:
            ALPHA = -ONE
         END IF

         ITER = 0
   30    CONTINUE
            A = EXP (ALPHA * TM)
            B = EXP (ALPHA)
            C = A - ONE
            D = B - ONE                   ! DALPHA = F(X) / F'(X) for F(X) = 0:


            DALPHA = ((C / D) - RHS) / (TM * (A / D) - B * C / (D * D))

            IF (BUNCHX1) THEN             ! Guard against zero from above:
               ALPHA = MIN (MAX (TOL, ALPHA - DALPHA),
     >                      TOOBIG)       ! ...and against overflow.
            ELSE                          ! Guard against zero from below:
               ALPHA = MIN (-TOL, ALPHA - DALPHA)
            END IF

            ITER = ITER + 1
            IF (ABS (DALPHA) .GT. TOL * ABS (ALPHA) .AND.
     >         ITER .LT. MXITER)
     >   GO TO 30

         IF (ITER .EQ. MXITER) THEN
            STOP 'CONDIS: iteration appears to diverge ... quitting.'
         END IF

C        Generate the distribution using the iteratively computed ALPHA:

         TERM = (XN - X1) / (EXP (ALPHA) - ONE)
         DO 40, I = 2, N - 1
            X (I) = X1 + (EXP ((I - 1) * (RNM1 * ALPHA)) - ONE) * TERM
   40    CONTINUE

      END IF

C     Make sure the specified points are exactly right:

      X (1) = X1
      X (M) = XM
      X (N) = XN

      RETURN
      END        ! End of CONDIS
