C+----------------------------------------------------------------------
C
      SUBROUTINE EXPDIS4 (MODE, XA, XB, DX, N, X, BETA, LUNOUT)
C
C  PURPOSE:
C
C        EXPDIS4 is the all-REAL version of EXPDIS2, for 64-bit systems.
C
C        EXPDIS2 (EXPonential DIStribution, 2nd form) generates N values
C     X(I) in the range [XA, XB] such that they form an exponential-type
C     distribution with a smallest interval equal to the given DX. Three
C     cases are provided:
C
C        MODE       Bunching towards ...         DX = ...
C        ------------------------------------------------------------
C          1               XB                  X(N) - X(N-1)
C          2            XA and XB       X(2) - X(1) and X(N) - X(N-1)
C          3               XA                  X(2) - X(1)
C
C     The two-sided bunching case is symmetric.  No attempt is made to
C     achieve a given largest interval, although the analysis would be
C     similar.        -------
C
C        See also the original module EXPDIS, which uses a given BETA to
C     produce a distribution in normalized form with no direct control
C     over the smallest interval.
C
C        For the UNsymmetric case of specified first and last intervals,
C     see HTDIS2, which gives the same result as EXPDIS2 if the intervals
C     are equal (although the implementations are quite different).
C
C  METHOD:
C
C        EXPDIS2 applies the analysis of EXPDIS to solve the "inverse"
C     problem:  "What BETA gives a specified DX?" as opposed to "What DX
C     does a given BETA produce?" by using a zero-finder to solve the
C     relevant nonlinear equation in BETA.
C
C  ARGUMENTS:
C     ARG    TYPE  I/O/S  DIM   DESCRIPTION
C     MODE     I     I     -    MODE=1 gives bunching towards XB.
C                               MODE=2 gives bunching at XA and XB.
C                               MODE=3 gives bunching towards XA.
C     XA       R     I     -    Left- and right-hand ends of the
C     XB       R     I     -    total interval to be distributed.
C                               Note that use of the interval [0., 1.]
C                               can give a RELATIVE distribution that
C                               may be readily applied to various other
C                               intervals.
C     DX       R     I     -    Desired length of the smallest interval
C                               in the distribution (see PURPOSE above).
C     N        I     I     -    Required number of points. N > 1.
C     X        R     O     N    Values varying nonuniformly from XA to
C                               XB, inclusive.
C     BETA     R     O     -    Value of the distribution parameter
C                               estimated here to give the desired DX.
C                               BETA is typically between 1.0001 and 2.
C                               BETA may be of no interest to the calling
C                               program but is returned just in case.
C     LUNOUT   I     I     -    Logical unit |LUNOUT| receives error messages;
C                               LUNOUT < 0 suppresses iteration history.
C
C  ERROR HANDLING:
C
C        Sufficient conditions for the existence of a solution have not
C     been determined.  Bad combinations of XA, XB, N, and DX are certainly
C     possible.  For example, if DX = (XB - XA) / (N-1) then the distribution
C     would be uniform, but this limiting case (infinite BETA) is not handled.
C     In fact, the chosen search interval for BETA is [1.00000001, 10.], which
C     is assumed to cover all likely applications.  If the zero-finder indicates
C     failure to find a solution, the input arguments are written to |LUNOUT|
C     and execution proceeds with BETA = 1.00000001D+0.
C
C  ENVIRONMENT:  FORTRAN 77, 64-bit arithmetic, with
C                IMPLICIT NONE, 8-character names, and ! comments
C
C  EXTERNAL REFERENCES:
C
C     ZERORC   General-purpose reverse-communication 1-D zero-finder,
C              REAL version.
C
C  HISTORY:
C
C  08/21/85    DAS    Original implementation of EXPDIS (after a formula
C                     found in the SCRAM2D flow solver of A.J.Kumar).
C  02/25/88    DAS    EXPDIS2 derived from the analysis of EXPDIS.
C  08/09/89    DAS    Calculation of BETA is in DOUBLE PRECISION now to
C                     to help requirement for initial dX of 1.E-7 on [0,1].
C                     If ZEROIN fails, keep going with very small BETA.
C  05/04/91    DAS    Description updated slightly in light of HTDIS2, etc.
C  12/01/94    DAS    Substituted ZERORC for ZEROIN to avoid COMMON and to
C                     avoid having both forms in one application.
C   "   "      DAS    EXPDIS4 adapted from EXPDIS2 for 64-bit systems.
C
C  AUTHOR: David Saunders, Sterling Software/NASA Ames, Mt. View, California
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   LUNOUT, MODE, N
      REAL
     >   BETA, DX, X (N), XA, XB

C     Local constants:

      INTEGER
     >   MAXFUN
      REAL
     >   ONE, TWO, THREE, BETAMAX, BETAMIN, TOL
      CHARACTER
     >   SUBNAME * 7

      PARAMETER
     >  (BETAMAX = 10.E+0,       ! Search interval for the value of BETA
     >   BETAMIN = 1.00000001E+0,! which produces the desired DX
     >   MAXFUN  = 30,           ! Max. # zero-finder function evaluations
     >   ONE     = 1.E+0,
     >   SUBNAME = 'EXPDIS4',
     >   THREE   = 3.E+0,
     >   TOL     = 0.E+0,        ! Ask ZERORC for full precision
     >   TWO     = 2.E+0)

C     Local variables:

      INTEGER
     >   I, ISTAT, LUNERR, NUMFUN
      REAL
     >   B, DBLEBETA, DELTA, FUN, EXPN, FI, GAMMA, HOLD (13),
     >   POWER, RANGE, RNM1
      LOGICAL
     >   ONESIDED

C     Procedures:

      EXTERNAL
     >   ZERORC

C     Execution:

C     Original set-up for DOUBLE PRECISION internal arithmetic is now REAL:

      DELTA = DX
      EXPN  = REAL (N)
      RANGE = XB - XA
      ONESIDED = MODE .NE. 2

      IF (ONESIDED) THEN
         POWER = (EXPN - TWO) / (EXPN - ONE)
      ELSE
         POWER = (EXPN - THREE) / (EXPN - ONE)
         RANGE = RANGE / TWO
      END IF

C     Initialize the reverse-communication zero-finding loop:

      ISTAT  = 2
      NUMFUN = MAXFUN

   10 CONTINUE

         CALL ZERORC (BETAMIN, BETAMAX, DBLEBETA, FUN, TOL, NUMFUN,
     >                SUBNAME, LUNOUT, HOLD, ISTAT)

         IF (ISTAT .LT. 0) THEN                     ! MAXFUN reached or ...
            IF (ISTAT .LT. -1) DBLEBETA = BETAMIN   ! ... unrecoverable

            LUNERR = ABS (LUNOUT)
            WRITE (LUNERR, 1000) ISTAT, XA, XB, DX, N, MODE, DBLEBETA

         ELSE IF (ISTAT .GT. 0) THEN                ! Evaluate the function

            GAMMA = ((DBLEBETA + ONE) / (DBLEBETA - ONE)) ** POWER
            FUN = (ONE - DBLEBETA * (GAMMA - ONE) / (GAMMA + ONE)) *
     >            RANGE - DELTA
            GO TO 10

C*****   ELSE ISTAT .EQ. 0          ! A zero has been found for Beta
         END IF


      BETA = DBLEBETA

C     Now we have the EXPDIS situation, except it would be indirect to
C     generate a relative distribution on [0, 1] then go to [XA, XB] -
C     better to do it directly.

C     Note that for F(I) = A ** P(I), where A is independent of I,
C     it is cheaper to use the form  F(I) = EXP (LOG(A) * P(I)):

      B = LOG ((DBLEBETA + ONE) / (DBLEBETA - ONE))
      RNM1 = B / REAL (N - 1)

      IF (MODE .EQ. 2) THEN

C        Symmetric bunching towards both XA and XB:

         DO I = 1, (N + 1) / 2
            FI = EXP (RNM1 * REAL (2*I-N-1))
            X (I) = RANGE * (ONE + DBLEBETA * (FI - ONE) / (FI + ONE))
     >              + XA
            X (N + 1 - I) = (XA + XB) - X (I)      ! XB - (X (I) - XA)
         END DO

      ELSE IF (MODE .EQ. 1) THEN

C        Bunch towards XB:

         DO I = 1, N
            FI = EXP (RNM1 * REAL (I - 1))
            X (I) = (RANGE * DBLEBETA) * (FI - ONE) / (FI + ONE) + XA
         END DO

      ELSE ! (MODE = 3)

C        Bunch towards XA:

         DO I = 1, N
            FI = EXP (RNM1 *  REAL (N - I))
            X (I) = RANGE * (ONE - DBLEBETA * (FI - ONE) / (FI + ONE))
     >              + XA
         END DO

      END IF

      RETURN

C     Formats:

 1000 FORMAT (///, ' EXPDIS4:  Bad return from ZERORC. ISTAT = ', I2,
     >        /, ' XA, XB, DX =', 1P, 3E15.6,
     >        /, ' N = ', I4, '   MODE = ', I1,
     >        /, ' Proceeding with BETA =', E15.6)
      END
