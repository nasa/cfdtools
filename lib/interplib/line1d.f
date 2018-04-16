C+----------------------------------------------------------------------
C
      SUBROUTINE LINE1D (NTABLE, XTABLE, YTABLE, NX, X, Y)
C
C ACRONYM: LINEar interpolation/extrapolation (1D)
C
C PURPOSE: For each X(I),  LINE1D gives a linearly-interpolated Y(I)
C          from the given 1D table, whose abscissas must be monotone
C          increasing or monotone decreasing.
C
C          This version permits extrapolation.
C
C          LINE1D is the subroutine form of FUNCTION TABLE1, for the
C          case when many linear interpolations in a given table are
C          required  (implying many external references to  TABLE1).
C
C METHOD:  Nothing more than linear interpolation  (2-point formula)
C          is attempted,  partly because multi-point formulas demand
C          special treatment at the ends of the table, and partly to
C          permit "safe" extrapolation, which is often due to round-
C          off or is otherwise intended. Furthermore, piecewise line
C          segments (as the table may well represent) MUST be treat-
C          ed linearly.
C
C ARGUMENTS:
C    ARG    DIM   TYPE I/O/S DESCRIPTION
C  NTABLE    -      I    I   Length of the table; NTABLE > 1 assumed
C  XTABLE  NTABLE   R    I   Abscissas in table, assumed monotonic -
C                            increasing or decreasing
C  YTABLE  NTABLE   R    I   Ordinates in table
C    NX      -      I    I   Number of interpolations required
C    X       NX     R    I   Abscissas at which interpolation is reqd.
C    Y       NX     R    O   Interpolated/extrapolated values
C
C EXTERNAL REFERENCES:
C    NAME     DESCRIPTION
C    INTERVAL Interpolation search for interval containing a point.
C
C ERROR HANDLING:
C     None - checking for monotonicity in XTABLE(*) would compromise
C     efficiency too much.  FORTRAN handles NTABLE <=0.
C
C ENVIRONMENT:  VAX/VMS FORTRAN 77
C
C DEVELOPMENT HISTORY:
C     DATE  INITIALS  DESCRIPTION
C   10/21/83   DAS    Initial design and code (sequential search; no
C                     extrapolation permitted) as FUNCTION TABLE1.
C   06/06/86   DAS    Extrapolation now permitted.  Adapted so-called
C                     "interpolation search" from CSEVAL (by Robert
C                     Kennelly) to handle monotonic decreasing as well
C                     as increasing case.
C   08/08/86   DAS    TRIAL = MIN (NTABLE-2,...), not MIN (NTABLE-1,...)
C   02/19/87   DAS    LINE1D adapted from TABLE1.
C   10/20/87   RAK    Abstract the interval search (with revisions) as
C                     a separate module.
C
C AUTHOR: David Saunders, Sterling Software, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   NTABLE, NX
      REAL
     >   X (NX), XTABLE (NTABLE), Y (NX), YTABLE (NTABLE)

C     Local variables:

      INTEGER
     >   I, LEFT
      REAL
     >   PM1

C     Execution:

      PM1  = SIGN (1.E+0, XTABLE (2) - XTABLE (1))
      LEFT = 1

      DO 10, I = 1, NX

C        Search for the best "left-hand" endpoint to use for interpolation.

         CALL INTERVAL (NTABLE, XTABLE, X (I), PM1, LEFT)

C        Linear interpolation or extrapolation.  Order does not matter.

         Y (I) = YTABLE (LEFT) + (YTABLE (LEFT+1) - YTABLE (LEFT)) *
     >      ((X (I) - XTABLE (LEFT)) /
     >      (XTABLE (LEFT+1) - XTABLE (LEFT)))
   10 CONTINUE

      RETURN
      END
