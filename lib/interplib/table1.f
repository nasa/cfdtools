C+----------------------------------------------------------------------
C
      FUNCTION TABLE1 (NTABLE, XTABLE, YTABLE, INDEX, X, IER)
C
C ACRONYM: TABLE look-up, 1 dimension (linear inter-/extra-polation)
C          -----          -
C PURPOSE: TABLE1 returns a linearly-interpolated value from a given
C          1-D table for the given abscissa. It also returns the in-
C          dex of the table abscissa nearest to the given  abscissa,
C          on the low side, as this can speed subsequent look-ups in
C          some contexts.    The table's abscissas are assumed to be
C          monotonic increasing or monotonic decreasing.
C
C          This version permits extrapolation.
C
C          NOTE:  The subroutine form of TABLE1, called LINE1D,  may
C          be more appropriate if an array of target Xs is involved.
C
C METHOD:  Nothing more than linear interpolation  (2-point formula)
C          is attempted, because multi-point formulas demand special
C          treatment at the ends of the table, and because it is de-
C          sired to permit "safe" extrapolation.  (This is often due
C          to round-off or is otherwise intended, and is unlikely to
C          cause grief.)
C
C ARGUMENTS:
C    ARG    DIM   TYPE I/O/S DESCRIPTION
C  NTABLE    -      I    I   Length of the table; NTABLE > 1 assumed
C  XTABLE  NTABLE   R    I   Abscissas in table, assumed monotonic -
C                            increasing or decreasing
C  YTABLE  NTABLE   R    I   Ordinates in table
C  INDEX     -      I   I/O  Input as an estimate of the interval in
C                            the table which contains the target X -
C                            normally INDEX=1 on the first call in a
C                            sequence, and left as set by TABLE1 for
C                            subsequent calls in the sequence;
C                            output points to the element before X
C                            (or exactly at X) in the table;
C                            if extrapolation occurred to the LEFT
C                            of the table, output INDEX=1;
C                            if extrapolation occurred to the RIGHT,
C                            output INDEX=NTABLE-1 (not NTABLE),
C                            meaning the last two table values were
C                            used to do the extrapolation.
C    X       -      R    I   Abscissa at which interpolation is reqd.
C   IER      -      I    O   Error return code - see ERROR HANDLING.
C                            (IER=0 means no error detected.)
C  TABLE1    -      R    O   Interpolated or extrapolated value.
C                            TABLE1=0.0 if an error is detected (rather
C                            than being left undefined).
C EXTERNAL REFERENCES:
C    NAME     DESCRIPTION
C    INTERVAL Interpolation search for interval containing a point.
C
C ERROR HANDLING:
C          Once extrapolation was permitted, it was tempting to remove
C          the error return argument altogether.  However, the follow-
C          ing input checks may save programmers some grief:
C             IER = 0 if no error was detected;
C                 = 1 if NTABLE < 2;
C                 = 2 if INDEX < 1;
C                 = 3 if INDEX > NTABLE-1.
C          There is no check for monotonicity in XTABLE(*).
C
C ENVIRONMENT:  VAX/VMS FORTRAN 77
C
C DEVELOPMENT HISTORY:
C     DATE  INITIALS  DESCRIPTION
C   10/21/83   DAS    Initial design and code (sequential search; no
C                     extrapolation permitted).
C   06/06/86   DAS    Extrapolation now permitted.  Adapted so-called
C                     "interpolation search" from CSEVAL (by Robert
C                     Kennelly) to handle monotonic decreasing as well
C                     as increasing case.
C   08/08/86   DAS    TRIAL = MIN (NTABLE-2,...), not MIN (NTABLE-1,...)
C   10/20/87   RAK    Abstract the interval search (with revisions) as
C                     a separate module.
C
C AUTHOR: David Saunders, Sterling Software, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   INDEX, NTABLE, IER
      REAL
     >   TABLE1, X, XTABLE (NTABLE), YTABLE (NTABLE)

C     Local variables:

      INTEGER
     >   LEFT
      REAL
     >   PM1

C     Execution:

      IER = 0
      TABLE1 = 0.E+0
      LEFT = INDEX

      IF (NTABLE .LT. 2) THEN
         IER = 1
         GO TO 99
      END IF

      IF (LEFT .LT. 1) THEN
         IER = 2
         GO TO 99
      END IF

      IF (LEFT .GT. NTABLE - 1) THEN
         IER = 3
         GO TO 99
      END IF

      PM1 = SIGN (1.E+0, XTABLE (LEFT+1) - XTABLE (LEFT))

C     Search for the best "left-hand" endpoint to use for interpolation.

      CALL INTERVAL (NTABLE, XTABLE, X, PM1, LEFT)

C     Linear interpolation or extrapolation.  Order does not matter.

      TABLE1 = YTABLE (LEFT)  + (YTABLE (LEFT+1) - YTABLE (LEFT)) *
     >   ((X - XTABLE (LEFT)) / (XTABLE (LEFT+1) - XTABLE (LEFT)))
      INDEX = LEFT

   99 CONTINUE
      RETURN
      END
