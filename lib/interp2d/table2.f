C+----------------------------------------------------------------------
C
      FUNCTION TABLE2 (NX, NY, IDIMX, XTBL, YTBL, FTBL, X, Y, IER)
C
C ACRONYM: TABLE look-up, 2 dimensions.
C          -----          -
C PURPOSE: TABLE2  returns a linearly-interpolated value from a given
C          2-D table for the given  (X,Y)  pair.  The table should be
C          rectangular (NX entries in each row and NY entries in each
C          column).  The coordinates need not be uniformly spaced but
C          they must be monotonic increasing or monotonic decreasing.
C
C          This version permits linear extrapolation.
C
C METHOD:  The efficient search used by TABLE1 is applied for each of
C          the two coordinates.  (Simple reuse of TABLE1 is less than
C          optimal, and is no longer necessary now that INTERVAL is a
C          distinct module.)
C
C ARGUMENTS:
C    ARG    DIM   TYPE I/O/S DESCRIPTION
C    NX      -      I    I   Number of columns in table (NX>=2)
C    NY      -      I    I   Number of rows    in table (NY>=2)
C  IDIMX     -      I    I   Declared first dimension of FTBL(*,*)
C   XTBL    NX      R    I   Table coordinates in the 1st dimension
C   YTBL    NY      R    I   Table coordinates in the 2nd dimension
C   FTBL IDIMX,NY   R    I   Table function values
C   X,Y      -      R    I   Point at which interpolation is reqd.
C   IER      -      I    O   Error return code - see ERROR HANDLING
C  TABLE2    -      R    O   FUNCTION  value is returned as that from
C                            linear interpolation or extrapolation.
C                            TABLE2=0. if an error is detected (rather
C                            than being left undefined).
C                     
C EXTERNAL REFERENCES:
C  INTERVAL   1-D search utility
C
C ERROR HANDLING:
C          IER = 0 if no error was detected.
C              = 1 if NX < 2;
C              = 2 if NY < 2.
C          There is no check for monotonicity in X(*) and Y(*).
C
C ENVIRONMENT:  VAX/VMS FORTRAN 77
C
C DEVELOPMENT HISTORY:
C     DATE  INITIALS  DESCRIPTION
C   10/02/85   DAS    Initial design and code, adapted from TABLE1.
C   06/06/86   DAS    Revised to suit revised TABLE1 (permitting
C                     extrapolation, with revised error handling).
C   10/28/87   DAS    Advent of INTERVAL as a distinct search module meant
C                     reuse of TABLE1 was unnecessarily indirect here.
C
C AUTHOR: David Saunders, Sterling Software, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   IDIMX, NX, NY, IER
      REAL
     >   TABLE2, X, Y, XTBL (NX), YTBL (NY), FTBL (IDIMX, NY)

C     Local constants:

      REAL
     >   ONE
      PARAMETER
     >  (ONE = 1.E+0)

C     Local variables:

      INTEGER
     >   I, J
      REAL
     >   P, Q

C     Procedures:

      EXTERNAL
     >   INTERVAL

C     Execution:

      TABLE2 = 0.E+0
      IER = 0

      IF (NX .LT. 2) THEN
         IER = 1
         GO TO 99
      END IF

      IF (NY .LT. 2) THEN
         IER = 2
         GO TO 99
      END IF

      I = 1
      CALL INTERVAL (NX, XTBL, X, SIGN (ONE, XTBL (2) - XTBL (1)), I)

      P = (X - XTBL (I)) / (XTBL (I+1) - XTBL (I))

      J = 1
      CALL INTERVAL (NY, YTBL, Y, SIGN (ONE, YTBL (2) - YTBL (1)), J)

      Q = (Y - YTBL (J)) / (YTBL (J+1) - YTBL (J))

C     Now I and J are subscripts to the "left" of X and Y, or at X or Y
C     exactly, unless we are extrapolating, in which case the formula
C     still applies using points 1 & 2 or N-1 & N in either direction:

      TABLE2 = (ONE - P) * (ONE - Q) * FTBL (I, J)   +
     >                 P * (ONE - Q) * FTBL (I+1, J) +
     >         (ONE - P) * Q         * FTBL (I, J+1) +
     >                 P * Q         * FTBL (I+1, J+1)
      
   99 RETURN
      END
