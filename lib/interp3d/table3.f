C+----------------------------------------------------------------------
C
      FUNCTION TABLE3 (NX, NY, NZ, IDIMX, IDIMY, XTBL, YTBL, ZTBL,
     >                 FTBL, X, Y, Z, IER)
C
C ACRONYM: TABLE look-up, 3 dimensions.
C          -----          -
C PURPOSE: TABLE3  returns a linearly-interpolated value from a given
C          3-D table for the given  (X,Y,Z)  combination.  The coordinates
C          need not be uniformly spaced but they must be monotonic
C          increasing or monotonic decreasing.
C
C          This version permits linear extrapolation.
C
C METHOD:  A linear combination of two lower order interpolations is
C          used here.  One lower order interpolation is done in the
C          (X,Y,k) plane and the other in the (X,Y,k+1) plane.  A
C          linear combination of these two values is then made.
C
C ARGUMENTS:
C    ARG    DIM      TYPE I/O/S DESCRIPTION
C    NX      -         I    I   Number of columns in table (NX>=2)
C    NY      -         I    I   Number of rows    in table (NY>=2)
C    NZ      -         I    I   Number of layers  in table (NZ>=2)
C  IDIMX     -         I    I   Declared first dimension of FTBL(*,*,*)
C  IDIMY     -         I    I   Declared second dimension of FTBL(*,*,*)
C  XTBL     NX         R    I   Table coordinates in the 1st dimension
C  YTBL     NY         R    I   Table coordinates in the 2nd dimension
C  ZTBL     NZ         R    I   Table coordinates in the 3rd dimension
C  FTBL IDIMX,IDIMY,NZ R    I   Table function values
C  X,Y,Z     -         R    I   Point at which interpolation is reqd.
C  IER       -         I    O   Error return code - see ERROR HANDLING
C  TABLE3    -         R    O   FUNCTION  value is returned as that from
C                               linear interpolation or extrapolation.
C                               TABLE3=0. if an error is detected (rather
C                               than being left undefined).
C                     
C EXTERNAL REFERENCES:
C  INTERVAL   1-D search utility
C
C ERROR HANDLING:
C          IER = 0 if no error was detected.
C              = 1 if NX < 2;
C              = 2 if NY < 2;
C              = 3 if NZ < 2.
C          There is no check for monotonicity in X(*), Y(*), and Z(*).
C
C ENVIRONMENT:  VAX/VMS FORTRAN 77
C
C DEVELOPMENT HISTORY:
C     DATE  INITIALS  DESCRIPTION
C   05/02/88   PTS    Initial design and code, adapted from TABLE2.
C
C AUTHOR: Phillip T. Snyder, Sterling Software, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   IDIMX, IDIMY, NX, NY, NZ, IER
      REAL
     >   TABLE3, X, Y, Z, XTBL (NX), YTBL (NY), ZTBL (NZ),
     >   FTBL (IDIMX, IDIMY, NZ)

C     Local constants:

      REAL
     >   ONE
      PARAMETER
     >  (ONE = 1.E+0)

C     Local variables:

      INTEGER
     >   I, J, K
      REAL
     >   P, Q, R, FK, FKP1

C     Procedures:

      EXTERNAL
     >   INTERVAL

C     Execution:

      TABLE3 = 0.E+0
      IER = 0

      IF (NX .LT. 2) THEN
         IER = 1
         GO TO 99
      END IF

      IF (NY .LT. 2) THEN
         IER = 2
         GO TO 99
      END IF

      IF (NZ .LT. 2) THEN
         IER = 3
         GO TO 99
      ENDIF

      I = 1
      CALL INTERVAL (NX, XTBL, X, SIGN (ONE, XTBL (2) - XTBL (1)), I)

      P = (X - XTBL (I)) / (XTBL (I+1) - XTBL (I))

      J = 1
      CALL INTERVAL (NY, YTBL, Y, SIGN (ONE, YTBL (2) - YTBL (1)), J)

      Q = (Y - YTBL (J)) / (YTBL (J+1) - YTBL (J))

      K = 1
      CALL INTERVAL (NZ, ZTBL, Z, SIGN (ONE, ZTBL (2) - ZTBL (1)), K)

      R = (Z - ZTBL (K)) / (ZTBL (K+1) - ZTBL (K))

C     Now I and J are subscripts to the "left" of X and Y, or at X or Y
C     exactly, unless we are extrapolating, in which case the formula
C     still applies using points 1 & 2 or N-1 & N in either direction:

      FK = (ONE - P) * (ONE - Q) * FTBL (I, J, K)   +
     >                 P * (ONE - Q) * FTBL (I+1, J, K) +
     >         (ONE - P) * Q         * FTBL (I, J+1, K) +
     >                 P * Q         * FTBL (I+1, J+1, K)
      
      FKP1 = (ONE - P) * (ONE - Q) * FTBL (I, J, K+1)   +
     >                 P * (ONE - Q) * FTBL (I+1, J, K+1) +
     >         (ONE - P) * Q         * FTBL (I, J+1, K+1) +
     >                 P * Q         * FTBL (I+1, J+1, K+1)
      
      TABLE3 = (1 - R) * FK + R * FKP1

   99 RETURN
      END
