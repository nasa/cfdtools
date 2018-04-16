C+----------------------------------------------------------------------
C
      SUBROUTINE POWER2 (L, M, LEFTOVER)
C
C PURPOSE: POWER2 determines whether given integer L is a positive power
C          of 2.  It finds M such that L = 2 ** M + LEFTOVER, where normal
C          usage will be looking for LEFTOVER = 0 and M > 0.  The degenerate
C          case L = 1 is indicated by M = 0 = LEFTOVER.  L <= 0 is considered
C          to be an error, indicated by M = -1 = LEFTOVER.
C
C METHOD:  One might be tempted to make use of the relation
C
C             M   =   log L   =   log L / log 2
C                        2           e       e
C
C          with some tolerance to allow for round-off error.
C          But a handful of divides by 2 seems preferable to use of
C          logarithms, and is the approach used here.
C
C ARGUMENTS:
C   ARG       TYPE I/O/S DESCRIPTION 
C    L         I     I   Integer in question - is it a power of 2?
C    M         I     O   M > 0  means L = 2 ** M + possible remainder;
C                        M = 0  means L = 1;
C                        M = -1 means L <= 0.
C  LEFTOVER    I     O   LEFTOVER = 0 means L is a power of 2 (possibly 2 ** 0)
C                        LEFTOVER > 0 means L is not a power of 2 - in fact,
C                                     L = 2 ** M + LEFTOVER;
C                        LEFTOVER = -1 if L <= 0.
C
C ENVIRONMENT:  VAX/VMS, FORTRAN 77.
C
C HISTORY:
C   05/14/86  Initial implementation (more explanation needed than expected)
C
C AUTHOR:  David Saunders, Sterling Software, Palo Alto, CA.
C
C-----------------------------------------------------------------------

C     Arguments:

      INTEGER
     >   L, M, LEFTOVER
      
C     Local variables:

      INTEGER
     >   TARGET

C     Execution:

      TARGET = L
      IF (TARGET .LE. 0) THEN
         M = -1
         LEFTOVER = -1
         GO TO 99
      END IF

      M = 0
   20 CONTINUE
         TARGET = TARGET / 2
         IF (TARGET .GT. 0) THEN
            M = M + 1
            GO TO 20
         END IF
         
      LEFTOVER = L - 2 ** M

   99 RETURN
      END
