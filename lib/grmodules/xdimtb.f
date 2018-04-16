C+------------------------------------------------------------------------------
C
      REAL FUNCTION XDIMTB (TEXT, LENTXT, LINE1, LINE2)
C
C ONE LINER:  Calculates length of a text block in inches. (Level 1 - 3)
C
C PURPOSE:
C    XDIMTB finds the length in inches of an imaginary rectangle enclosing a
C    block of text.  This routine may be used as an alternative to the DISSPLA
C    subroutine XSTORY, for handling data in character array format.
C
C    The character size and font are assumed to be set outside the routine
C    and hence constant within the indicated lines of text.
C
C    XDIMTB was introduced as when it was realized that the SMDLIB utility
C    TXTLEG needed such a function, both prior to calling TXTLEG (for 
C    positioning the legend) and within TXTLEG (for the optional box).  The
C    translation of TXTLEG to DISSPLA required the similar translation of 
C    XDIMTB.
C
C ARGUMENTS:
C     ARG      DIM     TYPE      I/O/S      DESCRIPTION
C     TEXT    (*) * (*)   C        I       Character array containing text.
C                                          (Terminating '$'s may be required,
C                                          depending on the usage of 
C                                          LENTXT (*).)  Elements LINE1:LINE2 
C                                          will be displayed in the current 
C                                          color (as opposed to the color 
C                                          associated with each line/symbol).
C     LENTXT     (*)      I        I       Array containing number of characters
C                                          in text line(s), if known.  (No 
C                                          trailing '$'s are needed in this 
C                                          case.)  Otherwise, pass 100 as the 
C                                          first (and only) value, and all lines
C                                          of text will be self counted and all
C                                          lines of text will be self counted
C                                          (requiring the trailing '$'s).
C    LINE1      -        I         I        First line of text to process.
C    LINE2      -        I         I        Last line of text to process.
C
C METHOD:
C    The size of the block is the length of the longest line in the block.
C    DISSPLA's XMESS utility, with its self counting option, is used for
C    finding the length of a string.
C
C ERROR HANDLING:  No error message if called at wrong level.
C
C EXTERNAL REFERENCES:
C    XMESS      Finds length of text string in inches
C
C ENVIRONMENT:  VMS/VAX; FORTRAN 77
C
C HISTORY:
C    04/19/89    M.D. Wong        Initial design and implementation for SMDLIB
C    06/07/89    M.D. Wong        Adapted for DISSPLA from SMDLIB.  Error 
C                                 handling removed (since internal common blocks
C                                 are avoided).
C    08/29/89    M.D. Wong        Updated to run on DISSPLA version 11.0.
C                                 (%REF taken out of call to function XMESS.)
C    02/26/90    M.D. Wong        Added LENTXT to argument list to avoid 
C                                 reliance on DISSPLA self counting option.
C    04/12/90    M.D. Wong        Enabled 100 to be passed as first and only
C                                 element of the LENTXT array to enable self 
C                                 counting option.
C
C AUTHOR:  Michael Wong, Sterling Software, Palo Alto, CA
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments
C     ---------

      CHARACTER TEXT (*) * (*)
      INTEGER   LENTXT (*), LINE1, LINE2

C     Local variables
C     ---------------

      INTEGER   I, LENGTH

C     Procedures
C     ----------

      REAL      XMESS
      EXTERNAL  XMESS

C     Execution
C     ---------

      XDIMTB = 0.
      DO 100, I = LINE1, LINE2
         IF (LENTXT (1) .EQ. 100) THEN
            LENGTH = LENTXT (1)
         ELSE
            LENGTH = LENTXT (I)
         END IF
         XDIMTB = MAX (XDIMTB, XMESS (TEXT (I), LENGTH))
  100 CONTINUE

      RETURN
      END
