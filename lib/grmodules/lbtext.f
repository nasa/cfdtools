C+------------------------------------------------------------------------------
C
      SUBROUTINE LBTEXT (TEXT, LENTXT, LINE1, LINE2, DELTAY, XCORNR, 
     >                   YCORNR, HITE)
C
C ONE LINER:
C     Displays a left-justified text block using DISSPLA utilities. (Level 2-3)
C
C PURPOSE:
C     LBTEXT composes a left-justified text block and displays it at a specified
C     position.  Character height and line spacing are assumed constant within
C     the indicated lines of text.
C
C     LBTEXT was developed as part of the TXTLEG alternative to the original
C     SMDLIB legend utility (TXTBLK).  This is the DISSPLA version. 
C
C ARGUMENTS:
C      ARG       DIM     TYPE    I/O/S     DESCRIPTION
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
C     LINE1       -       I        I       First line of text to display.
C     LINE2       -       I        I       Last line of text to display.
C     DELTAY      -       R        I       Space between lines in inches.
C     XCORNR      -       R        I       Horizontal position of lower left
C                                          hand corner of text block in inches.
C     YCORNR      -       R        I       Vertical position of lower left hand
C                                          corner of text block in inches.
C     HITE        -       R        I       Text character height in inches.
C
C
C ERROR HANDLING:  No check if LBTEXT is called at the wrong level.
C
C EXTERNAL REFERENCES:  DISSPLA utilities.
C
C ENVIRONMENT:  VAX/VMS; FORTRAN 77
C
C HISTORY:
C     04/17/89    M.D. Wong      Initial design and coding for SMDLIB.
C     06/06/89    M.D. Wong      Adapted from SMDLIB to DISSPLA.  
C     08/29/89    M.D. Wong      Updated to run on DISSPLA version 11.0
C                                (%REF taken out of call to LBTEXT.)
C     02/26/90    M.D. Wong      LENTXT added to argument list to avoid 
C                                reliance on DISSPLA's self counting option.
C     04/12/90    M.D. Wong      Enabled 100 to be passed as first and only
C                                element of the LENTXT array to envoke the
C                                self-counting option.
C
C AUTHOR:  Michael Wong, Sterling Software, Palo Alto, CA.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments
C     ---------

      CHARACTER  TEXT (*) * (*)
      INTEGER    LENTXT (*), LINE1, LINE2
      REAL       DELTAY, HITE, XCORNR, YCORNR

C     Local Variables
C     ---------------

      INTEGER    I, LENGTH
      REAL       GAP, Y 

C     Execution
C     ---------

C     Find distance between lines.

      GAP = DELTAY + HITE 

C     Start with upper left hand corner of text block.

      Y = YCORNR + (LINE2 - LINE1) * GAP

      DO 200 I = LINE1, LINE2
         IF (LENTXT (1) .EQ. 100) THEN
            LENGTH = LENTXT (1)
         ELSE
            LENGTH = LENTXT (I)
         END IF
         CALL MESSAG (TEXT (I), LENGTH, XCORNR, Y) 
         Y = Y - GAP
  200 CONTINUE

      RETURN
      END
