C+----------------------------------------------------------------------
C
      SUBROUTINE PLTYTL (XLOC, XLENMX, YTXT, HITTXT, TEXT)
C
C  PURPOSE:
C     PLTYTL centers a given character string at a specified vertical
C     location.  It reduces character size if necessary to get a fit.
C     Leading and trailing blanks are suppressed; embedded blanks are
C     left intact.  The output is in the current default font (upper
C     and lower case standard unless set otherwise outside this routine).
C
C  METHOD:    
C     *  Find first and last non-blank characters in TEXT.
C     *  Determine length of non-blank text. (Uses DISSPLA.)
C     *  If too long, adjust character size, redetermine length.
C     *  Determine starting point of centered TEXT.
C     *  Output TEXT.
C
C  ARGUMENTS:
C    ARG    TYPE  I/O/S   DIM   DESCRIPTION
C   XLOC     R      I      -    Left hand starting point in inches
C                               from the current physical origin for
C                               XLENMX.
C   XLENMX   R      I      -    Horizontal length of "box" allowed for
C                               output TEXT.  The height of the text
C                               is reduced until the length of the text
C                               does not exceed 0.9 * XLENMX.
C   YTXT     R      I      -    Vertical position on plot at which to
C                               output TEXT.
C   HITTXT   R      I      -    Character height to be used for TEXT.
C                               (A smaller value may have to be used.)
C   TEXT     C*     I      -    Text to be centered.
C
C  ENVIRONMENT:  VAX/VMS; FORTRAN 77; CA-DISSPLA
C
C  HISTORY:
C  04/22/82    PJT    Original design and coding.
C  11/22/82    RGL    Uses character variables and VAX extensions now.
C  01/19/83    RCL    Adapted as PLTYTL from CENTER. No longer sets
C                     TRIPLX character font.
C  09/01/89    MDW    Updated to run on DISSPLA version 11.0.
C  02/01/91    DAS    Removed the DO WHILEs and END DOs.
C
C  AUTHORS: P.TROSIN, R.LANGHI, R.LEFKOWITZ, INFORMATICS GENERAL
C
C-----------------------------------------------------------------------

C     Arguments:

      REAL       XLOC, XLENMX, YTXT, HITTXT
      CHARACTER  TEXT*(*)

C     Local constants:

      CHARACTER  BLANK*1
      PARAMETER (BLANK = ' ')

C  *  Find the first and last non-blank characters in TEXT:

      IBGN = 1
      IEND = LEN( TEXT )
  100 IF (TEXT (IBGN : IBGN) .EQ. BLANK .AND. IBGN .LT. IEND) THEN
         IBGN = IBGN + 1
         GO TO 100
      END IF

  200 IF (TEXT (IEND : IEND) .EQ. BLANK .AND. IEND .GT. IBGN) THEN
         IEND = IEND - 1
         GO TO 200
      END IF

      NCHARS = IEND - IBGN + 1

C  *  Adjust character height if text won't fit in space provided:

      XLEN = 1000.
      HITE = HITTXT / 0.9
  300 IF (XLEN .GT. 0.9 * XLENMX .AND. HITE .GT. 0.02) THEN
         HITE = 0.9 * HITE
         CALL HEIGHT (HITE)
         XLEN = XMESS (TEXT (IBGN : IEND), NCHARS)
         GO TO 300
      END IF

C  *  Output text:

      CALL MESSAG (TEXT (IBGN : IEND), NCHARS,
     >             XLOC + 0.5 * (XLENMX - XLEN), YTXT)
      CALL RESET ('HEIGHT')

      RETURN
      END
