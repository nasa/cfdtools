C+----------------------------------------------------------------------
C
      SUBROUTINE LSTFIL (LUNIN, LUNOUT, BUFFER)
C
C  PURPOSE: LSTFIL lists the indicated file in the sense of the TYPE
C           command of VMS.  The output may be the screen or another
C           file.  The file is assumed to be sequential/formatted.
C
C  METHOD:  Assume the file is already open and rewound.    Leave it
C           open (and not rewound) upon return.    Don't use a local
C           buffer -  the application is likely to have such scratch
C           available. Trailing blanks are suppressed - a frill per-
C           haps, but why not?
C
C  ENVIRONMENT: VAX/VMS - FORTRAN.
C
C  ARGUMENTS:
C  VAR    TYPE I/O/S    DIM    DESCRIPTION
C  LUNIN   I     I       -     Logical unit number for input file
C  LUNOUT  I     I       -     Logical unit number for output copy
C  BUFFER  C*(*) S       -     Buffer for one line of file - normally
C                              C*80, or maybe C*132.
C
C  FILES USED:  See arguments.
C
C  HISTORY:
C    April 1986  DAS/RAK  Adapted from ECHO - don't want the extra disk copy.
C    July  1989    DAS    Cosmetics in preparation for Mac translation.
C
C  AUTHORS:  David Saunders, Rob Kennelly, Sterling Software, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   LUNIN, LUNOUT
      CHARACTER
     >   BUFFER * (*)

C     Local constants:

      CHARACTER
     >   BLANK
      PARAMETER
     >  (BLANK = ' ')

C     Local variables:

      INTEGER
     >   LAST, MAXLEN

C     Execution:
      
C     Put out an initial blank line - saves the calling program from doing it:

      WRITE (LUNOUT, '(A)')

      MAXLEN = LEN (BUFFER)

   10 CONTINUE
         READ  (LUNIN, '(A)', END = 40) BUFFER

C        Find the last non-blank (if any).

         DO 20 LAST = MAXLEN, 1, -1
            IF (BUFFER (LAST : LAST) .NE. BLANK) GO TO 30
   20    CONTINUE
         LAST = 1

   30    WRITE (LUNOUT, '(1X, A)') BUFFER (1 : LAST)
         GO TO 10

   40 RETURN
      END
