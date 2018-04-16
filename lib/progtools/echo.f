C+----------------------------------------------------------------------
C
      SUBROUTINE ECHO (LUNINP, LUNOUT, LUNSCR)
C
C  PURPOSE: Copies normal input stream exactly to normal output file,
C           and creates a scratch file on disk to be used as input by
C           the calling program. (This bypasses possible difficulties
C           with rewinding the normal input file and means the echoed
C           data is unaffected by other activities requiring  inputs,
C           such as "UPDATE".)
C
C  ENVIRONMENT: CRAY/COS, VAX/VMS - FORTRAN.
C
C  ARGUMENTS:
C  VAR    TYPE I/O/S    DIM    DESCRIPTION
C  LUNINP  I     I       -     Logical unit number for input file
C  LUNOUT  I     I       -     Logical unit number for printed output
C  LUNSCR  I     I       -     Logical unit number for scratch file
C
C  LOCAL VARIABLES:
C  VAR    TYPE       DIM       DESCRIPTION
C  CARD    R         132       Storage for 1 record of input data
C
C  FILES USED:  See arguments.
C
C  AUTHOR:  Carol J. Green, Informatics General, Palo Alto, CA.
C
C  REVISIONS:
C    DATE      PERSON    STATEMENT OF CHANGES
C  04/26/82     CJG      Original coding
C  01/03/84     RAK      Read 132 columns of data, squeeze out trailing
C                        blanks.
C  02/13/84     RAK      Eliminated initial page feed.
C
C-----------------------------------------------------------------------


C     Variable declarations.
C     ----------------------

      IMPLICIT NONE

C     Parameters.

      INTEGER
     >   WIDTH
      CHARACTER
     >   BLANK
      PARAMETER
     >   (BLANK = ' ', WIDTH = 132)

C     Variables.

      INTEGER
     >   LAST, LUNINP, LUNOUT, LUNSCR
      CHARACTER
     >   LINE * (WIDTH)


C     Executable statements.
C     ----------------------

   10 CONTINUE
         READ  (LUNINP, '(A)', END = 40) LINE

C        Find the last non-blank (if any).

         DO 20 LAST = WIDTH, 1, -1
            IF (LINE (LAST : LAST) .NE. BLANK) GO TO 30
   20    CONTINUE
         LAST = 1

C        Print/output a line.

   30    CONTINUE
         WRITE (LUNOUT, '(1H , (A))') LINE (1 : LAST)
         WRITE (LUNSCR, '(A)') LINE (1 : LAST)
         GO TO 10


   40 CONTINUE
      REWIND (LUNSCR)
      RETURN
      END
