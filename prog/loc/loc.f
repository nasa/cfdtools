C+----------------------------------------------------------------------
C
      PROGRAM LOC
C
C     Description and usage:
C
C           A quick-and-dirty filter for counting lines of FORTRAN code
C        in an ASCII file.  We count everything EXCEPT blanks and lines
C        which begin with 'C' in column 1.  Total lines and lines of code
C        are reported to the terminal, and the user is prompted for another
C        source file.
C
C
C     Input format:
C
C           Any ASCII file.  Only the first 72 columns are examined.
C
C
C     External files:
C
C        Unit    I/O/S  Description
C         4      I      Source code file to be filtered.
C         5      I      User keyboard.
C         6        O    User terminal screen.
C
C
C     Environment:  Digital VAX-11/780 VMS/V4.1 FORTRAN.
C
C
C     Notes:
C
C        (1)  IMPLICIT NONE is non-standard, and some of the variable
C             names are longer than six characters.
C
C        (2)  We don't trace INCLUDE files.
C
C        (3)  VT100 escape codes are ANSI standard (I hope).
C
C        (4)  Use of CHAR function in PARAMETER statement is non-standard.
C
C
C     Author:  Robert Kennelly, Informatics General Corporation.
C
C
C     Development history:
C
C         1 Oct. 1985    RAK    Initial design and coding.
C         5 Mar. 1986    RAK    Added READONLY keyword in the OPEN to
C                               facilitate snooping in other people's
C                               directories.  "Clarity coefficient"
C                               displayed.  Playing with VT100 emphasis
C                               features for main results and errors.
C                               Added explicit CLOSE (was previously
C                               unable to look at same file twice).
C                               Use OPENER to prompt for and open the
C                               input file.
C        30 May  1986    DAS    Took advantage of revised OPENER (which
C                               now returns if filename is left blank).
C        20 May  1988    DAS    Introduced GETLINE to deal with trailing
C                               '!'-type comments properly; allowed '*'
C                               as well as 'C', 'c' for column 1 comments.
C        28 Jan. 2000    DAS    I5 format didn't allow for really large
C                               line counts.
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT NONE

C     Parameters.

      INTEGER
     >   KEYBRD, SCREEN, SOURCE
      PARAMETER
     >   (KEYBRD = 5,
     >    SCREEN = 6,
     >    SOURCE = 4)

      CHARACTER
     >   BLANK * 1, BLINK * 4, BOLD * 4, ESC, INVERSE * 4, RESET * 4
      PARAMETER
     >   (ESC     = CHAR (27),
     >    RESET   = ESC // '[0m',
     >    BOLD    = ESC // '[1m',
     >    BLINK   = ESC // '[5m',
     >    INVERSE = ESC // '[7m',
     >    BLANK   = ' ')

C     Variables.

      LOGICAL
     >   GETLOST
      INTEGER
     >   GROSS, IOSTAT, LAST, NET
      CHARACTER
     >   BUFFER * 72

C     Procedures.

      EXTERNAL
     >   GETLINE, OPENER

C     Execution.
C     ----------

C     Prompt for a filename to examine (we only pass in a substring of
C     BUFFER so that OPENER/READC will leave the cursor on the prompt line).
C     Note that OPENER provides termination if ^Z is entered.

   10 CONTINUE
         WRITE (SCREEN, 1000)
         BUFFER = BLANK
         CALL OPENER (SCREEN, BOLD // 'Source file? ' // RESET,
     >      KEYBRD, BUFFER (1 : 59), SOURCE, 'OLD')

         IF (BUFFER .EQ. BLANK) STOP ' '

C        Count total lines, and those which are FORTRAN.
C        -----------------------------------------------

         GROSS = 0
         NET = 0
   20    CONTINUE

            CALL GETLINE (SOURCE, '!', BUFFER, LAST, IOSTAT)
            IF (IOSTAT .LT. 0) GO TO 30                     ! End-of-file

            GROSS = GROSS + 1
            IF (LAST .GT. 0 .AND.
     >          BUFFER (1 : 1) .NE. 'C' .AND.
     >          BUFFER (1 : 1) .NE. 'c' .AND.
     >          BUFFER (1 : 1) .NE. '*') NET = NET + 1
            GO TO 20

   30    CONTINUE


C        Report findings.
C        ----------------

         IF (GROSS .GT. 0) THEN
            WRITE (SCREEN, 1000)
            WRITE (SCREEN, 1010) 'Total lines   =', GROSS
            WRITE (SCREEN, 1010) 'FORTRAN lines =', NET

C           Derived quantities.

            WRITE (SCREEN, 1000)
            WRITE (SCREEN, 1020)
     >         INVERSE,
     >         REAL (NET) / REAL (GROSS) * 100.0,
     >         ' = Purity (%) ' //
     >         RESET // ' (#FORTRAN / #total * 100)'

            IF (NET .GT. 0) THEN
               WRITE (SCREEN, 1020)
     >            INVERSE,
     >            REAL (GROSS - NET) / REAL (NET),
     >            ' = Clarity    ' //
     >            RESET // ' (#comments / #FORTRAN)'
            ELSE

C              No code - gripe.

               WRITE (SCREEN, 1000)
               WRITE (SCREEN, 1030)
     >            BLINK //
     >            '!@#$%  No FORTRAN:  clarity high, content low!' //
     >            RESET
            END IF
         ELSE

C           No source - gripe.

            WRITE (SCREEN, 1000)
            WRITE (SCREEN, 1030)
     >         BLINK //
     >         '!@#$%  Source file empty:  very pure indeed!' //
     >         RESET
         END IF

         WRITE (SCREEN, 1000)
         CLOSE (UNIT = SOURCE)
         GO TO 10


C     Formats.
C     --------

 1000 FORMAT (1X)
 1010 FORMAT (3X, A, I7)
 1020 FORMAT (3X, A, F6.2, A)
 1030 FORMAT (3X, A)

      END
