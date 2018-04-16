C+----------------------------------------------------------------------
C
      SUBROUTINE TOKENS (STRING, NUMBER, LIST)
C
C
C     Description and usage:
C
C           An aid to parsing input data.  The individual "tokens" in a
C        character string are isolated, converted to uppercase, and stored
C        in an array.  Here, a token is a group of significant, contiguous
C        characters.  The following are NON-significant, and hence may
C        serve as separators:  blanks, horizontal tabs, commas, colons,
C        and equal signs.  See SCANNR for details.  (Use TOKEN2 and SCAN2
C        if you need a variable list of delimiters.)
C
C           Processing continues until the requested number of tokens have
C        been found or the end of the input string is reached.
C
C
C     Parameters:
C
C        Name    Dimension  Type  I/O/S  Description
C        STRING              C    I      Input character string to be analyzed.
C        NUMBER              I    I/O    Number of tokens requested (input) and
C                                        found (output).
C        LIST    NUMBER      C      O    Array of tokens, changed to uppercase.
C
C
C     External references:
C
C        Name    Description
C        SCANNR  Finds positions of the first and last significant characters.
C        UPCASE  Converts a string to uppercase.
C
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77).
C
C
C     Notes:
C
C        (1)  IMPLICIT NONE is non-standard.
C
C
C     Author:  Robert Kennelly, Informatics General Corporation.
C
C
C     Development history:
C
C        16 Jan. 1984    RAK    Initial design and coding.
C        16 Mar. 1984    RAK    Revised header to reflect full list of
C                               separators, repaired faulty WHILE clause
C                               in "10" loop.
C        18 Sep. 1984    RAK    Change elements of LIST to uppercase one
C                               at a time, leaving STRING unchanged.
C        12 Mar. 1986    RAK    Cross-referenced TOKEN2 variation.
C
C-----------------------------------------------------------------------


C     Variable declarations.
C     ----------------------

      IMPLICIT NONE

C     Parameters.

      CHARACTER
     >   BLANK
      PARAMETER
     >   (BLANK = ' ')

C     Variables.

      INTEGER
     >   COUNT, FIRST, I, LAST, MARK, NUMBER
      CHARACTER
     >   STRING * (*), LIST (NUMBER) * (*)

C     Procedures.

      EXTERNAL
     >   UPCASE, SCANNR


C     Executable statements.
C     ----------------------

C     WHILE there are tokens to find, loop UNTIL enough have been found.

      FIRST = 1
      LAST = LEN (STRING)

      COUNT = 0
   10 CONTINUE

C        Get delimiting indices of next token, if any.

         CALL SCANNR (STRING, FIRST, LAST, MARK)
         IF (LAST .GT. 0) THEN
            COUNT = COUNT + 1

C           Pass token to output string array, then change case.

            LIST (COUNT) = STRING (FIRST : MARK)
            CALL UPCASE (LIST (COUNT))
            FIRST = MARK + 2
            IF (COUNT .LT. NUMBER) GO TO 10

         END IF


C     Fill the rest of LIST with blanks and set NUMBER for output.

      DO 20 I = COUNT + 1, NUMBER
         LIST (I) = BLANK
   20 CONTINUE

      NUMBER = COUNT


C     Termination.
C     ------------

      RETURN
      END
