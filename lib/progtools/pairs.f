C+----------------------------------------------------------------------
C
      SUBROUTINE PAIRS (STRING, NUMBER, KEYS, VALUES)
C
C
C     Description and usage:
C
C           An aid to parsing input data in keyword format.  A character
C        string is broken up into keywords and corresponding values.  The
C        "keywords" and "values" are simply alternating fields of significant,
C        contiguous characters, converted to uppercase.  The following
C        characters are NON-significant, and hence may serve as separators:
C           (a)  blanks,
C           (b)  horizontal tabs,
C           (c)  commas,
C           (d)  colons, and
C           (e)  equal signs.
C        See SCANNR for details.  Processing continues until the requested
C        number of (KEYS, VALUES) pairs have been found or the end of the
C        input string is reached.  A trailing key without a corresponding
C        value will be recognized and included in the count of pairs returned.
C        (Entries requested, but not found, will be set to the blank string.)
C
C
C     Parameters:
C
C        Name    Dimension  Type  I/O/S  Description
C        STRING              C    I      Input string.
C        NUMBER              I    I/O    Number of (KEYS,VALUES) pairs
C                                        requested (input) and found (output).
C        KEYS    NUMBER      C      O    Array of keywords, left-justified and
C                                        converted to uppercase.  "Leftover"
C                                        array entries are filled with blanks.
C        VALUES  NUMBER      C      O    Array of values in the form of
C                                        character strings, left-justified and
C                                        converted to uppercase.  "Leftover"
C                                        array entries are filled with blanks.
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
C        21 Feb. 1984  RAK/DAS    Initial design and coding based on TOKENS.
C        16 Mar. 1984    RAK      Cosmetic changes.
C         1 July 1985    RAK      Don't change input string to uppercase,
C                                 change keys/values individually instead.
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

      LOGICAL
     >   ODD
      INTEGER
     >   COUNT, I, FIRST, LAST, MARK, NUMBER
      CHARACTER
     >   KEYS (NUMBER) * (*), STRING * (*), VALUES (NUMBER) * (*)

C     Procedures.

      EXTERNAL
     >   SCANNR, UPCASE


C     Executable statements.
C     ----------------------

C     WHILE there are tokens left, loop UNTIL enough have been found.

      FIRST = 1
      LAST = LEN (STRING)

      ODD = .FALSE.
      COUNT = 0
   10 CONTINUE
         CALL SCANNR (STRING, FIRST, LAST, MARK)
         IF (LAST .GT. 0) THEN

C           The tokens go alternately into KEYS and VALUES.

            ODD = .NOT.ODD         
            IF (ODD) THEN
               COUNT = COUNT + 1
               KEYS (COUNT) = STRING (FIRST:MARK)
               CALL UPCASE (KEYS (COUNT))
            ELSE
               VALUES (COUNT) = STRING (FIRST:MARK)
               CALL UPCASE (VALUES (COUNT))
            END IF

            FIRST = MARK + 2
            IF (COUNT .LT. NUMBER .OR.
     >         (COUNT .EQ. NUMBER .AND. ODD)) GO TO 10

         END IF


C     Fill the rest of the output arrays with blanks and set NUMBER for output.

      IF (ODD) VALUES (COUNT) = BLANK
      DO 20 I = COUNT + 1, NUMBER
         KEYS (I) = BLANK
         VALUES (I) = BLANK
   20 CONTINUE

      NUMBER = COUNT


C     Termination.
C     ------------

      RETURN
      END
