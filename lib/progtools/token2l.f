C+------------------------------------------------------------------------------
C
      SUBROUTINE TOKEN2L (STRING, SEPS, NUMBER, LIST)
C
C        TOKEN2L is a variation of TOKEN2 that does NOT convert to uppercase.
C        Most of the description from TOKEN2 remains the same.
C
C     Description and usage:
C
C           An aid to parsing input data.  The individual "tokens" in a
C        character string are isolated, converted to uppercase, and stored
C        in an array.  Here, a token is a group of significant, contiguous
C        characters.  The separators (e.g., blank and comma) are supplied
C        by the calling program.  See SCAN2 for details.
C
C           Processing continues until the requested number of tokens have
C        been found or the end of the input string is reached.
C
C     Arguments:
C
C        Name    Dimension  Type  I/O/S  Description
C        STRING    *         C    I      Input character string to be analyzed.
C        SEPS      *         C    I      String of separators, left-justified.
C        NUMBER              I    I/O    Number of tokens requested (input) and
C                                        found (output).
C        LIST    NUMBER      C      O    Array of tokens, with case unchanged.
C
C     Procedures:
C
C        Name    Description
C        SCAN2   Finds positions of the first and last significant characters.
C
C     Author:  Robert Kennelly, Informatics General Corporation (TOKEN2).
C
C     History:
C
C        12 Mar. 1986   R.A.Kennelly  TOKEN2 variation of TOKENS to take
C                                     advantage of SCAN2's added flexibility
C                                     over SCANNR.
C        06 Sep. 2014   D.A.Saunders  TOKEN2L allows lowercase characters in
C                                     the returned tokens, unlike TOKEN2.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      CHARACTER (1), PARAMETER :: BLANK = ' '

C     Local variables:

      INTEGER   :: COUNT, FIRST, I, LAST, MARK, NUMBER
      CHARACTER :: STRING * (*), SEPS * (*), LIST (NUMBER) * (*)

C     Procedures:

      EXTERNAL :: SCAN2

C     Execution:

      FIRST = 1
      LAST  = LEN (STRING)
      COUNT = 0

      DO  ! Until no more tokens or the indicated number have been found.

         CALL SCAN2 (STRING, SEPS, FIRST, LAST, MARK)  ! Isolate next token

         IF (MARK == 0) EXIT

         COUNT = COUNT + 1
         LIST (COUNT) = STRING (FIRST : MARK)  ! Do NOT convert to uppercase
         IF (COUNT == NUMBER) EXIT

         FIRST = MARK + 2
      END DO

C     Fill the rest of LIST with blanks and set NUMBER for output.

      LIST (COUNT + 1 : NUMBER) = BLANK
      NUMBER = COUNT

      END SUBROUTINE TOKEN2L
