C+----------------------------------------------------------------------
C
      SUBROUTINE SCAN3 (STRING, SEPS, FIRST, LAST, MARK)
C
C     One-liner:
C
C        Backward-search variant of SCAN2 - find last token in a (sub)string
C
C     Description and usage:
C
C           Looks for the LAST non-blank field ("token") in a (sub)string by
C        searching from right to left.  The fields are of arbitrary length
C        and separated by any of a set of user-specified separators (e.g.,
C        blanks or commas).  This routine may be conveniently used within
C        a loop to process an entire line of text BACKwards.  (However, the
C        need for "TOKEN3" to go with SCAN3 appears limited - just use SCAN3
C        on the substring defined by FIRST and LAST = MARK-2 if you want further
C        tokens in reverse order.)
C
C           To clarify:  Given STRING and positions FIRST, LAST,
C
C        SCAN3 locates RIGHT-most token as STRING (MARK : LAST)  (LAST updated);
C        SCAN2    "    LEFT   "    "    "  STRING (FIRST : MARK) (FIRST  "   ").
C
C           In both cases, MARK = 0 if no token was found.
C
C     Arguments:
C
C        Name    Dimension  Type  I/O/S  Description
C        STRING    *         C    I      Text string containing data to be
C                                        scanned.
C        SEPS      *         C    I      String containing desired separators.
C                                        Each character in SEPS is treated
C                                        as a possible token delimiter.
C        FIRST               I    I      Index of first character of STRING
C                                        of interest.  Not updated here.
C        LAST                I    I/O    Input as last character of STRING
C                                        of interest.  Output is right-most
C                                        end of right-most token found in
C                                        specified (sub)string.  Unchanged if
C                                        no token was found (MARK = 0).
C        MARK                I      O    Points to left-most character of the
C                                        first token found searching BACKward
C                                        in the specified substring. MARK = 0
C                                        means that no token was found.
C
C
C     Environment:  VAX/VMS FORTRAN 77
C
C
C     Notes:
C
C        (1)  IMPLICIT NONE is non-standard.
C
C        (2)  Error checking:  None.  The calling program is expected
C             to set FIRST and LAST correctly on the first call for a
C             given STRING (probably at 1 and LEN (STRING)) and update
C             LAST as MARK-2 for further searches of the same string.
C
C
C     Author:  Robert Kennelly/David Saunders, Sterling Software, Palo Alto.
C
C
C     Development history:
C
C      4 Mar. 1986   RAK   SCAN2 derived from SCANNR.
C     21 Apr. 1987   DAS   SCAN3 derived from SCAN2 for the backward search
C                          case.  Simplified with one-or-two-trailing tokens
C                          per line in mind.
C
C-----------------------------------------------------------------------


      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   FIRST, LAST, MARK
      CHARACTER
     >   SEPS * (*), STRING * (*)

C     Local variables.

      INTEGER
     >   I, J, TAIL
      LOGICAL
     >   FOUND


C     Execution.
C     ----------

      MARK = 0
      TAIL = LAST
      FOUND = .FALSE.

C     Look at each character in STRING (FIRST : TAIL) from right to left
C     until a token is found or until FIRST is encountered.  Drop through
C     with MARK at 0 if TAIL < FIRST.

      DO 30, I = TAIL, FIRST, -1

C        Is the character a separator?

         DO 10, J = 1, LEN (SEPS)
            IF (STRING (I : I) .EQ. SEPS (J : J)) GO TO 20
   10    CONTINUE

C           Not a separator.  Check for the beginning of the right-most token.

            IF (.NOT. FOUND) THEN
               LAST = I
               FOUND = .TRUE.
            END IF
            GO TO 30

   20    CONTINUE

C           We found a separator.

            IF (FOUND) THEN

C              We just passed the "leading edge" of the desired token.
C              Bail out with a successful search.

               MARK = I + 1
               GO TO 99
            END IF
   30 CONTINUE

C     We reached the FIRST character.  Either it is part of a token, or
C     MARK is still zero.

      IF (FOUND) MARK = FIRST


C     Termination.

   99 RETURN
      END
