C+------------------------------------------------------------------------------
C
      SUBROUTINE SCAN4 (STRING, DELIMITER, FIRST, LAST, MARK)
C
C
C     Description and usage:
C
C           This variant of SCAN2 handles the case of blanks embedded in tokens
C        defined by some other delimiter or pair of delimiters.  For instance,
C        the variable names in an ASCII Tecplot file might look like this:
C
C                      VARIABLES = "X, m" "Y, m" "Z, m"
C
C        SCAN2 with '", ' as delimiters would find 6 names here, not 3, so a
C        different strategy is required.
C
C           Like SCAN2, SCAN4 examines substring STRING (FIRST : LAST), which
C        may of course be the entire string (in which case just call SCAN4
C        with FIRST = 1 and LAST = LEN (STRING) ).  The indices returned are
C        relative to STRING itself, not the substring.
C
C           Certain leading delimiters are assumed to be paired with natural
C        trailing delimiters as follows:
C
C                                 ()  []  {}
C
C     Arguments:
C
C        Name    Dimension  Type  I/O/S  Description
C        STRING      *       C    I      Text string containing data to be
C                                        scanned.
C        DELIMITER           C    I      Single character a pair of which
C                                        defines the next token in STRING.
C                                        See also possible pairing above.
C        FIRST               I    I/O    Index of first character of interest
C                                        in STRING.  FIRST >= 1 is assumed on
C                                        input; output is the index of the
C                                        start of the next token, if MARK is
C                                        not set to 0.
C        LAST                I    I      Index of last character of interest
C                                        in STRING.  LAST <= LEN (STRING) is
C                                        assumed.  LAST is NOT updated.
C                                        Using LEN_TRIM (STRING) is normally
C                                        appropriate under Fortran 90.
C        MARK                I      O    Points to the end of the first token
C                                        found in the specified portion of
C                                        STRING.  Use FIRST = MARK + 2 for
C                                        further searches in the same string.
C                                        MARK is set to 0 if no token was found.
C                                        MARK < 0 signals the special case of
C                                        a pair of adjacent delimiters, which
C                                        the application may choose to treat
C                                        as a blank token, perhaps.  The user
C                                        should reassign MARK = -MARK to
C                                        recover the normal pointer.
C
C     Notes:
C
C        (1)  FIRST >= 1 and LAST <= LEN (STRING) are assumed upon entry,
C             because these are the obvious bounds needed by the calling
C             program on the first scan of a given string - no need to
C             evaluate LEN (STRING) here with every call, while use of
C             MIN and MAX tends to disguise programming errors.
C
C        (2)  Entering with FIRST > LAST is an eventual normal consequence
C             of setting FIRST = MARK + 2 for further scans of the same string.
C             The zero DO loop iteration count feature of the language is
C             taken advantage of here.
C
C        (5)  Zeroing all of FIRST, LAST, and MARK (rather than just MARK)
C             when no token is found may, in retrospect, be undesirable in
C             that information is lost (esp. LAST).  Therefore, unlike SCAN2,
C             SCAN4 does NOT zero LAST if no token is found.  In fact it
C             does not change it at all now, because different kinds of
C             delimiters may be used for different keywords in the same
C             string.
C
C
C     History:
C
C         4 Mar 1986  Robert    SCAN2 variation of SCANNR, which is hard-coded
C                     Kennelly  (for historical reasons and a modest speed
C                     Sterling  advantage) for separators BLANK, TAB,
C                     Software  COMMA, COLON, and EQUAL.
C
C         5 May  1988   DAS     SCAN2's reverse search is used to find LAST;
C                    Sterling   MAX, MIN, and LEN have been eliminated.
C
C         5 May  2005   DAS     SCAN4 adaptation of SCAN2.  Note the MARK < 0
C                               possibility.
C
C        29 July 2006   DAS     Do NOT change LAST; extend ".." case to any of
C                               (..), [..], and {..}.
C
C     Author: David Saunders, ELORET/NASA Ames Research Ctr., Moffett Field, CA
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   FIRST, LAST, MARK
      CHARACTER
     >   DELIMITER * 1, STRING * (*)

C     Local variables.

      INTEGER
     >   HEAD, I, TAIL
      LOGICAL
     >   FOUND
      CHARACTER
     >   DELIM2 * 1, TARGET * 1

C     Execution.

      SELECT CASE (DELIMITER)
      CASE ('(')
         DELIM2 = ')'
      CASE ('[')
         DELIM2 = ']'
      CASE ('{')
         DELIM2 = '}'
      CASE DEFAULT
         DELIM2 = DELIMITER
      END SELECT

      TARGET = DELIMITER
      HEAD   = FIRST
      TAIL   = LAST

      FIRST  = 0
      MARK   = 0
      FOUND  = .FALSE.

C     Look at each character in STRING (HEAD : TAIL) from left to right
C     until a token is found.  (Drop through if LAST > FIRST.)

      DO I = HEAD, TAIL

         IF (STRING (I : I) == TARGET) THEN

            IF (.NOT. FOUND) THEN ! We've found the start of the first token
               FOUND  = .TRUE.
               TARGET = DELIM2
               FIRST  = I + 1
            ELSE                  ! We've found the end of the first token
               MARK   = I - 1
               EXIT
            END IF

         END IF ! Else keep searching for a delimiter

      END DO

C     An adjacent pair of delimiters might be considered to mean a blank token.
C     Allow the application to detect this:

      IF (FIRST > HEAD .AND. MARK == FIRST - 1) MARK = -MARK

      END SUBROUTINE SCAN4
