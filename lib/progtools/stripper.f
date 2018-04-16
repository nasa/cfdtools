C+----------------------------------------------------------------------
C
      SUBROUTINE STRIPPER (STRING, LENSTR, DELIM)
C
C
C     Description and usage:
C
C        STRIPPER scans a string for enclosing delimiters and passes back
C     the contents (possibly left-shifted one character) along with its
C     length (possibly reduced by 2).  The length of the string on input is
C     assumed to be positive (see the argument description).  If a pair of
C     delimiters is stripped out, the last two characters of the original
C     string are blanked.  Any blanks or tabs in front of the first delimiter
C     are ignored (and retained).  A returned length of zero is possible.
C
C
C        Examples:
C
C        (1)          '12345'           is returned as
C                     12345             (with length reduced by 2).
C                          ^^
C
C        (2)          " "               (a quoted blank) becomes
C                                       (3 blanks with length returned as 1).
C                     ^^^
C
C        (3)          ""                becomes
C                                       (2 blanks but length is returned as 0).
C                     ^^
C
C           STRIPPER was prompted by the need to distinguish blank or
C        apparently numeric text from missing or truly numeric data by "quoting"
C        the text.  The quotes need to be checked for and eliminated if present,
C        with a left shift in place, and a length reduced by 2 (normally).
C
C
C     Arguments:
C
C        Name    Dimension  Type  I/O/S  Description
C        STRING      *       C    I/O    Text string to be searched.  If found,
C                                        one pair of delimiters is eliminated
C                                        on output, and remaining characters
C                                        are left-shifted one position.
C        LENSTR              I    I/O    Length of text string.  LENSTR > 0 is
C                                        assumed as input since the string's
C                                        last significant character should have
C                                        been determined by a prior call to
C                                        GETLINE, SCANx, or equivalent.  If
C                                        enclosing delimiters are found, LENSTR
C                                        is reduced by two (meaning LENSTR = 0
C                                        is a possible output).
C        DELIM       *       C    I      String containing delimiters to be
C                                        searched for, probably ' or " or both.
C                                        If a matching pair is found, the search
C                                        terminates - no further delimiters are
C                                        considered.  A valid pair consists of
C                                        the character STRING (LENSTR:LENSTR)
C                                        and the first SIGNIFICANT character in
C                                        STRING.  (Leading blanks or tabs are
C                                        ignored.)
C
C     External References:
C
C        SCAN2   Used here to find the first significant character.
C
C
C     Environment:  Digital VAX-11/785 VMS FORTRAN (FORTRAN 77).
C                   IMPLICIT NONE is non-standard.
C
C     History:
C
C        14 Mar 1990    R.A.Kennelly    Initial design.
C        15 Mar 1990    M.D.Wong        Initial implementation.
C        02 Oct 1991    D.A.Saunders    Allowed LENSTR=0 upon return,
C                                       where 1 was forced originally.
C
C     Author:  Michael Wong, NASA Ames/Sterling Software, Palo Alto, CA
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   LENSTR
      CHARACTER
     >   DELIM * (*), STRING * (*)

C     Local constants.

      CHARACTER
     >   BLANK
      PARAMETER
     >  (BLANK = ' ')

C     Local variables.

      CHARACTER
     >   SEPS * 2
      INTEGER
     >   HEAD, I, J, MARK, NDELIM, TAIL

C     Execution.


      HEAD = 1
      TAIL = LENSTR
      SEPS = BLANK // CHAR (9)      ! Tab and blank are considered
                                    ! insignificant as leading characters.
      NDELIM = LEN (DELIM)

C     Compare each element in DELIM with first and last significant characters
C     in STRING.  The search is terminated if they both match.

      DO 20, J = 1, NDELIM
         IF (STRING (TAIL : TAIL) .EQ. DELIM (J : J)) THEN

C           The last character is a possible delimiter, so continue search
C           with the first significant character, located by SCAN2.

            CALL SCAN2 (STRING, SEPS, HEAD, TAIL, MARK)

            IF (STRING (HEAD : HEAD) .EQ. DELIM (J : J) .AND.
     >         HEAD .NE. TAIL) THEN        ! Guard against single delimiter.

C              Adjust string length, allowing 0 (where 1 was once forced).

               LENSTR = TAIL - 2

C              Shift left one character, in place, and blank out last two
C              characters.  (The loop does nothing in the degenerate case.)

               DO 10, I = HEAD, LENSTR
                  STRING (I : I) = STRING (I + 1 : I + 1)
   10          CONTINUE
               STRING (TAIL - 1 : TAIL) = BLANK

               GO TO 99
            END IF
         END IF
   20 CONTINUE

   99 RETURN
      END
