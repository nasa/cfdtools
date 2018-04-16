C+----------------------------------------------------------------------
C
      SUBROUTINE SCANNR (STRING, FIRST, LAST, MARK)
C
C
C     Description and usage:
C
C           Looks for significant fields ("tokens") in a string, where the
C        fields are of arbitrary length and separated by blanks, tabs, commas,
C        colons, or equal signs.  The position of the end of the first token
C        is also returned so that this routine may be conveniently used within
C        a loop to process an entire line of text.
C
C           The procedure examines a substring, STRING (FIRST : LAST), which
C        may of course be the entire string (in which case just call SCANNR
C        with FIRST = 1 and LAST = LEN (STRING) ).  The indices returned
C        are relative to STRING itself, not the substring.
C
C
C     Arguments:
C
C        Name    Dimension  Type  I/O/S  Description
C        STRING    *         C    I      Text string containing data to be
C                                        scanned.
C        FIRST               I    I/O    Index of first character of interest
C                                        in STRING.  FIRST >= 1 is assumed.
C                                        Output is index of the beginning of
C                                        the first token, or 0 if no token
C                                        was found.
C        LAST                I    I/O    Index of last character of interest
C                                        in STRING.  LAST <= LEN (STRING) is
C                                        assumed.  Output is index of the end
C                                        of the last non-separator, or 0 if no
C                                        token was found.
C        MARK                I      O    Points to the end of the first token
C                                        found in the specified portion of
C                                        STRING.  Use FIRST = MARK + 2 for
C                                        further searches in the same string.
C                                        MARK is set to 0 if no token was found.
C
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77).
C
C
C     Notes:
C
C        (1)  IMPLICIT NONE is non-standard.  Constant HT (Tab) is defined
C             in a non-standard way:  the CHAR function is not permitted
C             in a PARAMETER declaration (OK on VAX, though).  For Absoft
C             FORTRAN 77 on 68000 machines, use HT = 9.  In other cases, it
C             may be best to declare HT as a variable and assign HT = CHAR (9).
C
C        (2)  The pseudo-recursive structure of the original version has been
C             abandoned because the VAX compiler treated the SOLID statement
C             function as an in-line subroutine, with substantial penalty
C             in speed (factor of almost 3!).  The single-loop form used
C             later was almost as fast (especially for lines with only a few
C             tokens), and was more compact and easy to change if a different
C             set of delimiters was required.  However, ...
C
C        (3)  This version adjusts LAST using a BACKward search, which is
C             degenerate for all tokens after the first.  Repeated forward
C             scanning to locate LAST (in a calling program's loop over tokens)
C             amounts to an O(4N**2) operation count for N tokens in a string
C             versus the O(2N) that it should be.  (The 2 is in there because
C             the non-token fields are as numerous as the tokens and take a
C             similar time to scan.)  The price paid is some repeated code.
C
C        (4)  FIRST >= 1 and LAST <= LEN (STRING) are assumed upon entry,
C             because these are the obvious bounds needed by the calling
C             program on the first scan of a given string - no need to
C             evaluate LEN (STRING) here with every call, while use of
C             MIN and MAX tends to disguise programming errors.
C
C        (5)  Entering with FIRST > LAST is an eventual normal consequence
C             of setting FIRST = MARK + 2 for further scans of the same string.
C             The zero DO loop iteration count feature of the language is
C             taken advantage of here - the result is simply to drop
C             out the bottom with FIRST = LAST = MARK = 0.
C
C        (6)  Zeroing all of FIRST, LAST, and MARK (rather than just MARK)
C             when no token is found may, in retrospect, be undesirable in
C             that information is lost (esp. LAST), but there are too many
C             applications to change this now.
C
C        (7)  The variety of separators recognized limits the usefulness of
C             this routine somewhat.  The intent is to facilitate handling
C             such tokens as keywords or numerical values.  In other
C             applications, it might be necessary for ALL printing characters
C             to be significant.  Use SCAN2 (from the same author) if user-
C             supplied delimiters are appropriate.
C
C        (8)  Note that "null" characters are not treated here.  This should
C             not be a problem in that FORTRAN READs of short records, like
C             assignment of short strings, pads with blanks.
C
C
C     Author:  Robert Kennelly, Informatics General Corporation.
C
C
C     Development history:
C
C        29 Dec. 1984    RAK    Initial design and coding, (very) loosely
C                               based on SCAN_STRING by Ralph Carmichael.
C        25 Feb. 1984    RAK    Added ':' and '=' to list of separators.
C        16 Apr. 1985    RAK    Defined SOLID in terms of variable DUMMY
C                               (previous re-use of STRING was ambiguous).
C         3 Mar. 1986    RAK    Restructured, without SOLID.  Checks for
C                               "state transitions" while executing a
C                               single forward DO loop.  Protect against
C                               funny input (FIRST > LAST).
C        13 May  1988    DAS    Introduced backward search for LAST;
C                               eliminated MIN, MAX and LEN.  (See Notes.)
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   FIRST, LAST, MARK
      CHARACTER
     >   STRING * (*)

C     Local constants.

      CHARACTER
     >   BLANK, EQUAL, COLON, COMMA, HT
      PARAMETER
     >   (BLANK = ' ',
     >    EQUAL = '=',
     >    COLON = ':',
     >    COMMA = ',',
     >    HT = CHAR (9))

C     Local variables.

      INTEGER
     >   HEAD, I, J, TAIL
      LOGICAL
     >   FOUND


C     Execution.
C     ----------

      HEAD = FIRST
      TAIL = LAST

      FIRST = 0
      LAST = 0
      MARK = 0
      FOUND = .FALSE.

C     Look at each character in STRING (HEAD : TAIL) from left to right
C     until a token is found.  (Drop through if LAST > FIRST.)

      DO 30, I = HEAD, TAIL

C        Is the character a separator?

         IF (STRING (I : I) .EQ. BLANK) GO TO 20
         IF (STRING (I : I) .EQ. HT   ) GO TO 20
         IF (STRING (I : I) .EQ. COMMA) GO TO 20
         IF (STRING (I : I) .EQ. EQUAL) GO TO 20
         IF (STRING (I : I) .EQ. COLON) GO TO 20

C           Not a separator.  Check for the beginning of the FIRST token.

            IF (.NOT. FOUND) THEN
               FIRST = I
               FOUND = .TRUE.
            END IF
            GO TO 30

   20    CONTINUE

C           We found a separator.

            IF (FOUND) THEN

C              We just passed the "trailing edge" of the first token.

               MARK = I - 1
               GO TO 40
            END IF
   30 CONTINUE

C     We reached the last character while still seeking the first token.
C     Either this is part of the first token, or MARK and LAST are still zero.

      IF (FOUND) THEN
         MARK = TAIL
         LAST = TAIL
      END IF
      GO TO 99


   40 CONTINUE

C     We found the first token but haven't reached the end yet to adjust LAST.

      DO 60, I = TAIL, MARK + 1, -1

C        Keep searching backwards as long as we have a separator.

         IF (STRING (I : I) .EQ. BLANK) GO TO 60
         IF (STRING (I : I) .EQ. HT   ) GO TO 60
         IF (STRING (I : I) .EQ. COMMA) GO TO 60
         IF (STRING (I : I) .EQ. EQUAL) GO TO 60
         IF (STRING (I : I) .EQ. COLON) GO TO 60

C           Not a separator, so we've found LAST.

            LAST = I
            GO TO 99
   60 CONTINUE

C     Dropped through, so there were no more tokens.

      LAST = MARK
     
C     Termination.
C     ------------

   99 RETURN
      END
