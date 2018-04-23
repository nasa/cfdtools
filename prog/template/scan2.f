C+----------------------------------------------------------------------
C
      SUBROUTINE SCAN2 (STRING, SEPS, FIRST, LAST, MARK)
C
C
C     Description and usage:
C
C           Looks for non-blank fields ("tokens") in a string, where the
C        fields are of arbitrary length and separated by any of a set of
C        user-specified separators (e.g., blanks or commas).  The position
C        of the end of the first token is also returned so that this
C        routine may be conveniently used within a loop to process an
C        entire line of text.
C
C           The procedure examines a substring, STRING (FIRST : LAST), which
C        may of course be the entire string (in which case just call SCAN2
C        with FIRST = 1 and LAST = LEN (STRING) ).  The indices returned
C        are relative to STRING itself, not the substring.
C
C
C     Arguments:
C
C        Name    Dimension  Type  I/O/S  Description
C        STRING    *         C    I      Text string containing data to be
C                                        scanned.
C        SEPS      *         C    I      String containing the separators.
C                                        Each character in SEPS counts as a
C                                        token delimiter.
C        FIRST               I    I/O    Index of first character of interest
C                                        in STRING.  FIRST >= 1 is assumed.
C                                        Output is index of the beginning of
C                                        the first non-separator, or 0 if no
C                                        token was found.
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
C        (1)  IMPLICIT NONE is non-standard.
C
C        (2)  FIRST >= 1 and LAST <= LEN (STRING) are assumed upon entry,
C             because these are the obvious bounds needed by the calling
C             program on the first scan of a given string - no need to
C             evaluate LEN (STRING) here with every call, while use of
C             MIN and MAX tends to disguise programming errors.
C
C        (3)  Entering with FIRST > LAST is an eventual normal consequence
C             of setting FIRST = MARK + 2 for further scans of the same string.
C             The zero DO loop iteration count feature of the language is
C             taken advantage of here - the result is simply to drop
C             out the bottom with FIRST = LAST = MARK = 0.
C
C        (4)  This version adjusts LAST using a BACKward search, which is
C             degenerate for all tokens after the first.  (Doing it in
C             the one forward loop means unnecessary repeated tokenizing
C             to find the end.)
C
C        (5)  Zeroing all of FIRST, LAST, and MARK (rather than just MARK)
C             when no token is found may, in retrospect, be undesirable in
C             that information is lost (esp. LAST), but there are too many
C             applications to change this now.
C
C
C     Author:  Robert Kennelly, Informatics General Corporation.
C
C
C     Development history:
C
C         4 Mar 1986    RAK    Variation of SCANNR, which is hard-coded
C                              (for historical reasons and a modest speed
C                              advantage) for separators BLANK, TAB,
C                              COMMA, COLON, and EQUAL.
C
C         5 May  1988   DAS    Reverse search used to find LAST; MAX, MIN,
C                              and LEN also eliminated (see Notes).
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   FIRST, LAST, MARK
      CHARACTER
     >   SEPS * (*), STRING * (*)

C     Local variables.

      INTEGER
     >   HEAD, I, J, NSEPS, TAIL
      LOGICAL
     >   FOUND


C     Execution.
C     ----------

      NSEPS = LEN (SEPS)
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

         DO 10, J = 1, NSEPS
            IF (STRING (I : I) .EQ. SEPS (J : J)) GO TO 20
   10    CONTINUE

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

C        Keep searching as long as we have a separator.

         DO 50, J = 1, NSEPS
            IF (STRING (I : I) .EQ. SEPS (J : J)) GO TO 60
   50    CONTINUE

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
