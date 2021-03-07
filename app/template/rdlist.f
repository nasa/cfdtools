C+----------------------------------------------------------------------
C
      SUBROUTINE RDLIST (SCREEN, PROMPT, KEYBRD, NUMBER, LIST)
C
C
C     Description and usage:
C
C        RDLIST prompts for and reads a list of integers.  The list is of
C     indefinite length, and may contain sub-ranges which are expanded here.
C
C        Sample input (contrived to show the possibilities):
C
C     -1:1, 5, 7-10 4 -3  20:18       (read as a character string and parsed)
C
C        Corresponding output:
C
C     -1, 0, 1, 5, 7, 8, 9, 10, 4, -3, 20, 19, 18       (as an integer array)
C
C        Likely applications would be to entering "first, last, & increment,"
C     or specifying which transducers to display results for in an experiment,
C     or (the initial requirement) selecting a sequence of items from a menu.
C
C        The standard delimiters (comma, blank, or tab) are assumed between
C     normal list entries, while either a colon or a hyphen is permitted for
C     specifying a sub-range.  Negative increments are permitted in sub-range
C     specifications at no cost.  Blanks should not be mixed with hyphens in
C     sub-ranges because of ambiguity with negative integers.  Colons should
C     be used similarly.  Thus -5 -3 means the list -5, -3 and not the list
C     -5, -4, -3, while -5: -3 or -5 : -3 or -5 :3 are treated as invalid.
C     The valid forms would be -5:-3 or -5--3 (both meaning -5, -4, -3).
C
C        (The hyphen is too natural in the common positive integer cases for
C     it not to be allowed in addition to the colon, but it forces sub-ranges
C     to have no embedded blanks.  No doubt colons could be treated differently
C     from hyphens, but not in this version yet.)
C
C
C     Arguments:
C
C        Name    Dimension  Type  I/O/S  Description
C        SCREEN              I    I      Logical unit number to which
C                                        the prompt is written.
C        PROMPT              C    I      Prompt string to be written to
C                                        unit SCREEN.
C        KEYBRD              I    I      Logical unit number from which
C                                        data is to be read.
C        NUMBER              I    I/O    Input with the maximum number of
C                                        integers provided for.  On output,
C                                        NUMBER > 0 is the no. of values read;
C                                        NUMBER = 0 corresponds to <CR>;
C                                        NUMBER < 0 corresponds to EOF (^Z).
C        LIST    NUMBER      I      O    Array of integers entered (if any).
C
C
C     Implementation Notes:
C
C        The prompting part could well be done in-line, but READER is likely
C     to be used elsewhere by the application, so using it here too makes sense.
C
C        The SCAN2 utility makes virtually any manipulation of the input list
C     straightforward, so there is no need for arbitrary restrictions.  However,
C     assuming the standard delimiters as indicated above keeps the argument
C     list down.  (So does use of NUMBER <= 0 in place of READER-like <CR> and
C     <EOF> arguments.)
C
C        Disallowing blanks in sub-ranges has the virtue of making both single
C     integers and sub-ranges "tokens" in the usual sense.  In fact, TOKEN2
C     could well have been applied here, but the arbitrary limits on number of
C     tokens and maximum length of each one are eschewed in favor of delimiting
C     tokens in the buffer directly, since SCAN2 would still have to be applied
C     to TOKEN2 output.
C
C        Handling of more than one line's worth of data is not attempted in this
C     version, but a trailing hyphen (say) could easily serve as a flag.
C
C        Error handling is mostly accomplished by trapping bad internal reads
C     when the characters are converted to integers.  The offending substring
C     is displayed, and reentering is invited.  Overflow of the list is also
C     trapped, with options to proceed with a truncated list or to reenter.
C
C
C     Procedures:
C
C        READER  Prompting utility. Its READS option is used here to buffer
C                the input, and READY is handy when bad inputs are detected.
C        SCAN2   Finds positions of the first and last significant characters.
C
C
C     Environment:  DEC VAX/VMS FORTRAN, including:
C
C        >  IMPLICIT NONE
C        >  Trailing ! comments
C        >  Some 7-character variable names
C
C
C     Author:  David Saunders, Sterling Software, Palo Alto, CA.
C
C
C     08/30/88  DAS  Adapted TOKEN2 (Kennelly) as RDLIST for the INTEGER case,
C                    with sub-range handling being its main raison d'etre.
C                    Forcing of contiguous sub-ranges (no embedded blanks)
C                    to allow use of the hyphen is troubling, since the
C                    colon avoids ambiguity.  But the algorithm (dealing
C                    with one token at a time) would be affected seriously
C                    if colons were treated differently from hyphens.
C     02/20/89  DAS  Achieved limited cursor control via variable buffer.
C     07/23/04   "   Raised MAXBUF from 80 to 132 (likewise in READER).
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      CHARACTER
     >   PROMPT * (*)
      INTEGER
     >   KEYBRD, LIST (*), NUMBER, SCREEN

C     Local constants.

      INTEGER
     >   MAXBUF
      CHARACTER
     >   BLANK * 1, HYPHEN * 1

      PARAMETER
     >  (BLANK  = ' ',
     >   HYPHEN = '-',
     >   MAXBUF = 132)

C     Local variables.

      INTEGER
     >   COUNT, FIRST, I1, I2, INC, LAST, M, MARK, STATUS
      LOGICAL
     >   CR, EOF, PROCEED
      CHARACTER
     >   BUFFER * (MAXBUF), DELIMS * 2, FORMI * 6, SEPS * 3

C     Procedures.

      EXTERNAL
     >   READS, READY, SCAN2

C     Execution.


      DELIMS = ':' // HYPHEN            ! For defining sub-ranges
      SEPS = BLANK // ',' // CHAR (9)   ! Token separators, including tab
      FORMI = '(In)'                    ! Assume I9 format (max.) suffices
                                        ! (9 digits +ve; 8 digits -ve)

C     (Re)prompt user for input:

  200 CONTINUE
      BUFFER = BLANK                    ! Avoids referencing partial buffer.
      LAST = MIN (NUMBER * 10, MAXBUF)  ! Shorten buffer if possible to use
      COUNT = 0                         ! READER's limited cursor control.

      CALL READS (SCREEN, PROMPT, KEYBRD, BUFFER (1 : LAST), CR, EOF)

      IF (EOF) GO TO 900
      IF (CR)  GO TO 910

C     Convert character list to integers one token at a time:

      FIRST = 1
      LAST = MAXBUF
  300 CONTINUE

C        Delimit a token, if any:

         CALL SCAN2 (BUFFER, SEPS, FIRST, LAST, MARK)
         IF (MARK .GT. 0) THEN

C           Before checking for a sub-range, bypass possible negative sign:

            IF (BUFFER (FIRST : FIRST) .NE. HYPHEN) THEN
               I1 = 1
            ELSE
               I1 = 2            ! Negative integer encountered
            END IF
            I2 = MARK - FIRST + 1

C           Locate the first (and possibly only) integer:

            CALL SCAN2 (BUFFER (FIRST : MARK), DELIMS, I1, I2, M)

C           Check for trailing : or - characters (or nothing but):

            IF (I2 .NE. MARK - FIRST + 1) GO TO 800

C           Update the variable '(In)' format according to the no. of digits:

            WRITE (FORMI (3 : 3), 1001, IOSTAT = STATUS) M
            IF (STATUS .GT. 0) GO TO 800

C           Translate the (sub)token into integer form:

            COUNT = COUNT + 1
            IF (COUNT .GT. NUMBER) GO TO 810

            READ (BUFFER (FIRST : M + FIRST - 1), FORMI,
     >         IOSTAT = STATUS) LIST (COUNT)
            IF (STATUS .GT. 0) GO TO 800

C           If we're not at the end of the token yet, it's a sub-range:

            IF (M .LT. I2) THEN
               I1 = M + 2

C              Adjust the '(In)' format:

               WRITE (FORMI (3 : 3), 1001, IOSTAT = STATUS) I2 - I1 + 1
               IF (STATUS .GT. 0) GO TO 800

C              Translate from character to integer. Reuse of I2 is convenient:

               READ (BUFFER (I1 + FIRST - 1 : MARK), FORMI,
     >            IOSTAT = STATUS) I2
               IF (STATUS .GT. 0) GO TO 800

C              Expand the rest of the sub-range:

               INC = I2 - LIST (COUNT)
               IF (INC .NE. 0) THEN
                  INC = SIGN (1, INC)
                
                  DO 400, I1 = LIST (COUNT) + INC, I2, INC
                     COUNT = COUNT + 1
                     IF (COUNT .GT. NUMBER) GO TO 810

                     LIST (COUNT) = I1
  400             CONTINUE
               ELSE
C                 Sub-range is degenerate - not fatal - nothing to do:
               END IF

            END IF

C           Look for more in the list:

            FIRST = MARK + 2
            GO TO 300

         END IF

      GO TO 910


C     Error handling.

  800 WRITE (SCREEN, 1000) 'Invalid entry:  ', BUFFER (FIRST : MARK),
     >   'Try again.'
      GO TO 200

  810 COUNT = COUNT - 1
      WRITE (SCREEN, 1002) LIST (COUNT)
      PROCEED = .TRUE.
      CALL READY (SCREEN,
     >   'Proceed with truncated list? (<CR>=Yes; ^Z or No=reenter): ',
     >   KEYBRD, PROCEED, CR, EOF)
      IF (PROCEED .AND. .NOT.EOF) GO TO 910
      GO TO 200


C     Termination:

  900 NUMBER = -1         ! <EOF> = ^Z
      GO TO 999

  910 NUMBER = COUNT      ! <CR> (COUNT = 0) or valid list (COUNT > 0)

  999 RETURN

C     Formats.

 1000 FORMAT (1X, A)
 1001 FORMAT (I1)
 1002 FORMAT ('0List overflow. Last valid entry: ', I9)
      END
