C+----------------------------------------------------------------------
C
      SUBROUTINE RDREALS (SCREEN, PROMPT, KEYBRD, NUMBER, REALS)
C
C
C     Description and usage:
C
C        RDREALS reads an indefinite number of real values from one line or
C     record, with an optional prompt.  Standard delimiters (comma, blank,
C     or tab) are assumed between values.
C
C        One motivation here is to overcome an inherent limitation in the
C     earlier READER utility (which can accept just ONE value per prompt)
C     while still offering READER's error handling (reprompting) and its
C     flexibility with respect to null inputs (<CR> only).  [Standard list-
C     directed reads cannot handle indefinite lists.]  RDREALS also offers
C     more flexible cursor control.
C
C        RDREALS is appropriate for entering multiple values into a single
C     array.  Example: prompting for an arbitrary number of phase angles in
C     an oscillation.  Higher level routine RDTUPLES should be used for
C     entering (possibly many) pairs or triples.
C
C        RDLIST is available for the [in]definite-list-of-integers case.
C     There is no integer analog of RDTUPLES.
C
C
C     Arguments:
C
C        Name    Dimension  Type  I/O/S  Description
C        SCREEN              I    I      Logical unit number to which the
C                                        prompt and diagnostics are written.
C        PROMPT              C    I      Prompt string with carriage control
C                                        in PROMPT (1:1).  (Use ' ', '0', '1',
C                                        or '$' in the usual way.)
C                                        If PROMPT is blank it is suppressed.
C        KEYBRD              I    I      Logical unit number from which data
C                                        is to be read.  (May be a disk file if
C                                        called by RDTUPLES in indirect mode.)
C        NUMBER              I    I/O    Input with the maximum number of
C                                        values provided for.  On output,
C                                        NUMBER > 0 is the no. of values read;
C                                        NUMBER = 0 corresponds to null (CR);
C                                        NUMBER =-1 corresponds to end-of-data;
C                                        NUMBER =-2 corresponds to bad data if
C                                                   PROMPT is blank (else a
C                                                   reprompt is issued here);
C                                        NUMBER =-3 means a system-dependent
C                                                   input error was encountered.
C        REALS   NUMBER      R      O    Array of reals entered (if any).
C
C
C     Implementation Notes:
C
C        The prompting part is done in-line because READER's only suitable
C     entry point, READS, does not offer the desired cursor control.
C
C        Handling of more than one line's worth of data (132 characters) is not
C     attempted in this version, but a trailing hyphen (say) could serve as a
C     flag.
C
C        RDTUPLES makes use of RDREALS, and this application required three
C     minor features here that are not expected to be made use of directly:
C
C        >  the prompt may be suppressed;
C        >  indirect input may be REPORTed (but not DONE except upon further
C           calls from RDTUPLES).  An internal COMMON was reluctantly used to
C           return the @filename string (preferable to another argument);
C        >  bad data leads to termination in indirect mode (whereas reprompting
C           is appropriate interactively).  A blank PROMPT is taken to mean
C           indirect and hence an error return if an invalid value is read.
C
C
C     Procedures:
C
C        GETLINE  Reads one record, with option to strip off trailing comments.
C        SCAN2    Finds positions of the first and last significant characters.
C
C
C     Environment:  DEC VAX/VMS FORTRAN, including:
C
C        >  IMPLICIT NONE
C        >  Trailing ! comments
C        >  Some 8-character variable names
C
C
C     Author:  David Saunders, NASA Ames/Sterling Software, Palo Alto, CA.
C
C     02/20/89  DAS  Initial implementation, derived from RDLIST for the
C                    case of real values, with improved cursor control.
C     03/01/89  DAS  (With Robert Kennelly:)  Refinements indicated above
C                    to facilitate the RDTUPLES application.
C     02/21/90  DAS  IRIS 4D requires more cumbersome carriage control:
C                    print 'prompt' from PROMPT = '$prompt' with FORMAT
C                    (1X, A, $), else just use FORMAT (A) for entire PROMPT.
C     10/20/99  DAS  Fortran 90 version (avoiding $ carriage control).
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      CHARACTER
     >   PROMPT * (*)
      INTEGER
     >   KEYBRD, NUMBER, SCREEN
      REAL
     >   REALS (*)

C     Internal COMMON (required by the RDTUPLES application of RDREALS). (Ugh!)

      INTEGER, PARAMETER ::
     >   MAXBUF = 132       ! Arbitrary limit on input record length
      CHARACTER
     >   BUFFER * (MAXBUF)
      COMMON /TUPLES/
     >   BUFFER

C     Local constants:

      CHARACTER, PARAMETER ::
     >   BLANK    * 1 = ' ',
     >   COMMENT  * 1 = '!',
     >   FORMR    * 7 = '(F25.0)',  ! Arbitrary limit on numerical token length
     >   INDIRECT * 1 = '@'

C     Local variables:

      INTEGER
     >   COUNT, FIRST, LAST, MARK, STATUS
      LOGICAL
     >   INTERACT
      CHARACTER
     >   SEPS * 3

C     Procedures:

      EXTERNAL
     >   GETLINE, SCAN2

C     Execution:

      SEPS = ' ,' // CHAR (9)   ! Token separators are blank, comma, or tab
      INTERACT = PROMPT /= BLANK

  200 CONTINUE
      COUNT = 0
      IF (INTERACT) THEN
         IF (PROMPT (1 : 1) == '$') THEN
            WRITE (SCREEN, 1020, ADVANCE = 'NO') PROMPT (2 : )
         ELSE
            WRITE (SCREEN, 1020) PROMPT (2 : )
         END IF
      END IF

      CALL GETLINE (KEYBRD, COMMENT, BUFFER, LAST, STATUS)

      IF (STATUS > 0) THEN                     ! Read error
         IF (INTERACT) GO TO 200
         GO TO 800
      END IF
      IF (STATUS < 0) GO TO 900                ! ^Z = EOF = "quit"
      IF (LAST  == 0) GO TO 910                ! CR (or insignificant line)

C     Convert character list to reals one token at a time:

      FIRST = 1
  300 CONTINUE

C        Delimit a token, if any:

         CALL SCAN2 (BUFFER (1 : LAST), SEPS, FIRST, LAST, MARK)

         IF (MARK > 0) THEN

C           The following is relevant to RDTUPLES only.

            IF (COUNT == 0) THEN                              ! Indirection?
               IF (BUFFER (FIRST : FIRST) .EQ. INDIRECT) THEN ! E.g. @XY.DAT
                  BUFFER (FIRST : FIRST) = BLANK         ! Lets RDTUPLES look
                  BUFFER (1 : 1) = INDIRECT              ! in a fixed place.
                  NUMBER = 1           ! Helps termination test in RDTUPLES.
                  GO TO 999
               END IF
            END IF

            COUNT = COUNT + 1

C           Convert the token to a real value (internal read).

            READ (BUFFER (FIRST : MARK), FORMR, IOSTAT = STATUS)
     >         REALS (COUNT)
            IF (STATUS > 0) GO TO 810

            FIRST = MARK + 2
            IF (COUNT .LT. NUMBER) GO TO 300    ! Look for more if there's room,
                                                ! else ignore further entries.
         END IF

      GO TO 910

C     Error handling:

  800 WRITE (SCREEN, 1000) ' *** System error.  IOSTAT: ', STATUS
      IF (INTERACT) GO TO 200
      NUMBER = -3
      GO TO 999

  810 WRITE (SCREEN, 1010) ' *** Bad value found: ',
     >   BUFFER (FIRST : MARK)
      IF (INTERACT) GO TO 200
      NUMBER = -2
      GO TO 999

C     Termination:

  900 NUMBER = -1       ! EOF or ^Z with COUNT = 0
      GO TO 999

  910 NUMBER = COUNT    ! CR, insignificant line, or '@...' with COUNT = 0,
                        ! else valid list with COUNT > 0

  999 RETURN

C     Formats:

 1000 FORMAT (A, I6)
 1010 FORMAT (A, A)
 1020 FORMAT (1X, A)

      END SUBROUTINE RDREALS
