C+----------------------------------------------------------------------
C
      SUBROUTINE SELECT (PROMPT, NLINES, MENU, SUPRES, LUNCRT, LUNKBD,
     >                   SELNUM, SELNAM, QUIT)
C
C     Purpose:
C
C           SELECT modularizes the situation of choosing from a menu by
C        either NUMBER or NAME.  Its innovative (?) feature lies in its
C        parsing of the supplied menu text to identify valid responses,
C        thereby avoiding an additional dictionary as an argument. This
C        makes use of SELECT easy, at the expense of some complexity here.
C
C           If choice by NUMBER ALONE is adequate, it is probably better
C        to display the menu in-line and just prompt with, say, READER.
C        If device-dependent cursor control utilities are available and
C        acceptable from a portability point of view, a much different
C        approach could be taken.  This approach is as system-independent
C        as possible.
C
C           SELECT displays the menu unless the calling program suppresses
C        it, in which case the program user can still request the menu to
C        be shown if necessary.  A default selection must be defined by
C        the input values of SELNUM and SELNAM.  Valid keyboard entries
C        include all reasonable tokens or part-tokens from the menu text,
C        so long as they define a selection uniquely.  An ambiguous or
C        unidentified response leads to a diagnostic and a reprompt.
C
C           A numeric input is accepted right away if it is an in-range
C        item number, although the first tokens of each menu line are still
C        checked for possible duplication by a bad application (as well as
C        to identify the menu line number corresponding to the input).
C        (Descriptive menu text may contain numerical tokens which could
C        clash with the intended item number, so the normal search of all
C        tokens of every line is inappropriate for a numeric input.)
C
C           Two examples should clarify usage.
C
C           This version shows the user prompt and the default selection
C        in inverse video.  It also provides for having no default, via
C        a blank SELNAM on input, meaning both CR and CTRL Z mean quit
C        (handy for some applications).
C
C     Sample PROMPT:
C
C           'Select operating mode.'
C
C        The prompt seen on the screen will then be
C
C           Select operating mode.
C           ? = menu; CTRL Z = quit; <CR> = DISPLAY:
C
C     Corresponding sample menu text:
C
C           '0 = DISPLAY only (plot and/or tabulate - no changes)',
C           '1 = RECTIFY leading edge definition',
C           '2 = NORMALIZE or denormalize coordinates',
C           '3 = REDISTRIBUTE the abscissas',
C            :   :   :   :   :   :   :   :   :
C
C           (The calling program has this text in a DATA statement.)
C
C           Note that an item NUMBER is expected to be the FIRST token in
C        each menu line, and the corresponding NAME is taken as the SECOND
C        token, which should be unique to that item.  Tokens are separated
C        by any of the delimiters in SEPS below.  These include '.' so that
C        a menu could look like the following example, and so that decimal
C        points in the reply are ignored as probably unintentional.
C
C                Input MENU                       Output SELNUM, SELNAM
C
C           '   -1.  Return to previous menu.',              -1  'RETURN'
C           '    0.  Proceed.',                               0  'PROCEED'
C           '    1.  Display graphics parameters.',           1  'DISPLAY'
C                :   :   :   :   :   :   :   :                :   :   :
C
C        (Use of SELNAM in the calling program may not be appropriate with
C        this example because the mnemonics have not been well chosen.)
C        Other possible forms of item number are 1... and (1).
C
C           SYNONYMS may also be selected from the item descriptions as long
C        as they are unique across items.  An example of how this could be
C        handy is in item 2 above: 'denorm' could be a perfectly valid entry,
C        which would be interpreted as a request for operating mode 2.
C
C           Whatever valid entry (name or number) is identified, the returned
C        value of SELNAM will always be the corresponding SECOND token on the
C        line, fully expanded, in upper case, and left-justified. SELNAM should
C        therefore allow for the appropriate number of characters.
C
C     Arguments:
C
C        Name    Dimension  Type  I/O/S  Description
C        PROMPT      -        C     I    Prompt to appear first on the
C                                        screen (to which additional
C                                        standard text will be appended -
C                                        see example above).
C        NLINES      -        I     I    Number of lines of menu text.
C        MENU     NLINES      C     I    Menu, with each line in the form
C                                        shown by the above example - an
C                                        integer for token 1; a key name
C                                        for token 2; and possible synonym
C                                        names among further tokens.  Avoid
C                                        'HELP' for token 2 - it will be
C                                        interpreted as a request for the menu.
C        SUPRES      -        L     I    SUPRES=T means don't display the
C                                        menu (unless the user asks for it).
C        LUNCRT      -        I     I    Logical unit number for the screen.
C        LUNKBD      -        I     I    Logical unit number for the keyboard.
C        SELNUM      -        I    I/O   Input with default item number.
C                                        Output with item number selected from
C                                        the FIRST tokens on the menu lines.
C                                        SELNUM is NOT NECESSARILY within the
C                                        range [1:NLINES]. In the above example,
C                                        SELNUM = 0 is perfectly valid.
C        SELNAM      -        C    I/O   Input with default item name (which
C                                        should match the SECOND token of the
C                                        item referred to by SELNUM on input).
C                                        However a blank input (no default)
C                                        means both CR and CTRL Z will be
C                                        interpreted as meaning "QUIT."
C                                        Output with key name selected from
C                                        the SECOND tokens on the menu lines,
C                                        guaranteed upper case/left-justified.
C                                        May be preferred in calling program
C                                        over SELNUM, for mnemonic reasons.
C                                        Required length of SELNAM should not
C                                        exceed the parameter constant MXNAME
C                                        hard-coded below.
C        QUIT        -        L     O    QUIT = .TRUE. means Control-Z (EOF)
C                                        was entered and other outputs are
C                                        probably undefined.  See also
C                                        "SELNAM = blank" description.
C
C     Implementation notes:
C
C     1.    It would be convenient for the user and the programmer here not
C        to have to return a name as well as a number, because then the key
C        name would not have be token number TWO precisely, whereas if an
C        item number is entered and a corresponding name IS to be returned,
C        it must be in a known place, chosen as the second token.  The only
C        rationale is that the calling program may prefer to be mnemonic in
C        its use of the returned NAME rather than the returned NUMBER.  Here,
C        we provide for that choice, in the spirit of the lower-level LOOKUP
C        utility.
C
C     2.    Errors are handled in two main ways: they are either programmer
C        errors or user input errors.  If the menu text is found to break
C        the rules of the game, the application programmer has erred, and
C        execution halts here with a diagnostic.  If the selection from the
C        menu is too short to be unique, or is out of range or is otherwise
C        invalid, another selection is solicited as part of this routine's
C        purpose of ensuring meaningful inputs.
C
C           A third class of error consists of those that could conceivably
C        happen (such as two menu lines with the same leading integer, or
C        blank menu lines) but are more trouble than it's worth to handle.
C        Hooks are left in these cases, just in case someone has the urge
C        to make the routine more bullet-proof than it really needs to be.
C
C     3.    The SCAN2 version of SCANNR turned out to be the appropriate
C        utility here. TOKENS and LOOKUP could probably have been applied,
C        but the exact match for SELNUM, and the desired delimiters, would
C        cause difficulties.  Also, it was nice to avoid the local storage
C        that TOKENS would need to tokenize one menu line at a time.
C
C     Procedures:
C
C        NUMBER  Logical function for deciding if a string represents a number.
C        SCAN2   Finds first/last characters of next token; variable delimiters.
C        UPCASE  Converts text to upper case.
C
C     Environment:
C
C        VAX/VMS and IRIS 4D; FORTRAN 77 with some extensions:
C           >  IMPLICIT NONE
C           >  (A, $) carriage control
C           >  Escape sequences for inverse video
C
C     History:
C
C        09/24/86  DAS  Initial design and coding, for program PROFILE,
C                       using ideas from LOOKUP, etc., by R.A.Kennelly.
C        09/27/86  DAS  Switched from SCANNR to SCAN2; simplified defaults.
C        10/09/86  DAS  Added SUPRES argument for more flexibility, and
C                       inverse video at slight penalty to portability.
C        02/20/87  DAS  Provided for calling with SELNAM = "blank", as
C                       it is not always appropriate to have a default.
C                       In this case, both CR and CTRL Z mean "quit."
C        06/10/88  DAS  Use LAST = LEN (SELNAM), not MXNAME, since SCAN2
C                       was modified.
C        08/12/89  DAS  SCAN2 call on SELNAM uses SEPS now, not BLANK,
C                       to save scanning a default menu line externally.
C        02/07/90  DAS  Removed concatenation from I/O lists, and went
C                       to (A,$) in place of ('$',A) (both for IRIS 4D).
C                       Also moved escape sequences for inverse video
C                       out of parameter constants, where CHAR (27) and
C                       concatenation probably aren't portable.
C        03/19/97  DAS  The check for ambiguous inputs when a numeric
C                       selection is made can encounter unintended
C                       matches within the body of the menu text.
C                       Therefore, only FIRST tokens on each line should
C                       be checked for validity and ambiguity if the
C                       response is numeric.  (Is that clear?)
C
C     Author:
C
C        David Saunders, Sterling Software/NASA Ames, Mountain View, CA.
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   NLINES, LUNCRT, LUNKBD, SELNUM
      LOGICAL
     >   SUPRES, QUIT
      CHARACTER
     >   PROMPT * (*), MENU (NLINES) * (*), SELNAM * (*)

C     Local constants:

      INTEGER
     >   MXNAME
      CHARACTER
     >   BLANK, DFLTXT * 33, QUITXT * 28, SEPS * 7

      PARAMETER
     >  (BLANK   = ' ',
     >   DFLTXT  = ' ? = menu; CTRL Z = quit; <CR> = ',
     >   QUITXT  = ' ? = menu; CTRL Z or <CR> = ',
     >   MXNAME  = 16,
     >   SEPS    = ' =:.()/')

C     Local variables:

      INTEGER
     >   FIRST, I, I1, I2, LAST, LENMNU, LINE, MARK, N, NCHARS,
     >   NHITS, N1, N2, STATUS
      LOGICAL
     >   INTEGER, LISTIT, NODEFLT
      CHARACTER
     >   ESCAPE * 1, HELP * 4, INVERSE * 4, REPLY * (MXNAME), RESET * 4,
     >   TOKEN * (MXNAME)

C     Procedures:

      LOGICAL
     >   NUMBER
      EXTERNAL
     >   NUMBER, SCAN2, UPCASE

C     Storage:

      DATA
     >   HELP /'HELP'/


C     Execution.
C     ----------

C     Set up for inverse video (formerly done in VAX-dependent PARAMETERs).
C     2 formats and 3 I/O lists will need changing if this is removed.

      ESCAPE  = CHAR (27)
      INVERSE = ESCAPE // '[7m'
      RESET   = ESCAPE // '[0m'

C     Delimit the input default value for the selection name (may be blank).

      I1 = 1
      LAST = LEN (SELNAM)
      CALL SCAN2 (SELNAM, SEPS, I1, LAST, I2)

      NODEFLT = LAST .EQ. 0
      LISTIT = .NOT.SUPRES

 300  CONTINUE

C     Display menu only if calling program or program user asks for it.

      IF (LISTIT) THEN
         WRITE (LUNCRT, 1003) BLANK, MENU
      END IF

  310 CONTINUE

C     Issue supplied prompt plus standard default/quit/help info.

      WRITE (LUNCRT, 1005) INVERSE, PROMPT, RESET
 
      IF (NODEFLT) THEN
         WRITE (LUNCRT, 1006) QUITXT, INVERSE, 'quit', RESET, ': '
      ELSE
         WRITE (LUNCRT, 1006) DFLTXT, INVERSE, SELNAM (I1 : I2), RESET,
     >      ': '
      END IF

C     Accept the response, which should not exceed MXNAME characters.

      REPLY = BLANK
      READ (LUNKBD, 1001, IOSTAT=STATUS) REPLY

      IF (STATUS .GT. 0) THEN
         WRITE (LUNCRT, 1002) 'Input error -try again.'
         GO TO 310
      END IF

C     Check for CTRL Z or carriage return.

      IF (STATUS .LT. 0) THEN
         QUIT = .TRUE.
         GO TO 999
      END IF

      IF (REPLY .EQ. BLANK) THEN
         QUIT = NODEFLT
         GO TO 999
      END IF

      QUIT = .FALSE.

C     Response must be an item number, an item name, or a request to
C     display the menu, or it could be an invalid input.
C     Tokenize the response, and convert it to upper case, ignoring any
C     decimal point in what is presumably meant to be an item number.

      N1 = 1
      LAST = MXNAME
      CALL SCAN2 (REPLY, ' .', N1, LAST, N2)
      CALL UPCASE (REPLY (N1 : N2))
      N = N2 - N1 + 1

      IF (REPLY (N1 : N2) .EQ. '?' .OR.
     >    REPLY (N1 : N2) .EQ. HELP (1 : N)) THEN

C        (Some applications may find 'HELP' here clashes with an item name.
C        Just allowing for '?' here would be the way to go then.)

C        Display the menu, then reprompt.

         LISTIT = .TRUE.
         GO TO 300
      END IF

C     Try to match the reply with some token in the menu.  If the reply
C     is an item NUMBER, the match must be EXACT ('1' must not be claimed
C     to be a short form of '10', for instance).
C
C     The strategy is to isolate each token in the menu text and compare
C     its first N characters with the N significant characters of the
C     response string.  This has to be done in upper case, and since we
C     don't want to change the menu text in-place, we have to transfer
C     each token to a buffer of length MXNAME and UPCASE that.
C
C     Then ambiguities must be dealt with by searching the rest of the
C     menu for a second hit, in which case another response is solicited.

      INTEGER = NUMBER (REPLY (N1 : N2))
      LENMNU = LEN (MENU (1))
      NHITS = 0

C     For each line in the menu ...

      DO 600, I = 1, NLINES
         FIRST = 1
         LAST = LENMNU

C        For each token in the line (or just the FIRST if numeric) ...

  500    CONTINUE
            CALL SCAN2 (MENU (I), SEPS, FIRST, LAST, MARK)

            IF (MARK .GT. 0) THEN
C              Ignore tokens shorter than the reply.

               NCHARS = MARK - FIRST + 1
               IF (NCHARS .GE. N) THEN
                  IF (.NOT. INTEGER .OR.
     >               (INTEGER .AND. NCHARS .EQ. N)) THEN

                     TOKEN (1 : N) = MENU (I) (FIRST : MARK)
                     CALL UPCASE (TOKEN (1 : N))

                     IF (REPLY (N1 : N2) .EQ. TOKEN (1 : N)) THEN
                        NHITS = NHITS + 1
                        IF (NHITS .EQ. 1) THEN
C                          Ignore rest of this line, but go on to next line.
                           LINE = I
                           GO TO 600
                        ELSE
C                          Two hits:
                           WRITE (LUNCRT, 1002)
     >                        'Ambiguous selection - could be:',
     >                        MENU (LINE),
     >                        'or:',
     >                        MENU (I),
     >                        BLANK,
     >                        'Try again.'
                           GO TO 310
                        END IF
                     END IF
                  END IF
               END IF

               IF (.NOT. INTEGER) THEN
                  FIRST = MARK + 2
                  GO TO 500
               ELSE
                  ! Scan only FIRST tokens for numeric responses
               END IF

            END IF
  600 CONTINUE

      IF (NHITS .EQ. 0) THEN
         WRITE (LUNCRT, 1002) 'No such choice - try again.'
         GO TO 310
      END IF

C     Return valid choice in fully expanded/upper case form of tokens 1 and 2.

      FIRST = 1
      LAST = LENMNU
      CALL SCAN2 (MENU (LINE), SEPS, FIRST, LAST, MARK)
      READ (MENU (LINE) (FIRST : MARK), 1004) SELNUM

      FIRST = MARK + 2
      CALL SCAN2 (MENU (LINE), SEPS, FIRST, LAST, MARK)
      CALL UPCASE (MENU (LINE) (FIRST : MARK))
      SELNAM = MENU (LINE) (FIRST : MARK)


C     Termination.
C     ------------

  999 RETURN

C     Formats.
C     --------

 1001 FORMAT (A)
 1002 FORMAT (1X, A)
 1003 FORMAT (4X, A)
 1004 FORMAT (BN, I4)
 1005 FORMAT (/, 1X, 3A)
 1006 FORMAT (5A, $)

      END
