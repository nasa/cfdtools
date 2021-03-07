C+------------------------------------------------------------------------------
C
      SUBROUTINE READER (CASE, SCREEN, PROMPT, KEYBRD, CR, EOF,
     >                   INT4, REAL4, REAL8, CHARS, STRING, YES)
C
C
C     Description and usage:
C
C           READER provides a set of utilities offering some of the
C        flexibility of list-directed input to interactive FORTRAN programs,
C        while permitting use of <CR> (i.e. Carriage Return) as a legal,
C        identifiable response.  This is impossible with a list-directed
C        read, which ignores <CR> and continues to wait for input.  In
C        case of bad input (presumably due to user error), a warning is
C        displayed and the prompt repeated.  "Quit" (as opposed to "default")
C        is provided for by the "End-of-file" argument EOF.  (See Notes.)
C        If CR or EOF is .TRUE., the argument value returned is unchanged,
C        so these flags should be checked by the calling routine.  (However,
C        assigning a default value before calling READER can save checking
C        the CR argument.)
C
C           This version avoids the original's ENTRY points for the different
C        data types because of f90 compiler difficulties on an SGI system.
C
C           READER itself should not be called directly. Rather, six ancillary
C        routines provide for reading INTEGER, REAL, DOUBLE PRECISION, STRING,
C        CHARACTER, or YES/NO data, and each makes the appropriate call to
C        READER.  The names are READx, where x = I, R, D, S, C, or Y.  These
C        ancillary routines appear at the end of the READER source module.
C
C           READC returns all the non-blank characters found, converted to
C        upper case if alphabetic, and packed from the left into the output
C        string with blanks removed.  The READS option merely returns the
C        string as entered from the keyboard, without modification.  READC
C        is normally appropriate for entering a single item or token (such as
C        an identifier or a file name), while READS is appropriate for literal
C        strings such as plot titles.  HOWEVER: since Unix file names are case
C        sensitive, READS is actually the better choice for file name prompts
C        on Unix systems.
C
C           All calling sequences are identical except for the type of the
C        value to be returned.
C
C           Sample usage:
C
C           :      :                              (^D under Unix)
C       210 NPTS = 100                              |
C           CALL READI (LUNCRT,                     |
C          >   'Enter number of points.  <CR>=100; ^Z=quit: ',
C          >   LUNKBD, NPTS, CR, EOF)
C           IF (EOF) GO TO 999
C           IF (NPTS .LT. 1 .OR. NPTS .GT. MXPTS) GO TO 210
C           :      :
C
C
C     Arguments:
C
C        Name    Dimension  Type  I/O/S  Description
C        CASE                C    I      'I', 'R', 'D', 'S', 'C', or 'Y'
C                                        (but see use of READI, etc., above)
C        SCREEN              I    I      Logical unit number to which
C                                        the prompt is written (often 6).
C        PROMPT    *         C    I      Prompt string to be written to
C                                        unit SCREEN.
C        KEYBRD              I    I      Logical unit number from which
C                                        data is to be read (often 5).
C
C
C   ---> Only ONE of the following six output choices is to be used:
C
C        INT4                I      O    INTEGER quantity to be returned.
C        REAL4               R      O    REAL quantity to be returned.
C        REAL8               D      O    DOUBLE PRECISION quantity to be
C                                        returned.
C        STRING              C      O    A CHARACTER string, as entered.
C        CHARS               C      O    A CHARACTER string, converted to
C                                        upper case if alphabetic, packed
C                                        and left-justified.
C        YES                 L      O    Logical flag set .TRUE. for a 'YES'
C                                        response (first character 'Y' or 'y')
C                                        and .FALSE. for 'N' or 'n'.  Other
C                                        inputs force a reprompt.
C
C        CR                  L      O    CR = .TRUE. if a null value (i.e.
C                                        Carriage Return only) was read,
C                                        else CR = .FALSE.
C        EOF                 L      O    EOF = .TRUE. if "End-of-File"
C                                        was detected (^Z under VMS; ^D
C                                        under Unix), else EOF = .FALSE.
C                                        
C
C     External devices:  See arguments SCREEN and KEYBRD.
C
C
C     External procedures:  UPCASE is used by the READC and READY options
C                           in place of the original in-line code - see Notes.
C
C
C     Environment:  DEC/OpenVMS Fortran 90
C                   SGI/IRIX    Fortran 90
C
C     Method:
C
C           The user is prompted for input, with the cursor left at the
C        end of the prompt if there is room for the expected response.  The
C        reply is read into a buffer as a character string and, except for
C        the READS and READC entries, left-justified with blanks squeezed
C        out and re-read internally in the requested format.  For READC, the
C        same repacking is done, then the result is converted to upper case.
C        No repacking or conversion is done in the case of READS.
C
C
C     Notes:
C
C        (1)  Conversions to upper case were originally done in-line by
C             adding the appropriate offset (= ICHAR ('A') - ICHAR ('a')).
C             However this has been shown to fail on some systems and has
C             been replaced by use of the UPCASE utility, which cannot fail.
C
C        (2)  The length of the string returned by READS and READC is set by
C             the calling program - excess characters read will be truncated
C             from the right.  No error flag is set for this condition.  The
C             maximum length for CHARS and STRING is MAXBUF.
C
C        (3)  The '$' carriage control character used for short prompts is
C             not compatible with Fortran 90 - this version uses ADVANCE='NO'
C             instead.
C
C
C     Author:  Robert Kennelly, NASA Ames/Sterling Software, Palo Alto, CA.
C
C
C     Development history:
C
C        27 Sep. 1983    RAK    Initial design and coding.
C         6 Oct. 1983    RAK    Extended to multiple data types.
C         2 Nov. 1983  RAK/LJC  Added READS for strings without modification.
C        16 Mar. 1984  RAK/LJC  Included "BN" specifier in formats to avoid
C                               unneccessary packing.
C        14 June 1984  DAS/RAK  Corrected header, TEST always defined, added
C                               CC, eliminated redundant trap on CASE, used
C                               shorter flags, e.g. 'R' instead of 'REAL'.
C         7 Sep. 1984    RAK    Leave STRING unchanged when CR entered,
C                               thus consistent with the other modes.
C        19 Oct. 1984    RAK    Initialize BUFFER prior to READ.
C        14 Dec. 1984    RAK    Make sure length of CHARS is not exceeded
C                               during repacking/case-conversion.
C        19 Dec. 1985  RAK/DBS  Addition of READY entry, and general 
C                               streamlining of code.
C        30 Dec. 1985  RAK/DAS  Edited header.
C        09 Sep. 1988    DAS    Sample usage added above; other cosmetics.
C        05 May  1989  DAS/RGL  Unix- and VMS-compatible now: ^Z/^D usage
C                               documented; STOP 'READER' message <=8 chars.;
C                               OFFSET not defined as a PARAMETER constant.
C                               (UNICOS displays '$' in column 1 at time of
C                               writing, but this may go away...)
C        01 Feb. 1990    DAS    '$' carriage control changed to (A,$) form
C                               (see Notes above; thanks to Scott Thomas).
C                               READS recommended over READC for file names
C                               on Unix systems.  Conversions to upper case
C                               now done via UPCASE - see Notes above.
C        23 July 1997     "     Replaced (A,$) with ADVANCE='NO' for Fortran 90.
C        24 July 1997     "     Eliminated ENTRY points as a work-around for an
C                               SGI compiler.  (See 6 ancillary routines below.)
C        19 May  1999     "     The CASE argument had not been described.
C        20 May  1999     "     Dexter Hermstad solved a control-D (EOF) problem
C                               encountered with IRIX 6 f90, making this version
C                               SGI-specific.
C        23 July 2004     "     Raised all the 80s to 132 for RDLIST purposes.
C        13 Feb  2006     "     IMPLICIT NONE was missing from the ancillary
C                               routines introduced to avoid ENTRY points.
C                               The missing declarations may overcome problems
C                               encountered with program CAVITY_MAP when it is
C                               compiled with the pgf90 compiler.
C-------------------------------------------------------------------------------


C     DECLARATIONS.
C     -------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   INT4, KEYBRD, SCREEN
      REAL
     >   REAL4
      DOUBLE PRECISION
     >   REAL8
      CHARACTER
     >   CASE * 1, CHARS * (*), PROMPT * (*), STRING * (*)
      LOGICAL
     >   YES, CR, EOF

C     Local constants.  Note that the field lengths in the format strings
C     are hardcoded to MAXBUF.

      INTEGER
     >   MAXBUF, MAXLIN
      CHARACTER
     >   BLANK, FORMI * 10, FORMR * 12
      PARAMETER
     >  (BLANK  = ' ',
     >   FORMR  = '(BN, F132.0)',
     >   FORMI  = '(BN, I132)',
     >   MAXBUF = 132,
     >   MAXLIN = 132)

C     Local variables.

      LOGICAL
     >   ADVANCE
      INTEGER
     >   FILL, SCAN, STATUS, TEST
      CHARACTER
     >   BUFFER * (MAXBUF)

C     Procedures.

      EXTERNAL
     >   UPCASE


C     EXECUTION.
C     ----------

C     Set a typical maximum length for response to the non-character
C     entries.  This value and the length of the prompt string are used
C     to determine where the cursor is to be positioned after prompting.

      TEST = 16

      IF (CASE .EQ. 'S') THEN
         TEST = LEN (STRING)
      ELSE IF (CASE .EQ. 'C') THEN
         TEST = LEN (CHARS)
      ELSE IF (CASE .EQ. 'Y') THEN
         TEST = 5
      END IF

      ADVANCE = TEST + LEN (PROMPT) .GT. MAXLIN

C     Prompt user for input.
C     ----------------------

   30 CONTINUE
         BUFFER = BLANK

C        Issue prompt, then read up to MAXBUF characters from the keyboard.

         IF (ADVANCE) THEN
            WRITE (SCREEN, '(1X, A)') PROMPT
         ELSE
            WRITE (SCREEN, '(1X, A)', ADVANCE='NO') PROMPT
         END IF

         READ (KEYBRD, '(A)', IOSTAT = STATUS) BUFFER

C        Re-try on errors, exit on End-of-File, or continue.

         IF (STATUS .GT. 0) THEN
            WRITE (SCREEN, 1000)
            GO TO 30

         ELSE IF (STATUS .LT. 0) THEN ! SGI f90 won't allow further ^Ds
            EOF = .TRUE.              ! unless the keyboard is reopened
            CR = .FALSE.
            OPEN (KEYBRD, FILE='/dev/tty', STATUS='OLD') ! Dexter's fix
            GO TO 990

         ELSE
            EOF = .FALSE.
            CR = (BUFFER .EQ. BLANK)
         END IF


C        Process the data in the buffer, if any.
C        ---------------------------------------

         IF (.NOT. CR) THEN

            IF (CASE .EQ. 'I') THEN
               READ (BUFFER, FORMI, IOSTAT = STATUS) INT4

            ELSE IF (CASE .EQ. 'R') THEN
               READ (BUFFER, FORMR, IOSTAT = STATUS) REAL4

            ELSE IF (CASE .EQ. 'D') THEN
               READ (BUFFER, FORMR, IOSTAT = STATUS) REAL8

            ELSE IF (CASE .EQ. 'S') THEN

C              For literal strings, copy the buffered input directly to
C              the output variable.  If CR is true, we leave STRING alone
C              so that any default value set in the calling routine will
C              remain intact.  No error checking is required.

               STRING = BUFFER

            ELSE 

C              CASE is either 'C' or 'Y' - scan the input buffer, starting
C              from the left, to pack the data and count the non-blanks.

               FILL = 0
               DO SCAN = 1, MAXBUF
                  IF (BUFFER (SCAN : SCAN) .NE. BLANK) THEN
                     FILL = FILL + 1
                     BUFFER (FILL : FILL) = BUFFER (SCAN : SCAN)
                  END IF
               END DO

C              Convert to upper case.

               CALL UPCASE (BUFFER (1 : FILL))

               IF (CASE .EQ. 'C') THEN

                  CHARS = BUFFER (1 : FILL)

               ELSE

C                 CASE must be 'Y' - all we need is the first character.

                  IF (BUFFER (1 : 1) .EQ. 'Y') THEN
                     YES = .TRUE.
                  ELSE IF (BUFFER (1 : 1) .EQ. 'N') THEN
                     YES = .FALSE.
                  ELSE

C                    Buffer value is invalid - generate an ersatz input error
C                    to force a re-prompt.

                     STATUS = 999
                  END IF
               END IF
            END IF 

C           Re-try on errors, exit on End-of-File, else continue.
C           (Dropping through for 'C', 'S', 'Y' cases is OK and saves code.)

            EOF = (STATUS .LT. 0)
            IF (STATUS .GT. 0) THEN
               WRITE (SCREEN, 1000)
               GO TO 30

            END IF

      END IF


C     TERMINATION.
C     ------------

  990 RETURN

C     FORMATS.
C     --------

 1000 FORMAT (' Input error!  Please try again.')

      END SUBROUTINE READER


C     Ancillary subroutines used to avoid entry points in READER.
C     -----------------------------------------------------------
C
C     An SGI compiler also objects to omitted arguments (.., , ,),
C     so dummy arguments are passed as a work-around.

      SUBROUTINE READI (SCREEN, PROMPT, KEYBRD, INT4,   CR, EOF)

      IMPLICIT   NONE

C     Arguments:

      INTEGER    INT4, KEYBRD, SCREEN
      CHARACTER  PROMPT * (*)
      LOGICAL    CR, EOF

C     Local variables:

      REAL       R4, R8
      LOGICAL    Y

      CALL READER ('I', SCREEN, PROMPT, KEYBRD, CR, EOF,
     >             INT4, R4, R8, 'C', 'S', Y)
      END

      SUBROUTINE READR (SCREEN, PROMPT, KEYBRD, REAL4,  CR, EOF)

      IMPLICIT   NONE

C     Arguments:

      INTEGER    KEYBRD, SCREEN
      REAL       REAL4
      CHARACTER  PROMPT * (*)
      LOGICAL    CR, EOF

C     Local variables:

      INTEGER    I4
      REAL       R8
      LOGICAL    Y

      CALL READER ('R', SCREEN, PROMPT, KEYBRD, CR, EOF,
     >             I4, REAL4, R8, 'C', 'S', Y)
      END

      SUBROUTINE READD (SCREEN, PROMPT, KEYBRD, REAL8,  CR, EOF)

      IMPLICIT   NONE

C     Arguments:

      INTEGER    KEYBRD, SCREEN
      DOUBLE PRECISION   REAL8
      CHARACTER  PROMPT * (*)
      LOGICAL    CR, EOF

C     Local variables:

      INTEGER    I4
      REAL       R4
      LOGICAL    Y

      CALL READER ('D', SCREEN, PROMPT, KEYBRD, CR, EOF,
     >             I4, R4, REAL8, 'C', 'S', Y)
      END

      SUBROUTINE READC (SCREEN, PROMPT, KEYBRD, CHARS,  CR, EOF)

      IMPLICIT   NONE

C     Arguments:

      INTEGER    KEYBRD, SCREEN
      CHARACTER  CHARS * (*), PROMPT * (*)
      LOGICAL    CR, EOF

C     Local variables:

      INTEGER    I4
      REAL       R4, R8
      LOGICAL    Y

      CALL READER ('C', SCREEN, PROMPT, KEYBRD, CR, EOF,
     >             I4, R4, R8, CHARS, 'S', Y)
      END

      SUBROUTINE READS (SCREEN, PROMPT, KEYBRD, STRING, CR, EOF)

      IMPLICIT   NONE

C     Arguments:

      INTEGER    KEYBRD, SCREEN
      CHARACTER  PROMPT * (*), STRING * (*)
      LOGICAL    CR, EOF

C     Local variables:

      INTEGER    I4
      REAL       R4, R8
      LOGICAL    Y

      CALL READER ('S', SCREEN, PROMPT, KEYBRD, CR, EOF,
     >             I4, R4, R8, 'C', STRING, Y)
      END

      SUBROUTINE READY (SCREEN, PROMPT, KEYBRD, YES,    CR, EOF)

      IMPLICIT   NONE

C     Arguments:

      INTEGER    KEYBRD, SCREEN
      CHARACTER  PROMPT * (*)
      LOGICAL    YES, CR, EOF

C     Logical variables:

      INTEGER    I4
      REAL       R4, R8

      CALL READER ('Y', SCREEN, PROMPT, KEYBRD, CR, EOF,
     >             I4, R4, R8, 'C', 'S', YES)
      END
