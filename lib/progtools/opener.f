C+----------------------------------------------------------------------
C
      SUBROUTINE OPENER (LUNCRT, PROMPT, LUNKBD, FNAME, LUNFIL, FSTAT)
C
C
C  PURPOSE:
C
C        OPENER modularizes the common situation of prompting for a file
C     name, opening the file, and reprompting if the file is not found.
C     Isolating any system dependencies here (such as VAX/VMS's READONLY
C     extension) enhances the portability of typical applications.
C
C        This version is restricted to sequential files, either formatted
C     or unformatted, old or new, with a few of the other occasionally
C     desirable attributes.
C
C        This version also permits proceeding anyway if a file specified
C     as old is not found.  This option required FSTAT to be used as an
C     OUTPUT as well as an input; all other options use it as input only.
C
C        System-dependent feature:
C
C        If there is trouble opening the file, the program user has the
C     option to execute an operating-system command - probably to look in
C     some directory to see what the filename should be - then try again.
C
C
C  ARGUMENTS:
C
C   NAME     DIM    TYPE I/O/S DESCRIPTION
C  LUNCRT     -      I     I   Logical unit for screen prompts.
C  PROMPT     *      C     I   Prompt (possibly indicating default filename).
C                              May be blank to suppress the prompt (but FNAME
C                              must be non-blank in this case).
C  LUNKBD     -      I     I   Logical unit for keyboard entries.
C  FNAME      *      C    I/O  Name of file to be opened - see METHOD/NOTES.
C                              May be blank to permit termination of an
C                              indefinite loop over file names (carriage
C                              return response to prompt, no open attempted,
C                              and a check in the calling program to see if
C                              FNAME is still blank).  A non-blank FNAME on
C                              input will be treated as the default file name.
C  LUNFIL     -      I     I   Logical unit for file to be opened.
C  FSTAT      *      C   I[/O] INPUT:  String of file status attributes,
C                              separated by commas, colons, or blanks.
C                              Unique abbreviations suffice, in upper or
C                              lower case.  The more common possibilities
C                              provided for appear below.
C                              OUTPUT:  None, with one exception: if the
C                              keyword 'IfPresent' is included in the input 
C                              FSTAT string, then FSTAT = 'MISSING' (upper
C                              case) is returned as output if the file is
C                              not found.  ('Old' is implied by 'IfPresent'
C                              here.)  FSTAT must be a character VARIABLE
C                              in this one case; a CONSTANT is fine in all
C                              other cases.
C                              
C  FSTAT examples:
C
C     'NEW, LIST, 160'      New formatted file with up to 160 characters per
C                           record, and implied (not explicit) carriage control
C
C     'old:binary:write'    Existing unformatted file where READONLY access (the
C                           default) is not enough
C
C     'IfPresent'           Suppresses reprompting if the specified file is
C                           not found.
C
C
C  FSTAT token summary:
C
C     (Meaningful combinations should be concatenated.)
C
C  Attribute        Description          Corresponding OPEN keyword   Default
C
C     'OLD'        File already exists            STATUS              'UNKNOWN'
C     'NEW'        New file desired                "  "                 "   "
C     'SCRATCH'    New file; deleted upon closing  "  "                 "   "
C
C     'BINARY'     Unformatted file                FORM              'FORMATTED'
C     'UNFORMATTED'     "      "                   "  "                 "   "
C
C     'NONE'       No implied carriage ctrl     "  "  "    Default for unfrmttd.
C
C     'WRITE'      Write (and read) access        READONLY           The default
C                  is 'READONLY' if STATUS='OLD', else it is 'WRITE'=read/write.
C
C     'nnn'        Max. record length             RECL         <System default> 
C
C     'IfPresent'  See FSTAT description and examples.
C
C  Further FSTAT Notes:
C
C     Integers nnn for RECL refer to bytes (characters) for formatted files
C     (where 132 is normally the upper limit), or longwords (4-byte units)
C     for unformatted files (where the default is system-dependent).
C
C     READONLY is normally needed for reading files of another owner.
C
C     The defaults are also legitimate (if redundant) input values.
C
C
C  METHOD:
C
C     (1) Default the file attributes.  (The READONLY one will be adjusted
C         later if the file status is new or unknown.)
C
C     (2) Decode the subset of attributes indicated in FSTAT, one token at
C         a time - just itemize the finite number of cases in a way that
C         could be extended.  Any unknown attributes are considered fatal
C         programmer errors - stop with a diagnostic.
C
C     (3) Prompt for the name of the file to be opened (unless the prompt
C         is blank - handy for opening files with fixed names).
C
C     (4) If a carriage return is entered, and FNAME was passed to OPENER
C         as all blanks, return immediately without opening a file.  The
C         calling program can detect this case by checking if FNAME is
C         still blank - handy for performing an indefinite loop over
C         multiple files.
C
C     (5) If a carriage return is entered, and FNAME was NOT all blanks,
C         then FNAME is assumed to represent a default filename (presumably
C         indicated as part of the prompt).
C
C     (6) If "EOF" is entered (^Z under VMS; ^D under Unix), this is assumed
C         to mean "STOP."  Perhaps this should be indicated to the calling
C         program via a returned argument value.  However, the original
C         design had the stop occurring here, and there is no good, upwardly
C         compatible way of changing that.
C
C     (7) If the keyword 'IfPresent' has been input in FSTAT, use INQUIRE to
C         detect existence, and return with FSTAT = 'MISSING' if appropriate.
C
C     (8) Attempt to open the file:
C
C         IF the file cannot be opened THEN
C            Inform the user.
C            Prompt for either the correct file or an operating
C               system command (probably for a directory listing).
C            IF the first character found is a '$' THEN
C               Execute the associated command (system-dependent)
C               Go back to (3)
C            ELSE
C               Go back to (8)
C            END IF
C         END IF
C
C  EXTERNAL REFERENCES:
C
C     LOOKUP        Dictionary lookup: permits abbreviations
C     NUMBER        Identifies RECL parameter
C     READER        Prompting utility
C     SCANNR        Identifies tokens in FSTAT string
C     SYSTEM        IRIS utility analogous to VMS's LIB$SPAWN
C     UPCASE        Upper-casifies FSTAT prior to SCANNR/LOOKUP
C
C  ERROR HANDLING:
C
C        Normally, if a file cannot be opened, a message to that effect is
C     sent to the screen.  The user then has three options:  re-enter the
C     file, which is consequently opened; type a command to list directory
C     contents, whereupon the prompt/response cycle is reinitiated; or quit.
C
C        Alternatively, a missing file may be handled by the calling program.
C     See the FSTAT description for more.
C
C        File attributes found in FSTAT (e.g. 'NULL' for the BLANK keyword)
C     which are not (yet) handled by OPENER are fatal: OPENER will STOP.
C
C  NOTES:
C
C        The declared length of the name of the file to be opened should be
C     enough to accommodate a possibly lengthy directory specification.  A
C     generous size, such as CHARACTER * 60, is suggested.
C
C  SYSTEM DEPENDENCIES:
C
C     (1) IMPLICIT NONE is nonstandard.
C     (2) READONLY is a VAX/VMS extension - removed in this version.
C     (3) CALL SYSTEM is IRIS-specific.
C
C  ENVIRONMENT:  IRIS/IRIX, FORTRAN 77
C
C  HISTORY:
C
C  02/27/86  R.G.Langhi  Initial design and code (formatted files only).
C
C  05/24/86  D.Saunders  Provided for defaulting filename, and for
C                        returning with FNAME=' ' and no file opened
C                        to handle the indefinite-loop-over-files case;
C                        documentation clarified (sequential/formatted).
C
C  08/28/86  RGL         Added unformatted sequential capability.
C
C  05/28/87  RGL         Bug:  the default filename was lost after an
C                        erroneous response to the prompt; need to keep
C                        a copy locally.  On cancellation, FNAME is
C                        returned now with the default instead of blank.
C                        (Note:  The default filename MAY be blank.)
C
C  05/10/88  DAS         Generalized FSTAT usage to indicate more than
C                        just old/new and formatted/unformatted, giving
C                        CARRIAGECONTROL, READONLY, and RECL control too;
C                        suppressed prompt if it is blank, so that OPENER
C                        can still be used with fixed file names.
C
C  02/15/90  DAS         Generalized FSTAT further to permit the calling
C                        program to proceed without the specified (old)
C                        file (as opposed to OPENER's reprompting, which
C                        remains the usual option).
C                        Also: changed from READC to READS to accommodate
C                        the case-sensitive Unix community.
C
C  02/22/90  DAS         IRIS 4D version: suppressed LIB$SPAWN feature
C                        (is there an IRIX equivalent?); ^D instead of ^Z
C                        indicated in the prompts.
C
C  03/14/90  DLH/DAS     Dexter found "call system" as a substitute for
C                        LIB$SPAWN.
C
C  08/28/90  DAS         Unformatted files needed CC='NONE'.
C
C  06/10/91  DAS         NUMBER was declared INTEGER - should be LOGICAL.
C
C  10/15/99  DAS         Removed CARRIAGECONTROL & READONLY for SGI/IRIX f90.
C
C  AUTHOR:  Ronald Langhi, NASA Ames/Sterling Software, Palo Alto, CA
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   LUNCRT, LUNKBD, LUNFIL
      CHARACTER
     >   FNAME * (*), FSTAT * (*), PROMPT * (*)

C     Local constants:

      INTEGER
     >   MXCHAR, MXWORD
      CHARACTER
     >   BLANK * 1

      PARAMETER
     >  (BLANK  = ' ',
     >   MXCHAR = 11,
     >   MXWORD = 13)

C     Local variables:

      INTEGER
     >   ENTRY, FIRST, LAST, MARK, RECLEN
      CHARACTER
     >   ATTRIB * 80, CC * 7, DICT (MXWORD) * (MXCHAR), FORIG * 80,
     >   FORM * 11, KEY * (MXCHAR), STATUS * 7
      LOGICAL
     >   CR, DFRECL, ENQUIRE, EOF, PRESENT, RDONLY

C     Procedures:

      LOGICAL
     >   NUMBER
      EXTERNAL
     >   LOOKUP, NUMBER, SCANNR, READS, SYSTEM, UPCASE

C     Storage:

      DATA
     >   DICT
     >      /'BINARY', 'FORMATTED', 'FORTRAN', 'LIST', 'NEW', 'NONE',
     >       'OLD', 'READONLY', 'SCRATCH', 'UNFORMATTED', 'UNKNOWN',
     >       'WRITE', 'IFPRESENT'/
C     The dictionary should be in upper case.  It need not be alphabetized.


C     Execution:

C     Save the default input filename in case of an erroneous response 
C     to the prompt:

      FORIG = FNAME

C     Default the file attributes so that corresponding lookup hits can
C     be ignored (and input string FSTAT can be short):

      FORM    = DICT (2)
      CC      = DICT (3)
      STATUS  = DICT (11)
      RDONLY  = .TRUE.
      DFRECL  = .TRUE.
      ENQUIRE = .FALSE.

C     Ensure that the attributes text is upper case:

      ATTRIB = FSTAT
      CALL UPCASE (ATTRIB)

C     Start of loop over tokens in attributes text:

      FIRST = 1
      LAST = LEN (FSTAT)
   50 CONTINUE

         CALL SCANNR (ATTRIB, FIRST, LAST, MARK)
         IF (MARK .EQ. 0) GO TO 70

         KEY = ATTRIB (FIRST : MARK)
         CALL LOOKUP (MXWORD, DICT, .FALSE., KEY, ENTRY)

         IF (ENTRY .GT. 0) THEN

C           There is only a modest number of possibilities.
C           Avoid replicating the dictionary text by working with
C           subscripts rather than text.  Adding keywords at the
C           end of the dictionary will not affect this code, which
C           is why LOOKUP's non-alphabetic option is used above.

            IF (ENTRY .EQ. 1 .OR. ENTRY .EQ. 10) THEN
               FORM = DICT (10)
               CC = DICT (6)
            END IF
            IF (ENTRY .EQ. 4 .OR.
     >          ENTRY .EQ. 6)    CC = DICT (ENTRY)
            IF (ENTRY .EQ. 12)   RDONLY = .FALSE.
            IF (ENTRY .EQ. 5 .OR.
     >          ENTRY .EQ. 7 .OR.
     >          ENTRY .EQ. 9)    STATUS = DICT (ENTRY)
            IF (ENTRY .EQ. 13)   ENQUIRE = .TRUE.

         ELSE IF (NUMBER (ATTRIB (FIRST : MARK))) THEN
            DFRECL = .FALSE.
            READ (ATTRIB (FIRST : MARK), '(BN, I11)') RECLEN
         ELSE
            GO TO 810
         END IF

         FIRST = MARK + 2
         IF (FIRST .LE. LAST)
     >GO TO 50

   70 CONTINUE

C     Check for inconsistencies:

      IF (STATUS .NE. DICT (7)) RDONLY = .FALSE.


  100 CONTINUE

C     Start of loop over retries at opening the file:

      IF (PROMPT .EQ. BLANK) GO TO 300

C     The prompt should indicate any default file name here.
C     Use READS instead of the original READC, for Unix reasons.

      CALL READS (LUNCRT, PROMPT, LUNKBD, FNAME, CR, EOF)

  200 CONTINUE

      IF (CR .AND. FNAME .EQ. BLANK)  GO TO 999
      IF (EOF) GO TO 800

  300 CONTINUE

C     Try to open the file, but test for its existence first in one case:

      IF (ENQUIRE) THEN
         INQUIRE (FILE=FNAME, EXIST=PRESENT)
         IF (.NOT. PRESENT) THEN
            FSTAT = 'MISSING'
            GO TO 999
         END IF
      END IF

C     Oddball READONLY keyword, and uncertain system default for RECL
C     force four possibilities here:

      IF (RDONLY .AND. DFRECL) THEN

         OPEN (UNIT=LUNFIL, FILE=FNAME, FORM=FORM, STATUS=STATUS,
     >      ERR=400)

      ELSE IF (RDONLY .AND. .NOT.DFRECL) THEN

         OPEN (UNIT=LUNFIL, FILE=FNAME, FORM=FORM, STATUS=STATUS,
     >      RECL=RECLEN, ERR=400)

      ELSE IF (.NOT.RDONLY .AND. DFRECL) THEN

         OPEN (UNIT=LUNFIL, FILE=FNAME, FORM=FORM, STATUS=STATUS,
     >      ERR=400)

      ELSE
C****    IF (.NOT.RDONLY .AND. .NOT.DFRECL) THEN

         OPEN (UNIT=LUNFIL, FILE=FNAME, FORM=FORM, STATUS=STATUS,
     >      RECL=RECLEN, ERR=400)

      END IF

      GO TO 999


C     OPEN error handling:

  400 CONTINUE

      WRITE (LUNCRT, 1000) 'Error opening following file:', FNAME
      FNAME = FORIG

      IF (FNAME .EQ. BLANK) THEN
         CALL READS (LUNCRT,
     >      'Try again, look in directory ("$ls ..."), cancel open ' //
     >      '(CR), or stop (^D): ', LUNKBD, FNAME, CR, EOF)
      ELSE
         CALL READS (LUNCRT,
     >      'Try again, look in directory ("$ls ..."), open default ' //
     >      '(CR), or stop (^D): ', LUNKBD, FNAME, CR, EOF)
      END IF

C     Either another attempt at the file name was entered ...

      IF (FNAME (1 : 1) .NE. '$') GO TO 200

C     ... or an operating-system command was requested (system-dependent) ...

      CALL SYSTEM (FNAME (2 :))

C     ... and start over:

      FNAME = FORIG
      GO TO 100


  800 WRITE (LUNCRT, 1000) 'Stopping as requested.'
      GO TO 990

  810 WRITE (LUNCRT, 1000) 'OPENER: Bad file attribute in FSTAT:',
     >   FSTAT (FIRST : MARK), 'Aborting.'
C**** GO TO 990

  990 STOP ' '

  999 RETURN

C     Formats:

 1000 FORMAT (1X, A)

      END
