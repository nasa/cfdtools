C+------------------------------------------------------------------------------
C
      SUBROUTINE PLOTDEV (MODE, NAME, LUNCRT, LUNKBD, LUNDEV, EXPLAIN,
     >   NLIST, LIST, TYPELIST, PARMLIST, DEVINDEX)
C
C
C PURPOSE:
C
C        PLOTDEV prompts for and/or initializes the current graphics output
C     device, for applications using the CA-DISSPLA graphics package.  In
C     the case of metafile output, the basic file name is application-dependent
C     and hence indicated as an argument.  The appropriate file extension is
C     appended here.
C
C        This version is more elaborate than the original SMDLIB version in
C     order to permit switching devices repeatedly, as needed for flexible
C     previewing and hard-copying of multiple plot frames in a single run.
C
C        This version also provides a means of determining the graphics
C     output devices available to a given application at a given site.
C     Doing this via a configuration file in a choice of standard locations
C     is more easily extensible than hard-coding the lists of screen
C     devices and metafiles in the application.  Device attributes such
C     as X, Y scale factors needed for true inches are also handled here.
C
C
C USAGE:
C
C        Suppression of devices and selection of default devices (the first
C     in the preview/metafile lists) is determined by the configuration file
C     associated with each application.  The format of this file is indicated
C     under the DEVINDEX argument description below.
C
C        PLOTDEV's various modes are of two basic kinds: one kind sets up
C     the available list(s) of devices; the second kind picks one device
C     from a list and/or (re)initializes the device.
C
C        There are also two kinds of ways an application might do the
C     switching between previewing and hard-copy: the simple approach
C     provides for previewing one frame at a time until the user either
C     quits or requests hard-copy, in which case any input data files are
C     rewound and all frames are replotted to the metafile; the more
C     elaborate approach avoids rewinding the data by providing for
C     plotting of each frame to either one or two output devices, with
C     a "hard-copy-all-the-rest" option.  Thus an application may have
C     as few as two calls to PLOTDEV, or as many as five.
C
C        When initializing a device, the first use of PLOTDEV raises the
C     DISSPLA level from 0 to 1.  (An exception is for CGM metafile
C     initialization, which leaves the level at 0.)  In this case if
C     PLOTDEV is called at DISSPLA level 1, the level will remain at 1.
C
C
C ARGUMENTS:
C
C     ARG    TYPE/DIM   I/O/S  DESCRIPTION
C
C     MODE     C*1      I      Controls the mode of operation as follows:
C
C                |             'P': Read list of PREVIEW devices supported
C                |                  for this application at this site, from
C         Read one device           the file NAME.config in a standard place.
C       list.  NAME applies    'M': Read list of METAFILES supported from the
C       to the config. file.        configuration file.
C                |             'A': Reads list of ALL active devices from the
C                |                  configuration file (as needed by simpler
C                |                  applications).
C
C                |             'S': SELECT output device (preview or metafile)
C                |                  via a prompt; do not initialize the device.
C    Select and/or init. one   'I': INITIALIZE device according to the input
C    device. NAME applies to        device code index, without prompting.
C      the metafile if any.    'B': BOTH prompt for and initialize a device
C                |                  from the given LIST.  DEVINDEX is returned
C                |                  (see its description).
C                              
C     NAME     C*(*)    I      Name used to establish the input configuration
C                              file name and the output metafile name(s).
C                              The application program's name may serve in
C                              both cases, or the output file name may differ
C                              if NAME is used as indicated alongside the
C                              description of MODE.  Example: if NAME = 'qplot'
C                              and MODE = 'P', 'M', or 'A', the configuration
C                              file must be qplot.config in one of the standard
C                              locations shown below.
C                              For the other class of MODEs, NAME = 'qplot'
C                              would define the output PostScript file as
C                              qplot.ps.  But NAME here might instead be derived
C                              from the name of the application program's input
C                              data file, or from a user prompt.
C                               
C     LUNCRT     I      I      Unit number for screen output
C     LUNKBD     I      I      Unit number for keyboard input
C     LUNDEV     I      I      Unit number for device output (or config. file)
C
C     EXPLAIN  C*(*)    I      Explanatory text prefacing menu defined by LIST
C                              (applicable to MODE = 'S' or 'B').
C
C     NLIST      I      I/O    No. of device codes in LIST.  For MODE = 'P',
C                              'M' or 'A', NLIST should be input with the
C                              size of the array passed as LIST(*); this will
C                              be updated with the number of entries in the
C                              returned list.  For MODE = 'S' or 'I' or 'B',
C                              input the NLIST found from a prior call.
C
C     LIST   I(NLIST)   I/O    Subset of device codes defining desired menu.
C                              The first element should be the default code.
C                              Output for MODE = 'P', 'M', or 'A'; input for
C                              'S', 'I', or 'B'.  See DEVINDEX for options.
C
C   TYPELIST C(NLIST)*3   O    3-letter mnemonic device code corresponding to
C                              LIST(*) from config. file.  Output as for LIST.
C
C   PARMLIST R(NLIST,*)   O    Parameters associated with the output device.
C                              An arbitrary number per device is handled
C                              by scanning the config. file until a non-numeric
C                              token is encountered.  So far:
C                              PARMLIST (1) = X scale factor for true inches;
C                              PARMLIST (2) = Y   "     "     "     "     "
C
C   DEVINDEX     I      I/O    Device code index into the LIST arrays
C                              (input if MODE = 'I'; output if MODE = 'S' or
C                              'B'; ignored otherwise).  See definitions below.
C
C                              DEVINDEX = 0 output means no device was selected.
C
C     The configuration file format is as follows, with undesired options
C     commented out via '!'.  This should match the DATA statement below
C     (too much trouble to avoid the latter).
C
C     Standard locations searched    Samples
C
C     graph$:[NAME]                  graph$:[qplot]qplot.config
C     sys$system:                    sys$system:qplot.config
C     /usr/local/src/graphics/NAME   /usr/local/src/graphics/qplot/qplot.config
C
C     ! qplot.config file on node RALph
C     !
C     ! Type  Dev.code  Xscale  Yscale  ! Device/metafile
C     ! -------------------------------------------------
C     TEK       1       1.0     1.0     ! Tektronix 4014, default preview device
C     TEK       2       1.0     1.0     ! Tektronix 4105/4107
C     TEK       3       1.0     1.0     ! Tektronix 4109
C     ! TEK     4       1.0     1.0     ! Tektronix 4115
C     VT        5       1.0     1.0     ! VT240
C     ! LN      6       1.0     1.0     ! LN03 Plus
C     PS        7    0.996933  1.00719  ! PostScript, def. metaf.; Apple scales
C     CGM       8       1.0     1.0     ! CGM
C     ! DIS     9       1.0     1.0     ! DISSPOP
C     DIP      10       1.052   1.043   ! DIP    (QMS 1200 scale factors)
C     ! X      11       1.0     1.0     ! X Window System
C     ! <to be extended>
C
C EXTERNAL REFERENCES:
C     GETLINE    Gets one line of text, skipping over '!' lines.
C     NUMBER     Identifies (non)numeric tokens.
C     READER     Prompting utility.
C     SCANNR     Used here to trim trailing blanks.
C     UPCASE     Uppercase utility.
C     Also:      Numerous CA-DISSPLA device interfaces.
C
C ENVIRONMENT:   VAX/VMS, FORTRAN 77;  also IRIX
C                CA-DISSPLA, Version 11.0-9003
C
C HISTORY:
C
C     05/15/89   PLOTDEV (formerly SETDEV) was adapted from various ideas found
C     DAS/MDW    in SMDLIB and demo programs to solve the following problems:
C                1: IDEV had to be returned so that the screen could be cleared
C                   for SOME devices.
C                2: Several metafiles had to be handled, and these should be
C                   named after the application program.
C                3: A mechanism was needed for suppressing devices known to
C                   SMDLIB but irrelevant to the environment without losing
C                   information.  A hard-coded mapping was used as a reasonably
C                   maintainable compromise.
C
C     10/27/89   Adapted SMDLIB version to initialize devices supported by the
C     MDW        DISSPLA graphics package, including metafiles such as DIP, CGM,
C                PostScript and DISSPOP plus the option of choosing additional
C                Tektronix devices through further prompting.
C                Introduced argument MODE to separate prompting and initializ-
C                ation, and arguments LIST, NLIST to suppress devices from the
C                menu at the application level.
C
C     01/04/90   Deleted mythical GKS option.  Stripped trailing blanks from
C     MDW        device prompt.
C
C     01/23/90   Enabled hardware characters for PostScript and CGM output.
C     MDW/DAS    EXPLAIN argument found desirable to specify terminals and
C                metafiles in separate calls to PLOTDEV.
C
C     11/20/90   A correction.  Now passing 60 instead of 100 through the
C     MDW        argument list in CGM metafile calls.
C
C     01/23/91   Added configuration file scheme to keep changes to the
C     DAS        available devices from affecting applications.  Here was the
C                logical place, at some cost in explaining the different modes.
C                The argument list was reordered, and the two TYPE arguments
C                were added in the hope of keeping hard-coding of device codes
C                in application programs to a minimum.
C
C     02/19/91   Added PARMLIST argument to handle the X, Y scale factor
C     DAS        problem; changed IDEV to DEVINDEX (the relevant index into what
C                are now THREE lists); eliminated corresponding TYPE argument.
C
C     03/04/91   PostScript's pen-width was not being initialized when /PS was
C     DAS/DBS    specified on QPLOT's command line (to suppress prompting).
C                Would passing it IN via PARMLIST (3) (say) make sense?
C                For now, a data statement keeps the application simpler, but
C                forces it to invoke the menu to override the default.
C
C     03/22/91   The prompt for type of CGM file was in the wrong place, and
C     DAS        the <NAME>.cgm file was not being set up properly.
C
C     04/30/91   Had to introduce UPCASE when tried on the IRIS, which also
C     DAS        didn't like concatenation in a passed CHARACTER argument
C                (.pop and .dip file names).
C                Lots of trouble with the output PostScript file name on the
C                IRIS.  (Frames after the first went to fort.3, where LUNDEV
C                was 3.)  Had to avoid any kind of initialization if the
C                device did not change from the previous call (MODE = 'I').
C
C     05/31/91   DIP and DISSPOP metafile names were using self-counting 100s
C     DAS        where NCHARS should have been.
C
C     06/12/91   Eliminated DIP for Cray version (%REF (STRING) won't compile).
C                "Standard" location for *.config file had to be changed.
C                READONLY keyword suppressed in the OPEN for *.config.
C                IOMGR's weak way of dealing with file names is word-length
C                dependent.  LENWORD=4 or 8 and A4 or A8 are still hard-coded.
C
C     01/04/92   Metafile name(s) may not be related to the application program
C     DAS        name - made this clear in describing the two uses of NAME.
C
C AUTHORS:       David Saunders/Michael Wong
C                NASA Ames/Sterling Software, Palo Alto, CA.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:
C     ----------

      INTEGER
     >   DEVINDEX, LUNCRT, LUNKBD, LUNDEV, NLIST, LIST (*)

      REAL
     >   PARMLIST (NLIST, *)

      CHARACTER
     >   EXPLAIN * (*), NAME * (*), MODE * 1, TYPELIST (*) * 3

C     Local constants:
C     ----------------

      INTEGER
     >   LENSTR, LENWORD, MAXDEV, MAXTXT
      CHARACTER
     >   BLANK * 1
      PARAMETER
     >  (BLANK  = ' ',
     >   LENSTR = 60,   ! 60 >= MAXTXT + 27  (See prompt for menu item.)
                        ! NOTE: STRING is also used for the PostScript file
                        ! name set up, and thus should be at least 4*15 chars.
     >   LENWORD= 4,    ! Wordlength in bytes (8 for CRAY, 4 for DEC, SGI)
     >   MAXDEV = 10,
     >   MAXTXT = 30)   ! Length of each menu item

C     Local variables:
C     ----------------

      INTEGER
     >   FIRST, I, IDEV, IDEVLAST, IOS, ITEM, LAST, MARK, NCHARS,
     >   NWORDS, PSCALLS, IBUFF (16)
      REAL
     >   PWIDTH, XPAGE, YPAGE
      LOGICAL
     >   CR, EOF, PRE
      CHARACTER
     >   ALLDEV (MAXDEV) * (MAXTXT), CGMTYPE * 1, IFMT * 4, RFMT * 10,
     >   STRING * (LENSTR), TYP * 3

C     Procedures:
C     -----------

      LOGICAL
     >   NUMBER
      EXTERNAL
     >   NUMBER

C     Storage:
C     --------

      DATA
     +   ALLDEV /
     1   'Tektronix (monochrome) (4014) ',
     2   'Tektronix (color) (4105, 4107)',
     3   'Tek 4109                      ',
     4   'Tek 4115                      ',
     5   'VT240                         ',
     6   'LN03 Plus                     ',
     7   'PostScript metafile           ',
     8   'CGM metafile                  ',
     9   'DISSPOP metafile              ',
     +   'DIP metafile                  '/

      DATA
     >   CGMTYPE  /'B'/,
     >   IDEVLAST /-1/,     ! Initialize it to an invalid device code
     >   PSCALLS  /0/,
     >   PWIDTH   /.007/,
     >   IFMT     /'(In)'/,
     >   RFMT     /'(BN,Fnn.0)'/

      SAVE
     >   ALLDEV, CGMTYPE, IDEVLAST, IFMT, PSCALLS, PWIDTH, RFMT

C     Execution:
C     ----------

C     Either process the configuration file ...
C     -----------------------------------------

      IF (MODE .EQ. 'P' .OR. MODE .EQ. 'M' .OR. MODE .EQ. 'A') THEN

C        Look for the file NAME.config in one of several standard locations:

         DO 100, I = 1, 4
            IF (I .EQ. 1) STRING = '/codes/graphics/' // NAME // '/'
            IF (I .EQ. 2) STRING = '/codes/aero/'     // NAME // '/'
            IF (I. EQ. 3) STRING = 'graph$:[' // NAME // ']'
            IF (I .EQ. 4) STRING = 'sys$system:'
            LAST = INDEX (STRING, BLANK)
            STRING (LAST :) = NAME // '.config'
            LAST = LAST + LEN (NAME) + 6
            OPEN (UNIT = LUNDEV, FILE = STRING (1 : LAST),
C****>            STATUS = 'OLD', IOSTAT = IOS, READONLY)
     >            STATUS = 'OLD', IOSTAT = IOS)
            IF (IOS .EQ. 0) GO TO 120          ! Success
  100    CONTINUE

         WRITE (LUNCRT, '(/, A, A, A)')
     >     ' PLOTDEV: Unable to open ', NAME, '.config file. Aborting.',
     >     BLANK
         GO TO 990

  120    CONTINUE

C        Read the config. file till EOF.  There are so few choices for
C        preview (screen) devices that they are itemized here rather
C        requiring another column (P or M) to distinguish them from metafiles.

         ITEM = 0
  140    CONTINUE

            CALL GETLINE (LUNDEV, '!', STRING, LAST, IOS)

            IF (IOS .GT. 0) THEN
               WRITE (LUNCRT, '(A)')
     >            ' PLOTDEV: Error reading config. file.  Aborting.',
     >            BLANK
               GO TO 990
            END IF

            IF (IOS  .LT. 0) GO TO 160         ! EOF
            IF (LAST .EQ. 0) GO TO 140         ! Empty or suppressed line

C           Isolate the first token:
            FIRST = 1
            CALL SCANNR (STRING, FIRST, LAST, MARK)
            TYP = STRING (FIRST : MARK)
            CALL UPCASE (TYP)
            PRE = TYP .EQ. 'TEK' .OR.
     >            TYP .EQ. 'VT'  .OR.
     >            TYP .EQ. 'X'

            IF (MODE .EQ. 'A' .OR.
     >         (MODE .EQ. 'P' .AND. PRE) .OR.
     >         (MODE .EQ. 'M' .AND. .NOT. PRE)) THEN    ! Add it to the list.
               ITEM = ITEM + 1
               TYPELIST (ITEM) = TYP

C              The second token is the integer device code:

               FIRST = MARK + 2
               CALL SCANNR (STRING, FIRST, LAST, MARK)
               WRITE (IFMT (3 : 3), '(I1)') MARK - FIRST + 1
               READ (STRING (FIRST : MARK), IFMT) LIST (ITEM)

C              Remaining numeric tokens are device parameters.
C              E.g.: X, Y scale factors for true inches.

               I = 0
  150          CONTINUE

                  FIRST = MARK + 2
                  CALL SCANNR (STRING, FIRST, LAST, MARK)
                  IF (NUMBER (STRING (FIRST : MARK))) THEN
                     I = I + 1
                     WRITE (RFMT (6 : 7), '(I2)') MARK - FIRST + 1
                     READ (STRING (FIRST : MARK), RFMT)
     >                  PARMLIST (ITEM, I)
                     GO TO 150    ! No need to check for exceeding LAST.
                  END IF

            END IF
            IF (ITEM .LT. NLIST)
     >   GO TO 140

  160    NLIST = ITEM
         CLOSE (UNIT = LUNDEV)
      END IF                  ! List is set up - return.


C     ... or select and/or initialize a device:
C     -----------------------------------------

  400 IF (MODE .EQ. 'S' .OR. MODE .EQ. 'B') THEN

C        Prompt for and assign code for device or metafile.

         WRITE (LUNCRT, 1001) BLANK, EXPLAIN, BLANK
         WRITE (LUNCRT, 1002)
     >      (ITEM, ALLDEV (LIST (ITEM)), ITEM = 1, NLIST),
     >      NLIST + 1, '(or EOF) None of the above'
         WRITE (LUNCRT, 1001)

C        Strip trailing blanks from description of default device (LIST(1)):

         FIRST = 1
         LAST  = LEN (ALLDEV (LIST (1)))
         CALL SCANNR (ALLDEV (LIST (1)), FIRST, LAST, MARK)

         STRING (1:24) = 'Item number? (Def = 1 = '
         STRING (25:24+LAST) = ALLDEV (LIST (1)) (1:LAST)
         STRING (25+LAST:27+LAST) = '): '
         ITEM = 1
         CALL READI (LUNCRT, STRING (1:27+LAST), LUNKBD, ITEM, CR, EOF)

         IF (EOF .OR. ITEM .EQ. NLIST + 1) THEN  ! Quit.
            DEVINDEX = 0
            GO TO 999
         END IF
         IF ( ITEM .LE. 0 .OR. ITEM .GT. NLIST) THEN ! Must have been a mistake.
            GO TO 400
         END IF
         DEVINDEX = ITEM

         IF (TYPELIST (ITEM) .EQ. 'PS') THEN
            PWIDTH = 7.
            CALL READR (LUNCRT,
     >         'Pen width in thousandths of an inch? (<CR> = 7): ',
     >         LUNKBD, PWIDTH, CR, EOF)
            PWIDTH = PWIDTH * .001

         ELSE IF (TYPELIST (ITEM) .EQ. 'CGM') THEN
            CALL READC (LUNCRT, 'CGM metafile type? ' //
     >         '(B(inary), C(har.), T(ext); <CR>=B): ',
     >         LUNKBD, CGMTYPE, CR, EOF)
         END IF

      END IF

      IF (MODE .EQ. 'I' .OR. MODE .EQ. 'B') THEN

         IDEV = LIST (DEVINDEX)
         IF (IDEV .EQ. IDEVLAST) GO TO 999  ! No need to do anything.
                                            ! Too much indenting to avoid GO TO.
         IDEVLAST = IDEV

C        Initialize the selected device or metafile.
C        First, TURN OFF hardware characters for most devices, and
C        set the device configuration to primary I/O (normal default):

         IF (IDEV .NE. 7 .AND. IDEV .NE. 8) THEN
            CALL RESET ('HWCHAR')
            CALL IOMGR (0, -102)
         END IF

         IF (IDEV .EQ. 1) THEN                         ! Tektronix 4014 mono

            CALL TK4014 (960, 0)

         ELSE IF (IDEV .EQ. 2) THEN                    ! Tektronix 4105/7 color

            CALL PTK41

         ELSE IF (IDEV .EQ. 3) THEN                    ! Tektronix 4109

            CALL TK41DO (1, 4109)

         ELSE IF (IDEV .EQ. 4) THEN                    ! Tektronix 4115

            CALL TK41DO (1, 4115)

         ELSE IF (IDEV .EQ. 5) THEN                    ! VT240

            CALL VT240

         ELSE IF (IDEV .EQ. 6) THEN                    ! LN03 Plus

            CALL LN01TK (LUNDEV)

         ELSE IF (IDEV .EQ. 7) THEN                    ! PostScript

            PSCALLS = PSCALLS + 1
            IBUFF (1) = 16            ! Buffer length
            CALL IOMGR (IBUFF, -1)    ! Initialize I/O system

C           Set the I/O configuration:

            CALL IOMGR (5, -102)      ! Direct data file output

C           Set up new user-supplied file name in IBUFF (1 : 15) (awkward):

            STRING = NAME
            NCHARS = MIN (LEN (NAME) + 3, 15 * LENWORD)
            STRING (NCHARS - 2 : NCHARS) = '.ps'
            NWORDS = (NCHARS + LENWORD - 1) / LENWORD

	    DO 500, I = 1, 16
               IBUFF (I) = 0
  500       CONTINUE

            DO 510, I = 1, NWORDS
               READ (STRING ((I - 1) * LENWORD + 1 : I * LENWORD),
     >            '(A4)') IBUFF (I)
  510       CONTINUE

            IBUFF (16) = NCHARS
            CALL IOMGR (IBUFF, -103)

C           Set the file mode:

            IF (PSCALLS .EQ. 1) THEN
               IBUFF (1) = 3    ! No overwrite (= new version where applicable)
            ELSE
               IBUFF (1) = 0    ! Append
            END IF

            CALL IOMGR (IBUFF, -104)

C           Set the logical unit number:

            CALL IOMGR (LUNDEV, -110)

C           Set the file carriage control and padding:

            IBUFF (1) = 0             ! Carriage control off
            IBUFF (2) = 0             ! Shouldn't be needed (cc type = space)
            IBUFF (3) = 0             ! Pad with spaces, not nulls
            CALL IOMGR (IBUFF, -111)

C           Plot area:                ! Switch between portrait and landscape
                                      ! at the application level?
            XPAGE = 7.99              ! Recommended for Apple LaserWriter.
            YPAGE = 10.78

            CALL PSCRPT (XPAGE, YPAGE, PWIDTH)

            IF (PSCALLS .EQ. 1) THEN  ! Reset to append, else we get a new
               IBUFF (1) = 0          ! PostScript file for each plot on VAX.
               CALL IOMGR (IBUFF, -104)
            END IF

C           Specify hardware characters for all available alphabets,
C           with reasonable tolerances:

            CALL HWCHAR (20., 20., 90., 1, 'HEBREW')

         ELSE IF (IDEV .EQ. 8) THEN                    ! CGM metafile

            STRING = NAME
            NCHARS = LEN (NAME) + 4
            STRING (NCHARS - 3 : NCHARS) = '.cgm'

            IF (CGMTYPE .EQ. 'B') CALL CGMBO (STRING, NCHARS, 0)
            IF (CGMTYPE .EQ. 'C') CALL CGMCO (STRING, NCHARS, 0)
            IF (CGMTYPE .EQ. 'T') CALL CGMTO (STRING, NCHARS, 0)

            CALL HWCHAR (20., 20., 90., 1, 'HEBREW')

         ELSE IF (IDEV .EQ. 9) THEN                    ! DISSPOP metafile

            STRING = NAME
            NCHARS = LEN (NAME) + 4
            STRING (NCHARS - 3 : NCHARS) = '.pop'
            CALL COMPRS
            CALL SETCPR (LUNDEV, 0, 0, 0)
            CALL POPNAM (STRING, NCHARS)

         ELSE IF (IDEV .EQ. 10) THEN                   ! DIP metafile

            STRING = NAME
            NCHARS = LEN (NAME) + 4
            STRING (NCHARS - 3 : NCHARS) = '.dip'
C*****      CALL DIP (LUNDEV, %REF (STRING), NCHARS)
            WRITE (LUNCRT, 1001) 'PLOTDEV: DIP is no longer supported.'
            GO TO 990

         END IF
      END IF

      GO TO 999

  990 STOP ' '

  999 RETURN

C     Formats:
C     --------

 1001 FORMAT (1X, A)
 1002 FORMAT (3X, I2, ': ', A)

      END
