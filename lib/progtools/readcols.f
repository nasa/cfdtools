C+----------------------------------------------------------------------
C
      SUBROUTINE READCOLS (LREAD, LWRITE, XCOL, YCOL, MAXPT, POINTS,
     >   MAXCUR, CURVE, X, Y, BUFFERS, LAST, FINIS, ERROR)
C
C     One-liner:   READ two COLumnS of numerical data
C
C     Description:
C
C           READCOLS reads one set of X, Y pairs until either a curve-ending
C        keyword or EOF is encountered.  To read several sets of data, this
C        routine must be called from within a loop which should terminate
C        normally when flag FINIS is .TRUE. (signalling end of file for the
C        current data file).  The other return conditions are many - see the
C        ERROR argument description and "Notes" below.
C
C           READCOLS was originally written for use by QPLOT, but it may
C        also be useful in other applications involving multi-column files.
C
C           This version returns whenever the first token on a line is
C        found to be non-numeric (presumably a keyword).  All handling of
C        keywords is left to the calling program.  (END CURVE and END
C        FRAME were originally handled in READCOLS, but this was eliminated
C        with a view to implementing a keyword control scheme as opposed to
C        the original namelist scheme of QPLOT.)  This version also avoids
C        the backspacing that namelist handling had demanded.
C
C
C     Arguments:
C
C        Name    Dimension  Type  I/O/S  Description
C        LREAD               I    I      Logical unit number for input.
C        LWRITE              I    I      Logical unit number for output.
C                                        Negative suppresses warnings but
C                                        fatal errors are written to unit
C                                        ABS (LWRITE) unless LWRITE = 0
C                                        in which case NO output is written.
C        XCOL                I    I      Input column number from which X
C                                        data is to be read.  If XCOL=0,
C                                        then X(I) = I for I=1,2,3,....
C        YCOL                I    I/O    Input column number from which Y
C                                        data is to be read.  If YCOL=0,
C                                        then Y(I) = I for I=1,2,3,....
C        MAXPT               I    I      Maximum total number of points to
C                                        be read.
C        POINTS  MAXCUR      I    I/O    Number of points in each curve.
C                                        This array indexes the X and Y
C                                        data arrays.
C        MAXCUR              I    I      Maximum number of curves.
C        CURVE               I    I/O    Current number of curves.  Enter 0 for
C                                        first call.
C        X       MAXPT       R    I      Packed array of abscissas.
C        Y       MAXPT       R    I      Packed array of ordinates.
C        BUFFERS   2       C*(*)  I/O    IF LAST > 0, BUFFERS (1) (1:LAST)
C                                        contains the last data record read
C                                        (both on input and on output).
C                                        On output BUFFERS (2) stores formatted
C                                        point and curve information, which may
C                                        be used for upper level error handling.
C        LAST                I    I/O    See BUFFERS.  LAST = 0 means the
C                                        buffer is not meaningful.
C        FINIS               L      O    End-of-input flag (legal end of file).
C        ERROR               I      O    Error flag which indicates an input 
C                                        error was encountered. It assumes a
C                                        value as follows:
C                                        ERROR = 0  Normal (including FINIS=T).
C                                               -1  Non-correctable error
C                                                   reported (some software
C                                                   limit has been exceeded).
C                                               +1  Keyword caused termination 
C                                                   of data read.
C                                               +2  Invalid numeric data read.
C
C     Error Handling:
C
C           A message is printed on unit LWRITE for any error condition
C        (if LWRITE > 0), and on unit ABS (LWRITE) for fatal errors.
C        This permits the calling program to suppress warning messages.
C        Processing halts for READ errors judged fatal.  See LWRITE and
C        ERROR descriptions above, and NOTES below.
C
C
C     External References:
C
C        ALPHA   Identifies non-numeric characters.
C        EXIT    Stop and return status flag to system (VAX/VMS FORTRAN).
C        IOCHEK  Returns LOGICAL flags describing READ status.
C        GETLINE Reads one line, suppressing trailing blanks/comments
C        TOKENS  Finds tokens in text string and returns them in an array.
C
C
C     External Files:
C
C        Unit    I/O/S  Description
C        LREAD   I      Input file 
C        LWRITE    O    Output file ("line printer").
C
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN (FORTRAN 77).
C
C
C     Notes:
C
C        (1)  IMPLICIT NONE is non-standard.
C
C        (2)  Some symbols are up to eight characters long.
C
C        (3)  CALL EXIT (N) is VAX-specific.  It causes program termination,
C             closes all files, and returns control to the operating system.
C             Its argument is then available to the system.  N = 3 is defined
C             here to mean abnormal termination due to fatal read error.
C
C        (4)  Since this routine is intended to be embedded within a loop
C             over sets of data (several frames each consisting of several
C             curves, in plotting terminology), there are two sets of logical
C             control flags floating around.  It is important to keep them
C             straight:
C
C                FINIS  Indicates a normal EOF for the current data file.
C                ERROR  Tells the calling program to use what it has (if
C                       CURVE > 0) and then quit if ERROR .NE. 0 or +1.
C                DONE   Indicates when the internal loop over points is
C                       complete. (Local)
C                FATAL  General bail-out for unrecoverable read errors.
C                       (Local)
C
C        (5)  Some of the details concerning the data packing are handled
C             at this level in an attempt to simplify the calling routine
C             (cf. BASE, MAXCUR).
C
C
C     History:
C
C        31 Dec. 1982    RAK    Initial design and coding.
C        21 Feb. 1983    RAK    "End-of-curve" mark may begin in any column.
C        13 June 1983    RAK    Added "end-of-frame" input option.  FINIS
C                               output replaces EOF flag.
C        21 July 1983    RAK    Updated call to IOCHEK (added CONVER).
C        23 Jan. 1984    RAK    Extensively revised, but old data sets are
C                               upward compatible except that "END CURVE"
C                               and "END FRAME" marks are less free.
C                               Multiple-column input format and the packed
C                               data structure for the X and Y arrays are
C                               the significant changes.  A bug in the
C                               handling of the MAXCURth curve was repaired.
C        18 June 1984    RAK    Reduced BACKSPACE-ing by using TOKENS and BN
C                               format (much faster now on large files).
C                               List-directed READs eliminated.  Maximum
C                               number of columns and data field width
C                               changed.  Some GO TOs used to reduce level
C                               of nesting.  Dropped warning on blank lines.
C                               Reordered arguments.
C         5 Feb. 1985    RAK    Moved initialization of POINT to beginning
C                               of routine (exit processing was sometimes
C                               wrong when array bounds were violated).
C        18 May  1987  RAK/DAS  GETLINE introduced for more efficient reads
C                               (VAX-dependent, but easily modified if
C                               required).  Revise FORMATs to use /' ' for
C                               pushing a blank line (not 0).  Print
C                               warnings only if LWRITE > 0 and fatal
C                               error messages only if LWRITE <> 0.
C        31 July 1987    RAK    Increased MAXCOL to 30.  Pass EXCLAM to
C                               GETLINE as comment character.
C         3 Aug. 1987    RAK    Increased MAXTOK from 24 to 40 (enough
C                               for Cray DOUBLE PRECISION, plus a little).
C         2 Jan. 1990  MDW/RAK  Interpretation of unexpected alpha or numeric 
C                               input passed to to a higher level.  Ensured the
C                               reading of a single column of Y-data even when
C                               YCOL set incorrectly. 
C         5 Mar. '90  M.D.Wong  Keyword dictionary now passed as argument to
C                               terminate curve read.
C        16 Nov. '91 D.Saunders Major revision with indirection and a keyword
C                               input scheme in mind:
C                               1. Avoid all backspacing (except at a higher
C                                  level for a namelist) by passing a meaningful
C                                  BUFFER (1) (1:LAST) in (and out).
C                               2. Reduce end-of-curve handling at this level
C                                  to simply identifying a nonnumeric 1st token.
C                                  The dictionary argument is thus redundant.
C                                  (Using LOOKUP on the first token of every
C                                  line was inefficient anyway: ALPHA can detect
C                                  a nonnumeric string more cheaply.)
C                               3. Handling of arbitrary text in column 1 with
C                                  XCOL, YCOL > 1 appears impossible for a
C                                  keyword scheme - no way of telling if the
C                                  text is a misspelled keyword.  However,
C                                  the special case of *** in column 1 is
C                                  allowed here.
C
C     Author:  Robert Kennelly, NASA Ames/Sterling Software, Palo Alto, CA.
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT NONE

C     Arguments.

      LOGICAL
     >   FINIS
      INTEGER
     >   LAST, LREAD, LWRITE, MAXPT, MAXCUR, POINTS (MAXCUR), XCOL, YCOL
      REAL
     >   X (MAXPT), Y (MAXPT)
      CHARACTER
     >   BUFFERS (2) * (*)

C     Local constants.

      INTEGER
     >   MAXCOL, MAXTOK
      CHARACTER
     >   EXCLAM * 1
      PARAMETER
     >  (EXCLAM = '!',
     >   MAXCOL = 30,   ! Highest accessible column number.
     >   MAXTOK = 40)   ! Largest data field width which may be read.

C     Local variables.

      INTEGER
     >   BASE, CURVE, ERROR, INSERT, IOS, NUMCOL, POINT, XYMAX
      LOGICAL
     >   ALPHA, CONVER, DONE, EOF, FATAL, LDSYN, NUMERIC, OK
      CHARACTER
     >   CONCAT * (2 * MAXTOK), LIST (0:MAXCOL) * (MAXTOK), XYFMT * 11

C     Procedures.

      EXTERNAL
     >   ALPHA, GETLINE, IOCHEK, TOKENS

C     Storage.

      DATA
     >   XYFMT /'(BN,2Fnn.0)'/
      SAVE
     >   BASE, XYFMT


C     Execution.
C     ----------

      ERROR = 0
      POINT = 0
      CURVE = CURVE + 1
      FINIS = .FALSE.
      FATAL = .FALSE.

C     Variable format in portable form:

      WRITE (XYFMT (7:8), '(I2)') MAXTOK

C     (Re)set local pointer to beginning of storage for this curve.

      IF (CURVE .EQ. 1) THEN
         BASE = 1
      ELSE
         BASE = BASE + POINTS (CURVE - 1)
      END IF

C     We can't check for too many curves (or pts.) till we're sure we have data.

C     Check specified column numbers.

      XYMAX = MAX (XCOL, YCOL)
      IF (XYMAX .GT. MAXCOL) THEN

C        Too many columns requested for present array dimensions.

         IF (LWRITE .GT. 0) THEN
            WRITE (LWRITE, 1030)
            WRITE (LWRITE, 1040) XCOL, YCOL
         END IF

         ERROR = -1
         GO TO 800
      END IF

      IF (XCOL .LE. 0 .AND. YCOL .LE. 0) THEN

C        Nonsense input values.

         IF (LWRITE .GT. 0) THEN
            WRITE (LWRITE, 1050) CURVE
            WRITE (LWRITE, 1040) XCOL, YCOL
         END IF

         ERROR = -1
         GO TO 800
      END IF


C     Look for a curve.
C     -----------------

      DONE = .FALSE.
  200 CONTINUE
         POINT = POINT + 1

         IF (POINT .GT. 1 .OR. LAST .EQ. 0) THEN   ! Read another line.

            CALL GETLINE (LREAD, EXCLAM, BUFFERS (1), LAST, IOS)
            CALL IOCHEK (IOS, OK, EOF, LDSYN, CONVER)

         ELSE             ! The calling program passes a valid buffer initially.
            OK = .TRUE.   ! This is simpler than avoiding the next test.
         END IF

         IF (.NOT. OK) THEN
            POINT = POINT - 1
            DONE = .TRUE.

            IF (EOF) THEN     ! We're done (EOF can terminate last curve).
               FINIS = .TRUE.
            ELSE              ! Fatal read error.
               FATAL = .TRUE.
            END IF

         ELSE

C           Analyze the input line.
C           -----------------------

            IF (LAST .EQ. 0) THEN   ! Skip a blank line.  
               POINT = POINT - 1
               GO TO 200
            END IF

            NUMCOL = XYMAX
            CALL TOKENS (BUFFERS (1) (1:LAST), NUMCOL, LIST (1))

            IF (ALPHA (LIST (1) (1:3))) THEN   ! First 3 chars. are not numeric.

C              The token is MOST LIKELY a control keyword.  Return to caller.

C              SPECIAL CASE:  Allow a bunch of ***s in column 1 which will
C                             be ignored if XCOL, YCOL are not 1.

               NUMERIC = LIST (1) (1:3) .EQ. '***'

               IF (.NOT. NUMERIC) THEN   ! Signal some kind of keyword
                  ERROR = +1
                  POINT = POINT - 1               
                  DONE  = .TRUE.
               END IF
C
C              DISTURBING REALIZATION:
C
C              The possibility of having arbitrary text in column 1 along with
C              other columns of numeric data specified by XCOL and YCOL cannot
C              be handled here.  Even at the higher level, there is no way of
C              distinguishing this case from the case of a misspelled keyword
C              followed by (possibly numeric) keyword value(s).
              

            ELSE    ! The first token cannot be a keyword.

               NUMERIC = .TRUE.

            END IF

            IF (NUMERIC) THEN

C              Are there enough columns for the numeric data?

               IF (NUMCOL .LT. XYMAX) THEN
              
                  IF (XCOL .EQ. 0 .AND. NUMCOL .EQ. 1) THEN

C                    Special case for a single column of Y-data where
C                    QPLOT's default YCOL = 2 is corrected here.  [Note
C                    that a single-column row among 2-column rows in
C                    this case would give spurious results.  DAS, 11/16/91]

                     YCOL = 1
                  ELSE
                     IF (LWRITE .GT. 0) THEN
                        WRITE (LWRITE, 1060) 
     >                     CURVE, NUMCOL, BUFFERS (1) (1:60)
                        WRITE (LWRITE, 1020) POINT, CURVE
                     END IF
                     POINT = POINT - 1
                     ERROR = -1
                     DONE = .TRUE.
                  END IF
               END IF
            END IF


            IF (.NOT. DONE) THEN

               IF (POINT .EQ. 1) THEN         ! Don't test curve # every pt. #.
                  IF (CURVE .GT. MAXCUR) THEN

C                     We've reached the curve limit but not the end of the data.

                      IF (LWRITE .GT. 0) WRITE (LWRITE, 1010) CURVE
                      POINT = POINT - 1
                      ERROR = -1
                      GO TO 800
                  END IF
               END IF

C              Read the data internally in floating point format.
C              --------------------------------------------------

               INSERT = BASE + POINT - 1
               IF (INSERT .LE. MAXPT) THEN

C                 One or the other of XCOL, YCOL may be zero.

                  IF (XCOL * YCOL .EQ. 0)
     >               WRITE (LIST (0), XYFMT) FLOAT (POINT)

C                 "Read" both data points at once to save error handling.

                  CONCAT = LIST (XCOL) // LIST (YCOL)
                  READ (CONCAT, XYFMT, IOSTAT = IOS)
     >               X (INSERT), Y (INSERT)
                  CALL IOCHEK (IOS, OK, EOF, LDSYN, CONVER)
                  IF (.NOT. OK) THEN
                     IF (CONVER) THEN
                        ERROR = +2
                        WRITE (BUFFERS (2), 1025) POINT, CURVE
                     ELSE
                        FATAL = .TRUE.
                     END IF
                     POINT = POINT - 1
                     DONE = .TRUE.
                  END IF
               ELSE

C                 We have run out of array space.

                  IF (LWRITE .GT. 0) THEN
                     WRITE (LWRITE, 1070)
                     WRITE (LWRITE, 1020) POINT, CURVE
                  END IF
                  POINT = POINT - 1
                  ERROR = -1
                  DONE = .TRUE.
               END IF
            END IF
         END IF

         IF (.NOT. DONE) GO TO 200     ! Look for more data


C     Termination.
C     ------------

  800 CONTINUE
      IF (FATAL) THEN

C        Abnormal termination due to fatal read error - bail out!

         IF (LWRITE .NE. 0) THEN
            WRITE (ABS (LWRITE), 1080) IOS
            WRITE (ABS (LWRITE), 1020) POINT + 1, CURVE
         END IF
         CALL EXIT (3)
      ELSE

C        Check whether last attempt actually read any data and adjust
C        counters accordingly.

         IF (POINT .GT. 0) THEN
            POINTS (CURVE) = POINT
         ELSE
            CURVE = CURVE - 1
            IF (CURVE .GT. 0) BASE = BASE - POINTS (CURVE)
         END IF
      END IF

      RETURN


C     Formats.
C     --------

 1010 FORMAT (/' READCOLS: Trouble - too many curves for array size.'/
     >   11X, '(curve number ', I5, ')')
 1020 FORMAT (11X, '(point number ', I5, ', curve number ', I5, ')')
 1025 FORMAT ('(point number ', I5, ', curve number ', I5, ')')
 1030 FORMAT (/' READCOLS: Trouble - too many data columns requested.')
 1040 FORMAT (11X, '(XCOL = ', I3, ', YCOL = ', I3, ')')
 1050 FORMAT (/' READCOLS: Trouble - illegal input column data for ',
     >   'curve number ', I5)
 1060 FORMAT (/' READCOLS: Trouble - not enough columns of input data ',
     >   'for curve'/
     >   11X, 'number ', I5, '.  Only ', I3, ' column(s) were found !'/
     >   11X, 'The questionable data began with:'/
     >   11X, '>>', A, '<<')
 1070 FORMAT (/' READCOLS: Trouble - too many data points.  Terminate '
     >   'this frame.')
 1080 FORMAT (/' READCOLS: Trouble - fatal READ error, IO status  = ',
     >   I6, '.')

      END
