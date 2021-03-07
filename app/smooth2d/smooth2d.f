C+----------------------------------------------------------------------
C
      PROGRAM SMOOTH2D
C
C  PURPOSE:
C
C     SMOOTH2D is a driver for 2-D interpolation/smoothing routines. It
C     applies the selected methods to one or more datasets read from a
C     3-column-type file.  The dataset(s) may or may not be rectangular
C     but some of the methods assume rectangularity and are expected to
C     be used appropriately.  Some also just require pseudorectangularity.
C
C     SMOOTH2D evaluates the fit at points which may be generated or read
C     from a separate file (in a similar format except only two columns
C     are needed - any additional columns are ignored).
C
C     Results are written to a file in the same form as the original data.
C
C     If (X,Y)s are generated, one slice or more in each direction is
C     provided for, and the rectangular grid of interpolated data may be
C     ordered in the output file by either rows or columns to assist 1-D
C     plotting of slices.  This format is also suited to the 2-D graphics
C     utility FMAP.
C
C     Dataset format:            <Reads are list-directed>
C
C     Title                      <Up to 80 characters>
C     NX  NY   or   NPTS         <2 integers = (pseudo)rectangular; 1 = random>
C     X(1)     Y(1)     F(1)     <If (pseudo)rectangular, the table may be
C     X(2)     Y(2)     F(2)      entered either by columns or by rows.>
C       :        :        :
C       :        :        :
C     X(N)     Y(N)     F(N)     <where N is either NX * NY or NPTS>
C     NX  NY   or   NPTS         <repeat for further datasets in same file>
C       :        :        :
C
C     An optional format for binary, pseudorectangular datasets is indicated
C     by the following code segment (basically all Xs by rows or columns, then
C     all Ys the same way, and all in one unformatted record).  The
C     corresponding function values are read similarly from a separate file:
C
C     READ (LUNGRD, IOSTAT=IOS) NX, NY
C     NPTS = NX * NY
C     READ (LUNGRD, IOSTAT=IOS) (X(I), I=1,NPTS), (Y(I), I=1,NPTS)
C
C     READ (LUNFUN, IOSTAT=IOS) NX, NY
C     READ (LUNFUN, IOSTAT=IOS) (F(I), I=1,NPTS)
C
C  METHOD:
C
C     Open file for normal evaluated results (smooth2d.out).
C     Prompt for and open the input data file; read its title.
C     Prompt for descriptive text to go to results file (default = title).
C     Prompt for generating (X,Y)s for evaluation, or reading them.
C     IF <read points> THEN
C        Prompt for the evaluation file.  ! Same for all datasets.
C        Open and read it.
C     END IF
C     Prompt for scaling the data so that distances between pts. are meaningful.
C     PICK_METHOD:
C        Select a 2-D method.
C        IF <CR> or ^Z GO TO <DONE>
C        Rewind the input data file; skip its title.
C        GET_DATASET:
C           Try to read a dataset (formatted or unformatted; pseudorect. or not)
C           IF <EOF> GO TO PICK_METHOD.
C           Determine data range (partly for scaling, partly for grid gen.).
C
C           IF <generate points> THEN
C              IF <first dataset and first method> THEN
C                 Display data range.
C                 Prompt for low, high, and increment in X and in Y.
C              END IF
C              Generate rectangular evaluation mesh (possibly just a slice),
C              ordered by rows or columns according to a first-time prompt.
C           END IF
C
C           IF <first dataset> THEN
C              Prompt for smoothing parameters (if applicable).
C           END IF
C           Fit requested surface.                |  These two steps may be
C           Evaluate fit at requested points.     |  combined by some methods.
C           Write results to output file.
C        GO TO GET_DATASET.
C     <DONE>
C
C  PARAMETER CONSTANTS:
C     NAME  TYPE   DESCRIPTION
C    MXMENU   I    Limit on no. of smoothing options
C    MXDPTS   I    Limit on no. of points in the input dataset (most methods)
C    MXIPTS   I    Limit on no. of points at which to interpolate
C    MXIWRK   I    Amount of integer workspace provided
C    MXNSQR   I    Limit on no. of data points for methods involving N*N matrix
C    MXRWRK   I    Amount of real workspace provided
C
C  SIGNIFICANT LOCAL VARIABLES:
C     NAME     DIM     TYPE DESCRIPTION
C    X,Y,F    MXDPTS    R   One set of data
C    XEVAL,   MXIPTS    R   Points for evaluating smoothed data
C    YEVAL
C    FEVAL    MXIPTS    R   Evaluated results
C    FX,FY,FXY NX*NY    R   Partial derivatives used by BIMOND3
C    IWORK    MXIWRK    I   Integer workspace needed (Akima; PBHMD/PBHEV)
C    RWORK    MXRWRK    R   Real workspace needed (most methods)
C
C  FILES USED:
C     LUN     I/O/S   DESCRIPTION
C    LUNCRT     O     For prompts to user, and for error messages
C    LUNDAT     I     Data to be smoothed (one or more datasets);
C                     X, Y, and F triples if formatted, else just the
C                     X, Y coordinates for the unformatted case
C                     (see LUNFUN in this case)
C    LUNFUN     I     Unformatted function data file, if specified
C    LUNKBD     I     For user responses to prompts
C    LUNOUT     O     For evaluated results (same as for LUNDAT if
C                     formatted input; as for LUNFUN if unformatted)
C    LUNXYS     I     For input coordinates requiring evaluation;
C                     either as for first two columns of formatted
C                     data file, or as for the unformatted grid option
C
C  PROCEDURES:
C
C    BOUNDS  NUMODULES  Finds min. and max. of 1-D or 2-D array
C    EVAL624 [.TOMS]    Evaluates triangulation of TRMESH (in TOMS624.FOR)
C    FACETS  INTERP2D   Bilinear interpolation (triangular facets)
C    GETLINE PRMODULES  Used to permit !-type trailing comment on NPTS line
C    HARDY   INTERP2D   Hardy's method (N*N dense system)
C    INMESH  INTERP2D   Locates given point within a pseudorectangular mesh
C    IQHSCV  INTERP2D   The 2-D Akima method (adapted from IMSL routines)
C    LCSFIT2D   "       Local bicubic spline method
C    OPENER  PRMODULES  File prompting/opening utility
C    PBHEV   [.BIMOND3] Piecewise bicubic Hermite evaluation
C    PBHMD   [.BIMOND3] Monotonic piecewise bicubic Hermite fit
C    PLBICUBE INTERP2D  Parametric form of LCSFIT2D
C    R1MACH  VMSLIB     Single precision machine-dependent constants
C    READER  PRMODULES  Prompting utility (single quantity)
C    REORDR  [.TOMS]    Reorders scattered data by X (in TOMS624.FOR)
C    SECOND  PRMODULES  Returns CPU time used by calling program so far
C    TABLE2  INTERP2D   Bilinear interpolation (table look-up)
C    TOKENS  PRMODULES  Needed to distinguish NX and NY from NPTS on line 2
C    TPSPLN  INTERP2D   Thin-plate spline, similar to Hardy's method
C    TRMESH  [.TOMS]    Triangulates scattered data (in TOMS624.FOR)
C
C  ENVIRONMENT:  VAX/VMS FORTRAN including
C    > the '$' carriage control for keeping cursor at end of prompt;
C    > some VAX/VMS dependencies in the OPENER utility.
C
C  HISTORY:
C    DAS   02/12/86   Initial adaptation from SMOOTH (Akima, Hardy,
C                     and TPSPLN).
C    DAS   02/14/86   Developed FACETS and installed it, after above
C                     methods seemed to behave dismally.
C    DAS   02/19/86   Allowed for <multiple methods/one dataset> in one
C                     run; added DETAILS for self-descriptive output.
C    DAS   03/14/86   Provided for generating grid of evaluation points.
C    DAS   03/21/86   Provided for scaling the Xs to match the Ys so that
C                     the "distances" used in some methods are meaningful.
C    DAS   09/19/86   MXNSQR limits certain methods to fewer data points
C                     than others, since some involve solution of N*N dense
C                     systems.
C    DAS   01/20/88   Installed BIMOND3=PBHMD/PBHEV along with RDTABLE
C                     derived from part of BIMOND3).  Reading of regular
C                     or scattered data (NX and NY, or just NPTS resp. on
C                     line 2) is automated using ERR= keyword.  Also R1MACH
C                     (which came with BIMOND3 and is called by PBHMD) is
C                     now called by SMOOTH2D because it needs "BIG" for its
C                     own purposes (data range finding).
C    DAS   03/03/88   ERR= with list-directed reads isn't good enough:
C                     looking for NX and NY where just NPTS is present
C                     means NY is set to (the truncated value of) X(1).
C                     Switched to counting tokens instead.
C    DAS   06/12/89   Provided for multiple datasets in one input file.
C    DAS   03/09/90   MXNSQR raised reluctantly from 100 to 300 for one user.
C    DAS   03/20/90   The Akima method was not saving its results when
C                     some points (probably only a few) were not evaluated
C                     properly (IER=999; extrapolation bug in the method).
C                     Results are now saved anyway - bulk of them may be OK.
C    DAS   08/21/90   Added INMESH/AKIMA combination for case of pseudo-
C                     rectangular dataset.  Also provided unformatted option
C                     for both the data file and the file of target (X,Y)s.
C                     Output results may also be unformatted now, and have
C                     either NXE, NYE, or just NEVAL if regular or not,
C                     although the unformatted case assumes regularity.
C                     Had to eliminate RDTABLE (for truly rectangular data)
C                     in order to handle pseudorectangular datasets.  (Test
C                     in-line for true rectangularity now.)  Raised MXDPTS,
C                     MXIPTS to 201*101 for a wing surface application.
C    DAS   04/09/91   Added LCSINT2D method for data consisting of rows
C                     with constant Y along each row.
C     "    04/18/91   LCSINT2D had independent edge definition arguments added
C                     (but they are not made use of here).  It also needed
C                     first/last indices for the evaluation arrays, and
C                     different interpolation controls for "rows" & "columns"
C                     and "leading" & "trailing" edges.
C     "    11/28/93   Added LCSFIT2D.
C     "    11/29/93   Added PLBICUBE.  Parametric cases are different enough
C                     that some data options such as rescaling may not be valid.
C     "    05/23/94   Installed TOMS Algorithm 624 (REORDR, TRMESH, EVAL624).
C
C AUTHOR:  David Saunders, NASA Ames/Sterling Software, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT    NONE

C ... Local constants:

      INTEGER
     >   LUNCRT, LUNDAT, LUNFUN, LUNKBD, LUNOUT, LUNXYS,
     >   MXDPTS, MXIPTS, MXIWRK, MXMENU, MXNSQR, MXRWRK
      REAL
     >   EPS
      CHARACTER
     >   BLANK * 1, BLANKS * 80

      PARAMETER
     >  (LUNDAT=1, LUNFUN=2, LUNXYS=3, LUNOUT=4, LUNKBD=5,
     >   LUNCRT=6, MXDPTS=201*201, MXIPTS=201*201, MXMENU=11,
     >   MXNSQR=300, MXIWRK=31 * MXDPTS + MXIPTS,
     >   MXRWRK=MAX (6 * MXDPTS, MXNSQR * (MXNSQR + 1)),
     >   EPS=1.E-6, BLANK=' ', BLANKS='  ')
C        Why did BLANKS have '+ ' originally???

C ... Local variables:

      INTEGER
     >   DATASET, I, IC, IER, IEVAL, IEX, IEXTRAP, IFOUND, IMAT, ITEMP,
     >   J, JC, JEVAL, JMAT, JTEMP, LAST, METHOD, NEVAL, NPTS, NTOKENS,
     >   NX, NXE, NY, NYE, IWORK(MXIWRK)

      REAL
     >   COND, DX, DY, P, Q, SCALE, SHAPE, SX, SY, TENS, TIME1, TIME2,
     >   UMAX, UMIN, VMAX, VMIN, XI, XMIN, XMAX, Y1J, YJ, YMIN, YMAX,
     >   FX(MXDPTS), FY(MXDPTS), FXY(MXDPTS), FEVAL(MXIPTS),
     >   FXEVAL(MXIPTS), FYEVAL(MXIPTS), FZEVAL(MXIPTS), RWORK(MXRWRK),
     >   UD(MXDPTS), VD(MXDPTS),
     >   X(MXDPTS), Y(MXDPTS), F(MXDPTS),
     >   XCELL(4), YCELL(4), FCELL(4), FMATRIX(4,4,3),
     >   XEVAL(MXIPTS), YEVAL(MXIPTS), XEORIG(MXIPTS),
     >   XYZU(3), XYZV(3), YVECTOR(MXDPTS/2 + 1) ! Overestimate...

      LOGICAL
     >   BINARYD, BINARYI, BYROWS, CR, EOF, FIRSTFIT,
     >   FIRSTPASS, FIRSTSET, GENXYS, MONINX, MONINY, PARAMETRIC,
     >   RECTREQD, REGULAR, REGULARI, SCALED, SEMIRECT, YES

      CHARACTER
     >   COLMTHD*1, DATAFILE*80, DETAILS*50, LEDGMTHD*1,
     >   LINE*60, LIST(2)*5, OUTFILE*20, MENU(0:MXMENU)*63,
     >   ROWMTHD*1, TEDGMTHD*1, TITLE*80

C ... COMMON blocks:

      COMMON /PBHCOM/ BIG, LOUT, LTTY, IPRINT     ! Needed by PBHMD/PBHEV
      REAL        BIG
      INTEGER     LOUT, LTTY, IPRINT

C ... Procedures:

      REAL
     >   R1MACH, TABLE2
      LOGICAL
     >   INMESH

      EXTERNAL
     >   BOUNDS, FACETS, GETLINE, HARDY, IQHSCV, INMESH, LCSFIT2D,
     >   LCSINT2D, OPENER, PBHEV, PBHMD, PLBICUBE,
     >   READC, READI, READR, READS, READY, R1MACH, SECOND, TABLE2,
     >   TOKENS, TPSPLN

C ... Storage:

      DATA        MENU   /
     >' 0. Help                                                       ',
     >' 1. AKIMA   : 5th deg. polynom. in each triangle; C1 continuity',
     >' 2. HARDY   : basis = SQRT ((X-XD)**2 + (Y-YD)**2 + shape**2)  ',
     >' 3. TPSPLN  : basis = LOG (D) * D**2;  D**2 = ((X-XD)**2 + ...)',
     >' 4. FACETS  : bilinear - nearest 2 pts. + 3rd to enclose target',
     >' 5. TABLE2  : bilinear - rectangular table look-up             ',
     >' 6. BIMOND  : monotonic bicubic Hermite - rectangular (Fritsch)',
     >' 7. INMESH  : quasi-rectangular search + AKIMA on 1 quadrilatl.',
     >' 8. LCSINT2D: semi-rectangular data - local 1-D spline methods.',
     >' 9. LCSFIT2D: quasi-rectangular data - local bicubic spline    ',
     >' 10.PLBICUBE: parametric form of LCSFIT2D - local bicubic      ',
     >' 11.TOMS624 : TOMS Algorithm 624 by R. Renka (scattered data)  '/

C ... NOTE: MENU(5:12) shows up in "details" portion of output.

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C                         Initialization
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

C     The PBHMD/PBHEV modules extracted from BIMOND3 use an internal COMMON:

CCCCC BIG = R1MACH (2)         ! Largest single precision magnitude
      BIG = 1.E+36             ! Largest single precision magnitude
      LOUT = LUNCRT            ! Diagnostics only if IPRINT = -1
      LTTY = LUNCRT            ! Actually used by BIMOND3, not PBHMD/PBHEV
      IPRINT = -1              ! Suppresses everything except error messages;
                               ! Use 0, 1, 2 for details of computations

C ... Prompt for and open the input data file:

      WRITE (LUNCRT, '(A)')
     >   BLANK,
     >   ' SMOOTH2D: Interpolation in 2-space by a choice of methods.',
     >   '           (No 2-D smoothing techniques are available yet).',
     >   BLANK,
     >   ' Some methods expect quasi-, semi-, or true rectangularity.',
     >   ' Use them appropriately.',
     >   ' Parametric methods cannot evaluate at specified (X,Y)s -',
     >   ' only at generated (U,V)s.',
     >   BLANK

      YES = .TRUE.
      CALL READY (LUNCRT, 'Is the input dataset formatted? ' //
     >   '(<CR>=Yes=(X,Y,F) triples; No=unformatted) ',
     >   LUNKBD, YES, CR, EOF)
      BINARYD = .NOT. YES

      OUTFILE = 'smooth2d.out'

      IF (.NOT. BINARYD) THEN

C ...    Open formatted file for output results:

         CALL OPENER ( LUNCRT, BLANK, LUNKBD, OUTFILE, LUNOUT,
     >                 'UNKNOWN')

         DATAFILE = 'smooth2d.dat'
         CALL OPENER (LUNCRT,
     >      'Enter the input data file name [default: smooth2d.dat]',
     >      LUNKBD, DATAFILE, LUNDAT, 'OLD')

         READ (LUNDAT, 1001, ERR=800) TITLE
         WRITE (LUNCRT, 1002)
     >      'Input dataset title:', TITLE(1:LEN_TRIM(TITLE))
         CALL READS (LUNCRT,
     >      'Enter title for evaluated results [default: as above]',
     >      LUNKBD, TITLE, CR, EOF)
         IF (EOF) GO TO 980

      ELSE

C ...    Open unformatted file for output results:

         CALL OPENER (LUNCRT, BLANK, LUNKBD, OUTFILE, LUNOUT,
     >      'UNKNOWN:UNFORMATTED')

         DATAFILE = BLANK
         CALL OPENER (LUNCRT,
     >      'Enter the input (X,Y) file name: ',
     >      LUNKBD, DATAFILE, LUNDAT, 'OLD:UNFORMATTED')

         DATAFILE = BLANK
         CALL OPENER (LUNCRT,
     >      'Enter the input function file name: ',
     >      LUNKBD, DATAFILE, LUNFUN, 'OLD:UNFORMATTED')

      END IF

C ... Generate or read evaluation point coordinates?
C     (Same points for all datasets.)

      WRITE (LUNCRT, '(/, A)')
     >   ' Do you want to generate points at which to interpolate?'
      GENXYS = .TRUE.
      CALL READY (LUNCRT,
     >   'Yes=generate; No=read them from a file instead; <CR>=Y: ',
     >   LUNKBD, GENXYS, CR, EOF)
      IF (EOF) GO TO 980

      IF (.NOT. GENXYS) THEN

         DATAFILE = 'xy.dat'
         YES = .TRUE.
         CALL READY (LUNCRT, 'Are the target pts. formatted? ' //
     >      '(<CR>=Yes=formatted pairs; No=unformatted) ',
     >      LUNKBD, YES, CR, EOF)
         BINARYI = .NOT. YES

         IF (.NOT. BINARYI) THEN           ! Formatted file

            CALL OPENER (LUNCRT,
     >         'Target point file name?  [default is xy.dat]',
     >         LUNKBD, DATAFILE, LUNXYS, 'OLD')

            READ (LUNXYS, 1001, ERR=840)         ! Skip the title

C ...       Allow for NPTS or NX and NY:

            CALL GETLINE (LUNXYS, '!', LINE, LAST, IER)
            IF (IER .LT. 0) GO TO 700                     ! EOF
            IF (IER .NE. 0  .OR.  LAST .EQ. 0) GO TO 810  ! Error

            NTOKENS = 2
            CALL TOKENS (LINE (1:LAST), NTOKENS, LIST)
            IF (NTOKENS .EQ. 0  .OR.  NTOKENS .GT. 2) GO TO 810

            READ (LIST (1), '(BN, I5)', ERR=810) NXE

            REGULARI = NTOKENS .EQ. 2
            IF (REGULARI) THEN                    ! (Pseudo)rectangular data
               READ (LIST (2), '(BN, I5)', ERR=810) NYE
               NEVAL = NXE * NYE
            ELSE                                  ! Random data
               NEVAL = NXE
            END IF

            IF (NEVAL .LT. 1  .OR.  NEVAL .GT. MXIPTS) GO TO 850

            IFOUND = 0
            DO 50 I = 1, MXIPTS
               READ (LUNXYS, *, ERR=860, END=60) XEVAL(I), YEVAL(I)
               IFOUND = I
   50       CONTINUE
 
   60       IF (IFOUND .NE. NEVAL) THEN
               WRITE (LUNCRT, 1004)
     >            ' WARNING: # target pts. found: ', IFOUND,
     >            '          Number expected:     ', NEVAL,
     >            '          Proceeding...'
               NEVAL = IFOUND
            END IF

         ELSE          ! Unformatted target points (pseudorectangular only)

            CALL OPENER (LUNCRT,
     >         'Target point file name?  [default is xy.dat]',
     >         LUNKBD, DATAFILE, LUNXYS, 'OLD:UNFORMATTED')

            READ (LUNXYS, ERR=850) NXE, NYE
            NEVAL = NXE * NYE
            IF (NEVAL .LT. 1  .OR.  NEVAL .GT. MXIPTS) GO TO 850

            READ (LUNXYS, ERR=860) (XEVAL(I), I=1, NEVAL),
     >                             (YEVAL(I), I=1, NEVAL)
         END IF

      ELSE
C ...    Generate the target points once the data range is known.
      END IF


C ... SCALED = .TRUE. shall mean each dataset is already well scaled:

      WRITE (LUNCRT, '(/, A)')
     >   ' Most methods need similar X & Y scaling.'
      YES = .FALSE.
      CALL READY (LUNCRT,
     >   'Do you want SMOOTH2D to do any scaling? (Y/N <CR>=N): ',
     >   LUNKBD, YES, CR, EOF)
      IF (EOF) GO TO 980
      SCALED = .NOT. YES


C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C                    Top of loop over different methods
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      FIRSTFIT = .TRUE.
  100 CONTINUE

         DETAILS = BLANK
         WRITE (LUNCRT, 1001)
         CALL READI (LUNCRT,
     >      'Select a 2-D method, 0 for help, or <CR> when finished: ',
     >      LUNKBD, METHOD, CR, EOF)
         IF (CR .OR. EOF) GO TO 990

         IF (METHOD .LE. 0 .OR. METHOD .GT. MXMENU) THEN
            WRITE (LUNCRT, 1001) MENU
            GO TO 100
         END IF

         RECTREQD = METHOD .EQ. 5 .OR. METHOD .EQ. 6   ! TABLE2 or BIMOND
         SEMIRECT = METHOD .EQ. 8                      ! LCSINT2D
         PARAMETRIC = METHOD .EQ. 10                   ! PLBICUBE

         REWIND LUNDAT
         IF (.NOT. BINARYD) THEN
            READ (LUNDAT, 1001)    ! Skip the title
         ELSE
            REWIND LUNFUN
         END IF


C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C        Top of loop over (unknown no. of) datasets in the input file
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

         DATASET = 0
         FIRSTSET = .TRUE.
  200    CONTINUE

            FIRSTPASS = FIRSTSET .AND. FIRSTFIT
            DATASET = DATASET + 1

C ...       Read a dataset (2 formatted cases; 1 unformatted case):

            IF (.NOT. BINARYD) THEN          ! Formatted data

C ...          Allow for NPTS or NX and NY:

               CALL GETLINE (LUNDAT, '!', LINE, LAST, IER)
               IF (IER .LT. 0) GO TO 700                     ! EOF
               IF (IER .NE. 0  .OR.  LAST .EQ. 0) GO TO 810  ! Error

               NTOKENS = 2
               CALL TOKENS (LINE (1:LAST), NTOKENS, LIST)
               IF (NTOKENS .EQ. 0  .OR.  NTOKENS .GT. 2) GO TO 810

               READ (LIST (1), '(BN, I5)', ERR=810) NX

               REGULAR = NTOKENS .EQ. 2
               IF (REGULAR) THEN                     ! (Pseudo)rectangular data
                  READ (LIST (2), '(BN, I5)', ERR=810) NY
                  IF (FIRSTPASS)
     >               WRITE (LUNCRT, '(/, '' NX, NY found: '', 2I5)')
     >                  NX, NY
                  IF (NX .LT. 2  .OR.  NY .LT. 2) GO TO 815

                  NPTS = NX * NY
               ELSE                                  ! Random data
                  NPTS = NX
               END IF

               IF (NPTS .GT. MXDPTS) GO TO 820

               IFOUND = 0
               DO 210 I = 1, NPTS   ! Can't use MXDPTS if multiple datasets
                  READ (LUNDAT, *, ERR=830, END=220) X(I), Y(I), F(I)
                  IFOUND = I
  210          CONTINUE

  220          IF (IFOUND .NE. NPTS) THEN
                  WRITE (LUNCRT, 1004)
     >               ' WARNING: # data pts. found: ', IFOUND,
     >               '          Number expected:   ', NPTS,
     >               '          Proceeding with:   ', IFOUND,
     >               BLANK
                  NPTS = IFOUND
               END IF

            ELSE                ! Unformatted data (pseudorectangular only)

               REGULAR = .TRUE.
               READ (LUNDAT, END=700, ERR=870) NX, NY
               NPTS = NX * NY
               IF (NPTS .GT. MXDPTS) GO TO 820

               READ (LUNDAT, ERR=875) (X(I), I=1,NPTS),
     >                                (Y(I), I=1,NPTS)

               READ (LUNFUN, ERR=880) ITEMP, JTEMP
               IF (ITEMP .NE. NX .OR. JTEMP .NE. NY) GO TO 885

               READ (LUNFUN, ERR=887) (F(I), I=1,NPTS)

            END IF

C ...       Trap non-rectangular data if method expects it.
C ...       Assume storage by rows of constant Y.

            IF (RECTREQD) THEN
               IF (.NOT. REGULAR) GO TO 890

               DO 235, J = 1, NY
                  JTEMP = (J - 1) * NX
                  Y1J = Y (JTEMP + 1)
                  YVECTOR (J) = Y1J ! Rectangular methods expect Y(1:NY) only...
                  DO 230, I = 1, NX
                     ITEMP = JTEMP + I
                     IF (X (ITEMP) .NE. X (I)) GO TO 890
                     IF (Y (ITEMP) .NE. Y1J)   GO TO 890
  230             CONTINUE
  235          CONTINUE
            END IF

            IF (SEMIRECT) THEN
               IF (.NOT. REGULAR) GO TO 900

               DO 245, J = 1, NY
                  JTEMP = (J - 1) * NX
                  Y1J = Y (JTEMP + 1)
                  YVECTOR (J) = Y1J ! Rectangular methods expect Y(1:NY) only...
                  DO 240, I = 1, NX
                     ITEMP = JTEMP + I
                     IF (Y (ITEMP) .NE. Y1J) GO TO 900
  240             CONTINUE
  245          CONTINUE
            END IF

C ...       Determine data range:

            XMIN = X(1)
            XMAX = XMIN
            CALL BOUNDS (NPTS, 1, NPTS, X, XMIN, XMAX)

            YMIN = Y(1)
            YMAX = YMIN
            CALL BOUNDS (NPTS, 1, NPTS, Y, YMIN, YMAX)

            IF (YMIN .EQ. YMAX  .OR.  XMIN .EQ. XMAX) GO TO 835

            SCALE = (YMAX - YMIN) / (XMAX - XMIN)

C ...       Allow for generating only one evaluation grid for all datasets:

            IF (GENXYS) THEN

               IF (FIRSTPASS) THEN

                  WRITE (LUNCRT, '(/, (A, 1P, 2E12.3))')
     >               ' Data range: Xmin, Xmax = ', XMIN, XMAX,
     >               '             Ymin, Ymax = ', YMIN, YMAX

                  IF (PARAMETRIC) THEN

                    UMIN = 0.E+0
                    UMAX = 1.E+0

                    VMIN = 0.E+0
                    VMAX = 1.E+0

                    WRITE (LUNCRT, '(/, (A, 1P, 2E12.3))')
     >                 ' U, V range: Umin, Umax = ', UMIN, UMAX,
     >                 '             Vmin, Vmax = ', VMIN, VMAX
                  END IF

                  WRITE (LUNCRT, 1001) ' Enter definition of grid' //
     >               ' points at which to interpolate.'

C ...             Allow one or more slices each direction:

  250             WRITE (LUNCRT, '(A)', ADVANCE='NO')
     >               ' Min, max, inc for X (or U): '
                  READ  (LUNKBD, *, ERR=250, END=980) XMIN, XMAX, DX
                  IF (DX .EQ. 0.) DX = 1.  ! XMIN = XMAX is OK (one slice)

  260             WRITE (LUNCRT, '(A)', ADVANCE='NO')
     >               ' Min, max, inc for Y (or V): '
                  READ  (LUNKBD, *, ERR=260, END=980) YMIN, YMAX, DY
                  IF (DY .EQ. 0.) DY = 1.
               END IF

               NXE = 1 + INT ((XMAX-XMIN + DX * 0.01) / DX)
               NYE = 1 + INT ((YMAX-YMIN + DY * 0.01) / DY)
               NEVAL = NXE * NYE

               IF (NEVAL .GT. MXIPTS) THEN
                  WRITE (LUNCRT, 1004)
     >            ' Sorry - too many interpolation points: ', NEVAL,
     >            ' Try again.  The current limit is:      ', MXIPTS
                  GO TO 250
               END IF

               IF (FIRSTPASS) THEN
                  BYROWS = .TRUE.
                  CALL READY (LUNCRT, 'Save results by rows?' //
     >               ' (Y/N; No=columns) (<CR>=Yes=rows): ',
     >               LUNKBD, BYROWS, CR, EOF)
                  IF (EOF) GO TO 980
               END IF

               NEVAL = 0

               IF (BYROWS) THEN

                  DO 270 J = 1, NYE
                     YJ = YMIN + (J-1)*DY
                     IF (J .EQ. NYE) YJ = YMAX
                     DO 265 I = 1, NXE
                        NEVAL = NEVAL + 1
                        XEVAL(NEVAL) = XMIN + (I-1)*DX
                        YEVAL(NEVAL) = YJ
  265                CONTINUE
                     XEVAL(NEVAL) = XMAX
  270             CONTINUE

               ELSE  ! Generate grid by columns:

                  DO 280 I = 1, NXE
                     XI = XMIN + (I-1)*DX
                     IF (I .EQ. NXE) XI = XMAX
                     DO 275 J = 1, NYE
                        NEVAL = NEVAL + 1
                        XEVAL(NEVAL) = XI
                        YEVAL(NEVAL) = YMIN + (J-1)*DY
  275                CONTINUE
                     YEVAL(NEVAL) = YMAX
  280             CONTINUE

               END IF
            END IF


            IF (.NOT. SCALED) THEN

C ...          Scale the data so that distances make sense.
C              Save original values to avoid unscaling later.

               DO 310 I = 1, NPTS
                  X(I) = X(I) * SCALE
  310          CONTINUE

               DO 320 I = 1, NEVAL
                  XEORIG(I) = XEVAL(I)
                  XEVAL(I)  = XEVAL(I) * SCALE
  320          CONTINUE

            END IF


C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C              Apply current 2-D method to current dataset
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

            CALL SECOND (TIME1)

            GO TO (410, 420, 430, 440, 450, 460, 470, 480, 490, 500,
     >             510) METHOD


  410       CONTINUE   ! *** Akima's 2-D method ***

C ...       Interpolates data with a 5th-degree polynomial in each triangle
C           of a triangulation of the x-y projection of the surface.  Gives
C           continuity of function and first partial derivatives everywhere.
C           The data points may be irregularly distributed.

            IF (FIRSTSET) THEN
               IEXTRAP = 0
               YES = .TRUE.
               CALL READY (LUNCRT,
     >            'Do you want to permit extrapolation? (<CR>=Y): ',
     >            LUNKBD, YES, CR, EOF)
               IF (EOF) GO TO 980
               IF (.NOT. YES) THEN
                  IEXTRAP = 1
                  WRITE (LUNCRT, 1001)
     >      ' A value of 9999. means the target point was out of range.'
               END IF
           
  415          TENS = 0.0E+0
               CALL READR (LUNCRT,
     >            'Enter tension parameter in range [0,1] (<CR>=0.): ',
     >            LUNKBD, TENS, CR, EOF)
               IF (EOF) GO TO 980
               IF (TENS .LT. 0.E+0  .OR.  TENS .GT. 1.E+0) GO TO 415

               WRITE (DETAILS, '(''Tension ='', F7.3)') TENS
               WRITE (LUNCRT, 1001) '      Working ...'
            END IF

            CALL IQHSCV (X, Y, F, NPTS, XEVAL, YEVAL, FEVAL, NEVAL,
     >                   IWORK, RWORK, IER, IEXTRAP, TENS)
            IF (IER .EQ. 999) THEN
               WRITE (LUNCRT, 1002)
     >            BLANK,
     >            'Warning: Results contain some bad function values:',
     >            BLANK,
     >            '   9999. means out-of-range (X,Y)s but extrapolation'
     >            // ' was turned off;',
     >            '   9998. means requested extrapolation failed at'
     >            // ' some (X,Y)s (known bug);',
     >            '   9997. means tension term failure (unlikely).',
     >            BLANK
            ELSE IF (IER .NE. 0) THEN
               GO TO 910
            ELSE IF (FIRSTSET) THEN
               WRITE (LUNCRT, 1001) BLANKS
            END IF
            GO TO 600


  420       CONTINUE   ! *** Hardy's method ***

C ...       Bivariate interpolation: linear combination of basis functions
C           of the form
C
C              PHI(K,X,Y) = SQRT ((X-XD(K))**2 + (Y-YD(K))**2 + R**2)
C
C           where shape parameter R should be in the range [0,~2].
C           The data points may be irregular, but no two may be the same.

            IF (NPTS .GT. MXNSQR) GO TO 825

  425       IF (FIRSTSET) THEN
               SHAPE = 0.E+0
               CALL READR (LUNCRT,
     >            'Enter shape parameter in range [0,~2] (<CR>=0.>): ',
     >            LUNKBD, SHAPE, CR, EOF)
               IF (EOF) GO TO 980
               IF (SHAPE .LT. 0.E+0  .OR.  SHAPE .GT. 10.0) GO TO 425

               WRITE (LUNCRT, 1001) '      Working ...'
            END IF

            CALL HARDY (X, Y, F, NPTS, XEVAL, YEVAL, FEVAL, NEVAL,
     >                  IWORK, RWORK, SHAPE, IER, COND)
            IF (IER .NE. 0) GO TO 920

            IF (FIRSTSET) THEN
               WRITE (LUNCRT, 1012) ' Matrix condition number:', COND
            END IF

            WRITE (DETAILS, '(''Shape param:'', F7.3,
     >         ''   COND(Matrix):'', 1PE12.3)') SHAPE, COND
            GO TO 600


  430       CONTINUE   ! *** Thin-plate spline (similar to Hardy's method) ***

C ...       Bivariate interpolation: linear combination of basis functions
C           of the form
C
C              PHI(K,X,Y) = D(K)**2 * LOG (D(K))    where
C              D(K)**2 = (X-XD(K))**2 + (Y-YD(K))**2
C
C           The data points may be irregular, but no two may be the same.

            IF (NPTS .GT. MXNSQR) GO TO 825
            IF (FIRSTSET) WRITE (LUNCRT, 1001) '      Working ...'

            CALL TPSPLN (X, Y, F, NPTS, XEVAL, YEVAL, FEVAL, NEVAL,
     >                   IWORK, RWORK, IER, COND)
            IF (IER .NE. 0) GO TO 930

            IF (FIRSTSET) THEN
               WRITE (LUNCRT, 1012) ' Matrix condition number:', COND
            END IF
            WRITE (DETAILS, '(''COND(Matrix):'', 1PE12.3)') COND
            GO TO 600


  440       CONTINUE   ! *** FACETS: Use nearest 3 pts defining a plane... ***

            IF (FIRSTSET) WRITE (LUNCRT, 1001) '      Working ...'

            CALL FACETS (X, Y, F, NPTS, .TRUE.,
     >                   XEVAL, YEVAL, FEVAL, NEVAL, IWORK, IER)
            IF (IER .NE. 0) GO TO 940

            IF (FIRSTSET) WRITE (LUNCRT, 1001) BLANKS
            GO TO 600


  450       CONTINUE   ! *** 2-D Table Look-up ***

C ...       4-point formula for bilinear interpolation.

            IF (FIRSTSET) WRITE (LUNCRT, 1001) '      Working ...'

            DO 455 I = 1, NEVAL
               FEVAL(I) = TABLE2 (NX, NY, NX, X, YVECTOR, F,
     >                            XEVAL(I), YEVAL(I), IER)
  455       CONTINUE

            IF (FIRSTSET) WRITE (LUNCRT, 1001) BLANKS
            GO TO 600


  460       CONTINUE   ! *** BIMOND3: ***

C ...       Monotonic piecewise bicubic Hermite spline ***

C           No attempt to handle supplied boundary derivatives yet:

            DO 462, I = 1, 2 * MAX (NX, NY)
               IWORK(I) = 0
  462       CONTINUE

            IF (FIRSTSET) THEN
               MONINX = .FALSE.
               CALL READY (LUNCRT,
     >         'Do you want fit to match monotonicity in X? (<CR>=N): ',
     >         LUNKBD, MONINX, CR, EOF)
               IF (EOF) GO TO 980

               MONINY = .FALSE.
               CALL READY (LUNCRT,
     >         'Do you want fit to match monotonicity in Y? (<CR>=N): ',
     >         LUNKBD, MONINY, CR, EOF)
               IF (EOF) GO TO 980

  464          IEX = 0
               CALL READI (LUNCRT,
     >           'Extrapolate? <CR>=0=No; 1=in X; 2=in Y; 3=in X & Y: ',
     >         LUNKBD, IEX, CR, EOF)
               IF (EOF) GO TO 980
               IF (IEX .LT. 0 .OR. IEX .GT. 3) GO TO 464

               WRITE (LUNCRT, 1001) '      Working ...'
            END IF

C           Fit the interpolant:

            CALL PBHMD (NX, NY, NX, X, YVECTOR, F, IWORK, IWORK, MONINX,
     >                  MONINY, .TRUE., SX, SY, FX, FY, FXY, IER,
     >                  MXDPTS, RWORK)

C           Evaluate the fit at specified points:

            J = 0        ! ISTART argument for PBHEV.

            DO 466, I = 1, NEVAL
               CALL PBHEV (NX, NY, NX, NY, 1, J, X, YVECTOR, F, FX, FY,
     >                     FXY, XEVAL(I), YEVAL(I), FEVAL(I), IEX)
  466       CONTINUE

            IF (FIRSTSET) WRITE (LUNCRT, 1001) BLANKS
            GO TO 600


  470       CONTINUE   ! *** Hybrid method: search pseudo-rectangular
                       !     mesh for relevant cell; apply Akima's
                       !     method to the resulting 4 points ***

            IF (FIRSTSET) THEN
               TENS = 0.0E+0
               IEXTRAP = 1         ! Extrapolation is not meaningful.
               WRITE (LUNCRT, '(/, (a))')
     >     ' A value of 9999. means the target point was out of range.',
     >            '      Working ...'
            END IF

            DO 475 I = 1, NEVAL

               IF (INMESH (XEVAL (I), YEVAL (I), NX, NX, NY,
     >                     X, Y, IC, JC)) THEN

                  ITEMP = (JC - 1) * NX + IC
                  XCELL (1) = X (ITEMP)
                  YCELL (1) = Y (ITEMP)
                  FCELL (1) = F (ITEMP)
                  XCELL (2) = X (ITEMP - 1)
                  YCELL (2) = Y (ITEMP - 1)
                  FCELL (2) = F (ITEMP - 1)
                  ITEMP = ITEMP - NX
                  XCELL (3) = X (ITEMP)
                  YCELL (3) = Y (ITEMP)
                  FCELL (3) = F (ITEMP)
                  XCELL (4) = X (ITEMP - 1)
                  YCELL (4) = Y (ITEMP - 1)
                  FCELL (4) = F (ITEMP - 1)

                  CALL IQHSCV (XCELL, YCELL, FCELL, 4,
     >                         XEVAL (I), YEVAL (I), FEVAL (I), 1,
     >                         IWORK, RWORK, IER, IEXTRAP, TENS)
                  IF (IER .NE. 0) GO TO 910

               ELSE
                  FEVAL (I) = 9999.      ! Can't extrapolate
               END IF

  475       CONTINUE

            GO TO 600


  480       CONTINUE   ! *** Local 1-D cubic spline method for data in rows:
                       !     interpolate across rows at points of equal
                       !     RELATIVE position. ***

            IF (FIRSTSET) THEN
               LEDGMTHD = 'L'
               WRITE (LUNCRT, '(/, A)')
     >            ' Interpolation method for "leading" edge?'
               CALL READC (LUNCRT,
     >            '(Linear; Monotonic spline; Bessel-type ' //
     >            '(looser fit); <CR>=L) ', LUNKBD, LEDGMTHD, CR, EOF)

               TEDGMTHD = 'L'
               WRITE (LUNCRT, '(/, A)')
     >            ' Interpolation method for "trailing" edge?'
               CALL READC (LUNCRT, '(Same choices; <CR> = L) ',
     >            LUNKBD, TEDGMTHD, CR, EOF)

               ROWMTHD = 'B'
               WRITE (LUNCRT, 1001)' Interpolation method for "rows"?'
               CALL READC (LUNCRT, '(Same choices; <CR> = B) ',
     >            LUNKBD, ROWMTHD, CR, EOF)

               COLMTHD = 'B'
               WRITE(LUNCRT, 1001)' Interpolation method across "rows"?'
               CALL READC (LUNCRT, '(Same choices; <CR> = B) ',
     >            LUNKBD, COLMTHD, CR, EOF)

               WRITE (LUNCRT, 1001) '      Working ...'
            END IF

C           Too cumbersome to support the independent edge definitions here...

            CALL LCSINT2D (NX, NY, NX, NY, X, YVECTOR, F,
     >                     0, 0, 0, 0, 0, 0,
     >                     LEDGMTHD, TEDGMTHD, ROWMTHD, COLMTHD,
     >                     NXE, NYE, 1, NXE, 1, NYE, XEVAL,YEVAL,FEVAL)

            IF (FIRSTSET) WRITE (LUNCRT, 1001) BLANKS

            GO TO 600


  490       CONTINUE   ! *** LCSFIT2D: Local bicubic spline for regular data:

            WRITE (LUNCRT, 1001) '      Working ...'

            IEVAL = 1
            JEVAL = 1
            IMAT = 0   ! I.e., an invalid point first time in
            JMAT = 0

            DO I = 1, NEVAL
               CALL LCSFIT2D (0, NX, NY, 1, NX, 1, NY, X, Y, F,
     >                        XEVAL (I), YEVAL (I), IEVAL, JEVAL,
     >                        IMAT, JMAT, FMATRIX, EPS,
     >                        FEVAL (I), FXEVAL (1), FYEVAL (1), IER)
               IF (IER .NE. 0) THEN
                  WRITE (LUNCRT, 1004) ' Out-of-range point number: ', I
                  FEVAL (I)  = 9999.
               END IF
            END DO

            IF (FIRSTSET) WRITE (LUNCRT, 1001) BLANKS

            GO TO 600


  500       CONTINUE   ! *** PLBICUBE: Parametric local bicubic spline:

            WRITE (LUNCRT, 1001) '      Working ...'

            IEVAL = 1
            JEVAL = 1
            UD (1) = -1.      ! Forces parameterization of the surface

            DO I = 1, NEVAL
               CALL PLBICUBE (0, NX, NY, 1, NX, 1, NY, X, Y, F, UD, VD,
     >                        XEVAL (I), YEVAL (I), IEVAL, JEVAL,
     >                        IMAT, JMAT, FMATRIX, EPS, P, Q,
     >                        FXEVAL (I), FYEVAL (I), FZEVAL (I),
     >                        XYZU, XYZV, IER)
               IF (IER .NE. 0) THEN
                  WRITE (LUNCRT, 1004) ' Out-of-range point number: ', I
                  FXEVAL (I) = 9999.
                  FYEVAL (I) = 9999.
                  FZEVAL (I) = 9999.
               END IF
            END DO

            IF (FIRSTSET) WRITE (LUNCRT, 1001) BLANKS

            GO TO 600


  510       CONTINUE   ! *** Renka's method (TOMS Algorithm 624) ***

            WRITE (LUNCRT, '(A)')
     >         'An early TOMS 624 source has been lost.  The current',
     >         'TOMS 624 source is available as /SMOOTH2D/TOMS624.f',
     >         'but some effort is needed to make the switch.'

C ...       Triangulates scattered data and interpolates it at specified
C           points with first-derivative continuity.  Benefits from
C           sorting the data in X.

C           YES = .TRUE.
C           CALL READY (LUNCRT,
C    >         'Sort the data in-place by X? (Recommended; <CR>=Y): ',
C    >         LUNKBD, YES, CR, EOF)
C           IF (EOF) GO TO 980

C           IF (YES) THEN
C              CALL REORDR (NPTS, 3, X, Y, F, IWORK)
C           END IF
           
C           Triangulate the data:

C           CALL TRMESH (NPTS, X, Y, IWORK (NPTS + 1), IWORK (1), IER)
C ...                                IADJ (6*NPTS)     IEND (NPTS)

C           Interpolate, specifying that derivatives be generated once
C           in RWORK for the first target point and reused for the rest.

C           CALL EVAL624 (NPTS, X, Y, F, IWORK (NPTS + 1), IWORK (1),
C    >                    NEVAL, XEVAL, YEVAL, 2, RWORK, FEVAL, IER)
C ...                                             ZXZY (2, NPTS)
C           IF (IER .LT. 0) GO TO 950
C           IF (IER .GT. 0) THEN
C              WRITE (LUNCRT, '(/, A, I5, A, /)')
C    >            ' *** Warning:', IER, ' extrapolations performed.'
C           ELSE IF (FIRSTSET) THEN
C              WRITE (LUNCRT, 1001) BLANKS
C           END IF
            GO TO 600


  600       CONTINUE

C ...       End of case statement over smoothing methods.

            CALL SECOND (TIME2)

            TIME2 = TIME2 - TIME1
            IF (FIRSTSET) THEN
               WRITE (LUNCRT, 1006) ' Elapsed CPU secs:', TIME2
               WRITE (LUNCRT, 1001) '      Saving results ...'
            END IF

C ...       Save evaluations:

            IF (FIRSTPASS .AND. .NOT. BINARYD)
     >         WRITE (LUNOUT, 1002) TITLE

            IF (GENXYS) THEN
               IF (BYROWS) THEN
                  IF (.NOT. BINARYD) THEN
                     WRITE (LUNOUT, 1008) NXE, NYE, MENU(METHOD)(5:12),
     >                                    DETAILS
                  ELSE
                     WRITE (LUNOUT) NXE, NYE
                  END IF
               ELSE
                  IF (.NOT. BINARYD) THEN
                     WRITE (LUNOUT, 1008) NYE, NXE, MENU(METHOD)(5:12),
     >                                    DETAILS
                  ELSE
                     WRITE (LUNOUT) NYE, NXE
                  END IF
               END IF
            ELSE
               IF (.NOT. BINARYD) THEN
                  IF (REGULARI) THEN
                     WRITE (LUNOUT, 1008) NXE, NYE, MENU(METHOD)(5:12),
     >                                    DETAILS
                  ELSE
                     WRITE (LUNOUT, 1010) NEVAL, MENU(METHOD)(5:12),
     >                                    DETAILS
                  END IF
               ELSE     ! Pseudorectangular only for unformatted case.
                  WRITE (LUNOUT) NXE, NYE
               END IF
            END IF

            IF (.NOT. BINARYD) THEN
               IF (PARAMETRIC) THEN
                  WRITE (LUNOUT, 1016)
     >               (FXEVAL(I), FYEVAL(I), FZEVAL(I), I=1, NEVAL)
               ELSE
                  IF (SCALED) THEN
                      WRITE (LUNOUT, 1016)
     >                  (XEVAL(I), YEVAL(I), FEVAL(I), I=1, NEVAL)
                  ELSE
                     WRITE (LUNOUT, 1016)
     >                  (XEORIG(I), YEVAL(I), FEVAL(I), I=1, NEVAL)
                  END IF
               END IF
            ELSE
               WRITE (LUNOUT) (FEVAL(I), I=1, NEVAL)
            END IF

            IF (FIRSTSET) WRITE (LUNCRT, 1001) BLANKS


            FIRSTSET = .FALSE.
         GO TO 200

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C              End of loop over datasets in the input file
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


  700    CONTINUE

C ...    Come here upon EOF encountered looking for another dataset.

         FIRSTFIT = .FALSE.
      GO TO 100      

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C           End of loop over methods applied to each dataset
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C                           Error handling:
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  800 WRITE (LUNCRT, 1001) ' Error reading the input dataset title.'
      GO TO 999
  810 WRITE (LUNCRT, 1001)
     >   ' Error reading the number of input points.',
     >   ' Reminder: line 1 in the file should be a title;',
     >   '           line 2 should contain NX and NY (regular data),',
     >   '                  or             NPTS    (irregular data).',
     >   '           Begin any trailing comment on line 2 with "!".'
      GO TO 999
  815 WRITE (LUNCRT, 1004) ' Too few data points found.  NX:', NX,
     >                     '                             NY:', NY
      GO TO 999
  820 WRITE (LUNCRT, 1004) ' Too many data points found: ', NPTS,
     >                     ' Maximum allowed:            ', MXDPTS
      GO TO 999
  825 WRITE (LUNCRT, 1004)
     >   ' Too many points for NxN matrix method: ', NPTS,
     >   ' Maximum N allowed:                     ', MXNSQR
      GO TO 999
  830 WRITE (LUNCRT, 1001) ' Error reading the input (X,Y,F) triples.'
      GO TO 999
  835 WRITE (LUNCRT, 1001) ' Dataset is only 1-D ... aborting.'
      GO TO 999
  840 WRITE (LUNCRT, 1001) ' Error reading the evaluation pts. title.'
      GO TO 999
  850 WRITE (LUNCRT, 1001) ' Error in the no. of evaluation points:'
      WRITE (LUNCRT, 1004) ' Minimum required =', 1,
     >                     ' Maximum allowed  =', MXIPTS,
     >                     ' Number input     =', NEVAL
      GO TO 999
  860 WRITE (LUNCRT, 1001) ' Error reading (X,Y) evaluation coords.'
      GO TO 999
  870 WRITE (LUNCRT, 1001) ' Error reading NX, NY in grid file.'
      GO TO 999
  875 WRITE (LUNCRT, 1001) ' Error reading (X,Y) data.'
      GO TO 999
  880 WRITE (LUNCRT, 1001) ' Error reading NX, NY in function file.'
      GO TO 999
  885 WRITE (LUNCRT, 1001)
     >   ' Function array dimensions do not match grid dimensions.'
      GO TO 999
  887 WRITE (LUNCRT, 1001) ' Error reading function data.'
      GO TO 999
  890 WRITE (LUNCRT, 1001) ' This method requires rectangular data.'
      GO TO 100
  900 WRITE (LUNCRT, 1001) ' This method requires rows of data.'
      GO TO 100
  910 WRITE (LUNCRT, 1004) ' Bad return from IQHSCV.  IER: ', IER
      GO TO 200
  920 WRITE (LUNCRT, 1004) ' Bad return from HARDY.  IER: ', IER
      GO TO 200
  930 WRITE (LUNCRT, 1004) ' Bad return from TPSPLN.  IER: ', IER
      GO TO 200
  940 WRITE (LUNCRT, 1001)
     >   ' FACETS failed on some point(s), set to 9999.',
     >   ' Probable cause: interpolation point outside data range.'
      GO TO 600
  950 WRITE (LUNCRT, 1004) ' Bad return from EVAL624. IER: ', IER
      GO TO 200

C ... Premature termination:

  980 WRITE (LUNCRT, 1001) ' Abnormal termination.'
      GO TO 999

C ... Normal termination:

  990 WRITE (LUNCRT, 1001) ' *** Results are in smooth2d.out. ***'

  999 CONTINUE

C ... Formats:

 1001 FORMAT (A)
 1002 FORMAT (1X, A)
 1004 FORMAT (A, I5)
 1006 FORMAT (A, F7.2)
 1008 FORMAT (2I4, ' ! ', A, 2X, A)
 1010 FORMAT (1I4, ' ! ', A, 2X, A)
 1012 FORMAT (A, 1P, 2E12.3)
 1016 FORMAT (1X, 1P, 3E16.7)

      END
