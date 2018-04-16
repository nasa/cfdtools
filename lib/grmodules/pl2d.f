C+----------------------------------------------------------------------
C
      SUBROUTINE PL2D ( LUNCRT, LUNKBD, LUNDIP, DIPFIL,
     >                  NTITLE, TITLE, XLABEL, YLABEL,
     >                  NCURVE, NPTS, MXPTS, X, Y,
     >                  PLOT, LINE, SYMBOL )
C
C  ACRONYM: PLot in 2-D
C           --      - -
C
C  PURPOSE:
C
C    PL2D  plots one or more curves in two dimensions  on any of several
C    plotting devices using the DISSPLA graphics package.  It constructs
C    the entire plot (no preliminary set up assumed) and is intended for
C    interactive applications.  Multiple plots may be viewed on a screen
C    in one run, and optional hard-copy may  be generated along the way.
C    Plot type,  line type,  and symbols are selected at a higher level,
C    while suitable step sizes for axes are either user-selected here or
C    chosen automatically.
C
C    The available devices are:
C
C       1 - Any graphics terminal which emulates a Tektronix 4014
C           (e.g. DEC VT640, SELANAR HiREZ 100XL, GraphOn GO-250)
C       2 - Device-independent plot (DIP) file
C
C  METHOD:
C
C    (0) Do initialization, which includes computing physical axis
C        lengths, determining the number of characters in the .DIP file
C        name, and converting character line type codes to integer.
C    (1) Open the plotting device.
C        If ( first pass ) initialize plot window based on input data.
C        If ( not replotting ) then
C    (2)    Prompt for a plot window.
C           Set up the plot axes.
C        End if.
C    (3) Suppress DISSPLA messages, border, and grace margin.
C        Place the screen in graphics mode.
C        Set the page size, and define the subplot area.
C        Write axis labels and plot title(s).
C        Initialize the graph.
C        Loop over number of curves:
C           Draw curve(n).
C        Terminate the plot.
C    (4) Write to a .DIP file (replotting=true, return to 1), or
C        choose new axis limits (return to 2), or
C        clear the screen and place it in transparent mode and return
C        to the calling routine.
C
C  ARGUMENTS:
C    NAME     DIM   TYPE I/O/S DESCRIPTION
C   LUNCRT     -      I    I   Logical unit no. for prompts & screen plots.
C   LUNKBD     -      I    I   Logical unit no. for keyboard entries.
C   LUNDIP     -      I    I   Logical unit no. for device-independent
C                              plot (DIP) files (local to Ames).
C   DIPFIL     -      C    I   DIP file name.
C   NTITLE     -      I    I   Number of title lines.
C   TITLE   NTITLE    C    I   Plot title(s).
C   XLABEL     -      C    I   Horizontal axis label.
C   YLABEL     -      C    I   Vertical    "     "  .
C   NCURVE     -      I    I   No. of curves to be plotted in this frame.
C   NPTS    NCURVE    I    I   Number of points in a given curve.
C   MXPTS      -      I    I   Max. number of points in any one curve; i.e.,
C                              declared first dimension of X and Y arrays.
C   X, Y NPTS,NCURVE  R    I   Curve(s) to be plotted.
C   PLOT       -      C    I   Plot type (see NOTES).
C   LINE    NCURVE    C    I   Line type for each curve (see NOTES).
C   SYMBOL  NCURVE    I    I   Symbol type (see NOTES).
C
C  PARAMETER CONSTANTS:
C    NAME    TYPE   DESCRIPTION
C   AUTO      R     Flag to signal automatic scaling, etc. (see NOTES).
C   XPAGE     R     Horizontal page size (inches).
C   YPAGE     R     Vertical     "    "     "    .
C   DOLLAR    C     Dollar sign (DISSPLA end-of-string delimiter).
C
C  SIGNIFICANT LOCAL VARIABLES:
C    NAME    DIM   TYPE DESCRIPTION
C   DEFWID    -     R   Default horizontal axis physical length (inches).
C   DEFHGT    -     R   Default vertical     "     "        "      "    .
C   LTYPE   MXLTYP  C   Line type names.
C   XBOUND,   2     R   Data extremes - needed for axis layout.
C   YBOUND
C
C  EXTERNAL REFERENCES:
C    NAME    DESCRIPTION                        SOURCE
C   AXLIM    Prompts for axis limits.           Appended to this routine.
C   LAYOUT   Plot axis setup.                   GRMODULES
C   READER   Prompting utility.                 PRMODULES
C   SCLEAN   Screen clearing utility.           VT640
C   TERMIN   Adds a marker to end of string.    PRMODULES
C   DISSPLA graphics library.
C
C  NOTES:
C    (1) GRMODULES = RAL::GRAPH:[GRAPHICS.GENERALIB]GRMODULES.OLB
C        PRMODULES = RAL::PROG:[PROGTOOLS.FORTLIB]PRMODULES.OLB
C        VT640     = RAL::GRAPH:[GRAPHICS.VT640]VT640LIB.OLB
C
C    (2) IMPLICIT NONE is non-standard.
C
C    (3) Plot types are:  LINEAR, LOGLOG, LOGX, LOGY, POLAR, or SCALE.
C        For a full set of synonyms see subroutine LAYOUT.
C
C    (4) Line types are:
C        DEFAULT   = symbols plus solid line
C        SYMBOLS   = symbols alone
C        SOLID     = solid line
C        DOT       = dots
C        DASH      = short dashes
C        CHAINDOT  = dash-dot
C        CHAINDASH = long-short dashes
C        THICK     = heavy stroke
C        LONGDASH  = long dashes
C
C    (5) Symbol types apply only for those line types which use symbols.
C        Legal values are [0,18], but with automatic wraparound if > 18.
C        Refer to a DISSPLA manual for the key.
C
C    (6) AUTO must match the value of FLAG in subroutine LAYOUT.
C
C    (7) Diagnostics from LAYOUT have been suppressed by passing
C        a negative-valued logical unit number.
C
C  ENVIRONMENT:  VAX/VMS, FORTRAN 77 with DISSPLA
C
C  DEVELOPMENT HISTORY:
C     DATE   INITIALS   DESCRIPTION
C   10/01/86    RGL     Initial design and code.
C   04/01/87    RGL     Generalized for library inclusion.
C   08/31/89    MDW     Updated to run on DISSPLA version 11.0.
C   04/19/90    RGL     Will plot even when NCURVE exceeds dimensions.
C
C  AUTHOR:  Ronald G. Langhi, Sterling Software, Palo Alto, CA
C
C-----------------------------------------------------------------------


C     Declarations
C     ------------

      IMPLICIT NONE

C ... Arguments:

      INTEGER
     >   LUNCRT, LUNKBD, LUNDIP, NCURVE, MXPTS, NPTS(NCURVE), NTITLE,
     >   SYMBOL(NCURVE)
      REAL
     >   X(MXPTS,NCURVE), Y(MXPTS,NCURVE)
      CHARACTER
     >   DIPFIL*(*), LINE(NCURVE)*(*), PLOT*(*), TITLE(*)*(*),
     >   XLABEL*(*), YLABEL*(*)

C ... Constants:

      INTEGER
     >   ICPS, MXLTYP
      REAL
     >   AUTO, ONE, XPAGE, YPAGE, ZERO
      CHARACTER
     >   BLANK*1, DOLLAR*1
      PARAMETER
     >   (AUTO=999.E+0, BLANK=' ', DOLLAR='$', ICPS=960,
     >   MXLTYP=9, XPAGE=11.0, YPAGE=8.5, ONE=1.E+0, ZERO=0.E+0)

C ... Variables:

      INTEGER
     >   I, MODE, N, NCHARS, NCLOCL, NEXT, TYPE
      REAL
     >   DEFHGT, DEFWID, HEIGHT, WIDTH, XBOUND(2), YBOUND(2),
     >   XBGN, XEND, XINC, YBGN, YEND, YINC
      LOGICAL
     >   CR, EOF, FIRST, REPLOT
      CHARACTER
     >   DIPWRK*80, LTYPE(MXLTYP)*9, PLTWRK*40

C ... Commons:

      INTEGER
     >   MXCURV
      PARAMETER
     >   (MXCURV=1501)

      INTEGER
     >   LINTYP(MXCURV)

      COMMON /SPECS/ LINTYP

C ... Data:

      DATA
     >   LTYPE /'DEFAULT  ',
     >          'SYMBOLS  ',
     >          'SOLID    ',
     >          'DOT      ',
     >          'DASH     ',
     >          'CHAINDOT ',
     >          'CHAINDASH',
     >          'THICK    ',
     >          'LONGDASH '/

C ... Initialization:


      NCLOCL = NCURVE

      IF ( NCURVE.GT.MXCURV ) THEN
         WRITE (LUNCRT,'(''0PL2D:  NCURVE ='',I4,
     >      '' is greater than the local maximum'',I4)') NCURVE, MXCURV
         NCLOCL = MXCURV
      END IF

      DEFHGT = YPAGE - 2.0
      DEFWID = XPAGE - 2.0

      HEIGHT = AUTO
      WIDTH  = AUTO

      PLTWRK = PLOT

      IF ( DIPFIL.EQ.BLANK ) THEN
         DIPWRK = 'PL2D.DIP'
      ELSE
         DIPWRK = DIPFIL
      END IF

C ... Count the number of characters in the .DIP file name:

      DO 60, NCHARS = LEN (DIPWRK), 1, -1
         IF ( DIPWRK(NCHARS:NCHARS) .NE. BLANK ) GO TO 80
   60 CONTINUE

   80 CONTINUE

C ... Convert character line types to integer:

      DO 140, N = 1, NCLOCL
         LINTYP(N) = 0

         DO 100, I = 1, MXLTYP
            IF ( LINE(N) .EQ. LTYPE(I) ) THEN
               LINTYP(N) = I
               GO TO 120
            END IF
  100    CONTINUE

  120    CONTINUE
         LINTYP(N) = MAX (1, LINTYP(N))
  140 CONTINUE

      FIRST = .TRUE.
      REPLOT = .FALSE.
      MODE = 1

  200 CONTINUE

C ... Open the plotting device:

      IF ( MODE.EQ.1 ) THEN

C ...    Tektronix 4014 emulation:

         CALL TK4014 ( ICPS, 1 )

      ELSE IF ( MODE.EQ.2 ) THEN

C ...    DIP file:

         CALL DIP ( LUNDIP, DIPWRK, NCHARS )

      END IF

      IF ( FIRST ) THEN

C ...    Find the extreme axes values over all curves:

         XBGN = X(1,1)
         XEND = X(1,1)
         YBGN = Y(1,1)
         YEND = Y(1,1)

         DO 400, N = 1, NCLOCL
            DO 300, I = 1, NPTS(N)
               XBGN = MIN (X(I,N), XBGN)
               XEND = MAX (X(I,N), XEND)
               YBGN = MIN (Y(I,N), YBGN)
               YEND = MAX (Y(I,N), YEND)
  300       CONTINUE
  400    CONTINUE

         XINC = AUTO
         YINC = AUTO
      END IF

      IF ( .NOT.REPLOT ) THEN

C ...    Prompt for a plot window:

         WRITE (LUNCRT,'(''0Data range:''/
     >                   '' Xmin, Xmax = '', 2F11.4/
     >                   '' Ymin, Ymax = '', 2F11.4)')
     >                   XBGN, XEND, YBGN, YEND

  500    WRITE (LUNCRT,'(''0Specify a plot window (CR to'',
     >                   '' retain the present one):'')')

         CALL AXLIM ('Enter XMIN[, XMAX[, XINC]]:  ', XBGN, XEND, XINC )
         CALL AXLIM ('Enter YMIN[, YMAX[, YINC]]:  ', YBGN, YEND, YINC )

         XBOUND(1) = XBGN
         XBOUND(2) = XEND
         YBOUND(1) = YBGN
         YBOUND(2) = YEND

C ...    Set up the plot axes:

         CALL LAYOUT ( -LUNCRT, 2, XBOUND, YBOUND,
     >                 XBGN, XEND, XINC, YBGN, YEND, YINC,
     >                 HEIGHT, WIDTH, DEFHGT, DEFWID, PLTWRK )

      END IF

C ... Suppress DISSPLA messages, plot border, and grace margin:

      CALL SETDEV ( 0, 0 )
      CALL NOBRDR
      CALL GRACE ( ZERO )

      IF ( MODE.NE.2 ) THEN

C ...    Clear the screen and place it into graphics mode:

         CALL SCLEAN ( 2 )
         CALL HWSCAL ( 'SCREEN' )
      ELSE
         CALL HWSCAL ( 'NONE' )
      END IF

      CALL HWROT ( 'AUTO' )

C ... Set the page size:

      CALL PAGE ( XPAGE, YPAGE )

C ... Define the subplot area:

      CALL AREA2D ( WIDTH, HEIGHT )

C ... Mark labels and titles as "self-counting" for DISSPLA
C ... before writing them on the plot:

      CALL TERMIN ( XLABEL, DOLLAR )
      CALL TERMIN ( YLABEL, DOLLAR )

      CALL XNAME ( XLABEL, 100 )
      CALL YNAME ( YLABEL, 100 )

      IF ( NTITLE.GT.0 ) THEN
         DO 600, N = 1, NTITLE
            CALL TERMIN ( TITLE(N), DOLLAR )
            CALL HEADIN ( TITLE(N), 100, ONE, NTITLE )
  600    CONTINUE
      END IF

C ... Initialize the graph:

      CALL GRAF ( XBGN, XINC, XEND, YBGN, YINC, YEND )

C ... Draw one or more curves:

      DO 700, N = 1, NCLOCL
         IF ( LINTYP(N) .LE. 2 ) THEN

C           Curve with symbols - call MARKER to set the symbol.
C           ---------------------------------------------------

            CALL MARKER ( SYMBOL(N) )

            CALL CURVE ( X(1,N), Y(1,N), NPTS(N), 3 - 2 * LINTYP(N) )
         ELSE

C           Curve without symbols - DISSPLA calls MYSPEC for line type.
C           -----------------------------------------------------------

            CALL MYSPEC ( LINTYP(N) )

            CALL CURVE ( X(1,N), Y(1,N), NPTS(N), 0 )

C ...       Reset DISSPLA line type parameters.
C ...       (All broken-line patterns are reset by 'DOT'.)

            CALL RESET ('DOT')
            CALL RESET ('THKCRV')
         END IF
  700 CONTINUE

C ... Terminate the plot:

      CALL ENDPL ( 0 )

      FIRST  = .FALSE.
      REPLOT = .FALSE.

      NEXT = 3
      WRITE (LUNCRT,'(''0 1 - Write to a .DIP file''/
     >                ''  2 - Change the plot window'')')
      CALL READI ( LUNCRT, ' 3 - End this plot (default):  ',
     >             LUNKBD, NEXT, CR, EOF )

      IF ( NEXT.EQ.1 ) THEN

C ...    Write the current screen plot to a .DIP file:

         REPLOT = .TRUE.
         MODE = 2
         GO TO 200

      ELSE IF ( NEXT.EQ.2 ) THEN

C ...    Choose new axis limits.

         IF ( MODE.EQ.2 ) THEN

C ...       Change from .DIP file to screen display:

            MODE = 1
            GO TO 200
         ELSE

C ...       Already in screen mode:

            GO TO 500
         END IF

      ELSE

C ...    Clear the screen and return to transparent mode:

         CALL SCLEAN ( 3 )
      END IF

  999 RETURN

      END
C+----------------------------------------------------------------------
C
      SUBROUTINE AXLIM ( PROMPT, HEAD, TAIL, STEP )
C
C  ACRONYM: AXis LIMits
C
C  PURPOSE: AXLIM prompts for and reads one set of axis limits,
C           allowing fewer than the full set to be input.
C
C  METHOD:  1) Prompt for begin[,end[,increment]].
C           2) Read response as a character string.
C           3) Determine the number of tokens in the string.
C           4) Using internal reads, assign begin, end, and STEPement
C              in accordance with the number of tokens found.
C
C  ARGUMENTS:
C    NAME    DIM    TYPE I/O/S  DESCRIPTION
C   PROMPT    -      C     I    Message for screen prompt.
C   HEAD      -      R     O    Beginning of axis sub-range.
C   TAIL      -      R     O    End       "   "    "    "  .
C   STEP      -      R     O    Axis increment.
C
C  PARAMETER CONSTANTS:
C    NAME    TYPE   DESCRIPTION
C   LUNCRT    I     Logical unit for screen output.
C   LUNKBD    I        "     "    "  keyboard input.
C   MXNUM     I     Maximum number of tokens in string.
C
C  EXTERNAL REFERENCES:
C    NAME    DESCRIPTION AND SOURCE
C   READER   Prompting utility.
C   TOKENS   Parses a character string.
C
C  ERROR HANDLING:  An axis of zero length is trapped, and the
C                   prompt for axis limits is re-issued.
C
C  SYSTEM DEPENDENCIES:
C   (1) IMPLICIT NONE is a VAX extension.
C
C  ENVIRONMENT:  VAX/VMS, FORTRAN 77
C
C  DEVELOPMENT HISTORY:
C     DATE   INITIALS   DESCRIPTION 
C   06/10/86    RGL     Initial design and code.
C   10/01/86    RGL     Adapted from RDSET.
C
C  AUTHOR:  Ronald G. Langhi, Sterling Software, Palo Alto, CA
C
C-----------------------------------------------------------------------


C     Declarations.
C     ------------

      IMPLICIT NONE

C     Arguments:

      REAL
     >   HEAD, TAIL, STEP
      CHARACTER
     >   PROMPT*(*)

C     Constants:

      INTEGER
     >   LUNCRT, LUNKBD, MXNUM
      CHARACTER
     >   BLANK*1
      PARAMETER
     >   (LUNCRT=6, LUNKBD=5, MXNUM=3, BLANK=' ')

C     Variables:

      INTEGER
     >   NUMBER
      CHARACTER
     >   LIMITS*30, VALUE(MXNUM)*10
      LOGICAL
     >   CR, EOF


C     Begin execution.
C     ---------------

  100 CONTINUE

C        Prompt for one set of axis limits:

         LIMITS = BLANK

         CALL READS ( LUNCRT, PROMPT, LUNKBD, LIMITS, CR, EOF )

         IF ( CR .OR. EOF ) GO TO 999

         NUMBER = 3

         CALL TOKENS ( LIMITS, NUMBER, VALUE )

         IF ( NUMBER.GT.0 ) READ (VALUE(1),*,ERR=800) HEAD
         IF ( NUMBER.GT.1 ) READ (VALUE(2),*,ERR=800) TAIL
         IF ( NUMBER.GT.2 ) READ (VALUE(3),*,ERR=800) STEP

C        Error checking:

         IF ( HEAD.EQ.TAIL ) THEN
            WRITE (LUNCRT,'('' This leaves no axis at all. '',
     >                      '' Try again.'')')
            GO TO 100
         END IF

         IF ( NUMBER.EQ.0 ) GO TO 999
      GO TO 100

  800 WRITE (LUNCRT,'(''0Input error!  Try again.''/)')
      GO TO 100

  999 RETURN

      END
