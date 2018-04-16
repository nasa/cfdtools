C+----------------------------------------------------------------------
C
      SUBROUTINE PL3D ( LUNCRT, LUNKBD, LUNDIP, DIPFIL,
     >                  NTITLE, TITLE, XLABEL, YLABEL, ZLABEL,
     >                  NCURVE, NPTS, MXPTS, X, Y, Z,
     >                  LINE, SYMBOL )
C
C  ACRONYM: PLot in 3-Dimensions
C           --      - -
C
C  PURPOSE:
C
C    PL3D is the approximate 3-d analogue of PL2D but is less general.
C
C    PL3D plots one or more curves in three dimensions on any of several
C    plotting devices using the DISSPLA graphics package.  It constructs
C    the entire plot (no preliminary set up assumed), with equal scaling
C    of the axes enforced, and is intended for interactive applications.
C    Multiple plots  may be viewed on a screen in one run, and  optional
C    hard-copy may be generated along the way.   Plot type is linear  on
C    all axes;  line type  and  symbols are selected  at a higher level,
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
C    (2)    Prompt for a plot window and a view angle.
C           Error-check step sizes selected by the user.
C        End if.
C    (3) Suppress DISSPLA messages, border, and grace margin.
C        Place the screen in graphics mode.
C        Set the page size, and define the subplot area.
C        Write axis labels and plot title(s).
C        Force equal scaling on all axes.
C        Define the 3-d plot volume.
C        Calculate the radius of the viewpoint.
C        Set the viewpoint.
C        Initialize the graph, with automatic selection of any
C        step size not specified by the user.
C        Loop over number of curves:
C           Draw curve(n).
C        Terminate the plot.
C    (4) Write to a .DIP file (replotting=true, return to 1), or
C        choose new axis limits (return to 2), or
C        clear the screen and place it in transparent mode.
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
C   XLABEL     -      C    I   Horizontal axis label (base plane).
C   YLABEL     -      C    I   Vertical    "     "      "    "   .
C   ZLABEL     -      C    I   Vertical    "     "   (normal to base).
C   NCURVE     -      I    I   No. of curves to be plotted in this frame.
C   MXPTS      -      I    I   Max. number of points in any one curve; i.e.,
C                              declared first dimension of X, Y, and Z arrays.
C   NPTS    NCURVE    I    I   Number of points in a given curve.
C   X,Y,Z NPTS,NCURVE R    I   Curve(s) to be plotted.
C   LINE              C    I   Line type for each curve (see NOTES).
C   SYMBOL            I    I   Symbol type (see NOTES).
C
C  PARAMETER CONSTANTS:
C    NAME    TYPE   DESCRIPTION
C   AUTO      R     Flag to signal automatic scaling, etc.
C   XPAGE     R     Horizontal page size (inches).
C   YPAGE     R     Vertical    "    "       "   .
C   DOLLAR    C     Dollar sign (DISSPLA end-of-string delimiter).
C
C  SIGNIFICANT LOCAL VARIABLES:
C    NAME    DIM   TYPE DESCRIPTION
C   XAXIS     -     R   Horizontal (base plane) axis length (inches).
C   YAXIS     -     R   Vertical      "    "     "      "      "    .
C   LTYPE  MXLTYP   C   Line type names.
C
C  EXTERNAL REFERENCES:
C    NAME    DESCRIPTION                        SOURCE
C   AXLIM    Prompts for axis limits.           PRMODULES (part of PL2D)
C   READER   Prompting utility.                 PRMODULES
C   SCLEAN   Screen clearing utility.           VT640
C   TERMIN   Adds a marker to end of string.    PRMODULES
C   Also:    DISSPLA graphics library.
C
C  NOTES:
C    (1) PRMODULES = RAL::PROG:[PROGTOOLS.FORTLIB]PRMODULES.OLB
C        VT640     = RAL::GRAPH:[GRAPHICS.VT640]VT640LIB.OLB
C
C    (2) IMPLICIT NONE is non-standard.
C
C    (3) Line types are:
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
C    (4) Symbol types apply only for those line types which use symbols.
C        Legal values are [0,18] with automatic wraparound if > 18.
C        Refer to a DISSPLA manual for the key.
C
C  ENVIRONMENT:  VAX/VMS, FORTRAN 77
C
C  DEVELOPMENT HISTORY:
C     DATE   INITIALS   DESCRIPTION
C   10/01/86    RGL     Initial design and code.
C   04/01/87    RGL     Generalized somewhat for library inclusion.
C                       (PL3D is not as easy to generalize as PL2D.)
C   08/31/89    MDW     Updated to run on DISSPLA version 11.0.
C   04/19/90    RGL     Will plot even when NCURVE exceeds dimensions.
C
C  AUTHOR:  Ronald G. Langhi, Sterling Software, Palo Alto, CA
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C ... Arguments:

      INTEGER
     >   LUNCRT, LUNKBD, LUNDIP, MXPTS, NCURVE, NPTS(NCURVE), NTITLE,
     >   SYMBOL(NCURVE)
      REAL
     >   X(MXPTS,NCURVE), Y(MXPTS,NCURVE), Z(MXPTS,NCURVE)
      CHARACTER
     >   DIPFIL*(*), LINE(NCURVE)*(*), TITLE(*)*(*),
     >   XLABEL*(*), YLABEL*(*), ZLABEL*(*)

C ... Constants:

      INTEGER
     >   ICPS, MXLTYP
      REAL
     >   AUTO, ONE, XPAGE, YPAGE, ZERO
      CHARACTER
     >   BLANK*1, DOLLAR*1
      PARAMETER
     >   (AUTO=999.E+0, BLANK=' ', DOLLAR='$', ICPS=960, MXLTYP=9,
     >   ONE=1.E+0, XPAGE=11.0, YPAGE=8.5, ZERO=0.E+0)

C ... Variables:

      INTEGER
     >   I, MODE, N, NCHARS, NCLOCL, NEXT
      REAL
     >   DXDATA, XAXIS, XBGN, XCHECK, XEND, XINC,
     >   DYDATA, YAXIS, YBGN, YCHECK, YEND, YINC,
     >   DZDATA, ZAXIS, ZBGN, ZCHECK, ZEND, ZINC,
     >   ANGLAT, ANGLON, DFACT, RFACT, RPLOT, RVIEW
      LOGICAL
     >   CR, EOF, FIRST, REPLOT
      CHARACTER
     >   DIPWRK*80, LTYPE(MXLTYP)*9

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
         WRITE (LUNCRT,'(''0PL3D:  NCURVE ='',I4,
     >      '' is greater than the local maximum'',I4)') NCURVE, MXCURV
         NCLOCL = MXCURV
      END IF

      XAXIS = XPAGE - 2.0
      YAXIS = YPAGE - 2.0

      IF ( DIPFIL.EQ.BLANK ) THEN
         DIPWRK = 'PL3D.DIP'
      ELSE
         DIPWRK = DIPFIL
      END IF

C ... Count the number of characters in the .DIP filename:

      DO 20, NCHARS = LEN (DIPWRK), 1, -1
         IF ( DIPWRK(NCHARS:NCHARS) .NE. BLANK ) GO TO 40
   20 CONTINUE

   40 CONTINUE

C ... Convert character line types to integer:

      DO 100, N = 1, NCLOCL
         LINTYP(N) = 0

         DO 60, I = 1, MXLTYP
            IF ( LINE(N) .EQ. LTYPE(I) ) THEN
               LINTYP(N) = I
               GO TO 80
            END IF
   60    CONTINUE

   80    CONTINUE
         LINTYP(N) = MAX (1, LINTYP(N))
  100 CONTINUE

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

C ...    Find the extreme axes values over all curves,
C ...    and initialize the plot window:

         XBGN = X(1,1)
         XEND = X(1,1)
         YBGN = Y(1,1)
         YEND = Y(1,1)
         ZBGN = Z(1,1)
         ZEND = Z(1,1)

         DO 400, N = 1, NCLOCL
            DO 300, I = 1, NPTS(N)
               XBGN = MIN (X(I,N), XBGN)
               XEND = MAX (X(I,N), XEND)
               YBGN = MIN (Y(I,N), YBGN)
               YEND = MAX (Y(I,N), YEND)
               ZBGN = MIN (Z(I,N), ZBGN)
               ZEND = MAX (Z(I,N), ZEND)
  300       CONTINUE
  400    CONTINUE

         XINC = AUTO
         YINC = AUTO
         ZINC = AUTO
      END IF

      IF ( REPLOT ) GO TO 600

C ... Prompt for a plot window:

  500 CONTINUE

      WRITE (LUNCRT,'(''0Data range:''/
     >                '' Xmin, Xmax = '', 2F11.4/
     >                '' Ymin, Ymax = '', 2F11.4/
     >                '' Zmin, Zmax = '', 2F11.4)')
     >                XBGN, XEND, YBGN, YEND, ZBGN, ZEND

      WRITE (LUNCRT,'(''0Specify a plot window (CR to'',
     >                '' retain the present one):'')')

      CALL AXLIM ( 'Enter XMIN [, XMAX[, XINC]]:  ', XBGN, XEND, XINC )
      CALL AXLIM ( 'Enter YMIN [, YMAX[, YINC]]:  ', YBGN, YEND, YINC )
      CALL AXLIM ( 'Enter ZMIN [, ZMAX[, ZINC]]:  ', ZBGN, ZEND, ZINC )

C ... Error-check step sizes selected by the user:

      IF ( XINC.NE.AUTO ) THEN
         XCHECK = (XEND - XBGN) / XINC

         IF ( XCHECK.LT.ZERO ) THEN
            WRITE (LUNCRT,'(''0X-step and X-axis point in'',
     >                      '' opposite directions.'')')
            GO TO 500
         ELSE IF ( XCHECK.LT.ONE ) THEN
            WRITE (LUNCRT,'(''0X-step exceeds X-axis length.'')')
            GO TO 500
         END IF
      END IF

      IF ( YINC.NE.AUTO ) THEN
         YCHECK = (YEND - YBGN) / YINC

         IF ( YCHECK.LT.ZERO ) THEN
            WRITE (LUNCRT,'(''0Y-step and Y-axis point in'',
     >                      '' opposite directions.'')')
            GO TO 500
         ELSE IF ( YCHECK.LT.ONE ) THEN
            WRITE (LUNCRT,'(''0Y-step exceeds Y-axis length.'')')
            GO TO 500
         END IF
      END IF

      IF ( ZINC.NE.AUTO ) THEN
         ZCHECK = (ZEND - ZBGN) / ZINC

         IF ( ZCHECK.LT.ZERO ) THEN
         WRITE (LUNCRT,'(''0Z-step and Z-axis point in'',
     >                      '' opposite directions.'')')
            GO TO 500
         ELSE IF ( ZCHECK.LT.ONE ) THEN
            WRITE (LUNCRT,'(''0Z-step exceeds Z-axis length.'')')
            GO TO 500
         END IF
      END IF

C ... Prompt for a view angle:

      ANGLAT = 45.E+0
      ANGLON = 45.E+0
      DFACT  = ZERO
      WRITE (LUNCRT,'(''0Choose a view angle by entering'')')
      CALL READR ( LUNCRT, 'Latitude (+/-180 or <CR> for default)  : ',
     >             LUNKBD, ANGLAT, CR, EOF )
      CALL READR ( LUNCRT, 'Longitude (+/-90 or <CR> for default)  : ',
     >             LUNKBD, ANGLON, CR, EOF )
      CALL READR ( LUNCRT, 'Distortion factor (range [0,1] or <CR>): ',
     >             LUNKBD, DFACT, CR, EOF )

  600 CONTINUE

C ... Suppress DISSPLA messages, plot border, and grace margin:

      CALL SETDEV ( 0, 0 )
      CALL NOBRDR
      CALL GRACE ( ZERO )

      IF ( MODE.NE.3 ) THEN
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

      CALL AREA2D ( XAXIS, YAXIS )

C ... Use integer numbering on the axes:

      CALL INTAXS

C ... Set the numbers perpendicular to each axis:

      CALL XAXANG ( 90. )
      CALL YAXANG (  0. )
      CALL ZAXANG ( 90. )

C ... Mark labels and titles as "self-counting" for DISSPLA
C ... before writing them on the plot:

      IF ( NTITLE.GT.0 ) THEN
         DO 700, N = 1, NTITLE
            CALL TERMIN ( TITLE(N), DOLLAR )
            CALL HEADIN ( TITLE(N), 100, ONE, NTITLE )
  700    CONTINUE
      END IF

      CALL TERMIN ( XLABEL, DOLLAR )
      CALL TERMIN ( YLABEL, DOLLAR )
      CALL TERMIN ( ZLABEL, DOLLAR )

C ... Force equal scaling on all axes:

      CALL AXES3D ( XLABEL, 100, YLABEL, 100,
     >              ZLABEL, 100, ZERO, ZERO, ZERO )

C ... Define the 3-d plot volume:

      CALL VOLM3D ( ZERO, ZERO, ZERO )

C ... Calculate the radius of the viewpoint:

      IF ( DFACT.LT.ZERO ) DFACT = ZERO
      IF ( DFACT.GT.ONE  ) DFACT = ONE

      IF ( DFACT.EQ.ZERO ) THEN
         RFACT = 100.E+0
      ELSE
         RFACT = (ONE - DFACT) / DFACT
      END IF

      DXDATA = XEND - XBGN
      DYDATA = YEND - YBGN
      DZDATA = ZEND - ZBGN

      RPLOT = 0.5E+0 * SQRT (DXDATA**2 + DYDATA**2 + DZDATA**2)
      RVIEW = RPLOT + RFACT * RPLOT

C ... Set the viewpoint:

      CALL VUANGL ( ANGLAT, ANGLON, RVIEW )

C ... Initialize the graph:

      IF ( XINC.EQ.AUTO ) THEN
         IF ( YINC.EQ.AUTO ) THEN
            IF ( ZINC.EQ.AUTO ) THEN
               CALL GRAF3D ( XBGN, 'SCALE', XEND,
     >                       YBGN, 'SCALE', YEND,
     >                       ZBGN, 'SCALE', ZEND )
            ELSE
               CALL GRAF3D ( XBGN, 'SCALE', XEND,
     >                       YBGN, 'SCALE', YEND,
     >                       ZBGN,  ZINC,   ZEND )
            END IF
         ELSE
            IF ( ZINC.EQ.AUTO ) THEN
               CALL GRAF3D ( XBGN, 'SCALE', XEND,
     >                       YBGN,  YINC,   YEND,
     >                       ZBGN, 'SCALE', ZEND )
            ELSE
               CALL GRAF3D ( XBGN, 'SCALE', XEND,
     >                       YBGN,  YINC,   YEND,
     >                       ZBGN,  ZINC,   ZEND )
            END IF
         END IF
      ELSE
         IF ( YINC.EQ.AUTO ) THEN
            IF ( ZINC.EQ.AUTO ) THEN
               CALL GRAF3D ( XBGN,  XINC,   XEND,
     >                       YBGN, 'SCALE', YEND,
     >                       ZBGN, 'SCALE', ZEND )
            ELSE
               CALL GRAF3D ( XBGN,  XINC,   XEND,
     >                       YBGN, 'SCALE', YEND,
     >                       ZBGN,  ZINC,   ZEND )
            END IF
         ELSE
            IF ( ZINC.EQ.AUTO ) THEN
               CALL GRAF3D ( XBGN,  XINC,   XEND,
     >                       YBGN,  YINC,   YEND,
     >                       ZBGN, 'SCALE', ZEND )
            ELSE
               CALL GRAF3D ( XBGN, XINC, XEND,
     >                       YBGN, YINC, YEND,
     >                       ZBGN, ZINC, ZEND )
            END IF
         END IF
      END IF

C ... Draw one or more curves:

      DO 800, N = 1, NCLOCL
         CALL CURV3D ( X(1,N), Y(1,N), Z(1,N), NPTS(N), 0 )
  800 CONTINUE

C ... Terminate the plot:

      CALL ENDPL ( 0 )

      FIRST = .FALSE.
      REPLOT     = .FALSE.

      NEXT = 3
      WRITE (LUNCRT,'(''01 - Write to a .DIP file''/
     >                '' 2 - Change the plot window'')')
      CALL READI ( LUNCRT, '3 - End this plot (default):  ',
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
