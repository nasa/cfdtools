C+------------------------------------------------------------------------------
C
      PROGRAM INTEGRATE
C
C     PURPOSE:
C
C        INTEGRATE drives a choice of numerical integration routines.
C     Originally, selected methods were applied to one or more sets of data
C     read from one file (stored by rows or by columns in free format).
C     Now, one dataset is treated, and it is read using the much newer
C     table_io module at the expense of some mostly-unused features.
C
C        Calculation of moment about a user-supplied abscissa is included
C     as an option.  Results are written to the screen and to a printable
C     file.
C
C     ORIGINAL DATA FORMAT:
C
C        Rows:
C           Title <up to 80 characters>
C           NX
C           X(1)  X(2)  X(3)  ...  X(?)     <May exceed 1 line>
C           ...  X(NX)
C           Y(1)  Y(2)  Y(3)  ...  Y(?)     <Starts on new line>
C           ...  Y(NX)
C           <Repeat from NX for more data sets>
C
C        Columns:
C           Title <up to 80 characters>
C           NX
C           X(1)  Y(1)
C           X(2)  Y(2)
C           ...   ...
C           X(NX)  Y(NX)
C           <Repeat from NX for more data sets>
C
C     CURRENT DATA FORMAT:
C
C           Any set of columns with optional header lines that are ignored
C           and identified as not purely numeric.  Column 1 is taken to be
C           the abscissa, while any other column may be specified as the
C           ordinate.
C
C     METHOD:
C
C        INITIALIZATION:
C           Open the output file.
C           Prompt for the name of the input data file, and read/store it.
C           Prompt for the number of the column to be integrated against
C           column 1.  Some methods assume full abscissa range; some don't.
C
C        SELECTION:
C           Select an integration routine (else quit).
C           Integrate the (X,Y) distribution with the chosen method.
C           Output evaluated area to both screen and formatted file.
C           GO TO SELECTION.
C
C     PARAMETER CONSTANTS:
C
C        LENMENU    I    Maximum length of menu character strings.
C        MXMENU     I       "     "  "  integration options (including help).
C        MXPTS      I       "     "  "  data points per curve.
C
C     SIGNIFICANT LOCAL VARIABLES:
C
C        NAME       DIM     TYPE DESCRIPTION
C        COEFS    MXPTS,3    R   Spline coefficients.
C        X         MXPTS     R   One set of data abscissas (monotonic, but
C                                not necessarily increasing).  (Routine AREAXY
C                                is an exception - it expects a polygon.)
C        Y         MXPTS     R   Corresponding ordinates.
C        YMOM      MXPTS     R   Derived ordinates for moment calculation.
C
C     FILES USED:
C
C        LUNCRT     O     Terminal prompts and error messages.
C        LUNDAT     I     Data to be integrated (one or more datasets).
C        LUNKBD     I     Keyboard entries.
C        LUNPRT     O     Printed output, including debug.
C
C     EXTERNAL REFERENCES:
C
C        AREAXY   INTEGRATE  Calculates area of irregular simple polygon.
C        AVINT    INTEGRATE  Integrates by averaging quadratic functions.
C        COUNTR    FORTLIB   Counts the number of values on a line.
C        CSFIT    INTERPLIB  Fits conventional cubic spline.
C        CSQUAD   INTEGRATE  Cubic spline quadrature (rather specialized).
C        FSCOEF   INTERPLIB  Generates coefficients for FSQUAD.
C        FSQUAD   INTEGRATE  Integrates given finite Fourier series on [X1,X2].
C        LCSAREAS INTEGRATE  Variant of LCSQUAD for saving sub-areas.
C        LCSQUAD  INTEGRATE  Storage-efficient local cubic spline quadrature.
C        READER    FORTLIB   Prompting utility.
C        RVERSE   NUMODULES  Reverses the order of array elements.
C        SCTRPZ   INTEGRATE  Integrates by composite trapezoidal rule.
C        XFORMX   NUMODULES  Transforms [a,b] to [p,q].
C
C     ENVIRONMENT:  Fortran 90
C
C     HISTORY:
C        RGL   11/19/84   Initial design/code adapted from program SMOOTH.
C        RGL   11/29/84   Added calculation of moment about chosen abscissa.
C        RGL   05/21/85   Added a title for conformity to standard format.
C        DAS   09/16/86   CSQUAD in place of DCSQDU for use with CSFIT (mainly
C                         to test CSQUAD, extracted from FLO6QNM).  DCSQDU is
C                         retained with IQHSCU.  Bug in calls to DCSQDU fixed.
C                         Wrong values of IEND in calls to CSFIT fixed.
C        DAS   03/04/87   Added AREAXY here (rather than a separate driver).
C        DAS   03/13/87   Added FSCOEF + FSQUAD, to test FSQUAD.  (Pointless for
C                         full-range integration - do we need anything else?)
C        DAS   08/24/87   Installed MSFIT/CSQUAD option.
C        DAS   08/25/87   Wrong!  MSQUAD was needed - CSQUAD is strictly for
C                         conventional CSFIT-type splines.
C        DAS   06/10/92   Substituted LCSFIT for MSFIT/MSQUAD.
C        DAS   05/29/02   Eliminated IMSL references; some F90 upgrades;
C                         needed 64-bit version to test C versions of LCSQUAD.
C        DAS   06/13/02   Added LCSAREAS in order to test it; provided choice
C                         of spline method for LCSQUAD as for LCSAREAS.
C        DAS   01/10/16   Dinesh Prabhu needed tens of thousands of data points.
C        DAS   11/13/16   Made use of the table_io module.  Changes are kept
C                         to a minimum by copying data to the original arrays.
C        DAS   05/17/18   A file with nearly 2 million wavelengths required more
C                         digits for printing NROWS.  Why does LCSQUAD give
C                         NaN for the total integral?
C
C     AUTHOR:  Ronald Langhi, Informatics General Corp., Palo Alto, CA
C              David Saunders, AMA, Inc. at NASA Ames Research Center.
C
C-------------------------------------------------------------------------------

C     MODULES:

      USE TABLE_TYPE_MODULE  ! Both are in table_io.f90
      USE TABLE_IO

      IMPLICIT NONE

C     Constants:

      INTEGER, PARAMETER ::
     +   LENMENU = 60,    ! Length of menu strings
     +   LUNCRT  = 6,     ! Screen info and prompts
     +   LUNDAT  = 1,     ! Input data table
     +   LUNKBD  = 5,     ! Keyboard inputs
     +   LUNPRT  = 7,     ! Printable output
     +   MXMENU  = 7      ! Number of method choices

      REAL, PARAMETER ::
     +   ONE = 1., ZERO = 0.

C     Variables:

      INTEGER   :: DATASET, I, IEND, IER, IORD, ITYPE, LUN,
     +             NTOKS, N, NX

      REAL      :: AREA, CENTER, MOMENT, RANGE(0:3), TWOPI, XA, XB

CCC   REAL      :: A(0:MXPTS/2), B(0:MXPTS/2), AREAS(MXPTS),
CCC  +          :: COEFS(MXPTS,3),
CCC  +          :: X(MXPTS+1), XNORM(MXPTS), Y(MXPTS+1), YMOM(MXPTS)

      REAL, ALLOCATABLE, DIMENSION (:) ::
     +             A, B, AREAS, X, XNORM, Y, YNORM, YMOM

      REAL, ALLOCATABLE, DIMENSION (:,:) ::
     +             COEFS

      LOGICAL   :: BYCOLS, DEFAULT, EOF, FIRSTDATA, FOURIER, INCREASE,
     +             OK, POLYGON, QUIT, SPLINE_TYPE, XRANGE, YESMOM

      CHARACTER :: DATAFILE*64, LINE*80, MENU(0:MXMENU)*(LENMENU),
     +             METHOD*1

      TYPE (TABLE_TYPE) :: TABLE

C     Procedures:

      REAL      :: AVINT

      EXTERNAL  :: AREAXY, AVINT, COUNTR, CSFIT, CSQUAD,
     +             FSCOEF, FSQUAD, LCSAREAS, LCSQUAD, READS, READI,
     +             READR, RVERSE, SCTRPZ, XFORMX

C     Storage:

      DATA         MENU
     +  /'  0: Help',
     +   '  1: Integration by composite trapezoidal rule (SCTRPZ)',
     +   '  2: Integration of conventional spline (CSFIT/CSQUAD)',
     +   '  3: Integration of local cubic spline (LCSQUAD, any [a, b])',
     +   '  4: Integration of local cubic spline (LCSAREAS, sub-areas)',
     +   '  5: Integration by averaging of quadratic functions (AVINT)',
     +   '  6: Area of irregular simple polygon (AREAXY)',
     +   '  7: Fourier series method (FSCOEF/FSQUAD)'/


C     Execution:

C  *  Prompt for, read and store the input data:

   70 CONTINUE
         DATAFILE = 'integrate.dat'

         CALL READS (LUNCRT,
     +      'Enter data file name.  <CR> defaults to integrate.dat: ',
     +      LUNKBD, DATAFILE, DEFAULT, QUIT)
         IF (QUIT) GO TO 900

         TABLE%FILENAME = DATAFILE
         CALL TABLE_IO_READ_REAL (LUNDAT, TABLE, IER)

         IF (IER /= 0) THEN
            WRITE (LUNCRT, 1010)
     +      ' Unable to open this file.  Check spelling and try again.'
            GO TO 70
         END IF

      WRITE (LUNCRT, '(3X, A, I9)')
     +   'NHEADER:', TABLE%NHEADER,
     +   '# ROWS: ', TABLE%NROWS,
     +   '# COLS: ', TABLE%NCOLS

      DATASET = 1  ! Leave a hook for restoring multiple datasets
      IORD = 2
      IF (TABLE%ncols > 2) THEN
         CALL READI (LUNCRT,
     +   'Abscissa column is 1; what is the ordinate column? [CR = 2] ',
     +   LUNKBD, IORD, DEFAULT, QUIT)
      END IF

C  *  Prompt for abscissa location about which to calculate moment:

      CALL READR (LUNCRT,
     +   'Enter moment center abscissa; <CR> suppresses moments: ',
     +   LUNKBD, CENTER, DEFAULT, QUIT)
      YESMOM = .NOT. DEFAULT

C  *  Open the printable output file:

      OPEN (UNIT=LUNPRT, FILE='integrate.out', STATUS='UNKNOWN')

      WRITE (LUNPRT, 1010) ' Results from Program INTEGRATE'
      WRITE (LUNPRT, 1020) ' Data file:  ', TRIM (DATAFILE)
      IF (TABLE%NHEADER > 0) THEN
         WRITE (LUNPRT, 1010) ' File header:'
         DO LUN = LUNCRT, LUNPRT, LUNPRT - LUNCRT
            WRITE (LUN, '(1X, A)')
     +            (TRIM (TABLE%HEADER(I)), I = 1, TABLE%NHEADER)
         END DO
      ELSE
         WRITE (LUNPRT, '(1X, A)') TRIM (DATAFILE)
      END IF

CCC   REAL      :: A(0:MXPTS/2), B(0:MXPTS/2), AREAS(MXPTS),
CCC  +          :: COEFS(MXPTS,3), RANGE(0:3),
CCC  +          :: X(MXPTS+1), XNORM(MXPTS), Y(MXPTS+1), YMOM(MXPTS)

      NX = TABLE%NROWS
      ALLOCATE (A(0:NX/2), B(0:NX/2), AREAS(NX), COEFS(NX,3),
     +          X(NX+1), XNORM(NX), Y(NX+1), YMOM(NX))
      X(1:NX) = TABLE%VALUES(1,:)
      Y(1:NX) = TABLE%VALUES(IORD,:)
      FIRSTDATA = .TRUE.

  200 CONTINUE

C  *     Prompt for choice of integration scheme:

         WRITE (LUNCRT, 1010)
         CALL READI (LUNCRT,
     +      'Select integration method.  <CR> = done; 0 = help: ',
     +      LUNKBD, ITYPE, DEFAULT, QUIT)
         IF (DEFAULT .OR. QUIT) GO TO 900

         IF (ITYPE < 1 .OR. ITYPE > MXMENU) THEN
            WRITE (LUNCRT, 1010) MENU
            GO TO 200
         END IF

         INCREASE    = ITYPE == 0 ! Originally an IMSL method
         POLYGON     = ITYPE == 6
         FOURIER     = ITYPE == 7
         SPLINE_TYPE = ITYPE == 3 .OR. ITYPE == 4
         XRANGE      = ITYPE == 3 .OR. ITYPE == 5
         WRITE (LUNPRT, 1020) ' Method: ', MENU(ITYPE)(5:LENMENU)
         WRITE (LUNCRT, 1020) ' Method: ', MENU(ITYPE)(5:LENMENU)

C  *        Some routines require monotonically increasing abscissas:

            IF (.NOT. POLYGON) THEN

               IF (INCREASE) THEN
                  IF (X(1) > X(NX)) THEN
                     CALL RVERSE (NX, X, X)
                     CALL RVERSE (NX, Y, Y)
                  END IF
               END IF

C  *           Compute ordinates for moment calculation?

               IF (YESMOM .AND. .NOT. FOURIER) THEN
                  DO I = 1, NX
                     YMOM(I) = (CENTER - X(I)) * Y(I)
                  END DO
               END IF

C              Allow an arbitrary range of integration:

               IF (XRANGE) THEN
                  XA = X (1)
                  XB = X (NX)
                  WRITE (LUNCRT, '(/, A, 1P, 2E14.6)')
     +               ' Data range: ', XA, XB
                  CALL READR (LUNCRT,
     +               'Lower integration limit? <CR> = X(1):  ',
     +               LUNKBD, XA, DEFAULT, QUIT)
                  CALL READR (LUNCRT,
     +               'Upper integration limit? <CR> = X(NX): ',
     +               LUNKBD, XB, DEFAULT, QUIT)
               END IF
            END IF

            IF (SPLINE_TYPE) THEN
               METHOD = 'M'
               CALL READC (LUNCRT,
     +           'Type? (M=monotonic, B=loose fit, L=linear; <CR>=M): ',
     +            LUNKBD, METHOD, DEFAULT, QUIT)
            END IF

C  *        Evaluate integral(s) by the selected method -- may require
C           additional prompts, suppressed after the first dataset:

            GO TO (410, 420, 430, 440, 450, 460, 470) ITYPE

  410       CONTINUE ! *** Composite Trapezoidal Rule, SCTRPZ ***

               CALL SCTRPZ (NX, X, Y, AREA)
               IF (YESMOM) CALL SCTRPZ (NX, X, YMOM, MOMENT)

               GO TO 600

  420       CONTINUE ! *** Cubic Spline Quadrature, CSQUAD with CSFIT ***

               IF (FIRSTDATA) THEN

C  *              Prompt for the spline end conditions:

                  WRITE (LUNCRT, 1015)
     +            ' You must choose endpoint conditions as follows:',
     +            '  1 = match cubic defined by first/last 4 points;',
     +            '  2 = natural spline end conditions;',
     +            '  3 = periodic end conditions.',
     +            ' '

  425             CALL READI (LUNCRT, 'Select 1, 2, or 3: ',
     +                        LUNKBD, IEND, DEFAULT, QUIT)
                  IF (DEFAULT .OR. IEND < 1 .OR. IEND > 3) GO TO 425

                  IF (IEND == 1) THEN
                     LINE = ' Matching cubic at ends'
                     IEND = 0
                  ELSE IF (IEND == 2) THEN
                     LINE = ' Natural spline end conditions'
                  ELSE IF (IEND == 3) THEN
                     LINE = ' Periodic end conditions'
                     IEND = 4
                  END IF

                  WRITE (LUNCRT, 1010) LINE
                  WRITE (LUNPRT, 1010) LINE
               END IF

               CALL CSFIT (NX, X, Y, IEND, ZERO, IEND, ZERO,
     +                     COEFS(1,1), COEFS(1,2), COEFS(1,3), IER)
               IF (IER /= 0) GO TO 820

               CALL CSQUAD (NX, X, Y, ZERO, COEFS(1,2), COEFS(1,2))
               AREA = COEFS(NX,2)

               IF (YESMOM) THEN

                  CALL CSFIT (NX, X, YMOM, IEND, ZERO, IEND, ZERO,
     +                        COEFS(1,1), COEFS(1,2), COEFS(1,3), IER)
                  IF (IER /= 0) GO TO 820

                  CALL CSQUAD (NX, X, YMOM, ZERO, COEFS(1,2),
     +                         COEFS(1,2))
                  MOMENT = COEFS(NX,2)

               END IF

               GO TO 600

  430       CONTINUE ! *** Local Cubic Spline Quadrature (LCSQUAD) ***

               CALL LCSQUAD (NX, X, Y, XA, XB, METHOD, AREA)

               IF (YESMOM)
     +            CALL LCSQUAD (NX, X, YMOM, XA, XB, METHOD, MOMENT)

               GO TO 600

  440       CONTINUE ! *** Sub-areas variant of LCSQUAD (LCSAREAS) ***

               CALL LCSAREAS (NX, X, Y, METHOD, ZERO, AREAS)
               AREA = AREAS(NX)

               IF (NX < 200) THEN
                  WRITE (LUNCRT, '(3ES16.7)')
     +               (X(I), Y(I), AREAS(I), I = 1, NX)
               ELSE
                  WRITE (LUNCRT, '(3ES16.7)')
     +               (X(I), Y(I), AREAS(I), I = 1, 10)
                  WRITE (Luncrt, '(A)') ':::::::::::::::::::::::::::::'
                  WRITE (LUNCRT, '(3ES16.7)')
     +               (X(I), Y(I), AREAS(I), I = NX-10, NX)
               END IF

               IF (YESMOM) THEN
                  CALL LCSAREAS (NX, X, YMOM, METHOD, ZERO, AREAS)
                  MOMENT = AREAS(NX)
               END IF

               GO TO 600

  450       CONTINUE ! *** Averaging of Quadratic Functions, AVINT ***

               AREA = AVINT (X, Y, NX, XA, XB)
               IF (YESMOM) MOMENT = AVINT (X, YMOM, NX, XA, XB)
               GO TO 600

  460       CONTINUE ! *** Area of irregular, simple polygon, AREAXY ***

               CALL AREAXY (NX, X, Y, AREA)
               GO TO 600

  470       CONTINUE ! *** Fourier series method (FSCOEF/FSQUAD) ***

               N = (NX - 1) / 2
               CALL READI (LUNCRT,
     +            'Enter N to use Fourier coefficients (0:N): ',
     +            LUNKBD, N, DEFAULT, QUIT)
               IF (N < 0 .OR. N > (NX -1) / 2) GO TO 460

               DO I = LUNCRT, LUNPRT, LUNPRT - LUNCRT
                  WRITE (I, 1040)
     +               ' Fourier coefs. (0:N) used with N = ', N
               END DO

               WRITE (LUNCRT, 1010, ADVANCE='NO')
     +            ' Enter range of integration [X1:X2] (data units): '
               READ (LUNKBD, *) RANGE(1), RANGE(2)
               IF (RANGE(1) == X(1) .AND. RANGE(2) == X(NX))
     +            WRITE (LUNCRT, 1020)
     +            ' Warning: Full range integration requires just a(0)',
     +            ' - no sums needed.'

C  *           Generate a finite Fourier series for the given data:

               DO I = 0, N
                  CALL FSCOEF (3, NX-1, Y, I, A(I), B(I), IER)
               END DO

C  *           Transform the abscissas to [0, 2*pi] as FSCOEF assumes:

               TWOPI = 4.E+0 * ASIN (ONE)

               CALL XFORMX (NX, X, ZERO, TWOPI, XNORM)

C  *           XFORMX uses first and last values to calc. transformation, so:

               RANGE(0) = X(1)
               RANGE(3) = X(NX)
               CALL XFORMX (4, RANGE, ZERO, TWOPI, RANGE)

C  *           Integrate the transformed data (COEFS for scratch):

               CALL FSQUAD (3, N, A, B, RANGE(1), RANGE(2), COEFS,
     +                      AREA)
               AREA = AREA * (X(NX) - X(1)) / TWOPI
               GO TO 600

  600       CONTINUE

C  *        Output results of integration:

            IF (YESMOM .AND. .NOT. POLYGON .AND. .NOT. FOURIER) THEN
               IF (DATASET == 1) THEN
                  DO LUN = LUNCRT, LUNPRT, LUNPRT - LUNCRT
                     WRITE (LUN, 1020)
     +                  ' Dataset  # Pts.               Integral',
     +                  '   Moment Center          Moment'
                  END DO
               END IF
               DO LUN = LUNCRT, LUNPRT, LUNPRT - LUNCRT
                  WRITE (LUN, 1030) DATASET, NX, AREA, CENTER, MOMENT
               END DO
            ELSE
               IF (DATASET == 1) THEN
                  DO LUN = LUNCRT, LUNPRT, LUNPRT - LUNCRT
                     WRITE (LUN, 1010)
     +                  ' Dataset  # Pts.               Integral'
                  END DO
               END IF
               DO LUN = LUNCRT, LUNPRT, LUNPRT - LUNCRT
                  WRITE (LUN, 1030) DATASET, NX, AREA
               END DO
            END IF

C  *        Read NX for next case, if any:

CCC         READ (LUNDAT, *, END=650) NX
CCC         FIRSTDATA = .FALSE.
CCC         CYCLE

  650       EOF = .TRUE.

CCC      END DO ! Next dataset

C  *  Go back for another integration option:

      GO TO 200

C  *  Error handling:

  800 WRITE (LUNCRT, 1010) ' Error reading input data file.'
      STOP
  805 WRITE (LUNCRT, 1010) ' Data file contains only one column.'
      STOP
  810 WRITE (LUNCRT, 1040) ' Too many data points found. NX = ', NX
      STOP
  820 WRITE (LUNCRT, 1010) ' Error in fitting or integrating spline.'
      STOP

  900 CONTINUE

      WRITE (LUNCRT, 1015) ' ***  integrate.out  contains a  ***',
     +                     ' ***  copy of the screen output  ***'

 1010 FORMAT (A)
 1015 FORMAT (/, (A))
 1020 FORMAT (/, 2A)
 1030 FORMAT (I5, I9, ES25.15, 2ES16.7)
 1040 FORMAT (/, A, I4)

      END PROGRAM INTEGRATE
