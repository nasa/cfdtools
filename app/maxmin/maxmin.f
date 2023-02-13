C+----------------------------------------------------------------------
C
      PROGRAM MAXMIN
C
C PURPOSE:
C
C     MAXMIN drives a selection of methods of estimating derivative
C     information for one or more 1-D datasets ("curves").   It was
C     originally intended to test an algorithm for locating all  of
C     the maxima and minima contained in experimental measurements.
C     This searching for turning points is now optional  -  looking
C     at derivatives has become the main purpose, either to examine
C     algorithm behavior or to use the values elsewhere.
C
C     Results are presented in QPLOTable form.
C
C     This version concerns itself with 1st and 2nd derivatives.
C
C METHOD:  
C
C     MAXMIN is basically SMOOTH (which fits curves to data points)
C     adapted to look at 1st and 2nd derivatives, with an option to
C     calculate the locations of turning points by  locating  zeros
C     in the 1st derivatives.
C
C     Dataset format:
C
C        TITLE        <Default for plot frame titles>
C        N            <No. of pts. for first curve>
C        X    F       <Coordinates for first curve;
C        X    F        X ought to be monotonically increasing,
C        X    F        even though some methods don't require this;
C        :    :        distinctions could be made if necessary.>
C        :    :
C        N            <No. of pts. for next curve>
C        X    F
C        X    F
C        :    :
C        :    :       <and so on>
C
C
C     Plotting considerations:
C
C        Plotting function values, derivatives, and locations of turning
C        points is problematic because of possible scaling conflicts.
C        Any of these may be suppressed in the saved results.  (The data
C        points and fitted function values are linked - either both are
C        shown or neither is shown.)  All quantities show up on the same
C        frame (with a new frame for each method, and results for all
C        curves in the dataset on the same frame).
C
C     Method outline:
C
C        Prompt for dataset and plot title.
C        Prompt for suppressing the turning-point estimations (meaning
C           look at derivatives only).
C        Prompt for choice of smoothing routine.
C        WHILE <valid method selected> DO
C           Prompt for such things as smoothing parameters.
C           Initialize QPLOT file.
C           Read no. of data points in first dataset (curve).
C           WHILE <not EOF> DO
C              Read a dataset.
C              Fit requested smoothing curve.
C              Generate uniformly-spaced abscissas in data range.
C              Evaluate curve at each of these.
C              Also evaluate first derivative of curve at each of these.
C              Write QPLOT info for up to five curves (data, smoothed curve,
C                 first derivative curve, plus, if requested, any maxima
C                 as one curve, and any minima as another).
C              Read no. of points for next dataset.
C           END WHILE
C           Rewind data file.
C           Prompt for another method.
C        END WHILE
C
C PARAMETER CONSTANTS:
C    PARAM   TYPE   DESCRIPTION
C   MXEVAL     I    Number of points to evaluate smoothed curves at (also
C                   used for discrete search for maxima and minima)
C   MXFCOF     I    Limit on size of discrete Fourier series defined by N
C   MXMENU     I    Limit on number of smoothing options
C   MXNDEG     I    Limit on degree of polynomial fit
C   MXPTS      I    Limit on number of data points per curve
C   MXTPTS     I    Limit on number of turning points expected per curve
C
C SIGNIFICANT LOCAL VARIABLES:
C    VAR       DIM     TYPE DESCRIPTION
C     X       MXPTS     R   One set of data abscissas
C     F       MXPTS     R   One set of data ordinates
C    XNORM    MXPTS     R   Space for normalized abscissas (Fourier methods)
C    FNORM    MXPTS     R   Space for transformed ordinates (FSARBU/FSARBN)
C                           or smoothed function values (LSFIT1)
C     XK      MXPTS     R   Knot locations for ICSFKU, ICSVKU
C    WGHTS    MXPTS     R   Weights needed by spline routines
C    FSM      MXPTS     R   Smoothed result from ICSSCU, ICSFKU, ICSVKU
C     A       0:MXFCOF  R   Fourier series coefficients
C     B       1:MXFCOF  R   Fourier series coefficients
C    COEFS    MXPTS,3   R   Spline (or simple polynomial) coefficients
C    WORK     MXWORK    R   Other work-space needed by ICSSCU, etc.
C    XEVAL    MXEVAL    R   Abscissas for evaluating smoothed curves at
C    FEVAL    MXEVAL    R   Values of smoothed curves at XEVAL(*)
C    FPEVAL   MXEVAL    R   Corresp. values of 1st derivative
C    FPPEVAL  MXEVAL    R   Corresp. values of 2nd derivative (scratch only)
C    XENORM   MXEVAL    R   Normalized form of XEVAL(*) for Fourier methods
C    MAXIM    MXTPTS    L   MAXIM(J)=T if the Jth turning point found
C                           is a maximum, else it is a minimum
C    XTPTS    MXTPTS    R   Turning point abscissae
C
C FILES USED:
C     LUN     I/O/S   DESCRIPTION
C    LUNCRT     O     For prompts to user, and error messages
C    LUNDAT     I     Data to be analyzed
C    LUNKBD     I     For user responses to prompts
C    LUNOUT     O     For printed output (coefficients, etc.)
C    LUNPLT     O     Several sets of curves in QPLOT format
C
C EXTERNAL REFERENCES:
C   MODULE     SOURCE      DESCRIPTION
C   CSDVAL   INTERPLIB  Cubic spline evaluation (y, y', and y")
C   DCSEVU      IMSL    Calculates 1st derivatives of given cubic spline
C   DSTRIB   INTERPLIB  Generates abscissas for evaluating smoothed curves at
C   FD12K    PROFILLIB  3-point finite difference derivatives
C   FSARBE   INTERPLIB  Evaluates Fourier series calculated by FSARBU/FSARBN
C   FSARBN   INTERPLIB  Fourier coefs for arbitrary nonuniformly-spaced data
C   FSARBU   INTERPLIB  Fourier coefs for arbitrary uniformly-spaced data
C   FSCALC   INTERPLIB  Repackaged FSCOEF/FSDVAL: all coefs./data pt-only evals.
C   FSCOEF   INTERPLIB  Calculates jth Fourier coefficients explicitly
C   FSDVAL   INTERPLIB  Evaluates partial Fourier series, including derivs.
C   ICSEVU      IMSL    Evaluation of cubic splines output by ICS-routines.
C   ICSFKU      IMSL    Least squares spline - fixed knots
C   ICSSCU      IMSL    Simple smoothing spline
C   ICSVKU      IMSL    Least squares spline - variable knots
C   IQHSCU      IMSL    Quasi-Hermite splines
C   LCSFIT   INTERPLIB  Local cubic spline method with monotonic option
C   LSFIT1   INTERPLIB  Local least squares/spline hybrid method
C   MSFIT    INTERPLIB  Monotonic cubic spline fit
C   PNFIT    INTERPLIB  Fits (one) polynomial in least squares sense
C   PNDVAL   INTERPLIB  Evaluates polynomial and derivs. at arbitrary abscissas
C   SECOND     VMSLIB   Measures CPU time
C   XFORMX   NUMODULES  Transforms interval [X(1),X(N)] to specified [a,b]
C   ZEROS    INTERPLIB  Estimates the zeros of 1st derivatives (for max/mins)
C
C HISTORY:
C
C D.A.Saunders  06/27/83  Added max/min code to SMOOTH
C D.A.S.        07/09/83  Replaced FSCALC with FSCOEF
C C.L.Hooper    08/15/83  Added FSARBU and FSARBE
C D.A.S.        08/19/83  Added FSARBN
C C.L.H.        08/24/83  Added IQHSCU
C D.A.S.        06/10/87  Generalized for looking at derivatives (with
C                         max/mins optional)
C   "           06/12/87  Installed finite difference derivative option
C   "           06/16/87  Installed FSCALC and provided for 2nd derivs.
C   "           11/05/87  Installed MSFIT/CSDVAL
C   "           11/24/87  Installed LSFIT1/CSDVAL
C   "           08/26/88  Bug fix involving NFPP=0.  Linked display of
C                         data with display of fit so that both may be
C                         suppressed if f' and f" are the main interest.
C   "           07/12/89  FSARBE call needs B(0), not B(1) now.
C   "           03/13/90  Installed LCSFIT (no 2nd derivatives though).
C   "           03/06/17  Dinesh wanted to try locating maxima and
C                         minima within a very large radiation dataset.
C   "           02/12/23  Minor edits after a cfdtools module version
C                         was found to misbehave.
C
C AUTHOR:  David Saunders, NASA Ames/Sterling Software, June/July 1983.
C          Now with AMA, Inc. at NASA, ARC.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C  *  Constants:

      INTEGER
     +   LUNCRT, LUNDAT, LUNKBD, LUNOUT, LUNPLT, MXEVAL, MXFCOF,
     +   MXMENU, MXNXK, MXNDEG, MXPTS, MXWORK, MXTPTS
      PARAMETER
     +  (LUNCRT = 6,
     +   LUNDAT = 1,
     +   LUNKBD = 5,
     +   LUNOUT = 7,
     +   LUNPLT = 4,
     +   MXEVAL = 200000,
     +   MXFCOF = 256,
     +   MXMENU = 13,
     +   MXNXK  = 28,
     +   MXNDEG = 8,
     +   MXPTS  = 100000,
     +   MXWORK = MAX (7 * MXPTS + 14,
     +                 MXPTS * (MXNDEG + 4),
     +                 (MXPTS * (MXNXK + 6))),
     +   MXTPTS = 256)

      REAL
     +   ZERO
      PARAMETER
     +  (ZERO = 0.E+0)

      CHARACTER
     +   BLANK * 1
      PARAMETER
     +  (BLANK = ' ')

C  *  Variables:

      INTEGER
     +   I, IER, INCR, INCR0, ION, J, M, METHOD, NDEG, NEVAL, NEARPTS,
     +   NFCOF, NFP, NFPP, NTPTS, NX, NXK

      REAL
     +   A(0:MXFCOF), B(0:MXFCOF), COEFS(MXPTS,3), CPU, DX, EPS, ERROR,
     +   F(MXPTS), FEVAL(MXEVAL), FNORM(MXPTS), FPEVAL(MXEVAL),
     +   FPPEVAL(MXEVAL), FSM(MXPTS), RMSDEV, SCALE, SM, TIME1, TIME2,
     +   TWOPI, WGHTS(MXPTS), WORK(MXWORK), X(MXPTS), XEVAL(MXEVAL),
     +   XK(MXPTS), XENORM(MXEVAL), XNORM(MXPTS), XTPTS(MXTPTS)

      LOGICAL
     +   CLOSED, CYCLIC, DEFAULT, EOF, FOURIER, MAXIM(MXTPTS), MAXMINS,
     +   OLDXS, QUIT, SHOFIT, SHOFP, SHOFPP

      CHARACTER
     +   DATAFILE * 96, MENU(0 : MXMENU) * 61, SMETHOD * 8, MODE * 1,
     +   SUBTITLE * 96, TITLE * 96

C  *  Data:

      DATA
     +   MENU
     + /'  0: Help',
     +  '  1: ICSSCU (spline with smoothing param; SM=0 is like CSFIT)',
     +  '  2: ICSFKU (fixed knot least squares spline)',
     +  '  3: ICSVKU (variable knot least squares spline)',
     +  '  4: PNFIT  (least squares polynomial)',
     +  '  5: FSCOEF (Fourier series by recursion - uniform/periodic)',
     +  '  6: FSARBU (F.S. by periodic repetition - uniform/arbitrary)',
     +  '  7: FSARBN (F.S. by periodic repetition - nonuniform data)',
     +  '  8: IQHSCU (quasi-Hermite spline)',
     +  '  9: FD12K  (derivatives by 3-pt. finite differences)',
     +  ' 10: FSCALC (as for FSCOEF but ALL coefs; data-pt-only evals)',
     +  ' 11: MSFIT  (monotonic cubic spline; with CSDVAL)',
     +  ' 12: LSFIT1 (local least squares/cubic spline hybrid method)',
     +  ' 13: LCSFIT (local cubic spline; low storage; f, f''; no f")'/

      DATA
     +   WGHTS /MXPTS * 1.E+0/

C  *  Execution:

      TWOPI = 4.E+0 * ASIN (1.E+0)

      WRITE (LUNCRT, 1001)
     +   ' Program MAXMIN provides a look at f(x), f''(x), and f"(x),',
     +   ' and at maxima and minima if desired.', BLANK

C  *  Open QPLOTable and PRINTable output files:

      OPEN (UNIT=LUNPLT, FILE='maxmin.plt', STATUS='UNKNOWN')
      OPEN (UNIT=LUNOUT, FILE='maxmin.out', STATUS='UNKNOWN')

      CALL OPENER (LUNCRT, 'Data file name? ',
     +              LUNKBD, DATAFILE, LUNDAT, 'OLD')

      READ  (LUNDAT, 1001) TITLE
      WRITE (LUNCRT, 1002) BLANK, 'Dataset title found:', TITLE,
     +   BLANK,
     +   'Plot frame SUBtitles are derived from the methods selected.'
      CALL READS (LUNCRT, 'Enter a plot title (<CR>=dataset title):',
     +             LUNKBD, TITLE, DEFAULT, QUIT)

      WRITE (LUNOUT, 1002)
     +   'Program MAXMIN: A look at derivatives, etc.', BLANK, TITLE

  130 CONTINUE

C  *  Select quantities to be computed/plotted:

      MAXMINS = .FALSE.
      CALL READY (LUNCRT,
     +   'Do you want to estimate maxima/minima? (Y/N; <CR>=No):   ',
     +   LUNKBD, MAXMINS, DEFAULT, QUIT)

      SHOFIT = .FALSE.
      CALL READY (LUNCRT,
     +   'Plot data and fit(s)? (No = just derivatives; <CR>=No):  ',
     +   LUNKBD, SHOFIT, DEFAULT, QUIT)

      NFP = 0
      SHOFP = .FALSE.
      CALL READY (LUNCRT,
     +   'Do you want to plot the 1st derivatives? (Y/N; <CR>=No): ',
     +   LUNKBD, SHOFP, DEFAULT, QUIT)

      NFPP = 0
      SHOFPP = .FALSE.
      CALL READY (LUNCRT,
     +   'Do you want to plot the 2nd derivatives? (Y/N; <CR>=No): ',
     +   LUNKBD, SHOFPP, DEFAULT, QUIT)

      IF (.NOT. (MAXMINS .OR. SHOFIT .OR. SHOFP .OR. SHOFPP)) THEN
         WRITE (LUNCRT, 1001)
     +      ' Bad move: All options suppressed - try again.'
         GO TO 130
      END IF

  150 NEVAL = MXEVAL / 2
      WRITE (LUNCRT, 1002)
     +   'The following is appropriate for some methods only...'
      CALL READI (LUNCRT,
     +   'Enter number of pts. at which to evaluate function/derivs: ',
     +   LUNKBD, NEVAL, DEFAULT, QUIT)
      IF (NEVAL .GT. MXEVAL .OR. NEVAL .LT. 5) THEN
         WRITE (LUNCRT, 1004)
     +      ' *** Sorry - too many or too few.  Limit is ', MXEVAL
         GO TO 150
      END IF

  200 CONTINUE

C  *     Select a smoothing scheme:

         METHOD = 0
         CALL READI (LUNCRT,
     +      'Enter analysis method, ^D to quit, or <CR> for Help: ',
     +      LUNKBD, METHOD, DEFAULT, QUIT)
         IF (QUIT) GO TO 999
         IF (METHOD.LE.0 .OR. METHOD.GT.MXMENU) THEN
            WRITE (LUNCRT, 1001) MENU
            GO TO 200
         END IF

         FOURIER = .FALSE.
         OLDXS = .FALSE.

         GO TO (210,220,230,240,250,260,270,280,290,300,310,320,330)
     +      METHOD

  210       CALL READR (LUNCRT,
     +         'Enter smoothing parameter SM (SM*NX is used): ',
     +         LUNKBD, SM, DEFAULT, QUIT)
            IF (SM.LT.0.E+0 .OR. SM.GT.10.E+0) GO TO 210
            SUBTITLE = 'ICSSCU: SM = z.zzzzz'
            WRITE (SUBTITLE(14:20), '(F7.5)') SM
            GO TO 400

  220       CALL READI (LUNCRT, 'Enter knot-choosing increment: ',
     +         LUNKBD, INCR0, DEFAULT, QUIT)
            SUBTITLE = 'ICSFKU:  INCR0 = nnn'
            WRITE (SUBTITLE(18:20), '(I3)') INCR0
            GO TO 400

  230       WRITE (LUNCRT, 1001) ' Use of ICSVKU not implemented yet.'
            GO TO 200

  240       CALL READI (LUNCRT, 'Enter degree of polynomial: ',
     +         LUNKBD, NDEG, DEFAULT, QUIT)
            IF (NDEG.LT.1 .OR. NDEG.GT.MXNDEG) GO TO 240
            SUBTITLE = 'PNFIT: Degree = n'
            WRITE (SUBTITLE(17:17), '(I1)') NDEG
            GO TO 400

  250       CALL READI (LUNCRT,
     +         'Enter N defining number of Fourier coefs: ',
     +         LUNKBD, NFCOF, DEFAULT, QUIT)
            IF (NFCOF.LT.1 .OR. NFCOF.GT.MXFCOF) GO TO 250
            SUBTITLE = 'FSCOEF:  N = nnn'
            WRITE (SUBTITLE(14:16), '(I3)') NFCOF
            FOURIER = .TRUE.
            GO TO 400

  260       CALL READI (LUNCRT,
     +         'Enter N defining number of Fourier coefs: ',
     +         LUNKBD, NFCOF, DEFAULT, QUIT)
            IF (NFCOF.LT.1 .OR. NFCOF.GT.MXFCOF) GO TO 260
            SUBTITLE = 'FSARBU:  N = nnn'
            WRITE (SUBTITLE(14:16), '(I3)') NFCOF
            FOURIER = .TRUE.
            GO TO 400

  270       CALL READI (LUNCRT,
     +         'Enter N defining number of Fourier coefs: ',
     +         LUNKBD, NFCOF, DEFAULT, QUIT)
            IF (NFCOF.LT.1 .OR. NFCOF.GT.MXFCOF) GO TO 270
            SUBTITLE = 'FSARBN:  N = nnn'
            WRITE (SUBTITLE(14:16), '(I3)') NFCOF
            FOURIER = .TRUE.
            GO TO 400

  280       SUBTITLE = 'IQHSCU'
            GO TO 400

  290       SUBTITLE = 'FD12K '
            OLDXS = .TRUE.
            GO TO 400

  300       SUBTITLE = 'FSCALC'
            FOURIER = .TRUE.
            OLDXS = .TRUE.
            GO TO 400

  310       SUBTITLE = 'MSFIT / CSDVAL'
            GO TO 400

  320       CONTINUE           ! LSFIT1 method
            CLOSED = .FALSE.
            MODE = 'M'
            CALL READC (LUNCRT,
     +         'Hermite or monotone spline? (H/M; <CR>=Mono) ',
     +         LUNKBD, MODE, DEFAULT, QUIT)

            IF (MODE .EQ. 'M') THEN
               SUBTITLE = 'LSFIT1 (monotone)'
               SMETHOD = 'MONOTONE'
            ELSE
               SUBTITLE = 'LSFIT1 (Hermite)'
               SMETHOD = 'HERMITE'
            END IF

  322       NDEG = 3
            CALL READI (LUNCRT,
     +         'Degree of local lst.sqrs. polyn. fits? <CR>=3: ',
     +         LUNKBD, NDEG, DEFAULT, QUIT)
            IF (NDEG .LT. 1 .OR. NDEG .GT. 4) GO TO 322

  324       NEARPTS = NDEG + 3
            CALL READI (LUNCRT,
     +         'Number of pts. in each local fit? <CR>=deg.+3: ',
     +         LUNKBD, NEARPTS, DEFAULT, QUIT)
            IF (NEARPTS .LT. NDEG+1 .OR. NEARPTS .GT. 10)
     +         GO TO 324                  ! LSFIT1 has 10 hard-coded

            WRITE (SUBTITLE (20:30), 1008) 'Degree:', NDEG
            WRITE (SUBTITLE (35:48), 1008) 'Neighbors:', NEARPTS

            CYCLIC = .FALSE.
            CALL READY (LUNCRT,
     +          'Impose periodic boundary conditions? (<CR>=No) ',
     +          LUNKBD, CYCLIC, DEFAULT, QUIT)

            GO TO 400

  330       CONTINUE           ! LCSFIT method
            MODE = 'M'
            CALL READC (LUNCRT,
     +         'Monotone spline?  or loose fit? (M/L; <CR>=Mono) ',
     +         LUNKBD, MODE, DEFAULT, QUIT)

            IF (MODE .EQ. 'M') THEN
               SUBTITLE = 'LCSFIT (monotone)'
               SMETHOD = 'MONOTONE'
            ELSE
               SUBTITLE = 'LCSFIT (loose)'
               SMETHOD = 'BESSEL'
            END IF

            GO TO 400

  400    CONTINUE

C  *     Initialize a new frame in the plot file:

C *****  Note: Line type control assumes 1 or 2 data sets, each with
C *****        non-zero number of maxima and of minima...

         WRITE (LUNPLT, 1001) TITLE, SUBTITLE, 'X', '$'

C  *     Read number of data points NX in first dataset. Skip title.

         READ (LUNDAT, *, ERR=800) NX
         EOF = .FALSE.

         DO 700 WHILE (.NOT.EOF)

C  *        Check for too much data before the damage is done:

            IF (NX.GT.MXPTS) GO TO 810

C  *        Read one set of data.  Separate reads ignore further columns:

            DO I = 1, NX
               READ (LUNDAT, *) X(I), F(I)
            END DO

            IF (OLDXS) THEN
               NEVAL = NX
               DO I = 1, NX
                  XEVAL(I) = X(I)
               END DO
            ELSE

C  *           Generate uniformly-spaced abscissas in data range for
C              evaluating fitted curves and/or derivatives:

               CALL DSTRIB (0, 1, 1, NEVAL, X(1), X(NX), XEVAL)
            END IF

C  *        Apply selected method.

            CALL SECOND (TIME1)

            GO TO (410,420,430,440,450,460,470,480,490,500,510,520,530)
     +         METHOD

  410          CONTINUE

C  *           ICSSCU: Cubic spline with smoothing parameter:

               CALL ICSSCU (X, F, WGHTS, NX, SM*NX, FSM,
     +                      COEFS, MXPTS, WORK, IER)
               IF (IER.NE.0) GO TO 820

               IF (SHOFIT) THEN
                  CALL ICSEVU (X, FSM, NX, COEFS, MXPTS, XEVAL,
     +                         FEVAL, NEVAL, IER)
                  IF (IER.NE.0) GO TO 880
               END IF

               IF (SHOFP ) NFP = NEVAL
               IF (SHOFPP) NFPP = NEVAL
               IF (SHOFP .OR. SHOFPP) THEN
                  CALL DCSEVU (X, FSM, NX, COEFS, MXPTS, XEVAL,
     +                         FPEVAL, NFP, FPPEVAL, NFPP, IER)
                  IF (IER.NE.0) GO TO 890
               END IF

               GO TO 600

  420          CONTINUE

C  *           ICSFKU: Least squares spline with fixed knots.
C              ICSFKU will accept a maximum of 28 knots.
C              Pick every "nth" data point as knots:

               INCR = INCR0
               NXK  = (NX + (INCR-1)) / INCR

               DO WHILE (NXK.GT.MXNXK)
                  INCR = INCR + 1
                  NXK  = (NX + (INCR-1)) / INCR
               END DO

               XK(1) = X(1)
               J = 1 + INCR

C  *           Picking EVERY X seems to cause ICSFKU to fail, so:

               EPS = (X(NX) - X(NX-1)) / 10.E+0

               DO I = 2, NXK - 1
                  XK(I) = X(J) + EPS
                  J = J + INCR
               END DO
               XK(NXK) = X(NX)

               CALL ICSFKU (X, F, NX, 0, XK, NXK, FSM, COEFS,
     +                      MXPTS, ERROR, WORK, IER)
               IF (IER.NE.0) GO TO 830

               IF (SHOFIT) THEN
                  CALL ICSEVU (XK, FSM, NXK, COEFS, MXPTS, XEVAL,
     +                         FEVAL, NEVAL, IER)
                  IF (IER.NE.0) GO TO 880
               END IF

               IF (SHOFP ) NFP = NEVAL
               IF (SHOFPP) NFPP = NEVAL
               IF (SHOFP .OR. SHOFPP) THEN
                  CALL DCSEVU (XK, FSM, NXK, COEFS, MXPTS, XEVAL,
     +                         FPEVAL, NFP, FPPEVAL, NFPP, IER)
                  IF (IER.NE.0) GO TO 890
               END IF

               GO TO 600

  430          CONTINUE

C  *           ICSVKU: Least squares spline with variable knots.
C              Cannot get here, but...

               GO TO 200
               
  440          CONTINUE

C  *           PNFIT/PNDVAL: Polynomial by linear least squares:

               CALL PNFIT (NX, X, F, NDEG, .FALSE.,
     +                     MXWORK, WORK, COEFS, RMSDEV, IER)
               IF (IER.NE.0) GO TO 850

               CALL PNDVAL (NDEG, COEFS, NEVAL, XEVAL, FEVAL,
     +                      FPEVAL, FPPEVAL)

               IF (SHOFP ) NFP = NEVAL
               IF (SHOFPP) NFPP = NEVAL
               GO TO 600

  450          CONTINUE

C  *           FSCOEF (Fourier analysis for uniform, periodic data
C              by efficient recursion; [0,2*pi] data range assumed;
C              necessary transformations made here).

               M = NX - 1
               DO J = 0, NFCOF
                  CALL FSCOEF (3, M, F, J, A(J), B(J), IER)
                  IF (IER.NE.0) GO TO 870
               END DO

C  *           FSCOEF assumes abscissas are in the range [0,2*pi]:

               CALL XFORMX (NEVAL, XEVAL, ZERO, TWOPI, XENORM)

               ION = 0
               IF (SHOFIT) ION = 100
               IF (SHOFP ) ION = ION + 10
               IF (SHOFP ) NFP = NEVAL
               IF (SHOFPP) ION = ION + 1
               IF (SHOFPP) NFPP= NEVAL

               CALL FSDVAL (3, NEVAL, XENORM, NFCOF, A, B, WORK,
     +                      ION, FEVAL, FPEVAL, FPPEVAL)

C  *           Transform the derivatives back to original units:

               SCALE = TWOPI / (X(NX) - X(1))
               IF (SHOFP) THEN
                  DO I = 1, NEVAL
                     FPEVAL(I) = FPEVAL(I) * SCALE
                  END DO
               END IF
               IF (SHOFPP) THEN
                  DO I = 1, NEVAL
                     FPPEVAL(I) = FPPEVAL(I) * SCALE**2
                  END DO
               END IF

               GO TO 600

  460          CONTINUE

C  *           FSARBU: Fourier analysis for arbitrary (but uniformly dis-
C              tributed) data, by transformation and periodic continuation:
C
C              Compute coefficients appropriate in transformed space:

               M = NX - 1
               DO I = 0, NFCOF
                  A(I) = ZERO
               END DO
               B(0) = ZERO

               CALL FSARBU (M, F, FNORM, NFCOF, B(1), IER)
C ***          CALL FSARBU (M, F, F, NFCOF, B(1), IER) works but Fs are lost
               IF (IER.NE.0) GO TO 870

               ION = 0
               IF (SHOFIT) ION = 100
               IF (SHOFP ) ION = ION + 10
               IF (SHOFP ) NFP = NEVAL
               IF (SHOFPP) ION = ION + 1
               IF (SHOFPP) NFPP= NEVAL

               CALL FSARBE (NEVAL, XEVAL, M, X, F, NFCOF, 
     +                      B(0), WORK, ION, FEVAL, FPEVAL, FPPEVAL)
               GO TO 600

  470          CONTINUE

C  *           FSARBN: Fourier analysis for arbitrary, nonuniformly dis-
C              tributed data, by transformation and periodic continuation:
C
C              Compute coefficients appropriate in transformed space:

               M = NX - 1
               DO I = 0, NFCOF
                  A(I) = ZERO
               END DO
               B(0) = ZERO

               CALL FSARBN (M, F, X, FNORM, XNORM, NFCOF, B(1), IER)
               IF (IER.NE.0) GO TO 870

               ION = 0
               IF (SHOFIT) ION = 100
               IF (SHOFP ) ION = ION + 10
               IF (SHOFP ) NFP = NEVAL
               IF (SHOFPP) ION = ION + 1
               IF (SHOFPP) NFPP= NEVAL

               CALL FSARBE (NEVAL, XEVAL, M, X, F, NFCOF, 
     +                      B(0), WORK, ION, FEVAL, FPEVAL, FPPEVAL)
               GO TO 600

  480          CONTINUE

C  *           IQHSCU: Quasi-Hermite spline

               CALL IQHSCU (X, F, NX, COEFS, MXPTS, IER)
               IF (IER.NE.0) GO TO 910

               IF (SHOFIT) THEN
                  CALL ICSEVU (X, F, NX, COEFS, MXPTS, XEVAL,
     +                         FEVAL, NEVAL, IER)
                  IF (IER.NE.0) GO TO 880
               END IF

               IF (SHOFP ) NFP = NEVAL
               IF (SHOFPP) NFPP = NEVAL
               IF (SHOFP .OR. SHOFPP) THEN
                  CALL DCSEVU (X, F, NX, COEFS, MXPTS, XEVAL,
     +                         FPEVAL, NFP, FPPEVAL, NFPP, IER)
                  IF (IER.NE.0) GO TO 890
               END IF

               GO TO 600

  490          CONTINUE

C  *           FD12K: 3-point finite differences:

               CALL FD12K (NX, X, F, FPEVAL, FPPEVAL, FPPEVAL)

               IF (SHOFP ) NFP = NEVAL
               IF (SHOFPP) NFPP = NEVAL
               GO TO 600

  500          CONTINUE

C  *           FSCALC: Fourier analysis for uniform, periodic data by
C              efficient recursion; ALL coefficients; evals. at data pts only

C              For consistency with FSCOEF, etc., assume F(NX) = F(1)
C              is present (not necessary in other applications).

               M = NX - 1
               NEVAL = M
               NFCOF = (M / 2) - 1

C  *           First, calculate the coefficients:

               MODE = 'A'
               CALL FSCALC (MODE, M, F, A, B, WORK, DX, FPEVAL, IER)
               IF (IER. NE. 0) GO TO 920

               IF (SHOFIT) THEN
                  MODE = 'S'
                  CALL FSCALC (MODE, M, F, A, B, WORK, DX, FEVAL, IER)
                  IF (IER. NE. 0) GO TO 920
               END IF

               DX = X(2) - X(1)

               IF (SHOFP) THEN
                  NFP = NEVAL
                  MODE = '1'
                  CALL FSCALC (MODE, M, F, A, B, WORK, DX, FPEVAL, IER)
                  IF (IER. NE. 0) GO TO 920
               END IF

               IF (SHOFPP) THEN
                  NFPP = NEVAL
                  MODE = '2'
                  CALL FSCALC (MODE, M, F, A, B, WORK, DX, FPPEVAL,IER)
                  IF (IER. NE. 0) GO TO 920
               END IF

               GO TO 600

  510          CONTINUE

C  *           MSFIT: Monotonic cubic spline:

               CALL MSFIT (NX, X, F, .FALSE., COEFS(1,1), COEFS(1,2),
     +                     COEFS(1,3), IER)
               IF (IER.NE.0) GO TO 930

               CALL CSDVAL (NX, X, F, NEVAL, XEVAL, COEFS(1,1),
     +                      COEFS(1,2), COEFS(1,3), FEVAL, FPEVAL,
     +                      FPPEVAL)

               IF (SHOFP ) NFP = NEVAL
               IF (SHOFPP) NFPP = NEVAL
               GO TO 600

  520          CONTINUE

C  *           LSFIT1: Local least squares / spline hybrid method.
C              Use existing FNORM array for smoothed function values.

               CALL LSFIT1 (NX, X, F, NDEG, NEARPTS, SMETHOD, CYCLIC,
     +                      FNORM, COEFS(1,1), COEFS(1,2), COEFS(1,3),
     +                      IER)
               IF (IER.NE.0) GO TO 940

               CALL CSDVAL (NX, X, FNORM, NEVAL, XEVAL, COEFS(1,1),
     +                      COEFS(1,2), COEFS(1,3), FEVAL, FPEVAL,
     +                      FPPEVAL)

               IF (SHOFP ) NFP = NEVAL
               IF (SHOFPP) NFPP = NEVAL
               GO TO 600

  530          CONTINUE

C  *           LCSFIT: Local cubic spline method (low-storage; monotone option).

               CALL LCSFIT (NX, X, F, .TRUE., SMETHOD, NEVAL, XEVAL,
     +                      FEVAL, FPEVAL)

               IF (SHOFP ) NFP = NEVAL        ! Leave NFPP = 0
               GO TO 600

  600       CONTINUE

C  *        Save the results for this method:

            CALL SECOND (TIME2)
            CPU = TIME2 - TIME1
            WRITE (LUNCRT, 1007) 'CPU seconds for current method: ',
     +         CPU

            IF (FOURIER) THEN
               WRITE (LUNOUT, 1004) SUBTITLE (1:6) // ': NC = ', NFCOF
               WRITE (LUNOUT, 1001) ' A(0:NC):'
               WRITE (LUNOUT, 1005) (A(I), I=0, NFCOF)
               WRITE (LUNOUT, 1001) ' B(0:NC):'
               WRITE (LUNOUT, 1005) (B(I), I=0, NFCOF)
            END IF

            IF (SHOFIT) THEN

C  *           Show the data points and the fitted function values:

               WRITE (LUNPLT, 1001) ' $OPTIONS',
     +                              ' LINE = ''SYMBOLS'',',
     +                              ' LEGEND = ''Input data'',',
     +                              ' $END'
               WRITE (LUNPLT, 1003) (X(I), F(I), I = 1, NX)
               WRITE (LUNPLT, 1001) 'END CURVE'

               IF (SUBTITLE(1 : 5) .NE. 'FD12K') THEN
                  WRITE (LUNPLT, 1001) ' $OPTIONS',
     +                                 ' LINE = ''SOLID'',',
     +                                 ' LEGEND = ''Curve fit'',',
     +                                 ' $END'
                  WRITE (LUNPLT, 1003) (XEVAL(I), FEVAL(I), I=1,NEVAL)
                  WRITE (LUNPLT, 1001) 'END CURVE'
               END IF
            END IF

            IF (MAXMINS) THEN

C  *           Estimate turning points by locating zeros of 1st derivative:

               CALL ZEROS (NEVAL, XEVAL, FPEVAL, MXTPTS, NTPTS, XTPTS,
     +                     MAXIM, IER)
               IF (IER.NE.0) GO TO 900

C  *           Separate the maxima and minima as two curves (to get all
C              symbols for each type of turning point the same):

               WRITE (LUNPLT, 1001) ' $OPTIONS',
     +                              ' LINE = ''SYMBOLS'',',
     +                              ' LEGEND = ''Maxima'',',
     +                              ' $END'
               DO I = 1,NTPTS
                  IF (MAXIM(I)) WRITE (LUNPLT, 1003)
     +               XTPTS(I), FEVAL(1)
               END DO
               WRITE (LUNPLT, 1001) 'END CURVE'

               WRITE (LUNPLT, 1001) ' $OPTIONS',
     +                              ' LINE = ''SYMBOLS'',',
     +                              ' LEGEND = ''Minima'',',
     +                              ' $END'
               DO I = 1,NTPTS
                  IF (.NOT.MAXIM(I)) WRITE (LUNPLT, 1003)
     +               XTPTS(I), FEVAL(1)
               END DO
               WRITE (LUNPLT, 1001) 'END CURVE'
            END IF

C  *        Plot derivatives on the same frame as the data.
C           (Is there a better arrangement?)

            IF (SHOFP .AND. NFP .GT. 0) THEN
               WRITE (LUNPLT, 1001) ' $OPTIONS',
     +                              ' LINE = ''DASH'',',
     +                              ' LEGEND = ''1st derivative'',',
     +                              ' $END'
               WRITE (LUNPLT, 1003) (XEVAL(I), FPEVAL(I), I=1, NEVAL)
               WRITE (LUNPLT, 1001) 'END CURVE'
            END IF

            IF (SHOFPP .AND. NFPP .GT. 0) THEN
               WRITE (LUNPLT, 1001) ' $OPTIONS',
     +                              ' LINE = ''DOT'',',
     +                              ' LEGEND = ''2nd derivative'',',
     +                              ' $END'
               WRITE (LUNPLT, 1003) (XEVAL(I), FPPEVAL(I), I=1,NEVAL)
               WRITE (LUNPLT, 1001) 'END CURVE'
            END IF

C  *        Read NX for next case, if any:

            READ (LUNDAT, *, END=690, ERR=990) NX
            GO TO 700
C
  690       EOF = .TRUE.
  700    CONTINUE

C  *     Go back for another curve option:

         WRITE (LUNPLT, 1001) 'END FRAME'
         REWIND LUNDAT
         READ (LUNDAT, *)
      GO TO 200

C  *  Error handling:

  800 WRITE (LUNCRT, 1002) 'Error on first (re)read from data file.'
      GO TO 999
  810 WRITE (LUNCRT, 1004) 'Too many data points found. NX = ', NX
      GO TO 999
  820 WRITE (LUNCRT, 1004) 'Error in ICSSCU. IER =', IER
      GO TO 999
  830 WRITE (LUNCRT, 1004) 'Error in ICSFKU. IER =', IER
      GO TO 999
  850 WRITE (LUNCRT, 1004) 'Error in PNFIT.  IER =', IER
      GO TO 999
  870 WRITE (LUNCRT, 1001) 'Error in FSCOEF, FSARBU, or FSARBN.'
      WRITE (LUNCRT, 1004) 'NFCOF = ', NFCOF,
     +                     'IER   = ', IER,
     +                     'M     = ', M
      GO TO 999
  880 WRITE (LUNCRT, 1004) 'Error in ICSEVU. IER =', IER
      GO TO 999
  890 WRITE (LUNCRT, 1004) 'Error in DCSEVU. IER =', IER
      GO TO 999
  900 WRITE (LUNCRT, 1004) 'Error in ZEROS.  IER =', IER
      GO TO 999
  910 WRITE (LUNCRT, 1004) 'Error in IQHSCU. IER =', IER
      GO TO 999
  920 WRITE (LUNCRT, 1004) 'Error in FSCALC. IER =', IER
      WRITE (LUNCRT, 1006) ' MODE = ', MODE
      GO TO 999
  930 WRITE (LUNCRT, 1004) 'Error in MSFIT. IER =', IER
      GO TO 999
  940 WRITE (LUNCRT, 1004) 'Error in LSFIT1. IER =', IER
      GO TO 999

  990 WRITE (LUNCRT, 1002) 'Error looking for next dataset.'
      GO TO 999

  999 CONTINUE
      WRITE (LUNCRT, 1001) ' ',
     +                     ' *** Plottable data:  maxmin.plt ***',
     +                     ' *** Coefficients  :  maxmin.out ***'

CCC   STOP ' '  ! Avoid system-dependent behavior

C  *  Formats:

 1001 FORMAT (A)
 1002 FORMAT (1X, A)
 1003 FORMAT (2 (1X, ES14.6))
 1004 FORMAT (/, 1X, A, I7)
 1005 FORMAT (5 (1X, ES14.6))
 1006 FORMAT (A, A)
 1007 FORMAT (1X, A, F7.2)
 1008 FORMAT (A, I4)

      END
