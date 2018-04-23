C+------------------------------------------------------------------------------
C
      PROGRAM DISTRIBUTE
C
C  PURPOSE:
C        DISTRIBUTE is a driving program for the available 1-D grid point
C     distribution utilities.  It allows generation of one or more distri-
C     butions on the specified interval(s) [a, b].  Results may be saved
C     in one or both of two ways: one for plotting via QPLOT, and one that
C     is more convenient for use by programs such as SMOOTH or PROFILE or
C     some other grid generator.
C
C  METHOD:
C        The indicated output file(s) are opened, then an indefinite loop
C     is entered over the desired cases, which are selected one at a time
C     from a menu.  Some selections may require additional interactive
C     inputs.  The interval [a, b] may be changed for each case.  Defaults
C     are [0, 1] initially, and the previous interval subsequently.
C
C        All distributions generated in one run of DISTRIBUTE will appear on
C     the same plot frame (if plotting is requested).  Unconnected symbols are
C     the normal choice of line type, though the symbols may be suppressed.
C     The plot legend is used (with internal writes in some cases) to make
C     each plotted distribution as self-descriptive as possible.  Other plot
C     parameters such as max/min values are defaulted - QPLOT handles them.
C               
C  NOTES:
C     o  QPLOT is a general purpose plotting package developed by Robert
C        Kennelly (Sterling/Aerodynamics Division, NASA Ames).
C
C  FILES USED:
C     LUNCRT   Screen prompts
C     LUNKBD   Keyboard entries
C     LUNDAT   Output file in form usable by SMOOTH, etc. (one column with
C              title and "NX" for lines 1 and 2)
C     LUNPLT   Output file in QPLOT format (using YCOL=0 option)
C
C  EXTERNAL REFERENCES:
C     The following are now down a level in DISTRIB:
C     ARBDIS   Arbitrary-shape distributions
C     BLGRID   Hybrid geometric + Vinokur distribution for boundary layers
C     CONDIS   Exponential-type distribution with specified interior pt.
C     DSTRIB   Uniform and several sinusoidal-type distributions
C     EXPDIS   Exponential-type distribution with specified stretching param.
C     EXPDIS2  Variant of EXPDIS with first or last dX specified
C     GEODIS   Modified geometric distributions (1-sided)
C     GEODIS2  Modified geometric distributions (2-sided)
C     HTDIS2   Hyperbolic tangent method (first and last slope or dX specified)
C     Also:
C     READER   Prompting utility (entry pts. READI, etc.)
C     SCANNR   Needed to find where to add text to legend
C     TOGGLE   Utility for enabling/disabling options
C
C  ENVIRONMENT:  DEC VMS, FORTRAN 77 + minor extensions.
C
C  HISTORY:
C     03/07/88   DAS   Initial implementation (derived from EVALUATE).
C     03/18/88   DAS   Installed GEODIS2.
C     07/21/89   DAS   Provided for plotting dX vs. X (as well as pt. # vs. X).
C     07/24/89   DAS   Installed HTDIS2.
C     08/12/89   DAS   Extracted bulk of code as reusable routine DISTRIB.
C     08/28/89   DAS   ARBDIS workspace requirements reduced thanks to LCSFIT.
C     11/19/90   DAS   Installed HTDIS (in DISTRIB).
C     04/24/91   DAS   DISTRIB now expects a scale factor in WK (1) for the
C                      methods involving dX > 0. in RPARAM (*).
C     05/02/91   DAS   Added output of maximum growth factor.
C     05/07/91   DAS     "     "     " max., min., average dX (screen only);
C                      permitted suppressing both forms of saved results in
C                      case the screen output is all that is of interest.
C     10/26/91   DAS   Took advantage of QPLOT's integer axis-numbering and
C                      ESCAPE options (to allow for [0, 1], etc. in titles).
C     09/16/00   DAS   Installed BLGRID.
C     03/16/03   DAS   Eliminated CARRIAGECONTROL='LIST' for SGI F90.
C     02/01/14   DAS   Replaced ISAMIN, ISAMAX with MINVAL, MAXVAL.
C
C  AUTHOR:   David Saunders, Sterling Software/NASA Ames, Palo Alto, CA.
C
C-------------------------------------------------------------------------------

      IMPLICIT  NONE

C     Local constants:

      INTEGER
     >   LUNCRT, LUNDAT, LUNKBD, LUNPLT, MXC, MXP, MXPTS, MXWORK
      REAL
     >   HALF, ONE, ZERO
      CHARACTER
     >   BLANK * 1

      PARAMETER
     >  (BLANK  = ' ',
     >   HALF   = 0.5E+0,
     >   LUNCRT = 6,
     >   LUNDAT = 1,
     >   LUNKBD = 5,
     >   LUNPLT = 2,
     >   MXC    = 300,                 ! Max. # control points for ARBDIS
     >   MXP    = 5,                   ! Max. # params. needed by some methods
     >   MXPTS  = 1000,                ! Max. # points in one distribution
     >   MXWORK = 2 * MXC + 5 * MXPTS, ! For ARBDIS
     >   ONE    = 1.E+0,
     >   ZERO   = 0.E+0)

C     Local variables:

      INTEGER
     >   FIRST, I, IER, IG, IMX, IMN, IOS, IP (MXP), LAST, MARK, METHOD,
     >   NPTS
      REAL
     >   DX (MXPTS), DXMEAN, GROWTH (MXPTS), P (MXP), WORK (MXWORK),
     >   X (MXPTS), XA, XB
      LOGICAL
     >   ANOTHER, DEFAULT, DXVSX, FIRSTCASE, ONOFF (2), QUERY, QUIT,
     >   SAME, SHOWXS, SYMBOLS, YESDAT, YESPLT
      CHARACTER
     >   ANSWER * 1, FILENAME * 50, LEGEND * 80, OUTPUTS (2) * 46,
     >   START * 1, SUBTITLE * 70, TITLE * 70

C     Procedures:

      EXTERNAL
     >   DISTRIB, READC, READI, READR, READS, READY, SCANNR, TOGGLE

C     Storage:

      DATA
     >   OUTPUTS
     >     /'1:  Save results in QPLOT form?',
     >      '2:  Save results in one-column form for reuse?'/

C     Execution:


      WRITE (LUNCRT, '(/, (A))')
     >   ' Program DISTRIBUTE generates a variety of 1-D distributions.'
     >  ,' Results may be saved in two ways (either, neither, or both):'

      ONOFF (1) = .TRUE.
      ONOFF (2) = .FALSE.

      CALL TOGGLE (2, OUTPUTS, ONOFF, LUNCRT, LUNKBD)

      YESPLT = ONOFF (1)
      YESDAT = ONOFF (2)
         
      TITLE = 'Grid point distribution'

      IF (YESPLT) THEN

C        Open a QPLOT file and set up the plot frame:

         OPEN (UNIT = LUNPLT, FILE = 'distribute.plt',
     >         STATUS = 'UNKNOWN', IOSTAT = IOS)

         IF (IOS .NE. 0) THEN
            WRITE (LUNCRT, '(/, A)')
     >         ' Could not open distribute.plt.'
            GO TO 999
         END IF

         ANSWER = 'D'
         CALL READC (LUNCRT,
     >      'Plot dX vs. X or Ordinal # vs. X?  [D/O; <CR> = D] ',
     >      LUNKBD, ANSWER, DEFAULT, QUIT)
         DXVSX = ANSWER .EQ. 'D'

         CALL READS (LUNCRT,
     >      'Enter a plot title.  [<CR> = Grid point distribution]',
     >      LUNKBD, TITLE, DEFAULT, QUIT)

         SUBTITLE = '" "'    ! QPLOT skips blank lines
         I = 3
         CALL READS (LUNCRT,
     >      'Enter a subtitle.  [<CR> = none]',
     >      LUNKBD, SUBTITLE, DEFAULT, QUIT)
         IF (.NOT. DEFAULT) I = LEN (SUBTITLE)

         WRITE (LUNPLT, 1001) TITLE, SUBTITLE (1 : I), 'X'
         IF (DXVSX) THEN
            WRITE (LUNPLT, 1001) 'dX'
         ELSE
            WRITE (LUNPLT, 1001) 'Point number'
         END IF

         SYMBOLS = .TRUE.
         CALL READY (LUNCRT,
     >      'Do you want symbols at the plotted points?  [<CR> = Yes] ',
     >      LUNKBD, SYMBOLS, DEFAULT, QUIT)

         SHOWXS = .FALSE.
         IF (DXVSX)
     >      CALL READY (LUNCRT,
     >      'Do you want the Xs marked along the X axis?  [<CR> = No] ',
     >      LUNKBD, SHOWXS, DEFAULT, QUIT)
      END IF

      IF (YESDAT) THEN

  150    FILENAME = 'distribute.dat'
         CALL READS (LUNCRT,
     >      'Enter filename for generated grid points, or default it: ',
     >      LUNKBD, FILENAME, DEFAULT, QUIT)
         IF (QUIT) GO TO 999

         OPEN (UNIT=LUNDAT, FILE=FILENAME, STATUS='NEW', IOSTAT=IOS)
         IF (IOS .NE. 0) GO TO 150

      END IF



C ------------------- Begin loop over distributions --------------------

      NPTS = 51
      XA = ZERO
      XB = ONE
      QUERY = .TRUE.        ! Always prompt for method and parameters
      START = 'A'           ! Auto-start iterative methods
      FIRSTCASE = .TRUE.

  200 CONTINUE

         IF (FIRSTCASE) THEN

  220       CALL READI (LUNCRT,
     >         'How many points in the distribution?  Default = 51: ',
     >         LUNKBD, NPTS, DEFAULT, QUIT)
            IF (NPTS .LE. 1 .OR. NPTS .GT. MXPTS) GO TO 220

            CALL READR (LUNCRT,
     >         'Enter lower end of interval [a, b].  Default = 0.: ',
     >         LUNKBD, XA, DEFAULT, QUIT)

            CALL READR (LUNCRT,
     >         'Enter upper end of interval [a, b].  Default = 1.: ',
     >         LUNKBD, XB, DEFAULT, QUIT)

            WRITE (LUNCRT, '(/, A)') ' Distribution options:'
            METHOD = 12          ! 1 for DSTRIB is not unique.

         ELSE        ! Further cases will commonly use same N and [a, b].

            SAME = .TRUE.
            CALL READY (LUNCRT,
     >         'Same interval and same # points?  (Y/N; <CR>=Yes): ',
     >         LUNKBD, SAME, DEFAULT, QUIT)
            IF (QUIT) GO TO 900

            IF (.NOT. SAME) THEN

  240          CALL READI (LUNCRT,
     >            'How many points?  Default = same: ',
     >            LUNKBD, NPTS, DEFAULT, QUIT)
               IF (QUIT) GO TO 900
               IF (NPTS .LE. 1 .OR. NPTS .GT. MXPTS) GO TO 240

               CALL READR (LUNCRT,
     >            'Lower end of interval [a, b]?  Default = same: ',
     >            LUNKBD, XA, DEFAULT, QUIT)

               CALL READR (LUNCRT,
     >            'Upper end of interval [a, b]?  Default = same: ',
     >            LUNKBD, XB, DEFAULT, QUIT)
            END IF
         END IF


C        Select a method and apply it:
C        -----------------------------

         WORK (1) = ZERO  ! Since no scaling of dX > 0. values is done
                          ! in this application.

         CALL DISTRIB (METHOD, QUERY, NPTS, XA, XB, LUNCRT, LUNKBD,
     >      IP, P, START, MXWORK, WORK, X, LEGEND, IER)
         IF (IER .EQ. -1) GO TO 900                       !Quit


C        Determine the maximum growth factor and the max. & min. dX:

         DO 300, I = 1, NPTS - 1
            DX (I) = X (I + 1) - X (I)
  300    CONTINUE

         GROWTH (1) = ZERO
         DO 310, I = 2, NPTS - 1
            GROWTH (I) = DX (I) / DX (I - 1)
            IF (GROWTH (I) .LT. ONE / GROWTH (I))
     >         GROWTH (I) = -ONE / GROWTH (I)
  310    CONTINUE

         IMX = MAXVAL (DX(1:NPTS-1))
         IMN = MINVAL (DX(1:NPTS-1))
         IG  = MAXVAL (GROWTH(1:NPTS-1))
         DXMEAN = (XB - XA) / REAL (NPTS-1)

         WRITE (LUNCRT, 1007) DX (IMX), IMX, DX (IMN), IMN, DXMEAN,
     >                        GROWTH (IG), IG

C        Save results:
C        -------------

         IF (YESPLT) THEN      ! Save QPLOTable results:

C           Add max. growth factor to the legend.  Assume enough room:

            FIRST = 1
            LAST = LEN (LEGEND)
            CALL SCANNR (LEGEND, FIRST, LAST, MARK)
            WRITE (LEGEND (LAST + 3 :), 1008) GROWTH (IG), IG, NPTS

            WRITE (LUNPLT, 1001) ' $OPTIONS'
            WRITE (LUNPLT, 1003) ' ESCAPE=', '{}'  ! To allow [a, b].

            IF (DXVSX) THEN    ! Plot dX at the mid-point of each interval:

               IF (SYMBOLS) THEN
                  WRITE (LUNPLT, 1003) ' LINE=', 'SYMBOLS'
               ELSE
                  WRITE (LUNPLT, 1003) ' LINE=', 'SOLID'
               END IF
               WRITE (LUNPLT, 1003) ' LEGEND=', LEGEND
               WRITE (LUNPLT, 1001) ' $END'
               WRITE (LUNPLT, 1006)
     >            (HALF * (X (I) + X (I+1)), DX (I), I = 1, NPTS - 1)

               IF (SHOWXS) THEN  ! Show the Xs as '+'s along the X-axis:

                  WRITE (LUNPLT, 1001) 'END CURVE'
                  WRITE (LUNPLT, 1001) ' $OPTIONS'
                  WRITE (LUNPLT, 1003) ' LINE=', 'SYMBOLS'
                  WRITE (LUNPLT, 1001) ' SYMBOL=3,', ' $END'
                  WRITE (LUNPLT, 1006) (X (I), ZERO, I = 1, NPTS)
               END IF

            ELSE  ! Just save the Xs in 1 column - QPLOT supplies the ordinals.

               WRITE (LUNPLT, 1001) ' YCOL=0,'
               WRITE (LUNPLT, 1003) ' YNUMBERS=', 'INTEGER'
               IF (SYMBOLS) WRITE (LUNPLT, 1003) ' LINE=', 'SYMBOLS'
               WRITE (LUNPLT, 1003) ' LEGEND=', LEGEND
               WRITE (LUNPLT, 1001) ' $END'
               WRITE (LUNPLT, 1005) (X (I), I = 1, NPTS)
            END IF
            WRITE (LUNPLT, 1001) 'END CURVE'
         END IF

         IF (YESDAT) THEN    ! Save abscissas for probable reuse:
            IF (FIRSTCASE) WRITE (LUNDAT, 1001) TITLE
            WRITE (LUNDAT, 1004) NPTS
            WRITE (LUNDAT, 1005) (X (I), I = 1, NPTS)
         END IF

         ANOTHER = .TRUE.
         WRITE (LUNCRT, 1001)
         CALL READY (LUNCRT, 'Another case?  [Y/N/EOF; <CR> = Yes] ',
     >      LUNKBD, ANOTHER, DEFAULT, QUIT)
         IF (QUIT .OR. .NOT. ANOTHER) GO TO 900

         FIRSTCASE = .FALSE.
      GO TO 200

C ------------------- End of loop over distributions -------------------


  900 CONTINUE

      IF (YESPLT) WRITE (LUNCRT, '(/, A)')
     >   '    QPLOTable results:       distribute.plt'
      IF (YESDAT) WRITE (LUNCRT, 1002)
     >   '    Distribution abscissas:  ', FILENAME

  999 CONTINUE

C     Formats:

 1001 FORMAT (A)
 1002 FORMAT (A, A)
 1003 FORMAT (A, '''', A, ''',')
 1004 FORMAT (I4)
 1005 FORMAT (G15.7)
 1006 FORMAT (2G15.7)
 1007 FORMAT (/, ' Largest  dX: ', 1P, E12.6, ' at point', I4, /,
     >        ' Smallest dX: ', E12.6, ' at point', I4, /,
     >        ' Average  dX: ', E12.6, 0P, /
     >        ' Largest growth factor: ', F6.2, ' around point', I4)
 1008 FORMAT ('Gmax:', F6.2, ' at pt.', I4, ' of', I4)

      END
