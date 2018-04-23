C+------------------------------------------------------------------------------
C
      PROGRAM WINGSECTIONS
C
C     One-liner:  Prepares a wing surface as 2D sections in B-spline form
C     ----------
C
C     Purpose:
C     --------
C
C        WINGSECTIONS completes a simplified form of wing surface definition
C     from previously-prepared chord/thickness data and associated defining
C     sections.  The input sections are in B-spline form, and the results are
C     output as an expanded file of B-spline sections along with a summary
C     tabulation and an expanded file of XY coordinates for possible use by
C     a flow solver.
C
C     Method:
C     -------
C
C        The simplifying assumption here is that all sections can be represented
C     adequately by B-spline curves with the SAME KNOTS.  (See a version of the
C     SMOOTH program for methods of determining such B-spline representations
C     of related XY datasets.)
C
C        As long as two sections have the same knots, lofting of intermediate
C     sections can be achieved by lofting the control points.  This allows
C     standard algebraic methods to produce a compact surface representation
C     that is readily adapted by NURBS-based CAD/CAM systems as a full 3-D
C     surface.  Ideally, the same procedure could provide an alternative
C     input scheme for flow solver applications, although a CAD system is
C     better suited to geometry refinements such as detailed wing tip design.
C
C        Initially, WINGSECTIONS uses a slightly-modified version of the
C     geometry summary table output by R22OPT for the thickness/chord data.
C     For each station with "LOFT" = 1, a defining section in B-spline curve
C     form is expected.  The defining sections may be normalized (or not).
C     They should be concatenated as one input file.
C
C     Sample thickness/chord data file:
C     ---------------------------------
C
C     OAW-1E Geometry (Squared-tip form)                          <Title>
C     KG         Z       XLE       YLE     CHORD  ...    MAX T  ...  LOFT
C      1 -203.4000  10.87916   0.00000  21.80275  ...  1.27510  ...     0
C      2 -197.1225  10.34812   0.00000  23.46226  ...  1.57230  ...     0
C      3 -190.8450   9.81707   0.00000  25.12178  ...  1.87720  ...     1
C      4 -187.6130   9.54366   0.00000  25.97618  ...  2.06160  ...     0
C      . .........   .......   .......  ........  ...  .......  ...     .
C     43  -61.5653   0.03096   0.00000  55.70325  ...  8.88560  ...     0
C     44  -58.3333   0.00000   0.00000  55.80000  ...  8.92800  ...     1
C      . .........   .......   .......  ........  ...  .......  ...     .
C     58   58.3333   0.00000   0.00000  55.80000  ...  8.92800  ...     1
C     59   61.5653   0.03096   0.00000  55.70325  ...  8.90020  ...     0
C      . .........   .......   .......  ........  ...  .......  ...     .
C     98  187.6130   9.54366   0.00000  25.97618  ...  1.33900  ...     0
C     99  190.8450   9.81707   0.00000  25.12178  ...  1.11590  ...     1
C    100  197.1225  10.34812   0.00000  23.46226  ...  0.88500  ...     0
C    101  203.4000  10.87916   0.00000  21.80275  ...  0.36720  ...     0
C
C     Corresponding defining sections file:
C     -------------------------------------
C
C     OAW-1E-P2-MOD2 Tip-cap Section                                <Title>
C       1.000000000000000       2.000000000000000       5.000000000000000
C       21.00000000000000       21.00000000000000
C
C     0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00
C     0.0000000000000000E+00  0.0000000000000000E+00  9.1353959917732475E-02
C     0.2411894763703840      0.4024486663234111      0.5634075800663552
C     0.7232546095157547      0.8761017106844363      0.9779216600852452
C     1.022438833587189       1.065594150116754       1.154251426767800
C     1.286308448870403       1.428047814045446       1.568080209559026
C     1.710104876249171       1.852958747840903       1.984070742665636
C     2.067673456116899       2.067673456116899       2.067673456116899
C     2.067673456116899       2.067673456116899
C
C      1.000000000000000      0.9766573616010096      0.9167987120914935
C     0.8182593557168120      0.6759316525938030      0.5204441164611799
C     0.3669437255341476      0.2052794078141187      0.1237065524658922
C    -5.3327239734320598E-03 -1.3476529667296819E-03  7.8166813309687295E-02
C     0.1828410957291661      0.2952688069166972      0.4437506124776002
C     0.5805998573300908      0.7177141130728990      0.8390485132513349
C     0.9283616157847692      0.9785103137567076       1.000000000000000
C
C     4.7510568983852863E-03  4.2093382707367472E-03  5.6186343982026685E-03
C    -5.3332570885323248E-04 -9.6907468475426527E-03 -1.7661125500254172E-02
C    -3.1396893930427520E-02 -2.6467468887536488E-02 -2.1503660335958938E-02
C    -1.7591147529822391E-02  1.5159539809190553E-02  4.0424128341521668E-02
C     4.8592957581737368E-02  4.6747271055867917E-02  4.4292514180655268E-02
C     4.6004845731045817E-02  3.6009187829119820E-02  2.3870693882375195E-02
C     1.4528201781536322E-02  9.4037158124965926E-03  7.2510391473770142E-03
C     OAW 60-085-16 (OAW-1 Cabin section)            <Next defining section>
C     1.000000000000000       2.000000000000000       5.000000000000000
C     21.00000000000000       21.00000000000000
C
C     0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00
C     0.0000000000000000E+00  0.0000000000000000E+00  9.1353959917732475E-02
C     0.2411894763703840      0.4024486663234111      0.5634075800663552
C     0.7232546095157547      0.8761017106844363      0.9779216600852452
C     1.022438833587189       1.065594150116754       1.154251426767800
C     1.286308448870403       1.428047814045446       1.568080209559026
C     1.710104876249171       1.852958747840903       1.984070742665636
C     2.067673456116899       2.067673456116899       2.067673456116899
C     2.067673456116899       2.067673456116899                <Same knots>
C
C     1.000000000000000       0.9766573684613700      0.9167986913545925
C     0.8182593719759367      0.6759316192254129      0.5204440817792035
C     0.3669437646467861      0.2052793407910668      0.1237065913369712
C    -5.3327403303544115E-03 -1.3476487841357643E-03  7.8166807256873784E-02
C     0.1828411143917873      0.2952687381718988      0.4437506889818151
C     0.5805997456052916      0.7177141687293529      0.8390484488001356
C     0.9283616554188250      0.9785102941149742       1.000000000000000
C
C    -2.4999999441206455E-03 -1.5764163893617429E-03 -3.8547663530724645E-04
C    -1.0803315819805590E-02 -3.8935721655883260E-02 -5.1152378185592381E-02
C    -5.1869694038945192E-02 -5.0224055012760078E-02 -5.0413761744842232E-02
C    -3.5318966792019309E-02  3.0945843199744209E-02  7.4389077816659148E-02
C     9.5813309764173755E-02  0.1101315601432895      0.1102606125877267
C     9.1098402348239831E-02  6.2800820726638147E-02  3.6244110229431285E-02
C     1.8563923826951780E-02  8.3565019335236754E-03  2.4999999441206455E-03
C     <The next section follows similarly.>
C
C     Points to note:
C     ---------------
C
C        > The 5-word header of each curve must be readable list-directed.
C          It is offset for this reason.  The other blank lines are for
C          human readability only.
C
C        > A compact representation is the goal here.  (Otherwise, a vast
C          collection of points defining the wing section could be supplied
C          to the CAD system.)
C
C     Files used:
C     -----------
C
C     LUNCHORDTHICK   Input chord & thickness data from R22OPT's tabular summary
C     LUNCRT          Screen prompts, error messages, etc.
C     LUNKBD          Keyboard entries
C     LUNCOORDINATES  Discretized form of completed wing, output as R22OPT input
C     LUNINSECTIONS   Input defining sections (B-spline curves matching
C                     the "LOFT=1" stations, in DT_NURBS format)
C     LUNOUTSECTIONS  Complete wing sections in DT_NURBS curve form
C     LUNPRINT        Tabular summary of results
C
C     Environment:
C     ------------
C
C     VAX/VMS; FORTRAN 77, with...
C     > IMPLICIT NONE
C     > Trailing ! comments
C     > Names longer than 6 characters
C
C     History:
C     --------
C
C     02/23/93  D.A.Saunders  Initial implementation.
C     02/27/93    "     "     Added discretized output file for R22OPT, etc.
C     03/07/93    "     "     Added more commentary and checked for mismatches.
C     04/14/93    "     "     Allowed for 72 pts. per surface.
C     06/16/93    "     "     Allowed for 75 control pts. per wing section.
C
C     Author:   David Saunders, Sterling Software/NASA Ames Research Center
C     -------                                     Moffett Field, CA
C
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Procedures:
C     -----------

      EXTERNAL
     >   BSAF2XY,        ! Discretizes a B-spline-type airfoil
     >   BSAFGEOM,       ! Finds B-spline airfoil geometry properties
     >   BSAFNORM,       ! Imposes desired chord, thickness, leading edge
     >   BSCRVRD,        ! Reads one or more B-spline curves from one file
     >   DTRDWR,         ! Saves the output B-spline sections
     >   LINTRP,         ! Does the standard spanwise lofting (linear)
     >   OPENER,         ! File-opening utility
     >   READER          ! Prompting utility

C-------------------------------------------------------------------------------

C     Local constants:
C     ----------------

      INTEGER
     >   LUNCRT, LUNKBD, LUNCHORDTHICK, LUNCOORDINATES, LUNINSECTIONS,
     >   LUNOUTSECTIONS, LUNPRINT, MAXC, MAXCPTS, MAXORDER, MAXSECTIONS,
     >   NL, NU, NXY

      DOUBLE PRECISION
     >   ONE, PI, ZERO

      CHARACTER
     >   BLANK * 1

      PARAMETER
     >  (LUNCRT         = 6,
     >   LUNKBD         = 5,
     >   LUNCHORDTHICK  = 1,       ! See descriptions above
     >   LUNINSECTIONS  = 2,
     >   LUNOUTSECTIONS = 3,
     >   LUNPRINT       = 4,
     >   MAXCPTS        = 75,      ! Max. # control pts. per B-spline section
     >   MAXORDER       = 6,       ! Max. B-spline degree + 1
     >   MAXC           = 5 + MAXORDER + MAXCPTS * 3,
                                   ! Max. # words per nonrational XY "DT" curve
     >   MAXSECTIONS    = 101,     ! Max. # span stations handled
     >   NL             = 72,      ! # points on discretized surfaces
     >   NU             = 72,      !   (lower, upper, and wrap-around)
     >   NXY            = NL + NU - 1,
     >   ONE            = 1.D+0,
     >   PI             = 3.141592654D+0,
     >   ZERO           = 0.D+0,
     >   BLANK          = ' ')

C     Local variables:
C     ----------------

      INTEGER
     >   I, IER, IX1, IY1, K, K1, K2, KDEFINED (MAXSECTIONS),
     >   KG, KL, KORDER, KR, LOFT, NC, NCPTS, NDEFINED, NSECTIONS

      DOUBLE PRECISION
     >   AREATEMP, CAMBER, CHORDTEMP, DELLE, DELTE, DUM, THICKTEMP, TLE,
     >   XCAMBER, XLETEMP, XTEL, XTEU, XTHICKTEMP, YLETEMP, YTEL, YTEU,
     >   AREA (MAXSECTIONS), BSPLINE (MAXC, MAXSECTIONS),
     >   CHORD (MAXSECTIONS), RADIUS (MAXSECTIONS), T (NXY),
     >   THICKNESS (MAXSECTIONS), TWIST (MAXSECTIONS), X (NXY),
     >   XLE (MAXSECTIONS), XTHICK (MAXSECTIONS), Y (NXY),
     >   YLE (MAXSECTIONS), Z (MAXSECTIONS)

      LOGICAL
     >   CR, EOF

      CHARACTER
     >   HEADERS * 1, INFILE * 40, TITLE * 80, TITLEOUT * 80


C     Execution:
C     ----------

      CALL OPENER (LUNCRT, 'Chord/thickness data file?  ',
     >             LUNKBD, INFILE, LUNCHORDTHICK, 'OLD')

      CALL OPENER (LUNCRT, 'Defining B-spline sections file?  ',
     >             LUNKBD, INFILE, LUNINSECTIONS, 'OLD')

      CALL OPENER (LUNCRT, BLANK, LUNKBD, 'wingsections.bsplines',
     >             LUNOUTSECTIONS, 'NEW:LIST')

      CALL OPENER (LUNCRT, BLANK, LUNKBD, 'wingsections.xy',
     >             LUNCOORDINATES, 'NEW:LIST')

      CALL OPENER (LUNCRT, BLANK, LUNKBD, 'wingsections.summary',
     >             LUNPRINT, 'NEW:LIST')

C     Read the chord/thickness data for all stations (till EOF):

      READ (LUNCHORDTHICK, 1010) TITLE, HEADERS
      TITLEOUT (23:) = TITLE     ! See insertion of KG & Z below

      WRITE (LUNCRT, 1015) TITLEOUT
      CALL READS (LUNCRT,
     >   'Title for saved results?  [Above input title]',
     >   LUNKBD, TITLEOUT (23:), CR, EOF)
      IF (EOF) GO TO 999

      NDEFINED = 0
  200 CONTINUE

         READ (LUNCHORDTHICK, *, END=300) KG, Z (KG), XLE (KG),
     >      YLE (KG), CHORD (KG), DUM, DUM, DUM, THICKNESS (KG),
     >      DUM, DUM, DUM, DUM, LOFT

         IF (LOFT .GT. 0) THEN
            NDEFINED = NDEFINED + 1
            KDEFINED (NDEFINED) = KG
         END IF

         IF (KG .LE. MAXSECTIONS)
     >GO TO 200           ! Next station

      GO TO 910           ! Too many stations


  300 CONTINUE

      NSECTIONS = KG
      WRITE (LUNCRT, 1010)

C     Read the defining sections, in DT_NURBS format.

      DO K = 1, NDEFINED
         KG = KDEFINED (K)

         WRITE (LUNCRT, '(A, I3)')
     >      '    Processing defining station ', KG

         CALL BSCRVRD (LUNINSECTIONS, TITLE, MAXC, BSPLINE (1, KG),
     >                 NC, IER)
         IF (IER. NE. 0) GO TO 920

C        The defining section may be normalized.  Make sure it's not:

C        Determine what we've got:

         CALL BSAFGEOM (BSPLINE (1, KG), XLETEMP, YLETEMP, TLE,
     >                  XTEL, YTEL, XTEU, YTEU, CHORDTEMP,
     >                  THICKTEMP, XTHICKTEMP, CAMBER, XCAMBER,
     >                  AREATEMP, RADIUS (KG), -LUNCRT, IER)
         IF (IER .NE. 0) GO TO 930

C        Shift/scale to get what we want:

         CALL BSAFNORM (BSPLINE (1, KG), XLETEMP, YLETEMP,
     >                  CHORDTEMP, THICKTEMP, XLE (KG), YLE (KG),
     >                  CHORD (KG), THICKNESS (KG))

C        Perfectionist touch: do it again (a) for greater precision,
C        and (b) because defining stations involved in extrapolation
C        are processed once more often than those that are not.

         CALL BSAFGEOM (BSPLINE (1, KG), XLETEMP, YLETEMP, TLE,
     >                  XTEL, YTEL, XTEU, YTEU, CHORDTEMP,
     >                  THICKTEMP, XTHICKTEMP, CAMBER, XCAMBER,
     >                  AREATEMP, RADIUS (KG), -LUNCRT, IER)
         IF (IER .NE. 0) GO TO 930

         CALL BSAFNORM (BSPLINE (1, KG), XLETEMP, YLETEMP,
     >                  CHORDTEMP, THICKTEMP, XLE (KG), YLE (KG),
     >                  CHORD (KG), THICKNESS (KG))
      END DO


      KORDER = INT (BSPLINE (3, KG)) ! All curves are expected to
      NCPTS  = INT (BSPLINE (4, KG)) ! have the same values here
      IX1 = 6 + NCPTS + KORDER       ! Index of first control pt. X
      IY1 = IX1 + NCPTS              ! .......................... Y


C     Trap curves that don't have the same knots as required:

      KL = KG
      DO K = 1, NDEFINED - 1
         KG = KDEFINED (K)
         DO I = 6, IX1 - 1
            IF (BSPLINE (I, KG) .NE. BSPLINE (I, KL)) GO TO 960
         END DO
      END DO


C     Perform the spanwise lofting.  Modifying the control points
C     is equivalent to modifying the sections because the knots are
C     assumed to be the same for each section.

      DO K = 1, NDEFINED - 1
         KL = KDEFINED (K)
         KR = KDEFINED (K + 1)

C        Allow for extrapolation to the tips:

         IF (K .EQ. 1 .AND. KL .NE. 1) THEN
            K1 = 1
         ELSE
            K1 = KL + 1
         END IF

         IF (K + 1 .EQ. NDEFINED .AND. KR .NE. NSECTIONS) THEN
            K2 = NSECTIONS
         ELSE
            K2 = KR - 1
         END IF

C        Loft any missing stations linearly initially, and adjust later.
C        The X and Y ctl. pt. coords. are contiguous, so do both in one shot.

         WRITE (LUNCRT, 1010)

         DO KG = K1, K2

            WRITE (LUNCRT, '(A, I3)') '    Lofting station ', KG

            DO I = 1, IX1 - 1    ! Copy header + knots
               BSPLINE (I, KG) = BSPLINE (I, KL)
            END DO

            CALL LINTRP (2*NCPTS, BSPLINE (IX1, KL), BSPLINE (IX1, KR),
     >                   Z (KL), Z (KR), Z (KG), BSPLINE (IX1, KG))
         END DO

C        Impose the desired chord and maximum thickness.

         WRITE (LUNCRT, 1010)

         DO KG = K1, K2

            WRITE (LUNCRT, '(A, I3)') '    Retapering station ', KG

C           Determine what we've got:

            CALL BSAFGEOM (BSPLINE (1, KG), XLETEMP, YLETEMP, TLE,
     >                     XTEL, YTEL, XTEU, YTEU, CHORDTEMP,
     >                     THICKTEMP, XTHICKTEMP, CAMBER, XCAMBER,
     >                     AREATEMP, RADIUS (KG), -LUNCRT, IER)
            IF (IER .NE. 0) GO TO 930

C           Shift/scale to get what we want:

            CALL BSAFNORM (BSPLINE (1, KG), XLETEMP, YLETEMP,
     >                     CHORDTEMP, THICKTEMP, XLE (KG), YLE (KG),
     >                     CHORD (KG), THICKNESS (KG))
         END DO

      END DO


C     Verify the final chord, thickness, area, etc.:

      WRITE (LUNCRT, 1010)

      DO KG = 1, NSECTIONS

         WRITE (LUNCRT, '(A, I3)') '    Verifying station ', KG

         CALL BSAFGEOM (BSPLINE (1, KG), XLE (KG), YLE (KG), TLE,
     >                  XTEL, YTEL, XTEU, YTEU, CHORD (KG),
     >                  THICKNESS (KG), XTHICK (KG), CAMBER,
     >                  XCAMBER, AREA (KG), RADIUS (KG), -LUNCRT, IER)
         IF (IER .NE. 0) GO TO 930

      END DO

C     Tabulate the results for printing:

      WRITE (LUNPRINT, 1010) TITLEOUT (23:)
      WRITE (LUNPRINT, 1030)
      WRITE (LUNPRINT, 1040)
     >   (KG, Z (KG), XLE (KG), YLE (KG), CHORD (KG), THICKNESS (KG),
     >    XTHICK (KG), AREA (KG), KG = 1, NSECTIONS)

C     Save all sections in B-spline form:

      WRITE (LUNCRT, 1010)
     >   BLANK, ' Saving all sections as B-spline curves.'

      DO KG = 1, NSECTIONS

         WRITE (TITLEOUT (1:22), 1050) KG, Z (KG)  ! IGES limit is 72 chars.

         CALL DTRDWR ('W', LUNOUTSECTIONS, TITLEOUT, MAXC,
     >                BSPLINE (1, KG), NC, IER)
         IF (IER .NE. 0) GO TO 940

      END DO

C     ... and in discretized form for R22OPT, TranAir, etc.

      WRITE (LUNCRT, 1010)
     >   BLANK, ' Saving all sections in discretized form.'

      DO KG = 1, NSECTIONS

C        Determine reasonable discretization increments for the leading
C        and trailing edges.  This is empirical and sensitive: we don't
C        want too much bunching because the discretization is then too
C        coarse at mid-chord.  The square root of a function of the
C        circumference of the leading edge circle is used here, and bounded.
C        The calibration point is (PI * RHO = 0.1, DELTA ARC = 0.003).

         DELLE = SQRT (PI * RADIUS (KG) * 9.D-5)
         DELLE = MIN (MAX (DELLE, CHORD (KG) * 1.D-3),
     >                            CHORD (KG) * 4.D-3)
         DELTE = DELLE * 1.5D+0
         DELTE = MIN (MAX (DELTE, CHORD (KG) * 2.D-3),
     >                            CHORD (KG) * 4.D-3)

!!!         CALL BSAF2XY (BSPLINE (1, KG), NU, NL, -100.D+0 * DELLE,
!!!     >                 -100.D+0 * DELTE, T, X, Y, LUNCRT, IER)
         CALL BSAF2XY (BSPLINE (1, KG), NU, NL, -0.2D+0, -0.4D+0,
     >                 T, X, Y, LUNCRT, IER)
         IF (IER .NE. 0) GO TO 950

         WRITE (LUNCOORDINATES, 1060)
         WRITE (LUNCOORDINATES, 1070)
     >      KG, Z (KG), X (NL), Y (NL), CHORD (KG), ONE, ZERO, 1, 0
         WRITE (LUNCOORDINATES, 1080)
         WRITE (LUNCOORDINATES, 1090) ZERO, NU, NL
         WRITE (LUNCOORDINATES, 1100)
         WRITE (LUNCOORDINATES, 1110) (X (I), Y (I), I = NL, NXY)
         WRITE (LUNCOORDINATES, 1120)
         WRITE (LUNCOORDINATES, 1110) (X (I), Y (I), I = NL, 1, -1)

      END DO


      WRITE (LUNCRT, 1015) BLANK,
     >   'Tabular results file:       wingsections.summary',
     >   'All B-spline sections:      wingsections.bsplines',
     >   'Discretized sections:       wingsections.xy'

      GO TO 999


C     Error handling:
C     ---------------

  910 WRITE (LUNCRT, '(A)') ' *** Too many input stations.'
      GO TO 999

  920 WRITE (LUNCRT, 1020) 'BSCRVRD', IER
      GO TO 999

  930 WRITE (LUNCRT, 1020) 'BSAFGEOM', IER
      GO TO 999

  940 WRITE (LUNCRT, 1020) 'DTRDWR', IER
      GO TO 999

  950 WRITE (LUNCRT, 1020) 'BSAF2XY', IER
      GO TO 999

  960 WRITE (LUNCRT, 1025) KG
      GO TO 999

  999 STOP ' '

C     Formats:
C     --------

 1010 FORMAT (A)
 1015 FORMAT (1X, A)
 1020 FORMAT (' *** ', A, ' problem:  IER =', I4)
 1025 FORMAT (' *** The knots are mismatched for the last defining',
     >        ' section', /, ' *** and section ', I3, ': aborting.')
 1030 FORMAT ('! KG           Z       XLE       YLE     CHORD',
     >        '     MAX T X @ MAX T      AREA')
 1040 FORMAT (I4, F12.5, 5F10.5, F10.4)
 1050 FORMAT ('K: ', I3, '  Z:', F10.4, '  ')
 1060 FORMAT (' KG           Z         XLE         YLE       CHORD',
     >        '     THICK     TWIST NEW CRNK')
 1070 FORMAT (I3, 4F12.6, 2F10.6, 2I3)
 1080 FORMAT ('YSYM      NU        NL')
 1090 FORMAT (F2.0, 2I10)
 1100 FORMAT ('UPPER SURFACE')
 1110 FORMAT (2F11.6)
 1120 FORMAT ('LOWER SURFACE')

      END
