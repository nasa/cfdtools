C+------------------------------------------------------------------------------
C
      PROGRAM BSPROFILE
C
C     Purpose:
C     --------
C
C        BSPROFILE manipulates and/or displays a B-spline curve representing
C     an airfoil in wrap-around form.  It provides some of the functionality
C     of program PROFILE, which deals strictly with discretized (X,Y) airfoil
C     datasets but should still be of interest to a BSPROFILE user.
C
C        A given B-spline curve may be discretized for plotting purposes, or
C     for compatibility with other applications.  The plot includes the
C     control polygon and curve evaluations corresponding to the knots.
C
C        The airfoil is expected to consist of a single curve from the lower
C     trailing edge, round the leading edge (blunt or sharp), to the upper
C     trailing edge.  It need not be closed at the trailing edge.
C
C        This version can also discretize an arbitrary curve, which is
C     distinguished from the wrap-around airfoil case by whether the first
C     and last control point Xs (almost) match or not.
C
C        Results are in QPLOTable form.  Tailoring to some other plot package
C     may be required - possibly externally via a translator.  A modified
C     B-spline airfoil dataset may also be output.
C
C        BSPROFILE also provides for converting a discretized curve into a
C     B-spline curve, either by standard interpolation (where the number of
C     control points is determined by the number of data points), or by the
C     standard linear least squares method (to employ a specified, smaller,
C     number of control points), or by a nonlinear least squares method which
C     improves the curve fit attainable with a given number of control points
C     by iteratively adjusting the knots as well.  This more elaborate method
C     was developed at Ames to keep the number of control points to a minimum
C     for design-by-optimization purposes.
C
C        A curve fit may also be performed using a specified set of knots,
C     since multiple sections of a wing (say) with common knot vectors are
C     most convenient for a CAD system to convert to NURBS surface form.
C     Knot commonality also permits spanwise lofting of sections by applying
C     the lofting to the control points.
C
C        BSPROFILE is built upon the DT_NURBS library of NURBS utilities,
C     which is the only such package known to be in the public domain and
C     written in FORTRAN 77.  A collection of higher-level utilities, mostly
C     airfoil-related, has been developed in the Applied Aerodynamics Division
C     at Ames as the AA_NURBS library.  Some of these are used by BSPROFILE.
C     Along with a modest subset of the DT_NURBS package, they have also been
C     employed in CFD-based design-by-optimization applications with mixed
C     results.  Whether use of B-spline control points as design variables
C     is preferable to using conventional perturbing shape functions remains
C     an open question.  Either way, BSPROFILE enables an aerodynamicist to
C     interface with a CAD system in what is nowadays the preferred manner.
C     
C
C     Data format:
C     ------------
C
C        DT_NURBS employs a compact "C-array" format for 2-space B-spline
C     curves as follows:
C
C        C(1) = 1 (# parametric variables);
C        C(2) = 2 (for X and Y);
C        C(3) = polynomial order k (e.g., 4 for cubics);
C        C(4) = # control points, N;
C        C(5) = pointer for efficient evaluations;
C        C(6 : 5 + N + k) = knot vector,
C        followed by the X coords. of N control pts.,
C        followed by the Y coords. of N control pts.
C
C        Initially, BSPROFILE reads and writes B-spline curves in this format,
C     which is far more readable and compact than an IGES representation and
C     may well be preferred for a design application.  Direct I/O of IGES
C     representations may become an option, but for now, program DTIGES serves
C     to convert a curve from DT_NURBS to IGES format and vice versa.
C
C        The curve-fitting options expect (X,Y) coordinates forming a single
C     curve as follows:
C
C        Title  ! Descriptive text
C        N      ! # points
C        X   Y  ! Coordinates
C        :   :
C
C        If a knot vector is specified for the fit, it is extracted from a
C     B-spline curve file in the DT_NURBS format.  A file of values of the
C     parametric variable U corresponding to the data points is also required
C     as input.  It should be in a single-column version of the (X,Y) data
C     format, as output from the original fitting process which produced the
C     knot vector being reused.
C
C
C     Environment:
C     ------------
C
C        DEC/OpenVMS, SGI/IRIX
C        64-bit FORTRAN 77 with these Fortran 90 extensions:
C
C           ADVANCE='NO' in place of $ edit descriptors
C           IMPLICIT NONE
C           Names exceeding 6 characters
C           ! comments
C
C
C     History:
C     --------
C
C     05/19/92  D.A.Saunders  Adapted from TESTAFGEOM.
C     05/29/92    "    "      Added area calculation.
C     06/01/92    "    "      Added option to normalize the chord.
C     09/19/92    "    "      Made the disretization variable, for possible
C                             use of coordinates in a flow solver.
C     12/21/92    "    "      Provided for discretizing an arbitrary 2-space
C                             curve (as needed for a camber line).
C     02/26/93    "    "      Provided for adjusting chord, leading edge,
C                             and thickness of an airfoil curve.
C     03/31/93    "    "      '$' carriage control failed on an IRIS - use
C                             READS for the input file name.
C     04/05/93    "    "      Provided for adjusting the trailing edge with
C                             a shear to compensate for a moved leading edge.
C     04/26/93    "    "      Provided for discretization according to
C                             arc-length.  Introduced display of the leading
C                             edge radius of curvature (also used to control
C                             the discretization).
C     06/23/93    "    "      Radius of curvature is units-dependent: needed
C                             chord in the empirical formula for unnormalized
C                             airfoils.
C     01/28/94    "    "      Added knot insertion option (mainly to fix a
C                             control point at an airfoil leading edge).
C     02/11/94    "    "      Provided for standard linear transformations.
C     July '97    "    "      Version 2.0: Installed B-spline curve fitting
C                             options originally implemented in a version of
C                             SMOOTH; provided for multiple grid distributions.
C     Apr 2004    "    "      Handle up to 256 points per curve; TOL = 0.01 was
C                             too small to identify space shuttle wing sections
C                             as airfoils - make it 0.02.
C     08/13/13    "    "      Resurrecting BSPROFILE on an Intel x86 64-bit
C                             processor gave airfoil area trouble, though the
C                             nonlinear B-spline curve fit of a 12% oblique wing
C                             airfoil appears to work as earlier (degree 4, 16
C                             control points, 10,000 iterations, 500 points per
C                             surface evaluation of the fit).  The machine-
C                             dependent constants are set for an IRIS 4D though,
C                             so maybe that's the reason (no time to pursue it).
C
C     Author:  David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Procedures:
C     -----------

      EXTERNAL
     >   BSAF2XY,   ! Discretizes a wrap-around airfoil-type B-spline curve
     >   BSAFGEOM,  ! Determines geometric properties of "   "   "   "   "
     >   BSAFNORM,  ! Adjusts chord, leading edge, and thickness of an airfoil
     >   BSARCMAP,  ! Finds the T distribn. corresp. to a reqd. arc discretizn.
     >   BSCRVMOD,  ! Modifies a B-spline curve in-place (scale, shift, rotate)
     >   BSFIT,     ! Nonlinear least squares B-spline curve fit
     >   BSINSERT,  ! Inserts a knot with specified multiplicity
     >   BSZEROMN,  ! Finds a zero, minimum, or maximum of a B-spline curve
     >   CHORDS2D,  ! Cumulative chord length utility
     >   DSTRIB,    ! Basic 1-D grid distributions
     >   DTNSI,     ! Interpolatory B-spline curve fit
     >   DTRDWR,    ! Reads/writes one NURBS curve in DT_NURBS format
     >   DTSPVL,    ! Evaluates a NURBS curve at a specified T
     >   FOILGRD,   ! Linear + quadratic + sine + cosine 1-D grid utility
     >   GETDIS,    ! 1-D grid distribution prompting utility
     >   HTDIS4,    ! Vinokur distribution utility
     >   READER,    ! Prompting utility (entry points READI, READS, etc.)
     >   TOGGLE     ! Enable/disable utility

C-------------------------------------------------------------------------------

C     Constants:
C     ----------

      INTEGER
     >   LUNCRT, LUNDISK, LUNKBD, LUNLOG, LUNPLT,
     >   MAXC, MAXCPTS, MAXFIT, MAXGRID, MAXORDER, MAXPLT, MAXWORK,
     >   NFIT, NMENU, NMOD
      REAL
     >   BIG, HALF, ONE, PERCENT, TOL, ZERO
      CHARACTER
     >   PFILE * 14, SFILE * 17, UFILE * 17
      PARAMETER
     >  (PFILE   = 'bsprofile.plot',
     >   SFILE   = 'bsprofile.bspline',
     >   UFILE   = 'bsprofile.uvalues',
     >   LUNCRT  = 6,     ! Screen
     >   LUNDISK = 1,     ! Various in/out disk files
     >   LUNKBD  = 5,     ! Keyboard
     >   LUNLOG  = 2,     ! Printable output
     >   LUNPLT  = 3,     ! Plottable results
     >   MAXCPTS = 500,   ! Maximum # control points handled
     >   MAXFIT  = 1000,  ! Maximum # (X,Y) points able to be curve-fitted
     >   MAXGRID = 500,   ! Maximum # grid points imposable on a curve <= MAXPLT
     >   MAXORDER= 11,    ! Maximum value of (degree + 1) handled
     >   MAXC    = 5 + MAXORDER + 3 * MAXCPTS,
     >   MAXPLT  = 1000,  ! Maximum # discretized points plottable >= MAXGRID
     >   MAXWORK = MAXCPTS*(MAXORDER + 1) + 2*MAXORDER,  ! DTLSA requirement
     >   NFIT    = 5,     ! # curve-fitting options
     >   NMENU   = 5,     ! # main options
     >   NMOD    = 8,     ! # shift/scale/rotate transformations
     >   BIG     = 1.E+10,
     >   HALF    = 0.5,
     >   ONE     = 1.0,
     >   PERCENT = 0.01,
     >   TOL     = 0.02,  ! Fraction of data range used to identify an airfoil
     >   ZERO    = 0.0)

C     Variables:
C     ----------

      INTEGER
     >   I, I1, I2, ICC, IEND, IER, IMENU, IOSTAT, ISTART, IX, IY,
     >   LASTCPT, LUNERR, MODE, MULTI, NC, NCPTS, NCTL, NDEG, NEVAL,
     >   NIT, NK, NKNOTS, NL, NP, NU, NXY
      REAL
     >   ARC (MAXGRID), BSU(MAXFIT), C (MAXC), P (4),
     >   T (MAXPLT), TNORM (MAXGRID), X (MAXPLT), XFIT (MAXFIT),
     >   XYFIT (MAXFIT, 2), Y (MAXPLT), YFIT (MAXFIT), XY (2),
     >   WORK (MAXWORK),
     >   AREA, CAMBER, CHORD, CHORDNEW, CHORDOLD, DELTAT,
     >   DSQ, DUMMY, RADIUS, SLOPE, TEVAL, THICK, THICKNEW, TINSERT,
     >   TLE, TTOTAL, XCAMBER, XI, XLE, XLENEW, XMAX, XMIN, XTEL, XTEU,
     >   XTHICK, YLE, YLENEW, YMAX, YMIN, YTEL, YTELOLD, YTEU, YTEUOLD
      LOGICAL
     >   ADJUST, CR, ENABLED (NMENU), EOF, FIT, FIXLE, FIXTE, GRID,
     >   INSERT, INTERPOLATE, MODIFY, NEWSPLINE, NONLINEAR, PLOT,
     >   RESTART, SAMEKNOTS, WRAP, YES
      CHARACTER
     >   FILENAME * 32, DATASET * 32, FITMENU (NFIT) * 58,
     >   MENU (NMENU) * 57, MODMENU (NMOD) * 28, TITLE * 80

C     Storage:
C     --------

      EQUIVALENCE
     >  (XFIT, XYFIT (1, 1)), (YFIT, XYFIT (1, 2))  ! To suit DTNSI & BSFIT

      DATA MENU
     >  /'1:  Display a B-spline curve and its various properties. ',
     >   '2:  Scale, shift, and/or rotate a B-spline curve.        ',
     >   '3:  Add a [multiple?] knot (e.g. to fix a pt. on curve). ',
     >   '4:  Impose a mesh on a B-spline curve (several options). ',
     >   '5:  Fit a B-spline curve to (X,Y) pts. (several options).'/

      DATA ENABLED
     >  /.TRUE.,
     >   .FALSE.,
     >   .FALSE.,
     >   .FALSE.,
     >   .FALSE./

      DATA FITMENU
     >  /'1:  Interpolation?                                        ',
     >   '2:  Linear least squares (control points only)?           ',
     >   '3:  Nonlinear least squares (knots optimized as well)?    ',
     >   '4:  Restart a partially converged nonlinear iteration?    ',
     >   '5:  Fit a B-spline curve using knots from a related curve?'/

      DATA MODMENU
     >  /'1:  Scale X?',
     >   '2:  Scale Y?',
     >   '3:  Scale X and Y similarly?',
     >   '4:  Shift X?',
     >   '5:  Shift Y?',
     >   '6:  Scale X + shift?',
     >   '7:  Scale Y + shift?',
     >   '8:  Rotate curve?'/

      DATA NEWSPLINE /.FALSE./

C     Execution:
C     ----------

      WRITE (LUNCRT, 1009) 'Welcome to BSPROFILE Version 2.0, 07/25/97.'

C     Determine the main option(s):

      CALL TOGGLE (NMENU, MENU, ENABLED, LUNCRT, LUNKBD)

      PLOT   = ENABLED (1)
      MODIFY = ENABLED (2)
      INSERT = ENABLED (3)
      GRID   = ENABLED (4)
      FIT    = ENABLED (5)

      IF (FIT) THEN
         IF (.NOT. GRID) PLOT = .TRUE.  ! Ensure some display of the result
      END IF

      IF (.NOT. PLOT) NL = 0  ! Else enter NL >= 0 below

C     Open a plottable output file regardless because most results go there:

      OPEN (UNIT=LUNPLT, FILE=PFILE, STATUS='UNKNOWN')


C     Deal with the curve-fitting options first,
C     because further options then apply to the result.
C     -------------------------------------------------

      IF (FIT) THEN

C        Determine curve-fit options:

         WRITE (LUNCRT, 1009) 'Curve-fit options:'
         WRITE (LUNCRT, 1011) FITMENU, ' '
   20    IMENU = 3
         CALL READI (LUNCRT, 'Option? (<CR> = 3): ',
     >               LUNKBD, IMENU, CR, EOF)
         IF (EOF) GO TO 999
         IF (IMENU .LT. 1 .OR. IMENU .GT. NFIT) GO TO 20

         NEWSPLINE   = .TRUE.
         INTERPOLATE = IMENU .EQ. 1  ! Redefined if # ctl. pts. = # data pts.
         NONLINEAR   = IMENU .EQ. 3
         RESTART     = IMENU .EQ. 4
         SAMEKNOTS   = IMENU .EQ. 5

         IF (RESTART .OR. SAMEKNOTS) THEN  ! Further input files are needed

   30       IF (RESTART) THEN
               CALL READS (LUNCRT,
     >            'B-spline file to restart iteration from? ',
     >            LUNKBD, DATASET, CR, EOF)
            ELSE
               CALL READS (LUNCRT,
     >            'B-spline file with knots to be reused? ',
     >            LUNKBD, DATASET, CR, EOF)
            END IF

            IF (EOF) GO TO 999
            IF (CR)  GO TO 30

            OPEN (UNIT=LUNDISK, FILE=DATASET, STATUS='OLD', ERR=30)

            CALL DTRDWR ('R', LUNDISK, TITLE, MAXC, C, NC, IER)
            IF (IER .NE. 0) THEN
               WRITE (LUNCRT, 1010) 'DTRDWR (reusing)', IER
               GO TO 999
            END IF

            CLOSE (UNIT=LUNDISK)

            NDEG   = NINT (C (3)) - 1
            NCPTS  = NINT (C (4))
            NKNOTS = NCPTS - (NDEG + 1)  ! # interior knots
            NCTL   = NCPTS               ! # requested = previous #

   40       CALL READS (LUNCRT,
     >         'Corresponding file of U values <-> data points: ',
     >          LUNKBD, DATASET, CR, EOF)
            IF (EOF) GO TO 999
            IF (CR)  GO TO 40

            OPEN (UNIT=LUNDISK, FILE=DATASET, STATUS='OLD', ERR=40)

            READ (LUNDISK, *)  ! Skip title
            READ (LUNDISK, *, IOSTAT=IOSTAT) NU
            IF (IOSTAT .NE. 0) THEN
               WRITE (LUNCRT, '(/, A, I6)')
     >            ' *** Bad number of U values. IOSTAT:', IOSTAT
               GO TO 999
            END IF
            NU = MIN (NU, MAXFIT)  ! Any discrepancy is detected below

            READ (LUNDISK, *, IOSTAT=IOSTAT) (BSU (I), I = 1, NU)
            IF (IOSTAT .NE. 0) THEN
               WRITE (LUNCRT, '(/, A, I6)')
     >            ' *** Error reading U values. IOSTAT:', IOSTAT
               GO TO 999
            END IF

            CLOSE (UNIT=LUNDISK)
         END IF

C        Open and read the file of (X,Y) points:

  100    CALL READS (LUNCRT, '(X,Y) curve file? ', LUNKBD, DATASET,
     >               CR, EOF)
         IF (EOF) GO TO 999
         IF (CR)  GO TO 100

         OPEN (UNIT=LUNDISK, FILE=DATASET, STATUS='OLD', ERR=100)

         READ (LUNDISK, '(A)') TITLE
         READ (LUNDISK, *, IOSTAT=IOSTAT) NXY
         IF (IOSTAT .NE. 0) THEN
            WRITE (LUNCRT, '(/, A, A)') ' *** Error reading # (X,Y)',
     >         ' points following title - aborting.'
            GO TO 999
         END IF

         IF (NXY .LT. 2 .OR. NXY .GT. MAXFIT) THEN
            WRITE (LUNCRT, '(/, 2(A, I10))')
     >         ' *** Bad # (X,Y) points:', NXY, '     Limit:', MAXFIT
            GO TO 999
         END IF

         IF (RESTART .OR. SAMEKNOTS) THEN
            IF (NU .NE. NXY) THEN
               WRITE (LUNCRT, '(/, (A, I10))')
     >            ' *** # restart U values: ', NU,
     >            '     # (X,Y) data points:', NXY,
     >            '     Aborting.'
               GO TO 999
            END IF
         END IF

C        Read (X,Y)s individually to allow extra columns or comments:

         DO I = 1, NXY
            READ (LUNDISK, *, IOSTAT=IOSTAT) XFIT (I), YFIT (I)
            IF (IOSTAT .NE. 0) THEN
               WRITE (LUNCRT, '(/, 2(A, I6))')
     >            ' *** Error reading (X,Y) for point #', I,
     >            '.  IOSTAT:', IOSTAT
               GO TO 999
            END IF
         END DO

         CLOSE (UNIT=LUNDISK)
         WRITE (LUNCRT, '(/, A, I5)') ' # (X,Y)s read:', NXY

         IF (SAMEKNOTS) THEN

            NIT    = 0       ! Else the knots will change
            ISTART = 1       ! BSFIT restart flag

         ELSE IF (RESTART) THEN

            ISTART = 1

         ELSE

            ISTART = 0
            NDEG   = 4

            CALL READI (LUNCRT, 'Degree? (1 to 10; <CR>=4): ',
     >                  LUNKBD, NDEG, CR, EOF)
            IF (EOF) GO TO 999

            IF (INTERPOLATE) THEN

               NCPTS = NXY + 2 * (NDEG / 2)  ! Because of knot choice -
                                             ! see DTCGEN
            ELSE  ! IMENU = 2 or 3

  110          NCTL = 0
               CALL READI (LUNCRT,
     >'# control pts.? (Degree+1:NXY-1=least sqrs.; NXY=interpolate): ',
     >            LUNKBD, NCTL, CR, EOF)
               IF (EOF) GO TO 999

C              NXY here forces use of BSFIT, not DTNSI, presumably for
C              comparison purposes.

               IF (NCTL .LT. NDEG + 1 .OR. NCTL .GT. NXY) GO TO 110
               NCPTS = NCTL

            END IF

         END IF


         IF (NONLINEAR .OR. RESTART) THEN

  120       NIT = 2000

            IF (NONLINEAR) THEN
               CALL READI (LUNCRT,
     >     '# best-knot itns.? (0 = linear least squares; <CR>=2000): ',
     >                     LUNKBD, NIT, CR, EOF)
            ELSE
               CALL READI (LUNCRT,
     >                  '# further best-knot iterations? (<CR>=2000): ',
     >                     LUNKBD, NIT, CR, EOF)
            END IF

            IF (EOF) GO TO 999
            IF (NIT .LT. 0) GO TO 120

            LUNERR = LUNCRT
            IF (NIT .GT. 0) THEN
               YES = .FALSE.
               CALL READY (LUNCRT,
     >       'Show the iterations? (Y/N; No=just the last; <CR>=No): ',
     >                     LUNKBD, YES, CR, EOF)
               IF (.NOT. YES) THEN
                  LUNERR = -LUNCRT
               END IF
            ELSE  ! Equivalent of IMENU = 2 (DTLSA)
            END IF

         END IF


C        Perform the curve fit.
C        ----------------------

         IF (INTERPOLATE) THEN  ! "Natural spline" interpolation method

C           Use normalized chord lengths for the interpolation nodes.
C           DTNSI sets the knots to be these nodes if NDEG is odd,
C           else the knots are half way between the nodes if NDEG is even.

            CALL CHORDS2D (NXY, XFIT, YFIT, .TRUE., TTOTAL, BSU)

            ICC = -1    ! Initial/continue code; -1 = first dependent variable;
                        ! DTNSI returns it as |ICC|, or 0 for an error

            DO I = 1, 2 ! X and Y

               CALL DTNSI (NXY, BSU, XYFIT (1, I), NDEG, ICC,
     >                     WORK, MAXWORK, MAXC, C, NC, IER)
               IF (IER. LT. 0) THEN
                  WRITE (LUNCRT, 1010) 'DTNSI', IER
                  GO TO 999
               END IF

               ICC = ICC + 1  ! So as to APPEND control pt. coords. in C (*)
            END DO

         ELSE  ! Linear (NIT = 0) or nonlinear least squares methods

            YES = .TRUE.
            CALL READY (LUNCRT,
     >         'Constrain the end points? (Y/N; <CR>=yes): ',
     >         LUNKBD, YES, CR, EOF)
            IF (EOF) GO TO 999

            IF (YES) THEN
               IEND = 1
            ELSE
               IEND = 0
            END IF
 
            IF (LUNERR .LT. 0 .AND. NIT .GT. 5)
     >         WRITE (LUNCRT, 1011) 'Working ...'

            CALL BSFIT (NXY, BSU, XYFIT, MAXFIT, 2, NDEG, 0, DUMMY,
     >                  ISTART, IEND, NIT, NCPTS, MAXWORK, WORK,
     >                  MAXC, C, NC, LUNERR, IER)

C           DTLSA/BSFIT give diagnostics; 14 = non-monotonic Us

            IF (IER. LT. 0 .OR. IER .GE. 14) THEN
               WRITE (LUNCRT, 1010) 'BSFIT', IER
               GO TO 999
            END IF

            IF (NCPTS .EQ. NXY) INTERPOLATE = .TRUE.  ! For plotting purposes

         END IF

C        Save [optimized?] U (or T) values for possible reuse?

         IF (.NOT. SAMEKNOTS) THEN

            OPEN (UNIT=LUNDISK, FILE=UFILE, STATUS='UNKNOWN')

            WRITE (LUNDISK, 1001) 'U values from B-spline curve fit'
            WRITE (LUNDISK, '(I4)') NXY
            WRITE (LUNDISK, '(E23.16)') (BSU (I), I = 1, NXY)

            CLOSE (UNIT=LUNDISK)

         END IF

         FILENAME = SFILE  ! Used as for no-fit case below


      ELSE  ! No curve-fitting involved

C        Prompt for an input B-spline curve file and read it:
C        ----------------------------------------------------

  150    CALL READS (LUNCRT, 'B-spline curve file? ',
     >               LUNKBD, FILENAME, CR, EOF)
         IF (EOF) GO TO 999
         IF (CR)  GO TO 150

         OPEN (UNIT=LUNDISK, FILE=FILENAME, STATUS='OLD', ERR=150)

         CALL DTRDWR ('R', LUNDISK, TITLE, MAXC, C, NC, IER)
         IF (IER .NE. 0) THEN
            WRITE (LUNCRT, 1010) 'DTRDWR (reading)', IER
            GO TO 999
         END IF

         CLOSE (UNIT=LUNDISK)

      END IF


C     Most of the following applies whether there was a curve fit or not.
C     -------------------------------------------------------------------

      NDEG  = NINT (C (3)) - 1
      NCPTS = NINT (C (4))

      IX = 6 + NDEG + NCPTS   ! Offset of X coordinates of control pts.
      IY = IX + NCPTS         ! ......... Y ...........................

      IF (.NOT. FIT) WRITE (LUNCRT, 1006) 'Degree: ', NDEG,
     >   '    Control points: ', NCPTS

C     Determine the (control point) data range in order to distinguish
C     between a probable airfoil and some arbitrary curve:

      XMIN =  BIG
      XMAX = -BIG
      YMIN =  BIG
      YMAX = -BIG

      DO I = 1, NCPTS
         XMIN = MIN (XMIN, C (I + IX))
         XMAX = MAX (XMAX, C (I + IX))
         YMIN = MIN (YMIN, C (I + IY))
         YMAX = MAX (YMAX, C (I + IY))
      END DO

C     Compare first and last control points:

      DSQ  = (C (1 + IX) - C (NCPTS + IX)) ** 2 +
     >       (C (1 + IY) - C (NCPTS + IY)) ** 2
      WRAP = SQRT (DSQ) .LT. TOL * MAX (XMAX - XMIN, YMAX - YMIN)

      WRITE (LUNPLT, 1001) TITLE, '""', 'X', 'Y', '""'


C     Transform the curve?
C     --------------------

      IF (MODIFY) THEN

         WRITE (LUNCRT, 1009) 'Shift/scale/rotate options:'
         WRITE (LUNCRT, 1011) MODMENU, ' '
  200    IMENU = 0

         CALL READI (LUNCRT,
     >      'Option? (<CR> = done): ', LUNKBD, IMENU, CR, EOF)

         IF (IMENU .LT. 0 .OR. IMENU .GT. NMOD) GO TO 200
         IF (IMENU .EQ. 0) THEN    ! Drop through
            WRITE (LUNCRT, '(A)')  ! SGI seems to leave the cursor hanging
         ELSE    ! EOF misbehaves on SGI - avoid using it
            NEWSPLINE = .TRUE.
         END IF

         IF (IMENU .EQ. 1) THEN       ! Scale X

            WRITE (LUNCRT, 1007, ADVANCE='NO') 'X scale factor: '
            READ  (LUNKBD, *) P (1)
            CALL BSCRVMOD ('SCALEX', P, C)

         ELSE IF (IMENU .EQ. 2) THEN  ! Scale Y

            WRITE (LUNCRT, 1007, ADVANCE='NO') 'Y scale factor: '
            READ  (LUNKBD, *) P (1)
            CALL BSCRVMOD ('SCALEY', P, C)

         ELSE IF (IMENU .EQ. 3) THEN  ! Scale X and Y similarly?

            WRITE (LUNCRT, 1007, ADVANCE='NO') 'Scale factor: '
            READ  (LUNKBD, *) P (1)
            CALL BSCRVMOD ('SCALEXY', P, C)

         ELSE IF (IMENU .EQ. 4) THEN  ! Shift X

            WRITE (LUNCRT, 1007, ADVANCE='NO') 'X shift: '
            READ  (LUNKBD, *) P (1)
            CALL BSCRVMOD ('SHIFTX', P, C)

         ELSE IF (IMENU .EQ. 5) THEN  ! Shift Y

            WRITE (LUNCRT, 1007, ADVANCE='NO') 'Y shift: '
            READ  (LUNKBD, *) P (1)
            CALL BSCRVMOD ('SHIFTY', P, C)

         ELSE IF (IMENU .EQ. 6) THEN  ! Scale X + shift

            WRITE (LUNCRT, 1007, ADVANCE='NO') 'X scale and shift: '
            READ  (LUNKBD, *) P (2), P (1)
            CALL BSCRVMOD ('LINEARX', P, C)

         ELSE IF (IMENU .EQ. 7) THEN  ! Scale Y + shift

            WRITE (LUNCRT, 1007, ADVANCE='NO') 'Y scale and shift: '
            READ  (LUNKBD, *) P (2), P (1)
            CALL BSCRVMOD ('LINEARY', P, C)

         ELSE IF (IMENU .EQ. 8) THEN  ! Rotate curve

            WRITE (LUNCRT, 1007, ADVANCE='NO')
     >         'Angle?  (Degrees; +ve is anticlockwise): '
            READ  (LUNKBD, *) P (1)
            WRITE (LUNCRT, 1007, ADVANCE='NO')
     >         'Center of rotation, (X, Y): '
            READ  (LUNKBD, *) P (2), P (3)
            CALL BSCRVMOD ('ROTATE', P, C)

         END IF

         IF (IMENU .NE. 0) GO TO 200

      END IF


C     Airfoil curve options:
C     ----------------------

      IF (WRAP) THEN

C        Knot insertion:
C        ---------------

         IF (INSERT) THEN

            NEWSPLINE = .TRUE.

C           Default to fixing a control point at the leading edge.
C           Find the appropriate value of the parametric variable:

            CALL BSZEROMN ('MIN', 1, C, C (6), C (IX), TINSERT, XY,
     >                     LUNCRT, IER)

            IF (IER .NE. 0) THEN
               WRITE (LUNCRT, 1010) 'BSZEROMN', IER
               GO TO 999
            END IF

            CALL READD (LUNCRT,
     >         'Knot value to insert?  (<CR> = Leading edge value): ',
     >         LUNKBD, TINSERT, FIXLE, EOF)
            IF (EOF) GO TO 999

            MULTI = NDEG
            CALL READI (LUNCRT, 
     >         'Multiplicity?  (<CR> = Degree): ', LUNKBD, MULTI,
     >         CR, EOF)
            IF (EOF) GO TO 999 
            MULTI = MAX (1, MIN (MULTI, NDEG))
            FIXLE = FIXLE .AND. MULTI .EQ. NDEG

            CALL BSINSERT (MAXC, C, TINSERT, MULTI, LASTCPT, C, IER)

            IF (IER .NE. 0) THEN
               WRITE (LUNCRT, 1010) 'BSINSERT', IER
               GO TO 999
            END IF

            IF (LASTCPT .GT. 0) THEN
               NCPTS = C (4)
               IX = 6 + NDEG + NCPTS
               IY = IX + NCPTS

               IF (FIXLE) THEN
                  WRITE (LUNCRT, 1006)
     >               'Leading edge control point number: ', LASTCPT
               ELSE
                  WRITE (LUNCRT, 1006)
     >               'Last control point number inserted: ', LASTCPT
               END IF
            ELSE
               NEWSPLINE = .FALSE.
               WRITE (LUNCRT, 1006)
     >            'Knot insertion suppressed as invalid.'
            END IF

         END IF


C        Discretization for plot purposes:
C        ---------------------------------

         IF (PLOT) THEN
  250       NL = MAXPLT / 2
            CALL READI (LUNCRT,
     >         '# points per surface to plot? (500 max.; <CR>=200): ',
     >         LUNKBD, NL, CR, EOF)
            IF (EOF) GO TO 999

            IF (NL .GT. MAXPLT / 2) THEN
               WRITE (LUNCRT, '(A, I3)') ' Present limit: ', MAXPLT / 2
               GO TO 250
            END IF
         END IF

         NU = NL

C        Determine the standard geometric properties:

         CALL BSAFGEOM (C, XLE, YLE, TLE, XTEL, YTEL, XTEU, YTEU,
     >                  CHORD, THICK, XTHICK, CAMBER, XCAMBER, AREA,
     >                  RADIUS, LUNCRT, IER)

         IF (IER .NE. 0) THEN
            WRITE (LUNCRT, 1010) 'BSAFGEOM (1)', IER
CCC         GO TO 999
         END IF


C        Adjustment of chord and/or thickness:
C        -------------------------------------

         ADJUST = .TRUE.
         CALL READY (LUNCRT,
     >      'Adjust chord and/or maximum thickness?  [Yes] ',
     >      LUNKBD, ADJUST, CR, EOF)

         IF (ADJUST) THEN

            NEWSPLINE = .TRUE.

            CHORDOLD = CHORD
            YTELOLD = YTEL
            YTEUOLD = YTEU

            CHORDNEW = ONE
            CALL READD (LUNCRT, 'Desired chord?   [1.] ',
     >                  LUNKBD, CHORDNEW, CR, EOF)

            XLENEW = ZERO
            CALL READD (LUNCRT, 'Leading edge X?  [0.] ',
     >                  LUNKBD, XLENEW, CR, EOF)

            YLENEW = ZERO
            CALL READD (LUNCRT, 'Leading edge Y?  [0.] ',
     >                  LUNKBD, YLENEW, CR, EOF)

  300       THICKNEW = THICK
            CALL READD (LUNCRT, 'Maximum thickness?  [As is (NOT %)] ',
     >                  LUNKBD, THICKNEW, CR, EOF)
            IF (THICKNEW / CHORDNEW .GT. ONE) GO TO 300

            CALL BSAFNORM (C, XLE, YLE, CHORD, THICK,
     >                     XLENEW, YLENEW, CHORDNEW, THICKNEW)

            CALL BSAFGEOM (C, XLE, YLE, TLE, XTEL, YTEL, XTEU, YTEU,
     >                     CHORD, THICK, XTHICK, CAMBER, XCAMBER, AREA,
     >                     RADIUS, LUNCRT, IER)

            IF (IER .NE. 0) THEN
               WRITE (LUNCRT, 1010) 'BSAFGEOM (2)', IER
CCC            GO TO 999
            END IF

            IF (ABS (CHORD - CHORDOLD) .LT. CHORDOLD * 0.05) THEN
               FIXTE = .TRUE.
               CALL READY (LUNCRT,
     >            'Restore original trailing edge point(s)?  [Yes] ',
     >            LUNKBD, FIXTE, CR, EOF)
            ELSE
               FIXTE = .FALSE.
            END IF

            IF (FIXTE) THEN

C              The trailing edge has presumably moved slightly as a result of
C              adjusting to the B-spline's true leading edge.  The original
C              trailing edge control points need to be restored in some way
C              that implies as minimal a shape distortion as is practical.
C              Even a graduated vertical shear of the control points is
C              awkward.  (Rotating about the leading edge then realigning
C              is probably more awkward yet.)

C              (Almost) avoid changing the maximum thickness by leaving alone
C              any points to the left of the corresponding abscissa.

C              Lower surface:

               SLOPE = (YTELOLD - YTEL) / (XTEL - XTHICK)

               DO I = 1, NCPTS / 2
                  XI = C (IX + I)
                  IF (XI .GT. XTHICK) C (IY + I) =
     >                                C (IY + I) + (XI - XTHICK) * SLOPE
               END DO
                     
C              Upper surface:

               SLOPE = (YTEUOLD - YTEU) / (XTEU - XTHICK)

               DO I = NCPTS / 2, NCPTS
                  XI = C (IX + I)
                  IF (XI .GT. XTHICK) C (IY + I) =
     >                                C (IY + I) + (XI - XTHICK) * SLOPE
               END DO

               CALL BSAFGEOM (C, XLE, YLE, TLE, XTEL, YTEL, XTEU, YTEU,
     >                        CHORD, THICK, XTHICK, CAMBER, XCAMBER,
     >                        AREA, RADIUS, LUNCRT, IER)

               IF (IER .NE. 0) THEN
                  WRITE (LUNCRT, 1010) 'BSAFGEOM (3)', IER
CCC               GO TO 999
               END IF

            END IF

         END IF

C        The following attempts to align the numbers in the PostScript plot:

         WRITE (LUNPLT, 1005)
     >      'Leading edge X:          ', XLE,
     >      'Leading edge Y:          ', YLE,
     >      'Leading edge radius:   ', RADIUS,
     >      'Lower trailing edge X:', XTEL,
     >      'Lower trailing edge Y:', YTEL,
     >      'Upper trailing edge X: ', XTEU,
     >      'Upper trailing edge Y: ', YTEU,
     >      'Chord:                          ', CHORD,
     >      'Maximum thickness:   ', THICK,
     >      'Corresponding X:        ', XTHICK,
     >      'Maximum camber:      ', CAMBER,
     >      'Corresponding X:        ', XCAMBER,
     >      'Area:                            ', AREA

         WRITE (LUNPLT, 1001) '""'
         WRITE (LUNPLT, 1012) 'Degree: ', NDEG,
     >      '    Control points: ', NCPTS
         WRITE (LUNPLT, 1001)
         WRITE (LUNPLT, 1001)
     >      ' $OPTIONS PLOT=''SCALE'', WIDTH=6., $END'


C        Show curve-fit data?

         IF (FIT) THEN

            WRITE (LUNPLT, 1001) ' $OPTIONS'
            WRITE (LUNPLT, 1002) ' LEGEND=', 'Curve-fit data points',
     >         ' LINE=', 'SYMBOLS'

            WRITE (LUNPLT, 1001) ' SYMBOL=3', ' $END'
            WRITE (LUNPLT, 1003) (XFIT (I), YFIT (I), I = 1, NXY)

         END IF

         IF (NL .GT. 0) THEN

C           Discretize the B-spline curve for plotting purposes.
C           BSAF2XY is not as precise in terms of arc-length increments
C           as use of BSARCMAP would be, but it's adequate for airfoils.

            CALL BSAF2XY (C, NU, NL, T, X, Y, LUNCRT, IER)

            IF (IER .NE. 0) THEN
               WRITE (LUNCRT, 1010) 'BSAF2XY', IER
               GO TO 999
            END IF

            WRITE (LUNPLT, 1001) ' $OPTIONS LINE=''SOLID'' $END'
            WRITE (LUNPLT, 1003) (X (I), Y (I), I = 1, NL + NU - 1)
            WRITE (LUNPLT, 1004)

         END IF


C        Precise grid along the airfoil surfaces:
C        ----------------------------------------

  400    IF (GRID) THEN

C           Get the surface distribution details:

            NL = MAXPLT / 2  ! Max. is input to GETDIS
            NP = 4

            CALL GETDIS (LUNCRT, LUNKBD, MODE, NP, P, NL)

            IF (MODE .EQ. -99) THEN  ! Quit
               GRID = .FALSE.
               GO TO 400
            END IF

            NU = NL

C           Generate a normalized distribution on [0, 1] for LE to TE:

            IF (MODE .LE. 3) THEN  ! Uniform or sinusoidal

               CALL DSTRIB (MODE, NP, P, NU, ZERO, ONE, TNORM)

            ELSE IF (MODE .EQ. 4) THEN  ! Vinokur

               IF (P (1) .LT. ZERO) P (1) = -PERCENT * P (1)
               IF (P (2) .LT. ZERO) P (2) = -PERCENT * P (2)

               CALL HTDIS4 (.TRUE., ZERO, ONE, P (1), P (2), NU, TNORM,
     >                      -LUNCRT, IER)

               IF (IER .NE. 0) THEN
                  WRITE (LUNCRT, 1010) 'HTDIS4 on [0, 1]', IER
                  GO TO 999
               END IF

            ELSE IF (MODE .EQ. 5) THEN  ! Linear + quadratic + sine + cosine

               CALL FOILGRD (NU, ZERO, ONE, P (1), P (2), P (3), P (4),
     >                       TNORM)
            END IF

C           Determine the distribution of T that produces a distribution
C           of the arc-length corresponding to TNORM.  Upper surface ...

            CALL BSARCMAP (C, TLE, C (6 + NCPTS), NU, TNORM, ARC,
     >                     T (NL), LUNCRT, IER)

            IF (IER .NE. 0) THEN
               WRITE (LUNCRT, 1010) 'BSARCMAP (upper)', IER
               GO TO 999
            END IF
            

C           Reuse the normalized distribution for the lower surface:

            DO I = 1, NL
               T (I) = ONE - TNORM (I)
            END DO

            DO I = 1, NL
               TNORM (I) = T (NL + 1 - I)
            END DO

            CALL BSARCMAP (C, C (6 + NDEG), TLE, NL, TNORM, ARC, T (1),
     >                     LUNCRT, IER)

            IF (IER .NE. 0) THEN
               WRITE (LUNCRT, 1010) 'BSARCMAP (lower)', IER
               GO TO 999
            END IF
            

C           Evaluate the curve at the Ts corresp. to the specified arc values:

            NEVAL = NL + NU - 1
            DO I = 1, NEVAL

               CALL DTSPVL (T (I), C, WORK, MAXWORK, XY, IER)

               IF (IER .NE. 0) THEN
                  WRITE (LUNCRT, 1010) 'DTSPVL', IER
                  GO TO 999
               END IF

               X (I) = XY (1)
               Y (I) = XY (2)
            END DO

            X (NEVAL) = C (NCPTS + IX)      ! Avoids T > last knot by round-off
            Y (NEVAL) = C (NCPTS + IY)

            WRITE (LUNPLT, 1001)
     >         ' $OPTIONS LINE=''CONNECTED'',',
     >         ' LEGEND = ''Arc-length-based grid'',',
     >         ' $END'
            WRITE (LUNPLT, 1003) (X (I), Y (I), I = 1, NEVAL)
            WRITE (LUNPLT, 1004)
         END IF

      ELSE

C        Arbitrary curve case:
C        ---------------------

         IF (INSERT) THEN

            NEWSPLINE = .TRUE.
            TINSERT = HALF
            CALL READD (LUNCRT,
     >         'Knot value to insert?  [0.5] ',
     >         LUNKBD, TINSERT, CR, EOF)
            IF (EOF) GO TO 999

            MULTI = 1
            CALL READI (LUNCRT, 
     >         'Multiplicity?  [1] ', LUNKBD, MULTI, CR, EOF)
            IF (EOF) GO TO 999 
            MULTI = MAX (1, MIN (MULTI, NDEG))

            CALL BSINSERT (MAXC, C, TINSERT, MULTI, LASTCPT, C, IER)

            IF (IER .NE. 0) THEN
               WRITE (LUNCRT, 1010) 'BSINSERT', IER
               GO TO 999
            END IF

            IF (LASTCPT .GT. 0) THEN
               NCPTS = C (4)
               WRITE (LUNCRT, 1006)
     >            'Last control point number inserted: ', LASTCPT
            ELSE
               NEWSPLINE = .FALSE.
               WRITE (LUNCRT, 1006)
     >            'Knot insertion suppressed as invalid.'
            END IF

         END IF

         IF (PLOT) THEN

  500       NEVAL = 500
            CALL READI (LUNCRT,
     >         'How many points to be plotted?  0 is OK.  [500] ',
     >         LUNKBD, NEVAL, CR, EOF)
            IF (EOF) GO TO 999

            IF (NEVAL .GT. MAXPLT) THEN
               WRITE (LUNCRT, '(A, I3)') ' Present limit: ', MAXPLT
               GO TO 500
            END IF

            PLOT = NEVAL .GT. 0

         END IF

         IF (PLOT) THEN

            DELTAT = (C (6 + NCPTS) - C (6 + NDEG)) / (NEVAL - 1)

            DO I = 1, NEVAL - 1

               TEVAL = C (6 + NDEG) + DELTAT * (I - 1)

               CALL DTSPVL (TEVAL, C, WORK, MAXWORK, XY, IER)

               IF (IER .NE. 0) THEN
                  WRITE (LUNCRT, 1010) 'DTSPVL', IER
                  GO TO 999
               END IF

               X (I) = XY (1)
               Y (I) = XY (2)
            END DO

            X (NEVAL) = C (NCPTS + IX)      ! Avoids T > last knot by round-off
            Y (NEVAL) = C (NCPTS + IY)

            WRITE (LUNPLT, 1001) '""'
            WRITE (LUNPLT, 1012) 'Degree: ', NDEG,
     >         '    Control points: ', NCPTS

            WRITE (LUNPLT, 1001)
            WRITE (LUNPLT, 1001)
     >         ' $OPTIONS WIDTH=6., $END',
     >         ' $OPTIONS LINE=''SOLID'',',
     >         ' LEGEND=''B-spline curve'', $END'

            WRITE (LUNPLT, 1003) (X (I), Y (I), I = 1, NEVAL)
            WRITE (LUNPLT, 1004)

         END IF

         IF (GRID) THEN
            WRITE (LUNCRT, 1001)
     >         ' Arc-length-based discretization not implemented yet.'
         END IF

      END IF


C     Save [a new version of] the B-spline curve?
C     -------------------------------------------

      IF (NEWSPLINE) THEN

C        Open as "new" if possible - normally the VMS preference:

         I1 = INDEX (FILENAME, ']') + 1
         I2 = INDEX (FILENAME, ';') - 1
         IF (I2 .LT. 0) I2 = INDEX (FILENAME, ' ') - 1
         IF (I2 .LT. 0) I2 = LEN (FILENAME)

         OPEN (UNIT=LUNDISK, FILE=FILENAME (I1 : I2), STATUS='NEW',
     >         IOSTAT=IOSTAT)

         IF (IOSTAT .NE. 0) THEN      ! Unix/single-version?
            FILENAME = SFILE
            I1 = 1
            I2 = LEN (SFILE)  ! For end-report of output files

            OPEN (UNIT=LUNDISK, FILE=SFILE, STATUS='UNKNOWN',
     >            IOSTAT=IOSTAT)
            IF (IOSTAT .NE. 0) THEN   ! Unrecoverable
               WRITE (LUNCRT, 1006)
     >            'Error opening output B-spline file.  Status: ',
     >            IOSTAT
               GO TO 999
            END IF
         END IF

         CALL DTRDWR ('W', LUNDISK, TITLE, MAXC, C, NC, IER)

         IF (IER .NE. 0) THEN
            WRITE (LUNCRT, 1010) 'DTRDWR (writing)', IER
            GO TO 999
         END IF

         CLOSE (UNIT=LUNDISK)

      END IF

C     Show the control points ...

      WRITE (LUNPLT, 1001) ' $OPTIONS'
      WRITE (LUNPLT, 1002) ' LEGEND=', 'Control polygon', ' LINE=','DOT'
      WRITE (LUNPLT, 1001) ' SYMBOL=5,', ' $END'
      WRITE (LUNPLT, 1003) (C (I + IX), C (I + IY), I=1, NCPTS)
      WRITE (LUNPLT, 1004)

C     ... and the knots:

      X (1) = C (1 + IX)
      Y (1) = C (1 + IY)

      NK = NCPTS - NDEG + 1  ! Avoid end multiplicities

      DO I = 2, NK - 1

         CALL DTSPVL (C (I + 5 + NDEG), C, WORK, MAXWORK, XY, IER)

         IF (IER .NE. 0) THEN
            WRITE (LUNCRT, 1010) 'DTSPVL', IER
            GO TO 999
         END IF

         X (I) = XY (1)
         Y (I) = XY (2)

      END DO

      X (NK) = C (NCPTS + IX)
      Y (NK) = C (NCPTS + IY)

      WRITE (LUNPLT, 1001) ' $OPTIONS'
      WRITE (LUNPLT, 1002) ' LEGEND=', 'Spline evaluated at the knots',
     >                     ' LINE=', 'SYMBOLS'
      WRITE (LUNPLT, 1001) ' SYMBOL=16,', ' $END'
      WRITE (LUNPLT, 1003) (X (I), Y (I), I = 1, NK)


C     Show the quality of the data fit?

      IF (FIT) THEN

         IF (.NOT. INTERPOLATE) THEN

            X (1) = C (1 + IX)
            Y (1) = C (1 + IY)

            DO I = 2, NXY - 1

               CALL DTSPVL (BSU (I), C, WORK, MAXWORK, XY, IER)

               IF (IER .NE. 0) THEN
                  WRITE (LUNCRT, 1010) 'DTSPVL', IER
                  GO TO 999
               END IF

               X (I) = XY (1)
               Y (I) = XY (2)

            END DO

            X (NXY) = C (NCPTS + IX)
            Y (NXY) = C (NCPTS + IY)

            WRITE (LUNPLT, 1001) ' $OPTIONS'
            WRITE (LUNPLT, 1002) ' LEGEND=',
     >         'Spline evaluated at the data point parameter values',
     >         ' LINE=', 'SYMBOLS'
            WRITE (LUNPLT, 1001) ' SYMBOL=4,', ' $END'
            WRITE (LUNPLT, 1003) (X (I), Y (I), I = 1, NXY)

         END IF

      END IF


      IF (NEWSPLINE) THEN
         WRITE (LUNCRT, 1009)
     >      'Output B-spline curve file:  ', FILENAME (I1 : I2)
      ELSE
         WRITE (LUNCRT, 1001)
      END IF

      WRITE (LUNCRT, 1008) ' Plottable file:              ', PFILE

      IF (FIT) THEN
         IF (.NOT. SAMEKNOTS)
     >      WRITE (LUNCRT, 1008) ' File of reusable U values:   ', UFILE
      END IF

  999 WRITE (LUNCRT, 1001)

C     Avoid STOP problems by omitting it altogether.

C     Formats:

 1001 FORMAT (A)
 1002 FORMAT (A, '''', A, ''',')
 1003 FORMAT (2F24.16)
 1004 FORMAT ('END CURVE')
 1005 FORMAT (A, 2X, F24.16)
 1006 FORMAT (/, 1X, A, I3, A, I4)
 1007 FORMAT (1X, A)
 1008 FORMAT (A, A)
 1009 FORMAT (/, 1X, A, A)
 1010 FORMAT (/, ' *** Bad return from ', A, ':  IER = ', I3)
 1011 FORMAT (4X, A)
 1012 FORMAT (A, I3, A, I4)

      END
