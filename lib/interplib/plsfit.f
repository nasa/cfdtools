C+----------------------------------------------------------------------
C
      SUBROUTINE PLSFIT (NDATA, X, Y, TBEGIN, TEND, NEVAL, XEVAL, YEVAL,
     &   NEW, CLOSED, METHOD, DISTRIB, IER)
C
C     One-liner:  Storage-efficient parametric piecewise cubic fit (2-space)
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        PLSFIT (Parametric Local Spline FIT) is intended for general
C     curve drawing applications where it is not necessary to provide
C     spline coefficients as output. It is also suited to interpolating
C     the same data repeatedly at maximum efficiency (as in interactive
C     graphics, perhaps). Application to grid generation is another
C     possibility (as in redistributing X and Y vs. T where preferred
C     values of T are determined externally).
C
C        "Local" piecewise methods typically produce continuity in the
C     function and its first derivative across the data points, but the
C     second derivatives are not guaranteed to be continuous.  With
C     reasonably smooth data, this limitation should not be a problem.
C
C        Since it avoids storing spline coefficients, PLSFIT requires
C     very little working storage, yet its speed is on a par with other
C     methods in a single pass. The routine was designed for evaluation
C     of either many points at once (preferred) or one point at a time.
C     Computational efficiency is achieved by retaining some internal
C     quantities between calls for re-use where appropriate.
C
C        The primary method employed, references (1)-(3), uses monotone
C     piecewise cubic interpolation for X and Y vs. T, so the resulting
C     curve follows the data closely. In fact, the fit is sometimes too
C     tight, especially where X or Y appears constant over a region
C     which is intended to be curved. (A circle with unfortunate choice
C     of data points can develop flat spots.)  A partial solution is
C     offered by the "Bessel" option, which produces a more rounded
C     interpretation of the data. Both techniques are based on matching
C     the first (but not second) derivatives of the cubics on adjacent
C     intervals. The structure of PLSFIT permits easy addition of other
C     local, 3-point methods for these derivatives.
C
C        Output arrays XEVAL and YEVAL are filled by sampling the fit at
C     intervals along the arc's length, which is estimated by summing
C     the chords between adjacent points. For curve drawing (NEVAL large),
C     this version of PLSFIT offers only equal increments of arclength,
C     but other options may be added in the future (a method based on
C     curvature, for example). Other applications may be require calling
C     the subroutine repeatedly with NEVAL = 1 to evaluate the fit at
C     externally-generated values of T.  The evaluation can proceed
C     efficiently "backwards" (decreasing T) if desired.
C
C
C     Curve-drawing example:
C     ----------------------
C
C        PLSFIT usage, with extra-careful prescaling to avoid possible
C     overflow or underflow in chord length calculations, might look like
C     this for a curve-drawing application:
C
C        CALL RANGER (..., X, XSCALE, XSHIFT)   ! Scan input arrays
C        CALL RANGER (..., Y, YSCALE, YSHIFT)
C
C        CALL RESCALE (..., X, XSCALE, XSHIFT)  ! Transform data to [0, 1]
C        CALL RESCALE (..., Y, YSCALE, YSHIFT)
C
C        CALL PROTECT (..., DISTINCT)           ! Check for zero chords
C        IF (.NOT.DISTINCT) GO TO 900
C
C        NEW     = .TRUE.                       ! Initialize PLSFIT parameters
C        CLOSED  = .FALSE.
C        DISTRIB = 'U'                          ! Uniform (upper case)
C        METHOD  = 'M'                          ! Monotone ( "     " )
C        NFIT    = 500                          ! Some reasonable number
C        TBEGIN  = -1.0                         ! Interpolate entire curve
C        TEND    = -1.0
C
C        CALL PLSFIT (NXY, X, Y, TBEGIN, TEND, NFIT, XFIT, YFIT, NEW,
C       &   CLOSED, METHOD, DISTRIB, IER)
C        IF (IER .NE. 0) GO TO 910
C
C        Restore result using the inverse transformation:
C
C        CALL RESCALE (..., XFIT, ONE / XSCALE, -XSHIFT / XSCALE)
C        CALL RESCALE (..., YFIT, ONE / YSCALE, -YSHIFT / YSCALE)
C        (May need to restore the data also.)
C
C     Utilities RANGER, RESCALE, and PROTECT are available from the author.
C     (More recent utilities GETSCALE and USESCALE may be preferable now.)
C     They should only be needed if the magnitude and spacing of the data
C     points are unsuitable for ordinary single-precision calculations (all
C     too frequently the case, alas!). CHORD, which is used repeatedly by
C     PLSFIT, may also be invoked by the calling routine to provide initial
C     and final arclengths associated with particular array indices:
C
C        TBEGIN = CHORD (X, Y, 1, IBEGIN)
C        TEND   = TBEGIN + CHORD (X, Y, IBEGIN, IEND)
C
C     If the entire curve is to be interpolated at internally-generated
C     Ts, it is more efficient to pass -1.0s (indicating default choice
C     of arclength) and let PLSFIT compute (and re-use) the total arclength.
C
C     
C     Grid generation example:
C     ------------------------
C
C        A grid generation application may need to compute total arclength,
C     redistribute the intermediate values of T, and interpolate to these
C     by calling PLSFIT once per value of T.  Again, using PLSFIT to do
C     the initial calculation of TEND is preferable to direct use of CHORD,
C     because PLSFIT needs to do other initialization on its first call
C     if unnecessary recalculation of TEND is to be avoided:
C
C        Determine total arclength:
C
C        NEW     = .TRUE.                  ! Initialize PLSFIT parameters;
C        NEVAL   = 1                       ! force one pass through the
C                                          ! search to initialize it too.
C        TBEGIN  = -1.0                    ! Causes TEND to be returned;
C        TEND    = -1.0                    ! TBEGIN returns as zero.
C        ..............
C
C        CALL PLSFIT (NXY, X, Y, TBEGIN, TEND, NEVAL, XEVAL, YEVAL, NEW,
C       &   CLOSED, METHOD, DISTRIB, IER)
C        IF (IER .NE. 0) GO TO 910
C
C        Redistribute the total arclength:
C
C        CALL DISTRIB (..., NEVAL, ZERO, TEND, ..., T, ...)
C
C        Evaluate the redistributed X and Y values one pair at a time
C        to avoid PLSFIT's internally-generated Ts:
C
C        NEW = .FALSE.
C        DO 300, I = 1, NEVAL
C           CALL PLSFIT (NXY, X, Y, T (I), T (I), 1, XEVAL (I), ...)
C    300 CONTINUE
C        ..............
C
C
C     Arguments:
C     ----------
C
C     Name    Type/Dimension  I/O/S  Description
C     NDATA   I               I      Length of X, Y input data arrays.
C
C     X,      R (NDATA)       I      Input data coordinates. Successive
C     Y                              data points must be distinct in the
C                                    sense that DX ** 2 + DY ** 2 > 0, but
C                                    no check is performed at this level.
C                                    If CLOSED, then first and last points
C                                    must agree (not checked here).
C
C     TBEGIN  R               I/O    Specifies the arclength at which to
C                                    start interpolation. If NEVAL = 1,
C                                    TBEGIN is the only point at which
C                                    the fit will be evaluated. If both
C                                    TBEGIN and TEND are negative, a
C                                    default value of zero will be used
C                                    for TBEGIN.  In this case, TBEGIN
C                                    is returned as zero and TEND is
C                                    returned as the total arclength.
C
C     TEND    R               I/O    The arclength at which to end the
C                                    interpolation, unless both TBEGIN
C                                    and TEND are negative, in which case
C                                    the total chord length will be used
C                                    for TEND. If NEVAL = 2, TEND is the
C                                    second (and last) at which the fit
C                                    will be evaluated.  See TBEGIN and
C                                    NEVAL.
C
C     NEVAL   I               I      Number of points at which the fit is
C                                    to be evaluated.  NEVAL >= 1.
C                                    NEVAL > 1 means PLSFIT generates its
C                                    own values of T according to DISTRIB;
C                                    NEVAL = 1 with TBEGIN, TEND < 0. enables
C                                    PLSFIT to calculate total arclength
C                                    while initializing itself for further
C                                    interpolation at externally-supplied
C                                    values of T.  See the grid-generation
C                                    example above.
C
C     XEVAL,  R (NEVAL)         O    Interpolated coordinates.
C     YEVAL
C
C     NEW     L               I      If control flag NEW is .TRUE., the
C                                    search for a bracket starts from
C                                    scratch, otherwise locally-saved
C                                    search and fit information will be
C                                    assumed to be correct. If calling
C                                    PLSFIT from within a loop, set
C                                    NEW = .FALSE. after the first call.
C
C     CLOSED   L              I      Logical flag indicating that periodic
C                                    boundary conditions are to be used.
C                                    (The curve should wrap around smoothly
C                                    on itself.). The calling routine must
C                                    ensure that the ends agree.
C
C     METHOD   C*1            I      The type of fit to be used:
C                                    'M' means monotonic piecewise cubics;
C                                    'B' means non-monotonic "Bessel"-type
C                                    piecewise cubics (looser fit).
C                                    METHOD must be uppercase.
C
C     DISTRIB  C*1            I      Controls the distribution of the
C                                    interpolated points along the curve.
C                                    Irrelevant if NEVAL = 1, which is
C                                    appropriate for externally-generated
C                                    values of arclength for interpolation.
C                                    The only choice at present is 'U' for
C                                    uniform arclength increments. DISTRIB
C                                    must be uppercase.
C
C     IER      I                O    Error flag.  The possibilities are:
C                                       0:  No problems
C                                      +1:  NDATA < 2 (all cases), or
C                                           NDATA < 3 (periodic case)
C                                      +2:  Ends don't match (periodic)
C                                      +3:  Bad METHOD requested
C                                      +4:  Bad DISTRIB requested
C
C     Significant local variables:
C     ----------------------------
C
C     MEMORY         Indicates that search information and cubic coefficients
C                    are correct for the current point.
C
C     NEWFLAG        A local copy of NEW which is set .FALSE. during the
C                    first entry to the interval search routine, ARCSRCH,
C                    which SAVEs several local variables from one call to
C                    the next.
C
C     LEFT, RIGHT    Current endpoints of the bracketing interval.
C     TLEFT, TRIGHT  Corresponding cumulative arclengths.
C
C     IND            Index array which keeps track of the endpoints of
C                    the intervals preceding and following the bracketing
C                    interval, with wrap-around at the extremes.
C
C     H, DEL         Arclength and forward difference derivative arrays.
C
C     BX, CX, DX     Coefficients of X-cubic on the bracketing interval.
C     BY, CY, DY     Ditto, for Y-cubic.
C
C     Procedures:
C     -----------
C
C     ARCSRCH   Interpolation search along chord of a parametric curve.
C     BESSEL    First derivative using central 3-point formula.
C     BRODLIE   First derivative (central), adjusted for monotonicity.
C     BUTLAND   First derivative (non-central), adjusted for monotonicity.
C     CHORD     Summed chord-lengths for X-Y curve over range of indices.
C     THREEPT   First derivative using non-central 3-point formula.
C
C     Environment:  Digital VAX-11/780, VMS FORTRAN 77
C     ------------  Apple Macintosh, Absoft MacFORTRAN/020
C
C     Notes:
C     ------
C
C     (1)  IMPLICIT NONE, 8-character symbolic names, and "!" as comment
C          character are not (yet) standard.
C
C     (2)  Since many of the calculations must be repeated at both ends
C          of an interval, and for both X and Y data, the various finite
C          difference quantities used are stored as arrays. The following
C          "map" of a typical interior interval and its neighbors should
C          help in understanding the notation.  The local array indices
C          are all numbered relative to the left-hand end of the interval
C          which brackets the point to be evaluated.
C
C                                  LEFT       RIGHT
C
C          Point         -1          0         +1          +2
C
C          Data           x -------- x -------- x --------- x
C
C          Interval      -1          0         +1
C
C     (3)  Within PLSFIT, we have allowed for extrapolation to simplify
C          things, but it's not intended for general use (especially for
C          the periodic case where we don't attempt to wrap around past
C          the last interval).
C
C     (4)  Error-checking philosophy: several cheap tests are performed
C          upon entry (see IER description, above), but the user must
C          guarantee that successive points are distinct to avoid division
C          by arclength intervals of length zero. Such a test is best
C          performed at a higher level to avoid the considerable overhead
C          of scanning both data arrays on each entry.
C
C     (5)  Note that if some function of the data is to be plotted,
C          e.g., logarithm, then a smoother curve will be obtained by
C          passing transformed data, fitting at equal increments, then 
C          reverse-transforming the result.
C
C     Bibliography:
C     -------------
C
C     (1)  Brodlie, K. W.  A Review of Methods for Curve and Function
C             Drawing, in Mathematical Methods in Computer Graphics and
C             Design, ed. K. W. Brodlie.  London: Academic Press, 1980.
C             Pp. 1-37.  (See esp. the discussion, pp. 33-37)
C
C     (2)  Fritsch, F. N., and J. Butland.  A Method for Constructing
C             Local Monotone Piecewise Cubic Interpolants.  SIAM J. Sci.
C             Stat. Comput., Vol. 5, No. 2 (June 1984).  Pp. 300- 304.
C
C     (3)  Fritsch, F. N., and R. E. Carlson.  Monotone Piecewise Cubic
C             Interpolation.  SIAM J. Num. Anal., Vol. 17, No. 2
C             (April 1980).  Pp. 238-246.
C
C     (4)  Sedgewick, R.  Algorithms.  Reading: Addison-Wesley, 1983.
C             Chap. 14.  (Interpolation search)
C
C     Author:  Robert Kennelly, ex Sterling Software, now RAC Branch
C     -------
C              Mail Stop 227-2
C              NASA-Ames Research Center
C              Moffett Field, CA  94035
C
C              Phone (415) 604-5860
C
C     History:
C     --------
C
C     27 Feb. 1987    RAK    Initial design and coding.
C     16 Apr. 1987    RAK    Repaired handling of the degenerate linear
C                            case - DELX and DELY were not being set.
C     22 Apr. 1987    RAK    Fixed typo in linear case (DY not set).
C     25 Mar. 1988    RAK    Modularized ARCSRCH, with repairs (now works
C                            bidirectionally). Both TBEGIN and TEND have
C                            to be negative to signal use of defaults
C                            (zero and total chord length). Pass total
C                            arclength to ARCSRCH on first entry to
C                            avoid redundant calculation.
C     12 Apr. 1988    RAK    Set MEMORY = .TRUE. if possible when PLSFIT
C                            is called repeatedly with the same data.
C     23 Aug. 1989    DAS    Application to grid generation revealed that
C                            NEVAL = 0 (for arclength calculation) didn't
C                            initialize the search properly - had to go
C                            to NEVAL = 1.  Descriptions above revised
C                            to cover both graphics and grid generation.
C     20 June 1991    DAS    THREEPT was renamed BUTLAND (monotonic);
C                            THREEPT (pure 3-pt. formula) is now used
C                            for the loose fit at the boundaries.
C                            
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Arguments.

      LOGICAL
     &   CLOSED, NEW
      INTEGER
     &   IER, NDATA, NEVAL
      REAL
     &   TBEGIN, TEND, X (NDATA), XEVAL (NEVAL), Y (NDATA),
     &   YEVAL (NEVAL)
      CHARACTER
     &   DISTRIB * 1, METHOD * 1

C     Local constants.

      REAL
     &   ZERO, ONE, TWO, THREE
      PARAMETER
     &  (ZERO  = 0.0E+0,
     &   ONE   = 1.0E+0,
     &   TWO   = 2.0E+0,
     &   THREE = 3.0E+0)

C     Local variables.

      LOGICAL
     &   MEMORY, MONO, NEWFLAG
      INTEGER
     &   IEVAL, IND (-1:2), J, LEFT, RIGHT
      REAL
     &   BX (0:1), BY (0:1), CX, CY, DT, DELX (-1:1), DELY (-1:1),
     &   DX, DY, H (-1:1), RH, TEVAL, TINC, TLEFT, TRIGHT, TTOTAL

C     Procedures.

      REAL
     &   BESSEL, BRODLIE, BUTLAND, CHORD, THREEPT
      EXTERNAL
     &   BESSEL, BRODLIE, BUTLAND, CHORD, THREEPT

C     Storage.

      SAVE
     &   BX, BY, CX, CY, DX, DY, LEFT, RIGHT, TLEFT, TRIGHT, TTOTAL

C     Error checking and initialization.
C     ----------------------------------

      IER = 0
      IF (NDATA .LT. 2) THEN
         IER = +1
         GO TO 990
      END IF

      IF (CLOSED) THEN

C        At least three points are required, and the first and last
C        points must be identical.

         IF (NDATA .LT. 3) IER = +1
         IF (X (1) .NE. X (NDATA) .OR. Y (1) .NE. Y (NDATA)) IER = +2
      END IF

C     We'll need the total arclength to initialize the search efficiently.

      IF (NEW) THEN
         TTOTAL = CHORD (X, Y, 1, NDATA)
         TRIGHT = TTOTAL
      END IF

C     The whole curve will be interpolated by default if the calling
C     routine passes negative quantities for both TBEGIN and TEND. 

      IF (TBEGIN .LT. ZERO .AND. TEND .LT. ZERO) THEN
         TBEGIN = ZERO
         TEND   = TTOTAL
      END IF

C     Check the requested interpolation method and point distribution.

      MONO = METHOD .EQ. 'M'
      IF (.NOT. MONO .AND. METHOD .NE. 'B') IER = +3

      IF (NEVAL .GT. 1) THEN           ! NEVAL = 1 for T supplied externally.
         IF (DISTRIB .EQ. 'U') THEN    ! Uniform values of T generated here.
            TINC = (TEND - TBEGIN) / REAL (NEVAL - 1)
         ELSE
            IER = +4
         END IF
      END IF

C     Bail out if any error was found. (Only the last found is reported.)
C     Original test on NEVAL = 0 here was invalid - have to go through ARCSRCH.

      IF (IER .NE. 0) GO TO 990


C     Initialize bracket quantities. Note that when NEW = .TRUE. the
C     MEMORY flag is never set, thus ARCSRCH (below) gets initialized
C     properly, as long as it is entered at least once (NEVAL >= 1).

      TEVAL   = TBEGIN
      NEWFLAG = NEW

      IF (NEWFLAG) THEN
         MEMORY = .FALSE.
      ELSE
C        We can save a lot of time when PLSFIT is being called from within
C        a loop by setting MEMORY if possible. The out-of-range checking
C        relies on the fact that RIGHT = LEFT + 1 after the last return.
C        Cater to the more likely case of TEVAL in the previous, interior
C        interval.

         MEMORY = (TEVAL .GE. TLEFT) .AND. (TEVAL .LT. TRIGHT)
         IF (.NOT. MEMORY) THEN
            MEMORY = (LEFT  .EQ. 1)     .AND. (TEVAL .LT. TRIGHT) .OR.
     &               (RIGHT .EQ. NDATA) .AND. (TEVAL .GE. TLEFT)
         END IF
      END IF

C     Loop over evaluation points.
C     ----------------------------

      IEVAL = 0
   10 CONTINUE
         IEVAL = IEVAL + 1

         IF (MEMORY) GO TO 70        ! Skip the bulk of the computation.


C        Interpolation search for bracketing interval.
C        ---------------------------------------------

         CALL ARCSRCH (NDATA, X, Y, TEVAL, LEFT, TLEFT, RIGHT,
     &      TRIGHT, NEWFLAG)

C         -------------------------------------------------------------
C        |                                                             |
C        |   1 <= LEFT < RIGHT = LEFT + 1 <= NDATA                     |
C        |                                                             |
C         -------------------------------------------------------------

C        Compute derivatives by finite-differences.
C        ------------------------------------------

         IF (NDATA .GT. 2) THEN

C           Three cases are handled together: (1) both endpoints are
C           interior to the data range, or the interval is at the
C           beginning/end of the range with either (2) periodic,
C           or (3) free end conditions.  For special cases (2) and (3),
C           the initial calculations will be overridden below.

            IND (-1) = LEFT - 1
            IND ( 0) = LEFT
            IND (+1) = RIGHT
            IND (+2) = RIGHT + 1

C           Patch index array for periodic case (2). (Later, the free end
C           case will overridden again, following derivative calculation.)

            IF (LEFT .EQ. 1) THEN

C              Left side wrap-around boundary condition.

               IND (-1) = NDATA - 1
            ELSE IF (RIGHT .EQ. NDATA) THEN

C              Right side.

               IND (2) = 2
            END IF

C           Interval and derivative approximations.
C           ---------------------------------------

C           Eliminate possible division by zero due to cancellation while
C           subtracting by computing the chord from LEFT to RIGHT explicitly.

            H (-1) = CHORD (X, Y, IND (-1), IND (-1) + 1)
            H ( 0) = CHORD (X, Y, LEFT, RIGHT)
            H (+1) = CHORD (X, Y, IND (2) - 1, IND (2))

            DO 40, J = -1, +1
               RH       = ONE / H (J)
               DELX (J) = (X (IND (J + 1)) - X (IND (J))) * RH
               DELY (J) = (Y (IND (J + 1)) - Y (IND (J))) * RH
   40       CONTINUE

C           Select interpolation scheme.
C           ----------------------------

C           Compute adjusted X and Y derivatives at both left- and
C           right-hand endpoints of the interval.

            IF (MONO) THEN

C              Monotone - use Brodlie modification of Butland's
C              formula to adjust the derivatives at the knots.

               DO 50, J = 0, +1
                  BX (J) = BRODLIE (J, H, DELX)
                  BY (J) = BRODLIE (J, H, DELY)
  50           CONTINUE

            ELSE                             ! IF (METHOD .EQ. 'B') THEN

C              Bessel - use central difference formula at the knots.

               DO 60, J = 0, +1
                  BX (J) = BESSEL (J, H, DELX)
                  BY (J) = BESSEL (J, H, DELY)
  60           CONTINUE
            END IF

C           Patch initial/final derivatives if not periodic, case (3).

            IF (.NOT. CLOSED) THEN
               IF (LEFT .EQ. 1) THEN
                  IF (.NOT. MONO) THEN
                     BX (0) = THREEPT (0, H, DELX)
                     BY (0) = THREEPT (0, H, DELY)
                  ELSE
                     BX (0) = BUTLAND (0, H, DELX)
                     BY (0) = BUTLAND (0, H, DELY)
                  END IF
               ELSE IF (RIGHT .EQ. NDATA) THEN
                  IF (.NOT. MONO) THEN
                     BX (1) = THREEPT (1, H, DELX)
                     BY (1) = THREEPT (1, H, DELY)
                  ELSE
                     BX (1) = BUTLAND (1, H, DELX)
                     BY (1) = BUTLAND (1, H, DELY)
                  END IF
               END IF
            END IF

C           Compute the remaining cubic coefficients for X and Y relative
C           to the left-hand endpoint.

            RH = ONE / H (0)
            CX = (THREE * DELX (0) - TWO * BX (0) - BX (1)) * RH
            CY = (THREE * DELY (0) - TWO * BY (0) - BY (1)) * RH
            DX = (-TWO * DELX (0) + BX (0) + BX (1)) * RH ** 2
            DY = (-TWO * DELY (0) + BY (0) + BY (1)) * RH ** 2

         ELSE                                   ! IF (NDATA .EQ. 2) THEN

C           Degenerate case (linear).
C           -------------------------

            H (0)  = TRIGHT - TLEFT
            RH     = ONE / H (0)
            BX (0) = (X (RIGHT) - X (LEFT)) * RH
            BY (0) = (Y (RIGHT) - Y (LEFT)) * RH
            CX     = ZERO
            CY     = ZERO
            DX     = ZERO
            DY     = ZERO
         END IF

C        Evaluate the cubics for XEVAL and YEVAL.
C        ----------------------------------------

   70    CONTINUE

         DT            = TEVAL - TLEFT
         XEVAL (IEVAL) = X (LEFT) + DT * (BX (0) + DT * (CX + DT * DX))
         YEVAL (IEVAL) = Y (LEFT) + DT * (BY (0) + DT * (CY + DT * DY))

C        Choose next evaluation point and loop back.
C        -------------------------------------------

         IF (IEVAL .LT. NEVAL) THEN           ! Skips this if NEVAL = 1.

C           Uniform spacing - this fast, simple update should suffice
C           for graphics despite roundoff accumulation, especially since
C           extrapolation is allowed. (NOTE: Add test for DISTRIB here if
C           other distribution options are added.)

            TEVAL  = TEVAL + TINC
            MEMORY = (TEVAL .LT. TRIGHT) .AND. (TEVAL .GE. TLEFT)
            GO TO 10

         END IF

C     Termination.
C     ------------

  990 CONTINUE
      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE ARCSRCH (N, X, Y, T, LEFT, TLEFT, RIGHT, TRIGHT, NEW)
C
C     One-liner: Interpolation search along chord of a parametric curve.
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        A modular implementation of the arclength interval search
C     originally written for PLSFIT (Parametric Local Spline FIT). An
C     arclength interval is sought which contains some specified value,
C     T, or which is at least the nearest interval if T is off-scale.
C     The condition for a bracket is: TLEFT <= T <= T (LEFT+1). A logical
C     flag, NEW, must be set .TRUE. on the first call for a new set of
C     data, and the total arclength must be supplied as the initial value
C     of input variable TRIGHT. Subsequent calls with the same data make
C     use of local SAVEd variables to avoid unnecessary recalculation.
C
C        There is minimal error checking in this low-level routine. The
C     calling program is assumed to have verified that N >= 2. Efficiency
C     will benefit from passing the best estimate available, usually just
C     the result of the last call.
C
C        This is not a fully independent utility - a number of quantities
C     are shared between PLSFIT and ARCSRCH. Speed is thereby enhanced at
C     the expense of generality and ease-of-use. With some care, ARCSRCH
C     could be transplanted into another setting (see PLSFIT for usage).
C
C        The interpolation search was adapted from ideas in Sedgewick's
C     book, referenced below.
C
C     Arguments:
C     ----------
C
C     Name     Dimension  Type  I/O/S  Description
C     N                    I    I      Number of points in arrays X, Y;
C                                      must be >= 2; no check performed.
C
C     X, Y     N           R    I      Array of distinct points defining
C                                      the set of intervals to be examined;
C                                      not checked.
C
C     T                    R    I      The distance along the chord for
C                                      which a bracketing interval is sought.
C                                      Normally in the range zero to total
C                                      chord length, but may lie outside.
C
C     NOTE: If NEW = .TRUE., only the value of TRIGHT needs to be supplied.
C
C     LEFT                 I    I/O    Input: estimate of index of left
C                                      endpoint of the interval containing
C                                      the specified point. Must be in the
C                                      range [1, N-1]; not error checked.
C                                      LEFT < RIGHT is assumed.
C
C                                      Output: index of the largest array
C                                      value <=, or sometimes <, specified
C                                      point - see Notes. Special case for
C                                      data out of range: returns left
C                                      endpoint of closest interval.
C
C     TLEFT                R    I/O    Arclength up to index LEFT.
C
C     RIGHT                I    I/O    Input: estimate of index of right
C                                      endpoint of the interval containing
C                                      the specified point. Must be in the
C                                      range [2, N]; not error checked.
C                                      RIGHT > LEFT is assumed.
C
C                                      Output: index of the largest array
C                                      value >, or sometimes >=, specified
C                                      point - see Notes. Special case for
C                                      data out of range: return right
C                                      endpoint of closest interval.
C
C     TRIGHT               R    I/O    Arclength up to index RIGHT. NOTE:
C                                      when NEW = .TRUE., TRIGHT must be
C                                      the total arclength of the curve.
C                                      This trick permits PLSFIT & ARCSRCH
C                                      to avoid unnecessary recalculation.
C
C     NEW                  L    I/O    Must be set .TRUE. when ARCSRCH is
C                                      first called for any dataset.
C
C                                      Set to .FALSE. on output.
C
C     Environment:  Digital VAX-11/780, VMS FORTRAN
C     ------------  Apple Macintosh, Absoft MacFORTRAN/020 v2.3
C
C     Notes:
C     ------
C
C     (1)  IMPLICIT NONE and 8-character symbols are not (yet) standard.
C
C     (2)  The arclength calculations are approximations based on straight
C          line segments between the data points.
C
C     (3)  The algorithm is designed to return bracket endpoints such that
C          TLEFT <= T < TRIGHT = T (LEFT+1), with strict inequality on the
C          right, but roundoff errors in the chord calculations sometimes
C          interfere. We protect the routine from failures such as collapse
C          of the interval to zero, i.e., LEFT = RIGHT, by checking ahead
C          when the endpoints are adjusted during the main loop, and simply
C          accept the occasional odd result. Note also that the arclengths
C          of the interval's endpoints should be regarded as being slightly
C          fuzzy, perhaps a few parts in 10**7 for typical single precision.
C          Neither condition is likely to be a problem in typical graphics
C          applications, and the calling routine can easily check the result
C          if it is critical.
C
C     Bibliography:
C     -------------
C
C     (1) Sedgewick, R.  Algorithms.  Reading: Addison-Wesley, 1983.
C            (Chap. 14)
C
C     Author:  Robert Kennelly, Sterling Software/NASA Ames, Palo Alto, CA.
C     -------
C
C     History:
C     --------
C
C     25 Mar. 1988    RAK    Adapted from INTERVAL and PLSFIT.
C     24 Aug. 1988    RAK    Added protection against roundoff-induced
C                            failure in main loop (LEFT = RIGHT). Revised
C                            loop for single termination (pure WHILE).
C                            Take advantage if TEMP2 is already known.
C                            Retain URTRIGHT instead of re-calculating.
C     10 July 1989    DAS    URTRIGHT was integer by mistake.
C
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     &   ZERO
      PARAMETER
     &  (ZERO = 0.0E+0)

C     Arguments.

      LOGICAL
     &   NEW
      INTEGER
     &   LEFT, N, RIGHT
      REAL
     &   TLEFT, TRIGHT, X (N), Y (N), T

C     Local variables.

      INTEGER
     &   LENGTH, NLESS1, TRIAL
      REAL
     &   LSENTRY, RSENTRY, TEMP1, TEMP2, URTRIGHT

C     Storage.

      SAVE
     &   NLESS1, LSENTRY, RSENTRY, URTRIGHT

C     Procedures.

      REAL
     &   CHORD
      EXTERNAL
     &   CHORD

C     Execution.
C     ----------

      IF (NEW) THEN

C        Initialization for new data set. Note special use of TRIGHT here,
C        to avoid having to repeat the relatively expensive summation over
C        the entire curve.

         NEW      = .FALSE.
         URTRIGHT = TRIGHT
         NLESS1   = N - 1

         LSENTRY  = CHORD (X, Y, 1, 2)
         RSENTRY  = TRIGHT - CHORD (X, Y, NLESS1, N)

C        The following will be appropriate only if the "simplification"
C        below doesn't apply, but we initialize everything here to get it
C        over with.

         LEFT   = 2
         TLEFT  = LSENTRY
         RIGHT  = NLESS1
         TRIGHT = RSENTRY
      END IF

C     Simplify things by disposing of two important special cases so that
C     TLEFT and TRIGHT can really bracket T. As a byproduct, this also
C     takes care of the N = 2, 3 cases (one or two intervals).

      IF (T .LT. LSENTRY) THEN       ! First interval applies
         LEFT   = 1
         TLEFT  = ZERO
         RIGHT  = 2
         TRIGHT = LSENTRY
         GO TO 990
      ELSE IF (T .GE. RSENTRY) THEN  ! Last interval applies
         LEFT   = NLESS1
         TLEFT  = RSENTRY
         RIGHT  = N
         TRIGHT = URTRIGHT
         GO TO 990
      END IF

C      ----------------------------------------------------------------
C     |                                                                |
C     |   N > 3                                                        |
C     |                                                                |
C     |   1 <= LEFT < RIGHT <= N                                       |
C     |                                                                |
C     |   LSENTRY <= T < RSENTRY                                       |
C     |                                                                |
C      ----------------------------------------------------------------

C     Refine bracket estimate, checking in particular whether the current
C     values are already correct.

      IF (T .GE. TLEFT) THEN
         IF (T .LT. TRIGHT) THEN

C           T is somewhere in the original interval - are we done?

            IF (RIGHT - LEFT .LE. 1) GO TO 990
         ELSE

C           We'll look farther to the right. Evidently LEFT was < N - 2.

            LEFT   = RIGHT
            TLEFT  = TRIGHT
            RIGHT  = NLESS1
            TRIGHT = RSENTRY
         END IF
      ELSE

C        Look to the left of the guess. Evidently LEFT was > 2.

         RIGHT  = LEFT
         TRIGHT = TLEFT
         LEFT   = 2
         TLEFT  = LSENTRY
      END IF

C      ----------------------------------------------------------------
C     |                                                                |
C     |   2 <= LEFT < RIGHT <= N - 1                                   |
C     |                                                                |
C     |   LSENTRY <= TLEFT <= T < TRIGHT <= RSENTRY                    |
C     |                                                                |
C      ----------------------------------------------------------------

C     The interval length must decrease each search iteration. Terminate
C     when the interval length drops to 1.

   10 CONTINUE
         LENGTH = RIGHT - LEFT
         IF (LENGTH .GT. 1) THEN

C           The trial value is a "linear" estimate of the left-hand
C           endpoint of the interval bracketing the target T, with
C           protection against round-off (which can affect convergence).

            TRIAL = MIN (RIGHT - 1, LEFT + MAX (0, INT (REAL (LENGTH) *
     &         (T - TLEFT) / (TRIGHT - TLEFT))))

C            ----------------------------------------------------------
C           |                                                          |
C           |   2 <= LEFT <= TRIAL < RIGHT <= N - 1                    |
C           |                                                          |
C            ----------------------------------------------------------

C           Compute the cumulative arclength up to TRIAL as cheaply as
C           possible using previous results. (The search runs equally well
C           for an increasing or decreasing sequence of T's.)

            IF ((TRIAL - LEFT) .LE. (RIGHT - TRIAL)) THEN
               TEMP1 = TLEFT + CHORD (X, Y, LEFT, TRIAL)
            ELSE
               TEMP1 = TRIGHT - CHORD (X, Y, TRIAL, RIGHT)
            END IF

C           Similar trick for arclength up to TRIAL + 1 since we may well
C           already know it, e.g., just before termination.

            IF (RIGHT .EQ. TRIAL + 1) THEN
               TEMP2 = TRIGHT
            ELSE
               TEMP2 = TEMP1 + CHORD (X, Y, TRIAL, TRIAL + 1)
            END IF

C           Adjust the endpoints to reduce LENGTH as much as possible, but
C           not less than 1.

            IF (T .GE. TEMP2) THEN

C              Increase LEFT carefully to avoid overshoot (which can
C              occur due to roundoff error in chord calculations).

               IF (TRIAL + 1 .LT. RIGHT) THEN
                  LEFT   = TRIAL + 1
                  TLEFT  = TEMP2
               ELSE
                  LEFT   = RIGHT - 1
                  TLEFT  = TRIGHT - CHORD (X, Y, LEFT, RIGHT)
               END IF
            ELSE IF (T .LT. TEMP1) THEN

C              Decrease RIGHT.

               IF (TRIAL .GT. LEFT) THEN
                  RIGHT  = TRIAL
                  TRIGHT = TEMP1
               ELSE
                  RIGHT  = LEFT + 1
                  TRIGHT = TLEFT + CHORD (X, Y, LEFT, RIGHT)
               END IF
            ELSE

C              Adjust both LEFT and RIGHT. We're done since T is in the
C              interval [T (TRIAL), T (TRIAL+1)), but defer termination
C              until LENGTH is tested at the top of the loop.

               LEFT   = TRIAL
               TLEFT  = TEMP1
               RIGHT  = TRIAL + 1
               TRIGHT = TEMP2
            END IF
            GO TO 10

         END IF

C     Termination.
C     ------------

  990 CONTINUE
      RETURN
      END
