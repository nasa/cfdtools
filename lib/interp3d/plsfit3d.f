C+----------------------------------------------------------------------
C
      SUBROUTINE PLSFIT3D (NDATA, X, Y, Z, TBEGIN, TEND, NEVAL, XEVAL,
     &   YEVAL, ZEVAL, NEW, CLOSED, METHOD, DISTRIB, IER)
C
C     One-liner:  Parametric linear or cubic spline fit to arc in 3-space
C     ----------  (space efficient)
C
C     Description and usage:
C     ----------------------
C
C        PLSFIT3D is the XYZ analog of PLSFIT.  This version retains all
C     of the original features, with the addition of the linear option.
C     Its initial application is to grid generation rather than graphics.
C     See PLSFIT for a full description of this low-storage approach to
C     piecewise polynomial interpolation.
C
C        See PLSCURVE if you need first partial derivatives.
C
C        Application of PLSFIT3D to grid generation might appear as
C     follows, with extra-careful prescaling to avoid unnecessary
C     numerical effects attributable to the data units:
C
C
C        Determine normalization factors:
C
C        CALL GETSCALE ('G', 3, NDATA, X, Y, Z, SCALE, SHIFT, IER)   
C
C        <Possibly check for distinct points here (not shown).>
C
C        Normalize the original curve coordinates:
C
C        CALL USESCALE ('N', 3, NDATA, X, Y, Z, SCALE, SHIFT, IER)   
C
C        Use PLSFIT3D's option to compute total arclength only:
C
C        NEVAL  = 1               ! These four values give the
C        TBEGIN = -ONE            ! efficient way to initialize
C        TEND   = -ONE            ! PLSFIT3D.
C        NEW    = .TRUE.
C        ...............
C
C        CALL PLSFIT3D (NDATA, X, Y, Z, TBEGIN, TEND, NEVAL, XEVAL, ...)
C
C        Distribute the total arclength in some specified way:
C
C        CALL DISTRIB (..., NEVAL, ZERO, TEND, ...., T, ...)
C
C        Evaluate the interpolated curve coordinates one point at a time to
C        avoid the PLSFIT-inherited graphical prejudice towards internally-
C        generated Ts:
C
C        NEW = .FALSE.
C        DO 300, I = 1, NEVAL
C           CALL PLSFIT3D (NDATA, X, Y, Z, T (I), T (I), 1, XEVAL (I), ...)
C    300 CONTINUE
C
C        Restore the fitted result to the original units:
C
C        CALL USESCALE ('D', 3, NEVAL, XEVAL, Y..., SCALE, SHIFT, IER)   
C
C        (May need to restore the original data also.)
C
C
C        Utilities GETSCALE and USESCALE are available from the authors.
C     CHORD3D, which is used repeatedly by PLSFIT3D, may also be invoked
C     by the calling routine to provide initial and final arclengths
C     associated with particular array indices:
C
C        TBEGIN = CHORD3D (X, Y, Z, 1, IBEGIN)
C        TEND   = TBEGIN + CHORD3D (X, Y, Z, IBEGIN, IEND)
C
C     If the entire curve is to be interpolated at internally generated
C     Ts, it is most efficient to pass -1.0s (indicating default choice of
C     arclength) and let PLSFIT3D compute (and re-use) the total arclength.
C
C     Arguments:
C     ----------
C
C     Name    Type/Dimension  I/O/S  Description
C     NDATA   I               I      Length of X, Y, Z input data arrays.
C
C     X,      R (NDATA)       I      Input curve coordinates. Successive
C     Y,                             data points must be distinct (non-zero
C     Z                              distance between them), but no check
C                                    is performed at this level.
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
C                                    is returned as zero, and TEND is
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
C                                    NEVAL > 1 means PLSFIT3D generates its
C                                    own values of T according to DISTRIB;
C                                    NEVAL = 1 with TBEGIN, TEND < 0. enables
C                                    PLSFIT3D to calculate total arclength
C                                    while initializing itself for further
C                                    interpolation at externally-supplied
C                                    values of T.  See the grid-generation
C                                    example above.
C
C     XEVAL,  R (NEVAL)         O    Interpolated coordinates.
C     YEVAL,
C     ZEVAL
C
C     NEW     L               I      If control flag NEW is .TRUE., the
C                                    search for a bracket starts from
C                                    scratch, otherwise locally-saved
C                                    search and fit information will be
C                                    assumed to be correct. If calling
C                                    PLSFIT3D from within a loop, set
C                                    NEW = .FALSE. after the first call.
C
C     CLOSED   L              I      Logical flag indicating that periodic
C                                    boundary conditions are to be used.
C                                    (The curve should wrap around smoothly
C                                    on itself.) The calling routine must
C                                    ensure that the ends agree.
C
C     METHOD   C*1            I      The type of parametric fit to be used:
C                                    'M' means monotonic piecewise cubics;
C                                    'B' means non-monotonic "Bessel"-type
C                                        piecewise cubics (looser fit);
C                                    'L' means piecewise linear (only way
C                                        to guarantee no overshoots).
C                                    METHOD must be uppercase.
C
C     DISTRIB  C*1            I      Controls internal distribution of the
C                                    interpolated points along the curve.
C                                    Irrelevant if NEVAL = 1, which is
C                                    appropriate for externally distributed
C                                    values of arclength for interpolation.
C                                    The only choice at present is 'U' for
C                                    UNIFORM (equal arclength increments).
C                                    DISTRIB must be uppercase.
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
C                    first entry to the interval search routine, ASRCH3D,
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
C     BY, CY, DY     Ditto, for the Y-cubic...
C     BZ, CZ, DZ     ... and the Z-cubic.
C
C     Procedures:
C     -----------
C
C     ASRCH3D   Interpolation search along chord of a parametric curve.
C     BESSEL    First derivative (central 3-point formula).
C     BRODLIE   First derivative (central), adjusted for monotonicity.
C     BUTLAND   First derivative (non-central), adjusted for monotonicity.
C     CHORD3D   Summed chord-lengths for X-Y-Z curve over range of indices.
C     THREEPT   First derivative (non-central 3-point formula).
C
C     Environment:  FORTRAN 77 with minor extensions
C     ------------
C
C     Notes:
C     ------
C
C     (1)  IMPLICIT NONE, 8-character symbolic names, and "!" as comment
C          character are not (yet) standard.
C
C     (2)  Since many of the calculations must be repeated at both ends
C          of an interval, for each of X, Y, and Z, the various finite
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
C     (3)  Within PLSFIT3D, we have allowed for extrapolation to simplify
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
C     27 Feb. 1987    RAK    Initial design and coding of PLSFIT.
C     24 Aug. 1989  DAS/RAK  PLSFIT3D adapted from PLSFIT.
C     20 June 1991    DAS    THREEPT (monotonic) was renamed BUTLAND;
C                            THREEPT (pure 3-pt. formula) is now used
C                            for the loose fit at the ends.
C     05 Dec. 1993    DAS    ASRCH3D separated because it is now needed
C                            by PLSCURVE as well.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     &   IER, NDATA, NEVAL
      REAL
     &   TBEGIN, TEND, X (NDATA), XEVAL (NEVAL), Y (NDATA),
     &   YEVAL (NEVAL), Z (NDATA), ZEVAL (NEVAL)
      LOGICAL
     &   CLOSED, NEW
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
     &   BX (0:1), BY (0:1), BZ (0:1), CX, CY, CZ, DT, DELX (-1:1),
     &   DELY (-1:1), DELZ (-1:1), DX, DY, DZ, H (-1:1), RH, TEVAL,
     &   TINC, TLEFT, TRIGHT, TTOTAL

C     Procedures.

      REAL
     &   BESSEL, BRODLIE, BUTLAND, CHORD3D, THREEPT
      EXTERNAL
     &   BESSEL, BRODLIE, BUTLAND, CHORD3D, THREEPT, ASRCH3D

C     Storage.

      SAVE
     &   BX, BY, BZ, CX, CY, CZ, DX, DY, DZ, LEFT, RIGHT, TLEFT, TRIGHT,
     &   TTOTAL

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
         IF (X (1) .NE. X (NDATA) .OR. Y (1) .NE. Y (NDATA) .OR.
     &       Z (1) .NE. Z (NDATA)) IER = +2
      END IF

C     We'll need the total arclength to initialize the search efficiently.

      IF (NEW) THEN
         TTOTAL = CHORD3D (X, Y, Z, 1, NDATA)
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
      IF (.NOT. MONO .AND. METHOD .NE. 'B' .AND. METHOD .NE. 'L')
     &   IER = +3

      IF (NEVAL .GT. 1) THEN           ! NEVAL = 1 for T supplied externally.
         IF (DISTRIB .EQ. 'U') THEN    ! Uniform values of T generated here.
            TINC = (TEND - TBEGIN) / REAL (NEVAL - 1)
         ELSE
            IER = +4
         END IF
      END IF

C     Bail out if any error was found. (Only the last found is reported.)
C     Original test on NEVAL = 0 here was invalid - have to go through ASRCH3D.

      IF (IER .NE. 0) GO TO 990


C     Initialize bracket quantities. Note that when NEW = .TRUE. the
C     MEMORY flag is never set, thus ASRCH3D (below) always gets
C     initialized properly.

      TEVAL   = TBEGIN
      NEWFLAG = NEW

      IF (NEWFLAG) THEN
         MEMORY = .FALSE.
      ELSE
C        We can save time when PLSFIT3D is being called from within
C        a loop by setting MEMORY if possible.  The out-of-range checking
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

         IF (MEMORY) GO TO 70         ! Skip the bulk of the computation.


C        Interpolation search for bracketing interval.
C        ---------------------------------------------

         CALL ASRCH3D (NDATA, X, Y, Z, TEVAL, LEFT, TLEFT, RIGHT,
     &      TRIGHT, NEWFLAG)

C         -------------------------------------------------------------
C        |                                                             |
C        |   1 <= LEFT < RIGHT = LEFT + 1 <= NDATA                     |
C        |                                                             |
C         -------------------------------------------------------------

C        Compute derivatives by finite-differences.
C        ------------------------------------------

         IF (NDATA .GT. 2 .AND. METHOD .NE. 'L') THEN

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

            H (-1) = CHORD3D (X, Y, Z, IND (-1), IND (-1) + 1)
            H ( 0) = CHORD3D (X, Y, Z, LEFT, RIGHT)
            H (+1) = CHORD3D (X, Y, Z, IND (2) - 1, IND (2))

            DO 40, J = -1, +1
               RH       = ONE / H (J)
               DELX (J) = (X (IND (J + 1)) - X (IND (J))) * RH
               DELY (J) = (Y (IND (J + 1)) - Y (IND (J))) * RH
               DELZ (J) = (Z (IND (J + 1)) - Z (IND (J))) * RH
   40       CONTINUE

C           Select interpolation scheme.
C           ----------------------------

C           Compute adjusted first derivatives at both left- and
C           right-hand endpoints of the interval.

            IF (MONO) THEN

C              Monotone - use Brodlie modification of Butland's
C              formula to adjust the derivatives at the knots.

               DO 50, J = 0, +1
                  BX (J) = BRODLIE (J, H, DELX)
                  BY (J) = BRODLIE (J, H, DELY)
                  BZ (J) = BRODLIE (J, H, DELZ)
  50           CONTINUE

            ELSE IF (METHOD .EQ. 'B') THEN

C              Bessel - use central difference formula at the knots.

               DO 60, J = 0, +1
                  BX (J) = BESSEL (J, H, DELX)
                  BY (J) = BESSEL (J, H, DELY)
                  BZ (J) = BESSEL (J, H, DELZ)
  60           CONTINUE

            END IF      ! (METHOD .EQ. 'L') is handled with NDATA .EQ. 2

C           Patch initial/final derivatives if not periodic, case (3).

            IF (.NOT.CLOSED) THEN
               IF (LEFT .EQ. 1) THEN
                  IF (.NOT. MONO) THEN
                     BX (0) = THREEPT (0, H, DELX)
                     BY (0) = THREEPT (0, H, DELY)
                     BZ (0) = THREEPT (0, H, DELZ)
                  ELSE
                     BX (0) = BUTLAND (0, H, DELX)
                     BY (0) = BUTLAND (0, H, DELY)
                     BZ (0) = BUTLAND (0, H, DELZ)
                  END IF
               ELSE IF (RIGHT .EQ. NDATA) THEN
                  IF (.NOT. MONO) THEN
                     BX (1) = THREEPT (1, H, DELX)
                     BY (1) = THREEPT (1, H, DELY)
                     BZ (1) = THREEPT (1, H, DELZ)
                  ELSE
                     BX (1) = BUTLAND (1, H, DELX)
                     BY (1) = BUTLAND (1, H, DELY)
                     BZ (1) = BUTLAND (1, H, DELZ)
                  END IF
               END IF
            END IF

C           Compute the remaining cubic coefficients relative
C           to the left-hand endpoint.

            RH = ONE / H (0)
            CX = (THREE * DELX (0) - TWO * BX (0) - BX (1)) * RH
            CY = (THREE * DELY (0) - TWO * BY (0) - BY (1)) * RH
            CZ = (THREE * DELZ (0) - TWO * BZ (0) - BZ (1)) * RH
            DX = (-TWO * DELX (0) + BX (0) + BX (1)) * RH ** 2
            DY = (-TWO * DELY (0) + BY (0) + BY (1)) * RH ** 2
            DZ = (-TWO * DELZ (0) + BZ (0) + BZ (1)) * RH ** 2

         ELSE              ! IF (NDATA .EQ. 2 .OR. METHOD .EQ. 'L') THEN

C           Degenerate case (linear).
C           -------------------------

            H (0)  = TRIGHT - TLEFT
            RH     = ONE / H (0)
            BX (0) = (X (RIGHT) - X (LEFT)) * RH
            BY (0) = (Y (RIGHT) - Y (LEFT)) * RH
            BZ (0) = (Z (RIGHT) - Z (LEFT)) * RH
            CX     = ZERO
            CY     = ZERO
            CZ     = ZERO
            DX     = ZERO
            DY     = ZERO
            DZ     = ZERO
         END IF


   70    CONTINUE

C        Evaluate the cubics.
C        --------------------

         DT            = TEVAL - TLEFT
         XEVAL (IEVAL) = X (LEFT) + DT * (BX (0) + DT * (CX + DT * DX))
         YEVAL (IEVAL) = Y (LEFT) + DT * (BY (0) + DT * (CY + DT * DY))
         ZEVAL (IEVAL) = Z (LEFT) + DT * (BZ (0) + DT * (CZ + DT * DZ))

C        Choose next evaluation point and loop back.
C        -------------------------------------------

         IF (IEVAL .LT. NEVAL) THEN

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
