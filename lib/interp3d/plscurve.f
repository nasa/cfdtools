C+----------------------------------------------------------------------
C
      SUBROUTINE PLSCURVE (NDATA, X, Y, Z, NEW, CLOSED, T, TTOTAL,
     >                     IEVAL, XEVAL, YEVAL, ZEVAL, DERIVS)
C
C     One-liner:  Parametric local spline for a 3-space curve, with derivatives
C     ----------  (space efficient)
C
C     Description and usage:
C     ----------------------
C
C        PLSCURVE is adapted from PLSFIT3D to include return of the
C     partial derivatives of X, Y, and Z with respect to T for use in
C     curve/surface intersection calculations.  It is simpler in that
C     it handles only one value of T per call and assumes a monotonic cubic
C     fit for each of X, Y, and Z vs. T.
C
C        To match surface utilities such as PLBICUBE involving normalized
C     parametric variables U and V, PLSCURVE returns a normalized T so that
C     all of the elements of the nonlinear solution sought by an intersection
C     routine such as INTSEC4 are in the interval [0, 1].  The derivatives
C     are adjusted accordingly.
C
C        Initial scaling of the X/Y/Z data may also be advisable - see
C     PLSFIT3D for an outline.
C
C        The "NEW" argument provides for efficient repeated use on a given
C     curve dataset.  No need for working with more than one curve at a
C     time could be conceived of at the time of writing.  Thus the argument
C     list and error handling (none) are kept to a minimum.  (The periodic
C     case is not likely to be used either, but it was simpler to leave it
C     intact than to eliminate it.)
C
C     Arguments:
C     ----------
C
C     Name    Type/Dimension  I/O/S  Description
C
C     NDATA   I               I      Length of X, Y, Z input data arrays.
C                                    NDATA >= 2 (3 if CLOSED).
C
C     X,      R (NDATA)       I      Input curve coordinates. Successive
C     Y,                             data points must be distinct (non-zero
C     Z                              distance between them), but no check
C                                    is performed at this level.
C                                    If CLOSED, then first and last points
C                                    must agree (not checked here).
C
C     NEW     L               I      If control flag NEW is .TRUE., the
C                                    search for a bracket starts from
C                                    scratch, otherwise locally-saved
C                                    search and fit information including
C                                    total chord length are assumed reusable.
C
C     CLOSED  L               I      Logical flag indicating that periodic
C                                    boundary conditions are to be used.
C                                    (The curve should wrap around smoothly
C                                    on itself.)  The calling routine must
C                                    ensure that the ends agree.
C
C     T       R               I      Target value of the parametric variable
C                                    in [0, 1] (normalized chord length).
C
C     TTOTAL  R               I/O    Total chord length of curve, which may
C                                    be needed to determine the target T on
C                                    the initial call.  Therefore, enter
C                                    TTOTAL if it is known from a call to
C                                    CHORD3D or CHORDS3D, else enter TTOTAL = 0.
C                                    This is independent of the NEW usage.
C
C     IEVAL   I                 O    Points to the interval containing the
C                                    target point.  1 <= IEVAL < NDATA.
C
C     XEVAL,  R                 O    Interpolated coordinates at T.
C     YEVAL,
C     ZEVAL
C
C     DERIVS  R (3)             O    Partial derivatives of X, Y, Z with
C                                    respect to (normalized) T.
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
C     ASRCH3D   Interpolation search along chord of a parametric curve
C     BRODLIE   First derivative (central, adjusted for monotonicity)
C     BUTLAND   First derivative (non-central, adjusted for monotonicity)
C     CHORD3D   Summed chord-lengths for X-Y-Z curve over range of indices
C
C     Environment:  FORTRAN 77 with minor extensions
C     ------------
C
C     History:
C     --------
C
C     12/05/93  DAS  Initial adaptation of PLSFIT3D for curve/surface
C                    intersection calculations.  See Robert Kennelly's
C                    PLSFIT for further details about local cubic splines.
C
C     12/09/93   "   Added IEVAL argument because it's sometimes handy.
C
C     03/25/96   "   Switched from "loose" to monotonic fits to improve
C                    robustness of surface/line intersections (INTSEC4).
C
C     Author:  David Saunders, Sterling Software/NASA Ames, Mt. View, CA
C     -------
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   NDATA, IEVAL
      REAL
     >   X (NDATA), Y (NDATA), Z (NDATA), T, TTOTAL,
     >   XEVAL, YEVAL, ZEVAL, DERIVS (3)
      LOGICAL
     >   CLOSED, NEW

C     Local constants.

      REAL
     >   ZERO, ONE, TWO, THREE
      PARAMETER
     >  (ZERO  = 0.0E+0,
     >   ONE   = 1.0E+0,
     >   TWO   = 2.0E+0,
     >   THREE = 3.0E+0)

C     Local variables.

      LOGICAL
     >   MEMORY, NEWFLAG
      INTEGER
     >   IND (-1:2), J, LEFT, RIGHT
      REAL
     >   BX (0:1), BY (0:1), BZ (0:1), CX, CY, CZ, DT, DELX (-1:1),
     >   DELY (-1:1), DELZ (-1:1), DX, DY, DZ, H (-1:1), RH, TEVAL,
     >   TLEFT, TRIGHT

C     Procedures.

      REAL
     >   BRODLIE, BUTLAND, CHORD3D
      EXTERNAL
     >   BRODLIE, BUTLAND, CHORD3D, ASRCH3D

C     Storage.

      SAVE
     >   BX, BY, BZ, CX, CY, CZ, DX, DY, DZ, LEFT, RIGHT, TLEFT, TRIGHT


C     Initialization.
C     ---------------

C     We'll need the total arc-length to initialize the search efficiently.

      IF (NEW) THEN
         IF (TTOTAL .EQ. ZERO) TTOTAL = CHORD3D (X, Y, Z, 1, NDATA)
         TRIGHT = TTOTAL
      END IF

      TEVAL = T * TTOTAL  ! Since input T is normalized by total length

C     Initialize bracket quantities. Note that when NEW = .TRUE. the
C     MEMORY flag is never set, thus ASRCH3D (below) always gets
C     initialized properly.

      NEWFLAG = NEW

      IF (NEWFLAG) THEN
         MEMORY = .FALSE.
      ELSE
C        We can save time when PLSCURVE is being called from within
C        a loop by setting MEMORY if possible.  The out-of-range checking
C        relies on the fact that RIGHT = LEFT + 1 after the last return.
C        Cater to the common case of TEVAL in the previous, interior
C        interval.

         MEMORY = (TEVAL .GE. TLEFT) .AND. (TEVAL .LT. TRIGHT)
         IF (.NOT. MEMORY) THEN
            MEMORY = (LEFT  .EQ. 1)     .AND. (TEVAL .LT. TRIGHT) .OR.
     >               (RIGHT .EQ. NDATA) .AND. (TEVAL .GE. TLEFT)
         END IF
      END IF


      IF (.NOT. MEMORY) THEN

C        Interpolation search for bracketing interval.
C        ---------------------------------------------

         CALL ASRCH3D (NDATA, X, Y, Z, TEVAL, LEFT, TLEFT, RIGHT,
     >      TRIGHT, NEWFLAG)

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

            H (-1) = CHORD3D (X, Y, Z, IND (-1), IND (-1) + 1)
            H ( 0) = CHORD3D (X, Y, Z, LEFT, RIGHT)
            H (+1) = CHORD3D (X, Y, Z, IND (2) - 1, IND (2))

            DO 40, J = -1, +1
               RH       = ONE / H (J)
               DELX (J) = (X (IND (J + 1)) - X (IND (J))) * RH
               DELY (J) = (Y (IND (J + 1)) - Y (IND (J))) * RH
               DELZ (J) = (Z (IND (J + 1)) - Z (IND (J))) * RH
   40       CONTINUE

C           Hermite-type interpolation scheme.
C           ----------------------------------

C           Compute adjusted first derivatives at both left- and
C           right-hand endpoints of the interval.

            DO 60, J = 0, +1
               BX (J) = BRODLIE (J, H, DELX)
               BY (J) = BRODLIE (J, H, DELY)
               BZ (J) = BRODLIE (J, H, DELZ)
  60        CONTINUE

C           Patch initial/final derivatives if not periodic, case (3).

            IF (.NOT. CLOSED) THEN
               IF (LEFT .EQ. 1) THEN
                  BX (0) = BUTLAND (0, H, DELX)
                  BY (0) = BUTLAND (0, H, DELY)
                  BZ (0) = BUTLAND (0, H, DELZ)
               ELSE IF (RIGHT .EQ. NDATA) THEN
                  BX (1) = BUTLAND (1, H, DELX)
                  BY (1) = BUTLAND (1, H, DELY)
                  BZ (1) = BUTLAND (1, H, DELZ)
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

         ELSE              ! IF (NDATA .EQ. 2) THEN

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

      END IF


      IEVAL = LEFT

C     Evaluate the cubics for X, Y, and Z.
C     ------------------------------------

      DT    = TEVAL - TLEFT
      XEVAL = X (LEFT) + DT * (BX (0) + DT * (CX + DT * DX))
      YEVAL = Y (LEFT) + DT * (BY (0) + DT * (CY + DT * DY))
      ZEVAL = Z (LEFT) + DT * (BZ (0) + DT * (CZ + DT * DZ))


C     Evaluate the derivatives w.r.t. T = TEVAL / TTOTAL.
C     ---------------------------------------------------

      DERIVS (1) = (BX (0) + DT * (TWO * CX + DT * THREE * DX)) * TTOTAL
      DERIVS (2) = (BY (0) + DT * (TWO * CY + DT * THREE * DY)) * TTOTAL
      DERIVS (3) = (BZ (0) + DT * (TWO * CZ + DT * THREE * DZ)) * TTOTAL

C     Termination.
C     ------------

      RETURN
      END
