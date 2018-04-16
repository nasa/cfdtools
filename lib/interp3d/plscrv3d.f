C+------------------------------------------------------------------------------
C
      SUBROUTINE PLSCRV3D (NDATA, X, Y, Z, T, METHOD, NEW, CLOSED,
     >                     TEVAL, IEVAL, XEVAL, YEVAL, ZEVAL, DERIVS)
C
C     Description:
C
C        PLSCRV3D interpolates X, Y, Z along a curve at a specified TEVAL using
C     piecewise linear or local cubic spline techniques.  It is an adaptation of
C     the earlier PLSFIT3D and PLSCURVE intended to avoid internal arc length
C     calculations by assuming that the cumulative chord lengths associated with
C     the data points are supplied by the application. (PLSFIT3D avoided storing
C     all the Ts as part of its space-efficient graphics-oriented approach, but
C     in practical grid generation applications the Ts are needed anyway.)
C
C        PLSCRV3D restores the METHOD options that were dropped by PLSCURVE,
C     in anticipation of a requirement for piecewise linear results.  The option
C     to return derivatives for intersection calculations influenced the choice
C     between one and many target Ts: a single T wins.  Use of normalized Ts is
C     recommended for intersection calculations, so T is known to be in [0, 1],
C     but other applications may prefer not to normalize, and this is permitted.
C     T(*) may also be in descending order.
C
C        Arithmetic for unwanted derivatives is avoided by calls with input
C     DERIVS(1) = -999.
C
C        See Robert Kennelly's original 2-space PLSFIT for more on compact local
C     spline methods.
C
C     Environment:  Fortran 90
C
C     History:
C
C     09/22/97  DAS  Initial adaptation of PLSFIT3D/PLSCURVE with more efficient
C                    curve/surface intersections and wing surface grids in mind,
C                    especially NDATA = 2 and piecewise linear cases.
C     11/17/97   "   Allowed for non-normalized input Ts.
C     05/06/98   "   Minor Fortran 90 updates.
C     07/25/06   "   Allowed for descending T(*) order.
C
C     Author:  David Saunders, Sterling Software and
C                              ELORET/NASA Ames Research Ctr., Moffett Field, CA
C
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   NDATA               ! Length of X, Y, Z input data arrays.

      REAL, INTENT (IN) ::
     >   X(NDATA),           ! Input curve coordinates. Successive
     >   Y(NDATA),           ! data points must be distinct.
     >   Z(NDATA)            ! If CLOSED, then first and last points
                             ! must agree (not checked here).
      REAL, INTENT (IN) ::
     >   T(NDATA)            ! Parametric variable values associated
                             ! with the data points, usually normalized
                             ! cumulative chord lengths, but they need
                             ! not be normalized or increasing except for
                             ! intersection calculations via INTSEC4 or -5.
      CHARACTER, INTENT (IN) ::
     >   METHOD*1            ! The type of parametric fit to be used:
                             !    'M' means monotonic piecewise cubics;
                             !    'B' means non-monotonic "Bessel"-type
                             !        piecewise cubics (looser fit);
                             !    'L' means piecewise linear (only way
                             !        to guarantee no overshoots).
      LOGICAL, INTENT (IN) ::
     >   NEW                 ! If NEW = T, the search for a bracket starts
                             ! from scratch, otherwise locally-saved search
                             ! and fit information will be assumed to be
                             ! correct. If calling PLSCRV3D from within a
                             ! loop, set NEW = F after the first call.
      LOGICAL, INTENT (IN) ::
     >   CLOSED              ! CLOSED = T means periodic boundary conditions
                             ! are to be used. The curve should wrap around
                             ! smoothly on itself. The calling routine must
                             ! ensure that the end points agree.
      REAL, INTENT (IN) ::
     >   TEVAL               ! Target value of T at which to evaluate X,Y,Z.

      INTEGER, INTENT (INOUT) ::
     >   IEVAL               ! Input with the data index at which to start
                             ! searching for TEVAL; output with actual index.

      REAL, INTENT (OUT) ::
     >   XEVAL, YEVAL, ZEVAL ! Interpolated coordinates.

      REAL, INTENT (INOUT) ::
     >   DERIVS(3)           ! Partial derivatives of X, Y, Z with respect to T.
                             ! Enter DERIVS(1) = -999. to suppress these,
                             ! else  DERIVS(1) = 0. (say) on the first call
                             ! and something other than -999. thereafter.
C     Procedures:

      REAL, EXTERNAL ::
     >   BESSEL,             ! 1st deriv. (central 3-point formula)
     >   BRODLIE,            ! 1st deriv. (central) adjusted for monotonicity
     >   BUTLAND,            ! 1st deriv. (non-central)  "    "    "
     >   THREEPT             ! 1st deriv. (non-central 3-point formula)
      EXTERNAL
     >   INTERVAL            ! Interpolatory 1-D search utility
C-------------------------------------------------------------------------------

C     Local constants:

      REAL, PARAMETER :: ONE = 1., TWO = 2., THREE = 3.

C     Local variables:

      LOGICAL    DERIV, LINEAR, LOOSE, MEMORY, TWOPTS
      INTEGER    IND (-1:2), J, LEFT, RIGHT
      REAL       ARROW, BX (0:1), BY (0:1), BZ (0:1), CX, CY, CZ, DT,
     >           DELX (-1:1), DELY (-1:1), DELZ (-1:1), DX, DY, DZ,
     >           H (-1:1), RH, TBYARROW, TLEFT, TRIGHT

C     Storage:

      SAVE       ARROW, BX, BY, BZ, CX, CY, CZ, DX, DY, DZ, LEFT, RIGHT,
     >           TLEFT, TRIGHT, DERIV, LINEAR, LOOSE, TWOPTS

C     Execution:

C     The 2-data-point case is common enough for intersections that it
C     is treated specially.  Minor code duplication is OK.

      IF (NEW) THEN

         DERIV  = DERIVS (1) /= -999.
         TWOPTS = NDATA == 2

         IF (TWOPTS) THEN

            RH     = ONE / (T (2) - T (1))
            BX (0) = (X (2) - X (1)) * RH
            BY (0) = (Y (2) - Y (1)) * RH
            BZ (0) = (Z (2) - Z (1)) * RH

            DT    = TEVAL - T (1)
            XEVAL = X (1) + DT * BX (0)
            YEVAL = Y (1) + DT * BY (0)
            ZEVAL = Z (1) + DT * BZ (0)

            IF (DERIV) THEN
               DERIVS (1) = BX (0)
               DERIVS (2) = BY (0)
               DERIVS (3) = BZ (0)
            END IF

            GO TO 999

         ELSE
            MEMORY = .FALSE.
            LINEAR = METHOD == 'L'
            LOOSE  = METHOD == 'B' ! Else monotonic
            ARROW  = SIGN (ONE, T (2) - T (1))
         END IF

      ELSE IF (TWOPTS) THEN

         DT    = TEVAL - T (1)
         XEVAL = X (1) + DT * BX (0)
         YEVAL = Y (1) + DT * BY (0)
         ZEVAL = Z (1) + DT * BZ (0)

         IF (DERIV) THEN  ! INTSEC5 overwrites previous derivatives
            DERIVS (1) = BX (0)
            DERIVS (2) = BY (0)
            DERIVS (3) = BZ (0)
         END IF

         GO TO 999

      ELSE ! NEW = F and NDATA > 2

C        We can save time when PLSCRV3D is being called from within a loop
C        by setting MEMORY if possible.  The out-of-range checking relies on
C        the fact that RIGHT = LEFT + 1 after the last return.  Cater to the
C        more likely case of TEVAL in the previous (same) interior interval.

         TBYARROW = TEVAL * ARROW
         MEMORY   = TBYARROW >= TLEFT * ARROW .AND.
     >              TBYARROW < TRIGHT * ARROW

         IF (.NOT. MEMORY) THEN
            MEMORY = (LEFT  == 1)     .AND. (TBYARROW <  TRIGHT * ARROW)
     >               .OR.
     >               (RIGHT == NDATA) .AND. (TBYARROW >= TLEFT  * ARROW)
         END IF

      END IF

      IF (.NOT. MEMORY) THEN

C        Interpolation search for bracketing interval:

         LEFT = IEVAL

         CALL INTERVAL (NDATA, T, TEVAL, ARROW, LEFT)

         IEVAL  = LEFT
         TLEFT  = T (LEFT)
         RIGHT  = LEFT + 1
         TRIGHT = T (RIGHT)

C         -------------------------------------------
C        |   1 <= LEFT < RIGHT = LEFT + 1 <= NDATA   |
C         -------------------------------------------

         IF (LINEAR) THEN ! Piecewise linear, NDATA > 2

            RH     = ONE / (TRIGHT - TLEFT)
            BX (0) = (X (RIGHT) - X (LEFT)) * RH
            BY (0) = (Y (RIGHT) - Y (LEFT)) * RH
            BZ (0) = (Z (RIGHT) - Z (LEFT)) * RH

         ELSE ! Piecewise cubic

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
C           case will be overridden again, following derivative calculation.)

            IF (LEFT == 1) THEN
               IND (-1) = NDATA - 1 ! Left side wrap-around boundary condition
            ELSE IF (RIGHT == NDATA) THEN
               IND (2) = 2          ! Right side
            END IF

C           Interval and derivative approximations:

            H (-1) = T (IND (-1) + 1) - T (IND (-1))
            H ( 0) = T (RIGHT) - T (LEFT)
            H (+1) = T (IND (2)) - T (IND (2) - 1)

            DO J = -1, +1
               RH       = ONE / H (J)
               DELX (J) = (X (IND (J + 1)) - X (IND (J))) * RH
               DELY (J) = (Y (IND (J + 1)) - Y (IND (J))) * RH
               DELZ (J) = (Z (IND (J + 1)) - Z (IND (J))) * RH
            END DO

C           Select the interpolation scheme, and compute [adjusted] first
C           derivatives at both left- and right-hand endpoints of the interval:

            IF (LOOSE) THEN ! Bessel - use the central formula at 0, +1

               BX (0) = BESSEL (0, H, DELX)
               BX (1) = BESSEL (1, H, DELX)
               BY (0) = BESSEL (0, H, DELY)
               BY (1) = BESSEL (1, H, DELY)
               BZ (0) = BESSEL (0, H, DELZ)
               BZ (1) = BESSEL (1, H, DELZ)

            ELSE ! Monotone - use the Brodlie modification of Butland's formula

               BX (0) = BRODLIE (0, H, DELX)
               BX (1) = BRODLIE (1, H, DELX)
               BY (0) = BRODLIE (0, H, DELY)
               BY (1) = BRODLIE (1, H, DELY)
               BZ (0) = BRODLIE (0, H, DELZ)
               BZ (1) = BRODLIE (1, H, DELZ)

            END IF

C           Patch initial/final derivatives if not periodic, case (3):

            IF (.NOT. CLOSED) THEN
               IF (LEFT == 1) THEN
                  IF (LOOSE) THEN
                     BX (0) = THREEPT (0, H, DELX)
                     BY (0) = THREEPT (0, H, DELY)
                     BZ (0) = THREEPT (0, H, DELZ)
                  ELSE
                     BX (0) = BUTLAND (0, H, DELX)
                     BY (0) = BUTLAND (0, H, DELY)
                     BZ (0) = BUTLAND (0, H, DELZ)
                  END IF
               ELSE IF (RIGHT == NDATA) THEN
                  IF (LOOSE) THEN
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

         END IF

      END IF

C     Evaluate the identified polynomials:

      DT = TEVAL - TLEFT

      IF (LINEAR) THEN

         XEVAL = X (LEFT) + DT * BX (0)
         YEVAL = Y (LEFT) + DT * BY (0)
         ZEVAL = Z (LEFT) + DT * BZ (0)

         IF (DERIV) THEN
            DERIVS (1) = BX (0)
            DERIVS (2) = BY (0)
            DERIVS (3) = BZ (0)
         END IF

      ELSE ! Parametric cubics

         XEVAL = X (LEFT) + DT * (BX (0) + DT * (CX + DT * DX))
         YEVAL = Y (LEFT) + DT * (BY (0) + DT * (CY + DT * DY))
         ZEVAL = Z (LEFT) + DT * (BZ (0) + DT * (CZ + DT * DZ))

         IF (DERIV) THEN
            DERIVS (1) = BX (0) + DT * (TWO * CX + DT * THREE * DX)
            DERIVS (2) = BY (0) + DT * (TWO * CY + DT * THREE * DY)
            DERIVS (3) = BZ (0) + DT * (TWO * CZ + DT * THREE * DZ)
         END IF

      END IF

  999 RETURN

      END SUBROUTINE PLSCRV3D
