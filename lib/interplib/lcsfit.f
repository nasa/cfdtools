C+----------------------------------------------------------------------
C
      SUBROUTINE LCSFIT (NDATA, X, Y, NEW, METHOD, NEVAL, XEVAL, YEVAL,
     &   YPEVAL)
C
C     Two-liner:  Storage-efficient local cubic spline fit (2-space)
C     ----------  (monotonic and piecewise linear options too)
C
C     Description and usage:
C     ----------------------
C
C        LCSFIT is the non-parametric analog of PLSFIT (parametric).
C     It is intended for spline applications which do not require the
C     spline coefficients as output.  It is efficient for repeated
C     calls with the same data, so repeated use with NEVAL = 1 may be
C     preferable to storing vectors of results.
C
C        LCSFIT offers monotonic spline and piecewise linear options
C     also.  And it returns an interpolated first derivative along
C     with the function value.  (The second derivative is omitted
C     because Y" is not guaranteed to be continuous by local methods.)
C
C        See PLSFIT for more details on local methods.  As with most
C     numerical methods, scaling of the data to the unit interval (and
C     unscaling of the result) is recommended to avoid unnecessary
C     effects attributable to the data units.  Utilities GETSCALE and
C     USESCALE from the present authors are appropriate.  The data
C     abscissas should be distinct and either ascending or descending.
C     PROTECT is available to check this.  Extrapolation is permitted
C     (mainly in case of round-off; it is normally inadvisable).
C
C        The CSFIT/CSEVAL or CSDVAL pair are probably preferable if
C     efficiency is not an issue, since CSFIT gives Y" continuity.
C
C     Arguments:
C     ----------
C
C     Name    Type/Dimension  I/O/S  Description
C     NDATA   I               I      Length of X, Y input data arrays.
C
C     X,      R (NDATA)       I      Input data coordinates.  The Xs
C     Y                              must be distinct and monotonic,
C                                    either ascending or descending.
C                                    (No check here.) 
C
C     NEW     L               I      If control flag NEW is .TRUE., the
C                                    search for a bracket starts from
C                                    scratch, otherwise locally-saved
C                                    search and fit information will be
C                                    assumed to be correct. If calling
C                                    LCSFIT from within a loop, set
C                                    NEW = .FALSE. after the first call.
C
C     METHOD   C*1            I      (Uppercase) Type of fit to be used:
C                                    'M' means Monotonic piecewise cubics;
C                                    'B' means non-monotonic "Bessel"-type
C                                        piecewise cubics (looser fit);
C                                    'L' means piecewise Linear fit;
C                                    'C' means Cyclic (periodic) end
C                                        conditions: loose fit assumed.
C
C     NEVAL   I               I      Number of interpolations requested.
C                                    NEVAL >= 1.  One call per result
C                                    (NEVAL = 1) may save storage, and is
C                                    not too inefficient as long as NEW
C                                    is set to .FALSE. after the first.
C
C     XEVAL   R (NEVAL)       I      Abscissa(s) to interpolate to.  These
C                                    are normally in the data range, but
C                                    extrapolation - probably due to
C                                    round-off - is not prevented.
C
C     YEVAL   R (NEVAL)       O      Interpolated function value(s).
C
C     YPEVAL  R (NEVAL)       O      Interpolated 1st derivative value(s).
C                                    Pass the same storage as for YEVAL
C                                    if no derivatives are required.
C
C     Significant local variables:
C     ----------------------------
C
C     MEMORY         Indicates that coefficients are correct for the
C                    current point.
C
C     H, DEL         Delta X and forward difference derivative arrays.
C
C     B, C, D        Coefficients of cubic on the bracketing interval.
C
C     Procedures:
C     -----------
C
C     INTERVAL  1-D "interpolation" search.
C     BESSEL    First derivative (central 3-point formula).
C     BRODLIE   First derivative (central), adjusted for monotonicity.
C     BUTLAND   First derivative (non-central), adjusted for monotonicity.
C     THREEPT   First derivative (non-central 3-point formula).
C
C     Environment:  FORTRAN 90
C     ------------
C
C     Error handling:  None
C     ---------------
C
C     Notes:
C     ------
C
C     (1)  Since many of the calculations must be repeated at both ends
C          of an interval, the various finite difference quantities used
C          are stored as arrays. The following "map" of a typical interior
C          interval and its neighbors should help in understanding the
C          notation.  The local array indices are all numbered relative
C          to the left-hand end of the interval which brackets the point
C          to be evaluated.
C
C                                  LEFT       RIGHT
C
C          Point         -1          0         +1          +2
C
C          Data           x -------- x -------- x --------- x
C
C          Interval      -1          0         +1
C
C
C     Author: Robert Kennelly, Sterling Software/NASA Ames  (PLSFIT)
C     -------
C
C     History:
C     --------
C
C     27 Feb. 1987  R.A.Kennelly  Initial implementation of PLSFIT.
C     23 Aug. 1989  D.A.Saunders  LCSFIT adapted as non-parametric form,
C                                 for embedding in other utilities where
C                                 minimizing work-space is desirable.
C     20 June 1991    "    "      THREEPT (monotonic) renamed BUTLAND;
C                                 THREEPT (pure 3-pt. formula) now used
C                                 for nonmonotonic end-point handling;
C                                 METHOD='C' case belatedly added, as
C                                 needed by PLSINTRP for closed curves.
C     23 July 1991    "    "      The tests for being in the same interval
C                                 as before were not allowing for the
C                                 descending-Xs case.
C     06 May  1998    "    "      Minor Fortran 90 updates.
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER,   INTENT (IN)  :: NDATA, NEVAL
      REAL,      INTENT (IN)  :: X (NDATA), Y (NDATA), XEVAL (NEVAL)
      REAL,      INTENT (OUT) :: YEVAL (NEVAL), YPEVAL (NEVAL)
      LOGICAL,   INTENT (IN)  :: NEW
      CHARACTER, INTENT (IN)  :: METHOD * 1

C     Local constants:

      REAL, PARAMETER :: ZERO = 0., ONE = 1., TWO = 2., THREE = 3.

C     Local variables:

      LOGICAL
     &   CYCLIC, LINEAR, MEMORY, MONO
      INTEGER
     &   IEVAL, J, K, LEFT, RIGHT
      REAL
     &   ARROW, B (0:1), C, DELY (-1:1), D, DX, H (-1:1), XBYARROW, XE

C     Procedures:

      REAL, EXTERNAL ::
     &   BESSEL, BRODLIE, BUTLAND, THREEPT

C     Storage:

      SAVE
     &   ARROW, B, C, D, LEFT, RIGHT

C     Execution:
C     ----------

      MONO   = METHOD == 'M'
      CYCLIC = METHOD == 'C'
      LINEAR = METHOD == 'L'

      IF (CYCLIC) THEN
         IF (Y (NDATA) /= Y (1)) STOP 'LCSFIT: End points must match.'
      END IF

C     Initialize search or avoid it if possible:

      IF (NEW) THEN
         MEMORY = .FALSE.
         ARROW  = SIGN (ONE, X (2) - X (1))
         LEFT   = 1
      END IF

      IEVAL = 1
      XE = XEVAL (1)
      XBYARROW = XE * ARROW

      IF (.NOT. NEW) THEN
      
C        We can save a lot of time when LCSFIT is being called from within
C        a loop by setting MEMORY if possible. The out-of-range checking
C        relies on the fact that RIGHT = LEFT + 1 after the last return.
C        Cater to the more likely case of XE in the previous, interior
C        interval.

         MEMORY = XBYARROW >= X (LEFT)  * ARROW .AND.
     &            XBYARROW <  X (RIGHT) * ARROW

         IF (.NOT. MEMORY) THEN
            MEMORY =
     &         LEFT  == 1     .AND. XBYARROW <  X (RIGHT) * ARROW
     &         .OR.
     &         RIGHT == NDATA .AND. XBYARROW >= X (LEFT)  * ARROW
         END IF

      END IF

      IF (MEMORY) GO TO 70 ! Skip the bulk of the computation

C     Loop over evaluation points requiring a new search:
C     ---------------------------------------------------

   10 CONTINUE

C        Interpolation search for bracketing interval:
C        ---------------------------------------------

         CALL INTERVAL (NDATA, X, XE, ARROW, LEFT)

         RIGHT = LEFT + 1

C         -------------------------------------------
C        |                                           |
C        |   1 <= LEFT < RIGHT = LEFT + 1 <= NDATA   |
C        |                                           |
C         -------------------------------------------

C        Compute derivatives by finite-differences:
C        ------------------------------------------

         IF (NDATA > 2 .AND. .NOT. LINEAR) THEN

C           Interval and derivative approximations:
C           ---------------------------------------

C           The following duplicates more code than PLSFIT's approach,
C           but eliminates some indirection - no need to wrap-around here.
C           Handle the end conditions first to minimize testing LEFT, RIGHT.

            IF (LEFT == 1) THEN

               H (0) = X (2) - X (1)
               DELY (0) = (Y (2) - Y (1)) / H (0)
               H (1) = X (3) - X (2)
               DELY (1) = (Y (3) - Y (2)) / H (1)

               IF (CYCLIC) THEN ! Loose fit assumed
                  H (-1) = X (NDATA) - X (NDATA - 1)
                  DELY (-1) = (Y (NDATA) - Y (NDATA - 1)) / H (-1)
                  B (0) = BESSEL (0, H, DELY)
                  B (1) = BESSEL (1, H, DELY)
               ELSE
                  IF (MONO) THEN
                     B (0) = BUTLAND (0, H, DELY)
                     B (1) = BRODLIE (1, H, DELY)
                  ELSE
                     B (0) = THREEPT (0, H, DELY)
                     B (1) = BESSEL  (1, H, DELY)
                  END IF
               END IF

            ELSE IF (RIGHT == NDATA) THEN

               H(-1) = X (LEFT) - X (LEFT-1)
               DELY(-1) = (Y (LEFT) - Y (LEFT-1)) / H (-1)
               H (0) = X (RIGHT) - X (LEFT)
               DELY (0) = (Y (RIGHT) - Y (LEFT))  / H (0)

               IF (CYCLIC) THEN
                  H (1) = X (2) - X (1)
                  DELY (1) = (Y (2) - Y (1)) / H (1)
                  B (0) = BESSEL (0, H, DELY)
                  B (1) = BESSEL (1, H, DELY)
               ELSE

                  IF (MONO) THEN
                     B (0) = BRODLIE (0, H, DELY)
                     B (1) = BUTLAND (1, H, DELY)
                  ELSE
                     B (0) = BESSEL  (0, H, DELY)
                     B (1) = THREEPT (1, H, DELY)
                  END IF
               END IF

            ELSE

               K = LEFT
               DO J = -1, +1
                  H (J)    =  X (K) - X (K-1)
                  DELY (J) = (Y (K) - Y (K-1)) / H (J)
                  K = K + 1
               END DO

C              Select interpolation scheme:
C              ----------------------------

C              Compute (possibly adjusted) first derivatives at both
C              left- and right-hand endpoints of the interval.

               IF (MONO) THEN

C                 Monotone - use Brodlie modification of Butland's
C                 formula to adjust the derivatives at the knots.

                  B (0) = BRODLIE (0, H, DELY)
                  B (1) = BRODLIE (1, H, DELY)

               ELSE ! METHOD = 'B'

C                 Bessel - use central difference formula at the knots.

                  B (0) = BESSEL (0, H, DELY)
                  B (1) = BESSEL (1, H, DELY)

               END IF

            END IF

C           Compute the remaining cubic coefficients.

            C = (THREE * DELY (0) - TWO * B (0) - B (1)) / H (0)
            D = ( -TWO * DELY (0) + B (0) + B (1)) / H (0) ** 2

         ELSE ! NDATA = 2 .OR. METHOD = 'L'

C           Degenerate case (linear).
C           -------------------------

            B (0) = (Y (RIGHT) - Y (LEFT)) / (X (RIGHT) - X (LEFT))
            C     = ZERO
            D     = ZERO

         END IF

C        Evaluate the cubic (derivative first in case only YEVAL is reqd.):
C        ------------------------------------------------------------------

   70    CONTINUE ! Start of same-interval loop inside new-interval loop

            DX = XE - X (LEFT)
            YPEVAL (IEVAL) = B (0) + DX * (TWO * C + DX * THREE * D)
            YEVAL  (IEVAL) = Y (LEFT) + DX * (B (0) + DX * (C + DX * D))

C           The next evaluation point may be in the same interval:
C           ------------------------------------------------------

            IF (IEVAL < NEVAL) THEN ! Skips this if NEVAL = 1

               IEVAL = IEVAL + 1
               XE = XEVAL (IEVAL)
               XBYARROW = XE * ARROW
               IF (XBYARROW >= X (LEFT)  * ARROW  .AND.
     &             XBYARROW <  X (RIGHT) * ARROW) GO TO 70
            
               GO TO 10 ! Else much more work required.

            END IF

C     Termination:
C     ------------

      RETURN

      END SUBROUTINE LCSFIT
