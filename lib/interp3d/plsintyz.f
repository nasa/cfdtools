C+----------------------------------------------------------------------
C
      SUBROUTINE PLSINTYZ (NDATA, X, Y, Z, I1, I2, METHOD, NEVAL, XEVAL,
     >                     YEVAL, ZEVAL)
C
C     Two-liner:  Storage-efficient parametric local cubic spline
C     ----------  interpolation in 3-space of Y and Z for given X
C
C     Description and usage:
C     ----------------------
C
C        PLSINTYZ is the (X,Y,Z)-space version of PLSINTRP ((X,Y)-space.)
C     PLSINTYZ determines Y and Z values along a curve corresponding to
C     given Xs which must be confined to a monotonic sub-section of the
C     curve.  It applies the local spline technique of LCSFIT, which
C     avoids storage of spline coefficients - they are computed as needed.
C
C        The method involves identifying the value of T (within the
C     monotonic sub-curve defined by subscripts I1, I2) which produces
C     the target X, then evaluating Y and Z for that T.  Locating T is a
C     zero-finding problem, done here with a Newton iteration rather
C     than with an external utility, to avoid the communication problems
C     (via COMMON blocks) typical of such utilities which expect a
C     single-argument routine defining the function to be zeroed.
C
C        As with LCSFIT, a piecewise linear option is included.  The
C     interpolated first derivatives of PLSINTRP are omitted here.
C     Remember that local methods guarantee only first-derivative
C     continuity: there is no control over second derivative.
C
C        See PLSFIT for more details on local methods.  As with most
C     numerical methods, scaling of the data to the unit interval (and
C     unscaling of the result) is recommended to avoid unnecessary
C     effects attributable to the data units.  Utilities GETSCALE and
C     USESCALE from the present authors are appropriate.  Extrapolation
C     is permitted (mainly in case of round-off; it is normally not
C     advisable).
C
C     Arguments:
C     ----------
C
C     Name    Type/Dimension  I/O/S  Description
C     NDATA   I               I      Length of X, Y input data arrays.
C
C     X,      R (NDATA)       I      Input data coordinates.  The Xs
C     Y                              must be distinct but need not be
C     Z                              monotonic except from I1 to I2,
C                                    where they may be either ascending
C                                    or descending.
C
C     I1,    I                I      Data indices defining a sub-curve
C     I2                             in X, Y, Z (*) where X is monotonic.
C
C     METHOD   C*1            I      (Uppercase) Controls fits to be used:
C                                    'M' means Monotonic piecewise cubics;
C                                    'B' means non-monotonic "Bessel"-type
C                                        piecewise cubics (looser fit);
C                                    'L' means piecewise Linear fit;
C                                    'C' means the curve is to be Closed
C                                        smoothly: points 1 and NDATA must
C                                        match, and loose fits are used.
C
C     NEVAL   I               I      Number of interpolations requested.
C                                    NEVAL >= 1.
C
C     XEVAL   R (NEVAL)       I      Abscissa(s) to interpolate to.  These
C                                    are normally in the data range, but
C                                    extrapolation - probably due to
C                                    round-off - is not prevented.
C
C     YEVAL   R (NEVAL)       O      Interpolated Y value(s).
C
C     ZEVAL   R (NEVAL)       O      Interpolated Z value(s).
C
C     Procedures:
C     -----------
C
C     CHORD3D   Safeguarded chord-length function
C     INTERVAL  1-D "interpolation" search
C     LCSFIT    Local cubic spline method for monotonic data
C
C
C     Environment:  VAX/VMS; FORTRAN 77
C     ------------
C
C     IMPLICIT NONE, 8-character symbolic names, and "!" as comment
C     character are not (yet) standard.
C
C
C     Author: David Saunders, NASA Ames/Sterling Software, Palo Alto, CA.
C     -------
C
C     History:
C     --------
C
C     06/27/91  D.A.Saunders  PLSINTYZ adapted from PLSINTRP.  Debug
C                             printout to unit 6 if convergence fails
C                             has been retained: it should point to
C                             programmer error in the application.
C     07/15/91    "    "      Changed TOL to 6.E-7 as for PLSINTRP.
C     07/24/91    "    "      The test for the same interval was over-
C                             looking the descending Xs case.  And the
C                             Newton iteration can take advantage of
C                             LCSFIT's "NEW" argument.
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   I1, I2, NDATA, NEVAL
      REAL
     >   X (NDATA), XEVAL (NEVAL), Y (NDATA), YEVAL (NEVAL),
     >   Z (NDATA), ZEVAL (NEVAL)
      CHARACTER
     >   METHOD * 1

C     Local constants.

      INTEGER
     >   MXITER
      REAL
     >   EPS, ONE, TOL, ZERO
      PARAMETER
     >  (MXITER = 20,
     >   EPS    = 1.E-30,    ! Saveguards dX / dT
     >   ONE    = 1.E+0,
     >   TOL    = 6.E-7,     ! Relative tolerance for X by Newton iteration
     >   ZERO   = 0.0E+0)    ! (suited to VAX single precision)

C     Local variables.

      INTEGER
     >   IEVAL, ITER, J1, J2, LEFT, N, NT
      REAL
     >   ARROW, T (4), TI, TMAX, TMIN, TOLER, XE, XI, XPI, YI
      LOGICAL
     >   ENDPT, NEW

C     Procedures.

      REAL
     >   CHORD3D
      EXTERNAL
     >   CHORD3D, INTERVAL, LCSFIT

C     Execution.
C     ----------

      T (1) = ZERO
      IEVAL = 1
      XE = XEVAL (1)
      ARROW = SIGN (ONE, X (I1 + 1) - X (I1))
      N = I2 - I1 + 1

C     Loop over evaluation points.
C     ----------------------------

   10 CONTINUE

C        End points can cause trouble if dX / dT ~ 0.  Avoid them:

         IF (XE .EQ. X (I1)) THEN
            ENDPT = .TRUE.
            YEVAL (IEVAL) = Y (I1)
            ZEVAL (IEVAL) = Z (I1)

         ELSE IF (XE .EQ. X (I2)) THEN
            ENDPT = .TRUE.
            YEVAL (IEVAL) = Y (I2)
            ZEVAL (IEVAL) = Z (I2)

         ELSE
            ENDPT = .FALSE.

C           Interpolation search for interval containing target X:

            CALL INTERVAL (N, X (I1), XE, ARROW, LEFT)

            LEFT = LEFT + I1 - 1          ! I1 <= LEFT < LEFT + 1 <= I2

C           LCSFIT requires only the neighboring 4 points (3 at a boundary).

            J1 = MAX (LEFT - 1, 1)        ! These MAY go outside I1, I2
            J2 = MIN (LEFT + 2, NDATA)    ! for curve fit (but not TI).
            NT = J2 - J1 + 1

            T (2) = CHORD3D (X, Y, Z, J1, J1 + 1)
            T (3) = CHORD3D (X, Y, Z, J1 + 1, J1 + 2) + T (2)
            T (4) = CHORD3D (X, Y, Z, J2 - 1, J2) + T (3)    ! Not used at ends

C           Nasty potential for getting (X,T) on the wrong sub-curve unless...

            TMIN = ZERO                   ! Normally
            TMAX = T (NT)
            IF (J2 .GT. I2) THEN          ! Points in use have wrapped;
               TMAX = T (3) - ABS (X (I2) - XE)     ! dX < dT; avoid I2.
            ELSE IF (J1 .LT. I1) THEN
               TMIN = T (2) + ABS (X (I1) - XE)     ! Avoid I1
            END IF

C           Newton iteration for T that corresponds to the target X.
C           Find where F (T) = X (T) - XE = 0.   Note that F'(T) = X'(T).

            TI = MAX (T (2), TMIN)  ! T (2) corresponds to X (LEFT) except
                                    ! when LEFT = 1.  Can't exceed TMAX.

   20       ITER = 0
            TOLER = (ONE + ABS (XE)) * TOL  ! Relative test, & protects XE = 0.
            NEW = .TRUE.

   30       CONTINUE
               ITER = ITER + 1

               CALL LCSFIT (NT, T (1), X (J1), NEW, METHOD, 1,TI,XI,XPI)

               NEW = .FALSE.
               IF (ABS (XPI) .LT. EPS) XPI = SIGN (EPS, XPI)

               IF (ABS (XI - XE) .GT. TOLER) THEN
                  IF (ITER .LT. MXITER) THEN
                     TI = MAX (MIN (TI - (XI - XE) / XPI, TMAX), TMIN)
                     GO TO 30

                  ELSE ! This should help uncover errors in the application:

                     WRITE (6, '(A,/,A,I2,1P,4E14.6)')
     >                  '0** PLSINTYZ: Bad input data likely.',
     >                  '0ITER,XE,XI,TI,XPI: ', ITER,XE,XI,TI,XPI
                     WRITE (6, '(1X,A,5F9.5)')
     >                  'XE,X(J1:J2): ', XE,X(J1),X(J1+1),X(J1+2),X(J2),
     >                  'TI,T(1:4):   ', TI,T,
     >                  'XE,Y(J1:J2): ', XE,Y(J1),Y(J1+1),Y(J1+2),Y(J2),
     >                  'XE,Z(J1:J2): ', XE,Z(J1),Z(J1+1),Z(J1+2),Z(J2)
                     WRITE (6, '(1X,A,5I9)')
     >                  'NDATA,I1,I2,J1,J2: ', NDATA,I1,I2,J1,J2
                     WRITE (6, '(1X, A1,/,(1X,8F9.5))') 'X', X
                     WRITE (6, '(1X, A1,/,(1X,8F9.5))') 'Y', Y
                     WRITE (6, '(1X, A1,/,(1X,8F9.5))') 'Z', Z

                  END IF
               END IF

C           Evaluate Y and Z at this T:

            NEW = .TRUE.

            CALL LCSFIT (NT, T (1), Y (J1), NEW, METHOD, 1, TI, YI, YI)
            YEVAL (IEVAL) = YI

            CALL LCSFIT (NT, T (1), Z (J1), NEW, METHOD, 1, TI, YI, YI)
            ZEVAL (IEVAL) = YI

         END IF

C        The next evaluation point may be in the same interval.
C        ------------------------------------------------------

         IF (IEVAL .LT. NEVAL) THEN      ! Skips this if NEVAL = 1.
            IEVAL = IEVAL + 1
            XE = XEVAL (IEVAL)
            IF (ENDPT) GO TO 10          ! LEFT is not valid in this case

            IF (XE * ARROW .GE. X (LEFT)     * ARROW  .AND.
     >          XE * ARROW .LT. X (LEFT + 1) * ARROW) THEN
               IF (XE .NE. X (I1) .AND. XE .NE. X (I2)) GO TO 20  ! Can't avoid
            END IF

            GO TO 10                     ! Else more work required.

         END IF

C     Termination.
C     ------------

      RETURN
      END
