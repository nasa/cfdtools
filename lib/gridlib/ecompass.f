C+----------------------------------------------------------------------
C
      SUBROUTINE ECOMPASS (XC, YC, A, B, NPTS, X, Y, T,
     >                     XINT, YINT, TINT, INDEX, IER)
C
C     One-liner:  Locates the intersection of an ellipse and a curve
C
C     Description and usage:
C     ----------------------
C
C        ECOMPASS is the "ellipse" version of the circle/curve module
C     named COMPASS.  It locates where a given semi-ellipse meets a
C     (somewhat) arbitrary curve defined by discrete points in 2-space.
C     (See the warning below about the curve.)  It was prompted by a
C     need to adjust the radial extent of an O-grid calculated by a
C     hyperbolic PDE method, which marches out from the body surface
C     and stops at an indefinite distance from the surface.
C
C        Spline-extrapolation of the line is possible, but will normally
C     be guarded against by the application.
C
C        Warning:  the points along the curve are assumed to have
C     distances from (XC, YC) varying monotonically.  If this is not
C     true, an apparently valid solution may not be unique.  On the
C     other hand, non-convergence of the Newton iteration should be
C     a sure sign that there is NO solution.
C
C     Arguments:
C     ----------
C
C     Name    Type/Dimension  I/O/S  Description
C     XC,     R               I      Coordinates of the center of the
C     YC                             ellipse.
C     A,      R               I      X and Y axis lengths of the ellipse.
C     B
C     NPTS    I               I      Number of points defining the curve.
C     X,      R (NPTS)        I      Input curve coordinates.  The points
C     Y                              must be distinct but need not be
C                                    monotonic in X.  Their DISTANCES from
C                                    the ellipse center are assumed to be
C                                    monotonic.
C     T       R (NPTS)        I      Cumulative chord-lengths corresponding
C                                    to X and Y, increasing monotonically.
C     XINT,   R                 O    Coordinates of desired point of
C     YINT                           intersection.
C     TINT    R                 O    Corresponding chord-length value.
C     INDEX   I                 O    (XINT, YINT) was determined to lie
C                                    in the X,Y interval defined by
C                                    INDEX and INDEX + 1.  Special cases:
C                                    (1) INDEX < 0 means an exact hit:
C                                        the point matches point I = -INDEX.
C                                        The application may need to know
C                                        this if it would otherwise insert
C                                        (XINT, YINT) into X and Y.
C                                    (2) INDEX = 0 means the point is off
C                                        the end to the "left."
C                                    (3) INDEX = NPTS means it is off
C                                        the end to the "right."
C     IER     I                 O    Error return code:
C                                    IER = 0  means INDEX is valid;
C                                        = -1 means the Newton iteration
C                                             failed: the ellipse does not
C                                             intersect the curve at all.
C
C     Method:
C     -------
C
C        Since the initial application is working with chord-lengths to
C     redistribute points along the curve, no attempt is made to avoid
C     storing them: they are expected as input.  (Use module CHORD.)
C
C        This allows representation of the curve as X = X(T) and Y = Y(T).
C     T(NPTS) is taken as the initial guess for the Newton iteration used
C     to solve the equation
C
C        F(T) = ((X(T) - XC) / A) ** 2 + ((Y(T) - YC) / B) ** 2 - 1 = 0
C
C     via the local spline method of LCSFIT, which supplies X'(T) and
C     Y'(T) as needed for F'(T).
C
C     Procedures:
C     -----------
C
C     INTERVAL  Search utility used to identify INDEX
C     LCSFIT    Local cubic spline method for monotonic data
C
C     Environment:
C     ------------
C
C     VAX/VMS; FORTRAN 77
C     IMPLICIT NONE, 8-character symbolic names, and "!" as comment
C     character are not (yet) standard.
C
C     Author: David Saunders, NASA Ames/Sterling Software, Palo Alto, CA.
C     -------
C
C     History:
C     --------
C
C     07/23/91  D.A.Saunders  ECOMPASS adapted from COMPASS when the
C                             HYPEROH grid generator application was
C                             modified to generate elliptical outer
C                             boundaries.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   NPTS, INDEX, IER
      REAL
     >   XC, YC, A, B, X (NPTS), Y (NPTS), T (NPTS), XINT, YINT, TINT

C     Local constants.

      INTEGER
     >   MXITER
      REAL
     >   ONE, TWO, TOL, ZERO
      LOGICAL
     >   NEW
      CHARACTER
     >   METHOD * 1

      PARAMETER
     >  (MXITER = 12,
     >   METHOD = 'B',       ! "Bessel" = "loose" fit should always suffice
     >   NEW    = .TRUE.,    ! Alternating between X and Y precludes efficiency
     >   TOL    = 1.E-6,     ! Tolerance for F ~ 0 by Newton iteration
     >   ONE    = 1.E+0,     !    (suited to VAX single precision;
     >   TWO    = 2.E+0,     !    shouldn't need to be relative)
     >   ZERO   = 0.0E+0)

C     Local variables.

      INTEGER
     >   ITER
      REAL
     >   DX, DY, FI, FPI, RA, RB, TI, TOLER, XI, XPI, YI, YPI

C     Procedures.

      EXTERNAL
     >   INTERVAL, LCSFIT

C     Execution.
C     ----------

      IER = 0
      RA = ONE / A
      RB = ONE / B

C     Newton iteration for the T that corresponds to (X, Y) on the ellipse:

      TI = T (NPTS)
      ITER = 0

   20 CONTINUE
         ITER = ITER + 1

         CALL LCSFIT (NPTS, T, X, NEW, METHOD, 1, TI, XI, XPI)
         CALL LCSFIT (NPTS, T, Y, NEW, METHOD, 1, TI, YI, YPI)

         DX = XI - XC
         DY = YI - YC
         FI = (DX * RA) ** 2 + (DY * RB) **2 - ONE

         IF (ABS (FI) .GT. TOL) THEN
            IF (ITER .LT. MXITER) THEN
               FPI = ((DX * RA) * RA * XPI +
     >                (DY * RB) * RB * YPI) * TWO
C****          TI = MAX (MIN (TI - FI / FPI, T (NPTS)), ZERO)
               TI = TI - FI / FPI       ! Allow extrapolation.  Even then,
                                        ! unreasonable cases can fail...
               GO TO 20

            ELSE
               IER = -1
               GO TO 99
            END IF

         END IF

      XINT = XI
      YINT = YI
      TINT = TI

C     Identify the interval containing (XINT, YINT) in case a point
C     is to be inserted in X and Y.

      CALL INTERVAL (NPTS, T, TI, ONE, INDEX)

C     Check for an exact hit (which could affect the application):

      IF (TI .GT. T (NPTS)) THEN
         INDEX = NPTS
      ELSE IF (TI .LT. ZERO) THEN
         INDEX = 0
      ELSE IF (TI .EQ. T (INDEX)) THEN
         INDEX = -INDEX
      ELSE IF (TI .EQ. T (NPTS)) THEN
         INDEX = -NPTS           ! INTERVAL gives 1 <= INDEX <= NPTS - 1
      END IF

   99 RETURN
      END
