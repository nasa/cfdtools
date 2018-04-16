C+----------------------------------------------------------------------
C
      SUBROUTINE COMPASS (XC, YC, R, NPTS, X, Y, XINT, YINT, INDEX, IER)
C
C     One-liner:  Locates the intersection of a circle and a curve
C
C     Description and usage:
C     ----------------------
C
C        COMPASS serves one function of a geometry compass: it locates
C     where a given circle meets a (somewhat) arbitrary curve defined by
C     discrete points in 2-space.  (See warning below.)  It was prompted
C     by a need to adjust the radial extent of an O-grid calculated by a
C     hyperbolic PDE method, which marches out from the body surface and
C     stops at an indefinite distance from the surface.
C
C        Spline-extrapolation of the line is possible, but will normally
C     be guarded against by the application.
C
C        Warning:  the points along the curve from 1 : NPTS are assumed
C     to be in the order of increasing distance from (XC, YC).  Spurious
C     results are possible if this is not true.
C
C     Arguments:
C     ----------
C
C     Name    Type/Dimension  I/O/S  Description
C     XC,     R               I      Coordinates of the center of the
C     YC                             circle.
C     R       R               I      Radius of the circle.
C     NPTS    I               I      Number of points defining the curve.
C     X,      R (NPTS)        I      Input curve coordinates.  The points
C     Y                              must be distinct but need not be
C                                    monotonic in X.  Their DISTANCES from
C                                    the circle center are assumed to be
C                                    monotonically increasing.
C     XINT,   R                 O    Coordinates of desired point of
C     YINT                           intersection.
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
C                                             failed: the circle does not
C                                             intersect the curve at all.
C
C     Method:
C     -------
C
C        To avoid storing more than a few arc-lengths, the relevant interval
C     is initially isolated by a backward sequential search (influenced by
C     the initial application, where the last point is likely to have been
C     inserted by LINEAR extrapolation for the very purpose of avoiding
C     nonlinear extrapolation).
C
C        A Newton iteration is then performed to solve the equation
C
C     F(T)  =  (X(T) - XC) ** 2  +  (Y(T) - YC) ** 2  -  R ** 2  =  0
C
C     using X = X(T), Y = Y(T) and their derivatives as given by the local
C     spline method of LCSFIT, which needs only the Ts of at most 4 data
C     points in the neighborhood.  (F'(T) involves X'(T) and Y'(T).)
C
C     Procedures:
C     -----------
C
C     CHORD     Safeguarded chord-length function
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
C     07/23/91  D.A.Saunders  COMPASS adapted from PLSINTRP in aid of
C                             the HYPEROH grid generator application.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   NPTS, INDEX, IER
      REAL
     >   XC, YC, R, X (NPTS), Y (NPTS), XINT, YINT

C     Local constants.

      INTEGER
     >   MXITER
      REAL
     >   ONE, TOL, TWO, ZERO
      LOGICAL
     >   NEW
      CHARACTER
     >   METHOD * 1

      PARAMETER
     >  (MXITER = 12,
     >   METHOD = 'B',       ! "Bessel" = "loose" fit should always suffice
     >   NEW    = .TRUE.,    ! Alternating between X and Y precludes efficiency
     >   ONE    = 1.E+0,
     >   TOL    = 1.E-6,     ! Tolerance for F ~ 0 by Newton iteration,     
     >   TWO    = 2.E+0,     !    suited to VAX single precision; division
     >   ZERO   = 0.0E+0)    !    by R**2 avoids the usual relative test

C     Local variables.

      INTEGER
     >   I, ITER, I1, I2, NT
      REAL
     >   DX, DY, FI, FPI, RR, T (4), TI, XI, XPI, YI, YPI

C     Procedures.

      REAL
     >   CHORD
      EXTERNAL
     >   CHORD, LCSFIT

C     Execution.
C     ----------

      IER = 0
      RR = ONE / R
      I = NPTS

C     Isolate the interval containing the point of intersection in a way
C     that avoids storing all the arc-lengths:

   10 CONTINUE
         DX = X (I) - XC
         DY = Y (I) - YC
         FI = (DX * RR) ** 2 + (DY * RR) ** 2 ! Avoid squaring large nos.
         IF (FI .GT. ONE) THEN
            I = I - 1
            IF (I .GT. 0) GO TO 10
         END IF

      INDEX = I

C     Check for an exact hit (which could affect the application):

      IF (FI .EQ. ONE) THEN
         INDEX = -I
         XINT = X (I)
         YINT = Y (I)
         GO TO 99
      END IF

C     To locate the desired point more exactly on an arbitrary curve,
C     we need to solve for the arc length T within the identified interval.
C     We take advantage of the local spline method employed by LCSFIT:
C     it needs at most 4 points to define the curve in this region.

      I1 = MAX (1, I - 1)
      I2 = MIN (NPTS, I + 2)
      NT = I2 - I1 + 1

      T (1) = ZERO
      T (2) = CHORD (X, Y, I1, I1 + 1)
      T (3) = CHORD (X, Y, I1 + 1, I1 + 2) + T (2)
      T (4) = CHORD (X, Y, I2 - 1, I2) + T (3)    ! Not used at ends

C     Newton iteration for the T that corresponds to the target distance R.

      TI = T (2)
      ITER = 0

   20 CONTINUE
         ITER = ITER + 1

         CALL LCSFIT (NT, T (1), X (I1), NEW, METHOD, 1, TI, XI, XPI)
         CALL LCSFIT (NT, T (1), Y (I1), NEW, METHOD, 1, TI, YI, YPI)

         DX = XI - XC
         DY = YI - YC
         FI = (DX * RR) ** 2 + (DY * RR) ** 2 - ONE   ! Avoid large nos.

         IF (ABS (FI) .GT. TOL) THEN
            IF (ITER .LT. MXITER) THEN
               FPI = ((DX * RR) * XPI + (DY * RR) * YPI) * (TWO * RR)
C****          TI = MAX (MIN (TI - FI / FPI, T (NT)), ZERO)
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

   99 RETURN
      END
