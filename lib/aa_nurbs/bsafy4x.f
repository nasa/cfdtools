C+------------------------------------------------------------------------------
C
      SUBROUTINE BSAFY4X (C, TLE, SURFACE, NX, X, Y, LUNERR, IER)
C
C ONE-LINER: B-Spline AirFoil evaluation of Y for given X(s) (one surface)
C            - -      -  -                  - -         -
C PURPOSE:
C
C        BSAFY4X discretizes ONE surface of an airfoil, given in DT_NURBS
C     B-spline curve form, at specified X abscissa(s) as opposed to given
C     (or generated) values of parametric variable T.  See also BSAF2XY.
C
C METHOD:
C
C        The B-spline curve should represent the airfoil in clockwise
C     wrap-around form, from lower trailing edge to upper trailing edge.
C     Unless TLE is input (not 999.), general purpose utility BSZEROMN
C     is used to locate it as the value of the parametric variable T
C     corresponding to the leading edge.  Then the interval [T1, TLE]
C     or [TLE, TN], depending on SURFACE, is used with BSEVAL to estimate
C     the Y for each X, where T1 and TN are the first and last knots.
C
C        The leading edge point presents a problem: estimating its Y from
C     its X via BSEVAL could well fail because of round-off difficulties.
C     Therefore, if the smallest X is arbitrarily close to the leading edge
C     as determined from TLE, then the corresponding Y is determined from
C     TLE, not this X.
C
C ENVIRONMENT:
C     VAX/VMS; FORTRAN 77, with...
C     > IMPLICIT NONE
C     > Trailing ! comments
C     > Names up to 8 characters
C
C HISTORY:
C     04/21/93  D.A.Saunders  Adaptation of BSAF2XY.
C     02/01/94       "        Dispensed with DOUBLE PRECISION.  Compile with
C                             /REAL_LENGTH=64 or -r8 switches as necessary.
C
C AUTHOR: David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments
C     ---------

      REAL    C (*)     ! (I)   DT_NURBS B-spline curve vector:
                        !       C(1) = 1 (# parametric variables);
                        !       C(2) = 2 (for X and Y);
                        !       C(3) = polynomial order k (4 for cubics);
                        !       C(4) = # control points, N;
                        !       C(5) = pointer for efficient evaluations;
                        !       C(6 : 5 + N + k) = knot vector,
                        !       followed by X coords. of control pts.,
                        !       followed by Y coords. of control pts.

      REAL    TLE       ! (I/O) The parametric variable value corresponding
                        !       to the leading edge.  Enter 999. if it is
                        !       unknown, in which case it will be calculated
                        !       and returned.

      CHARACTER SURFACE * (*)   ! (I)   U[pper] or L[ower]: only character 1
                                !       is checked; it must be UPPER case.

      INTEGER NX        ! (I)   Number of Xs at which to evaluate Y.  NX >= 1.

      REAL    X (NX)    ! (I)   The target abscissas, which may be
                        !       either ascending or descending, and
                        !       may or may not span a whole surface.

      REAL    Y (NX)    ! (O)   The ordinates estimated for the indicated
                        !       surface at each X, in the same order.

      INTEGER LUNERR    ! (I)   Logical unit for error messages.

      INTEGER IER       ! (O)   Success/error code:
                        !        0 means no problem was encountered;
                        !       12 means <to be completed>

C     Procedures
C     ----------

      EXTERNAL BSEVAL   ! Evaluates a 2-D B-spline curve (Y for X here).
      EXTERNAL BSZEROMN ! Finds a zero, min, or max of a B-spline curve.
      EXTERNAL DTSPVL   ! Evaluates a B-spline curve at the given T.

C-------------------------------------------------------------------------------

C     Local constants
C     ---------------

      INTEGER
     >   MAXK, NWORK
      REAL
     >   EPS, FLAG, ZERO
      PARAMETER
     >  (MAXK    = 7,              ! Local storage handles polynomial deg. <= 6
     >   NWORK   = 5 * MAXK - 2,   ! This is why the degree is limited
     >   EPS     = 1.E-6,          ! Tolerance for XMIN vs. XLE relative to XAVG
     >   FLAG    = 999.E+0,        ! Tells whether a good TLE is input or not
     >   ZERO    = 0.E+0)

C     Local variables
C     ---------------

      INTEGER
     >   I1, I2, IMIN
      REAL
     >   DUM, T1, T2, TA, TB, WORK (NWORK), XY (2)

C     Execution
C     ---------

      I2 = 5 + INT (C (3) + C (4))
      T2 = C (I2)                      ! Last knot
      T1 = C (6)                       ! First knot

      IF (TLE .EQ. FLAG) THEN

C        Find the leading edge value of T by minimizing X w.r.t. T:

         CALL BSZEROMN ('MIN', 1, C, T1, T2, TLE, XY, LUNERR, IER)
         IF (IER .NE. 0) GO TO 800

      ELSE

C        We still need X and Y at the leading edge as functions of TLE:

         CALL DTSPVL (TLE, C, WORK, NWORK, XY, IER)
         IF (IER .NE. 0) GO TO 810

      END IF

C     Avoiding the leading edge X among the main evaluations makes
C     for some awkward gyrations:

      IF (X (1) .LT. X (NX)) THEN
         IMIN = 1
      ELSE
         IMIN = NX
      END IF

      IF (ABS (X (IMIN) - XY (1)) .LT. EPS * (X (1) + X (NX))) THEN
         Y (IMIN) = XY (2)      ! Leading edge Y
         IF (IMIN .EQ. 1) THEN
            I1 = 2
            I2 = NX
         ELSE
            I1 = 1
            I2 = NX - 1
         END IF
      ELSE
         I1 = 1
         I2 = NX
      END IF

C     Distinguish the surfaces:

      IF (SURFACE (1:1) .EQ. 'L') THEN
         TA = T1
         TB = TLE
      ELSE
         TA = TLE
         TB = T2
      END IF


C     Evaluate Y at each X:

      CALL BSEVAL ('Y', C, .TRUE., I2 - I1 + 1, X (I1), TA, TB, DUM,
     >             Y (I1), LUNERR, IER)
      IF (IER. NE. 0) GO TO 820

      GO TO 999


C     Error handling.

  800 WRITE (LUNERR, 1010) 'BSZEROMN', IER
      GO TO 999
      
  810 WRITE (LUNERR, 1010) 'DTSPVL', IER
      GO TO 999

  820 WRITE (LUNERR, 1010) 'BSEVAL', IER


  999 RETURN

C     Formats.

 1010 FORMAT ('0BSAFY4X: Bad return from ', A, '.  IER: ', I4)

      END
