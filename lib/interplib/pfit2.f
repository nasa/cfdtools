C+----------------------------------------------------------------------
C
      SUBROUTINE PFIT2 (NGEOM, XGEOM, YGEOM, TGEOM, PSCOFS, CLOSED, IER)
C
C  ACRONYM:   Parametric spline FIT in 2 dimensions
C             -                 ---    -
C  PURPOSE:   PFIT2 fits a 2-D interpolating parametric cubic spline to
C             an open or closed curve represented by the given  set  of
C             discrete points.  The parametric variable T is cumulative
C             chord length - the usual approximation to arc length.
C
C             See PEVAL2 for evaluation of the spline at given value(s)
C             of T.
C
C  METHOD:    TGEOM(*) is set up here with cumulative chord length, and
C             returned for possible use by the calling program.
C
C             The fitted splines for  X vs. T and  Z vs. T are standard
C             interpolating cubic splines with natural  end  conditions
C             (Y" = 0), except when the curve is to be closed smoothly,
C             in which case boundary conditions are periodic.
C             
C  ARGUMENTS:
C   ARG     TYPE I/O/S  DIM      DESCRIPTION
C  NGEOM      I    I     -       Number of data points to be fitted.
C  XGEOM,     R    I   NGEOM     Data point coordinates.
C   YGEOM
C  TGEOM      R   S/O  NGEOM     TGEOM(1)=0.; TGEOM(J) = TGEOM(J-1) +
C                                distance between points J-1 and J.
C  PSCOFS     R   S/O  NGEOM*6   PSCOFS(*,L,M) are the spline coefs. of
C                                type L (1,2,3 for B,C,D) for dimension
C                                M (1,2 for X,Y) vs. T.
C  CLOSED     L    I     -       .TRUE. means the data points represent
C                                a closed curve. In this case, the call-
C                                ing program MUST duplicate the first
C                                data point (XGEOM(1), YGEOM(1)) as
C                                (XGEOM(NGEOM), YGEOM(NGEOM)), and
C                                smoothness in the spline will be sought
C                                here as at all the other data points;
C                                .FALSE. means the end-points are treat-
C                                ed as distinct.
C  IER          I    O   -       IER = 0: No errors detected;
C                                IER = 1: Too few data points:
C                                         (NGEOM < 2 for open curve or
C                                          NGEOM < 3 for closed curve);
C                                IER = 4: Closed curve case only:
C                                         First and last data points do
C                                         not match.
C  PROCEDURES:
C     CSFIT:  Fits conventional interpolating cubic spline
C
C  ENVIRONMENT:  FORTRAN 77 + minor extensions
C
C  HISTORY:
C     03/29/94  D.A.Saunders  Adapted from PFIT3 for the airfoil case.
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER    NGEOM, IER
      REAL       XGEOM (NGEOM), YGEOM (NGEOM), TGEOM (NGEOM),
     >           PSCOFS (NGEOM, 3, 2)
      LOGICAL    CLOSED

C     Local constants:

      REAL       ZERO
      PARAMETER (ZERO = 0.E+0)

C     Local variables:

      INTEGER    IENDL, J

C     Procedures:

      EXTERNAL   CSFIT

C     Execution:

C     Too few data points?

      IF (NGEOM .LT. 2 .OR. (CLOSED .AND. NGEOM .LT. 3)) THEN
         IER = 1
         GO TO 999
      END IF

C     Fit natural splines or splines with periodic end conditions:

      IF (CLOSED) THEN
         IENDL = 4
      ELSE
         IENDL = 0
      END IF

C     Set up parametric variable values:

      TGEOM (1) = ZERO
      DO J = 2, NGEOM
         TGEOM (J) = TGEOM (J - 1) +
     >      SQRT ((XGEOM (J) - XGEOM (J - 1)) ** 2 + 
     >            (YGEOM (J) - YGEOM (J - 1)) ** 2)
      END DO

C     Fit spline for XGEOM vs. TGEOM:

      CALL CSFIT (NGEOM, TGEOM, XGEOM, IENDL, ZERO, 0, ZERO,
     >   PSCOFS (1, 1, 1), PSCOFS (1, 2, 1), PSCOFS (1, 3, 1), IER)
      IF (IER .NE. 0) GO TO 999

C     ... and for YGEOM vs. TGEOM:

      CALL CSFIT (NGEOM, TGEOM, YGEOM, IENDL, ZERO, 0, ZERO,
     >   PSCOFS (1, 1, 2), PSCOFS (1, 2, 2), PSCOFS (1, 3, 2), IER)

 999  RETURN
      END
