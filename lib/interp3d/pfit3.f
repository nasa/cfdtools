C+---------------------------------------------------------------------
C
      SUBROUTINE PFIT3 ( NGEOM, XGEOM, YGEOM, ZGEOM, TGEOM, PSCOFS,   
     +                   CLOSED, IER )
C
C  ACRONYM:   Parametric spline FIT in 3 dimensions
C             -                 ---    -
C  PURPOSE:   PFIT3 fits a 3-D interpolating parametric cubic spline to
C             an open or closed curve represented by the given  set  of
C             discrete points.  The parametric variable T is cumulative
C             chord length - the usual approximation to arc length.
C
C             Evaluation of the spline for a given value or values of T  
C             is provided by subroutine PEVAL3.
C
C  METHOD:    TGEOM(*) is set up here with cumulative chord length, and
C             returned for possible use by the calling program.
C
C             The fitted splines for  X vs. T, Y vs. T, and Z vs. T are
C             standard interpolating cubic splines with the natural end
C             conditions (Y"=0),  except when the curve is to be closed
C             smoothly, in which case boundary conditions are periodic.
C             
C  ARGUMENTS:
C   ARG     TYPE I/O/S  DIM      DESCRIPTION
C  NGEOM      I    I     -       Number of data points to be fitted.
C  XGEOM,     R    I   NGEOM     Coordinates of the data points to be
C   YGEOM,                       fitted with a parametric spline.
C   ZGEOM
C  TGEOM      R   S/O  NGEOM     TGEOM(1)=0.; TGEOM(J) = TGEOM(J-1) +
C                                distance between points J-1 and J.
C  PSCOFS     R   S/O  NGEOM,3,3 PSCOFS(*,L,M) are the spline coefs. of
C                                type L (1,2,3 for B,C,D) for dimension
C                                M vs T (1,2,3 for X,Y,Z).  Note that
C                                the calling program just needs work
C                                space of length NGEOM*9 for PSCOFS(*).
C  CLOSED     L    I     -       .TRUE. means the data points represent
C                                a closed curve. In this case, the call-
C                                ing program MUST duplicate the first
C                                data point (XGEOM(1),YGEOM(1),ZGEOM(1))
C                                as (XGEOM(NGEOM),YGEOM(NGEOM),ZGEOM(NGEOM)),
C                                and smoothness in the spline will be sought
C                                here as at all the other data points;
C                                .FALSE. means the end-points are treat-
C                                ed as distinct, and are not joined
C                                (though in fact they MAY be the same
C                                point, as in the case of an airfoil
C                                with a sharp trailing edge).
C  IER          I    O     -     Error return code:
C                                IER = 0: No errors were detected;
C                                IER = 1: Too few data points:
C                                         (NGEOM < 2 for open curve or
C                                          NGEOM < 3 for closed curve);
C                                IER = 4: Closed curve case only:
C                                         First and last data points do
C                                         not match.
C  ERROR HANDLING:
C   See IER description.  It is the user's responsibility to check IER
C   on return from PFIT3, and possibly stop before calling PEVAL3.
C
C  EXTERNAL REFERENCES:
C   CSFIT:  Fits interpolating cubic spline (many end condition choices)
C
C  ENVIRONMENT:  VAX/VMS FORTRAN
C
C  AUTHOR:  David Saunders, Sterling Software, Palo Alto, CA.
C
C  DEVELOPMENT HISTORY:
C       DATE   INITIALS  DESCRIPTION
C     09/18/86   DAS     Adapted from PSFIT, which deals with just X, Y.
C     05/11/87   DAS     Bug: Third fit had Y where it should have had Z.
C
C----------------------------------------------------------------------

      IMPLICIT  NONE

C  *  Arguments:

      INTEGER   NGEOM, IER
      REAL      XGEOM(NGEOM), YGEOM(NGEOM), ZGEOM(NGEOM), TGEOM(NGEOM),
     +          PSCOFS(NGEOM,3,3)
      LOGICAL   CLOSED

C  *  Local constants:

      REAL      ZERO
      PARAMETER ( ZERO = 0.E+0 )

C  *  Local variables:

      INTEGER   IENDL, J

C  *  Procedures:

      EXTERNAL  CSFIT

C  *  Execution:

C  *  Too few data points?

      IF ( NGEOM.LT.2 .OR. (CLOSED.AND.NGEOM.LT.3) ) THEN
         IER = 1
         GO TO 999
      END IF

C  *  Fit natural splines or splines with periodic end conditions:

      IF ( CLOSED ) THEN
            IENDL = 4
         ELSE
            IENDL = 0
      END IF

C  *  Set up parametric variable values:

      TGEOM(1) = ZERO
      DO 200, J = 2, NGEOM
         TGEOM(J) = TGEOM(J-1) + SQRT ( ( XGEOM(J) - XGEOM(J-1) ) ** 2 + 
     +                                  ( YGEOM(J) - YGEOM(J-1) ) ** 2 +
     +                                  ( ZGEOM(J) - ZGEOM(J-1) ) ** 2 )
 200  CONTINUE

C  *  Compute spline for XGEOM vs. TGEOM:

      CALL CSFIT ( NGEOM, TGEOM, XGEOM, IENDL, ZERO, 0, ZERO,
     +             PSCOFS(1,1,1), PSCOFS(1,2,1), PSCOFS(1,3,1), IER )
      IF ( IER.NE.0 ) GO TO 999

C  *  Compute spline for YGEOM vs. TGEOM:

      CALL CSFIT ( NGEOM, TGEOM, YGEOM, IENDL, ZERO, 0, ZERO,
     +             PSCOFS(1,1,2), PSCOFS(1,2,2), PSCOFS(1,3,2), IER )
      IF ( IER.NE.0 ) GO TO 999

C  *  Compute spline for ZGEOM vs. TGEOM:

      CALL CSFIT ( NGEOM, TGEOM, ZGEOM, IENDL, ZERO, 0, ZERO,
     +             PSCOFS(1,1,3), PSCOFS(1,2,3), PSCOFS(1,3,3), IER )
C
 999  RETURN
      END
