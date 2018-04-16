C+---------------------------------------------------------------------
C
      SUBROUTINE PSTVAL (NT, T, X, Y, XP, YP, XPP, YPP, NGEOM, XGEOM,
     >                   YGEOM)
C
C  ACRONYM:    Parametric Spline eVALuation for given T(s)
C              -          -      ----                 -
C  PURPOSE:    PSTVAL is an alternative to PSEVAL for the case when the
C              curve fitted by PSFIT is to be evaluated (X and Y) for a
C              given value or values of parametric variable T.  This is
C              in contrast to evaluation for a given X.
C
C              Evaluation for given T is much more straightforward than
C              for given X,  and may be appropriate for graphical work.
C
C              Provision is made for convenient evaluations at the data
C              points (no need to pass in the corresponding Ts).
C
C              The optional derivatives with respect to T will probably
C              be used for curvature calculations.  See also CURV2D and
C              XDERIVS.
C
C  METHOD:     The splines XGEOM vs. TGEOM and YGEOM vs. TGEOM, already
C              calculated by PSFIT, are evaluated at each of the values
C              of T supplied, or at TGEOM (1 : -NT) if NT < 0 is input.
C              First and second derivatives with respect to  T  (not X)
C              are also calculated.  As with PSEVAL, these are computed
C              in the order x", x', x and y", y', y.   If x" and y" are
C              not needed, pass the same arrays as for x' and y'. If x'
C              and y' are not needed either, pass x and y appropriately
C              to save storage.
C
C  ARGUMENTS:
C   ARG    TYPE  I/O/S   DIM   DESCRIPTION
C   NT      I      I      -    Number of evaluations required.  NT < 0
C                              is taken to mean the evaluations are to
C                              be at the data points 1:-NT.  |NT|>=1.
C   T       R      I      NT   Arbitrary values of parametric variable
C                              T for which evaluations are desired.  If
C                              NT < 0, the existing values of T for the
C                              data points (TGEOM(*)) stored internally
C                              by PSFIT are used - no need to rederive
C                              them in the application.
C   X       R      O      NT   Corresponding values of X on the curve.
C                              May be same locations as T(*).
C   Y       R      O      NT   Corresponding values of Y on the curve.
C   XP      R      O      NT   First derivatives, DX/DT.
C   YP      R      O      NT   First derivatives, DY/DT.
C   XPP     R      O      NT   Second derivatives of X and Y w.r.t. T.
C   YPP     R      O      NT   (See above if these are not needed.)
C   NGEOM   R      I      -    Number of data points fitted by PSFIT.
C   XGEOM   R      I    NGEOM  Geometry abscissas fitted by PSFIT.
C   YGEOM   R      I    NGEOM  Geometry ordinates fitted by PSFIT.
C
C  INTERNAL COMMON BLOCK:
C   /PSWORK/ Work arrays, set up by PSFIT.
C            Dimensions are set to MAXPS by PARAMETER statement.
C   VAR    TYPE I/O/S DIM  DESCRIPTION
C   TGEOM   R     I MAXPS  Values of the paramateric variable at the
C                          data points, set up by PSFIT.
C   PSCOFS  R     I MAXPS  Cubic spline coefficients from CSFIT.
C                   *3*2   An array of size MAXPS*3 is required for
C                          XGEOM and another for YGEOM.
C   XLOC    R     S   -    Copy of X, abscissa at which PSEVAL is seeking
C                          an ordinate.
C   NCURV   I     O   -    Copy of NGEOM.
C
C  ERROR HANDLING: None.
C
C  EXTERNAL REFERENCES:
C   CSDVAL  Evaluates a 1-dimensional cubic spline along with its first
C           and second derivatives.
C
C  ENVIRONMENT:  VAX/VMS, FORTRAN 77
C
C  HISTORY:
C  Sept. 85  DAS  Adapted from PSEVAL for the simpler situation,
C                 in case it proves handy some time.
C  April 86  DAS  Allowed for caller to overwrite T(*) with X(*).
C  10/16/91  DAS  Derivatives returned are with respect to T, not X,
C                 now (to facilitate curvature calculations).  Gave
C                 up on allowing T to be overwritten with X, to
C                 avoid 2*NT calls to CSDVAL where two calls will do.
C  10/25/91  DAS  Derivatives at the data points required the internal
C                 TGEOM values.  Use of TSUBJ for all of them would
C                 be ugly.  Thus NT < 0 now avoids having to input Ts.
C  04/01/95  DAS  MAXPS had not been raised to 1001 to match PSFIT.
C
C  AUTHOR: David Saunders, Sterling Software/NASA Ames, Moffett Field, CA
C
C----------------------------------------------------------------------

      IMPLICIT NONE

C  *  Arguments:

      INTEGER
     >   NT, NGEOM
      REAL
     >   T (*), X (*), XP (*), XPP (*), Y (*), YP (*), YPP (*),
     >   XGEOM (NGEOM), YGEOM (NGEOM)

C  *  Internal COMMON:

      INTEGER
     >   MAXPS
      PARAMETER
     >  (MAXPS = 1001)
      INTEGER
     >   NCURV
      REAL
     >   PSCURV, TGEOM, PSCOFS, XLOC
      COMMON /PSWORK/
     >   PSCURV (MAXPS), TGEOM (MAXPS), PSCOFS (MAXPS,3,2), XLOC, NCURV

C  *  Local variables:

      INTEGER
     >   NPTS

C  *  Execution:

      NPTS = ABS (NT)

      IF (NT .GT. 0) THEN

C  *     Evaluate the XGEOM vs. TGEOM spline at each input T:

         CALL CSDVAL (NGEOM, TGEOM, XGEOM, NPTS, T,
     >                PSCOFS (1,1,1), PSCOFS (1,2,1), PSCOFS (1,3,1),
     >                X, XP, XPP)

C  *     ... and the YGEOM vs. TGEOM spline at each input T:

         CALL CSDVAL (NGEOM, TGEOM, YGEOM, NPTS, T,
     >                PSCOFS (1,1,2), PSCOFS (1,2,2), PSCOFS (1,3,2),
     >                Y, YP, YPP)

      ELSE  ! NT < 0 means evaluate at the original data points:

C        First, the X derivatives:

         CALL CSDVAL (NGEOM, TGEOM, XGEOM, NPTS, TGEOM,
     >                PSCOFS (1,1,1), PSCOFS (1,2,1), PSCOFS (1,3,1),
     >                X, XP, XPP)

C  *     ... and now the Y derivatives:

         CALL CSDVAL (NGEOM, TGEOM, YGEOM, NPTS, TGEOM,
     >                PSCOFS (1,1,2), PSCOFS (1,2,2), PSCOFS (1,3,2),
     >                Y, YP, YPP)

      END IF

      RETURN
      END
