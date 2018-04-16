C+---------------------------------------------------------------------
C
      SUBROUTINE PSEVAL (NX, X, Y, YP, YPP, T1, T2, NGEOM, XGEOM,
     >                   YGEOM, TOL, IER)
C
C  ACRONYM:    Parametric Spline EVALuation (Y for given X, not T)
C              -          -      ----
C  PURPOSE:
C                 PSEVAL will normally be called after subroutine PSFIT
C              has fit a parametric cubic spline to a given  curve,  to
C              evaluate the monotonic sub-curve defined by  T1 and  T2,
C              at arbitrary abscissas which are assumed to be  in-range
C              for this sub-curve.  PSEVAL computes ordinates and first
C              and second derivatives.  (See NOTES below if derivatives
C              are not needed, and for tips on choosing T1, T2.)
C
C                 Use PSTVAL if evaluations for given T, not X, are all
C              that is needed.
C
C  METHOD:
C                 A zero finder is used to approximate the value of the
C              parametric variable corresponding to the given X  on the
C              subcurve defined by T1 and T2.   The corresponding Y  is
C              theneasily evaluated, as are the derivatives.   Outline:
C
C              FOR each abscissa X=X(J), J=1:NX, DO
C               * Use a zero finder to find the zero of the function
C                 represented by routine GETT in the range  [T1,T2].
C                 This value, TBEST,  is such that evaluation of the
C                 spline from PSFIT on  TGEOM, XGEOM  at TBEST gives
C                 an XTEST arbitrarily close to the given X.
C               * Evaluate the spline from PSFIT on TGEOM, YGEOM for
C                 the corresponding ordinate, Y(J).
C               * Compute 1st and 2nd derivatives by the chain rule.
C
C NOTES:
C         (1)  PSEVAL calculates derivatives whether they are needed or
C              not, for reasons of simplicity.  Undesirable storage may
C              be avoided, however, since the output values are assign-
C              ed in the order YPP(*), YP(*), Y(*).  For example, if no
C              derivatives are required, use the same array for each of
C              these arguments; if just ordinates and first derivatives
C              are desired  (but no second derivatives),  pass the same
C              array for YPP as for YP.
C
C         (2)  Derivatives with respect to  X  require derivatives with
C              respect to the parametric variable,  so CSDVAL is needed
C              rather than CSEVAL.  GETT does not need derivatives, but
C              it was considered cleaner to use CSDVAL in both places.
C
C         (3)  FUNCTION TSUBJ may be used to determine T1 and T2. TSUBJ
C              is included in the same source module as  PSEVAL & GETT.
C              Note that the end-points of sub-curves are not necessar-
C              ily at original data points.   (Consider, say,  a modest
C              number of points roughly in the form of a  circle.   The
C              fitted curve may "bulge" at the sides...)
C
C  ARGUMENTS:
C   ARG    TYPE  I/O/S   DIM   DESCRIPTION
C   NX      I      I      -    Number of evaluations required.  NX>=1.
C   X       R      I      NX   Arbitrary abscissas (except all must be
C                              valid for the indicated sub-curve).
C   Y       R      O      NX   Ordinates corresponding to X(*).
C   YP      R      O      NX   First derivatives, DY/DX (see NOTES).
C   YPP     R      O      NX   Second derivatives w.r.t. X (see NOTES).
C   T1,     R      I      -    Values of the parametric variable that
C    T2     R      I      -    delimit the desired monotonic sub-curve.
C                              T1 < T2 is required.  See NOTES.  Using
C                              TSUBJ ( J ) to define the value of T at
C                              data point J will commonly suffice, but
C                              the X-extent of sub-curves is not nece-
C                              ssarily bounded by the original points.
C   NGEOM   R      I      -    Number of data points fitted by PSFIT.
C   XGEOM   R      I    NGEOM  Geometry abscissas fitted by PSFIT.
C   YGEOM   R      I    NGEOM  Geometry ordinates fitted by PSFIT.
C   TOL     R      I      -    Convergence criterion for zero finder, which
C                              returns TBEST as its best approximation to
C                              the zero of function GETT(T).  TOL is applied
C                              relative to the data range as given by the
C                              total arc length.  For normal applications
C                              and 32-bit arithmetic, TOL=1.E-6 is suggested;
C                              for 64-bit arithmetic, TOL~1.E-14.
C   IER     I      O      -    Error return code from ZEROIN:
C                                0 if no error;
C                                1 if F(T1), F(T2) are not of opposite sign;
C                                2 if T1 >= T2.
C
C  ERROR HANDLING: Nonzero IER from ZEROIN will force early return.
C
C  PROCEDURES:
C   CSDVAL  Evaluates a 1-dimensional cubic spline along with its first
C           and second derivatives.  (Also used by GETT.)
C   GETT    Real function passed as argument to ZEROIN.  GETT is in
C           the same source module as PSEVAL, and should not concern
C           the user unless the following COMMON block is changed.
C   ZEROIN  Finds zero of a function which changes sign in a given interval.
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
C  ENVIRONMENT:  VAX/VMS, FORTRAN 77
C
C  HISTORY:
C  Jan.  83  R.Lefkowitz  Original design and coding.
C  Sept. 84    "    "     CSDVAL used in place of IMSL's ICSEVU, DCSEVU.
C  Sept. 85  D.Saunders   Description of T1, T2 changed to reflect move
C                         to cumulative chord length from point index as
C                         the parametric variable. No code changes needed.
C  April 86    "    "     Calling sequence to ZEROIN changed (no LUN now).
C  Oct.  87    "    "     EPS=1.E-12 was too extreme - try 1.E-8. (What IS
C                         the proper fix here for zero DX/DT?)
C  Oct.  91    "    "     TOL is applied relative to total arc-length now.
C                         So is EPS (back to 1.E-10, with EPS**3 avoided).
C                         Other cosmetics.
C
C  AUTHOR: Rosalie Lefkowitz, Sterling Software/NASA Ames, Moffett Field, CA
C
C----------------------------------------------------------------------

      IMPLICIT NONE

C  *  Arguments:

      INTEGER
     >   NX, NGEOM, IER
      REAL
     >   T1, T2, TOL, X (NX), Y (NX), YP (NX), YPP (NX), XGEOM (NGEOM),
     >   YGEOM (NGEOM)

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

C  *  Local constants:

      REAL
     >   EPS
      PARAMETER
     >  (EPS = 1.E-10)  ! EPS is a limit on |DX/DT| to guard against
                        ! division by zero.

C  *  Local variables:

      INTEGER
     >   J
      REAL
     >   EPSIL, TBEST, TOLER, XD1, XD2, YD1, YD2, YJ

C  *  Procedures:

      REAL
     >   GETT, ZEROIN
      EXTERNAL
     >   CSDVAL, GETT, ZEROIN

C  *  Execution:

      TOLER = TOL * TGEOM (NGEOM)
      EPSIL = EPS * TGEOM (NGEOM)

      DO 100, J = 1, NX

C  *     Estimate parametric variable value TBEST corresponding to X (J),
C        using spline for TGEOM, XGEOM.

         XLOC = X(J)

         TBEST = ZEROIN (T1, T2, GETT, TOLER, IER)
         IF (IER .NE. 0) GO TO 999

C  *     Evaluate spline TGEOM, XGEOM for derivatives at TBEST:

         CALL CSDVAL (NGEOM, TGEOM, XGEOM, 1, TBEST, PSCOFS (1,1,1),
     >                PSCOFS (1,2,1), PSCOFS (1,3,1), XLOC, XD1, XD2)

C  *     Evaluate spline TGEOM, YGEOM for Y and its derivatives at TBEST:

         CALL CSDVAL (NGEOM, TGEOM, YGEOM, 1, TBEST, PSCOFS (1,1,2),
     >                PSCOFS (1,2,2), PSCOFS (1,3,2), YJ, YD1, YD2)

C  *     Evaluate derivatives using chain rule.  (Do it in the order
C        y", y', y because derivatives may not be needed - see NOTES.)

         IF (ABS (XD1) .LT. EPSIL) XD1 = SIGN (EPSIL, XD1)

C****    YPP(J) = (XD1 * YD2 - YD1 * XD2) / XD1 ** 3
         YPP(J) = (YD2 - YD1 * XD2 / XD1) / XD1 ** 2
         YP(J)  = YD1 / XD1
         Y(J)   = YJ
  100 CONTINUE

  999 RETURN
      END
C+---------------------------------------------------------------------
C
      FUNCTION GETT (T)
C
C  ACRONYM:  GET T corresponding to given X
C            --- -
C  PURPOSE:
C             GETT is used by PSEVAL as part of evaluating a parametric
C          spline as fit by PSFIT, for some X.  It evaluates
C                              F = XLOC - XTEST
C          where XLOC is an arbitrary X, and XTEST is the ordinate cor-
C          responding to the abscissa T of the spline for TGEOM & XGEOM
C          (TGEOM & PSCURV actually, since these must be passed through
C          COMMON).
C
C  METHOD:
C             * Call CSDVAL to evaluate spline for TGEOM & PSCURV at T,
C               with XTEST as output.
C             * Set GETT = XLOC - XTEST.
C
C             Eventually the zero-finder ZEROIN will find a T such that
C          this difference is arbitrarily small, meaning a value of the
C          parametric variable for the given abscissa has been found.
C
C  ARGUMENTS:
C   ARG    TYPE  I/O/S   DIM   DESCRIPTION
C    T      R    I        -    Used as test abscissa of spline for TGEOM
C                              & PSCURV.  Desired ordinate corresponding
C                              to T is XLOC.
C
C  PROCEDURES:
C   CSDVAL  Spline evaluation with optional derivatives
C
C  INTERNAL COMMON BLOCK:
C           /PSWORK/ Work space set up by PSFIT/PSEVAL.
C           Row dimensions are set to MAXPS by PARAMETER statement.
C   VAR    TYPE I/O/S DIM  DESCRIPTION
C   PSCURV  R     I MAXPS  Copy of original data abscissas, XGEOM as
C                          passed to PSFIT, for FUNCTION GETT, which
C                          is an argument of function  ZEROIN and is
C                          itself limited to a single argument.
C   TGEOM   R     I MAXPS  Values of the parametric variable at the
C                          original data points.
C   PSCOFS  R     I MAXPS  Cubic spline coefficients from CSFIT.
C                   *3*2   An array of size MAXPS*3 is required for
C                          XGEOM and another for YGEOM.
C   XLOC    R     I   -    Copy of abscissa at which PSEVAL is seeking
C                          an ordinate.
C   NCURV   I     I   -    Number of elements defining the curve.
C
C  ENVIRONMENT:  VAX/VMS, FORTRAN 77
C
C  HISTORY:
C     Jan.  83    RCL    Original design and coding.
C     Sept. 84    RCL    Now uses CSDVAL instead of IMSL's ICSEVU.
C
C  AUTHOR: Rosalie Lefkowitz, Sterling Software/NASA Ames, Moffett Field, CA
C
C----------------------------------------------------------------------

      IMPLICIT NONE

C  *  Arguments:

      REAL
     >   GETT, T

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

      REAL
     >   DUM, XTEST

C  *  Execution:

C  *  Evaluate spline TGEOM, PSCURV at T to get trial X, XTEST:

      CALL CSDVAL (NCURV, TGEOM, PSCURV, 1, T, PSCOFS (1,1,1),
     >             PSCOFS (1,2,1), PSCOFS (1,3,1), XTEST, DUM, DUM)

C  *  Compare XTEST with the target XLOC:

      GETT = XLOC - XTEST

      RETURN
      END
C+---------------------------------------------------------------------
C
      FUNCTION TSUBJ (J)
C
C  TWO-LINER: Find T (J), where T is the parametric variable defined at
C             points (X (J), Y (J)) (probably cumulative chord length).
C
C  PURPOSE:
C             TSUBJ returns the value of the parametric variable set up
C          by PSFIT for the given point J.   It may be used to find the
C          interval of this variable needed by PSEVAL.
C
C  METHOD:
C             Simply use the internal COMMON set up by PSFIT.  The user
C          is thereby isolated from this COMMON.
C
C  ARGUMENTS:
C   ARG  TYPE  I/O/S  DIM  DESCRIPTION
C    J    I      I     -   Subscript of a point known to PSFIT at which
C                          the value of the parametric variable set  up
C                          by PSFIT is required.
C
C  INTERNAL COMMON BLOCK:
C           /PSWORK/ Work space set up by PSFIT/PSEVAL, q.v.
C
C  ENVIRONMENT:  VAX/VMS, FORTRAN 77
C
C  HISTORY:
C     Sept. 85  DAS  Introduced when cumulative chord length was
C                    substituted for index number as the parametric
C                    variable used by PSFIT/PSEVAL.
C
C  AUTHOR: David Saunders, Sterling Software/NASA Ames, Moffett Field, CA
C
C----------------------------------------------------------------------

      IMPLICIT NONE

C  *  Arguments:

      INTEGER
     >   J
      REAL
     >   TSUBJ

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

C  *  Execution:

      TSUBJ = TGEOM (J)

      RETURN
      END
