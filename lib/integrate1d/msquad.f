C+----------------------------------------------------------------------
C
      SUBROUTINE MSQUAD (N, X, F, OFFSET, B, C, D, FINT)
C
C
C     Description and usage:
C
C           MSQUAD (Monotonic Spline QUADrature) integrates a function
C        using a previously calculated cubic spline representation.  In
C        contrast to CSQUAD, it uses the obvious quartic formula in each
C        interval, so that any cubic spline coefficients (including those
C        from, say, MSFIT) apply.  (CSQUAD takes advantage of the special
C        relationship between the coefficients of a conventional cubic
C        spline from, say, CSFIT, to derive the integral from just one of
C        the three sets of coefficients.  CSQUAD is not appropriate for
C        local-spline methods.)

C           The (definite) integral up to each knot is returned, with an
C        optional constant added to the value of the integral array at the
C        first point.  Output array FINT may share storage with any of the
C        spline coefficient arrays B, C, D if desired.
C
C
C     Arguments:
C
C        Name    Dimension  Type  I/O/S  Description
C         N                  I    I      Number of knots.
C         X        N         R    I      Abscissas (spline knots Xi).
C         F        N         R    I      Function values.
C         OFFSET             R    I      Value assigned to the first element
C                                        of FINT (will be added to all).
C         B, C, D  N         R    I      Spline coefficients from MSFIT (or
C                                        CSFIT or others) in the order 1st,
C                                        2nd, and 3rd derivative-like resp.
C         FINT     N         R      O    Integral of F up to each knot.
C                                        FINT (1) = OFFSET; FINT (N) =
C                                        OFFSET + integral from X (1) to X (N).
C
C
C     Standards violations:  IMPLICIT NONE is non-standard.
C
C
C     Environment:  DEC VAX-11/780   VMS/V4.4   FORTRAN
C
C
C     Authors:  David Saunders/Robert Kennelly, Sterling Software/NASA-Ames
C
C
C     History:
C
C        25 Aug. 1987  DAS/RAK  Adapted from CSQUAD when it was realized
C                               that CSQUAD does not apply to MSFIT's coefs.
C
C-----------------------------------------------------------------------


C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     >   HALF, THIRD, FOURTH
      PARAMETER
     >  (HALF   = 1.0E+0 / 2.0E+0,
     >   THIRD  = 1.0E+0 / 3.0E+0,
     >   FOURTH = 1.0E+0 / 4.0E+0)

C     Arguments.

      INTEGER
     >   N
      REAL
     >   B (N), C (N), D (N), F (N), FINT (N), OFFSET, X (N)

C     Local variables.

      INTEGER
     >   I
      REAL
     >   H


C     Execution.
C     ----------

C     Compute integral on each subinterval.  Work backwards to avoid
C     overwriting coefficients if one of those vectors was also passed
C     in as FINT.

      DO 10, I = N - 1, 1, -1
         H = X (I + 1) - X (I)
         FINT (I + 1) = H * (F (I) + H * (HALF * B (I) +
     >      H * (THIRD * C (I) + H * FOURTH * D (I))))
   10 CONTINUE

C     Sum the contributions from each subinterval.

      FINT (1) = OFFSET
      DO 20, I = 2, N
         FINT (I) = FINT (I) + FINT (I - 1)
   20 CONTINUE


C     Termination.
C     ------------

      RETURN
      END
