C+---------------------------------------------------------------------
C
      SUBROUTINE CSDVAL (NX, X, Y, NU, U, B, C, D, S, SP, SPP)
C
C     One-liner: Cubic spline evaluation on an array, with derivatives.
C     ----------
C
C     Purpose:
C     --------
C
C        CSDVAL evaluates the 1-dimensional cubic spline function
C
C     S(U) = Y(I) + B(I)*(U - X(I)) + C(I)*(U - X(I))**2 + D(I)*(U - X(I))**3
C
C     along with its first and second derivatives
C
C     SP(U)  =      B(I)          + 2*C(I)*(U - X(I))  + 3*D(I)*(U - X(I))**2
C     SPP(U) =    2*C(I)          + 6*D(I)*(U - X(I))
C
C     using Horner's Rule for polynomials.  If U < X(1), then I = 1 is used;
C     if U > X(NX) then I = NX-1 is used. The data must be monotonic, but
C     may be increasing or decreasing.  (The coefficients associated with
C     the NXth point are not required.)  Most of the effort is in deciding
C     which set of spline coefficients to use - an "interpolation search"
C     is employed which should be highly efficient if the knots are fairly
C     uniformly spaced (commonly the case).
C
C        CSDVAL was developed as a companion routine to CSFIT, but it should
C     work with any method of generating cubic spline coefficients.  CSEVAL
C     should be used if derivatives are not required.  See also NOTES.
C
C     Arguments:
C     ----------
C
C     Name     Dimension  Type  I/O/S  Description
C     NX                   I    I      Number of data points defining the
C                                      spline ("knots"); >= 2.
C
C     X        NX          R    I      Abscissas of knots. Must be monotone
C                                      increasing or decreasing (consistent
C                                      with the data used to fit the spline).
C
C     Y        NX          R    I      Ordinates of knots.
C
C     NU                   I    I      Number of points at which to evaluate
C                                      the spline.
C
C     U        NU          R    I      Abscissas at which to evaluate spline.
C
C     B,C,D    NX          R    I      Spline coefficients, e.g. as computed
C                                      by CSFIT.
C
C     S       NU           R      O    Spline values at points U(I).
C
C     SP      NU           R      O    1st derivative values at points U(I).
C
C     SPP     NU           R      O    2nd derivative values at points U(I).
C
C     External modules:
C     -----------------
C
C     INTERVAL  Interpolation search for interval containing a point.
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN 77 V4.7
C     ------------
C
C     Notes:
C     ------
C
C     (1)  IMPLICIT NONE and eight character symbols are not (yet) standard.
C
C     (2)  If derivatives are not desired, use CSEVAL instead. Or, take
C          advantage of the fact that the returned values are computed in
C          the order  SPP, SP, S so that if only first derivative is needed,
C          for example, the same array can be passed for SPP as for SP.
C
C     Bibliography:
C     -------------
C
C     (1)  Forsythe, G., M. Malcolm, and C. Moler.  Computer Methods for
C             Mathematical Computations.  Englewood Cliffs: Prentice-Hall,
C             1977.  (Chapter 4)
C
C     Author:  Robert Kennelly, Informatics General Corp.
C     -------
C
C     History:
C     --------
C
C      8/27/84    RGL    Adapted from CSEVAL to provide derivatives.
C     10/20/87    RAK    Abstract the interval search (with revisions) as
C                        a separate module.
C     06/09/91    DAS    Ensured S (U) = Y (NX) if U = X (NX).  (This
C                        is not guaranteed by CSDVAL's use of the cubic
C                        in the (NX-1)th interval for such a U, except
C                        if the arithmetic is exact, which it isn't.)
C
C----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     &   ONE, TWO, THREE, SIX
      PARAMETER
     &  (ONE   = 1.0E+0,
     &   TWO   = 2.0E+0,
     &   THREE = 3.0E+0,
     &   SIX   = 6.0E+0)

C     Arguments.

      INTEGER
     &   NU, NX
      REAL
     &   B (NX), C (NX), D (NX), S (NU), SP (NU), SPP (NU), U (NU),
     &   X (NX), Y (NX)

C     Local variables.

      INTEGER
     &   IU, LEFT
      REAL
     &   ARROW, DX, XRIGHT

C     Execution.
C     ----------

      ARROW = SIGN (ONE, X (2) - X (1))
      XRIGHT = X (NX)
      LEFT = 1

      DO 10, IU = 1, NU

C        Search for the best "left-hand" endpoint to use for interpolation.

         CALL INTERVAL (NX, X, U (IU), ARROW, LEFT)

C        Evaluate the spline.  Note that if U is off-scale on the left,
C        the coefficients of the first interval are used for extrapolation,
C        while if U is greater than X(NX), the (NX-1)th set is used.
C        This means when U = X (NX), roundoff could mean S (U) is not
C        exactly Y (NX).  Hence the test for this special case.
C        Note also the order, in case y" is not really needed - see NOTES.

         DX = U (IU) - X (LEFT)

         SPP (IU) = TWO * C (LEFT) +
     &       (DX * (SIX * D (LEFT)))

         SP  (IU) = B (LEFT) +
     &       (DX * (TWO * C (LEFT) +
     &       (DX * (THREE * D (LEFT)))))

         S   (IU) = Y (LEFT) +
     &       (DX * (B (LEFT) +
     &       (DX * (C (LEFT) +
     &       (DX * (D (LEFT)))))))

         IF (U (IU) .EQ. XRIGHT) S (IU) = Y (NX)
   10 CONTINUE

C     Termination.
C     ------------

      RETURN
      END
