C+---------------------------------------------------------------------
C
      SUBROUTINE CSEVAL (NDATA, X, Y, NU, U, B, C, D, S)
C
C     One-liner: Cubic spline evaluation on an array.
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        CSEVAL evaluates the 1-dimensional cubic spline function
C
C     S(U) = Y(I) + B(I)*(U - X(I)) + C(I)*(U - X(I))**2 + D(I)*(U - X(I))**3
C
C     using Horner's Rule for polynomials. Normally, if U < X(1), then I = 1
C     is used, and if U > X(NX) then I = NX - 1 is used.  However, this version
C     treats periodic data properly - see the usage of NDATA.
C
C        The data must be monotonic, but may be increasing or decreasing.
C     (The coefficients associated with the NXth point are not required.)
C     Most of the effort is in deciding which set of spline coefficients to
C     use.  An "interpolation" search is employed which should be highly
C     efficient if the knots are fairly uniformly spaced (commonly the case).
C
C        CSEVAL was developed as a companion routine to CSFIT, but it should
C     work with any method of generating cubic spline coefficients.  CSDVAL
C     should be used if derivatives of the spline are required.
C
C     Arguments:
C     ----------
C
C     Name     Dimension  Type  I/O/S  Description
C     NDATA                I    I      NX = |NDATA| is the number of data points
C                                      defining the spline ("knots");
C                                      use NDATA = -NX to signify the periodic
C                                      data case; NX >= 2.
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
C     External modules:
C     -----------------
C
C        INTERVAL   Interpolation search for interval containing a point.
C
C     Environment:  FORTRAN 77 + IMPLICIT NONE
C     ------------
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
C      7/27/84    RAK    Original design and coding, based in part on
C                        SEVAL from Forsythe, Malcolm, and Moler, and
C                        using a search method adapted from Sedgewick.
C                        An earlier rewrite of SEVAL by Trosin included
C                        conversion from a FUNCTION to a SUBROUTINE to
C                        allow more than one evaluation per call.
C     10/19/87    RAK    Abstracted the interval search (with revisions)
C                        as a separate module.
C     06/09/91    DAS    Ensured S (U) = Y (NX) if U = X (NX).  (This
C                        is not guaranteed by CSEVAL's use of the cubic
C                        in the (NX-1)th interval for such a U, except
C                        if the arithmetic is exact, which it isn't.)
C     08/16/97     "     Handled the periodic case via NX < 0 kludge.
C
C----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     &   ONE
      PARAMETER
     &  (ONE = 1.0E+0)

C     Arguments.

      INTEGER
     &   NDATA, NU
      REAL
     &   B (*), C (*), D (*), S (NU), U (NU), X (*), Y (*)

C     Local variables.

      INTEGER
     &   IU, LEFT, NX
      REAL
     &   ARROW, DX, PERIOD, XEVAL, XLEFT, XRIGHT

C     Execution.
C     ----------

      ARROW  = SIGN (ONE, X (2) - X (1))
      NX     = ABS (NDATA)
      XRIGHT = X (NX)
      LEFT   = 1

      IF (NDATA .GT. 0) THEN  ! Non-periodic case

         DO IU = 1, NU

            XEVAL = U (IU)

C           Search for the best "left-hand" endpoint to use for interpolation.

            CALL INTERVAL (NX, X, XEVAL, ARROW, LEFT)

C           Evaluate the spline.  Note that if U is off-scale on the left,
C           the coefficients of the first interval are used for extrapolation,
C           while if U is greater than X(NX), the (NX-1)th set is used.
C           This means when U = X (NX), roundoff could mean S (U) is not
C           exactly Y (NX).  Hence the test for this special case.

            DX = XEVAL - X (LEFT)

            S (IU) = Y (LEFT) +
     &        (DX * (B (LEFT) +
     &        (DX * (C (LEFT) +
     &        (DX * (D (LEFT)))))))

            IF (XEVAL .EQ. XRIGHT) S (IU) = Y (NX)

         END DO

      ELSE  ! Periodic case

         IF (ARROW .EQ. ONE) THEN
            XLEFT = X (1)
         ELSE
            XLEFT = XRIGHT
            XRIGHT = X (1)
         END IF

         PERIOD = XRIGHT - XLEFT  ! Definitely positive

         DO IU = 1, NU

            XEVAL = U (IU)

            IF (XEVAL .LT. XLEFT) THEN

               XEVAL = XRIGHT - MOD (XLEFT - XEVAL, PERIOD)

            ELSE IF (XEVAL .GT. XRIGHT) THEN

               XEVAL = XLEFT + MOD (XEVAL - XRIGHT, PERIOD)

            END IF

            CALL INTERVAL (NX, X, XEVAL, ARROW, LEFT)

            DX = XEVAL - X (LEFT)

            S (IU) = Y (LEFT) +
     &        (DX * (B (LEFT) +
     &        (DX * (C (LEFT) +
     &        (DX * (D (LEFT)))))))

            IF (XEVAL .EQ. X (NX)) S (IU) = Y (NX)

         END DO

      END IF

C     Termination.
C     ------------

      RETURN
      END
