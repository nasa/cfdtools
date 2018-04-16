C+----------------------------------------------------------------------
C
      SUBROUTINE MSFIT (N, X, Y, CYCLIC, B, C, D, IER)
C
C  PURPOSE:
C     MSFIT computes the coefficients B(I), C(I), and D(I), I=1,2,...,N
C     for the MONOTONE cubic interpolating spline defined by
C
C     S(X) = Y(I) + B(I)*(X-X(I)) + C(I)*(X-X(I))**2 + D(I)*(X-X(I))**3
C
C     for  X(I) <= X <= X(I+1).  It is derived from CSPLIN by J. Cordova
C     (Sterling Software, Oct. '86), which implements the algorithm of
C     Fritsch and Butland in SIAM J. SCI. STAT. COMPUT., VOL. 5, NO. 2,
C     June, 1984. MSFIT separates fitting the spline from evaluating it.
C     Use CSEVAL to evaluate the fit.
C
C     This class of spline ensures continuity of the first derivative at
C     the knots, but not of the second derivative.
C
C     This version provides a periodic (cyclic) option, needed to close
C     a curve smoothly as part of a parametric fit - see PSFIT.  (But
C     PLSFIT or PLSINTRP may be preferable for parametric fits.)
C
C     See also the LCSFIT implementation (Y vs. X derivation of PLSFIT's
C     X vs. T, Y vs. T).  These both calculate the spline coefficients
C     as they are needed, to avoid storing them.
C
C  ARGUMENTS:
C     NAME   DIM TYPE I/O/S DESCRIPTION
C     N       -   I     I   Number of data points, or knots.  N >= 2.
C     X       N   R     I   Abscissas of the knots, strictly increasing
C                           or strictly decreasing.
C     Y       N   R     I   Ordinates of the knots; Y(N)=Y(1) if CYCLIC.
C     CYCLIC  -   L     I   .TRUE. means periodic case; Y(N)=Y(1) reqd.;
C                           .FALSE. otherwise.
C     B,C,D   N   R     O   Spline coefficients (see PURPOSE and NOTES).
C     IER     -   I     O   0: No errors were detected;
C                           1: Too few data points; N < 2;
C                           4: Cyclic mode requested but Y(N) .NE. Y(1).
C  PROCEDURES:
C     PCHST   Function returning +1, -1, or 0. (No obvious way to make it
C             a statement function.  Source is in same module as MSFIT.)
C
C  NOTES:
C         Y(I) = S (X(I))
C         B(I) = S' (X(I))
C         C(I) = S'' (X(I))/2
C         D(I) = S''' (X(I))/6
C
C     Note that C (*) and D (*) are used for intermediate results before
C     their final definitions, to avoid local storage.
C
C  ENVIRONMENT:  VAX/VMS; FORTRAN 77
C                IMPLICIT NONE is an extension of FORTRAN 77.
C
C  HISTORY:
C  10/26/86  J. Cordova   Fit + evaluation, as routine CSPLIN
C  12/23/86  D. Saunders  Fit only, as MSFIT, for use with CSEVAL; avoid
C                         local storage; CSFIT-type nomenclature
C  01/14/87  D. Saunders  Introduced CYCLIC argument, for use by PSFIT.
C  02/27/87    "    "     Revised nomenclature following Kennelly's
C                         PLSFIT and the PCHIP package of Fritsch, et al.
C  10/23/87    "    "     Descending abscissas permissible now in light
C                         of revised CSEVAL.
C  06/21/91    "    "     References made to later related subroutines.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   N, IER
      REAL
     >   X (N), Y (N), B (N), C (N), D (N)
      LOGICAL
     >   CYCLIC

C     Local constants:

      REAL
     >   ZERO, ONE, TWO, THREE, THIRD
      PARAMETER
     >  (ZERO  = 0.E+0,
     >   ONE   = 1.E+0,
     >   TWO   = 2.E+0,
     >   THREE = 3.E+0,
     >   THIRD = ONE/THREE)

C     Local variables:

      INTEGER
     >   I, ILAST
      REAL
     >   ALPHA, BMAX, H, TOP, W1, W2

C     Procedures:

      REAL
     >   PCHST
      EXTERNAL
     >   PCHST

C     Execution:

      IER = 0
      IF (N .LT. 2) THEN
         IER = 1
         GO TO 999
      END IF

      IF (CYCLIC) THEN
         IF (Y (N) .NE. Y (1)) THEN
            IER = 4
            GO TO 999
         END IF         
      END IF         

      DO 100, I = 1, N-1
         C (I) = X (I+1) - X (I)
         D (I) = (Y (I+1) - Y (I)) / C (I)
  100 CONTINUE

C     Linear interpolation for the N = 2 case.

      IF (N .EQ. 2) THEN
         B (1) = D (1)
         B (2) = D (1)
         C (1) = ZERO
         C (2) = ZERO
         D (1) = ZERO
         D (2) = ZERO
         GO TO 999
      END IF

C     Monotone cubic spline for N >= 3 case.

      IF (CYCLIC) THEN

C        Periodic case.  B (N) is like other B (*) if ...

         ILAST = N
         C (N) = C (1)
         D (N) = D (1)
      ELSE

C        First boundary point: 3-pt formula modified to preserve shape.

         ILAST = N - 1
         H = C (1) + C (2)
         W2 = -C (1) / H
         W1 = ONE - W2
         B (1) = W1 * D (1) + W2 * D (2)

         IF (PCHST (B (1), D (1)) .LE. ZERO) THEN
            B (1) = ZERO
         ELSE IF (PCHST (D (1), D (2)) .LT. ZERO) THEN
            BMAX = THREE * D (1)
            IF (ABS (B (1)) .GT. ABS (BMAX)) B (1) = BMAX
         END IF

C        Last boundary point, as for first.

         H = C (N-2) + C (N-1)
         W1 = -C (N-1) / H
         W2 = ONE - W2
         B (N) = W1 * D (N-2) + W2 * D (N-1)

         IF (PCHST (B (N), D (N-1)) .LE. ZERO) THEN
            B (N) = ZERO
         ELSE IF (PCHST (D (N-2), D (N-1)) .LT. ZERO) THEN
            BMAX = THREE * D (N-1)
            IF (ABS (B (N)) .GT. ABS (BMAX)) B (N) = BMAX
         END IF
      END IF

C     Interior points (Brodlie modification of Butland formula).
C     (This loop is vectorizable using SIGN for TOP and adding an
C     epsilon to the denominator, but at a cost in clarity.)

      DO 200, I = 2, ILAST
         TOP = D (I-1) * D (I)
         IF (TOP .GT. ZERO ) THEN          ! Slopes are same sign
            ALPHA = THIRD * (ONE + C (I) / (C (I-1) + C (I)))
            B (I) = TOP / (ALPHA * D (I) + (ONE - ALPHA) * D (I-1))
         ELSE
            B (I) = ZERO                   ! Make this point an extremum
         END IF
  200 CONTINUE

      IF (CYCLIC) THEN
         B (1) = B (N)
      END IF         

C     Quadratic and cubic coefficients for intervals 1:N-1.

      DO 300, I = 1, N-1
         H = ONE / C (I)
         C (I) = (THREE * D (I) - TWO * B (I) - B (I+1)) * H
         D (I) = (B (I) + B (I+1) - TWO * D (I)) * H * H
  300 CONTINUE

      IF (CYCLIC) THEN

C        Coefs. for interval N are same as coefs. for interval 1,
C        which in turn use interval 2.

         H = ONE / C (N)
         C (N) = (THREE * D (N) - TWO * B (N) - B (2)) * H
         D (N) = (B (N) + B (2) - TWO * D (N)) * H * H
      ELSE

C        Coefs. for X >= X (N) are related to derivatives at X (N),
C        using coefs. for interval N-1.

         H = X (N) - X (N-1)
         B (N) = B (N-1) + C (N-1) * TWO * H + D (N-1) * THREE * H * H
         C (N) = C (N-1) + D (N-1) * THREE * H
         D (N) = D (N-1)
      END IF

  999 RETURN
      END
C+----------------------------------------------------------------------
C
      FUNCTION PCHST (ARG1, ARG2)
C
C  PURPOSE:
C     PCHST returns +1., -1., or 0. according to the signs of its two
C     arguments. It is convenient for monotonic spline routine MSFIT.
C     Nomenclature is that of PCHIP (Fritsch, et al.).
C
C-----------------------------------------------------------------------

      REAL
     >   ARG1, ARG2, ONE, PCHST, ZERO

      PARAMETER
     >   (ONE = 1.E+0, ZERO = 0.E+0)

      IF (ARG1 * ARG2 .EQ. ZERO) THEN
         PCHST = ZERO
      ELSE
         PCHST = SIGN (ONE, ARG1) * SIGN (ONE, ARG2)
      END IF

      END
