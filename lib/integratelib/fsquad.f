C+----------------------------------------------------------------------
C
      SUBROUTINE FSQUAD ( MODE, N, A, B, X1, X2, WRK, AREA )
C
C ACRONYM: Fourier Series QUADrature
C          -       -      ----
C PURPOSE: FSQUAD integrates a function given in the finite Fourier ser-
C          ies form
C
C             f(x) = a(0)/2 + SUM <j=1 to N> a(j)*cos(jx) + b(j)*sin(jx)
C
C          over the interval [X1,X2] (radians).    Note that b(0) = 0 is
C          required as input.  Note also that if [X1,X2] = [0,2*pi] then
C          the integral is just a(0) * pi and a lot of calculation would
C          be wasted.
C
C          This version handles odd or even series efficiently.
C
C METHOD:  Given the recurrence-relation method used by FSEVAL (far from
C          obvious),  the indefinite integral of a finite Fourier series
C          can be interpreted as a related series in its own right. This
C          permits reuse of FSEVAL for each of two evaluations at X1 and
C          X2, with a subtraction giving the desired "area."
C
C PARAMETERS:
C    ARG    DIM  TYPE I/O/S DESCRIPTION
C    MODE    -     I    I   MODE = 1 means series is odd:  all a(j) = 0.
C                           MODE = 2 means series is even: all b(j) = 0.
C                           MODE = 3 means series is general.
C    N       -     I    I   Defines number of terms in the series, >= 0.
C    A      0:N    R    I   Coefficients a(j), j = 0:N.
C    B      0:N    R    I   Coefficients b(j), j = 0:N, with b(0) = 0.
C    X1      -     R    I   Lower limit for definite integral (radians).
C    X2      -     R    I   Upper limit for definite integral (radians).
C                           See note above if [X1,X2] = [0,2*pi].
C    WRK   2*(N+1) R    S   Stores the modified coefficients needed for
C                           the derived series.
C    AREA    -     R    O   Desired definite integral.
C
C EXTERNAL REFERENCES:  FSEVAL  (sum of series by recurrence relation)
C
C ENVIRONMENT: FORTRAN 77
C
C DEVELOPMENT HISTORY:
C       DATE   INITIALS   DESCRIPTION
C    03/13/87    DAS      Adapted from FSDVAL (same idea).
C
C AUTHOR: David Saunders, Sterling Software, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER  MODE, N
      REAL     A(0:N), B(0:N), X1, X2, WRK(0:N,2), AREA

C     Local variables:

      INTEGER  I, NUMODE
      REAL     Q1, Q2

C     Set up the new series representing the indefinite integral.

      WRK(0,1) = 0.E+0
      WRK(0,2) = 0.E+0
      NUMODE = 3 - MOD ( MODE, 3 )

      IF ( NUMODE .NE. 1 ) THEN  ! New series is not odd, so set up new a(*):

         DO 20 I = 1, N
            WRK(I,1) = -B(I) / I
   20    CONTINUE

      END IF

      IF ( NUMODE .NE. 2 ) THEN  ! New series is not even, so set up new b(*):

         DO 30 I = 1, N
            WRK(I,2) = A(I) / I
   30    CONTINUE

      END IF

      CALL FSEVAL ( NUMODE, 1, X1, N, WRK(0,1), WRK(0,2), Q1 )
      CALL FSEVAL ( NUMODE, 1, X2, N, WRK(0,1), WRK(0,2), Q2 )

      AREA = 0.5E+0 * A(0) * ( X2 - X1 ) + ( Q2 - Q1 )

      RETURN
      END
