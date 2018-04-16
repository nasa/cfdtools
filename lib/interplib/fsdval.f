C+----------------------------------------------------------------------
C
      SUBROUTINE FSDVAL ( MODE, NX, X, N, A, B, WRK, ION, F, FP, FPP )
C
C ACRONYM: Fourier Series partial sum Derivative eVALuation
C          -       -                  -           ---
C PURPOSE: FSDVAL  evaluates  the  finite Fourier series  represented by
C
C             f(x) = a(0)/2 + SUM <j=1 to N> a(j)*cos(jx) + b(j)*sin(jx)
C
C          for each abscissa x given.  It also evaluates the 1st and 2nd
C          derivatives at each x.  (Use FSEVAL if you don't need these.)
C          Note that b(0) = 0 is required as input.
C
C          This version handles odd or even series efficiently.  It also
C          permits turning off calculations for any of F, FP, and FPP.
C
C METHOD:  Given the recurrence-relation method used by FSEVAL (far from
C          obvious), it is simplest to interpret the derivatives as fin-
C          ite Fourier series in their own right, with coefficients dir-
C          ectly related to those of the original series.   This permits
C          reuse of FSEVAL for each derivative, at the expense of 2N + 1
C          words of work-space needed to set up the derived coefs. (Note
C          that a(0) = 0. for each of the derived series.)
C
C PARAMETERS:
C    ARG    DIM  TYPE I/O/S DESCRIPTION
C    MODE    -     I    I   MODE = 1 means series is odd:  all a(j) = 0.
C                           MODE = 2 means series is even: all b(j) = 0.
C                           MODE = 3 means series is general.
C    NX      -     I    I   Number of abscissas at which  Fourier series
C                           is to be evaluated.  NX >= 1.
C    X      NX     R    I   Abscissas requiring evaluation at (radians).
C    N       -     I    I   Defines number of terms in the series, >= 0.
C    A      0:N    R    I   Coefficients a(j), j = 0:N.
C    B      0:N    R    I   Coefficients b(j), j = 0:N, with b(0) = 0.
C    WRK   2*(N+1) R    S   Stores the modified coefficients needed  for
C                           the derivative series.
C    ION     -     I    I   3-digit decimal integer indicating which  of
C                           F, FP, FPP are requested (left-to-right dig-
C                           its respectively). E.g., ION = 110 turns off
C                           the 2nd-derivative calculations: 1 means on.
C    F      NX     R    O   F(I) = Fourier series evaluated at X(I).
C    FP     NX     R    O   FP(I) = 1st derivative at X(I).
C    FPP    NX     R    O   FPP(I) = 2nd derivative at X(I).
C
C EXTERNAL REFERENCES:  FSEVAL  (used once for each of F, FP, and FPP).
C
C ENVIRONMENT: FORTRAN 77
C
C DEVELOPMENT HISTORY:
C       DATE   INITIALS   DESCRIPTION
C    07/08/83    DAS      Original design and code.
C    03/04/87    DAS      B(0) = 0 is expected now.
C
C AUTHOR: David Saunders, Informatics General, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C ... Arguments:

      INTEGER  ION, MODE, N, NX
      REAL     A(0:N), B(0:N), WRK(0:N,2), X(NX), F(NX), FP(NX), FPP(NX)

C ... Local variables:

      INTEGER  I, I0
      LOGICAL  EVEN, ODD
      
C ... Execution:

      I0 = ION/100
      IF ( I0 .GT. 0 ) CALL FSEVAL ( MODE, NX, X, N, A, B, F )

      ODD  = MODE.EQ.1
      EVEN = MODE.EQ.2
      WRK(0,1) = 0.E+0
      WRK(0,2) = 0.E+0

      IF ( ( ION - I0*100 )/10 .GT. 0 ) THEN

C ...    Set up coefficients of series representing 1st derivative:

         DO 20 I = 1, N
            IF ( .NOT. EVEN ) WRK(I,1) =  I * B(I)
            IF ( .NOT. ODD  ) WRK(I,2) = -I * A(I)
   20    CONTINUE

         CALL FSEVAL ( 3-MOD(MODE,3), NX, X, N, WRK(0,1), WRK(0,2), FP )

      END IF

      IF ( MOD ( ION, 10 ) .GT. 0 ) THEN

C ...    Set up coefficients of series representing 2nd derivative:

         DO 30 I = 1, N
            IF ( .NOT. ODD  ) WRK(I,1) = (-I*I) * A(I)
            IF ( .NOT. EVEN ) WRK(I,2) = (-I*I) * B(I)
   30    CONTINUE

         CALL FSEVAL ( MODE, NX, X, N, WRK(0,1), WRK(0,2), FPP )

      END IF

      RETURN
      END
