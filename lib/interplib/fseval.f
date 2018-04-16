C+----------------------------------------------------------------------
C
      SUBROUTINE FSEVAL ( MODE, NX, X, N, A, B, F )
C
C ACRONYM: Finite Fourier Series EVALuation
C                 -       -      ----
C PURPOSE: FSEVAL  evaluates  the  finite Fourier series  represented by
C
C          f(x) = a(0)/2 + SUM <j=1 to N> a(j)*cos(jx) + b(j)*sin(jx)
C
C          for each abscissa x given.    It is required that b(0) = 0 be
C          given to FSEVAL.
C
C          This version handles series that are odd,  even,  or neither,
C          efficiently. If derivatives of the series are needed too, use
C          FSDVAL.
C
C METHOD:  Acton's "Numerical Methods That Work" (1970), pp.11-12, shows
C          how to evaluate sums of orthogonal functions efficiently with
C          the use of recurrence relations. The above sum degenerates to
C
C          f(x) = a(0)/2 + A(1)*cos(x) - A(2) + B(1)*sin(x)        where
C
C          A(j) = 2*cos(x) * A(j+1)  -  A(j+2)  +  a(j)              for
C
C                 j = N, N-1, ..., 1,  with  A(N+1) = A(N+2) = 0.
C
C          The B(j) are defined analogously.  Thus only one sine and one
C          cosine evaluation are really needed, regardless of N.
C
C          MODE was introduced so that FSEVAL can be used to compute the
C          coefficients of finite Fourier series   (which are themselves
C          sums of cosine or sine terms  when simple quadrature formulae
C          are used to estimate the integrals involved).
C
C PARAMETERS:
C    ARG    DIM  TYPE I/O/S DESCRIPTION
C    MODE    -     I    I   MODE = 1 means series is odd:  all a(j) = 0.
C                                = 2 means series is even: all b(j) = 0.
C                                = 3 means series is general.
C    NX      -     I    I   Number of abscissas at which series is to be
C                           evaluated.  NX >= 1.
C    X      NX     R    I   Abscissas requiring evaluation at (radians).
C    N       -     I    I   Defines number of terms in the series, >= 0.
C    A      0:N    R    I   Coefficients a(j), j = 0:N.
C    B      0:N    R    I   Coefficients b(j), j = 0:N,  where b(0) = 0.
C    F      NX     R    O   F(I) = Fourier series evaluated at X(I).
C
C ENVIRONMENT: FORTRAN 77
C
C DEVELOPMENT HISTORY:
C       DATE   INITIALS   DESCRIPTION
C    07/08/83    DAS      Initial design and code.
C    06/05/86    RAK      Removed non-standard tabs and DO/END DO loops,
C                         and added non-standard IMPLICIT NONE.
C    03/03/87    DAS/RAK  Should work for N=0 - it does now.
C
C AUTHOR: David Saunders, Informatics General, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER  I, J, MODE, N, NX
      REAL     A(0:N), B(0:N), F(NX), X(NX)
      REAL     AJ, AJP1, AJP2, COSX, EVENS, ODDS

      AJ = 0.E+0                                ! For the N = 0 case.

      DO 30, I = 1, NX
         COSX = 2.E+0 * COS( X(I) )

         IF ( MODE.EQ.1 ) THEN

C ...       Series is odd - no cosine terms, so:

            EVENS = 0.E+0
         ELSE
            AJP2 = 0.E+0
            AJP1 = 0.E+0
            DO 10, J = N, 1, -1
               AJ   = COSX * AJP1  -  AJP2  +  A(J)
               AJP2 = AJP1
               AJP1 = AJ
   10       CONTINUE
            EVENS = ( A(0) + AJ*COSX )*0.5E+0  -  AJP2
         END IF

         IF ( MODE.EQ.2 ) THEN

C ...         Series is even - no sine terms, so:

            ODDS = 0.E+0
         ELSE
            AJP2 = 0.E+0
            AJP1 = 0.E+0

            DO 20, J = N, 1, -1
               AJ   = COSX * AJP1  -  AJP2  +  B(J)
               AJP2 = AJP1
               AJP1 = AJ
   20       CONTINUE
            ODDS = AJ * SIN( X(I) )
         END IF

         F(I) = ODDS + EVENS
   30 CONTINUE

      RETURN
      END
