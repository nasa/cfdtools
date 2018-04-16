C+----------------------------------------------------------------------
C
      SUBROUTINE FSCOEF ( MODE, M, F, J, AJ, BJ, IER )
C
C ACRONYM: Finite Fourier Series COEFficient calculation
C                 -       -      ----
C PURPOSE: FSCOEF  calculates the  jth  coefficients or harmonics in the
C          finite Fourier series represented by
C
C             f(x) = a(0)/2 + SUM <j=1 to N> a(j)*cos(jx) + b(j)*sin(jx)
C
C          where 2*N+1 <= M and M is the number of data points which the
C          series is representing.  Also, b(0) = 0 is assigned here.
C
C METHOD:  Assume that f(x) is available as values  f(i) = f( x(i) )  at
C          uniformly-spaced points  x(i) = 0, h, 2h, ..., mh (h=2*PI/m).
C          Then the coefficients a(j), b(j) are given by
C
C             a(j) = (2/m) * SUM <k=0 to m-1> f(k)*cos(k*jh)   j=0,1,...
C
C             b(j) = (2/m) * SUM <k=1 to m-1> f(k)*sin(k*jh)   j=1,2,...
C
C          for j <= (m-1)/2, with b(j) = 0.
C
C          These are the same as what the trapezoidal rule on  the  f(i)
C          gives for the integrals associated with the infinite  Fourier
C          series,  if  f(m) = f(0)  is assumed.   However, the sums are
C          carried out over the full cycle with no symmetry assumptions,
C          and f(m) is not used.    Passing the average of the values at
C          x=0 and x=2*PI for f(0) might be a good idea.
C
C          Evaluation of  a(j)  and  b(j) is done here by observing that
C          the sums can be interpreted as (even and odd)  finite Fourier
C          series in their own right. Thus FSEVAL is used appropriately.
C
C PARAMETERS:
C    ARG    DIM  TYPE I/O/S DESCRIPTION
C    MODE    -     I    I   MODE = 1 means series is odd:  all a(j) = 0.
C                           MODE = 2 means series is even: all b(j) = 0.
C                           MODE = 3 means series is general.
C    M       -     I    I   Defines the discretization of f(x),  at uni-
C                           formly spaced points x=0,h,...,mh (h=2*PI/m)
C    F     0:M-1   R    I   Discretized function.  F(M) is not used.
C    J       -     I    I   Defines the coefs. required. J <= (M-1)/2.
C    AJ      -     R    O   Coefficient a(j) for given j; AJ=0 if MODE=1
C    BJ      -     R    O   Coefficient b(j) for given j; BJ=0 if MODE=2
C    IER     -     I    O   IER = 0 means no input errors were noticed.
C                           IER = 2 means J is not meaningful (J < 0, or
C                                   J > (M-1)/2.
C ENVIRONMENT: FORTRAN 77
C
C DEVELOPMENT HISTORY:
C       DATE   INITIALS   DESCRIPTION
C    07/09/83    DAS      Initial implementation (b(0) an exception).
C    03/04/87    DAS      BJ = 0 is returned now if J = 0.
C
C AUTHOR: David Saunders, Informatics General, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      REAL       AJ, BJ, F(0:M-1)
      INTEGER    IER, J, M, MODE
      
      PARAMETER  ( PI = 3.1415927 )


      IF ( J. LT. 0  .OR.  J .GT. (M-1)/2 ) GO TO 90

      C  = 2.E+0/M
      H  = C*PI
      HJ = H*J

      IF ( MODE.EQ.1 ) THEN

C ...    Series is odd - no cosine terms, so:

         AJ = 0.E+0

      ELSE

         F0 = F(0)
         F(0) = F0 + F0

         CALL FSEVAL ( 2, 1, HJ, M-1, F, DUMMY, AJ )

         AJ = C*AJ
         F(0) = F0
      END IF

      IF ( MODE.EQ.2 ) THEN

C ...    Series is even - no sine terms, so:

         BJ = 0.E+0

      ELSE

         CALL FSEVAL ( 1, 1, HJ, M-1, DUMMY, F, BJ )
         BJ = C*BJ

      END IF

      IER = 0
      GO TO 99

   90 IER = 2

   99 RETURN
      END
