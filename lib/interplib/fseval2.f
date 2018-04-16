C+----------------------------------------------------------------------
C
      SUBROUTINE FSEVAL2 (MODE, NEVAL, XEVAL, ORIGIN, PERIOD, N, A, B,
     >                    FEVAL)
C
C ONE-LINER:  Finite Fourier series evaluation (any coordinate system)
C
C PURPOSE:
C     FSEVAL2  evaluates  the  finite Fourier series  represented by
C                         N             _              _
C        f(x) = a(0)/2 + SUM  a(j)*cos(jx) + b(j)*sin(jx)     where
C                        j=1
C               _
C               x = 2 pi (x - origin) / period
C
C     for each abscissa x given.  The units are arbitrary, but XEVAL,
C     ORIGIN, PERIOD should be related to the data used to calculate
C     the Fourier coefficients.
C
C     This version handles series that are odd, even, or neither.
C
C METHOD:
C     Acton's "Numerical Methods That Work" (1970), pp.11-12, shows
C     how to evaluate sums of orthogonal functions efficiently with
C     the use of recurrence relations. The above sum degenerates to
C
C        f(x) = a(0)/2 + A(1)*cos(x) - A(2) + B(1)*sin(x)     where
C
C        A(j) = 2*cos(x) * A(j+1)  -  A(j+2)  +  a(j)           for
C
C               j = N, N-1,..., 1,  with  A(N+1) = A(N+2) = 0.
C
C     The B(j) are defined analogously.  Thus only one sine and one
C     cosine evaluation are really needed, regardless of N.
C
C     MODE was introduced so that FSEVAL can be used to compute the
C     coefficients of finite Fourier series   (which are themselves
C     sums of cosine or sine terms  when simple quadrature formulae
C     are used to estimate the integrals involved and the abscissas
C     are uniformly spaced).
C
C ARGUMENTS:
C    ARG    DIM  TYPE I/O/S DESCRIPTION
C    MODE    -     I    I   MODE = 1 means series is odd:  all a(j) = 0.
C                                = 2 means series is even: all b(j) = 0.
C                                = 3 means series is general.
C    NEVAL   -     I    I   Number of abscissas at which series is to be
C                           evaluated.  NEVAL >= 1.
C    XEVAL  NEVAL  R    I   Abscissas requiring evaluation.
C    ORIGIN  -     R    I   Should be X(1) of original dataset.
C    PERIOD  -     R    I   Normally the original data range (but may be
C                           twice that when used with FSERIES and MODE=1).
C    N       -     I    I   Defines number of terms in the series, >= 0.
C    A      0:N    R    I   Coefficients a(j), j = 0:N.
C    B      0:N    R    I   Coefficients b(j), j = 0:N.  b(0) is ignored.
C    FEVAL  NEVAL  R    O   FEVAL(I) = Fourier series evaluated at XEVAL(I).
C
C ENVIRONMENT: VAX/VMS FORTRAN; MacFORTRAN;
C              IMPLICIT NONE and '!' comments are nonstandard.
C
C HISTORY:
C    07/08/83    DAS      Initial design and code of FSEVAL.
C    04/23/86    RAK      MacFORTRAN version: no unnumbered DOs.
C    06/27/86    RAK      Added X1 and XN for periods other than 2 * pi.
C                         Use IMPLICIT NONE, etc.  Be sure that X1 does
C                         not equal XN - we don't check.
C    03/02/87    RAK      Re-ordered and re-named XMIN, XMAX to match
C                         FSFIT. Permit N = 0 (just ignore B(0) element).
C    07/15/89    DAS      XMIN, XMAX replaced by ORIGIN and PERIOD to
C                         suit FSERIES (which displaces FSFIT now on Mac).
C                         VAX version named FSEVAL2 to avoid reworking
C                         everything using the original FSEVAL (0 : 2 pi).
C
C AUTHOR: David Saunders, Sterling Software/NASA Ames, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   MODE, N, NEVAL
      REAL
     >   A (0:N), B (0:N), FEVAL (NEVAL), ORIGIN, PERIOD, XEVAL (NEVAL)

C     Local constants.

      REAL
     >   ZERO, HALF, TWO, PI
      PARAMETER
     >  (ZERO = 0.0E+0,
     >   HALF = 0.5E+0,
     >   TWO  = 2.0E+0,
     >   PI   = 3.141592653589793E+0)

C     Local variables.

      INTEGER
     >   I, J
      REAL
     >   AJ, AJP1, AJP2, COSX, EVENS, FREQ, ODDS, X1


C     Execution.
C     ----------

      FREQ = TWO * PI / PERIOD
      X1 = ORIGIN
      AJ = ZERO                                   ! For N = 0 case.

      DO 30, I = 1, NEVAL
         COSX = TWO * COS ((XEVAL (I) - X1) * FREQ)
         IF (MODE .EQ. 1) THEN

C           Series is odd.

            EVENS = ZERO
         ELSE
            AJP2 = ZERO
            AJP1 = ZERO
            DO 10, J = N, 1, -1
               AJ   = COSX * AJP1 - AJP2 + A (J)
               AJP2 = AJP1
               AJP1 = AJ
   10       CONTINUE
            EVENS = (A (0) + AJ * COSX) * HALF - AJP2
         END IF

         IF (MODE .EQ. 2) THEN

C           Series is even.

            ODDS = ZERO
         ELSE
            AJP2 = ZERO
            AJP1 = ZERO
            DO 20, J = N, 1, -1
               AJ   = COSX * AJP1  -  AJP2  +  B (J)
               AJP2 = AJP1
               AJP1 = AJ
   20       CONTINUE
            ODDS = AJ * SIN ((XEVAL (I) - X1) * FREQ)
         END IF
         FEVAL (I) = ODDS + EVENS
   30 CONTINUE

      RETURN
      END
