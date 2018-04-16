C+----------------------------------------------------------------------
C
      SUBROUTINE FSFIT ( M, X, F, N, MAXWRK, WRK, A, B, RMSDEV, IER )
C
C  ACRONYM: Fourier Series: FITting of first few terms by least squares
C           -       -       ---
C  PURPOSE: FSFIT  fits  a partial Fourier series in the linear least
C           squares sense to the given data.    The series is:
C
C           f(x) = a(0)/2 + SIGMA <n=1:N> a(n)*cos(nx) + b(n)*sin(nx)
C
C           with b(0) = 0 (assigned here).
C
C           Abscissas are assumed to be in radians, but no assumption
C           is made about the their range or uniformity.   Meaningful
C           applications should span an interval of length PI or 2*PI
C           adequately, however.
C
C           Much more efficient methods are available for uniform and
C           periodic data, giving all M of the possible coefficients.
C           An arbitrary limit of N=15 is imposed here,  meaning that
C           in general only the first few of many terms in the series
C           are calculable using FSFIT.
C
C  METHOD:  FSFIT sets up the particular overdetermined system that a
C           partial Fourier series leads to, then uses a general lin-
C           ear least squares solver to calculate the coefficients of
C           the optimal partial series.  No attempt is made to handle
C           odd or even series specially.
C
C  NOTES:   Fitting a(0:N) and b(1:N) in the least squares sense does
C           NOT amount to the same as calculating these  coefficients
C           directly using explicit integrals.   It is up to the user
C           to decide which is best suited to his problem.  If actual
C           harmonics are not of particular interest,  but rather the
C           dataset is required to be smoothed in some way and though
C           noisy, is known to be basically cyclic in nature, then it
C           may be appropriate to use FSFIT.   If the data points are
C           not uniform, FSFIT MAY be an answer for the first handful
C           of coefficients. Plots or tabulations will be facilitated
C           by FSEVAL, which IS efficient.
C
C           Note that calling FSFIT twice with different values of N,
C           for the same input data, will NOT give the same values of
C           coefficients 0 through MIN ( N1, N2 ) in both cases.  The
C           situation is more like fitting sums of powers of X  (that
C           is, polynomials) - the coefficients from polynomials with
C           different degrees for the same data are not related.
C
C  PARAMETERS:
C    ARG    TYPE  I/O/S   DIM       DESCRIPTION
C    M       I      I      -        Number of input data points.
C    X       R      I      M        Input abscissas (radians, but not
C                                   necessarily uniform).
C    F       R      I      M        Data values to be fitted/smoothed
C                                   - presumably basically periodic.
C    N       I      I      -        Defines the partial series.  Must
C                                   have 2*N+1 <= M to give an  over-
C                                   determined problem.   Values of N
C                                   greater than 15 are not permitted
C                                   to protect over-zealous users.
C                                   Must have N >= 0.
C    MAXWRK  I      I      -        Max. work-space provided.  Should
C                                   be at least M*(2*N+3).
C    WRK     R      S    MAXWRK     Work-space for  HDESOL.
C    A       R      O     0:N       Fourier coefs. a(0), a(1),..,a(N)
C    B       R      O     0:N       Fourier coefs. b(0), b(1),..,b(N)
C                                   where b(0) = 0 always.
C    RMSDEV  R      O      -        SQRT(Optimal sum of squares / M).
C    IER     I      O      -        0 means no errors;
C                                   4 means no fit (N < 1 or > 20 );
C                                   5 means no fit (M < 2*N+1);
C                                   6 means no fit (MAXWRK too small).
C  EXTERNAL REFERENCES:
C    HDESOL     Linear least squares by Householder decomposition.
C
C  ENVIRONMENT:  VAX/VMS FORTRAN
C
C  DEVELOPMENT HISTORY:
C    DATE   INITIALS   DESCRIPTION
C  06/20/83   DAS      Initial design and code.
C  12/23/85   DAS      Revised documentation - more qualifications.
C  03/03/87   DAS/RAK  N=0 should be a valid input - it is now.
C  05/18/04   DAS      Raised the limit from 15 to 20 terms.
C
C  AUTHOR:  David Saunders, Informatics General, Palo Alto, CA.
C
C-----------------------------------------------------------------------

C  *  Arguments:

      INTEGER           M, N, MAXWRK, IER
      REAL              X(M), F(M), WRK(MAXWRK), A(0:N), B(0:N), RMSDEV

C  *  Local constants:

      INTEGER           MAXN
      REAL              PI
      PARAMETER         ( MAXN = 20, PI = 3.14159265358979323 )

C  *  Local variables:

      INTEGER           I, IRHS, ISOLN, J, L

C  *  Externals:

      EXTERNAL          HDESOL
      INTRINSIC         COS, SIN, SQRT

C  *  Execution:

C  *  Check for the more obvious erroneous calls:

      IF ( N .LT. 0  .OR.  N .GT. MAXN ) GO TO 940
      IF ( M .LT. 2*N + 1 ) GO TO 950
      IF ( MAXWRK .LT. M * ( 2*N+3 ) ) GO TO 960

C  *  Set up the overdetermined (or square) system by rows:

      DO 250 I = 1, M
         L = M*N + I
         WRK(I) = 0.5E+0

         DO 200 J = 1, N
            WRK(M*J+I) = COS( J*X(I) )
            WRK(M*J+L) = SIN( J*X(I) )
  200    CONTINUE

C  *     The RHS goes in as a final column of the overdetermined system:

         IRHS = M*(2*N+1) + I
         WRK(IRHS) = F(I)
  250 CONTINUE

      ISOLN = IRHS+1

C  *  Solve the overdetermined system by Householder factorization.
C     Fourier coefs a0,a1,...,b1,b2,... are returned starting in WRK(ISOLN):

      CALL HDESOL ( M, M, 2*N+2, WRK(1), WRK(ISOLN), RMSDEV )

      A(0) = WRK(ISOLN)
      B(0) = 0.E+0

      DO 300 I = 1, N
         A(I) = WRK(I+ISOLN)
         B(I) = WRK(I+ISOLN+N)
  300 CONTINUE

      RMSDEV = SQRT( RMSDEV / M )

      IER = 0
      GO TO 999

C  *  Error handling:

 940  CONTINUE
C  *  Requested partial sum has too many terms to be a reasonable mathematical
C     model fitted by these techniques ( i.e., N > MAXN ), else N < 0:
      IER = 4
      GO TO 999

 950  CONTINUE
C  *  N is too big relative to M for the mathematical model to be fitted
C     this way:
      IER = 5
      GO TO 999

 960  CONTINUE
C  *  Not enough work-space provided:
      IER = 6

 999  RETURN
      END
