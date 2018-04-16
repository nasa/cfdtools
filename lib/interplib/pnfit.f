C+----------------------------------------------------------------------
C
      SUBROUTINE PNFIT (MPTS, X, Y, NDEG, THRU00, MAXWRK, WRK,
     +                  COEFS, RMSDEV, IER)
C
C  PURPOSE: PNFIT  fits a polynomial of requested degree in the least
C           squares sense to the given data.   It should be used with
C           PNEVAL  when a polynomial is an appropriate  mathematical
C           model for the data.   This version has an option to force
C           the curve through the origin.
C
C  METHOD:  PNFIT sets up the particular overdetermined system that a
C           polynomial  (as opposed to some other linear model) leads
C           to.   It then uses a general linear least squares  solver
C           to calculate the coefficients of the optimal  polynomial.
C
C  ARGUMENTS:
C    ARG    TYPE  I/O/S   DIM    DESCRIPTION
C    MPTS    I      I      -     Number of data points (>=NDEG+1)
C    X       R      I    MPTS    Abscissas of data points
C    Y       R      I    MPTS    Ordinates of data points
C    NDEG    I      I      -     Requested order of fit (>=0)
C    THRU00  L      I      -     .TRUE. forces curve through (0,0)
C    MAXWRK  I      I      -     Max. work-space provided; should be 
C                                at least MPTS*(NDEG+4)
C    WRK     R      S   MAXWRK   Work-space for HDESOL
C    COEFS   R      O   0:NDEG   Polynomial coefficients in the order
C                                C(0), C(1),...,C(NDEG) expected by
C                                PNEVAL; C(0) = 0. if THRU00 is TRUE
C    RMSDEV  R      O      -     SQRT (Minimum sum of squares / MPTS)
C    IER     I      O      -     0 means no errors;
C                                4 means no fit (NDEG < 0 or > 10 );
C                                5 means no fit (MPTS < NDEG+1);
C                                6 means no fit (MAXWRK too small)
C
C  PROCEDURES:
C    HDESOL     Linear least squares by Householder decomposition.
C
C  NOTES:
C         * Please note that this routine does NOT set up the "Normal
C           Equations" and solve the resulting (square) system, which
C           can be poorly conditioned, particularly where polynomials
C           are involved.   Its general purpose solver works with the
C           overdetermined system directly instead.   This takes more
C           storage (proportional to MPTS), but requires roughly half
C           the machine precision for the same accuracy.   High-order
C           polynomials are generally inadvisable regardless  of  the
C           method used to fit them, so PNFIT sets an arbitrary upper
C           limit on the degree of fit that it will handle,  not  for
C           work-space reasons but to protect the uninformed user. An
C           example may help understand why: suppose a data abscissa,
C           x, is of order 10; this leads to the quantity  x**n for a
C           polynomial of degree n; if n is just 6, the data has been
C           "blown up" to 10**6, even for the present implementation;
C           using the Normal Equations would involve 10**12 in such a
C           case.   Fitting higher-order polynomials involves solving
C           a poorly-conditioned problem because the higher powers of
C           x are very nearly linearly dependent.   Use PNFIT wisely!
C
C         * If tabulations of deviations  are  required,  try  LSPOLY
C           instead.  But PNFIT/PNEVAL are more flexible.
C
C         * Beware of the order of the coefficients.
C
C  ENVIRONMENT:  FORTRAN 77
C
C  HISTORY:
C
C  May 1983  DAS  Adapted from LSPOLY to go with PNEVAL.
C  Apr 1997  DAS  Removed THRU00 test from RHS loop; updated style.
C  Jan 2004  DAS  MAXWRK >= NPTS*(NDEG+3), not NPTS*(NDEG+4).
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C-----------------------------------------------------------------------

C  *  Arguments:

      INTEGER
     >   MPTS, NDEG, MAXWRK, IER
      REAL
     >   X(MPTS), Y(MPTS), WRK(MAXWRK), COEFS(0:NDEG), RMSDEV
      LOGICAL
     >   THRU00

C  *  Local constants:

      INTEGER
     >   MAXDEG
      PARAMETER
     >  (MAXDEG = 10)

C  *  Local variables:

      INTEGER
     >   I, IFREE, N, NCOEFS

C  *  Externals:

      EXTERNAL
     >   HDESOL
      INTRINSIC
     >   SQRT

C  *  Execution:

C  *  Check for the more obvious erroneous calls:

      IF (THRU00) THEN
         NCOEFS = NDEG
         IFREE  = 1
      ELSE
         NCOEFS = NDEG + 1
         IFREE  = MPTS + 1
      END IF

      IF (NCOEFS .LT. 1  .OR.  NDEG .GT. MAXDEG) GO TO 940
      IF (MPTS .LT. NCOEFS) GO TO 950
      IF (MAXWRK .LT. MPTS * (NDEG + 3)) GO TO 960

C  *  Set up the overdetermined (or square) system by columns:

      DO N = 1, NDEG     ! NDEG = 0 is OK
         DO I = 1, MPTS
            WRK(IFREE) = X(I) ** N
            IFREE = IFREE + 1
         END DO
      END DO

      IF (THRU00) THEN
         DO I = 1, MPTS
            WRK(IFREE) = Y(I)
            IFREE = IFREE + 1
         END DO
      ELSE
         DO I = 1, MPTS
            WRK(I) = 1.E+0
            WRK(IFREE) = Y(I)
            IFREE = IFREE + 1
         END DO
      END IF

C  *  Solve the overdetermined system by orthogonal (QR) factorization:

      CALL HDESOL (MPTS, MPTS, NCOEFS+1, WRK(1), WRK(IFREE), RMSDEV)

      RMSDEV = SQRT (RMSDEV / MPTS)

C  *  Desired coefficients are returned in  WRK(IFREE), WRK(IFREE+1),...

      IF (THRU00) THEN
         COEFS(0) = 0.E+0
         N = 1
      ELSE
         N = 0
      END IF

      DO I = 1, NCOEFS
         COEFS(N) = WRK(I-1+IFREE)
         N = N + 1
      END DO

      IER = 0
      GO TO 999

C  *  Error handling:

  940 CONTINUE
C  *  Requested degree of fit not meaningful:
      IER = 4
      GO TO 999

  950 CONTINUE
C  *  Too few data points for this degree of fit:
      IER = 5
      GO TO 999

  960 CONTINUE
C  *  Not enough work-space provided:
      IER = 6

  999 RETURN
      END
