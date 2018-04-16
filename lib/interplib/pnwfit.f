C+----------------------------------------------------------------------
C
      SUBROUTINE PNWFIT (NPTS, X, Y, WEIGHT, NDEG, MAXWRK, WRK,
     +                   COEFS, RMSDEV, IER)
C
C  ACRONYM: PolyNomial Weighted least squares FIT
C           -   -      -                      ---
C
C  PURPOSE: PNWFIT fits a polynomial of requested degree in the weighted
C           least squares sense to the given data.  (Each equation K of
C           the overdetermined system is weighted by WEIGHT(K).)
C
C  METHOD:  PNWFIT sets up the weighted overdetermined system that a
C           polynomial (as opposed to some other linear model) leads
C           to.  It then uses a general linear least squares solver
C           to factorize the (over)square matrix (QR) and calculate the
C           coefficients of the optimal polynomial.
C
C  NOTES:   * Use PNFIT if no meaningful weighting is known.  PNFIT has
C             the option to force the polynomial through the origin also.
C           * Use PNEVAL or PNDVAL to evaluate the fitted polynomial.
C           * The RMSDEV argument may not be as meaningful here as it is
C             for PNFIT, unless perhaps the weights are of order 1.
C
C  ARGUMENTS:
C    ARG    TYPE  I/O/S   DIM    DESCRIPTION
C    NPTS    I      I      -     Number of data points (>=NDEG+1).
C    X       R      I    NPTS    Abscissas of data points.
C    Y       R      I    NPTS    Ordinates of data points.
C    WEIGHT  R      I    NPTS    Weights for each eqn. of overdet. system.
C    NDEG    I      I      -     Requested order of fit (>=0).
C    MAXWRK  I      I      -     Max. work-space provided. Should be 
C                                at least NPTS*(NDEG+3).
C    WRK     R      S   MAXWRK   Work-space for HDESOL.
C    COEFS   R      O   0:NDEG   Polynomial coefficients in the order
C                                C(0), C(1),...,C(NDEG) expected by PNEVAL.
C    RMSDEV  R      O      -     SQRT( Minimum weighted sum of squares / NPTS ).
C    IER     I      O      -     0 means no errors;
C                                4 means no fit (NDEG < 0 or > 10 );
C                                5 means no fit (NPTS < NDEG+1);
C                                6 means no fit (MAXWRK too small).
C
C  PROCEDURES:
C    HDESOL     Linear least squares by Householder decomposition.
C
C  ENVIRONMENT:  VAX/VMS FORTRAN  (IMPLICIT NONE is the only extension used.)
C
C  HISTORY:
C    Nov. 1987  D.A.Saunders  Adapted from PNFIT for weighted case.
C    Jul. 1989    "     "     Cosmetics, in preparation for Mac translation.
C    Jan. 2004    "     "     MAXWRK >= NPTS*(NDEG+3), not NPTS*(NDEG+4).
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C  *  Arguments:

      INTEGER
     >   NPTS, NDEG, MAXWRK, IER
      REAL
     >   WEIGHT (NPTS), X (NPTS), Y (NPTS), WRK (MAXWRK),
     >   COEFS (0:NDEG), RMSDEV

C  *  Local constants:

      INTEGER
     >   MAXDEG
      PARAMETER
     >  (MAXDEG = 10)

C  *  Local variables:

      INTEGER
     >   I, IFREE, N, NCOEFS

C  *  Execution:

C  *  Check for the more obvious erroneous calls:

      NCOEFS = NDEG + 1
      IF (NCOEFS .LT. 1  .OR.  NDEG .GT. MAXDEG) GO TO 940
      IF (NPTS .LT. NCOEFS) GO TO 950
      IF (MAXWRK .LT. NPTS * (NDEG + 3)) GO TO 960

C  *  Set up the overdetermined (or square) system by columns:

      IFREE = NPTS + 1
      IF (NDEG .LT. 1) GO TO 325

         DO 320 N = 1, NDEG
            DO 310 I = 1, NPTS
               WRK (IFREE) = WEIGHT (I) * X (I) ** N
               IFREE = IFREE + 1
  310       CONTINUE
  320    CONTINUE

  325 DO 330 I = 1, NPTS
         WRK (I) = WEIGHT (I)
         WRK (IFREE) = WEIGHT (I) * Y (I)
         IFREE = IFREE + 1
  330 CONTINUE

C  *  Solve the overdetermined system by Householder (QR) factorization:

      CALL HDESOL (NPTS, NPTS, NCOEFS+1, WRK (1), WRK (IFREE), RMSDEV)

      RMSDEV = SQRT (RMSDEV / NPTS)

C  *  Desired coefficients are returned in  WRK (IFREE), WRK (IFREE+1),...

      N = 0
      DO 340 I = 1, NCOEFS
         COEFS (N) = WRK (I - 1 + IFREE)
         N = N + 1
  340 CONTINUE

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
