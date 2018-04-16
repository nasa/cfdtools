C+----------------------------------------------------------------------
C
      SUBROUTINE WAGFIT (NX, X, Y, N, MAXWRK, WRK, C, RMSDEV, IER)
C
C  ACRONYM: WAGner series: FITting of first few terms by least squares
C           ---            ---
C  PURPOSE: WAGFIT  fits a linear combination of the  first  N  Wagner
C           functions in the least squares sense to the given dataset,
C           which is probably an airfoil surface or a  thickness  dis-
C           tribution or a meanline (camber) distribution. It may also
C           be a distribution of airfoil perturbations.   It must span
C           the interval [0.,1.], with Y(1) = Y(NX) = 0.
C
C           This version requires that the "ramp" function  needed  to
C           retain nonzero Y at X=1. be dealt with in the calling pro-
C           gram, for reasons of symmetry:  evaluation of the fit must
C           handle this nonzero case,  so  setting up the fit to allow
C           for it should be done in the same place.
C           
C           (WAGFIT originally included the ramp term  as  an  (N+1)th
C           coefficient to compute, but the preferred result is to re-
C           tain any original nonzero Y(NX) exactly.)
C
C  METHOD:  WAGFIT sets up the particular overdetermined system led to
C           by a Wagner series, one column at a time, then uses a gen-
C           eral linear least squares solver  to compute optimal coef-
C           ficients and a goodness-of-fit measure.
C
C  NOTES:   Wagner functions have found application in airfoil design,
C           as combinations of them can produce  airfoil-like  curves.
C           See subroutine BEVAL for details of the  Wagner functions,
C           and for evaluating the linear combination determined here.
C
C  ARGUMENTS:
C    ARG    TYPE  I/O/S   DIM    DESCRIPTION
C    NX      I      I      -     Number of data points being fitted
C    X       R      I      NX    Abscissas of data points, in [0.,1.]
C    Y       R      I      NX    Corresponding ordinates: Y(1)=Y(NX)=0
C    N       I      I      -     Wagner functions defined by  1:N  are
C                                to be fitted as a linear combination,
C                                1 <= N <= NX. Values of N > 20 aren't
C                                permitted  (to  protect  over-zealous
C                                users from themselves).
C    MAXWRK  I      I      -     Maximum work-space provided.  Must be
C                                be at least NX*(N+2).
C    WRK     R      S    MAXWRK  Work-space for least squares solver
C    C       R      O      N     Computed coefs. C(I) applies to Wagner
C                                function I for I = 1:N.
C    RMSDEV  R      O      -     SQRT ( Optimal sum of squares / NX )
C    IER     I      O      -     0 means no errors;
C                                4 means no fit (N < 1 or > 20 );
C                                5 means no fit (N > NX);
C                                6 means no fit (MAXWRK too small).
C  EXTERNAL REFERENCES:
C    BEVAL      Evaluates a Wagner function for given N, X(*)
C    HDESOL     Linear least squares by Householder decomposition
C
C  ENVIRONMENT:  VAX/VMS FORTRAN
C
C  HISTORY:
C    03/05/85   DAS   Initial implementation, including "ramp" term.
C    07/17/86   DAS   Nonzero Y(NX) must be handled outside now.
C    04/05/92   DAS   Meaning of IER=5 above had not been updated.
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C ... Arguments:

      INTEGER    N, NX, MAXWRK, IER
      REAL       C (N), RMSDEV, WRK (MAXWRK), X (NX), Y (NX)

C ... Local declarations:

      INTEGER    I, J, MAXN
      REAL       P (2)
      PARAMETER  (MAXN = 20)

C ... Externals:

      INTRINSIC  SQRT
      EXTERNAL   BEVAL, HDESOL

C ... Execution:

C ... Check for the more obvious erroneous calls:

      IF (N .LT. 1  .OR.  N .GT. MAXN) GO TO 940
      IF (N .GT. NX) GO TO 950
      IF (MAXWRK .LT. NX * (N + 2)) GO TO 960

C ... Set up the over-determined (or square) system by columns:

      P (2) = 1.E+0
      I = 1

      DO 200 J = 1, N
         P (1) = J

         CALL BEVAL ('WAGNER', 2, P, .FALSE., NX, X, WRK (I))

         I = I + NX
  200 CONTINUE

      J = I - 1

C ... Now for the right-hand-side:

      DO 300 I = 1, NX
         WRK (I+J) = Y (I)
  300 CONTINUE

      J = NX + J

C ... Solve the overdetermined system by orthogonal factorization.
C     The required coefs. are returned starting in WRK (J+1):

      CALL HDESOL (NX, NX, N+1, WRK (1), WRK (J+1), RMSDEV)

      DO 400 I = 1, N
         C (I) = WRK (J+I)
  400 CONTINUE

      RMSDEV = SQRT (RMSDEV / NX)

      IER = 0
      GO TO 999

C ... Error handling:

 940  CONTINUE
C ... Requested sum has too many terms to be a reasonable mathematical
C     model fitted by these techniques (i.e., N > MAXN), else N < 1:

      IER = 4
      GO TO 999

 950  CONTINUE
C ... N is too big relative to NX for the mathematical model to be fit
C     this way:

      IER = 5
      GO TO 999

 960  CONTINUE
C     Not enough work-space provided:

      IER = 6

 999  RETURN
      END
