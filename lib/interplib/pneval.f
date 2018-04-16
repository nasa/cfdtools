C+----------------------------------------------------------------------
C
      SUBROUTINE PNEVAL ( NDEG, COEFS, NPTS, X, Y )
C
C ACRONYM:  PolyNomial EVALuation
C           -   -      ----
C PURPOSE:  PNEVAL evaluates the polynomial curve defined by NDEG and
C           COEFS(*), at the abscissas supplied in X(*).  (Use PNDVAL
C           instead if you need 1st or 2nd derivatives as well.)
C
C METHOD:   Horner's rule (valid for NDEG>=0), for each abscissa.
C     
C PARAMETERS:
C    ARG    DIM      TYPE I/O/S DESCRIPTION 
C   NDEG     -        I     I   Degree of polynomial (>=0)
C   COEFS  0:NDEG     R     I   Coefficients in the order
C                               C(0), C(1), ..., C(NDEG)
C   NPTS     -        I     I   No. of points to evaluate curve at (>=1)
C   X       NPTS      R     I   Abscissas
C   Y       NPTS      R     O   Ordinates on the curve
C    
C AUTHOR: David Saunders, Informatics, May 1983.
C         DAS, 06/12/86:  Double precision accumulation of sums added -
C                         it helps extreme cases (e.g. X=year number)
C                         and does not slow normal cases signficantly.
C         DAS, 12/11/07:  64-bit version shouldn't have SNGL in it.
C                         Compile with -r8 or equivalent.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER           NDEG, NPTS
      REAL              COEFS(0:NDEG), X(NPTS), Y(NPTS)

      INTEGER           I, J
      REAL              SUM, XI

      DO I = 1, NPTS
         SUM = 0.
         XI  = X(I)
         DO J = NDEG, 0, -1
            SUM = SUM*XI + COEFS(J)
         END DO
         Y(I) = SUM
      END DO

      END
