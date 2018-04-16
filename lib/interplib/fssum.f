C+----------------------------------------------------------------------
C
      SUBROUTINE FSSUM (MODE, NEVAL, XEVAL, NX, X, F, N, A, B, FEVAL)
C
C ONE-LINER:  Fourier series summation (irregular data fit evaluation).
C
C PURPOSE:
C     FSSUM evaluates the Fourier series-type function determined by
C     FSERIES, at the given abscissas.
C
C METHOD:
C     FSEVAL2 does the bulk of the work efficiently, but FSSUM is needed
C     to hide the addition of the "ramp" term subtracted by FSERIES prior
C     to its calculation of the Fourier coefficients.
C
C ARGUMENTS:
C    ARG    DIM  TYPE I/O/S DESCRIPTION
C    MODE    -     I    I   MODE = 1 means treat data as spanning half the
C                                    underlying cycle:  sine terms only;
C                           MODE = 3 means full-period data: sines + cosines.
C   NEVAL    -     I    I   Number of abscissas at which the function is
C			    to be evaluated.  NEVAL >= 1.
C   XEVAL  NEVAL   R    I   Abscissas in the original coordinate system (not
C			    necessarily radians;  no ordering required, but
C                           should be in the original range.
C     NX     -     I    I   The number of data points fitted by FSERIES.
C    X,F     NX    R    I   Abscissas and ordinates of the original data.
C                           Only the end-points are used here.
C     N      -     I    I   Fourier coefficient pairs 0:N will be applied.
C     A     0:N    R    I   Cosine coefficients (all zero if MODE = 1).
C     B     0:N    R    I   Sine coefficients. (B(0) is always zero.)
C   FEVAL  NEVAL   R    O   Desired interpolations.
C
C NOTE:  XEVAL may be passed as X (to smooth the original data), and
C        FEVAL may be passed as F (to overwrite the original data).
C
C EXTERNAL REFERENCES:
C   FSEVAL2     Sums sine and cosine series efficiently (arbitrary coords.)
C
C ENVIRONMENT: VAX/VMS FORTRAN; MacFORTRAN;
C              IMPLICIT NONE is nonstandard.
C
C HISTORY:
C  07/15/89   D.A.Saunders   Adapted earlier ideas as needed for FSERIES.
C
C AUTHORS: David Saunders, Sterling Software/NASA Ames, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   MODE, NEVAL, NX, N
      REAL
     >   XEVAL (NEVAL), FEVAL (NEVAL), X (NX), F (NX), A (0:N), B (0:N)

C     Local variables:

      INTEGER
     >   I
      REAL
     >   F1, PERIOD, SLOPE, X1

C     Execution:

      F1 = F (1)
      X1 = X (1)
      PERIOD = X (NX) - X1
      SLOPE = (F (NX) - F1) / PERIOD
      IF (MODE .EQ. 1 ) PERIOD = PERIOD + PERIOD

C     Evaluate the Fourier series represented by the a(j)s and b(j)s:

      CALL FSEVAL2 (MODE, NEVAL, XEVAL, X1, PERIOD, N, A, B, FEVAL)

C     Restore the "ramp" subtracted by FSERIES:

      DO 100, I = 1, NEVAL
         FEVAL (I) = FEVAL(I) + F1 + (XEVAL (I) - X1) * SLOPE
  100 CONTINUE

      RETURN
      END
