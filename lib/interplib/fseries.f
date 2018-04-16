C+----------------------------------------------------------------------
C
      SUBROUTINE FSERIES (MODE, NX, X, F, G, N, A, B, IER)
C
C ONE-LINER: Fourier series coefs. for irregular data, full or half period
C       
C PURPOSE:
C     FSERIES computes the coefficients A(J), B(J), J=0:N representing
C     the finite Fourier series for an arbitrary dataset, which may be
C     indicated as spanning either half the underlying cycle (MODE=1,
C     meaning sine terms only) or the full cycle (MODE=3, meaning both
C     sines and cosines).  (Cosines only does not apply here; these
C     choices for MODE match those of FSEVAL/FSEVAL2.)
C
C     FSERIES offers two kinds of "fitting" for irregular datasets, with
C     modest values of NX and N expected since much more efficient methods
C     are available for regular data.
C
C METHOD:
C     The composite trapezoidal rule is used for the integrations, since
C     irregular abscissas prevent use of efficient methods.  Periodic
C     continuation is achieved as well as possible by subtracting the line
C     connecting the first and last points - the coefficients actually apply
C     to these transformed values, G(I):
C                           N             _               _
C        G(x) = a(0)/2  +  SUM  a(j) cos jx  +  b(j) sin jx
C           _              j=1
C     where x is x transformed to [0, pi] or [0, 2 pi] depending on MODE.
C
C     Subroutine FSSUM uses FSEVAL2 to evaluate this representation efficiently
C     at given abscissas, then add back the "ramp" term.
C     
C ARGUMENTS:
C    ARG    DIM  TYPE I/O/S DESCRIPTION
C    MODE    -     I    I   MODE = 1 means treat data as spanning half the
C                                    underlying cycle:  sine terms only;
C                           MODE = 3 means full-period data: sines + cosines.
C     NX     -     I    I   The number of data points.  NX >= 3.
C     X      NX    R    I   Data abscissas, presumed irregular but monotonic.
C     F      NX    R    I   Data ordinates corresponding to X, not necessarily
C                           periodic.
C     G      NX    R    S   Work-space for the transformed ordinates.
C     N      -     I    I   Fourier coefficient pairs 0:N will be calculated.
C                           N <= (NX - 1) / 2.
C     A     0:N    R    O   Cosine coefficients (all zero if MODE = 1).
C     B     0:N    R    O   Sine coefficients. (B(0) is always zero.)
C    IER     -     I    O   IER = 0 means no input errors were noticed;
C                           IER = 1 means NX is too small: NX < 3.
C                           IER = 2 means N is too large: N > (NX - 1)/2.
C
C EXTERNAL REFERENCES: None
C
C ENVIRONMENT: VAX/VMS FORTRAN; MacFORTRAN;
C              IMPLICIT NONE is nonstandard.
C
C HISTORY:
C  07/15/89   D.A.Saunders   Adapted from earlier ideas, mainly for the
C                            Mac to provide two basic capabilities compactly.
C
C AUTHORS: David Saunders, Sterling Software/NASA Ames, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   MODE, NX, N, IER
      REAL
     >   X (NX), F (NX), G (NX), A (0:N), B (0:N)

C     Local constants:

      REAL
     >   ONE, ZERO, PI
      PARAMETER
     >  (ONE  = 1.E+0,
     >   ZERO = 0.E+0,
     >   PI   = 3.141592653589793E+0)

C     Local variables:

      INTEGER
     >   I, J
      REAL
     >   ENDS, FREQ, F1, RRANGE, SLOPE, SUM, X1

C     Execution:

      IER = 1
      IF (NX .LE. 2) GO TO 999
      IER = 2
      IF (N .LT. 0 .OR. N .GT. (NX - 1) / 2) GO TO 999
      IER = 0

      F1 = F (1)
      X1 = X (1)
      RRANGE = ONE / (X (NX) - X1)
      SLOPE = (F (NX) - F1) * RRANGE

C     Subtract the "ramp" term to avoid discontinuities:

      DO 100, I = 1, NX
         G (I) = F (I) - F1 - (X (I) - X1) * SLOPE
  100 CONTINUE

C     Apply trapezoidal rule to integral for each coefficient:

      FREQ = PI * RRANGE
      IF (MODE .NE. 1) THEN
         FREQ = FREQ + FREQ
         ENDS = (X (2) - X1) * G (1) + (X (NX) - X (NX-1)) * G (NX)

         DO 300, J = 0, N
            SUM = ZERO
            DO 200, I = 2, NX-1
               SUM = (X (I+1) - X (I-1)) * G(I) *
     >               COS (J * (X (I) - X1) * FREQ) + SUM
  200       CONTINUE
            A (J) = (SUM + ENDS) * RRANGE
  300    CONTINUE

      ELSE
         DO 400, J = 0, N
            A (J) = ZERO
  400    CONTINUE
      END IF

      B (0) = ZERO
      DO 600, J = 1, N
         SUM = ZERO
         DO 500, I = 2, NX-1
            SUM = (X (I+1) - X (I-1)) * G(I) *
     >            SIN (J * (X (I) - X1) * FREQ) + SUM
  500    CONTINUE
         B (J) = SUM * RRANGE
  600 CONTINUE

  999 RETURN
      END
