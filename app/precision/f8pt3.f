C+----------------------------------------------------------------------
C
      SUBROUTINE F8PT3 (N, X, F)
C
C  PURPOSE:
C        F8PT3 calculates the 1-D function of Example 8.3 in Practical
C     Optimization, for exercising program PRECISION on.
C
C        Adjust the constant NDIGITS below to suit.
C
C  HISTORY:
C     05/10/91  D.A.Saunders  Second test function for PRECISION.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   N
      REAL
     >   X (N), F

C     Local constants:

      INTEGER
     >   NDIGITS
      PARAMETER
     >  (NDIGITS = 7) ! NDIGITS is the (reduced) number of significant
                      ! digits achieved by internal I/O

C     Local variables:

      CHARACTER
     >   DIGITS * (NDIGITS + 6)

      LOGICAL
     >   FIRST

      DATA FIRST /.TRUE./
         SAVE FIRST

C     Execution:

      IF (FIRST) THEN
         FIRST = .FALSE.
         WRITE (6, '(//, A)') ' F (X) = EXP (X) + X ** 3 - 3X - 1.1',
     >      ' See Example 8.3 of Practical Optimization, p. 336.'
         WRITE (6, '(//, A, I1)')
     >      ' Number of significant digits in the function values: ',
     >      NDIGITS
     >     , ' NO ARTIFICIAL TRUNCATION'
      END IF

      F = EXP (X (1)) + X (1) ** 3  - 3.E+0 * X (1) - 1.1E+0

C***      WRITE (DIGITS, '(E<NDIGITS + 6>.<NDIGITS>)') F
C***      READ  (DIGITS, '(E<NDIGITS + 6>.0)') F

      RETURN
      END
