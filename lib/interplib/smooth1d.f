C+----------------------------------------------------------------------
C
      SUBROUTINE SMOOTH1D (I1, I2, S, V)
C
C     SMOOTH1D smooths a vector in-place by an explicit method with
C     hard-coded parameters.  End elements are unchanged.
C
C     06/23/97  D.Saunders  GRIDGEN-type explicit smoothing to cope with
C                           spikes in edge Phi, Psi.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN)    :: I1, I2   ! Index range to be smoothed
      REAL,    INTENT (IN)    :: S(I1:I2) ! Arc-lengths associated with V(I1:I2)
      REAL,    INTENT (INOUT) :: V(I1:I2) ! Vector being smoothed; elements
                                          ! I1 & I2 are left unchanged
C     Local constants:

      INTEGER, PARAMETER :: NITERS = 5
      REAL,    PARAMETER :: ALPHA = 0.5, ONE = 1.

C     Local variables:

      INTEGER  I, ITER
      REAL     AI, BI, CI, HL, HR, TERM
      REAL     W(I1:I2) ! Work-space for previous iterate (to vectorize)

C     Execution:

      DO ITER = 1, NITERS

         W(I1:I2) = V(I1:I2)

         DO I = I1 + 1, I2 - 1
            HL   = S(I) - S(I-1)
            HR   = S(I+1) - S(I)
            TERM = ALPHA * ((MIN (HL, HR) ** 2) / (HL + HR))
            AI   = TERM / HL
            CI   = TERM / HR
            BI   = ONE - (AI + CI)
            V(I) = AI * W(I-1) + BI * W(I) + CI * W(I+1)
         END DO

      END DO

      END SUBROUTINE SMOOTH1D
