!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE SMOOTH2D (NITERS, NI, NJ, U, V, F)
!
!     SMOOTH2D smooths a function of 2 variables defined on a structured grid,
!     in-place, using an iterative explicit method with hard-coded parameters.
!     Interior edge values are smoothed by a 1-dimensional method (as opposed
!     to an edge version of the 2-D method) to allow application to grid blocks
!     that are being processed independently but must produce identical results
!     along common boundaries.  This edge smoothing may be suppressed by passing
!     NITERS < 0.  Corner values are never changed.
!
!     02/15/05  D.Saunders  Initial analog of SMOOTH1D for smoothing the arc-
!                           length perturbations calculated by linear methods
!                           of tailoring outer boundaries of hypersonic flow
!                           grids to the shock iso-surface.
!     02/18/05    "    "    The arc lengths need not be normalized, and indeed
!                           may be better left unnormalized, depending on the
!                           application; NITERS < 0 suppresses edge smoothing.
!     07/05/05    "    "    Edges NI and NJ were being smoothed after the
!                           interior smoothing.  Better to smooth all edges, if
!                           any, before smoothing the interior each iteration.
!     07/06/05    "    "    The previous change made negligible difference for
!                           20 iterations.  However, assuming the surface grid
!                           is far from square, and u and v are unnormalized,
!                           weighting the combination of i and j direction
!                           smoothings should help avoid skewing at the
!                           boundaries that simple averaging introduces.
!                           Parabolic decay from 1 at mid arc length to 0
!                           towards 0 at the edges helps extreme cases.
!     10/13/06    "    "    The parabolic weighting was right for normalized
!                           arc lengths only!
!
!     David Saunders, ELORET/NASA Ames Research Center, Moffett Field, CA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IMPLICIT NONE

!     Arguments:

      INTEGER, INTENT (IN)    :: NITERS   ! |NITERS| = # iterations;
                                          ! NITERS < 0 leaves all edges alone
      INTEGER, INTENT (IN)    :: NI, NJ   ! Surface array dimensions
      REAL,    INTENT (IN)    :: U(NI,NJ) ! Arc lengths in the I and J dirns.,
      REAL,    INTENT (IN)    :: V(NI,NJ) !    normalized or not; e.g., PARAM2D
      REAL,    INTENT (INOUT) :: F(NI,NJ) ! Surface function being smoothed;
                                          ! corner elements are left unchanged
!     Local constants:

      REAL, PARAMETER :: ALPHA = 0.5, BETA = 0.5, FOUR = 4., ONE = 1.

!     Local variables:

      INTEGER :: I, ITER, J
      LOGICAL :: SMOOTH_EDGES
      REAL    :: AI, BI, CI, AJ, BJ, CJ, DL, DR, FI, FJ, HL, HR
      REAL    :: P, Q, TERMI, TERMJ, WI, WJ
      REAL, ALLOCATABLE :: G(:,:)

!     Execution:

      SMOOTH_EDGES = NITERS > 0

      ALLOCATE (G(NI,NJ))

      DO ITER = 1, ABS (NITERS)

         G = F

         IF (SMOOTH_EDGES) THEN

!           I = 1 edge:

            DO J = 2, NJ - 1
               DL     = V(1,J) - V(1,J-1)
               DR     = V(1,J+1) - V(1,J)
               TERMJ  = BETA * ((MIN (DL, DR) ** 2) / (DL + DR))
               AJ     = TERMJ / DL
               CJ     = TERMJ / DR
               BJ     = ONE - (AJ + CJ)
               G(1,J) = AJ * F(1,J-1) + BJ * F(1,J) + CJ * F(1,J+1)
            END DO

            I = NI ! edge:

            DO J = 2, NJ - 1
               DL     = V(I,J) - V(I,J-1)
               DR     = V(I,J+1) - V(I,J)
               TERMJ  = BETA * ((MIN (DL, DR) ** 2) / (DL + DR))
               AJ     = TERMJ / DL
               CJ     = TERMJ / DR
               BJ     = ONE - (AJ + CJ)
               G(I,J) = AJ * F(I,J-1) + BJ * F(I,J) + CJ * F(I,J+1)
            END DO

!           J = 1 edge:

            DO I = 2, NI - 1
               HL     = U(I,1) - U(I-1,1)
               HR     = U(I+1,1) - U(I,1)
               TERMI  = ALPHA * ((MIN (HL, HR) ** 2) / (HL + HR))
               AI     = TERMI / HL
               CI     = TERMI / HR
               BI     = ONE - (AI + CI)
               G(I,1) = AI * F(I-1,1) + BI * F(I,1) + CI * F(I+1,1)
            END DO

            J = NJ ! edge:

            DO I = 2, NI - 1
               HL     = U(I,J) - U(I-1,J)
               HR     = U(I+1,J) - U(I,J)
               TERMI  = ALPHA * ((MIN (HL, HR) ** 2) / (HL + HR))
               AI     = TERMI / HL
               CI     = TERMI / HR
               BI     = ONE - (AI + CI)
               G(I,J) = AI * F(I-1,J) + BI * F(I,J) + CI * F(I+1,J)
            END DO    

         END IF
 
!        Interior points:

         DO J = 2, NJ - 1
            DO I = 2, NI - 1
               P      = U(I,J) / U(NI,J)
               HL     = U(I,J) - U(I-1,J)
               HR     = U(I+1,J) - U(I,J)
               WI     = FOUR * P * (ONE - P) ! Parabolic variation between 0 & 1
               TERMI  = ALPHA * ((MIN (HL, HR) ** 2) / (HL + HR))
               AI     = TERMI / HL
               CI     = TERMI / HR
               BI     = ONE - (AI + CI)
               FI     = AI * F(I-1,J) + BI * F(I,J) + CI * F(I+1,J)

               Q      = V(I,J) / V(I,NJ)
               DL     = V(I,J) - V(I,J-1)
               DR     = V(I,J+1) - V(I,J)
               WJ     = FOUR * Q * (ONE - Q) ! Parabolic variation between 0 & 1
               TERMJ  = BETA * ((MIN (DL, DR) ** 2) / (DL + DR))
               AJ     = TERMJ / DL
               CJ     = TERMJ / DR
               BJ     = ONE - (AJ + CJ)
               FJ     = AJ * F(I,J-1) + BJ * F(I,J) + CJ * F(I,J+1)

               G(I,J) = (WI * FI + WJ * FJ) / (WI + WJ)
            END DO
         END DO

         F = G

      END DO ! Next iteration

      DEALLOCATE (G)

      END SUBROUTINE SMOOTH2D
