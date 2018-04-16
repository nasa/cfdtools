!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE CHANGEN1D (I1, I2, T1, IA, IB, T2)
!
!        CHANGEN1D redistributes abscissas T1(I1:i2) to a different number
!     of abscissas T2(IA:IB) with the same end points and the same relative
!     spacing, meaning if the number of spaces is doubled, the spaces are all
!     halved (say). The abscissas might be "x" or "t" (times in a time history).
!
!        The formulation comes from the linear transformation that converts
!     x in [a,b] to y in [p,q], namely y = ((q-p)/(b-a))x + (bp-aq)/(b-a).
!     Here, we want to transform j in [ia,ib] to i in [i1,i2], giving a
!     fractional index or a data interval in which to interpolate linearly for
!     each output index j (not as obvious as one would think).
!
!     12/22/97  DAS  3-space grid line utility CHANGEN wrapped around PLSCRV3D.
!     10/23/02   "   2-space variant: CHANGEN2D.
!     06/23/14   "   1-space variant: CHANGEN1D, for a time history application.
!                    1:m --> 1:n is the likely usage, but need not be.
!
!     Author:  David Saunders, ERC, Inc. at NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IMPLICIT NONE

!     Arguments:

      INTEGER, INTENT (IN)  :: I1, I2  ! First & last data points
      REAL,    INTENT (IN)  :: T1(I2)  ! Data points in I1:I2
      INTEGER, INTENT (IN)  :: IA, IB  ! First & last redistributed points
      REAL,    INTENT (OUT) :: T2(IB)  ! Redistributed coordinates in IA:IB

!     Local constants:

      REAL, PARAMETER :: ONE = 1.

!     Local variables:

      INTEGER :: I, IP
      REAL    :: P, R, RI1, RIP

!     Execution:

      R   = REAL (I2 - I1) / REAL (IB - IA)
      RI1 = REAL (I1)

      DO I = IA, IB - 1
         RIP   = RI1 + R * REAL (I - IA)
         IP    = INT (RIP)
         P     = RIP - REAL (IP)
         T2(I) = (ONE - P) * T1(IP) + P * T1(IP+1)
      END DO

      T2(IB) = T1(I2)

      END SUBROUTINE CHANGEN1D
