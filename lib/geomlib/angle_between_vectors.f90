!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine angle_between_vectors (v1, v2, angle)
!
!  Calculate the angle in degrees between two vectors in 3-space.
!
!  History:
!     02/15/2018  D.A.Saunders  Initial implementation.
!     05/05/2018    "     "     The arc cosine needs to be protected against
!                               arguments slightly outside +/-1.
!
!  Author:  David Saunders, AMA, Inc. at NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use trigd
   implicit none

!  Arguments:

   real,    intent (in)  :: v1(3), v2(3)  ! Given vectors
   real,    intent (out) :: angle         ! Desired angle in [0, 180] degrees

!  Local variables:

   real :: cosine

!  Execution:

   cosine = max (-1., min (1., dot_product (v1, v2) / &
                         sqrt (dot_product (v1, v1) * dot_product (v2, v2))))
   angle  = acosd (cosine)

   end subroutine angle_between_vectors
