!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine angle2d (n, i, x, y, angle)
!
!  For interior point i on a 2-space curve defined by x(i), y(i), i = 1, n,
!  compute the angle between vectors defined by points i, i+/-1.
!  One potential use is to detect oscillations in the curve via small angles.
!
!  History:
!     06/26/2014  D.A.Saunders  Initial implementation, for FIAT_Opt.
!
!  Author:  David Saunders, ERC, Inc. at NASA Ames Research Center, CA
!     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)  :: n           ! Number of data points
   integer, intent (in)  :: i           ! Interior point of interest
   real,    intent (in)  :: x(n), y(n)  ! Data point coordinates
   real,    intent (out) :: angle       ! Desired angle in [0, 180] degrees

!  Local variables:

   real :: v1(2), v2(2)

!  Execution:

   v1(1) = x(i-1) - x(i)
   v2(1) = x(i+1) - x(i)
   v1(2) = y(i-1) - y(i)
   v2(2) = y(i+1) - y(i)

   angle = acosd (      dot_product (v1, v2) / &
                  sqrt (dot_product (v1, v1) * dot_product (v2, v2)))

   end subroutine angle2d
