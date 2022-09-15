!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine quad_area (xyz1, xyz2, xyz3, xyz4, area)
!
!  Purpose:
!     Calculate the area of the quad cell defined by 3-space points 1, 2, 3, 4
!     in circular order, as the sum of two triangle areas formed by connecting
!     points 1 and 3.  Each triangle area is provided by half the magnitude of
!     a cross product.
!
!  History:
!     04/08/2022  D.A.Saunders  Preferable to the area4 function from 1989.
!
!  Author:  David Saunders, AMA, Inc. at NASA Ames Research Center.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   real, intent (in), dimension (3) :: xyz1, xyz2, xyz3, xyz4
   real, intent (out) :: area

!  Local constants:

   real, parameter :: half = 0.5

!  Local variables:

   real :: A1, A2
   real, dimension (3) :: v1, v2, cp

!  Execution:

   v1(:) = xyz1(:) - xyz2(:)
   v2(:) = xyz3(:) - xyz2(:)
   call cross_product (v1, v2, cp)
   A1 = sqrt (dot_product (cp, cp))

   v1(:) = xyz3(:) - xyz4(:)
   v2(:) = xyz1(:) - xyz4(:)
   call cross_product (v1, v2, cp)
   A2 = sqrt (dot_product (cp, cp))

   area = half*(A1 + A2)

   end subroutine quad_area
