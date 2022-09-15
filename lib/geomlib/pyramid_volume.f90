!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine pyramid_volume (p1, p2, p3, p4, p5, volume)

!  Purpose:
!
!     Calculate the volume of one pyramid with a quadrilateral base as:
!
!                     V = (base area * height) / 3
!
!     The five pyramid vertices should start with the base points in circular
!     order.  The height of the 5th pt. is calculated by orthogonal projection.
!
!  History:
!
!     03/12/2022  D.A.Saunders  Part of computing a hex cell volume.
!
!  Author:  David Saunders, AMA, Inc. at NASA Ames Research Center, CA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   real, intent (in), dimension (3) :: p1, p2, p3, p4, p5  ! Vertex (x,y,z)s
   real, intent (out) ::               volume              ! Pyramid volume

!  Local constants:

   real, parameter :: one_third = 1./3.

!  Local variables:

   integer :: ier
   real    :: area, heightsq, p, q, xfoot(3)

!  Procedures:

   real, external :: area4  ! Area of a quadrilateral in 3-space
   external :: project4     ! Distance ofrom a point to a quadrilateral

!  Execution:

   area = area4 (p1(1), p1(2), p1(3), p2(1), p2(2), p2(3), &
                 p3(1), p3(2), p3(3), p4(1), p4(2), p4(3))

   call project4 (p1, p2, p3, p4, p5, xfoot, heightsq, p, q, ier)

!!!   if (ier /= 0) write (*, '(a, i2)') &
!!!      'Pyramid_volume: Project4 ier =', ier  ! Keep going

   volume = one_third * area * sqrt (heightsq)

   end subroutine pyramid_volume
