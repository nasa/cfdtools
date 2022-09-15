!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine tet_vol (p1, p2, p3, p4, vol)

!  Calculate the volume of a tetrahedron by the most concise of formulas, found
!  in Wolfram, for example:
!
!     6 V = | a . (b x c) |  where a, b, c are vectors from vertex 1, 2, 3 to
!                            vertex 4.
!
!  03/14/2022  D.A.Saunders  ! Should be preferable to calculating vertical
!                            ! heights of pyramids with quadrilateral bases
!                            ! originally used for hexahedral cell volumes.
!  03/16/2022    "     "     ! Puzzling test results were explained by spurious
!                            ! squaring of the dot product above. Thanks, Chris!
!
!  Author:  David Saunders, AMA, Inc. at MASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   real, intent (in), dimension (3) :: p1, p2, p3, p4  ! (x,y,z)s of vertices
   real, intent (out) :: vol                           ! Tetrahedron volume

!  Local constants:

   real, parameter :: sixth = 1./6.

!  Local variables:

   real, dimension (3) :: v1, v2, v3, v4

!  Execution:

   v1(:) = p4(:) - p1(:)
   v2(:) = p4(:) - p2(:)
   v3(:) = p4(:) - p3(:)

   call cross_product (v2, v3, v4)

   vol = sixth * dot_product (v1, v4)

   end subroutine tet_vol
