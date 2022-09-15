!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine hex_volume (p1, p2, p3, p4, p5, p6, p7, p8, volume)

!  Purpose:
!
!     Calculate the volume of one hexahedral cell by treating it as three
!     pyramids with quadrilateral bases.  This is in contrast to the method of
!     Antony Jameson, which treats six such pyramids meeting at the cell center.
!     The hex vertices should follow the Fluent convention: 1234 in circular
!     order on one face and 5678 on the opposite face.
!
!  History:
!
!     03/12/2022  D.A.Saunders  Initial implementation as needed for moving
!                               cell-centered function data to the vertices.
!     03/16/2022    "     "     The Newton iteration in pyramid_volume may not
!                               converge normally and is surely more expensive
!                               than the concise formulation of tet_vol as
!                               used 6 times by hex_vol, q.v.
!
!  Author:  David Saunders, AMA, Inc. at NASA Ames Research Center, CA.
!     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   real, intent (in), dimension (3) :: p1, p2, p3, p4, &  ! Vertex (x,y,z)s
                                       p5, p6, p7, p8     ! in Fluent order
   real, intent (out) :: volume                           ! Cell volume

!  Local variables:

   real :: vol1, vol2, vol3

!  Procedures:

   external :: pyramid_volume

!  Execution:

   call pyramid_volume (p1, p2, p3, p4, p8, vol1)
   call pyramid_volume (p1, p2, p5, p6, p8, vol2)
   call pyramid_volume (p2, p3, p7, p6, p8, vol3)

   volume = vol1 + vol2 + vol3

   end subroutine hex_volume
