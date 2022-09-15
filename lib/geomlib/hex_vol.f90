!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine hex_vol (p1, p2, p3, p4, p5, p6, p7, p8, vol)

!  Purpose:
!
!     Calculate the volume of one hexahedral cell by treating it as six
!     tetrahedra.  This is in contrast to the method of Antony Jameson, which
!     treats six quadrilateral-base pyramids meeting at the cell center.
!     The hex vertices should follow the Fluent convention: 1234 in circular
!     order on one face and 5678 on the opposite face.
!
!     Consider a cubic cell:  dividing it into two equal halves with a diagonal
!     slice through one face to the opposite face helps visualize two sets of
!     three tetrahedra.  The order of the vertices passed to tet_vol seems to
!     matter, and the abs has been omitted from the tet_vol result in order for
!     hex_vol to provide a second opinion in the GSMOOTH application, where
!     negative volumes are looked for as a sign of bad cells.
!
!  History:
!
!     03/12/2022  D.A.Saunders  Implementation of hex_volume as three quadri-
!                               lateral-base pyramids, as needed for moving
!                               cell-centered function data to the vertices.
!     03/15/2022    "     "     Hex_vol adapted from hex_volume to make use of
!                               tet_vol (6 times) rather than pyramid_vol, which
!                               has to compute a perpendicular distance via a
!                               Newton iteration that may not converge normally.
!     03/16/2022    "     "     The fourth tet vertex order p7, p8, p5, p3
!                               produces a negative volume (unclear why).  The
!                               order of all vertices below is empirically not
!                               in need of an abs on the tet_vol results when
!                               the cell is a perfect cube, so the abs has been
!                               omitted in tet_vol for the reason stated above.
!
!  Author:  David Saunders, AMA, Inc. at NASA Ames Research Center, CA.
!     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   real, intent (in), dimension (3) :: p1, p2, p3, p4, &  ! Vertex (x,y,z)s
                                       p5, p6, p7, p8     ! in Fluent order
   real, intent (out) :: vol                              ! Cell volume

!  Local variables:

   real :: vol1, vol2, vol3, vol4, vol5, vol6

!  Procedures:

   external :: tet_vol

!  Execution:

!  Consider two halves formed by a diagonal slice if it were a cube:

   call tet_vol (p1, p2, p3, p5, vol1)
   call tet_vol (p2, p3, p5, p7, vol2)
   call tet_vol (p2, p5, p6, p7, vol3)

   call tet_vol (p7, p5, p8, p3, vol4)
   call tet_vol (p8, p5, p1, p7, vol5)
   call tet_vol (p8, p1, p4, p3, vol6)

   vol = vol1 + vol2 + vol3 + vol4 + vol5 + vol6

!! write (*, '(a, 7es16.8)') &
!!    'vol123456tot:', vol1, vol2, vol3, vol4, vol5, vol6, vol

   end subroutine hex_vol
