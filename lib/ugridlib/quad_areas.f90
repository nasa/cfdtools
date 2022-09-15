!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine quad_areas (nnode, nquad, xyz, conn, area, area_total)
!
!  Calculate the areas of the cells of one zone of an unstructured quad surface
!  grid.  The cells are expected to be defined in clockwise order via pointers
!  into a list of xyz coordinates.  (Counterclockwise is permitted too.)  Each
!  quad is treated as two triangles by connecting vertices 1 and 3.  A cross
!  product provides the area of the two resulting triangles.
!
!  Also return, for each vertex, the summed areas of quad cells sharing that
!  vertex, as needed for area-weighted interpolation of function data from
!  cell centers to cell vertices.
!
!  Programming note retained for posterity from the original project3 approach:
!
!  03/17/2022  D.A.Saunders  Initial implementation of TET_VOLUMES.
!  04/27/2010    "     "     Analogous adaptation of QUAD_AREAS..
!
!  Author:  David Saunders, AMA, Inc. at NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)  :: nnode               ! # vertices/nodes in the zone
   integer, intent (in)  :: nquad               ! # quad cells in the zone
   real,    intent (in)  :: xyz(3,nnode)        ! (x,y,z)s of the vertices
   integer, intent (in)  :: conn(4,nquad)       ! Vertex pointers into xyz(:,*)
   real,    intent (out) :: area(nquad)         ! Areas of the quad cells
   real,    intent (out) :: area_total(nnode)   ! Summed area surrounding each
                                                ! vertex/node
!  Local constants:

   real, parameter :: zero = 0.

!  Local variables:

   integer :: i1, i2, i3, i4, iquad

!  Execution:

   area_total(:) = zero

   do iquad = 1, nquad
      i1 = conn(1,iquad)
      i2 = conn(2,iquad)
      i3 = conn(3,iquad)
      i4 = conn(4,iquad)

      call quad_area (xyz(:,i1), xyz(:,i2), xyz(:,i3), xyz(:,i4), area(iquad))

      area_total(i1) = area_total(i1) + area(iquad)
      area_total(i2) = area_total(i2) + area(iquad)
      area_total(i3) = area_total(i3) + area(iquad)
      area_total(i4) = area_total(i4) + area(iquad)
   end do

   end subroutine quad_areas
