!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine tri_areas (nnode, ntri, xyz, conn, area, area_total)
!
!  Calculate the areas of the triangles of one zone of a 3-space surface
!  triangulation.  Half the magnitude of the cross product of two sides gives
!  the area of each triangle (see earlier subroutine area3.f).
!
!  Also return, for each vertex, the summed areas of triangles sharing that
!  vertex, as needed for area-weighted interpolation of function data from
!  cell centers to cell vertices.
!
!  03/17/2010  D.A.Saunders  Initial implementation.
!
!  Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)  :: nnode              ! # vertices/nodes in the zone
   integer, intent (in)  :: ntri               ! # triangles in the zone
   real,    intent (in)  :: xyz(3,nnode)       ! (x,y,z)s of the vertices
   integer, intent (in)  :: conn(3,ntri)       ! Vertex pointers into xyz(:,*)
   real,    intent (out) :: area(ntri)         ! Triangle areas
   real,    intent (out) :: area_total(nnode)  ! Summed areas surrounding each
                                               ! vertex/node
!  Local constants:

   real, parameter :: half = 0.5, zero = 0.

!  Local variables:

   integer :: i1, i2, i3, itri
   real    :: x1, x2, x3, y1, y2, y3, z1, z2, z3

!  Execution:

   area_total(:) = zero

   do itri = 1, ntri
      i1 = conn(1,itri);  x1 = xyz(1,i1);  y1 = xyz(2,i1);  z1 = xyz(3,i1)
      i2 = conn(2,itri);  x2 = xyz(1,i2);  y2 = xyz(2,i2);  z2 = xyz(3,i2)
      i3 = conn(3,itri);  x3 = xyz(1,i3);  y3 = xyz(2,i3);  z3 = xyz(3,i3)

      area(itri) = half * &
         sqrt ( (y1*(z2 - z3) - z1*(y2 - y3) + y2*z3 - y3*z2) ** 2 +           &
                (z1*(x2 - x3) - x1*(z2 - z3) + z2*x3 - z3*x2) ** 2 +           &
                (x1*(y2 - y3) - y1*(x2 - x3) + x2*y3 - x3*y2) ** 2 )

      area_total(i1) = area_total(i1) + area(itri)
      area_total(i2) = area_total(i2) + area(itri)
      area_total(i3) = area_total(i3) + area(itri)
   end do

   end subroutine tri_areas
