!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine tri_normal_and_area (nnode, ntri, tri_xyz, conn, itri, area,     &
                                   unit_normal)

!  Calculate the area and unit normal of one element of a surface triangulation.
!  Both outputs require the same cross product vector, so it makes sense to
!  calculate both together.  Collapsed triangles should be trapped at a higher
!  level.
!
!  Note that the projected area defined by a unit look direction vector d is
!  |d.n| A, where A is the full area, and if d.n is negative for outward
!  pointing normal n, then the triangle can be seen from that direction.
!  Otherwise it is hidden if the surface containing the triangle is a closed
!  convex surface.
!
!  Consistent with companion routine "triangulate_patches", vertices 1, 2, 3 of
!  the triangle are assumed to be such that the normal should be derived from
!  the cross product w = (v2 - v1) x (v3 - v1).  Then the unit normal vector
!  is n = w / |w| where |w| = sqrt (w.w) and the triangle area is 0.5 |w|.
!
!  07/19/06  D.A.Saunders  Initial modularization of unit normal calculation.
!  07/23/08    "     "     Added area as a second output since it requires much
!                          the same arithmetic.
!
!  Author:  David Saunders, ELORET/NASA Ames Research Center, Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)  :: nnode             ! # (x,y,z) points in tri_xyz(:,*)
   integer, intent (in)  :: ntri              ! # triangles defined in conn(:,*)
   real,    intent (in)  :: tri_xyz(3,nnode)  ! (x,y,z) coordinates pointed to
   integer, intent (in)  :: conn(3,ntri)      ! conn(1:3,itri) point to vertices
   integer, intent (in)  :: itri              ! ... 1, 2, 3 of triangle itri
   real,    intent (out) :: area              ! Triangle area
   real,    intent (out) :: unit_normal(3)    ! Desired unit normal vector

!  Local constants:

   real, parameter :: half = 0.5

!  Local variables:

   integer :: i1, i2, i3
   real    :: v1(3), v2(3), w(3), wsize

!  Execution:

   i1 = conn(1,itri);  i2 = conn(2,itri);  i3 = conn(3,itri)

   v1(1) = tri_xyz(1,i2) - tri_xyz(1,i1);  v2(1) = tri_xyz(1,i3) - tri_xyz(1,i1)
   v1(2) = tri_xyz(2,i2) - tri_xyz(2,i1);  v2(2) = tri_xyz(2,i3) - tri_xyz(2,i1)
   v1(3) = tri_xyz(3,i2) - tri_xyz(3,i1);  v2(3) = tri_xyz(3,i3) - tri_xyz(3,i1)

!! call cross (v1, v2, w)

   w(1)  =   v1(2) * v2(3) - v1(3) * v2(2)
   w(2)  = -(v1(1) * v2(3) - v1(3) * v2(1))
   w(3)  =   v1(1) * v2(2) - v1(2) * v2(1)

   wsize =    sqrt (dot_product (w, w))
   area  =    half * wsize
   unit_normal = w / wsize

   end subroutine tri_normal_and_area
