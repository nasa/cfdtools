!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine tet_volumes (nnode, ntet, xyz, conn, volume, volume_total)
!
!  Calculate the volumes of the tetrahedra of one zone of an unstructured volume
!  grid.  One third of a face area times the altitude to the fourth vertex is
!  employed for the volume calculations, in favor of a determinant formula.
!
!  Also return, for each vertex, the summed volumes of tetrahedra sharing that
!  vertex, as needed for volume-weighted interpolation of function data from
!  cell centers to cell vertices.
!
!  Programming Note:
!
!  Determinant formulas tend to become ill-conditioned for extreme cases.
!  See www.eecs.berkeley.edu/~wkahan/VtetLang.pdf for an excellent treatment
!  of the issues in calculating areas and volumes accurately yet efficiently.
!  This routine avoids one determinant formula but employs the area formula
!  inherited from TRI_AREAS, so it is neither one thing nor the other.  Since
!  the intended application is for weighted interpolation of cell-centered
!  functions to vertices on computational grids where the cells should have
!  sensible aspect ratios, the potential weakness of the area formula should
!  not be an issue, while the convenience of reusing the existing PROJECT3 was
!  considered preferable to evaluating a second, more elaborate, determinant.
!  If performance becomes a problem (two external references per cell), consider
!  the 4 x 4 or 3 x 3 determinant formulas readily available on the web.
!
!  03/17/2010  D.A.Saunders  Initial implementation of TRI_AREAS.
!  04/27/2010    "     "     Analogous adaptation of TET_VOLUMES.
!
!  Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)  :: nnode               ! # vertices/nodes in the zone
   integer, intent (in)  :: ntet                ! # tetrahedra in the zone
   real,    intent (in)  :: xyz(3,nnode)        ! (x,y,z)s of the vertices
   integer, intent (in)  :: conn(4,ntet)        ! Vertex pointers into xyz(:,*)
   real,    intent (out) :: volume(ntet)        ! Triangle volumes
   real,    intent (out) :: volume_total(nnode) ! Summed volume surrounding each
                                                ! vertex/node
!  Local constants:

   real, parameter :: sixth = 1./6., zero = 0.

!  Local variables:

   integer :: i1, i2, i3, i4, itet
   real    :: dsq, four_area_squared, p, q, xfoot(3)
   real    :: x1, x2, x3, y1, y2, y3, z1, z2, z3

!  Execution:

   volume_total(:) = zero

   do itet = 1, ntet
      i1 = conn(1,itet);  x1 = xyz(1,i1);  y1 = xyz(2,i1);  z1 = xyz(3,i1)
      i2 = conn(2,itet);  x2 = xyz(1,i2);  y2 = xyz(2,i2);  z2 = xyz(3,i2)
      i3 = conn(3,itet);  x3 = xyz(1,i3);  y3 = xyz(2,i3);  z3 = xyz(3,i3)
      i4 = conn(4,itet)

      four_area_squared = &
         (y1*(z2 - z3) - z1*(y2 - y3) + y2*z3 - y3*z2) ** 2 +                  &
         (z1*(x2 - x3) - x1*(z2 - z3) + z2*x3 - z3*x2) ** 2 +                  &
         (x1*(y2 - y3) - y1*(x2 - x3) + x2*y3 - x3*y2) ** 2

!     Squared distance from vertex 4 to this face via orthogonal projection:

      call project3 (xyz(:,i1), xyz(:,i2), xyz(:,i3), xyz(:,i4), xfoot, dsq,   &
                     p, q)

      volume(itet) = sixth * sqrt (four_area_squared * dsq)

      volume_total(i1) = volume_total(i1) + volume(itet)
      volume_total(i2) = volume_total(i2) + volume(itet)
      volume_total(i3) = volume_total(i3) + volume(itet)
      volume_total(i4) = volume_total(i4) + volume(itet)
   end do

   end subroutine tet_volumes
