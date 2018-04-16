!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine solid_angle_tri (nnodes, ntris, xyz, conn, itri, Pxyz, &
                               area, domega)
!  Description:
!
!     For cell itri of a triangulated surface zone, estimate the solid angle
!     subtended at a given point Pxyz.
!
!     This was prompted by a need to integrate radiance vs. solid angle in
!     order to estimate radiative heat flux at a point more correctly than
!     with the tangent-slab approximation.
!
!  Technique:
!
!     Wikipedia's Solid Angle web page shows that the solid angle for an
!     arbitrary oriented surface element S subtended at a point P is
!
!                    omega = integral over S of (r.n/r^3) ds
!
!     where r is the vector position of a point on S w.r.t. P and n is the
!     local unit normal vector.  We simply average the three values of this
!     integrand calculated at the cell vertices and multiply by the cell
!     area, which is presumed to be small.  Existing utilities provide that
!     cell area and the unit normals.
!
!     Unlike the structured grid case of quadrilaterals, we can't (easily)
!     identify cell neighbors in order to estimate a refined unit normal at
!     each cell vertex.  Therefore, we simply employ the same normal for
!     each vertex.
!
!  History:
!
!     Feb. 25, 2014  D.A.Saunders  Adaptation of solid_angle_quad.
!     Feb. 26, 2014    "     "     Work with three r vectors instead of
!                                  one averaged r vector.
!     Mar. 25, 2014    "     "     Use Pxyz(1:3) instead of px, py, pz.
!     Mar. 26, 2014    "     "     Radiance mentioned above (W.sr^-1.m^-2)
!                                  was incorrectly referred to as radiant
!                                  intensity (W.sr^-1).
!
!  Author:  David Saunders, ERC, Inc./NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)  :: nnodes         ! # x/y/z points/nodes
   integer, intent (in)  :: ntris          ! # triangular elements
   real,    intent (in)  :: xyz(3,nnodes)  ! Node coordinates
   integer, intent (in)  :: conn(3,ntris)  ! Vertex indices in xyz(:,*)
   integer, intent (in)  :: itri           ! Indicated triangular cell
   real,    intent (in)  :: Pxyz(3)        ! Coordinates of P
   real,    intent (out) :: area           ! Cell area as a by-product
   real,    intent (out) :: domega         ! Desired solid angle element

!  Local constants:

   real, parameter :: third = 1.0/3.0, zero = 0.0

!  Local variables;

   integer :: inode, iv
   real    :: rsq, rvector(3), sum, unormal(3)

!  Procedures:

   external :: tri_normal_and_area

!  Execution:

!  Compute a unit normal to the triangular element, and its area:

   call tri_normal_and_area (nnodes, ntris, xyz, conn, itri, area, unormal)

!  Compute the averaged r vector and apply it to the common unit normal.
!  Later: no - this approaches the right answer from above with grid refinement.
!  Instead, average three values of the integrand, even though they all use the
!  same unit normal.  This approaches the right answer from below as expected.

   sum = zero

   do iv = 1, 3
      inode = conn(iv,itri)
      rvector(:) = xyz(:,inode) - Pxyz(:)
      rsq = dot_product (rvector, rvector)
      sum = sum + dot_product (rvector, unormal) / (rsq * sqrt (rsq))
   end do

   domega = third * sum * area

   end subroutine solid_angle_tri
