!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine solid_angle_quad (ni, nj, x, y, z, i, j, px, py, pz, domega)
!
!  Description:
!
!     For the quadrilateral cell (i:i+1,j:j+1) of a structured surface patch,
!     estimate the solid angle subtended at a given point P(px,py,pz).
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
!     local unit normal vector.  We simply average the four values of this
!     integrand calculated at the four cell vertices and multiply by the quad.
!     cell area, which is presumed to be small.  Existing utilities provide
!     that cell area and the unit normals.
!
!  History:
!
!     Feb. 19, 2014  D.A.Saunders  Initial implementation.
!     Feb. 21, 2014    "     "     The surface normal utility handles target
!                                  points at a patch boundary now, so there
!                                  is no need to safeguard it here.
!     Mar. 26, 2014    "     "     Radiance mentioned above (W.sr^-1.m^-2)
!                                  was incorrectly referred to as radiant
!                                  intensity (W.sr^-1).
!
!  Author:  David Saunders, ERC, Inc./NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)  :: ni, nj                 ! (Packed) patch array dims.
   real, intent (in), dimension (ni,nj) :: x, y, z ! Patch coordinates
   integer, intent (in)  :: i, j  ! i < ni, j < nj ! Indicated cell, lower left
   real,    intent (in)  :: px, py, pz             ! Coordinates of P
   real,    intent (out) :: domega                 ! Desired solid angle element

!  Local constants:

   real, parameter :: fourth = 0.25, p = 0.0, q = 0.0, zero = 0.0

!  Local variables;

   integer :: iv, jv
   real    :: rsq, sum, rvector(3), unormal(3)

!  Procedures:

   real     :: area4
   external :: area4           ! Quad cell area
   external :: surface_normal  ! Careful quad. cell unit normal at (p, q)

!  Execution:

!  Compute the integrand at each vertex of the given quad. cell:

   sum = zero

   do jv = j, j + 1

      do iv = i, i + 1
         rvector(1) = x(iv,jv) - px
         rvector(2) = y(iv,jv) - py
         rvector(3) = z(iv,jv) - pz

!        p, q are always zero because we're always at a cell vertex:

         call surface_normal (ni, nj, x, y, z, iv, jv, p, q, unormal)

!!!      write (*, '(a, 2i4, 3f12.8, 2x, 3f12.8)') &
!!!         'iv,jv,r,un:', iv,jv,rvector,unormal

         rsq = dot_product (rvector, rvector)
         sum = sum + dot_product (rvector, unormal) / (rsq * sqrt (rsq))
      end do
   end do

!  Integrand average x cell area:

   domega = fourth * sum * area4 (x(i,j), y(i,j), z(i,j),             &
                                  x(i+1,j), y(i+1,j), z(i+1,j),       &
                                  x(i+1,j+1), y(i+1,j+1), z(i+1,j+1), &
                                  x(i,j+1), y(i,j+1), z(i,j+1))

   end subroutine solid_angle_quad
