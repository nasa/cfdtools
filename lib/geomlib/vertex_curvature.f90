!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine vertex_curvature (n, s, curvature, iv, nv, exponent, fraction)
!
!  Purpose:
!
!     Curvature distributions associated with ~rectilinear "curves" such as in
!  portions of a generatrix for an atmospheric entry capsule with vertices that
!  are not to be rounded are not handled well by standard curvature-based grid
!  point redistribution.  Such curvature distributions (which can be assumed to
!  come from geometry normalized to the unit square to get scaling out of the
!  picture) typically have huge spikes at each vertex point, and either zero or
!  small curvature on either side.  The technique of subroutine CURVDIS simply
!  cannot "see" such 1-point spikes, so the redistributed spacing is essentially
!  uniform where clustering towards each vertex is really desirable.
!
!     This utility is intended to make a given vertex more "visible" to the
!  CURVDIS technique by artificially broadening the curvature data on either
!  side of the vertex, without touching the underlying geometry.  The indicated
!  curvature data points are modified in-place by adding certain "sine bumps"
!  available from earlier airfoil shape optimization work.
!
!     The shape functions are appropriately denormalized so they match the
!  (scaled) height of the curvature spike above the curvature values at the
!  points iv +/- nv, where they go to zero.  Arguments nv and exponent both
!  control the broadening of the curvature spike:  nv should be fairly small
!  to localize the eventual clustering, and exponent should be at least 2. to
!  produce a shape with the intended curvature at both ends for smooth blending
!  with the unaltered curvature points.
!
!     Note that curvature may be positive or negative at this level, so spikes
!  with large negative curvature are allowed for.
!
!  History:
!
!     10/31/13  D.A.Saunders  Initial implementation, for CAPSULE_GRID,
!                             prompted by the Galileo probe.
!     11/01/13    "     "     Added the fraction argument, since the spike
!                             height is artificially high.
!
!  Author: David Saunders, ERC, Inc./NASA Ames Research Center, Moffett Fld., CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)    :: n             ! # data pts. in s(:) & curvature(:)
   real,    intent (in)    :: s(n)          ! Arc length along the curve
   real,    intent (inout) :: curvature(n)  ! Curvature before and after the
                                            ! treatment of the indicated vertex
   integer, intent (in)    :: iv            ! Index where a vertex has been IDd
   integer, intent (in)    :: nv            ! Only curvature(iv-nv+1:iv+nv-1)
                                            ! elements are modified; nv ~6:10
   real,    intent (in)    :: exponent      ! Exponent applied to COSL and COSR
                                            ! shape functions; exponent ~2.0
   real,    intent (in)    :: fraction      ! Multiplier applied to the spike
                                            ! height, which is artificially
                                            ! high; fraction ~0.2 - 0.25
!  Local constants:

   logical :: add = .false.  ! BEVAL won't be adding to earlier BEVAL calls

!  Local variables:

   integer :: i, ileft
   real    :: f1, f2, h1, h2, sleft, width, p(2), bnorm(nv+1), snorm(nv+1)

!  Execution:

!  We apply a scaled COSR left of the vertex:

   ileft = iv - nv
   sleft = s(ileft)
   width = s(iv) - sleft
   h1    = curvature(iv) - curvature(ileft)
   h2    = curvature(iv) - curvature(iv + nv)

   do i = 2, nv + 1
      snorm(i) = (s(i + ileft - 1) - sleft) / width
   end do

   p(1) = exponent   ! Controls width/shape of shape function
   f1   = fraction   ! Applies to left shape function only
   p(2) = f1*h1      ! Multiplier of unit-height shape function

   call beval ('COSR', 2, p, add, nv, snorm(2), bnorm(2))

   do i = 2, nv + 1
      curvature(ileft + i - 1) = curvature(ileft) + bnorm(i)
   end do

!  Likewise, apply a scaled COSL to the right of the vertex:

   sleft = s(iv)
   width = s(iv + nv) - s(iv)
   f2    = (h2 - h1*(1. - f1))/h2  ! Awkwardness if h2 /= h1
   p(2)  = f2*h2

   do i = 2, nv
      snorm(i) = (s(iv + i - 1) - sleft) / width
   end do

   call beval ('COSL', 2, p, add, nv-1, snorm(2), bnorm(2))

   do i = 2, nv
      curvature(iv + i - 1) = curvature(iv + nv) + bnorm(i)
   end do

   end subroutine vertex_curvature
