!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine nearest_quad_point (xvert, xtarget, xnear, coefs, dsq)

!  See subroutine NEAREST_CELL_POINT for a thorough description.  This initial
!  version of the quadrilateral cell case simply adapts the nonlinear method of
!  searching structured surface grids but returns coefs(1:4) in the form
!  described for the linear method that may or may not be implemented here.
!
!  The virtue of the nonlinear case is that constraining p, q to [0, 1] during
!  the Newton iteration produces a result that is not outside the cell, as reqd.
!
!  Vertex numbering (Fluent convention):
!
!                       4 o-------------------o 3
!                        /                   /
!                       /                   /
!                      /                   /
!                   1 o-------------------o 2
!
!  06/13/13  David Saunders  Initial implementation.
!
!  Author:  David Saunders, ERC Inc./NASA Ames Research Center, Moffett Fld., CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   real, intent (in)  :: xvert(3,4)  ! Cell vertex coordinates in Fluent order
   real, intent (in)  :: xtarget(3)  ! Target pt. being tested against this cell
   real, intent (out) :: xnear(3)    ! Point found closest to the target that is
                                     ! not outside this cell
   real, intent (out) :: coefs(4)    ! Interpolation coefficients associated
                                     ! with xnear, in [0, 1] and summing to 1
   real, intent (out) :: dsq         ! The corresponding shortest distance
                                     ! between xtarget and xnear, squared
!  Local constants:

   integer, parameter :: niter_max = 15
   real,    parameter :: fourth = 0.25, half = 0.5, one = 1.0, zero = 0.0
   real,    parameter :: eps = 1.e-25, eps_stop = 1.e-8

!  Local variables:

   integer :: iter
   real :: pq, pm1, qm1, vn
   real, dimension (2) :: dp, pl, pold
   real, dimension (3) :: a, an, b, bn, dx, n, vf, vt, xf, &
                          x21, x41, x32, x3142

!  Execution:

   pl(:) = half;  pq = fourth

   x21(:)   = xvert(:,2) - xvert(:,1)
   x41(:)   = xvert(:,4) - xvert(:,1)
   x32(:)   = xvert(:,3) - xvert(:,2)
   x3142(:) = x32(:)     - x41(:)
   xf(:)    = xvert(:,1) + pl(1)*x21(:) + pl(2)*x41(:) + pq*x3142(:)

   Newton:  do iter = 1, niter_max

      pold(:) = pl(:)

      a(:) = x21(:) + pl(2)*x3142(:)  ! Tangent vector in q direction
      b(:) = x41(:) + pl(1)*x3142(:)  ! Tangent vector in p direction

      n(1) = a(2)*b(3) - a(3)*b(2)    ! Cross product of a and b
      n(2) = a(3)*b(1) - a(1)*b(3)
      n(3) = a(1)*b(2) - a(2)*b(1)
      n(:) = n(:) / max (eps, sqrt (n(1)**2 + n(2)**2 + n(3)**2))  ! Unit normal

      vf(:) = xtarget(:) - xf(:)      ! Vector from xf(:) to the target
      vn    = vf(1)*n(1) + vf(2)*n(2) + vf(3)*n(3)
      vt(:) = vf(:) - vn*n(:)         ! Projection of vf(:) onto the quad

      an(1) = a(2)*n(3) - a(3)*n(2)
      an(2) = a(3)*n(1) - a(1)*n(3)
      an(3) = a(1)*n(2) - a(2)*n(1)

      bn(1) = b(2)*n(3) - b(3)*n(2)
      bn(2) = b(3)*n(1) - b(1)*n(3)
      bn(3) = b(1)*n(2) - b(2)*n(1)

!     Solve for delta p and q.  See the structured surface grid form for more.

      vn    = a(1)*bn(1) + a(2)*bn(2) + a(3)*bn(3)
      vn    = sign (max (eps, abs (vn)), vn)
      dp(1) = (vt(1)*bn(1) + vt(2)*bn(2) + vt(3)*bn(3))/vn

      vn    = b(1)*an(1) + b(2)*an(2) + b(3)*an(3)
      vn    = sign (max (eps, abs (vn)), vn)
      dp(2) = (vt(1)*an(1) + vt(2)*an(2) + vt(3)*an(3))/vn

      pl(:) = min (one, max (zero, pl(:) + dp(:)))
      dp(:) = abs (pl(:) - pold(:))

      pq    = pl(1)*pl(2)
      xf(:) = xvert(:,1) + pl(1)*x21(:) + pl(2)*x41(:) + pq*x3142(:)

      if (maxval (dp(:)) < eps_stop) exit

   end do Newton

   xnear(:) = xf(:)
   dx(:)    = xf(:) - xtarget(:)
   dsq      = dx(1)**2 + dx(2)**2 + dx(3)**2

   pm1 = one - pl(1)
   qm1 = one - pl(2)

   coefs(1) = qm1*pm1
   coefs(2) = qm1*pl(1)
   coefs(3) = pl(2)*pl(1)
   coefs(4) = pl(2)*pm1

   end subroutine nearest_quad_point
