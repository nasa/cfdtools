!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine nearest_brick_point (ni, nj, nk, x, y, z, ic, jc, kc, xtarget, &
                                   xnear, coefs, dsq)
!
!  This is a structured-volume-grid variant of nearest_hex_point. It can be used
!  outside of the ADT (alternating digital tree) package that has its equivalent
!  in-line.  It includes the transfer of coordinates from separate X/Y/Z arrays
!  to the more convenient xvert(3,8), then follows nearest_hex_point except that
!  the vertex ordering is better suited to Fortran rather than being that of the
!  Fluent convention.  See the revised diagram below.
!
!  Input are the coordinates of some structured grid block, indices for some
!  cell in that block, and the coordinates of some target point.
!
!  Output is the point NOT OUTSIDE the indicated cell that is nearest to the
!  target point, and an associated squared distance.  This is accomplished by
!  constraining the underlying p/q/r coefficients to [0, 1] during the Newton
!  iteration involved.
!
!  Vertex numbering (Fortran convention, NOT Fluent convention):
!
!                       7 o-------------------o 8
!                        /|                  /|
!                       / |                 / |
!                      /  |                /  |
!                   5 o-------------------o 6 |
!                     |   |               |   |
!                     |   |               |   |
!                     | 3 o---------------|---o 4
!                     k  /                |  /
!                     | j                 | /
!                     |/                  |/
!                   1 o---------i---------o 2
!
!  06/10/13  D.A.Saunders  Initial implementation of nearest_hex_point.
!  07/30/13    "      "    Adaptation for hybrid method of FLOW_INTERP as
!                          nearest_brick_point.
!  08/05/13    "     "     Using ADT subroutine TERMINATE was a mistake once all
!                          ADT variants were combined into one package.
!
!  Author:  David Saunders, ERC Inc./NASA Ames Research Center, Moffett Fld., CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)  :: ni, nj, nk ! Grid block dimensions
   real,    intent (in), dimension (ni,nj,nk) :: x, y, z ! Block coordinates
   integer, intent (in)  :: ic, jc, kc ! Lower-left indices of the intended cell
   real,    intent (in)  :: xtarget(3) ! Target pt. being tested vs. this cell
   real,    intent (out) :: xnear(3)   ! Point found nearest to the target that
                                       ! is not outside this cell
   real,    intent (out) :: coefs(8)   ! Interpolation coefficients associated
                                       ! with xnear, in [0, 1] and summing to 1
   real,    intent (out) :: dsq        ! The corresponding shortest distance
                                       ! between xtarget and xnear, squared
!  Local constants:

   integer, parameter :: niter_max = 15
   real,    parameter :: half = 0.5, one = 1.0, zero = 0.0
   real,    parameter :: eps_stop = 1.e-8

!  Local variables:

   integer :: i, ierr, iter, j, k
   real :: p, q, r, pm1, qm1, rm1
   real, dimension (3) :: dp, dx, pl, pold, xp, xq, xr, xt, xpq, xqr, xrp, xpqr
   real, dimension (3,3) :: ajac
   real, dimension (3,8) :: xvert

!  Execution:

!  The following iteration is adapted from that of TRILINT.  Here,
!  p, q, r are forced to remain inside the unit cube, whereas TRILINT
!  allows convergence to a point outside the current cell as part of
!  the (misguided) RIPPLE3D strategy for determining inside/outside status.

   i = ic;  j = jc;  k = kc        ! Avoid excessive argument references

!! xvert(1,:) = x(i:i+1,j:j+1,k:k+1)  ! Doesn't compile.  :(

   xvert(1,1:2) = x(i:i+1,j,k)
   xvert(1,3:4) = x(i:i+1,j+1,k)
   xvert(1,5:6) = x(i:i+1,j,k+1)
   xvert(1,7:8) = x(i:i+1,j+1,k+1)
   xvert(2,1:2) = y(i:i+1,j,k)
   xvert(2,3:4) = y(i:i+1,j+1,k)
   xvert(2,5:6) = y(i:i+1,j,k+1)
   xvert(2,7:8) = y(i:i+1,j+1,k+1)
   xvert(3,1:2) = z(i:i+1,j,k)
   xvert(3,3:4) = z(i:i+1,j+1,k)
   xvert(3,5:6) = z(i:i+1,j,k+1)
   xvert(3,7:8) = z(i:i+1,j+1,k+1)

   xt(:)   = xvert(:,1) - xtarget(:)
   xp(:)   = xvert(:,2) - xvert(:,1)
   xq(:)   = xvert(:,3) - xvert(:,1)
   xr(:)   = xvert(:,5) - xvert(:,1)
   xpq(:)  = xvert(:,4) - xvert(:,3) - xp(:)
   xrp(:)  = xvert(:,6) - xvert(:,5) - xp(:)
   xqr(:)  = xvert(:,7) - xvert(:,5) - xq(:)
   xpqr(:) = xvert(:,8) - xvert(:,7) - xvert(:,6) + xvert(:,5) - xpq(:)

   pl(:) = half

   Newton:  do iter = 1, niter_max

      pold(:) = pl(:)

!     Jacobian (i, j) = partial df(i)/dp(j)

      ajac(:,1) = pl(2)*(pl(3)*xpqr(:) + xpq(:)) + pl(3)*xrp(:) + xp(:)
      ajac(:,2) = pl(3)*(pl(1)*xpqr(:) + xqr(:)) + pl(1)*xpq(:) + xq(:)
      ajac(:,3) = pl(1)*(pl(2)*xpqr(:) + xrp(:)) + pl(2)*xqr(:) + xr(:)

!     RHS elements are f (1), f (2), f (3):

      dp(:) = pl(1)*ajac(:,1) + pl(2)*(pl(3)*xqr(:) + xq(:)) + pl(3)*xr(:) + &
              xt(:)

!     The RHS elements can go to zero only if the target is inside the
!     cell because p,q,r are forced to stay inside the unit cube.

!!!   dsq = dp(1)**2 + dp(2)**2 + dp(3)**2

!!!   write (*, '(a, i4, a, i2, a, es9.2, a, 3e9.2, a, 3f12.9)') &
!!!     ' ii', ii, ' it', iter, ' dsq', dsq, ' xt', xtarget(:), ' pqr', pl(:)

      call lusolve (3, 3, ajac, dp, ierr)       ! Solve J dp = f

      if (ierr /= 0) then
         write (*, '(a)') ' NEAREST_BRICK_POINT:  Singular matrix.'
         dsq = 1.e+10
         coefs(:) = zero
         exit
      end if

      pl(:) = min (one, max (zero, pl(:) - dp(:)))
      dp(:) = abs (pl(:) - pold(:))

      if (maxval (dp(:)) < eps_stop) exit

   end do Newton

!  The deltas from the target (dp(:) before the solve) are one iteration
!  behind, so update the squared distance:

   ajac(:,1) = pl(2)*(pl(3)*xpqr(:) + xpq(:)) + pl(3)*xrp(:) + xp(:)
   dx(:) = pl(1)*ajac(:,1) + pl(2)*(pl(3)*xqr(:) + xq(:)) + pl(3)*xr(:) + xt(:)
   xnear(:) = xtarget(:) + dx(:)
   dsq = dot_product (dx, dx)

   p = pl(1);  pm1 = one - p
   q = pl(2);  qm1 = one - q
   r = pl(3);  rm1 = one - r

   coefs(1) = (rm1*qm1)*pm1
   coefs(2) = (rm1*qm1)*p
   coefs(3) = (rm1*q)*pm1
   coefs(4) = (rm1*q)*p

   coefs(5) = (r*qm1)*pm1
   coefs(6) = (r*qm1)*p
   coefs(7) = (r*q)*pm1
   coefs(8) = (r*q)*p

   end subroutine nearest_brick_point
