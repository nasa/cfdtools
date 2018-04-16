!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine nearest_hex_point (xvert, xtarget, xnear, coefs, dsq)

!  See subroutine NEAREST_CELL_POINT for a thorough description.  This initial
!  version of the hexahedral cell case simply adapts the nonlinear method used
!  for searching structured volume grids but returns coefs(1:8) in the form
!  described for the linear method that may or may not be implemented here.
!
!  The virtue of the nonlinear case is that constraining p, q, r to [0, 1]
!  during the Newton iteration produces a result that is not outside the cell,
!  as required.
!
!  Vertex numbering (Fluent convention):
!
!                       8 o-------------------o 7
!                        /|                  /|
!                       / |                 / |
!                      /  |                /  |
!                   5 o-------------------o 6 |
!                     |   |               |   |
!                     |   |               |   |
!                     | 4 o---------------|---o 3
!                     |  /                |  /
!                     | /                 | /
!                     |/                  |/
!                   1 o-------------------o 2
!
!  06/10/13  D.A.Saunders  Initial implementation.
!  08/05/13    "     "     Using ADT subroutine TERMINATE was a mistake once all
!                          ADT variants were combined into one package.
!
!  Author:  David Saunders, ERC Inc./NASA Ames Research Center, Moffett Fld., CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   real, intent (in)  :: xvert(3,8)  ! Cell vertex coordinates in Fluent order
   real, intent (in)  :: xtarget(3)  ! Target pt. being tested against this cell
   real, intent (out) :: xnear(3)    ! Point found closest to the target that is
                                     ! not outside this cell
   real, intent (out) :: coefs(8)    ! Interpolation coefficients associated
                                     ! with xnear, in [0, 1] and summing to 1
   real, intent (out) :: dsq         ! The corresponding shortest distance
                                     ! between xtarget and xnear, squared
!  Local constants:

   integer, parameter :: niter_max = 15
   real,    parameter :: half = 0.5, one = 1.0, zero = 0.0
   real,    parameter :: eps_stop = 1.e-8

!  Local variables:

   integer :: ierr, iter
   real :: p, q, r, pm1, qm1, rm1
   real, dimension (3) :: dp, dx, pl, pold, xp, xq, xr, xt, xpq, xqr, xrp, xpqr
   real, dimension (3,3) :: ajac

!  Execution:

!  The following iteration is adapted from that of TRILINT.  Here,
!  p, q, r are forced to remain inside the unit cube, whereas TRILINT
!  allows convergence to a point outside the current cell as part of
!  the (misguided) RIPPLE3D strategy for determining inside/outside status.

   xt(:)   = xvert(:,1) - xtarget(:)
   xp(:)   = xvert(:,2) - xvert(:,1)
   xq(:)   = xvert(:,4) - xvert(:,1)
   xr(:)   = xvert(:,5) - xvert(:,1)
   xpq(:)  = xvert(:,3) - xvert(:,4) - xp(:)
   xrp(:)  = xvert(:,6) - xvert(:,5) - xp(:)
   xqr(:)  = xvert(:,8) - xvert(:,5) - xq(:)
   xpqr(:) = xvert(:,7) - xvert(:,8) - xvert(:,6) + xvert(:,5) - xpq(:)

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

!!!   write (6, '(a, i4, a, i2, a, es9.2, a, 3e9.2, a, 3f12.9)') &
!!!     ' ii', ii, ' it', iter, ' dsq', dsq, ' xt', xtarget(:), ' pqr', pl(:)

      call lusolve (3, 3, ajac, dp, ierr)       ! Solve J dp = f

      if (ierr /= 0) then
         write (*, '(a)') ' NEAREST_HEX_POINT:  Singular matrix.'
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
   coefs(3) = (rm1*q)*p
   coefs(4) = (rm1*q)*pm1

   coefs(5) = (r*qm1)*pm1
   coefs(6) = (r*qm1)*p
   coefs(7) = (r*q)*p
   coefs(8) = (r*q)*pm1

   end subroutine nearest_hex_point
