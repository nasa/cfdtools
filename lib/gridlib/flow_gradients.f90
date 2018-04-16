!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine flow_gradients (ni, nj, nk, lxyz, x, y, z, f, fx, fy, fz)
!
!  For all points of one structured grid block, calculate the partial first
!  derivatives of a single function with respect to x, y, and/or z.
!
!  See also flow_gradients_nf for the multiple-functions case.
!
!  NOTE:  This routine differs from those in typical flow solvers, for which it
!  is unsuited, because the grid coordinate partial derivatives (metrics) aren't
!  separated from the derivatives of the function on the computational grid.
!  Do not use it if the function is varying iteratively on a fixed grid.
!
!  This implementation avoids the traditional use of Cramer's rule, at some
!  cost, but that cost is intended to be a one-time thing for a given grid.
!  [Cramer's rule is notoriously poor as a means of solving a linear system.]
!
!  At each grid point, if xi, xj, xk = partial dx/di, dx/dj, dx/dk on the usual
!  rectangular/unit-spaced computational grid, and similarly for derivatives of
!  y, z and f, then the system to be solved for the x/y/z derivatives of f is
!  (from the chain rule):
!
!                         |xi  yi  zi|  |fx|     |fi|
!                         |          |  |  |     |  |
!                         |xj  yj  zj|  |fy|  =  |fj|
!                         |          |  |  |     |  |
!                         |xk  yk  zk|  |fz|     |fk|
!
!  The function values are assumed to be vertex-centered, not cell-centered, so
!  1-sided differencing is needed at the boundaries, which spoils vectorization.
!  (Three-point differencing is used everywhere.)  LU factorization of 3 x 3
!  matrices and partial pivoting for stability may be overkill for good grids,
!  but why not do it as well as possible for any grid?  This would not be the
!  approach to use in a typical PDE solver, but the functionality it provides
!  for simulated shadowgraph calculations is implemented in somewhat reusable
!  form in case other applications emerge.
!
!  01/28/2010  D. A. Saunders  Utility needed for calculating shadowgraphs from
!                              CFD solutions.
!
!  Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)  :: ni, nj, nk    ! Grid block dimensions
   logical, intent (in)  :: lxyz(3)       ! Allows suppression of fx/fy/fz
   real,    intent (in)  :: x(ni,nj,nk)   ! Grid coordinates
   real,    intent (in)  :: y(ni,nj,nk)
   real,    intent (in)  :: z(ni,nj,nk)
   real,    intent (in)  :: f(ni,nj,nk)   ! Function values
   real,    intent (out) :: fx(ni,nj,nk)  ! Partial df/dx if lxyz(1)
   real,    intent (out) :: fy(ni,nj,nk)  ! Partial df/dy if lxyz(2)
   real,    intent (out) :: fz(ni,nj,nk)  ! Partial df/dz if lxyz(3)

!  Local constants:

   real, parameter :: three = 3., four = 4.

!  Local variables:

   integer :: i, i1, i2, ier, j, j1, j2, k, k1, k2
   real    :: A(3,3), b(3)

!  Procedures:

   external :: lusolve  ! Linear system solver for one right-hand side

!  Execution:

!  Note that the halves resulting from uniform unit spacing in i/j/k space
!  cancel on the left- and right-hand sides of each 3 x 3 system.  Likewise,
!  at upper boundaries, the missing minus sign appears on both the left- and
!  right-hand sides, so it is not needed.

   do k = 1, nk

      if (k == 1) then
         k1 = 2;     k2 = 3
      else if (k == nk) then
         k1 = nk-1;  k2 = k1-1
      else
         k1 = 0
      end if

      do j = 1, nj

         if (j == 1) then
            j1 = 2;     j2 = 3
         else if (j == nj) then
            j1 = nj-1;  j2 =j1-1 
         else
            j1 = 0
         end if

         do i = 1, ni

            if (i == 1) then
               i1 = 2;     i2 = 3
            else if (i == ni) then
               i1 = ni-1;  i2 =i1-1 
            else
               i1 = 0
            end if

            if (i1 == 0) then
               A(1,1) = x(i+1,j,k) - x(i-1,j,k)
               A(1,2) = y(i+1,j,k) - y(i-1,j,k)
               A(1,3) = z(i+1,j,k) - z(i-1,j,k)
               b(1)   = f(i+1,j,k) - f(i-1,j,k)
            else
               A(1,1) = -three*x(i,j,k) + four*x(i1,j,k) - x(i2,j,k)
               A(1,2) = -three*y(i,j,k) + four*y(i1,j,k) - y(i2,j,k)
               A(1,3) = -three*z(i,j,k) + four*z(i1,j,k) - z(i2,j,k)
               b(1)   = -three*f(i,j,k) + four*f(i1,j,k) - f(i2,j,k)
            end if

            if (j1 == 0) then
               A(2,1) = x(i,j+1,k) - x(i,j-1,k)
               A(2,2) = y(i,j+1,k) - y(i,j-1,k)
               A(2,3) = z(i,j+1,k) - z(i,j-1,k)
               b(2)   = f(i,j+1,k) - f(i,j-1,k)
            else
               A(2,1) = -three*x(i,j,k) + four*x(i,j1,k) - x(i,j2,k)
               A(2,2) = -three*y(i,j,k) + four*y(i,j1,k) - y(i,j2,k)
               A(2,3) = -three*z(i,j,k) + four*z(i,j1,k) - z(i,j2,k)
               b(2)   = -three*f(i,j,k) + four*f(i,j1,k) - f(i,j2,k)
            end if

            if (k1 == 0) then
               A(3,1) = x(i,j,k+1) - x(i,j,k-1)
               A(3,2) = y(i,j,k+1) - y(i,j,k-1)
               A(3,3) = z(i,j,k+1) - z(i,j,k-1)
               b(3)   = f(i,j,k+1) - f(i,j,k-1)
            else
               A(3,1) = -three*x(i,j,k) + four*x(i,j,k1) - x(i,j,k2)
               A(3,2) = -three*y(i,j,k) + four*y(i,j,k1) - y(i,j,k2)
               A(3,3) = -three*z(i,j,k) + four*z(i,j,k1) - z(i,j,k2)
               b(3)   = -three*f(i,j,k) + four*f(i,j,k1) - f(i,j,k2)
            end if

            call lusolve (3, 3, A, b, ier)  ! RHS is overwritten

!           Assume nonsingularity for now.

            if (lxyz(1)) fx(i,j,k) = b(1)
            if (lxyz(2)) fy(i,j,k) = b(2)
            if (lxyz(3)) fz(i,j,k) = b(3)

         end do

      end do

   end do

   end subroutine flow_gradients
