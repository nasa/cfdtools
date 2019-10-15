!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine flow_gradients_2d (ni, nj, nf, lxyz, x, y, f, fx, fy)
!
!  This is a 2D analogue of flow_gradients_nf, q.v.  The descriptions have been
!  adapted accordingly.
!
!  For all points of one structured grid block, calculate the partial first
!  derivatives of one or more functions with respect to x and y.
!
!  NOTE:  This routine differs from those in typical flow solvers, for which it
!  is unsuited, because the grid coordinate partial derivatives (metrics) aren't
!  separated from the derivatives of the functions on the computational grid.
!  Do not use it if each function is varying iteratively on a fixed grid.
!
!  This implementation avoids the traditional use of Cramer's rule, at some
!  cost, but that cost is intended to be a one-time thing for a given grid.
!  [Cramer's rule is notoriously poor as a means of solving a linear system.]
!
!  At each grid point, if xi, xj = partial dx/di, dx/dj on the usual rectangular
!  unit-spaced computational grid, and similarly for derivatives of y and f,
!  then the system to be solved for the x/y derivatives of f is (from the chain
!  rule):
!
!                         |xi  yi|  |fx|     |fi|
!                         |xj  yj|  |fy|  =  |fj|
!
!  The function values are assumed to be vertex-centered, not cell-centered, so
!  1-sided differencing is needed at the boundaries, which spoils vectorization.
!  (Three-point differencing is used everywhere.)  LU factorization of 2 x 2
!  matrices and partial pivoting for stability may be overkill for good grids,
!  but why not do it as well as possible for any grid?  This would not be the
!  approach to use in a typical PDE solver, but the functionality it provides
!  (prompted by shadowgraph calculations) is implemented in somewhat reusable
!  form in case other applications emerge.
!
!  01/28/2010  D. A. Saunders  Utility prompted by calculating shadowgraphs from
!                              CFD solutions (but the single-function form is
!                              all that's needed, and it's less inefficient).
!  02/10/2010     "      "     Multi-function form of FLOW_GRADIENTS (3D).
!  08/13/2018     "      "     2D form of FLOW_GRADIENTS_NF.
!
!  Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA
!           Now with AMA, Inc. at NASA Ames.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)  :: ni, nj        ! Grid block dimensions
   integer, intent (in)  :: nf            ! Number of functions >= 1
   logical, intent (in)  :: lxyz(2)       ! Allows suppression of fx/fy
   real,    intent (in)  :: x(ni,nj)      ! Grid coordinates
   real,    intent (in)  :: y(ni,nj)
   real,    intent (in)  :: f(nf,ni,nj)   ! Function values
   real,    intent (out) :: fx(nf,ni,nj)  ! Partial df/dx if lxyz(1)
   real,    intent (out) :: fy(nf,ni,nj)  ! Partial df/dy if lxyz(2)

!  Local constants:

   real, parameter :: three = 3., four = 4.

!  Local variables:

   integer :: i, i1, i2, ier, j, j1, j2, n, ip(2)
   real    :: A(2,2)

   real, allocatable :: b(:,:)

!  Procedures:

   external :: decomp, solve  ! Linear system solver for 2+ right-hand sides

!  Execution:

   allocate (b(2,nf))

!  Note that the halves resulting from uniform unit spacing in i/j space
!  cancel on the left- and right-hand sides of each 2 x 2 system.  Likewise,
!  at upper boundaries, the missing minus sign appears on both the left- and
!  right-hand sides, so it is not needed.

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
            A(1,1) = x(i+1,j) - x(i-1,j)
            A(1,2) = y(i+1,j) - y(i-1,j)
            b(1,:) = f(:,i+1,j) - f(:,i-1,j)
         else
            A(1,1) = -three*x(i,j) + four*x(i1,j) - x(i2,j)
            A(1,2) = -three*y(i,j) + four*y(i1,j) - y(i2,j)
            b(1,:) = -three*f(:,i,j) + four*f(:,i1,j) - f(:,i2,j)
         end if

         if (j1 == 0) then
            A(2,1) = x(i,j+1) - x(i,j-1)
            A(2,2) = y(i,j+1) - y(i,j-1)
            b(2,:) = f(:,i,j+1) - f(:,i,j-1)
         else
            A(2,1) = -three*x(i,j) + four*x(i,j1) - x(i,j2)
            A(2,2) = -three*y(i,j) + four*y(i,j1) - y(i,j2)
            b(2,:) = -three*f(:,i,j) + four*f(:,i,j1) - f(:,i,j2)
         end if

         call decomp (2, 2, A, ip)       ! LU factorization

!        Assume nonsingularity for now.

         do n = 1, nf

            call solve (2, 2, A, b(:,n), ip)  ! RHS is overwritten

            if (lxyz(1)) fx(n,i,j) = b(1,n)
            if (lxyz(2)) fy(n,i,j) = b(2,n)

         end do

      end do

   end do

   deallocate (b)

   end subroutine flow_gradients_2d
