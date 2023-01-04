!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine radial_basis_1d_weights (n, x, y, shape, w, ier)

!  Description:
!
!     For an (x,y) dataset in one dimension, calculate weights w such that the
!     function y(x) = sum w(i) g(||x - x(i)||) interpolates the data points.
!     Following the Wikipedia radial basis function example, we choose
!     g(r) = exp (-(sr)**2) where s is a shape parameter such as 2 or 3.
!     The smooth function that the data points are assumed to represent can then
!     be interpolated at any x.  This is NOT a smoothing technique for noisy
!     data because the data points are interpolated exactly by construction.
!
!  Method:
!
!     The n x n linear system obtained by forcing the weighted sum of radial
!     basis terms to pass through the data points is solved by LU decomposition
!     with partial pivoting, and a check for apparent singularity that could
!     be a consequence of ill-conditioning.
!
!     To evaluate the computed interpolation function at any in-range x(s), see
!     companion subroutine radial_basis_1d_eval.
!
!  Reference:
!
!     Wikipedia, Radial basis function interpolation
!
!  History:
!
!     12/28/2022  D.A.Saunders  Initial implementation, out of curiosity.
!
!  Author:  David Saunders, AMA, Inc. at NASA Ames Research Center, CA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer,             intent (in)  :: n      ! Number of data points
   real, dimension (n), intent (in)  :: x, y   ! Data point coordinates
   real,                intent (in)  :: shape  ! Shape parameter > 0 (say 3.)
   real, dimension (n), intent (out) :: w      ! Computed weights for data pts.
   integer,             intent (out) :: ier    ! 0 => no singularity detected;
                                               ! 1 => apparent singularity
!  Local variables:

   integer :: i, j
   real    :: rsq, ssq
   real, allocatable :: A(:,:)

!  Procedures:

   external :: lusolve  ! LU decomposition and solution for one system

!  Execution:

   ssq = shape**2
   allocate (A(n,n))

   do j = 1, n
      do i = 1, n
         rsq = (x(j) - x(i))**2
         A(i,j) = exp (-ssq*rsq)
      end do
      w(j) = y(j)
   end do

   call lusolve (n, n, A, w, ier)

   deallocate (A)

   end subroutine radial_basis_1d_weights
