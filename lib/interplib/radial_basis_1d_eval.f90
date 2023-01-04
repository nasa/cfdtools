!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine radial_basis_1d_eval (ndat, xdat, w, shape, neval, xeval, yeval)

!  Description:
!
!     Evaluate the radial basis function-based function for which data point
!     weights have been computed by radial_basis_1d_weights, at one or more
!     abscissas xeval(:).  yeval(x) = sum w(i) g(||x - xdat(i)||), i = 1:ndat
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

   integer,                 intent (in)  :: ndat   ! Number of data points
   real, dimension (ndat),  intent (in)  :: xdat   ! Data point abscissas
   real, dimension (ndat),  intent (in)  :: w      ! Data point weights
   real,                    intent (in)  :: shape  ! Shape parameter used for w
   integer,                 intent (in)  :: neval  ! # abscissas to evaluate at
   real, dimension (neval), intent (in)  :: xeval  ! Evaluation abscissas
   real, dimension (neval), intent (out) :: yeval  ! Interpolations at xeval(:)

!  Local variables:

   integer :: i, ieval
   real    :: rsq, ssq, sum

!  ExecutionL

   ssq = shape**2

   do ieval = 1, neval
      sum = 0.
      do i = 1, ndat
         rsq = (xeval(ieval) - xdat(i))**2
         sum = sum + w(i) * exp (-ssq*rsq)
      end do
      yeval(ieval) = sum
   end do

   end subroutine radial_basis_1d_eval
