!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine lcsfit2 (ndata, x, y, method, neval, xeval, yeval)
!
!  One-liner:  Variant of LCSFIT that normalizes the data first
!
!  Description:
!
!     LCSFIT2 was prompted by the need to redistribute multiple 1-D
!     profiles of flow field data with potentially huge units (species
!     number densities in m^-3 and temperatures in Kelvin).  It simply
!     normalizes the data, performs the interpolations, then denormalizes
!     the results.  The input data are NOT normalized in place.  For this
!     reason, no attempt is made to avoid normalizing the abscissas after
!     a call for the first of multiple profiles associated with the same
!     coordinates.
!
!     Note that "x" is commonly arc length along a curve in 2- or 3-space
!     and "y" is a function rather than a coordinate in the usual sense.
!
!     The NEW (data) argument of LCSFIT has been removed because the NEVAL
!     input is unlikely to be 1.  Likewise, the first derivatives provided
!     by LCSFIT are suppressed as almost certainly redundant.
!
!     See LCSFIT for further details of the local cubic spline method it
!     uses (specially the monotonic option).
!
!  History:
!
!     06/06/1998  R.A.Kennelly,  Version of LCSFIT started with for LCSFIT2.
!                 D.A.Saunders
!     11/22/2013  D.A.Saunders   Initial complement to LCSFIT, intended for
!                                data with abnormally large or small units.
!
!  Author:  David Saunders, ERC Inc./NASA Ames Research Center, Moffett Fld., CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer,   intent (in)  :: ndata               ! # data points
   real,      intent (in)  :: x(ndata), y(ndata)  ! Data coordinates, with x(:)
                                                  ! monotonically increasing or
                                                  ! decreasing
   character, intent (in)  :: method*1            ! Uppercase fit type:
                                                  ! 'M', 'B', 'L', 'C' as for
                                                  ! LCSFIT, q.v.
   integer,   intent (in)  :: neval               ! # interpolations requested;
                                                  ! neval >= 1
   real,      intent (in)  :: xeval(neval)        ! Abscissa(s) to interpolate
                                                  ! to, normally within the
                                                  ! range of x(:)
   real,      intent (out) :: yeval(neval)        ! Interpolated function values

!  Local constants:

   integer, parameter :: ndim = 2       ! Meaning x and y are normalized
   logical, parameter :: new  = .true.  ! Always new data assumed for LCSFIT

!  Local variables:

   integer                 :: ier
   real                    :: scale(ndim), shift(ndim)
   real, dimension (ndata) :: xdata, ydata    ! Avoid overwriting input data
   real, dimension (neval) :: xevaln, ypeval  ! Unused 1st derivs. at xevaln(:)

!  Procedures:


!  Execution:

   xdata(:) = x(:)
   ydata(:) = y(:)

!  Normalize the data x(:) & y(:) independently; they are probably s(:) & f(:):

   call getscale ('I', ndim, ndata, xdata, ydata, ydata, scale, shift, ier)
   call usescale ('N', ndim, ndata, xdata, ydata, ydata, scale, shift, ier)

!  Apply the same normalization to xeval(:) as has been applied to xdata(:):

   xevaln(:) = xeval(:)

   call usescale ('N', 1, neval, xevaln, xevaln, xevaln, scale, shift, ier)

!  Interpolate in normalized space:

   call lcsfit (ndata, xdata, ydata, new, method, neval, xevaln, yeval, ypeval)

!  Denormalize the interpolations:

   call usescale ('D', 1, neval, yeval, yeval, yeval, scale(ndim), &
                  shift(ndim), ier)

   end subroutine lcsfit2
