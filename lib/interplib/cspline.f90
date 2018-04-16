!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   module cspline_module
!
!  This module for conventional 1-D cubic spline interpolation is intended to
!  reduce clutter in an application for the common case of doing all of the
!  interpolations of a given dataset in a single call.  Moreover, as long as
!  datasets are not interleaved, further interpolations may be performed in
!  further calls without having to recalculate the coefficients, yet those
!  coefficients are hidden from the application.  End-condition controls may
!  also be initialized once and for all for the common case where one set
!  suffices in a given application (or changed as necessary, but effectively
!  hidden from where the action is).
!
!  The point is to make the coefficients private to the module.  Otherwise,
!  this is just a higher level interface permitted by Fortran 90 to simplify
!  use of some well-proven FORTRAN 77 subroutines.
!
!  By conventional spline, we mean the kind involving solution of a tridiagonal
!  system for the coefficients associated with one (x, y) dataset.  First and
!  second derivative continuity is enforced at the data points ("knots"), which
!  means geometric curvature is also continuous, if not necessarily SMOOTHLY so.
!
!  A choice of end conditions is provided, including periodic as might be needed
!  for a closed curve.  See the underlying Forsythe, Malcolm, and Moler-based
!  subroutines for all details, including the options for iendl, derivl, etc.
!
!  Abscissas may be increasing or decreasing, as long as they are monotonically
!  so.
!
!  Note that if spline coefficient arrays were included in a derived data type,
!  interpolation of different datasets could be interleaved. This implementation
!  chooses to address the simpler situation where datasets are not interleaved,
!  thus avoiding pointer arrays.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

   public

      integer :: iendl    ! Left  end condition control; iendl = 0|1|2|3|4
      real    :: derivl   ! Left  end derivative used if iendl =   1|2|3
      integer :: iendr    ! Right end condition control; iendr = 0|1|2|3|4
      real    :: derivr   ! Right end derivative used if iendr =   1|2|3

   real, allocatable, private, save, dimension (:) :: b, c, d  ! Spline coefs.

   public :: cspline      ! For the indicated dataset and target abscissas,
                          ! calculate spline coefficients (if new data) and
                          ! perform the interpolations

   external :: csfit      ! Underlying spline fitting ...
   external :: cseval     ! ... and evaluating routines, q.v.

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine cspline (new, ndata, x, y, neval, xeval, yeval, ier)

!     Streamlined interface to CSFIT and CSEVAL.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      logical, intent (in)  :: new          ! T => calculate spline coefficients
                                            !      and interpolate at xeval(:);
                                            ! F => perform more interpolations

      integer, intent (in)  :: ndata        ! Number of data points
      real,    intent (in)  :: x(ndata)     ! Data abscissas
      real,    intent (in)  :: y(ndata)     ! Data ordinates

      integer, intent (in)  :: neval        ! Number of target abscissas
      real,    intent (out) :: xeval(neval) ! Target abscissas
      real,    intent (out) :: yeval(neval) ! Spline interpolations

      integer, intent (out) :: ier          ! 0 means no error; see CSFIT for
                                            ! other possibilities, mostly to do
                                            ! with small ndata and out-of-range
                                            ! settings of iendl, iendr above
!     Local variables:

!     Execution:

      if (new) then

         if (allocated (b)) deallocate (b, c, d)

         allocate (b(ndata), c(ndata), d(ndata))

         call csfit (ndata, x, y, iendl, derivl, iendr, derivr, b, c, d, ier)

         if (ier /= 0) go to 99  ! Single-exit philosophy

      end if

      call cseval (ndata, x, y, neval, xeval, b, c, d, yeval)

  99  return

      end subroutine cspline
   
   end module cspline_module

