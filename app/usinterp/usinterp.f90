!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program usinterp
!
!  Purpose:
!
!     Read a list of data points with any number of dimensions and any number of
!     function values at each data point as a column dataset with any number of
!     header lines that are ignored.  Also read a similar table of coordinates
!     at which to calculate estimates of those functions.  This serves as a
!     driver for the optimal interpolation module originally made available by
!     Alexander Barth (University of Liege, Belgium) and extended by David
!     Saunders at NASA Ames Research Center.  Interpolating/smoothing of an
!     unstructured dataset is intended.
!
!  Method:
!
!     The optimal interpolation scheme is designed to minimize interpolation
!     errors under certain statistical assumptions about the function samples.
!     Experience with smoothly varying test datasets shows remarkably accurate
!     results that may be artificially good.  This driver is part of setting up
!     and running further test cases, including working with noisy datasets.
!
!  Unstructured Input Dataset Format:
!
!     Header records (zero or more lines that are ignored)
!     x1 x2 ... xn  f1 f2 .. fm  ! First data point: n dimensions, m functions
!     x1 x2 ... xn  f1 f2 .. fm  ! Further data points, any order
!      :  : ...  :   :  :     :  ! n >= 1, m >= 1
!      :  : ...  :   :  :     :
!     x1 x2 ... xn  f1 f2 .. fm
!
!  Target Coordinates Dataset Format:
!
!     Header records (one or more lines, ignored)
!     x1 x2 ... xn
!     x1 x2 ... xn
!      :  : ...  :  ! Any number of target points, in some order or not
!      :  : ...  :
!     x1 x2 ... xn
!
!  Input Control File Format (on Standard Input):
!
!     USINTERP Control File
!     nd                    ! Number of coordinate dimensions; nd >= 1
!     nf                    ! Number of functions at each data point
!     unstructured.dat      ! Input dataset
!     target.dat            ! Target coordinates dataset
!     usinterp.dat          ! Output results
!     (3es16.8, 4es15.7)    ! Run-time output format
!     m                     ! Number of nearest neighbors to interpolate within
!     normalize             ! T|F; T => normalize the input coordinates
!     verbose               ! T|F; T => diagnostic output is written to fort.50
!     cf1 cf2 ...           ! nd correlation functions: each is g|l|G|L
!     cl1 cl2 ...           ! nd correlation lengths in (0, 1]
!     ovar1 ovar2 ...       ! nf observation variances >= 0.
!     adapt                 ! -1.|0.|a < 1.; -1. => simple inverse dist. weights
!
!  Output Dataset Format:
!
!     Same as the target coordinates, with interpolated function(s) appended
!     as in the input data format, no header lines, and interpolated function
!     variances as nf additional functions:
!
!     x1 x2 ... xn  f1 f2 .. fm v1 v2 .. vm
!     x1 x2 ... xn  f1 f2 .. fm v1 v2 .. vm
!      :  : ...  :  ! Target coords., interpolated fns., interpolated variances
!      :  : ...  :
!     x1 x2 ... xn  f1 f2 .. fm v1 v2 .. vm
!
!  Notes:
!
!     1. Best choices for the multiple controls tend to require trial and error
!        because of inadequate documentation.
!     2. The normalization option converts the working coordinates to [0, 1]
!        in case of disparate raw units; unnormalized coordinates are output.
!     3. Forcing some nonzero contribution from the furthest of m neighbors to
!        a given target point may or may not be a good idea.  Increasing the
!        correlation length(s) should have a similar effect.
!     4. See the available test_optiminterp.f90 for testing on a dataset that
!        is constructed from a trigonometric function, with an option to impose
!        Gaussian noise on the data point function values.
!     5. The observation variance inputs are assumed to be a measure of the
!        uncertainty in the function data.  Originally, they could differ at
!        each data point for a single function.  Handling of multiple functions
!        led to one variance per function type instead.  Experience suggests
!        that specifying zeros for the variances is preferable, but this is
!        not well understood.
!     6. The adaptation option is not working as intended because the trig.
!        squared test case suggests that large correlation lengths produce the
!        best results.
!
!  History:
!
!     01/31/2022  D.A.Saunders  Initial implementation to enable more testing.
!     02/11/2022    "     "     [After a severed internet cable hiatus:]
!                               trig-squared test case for Optimal_Interpolation
!                               reproduced here from tabular datasets.
!     02/27/2022    "     "     Perform a second set of interpolations at the
!                               data point coordinates and show the sum of
!                               squared differences from the data functions.
!
!  Author: David Saunders, AMA, Inc. at NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use table_type_module      ! Both are in table_io.f90
   use table_io
   use optimal_interpolation  ! Unstructured data interpolation/smoothing

   implicit none

!  Constants:

   integer, parameter :: &
      lundata  = 1,      &  ! Input dataset
      luntarg  = 2,      &  ! Input target coordinates
      lunout   = 3,      &  ! Output interpolations
      lunctl   = 5,      &  ! Keyboard inputs
      luncrt   = 6          ! Screen prompts and diagnostics

!  Local variables:

   integer :: &
      i, ios, j, l, m, nf, nd, ndata, ninterp

   real :: &
      adapt, sumsq

   real, allocatable, dimension (:) :: &
      cl, ovar                         ! Correlation length(s) and
                                       ! observation variances

   logical :: &
      normalize, verbose

   character (1), allocatable, dimension (:) :: &
      cf                                        ! Correlation function(s)

   type (table_type) :: &
      table_data, table_interp, table_target, table_interp_at_data_pts

!  Execution:

!  Process the control file:

   read (lunctl, '(a)')  ! Skip header
   read (lunctl, *) nd   ! Dimension/number of independent coordinates
   read (lunctl, *) nf   ! Number of functions at each data point

   read (lunctl, *) table_data%filename
   call table_io_read_real (lundata, table_data, ios)
   if (ios /= 0) go to 99

   ndata = table_data%nrows
   if (table_data%ncols /= nd + nf) then
      write (luncrt, '(a, i7)') &
         'Expected number of data columns:', nd + nf, &
         'Number of data columns read:    ', table_data%ncols
      go to 99
   end if

   read (lunctl, *) table_target%filename
   call table_io_read_real (luntarg, table_target, ios)
   if (ios /= 0) go to 99

   ninterp = table_target%nrows
   if (table_target%ncols /= nd) then
      write (luncrt, '(a, i7)') &
         'Expected number of target columns:', nd, &
         'Number of data columns read:      ', table_target%ncols
      go to 99
   end if

   read (lunctl, *) table_interp%filename
   read (lunctl, '(a)') table_io_format
   l = index (table_io_format, '!')
   if (l > 0) table_io_format(l:) = ' '; l = len_trim (table_io_format)
   table_interp%nrows   = ninterp
   table_interp%ncols   = table_data%ncols + nf ! Include interpolated variances
   table_interp%nheader = 0

   allocate (cf(nd), cl(nd), ovar(nf))

   read  (lunctl, *) m  ! # nearest neighbors to interpolate from at each target
   read  (lunctl, *) normalize
   read  (lunctl, *) verbose
   read  (lunctl, *) cf(:)   ! Correlation functions, one per dimension
   read  (lunctl, *) cl(:)   ! Correlation lengths,    "   "   "   "
   read  (lunctl, *) ovar(:) ! Measures of function uncertainty (?), 1 per fn.
   read  (lunctl, *) adapt   ! Explain

   if (verbose) then  ! The optimal interpolation package writes to unit 50
      rewind (lunctl);  call transcribe (lunctl, 50, ios)
   end if

   close (lunctl)

!  Allocate storage for the interpolation results:

   allocate (table_interp%values(nd+nf+nf,ninterp))

   do i = 1, ninterp
      table_interp%values(1:nd,i) = table_target%values(:,i)
   end do

!  More storage for further interpolations at the data coordinates:

   table_interp_at_data_pts%nheader = 1
   allocate (table_interp_at_data_pts%header(1))
   allocate (table_interp_at_data_pts%values(nd+nf+nf,ndata))
   table_interp_at_data_pts%values(1:nd,:) = table_data%values(1:nd,:)

!  Perform the interpolations:

   call optiminterp (table_data%values(1:nd,:), &
                     table_data%values(nd+1:nd+nf,:), ovar, cl, adapt, cf, m, &
                     table_target%values, normalize, verbose, &
                     table_interp%values(nd+1:nd+nf,:), &
                     table_interp%values(nd+nf+1:nd+nf+nf,:))

!! call optiminterp (x, f, var, cor_len, adapt, fun, m, xi, normalize, &
!!                   verbose, fi, vari)

   call table_io_write_real (lunout, table_interp, table_io_format_sep, ios)
   if (ios /= 0) go to 99

!  Perform a second set of interpolations at the data points and compute
!  a measure of their goodness:

   call optiminterp (table_data%values(1:nd,:), &
                     table_data%values(nd+1:nd+nf,:), ovar, cl, adapt, cf, m, &
                     table_interp_at_data_pts%values(1:nd,:), &
                     normalize, verbose, &
                     table_interp_at_data_pts%values(nd+1:nd+nf,:), &
                     table_interp_at_data_pts%values(nd+nf+1:nd+nf+nf,:))

   table_interp_at_data_pts%nrows = ndata
   table_interp_at_data_pts%ncols = nd + nf + nf
   table_interp_at_data_pts%filename  = 'interpolations_at_data_pts.dat'
   table_interp_at_data_pts%header(1) = 'Squared difference sums(s):'

   do j = 1, nf
      sumsq = 0.
      do i = 1, ndata
         sumsq = sumsq + (table_interp_at_data_pts%values(nd+j,i) - &
                          table_data%values(nd+j,i))**2
      end do
      l = len_trim (table_interp_at_data_pts%header(1))
      write (table_interp_at_data_pts%header(1)(l+1:), '(es16.8)') sumsq
   end do

   if (verbose) write (50, '(/, a)') table_interp_at_data_pts%header(1)(1:l+16)

   call table_io_write_real (lunout, table_interp_at_data_pts, &
                             table_io_format_sep, ios)

99 continue

   end program usinterp
