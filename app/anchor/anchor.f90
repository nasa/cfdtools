!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   program anchor

!  INTRODUCTION:
!
!     Anchoring refers to determining corrections for low-fidelity calculations
!  from a limited number of more expensive higher-fidelity calculations.  Where
!  the high-fidelity "observations" are available, the corrected results should
!  match, while far away from them the adjustments should be essentially zero.
!
!     A related issue is determining where to compute the high-fidelity results.
!  The short answer is to add them where gradients in the preliminary estimates
!  are highest.  (Actually, in 1- or 2-D, adding fidelity in the high-curvature
!  regions makes most sense. These correspond roughly to high 2nd derivatives.)
!  This aspect is not covered here (yet).
!
!     This application is prompted by the need to correct Euler calculations of
!  aerodynamic lift and drag coefficients using a limited number of more costly
!  Navier-Stokes solutions.   It is implemented in a general enough way that it
!  may find other uses.
!
!     The data points are treated as unstructured, and the interpolation method
!  is a form of Kriging known as Optimal Interpolation.  Since both the low and
!  high fidelity datasets are being interpolated, it may seem that the datasets
!  should have independent variance estimates. However, only low-fidelity error
!  estimates are meaningful because it is actually DIFFERENCES between high and
!  low fidelity data that are involved in the second set of interpolations. The
!  same estimates are applied to all points of a given dataset, meaning they're
!  crude at best.    Provision is made for specifying different error variances
!  for the different functions (originally assumed to be comparable).
!
!  ASSUMPTIONS:
!
!     o  Allow for any number of coordinate dimensions and any number of
!        functions at each point.
!     o  Read and write one data point per line; count the number of input data
!        points by reading to EOF.
!     o  Outputs are assumed to be desired as a rectangular table defined by a
!        list of coordinates in each dimension (original design; this version
!        can alternatively read an unstructured list of target points at which
!        to interpolate, in which case derivative outputs are not an option).
!        For rectangular-table mode, the table coordinates are read from the
!        control file, with each dimension starting on a new line.
!        For the arbitrary-list-of-points-at-which-to-interpolate mode, those
!        target points are read one tuple per line until EOF from the file
!        indicated on the control line following some meaningful rectangular
!        table definition (acting as a place holder, but ignored after being
!        read).  This line is skipped if the table definition is to be used,
!        but it must be present.
!     o  Results are tabulated in Fortran order (first coordinate varying first,
!        etc.) if they are in rectangular table form, else they are tabulated
!        in the input target point order.  Plottable results are also written
!        in Tecplot format with POINT order.  The file name is derived from the
!        tabulation name by adding ".Tec.dat".  This file includes any high-
!        fidelity points as a second zone.
!
!  OUTLINE:
!
!     o  A short control file specifies file names and interpolation parameters.
!        It is read from standard input.
!     o  Interpolate the lo-fi functions at the hi-fi coordinates and at the
!        target coordinates (concatenated to do both in one call).
!     o  Interpolate the  (hi-fi - lo-fi) function differences at the target
!        output coordinates.
!     o  Add the interpolated corrections to the lo-fi functions evaluated at
!        the target coordinates.
!
!  NO-ANCHORING OPTION:
!
!     If the high fidelity file name is "none", only the initial interpolations
!     of the low-fi data are performed.  This provides a way of testing control
!     inputs on one set of interpolations, and/or performing Kriging-type
!     interpolation within a single dataset (the "low-fi" data) without the
!     "anchoring" idea in the picture.
!
!  SAMPLE CONTROL FILE:
!
!     --------------------------------------------------------------------------
!     ANCHOR Control File
!     --------------------------------------------------------------------------
!     2  3              ! # coordinates and # functions
!     Mach Alpha        ! coordinate names
!     CL CD CM          ! function names
!     cbaero.dat        ! low fidelity data input file
!     overflow.dat      ! high fidelity data input file or "none"
!     anchored.dat      ! output results file (structured table or just a list)
!     derivatives.dat   ! output derivatives (plottable) or "none"
!     1                 ! target point mode: 1 = rectangular table; 2 = list
!     14  52            ! dimensions of output table, mode = 1; coords. follow:
!     0.5 0.7 0.9 0.95 1.05 1.1 1.2 1.6 2.0 2.5 3.0 4.0 6.0 8.0
!     0 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100 105 110 115
!     120 125 130 132 134 136 138 140 142 144 146 148 150 152 154 156 158 160
!     162 164 166 168 170 172 174 176 176 180
!     target.dat        ! file name of target pts., mode = 2; ignore if mode = 1
!     --------------------------------------------------------------------------
!     Low-fidelity Data Interpolation Controls
!     --------------------------------------------------------------------------
!     40                ! # nearest neighbors used for each interpolation
!     4.   40.          ! data correlation lengths in coordinate space
!     G    G            ! correlation functions: G = Gaussian, L = Linear
!     0.1               ! adaptation control in [0, ~0.5]; 0 = no adaptation
!     1.e-6 1.e-6 1.e-6 ! function error variances >= 0
!     --------------------------------------------------------------------------
!     High-fidelity Data Interpolation Controls
!     --------------------------------------------------------------------------
!     40                ! # nearest neighbors used for each interpolation
!     4.   40.          ! data correlation lengths in coordinate space
!     G    G            ! correlation functions: G = Gaussian, L = Linear
!     0.1               ! adaptation control in [0, ~0.5]; 0 = no adaptation
!     0.0   0.0   0.0   ! function error variances >= 0
!
!  HISTORY:
!
!    01/17/07- D.A.Saunders  Initial design ...
!    01/19/07                ... and implementation.
!    01/26/07    "     "     Replaced PLOT3D outputs with Tecplot format;
!                            separated the lo-fi & hi-fi interpolation controls.
!    02/06/07    "     "     Skip the anchoring procedures if the hi-fi dataset
!                            name matches the lo-fi dataset name, to allow the
!                            interpolation method to be tested more easily.
!    02/07/07    "     "     Linear correlation function option appears to work
!                            better than Gaussian for CL, CD data.
!    02/13/07    "     "     The Optimal Interpolation package now provides for
!                            a choice of correlation functions, and the scheme
!                            for dealing with wide data-spacing variations by
!                            adapting the correlation lengths is no longer hard-
!                            coded.
!    03/14/07    "     "     Added calculation of derivatives of the structured
!                            interpolations via finite differencing. The hope is
!                            that they will highlight where high-fidelity data
!                            points are desirable, since the error variances
!                            produced by the optimal interpolation procedure
!                            seem to indicate little other than where the data
!                            points are sparse (which could be intentional if
!                            the functions are ~constant or linear in those
!                            regions).
!    03/15/07    "     "     Omitting the function data from the plottable file
!                            of derivatives and curvature was unwise - we may
!                            want to plot (say) CD vs. Mach and Alpha but color
!                            the surface with contours of (say) d^2CD/dMach^2.
!    05/19/07    "     "     Added the option to interpolate at an arbitrary
!                            list of points read from a file.  Some rectangular
!                            table definition should remain in the control file
!                            as a place holder, but it will be ignored after
!                            being read if "output mode" = 2.  In this case, the
!                            target list of coordinates is read one tuple per
!                            line from the file named on a new control line
!                            following the table coordinates (for which each
!                            dimension should start on a new line).
!    11/14/08    "     "     The optimal interpolation package now handles
!                            different observation error variances for the
!                            different functions.
!    08/31/16    "     "     The OI package now has a "verbose" argument, and
!                            an option to use a simple inverse-distance-weighted
!                            method (via a negative "adapt").
!    09/20/22    "     "     A normalize argument has been added to optiminterp,
!                            among other changes therein.  See also the more
!                            recent USINTERP for interpolating scattered data.
!
!  Author:  David Saunders, ELORET/NASA Ames Research Center, Moffett Field, CA
!           Now with AMA, Inc. at NASA ARC.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use optimal_interpolation ! See optimal_interpolation.f90

   implicit none

!  Constants:

   integer, parameter ::  &
      luncrt = 6,         &  ! For diagnostics
      lunctl = 5,         &  ! Control file (standard input)
      lunin  = 1,         &  ! Input data files
      lunout = 2             ! Tabulated results and plottable form of them

   logical, parameter ::  &
      normalize = .true., &  ! Wise when the independent variables are disparate
      verbose = .false.      ! The input scheme doesn't lend itself to adding
                             ! this new control 
!  Variables:

   integer :: &
      i, inc, ios, j, k, l, mode_target_points, n, ncombined, ndim, nf, nhifi, &
      nlofi, ntarget, num_neighbors_lofi, num_neighbors_hifi

   integer, allocatable, dimension (:) :: &
      length_dimensions, &  ! # digits needed to print each dimension
      table_dimensions      ! # output table entries in each coord. direction

   real :: &
      adapt_lofi, adapt_hifi

   real, allocatable, dimension (:) :: &
      cor_lengths_lofi,   &  ! Correlation lengths in each dimension for ...
      cor_lengths_hifi,   &  ! ... the low-fidelity and high-fidelity data
      variances_lofi,     &  ! Function error variances, 1/function (lo-fi data)
      variances_hifi,     &  ! .................................... (hi-fi data)
      x_table                ! Rectangular coordinates for interpolated results

   real, allocatable, dimension (:,:) :: &
      x_lofi, f_lofi,     &  ! Low fidelity data to be corrected
      x_hifi, f_hifi,     &  ! Higher fidelity data
      x_interp,           &  ! At which lo-fi data are to be interpolated
      f_lo_interp,        &  ! Corresponding interpolations
      f_hi_minus_lo,      &  ! Differences at hi-fi points
      x_target,           &  ! Output table coordinates (expanded)
      f_target,           &  ! Desired table of corrected function values
      variances_interp       ! Interpolated function error variances, minimized

   real, allocatable, dimension (:,:,:) :: &
      f1, f2, fk            ! Corresponding 1st & 2nd derivatives & curvature-
                            ! like quantity
   logical :: &
      derivatives, skip_anchoring

   character :: &
      filename_lofi * 64, filename_hifi * 64, filename_derivatives * 64,       &
      filename_results * 64, filename_targets * 64

   character, allocatable, dimension (:) :: &
      coordinate_names * 16, function_names * 16, &
      cor_functions_lofi * 1, cor_functions_hifi * 1

!  Execution:
!  ----------

!  Read controls, open files, allocate work=space, and read input data:

   call initialize ()  ! Internal procedure below

   if (ios /= 0) go to 99

!  Option to test the interpolation scheme with a single set of calculations:
!  --------------------------------------------------------------------------

   if (skip_anchoring) then ! Keep it simple

      allocate (f_target(nf,ntarget), variances_interp(nf,ntarget))

      call optiminterp (x_lofi, f_lofi, variances_lofi, cor_lengths_lofi,      &
                        adapt_lofi, cor_functions_lofi, num_neighbors_lofi,    &
                        x_target, normalize, verbose, f_target,                &
                        variances_interp)

      if (derivatives) then

         allocate (f1(nf,ntarget,ndim), f2(nf,ntarget,ndim), &
                   fk(nf,ntarget,ndim))

         call calculate_derivatives ()  ! Internal procedure

      end if

      call save_results ()  ! Internal procedure

      deallocate (cor_lengths_lofi, cor_lengths_hifi, cor_functions_lofi,      &
                  cor_functions_hifi, variances_lofi, variances_hifi,          &
                  variances_interp, x_table, x_lofi, f_lofi,                   &
                  x_target, f_target, coordinate_names, function_names,        &
                  length_dimensions, table_dimensions)

      if (derivatives) deallocate (f1, f2, fk)

      go to 99

   end if

!  Normal anchoring case:
!  ----------------------

!  Append the target coordinates to the hi-fi coordinates so that all the
!  interpolations of the lo-fi data can be done in one call:

   ncombined = nhifi + ntarget

   allocate (x_interp(ndim,ncombined))

   do n = 1, ndim
      x_interp(n,1:nhifi) = x_hifi(n,:)
      x_interp(n,nhifi+1:ncombined) = x_target(n,:)
   end do

   allocate (f_lo_interp(nf,ncombined), variances_interp(nf,ncombined))

!  Interpolate the low-fi data at the hi-fi points and the table points:

   call optiminterp (x_lofi, f_lofi, variances_lofi, cor_lengths_lofi,         &
                     adapt_lofi, cor_functions_lofi, num_neighbors_lofi,       &
                     x_interp, normalize, verbose, f_lo_interp,                &
                     variances_interp)

!  Corrections at the high-fidelity coordinates:

   allocate (f_hi_minus_lo(nf,nhifi))

   do n = 1, nf
      f_hi_minus_lo(n,:) = f_hifi(n,:) - f_lo_interp(n,1:nhifi)
   end do

   open  (9, file='f_hi_minus_lo.dat', status='unknown')
   write (9, '(a)') 'TITLE = "High- minus low-fidelity data at hifi points"'
   do n = 1, nhifi
      write (9, '(1p, 10e13.5)') x_hifi(:,n), f_hi_minus_lo(:,n)
   end do
   close (9)

!  Interpolate the corrections at the output table coordinates:

   deallocate (variances_interp);  allocate (variances_interp(nf,ntarget))

   allocate (f_target(1:nf,1:ntarget))

   call optiminterp (x_hifi, f_hi_minus_lo, variances_hifi, cor_lengths_hifi,  &
                     adapt_hifi, cor_functions_hifi, num_neighbors_hifi,       &
                     x_target, normalize, verbose, f_target, variances_interp)

   open  (9, file='f_corrections_at_target.f', status='unknown')
   write (9, '(i1)') 1
   write (9, '(4i4)') ntarget, 1, 1, nf
   do n = 1, nf
      write (9, '(1p, 10e13.5)') f_target(n,:)
   end do
   close (9)

!  Apply the interpolated corrections to the lo-fi interpolations, table coords.

   do n = 1, nf
      f_target(n,:) = f_lo_interp(n,nhifi+1:ncombined) + f_target(n,:)
   end do

   if (derivatives) then

      allocate (f1(nf,ntarget,ndim), f2(nf,ntarget,ndim), fk(nf,ntarget,ndim))

      call calculate_derivatives ()  ! Internal procedure

   end if

!  Tabulate results and save them in plottable form:

   call save_results ()

   deallocate (cor_lengths_lofi, cor_lengths_hifi, cor_functions_lofi,         &
               cor_functions_hifi, variances_lofi, variances_hifi,             &
               variances_interp, x_table, x_lofi, f_lofi,                      &
               x_hifi, f_hifi, x_interp, f_lo_interp, f_hi_minus_lo,           &
               x_target, f_target, coordinate_names, function_names,           &
               length_dimensions, table_dimensions)

   if (derivatives) deallocate (f1, f2, fk)

99 continue ! Avoid system-dependent "stop" behavior.

!  Internal procedures for program anchor, in order of execution:

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine initialize ()  ! From standard input

!     Most variables are inherited from the host by the internal procedures.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

!     Execution:

      do l = 1, 3
         read (lunctl, '(a)', iostat=ios)  ! Skip header lines
         if (ios /= 0) then
            write (luncrt, '(/, a)') ' Usage:  anchor < xxx.inp'
            go to 99  ! Single return philosophy
         end if
      end do

      read (lunctl, *, iostat=ios) ndim, nf
      if (ios /= 0) then
         write (luncrt, '(/, a)') 'Trouble reading ndim & nf.'
         go to 99
      end if

      allocate (coordinate_names(ndim), function_names(2*nf), stat=ios)

!     (1 extra function name for the error variances.)

      if (ios /= 0) then
         write (luncrt, '(/, a)') 'Trouble allocating variable names.'
         go to 99
      end if

      read (lunctl, *, iostat=ios) coordinate_names
      if (ios /= 0) then
         write (luncrt, '(/, a)') 'Trouble reading coordinate names.'
         go to 99
      end if

      read (lunctl, *, iostat=ios) function_names(1:nf)
      if (ios /= 0) then
         write (luncrt, '(/, a)') 'Trouble reading function names.'
         go to 99
      end if
      do l = 1, nf  ! Add interpolated function variances to plot file
         function_names(nf+l) = 'Errvar_' // function_names(l)
      end do

      read (lunctl, *, iostat=ios) filename_lofi
      if (ios /= 0) then
         write (luncrt, '(/, a)') 'Trouble reading low-fidelity file name.'
         go to 99
      end if

      read (lunctl, *, iostat=ios) filename_hifi 
      if (ios /= 0) then
         write (luncrt, '(/, a)') 'Trouble reading high-fidelity file name.'
         go to 99
      end if

      skip_anchoring = filename_hifi(1:4) == 'none'

      read (lunctl, *, iostat=ios) filename_results
      if (ios /= 0) then
         write (luncrt, '(/, a)') 'Trouble reading output results file name.'
         go to 99
      end if

      read (lunctl, *, iostat=ios) filename_derivatives
      if (ios /= 0) then
         write (luncrt, '(/, a)') 'Trouble reading output derivatives filename.'
         go to 99
      end if

      derivatives = filename_derivatives(1:4) /= 'none'

      read (lunctl, *, iostat=ios) mode_target_points
      if (ios /= 0) then
         write (luncrt, '(/, a)') 'Trouble reading target point mode (1 or 2).'
         go to 99
      end if

!     A rectangular table definition should be present as a place-holder even if
!     mode_target_point = 2:

      allocate (table_dimensions(ndim), length_dimensions(ndim))

      read (lunctl, *, iostat=ios) table_dimensions
      if (ios /= 0) then
         write (luncrt, '(/, a)') 'Trouble reading output table dimensions.'
         go to 99
      end if

      ntarget = sum (table_dimensions)

      allocate (x_table(ntarget))

      i = 1;  l = 0
      do n = 1, ndim  ! Pack the rectangular coordinates
         l = l + table_dimensions(n)
         read (lunctl, *, iostat=ios) x_table(i:l)
         if (ios /= 0) then
            write (luncrt, '(/, a)') 'Trouble reading rectangular coordinates.'
            go to 99
         end if
         i = l + 1
      end do

      read (lunctl, *, iostat=ios) filename_targets
      if (ios /= 0) then
         write (luncrt, '(/, a)') 'Trouble reading target points filename.'
         go to 99
      end if

      do l = 1, 3
         read (lunctl, '(a)', iostat=ios)  ! Skip more header lines
         if (ios /= 0) then
            write (luncrt, '(/, a)') 'Trouble skipping lo-fi interp. headers'
            go to 99
         end if
      end do

      read (lunctl, *, iostat=ios) num_neighbors_lofi
      if (ios /= 0) then
         write (luncrt, '(/, a)') 'Trouble reading # nearest lo-fi neighbors.'
         go to 99
      end if

      allocate (cor_lengths_lofi(ndim), variances_lofi(nf))

      read (lunctl, *, iostat=ios) cor_lengths_lofi
      if (ios /= 0) then
         write (luncrt, '(/, a)') 'Trouble reading lo-fi correlation lengths.'
         go to 99
      end if

      allocate (cor_functions_lofi(ndim))

      read (lunctl, *, iostat=ios) cor_functions_lofi
      if (ios /= 0) then
         write (luncrt, '(/, a)') 'Trouble reading lo-fi correlation functions.'
         go to 99
      end if

      read (lunctl, *, iostat=ios) adapt_lofi
      if (ios /= 0) then
         write (luncrt, '(/, a)') 'Trouble reading lo-fi adaptation control.'
         go to 99
      end if

      read (lunctl, *, iostat=ios) variances_lofi
      if (ios /= 0) then
         write (luncrt, '(/, a)') 'Trouble reading lo-fi function variances.'
         go to 99
      end if

      do l = 1, 3
         read (lunctl, '(a)', iostat=ios)  ! Skip more header lines
         if (ios /= 0) then
            write (luncrt, '(/, a)') 'Trouble skipping hi-fi interp. headers'
            go to 99
         end if
      end do

      read (lunctl, *, iostat=ios) num_neighbors_hifi
      if (ios /= 0) then
         write (luncrt, '(/, a)') 'Trouble reading # nearest hi-fi neighbors.'
         go to 99
      end if

      allocate (cor_lengths_hifi(ndim), variances_hifi(nf))

      read (lunctl, *, iostat=ios) cor_lengths_hifi
      if (ios /= 0) then
         write (luncrt, '(/, a)') 'Trouble reading hi-fi correlation lengths.'
         go to 99
      end if

      allocate (cor_functions_hifi(ndim))

      read (lunctl, *, iostat=ios) cor_functions_hifi
      if (ios /= 0) then
         write (luncrt, '(/, a)') 'Trouble reading hi-fi correlation functions.'
         go to 99
      end if

      read (lunctl, *, iostat=ios) adapt_hifi
      if (ios /= 0) then
         write (luncrt, '(/, a)') 'Trouble reading hi-fi adaptation control.'
         go to 99
      end if

      read (lunctl, *, iostat=ios) variances_hifi
      if (ios /= 0) then
         write (luncrt, '(/, a)') 'Trouble reading hi-fi function variances.'
         go to 99
      end if

      close (lunctl)

!     End of control file processing.
!     Now open and read the data files.

!     Low-fidelity data:
!     ------------------

      l = len_trim (filename_lofi)
      open (lunin, file=filename_lofi(1:l), iostat=ios)
      if (ios /= 0) then
         write (luncrt, '(/, a)') 'Trouble opening ', filename_lofi(1:l)
         go to 99
      end if

!     Count the lines:

      n = 0
      do ! Until EOF
         read (lunin, '(a)', iostat=ios)
         if (ios /= 0) exit
         n = n + 1
      end do

      rewind (lunin)
      nlofi = n
      allocate (x_lofi(ndim,n), f_lofi(nf,n))

      do l = 1, n
         read (lunin, *, iostat=ios) x_lofi(:,l), f_lofi(:,l)
         if (ios /= 0) then
            write (luncrt, '(/, a)') 'Trouble reading lo-fi data.'
            go to 99
         end if
      end do

      close (lunin)

!     High-fidelity data:
!     -------------------

      if (skip_anchoring) then
         nhifi = 0
      else
         l = len_trim (filename_hifi)
         open (lunin, file=filename_hifi(1:l), iostat=ios)
         if (ios /= 0) then
            write (luncrt, '(/, a)') 'Trouble opening ', filename_hifi(1:l)
            go to 99
         end if

         n = 0  ! Count the lines
         do ! Until EOF
            read (lunin, '(a)', iostat=ios)
            if (ios /= 0) exit
            n = n + 1
         end do

         rewind (lunin)
         nhifi = n
         allocate (x_hifi(ndim,n), f_hifi(nf,n))

         do l = 1, n
            read (lunin, *, iostat=ios) x_hifi(:,l), f_hifi(:,l)
            if (ios /= 0) then
               write (luncrt, '(/, a)') 'Trouble reading hi-fi data.'
               go to 99
            end if
         end do

         close (lunin)
      end if

!     Deal with two cases for the target points to interpolate at:
!     ------------------------------------------------------------

      select case (mode_target_points)

         case (1) ! Rectangular table

!           Expand the rectangular table of output coordinates, needed in packed
!           form for the first set of interpolations of the lo-fi data:

            ntarget = product (table_dimensions)

            allocate (x_target(ndim,ntarget))

            select case (ndim)

               case (1)
                  x_target(1,:) = x_table(:)

               case (2)
                  call expand_2d (ndim, &
                                  table_dimensions(1), table_dimensions(2),    &
                                  x_table, x_target)
               case (3)
                  call expand_3d (ndim, &
                                  table_dimensions(1), table_dimensions(2),    &
                                  table_dimensions(3), x_table, x_target)
               case (4)
                  call expand_4d (ndim, &
                                  table_dimensions(1), table_dimensions(2),    &
                                  table_dimensions(3), table_dimensions(4),    &
                                  x_table, x_target)
               case (5)
                  call expand_5d (ndim, &
                                  table_dimensions(1), table_dimensions(2),    &
                                  table_dimensions(3), table_dimensions(4),    &
                                  table_dimensions(5), x_table, x_target)
               case (6)
                  call expand_6d (ndim, &
                                  table_dimensions(1), table_dimensions(2),    &
                                  table_dimensions(3), table_dimensions(4),    &
                                  table_dimensions(5), table_dimensions(6),    &
                                  x_table, x_target)
               case default
                  write (luncrt, '(/, a, a)') &
                     ' Too many dimensions: extend table_coordinates.f90,',    &
                     ' anchor.f90, and table_derivatives.f90.'
                  ios = 99
                  go to 99

            end select

         case (2) ! Arbitrary list of target points, one per line

            l = len_trim (filename_targets)
            open (lunin, file=filename_targets(1:l), iostat=ios)
            if (ios /= 0) then
               write (luncrt,'(/, a)') 'Trouble opening ', filename_targets(1:l)
               go to 99
            end if

            n = 0  ! Count the lines
            do ! Until EOF
               read (lunin, '(a)', iostat=ios)
               if (ios /= 0) exit
               n = n + 1
            end do

            rewind (lunin)
            ntarget = n
            allocate (x_target(ndim,ntarget))

            do l = 1, n
               read (lunin, *, iostat=ios) x_target(:,l)
               if (ios /= 0) then
                  write (luncrt, '(/, a)') 'Trouble reading target coords.'
                  go to 99
               end if
            end do

            close (lunin)

            derivatives = .false. ! They assume a rectangular table of functions
            table_dimensions(:) = 1
            table_dimensions(1) = ntarget

         case default

            write (luncrt, '(/, a, i6)') &
               ' Invalid input for target point mode: ', mode_target_points
            ios = 99

      end select

!     Determine the numbers of digits needed to print the table dimensions:

      do n = 1, ndim
         call ndigits (table_dimensions(n), length_dimensions(n))
      end do

   99 return

      end subroutine initialize

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine calculate_derivatives ()

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      select case (ndim)

        case (1)

           call derivs_1d (ntarget, nf, x_target, f_target, f1, f2, fk)

        case (2)

           call derivs_2d (ndim, table_dimensions(1), table_dimensions(2),     &
                           nf, x_target, f_target, f1, f2, fk)
        case (3)

           call derivs_3d (ndim, table_dimensions(1), table_dimensions(2),     &
                           table_dimensions(3), nf, x_target, f_target,        &
                           f1, f2, fk)
        case (4)

           call derivs_4d (ndim, table_dimensions(1), table_dimensions(2),     &
                           table_dimensions(3), table_dimensions(4),           &
                           nf, x_target, f_target, f1, f2, fk)
        case (5)

           call derivs_5d (ndim, table_dimensions(1), table_dimensions(2),     &
                           table_dimensions(3), table_dimensions(4),           &
                           table_dimensions(5), nf, x_target, f_target,        &
                           f1, f2, fk)

!!      case (6) ! Too many subscripts in derivs_6d - deal with it if we have to

!!         call derivs_6d (ndim, table_dimensions(1), table_dimensions(2),     &
!!                         table_dimensions(3), table_dimensions(4),           &
!!                         table_dimensions(5), table_dimensions(6),           &
!!                         nf, x_target, f_target, f1, f2, fk)
      end select

      end subroutine calculate_derivatives

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine save_results ()

!     Write a tabulation of the anchored results, and a plottable form of it.
!     The plottable file has the error variances as an extra function.
!     If anchoring is being skipped, use the tabulation file name as a Tecplot
!     file name.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      integer,   parameter :: numdigits = 13        ! Width of tabulated columns
      character, parameter :: indices * 7 = 'IJKLMN?'

      integer   :: len, k, l, m
      character :: buffer * 256, c * 1, iformat * 4

!     Execution:

      l = len_trim (filename_results)
      open (lunout, file=filename_results(1:l), status='unknown')

      if (skip_anchoring) then ! No tabulation

         write (lunout, '(a)') 'TITLE = "Optimal Interpolation"'

      else ! Tabulate as well as save in Tecplot form

!        Make the tabulation self-descriptive:

         write (lunout, '(a)') 'Anchored results'
         write (lunout, '(i3, a)') ndim, '  ! # dimensions'
         write (lunout, '(i3, a)')   nf, '  ! # functions'
         write (lunout, *) table_dimensions(:)  ! array dimensions

         buffer = ' ';  l = 0
         do n = 1, ndim
            len = len_trim (coordinate_names(n))
            l = l + numdigits
            i = l - len + 1
            buffer(i:l) = coordinate_names(n)(1:len)
         end do

         do n = 1, nf
            len = len_trim (function_names(n)) 
            l = l + numdigits
            i = l - len + 1
            buffer(i:l) = function_names(n)(1:len)
         end do
         
         write (lunout, '(a)') buffer(1:l)

         do i = 1, ntarget
            write (lunout, '(1p, 10e13.5)') x_target(:,i), f_target(:,i)
         end do

         close (lunout)

!        Derive a Tecplot file name from the tabulated results file name:

         l = len_trim (filename_results)
         i = index (filename_results(1:l), '.')
         if (i == 0) then
            i = l + 1;  filename_results(i:i) = '.'
         end if
         i = i + 1

         filename_results(i:)  = ' ';  l = i + 6
         filename_results(i:l) = 'Tec.dat'
         open  (lunout, file=filename_results(1:l), status='unknown')

         write (lunout, '(a)') 'TITLE = "Anchored data"'

      end if

      buffer = ' ';  i = 2
      do n = 1, ndim + nf + nf  ! The extra nf is for the error variances
         buffer(i:i) = '"'
         if (n <= ndim) then
            len = len_trim (coordinate_names(n))
         else
            len = len_trim (function_names(n-ndim))
         end if
         l = i + len;  i = i + 1
         if (n <= ndim) then
            buffer(i:l) = coordinate_names(n)(1:len)
         else
            buffer(i:l) = function_names(n-ndim)(1:len)
         end if
         l = l + 1;  buffer(l:l+1) = '",'
         i = l + 3
      end do
      write (lunout, '(2a)') 'VARIABLES =', buffer(1:l)

      if (skip_anchoring) then
         write (lunout,  '(a)') 'ZONE T = "Interpolated data"'
      else
         write (lunout,  '(a)') 'ZONE T = "Anchored data"'
      end if

!     Does this work for more than 3 dimensions??

      if (ndim > 6) then
         write (luncrt, '(/, a, i3)') ' Too many dimensions? ', ndim
         write (luncrt, '(a)') ' Check output Tecplot zone dimensions.'
      end if

      iformat = '(in)'
      i = 1;  buffer = ' '

      do n = 1, ndim
         l = min (n, 7)
         c = indices(l:l)
         buffer(i:i) = c;  i = i + 1
         buffer(i:i) = '='
         l = length_dimensions(n)
         write (iformat(3:3), '(i1)') l
         write (buffer(i+1:i+l), iformat) table_dimensions(n)
         i = i + l + 1;  buffer(i:i+1) = ', '
         i = i + 2
      end do

      write (lunout, '(a)') buffer(1:i-3)
      write (lunout, '(a)') 'DATAPACKING=POINT'

      do i = 1, ntarget
         write (lunout, '(1p, 15e13.5)') &
            x_target(:,i), f_target(:,i), variances_interp(:,i)
      end do

      buffer(3:) = ' '  ! Leave I= intact

      if (skip_anchoring) then

         write (lunout,  '(a)') 'ZONE T = "Data points"'
         call ndigits (nlofi, l)
         write (iformat(3:3), '(i1)') l
         write (buffer(3:2+l), iformat) nlofi
         i = 3 + l;  buffer(i:i+1) = ', '
         i = i + 2

         do n = 2, ndim
            l = min (n, 7)
            c = indices(l:l)
            buffer(i:i) = c;  i = i + 1
            buffer(i:i+3) = '=1, '
            i = i + 4
         end do

         write (lunout, '(a)') buffer(1:i-3)

         do i = 1, nlofi
            write (lunout, '(1p, 10e13.5)') &
               x_lofi(:,i), f_lofi(:,i), variances_lofi(:)
         end do

      else ! Anchoring case:  just show the hi-fi data

         write (lunout,  '(a)') 'ZONE T = "High fidelity data"'
         call ndigits (nhifi, l)
         write (iformat(3:3), '(i1)') l
         write (buffer(3:2+l), iformat) nhifi
         i = 3 + l;  buffer(i:i+1) = ', '
         i = i + 2

         do n = 2, ndim
            l = min (n, 7)
            c = indices(l:l)
            buffer(i:i) = c;  i = i + 1
            buffer(i:i+3) = '=1, '
            i = i + 4
         end do

         write (lunout, '(a)') buffer(1:i-3)

         do i = 1, nhifi
            write (lunout, '(1p, 10e13.5)') &
               x_hifi(:,i), f_hifi(:,i), variances_hifi(:)
         end do

      end if

      close (lunout)

!     Derivatives as well?

      if (derivatives) then

         open (lunout, file=filename_derivatives, status='unknown')

         write (lunout, '(a)') 'TITLE = "Derivatives and curvatures"'

         buffer = ' ';  i = 2
         do n = 1, ndim
            buffer(i:i) = '"'
            len = len_trim (coordinate_names(n))
            l = i + len;  i = i + 1
            buffer(i:l) = coordinate_names(n)(1:len)
            l = l + 1;  buffer(l:l+1) = '",'
            i = l + 3
         end do

!        Include the functions, as plotting them from another file is awkward.

         do n = 1, nf
            buffer(i:i) = '"'
            len = len_trim (function_names(n))
            l = i + len;  i = i + 1
            buffer(i:l) = function_names(n)(1:len)
            l = l + 1;  buffer(l:l+1) = '",'
            i = l + 3
         end do

!        Should check for exceeding buffer length (but don't).

         do m = 1, 3 ! f', f", fk
            if (m == 1) then
               c = '1'
            else if (m == 2) then
               c = '2'
            else
               c = 'k'
            end if
            do n = 1, nf
               do k = 1, ndim
                  buffer(i:i) = '"'
                  len = len_trim (function_names(n))
                  l = i + len;  i = i + 1
                  buffer(i:l) = function_names(n)(1:len)
                  l = l + 1
                  buffer(l:l) = c
                  l = l + 1
                  write (buffer(l:l), '(i1)') k
                  l = l + 1;  buffer(l:l+1) = '",'
                  i = l + 3
               end do
            end do
         end do

         write (lunout, '(2a)') 'VARIABLES =', buffer(1:l)

         if (skip_anchoring) then
            write (lunout, '(a)') 'ZONE T = "Interpolated data and derivatives"'
         else
            write (lunout, '(a)') 'ZONE T = "Anchored data and derivatives"'
         end if

!        Does this work for more than 3 dimensions??

         i = 1;  buffer = ' '

         do n = 1, ndim
            l = min (n, 7)
            c = indices(l:l)
            buffer(i:i) = c;  i = i + 1
            buffer(i:i) = '='
            l = length_dimensions(n)
            write (iformat(3:3), '(i1)') l
            write (buffer(i+1:i+l), iformat) table_dimensions(n)
            i = i + l + 1;  buffer(i:i+1) = ', '
            i = i + 2
         end do

         write (lunout, '(a)') buffer(1:i-3)
         write (lunout, '(a)') 'DATAPACKING=BLOCK'

         do i = 1, ndim
            write (lunout, '(1p, 10e13.5)') x_target(i,:)
         end do
         do n = 1, nf
            write (lunout, '(1p, 10e13.5)') f_target(n,:)
         end do
         do n = 1, nf
            do i = 1, ndim
               write (lunout, '(1p, 10e13.5)') f1(n,:,i)
            end do
         end do
         do n = 1, nf
            do i = 1, ndim
               write (lunout, '(1p, 10e13.5)') f2(n,:,i)
            end do
         end do
         do n = 1, nf
            do i = 1, ndim
               write (lunout, '(1p, 10e13.5)') fk(n,:,i)
            end do
         end do

         close (lunout)

      end if

      end subroutine save_results

   end program anchor
