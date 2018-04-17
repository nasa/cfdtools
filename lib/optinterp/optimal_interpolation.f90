!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   module optimal_interpolation

!  Fortran 90 module for n-dimensional optimal interpolation.
!  Released under the BSD 2-Clause License

!  Copyright (c) 2005, Alexander Barth <a.barth@ulg.ac.be>,
!  <barth.alexander@gmail.com> 
!  All rights reserved.
!
!  Redistribution and use in source and binary forms, with or without
!  modification, are permitted provided that the following conditions are
!  met:
!
!  1. Redistributions of source code must retain the above copyright
!  notice, this list of conditions and the following disclaimer.
!
!  2. Redistributions in binary form must reproduce the above copyright
!  notice, this list of conditions and the following disclaimer in the
!  documentation and/or other materials provided with the distribution.
!
!  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
!  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
!  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
!  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
!  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
!  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
!  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!  Author:       Alexander Barth
!  Dependencies: LAPACK (dsyev and ancillary routines) *** No longer needed. ***
!
!  Extensions made at NASA Ames Research Center by David Saunders, ELORET Corp.:
!
!     This package provides a powerful means of treating unstructured data.  It
!     belongs in the family of methods best known as Kriging, which determine
!     the linear combination of neighboring data points that minimizes the
!     interpolation error in some sense under certain statistical assumptions.
!
!     The method can be used for smoothing noisy data points as well as for
!     interpolating the points precisely under reasonable circumstances.  Ill-
!     conditioned cases are handled by returning shortest-length solutions
!     according to a tolerance on the size of diagonal elements in the matrix
!     eigenvalue decomposition, and these solutions produce some smoothing.
!
!     (Later:  The 01/29/07 fudge of the diagonal elements should suffice to
!     protect the Cholesky method now adopted.)
!
!     The initial application at ARC is for interpolating unstructured high
!     fidelity data corrections on to lower-fidelity aerodynamic data, and for
!     showing where further high-fidelity data points are desirable.
!
!     This version normalizes and denormalizes the data in-place, so save the
!     inputs at the higher level if round-off differences are significant.
!     (Actually, it no longer normalizes the function data - only their
!     error variances, now one per function for all observations.)
!
!     Does it matter that each optimization uses just a subset of the
!     observations, yet the normalization uses all the observations?
!
!     Provision for more than just the original Gaussian correlation function
!     has been made, and a partial answer to using one correlation length (per
!     dimension) for what might be widely-varying spacing in the data has been
!     incorporated.  This allows the correlation lengths to be rescaled (in the
!     same proportion) for each interpolation in order to force the furthest
!     neighbor to make a specified contribution to the right-hand-side.  This
!     option can also be suppressed.  It was prompted by trying to deal with
!     aerodynamic coefficients that are insensitive to high Mach numbers (where
!     the data can therefore be more sparse) while also coping (in the same
!     dataset) with transonic Mach numbers where much higher resolution is
!     needed in the data points.
!
!  11/14/06  David Saunders  Extended for multiple functions per observation.
!                            All functions use the same error variance at
!                            all observations.
!  11/15/06    "      "      Dispensed with (wp) = 4 or 8.  Use a compiler
!                            switch to give the 64-bit arithmetic that is
!                            highly recommended.  Switched to free format
!                            (.f90 extension) and made cosmetic changes.
!                            Introduced normalization and denormalization.
!  11/17/06    "      "      Avoid forming the pseudo-inverse explicitly.
!                            This speeds the method by about 30% if m = 30.
!  01/23/07    "      "      Guard against zero standard deviations from flat
!                            functions, since we normalize by them.
!  01/29/07    "      "      Guard against singularity by adding (10 + m)eps to
!                            the diagonal elements.
!  01/30/07    "      "      Adjust the correlation lengths for each target pt.
!                            so that the mth nearest observation contributes a
!                            finite amount.  This should cope with observations
!                            that are widely differing distances apart.  It
!                            represents a crude attempt at optimizing what are
!                            largely unknown quantities (correlation lengths).
!  02/02/07    "      "      Don't adjust if the mth neighbor is closer than the
!                            tolerance.
!  02/12/07    "      "      Provide a choice of correlation functions, and for
!                            suppressing the adaptive correlation length option.
!  02/14/07    "      "      Anisotropic correlations aren't working well.
!                            Arrange for evaluating them on a [0, 2] grid in
!                            each (normalized) coordinate direction if the first
!                            correlation length is negative.  Its absolute value
!                            will be used.
!  05/09/07    "      "      Experience suggests that the systems solved are
!                            invariably positive definite.  Use of the Cholesky
!                            solver from LINPACK cleans things up and speeds
!                            the method by another ~30% for m = 30 neighbors.
!  11/14/08    "      "      Normalizing the function data is actually redundant
!                            but calculating their variances and normalizing
!                            them is not.  The error variances at the target
!                            points (and the observation variances) still need
!                            to be denormalized.  Allowing for different
!                            variances for each function meant gvar(:) had to
!                            become gvar(:,:) in subroutine optiminterp.
!                            Also, allowing for different error variances at
!                            each observation appears impractical.  Therefore,
!                            change that to the same error variance at all
!                            observation points for a given function. Since the
!                            different function error variances are normalized
!                            by different function data variances, we cannot get
!                            away with the same left hand side for all functions
!                            at each interpolation: we have to factorize the
!                            LHS matrix nf times per target point.
!  04/13/16    "      "      Application to a radio blackout problem prompted a
!                            kludged option to use inverse-distance weighting
!                            that may benefit development of other applications.
!                            (Testing on function data that are all identical,
!                            for instance, should produce identical interpolated
!                            function values, but this seems not to be true for
!                            the optimal interpolations.)  The "adapt" argument
!                            to optiminterp, if negative, invokes this option.
!                            Installing the option here allows reuse of the
!                            search for m nearest neighbors and of the scaling/
!                            unscaling that has been automated here.  A new
!                            argument, verbose, has also been added to the
!                            subroutine optiminterp to help understand what is
!                            going on with diagnostic output to unit 50.
!  05/09/16    "      "      The normalized coordinate printout (verbose mode)
!                            now includes the function data as well.  Tecplot
!                            can triangulate these coordinates and make contour
!                            plots.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

   private

!  Constants used by the package:

   real, parameter :: one = 1., zero = 0.

!  Callable routine for interpolating one dataset to one grid:

   public :: optiminterp

   contains

!     --------------------------------------------------------------------------

      subroutine mean_and_sd (ne, np, ie, x, xmean, xsd)

!     Calculate the mean and standard deviation for element ie of packed array
!     x(1:ne,1:np) using a single pass through the data.

!     11/15/06  D. Saunders  Adaptation of earlier work.

      implicit none

!     Arguments:

      integer, intent (in)  :: ne, np    ! # elements per point & # points
      integer, intent (in)  :: ie        ! Element x/y/.. or f1/f2/.. to average
      real,    intent (in)  :: x(ne,np)  ! Packed data
      real,    intent (out) :: xmean     ! Mean of indicated quantity
      real,    intent (out) :: xsd       ! Standard deviation of that quantity

!     Local variables:

      integer :: i, j
      real    :: ssq, term, xbar

!     Execution:

      xbar = zero;  ssq = zero;  i = ie

      do j = 1, np
         term = (x(i,j) - xbar) / real (j)            ! Permits just 1 divide
         xbar = xbar + term
         ssq  = ssq  + (term * real (j)) * (term * real (j-1))   ! May be better
      end do                                          ! than term**2*r(j)*r(j-1)

      xmean = xbar;    xsd = sqrt (ssq / real (np))
      if (xsd < 1.e-6) xsd = one                      ! Since we normalize by it

      end subroutine mean_and_sd

!     --------------------------------------------------------------------------

      subroutine norm_denorm (mode, ne, np, ie, xmean, xsd, x)

!     Normalize or denormalize the indicated data, in-place.

      implicit none

!     Arguments:

      integer, intent (in)    :: mode     ! 1 => normalize; 2 => denormalize
      integer, intent (in)    :: ne, np   ! # elements per point & # points
      integer, intent (in)    :: ie       ! Element x/y/.. or f1/f2/.. to adjust
      real,    intent (in)    :: xmean    ! Mean of indicated quantity
      real,    intent (in)    :: xsd      ! Standard deviation of that quantity
      real,    intent (inout) :: x(ne,np) ! Packed data

!     Local variable:

      real :: rsd

!     Execution:

      select case (mode)

      case (1)

         rsd = one / xsd;  x(ie,:) = (x(ie,:) - xmean) * rsd

      case (2)

         x(ie,:) = x(ie,:) * xsd + xmean

      end select

      end subroutine norm_denorm

!     --------------------------------------------------------------------------

      subroutine select_nearest (m, x, ox, param, adapt, fun, lunpr, indices)

!     Select the m observations from ox(1:nd,1:on) closest to point x(1:nd).

!     01/30/07  D. Saunders  Adjust the correlation lengths (via param(:)) to
!                            ensure that the mth closest point contributes some
!                            minimum amount (hard-coded).
!     02/12/07   "      "    The adaptive option gets awkward with more than
!                            one correlation function (isotropic or not).

      implicit none

!     Arguments:

      integer, intent (in)    :: m
      real,    intent (in)    :: x(:), ox(:,:)
      real,    intent (inout) :: param(:)    ! May be rescaled if adapt > 0.
      real,    intent (in)    :: adapt       ! 0. <= adapt << 1.;
                                             ! < 0. invokes an option not to
                                             ! use optimal interpolation in
                                             ! favor of using inverse-distance
                                             ! weighting
      character, intent (in)  :: fun(:)*1    ! 'G' = Gaussian correlation;
                                             ! 'L' = linear;
                                             ! one per dimension;
                                             ! lower case => isotropic
      integer, intent (in)    :: lunpr       ! > 0 turns on diagnostics to lunpr
      integer, intent (out)   :: indices(m)

!     Local variables:

      integer :: i, mth
      integer :: max_pannier(1)   ! See sort utility
      logical :: verbose
      real    :: d(size (ox, 2)), distance_sq(m), dmsq, dsq_bound, fmin

!     Execution:

      verbose = lunpr > 0

!     Calculate a measure of (squared) normalized distance to each observation:

      do i = 1, size (ox, 2)
!!!     d(i) = sum (((x(:) - ox(:,i)) * param(:))**2)  ! This denormalizes!
        d(i) = sum (((x(:) - ox(:,i)))**2)
      end do

      call sort (d, m, indices)

      distance_sq(:) = d(indices)

      if (verbose) then
         write (lunpr, '(a, (9es14.7))') ' dsq(1:ndata):', d
         write (lunpr, '(a, (5es16.8))') ' param0:', param
         write (lunpr, '(a, (40i6))')    ' indices: ', indices
         write (lunpr, '(a, (9es14.7))') ' initial dnorm_sqrd(*):', distance_sq
      end if

!     Rescale cor. lengths so that the mth closest point is not worthless?

      if (adapt > zero) then ! Assume adapt << 1.

         fmin = adapt
         max_pannier = maxloc (distance_sq);  mth = max_pannier(1)
         dmsq = distance_sq(mth)

         select case (fun(1))

            case ('g') ! Gaussian, isotropic

               dsq_bound = -log (fmin)
               if (dmsq > dsq_bound) then
                  param = param * sqrt (dsq_bound / dmsq)
                  if (verbose) write (lunpr, *) 'scaled param:', param
               end if

            case ('l') ! Linear, isotropic

               dsq_bound = (one - fmin)**2
               if (dmsq > dsq_bound) then
                  param = param * sqrt (dsq_bound / dmsq)
                  if (verbose) write (lunpr, *) 'scaled param:', param
               end if

            case ('G', 'L')

               do i = 1, size (x)

                  select case (fun(i))

                     case ('G') ! Gaussian, anisotropic

                        dsq_bound = -log (fmin)
                        if (distance_sq(i) > dsq_bound) then
                           param(i) = param(i) * sqrt (dsq_bound/distance_sq(i))
                           if (verbose) &
                              write (lunpr, *) 'scaled param i:', i, param(i)
                        end if

                     case ('L') ! Linear, anisotropic

                        dsq_bound = (one - fmin)**2
                        if (distance_sq(i) > dsq_bound) then
                           param(i) = param(i) * sqrt (dsq_bound/distance_sq(i))
                           if (verbose) &
                              write (lunpr, *) 'scaled param i:', i, param(i)
                        end if

                  end select

               end do

         end select

      end if

      end subroutine select_nearest

!     --------------------------------------------------------------------------

      subroutine sort (d, m, pannier)

!     Return the indices of the m smallest elements in d(:).
!     The algorithm is succinctly coded, but would a heap sort be faster?

      implicit none

!     Arguments:

      real,    intent (in)  :: d(:)
      integer, intent (in)  :: m
      integer, intent (out) :: pannier(m)  ! Pannier = French for basket

!     Local variables:

      integer :: i, max_pannier(1)

!     Execution:

      do i = 1, m
         pannier(i) = i
      end do

      max_pannier = maxloc (d(pannier))   ! maxloc returns a vector

      do i = m + 1, size (d)
         if (d(i) < d(pannier(max_pannier(1)))) then
            pannier(max_pannier(1)) = i
            max_pannier = maxloc (d(pannier))
         end if
      end do

      end subroutine sort

!     --------------------------------------------------------------------------

      subroutine observation_covariance (ovar, indices, R)

      implicit none

      real,    intent (in)  :: ovar(:)          ! Observation variances
      integer, intent (in)  :: indices(:)       ! List of m indices into ovar(:)
      real,    intent (out) :: R(size(indices)) ! Originally m x m, now just m

      integer :: i

!!!   R = 0.  No need to store these or use them

      do i = 1, size (indices)
         R(i) = ovar(indices(i))
      end do

      end subroutine observation_covariance

!     --------------------------------------------------------------------------

      function background_covariance (x1, x2, param, fun) result (c)

!     Original               Only Gaussian correlation was available.
!     02/09/07  D. Saunders  Allowed for anisotropic correlation.
!
!     All the correlation functions should be 1 at distance d = 0.
!     The correlation lengths determine where the functions drop to zero (or
!     how quickly they asymptote to zero).
!     Since the Cartesian product of two 1-D functions is in general NOT the
!     same as using an isotropic function of distance, the higher level should
!     check for a request to use the same kind of correlation in all dimensions
!     and adjust the function type accordingly.

      implicit none

      real, intent (in)      :: x1(:), x2(:) ! Two pts. distance d apart
      real, intent (in)      :: param(:)     ! 1/cor_lengths before normalizing
      character, intent (in) :: fun(:)*1     ! 'G' = Gaussian; 'L' = linear;
                                             ! one per dimension; if all are
                                             ! dimensions are given the same
                                             ! function, use 'g' or 'l' instead
                                             ! to work with total distance
                                             ! instead of distance components
      real                   :: c            ! Returned result

!     Local variables:

      integer :: i
      real    :: d(size (x1))

!     Execution:

      d = (x1 - x2)*param

      select case (fun(1))

         case ('g') ! Gaussian, isotropic

            c =  exp (-sum (d**2))

         case ('l') ! Linear, isotropic

            c = max (one - sqrt (sum (d**2)), zero)

         case ('G', 'L')

            c = one

            do i = 1, size (x1)

               select case (fun(i))

                  case ('G') ! Gaussian, anisotropic

                     c = c * exp (-d(i)**2)

                  case ('L') ! Linear, anisotropic

                     c = c * max (one - abs (d(i)), zero)

               end select

            end do

      end select

      end function background_covariance

!     --------------------------------------------------------------------------

      subroutine plot_cor_function (nd, param, fun)

!     Evaluate the correlation function(s) on [0,2] in each coordinate direction
!     for plotting purposes (prompted by anisotropic cases).
!     Initially, the nd = 2 case (only) is coded.  Add others as needed.

      implicit none

!     Arguments:

      integer,   intent (in) :: nd       ! # dimensions
      real,      intent (in) :: param(:) ! Cor. length inverses, normalized
      character, intent (in) :: fun(:)*1 ! 'G' = Gaussian correlation;
                                         ! 'L' = linear;
                                         ! one per dimension;
                                         ! lower case => isotropic
!     Local constants:

      integer, parameter :: lunplot = 49 ! Sorry to hard-code it
      integer, parameter :: neval   = 65 ! Adjust to suit

!     Local variables:

      integer :: i, j, ieval, ntotal
      real    :: dx, x1(nd), x2(nd)
      real, allocatable :: feval(:), xeval(:,:)

!     Execution:

      open (lunplot, file='correlation_function.dat', status='unknown')

      x1(:)  = zero;  ieval = 1
      dx     = (one + one) / real (neval - 1)
      ntotal = neval**nd

      allocate (feval(ntotal))

      select case (nd)

         case (2)

         allocate (xeval(nd,ntotal))

         do j = 1, neval
            do i = 1, neval
               xeval(1,ieval) = dx * real (i - 1)
               xeval(2,ieval) = dx * real (j - 1)
               feval(ieval)   = background_covariance (x1, xeval(:,ieval),     &
                                param, fun)
               ieval = ieval + 1
            end do
         end do

         write (lunplot, '(a)') &
            'TITLE = "Correlation function(s)"', &
            'VARIABLES = "X1_norm", "X2_norm", "Correlation function"', &
            'ZONE T = "Correlation function(s)"', &
            'DATAPACKING=POINT'
         write (lunplot, '(2(a, i2))') 'I=', neval, ', J=', neval
         write (lunplot, '(3es13.5)') (xeval(:,i), feval(i), i = 1, ntotal)

         case default

      end select

      deallocate (xeval, feval)

      close (lunplot)

      end subroutine plot_cor_function

!     --------------------------------------------------------------------------

      subroutine optiminterp (ox, of, ovar, cor_len, adapt, fun, m, gx,        &
                              verbose, gf, gvar)

!     Main optimal interpolation routine, with automated scaling of the data
!     in-place.  (Hence the (inout) intents.)
!
!     Original:              Gaussian correlation function, isotropic (only).
!     02/09/07  D. Saunders  Added correlation function choices.
!     02/12/07   "      "    Option to adapt correlation lengths to sparse data.
!     02/14/07   "      "    Option to evaluate correlation functions for plots.

      implicit none

!     Arguments:

      real,    intent (inout) :: ox(:,:), of(:,:) ! Observations; dimensions are
                                                  ! (1:nd,1:np) and (1:nf,1:np)
      real,    intent (inout) :: ovar(:)          ! Observation error variances,
                                                  ! (1:nf) (originally 1:np)
      real,    intent (inout) :: cor_len(:)       ! Correlation lengths, (1:nd);
                                                  ! cor_len(1) < 0 produces a
                                                  ! plottable dataset of the
                                                  ! correlation function(s)
      real,    intent (in)    :: adapt            ! Minimum contribution of each
                                                  ! correlation function for the
                                                  ! mth nearest neighbor to
                                                  ! contribute;  0. suppresses
                                                  ! this option; 0. < adapt < 1.
                                                  ! raises correlation length(s)
                                                  ! if necessary for the current
                                                  ! target point; < 0. invokes
                                                  ! inverse-distance weighting
      character (1), intent (in)  :: fun(:)       ! 'G' = Gaussian correlation;
                                                  ! 'L' = linear;
                                                  ! one per dimension
      integer, intent (in)    :: m                ! # nearest observations used
      real,    intent (inout) :: gx(:,:)          ! Target grid coords.  (nd,gn)
      logical, intent (in)    :: verbose          ! For diagnostics to unit 50
      real,    intent (out)   :: gf(:,:),       & ! Interpolated func(s) (nf,gn)
                                 gvar(:,:)        ! and error variances  (nf,gn)

!     Local constants:

      integer, parameter :: lunpr = 50  ! For diagnostic output via verbose arg.

!     Local variables:

      integer  :: i, info, ipdone, j, j1, j2, k, lun, n, nd, nf, np
      integer  :: gn, indices(m)
      real     :: A(m,m), Asave(m,m), D(m), R(m), PH(m), RHSsave(m), W(m)
      real     :: fn(m),  param(size (ox, 1)), param0(size (ox, 1))
      real     :: oxmean(size (ox, 1)), oxsd(size (ox, 1))
      real     :: ofmean(size (of, 1)), ofsd(size (of, 1))
      real     :: ovar_norm(size (of, 1))
      real     :: add_to_diagonal, fudge

      real     :: dsq, dx(size (ox, 1)), sumW

      logical   :: isotropic, show_cor_function
      character :: func(size (ox, 1)) * 1

      logical, save :: first_call = .true.

!     Execution:

      if (verbose) then
         lun = lunpr
      else
         lun = 0
      end if

      fudge = real (10 + m) * epsilon (fudge)  ! See DACE Kriging Toolbox guide
      show_cor_function = cor_len(1) < zero
      if (show_cor_function) cor_len(1) = -cor_len(1)

      nd = size (ox, 1)  ! # coordinates
      np = size (ox, 2)  ! # observation points
      nf = size (of, 1)  ! # functions at each observation point
      gn = size (gx, 2)  ! # target (packed grid) points to interpolate at

      func = fun;  isotropic = .true.

!     Normalize the observation coordinates in-place:

      do j = 1, nd
         if (func(j) /= func(1)) isotropic = .false.
         call mean_and_sd (nd, np, j, ox, oxmean(j), oxsd(j))
         call norm_denorm ( 1, nd, np, j, oxmean(j), oxsd(j), ox)
      end do

      if (verbose) then  ! Write normalized coords. as columns for plotting
         if (first_call) then
            first_call = .false.
            write (lunpr, '(/, a)') &
               ' Normalized data coordinates 1:ndim and function data:'
            do j = 1, np
               write (lunpr, '(20es13.5)') ox(:,j), of(:,j)
            end do
         end if
      end if

      if (isotropic) then
         select case (func(1))
            case ('G')
               func(:) = 'g'
            case ('L')
               func(:) = 'l'
            case default
               func(:) = 'g'
         end select
      end if

      param0(:) = oxsd(:) / cor_len(:)  ! 1. / cor_len prior to normalizing
                                        ! May be rescaled for each target point

!     Normalize the observation function data in-place. [No: variances only.]

      do j = 1, nf
         call mean_and_sd (nf, np, j, of, ofmean(j), ofsd(j))
!!!      call norm_denorm ( 1, nf, np, j, ofmean(j), ofsd(j), of)
      end do

      ovar_norm(:) = ovar(:) / ofsd(:)**2

      if (verbose) then
         write (lunpr, '(/, a, 3es13.5)') ' oxmean:   ', oxmean
         write (lunpr, '(   a, 3es13.5)') ' oxsd:     ', oxsd
         write (lunpr, '(/, a,12es13.5)') ' ofmean:   ', ofmean
         write (lunpr, '(   a,12es13.5)') ' ofsd:     ', ofsd
         write (lunpr, '(/, a,12es13.5)') ' ovar in:  ', ovar
         write (lunpr, '(   a,12es13.5)') ' ovar norm:', ovar_norm
      end if

!     Normalize the target coordinates in-place by the observation mean & sd.

      do j = 1, nd
         call norm_denorm (1, nd, gn, j, oxmean(j), oxsd(j), gx)
      end do

      if (show_cor_function) then
          call plot_cor_function (nd, param0, func)
      end if

!     For each target grid point:

      do i = 1, gn

         if (verbose) then
            write (lunpr, '(/, a, 6es13.5)') ' normalized target:', gx(:,i)
         end if

         param = param0  ! Possibly rescaled for each target point

!        Determine the indices of the m nearest observations.
!        If specified, adjust correlation lengths so that the furthest neighbor
!        contributes some minimum amount to the RHS vector in normalized space.

         call select_nearest (m, gx(:,i), ox, param, adapt, func, lun, indices)

         if (adapt < zero) then  ! Use simple inverse distance weighting:

            do k = 1, m
               dx(:) = gx(:,i) - ox(:, indices(k))
               dsq = dot_product (dx, dx)
               if (abs (dsq) < 1.e-6) then
                  W(:) = 0.; W(k) = 1.
                  exit
               end if
               W(k) = 1.0 / sqrt (dsq)
            end do
            sumW = sum (W)
            if (verbose) then
               write (lunpr, '(9x, a)') 'Target normlzd. V&h,  neighbor V&h,  W'
               write (lunpr, '(7x, 2es16.8, 2x, 2es16.8, es20.8)') &
                  (gx(:,i), ox(:, indices(k)), W(k)/sumW, k = 1, m)
            end if

            do n = 1, nf
               do k = 1, m
                  fn(k) = of(n,indices(k))
               end do
               gf(n,i) = dot_product (W, fn) / sumW
               if (verbose) write (lunpr, '(a, 2i3, a, es13.5)') &
                  ' finterp', n, i, ':', gf(n,i)
               gvar(n,i) = zero
            end do

            cycle
         end if

!        Form the error covariance matrix background field and the RHS vector:

         do j2 = 1, m
            do j1 = 1, j2  ! Upper triangle only
               A(j1,j2) = &
                  background_covariance (ox(:,indices(j1)), ox(:,indices(j2)), &
                                         param, func)
            end do
            PH(j2) = background_covariance (gx(:,i), ox(:,indices(j2)),        &
                                            param, func)
         end do

!        Allowing for different observation error variances per function means
!        a different LHS for each function.  Originally, these error variances
!        were assumed to be the same for each function but (possibly) different
!        at different observation points. That has been abandoned now at a cost.

         if (nf > 1) then
            Asave = A;  RHSsave = PH
         end if

!        For each function at this target point:

         do n = 1, nf

!           Form the error covariance matrix of the observations:

!!!         call observation_covariance (ovar, indices, R)  ! See comment above

!           Covariance matrix of the "innovation":

!!!         add_to_diagonal = ovar_norm(n) + fudge  ! Guard against singularity
            add_to_diagonal = fudge

            do j = 1, m
!!!            A(j,j) = A(j,j) + R(j) + fudge
               A(j,j) = A(j,j) + add_to_diagonal
            end do

            if (verbose) then
               do j1 = 1, m
                  write (lunpr, '(es13.5, 2x, 20es13.5)') PH(j1), A(j1,j1:m)
               end do
            end if

!           Cholesky factorization of A:

            call cholesky_factorization (m, A, info)

            if (info > 0) then  ! Should never happen
               write (*, *) ' Non-positive definite leading minor ', info
               write (*, *) ' Target point and function #: ', i, n
               stop
            end if

!           Complete the solution by Cholesky decomposition.
!           The RHS is overwritten.

            W(:) = PH(:)

            call cholesky_solution (m, A, W)

!           "Compute the analysis" (interpolate this function):

            do k = 1, m
               fn(k) = of(n,indices(k))
            end do

            gf(n,i) = dot_product (W, fn)

            if (verbose) then
               write (lunpr, '(a, (12es13.5))') ' W:', W
               write (lunpr, '(a, 2i3, a, es13.5)') &
                  ' finterp', n, i, ':', gf(n,i)
            end if

!           Compute the error variance of the analysis for this function:

            gvar(n,i) = one - dot_product (W, PH)

            if (nf > 1) then
               A = Asave;  PH = RHSsave
            end if

         end do  ! Next function

      end do  ! Next target point

!     Denormalize certain inputs and outputs in-place:

      do j = 1, nd  ! Observation coordinates
         call norm_denorm (2, nd, np, j, oxmean(j), oxsd(j), ox)
      end do

      do j = 1, nf  ! Observation functions [no need]
!!!      call norm_denorm (2, nf, np, j, ofmean(j), ofsd(j), of)
      end do

      do j = 1, nd  ! Target coordinates
         call norm_denorm (2, nd, gn, j, oxmean(j), oxsd(j), gx)
      end do

      do j = 1, nf  ! Interpolated functions [no need]
!!!      call norm_denorm (2, nf, gn, j, ofmean(j), ofsd(j), gf)
      end do

      do i = 1, gn  ! Interpolated error variances
         gvar(:,i) = gvar(:,i) * ofsd(:)**2
      end do

      if (show_cor_function) cor_len(1) = -cor_len(1)

      end subroutine optiminterp

!     --------------------------------------------------------------------------

   end module optimal_interpolation
