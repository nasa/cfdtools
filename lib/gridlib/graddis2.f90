!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine graddis2 (ndat, xdat, ydat, sdat, fdat, n, power, &
                        ismooth, lunout, x, y, s, ier)
!
!  GRADDIS2 is a gradient-based analogue of curvature-based utility CURVDIS2.
!  It needs a profile of some function along with the 2-space geometric curve
!  and its cumulative arc (chord) lengths.
!
!  The gradient of |f| is taken to be |df/ds| were s is arc length from point 1.
!  Three-point derivatives are used (1-sided at the end points).
!
!  GRADDIS2 applies geometric normalization to the coordinates and normalizes
!  the function data before calling GRADDIS in a loop over the power argument
!  until convergence is achieved, then denormalizes results before returning.
!
!  This means the coordinate with the larger data range is transformed to [0,1],
!  and the other coordinate is scaled by the same factor to preserve shape,
!  along with being shifted to the origin.  The function data are also shifted
!  and normalized to [0,1] before gradients are calculated.  Larger gradient
!  magnitudes will lead to more point clustering just as for curvature, so the
!  same technique is employed at the lower level on different data.
!
!  Note that the data points are NOT normalized/denormalized in place.
!
!  10/23/13  D.A.Saunders  CURVDIS2 starting point.
!  11/20/13    "      "    Initial adaptation of GRADDIS2 from CURVDIS2.
!  11/22/13    "      "    The advent of LCSFIT2 prompted NOT normalizing and
!                          denormalizing in place: use automatic arrays instead.
!  12/03/13    "      "    GRADDIS now returns the updated arc lengths, not just
!                          the updated coordinates, because recovering them from
!                          the latter (for function interpolations) is not
!                          identical to using the arcs that produced those new
!                          coordinates.
!
!  Author:  David Saunders, ERC, Inc. at NASA Ames Research Center, CA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)  :: ndat        ! Number of input curve points
   real,    intent (in)  :: xdat(ndat)  ! X coordinates of input curve
   real,    intent (in)  :: ydat(ndat)  ! Y coordinates of input curve
   real,    intent (in)  :: sdat(ndat)  ! Corresponding arc lengths, probably
                                        ! needed for subsequent fn. interplns.
   real,    intent (in)  :: fdat(ndat)  ! Some function such as temperature,
                                        ! used for a df/ds-based shape function
   integer, intent (in)  :: n           ! Number of output curve points
   real,    intent (in)  :: power       ! Exponent for spacing function;
                                        ! controls the point clustering.
                                        ! POWER is in [0., 1.] (or higher);
                                        ! 0.0 produces uniform spacing;
                                        ! 1.0 maximizes the curvature effect;
                                        ! 1.0 is suggested
   integer, intent (in)  :: ismooth     ! 0 => no smoothing;
                                        ! 1 => smooth shape function only;
                                        ! 2 => smooth redistributed arcs only;
                                        ! 3 => perform both sets of smoothings;
                                        ! 1 is suggested
   integer, intent (in)  :: lunout      ! Logical unit for showing convergence
                                        ! history; LUNOUT < 0 suppresses it
   real,    intent (out) :: x(n), y(n)  ! X & Y coordinates of output curve
   real,    intent (out) :: s(n)        ! Corresponding new arc lengths
   integer, intent (out) :: ier         ! 0 => no errors;
                                        ! 1 => POWER is outside [0., 2.];
                                        ! 2 => failure in ARBDIS
!  Local constants:

   integer, parameter :: ndim = 2  ! Any "z" argument is ignored by get/usescale
   real,    parameter :: zero = 0.

!  Local variables:

   integer :: ier_loc
   real    :: pwr, xscale(ndim), xshift(ndim)
   real    :: fscale, fshift, sshift, total
   real    :: xdata(ndat), ydata(ndat), &     ! Avoided touching inputs
              sdata(ndat), fdata(ndat), &
              fpdat(ndat), fpp(ndat)          ! Use fpp for unneeded fk,
                                              ! not fp for unneeded fpp
!  Procedures:

   external :: getscale  ! Get normalizing shift/scale factors
   external :: usescale  ! Apply normalizing/denormalizing factors
   external :: fd12k     ! Finite difference derivatives [& curvature, unused]
   external :: graddis   ! Redistribute points according to |df/ds|

!  Execution:

   xdata(:) = xdat(:)
   ydata(:) = ydat(:)

   call getscale ('G', ndim, ndat, xdata, ydata, ydata, xscale, xshift, ier_loc)
   call usescale ('N', ndim, ndat, xdata, ydata, ydata, xscale, xshift, ier_loc)

   sdata(:) = sdat(:)
   sshift   = zero  ! Arc lengths already start at zero; just scale arcs

   call usescale ('N', 1, ndat, sdata, sdata, sdata, xscale, sshift, ier_loc)

   fdata(:) = fdat(:)

   call getscale ('N', 1, ndat, fdata, fdata, fdata, fscale, fshift, ier_loc)
   call usescale ('N', 1, ndat, fdata, fdata, fdata, fscale, fshift, ier_loc)

   call fd12k (ndat, sdata, fdata, fpdat, fpp, fpp) ! fpdat = f(norm) gradient

   fpdat(:) = abs (fpdat(:))

   pwr = power
   do  ! Until a low-enough exponent allows convergence

      call graddis (ndat, xdata, ydata, sdata, fpdat, n, pwr, &
                    ismooth, lunout, x, y, s, ier)
      if (ier == 0) exit

      write (*, '(a, i3, f5.1)') &
         'GRADDIS2: Trouble in GRADDIS.  ier, power:', ier, pwr
      if (pwr == zero) then  ! Should never fail with zero
         write (*, '(a)') 'Aborting.'
         stop
      end if
      pwr = max (zero, pwr - 0.1)
      write (*, '(a, f5.1)') 'Retrying with lower power:', pwr
   end do

   call usescale ('D', ndim, n, x, y, y, xscale, xshift, ier_loc)
   call usescale ('D', 1,    n, s, s, s, xscale, sshift, ier_loc)

   end subroutine graddis2
