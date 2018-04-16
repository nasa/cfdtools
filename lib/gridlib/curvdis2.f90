!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine curvdis2 (ndat, xdat, ydat, n, power, ismooth, lunout, &
                        arcs_only, x, y, ier)
!
!  CURVDIS2 has the same calling sequence as CURVDIS, but it applies geometric
!  normalization to the data first, then denormalizes results before returning.
!  This means the coordinate with the larger data range is transformed to [0,1],
!  and the other coordinate is scaled by the same factor to preserve shape,
!  along with being shifted to the origin, before curvatures are calculated.
!  The purpose is to avoid extreme curvature values from tiny geometries that
!  can give undesirable curvature-based point distributions.
!
!  Later extension: In order to test the behavior of redistributing arc lengths
!  of a nongeometric curve treated much as it appears in an x/y plot, an option
!  to perform nongeometric normalization has been retrofitted in a way that
!  does not affect existing applications:  Enter ier = -99 to invoke this
!  option.
!
!  Note that the data points are NOT normalized/denormalized in place.
!
!  11/08/12  D.A.Saunders  Application of existing utilities for use by
!                          CAPSULE_GRID.
!  10/23/13    "      "    In case of CURVDIS failure, reduce the exponent used
!                          by 0.1 until it converges (to uniform spacing in the
!                          worst case of reducing it to 0.).
!  11/25/13    "      "    The arcs_only case should just unscale the arc
!                          lengths (no shift from starting at 0.).
!
!  Author:  David Saunders, ERC, Inc. at NASA Ames Research Center, CA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)  :: ndat       ! Number of input curve points
   real,    intent (in)  :: xdat(ndat) ! X coordinates of input curve
   real,    intent (in)  :: ydat(ndat) ! Y coordinates of input curve
   integer, intent (in)  :: n          ! Number of output curve points
   real,    intent (in)  :: power      ! Exponent for spacing function;
                                       ! controls the point clustering.
                                       ! POWER must be in [0., 1.];
                                       ! 0.0 produces uniform spacing;
                                       ! 1.0 maximizes the curvature effect;
                                       ! 0.5 is suggested
   integer, intent (in)  :: ismooth    ! 0 => no smoothing;
                                       ! 1 => smooth shape function only;
                                       ! 2 => smooth redistributed arcs only;
                                       ! 3 => perform both sets of smoothings
   integer, intent (in)  :: lunout     ! Logical unit for showing convergence
                                       ! history; LUNOUT < 0 suppresses it
   logical, intent (in)  :: arcs_only  ! T => skip calculations of X & Y(1:N)
                                       ! but return revised arcs as X(1:N)
   real,    intent (out) :: x(n), y(n) ! X & Y coordinates of output curve,
                                       ! unless ARCS_ONLY = T, in which case
                                       ! just the arc-lengths are returned as
                                       ! X(1:N)
   integer, intent (inout) :: ier      ! On input:
                                       !  -99 => apply nongeometric normalizn.
                                       !         as explained above;
                                       !    anything else => geometric nrmlzn.
                                       ! On output:
                                       !    0 => no errors;
                                       !    1 => POWER is out of range;
                                       !    2 => failure in ARBDIS
!  Local constants:

   integer, parameter :: ndim = 2  ! Any "z" argument is ignored by get/usescale
   real,    parameter :: zero = 0.

!  Local variables:

   integer :: ier_local
   real    :: xdata(ndat), ydata(ndat)
   real    :: pwr, scale(ndim), shift(ndim)
   character (1) :: nrm

!  Execution:

   xdata(:) = xdat(:)
   ydata(:) = ydat(:)

   if (ier == -99) then
       nrm = 'N'
   else
       nrm = 'G'
   end if

   call getscale (nrm, ndim, ndat, xdata, ydata, ydata, scale, shift, ier_local)
   call usescale ('N', ndim, ndat, xdata, ydata, ydata, scale, shift, ier_local)

   pwr = power
   do  ! Until a low-enough exponent allows convergence

      call curvdis (ndat, xdata, ydata, n, pwr, ismooth, lunout, arcs_only, &
                    x, y, ier)
      if (ier == 0) exit

      write (*, '(a, i3, f5.1)') &
         'CURVDIS2 trouble.  ier, power:', ier, pwr
      if (pwr == zero) then  ! Should never fail with zero
         write (*, '(a)') 'Aborting.'
         stop
      end if
      pwr = max (zero, pwr - 0.1)
      write (*, '(a, f5.1)') 'Retrying with lower power:', pwr
   end do

   if (arcs_only) then  ! x(:) is really arc lengths starting at 0.
      y(:) = x(:)       ! Avoid undefined y(:) when usescale denormalizes
      shift(:) = zero
   end if

   call usescale ('D', ndim, n, x, y, ydata, scale, shift, ier_local)

   end subroutine curvdis2
