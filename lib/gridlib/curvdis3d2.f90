!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine curvdis3d2 (ndat, xdat, ydat, zdat, n, power, ismooth, lunout, &
                          arcs_only, x, y, z, ier)
!
!  CURVDIS3D2 is the 3-space analogue of CURVDIS2.  It has the same calling
!  sequence as CURVDIS3D, but it applies geometric normalization to the data
!  first, then denormalizes results before returning.  It also iterates until
!  the lower-level routine converges successfully, by lowering power by 0.1 as
!  often as necessary.
!
!  The coordinate with the largest data range is transformed to [0,1], and the
!  other coordinates are scaled by the same factor to preserve shape, along with
!  being shifted to the origin, before curvatures are calculated.  The purpose
!  is to avoid extreme curvature values from tiny geometries that can give
!  undesirable curvature-based point distributions.
!
!  Note that the data points are NOT normalized/denormalized in place.
!
!  11/25/13  D.A.Saunders  Version of CURVDIS2 used to derive CURVDIS3D2.
!  11/29/13    "      "    Higher level driving routine for CURVDIS3D as part
!                          of updating the latter to match CURVDIS refinements.
!
!  Author:  David Saunders, ERC, Inc. at NASA Ames Research Center, CA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)  :: ndat       ! Number of input curve points
   real,    intent (in)  :: xdat(ndat) ! X coordinates of input curve
   real,    intent (in)  :: ydat(ndat) ! Y coordinates of input curve
   real,    intent (in)  :: zdat(ndat) ! Z coordinates of input curve
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
                                       ! 3 => perform both sets of smoothings;
                                       ! 3 or 1 are suggested
   integer, intent (in)  :: lunout     ! Logical unit for showing convergence
                                       ! history; LUNOUT < 0 suppresses it
   logical, intent (in)  :: arcs_only  ! T => skip calculations of X & Y(1:N)
                                       ! but return revised arcs as X(1:N)
   real,    intent (out) :: x(n), &    ! X-Y-Z coordinates of output curve,
                            y(n), &    ! unless ARCS_ONLY = T, in which case
                            z(n)       ! just the arc-lengths are returned as
                                       ! X(1:N)
   integer, intent (out) :: ier        ! 0 => no errors;
                                       ! 1 => POWER is out of range;
                                       ! 2 => failure in ARBDIS
!  Local constants:

   integer, parameter :: ndim = 3
   real,    parameter :: zero = 0.

!  Local variables:

   integer :: ier_local
   real    :: xdata(ndat), ydata(ndat), zdata(ndat)
   real    :: pwr, scale(ndim), shift(ndim)

!  Procedures:

   external :: getscale  ! Get normalizing shift/scale factors
   external :: usescale  ! Apply normalizing/denormalizing factors
   external :: curvdis3d ! Redistribute points according to curvature

!  Execution:

   xdata(:) = xdat(:)
   ydata(:) = ydat(:)
   zdata(:) = zdat(:)

   call getscale ('G', ndim, ndat, xdata, ydata, zdata, scale, shift, ier_local)
   call usescale ('N', ndim, ndat, xdata, ydata, zdata, scale, shift, ier_local)

   pwr = power
   do  ! Until a low-enough exponent allows convergence

      call curvdis3d (ndat, xdata, ydata, zdata, n, pwr, ismooth, lunout, &
                      arcs_only, x, y, z, ier)
      if (ier == 0) exit

      write (*, '(a, i3, f5.1)') &
         'CURVDIS3D2: Troublein CURVDIS3D.  ier, power:', ier, pwr
      if (pwr == zero) then  ! Should never fail with zero
         write (*, '(a)') 'Aborting.'
         stop
      end if
      pwr = max (zero, pwr - 0.1)
      write (*, '(a, f5.1)') 'Retrying with lower power:', pwr
   end do

   if (arcs_only) then  ! x(:) is really arc lengths starting at 0.
      y(:) = x(:)       ! Avoid undefined y(:) & z(:) when usescale denormalizes
      z(:) = x(:)
      shift(:) = zero
   end if

   call usescale ('D', ndim, n, x, y, z, scale, shift, ier_local)

   end subroutine curvdis3d2
