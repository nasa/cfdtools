!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine intsec8 (nnode, ntri, conn, coord,                               &
                       nline, xline, yline, zline, tline, lcs_method,          &
                       itri, pint, qint, rint, tint, xyzint, dsq)
!
!     This variant of INTSEC7 applies to a triangulated surface rather than a
!  structured surface of quadrilateral cells.
!
!     The following documentation is slightly adapted from that in INTSEC7.
!
!  Purpose:
!
!     Calculate the intersection of a 3-space line or curve defined by two or
!  more points with a triangulated surface defined by a list of triangles.
!  If the curve and surface do not actually meet, the curve may be extrapolated
!  (tint outside [0, 1]), but the surface will not be extrapolated.  The point
!  on the surface nearest to the [possibly extrapolated] line will always be
!  determined.  The shortest squared distance returned indicates whether a true
!  intersection was found or not (although a tangential meeting is possible and
!  not easily distinguished).
!
!  Strategy:
!
!     The two-point line case is explained most readily.  Any point on P1 P2
!  may be represented as P(t) = (1 - t) P1 + t P2.  The surface triangle closest
!  to this point may be determined efficiently via an ADT search of the surface
!  grid.  A 1-D minimization of the squared distance w.r.t. t solves the
!  problem.  In the case of a curve, the local spline method of LCSFIT can be
!  used similarly to calculate (x,y,z) as functions of normalized arc length t.
!
!  Search tree construction:
!
!     allocate (conn(3,ntri))   ! Triangle connectivity info.
!     allocate (coord(3,nnode)) ! (x,y,z)s of nodes pointed to be conn(1:3,*)
!
!     call build_unstructured_surface_adt (nnode, ntri, conn, coord)
!
!  History:
!
!  06/13/05  DAS  Initial implementation of INTSEC6 for the 2-point line case
!                 and a structured surface.
!  08/25/06   "   nline >= 2 INTSEC6 case.  The fact that the arc lengths are
!                 normalized still allows tline(1) and tline(n) to be used to
!                 define the search interval (because we know those arcs are
!                 0 and 1), but two other arguments would have been better.
!                 The 0 and 1 have to be substituted here, then ta and tb are
!                 restored before the return.
!                 Requiring tline(:) as input (as opposed to deriving arcs from
!                 x/y/zline(:) here) is considered better since the application
!                 is likely to need those arc lengths anyway.
!  04/20/07   "   INTSEC7 adjustment of INTSEC6 allows more than one variant of
!                 the ADT package to be used in a single applications.
!  02/25/13   "   Adaptation of INTSEC7 for the triangulated surface case.
!  08/05/13   "   All ADT variants have been merged as adt_utilities (generic
!                 build and search interfaces) for distribution purposes.  This
!                 makes INTSEC7 a redundant substitute for INTSEC6 now.
!
!  Author:  David Saunders, ERC, Inc. at NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use adt_utilities         ! All variants of the ADT build & search utilities

   implicit none

!  Arguments:

   integer, intent (in)  :: nnode          ! # (x,y,z)s packed in coord(1:3,*)

   integer, intent (in)  :: ntri           ! # triangles in the surface list

   integer, intent (in)  :: conn(3,ntri)   ! Node pointers into coord(:,1:nnode)

   real,    intent (in)  :: coord(3,nnode) ! (x,y,z)s of triangulated nodes

   integer, intent (in)  :: nline          ! # points defining the line, >= 2

   real,    intent (in), dimension (nline) :: &
            xline, yline, zline            ! Line coordinates

   real,    intent (inout) :: tline(nline) ! Input with normalized arc lengths
                                           ! if nline > 2.  Use tline(1) and
                                           ! tline(nline) to indicate the range
                                           ! of t (possibly beyond [0, 1]) that
                                           ! should contain the intersection;
                                           ! see the 08/25/06 history above
   character, intent (in) :: lcs_method*1  ! 'L', 'M', or 'B' as for LCSFIT
                                           ! (if nline > 2)
   integer, intent (out)  :: itri          ! conn(:,itri) points to the best
                                           ! surface triangle found
   real,    intent (out)  :: pint,qint,rint! Corresponding interpolation coefs.:
                                           ! xi = p x1 + q x2 + r x3
   real,    intent (out)  :: tint          ! Normalized arc length along the
                                           ! line at the best point found
   real,    intent (out)  :: xyzint(3)     ! Surface point nearest to the
                                           ! intended intersection point
   real,    intent (out)  :: dsq           ! Corresponding squared distance
                                           ! to the nearest line point

!  Local constants:

   integer,   parameter :: lunout = -6     ! Suppress FMINRC iteration printout
   integer,   parameter :: nfmax  = 50     ! Limit on # function evaluations;
   real,      parameter :: one    = 1.     ! ~10 iterations appear typical
   real,      parameter :: zero   = 0.
   logical,   parameter :: false  = .false.
   logical,   parameter :: true   = .true.
   character, parameter :: caller * 7 = 'INTSEC8'

!  Local variables:

   integer :: ieval, istat, lunerr, numfun
   real    :: t, ta, tb, tol
   logical :: new_plscrv3d, two_points

!  Execution:

!  A minimum can be found to within sqrt (machine epsilon), but avoid the sqrt:

   if (epsilon (tol) < 1.e-10) then
      tol = 1.e-8
   else
      tol = 1.e-4
   end if

   two_points = nline == 2

   ta = tline(1)       ! 0 and 1 ...
   tb = tline(nline)   ! ... unless extrapolation is being permitted

   if (.not. two_points) then
      tline(1)     = zero
      tline(nline) = one
      ieval        = nint (0.5 * (ta + tb) * real (nline))
      ieval        = min (nline - 1, max (1, ieval))
      new_plscrv3d = true
   end if

   numfun = nfmax      ! Limit; FMINRC typically takes about 6 iterations
   lunerr = abs (lunout)
   istat = 2           ! Initialize the minimization

10 continue

      call fminrc (ta, tb, t, dsq, tol, numfun, caller, lunout, istat)

      if (istat < -1) then

         write (lunerr, '(/, 2a)') caller, ': FMINRC fatal error'
         stop

      else if (istat < 0) then ! Iteration limit; may be usable

         write (lunerr, '(/, 2a)') caller, ': Iteration limit.'

      else if (istat > 0) then ! Evaluate the objective function at t

         call objective   ! Internal procedure below
         go to 10

      else ! istat = 0 (success).

      end if

!  Ensure that everything matches t(best), not t(last):

   call objective ()

   tint = t

   if (.not. two_points) then
      tline(1)     = ta
      tline(nline) = tb
   end if

!  Internal procedure for subroutine intsec8:

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine objective ()

!     Calculate the squared distance from the line point defined by t to the
!     nearest triangle of the unstructured surface grid.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real :: tm1, xyztarget(3)
      real, save :: derivs = -999. ! -999. suppresses partial deriv. calcs.

      if (two_points) then

         tm1 = one - t

         xyztarget(1) = tm1 * xline(1) + t * xline(2)
         xyztarget(2) = tm1 * yline(1) + t * yline(2)
         xyztarget(3) = tm1 * zline(1) + t * zline(2)

      else  ! Spline interpolation at t:

         call plscrv3d (nline, xline, yline, zline, tline, lcs_method,         &
                        new_plscrv3d, false, t, ieval, xyztarget(1),           &
                        xyztarget(2), xyztarget(3), derivs)

         new_plscrv3d = false  ! After the first call
      end if

!!    write (*, '(a, 1p, 3es19.11)') ' xyztarget:', xyztarget

      call search_adt (xyztarget, itri, pint, qint, rint, dsq, true, nnode,    &
                       ntri, conn, coord, xyzint)

!!    write (*, '(a, 5es19.11)') ' p, q, r, t, dsq:', pint, qint, rint, t, dsq

      end subroutine objective

   end subroutine intsec8
