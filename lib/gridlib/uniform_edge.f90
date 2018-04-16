!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine uniform_edge (nd, xd, yd, zd, md, method, nu, xu, yu, zu)

!  For a list of points in 3-space representing a line or edge, interpolate a
!  (roughly) uniform point distribution.  The local spline interpolation of x,
!  y, z with respect to approximate arc length performed by PLSCRV3D (q.v.) is
!  applied to what might be points obtained by probing a dataset graphically.
!
!  If there are only two data points, the output points are truly uniform, else
!  the best we can do is evaluate the spline at uniform chord lengths.  There is
!  no attempt (yet) to ensure capturing of interior data points - just the end
!  points are included exactly.
!
!  Provision has been made for ensuring linear interpolation where desired, but
!  it's awkward.  The number of data points is assumed to be small, though, so
!  it shouldn't be too burdensome to supply a flag for each data interval.
!  Thus, the interpolation will either be the type specified by "method" or
!  linear if the flag for the current interval requests it:  md(id) = 1 means do
!  linear interpolation in the data interval id : id + 1, id = 1 : nd - 1.
!
!  History:
!
!     06/10/2008  D.A.Saunders  Initial implementation, prompted by the need to
!                               represent a cloud of laser-scanned data points
!                               from a wedge of a recovered heat shield as a set
!                               of structured surface patches.
!
!  Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer,   intent (in)                  :: nd         ! # data points >= 2
   real,      intent (in),  dimension (nd) :: xd, yd, zd ! Data coordinates
   integer,   intent (in),  dimension (nd) :: md         ! See description above
   character, intent (in)                  :: method * 1 ! 'B' = "loose" fit
                                                         ! 'M' = tight/monotonic
                                                         ! 'L' = linear
   integer,   intent (in)                  :: nu         ! # ~uniform points
   real,      intent (out), dimension (nu) :: xu, yu, zu ! Interpolated points

!  Local constants:

   real,    parameter :: derivs = -999.  ! Suppresses derivatives
   real,    parameter :: arrow  =  1.0   ! Arc lengths increase from zero
   logical, parameter :: true   = .true.
   logical, parameter :: false  = .false.

!  Local variables:

   integer   :: id, iu
   real      :: dt, td(nd), tdtotal, tu
   character :: type * 1

!  Execution:

   call chords3d (nd, xd, yd, zd, false, tdtotal, td)  ! Don't normalize

   dt = tdtotal / real (nu - 1)

   xu(1) = xd(1);  yu(1) = yd(1);  zu(1) = zd(1)

   id = 1

   do iu = 2, nu - 1

      tu = dt * real (iu - 1)
     
      call interval (nd, td, tu, arrow, id)  ! id points to the left data index

      if (md(id) == 1) then
         type = 'L'
      else
         type = method
      end if
!                                              new   not closed
      call plscrv3d (nd, xd, yd, zd, td, type, true, false, tu, id, &
                     xu(iu), yu(iu), zu(iu), derivs)
   end do

   xu(nu) = xd(nd);  yu(nu) = yd(nd);  zu(nu) = zd(nd)

   end subroutine uniform_edge
