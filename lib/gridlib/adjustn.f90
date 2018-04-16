!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine adjustn (ndata, i1, i2, x1, y1, z1, ia, ib, x2, y2, z2, method)

!  ADJUSTN is a slight generalization of CHANGEN to allow for data points
!  outside the range of interest, as needed for treating periodic data.
!
!  Both redistribute the points along a curve represented by x1, y1, z1 so that
!  data points between i1 and i2 are converted to interpolated points ia:ib,
!  where those indices imply the changed number of points.
!
!  Unlike CHANGEN, ADJUSTN can use data points i1 - 1 and i2 + 1 if they are
!  present: data arrays indexed 1:ndata are passed to PLSCRV3D, not i1:i2.
!
!  Normally, the character of the input distribution is preserved in the new
!  spacing, but the output points are forced to be uniform (or as near to
!  uniform as use of uniform accumulated chord lengths allows) if a lowercase
!  "method" argument is provided instead of the uppercase choice recognized
!  by the local cubic spline interpolation scheme of PLSCRV3D.
!
!  12/22/97  DAS  3-space CHANGEN grid line utility wrapped around PLSCRV3D.
!  01/11/08   "   Added option to make the output distribution uniform.
!                 See expanded use of "method".
!  01/14/08   "   ADJUSTN adapted from CHANGEN as explained above, with a new
!                 name, since adding an argument would affect existing uses.
!
!  Author:  David Saunders, ELORET Corporation/NASA Ames, Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer,   intent (in)  :: ndata      ! Number of data points, 1:ndata
   integer,   intent (in)  :: i1, i2     ! Pts. i1:i2 are redistributed as ia:ib
   real,      intent (in)  :: x1(ndata), y1(ndata), z1(ndata)  ! Data points
   integer,   intent (in)  :: ia, ib     ! First & last redistributed points
   real,      intent (out) :: x2(ib), y2(ib), z2(ib)   ! Redistributed coords.
   character, intent (in)  :: method * 1 ! Type of fit to be used by PLSCRV3D:
                                         ! M(onotonic) | B ("loose") | L(inear);
                                         ! lowercase m | b | l produce output
                                         ! pts. that are (essentially) uniform
!  Procedures:

   external chords3d,   &                ! Arc-length utility
            plscrv3d,   &                ! 3-space local cubic spline utility
            upcase                       ! Character string utility

!  Local constants:

   real, parameter    :: derivs = -999.  ! Suppresses them
   real, parameter    :: one    = 1.
   logical, parameter :: norm   = .false.! No need to normalize
   logical, parameter :: closed = .false.

!  Local variables:

   integer   :: i, ieval, ip
   real      :: t(ndata), dt, p, r, ri1, rip, teval, total
   logical   :: new, uniform
   character :: mupper*1

!  Execution:

   mupper = method

   call upcase (mupper)

   uniform = mupper /= method

   call chords3d (ndata, x1, y1, z1, norm, total, t)

   ieval  = i1
   new    = .true.

   x2(ia) = x1(i1)
   y2(ia) = y1(i1)
   z2(ia) = z1(i1)

   if (uniform) then

      dt = (t(i2) - t(i1)) / real (ib - ia)

      do i = ia + 1, ib - 1

         teval = t(i1) + dt * real (i - ia)

         call plscrv3d (ndata, x1, y1, z1, t, mupper, new, closed, &
                        teval, ieval, x2(i), y2(i), z2(i), derivs)
         new = .false.

      end do

   else  ! Preserve original form of spacing as much as possible

      r   = real (i2 - i1) / real (ib - ia)
      ri1 = real (i1)

      do i = ia + 1, ib - 1

         rip   = ri1 + r * real (i - ia)
         ip    = int (rip)
         p     = rip - real (ip)
         teval = (one - p) * t(ip) + p * t(ip+1)

         call plscrv3d (ndata, x1, y1, z1, t, mupper, new, closed, &
                        teval, ieval, x2(i), y2(i), z2(i), derivs)
         new = .false.

      end do

   end if

   x2(ib) = x1(i2)
   y2(ib) = y1(i2)
   z2(ib) = z1(i2)

   end subroutine adjustn
