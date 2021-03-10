!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine curve_diffs (n1, x1, y1, n2, x2, y2, shortest_diff)

!  Description:
!
!     This utility calculates a measure of the distance between two related
!  curves in 2-space that should be preferable to the standard measure of
!  vertical distances at the defining abscissas, which can be misleading in the
!  high-gradient regions.  It is the underlying curves, not just the discretized
!  points on them, that should be of interest.  The curves are treated as
!  piecewise linear because of how the shortest-distance utility used works.
!  They are likely to be time histories.
!
!     For each of the points on curve 1, the shortest distance to the curve
!  defined by the given points on curve 2 is determined. A utility first written
!  for an edge grid form of earlier surface and volume grid searching methods
!  makes the method straightforward.
!
!     The sign assigned to shortest_diff(i) is that of the cross product of the
!  chord length at point i on curve 1 with the vector between that point and the
!  nearest point found on curve 2.  If the curves cross at a curve 1 point, the
!  shortest difference is zero, but this is unlikely in practice.
!
!     All calculations are performed on normalized forms of the curves, in order
!  to handle arbitrary units for x and y.  The normalizations do NOT overwrite
!  the input data.  The units of the differences are somewhat problematic.
!  Since the differences curve should not be affected (much) by the resolution
!  of the abscissas of curve 1, the units of shortest_diffs are denormalized by
!  the curve 1 shift and scale, meaning they are roughly the units of y1.
!  (Both curves are normalized with one shift and scale, namely those from
!  curve 1, so as not to distort the relative shapes.)
!
!  History:
!
!     01/19/2021  D.A.Saunders  Initial implementation, prompted by FIAT_Opt.
!     01/20/2021    "     "     Fixed initially too-large residuals.
!
!  Author:  David Saunders, AMA, Inc. at NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)  :: n1                  ! # points defining curve 1
   real,    intent (in)  :: x1(n1), y1(n1)      ! 2-space coordinates of curve 1
   integer, intent (in)  :: n2                  ! # points defining curve 2
   real,    intent (in)  :: x2(n2), y2(n2)      ! 2-space coordinates of curve 2
   real,    intent (out) :: shortest_diff(n1)   ! Shortest differences from all
                                                ! curve 1 points to curve 2, in
                                                ! the units of y1
!  Local constants:

   integer,       parameter :: ndim      =  2   ! Everything is 2-space
   real,          parameter :: one       =  1.
   character (1), parameter :: norm_mode = 'I'  ! Independent x/y normalizations
   character (1), parameter :: norm      = 'N'  ! Normalize
   character (1), parameter :: denorm    = 'D'  ! Denormalize

!  Local variables:

   integer :: i1, i2, i2best, ier, it
   real    :: dsq, dsqmin, p, q
   real, dimension (ndim) :: scale1, shift1, p1, p2, pbest, pnear, x1y1, &
                             xybest, v1, v2, yd, zd
   real, allocatable, dimension (:) :: x1norm, x2norm, y1norm, y2norm

!  Procedures:

   external ::  &
      getscale, &    ! Gets scale and shift factors to normalize coordinates
      usescale, &    ! Applies scale and shift factors to normalize/denormalize
      nearest_edge_pt_2d  ! ... from each normalized crv1 pt. to normalized crv2

!  Execution:

!  Normalize the curves, x and y being treated independently:

   call getscale (norm_mode, ndim, n1, x1, y1, zd, scale1, shift1, ier)
   allocate (x1norm(n1), y1norm(n1))
   x1norm(:) = x1(:)
   y1norm(:) = y1(:)
   call usescale (norm, ndim, n1, x1norm, y1norm, zd, scale1, shift1, ier)

!! write (21, '(2es16.8)') (x1norm(i1), y1norm(i1), i1 = 1, n1)

   allocate (x2norm(n2), y2norm(n2))
   x2norm(:) = x2(:)
   y2norm(:) = y2(:)
   call usescale (norm, ndim, n2, x2norm, y2norm, zd, scale1, shift1, ier)

!! write (22, '(2es16.8)') (x2norm(i2), y2norm(i2), i2 = 1, n2)

!  For each point on normalized c1, find the nearest point on normalized c2:

   do i1 = 1, n1
      dsqmin  = 1.e+10
      x1y1(1) = x1norm(i1);  x1y1(2) = y1norm(i1)
      do i2 = 1, n2 - 1  ! For each interval along normalized curve 2
         p1(1) = x2norm(i2);  p2(1) = x2norm(i2+1)
         p1(2) = y2norm(i2);  p2(2) = y2norm(i2+1)

         call nearest_edge_pt_2d (p1, p2, x1y1, pnear, p, q, dsq)

         if (dsq < dsqmin) then
             dsqmin = dsq
!!           i2best = i2
             xybest = pnear
         end if
      end do
      it = min (i1+1, n1)
      v1(1) = x1norm(it) - x1norm(it-1)  ! Vector from i1/crv. 1 (w/ care at n1)
      v1(2) = y1norm(it) - y1norm(it-1)
      v2(1) = xybest(1)  - x1norm(i1) ! Vector to closest point on curve 2
      v2(2) = xybest(2)  - y1norm(i1)
      shortest_diff(i1)  = sqrt (dsqmin) * sign (one, v1(1)*v2(2) - v1(2)*v2(1))
   end do

!! write (23, '(2es16.8)') (x1(i1), shortest_diff(i1), i1 = 1, n1)

   shortest_diff(:) = shortest_diff(:)*scale1(2)

   end subroutine curve_diffs
