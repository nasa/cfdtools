!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine line_surface (linenum, npatches, surface_patches, nquad, conn,   &
                            nline, xline, yline, zline, sline, lcs_method,     &
                            ns, s1, s2)
!  Purpose:
!
!     LINE_SURFACE is a specialized driver of the lower-level INTSEC6 utility
!  for the case of finding the two possible intersections between a 3-space
!  curve and a closed convex surface.  The curve is defined by 2 or more points
!  and the surface is a multipatch structured surface grid such as the inner
!  or outer boundary of a CFD grid for hypersonic flow calculations on a convex
!  body such as a capsule.

!     This version implements the 2-point line case, but hooks are provided as
!  for INTSEC6 in case a more general curve is called for.
!
!  Original Strategy:
!
!     INTSEC6 solves the intersection problem for a given range of parametric
!  variable t defining distance along the curve.  The curve (most likely a
!  straight line) is assumed to intersect the surface either not at all or at
!  two points.  The possibility of coincident points is extremely unlikely in
!  the intended application to a capsule being observed during atmospheric
!  entry.  (Returning a value of ns = 1 is possible, but a derivative-based
!  method would be a more reliable way to calculate a tangency point defined
!  by a flat plane rather than a straight line.)
!
!     The input arc lengths in tline(*) should not be normalized, because the
!  line length is used to scale the distance tolerance that determines whether
!  intersections exist or not.
!
!     If no intersection is found (as indicated by INTSEC6 argument dsq > tol),
!  then ns will be set to 0 and we're done.
!
!     If a first intersection is found at normalized line location t1, then the
!  intervals [0, t1] and [t1, 1] are searched for the second intersection that
!  is presumed to exist.  One of these searches should recalculate t1, while the
!  other will be the second solution t2.  [In practice, the t1 boundaries are
!  fudged a little because looking for a minimum right at an end-point is not a
!  good idea for convergence reasons.]
!
!     Denormalized locations s1 and s2 are returned in ascending order along
!  with ns = 2.  Note that recovering (x,y,z)s is left to the calling program,
!  while returning surface details at the two intersections would require many
!  more arguments that may never be needed in practice.
!
!  Revised Strategy:
!
!     Unfortunately, the above scheme cannot be relied upon not the find the
!  same minimum THREE times.  It works some of the time but not always.
!
!     What is really needed is a better way of determining the two intervals
!  normally containing one minimum each.  Strategy:  Maximize a measure of
!  being inside the convex surface.  Assuming the surface patches are right-
!  handed, such a measure for point P(t) on the line is the dot product of the
!  surface normal at the foot of the perpendicular F(t) on the surface with
!  the vector (F - P).  If this is negative, the point is outside the surface,
!  so if the maximum as a function of t is negative, the entire line is outside
!  the surface.  Otherwise, the maximum is positive and its location provides
!  a good definition of two intervals in which to look for one intersection.
!  The possibility of tangency remains, so if the maximized dot product is
!  below some slightly positive tolerance, don't look for any intersections.
!
!
!  Further Notes:
!
!     See FURTHEST_POINT for details of maximizing the distance of a point
!  on a line to a convex structured surface.
!
!     See INTSEC6 for details of how a 3-variable line-surface intersection
!  problem is treated with a 1-variable non-derivative minimization method using
!  an efficient tree search package to locate the nearest surface point to any
!  line point that is typically off the surface except at an intersection where
!  the squared distance being minimized goes to within a tolerance of zero.
!
!     The following search tree construction outline is reproduced from INTSEC6:
!
!  nquad = 0
!  do ib = 1, npatches
!    nquad = (surface_patches(ib)%ni - 1) * (surface_patches(ib)%nj - 1) + nquad
!  end do
!
!  allocate (conn(3,nquad)) ! For patch # and (i,j)
!
!  call build_adt (npatches, surface_patches, nquad, conn)
!
!  History:
!
!  06/13/05  DAS  Initial implementation of INTSEC6 for the 2-point line case,
!                 with a single intersection assumed.
!  03/01/06   "   Initial implementation of LINE_SURFACE (still just the 2-point
!                 line case) to help with radiative heating estimates derived
!                 from CFD solutions for the Stardust capsule.  This original
!                 strategy described above cannot be made reliable.
!  03/24/06   "   Revised strategy using new FURTHEST_POINT utility.
!  09/17/20   "   Hayabusa 2 produced puzzling radiation hot spots with NEQAIR'S
!                 black body BC used correctly.  We believe that the aft body's
!                 concavities produced flawed line-surface intersections from
!                 the method of intsec6 that assumes a convex surface.  This
!                 version invokes intsec9, which can trap an intersection
!                 squared distance not essentially zero and use a back-up method
!                 to isolate the true intersection.  Added the linenum argument
!                 for improved diagnostics.
!
!  Author:  David Saunders, ELORET/NASA Ames Research Center, Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use grid_block_structure  ! For the derived data type used by the ADT package

   implicit none

!  Arguments:

   integer, intent (in)  :: linenum       ! Line number for disgnostic purposes

   integer, intent (in)  :: npatches      ! # surface grid patches

   type (grid_type), intent (in) :: &
            surface_patches (npatches)    ! Surface grid; may be a volume grid
                                          ! of which only k = 1 is used
   integer, intent (in)  :: nquad         ! # surface quads., all patches

   integer, intent (in)  :: conn(3,nquad) ! Patch # and (i,j) for each quad.

   integer, intent (in)  :: nline         ! # points defining the line, >= 2
                                          ! (incomplete for more than 2)
   real,    intent (in), dimension (nline) :: &
            xline, yline, zline, sline    ! Line coordinates and UNnormalized
                                          ! arc lengths
   character, intent (in) :: lcs_method*1 ! 'L', 'M', or 'B' as for LCSFIT
                                          ! (if nline > 2)
   integer, intent (out)  :: ns           ! # intersection points found (0 or 2
                                          ! normally; ns = 1 indicates tangency)
   real,    intent (out)  :: s1, s2       ! UNnormalized arc lengths along the
                                          ! line of the intersections if ns = 2;
                                          ! ns = 1 indicates ~tangency at s1
!  Local constants:

   real, parameter :: zero = 0., one = 1.

!  Local variables:

   integer :: iquad, lunerr, numfun
   real    :: dotprodmax, dsq, pint, qint, sa, sb, &
              tint0, tint1, tint2, tol
   real    :: xyzint(3)
   real, allocatable :: tline(:)

!  Execution:

!  Establish a tolerance for the distance of the point on the line that is
!  closest to the surface, and scale it by the line length:

   if (epsilon (tol) < 1.e-10) then
      tol = 1.e-8
   else
      tol = 1.e-4
   end if

   sa  = sline(1);  sb = sline(nline)
   tol = (sb - sa) * tol

   allocate (tline(nline))

   tline(:) = (sline(:) - sa) / (sb - sa) ! Normalized arcs

!  Look for the point on the line that is "most inside" the surface.
!  If it is not inside the surface at all, we are done.

   call furthest_point (npatches, surface_patches, nquad, conn,                &
                        nline, xline, yline, zline, tline, lcs_method,         &
                        iquad, pint, qint, tint0, xyzint, dotprodmax)

!! write (6, '(a, i7, a, 2f12.8, a, f12.8, a, 3f13.8, a, 1p, e14.6)')          &
!!    ' iquad:', iquad, '  (p,q):', pint, qint, '  t:', tint0, '  xyzint:',    &
!!    xyzint, '  dpmax:', dotprodmax

   if (dotprodmax < zero) then ! No intersection; we're done
      ns = 0
      go to 99
   else if (dotprodmax < tol**2) then
      ns = 1
      s1 = sa + (sb -sa) * tint0
      go to 99
   end if

!  Look for an intersection in [0, tint0].

   tline(2) = tint0 ! For a 2-point line; rethink this for a general curve

   call intsec9 (linenum, npatches, surface_patches, nquad, conn,              &
                 nline, xline, yline, zline, tline, lcs_method,                &
                 iquad, pint, qint, tint1, xyzint, dsq)

!! if (dsq > 1.e-6) write (52, '(a)') '???'
!! write (52, '(a, i6, es16.8)') 'line, dsq1:', linenum, dsq

!  Look for a second intersection in [tint0, 1]:

   tline(1) = tint0;  tline(2) = one

   call intsec9 (linenum, npatches, surface_patches, nquad, conn,              &
                 nline, xline, yline, zline, tline, lcs_method,                &
                 iquad, pint, qint, tint2, xyzint, dsq)

!! if (dsq > 1.e-6) write (52, '(a)') '%%%'
!! write (52, '(a, i6, es16.8)') 'line, dsq2:', linenum, dsq

!  Sort the two values of t:

   if (tint2 < tint1) then
      tint0 = tint1
      tint1 = tint2
      tint2 = tint0
   end if

   s1 = sa + (sb -sa) * tint1
   s2 = sa + (sb -sa) * tint2
   ns = 2

99 deallocate (tline)

   end subroutine line_surface
