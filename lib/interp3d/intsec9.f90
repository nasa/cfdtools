!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine intsec9 (linenum, npatches, surface_patches, nquad, conn, &
                       nline, xline, yline, zline, tline, lcs_method,   &
                       iquad, pint, qint, tint, xyzint, dsq)
!  Purpose:
!
!     This is an adaptation of intsec6 that has a recovery method for occasions
!  when the minimized squared distance between the line and (structured) surface
!  is not essentially zero.  The recovery method simply discretizes the line
!  finely, finds the squared distance of each point from the surface (quite
!  cheap because of the ADT search method) and repeats the minimization method
!  on a much smaller interval virtually guaranteed to contain the intersection.
!  See Strategy below for turning what is really a 3-variable nonlinear set of
!  equations into a 1-variable minimization.  That strategy assumes that the
!  surface is convex.  This adaptation should handle surfaces that contain
!  concavities.
!
!     The following description is mostly from intsec6:
!
!     Calculate the intersection of a 3-space line or curve defined by two or
!  more points with a structured surface defined by one or more surface patches.
!  If the line and curve do not actually meet, the curve may be extrapolated
!  (tint outside [0, 1]), but the surface will not be extrapolated.  The point
!  on the surface nearest to the [possibly extrapolated] line will always be
!  determined if the surface is convex.  The shortest squared distance returned
!  indicates whether a true intersection has been found or not (although a
!  tangential meeting is possible and not easily distinguished).  If there is a
!  dip in the surface (such as in the wake part of a flow field shock boundary
!  for an atmospheric entry vehicle, or a concavity in the aft body part of such
!  a vehicle's surface, a back-up method of finding the true intersection is
!  needed, and that is what this intsec9 variation does.
!
!  Strategy (Reliable for a Convex Surfaxe):
!
!     The two-point line case is explained most readily.  Any point on P1 P2
!  may be represented as P(t) = (1 - t) P1 + t P2.  The surface quadrilateral
!  closest to this point may be determined efficiently via an ADT search of the
!  surface grid.  A 1-D minimization of the squared distance w.r.t. t solves the
!  problem.  In the case of a curve, the local spline method of LCSFIT can be
!  used similarly to calculate (x,y,z) as functions of normalized arc length t.
!
!  Outline of search tree construction:
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
!  06/13/05  DAS  Initial intsec6 implementation for the 2-point line case;
!                 leave hooks for the 3+-pt. line case, piecewise linear or not.
!  08/25/06   "   Completed the 3+-point case.  The fact that the arc lengths
!                 are normalized still allows tline(1) and tline(n) to be used
!                 to define the search interval (because we know those arcs are
!                 0 and 1), but two other arguments would have been better.
!                 The 0 and 1 have to be substituted here, then ta and tb are
!                 restored before the return.
!                 Requiring tline(:) as input (as opposed to deriving arcs from
!                 x/y/zline(:) here) is considered better since the application
!                 is likely to need those arc lengths anyway.
!  08/05/13   "   All ADT variants have been merged as adt_utilities (generic
!                 build and search interfaces) for distribution purposes.
!  09/17/20   "   Intsec9 variant of intsec6, with the backup-method first
!                 devised for program SLOS incorporated here instead of in the
!                 application progam.  Invoke intsec6 twice rather than making
!                 an internal procedure version of it as initially tried.  The
!                 option to treat multipoint curves got in the way, but it has
!                 been retained with this approach.  No: discretizing the line
!                 for the recovery method should use PLSCRV3D.
!                 Add the linenum argument for the diagnostics.
!  09/18/20   "   Last bug: The recovery method can't assume that the tline(:)
!                 argument goes from 0 to 1.  Also, 128 points may be too few
!                 for the crude search, so double it.
!  10/08/20   "   Comment out excess recovery method printing, and ease the
!                 dtolsq tolerance by making it 100 times bigger.
!  01/15/21   "   Fixed a few typos, and removed the objective internal pro-
!                 cedure that was unwittingly carried over from intsec6
!
!  Author:  David Saunders, ELORET, Inc. at NASA Ames Research Center, CA.
!           Later with AMA, Inc. at NASA ARC.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use grid_block_structure  ! For derived data type used by the ADT package
   use adt_utilities         ! All variants of the ADT build & search utilities

   implicit none

!  Arguments:

   integer, intent (in)  :: linenum        ! Line nnumber for diagnostics

   integer, intent (in)  :: npatches       ! # surface grid patches

   type (grid_type), intent (in) :: &
            surface_patches (npatches)     ! Surface grid; may be a volume grid
                                           ! of which only k = 1 is used
   integer, intent (in)  :: nquad          ! # surface quads., all patches

   integer, intent (in)  :: conn(3,nquad)  ! Patch # and (i,j) for each quad.

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
   integer, intent (out)  :: iquad         ! conn(:,iquad) points to the best
                                           ! surface cell found
   real,    intent (out)  :: pint, qint    ! Corresponding interpolation coefs.
                                           ! in the unit square
   real,    intent (out)  :: tint          ! Normalized arc length along the
                                           ! line at the best point found
   real,    intent (out)  :: xyzint(3)     ! Surface point nearest to the
                                           ! intended intersection point
   real,    intent (out)  :: dsq           ! Corresponding squared distance
                                           ! to the nearest line point

!  Local constants:

   integer,   parameter :: lunout = -6     ! Suppress FMINRC iteration printout
   integer,   parameter :: nfmax  = 50     ! Limit on # function evaluations
   real,      parameter :: one    = 1.
   real,      parameter :: zero   = 0.
   logical,   parameter :: false  = .false.
   logical,   parameter :: true   = .true.

!  Local variables:

   integer :: ieval, istat, lunerr, numfun
   real    :: t, ta, tb, tol, dtolsq, xyztarget(3)
   logical :: new_plscrv3d, two_points
   logical, save :: first = .true.

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
   dtolsq = 1000.*(tol*(max (tb - ta, one)))**2

   if (first) then
       first = false
       write (*, '(a, es16.8)') 'dtolsq:', dtolsq
   end if

!  The following is vulnerable to a local minimum if the surface is not convex;
!  it returns a minimized dsq:

   call intsec6 (npatches, surface_patches, nquad, conn,        &
                 nline, xline, yline, zline, tline, lcs_method, &
                 iquad, pint, qint, tint, xyzint, dsq)

   if (dsq > dtolsq) then
      write (*, '(a, i6, es13.6)') &
         ' *** Intersection trouble. Line #, dsqmin:', linenum, dsq

      call bad_local_min_recovery ()  ! Brute force recovery method

   end if

!  Internal procedure for subroutine intsec9:

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine bad_local_min_recovery ()

!     Crude search for a best small interval, followed by efficient refinement.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer, parameter :: nu = 256 ! # uniform points for fine discretization
      integer, parameter :: nl = 2   ! Always a 2-point line for intsec6 here

      integer :: i, iuse, ilocmin(1)
      real    :: dtu, distsq(nu), sl(nl), su(nu), xu(nu), yu(nu), zu(nu)

!     Discretize the total arc length to a reasonable fineness.
!     Remember they are NORMALIZED arc lengths.

      dtu = (tline(nline) - tline(1))/real (nu-1)
      do i = 1, nu
         su(i)        = tline(1) + dtu*real (i-1)
         xu(i)        = xline(1) + (xline(nline) - xline(1))*su(i)
         xyztarget(1) = xu(i)
         yu(i)        = yline(1) + (yline(nline) - yline(1))*su(i)
         xyztarget(2) = yu(i)
         zu(i)        = zline(1) + (zline(nline) - zline(1))*su(i)
         xyztarget(3) = zu(i)

         call search_adt (xyztarget, iquad, pint, qint, distsq(i), true, &
                          npatches, surface_patches, nquad, conn, xyzint)

!!!      write (*, '(a, i4, es14.6, 2x, 3es14.6, 2x, es14.6, 2x, 3es14.6)') &
!!!         'i, s(i), dsq(i):', i, su(i), xyztarget(:), distsq(i), xyzint(:)
      end do

      ilocmin = minloc (distsq)
      i       = ilocmin(1)
      write (*, '(a, i4)') 'Index of minimum distance:', i
      iuse    = max (1, i-1)
      sl(1)   = su(iuse)
      iuse    = min (nu, i+1)
      sl(2)   = su(iuse)

      call intsec6 (npatches, surface_patches, nquad, conn,     &
                    nline, xline, yline, zline, sl, lcs_method, &
                    iquad, pint, qint, tint, xyzint, dsq)

      if (dsq < dtolsq) then
         write (*, '(a, i6, 4es13.6)') &
            '     Recovered.            Line #, sa, sb, dsq, s:', &
            linenum, sl(:), dsq, tint
      else
         write (*, '(a, i6, 4es13.6)') &
            ' *** Intersection failure. Line #, sa, sb, dsq, s:', &
            linenum, sl(:), dsq, tint

         ! Keep going
      end if

      end subroutine bad_local_min_recovery

   end subroutine intsec9
