!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine usintsec (linenum, nnode, ntri, conn, coord,             &
                        nline, xline, yline, zline, tline, lcs_method, &
                        itri, pint, qint, rint, tint, xyzint, dsq)
!  Purpose:
!
!     This is the unstructured surface adaptation of structured grid utility
!     intsec9, for intersecting a 3-space line or curve defined by two or more
!     points with a triangulated surface defined as a list of elements.  Like
!     intsec9, it has a recovery method for occasions when the minimized squared
!     distance between the line and the surface is not essentially zero.  The
!     recovery method simply discretizes the line finely, finds the squared
!     distance of each point from the surface (quite cheap because of the ADT
!     search method) and repeats the minimization method on a much smaller
!     interval virtually guaranteed to contain the intersection.
!
!     See Strategy below for turning what is really a 3-variable nonlinear set
!     of equations into a 1-variable minimization.  That strategy assumes that
!     the surface is convex. This adaptation should handle surfaces that contain
!     concavities.
!
!  Strategy (Reliable for a Convex Surfaxe):
!
!     The two-point line case is explained most readily.  Any point on P1 P2 may
!     be represented as P(t) = (1 - t) P1 + t P2.  The surface element closest
!     to this point may be determined efficiently via an ADT search of the surf-
!     ace grid.  A 1-D minimization of the squared distance w.r.t. t solves the
!     problem.  In the case of a curve, the local spline method of LCSFIT can be
!     used similarly to calculate (x,y,z) as functions of normalized arc length
!     t.
!
!  History:
!
!     01/15/21  DAS  Date of intsec9 source from which this is derived.
!     03/05/22   "   Initial adaptation of usintsec from intsec9, using the
!                    existing intsec8 (twice) in place of intsec6.
!     09/30/22   "   bad_local_min_recovery was using tline instead of sl(1:2)
!                    in the recovery intsec8 call.
!
!  Author:  David Saunders, AMA,, Inc. at NASA Ames Research Center, CA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use adt_utilities         ! All variants of the ADT build & search utilities

   implicit none

!  Arguments:

   integer, intent (in)  :: linenum        ! Line number for diagnostics

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
   integer,   parameter :: nfmax  = 50     ! Limit on # function evaluations
   real,      parameter :: one    = 1.
   real,      parameter :: zero   = 0.
   logical,   parameter :: false  = .false.
   logical,   parameter :: true   = .true.

!  Local variables:

   integer :: i, ieval, istat, lunerr, numfun
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

   call intsec8 (nnode, ntri, conn, coord,                      &
                 nline, xline, yline, zline, tline, lcs_method, &
                 itri, pint, qint, rint, tint, xyzint, dsq)

   if (dsq > dtolsq) then
      write (*, '(a, i6, es13.6)') &
         ' *** Usintsec trouble. Line #, dsqmin:', linenum, dsq
      write (*, '(a, /, (i4, 4es16.7))') ' i, x, y, z, t: ', &
         (i, xline(i), yline(i), zline(i), tline(i), i = 1, nline)

      call bad_local_min_recovery ()  ! Brute force recovery method

   end if

!  Internal procedure for subroutine usintsec::

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine bad_local_min_recovery ()

!     Crude search for a best small interval, followed by efficient refinement.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer, parameter :: nu = 256 ! # uniform points for fine discretization
      integer, parameter :: nl = 2   ! Always a 2-point line for intsec8 here

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

         call search_adt (xyztarget, itri, pint, qint, rint, distsq(i), true, &
                          nnode, ntri, conn, coord, xyzint)

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

      call intsec8 (nnode, ntri, conn, coord,                      &
                    nline, xline, yline, zline, sl, lcs_method, &
                    itri, pint, qint, rint, tint, xyzint, dsq)

      if (dsq < dtolsq) then
         write (*, '(a, i6, 4es13.6)') &
            '     Recovered.            Line #, sa, sb, dsq, s:', &
            linenum, sl(:), dsq, tint
      else
         write (*, '(a, i6, 4es13.6)') &
            ' *** USINTSEC: Intersection failure. Line #, sa, sb, dsq, s:', &
            linenum, sl(:), dsq, tint

         ! Keep going
      end if

      end subroutine bad_local_min_recovery

   end subroutine usintsec
