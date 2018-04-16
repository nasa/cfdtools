!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine intsec2d (nblock, grid, ncell, j, conn, &
                        nline, xline, yline, tline, lcs_method, &
                        icell, pint, qint, tint, xyint, dsq)
!  Purpose:
!
!     Calculate the intersection of a 2-space line or curve defined by two or
!  more points with a J line defined by the (:,J) points of a 2D structured
!  volume grid (1 or more blocks) where all blocks are indexed consistently
!  and J is probably 1 or %nj.
!
!  If the line/curve and grid curve do not actually meet, the line/curve may be
!  extrapolated (tint outside [0, 1]), but the grid J curve will not be extrap-
!  olated.  The point on the grid curve nearest to the [possibly extrapolated]
!  line will always be determined.  The shortest squared distance returned
!  indicates whether a true intersection was found or not (although a tangential
!  meeting is possible and not easily distinguished).
!
!  Strategy:
!
!     The two-point line case is explained most readily.  Any point on P1 P2
!  may be represented as P(t) = (1 - t) P1 + t P2.  The cell of the grid J curve
!  closest to this point may be determined efficiently via an ADT search of that
!  grid line.  A 1-D minimization of the squared distance w.r.t. t solves the
!  problem.  In the case of a curve, the local spline method of LCSFIT can be
!  used similarly to calculate (x,y) as functions of normalized arc length t.
!
!  Outline of search tree construction:
!
!  ncell = 0
!  do ib = 1, nblock
!     ncell = (grid(ib)%ni - 1) + ncell
!  end do
!
!  allocate (conn(2,ncell)) ! For block # and cell i (1 through %ni - 1)
!
!  call build_adt (nblock, grid, ncell, j, conn)
!
!  History:
!
!  08/05/13  DAS  Version of INTSEC6 (3-space structured surface intersection)
!                 from which this 2-space analogue has been derived.
!  03/31/15   "   INTSEC2D adapted from INTSEC6 and INTSEC2T.
!
!  Author:  David Saunders, ERC, Inc. at/NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use grid_block_structure  ! For derived data type used by the ADT package
   use adt_utilities         ! All variants of the ADT build & search utilities

   implicit none

!  Arguments:

   integer, intent (in)   :: nblock        ! blocks in the 2D structured grid

   type (grid_type), intent (in) :: grid (nblock)  ! 2D structured grid

   integer, intent (in)   :: ncell         ! # 2-pt. cells along a grid J line

   integer, intent (in)   :: j             ! J is probably 1 or %nj, defining a
                                           ! boundary curve formed by all blocks

   integer, intent (in)   :: conn(2,ncell) ! Block # and i for each 2-pt. cell

   integer, intent (in)   :: nline         ! # points defining the line, >= 2

   real,    intent (in), dimension (nline) :: &
            xline, yline                   ! Line coordinates

   real,    intent (in)   :: tline(nline)  ! Normalized arc lengths if nline > 2
                                           ! else 0 and 1 are set here.
   character, intent (in) :: lcs_method*1  ! 'L', 'M', or 'B' as for LCSFIT
                                           ! (if nline > 2)
   integer, intent (out)  :: icell         ! conn(:,icell) points to the best
                                           ! grid J line cell found
   real,    intent (out)  :: pint, qint    ! Corresponding interpolation coefs.
                                           ! applied to pts. icell, icell + 1
   real,    intent (out)  :: tint          ! Normalized arc length along the
                                           ! line at the best point found
   real,    intent (out)  :: xyint(2)      ! Grid J line point nearest to the
                                           ! intended intersection point
   real,    intent (out)  :: dsq           ! Corresponding squared distance
                                           ! to the nearest line point, >= 0.

!  Local constants:

   integer,       parameter :: lunout = -6 ! < 0 suppresses FMINRC printout
   integer,       parameter :: nfmax  = 50 ! Limit on # function evaluations
   real,          parameter :: one    = 1.
   logical,       parameter :: true   = .true.
   character (8), parameter :: caller = 'INTSEC2D'

!  Local variables:

   integer :: istat, lunerr, numfun
   real    :: t, ta, tb, tol
   logical :: two_points

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
   numfun = nfmax      ! Limit; FMINRC typically takes about 6 iterations
   lunerr = abs (lunout)
   istat  = 2          ! Initialize the minimization

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

!  Internal procedure for subroutine intsec6:

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine objective ()

!     Calculate the squared distance from the line point defined by t to the
!     nearest cell of the curve formed by line J of the structured 2D grid.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real    :: deriv, r, rm1, xytarget(2)
      integer :: neval = 1

      if (two_points) then

         r   = (t - ta) / (tb - ta)
         rm1 = one - r
         xytarget(1) = rm1 * xline(1) + r * xline(2)
         xytarget(2) = rm1 * yline(1) + r * yline(2)

      else  ! Spline interpolation at t:

         call lcsfit (nline, tline, xline, true, lcs_method, neval, t, &
                      xytarget(1), deriv)
         call lcsfit (nline, tline, yline, true, lcs_method, neval, t, &
                      xytarget(2), deriv)
      end if

!!!   write (*, '(a, 2es19.11, a, 3es19.11)') &
!!!      ' xytarget:', xytarget, ' ta, tb, t:      ', ta, tb, t

      call search_adt (xytarget, icell, pint, qint, dsq, true, nblock, &
                       grid, ncell, j, conn, xyint)

!!!   write (*, '(a, 3i5)') 'icell,ib,i:', icell, conn(:,icell)
!!!   write (*, '(a, 2es19.11, a, 3es19.11)') &
!!!      ' xyint:   ', xyint, ' pint, qint, dsq:', pint, qint, dsq

      end subroutine objective

   end subroutine intsec2d
