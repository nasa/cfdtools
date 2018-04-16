!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine furthest_point (npatches, surface_patches, nquad, conn,          &
                              nline, xline, yline, zline, tline, lcs_method,   &
                              iquad, pint, qint, tint, xyzint, dotprodmax)
!  Purpose:
!
!     Determine the point on a 3-space line or curve defined by two or more
!  points that is "most inside" a convex structured surface in some sense.
!  This can help locate the two intersections between the line and the surface
!  (see LINE_SURFACE).
!
!  Strategy:
!
!     The two-point line case is explained most readily.  Any point on P1 P2
!  may be represented as P(t) = (1 - t) P1 + t P2.  The surface quadrilateral
!  closest to this point may be determined efficiently via an ADT search of the
!  surface grid.  If F is the surface point nearest to P, a 1-D maximization of
!  the dot product of the outward surface unit normal at F with the vector F - P
!  gives the point P that is furthest from the surface but inside it if that
!  dot product is positive, or the point on the line nearest to the surface but
!  outside it if the dot product is negative.
!
!  In the case of a curve, the local spline method of LCSFIT can be used
!  similarly to calculate (x,y,z) as functions of normalized arc length t.
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
!  03/24/06  DAS  Initial adaptation of INTSEC6 for the 2-point line case;
!                 leave hooks for the 3+-point line case.
!  08/05/13   "   All ADT variants have been merged as adt_utilities (generic
!                 build and search interfaces) for distribution purposes.
!
!  Author:  David Saunders, ELORET/NASA Ames Research Center, Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use grid_block_structure  ! For derived data type used by the ADT package
   use adt_utilities         ! All variants of the ADT build & search routines

   implicit none

!  Arguments:

   integer, intent (in)  :: npatches      ! # surface grid patches

   type (grid_type), intent (in) :: &
            surface_patches (npatches)    ! Surface grid; may be a volume grid
                                          ! of which only k = 1 is used
   integer, intent (in)  :: nquad         ! # surface quads., all patches

   integer, intent (in)  :: conn(3,nquad) ! Patch # and (i,j) for each quad.

   integer, intent (in)  :: nline         ! # points defining the line, >= 2

   real,    intent (in), dimension (nline) :: &
            xline, yline, zline, tline    ! Line coordinates, and work-space
                                          ! for normalized arc lengths.  Use
                                          ! tline(1) and tline(nline) to
                                          ! indicate the range of t (possibly
                                          ! beyond [0, 1]) within which to
                                          ! perform the minimization
   character, intent (in) :: lcs_method*1 ! 'L', 'M', or 'B' as for LCSFIT
                                          ! (if nline > 2)
   integer, intent (out)  :: iquad        ! conn(:,iquad) points to the best
                                          ! surface cell cell found
   real,    intent (out)  :: pint, qint   ! Corresponding interpolation coefs.
                                          ! in the unit square
   real,    intent (out)  :: tint         ! Normalized arc length along the
                                          ! line at the best point found
   real,    intent (out)  :: xyzint(3)    ! Surface point nearest to the
                                          ! intended intersection point
   real,    intent (out)  :: dotprodmax   ! Corresponding dot product as
                                          ! described above; a negative
                                          ! value means the line does not
                                          ! intersect the surface

!  Local constants:

   integer,   parameter :: lunout = -6    ! Suppress FMINRC iteration printout
   integer,   parameter :: nfmax  = 50    ! Limit on # function evaluations
   character, parameter :: caller * 14 = 'FURTHEST_POINT'

!  Local variables:

   integer :: istat, lunerr, numfun
   real    :: dotprodmin, t, ta, tb, tol

!  Execution:

!  A minimum can be found to within sqrt (machine epsilon), but avoid the sqrt:

   if (epsilon (tol) < 1.e-10) then
      tol = 1.e-8
   else
      tol = 1.e-4
   end if

   ta = tline(1)       ! 0 and 1 ...
   tb = tline(nline)   ! ... unless extrapolation is being permitted
   numfun = nfmax      ! Limit; FMINRC typically takes about 6 iterations
   lunerr = abs (lunout)
   istat = 2           ! Initialize the minimization

10 continue

      call fminrc (ta, tb, t, dotprodmin, tol, numfun, caller, lunout, istat)

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
   dotprodmax = -dotprodmin

!  Internal procedure for subroutine furthest_point:

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine objective ()

!     Calculate the squared distance from the line point P defined by t to the
!     nearest cell of the structured surface grid.  Dot product the unit normal
!     to the surface at point F with the vector P - F as a measure of being
!     inside the surface or not.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real,    parameter :: one = 1.0
      logical, parameter :: true = .true.

      integer :: n, ni, nj, i, j
      real    :: dsq, tm1, unit_normal(3), xyztarget(3)

      if (nline == 2) then

         tm1 = one - t

         xyztarget(1) = tm1 * xline(1) + t * xline(2)
         xyztarget(2) = tm1 * yline(1) + t * yline(2)
         xyztarget(3) = tm1 * zline(1) + t * zline(2)

      else ! Install LCSFIT if it's ever needed

         write (lunerr, '(/, 2a)') caller, ' needs to be extended for 3+ pts.'
         stop

      end if

      call search_adt (xyztarget, iquad, pint, qint, dsq, true, npatches,      &
                       surface_patches, nquad, conn, xyzint)

      n  = conn(1,iquad)  ! Patch #
      i  = conn(2,iquad)  ! Corner indices of relevant cell
      j  = conn(3,iquad)
      ni = surface_patches(n)%ni
      nj = surface_patches(n)%nj

      call surface_normal (ni, nj, surface_patches(n)%x, surface_patches(n)%y, &
                           surface_patches(n)%z, i, j, pint, qint, unit_normal)
      xyztarget(:) = xyztarget(:) - xyzint(:)

      dotprodmin = dot_product (xyztarget, unit_normal)

      end subroutine objective

   end subroutine furthest_point
