!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine inside_outside (npatches, surface_patches, nquad, conn, &
                              px, py, pz, dsq, inside)
!  Purpose:
!
!     Determine whether the given point in 3-space appears to be "inside" or
!  "outside" a structured surface grid.  By definition, "outside" means the
!  side that the local surface normal points into if the point is close to
!  the surface.  For a typical right-handed surface grid on a convex body,
!  the normals point outside the body.  Some bodies such as axisymmetric
!  bodies of revolution, while convex in the rotational direction, may have
!  surface concavities in the axial direction.  The aft body of the Mars
!  Science Laboratory aeroshell is an example which has prompted this utility.
!  Its aft-body cone angles are not monotonically increasing.
!
!     Thorough estimation of radiative heat flux at a body point during entry
!  into an atmosphere can be done by constructing lines of sight from the
!  body point in many directions covering a hemisphere centered tangentially
!  at the body point.  Concavities mean some of those hemisphere lines may
!  encounter the body (inner volume grid surface) before they reach the outer
!  volume grid surface.  Such lines need to be clipped appropriately (or maybe
!  have their flow data partially zeroed out) before the radiation calculations
!  are performed along each line.
!
!     The intended use in this application is to check each hemisphere line
!  starting at the body point end, and locate the first line point, if any,
!  that appears to be inside the body and not outside it.  Flow data at such a
!  point and beyond can be suppressed because it cannot contribute to the
!  radiation at the body point.
!
!  Strategy:
!
!     For a given point P being checked, the nearest surface point F may be
!  determined efficiently via an ADT search.  If the dot product of F - P
!  with the unit normal at F is positive, then P should be inside the surface.
!  This is clear if P is near the surface because F will then be the foot of
!  the normal from P to the surface.  As P moves further from the surface,
!  F - P may no longer be normal to the surface (it just defines the shortest
!  distance) but at the time of writing it is believed that the dot product
!  will still answer the inside/outside question.
!
!  Outline of preliminary search tree construction:
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
!  Then search_adt accesses the search tree through internal module ADT_DATA,
!  and returns details of the (interpolated) surface point nearest to each
!  target point.
!
!  History:
!
!  03/24/06  DAS  Initial FURTHEST_POINT utility from which INSIDE_OUTSIDE
!                 has been adapted.
!  02/17/15   "   Initial implementation of INSIDE_OUTSIDE.
!
!  Author:  David Saunders, ERC, Inc. at NASA Ames Research Center, Field, CA
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

   real,    intent (in)  :: px, py, pz    ! Current point being checked

   real,    intent (out) :: dsq           ! Shortest squared distance from
                                          ! the point to the surface
   logical, intent (out)  :: inside       ! T if the point appears to be
                                          ! "inside" the surface
!  Local constants:

   real,    parameter :: zero = 0.
   logical, parameter :: true = .true.

!  Local variables:

   integer :: i, iquad, j, n, ni, nj
   real    :: pint, qint, unit_normal(3)
   real    :: xyzint(3), xyztarget(3)

!  Execution:

   xyztarget(1) = px
   xyztarget(2) = py
   xyztarget(3) = pz

   call search_adt (xyztarget, iquad, pint, qint, dsq, true, npatches, &
                    surface_patches, nquad, conn, xyzint)

   n  = conn(1,iquad)  ! Patch #
   i  = conn(2,iquad)  ! Corner indices of relevant cell
   j  = conn(3,iquad)
   ni = surface_patches(n)%ni
   nj = surface_patches(n)%nj

   call surface_normal (ni, nj, surface_patches(n)%x, surface_patches(n)%y, &
                        surface_patches(n)%z, i, j, pint, qint, unit_normal)

   xyztarget(:) = xyzint(:) - xyztarget(:)

   inside = dot_product (xyztarget, unit_normal) > zero

   end subroutine inside_outside
