!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine quad_to_tri (header, zone, ios)
!
!  Purpose:
!
!     Convert one quad-cell surface zone to an equivalent triangulation, as a
!     new single zone, in place.
!
!     The input surface is expected to be in Fluent format, with the vertices
!     of each quad cell in circular order.  Vertices 1 and 3 of each quad are
!     joined to form the two triangles derived from it.  Only the header and
!     the zone connectivity pointers are affected.
!
!  History:
!
!     04/01/2022  D.A.Saunders  Initial design, with USLOS in mind (inner and
!                               outer boundaries).
!     04/12/2022    "      "    The test case produced triangle normals pointing
!                               into the body surface, so the index orders have
!                               now been reversed for each pair of triangles.
!     04/15/2022    "      "    After switching to triangulating all zones of a
!                               quad surface, switched back to a single zone.
!                               Concatenating all zones into header%xyz and
!                               %conn is now done with new combine_zones in
!                               triangulation_io.f90.
!     04/29/2022    "      "    Nasty one: we can't change header%nvertices
!                               until the last zone is done.  This means it has
!                               to be changed at the higher level.
!
!  Author:  David Saunders, AMA, Inc. at NASA Ames Research Center, CA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use tri_header_structure  ! See triangulation_io for both
   use tri_zone_structure

   implicit none

!  Arguments:

   type (tri_header_type), intent (inout) :: header  ! Input with quad surface;
                                                     ! output as a triangulation
   type (tri_type),        intent (inout) :: zone    ! Quad zone to convert

   integer,                intent (out)   :: ios     ! 0 = no error

!  Local variables:

   integer :: ie, ne, nelements, nnodes, nvertices
   integer, allocatable :: icontemp(:,:)

!  Execution:

   nvertices = header%nvertices
   if (nvertices /= 4) then
      write (*, '(a, i6)') '*** quad_to_tris: nvertices is not 4:', nvertices
      ios = 1
      go to 99
   end if

   nnodes    = zone%nnodes
   nelements = zone%nelements

   allocate (icontemp(nvertices,nelements))
   icontemp(:,:) = zone%conn(:,:)
   deallocate (zone%conn)
   allocate (zone%conn(3,2*nelements))

   ie = 0
   do ne = 1, nelements
      ie = ie + 1
      zone%conn(1:3,ie) = icontemp(3:1:-1,ne);  ie = ie + 1
      zone%conn(1:2,ie) = icontemp(4:3:-1,ne)
      zone%conn(3,  ie) = icontemp(1,     ne)
   end do

   deallocate (icontemp)

   zone%nelements    = 2*nelements
   zone%zone_type    = 'FETRIANGLE'
   zone%element_type = 'TRI'
   header%zone_type  = zone%zone_type
!! header%nvertices  = 3  ! This means all zones have to be the same type
!  No: this affects the next zone; do it in the calling program.

99 continue

   end subroutine quad_to_tri
