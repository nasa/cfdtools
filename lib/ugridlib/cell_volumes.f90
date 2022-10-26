!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine cell_volumes (nnode, ncell, nvert, xyz, conn, volume, &
                            volume_total, ier)
!
!  Calculate the volumes of the cells of one zone of a 3-space unstructured
!  volume mesh.  Tet or hex cells are treated initially; others may need to be
!  handled.  All cells of one zone are assumed to have the same cell type.
!
!  Also return, for each vertex, the summed volumes of cells sharing that
!  vertex, as needed for volume-weighted interpolation of function data from
!  cell centers to cell vertices for ADT searching.  Note that this doesn't
!  handle the possibility of nodes shared by adjacent zones.
!
!  The hex cell vertex ordering is assumed to be that of Fluent: counterclock-
!  wise for a base quad. face, 1234, and similarly for the opposite face, 5678.
!
!  03/17/2010  D.A.Saunders  Initial implementation of tri_areas.
!  03/13/2022    "      "    Analogous cell_volumes (tet and hex cases first).
!  03/17/2022    "      "    Switched the hexahedron method from the three
!                            quad-based pyramids of tet_volume to the six
!                            tetrahedra of hex_vol.
!
!  Author:  David Saunders, AMA, Inc./NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)  :: nnode                ! # nodes in the zone
   integer, intent (in)  :: ncell                ! # cells in the zone
   integer, intent (in)  :: nvert                ! # vertices in one cell, 4|8
   real,    intent (in)  :: xyz(3,nnode)         ! (x,y,z)s of all nodes
   integer, intent (in)  :: conn(nvert,ncell)    ! Vertex pointers into xyz(:,*)
   real,    intent (out) :: volume(ncell)        ! Cell volumes
   real,    intent (out) :: volume_total(nnode)  ! Summed volumes surrounding
                                                 ! each node
   integer, intent (out) :: ier                  ! 0|1 = no error|bad nvert

!  Local constants:

   real, parameter :: zero = 0.

!  Local variables:

   integer :: i1, i2, i3, i4, i5, i6, i7, i8, icell

!  Procedures:

   external :: &
      tet_volumes, &  ! Analogue of cell_volumes for the all-tet case
      hex_vol         ! Volume of a hex cell via six tet volumes

!  Execution:

   ier = 0
   volume_total(:) = zero

   select case (nvert)

     case (4)  ! Tetrahedral cells; this had already been implemented

        call tet_volumes (nnode, ncell, xyz, conn, volume, volume_total)

     case (8)  ! Hexahedral cells

        do icell = 1, ncell
           i1 = conn(1,icell)
           i2 = conn(2,icell)
           i3 = conn(3,icell)
           i4 = conn(4,icell)
           i5 = conn(5,icell)
           i6 = conn(6,icell)
           i7 = conn(7,icell)
           i8 = conn(8,icell)

           call hex_vol (xyz(:,i1), xyz(:,i2), xyz(:,i3), xyz(:,i4), &
                         xyz(:,i5), xyz(:,i6), xyz(:,i7), xyz(:,i8), &
                         volume(icell))

           volume_total(i1) = volume_total(i1) + volume(icell)
           volume_total(i2) = volume_total(i2) + volume(icell)
           volume_total(i3) = volume_total(i3) + volume(icell)
           volume_total(i4) = volume_total(i4) + volume(icell)
           volume_total(i5) = volume_total(i5) + volume(icell)
           volume_total(i6) = volume_total(i6) + volume(icell)
           volume_total(i7) = volume_total(i7) + volume(icell)
           volume_total(i8) = volume_total(i8) + volume(icell)
        end do


     case default

        write (*, '(a, i4)') &
           '*** Cell_volumes: unhandled # vertices per cell:', nvert
        ier = 1

   end select

   end subroutine cell_volumes
