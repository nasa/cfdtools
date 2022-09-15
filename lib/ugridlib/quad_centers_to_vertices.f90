!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine quad_centers_to_vertices (nnode, ncell, nvert, nf, xyz, conn, &
                                       fcenter, area, area_total, fnode)
!
!  Interpolate cell-centered function values to the vertices of a 3-space
!  unstructured quad surface zone via area-weighted averaging of the functions
!  of the cells sharing common vertices.
!
!  Ideally, the input function array would be updated with new dimensions in
!  place, but Fortran 90 seems not to handle real, pointer :: f(:,:) upon input.
!  Therefore, it is up to the application to do any in-place substitution.
!
!  03/17/2010  D.A.Saunders  Adaptation of triangle techniques from program
!                            ADB2ADSI as tri_centers_to_vertices.
!  03/18/2022    "      "    Adaptation of tri_centers_to_vertices as
!                            vol_centers_to_vertices.
!  04/08/2022    "      "    Adaptation of vol_centers_to_vertices as
!                            quad_centers_to_vertices.
!
!  Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA
!           Now with AMA, Inc. at NASA ARC.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)  :: nnode             ! # vertices in the zone
   integer, intent (in)  :: ncell             ! # cells in the zone
   integer, intent (in)  :: nvert             ! # vertices per cell
   integer, intent (in)  :: nf                ! # functions to interpolate
   real,    intent (in)  :: xyz(3,nnode)      ! (x,y,z)s of all nodes
   integer, intent (in)  :: conn(nvert,ncell) ! Pointers into xyz(:,*)
   real,    intent (in)  :: fcenter(nf,ncell) ! Cell-centered function values
   real,    intent (out) :: area(ncell)       ! Cell areas
   real,    intent (out) :: area_total(nnode) ! Area sums at each vertex
   real,    intent (out) :: fnode(nf,nnode)   ! Vertex-centered function values

!  Local constants:

   real, parameter :: zero = 0.

!  Local variables:

   integer :: icell, inode, iv, ivert

!  Execution:

   call quad_areas (nnode, ncell, xyz, conn, area, area_total)

   fnode(:,:) = zero
   do icell = 1, ncell
      do ivert = 1, nvert
         iv = conn(ivert,icell)
         fnode(:,iv) = fnode(:,iv) + area(icell)*fcenter(:,icell)
      end do
   end do

   do inode = 1, nnode
      fnode(:,inode) = fnode(:,inode) / area_total(inode)
   end do

99 return

   end subroutine quad_centers_to_vertices
