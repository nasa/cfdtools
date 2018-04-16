!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine tri_centers_to_vertices (nnode, ntri, nf, area, area_total,      &
                                       conn, fcenter, fnode)
!
!  Interpolate cell-centered function values to the vertices of a 3-space
!  triangulated surface zone via area-weighted averaging of the functions of
!  the triangles sharing common vertices.
!
!  Ideally, the input function array would be updated with new dimensions in
!  place, but Fortran 90 seems not to handle real, pointer :: f(:,:) upon input.
!  Therefore, it is up to the application to do any in-place substitution.
!
!  03/17/2010  D.A.Saunders  Adaptation of techniques from program ADB2ADSI.
!
!  Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)  :: nnode             ! # vertices in the triangulation
   integer, intent (in)  :: ntri              ! # triangles in the zone
   integer, intent (in)  :: nf                ! # functions to interpolate
   real,    intent (in)  :: area(ntri)        ! Triangle areas
   real,    intent (in)  :: area_total(nnode) ! Area sums at each vertex
   integer, intent (in)  :: conn(3,ntri)      ! Pointers into xyz(:,*)
   real,    intent (in)  :: fcenter(nf,ntri)  ! Cell-centered function values
   real,    intent (out) :: fnode(nf,nnode)   ! Vertex-centered function values

!  Local constants:

   real, parameter :: zero = 0.

!  Local variables:

   integer :: i1, i2, i3, inode, itri

!  Execution:

   fnode(:,:) = zero

   do itri = 1, ntri
      i1 = conn(1,itri);  fnode(:,i1) = fnode(:,i1) + area(itri)*fcenter(:,itri)
      i2 = conn(2,itri);  fnode(:,i2) = fnode(:,i2) + area(itri)*fcenter(:,itri)
      i3 = conn(3,itri);  fnode(:,i3) = fnode(:,i3) + area(itri)*fcenter(:,itri)
   end do

   do inode = 1, nnode
      fnode(:,inode) = fnode(:,inode) / area_total(inode)
   end do

   end subroutine tri_centers_to_vertices
