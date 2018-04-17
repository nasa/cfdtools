!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine angle_between_planes (p1, p2, p3, q1, q2, q3, fixed_p, angle)

!  Determine the angle between two planes, where each plane is defined by the
!  vertices of one triangle.  The triangles may or may not have a common edge.
!
!  The desired angle is the angle between vectors normal to the triangles.
!  Thus an existing unit normal utility is made use of.  If the planes are
!  parallel, the result could be either 0 or 180 degrees, depending on the
!  order of the vertices in the two triangles (application-dependent).
!
!  Provision is made for not recomputing the first normal if only the second
!  triangle is changing between calls.
!
!  It is up to the application to avoid passing a degenerate triangle (with
!  two vertices identical, meaning the angle is undefined).
!
!  10/27/2014  D.A.Saunders  Adaptation of tri_normal_and_area utility.
!
!  Author:  David Saunders, ERC, Inc./NASA Ames Research Center, Moffett Field.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use trigd
   implicit none

!  Argument:

   real, dimension (3), intent (in)  :: p1, p2, p3   ! Triangle 1 vertices
   real, dimension (3), intent (in)  :: q1, q2, q3   ! Triangle 2 vertices
   logical,             intent (in)  :: fixed_p      ! T if the normal for
                                                     ! triangle p1-p2-p3 is
                                                     ! known from a prior call
   real,                intent (out) :: angle        ! Desired angle (degrees)

!  Local constants:

   integer, parameter :: nnode = 6  ! 3 vertices for 2 triangles
   integer, parameter :: ntri  = 2  ! 2 triangles

!  Local variables:

   integer, save    :: conn(3,ntri)      ! Pointers to 3 vertices per triangle
   real,    save    :: tri_xyz(3,nnode)  ! Coordinates for 2 x 3 vertices
   real             :: areap, areaq      ! Unneeded triangle areas
   real,    save    :: upsq              ! |uq|**2
   real,    save    :: up(3)             ! Unit normal for triangle p
   real             :: uq(3)             ! Unit normal for triangle q

!  Execution:

   if (.not. fixed_p) then  ! First triangle
      tri_xyz(:,1) = p1(:)
      tri_xyz(:,2) = p2(:)
      tri_xyz(:,3) = p3(:)
      conn(1,1)    = 1
      conn(2,1)    = 2
      conn(3,1)    = 3
      conn(1,2)    = 4
      conn(2,2)    = 5
      conn(3,2)    = 6

      call tri_normal_and_area (nnode, ntri, tri_xyz, conn, 1, areap, up)

      upsq = dot_product (up, up)
   end if

   tri_xyz(:,4) = q1(:)
   tri_xyz(:,5) = q2(:)
   tri_xyz(:,6) = q3(:)

   call tri_normal_and_area (nnode, ntri, tri_xyz, conn, 2, areaq, uq)

   angle = acosd (dot_product (up, uq) / sqrt (upsq * dot_product (uq, uq)))

   end subroutine angle_between_planes
