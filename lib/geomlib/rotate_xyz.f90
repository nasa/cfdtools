!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine rotate_xyz (n, xyz, angle, p, q)
!
!     ROTATE_XYZ is a variant of ROTATE3D for the case where x/y/z are stored
!  as triples, xyz(1:3,:).  As with the earlier routine, it rotates the given
!  point(s) in 3-space about the straight line joining points P (Px, Py, Pz)
!  and Q (Qx, Qy, Qz) by the indicated angle (degrees, right-hand rule).
!  Input coordinates are overwritten.  For consistency, P and Q are also
!  entered as 3-vectors.
!
!  HISTORY: 06/14/00  DAS  ROTATE3D adapted from ROTATE2D and ROTATE program.
!           02/28/14   "   ROTATE_XYZ adapted from ROTATE3D.
!
!  AUTHOR: David Saunders, ERC, Inc./NASA Ames Research Center, CA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Arguments:

   integer, intent (in) :: n          ! Number of points to rotate
   real, intent (inout) :: xyz(3,n)   ! Point coordinates, overwritten
   real, intent (in)    :: angle      ! Angle in degrees; RH rule (thumb P -> Q)
   real, intent (in)    :: p(3), q(3) ! Points defining the rotation axis

!  Local constants:

   real, parameter :: one = 1.0, zero = 0.0

!  Local variables:

   integer :: i, j, k
   real    :: cth, sth, vth, d, vec(3), vecrot(3), rot(3,3)

!  Intrinsics:

   real    :: cosd, sind

!  Execution:

   cth = cosd (angle)
   sth = sind (angle)
   vth = one - cth

!  Derive a unit vector translated to the origin from the vector PQ:

   vec(:) = q(:) - p(:)
   d      = sqrt (dot_product (vec, vec))
   vec(:) = vec(:) / d

!  Rotation matrix ("Introduction to Robotics" by John Craig, 1986):

   rot(1,1) = vec(1) * vec(1) * vth + cth
   rot(1,2) = vec(1) * vec(2) * vth - vec(3) * sth
   rot(1,3) = vec(1) * vec(3) * vth + vec(2) * sth

   rot(2,1) = vec(1) * vec(2) * vth + vec(3) * sth
   rot(2,2) = vec(2) * vec(2) * vth + cth
   rot(2,3) = vec(3) * vec(2) * vth - vec(1) * sth

   rot(3,1) = vec(1) * vec(3) * vth - vec(2) * sth
   rot(3,2) = vec(2) * vec(3) * vth + vec(1) * sth
   rot(3,3) = vec(3) * vec(3) * vth + cth

!  Translate each point by P's coordinates, rotate about the unit axis
!  at the origin, then translate back:

   do k = 1, n
      vec(:) = xyz(:,k) - p(:)

      do i = 1, 3
         vecrot(i) = zero
         do j = 1, 3
            vecrot(i) = vecrot(i) + rot(i,j) * vec(j)
         end do
      end do

      xyz(:,k) = vecrot(:) + p(:)
   end do

   end subroutine rotate_xyz
