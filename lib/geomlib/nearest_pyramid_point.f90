!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine nearest_pyramid_point (xvert, xtarget, xnear, coefs, dsq)
!
!  To be completed as needed.  See nearest_tet_point for an analogue.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   real, intent (in)  :: xvert(3,8)  ! Cell vertex coordinates in Fluent order
   real, intent (in)  :: xtarget(3)  ! Target pt. being tested against this cell
   real, intent (out) :: xnear(3)    ! Point found closest to the target that is
                                     ! not outside this cell
   real, intent (out) :: coefs(8)    ! Interpolation coefficients associated
                                     ! with xnear, in [0, 1] and summing to 1
   real, intent (out) :: dsq         ! The corresponding shortest distance
                                     ! between xtarget and xnear, squared
!  Execution:

   xnear(:) = 0.                     ! Avoid compiler warnings
   coefs(:) = 0.
   dsq      = 0.

   end subroutine nearest_pyramid_point
