!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program cone_angles
!
!  For a given number of points along an edge of a hemispherical quadrant, show
!  the cone-angles associated with each of the points.
!
!  David Saunders  08/07/2018  Initial Implementation.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Constants:

   integer, parameter :: lunkbd = 5
   integer, parameter :: luncrt = 6

!  Variables:

   integer :: i, nedge
   real    :: angle, delta

!  Execution:

   write (luncrt, '(a)', advance='no') '# hemisphere quadrant edge points: '
   read  (lunkbd, *) nedge
   delta = 90./real (nedge - 1)

   do i = 1, nedge
      angle = delta*real (i - 1)
      write (luncrt, '(i4, f10.6)') i, angle
   end do

   end program cone_angles
