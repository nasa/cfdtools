   program test_angle2d

   implicit none

   integer :: i
   real    :: angle
   real    :: x(3), y(3)

   write (*, '(a)', advance='no') 'Enter x(1:3): '
   read  (*, *) x(:)
   write (*, '(a)', advance='no') 'Enter y(1:3): '
   read  (*, *) y(:)

   call angle2d (3, 2, x, y, angle)

   write (*, '(f12.7)') angle

   end program test_angle2d
