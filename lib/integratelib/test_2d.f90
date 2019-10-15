   program test_2d

!  Test 2-D quadrature on a generated structured surface.

   implicit none

   integer, parameter :: nf = 2
   integer :: i, j, n, ni, nj
   real    :: volume(nf)
   real, allocatable, dimension (:,:) :: x, y
   real, allocatable, dimension (:,:,:) :: f
   character :: filename * 32

   write (6, '(/, a)', advance='no') ' f = f(x,y) file name: '
   read  (5, *) filename
   open  (1, file=filename, status='old')  ! From GEN2D in SMOOTH2D directory
   read  (1, *)
   read  (1, *) ni, nj

   allocate (x(ni,nj), y(ni,nj), f(nf,ni,nj))

   read  (1, *) ((x(i,j), y(i,j), f(1,i,j), i = 1, ni), j = 1, nj)

   close (1)

   do n = 2, nf
      f(n,:,:) = f(1,:,:)
   end do

   call trapezoidal_quadrature_2d (ni, nj, nf, x, y, f, volume)

   write (6, '(/, a, 1p, 3e15.7)') ' Volumes:', volume(:)

   end program test_2d
