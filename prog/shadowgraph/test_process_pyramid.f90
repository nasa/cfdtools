!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   program test_process_pyramid

!  05/27/2010  David Saunders

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

   integer, parameter :: lunkbd = 5, luncrt = 6

   integer :: i, npts, inode(5)
   real    :: x(5), y(5), z(5), f(5), fx(5), fz(5), xt, zt
   real    :: ypts(2), fpts(3,2)

!  Execution:

   write (luncrt, '(a)') 'The quad face must be face 1.'

   do i = 1, 5
      write (luncrt, '(a, i2, a)', advance='no') &
         'x/y/z, f(1:3) for pyramid vertex', i, ': '
      read  (lunkbd, *, end=99) x(i), y(i), z(i), f(i), fx(i), fz(i)
      inode(i) = i
   end do

   do ! Until ^D
      write (luncrt, '(a)', advance='no') 'x, z for line _|_ y = 0 plane: '
      read  (lunkbd, *, end=99) xt, zt

      call process_pyramid (5, inode, x, y, z, f, fx, fz, xt, zt, &
                            npts, ypts, fpts)

      write (luncrt, '(a, i2, 1p, 2e14.6, 2(1x, 3e14.6))') &
         '   n, y, f: ', npts, ypts, fpts
   end do

99 continue

   end program test_process_pyramid
