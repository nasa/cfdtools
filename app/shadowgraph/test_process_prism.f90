!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   program test_process_prism

!  05/27/2010  David Saunders

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

   integer, parameter :: lunkbd = 5, luncrt = 6

   integer :: i, npts, inode(6)
   real    :: x(6), y(6), z(6), f(6), fx(6), fz(6), xt, zt
   real    :: ypts(2), fpts(3,2)

!  Execution:

   write (luncrt, '(a)') 'Follow the FieldView node numbering convention.'

   do i = 1, 6
      write (luncrt, '(a, i2, a)', advance='no') &
         'x/y/z, f(1:3) for prism vertex', i, ': '
      read  (lunkbd, *, end=99) x(i), y(i), z(i), f(i), fx(i), fz(i)
      inode(i) = i
   end do

   do ! Until ^D
      write (luncrt, '(a)', advance='no') 'x, z for line _|_ y = 0 plane: '
      read  (lunkbd, *, end=99) xt, zt

      call process_prism (6, inode, x, y, z, f, fx, fz, xt, zt, &
                          npts, ypts, fpts)

      write (luncrt, '(a, i2, 1p, 2e14.6, 2(1x, 3e14.6))') &
         '   n, y, f: ', npts, ypts, fpts
   end do

99 continue

   end program test_process_prism
