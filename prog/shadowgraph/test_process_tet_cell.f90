!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   program test_process_tet_cell

!  05/27/2010  David Saunders

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

   integer, parameter :: lunkbd = 5, luncrt = 6

   integer :: i, npts, inode(4)
   real    :: x(4), y(4), z(4), f(4), fx(4), fz(4), xt, zt
   real    :: ypts(2), fpts(3,2)

!  Execution:

   do i = 1, 4
      write (luncrt, '(a, i2, a)', advance='no') &
         'x/y/z, f(1:3) for tet vertex', i, ': '
      read  (lunkbd, *, end=99) x(i), y(i), z(i), f(i), fx(i), fz(i)
      inode(i) = i
   end do

   do ! Until ^D
      write (luncrt, '(a)', advance='no') 'x, z for line _|_ y = 0 plane: '
      read  (lunkbd, *, end=99) xt, zt

      call process_tet_cell (4, inode, x, y, z, f, fx, fz, xt, zt, &
                             npts, ypts, fpts)

      write (luncrt, '(a, i2, 1p, 2e14.6, 2(1x, 3e14.6))') &
         '   n, y, f: ', npts, ypts, fpts
   end do

99 continue

   end program test_process_tet_cell
