!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      program combine
!
!     Form a linear combination of a pair of (x,y) datasets.  In order to deal
!     with differing data ranges, the two sets of abscissas are normalized to
!     [0, 1] prior to forming the combination.  The result is denormalized to
!     a range specified interactively.  In order to deal with differing numbers
!     of points and possibly differing distributions, the output number of
!     uniformly distributed points is also prompted for.  Time histories are
!     the intended application.
!
!     Data format:
!
!     Title
!     n
!     x1  y1
!     x2  y2
!     :   :
!     xn  yn
!
!     11/02/00  DAS  Initial implementation, for starting trajectories.
!     08/26/02   "   Don't require roughly the same data range.
!
!     David Saunders, ELORET/NASA Ames Research Center, Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      integer, parameter ::
     >   lundat = 1, lunkbd = 5, luncrt = 6

      real, parameter ::
     >   one = 1., zero = 0.
 
      logical, parameter ::
     >   new = .true., normalize = .true.

      character, parameter ::
     >   method * 1 = 'M'

      integer
     >   i, n, n1, n2

      real
     >   dx, w1, w2, xa, xb

      real, allocatable, dimension (:) ::
     >   x, x1, x2, y, y1, y1int, y2, y2int, yp

      character
     >   filename * 48

!     Execution:

      write (luncrt, '(/, a, /)')
     >   ' Linear combination of two datasets.'

      write (luncrt, '(a)', advance='no')
     >   ' Primary dataset?   '
      read  (lunkbd, '(a)') filename
      open  (lundat, file=filename, status='old')
      read  (lundat, '(a)') ! Title
      read  (lundat, *) n1
      allocate (x1(n1), y1(n1))
      read  (lundat, *) (x1(i), y1(i), i = 1, n1)
      close (lundat)

       write (luncrt, '(a)', advance='no')
     >    ' Secondary dataset? '
      read  (lunkbd, '(a)') filename
      open  (lundat, file=filename, status='old')
      read  (lundat, '(a)') ! Title
      read  (lundat, *) n2
      allocate (x2(n2), y2(n2))
      read  (lundat, *) (x2(i), y2(i), i = 1, n2)
      close (lundat)

      write (luncrt, '(/, 15x, a)')
     >   'n              x1              xn'
      write (luncrt, '(a, i3, 1p, 2e16.7)')
     >   ' Primary:    ', n1, x1(1), x1(n1),
     >   ' Secondary:  ', n2, x2(1), x2(n2)
      write (luncrt, '(a)', advance='no') ' Combo:      '
      read  (lunkbd, *)  n, xa, xb

      write (luncrt, '(/, a)', advance='no')
     >   ' Output file name?  '
      read  (lunkbd, '(a)') filename
      open  (lundat, file=filename, status='new')

      write (luncrt, '(a)', advance='no')
     >   ' Weight to apply to the primary dataset (e.g., 0.5): '
      read  (lunkbd, *) w1
      w2 = one - w1

      write (lundat, '(a, 2f6.3)') 'Combination coefficients: ', w1, w2
      write (lundat, '(i5)') n

      allocate (x(n), y(n), y1int(n), y2int(n), yp(n))

!     Normalize the input abscissas, and generate the target abscissas:

      call nmlizx (normalize, n1, x1, x1, x1(1), x1(n1))
      call nmlizx (normalize, n2, x2, x2, x2(2), x2(n2))

      dx = one / real (n - 1)
      x(1) = zero
      do i = 2, n - 1
         x(i) = dx * real (i - 1)
      end do
      x(n) = one

!     Interpolate the normalized curves at the target abscissas:

      call lcsfit (n1, x1, y1, new, method, n, x, y1int, yp)
      call lcsfit (n2, x2, y2, new, method, n, x, y2int, yp)

!     Form the linear combination and denormalize the abscissas:

      y = w1 * y1int + w2 * y2int

      call nmlizx (.not. normalize, n, x, x, xa, xb)

      write (lundat, '(1p, 2e16.7)') (x(i), y(i), i = 1, n)
      close (lundat)

      end program combine
