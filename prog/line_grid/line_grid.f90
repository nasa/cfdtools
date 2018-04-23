!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   program line_grid

!  Generate a 3-space grid between specified end points.  Write results either
!  as columns or for use as a PLOT3D file.
!
!  02/03/05  David Saunders  Initial implementation to generate target points at
!                            which to interpolate CFD data for comparison with
!                            wind tunnel boundary layer measurements.  Uniform
!                            distribution only, initially.
!  03/02/05    "      "      Kludged in a way of extrapolating, as needed for
!                            generating a normal from adjacent radial grid pts.
!  06/30/16    "      "      2D/3D choice now; higher-precision output.
!
!  David Saunders, ELORET/NASA Ames Research Center, Moffett Field, CA
!                  Now with AMA, Inc. at NASA ARC.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Constants:

   integer, parameter :: &
      luncrt = 6,        &
      lunkbd = 5,        &
      lunout = 1

   real, parameter :: &
      one = 1.0

!  Variables:

   integer :: &
      i, ios, n, nd

   real :: &
      dr, r, rmultiple

   real, dimension (3) :: &
      xa, xb

   real, dimension (:,:), allocatable :: &
      xi

   logical :: &
      extrap, first, rows, threed

   character :: &
      answer*1, filename*64

!  Execution:

   write (luncrt, '(/, (a))') &
      ' Generate a straight-line uniform grid segment.', &
      ' Enter ^D to terminate the following loop over line grids.'

   write (luncrt, '(/, a)', advance='no') ' 2D or 3D?  [2|3]: '
   read  (lunkbd, *) nd
   threed = nd == 3
   write (luncrt, '(/, a)', advance='no') ' Output file name: '
   read  (lunkbd, '(   a)') filename
   open  (lunout, file=filename, status='unknown')

   first = .true.
   do  ! Until ^D

      write (luncrt, '(/, a)', advance='no') ' 1st end point coordinates: '
      read  (lunkbd, *, iostat=ios) xa(1:nd)
      if (ios /= 0) exit
      write (luncrt, '(   a)', advance='no') ' 2nd end point coordinates: '
      read  (lunkbd, *, iostat=ios) xb(1:nd)
      if (ios /= 0) exit

      if (first) then

         write (luncrt, '(a)', advance='no') &
            ' Extrapolate beyond point 2? [y|n]: '
         read  (lunkbd, '(a)') answer
         extrap = answer == 'y' .or. answer == 'Y'

         if (extrap) then
            write (luncrt, '(a)', advance='no') &
               ' Multiple of distance from point 1 to point 2 [+/-]: '
            read  (lunkbd, *) rmultiple
            do i = 1, nd
               dr = xb(i) - xa(i)
               xb(i) = xa(i) + dr * rmultiple
            end do
            rows = .false.
            write (luncrt, '(a, 3es24.16)') ' Adjusted point 2: ', xb(1:nd)
         end if

         write (luncrt, '(a)', advance='no') ' # points in grid (uniform): '
         read  (lunkbd, *) n
         dr = one / real (n - 1)

         allocate (xi(nd,n))

         if (.not. extrap) then
            write (luncrt, '(a)', advance='no') ' Rows or columns? [r|c]: '
            read  (lunkbd, '(a)') answer
            rows = answer == 'r' .or. answer == 'R'
         end if

      end if

      do i = 1, n
         r = dr * real (i - 1)
         xi(:,i) = (one - r) * xa(1:nd) + r * xb(1:nd)
      end do
      xi(:,n) = xb(1:nd)  ! Make sure of it

      if (.not. rows) then
         if (threed) then
            write (lunout, '(3es24.16)') xi(:,:)
         else
            write (lunout, '(2es24.16)') xi(:,:)
         end if
      else
         write (lunout, '(5es24.16)') xi(1,:)
         write (lunout, '(5es24.16)') xi(2,:)
         if (threed) write (lunout, '(5es24.16)') xi(3,:)
      end if

      first = .false.
   end do

   close (lunout)
   deallocate (xi)

   end program line_grid
