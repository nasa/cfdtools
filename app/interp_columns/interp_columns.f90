!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      program interp_columns
!
!     PURPOSE:
!
!        Interpolate all columns of a dataset at a specified abscissa or set
!        of uniform abscissas (or arbitrary abscissas from an auxiliary file).
!        The column to treat as "x" is prompted for, and should be monotonically
!        increasing or decreasing.
!
!        Monotonic local cubic splines are used for the interpolations.
!        List-directed output to a buffer deals with problems involving output
!        record length and the number of significant digits.
!
!     HISTORY:
!
!        10/03/01  DAS  Earlier REGULARIZE program (first column = "x" = Time).
!        06/12/07   "   Adaptation for something other than a time history.
!        06/20/07   "   Added the option to read abscissas from a file.
!        04/29/11   "   Formats changed from e15.6 to e16.8 (FIAT_Opt).
!        07/08/13   "   Dinesh Prabhu needed to interpolate BLAYER output (68
!                       columns): raised the buffer length from 500 to 2048.
!
!     AUTHOR: David Saunders, ELORET Corporation/NASA Ames Research Center, CA
!                             Now with ERC, Inc. at ARC.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

!     Constants:

      integer, parameter :: &
         lunin = 1, lunin2 = 2, lunout = 3, lunkbd = 5, luncrt = 6, &
         maxlen = 2048

!     Variables:

      integer :: &
         first, i, ios, j, last, mark, ncols, nlines, nabscissas, xcol

      real :: &
         x1, x2, xinc

      real, allocatable :: &
         abscissas(:), data_in(:,:), data_out(:,:), derivs(:)

      logical :: &
         read_abscissas

      character :: &
         buffer * (maxlen), filein * 64, filein2 * 64, fileout * 64

!     Procedures:

      logical :: &
         number

      external :: &
         scannr,  &  ! Token-identifying utility used for counting columns
         lcsfit,  &  ! Monotonic spline utility
         number      ! Distinguishes a numeric token from a text token

!     Execution:

      write (luncrt, '(/, a)', advance='no') ' Input file:   '
      read  (lunkbd, *) filein
      open  (lunin,  file=filein,  status='old')

      write (luncrt, '(a)', advance='no') ' Output file:  '
      read  (lunkbd, *) fileout
      open  (lunout, file=fileout, status='unknown', recl=maxlen)

      write (luncrt, '(a)', advance='no') ' Column to treat as "x": '
      read  (lunkbd, *) xcol

      write (luncrt, '(a)', advance='no') &
         ' File containing target abscissas ("none" = uniform via a prompt): '
      read  (lunkbd, *) filein2
      read_abscissas = filein2(1:4) /= 'none'

      if (read_abscissas) then

         call read_targets ()  ! Local procedure below

         if (nabscissas <= 0) go to 99

      else  ! Uniform abscissas (or just one):

         write (luncrt, '(a)', advance='no') ' xfirst, xlast, xinc: '
         read  (lunkbd, *) x1, x2, xinc

!        How many output lines?

         if (x1 == x2) then
            nabscissas = 1
            xinc = 1.
         else if (xinc == 0.) then
            nabscissas = 2
            xinc = x2 - x1
         else
            nabscissas = 1 + nint ((x2 - x1) / xinc)
         end if

         allocate (abscissas(nabscissas))

         do i = 1, nabscissas
            abscissas(i) = x1 + xinc * real (i - 1)
         end do

      end if

!     Count the number of columns in the first line starting with a number.
!     All numeric lines are assumed to contain this number of columns, but
!     non-numeric lines may be present and will be ignored.

      do ! Until the first numeric line is found

         read  (lunin, '(a)') buffer
         last  = len_trim (buffer)
         if (last < 1) cycle

         first = 1
         call scannr (buffer, first, last, mark)

         if (number (buffer(first:mark))) exit

      end do

      ncols = 0
      first = 1

      do ! Until end of line

         call scannr (buffer, first, last, mark)
         if (first > 0) ncols = ncols + 1
         first = mark + 2
         if (first > last) exit

      end do

!     Count the number of numeric input lines:

      nlines = 1

      do ! Until EOF

         read (lunin, '(a)', iostat=ios) buffer
         if (ios < 0) exit

         last  = len_trim (buffer)
         if (last < 1) cycle

         first = 1
         call scannr (buffer, first, last, mark)

         if (number (buffer(first:mark))) nlines = nlines + 1

      end do

      rewind (lunin)

!     Store all input numeric lines:

      write (luncrt, '(/, a, 2i6)') &
         ' # numeric lines and # columns found:', nlines, ncols

      allocate (data_in(nlines,ncols))
      allocate (data_out(nabscissas,ncols), derivs(nabscissas))

      nlines = 0

      do ! Until EOF

         read (lunin, '(a)', iostat=ios) buffer
         if (ios < 0) exit

         last  = len_trim (buffer)
         if (last < 1) cycle

         first = 1
         call scannr (buffer, first, last, mark)

         if (number (buffer(first:mark))) then
            nlines = nlines + 1
            read  (buffer(1:last), *, iostat=ios) data_in(nlines,1:ncols)
            if (ios /= 0) then
               write (luncrt, '(/, a, i6)') &
                  ' Trouble reading data.  Numeric line #:', nlines
               write (luncrt, '(2a)') ' Line: ', buffer(1:last)
               go to 99
            end if
         end if

      end do

      close (lunin)

!     Interpolate each column at the specified values of the specified column:

!!!   write (luncrt, '(/, a, /, (8es12.3))') ' Abscissas:', data_in(:,xcol)
!!!   write (luncrt, '(/, a, /, (8es12.3))') ' Ordinates:', data_in(:,1)

      do j = 1, ncols

         call lcsfit (nlines, data_in(1,xcol), data_in(1,j), .true., 'M',      &
                      nabscissas, abscissas, data_out(1,j), derivs)
      end do

!     Save results.  The buffer works around line-length problems.

      if (nabscissas < 20) write (luncrt, '(a)')
      do i = 1, nabscissas
         write (lunout, '(99es16.8)') data_out(i,:)
         if (nabscissas < 20) write (luncrt, '(99es16.8)') data_out(i,:)
      end do

      close (lunout)

   99 continue

!     Local procedure for program interp_columns:

      contains

!        -----------------------------------------------------------------------
         subroutine read_targets ()
!        -----------------------------------------------------------------------

         nabscissas = 0

         open (lunin2, file=filein2, status='old', iostat=ios)
         if (ios /= 0) go to 99

         do ! Until EOF

            read (lunin2, '(a)', iostat=ios) buffer
            if (ios < 0) exit

            last = len_trim (buffer)
            if (last < 1) cycle

            first = 1

            do ! Until no more tokens on this line

               call scannr (buffer, first, last, mark)

               if (mark == 0) exit

               nabscissas = nabscissas + 1
               first = mark + 2

            end do

         end do

         if (nabscissas > 0) then

            write (luncrt, '(/, a, i6)') ' # abscissas found: ', nabscissas

            allocate (abscissas(nabscissas))

            rewind (lunin2)
            read (lunin2, *, iostat=ios) abscissas

            if (ios /= 0) then
               write (luncrt, '(/, a)')  ' Trouble reading target abscissas.'
               nabscissas = 0
            end if

         end if

         close (lunin2)

      99 return

         end subroutine read_targets

      end program interp_columns
