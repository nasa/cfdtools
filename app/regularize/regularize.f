********************************************************************************
*
      program regularize
*
*     PURPOSE:
*
*        Regularize all columns of a dataset by interpolating them at uniform
*        steps within the first column.  A step size is prompted for, and the
*        first line is assumed to start at or close to zero.
*
*        The intended application is to time histories from an ODE solver such
*        as a trajectory package.  Plotting symbols at every nth second on top
*        of the full curves can help interpret results.
*
*        Monotonic splines are used for the interpolations.  List-directed
*        output to a buffer deals with problems involving output record length
*        and the number of significant digits.
*
*     HISTORY:
*
*        10/02/01  DAS  Initial implementation for dealing with irregular time
*                       steps saved by a trajectory package.
*        10/03/01   "   Had to remove line feeds (preferable to formatted I/O).
*
*     AUTHOR: David Saunders, ELORET/NASA Ames Research Ctr., Moffett Field, CA
*
********************************************************************************

      implicit none

*     Constants:

      integer, parameter ::
     >   lunin = 1, lunout = 2, lunkbd = 5, luncrt = 6,
     >   maxlen = 500

*     Variables:

      integer
     >   first, i, ios, j, last, mark, ncols, nlines, nuniform

      real
     >   step

      real, allocatable ::
     >   data_in(:,:), data_out(:,:), derivs(:), uniform(:)

      character
     >   buffer * (maxlen), filein * 48, fileout * 48

*     Procedures:

      external
***  >   cleantext,  ! Utility for removing line-feeds, commas, etc.
     >   scannr,     ! Token-identifying utility used for counting columns
     >   lcsfit      ! Monotonic spline utility

*     Execution:

      write (luncrt, '(/, a)', advance='no') ' Input file:   '
      read  (lunkbd, *) filein
      open  (lunin,  file=filein,  status='old')

      write (luncrt, '(a)', advance='no') ' Output file:  '
      read  (lunkbd, *) fileout
      open  (lunout, file=fileout, status='new', recl=maxlen)

      write (luncrt, '(a)', advance='no') ' Uniform step: '
      read  (lunkbd, *) step

*     Count the number of columns in the first line (& all lines assumed):

      read  (lunin, '(a)') buffer

      last  = len_trim (buffer)
      ncols = 0
      first = 1

      do ! Until end of line
         call scannr (buffer, first, last, mark)
         if (first > 0) ncols = ncols + 1
         first = mark + 2
         if (first > last) exit
      end do

*     Count the number of input lines with a quick (?) scan:

      nlines = 1
      do ! Until EOF
         read (lunin, *, iostat=ios)
         if (ios < 0) exit
         nlines = nlines + 1
      end do

      rewind (lunin)


*     Store all input columns:

      allocate (data_in(nlines,ncols))

      do i = 1, nlines
         read  (lunin, *) data_in(i,1:ncols)
      end do

      close (lunin)

*     How many output lines?

      nuniform = int (data_in(nlines,1) / step)

      allocate (data_out(nuniform,ncols),
     >          uniform(nuniform), derivs(nuniform))

*     Target abscissas:

      do i = 1, nuniform
         uniform(i) = step * real (i)
      end do

*     Interpolate each column after the first at uniform steps in the first:

      do j = 2, ncols

         call lcsfit (nlines, data_in(1,1), data_in(1,j), .true., 'M',
     >                nuniform, uniform, data_out(1,j), derivs)
      end do

*     Save results.  The buffer works around line-length problems.

      do i = 1, nuniform
***      write (lunout, *) uniform(i), data_out(i,2:ncols) ! Too long on Intel?
         write (lunout, '(1p, 99e15.6)') uniform(i), data_out(i,2:ncols)

***      j = len_trim (buffer)

***      call cleantext (buffer, buffer, j) ! Strip out line feeds, etc.

***      write (lunout, '(a)') buffer(1:j)
      end do

      close (lunout)

      end program regularize
