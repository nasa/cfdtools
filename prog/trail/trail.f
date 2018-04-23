************************************************************************
*
      program trail
*
*     Strip trailing blanks from a text file.
*
*     12/08/00  David Saunders  Initial implementation for Unix systems,
*                               in anticipation of losing XDECK on a VMS
*                               system where file version numbers are so
*                               convenient.
*     04/19/06    "      "      Flag lines with non-blanks beyond the
*                               specified column (probably 80 or 132).
*
************************************************************************

      implicit none

*     Constants:

      integer, parameter ::
     >   lunkbd = 5, luncrt = 6, lunin = 7, lunout = 8

*     Variables:

      integer
     >  icol, iline, ios, last

      character
     >  filename * 64, line * 132

*     Execution:

      ios = 1
      do while (ios /= 0)
         write (luncrt, '(/, a)', advance = 'no') ' Input file:  '
         read  (lunkbd, '(a)') filename
         open  (lunin, file = filename, status = 'old', iostat = ios)
      end do

      ios = 1
      do while (ios /= 0)
         write (luncrt, '(a)', advance = 'no') ' Output file: '
         read  (lunkbd, '(a)') filename
         open  (lunout, file = filename, status = 'new', iostat = ios)
      end do

      write (luncrt, '(a)', advance='no')
     >   ' Flag lines extending beyond which column? [80|132|999] '
      read (lunkbd, *) icol

      iline = 0
      do ! Until EOF
         read  (lunin, '(a)', iostat = ios) line
         if (ios < 0) exit

         iline = iline + 1
         last  = len_trim (line)

         if (last > icol)
     >      write (luncrt, '(a, i7)') ' Too long at line #', iline

         write (lunout, '(a)') line(1:last)
      end do

      end program trail
