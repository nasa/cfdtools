************************************************************************
*
      program extract_lines
*
*     Extract every nth line from a text file consisting of a title
*     followed by groups of n lines, until EOF.
*
*     03/13/00  D.Saunders  Initial implementation for dealing with
*               ELORET      trajectory results from HAVOC/Excel file.
*               NASA Ames
*     08/03/02      "       If i1 = 0, no title is assumed.
*
************************************************************************

      implicit none

*     Constants:

      integer, parameter ::
     >   luni = 1, luno = 2, lunkbd = 5, luncrt = 6

*     Variables:

      integer
     >   i, i1, ios, n, nchars

      character
     >   filename * 80, line * 340

*     Execution:

      write (luncrt, '(/, a)', advance='no') ' Input file?   '
      read  (lunkbd, '(a)') filename
      open  (luni, file=filename, status='old')

      write (luncrt, '(a)', advance='no')    ' Output file?  '
      read  (lunkbd, '(a)') filename
      open  (luno, file=filename, status='new')

      write (luncrt, '(a, a)') ' Any title and every nth line',
     >   ' starting from line i1 will be extracted. '
      write (luncrt, '(a)', advance='no')
     >   ' i1 = 0 means no title.  Enter n & i1: '
      read  (lunkbd, *) n, i1

      if (i1 > 0) then
         read  (luni, '(a)') line ! Title
         nchars = len_trim (line)
         write (luno, '(a)') line(1:nchars)
      end if

      do i = 1, i1 - 1 ! Skip
         read  (luni, '(a)', iostat=ios)
         if (ios < 0) exit
      end do

      do ! Until EOF

         read  (luni, '(a)', iostat=ios) line ! nth line of current group
         if (ios < 0) exit

         nchars = len_trim (line)
         write (luno, '(a)') line(1:nchars)

         do i = 1, n - 1 ! Skip
            read  (luni, '(a)', iostat=ios)
            if (ios < 0) exit
         end do
         if (ios < 0) exit

      end do

      close (luni)
      close (luno)

      end program extract_lines
