!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program sort_rows
!
!  Description:
!
!     This is a slightly specialized sorting utility that should still prove
!     handy in a variety of situations.  The prompting application is to a list
!     of files produced by searching with "grep" for a certain string in related
!     files (possibly in related directories).  If those files are effectively
!     numbered in a manner that does not include leading 0s for the shorter
!     numbers, then grep will not list them in true numerical order.  E.g.:
!
!        Case2/LINE-1/neqair.out:Total radiative heating =     0.381727 W/cm2
!        Case2/LINE-10/neqair.out:Total radiative heating =     0.456789 W/cm2
!        Case2/LINE-11/neqair.out:Total radiative heating =     0.123456 W/cm2
!        Case2/LINE-12/neqair.out:Total radiative heating =     0.234567 W/cm2
!             :       :      :      :
!        Case2/LINE-19/neqair.out:Total radiative heating =     0.345678 W/cm2
!        Case2/LINE-2/neqair.out:Total radiative heating =     0.765432 W/cm2
!        Case2/LINE-20/neqair.out:Total radiative heating =     0.654321 W/cm2
!             :       :      :      :
!
!     This list can be treated as a table of 6 columns.
!
!     The desired order is most likely to be:
!
!        Case2/LINE-1/neqair.out:Total radiative heating =     0.381727 W/cm2
!        Case2/LINE-2/neqair.out:Total radiative heating =     0.765432 W/cm2
!             :       :      :      :
!        Case2/LINE-9/neqair.out:Total radiative heating =     0.345678 W/cm2
!        Case2/LINE-10/neqair.out:Total radiative heating =     0.456789 W/cm2
!        Case2/LINE-12/neqair.out:Total radiative heating =     0.234567 W/cm2
!             :       :      :      :
!
!     To achieve this order, the user would enter LINE- as the prefix prompted
!     for.  The utility will then read the full list as a table and search row 1
!     for the prefix in some column.  That column is then used to reorder the
!     rows as suggested above.  Not finding the prefix on row 1 is considered an
!     error, as is not finding it in the same column for some other row.
!
!     Existing subroutine get_coordinates is reused to isolate the number that
!     follows the prefix in each row.  Those numbers are sorted via an index
!     list, and the rows are written in the sorted order.  The numbers are most
!     likely to be integers, but reals are also allowed for, such as 12.345.
!
!     (Later:)  Since get_coordinates searches for numbers from the right, the
!     original implementation could not sort something like this ...
!
!        ant1_13.attenuation_v3.30_h38.00.dat
!        ant1_13.attenuation_v3.60_h39.00.dat
!        ant1_13.attenuation_v3.80_h40.00.dat
!        ant1_13.attenuation_v4.00_h41.00.dat
!        ant1_13.attenuation_v4.30_h42.00.dat
!        ant1_13.attenuation_v3.26_h42.08.dat
!        ant1_13.attenuation_v4.70_h43.00.dat
!        ant1_13.attenuation_v3.70_h43.00.dat
!
!     ... except in the order of h, while sorting in order of v may be desired.
!     This has been handled now.
!
!  History:
!
!     02/10/17  D.A.Saunders  Initial design and implementation minus sorting.
!     02/13/17    "      "    Completed the sorting part after extending the
!                             table_io module so that the "read as alpha" option
!                             can work with input lines that are all alpha-
!                             numeric (meaning no leading header lines).
!     03/16/17    "      "    Handled more than one "coordinate" (number) in the
!                             column containing the prefix.
!
!  Author:  David Saunders, AMA, Inc./NASA Ames Research Center, Mnt. View, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use table_type_module  ! Both are in table_io.f90
   use table_io

   implicit none

!  Constants:

   integer, parameter :: &
      lunin  = 1,        &  ! Input file list/table
      lunout = 2,        &  ! Reordered output list
      lunkbd = 5,        &  ! Keyboard inputs
      luncrt = 6,        &  ! Screen prompts and diagnostics
      maxbuf = 512,      &  ! Should match table_io module
      ndim   = 1            ! One number per row is looked for

!  Variables:

   integer :: i, i1, i2, icol, inc, ios, irow, l, lpre, ncols, nrows

   real :: coords(ndim)  ! An array argument is expected by get_coordinates
   logical :: cr, eof    ! Handle carriage return or ctrl D
   logical :: descending
   character (maxbuf)   :: buffer
   character (132)      :: prefix

   integer, allocatable :: ilist(:)
   real,    allocatable :: rlist(:)

   type (table_type) :: table

!  Execution:

   call reads (luncrt, 'List or table to be reordered: ', &
               lunkbd, table%filename, cr, eof)
   if (cr .or. eof) go to 99

!  Read the table as character tokens:

   call table_io_read_alpha (lunin, table, ios)
   if (ios /= 0) go to 99

   nrows = table%nrows
   ncols = table%ncols

   write (luncrt, '(a, 2i6)') ' # rows and columns found:', nrows, ncols

   call reads (luncrt, 'Prefix for numbers to be reordered: ', &
               lunkbd, prefix, cr, eof)
   if (cr .or. eof) go to 99

   lpre = len_trim (prefix)

!  Locate the column containing the prefix:

   do icol = 1, ncols
      l = len_trim (table%tokens(icol,1))
      i = index (table%tokens(icol,1)(1:l), prefix(1:lpre))
      if (i > 0) then
         l = icol
         exit
      end if
   end do

   if (i <= 0) then
      write (luncrt, '(2a)') ' No such prefix found: ', prefix(1:lpre)
      go to 99
   end if

   icol = l
   write (luncrt, '(a, i5)') ' Column containing this prefix:', icol

   descending = .false.
   i1  = 1
   i2  = nrows
   inc = 1
   call ready (luncrt, 'Output in descending order? [y|n; <cr>=n] ', &
               lunkbd, descending, cr, eof)
   if (eof) go to 99

   if (descending) then
      i1  = nrows
      i2  = 1
      inc = -1
   end if

   allocate (rlist(nrows))

!  Search this col. of each row for the prefix and decode the associated number:

   do irow = 1, nrows
      l = len_trim (table%tokens(icol,irow))
      call skip_any_others ()  ! Since get_coordinates searches from the right
      call get_coordinates (table%tokens(icol,irow)(1:l), ndim, coords, ios)
      if (ios /= 0) then
         write (luncrt, '(a, i6)') &
            ' Number missing from this column; row #:', irow
         go to 99
      end if
      rlist(irow) = coords(1)
   end do

!  Set up the output file:

   call reads (luncrt, 'File name for sorted table: ', &
               lunkbd, table%filename, cr, eof)
   if (cr .or. eof) go to 99

   open (lunout, file=table%filename, status='unknown', iostat=ios)
   if (ios /= 0) then
      l = len_trim (table%filename)
      write (luncrt, '(2a)') 'Trouble opening ', table%filename(1:l)
   end if

!  Sort the list in-place, and an index list likewise:

   allocate (ilist(nrows))
   do i = 1, nrows
      ilist(i) = i
   end do

   call hsortri (rlist, ilist, nrows)

!  Write the rows in the sorted order.  We can't use table_io_write_alpha.

   do i = i1, i2, inc
      irow = ilist(i)
      call pack_row ()  ! Internal procedure below
   end do

   close (lunout)

99 continue

!  Internal subroutines for program sort_rows:

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine skip_any_others ()

!     Handle a string such as ant1_13.attenuation_v3.30_h38.00.dat where v, not
!     h, is the specified prefix, given that get_coordinates searches from the
!     right and would return 38.0 if the full string were treated.  In this
!     example, l is 36 upon input, and will be returned as 25 if 'v' or '_v' is
!     the intended prefix (or 32 for 'h' or '_h').

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: i, j, last
      logical, external :: number  ! Utility for distinguishing digits

      i = index (table%tokens(icol,irow)(1:l), prefix(1:lpre))  ! Like above
      last = l
      do j = i + lpre, last
         if (.not. number (table%tokens(icol,irow)(j:j))) then
            if (table%tokens(icol,irow)(j:j) /= '.') then  ! Nasty one
               l = j - 1
               exit
            end if
         end if
      end do

      end subroutine skip_any_others

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine pack_row ()  ! Pack the tokens for irow and write the row.
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: ibuf

      ibuf = 1
      do icol = 1, ncols
         l = len_trim (table%tokens(icol,irow))
         buffer(ibuf:ibuf+l-1) = table%tokens(icol,irow)(1:l)
         ibuf = ibuf + l
         buffer(ibuf:ibuf) = ' '
         ibuf = ibuf + 1
      end do
      write (lunout, '(a)') buffer(1:ibuf-2)

      end subroutine pack_row

   end program sort_rows
