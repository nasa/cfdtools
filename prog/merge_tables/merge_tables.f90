!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program merge_tables
!
!  Outline:
!
!     From any number of tables of equal length but variable numbers of columns,
!     produce a new table by combining the indicated columns in the specified
!     order.  The intent is to be able to automate a lot of cutting and pasting.
!     The columns of a single table may also be reordered.
!
!  Control Format (illustrated with an example that merges parts of 3 tables):
!
!     table1.dat
!     table2.dat
!     table3.dat
!     1  1  1  2  2  3
!     2  4  3  8 12  1
!     newtable.dat
!
!  Usage:
!
!     merge_tables < merge.inp  (say)
!
!  Notes:
!
!     o  If column headers are present in the first 1 or more lines, those lines
!        should contain as many tokens as there are (fully) numeric columns.
!        Otherwise, column headers cannot be carried along to the new table.
!     o  If this is true, the entire table could be non-numeric, but expected
!        usage is on numeric tables with column headers.
!     o  If it is not true, initially no attempt is made to ignore headers.
!        Therefore all lines of all tables are actually treated as text, not
!        partly text and partly numeric.
!     o  The first indefinite list of integers consists of table numbers.
!     o  The second indefinite list (same length) contains column numbers.
!        These two lists define the output columns.
!
!  Sample Input Tables:
!
!     o  The idea is to combine certain columns from tables that could have
!        been extracted from (say) a trajectory time history as CFD points
!        analyzed by a flow solver, extracted from the resulting boundary
!        layer datasets (see EXTRACT_BLAYER_DATA), and "grep"ed from (say)
!        a radiative heating solver.  The result in this case is another
!        table from which certain columns can be curve fit to produce full
!        trajectory time histories at some body point for a material reponse
!        solver.  E.g., a "CFD point" table might look like this:
!
!    t,s   V,km/s p,Pa    rho,kg/m^3 T,K    pstag,Pa M  qconv,W/cm^2 qrad,W/cm^2
!    47.00 11.944 3.90e+00 1.18e-04  171.84 1.64e+04 56.4  395.8  311.5
!    50.00 11.740 1.60e+01 4.94e-04  169.06 6.61e+04 55.8  788.5 1268.2
!    52.10 11.309 4.12e+01 1.25e-03  172.21 1.55e+05 53.3 1126.1 1967.8
!    54.10 10.455 9.28e+01 2.70e-03  179.59 2.87e+05 48.3 1283.9 1208.2
!    56.60  8.679 2.15e+02 5.90e-03  190.54 4.33e+05 39.0  992.2  465.8
!    58.70  6.869 3.69e+02 9.72e-03  198.59 4.45e+05 30.3  508.9  133.5
!    61.00  5.064 5.70e+02 1.45e-02  206.07 3.62e+05 22.0  248.2   19.4
!
!     o  Note the single-token column headers.  There may be more than one
!        header line, but all tables being merged must contain the same total
!        number of lines.  Header and numeric lines are all treated the same
!        way, as text, and the numbers of each should match across tables.
!
!     o  Simply numbering the unwanted columns, as from a "grep", may serve
!        this purpose.  For instance, only the radiative heat flux column 5
!        in the following example needs to have a meaningful header:
!
!    1                                 2         3       4 qrad,W/cm^2  6
!    G12-t47.0/LINE-1/neqair.out:Total radiative heating = 2.099655E-01 W/cm2
!    G12-t50.0/LINE-1/neqair.out:Total radiative heating = 6.122821E-01 W/cm2
!    G12-t52.1/LINE-1/neqair.out:Total radiative heating = 1.643167E+00 W/cm2
!    G12-t54.1/LINE-1/neqair.out:Total radiative heating = 2.978185E+00 W/cm2
!    G12-t56.6/LINE-1/neqair.out:Total radiative heating = 4.252228E+00 W/cm2
!    G12-t58.7/LINE-1/neqair.out:Total radiative heating = 5.858140E+00 W/cm2
!    G12-t61.0/LINE-1/neqair.out:Total radiative heating = 6.259667E+00 W/cm2
!
!  Merging Strategy:
!
!     o  Memory these days is plentiful, so we can work with a large limit on
!        the maximum size of any token, and we can store the input data as a
!        3-D array of all tokens from all columns of all tables.  This greatly
!        simplifies the merging, at the expense of allocating as much storage
!        for each table as is needed for the largest table.
!     o  If all tables are not of equal length (same number of rows), the
!        program aborts.  The same is true  if any input table does not have
!        the same number of columns in all of its rows.
!     o  Fortran 90 does not appear capable of allocating arrays of variable-
!        length strings.  This is why some upper limit on the maximum length of
!        a token is necessary.  (It can be larger than is ever likely needed.)
!     o  Upon allocation of enough storage, each table is reread line by line,
!        with each line tokenized as we go.  All tokens are treated as text.
!        There is no need to decode them as integers or reals.
!     o  The rest is now easy.  The output columns are no wider than the widest
!        token in each columm, separated by a single blank.
!
!  History:
!
!     08/24/2014  D.A.Saunders  Initial design, as part of automating the
!                               preparation of flow solver data and radiation
!                               solver data in the form expected by a material
!                               response solver.  At NASA ARC, these are DPLR,
!                               NEQAIR, and FIAT respectively.  A trajectory
!                               solver (Traj at ARC) is also in the picture
!                               if time histories are being treated.
!     09/05/2014    "     "     Completed initial coding after a hiatus.
!     09/06/2014    "     "     TOKEN4 mishandles the = character here.
!                               Introduced a simpler TOKEN2L variant of TOKEN2
!                               to avoid the latter's uppercasing.
!     10/03/2014    "     "     Considered making use of the (later) table_io
!                               module, but its separate handling of header
!                               lines is problematic.  Requiring single-token
!                               headers matching the columns allows the column
!                               headers to be carried along into the merged
!                               table, while there is no good solution for
!                               arbitrary numbers of arbitrary header lines
!                               as handled by table_io, with its completely
!                               different data structure.  The table format
!                               expected here has been clarified with examples.
!                               See PREPARE_FIAT_DATA for a different table
!                               application that uses table_io effectively.
!
!  Author:  David Saunders, ERC, Inc./NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Local constants:

   integer, parameter :: &
      lundat    = 1,     &    ! For input tables
      lunout    = 2,     &    ! For the output table
      lunctl    = 5,     &    ! For control inputs (standard in)
      luncrt    = 6,     &    ! For diagnostics
      lunoffset = 10,    &    ! For an indefinite list of logical units
      maxbuffer = 1024,  &    ! Limit on the length of any table line, in or out
      maxcols   = 256,   &    ! Limit on output # columns (the integer lists)
      maxtables = 16,    &    ! Limit on # input tables (for logical unit array)
      maxtoken  = 64,    &    ! Limit on the length of any token
      maxtokens = 128         ! Limit on # columns in an input table

   character (1) :: &
      blank = ' '

!  Variables:

   integer :: &
      i, ios, itable, j, k, l, ncols, ncols2, ncolsin, ncolsout, &
      nrows, nrowsall, ntable
   integer :: &
      list_columns(maxcols), list_tables(maxcols), luntable(maxtables)
   character (2) :: &
      seps
   character (maxbuffer) :: &
      buffer
   character (maxtoken), allocatable :: &
      tokens(:,:,:)

!  Procedures:

   logical  :: alpha
   external :: alpha  ! Distinguishes text from numeric input

!  Execution:

   seps = blank // char (9)  ! Blank | tab may delimit table columns (not comma)

   call read_controls_and_scan_data ()  ! Internal procedure below
   if (ios /= 0) go to 99

   allocate (tokens(ncolsin,nrows,ntable))  ! For all tokens of all tables

   call read_all_tables ()     ! Rewind each table file and read all tokens

   call write_merged_table ()  ! ... and wrap up

99 continue

!  Internal procedures for program merge_tables:

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine read_controls_and_scan_data ()
!
!     The control file is processed by reading file names and opening each file
!     on its own unit number until an integer list is encountered.  This list of
!     table numbers is decoded, as is the following list of table columns.
!     Finally, the output file name is read and the file is opened.
!     Error checking such as detecting tables of different length is also done,
!     and the largest number of columns in any input table is identified.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ntable = 0  ! # input tables found

      do  ! Scan tables until a fully numeric line is found (first integer list)

         read (lunctl, '(a)', iostat=ios) buffer
         if (ios /= 0) exit   ! Unexpected EOF

         l = len_trim (buffer)
         write (luncrt, '(3x, a)') buffer(1:l)

         if (alpha (buffer(1:l))) then      ! Expect a table file name
            ntable = ntable + 1
            buffer = adjustl (buffer(1:l))  ! Eliminate any leading blank(s)
            l = len_trim (buffer(1:l))
            luntable(ntable) = lunoffset + ntable

            open (luntable(ntable), file=buffer(1:l), status='old', iostat=ios)
            if (ios /= 0) then
               write (luncrt, '(2a)') '*** Table not found: ', buffer(1:l)
               exit
            end if

            call scan_table ()  ! Internal procedure below
            if (ios /= 0) exit

            if (ntable == 1) then
               nrowsall = nrows
               ncolsin  = ncols
            else
               if (nrows /= nrowsall) then
                  write (luncrt, '(3a, 2i6)') &
               '*** Table ', buffer(1:l), ' row count differs:', nrows, nrowsall
                  ios = 1
                  exit
               end if
               ncolsin = max (ncolsin, ncols)
            end if

         else  ! Expect two equal-length lists of integers ...

            ncolsout = maxcols  ! Input with limit; output with # found

            call decode_ilist (buffer(1:l), ncolsout, list_tables, ios)

            if (ios /= 0) then
               write (luncrt, '(a)') &
                  '*** Apparently a bad index in this control line:', &
                  buffer(1:l)
               exit
            end if

            read (lunctl, '(a)', iostat=ios) buffer
            if (ios /= 0) exit

            l = len_trim (buffer)
            write (luncrt, '(3x, a)') buffer(1:l)

            ncols2 = maxcols

            call decode_ilist (buffer(1:l), ncols2, list_columns, ios)

            if (ios /= 0) then
               write (luncrt, '(a)') &
                  '*** Apparently a bad index in this control line:', &
                  buffer(1:l)
               exit
            end if

            if (ncols2 /= ncolsout) then
               write (luncrt, '(a, 2i5)') &
                  '*** Column list does not match table list:', ncolsout, ncols2
               ios = 1
               exit
            end if

!           ... followed by the output file name:

            read (lunctl, '(a)', iostat=ios) buffer
            if (ios /= 0) exit

            l = len_trim (buffer)
            write (luncrt, '(3x, a)') buffer(1:l)

            open (lunout, file=buffer(1:l), status='unknown', iostat=ios)

            if (ios /= 0) then
               write (luncrt, '(2a)') &
                  '*** Bad output file name? ', buffer(1:l)
            end if

            exit

         end if

      end do  ! End of control file processing

      if (ios == 0) then
         write (luncrt, '(/, (3x, a, i6))') &
            '# tables found:           ', ntable,   &
            '# rows in all tables:     ', nrows,    &
            '# columns in widest table:', ncolsin,  &
            '# columns in output table:', ncolsout
      end if

      end subroutine read_controls_and_scan_data

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine scan_table ()
!
!     Count the number of rows and columns in the current table.
!     Abort if the table does not have a fixed number of columns in each row.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local variables:

      integer :: ncolsallrows

!     Execution:

      nrows = 0

      do  ! Until EOF
         read (luntable(ntable), '(a)', iostat=ios) buffer
         if (ios < 0) exit

         nrows = nrows + 1

         call token_count (trim (buffer), seps, ncols)

         if (nrows == 1) then
            ncolsallrows = ncols
         else
            if (ncols /= ncolsallrows) then
               write (luncrt, '(a, i3, a, 2i4)') &
                  '***Table', ntable, ' has variable numbers of columns:', &
                  ncolsallrows, ncols
               ios = 1
               go to 99
            end if
         end if
      end do

      ios = 0
 99   return

      end subroutine scan_table

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine read_all_tables ()  ! Rewind each table file & read all tokens
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local variables:

      integer :: j, lun, n, ncols

!     Execution:

!     Use TOKEN2L rather than TOKEN2 to avoid TOKEN2's use of UPCASE.

      do n = 1, ntable
         lun = luntable(n)
         rewind (lun)
         do j = 1, nrows
            read (lun, '(a)') buffer
            ncols = ncolsin   ! At most
            call token2l (trim (buffer), seps, ncols, tokens(:,j,n))
         end do
         close (lun)
      end do

      end subroutine read_all_tables

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine write_merged_table ()  ! ... and wrap up
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local variables:

      integer :: i, i1, i2, ier, j, l, m, n, width
      integer, allocatable :: colwidth(:)

!!!!  allocate (tokens(ncolsin,nrows,ntable))  ! For reference

      allocate (colwidth(ncolsout))

!     Determine a maximum width for each output column:

      do i = 1, ncolsout
         width = 0
         n = list_tables(i)   ! Table # for output column i
         m = list_columns(i)  ! Column from table n for output column i
         do j = 1, nrows
            width = max (width, len_trim (tokens(m,j,n)))
         end do
         colwidth(i) = width
      end do

      write (luncrt, '(a, 128i3)') '   Output column widths:', colwidth(:)

!     Filling a buffer for each line is easier than constructing the right
!     format for writing all the columns of a line.

      ier = 0

      do j = 1, nrows
         buffer = blank
         i1 = 1  ! Left-justify to avoid any leading blanks that Excel
         i2 = 0  ! seems to treat as an extra column when pasting

         do i = 1, ncolsout
            n = list_tables(i)
            m = list_columns(i)
            l = colwidth(i)
            i2 = i1 + l
            if (i2 > maxbuffer) then
               if (ier == 0) then
                   ier = 1
                   write (luncrt, '(a, i4)') &
                     '*** Buffer too short for output table. Trouble column:', i
               end if
            else
               buffer(i1:i2) = tokens(m,j,n)(1:l+1)  ! At least 1 blank between
            end if
            i1 = i2 + 1
         end do

         if (ier /= 0) then
            write (luncrt, '(a, 2i5)') &
               '    Current and needed buffer length:', maxbuffer, i2 - 1
            exit
         end if

         write (lunout, '(a)') buffer(1:i2-1)     ! ... but not on the end
      end do

      deallocate (tokens, colwidth)
      close (lunout)

      end subroutine write_merged_table

   end program merge_tables
