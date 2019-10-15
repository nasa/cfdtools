!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   module table_type_module

!  This module defines a derived data type representing columns of numerical data with optional header line(s) -- a text table.
!  Header lines are identified as being alphanumeric, while non-header lines are not.  All numeric lines are expected to have the
!  same number of tokens.  Embedded header lines are presently disallowed, but the benefit of modularizing this sort of thing is
!  that desirable extensions should be able to be made as the needs arise.
!
!  Table contents may initially be stored as character data.  For some applications, such as merging of tables, the numerical values
!  may not be needed.  It may be worth marking each column as integer or real or neither (= character) and treat them appropriately
!  at the higher level.  The most common case of real-number tables has been addressed first, with a convenient data structure and
!  the associated I/O utilities.  With the run-time formatting option, integer columns may now be stored as reals and written as
!  integers.  This is implemented by looking for 'i' in table_io_format if table_io_format_sep is passed to table_io_write_real as
!  its separator argument.

!  Do these max* limits need to be public?  That is unclear at the time of writing.  Choose names that should allow it.

   integer,       parameter :: max_table_file_name     = 192  ! Limit on the length of the table filename
   integer,       parameter :: max_table_buffer        = 512  ! Limit on the length of a numeric line or row in the table
   integer,       parameter :: max_table_format        = 128  ! Limit on the length of run-time format string table_io_format
   integer,       parameter :: max_table_header_length = 512  ! Limit on the length of a header line
   integer,       parameter :: max_table_token_length  = 64   ! Limit on the length of a numeric token (width of a numeric column)
   integer,       parameter :: max_table_columns       = 192  ! Limit on the number of numeric columns in a table
   character (1), parameter :: table_io_format_sep     = 'F'  ! Pass this separator to table_io_write_* to use run-time formatting

   character (max_table_buffer) :: table_io_buffer      ! All table lines are buffered by the accompanying I/O utilities
   character (max_table_format) :: table_io_format      ! Run-time output format like '(10es20.12)' invoked via table_io_format_sep

   type table_type
      integer                                      :: nheader       ! # header lines; nheader >= 0
      integer                                      :: nrows         ! # numeric lines or rows; nrows >= 1
      integer                                      :: ncols         ! # numeric columns >= 1; each numeric row has ncols columns
      integer, dimension (max_table_columns)       :: column_type   ! 1 = integer, 2 = real, 3 = character
      integer, dimension (max_table_columns)       :: column_width  ! Maximum width of each input column; may have other uses
      real,                                pointer :: values (:,:)  ! For the common case of all-real-value tables
      character (max_table_file_name)              :: filename      ! Name of the file to use if the table is to be read or written
      character (max_table_header_length), pointer :: header (:)    ! Header lines to be allocated as needed
      character (max_table_token_length),  pointer :: tokens (:,:)  ! Character tokens(i,j) = token or column i on non-header row j
   end type table_type

!  Storage Note:  These are line-oriented utilities, so the tokens on a line/row are stored contiguously (i = 1:ncols for row j).
!                 Conversely, matrix notation A(i,j) means row i, column j, so COLUMN elements are contiguous in Fortran.

   end module table_type_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   module table_io

!  This module packages I/O utilities for common data tables such as trajectory time histories.  Depending on the application, a
!  given table may be treated as numeric (with optional header text line(s)) or as alphanumeric text.  Either way, numeric lines
!  are distinguished from header lines by containing nothing that could not be part of a numerical value.
!
!  If the dataset is strictly alphanumeric, and read as such, it is checked for having all rows matching in terms of number of
!  columns, and if so it is treated as having no header lines.
!
!  Restrictions:
!
!     o  Header lines must be at the top of a table: embedded header lines and trailing comments are presently treated as errors.
!     o  A valid table must have the same number of tokens (items) on each non-header line.
!     o  Numeric column separators may be blank, tab or comma characters, or any other sensible character.
!     o  The alphanumeric case should NOT contain commas if they are not intended to be column delimiters.
!
!  Other Notes:
!
!     o  When written as character tokens (table_io_write_alpha) with separator = blank, table columns are left-justified and as
!        wide as the widest token in the column and separated by at least a single space.
!     o  When written as real values (table_io_write_real) with "sep" /= blank, some blanks may accompany the indicated separator.
!     o  Both of these options may use a run-time-variable format, invoked by passing "sep" as table_io_format_sep along with a
!        sensible format string assigned to table_io_format.
!     o  If the table_io_format string contains any 'i' character(s), and table_io_format_sep is passed to table_io_write_real,
!        then the indicated column(s) are written as integers although they are stored as reals.
!     o  Thorny problem:  how to allow commenting out data lines?
!     o  What else?
!
!  History:
!
!     09/16/14-  David Saunders  Initial implementation following work that involved too much cutting and pasting into Excel.
!     09/18/14                   Quite a few issues arise, though, so only the basic requirements are implemented to begin with.
!     10/02/14     "      "      Provided for writing some real columns as integers if i' or 'I' appear in table_io_format.
!     01/20/15     "      "      Awkwardness:  if at least one column, stored as real, is desired to be written as integer via a
!                                format like 'i4', the strategy implemented here means any accompanying column intended to end
!                                with a '.' such as f6.0 will also be adjusted and appear as i5 in this example.  Compromise:
!                                keep a line-ending '.' not associated with an integer format, but any interior 'f6.0'-type
!                                columns will still be printed as integers, which may not or may not be an issue.
!                                Also: an option to thin a table in-place has been added.
!     11/18/16     "      "      A table with 369 characters required some larger limits above.
!     02/11/17     "      "      The alphanumeric case (table_io_read_alpha) required allowing for the case of all alphanumeric
!                                that match in terms of number of columns, with no header records assumed.
!     11/20/18     "      "      Displaying a (real) table on the screen required suppression of the open and close in subroutine
!                                table_io_write_real.  Do the same in table_io_write_alpha.  While we're at it, add table_io_copy.
!
!  Author:  David Saunders, ERC, Inc.,/NASA Ames Research Center, Moffett Field, CA (now with AMA, Inc. at NASA ARC).

!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use table_type_module  ! Defines various string limits and the "table_type" derived data type for a table

   implicit none

   private

!  Constants used by the package:

   logical,   parameter :: false = .false.,        &
                           true  = .true.
   character (1), parameter :: blank = ' '
   character (3), parameter :: seps  = blank // char (9) // ',' ! Blank or tab or comma

!  Internal variables used by the package:

   integer   :: inumber, line, nheader, ncols, nrows            ! Line counter, and internal copies of table% components
   real      :: rnumber

!  Parsing utilities used by the package:

   logical   :: alpha
   external  :: alpha, scan2l, token_count

!  Utilities provided for applications:

   public :: table_io_scan               ! Determines table dimensions as may be needed in an application for other variables
   public :: table_io_read_alpha         ! Reads  a table for which all numeric columns are treated as alphanumeric text
   public :: table_io_write_alpha        ! Writes a table for which all numeric columns are treated as alphanumeric text
   public :: table_io_read_real          ! Reads  a table for which all numeric columns are treated as reals
   public :: table_io_write_real         ! Writes a table for which all numeric columns are stored  as reals; may write as integer
!! public :: table_io_read_mixed         ! Reads  a table for which the numeric column data types are defined by column_type(:)
!! public :: table_io_write_mixed        ! Writes a table for which the numeric column data types are defined by column_type(:)
   public :: table_io_thin               ! Repack a table with row 1 + every nth following row, as may be needed for a time history
   public :: table_io_deallocate         ! Deallocates table arrays
   public :: table_io_copy               ! Copies a table to a new or existing table

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine table_io_scan (lun, table, ios)

!     Open the indicated table on the indicated unit, determine its number of header lines, whether all numeric rows have the same
!     number of entries, and [some day?] which columns appear to be integers and which appear to be reals (although the way they are
!     treated may be up to the application).  Trap what appears to be a defective table.
!
!     The table is left open and rewound upon return if it appears good.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer,           intent (in)         :: lun      ! Logical unit for the file being read;
                                                         ! opened here upon entry; rewound here before the return
      type (table_type), intent (inout)      :: table    ! Data structure for a numerical table with optional header line(s),
                                                         ! input with filename set; output with nheader, nrows, ncols (only) set
      integer,           intent (out)        :: ios      ! 0 means no error was detected;
                                                         ! 1 means the table appears to be invalid (embedded header line, or
                                                         ! not a constant # items per row)
!     Local variables:

      integer :: lbuf, lname, ncolsallrows, nlines

!     Execution:

      lname = len_trim (table%filename);  open (lun, file=table%filename(1:lname), status='old', iostat=ios)

      if (ios /= 0) then
         write (*, '(2a)') '*** Unable to open ', table%filename(1:lname)
         go to 99
      end if

      nheader = 0;  nrows = 0;  nlines = 0

      do  ! Until EOF
         read (lun, '(a)', iostat=ios) table_io_buffer

         if (ios < 0) then
             ios = 0
             rewind (lun)
             exit
         end if

         lbuf = len_trim (table_io_buffer)
         nlines = nlines + 1

         if (alpha (table_io_buffer(1:lbuf))) then  ! Appears to be a header line
            if (nrows == 0) then
               nheader = nheader + 1
            else  ! What about commented-out lines or trailing comments??  Don't deal with those yet.
               write (*, '(3a, i7)') '*** Non-numeric data apparently embedded in table ', table%filename(1:lname), ', line', nlines
               ios = 1
               go to 99
            end if
         else  ! Strictly numeric line: count how many columns it has

            nrows = nrows + 1

            call token_count (table_io_buffer(1:lbuf), seps, ncols)

            if (nrows == 1) then
                ncolsallrows = ncols
            else
               if (ncols /= ncolsallrows) then
                  write (*, '(3a, 2i5)') '*** Table ', table%filename(1:lname), ' has variable numbers of columns: ', &
                    ncolsallrows, ncols, '    See', blank, 'line ', nlines
                  ios = 1
                  go to 99
               end if
            end if
         end if
      end do

      table%nheader = nheader
      table%nrows   = nrows
      table%ncols   = ncols

 99   continue

      end subroutine table_io_scan

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine table_io_rescan (lun, table, ios)

!     This variant of the original scanning routine is prompted by the case of using table_io_read_alpha on a table that appears to
!     be all header records because the columns of data aren't necessarily strictly numeric (but are to be read as alphanumeric).
!     There may be a need to work with alphanumeric columns, as in the case of the author's SORT_ROWS utility.
!
!     The strategy (in read_alpha) is to interpret nheader > 0 and nrows = 0 to mean there are no header lines and the data rows
!     are alphanumeric, not strictly numeric.  We still need to check that all rows have the same number of columns.
!
!     Rewind the table, count the tokens on line 1, and make sure the remaining lines contain the same number of tokens.  If so,
!     rewind again and return for the read_alpha calling routine to continue.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer,           intent (in)         :: lun      ! Logical unit for the file being read;
                                                         ! opened here upon entry; rewound here before the return
      type (table_type), intent (inout)      :: table    ! Data structure for what is now considered to be a no-header alpha table;
                                                         ! output with nheader = 0, and nrows, ncols set
      integer,           intent (out)        :: ios      ! 0 means no error was detected;
                                                         ! 1 means the table appears to be invalid (not a constant # items per row)
!     Local variables:

      integer :: lbuf, ncolsallrows

!     Execution:

      nrows = 0

      do  ! Until EOF
         read (lun, '(a)', iostat=ios) table_io_buffer

         if (ios < 0) then
             ios = 0
             rewind (lun)
             exit
         end if

         lbuf  = len_trim (table_io_buffer)
         nrows = nrows + 1

         call token_count (table_io_buffer(1:lbuf), seps, ncols)

         if (nrows == 1) then
             ncolsallrows = ncols
         else
            if (ncols /= ncolsallrows) then
               write (*, '(3a, 2i5)') '*** Table ', trim (table%filename), ' has variable numbers of columns: ', &
                 ncolsallrows, ncols, '    See', blank, 'line ', nrows
               ios = 1
               go to 99
            end if
         end if
      end do

      table%nheader = 0
      table%nrows   = nrows
      table%ncols   = ncols

 99   continue

      end subroutine table_io_rescan

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine table_io_read_alpha (lun, table, ios)  ! Scan a table then (if all is well) read it as character tokens (only).

!     Afterthought:  The table data need not be numeric in this case, but it will look like all header lines if (all) the rows are
!     alphanumeric.  Therefore, if the initial scan finds nheader > 0 and nrows = 0, rescan the able for the possibility of all
!     rows being alphanumeric with equal numbers of tokens, and no header lines.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer,           intent (in)    :: lun    ! The table is opened here, read as header and tokens, and closed upon return
      type (table_type), intent (inout) :: table  ! %filename should be set on input;
                                                  ! %header(:) and %tokens(:,:) are set upon return
      integer,           intent (out)   :: ios    ! 0 means no problem

!     Execution:

      call table_io_scan (lun, table, ios)  ! Determine the table's number of header lines, rows, and columns
      if (ios /= 0) then
         write (*, '(2a)') '*** Trouble scanning table ', trim (table%filename)
         go to 99
      end if

      nheader = table%nheader
      if (nheader > 0 .and. nrows /= 0) then
         allocate (table%header(nheader))
         do line = 1, nheader
            read (lun, '(a)', iostat=ios) table%header(line)  ! Separate read for each header line
            if (ios /= 0) then
               write (*, '(3a, i4)') '*** Trouble reading table ', trim (table%filename), ' at header line', line
               go to 99
            end if
         end do
      else if (nheader > 0 .and. nrows == 0) then  ! Retrofitted option described above
         call table_io_rescan (lun, table, ios)    ! Check that all rows have the same number of columns, and assume nheader = 0
         if (ios /= 0) then
            write (*, '(2a)') '*** Trouble rescanning table ', trim (table%filename)
            go to 99
         end if
      end if

      ncols = table%ncols  ! No need to initialize it for every line; table_io_[re]scan should already have checked all lines
      nrows = table%nrows
      allocate (table%tokens(ncols,nrows), stat=ios)
      if (ios /= 0) then
         write (*, '(3a, i6)') '*** Trouble allocating table ', trim (table%filename), '.  nrows, ncols:', nrows, ncols
         go to 99
      end if

      do line = 1, nrows
         read (lun, '(a)', iostat=ios) table_io_buffer
         if (ios /= 0) then
            write (*, '(3a, i6)') '*** Trouble reading table ', trim (table%filename), ' at numeric line', line
            go to 99
         end if
         call token2l (trim (table_io_buffer), seps, ncols, table%tokens(:,line))  ! No need to upcase these values
      end do

      close (lun)

 99   return

      end subroutine table_io_read_alpha

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine table_io_read_real (lun, table, ios)  ! Read a table as character tokens then decode numeric tokens as reals.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer,           intent (in)    :: lun    ! The table is ready to read from the top, and closed upon return
      type (table_type), intent (inout) :: table  ! filename, nheader, nrows, and ncols have been set on input;
                                                  ! header(:) and tokens(:,:) and values(:,:) are set upon return
      integer,           intent (out)   :: ios    ! 0 means no problem

!     Local variables:

      integer :: icol, l

!     Execution:

      call table_io_read_alpha (lun, table, ios)  ! Returns any header line(s) and all numeric entries as character tokens

      if (ios /= 0) then
         write (*, '(a)') '*** Bad return from table_io_read_alpha'
         go to 99
      end if

      nrows = table%nrows
      ncols = table%ncols

      allocate (table%values(ncols,nrows), stat=ios)
      if (ios /= 0) then
         write (*, '(a)') '*** Trouble allocating table%values:', nrows, ncols
         go to 99
      end if

      do line = 1, nrows
         do icol = 1, ncols
            l = len_trim (table%tokens(icol,line))
            call decode_number (table%tokens(icol,line)(1:l), inumber, rnumber, ios)
            if (ios >= 1) then  ! 1 = real; 2 = integer or real
               table%values(icol,line) = rnumber               ! Exploit ios = 2 (integer) some day for distinguishing integer cols?
            else
               write (*, '(5a, i5)') '*** Bad number encountered: ', table%tokens(icol,line)(1:l), ' in table ', &
                  trim (table%filename), '.  Line number:', line
               ios = 1
               go to 99
            end if
         end do
      end do

      ios = 0

 99   return

      end subroutine table_io_read_real

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine table_io_write_alpha (lun, table, sep, ios)

!     Write a table from its form as an array of ASCII tokens.  The columns are separated with the indicated character.
!     If the separator is a blank, then the columns are left-justified and each is of a uniform width (possibly meaning more than
!     one blank between some columns).  Otherwise, strictly one separator is used.
!
!     Run-time-variable formatting is also an option, invoked by passing table_io_format_sep as the "sep" argument and assigning
!     a sensible format string, such as:  table_io_format = '(20a10)'.  A separate write is used for each line, so that the 20 in
!     the example may safely be more than enough for one line.
!
!     No deallocations are performed after the table is written.  See the lun argument definition.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer,           intent (in)    :: lun    ! The output file is opened here, written to, and closed upon return
                                                  ! unless lun = standard output (screen), when the open and close are suppressed
      type (table_type), intent (inout) :: table  ! filename, nheader, nrows, and ncols have been set on input, along with
                                                  ! header(:) and tokens(:,:)
      character (1),     intent (in)    :: sep    ! Desired column separator, probably blank, tab, or comma
      integer,           intent (out)   :: ios    ! 0 means no problem

!     Local constants:

      integer, parameter :: std_out = 6           ! Suppress the open and close if we're writing to the screen

!     Local variables:

      integer :: i, i1, i2, j, l, lform, lname, width
      logical :: screen_output

!     Execution:

      lname = len_trim (table%filename)
      screen_output = lun == std_out
      if (.not. screen_output) then
         open (lun, file=table%filename(1:lname), status='unknown', iostat=ios)
         if (ios /= 0) then
            write (*, '(2a)') '*** Unable to open table file ', table%filename(1:lname)
            go to 99
         end if
      end if

      if (table%nheader > 0) then
         do line = 1, table%nheader
            write (lun, '(a)', iostat=ios) trim (table%header(line))
            if (ios /= 0) then
               write (*, '(3a, i4)') '*** Trouble writing table ', trim (table%filename), ' at header line', line
               go to 99
            end if
         end do
      end if

      ncols = table%ncols  ! No need to initialize it for every line, because table_io_scan should already have checked all lines
      nrows = table%nrows

      select case (sep)

         case (blank)  ! Columns will be aligned/left-justified

            do i = 1, ncols  ! Make columns no wider than necessary
               width = 0
               do j = 1, nrows
                  width = max (width, len_trim (table%tokens(i,j)))
               end do
               table%column_width(i) = width
            end do

            do j = 1, nrows           ! Fill a buffer for each line, but trap going off the end of the buffer
               table_io_buffer = blank
               i1 = 1
               i2 = 0
               do i = 1, ncols
                  l = table%column_width(i)
                  i2 = i1 + l
                  if (i2 > max_table_buffer) then
                     if (ios == 0) then  ! Avoid more than one error message
                         ios =  1
                         write (*, '(3a, 2i6)') '*** Buffer too short for table ', trim (table%filename), '; row, column:', j, i
                     end if
                  else
                     table_io_buffer(i1:i2) = table%tokens(i,j)(1:l+1)  ! At least 1 blank between columns ...
                  end if
                  i1 = i2 + 1
               end do
               if (ios /= 0) go to 90

               write (lun, '(a)') table_io_buffer(1:i2-1)               ! ... but not on the end
            end do

         case (table_io_format_sep)  ! Run-time control of output formatting

            lform = len_trim (table_io_format)

            do j = 1, nrows
               write (lun, table_io_format(1:lform), iostat=ios) table%tokens(:,j)
               if (ios /= 0) then
                  write (*, '(5a, i5)') &
                     '*** Trouble writing ', table%filename(1:lname), ' with format ', table_io_format(1:lform), ' at line', j
                  go to 99
               end if
            end do

         case default  ! Non-blank separator: strictly 1 between columns - a little simpler

            do j = 1, nrows           ! Fill a buffer for each line, but trap going off the end of the buffer
               i1 = 1
               i2 = 0
               do i = 1, ncols
                  l = len_trim (table%tokens(i,j))
                  i2 = i1 + l
                  if (i2 > max_table_buffer) then
                     if (ios == 0) then  ! Avoid more than one error message
                         ios =  1
                         write (*, '(3a, 2i6)') '*** Buffer too short for table ', trim (table%filename), '; row, column:', j, i
                     end if
                  else
                     table_io_buffer(i1:i2-1) = table%tokens(i,j)(1:l)
                     table_io_buffer(i2:i2)   = sep
                  end if
                  i1 = i2 + 1
               end do
               if (ios /= 0) go to 90

               write (lun, '(a)') table_io_buffer(1:i2-1)               ! Suppress the last separator
            end do

      end select

      if (.not. screen_output) close (lun)
      go to 99

 90   write (*, '(a, 2i5)') '    Current and needed buffer length:', max_table_buffer, i2 - 1

 99   return

      end subroutine table_io_write_alpha

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine table_io_write_real (lun, table, sep, ios)

!     Write a table as reals from %values(:,:).  The hard-coded formatting options may be an issue, but run-time-variable
!     formatting is also an option, invoked by passing table_io_format_sep as the "sep" argument and assigning a sensible
!     format string, such as:   table_io_format = '(f5.2, es9.2, f7.3, 20es16.8)'
!
!     A separate write is used for each run-time formatted line, so that the 20 in the example may safely be more than enough
!     for one line.
!
!     Later option, if sep = table_io_format_sep:  use of 'i' or 'I' editing in table_io_format causes the indicated column(s)
!     to be written as integers, even though they're stored as reals.  E.g.:  table_io_format = '(i1, f7.1, 4es16.8, i2)' will
!     write the nearest integer for two of the columns.
!
!     No deallocations are performed after the table is written.  See the lun argument definition.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer,           intent (in)    :: lun    ! The table file is opened here as "unknown", written to, then closed upon return
                                                  ! unless lun = standard output (screen), when the open and close are suppressed
      type (table_type), intent (inout) :: table  ! In this version, %tokens(:,:) are not accessed; most other table components are
      character (1),     intent (in)    :: sep    ! Desired column separator, probably blank, tab, or comma, but see also the use
                                                  ! of sep = table_io_format_sep
      integer,           intent (out)   :: ios    ! 0 means no problem

!     Local constants:

      integer, parameter :: std_out = 6           ! Omit the open and close if lun = standard output

!     Local variables:

      integer        :: i, j, lname
      integer, save  :: lform  ! Internal procedure uses it
      logical        :: screen_output
      character (11) :: format1
      character (24) :: format2

!     Execution:

      lname = len_trim (table%filename)
      screen_output = lun == std_out
      if (.not. screen_output) then
         open (lun, file=table%filename(1:lname), status='unknown', iostat=ios)
         if (ios /= 0) then
            write (*, '(2a)') '*** Unable to open table file ', table%filename(1:lname)
            go to 99
         end if
      end if

      do line = 1, table%nheader
         write (lun, '(a)') trim (table%header(line))
      end do

      nrows = table%nrows
      ncols = table%ncols

      select case (sep)

         case (blank)

            format1 = '(   es16.8)'
            write (format1(2:4), '(i3)') ncols

            write (lun, format1) table%values(:,:)  ! For now, a single write, but future enhancements may require fancier footwork.

         case (table_io_format_sep)  ! Run-time control of output formatting, possibly with some columns output as integers

            lform = len_trim (table_io_format)

            if (index (table_io_format(1:lform), 'i') == 0 .and. index (table_io_format(1:lform), 'I') == 0) then  ! All reals

               do j = 1, nrows
                  write (lun, table_io_format(1:lform), iostat=ios) table%values(:,j)
                  if (ios /= 0) then
                     write (*, '(5a, i5)') &
                        '*** Trouble writing to ', table%filename(1:lname), ' with format ', table_io_format(1:lform), &
                        ' at numeric line', j
                     go to 99
                  end if
               end do

            else  ! Some tricky gyrations are modularized as an internal procedure to handle output of [some] reals as integers

               call table_io_write_i_or_r ()
               if (ios /= 0) go to 99

            end if

         case default  ! Non-blank separator

            format2 = '(es16.8,   (a1, es16.8))'
            write (format2(9:11), '(i3)') ncols - 1

            do j = 1, nrows
               write (lun, format2) table%values(1,j), (sep, table%values(i,j), i = 2, ncols)
            end do

      end select

      if (.not.screen_output) close (lun)
      ios = 0

 99   return

!     Internal procedure for subroutine table_io_write_real:

      contains

!        -----------------------------------
         subroutine table_io_write_i_or_r ()  ! Write some columns as integers even though they're stored as reals.
!        -----------------------------------

!        Strategy:  Adjust the run-time format provided in table_io_format(1:lform) by replacing any 'i' or 'I' editing token
!                   such as 2i4 with 2f5.0; for each table line, write with the adjusted format to a buffer, then convert any
!                   '. ' substrings in the buffer to ' ' before writing the adjusted buffer to the file.  If the last non-blank
!                   character in the buffer is '.', that is also suppressed, but only if it is the result of an integer format.
!                   Avoiding interior conversion of (say) f5.0 to i4 is awkward and not attempted here.
!
!                   No doubt, some input run-time format could thwart such a strategy, but this is unlikely.  Replacing what is
!                   written with nint (value) and using the given format as is would be another approach, but handling leading
!                   multipliers such as 4i3 appears more awkward.

!        Local variables:

         integer :: i1, i2, lbuf, lnew, nshift
         logical :: keep_last_period
         character (4) :: iformat
         character (max_table_format) :: rformat  ! Avoid adjusting the input format in-place

!        Execution:

         lnew    = lform
         nshift  = 0
         iformat = '(i1)'
         rformat = table_io_format
         keep_last_period = rformat(lform-2:lform) == '.0)'

         call upcase (rformat(1:lform))

         do i = 2, lform - 1
            if (table_io_format(i:i) == 'i' .or. &
                table_io_format(i:i) == 'I') then  ! Is it I1, I10, or what?
                i1 = i + nshift                    ! Equivalent index in the copy that is being right-shifted
                rformat(i1:i1) =  'F'              ! We'll replace I1 with F2.0, I9 with F10.0, I10 with F11.0, and so on
                i2 = i1 + 1
                do  ! Until ',' or ')' is found
                    if (rformat(i2:i2) == ',' .or. rformat(i2:i2) == ')') exit
                    i2 = i2 + 1
                end do
                if (rformat(i1:i1+1) == 'F9') then  ! Special awkward case: right shift by 3, not 2
                    nshift = nshift + 3
                    lnew   = lnew   + 3
                    rformat(i2+3:lnew) = rformat(i2:lnew-3)
                    rformat(i1+1:i2+2) = '10.0'
                else
                    call decode_number (rformat(i1+1:i2-1), inumber, rnumber, ios)
                    if (ios /= 2) then
                       write (*, '(2a)') '*** table_io_write_i_or_r: bad integer spec in run-time format ', table_io_format(1:lform)
                       ios = 1
                       go to 99
                    end if
                    write (iformat(3:3), '(i1)') i2-i1-1  ! # digits in the i format
                    nshift = nshift + 2
                    lnew   = lnew   + 2
                    rformat(i2+2:lnew) = rformat(i2:lnew-2)
                    write (rformat(i1+1:i2-1), iformat) inumber + 1
                    rformat(i2:i2+1) = '.0'
                end if
            end if
         end do

!        Now write each table line to a buffer with the modified all-real format, suppress any unwanted decimal points,
!        then write the left-shifted buffer to the output file:

         table_io_buffer = blank
         do j = 1, nrows
            write (table_io_buffer, rformat(1:lnew), iostat=ios) table%values(:,j)
            if (ios /= 0) then
               write (*, '(5a, i5)') &
                  '*** Trouble writing to ', table%filename(1:lname), ' with format ', rformat(1:lnew), ' at numeric line', j
               go to 99
            end if
            nshift = 0
            lbuf = len_trim (table_io_buffer)

            do i = 2, len_trim (table_io_buffer) - 2
               i1 = i - nshift
               if (table_io_buffer(i1:i1+1) == '. ') then
                   nshift = nshift + 1
                   lbuf   = lbuf   - 1
                   table_io_buffer(i1:lbuf) = table_io_buffer(i1+1:lbuf+1)  ! Shift left 1 character
               end if
            end do

            if (.not. keep_last_period) then
               if (table_io_buffer(lbuf:lbuf) == '.') lbuf = lbuf - 1
            end if

            write (lun, '(a)') table_io_buffer(1:lbuf)
         end do

 99      return

         end subroutine table_io_write_i_or_r

      end subroutine table_io_write_real

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine table_io_thin (table, ithin)

!     Thin a table by repacking every nth row, as may be needed for a time history.
!     The rows are packed in place and the row count is updated.  The row/col arrays (data proper) are NOT deallocated and
!     reallocated (although they could be if the need arises).

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      type (table_type), intent (inout) :: table  ! All columns are thinned;
                                                  ! the number of rows changes
      integer,           intent (in)    :: ithin  ! If ithin = 2, rows 1,3,5, ..
                                                  ! are repacked as rows 1,2, ..
!     Local variables:

      integer :: jin, jout

!     Execution:

      jout = 0
      do jin = 1, table%nrows, ithin
         jout = jout + 1
         table%values(:,jout) = table%values(:,jin)
      end do

      table%nrows = jout

      end subroutine table_io_thin

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine table_io_deallocate (table, ios)  ! We may need to arrange to deallocate tokens | values | both, but not yet.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      type (table_type), intent (inout) :: table
      integer,           intent (out)   :: ios     ! 0 means no problem

!     Execution:

      if (table%nheader > 0) then
         deallocate (table%header, stat=ios)
         if (ios /= 0) go to 90
      end if

      deallocate (table%tokens, stat=ios)
      if (ios /= 0) go to 90

      deallocate (table%values, stat=ios)
      if (ios /= 0) go to 90

      go to 99

   90 write (*, '(2a)') '*** Trouble deallocating table ', trim (table%filename)

   99 return

      end subroutine table_io_deallocate

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine table_io_copy (table1, table2, new, real, ios)  ! Copy table1 to new table2, which is allocated here appropriately.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      type (table_type), intent (in)  :: table1  ! Input table (treated as alphnumeric or real, according to the real argument)
      type (table_type), intent (out) :: table2  ! Output copy of the input table
      logical,           intent (in)  :: new     ! T means allocate table2 arrays; F means they're allocated already
      logical,           intent (in)  :: real    ! The application needs to know if table%values(:,:) is present
      integer,           intent (out) :: ios     ! 0 means no problem

!     Local variables:

      integer :: ncolumns, nheader, nrows

!     Execution:

      nheader = table1%nheader;  table2%nheader = nheader
      nrows   = table1%nrows;    table2%nrows   = nrows
      ncols   = table1%ncols;    table2%ncols   = ncols

      if (nheader > 0) then
         if (new) then
            allocate (table2%header(nheader), stat=ios)
            if (ios /= 0) then
               write (*, '(a)') '*** Trouble allocating table2%header:', nheader
               go to 99
            end if 
         end if
         table2%header(:) = table1%header(:)
      end if

      if (new) then
         allocate (table2%tokens(ncols,nrows), stat=ios)
         if (ios /= 0) then
            write (*, '(a)') '*** Trouble allocating table2%tokens:', nrows, ncols
            go to 99
         end if
      end if
      table2%tokens(:,:) = table1%tokens(:,:)

      if (real) then
!!!      if (allocated (table2%values)) then  ! Not allowed for pointers; only for allocatable arrays
         if (new) then
            allocate  (table2%values(ncols,nrows), stat=ios)
            if (ios /= 0) then
               write (*, '(a)') '*** Trouble allocating table2%values:', nrows, ncols
               go to 99
            end if
         end if
         table2%values(:,:) = table1%values(:,:)
      end if

 99   return

      end subroutine table_io_copy

   end module table_io
