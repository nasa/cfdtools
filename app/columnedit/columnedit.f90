!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program columnedit
!
!  Description:
!
!     COLUMNEDIT provides for replacing one column of a text file with any
!     column from another text file, or for inserting a new column, deleting
!     a column, or extracting a column.  Here, "column" is used in the sense
!     of tokens delimited by white space (blanks or tabs).
!
!     This version has the handy option to extract more than one column at a
!     time.  This is implemented in a way that minimizes impact to the one-
!     column-per-run original design.  The prompting approach here is more
!     convenient than the control files required by related utilities such as
!     EXTRACT_COLUMNS and MERGE_TABLES.
!
!  Limitations:
!
!     > Original column spacing may be lost - recover it with an editor.
! 
!  History:
!
!     12/11/96  D.Saunders  Initial implementation to facilitate entering
!                           nonzero design variable values in the control
!                           file of a wing/body design code.
!     05/21/99      "       Minor Fortran 90 revisions.
!     06/13/00      "       How come no "extract column" option?  (Added it.)
!     08/24/00      "       Gnuplot files can be more than 132 characters long.
!     06/30/04      "       DO LINENO = 1, 99999 was a poor way of reading to
!                           end-of-file.  Use no upper limit.
!     12/28/10      "       Free-format version; considered allowing more than
!                           one operation to a file in one run, but it's not
!                           worth the trouble (which would include storing the
!                           entire file as a 2-D array of tokens).
!     11/16/16      "       Max. column width of 32 didn't handle long file
!                           names in columns.  Make it 80.
!                           Also, use TOKEN2L in place of TOKEN2 to avoid
!                           any uppercasing.
!     02/07/17      "       Made use of RDLIST to enable extraction of any
!                           number of columns, not just one.
!     02/08/17      "       List-directed reading of tokens for multicolumn
!                           extraction doesn't handle commas in header tokens.
!     07/18/18      "       Handle bad file names properly.
!
!  Author:  David Saunders, Sterling Software/ELORET Corp./ERC, Inc.
!                           NASA Ames Research Center, Moffett Field, CA
!     Now with AMA, Inc. at NASA Ames.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Constants:

   integer, parameter :: &
      luncrt    = 6,     &
      lunkbd    = 5,     &
      lun1      = 1,     &  ! File to be modified
      lun2      = 2,     &  ! File providing an inserted/replaced column
      lunout    = 3,     &  ! Modified file
      maxline   = 512,   &  ! Max. length of a line + 1
      maxcolumn = 200,   &  ! Max. # columns handled
      maxtoken  = 200       ! Max. width of any column (token)

   character, parameter :: &
      blank*1 = ' '

!  Variables:

   integer :: &
      count, first, icol, icol1, icol2, ios, last, lineno, mark, ncol1, ncol2
   logical :: &
      cr, delete, eof, extract, extract_many, insert, replace
   character :: &
      answer*1, filename*48, line1*(maxline), line2*(maxline), &
      tokens1(maxcolumn)*(maxtoken), tokens2(maxcolumn)*(maxtoken), white*2

!  Execution:

   white(1:1) = blank  ! White space = blank or tab
   white(2:2) = char (9)

   ios = 1
   write (luncrt, '(/, a)', advance='no') ' File to be processed: '
   do while (ios /= 0)
      read (lunkbd, '(a)', iostat=ios) filename
      open (lun1, file=filename, status='old', iostat=ios)
      if (ios < 0) go to 99  ! Quit
      if (ios > 0) write (luncrt, '(a)') ' File not opened; try again.'
   end do

10 answer = 'E'
   call readc (luncrt, &
      'Replace|insert|delete|extract a column|extract more than one? ' // &
      '(r|i|d|e|m; [e]): ', lunkbd, answer, cr, eof)
   if (eof) go to 99

   replace = answer == 'R'
   insert  = answer == 'I'
   delete  = answer == 'D'
   extract = answer == 'E'
   extract_many = answer == 'M'

   if (replace) then
      write (luncrt, '(a)', advance='no') ' Number of column to replace: '
   else if (insert) then
      write (luncrt, '(a)', advance='no') ' Number of column to insert: '
   else if (delete) then
      write (luncrt, '(a)', advance='no') ' Number of column to delete: '
   else if (extract) then
      write (luncrt, '(a)', advance='no') ' Number of column to extract: '
   else if (extract_many) then
      call more_than_one_column ()  ! It's cleaner to treat this separately
      go to 99
   else
      write (luncrt, '(a)') ' Invalid choice.  Try again.'
      go to 10
   end if

   ios = 1
   do while (ios /= 0)
      read (lunkbd, *, iostat=ios) icol1
      if (ios < 0) go to 99  ! Quit
      if (ios > 0) write (luncrt, '(a)') ' Bad entry; try again.'
   end do

   if (.not. (delete .or. extract)) then
      write (luncrt, '(a)', advance='no') ' File with desired new column: '
      ios = 1
      do while (ios /= 0)
         read (lunkbd, '(a)', iostat=ios) filename
         open (lun2, file=filename, status='old', iostat=ios)
         if (ios < 0) go to 99
         if (ios > 0) write (luncrt, '(a)') ' File not opened; try again.'
      end do

      write (luncrt, '(a)', advance='no') ' Number of column to transcribe: '
      ios = 1
      do while (ios /= 0)
         read (lunkbd, *, iostat=ios) icol2
         if (ios < 0) go to 99
         if (ios > 0) write (luncrt, '(a)') ' Bad entry; try again.'
      end do
   end if

   ios = 1
   write (luncrt, '(a)', advance='no') ' Output file name: '
   do while (ios /= 0)
      read (lunkbd, '(a)', iostat=ios) filename
      open (lunout, file=filename, status='new', iostat=ios)
      if (ios < 0) go to 99  ! Quit
      if (ios > 0) then
         write (luncrt, '(a)', advance='no') &
            ' File already exists. Overwrite it? [y|n]: '
         read (lunkbd, '(a)', iostat=ios) answer
         if (ios < 0) go to 99  ! Quit
         if (answer == 'y') then
            open (lunout, file=filename, status='unknown')
            exit
         end if
         ios = 1
      end if
   end do 

!  Process file(s) line-by-line until we hit an EOF:

   lineno = 0

   do !!!! lineno = 1, 99999  ! No: may not be high enough!

!     Read and tokenize a line from the first file:

      read (lun1, '(a)', iostat=ios) line1
      if (ios < 0) exit ! EOF

      lineno = lineno + 1
      if (ios /= 0) go to 70

      if (extract) then
         ncol1 = icol1
      else
         ncol1 = maxcolumn
      end if

      call token2l (line1, white, ncol1, tokens1)

      if (replace .or. delete .or. extract) then
         if (ncol1 < icol1) go to 90
      else if (insert) then
         if (ncol1 + 1 < icol1) go to 90
      end if
    
      if (delete) then  ! Shift columns left

         do icol = icol1 + 1, ncol1
            tokens1(icol - 1) = tokens1(icol)
         end do
         ncol1 = ncol1 - 1

      else if (extract) then  ! Copy the token to column 1 for convenience

         tokens1(1) = tokens1(icol1)
         ncol1 = 1

      else if (replace .or. insert) then

!        Read and tokenize a line from the second file:

         read (lun2, '(a)', iostat=ios) line2
         if (ios /= 0) go to 80

         ncol2 = maxcolumn

         call token2l (line2, white, ncol2, tokens2)

         if (ncol2 < icol2) go to 91
    
         if (replace) then

            tokens1(icol1) = tokens2(icol2)

         else if (insert) then  ! Shift right ...

            do icol = ncol1, icol1, -1
               tokens1(icol + 1) = tokens1(icol)
            end do

            tokens1(icol1) = tokens2(icol2)  ! ... and insert
            ncol1 = ncol1 + 1

         end if

      end if

!     Compose the output line in line2:

      count = 1
      do icol = 1, ncol1
         first = 1
         last  = maxtoken

         call scan2 (tokens1(icol), blank, first, last, mark)

         line2(count : count + mark - 1) = tokens1(icol)(1 : mark)
         count = count + mark
         line2(count : count) = blank
         count = count + 1
      end do
 
!     Output the updated line:

      write (lunout, '(a)') line2(1 : count - 2)

   end do ! Next input line

   write (luncrt, '(/, 1x, i6, a, a, /)') &
      lineno, ' lines processed into file ', trim (filename)

   go to 99


!  Error handling:

70 write (luncrt, '(/, a, i7)') &
      ' Error reading from file 1 at line #', lineno
   write (luncrt, '(1x, a)') 'Offending line:', line1(1 : 131)
   go to 99

80 write (luncrt, '(/, a, i7)') &
      ' Error reading from file 2 at line #', lineno
   write (luncrt, '(1x, a)') 'Offending line:', line2(1 : 131)
   go to 99

90 write (luncrt, '(/, a, i7)') &
      ' First file has too few columns in line #', lineno
   write (luncrt, '(1x, a)') 'Offending line:', line1(1 : 131)
   go to 99

91 write (luncrt, '(/, a, i7)') &
      ' Second file has too few columns in line #', lineno
   write (luncrt, '(1x, a)') 'Offending line:', line2(1 : 131)

99 continue

   contains

!     Internal procedure for program columnedit:

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine more_than_one_column ()  ! Clean separation for special case
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local variables:

      integer :: i, l, ncols
      integer :: list(maxcolumn)  ! List of column numbers to extract

!     Execution:

!     Counting the tokens in line 1 (and all other lines) simplifies things:

      read (lun1, '(a)', iostat=ios) line1
      if (ios /= 0) then
         write (luncrt, '(a)') ' Error reading first line.'
         go to 99
      end if

      call token_count (line1, white, ncol1)
      if (ncol1 <= 0) then
         write (luncrt, '(a)') ' No columns found on line 1.'
         go to 99
      else
         if (ncol1 > maxcolumn) write (luncrt, '(a, i5)') &
            ' # columns found on line 1:', ncol1, &
            ' # columns handled:        ', maxcolumn
         ncol1 = min (ncol1, maxcolumn)
      end if

      write (luncrt, '(a, i5)') ' # columns expected on all lines:', ncol1
      ncol2 = ncol1
      call rdlist (luncrt, 'Column numbers to extract, in desired order: ', &
                   lunkbd, ncol2, list)
      if (ncol2 <= 0) go to 99

      rewind (lun1)

      write (luncrt, '(a)', advance='no') ' Output file name: '
      read (lunkbd, '(a)', end=99) filename
      open (lunout, file=filename, status='new', iostat=ios)
      if (ios /= 0) then
         write (luncrt, '(/, 2a)') &
            ' Output file already exists: ', trim (filename)
         go to 99
      end if

      lineno = 0
      do  ! Until EOF
!!!      read (lun1, *, iostat=ios) tokens1(1:ncol1)  ! Bad if header commas
         read (lun1, '(a)', iostat=ios) line1
         if (ios < 0) exit

         lineno = lineno + 1
         if (ios /= 0) exit

!        Tokenize the line:

         ncols = ncol1
         call token2l (line1, white, ncols, tokens1)
         if (ncols /= ncol1) then
            write (luncrt, '(a, 2i5)') &
               ' Bad # columns; found/expected:', ncols, ncol1
            ios = 1
            exit
         end if

!        Compose the output line and write it:

         count = 1
         do icol = 1, ncol2
            i = list(icol)
            l = len_trim (tokens1(i))
            line2(count:count+l-1) = tokens1(i)(1:l)
            count = count + l
            line2(count:count) = blank
            count = count + 1
         end do

         write (lunout, '(a)') line2(1:count-2)

      end do

      if (ios > 0) then
         write (luncrt, '(a, i7)') ' Read error; line number:', lineno
         go to 99
      end if

      write (luncrt, '(/, 1x, i6, a, a, /)') &
         lineno, ' lines processed into file ', trim (filename)

   99 return

      end subroutine more_than_one_column

   end program columnedit
