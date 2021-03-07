!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program filter_rows
!
!  Description:
!
!     This is a slightly specialized utility that may prove handy in a variety
!     of situations.  It looks for lines common to two input text files and
!     writes a third file with any duplicates removed.  The option to include
!     or exclude any non-duplicates from the second file is provided.
!
!     The application prompting this functionality is to two lists of files
!     produced by something like "find . -name PP" which lists recursively all
!     the files at or below the working directory that are named PP.  E.g., in
!     a directory containing directories for multiple sets of flight conditions,
!     something like this may show up:
!
!        ./v7.0/h64/a162/PP
!        ./v7.0/h64/a162/Baldwin-Lomax/PP
!        ./v7.0/h64/a162/RCG/PP
!        ./v7.0/h72/a162/PP
!        ./v10.95/h73.6/a162.5/PP
!        ::::::::
!
!     An earlier list may contain a subset of this list and any duplicates need
!     to be deleted for some reason such as for archiving the directories.
!
!  Restrictions:
!
!     Using the table_io module simplifies the coding, but it means that the
!     lists should contain only one column.  Any embedded blanks or other
!     delimiters will be read as more than one column, complicating string
!     comparisons.  Sorry ...
!
!   History:
!
!     07/06/17  D.A.Saunders  Initial implementation to assist archiving of
!                             flow calculations.
!     07/07/17    "      "    Got the "include nonduplicates" option going.
!
!  Author:  David Saunders, AMA, Inc./NASA Ames Research Center, Mtn. View, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use table_type_module  ! Both are in table_io.f90
   use table_io

   implicit none

!  Constants:

   integer, parameter :: &
      lun1   = 1,        &  ! Input list 1
      lun2   = 2,        &  ! Input list 2
      lun3   = 3,        &  ! Output list 3
      lunkbd = 5,        &  ! Keyboard inputs
      luncrt = 6,        &  ! Screen prompts and diagnostics
      maxbuf = 512          ! Should match table_io module

!  Variables:

   integer           :: ios, l1, l2
   integer           :: irow1,  irow2,  irow3,  &
                        ncols1, ncols2, ncols3, &
                        nrows1, nrows2, nrows3
   logical           :: cr, eof              ! Handle carriage return or ctrl D
   logical           :: duplicate, include_nonduplicates
   type (table_type) :: table1, table2, table3

!  Execution:


   call reads (luncrt, &
               'First list which may have lines common to the second list: ', &
               lunkbd, table1%filename, cr, eof)
   if (cr .or. eof) go to 99

!  Read the list as character tokens:

   call table_io_read_alpha (lun1, table1, ios)
   if (ios /= 0) go to 99

   nrows1 = table1%nrows
   ncols1 = table1%ncols

   write (luncrt, '(a, 2i6)') ' # rows and columns found:', nrows1, ncols1

   if (ncols1 /= 1) then
      write (luncrt, '(a)') 'Only single-column lists are handled -- sorry.'
      go to 99  ! Single stop philosophy
   end if

   call reads (luncrt, &
               'Second list containing possible duplicates: ', &
               lunkbd, table2%filename, cr, eof)
   if (cr .or. eof) go to 99

!  Read the list as character tokens:

   call table_io_read_alpha (lun2, table2, ios)
   if (ios /= 0) go to 99

   nrows2 = table2%nrows
   ncols2 = table2%ncols

   write (luncrt, '(a, 2i6)') ' # rows and columns found:', nrows2, ncols2

   if (ncols2 /= 1) then
      write (luncrt, '(a)') 'Only single-column lists are handled.  Sorry.'
      go to 99
   end if

   include_nonduplicates = .false.
   call ready (luncrt, 'Retain nonduplicates from file 2? [y|n; <cr>=n] ', &
               lunkbd, include_nonduplicates, cr, eof)
   if (eof) go to 99

   if (include_nonduplicates) then
      nrows3 = nrows1 + nrows2   ! Max.
   else
      nrows3 = nrows1            ! Max.
   end if

   table3%nrows           = nrows3
   table3%ncols           = 1
   table3%nheader         = 0
   table3%column_type(1)  = 3    ! Character data
   table3%column_width(1) = max_table_token_length
   table3%filename        = 'filtered_list.txt'

   call reads (luncrt, 'Output list name? [<cr>=filtered_list.txt] ', &
               lunkbd, table3%filename, cr, eof)
   if (eof) go to 99

   allocate (table3%tokens(1,nrows3))

!  Process the two lists, one line of list 1 at at time:

   nrows3 = 0;  l2 = 0  ! In case the first l1 exceeds the limit
   do irow1 = 1, nrows1
      duplicate = .false.
      l1 = len_trim (table1%tokens(1,irow1))
      if (l1 == max_table_token_length) go to 90  ! May have been truncated
      do irow2 = 1, nrows2
         l2 = len_trim (table2%tokens(1,irow2))
         if (l2 == max_table_token_length) go to 90  ! May have been truncated
         if (l2 == l1) then
            if (table2%tokens(1,irow2)(1:l1) == &
                table1%tokens(1,irow1)(1:l1)) then
               duplicate = .true.
               exit
            end if
         end if
      end do
      if (.not. duplicate) then
         nrows3 = nrows3 + 1
         table3%tokens(1,nrows3) = table1%tokens(1,irow1)
      end if
   end do

!  If nonduplicates from list 2 are to be included, tack them on the end:

   if (include_nonduplicates) then
      do irow2 = 1, nrows2
         duplicate = .false.
         l2 = len_trim (table2%tokens(1,irow2))
         do irow1 = 1, nrows1
            l1 = len_trim (table1%tokens(1,irow1))
            if (l1 == l2) then
               if (table2%tokens(1,irow2)(1:l1) == &
                   table1%tokens(1,irow1)(1:l1)) then
                  duplicate = .true.
                  exit
               end if
            end if
         end do
         if (.not. duplicate) then
            nrows3 = nrows3 + 1
            table3%tokens(1,nrows3) = table2%tokens(1,irow2)
         end if
      end do
   end if

   table3%nrows = nrows3
   write (luncrt, '(a, i6)') ' # output lines:', nrows3

   call table_io_write_alpha (lun3, table3, ' ', ios)
   go to 99

90 write (luncrt, '(2a, 3i5)') &
    ' List entry may exceed table_io.f90 limit.', &
    ' l1, l2, max_table_token_length:', &
      l1, l2, max_table_token_length

99 continue  ! Avoid system-dependent stop behavior

   end program filter_rows
