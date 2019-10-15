!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program table_arithmetic
!
!  Description:
!     This aid to Excel novices is intended to ease the burden of compiling
!     tables of results from multiple related analyses.  It was prompted by
!     the need to adjust units (e.g., convert W/m^2 to W/cm^2) (scale and/or
!     shift operations) and to form ratios of radiative to convective heat
!     flux at a given body point on an atmospheric entry vehicle (implying
!     input of a second table and operations on matching columns).
!
!     Some options involve only one table, and multiple operations can be
!     performed on any of the columns in one run.  Output number formats
!     can be specified at run time (e.g., an input E format may be preferred
!     as an output F format, for preparing a slide).  Reals desired to be
!     treated as integers are handled with a format such as i3 for column
!     1 in the example below.
!
!     Any second table, at least initially, is expected to match the first
!     in terms of row counts (cf. matrix multiplication, which requires
!     certain dimensions to match).
!
!     See also the author's much earlier COLUMNEDIT for extracting, inserting,
!     or replacing table columns.  The TABLE_IO module and the RDLIST subroutine
!     are key building blocks, without which contemplating this type of utility
!     would probably be foolhardy.  (RDLIST allows entry of column lists, say,
!     via any reasonable shorthand, such as 2:8 for the sample table below.)
!
!  Table Inputs:
!     Any initial lines that are not purely numeric are treated as header
!     information by the table_io package, and ignored.  Any header records
!     from table 1 are transferred to the output table.  An example of a
!     radiative heat flux table follows:
!
!         t,s  BP 1   BP 2   BP 3   BP 4   BP 5   BP 6   BP 7
!         50 0.2887 0.1785 0.1341 0.1383 0.1005 0.1362 0.1946
!         59 0.7429 0.4754 0.3245 0.2612 0.2134 0.2669 0.4218
!         78 2.9063 1.8416 1.1830 0.9051 0.8218 1.1402 0.4218
!         90 4.3462 2.5807 1.6881 1.2732 1.1256 1.5980 1.9729
!         95 4.6270 2.8072 1.9010 1.4455 1.2590 1.6965 3.5093
!        104 3.6444 2.2715 1.8328 1.5137 1.3239 1.5524 3.8476
!        112 1.3757 0.8778 0.7848 0.6097 0.5070 0.5730 3.4849
!        121 0.3423 0.2136 0.1826 0.1289 0.0978 0.1198 1.4563
!        130 0.0918 0.0500 0.0600 0.0600 0.0492 0.0600 0.3987
!
!      All rows of specified columns (except headers) are operated on.
!
!  Programmer Note:
!      If a menu is added to, be sure to adjust the appropriate prompt for a
!      menu choice, where ^D and maxmenu* are both indicated.
!
!  History:
!     11/17/2018  D.A.Saunders  Initial design.
!     11/19/2018    "      "    Initial coding and testing completed.
!     11/20/2018    "      "    Displaying a table on the screen required a
!                               a change in table_io.f90, which now has a
!                               table_io_copy option that allows starting
!                               over with either table.  Also, table 2 may
!                               now be scaled or shifted.
!     08/15/2018    "      "    The requirement that a second table needs the
!                               same number of columns as the first was more
!                               restrictive than it needed to be.
!
!  Author:  David Saunders, AMA, Inc. at NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use table_type_module  ! Both are in table_io.f90
   use table_io

   implicit none

!  Constants:

   integer, parameter :: &
      lunin1   = 1,      &  ! Input file/table
      lunin2   = 2,      &  ! Optional second input table
      lunout   = 3,      &  ! Output table
      lunkbd   = 5,      &  ! Keyboard inputs
      luncrt   = 6,      &  ! Screen prompts and diagnostics
      lmenu1   = 31,     &  ! Length of longest menu1 entry
      lmenu2   = 33,     &  !   "    "    "     menu2   "
      maxbuf   = 512,    &  ! Should match table_io module
      maxmenu1 = 4,      &  ! Number of menu 1 entries
      maxmenu2 = 6          !   "    "    "  2    "

   real, parameter :: &
      eps = 1.e-32          ! Guard against divides by zero

!  Variables:

   integer :: &
      i, icol1, icol2, ios, j, mcols, menu_choice, ncols1, ncols, ncols2, &
      nrows1, nrows2

   real :: &
      scale, shift

   logical :: &
      cr, eof

   character (lmenu1) :: &
      menu1(-2:maxmenu1 + 1)  ! Allow for menu items -2, -1, 0, and "Done"

   character (lmenu2) :: &
      menu2(-2:maxmenu2 + 1)

   character (max_table_file_name) :: &
      table1_name, table2_name

   integer, allocatable :: &
      icols1(:), icols2(:)

   type (table_type) :: &
      table1, table2, table1_original, table2_original

!  Storage:

   data menu1 / &
      '  -2: Start over with table 1',     &
      '  -1: Undo the last operation',     &
      '   0: Review table 1',              &
      '   1: Scale table 1 column(s)',     &
      '   2: Shift table 1 column(s)',     &
      '   3: Reverse table 1 row order',   &
      '   4: Read a second table',         &
      '   5: Done'/

   data menu2 / &
      '  -2: Start over with table 2',     &
      '  -1: Undo the last operation',     &
      '   0: Review table 2',              &
      '   1: Add table 1 & 2 column(s)',   &
      '   2: Multiply corresp. col(s)',    &
      '   3: Calculate table 1 - table 2', &
      '   4: Calculate table 1 / table 2', &
      '   5: Scale table 2 column(s)',     &
      '   6: Shift table 2 column(s)',     &
      '   7: Done'/

!  Execution:

   write (luncrt, '(/, (a))') &
      ' N.B.: Computed results overwrite table 1 in memory, but', &
      '       the output table name must differ from the input name.'
   call reads (luncrt, 'Name of first (only?) table: ', &
               lunkbd, table1_name, cr, eof)
   if (cr .or. eof) go to 99

   table1%filename = table1_name  ! We avoid overwriting input files
   table2_name     = table1_name  ! Avoid a possible undefined name

!  Read table 1 as real tokens:

   call table_io_read_real (lunin1, table1, ios)
   if (ios /= 0) go to 99

   nrows1 = table1%nrows
   ncols1 = table1%ncols

   call table_io_copy (table1, table1_original, .true., .true., ios)
   if (ios /= 0) then
      write (luncrt, '(a, 2i6)') '*** Trouble cloning table 1:', nrows1, ncols1
      go to 99
   end if

   write (luncrt, '(a, 2i6)') ' # data rows and columns found:', nrows1, ncols1

   allocate (icols1(ncols1))
   ncols = ncols1  ! At most
   call rdlist (luncrt, 'Column number(s) to operate on (e.g., 2:8): ', &
                lunkbd, ncols, icols1)
   if (ncols < 0) then  ! ^D option
      ios = -1
      go to 99
   end if

!  All operations performed on table 1 are performed in-place:

   write (luncrt, '(/, (a))') menu1(:)
   menu_choice = 0
   do while (menu_choice /= maxmenu1 + 1)

      call readi (luncrt, 'Menu 1 choice? [^D or 5 means no more]: ', &
                  lunkbd, menu_choice, cr, eof)
      if (eof) exit

      call operate_on_1 ()
      if (ios < 0) exit

      if (menu_choice == maxmenu1) then  ! We'll read a second table

         call table_io_read_real (lunin2, table2, ios)
         if (ios /= 0) go to 99

         nrows2 = table2%nrows
         ncols2 = table2%ncols

         write (luncrt, '(a, 2i6)') ' # data rows and columns found:', &
            nrows2, ncols2

         if (nrows2 /= nrows1) then
            write (luncrt, '(a, /, (2i6))') &
               ' *** Table row mismatch; row and column counts:', &
               nrows1, ncols1, nrows2, ncols2
            ios = -1
            go to 99
         end if

         call table_io_copy (table2, table2_original, .true., .true., ios)
         if (ios /= 0) then
            write (luncrt, '(a, 2i6)') &
               '*** Trouble cloning table 2:', nrows2, ncols2
            go to 99
         end if

         allocate (icols2(ncols2))
         mcols = ncols2  ! At most
         call rdlist (luncrt, 'Column number(s) to operate on (e.g., 2:8): ', &
                      lunkbd, mcols, icols2)
         if (mcols < 0) then  ! ^D option
            ios = -1
            go to 99
         end if

         if (mcols /= ncols) then
            write (luncrt, '(a, 2i6)') &
               '*** Likely mismatch in numbers of columns to operate on: ', &
               ncols, mcols
            ios = -1
            go to 99
         end if

         write (luncrt, '(/, (a))') menu2(:)
         menu_choice = 0
         do while (menu_choice /= maxmenu2 + 1)
            call readi (luncrt, 'Menu 2 choice? [^D or 7 means no more]: ', &
                     lunkbd, menu_choice, cr, eof)
            if (eof) exit

            call operate_on_1_and_2 ()
            if (ios < 0) exit

         end do  ! Loop over operations on table 1 and 2 exit end if
      end if

   end do  !  Loop over operations on table 1 alone

!  Save results?

   cr = .true.
   do while (cr)  ! Avoid losing results with an unwitting carriage return
      call reads (luncrt, 'Output table name? [^D = do not save results] ', &
                  lunkbd, table1%filename, cr, eof)
      if (eof) go to 99
      if (table1%filename == table1_name .or. &
          table1%filename == table2_name) then
         write (luncrt, '(a)') ' Overwriting an input table is not permitted.'
         cr = .true.
      end if
   end do

   ios = 1
   do while (ios /= 0)
      table_io_format = '(i3, 10f7.4)'
      call reads (luncrt, &
                  'Output table format? E.g.: (i3, 10f7.4) with parens.: ', &
                  lunkbd, table_io_format, cr, eof)
      if (eof) exit

      call table_io_write_real (lunout, table1, table_io_format_sep, ios)
      if (ios /= 0) rewind (lunout)  ! Avoid duplicating header record(s)
   end do

99 continue

!  Internal procedures for program table_arithmetic:

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine operate_on_1 ()  ! Table 1 operations
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      select case (menu_choice)

         case (-2)
            call table_io_copy (table1_original, table1, .false., .true., ios)
            if (ios /= 0) then
               write (luncrt, '(a)') '*** Trouble restoring table1.'
            end if

         case (-1)
            write (luncrt, '(a)') 'Sorry: undoing has not yet been implemented.'

         case (0)  ! Review current table1 values
            call table_io_write_real (luncrt, table1, ' ', ios)

         case (1)
            call readr (luncrt, 'Scale factor: ', lunkbd, scale, cr, eof)
            do i = 1, ncols
               icol1 = icols1(i)
               table1%values(icol1,:) = table1%values(icol1,:)*scale
            end do

         case (2)
            call readr (luncrt, 'Shift: ', lunkbd, shift, cr, eof)
            do i = 1, ncols
               icol1 = icols1(i)
               table1%values(icol1,:) = table1%values(icol1,:) + shift
            end do

         case (3)  ! Reverse the row order of specified columns
            do i = 1, ncols
               icol1 = icols1(i)
               call rverse (nrows1, table1%values(icol1,:), &
                                    table1%values(icol1,:))
            end do

         case (4)  ! Read a second table at the higher level
             cr = .true.
             do while (cr)
                call reads (luncrt, 'Name of second table: ', &
                            lunkbd, table2_name, cr, eof)
             end do
             if (eof) then
                ios = -1
                go to 99
             end if

             table2%filename = table2_name
             ios = 0

         case (5)  ! Done
            ios = -1

         case default
            ios = 0
            write (luncrt, '(a, i6)') 'Invalid choice:', menu_choice

      end select

99    continue

      end subroutine operate_on_1

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine operate_on_1_and_2 ()  ! Table 1 & 2 operations
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ios = 0  ! Allows further menu choices

      select case (menu_choice)

         case (-2)
            call table_io_copy (table2_original, table2, .false., .true., ios)
            if (ios /= 0) then
               write (luncrt, '(a)') '*** Trouble restoring table2.'
            end if

         case (-1)
            write (luncrt, '(a)') 'Sorry: undoing has not yet been implemented.'

         case (0)  ! Review current table2 values
            call table_io_write_real (luncrt, table2, ' ', ios)

         case (1)  ! Add corresponding table entries

            do i = 1, ncols  ! = mcols for table2 (but possibly another list)
               icol1 = icols1(i)
               icol2 = icols2(i)
               table1%values(icol1,:) = table1%values(icol1,:) + &
                                        table2%values(icol2,:)
            end do

         case (2)  ! Multiply corresponding table entries
            do i = 1, ncols
               icol1 = icols1(i)
               icol2 = icols2(i)
               table1%values(icol1,:) = table1%values(icol1,:) * &
                                        table2%values(icol2,:)
            end do

         case (3)  ! Subtract table 2 entries from those of table 1
            do i = 1, ncols
               icol1 = icols1(i)
               icol2 = icols2(i)
               table1%values(icol1,:) = table1%values(icol1,:) - &
                                        table2%values(icol2,:)
            end do

         case (4)  ! Divide table 1 entries by those of table 2
            do i = 1, ncols
               icol1 = icols1(i)
               icol2 = icols2(i)
               do j = 1, nrows1
                  if (abs (table2%values(icol2,j)) < eps) &
                     table2%values(icol2,j) = sign (eps, table2%values(icol2,j))
               end do
               table1%values(icol1,:) = table1%values(icol1,:) / &
                                        table2%values(icol2,:)
            end do

         case (5)
            call readr (luncrt, 'Scale factor: ', lunkbd, scale, cr, eof)
            do i = 1, ncols
               icol2 = icols2(i)
               table2%values(icol2,:) = table2%values(icol2,:)*scale
            end do

         case (6)
            call readr (luncrt, 'Shift: ', lunkbd, shift, cr, eof)
            do i = 1, ncols
               icol2 = icols2(i)
               table2%values(icol2,:) = table2%values(icol2,:) + shift
            end do

         case (7)  ! Done
            ios = maxmenu2 + 1

         case default
            ios = 0
            write (luncrt, '(a, i6)') 'Invalid choice:', menu_choice

      end select

99    continue

      end subroutine operate_on_1_and_2

   end program table_arithmetic
