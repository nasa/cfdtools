
   program test_table_io

!  09/17/14  David Saunders  Initial testing of a package that may evolve.

   use table_type_module  ! Both are in table_io.f90
   use table_io

   implicit none

!  Constants:

   integer, parameter :: lunin  = 1  ! Input table
   integer, parameter :: lunout = 2  ! Output table

!  Variables:

   integer :: ios
   character (1) :: sep
   type (table_type) :: tablein  ! The same table is written for test purposes

!  Execution:

   write (*, '(a)', advance='no') 'Input table file: '
   read  (*, '(a)') tablein%filename

   call table_io_read_real (lunin, tablein, ios)

   write (*, '(a, i3)') 'After reading, ios =', ios
   write (*, '(3x, a, i3)') &
      'nheader:', tablein%nheader,  &
      '# rows: ', tablein%nrows,    &
      '# cols: ', tablein%ncols
   if (ios /= 0) go to 99

   write (*, '(a)', advance='no') 'Output table file (alpha): '
   read  (*, '(a)') tablein%filename
   write (*, '(a)', advance='no') 'Separator to use:  '
   read  (*, '(a)') sep

   call table_io_write_alpha (lunout, tablein, sep, ios)

   write (*, '(a, i3)') 'After writing, ios =', ios
   if (ios /= 0) go to 99

   write (*, '(a)', advance='no') 'Output table file (real):  '
   read  (*, '(a)') tablein%filename

   call table_io_write_real (lunout, tablein, sep, ios)

   write (*, '(a, i3)') 'After writing, ios =', ios
   if (ios /= 0) go to 99

   write (*, '(a)', advance='no') 'Output table file (run-time-formatted, alpha):  '
   read  (*, '(a)') tablein%filename
   write (*, '(a)', advance='no') 'Run-time alpha format: '
   read  (*, '(a)') table_io_format

   call table_io_write_alpha (lunout, tablein, table_io_format_sep, ios)

   write (*, '(a, i3)') 'After writing, ios =', ios
   if (ios /= 0) go to 99

   write (*, '(a)', advance='no') 'Output table file (run-time-formatted, real):  '
   read  (*, '(a)') tablein%filename
   write (*, '(a)', advance='no') 'Run-time real format: '
   read  (*, '(a)') table_io_format

   call table_io_write_real (lunout, tablein, table_io_format_sep, ios)

   write (*, '(a, i3)') 'After writing, ios =', ios

99 continue

   end program test_table_io
