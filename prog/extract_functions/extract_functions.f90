!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program extract_functions
!
!  Description:
!
!     EXTRACT_FUNCTIONS extracts the indicated list of functions from a PLOT3D
!  function file, probably to avoid redundant flow interpolation for some other
!  application.  The circumstance that prompted it involves working with volume
!  datasets compatible with BLAYER when only pressure, temperature and electron
!  number density are required for a RADIO_BLACKOUT program.
!
!     The input and output function file names are entered via prompts, and the
!  list of function numbers to extract is entered via a control file.  The form
!  of the input file is determined automatically, and the output file takes the
!  the same form (ASCII or unformatted).  There is no need to deal with the
!  associated grid file.
!
!  Control file ('extract_functions.inp'):
!
!     This should contain a list of integers representing the desired functions
!  and their output order.  The list can take any of the forms handled by the
!  RDLIST utility and it should be on a single line.  Examples:
!
!     2 3 14                    (p, T and electron density from a BLAYER input)
!
!     2 3 15:17                 (p, T and velocity components likewise)
!
!     15-17 1:3                 (any sensible list, not necessarily ascending)
!
!  History:
!
!     03/04/2016  D.A.Saunders  Initial adaptation of earlier utilities.
!     03/05/2016    "     "     Introduced f_io.f90 and made use of it for
!                               much more efficient function file I/O.
!
!  Author:  David Saunders, AMA, Inc./NASA Ames Research Ctr., Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_block_structure  ! Derived data type for one grid block
   use xyzq_io_module        ! I/O package for PLOT3D-type files
   use    f_io_module        ! I/O package for PLOT3D-type fn. files (see above)

   implicit none

!  Constants:

   integer, parameter :: &
      ignore     = -99,  &         ! Disables grid/fn. file block count check
      lunctl     = 1,    &         ! Control file of function numbers
      lunf_in    = 2,    &         ! Input  function file
      lunf_out   = 3,    &         ! Output function file
      lunkbd     = 5,    &         ! Keyboard inputs
      luncrt     = 6,    &         ! Screen prompts and diagnostics
      maxname    = 132             ! Limit on file name lengths

   character (11), parameter :: &
      format = 'unformatted'

!  Variables:

   integer :: &
      i1, ib, ios, l, m, mi, mj, mk, n, nblocks, numf_in, numf_out

   integer, allocatable, dimension (:) :: &
      ifun

   logical :: &
      cell_centered, cr, eof, formatted

   character (maxname) :: &
      filename

!  Derived data types:

   type (grid_type), pointer, dimension (:) :: &
      f_in, one_block

!  Execution:

   ios = 1
   do while (ios /= 0)  ! Check existence and format

      filename = 'volume_for_blayer.fu'
      call reads (luncrt, &
         'Input PLOT3D-type function file [<cr> = volume_for_blayer.fu]: ',  &
         lunkbd, filename, cr, eof)
      if (eof) go to 99

      l = len_trim (filename)
      call determine_grid_form (filename(1:l), lunf_in, formatted, ios)

   end do

   i1 = 1;  if (formatted) i1 = 3

   open (lunf_in, file=filename, form=format(i1:11), status='old')

!  Read the function file header.  Not having read a grid first makes trouble:

   if (formatted) then
      read (lunf_in, *, iostat=ios) nblocks
   else
      read (lunf_in,    iostat=ios) nblocks
   end if
   if (ios /= 0) go to 99

   rewind (lunf_in)
   allocate (f_in(nblocks))

   call q_header_io (1, lunf_in, formatted, nblocks, numf_in, f_in, ios)
   if (ios /= 0) go to 99

!  Open and read the list of output function numbers:

   allocate (ifun(numf_in))

   numf_out = numf_in  ! Max. length is numf

   open (lunctl, file='extract_functions.inp', status='old', iostat=ios)
   if (ios /= 0) then
      write (luncrt, '(a)') '*** Cannot find extract_functions.inp'
      go to 99
   end if

   call rdlist (luncrt, ' ', lunctl, numf_out, ifun)
   if (numf_out <= 0) then
      write (luncrt, '(a, i6)') &
         '*** Output function list is bad.  numf_out:', numf_out
      go to 99
   end if

!  Get the output file name, open the file, and write its header:

   filename = 'extract.fu'
   if (formatted) filename(10:10) = ' '

   call reads (luncrt, 'Output file name: ', lunkbd, filename, cr, eof)

   open (lunf_out, file=filename, form=format(i1:11), status='unknown')

   call q_header_io (2, lunf_out, formatted, nblocks, numf_out, f_in, ios)
   if (ios /= 0) go to 99

!  Process one block at a time:

   allocate (one_block(1))

   do ib = 1, nblocks

      call f_allocate (f_in(ib), numf_in, ios)
      if (ios /= 0) go to 99

      mi = f_in(ib)%mi;  mj = f_in(ib)%mj;  mk = f_in(ib)%mk
      call f_block_io (1, lunf_in, formatted, numf_in, mi, mj, mk, &
                       f_in(ib)%q, ios)
     if (ios /= 0) go to 99

      one_block(1)%mi = mi
      one_block(1)%mj = mj
      one_block(1)%mk = mk

      call f_allocate (one_block(1), numf_out, ios)
      if (ios /= 0) go to 99

      do n = 1, numf_out
         m = ifun(n)
         one_block(1)%q(:,:,:,n) = f_in(ib)%q(:,:,:,m)
      end do

      call f_block_io (2, lunf_out, formatted, numf_out, mi, mj, mk, &
                       one_block(1)%q, ios)
      if (ios /= 0) go to 99

      deallocate (one_block(1)%q, f_in(ib)%q)
   end do

   close (lunf_in)
   close (lunf_out)

99 continue  ! Avoid system-dependent STOP behavior.

   end program extract_functions
