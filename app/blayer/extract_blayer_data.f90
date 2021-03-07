!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program extract_blayer_data
!
!  Outline:
!
!     This utility automates tabulation of certain flow solver results for CFD
!     points along a trajectory in readiness for the curve fits and padding to
!     full time histories as needed for TPS sizing at some body point(s).
!
!     For one body point indicated by (iblock, i[, j]), and any number of BLAYER
!     output files, extract the value(s) indicated by a list of variable numbers
!     in any order and all on one line.  Here, x and y count as the first two
!     variables, and s or z is the third, for 2-D or 3-D data respectively.
!
!     The dimensionality is determined automatically from the first BLAYER file
!     read.  All should be the same dimension, and will presumably be entered in
!     some sensible order.  Results are written to standard output.
!
!     This version also looks for tauwx, tauwy[, tauwz] and, if present, adds
!     |tauw| on the end of the requested list of variables.
!
!  Control file example:
!
!     1 163                         ! (iblock, i) for a 2-D case
!     5 26 6 58 32 7 27 28          ! pw, qconv, Tw, CH, Hedge, Hw, tauwx, tauwy
!     MHL-t52.0/blayer.dat          ! Any number of BLAYER output files, to EOF
!     MHL-t55.3/blayer.dat
!     MHL-t57.9/blayer.dat
!     MHL-t61.2/blayer.dat
!     MHL-t64.2/blayer.dat
!     MHL-t65.8/blayer.dat
!     MHL-t70.0/blayer.dat
!
!  Corresponding output (with variable names from the first BLAYER file):
!
!     pw (Pa)          qw (W/m^2)       Tw (K)           ......  tauwy (Pa)
!     1.015842438E+04  4.485145706E+06  3.107182014E+03  ......  8.951380010E-30
!     3.740018736E+04  7.853650464E+06  3.574295442E+03  ......  2.430156551E-29
!     9.548843572E+04  1.160357837E+07  3.940674287E+03  ......  5.841697623E-29
!     2.367064900E+05  1.546494660E+07  4.234085869E+03  ......  1.423226761E-28
!     3.591160591E+05  1.251104156E+07  4.015558173E+03  ......  2.195770009E-28
!     3.768805570E+05  9.299941983E+06  3.728572697E+03  ......  2.362666022E-28
!     2.860373898E+05  2.720219644E+06  2.742038501E+03  ......  2.108604270E-28
!
!  Actually, the variable names now have embedded blanks removed:
!
!     pw,Pa            qw,W/m^2         Tw,K             ......  tauwy,Pa
!
!  History:
!
!     08/06/2014  D.A.Saunders  Initial implementation.
!     09/04/2014    "     "     Convert variable names to single-token names
!                               with MERGE_TABLES in mind.
!     09/09/2014    "     "     |tauw| is now added as a last column if the
!                               shear stress components are present.
!
!  Author:  David Saunders, ERC, Inc./NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_header_structure      ! See Tecplot_io.f90 for all of these modules
   use grid_block_structure
   use tecplot_io_module

   implicit none

!  Constants:

   integer, parameter :: &
      lunin      = 1,    &        ! For reading the BLAYER file(s)
      lunctl     = 5,    &        ! Input control file
      luncrt     = 6,    &        ! Standard output
      name_limit = 32,   &        ! Must match the value in Tecplot_io_module
      maxnv      = 128,  &        ! Limit on no. of variables/values to extract
      maxsize    = 9999, &        ! For trapping a likely bad surface index
      width      = 17             ! If we use es17.9 format for output columns
   character (1), parameter :: &
      blank      = ' '

!  Variables:

   integer :: &
      i, i1, i2, indices(3), ios, isurf, ivnumbers(maxnv), izone, j, jsurf, l, &
      n, ndim, nindices, nfile, nv, nvextract, nvwrite
   integer :: &
      indices_tau(3)
   real :: &
      rnumber, tau, values(maxnv) ! x/y/z aren't in xyzq(:)%q
   logical :: &
      append_tau
   character (4) :: &
      vnumbers(maxnv)             ! For tokenizing the variable numbers
   character (4) :: &
      tau_name
   character (132) :: &
      buffer
   character (maxnv*width) :: &
      column_headers

!  Derived data types:

   type (grid_header) :: &
      header
   type (grid_type), pointer, dimension (:) :: &
      xyzq

!  Procedures:

   logical :: number  ! Determines if a string is a number or not

!  Execution:
!  !!!!!!!!!!

   write (luncrt, '(a)')

!  We don't want to backspace or rewind standard input, so each control file
!  line must be handled in order.

   read (lunctl, '(a)') buffer  ! iblock, isurf[, jsurf][  ! Trailing comment]

   nindices = 3;  read (buffer, *, iostat=ios) indices(1:3)

   if (ios /= 0) then
      nindices = 2;  read (buffer, *, iostat=ios) indices(1:2)
      if (ios /= 0) then
         write (luncrt, '(a)') 'Trouble reading 2 or 3 indices on line 1:', &
           trim (buffer)
         go to 99
      end if
   else  ! Allow for list-directed reading not giving an error for a non-integer
      if (indices(3) <= 0 .or. indices(3) > maxsize) nindices = 2
   end if

!  For the list of variable numbers, we also allow for a trailing comment:

   read (lunctl, '(a)') buffer  ! E.g.:   5  26  6  ... 28  ! pw, qw, ..., tauwy

   nv = maxnv
   call token2 (buffer, blank, nv, vnumbers)  ! Gives an array of strings

!  Suppress any trailing text:

   do i = 1, nv
      nvextract = i
      call decode_number (vnumbers(i), ivnumbers(i), rnumber, ios)
      if (ios /= 2) then  ! Not an integer
         nvextract = i - 1
         exit
      end if
   end do

!! write (luncrt, '(a, i4)') ' Number of variables found to extract:', nvextract

!  Process the list of BLAYER output files until EOF:

   header%formatted = .true.  ! Not sure how to distinguish a Tecplot *.plt file
   header%ndim      = -1      ! Tells tecplot_io to determine ndim
   nfile = 0

   do  ! Until EOF (or error)

      read (lunctl, '(a)', iostat=ios) buffer
      if (ios /= 0) exit

      buffer = adjustl (buffer)      ! In case of leading blanks
      l = index (buffer, blank) - 1
      header%filename = buffer(1:l)  ! First token only
      nfile = nfile + 1

!!    ios = 1  ! Verbose mode (mucks up the desired table on standard output)

      call tecplot_read (lunin, header, xyzq, ios)  ! Easiest to read whole file
      if (ios /= 0) exit

      if (nfile == 1) then  ! Print the variable names (fixed-width columns)
         if (header%ndim == 3) then
            if (nindices /= 3) then
               write (luncrt, '(a)') &
                  'BLAYER file is 3-D, but only (iblock,i) was read.'
               go to 99
            end if
            jsurf = indices(3)
         else
            jsurf = 1
            nindices = 2
         end if

         ndim  = nindices
         izone = indices(1)
         isurf = indices(2)
         column_headers = blank
         i2 = 2  ! Moves the names over 2 spaces, where positive values start
         do i = 1, nvextract
            i1 = i2 + 1
            i2 = i*width + 2
            column_headers(i1:i2) = header%varname(ivnumbers(i))
         end do

         call look_for_tau ()  ! Option to append |tauw| or |taue|

         if (append_tau) then
            i1 = i2 + 1
            i2 = i2 + 9
            column_headers(i1:i1+3) = tau_name
            column_headers(i1+4:i2) = ' (Pa)'
         end if

         call noblanks (column_headers(1:i2))  ! Single-token column headers

         write (luncrt, '(a)') trim (column_headers)
      end if

      do i = 1, nvextract
         j = ivnumbers(i)
         if (j == 1) then
            values(i) = xyzq(izone)%x(isurf,jsurf,1)
         else if (j == 2) then
            values(i) = xyzq(izone)%y(isurf,jsurf,1)
         else if (j == 3 .and. ndim == 3) then
            values(i) = xyzq(izone)%z(isurf,jsurf,1)
         else
            values(i) = xyzq(izone)%q(j-ndim,isurf,jsurf,1)
         end if
      end do

      nvwrite = nvextract

      if (append_tau) then
         tau = 0.
         do n = 1, ndim
            j = indices_tau(n)
            tau = tau + xyzq(izone)%q(j-ndim,isurf,jsurf,1)**2
         end do
         tau = sqrt (tau)
         nvwrite = nvextract + 1
         values(nvwrite) = tau
      end if

      write (luncrt, '(128es17.9)') values(1:nvwrite)

      deallocate (header%varname, stat=ios)
      if (ios /= 0) go to 99

      call deallocate_blocks (1, header%nblocks, ndim, header%numq, xyzq, ios)
      if (ios /= 0) go to 99

      deallocate (xyzq, stat=ios)
      if (ios /= 0) go to 99

   end do  ! Next BLAYER file

99 continue

!  Internal procedure for program EXTRACT_BLAYER_DATA

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine look_for_tau ()  ! Determine whether tau?x/y[,z] are present.
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: i, i1, j, n

!     Are there ndim occurrences of "tau" in the variable names?

      n = 0
      do i = 1, nvextract
         j = ivnumbers(i)
         i1 = index (header%varname(j), 'tau')
         if (i1 > 0) then
            n = n + 1
            indices_tau(n) = j
            if (n == ndim) then
               tau_name = header%varname(j)(1:4)  ! tauw | taue; not both
               exit
            end if
         end if
      end do

      append_tau = n == ndim

      if (append_tau) then
         if (nvextract == maxnv) then
            append_tau = .false.
            write (luncrt, '(a, i4)') &
               '*** No room to append tau variable; maxnv:', maxnv
         end if
      end if

      end subroutine look_for_tau

   end program extract_blayer_data
