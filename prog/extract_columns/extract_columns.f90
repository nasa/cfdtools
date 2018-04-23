!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program extract_columns
!
!  PURPOSE:
!
!     Extract specified columns from a text file consisting of multiple
!     columns of floating point numbers, possibly needed scaling and/or
!     shifting.  The scaling is of the form x' <- scale * x + shift.
!
!     Enhancement:  This version performs the following transformation:
!
!        x' <- scale * (x ** power) + shift   where power is usually 1.
!
!     This version handles embedded text lines or lines with fewer than
!     the specified number of columns (possibly header lines of integers)
!     by either transmitting them intact or suppressing them as requested.
!
!     Converting to absolute values is done with a kludge:  use a negative
!     column number.
!
!     Simpler Option:  Provide for rearranging columns interactively if the
!                      control file is not found, in which case none of the
!                      shifting/scaling is attempted.
!
!  SAMPLE CONTROL FILE             ! Control file name: 'extract_columns.inp'
!
!     INPUT FILE        OUTPUT FILE
!     traj.plt          SHARP-V7.trajectory
!
!     # COLUMNS INPUT   # COLUMNS OUTPUT   TRANSMIT ODD LINES?
!     26                4                  T
!
!     COLUMN NUMBERS TO EXTRACT
!     1  9  6           25
!
!     POWERS
!     1. 1. 1. 1.
!
!     SCALES
!     1. 1. 0.020885472 1.
!
!     SHIFTS
!     0. 0. 0.          0.
!
!     OUTPUT FORMAT                ! Exactly as for F90, with parens.
!     (F5.0, F9.5, 1P, E10.3, E15.7)
!
!
!  HISTORY:
!     
!     07/16/01  D.Saunders  Initial implementation for dealing with
!                           trajectory results from Traj.
!     08/23/01    "   "     Transmit or suppress text lines or lines that
!                           generate errors when read as "ncol_out" reals.
!     09/27/02    "   "     Mission outputs required 500+ chars/line and
!                           ~40 columns input.  Kludged a way of ignoring
!                           the output sign via a negative column input.
!     09/10/03    "   "     Introduced POWERS to enable conversion of
!                           temperatures (F) to heat fluxes (W/cm^2),
!                           albeit in two runs of this utility.
!     03/21/11   "    "     Added an interactive action (no control file found)
!                           for such simple requirements as reversing the order
!                           of a pair of columns.  Allow for extracting any
!                           subset of columns in any order, interactively.
!
!  AUTHOR:  David Saunders, ELORET Corp./NASA Ames Research Center, CA
!                           (Now ERC, Inc./ARC)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Constants:

   integer, parameter :: &
      lunin = 1, lunout = 2, lunctl = 3, lunkbd = 5, luncrt = 6, ncol_max = 40

   real, parameter :: &
      ONE = 1.

!  Variables:

   integer :: &
      i, ios, j, jcol, nchar, ncol_in, ncol_out, nform, nodd

   integer, dimension (ncol_max) :: &
      icol_out

   real, dimension (ncol_max) :: &
      power, result, scale, shift, value

   logical :: &
      interactive, odd, transmit

   character :: &
      buffer * 512, filein * 128, fileout * 128, format_string * 132

!  Procedure:

   logical  :: alpha ! String function for identifying alphabetic text
   external :: alpha

!  Execution:

   open (lunctl, file='extract_columns.inp', status='old', iostat=ios)

   interactive = ios /= 0

   if (interactive) then

      call simple_case ()  ! Internal procedure below

      go to 99  ! Single STOP philosophy

   end if

   read (lunctl, *)
   read (lunctl, *) filein, fileout

   read (lunctl, *)
   read (lunctl, *)
   read (lunctl, *) ncol_in, ncol_out, transmit

   read (lunctl, *)
   read (lunctl, *)
   read (lunctl, *) icol_out(1:ncol_out)

   read (lunctl, *)
   read (lunctl, *)
   read (lunctl, *) power(1:ncol_out)

   read (lunctl, *)
   read (lunctl, *)
   read (lunctl, *) scale(1:ncol_out)

   read (lunctl, *)
   read (lunctl, *)
   read (lunctl, *) shift(1:ncol_out)

   read (lunctl, *)
   read (lunctl, *)
   read (lunctl, '(a)') format_string

   nform = len_trim (format_string)

   close (lunctl)

   open  (lunin,  file=filein,  status='old')
   open  (lunout, file=fileout, status='new')

   nodd = 0

   do ! Until EOF

      read  (lunin, '(a)', iostat=ios) buffer
      if (ios < 0) exit

      nchar = len_trim (buffer)

      odd = alpha (buffer(1:nchar)) ! T probably means alphabetic text

      if (.not. odd) then ! Try to read as columns

         read  (buffer(1:nchar), *, iostat=ios) value(1:ncol_in)

         odd = ios /= 0 ! Too few columns, or un-number-like text

      end if

      if (odd) then
         nodd = nodd + 1
         if (transmit) write (lunout, '(a)') buffer(1:nchar) ! Else suppress
      else
         do i = 1, ncol_out
            j = icol_out(i)
            jcol = abs (j)

            if (power(i) == ONE) then
               result(i) = value(jcol) * scale(i) + shift(i)
            else
               result(i) = &
                  (value(jcol) ** power(i)) * scale(i) + shift(i)
            end if

            if (j < 0) result(i) = abs (result(i))  ! Kludge for |accln.| in Gs
         end do

         write (lunout, format_string(1:nform)) result(1:ncol_out)
      end if

   end do

   close (lunin)
   close (lunout)

   if (nodd > 0) then
      if (transmit) then
         write (*, '(/, a, i6, a, /)') &
            ' Note: ', nodd, ' lines were transmitted as plain text.'
      else
         write (*, '(/, a, i6, 2a, /)') &
            ' Note: ', nodd, ' lines were suppressed as short or non-numeric.'
      end if
   end if

99 continue
 
!  Internal procedure:

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine simple_case ()

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Simply prompt for a list of columns in the desired output order.
!     Omit shifting/scaling/reformatting; transmit columns verbatim.
!     Rows with unexpected numbers of columns will cause an abort.

!     Local constants:

      integer, parameter :: nchars_max = 30  ! Limit on width of any column

!     Local variables:

      integer :: i, line, nwords
      character, allocatable, dimension (:) :: words*(nchars_max)

!     Execution:

      ios = 1
      do while (ios /= 0)
         write (luncrt, '(a)', advance='no') 'Input file: '
         read  (lunkbd, '(a)', iostat=ios) filein
         if (ios < 0) go to 99
         open  (lunin, file=filein, status='old', iostat=ios)
      end do

      write (luncrt, '(a)', advance='no') 'Expected number of input columns: '
      read  (lunkbd, *) ncol_in
      ncol_out = ncol_max

      call rdlist (luncrt, 'Output column(s) in desired order:', lunkbd, &
                   ncol_out, icol_out)
      if (ncol_out <= 0) go to 99

      write (luncrt, '(a)', advance='no') 'Output file: '
      read  (lunkbd, '(a)', iostat=ios) fileout
      if (ios < 0) go to 99
      open  (lunout, file=fileout, status='unknown', iostat=ios)

      allocate (words(ncol_in))

      line = 0
      do ! Until EOF

         read  (lunin, '(a)', iostat=ios) buffer
         if (ios < 0) exit

         line = line + 1
         nwords = ncol_in

         call tokens (buffer, nwords, words)

         if (nwords /= ncol_in) then
            write (luncrt, '(a, 2i10, /, a, i10)') &
               'Unexpected number of columns found:', nwords, ncol_in, &
               'Offending line:', line
            go to 99
         end if

         write (lunout, '(40(1x,a))') &
            (trim (words(icol_out(i))), i = 1, ncol_out)

      end do  ! Next input line

      close (lunin)
      close (lunout)

99    return

      end subroutine simple_case

   end program extract_columns
