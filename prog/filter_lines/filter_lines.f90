!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program filter_lines
!
!  Description:
!
!     Transcribe any initial non-numeric header lines from a multi-column file
!     (if specified, else ignore them) then filter out any further non-numeric
!     lines.  Optionally, comment the header line(s) and/or append the same
!     line(s) as a trailer.
!
!  Controls (namelist on standard input; sample shown):
!
!     $filter
!     infile = 'freestream.dat'
!     outfile = 'freestream-clean.dat'
!     keep_top_header = T
!     add_bottom_header = F
!     header_comment = '#'
!     $end
!
!  Notes:
!
!     o  This is for people who don't know awk or similar languages.
!     o  Namelist control allows for future extensions without affecting
!        existing control files or scripts.
!     o  If header_comment = ' ', any header lines are left untouched.
!     o  Transcribed numeric lines are copied exactly as is (no reformatting).
!
!  History:
!     08/21/2014  D.A.Saunders  Initial design, prompted by the need to produce
!                               trajectory time histories of high-fidelity flow
!                               data, starting with output from NASA ARC's Traj
!                               program and selected CFD point results.
!                               (Traj repeats its headers every 20 lines, which
!                               is inconvenient when pasting a trajectory into
!                               a spreadsheet.)
!     08/23/2014    "     "     Initial implementation.
!
!  Author:  David Saunders, ERC, Inc./NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Constants:

   integer, parameter :: &
      lunin  = 1,        &    ! Input data file, one or more columns
      lunout = 2,        &    ! Filtered output file
      lunctl = 5,        &    ! Standard input namelist control file
      luncrt = 6              ! Output diagnostics

!  Variables:

   integer          :: i, ios, l, nheader
   logical          :: add_bottom_header, keep_top_header, comment
   character (1)    :: header_comment
   character (80)   :: infile, outfile
   character (2048) :: buffer ! Consistent with INTERP_COLUMNS

!  Namelist controls:

   namelist /filter/ &
      infile, outfile, keep_top_header, add_bottom_header, header_comment

!  Procedures:

   external :: alpha  ! Identifies a string as (probably) text or not
   logical  :: alpha

!  Execution:

   keep_top_header   = .true.
   add_bottom_header = .false.
   header_comment    = '#'

   read (lunctl, nml=filter, iostat=ios)
   if (ios /= 0) then
      write (luncrt, '(a)') 'Trouble reading namelist "filter".'
      go to 99
   end if

   comment = header_comment /= ' '
   open (lunin, file=infile, status='old', iostat=ios)
   if (ios /=0) then
      write (luncrt, '(2a)') 'Cannot find input file ', trim (infile)
      go to 99
   end if

   open (lunout, file=outfile, status='unknown', iostat=ios)

!  Count any header lines:

   nheader = 0
   do  ! Until a numeric line is found
      read (lunin, '(a)', iostat=ios) buffer
      if (ios < 0) then
         if (nheader == 0) then
            write (luncrt, '(a)') 'Input file appears to be empty.'
         else
            write (luncrt, '(a)') 'Input file appears to be all text.'
         end if
         go to 90  ! Add bottom header? Clean up.
      end if
      l = len_trim (buffer)
      if (.not. alpha (buffer(1:l))) exit

      nheader = nheader + 1
      if (keep_top_header) then
         if (comment) buffer(1:1) = header_comment
         write (lunout, '(a)') buffer(1:l)
      end if
   end do

   write (lunout, '(a)') buffer(1:l)  ! First numeric line

   do  ! Until EOF
      read (lunin, '(a)', iostat=ios) buffer
      if (ios < 0) exit

      l = len_trim (buffer)
      if (.not. alpha (buffer(1:l))) write (lunout, '(a)') buffer(1:l)
   end do

90 continue  ! End processing

   if (add_bottom_header) then  ! Preferable to storing unknown # header lines
      rewind (lunin)
      do i = 1, nheader
         read (lunin, '(a)') buffer
         if (comment) buffer(1:1) = header_comment
         write (lunout, '(a)') trim (buffer)
      end do
   end if

   close (lunin)
   close (lunout)

99 continue

   end program filter_lines
