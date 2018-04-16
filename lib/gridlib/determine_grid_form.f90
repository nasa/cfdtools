!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine determine_grid_form (grid_file_name, lun, formatted, ios)
!
!  Open a multiblock grid or function file (PLOT3D-type) and attempt to read the
!  number of blocks, thus determining whether it seems to be unformatted or not.
!  Close the file upon return.  Nonexistence of the named file is also trapped.
!
!  Note that use of the INQUIRE statement with the FORM or FORMATTED keyword
!  does not solve the problem because unless the file is open already, the
!  result is 'UNDEFINED' or 'UNKNOWN'.  This seems to be a Catch-22 situation.
!
!  Using this utility in a while loop checking ios allows retries in case of bad
!  file names as shown below.  Here, "reads" is part of a prompting utility
!  "reader" available from the present author, "cr" means "carriage return" or
!  "default", and entering control-D is the standard "end-of-file" or "quit"
!  indicator under Unix):
!
!     integer,   parameter :: lunin = 1, lunkbd = 5, luncrt = 6
!     character, parameter :: format * 11 = 'unformatted'
!
!     integer   :: i1, ios
!     logical   :: cr, eof, formatted
!     character :: filename * 64
!
!     :::::::
!
!     ios = 1
!     do while (ios /= 0)  ! Check existence and format
!
!        filename = 'both-halves.cvertex.gu'
!        call reads (luncrt, &
!                    'Baseline grid [<CR> = both-halves.cvertex.gu]: ',  &
!                    lunkbd, filename, cr, eof)
!        if (eof) go to 99  ! Abort, or whatever
!
!        call determine_grid_form (filename, lunin, formatted, ios)
!
!     end do
!
!     i1 = 1;  if (formatted) i1 = 3
!
!     open (lunin, file=filename, form=format(i1:11), status='OLD')
!
!     :::::::::
!
!  History:
!
!     Nov. 2005  D.A.Saunders  Original "determine_form" internal procedure as
!                              part of prepare_rapid_analysis.f90.
!     01/23/2008    "     "    Reusable "determine_grid_form" utility.
!     10/24/2016    "     "    A certain CAD-related file had 10493 blocks!
!
!  Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   character, intent (in)  :: grid_file_name * (*)
   integer,   intent (in)  :: lun
   logical,   intent (out) :: formatted
   integer,   intent (out) :: ios      ! ios /= 0 suggests a bad file name

!  Local constants:

   integer, parameter :: nb_suspicious = 20000 ! Trap questionable block count

!  Local variables:

   integer :: l, nbf, nbu
   logical :: doubtful

   l = len_trim (grid_file_name)

!  First, try opening the file as formatted (the default):

   open (lun, file=grid_file_name(1:l), status='old', iostat=ios)

   if (ios /= 0) then
      write (*, '(2a)') ' Unable to open ', grid_file_name(1:l)
      go to 99
   end if

   read (lun, *, iostat=ios) nbf

   doubtful = ios /= 0
   if (.not. doubtful) doubtful = nbf <= 0 .or. nbf > nb_suspicious

   if (doubtful) then
!!!   write (*,*) 'Formatted read gives nb & ios = ', nbf, ios
!!!   write (*,*) 'Try again as unformatted.'
      close (lun)
      open (lun, file=grid_file_name(1:l), status='old', form='unformatted')
      formatted = .false.
      read  (lun, iostat=ios) nbu
!!!   write (*, *) 'unformatted read gives nb & ios = ', nbu, ios
      if (ios == 0) then
         write (*, *) 'The file appears to be unformatted.'
      else
         write (*, *) &
         '# blocks found from formatted and unformatted reads: ', nbf, nbu, &
         'Bad file name?'
      end if
   else
      formatted = .true.
      write (*, *) 'The file appears to be formatted.'
!!!   write (*, *) 'nb, ios:', nbf, ios
   end if

   close (lun)

99 return

   end subroutine determine_grid_form
