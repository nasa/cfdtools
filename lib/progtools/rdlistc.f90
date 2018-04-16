!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine rdlistc (lun, delimiters, mode, maxperline, nstrings, strings)

!  Description:
!
!     Read an indefinite list of character strings from the indicated file,
!     using the indicated string delimiter(s).
!
!     Three likely situations are handled according to the "mode" argument:
!
!        mode = 1 means all strings are expected to be on a single next line;
!        mode = 2 means one or more lines of input terminated by a blank line
!                 or by a line beginning with "END" (case-insensitive),
!                 or by end-of-file.
!
!     The blank line option and the "maxperline" argument were forced by the
!     initial application to NEQAIR.  The "end of input" option may be handy
!     for other applications.
!
!  History:
!
!     05/09/08   D.A.Saunders   Initial implementation, for enhancing NEQAIR.
!
!  Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer,   intent (in)    :: &
      lun                       ! Logical unit number, ready to read from

   character, intent (in)    :: &
      delimiters * (*)          ! String separator(s) such as ' ' or ', ' or
                                ! perhaps including the tab character
   integer,   intent (in)    :: &
      mode                      ! See above; tells how to end reading

   integer,   intent (in)    :: &
      maxperline                ! Allows trailing comments to be ignored;
                                ! use a large number if there's no real limit
   integer,   intent (inout) :: &
      nstrings                  ! Input with the maximum allowed for in the
                                ! calling program; output with the number found;
                                ! 0 <= output nstrings <= input nstrings
   character, intent (out)   :: &
      strings (*) * (*)         ! Strings found, left-justified and possibly
                                ! truncated if longer than allowed for
!  Local variables:

   integer   :: first, fplus2, ios, last, lenbuf, mark, nfound, nthisline
   character :: buffer * 132

!  Execution:

   nfound = 0

   read (lun, '(a)', iostat=ios) buffer

   do while (ios == 0)

      lenbuf = len_trim (buffer)

      if (lenbuf == 0) exit    ! Blank line termination - evidently not mode 1

!     Isolate tokens, but no more than maxperline on this line:

      nthisline = 0;  first = 1;  last = lenbuf

      do while (first     <= last       .and.  &
                nthisline <  maxperline .and.  &
                nfound    <  nstrings)

         call scan2 (buffer, delimiters, first, last, mark)

         nthisline = nthisline + 1;  nfound = nfound + 1

         strings(nfound) = buffer(first:mark)  ! Pads with blanks or truncates

         if (mode == 2) then
            if (nthisline == 1) then  ! Is this an "end of data" marker?
               fplus2 = first + 2
               if (fplus2 <= mark) then
                  call upcase (buffer(first:fplus2))
                  if (buffer(first:fplus2) == 'END') then
                     nfound = nfound - 1
                     ios = 1  ! Don't want to keep reading
                     exit
                  end if
               end if
            end if
         end if

         first = mark + 2

      end do

      if (mode == 1) exit

      if (ios == 0) read (lun, '(a)', iostat=ios) buffer

   end do

   nstrings = nfound

   end subroutine rdlistc
