!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine read_to_line_feed (lun, buffer, nchars, ios)

!        Read character data one byte at a time from the current position in an
!     unformatted file opened with access='stream' until a line feed character
!     is detected.  The characters preceding that line feed are returned packed
!     as buffer(1:nchars).  If no line feed is found, nchars is returned as the
!     negative of the buffer length and it is up to the application to deal with
!     this possibility (or not).  Using a buffer more than long enough is highly
!     recommended.
!
!        This utility was prompted by the need to read binary Virtual Toolkit
!     (VTK) legacy files written by the flow solver US3D.
!
!  History:
!
!     07/07/2011  D.A.Saunders  The author's first of what might be a growing
!                               set of stream I/O utilities prompted by Intel's
!                               extensions to the Fortran 95 standard.
!     07/08/2011    "     "     Blank the buffer in case part of a prior call
!                               shows up in an error message.
!
!  Author:  David Saunders, ERC, Inc./NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer,   intent (in)  :: lun          ! Logical unit # to read from
   character, intent (out) :: buffer*(*)   ! For packing the characters that
                                           ! precede the first line feed read
   integer,   intent (out) :: nchars       ! The number of characters found,
                                           ! not counting the line feed;
                                           ! nchars < 0 means no line feed
                                           ! was found before the buffer was
                                           ! filled; nchars is returned as
                                           ! -len (buffer) in this case
   integer,   intent (out) :: ios          ! ios = 0 if no read error occurred;
                                           ! ios > 0 if  a read error occurred;
                                           ! ios < 0 if end-of-file was read

!  Local constants:

   character, parameter :: LF*1 = char (10)

!  Local variables:

   integer   :: i, limit
   character :: c*1

!  Execution:

   nchars = 0
   limit  = len (buffer)
   buffer = ' '

   do i = 1, limit + 1
      read (lun, iostat=ios) c
      if (ios /= 0)        exit
      if (c  == LF)        exit  ! Deals with nchars = 0 & nchars = limit cases
      if (i > limit) then
         nchars = -limit;  exit
      end if

      nchars = i
      buffer(i:i) = c
   end do

   end subroutine read_to_line_feed
