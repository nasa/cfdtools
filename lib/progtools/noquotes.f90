!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine noquotes (string, length)

!        Left-justify the input character string in-place, and suppress possible
!     leading and trailing single or double quotes.  Return the adjusted length.
!     Any leading quote is assumed to be matched by a trailing quote.
!
!        This was prompted by misbehavior of list-directed reading of file names
!     when quotes are omitted:  a slash (/) behaves as end-of-information, so
!     reading with '(a)' format is needed, and that does not left-justify the
!     way list-directed reading does.
!
!  History:
!
!     04 Nov. 2009  D. Saunders  Yet another string utility to work around
!                                odd behavior of list-directed reading of
!                                unquoted path names.
!
!  Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   character, intent (inout)  :: string * (*)   ! String to be adjusted in-place
   integer,   intent (out)    :: length         ! Adjusted trimmed length

!  Local variables:

   integer :: l

!  Execution:

   string = adjustl (string)  ! Left justify

   l = len_trim (string)

   if (string(1:1) == '''' .or. string(1:1) == '"') then
       string(1:l-2) = string(2:l-1)
       string(l-1:)  = ' '
       l = l-2
   end if

   length = l

   end subroutine noquotes
