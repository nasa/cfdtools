!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine ndigits (number, length)

!  Description:
!
!     For the given integer, determine its number of digits as needed in order
!     to print the integer with a run-time-generated format.  For instance,
!     the lengths of the numbers 1001 and -999 are both returned as 4.
!
!  History:
!
!     05/19/08  D.A.Saunders  Initial implementation, for writing Tecplot files.
!
!  Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer,   intent (in)    :: &
      number                    ! Integer, possibly negative

   integer,   intent (out)   :: &
      length                    ! # digits, possibly including a minus sign

!  Local variables:

   integer   :: l, len
   character :: buffer * 16

!  Execution:

!  Take the easy way out, assuming list-directed internal writes are legal:

   write (buffer, *) number
   length = len_trim (buffer)

!  It should be left-justified, but just in case:

   len  = length
   do l = 1, len
      if (buffer(l:l) == ' ') length = length - 1
   end do

   end subroutine ndigits
