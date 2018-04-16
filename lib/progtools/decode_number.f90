!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine decode_number (string, inumber, rnumber, ios)

!        DECODE_NUMBER modularizes internal reading of strings that supposedly
!     represent integers or reals.  The calling program should check the output
!     status argument (ios) first in case of a bad input.  If it's OK, the
!     context should allow picking inumber or rnumber as the value to use.
!     If a string could be interpreted as real as well as integer, both those
!     arguments will be set as a matter of convenience, at negligible cost
!     since the "real" intrinsic can be used rather than a second internal read.
!
!  Limitations:
!
!         Applications with mixed floating-point precision numbers may not be
!     able to use this utility.  Compiler switches -r8 or -r4 (or equivalent)
!     are assumed to be employed.  Complex numbers are not treated either.
!
!         If a compiler does not handle list-directed internal reads, it would
!     be straightforward to construct a format string here based on the length
!     of the string.
!
!  History:
!
!     23 Jan. 2009  D. Saunders  One more string utility to save coding lots of
!                                internal reads in a free-form input scheme.
!
!  Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   character, intent (in)  :: string * (*)   ! String to be decoded
   integer,   intent (out) :: inumber        ! Integer value to use (see ios)
   real,      intent (out) :: rnumber        ! Real value to use (see ios)
   integer,   intent (out) :: ios            ! 0 => neither integer nor real
                                             ! 1 => real (only): inumber not set
                                             ! 2 => integer or real (both set)
!  Local constants:

   character, parameter :: not_in_integer * 5 = '.dDeE'

!  Local variables:

   integer :: in, is
   logical :: not_integer

!  Execution:

   not_integer = .false.

   do is = 1, len (string)

      do in = 1, 5
         if (string(is:is) == not_in_integer(in:in)) not_integer = .true.
      end do

   end do

   if (not_integer) then
      read (string, *, iostat=ios) rnumber
      if (ios == 0) then
         ios = 1
      else
         ios = 0  ! Not real either - must be a bad string
      end if
   else
      read (string, *, iostat=ios) inumber
      if (ios == 0) then
         ios = 2  ! Could still be real
         rnumber = real (inumber)
      else
         ios = 0  ! Bad input
      end if
   end if

   end subroutine decode_number
