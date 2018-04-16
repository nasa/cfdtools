!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine token_count (string, seps, number)

!        TOKEN_COUNT is a reduced form of TOKEN2 for the case of determining
!     the number of items in an indefinite list so that they can then be read
!     simply via a list-directed internal read.  The application needs to know
!     the data type represented by the indefinite list (integers? reals?), but
!     at this level the items are just character tokens.
!
!  History:
!
!     15 Jan. 2009  D. Saunders  TOKEN_COUNT adapted from TOKEN4, as prompted
!                                by its handling of generalized tokens that
!                                are really delimited lists.
!
!  Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   character, intent (in)  :: string * (*)   ! String to be parsed for tokens
   character, intent (in)  :: seps   * (*)   ! Delimiter(s) for simple tokens
   integer,   intent (out) :: number         ! Number of tokens found

!  Local variables:

   integer  :: first, last, mark

!  Procedures:

   external :: scan2

!  Execution:

   first  = 1
   last   = len (string)  ! The string is most likely a substring in the calling
   number = -1            ! so len_trim is redundant
   mark   = 1

   do while (mark > 0)

      number = number + 1

!     Find the delimiting indices of the next token, if any.

      call scan2 (string, seps, first, last, mark)

      first = mark + 2

   end do

   end subroutine token_count
