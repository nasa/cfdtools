!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine decode_rlist (string, nlist, rlist, ier)
!
!  Description:
!
!     From the given character string, expected to contain an indefinite list
!     of reals, determine how many there are (up to the limit indicated by the
!     input value of nlist) and return the corresponding real array.  If a bad
!     number is detected, abort with ier /= 0.
!
!     Token separators may be any number of blanks, commas, or tabs.
!
!  History:
!
!     09/03/2014  D.A.Saunders  Simultaneous adaptation of DECODE_ILIST.
!
!  Author:  David Saunders, ERC, Inc./NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   character (*), intent (in)    :: string    ! String to be decoded
   integer,       intent (inout) :: nlist     ! Input  with # allowed for;
                                              ! output with # integers found,
                                              ! not exceeding the input #
   real,          intent (out)   :: rlist(*)  ! Desired list of reals
   integer,       intent (out)   :: ier       ! 0 => no problem detected;
                                              ! 1 => bad number in the list
!  Local variables:

   integer       :: first, inumber, last, mark, nmost
   real          :: rnumber
   character (3) :: seps

!  Procedures:

   external :: scan2          ! Isolates the next token in a string
   external :: decode_number  ! Converts character token to integer|real|both

!  Execution:

   seps  = ' ,' // char (9)  ! Blank, comma, or tab delimiters
   nmost = nlist  
   nlist = 0
   ier   = 0
   first = 1
   last  = len (string)

   do  ! Until nmost good values have been found, or fewer if the list is short

      call scan2 (string, seps, first, last, mark)  ! Isolate the next token

      if (mark == 0) exit  ! End of string

      call decode_number (string(first:mark), inumber, rnumber, ier)

      if (ier == 0) then  ! 1 => real (only); 0 => neither integer nor real
          ier = 1         ! 2 => integer or real
          exit
      end if

      ier = 0
      nlist = nlist + 1
      rlist(nlist) = rnumber
      first = mark + 2
      if (nlist >= nmost) exit

   end do

   end subroutine decode_rlist
