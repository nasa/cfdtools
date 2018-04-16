!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine numbered_name_2d (prefix1, number1, prefix2, number2, name, &
                                length)
!
!  This utility applies the earlier numbered_name to construct a name from an
!  integer pair.  It was prompted by a need for unique names defined by grid
!  (block, face) pairs, such as 'Block 11 Face 6'.  Embedded blanks as in this
!  example are permitted.  Here, they would be obtained with a call like this:
!
!    call numbered_name_2d ('Block ', ib, ' Face ', iface, group_name, length)
!
!  The following description from numbered_name applies similarly:
!
!  The character argument name is assumed to be long enough for concatenating
!  the prefixes with any default (4-byte) integers that the calling program is
!  expected to encounter.  NO ERROR HANDLING IS PERFORMED.  Rather, the largest
!  possible number (in magnitude, positive or negative) is written by each call
!  to numbered_name if the name argument is not long enough.
!
!  Upon return, name(1:length) could be written with ('a') format as an
!  alternative to constructing a run-time format from just the active length.
!  The active string may also be needed for other purposes (see numbered_name).
!
!  11/26/12  D.A.S.  Adaptation of numbered_name as described above.
!
!  Author:  David Saunders, ERC, Inc. at NASA Ames Research Center, California.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   character (*), intent (in)  :: prefix1  ! 1st character portion of the name
   integer,       intent (in)  :: number1  ! 1st integer, possibly negative
   character (*), intent (in)  :: prefix2  ! 2nd character portion of the name
   integer,       intent (in)  :: number2  ! 2nd integer, possibly negative
   character (*), intent (out) :: name     ! Desired name, left justified
   integer,       intent (out) :: length   ! name(1:length) is significant;
                                           ! any trailing characters are blank
!  Local variables:

   integer :: length1, length2

!  Execution:

   call numbered_name (prefix1, number1, name, length1)
   call numbered_name (prefix2, number2, name(length1+1:), length2)

   length = length1 + length2

   end subroutine numbered_name_2d
