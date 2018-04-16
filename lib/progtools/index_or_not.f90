!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine index_or_not (string, indices)

!        INDEX_OR_NOT is prompted by applications involving definition of body
!     points within a surface grid via a list of either integers (indices) or
!     real numbers.  It should be possible to distinguish the type of list
!     automatically, but there are pitfalls.  For instance, the presence of a
!     decimal points (period) in the string is insufficient evidence of a real
!     number list if trailing comments are allowed for.  Therefore, only the
!     first token is treated, and this is isolated here rather than at the
!     higher level (although the string may already be just one token).
!
!        Another pitfall would be entry of "0 0" to mean the nose coordinates
!     of a typical aerospace vehicle in an axisymmetric grid, meaning "0. 0.".
!     Therefore, if an apparent index is "0" it is instead interpreted as a
!     real, whereas "1" or any other positive integer is interpreted as an
!     index.  Always include the decimal point in the string if you want it
!     to be interpreted as a real (with 0 meaning 0. the one exception).
!
!        Tabs rather than blanks could be allowed for, but presently are not.
!
!  History:
!
!     16 Aug. 2014  D. Saunders  One more string utility prompted by the
!                                author's LINES_OF_SIGHT[_2D] programs, where
!                                the desirability of allowing trailing comments
!                                in a list of body points caused a wrong
!                                interpretation when merely the presence of a
!                                period or not was thought to suffice.
!
!  Author:  David Saunders, ERC, Inc./NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   character, intent (in)  :: string * (*)   ! String of which only the first
                                             ! token is treated; blanks, not
                                             ! tabs, are assumed delimiters
   logical,   intent (out) :: indices        ! F if the first token contains a
                                             ! decimal point or just '0';
                                             ! T otherwise; can't use the name
                                             ! "index" because of its use as
                                             ! an intrinsic
!  Local variables:

   integer :: l
   character (32) :: token  ! Surely more than long enough for a real number

!  Execution:

   token = adjustl (string)  ! Handle any leading blank(s)
   l = index (token, ' ') - 1
   if (l < 0) l = len (token)  ! Highly unlikely

   indices = index (token(1:l), '.') == 0

   if (l == 1) then
      if (token(1:1) == '0') indices = .false.  ! Special case (see above)
   end if

   end subroutine index_or_not
