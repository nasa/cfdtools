!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine numbered_name (prefix, number, name, length)
!
!  Convert the given prefix string and integer to a left-justified character
!  string as a numbered name such as 'Block12' or 'Block 12', and return the
!  active length of that name.  The character argument name is assumed to be
!  long enough for concatenating the prefix with any default (4-byte) integer
!  the calling program is expected to encounter.  For instance, if a multiblock
!  grid is certain to contain less than 10,000 blocks, the portion of name that
!  is output with the block number need not be any longer than 4 characters.
!  NO ERROR HANDLING IS PERFORMED.  Rather, the largest possible number (in
!  magnitude, positive or negative) is written if the name argument is not long
!  enough.
!
!  See also the earlier int_to_char subroutine if no prefix is in the picture.
!
!  Upon return, name(1:length) could be written with ('a') format as an
!  alternative to constructing a run-time format from just the active length.
!  The active string may also be needed for other purposes (see History).
!
!  09/28/12  D.A.S.  Initial implementation of int_to_char, for constructing
!                    "group" names from grid block numbers during HDF5 I/O.
!  10/12/12    "     There is no need not to put 'Block' in front of the block
!                    number for the group name, so allow for any prefix in a
!                    more general utility.
!
!  Author:  David Saunders, ERC, Inc. at NASA Ames Research Center, California.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   character (*), intent (in)  :: prefix  ! Character portion of the output name
   integer,       intent (in)  :: number  ! Integer, possibly negative
   character (*), intent (out) :: name    ! Prefix // number, left justified
   integer,       intent (out) :: length  ! name(1:length) is significant;
                                          ! any trailing characters are blank
!  Local variables:

   integer       :: i1, l, nabs, nchar_available, nchar_required
   character (4) :: iformat = '(in)'

!  Execution:

   name = prefix        ! Pads with blanks
   l = len (prefix) + 1 ! Avoid len_trim in case a trailing blank is desired
   nabs = abs (number)

   if (nabs == 0) then  ! We have to avoid the logarithm, so handle it first
      name(l:l) = '0'
      length = l
   else
      nchar_available = len (name) - (l - 1)
      i1 = l
      if (number < 0) then
         nchar_available = nchar_available - 1
         i1 = i1 + 1
         name(i1:i1) = '-'
      end if
      nchar_required = int (log10 (real (nabs))) + 1
      if (nchar_required > nchar_available) then  ! Unexpected, but safeguard it
         nabs = 10**nchar_available - 1
         nchar_required = nchar_available
      end if
      write (iformat(3:3), '(i1)') nchar_required
      length = i1 + nchar_required - 1
      write (name(i1:length), iformat) nabs
   end if

   end subroutine numbered_name
