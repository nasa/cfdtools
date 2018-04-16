!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine int_to_char (number, string, length)
!
!  Convert the given integer to a left-justified character string and return
!  the active length of that string.  The input character variable is assumed
!  to be long enough for any default (4-byte) integer the calling program is
!  expected to encounter.  For instance, if a multiblock grid is certain to
!  contain less than 10,000 blocks, the input string for a block number need
!  not be any longer than 4 bytes.  NO ERROR HANDLING IS PERFORMED.  Rather,
!  the largest possible number (in magnitude, positive or negative) is written
!  if the string argument is not long enough.
!
!  Upon return, string(1:length) could be written with ('a') format as an
!  alternative to constructing a run-time format from just the active length.
!  The active string may also be needed for other purposes (see History).
!
!  09/28/12  D.A.S.  Initial implementation, for constructing "group" names
!                    from grid block numbers during HDF5 I/O.
!
!  Author:  David Saunders, ERC, Inc. at NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer,       intent (in)  :: number  ! Integer, possibly negative
   character (*), intent (out) :: string  ! Equivalent string, left justified
   integer,       intent (out) :: length  ! string(1:length) is significant;
                                          ! any trailing characters are blank
!  Local variables:

   integer       :: i1, nabs, nchar_available, nchar_required
   character (4) :: iformat = '(in)'

!  Execution:

   nabs   = abs (number)
   string = ' '

   if (nabs == 0) then  ! We have to avoid the logarithm, so handle it first
      string(1:1) = '0'
      length = 1
   else
      nchar_available = len (string)
      i1 = 1
      if (number < 0) then
         nchar_available = nchar_available - 1
         i1 = 2
         string(1:1) = '-'
      end if
      nchar_required = int (log10 (real (nabs))) + 1
      if (nchar_required > nchar_available) then  ! Unexpected, but safeguard it
         nabs = 10**nchar_available - 1
         nchar_required = nchar_available
      end if
      write (iformat(3:3), '(i1)') nchar_required
      length = i1 + nchar_required - 1
      write (string(i1:length), iformat) nabs
   end if

   end subroutine int_to_char
