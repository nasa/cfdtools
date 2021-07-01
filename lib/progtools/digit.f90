!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   function digit (char)
!
!  Description:
!
!     if (digit (c)) then
!        single-character c is an integer from 0 through 9.
!
!  History:
!
!     06/09/2021  D.A.Saunders  Belated implementation when long-used functions
!                               number and alpha proved inadequate for parsing a
!                               Fortran format string.
!
!  Author:  David Saunders, AMA, Inc. at NASA Ames Research Center, CA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Argument:

   logical                    :: digit  ! Returned as T if char is one of 0-9
   character (1), intent (in) :: char   ! Character being tested

!  Constant:

   character (10), parameter :: digits = '0123456789'

!  Execution:

   digit = index (digits, char) /= 0

   end function digit
