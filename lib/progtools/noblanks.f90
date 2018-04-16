!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine noblanks (string)

!     This is a specialized string utility prompted by the need to convert
!     column headers of an inconvenient form into single-token column headers:
!
!        Input string:    pw (Pa)      qw (W/m^2)     Tw (K)  ...
!                  or:    pw [Pa]      qw [W/m^2]     Tw [K]  ...
!
!        Output string:   pw,Pa        qw,W/m^2       Tw,K
!
!     This allows column headers to be cut and pasted along with column values,
!     because the numbers of tokens now match.
!
!     Note that unitless quantities pose no implementation problem: they are
!     simply left in-place.
!
!     No error checking is done; use on sensible data only.
!
!  Assumptions:
!
!     o  The string contains (mostly) variable name/units pairs, with one
!        space between the elements of each pair.
!     o  There is no need to compact the output tokens by separating them
!        with just a single space.  Therefore, each output token simply starts
!        in the same string position as the associated input item.
!
!  History:
!
!     23 Aug. 2014  D.A.Saunders  Adaptation of subroutine noquotes, for a
!                                 set of procedures involving spreadsheets.
!                                 Variable names from program BLAYER containing
!                                 embedded blanks prompted this functionality.
!
!  Author:  David Saunders, ERC, Inc./NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   character, intent (inout)  :: string * (*)  ! String to be adjusted in-place

!  Local variables:

   integer       :: i1, i2, l
   character (1) :: right

!  Execution:

   l  = len_trim (string)
   i2 = 1

   do  ! Until i2 = l
      i1 = i2 + 1
      if (i1 > l) exit

      if (string (i1:i1) == '(') then
         right = ')'
      else if (string (i1:i1) == '[') then
         right = ']'
      else
         i2 = i1
         cycle
      end if

      string (i1-1:i1-1) = ','
      i2 = i1 + 1

      do  ! Until matching paren or bracket is found
         if (string(i2:i2) == right) exit
         if (i2 > l) exit
         i2 = i2 + 1
      end do

      string(i1:i2-2) = string(i1+1:i2-1)
      string(i2-1:i2) = ' '
   end do

   end subroutine noblanks
