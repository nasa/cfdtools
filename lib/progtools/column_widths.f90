!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine column_widths (format_string, ncol, width, ier)
!
!  Description:
!
!     Determine the apparent widths of the items intended to be printed in one
!     row by the given format string, such as (f5.1, 2f10.6, 3es14.6, 3f9.3),
!     where the parentheses are expected.  The purpose is to enable right-justi-
!     fication (or centering?) of column headers in a printed table.
!
!     Handling of Fortran's option to indicate repeats appears awkward enough
!     that it is NOT attempted here.  Thus, something like this is NOT handled:
!                       (f5.1, 3(f10.6, 2es14.6))
!
!    Argument width(:) is a pointer so that it can be allocated here after the
!    number of columns has been determined from the format string.
!
!  History:
!
!     06/08/2021  D.A.Saunders  Initial design and implementation, prompted by
!                               the BODY_POINT_DATA utility.
!     06/09/2021    "      "    Introduced new "digit" function after finding
!                               function "number" unable to isolate multipliers
!                               because it can't reject E or F.
!
!  Author:  David Saunders, AMA, Inc. at NASA Ames Research Center, CA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   character*(*), intent (in)  :: format_string  ! Fortran format string,
                                                 ! including leading and
                                                 ! trailing parentheses
   integer,          intent (out) :: ncol        ! Apparent # columns
   integer, pointer, intent (out) :: width(:)    ! Apparent column widths,
                                                 ! allocated here
   integer,          intent (out) :: ier         ! 0 => no error detected;
                                                 ! 1 => flawed format?
!  Local constants:

   integer, parameter :: mxcol = 1000  ! Workspace overkill, but at no cost
   integer, parameter :: mxwidth = 10  ! Can any token exceed this # chars.?
!! character (2), parameter :: seps = ' ,'  ! Preferable, can't use token2

!  Local variables:

   integer :: i, idot, ino, inumber, ios, j, k, l, m, n, nx, ntoken
   integer, allocatable :: iwidth(:)  ! Enables setting widths before
                                      ! actual # columns is known
   real    :: rno, rnumber
   character (mxwidth), allocatable :: token(:)

!  Procedures:

   logical, external :: digit  ! T if the character is in the range 0-9
   external :: decode_number   ! Character string --> integer|real
   external :: tokens  ! Isolates tokens as a list, and up-cases them

!  Execution:

   ier = 1  ! Until done with no error
   l = len_trim (format_string)
   if (format_string(1:1) /= '(' .or. &
       format_string(l:l) /= ')') then
       write (*, '(2a)') 'Format string must include (): ', format_string(1:l)
       go to 99
   end if

   ntoken = 1
   do i = 2, l - 1  ! Count 1 or more tokens in the format string
      if (format_string(i:i) == ',') ntoken = ntoken + 1
   end do
   allocate (token(ntoken), iwidth(mxwidth))

!! call token2 (format_string(2:l-1), seps, ntoken, token)  ! Doesn't up-case
   call tokens (format_string(2:l-1), ntoken, token)        ! Up-cases

!  Parse each token for a possible multiplier at the beginning:

   nx   = 0  ! Awkwardness caused by mX in the format, affecting next token
   ncol = 0
   do j = 1, ntoken
      k = len_trim (token(j))
      do m = 1, k ! Until leading string is not an integer
         n = m
         if (.not. digit (token(j)(m:m))) exit
      end do
      if (n > 1) then  ! Apparent multiplier
         call decode_number (token(j)(1:n-1), inumber, rnumber, ios)
         if (ios /= 2) then  ! 0 => not a number; 1 => real number (only)
            write (*, '(2a)') 'Expected integer not found in ', token(j)(1:n-1)
            go to 99
         end if
      else
         inumber = 1
      end if
      ncol = ncol + inumber

      select case (token(j)(n:n))

         case ('A', 'I', 'L')  ! Not a real
            call decode_number (token(j)(n+1:k), ino, rno, ios)
            if (ios /= 2) then  ! 0 => not a number; 1 => real number (only)
               write (*, '(2a)') 'Expected integer missing in ', token(j)(n+1:k)
               go to 99
            end if
            iwidth(ncol-inumber+1:ncol) = ino
            iwidth(ncol-inumber+1) = iwidth(ncol-inumber+1) + nx
            nx = 0

         case ('E', 'F')  ! A real
            if (token(j)(n+1:n+1) == 'S') n = n + 1
            idot = index (token(j)(n+1:k), '.')
            call decode_number (token(j)(n+1:n+idot-1), ino, rno, ios)
            if (ios /= 2) then  ! 0 => not a number; 1 => real number (only)
               write (*, '(2a)') 'Expected integer absent in ', &
                                 token(j)(n+1:n+idot-1)
               go to 99
            end if
            iwidth(ncol-inumber+1:ncol) = ino
            iwidth(ncol-inumber+1) = iwidth(ncol-inumber+1) + nx
            nx = 0

         case ('X')
            ncol = ncol - inumber  ! Not an actual column
            nx   = inumber         ! Next column is this much wider

         case default
            write (*, '(2a)') 'Unhandled formatting: ', token(j)(n:n), &
               'Abandoning column justification.'
            go to 99
      end select

   end do  ! Next token

   allocate (width(ncol))
   width(:) = iwidth(1:ncol)
   deallocate (iwidth)
   ier = 0

99 return

   end subroutine column_widths
