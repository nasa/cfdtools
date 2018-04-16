!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine insertxy (ndim, n, xin, yin, x, y)
!
!  Insert real point (xin, yin) into arrays x(:) and y(:) where x(:) is ordered,
!  increasing or decreasing.  If xin is already present in x(:), the insertion
!  is suppressed.  The input value of n is assumed to be at least 2.
!
!  06/02/2014  D.A.Saunders  Initial implementation of a utility that should
!                            have been available long ago.
!
!  Author:  David Saunders, ERC, Inc./NASA Ames Research Center, Moffett Fld, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)    :: ndim      ! Space available in x(:) and y(:)
   integer, intent (inout) :: n         ! Input with space used so far; normally
                                        ! incremented by 1 (else not changed)
   real,    intent (in)    :: xin, yin  ! Point being inserted
   real,    intent (inout) :: x(ndim)   ! Abscissas, ordered from 1:n, either
                                        ! increasing or decreasing
   real,    intent (inout) :: y(ndim)   ! Corresponding ordinates

!  Constants:

   real, parameter :: one = 1.

!  Variables:

   integer :: left, nshift
   real    :: arrow

!  Procedures:

   external :: interval  ! Efficient 1-D search utility

!  Execution:

   arrow = sign (one, x(2) - x(1))
   left  = 1
   call interval (n, x, xin, arrow, left)

!  Pointer left cannot be n because of how interval works, so:

   if (xin*arrow >= x(n)*arrow) left = left + 1

   if (xin /= x(left)) then
      if (n < ndim) then
         nshift = n - left
         if (nshift > 0) then  ! Make room by shifting to the "right":
            x(left+2:left+1+nshift) = x(left+1:left+nshift)
            y(left+2:left+1+nshift) = y(left+1:left+nshift)
         end if
         x(left+1) = xin
         y(left+1) = yin
         n = n + 1
      end if
   end if

   end subroutine insertxy
