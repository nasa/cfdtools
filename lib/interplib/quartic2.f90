!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine quartic2 (nx, x, mode, x1, y1, yp1, x2, y2, yp2, ypp, y)
!
!  Description:   Evaluate the quartic defined by:
!
!                 mode = 1:  (x,y,y',y") at left end and (x,y,y')    at right
!                 mode = 2:  (x,y,y')    at left end and (x,y,y',y") at right
!
!                 The quartic is of this form, where h = x - x1:
!
!                 Y(X)  =  Y1 + H (C1 + H (C2 + H (C3 + H C4)))
!
!  HISTORY:  11/14/01  DAS  QUARTIC adapted from QUINTIC (y" at right end only).
!            06/28/13   "   QUARTIC2 adapted from QUARTIC to look like CUBINTRP.
!
!  AUTHOR:   David Saunders, ERC, Inc./NASA Ames, Moffett Field, CA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)  :: nx       ! # abscissas at which to evaluate quartic
   real,    intent (in)  :: x(nx)    ! Evaluation abscissas
   integer, intent (in)  :: mode     ! 1 => ypp applies to (x1, y1) end point;
                                     ! 2 => ypp applies to (x2, y2) end point
   real,    intent (in)  :: &
      x1, y1, yp1, x2, y2, yp2       ! Coordinates & 1st derivates at end pts.
   real,    intent (in)  :: ypp      ! 2nd derivative at end defined by mode
   real,    intent (out) :: y(nx)    ! Evaluations of the quartic at x(1:nx)

!  Local constants:

   integer, parameter :: ndim = 3    ! Declared dimension of matrix A

!  Local variables:

   integer :: i, ier
   real    :: A(ndim,ndim), coefs(4), h, h2, h3

!  Execution:

   h  = x2 - x1
   h2 = h*h
   h3 = h*h2

   select case (mode)

      case (1)  ! y" = ypp applies to the left end

         a(1,1) = h3
         a(2,1) = h2 * 3.
         a(1,2) = h2 * h2
         a(2,2) = h3 * 4.

         coefs(1) = yp1
         coefs(2) = ypp * 0.5
         coefs(3) = y2  - y1  - yp1*h - coefs(2)*h2
         coefs(4) = yp2 - yp1 - ypp*h

         call lusolve (2, ndim, A, coefs(3), ier)

      case (2)  ! y" = ypp applies to the right end

         a(1,1) = h2
         a(2,1) = h + h
         a(3,1) = 2.
         a(1,2) = h3
         a(2,2) = h2 * 3.
         a(3,2) = h  * 6.
         a(1,3) = h2 * h2
         a(2,3) = h3 * 4.
         a(3,3) = h2 * 12.

         coefs(1) = yp1
         coefs(2) = y2  - y1  - yp1*h
         coefs(3) = yp2 - yp1
         coefs(4) = ypp

         call lusolve (3, ndim, A, coefs(2), ier)

   end select

   do i = 1, nx
      h = x(i) - x1
      y(i) = y1 + h*(coefs(1) + h*(coefs(2) + h*(coefs(3) + h*coefs(4))))
   end do

   end subroutine quartic2
