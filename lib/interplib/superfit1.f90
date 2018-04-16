!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine superfit1 (n, x, y, p, a, b, yc, rsq)

!     SUPERFIT1 performs the linear least squares portion of fitting a super-
!  ellipse to data points where the major axis a is known, the center in the x
!  direction is assumed to be zero, and the unknowns are the minor axis b, the
!  center in the y direction yc, and the exponent p.  SUPERFIT2 computes the
!  best exponent p by performing a 1-dimensional nonlinear minimization of the
!  squared residual rsq from SUPERFIT1 via calls with the iterates for p.  If
!  fitting of a standard ellipse is desired, simply use SUPERFIT1 on its own
!  with p input as 2.
!
!     This capability was prompted by a need to represent the cross-section of a
!  recessed arc-jet test article as a smooth discretized curve.  The data points
!  (x,y) here for the initial application are an averaged section derived from
!  laser scan measurements by other numerical techniques.
!
!     The data points may cover either one or two quadrants, but they must
!  represent either the top half of the [super]ellipse or the bottom half, but
!  not both.  If they are on the bottom half, the output b will be negative.
!  |x(i)| <= a is assumed for all i = 1:n.
!
!     Use SUPEREVAL to discretize the fitted curve.
!
!  07/29/2008  D.A.Saunders  Initial implementation.
!
!  Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)  :: n        ! Number of data points

   real,    intent (in)  :: x(n), &  ! Data coordinates on either the upper or
                            y(n)     ! lower half (or quarter) of the expected
                                     ! [super]ellise, but not both

   real,    intent (in)  :: p        ! Exponent of the superellipse

   real,    intent (in)  :: a        ! Major axis, with center x = 0 assumed

   real,    intent (out) :: b        ! Minor axis calculated here; b < 0
                                     ! if the lower portion of the curve is
                                     ! represented by the data points

   real,    intent (out) :: yc       ! Center y coordinate calculated here

   real,    intent (out) :: rsq      ! Squared residual from the linear least
                                     ! squares fit performed here
!  Procedures:

   external  hdesol                  ! Linear least least squares via QR
                                     ! factorization (Householder transformns.)
!  Local constants:

   real, parameter   :: one = 1.

!  Local variables:

   integer           :: i
   real              :: oneoverp
   real, allocatable :: Ab(:,:)      ! Augmented matrix A | b  for Ax ~ b
   real, allocatable :: w(:)         ! Work-space for HDESOL

!  Execution:

   oneoverp = one / p;                     allocate (Ab(n,3), w(n))

   do i = 1, n
      Ab(i,1) = (one - (abs (x(i) / a)) ** p) ** oneoverp
      Ab(i,2) =  one
      Ab(i,3) =  y(i)
   end do

   call hdesol (n, n, 3, Ab, w, rsq);      deallocate (Ab)

   b = w(1);  yc = w(2);                   deallocate (w)

   end subroutine superfit1
