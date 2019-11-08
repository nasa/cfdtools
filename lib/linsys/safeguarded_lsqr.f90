!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine safeguarded_lsqr (m, n, A, b, lambda, x, ssqmin, ier)
!
!  Description:
!
!     Solve the square or rectangular system A x = b for the case where A is
!     (close to) singular by augmenting the system with additional equations
!     lambda I ~ 0 that have the effect of ensuring full rank while promoting a
!     solution x with norm as small as possible, for any lambda > 0.  This is
!     equivalent to working with the form (A'A + lambda^2 I) x = A'b, but is
!     solved here more simply (and accurately) with the QR factorization of
!     A augmented by lambda I.  A good choice for lambda should be sqrt (machine
!     epsilon).
!
!  Programmer note:
!
!     This capability was prompted by singularities observed in the boundary
!     layer regions of computational fluid dynamics volume grids.  The number
!     of variables n is not expected to be large, so dense-matrix techniques
!     are employed.  Local use of a copy of A with room for the augmentation is
!     considered preferable to having the calling program provide the matrix
!     with room for the augmentation.
!
!  History:
!
!     10/15/19  D.A.S.  Introduced as a work-around for occasional singularities
!                       encountered with the structured volume grid form of the
!                       ADT (alternating digital tree) search method, after it
!                       was found that an alternative KDTREE-based hybrid method
!                       in FLOW_INTERP, applied to cell centroids, can choose
!                       seriously wrong cells near shock boundaries.  Michael
!                       Saunders advised this Levenberg-Marquardt-type approach.
!
!  Author:  David Saunders, AMA, Inc. at NASA Ames Research Center, CA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)  :: m, n    ! Row/column dimensions of A, which may be
                                    ! square (m = n) or have m < n or m > n
   real,    intent (in)  :: A(m,n)  ! Matrix found to be (close to) singular
   real,    intent (in)  :: b(m)    ! RHS vector; lost if b(:) is reused as x(:)
   real,    intent (in)  :: lambda  ! Small positive multiplier applied to I
   real,    intent (out) :: x(n)    ! Desired least-squares solution; passing
                                    ! b(:) as x(:) is OK
   real,    intent (out) :: ssqmin  ! Residual || b - Ax ||**2
   integer, intent (out) :: ier     ! 0 means no error found;
                                    ! > 0 means an allocation error
!  Local variables:

   integer :: i
   real, allocatable :: Aplus(:,:), w(:)

!  Execution:

!  The venerable HDESOL implementation augments the intended overdetermined
!  system with the RHS column for efficiency reasons.

   allocate (Aplus(m+n,n+1), w(m+n), stat=ier)
   if (ier /= 0) go to 99

   Aplus(1:m,1:n) = A(:,:)
   Aplus(1:m,n+1) = b(:)
   do i = 1, n
      Aplus(m+i,i) = lambda
      Aplus(m+i,n+1) = 0.
   end do

   call hdesol (m+n, m+n, n+1, Aplus, w, ssqmin)

   x(:) = w(1:n)

   deallocate (Aplus, w)

99 continue

   end subroutine safeguarded_lsqr
