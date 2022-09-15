!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine symevdsolve (n, V, D, tol, b)
!
!  Purpose:
!
!     Solve the symmetric linear system Ax = b given the eigenvalue decompostion
!     A = VDV' from subroutine symevd.  An ill-conditioned A is handled by sub-
!     stituting zero for the inverse of any eigenvalue with absolute value less
!     than the given tolerance.  This can be shown (thanks, Michael Saunders) to
!     produce a residual vector b - Ax with norm not exceeding that of RHS b,
!     meaning that the solution from the "pseudo" inverse is guaranteed to be
!     bounded.  Furthermore, this pseudoinverse Ak of rank k < n can be shown to
!     minimize || A - Ak || among rank k approximations.
!
!     The solution x overwrites input right-hand-side vector b.  Arguments D and
!     V are not changed.
!
!  History:
!
!     01/17/2022  D.A.Saunders  Companion to symmetric eigenvalue decomposition
!                               routine symevd, q.v.
!
!  Author:  David Saunders, AMA, Inc. at NASA Ames Research Center, CA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)    :: n      ! Dimension of the system being solved
   real,    intent (in)    :: V(n,n) ! Matrix of orthonormal eigenvectors of A
                                     ! as returned in argument A of symevd
   real,    intent (in)    :: D(n)   ! Eigenvalues in ascending order as output
                                     ! by symevd
   real,    intent (in)    :: tol    ! Eigenvalues with magnitude less than this
                                     ! value (say 1.e-5) are not inverted here;
                                     ! zeros are substituted for the inverses
   real,    intent (inout) :: b(n)   ! RHS vector on input; solution x on output

!  Local variables:

   integer :: i
   real    :: Dinverse(n)

!  Execution:

!  Form the adjusted inverse of D that stabilizes solution of Ax = VDV'x = b.
!  The inverse of (orthogonal) V is V'.

   do i = 1, n
      if (abs (D(i)) > tol) then
         Dinverse(i) = 1.0/D(i)
      else
         Dinverse(i) = 0.0
      end if
   end do

!  For Ax = VDV'x = b, the soln. is x = VEV'b  where E = [modified?] D inverse.

   b = matmul (V, Dinverse * matmul (transpose (V), b))

   end subroutine symevdsolve
