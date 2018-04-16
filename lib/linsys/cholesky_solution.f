c-------------------------------------------------------------------------------

      subroutine cholesky_solution (n, a, b)

c     This is a Fortran 90 update of the Linpack routine sposl (a,lda,n,b).
c     It replaces BLAS utilities with language instrinsics, and assumes the
c     matrix is stored in a(n,n) instead of the more general a(lda,n).

      implicit none

c     Arguments:

      integer, intent (in)    :: n
      real,    intent (in)    :: a(n,n)
      real,    intent (inout) :: b(n)

c     The following documentation from sposl still applies:

c     sposl solves the real symmetric positive definite system
c     a * x = b
c     using the factors computed by spoco or spofa.
c
c     on entry
c
c        a       real(n, n)
c                the output from cholesky_factorization
c
c        n       integer
c                the order of the matrix  a .
c
c        b       real(n)
c                the right hand side vector.
c
c     on return
c
c        b       the solution vector  x .
c
c
c     linpack version dated 08/14/78.
c     Cleve Moler, University of New Mexico, Argonne National Lab.
c
c     05/14/07     David Saunders    Fortran 90 version better suited to a(n,n).
c
c-------------------------------------------------------------------------------

c     Local variables:

      integer :: k

c     Execution:

c     Solve trans(r)*y = b

      b(1) = b(1) / a(1,1)

      do k = 2, n
         b(k) = (b(k) - dot_product (a(1:k-1,k), b(1:k-1))) / a(k,k)
      end do

c     Solve r*x = y

      do k = n, 1, -1
         b(k) = b(k) / a(k,k)
cccc     t = -b(k)
cccc     call saxpy(k-1,t,a(1,k),1,b(1),1)
         b(1:k-1) = b(1:k-1) - b(k)*a(1:k-1,k)
      end do

      end subroutine cholesky_solution
