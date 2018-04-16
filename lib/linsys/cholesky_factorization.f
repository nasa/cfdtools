c-------------------------------------------------------------------------------

      subroutine cholesky_factorization (n, a, info)

c     This is a Fortran 90 update of the Linpack routine spofa (a,lda,n,info).
c     It replaces BLAS utilities with language intrinsics, and assumes the
c     matrix is stored in a(n,n) instead of the more general a(lda,n).

      implicit none

c     Arguments:

      integer, intent (in)    :: n
      real,    intent (inout) :: a(n,n)
      integer, intent (out)   :: info

c     The following documentation from spofa still applies:
c
c     spofa factors a real symmetric positive definite matrix.
c
c     spofa is usually called by spoco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for spoco) = (1 + 18/n)*(time for spofa) .
c
c     on entry
c
c        n       integer
c                the order of the matrix  a .
c
c        a       real(lda, n)
c                the symmetric matrix to be factored.  only the
c                diagonal and upper triangle are used.
c
c     on return
c
c        a       an upper triangular matrix  r  so that  a = trans(r)*r
c                where  trans(r)  is the transpose.
c                the strict lower triangle is unaltered.
c                if  info .ne. 0 , the factorization is not complete.
c
c        info    integer
c                = 0  for normal return.
c                = k  signals an error condition.  the leading minor
c                     of order  k  is not positive definite.
c
c     linpack version dated 08/14/78.
c     Cleve Moler, University of New Mexico, Argonne National Lab.
c
c     05/14/07   David Saunders  Fortran 90 version better suited to a(n,n).
c
c-------------------------------------------------------------------------------

c     Local constants:

      real, parameter :: zero = 0.0e0

c     Local variables:

      real    :: s, t
      integer :: j, k

c     Execution:

c     Avoid the j > 1 test inside the outer loop:

      if (a(1,1) <= zero) then

         info = 1

      else

         info = 0
         a(1,1) = sqrt (a(1,1))

         do j = 2, n

            t = a(1,j) / a(1,1)  ! Avoid the dot product when k = 1
            a(1,j) = t
            s = t*t

            do k = 2, j - 1
               t = (a(k,j) - dot_product (a(1:k-1,k), a(1:k-1,j))) /
     >              a(k,k)
               a(k,j) = t
               s = s + t*t
            end do

            s = a(j,j) - s

            if (s <= zero) then
               info = j
               exit
            end if

            a(j,j) = sqrt (s)

         end do

      end if

      end subroutine cholesky_factorization
