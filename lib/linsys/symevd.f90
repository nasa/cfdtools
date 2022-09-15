!-------------------------------------------------------------------------------
!
   subroutine symevd (n, A, D, ier)
!
!  Purpose:
!
!     Calculate the eigenvalue decomposition VDV' of symmetrix matrix A.
!     Originally, this repackaging of the LAPACK routine included the possible
!     adjustment that stabilizes solution of ill-conditioned systems.  But that
!     meant returning the [possibly modified] inverse of D, and the adjustment
!     for a given tolerance really belongs in the solve step.  Now, this level
!     above dsyev merely relieves the application of determining and allocating
!     the necessary work-space, and omits the option to enter the lower triangle
!     in A.
!
!  History:
!
!     c. 2005     Alexander Barth  Part of a key portion of an n-dimensional
!                                  optimal interpolation module.
!     01/17/2022  D.A.Saunders     Initial remodularization following futile
!                                  attempt to substitute the Cholesky factor-
!                                  ization for what turns out to be possibly
!                                  badly-conditioned systems.  [With smooth
!                                  analytic function data, the interpolations
!                                  can be extraordinarily accurate, but adding
!                                  some noise throws off the Cholesky approach,
!                                  which does not deal with ill-conditioning.]
!     01/28/2022     "     "       Passing rwork in the second call to dsyev
!                                  instead of work eluded detection in test
!                                  driving program nbyn4.  (See NBYN.)
!
!  Author:  David Saunders, AMA, Inc. at NASA Ames Research Center, CA.
!  
!-------------------------------------------------------------------------------

   implicit none

!  Arguments:

   integer, intent (in)    :: n         ! Order of sysmmetric matrix A; n >= 1
   real,    intent (inout) :: A(n,n)    ! Upper triangle of A on input;
                                        ! orthonormal eigenvector factor V on
                                        ! output, by columns
   real,    intent (out)   :: D(n)      ! The vector of eigenvalues in ascending
                                        ! order
   integer, intent (out)   :: ier       ! = 0 => successful decomposition
                                        ! > 0 => convergence failure

!  Local variables:

!! integer :: i
   integer :: lwork
   real    :: rwork, dummy
   real, allocatable, dimension (:) :: work

!  Execution:

!  Determine the necessary symmetric eigenvalue decomposition work-space.

   call dsyev ('V','U', n, dummy, n, dummy, rwork, -1, ier)

   lwork = ceiling (rwork)  ! Smallest integer >= rwork

   allocate (work(lwork))

!  Calculate the eigenvalues and eigenvectors:

   call dsyev ('V','U', n, A, n, D, work, lwork, ier)

!! write (6, '(a)') ' Eigenvalues:'
!! write (6, '(1x, 6es24.14)') D(:)
!! write (6, '(a)') ' Eigenvectors:'
!! do i = 1, n
!!    write (6, '(1x, 6es24.14)') A(i,:)
!! end do

   deallocate (work)

   end subroutine symevd
