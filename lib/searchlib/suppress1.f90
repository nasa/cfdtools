!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine suppress1 (nin, v1, nout)
!
!  Filter out repeated abscissas from a sorted real variable array v1.
!  The input array is modified in place if any duplicates are present in
!  v1(1:nin).  Abscissa comparisons employ a tolerance of 10*machine epsilon.
!
!  05/18/2018  D.A.Saunders  Adaptation of UNIQUE_RI utility that avoids the
!                            integer list that would normally contain 1:nin.
!                            See also SUPPRESS2 for the common case of removing
!                            duplicates from a pair of related arrays.
!                            Prompted by extensions to NEQAIR_Integration.
!
!  Author:  David Saunders, AMA, Inc./NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)    :: nin      ! # abscissas input with possible repeats
   real,    intent (inout) :: v1(nin)  ! Sorted abscissas: possible repeats
                                       ! upon input; none upon output
   integer, intent (out)   :: nout     ! # abscissas upon output, with any
                                       ! repeats suppressed
!  Local variables:

   integer :: i, ikeep, nkeep
   real    :: tol
   integer, allocatable :: nlist(:)       ! New list of elements to keep

!  Execution:

   tol = 10. * epsilon (tol)

!  Avoid multiple moves by performing two passes:

   allocate (nlist(nin))

   nkeep = 1
   do i = 2, nin
      if (abs (v1(i) - v1(i-1)) > tol) then
         nkeep = nkeep + 1
         nlist(nkeep) = i
      end if
   end do

   do i = 2, nkeep
      ikeep = nlist(i)
      v1(i) = v1(ikeep)
   end do
   nout = nkeep

   deallocate (nlist)

   end subroutine suppress1
