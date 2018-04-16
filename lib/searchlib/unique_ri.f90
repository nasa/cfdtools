!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine unique_ri (nin, rlist, ilist, nout)
!
!  Filter out repeated abscissas from a sorted real list accompanied by an
!  integer list that indicates how any ordinate arrays should be matched with
!  the sorted abscissas.  See heap sort subroutine hsortri for performing the
!  preliminary sorting of the original list (but not eliminating duplicates).
!  Both input arrays are modified in place if any duplicates are present in
!  rlist(1:nin).  Abscissa comparisons employ a tolerance 10 * machine epsilon.
!
!  01/13/2010  D.A.Saunders  Initial implementation with 1-D quadrature in mind.
!
!  Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)    :: nin         ! Length of list upon input
   real,    intent (inout) :: rlist(nin)  ! Sorted abscissas: possible repeats
                                          ! upon input; none upon output
   integer, intent (inout) :: ilist(nin)  ! Associated ordinal list, likewise
   integer, intent (out)   :: nout        ! Length of filtered list upon output

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
      if (abs (rlist(i) - rlist(i-1)) > tol) then
         nkeep = nkeep + 1
         nlist(nkeep) = i
      end if
   end do

   do i = 2, nkeep
      ikeep = nlist(i)
      rlist(i) = rlist(ikeep)
      ilist(i) = ilist(ikeep)
   end do

   nout = nkeep

   deallocate (nlist)

   end subroutine unique_ri
