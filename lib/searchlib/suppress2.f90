!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine suppress2 (nin, v1, v2, nout)
!
!  Filter out repeated abscissas from a sorted real variable array v1.
!  Associated ordinates v2(:) are kept in sync with v1(:).  The input arrays
!  are modified in place if any duplicates are present in v1(1:nin).
!  For 3 or more related arrays, UNIQUE_RI may be better (or a SUPPRESS3, say,
!  might be preferred.)
!
!  Abscissa comparisons employ a tolerance of 10*machine epsilon.
!
!  05/18/2018  D.A.Saunders  Companion to SUPPRESS1, q.v. Avoid the integer list
!                            of UNIQUE_RI that would normally contain 1:nin.
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
   real,    intent (inout) :: v2(nin)  ! Ordinates associated with v1(:)
   integer, intent (out)   :: nout     ! # array elements upon output, with any
                                       ! repeats found in v1(:) suppressed
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
      v2(i) = v2(ikeep)
   end do
   nout = nkeep

   deallocate (nlist)

   end subroutine suppress2
