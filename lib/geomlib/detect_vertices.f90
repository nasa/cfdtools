!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine detect_vertices (n, curvature, nvmax, iv, nv)
!
!  Purpose:
!
!     For the given curvature distribution, which can be assumed to apply to a
!  geometry that has been normalized to avoid scaling effects, identify apparent
!  sharp corners (vertices) by looking for one-point spikes typical of what
!  finite-difference derivatives produce.  Heuristics are unavoidable.
!
!  History:
!
!     11/01/13  D.A.Saunders  Initial implementation, for CAPSULE_GRID.
!
!  Author: David Saunders, ERC, Inc./NASA Ames Research Center, Moffett Fld., CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)    :: n             ! # data pts. in curvature(:)
   real,    intent (inout) :: curvature(n)  ! Curvature data
   integer, intent (in)    :: nvmax         ! Maximum # vertices allowed for
   integer, intent (out)   :: iv(nvmax)     ! Indices apparently <-> vertices
   integer, intent (out)   :: nv            ! Actual # vertices detected, >= 0

!  Local constants:

   real, parameter ::    &
      background = 20.,  &  ! Curvature below which a spike neighbor must be
      spike      = 200.     ! Curvature above which a spike must be

!  Local variables:

   integer :: i

!  Execution:

   nv = 0

   do i = 2, n - 1
      if (abs (curvature(i)) > spike) then
         if (abs (curvature(i-1)) < background) then
            if (abs (curvature(i+1)) < background) then
               nv = nv + 1
               iv(nv) = i
               if (nv == nvmax) exit
            end if
         end if
      end if
   end do

   end subroutine detect_vertices
