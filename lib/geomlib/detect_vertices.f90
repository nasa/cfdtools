!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine detect_vertices (lunout, n, curvature, nvmax, iv, nv)
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
!     05/09/18     "     "    Background = 20 and spike = 200 were missing
!                             (normalized) vertices in the 100-150 curvature
!                             range for a Mars Sample Return Lander geometries.
!                             Also, change the inner two tests from and to or.
!     05/10/18     "     "    Lower spike to 50. so as not to miss a shoulder
!                             cusp.
!
!  Author: David Saunders, ERC, Inc./NASA Ames Research Center, Moffett Fld., CA
!                 Now with AMA, Inc. at ARC.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)    :: lunout        ! lunout > 0 produces diagnostics
   integer, intent (in)    :: n             ! # data pts. in curvature(:)
   real,    intent (inout) :: curvature(n)  ! Curvature data
   integer, intent (in)    :: nvmax         ! Maximum # vertices allowed for
   integer, intent (out)   :: iv(nvmax)     ! Indices apparently <-> vertices
   integer, intent (out)   :: nv            ! Actual # vertices detected, >= 0

!  Local constants:

   real, parameter ::    &
      background = 15.,  &  ! Curvature below which a spike neighbor must be
      spike      = 50.      ! Curvature above which a spike must be

!  Local variables:

   integer :: i

!  Execution:

   nv = 0

   do i = 2, n - 1
      if (abs (curvature(i)) > spike) then
         if (abs (curvature(i-1)) < background .or. &
             abs (curvature(i+1)) < background) then
             nv = nv + 1
             iv(nv) = i
             if (lunout > 0) write (lunout, '(a, i5)') ' Vertex found at i =', i
             if (nv == nvmax) exit
         end if
      end if
   end do

   end subroutine detect_vertices
