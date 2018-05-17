!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine handle_flatness (lunout, n, s, curvature, nfl, ifl, ifr)
!
!  Purpose:
!
!     In conjunction with companion utility DETECT_FLATNESS, we seek to blend
!  redistributed generatrix points further into regions of low/zero curvature
!  than is possible without artificial broadening of the curvature-based shape
!  function seen by the general-purpose routine ARBDIS used by CURVDIS.
!
!  The input curvature data points are modified in-place by adding certain
!  "sine bumps" available from earlier airfoil shape optimization work.
!
!     The shape functions are appropriately denormalized so they match the
!  height of the non-flat curvature region next to an end of each flat region.
!
!  History:
!
!     06/27/16  D.A.Saunders  Initial adaptation of VERTEX_CURVATURE, for
!                             CAPSULE_GRID purposes.
!     06/29/16    "     "     An Orion on a sting model has three flat segments
!                             such that the middle segment has all of its length
!                             broadened!  We need to limit such situations to no
!                             more than half the segment if the other end may be
!                             broadened as well.
!     05/11/18    "     "     Added lunout argument to allow a verbose mode in
!                             CURVDIS.
!
!  Author: David Saunders, AMA, Inc./NASA Ames Research Center, Moffett Fld., CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)    :: lunout        ! lunout > 0 produces diagnostics
   integer, intent (in)    :: n             ! # data pts. in s(:) & curvature(:)
   real,    intent (in)    :: s(n)          ! Arc length along the curve
   real,    intent (inout) :: curvature(n)  ! Curvature before and after the
                                            ! treatment of flat regions
   integer, intent (in)    :: nfl           ! Actual # flat regs. detected, >= 0
   integer, intent (in)    :: ifl(nfl)      ! Left and right index pairs that
   integer, intent (in)    :: ifr(nfl)      ! appear to bracket flat regions

!  Local constants:

   integer, parameter :: npmin = 5      ! Min. # pts. to broaden over
   integer, parameter :: npmax = 40     ! Max. ...
   real,    parameter :: exponent = 2.  ! Controls with/shape of blend fn.
   real,    parameter :: half = 0.5
   logical, parameter :: add = .false.  ! BEVAL won't be adding to earlier calls

!  Local variables:

   integer :: i, ie, il, ir, iregion, nb, nf
   real    :: blend(n), p(2), snorm(n), w  ! More than enough

!  Execution:

   p(1) = exponent

   do iregion = 1, nfl

!     Handle the left-hand end of this flat region:

      ie = ifl(iregion)           ! Left end point of current flat region
      if (ie >= npmin) then       ! Else too near the left end of the data
         ir = ifr(iregion)        ! Right end of flat region
         nf = ir - ie + 1         ! # pts. in the flat region
         nb = min (nf, npmax)     ! Limit the blending region
         if (iregion < nfl) then  ! Avoid possibly bumping into next broadening
            if (ir < n - npmin) nb = min (nb, nf/2)
         end if
         ir = ie + nb - 1         ! Adjusted end of flat region to treat
         if (ir  - ie >= npmin) then
            p(2) = curvature(ie-2) - curvature(ir)  ! Scaled height of shape fn.
            w = s(ir) - s(ie)                       ! Width of blend region
            do i = ie, ir
               snorm(i) = (s(i) - s(ie))/w
            end do
            call beval ('COSL', 2, p, add, ir-ie, snorm(ie), blend(ie))
            if (lunout > 0) then
               write (lunout, '(a, i2, a)') 'Flat region', iregion, ', left'
               write (lunout, '(i4, 2es16.8)') &
                  (i, curvature(i), blend(i), i = ie, ir-1)
            end if
            curvature(ie:ir-1) =  curvature(ie:ir-1) + blend(ie:ir-1)
            curvature(ie-1)    = (curvature(ie-2) + curvature(ie))*half
            if (lunout > 0) &
               write (*, '(i4, es16.8)') (i, curvature(i), i = ie-2, ir)
         end if
      end if

!     Likewise for the right-hand end:

      ie = ifr(iregion)      ! Right end point of current flat region
      if (n - ie < npmin) cycle

      il = ifl(iregion)      ! Left end of flat region
      nf = ie - il + 1       ! # pts. in the flat region
      nb = min (nf, npmax)   ! Limit the blending region
      if (iregion > 1) then  ! Avoid possibly bumping into next broadening
         if (il >= npmin) nb = min (nb, nf/2)
      end if
      il = ie - nb + 1       ! Adjusted end of flat region to treat
      if (ie  - il < npmin) cycle

      p(2) = curvature(ie+2) - curvature(il)  ! Scaled height of shape fn.
      w = s(ie) - s(il)                       ! Width of blend region
      do i = il, ie
         snorm(i) = (s(i) - s(il))/w
      end do
      call beval ('COSR', 2, p, add, ie-il, snorm(il+1), blend(il+1))
      if (lunout > 0) then
         write (lunout, '(a, i2, a)') 'Flat region', iregion, ', right'
         write (lunout, '(i4, 2es16.8)') &
            (i, curvature(i), blend(i), i = il+1, ie)
      end if
      curvature(il+1:ie) =  curvature(il+1:ie) + blend(il+1:ie)
      curvature(ie+1)    = (curvature(ie) + curvature(ie+2))*half
      if (lunout > 0) write (*, '(i4, es16.8)') (i, curvature(i), i = il, ie+2)
   end do

   if (lunout > 0) then
      write (lunout, '(a)') ' HANDLE_FLATNESS:  Arc length + adjusted curvature'
      write (lunout, '(i4, 2es16.8)') (i, s(i), curvature(i), i = 1, n)
   end if

   end subroutine handle_flatness
