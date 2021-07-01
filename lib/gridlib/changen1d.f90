!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine changen1d (i1, i2, s1, ia, ib, s2)
!
!  Description:
!
!     CHANGEN1D redistributes abscissas s1(i1:i2) to a different number
!     of abscissas s2(ia:ib) with the same end points and the same relative
!     spacing, meaning if the number of spaces is doubled, the spaces are all
!     halved (say). The abscissas might be "x" or "t" (times in a time history),
!     or normalized arc lengths on [0, 1].
!
!     The formulation comes from the linear transformation that converts
!     x in [a,b] to y in [p,q], namely y = ((q-p)/(b-a))x + (bp-aq)/(b-a).
!     Here, we want to transform j in [ia,ib] to i in [i1,i2], giving a
!     fractional index or a data interval in which to interpolate linearly for
!     each output index j (not as obvious as one would think).  See History.
!
!  History:
!
!     12/22/97  DAS  3-space grid line utility CHANGEN wrapped around PLSCRV3D.
!     10/23/02   "   2-space variant: CHANGEN2D.
!     06/23/14   "   1-space variant: CHANGEN1D, for a time history application.
!                    1:m --> 1:n is the likely usage, but need not be.
!     06/23/21   "   Added the refinement from densify_grid_block that fixes
!                    the weakness of the algebraic method.  (E.g., if we double
!                    the number of cells, it amounts to inserting cell mid-pts.,
!                    which means cell growth rates become stair-stepped.)
!
!  Author:  David Saunders, ERC, Inc. at NASA Ames Research Center, CA
!           Now with AMA, Inc. at NASA ARC.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)  :: i1, i2  ! First & last data points
   real,    intent (in)  :: s1(i2)  ! Data points in i1:i2
   integer, intent (in)  :: ia, ib  ! First & last redistributed points
   real,    intent (out) :: s2(ib)  ! Redistributed coordinates in ia:ib,
                                    ! with s2(ia) = s1(i1) and s2(ib) = s1(i2)

!  Local constants:

   integer, parameter :: lunlog = 6
   real,    parameter :: half = 0.5, one = 1., zero = 0.

!  Local variables:

   integer :: i, ier, ip, mi, ni
   real    :: p, r, ri1, rip
   real, allocatable, dimension (:) :: smid_in, spacing_in, work

!  Execution:

   r   = real (i2 - i1) / real (ib - ia)
   ri1 = real (i1)

   do i = ia, ib - 1
      rip   = ri1 + r * real (i - ia)
      ip    = int (rip)
      p     = rip - real (ip)
      s2(i) = (one - p) * s1(ip) + p * s1(ip+1)
   end do

   s2(ib) = s1(i2)

!  The above does not guarantee preservation of the smoothly
!  varying cell growth rates/cell spacings that are presumed
!  to exist in the input grid.  Use the input cell spacings as
!  the shape function ARBDIS employs to impose a desired type
!  of point distribution (as in CURVDIS).  Moreover, use the
!  above algebraic result as a starting guess:

   allocate (smid_in(i2+1), spacing_in(i2+1), work(5*ib))

   do i = i1 + 1, i2
      smid_in(i)    = (s1(i) + s1(i-1))*half
      spacing_in(i) =  s1(i) - s1(i-1)
   end do

!  Avoid extrapolation off the ends of the shape function:

   smid_in(i1) = zero;  smid_in(i2+1) = s1(i2)
   call linear_interp (smid_in(i1+1), spacing_in(i1+1), &
                       smid_in(i1+2), spacing_in(i1+2), &
                       smid_in(i1),   spacing_in(i1))
   call linear_interp (smid_in(i2-1), spacing_in(i2-1), &
                       smid_in(i2),   spacing_in(i2),   &
                       smid_in(i2+1), spacing_in(i2+1))

   mi = ib - ia + 1
   ni = i2 - i1 + 1

   call arbdis (mi, zero, s1(i2), ni+1, smid_in(i1), &
                spacing_in(i1), 'I', -lunlog, work, s2, ier)
   if (ier /= 0) then
      write (*, '(a, i2, /, a, 3i5)') &
         ' *** CHANGEN1D: Bad ier from ARBDIS:', ier
   end if

   deallocate (smid_in, spacing_in, work)

   return

!  Internal procedure:

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine linear_interp (x1, y1, x2, y2, xinterp, yinterp)  ! Obvious ...
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      real, intent (in)  :: x1, y1, x2, y2, xinterp
      real, intent (out) :: yinterp

!     Local variables:

      real :: slope

!     Execution:

      slope   = (y2 - y1)/(x2 - x1)
      yinterp = y1 + slope*(xinterp - x1)

      end subroutine linear_interp

      END SUBROUTINE CHANGEN1D
