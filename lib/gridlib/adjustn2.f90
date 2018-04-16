!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine adjustn2 (ndata, i1, i2, x1, y1, z1, ia, ib, snorm, &
                        x2, y2, z2, method)

!  ADJUSTN2 is a variant of ADJUSTN (itself a version of CHANGEN to allow for
!  data points outside the range of interest) intended for changing either the
!  number of points on a 3-space line or the relative spacing, or both.  The
!  desired relative spacing snorm(1:ib-ia+1) is an input as opposed to being
!  derived from the given curve in ADJUSTN.
!
!  03/07/13  DAS  Adaptation of ADJUSTN for CAPSULE_GRID purposes.
!
!  Author:  David Saunders, ERC, Inc./NASA Ames Research Ctr., Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer,   intent (in)  :: ndata      ! Number of data points, 1:ndata
   integer,   intent (in)  :: i1, i2     ! Pts. i1:i2 are redistributed as ia:ib
   real,      intent (in)  :: x1(ndata), y1(ndata), z1(ndata)  ! Data points
   integer,   intent (in)  :: ia, ib     ! First & last redistributed points
   real,      intent (in)  :: snorm(ib-ia+1)  ! Desired relative spacing
   real,      intent (out) :: x2(ib), y2(ib), z2(ib)   ! Redistributed coords.
   character, intent (in)  :: method * 1 ! Type of fit to be used by PLSCRV3D:
                                         ! M(onotonic) | B ("loose") | L(inear)
!  Procedures:

   external :: chords3d,   &             ! Arc-length utility
               plscrv3d                  ! 3-space local cubic spline utility

!  Local constants:

   real,    parameter :: derivs = -999.  ! Suppresses them
   logical, parameter :: norm   = .false.! No need to normalize
   logical, parameter :: closed = .false.

!  Local variables:

   integer :: i, ieval
   real    :: t(ndata), teval, total
   logical :: new

!  Execution:

   call chords3d (ndata, x1, y1, z1, norm, total, t)

   total  = t(i2) - t(i1)
   ieval  = i1
   new    = .true.

   x2(ia) = x1(i1)
   y2(ia) = y1(i1)
   z2(ia) = z1(i1)

   do i = ia + 1, ib - 1

      teval = t(i1) + snorm(i-ia+1)*total

      call plscrv3d (ndata, x1, y1, z1, t, method, new, closed, &
                     teval, ieval, x2(i), y2(i), z2(i), derivs)
      new = .false.

   end do

   x2(ib) = x1(i2)
   y2(ib) = y1(i2)
   z2(ib) = z1(i2)

   end subroutine adjustn2
