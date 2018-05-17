!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine detect_flatness (lunout, n, curvature, nflmax, ifl, ifr, nfl)
!
!  Purpose:
!
!     For the given curvature distribution, which can be assumed to apply to a
!  geometry that has been normalized to avoid scaling effects, identify apparent
!  flat regions as part of blending them better into neighboring regions of
!  high curvature by broadening the non-flat region artificially.  The goal is
!  to put more points on the conical regions of a typical capsule generatrix
!  near the outer shoulder (in particular) as has been found desirable to
!  improve resolution of the surface heat flux spikes that are likely to be
!  present in real gas flow solutions.
!
!     The logic below could probably be more elegant, but it appears to be
!  robust for the cases of Stardust, Sample Return Capsule, and a scale model of
!  Orion on a sting.
!
!  History:
!
!     06/27/16  D.A.Saunders  Initial adaptation of the earlier detect_vertices,
!                             for CAPSULE_GRID.
!     06/29/16    "     "     Raised flatmax from 0.1 to 0.5 to allow the same
!                             idea to work on the forebody of Orion-type cases.
!     05/10/18    "     "     See the comment for flatmax below.
!
!  Author: David Saunders, AMA, Inc./NASA Ames Research Center, Moffett Fld., CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)    :: lunout        ! lunout > 0 produces diagnostics
   integer, intent (in)    :: n             ! # data pts. in curvature(:)
   real,    intent (inout) :: curvature(n)  ! Curvature data
   integer, intent (in)    :: nflmax        ! Max. # flat regions  allowed for
   integer, intent (out)   :: ifl(nflmax)   ! Left and right index pairs that
   integer, intent (out)   :: ifr(nflmax)   ! appear to bracket flat regions
   integer, intent (out)   :: nfl           ! Actual # flat regs. detected, >= 0

!  Local constants:

   real, parameter ::    &
      flatmax = 0.9         ! Curvature value below which we say "flat";
                            ! normalized Orion spherical section has 0.4166667;
                            ! A Mars capsule had a spherical section varying
                            ! between 0.493 and 0.502, causing consternation
                            ! when flatmax was 0.5

!  Local variables:

   integer :: i, j

!  Execution:

   if (lunout > 0) write (lunout, '(a, i5)') ' DETECT_FLATNESS # pts. n:', n
   nfl = 0
   i   = 1

   do  ! Until i > n
      if (lunout > 0) write (lunout, '(i4, es18.6)') i, curvature(i)
      if (abs (curvature(i)) < flatmax) then  ! Determine the extent of flatness
         nfl = nfl + 1
         if (lunout > 0) write (lunout, '(a, i4)') '   nfl:', nfl
         ifl(nfl) = i
         ifr(nfl) = i
         do j = i+1, n
            if (lunout > 0) write (lunout, '(i7, es18.6)') j, curvature(j)
            if (abs (curvature(j)) > flatmax) then
               ifr(nfl) = j - 1
               i = j + 1
               exit
            else if (j == n) then
               ifr(nfl) = n
               i = n + 1
               exit
            end if
            i = j
         end do
      else
         i = i + 1
      end if
      if (i > n) exit
      if (nfl == nflmax) exit
   end do

   if (lunout > 0) then
      write (lunout, '(a)') '#  i   ifl   ifr'
      write (lunout, '(i4, 2i6)') (i, ifl(i), ifr(i), i = 1, nfl)
   end if

   end subroutine detect_flatness
