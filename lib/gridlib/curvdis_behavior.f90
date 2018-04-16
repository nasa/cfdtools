!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   program curvdis_behavior

!  Evaluate the curvature-based shape function constructed by CURVDIS at each
!  data point for use by ARBDIS to generate 1-D grid point distributions.  In
!  particular, show how it varies with the exponent or power p for a given value
!  of curvature, k.  The shape function, for p in [0. 1] and |k| >= 0, is:
!
!                                                1
!                         shape (p, k)  =  ------------
!                                          (1 + |k|)**p
!
!  Note that curvature is extremely dependent on the data scaling.  It is highly
!  recommended that CURVDIS2 be used, not CURVDIS directly.  CURVDIS2 normalizes
!  the data to the unit square prior to the grid point generation, and handles
!  non-convergence by reducing the power by 0.1 and trying again.  Normalized
!  curvatures should not be as extreme as those associated with subscale geom-
!  etries such as capsules tested in shock tunnels.  Normalized radii may be as
!  small as, say, 0.01, meaning normalized curvatures may be as large as 100 in
!  magnitude, but not in the thousands.
!
!  10/24/13  D.A.Saunders  Belated implementation prompted by efforts to improve
!                          the blending between straight-line segments with zero
!                          curvature and circular segments typical of reentry
!                          capsules.
!
!  Author:  David Saunders, ERC, Inc. at NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Constants:

   integer, parameter :: &
      luncrt = 6, &
      lunkbd = 5, &
      lunout = 1
   real, parameter :: &
      dp  = 0.05, &          ! Discretization interval of p for a given k
      one = 1.

!  Variables:

   integer :: i, neval
   real    :: k, p

   real, allocatable, dimension (:) :: power, shape

!  Execution:

   open (lunout, file='shape_function.dat', status='unknown')

   write (lunout, '(a)') &
      'TITLE = "Curvature-based shape function"', &
      'VARIABLES =', &
      '"Power"', &
      '"Shape"'

   neval = 1 + nint (one  /dp)
   allocate (power(neval), shape(neval))

!  Prompt for a curvature value, evaluate the shape function on [0, 1],
!  and save the values as a zone for Tecplot plotting.

   do  ! Until ^D (EOF)

      write (luncrt, '(a)', advance='no') &
         'Normalized curvature value, or ^D to quit: '
      read  (lunkbd, *, end=90) k
      k = abs (k)

      do i = 1, neval
         p = dp * real (i - 1);  if (i == neval) p = one  ! Be sure of it
         power(i) = p
         shape(i) = one / (one + k)**p
      end do

      write (lunout, '(a, f6.1, a)') 'ZONE T="k =', k, '"'
!!    write (lunout, '(a, i3, a)') 'I=', neval, ', J=1, K=1, ZONETYPE=Ordered'
      write (lunout, '(f6.3, es16.8)') (power(i), shape(i), i = 1, neval)

   end do

90 close (lunout)
   deallocate (power, shape)

   end program curvdis_behavior
