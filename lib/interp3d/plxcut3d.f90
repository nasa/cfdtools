!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine plxcut3d (n, x, y, z, t, ix, calct, xcut, tol, tcut, lunout, ier)
!
!  One-liner: Parametric local spline: Find t for given x on 3-space curve
!
!  Description:
!
!     PLXCUT3D estimates the value of the parametric variable t where a 3-space
!     curve has x = xcut.  It is the 3-space analog of 2-space PLXCUT.
!
!     A safeguarded Newton iteration is used to solve the nonlinear equation
!                         f(t) = x(t) - xcut = 0.
!
!  History:
!     03/24/95  D.A.Saunders  Initial adaptation of INTSEC2 as PLXCUT (2-space).
!     06/25/21    "     "     Adaptation of PLXCUT for the 3-space case, as
!                             needed for a CAPSULE_GRID extension.
!
!  Author: David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
!          Now with AMA, Inc. at NASA Ames Research Center.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)    :: n       ! # points on curve; n >= 2

   real,    intent (in)    :: x(n), y(n), z(n) ! Curve coordinates

   real,    intent (inout) :: t(n)    ! Cumulative chord lengths, either input
                                      ! or computed here; see calct

   integer, intent (inout) :: ix      ! Estimate of index near the point of
                                      ! intersection; ix = 0 ==> use mid index;
                                      ! Actual index nearest to the point of
                                      ! intersection in INTERVAL's sense (that
                                      ! is, between 1 and n - 1 inclusive).

   logical, intent (in)    :: calct   ! F means the calling program has
                                      ! supplied t(:)

   real,    intent (in)    :: xcut    ! x for which t is desired

   real,    intent (in)    :: tol     ! Tolerance on ||f|| used relative to the
                                      ! the curve length; try 1.E-6 or 1.E-14
                                      ! for 32- or 64-bit arithmetic

   real,    intent (out)   :: tcut    ! Estimated t corresponding to x = xcut

   integer, intent (in)    :: lunout  ! Logical unit for displaying iterations,
                                      ! which are suppressed if lunout < 0;
                                      ! |lunout| is used for error messages

   integer, intent (out)   :: ier     ! 0 means no problem was encountered;
                                      ! 1 means the safeguarded iteration failed
                                      !   to satisfy TOL but did satisfy 10*TOL;
                                      ! 2 means the iteration limit was reached
                                      !   without (close) convergence;
                                      ! 3 means the step-halving failed - somehow
                                      !   the iteration diverged;
                                      ! 4 means the curve was probably parallel
                                      !   to x = xcut
!  Procedures:

   external :: chords3d   ! Chord-length utility
   external :: interval   ! Search utility, also used by LCSFIT
   external :: lcsfit     ! Local spline utility; "loose" fits are used

!  Local constants:

   integer, parameter :: itmax = 15   ! Outer iteration limit; halving the step
                                      ! is the inner iteration, using tol
   real,    parameter :: half = 0.5, one = 1., zero = 0.
   character (1)      :: method = 'B' ! "Bessel" = "loose" spline fits

!  Local variables:

   integer :: iter, lunerr
   real    :: dt, fi, ti, tlast, xi, xp
   real    :: alpha, fnorm, fnorm0, toler, ttotal
   logical :: new

!  Execution:

   lunerr = abs (lunout)
   if (ix == 0) ix = (n + 1) / 2

!  Set up the cumulative chord lengths unless they are supplied by the higher
!  level for efficiency reasons:

   if (calct) then
      call chords3d (n, x, y, z, .false., ttotal, t)
   else
      ttotal = t(n)
   end if

!  Perform a safeguarded Newton iteration:

   ier    = 0
   iter   = 0
   ti     = t(ix)
   new    = .true.
   alpha  = zero           ! For printout purposes only
   fnorm0 = 1.e+20         ! I.e., big to avoid iteration 0 test
   toler  = tol*ttotal

300 continue

!     Evaluate x, dx/dt at the current t:

      call lcsfit (n, t, x, new, method, 1, ti, xi, xp)

      new   = .false.
      fi    = xi - xcut
      fnorm = abs (fi)

!     Halve the step until |f| is reduced (except first time through):

      if (fnorm > fnorm0) then
         if (alpha > tol) then   ! Presume tol > 2 * machine epsilon
             alpha = half*alpha
             ti    = tlast - alpha*dt
             go to 300
         end if
         go to 830
      end if

      if (lunout > 0) then
         write (lunout, '(a, i3, a, es9.2, a, es9.2, a, es16.8)') &
            ' PLXCUT3D:', iter, '  |f|:', fnorm, '  step:', alpha, '  t:', ti
      end if

      if (fnorm < toler) go to 900  ! Success


      iter = iter + 1
      if (iter == itmax) go to 810

      fnorm0 = fnorm

!     Newton step:

      if (xp == zero) go to 840

      dt    = fi/xp
      tlast = ti
      ti    = ti - dt

      alpha = one

   go to 300                           ! Do another iteration


!  Error handling:

810 ier = 1
    write (lunerr, '(a)') 'Iteration limit reached.'
    if (fnorm < 10.*toler) go to 900

    ier = 2
    go to 999

830 ier = 3
    write (lunerr, '(a)') 'Step-halving iteration failed.'
    go to 999

840 ier = 4
    write (lunerr, '(a)') 'Curve parallel to cut?'
    go to 999

900 continue  ! Wrap up, having converged or almost converged

    tcut = ti
    call interval (n, t, ti, one, ix)

999 return

    end subroutine plxcut3d
