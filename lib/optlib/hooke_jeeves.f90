!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  function hooke (nvars, blo, bup, startpt, endpt, rho, step1, eps, itermax,   &
                  fun, f)
!
!    HOOKE seeks a minimizer of a scalar function of several variables.
!    This version allows simple constraints (lower and upper variable bounds).
!
! Summary:
!
!    This version has been adapted at NASA from the version provided at
!      http://people.sc.fsu.edu/~jburkardt/f_src/toms178/toms178.html
!    in order to provide a non-derivative alternative to the gradient-based
!    quasi-Newton method QNMDIF2 and its QNMDRIVER3 framework, including use
!    of the same subroutine fun (...) for controlling the function evaluations
!    during all phases of a sequence of minimizations.  The framework allows
!    all of the application-specifics to be confined to subroutine fun.
!    This version also provides for taking advantage of what should be a good
!    starting guess as a sequence of related optimizations proceeds.
!
! Recent History:
!
!    Drosos Kourounis, Stanford University, suggested this non-derivative method
!    in Aug., 2011 for application to Mars Science Laboratory TPS heating data.
!
!    09/16/11   David Saunders, ERC, Inc./NASA Ames Research Center: initial
!                               adaptation to fun (n, x, f, ncall, fail) form.
!    09/19/11     "      "      Added lower and upper bounds on the variables.
!    09/26/11     "      "      Discussion with Drosos led to the idea of
!                               shortening the initial step length if the
!                               starting guess is believed to be good.  This
!                               requires new argument STEP1.  See usage below.
! Discussion:
!
!    This routine finds a point X where the nonlinear objective function
!    F(X) has a local minimum.  X is an N-vector and F(X) is a scalar.
!    The objective function F(X) is not required to be differentiable
!    or even continuous.  The program does not use or require derivatives
!    of the objective function.
!
!    The user supplies three things:
!    1) a subroutine that computes F(X),
!    2) an initial "starting guess" of the minimum point X,
!    3) values for the algorithm convergence parameters.
!
!    The program searches for a local minimum, beginning from the starting
!    guess, using the Direct Search algorithm of Hooke and Jeeves.
!
!    This program is adapted from the Algol pseudocode found in the paper by
!    Kaupe, and includes improvements suggested by Bell and Pike, and by
!    Tomlin and Smith.
!
!    The algorithm works by taking steps from one estimate of a minimum, to
!    another (hopefully better) estimate.  Taking big steps gets to the minimum
!    more quickly, at the risk of "stepping right over" an excellent point.
!    The stepsize is controlled by a user supplied parameter called RHO.  At
!    each iteration, the stepsize is multiplied by RHO  (0 < RHO < 1), so the
!    stepsize is successively reduced until RHO**m < EPS for some m that depends
!    on the scaling of the variables (they should be O(1)) and the choice of
!    EPS, which should not be smaller than sqrt (machine epsilon).
!
!    Small values of RHO correspond to big stepsize changes,
!    which make the algorithm run more quickly.  However, there
!    is a chance (especially with highly nonlinear functions)
!    that these big changes will accidentally overlook a
!    promising search vector, leading to nonconvergence.
!
!    Large values of RHO correspond to small stepsize changes,
!    which force the algorithm to carefully examine nearby points
!    instead of optimistically forging ahead.  This improves the
!    probability of convergence.
!
!    The stepsize is reduced until it is equal to (or smaller
!    than) EPS.  So the number of iterations performed by
!    Hooke-Jeeves is determined by RHO and EPS:
!
!      RHO^(number_of_iterations) = EPS    ! But see possible use of STEP1 below
!
!    In general it is a good idea to set RHO to an aggressively
!    small value like 0.5 (hoping for fast convergence).  Then,
!    if the user suspects that the reported minimum is incorrect
!    (or perhaps not accurate enough), the program can be run
!    again with a larger value of RHO such as 0.85, using the
!    result of the first minimization as the starting guess to
!    begin the second minimization.
!
!    If there are doubts about the result, the computed minimizer
!    can be used as the starting point for a second minimization attempt.
!
!    In this version, provision is made for taking advantage of what should be
!    a good starting guess in the middle of a sequence of related optimizations.
!    If new input STEP1 < RHO is entered, the sequence of successive stepsize
!    reductions  delta(i) <-- step**m |startpt(i))|  starts with step = STEP1
!    instead of the input RHO, thus reducing the number of iterations.
!    The effect on the number of function evaluations is not so clear.
!
!    To apply this method to data fitting, code your function to be
!    the sum of the squares of the errors (differences) between the
!    computed values and the measured values.  Then minimize it
!    using Hooke-Jeeves.
!
!    For example, you have 20 datapoints (T(i), Y(i)) and you want to
!    find A, B and C so that:
!
!      A*t*t + B*exp(t) + C*tan(t)
!
!    fits the data as closely as possible.  Then the objective function
!    F() to be minimized is just
!
!      F(A,B,C) = sum ( 1 <= i <= 20 )
!        ( y(i) - A*t(i)*t(i) - B*exp(t(i)) - C*tan(t(i)) )^2.
!
! Authors:
!
!    ALGOL original by Arthur Kaupe.
!    C version by Mark Johnson.
!    Fortran 90 version by John Burkardt, distributed as TOMS178, February 2008.
!
! References:
!
!    M Bell, Malcolm Pike,
!    Remark on Algorithm 178: Direct Search,
!    Communications of the ACM,
!    Volume 9, Number 9, September 1966, page 684.
!
!    Robert Hooke, Terry Jeeves,
!    Direct Search Solution of Numerical and Statistical Problems,
!    Journal of the ACM,
!    Volume 8, Number 2, April 1961, pages 212-229.
!
!    Arthur Kaupe,
!    Algorithm 178:
!    Direct Search,
!    Communications of the ACM,
!    Volume 6, Number 6, June 1963, page 313.
!
!    FK Tomlin, LB Smith,
!    Remark on Algorithm 178: Direct Search,
!    Communications of the ACM,
!    Volume 12, Number 11, November 1969, page 637-638.
!
! Parameters:
!
!    Input, integer (kind = 4) NVARS, the number of spatial dimensions.
!
!    Input, real (kind = 8) BLO(NVARS), the lower bounds on the variables;
!    Input, real (kind = 8) BUP(NVARS), the upper bounds on the variables;
!    enter (say) -1.e+30 or +1.e+30 if no useful bounds are available.
!
!    Input, real (kind = 8) STARTPT(NVARS), the user-supplied
!    initial estimate for the minimizer.
!
!    Output, real (kind = 8) ENDPT(NVARS), the estimate for the
!    minimizer, as calculated by this routine.
!
!    Output, real (kind = 8) F, the corresponding minimum function value found.
!
!    Input, real (kind = 8) RHO, a user-supplied convergence parameter
!    which should be set to a value between 0.0 and 1.0.  Larger values
!    of RHO give greater probability of convergence on highly nonlinear
!    functions, at a cost of more function evaluations.  Smaller values
!    of RHO reduce the number of evaluations and the program running time,
!    but increases the risk of missing better minimum estimates.
!
!    Input, real (kind = 8) STEP1, <= RHO as described above.
!
!    Input, real (kind = 8) EPS, the criterion for halting the search for a
!    minimum.  When the algorithm begins to make less and less progress on
!    each iteration, it checks the halting criterion: if the stepsize is below
!    EPS, terminate the iteration and return the current best estimate of the
!    minimum.  Larger values of EPS (such as 1.0e-4) give quicker running time,
!    but a !    less accurate estimate of the minimum.  Smaller values of EPS
!    (such as 1.0e-7) give longer running time, but a more accurate estimate of
!    the minimum.
!
!    Input, integer (kind = 4) ITERMAX, a limit on the number of iterations.
!
!    Input, logical (kind = 4) VERBOSE, turns on print of x(:) and f each itn.
!
!    Input, external FUN, the name of the function subroutine, which should
!    have the form compatible with QNMDIF2 and QNMDRIVER3 for performing a
!    sequence of minimizations:
!
!      subroutine fun (n, x, f, ncall, fail)
!      integer (kind = 4) n
!      real (kind = 8) x(n)
!      real (kind = 8) f
!      integer (kind = 4) ncall
!      logical (kind = 4) fail
!
!      Usage of ncall:
!
!         -5 => first function call for first minimization (reset to -3 in fun);
!         -3 => first function call for subsequent minimizations;
!         -4 => last   "    "    "    "  the current minimization;
!         -2 => end of a search iteration;
!          0 => within a search iteration;
!         99 => stop the sequence of minimizations
!
!    Output, integer (kind = 4) HOOKE, the number of iterations taken.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  implicit none

! Arguments:

  integer (kind = 4) :: hooke        ! Return value is the number of iterations
  integer (kind = 4) :: nvars        ! Assumed known from an external first call
  real    (kind = 8) :: blo(nvars)   ! Lower bounds on the variables, or -"big"
  real    (kind = 8) :: bup(nvars)   ! Upper bounds on the variables, or +"big"
  real    (kind = 8) :: startpt(nvars) ! See descriptions above
  real    (kind = 8) :: endpt(nvars)
  real    (kind = 8) :: rho
  real    (kind = 8) :: step1
  real    (kind = 8) :: eps
  integer (kind = 4) :: itermax
  external           :: fun          ! Function subroutine described above
  real    (kind = 8) :: f            ! Input with the initial function value;
                                     ! output with the minimum function value
! Procedures:

  real    (kind = 8) , external :: best_nearby

! Local constants:

  real   (kind = 8), parameter :: very_big = 1.e+30
  real   (kind = 8), parameter :: zero = 0.e+0

! Local variables:

  real    (kind = 8) delta(nvars)
  logical (kind = 4) fail
  real    (kind = 8) fbefore
  integer (kind = 4) funevals
  integer (kind = 4) i
  integer (kind = 4) iters
  integer (kind = 4) j
  integer (kind = 4) keep
  integer (kind = 4) ncall
  real    (kind = 8) newf
  real    (kind = 8) newx(nvars)
  real    (kind = 8) steplength
  real    (kind = 8) tmp
  real    (kind = 8) xbefore(nvars)

! Execution:

  do i = 1, nvars
    newx(i)    = max (blo(i), min (bup(i), startpt(i)))
    xbefore(i) = newx(i)
  end do

  do i = 1, nvars
    if (startpt(i) == zero) then
      delta(i) = rho
    else
      delta(i) = rho * abs (startpt(i))
    end if
  end do

  steplength = step1  ! Not rho now
  funevals = 0
  iters = 0

!!fbefore = f (newx, nvars)  ! The initial function evaluation is external now

  funevals = funevals + 1
  fbefore  = f
  newf = fbefore

  do while (iters < itermax .and. eps < steplength)

    iters = iters + 1

!   Find best new point, one coordinate at a time.

    do i = 1, nvars
      newx(i) = xbefore(i)
    end do

    newf = best_nearby (nvars, blo, bup, delta, newx, fbefore, fun, funevals)

!   If we made some improvements, pursue that direction.

    keep = 1

    do while (newf < fbefore .and. keep == 1)

      do i = 1, nvars

!       Arrange the sign of DELTA.

        if (newx(i) <= xbefore(i)) then
          delta(i) = -abs (delta(i))
        else
          delta(i) =  abs (delta(i))
        end if

!       Now, move further in this direction.

        tmp = xbefore(i)
        xbefore(i) = newx(i)
        newx(i) = max (blo(i), min (bup(i), newx(i) + newx(i) - tmp))

      end do

      fbefore = newf
      newf = best_nearby (nvars, blo, bup, delta, newx, fbefore, fun, funevals)

!     If the further (optimistic) move was bad...

      if (fbefore <= newf) exit  ! No more in this direction

!     Make sure that the differences between the new and the old points
!     are due to actual displacements; beware of roundoff errors that
!     might cause NEWF < FBEFORE.

      keep = 0

      do i = 1, nvars
        if (0.5e+00 * abs (delta(i)) < abs (newx(i) - xbefore(i))) then
          keep = 1
          exit
        end if
      end do

    end do

    if (eps <= steplength .and. fbefore <= newf) then
      steplength = steplength * rho
      do i = 1, nvars
        delta(i) = delta(i) * rho
      end do
    end if

!   Repeat the final evaluation from this iteration, for possible plotting, etc.

    ncall = -2

    if (fbefore <= newf) then

       call fun (nvars, xbefore, fbefore, ncall, fail)
       if (fail) fbefore = very_big
    else

       call fun (nvars, newx, newf, ncall, fail)
       if (fail) newf = very_big
    end if

  end do  ! Next iteration

  f = newf

  do i = 1, nvars
    endpt(i) = xbefore(i)
  end do

  hooke = iters

  end function hooke

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  function best_nearby (nvars, blo, bup, delta, point, prevbest, fun, funevals)
!
!    BEST_NEARBY looks for a better nearby point, one coordinate at a time.
!
! Modified:
!
!    John Burkardt   12 February 2008  Version distributed as TOMS178.
!    David Saunders  29 August   2011  Made the function QNMDIF2-compatible in
!                                      the form call fun (n, x, f, ncall, fail).
!      "      "      19 Septemb. 2011  Introduced simple variable bounds.
! Author:
!
!    ALGOL original by Arthur Kaupe.
!    C version by Mark Johnson.
!    Fortran 90 version by John Burkardt.
!
! References:
!
!    As for FUNCTION HOOKE.
!
! Parameters:
!
!    Input, integer (kind = 4) NVARS, the number of variables.
!
!    Input, real (kind = 8) BLO(NVARS), the lower bounds on the variables.
!    Input, real (kind = 8) BUP(NVARS), the upper bounds on the variables.

!    Input, real (kind = 8) DELTA(NVARS), the size of a step in each direction.
!
!    Input/output, real (kind = 8) POINT(NVARS); on input, the current
!    candidate.  On output, the value of POINT may have been updated.
!
!    Input, real (kind = 8) PREVBEST, the minimum value of the function seen
!    so far.
!
!    Input, external FUN, the name of the function subroutine, which should
!    have the form compatible with QNMDIF2 and QNMDRIVER3 for performing a
!    sequence of minimizations:
!
!      subroutine fun (n, x, f, ncall, fail)
!      integer (kind = 4) n
!      real (kind = 8) x(n)
!      real (kind = 8) f
!      integer (kind = 4) ncall
!      logical (kind = 4) fail
!
!    Input/output, integer (kind = 4) FUNEVALS, the number of function
!    evaluations.
!
!    Output, real (kind = 8) BEST_NEARBY, the minimum value of the function
!    seen after checking the nearby neighbors.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  implicit none

! Arguments:

  real    (kind = 8) best_nearby
  integer (kind = 4) nvars
  real    (kind = 8) blo(nvars)
  real    (kind = 8) bup(nvars)
  real    (kind = 8) delta(nvars)
  real    (kind = 8) point(nvars)
  real    (kind = 8) prevbest
  external :: fun
  integer (kind = 4) funevals

! Local constants:

  real   (kind = 8), parameter :: very_big = 1.e+30

! Local variables:

  real    (kind = 8) ftmp
  integer (kind = 4) i
  real    (kind = 8) minf
  integer (kind = 4) ncall
  real    (kind = 8) z(nvars)
  logical (kind = 4) fail

! Execution:

  ncall = 0  ! Comparable to QNMDIF2's line search situation

  minf = prevbest

  do i = 1, nvars
    z(i) = point(i)
  end do

  do i = 1, nvars

    z(i) = max (blo(i), min (bup(i), point(i) + delta(i)))

!!  ftmp = f (z, nvars)

    call fun (nvars, z, ftmp, ncall, fail)

    funevals = funevals + 1

    if (fail) ftmp = very_big

    if (ftmp < minf) then

      minf = ftmp

    else

      delta(i) = -delta(i)
      z(i) = max (blo(i), min (bup(i), point(i) + delta(i)))

!!    ftmp = f (z, nvars)

      call fun (nvars, z, ftmp, ncall, fail)

      funevals = funevals + 1

      if (fail) ftmp = very_big

      if (ftmp < minf) then
        minf = ftmp
      else
        z(i) = point(i)
      end if

    end if

  end do

  do i = 1, nvars
    point(i) = z(i)
  end do

  best_nearby = minf

  end function best_nearby
