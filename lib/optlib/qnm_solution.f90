!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine qnm_solution (fun, eachstep, fail)
!
!  This is a subroutine form of the essence of the earlier program QNMDRIVER2,
!  which is a reusable framework for driving the QNMDIF2 quasi-Newton method of
!  unconstrained optimization with all the application specifics confined to the
!  objective function routine FUN (N, X, F, NCALL, FAIL) and the "EACHSTEP"
!  routine with its long calling sequence, q.v.  In fact, QNMDRIVER2 now calls
!  QNM_SOLUTION as the reusable portion that could find application in other
!  programs requiring from one to many optimizations in one run, where the
!  optimizations aren't necessarily related in some way.  (For doing a sequence
!  of related optimizations, see QNMDRIVER3.)
!
!  History:
!
!     07/03/13  D.A.Saunders  Subroutine form of the QNMDRIVER2 program used
!                             by the first (unconstrained) FIAT_Opt scheme,
!                             which moved xscale(:), xshift(:) and h(:) to
!                             common /QNMCTL/ where FUN can see them.
!     07/11/13    "      "    Time to abandon the common block: use a module.
!
!  Author:  David Saunders, ERC, Inc./NASA Ames Research Center, Moffett Fld. CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use qnm_module  ! Replaces the original common block for numerous controls

   implicit none

!  Arguments:

   external :: fun ! Objective function routine with the above calling sequence
                   ! that uses the ncall argument to determine its action:
                   !
                   ! ncall = -4 means QNMDIF2 is finished; FUN is being called
                   !            one more time to reevaluate the objective at
                   !            the optimal Xs returned by QNMDIF2;
                   ! ncall = -3 means initialize a new optimization;
                   ! ncall = -2 means FUN is being called at the end of a step;
                   ! ncall = -1 means FUN is being called during a local search;
                   ! ncall =  0 means FUN is being called during a line search;
                   ! ncall >  0 means FUN is being called to estimate element
                   !            ncall of the gradient;

   external :: eachstep ! A subroutine of the form expected by QNMDIF2 for its
                        ! USER argument.  If QNM_SOLUTION is being called many
                        ! times in one run, EACHSTEP may do no more than call
                        ! FUN at the end of a line search (ncall = -2), but this
                        ! is up to the application for monitoring progress, etc.

   logical, intent (out) :: fail ! It's up to the application to know what to do

!  General purpose utilities:

   external :: cendif2  ! Supervises calculation of good finite difference
                        ! intervals & initializes the gradient vector and
                        ! the diagonal factor of the Hessian approximation
   external :: qnmdif2  ! Quasi-Newton unconstrained minimization utility
                        ! using finite difference gradient approximations
!  Local variables:

   integer :: &
      istart, n, ncall, nfcen, nldim, ntype
   integer, dimension (maxvar) :: &
      ierror, nfdel
   real :: &
      f, gtp, gnm
   real, dimension (maxvar) :: &
      d, ffwd, g, oldg, p, pwork, sh, urh, x, y, z
   real :: &
       l(maxvar*(maxvar - 1)/2)
   logical :: &
      unitl, conv, lfail(maxvar)

!  Execution:
!  !!!!!!!!!!

!  Default certain arguments to QNMDIF2.  These can be overridden within FUN.

   contin = .true.         ! See qnm_module or QNMDIF2 for the descriptions
   epsmch = epsilon (epsmch)
   eta    = 0.1
   igrad  = 0
   local  = .false.
   lunlog = 6
   nallow = 20000
   rescue = .false.
   zeta   = 10.
   nftotl = 0

   ncall  = -3             ! Initial function evaluation for this minimization
   n      = maxvar         ! FUN (ncall = -3) sets # vars. n to the right value

   call fun (n, x, f, ncall, fail)

   if (fail) go to 90

   nftotl = 1

   if (niter > 0) then

!     Estimate good finite differencing intervals, and set up the initial
!     gradient and diagonal elements of the Hessian-like matrix.

      call cendif2 (n, x, f, epsobj, h, g, d, nfcen, nallow, zeta, lunlog,     &
                    fun, urh, nfdel, ierror)

      nftotl = nftotl + nfcen
      ntype  = 1              ! A central diff. gradient is input to QNMDIF2
      unitl  = .false.
      istart = 0
      nldim  = n*(n - 1)/2

      l(1:nldim) = 0.

!     Perform the minimization:

      call qnmdif2 (n, nldim, nftotl, niter, ntype, lunlog, x, f, lowbound, g, &
                    h, l, d, eta, tol, stepmx, epsmch, epsobj, unitl, local,   &
                    conv, contin, rescue, print, fun, eachstep, oldg, p, ffwd, &
                    sh, y, z, pwork, lfail, gtp, gnm, zeta, nallow, igrad,     &
                    istart)

!     Repeat the optimal function evaluation, to wrap up:

      ncall = -4

      call fun (n, x, f, ncall, fail)

      if (fail) go to 90

   end if

   go to 99

90 write (lunlog, '(/, a, i4)') &
      'QNM_SOLUTION:  Objective function evaluation failed.  NCALL:', ncall

99 return

   end subroutine qnm_solution
