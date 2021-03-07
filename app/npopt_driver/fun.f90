!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine fun (n, x, f, ncall, fail)
!
!  Subroutines FUN, OBJECTIVE, LINCON, NLCON and SUMMARY should contain all of
!  the application-specifics for performing a sequence of nonlinear constrained
!  optimizations using the NPOPT_DRIVER framework.  Of these, only NLCON is not
!  an internal procedure of FUN, because it is called by both FUN and CONFUN.
!
!  See further usage details below.
!
!  Present Application:
!
!  This multidimensional Rosenbrock function from Wikipedia requires that the
!  number of variables n is even:
!
!                            n/2
!     f(x , x , ..., x )  =  Sum  [100 (x    **2 - x  )**2 + (x     - 1)**2]
!        1   2        n      i=1         2i-1       2i         2i-1
!
!  Control Files:
!
!     namelist.inputs    If present, this file can override the defaults for
!                        each optimization.  If EOF is encountered before the
!                        sequence has been completed, the last controls read
!                        apply to any further optimizations.
!
!     constraint.inputs  In much the same way, if this file is present, it
!                        should contain one or more pairs of namelist-like
!                        tables.  See read_constraints below for details.
!                        If this file is not found, unconstrained optimization
!                        proceeds (with optional bounds on the variables).
!
!  Namelist Inputs:
!
!     epsobj             Estimate of the relative precision of 1 + |obj. fn.|;
!                        e.g., if F(x) has 6 correct significant digits, 1.e-6
!                        is appropriate if F(x) is large, while if F(x) ~ 1.e-4,
!                        then epsobj = 1.e-10 would be appropriate
!     eta                Linesearch tolerance in [0, 1); loose: 0.9; tight: 0.1;
!                        0.9 is the default, suited to expensive objectives
!     major_print_level  The default value of 10 prints one line per major
!                        iteration plus details of the final solution
!     maxn               Total # variables not to be exceeded during any
!                        minimization; must not exceed the hard-coded value of
!                        parameter maxvar in NPOPT_DRIVER (module opt_mod), but
!                        may well be smaller to end a sequence earlier
!     minor_iteration_limit  applies to each QP subproblem
!     mode_gradients     2 => forward diffs. for objective function gradients;
!                        3 => central diffs. (nonlinear constraint grads. too)
!     niter              Number of optimization iterations; niter >= 0; if
!                        niter = 0, there is no need to set noptimizations = 0
!     noptimizations     Number of minimizations to perform (counting the first)
!     step_limit         Limit on variable changes between iterations;
!                        step_limit = 2 (default) allows 200% changes
!     testing            T turns on printing of any diagnostics to lunlog
!     tol_linear         Feasibility tolerance for satisfying linear constraints
!                        and variable bounds at all iterations; the default is
!                        sqrt (epsobj)
!     tol_nonlinear      The largest nonlinear constraint violation that is
!                        acceptable at an optimal point: default: sqrt (epsobj)
!     tol_opt            Optimality tolerance in [machine eps, 1]; the default
!                        is sqrt (epsobj); tol_opt = 1.e-4 (say) means if NPOPT
!                        succeeds, the final objective function should have
!                        approximately 4 correct significant digits
!  Constraint Inputs:
!
!     See read_constraints procedure below.
!
!  History:
!
!  08/23/94  D.A.S.  Initial Rosenbrock fn. (n = 2, QNMDIF, unconstrained).
!  02/18/95  D.A.S.  NCALL argument added for QNMDIF2 variant.
!  05/02/12- D.A.S.  Multidimensional form for demonstrating a sequence of
!  05/30/12          constrained minimizations.  Usage of NCALL is revised.
!  06/01/12  D.A.S.  The FIAT_Opt application showed that variable bounds
!                    and scale factors and differencing intervals may be
!                    set during initialization, before the allocations of
!                    NPOPT work-space.  So we use holding arrays for them
!                    and resort to (maxvar) dimensions, not [unknown] (n).
!  04/25/12  D.A.S.  The major iteration limit for optimizations >= 3 became
!                    the actual iteration count for optimization 2 because of
!                    mishandling of end_of_namelists when only one was present.
!
!  David Saunders, ERC, Inc./Aerothermodynamics Branch, NASA Ames Research Cntr.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Global variables:
!  -----------------

   use opt_mod   ! Optimization variables needed for any NPOPT minimization

   implicit none

!  Arguments:
!  ----------

   integer, intent (inout) :: n      ! No. of optimization variables; see below
   real,    intent (inout) :: x(*)   ! Variables at which to evaluate f
   real,    intent (out)   :: f      ! Function value (unless ncall = -6 or -5)
   integer, intent (inout) :: ncall  ! Input apart from 999 exception (below)
   logical, intent (out)   :: fail   ! T means too many variables were indicated
                                     ! or initialization or evaluation of the
                                     ! function failed in some other way
!  Clarifications:
!  ---------------
!
!     FUN's main purpose is to apply the current optimization variables and
!     evaluate the objective function (which is further isolated in subroutine
!     OBJECTIVE).  FUN uses its ncall argument to control outputs to the log
!     file and to distinguish start- and end-of-optimization situations.
!
!     FUN is called in the main program before NPOPT is called (to perform
!     initializations but no evaluation of objective or constraints) and after
!     NPOPT (to repeat the optimal function and constraint evaluations and
!     perform any wrap-up).
!
!     FUN is also called by OBJFUN unless there are nonlinear constraints, in
!     which case CONFUN evaluates the objective as well as the constraints (to
!     avoid likely repetition in OBJFUN, because NPOPT calls CONFUN first).
!
!     FUN is thus also called by CONFUN if nc_nonlinear > 0.
!
!  Usage of ncall:
!  ---------------
!
!     Both ncall = -6 and -5 are initialization calls preceding an NPOPT call,
!     and neither the objective function nor any constraints should be evaluated
!     in these cases because such evaluations would be repeated when ncall = -3.
!
!     -6 => first FUN call outside NPOPT, before the first minimization
!     -5 => first FUN call outside NPOPT, before any further minimizations;
!     -3 => first FUN call within NPOPT (from CONFUN|OBJFUN), all minimizations;
!     -4 => last (repeat) function call for current minimization;
!     -2 => end of line search;
!     -1 => inactive
!      0 => within line search;
!     >0 => element of gradient;
!    999 on output => stop the sequence of minimizations
!
!  Usage of n:
!  -----------
!
!     If ncall = -6 (very first call to FUN), n should be input with the
!     maximum number of variables allowed for (i.e., maxvar in NPOPT_DRIVER),
!     then output with the actual number of variables to use.
!
!     If ncall = -5 (preceding any further minimization), n should be input as
!     the value from the previous minimization (from which the next value of n
!     may be derived, a value that should also be checked against the original
!     limit, maxvar).
!
!     If ncall = anything else, n is an input that is neither checked nor
!     changed.
!
!  Scaling of the variables, differencing intervals, and upper/lower bounds:
!  -------------------------------------------------------------------------
!
!     xscale(:)  These should be used to make each variable x(:) seen by the
!     xshift(:)  optimizer to be of order 1, roughly (not critical).
!                x(i) * xscale(i) + xshift(i) = physical variable (i).
!                x(i) = (physical variable(i) - xshift(i)) / xscale(i).
!
!     h(:)       Finite differencing intervals for the SCALED variables
!                x(:) seen by NPOPT.  The OBJFUN/CONFUN gradient estimates
!                use h(i) where h(i) * xscale(i) = physical interval (i),
!                         and   h(i) = physical h(i) / xscale(i).
!                Some applications benefit from defining physical h(i) as a
!                fraction of physical variable v(i) if v(i) > 0.
!
!     snlcon(:)  Scale factors applied to the physical nonlinear constraint
!                values in order to produce constraint gradients of order 1.
!
!     bl/bu(:)   Lower and upper bounds on the variables, linear constraints,
!                and nonlinear constraints, should be packed in that order in
!                scaled form before NPOPT sees them.

!  Local constants:
!  ----------------

   real, parameter :: &
      one = 1., zero = 0.

!  Local variables:
!  ----------------

   logical, save :: &
      no_constraints, testing

!  Execution:
!  ----------

   fail = .false.  ! Normally;  set it true if there's a problem

!  Initialization call to FUN outside NPOPT?
!  -----------------------------------------

   if (ncall <= -5) then  ! I.e., -6 or -5: initialize but don't evaluate f

      call initialize_minimization ()  ! Internal procedure below

      if (ncall == 999) then  ! No more minimizations
         if (num_evaluations > 0) then
            write (lunlog, '(/, a, es44.12)') 'Final objective function:', f
         end if
      end if
      go to 99

   end if

!  For every objective function evaluation via OBJFUN or CONFUN:
!  -------------------------------------------------------------

!  Apply the current variables:

!  (Nothing to do in the case of Rosenbrock's function.)

!  Evaluate the objective:

   call objective ()

   if (ncall == -3) then
      obj0  = f                 ! Initial obj. for this minimization
      neval = 0
   end if

   neval           = neval + 1            ! For current minimization
   num_evaluations = num_evaluations + 1  ! Across all minimizations

!  Save/summarize best results so far?

   if (ncall < 0) then  ! Start or end of an iteration

      call summary ()   ! Print the current status; maybe save results

   end if

99 return

!  Internal procedures for subroutine fun:

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine initialize_minimization ()

!     Set up a minimization, either the very first of the run (ncall = -6) or
!     (ncall = -5) another minimization related in some way to the previous one.
!     Some applications won't know when to stop the sequence ahead of time.
!     The number of optimization variables may also change.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local variables:
!     ----------------

      integer       :: i, ios, j
      integer, save :: maxiter, maxn
      logical, save :: end_of_namelists

!     Namelist to be read before each minimization; adapt as needed:
!     --------------------------------------------------------------

      namelist /controls/ &
         epsobj, eta, major_print_level, maxn, minor_iteration_limit, &
         mode_gradients, niter, noptimizations, step_limit, testing,  &
         tol_linear, tol_nonlinear, tol_opt

!     Execution:
!     ----------

!     Very first initialization?
!     --------------------------

      if (ncall == -6) then  ! General-purpose one-time initializations go here

         open  (lunlog, file='Rosenbrock.out', status='unknown')
         write (lunlog, '(/, a)') 'Minimization(s) by NPOPT:'

!        Look for an [optional] namelist control file:

         open (lunnml, file='namelist.inputs', status='old', iostat=ios)
         end_of_namelists = ios /= 0

!        Set some NPOPT defaults and zero some counts:

         bigbnd        = 1.e10  ! +/-bigbnd is substituted for +/-999. inputs
         epsobj        = 5.* epsmch  ! For high-precision F(x) of O(1)
         eta           = 0.2    ! Linesearch tol., [0,1]; loose: 0.9; tight: 0.1
         maxn          = n      ! Should be input with the limit allowed for
         niter         = 50     ! Copied to major_iteration_limit then reused
         maxiter       = niter  ! For initializing minimizations after the first
         major_iteration_limit = niter
         minor_iteration_limit = 1000  ! For each QP subproblem
         major_print_level     = 10    ! 1 line per iteration + final details
         minor_print_level     = 0
         nprint_frequency      = 100                   ! What does this control?
!!       penalty_parameter     = zero                  ! SNOPT only?
         step_limit            = 2.    ! 2 allows a 200% change in the variables
         testing               = .false.  ! Suppress any application diagnostics
         tol_opt               = sqrt (epsobj)  ! 1.e-8 => ~8 digits in last f
         tol_linear            = tol_opt
         tol_nonlinear         = tol_opt

         mode_gradients        = 3  ! Central differencing
         num_evaluations       = 0
         num_iterations        = 0
         num_optimizations     = 0
         noptimizations        = 1  ! The number to do, not the number done

      end if

!     For all initializations:
!     ------------------------

      if (end_of_namelists) then
         niter = maxiter
      else
         niter = -1  ! To detect whether a new value is input or not
         read  (lunnml, nml=controls, iostat=ios)
         if (niter == -1) then
            niter = maxiter
         else
            maxiter = niter
         end if
         end_of_namelists = ios /= 0
      end if

      if (num_optimizations == noptimizations) go to 90  ! # done = # requested

      major_iteration_limit = niter
      maxiter = niter             ! For initializing next minimization

      write (lunlog, '(a)')       ! Show the namelist whether it changed or not
      write (lunlog, nml=controls)

      niter   = 0                 ! Reused as counter for any one minimization
      epsm6th = epsobj**(-1./6.)  ! Applied to h(:) if mode_gradients is 3
      tol_opt = sqrt (epsobj)     ! In case it is changing during the sequence

!     Initialize reading of linear and nonlinear constraints?
!     -------------------------------------------------------

      if (ncall == -6) then  ! After namelist overrides in case testing is now T
         call read_constraints (0, ios)
         if (ios /= 0) then
            fail = .true.
            go to 90
         end if
      end if

!     Application-specific initializations:
!     -------------------------------------

      if (ncall == -6) then

         maxn =  6    ! Override the default for this demonstration case
         n    =  2    ! Size of the first problem in this case
         x(1) = -1.2  ! Starting guesses for Rosenbrock's function
         x(2) =  1.0

      else  ! ncall = -5; set up a next minimization related to the previous one

         n = n + 2    ! This sequence passes through the even values of n
         if (n > maxn) go to 90

         major_iteration_limit = 25 * n

         do i = 1, n / 2  ! Starting guesses
            j = i + i
            x(j-1) = -1.2
            x(j)   =  1.0
         end do

      end if

      blv(1:n)    = -2.
      buv(1:n)    = +2.
      h(1:n)      = 1.e-8  ! For forward differencing
      xscale(1:n) = 1.
      xshift(1:n) = 0.

!     Set up NPOPT's optional inputs (read before first CONFUN/OBJFUN call):
!     ----------------------------------------------------------------------

      open (unit=lnpopt, file='optional.NPOPT.inputs', status='unknown')

      write (lnpopt, '(a)') &
         'Begin', &
         'Nonderivative Line Search'
      write (lnpopt, '(a, i6)') &
         'Major Iteration Limit          ', major_iteration_limit,  &
         'Minor Iteration Limit          ', minor_iteration_limit,  &
         'Major Print Level              ', major_print_level,      &
         'Minor Print Level              ', minor_print_level,      &
         'Print Frequency                ', nprint_frequency,       &
         'Verify Level                   ', -1,                     &
         'Derivative Level               ',  3
      write (lnpopt, '(a, es12.3)') &
         'Function Precision             ', epsobj,            &
         'Infinite Bound                 ', bigbnd,            &
         'Step Limit                     ', step_limit,        &
         'Linear Feasibility Tolerance   ', tol_linear,        &
         'Nonlinear Feasibility Tolerance', tol_nonlinear,     &
         'Optimality Tolerance           ', tol_opt,           &
         'Linesearch Tolerance           ', eta,               &
!!       'Penalty Parameter              ', penalty_parameter, &
         'End'

      rewind (lnpopt)

      call npprms (lnpopt, inform);  close (lnpopt)

      if (inform /= 0) then
         write (lunlog, '(a)') &
            'An error has been detected in the NPOPT options file.'
         fail = .true.
         go to 90
      end if

!     Count the next set of constraints (or recount the last set):

      call read_constraints (1, ios)

      if (testing) write (lunlog, '(/, a, 2i5)') &
         'Numbers of linear & nonlinear constraints:', nc_linear, nc_nonlinear

!     Allocate work-space for NPOPT:
!     ------------------------------

      if (ncall == -5) then
         deallocate (ilcon, inlcon, jnlcon, istate, iwvec)
         deallocate (bl, bllin, blnlin, bu, bulin, bunlin, cjacbak, clambda,   &
                     conlin, cprint, cvec, ffwd, grad, snlcon, tlcon, xnlcon,  &
                     wvec)
         deallocate (amat, cjac, rmat)
         deallocate (lctype, nlctype)
      end if

      nbounds = n + nc_linear + nc_nonlinear  ! For variable + constraint bounds
      nrowa   = max (nc_linear, 1)            ! Max. # linear constraints ...
      nrowj   = max (nc_nonlinear, 1)         ! ... & nonlinear constraints
      nrowr   = max (n, 1)                    ! Row dim. of Cholesky factor R
                                              ! of the Hessian of the Lagrangian
      leniw   = 3*n + nrowa + 2*nrowj         ! Length of integer workspace ...
      lenw    = n*(2*n + nrowa + 2*nrowj + 20) &  ! .... and real workspace
                + 11*nrowa + 21*nrowj
      if (n == 2 .and. nc_linear == 1) lenw = lenw + n  ! Else NPOPT aborts

      allocate (ilcon(nrowa), inlcon(nrowj), jnlcon(nrowj),                    &
                istate(nbounds), iwvec(leniw))
      allocate (bl(nbounds), bllin(nrowa), blnlin(nrowj),                      &
                bu(nbounds), bulin(nrowa), bunlin(nrowj),                      &
                cjacbak(nrowj), clambda(nbounds), conlin(nrowa),               &
                cprint(nrowj), cvec(nrowj), ffwd(n), grad(n),                  &
                snlcon(nrowj), tlcon(nrowa), xnlcon(nrowj), wvec(lenw))
      allocate (amat(nrowa,n), cjac(nrowj,n), rmat(nrowr,nrowr))
      allocate (lctype(nrowa), nlctype(nrowj))

      confun_does_obj = nc_nonlinear > 0

!     Read the next set of constraints (or reread the last set):

      call read_constraints (2, ios)
      if (ios /= 0) then
         fail = .true.
         go to 90
      end if

!     Set up the linear constraint matrix if there is one:
!     ----------------------------------------------------

      if (nc_linear > 0) call lincon (0)  ! Internal procedure below

!     Pack the bounds for the variables and constraints as expected by NPOPT:
!     -----------------------------------------------------------------------

      call pack_bounds ()                 ! Internal procedure

!     More application-specifics may go here:
!     ---------------------------------------

      go to 99

 90   ncall = 999  ! No more minimizations

 99   return

      end subroutine initialize_minimization

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine objective ()
!
!     OBJECTIVE is called by FUN.  It keeps the objective function details in
!     one place.
!
!     This Application:
!
!     Multidimensional form of Rosenbrock's function (see Wikipedia).
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Global variables:

      ! None for this demonstration test case

      implicit none

!     Arguments:

!     Local variables:

      integer :: i, j

!     Execution:

      f = 0.
      do i = 1, n / 2
         j = i + i
         f = f + 100.* (x(j-1)**2 - x(j))**2 + (x(j-1) - 1.)**2
      end do

      end subroutine objective

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine pack_bounds ()
!
!     Pack the bounds for the variables and linear and nonlinear constraints in
!     that order, as expected by NPOPT.  These may well have been set during the
!     initialization before work-space for NPOPT has been allocated.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local variables:

      integer :: l, m

!     Execution:

      bl(1:n) = blv(1:n)
      bu(1:n) = buv(1:n)

      l = n + nc_linear
      if (nc_linear > 0) then
         bl(n+1:l) = bllin(1:nc_linear)
         bu(n+1:l) = bulin(1:nc_linear)
      end if

      m = l + nc_nonlinear
      if (nc_nonlinear > 0) then
         bl(l+1:m) = blnlin(1:nc_nonlinear)
         bu(l+1:m) = bunlin(1:nc_nonlinear)
      end if

      end subroutine pack_bounds

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine read_constraints (mode, ios)
!
!     Purpose:
!
!        Constraint tables of the formats shown can be scanned (mode = 1) and
!        read (mode = 2) with a pair of calls to this routine separated by
!        allocations of the relevant arrays found in module opt_mod.
!
!     Format of 'constraint.inputs':
!
!        $linear_constraints
!        i       lctype     bl     bu    tlcon  ilcon
!        1   'SUM     '     2.     2.     999.    999
!        2   '????????'     0.   999.     999.    999
!        $end [linear_constraints]
!   
!        Park any inactive linear constraints here.
!   
!        $nonlinear_constraints
!        i      nlctype     bl     bu   xnlcon inlcon jnlcon snlcon
!        1   'PRODUCT '     1.     1.     999.    999    999     1.
!        2   '????????'  -999.     0.     999.    999    999     1.
!        $end [nonlinear_constraints]
!   
!        Park any inactive nonlinear constraints here.
!   
!        $linear_constraints
!        i       lctype     bl     bu    tlcon  ilcon
!        ::::::::::::::::::::::::::::::::::::::::::::    Any more as needed
!        $end
!   
!        $nonlinear_constraints
!        i      nlctype     bl     bu   xnlcon inlcon jnlcon snlcon
!        ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!        $end
!   
!        ::::::::::::::::::::::::::::::::::::::::::::
!
!     Guidelines:
!
!        Pairs of linear and nonlinear constraint sets are delimited as though
!        they were namelists.  A pair is looked for as part of initializing each
!        minimization.  A line is considered active if it appears between
!        $[non]linear_constraints and $end.  (Truly) blank lines anywhere
!        are permitted.
!
!        Inactive constraints can be retained after an $end line.
!
!        If fewer pairs of constraint sets appear than the number of optimiza-
!        tions performed, the last set present will be reused.
!
!        Constraint names lctype and nlctype are limited to 8 characters, left-
!        justified.  These names should be handled by application-specific sub-
!        routines LINCON and NLCON.
!
!        Bound inputs of +/-999. are set to the value of namelist input bigbnd;
!        other inputs of 999[.] are suggested to indicate that the available
!        integer or real controls (ilcon, tlcon, i/jnlcon, xnlcon) are unused.
!
!     Strategy:
!
!        The numbers of constraints are counted as opposed to being indicated
!        explicitly.  This is typically done by counting lines then rewinding
!        the file.  However, for a sequence of such scans and reads, it appears
!        preferable to avoid repeated rewinding by reading ALL constraint sets
!        into a character array and accessing that directly.
!
!     History:
!
!        05/25/12  D.A.Saunders  Initial implementation.
!
!     Author:  David Saunders, ERC, Inc./Aerothermodynamics Branch, NASA Ames.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Global variables:

      use opt_mod     ! No need - this internal procedure inherits the module

!     Arguments:

      integer, intent (in)  :: mode  ! 0 => read and store all sets of
                                     !      constraints (only)
                                     ! 1 => count the numbers of constraints in
                                     !      the next set pair; if no more sets
                                     !      were present, reuse the last sets
                                     ! 2 => read (internally) the indicated pair
                                     !      of constraints after allocating the
                                     !      storage
      integer, intent (out) :: ios   ! 0 => no problem detected;
                                     ! 1 => read error - probably time to abort
!     Local constants:

      integer,   parameter :: maxlength    = 80   ! For any meaningful line
      character, parameter :: header*4     = 'lcon'
      character, parameter :: endset*2     = '$e'
      character, parameter :: startlcon*2  = '$l'
      character, parameter :: startnlcon*2 = '$n'
      character, parameter :: filename*17  = 'constraint.inputs'
      character, parameter :: subname*18   = 'READ_CONSTRAINTS: '

!     Local variables:

      integer                      :: i, ic, j, l ! Short index variables
      integer                      :: ie, il, in  ! For detecting '$e','$l','$n'
      integer                      :: ih          ! For detecting header lines
      integer,   save              :: lastlin     ! Line # of last $linear...
      integer,   save              :: lastnlin    ! Line # of last $nonlinear...
      integer,   save              :: line        ! Next line to read
      integer,   save              :: line_keep   ! Mode 2 starts as for Mode 1
      integer,   save              :: nlines      ! # meaningful lines found
      real                         :: boundlo     ! Obvious
      real                         :: boundup
      logical                      :: counting    ! T if line is meaningful
      logical                      :: store       ! Saves repeating buffer copy
      character                    :: buffer*(maxlength)
      character, save, allocatable :: constraints(:)*(maxlength)

!     Execution:

      select case (mode)

         case (0)  ! Load all sets of constraints into memory

            open (luncon, file=filename, status='old', iostat=ios)

            if (ios /= 0) then
               write (lunlog, '(3a)') subname, 'Unable to open ', filename
               write (lunlog, '(a)') 'Proceeding in unconstrained mode.'
               no_constraints = .true.
               nc_linear      = 0
               nc_nonlinear   = 0
               ios = 0
               go to 99  ! Single-return philosphy
            end if

!           Count the number of meaningful lines in the file:

            nlines = 0
            counting = .false.

            do  ! Until EOF
               read (luncon, '(a)', iostat=ios) buffer
               if (ios < 0) exit

               ih = 0;  ie = 0;  il = 0;  in = 0
               l = len_trim (buffer)
               if (l > 0) then
                  ih = index (buffer(1:l), header)
                  ie = index (buffer(1:l), endset)
                  il = index (buffer(1:l), startlcon)
                  in = index (buffer(1:l), startnlcon)
               end if

               if (ih > 0) then
                  ! Skip header lines
               else if (il > 0 .or. in > 0) then
                  nlines = nlines + 1
                  counting = .true.
               else if (ie > 0) then
                  nlines = nlines + 1
                  counting = .false.
               else if (counting) then
                  nlines = nlines + 1
               end if
            end do

            if (testing) write (lunlog, '(a, i5)') &
               'Number of meaningful constraint input lines:', nlines

!           Now reread the lines and store them:

            rewind (luncon)
            allocate (constraints(nlines))
            nlines = 0
            counting = .false.

            do  ! Until EOF
               read (luncon, '(a)', iostat=ios) buffer
               if (ios < 0) then
                   ios = 0
                   exit
               end if

               ih = 0;  ie = 0;  il = 0;  in = 0
               l = len_trim (buffer)
               if (l > 0) then
                  ih = index (buffer(1:l), header)
                  ie = index (buffer(1:l), endset)
                  il = index (buffer(1:l), startlcon)
                  in = index (buffer(1:l), startnlcon)
               end if

               store = .false.
               if (ih > 0) then
                  ! Skip header lines
               else if (il > 0) then
                  nlines   = nlines + 1
                  counting = .true.
                  store    = .true.
                  lastlin = nlines
               else if (in > 0) then
                  nlines   = nlines + 1
                  counting = .true.
                  store    = .true.
                  lastnlin = nlines
               else if (ie > 0) then
                  nlines   = nlines + 1
                  counting = .false.
                  store    = .true.
               else if (counting .and. l > 0) then
                  nlines = nlines + 1
                  store  = .true.
               end if
               if (store) constraints(nlines) = buffer
            end do

            close (luncon)

            no_constraints = nlines == 0

            if (testing) write (lunlog, '(/, (a))') &
               (trim (constraints(i)), i = 1, nlines)

            line = 0     ! Next line to read

         case (1)  ! Count the constraints in the next set (or reuse last set)

            if (no_constraints) then
               ios = 0
               go to 99
            end if

            line = line + 1  ! Non-meaningful lines have been suppressed
            if (line > nlines) line =  lastlin  ! Reread the last set
            line_keep = line   ! For mode = 2
            nc_linear = 0
            do  ! Until $end
               line = line + 1
               ie   = index (constraints(line), endset)
               if (ie > 0) exit
               nc_linear = nc_linear + 1
            end do

            line = line + 1
            if (line > nlines) line =  lastnlin
            nc_nonlinear = 0
            do  ! Until $end
               line = line + 1
               ie   = index (constraints(line), endset)
               if (ie > 0) exit
               nc_nonlinear = nc_nonlinear + 1
            end do

!           The calling program should now allocate constraint storage.

         case (2)  ! Read (internally) the indicated set(s) of constraints

            if (no_constraints) then
               ios = 0
               go to 99
            end if

            line = line_keep
            nc_linear = 0
            do  ! Until $end
               line = line + 1
               ie   = index (constraints(line), endset)
               if (ie > 0) exit
               nc_linear = nc_linear + 1
               i = nc_linear
               read (constraints(line), *, iostat=ios) &
                  ic, lctype(i), boundlo, boundup, tlcon(i), ilcon(i)
               if (ios /= 0) then
                  write (lunlog, '(a)') &
                     'Error reading linear constraints:', constraints(line)
                  go to 90
               end if
               if (int (boundlo) == -999) boundlo = -bigbnd
               bllin(i) = boundlo
               if (int (boundup) ==  999) boundup =  bigbnd
               bulin(i) = boundup
            end do

            line = line + 1
            if (line > nlines) line =  lastnlin  ! Reread the last set 
            nc_nonlinear = 0
            do  ! Until $end
               line = line + 1
               ie   = index (constraints(line), endset)
               if (ie > 0) exit
               nc_nonlinear = nc_nonlinear + 1
               j = nc_nonlinear
               read (constraints(line), *, iostat=ios) &
                  ic, nlctype(j), boundlo, boundup, xnlcon(j), inlcon(j), &
                  jnlcon(j), snlcon(j)
               if (ios /= 0) then
                  write (lunlog, '(a)') &
                     'Error reading nonlinear constraints:', constraints(line)
                  go to 90
               end if
               if (int (boundlo) == -999) then
                  boundlo = -bigbnd
               else
                  boundlo = boundlo * snlcon(j)
               end if
               blnlin(j) = boundlo
               if (int (boundup) ==  999) then
                  boundup =  bigbnd
               else
                  boundup = boundup * snlcon(j)
               end if
               bunlin(j) = boundup
            end do

         case default  ! Programming error to get here

      end select

      go to 99

 90   ios = 1  ! Trouble
 99   return

      end subroutine read_constraints

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine summary ()
!
!     Summarize the initial and optimized result and save.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Global variables:

      ! None for this demonstration test case

      implicit none

!     Local constants:

      logical, parameter :: print = .true.

!     Execution:

      if (ncall == -2) then
         num_iterations = num_iterations + 1  ! Total across all minimizations
         niter          = niter          + 1  ! Total for this minimization
      end if

      if (ncall /= -4) then
         if (num_optimizations == 0) write (lunlog, 90) niter
         if (num_optimizations >= 1) write (lunlog, 91) niter, num_iterations
         write (lunlog, 92) &
            obj0, f
!!          obj0, f, &
!!          obj0, f     ! Print other items of interest here

!        Linear constraints can be handled with an internal procedure:

         if (nc_linear > 0) call lincon (1)  ! More readable print than NPOPT's

!        Nonlinear constraints cannot, because CONFUN also calls NLCON:

         if (nc_nonlinear > 0) call nlcon (n, x, print, cprint) ! Evaluate/print

      else  ! ncall == -4 -- last repeat evaluation for this minimization

         num_optimizations = num_optimizations + 1
         write (lunlog, '(/, a, i6, i7, /, a, i6, i7, /, a, i4)') &
            '# iterations & objective evaluations, this minimization:',   &
            niter, neval, &
            'Total # iterations & objective evaluations, this run:   ',   &
            num_iterations, num_evaluations, &
            '# minimizations performed:', num_optimizations
      end if

90    format (/, 'STATUS AT ITERATION', I4)
91    format (/, 'STATUS AT ITERATION', I4, '  (TOTAL ITERATIONS:', I5, ')')
92    format (42X, ' Initial', 12X, 'Current', /, &
         ' Objective function            ', 2es19.10)
!!       ' Objective function            ', 2es19.10, /, &
!!       ' Some other quantity           ', 2es19.10)

      end subroutine summary

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine lincon (mode)
!
!     Set up the linear constraint matrix if mode = 0, or print the status of
!     the linear constraints at each major iteration if mode = 1.  NPOPT ensures
!     that the linear constraints are always satisfied, but this printout shows
!     that in more readable form.  It is application-specific.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in) :: mode  ! See above

!     Local variables:

      integer :: j, l

!     Execution:

      select case (mode)

         case (0)  ! Set up the linear constraint matrix

            amat(:,:) = zero  ! (1:nc_linear,1:n)

            do l = 1, nc_linear

               select case (lctype(l))

                  case ('SUM     ')  ! Sum of variables = n (trivial example)

                     do j = 1, n
                        amat(j,1) = 1.
                     end do

                  case default

               end select

            end do

         case (1)  ! Print the linear constraint status

            conlin = matmul (amat, x(1:n))

            write (lunlog, '(/, 2a)') &
               ' #  LCTYPE     TLCON  ILCON      BL      BU           AMAT * V'

            do l = 1, nc_linear
               write (lunlog, '(i2, 2x, a, f8.2, i7, 2f8.2, es19.10)') &
                 l, lctype(l), tlcon(l), ilcon(l), bllin(l), bulin(l), conlin(l)
            end do

         case default  ! This would be a programming error

      end select

      end subroutine lincon

   end subroutine fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine nlcon (nvpar, vpar, print, cpar)
!
!  NLCON evaluates all nonlinear constraints, CPAR(*), with the option to
!  tabulate UNscaled results.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Global variables:

   use opt_mod

   implicit none

!  Arguments:

   integer, intent (in)  :: nvpar       ! The number of optimization variables
   real,    intent (in)  :: vpar(nvpar) ! Current opimization variable values
   logical, intent (in)  :: print       ! T means tabulate UNscaled constraints
   real,    intent (out) :: cpar(*)     ! Nonlinear constraint values, scaled by
                                        ! SNLCON(:) which should be set to give
                                        ! constraint gradients of order 1
!  Local constants:

   real, parameter :: zero = 0.

!  Local variables:

   integer   :: i, k, m
   real      :: boundl, boundu, scale, time, value
   character :: flag*3

!  Execution:

   k = nvpar + nc_linear ! Offset for the bounds

   if (print) write (lunlog, 80)

   do m = 1, nc_nonlinear

      select case (nlctype(m))

         case ('PRODUCT ')   ! Simple example
            value = product (vpar(1:nvpar))  ! Physical-space value

         case ('SUM     ')
            value = sum (vpar(1:nvpar))

         case default
            write (lunlog, '(/, 2a)') ' NLCON: bad constraint type = ', &
               nlctype(m)
            stop

      end select

      scale   = snlcon(m)
      cpar(m) = value * scale

      if (print) then
         k = k + 1
         boundl = bl(K) / scale
         boundu = bu(k) / scale

         if (value < boundl .or. value > boundu) THEN
            flag = '***'
         else
            flag = '   '
         end if

         write (lunlog, 85) m, nlctype(m), value, boundl, boundu, &
                            xnlcon(m), inlcon(m), flag
      end if

   end do

80 format (/, ' #  NLCTYPE      VALUE    LOWER BOUND UPPER BOUND', &
           '    XNLCON  INLCON  VIOLATION?')
85 format (i2, 2x, a, f13.6, 2f12.5, f10.4, i8, 5x, a3)

   end subroutine nlcon
