
!!!!!!                    Program NPOPT_DRIVER modules                    !!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   module opt_mod               ! Variables needed by NPOPT for any problem but
                                ! set up in an initialization call to FUN
      integer, parameter :: &
         lunlog = 1,        &   ! For main program & FUN; NPOPT writes to unit 6
         lunnml = 2,        &   ! For [optional] namelist.inputs file
         luncon = 3,        &   ! For constraint.inputs file
         lnpopt = 4,        &   ! Temporary file for NPOPT's keyword controls
         maxvar = 200           ! Limit on the number of optimization variables

!     NOTE:  The optimization variables v(:) cannot be allocatable because they
!            are an argument to FUN which may determine how many there are at
!            the start of each optimization. Some reasonable limit is OK because
!            the dense matrix method of NPOPT is not intended for big problems.

      integer :: &
         inform,              & ! NPOPT's returned status
         lencw, leniw, lenw,  & ! Work-space lengths set in initialize_min.
         max_iteration_total, & ! Limit on sum of all minor (QP) iterations
         major_iteration_limit, & ! Limit on the number of NPOPT main iterations
         minor_iteration_limit, & ! Limit for each QP subproblem; set it large
         major_print_level,   & ! 10 gives one line per major itn. + final soln.
         minor_print_level,   & ! 0 suppresses all QP subproblem printing
         mode_gradients,      & ! 2 = 2-point forward differencing; 3 = central
         nbounds,             & ! Derived total for variables and constraints
         nc_linear,           & ! Input number of linear constraints (0 for now)
         nc_nonlinear,        & ! Input number of nonlinear constraints
         neval,               & ! Number of objective fn. evaluations (1 opt.)
         niter,               & ! Number of major iterations (1 optimization)
         noptimizations,      & ! Requested number of calls to NPOPT
         nprint_frequency,    & ! What does this affect?
         nrowa, nrowj, nrowr, & ! Matrix row dimensions set by initialize_min.
         num_evaluations,     & ! FUN should increment these appropriately
         num_iterations,      & ! Sum of niter over all optimizations
         num_optimizations,   & ! Number of optimizations performed so far
         verify_level           ! Allows checking of derivatives by DNOPT

      real :: &
         bigbnd,            &   ! Value replacing +/-999. inputs as "infinite"
         epsmch,            &   ! Machine epsilon computed via epsilon (epsmch)
         epsm6th,           &   ! Value applied to 2-point h to get 3-point h
         epsobj,            &   ! Min. significant change in objective function
         eta,               &   ! Line search tol., [0,1]; 0.2 if expensive obj.
         h_all_x,           &   ! Finite diff. intervals for all scaled x vars.
         obj0,              &   ! Initial objective function for this NPOPT call
         obj_confun,        &   ! CONFUN evaluates objective if nc_nonlinear > 0
         obj_scale,         &   ! Multiplier applied to the objective function
         penalty_parameter, &   ! > 0. may help start a difficult optimization
         step_limit,        &   ! E.g., 2.0 allows a 200% change in a variable
         tol_linear,        &   ! Largest abs. violation at a "feasible" point
         tol_nonlinear,     &   ! Largest acceptable nonlinear constraint violn.
         tol_opt                ! No. of correct figures sought in the optimized
                                ! objective
      real, dimension (maxvar) :: &
         blv,               &   ! Holding space for variable bounds, which may
         buv,               &   ! be set during initializn. before bl/bu allocn.
         h,                 &   ! Intervals for differencing w.r.t. opt. vars.
         v,                 &   ! Optimzn. vars., passed to FUN as its "x" arg.
         xscale,            &   ! Scale and shift inputs making the ...
         xshift                 ! ... variables seen by NPOPT to be ~ 1

      logical :: &
         confun_does_obj        ! T if nc_nonlinear > 0

      integer, allocatable, dimension (:) :: &
         ilcon,             &   ! Input index usable for each linear constraint
         inlcon,            &   ! Two indices used (or not) for each ...
         jnlcon,            &   ! ... nonlinear constraint
         istate,            &   ! Need not be initialized for a cold NPOPT start
         iwvec                  ! Work-space passed to NPOPT

      real, allocatable, dimension (:) :: &
         bl, bu,            &   ! Bounds for vars. + lin. + nonlin. constraints
         bllin, bulin,      &   ! Separate copy of linear constraint bounds ...
         blnlin, bunlin,    &   ! ... and nonlinear constraint bounds
         cjacbak,           &   ! Work-space for constraint gradients
         clambda,           &   ! Work-space for NPOPT, not input for cold start
         conlin,            &   ! Work-space for evaluating linear constraints
         cprint,            &   ! Work-space for evaluating nonlinear constrnts.
         cvec,              &   ! Work-space for NPOPT's output argument C
         ffwd,              &   ! Work-space for objective fn. finite diffs.
         grad,              &   ! Work-space for NPOPT's output obj. grad. arg.
         snlcon,            &   ! Inputs scaling nonlinear constr. grads. to ~ 1
         tlcon,             &   ! Inputs used (or not) for linear constraints
         xnlcon,            &   ! Inputs used (or not) for nonlinear constraints
         wvec                   ! Work-space for NPOPT; no need to initialize

      real, allocatable, dimension (:,:) :: &
         amat,              &   ! Input linear constraint matrix for NPOPT
         cjac,              &   ! Work-space for nonlinear constraint gradients
         rmat                   ! Work-space for NPOPT (Cholesky factor)

      character (8), allocatable, dimension (:) :: &
         cwvec,             &   ! Work-space for DNOPT; no need to initialize
         lctype,            &   ! For specifying linear constraints by name
         nlctype                ! For specifying nonlinear constraints by name

   end module opt_mod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program npopt_driver
!
!  Prologue:
!
!        This adaptation of QNMDRIVER3 (sequence of unconstrained optimizations
!     using QNMDIF2) and the HEAT_SHIELD constrained optimization application
!     drives a sequence of constrained optimizations using the NPOPT sequential
!     quadratic programming package.  Analytic gradients are assumed not to be
!     available: both schemes employ finite difference approximations to first
!     derivatives, so the user needs to supply only function values for the
!     objective and any nonlinear constraints.  Linear constraints are also
!     supported, along with simple bounds on the variables.
!
!        A sequence of related optimizations may arise from applying the same
!     variables and objective function to a series of similar datasets, or from
!     applying increasing numbers of variables to the same dataset.
!
!        NPOPT expects subroutines OBJFUN and CONFUN as two of its arguments.
!     A complication is that, when nonlinear constraints are present, NPOPT
!     calls CONFUN before it calls OBJFUN, so it is desirable to arrange for
!     CONFUN to evaluate the objective function too (since the constraints and
!     the objective are likely to need the same underlying calculations).
!     This has been implemented in a way that keeps both OBJFUN and CONFUN
!     independent of the application, through calls to a "FUN" routine that
!     is actually in the same form as that used by the unconstrained scheme.
!
!        Ideally, all of the problem specifics are confined to fun.f90, of which
!     a straightforward example accompanies this framework.  Both FUN and CONFUN
!     call a lower-level NLCON to evaluate/print any nonlinear constraints, and
!     it makes sense for NLCON to be part of fun.f90 following the internal
!     procedures of FUN (which include linear constraint subroutine LINCON).
!
!        The calling sequences are:
!
!        SUBROUTINE  CONFUN (MODE, NC, NV, NROWJ, NEEDC, V, C, CJAC, NSTATE)
!        SUBROUTINE  OBJFUN (MODE, NV, V, OBJ, GRAD, NSTATE)
!        SUBROUTINE     FUN (NV, V, OBJ, NCALL, FAIL)
!        SUBROUTINE   NLCON (NV, V, PRINT, C)
!
!        All of these subroutines use the OPT_MOD module above, which contains
!     the numerous application-independent variables required for any sequence
!     of optimizations using NPOPT.  It should be possible to confine further
!     application-specific modules to FUN and NLCON.
!
!        Note that the optimization variables themselves are variously named
!     either V or X for historical reasons.
!
!        See a sample FUN routine for usage of the NCALL argument to FUN, and
!     for its use of 'namelist.inputs' and 'constraint.inputs' control files.
!
!  Abbreviated History:
!
!     QNMDIF is an unconstrained quasi-Newton method using finite-difference
!     derivative approximations and modified Cholesky factors of the Hessian-
!     like matrix to ensure descent directions at each step.  QNMDIF2 differs
!     only in being full argument-driven (no local arrays).  Reference:
!     National Physical Laboratory Report NAC 11, "The Implementation of Two
!     Revised Quasi-Newton Algorithms for Unconstrained Optimization" by
!     P.E.Gill, W.Murray, and R.A.Pitfield (1972).
!
!     NPOPT has the same calling sequence as the better-known NPSOL, from which
!     it differs only in the details of the QP solver used for the SQP scheme
!     (sequential quadratic programming algorithm).  The version employed at
!     NASA Ames Research Center is Version 1.0-10, 1995.  Reference:
!     Systems Optimization Laboratory Technical Report SOL 86-2, "User's Guide
!     for NPSOL (Version 4.0): A Fortran Package for Nonlinear Programming" by
!     Philip E. Gill, Walter Murray, Michael A. Saunders and Margaret H. Wright
!     (Department of Operations Research, Stanford University, 1986).
!
!     Robert Kennelly, James Reuther, and David Saunders contributed to various
!     applications of these optimizers at NASA Ames Research Center.
!
!     05/12/92  D.A.Saunders  Framework developed as QNMDRIVER (unconstrained),
!                             for initial application to B-spline mysteries.
!     08/23/94    "    "      QNMDRIVER2 version for testing mods. to CENDIF2
!                             and QNMDIF2.  Specifics are confined to "FUN".
!     05/03/11    "    "      QNMDRIVER3 variant performs a sequence of related
!                             optimizations.  This required NCALL = -5 usage.
!     05/02/12 -  "    "      Adaptation of QNMDRIVER3 and HEAT_SHIELD as the
!     05/30/12                NPOPT_DRIVER framework for a sequence of (one or
!                             more) constrained optimizations.  This required
!                             use of NCALL = -6.
!     06/01/12    "    "      Introduced holding space for variable bounds. They
!                             may need to be set during initialization before
!                             NPOPT work-space has been allocated.  For the same
!                             reason, other arrays like h and xscale that are
!                             associated with the variables are too awkward to
!                             make allocatable, so use maxvar as for x/v.
!     06/22/16    "    "      Slight changes to opt_mod to match what was done
!                             in DNOPT_DRIVER recently.
!
!  Author:  David Saunders, ERC, Inc./Aerothermodynamics Branch, NASA Ames, CA.
!                 Now with  AMA, Inc.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use opt_mod  ! Optimization variables common to any NPOPT problem

   implicit none

!  Local variables:

   integer :: &
      i, major_iters, n, ncall
   real :: &
      f, obj_final       ! NPOPT outputs obj_final and major_iters, but FUN is
   logical :: &          ! most likely displaying these another way
      fail
   character :: &
      progname*14 = 'NPOPT_DRIVER: '

!  Procedures:

   external :: &
      confun,  &         ! Nonlinear constraint function, application-specific
      fun,     &         ! Should be able to contain most application-specifics
      objfun,  &         ! Objective function routine; like CONFUN, it calls FUN
      npopt              ! Nonlinear constrained optimization (Stanford Univ.)


!  Execution:
!  ----------

   epsmch = epsilon (epsmch)

!  Loop over one or more minimizations, relying on FUN to flag when the loop
!  is finished (NCALL = 999).  Unlike QNMDIF2, which expects the objective
!  function to have been evaluated at the starting guess before it is called,
!  NPOPT itself performs the first call(s) to CONFUN and/or OBJFUN after it
!  does its own initialization.  These in turn will call FUN appropriately.
!  Therefore, in order to confine the application-specifics to FUN, we still
!  make its initial call outside of NPOPT so we can set up for NPOPT's own
!  initialization.  This initial FUN call (for each minimization) should avoid
!  evaluating the initial objective function and any constraints because those
!  are done within NPOPT (even if niter = 0).

   ncall = -6  ! First call outside NPOPT, before the first minimization
   n = maxvar  ! FUN should set n <= maxvar correctly before each NPOPT call

10 continue

      call fun (n, v, f, ncall, fail)  ! Initialize (only - no objective eval.)

      if (fail) then
         write (lunlog, '(/, 2a)') progname, 'Initialization problem in FUN.'
         go to 999
      end if

      if (ncall == 999) go to 999  ! No more minimizations

!     Perform the (current) minimization:

      call npopt (n, nc_linear, nc_nonlinear, nrowa, nrowj, nrowr, &
                  amat, bl, bu, confun, objfun, inform, major_iters, &
                  istate, cvec, cjac, clambda, obj_final, grad, rmat, &
                  v, iwvec, leniw, wvec, lenw)

      if (major_iters == 0) then
         write (lunlog,'(/, 2a)') progname, 'Terminating (0 iterations).'
         go to 999
      else
!        Repeat the optimal function evaluation, to wrap up this minimization.

         ncall = -4

         call fun (n, v, f, ncall, fail)

         if (fail) then
            write (lunlog,'(/, 2a)') &
               progname, 'Final objective function failed.'
            go to 999
         end if

      end if

      ncall = -5  ! First call outside NPOPT before each further minimization
      go to 10

999 continue

   close (lunlog)

   end program npopt_driver
