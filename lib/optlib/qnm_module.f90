!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This is a replacement for the original common block /qnmctl/ originally
!  employed by QNMDIF2.

   module qnm_module

      integer, parameter :: maxvar = 20  ! Adjust this limit on the number of
                                         ! optimization variables to suit the
                                         ! application

      integer :: igrad       ! 0 => QNMDIF2 decides when to use forward and
                             !      central finite differences;
                             ! 1 => QNMDIF2 should NOT switch from forward to
                             !      central differences
      integer :: lunlog      ! For output from QNMDIF2
      integer :: nallow      ! Max. # function evaluations allowed
      integer :: nftotl      ! # function evaluations, incremented by QNMDIF2
                             ! (and QNM_SOLUTION for CENDIF2 calls)
      integer :: niter       ! Limit on the # steps taken by QNMDIF2

      real    :: epsmch      ! Machine precision (EPSILON intrinsic)
      real    :: epsobj      ! Estimate of the function precision
      real    :: eta         ! In [0,1]; 1. => loose line searches; 0.1 => tight
      real    :: lowbound    ! Some lower bound on the objective at the minimum
      real    :: stepmx      ! Either some known bound on the distance to the
                             ! minimum, or a very large number for no effect
      real    :: tol         ! Often sqrt (epsobj), or larger; QNMDIF2 is done
                             ! when ||gradient|| < tol**(1/3) & last step < tol
      real    :: zeta        ! Factor by which h(:) are reduced if FUN fails
      real    :: h(maxvar)   ! Finite differencing intervals for the (scaled)
                             ! variables x(:).  Setting any h(i) < 0 means that
                             ! CENDIF2 will estimate h(i).  If likely values are
                             ! known, it is recommended that they be input in
                             ! units of the physical variables.  CENDIF2 and
                             ! QNMDIF2 will set h(i) = physical h(i)/xscale(i)
                             ! where  h(i)*xscale(i) = physical interval (i).
      real :: xscale(maxvar) ! Should be used to make each variable x(i) seen
                             ! by QNMDIF2 to be O(1), roughly (not critical).
                             ! x(i)*xscale(i) + xshift(i) = physical var. i;
                             ! x(i) + (physical var.(i) - xshift(i))/xscale(i)
      real :: xshift(maxvar) ! See xscale(:)

      logical :: contin      ! T tells QNMDIF2 to update the gradient and the
                             ! Hessian approximation if max. # iters. is reached
      logical :: local       ! F => convergence to a saddle point is not likely;
                             ! T => QNMDIF2 does a local search after satisfying
                             ! the convergence criteria, in case a saddle point
                             ! has been arrived at rather than a true minimum
      logical :: print       ! F suppresses most printing from QNMDIF2
      logical :: rescue      ! F is the norm.  See QNMDIF2 for how to start with
                             ! a local search for a difficult problem.
   end module qnm_module
