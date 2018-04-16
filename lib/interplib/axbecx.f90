!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine AxBeCx (mode, n, x, y, A, B, C, ssqmin, ier)
!
!  Purpose:
!
!     AxBeCx is a generalization of AeBxplusC, q.v., although the technique is
!     different and the constant term is omitted for starting guess reasons.
!
!     For an (x,y) dataset, calculate the nonlinear least squares fit of an
!     exponential-type function as follows:
!
!             Mode 1:  y = Ax^Be^(C*x)       ! x > 0
!
!             Mode 2:  y = Ax^Be^(C/x)       ! x > 0
!
!     A > 0 is assumed in both cases.  See subroutine AeBxplusC if the x^B term
!     is not wanted.
!
!  Method:
!
!     A starting guess for the coefficients could be obtained by calling the
!     related AeBxplusC with its input C = 0., giving A & C for when B = 0.
!     However, once B is treated as nonzero, the idea of using a pair of slope
!     estimates breaks down.  Instead, we omit the + D constant option and
!     compute estimates of A, B and C via linear least squares in log space.
!     Since the effect of deviations from the curve is distorted in log space,
!     we go further by optimizing further in real space, using a pair of nested
!     minimization iterations for the two nonlinear coefficients along with
!     updates of A, which can only be positive from the log-space estimate.
!
!     The robustness and efficiency of the FMIN-based 1-D minimizer is expected
!     to be preferable to use of a general-purpose n-dimensional optimizer.
!     Note that FMINRC2 is simply a copy of FMINRC to facilitate the nesting by
!     avoiding having to save and restore internal variables between calls.
!
!     There is no requirement for the data abscissas to be uniformly spaced, but
!     they are assumed to increase monotonically and to be non-negative.
!
!  History:
!
!     11/16/2010  D.A.Saunders  Initial dual-mode implementation of AeBxplusC.
!     11/17/2010    "    "      Log-space + nested iteration variant to include
!                               x^B term (but exclude the + D constant option).
!  Author:
!
!     David Saunders, ERC, Inc./NASA Ames Research Center, Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)    :: mode       ! 1 or 2 as outlined above
   integer, intent (in)    :: n          ! Number of data points
   real,    intent (in)    :: x(n), y(n) ! Data pt. coords.; x >= 0, monotonic
   real,    intent (out)   :: A, B, C    ! Coefs. of fitted y = Ax^Be^(C[/]x)
   real,    intent (out)   :: ssqmin     ! Minimized sum of squared deviations
   integer, intent (inout) :: ier        ! Input 0 suppresses iteration printing
                                         ! Output 0 => no error detected;
                                         !       -1 => 1-D minimizer fatal error
                                         !       -2 => too few data pts or x < 0
!  Constants:

   real,      parameter :: one    = 1.
   real,      parameter :: zero   = 0.
   real,      parameter :: scale1 = 2.   ! Applied to the B, C estimates to give
   real,      parameter :: scale2 = 0.5  ! the search intervals for reoptimizing
   integer,   parameter :: nfmax  = 50   ! Limit on # function evaluations
   character, parameter :: calleri * 7  = 'AxBeCxi'
   character, parameter :: callero * 7  = 'AxBeCxo'

!  Variables:

   integer :: istati, istato, lunerr, lunout, numfuni, numfuno
   real    :: Bhi, Blo, Chi, Clo, ssqi, ssqo, tol
   real    :: AA(n,2), coefso(n)

!  Procedures:

   external :: fminrc, fminrc2  ! Reverse-communication 1-D minimizer & a copy
   external :: hdesol           ! Linear least sqrs via orthogonal factorization

!  Execution:

!  Check that the inputs are reasonable:

   if (n < 4 .or. x(1) <= zero) then
      ier = -2
      go to 99
   end if

!  Calculate initial estimates in log space (linear least squares):

   call log_space_estimate ()

   lunout = -6
   if (ier /= 0) lunout = -lunout

   lunerr = abs (lunout)

   if (lunout > 0) then
      write (*, '(a, 1p, 4e17.8)') &
         'Log-space estimates for A, B, C, ssq:   ', A, B, C, ssqo
      call sumofsquares (ssqo)  ! In real space now
      write (*, '(a, 1p, e17.8)') &
         'Corresponding real space sum of squares:', ssqo
   end if

   Blo = scale1 * B  ! Search interval for optimized B; B is negative
   Bhi = scale2 * B

!  A minimum can be found to within sqrt (machine epsilon), but avoid the sqrt:

   if (epsilon (tol) < 1.e-10) then
      tol = 1.e-8
   else
      tol = 1.e-4
   end if

   numfuno = nfmax      ! Limit; FMINRC typically takes 12-14 iterations
   istato  = 2          ! Initialize the minimization
   ier     = 0          ! Upon exit, normally

   do while (istato > 0)

      call fminrc2 (Blo, Bhi, B, ssqo, tol, numfuno, callero, lunout, istato)

      if (istato < -1) then

         write (lunerr, '(/, 2a)') callero, ': fatal outer FMINRC error'
         ier = -1

      else if (istato < 0) then   ! Iteration limit; may be usable

         write (lunerr, '(/, 2a)') callero, ': Iteration limit.'
         istato = 0

      else  ! istato >= 0         ! Evaluate the objective function at B

         call computeC ()         ! Inner nonlinear C optimization for this A, B
         if (ier /= 0) go to 99

         call computeA ()         ! Linear least sqrs. soln. for A for this B, C

!!!      call sumofsquares (ssqo) ! The sum of squares is a by-product

         if (istato == 0) exit    ! Success
      end if

   end do

   ssqmin = ssqo

99 return

   contains

!     Internal procedures for subroutine AxBeCx:

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine log_space_estimate ()  ! Linear least squares estimate
                                        ! of log A, B and C
!     Local variables:

      integer :: i
      real    :: A1(n,4), coefs(n)

!     Execution:

!     Set up the overdetermined system for 3 coefs. with the RHS as column 4:

      if (mode == 1) then
         do i = 1, n
            A1(i,1) = one
            A1(i,2) = log (x(i))
            A1(i,3) = x(i)
            A1(i,4) = log (y(i))
         end do
      else ! mode = 2
         do i = 1, n
            A1(i,1) = one
            A1(i,2) = log (x(i))
            A1(i,3) = one /x(i)
            A1(i,4) = log (y(i))
         end do
      end if

      call hdesol (n, n, 4, A1, coefs, ssqo) ! Orthogonal factorzn. & soln.

      A = exp (coefs(1))
      B = coefs(2)
      C = coefs(3)

      end subroutine log_space_estimate

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine computeC ()

!     Inner nonlinear minimization, optimizing C for current values of A and B.

!     Local variables:

!     Execution:

      Clo     = scale1 * C
      Chi     = scale2 * C
      numfuni = nfmax      ! Limit
      istati  = 2          ! Initialize the minimization
      ier     = 0          ! Upon exit, normally

      do while (istati > 0)

         call fminrc (Clo, Chi, C, ssqi, tol, numfuni, calleri, lunout, istati)

         if (istati < -1) then

            write (lunerr, '(/, 2a)') calleri, ': fatal inner FMINRC error'
            ier = -1

         else if (istati < 0) then   ! Iteration limit; may be usable

            write (lunerr, '(/, 2a)') calleri, ': Iteration limit.'
            istati = 0

         else  ! istati >= 0         ! Evaluate the objective function at C

            call sumofsquares (ssqi) ! Sum of squared deviations being minimized
                                     ! with respect to C for current A and B
            if (istati == 0) exit    ! Success
         end if

      end do

      end subroutine computeC

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine computeA ()  ! Compute linear coefficient A for current B, C

!     Local variables:

      integer :: i

!     Execution:

!     Set up the overdetermined system with the RHS as column 2:

      if (mode == 1) then
         do i = 1, n
            AA(i,1) = x(i)**B * exp (C*x(i))
            AA(i,2) = y(i)
         end do
      else ! mode = 2
         do i = 1, n
            AA(i,1) = x(i)**B * exp (C/x(i))
            AA(i,2) = y(i)
         end do
      end if

      call hdesol (n, n, 2, AA, coefso, ssqo)  ! Orthogonal factzn. & soln.

      A = coefso(1)

      if (lunout > 0) write (*, '(a, 1p, e16.8)') 'A from llsq:', A

      end subroutine computeA

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine sumofsquares (ssq)  ! Sum of squared deviations

!     Arguments:

      real, intent (out) :: ssq

!     Local variables:

      integer :: i

!     Execution:

      ssq = zero
      if (mode == 1) then
         do i = 1, n
            ssq = (A * x(i)**B * exp (C*x(i)) - y(i))**2 + ssq
         end do
      else ! mode = 2
         do i = 1, n
            ssq = (A * x(i)**B * exp (C/x(i)) - y(i))**2 + ssq
         end do
      end if

      end subroutine sumofsquares

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   end subroutine AxBeCx
