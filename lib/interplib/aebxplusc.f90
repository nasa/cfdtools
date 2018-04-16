!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine AeBxplusC (mode, n, x, y, A, B, C, ssqmin, ier)
!
!  Purpose:
! 
!     For an (x,y) dataset, calculate the nonlinear least squares fit of an
!     exponential-type function as follows:
!
!             Mode 1:  y = Ae^(B*x) + C   ! x >= 0
!
!             Mode 2:  y = Ae^(B/x) + C   ! x >  0
!
!     A > 0 is assumed in both cases, and the constant term may be suppressed
!     from the curve fit in each case by entering C = 0.
!
!     See related subroutine AxBeCx if an x^B term is desired as well.
!
!  Method (mode 1):
!
!     A starting guess for the nonlinear coefficient B is provided by the method
!     of Geoffrey Rowe outlined in the Software Tech Briefs supplement to NASA
!     Tech Briefs, September 2010.  Then the A and C coefficients are calculated
!     with a linear least squares method that avoids the so-called "normal"
!     equations.  The published article highlights avoidance of an iteration as
!     a benefit of the method for estimating B, specially for real-time applica-
!     tions where a predictable computational time for a constant number of data
!     points is desirable.  However, any estimate of B is not necessarily the
!     best in some sense.  An iteration is still desirable to optimize B along
!     with A and C.  Since the non-derivate 1-D minimizer employed here performs
!     a virtually constant number of iterations, the computational time remains
!     quite consistent (and efficient) for a given dataset size.  Key idea:
!
!             slope = dy/dx = ABe^(Bx)
!         =>  slope1/slope2 = e^B(x1 - x2)
!         =>              B = ln (slope1/slope2) / (x1 - x2)
!
!     The two slope estimates are obtained from straight line fits within a pair
!     of data windows taken from the "left" portion of the dataset, with some
!     overlap. (Since they serve only to provide a reasonable search interval
!     for B, the particular choice is not critical.)
!
!     There is no requirement for the data abscissas to be uniformly spaced, but
!     they are assumed to increase monotonically and to be non-negative.
!
!  Method (mode 2):
!
!     Analogously, the estimate for B can be shown to be:
!
!             B = [ln (slope1/slope2) - 2 ln (x2/x1)] / (1/x1 - 1/x2)
!  History:
!
!     11/10/2010  D.A.Saunders  Initial design (mode 1).
!     11/12/2010    "    "      Initial implementation (mode 1).
!     11/15/2010    "    "      Added mode 2 option.
!     11/16/2010    "    "      Suppress calculation of C if C = 0. on input.
!
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
   real,    intent (out)   :: A, B       ! Coefs. of fitted y = Ae^(Bx) + C
   real,    intent (inout) :: C          ! C = 0. on input means exclude it
   real,    intent (out)   :: ssqmin     ! Minimized sum of squared deviations
   integer, intent (inout) :: ier        ! Input 0 suppresses iteration printing
                                         ! Output 0 => no error detected;
                                         !       -1 => 1-D minimizer fatal error
                                         !       -2 => too few data pts or x < 0
!  Constants:

   real,      parameter :: one    = 1.
   real,      parameter :: zero   = 0.
   real,      parameter :: scale1 = 5.   ! Applied to the estimate of B to give
   real,      parameter :: scale2 = 0.2  ! the search interval for optimizing it
   integer,   parameter :: nfmax  = 50   ! Limit on # function evaluations
   character, parameter :: caller * 9 = 'AeBxplusC'

!  Variables:

   integer :: istat, lunerr, lunout, numfun
   logical :: excludeC
   real    :: Bhi, Blo, tol

!  Procedures:

   external :: fminrc  ! Reverse-communication 1-D minimizer
   external :: hdesol  ! Linear least squares via orthogonal factorization

!  Execution:

!  Check that the inputs are reasonable:

   if (n < 3 .or. x(1) < zero) then
      ier = -2
      go to 99
   end if

   if (mode == 2) then
      if (x(1) == zero) then
         ier = -2
         go to 99
      end if
   end if

   excludeC = C == zero

   lunout = -6
   if (ier /= 0) lunout = -lunout

!  Calculate an estimate for B.  Handling small numbers of data points is the
!  main awkwardness for the heuristics in choosing a couple of windows.

   call estimateB ()

   Blo = scale1 * B  ! Search interval for optimized B; B is negative
   Bhi = scale2 * B

   if (mode == 2) then  ! Estimate is less robust
      Blo = Blo * 2.0
      Bhi = Bhi * 0.5
   end if

!  A minimum can be found to within sqrt (machine epsilon), but avoid the sqrt:

   if (epsilon (tol) < 1.e-10) then
      tol = 1.e-8
   else
      tol = 1.e-4
   end if

   lunerr = abs (lunout)
   numfun = nfmax      ! Limit; FMINRC typically takes 12-14 iterations
   istat  = 2          ! Initialize the minimization
   ier    = 0          ! Upon exit, normally

   do while (istat > 0)

      call fminrc (Blo, Bhi, B, ssqmin, tol, numfun, caller, lunout, istat)

      if (istat < -1) then

         write (lunerr, '(/, 2a)') caller, ': FMINRC fatal error'
         ier = -1

      else if (istat < 0) then ! Iteration limit; may be usable

         write (lunerr, '(/, 2a)') caller, ': Iteration limit.'
         istat = 0

      else  ! istat >= 0       ! Evaluate the objective function at B

         call computeAC ()     ! Linear least squares soln. for A, C for this B;
                               ! the sum of squared deviations is a by-product
         if (istat == 0) exit  ! Success
      end if

   end do

99 return

   contains

!     Internal procedures for subroutine AeBxplusC:

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine computeAC ()  ! Linear least squares calculation for given B

!     Local variables:

      integer :: i, ncp1
      real    :: A1(n,3), coefs(n)

!     Execution:

      if (excludeC) then
         ncp1 = 2  ! Only one coefficient to calculate
      else
         ncp1 = 3
      end if

!     Set up the overdetermined system with the RHS as column nc + 1:

      if (mode == 1) then
         if (excludeC) then
            do i = 1, n
               A1(i,1) = exp (B * x(i))
               A1(i,2) = y(i)
            end do
         else
            do i = 1, n
               A1(i,1) = exp (B * x(i))
               A1(i,2) = one
               A1(i,3) = y(i)
            end do
         end if
      else ! mode = 2
         if (excludeC) then
            do i = 1, n
               A1(i,1) = exp (B / x(i))
               A1(i,2) = y(i)
            end do
         else
            do i = 1, n
               A1(i,1) = exp (B / x(i))
               A1(i,2) = one
               A1(i,3) = y(i)
            end do
         end if
      end if

      call hdesol (n, n, ncp1, A1, coefs, ssqmin) ! Orthogonal factorzn. & soln.

      A = coefs(1)
      if (.not. excludeC) C = coefs(2)

      end subroutine computeAC

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine estimateB ()

      ! Define a pair of data windows towards at the lower end (higher slopes).
      ! Fit straight lines to the points in the two windows.
      ! Estimate B from the slopes of those lines and the window separation.

!     Local constants:

      real,    parameter :: half = 0.5, two = 2.
      real,    parameter :: wfraction = 0.2  ! Window size as fraction of n > 10

!     Local variables:

      integer :: i1w2, i2w2, nw
      real    :: c1, c2, logterm, slope1, slope2, x1, x2

!     Execution:

      select case (n)  ! n = # data points

         case (3)  ! Minimum to define two distinct windows

            nw   = 2  ! Window size
            i1w2 = 2
            i2w2 = 3

         case (4, 5)

            nw   = 3
            i1w2 = 2
            i2w2 = 4

         case (6:9)

            nw   = 3
            i1w2 = 3
            i2w2 = 5

         case default

            nw   = nint (real (n) * wfraction) + 2
            i1w2 = nint (real (nw) * half)
            i2w2 = i1w2 + nw - 1

      end select

      if (lunout > 0) &
         write (lunout, '(a, 4i6)') 'n,nw,i1w2,i2w2:', n, nw, i1w2, i2w2

      call linefit (nw, x(1),    y(1),    slope1, c1)
      call linefit (nw, x(i1w2), y(i1w2), slope2, c2)

      logterm = log (slope1 / slope2)
      x1 = xmean (nw, x(1))
      x2 = xmean (nw, x(i1w2))

      if (mode == 1) then
         B = logterm / (x1 - x2)
      else ! mode = 2: the estimate is not as robust
         B = (logterm - two * log (x2 / x1)) / (one/x1 - one/x2)
         if (B > zero) B = -B
      end if

      if (lunout > 0) write (lunout, '(a, 1p, 5e17.8)') &
         'x1,x2,s1,s2,B:', x1, x2, slope1, slope2, B

      end subroutine estimateB

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine linefit (npts, x, y, slope, yintercept)

!     Fit a straight line to (x,y) data via orthogonal decomposition.
!     Avoid the more general pnfit (polynomial fit) for this specialized need.

!     Arguments:

      integer, intent (in)  :: npts  ! # data points
      real,    intent (in)  :: x(npts), y(npts)    ! Data coordinates
      real,    intent (out) :: slope, yintercept   ! y ~ slope * x + yintercept

!     Local variables:

      integer :: i
      real    :: A1(npts,3), coefs(npts), ssq

!     Execution:

!     Set up the overdetermined system with the RHS as column 3:

      do i = 1, npts
         A1(i,1) = x(i)
         A1(i,2) = one
         A1(i,3) = y(i)
      end do

      call hdesol (npts, npts, 3, A1, coefs, ssq)  ! Orthogonal factzn. & soln.

      slope      = coefs(1)
      yintercept = coefs(2)

      end subroutine linefit

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real function xmean (n, x)     ! Arithmetic mean of n values

      integer, intent (in)  :: n     ! # points to average
      real,    intent (in)  :: x(n)  ! Values to average

      xmean = sum (x) / real (n)

      end function xmean

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   end subroutine AeBxplusC
