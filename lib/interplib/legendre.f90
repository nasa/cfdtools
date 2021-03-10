!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine legendre (fit, m, n, x, y, c, ssq, eval, neval, xeval, yeval, ier)

!     LEGENDRE calculates a least squares fit of a linear combination of the
!  Legendre polynomials P0(x):Pn(x) to an (x,y) dataset and/or evaluates the fit
!  at the given abscissas, or both.  The input dataset units may be arbitrary;
!  the abscissas are transformed to [-1, +1] here.  The data ordinates are NOT
!  scaled, so for the solution coefficients to be reasonably scaled, so should
!  those data ordinates.
!
!     This capability was prompted by a need for the 2-D analog of computing the
!  spherical harmonics for a 3-space surface (as a way of smoothing the surface)
!  and implemented as a possible learning exercise.  It supplements the many 1-D
!  interpolation and smoothing methods driven by the author's program SMOOTH.
!
!     The Legendre polynomials are defined by the following recurrence relation:
!
!       P (x) = 1      P (x) = x      nP (x) = (2n-1) x P   (x) - (n-1)P   (x)
!        0              1               n                n-1            n-2
!
!  for n = 2, 3, 4, ...
!
!     The least squares fit DOES include a coefficient for P (x) to allow for a
!  vertical shift to improve the fit.                       0
!
!     The fit is performed via the orthogonal factorization of subroutine HDESOL
!  which expects the right-hand-side vector to be entered as an extra column of
!  the left-hand-side matrix.
!
!  History:
!
!     06/30/2020  D.A.Saunders  Initial implementation.
!
!  Author:  David Saunders, AMA, Inc. at NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   logical, intent (in)  :: fit      ! T => calculate coefficients c(0:n);
                                     ! F => the coefficients are input

   integer, intent (in)  :: m        ! Number of data points; m >= n + 1

   integer, intent (in)  :: n        ! Defines number of polynomial terms 0:n

   real,    intent (in)  :: x(m), &  ! Data abscissas and ordinates, arbitrary
                            y(m)     ! units; scaled/shifted here prior to fit

   real,    intent (inout) :: c(0:n) ! Coefficients of best-fit polynomials for
                                     ! scaled/shifted data; ensuing evaluations
                                     ! work with abscissas similarly scaled and
                                     ! shifted here

   real,    intent (out) :: ssq      ! Sum of squared residuals (data units)

   logical, intent (in)  :: eval     ! T => evaluate the fit at xeval(:);
                                     ! F => don't perform any evaluations

   integer, intent (in)  :: neval    ! If (eval), perform neval evaluations

   real,    intent (in)  :: xeval(neval)  ! Abscissas at which to evaluate fit,
                                          ! scaled/shifted here as for x(:)

   real,    intent (out) :: yeval(neval)  ! Evaluated ordinates (if (eval))

   integer, intent (out) :: ier      ! 0 => No error detected;
                                     ! 1 => m < n + 1;
                                     ! 2 => n appears too big (n > 100)
!  Procedures:

   external  hdesol                  ! Linear least least squares via QR
                                     ! factorizn. (Householder transformations)
!  Local constants:

   integer, parameter :: nmax = 20   ! Reasonable limit on n
   real,    parameter :: one  = 1.
   real,    parameter :: u = -1., v = 1.  ! Working interval for Legendre polys.

!  Local variables:

   integer           :: j, jn, jnm1, jnm2, jt
   real              :: rmsdev
   real, save        :: xscale, xshift ! May be needed on a second "eval" call
   real              :: xmin, xmax     ! Data abscissa range
   real, allocatable :: Ab(:,:)        ! Augmented matrix A | b  for Ax ~ b
   real, allocatable :: P(:,:)         ! For polynomials j, j+1, j+2
   real, allocatable :: w(:)           ! Work-space for HDESOL
   real, allocatable :: xnorm(:)       ! For transforming x(:) to [-1, +1]

!  Execution:

   ier = 0

   if (fit) then

      if (m < n + 1) then
         write (*, *) ' *** Legendre error: m, n = ', m, n
         ier = 1
         go to 99
      end if

      if (n > nmax) then
         write (*, *) ' *** Likely Legendre error: n = ', n
         ier = 2
         go to 99
      end if

      xmin = minval (x(:))
      xmax = maxval (x(:))

      call getxform (xmin, xmax, u, v, xscale, xshift)

      allocate (xnorm(m), Ab(m,0:n+1), w(m), P(m,0:2))

      xnorm(:) = xscale*x(:) + xshift

      jnm2 = 0  ! Rotating indices for the recurrence relation
      jnm1 = 1
      jn   = 2

      P(:,0)  = one;  P(:,1)  = xnorm(:)
      Ab(:,0) = one;  Ab(:,1) = xnorm(:)

      do j = 2, n
         P(:,jn) = (real(2*j-1)*xnorm(:)*P(:,jnm1) - real(j-1)*P(:,jnm2))/ &
                    real(j)
         Ab(:,j) = P(:,jn)
         jt   = jnm2
         jnm2 = jnm1
         jnm1 = jn
         jn   = jt
      end do

      deallocate (xnorm)

      Ab(:,n+1) = y(:)  ! RHS vector

      call hdesol (m, m, n+2, Ab, w, ssq);  deallocate (Ab, P)

      c(0:n) = w(1:n+1);                    deallocate (w)
      rmsdev = sqrt (ssq/real(m))

   end if

   if (eval) then

      allocate (xnorm(neval), P(neval,0:2))
      xnorm(:) = xscale*xeval(:) + xshift

      jnm2 = 0  ! Rotating indices for the recurrence relation
      jnm1 = 1
      jn   = 2

      P(:,0) = one;  P(:,1) = xnorm(:)
      yeval(:) = c(0)*P(:,0) + c(1)*P(:,1)

      do j = 2, n
         P(:,jn) = (real(2*j-1)*xnorm(:)*P(:,jnm1) - real(j-1)*P(:,jnm2))/ &
                    real(j)
         yeval(:) = yeval(:)  + c(j)*P(:,jn)
         jt   = jnm2
         jnm2 = jnm1
         jnm1 = jn
         jn   = jt
      end do

      deallocate (xnorm)

   end if

99 return

   end subroutine legendre
