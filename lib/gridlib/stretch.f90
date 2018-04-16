!-Stretch--------------------------------------------------------------------
!
!  Routine to compute one-sided grid stretching on [0.0, 1.0].
!
!  Inputs:
!     1. nk    = number of points on the line
!     2. ds    = initial spacing (normalized to the unit interval)
!     3. beta  = initial guess for the stretching parameter (> 1.0);
!                try 1 + ds for small ds
!  Outputs:
!     1. beta  = value of the stretching parameter
!     2. s     = stretching function (ranges from 0 to 1)
!
!  Routines called: None
!
!  References:  Roberts, G. O. Computational Meshes for Boundary Layer
!               Problems, Proc. Second Int. Conf. Num. Methods Fluid Dyn.,
!               Lect. Notes Phys., vol. 8, Springer-Verlag, New York,
!               pp. 171-177 (1971).
!               Tannehill, J. C., Anderson, D. A. and Pletcher, R. H.,
!               Computational Fluid Mechanics and Heat Transfer,
!               second Edition, Taylor & Francis, Washington DC,
!               pp. 333-338 (1997).
!
!  Written by:  Dr. Dinesh K. Prabhu, Thermosciences Institute, NASA ARC
!
!  Date:        Nov. 18, 1996  Initial Fortran 90 version.
!               Mar. 22, 2004  Omit total length argument;
!                              pass normalized ds instead.
!----------------------------------------------------------------------------

SUBROUTINE Stretch (nk, ds, beta, s)

   IMPLICIT NONE

   !  Input arguments:

   INTEGER, INTENT(IN) :: nk
   REAL, INTENT(IN)    :: ds
   REAL, INTENT(INOUT) :: beta

   !  Output arguments:

   REAL, DIMENSION(nk), INTENT(INOUT) :: s

   !  Local constants:

   INTEGER, PARAMETER :: maxit = 25
   REAL, PARAMETER    :: one = 1.0, toler = 1.E-06

   !  Local variables:

   INTEGER :: iter, k
   REAL    :: dbeta, deta, eta, alfa, alfa2
   REAL    :: rbeta, betasq, term
   REAL    :: func, derv

   !  Execution:
   
   deta  = one/REAL(nk-1)
   alfa  = one-deta
   alfa2 = alfa+alfa
   dbeta = one

   DO iter = 1,maxit
!!!   write (6,2000) beta,dbeta
      rbeta  = (beta+one)/(beta-one)
      betasq = beta*beta-one
      term   = rbeta**alfa
      func   = ds+(beta*(term-one)/(term+one))-one
      derv   = ((rbeta**alfa2-one)- &
                (4.0*alfa*beta*term/betasq))/((term+one)**2)
      dbeta  = func/derv
      beta   = beta-dbeta
      IF (ABS(dbeta) < toler) EXIT
   END DO
   
!!!write (6,2000) beta,dbeta

   rbeta = (beta+one)/(beta-one)

   s(1) = 0.0
   DO k = 2,nk-1
      eta  = one-REAL(k-1)*deta
      term = rbeta**eta
      s(k) = one-beta*(term-one)/(term+one)
   END DO
   s(nk) = one

!!!2000 FORMAT(5X,'Stretch - beta, dbeta: ',1PE15.7,2X,'(', 1PE12.5,')')

END SUBROUTINE Stretch
