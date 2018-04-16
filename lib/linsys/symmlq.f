C+----------------------------------------------------------------------

      subroutine symmlq( n, b, r1, r2, v, w, x, y,
     $                   Aprod, Msolve, goodb, precon, shift,
     $                   nout , ncheck, itnlim, rtol,
     $                   istop, itn, Anorm, Acond, rnorm, ynorm )
 
      external           Aprod, Msolve
      integer            n, nout, ncheck, itnlim, istop, itn
      logical            goodb, precon
      real               shift, rtol, Anorm, Acond, rnorm, ynorm,
     $                   b(n), r1(n), r2(n), v(n), w(n), x(n), y(n)
*     ------------------------------------------------------------------
*
*     symmlq  is designed to solve the system of linear equations
*
*                Ax = b
*
*     where A is an n by n symmetric matrix and b is a given vector.
*     The matrix A is not required to be positive definite.
*     (If A is known to be definite, the method of conjugate gradients
*     might be preferred, since it will require about the same number of
*     iterations as symmlq but slightly less work per iteration.)
*
*
*     The matrix A is intended to be large and sparse.  It is accessed
*     by means of a subroutine call of the form
*
*                call Aprod ( n, x, y )
*
*     which must return the product y = Ax for any given vector x.
*
*
*     More generally, symmlq is designed to solve the system
*
*                (A - shift*I) x = b
*
*     where  shift  is a specified scalar value.  If  shift  and  b
*     are suitably chosen, the computed vector x may approximate an
*     (unnormalized) eigenvector of A, as in the methods of
*     inverse iteration and/or Rayleigh-quotient iteration.
*     Again, the matrix (A - shift*I) need not be positive definite.
*     The work per iteration is very slightly less if  shift = 0.
*
*
*     A further option is that of preconditioning, which may reduce
*     the number of iterations required.  If M = C C' is a positive
*     definite matrix that is known to approximate  (A - shift*I)
*     in some sense, and if systems of the form  My = x  can be
*     solved efficiently, the parameters precon and Msolve may be
*     used (see below).  When  precon = .true., symmlq will
*     implicitly solve the system of equations
*
*             P (A - shift*I) P' xbar  =  P b,
*
*     i.e.                  Abar xbar  =  bbar
*     where                         P  =  C**(-1),
*                                Abar  =  P (A - shift*I) P',
*                                bbar  =  P b,
*
*     and return the solution       x  =  P' xbar.
*     The associated residual is rbar  =  bbar - Abar xbar
*                                      =  P (b - (A - shift*I)x)
*                                      =  P r.
*
*     In the discussion below, eps refers to the machine precision.
*     eps is computed by symmlq.  A typical value is eps = 2.22e-16
*     for IBM mainframes and PCs using double-precision arithmetic.
*
*     Parameters
*     ----------
*
*     n       input      The dimension of the matrix A.
*
*     b(n)    input      The rhs vector b.
*
*     r1(n)   workspace
*     r2(n)   workspace
*     v(n)    workspace
*     w(n)    workspace
*
*     x(n)    output     Returns the computed solution  x.
*
*     y(n)    workspace
*
*     Aprod   external   A subroutine defining the matrix A.
*                        For a given vector x, the statement
*
*                              call Aprod ( n, x, y )
*
*                        must return the product y = Ax
*                        without altering the vector x.
*
*     Msolve  external   An optional subroutine defining a
*                        preconditioning matrix M, which should
*                        approximate (A - shift*I) in some sense.
*                        M must be positive definite.
*                        For a given vector x, the statement
*
*                              call Msolve( n, x, y )
*
*                        must solve the linear system My = x
*                        without altering the vector x.
*
*                        In general, M should be chosen so that Abar has
*                        clustered eigenvalues.  For example,
*                        if A is positive definite, Abar would ideally
*                        be close to a multiple of I.
*                        If A or A - shift*I is indefinite, Abar might
*                        be close to a multiple of diag( I  -I ).
*
*                        NOTE.  The program calling symmlq must declare
*                        Aprod and Msolve to be external.
*
*     goodb   input      Usually, goodb should be .false.
*                        If x is expected to contain a large multiple of
*                        b (as in Rayleigh-quotient iteration),
*                        better precision may result if goodb = .true.
*                        See Lewis (1977) below.
*                        When goodb = .true., an extra call to Msolve
*                        is required.
*
*     precon  input      If precon = .true., preconditioning will
*                        be invoked.  Otherwise, subroutine Msolve
*                        will not be referenced; in this case the
*                        actual parameter corresponding to Msolve may
*                        be the same as that corresponding to Aprod.
*
*     shift   input      Should be zero if the system Ax = b is to be
*                        solved.  Otherwise, it could be an
*                        approximation to an eigenvalue of A, such as
*                        the Rayleigh quotient b'Ab / (b'b)
*                        corresponding to the vector b.
*                        If b is sufficiently like an eigenvector
*                        corresponding to an eigenvalue near shift,
*                        then the computed x may have very large
*                        components.  When normalized, x may be
*                        closer to an eigenvector than b.
*
*     nout    input      A file number.
*                        If nout .gt. 0, a summary of the iterations
*                        will be printed on unit nout.
*
*     ncheck  input      If ncheck .gt. 0, an extra call of Aprod will be
*                        used to check if A is symmetric.  Also,
*                        if precon = .true., an extra call of Msolve
*                        will be used to check if M is symmetric.
*
*     itnlim  input      An upper limit on the number of iterations.
*
*     rtol    input      A user-specified tolerance.  symmlq terminates
*                        if it appears that norm(rbar) is smaller than
*                              rtol * norm(Abar) * norm(xbar),
*                        where rbar is the transformed residual vector,
*                              rbar = bbar - Abar xbar.
*
*                        If shift = 0 and precon = .false., symmlq
*                        terminates if norm(b - A*x) is smaller than
*                              rtol * norm(A) * norm(x).
*
*     istop   output     An integer giving the reason for termination...
*
*              -1        The matrix Abar appears to be a multiple of
*                        the identity.  The solution is a multiple of b.
*
*               0        b = 0, so the exact solution is x = 0.
*                        No iterations were performed.
*
*               1        Norm(rbar) appears to be less than
*                        the value  rtol * norm(Abar) * norm(xbar).
*                        The solution in  x  should be acceptable.
*
*               2        Norm(rbar) appears to be less than
*                        the value  eps * norm(Abar) * norm(xbar).
*                        This means that the residual is as small as
*                        seems reasonable on this machine.
*
*               3        Norm(Abar) * norm(xbar) exceeds norm(b)/eps,
*                        which should indicate that x has essentially
*                        converged to an eigenvector of A
*                        corresponding to the eigenvalue shift.
*
*               4        Acond (see below) has exceeded 0.1/eps, so
*                        the matrix Abar must be very ill-conditioned.
*                        x may not contain an acceptable solution.
*
*               5        The iteration limit was reached before any of
*                        the previous criteria were satisfied.
*
*               6        The matrix defined by Aprod does not appear
*                        to be symmetric.
*                        For certain vectors y = Av and r = Ay, the
*                        products y'y and r'v differ significantly.
*
*               7        The matrix defined by Msolve does not appear
*                        to be symmetric.
*                        For vectors satisfying My = v and Mr = y, the
*                        products y'y and r'v differ significantly.
*
*               8        An inner product of the form  x' M**(-1) x
*                        was not positive, so the preconditioning matrix
*                        M does not appear to be positive definite.
*
*                        If istop .ge. 5, the final x may not be an
*                        acceptable solution.
*
*     itn     output     The number of iterations performed.
*
*     Anorm   output     An estimate of the norm of the matrix operator
*                        Abar = P (A - shift*I) P,   where P = C**(-1).
*
*     Acond   output     An estimate of the condition of Abar above.
*                        This will usually be a substantial
*                        under-estimate of the true condition.
*
*     rnorm   output     An estimate of the norm of the final
*                        transformed residual vector,
*                           P (b  -  (A - shift*I) x).
*
*     ynorm   output     An estimate of the norm of the final
*                        transformed solution vector y = P x.
*
*
*
*     To change precision
*     -------------------
*
*     Globally change the words
*            real
*            saxpy, scopy, sdot, snrm2
*     to their equivalents, and change the definitions of the parameter
*     constants for zero, one, ... below.
*
*     ------------------------------------------------------------------
*
*
*     This routine is an implementation of the algorithm described in
*     the following references:
*
*     C.C. Paige and M.A. Saunders,  Solution of Sparse Indefinite
*          Systems of Linear Equations,
*          SIAM J. Numer. Anal. 12, 4, September 1975, pp. 617-629.
*
*     J.G. Lewis,  Algorithms for Sparse Matrix Eigenvalue Problems,
*          Report STAN-CS-77-595, Computer Science Department,
*          Stanford University, Stanford, California, March 1977.
*
*     Applications of symmlq and the theory of preconditioning
*     are described in the following references:
*
*     D.B. Szyld and O.B. Widlund,  Applications of Conjugate Gradient
*          Type Methods to Eigenvalue Calculations,
*          in R. Vichnevetsky and R.S. Steplman (editors),
*          Advances in Computer Methods for Partial Differential
*          Equations -- III, IMACS, 1979, 167-173.
*
*     D.B. Szyld,  A Two-level Iterative Method for Large Sparse
*          Generalized Eigenvalue Calculations,
*          Ph. D. dissertation, Department of Mathematics,
*          New York University, New York, October 1983.
*
*     P.E. Gill, W. Murray, D.B. Ponceleon and M.A. Saunders,
*          Preconditioners for indefinite systems arising in
*          optimization, Report SOL 90-8, Dept of Operations Research,
*          Stanford University, Stanford, CA, 1990.
*     ------------------------------------------------------------------
*
*
*     SYMMLQ development:
*            1972: First version.
*            1975: John Lewis recommended modifications to help with
*                  inverse iteration:
*                  1. Reorthogonalize v1 and v2.
*                  2. Regard the solution as x = x1  +  bstep * b,
*                     with x1 and bstep accumulated separately
*                     and bstep * b added at the end.
*                     (In inverse iteration, b might be close to the
*                     required x already, so x1 may be a lot smaller
*                     than the multiple of b.)
*            1978: Daniel Szyld and Olof Widlund implemented the first
*                  form of preconditioning.
*                  This required both a solve and a multiply with M.
*            1979: Implemented present method for preconditioning.
*                  This requires only a solve with M.
*            1984: Sven Hammarling noted corrections to tnorm and x1lq.
*                  symmlq added to NAG Fortran Library.
*     15 Sep 1985: Final F66 version.  symmlq sent to "misc" in netlib.
*     16 Feb 1989: First F77 version.
*
*     22 Feb 1989: Hans Mittelmann observed beta2 = 0 (hence failure)
*                  if Abar = const*I.  istop = -1 added for this case.
*
*     01 Mar 1989: Hans Mittelmann observed premature termination on
*                  ( 1  1  1 )     (   )                   ( 1  1    )
*                  ( 1  1    ) x = ( 1 ),  for which  T3 = ( 1  1  1 ).
*                  ( 1     1 )     (   )                   (    1  1 )
*                  T2 is exactly singular, so estimating cond(A) from
*                  the diagonals of Lbar is unsafe.  We now use
*                  L       or  Lbar         depending on whether
*                  lqnorm  or  cgnorm       is least.
*
*     03 Mar 1989: eps computed internally instead of coming in as a
*                  parameter.
*     07 Jun 1989: ncheck added as a parameter to say if A and M
*                  should be checked for symmetry.
*     20 Nov 1990: goodb added as a parameter to make Lewis's changes
*                  an option.  Usually b is NOT much like x.  Setting
*                  goodb = .false. saves a call to Msolve at the end.
*     20 Nov 1990: Residual not computed exactly at end, to save time
*                  when only one or two iterations are required
*                  (e.g. if the preconditioner is very good).
*                  Beware, if precon is true, rnorm estimates the
*                  residual of the preconditioned system, not Ax = b.
*
*
*     Michael A. Saunders            (na.saunders @ NA-net.stanford.edu)
*     Department of Operations Research                   (415) 723-1875
*     Stanford University
*     Stanford, CA 94305-4022
*     ------------------------------------------------------------------
*
*
*     Subroutines and functions
*
*     USER       Aprod, Msolve
*     BLAS       saxpy, scopy, sdot, snrm2
*
*
C-----------------------------------------------------------------------

*     Intrinsics and local variables
 
      intrinsic          abs, max, min, mod, sqrt
      real               sdot, snrm2
      real               alfa, b1, beta, beta1, bstep, cs,
     $                   cgnorm, dbar, delta, denom, diag,
     $                   eps, epsa, epsln, epsr, epsx,
     $                   gamma, gbar, gmax, gmin, gpert,
     $                   lqnorm, oldb, qrnorm, rhs1, rhs2,
     $                   s, sn, snprod, t, tnorm,
     $                   x1cg, x1lq, ynorm2, zbar, z
      integer            i
 
      real               zero         ,  one         ,  two
      parameter        ( zero = 0.0e+0,  one = 1.0e+0,  two = 2.0e+0 )
 
      character*16       enter, exit
      character*60       msg(-1:8)
 
      data               enter /' Enter symmlq.  '/,
     $                   exit  /' Exit  symmlq.  '/
 
      data               msg
     $ / 'The matrix Abar is a multiple of the identity',
     $   'The exact solution is  x = 0',
     $   'Requested accuracy achieved, as determined by rtol',
     $   'Reasonable accuracy achieved, given eps',
     $   'x has converged to an eigenvector',
     $   'Acond has exceeded 0.1/eps',
     $   'The iteration limit was reached',
     $   'Aprod  does not define a symmetric matrix',
     $   'Msolve does not define a symmetric matrix',
     $   'Msolve does not define a positive definite preconditioner' /
*     ------------------------------------------------------------------
 
 
*     Compute eps, the machine precision.  The call to saxpy is
*     intended to fool compilers that use extra-length registers.
 
      eps    = one
 
   10 eps    = eps / two
      x(1)   = eps
      y(1)   = one
      call saxpy ( 1, one, x, 1, y, 1 )
      if (y(1) .gt. one) go to 10
 
      eps    = eps * two
 
*     Print heading and initialize.
 
      if (nout .gt. 0) write(nout, 1000)
     $    enter, n, goodb, precon, itnlim, rtol, shift
      istop  = 0
      itn    = 0
      Anorm  = zero
      Acond  = zero
      rnorm  = zero
      ynorm  = zero
 
      do 50 i = 1, n
         x(i) = zero
   50 continue
 
*     Set up y for the first Lanczos vector v1.
*     y is really beta1 * P * v1  where  P = C**(-1).
*     y and beta1 will be zero if b = 0.
 
      call scopy ( n, b, 1, y , 1 )
      call scopy ( n, b, 1, r1, 1 )
      if ( precon ) call Msolve( n, r1, y )
      if ( goodb  ) then
         b1  = y(1)
      else
         b1  = zero
      end if
      beta1  = sdot  ( n, r1, 1, y, 1 )
 
*     See if Msolve is symmetric.
 
      if (ncheck .gt. 0  .and.  precon) then
         call Msolve( n, y, r2 )
         s      = sdot  ( n, y, 1, y, 1 )
         t      = sdot  ( n,r1, 1,r2, 1 )
         z      = abs( s - t )
         epsa   = (s + eps) * eps**0.33333
         if (z .gt. epsa) then
            istop = 7
            go to 900
         end if
      end if
 
*     Test for an indefinite preconditioner.
 
      if (beta1 .lt. zero) then
         istop = 8
         go to 900
      end if
 
*     If b = 0 exactly, stop with x = 0.
 
      if (beta1 .eq. zero) then
         go to 900
      end if
 
*     Here and later, v is really P * (the Lanczos v).
 
      beta1  = sqrt( beta1 )
      s      = one / beta1
      do 100 i = 1, n
         v(i)  = s * y(i)
  100 continue
 
*     See if Aprod  is symmetric.
 
      call Aprod ( n, v, y )
      if (ncheck .gt. 0) then
         call Aprod ( n, y, r2 )
         s      = sdot  ( n, y, 1, y, 1 )
         t      = sdot  ( n, v, 1,r2, 1 )
         z      = abs( s - t )
         epsa   = (s + eps) * eps**0.33333
         if (z .gt. epsa) then
            istop = 6
            go to 900
         end if
      end if
 
*     Set up y for the second Lanczos vector.
*     Again, y is beta * P * v2  where  P = C**(-1).
*     y and beta will be zero or very small if Abar = I or constant * I.
 
      call saxpy ( n, (- shift), v, 1, y, 1 )
      alfa   = sdot  ( n, v, 1, y, 1 )
      call saxpy ( n, (- alfa / beta1), r1, 1, y, 1 )
 
*     Make sure  r2  will be orthogonal to the first  v.
 
      z      = sdot  ( n, v, 1, y, 1 )
      s      = sdot  ( n, v, 1, v, 1 )
      call saxpy ( n, (- z / s), v, 1, y, 1 )
 
      call scopy ( n, y, 1, r2, 1 )
      if ( precon ) call Msolve( n, r2, y )
      oldb   = beta1
      beta   = sdot  ( n, r2, 1, y, 1 )
 
      if (beta .lt. zero) then
         istop = 8
         go to 900
      end if
 
*     Cause termination (later) if beta is essentially zero.
 
      beta   = sqrt( beta )
      if (beta .le. eps) then
         istop = -1
      end if
 
*     See if the local reorthogonalization achieved anything.
 
      denom  = sqrt( s ) * snrm2( n, r2, 1 )  +  eps
      s      = z / denom
      t      = sdot  ( n, v, 1, r2, 1 ) / denom
      if (nout .gt. 0  .and.  goodb) then
         write(nout, 1100) beta1, alfa, s, t
      end if
 
*     Initialize other quantities.
 
      cgnorm = beta1
      gbar   = alfa
      dbar   = beta
      rhs1   = beta1
      rhs2   = zero
      bstep  = zero
      snprod = one
      tnorm  = alfa**2
      ynorm2 = zero
      gmax   = abs( alfa )
      gmin   = gmax
 
      if ( goodb ) then
         do 200 i = 1, n
            w(i)  = zero
  200    continue
      else
         call scopy ( n, v, 1, w, 1 )
      end if
 
*     ------------------------------------------------------------------
*     Main iteration loop.
*     ------------------------------------------------------------------
 
*     Estimate various norms and test for convergence.
 
  300 Anorm  = sqrt( tnorm  )
      ynorm  = sqrt( ynorm2 )
      epsa   = Anorm * eps
      epsx   = Anorm * ynorm * eps
      epsr   = Anorm * ynorm * rtol
      diag   = gbar
      if (diag .eq. zero) diag = epsa
 
      lqnorm = sqrt( rhs1**2 + rhs2**2 )
      qrnorm = snprod * beta1
      cgnorm = qrnorm * beta / abs( diag )
 
*     Estimate  cond(A).
*     In this version we look at the diagonals of  L  in the
*     factorization of the tridiagonal matrix,  T = L*Q.
*     Sometimes, T(k) can be misleadingly ill-conditioned when
*     T(k+1) is not, so we must be careful not to overestimate Acond.
 
      if (lqnorm .le. cgnorm) then
         Acond  = gmax / gmin
      else
         denom  = min( gmin, abs( diag ) )
         Acond  = gmax / denom
      end if
 
*     See if any of the stopping criteria are satisfied.
*     In rare cases, istop is already -1 from above (Abar = const * I).
 
      if (istop .eq. 0) then
         if (itn    .ge. itnlim ) istop = 5
         if (Acond  .ge. 0.1/eps) istop = 4
         if (epsx   .ge. beta1  ) istop = 3
         if (cgnorm .le. epsx   ) istop = 2
         if (cgnorm .le. epsr   ) istop = 1
      end if
*     ==================================================================
 
*     See if it is time to print something.
 
      if (nout .le.  0)          go to 600
      if (n    .le. 40)          go to 400
      if (itn  .le. 10)          go to 400
      if (itn  .ge. itnlim - 10) go to 400
      if (mod(itn,10)  .eq.   0) go to 400
      if (cgnorm .le. 10.0*epsx) go to 400
      if (cgnorm .le. 10.0*epsr) go to 400
      if (Acond  .ge. 0.01/eps ) go to 400
      if (istop  .ne. 0)         go to 400
      go to 600
 
*     Print a line for this iteration.
 
  400 zbar   = rhs1 / diag
      z      = (snprod * zbar  +  bstep) / beta1
      x1lq   = x(1)  +  b1 * bstep / beta1
      x1cg   = x(1)  +  w(1) * zbar  +  b1 * z
 
      if (    itn     .eq. 0) write(nout, 1200)
      write(nout, 1300) itn, x1cg, cgnorm, bstep/beta1, Anorm, Acond
      if (mod(itn,10) .eq. 0) write(nout, 1500)
*     ==================================================================
 
 
*     Obtain the current Lanczos vector  v = (1 / beta)*y
*     and set up  y  for the next iteration.
 
  600 if (istop .ne. 0) go to 800
      s      = one / beta
 
      do 620 i = 1, n
         v(i)  = s * y(i)
  620 continue
 
      call Aprod ( n, v, y )
      call saxpy ( n, (- shift), v, 1, y, 1 )
      call saxpy ( n, (- beta / oldb), r1, 1, y, 1 )
      alfa   = sdot( n, v, 1, y, 1 )
      tnorm  = tnorm  +  alfa**2  +  two * beta**2
      call saxpy ( n, (- alfa / beta), r2, 1, y, 1 )
      call scopy ( n, r2, 1, r1, 1 )
      call scopy ( n, y, 1, r2, 1 )
      if ( precon ) call Msolve( n, r2, y )
      oldb   = beta
      beta   = sdot  ( n, r2, 1, y, 1 )
      if (beta .lt. zero) then
         istop = 6
         go to 800
      end if
      beta   = sqrt( beta )
 
*     Compute the next plane rotation for  Q.
 
      gamma  = sqrt( gbar**2 + oldb**2 )
      cs     = gbar / gamma
      sn     = oldb / gamma
      delta  = cs * dbar  +  sn * alfa
      gbar   = sn * dbar  -  cs * alfa
      epsln  = sn * beta
      dbar   =            -  cs * beta
 
*     Update  x.
 
      z      = rhs1 / gamma
      s      = z * cs
      t      = z * sn
 
      do 700 i = 1, n
         x(i)  = (w(i) * s   +   v(i) * t)  +  x(i)
         w(i)  =  w(i) * sn  -   v(i) * cs
  700 continue
 
*     Accumulate the step along the direction  b,
*     and go round again.
 
      bstep  = snprod * cs * z  +  bstep
      snprod = snprod * sn
      gmax   = max( gmax, gamma )
      gmin   = min( gmin, gamma )
      ynorm2 = z**2  +  ynorm2
      rhs1   = rhs2  -  delta * z
      rhs2   =       -  epsln * z
      itn    = itn   +  1
      go to 300
 
*     ------------------------------------------------------------------
*     End of main iteration loop.
*     ------------------------------------------------------------------
 
*     Move to the CG point if it seems better.
*     In this version of SYMMLQ, the convergence tests involve
*     only cgnorm, so we're unlikely to stop at an LQ point,
*     EXCEPT if the iteration limit interferes.
 
  800 if (cgnorm .le. lqnorm) then
         zbar   = rhs1 / diag
         bstep  = snprod * zbar  +  bstep
         ynorm  = sqrt( ynorm2  +  zbar**2 )
         rnorm  = cgnorm
         call saxpy ( n, zbar, w, 1, x, 1 )
      else
         rnorm  = lqnorm
      end if
 
      if ( goodb ) then
 
*        Add the step along  b.
 
         bstep  = bstep / beta1
         call scopy ( n, b, 1, y, 1 )
         if ( precon ) call Msolve( n, b, y )
         call saxpy ( n, bstep, y, 1, x, 1 )
      end if
 
*     ==================================================================
*     Display final status.
*     ==================================================================
  900 if (nout  .gt. 0) then
         write(nout, 2000) exit, istop, itn,
     $                     exit, Anorm, Acond,
     $                     exit, rnorm, ynorm
         write(nout, 3000) exit, msg(istop)
      end if
 
      return
 
*     ------------------------------------------------------------------
 1000 format(// 1p,    a, 5x, 'Solution of symmetric   Ax = b'
     $ / ' n      =', i7, 5x, 'goodb  =',   l4, 12x, 'precon =', l4
     $ / ' itnlim =', i7, 5x, 'rtol   =', e11.2, 5x, 'shift  =', e23.14)
 1100 format(/ 1p, ' beta1  =', e10.2, 3x, 'alpha1 =', e10.2
     $       / ' (v1,v2) before and after ', e14.2
     $       / ' local reorthogonalization', e14.2)
 1200 format(// 5x, 'itn', 7x, 'x1(cg)',
     $   10x, 'norm(r)', 5x, 'bstep', 7x, 'norm(A)', 3X, 'cond(A)')
 1300 format(1p, i8, e19.10, e11.2, e14.5, 2e10.2)
 1500 format(1x)
 2000 format(/ 1p, a, 6x, 'istop =', i3,   15x, 'itn   =', I8
     $       /     a, 6x, 'Anorm =', e12.4, 6x, 'Acond =', e12.4
     $       /     a, 6x, 'rnorm =', e12.4, 6x, 'ynorm =', e12.4)
 3000 format(      a, 6x, a )
*     ------------------------------------------------------------------
*     end of symmlq
      end
