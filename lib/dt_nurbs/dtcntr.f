      subroutine dtcntr (cntour, aux, tol, mc, work, nwork, c, ier)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc     DTCNTR follows contours using spline collocation.  It assumes
cc     an initial guess which is a cubic spline with double knots in
cc     the interior which is accurate at the endpoints.  It also
cc     assumes that the parametric domain is [0,1] and that there are
cc     n >= 2 dependent variables. Usage is as follows:
cc
cc     double precision aux(naux), tol, work(nwork), c(mc)
cc     call dtcntr (cntour, aux, tol, mc, work, nwork, c, ier)
cc
cc     where
cc
cc     cntour - a user callable function which evaluates the contouring
cc              function.  It should be of the form
cc
cc              subroutine cntour (t, aux, f, jac, ier)
cc              double precision t(n), aux(naux), f(n-1),
cc             +                 jac(n-1,n)
cc
cc              where
cc
cc              t   - the parameter values to evaluate at
cc              aux - an auxiliary array which gets passed on
cc              f   - on output, the vector of function values
cc              jac - on output, the n-1 x n Jacobian matrix
cc              ier - error flag; values other than 0 are returned
cc                    to the upper level program
cc
cc     aux    - an auxiliary array which gets passed on to cntour
cc     tol    - the tolerance to which the problem should be solved
cc     work   - the working storage array
cc     nwork  - the size of work array.  It should be >= 9n*mc
cc     c      - on input, the initial guess at the contour
cc              on output, the final solution
cc     ier    - success / error flag
cc              ier = 1;  contour is singular at solution
cc              ier = 2;  mc not large enough to satisfy tol
cc              ier = 0;  success - no errors encountered
cc              ier = -1;  c(3) <> 4
cc              ier = -3;  nwork too small
cc              ier = -6;  c(4) < 4
cc              ier = -8;  invalid knot sequence
cc              ier = -51;  c(1) <> 1
cc              ier = -52;  c(2) < 2
cc              ier = -61;  mc too small to get started
cc              ier = -150;  tol <= 0
cc              ier = -100+i;  error i returned from cntour
cc              ier = -200;  unexpected error from lower level
cc
cc     Thomas Grandine
cc     April, 1992
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer n, mc, nwork, ier
      integer need, maxint, maxbsp, isplit, iseq, ndeg, nbs, i, mode
      integer ibsval, ixknot, irhs, iwork, iblok, ifval, ijac, ipval
      double precision aux(*), tol, work(*), c(*)
      character*8 subnam
      parameter (subnam = 'DTCNTR')
      external cntour

      n = c(2)
c----- Determine input errors

      if (int(c(1)) .ne. 1) then
        ier = -51
        go to 9000
        endif
      if (n .lt. 2) then
        ier = -52
        go to 9000
        endif
      if (int(c(3)) .ne. 4) then
        ier = -1
        go to 9000
        endif
      need = 9 * n * mc
      if (nwork .lt. need) then
        ier = -3
        go to 9000
        endif
      nbs = c(4)
      if (nbs .lt. 4) then
        ier = -6
        go to 9000
        endif

c----- To check knots, verify 4 knots at 0 and 4 at one,
c      be sure there are an even number of them, and that
c      internal knots are all double.

      do 10 i = 1, 4
        if ( c(5+i) .ne. 0.d0  .or.  c(5+nbs+i) .ne. 1.d0 ) then
          ier = -8
          go to 9000
          endif
   10   continue
      if ( mod(nbs,2) .ne. 0 ) then
        ier = -8
        go to 9000
        endif

c----- This loop starts at 8 rather than 10 to ensure that
c      the first double knot is .gt. 0

      do 20 i = 8, 4 + nbs, 2
        if ( c(i) .ne. c(i+1) .or. c(i) .ge. c(i+2) ) then
          ier = -8
          go to 9000
          endif
   20   continue

      if ( mc .lt. 9 + nbs*(n+1) ) then
        ier = -61
        go to 9000
        endif
      if (tol .le. 0.0d0) then
        ier = -150
        go to 9000
        endif

c----- Determine amount of room in spline definition array

      ndeg = 3
      maxbsp = (mc - 6 - ndeg) / (n + 1)
      maxint = maxbsp / (ndeg - 1)

c----- Chop up the working storage section for passing along

      isplit = 1
      iblok = isplit + maxint
      irhs = iblok + (n * maxbsp) * (9 * n - 2)
      iseq = irhs + n * maxbsp
      ibsval = iseq + maxint + 1
      ixknot = ibsval + 3 * maxbsp
      ifval = ixknot + maxbsp + ndeg + 1
      ijac = ifval + n - 1
      ipval = ijac + (n - 1) * n
      iwork = ipval + 3 * n
      need = nwork - iwork + 1

c----- Proceed to the next routine with all of the arrays chopped up

      ier = 0
      call dterpt (1)
      call dtcnt1 (n, cntour, aux, tol, mc, maxbsp, maxint,
     *             work(isplit), work(iblok), work(irhs), work(iseq),
     *             work(ibsval), work(ixknot), work(ifval), work(ijac),
     *             work(ipval),
     *             work(iwork), need, c, ier)
      call dterpt (1)

 9000 continue
      if (ier .gt. 0) then
        call dterr (0, subnam, ier, 0)
        endif
      if (ier .lt. 0) then      
        c(1) = -1.0d0
        mode = 1
        if (ier .eq. -3) mode = 2
        if (ier .eq. -200) mode = 3
        call dterr (mode, subnam, ier, need)
        endif
      end
