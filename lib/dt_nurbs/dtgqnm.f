      subroutine dtgqnm (n, func, aux, xint, maxit, tol, prob, ldsol,
     *                   ncol, work, nwork, num, soln, ier)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc     This subroutine can be used to find all zeros of a nonlinear
cc     system of equations defined over a rectangular region.  The
cc     subroutine uses a quasi-Newton method to find the solutions.
cc     The calling sequence for the routine is the following:
cc
cc     call dtgqnm (n, func, aux, xint, maxit, tol, prob, ldsol, ncol,
cc                  work, nwork, num, soln, ier)
cc
cc     where the arguments have the following significance:
cc
cc     n     - the number of equations and the number of variables
cc     func  - the name of a user supplied subroutine which evaluates
cc             the nonlinear functions.  The calling sequence for the
cc             subroutine is
cc
cc             call func (n, x, aux, f, ier)
cc
cc             where
cc
cc             n   - the size of the system
cc             x   - the vector to evaluate at
cc             aux - the auxilliary array
cc             f   - the value of the function at x (output)
cc             ier - the error / success code (output)
cc                   a negative value will abort dtgqnm
cc
cc     aux   - the auxilliary array which is passed to func.  This array
cc             is not altered or touched by dtgqnm
cc     xint  - a 2 x n array containing the upper and lower bounds for
cc             the variables
cc     maxit - the maximum number of iterations allowed for each zero
cc     tol   - the tolerance to which the zeros should be found
cc     prob  - the probability of finding all of the zeros.  The closer
cc             this value is to 1, the greater the amount of computer
cc             time required.  Values larger than 0.95 should be used
cc             with caution.
cc     ldsol - the leading declared dimension of array soln.
cc     ncol  - the number of columns declared in the soln array.
cc     work  - a working storage array of length nwork
cc     nwork - the length of the working storage array.  This number
cc             should be at least 2 n ** 2 + 8 n
cc     num   - on output, the number of zeros found
cc     soln  - an ldsol x ncol array containing the solutions on output
cc     ier   - the success / error code for this program
cc             ier = 0; success (results computed)
cc             ier = 1; ncol zeros found before completion; increase the
cc                      size of ncol and try again
cc             ier = -1; n <= 0
cc             ier = -2; xint(1,i) >= xint(2,i) for some i
cc             ier = -3; maxit <= 0
cc             ier = -4; tol <= 0
cc             ier = -5; prob <= 0 or prob >= 1
cc             ier = -6; ldsol < n
cc             ier = -7; ncol <= 0
cc             ier = -8; nwork < 2 n ** 2 + 8 N
cc             ier = -100 + ifail; ifail returned by func
cc             ier = -200; unexpected error from lower level
cc
cc     Thomas Grandine
cc     August, 1990
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c----- Declare the parameters for this subprogram

      integer n, maxit, ldsol, ncol, nwork, num, ier
      double precision aux(*), xint(2,*), tol, work(*)
      double precision soln(ldsol,*), prob

c----- Determine clustering value

      double precision clustr
      parameter (clustr = 0.5d0)

c----- Declare the variables that will be needed in this case

      integer ia, isave, iguess, ires, iresn, iupd, idiff, isol, iwork
      integer nleft, info, iter, maxtry, ix, iy, tries, need
      double precision expect, resid, residn, rcond, scale, alfa, dist
      double precision check, sofar, sftol, ptol, update, skip, skipc
      real rand
      character*8 subnam

c----- Declare functions which are needed

      double precision ddot

c----- Declare the user function to be external

      external func

c----- Initialize subnam

      data subnam /'DTGQNM'/

c----- Check the input for errors

      if (n .le. 0) then
        ier = -1
        num = -1
        call dterr (1, subnam, ier, 0)
        return
        endif
      do 10 ix = 1,n
        if (xint(1,ix) .ge. xint(2,ix)) then
          ier = -2
          num = -1
          call dterr (1, subnam, ier, 0)
          return
          endif
 10     continue
      if (maxit .le. 0) then
        ier = -3
        num = -1
        call dterr (1, subnam, ier, 0)
        return
        endif
      if (tol .le. 0) then
        ier = -4
        num = -1
        call dterr (1, subnam, ier, 0)
        return
        endif
      if ((prob .le. 0.0d0) .or. (prob .ge. 1.0d0)) then
        ier = -5
        num = -1
        call dterr (1, subnam, ier, 0)
        return
        endif
      if (ldsol .lt. n) then
        ier = -6
        num = -1
        call dterr (1, subnam, ier, 0)
        return
        endif
      if (ncol .le. 0) then
        ier = -7
        num = -1
        call dterr (1, subnam, ier, 0)
        return
        endif
      need = 2 * n ** 2 + 8 * n
      if (nwork .lt. need) then
        ier = -8
        num = -1
        call dterr (2, subnam, ier, need)
        return
        endif

c----- Initialize random number generator and set up expected value

      call dtrpst (nwork)
      expect = (3.0d0 * prob - 1.0d0) / (2.0d0 * (1.0d0 - prob))
      maxtry = 2.0d0 * expect + 3.5d0
      skipc = log (clustr) / (maxtry - 1)
      sofar = 0.0d0
      ptol = sqrt (tol)

c----- Chop up the working storage into useful pieces

      ia = 1
      isave = ia + n ** 2
      iguess = isave + n ** 2
      ires = iguess + n
      iresn = ires + n
      iupd = iresn + n
      idiff = iupd + n
      isol = idiff + n
      iwork = isol + n
      nleft = nwork - iwork + 1

c----- Set up number of zeros found so far

      num = 0
      tries = 0
      check = - sqrt (tol)

c----- Start another attempt at solution

 100  continue
      tries = tries + 1
      if (tries .gt. maxtry) return
      iter = 0

c----- Set up the initial guess

      do 110 ix = 1,n
        call dtrpun (1, rand)
        work(isol+ix-1) = xint(1,ix) + rand * (xint(2,ix) - xint(1,ix))
 110    continue

c----- Modify guess away from the boundary

 120  continue
      call daxpy (n, check, 1.0d0, 0, work(isol), 1)
      call dtgqn1 (n, xint, work(isol))
      call func (n, work(isol), aux, work(ires), info)
      if (info .lt. 0) then
        ier = -100 + info
        call dterr (3, subnam, ier, 0)
        return
        endif
      resid = sqrt (ddot (n, work(ires), 1, work(ires), 1))
      sofar = max (sofar, resid)
      skip = max (exp (skipc * (tries - 1)), clustr)
      if (resid .gt. skip * sofar) go to 100

c----- Determine initial guess for Jacobian

      do 130 ix = 1,n
        work(isol+ix-1) = work(isol+ix-1) - check
        call func (n, work(isol), aux, work(iresn), info)
        if (info .lt. 0) then
          ier = -100 + info
          call dterr (3, subnam, ier, 0)
          return
          endif
        call daxpy (n, -1.0d0, work(iresn), 1, work(ires), 1)
        call dscal (n, 1.0d0 / check, work(ires), 1)
        call dcopy (n, work(ires), 1, work(isave+(ix-1)*n), 1)
        call dcopy (n, work(iresn), 1, work(ires), 1)
 130    continue
      call dcopy (n**2, work(isave), 1, work(ia), 1)

c----- Evaluate the function at the initial point

      resid = sqrt (ddot (n, work(ires), 1, work(ires), 1))

c----- Branch back point to start the main iteration

 200  continue
      sofar = max (sofar, resid)
      iter = iter + 1
      if (iter .gt. maxit) go to 100

c----- Solve for the update direction

      call dcopy (n, work(ires), 1, work(iupd), 1)
      call dterpt (1)
      call dtgele (work(ia), n, n, work(iupd), work(iwork), nleft,
     *             rcond, info)
      call dterpt (1)
      if (info .eq. -4) go to 120
      if (info .lt. 0) then
        ier = -200
        call dterr (5, subnam, ier, 0)
        return
        endif

c----- Determine maximum stepsize

      alfa = -1.0d0
      do 210 ix = 1,n
        if (work(iupd+ix-1) .eq. 0.0d0) go to 210
        if (work(iupd+ix-1) .lt. 0.0d0) then
          dist = xint(2,ix) - work(isol+ix-1)
        else
          dist = xint(1,ix) - work(isol+ix-1)
          endif
        if (abs (dist) .ge. tol) alfa = max (alfa,
     *                                       dist / work(iupd+ix-1))
 210    continue

c----- Throw points away if too far out of bounds

      if (alfa .gt. -0.5d0) go to 100

c----- Determine the update

 220  continue
      call dcopy (n, work(isol), 1, work(iguess), 1)
      call daxpy (n, alfa, work(iupd), 1, work(iguess), 1)
      call dtgqn1 (n, xint, work(iguess))

c----- Check for reduction in residual

      call func (n, work(iguess), aux, work(iresn), info)
      if (info .lt. 0) then
        ier = -100 + info
        call dterr (3, subnam, ier, 0)
        return
        endif
      residn = sqrt (ddot (n, work(iresn), 1, work(iresn), 1))
      if (residn .gt. resid) then
        alfa = 0.5d0 * alfa
        if (alfa .gt. check) go to 120
        go to 220
        endif

c----- Take this new step

      call dcopy (n, work(iguess), 1, work(isol), 1)
      update = sqrt (alfa ** 2 * ddot (n, work(iupd), 1, work(iupd), 1))
     *         * rcond

c----- Check to see if we are done

      sftol = max (1.0d0, sofar)
      if ((residn .lt. sftol * tol) .and. (update .lt. tol)) go to 300

c----- Determine difference in gradients

      call dcopy (n, work(iresn), 1, work(idiff), 1)
      call daxpy (n, -1.0d0, work(ires), 1, work(idiff), 1)

c----- Determine product of Jacobian and step

      do 230 ix = 1,n
        work(idiff+ix-1) = work(idiff+ix-1) - alfa *
     *                     ddot (n, work(isave+ix-1), n, work(iupd), 1)
 230    continue

c----- Perform the quasi-Newton update step

      scale = alfa * ddot (n, work(iupd), 1, work(iupd), 1)
      do 250 ix = 1,n
        do 240 iy = 1,n
          work(isave+n*iy+ix-n-1) = work(isave+n*iy+ix-n-1) +
     *             work(idiff+ix-1) * work(iupd+iy-1) / scale
 240      continue
 250    continue

c----- Copy over results and run again

      call dcopy (n**2, work(isave), 1, work(ia), 1)
      call dcopy (n, work(iresn), 1, work(ires), 1)
      resid = residn
      go to 200

c----- See if this is a new zero or not

 300  continue
      do 310 ix = 1,num
        call dcopy (n, work(isol), 1, work(idiff), 1)
        call daxpy (n, -1.0d0, soln(1,ix), 1, work(idiff), 1)
        scale = sqrt (ddot (n, work(idiff), 1, work(idiff), 1))
        if (scale .lt. ptol) go to 100
 310    continue

c----- New root found - record it

      num = num + 1
      call dcopy (n, work(isol), 1, soln(1,num), 1)
      maxtry = (expect * (num + 1.0d0) + (expect + 1.0d0)) *
     *         (num + 1.0d0) + 2.5d0
      if (num .lt. ncol) go to 100

c----- We have found the maximum number of roots

      ier = 1
      call dterr (0, subnam, ier, 0)

c----- This is the end of the program

      end

c----- This subroutine projects points back onto the boundary

      subroutine dtgqn1 (n, xint, guess)

      double precision xint(2,*), guess(*)
      integer n, ix

      do 10 ix = 1,n
        guess(ix) = max (xint(1,ix), min (xint(2,ix), guess(ix)))
 10     continue
      end

