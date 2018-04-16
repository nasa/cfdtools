      subroutine dtcnt1 (n, cntour, aux, tol, mc, maxbsp, maxint, split,
     *                   bloks, rhs, seq, bsval, xknot, fval, fjac,
     *                   pval,
     *                   work, nwork, c, ier)


c----- This is the routine which actually computes the solution

      integer n, nwork, ier, ndeg, nord, ix, iy, iz, numint, ncoll, il
      integer nsize, maxint, maxbsp, mc, nvar, nknot, nrow, ir, ic
      integer ncoef, colskp, colinc, sbs, ebs, mdiag, mdim
      double precision tol, work(*), c(*), gauss(0:3)
      double precision bloks(*), seq(*), bsval(maxbsp,*), aux(*)
      double precision error, xknot(*), fval(*), fjac(*), pval(n,*)
      double precision tau, rhs(*), split(*), rcond, scale
      double precision dnrm2, ddot
      external cntour

c----- Determine the initial knot set for the solution

      nord = c(3)
      ndeg = nord - 1
      ncoll = ndeg - 1

c----- Determine the breakpoint sequence

      nsize = c(4)
      numint = 1
      seq(1) = c(6)
      do 10 ix = 1,nsize
        if (c(6+ix) .ne. seq(numint)) then
          numint = numint + 1
          seq(numint) = c(6+ix)
          endif
 10     continue
      numint = numint - 1

c----- Specify the collocation points

      gauss(3) = 1.0d0
      gauss(2) = sqrt (3.0d0) / 3.0d0
      gauss(1) = -gauss(2)
      gauss(0) = -1.0d0

c----- Clean up the initial guess

      nsize = ncoll * numint + 2
      nknot = nord + nsize
      do 50 ix = 1,n*nsize
        c(5+nknot+ix) = max (min (c(5+nknot+ix), 1.0d0), 0.0d0)
 50     continue

c----- Initialize for linear system solve

 100  continue
      il = 1
      nsize = ncoll * numint
      ncoef = nsize + 2
      nvar = n * nsize
      nknot = nord + ncoef
      mdiag = min (ndeg * n - 1, nvar - 1)
      mdim = 3 * mdiag + 1
      iz = nvar * mdim
      call dcopy (iz, 0.0d0, 0, bloks, 1)

c----- Establish the banded matrix and right-hand side

      nrow = 0
      do 250 ix = 1,numint
        sbs = 0
        if (ix .eq. 1) sbs = 1
        colskp = (2 * ix - 3 + sbs) * n
        colinc = n * (mdim - 1)
        ebs = nord - sbs
        if (ix .eq. numint) ebs = ebs - 1
        sbs = (ix - 1) * ncoll + 1 + sbs
        do 240 iy = 1,ncoll
          tau = seq(ix) + (seq(ix+1) - seq(ix)) *
     *                    (gauss(iy) + 1.0d0) * 0.5d0
          call dtspdr (tau, 2, c, work, nwork, pval, n, ier)
          if (ier .lt. 0) then
            ier = -200
            return
            endif
          call dtbspl (c(6), nknot, tau, il, nord, 2, maxbsp,
     *                 work, nwork, bsval, ier)
          if (ier .lt. 0) then
            ier = -200
            return
            endif
          call cntour (pval, aux, fval, fjac, ier)
          if (ier .ne. 0) then
            ier = -100 + ier
            return
            endif
          nrow = nrow + 1
          call dcopy (n-1, fval, 1, rhs(nrow), 1)
          do 220 ir = 1,n-1
            do 210 ic = 1,n
              iz = (colskp + ic - 1) * mdim + nrow + ir - colskp - ic +
     *             2 * mdiag
              call daxpy (ebs, fjac((ic-1)*(n-1)+ir),
     *                    bsval(sbs,1), 1, bloks(iz), colinc)
 210          continue
 220        continue
          nrow = nrow + n - 1
          rhs(nrow) = ddot (n, pval(1,2), 1, pval(1,3), 1)
          do 230 ic = 1,n
            iz = (colskp + ic - 1) * mdim + nrow + 1 - colskp - ic +
     *           2 * mdiag
            call daxpy (ebs, pval(ic,3), bsval(sbs,2), 1,
     *                  bloks(iz), colinc)
            call daxpy (ebs, pval(ic,2), bsval(sbs,3), 1,
     *                  bloks(iz), colinc)
 230        continue
 240      continue
 250    continue

c----- Solve the linear system

      call D0GBLE (bloks, mdim, nvar, mdiag, mdiag, rhs, work, nwork,
     *             rcond, ier)
      ier = min (ier, 0)
      if (ier .eq. -4) ier = 1
      if (ier .lt. 0) ier = -200
      if (ier .ne. 0) return

c----- Restrict solution to the appropriate domain

      do 320 ix = 1,nsize
        do 310 iy = 1,n
          if (c(6+nknot+(iy-1)*ncoef+ix) .lt. rhs((ix-1)*n+iy))
     *        rhs((ix-1)*n+iy) = c(6+nknot+(iy-1)*ncoef+ix)
          if (c(6+nknot+(iy-1)*ncoef+ix) - rhs((ix-1)*n+iy) .gt. 1.0d0)
     *        rhs((ix-1)*n+iy) = c(6+nknot+(iy-1)*ncoef+ix) - 1.0d0
 310      continue
 320    continue

c----- Update the solution

      do 330 ix = 1,n
        call daxpy (nsize, -1.0d0, rhs(ix), n,
     *              c(7+nknot+(ix-1)*ncoef), 1)
 330    continue

c----- See if another Newton iteration is required

      error = dnrm2 (n * nsize, rhs, 1)
      error = error / sqrt (dble (n * nsize))
      if (error .gt. tol) go to 100

c----- Determine which intervals need to be split

      do 430 ix = 1,numint
        error = 0.0d0
        do 420 iy = 0,ncoll
          tau = seq(ix) + (seq(ix+1) - seq(ix)) * 0.5d0 *
     *                    (0.5d0 * (gauss(iy) + gauss(iy+1)) + 1.0d0)
          call dtspvl (tau, c, work, nwork, pval, ier)
          if (ier .lt. 0) then
            ier = -200
            return
            endif
          call cntour (pval, aux, fval, fjac, ier)
          if (ier .ne. 0) then
            ier = -100 + ier
            return
            endif
          do 410 iz = 1,n-1
            scale = dnrm2 (n, fjac(iz), n-1)
            if (scale .ne. 0.0d0) then
              error = max (error, abs (fval(iz)) / scale)
            else
              error = max (error, abs (fval(iz)))
              endif
 410        continue
 420      continue
        split(ix) = error
 430    continue

c----- Split the appropriate intervals

      nsize = 0
      do 530 ix = numint,1,-1
        if (split(ix) .ge. tol) then
          iz = (split(ix) / tol) ** (1.0d0 / nord)
          if (numint + iz .gt. maxint) then
            ier = 2
            return
            endif
          do 510 iy = numint+iz+1,ix+iz+1,-1
            seq(iy) = seq(iy-iz)
 510        continue
          do 520 iy = 1,iz
            seq(ix+iy) = ((iz - iy + 1) * seq(ix) + iy * seq(ix+iz+1)) /
     *                   (iz + 1)
            nsize = nsize + 1
            xknot(nsize) = seq(ix+iy)
            nsize = nsize + 1
            xknot(nsize) = seq(ix+iy)
 520        continue
          numint = numint + iz
          endif
 530    continue

c----- If no knots need to be added, then quit

      if (nsize .eq. 0) then
        ier = 0
        return
        endif

c----- Add the knots and refit the curve

      call dcopy (nknot, c(6), 1, xknot(nsize+1), 1)
      nknot = nsize + nknot
      call dtsort (xknot, nknot, ier)
      if (ier .lt. 0) then
        ier = -200
        return
        endif
      il = mc
      call dtoslo (c, 1, xknot, nknot, il, work, nwork, bloks, ier)
      if (ier .lt. 0) then
        ier = -200
        return
        endif
      call dcopy (il, bloks, 1, c, 1)
      go to 100
      end

