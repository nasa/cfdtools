      subroutine dtpcut (c, plane, tol, mxpts, iopt, work, nwork,
     *                   ncuts, cuts, icptr, ier)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc      The subroutine DTPCUT is used to find all intersections of
cc      a given surface with a plane.  The surface may have any number
cc      of dependent variables.  In this context, "plane" means any
cc      hyperplane in a space whose dimension is the same as the
cc      number of dependent variables.  The calling sequence for this
cc      subroutine has the following form:
cc
cc      call dtpcut (c, plane, tol, mxpts, iopt, work, nwork,
cc     *             ncuts, cuts, icptr, ier)
cc
cc      where the parameters mean
cc
cc      c      - the input spline surface
cc      plane  - the plane with which to intersect the surface
cc      tol    - the tolerance to which curves are obtained
cc      mxpts  - the maximum number of discrete points on a cut branch
cc      iopt   - the output curve option flag
cc               iopt = 1    generate curves in parameter space
cc               iopt = 2    generate curves in range space
cc               iopt = 3    generate both sets of curves
cc      work   - a work array of length nwork
cc      nwork  - the length of the work array
cc      ncuts  - the number of curve components of the planar cut
cc      cuts   - the output spline vectors for each planar cut component
cc      icptr  - the index array which points to the start, inside
cc               cuts, of each of the spline vectors of the solution
cc      ier    - the error control flag
cc               ier = 0     no errors detected
cc               ier = 1     tolerance may not be achieved; increase
cc                           mxpts or tol
cc               ier = 2     input plane contains the surface
cc               ier = -1    curve(1) <> 2
cc               ier = -2    cin(2) = 0
cc               ier = -3    invalid spline order;
cc               ier = -4    invalid number of B-spline coefficients
cc               ier = -5    invalid knot sequence
cc               ier = -6    tol <= 0
cc               ier = -7    insufficient working storage
cc               ier = -8    iopt <> 1 and iopt <> 2 and iopt <> 3
cc               ier = -9    tol too large
cc               ier = -100  unexpected error return from lower level
cc
cc      Thomas Grandine
cc      August, 1989
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c----- Declare the necessary variables for this program

      double precision c(*), work(*), cuts(*), plane(*), tol
      integer nwork, ncuts, iopt, icptr(*), ier, iw1, iw2, iw3, mxpts
      integer mc, ndep, nknots, ncoef, nleft, iw4, iw5, iw6, need, nw
      character*8 subnam
      data subnam /'DTPCUT  '/

c----- Check the input for validity

      call dterpt (0)
      call dtschk (c, ier)
      if (ier .ne. 0) go to 9900
      if (c(1) .ne. 2) then
        ier = -1
        go to 9900
        endif
      if (tol .le. 0.0d0) then
        ier = -6
        go to 9900
        endif
      if ((iopt .ne. 1) .and. (iopt .ne. 2) .and. (iopt .ne. 3)) then
        ier = -8
        go to 9900
        endif

c----- Chop up the working storage into appropriate chunks

      ndep = c(2)
      if (ndep .lt. 0) ndep = -ndep - 1
      ncoef = c(5) * c(6)
      nknots = c(3) + c(4) + c(5) + c(6)
      mc = 8 + nknots + ncoef
      nw = 2.0d0 * (c(3) ** 2 + c(3) + c(4) ** 2 + c(4) + 1.0d0)
      need = mc + 2 * nw + 148 + (ndep + 4) * mxpts + mxpts ** 2 + ndep
      if (nwork .lt. need) then
        ier = -7
        go to 9900
        endif
      iw1 = 1
      iw2 = iw1 + mc + nw
      iw3 = iw2 + mxpts
      iw4 = iw3 + mxpts
      if (iopt .eq. 1) iw5 = iw4 + 2 * mxpts
      if (iopt .eq. 2) iw5 = iw4 + ndep * mxpts
      if (iopt .eq. 3) iw5 = iw4 + max (2, ndep) * mxpts
      iw6 = iw5 + ndep + 1
      nleft = nwork - iw6 + 1

c----- Now call the lower level subroutine to do the computations

      call dtpcu1 (c, plane, tol, iopt, mxpts, work(iw1),
     *             work(iw2), work(iw3), work(iw4), work(iw5),
     *             work(iw6), nleft, ncuts, cuts, icptr, ier)

c----- Now perform the error processing step

9900  continue
      call dterpt (1)
      if (ier .lt. 0) then
        if (ier .eq. -7) call dterr (2, subnam, ier, need)
        if (ier .eq. -100) call dterr (5, subnam, ier, 0)
        if ((ier .ne. -7) .and. (ier .ne. -100)) then
          call dterr (1, subnam, ier, 0)
          endif
        endif
      if (ier .gt. 0) call dterr (0, subnam, ier, 0)
      end
      
c----- This next subroutine is used to drive the planar cutter

      subroutine dtpcu1 (c, plane, tol, iopact, mxpts, cp,
     *                   xknots, tpts, upts, vals, work, nwork,
     *                   ncuts, cuts, icptr, ier)

c----- Declare the variables for this subprogram

      integer ier, ix, iy, nz, ncuts, ider(7), icptr(*), iprm(2)
      integer nv, nh, ioff, iv, ih, nc, nwork, ndep, nleft, mxpts, ic
      integer iopt, iopact, nfit, iw1, iw2, iw3, iw4, iwarn, n1, n2
      double precision c(*), work(*), spts(2,10), epts(2,10), plane(*)
      double precision tpts(*), upts(mxpts,*), plormn, cp(*), vals(*)
      double precision xknots(*), cuts(*), arclen, bnd(4), t, fu, fv
      double precision tol, dtmcon, ptol, uv(2), dist, fuu, fuv, fvv

c----- Define the indexing operations

      ider(1) = 1
      ider(2) = 0
      ider(3) = 1
      ider(4) = 1
      ider(5) = 2
      ider(6) = 0
      ider(7) = 2
      iprm(1) = 2
      iprm(2) = 1

c----- Determine the number of dependent variables and reset c(2)

      ndep = c(2)
      if (ndep .lt. 0) ndep = - ndep - 1

c----- Copy and modify the surface into a contour array

      nh = c(5)
      nv = c(6)
      nc = nh * nv
      ioff = 8 + c(3) + c(4) + nv + nh
      call dcopy (ioff, c, 1, cp, 1)
      cp(2) = 1
      call dcopy (nc, 0.0d0, 0, cp(1+ioff), 1)
      do 20 ix = 1,ndep
        do 10 iy = 1,nc
          cp(ioff+iy) = cp(ioff+iy) + plane(ix) * c(ioff+nc*(ix-1)+iy)
10        continue
20      continue
      do 30 iy = 1,nc
        if (c(2) .gt. 0.0d0) then
          cp(ioff+iy) = cp(ioff+iy) - plane(ndep+1)
        else
          cp(ioff+iy) = cp(ioff+iy) - plane(ndep+1) * c(ioff+nc*ndep+iy)
          endif
30      continue

c----- Determine whether or not the plane contains the surface

      do 40 iy = 1,nc
        if (abs (cp(ioff+iy)) .gt. tol) go to 50
40      continue
      ier = 2
      return
50    continue

c----- Determine whether or not a permutation is necessary

      iv = 0
      ih = 0
      do 110 ix = 2,nh
        if (cp(ioff+ix) * cp(ioff+ix-1) .le. 0.0d0) ih = ih + 1
        if (cp(ioff+nc-ix+1) * cp(ioff+nc-ix+2) .le. 0.0d0) ih = ih + 1
110     continue
      do 120 ix = 1,nv-1
        if (cp(ioff+ix*nh+1) * cp(ioff+(ix-1)*nh+1) .le. 0.0d0)
     *             iv = iv + 1
        if (cp(ioff+ix*nh) * cp(ioff+(ix+1)*nh) .le. 0.0d0) iv = iv + 1
120     continue
      fu = 0.0d0
      fuu = 0.0d0
      fv = 0.0d0
      fvv = 0.0d0
      do 130 ix = 1,nh
        fu = fu + abs (cp(ioff+ix))
        fuu = fuu + abs (cp(1+ioff+nc-ix))
130     continue
      if (fu .lt. tol) call dcopy (nh, 0.0d0, 0, cp(ioff+1), 1)
      if (fuu .lt. tol) call dcopy (nh, 0.0d0, 0, cp(ioff+nc-nh+1), 1)
      do 140 ix = 1,nv
        fv = fv + abs (cp(ioff+(ix-1)*nh+1))
        fvv = fvv + abs (cp(ioff+ix*nh))
140     continue
      if (fv .lt. tol) call dcopy (nv, 0.0d0, 0, cp(ioff+1), nh)
      if (fvv .lt. tol) call dcopy (nv, 0.0d0, 0, cp(ioff+nh), nh)

c----- Permute the parameters to minimize number of turning points

      if ((fu .lt. tol) .or. (fuu .lt. tol)) iv = ih
      if ((fv .lt. tol) .or. (fvv .lt. tol)) ih = iv + 1
      if (ih .gt. iv) then
        call dcopy (ioff+nc, cp, 1, work, 1)
        nleft = nwork - ioff - nc
        call dtsprm (work, iprm, work(ioff+nc+1), nleft, cp, ier)
        if (ier .lt. 0) then
          ier = -100
          return
          endif
        endif
      call dtsutl (cp, 2, bnd, ier)
      if (ier .lt. 0) then
        ier = -100
        return
        endif

c----- Scale the parametric region to be [0,1]x[0,1]

      n1 = cp(3) + cp(5)
      n2 = cp(4) + cp(6)
      do 150 ix = 1,n1
        cp(8+ix) = (cp(8+ix) - bnd(1)) / (bnd(2) - bnd(1))
150     continue
      do 160 ix = 1,n2
        cp(8+n1+ix) = (cp(8+n1+ix) - bnd(3)) / (bnd(4) - bnd(3))
160     continue

c----- Determine the topology of the contours

      call dtpcu4 (cp, work, nwork, ncuts, spts, epts, ier)
      if (ier .lt. 0) then
        ier = -100
        return
        endif

c----- Set up what has to happen in order to fit both curves

      iopt = iopact
      if (iopt .eq. 3) iopt = 1
      if (iopt .eq. 1) then 
        nfit = 2
      else
        nfit = ndep
        endif

c----- Chop up the workspace for tracing the contours

      ic = 1
      icptr(ic) = 1
      iw1 = 1
      iw2 = iw1 + mxpts * (mxpts - 2)
      iw3 = iw2 + mxpts
      iw4 = iw3 + mxpts
      nleft = nwork - iw4 + 1
      iwarn = 0
      do 220 ix = 1,ncuts

c----- Determine which direction to trace contour by checking the cases

        plormn = 1.0d0

c----- Case 1:  Start at left boundary

        if (spts(1,ix) .eq. 0.0d0) then
          call dtnpdr (spts(1,ix), 1, ider(2), cp, work, nwork, t, ier)
          if (ier .lt. 0) then
            ier = -100
            return
            endif
          if (t .lt. 0.0d0) plormn = -1.0d0
          endif

c----- Case 2:  Start at right boundary

        if (spts(1,ix) .eq. 1.0d0) then
          call dtnpdr (spts(1,ix), 1, ider(2), cp, work, nwork, t, ier)
          if (ier .lt. 0) then
            ier = -100
            return
            endif
          if (t .gt. 0.0d0) plormn = -1.0d0
          endif

c----- Case 3:  Start at bottom boundary

        if (spts(2,ix) .eq. 0.0d0) then
          call dtnpdr (spts(1,ix), 1, ider(1), cp, work, nwork, t, ier)
          if (ier .lt. 0) then
            ier = -100
            return
            endif
          if (t .gt. 0.0d0) plormn = -1.0d0
          endif

c----- Case 4:  Start at top boundary

        if (spts(2,ix) .eq. 1.0d0) then
          call dtnpdr (spts(1,ix), 1, ider(1), cp, work, nwork, t, ier)
          if (ier .lt. 0) then
            ier = -100
            return
            endif
          if (t .lt. 0.0d0) plormn = -1.0d0
          endif

c----- Set the initial start parameters

        uv(1) = spts(1,ix)
        uv(2) = spts(2,ix)

c----- Set the initial point for the contour trace

        t = 0.0d0
        nz = 0
        call dtpcu7 (iopt, iv, ih, nz, ndep, upts, mxpts, uv, bnd, c,
     *               plane, work, nwork, vals, ier)
        if (ier .lt. 0) return
        tpts(nz) = t

c----- Determine the stepsize for this call to the rootfinder

210     call dtnpdr (uv, 1, ider(1), cp, work, nwork, fu, ier)
        if (ier .lt. 0) then
          ier = -100
          return
          endif
        call dtnpdr (uv, 1, ider(2), cp, work, nwork, fv, ier)
        if (ier .lt. 0) then
          ier = -100
          return
          endif

c----- Eliminate the case for single point contours

        if ((abs (fu) .lt. tol) .and. (abs (fv) .lt. tol)) then
          call dtpcu7 (iopt, iv, ih, nz, ndep, upts, mxpts, uv, bnd, c,
     *                 plane, work, nwork, vals, ier)
          tpts(nz) = 0.25d0
          call dtpcu7 (iopt, iv, ih, nz, ndep, upts, mxpts, uv, bnd, c,
     *                 plane, work, nwork, vals, ier)
          tpts(nz) = 0.5d0
          call dtpcu7 (iopt, iv, ih, nz, ndep, upts, mxpts, uv, bnd, c,
     *                 plane, work, nwork, vals, ier)
          tpts(nz) = 0.75d0
          t = 0.0d0
          dist = 1.0d0
          go to 215
          endif

c----- Resume establishing the stepsize to use

        call dtnpdr (uv, 1, ider(3), cp, work, nwork, fuv, ier)
        if (ier .lt. 0) then
          ier = -100
          return
          endif
        call dtnpdr (uv, 1, ider(5), cp, work, nwork, fuu, ier)
        if (ier .lt. 0) then
          ier = -100
          return
          endif
        call dtnpdr (uv, 1, ider(6), cp, work, nwork, fvv, ier)
        if (ier .lt. 0) then
          ier = -100
          return
          endif
        ptol = abs (fvv * fu**2 + fuu * fv**2 - 2.0d0 * fu * fv * fuv)
        if (ptol .ne. 0.0d0) then
          ptol = 2.0d0 * dtmcon(12) * (fu ** 2 + fv ** 2) ** 1.5 / ptol
        else
          ptol = dtmcon (2)
          endif
        ptol = 2.0d0 * min (ptol, 2.0d0) / mxpts
        arclen = t + ptol

c----- Call the rootfinder to get the next point

        call dtpcu2 (cp, uv, plormn, ptol, tol, ier)
        if (ier .lt. 0) then
          ier = -100
          return
          endif
        if (ier .eq. 5) then
          ier = 0
          go to 215
          endif
        t = arclen

c----- Add this point to the list of points so far

        call dtpcu7 (iopt, iv, ih, nz, ndep, upts, mxpts, uv, bnd, c,
     *               plane, work, nwork, vals, ier)
        if (ier .lt. 0) return
        tpts(nz) = t

c----- Determine if more points need to be added to the list

        dist = sqrt ((uv(1) - epts(1,ix)) ** 2 +
     *               (uv(2) - epts(2,ix)) ** 2)
        if ((dist .gt. 2.0d0 * ptol) .or. (nz .le. 3)) go to 210

c----- Add the final point for this contour

215     call dtpcu7 (iopt, iv, ih, nz, ndep, upts, mxpts, epts(1,ix),
     *               bnd, c, plane, work, nwork, vals, ier)
        if (ier .lt. 0) return
        tpts(nz) = t + dist

c----- Now solve the overdetermined system of equations

        call dtpcu3 (nz, upts, mxpts, nfit, tpts, work(iw1), tol,
     *               work(iw4), nleft, work(iw2), work(iw3),
     *               cuts(icptr(ic)), iwarn, ier)
        if (ier .lt. 0) then
          ier = -100
          return
          endif

c----- Fit the second curve in both fit cases (iopt = 3)

        if (iopact .eq. 3) then
          call dtsutl (cuts(icptr(ic)), 1, icptr(ic+1), ier)
          if (ier .lt. 0) then
            ier = -100
            return
            endif
          icptr(ic+1) = icptr(ic+1) + icptr(ic)
          ic = ic + 1
          call dtpcu8 (iv, ih, nz, ndep, upts, mxpts, c, plane,
     *                 work, nwork, vals, ier)
          if (ier .lt. 0) return
          call dtpcu3 (nz, upts, mxpts, ndep, tpts, work(iw1), tol,
     *                 work(iw4), nleft, work(iw2), work(iw3),
     *                 cuts(icptr(ic)), iwarn, ier)
          if (ier .lt. 0) then
            ier = -100
            return
            endif
          endif
 
c----- Determine the length of this vector and prepare for next section

        call dtsutl (cuts(icptr(ic)), 1, icptr(ic+1), ier)
        if (ier .lt. 0) then
          ier = -100
          return
          endif
        icptr(ic+1) = icptr(ic+1) + icptr(ic)
        ic = ic + 1
220     continue

c----- Set the error and warning flag and return

      ier = iwarn
      end

c----- This subroutine finds points on the contour

      subroutine dtpcu2 (cp, uv, plormn, ptol, tol, ier)

      integer mxiter
      parameter (mxiter = 20)
      integer ider(3), nwork, nused, iter, ier
      double precision cp(*), uv(*), ptol, tol, plormn, use(2)
      double precision f1, f2, det, fu, fv, fu2fv2, uvn(2), j11, j12

c----- Initialize the partial derivative indices

      ider(1) = 1
      ider(2) = 0
      ider(3) = 1

c----- Chop up the workspace for this subroutine

      nwork = 2.0d0 * (cp(3) ** 2 + cp(3) + cp(4) ** 2 + cp(4) + 1.0d0)
      nused = 9.0d0 + cp(3) + cp(4) + cp(5) + cp(6) + cp(5) * cp(6)

c----- Determine the initial parameter values

      use(1) = max (min (uv(1), 1.0d0), 0.0d0)
      use(2) = max (min (uv(2), 1.0d0), 0.0d0)

c----- Determine the initial guess

      call dtnpdr (use, 1, ider(1), cp, cp(nused), nwork, fu, ier)
      if (ier .lt. 0) return
      call dtnpdr (use, 1, ider(2), cp, cp(nused), nwork, fv, ier)
      if (ier .lt. 0) return
      fu2fv2 = sqrt (fu ** 2 + fv ** 2)
      if (fu2fv2 .ne. 0.0d0) then
        uvn(1) = uv(1) + ptol * plormn * fv / fu2fv2
        uvn(2) = uv(2) - ptol * plormn * fu / fu2fv2
      else
        uvn(1) = uv(1) + ptol * plormn
        uvn(2) = uv(2)
        endif
      iter = 0

c----- Perform the nonlinear Newton update by inverting the Jacobian

10    continue
      iter = iter + 1
      use(1) = max (min (uvn(1), 1.0d0), 0.0d0)
      use(2) = max (min (uvn(2), 1.0d0), 0.0d0)
      f1 = (use(1) - uv(1)) ** 2 + (use(2) - uv(2)) ** 2 - ptol ** 2
      call dtnpvl (use, 1, cp, cp(nused), nwork, f2, ier)
      if (ier .lt. 0) return
      if (sqrt (f1 ** 2 + f2 ** 2) .lt. tol) then
        uv(1) = use(1)
        uv(2) = use(2)
        return
        endif
      if (iter .gt. mxiter) then
        ier = 5
        return
        endif
      call dtnpdr (use, 1, ider(1), cp, cp(nused), nwork, fu, ier)
      if (ier .lt. 0) return
      call dtnpdr (use, 1, ider(2), cp, cp(nused), nwork, fv, ier)
      if (ier .lt. 0) return
      j11 = 2.0d0 * (use(1) - uv(1))
      j12 = 2.0d0 * (use(2) - uv(2))
      det = j11 * fv - j12 * fu

c----- Regularize if the Jacobian is singular

      if (abs (det) .lt. tol) then
        if (j11 .eq. 0.0d0) j11 = tol
        if (fv .eq. 0.0d0) fv = tol
        j11 = j11 + sign (1.0d0, j11)
        fv = fv + sign (1.0d0, fv)
        det = j11 * fv - j12 * fu
        endif

c----- Update the guess using the new step

      uvn(1) = use(1) - (fv * f1 - j12 * f2) / det
      uvn(2) = use(2) - (j11 * f2 - fu * f1) / det
      go to 10
      end

c----- The following subroutine fits curves to a set of points

      subroutine dtpcu3 (npts, pts, mxpts, ndep, parms, amat, tol,
     *                   work, nwork, rsd, rhs, cfit, iwarn, ier)

      integer mxpts, npts, ix, iy, il, ier, nrow, nintk
      integer nwork, ndep, iwarn
      double precision pts(mxpts,*), parms(*), amat(mxpts,*), cfit(*)
      double precision rsd(*), rhs(*), rpert, error, dmin, dmax, tol
      double precision work(*), region

c----- Handle the case with only two points

      if (npts .lt. 4) then
        cfit(1) = 1
        cfit(2) = ndep
        cfit(3) = 2
        cfit(4) = 2
        cfit(5) = 0.0d0
        cfit(6) = 0.0d0
        cfit(7) = 0.0d0
        cfit(8) = 1.0d0
        cfit(9) = 1.0d0
        do 5 ix = 1,ndep
          call dcopy (2, pts(1,ix), npts-1, cfit(8+2*ix), 1)
5         continue
        return
        endif

c----- Define the knot set for the preliminary no interior knot fit

      cfit(1) = 1
      cfit(2) = ndep
      cfit(3) = 4
      cfit(4) = 4
      cfit(6) = 0.0d0
      cfit(7) = 0.0d0
      cfit(8) = 0.0d0
      cfit(9) = 0.0d0
      cfit(10) = parms(npts)
      cfit(11) = cfit(10)
      cfit(12) = cfit(10)
      cfit(13) = cfit(10)

c----- Generate the linear system that needs to be solved

      do 10 ix = 2,npts-1
        call dtbspl (cfit(6), 8, parms(ix), il, 4, 0, mxpts,
     *               work, nwork, rsd, ier)
        if (ier .lt. 0) return
        call dcopy (4, rsd, 1, amat(ix-1,1), mxpts)
10      continue

c----- Establish the right hand side

      nrow = npts - 2
      call dcopy (nrow, pts(2,1), 1, rhs, 1)
      call daxpy (nrow, -pts(1,1), amat(1,1), 1, rhs, 1)
      call daxpy (nrow, -pts(npts,1), amat(1,4), 1, rhs, 1)

c----- Solve the least squares system using QR factorization

      call dtlsle (amat(1,2), mxpts, nrow, 2, rhs, rsd, work, nwork,
     *             rpert, ier)
      if (ier .lt. 0) return
      error = max (abs (dmin (nrow, rsd, 1)), dmax (nrow, rsd, 1))
      cfit(14) = pts(1,1)
      cfit(15) = rhs(1)
      cfit(16) = rhs(2)
      cfit(17) = pts(npts,1)

c----- Fill in the additional dependent variables using factorization

      do 20 ix = 2,ndep
        call dcopy (nrow, pts(2,ix), 1, rhs, 1)
        call daxpy (nrow, -pts(1,ix), amat(1,1), 1, rhs, 1)
        call daxpy (nrow, -pts(npts,ix), amat(1,4), 1, rhs, 1)
        call dtlssl (amat(1,2), mxpts, nrow, 2, rhs, rsd, work, nwork,
     *               rpert, ier)
        if (ier .lt. 0) return
        error = max (error, abs (dmin (nrow, rsd, 1)),
     *                           dmax (nrow, rsd, 1))
        cfit(10+4*ix) = pts(1,ix)
        cfit(11+4*ix) = rhs(1)
        cfit(12+4*ix) = rhs(2)
        cfit(13+4*ix) = pts(npts,ix)
20      continue

c----- Evaluate and check the error at the points

      nintk = sqrt (sqrt (error / tol))
      if (nintk .ge. npts - 4) then
        nintk = npts - 4
        iwarn = max (iwarn, 1)
        endif

c----- Determine additional knots to add

      cfit(4) = nintk + 4
      cfit(6) = 0.0d0
      cfit(7) = 0.0d0
      cfit(8) = 0.0d0
      cfit(9) = 0.0d0
      do 110 ix = 1,nintk
        region = ix * (npts - 1.0d0) / (nintk + 1.0d0)
        iy = region
        region = region - iy
        cfit(9+ix) = (1.0d0 - region) * parms(iy+1) +
     *               region * parms(iy+2)
110     continue
      cfit(10+nintk) = parms(npts)
      cfit(11+nintk) = cfit(10+nintk)
      cfit(12+nintk) = cfit(10+nintk)
      cfit(13+nintk) = cfit(10+nintk)

c----- Generate the matrix for the new linear system

      do 120 ix = 2,npts-1
        call dtbspl (cfit(6), 8+nintk, parms(ix), il, 4, 0, mxpts,
     *               work, nwork, rsd, ier)
        if (ier .lt. 0) return
        call dcopy (4+nintk, rsd, 1, amat(ix-1,1), mxpts)
120     continue

c----- Fill in the right hand side

      nrow = npts - 2
      call dcopy (nrow, pts(2,1), 1, rhs, 1)
      call daxpy (nrow, -pts(1,1), amat(1,1), 1, rhs, 1)
      call daxpy (nrow, -pts(npts,1), amat(1,4+nintk), 1, rhs, 1)

c----- Solve the least squares system as before, using QR

      call dtlsle (amat(1,2), mxpts, nrow, 2+nintk, rhs, rsd, work,
     *             nwork, rpert, ier)
      if (ier .lt. 0) return
      cfit(14+nintk) = pts(1,1)
      call dcopy (2+nintk, rhs, 1, cfit(15+nintk), 1)
      cfit(17+2*nintk) = pts(npts,1)

c----- Solve for each of the additional variables using factorization

      do 130 ix = 2,ndep
        call dcopy (nrow, pts(2,ix), 1, rhs, 1)
        call daxpy (nrow, -pts(1,ix), amat(1,1), 1, rhs, 1)
        call daxpy (nrow, -pts(npts,ix), amat(1,4+nintk), 1, rhs, 1)
        call dtlssl (amat(1,2), mxpts, nrow, 2+nintk, rhs, rsd, work,
     *               nwork, rpert, ier)
        if (ier .lt. 0) return
        cfit(10+4*ix+ix*nintk) = pts(1,ix)
        call dcopy (2+nintk, rhs, 1, cfit(11+4*ix+ix*nintk), 1)
        cfit(13+4*ix+(ix+1)*nintk) = pts(npts,ix)
130     continue
      end

c----- This is the subroutine which determines the topology workspace

      subroutine dtpcu4 (c0, work, nwork, ncuts, spts, epts, ier)

      double precision c0(*), work(*), spts(2,*), epts(2,*)
      integer nwork, ncuts, ier, needed, each, coutsz
      integer w1, w2, w3, w4, w5, w6, w7, w8

c----- Determine amount of work space and chop this thing up

      needed = max (3*c0(3)+c0(5),3*c0(4)+c0(6)) *
     *         max (c0(3)+2, c0(4)+2) + 3 * (2 * c0(3) + c0(5)) *
     *         (2 * c0(4) + c0(6)) + 2 * (3 * c0(3) + 3 * c0(4) +
     *          c0(5) + c0(6)) + 16
      call dtsutl (c0, 1, coutsz, ier)
      if (ier .lt. 0) return
      each = (nwork - needed - coutsz) / 8
      w1 = 1
      w2 = w1 + coutsz
      w3 = w2 + each
      w4 = w3 + each
      w5 = w4 + each
      w6 = w5 + each
      w7 = w6 + 2 * each
      w8 = w7 + 2 * each

c----- Call the lower level subroutine which has this chopped up right

      call dtpcu5 (c0, work(w8), needed, work(w1), work(w2), work(w3),
     *             work(w4), work(w5), work(w6), work(w7),
     *             ncuts, spts, epts, ier)
      end

c----- This is the subroutine which determines the topology

      subroutine dtpcu5 (c0, work, nwork, cout, ucuts, scrtch, cpts,
     *                   indx, trnpts, traces, ncuts, spts, epts, ier)

      double precision c0(*), work(*), spts(2,*), epts(2,*)
      double precision region(2,2), trimto(2,2), cout(*), ucuts(*)
      double precision scrtch(*), trnpts(2,*), bnd(2), cpts(*), ptol
      double precision dtmcon, wtol, traces(2,*), gradu, gradv
      integer nwork, ncuts, ncur, ier, nucuts, ix, iy, iz, ntpt
      integer tsofar, indx(*), ider(3), jder(2), tptr, iw, cntr
      logical tphere
      external dtpcu6
      common /savblk/ cntr

      data ider /1, 0, 1/
      data jder /0, 2/

c----- Initialize the variables

      ncuts = 0
      ncur = 0
      wtol = dtmcon (6) ** 0.4d0
      ptol = dtmcon (6) ** 0.7d0

c----- Find contours through the bottom boundary

      call dtsutl (c0, 2, region, ier)
      if (ier .lt. 0) return
      trimto(1,1) = region(1,1)
      trimto(1,2) = region(2,1)
      trimto(2,1) = region(1,2)
      trimto(2,2) = region(1,2)
      call dtstrm (c0, trimto, 2, work, nwork, cout, ier)
      if (ier .lt. 0) return
      call dtsutl (cout, 3, bnd, ier)
      if (ier .lt. 0) return
      if (bnd(1) * bnd(2) .le. 0.0d0) then
        call dtszer (cout, work, nwork, nucuts, ucuts, ix, scrtch, ier)
        if (ier .lt. 0) return
      else
        nucuts = 0
        endif

c----- Find contours through the top boundary

      trimto(2,1) = region(2,2)
      trimto(2,2) = region(2,2)
      call dtstrm (c0, trimto, 2, work, nwork, cout, ier)
      if (ier .lt. 0) return
      call dtsutl (cout, 3, bnd, ier)
      if (ier .lt. 0) return
      if (bnd(1) * bnd(2) .le. 0.0d0) then
        call dtszer (cout, work, nwork, iy, ucuts(nucuts+1), ix,
     *               scrtch, ier)
        if (ier .lt. 0) return
      else
        iy = 0
        endif
      nucuts = nucuts + iy

c----- Find all of the turning points

      cntr = 0
      call dtgqnm (2, dtpcu6, c0, region, 30, ptol, 0.9d0, 2, 30,
     *             work, nwork, ntpt, trnpts, ier)
c      write (*,*) 'Number function values = ',cntr
      if (ier .lt. 0) return

c----- Now sort this into an ordered list

      call dcopy (ntpt, trnpts, 2, ucuts(nucuts+1), 1)
      nucuts = nucuts + ntpt + 2
      ucuts(nucuts-1) = region(1,1)
      ucuts(nucuts) = region(2,1)
      call dtsort (ucuts, nucuts, ier)

c----- Remove the duplicates from the list

      do 10 ix = nucuts,2,-1
        if (abs (ucuts(ix) - ucuts(ix-1)) .lt. wtol) then
          call dcopy (nucuts-ix, ucuts(ix+1), 1, ucuts(ix), 1)
          nucuts = nucuts - 1
          endif
10        continue

c----- Sort the turning points into order

      call dcopy (ntpt, trnpts, 2, scrtch, 1)
      if (ntpt .ne. 0) call dtsrtn (scrtch, ntpt, 0, 1, indx, ier)
      call dcopy (ntpt, scrtch, 1, trnpts, 2)
      call dcopy (ntpt, trnpts(2,1), 2, scrtch, 1)
      do 20 ix = 1,ntpt
        trnpts(2,ix) = scrtch(indx(ix))
20      continue

c----- Loop through each of the panels to resolve the topology

      tsofar = 0
      trimto(2,1) = region(1,2)
      trimto(2,2) = region(2,2)
      do 160 ix = 1,nucuts

c----- Find all of the contours passing through this panel edge

        trimto(1,1) = ucuts(ix)
        trimto(1,2) = ucuts(ix)
        call dtstrm (c0, trimto, 2, work, nwork, cout, ier)
        if (ier .lt. 0) return
        iy = 5 + cout(3) + cout(4)
        iz = cout(4)
        do 110 iw = 1,iz
          if (abs (cout(iy+iw)) .lt. wtol) cout(iy+iw) = 0.0d0
110       continue
        call dtsutl (cout, 3, bnd, ier)
        if (ier .lt. 0) return
        if (bnd(1) * bnd(2) .le. 0.0d0) then
          call dtszer (cout, work, nwork, iz, cpts, iy, scrtch, ier)
          if (ier .lt. 0) return
        else
          iz = 0
          endif

c----- See if this panel edge corresponds to a turning point location

        tphere = ((trnpts(1,tsofar+1) .eq. ucuts(ix)) .and.
     *            (tsofar .lt. ntpt))

c----- Make sure that the turning points are correctly resolved

        if (tphere) then
          tsofar = tsofar + 1
          do 120 iy = iz,1,-1
            if (abs (cpts(iy) - trnpts(2,tsofar)) .lt. wtol) then
              call dcopy (iz-iy, cpts(iy+1), 1, cpts(iy), 1)
              iz = iz - 1
              endif
120         continue
          iz = iz + 1
          cpts(iz) = trnpts(2,tsofar)
          call dtsort (cpts, iz, ier)
          if (ier .lt. 0) return
          endif

c----- Loop through for each trace line intersection

        tptr = 1
        do 150 iy = 1,iz
          bnd(1) = ucuts(ix)
          bnd(2) = cpts(iy)

c----- If this is a new point at the bottom boundary, process it

          if (abs (cpts(iy) - region(1,2)) .le. wtol) then
            bnd(2) = region(1,2)
            call dtnpdr (bnd, 1, ider(1), c0, work, nwork, gradu, ier)
            if (ier .lt. 0) return
            call dtnpdr (bnd, 1, ider(2), c0, work, nwork, gradv, ier)
            if (ier .lt. 0) return
            if (ucuts(ix) .eq. region(1,1)) gradu = - gradv
            if (ucuts(ix) .eq. region(2,1)) gradu = gradv
            if (abs (gradu) .lt. wtol) go to 145
            if (gradu * gradv .gt. 0.0d0) then

c----- Contour terminates here

              if (ucuts(ix) .ne. region(1,1)) then
                ncuts = ncuts + 1
                spts(1,ncuts) = traces(1,1)
                spts(2,ncuts) = traces(2,1)
                epts(1,ncuts) = bnd(1)
                epts(2,ncuts) = bnd(2)
                ncur = ncur - 1
                call dcopy (2 * ncur, traces(1,2), 1, traces, 1)
                endif
            else

c----- New contour starts here

              if (ucuts(ix) .ne. region(2,1)) then
                call dcopy (2 * ncur, traces, -1, traces(1,2), -1)
                ncur = ncur + 1
                traces(1,1) = bnd(1)
                traces(2,1) = bnd(2)
                tptr = 2
                endif
              endif
            go to 150
            endif            

c----- If this is a new point at the top boundary, process it

          if (abs (cpts(iy) - region(2,2)) .le. wtol) then
            bnd(2) = region(2,2)
            call dtnpdr (bnd, 1, ider(1), c0, work, nwork, gradu, ier)
            if (ier .lt. 0) return
            call dtnpdr (bnd, 1, ider(2), c0, work, nwork, gradv, ier)
            if (ier .lt. 0) return
            if (ucuts(ix) .eq. region(1,1)) gradu = gradv
            if (ucuts(ix) .eq. region(2,1)) gradu = - gradv
            if (abs (gradu) .lt. wtol) go to 145
            if (gradu * gradv .lt. 0.0d0) then

c----- Contour terminates here

              if (ucuts(ix) .ne. region(1,1)) then
                ncuts = ncuts + 1
                spts(1,ncuts) = traces(1,ncur)
                spts(2,ncuts) = traces(2,ncur)
                epts(1,ncuts) = bnd(1)
                epts(2,ncuts) = bnd(2)
                ncur = ncur - 1
                endif
            else

c----- New contour starts here

              if (ucuts(ix) .ne. region(2,1)) then
                ncur = ncur + 1
                traces(1,ncur) = bnd(1)
                traces(2,ncur) = bnd(2)
                endif
              endif
            go to 150
            endif            

c----- If this is a turning point, decide what to do              

          if ((abs (cpts(iy) - trnpts(2,tsofar)) .lt. wtol)
     *        .and. tphere) then
            call dtnpdr (bnd, 1, ider, c0, work, nwork, gradu, ier)
            if (ier .lt. 0) return
            call dtnpdr (bnd, 1, jder, c0, work, nwork, gradv, ier)
            if (ier .lt. 0) return
            if (abs (gradu * gradv) .lt. wtol) then

c----- Contour forms single point, so register it as such

              ncuts = ncuts + 1
              spts(1,ncuts) = bnd(1)
              spts(2,ncuts) = bnd(2)
              epts(1,ncuts) = bnd(1)
              epts(2,ncuts) = bnd(2)
            else
              if (gradu * gradv .gt. 0.0d0) then

c----- Two branches of contour close on each other here

                if (ucuts(ix) .ne. region(1,1)) then
                  ncuts = ncuts + 1
                  spts(1,ncuts) = traces(1,tptr)
                  spts(2,ncuts) = traces(2,tptr)
                  epts(1,ncuts) = traces(1,tptr+1)
                  epts(2,ncuts) = traces(2,tptr+1)
                  call dcopy (2*(ncur-tptr-1), traces(1,tptr+2), 1,
     *                        traces(1,tptr), 1)
                  ncur = ncur - 2
                  endif
              else

c----- Two new branches of the contour start here

                if (ucuts(ix) .ne. region(2,1)) then
                  call dcopy (2*(ncur-tptr+1), traces(1,tptr), -1,
     *                        traces(1,tptr+2), -1)
                  ncur = ncur + 2
                  traces(1,tptr) = bnd(1)
                  traces(2,tptr) = bnd(2)
                  traces(1,tptr+1) = bnd(1)
                  traces(2,tptr+1) = bnd(2)
                  tptr = tptr + 2
                  endif
                endif
              endif
            go to 150
            endif

c----- Start a new contour if this is the initial panel boundary

          if (ix .eq. 1) then
            do 130 iw = 1,ntpt
              if ((abs (bnd(1) - trnpts(1,iw)) .lt. wtol) .and.
     *            (abs (bnd(2) - trnpts(2,iw)) .lt. wtol)) go to 150
130           continue
            ncur = ncur + 1
            traces(1,ncur) = bnd(1)
            traces(2,ncur) = bnd(2)
            endif

c----- End the contour if this is the final panel boundary

          if (ix .eq. nucuts) then
            do 140 iw = 1,ntpt
              if ((abs (bnd(1) - trnpts(1,iw)) .lt. wtol) .and.
     *            (abs (bnd(2) - trnpts(2,iw)) .lt. wtol)) go to 150
140           continue
            ncuts = ncuts + 1
            spts(1,ncuts) = traces(1,tptr)
            spts(2,ncuts) = traces(2,tptr)
            epts(1,ncuts) = bnd(1)
            epts(2,ncuts) = bnd(2)
            endif
145       continue
          tptr = tptr + 1
150       continue
160     continue

c----- Merge branches of contours together, if possible

      do 210 ix = ncuts,2,-1
        if ((spts(1,ix) .eq. spts(1,ix-1)) .and.
     *      (spts(2,ix) .eq. spts(2,ix-1))) then
          spts(1,ix-1) = epts(1,ix)
          spts(2,ix-1) = epts(2,ix)
          call dcopy (2*(ncuts-ix), spts(1,ix+1), 1, spts(1,ix), 1)
          call dcopy (2*(ncuts-ix), epts(1,ix+1), 1, epts(1,ix), 1)
          ncuts = ncuts - 1
          go to 210
          endif
        if ((spts(1,ix) .eq. epts(1,ix-1)) .and.
     *      (spts(2,ix) .eq. epts(2,ix-1))) then
          epts(1,ix-1) = epts(1,ix)
          epts(2,ix-1) = epts(2,ix)
          call dcopy (2*(ncuts-ix), spts(1,ix+1), 1, spts(1,ix), 1)
          call dcopy (2*(ncuts-ix), epts(1,ix+1), 1, epts(1,ix), 1)
          ncuts = ncuts - 1
          endif
210     continue
      end

c----- Zeros of this function are the turning points

      subroutine dtpcu6 (n, u, c, f, ier)

      integer nwork
      parameter (nwork = 100)
      integer n, ier, ider(2), cntr
      double precision u(*), c(*), f(*), work(nwork)
      common /savblk/ cntr
      data ider /0,1/

      call dtnpvl (u, 1, c, work, nwork, f(1), ier)
      call dtnpdr (u, 1, ider, c, work, nwork, f(2), ier)
      cntr = cntr + 1
      end

c----- This subroutine stores the computed points on the contour

      subroutine dtpcu7 (iopt, iv, ih, nz, ndep, upts, mxpts, uv, bnd,
     *                   c, plane, work, nwork, vals, ier)

      integer iopt, iv, ih, nz, ndep, mxpts, nwork, ier
      double precision upts(mxpts,*), uv(*), work(*), vals(*), c(*)
      double precision uvused(2), bnd(*), plane(*), lgm, ddot

      nz = nz + 1
      if (nz .ge. mxpts) then
        ier = -9
        return
        endif

c----- Store the parameter values of the point

      if (iopt .eq. 1) then
        if (iv .ge. ih) then
          upts(nz,1) = (bnd(2) - bnd(1)) * uv(1) + bnd(1)
          upts(nz,2) = (bnd(4) - bnd(3)) * uv(2) + bnd(3)
        else
          upts(nz,1) = (bnd(4) - bnd(3)) * uv(2) + bnd(3)
          upts(nz,2) = (bnd(2) - bnd(1)) * uv(1) + bnd(1)
          endif
      else

c----- Store the parameter values of the point in temporary space

        if (iv .ge. ih) then
          uvused(1) = (bnd(2) - bnd(1)) * uv(1) + bnd(1)
          uvused(2) = (bnd(4) - bnd(3)) * uv(2) + bnd(3)
        else
          uvused(1) = (bnd(4) - bnd(3)) * uv(2) + bnd(3)
          uvused(2) = (bnd(2) - bnd(1)) * uv(1) + bnd(1)
          endif

c----- Determine the actual point on the planar cut and store it

        call dtnpvl (uvused, 1, c, work, nwork, vals, ier)
        if (ier .lt. 0) then
          ier = -100
          return
          endif
        lgm = (plane(ndep+1) - ddot (ndep, plane, 1, vals, 1)) /
     *        ddot (ndep, plane, 1, plane, 1)
        call daxpy (ndep, lgm, plane, 1, vals, 1)
        call dcopy (ndep, vals, 1, upts(nz,1), mxpts)
        endif
      end

c----- Make the points actual space points in this subroutine

      subroutine dtpcu8 (iv, ih, nz, ndep, upts, mxpts, c, plane, work,
     *                   nwork, vals, ier)

      integer iv, ih, nz, ndep, mxpts, nwork, ier, ix
      double precision upts(mxpts,*), work(*), vals(*), c(*), plane(*)
      double precision lgm, ddot

c----- Loop through for each point

      do 10 ix = 1,nz
        call dtnpvl (upts(ix,1), mxpts, c, work, nwork, vals, ier)
        if (ier .lt. 0) then
          ier = -100
          return
          endif
        lgm = (plane(ndep+1) - ddot (ndep, plane, 1, vals, 1)) /
     *        ddot (ndep, plane, 1, plane, 1)
        call daxpy (ndep, lgm, plane, 1, vals, 1)
        call dcopy (ndep, vals, 1, upts(ix,1), mxpts)
10      continue
      end
