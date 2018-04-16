      subroutine D0SSX3 ( upts, vpts, spts, tpts, ixptr, fit,
     +                    tval, guess, work, nwork,
     +                    cuts, icptr, ncuts, ier )

      double precision upts(*), vpts(*), spts(*), tpts(*),
     +                 fit(*), tval(*), guess(*), work(*), cuts(*)
      integer ixptr(*), nwork, icptr(*), ncuts

      integer ix, iy, icc, ier, next, npts, osize, mcu

c----- Loop through for each cut and determine contour guess

      next = 1
      do 250 ix = 1,ncuts
        icptr(ix) = next
        npts = ixptr(ix+1) - ixptr(ix)

c----- Fill in single point contours

        if (npts .eq. 1) then
          cuts(next) = 1.0d0
          cuts(next+1) = 4.0d0
          cuts(next+2) = 1.0d0
          cuts(next+3) = 1.0d0
          cuts(next+4) = 1.0d0
          cuts(next+5) = 0.0d0
          cuts(next+6) = 1.0d0
          cuts(next+7) = upts(ixptr(ix))
          cuts(next+8) = vpts(ixptr(ix))
          cuts(next+9) = spts(ixptr(ix))
          cuts(next+10)= tpts(ixptr(ix))
          next = next + 11
          go to 250
          endif

c----- Determine the fitting variables

        tval(1) = 0.0d0
        do 210 iy = 2,npts
          tval(iy) = tval(iy-1) + sqrt (
     *         (upts(ixptr(ix)+iy-1) - upts(ixptr(ix)+iy-2)) ** 2 +
     *         (vpts(ixptr(ix)+iy-1) - vpts(ixptr(ix)+iy-2)) ** 2)
 210      continue
        call drscl (npts, tval(npts), tval, 1)
        call dcopy (npts, upts(ixptr(ix)), 1, fit, 1)
        call dcopy (npts, vpts(ixptr(ix)), 1, fit(npts+1), 1)
        call dcopy (npts, spts(ixptr(ix)), 1, fit(2*npts+1), 1)
        call dcopy (npts, tpts(ixptr(ix)), 1, fit(3*npts+1), 1)

c----- Determine initial fit to the contour

        icc = -1
        call D0NPMS (npts, tval, fit, npts, 3, 4, icc,
     +               work, nwork, guess, ier)
        if (ier .lt. 0) return

c----- Insert the extra knots that are required

        tval(npts+1) = 0.0d0
        tval(npts+2) = 1.0d0
        call dcopy (npts+2, tval, 1, tval(npts+3), 1)
        npts = 2 * npts + 4
        call dtsort (tval, npts, ier)
        if (ier .lt. 0) return
        osize = 5 * npts + 1
        call dtoslo (guess, 1, tval, npts, osize, work, nwork,
     *               cuts(next), ier)
        if (ier .lt. 0) return
        call dtsutl ( cuts(next), 1, mcu, ier )
        if (ier .lt. 0) return
        next = next + mcu

  250   continue

      icptr(ncuts+1) = next
      return
      end
