      subroutine D1SSXT (c, d, mc, md, tol, prob, mcuts, nwc,
     +                   cp, dp, upts, vpts, spts, tpts,
     *                   ixptr, aux, wsp, cr, work, nwork,
     +                   nprt, ncuts, cuts, icptr, ier)

c----- Declare the variables for this subprogram

      double precision c(*), d(*), tol, prob, cp(*), dp(*),
     +                 upts(*), vpts(*), spts(*), tpts(*),
     +                 aux(*), wsp(*), cr(*), work(*), cuts(*)
      double precision cbnd(4), dbnd(4)
      integer mcuts, ier, ix, ncuts, ixptr(*), icptr(*)
      integer mc, md, next, npts, osize, nwork, n1, n2
      integer j, k, iu, iv, is, it, ncoef, ncu, ncv, nku, nkv, nwc
      integer ipos, nprt, ic, ii, ie, ifi, itv, ig, iused, imv, isp
      logical snglr
      external D0SSX9

c----- Copy surfaces into arrays to be scaled

      call dcopy (mc, c, 1, cp, 1)
      call dcopy (md, d, 1, dp, 1)

c----- Scale the parametric regions to be [0,1] x [0,1]

      call dtsutl (cp, 2, cbnd, ier)
      if (ier .lt. 0) then
        ier = -200
        return
        endif
      n1 = cp(3) + cp(5)
      n2 = cp(4) + cp(6)
      do 110 ix = 1,n1
        cp(8+ix) = (cp(8+ix) - cbnd(1)) / (cbnd(2) - cbnd(1))
 110    continue
      do 120 ix = 1,n2
        cp(8+n1+ix) = (cp(8+n1+ix) - cbnd(3)) / (cbnd(4) - cbnd(3))
 120    continue

      call dtsutl (dp, 2, dbnd, ier)
      if (ier .lt. 0) then
        ier = -200
        return
        endif
      n1 = dp(3) + dp(5)
      n2 = dp(4) + dp(6)
      do 130 ix = 1, n1
        dp(8+ix) = (dp(8+ix) - dbnd(1)) / (dbnd(2) - dbnd(1))
  130   continue
      do 140 ix = 1, n2
        dp(8+n1+ix) = (dp(8+n1+ix) - dbnd(3)) / (dbnd(4) - dbnd(3))
  140   continue

c----- Split up work

       iu = 1
       ic = iu + 2*nprt
       ii = ic + 3*nprt
       it = ii + nprt
       ie = it + 6*nprt

c      (the part of work starting at ie should be 5*nprt long)

c----- Determine the topology of the contours

      call D0SSX2 (cp, dp, prob, work, nwork, aux, wsp(iu),
     +             wsp(ic), wsp(ii), wsp(it), wsp(ie), ncuts,
     +             nprt, upts, vpts, spts, tpts, ixptr, ier)
      if (ier .lt. 0) then

c----- If necessary, try reversing roles of u and v

        call dcopy ( mc, cp, 1, cr, 1 )
        ncu = cp(5)
        ncv = cp(6)
        nku = cp(3) + cp(5)
        nkv = cp(4) + cp(6)
        call dswap ( 2, cr(3), 2, cr(4), 2 )
        call dcopy ( nku, cp(9), 1, cr(9+nkv), 1 )
        call dcopy ( nkv, cp(9+nku), 1, cr(9), 1 )
        do 150 j = 1, int(abs(c(2)))
           ipos = 9 + nku + nkv + (j-1)*ncu*ncv
           call D0IPTR ( cr(ipos), ncu, ncv,
     +                   work, nwork, k, ier)
           if (ier .lt. 0) then
              ier = -200
              return
              endif
  150      continue
        call D0SSX2 (cr, dp, prob, work, nwork, aux, wsp(iu),
     +               wsp(ic), wsp(ii), wsp(it), wsp(ie), ncuts,
     +               nprt, upts, vpts, spts, tpts, ixptr, ier)
        if (ier .lt. 0) return
        npts = ixptr(ncuts+1)
        call dswap ( npts, upts, 1, vpts, 1 )
        endif

c----- Put surface guesses into cuts (reuse wsp)

      ifi = 1
      itv = ifi + 2*nprt
      ig  = itv + 4*nprt

      call D0SSX3 ( upts, vpts, spts, tpts, ixptr,
     +              wsp(ifi), wsp(itv), wsp(ig), work, nwork,
     +              cuts, icptr, ncuts, ier )
      if (ier .lt. 0) then
        ier = -200
        return
        endif

c----- Move guesses to bottom of cuts

      iused = icptr(ncuts+1) - 1
      imv = mcuts - iused
      do 160 j = iused, 1, -1
        cuts(j+imv) = cuts(j)
  160 continue
      do 170 j = 1, ncuts + 1
        icptr(j) = icptr(j) + imv
  170 continue

c----- Put surfaces at start of work for use by dtcntr and D0SSX9

      aux(1) = mc + 2
      call dcopy ( mc, cp, 1, aux(2), 1 )
      call dcopy ( md, dp, 1, aux(mc+2), 1 )

c----- Determine final contours from guesses

      next = 1
      snglr = .false.
      do 250 ix = 1, ncuts
        isp = icptr(ix) - next
        osize = icptr(ix+1) - icptr(ix)
        if ( isp .lt. osize ) then
           ier = -62
           return
           endif

        call dcopy ( osize, cuts(icptr(ix)), 1, cuts(next), 1 )
        icptr(ix) = next

c----- Determine the cut branch corresponding to this one
c      (note that cp is used for working storage)

        call dtcntr (D0SSX9, aux, tol, isp+osize,
     +               cp, nwc, cuts(next), ier)
        if (ier .eq. 1) snglr = .true.
        if ( ier .lt. 0 ) then
           ier = -200
           return
           endif
        if ( ier .eq. 2 ) then
           ier = -62
           return
           endif

c----- Scale the coefficients of cuts(next) back to original size

       ncoef = cuts(next+3)
       iu = next + 4 + cuts(next+2) + ncoef
       iv = iu + ncoef
       is = iv + ncoef
       it = is + ncoef
       do 220 j = 1, ncoef
         cuts(iu+j) = cuts(iu+j)*(cbnd(2) - cbnd(1)) + cbnd(1)
         cuts(iv+j) = cuts(iv+j)*(cbnd(4) - cbnd(3)) + cbnd(3)
         cuts(is+j) = cuts(is+j)*(dbnd(2) - dbnd(1)) + dbnd(1)
         cuts(it+j) = cuts(it+j)*(dbnd(4) - dbnd(3)) + dbnd(3)
  220  continue

c----- Go back and get the cut branch

        next = next + 5 + cuts(next+2) + 5 * cuts(next+3)
 250    continue
      icptr(ncuts+1) = next
      if (snglr) ier = 1
      end
