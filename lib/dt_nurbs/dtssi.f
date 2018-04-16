        subroutine dtssi(s1, p1, s2, p2, npmax, 
     *                   eps, gamma, fnu, mxit, mxha, mxmarq,
     *                   work, nwork, 
     *                   pts, npts, prm1, prm2, info, ier)
C  **
C  **           Given an initial point (s1(p1) = s2(p2)) find
C  **           a sequence of points lying on the intersection
C  **           of two surfaces.
C  **
C  **   input:
C  **
C  **           s1      first surface
C  **
C  **           p1      parameter value on first surface
C  **
C  **           s2      second surface
C  **
C  **           p2      parameter value on second surface
C  **
C  **           npmax   maximum number of points to be found
C  **
C  **           eps(1)  tolerance for accepting a new point
C  **
C  **           eps(2)  tolerance for testing closure 
C  **
C  **          eps(3)  tolerance for controling advance
C  **
C  **           eps(4 - 7), gamma, fnu, mxit, mxha, mxmarq
C  **                   various control parameters, see dtclspt and
C  **                   dtdist for descriptions.
C  **
C  **   workspace:
C  **
C  **           work    work array of length nwork
C  **
C  **           nwork   length of work, 
C  **                   nwork .ge. 5 * k ** 2 + 3 * k + 16
C  **                   where k = max( s1(3), s1(4), s2(3), s2(4) )
C  **
C  **   output:
C  **
C  **           pts     list of points in space that are on the 
C  **                   intersection curve
C  **
C  **           npts    number of points output, npts .le. npmax
C  **
C  **           prm1    corresponding parameter values in parameter
C  **                   domain of the first surface.
C  **
C  **           prm2    corresponding parameter values in parameter
C  **                   domain of the second surface.
C  **
C  **           info    convergence information, info(1) contains
C  **                   convergence information for the intersector,
C  **                   info(2), info(3) are outputs from dtclspt and 
C  **                   dtdist.
C  **                   info(1) = 10    npts .ge. npmax
C  **                           = 20    curve closed on itself
C  **                           = 30    computation stopped, no 
C  **                                   progress was being made
C  **                           = 40    boundary was reached
C  **
C  **           ier     success/error flag
C  **                        =  0    Success
C  **                        = -1    DTNPVL error
C  **                        = -3    DTNPVL error
C  **                        = -6    DTNPVL error
C  **                        = -40   if the number of independent
C  **                                variables of either surface 
C  **                                is not 2
C  **                        = -41   if the number of dependent 
C  **                                variables of either surface 
C  **                                is not 3 
C  **                        = -42   NWORK too small 
C  **                        = -50   DTNPVL error 
C  **                        = -51   DTNPVL error 
C  **                        = -52   DTNPVL error 
C  **                        = -60   Failed using DTSSI2 
C  **                        = -70   Failed in DTSSI1 
        double precision s1(*), s2(*), p1(*), p2(*), eps(*)
        double precision work(*), pts(3, *), prm1(2, *), prm2(2, *)
        double precision gamma, fnu
        integer info(*)
        integer npmax, mxit, mxha, mxmarq, nwork, npts, ier
C  **
C  **           internal variables
C  **
        double precision a1, b1, c1, d1, a2, b2, c2, d2
        double precision step, t2(3), q(3), pout1(2), pout2(2)
        double precision oldstp, ptmp1(2), ptmp2(2), d
        integer nrng,  ku1, kv1, ncu1, ncv1, nku1, nkv1
        integer iflag, ku2, kv2, ncu2, ncv2, nku2, nkv2
        integer ifill, ntmp, iyes, inc
        integer ntmax, ntimes, i, nhalf, j
C  **
C  **           initial error checking 
C  **           (initially, npts, nrng are used as temporaries)
C  **
        npts = s1(1)
        ier  = -40
        if(npts .ne. 2)go to 9000
        nrng = s2(1)
        if(npts .ne. nrng)go to 9000
C  **
C  **           define some useful integers and continue error checking
C  **
        nrng = s1(2)
        if(nrng .lt. 0)nrng = -nrng - 1
        ier = -41
        if(nrng .ne. 3)go to 9000
        nrng = s2(2)
        if(nrng .lt. 0)nrng = -nrng - 1
        if(nrng .ne. 3)go to 9000
        ku1  = s1(3)
        kv1  = s1(4)
        ncu1 = s1(5)
        ncv1 = s1(6)
        nku1 = ncu1 + ku1
        nkv1 = ncv1 + kv1
        ku2  = s2(3)
        kv2  = s2(4)
        ncu2 = s2(5)
        ncv2 = s2(6)
        nku2 = ncu2 + ku2
        nkv2 = ncv2 + kv2
C  **
C  **           check length of work
C  **
        inc = max(ku1, ku2, kv1, kv2)
        inc = 5 * inc ** 2 + 3 * inc + 16
        ier = -42
        if( inc .gt. nwork ) return
C  **
C  **           define boundaries of the parametric rectangles
C  **
        a1 = s1(9)
        b1 = s1(8 + nku1)
        c1 = s1(9 + nku1)
        d1 = s1(8 + nku1 + nkv1)
        a2 = s2(9)
        b2 = s2(8 + nku2)
        c2 = s2(9 + nku2)
        d2 = s2(8 + nku2 + nkv2)
C  **
C  **           now set up controls for future use
C  **
        ier = 0
C  **
C  **           are we starting at a parametric boundary (set iflag = 1)
C  **           or in the interior (set iflag = 0)?
C  **
        iflag = 0
        if(p1(1) .le. a1 .or. p1(1) .ge. b1 .or.
     *     p1(2) .le. c1 .or. p1(2) .ge. d1) iflag = 1 
        if(p2(1) .le. a2 .or. p2(1) .ge. b2 .or.
     *     p2(2) .le. c2 .or. p2(2) .ge. d2) iflag = 1 
C  **
C  **           initialize the fill-in flag, either the curve is
C  **           being extended into uncharted territory (ifill = 0)
C  **           in which case we can only go forward along the tangent
C  **           or we are trying to fill in more points between two
C  **           points that are already set (ifill = 1) in which case 
C  **           we use both points to fill-in.
C  **
        ifill = 0
        step = -2.
      oldstp = -2.
C  **
C  **           initialize output arrays
C  **
        npts = 1
        ntmp = npmax
        inc = 1
        call dtnpvl(p1, inc, s1, work, nwork, pts(1, 1), ier)
        if(ier .eq. 0)go to 7
        if(ier .lt. 0)go to 9000
7       continue
        do 10 j = 1, 2
        prm1(j, 1) = p1(j)
        prm2(j, 1) = p2(j)
10      continue
C  **
C  **           start filling in the curve
C  **
        ntmax = 10
        ntimes = 1
200     continue
        step = .5 * step
        if(npts .eq. npmax)go to 1000
C  **
C  **           Find a point on the tangent line(s) from which 
C  **           to start search for intersection point.
C  **
        call dtssi5(prm1(1, npts), pout1, s1, 
     *              prm2(1, npts), pout2, s2,
     *              step, ifill, work, nwork, q, ier)
        if( ier .lt. 0) ier = -70
        if( ier .lt. 0) go to 9000
C  **
C  **           Find intersection point.
C  **
        call dtssi2(s1, prm1(1, npts), s2, prm2(1, npts), q, 
     *             eps(4), gamma, fnu, mxit, mxha, mxmarq, work, nwork,
     *             pout1, pout2, info(2), ier)
        if(ier .ge. 0)go to 14
        if(ier .lt. 0)ntimes = ntimes + 1
        if(ntimes .lt. ntmax)go to 200
        if( ier .lt. 0) ier = -60
        if(ier .lt. 0)go to 9000
14      continue 
C  **
C  **           A candidate point has been found.  Before adding it 
C  **           to the list, check tolerances.
C  **
C  **           Find a temporary point q closer to the initial points
C  **           defined by prm1(1, npts).
C  **
        call dtnpvl(pout1, inc, s1, work, nwork, t2, ier)
        if(ier .lt. 0)go to 9000
        do 20 i = 1, 3
        q(i) = .5 * (t2(i) + pts(i, npts))
20      continue
C  **
C  **           Find new intersection point.
C  **
        call dtssi2(s1, prm1(1, npts), s2, prm2(1, npts), q, 
     *             eps(4), gamma, fnu, mxit, mxha, mxmarq, work, nwork,
     *             ptmp1, ptmp2, info(2), ier)
        if( ier .lt. 0) ier = -60
        if(ier .lt. 0)go to 9000
C  **
C  **           Check to see if the new point meets the chord-height
C  **           tolerance as required.
C  **
        call dtssi6(prm1(1, npts), ptmp1, pout1, 
     *              prm2(1, npts), ptmp2, pout2, eps(1), iyes)
        if(iyes .eq. 0)go to 30
C  **
C  **           All is within tolerance, add new point to the list
C  **
        ntimes = 1 
        step = 4. * step 
        npts = npts + 1
        do 40 j = 1, 2
        prm1(j, npts) = pout1(j)
        prm2(j, npts) = pout2(j)
40      continue
        call dtnpvl(prm1(1, npts), inc, s1
     *      , work, nwork, pts(1, npts), ier)
        if(ier .lt. 0)go to 9000
C  **
C  **           Go to the check for finish section
C  **
        ifill = 0
        go to 1000
30      continue
C  **
C  **           The mid points are not within tolerance.  Put the pout
C  **           and the ptmp points on the stack and try again.
C  **
C  **           If there is no more room on the stack, write over 
C  **           an old point.
C  **
        ifill = 1
        ntimes = ntimes + 1
        go to 200
1000    continue
C  **
C  **           Check to see if we are finished.
C  **
        info(1) = 10
        if(npts .ge. npmax)go to 9000
C  **
C  **           Are we back to the beginning?
C  **
        d = 0
        do 1010 i = 1, 3
        d = d + (pts(i, 1) - pts(i, npts)) ** 2
1010    continue
        info(1) = 20
        if (d .lt. eps(2))go to 9000
C  **
C  **           Are we making progress?
C  **
        if(npts .eq. 0)go to 1030
        d = 0
        do 1040 i = 1, 3
        d = d + (pts(i, npts) - pts(i, npts - 1)) ** 2
1040    continue
        info(1) = 30
        d = sqrt( d )
        if (d .lt. eps(3))go to 9000
1030    continue
C  **
C  **           Are we on the boundary?
C  **
        if(prm1(1, npts) .le. a1 .or. prm1(1, npts) .ge. b1)go to 1020
        if(prm1(2, npts) .le. c1 .or. prm1(2, npts) .ge. d1)go to 1020
        if(prm2(1, npts) .le. a2 .or. prm2(1, npts) .ge. b2)go to 1020
        if(prm2(2, npts) .le. c2 .or. prm2(2, npts) .ge. d2)go to 1020
C  **
C  **           No, continue to find points.
C  **
        go to 200
1020    continue
C  **
C  **           On the boundary, if iflag = 1 we are done.
C  **
        info(1) = 40
        if(iflag .eq. 1)go to 9000
        if(npts .le. 1)go to 1045
        nhalf = npts/2
        do 1050 i = 1, nhalf
                do 1060 j = 1, 3
                        work(1)               = pts(j, i)
                        pts(j, i)             = pts(j, npts + 1 - i)
                        pts(j, npts + 1 - i)  = work(1)
1060            continue
                do 1070 j = 1, 2
                        work(1)               = prm1(j, i)
                        prm1(j, i)            = prm1(j, npts + 1 - i)
                        prm1(j, npts + 1 - i) = work(1)
1070            continue
                do 1080 j = 1, 2
                        work(1)               = prm2(j, i)
                        prm2(j, i)            = prm2(j, npts + 1 - i)
                        prm2(j, npts + 1 - i) = work(1)
1080            continue
1050    continue
1045    continue
        step = -oldstp
        iflag = 1
        info(10) = 0
        go to 200
9000    return
        end
