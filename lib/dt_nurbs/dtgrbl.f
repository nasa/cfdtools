        Subroutine DTGRBL(cru, ndimu, v, nv
     *                   ,crv, ndimv, u, nu
     *                   ,work, nwork, surf, ier)
C  **
C  **           Subroutine for blending a network of curves into 
C  **           a surface.
C  **
C  **           Input variables:
C  **
        double precision cru(ndimu, *), v(*)
     *                  ,crv(ndimv, *), u(*)
     *                  ,work(*), surf(*)
        integer nv, nwork, ier
        integer i
C  **
C  **           Input:
C  **
C  **                   cru     Array holding the u-input curves.
C  **                           Each column of cru represents one
C  **                           input curve.
C  **
C  **                   ndimu   row dimension of cru.
C  **
C  **                   v       Array holding constant parameter 
C  **                           values for the input curves.
C  **
C  **                   nv      Number of curves (also, length of 
C  **                           v array).
C  **
C  **                   crv     Array holding the v-input curves.
C  **                           Each column of crv represents one
C  **                           input curve.
C  **
C  **                   ndimv   row dimension of crv.
C  **
C  **                   u       Array holding constant parameter 
C  **                           values for the input curves.
C  **
C  **                   nu      Number of curves (also, length of 
C  **                           u array).
C  **
C  **           Work space:
C  **
C  **                   work    work array of length nwork
C  **
C  **                   nwork   length of the work array.  
C  **                           nwork .ge. 
C  **                           maxc + nhold + nv + l1 + l2
C  **                           where   
C  **                           nhold = ncv * (3 * kv -2) 
C  **                                   + 2 * kv ** 2
C  **                           maxc  = 2 * ncv + kv + 5
C  **                           ncv   = nv + ndeg - 1
C  **                           kv    = ndeg + 1
C  **                           l1    = nv * (5 + ku + 2 * ncu)
C  **                           l2    = 8 + ku + kv + 
C  **                                   ncu + ncv + ncu * ncv
C  **
C  **           Output:
C  **
C  **                   surf    Array of spline parameters for the 
C  **                           surface.
C  **
C  **                   ier     Success/error flag
C  **
C  **   Calls:  DTCRBL, DTSPRM, DTCNPR, DTERR
C  **
C  **           Internal variables
C  **
        integer nrng, ku, kv, ncu, ncv, ncuv, lenu, locu
        integer lens0, lens, locw, nwork1, iswap(2), iermax
        data iswap / 2, 1 /
C  **
C  **           allocate workspace
C  **
        if (cru(2,1) .ne. crv(2,1)) go to 9008
        nrng = abs( cru(2,1))
        ku  = cru(3, 1)
        kv  = crv(3, 1)
        if (mod (ku, 2) .ne. 0 .or. mod (kv, 2) .ne. 0) go to 9009
        ncu = cru(4, 1)
        ncv = crv(4, 1)
        ncuv = ncu * ncv
        lenu = 5 + ku + ncu + nrng * ncu
        lens0 = 8 + ku + ncu + kv + ncv
        lens = lens0 + nrng * ncuv
        nwork1 = lens + nv * lenu + ncuv
        if (nwork .le. nwork1) go to 9007
        locu = 1 + lens
        locw = locu + nv * lenu
        nwork1 = nwork - locw + 1
        iermax = 0
C  **
C  **           first blend the v-curves in the u-direction
C  **           note: first coordinate is v, second is u.
C  **
        Call DTCRBL(crv, ndimv, u, nu, ku-1, work(locw), nwork1 
     *                   , surf, ier)
        if (ier .lt. 0) go to 9099
        iermax = max (iermax, ier)
C  **
C  **           copy surface to workspace and permute the parameter 
C  **           variables to be consistent with the next surface.
C  **
        Call DTSPRM(surf, iswap, work(locw), nwork1, work, ier)
        if (ier .ne. 0) go to 9099
C  **
C  **           second blend the u-curves in the v-direction
C  **           note: first coordinate is u, second is v
C  **
        Call DTCRBL(cru, ndimu, v, nv, kv-1, work(locw), nwork1
     *                   , surf, ier)
        if (ier .lt. 0) go to 9099
        iermax = max (iermax, ier)
C  **
C  **           add the two surfaces
C  **
        do 10 i = lens0+1, lens
          surf(i) = surf(i) + work(i)
10      continue
C  **
C  **           extract u-curves from first surface to use for 
C  **           building correction surface
C  **
        do 20 i = 1, nv
          Call DTCNPR(work, v(i), 2, work(locw), nwork1, 
     *        work(locu + (i-1) * lenu), ier)
          if (ier .ne. 0) go to 9099
20      continue
C  **
C  **           blend the extracted curves to get correction surface
C  **
        Call DTCRBL(work(locu), lenu, v, nv, kv-1, work(locw), nwork1
     *                   , work, ier)
        if (ier .lt. 0) go to 9099
        iermax = max (iermax, ier)
C  ** 
C  **           subtract the correction surface
C  ** 
        do 30 i = lens0+1, lens
          surf(i) = surf(i) - work(i)
30      continue
        ier = iermax
        return
C  **
C  **           error exits
C  **
C  **           workspace too small
9007    ier = -7
        call dterr (2, 'DTGRBL  ', ier, nwork1)
        return
C  **           incompatible types (poly vs rational) or range spaces
9008    ier = -8
        call dterr (1, 'DTGRBL  ', ier, 0)
        return
C  **           orders of curves not even (degrees not odd)
9009    ier = -9
        call dterr (1, 'DTGRBL  ', ier, 0)
        return
C  **           traceback from error detected in called routine
9099    ier = -99
        call dterr (4, 'DTGRBL  ', ier, 0)
        return
        end
