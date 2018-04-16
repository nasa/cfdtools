        Subroutine DTCRBL(cr, ndim, v, nv, ndeg, work, nwork
     *                   , surf, ier)
      external dtmcon
      double precision dtmcon
C  **
C  **           Subroutine for blending two or more curves to form
C  **           a surface.
C  **
C  **           Input variables:
C  **
        double precision cr(ndim, *), v(*), work(*), surf(*)
        integer nv, ndeg, nwork, ier
C  **
C  **           Input:
C  **
C  **                   cr      Array holding the input curves.
C  **                           Each column of cr represents one
C  **                           input curve.
C  **
C  **                   ndim    row dimension of cr.
C  **
C  **                   v       Array holding constant parameter 
C  **                           values for the input curves.
C  **
C  **                   nv      Number of curves (also, length of 
C  **                           v array).
C  **
C  **                   ndeg    Degree of the blending.
C  **
C  **           Work space:
C  **
C  **                   work    work array of length nwork
C  **
C  **                   nwork   length of the work array.  
C  **                           nwork .ge. maxc + nhold + nv
C  **                           where   
C  **                           nhold = ncv * (3 * kv -2) 
C  **                                   + 2 * kv ** 2 + 2 * kcv 
C  **                                   + 4 * kv + 9
C  **                           maxc  = 2 * ncv + kv + 5
C  **                           ncv   = nv + ndeg
C  **                           kv    = ndeg + 1
C  **
C  **           Output:
C  **
C  **                   surf    Array of spline parameters for the 
C  **                           surface.
C  **
C  **                   ier     Success/error flag
C  **                           =  0    success
C  **                           =  1    some ill-conditioning may be 
C  **                                   present.
C  **                           =  2    severe ill-conditioning may be 
C  **                                   present
C  **                           = -1    one or more splines are not 
C  **                                   curves
C  **                           = -2    curves not all same type (poly or
C  **                                   rat) or not in same dimension
C  **                                   image space
C  **                           = -3    curves not all same order
C  **                           = -4    curves not all same knot length
C  **                           = -5    curves not all same knot sequence
C  **                           = -6    insufficient work space
C  **                           = -7    nv less than two
C  **                           = -8    v values not increasing
C  **                           = -9    ndeg less than one
C  **                           =-10    internal check failed
C  **                           =-99    traceback from error in dtnsi
C  **
C  **   Calls: DTNSI, DTERR, DTMCON
C  **
C  **           Internal variables
C  **
        integer nrng, ku, kv, nku, nkv, ncu, ncv, locc, locs
        integer i, j, need, iermax, locv, locw, nwork1
        integer maxc, icc, nc, n
        character subnam*8
C  **
        data subnam /'DTCRBL'/
C  **
C  **           Input data checks
C  **
        if (ndeg .lt. 1) go to 9009
        if (nv .lt. 2) go to 9007
        if (cr(1,1) .ne. 1.0d0) go to 9001
        nrng = abs (cr(2,1))
        ku   = cr(3,1)
        ncu  = cr(4,1)
        nku  = ku + ncu
        do 10 i = 2, nv
          if (v(i) .le. v(i-1)) go to 9008
          if (cr(1,i) .ne. cr(1,1)) go to 9001
          if (cr(2,i) .ne. cr(2,1)) go to 9002
          if (cr(3,i) .ne. cr(3,1)) go to 9003
          if (cr(4,i) .ne. cr(4,1)) go to 9004
          do 5 j = 6, 5 + nku
            if (cr(j,i) .ne. cr(j,1)) go to 9005
5         continue
10      continue
C  **
C  **           Work out parameters for blending direction (v)
C  **
        kv = ndeg + 1
        ncv = nv + ndeg - 1
        if (mod (ndeg, 2) .eq. 0) ncv = ncv + 1
        nkv = kv + ncv
        maxc = 5 + nkv + ncv
C  **
C  **           Check size of work array
C  **   
        need = nv + maxc + ncv * (3 * kv - 2) + 2 * kv ** 2
     *         + 2 * ncv + 4 * kv + 9
        if (need .gt. nwork) go to 9006
C  **
C  **           Allocate pieces of work array
C  **
        locv = 1 + nv
        locw = locv + maxc
        nwork1 = nwork - nv - maxc
C  **
C  **           Now start filling in the surf array
C  **
        ier = 0
        surf(1) = 2
        surf(2) = cr(2,1)
        surf(3) = ku
        surf(4) = kv
        surf(5) = ncu
        surf(6) = ncv
        surf(7) = ku
        surf(8) = kv
C  **
C  **           Put in the u - knots.  Leave the v - knots until DTNSI 
C  **           computes them.
C  **
        do 100 i = 1, nku
          surf(8 + i) = cr(5 + i, 1)
100     continue
C  **
C  **           Now begin the construction of the coefficients.
C  **           The first loop builds the coefficients dependent 
C  **           variable by dependent variable.
C  **
C  **           The second loop is a loop through the curves.
C  **
        iermax = 0
        icc = -1
C  **
C  **           locc locates the zero position for the list of
C  **           coefficients for the current dependent variable in
C  **           the cr(*, j) curve spline array.
C  **
        locc = 5 + nku
C  **
C  **           locs locates the zero position for the matrix of
C  **           coefficients for the current dependent variable in 
C  **           the surface spline array.
C  **
        locs = 8 + nku + nkv
        do 300 n = 1, nrng
C  **
C  **           Fit a curve across each coefficient of the u-curves.
C  **
                do 200 i = 1, ncu
C  **
C  **                  Extract the coefficients.
C  **
                       do 120 j = 1, nv
                         work(j) = cr(locc + i, j)
120                    continue
C  **
C  **                  DTNSI fits a spline curve to work(1:nv) and puts
C  **                  it in work(locv:locw-1) using work(locw:...) as
C  **                  its workspace.
C  **
                       Call DTNSI (nv, v, work, ndeg, icc, work(locw),
     *                             nwork1, maxc, work(locv), nc, ier)
                       if (ier .lt. 0) go to 9099
                       iermax = max (iermax, ier)
                       if (nc .ne. maxc) go to 9010
C  **
C  **                  Put the computed coefficients into the surf array
C  **
                       do 180 j = 1, ncv
                         surf(locs + (j-1) * ncu + i) = 
     *                           work( locv + 4 + nkv + j)
180                    continue
200             continue
                locc = locc + ncu
                locs = locs + ncu * ncv
300     continue
C  **
C  **           Retrieve v - knots from last v - curve 
C  **
        do 400 i = 1, nkv
          surf(8 + nku + i) = work(locv + 4 + i)
400     continue
        ier = iermax
        if (ier .gt. 0) call dterr (0, subnam, ier, 0)
        return
C  **
C  **   Error Exits
C  **
C  **   One or more splines are not curves  (cr(1,*) .ne. 1)
9001    ier = -1
        go to 9900
C  **   Curves not all same type (polynomial or rational) or not in same
C  **   dimension image space  (cr(2,*) not all same)
9002    ier = -2
        go to 9900
C  **   Curves not all same order  (cr(3,*) not all same)
9003    ier = -3
        go to 9900
C  **   Curves not all same knot length  (cr(4,*) not all same)
9004    ier = -4
        go to 9900
C  **   Curves not all same knot sequence
9005    ier = -5
        go to 9900
C  **   Insufficient workspace
9006    ier = -6
        call dterr (2, subnam, ier, need)
        surf(1) = dtmcon(1)
        return
C  **   Less than two curves provided
9007    ier = -7
        go to 9900
C  **   v array values not strictly increasing
9008    ier = -8
        go to 9900
C  **   Requested degree is less than one
9009    ier = -9
        go to 9900
C  **   Curve produced by DTNSI does not have the expected spline length
9010    ier = -10
        call dterr (3, subnam, ier, 0)
        surf(1) = dtmcon(1)
        return
C  **   Traceback message after error occurred in DTNSI
9099    ier = -99
        call dterr (4, subnam, ier, 0)
        surf(1) = dtmcon(1)
C  **   Common exit for input data errors
9900    continue
        call dterr (1, subnam, ier, 0)
        surf(1) = dtmcon(1)
        return
        end
