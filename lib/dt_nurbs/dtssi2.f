        subroutine dtssi2(s1, p1, s2, p2, q, eps, gamma, fnu, 
     *                   mxit, mxha, mxmarq, work, nwork, 
     *                   pout1, pout2, info, ier) 
C  **
C  **           Subroutine projects a point q in space onto the 
C  **           intersection of the two surfaces s1 and s2.
C  **
C  **   input:
C  **
C  **           s1      first surface
C  **
C  **           p1      initial point on first surface 
C  **
C  **           s2      second surface
C  **
C  **           p2      initial point on second surface 
C  **
C  **           q       point in space
C  **
C  **           eps     array of tolerances (see dtdist)
C  **
C  **           gamma, fnu, mxit, mxha, mxmarq
C  **                   various control parameters. See dtclspt and
C  **                   dtdist for descriptions.
C  **
C  **   workspace:
C  **
C  **           work    work array of length nwork
C  **
C  **           nwork   length of work array, see dtssi
C  **
C  **   output:
C  **
C  **           pout1   parameter values of common point on the first 
C  **                   surface
C  ** 
C  **           pout2   parameter values of common point on the second
C  **                   surface 
C  **
C  **           info    convergence information - see dtclspt and       
C  **                   dtdist for descriptions
C  **
C  **           ier     success/error flag
C  **                   = 0     success
C  **                   = -60   failed to find a common point
C  **                   = other error return from dtclsp or dtdist
C  **
        double precision s1(*), p1(*), s2(*), p2(*), q(*)
        double precision eps(*), work(*), pout1(*), pout2(*)
        double precision gamma, fnu     
        integer info(2)
        integer mxit, mxha, mxmarq, nwork, ier
C  **
C  **           internal variables
C  **
        double precision d, pin1(2), pin2(2)
        integer niter
C  **
C  **           for each surface, find closest point to q
C  **
        call dtclsp(s1, q, p1, eps, gamma, fnu,
     *              mxit, mxha, mxmarq, work, nwork, pin1, d, niter,
     *              info, ier)
        call dtclsp(s2, q, p2, eps, gamma, fnu, 
     *              mxit, mxha, mxmarq, work, nwork, pin2, d, niter,
     *              info, ier) 
        if( ier .lt. 0) return
C  **
C  **           using pin1 and pin2 as initial quesses, use dtdist to
C  **           find common point.
C  **
        call dtdist (s1, s2, pin1, pin2, eps, gamma, fnu, 
     *              mxit, mxha, mxmarq, work, nwork, 
     *               pout1, pout2, d, niter, info, ier)
        if (ier .lt. 0 )return
        if( (d ** 2) .gt. eps(1))ier = -60
        return
        end
