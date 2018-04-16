        subroutine dtssi3(p1, t1, p2, t2, q, ier)
C  **
C  **           Find the point q closest to the two lines 
C  **           p1 + lambda * t1 and p2 + mu * t2 where t1, t2
C  **           are unit vectors.
C  **
C  **   input:
C  **
C  **           p1      first base point
C  **
C  **           t1      first unit direction vector
C  ** 
C  **           p2      second base point 
C  ** 
C  **           t2      second unit direction vector 
C  **
C  **   output:
C  **
C  **           q       closest point
C  **
C  **           ier     success/error flag
C  **                   =  0    success
C  **                   = -40   lines are parallel
C  **
        double precision p1(*), t1(*), p2(*), t2(*), q(*)
        integer ier
C  **
C  **           internal variables
C  **
        double precision z(3), a, b, d, r, r1, r2
        double precision ddot, dtmcon
        integer incx, incy, n, i
        incx = 1
        incy = 1
        n = 3
        do 10 i = 1, 3
10      z(i) = p2(i) - p1(i)
        r  = ddot(n, t1, incx, t2, incy)
        r1 = ddot(n, t1, incx, z,  incy)
        r2 = ddot(n, t2, incx, z,  incy)
        d = 1.d0 - r ** 2
        ier = -40
        if(d .le. dtmcon(5))return
        ier = 0
        a = ( r1 - r * r2) / d
        b = (-r2 + r * r1) / d
        do 20 i = 1, 3
        q(i) = .5 * (p1(i) + a * t1(i) + p2(i) + b * t2(i))
20      continue
        return
        end
