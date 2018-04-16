        subroutine dtssi6(p1, q1, r1, p2, q2, r2, eps, iyes)
C  **
C  **           Check to see if points are in tolerance
C  **
C  **
C  **   input:
C  **
C  **           p1      first point in the first parameter space
C  **
C  **           q1      second point in the first parameter space
C  **
C  **           r1      third point in the first parameter space
C  **
C  **           p2      first point in the second parameter space
C  **
C  **           q2      second point in the second parameter space
C  **
C  **           r2      third point in the second parameter space
C  **
C  **           eps     tolerance
C  **
C  **   output:
C  **
C  **           iyes    answer
C  **                   =  0    The second point is not within 
C  **                           tolerance of the line connecting 
C  **                           the first and third points in one 
C  **                           or the other parameter space.
C  **                   =  1    Both second points are within tolerance
C  **
        double precision p1(*), q1(*), r1(*), p2(*), q2(*), r2(*)
        double precision eps(*)
        integer iyes
C  **
C  **           internal variables
C  **
        double precision s(2), v(2), a, b, c
        double precision ddot, dtmcon
        integer ix, iy, n, i, ier
        iyes = 0
        ix = 1
        iy = 1
        n = 2
        do 10 i = 1, n
        s(i) = q1(i) - p1(i)
        v(i) = r1(i) - p1(i)
10      continue
        a = ddot(n, s, ix, v, iy)
        b = ddot(n, v, ix, v, ix)
        c = ddot(n, s, ix, s, iy)
        ier = -40
        if(b .le. dtmcon(5))return
        ier = 0
        a = c - (a ** 2) / b
        if(a .lt. 0.d0)a = 0.d0
C        a = dsqrt(a)
        if(a .gt. eps(1))return
        do 20 i = 1, n
        s(i) = q2(i) - p2(i)
        v(i) = r2(i) - p2(i)
20      continue
        a = ddot(n, s, ix, v, iy)
        b = ddot(n, v, ix, v, ix)
        c = ddot(n, s, ix, s, iy) 
        ier = -40
        if(b .le. dtmcon(5))return
        ier = 0
        a = c - (a ** 2) / b
        if(a .lt. 0.d0)a = 0.d0
C        a = dsqrt(a)
        if(a .gt. eps(1))return
        iyes = 1
        return
        end
