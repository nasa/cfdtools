        subroutine dtssi1(s1, p1, s2, p2, work, nwork, t, ier)
C  **
C  **           Compute tangent vector common to the two surfaces.
C  **           It is assumed that the surfaces intersect at
C  **           parameter pairs p1 and p2 (i.e., s1( p1 ) = s2( p2 ).
C  **
C  **   input:
C  **
C  **           s1      first surface
C  **
C  **           p1      first parameter pair
C  **
C  **           s2      second surface
C  **
C  **           p2      second parameter pair
C  **
C  **   workspace:
C  **
C  **           work    work array of length nwork
C  **
C  **           nwork   length of the work array, see dtssi
C  **
C  **   output:
C  **
C  **           t       common unit tangent vector
C  **
C  **           ier     success/error flag
C  **                   =  0    success
C  **                   = -40   normal vector is essentially zero
C  **                   = other error return from dtsrfn
C  **
        double precision s1(*), p1(*), s2(*), p2(*), work(*), t(*)
        integer nwork, ier
C  **
C  **           internal variables
C  **
        double precision vn1(3), vn2(3), cross
        double precision dtmcon
        integer inc, i
        ier = 0
        inc = 1
C  **
C  **           compute unit surface normals
C  **
        call dtsrfn(p1, inc, s1, work, nwork, vn1, cross, ier)
        if( ier .lt. 0 )return
        call dtsrfn(p2, inc, s2, work, nwork, vn2, cross, ier)
        if( ier .lt. 0 )return
C  **
C  **           compute cross product of normals to get tangent
C  **
        t(1) =  vn1(2) * vn2(3) - vn1(3) * vn2(2)
        t(2) = -vn1(1) * vn2(3) + vn1(3) * vn2(1)
        t(3) =  vn1(1) * vn2(2) - vn1(2) * vn2(1)
C  **
C  **           and normalize
C  **
        cross = t(1) ** 2 + t(2) ** 2 + t(3) ** 2 
        cross = dsqrt(cross)
        ier = -40
        if(cross .le. dtmcon(5))return
        ier = 0
        do 10 i = 1, 3
10      t(i) = t(i) / cross
        return
        end
