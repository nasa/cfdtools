        subroutine dtssi5(p1, q1, s1, p2, q2, s2, step, ifill,
     *             work, nwork, q, ier)
C  **
C  **           Find a point in space from which to start the projection
C  **
C  **   input:
C  **
C  **           p1      point in first parameter space
C  **
C  **           q1      point in first parameter space
C  **
C  **           s1      first surface
C  **
C  **           p2      point in second parameter space
C  **
C  **           q2      point in second parameter space
C  **
C  **           s2      second surface
C  **
C  **           step    a step-size
C  **
C  **           ifill   a control integer (see method)
C  **                   
C  **
C  **   workspace:
C  **
C  **           work    work array of length nwork
C  **
C  **           nwork   length of work, see dtssi for length
C  **
C  **   output:
C  **
C  **           q       a point in space
C  **
C  **           ier     success/error flag
C  **                   =  0    success
C  **                   = -70   failure in dtssi1
C  **                   = other failure in dtnpvl
C  **
C  **   method:
C  **
C  **           It is assumed that s1(p1) = s2(p2) and, if ifill = 1, 
C  **           s1(q1) = s2(q2).  
C  **
C  **           If ifill .ne. 1, a point q on the tangent common to 
C  **                          both s1 and s2 at s1(p1) is found.
C  **
C  **           If ifill = 1,  the point q on the intersection of 
C  **                          the tangent common to both surfaces 
C  **                          at s1(p1) with the tangent common at 
C  **                          s1(p2) is found.
C  **
C  **           If ifill .ne. 0, q1 and q2 are not used.  
C  **           They may be dummies.
C  **
        double precision p1(*), q1(*), s1(*), 
     *                   p2(*), q2(*), s2(*), 
     *                   work(*), q(*)
        double precision step
        integer ifill, nwork, ier
C  **
C  **           Internal variables
C  **
        double precision t1(3), t2(3), val(3, 2)
        integer inc
        inc = 1
        call dtnpvl(p1, inc, s1, work, nwork, val, ier)
        if(ifill .eq. 1)
     *  call dtnpvl(q1, inc, s1, work, nwork, val(1, 2), ier)
        ier = 0
        call dtssi1(s1, p1, s2, p2, work, nwork, t1, ier)
        if( ier .lt. 0) ier = -70
        if(ier .lt. 0)return
        if(ifill .eq. 1)
     *  call dtssi1(s1, q1, s2, q2, work, nwork, t2, ier)
        if( ier .lt. 0) ier = -70
        if(ier .lt. 0)return
        if(ifill .eq. 1)
     *  call dtssi3(val, t1, val(1, 2), t2, q, ier)
        if(ier .lt. 0)ifill = 0
        if(ifill .eq. 0)
     *  call dtssi4(val, t1, step, q)
        return
        end
