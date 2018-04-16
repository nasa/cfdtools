      subroutine d0ss12 ( c, d, p, dudv, ier )

c     given a point (u,v) on a surface defined by c and a point
c     (s,t) on a surface defined by d, this subroutine finds the
c     derivatives of s, t, and u with respect to v, assuming that
c     in range space the two points coincide.  u, v, s, and t are
c     stored in the array p.  on output, ds/dv and dt/dv are
c     stored in p(5) and p(6).

      double precision c(*), d(*), p(6), dudv
      integer ier
      integer nwork, ider1(2), ider2(2)
      integer ier1, ier2, ier3, ier4
      parameter ( nwork = 200 )
      double precision a(3,3), b(3), work(nwork), rcond

      data ider1 / 1, 0 /, ider2 / 0, 1 /

c     set up matrix  ( g    g   -f  ) and right hand side ( f   )
c                       s    t    u                          v

      call dtnpdr ( p(3), 1, ider1, d, work, nwork, a(1,1), ier1 )
      call dtnpdr ( p(3), 1, ider2, d, work, nwork, a(1,2), ier2 )
      call dtnpdr ( p,    1, ider1, c, work, nwork, a(1,3), ier3 )
      call dscal ( 3, -1.d0, a(1,3), 1 )

      call dtnpdr ( p,    1, ider2, c, work, nwork, b,      ier4 )

      if ( ier1 + ier2 + ier3 + ier4 .lt. 0 ) then
         ier = -200
         return
         endif

c     solve equations

      call dtgele ( a, 3, 3, b, work, nwork, rcond, ier )

      if ( ier .eq. -4 ) then
         ier = 10
         dudv = 1.d12

       else if ( ier .lt. 0 ) then
         ier = -200

       else
        p(5) = b(1)
        p(6) = b(2)
        dudv = b(3)
        endif

      return
      end
