      subroutine d0ss11 ( iopt, cp, dp, fpt, gpt, wtol, dir, ier )

      double precision cp(*), dp(*), fpt(2), gpt(2), wtol, dir
      integer ier, iopt

c     given 2 surfaces f(u,v) and g(s,t), and a point of intersection
c     in the interior of f but on an edge of g, d0ss11 determines
c     the direction of the contour with respect to u.  this means
c     finding the sign of either du/dt (if top or bottom edge of g)
c     or du/ds (if left or right edge of g).  this sign is then
c     reversed if top or right edge of g.

c     iopt    = 1 if we are finding sign of du/dt
c             = 2 if we are finding sign of du/ds
c     cp      spline array for first surface
c     dp      spline array for second surface
c     fpt     parameter values of point on first surface
c     gpt     parameter values of point on second surface
c     wtol    tolerance within which sign should be zero
c     dir     direction of contour with respect to u
c     ier     error flag

      integer nwork
      parameter ( nwork = 200 )
      double precision a(3,3), b(3), work(nwork), rcond, temp

      integer ider1(2), ider2(2), ier1, ier2, ier3, ier4

      data ider1 / 1, 0 /, ider2 / 0, 1 /

      if ( gpt(2) .lt. wtol .or. 1.d0 - gpt(2) .lt. wtol ) then
         iopt = 1
        else
         iopt = 2
      endif

      call dtnpdr ( fpt, 1, ider1, cp, work, nwork, a(1,1), ier1 )
      call dtnpdr ( fpt, 1, ider2, cp, work, nwork, a(1,2), ier2 )
      call dtnpdr ( gpt, 1, ider1, dp, work, nwork, a(1,3), ier3 )
      call dtnpdr ( gpt, 1, ider2, dp, work, nwork, b,      ier4 )

      if ( ier1 + ier2 + ier3 + ier4 .lt. 0 ) then
         ier = -200
         return
         endif

      if ( iopt .eq. 2 ) call dswap ( 3, a(1,3), 1, b, 1 )
      call dscal ( 3, -1.d0, a(1,3), 1 )

      call dtgele ( a, 3, 3, b, work, nwork, rcond, ier )
      if ( ier .lt. 0 ) then
         ier = -42
         return
         endif

      temp = 1.d0
      dir = sign ( temp, b(1) )
      if ( abs(b(1)) .lt. wtol ) dir = 0.d0

      if ( gpt(1) .lt. wtol  .or.  gpt(2) .lt. wtol ) then
        continue
       else
        dir = -dir
      endif

      return
      end
