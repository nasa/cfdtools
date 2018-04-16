      subroutine d0ss10 ( cp, dp, fpt, gpt, wtol, isign, ier )

      double precision cp(*), dp(*), fpt(2), gpt(2), wtol
      integer isign, ier

c     given 2 surfaces and the parameter values on each surface
c     of a point where they intersect, d0ss10 finds the derivative
c     of the contour on the first surface and returns its sign

c     cp      spline array for first surface
c     dp      spline array for second surface
c     fpt     parameter values of point on first surface
c     gpt     parameter values of point on second surface
c     wtol    tolerance within which sign should be zero
c     isign   sign of derivative du/dv
c     ier     error flag

      integer nwork
      parameter ( nwork = 200 )
      double precision a(3,3), b(3), work(nwork), rcond, temp

      integer ider1(2), ider2(2), ier1, ier2, ier3, ier4

      data ider1 / 1, 0 /, ider2 / 0, 1 /

      call dtnpdr ( fpt, 1, ider1, cp, work, nwork, a(1,1), ier1 )
      call dscal ( 3, -1.d0, a(1,1), 1 )
      call dtnpdr ( gpt, 1, ider1, dp, work, nwork, a(1,2), ier2 )
      call dtnpdr ( gpt, 1, ider2, dp, work, nwork, a(1,3), ier3 )

      call dtnpdr ( fpt, 1, ider2, cp, work, nwork, b, ier4 )

      if ( ier1 + ier2 + ier3 + ier4 .lt. 0 ) then
         ier = -200
         return
         endif

      call dtgele ( a, 3, 3, b, work, nwork, rcond, ier )
      if ( ier .lt. 0 ) then
         ier = -42
         return
         endif

      temp = 1.d0
      isign = sign ( temp, b(1) )
      if ( abs(b(1)) .lt. wtol ) isign = 0

      return
      end
