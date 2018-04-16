      subroutine D0SSX6 ( c, d, p, wtol, idir, ier )

c     this subroutine determines local behavior at a turning
c     point.  c and d are the two spline surface arrays.
c     p has the same meaning as in D0SSX5; that is, it contains
c     u, v, s, t, alpha, and beta, where 

c     u and v are the coordinates of the turning point on f,
c     s and t are the coordinates of the turning point on g,
c     alpha is ds/dv, and beta is dt/dv.

c     on output, idir is -1 if curve opens to the left
c                         0 if there is higher order singularity
c                         1 if curve opens to the right
c     and ier is an error flag

      double precision c(*), d(*), p(6), wtol, rcond
      integer i, ier, idir
      integer ier1, ier2, ier3, ier4, ier5, ier6, ier7
      integer nwork, ider(2)
      parameter (nwork = 200)
      double precision a(3,3), b(3), work(nwork),
     +                 gss(3), gst(3), gtt(3), fvv(3)

c     set up matrix  ( f   -g   -g   )
c                       u    s    t

      ider(1) = 1
      ider(2) = 0
      call dtnpdr ( p, 1, ider, c, work, nwork, a(1,1), ier1 )
      call dtnpdr ( p(3), 1, ider, d, work, nwork, a(1,2), ier2 )
      call dscal ( 3, -1.d0, a(1,2), 1 )
      ider(1) = 0
      ider(2) = 1
      call dtnpdr ( p(3), 1, ider, d, work, nwork, a(1,3), ier3 )
      call dscal ( 3, -1.d0, a(1,3), 1 )

c     set up right hand side

      ider(1) = 2
      ider(2) = 0
      call dtnpdr ( p(3), 1, ider, d, work, nwork, gss, ier4 )
      ider(1) = 1
      ider(2) = 1
      call dtnpdr ( p(3), 1, ider, d, work, nwork, gst, ier5 )
      ider(1) = 0
      ider(2) = 2
      call dtnpdr ( p(3), 1, ider, d, work, nwork, gtt, ier6 )
      call dtnpdr ( p, 1, ider, c, work, nwork, fvv, ier7 )

      if (ier1+ier2+ier3+ier4+ier5+ier6+ier7 .lt. 0) then
        ier = -200
        return
        endif
      do 10 i = 1,3
        b(i) = p(5)**2 * gss(i) + 2.d0*p(5)*p(6)*gst(i) +
     +         p(6)**2 * gtt(i) - fvv(i)
   10 continue

c     solve equations and set idir

      call dtgele ( a, 3, 3, b, work, nwork, rcond, ier )
      if ( ier .eq. -4 ) then
        ier = 0
        b(1) = 0.d0
        endif
      if ( ier .lt. 0 ) then
        ier = -200
        return
        endif
      idir = sign ( 1.d0, b(1) )
      if ( abs(b(1)) .lt. wtol ) idir = 0
      ier = 0

      return
      end
