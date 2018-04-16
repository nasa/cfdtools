      subroutine D0SSX5 ( n, p, aux, qf, ier )

      integer nwork
      parameter (nwork = 200)
      integer n, ier, i, ipos, ier1, ier2
      double precision p(*), aux(*), qf(*), work(nwork),
     +                 f(3), g(3), dudv

c     zeros of this function are solutions to the 4 by 4 set
c     of nonlinear equations:

c       f(u,v) - g(s,t) = 0
c       du/dv = 0

c     these give the turning points of the contours

c     the array p contains the variables in the order u, v,
c     s, t.  as a bonus, ds/dv and dt/dv are returned in
c     p(5) and p(6).

c     the array aux contains the c arrays for f and g.  the
c     c array for f starts at aux(2), while aux(1) gives the
c     starting position for the c array for g. 

      ipos = aux(1)
      call dtnpvl ( p, 1, aux(2), work, nwork, f, ier1 )
      call dtnpvl ( p(3), 1, aux(ipos), work, nwork, g, ier2 )

      if ( ier1 + ier2 .lt. 0 ) then
         ier = -200
         return
         endif
      
      do 10 i = 1, 3
         qf(i) = f(i) - g(i)
   10 continue

      call d0ss12 ( aux(2), aux(ipos), p, dudv, ier )
      qf(4) = dudv

      return
      end
