      subroutine D0SSX4 ( n, u, aux, f, ier )

      integer nwork
      parameter (nwork = 200)
      integer n, ier, i, ipos, ier1, ier2
      double precision u(*), aux(*), f(*), work(nwork),
     +                 cval(3), sval(3)

c     this function evaluates a space curve and a surface
c     and subtracts the surface from the curve.  thus, zeros
c     of this function are intersections of the curve and
c     surface.

c     u(1) and u(2) are the parameter values for the surface,
c     while u(3) is the parameter value for the curve.

c     the array aux contains the c arrays for both the
c     curve and surface.  aux(1) gives the position of
c     start of the c array for the curve, while the
c     c array for the surface starts at aux(2)

      ier = 0
      ipos = aux(1)
      call dtspvl ( u(3), aux(ipos), work, nwork, cval, ier1 )
      call dtnpvl ( u(1), 1, aux(2), work, nwork, sval, ier2 )

      if ( ier1 + ier2 .lt. 0 ) then
         ier = -200
         return
         endif

      do 10 i = 1, 3
         f(i) = cval(i) - sval(i)
   10 continue

      return
      end
