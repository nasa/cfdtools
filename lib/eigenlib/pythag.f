      real function pythag (a, b)
      real a,b
c
c     finds sqrt (a**2+b**2) without overflow or destructive underflow
c
c     06/02/05  D. Saunders  Generic intrinsics retrofitted.

      real p,r,s,t,u
      p = max (abs(a), abs(b))
      if (p .eq. 0.0e0) go to 20
      r = (min (abs(a), abs(b))/p)**2
   10 continue
         t = 4.0e0 + r
         if (t .eq. 4.0e0) go to 20
         s = r/t
         u = 1.0e0 + 2.0e0*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      end
