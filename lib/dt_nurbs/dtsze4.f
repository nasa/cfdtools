c----- This subroutine determines distances from points

      subroutine dtsze4 (x, c, tol, work, nwork, dist, ier)

      integer nwork, ier, ix, kord, ndeg, nleft
      double precision x, c(*), tol, work(*), dist

c----- Determine the first non-zero derivative

      kord = c(3)
      ndeg = kord - 1
      nleft = nwork - kord
      call dterpt (1)
      call dtspdr (x, ndeg, c, work(kord+1), nleft, work, 1, ier)
      call dterpt (1)
      if (ier .ne. 0) return
      do 10 ix = 1,ndeg
        if (abs (work(ix+1)) .gt. tol) go to 20
 10     continue

c----- All derivatives are zero, so fudge the derivative

      work(ix+1) = 1.0d0

c----- The distance should be dependent on the value of the derivative

 20   ndeg = ix
      dist = tol
      do 30 ix = 2,ndeg
        dist = dist * ix
 30     continue
      dist = 2.0d0 * (dist / abs (work(ndeg+1))) ** (1.0d0 / ndeg)
      end
