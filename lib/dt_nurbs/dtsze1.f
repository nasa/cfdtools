c----- This routine is used to evaluate a spline over an interval

      subroutine dtsze1 (c, xleft, xright, work, nwork, vleft, vright)

      integer nwork, low, high, order, offset, knots, endind, ix, iy
      double precision c(*), xleft, xright, work(*), vleft, vright
      double precision alfa1, alfa2

c----- Find the appropriate start and end spans

      order = c(3)
      low = 1
10    if (c(5 + low + order) .ge. xleft) go to 20
      low = low + 1
      go to 10
20    high = low
30    if (c(5 + high + order) .ge. xright) go to 40
      high = high + 1
      go to 30
40    offset = high - low
      endind = offset + order

c----- Copy over the knots that will be needed

      knots = 2 * order + offset
      do 110 ix = 1,knots
        work(ix) = c(ix + low + 4)
110     continue

c----- Copy over the coefficients that will be needed

      do 120 ix = 1,endind
        work(ix + knots) = c(ix + low + 4 + int (c(3)) + int (c(4)))
120     continue

c----- Evaluate the spline at the left end

      do 220 ix = 1,order - 1
        do 210 iy = 1,order - ix
          alfa1 = (work(order + iy) - xleft) /
     &            (work(order + iy) - work(iy + ix))
          alfa2 = 1.0d0 - alfa1
          work(iy + knots) = alfa1 * work(iy + knots) +
     &                       alfa2 * work(iy + knots + 1)
210       continue
220     continue

c----- Replace the knots on the low end

      do 310 ix = 2,order
        work(ix) = xleft
310     continue

c----- Evaluate the spline at the right end

      do 420 ix = 1,order - 1
        do 410 iy = order - ix,1,-1
          alfa1 = (work(order + iy + offset) - xright) /
     &            (work(order + iy + offset) - work(iy + offset + ix))
          alfa2 = 1.0d0 - alfa1
          work(iy + offset + knots + ix) =
     &             alfa1 * work(iy + offset + knots + ix - 1) +
     &             alfa2 * work(iy + offset + knots + ix)
410       continue
420     continue

c----- Determine the interval of the answer

      vleft = work(knots + 1)
      vright = vleft
      do 510 ix = 2,endind
        if (work(ix + knots) .lt. vleft) vleft = work(ix + knots)
        if (work(ix + knots) .gt. vright) vright = work(ix + knots)
510     continue
      end

