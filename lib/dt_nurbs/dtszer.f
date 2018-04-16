      subroutine dtszer (c, work, nwork, nroots, zeros, nints,
     *                   zints, ier)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc    The subprogram DTSZER computes all zeros of a spline function
cc    using the interval Newton method.  The arguments in the
cc    calling sequence have the following meaning:
cc
cc    c      - The input spline definition array.
cc    work   - A real array of length nwork used for working storage.
cc    nwork  - The length of the work array.
cc    nroots - The number of computed isolated zeros of the spline.
cc    zeros  - A real array containing the isolated zeros.
cc    nints  - The number of intervals over which the spline is zero.
cc    zints  - A real array containing left and right endpoints
cc             of the intervals.  zints(1,i) and zints(2,i) are
cc             the left and right endpoints of the ith interval.
cc    ier    - The error control flag.  See the BCSLIB user document
cc             for a description of the values of this flag.
cc
cc    Thomas Grandine
cc    November, 1988
cc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c----- Declare the variables that will be needed in this routine

      integer nwork, nroots, ier, kord, koeff, knots, ix, stack, sofar
      integer spadsz, worksz, staksz, workst, stakst, need, nints
      double precision c(*), work(*), zeros(*), xleft, xright, xmid
      double precision dleft, dright, tol, dtmcon, zints(2,*), fvalue
      double precision scale, small
      logical split
      character*8 subnam
      parameter (subnam = 'DTSZER')

c----- Check for input errors

      ier = 0
      nroots = 0
      nints = 0
      if (c(3) .le. 0.0d0) then
        ier = -1
        call dterr (1, subnam, ier, 0)
        return
        endif
      need = 2 * (5 + c(3) + 2 * c(4)) + 100
      if (nwork .lt. need) then
        ier = -3
        call dterr (2, subnam, ier, need)
        return
        endif
      if (c(4) .lt. c(3)) then
        ier = -6
        call dterr (1, subnam, ier, 0)
        return
        endif
      kord = c(3)
      koeff = c(4)
      knots = kord + koeff
      do 10 ix = 2,knots
        if ((c(4+ix) .gt. c(5+ix)) .or. ((ix .le. (koeff + 1)) .and.
     *      (c(4+ix) .ge. c(4+ix+kord)))) then
          ier = -8
          call dterr (1, subnam, ier, 0)
          return
          endif
10      continue
      if (c(1) .ne. 1.0d0) then
        ier = -51
        call dterr (1, subnam, ier, 0)
        return
        endif
      if (c(2) .ne. 1.0d0) then
        ier = -52
        call dterr (1, subnam, ier, 0)
        return
        endif

c----- Preprocess the intervals which are all zero

      sofar = 0
      scale = 0.0d0
      tol = dtmcon (6)
      do 20 ix = 1,koeff
        scale = max (scale, abs (c(5+knots+ix)))
        if (abs (c(5+knots+ix)) .le. tol) then
          sofar = sofar + 1
        else
          sofar = 0
          endif
        if (sofar .eq. kord) then
          if (c(5+ix) .ne. c(6+ix)) then
            nints = nints + 1
            zints(1,nints) = c(5+ix)
            zints(2,nints) = c(6+ix)
          else
            sofar = sofar - 1
            endif
          endif
        if (sofar .gt. kord) zints(2,nints) = c(6+ix)
20      continue

c----- Set up the stopping tolerances

      tol = 10.0d0 * dtmcon (5)
      scale = max (tol, scale * tol)
      small = sqrt (dtmcon (4))
      
c----- Chop up the working storage space

      spadsz = 2 + c(3) + 2 * c(4)
      worksz = max (spadsz, int (3 * c(3) + 2 * c(3) ** 2))
      staksz = nwork - spadsz - worksz
      workst = spadsz + 1
      stakst = workst + worksz

c----- Initialize the stack and the interval

      stack = 0
      xleft = c(5+kord)
      xright = c(6+koeff)
      do 30 ix = nints,1,-1
        call dtsze4 (zints(2,ix), c, scale, work(workst), worksz,
     *               dright, ier)
        if (ier .ne. 0) then
          ier = -200
          call dterr (5, subnam, ier, 0)
          return
          endif
        xmid = zints(2,ix) + dright
        call dtsze2 (xmid, xright, stack, work(stakst))
        call dtsze4 (zints(1,ix), c, scale, work(workst), worksz,
     *               dleft, ier)
        if (ier .ne. 0) then
          ier = -200
          call dterr (5, subnam, ier, 0)
          return
          endif
        xright = zints(1,ix) - dleft
30      continue

c----- Compute the analytic derivative of the input spline

      call dtspad (c, 1, work(workst), worksz, work, ier)
      if (ier .ne. 0) then
        ier = -200
        call dterr (5, subnam, ier, 0)
        return
        endif

c----- Handle the case when the program terminates

50    if (abs (xright - xleft) .le. tol) then
        xmid = 0.5d0 * (xleft + xright)
        xmid = max (xleft, min (xright, xmid))
60      call dterpt (1)
        call dtspvl (xmid, c, work(workst), worksz, fvalue, ier)
        call dterpt (1)
        if (ier .eq. -38) then
          xmid = (1.0e0 + dtmcon (5)) * xmid
          go to 60
          endif
        if (ier .ne. 0) then
          call dterr (1, 'DTSPVL', ier, 0)
          ier = -200
          call dterr (5, subnam, ier, 0)
          return
          endif
        if (abs (fvalue) .le. scale) then
          nroots = nroots + 1
          zeros(nroots) = xmid
          call dtsze3 (xleft, xright, stack, work(stakst))
          if (stack .eq. -1) return
          go to 50
          endif
        endif
      if ((xleft .ge. xright) .and. (stack .ne. -1)) then
        call dtsze3 (xleft, xright, stack, work(stakst))
        if (stack .eq. -1) return
        go to 50
        endif

c----- See if this interval contains any zeros

      call dtsze1 (c, xleft, xright, work(workst), worksz,
     *             dleft, dright)
      if (dleft .lt. scale) dleft = -1.0d0
      if (dright .gt. -scale) dright = 1.0d0
      if ((dleft * dright) .gt. 0.0d0) then
        call dtsze3 (xleft, xright, stack, work(stakst))
        if (stack .eq. -1) return
        go to 50
        endif

c----- Evaluate the function at this point

40    xmid = 0.5d0 * (xleft + xright)
45    call dterpt (1)
      call dtspvl (xmid, c, work(workst), worksz, fvalue, ier)
      call dterpt (1)
      if (ier .eq. -38) then
        xmid = (1.0d0 + dtmcon (5)) * xmid
        go to 45
        endif
      if (ier .ne. 0) then
        call dterr (1, 'DTSPVL', ier, 0)
        ier = -200
        call dterr (5, subnam, ier, 0)
        return
        endif

c----- If this value is small enough, then save this value

      if (abs (fvalue) .lt. scale) then
        nroots = nroots + 1
        zeros(nroots) = xmid
        call dtsze4 (xmid, c, scale, work(workst), worksz, dleft, ier)
        if (ier .ne. 0) then
          ier = -200
          call dterr (5, subnam, ier, 0)
          return
          endif
        dright = xmid + dleft
        dleft = xmid - dleft
        if ((2 * stack + 2) .gt. staksz) then
          ier = 6
          call dterr (0, subnam, ier, 0)
          return
          endif
        call dtsze2 (dright, xright, stack, work(stakst))
        if ((2 * stack + 2) .gt. staksz) then
          ier = 6
          call dterr (0, subnam, ier, 0)
          return
          endif
        xright = dleft
        go to 50
        endif

c----- Evaluate the derivative and take the Newton step

      call dtsze1 (work, xleft, xright, work(workst), worksz,
     *             dleft, dright)
      split = (dleft * dright .lt. 0.0d0)
      if (dleft .eq. 0.0d0) dleft = small
      if (dright .eq. 0.0d0) dright = -small
      dleft = xmid - fvalue / dleft
      dright = xmid - fvalue / dright
      if (dleft .gt. dright) then
        fvalue = dleft
        dleft = dright
        dright = fvalue
        endif

c----- Handle the case when the interval splits into two pieces

      if (split) then
        if ((2 * stack + 2) .gt. staksz) then
          ier = 6
          call dterr (0, subnam, ier, 0)
          return
          endif
        call dtsze2 (dright, xright, stack, work(stakst))
        xright = dleft
      else
        if (dleft .gt. xleft) xleft = dleft
        if (dright .lt. xright) xright = dright
        endif
      go to 50
      end

