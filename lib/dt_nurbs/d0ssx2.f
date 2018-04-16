      subroutine d0ssx2 (cp, dp, prob, work, nwork, aux, ucuts,
     *                   cpts, indx, trnpts, epts, ncuts, nprt,
     *                   upts, vpts, spts, tpts, icptr, ier)

C===============================================================================
C  PURPOSE:  This is the subroutine which determines the topology
C
C  USAGE:
C
C  INPUTS:  
C
C
C
C  WORKING STORAGE:
C
C
C  OUTPUT:  
C
C  ROUTINES CALLED:
C
C  PROGRAMMER:
C
C  CREATION DATE: 
C
C  UPDATE DATE:  6/4/94                     UPDATE VERSION #: 2.01
C  UPDATE NOTICE: Fix from P. Krausher pertaining to DTSSXT -200 problems.
C                 One line added from previous version.  Not completely fixed.
C
C  UPDATE DATE:  6/14/94                     UPDATE VERSION #: 2.02
C  UPDATE NOTICE: Fix from P. Krausher pertaining to DTSSXT -200 problems.
C                 Problem fixed.
C
C===============================================================================
c----- This is the subroutine which determines the topology

      double precision cp(*), dp(*), work(*), prob
      double precision upts(*), vpts(*), spts(*), tpts(*)
      double precision cpts(nprt,*), aux(*), ucuts(*)
      double precision trnpts(6,*), epts(5,*)
      double precision fpt(2), gpt(2), ptol
      double precision dtmcon, wtol, btol, bndary(2,4)
      double precision p(6), dudv, trimto(2,2)
      integer nwork, ncuts, nprt, m, mm
      integer ncur, nseg, ier, nucuts, ix, iy, iz, ntpt
      integer i, mc, md, nbepts, ntepts, k, isign, idir, lc, rc, lp
      integer tsofar, indx(*), icptr(*), tptr, iw, idup
      integer nepts, n1epts, nb2pts, nt2pts, nl2pts, nr2pts, esofar
      integer iopt, nlimit, j, ifnd, ilast, id, itry, ist, d0tbu1
      logical tphere, ephere, found, again
      external d0ssx4, d0ssx5

c----- Initialize the variables

      ncuts = 0
      ncur = 0
      nseg = 0
      wtol = 50.d0 * ( dtmcon (6) ** 0.4d0 )
      ptol = dtmcon (6) ** 0.8d0
      btol = 10.d0 * wtol
      call dtsutl (cp, 1, mc, ier)
      call dtsutl (dp, 1, md, ier)
      call dcopy ( 4, 0.d0, 0, bndary(1,1), 2 )
      call dcopy ( 4, 1.d0, 0, bndary(2,1), 2 )

c----- Find contours through the bottom boundary

      aux(1) = md + 2
      call dcopy ( md, dp, 1, aux(2), 1 )
      trimto(1,1) = 0.0d0
      trimto(1,2) = 1.0d0
      trimto(2,1) = 0.0d0
      trimto(2,2) = 0.0d0
      call dtstrm (cp, trimto, 2, work, nwork, aux(md+2), ier)
      if (ier .lt. 0) then
        ier = -200
        return
        endif
      call dtgqnm ( 3, d0ssx4, aux, bndary, 30, ptol, prob,
     +              5, nprt, work, nwork, nbepts, epts, ier )
      if ( ier .ne. 0 ) then
         if ( ier .eq. 1 ) ier = -41
         if ( ier .lt. 0 ) ier = -200
         return
         endif
      call dcopy ( nbepts, 0.d0, 0, epts(4,1), 5 )
      nepts = nbepts

c----- Find contours through the top boundary

      trimto(2,1) = 1.0d0
      trimto(2,2) = 1.0d0
      call dtstrm (cp, trimto, 2, work, nwork, aux(md+2), ier)
      if (ier .lt. 0) then
        ier = -200
        return
        endif
      call dtgqnm ( 3, d0ssx4, aux, bndary, 30, ptol, prob,
     +              5, nprt-nepts, work, nwork, ntepts,
     +              epts(1,nepts+1), ier )
      if ( ier .ne. 0 ) then
         if ( ier .eq. 1 ) ier = -41
         if ( ier .lt. 0 ) ier = -200
         return
         endif
      call dcopy ( ntepts, 1.d0, 0, epts(4,nepts+1), 5 )
      nepts = nepts + ntepts

      n1epts = nepts
      call dswap (n1epts, epts(1,1), 5, epts(3,1), 5)
      call dswap (n1epts, epts(2,1), 5, epts(4,1), 5)

c----- Find intersections with bottom edge of 2nd surface

       aux(1) = mc + 2
       call dcopy (mc, cp, 1, aux(2), 1)
       trimto(2,1) = 0.d0
       trimto(2,2) = 0.d0
       call dtstrm (dp, trimto, 2, work, nwork, aux(mc+2), ier)
       if (ier .lt. 0) then
         ier = -200
         return
         endif
       call dtgqnm ( 3, d0ssx4, aux, bndary, 30, ptol, prob,
     +               5, nprt-nepts, work, nwork, nb2pts,
     +               epts(1,nepts+1), ier )
      if ( ier .ne. 0 ) then
         if ( ier .eq. 1 ) ier = -41
         if ( ier .lt. 0 ) ier = -200
         return
         endif
       do 10 i = 1, nb2pts
          fpt(1) = epts(1,nepts+i)
          fpt(2) = epts(2,nepts+i)
          gpt(1) = epts(3,nepts+i)
          gpt(2) = 0.d0
          epts(4,nepts+i) = 0.d0
          iopt = 1
          call d0ss11 (iopt, cp, dp, fpt, gpt, wtol,
     +                 epts(5,nepts+i), ier)
          if (ier .lt. 0) return
   10  continue
       nepts = nepts + nb2pts

c----- Find intersections with top edge of 2nd surface

       trimto(2,1) = 1.d0
       trimto(2,2) = 1.d0
       call dtstrm (dp, trimto, 2, work, nwork, aux(mc+2), ier)
       if (ier .lt. 0) then
         ier = -200
         return
         endif
       call dtgqnm ( 3, d0ssx4, aux, bndary, 30, ptol, prob,
     +               5, nprt-nepts, work, nwork, nt2pts,
     +               epts(1,nepts+1), ier )
      if ( ier .ne. 0 ) then
         if ( ier .eq. 1 ) ier = -41
         if ( ier .lt. 0 ) ier = -200
         return
         endif
       do 20 i = 1, nt2pts
          fpt(1) = epts(1,nepts+i)
          fpt(2) = epts(2,nepts+i)
          gpt(1) = epts(3,nepts+i)
          gpt(2) = 1.d0
          epts(4,nepts+1) = 1.d0
          iopt = 1
          call d0ss11 (iopt, cp, dp, fpt, gpt, wtol,
     +                 epts(5,nepts+i), ier)
          if (ier .lt. 0) return
   20  continue
       nepts = nepts + nt2pts

c----- Find intersections with left edge of 2nd surface

       trimto(1,1) = 0.d0
       trimto(1,2) = 0.d0
       trimto(2,1) = 0.d0
       trimto(2,2) = 1.d0
       call dtstrm (dp, trimto, 2, work, nwork, aux(mc+2), ier)
       if (ier .lt. 0) then
         ier = -200
         return
         endif
       call dtgqnm ( 3, d0ssx4, aux, bndary, 30, ptol, prob,
     +               5, nprt-nepts, work, nwork, nl2pts,
     +               epts(1,nepts+1), ier )
      if ( ier .ne. 0 ) then
         if ( ier .eq. 1 ) ier = -41
         if ( ier .lt. 0 ) ier = -200
         return
         endif
       do 30 i = 1, nl2pts
          fpt(1) = epts(1,nepts+i)
          fpt(2) = epts(2,nepts+i)
          gpt(1) = 0.d0
          gpt(2) = epts(3,nepts+i)
          epts(3,nepts+i) = 0.d0
          epts(4,nepts+i) = gpt(2)
          iopt = 1
          call d0ss11 (iopt, cp, dp, fpt, gpt, wtol,
     +                 epts(5,nepts+i), ier)
          if (ier .lt. 0) return
   30  continue
       nepts = nepts + nl2pts

c----- Find intersections with right edge of 2nd surface

       trimto(1,1) = 1.d0
       trimto(1,2) = 1.d0
       call dtstrm (dp, trimto, 2, work, nwork, aux(mc+2), ier)
       if (ier .lt. 0) then
         ier = -200
         return
         endif
       call dtgqnm ( 3, d0ssx4, aux, bndary, 30, ptol, prob,
     +               5, nprt-nepts, work, nwork, nr2pts,
     +               epts(1,nepts+1), ier )
      if ( ier .ne. 0 ) then
         if ( ier .eq. 1 ) ier = -41
         if ( ier .lt. 0 ) ier = -200
         return
         endif
       do 40 i = 1, nr2pts
          fpt(1) = epts(1,nepts+i)
          fpt(2) = epts(2,nepts+i)
          gpt(1) = 1.d0
          gpt(2) = epts(3,nepts+i)
          epts(3,nepts+i) = 1.d0
          epts(4,nepts+i) = gpt(2)
          iopt = 1
          call d0ss11 (iopt, cp, dp, fpt, gpt, wtol,
     +                 epts(5,nepts+i), ier)
          if (ier .lt. 0) return
   40  continue
       nepts = nepts + nr2pts

c----- Remove from epts points that are on top or bottom edge
c      of first surface

       nlimit = nepts
       do 60 j = nlimit, n1epts+1, -1
          if ( (epts(1,j) .lt. wtol) .or.
     +         (1.d0 - epts(1,j) .lt. wtol) .or.
     +         (epts(2,j) .lt. btol) .or.
     +         (1.d0 - epts(2,j) .lt. btol) ) then
             nepts = nepts - 1
             do 50 i = j, nepts
                call dcopy (5, epts(1,i+1), 1, epts(1,i), 1)
   50        continue
          endif
   60  continue

c----- Sort the 2nd surface edge points by u and put u-coordinates
c      in ucuts

       if ( nepts .gt. 0 ) then
         call dcopy (nepts, epts, 5, work, 1)
         call dtsrtn ( work, nepts, 0, 1, indx, ier )
         call d0pr2x ( epts, 5, 5, nepts, indx, work, nwork, ier )
         call dcopy (nepts, epts(1,1), 5, ucuts, 1)
         endif
       nucuts = nepts

c----- Find all of the turning points

      aux(1) = mc + 2
      call dcopy ( mc, cp, 1, aux(2), 1 )
      call dcopy ( md, dp, 1, aux(mc+2), 1 )
      call dtgqnm (4, d0ssx5, aux, bndary, 30, ptol, prob,
     *             6, nprt, work, nwork, ntpt, trnpts, ier)
      if ( ier .ne. 0 ) then
         if ( ier .eq. 1 ) ier = -41
         if ( ier .lt. 0 ) ier = -200
         return
         endif

c----- For each turning point, find ds/dv and dt/dv again

      do 65 j = 1, ntpt
        call d0ss12 ( cp, dp, trnpts(1,j), work, ier )
        if (ier .lt. 0) return
   65   continue

c----- Sort the turning points into order

      if (ntpt .ne. 0) then
        call dcopy (ntpt, trnpts, 6, work, 1)
        call dtsrtn (work, ntpt, 0, 1, indx, ier)
        call d0pr2x ( trnpts, 6, 6, ntpt, indx, work, nwork, ier )
        endif

c----- Add in the averages of the turning points

      if (ntpt .ge. 2) then
        call dcopy (ntpt-1, trnpts, 6, ucuts(nucuts+1), 1)
        call daxpy (ntpt-1, 1.0d0, trnpts(1,2), 6, ucuts(nucuts+1), 1)
        call drscl (ntpt-1, 2.0d0, ucuts(nucuts+1), 1)
        nucuts = nucuts + ntpt - 1
        endif

c----- Now add turning points and the left and right boundaries.
c      Then sort this into an ordered list

      call dcopy (ntpt, trnpts, 6, ucuts(nucuts+1), 1)
      nucuts = nucuts + ntpt + 2
      ucuts(nucuts-1) = 0.0d0
      ucuts(nucuts) = 1.0d0
      call dtsort (ucuts, nucuts, ier)

c----- Remove all of the duplicates in this list

      do 70 ix = nucuts,2,-1
        if (ucuts(ix) - ucuts(ix-1) .lt. wtol) then
          call d0xcpy (nucuts-ix, ucuts, ix+1, 1, ix, 1)
          nucuts = nucuts - 1
          endif
   70   continue

c----- Line up and sort vertical turning points

      ilast = 1
      do 80 j = 1, ntpt + 1
         if ( j .le. ntpt ) then
           id = d0tbu1 ( trnpts(1,j), ucuts, nucuts, 1, 0, ifnd )
           if ( abs(trnpts(1,j) - ucuts(ifnd)) .lt. wtol )
     +        trnpts(1,j) = ucuts(ifnd)
           if ( ifnd .lt. nucuts ) then
             if ( abs(trnpts(1,j) - ucuts(ifnd+1)) .lt. wtol )
     +          trnpts(1,j) = ucuts(ifnd+1)
             endif
           if ( j .eq. 1 ) go to 80
           endif

         if ( trnpts(1,j) .gt. trnpts(1,j-1) .or.
     +        j .eq. ntpt + 1) then
             if ( ilast .gt. 1 ) then
                ist = j - ilast
                call dcopy ( ilast, trnpts(2,ist), 6, work, 1 )
                call dtsrtn ( work, ilast, 0, 1, indx, ier )
                call d0pr2x ( trnpts(1,ist), 6, 6, ilast, indx,
     +                      work, nwork, ier )
                endif
             ilast = 1
           else
             ilast = ilast + 1
             endif
  80     continue

c----- Line up and sort vertical edge points (from 2nd surface)

      ilast = 1
      do 90 j = 1, nepts + 1
         if ( j .le. nepts ) then
           id = d0tbu1 ( epts(1,j), ucuts, nucuts, 1, 0, ifnd )
           if ( abs(epts(1,j) - ucuts(ifnd)) .lt. wtol )
     +        epts(1,j) = ucuts(ifnd)
           if ( ifnd .lt. nucuts ) then
             if ( abs(epts(1,j) - ucuts(ifnd+1)) .lt. wtol )
     +          epts(1,j) = ucuts(ifnd+1)
             endif
           if ( j .eq. 1 ) go to 90
           endif

         if ( epts(1,j) .gt. epts(1,j-1) .or.
     +        j .eq. nepts+1 ) then
             if ( ilast .gt. 1 ) then
                ist = j - ilast
                call dcopy ( ilast, epts(2,ist), 5, work, 1 )
                call dtsrtn ( work, ilast, 0, 1, indx, ier )
                call d0pr2x ( epts(1,ist), 5, 5, ilast, indx,
     +                        work, nwork, ier )
                do 89 k = ist, ist + ilast - 2
                  if (abs(epts(2,k) - epts(2,k+1)) .lt. btol) then
                    nlimit = ntpt
                    do 88 m = nlimit, 1, -1
                      if ( trnpts(1,m) .eq. epts(1,k) .and.
     +                   abs(trnpts(2,m) - epts(2,k)) .lt. btol .and.
     +                   abs(trnpts(2,m) - epts(2,k+1)) .lt. btol )
     +                   then
                        ntpt = ntpt - 1
                        do 85 mm = m, ntpt
                          call dcopy (5, trnpts(1,mm+1), 1,
     +                                   trnpts(1,mm), 1 )
   85                     continue
                        endif
   88                 continue
                    endif
   89             continue
                endif
             ilast = 1
           else
             ilast = ilast + 1
             endif
  90     continue

c----- Loop through each of the panels to resolve the topology

      tsofar = 0
      esofar = 0
      trimto(2,1) = 0.0d0
      trimto(2,2) = 1.0d0
      icptr(1) = 1
      aux(1) = md + 2
      call dcopy ( md, dp, 1, aux(2), 1 )

      do 160 ix = 1,nucuts

c----- Find all of the contours passing through this panel edge

        trimto(1,1) = ucuts(ix)
        trimto(1,2) = ucuts(ix)
        call dtstrm (cp, trimto, 2, work, nwork, aux(md+2), ier)
        if (ier .lt. 0) then
           ier = -200
           return
           endif
        call dtgqnm ( 3, d0ssx4, aux, bndary, 30, ptol, prob,
     *                3, nprt, work, nwork, iz, cpts, ier )
        if ( ier .ne. 0 ) then
           if ( ier .eq. 1 ) ier = -41
           if ( ier .lt. 0 ) ier = -200
           return
           endif
        call d0iptr ( cpts, 3, nprt, work, nwork, k, ier )

c----- Process the case where the entire cut is a contour

c         if ((iy .eq. 1) .and. (scrtch(1) .eq. 0.0d0) .and.
c    *        (scrtch(2) .eq. 1.0d0)) then
c           ncuts = ncuts + 1
c           nseg = nseg + 1
c           upts(icptr(nseg)) = ucuts(ix)
c           vpts(icptr(nseg)) = 0.0d0
c           upts(icptr(nseg)+1) = ucuts(ix)
c           vpts(icptr(nseg)+1) = 1.0d0
c           icptr(nseg+1) = icptr(nseg) + 2
c           go to 160
c           endif
c       else
c         iz = 0
c         endif

c----- See if this panel edge corresponds to a turning point 
c      or 2nd surface edge location

        tphere = ((trnpts(1,tsofar+1) .eq. ucuts(ix)) .and.
     *            (tsofar .lt. ntpt))
        ephere = ((epts(1,esofar+1) .eq. ucuts(ix)) .and.
     *            (esofar .lt. nepts))

c----- Make sure that the edge points are correctly resolved

        if (ephere) then
          do 110 j = esofar + 1, nepts
            if ( epts(1,j) .gt. ucuts(ix) ) go to 111
            found = .false.
            do 105 k = 1, iz
              if ( abs(epts(2,j) - cpts(k,3)) .lt. btol )
     +           found = .true.
  105         continue
            if ( .not. found ) then
               iz = iz + 1
               cpts(iz,3) = epts(2,j)
               cpts(iz,1) = epts(3,j)
               cpts(iz,2) = epts(4,j)
               endif
  110       continue
  111       continue
          endif

c----- Make sure that the turning points are correctly resolved

        if (tphere) then
          do 120 j = tsofar + 1, ntpt
            if ( trnpts(1,j) .gt. ucuts(ix) ) go to 121
            found = .false.
            do 115 k = 1, iz
              if ( abs(trnpts(2,j) - cpts(k,3)) .lt. btol )
     +           found = .true.
  115         continue
            if ( .not. found ) then
               iz = iz + 1
               cpts(iz,3) = trnpts(2,j)
               cpts(iz,1) = trnpts(3,j)
               cpts(iz,2) = trnpts(4,j)
               endif
  120       continue
  121       continue
          endif

c----- Loop through for each trace line intersection and remove
c      duplicate zeros

        if (iz .ne. 0) then
          call dtsrtn (cpts(1,3), iz, 0, 1, indx, ier)
          call d0prmx (cpts(1,1), iz, indx, ier)
          call d0prmx (cpts(1,2), iz, indx, ier)
          idup = 0
          do 125 iy = iz, 2, -1
            if ( cpts(iy,3) - cpts(iy-1,3) .le. btol ) then
              idup = idup + 1
              do 124 iw = iy, iz + 1 - idup
                 call dcopy ( 3, cpts(iw,1), nprt,
     +                           cpts(iw-idup,1), nprt)
  124         continue
              endif
  125       continue
          iz = iz - idup
          endif
        tptr = 1
        lp = lc
        lc = 0
        rc = 0

        do 150 iy = 1,iz
          fpt(1) = ucuts(ix)
          fpt(2) = cpts(iy,3)
          gpt(1) = cpts(iy,1)
          gpt(2) = cpts(iy,2)

c----- If this is a point at the bottom boundary, process it

          if (abs (cpts(iy,3)) .le. btol) then
            esofar = esofar + 1
            fpt(2) = 0.0d0
            call d0ss10 ( cp, dp, fpt, gpt, wtol, isign, ier )
            if ( ier .lt. 0 ) return
            if ( isign .eq. 0) then
               p(1) = fpt(1)
               p(2) = fpt(2)
               p(3) = gpt(1)
               p(4) = gpt(2)
               call d0ss12 ( cp, dp, p, dudv, ier )
               if ( ier .eq. 10 ) go to 145
               call d0ssx6 ( cp, dp, p, wtol, isign, ier )
               endif
            if (isign .eq. 0) go to 145
            if (isign .lt. 0) then

c----- Contour terminates here

              if (ucuts(ix) .ne. 0.0d0) then
                call d0ssx7 (upts, vpts, spts, tpts, icptr, nseg,
     +                       1, nseg+1, fpt, gpt)
                ncur = ncur - 1
                ncuts = ncuts + 1
                rc = rc + 1
                endif
            else

c----- New contour starts here

              if (ucuts(ix) .ne. 1.0d0) then
                call d0ssx7 (upts, vpts, spts, tpts, icptr, nseg,
     +                       0, 1, fpt, gpt)
                ncur = ncur + 1
                nseg = nseg + 1
                tptr = 2
                lc = lc + 1
                endif
              endif
            go to 150
            endif            

c----- If this is a new point at the top boundary, process it

          if (abs (cpts(iy,3) - 1.0d0) .le. btol) then
            esofar = esofar + 1
            fpt(2) = 1.0d0
            call d0ss10 ( cp, dp, fpt, gpt, wtol, isign, ier )
            if ( ier .lt. 0 ) return
            if ( isign .eq. 0) then
               p(1) = fpt(1)
               p(2) = fpt(2)
               p(3) = gpt(1)
               p(4) = gpt(2)
               call d0ss12 ( cp, dp, p, dudv, ier )
               if ( ier .eq. 10 ) go to 145
               call d0ssx6 ( cp, dp, p, wtol, isign, ier )
               isign = -isign
               endif
            if (isign .eq. 0) go to 145
            if (isign .gt. 0) then

c----- Contour terminates here

              if (ucuts(ix) .ne. 0.0d0) then
                call d0ssx7 (upts, vpts, spts, tpts, icptr, nseg,
     +                       ncur, nseg+1, fpt, gpt)
                ncuts = ncuts + 1
                ncur = ncur - 1
                rc = rc + 1
                endif
            else

c----- New contour starts here

              if (ucuts(ix) .ne. 1.0d0) then
                call d0ssx7 (upts, vpts, spts, tpts, icptr, nseg,
     +                       0, ncur+1, fpt, gpt)
                ncur = ncur + 1
                nseg = nseg + 1
                lc = lc + 1
                endif
              endif
            go to 150
            endif            

          lc = lc + 1
          rc = rc + 1

c----- If this is a turning point, decide what to do              

          if ( tphere ) then

             itry = tsofar + 1
             if ( (abs(trnpts(1,itry) - fpt(1)) .lt. wtol) .and.
     +            (abs(trnpts(2,itry) - fpt(2)) .lt. btol) .and.
     +            (abs(trnpts(3,itry) - gpt(1)) .lt. btol) .and.
     +            (abs(trnpts(4,itry) - gpt(2)) .lt. btol) ) then

               tsofar = itry
               call d0ssx6 ( cp, dp, trnpts(1,tsofar), wtol,
     +                       idir, ier )

               if (idir .eq. 0) then

c----- Contour forms single point, so register it as such

                 ncuts = ncuts + 1
                 nseg = nseg + 1
                 upts(icptr(nseg)) = fpt(1)
                 vpts(icptr(nseg)) = fpt(2)
                 spts(icptr(nseg)) = gpt(1)
                 tpts(icptr(nseg)) = gpt(2)
                 icptr(nseg+1) = icptr(nseg) + 1
                 lc = lc - 1
                 rc = rc - 1
               else
                 if (idir .lt. 0) then

c----- Two branches of contour close on each other here

                   if (ucuts(ix) .ne. 0.0d0) then
                     call d0ssx7 (upts, vpts, spts, tpts, icptr, nseg,
     +                            tptr, nseg+1, fpt, gpt)
                     call d0ssx8 (upts, vpts, spts, tpts, icptr,
     +                            tptr, work)
                     call d0ssx7 (upts, vpts, spts, tpts, icptr, nseg,
     +                            tptr, nseg+1, fpt, gpt)
                     ncur = ncur - 2
                     ncuts = ncuts + 1
                     icptr(nseg) = icptr(nseg+1) - 1
                     nseg = nseg - 1
                     lc = lc - 1
                     rc = rc + 1
                     endif
                 else

c----- Two new branches of the contour start here

                   if (ucuts(ix) .ne. 1.0d0) then
                     call d0ssx7 (upts, vpts, spts, tpts, icptr, nseg,
     +                            0, tptr, fpt, gpt)
                     nseg = nseg + 1
                     call d0ssx7 (upts, vpts, spts, tpts, icptr, nseg,
     +                            0, tptr, fpt, gpt)
                     nseg = nseg + 1
                     ncur = ncur + 2
                     tptr = tptr + 2
                     rc = rc - 1
                     lc = lc + 1
                     endif
                   endif
                 endif
               go to 150
               endif
            endif

c----- If this is on the edge of the 2nd surface, decide what to do

          again = .false.
          if ( ephere ) then
  128        continue
             itry = esofar + 1
             if ( abs(epts(1,itry) - fpt(1)) .lt. wtol  .and.
     +            abs(epts(2,itry) - fpt(2)) .lt. btol ) then

c----- If this point on 1st surface corresponds to 2 points on the
c      edge of the 2nd, we will need to take care of two cuts

                if (again) then
                   ier = -40
                   return
                   endif
                again = .false.
                if ( abs(epts(1,itry+1) - fpt(1)) .lt. wtol  .and.
     +               abs(epts(2,itry+1) - fpt(2)) .lt. btol )
     +             again = .true.
             
                esofar = esofar + 1
                isign = epts(5,esofar)
                if ( isign .eq. 0 ) then
                   p(1) = fpt(1)
                   p(2) = fpt(2)
                   p(3) = gpt(1)
                   p(4) = gpt(2)
                   call d0ss12 ( cp, dp, p, dudv, ier )
                   if ( ier .eq. 10 ) then
                      ier = -42
                      return
                      endif
                   call d0ssx6 ( cp, dp, p, wtol, isign, ier )
                   if ( isign .eq. 0 ) then
                      ier = -42
                      return
                      endif
                   endif
                if ( isign .gt. 0 ) then

c----- Start a new contour in the interior

                   call d0ssx7 ( upts, vpts, spts, tpts, icptr,
     +                           nseg, 0, tptr, fpt, gpt )
                   nseg = nseg + 1
                   ncur = ncur + 1
                   tptr = tptr + 1
                   rc = rc - 1
                   if (again) go to 128
                   go to 150
 
                   else

c----- End a contour in the interior

                   call d0ssx7 ( upts, vpts, spts, tpts, icptr,
     +                           nseg, tptr, nseg+1, fpt, gpt )
                   ncuts = ncuts + 1
                   ncur = ncur - 1
                   lc = lc - 1
                   if (again) go to 128
                   go to 150
                   endif
                endif
             endif

c----- Start a new contour if this is the initial panel boundary

          if (ix .eq. 1) then
            do 130 iw = 1,ntpt
              if ((abs (fpt(1) - trnpts(1,iw)) .lt. wtol) .and.
     *            (abs (fpt(2) - trnpts(2,iw)) .lt. btol)) go to 150
 130          continue
            ncur = ncur + 1
            nseg = nseg + 1
            upts(icptr(ncur)) = fpt(1)
            vpts(icptr(ncur)) = fpt(2)
            spts(icptr(ncur)) = gpt(1)
            tpts(icptr(ncur)) = gpt(2)
            icptr(ncur+1) = icptr(ncur) + 1
            go to 150
            endif

c----- End the contour if this is the final panel boundary

          if (ix .eq. nucuts) then
            do 140 iw = 1,ntpt
              if ((abs (fpt(1) - trnpts(1,iw)) .lt. wtol) .and.
     *            (abs (fpt(2) - trnpts(2,iw)) .lt. btol)) go to 150
 140          continue
            call d0ssx7 (upts, vpts, spts, tpts, icptr, nseg,
     +                   tptr, nseg+1, fpt, gpt)
            ncuts = ncuts + 1
            ncur = ncur - 1
            go to 150
            endif
 145      continue
          if (tptr .gt. ncur) then
            ier = -41
c           return
            endif
          tptr = tptr + 1
          call d0ssx7 (upts, vpts, spts, tpts, icptr, nseg,
     +                   0, tptr, fpt, gpt)
          call d0icpy (nseg-tptr+2, icptr(tptr+1), 1, icptr(tptr), 1)
 150      continue
          if ( ix .gt. 1  .and.  lp .ne. rc ) then
             ier = -41
             return
             endif
 160    continue

c----- Merge branches of contours together, if possible
c      (Note that from here on, fpt and gpt are dummy arrays)

 200  continue
      do 220 ix = 1,ncuts-1
        do 210 iy = ix+1,ncuts
          if ((upts(icptr(ix)) .eq. upts(icptr(iy))) .and.
     *        (vpts(icptr(ix)) .eq. vpts(icptr(iy))) .and.
     *        (spts(icptr(ix)) .eq. spts(icptr(iy))) .and.
     *        (tpts(icptr(ix)) .eq. tpts(icptr(iy)))) then
            call d0ssx8 (upts, vpts, spts, tpts, icptr, iy, work)
            call d0ssx7 (upts, vpts, spts, tpts, icptr, ncuts,
     +                   iy, ncuts+1, fpt, gpt)
            icptr(ncuts+1) = icptr(ncuts+1) - 2
            call d0ssx7 (upts, vpts, spts, tpts, icptr, ncuts,
     +                   ix, ncuts+1, fpt, gpt)
            icptr(ncuts) = icptr(ncuts+1) - 1
            ncuts = ncuts - 1
            go to 200
            endif
          if ((upts(icptr(ix)) .eq. upts(icptr(iy+1)-1)) .and.
     *        (vpts(icptr(ix)) .eq. vpts(icptr(iy+1)-1)) .and.
     *        (spts(icptr(ix)) .eq. spts(icptr(iy+1)-1)) .and.
     *        (tpts(icptr(ix)) .eq. tpts(icptr(iy+1)-1))) then
            call d0ssx7 (upts, vpts, spts, tpts, icptr, ncuts,
     +                   iy, ncuts+1, fpt, gpt)
            icptr(ncuts+1) = icptr(ncuts+1) - 2
            call d0ssx7 (upts, vpts, spts, tpts, icptr, ncuts,
     +                   ix, ncuts+1, fpt, gpt)
            icptr(ncuts) = icptr(ncuts+1) - 1
            ncuts = ncuts - 1
            go to 200
            endif
          if ((upts(icptr(ix+1)-1) .eq. upts(icptr(iy))) .and.
     *        (vpts(icptr(ix+1)-1) .eq. vpts(icptr(iy))) .and.
     *        (spts(icptr(ix+1)-1) .eq. spts(icptr(iy))) .and.
     *        (tpts(icptr(ix+1)-1) .eq. tpts(icptr(iy)))) then
            call d0ssx7 (upts, vpts, spts, tpts, icptr, ncuts,
     +                   ix, ncuts+1, fpt, gpt)
            icptr(ncuts+1) = icptr(ncuts+1) - 2
            call d0ssx7 (upts, vpts, spts, tpts, icptr, ncuts,
     +                   iy-1, ncuts+1, fpt, gpt)
            icptr(ncuts) = icptr(ncuts+1) - 1
            ncuts = ncuts - 1
            go to 200
            endif
          if ((upts(icptr(ix+1)-1) .eq. upts(icptr(iy+1)-1)) .and.
     *        (vpts(icptr(ix+1)-1) .eq. vpts(icptr(iy+1)-1))) then
            call d0ssx8 (upts, vpts, spts, tpts, icptr, iy, work)
            call d0ssx7 (upts, vpts, spts, tpts, icptr, ncuts,
     +                   iy, ncuts+1, fpt, gpt)
            icptr(ncuts+1) = icptr(ncuts+1) - 2
            call d0ssx7 (upts, vpts, spts, tpts, icptr, ncuts, i
     +                  x, ncuts+1, fpt, gpt)
            icptr(ncuts) = icptr(ncuts+1) - 1
            ncuts = ncuts - 1
            go to 200
            endif
 210      continue
 220    continue
      end
