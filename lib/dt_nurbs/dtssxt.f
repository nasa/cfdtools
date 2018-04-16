      subroutine dtssxt (c, d, tol, prob, mcuts, work, nwork,
     *                   ncuts, cuts, icptr, ier)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc      The subroutine DTSSXT is used to find all intersections of
cc      two given surfaces.  The surfaces must each have 2 independent
cc      and 3 dependent variables.  The calling sequence for this
cc      subroutine has the following form:
cc
cc      call dtssxt (c, d, tol, prob, work, nwork,
cc     *             ncuts, cuts, icptr, ier)
cc
cc      where the parameters mean
cc
cc      c      - the first input spline surface
cc      d      - the second input spline surface
cc      tol    - the tolerance to which curves are obtained
cc      prob   - the probability of getting all cuts
cc      mcuts  - size of the array cuts
cc      work   - a work array of length nwork
cc      nwork  - the length of the work array, must be at least
cc                  36*mcuts + mc + md + 1, where mc and md
cc                  are the lengths of the arrays c and d
cc      ncuts  - the number of curve components of the planar cut
cc      cuts   - the output spline vectors for each planar cut component
cc      icptr  - the index array which points to the start, inside
cc               cuts, of each of the spline vectors of the solution
cc      ier    - the error control flag
cc               ier = 0     no errors detected
cc               ier = 1     tolerance may not be achieved; increase tol
cc               ier = -1    invalid spline order for c or d
cc               ier = -3    nwork too small
cc               ier = -6    invalid number of B-spline coefficients
cc               ier = -8    invalid knot sequence
cc               ier = -41   unable to resolve topology
cc               ier = -42   high-order singularity found
cc               ier = -51   c(1) <> 2 or d(1) <> 2
cc               ier = -52   c(2) or d(2) other than 3 or -4
c                ier = -61   mcuts < 200 (needed to get started)
cc               ier = -62   mcuts too small (process error)
cc               ier = -150  tol <= 0
cc               ier = -161  prob <= 0 or prob >= 1
cc               ier = -200  unexpected error return from lower level
c
c                (note:  if attempt to resolve topology fails for
c                        any reason, the roles of u and v are
c                        reversed and another attempt is made.
c                        if this also fails, the ier value applies
c                        to this second attempt.)
cc
cc      Developed from DTPCUT by Thomas Grandine and Fritz Klein
cc      June, 1992
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c----- Declare the necessary variables for this program

      double precision c(*), d(*), work(*), cuts(*), prob, tol
      double precision crange(2,3), drange(2,3)
      integer mcuts, nwork, ncuts, icptr(*), ier, icp, iup, ivp, iix
      integer idp, isp, itp, iwork, ntwork, nkmax, iaux, itwork, icr
      integer mc, md, nleft, need, ntc, ntd, nprt, nwc, j, mode
      character*8 subnam
      data subnam /'DTSSXT  '/

c----- Check the input for validity

      call dterpt (0)
      call dtschk (c, ier)
      if (ier .ne. 0) go to 9900
      call dtschk (d, ier)
      if (ier .ne. 0) go to 9900
      if (c(1) .ne. 2) then
        ier = -51
        go to 9900
        endif
      if (d(1) .ne. 2) then
        ier = -51
        go to 9900
        endif
      if ( (c(2) .ne. 3 .and. c(2) .ne. -4) .or.
     +     (d(2) .ne. 3 .and. d(2) .ne. -4) ) then
        ier = -52
        go to 9900
        endif
      if (tol .le. 0.0d0) then
        ier = -150
        go to 9900
        endif
      if ( prob .le. 0.d0 .or. prob .ge. 1.d0) then
        ier = -161
        go to 9900
        endif

c----- Check that mcuts is at least big enough to get started

       if (mcuts .lt. 200) then
         ier = -61
         go to 9900
         endif

c----- Quick check to see if bounding boxes intersect

      call dtsutl (c, 3, crange, ier)
      call dtsutl (d, 3, drange, ier)

      do 5 j = 1, 3
        if ( crange(1,j) .gt. drange(2,j) .or.
     +       drange(1,j) .gt. crange(2,j) )  then

           ncuts = 0
           ier = 0
           return
           endif
    5   continue

c----- Chop up the working storage into appropriate chunks

      call dtsutl (c, 1, mc, ier)
      call dtsutl (d, 1, md, ier)
      nkmax = max (3*c(3) + c(5), 3*c(4) + c(6))
      ntc = nkmax + 16 +
     +      2*abs(c(2))*(2*c(3) + c(5))*(2*c(4) + c(6)) +
     +      6*(c(3) + c(4)) + 2*(c(5) + c(6)) +
     +      max (
     +        (2*c(3) + c(5))*(c(3) + 1),
     +        (2*c(4) + c(6))*(c(4) + 1),
     +        2*nkmax + (2*c(3) + c(5))*(2*c(4) + c(6)) )
      nkmax = max (3*d(3) + d(5), 3*d(4) + d(6))
      ntd = nkmax + 16 +
     +      2*abs(d(2))*(2*d(3) + d(5))*(2*d(4) + d(6)) +
     +      6*(d(3) + d(4)) + 2*(d(5) + d(6)) +
     +      max (
     +        (2*d(3) + d(5))*(d(3) + 1),
     +        (2*d(4) + d(6))*(d(4) + 1),
     +        2*nkmax + (2*d(3) + d(5))*(2*d(4) + d(6)) )
      ntwork = max(ntc, ntd)
      need = 36*mcuts + mc + md + 1
      if (nwork .lt. need) then
        ier = -3
        go to 9900
        endif
      iaux = 1
      icp = iaux + mc + md + 1
      idp = icp + mc
      itwork = idp + md
      icr = itwork + ntwork
      iup = icr + mc
      nleft = nwork - iup + 1
      nprt = nleft/60
      ivp = iup + 10*nprt
      isp = ivp + 10*nprt
      itp = isp + 10*nprt
      iix = itp + 10*nprt
      iwork = iix + 3*nprt
      nwc = nwork - ( mc + md + 1 )

c----- Now call the lower level subroutine to do the computations

      call D1SSXT (c, d, mc, md, tol, prob, mcuts, nwc,
     +             work(icp), work(idp),
     +             work(iup), work(ivp), work(isp), work(itp),
     *             work(iix), work(iaux), work(iwork), work(icr),
     +             work(itwork), ntwork, nprt, ncuts, cuts, icptr, ier)

c----- Now perform the error processing step

 9900 continue
      call dterpt (1)
      if (ier .lt. 0) then
        ncuts = 0
        cuts(1) = -1.d0
        mode = 1
        if (ier .eq. -3) mode = 2
        if (ier .eq. -41  .or.  ier .eq. -42  .or.
     +      ier .eq. -62  .or.  ier .eq. -200)    mode = 3 
        call dterr (mode, subnam, ier, need)
        endif
      if (ier .gt. 0) call dterr (0, subnam, ier, 0)
      end
