!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine cutorder (nad, NADDIM, x, y, z, f, ICUTDIM, xyzf,
     >                     nnode, MNODE, ndg, nedge, MEDGE,
     >                     epsxyz, epsf, ierr)
!
!  Description:
!
!     Cutorder takes a collection of 2-point line segments (edges) and orders
!     them into polygonal curves.  The input edges are supplied as a point list,
!     xyzf, and a list of edge index pairs, ndg.  The edges are oriented.
!     A list of edges may be obtained, for example, by calling cutsurf, which
!     cuts a structured surface mesh with a given plane.
!
!     The points of the polygonal curve are stored into new arrays x, y, z, f.
!     If two points are close together then the first one is taken and the
!     second one is ignored, according to epsxyz/f as described in cutsurf.
!
!     The reordered points may be traversed by using nad(*) as follows:
!
!        nedge = 0 ! Normal initialization
!        nnode = 0
!
!        call cutsurf (...)  ! For a given cutting plane, giving packed segments
!
!        nad(1) = nedge      ! May be 0, which may need to be known by caller
!
!        if (nedge > 0) then
!
!           call cutorder (...) ! Merge segments to eliminate common end points
!
!           ncrv = nad(1)
!
!           nad(1) = 1 ! So the loop over disconnected curves works for curve 1
!
!           do line = 1, ncrv
!              l1 = nad(line)
!              l2 = nad(line+1) - 1
!              ::::::::::::::     ! Plot points l1:l2 or whatever
!           end do
!
!           nad(1) = ncrv ! Restore it as the number of curves >= 0
!
!        end if
!
!  History:
!
!     09-APR-1998  S.Thomas    Companion for cutsurf to ease plotting of cuts.
!     18-AUG-1998     "        Initial Fortran 90 translation.
!     20-AUG-1998  D.Saunders  Streamlined with cutsurf as a library utility.
!
!  Author:
!
!     Scott D. Thomas, Sterling Software, Inc.,  Contract NAS2-13210/6/5
!                      Raytheon STX Corporation, Contract NAS2-98080/53
!                      NASA Ames Research Center, Moffett Field, CA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

!  Arguments:

      integer, intent (in) ::
     >  NADDIM,                  ! Dimension of nad(*) >= 1 + max. # connected
     >                           ! line segments expected per set of xyzf(*,*)
     >  ICUTDIM,                 ! Max. # points handled as output from the
     >                           ! merging into connected curves done here
     >  nnode,                   ! See cutsurf - not used here
     >  MNODE,                   ! Max. # end points of line segments handled
     >  MEDGE,                   ! Max. # 2-pt. line segments handled per cut
     >  ndg(2,MEDGE),            ! Pointers defining the ends of each line
     >                           ! segment as packed in xyzf(*,1:nnode)
     >  nedge                    ! # 2-pt. line segments in xyzf(*,*), >= 0

      real, intent (in) ::
     >  xyzf(4,MNODE),           ! Packed (x,y,z,f) coordinates of end points
     >                           ! of 2-pt. line segments as from cutsurf
     >  epsxyz, epsf             ! Tolerances for coincident pts. as described
                                 ! in cutsurf

      real, intent (out),
     >  dimension (ICUTDIM) ::   ! Coordinates and scalar function values for
     >  x, y, z, f               ! one or more connected polygonal curves,
                                 ! accessed via nad(*) as indicated above

      integer, intent (out) ::
     >  nad(NADDIM),             ! Pointers into x, y, z, f; usage shown above;
     >                           ! nad(1) >= 0 is the number of connected
     >                           ! curves found
     >  ierr                     ! 0 = no error detected;
                                 ! 1 = NADDIM limit encountered;
                                 ! 2 = NADDIM < 2;
                                 ! 3 = nnode < 0;
                                 ! 4 = nedge < 0;
                                 ! 5 = ICUTDIM limit encountered

!     Local variables:

      integer
     >  ipoint(nedge),
     >  it, itbeg, line, m, m1, m2, n, n1, n2, nstart, ntouch

!     Execution:

      if (nedge < 0) then
        ierr = 4
      else if (nnode < 0) then
        ierr = 3
      else if (NADDIM < 2) then
        ierr = 2
      else
        ierr = 0
      end if

      if (ierr /= 0) go to 999

!     Mark all edges as available.

      ipoint = 0 ! 1:nedge

!     Mark redundant edges as unavailable.
!     Mark collapsed edges that touch other edges as unavailable.
!     This assumes that the point list xyzf is unique to within
!     epsxyz, epsf as explained in cutsurf.

      do n = 1, nedge
        n1 = ndg(1,n)
        n2 = ndg(2,n)
        do m = n + 1, nedge
          m1 = ndg(1,m)
          m2 = ndg(2,m)
          if (n1 == m1 .and. n2 == m2) then
            ipoint(n) = -1
            exit
          else if (n1 == n2) then
            if (n1 == m1 .or. n1 == m2) then
              ipoint(n) = -2
              exit
            end if
          end if
        end do
      end do

!     This algorithm is fairly straightforward but it has a lot of go to
!     statements - sorry.  Hope the comments help.

!     Initialize the point counter, it, and the line counter, line.

      it   = 0
      line = 0

   20 continue

!       Look for a place to start a new line.

        do n = 1, nedge
          if (ipoint(n) == 0) then
            nstart = n
            go to 30
          end if
        end do

        go to 80 ! No place to start - bail out.

   30   continue

!       Start a new line, taking both points of this segment.
!       Cutordap won't accept the second point of a collapsed edge.

        line = line + 1
        if (line + 1 > NADDIM) then
          ierr = 1
          go to 999
        end if

!       Append a segment.

        ipoint(nstart) = line

        call cutordap (it, ndg(1,nstart))

        if (ierr /= 0) go to 999

        itbeg = it

        call cutordap (it, ndg(2,nstart))

        if (ierr /= 0) go to 999

   40   continue

!         Look for a segment touching the point at x(i),y(i),z(i),f(i),
!         where i == it (end of this poly line).

          do n = 1, nedge
            if (ipoint(n) == 0) then
              n1 = ndg(1,n)
              if (abs (x(it) - xyzf(1,n1)) <= epsxyz .and.
     >            abs (y(it) - xyzf(2,n1)) <= epsxyz .and.
     >            abs (z(it) - xyzf(3,n1)) <= epsxyz .and.
     >            abs (f(it) - xyzf(4,n1)) <= epsf) then
                ntouch = n
                go to 50
              end if
            end if
          end do

!         No more segments join the end of the line, so go look for segments
!         to insert at the beginning of the line.  Inserting is more costly.

          go to 60

   50     continue

!         Append this segment, then go back and look for more.

          ipoint(ntouch) = line

          call cutordap (it, ndg(2,ntouch))

          if (ierr /= 0) go to 999

        goto 40

   60   continue

!          Look for a segment touching the point at x(i),y(i),z(i),f(i),
!          where i == itbeg (start of this poly line).

          do n = 1, nedge
            if (ipoint(n) == 0) then
              n2 = ndg(2,n)
              if (abs(x(itbeg) - xyzf(1,n2)) <= epsxyz .and.
     >            abs(y(itbeg) - xyzf(2,n2)) <= epsxyz .and.
     >            abs(z(itbeg) - xyzf(3,n2)) <= epsxyz .and.
     >            abs(f(itbeg) - xyzf(4,n2)) <= epsf) then
                ntouch = n
                go to 70
              end if
            end if
          end do

!         No more touching segments.  End this line and look for next line.

          nad(line+1) = it + 1

        goto 20

   70   continue

!         Insert this segment, then go back and look for more to insert.
!         It takes more work to insert than to append.

          ipoint(ntouch) = line

          call cutordin (itbeg, it, ndg(1,ntouch))

          if (ierr /= 0) go to 999

        go to 60

   80 continue

!     Normal termination.

      nad(1) = line

  999 return

!     Internal procedures for cutorder (cutordap and cutordin).

      contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine cutordap (it, ndg)
!
!     Append a point as x/y/z/f(it+1).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in)    :: ndg ! Pointer into xyzf(*,*)
      integer, intent (inout) :: it  ! Pointer into x, y, z, f(*)

!     Execution:

      ierr = 0

      if (it == 0) then

!       Append the first point if there is room for it.

        if (it + 1 > ICUTDIM) then
          ierr = 5
        else
          it    = it + 1
          x(it) = xyzf(1,ndg)
          y(it) = xyzf(2,ndg)
          z(it) = xyzf(3,ndg)
          f(it) = xyzf(4,ndg)
        end if

      else

!       Append a point if there is enough room, unless it is within epsilon
!       of the previous point.

        if (abs (xyzf(1,ndg) - x(it)) <= epsxyz .and.
     >      abs (xyzf(2,ndg) - y(it)) <= epsxyz .and.
     >      abs (xyzf(3,ndg) - z(it)) <= epsxyz .and.
     >      abs (xyzf(4,ndg) - f(it)) <= epsf) then
        else
          if (it + 1 > ICUTDIM) then
            ierr = 5
          else
            it    = it + 1
            x(it) = xyzf(1,ndg)
            y(it) = xyzf(2,ndg)
            z(it) = xyzf(3,ndg)
            f(it) = xyzf(4,ndg)
          end if
        end if
      end if

      end subroutine cutordap

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine cutordin (itbeg, it, ndg)
!
!     Insert a point at x/y/z/f(itbeg).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in)    :: itbeg ! Pointer into x, y, z, f for new pt.
      integer, intent (inout) :: it    ! Pointer into x, y, z, f(*) next pt.
      integer, intent (in)    :: ndg   ! Pointer into xyzf(*,*)

!     Execution:

      ierr = 0

      if (it == 0) then

!       Insert the first point if there is room for it.

        if (it + 1 > ICUTDIM) then
          ierr = 5
        else
          it    = it + 1
          x(it) = xyzf(1,ndg)
          y(it) = xyzf(2,ndg)
          z(it) = xyzf(3,ndg)
          f(it) = xyzf(4,ndg)
        end if

      else

!       Insert a point if there is enough room, unless it is within epsilon
!       of the beginning point.

        if (abs (xyzf(1,ndg) - x(itbeg)) <= epsxyz .and.
     >      abs (xyzf(2,ndg) - y(itbeg)) <= epsxyz .and.
     >      abs (xyzf(3,ndg) - z(itbeg)) <= epsxyz .and.
     >      abs (xyzf(4,ndg) - f(itbeg)) <= epsf) then
        else
          if (it + 1 > ICUTDIM) then
            ierr = 5
          else
            it = it + 1

!           Shift the last curve by one place and insert the new point.

            x(itbeg+1:it) = x(itbeg:it-1)
            y(itbeg+1:it) = y(itbeg:it-1)
            z(itbeg+1:it) = z(itbeg:it-1)
            f(itbeg+1:it) = f(itbeg:it-1)

            x(itbeg) = xyzf(1,ndg)
            y(itbeg) = xyzf(2,ndg)
            z(itbeg) = xyzf(3,ndg)
            f(itbeg) = xyzf(4,ndg)
          end if
        end if
      end if

      end subroutine cutordin

      end subroutine cutorder
