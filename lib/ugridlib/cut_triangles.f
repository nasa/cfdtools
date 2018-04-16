!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine cut_triangles (npts, ntris, iv, xyz, f, xyzfcut,
     >                          MNODE, MEDGE, nnode, ndgcut, nedge,
     >                          p, s, epsxyz, epsf, ierr)
!
!  Description:
!
!     Cut_triangles cuts a triangulated surface mesh with a specified plane.
!     A single scalar function is interpolated along with the Cartesian
!     coordinates representing the cut.  More precisely, cut_triangles
!     finds oriented line segments, or edges, which lie on the intersection
!     of an unstructured surface patch with a scalar function on it
!     (x,y,z,f) and a plane specified by a point and a surface normal
!     vector (p,s).  Edges are stored as points (xyzfcut) and oriented
!     edge numbers (ndgcut).  Cut_triangles may be called repeatedly to
!     accumulate edges if nnode and nedge are left alone between calls.
!
!     If two points are sufficiently close together then the first one is
!     taken and the second one is ignored.  Point p1 = (x1,y1,z1,f1) is
!     considered to be close to p2 = (x2,y2,z2,f2) if and only if:
!
!        |x1-x2| <= epsxyz  and  |y1-y2| <= epsxyz  and
!        |z1-z2| <= epsxyz  and  |f1-f2| <= epsf
!
!     The caller must set nedge = 0 to start a new cut, and should usually
!     set nnode = 0 at the same time.  The edges may be ordered into a 
!     collection of unbroken polygonal curves by using subroutine cutorder.
!
!  History:
!
!     01-APR-1998  S.Thomas    Initial implementation of cutsurf for quads.
!     09-APR-1998     "        Companion cutorder eases plotting of cuts.
!     18-AUG-1998     "        Initial Fortran 90 translation.
!     20-AUG-1998  D.Saunders  Cutsurf & cutorder refined as library utilities.
!     22-DEC-2001     "        Cut_triangles adapted from cutsurf.  It was
!                              simpler to retain f than to remove it.
!
!  Authors:
!
!     Scott D. Thomas, Sterling Software, Inc.,  Contract NAS2-13210/6/5
!                      Raytheon STX Corporation, Contract NAS2-98080/53
!                      NASA Ames Research Center, Moffett Field, CA.
!
!     David Saunders,  ELORET, NASA Ames.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

!  Arguments:

      integer, intent (in) ::
     >  npts,                    ! Number of points in xyz(1:3,*) and f(*)
     >  ntris,                   ! Number of triangles defined by iv(1:3,*)
     >  iv(3,ntris),             ! Pointers to triangle vertices in xyz(1:3,*)
     >  MNODE,                   ! Max. # end points of line segments handled
     >                           ! per cut; MNODE = 2 x MEDGE is safe
     >  MEDGE                    ! Max. # 2-pt. line segments handled per cut;
                                 ! MEDGE = ntris is an upper bound
      real, intent (in) ::
     >  xyz(3,npts),             ! Unstructured surface mesh coordinates and ...
     >  f(npts)                  ! ... associated scalar surface function

      real, intent (in) ::
     >  p(3),                    ! Point p(1:3) and vector s(1:3) from p ...
     >  s(3),                    ! ... define the cutting plane.  E.g., for a
     >                           ! cut at z = Z, use p = (0,0,Z) & s = (0,0,1.)
     >  epsxyz, epsf             ! Tolerances for coincident pts. such as the
                                 ! data ranges x 1.E-6; see usage above

      real, intent (inout) ::
     >  xyzfcut(4,MNODE)         ! Packed (x,y,z,f) coordinates of end points
                                 ! of 2-pt. line segments for this cut, with
                                 ! no coincident points according to epsxyz/f
      integer, intent (inout) ::
     >  ndgcut(2,MEDGE),         ! Pointers defining the ends of each line
     >                           ! segment as packed in xyzfcut(*,1:nnode)
     >  nnode,                   ! Current number of points on the cut and ...
     >  nedge,                   ! ... corresp. no. of edges (2-pt. line segs.)
     >  ierr                     !  0 = no error detected;
                                 ! -1 = surface-defining vector s = (0.,0.,0.);
                                 !  1 = MNODE limit encountered;
                                 !  2 = MEDGE   "     " ;
                                 !  3 = nnode < 0 on input;
                                 !  4 = nedge < 0  "  "
!     Local constants:

      real, parameter :: one = 1., zero = 0.

!     Local variables:

      integer
     >   i, j
      real
     >   a(4), b(4), c(4)

!     Execution:

!     S, the vector defining the cutting surface, must be positive.

      if (s(1)**2 + s(2)**2 + s(3)**2 <= zero) then
         ierr = -1
         go to 999
      end if

!     Treat each triangle:

      do i = 1, ntris

         j    = iv(1,i)
         a(1) = xyz(1,j)
         a(2) = xyz(2,j)
         a(3) = xyz(3,j)
         a(4) = f(j)

         j    = iv(2,i)
         b(1) = xyz(1,j)
         b(2) = xyz(2,j)
         b(3) = xyz(3,j)
         b(4) = f(j)

         j    = iv(3,i)
         c(1) = xyz(1,j)
         c(2) = xyz(2,j)
         c(3) = xyz(3,j)
         c(4) = f(j)

         call cutsurft (a, b, c)

         if (ierr /= 0) go to 999

      end do

  999 return

!  Internal procedures for cut_triangles:
!  cutsurft, cutsurf1, cutsurf2, cutsurf3

      contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine cutsurft (a, b, c)
!
!     Cut a surface triangle abc by the plane defined by a point P and a
!     direction vector S.  Put the resulting line segment, if any, into
!     xyzfcut, ndgcut but do not exceed MNODE, MEDGE.  Return ierr = 0
!     if abc does not touch the plane or if one or more line segments
!     are found in the cutting plane, else ierr /= 0 is as described in
!     cutsurf.  E.g., if another line segment would make nedge > MEDGE
!     then ierr is set to 2 followed by an immediate return.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      real, intent (in) :: a(4), b(4), c(4)

!     Local variables:

      integer numzer

      real da, db, dc, dmin, dmax, q(3)

!     Execution:

!     The inner product determines the signed distance of each point
!     from the plane.
!
      da = (a(1)-p(1))*s(1) + (a(2)-p(2))*s(2) + (a(3)-p(3))*s(3)
      db = (b(1)-p(1))*s(1) + (b(2)-p(2))*s(2) + (b(3)-p(3))*s(3)
      dc = (c(1)-p(1))*s(1) + (c(2)-p(2))*s(2) + (c(3)-p(3))*s(3)

      dmin = min (da, db, dc)
      dmax = max (da, db, dc)

      if (dmax < zero .or. dmin > zero) then

!       The triangle does not touch or cross the plane, so bail out.

        ierr = 0
        go to 999

      end if

      numzer = 0
      if (da == zero) numzer = numzer + 1
      if (db == zero) numzer = numzer + 1
      if (dc == zero) numzer = numzer + 1

!     Q = (b-a) x (c-a) is a vector normal to triangle abc which will be
!     used to assure that the orientation of each line segment matches the
!     orientation of the triangle, if possible.

      q(1) = (b(2)-a(2))*(c(3)-a(3)) - (b(3)-a(3))*(c(2)-a(2))
      q(2) = (b(3)-a(3))*(c(1)-a(1)) - (b(1)-a(1))*(c(3)-a(3))
      q(3) = (b(1)-a(1))*(c(2)-a(2)) - (b(2)-a(2))*(c(1)-a(1))

      if (numzer >= 2) then

!       One or more edges is/are on the plane; usually either one or three.

        if (da == zero .and. db == zero) then

          call cutsurf1 (a, b, q, s)

          if (ierr /= 0) go to 999
        end if

        if (db == zero .and. dc == zero) then

          call cutsurf1 (b, c, q, s)

          if (ierr /= 0) go to 999
        end if

        if (dc == zero .and. da == zero) then

          call cutsurf1 (c, a, q, s)

          if (ierr /= 0) go to 999
        end if

      else if (numzer == 1) then

!       One point is on the plane; abc may or may not cross the plane.

        if (sign (one, dmin) * dmax >= zero) then

!         The opposite edge does not cross the plane, so get a
!         zero-length line segment on the point that touches the plane.

          if (da == zero) then

            call cutsurf1 (a, a, q, s)

          else if (db == zero) then

            call cutsurf1 (b, b, q, s)

          else

            call cutsurf1 (c, c, q, s)

          end if

        else ! The opposite edge crosses the plane.

          if (da == zero) then

            call cutsurf2 (a, b, c, db, dc, q, s)

          else if (db == zero) then

            call cutsurf2 (b, c, a, dc, da, q, s)

          else

            call cutsurf2 (c, a, b, da, db, q, s)

          end if

        end if

      else

!       No points are on the plane so the triangle crosses at two edges.

        if (sign (one, db) * dc > zero) then      ! a & bc are on opposite sides

          call cutsurf3 (c, a, b, dc, da, db, q, s)

        else if (sign (one, dc) * da > zero) then ! b & ca  "   "   "

          call cutsurf3 (a, b, c, da, db, dc, q, s)

        else                                      ! c & ab  "   "   "

          call cutsurf3 (b, c, a, db, dc, da, q, s)

        end if

      end if

  999 return

      end subroutine cutsurft

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine cutsurf1 (a, b, q, s)
!
!     Get the line segment from point a to point b.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      real, intent (in) :: a(4), b(4), q(3), s(3)

!     Local variables:

      real det

!     Execution:

      det = s(1) * (q(2)*(b(3)-a(3)) - q(3)*(b(2)-a(2))) +
     >      s(2) * (q(3)*(b(1)-a(1)) - q(1)*(b(3)-a(3))) +
     >      s(3) * (q(1)*(b(2)-a(2)) - q(2)*(b(1)-a(1)))

      if (det < zero) then

        call cutsurfa (a, b)

      else

        call cutsurfa (b, a)

      end if

      end subroutine cutsurf1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine cutsurf2 (a, b, c, db, dc, q, s)
!
!     Get the line segment from point a to an interpolated point on bc.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      real, intent (in) :: a(4), b(4), c(4), db, dc, q(3), s(3)

!     Local variables:

      real d(4), det, sigma

!     Execution:

      sigma = abs (db) / (abs (db) + abs (dc))

      d(1)  = (one - sigma) * b(1) + sigma * c(1)
      d(2)  = (one - sigma) * b(2) + sigma * c(2)
      d(3)  = (one - sigma) * b(3) + sigma * c(3)
      d(4)  = (one - sigma) * b(4) + sigma * c(4)

      det   = s(1) * (q(2)*(d(3)-a(3)) - q(3)*(d(2)-a(2))) +
     >        s(2) * (q(3)*(d(1)-a(1)) - q(1)*(d(3)-a(3))) +
     >        s(3) * (q(1)*(d(2)-a(2)) - q(2)*(d(1)-a(1)))

      if (det < zero) then

        call cutsurfa (a, d)

      else

        call cutsurfa (d, a)

      end if

      end subroutine cutsurf2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine cutsurf3 (a, b, c, da, db, dc, q, s)
!
!  Get the line segment from a point on ab to a point on bc.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      real, intent (in) :: a(4), b(4), c(4), da, db, dc, q(3), s(3)

!     Local variables:

      real d(4), e(4), det, sigma

!     Execution:

      sigma = abs (da) / (abs (da) + abs (db))
      d(1)  = (one - sigma) * a(1) + sigma * b(1)
      d(2)  = (one - sigma) * a(2) + sigma * b(2)
      d(3)  = (one - sigma) * a(3) + sigma * b(3)
      d(4)  = (one - sigma) * a(4) + sigma * b(4)

      sigma = abs (db) / (abs (db) + abs (dc))
      e(1)  = (one - sigma) * b(1) + sigma * c(1)
      e(2)  = (one - sigma) * b(2) + sigma * c(2)
      e(3)  = (one - sigma) * b(3) + sigma * c(3)
      e(4)  = (one - sigma) * b(4) + sigma * c(4)

      det   = s(1) * (q(2)*(d(3)-e(3)) - q(3)*(d(2)-e(2))) +
     >        s(2) * (q(3)*(d(1)-e(1)) - q(1)*(d(3)-e(3))) +
     >        s(3) * (q(1)*(d(2)-e(2)) - q(2)*(d(1)-e(1)))

      if (det < zero) then

        call cutsurfa (e, d)

      else

        call cutsurfa (d, e)

      end if

      end subroutine cutsurf3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine cutsurfa (a, b)
!
!     Add segment ab to xyzfcut and ndgcut unless it is already there.
!     A second point within epsilon of another point will be ignored.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      real, intent (in) :: a(4), b(4)

!     Local variables:

      integer n, n1, n2

!     Execution:

      if (nedge < 0) then
        ierr = 4
      else if (nnode < 0) then
        ierr = 3
      else if (nedge == 0) then
        if (nedge + 1 > MEDGE) then
          ierr = 2
        else if (nnode + 2 > MNODE) then
          ierr = 1
        else
          ierr = 0

!         Add the first line segment.

          nnode = nnode + 2
          xyzfcut(1:4,nnode-1) = a
          xyzfcut(1:4,nnode)   = b

          nedge = nedge + 1
          ndgcut(1,nedge) = nnode-1
          ndgcut(2,nedge) = nnode
        end if

      else ! Current nedge > 0

!       Add the line segment if it is not there already.

        do n = 1, nedge

          n1 = ndgcut(1,n)
          n2 = ndgcut(2,n)

          if (abs (xyzfcut(1,n1) - a(1)) <= epsxyz .and.
     >        abs (xyzfcut(2,n1) - a(2)) <= epsxyz .and.
     >        abs (xyzfcut(3,n1) - a(3)) <= epsxyz .and.
     >        abs (xyzfcut(4,n1) - a(4)) <= epsf   .and.
     >        abs (xyzfcut(1,n2) - b(1)) <= epsxyz .and.
     >        abs (xyzfcut(2,n2) - b(2)) <= epsxyz .and.
     >        abs (xyzfcut(3,n2) - b(3)) <= epsxyz .and.
     >        abs (xyzfcut(4,n2) - b(4)) <= epsf) then

!           The line segment is already there: don't add it again - bail out.
!
            ierr = 0
            go to 999

          end if

        end do

!       It's not there already so add the line segment if it fits.

        if (nedge + 1 > MEDGE) then
          ierr = 2
        else
          nedge = nedge + 1

!         Add the first point if it is not there already.
!
          do n = 1, nnode
            if (abs (xyzfcut(1,n) - a(1)) <= epsxyz .and.
     >          abs (xyzfcut(2,n) - a(2)) <= epsxyz .and.
     >          abs (xyzfcut(3,n) - a(3)) <= epsxyz .and.
     >          abs (xyzfcut(4,n) - a(4)) <= epsf) then
              ndgcut(1,nedge) = n
              go to 10
            end if
          end do

          if (nnode + 1 > MNODE) then
            ierr = 1
            go to 999
          else
            nnode = nnode + 1
            xyzfcut(1:4,nnode) = a
            ndgcut(1,nedge) = nnode
          end if

   10     continue

!         Add the second point if it is not there already.

          do n = 1, nnode
            if (abs (xyzfcut(1,n) - b(1)) <= epsxyz .and.
     >          abs (xyzfcut(2,n) - b(2)) <= epsxyz .and.
     >          abs (xyzfcut(3,n) - b(3)) <= epsxyz .and.
     >          abs (xyzfcut(4,n) - b(4)) <= epsf) then
              ndgcut(2,nedge) = n
              go to 20
            end if
          end do

          if (nnode + 1 > MNODE) then
            ierr = 1
            go to 999
          else
            nnode = nnode + 1
            xyzfcut(1:4,nnode) = b
            ndgcut(2,nedge) = nnode
          end if

   20     continue
          ierr = 0

        end if

      end if

  999 return

      end subroutine cutsurfa

      end subroutine cut_triangles
