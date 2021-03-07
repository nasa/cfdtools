!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine process_hex_cell (ni, nj, nk, i, j, k, x, y, z, f, fx, fz,       &
                                xt, zt, npts, ypts, fpts)
!
!  Description:
!
!  Intersect the straight line parallel to the y axis through the point (xt, yt)
!  with the indicated hex cell of a structured grid block.  Interpolate vertex-
!  centered functions at the intersection points if any are found.  This number
!  is expected to be 2 or 0 in practice.
!
!  The intended application is to interpolating functions of refractive index n
!  at cell face/light ray intersections as needed for constructing shadowgraph/
!  schlieren-type images.
!
!  More details:
!
!  For the indicated structured grid cell in 3-space, intersect it with the
!  straight line perpendicular to the y = 0 plane through (xt, zt).  Each
!  cell face is triangulated with its shortest diagonal.  The intent is to
!  integrate the indicated function(s) along the line segments formed by the
!  points of intersection with the grid found by many calls for a given line
!  (or by multiple calls for a given cell, if trapezoidal integration suffices).
!
!  The number of intersections, npts, will normally be either 2 or 0, but
!  special cases are possible as follows (where "same" means to within a
!  small tolerance):
!
!     # intersections     interpretation
!
!     2 distinct          2 distinct triangles; use linearly interpolated fs
!     1                   shouldn't be possible
!     0                   line misses the cell
!     infinite            line is in the plane of a cell face _|_ y = 0
!     2+ not all same     edge point or corner point at one end; pick 1 pair
!
!  The "infinite" case is assumed not to occur in practice, since it can be
!  largely guarded against through the choice of the image plane resolution.
!  The reason is that there are too many pathological cases that would lead
!  to counting the same line segment in a cell face more than once if hex
!  cells are processed independently, as is the plan.
!
!  The calling program is expected to sort the ypts for each (xt, yt) line,
!  eliminate duplicate abscissas, and perform numerical quadrature for the
!  functions of interest.  The caller should also eliminate most grid cells
!  for a given line using bounding box techniques.  Actually, visiting each
!  grid cell only once is much more efficient if trapezoidal quadrature is
!  acceptable: simply increment the integration sums for the light rays found
!  to pass through the current cell; no sorting of line segments is necessary.
!
!  Jan. 2010  D. A. Saunders  Initial implementation for program SHADOWGRAPH.
!
!  Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)  :: ni, nj, nk   ! Structured grid block dimensions
   integer, intent (in)  :: i, j, k      ! "Lower left" indices of target cell
   real,    intent (in)  :: x(ni,nj,nk)  ! Grid coordinates
   real,    intent (in)  :: y(ni,nj,nk)
   real,    intent (in)  :: z(ni,nj,nk)
   real,    intent (in)  :: f(ni,nj,nk)  ! Associated functions are n (- 1) ...
   real,    intent (in)  :: fx(ni,nj,nk) ! ... and partial dn/dx ...
   real,    intent (in)  :: fz(ni,nj,nk) ! ... and partial dn/dz
   real,    intent (in)  :: xt, zt       ! Coordinates defining line _|_ y = 0
   integer, intent (out) :: npts         ! # line/cell intersns. found, 0 | 2
   real,    intent (out) :: ypts(2)      ! Corresp. ends of line segment if 2
   real,    intent (out) :: fpts(3,2)    ! Corresp. linearly interpolated fns.

!  Local variables:

   integer :: i1, i2, i3                 ! For triangle vertex indices
   integer :: istate                     ! For triangle/line intersection case
   integer :: ii, jj, kk                 ! For min or max index faces
   integer :: iface, itri                ! For loops over 6 faces/2 triangles
   integer :: iv(3,2)                    ! For the 3 vertices of 2 triangles
   real    :: xyz(3,4)                   ! For the 4 corners of a cell face
   real    :: fv(3,4)                    ! For 3 functions at 4 face vertices
   real    :: p, q, r, yt                ! Interpolation coefficients + yinterp

!  Procedures:

   external :: tri_line  ! Intersects straight line _|_ y = 0 with triangle

!  Execution:

   npts = 0

!  Traverse the six cell faces in the standard order in a way that allows
!  twelve triangles to be treated similarly:

   do iface = 1, 6

      select case (iface)

         case (1)  ! imin face

            ii = i
            xyz(1,1) =  x(ii,j,k)
            xyz(2,1) =  y(ii,j,k)
            xyz(3,1) =  z(ii,j,k)
            fv(1,1)  =  f(ii,j,k)
            fv(2,1)  = fx(ii,j,k)
            fv(3,1)  = fz(ii,j,k)
            xyz(1,2) =  x(ii,j+1,k)
            xyz(2,2) =  y(ii,j+1,k)
            xyz(3,2) =  z(ii,j+1,k)
            fv(1,2)  =  f(ii,j+1,k)
            fv(2,2)  = fx(ii,j+1,k)
            fv(3,2)  = fz(ii,j+1,k)
            xyz(1,3) =  x(ii,j,k+1)
            xyz(2,3) =  y(ii,j,k+1)
            xyz(3,3) =  z(ii,j,k+1)
            fv(1,3)  =  f(ii,j,k+1)
            fv(2,3)  = fx(ii,j,k+1)
            fv(3,3)  = fz(ii,j,k+1)
            xyz(1,4) =  x(ii,j+1,k+1)
            xyz(2,4) =  y(ii,j+1,k+1)
            xyz(3,4) =  z(ii,j+1,k+1)
            fv(1,4)  =  f(ii,j+1,k+1)
            fv(2,4)  = fx(ii,j+1,k+1)
            fv(3,4)  = fz(ii,j+1,k+1)

         case (2)  ! imax face

            ii = i + 1
            xyz(1,1) =  x(ii,j,k)
            xyz(2,1) =  y(ii,j,k)
            xyz(3,1) =  z(ii,j,k)
            fv(1,1)  =  f(ii,j,k)
            fv(2,1)  = fx(ii,j,k)
            fv(3,1)  = fz(ii,j,k)
            xyz(1,2) =  x(ii,j+1,k)
            xyz(2,2) =  y(ii,j+1,k)
            xyz(3,2) =  z(ii,j+1,k)
            fv(1,2)  =  f(ii,j+1,k)
            fv(2,2)  = fx(ii,j+1,k)
            fv(3,2)  = fz(ii,j+1,k)
            xyz(1,3) =  x(ii,j,k+1)
            xyz(2,3) =  y(ii,j,k+1)
            xyz(3,3) =  z(ii,j,k+1)
            fv(1,3)  =  f(ii,j,k+1)
            fv(2,3)  = fx(ii,j,k+1)
            fv(3,3)  = fz(ii,j,k+1)
            xyz(1,4) =  x(ii,j+1,k+1) 
            xyz(2,4) =  y(ii,j+1,k+1) 
            xyz(3,4) =  z(ii,j+1,k+1)
            fv(1,4)  =  f(ii,j+1,k+1)
            fv(2,4)  = fx(ii,j+1,k+1)
            fv(3,4)  = fz(ii,j+1,k+1)

         case (3)  ! jmin face

            jj = j
            xyz(1,1) =  x(i,jj,k)
            xyz(2,1) =  y(i,jj,k)
            xyz(3,1) =  z(i,jj,k)
            fv(1,1)  =  f(i,jj,k)
            fv(2,1)  = fx(i,jj,k)
            fv(3,1)  = fz(i,jj,k)
            xyz(1,2) =  x(i+1,jj,k)
            xyz(2,2) =  y(i+1,jj,k)
            xyz(3,2) =  z(i+1,jj,k)
            fv(1,2)  =  f(i+1,jj,k)
            fv(2,2)  = fx(i+1,jj,k)
            fv(3,2)  = fz(i+1,jj,k)
            xyz(1,3) =  x(i,jj,k+1)
            xyz(2,3) =  y(i,jj,k+1)
            xyz(3,3) =  z(i,jj,k+1)
            fv(1,3)  =  f(i,jj,k+1)
            fv(2,3)  = fx(i,jj,k+1)
            fv(3,3)  = fz(i,jj,k+1)
            xyz(1,4) =  x(i+1,jj,k+1)
            xyz(2,4) =  y(i+1,jj,k+1)
            xyz(3,4) =  z(i+1,jj,k+1)
            fv(1,4)  =  f(i+1,jj,k+1)
            fv(2,4)  = fx(i+1,jj,k+1)
            fv(3,4)  = fz(i+1,jj,k+1)

         case (4)  ! jmax face

            jj = j + 1
            xyz(1,1) =  x(i,jj,k)
            xyz(2,1) =  y(i,jj,k)
            xyz(3,1) =  z(i,jj,k)
            fv(1,1)  =  f(i,jj,k)
            fv(2,1)  = fx(i,jj,k)
            fv(3,1)  = fz(i,jj,k)
            xyz(1,2) =  x(i+1,jj,k)
            xyz(2,2) =  y(i+1,jj,k)
            xyz(3,2) =  z(i+1,jj,k)
            fv(1,2)  =  f(i+1,jj,k)
            fv(2,2)  = fx(i+1,jj,k)
            fv(3,2)  = fz(i+1,jj,k)
            xyz(1,3) =  x(i,jj,k+1)
            xyz(2,3) =  y(i,jj,k+1)
            xyz(3,3) =  z(i,jj,k+1)
            fv(1,3)  =  f(i,jj,k+1)
            fv(2,3)  = fx(i,jj,k+1)
            fv(3,3)  = fz(i,jj,k+1)
            xyz(1,4) =  x(i+1,jj,k+1)
            xyz(2,4) =  y(i+1,jj,k+1)
            xyz(3,4) =  z(i+1,jj,k+1)
            fv(1,4)  =  f(i+1,jj,k+1)
            fv(2,4)  = fx(i+1,jj,k+1)
            fv(3,4)  = fz(i+1,jj,k+1)

         case (5)  ! kmin face

            kk = k
            xyz(1,1) =  x(i,j,kk)
            xyz(2,1) =  y(i,j,kk)
            xyz(3,1) =  z(i,j,kk)
            fv(1,1)  =  f(i,j,kk)
            fv(2,1)  = fx(i,j,kk)
            fv(3,1)  = fz(i,j,kk)
            xyz(1,2) =  x(i+1,j,kk)
            xyz(2,2) =  y(i+1,j,kk)
            xyz(3,2) =  z(i+1,j,kk)
            fv(1,2)  =  f(i+1,j,kk)
            fv(2,2)  = fx(i+1,j,kk)
            fv(3,2)  = fz(i+1,j,kk)
            xyz(1,3) =  x(i,j+1,kk)
            xyz(2,3) =  y(i,j+1,kk)
            xyz(3,3) =  z(i,j+1,kk)
            fv(1,3)  =  f(i,j+1,kk)
            fv(2,3)  = fx(i,j+1,kk)
            fv(3,3)  = fz(i,j+1,kk)
            xyz(1,4) =  x(i+1,j+1,kk)
            xyz(2,4) =  y(i+1,j+1,kk)
            xyz(3,4) =  z(i+1,j+1,kk)
            fv(1,4)  =  f(i+1,j+1,kk)
            fv(2,4)  = fx(i+1,j+1,kk)
            fv(3,4)  = fz(i+1,j+1,kk)

         case (6)  ! kmax face

            kk = k + 1
            xyz(1,1) =  x(i,j,kk)
            xyz(2,1) =  y(i,j,kk)
            xyz(3,1) =  z(i,j,kk)
            fv(1,1)  =  f(i,j,kk)
            fv(2,1)  = fx(i,j,kk)
            fv(3,1)  = fz(i,j,kk)
            xyz(1,2) =  x(i+1,j,kk)
            xyz(2,2) =  y(i+1,j,kk)
            xyz(3,2) =  z(i+1,j,kk)
            fv(1,2)  =  f(i+1,j,kk)
            fv(2,2)  = fx(i+1,j,kk)
            fv(3,2)  = fz(i+1,j,kk)
            xyz(1,3) =  x(i,j+1,kk)
            xyz(2,3) =  y(i,j+1,kk)
            xyz(3,3) =  z(i,j+1,kk)
            fv(1,3)  =  f(i,j+1,kk)
            fv(2,3)  = fx(i,j+1,kk)
            fv(3,3)  = fz(i,j+1,kk)
            xyz(1,4) =  x(i+1,j+1,kk)
            xyz(2,4) =  y(i+1,j+1,kk)
            xyz(3,4) =  z(i+1,j+1,kk)
            fv(1,4)  =  f(i+1,j+1,kk)
            fv(2,4)  = fx(i+1,j+1,kk)
            fv(3,4)  = fz(i+1,j+1,kk)

      end select

!     If the target line is outside the face bounding box, skip the face.
!     Testing here avoids storing cell bounding boxes for skipping cells,
!     which would take more memory than the grid itself.

      if (xt > max (xyz(1,1), xyz(1,2), xyz(1,3), xyz(1,4))) cycle
      if (xt < min (xyz(1,1), xyz(1,2), xyz(1,3), xyz(1,4))) cycle
      if (zt > max (xyz(3,1), xyz(3,2), xyz(3,3), xyz(3,4))) cycle
      if (zt < min (xyz(3,1), xyz(3,2), xyz(3,3), xyz(3,4))) cycle

!     Triangulate the quad face using the shorter diagonal:

      call triangulate_quad_face (xyz, iv)

!     Intersect each triangle of this hex face with the target line:

      do itri = 1, 2

         i1 = iv(1,itri);  i2 = iv(2,itri);  i3 = iv(3,itri)

         call tri_line (xyz(1,i1), xyz(1,i2), xyz(1,i3), &
                        xt, yt, zt, p, q, r, istate)

!!!!     write (11, '(a, 4i5, 1p, 3e15.6)') 'istate, i, j, k, p, q, r', &
!!!!                                         istate, i, j, k, p, q, r

         if (istate == 1) then  ! Intersection found

            npts = npts + 1
            ypts(npts) = yt

            fpts(:,npts) = p*fv(:,i1) + q*fv(:,i2) + r*fv(:,i3)

            exit  ! Avoid possible duplicate intersection on common diagonal.
                  ! Note that if the quad is not planar, this could miss two
                  ! distinct intersections in the same quad face.  There is
                  ! no good answer, but the cell or its neighbor at that face
                  ! would then be nonconvex.  At a convex wall, we don't want
                  ! to count such a possibility anyway.
         end if

      end do  ! Next triangle of this face

      if (npts == 2) exit

   end do  ! Next hex cell face

   end subroutine process_hex_cell
