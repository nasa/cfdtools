!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine process_pyramid (nnode, inode, x, y, z, f, fx, fz, xt, zt,       &
                               npts, ypts, fpts)
!
!  Description:
!
!  Intersect the straight line parallel to the y axis through the point (xt, yt)
!  with the indicated pyramid cell of an unstructured grid.  Interpolate vertex-
!  centered functions at the intersection points if any are found.  This number
!  is expected to be 2 or 0 in practice.
!
!  The intended application is to interpolating functions of refractive index n
!  at cell face/light ray intersections as needed for constructing shadowgraph/
!  schlieren-type images.
!
!  More details (adapted from the original process_hex_cell variant):
!
!  For the indicated unstructured grid cell in 3-space, intersect it with the
!  straight line perpendicular to the y = 0 plane through (xt, zt).  Each cell
!  face is intersected until 2 intersections are found, if any.  The intent is
!  to integrate the indicated function(s) along the line segments formed by the
!  points of intersection with the grid found by many calls for a given line
!  (or by multiple calls for a given cell, if trapezoidal integration suffices).
!
!  The vertex numbering is assumed to follow the FieldView convention, with the
!  quadrilateral face of the pyramid defined by vertices 1, 2, 3, 4.
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
!     2+ not all same     edge point or vertex point at one end; pick 1 pair
!
!  The "infinite" case is assumed not to occur in practice, since it can be
!  largely guarded against through the choice of the image plane resolution.
!  The reason is that there are too many pathological cases that would lead
!  to counting the same line segment in a cell face more than once if grid
!  cells are processed independently, as is the plan.
!
!  See the preceding variant process_hex_cell for the structured grid case.
!
!  Jan. 2010  D. A. Saunders  Initial implementation of process_hex_cell.
!  May  2010     "      "     Adapted as process_pyramid for FUN3D grids.
!
!  Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)  :: nnode        ! # points/nodes in this grid
   integer, intent (in)  :: inode(5)     ! Indices defining this pyramid cell
   real,    intent (in)  :: x(nnode)     ! Grid point coordinates
   real,    intent (in)  :: y(nnode)   
   real,    intent (in)  :: z(nnode)   
   real,    intent (in)  :: f(nnode)     ! Associated functions are n (- 1) ...
   real,    intent (in)  :: fx(nnode)    ! ... and partial dn/dx ...
   real,    intent (in)  :: fz(nnode)    ! ... and partial dn/dz
   real,    intent (in)  :: xt, zt       ! Coordinates defining line _|_ y = 0
   integer, intent (out) :: npts         ! # line/cell intersns. found, 0 | 2
   real,    intent (out) :: ypts(2)      ! Corresp. ends of line segment if 2
   real,    intent (out) :: fpts(3,2)    ! Corresp. linearly interpolated fns.

!  Local variables:

   integer :: i, in                      ! Vertex index & node index
   integer :: i1, i2, i3                 ! Vertex indices
   integer :: iface                      ! For loop over 4 triangular faces
   integer :: istate                     ! For triangle/line intersection case
   integer :: itri                       ! For loop over triangles of quad face
   integer :: iv(3)                      ! For the 3 vertices of 1 triangle
   integer :: iv2(3,2)                   ! For triangulating the quad face
   real    :: xyz(3,4)                   ! For 3 coordinates at 3 or 4 vertices
   real    :: fv(3,4)                    ! For 3 functions   at 3 or 4 vertices
   real    :: p, q, r, yt                ! Interpolation coefficients + yinterp

!  Procedures:

   external :: tri_line  ! Intersects straight line _|_ y = 0 with triangle

!  Execution:

!  Process the cell faces in the FieldView Reference Manual order:

!                             Face   Nodes
!                             1      1 2 3 4  ! quad face
!                             2      2 3 5
!                             3      3 4 5
!                             4      4 5 1
!                             5      5 1 2
   npts = 0

!  Quadrilateral face:

   do i = 1, 4
      in       = inode(i)
      xyz(1,i) = x(in)
      xyz(2,i) = y(in)
      xyz(3,i) = z(in)
      fv(1,i)  = f(in)
      fv(2,i)  = fx(in)
      fv(3,i)  = fz(in)
   end do

!  If the target line is outside the face bounding box, skip the face.
!  Testing here avoids storing cell bounding boxes for skipping cells,
!  which would take more memory than the grid itself.

   if (xt > max (xyz(1,1), xyz(1,2), xyz(1,3), xyz(1,4))) go to 50  ! Pragmatic
   if (xt < min (xyz(1,1), xyz(1,2), xyz(1,3), xyz(1,4))) go to 50
   if (zt > max (xyz(3,1), xyz(3,2), xyz(3,3), xyz(3,4))) go to 50
   if (zt < min (xyz(3,1), xyz(3,2), xyz(3,3), xyz(3,4))) go to 50

!  Triangulate the quad using its shorter diagonal:

   call triangulate_quad_face (xyz, iv2)

!  Intersect each triangle of this quad face with the target line:

   do itri = 1, 2

      i1 = iv2(1,itri);  i2 = iv2(2,itri);  i3 = iv2(3,itri)

      call tri_line (xyz(1,i1), xyz(1,i2), xyz(1,i3), &
                     xt, yt, zt, p, q, r, istate)

      if (istate == 1) then  ! Intersection found
         npts = npts + 1
         ypts(npts) = yt
         fpts(:,npts) = p*fv(:,i1) + q*fv(:,i2) + r*fv(:,i3)
         exit  ! Avoid possible duplicate intersection on common diagonal
      end if

   end do  ! Next triangle of this face

50 continue

!  For the triangle faces of the pyramid, follow the process_tet_cell approach:

   do iface = 2, 5

      iv(1) = iface
      iv(2) = iv(1) + 1;  if (iv(2) > 5) iv(2) = 1
      iv(3) = iv(2) + 1;  if (iv(3) > 5) iv(3) = 1
      if (iface == 2) iv(3) = 5  ! Special awkward case (see above)

      do i = 1, 3
         in       = inode(iv(i))
         xyz(1,i) =  x(in)
         xyz(2,i) =  y(in)
         xyz(3,i) =  z(in)
         fv(1,i)  =  f(in)
         fv(2,i)  = fx(in)
         fv(3,i)  = fz(in)
      end do

!     If the target line is outside the face bounding box, skip the face.

      if (xt > max (xyz(1,1), xyz(1,2), xyz(1,3))) cycle
      if (xt < min (xyz(1,1), xyz(1,2), xyz(1,3))) cycle
      if (zt > max (xyz(3,1), xyz(3,2), xyz(3,3))) cycle
      if (zt < min (xyz(3,1), xyz(3,2), xyz(3,3))) cycle

!     Intersect this triangular face with the target line:

      call tri_line (xyz(1,1), xyz(1,2), xyz(1,3), xt, yt, zt, p, q, r, istate)

      if (istate == 1) then  ! Intersection found
         npts = npts + 1
         ypts(npts) = yt

         fpts(:,npts) = p*fv(:,1) + q*fv(:,2) + r*fv(:,3)

         if (npts == 2) exit
      end if

   end do  ! Next triangular face

99 return

   end subroutine process_pyramid
