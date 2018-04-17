!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine spherical_triangulation (n, tri_quadrant)
!
!  Purpose:
!
!     This utility triangulates a basic polar triangle defined by 3 equispaced
!  points on a unit-radius sphere: one at its north pole, two on its equator,
!  and (consequently) angles of pi/2 at each vertex.  Argument n specifies the
!  number of points along each edge of the basic triangle.  The intent is to
!  enable discretization of a hemisphere into sensibly-spaced points in order
!  to construct radiation lines of sight at point on a convex body.  The output
!  from this routine can be transformed to produce the four quadrants of some
!  discretized hemisphere at some orientation with some radius.  (Suppressing
!  duplicate edge nodes is an awkwardness left to the higher level ...)
!
!  Further details:
!
!     This is the spherical analogue of tessellating an equilateral triangle
!  in planar geometry with equal-sized equilateral triangular elements.  Note
!  that in the spherical case, though, the triangular elements are not all the
!  same: the elements at the vertices of the underlying basic polar triangle
!  necessarily have one pi/2 angle.  What can be said about the other elements?
!
!     The variation of the numbers of elements and points (nodes) is this:
!
!                   # nodes/edge   # elements   # nodes
!                        2             1            3
!                        3             4            6
!                        4             9           10
!                        5            16           15
!                        6            25           21
!                        :             :            :
!                        n        (n - 1)**2    n(n + 1)/2
!
!     Construction can be viewed as slicing the unit quarter sphere by planes
!  parallel to the equator that are uniformly spaced in latitude.  Starting
!  from the north pole, these latitude lines are divided into one more sector
!  with each lower latitude.
!
!  Conventions:
!
!     The x/y/z axes are right-handed with z pointing up through the north
!  pole of the underlying unit sphere:
!
!                                      z
!                                      |
!                                      |
!                                     / \
!                                   /     \
!                                 x         y
!
!     The nodes are constructed in the order shown below left for a planar
!  triangle, and the resulting elements are numbered correspondingly:
!
!                        Nodes                Elements
!                          1
!                                                 1
!                        2   3
!                                               2   3
!                      4   5   6
!                                             4   5   6
!                    7   8   9   10
!                                           7   8   9  10
!                  :   :   :   :   :              :
!
!     The storage convention for the x/y/z coordinates and the connectivity
!  information is that of an earlier I/O package for triangulated surfaces in
!  Tecplot format:
!
!     xyz(1:3,1:nnodes) are the coordinates of the grid points/nodes, with no
!                       repeated points;
!     conn(1:3,m)       are node indices into xyz(1:3,:) for triangular
!                       element m, with triangle vertices in clockwise order
!                       so that element normals all point outwards.
!
!  History:
!
!     Feb. 24, 2014  D.A.Saunders  Initial implementation, when it was belately
!                                  realized that working with lines of latitude
!                                  and longitude and structured surfaces is not
!                                  the right approach, as the nodes inevitably
!                                  cluster towards the north pole.
!     Feb. 25, 2014    "     "     Completed the awkward connectivity coding.
!
!  Author:  David Saunders, ERC, Inc./NASA Ames Research Center, Moffett Field
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use tri_zone_structure  ! See triangulation_io.f90
   use trigd

   implicit none

!  Arguments:

   integer, intent (in)          :: n             ! Requested number of
                                                  ! uniform edge points; n >= 2
   type (tri_type), intent (out) :: tri_quadrant  ! Triangulated quarter
                                                  ! hemisphere
!  Local constants:

   real, parameter :: ninety = 90., one = 1., zero = 0.

!  Local variables:

   integer :: i, j, n2, n3, ne, nelements, nl, nn, nnodes, np
   real    :: dphi, dtheta, phi, theta, r, z

!  Execution:

   nelements = (n - 1)**2;   tri_quadrant%nelements = nelements
   nnodes    = n*(n + 1)/2;  tri_quadrant%nnodes    = nnodes

   allocate (tri_quadrant%xyz(3,nnodes), tri_quadrant%conn(3,nelements))

!  Avoid computing the trig. functions more than once:

   dphi = ninety / real (n - 1)  ! Angular spacing of latitude lines

!  j = 1 (north pole):

   tri_quadrant%xyz(1,1) = zero
   tri_quadrant%xyz(2,1) = zero
   tri_quadrant%xyz(3,1) = one

   nn = 1  ! Node counter
   nl = 1  ! # longitude points on current latitude line

   do j = 2, n  ! Lines of latitude, starting nearest to the north pole
      phi = dphi * real (j-1)
      z   = cosd (phi)
      r   = sind (phi)
      dtheta = ninety / real (nl)
      nl  = nl + 1
      do i = 1, nl  ! Lines of longitude, starting at the x axis
         nn = nn + 1
         theta = dtheta * real (i-1)
         tri_quadrant%xyz(1,nn) = r * cosd (theta)
         tri_quadrant%xyz(2,nn) = r * sind (theta)
         tri_quadrant%xyz(3,nn) = z
      end do
   end do

!  Any discernible pattern in the numbering of the element vertices is not
!  easy to find!  If we color the top element black, and the first element of
!  each row black, and alternate the others in each row white and black, some
!  pattern emerges among the black cells and another among the white cells.
!  Once the first (black) cell of a row is established, the others follow.

   tri_quadrant%conn(1,1) = 1
   tri_quadrant%conn(2,1) = 2
   tri_quadrant%conn(3,1) = 3

   nn = 1  ! Node counter (= vertex 1 of each "black" cell)
   ne = 1  ! Element counter
   np = 0  ! # "white/black" cell pairs after the 1st black cell in current row

   do j = 2, n - 1   ! Rows of elements starting from row 2 near the north pole
      nn = nn + 1    ! Vertex 1  of 1st (black) element in this row
      ne = ne + 1    ! Element #  "   "   "   "   "   "   "   "   "
      np = np + 1    ! # white/black cell pairs for this j, after 1st black cell
      n2 = nn + np + 1   ! Vertex 2 of 1st element
      n3 = n2 + 1        ! Vertex 3 "   "    "
      tri_quadrant%conn(1,ne) = nn
      tri_quadrant%conn(2,ne) = n2
      tri_quadrant%conn(3,ne) = n3
      do i = 1, np  ! For each cell pair after black cell 1
         nn = nn + 1  ! Vertex 1  of current black cell
         ne = ne + 2  ! Element #  "   "   "   "   "
         tri_quadrant%conn(1,ne)   = nn      ! "Black" cell indices
         tri_quadrant%conn(2,ne)   = n3      ! Vertex 3 of previous black cell
         tri_quadrant%conn(3,ne)   = n3 + 1
         tri_quadrant%conn(1,ne-1) = n3      ! "White" cell indices
         tri_quadrant%conn(2,ne-1) = nn
         tri_quadrant%conn(3,ne-1) = nn - 1
         n3 = n3 + 1
      end do  ! Next cell pair
   end do  ! Next line of latitude

   end subroutine spherical_triangulation
