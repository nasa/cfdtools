!*****************************************************************************80
subroutine points_delaunay_naive_2d ( node_num, node_xy, maxtri, &
  element_num, element_node )
!
!! POINTS_DELAUNAY_NAIVE_2D is a naive Delaunay triangulation scheme.
!
!  Discussion:
!
!    This routine is only suitable as a demonstration code for small
!    problems.  Its running time is of order NODE_NUM**4.  Much faster
!    algorithms are available.
!
!    Given a set of nodes in the plane, a triangulation is set of
!    triples of distinct nodes, forming triangles, so that every
!    point within the convex hull of the set of nodes is either
!    one of the nodes, or lies on an edge of one or more triangles,
!    or lies within exactly one triangle.
!
!    A Delaunay triangulation is a triangulation with additional
!    properties.
!
!    NODE_NUM must be at least 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Joseph ORourke,
!    Computational Geometry,
!    Cambridge University Press,
!    Second Edition, 1998, page 187.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of the nodes.
!
!    Input, integer ( kind = 4 ) MAXTRI, the maximum number of triangles.
!
!    Output, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles in
!    the triangulation.
!
!    Output, integer ( kind = 4 ) ELEMENT_NODE(3,MAXTRI), the indices of
!    the triangle nodes.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ) maxtri
  integer ( kind = 4 ) node_num

  logical flag
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  real ( kind = 8 ) node_xy(dim_num,node_num)
  integer ( kind = 4 ) element_node(3,maxtri)
  integer ( kind = 4 ) element_num
  real ( kind = 8 ) xn
  real ( kind = 8 ) yn
  real ( kind = 8 ) z(node_num)
  real ( kind = 8 ) zn

  element_num = 0

  if ( node_num < 3 ) then
    return
  end if
!
!  Compute Z = X*X + Y*Y.
!
  z(1:node_num) = node_xy(1,1:node_num)**2 + node_xy(2,1:node_num)**2
!
!  For each triple (I,J,K):
!
  do i = 1, node_num - 2
    do j = i+1, node_num
      do k = i+1, node_num

        if ( j /= k ) then

          xn = ( node_xy(2,j) - node_xy(2,i) ) * ( z(k) - z(i) ) &
             - ( node_xy(2,k) - node_xy(2,i) ) * ( z(j) - z(i) )

          yn = ( node_xy(1,k) - node_xy(1,i) ) * ( z(j) - z(i) ) &
             - ( node_xy(1,j) - node_xy(1,i) ) * ( z(k) - z(i) )

          zn = ( node_xy(1,j) - node_xy(1,i) ) &
             * ( node_xy(2,k) - node_xy(2,i) ) &
             - ( node_xy(1,k) - node_xy(1,i) ) &
             * ( node_xy(2,j) - node_xy(2,i) )

          flag = ( zn < 0.0D+00 )

          if ( flag ) then
            do m = 1, node_num
              flag = flag .and. &
                ( ( node_xy(1,m) - node_xy(1,i) ) * xn &
                + ( node_xy(2,m) - node_xy(2,i) ) * yn &
                + ( z(m)   - z(i) )   * zn <= 0.0D+00 )
            end do
          end if

          if ( flag ) then
            if ( element_num < maxtri ) then
              element_num = element_num + 1
              element_node(1:3,element_num) = (/ i, j, k /)
            end if
          end if

        end if

      end do
    end do
  end do

  return
end
