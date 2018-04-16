!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine nonlinear_tri_interp (ndegree, nf, nnear, nnode, ntri, conn,     &
                                    xnode, fnode, itri, p, q, r, finterp)

!  Purpose:
!
!     Interpolate one or more function values at the target point defined by
!  linear interpolation coefficients p, q, r within the triangle itri of a
!  surface triangulation in 3-space, using higher order corrections to the
!  initial linear interpolations.
!
!  Outline:
!
!     The linear method can be found in earlier routines such as TRI_COEFS and
!  NEAREST_TRI_COEFS (which adjusts TRI_COEFS-type coefficients to be not out-
!  side the current triangle).  The linear coefficients should be input here,
!  and the target (x,y,z) is actually redundant.
!
!     Nonlinear corrections are obtained using linear least squares techniques
!  on data from points in the triangulation outside the triangle itri already
!  identified for linear interpolation.  This calculation minimizes a measure
!  of the combined differences of the interpolated function from the functions
!  at the neighboring data points.
!
!     The number of nearest neighbor points to use is a variable, but it should
!  be large enough to avoid possible degeneracy.  For small datasets, use of all
!  points in every interpolation may be appropriate: just input nnear >= nnode.
!  This would not be appropriate for very large triangulations, where location
!  of nearest neighbors becomes an issue (not implemented yet).
!
!  Reference:  "Interpolation from a cloud of points" by Timothy J. Baker
!              Proceedings, 12the International Meshing Roundtable,
!              Sandia National Laboratories, pp. 55-63, September 2003
!
!  History:
!
!     02/04/05  DAS  Initial implementation, with about a dozen triangulation
!                    points in mind.
!
!  Author:  David Saunders, ELORET/NASA Ames Research Center, Moffett Field, CA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in) :: ndegree ! 2 => quadratic corrections (m x 3 systems);
                                   ! 3 => cubic     corrections (m x 7 systems)
   integer, intent (in) :: nf      ! 1 or more functions to be interpolated at
                                   ! at the indicated target point
   integer, intent (in) :: nnear   ! # neighboring points outside the current
                                   ! triangle to include in the corrections
   integer, intent (in) :: nnode   ! # data points in the triangulation

   integer, intent (in) :: ntri    ! # triangles defined by conn(1:3, 1:ntri)

   integer, intent (in) :: conn(3,ntri)    ! Pointers to triangle vertices

   real,    intent (in) :: xnode(3,nnode)  ! (x,y,z)s pointed to by conn(1:3,:)

   real,    intent (in) :: fnode(nf,nnode) ! Corresponding function data

   integer, intent (in) :: itri    ! Triangle already located as the nearest to
                                   ! the current target, xinterp(:)
   real,    intent (in) :: p, q, r ! Corresponding linear interpolation coefs.

   real,    intent (out) :: finterp(nf)    ! Desired interpolated functions

!  Local variables:

   integer :: i, in, l, m, n, n1, n2, n3

   real    :: dsq, flinear, pi, qi, ri, residual, x0(3)

   real, allocatable :: A(:,:), b(:,:), u(:)

!  Execution:

   m = min (nnear, nnode - 3)  ! # neigbors used in the nonlinear corrections
   n = 4 * ndegree - 5         ! 3 for quadratic, 7 for cubic corrections

   n1 = conn(1,itri)
   n2 = conn(2,itri)
   n3 = conn(3,itri)

   allocate (A(m,n), b(m,nf), u(m))

!  Using data points other than those of the current triangle, set up the LHS:

   i = 0

   do in = 1, nnode

      if (in == n1) cycle
      if (in == n2) cycle
      if (in == n3) cycle

!     Interpolation coefs. for this node within the plane of triangle itri:

      call tri_coefs (xnode(1,n1), xnode(1,n2), xnode(1,n3), xnode(1,in), x0,  &
                      dsq, pi, qi, ri)
      i = i + 1

      select case (n) ! For LHS matrix

         case (3) ! Quadratic

            A(i,1) = pi * qi;  A(i,2) = qi * ri;  A(i,3) = ri * pi

         case (7) ! Cubic

            A(i,1) = pi*qi*qi;  A(i,2) = pi*ri*ri;  A(i,3) = qi*pi*pi
            A(i,4) = qi*ri*ri;  A(i,5) = ri*pi*pi;  A(i,6) = ri*qi*qi
            A(i,7) = pi*qi*ri

      end select

      do l = 1, nf ! Set up RHS element for each of the functions

         b(i,l) = fnode(l,in) -             & ! Function data value at this node
           (pi*fnode(l,n1) + qi*fnode(l,n2) + ri*fnode(l,n3)) ! Linear interp. f
      end do

   end do ! Next node point

!  Factorize the m x n matrix:

   call hdecom (m, m, n, A, u)   ! Q R factors 

!  For each function, set up the right-hand-side, solve, and interpolate:

   do l = 1, nf

      call hsolve (m, m, n, A, u, b(1,l), residual)

      flinear = p * fnode(l,n1) + q * fnode(l,n2) + r * fnode(l,n3)

      select case (n)

         case (3) ! Quadratic

            finterp(l) = flinear +                                            &
               b(1,l) * p * q + b(2,l) * q * r + b(3,l) * r * p

         case (7) ! Cubic

            finterp(l) = flinear +                                            &
               b(1,l) * p * q * q + b(2,l) * p * r * r + b(3,l) * q * p * p + &
               b(4,l) * q * r * r + b(5,l) * r * p * p + b(6,l) * r * q * q + &
               b(7,l) * p * q * r

      end select

   end do ! Next function

   deallocate (A, b, u)

   end subroutine nonlinear_tri_interp
