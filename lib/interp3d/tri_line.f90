!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine tri_line (xyz1, xyz2, xyz3, xt, yt, zt, p, q, r, istate)

!  For the 3-space plane defined by the indicated vertices of a triangle, find
!  the intersection with a line parallel to a coordinate axis.  For consistency
!  with a common CFD convention, the line is taken to be normal to the y = 0
!  plane of symmetry, i.e., the line is x = xt and z = zt, and we need to find
!  yt in the plane of the triangle.  For other lines, rotating the triangle
!  coordinates allows use of this choice.
!
!  If the point of intersection lies inside the triangle, all of the associated
!  linear interpolation coefficients p, q, r will lie in [0., 1.] and return
!  flag istate will be set to 1.  If the point is outside the triangle, it will
!  be set to 0 and the interpolation coefficients may still be of use.  If the
!  triangle is effectively perpendicular (edge on) to the y = 0 plane, the
!  (x, z) coordinates of the vertices will appear to be collinear, which will
!  be detected as a singular matrix in the coding.  In this case, if (xt, zt) is
!  also collinear with the triangle edge (and not off its ends), istate will be
!  returned as 2 since there are effectively two points of intersection with the
!  3-space triangle.  Otherwise, there are no intersection points and istate
!  will be returned as -1.  It is up to the application to handle the possible
!  cases appropriately.
!
!  The strategy is to solve the (x, z)-space problem for the triangle projected
!  onto the y = 0 plane (giving linear interpolation coefficients p, q, r) then
!  to use those coefficients to find yt in the 3-space plane of the triangle:
!
!                            x1 p + x2 q + x3 r = xt
!                            z1 p + z2 q + z3 r = zt
!                               p +    q +    r = 1
!
!  This 3x3 is reducible to a 2x2 and the solution uniquely defines p, q and r.
!  Then:
!                            yt = p y1 + q y2 + r y3
!
!  Linear system solver LUSOLVE could be used for the 2x2, but it is in-lined
!  here for efficiency (as done earlier in subroutine BILINT).
!
!  01/08/2011  D.A.Saunders  Variant of nearest-point least squares techniques
!                            for intersecting a line with a triangle as needed
!                            to construct shadowgraphs from CFD solutions.
!                            Working with a general 3-space line would be
!                            unnecessarily inefficient for the common case of
!                            lines perpendicular to the y = 0 symmetry plane.
!
!  Author:  David Saunders, ELORET/NASA Ames Research Center, Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   real, intent (in)  :: xyz1(3)   ! (x,y,z) coords. of vertex 1
   real, intent (in)  :: xyz2(3)   ! (x,y,z) coords. of vertex 2
   real, intent (in)  :: xyz3(3)   ! (x,y,z) coords. of vertex 3
   real, intent (in)  :: xt, zt    ! Coords. defining line _|_ y = 0
   real, intent (out) :: yt        ! Coord. of line/triangle intersection
                                   ! if it exists (istate = 1 case)
   real, intent (out) :: p, q, r   ! Linear interpolation coefficients
                                   ! associated with the intersection
                                   ! if it exists (istate = 0, 1 cases)
   integer, intent (out) :: istate ! Flag distinguishing the cases:
                                   !  2 => line is in the triangle plane;
                                   !  1 => intersection is inside  the triangle;
                                   !  0 => intersection is outside the triangle;
                                   ! -1 => line is parallel to triangle, offset
                                         

!  Local constants:

   real, parameter :: eps = 1.e-14, one = 1.0, zero = 0.0

!  Local variables:

   integer :: m
   real    :: A(2,2), b(2), Anorm, t, tol

!  Execution:

!  Set up the linear system that produces linear interpolation coefficients
!  for a point in the plane defined by the triangle's (x, z) coordinates:

   A(1,1) = xyz1(1) - xyz3(1);  A(1,2) = xyz2(1) - xyz3(1);  b(1) = xt - xyz3(1)
   A(2,1) = xyz1(3) - xyz3(3);  A(2,2) = xyz2(3) - xyz3(3);  b(2) = zt - xyz3(3)

   Anorm = max (abs (A(1,1)), abs (A(2,1)), abs (A(1,2)), abs (A(2,2)))
   tol   = eps * Anorm

!! call lusolve (2, 2, A, b, istate)    ! Solve A p = b

!  Solve the 2x2 system explicitly to avoid an external reference and to avoid
!  degenerate loops:

!  Perform the LU factorization, solving L y = b as we go:

!  Determine pivot element:

   m = 1
   t = A(1,1)
   if (abs (A(2,1)) > abs (t)) then
      A(1,1) = A(2,1)
      A(2,1) = t
      m = 2
   end if

   if (abs (A(1,1)) < tol) go to 90  ! Singular matrix

!  Store the multiplier:

   A(2,1) = -A(2,1) / A(1,1)

!  Apply the multiplier to current submatrix and RHS:

   t      = A(m,2)
   A(m,2) = A(1,2)
   A(1,2) = t
   A(2,2) = A(2,2) + A(2,1) * t

   t    = b(m)
   b(m) = b(1)
   b(1) = t
   b(2) = b(2) + A(2,1) * t

   if (abs (A(2,2)) < tol) go to 90  ! Singular matrix

!  Back substitution (solution of U p = y):

   b(2) =  b(2) / A(2,2)
   b(1) = (b(1) - A(1,2) * b(2)) / A(1,1)

   p = b(1)     ! Desired solution
   q = b(2)
   r = one - b(1) - b(2)

!  Third coordinate of the intersection of the line with the triangle plane:

   yt = p * xyz1(2) + q * xyz2(2) + r * xyz3(2)

   istate = 0  ! If the intersection is outside the triangle

   if (p >= zero) then
      if (p <= one) then
         if (q >= zero) then
            if (q <= one) then
               if (r >= zero) then
                  if (r <= one) istate = 1  ! Intersection not outside triangle
               end if
            end if
         end if
      end if
   end if

   go to 99  ! Single return philosophy

90 continue  ! Matrix was singular: the three (x, z) vertices are collinear

!  Is the line (x, z) also collinear?

   if (abs ((xyz1(1) - xt)*(xyz2(3) - zt) - (xyz2(1) - xt)*(xyz1(3) - zt)) <   &
       tol) then
      istate = 2  ! Line is in the plane of the triangle
   else
      istate = -1 ! Line is parallel to the triangle plane but offset
   end if

99 return

   end subroutine tri_line
