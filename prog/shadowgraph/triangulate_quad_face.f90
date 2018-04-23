!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine triangulate_quad_face (xyz, iv)
!
!  Set up two sets of triangle vertices from a quad face by using the
!  shorter diagonal.  This is common to all 6 faces of a hex cell.
!
!  Jan. 2010  D.A.Saunders  Initial implementation for program SHADOWGRAPH.
!
!  Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   real,    intent (in)  :: xyz(3,4)  ! (x,y,z) for 4 corners of a quad face
   integer, intent (out) :: iv(3,2)   ! (i,j,k) for vertices of 2 triangles

!  Local variables:

   real :: diag1, diag2

!  Execution:

   diag1 = (xyz(1,2) - xyz(1,3))**2 + &
           (xyz(2,2) - xyz(2,3))**2 + &
           (xyz(3,2) - xyz(3,3))**2 
   diag2 = (xyz(1,1) - xyz(1,4))**2 + &
           (xyz(2,1) - xyz(2,4))**2 + &
           (xyz(3,1) - xyz(3,4))**2

   if (diag1 <= diag2) then
      iv(1,1) = 1;  iv(2,1) = 2;  iv(3,1) = 3
      iv(1,2) = 2;  iv(2,2) = 4;  iv(3,2) = 3
   else
      iv(1,1) = 1;  iv(2,1) = 2;  iv(3,1) = 4
      iv(1,2) = 1;  iv(2,2) = 4;  iv(3,2) = 3
   end if

   end subroutine triangulate_quad_face
