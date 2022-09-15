!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine cross_product (a, b, c)

!  Calculate the cross product c = a x b for vectors in three-space.
!  The argument description is obvious.
!
!  03/15/22  D.A.Saunders Adapted from a trianulation_io.f90 internal procedure.
!
!  Author:  David Saunders, AMA, Inc. at NASA Ames Research Center, CA
!
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Arguments:

   real, intent (in)  :: a(3), b(3)
   real, intent (out) :: c(3)

!  Execution:

   c(1) = a(2)*b(3) - a(3)*b(2)
   c(2) = a(3)*b(1) - a(1)*b(3)
   c(3) = a(1)*b(2) - a(2)*b(1)

   end subroutine cross_product
