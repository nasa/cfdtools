!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine solid_angle_quadrature_3 (nnodes, ntris, nf, xyz, f, conn, &
                                        cell_centered, Pxyz, fquad)
!  Description:
!
!     Integrate one or more functions defined at either the vertices or the
!     cell centers of a triangulated surface with respect to the solid angle
!     subtended at a given point Pxyz.  If the function values are at the
!     vertices, the average of three is used with each elemental solid angle.
!
!     This was prompted by a need to integrate radiance vs. solid angle in
!     order to estimate radiative heat flux at a point more correctly than
!     with the tangent-slab approximation.
!
!  History:
!
!     Mar. 26, 2014  D.A.Saunders  Companion to the earlier solid_angle_tri.
!
!  Author:  David Saunders, ERC, Inc./NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)  :: nnodes         ! # x/y/z points/nodes/vertices
   integer, intent (in)  :: ntris          ! # triangular elements
   integer, intent (in)  :: nf             ! # functions to be integrated >= 1
   real,    intent (in)  :: xyz(3,nnodes)  ! Node coordinates
   real,    intent (in)  :: f(nf,*)        ! Discretized function values to be
                                           ! integrated; the second dimension
                                           ! is nnodes|ntris if cell_centered
                                           ! is F|T respectively
   integer, intent (in)  :: conn(3,ntris)  ! Vertex indices in xyz(:,*)
   logical, intent (in)  :: cell_centered  ! T if f(:,n) applies to element n;
                                           ! F if f(:,n)   "     "  vertex  n
   real,    intent (in)  :: Pxyz(3)        ! Coordinates of P
   real,    intent (out) :: fquad(nf)      ! Desired integral estimate(s)

!  Local constants:

   real, parameter :: third = 1.0/3.0, zero = 0.0

!  Local variables;

   integer :: itri
   real    :: darea, domega

!  Procedures:

   external :: tri_normal_and_area

!  Execution:

   fquad(:) = zero

   if (cell_centered) then

      do itri = 1, ntris
         call solid_angle_tri (nnodes, ntris, xyz, conn, itri, Pxyz, darea, &
                               domega)
         fquad(:) = fquad(:) + f(:,itri)*domega
      end do

   else  ! Average the three vertex-centered function values for each element

      do itri = 1, ntris
         call solid_angle_tri (nnodes, ntris, xyz, conn, itri, Pxyz, darea, &
                               domega)
         fquad(:) = fquad(:) + domega * &
                    (f(:,conn(1,itri)) + f(:,conn(2,itri)) + f(:,conn(3,itri)))
      end do

      fquad(:) = third * fquad(:)

   end if

   end subroutine solid_angle_quadrature_3
