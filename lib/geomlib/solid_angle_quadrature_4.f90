!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine solid_angle_quadrature_4 (ni, nj, nf, x, y, z, f, px, py, pz, &
                                        fquad)
!  Description:
!
!     Integrate one or more functions defined at the vertices of a structured
!     surface with respect to the solid angle subtended at a given point Pxyz.
!     The cell-centered option found in the earlier triangulated zone utility
!     appears to be redundant here.  One can always convert the grid vertices
!     to a grid of cell centroids.
!
!     A single zone/block/patch is treated here.
!
!
!  History:
!
!     Feb. 10, 2020  D.A.Saunders  Adaptation of solid_angle_quadrature_3, as
!                                  needed for data on a rectangular array of
!                                  "pixels" used to predict radiation intensity
!                                  observable during a hypersonic atmospheric
!                                  entry such as that of Hayabusa 2.
!
!  Author:  David Saunders, AMA, Inc./NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)  :: ni, nj         ! (Packed) patch array dimensions
   integer, intent (in)  :: nf             ! # functions to be integrated >= 1
   real,    intent (in), &
       dimension (ni,nj) :: x, y, z        ! Patch coordinates
   real,    intent (in)  :: f(nf,ni,nj)    ! Discretized function values to be
                                           ! integrated
   real,    intent (in)  :: px, py, pz     ! Coordinates of P
   real,    intent (out) :: fquad(nf)      ! Desired integral estimate(s)

!  Local constants:

   real, parameter :: fourth = 0.25, zero = 0.0

!  Local variables;

   integer :: i, j
   real    :: darea, domega

!  Procedures:

   external :: solid_angle_quad  ! Solid angle element at P for a given cell

!  Execution:

   fquad(:) = zero

   do j = 1, nj - 1
      do i = 1, ni - 1
         call solid_angle_quad (ni, nj, x, y, z, i, j, px, py, pz, domega)
         fquad(:) = fquad(:) + domega * &
                    (f(:,i,j) + f(:,i+1,j) + f(:,i,j+1) + f(:,i+1,j+1))
      end do
   end do

   fquad(:) = fourth * fquad(:)

   end subroutine solid_angle_quadrature_4
