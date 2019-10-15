!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine trapezoidal_quadrature_3d (ni, nj, nf, x, y, z, f, area, volume)
!
!  Description:
!
!     Integrate one or more functions over the area of a 3-space structured
!     surface grid (or grid block) by the trapezoidal method.  Each subvolume
!     is approximated as the average of four function values multiplied by the
!     surface cell area.
!
!     All functions are integrated and arrays are assumed to be packed.
!
!  History:
!
!     03/22/2007  D.A.Saunders  Initial implementation of 2-space form.
!     05/02/2017    "     "     3-space surface adaptation.
!
!  Author:  David Saunders, AMA, Inc. at NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)  :: ni, nj      ! Dimensions of (packed) 3-space grid
   integer, intent (in)  :: nf          ! # functions at each grid point >= 1
   real,    intent (in)  :: x(ni,nj)    ! Grid coordinates
   real,    intent (in)  :: y(ni,nj)
   real,    intent (in)  :: z(ni,nj)
   real,    intent (in)  :: f(nf,ni,nj) ! Function(s) at each grid point
   real,    intent (out) :: area        ! Estimated surface grid area
   real,    intent (out) :: volume(nf)  ! Integration(s) of f over (x,y).

!  Local constants:

   real, parameter :: fourth = 0.25, zero = 0.

!  Local variables:

   integer :: i, j
   real    :: cell_area

!  Procedures:

   real, external :: area4  ! Quadrilateral area utility

!  Execution:

   area      = zero
   volume(:) = zero

   do j = 1, nj - 1
      do i = 1, ni - 1
         cell_area = area4 (x(i,  j  ), y(i,  j  ), z(i,  j  ), &
                            x(i+1,j  ), y(i+1,j  ), z(i+1,j  ), &
                            x(i+1,j+1), y(i+1,j+1), z(i+1,j+1), &
                            x(i,  j+1), y(i,  j+1), z(i,  j+1))
         area      = area      + cell_area
         volume(:) = volume(:) + cell_area * fourth * &
                     (f(:,i,j) + f(:,i+1,j) + f(:,i,j+1) + f(:,i+1,j+1))
      end do
   end do

   end subroutine trapezoidal_quadrature_3d
