!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine trapezoidal_quadrature_2d (ni, nj, nf, x, y, f, area,
     >                                      volume)
!
!     Integrate one or more functions over the area of a 2-D structured grid by
!     the trapezoidal method.  The quadrilateral surface cells are NOT assumed
!     to be rectangles.  Each subvolume is approximated as the average of four
!     function values multiplied by the surface cell area.
!
!     All functions are integrated and arrays are assumed to be packed.
!
!     03/22/2007  D. Saunders  Initial implementation.
!     03/23/2007   "     "     Make total area another output.
!
!     Author:  David Saunders, ELORET/NASA Ames Research Center, Moffett Fld. CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

!     Arguments:

      integer, intent (in)  :: ni, nj      ! Dimensions of (packed) 2-space grid
      integer, intent (in)  :: nf          ! # functions at each grid point >= 1
      real,    intent (in)  :: x(ni,nj)    ! Grid coordinates
      real,    intent (in)  :: y(ni,nj)
      real,    intent (in)  :: f(nf,ni,nj) ! Function(s) at each grid point
      real,    intent (out) :: area        ! Estimated grid area
      real,    intent (out) :: volume(nf)  ! Integration(s) of f over (x,y).

!     Local constants:

      real, parameter :: fourth = 0.25, zero = 0.

!     Local variables:

      integer :: i, j, n
      real    :: cell_area, xp(5), yp(5)  ! Allow for closing list of vertices

!     Execution:

      area      = zero
      volume(:) = zero

      do j = 1, nj - 1

         do i = 1, ni - 1

            xp(1) = x(i,j)
            yp(1) = y(i,j)
            xp(2) = x(i+1,j)
            yp(2) = y(i+1,j)
            xp(3) = x(i+1,j+1)
            yp(3) = y(i+1,j+1)
            xp(4) = x(i,j+1)
            yp(4) = y(i,j+1)

            call areaxy (4, xp, yp, cell_area)  ! Polygon utility

            area      = area      + cell_area
            volume(:) = volume(:) + cell_area *
     >         (f(:,i,j) + f(:,i+1,j) + f(:,i,j+1) + f(:,i+1,j+1))

         end do

      end do

      volume(:) = fourth * volume(:)

      end subroutine trapezoidal_quadrature_2d
