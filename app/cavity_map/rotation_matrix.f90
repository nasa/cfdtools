!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine rotation_matrix (nblocks, grid, nlist, list, cm, R, lambda)

!  This routine is prompted by the need to transform the grid blocks defining an
!  arbitrary cavity or gouge in a surface to a normalized cavity for comparisons
!  with other cases.  For the gouge-related subset of grid blocks, the principal
!  axes needed to align the cavity with the working axes are taken to be defined
!  by the eigenvectors of the inertia matrix established here.  The eigenvectors
!  in turn are stacked as three unit column vectors of a rotation matrix R which
!  describes the principal axes relative to the working axes.
!
!  References:  http://kwon3d.com/ (Theories and Practices of Motion Analysis);
!               Introduction to Robotics Mechanics & Control, John J. Craig
!
!  06/02/05  D. A. Saunders  Initial implementation.
!  06/18/05     "     "      Inertias need to be about the center of mass.
!  06/20/05     "     "      Replaced crude use of unit masses at the vertices
!                            with cell centroids/surface areas/unit densities.
!
!  David Saunders, ELORET/NASA Ames Research Center, Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use grid_block_structure    ! Module defining grid_type (Tecplot_io.f90)
   use surface_patch_utilities ! Module containing centroids_areas, etc.

   implicit none

!  Arguments:

   integer, intent (in)          :: nblocks        ! Dimension of grid(:)
   type (grid_type), intent (in) :: grid(nblocks)  ! Grid blocks
   integer, intent (in)          :: nlist          ! Active length of list(:)
   integer, intent (in)          :: list(nlist)    ! Blocks to be processed
   real,    intent (out)         :: cm(3)          ! Center of mass of cavity
   real,    intent (out)         :: R(3,3)         ! Rotation matrix describing
                                                   ! the principal axes of the
                                                   ! cavity w.r.t. Ox, Oy, Oz;
   real,    intent (out)         :: lambda(3)      ! Eigenvalues, ascending

!  Consider application to a flat plate longer than it is wide, aligned with the
!  axes.  The smallest moment of inertia (<-> smallest eigenvalue) will be about
!  the length axis, the second smallest about the width axes, and the largest is
!  about the axis perpendicular to the plate.  This is the assumed order for the
!  cavity patches.  If the cavity is deeper than it is wide, the application may
!  have to swap columns 2 and 3 of the rotation matrix.

!  Local constants:

   real, parameter :: zero = 0.

!  Local variables:

   integer :: i, j, l, ib, ier, ni, nj
   real    :: area, cell_area
   real    :: x, y, z, xx, yy, zz, xy, yz, zx
   real    :: fv1(3), fv2(3)                      ! Work-space for eigen routine
   real    :: A(3,3)                              ! Inertia matrix; columns
                                                  ! correspond to eigenvectors
   type (grid_type), allocatable :: grid_cells(:) ! Ancillary fields for cell
                                                  ! centroids and surface areas
!  Execution:

! Adding centroids & areas as further fields of grid(*) is not an option if
! Tecplot_io is in use, so use a parallel data structure of the same type.

   allocate (grid_cells(nblocks))

!  Calculate the center of mass for the relevant surface cells (unit density):

   area = 0;  x = zero;  y = zero;  z = zero

   do l = 1, nlist
      ib = list(l)

      call centroids_areas (grid(ib), grid_cells(ib))

      do j = 1, grid_cells(ib)%nj
         do i = 1, grid_cells(ib)%ni
            cell_area = grid_cells(ib)%q(1,i,j,1)
            x    = cell_area * grid_cells(ib)%x(i,j,1) + x
            y    = cell_area * grid_cells(ib)%y(i,j,1) + y
            z    = cell_area * grid_cells(ib)%z(i,j,1) + z
            area = cell_area + area
         end do
      end do
   end do

   cm(1) = x / area;  cm(2) = y / area;  cm(3) = z / area

!  Generate the inertia matrix for rotations about the center of mass:

   xx = zero;  yy = zero;  zz = zero
   xy = zero;  yz = zero;  zx = zero

   do l = 1, nlist
      ib = list(l)
      do j = 1, grid_cells(ib)%nj
         do i = 1, grid_cells(ib)%ni
            x = grid_cells(ib)%x(i,j,1) - cm(1)
            y = grid_cells(ib)%y(i,j,1) - cm(2)
            z = grid_cells(ib)%z(i,j,1) - cm(3)
            area = grid_cells(ib)%q(1,i,j,1)
            xx = area * (y*y + z*z) + xx
            yy = area * (z*z + x*x) + yy
            zz = area * (x*x + y*y) + zz
            xy = area *  x*y + xy
            yz = area *  y*z + yz
            zx = area *  z*x + zx
         end do
      end do
      deallocate (grid_cells(ib)%x, grid_cells(ib)%y, grid_cells(ib)%z)
      deallocate (grid_cells(ib)%q)
   end do

   deallocate (grid_cells)

   A(1,1) =  xx;  A(1,2) = -xy;  A(1,3)  = -zx
   A(2,1) = -xy;  A(2,2) =  yy;  A(2,3)  = -yz
   A(3,1) = -zx;  A(3,2) = -yz;  A(3,3)  =  zz

   write (6, '(/, a, /, (1p, 3e15.6, 3x))')           &
      ' Inertia matrix:', (A(i,:), i = 1, 3)

!  Calculate the eigenvectors (and values) of the inertia matrix:

   call rs (3, 3, A, lambda, 1, R, fv1, fv2, ier)

   if (ier /= 0) then
      write (*, *) ' Eigenvector routine "rs" reports "tql2" error: ', ier
      stop
   end if

   end subroutine rotation_matrix
