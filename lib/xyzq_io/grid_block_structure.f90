!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   module grid_block_structure

!  Module for a structured grid block with associated flow field variables
!  and optional PLOT3D-type blanking of grid points.
!  The %mi/mj/mk variables were inherited from NASA LaRC, but all applica-
!  tions at NASA ARC have %mi = %ni, etc.

   implicit none

   type grid_type
      real,    dimension(:,:,:),   pointer :: x,y,z  ! Grid coordinates
      real,    dimension(:,:,:,:), pointer :: q      ! Fns. in (1:nf,:,:,:)
      integer, dimension(:,:,:),   pointer :: iblank ! 0 => suppress the pt

      real    :: xmin, xmax        ! Data range for the block in the x dir.
      real    :: ymin, ymax        ! ... and the y direction ...
      real    :: zmin, zmax        ! ... and the z direction

      integer :: ni, nj, nk        ! # grid vertices in i, j, k directions
      integer :: mi, mj, mk        ! # dependent variables in each dir.
   end type grid_type

   end module grid_block_structure

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
