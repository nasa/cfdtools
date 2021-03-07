!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   module grid_block_structure

!  Module for a grid block with associated flow field variables.

   type grid_type
      real,    dimension(:,:,:),   pointer :: x,y,z  ! grid coordinates
      real,    dimension(:,:,:,:), pointer :: q      ! dependent variables
      integer, dimension(:,:,:),   pointer :: iblank ! 0 => suppress the point
      real    :: xmin, xmax      ! Data range for the grid block in x direction
      real    :: ymin, ymax      ! ... and y direction ...
      real    :: zmin, zmax      ! ... and z direction
      integer :: ni, nj, nk      ! Numbers of nodes in the i, j, k directions
      integer :: mi, mj, mk      ! Numbers of dependent variables
   end type grid_type

   end module grid_block_structure

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
