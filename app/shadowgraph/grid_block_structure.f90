!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   module grid_block_structure  ! SHADOWGRAPH version

!  Module for a structured grid block with associated flow field variables.
!  This version is consistent with the xyzq_io package for reading flow species
!  densities as %q, but adds further arrays for derived functions as needed by
!  program SHADOWGRAPH.  Here, %mi = %ni, etc.

   type grid_type
      real, dimension(:,:,:),   pointer :: x,y,z     ! Grid coordinates
      real, dimension(:,:,:,:), pointer :: q         ! Input species densities
      real, dimension(:,:,:),   pointer :: r         ! Refractive index n (- 1)
      real, dimension(:,:,:),   pointer :: fx        ! Partial dn/dx
      real, dimension(:,:,:),   pointer :: fy        ! Partial dn/dy
      real, dimension(:,:,:),   pointer :: fz        ! Partial dn/dz

      integer, dimension(:,:,:),   pointer :: iblank ! 0 => suppress the point

      integer :: ni, nj, nk     ! # grid vertices in i, j, k directions
      integer :: mi, mj, mk     ! # dependent variables in each direction
      real    :: xmin, xmax, ymin, ymax, zmin, zmax  ! Grid block data range
   end type grid_type

   end module grid_block_structure

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
