!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program thin_grid_2d
!
!  Description:
!
!     THIN_GRID_2D extracts every nth and mth point from all blocks of a 2D
!  multiblock grid in PLOT2D form.
!
!  History:
!
!     08/28/12  D.A.Saunders  Adaptation of THIN_FLOW (easier than THIN_GRID
!                             which has the option to do different things to
!                             different blocks).
!
!  Author:  David Saunders, ERC, Inc. at NASA Ames Research Center, CA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_block_structure  ! For one 3D grid block; OK for 2D
   use xyq_io_module         ! 2D I/O utilities analogous to the 3D utilities

   implicit none

!  Constants:

   integer, parameter :: &
      lunin  = 1,        &
      lunout = 2,        &
      lunkbd = 5,        &
      luncrt = 6

   logical, parameter :: &
      true = .true.

!  Variables:

   integer :: &
      i, ib, ii, inc, ios, j, jj, jnc, nblocks, ni, nj, npts

   logical :: &
      formatted_in, formatted_out

!  Composite data types:

   type (grid_type), pointer, dimension (:) :: &
      gridin, gridout

!  Execution:
!  ----------

   call file_prompt_2d (lunin, 'input 2D grid', 'old', true, formatted_in, ios)
   if (ios /= 0) go to 99

!  Read the number of blocks, allocate them, and read dimensions:

   call xy_header_io (1, lunin, formatted_in, nblocks, gridin, ios)
   if (ios /= 0) go to 99

   write (luncrt, '(a)', advance='no') ' Increments to use for thinning i, j: '
   read  (lunkbd, *) inc, jnc

   call file_prompt_2d (lunout, 'output grid', 'unknown', true, formatted_out, &
                        ios)
   if (ios /= 0) go to 99

   allocate (gridout(nblocks))

   do ib = 1, nblocks
      gridout(ib)%ni = (gridin(ib)%ni + inc - 1) / inc
      gridout(ib)%nj = (gridin(ib)%nj + jnc - 1) / jnc
   end do

   call xy_header_io (2, lunout, formatted_out, nblocks, gridout, ios)
   if (ios /= 0) go to 99

!  Thin the blocks one at a time:

   do ib = 1, nblocks

      call xy_allocate (gridin(ib), ios)
      if (ios /= 0) go to 99

      ni = gridin(ib)%ni;  nj = gridin(ib)%nj;  npts = ni*nj

      call xy_block_io (1, lunin, formatted_in, npts, &
                        gridin(ib)%x, gridin(ib)%y, ios)
      if (ios /= 0) go to 99

      call xy_allocate (gridout(ib), ios)
      if (ios /= 0) go to 99

      jj = 0
      do j = 1, nj, jnc
         jj = jj + 1
         ii = 0
         do i = 1, ni, inc
            ii = ii + 1
            gridout(ib)%x(ii,jj,1) = gridin(ib)%x(i,j,1)
            gridout(ib)%y(ii,jj,1) = gridin(ib)%y(i,j,1)
         end do
      end do

      npts = gridout(ib)%ni * gridout(ib)%nj

      call xy_block_io (2, lunout, formatted_out, npts, &
                        gridout(ib)%x, gridout(ib)%y, ios)
      if (ios /= 0) go to 99
         
      deallocate (gridin(ib)%x, gridin(ib)%y, gridout(ib)%x, gridout(ib)%y)

   end do

99 continue

!  stop ! Avoid system dependencies.

   end program thin_grid_2d
