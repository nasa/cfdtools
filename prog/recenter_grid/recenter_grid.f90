!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      program recenter_grid
!
!     Description:
!
!        RECENTER_GRID converts all blocks of a multiblock grid to cell-centered
!     equivalents, including addition of single-layer halo cells whose faces lie
!     in the original block faces.
!
!     Procedures:
!
!        XYZQ_IO package  I/O utilities for PLOT3D grid and function files
!
!     History:
!
!        02/10/05  D. Saunders  Initial adaptation of THIN_GRID to enable
!                               inspection of interpolated cell-centered flow
!                               solutions output by LaRC's MORPH program.
!
!     Author:  David Saunders, Eloret/NASA Ames, Moffett Field, CA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Modules:

      use grid_block_structure
      use xyzq_io_module

      implicit none

!     Constants:

      integer, parameter :: &
         lunin  = 1,        &
         lunout = 2,        &
         lunkbd = 5,        &
         luncrt = 6

      logical, parameter :: &
         true = .true.

!     Variables:

      integer :: &
         i, ib, ii, inc, ios, j, jj, jnc, k, kk, knc, nblocks, ni, nj, nk, npts

      logical :: &
         formatted_in, formatted_out

!     Composite data types:

      type (grid_type), pointer, dimension (:) :: &
         gridin, gridout

!     Execution:
!     ----------

      call file_prompt (lunin,  -luncrt, 'input grid', 'old', true, &
                        formatted_in, ios)
      if (ios /= 0) go to 999

!     Read the header records and allocate work-space:

      call xyz_header_io (1, lunin, formatted_in, nblocks, gridin, ios)

      if (ios /= 0) go to 999

      call file_prompt (lunout, -luncrt, 'output grid', 'unknown', true, &
                        formatted_out, ios)
      if (ios /= 0) go to 999

      allocate (gridout(nblocks))

      do ib = 1, nblocks
         gridout(ib)%ni = gridin(ib)%ni + 1
         gridout(ib)%nj = gridin(ib)%nj + 1
         gridout(ib)%nk = gridin(ib)%nk + 1
      end do

      call xyz_header_io (2, lunout, formatted_out, nblocks, gridout, ios)

      if (ios /= 0) go to 999

!     Process the blocks one at a time:

      do ib = 1, nblocks

         call xyz_allocate (gridin(ib), ios)
         if (ios /= 0) go to 999

         ni = gridin(ib)%ni;  nj = gridin(ib)%nj;  nk = gridin(ib)%nk
         npts = ni * nj * nk

         call xyz_block_io (1, lunin, formatted_in, npts, &
                            gridin(ib)%x, gridin(ib)%y, gridin(ib)%z, ios)
         if (ios /= 0) go to 999

         call xyz_allocate (gridout(ib), ios)
         if (ios /= 0) go to 999

         call recenter_volume (gridin(ib), gridout(ib))

         npts = gridout(ib)%ni * gridout(ib)%nj * gridout(ib)%nk

         call xyz_block_io (2, lunout, formatted_out, npts, &
                            gridout(ib)%x, gridout(ib)%y, gridout(ib)%z, ios)
         if (ios /= 0) go to 999
         
         deallocate (gridin(ib)%x,  gridin(ib)%y,  gridin(ib)%z,  stat=ios)
         deallocate (gridout(ib)%x, gridout(ib)%y, gridout(ib)%z, stat=ios)

      end do

  999 continue

! *** stop ! Avoid system dependencies.

      end program recenter_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine recenter_volume (block_in, block_out)

!  For one volume grid block, convert vertex-centered coordinates to cell-
!  centered coordinates including one-layer halo cells.
!  Input and output block dimensions should be set upon entry.
!
!  12/03/04  Adaptation of version in RADIAL_INTERP, with no local I/O.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use grid_block_structure  ! External module

   implicit none

!  Arguments:

   type (grid_type), intent (in)    :: block_in  ! Vertex-cntrd. block, no halos
   type (grid_type), intent (inout) :: block_out ! Cell-centered block, w/ halos

!  Local constants:

   real, parameter :: eighth = 0.125

!  Local variables:

   integer :: i, j, k, il, ir, jl, jr, kl, kr, ni, nj, nk

!  Execution:

   ni = block_in%ni  ! Vertex-centered point counts (no halos)
   nj = block_in%nj
   nk = block_in%nk

   do k = 1, nk + 1
      kl = max (1, k - 1);  kr = min (nk, k)
      do j = 1, nj + 1
         jl = max (1, j - 1);  jr = min (nj, j)
         do i = 1, ni + 1
            il = max (1, i - 1);  ir = min (ni, i)
            block_out%x(i,j,k) = eighth * (                  &
               block_in%x(il,jl,kl) + block_in%x(ir,jl,kl) + &
               block_in%x(il,jr,kl) + block_in%x(ir,jr,kl) + &
               block_in%x(il,jl,kr) + block_in%x(ir,jl,kr) + &
               block_in%x(il,jr,kr) + block_in%x(ir,jr,kr) )
            block_out%y(i,j,k) = eighth * (                  &
               block_in%y(il,jl,kl) + block_in%y(ir,jl,kl) + &
               block_in%y(il,jr,kl) + block_in%y(ir,jr,kl) + &
               block_in%y(il,jl,kr) + block_in%y(ir,jl,kr) + &
               block_in%y(il,jr,kr) + block_in%y(ir,jr,kr) )
            block_out%z(i,j,k) = eighth * (                  &
               block_in%z(il,jl,kl) + block_in%z(ir,jl,kl) + &
               block_in%z(il,jr,kl) + block_in%z(ir,jr,kl) + &
               block_in%z(il,jl,kr) + block_in%z(ir,jl,kr) + &
               block_in%z(il,jr,kr) + block_in%z(ir,jr,kr) )
         end do
      end do
   end do

   end subroutine recenter_volume
