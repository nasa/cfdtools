!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine vertices_to_centers_2d (initialize, mode, nf, block_in, block_out)

!  Description:
!
!     This is a 2-space variant of VERTICES_TO_CENTERS, from which the following
!     description is adapted.
!
!     Transfer a grid and optional flow variables from cell vertices to cell
!     centers for one grid block.  The block may be a volume or surface grid
!     or just a single grid line.
!
!     DPLR-type treatment of "halo" or "ghost" cells is supported by mode 1
!     (such cells are on the block boundaries), but mode 2 also provides for
!     omitting halo cells.
!
!     An original intent of including centers-to-vertices options has been
!     abandoned.  It makes little sense for a cell-centered grid with no halo
!     cells, because true boundary vertices have been lost.  Even for function
!     data, the usual implementation of averaging neighbors as for v-to-c is
!     only approximate in the interior, and should involve extrapolation at
!     the boundaries but typically doesn't.
!
!     Note that any input dimension of 1 means NO halo cells in that index
!     direction, regardless.  For instance, a j surface input with mode = 1
!     still has block_out%nj = 1.
!
!     The input block is NOT modified in place.  Likely usage will loop over all
!     blocks with "initialize" input as T, so that an output file header can be
!     written before each block is processed in a second loop.  This requires
!     the array of modified blocks to be separate from the initial blocks.
!
!     Deallocates of input blocks should therefore be performed at a higher
!     level.  Output arrays should also be allocated before calling this
!     routine with initialize = F.
!
!     For now, processing the grid is not optional.  If processing function
!     data only is ever needed, nf < 0 could be the appropriate flag.
!
!  History:
!
!     01/06/10  D.A.Saunders  Date of 3-space routine from which this 2-space
!                             variant is adapted.
!
!  Author:  David Saunders, AMA, Inc. at NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_block_structure                  ! Derived data type for a grid blk.

   implicit none

!  Arguments:

   logical, intent (in) :: initialize        ! T  = set output dimensions only;
                                             ! F  = set output x,y[,q] (only)
   integer, intent (in) :: mode              ! 1  = vertices -> centers + halos;
                                             ! 2  = vertices -> centers/no halos
   integer, intent (in) :: nf                ! # flow variables >= 0
   type (grid_type), intent (in) :: block_in ! Input grid block
   type (grid_type), intent (inout) :: block_out
                                             ! If initialize = T, the revised
                                             ! block dimensions are set, else
                                             ! they are assumed to be set and
                                             ! just x,y[,q] are processed
!  Local constants:

   real, parameter :: fourth = 0.25

!  Local variables:

   integer :: i, il, ir, j, jl, jr, m, mi, mj, ni, nj

!  Execution:

   select case (mode)

      case (1) ! Vertices -> centers with halos

         ni = block_in%ni;  mi = ni + 1;  if (ni == 1) mi = 1 
         nj = block_in%nj;  mj = nj + 1;  if (nj == 1) mj = 1

         if (initialize) go to 90  ! Avoid excessive indenting

         do j = 1, mj
            jl = max (1, j - 1);  jr = min (nj, j)
            do i = 1, mi
               il = max (1, i - 1);  ir = min (ni, i)
               block_out%x(i,j,1) = fourth * (                &
                  block_in%x(il,jl,1) + block_in%x(ir,jl,1) + &
                  block_in%x(il,jr,1) + block_in%x(ir,jr,1) )
               block_out%y(i,j,1) = fourth * (                &
                  block_in%y(il,jl,1) + block_in%y(ir,jl,1) + &
                  block_in%y(il,jr,1) + block_in%y(ir,jr,1) )
            end do
         end do

         if (nf > 0) then

            do j = 1, mj
               jl = max (1, j - 1);  jr = min (nj, j)
               do i = 1, mi
                  il = max (1, i - 1);  ir = min (ni, i)
                  block_out%q(:,i,j,1) = fourth * (                  &
                     block_in%q(:,il,jl,1) + block_in%q(:,ir,jl,1) + &
                     block_in%q(:,il,jr,1) + block_in%q(:,ir,jr,1) )
               end do
            end do

         end if

      case (2) ! Vertices -> centers without halos

         ni = block_in%ni;  mi = ni - 1;  if (ni == 1) mi = 1
         nj = block_in%nj;  mj = nj - 1;  if (nj == 1) mj = 1
 
         if (initialize) go to 90  ! Avoid excessive indenting

         do j = 1, mj
            jl = j;  jr = j + 1;  if (mj == 1) jr = 1
            do i = 1, mi
               il = i;  ir = i + 1;  if (mi == 1) ir = 1
               block_out%x(i,j,1) = fourth * (                &
                  block_in%x(il,jl,1) + block_in%x(ir,jl,1) + &
                  block_in%x(il,jr,1) + block_in%x(ir,jr,1) )
               block_out%y(i,j,1) = fourth * (                &
                  block_in%y(il,jl,1) + block_in%y(ir,jl,1) + &
                  block_in%y(il,jr,1) + block_in%y(ir,jr,1) )
            end do
         end do

         if (nf > 0) then

            do j = 1, mj
               jl = j;  jr = j + 1;  if (mj == 1) jr = 1
               do i = 1, mi
                  il = i;  ir = i + 1;  if (mi == 1) ir = 1
                  block_out%q(:,i,j,1) = fourth * (                  &
                     block_in%q(:,il,jl,1) + block_in%q(:,ir,jl,1) + &
                     block_in%q(:,il,jr,1) + block_in%q(:,ir,jr,1) )
               end do
            end do

         end if

      case default

!        Deal with bad inputs at the higher level.

   end select

90 continue

   if (initialize) then
      block_out%ni = mi;  block_out%mi = mi
      block_out%nj = mj;  block_out%mj = mj
   end if

   end subroutine vertices_to_centers_2d
