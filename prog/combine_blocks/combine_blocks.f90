!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program combine_blocks
!
!  Description:
!
!     COMBINE_BLOCKS appends the blocks of a second grid to those of an
!  initial grid.  Accompanying flow solutions are optional.  This could be
!  convenient for overlaying corresponding blocks from two grids via PLOT3D,
!  which handles just one file at a time.  More recently, it has proved handy
!  for appending tile cavity blocks to smooth OML grids during damage/repair
!  studies for the Shuttle Return-to-Flight program.
!
!     Some specialization is pragmatic for cavity cases when solution files are
!  present.  Wing leading edge plug repair cases are also handled by this same
!  option, except that a second specialized procedure (PLUG_INTERP) should then
!  be used to adjust the outer boundary flow interpolations.
!
!     The input and output files are prompted for.  If a flow solution is
!  present with the first grid but not the second, then blocks corresponding
!  to the second grid file (presumably in a cavity) are added to the output
!  solution file and their flow fields are initialized from the k = 1 layer
!  of the flow for the first set of blocks as follows:
!
!         pressure is reduced by a factor of 10 (default; now prompted for);
!         velocity components are zeroed;
!         remaining state variables are retained.
!
!     The nearest point in the first grid's k = 1 surface layer is located
!  efficiently (for each volume point of the cavity blocks), via an ADT search.
!  This avoids dealing with interface files, and produces smoothly varying flow
!  within the cavity blocks which is presumably a little better than choosing
!  some constant flow throughout each cavity block.
!
!  Further clarification of flow-solution cases:
!
!     (1) If both grids are accompanied by solution files, the solutions may
!         be either vertex-centered or cell-centered - it doesn't matter.
!         The operation is quite general - just the number of flow variables
!         must be the same in each solution file.
!
!     (2) If no solution file accompanies the second grid, then the first
!         grid and solution should both be cell-centered as output by the
!         POSTFLOW procedure for DPLR solutions.  The second grid, on the
!         other hand, is most conveniently vertex-centered.  In this case,
!         it is converted here to cell-centered form (with halos) prior to
!         transcribing to the output grid, and a corresponding flow solution
!         is generated as described above.
!
!  History:
!
!     05/04/04  D. Saunders  Initial implementation of COMBINE_GRIDS.
!     12/02/04   "     "     COMBINE_BLOCKS allows for solution files too,
!                            with specialized treatment for solution blocks
!                            that do not accompany the second grid file.
!     09/12/05   "     "     Slight change to allow for variable # species.
!     10/22/05   "     "     Postflow for DPLR V3.04.1 doesn't output a second
!                            temperature unless ivib=1, so allow 9 flow vars.
!                            rather than the previous 10 in the specialized
!                            option of q_block_estimate (cavity block case).
!     02/27/08   "     "     Ideal gas flow is now handled by the specialized
!                            option for initializing the second set of blocks.
!     06/21/13   "     "     For lee-side cavities, the hard-coded factor of
!                            0.1 x pressure is probably excessive.  This value
!                            is still the default, but a prompt allows it to be
!                            overridden.  EOF when reading a scale factor is
!                            handled, so existing scripts should still be OK.
!     08/05/13   "     "     All variants of the ADT routines have been merged
!                            into a module (generic build & search calls).
!
!  Author:  David Saunders, ELORET/NASA Ames Research Center, Moffett Field, CA
!                           Now ERC, Inc. at NASA/ARC.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_block_structure  ! Derived data type for one grid block
   use xyzq_io_module        ! PLOT3D file I/O package

   implicit none

!  Constants:

   integer, parameter :: &
      lunin1g = 1,       &
      lunin1f = 2,       &
      lunin2g = 3,       &
      lunin2f = 4,       &
      lunkbd  = 5,       &
      luncrt  = 6,       &
      lunoutg = 7,       &
      lunoutf = 8

   logical, parameter :: &
      false = .false.,   &
      true  = .true.

!  Variables:

   integer :: &
      i, ib, ios, j, nblocks_1, nblocks_2, nblocks_out, npts, nvars_1, nvars_2

   logical :: &
      cavity_case, formatted_1, formatted_2, formatted_out, solution_1,        &
      solution_2

   character :: &
      response * 1

!  Composite data types:

   type (grid_type), pointer, dimension (:) :: &
      grid_1, grid_2, grid_out, surf_1

!  Execution:

   call file_prompt (lunin1g, -luncrt, 'input grid #1', 'old', true,           &
                     formatted_1, ios)
   if (ios /= 0) go to 99

   write (luncrt, '(a)', advance='no') ' Accompanying flow solution? (y/n) '
   read  (lunkbd, '(a)') response
   solution_1 = response == 'y'

   if (solution_1) then
      call file_prompt (lunin1f, -luncrt, 'input flow #1', 'old', false,       &
                        formatted_1, ios)
      if (ios /= 0) go to 99
   end if

   call file_prompt (lunin2g, -luncrt, 'input grid #2', 'old', true,           &
                     formatted_2, ios)
   if (ios /= 0) go to 99

   if (solution_1) then
      write (luncrt, '(a)', advance='no') ' Accompanying flow solution? (y/n) '
      read  (lunkbd, '(a)') response
      solution_2 = response == 'y'
   else
      solution_2 = false
   end if

   cavity_case = solution_1 .and. .not. solution_2

   if (solution_2) then
      call file_prompt (lunin2f, -luncrt, 'input flow #2', 'old', false,       &
                        formatted_2, ios)
      if (ios /= 0) go to 99
   end if

!  Read the input file headers:

   call xyz_header_io (1, lunin1g, formatted_1, nblocks_1, grid_1, ios)
   if (ios /= 0) go to 99

   if (solution_1) then
      call q_header_io (1, lunin1f, formatted_1, nblocks_1, nvars_1, grid_1,   &
                        ios)
      if (ios /= 0) go to 99
   end if 

   call xyz_header_io (1, lunin2g, formatted_2, nblocks_2, grid_2, ios)
   if (ios /= 0) go to 99

   if (solution_2) then

      call q_header_io (1, lunin2f, formatted_2, nblocks_2, nvars_2, grid_2,   &
                        ios)
      if (ios /= 0) go to 99

      if (nvars_2 /= nvars_1) then
         write (luncrt, '(a, 2i4)') ' Mismatched function counts: ',           &
            nvars_1, nvars_2
         go to 99
      end if

   end if

   if (cavity_case) then ! We will generate flow for grid 2

!     Store the surface data for grid 1 (k = 1):

      allocate (surf_1(nblocks_1), stat=ios)
      if (ios /= 0) then
         write (luncrt,'(a)') ' Trouble allocating surface array for grid 1.'
         go to 99
      end if

      do ib = 1, nblocks_1
         surf_1(ib)%ni = grid_1(ib)%ni
         surf_1(ib)%nj = grid_1(ib)%nj
         surf_1(ib)%nk = 1
         surf_1(ib)%mi = grid_1(ib)%mi
         surf_1(ib)%mj = grid_1(ib)%mj
         surf_1(ib)%mk = 1
      end do

      nvars_2 = nvars_1

   end if

!  Set up the output file(s):

   call file_prompt (lunoutg, -luncrt, 'output grid', 'unknown', true,         &
                     formatted_out, ios)
   if (ios /= 0) go to 99

   if (solution_1) then
      call file_prompt (lunoutf, -luncrt, 'output flow', 'unknown', false,     &
                        formatted_out, ios)
      if (ios /= 0) go to 99
   end if

   nblocks_out = nblocks_1 + nblocks_2

   allocate (grid_out(nblocks_out))

   do ib = 1, nblocks_1
      grid_out(ib)%ni = grid_1(ib)%ni
      grid_out(ib)%nj = grid_1(ib)%nj
      grid_out(ib)%nk = grid_1(ib)%nk
      if (solution_1) then
         grid_out(ib)%mi = grid_1(ib)%mi
         grid_out(ib)%mj = grid_1(ib)%mj
         grid_out(ib)%mk = grid_1(ib)%mk
      end if
   end do

   i = nblocks_1
   do ib = 1, nblocks_2
      i = i + 1
      grid_out(i)%ni = grid_2(ib)%ni ! If not the cavity case
      grid_out(i)%nj = grid_2(ib)%nj
      grid_out(i)%nk = grid_2(ib)%nk
      if (solution_2) then
         grid_out(i)%mi = grid_2(ib)%mi
         grid_out(i)%mj = grid_2(ib)%mj
         grid_out(i)%mk = grid_2(ib)%mk
      else if (solution_1) then     ! Cavity case: assume vertex-centered grid 2
         grid_out(i)%ni = grid_2(ib)%ni + 1 ! g & f are cell-centered with halos
         grid_out(i)%mi = grid_2(ib)%ni + 1
         grid_out(i)%nj = grid_2(ib)%nj + 1
         grid_out(i)%mj = grid_2(ib)%nj + 1
         grid_out(i)%nk = grid_2(ib)%nk + 1
         grid_out(i)%mk = grid_2(ib)%nk + 1
      else ! Grid only
      end if
   end do

!  Write the output file header(s):

   call xyz_header_io (2, lunoutg, formatted_out, nblocks_out, grid_out, ios)
   if (ios /= 0) go to 99

   if (solution_1) then
      call q_header_io (2, lunoutf, formatted_out, nblocks_out, nvars_1,       &
                        grid_out, ios)
      if (ios /= 0) go to 99
   end if


!  Transcribe the blocks one at a time, the first grid first:
!  ----------------------------------------------------------

   do ib = 1, nblocks_1

      call xyz_allocate (grid_1(ib), ios)
      if (ios /= 0) go to 99

      npts = grid_1(ib)%ni * grid_1(ib)%nj * grid_1(ib)%nk

      call xyz_block_io (1, lunin1g, formatted_1, npts,                        &
                         grid_1(ib)%x, grid_1(ib)%y, grid_1(ib)%z, ios)
      if (ios /= 0) go to 99

      call xyz_block_io (2, lunoutg, formatted_out, npts,                      &
                         grid_1(ib)%x, grid_1(ib)%y, grid_1(ib)%z, ios)
      if (ios /= 0) go to 99

      if (solution_1 .and. .not. solution_2) then  ! Save surface data

         call xyz_allocate (surf_1(ib), ios)
         if (ios /= 0) go to 99

         do j = 1, grid_1(ib)%nj
            do i = 1, grid_1(ib)%ni
               surf_1(ib)%x(i,j,1) = grid_1(ib)%x(i,j,1)
               surf_1(ib)%y(i,j,1) = grid_1(ib)%y(i,j,1)
               surf_1(ib)%z(i,j,1) = grid_1(ib)%z(i,j,1)
            end do
         end do

      end if

      deallocate (grid_1(ib)%x, grid_1(ib)%y, grid_1(ib)%z, stat=ios)

      if (ios /= 0) then
         write (luncrt, '(a, 2i5)') &
            ' Trouble deallocating first grid block #', ib, ios
         go to 99
      end if

      if (solution_1) then

         call q_allocate (grid_1(ib), nvars_1, ios)
         if (ios /= 0) go to 99
 
         call q_block_io (1, lunin1f, formatted_1, nvars_1,                    &
                          grid_1(ib)%mi, grid_1(ib)%mj, grid_1(ib)%mk,         &
                          grid_1(ib)%q, ios)
         if (ios /= 0) go to 99

         call q_block_io (2, lunoutf, formatted_out, nvars_1,                  &
                          grid_1(ib)%mi, grid_1(ib)%mj, grid_1(ib)%mk,         &
                          grid_1(ib)%q, ios)
         if (ios /= 0) go to 99

         if (.not. solution_2) then ! Save surface data

            call q_allocate (surf_1(ib), nvars_1, ios)
            if (ios /= 0) go to 99

            do j = 1, grid_1(ib)%mj
               do i = 1, grid_1(ib)%mi
                  surf_1(ib)%q(:,i,j,1) = grid_1(ib)%q(:,i,j,1)
               end do
            end do

         end if

         deallocate (grid_1(ib)%q, stat=ios)
         if (ios /= 0) then
            write (luncrt, '(a, 2i5)') &
               ' Trouble deallocating first flow block #', ib, ios
            go to 99
         end if

      end if

   end do ! Next block from the first file(s)

   close (lunin1g)
   if (solution_1) close (lunin1f)


!  Append the blocks from the second grid:
!  ---------------------------------------

   do ib = 1, nblocks_2

      call xyz_allocate (grid_2(ib), ios)
      if (ios /= 0) go to 99

      npts = grid_2(ib)%ni * grid_2(ib)%nj * grid_2(ib)%nk

      call xyz_block_io (1, lunin2g, formatted_2, npts,                        &
                         grid_2(ib)%x, grid_2(ib)%y, grid_2(ib)%z, ios)
      if (ios /= 0) go to 99

      if (.not. cavity_case) then ! Simply transcribe

         call xyz_block_io (2, lunoutg, formatted_out, npts,                   &
                            grid_2(ib)%x, grid_2(ib)%y, grid_2(ib)%z, ios)
         if (ios /= 0) go to 99

      else ! Cavity case: recenter the grid first

         i = ib + nblocks_1

         call xyz_allocate (grid_out(i), ios)
         if (ios /= 0) go to 99

         call recenter_volume (grid_2(ib), grid_out(i))

         npts = grid_out(i)%ni * grid_out(i)%nj * grid_out(i)%nk

         call xyz_block_io (2, lunoutg, formatted_out, npts,                   &
                            grid_out(i)%x, grid_out(i)%y, grid_out(i)%z, ios)
         if (ios /= 0) go to 99

      end if

      deallocate (grid_2(ib)%x, grid_2(ib)%y, grid_2(ib)%z, stat=ios)
      if (ios /= 0) then
         write (luncrt, '(a, 2i5)') &
            ' Trouble deallocating second grid block #', ib, ios
         go to 99
      end if

      if (solution_2) then ! Simply transcribe

         call q_allocate (grid_2(ib), nvars_2, ios)
         if (ios /= 0) go to 99

         call q_block_io (1, lunin2f, formatted_2, nvars_2,                    &
                          grid_2(ib)%mi, grid_2(ib)%mj, grid_2(ib)%mk,         &
                          grid_2(ib)%q, ios)
         if (ios /= 0) go to 99

         call q_block_io (2, lunoutf, formatted_out, nvars_2,                  &
                          grid_2(ib)%mi, grid_2(ib)%mj, grid_2(ib)%mk,         &
                          grid_2(ib)%q, ios)
         if (ios /= 0) go to 99

         deallocate (grid_2(ib)%q, stat=ios)

         if (ios /= 0) then
            write (luncrt, '(a, 2i5)') &
               ' Trouble deallocating second grid q for block #', ib, ios
            go to 99
         end if

      end if

      if (cavity_case) then ! Generate flow starting guesses

         call q_allocate (grid_out(i), nvars_1, ios)
         if (ios /= 0) go to 99

         call q_block_estimate (nblocks_1, nvars_1, surf_1, grid_out(i), ib,   &
                                ios)
         if (ios /= 0) go to 99

         call q_block_io (2, lunoutf, formatted_out, nvars_1,                  &
                          grid_out(i)%mi, grid_out(i)%mj, grid_out(i)%mk,      &
                          grid_out(i)%q, ios)
         if (ios /= 0) go to 99

         deallocate (grid_out(i)%x, grid_out(i)%y, grid_out(i)%z,              &
                     grid_out(i)%q, stat=ios)
         if (ios /= 0) then
            write (luncrt, '(a, 2i5)') &
               ' Trouble deallocating second grid x,y,z,q for block #', ib, ios
            go to 99
         end if

      end if

   end do ! Next block from the second file(s)

   close (lunin2g)
   if (solution_2) close (lunin2f)
   close (lunoutg)
   if (solution_1) close (lunoutf)

99 continue

! *** stop ! Avoid system dependencies.

   end program combine_blocks

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine q_block_estimate (nblocks_1, nvars, surf, block_2, iblock, ios)

!  Estimate a flow solution for one block that is assumed to be a cavity block
!  near some surface patch representing a smooth OML grid.  For every point of
!  the cell-centered input grid block, the nearest cell of the given surface
!  grid is located with an efficient ADT search, and the surface solution is
!  interpolated within that cell.  Velocity components are reset to zero, and
!  species densities are reduced by a factor of (say) 10 (now an input).
!
!  Species are assumed to appear first among the flow variables, followed by 5
!  more variables (3 velocity components and 2 temperatures).  However, 9 or
!  fewer variables are interpreted as including only 1 temperature, with less
!  than 5 variables being illegal.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_block_structure
   use adt_utilities         ! All variants of the ADT build & search utilities

   implicit none

!  Arguments:

   integer, intent (in)             :: nblocks_1        ! # surface patches
   integer, intent (in)             :: nvars            ! # flow variables;
                                                        ! see below
   type (grid_type), intent (in)    :: surf(nblocks_1)  ! Grid 1 surface data
   type (grid_type), intent (inout) :: block_2          ! Input with a grid;
                                                        ! output with a flow
   integer, intent (in)             :: iblock           ! 1 => set up ADT here
   integer, intent (out)            :: ios              ! 0 means no error

!  Local constants:

   real,    parameter :: default = 0.1    ! For scaling species densities
   real,    parameter :: one     = 1.0
   real,    parameter :: zero    = 0.0
   logical, parameter :: true    = .true. ! Tells search routine to initialize
                                          ! closest distance to a large number
!  Local variables:

   integer :: i, ib, ic, iquad, j, jc, k, nspecies
   real    :: dmin, dmax, dsqmin, p, pm1, q, qm1
   real    :: foot_coords(3), target_coords(3)

!  Variables saved for further calls:

   real,    save :: pscale
   integer, save :: nquad                               ! # cells being searched
   integer, allocatable, dimension (:,:), save :: conn  ! (1:3,*) is used by the
                                                        ! ADT package for block
!                                                       ! # & lower left indices
!  Execution:

   if (nvars <= 9) then  ! Assume only 1 temperature
      nspecies = nvars - 4
   else                  ! Assume 2 temperatures
      nspecies = nvars - 5
   end if

   if (nspecies < 1) then
      write (*, '(/, a, 2i4, i5)') &
         ' Q_BLOCK_ESTIMATE: Data error likely.  nvars, nspecies, block #:',   &
         nvars, nspecies, iblock
      ios = -1
      go to 99
   end if

   if (iblock == 1) then

!     Allow for variable scaling of the densities/pressures:

      pscale = default
      write (*, '(a)', advance='no') &
         ' Scale factor for densities/pressure? [EOF => 0.1] '
      read (*, *, iostat=ios) pscale

!     Count the surface patch quad. cells and allocate work-space for ADT pkg.:

      nquad = 0
      do ib = 1, nblocks_1
         nquad = (surf(ib)%ni - 1) * (surf(ib)%nj - 1) + nquad
      end do

      allocate (conn(3,nquad))

!     Build the surface patch search tree:

      call build_adt (nblocks_1, surf, nquad, conn) 

      write (*, '(/, a, /)') ' Block #   Min. distance   Max. distance'

   end if

   dmin = 1.e+30 ! Report some statistics to help catch mismatched units, maybe
   dmax = -dmin

   do k = 1, block_2%mk
      do j = 1, block_2%mj
         do i = 1, block_2%mi
            target_coords(1) = block_2%x(i,j,k)
            target_coords(2) = block_2%y(i,j,k)
            target_coords(3) = block_2%z(i,j,k)

            call search_adt (target_coords, iquad, p, q, dsqmin, true,         &
                             nblocks_1, surf, nquad, conn, foot_coords)

!           Finding a surface cell closer than some tolerance is not possible
!           in this application.  Just record the nearest and farthest cases.

            dmin = min (dsqmin, dmin)
            dmax = max (dsqmin, dmax)

            ib = conn(1,iquad) ! Surface patch #
            ic = conn(2,iquad) ! Lower left indices for best quad. cell found
            jc = conn(3,iquad)

            pm1 = one - p
            qm1 = one - q

            block_2%q(:,i,j,k) = qm1 * (pm1 * surf(ib)%q(:,ic,  jc,  1)  +     &
                                          p * surf(ib)%q(:,ic+1,jc,  1)) +     &
                                   q * (pm1 * surf(ib)%q(:,ic,  jc+1,1)  +     &
                                          p * surf(ib)%q(:,ic+1,jc+1,1))
            block_2%q(1:nspecies,i,j,k) = pscale * block_2%q(1:nspecies,i,j,k)
            block_2%q(nspecies+1:nspecies+3,i,j,k) = zero    ! Velocities
         end do
      end do
   end do

   dmin = sqrt (dmin)
   dmax = sqrt (dmax)

   write (*, '(i8, 2f16.8)') iblock, dmin, dmax
   ios = 0

99 return

   end subroutine q_block_estimate

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
   use xyzq_io_module        ! External module

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
