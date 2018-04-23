!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   program tet_interp

!  This is an adaptation of the earlier FLOW_INTERP, intended to interpolate 3-space flow field data from an unstructured (tetra-
!  hedral) volume grid on to a multiblock structured grid.  It is prompted by a need to interpolate a DSMC (Direct Simulation Monte
!  Carlo) solution containing (at least) pressure, temperature, and electron number density onto a hemisphere line of sight grid for
!  program Radio_Attenuation purposes.   It is more general than that:  the flow data may contain any number of variables and the
!  target grid may be any multiblock structured grid.
!
!  Function values are assumed to be vertex-centered as needed for trilinear interpolation within tetrahedral cells.  Deal with
!  cell-centered function values only if called for (via volume-weighted averaging?).
!
!  The ADT search scheme employed presently handles only one zone of cells.  Generalize this if called for< in adt_utilities.f90.
!
!  As with FLOW_INTERP, all of the flow variables present are interpolated.  If necessary, use EXTRACT_FUNCTIONS on the result to
!  extract some subset of variables.  Similarly, the target grid presumably matches the unstructured grid in some sense.  Any mis-
!  matches towards the outer boundary will not cause interpolation trouble in that the closest point of the best possible cell is
!  always found by the ADT (Alternating Digital Tree) search scheme employed.  Extrapolation is explicitly avoided.
!
!  Control inputs are all via prompts.
!
!  History:
!
!     08/05/13  DAS  FLOW_INTERP version from which TET_INTERP is being adapted.  The ADT variants had just been merged into a
!                    module with generic build and search calls.
!     08/15/16-  "   Initial adaptation of FLOW_INTERP, for vertex-centered tetrahedral data, with radio attenuation calculations
!     08/22/16       in mind at lower-density conditions than conventional CFD can handle.
!
!  Author:
!
!     David Saunders, AMA, Inc./NASA Ames Research Center, Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:
!  !!!!!!!!

   use tri_header_structure  ! Part of triangulation_io.f90
   use tri_zone_structure    ! Likewise
   use triangulation_io      ! Unstructured data file I/O
   use grid_block_structure  ! For on structured grid block
   use xyzq_io_module        ! For PLOT3D files

   implicit none

!  Local constants:
!  !!!!!!!!!!!!!!!!

   integer, parameter ::   &
      lunsol = 2,          & ! Input solution file
      luntar = 3,          & ! Input target grid
      lunf   = 4,          & ! Output flow solution
      lunkbd = 5,          & ! Keyboard entries
      luncrt = 6             ! Prompts and diagnostics

   character (11), parameter :: &
      format = 'unformatted'

!  Derived data types:
!  !!!!!!!!!!!!!!!!!!!

   type (tri_header_type)                   :: header
   type (tri_type),  pointer, dimension (:) :: solution
   type (grid_type), pointer, dimension (:) :: target

!  Local variables:
!  !!!!!!!!!!!!!!!!

   integer        :: i1, ib, ios, method, mi, mj, mk, n, nzones, nblocks_target, ni, nj, nk, npts, numf
   real           :: time_end, time_start
   logical        :: cleanup, formatted_out, formatted_sol, formatted_tar, initialize
   character  (1) :: answer
   character (80) :: filename

!  Execution:
!  !!!!!!!!!!

   write (luncrt, '(/, a)', advance='no') 'Input unstructured volume dataset: '
   read  (lunkbd, *) header%filename
   header%formatted = .true.  ! Unformatted tet. reads are not handled
   header%nvertices = 4       ! Tets., not triangles
   header%fileform  = 1       ! Vertex-centered is required

   write (luncrt, '(a)', advance='no') 'Input structured target volume grid: '
   read  (lunkbd, *) filename
   call determine_grid_form (filename, luntar, formatted_tar, ios);  if (ios /= 0) go to 999

   i1 = 1; if (formatted_tar) i1 = 3
   open (luntar, file=filename, status='old', form=format(i1:11), iostat=ios)
   call xyz_header_io (1, luntar, formatted_tar, nblocks_target, target, ios);  if (ios /= 0) go to 999

   write (luncrt, '(a)', advance='no') 'Output interpolated flow dataset: '
   read  (lunkbd, *) filename
   write (luncrt, '(a)', advance='no') 'Formatted (y) or unformatted (n)? '
   read  (lunkbd, *) answer
   formatted_out = answer == 'y'
   i1 = 1; if (formatted_out) i1 = 3
   open (lunf, file=filename, status='unknown', form=format(i1:11), iostat=ios)
   if (ios /= 0) then
      write (luncrt, '(/, 2a)') ' Unable to open output flow file: ', filename(1:len_trim(filename));  go to 999
   end if

   call vol_read (lunsol, header, solution, ios);  if (ios /= 0) go to 999   !  Read the entire input solution
   
   numf   = header%numf
   nzones = header%nzones
   if (nzones > 1) then
      write (luncrt, '(/, a)') 'More than one zone of tetrahedra requires generalization in adt_utilities.f90.'
      go to 999
   end if

   target(:)%mi = target(:)%ni                                               !  Write the interpolated flow solution file header
   target(:)%mj = target(:)%nj
   target(:)%mk = target(:)%nk
   call q_header_io (2, lunf, formatted_out, nblocks_target, numf, target, ios);  if (ios /= 0) go to 999

!  Perform the flow field interpolation one target block at a time:
!  ----------------------------------------------------------------

   call cpu_time (time_start)

   do ib = 1, nblocks_target

      initialize = ib == 1
      cleanup    = ib == nblocks_target

      call xyz_allocate (target(ib), ios)                                    ! Set up for and read the current target grid block
      ni = target(ib)%ni
      nj = target(ib)%nj
      nk = target(ib)%nk
      npts = ni * nj * nk

      if (ios /= 0) then
         write (luncrt, '(/, a, i4, a, i4)') ' Trouble allocating target grid.  Block #:', ib, '  I/O status:', ios
         write (luncrt, '(a, 3i5)') ' Dimensions: ', ni, nj, nk
         go to 999
      end if

      call xyz_block_io (1, luntar, formatted_tar, npts, target(ib)%x, target(ib)%y, target(ib)%z, ios)

      if (ios /= 0) then
         write (luncrt, '(/, a, i4, a, i4)') ' Trouble reading target grid.  Block #:', ib, '  I/O status:', ios
         write (luncrt, '(a, 3i5)') ' Dimensions: ', ni, nj, nk
         go to 999
      end if

      call q_allocate (target(ib), numf, ios)                                ! Allocate the corresponding interpolated flow block

      if (ios /= 0) then
         write (luncrt, '(/, a, i4, a, i4)') ' Trouble allocating interpolated flow space.  Block #:', ib, '  I/O status:', ios
         write (luncrt, '(a, 4i5)') ' Dimensions: ', ni, nj, nk, numf
         go to 999
      end if

!     Interpolate the flow at the points of the current target block.

      call interp_q (header, solution(1), ib, target(ib), initialize, cleanup, luncrt)
!     -------------

!     Save results for this target block:

      mi = target(ib)%mi
      mj = target(ib)%mj
      mk = target(ib)%mk

      call q_block_io (2, lunf, formatted_out, numf, mi, mj, mk, target(ib)%q, ios)

      if (ios /= 0) then
         write (luncrt, '(/, a, i4, a, i4)') ' Trouble writing interpolated flow solution.  Block #:', ib, '  I/O status:', ios
         write (luncrt, '(a, 4i5)') ' Dimensions and # functions: ', mi, mj, mk, numf
         go to 999
      end if

      deallocate (target(ib)%x, target(ib)%y, target(ib)%z, target(ib)%q)

   end do ! Next target block

   call cpu_time (time_end)

   write (luncrt, '(/, a, f9.2)') ' Total run time, seconds:', time_end - time_start

   close (luntar)
   close (lunf)

!  Done.

   999 continue  ! Avoid system-dependent STOP behavior

   end program tet_interp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine interp_q (header, solution, iblock_target, target_block, initialize, cleanup, luncrt)
!
!  This routine interpolates the given unstructured flow solution onto one block of a structured target grid.
!  It is the argument-driven, reusable portion of program tet_interp, q.v.
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:
!  !!!!!!!!

   use tri_header_structure  ! Part of triangulation_io.f90
   use tri_zone_structure    ! Likewise
   use triangulation_io      ! Unstructured data file I/O
   use grid_block_structure  ! Derived data type for one structured grid block
   use adt_utilities         ! ADT search package (all variants)
   use xyzq_io_module        ! For PLOT3D file I/O

   implicit none

!  Arguments:
!  !!!!!!!!!!

   type (tri_header_type), intent (inout) :: header         ! Input are numf, nvertices, etc.; xmin, xmax, etc. are assigned here
   type (tri_type),        intent (inout) :: solution       ! Solution grid and flow data, etc. input; single zone assumed for now
   integer,                intent (in)    :: iblock_target  ! Current target block number, for the running log
   type (grid_type),       intent (inout) :: target_block   ! Current target block, input with its grid defined (only),
                                                            ! and output with its flow field assigned via interpolation
                                                            ! or from nearest boundary point
   logical,                intent (in)    :: initialize     ! True means allocate solution-specific work-space (normally
                                                            ! when iblock_target = 1, but we allow for more than one set
                                                            ! of solution blocks in the same run)
   logical,                intent (in)    :: cleanup        ! True means deallocate solution-specific work-space before returning
   integer,                intent (in)    :: luncrt         ! Logical unit for running log and error messages

!  Local constants:
!  ----------------

   real,    parameter :: big  = 1.e+32, zero = 0.0
   logical, parameter :: true = .true.

!  Local variables:
!  !!!!!!!!!!!!!!!!

   integer :: i, i1, i2, i3, i4, icell, j, k, ni, ninside, nj, nk, nout
   logical :: profile
   real    :: dsqmin, dsqmax, time1, time2, unused, wall_distance
   real, dimension (3) :: delta_coords, interp_coords, target_coords, wall_coords
   real, dimension (4) :: pqrs

!  Variables saved for further calls:
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   character (36), save :: pformat
   integer,        save :: ncells, nnodes
   real,           save :: dsqtol

!  Execution:
!  !!!!!!!!!!

!  Pre-process the solution blocks?

   if (initialize) then

      nnodes = solution%nnodes     ! Only one zone assumed for now
      ncells = solution%nelements

!     Handle any number of functions in profile tabulations:

      pformat = '(i4, 5i11, 4f12.8, es15.7, nfes12.4)'
      write (pformat(28:29), '(i2)') header%numf

!     Set up a tolerance for inside/outside searches:

      header%xmin = big;  header%xmax = -big
      header%ymin = big;  header%ymax = -big
      header%zmin = big;  header%zmax = -big

      do i = 1, ncells
         header%xmin = min (header%xmin, solution%xyz(1,i));  header%xmax = max (header%xmax, solution%xyz(1,i))
         header%ymin = min (header%ymin, solution%xyz(2,i));  header%ymax = max (header%ymax, solution%xyz(2,i))
         header%zmin = min (header%zmin, solution%xyz(3,i));  header%zmax = max (header%zmax, solution%xyz(3,i))
      end do

      dsqtol = (10.*epsilon (dsqtol) * max (header%xmax - header%xmin, header%ymax - header%ymin, header%zmax - header%zmin))**2

!     Build an Alternating Digital Tree for all cells of all solution zones (just one, however, unless more are ever needed).

      call cpu_time (time1)

      call build_adt (nnodes, ncells, solution%conn, solution%xyz, unused)

      call cpu_time (time2)

      write (luncrt, '(/, a, 2i10, f8.2)') 'Search tree # pts., # cells, and tree build time, s:', nnodes, ncells, time2 - time1

   end if        ! End of initialization

   ninside = 0   ! # target block points found inside a solution cell
   nout    = 0   ! # such points found more than the distance tolerance from the nearest solution cell

   ni = target_block%ni
   nj = target_block%nj
   nk = target_block%nk

   profile = ni == 1 .and. nj == 1 ! A single radial line gets special treatment for wind tunnel data comparison or lines of sight

   if (profile) then
      write (luncrt, '(/, a, i4, /, 2a)') '! Profile #:', iblock_target, &
         '!  k       Tet.         i1         i2         i3         i4           p           q           r           s', &
         '  Wall distance  Interpolated flow'
   else
      write (luncrt, '(/, a, i4, a, 3i5, a, i8)') &
         ' Beginning target grid block', iblock_target, '.  Dimensions:', ni, nj, nk, '  Total:', ni * nj * nk
   end if

   dsqmax = zero

   call cpu_time (time1)

!  For each point of the current target block:

   do k = 1, nk
      do j = 1, nj
         do i = 1, ni

            target_coords(1) = target_block%x(i,j,k)
            target_coords(2) = target_block%y(i,j,k)
            target_coords(3) = target_block%z(i,j,k)

            call search_adt (target_coords, icell, pqrs, dsqmin, true, nnodes, ncells, solution%conn, solution%xyz, interp_coords)

            delta_coords(:) = target_coords(:) - interp_coords(:)

!!!         write (6, '(a, 4i4, a, i11, 4f12.8, a, 3es15.7)') &
!!!            ' ib,i,j,k:', iblock_target, i, j, k, ' tet, pqrs:', icell, pqrs(:), ' dxyz:', delta_coords

            dsqmax = max (dsqmin, dsqmax)

            if (dsqmin < dsqtol) then ! The best solution cell found was within the distance tolerance
               ninside = ninside + 1
            else
               nout = nout + 1
!!!            write (luncrt, '(a, 4i4, a, es13.5)') &
!!!               ' Tolerance exceeded for new (i,j,k,ib) =', i, j, k, iblock_target, '  min. distance:', sqrt (dsqmin)
            end if

!           Use the fractional coordinates to interpolate the flow solution to the target point:

            i1 = solution%conn(1,icell)
            i2 = solution%conn(2,icell)
            i3 = solution%conn(3,icell)
            i4 = solution%conn(4,icell)

            target_block%q(:,i,j,k) = pqrs(1)*solution%f(:,i1) + pqrs(2)*solution%f(:,i2) + &
                                      pqrs(3)*solution%f(:,i3) + pqrs(4)*solution%f(:,i4)
            if (profile) then
               if (k == 1) then
                  wall_coords(:) = pqrs(1)*solution%xyz(:,i1) + pqrs(2)*solution%xyz(:,i2) + &
                                   pqrs(3)*solution%xyz(:,i3) + pqrs(4)*solution%xyz(:,i4)
               end if
               wall_distance = sqrt ((wall_coords(1) - target_coords(1))**2 + &
                                     (wall_coords(2) - target_coords(2))**2 + &
                                     (wall_coords(3) - target_coords(3))**2)
               write (luncrt, pformat) k, icell, i1, i2, i3, i4, pqrs(:), wall_distance, target_block%q(:,1,1,k) 
            end if

         end do ! Next i for this target block
      end do ! Next j for this target block

!!!   if (.not. profile) then
!!!      write (luncrt, '(a, i4, a, 3i5, a, 2i9)') &
!!!         ' Target block:', iblock_target, '   point:', i-1,j-1,k, '  in, out:', ninside, nout
!!!   end if

   end do ! Next k for this target block


   if (.not. profile) then
      call cpu_time (time2)
      write (luncrt, '(a, i4, a, i9, i8, a, es13.5, a, f8.2)') &
         ' Finished  target grid block', iblock_target,  ';  in, out:', ninside, nout, &
         ';  largest distance from target:', sqrt (dsqmax), '  CPU seconds:', time2 - time1
   end if

!  Clean up?  (This allows for another set of solution blocks in the same run.)

   if (cleanup) then

   end if

   end subroutine interp_q
