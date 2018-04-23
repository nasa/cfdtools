!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   program surface_pad

!  For a structured surface grid (one or more patches), apply 1-dimensional
!  interpolations versus arc length to pad the number(s) of points.  Existing
!  program MSPLINE should be able to do this, but it seems to fail for surface
!  grids, and it preserves relative arc lengths where one might prefer roughly
!  uniform distributions as offered (only) here, at least initially.
!
!  Either or both index directions may be treated, patch by patch.  If the i
!  direction is specified, it is treated first.  Those results are further
!  interpolated in the j direction if so-specified.  (Program ADJUST_GRID can
!  transpose indices if some other order is desired.)
!
!  Two-dimensional surface interpolations are avoided, since the nonlinear ones
!  are dubious at best in the author's experience.
!
!  Modest numbers of data points are anticipated, for which existing subroutine
!  uniform_edge is appropriate (although its linear interpolation option for
!  specified index intervals is unlikely to be made use of here).
!
!  Input and output files are in PLOT3D format, ASCII or unformatted.
!
!  03/30/09  D.A.Saunders  Initial implementation to pad a mesh of wing leading
!                          edge sensor locations.
!  07/30/09    "     "     Densifying our finest Shuttle surface grid by factors
!                          of 4 in each direction showed that the new number of
!                          points for a multiple of m should be m*(n-1) + 1.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use grid_block_structure     ! Derived data type for one grid block
   use surface_patch_utilities  ! For transposing patches to treat j directions
   use xyzq_io_module           ! PLOT3D file I/O package

   implicit none

!  Local constants:

   integer, parameter :: &
      lunin  = 1,        &      ! Input grid
      lunout = 2,        &      ! Output grid
      lunkbd = 5,        &      ! Keyboard
      luncrt = 6                ! Screen

   logical, parameter :: &
      false  = .false.,  &
      true   = .true.

!  Local variables:

   integer :: &
      i, ib, idirection, ios, j, l, multiple(2), nblocks, nf, &
      npts_in(2), npts_out(2)
   integer, allocatable :: &
      md(:)                   ! Potentially forces linear interp. in places
   logical :: &
      cell_centered, cr, direction(2), eof, formatted_in, formatted_out
   character :: &
      index_direction(2)*1, interp_method(2)*1, prompt*64
   type (grid_type), pointer :: &
      grid_in(:), grid_out(:)

!  Execution:

   formatted_in  = true  ! Activates the prompt for ASCII or not

   call file_prompt (lunin,  -1, 'input grid',  'old',     true, &
                     formatted_in,  ios)
   if (ios /= 0) go to 99

   formatted_out = true

   call file_prompt (lunout, -1, 'output grid', 'unknown', true, &
                     formatted_out, ios)
   if (ios /= 0) go to 99

!  Storing all patches is lazy but memory is cheap these days ...

   call xyzq_read (lunin, -1, formatted_in, nblocks, nf, cell_centered, &
                   grid_in, ios)
   if (ios /= 0) go to 99

   index_direction(1) = 'i';  npts_in(1) = grid_in(1)%ni
   index_direction(2) = 'j';  npts_in(2) = grid_in(1)%nj

   do idirection = 1, 2

      direction(idirection) = true
      prompt = 'Interpolate in the ' // index_direction(idirection) // &
               ' direction?'
      l = len_trim (prompt) + 1
      call ready (luncrt, prompt(1:l), lunkbd, direction(idirection), &
                  cr, eof)
      if (eof) go to 99
      if (.not. direction(idirection)) cycle

      if (nblocks == 1) then
         npts_out(idirection) = 2*npts_in(idirection) - 1
         call readi (luncrt, 'Output number [<cr> = 2*n - 1]: ', &
                     lunkbd, npts_out(idirection), cr, eof)
      else
         multiple(idirection) = 2
         call readi (luncrt, 'Multiply point count by [<cr> = 2]: ', &
                     lunkbd, multiple(idirection), cr, eof)
      end if
      if (eof) go to 99

      interp_method(idirection) = 'M'  ! Monotonic cubic spline fit
      call readc (luncrt, 'Spline method [B|M|L; <cr> = M (monotonic)]: ', &
                  lunkbd, interp_method(idirection), cr, eof)
   end do

   allocate (grid_out(nblocks))

   if (direction(1)) then

      do ib = 1, nblocks

         npts_in(1) = grid_in(ib)%ni

         allocate (md(npts_in(1)));  md(:) = 3  ! 1 would force linear interp.

         if (nblocks > 1) npts_out(1) = multiple(1)* (npts_in(1) - 1) + 1

         grid_out(ib)%ni = npts_out(1)
         grid_out(ib)%nj = grid_in(ib)%nj
         grid_out(ib)%nk = 1

         allocate (grid_out(ib)%x(npts_out(1),grid_in(ib)%nj,1), &
                   grid_out(ib)%y(npts_out(1),grid_in(ib)%nj,1), &
                   grid_out(ib)%z(npts_out(1),grid_in(ib)%nj,1), stat=ios)

         if (ios /= 0) then
            write (luncrt, '(/, a, i4)') &
               '*** Trouble allocating output block', ib
            go to 99
         end if

         do j = 1, grid_in(ib)%nj

            call uniform_edge (npts_in(1),            grid_in(ib)%x(1,j,1),    &
                               grid_in(ib)%y(1,j,1),  grid_in(ib)%z(1,j,1),    &
                               md, interp_method(1),  npts_out(1),             &
                               grid_out(ib)%x(1,j,1), grid_out(ib)%y(1,j,1),   &
                               grid_out(ib)%z(1,j,1))
         end do

         deallocate (md)

         if (direction(2)) then  ! Transcribe interim result, deallocating

            call update_patch (grid_out(ib), 0, true, true, grid_in(ib))

         end if

      end do  ! Next block

   end if

   if (direction(2)) then

      do ib = 1, nblocks

         call transpose_patch (grid_in(ib), 0)

         npts_in(2) = grid_in(ib)%ni

         allocate (md(npts_in(2)));  md(:) = 3  ! 1 would force linear interp.

         if (nblocks > 1) npts_out(2) = multiple(2)* (npts_in(2) - 1) + 1

         grid_out(ib)%ni = npts_out(2)
         grid_out(ib)%nj = grid_in(ib)%nj
         grid_out(ib)%nk = 1

         allocate (grid_out(ib)%x(npts_out(2),grid_in(ib)%nj,1), &
                   grid_out(ib)%y(npts_out(2),grid_in(ib)%nj,1), &
                   grid_out(ib)%z(npts_out(2),grid_in(ib)%nj,1))

         do j = 1, grid_in(ib)%nj

            call uniform_edge (npts_in(2),            grid_in(ib)%x(1,j,1),    &
                               grid_in(ib)%y(1,j,1),  grid_in(ib)%z(1,j,1),    &
                               md, interp_method(2),  npts_out(2),             &
                               grid_out(ib)%x(1,j,1), grid_out(ib)%y(1,j,1),   &
                               grid_out(ib)%z(1,j,1))
         end do

         deallocate (md)

         call transpose_patch (grid_out(ib), 0)

      end do  ! Next block

   end if

!  Save the results:

   call xyz_write (lunout, formatted_out, nblocks, grid_out, ios)

   do ib = 1, nblocks
      deallocate (grid_in(ib)%x,  grid_in(ib)%y,  grid_in(ib)%z, &
                  grid_out(ib)%x, grid_out(ib)%y, grid_out(ib)%z)
   end do
   deallocate (grid_in, grid_out)

99 continue

   end program surface_pad
