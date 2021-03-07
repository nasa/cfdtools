!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program thin_grid
!
!  Description:
!
!     THIN_GRID extracts a subset of points from all blocks of a multiblock
!     grid.  This version has the option to do different things to different
!     blocks, as might be needed to make a grid point-to-point matched.
!
!     In order to avoid affecting existing scripts, the more general option is
!     implemented via an optional 'thin_grid.inp.2' control.  (There is no
!     'thin_grid.inp' but the '.2' is consistent with ancillary control files
!     for other grid utilities.)
!
!  Optional 'thin_grid.inp.2' Control File Format:
!
!     1 2 2 1    ! Halve the i & j density for block 1; leave nk alone
!     2 2 2 1    ! ..................................2................
!     3 1 1 1    ! Leave block 3 alone
!     .......    ! One line per block
!
!     (Later:)  In order to avoid duplicate points along surface patch
!     edges when setting up body points for radiation calculations, another
!     generalized capability has been implemented via thin_grid.inp.3 that
!     allows the initial index to be something other than 1 (same for all
!     blocks, since it's probably a surface grid where avoiding common edge
!     points makes some sense for body point purposes, as in the case of
!     a calculation on an asteroid-like shape).
!
!  Optional 'thin_grid.inp.3' Control File Format (cell-centered surface):
!
!     i1, inc  j1, jnc  k1, knc
!      3   4    3   4    1   1
!
!  Procedures:
!
!     XYZQ_IO package  I/O utilities for PLOT3D grid and function files
!
!  History:
!
!     04/29/04  D. Saunders  Initial adaptation of EXTRACT_BLOCKS.
!     01/29/11   "      "    More general functionality via thin_grid.inp.2.
!     06/03/15   "      "    Provided for avoiding duplicate edge points via
!                            thin_grid.inp.3.
!
!  Author:  David Saunders, ELORET Corp./NASA Ames, Moffett Field, CA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_block_structure
   use xyzq_io_module

   implicit none

!  Constants:

   integer, parameter :: &
      lunin  = 1,        &
      lunctl = 2,        &
      lunout = 3,        &
      lunkbd = 5,        &
      luncrt = 6

   logical, parameter :: &
      true = .true.

!  Variables:

   integer :: &
      i, i1, ib, ii, inc, ios, j, j1, jj, jnc, k, k1, kk, knc, nblocks, &
      ni, nj, nk, npts

   logical :: &
      formatted_in, formatted_out, i1_is_1, same_for_all_blocks

   integer, allocatable, dimension (:) :: &
      iincr, jincr, kincr

!  Derived data types:

   type (grid_type), pointer, dimension (:) :: &
      gridin, gridout

!  Execution:
!  ----------

   open (lunctl, file='thin_grid.inp.2', status='old', iostat=ios)
   same_for_all_blocks = ios /= 0

   i1 = 1;  j1 = 1;  k1 = 1
   if (same_for_all_blocks) then
      close (lunctl)
      open (lunctl, file='thin_grid.inp.3', status='old', iostat=ios)
      i1_is_1 = ios /= 0
      if (.not. i1_is_1) then
         read  (lunctl, '(a)')  ! Header
         read  (lunctl, *) i1, inc, j1, jnc, k1, knc
         close (lunctl)
      end if
   else
      i1_is_1 = true
   end if

!  Prompt for the input grid and read its header records:

   call file_prompt (lunin,  -luncrt, 'input grid', 'old', true, &
                     formatted_in, ios)
   if (ios /= 0) go to 99

   call xyz_header_io (1, lunin, formatted_in, nblocks, gridin, ios)
   if (ios /= 0) go to 99

!  Prompt for the output grid and set up its header records:

   call file_prompt (lunout, -luncrt, 'output grid', 'unknown', true, &
                     formatted_out, ios)
   if (ios /= 0) go to 99

   allocate (gridout(nblocks))
   allocate (iincr(nblocks), jincr(nblocks), kincr(nblocks))

   if (same_for_all_blocks) then

      if (i1_is_1) then
         write (luncrt, '(a)', advance='no') &
            ' Increments to use for thinning i, j and k: '
         read  (lunkbd, *) inc, jnc, knc
      end if

      iincr(:) = inc
      jincr(:) = jnc
      kincr(:) = knc

   else

      do ib = 1, nblocks
         read (lunctl, *, iostat=ios) i, iincr(ib), jincr(ib), kincr(ib)
         if (ios /= 0) then
            write (luncrt, '(a, i5)') &
               'Trouble reading thinning controls.  Block #:', ib
            go to 99
         end if
      end do
      close (lunctl)

   end if

   do ib = 1, nblocks
      inc = iincr(ib);  gridout(ib)%ni = (gridin(ib)%ni + inc - i1) / inc
      jnc = jincr(ib);  gridout(ib)%nj = (gridin(ib)%nj + jnc - j1) / jnc
      knc = kincr(ib);  gridout(ib)%nk = (gridin(ib)%nk + knc - k1) / knc
   end do

   call xyz_header_io (2, lunout, formatted_out, nblocks, gridout, ios)
   if (ios /= 0) go to 99

!  Thin the blocks one at a time:

   do ib = 1, nblocks

      call xyz_allocate (gridin(ib), ios)
      if (ios /= 0) go to 99

      ni = gridin(ib)%ni;  nj = gridin(ib)%nj;  nk = gridin(ib)%nk
      npts = ni * nj * nk

      call xyz_block_io (1, lunin, formatted_in, npts, &
                         gridin(ib)%x, gridin(ib)%y, gridin(ib)%z, ios)
      if (ios /= 0) go to 99

      call xyz_allocate (gridout(ib), ios)
      if (ios /= 0) go to 99

      inc = iincr(ib)
      jnc = jincr(ib)
      knc = kincr(ib)

      kk = 0
      do k = k1, nk, knc
         kk = kk + 1
         jj = 0
         do j = j1, nj, jnc
            jj = jj + 1
            ii = 0
            do i = i1, ni, inc
               ii = ii + 1
               gridout(ib)%x(ii,jj,kk) = gridin(ib)%x(i,j,k)
               gridout(ib)%y(ii,jj,kk) = gridin(ib)%y(i,j,k)
               gridout(ib)%z(ii,jj,kk) = gridin(ib)%z(i,j,k)
            end do
         end do
      end do

      npts = gridout(ib)%ni * gridout(ib)%nj * gridout(ib)%nk

      call xyz_block_io (2, lunout, formatted_out, npts, &
                         gridout(ib)%x, gridout(ib)%y, gridout(ib)%z, ios)
      if (ios /= 0) go to 99
         
      deallocate (gridin(ib)%x,  gridin(ib)%y,  gridin(ib)%z,  stat=ios)
      deallocate (gridout(ib)%x, gridout(ib)%y, gridout(ib)%z, stat=ios)

   end do

99 continue

! *** stop ! Avoid system dependencies.

   end program thin_grid
