!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      program reorder_blocks
!
!     Description:
!
!        REORDER_BLOCKS can reorder the blocks in a multiblock grid and an
!     optional PLOT3D-type function file.
!
!     Control file ('reorder_blocks.inp'):
!
!        This should contain a list of integers representing the desired block
!     order.  The list can span any number of lines except if the various short
!     forms handled by the RDLIST utility are employed, the list should be on a
!     single line.
!
!        1:99 104 103 102 100 101   (all on one line)   or

!        1 3 5 7 9 11 13 15
!        2 4 6 8 10 12 14           (any number of lines)
!
!     Procedures:
!
!        XYZQ_IO package  I/O utilities for PLOT3D grid and function files
!
!     History:
!
!        06/19/06  D. Saunders  Initial adaptation of EXTRACT_BLOCKS.
!
!     Author:  David Saunders, ELORET/NASA Ames Research Ctr., Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Modules:

      use grid_block_structure
      use xyzq_io_module

      implicit none

!     Constants:

      integer, parameter :: &
         lunctl     = 1,    &
         lunkbd     = 5,    &
         luncrt     = 6,    &
         lunxyz_in  = 7,    &
         lunq_in    = 8,    &
         lunxyz_out = 9,    &
         lunq_out   = 10

!     Variables:

      integer :: &
         i, ib, ii, ios, lunqi, lunqo, nblocks_in, nblocks_out, npts, num_q

      integer, allocatable, dimension (:) :: &
         iblock

      logical :: &
         cell_centered, formatted_in, formatted_out, qfile

      character :: &
         answer * 1

!     Composite data types:

      type (grid_type), pointer, dimension (:) :: &
         xyzq_in, xyzq_out

!     Execution:

      open (lunctl, file='reorder_blocks.inp', status='old', iostat=ios)

      if (ios /= 0) then
         write (luncrt, '(/, a)') &
            ' Unable to open reorder_blocks.inp control file.'
         go to 999
      end if

!     Allow for a shorthand form of the block order, but it must be on one line.

      ii = 0
      do ! Until EOF
         read (lunctl, *, iostat=ios)
         if (ios < 0) exit

         ii = ii + 1
      end do

      rewind (lunctl)

!     Find how many input blocks there are before reading the output order:

      call file_prompt (lunxyz_in, -luncrt, 'input grid', 'old', .true., &
                        formatted_in, ios)
      if (ios /= 0) go to 999

      write (luncrt, '(a)', advance='no') &
         ' Is there a function file too? [y/n]: '
      read  (lunkbd, *) answer

      qfile = answer == 'y' .or. answer == 'Y'

      if (qfile) then
         call file_prompt (lunq_in, -luncrt, 'input flow', 'old', .false., &
                           formatted_in, ios)
         if (ios /= 0) go to 999
         lunqi = lunq_in
         lunqo = lunq_out
      else
         lunqi = -1
         lunqo = -1
      end if

!     Read the input file header(s):

      call xyz_header_io (1, lunxyz_in, formatted_in, nblocks_in, xyzq_in, ios)

      if (ios /= 0) go to 999

      if (qfile) then
         call q_header_io (1, lunqi, formatted_in, nblocks_in, num_q, &
                           xyzq_in, ios)
         if (ios /= 0) go to 999
      end if

!     Display the block dimensions to reassure the user:

      write (luncrt, '(a)')
      write (luncrt, '(4(i6, 2X, 3i5))') &
         (ib, xyzq_in(ib)%ni, xyzq_in(ib)%nj, xyzq_in(ib)%nk, &
          ib = 1, nblocks_in)

!     Now that we know how many blocks should be listed, read the list:

      nblocks_out = nblocks_in

      allocate (iblock(nblocks_out))

      if (ii > 1) then ! Simple list-directed read (see above)

         read (lunctl, *, iostat=ios) iblock

         if (ios /= 0) then
            write (luncrt, '(/, a)') &
               ' Trouble reading desired block order (list-directed read).'
            go to 999
         end if

      else ! All blocks are listed on one line, possibly in shortened form:

         call rdlist (luncrt, ' ', lunctl, nblocks_out, iblock)
 
         if (nblocks_out /= nblocks_in) then
            write (luncrt, '(/, a, /, (a, i5))') &
               ' Trouble reading desired block order (RDLIST).', &
               ' # blocks expected:', nblocks_in, &
               ' # blocks in list: ', nblocks_out
            go to 999
         end if

      end if

      close (lunctl)

      write (luncrt, '(/, a, (16i5))') &
         ' Desired output block order:', iblock

!     Check the inputs for bad block numbers:

      do ii = 1, nblocks_out
         ib = iblock(ii)
         if (ib < 1 .or. ib > nblocks_in) then
            write (luncrt, '(a, i6, a, i5)') &
               ' Output block # is bad:', ib, ';  # blocks read:', nblocks_in
            go to 999
         end if
      end do

!     Set up the output file(s):

      call file_prompt (lunxyz_out, -luncrt, 'output grid', 'unknown', .true., &
                        formatted_out, ios)
      if (ios /= 0) go to 999

      if (qfile) then
         call file_prompt (lunq_out, -luncrt, 'output flow', 'unknown', &
                           .false., formatted_out, ios)
         if (ios /= 0) go to 999
      end if

      allocate (xyzq_out(nblocks_out))

      do i = 1, nblocks_out

         ib = iblock(i)
         xyzq_out(i)%ni = xyzq_in(ib)%ni
         xyzq_out(i)%nj = xyzq_in(ib)%nj
         xyzq_out(i)%nk = xyzq_in(ib)%nk

         if (qfile) then
            xyzq_out(i)%mi = xyzq_in(ib)%mi
            xyzq_out(i)%mj = xyzq_in(ib)%mj
            xyzq_out(i)%mk = xyzq_in(ib)%mk
         end if

      end do

      call xyz_header_io (2, lunxyz_out, formatted_out, nblocks_out, &
                          xyzq_out, ios)
      if (ios /= 0) go to 999

      if (qfile) then

         call q_header_io (2, lunq_out, formatted_out, nblocks_out, num_q, &
                           xyzq_out, ios)
         if (ios /= 0) go to 999
      end if

!     It's awkward to avoid storing all blocks.  Treat the grid first,
!     then treat the function file if it's present.

      do ib = 1, nblocks_in

         call xyz_allocate (xyzq_in(ib), ios)
         if (ios /= 0) go to 999

         npts = xyzq_in(ib)%ni * xyzq_in(ib)%nj * xyzq_in(ib)%nk

         call xyz_block_io (1, lunxyz_in, formatted_in, npts, &
                            xyzq_in(ib)%x, xyzq_in(ib)%y, xyzq_in(ib)%z, ios)
         if (ios /= 0) go to 999

      end do

      close (lunxyz_in)

      do ii = 1, nblocks_out

         ib = iblock(ii)

         npts = xyzq_in(ib)%ni * xyzq_in(ib)%nj * xyzq_in(ib)%nk

         call xyz_block_io (2, lunxyz_out, formatted_out, npts, &
                            xyzq_in(ib)%x, xyzq_in(ib)%y, xyzq_in(ib)%z, ios)
         if (ios /= 0) go to 999

         deallocate (xyzq_in(ib)%x, xyzq_in(ib)%y, xyzq_in(ib)%z, stat=ios)

         if (ios /= 0) then
            write (luncrt, '(a, 2i5)') &
               ' Trouble deallocating input grid block #', ib, ios
            go to 999
         end if

      end do

      close (lunxyz_out)

      if (qfile) then ! Read all blocks

         do ib = 1, nblocks_in

            call q_allocate (xyzq_in(ib), num_q, ios)
            if (ios /= 0) go to 999

            call q_block_io (1, lunq_in, formatted_in, num_q, &
                             xyzq_in(ib)%mi, xyzq_in(ib)%mj, xyzq_in(ib)%mk, &
                             xyzq_in(ib)%q, ios)
            if (ios /= 0) go to 999

        end do

        close (lunq_in)

!       Write the blocks in the desired order:

        do ii = 1, nblocks_out

            ib = iblock(ii)

            call q_block_io (2, lunq_out, formatted_out, num_q, &
                             xyzq_in(ib)%mi, xyzq_in(ib)%mj, xyzq_in(ib)%mk, &
                             xyzq_in(ib)%q, ios)
            if (ios /= 0) go to 999

            deallocate (xyzq_in(ib)%q, stat=ios)

            if (ios /= 0) then
               write (luncrt, '(a, 2i5)') &
                  ' Trouble deallocating input flow block #', ib, ios
               go to 999
            end if

         end do

      end if

      close (lunq_in)

  999 continue

! *** stop ! Avoid system dependencies.

      end program reorder_blocks
