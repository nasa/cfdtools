!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      program extract_blocks
!
!     Description:
!
!        EXTRACT_BLOCKS extracts a list of blocks from a multiblock grid and an
!     optional "q" file (or rather a PLOT3D-type function file).
!
!        This version can handle "iblanked" files at the cost of another prompt.
!
!     Control file ('extract_blocks.inp'):
!
!        This should contain a list of ascending block numbers, one per line.
!     E.g.:
!
!        1
!        3
!        7
!        :
!        15
!
!     Procedures:
!
!        XYZQ_IO package  I/O utilities for PLOT3D grid and function files
!
!     History:
!
!        04/06/04  D.A.Saunders  Initial adaptation of GRID_FACES.
!        07/02/10    "      "    Extended to handle OVERFLOW solutions.
!
!     Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA
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
         cell_centered, formatted_in, formatted_out, iblanked, qfile

      character :: &
         answer * 1

!     Composite data types:

      type (grid_type), pointer, dimension (:) :: &
         xyzq_in, xyzq_out

!     Execution:

      open (lunctl, file='extract_blocks.inp', status='old', iostat=ios)

      if (ios /= 0) then
         write (luncrt, '(/, a)') &
            ' Unable to open extract_blocks.inp control file.'
         go to 999
      end if

      nblocks_out = 0

      do ! Until EOF
         read (lunctl, *, iostat=ios) ib
         if (ios == 0) then
            nblocks_out = nblocks_out + 1
         else if (ios < 0) then
            exit
         else
            write (luncrt, '(a, i5)') &
               ' Bad integer in control file.  Line #: ', nblocks_out + 1
            go to 999
         end if
      end do

      rewind (lunctl)

      allocate (iblock(nblocks_out))

      read (lunctl, *) iblock

      close (lunctl)

      call file_prompt (lunxyz_in, -luncrt, 'input grid', 'old', .true., &
                        formatted_in, ios)
      if (ios /= 0) go to 999

      write (luncrt, '(a)', advance='no') ' Is there a "q" file too? [y|n]: '
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

      write (luncrt, '(a)', advance='no') ' Is the grid iblanked? [y|n]: '
      read  (lunkbd, *) answer

      iblanked = answer == 'y' .or. answer == 'Y'

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

!     Extract the specified blocks in the indicated ascending order:

      i = 1
      do ii = 1, nblocks_out

         ib = iblock(ii)

         if (ib > nblocks_in) then
            write (luncrt, '(a, 2i5)') &
               ' Specified block number is too big:', ib, nblocks_in
            cycle
         end if

         do while (i <= ib) ! Read an input block - no way of skipping

            if (iblanked) then
               call xyzi_allocate (xyzq_in(i), ios)
            else
               call xyz_allocate (xyzq_in(i), ios)
            end if

            if (ios /= 0) go to 999

            npts = xyzq_in(i)%ni * xyzq_in(i)%nj * xyzq_in(i)%nk

            if (iblanked) then
               call xyzi_block_io (1, lunxyz_in, formatted_in, npts, &
                                   xyzq_in(i)%x, xyzq_in(i)%y, xyzq_in(i)%z,   &
                                   xyzq_in(i)%iblank, ios)
            else
               call xyz_block_io (1, lunxyz_in, formatted_in, npts, &
                                  xyzq_in(i)%x, xyzq_in(i)%y, xyzq_in(i)%z, ios)
            end if

            if (ios /= 0) go to 999

            if (i < ib) then
               deallocate (xyzq_in(i)%x, xyzq_in(i)%y, xyzq_in(i)%z, stat=ios)
               if (ios /= 0) then
                  write (luncrt, '(a, 2i5)') &
                     ' Trouble deallocating x/y/z for input block #', i, ios
                  go to 999
               end if
               if (iblanked) deallocate (xyzq_in(i)%iblank, stat=ios)
               if (ios /= 0) then
                  write (luncrt, '(a, 2i5)') &
                     ' Trouble deallocating iblank for input block #', i, ios
                  go to 999
               end if
            end if

            if (qfile) then

               call q_allocate (xyzq_in(i), num_q, ios)
               if (ios /= 0) go to 999

               call q_block_io (1, lunq_in, formatted_in, num_q, &
                                xyzq_in(i)%mi, xyzq_in(i)%mj, xyzq_in(i)%mk, &
                                xyzq_in(i)%q, ios)
               if (ios /= 0) go to 999

               if (i < ib) then
                  deallocate (xyzq_in(i)%q, stat=ios)
                  if (ios /= 0) then
                     write (luncrt, '(a, 2i5)') &
                        ' Trouble deallocating input flow block #', i, ios
                     go to 999
                  end if
               end if

            end if

            i = i + 1

         end do ! Next input block to read

!        Transcribe the block:

         if (iblanked) then
            call xyzi_block_io (2, lunxyz_out, formatted_out, npts, &
                                xyzq_in(ib)%x, xyzq_in(ib)%y, xyzq_in(ib)%z,   &
                                xyzq_in(ib)%iblank, ios)
         else
            call xyz_block_io (2, lunxyz_out, formatted_out, npts, &
                               xyzq_in(ib)%x, xyzq_in(ib)%y, xyzq_in(ib)%z, ios)
         end if

         if (ios /= 0) go to 999
         
         deallocate (xyzq_in(ib)%x, xyzq_in(ib)%y, xyzq_in(ib)%z, stat=ios)
         if (ios /= 0) then
            write (luncrt, '(a, 2i5)') &
               ' Trouble deallocating x/y/z for input block #', ib, ios
            go to 999
         end if

         if (iblanked) then
            deallocate (xyzq_in(ib)%iblank, stat=ios)
            if (ios /= 0) then
               write (luncrt, '(a, 2i5)') &
                  ' Trouble deallocating iblank for input block #', ib, ios
               go to 999
            end if
         end if

         if (qfile) then
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
         end if

      end do

  999 continue

! *** stop ! Avoid system dependencies.

      end program extract_blocks
