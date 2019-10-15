!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      program extract_blocks_2d
!
!     Description:
!
!        This is a 2D analogue of EXTRACT_BLOCKS.
!
!        EXTRACT_BLOCKS_2D extracts a list of blocks from a multiblock grid and
!     an optional "q" file (or rather a PLOT3D-type function file).
!
!        This version can handle "iblanked" files at the cost of another prompt.
!
!     Control file:
!
!        This is dispensed with.  At the prompt, enter a list of desired blocks
!     in ascending order, on a single line.  E.g.:
!
!        1 3 7 : : 15
!
!     Procedures:
!
!        RDLIST          Read an indefinite list of integers
!        XYQ_IO package  I/O utilities for PLOT3D grid and function files
!
!     History:
!
!        04/06/04  D.A.Saunders  EXTRACT_BLOCKS adaptation of GRID_FACES.
!        09/14/18    "      "    EXTRACT_BLOCKS_2D variant for the case of
!                                extracting the forebody block 1 from the
!                                2-block output from COMPRESS2D.  Including an
!                                aft body in the CAPSULE_GRID step and its
!                                hyperbolic volume grid step is a workaround
!                                for the splay problem encountered with Gridgen
!                                on forebody-only cases.
!
!     Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA
!              Later with AMA, Inc. at NASA ARC.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Modules:

      use grid_block_structure
      use xyq_io_module

      implicit none

!     Constants:

      integer, parameter :: &
         lunkbd    = 5,     &
         luncrt    = 6,     &
         lunxy_in  = 7,     &
         lunq_in   = 8,     &
         lunxy_out = 9,     &
         lunq_out  = 10,    &
         mxextract = 32

!     Variables:

      integer :: &
         i, ib, ii, ios, lunqi, lunqo, nblocks_in, nblocks_out, npts, num_q

      integer, dimension (mxextract) :: &
         iblock

      logical :: &
         cell_centered, formatted_in, formatted_out, iblanked, qfile

      character (1) :: &
         answer

!     Derived data types:

      type (grid_type), pointer, dimension (:) :: &
         xyq_in, xyq_out

!     Execution:

      call file_prompt_2d (lunxy_in, 'input grid', 'old', .true., &
                           formatted_in, ios)
      if (ios /= 0) go to 999

      write (luncrt, '(a)', advance='no') ' Is there a "q" file too? [y|n]: '
      read  (lunkbd, *) answer

      qfile = answer == 'y' .or. answer == 'Y'

      if (qfile) then
         call file_prompt_2d (lunq_in, 'input flow', 'old', .false., &
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

      call xy_header_io (1, lunxy_in, formatted_in, nblocks_in, xyq_in, ios)
      if (ios /= 0) go to 999

      if (qfile) then
         call q_header_io_2d (1, lunqi, formatted_in, nblocks_in, num_q, &
                              xyq_in, ios)
         if (ios /= 0) go to 999
      end if

!     Display the block dimensions to reassure the user:

      write (luncrt, '(a)')
      write (luncrt, '(4(i6, 2X, 2i5))') &
         (ib, xyq_in(ib)%ni, xyq_in(ib)%nj, ib = 1, nblocks_in)

!     Prompt for the blocks to be extracted (in order):

      nblocks_out = mxextract
      call rdlist (luncrt, 'Blocks to extract (ascending order): ', &
                   lunkbd, nblocks_out, iblock)
      if (nblocks_out <= 0) go to 999

!     Set up the output file(s):

      call file_prompt_2d (lunxy_out, 'output grid', 'unknown', &
                           .true., formatted_out, ios)
      if (ios /= 0) go to 999

      if (qfile) then
         call file_prompt_2d (lunq_out, 'output flow', 'unknown', &
                              .false., formatted_out, ios)
         if (ios /= 0) go to 999
      end if

      allocate (xyq_out(nblocks_out))

      do i = 1, nblocks_out

         ib = iblock(i)
         xyq_out(i)%ni = xyq_in(ib)%ni
         xyq_out(i)%nj = xyq_in(ib)%nj
         xyq_out(i)%nk = xyq_in(ib)%nk

         if (qfile) then
            xyq_out(i)%mi = xyq_in(ib)%mi
            xyq_out(i)%mj = xyq_in(ib)%mj
            xyq_out(i)%mk = xyq_in(ib)%mk
         end if

      end do

      call xy_header_io (2, lunxy_out, formatted_out, nblocks_out, &
                         xyq_out, ios)
      if (ios /= 0) go to 999

      if (qfile) then

         call q_header_io_2d (2, lunq_out, formatted_out, nblocks_out, num_q, &
                              xyq_out, ios)
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
               call xyi_allocate (xyq_in(i), ios)
            else
               call xy_allocate (xyq_in(i), ios)
            end if

            if (ios /= 0) go to 999

            npts = xyq_in(i)%ni * xyq_in(i)%nj

            if (iblanked) then
               call xyi_block_io (1, lunxy_in, formatted_in, npts, &
                                  xyq_in(i)%x, xyq_in(i)%y, xyq_in(i)%iblank, &
                                  ios)
            else
               call xy_block_io (1, lunxy_in, formatted_in, npts, &
                                 xyq_in(i)%x, xyq_in(i)%y, ios)
            end if

            if (ios /= 0) go to 999

            if (i < ib) then
               deallocate (xyq_in(i)%x, xyq_in(i)%y, stat=ios)
               if (ios /= 0) then
                  write (luncrt, '(a, 2i5)') &
                     ' Trouble deallocating x/y for input block #', i, ios
                  go to 999
               end if
               if (iblanked) deallocate (xyq_in(i)%iblank, stat=ios)
               if (ios /= 0) then
                  write (luncrt, '(a, 2i5)') &
                     ' Trouble deallocating iblank for input block #', i, ios
                  go to 999
               end if
            end if

            if (qfile) then

               call q_allocate_2d (xyq_in(i), num_q, ios)
               if (ios /= 0) go to 999

               call q_block_io_2d (1, lunq_in, formatted_in, num_q, &
                                   xyq_in(i)%mi, xyq_in(i)%mj, xyq_in(i)%q, ios)
               if (ios /= 0) go to 999

               if (i < ib) then
                  deallocate (xyq_in(i)%q, stat=ios)
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
            call xyi_block_io (2, lunxy_out, formatted_out, npts, &
                               xyq_in(ib)%x, xyq_in(ib)%y, xyq_in(ib)%iblank, &
                               ios)
         else
            call xy_block_io (2, lunxy_out, formatted_out, npts, &
                              xyq_in(ib)%x, xyq_in(ib)%y, ios)
         end if

         if (ios /= 0) go to 999
         
         deallocate (xyq_in(ib)%x, xyq_in(ib)%y, stat=ios)
         if (ios /= 0) then
            write (luncrt, '(a, 2i5)') &
               ' Trouble deallocating x/y for input block #', ib, ios
            go to 999
         end if

         if (iblanked) then
            deallocate (xyq_in(ib)%iblank, stat=ios)
            if (ios /= 0) then
               write (luncrt, '(a, 2i5)') &
                  ' Trouble deallocating iblank for input block #', ib, ios
               go to 999
            end if
         end if

         if (qfile) then
            call q_block_io_2d (2, lunq_out, formatted_out, num_q, &
                                xyq_in(ib)%mi, xyq_in(ib)%mj, xyq_in(ib)%q, &
                                ios)
            if (ios /= 0) go to 999

            deallocate (xyq_in(ib)%q, stat=ios)
            if (ios /= 0) then
               write (luncrt, '(a, 2i5)') &
                  ' Trouble deallocating input flow block #', ib, ios
               go to 999
            end if
         end if

      end do

  999 continue

! *** stop ! Avoid system dependencies.

      end program extract_blocks_2d
