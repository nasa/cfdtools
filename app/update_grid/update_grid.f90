!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program update_grid
!
!  Description:
!
!     UPDATE_GRID overwrites the indicated block(s) of a multiblock grid.
!  Blocks being replaced should be specified in ascending order so that
!  only one block needs to be stored at a time.  The dimensions of the
!  block(s) being replaced cannot change.
!
!     UPDATE_GRID might be used in conjunction with GSMOOTH, which smooths
!  one block of a grid at a time.
!
!  Known weakness:
!
!     If more than one block is being updated from the same replacement file,
!  there is no attempt to suppress the prompt for that file (it is done from
!  within the I/O package) or to avoid repeated reads from the top of the
!  file to the desired blocks.  The anticipated usage is with multiple
!  single-block files being plugged into one multiblock file.  Substituting
!  patches in a surface grid is another possibility.
!
!  History:
!
!     02/20/00  D.A.Saunders  Initial adaptation of GSMOOTH.
!     10/20/03    "      "    Remind the user to enter block changes in order.
!     02/10/14    "      "    Replaced CFD_IO_PACKAGE with xyzq_io package.
!     04/12/16    "      "    Fixed a glitch in transfer of untouched blocks.
!
!  Author:  David Saunders, ERC, Inc./NASA Ames Research Center, CA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_block_structure   ! Derived data type for one grid block
   use xyzq_io_module         ! PLOT3D file I/O package

   implicit none

!  Constants:

   integer, parameter :: &
      lunold = 1,        &    ! Grid being updated
      lunnew = 2,        &    ! Grid with desired block(s)
      lunout = 3,        &    ! Updated output grid
      lunkbd = 5,        &
      luncrt = 6

!  Variables:

   integer :: &
      ib, ios, n, nblock, nblocknew, ni, nj, nk, nnew, nprevious, npts

   logical :: &
      cell_centered, formatted_old, formatted_new, formatted_out

!  Composite data types:

   type (grid_type), pointer, dimension (:) :: &
      oldgrid, newgrid

!  Execution:

   write (luncrt, '(/, (a))') &
      '   N.B.:  Update blocks in ascending order, and', &
      '          don''t try to change any block dimensions.'

!  Get the initial grid name and formatting:

   call file_prompt (lunold, -luncrt, 'initial grid', 'old', .true., &
                     formatted_old, ios)

!  Read the initial grid header:

   call xyz_header_io (1, lunold, formatted_old, nblock, oldgrid, ios)
   if (ios /= 0) go to 99

!  Get the output grid name and formatting:

   call file_prompt (lunout, -luncrt, 'output grid', 'unknown', .true., &
                     formatted_out, ios)
   if (ios /= 0) go to 99

!  The output grid header info is the same as the initial header:

   call xyz_header_io (2, lunout, formatted_out, nblock, oldgrid, ios)
   if (ios /= 0) go to 99

!  Display the block dimensions to help specification of substitutes:

   write (luncrt, '(/, (i6, 2x, 3i5))') &
      (n, oldgrid(n)%ni, oldgrid(n)%nj, oldgrid(n)%nk, n = 1, nblock)

   nprevious = 0

!  Look for blocks to update, in ascending order:

   do ! Until no more blocks to update

      write (luncrt, '(/, a)', advance='no') &
         ' Block number to replace [n > 0; -99 = no more]: '
      read  (lunkbd, *) n

      if (n > nprevious .and. n <= nblock) then

         ni = oldgrid(n)%ni
         nj = oldgrid(n)%nj
         nk = oldgrid(n)%nk

!        Transcribe blocks from input to output, up to the specified block:

         do ib = nprevious + 1, n

            write (luncrt, '(a, i4)') 'Reading old block', ib

            call xyz_allocate (oldgrid(ib), ios)

            npts = oldgrid(ib)%ni * oldgrid(ib)%nj * oldgrid(ib)%nk

            call xyz_block_io (1, lunold, formatted_old, npts, &
                               oldgrid(ib)%x, oldgrid(ib)%y, oldgrid(ib)%z, ios)

            if (ib < n) then  ! Leave output file ready for new block n

               write (luncrt, '(a, i4)') 'Writing old block', ib

               call xyz_block_io (2, lunout, formatted_out, npts, &
                               oldgrid(ib)%x, oldgrid(ib)%y, oldgrid(ib)%z, ios)
               deallocate     (oldgrid(ib)%x, oldgrid(ib)%y, oldgrid(ib)%z)

            end if

         end do

         nprevious = n

         call file_prompt (lunnew, -luncrt, 'replacement block grid', &
                           'unknown', .true., formatted_new, ios)
         if (ios /= 0) go to 99

!        Read the header of the file with the replacement for block n:

         write (luncrt, '(a)') 'Reading header of new block file.'

         call xyz_header_io (1, lunnew, formatted_new, nblocknew, newgrid, ios)
         if (ios /= 0) go to 99

         if (nblocknew == 1) then
            nnew = 1
         else
            nnew = 0
            do while (nnew < 1 .or. nnew > nblocknew)
               write (luncrt, '(a, i4, a)', advance='no') &
                  ' Replace block', n, '  with which block?: '
               read  (lunkbd, *) nnew
            end do
         end if

         if (newgrid(nnew)%ni /= ni .or. &
             newgrid(nnew)%nj /= nj .or. &
             newgrid(nnew)%nk /= nk) then

            write (luncrt, '(/, a, /, (a, 3i6))') &
               ' Mismatched block dimensions.',   &
               ' Old: ', ni, nj, nk,              &
               ' New: ', newgrid(nnew)%ni, newgrid(nnew)%nj, newgrid(nnew)%nk, &
               ' Proceeding.'
            cycle
         end if

!        Read the replacement for block n (no direct access):

         do ib = 1, nnew
            write (luncrt, '(a, i4)') 'Reading new block', ib

            call xyz_allocate (newgrid(ib), ios)

            npts = newgrid(ib)%ni * newgrid(ib)%nj * newgrid(ib)%nk

            call xyz_block_io (1, lunnew, formatted_new, npts, &
                               newgrid(ib)%x, newgrid(ib)%y, newgrid(ib)%z, ios)
            if (ib < nnew) then
               deallocate     (newgrid(ib)%x, newgrid(ib)%y, newgrid(ib)%z)
            end if
         end do

!        Write the replacement for block n:

         ib = nnew
         write (luncrt, '(a, i4)') 'Writing new block', ib

         call xyz_block_io (2, lunout, formatted_out, npts, &
                            newgrid(ib)%x, newgrid(ib)%y, newgrid(ib)%z, ios)
         deallocate        (newgrid(ib)%x, newgrid(ib)%y, newgrid(ib)%z)
         close (lunnew)

      else if (n <= 0) then ! -99

!        Transfer any remaining blocks from input to output:

         do ib = nprevious + 1, nblock

            call xyz_allocate (oldgrid(ib), ios)

            npts = oldgrid(ib)%ni * oldgrid(ib)%nj * oldgrid(ib)%nk

            write (luncrt, '(a, i4)') 'Reading old block', ib

            call xyz_block_io (1, lunold, formatted_old, npts, &
                               oldgrid(ib)%x, oldgrid(ib)%y, oldgrid(ib)%z, ios)

            write (luncrt, '(a, i4)') 'Writing old block', ib

            call xyz_block_io (2, lunout, formatted_out, npts, &
                               oldgrid(ib)%x, oldgrid(ib)%y, oldgrid(ib)%z, ios)
            deallocate        (oldgrid(ib)%x, oldgrid(ib)%y, oldgrid(ib)%z)

         end do

         close (lunout)

         exit

      else ! n <= nprevious .or. n > nblock

         write (luncrt, '(/, (3x, a))') &
            'Illegal block number.',    &
            'Not in ascending order?  Greater than nblock?', &
            'If necessary, save results (-99) and perform further updates', &
            'in another run.'

      end if

   end do ! Next replacement block

99 continue

   end program update_grid
