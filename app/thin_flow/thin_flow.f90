!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      program thin_flow
!
!     Description:
!
!        THIN_FLOW extracts a subset of points from all blocks of a multiblock
!     flow solution.
!
!     Procedures:
!
!        XYZQ_IO package  I/O utilities for PLOT3D grid and function files
!
!     History:
!
!        12/30/04  D. Saunders  Adaptation of THIN_GRID (easier than combining
!                               grid and/or flow cases).
!
!     Author:  David Saunders, ELORET/NASA Ames, Moffett Field, CA.
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
         i, ib, ii, inc, ios, j, jj, jnc, k, kk, knc, nblocks, nf, mi, mj, mk

      logical :: &
         formatted_in, formatted_out

!     Composite data types:

      type (grid_type), pointer, dimension (:) :: &
         gridin, gridout

!     Execution:
!     ----------

      call file_prompt (lunin,  -luncrt, 'input flow solution', 'old', true, &
                        formatted_in, ios)
      if (ios /= 0) go to 999

!     Read the number of blocks, allocate them, and read dimensions:

      if (formatted_in) then
         read (lunin, *, iostat=ios) nblocks
      else
         read (lunin,    iostat=ios) nblocks
      end if

      if (ios /= 0) then
         write (luncrt, '(/, a)') ' Trouble reading the number of blocks.'
         go to 999
      end if

      rewind (lunin)

      allocate (gridin(nblocks))

      call q_header_io (1, lunin, formatted_in, nblocks, nf, gridin, ios)

      if (ios /= 0) then
         write (luncrt, '(/, a)') ' Trouble reading function file header.'
         go to 999
      end if

      write (luncrt, '(a)', advance='no') &
         ' Increments to use for thinning i, j and k: '

      read  (lunkbd, *) inc, jnc, knc

      call file_prompt (lunout, -luncrt, 'output flow file', 'unknown', true, &
                        formatted_out, ios)
      if (ios /= 0) go to 999

      allocate (gridout(nblocks))

      do ib = 1, nblocks
         gridout(ib)%mi = (gridin(ib)%mi + inc - 1) / inc
         gridout(ib)%mj = (gridin(ib)%mj + jnc - 1) / jnc
         gridout(ib)%mk = (gridin(ib)%mk + knc - 1) / knc
      end do

      call q_header_io (2, lunout, formatted_out, nblocks, nf, gridout, ios)

      if (ios /= 0) go to 999

!     Thin the blocks one at a time:

      do ib = 1, nblocks

         call q_allocate (gridin(ib), nf, ios)
         if (ios /= 0) go to 999

         mi = gridin(ib)%mi;  mj = gridin(ib)%mj;  mk = gridin(ib)%mk

         call q_block_io (1, lunin, formatted_in, nf, mi, mj, mk,             &
                          gridin(ib)%q, ios)
         if (ios /= 0) go to 999

         call q_allocate (gridout(ib), nf, ios)
         if (ios /= 0) go to 999

         kk = 0
         do k = 1, mk, knc
            kk = kk + 1
            jj = 0
            do j = 1, mj, jnc
               jj = jj + 1
               ii = 0
               do i = 1, mi, inc
                  ii = ii + 1
                  gridout(ib)%q(:,ii,jj,kk) = gridin(ib)%q(:,i,j,k)
               end do
            end do
         end do

         call q_block_io (2, lunout, formatted_out, nf, &
                          gridout(ib)%mi, gridout(ib)%mj, gridout(ib)%mk,     &
                          gridout(ib)%q, ios)
         if (ios /= 0) go to 999
         
         deallocate (gridin(ib)%q, gridout(ib)%q)

      end do

  999 continue

! *** stop ! Avoid system dependencies.

      end program thin_flow
