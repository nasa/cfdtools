!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      program reflect_blocks_2d
!
!     Description:
!
!        This is the 2D analogue of REFLECT_BLOCKS, from which it is adapted.
!
!        REFLECT_BLOCKS_2D reflects the blocks of a multiblock grid in the plane
!     specified.  An optional function file may also be reflected.  Either both
!     halves or just the reflected half may be output.  An option to force
!     exact zeros at the symmetry plane is provided.  If a function file is
!     present, one or more quantities may have the sign changed as part of the
!     reflection.  Normally this would be just one velocity component.
!
!        The function file may represent either a vertex- or cell-centered
!     solution, with or without halo cells.  (Only the number of blocks is
!     required to match that of the grid.)
!
!     Procedures:
!
!        XYQ_IO package  I/O utilities for PLOT3D grid and function files
!        RDLIST          Utility for reading an indefinite list of integers
!
!     History:
!
!        02/23/06  D.A.Saunders  REFLECT_BLOCKS from which this is adapted.
!        03/25/15    "     "     REFLECT_BLOCKS_2D is prompted by a Comet Sample
!                                Return capsule for which some body-normal lines
!                                of sight cross the X-axis, causing trouble with
!                                line/surface intersections unless the whole
!                                volume grid is present, not just the top half.
!        04/12/15    "     "     Preserving handedness by reordering j caused
!                                another problem: we want j = 1 at the surface.
!                                Therefore, reverse the i indices instead.
!
!     Author:  David Saunders, ELORET/NASA Ames Research Cntr, Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Modules:

      use grid_block_structure
      use xyq_io_module

      implicit none

!     Constants:

      integer, parameter :: &
         lunkbd    = 5,    &
         luncrt    = 6,    &
         lunxy_in  = 7,    &
         lunq_in   = 8,    &
         lunxy_out = 9,    &
         lunq_out  = 10

!     Variables:

      integer :: &
         i, ib, ios, ir, ixy, j, k, m, n, nblocks_in, nblocks_out, &
         nchange_sign, ni, nj, npts, num_q

      integer, allocatable, dimension (:) :: &
         iq_change_sign

      real :: &
         tol

      logical :: &
         both_halves, clean_up_symmetry_plane, formatted_in, formatted_out, &
         qfile

      character :: &
         answer * 1

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
      end if

   20 write (luncrt, '(a)', advance='no') &
         ' What coordinate is to be reflected?  [1|2 => x|y]: '
      read  (lunkbd, *, iostat=ios) ixy
      if (ios /= 0) go to 20
      if (ixy < 1 .or. ixy > 2) go to 20

      write (luncrt, '(a)', advance='no') &
         ' Save both halves [y]?  [n = just save reflected half]: '
      read  (lunkbd, *) answer
      both_halves = answer == 'y' .or. answer == 'Y'

      write (luncrt, '(a)', advance='no') &
         ' Clean up the symmetry plane? [y|n]: '
      read  (lunkbd, *) answer
      clean_up_symmetry_plane = answer == 'y' .or. answer == 'Y'

      if (clean_up_symmetry_plane) then
         write (luncrt, '(a)', advance='no') &
            ' Tolerance for average absolute value of appropriate coordinate: '
         read  (lunkbd, *) tol
      end if

!     Read the input file header(s):

      call xy_header_io (1, lunxy_in, formatted_in, nblocks_in, xyq_in, ios)

      if (ios /= 0) go to 999

      if (qfile) then
         call q_header_io_2d (1, lunq_in, formatted_in, nblocks_in, num_q, &
                              xyq_in, ios)
         if (ios /= 0) go to 999
      else
         num_q = 0
      end if

      if (both_halves) then
         nblocks_out = nblocks_in + nblocks_in
      else
         nblocks_out = nblocks_in
      end if

!     Display the block dimensions to reassure the user:

      write (luncrt, '(a)')
      write (luncrt, '(4(i6, 2x, 2i5))') &
         (ib, xyq_in(ib)%ni, xyq_in(ib)%nj, ib = 1, nblocks_in)

!     Set up the output file(s):

      call file_prompt_2d (lunxy_out, 'output grid', 'unknown', .true., &
                           formatted_out, ios)
      if (ios /= 0) go to 999

      if (qfile) then
         call file_prompt_2d (lunq_out, 'output flow', 'unknown',        &
                              .false., formatted_out, ios)
         if (ios /= 0) go to 999

         allocate (iq_change_sign(num_q))

         nchange_sign = num_q

         call rdlist (luncrt, 'Function number(s) needing sign change: ',      &
                      lunkbd, nchange_sign, iq_change_sign)
      end if

      allocate (xyq_out(nblocks_out))

!     Set up the output block dimensions:

      do ib = 1, nblocks_in

         xyq_out(ib)%ni = xyq_in(ib)%ni
         xyq_out(ib)%nj = xyq_in(ib)%nj

         if (qfile) then
            xyq_out(ib)%mi = xyq_in(ib)%mi
            xyq_out(ib)%mj = xyq_in(ib)%mj
         end if

         if (both_halves) then

            ir = ib + nblocks_in
            xyq_out(ir)%ni = xyq_in(ib)%ni
            xyq_out(ir)%nj = xyq_in(ib)%nj

            if (qfile) then
               xyq_out(ir)%mi = xyq_in(ib)%mi
               xyq_out(ir)%mj = xyq_in(ib)%mj
            end if

         end if

      end do

      call xy_header_io (2, lunxy_out, formatted_out, nblocks_out, xyq_out, &
                         ios)
      if (ios /= 0) go to 999

      if (qfile) then

         call q_header_io_2d (2, lunq_out, formatted_out, nblocks_out, num_q, &
                              xyq_out, ios)
         if (ios /= 0) go to 999
      end if

!     Read the blocks one at a time, clean them up if specified, and write
!     them as is if both halves were requested:

      do ib = 1, nblocks_in

         call xy_allocate (xyq_in(ib), ios)
         if (ios /= 0) go to 999

         ni = xyq_in(ib)%ni;  nj = xyq_in(ib)%nj
         npts = ni * nj

         call xy_block_io (1, lunxy_in, formatted_in, npts, &
                           xyq_in(ib)%x, xyq_in(ib)%y, ios)
         if (ios /= 0) go to 999

         if (clean_up_symmetry_plane) then

            select case (ixy)

               case (1)

                  call clean_up (ib, ni, nj, tol, xyq_in(ib)%x)

               case (2)

                  call clean_up (ib, ni, nj, tol, xyq_in(ib)%y)

            end select

         end if

         if (qfile) then

            call q_allocate_2d (xyq_in(ib), num_q, ios)
            if (ios /= 0) go to 999

            call q_block_io_2d (1, lunq_in, formatted_in, num_q, &
                                xyq_in(ib)%mi, xyq_in(ib)%mj, xyq_in(ib)%q, ios)
            if (ios /= 0) go to 999

         end if

!        Transcribe the input block?

         if (both_halves) then

            call xy_block_io (2, lunxy_out, formatted_out, npts, &
                              xyq_in(ib)%x, xyq_in(ib)%y, ios)
            if (ios /= 0) go to 999

            if (qfile) then

               call q_block_io_2d (2, lunq_out, formatted_out, num_q, &
                                   xyq_in(ib)%mi, xyq_in(ib)%mj, xyq_in(ib)%q, &
                                   ios)
               if (ios /= 0) go to 999

            end if

         else
!           We don't want to reorder the original blocks by interleaving them
!           with the reflected blocks, so we wait until all blocks are read.
         end if

      end do ! Next input block to read

!     Perform the reflection:

      do ib = 1, nblocks_in

         ni = xyq_in(ib)%ni;  nj = xyq_in(ib)%nj
         npts = ni * nj

         select case (ixy)

            case (1)

               call reflect_block (npts, xyq_in(ib)%x)

            case (2)

               call reflect_block (npts, xyq_in(ib)%y)

         end select

         if (.not. qfile) then

!           Restore the original handedness:

            call fix_handedness (xyq_in(ib), num_q)  ! num_q = 0

            call xy_block_io (2, lunxy_out, formatted_out, npts, &
                              xyq_in(ib)%x, xyq_in(ib)%y, ios)
            if (ios /= 0) go to 999

            deallocate (xyq_in(ib)%x, xyq_in(ib)%y, stat=ios)

            if (ios /= 0) then
               write (luncrt, '(a, 2i5)') &
                  ' Trouble deallocating grid block #', ib, ios
               go to 999
            end if

         else ! We can't fix xyz handedness until we've dealt with the flow data

            if (nchange_sign > 0) then
               do j = 1, xyq_in(ib)%mj
                  do i = 1, xyq_in(ib)%mi
                     do m = 1, nchange_sign
                        n = iq_change_sign(m)
                        xyq_in(ib)%q(n,i,j,1) = -xyq_in(ib)%q(n,i,j,1)
                     end do
                  end do
               end do
            end if

            call fix_handedness (xyq_in(ib), num_q)

            call xy_block_io (2, lunxy_out, formatted_out, npts, &
                              xyq_in(ib)%x, xyq_in(ib)%y, ios)
            if (ios /= 0) go to 999

            deallocate (xyq_in(ib)%x, xyq_in(ib)%y, stat=ios)

            if (ios /= 0) then
               write (luncrt, '(a, 2i5)') &
                  ' Trouble deallocating grid block #', ib, ios
               go to 999
            end if

            call q_block_io_2d (2, lunq_out, formatted_out, num_q, &
                                xyq_in(ib)%mi, xyq_in(ib)%mj, xyq_in(ib)%q, ios)
            if (ios /= 0) go to 999

            deallocate (xyq_in(ib)%q, stat=ios)

            if (ios /= 0) then
               write (luncrt, '(a, 2i5)') &
                  ' Trouble deallocating flow block #', ib, ios
               go to 999
            end if

         end if

      end do ! Next reflected block to write

  999 continue

! *** stop ! Avoid system dependencies.

!     Internal procedures for program reflect_blocks:

      contains

!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
         subroutine clean_up (ib, ni, nj, tol, x)
!
!        Look for a block face that appears to be on the symmetry plane, and set
!        its given coordinates to zero.  The strategy is to find the average
!        absolute value of the indicated coordinate, and if this is below the
!        given tolerance then assume the face is meant to be on the "x" = 0
!        plane exactly.  "x" may be x or y.
!
!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         integer, intent (in) :: ib          ! Block number
         integer, intent (in) :: ni, nj      ! Block dimensions
         real,    intent (in) :: tol         ! Criterion for cleaning up a face
         real, intent (inout) :: x(ni,nj,1)  ! One coordinate for all points

         real, parameter :: zero = 0.

         integer :: i, j
         real    :: sum, tol100

         tol100 = tol * 100

         sum = zero
         do j = 1, nj
            sum = sum + abs (x(1,j,1))
         end do
         sum = sum / real (nj)
         if (sum < tol100) write (luncrt, '(a, f12.7, a, i4)') &
            ' Face 1 has average |coordinate| of', sum, '  for block', ib
         if (sum < tol) then
            write (luncrt, '(a)') '    Forcing face 1 to symmetry plane.'
            x(1,1:nj,1) = zero
         end if

         sum = zero
         do j = 1, nj
            sum = sum + abs (x(ni,j,1))
         end do
         sum = sum / real (nj)
         if (sum < tol100) write (luncrt, '(a, f12.7, a, i4)') &
            ' Face 2 has average |coordinate| of', sum, '  for block', ib
         if (sum < tol) then
            write (luncrt, '(a)') '    Forcing face 2 to symmetry plane.'
            x(ni,1:nj,1) = zero
         end if

         sum = zero
         do i = 1, ni
            sum = sum + abs (x(i,1,1))
         end do
         sum = sum / real (ni)
         if (sum < tol100) write (luncrt, '(a, f12.7, a, i4)') &
            ' Face 3 has average |coordinate| of', sum, '  for block', ib
         if (sum < tol) then
            write (luncrt, '(a)') '    Forcing face 3 to symmetry plane.'
            x(1:ni,1,1) = zero
         end if

         sum = zero
         do i = 1, ni
            sum = sum + abs (x(i,nj,1))
         end do
         sum = sum / real (ni)
         if (sum < tol100) write (luncrt, '(a, f12.7, a, i4)') &
            ' Face 4 has average |coordinate| of', sum, '  for block', ib
         if (sum < tol) then
            write (luncrt, '(a)') '    Forcing face 4 to symmetry plane.'
            x(1:ni,nj,1) = zero
         end if

         end subroutine clean_up

!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
         subroutine reflect_block (npts, x)

!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         integer, intent (in) :: npts
         real, intent (inout) :: x(npts) ! One coordinate for all points

         integer :: i

         do i = 1, npts
            x(i) = -x(i)
         end do

         end subroutine reflect_block

      end program reflect_blocks_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine fix_handedness (block, nvars)

!  Change the handedness of a grid block by reversing the i indices.
!  The block is updated in place so that subsequent I/O can remain efficient.
!
!  10/15/04   David Saunders, ELORET/NASA Ames   Initial implementation.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use grid_block_structure  ! External module

   implicit none

!  Arguments:

   type (grid_type), intent (inout) :: block
   integer,          intent (in)    :: nvars ! nvars > 0 means flip "q" too

!  Local variables:

   integer                          :: i, j, l, nip1, niby2
   real                             :: t
   real, allocatable, dimension (:) :: s

!  Execution:

   nip1  = block%ni + 1
   niby2 = nip1 / 2

   do j = 1, block%nj
      do i = 1, niby2
         l = nip1 - i
         t = block%x(i,j,1)
         block%x(i,j,1) = block%x(l,j,1)
         block%x(l,j,1) = t
         t = block%y(i,j,1)
         block%y(i,j,1) = block%y(l,j,1)
         block%y(l,j,1) = t
      end do
   end do

   if (nvars > 0) then

      allocate (s(nvars))

      do j = 1, block%nj
         do i = 1, niby2
            l = nip1 - i
            s = block%q(:,i,j,1)
            block%q(:,i,j,1) = block%q(:,l,j,1)
            block%q(:,l,j,1) = s
         end do
      end do

      deallocate (s)

   end if

   end subroutine fix_handedness
