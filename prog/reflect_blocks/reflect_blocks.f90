!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      program reflect_blocks
!
!     Description:
!
!        REFLECT_BLOCKS reflects the blocks of a multiblock grid in the plane
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
!        XYZQ_IO package  I/O utilities for PLOT3D grid and function files
!        RDLIST           Utility for reading an indefinite list of integers
!
!     History:
!
!        08/30/05  D.A.Saunders  Initial adaptation of EXTRACT_BLOCKS.
!        02/23/06    "     "     Belated realization that reflecting changes the
!                                handedness.  Restore the original.
!
!     Author:  David Saunders, ELORET/NASA Ames Research Cntr, Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Modules:

      use grid_block_structure
      use xyzq_io_module

      implicit none

!     Constants:

      integer, parameter :: &
         lunkbd     = 5,    &
         luncrt     = 6,    &
         lunxyz_in  = 7,    &
         lunq_in    = 8,    &
         lunxyz_out = 9,    &
         lunq_out   = 10

!     Variables:

      integer :: &
         i, ib, ios, ir, ixyz, j, k, m, n, nblocks_in, nblocks_out,            &
         nchange_sign, ni, nj, nk, npts, num_q

      integer, allocatable, dimension (:) :: &
         iq_change_sign

      real :: &
         tol

      logical :: &
         both_halves, clean_up_symmetry_plane, formatted_in, formatted_out,    &
         qfile

      character :: &
         answer * 1

!     Composite data types:

      type (grid_type), pointer, dimension (:) :: &
         xyzq_in, xyzq_out

!     Execution:

      call file_prompt (lunxyz_in, -luncrt, 'input grid', 'old', .true.,       &
                        formatted_in, ios)
      if (ios /= 0) go to 999

      write (luncrt, '(a)', advance='no') ' Is there a "q" file too? [y|n]: '
      read  (lunkbd, *) answer
      qfile = answer == 'y' .or. answer == 'Y'

      if (qfile) then
         call file_prompt (lunq_in, -luncrt, 'input flow', 'old', .false.,     &
                           formatted_in, ios)
         if (ios /= 0) go to 999
      end if

   20 write (luncrt, '(a)', advance='no') &
         ' What coordinate is to be reflected?  [1|2|3 => x|y|z]: '
      read  (lunkbd, *, iostat=ios) ixyz
      if (ios /= 0) go to 20
      if (ixyz < 1 .or. ixyz > 3) go to 20

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

      call xyz_header_io (1, lunxyz_in, formatted_in, nblocks_in, xyzq_in, ios)

      if (ios /= 0) go to 999

      if (qfile) then
         call q_header_io (1, lunq_in, formatted_in, nblocks_in, num_q,        &
                           xyzq_in, ios)
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
      write (luncrt, '(4(i6, 2x, 3i5))') &
         (ib, xyzq_in(ib)%ni, xyzq_in(ib)%nj, xyzq_in(ib)%nk, &
          ib = 1, nblocks_in)

!     Set up the output file(s):

      call file_prompt (lunxyz_out, -luncrt, 'output grid', 'unknown', .true., &
                        formatted_out, ios)
      if (ios /= 0) go to 999

      if (qfile) then
         call file_prompt (lunq_out, -luncrt, 'output flow', 'unknown',        &
                           .false., formatted_out, ios)
         if (ios /= 0) go to 999

         allocate (iq_change_sign(num_q))

         nchange_sign = num_q

         call rdlist (luncrt, 'Function number(s) needing sign change: ',      &
                      lunkbd, nchange_sign, iq_change_sign)
      end if

      allocate (xyzq_out(nblocks_out))

!     Set up the output block dimensions:

      do ib = 1, nblocks_in

         xyzq_out(ib)%ni = xyzq_in(ib)%ni
         xyzq_out(ib)%nj = xyzq_in(ib)%nj
         xyzq_out(ib)%nk = xyzq_in(ib)%nk

         if (qfile) then
            xyzq_out(ib)%mi = xyzq_in(ib)%mi
            xyzq_out(ib)%mj = xyzq_in(ib)%mj
            xyzq_out(ib)%mk = xyzq_in(ib)%mk
         end if

         if (both_halves) then

            ir = ib + nblocks_in
            xyzq_out(ir)%ni = xyzq_in(ib)%ni
            xyzq_out(ir)%nj = xyzq_in(ib)%nj
            xyzq_out(ir)%nk = xyzq_in(ib)%nk

            if (qfile) then
               xyzq_out(ir)%mi = xyzq_in(ib)%mi
               xyzq_out(ir)%mj = xyzq_in(ib)%mj
               xyzq_out(ir)%mk = xyzq_in(ib)%mk
            end if

         end if

      end do

      call xyz_header_io (2, lunxyz_out, formatted_out, nblocks_out, xyzq_out, &
                          ios)
      if (ios /= 0) go to 999

      if (qfile) then

         call q_header_io (2, lunq_out, formatted_out, nblocks_out, num_q, &
                           xyzq_out, ios)
         if (ios /= 0) go to 999
      end if

!     Read the blocks one at a time, clean them up if specified, and write
!     them as is if both halves were requested:

      do ib = 1, nblocks_in

         call xyz_allocate (xyzq_in(ib), ios)
         if (ios /= 0) go to 999

         ni = xyzq_in(ib)%ni;  nj = xyzq_in(ib)%nj;  nk = xyzq_in(ib)%nk
         npts = ni * nj * nk

         call xyz_block_io (1, lunxyz_in, formatted_in, npts, &
                            xyzq_in(ib)%x, xyzq_in(ib)%y, xyzq_in(ib)%z, ios)
         if (ios /= 0) go to 999

         if (clean_up_symmetry_plane) then

            select case (ixyz)

               case (1)

                  call clean_up (ib, ni, nj, nk, tol, xyzq_in(ib)%x)

               case (2)

                  call clean_up (ib, ni, nj, nk, tol, xyzq_in(ib)%y)

               case (3)

                  call clean_up (ib, ni, nj, nk, tol, xyzq_in(ib)%z)

            end select

         end if

         if (qfile) then

            call q_allocate (xyzq_in(ib), num_q, ios)
            if (ios /= 0) go to 999

            call q_block_io (1, lunq_in, formatted_in, num_q,                  &
                             xyzq_in(ib)%mi, xyzq_in(ib)%mj, xyzq_in(ib)%mk,   &
                             xyzq_in(ib)%q, ios)
            if (ios /= 0) go to 999

         end if

!        Transcribe the input block?

         if (both_halves) then

            call xyz_block_io (2, lunxyz_out, formatted_out, npts,             &
                               xyzq_in(ib)%x, xyzq_in(ib)%y, xyzq_in(ib)%z, ios)
            if (ios /= 0) go to 999

            if (qfile) then

               call q_block_io (2, lunq_out, formatted_out, num_q,             &
                               xyzq_in(ib)%mi, xyzq_in(ib)%mj, xyzq_in(ib)%mk, &
                               xyzq_in(ib)%q, ios)
               if (ios /= 0) go to 999

            end if

         else
!           We don't want to reorder the original blocks by interleaving them
!           with the reflected blocks, so we wait until all blocks are read.
         end if

      end do ! Next input block to read

!     Perform the reflection:

      do ib = 1, nblocks_in

         ni = xyzq_in(ib)%ni;  nj = xyzq_in(ib)%nj;  nk = xyzq_in(ib)%nk
         npts = ni * nj * nk

         select case (ixyz)

            case (1)

               call reflect_block (npts, xyzq_in(ib)%x)

            case (2)

               call reflect_block (npts, xyzq_in(ib)%y)

            case (3)

               call reflect_block (npts, xyzq_in(ib)%z)

         end select

         if (.not. qfile) then

!           Restore the original handedness:

            call fix_handedness (xyzq_in(ib), num_q)  ! num_q = 0

            call xyz_block_io (2, lunxyz_out, formatted_out, npts,             &
                            xyzq_in(ib)%x, xyzq_in(ib)%y, xyzq_in(ib)%z, ios)
            if (ios /= 0) go to 999

            deallocate (xyzq_in(ib)%x, xyzq_in(ib)%y, xyzq_in(ib)%z, stat=ios)

            if (ios /= 0) then
               write (luncrt, '(a, 2i5)') &
                  ' Trouble deallocating grid block #', ib, ios
               go to 999
            end if

         else ! We can't fix xyz handedness until we've dealt with the flow data

            if (nchange_sign > 0) then
               do k = 1, xyzq_in(ib)%mk
                  do j = 1, xyzq_in(ib)%mj
                     do i = 1, xyzq_in(ib)%mi
                        do m = 1, nchange_sign
                           n = iq_change_sign(m)
                           xyzq_in(ib)%q(n,i,j,k) = -xyzq_in(ib)%q(n,i,j,k)
                        end do
                     end do
                  end do
               end do
            end if

            call fix_handedness (xyzq_in(ib), num_q)

            call xyz_block_io (2, lunxyz_out, formatted_out, npts,             &
                            xyzq_in(ib)%x, xyzq_in(ib)%y, xyzq_in(ib)%z, ios)
            if (ios /= 0) go to 999

            deallocate (xyzq_in(ib)%x, xyzq_in(ib)%y, xyzq_in(ib)%z, stat=ios)

            if (ios /= 0) then
               write (luncrt, '(a, 2i5)') &
                  ' Trouble deallocating grid block #', ib, ios
               go to 999
            end if

            call q_block_io (2, lunq_out, formatted_out, num_q,                &
                             xyzq_in(ib)%mi, xyzq_in(ib)%mj, xyzq_in(ib)%mk,   &
                             xyzq_in(ib)%q, ios)
            if (ios /= 0) go to 999

            deallocate (xyzq_in(ib)%q, stat=ios)

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
         subroutine clean_up (ib, ni, nj, nk, tol, x)
!
!        Look for a block face that appears to be on the symmetry plane, and set
!        its given coordinates to zero.  The strategy is to find the average
!        absolute value of the indicated coordinate, and if this is below the
!        given tolerance then assume the face is meant to be on the "x" = 0
!        plane exactly.  "x" may be x, y, or z.
!
!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         integer, intent (in) :: ib          ! Block number
         integer, intent (in) :: ni, nj, nk  ! Block dimensions
         real,    intent (in) :: tol         ! Criterion for cleaning up a face
         real, intent (inout) :: x(ni,nj,nk) ! One coordinate for all points

         real, parameter :: zero = 0.

         integer :: i, j, k
         real    :: sum, tol100

         tol100 = tol * 100

         sum = zero
         do k = 1, nk ! i = 1 face
            do j = 1, nj
               sum = sum + abs (x(1,j,k))
            end do
         end do
         sum = sum / real (nj * nk)
         if (sum < tol100) write (luncrt, '(a, f12.7, a, i4)') &
            ' Face 1 has average |coordinate| of', sum, '  for block', ib
         if (sum < tol) then
            write (luncrt, '(a)') '    Forcing face 1 to symmetry plane.'
            x(1,1:nj,1:nk) = zero
         end if

         sum = zero
         do k = 1, nk ! i = ni face
            do j = 1, nj
               sum = sum + abs (x(ni,j,k))
            end do
         end do
         sum = sum / real (nj * nk)
         if (sum < tol100) write (luncrt, '(a, f12.7, a, i4)') &
            ' Face 2 has average |coordinate| of', sum, '  for block', ib
         if (sum < tol) then
            write (luncrt, '(a)') '    Forcing face 2 to symmetry plane.'
            x(ni,1:nj,1:nk) = zero
         end if

         sum = zero
         do k = 1, nk ! j = 1 face
            do i = 1, ni
               sum = sum + abs (x(i,1,k))
            end do
         end do
         sum = sum / real (nk * ni)
         if (sum < tol100) write (luncrt, '(a, f12.7, a, i4)') &
            ' Face 3 has average |coordinate| of', sum, '  for block', ib
         if (sum < tol) then
            write (luncrt, '(a)') '    Forcing face 3 to symmetry plane.'
            x(1:ni,1,1:nk) = zero
         end if

         sum = zero
         do k = 1, nk ! j = nj face
            do i = 1, ni
               sum = sum + abs (x(i,nj,k))
            end do
         end do
         sum = sum / real (nk * ni)
         if (sum < tol100) write (luncrt, '(a, f12.7, a, i4)') &
            ' Face 4 has average |coordinate| of', sum, '  for block', ib
         if (sum < tol) then
            write (luncrt, '(a)') '    Forcing face 4 to symmetry plane.'
            x(1:ni,nj,1:nk) = zero
         end if

         sum = zero
         do j = 1, nj ! k = 1 face
            do i = 1, ni
               sum = sum + abs (x(i,j,1))
            end do
         end do
         sum = sum / real (ni * nj)
         if (sum < tol100) write (luncrt, '(a, f12.7, a, i4)') &
            ' Face 5 has average |coordinate| of', sum, '  for block', ib
         if (sum < tol) then
            write (luncrt, '(a)') '    Forcing face 5 to symmetry plane.'
            x(1:ni,1:nj,1) = zero
         end if

         sum = zero
         do j = 1, nj ! k = nk face
            do i = 1, ni
               sum = sum + abs (x(i,j,nk))
            end do
         end do
         sum = sum / real (ni * nj)
         if (sum < tol100) write (luncrt, '(a, f12.7, a, i4)') &
            ' Face 6 has average |coordinate| of', sum, '  for block', ib
         if (sum < tol) then
            write (luncrt, '(a)') '    Forcing face 6 to symmetry plane.'
            x(1:ni,1:nj,nk) = zero
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

      end program reflect_blocks

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine fix_handedness (block, nvars)

!  Change the handedness of a grid block by reversing the j indices.
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

   integer                          :: i, j, k, l, njp1, njby2
   real                             :: t
   real, allocatable, dimension (:) :: s

!  Execution:

   njp1  = block%nj + 1
   njby2 = njp1 / 2

   do k = 1, block%nk
      do i = 1, block%ni
         do j = 1, njby2
            l = njp1 - j
            t = block%x(i,j,k)
            block%x(i,j,k) = block%x(i,l,k)
            block%x(i,l,k) = t
            t = block%y(i,j,k)
            block%y(i,j,k) = block%y(i,l,k)
            block%y(i,l,k) = t
            t = block%z(i,j,k)
            block%z(i,j,k) = block%z(i,l,k)
            block%z(i,l,k) = t
         end do
      end do
   end do

   if (nvars > 0) then

      allocate (s(nvars))

      do k = 1, block%nk
         do i = 1, block%ni
            do j = 1, njby2
               l = njp1 - j
               s = block%q(:,i,j,k)
               block%q(:,i,j,k) = block%q(:,i,l,k)
               block%q(:,i,l,k) = s
            end do
         end do
      end do

      deallocate (s)

   end if

   end subroutine fix_handedness
