!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      program rectify_grid
!
!     Description:
!
!        RECTIFY_GRID scans a grid for blocks that are not right-handed, and
!     swaps their (i,j) ordering for x, y, and z.  An optional "q" file (or
!     rather a PLOT3D-type function file) may be similarly reordered.
!
!        For the surface case, the user needs to indicate which patches are to
!     be rectified via any meaningful integer list entered at the keyboard on a
!     single line.
!
!        E.g.:  1 3:6 11-13
!
!        In this case, if an empty line is entered, the program will write a
!     function file containing components of the unit normals to each point.
!     Typically, contour plotting of one of these components should suffice to
!     tell whether the patch is right-handed (normal pointing out) or not.
!     Omitting any function file if normals are to be calculated is sensible.
!     Then, if necessary, the program can be rerun and the appropriate patches
!     to rectify can be entered interactively.
!
!     Control file ('rectify_grid.inp'):
!
!        RECTIFY_GRID control file
!        mygrid.g     Input grid
!        T            Formatted?
!        none         Accompanying function file, or 'none'
!        T
!        rectified.g  Rectified grid
!        T            Formatted?
!        none         Rectified function file
!
!     Procedures:
!
!        XYZQ_IO package  I/O utilities for PLOT3D grid and function files
!
!     History:
!
!        04/21/04  D. Saunders  Adaptation of EXTRACT_BLOCKS.
!        04/22/04       "       Cell volumes work except for surface patches,
!                               so prompt for the patches to flip.
!        07/15/06       "       Added output of surface normal components as an
!                               aid to identifying left-handed patches.
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
         lunctl = 1,        &
         luncrt = 6,        &
         lunkbd = 5,        &
         lunxyz = 2,        &
         lunq   = 3

      real, parameter ::    &
         one = 1., zero = 0.

      character, parameter :: &
         format * 11 = 'unformatted'

!     Variables:

      integer :: &
         i, i1, ib, ic, ios, j, jc, lunfun, nblocks, nflip, ni, nj, nk,        &
         nrectify, num_q

      integer, allocatable, dimension (:) :: &
         iblock

      real :: &
         eps, p, q

      real :: &
         un(3), &
         volume(3:3,3:3,3:3) ! For use on cells defined by the (2,2,2) corner

      logical :: &
         cell_centered, formatted, normals, qfile, surface

      logical, allocatable, dimension (:) :: &
         flip

      character :: &
         filename * 80

!     Derived data types:

      type (grid_type), pointer, dimension (:) :: &
         grid

!     Execution:

      open (lunctl, file='rectify_grid.inp', status='old', iostat=ios)

      if (ios /= 0) then
         write (luncrt, '(/, a)') &
            ' Unable to open rectify_grid.inp control file.'
         go to 999
      end if

      read (lunctl, *) ! Header
      read (lunctl, *) filename
      read (lunctl, *) formatted

      i1 = 1;  if (formatted) i1 = 3

      open (lunxyz, file=filename, form=format(i1:11), status='old', iostat=ios)

      if (ios /= 0) then
         write (luncrt, '(/, 2a)') &
            ' Unable to open input grid: ', filename(1:len_trim(filename))
         go to 999
      end if

      read (lunctl, *) filename
      read (lunctl, *) ! formatted  ! The I/O package assumes grid & q match
      qfile = filename(1:4) /= 'none'

      if (qfile) then
         lunfun = lunq
         open (lunq, file=filename, form=format(i1:11), status='old',iostat=ios)
         if (ios /= 0) then
            write (luncrt, '(/, 2a)') &
               ' Unable to open input function file: ', &
               filename(1:len_trim(filename))
            go to 999
         end if
      else
         lunfun = -lunq
         num_q  = 0
      end if

!     Allocate grid arrays and read them:

      call xyzq_read (lunxyz, lunfun, formatted, nblocks, num_q, &
                      cell_centered, grid, ios)
      if (ios /= 0) go to 999

      surface = grid(1)%nk == 1

      allocate (flip(nblocks))

      if (surface) then ! We need help from the user

         allocate (iblock(nblocks))

         nflip = nblocks

         call rdlist (luncrt, ' Enter a list of patch numbers to flip:', &
                      lunkbd, nflip, iblock)

         flip(:) = .false.

         normals = nflip == 0 ! No patches entered - user needs help

         do i = 1, nflip
            ib = iblock(i)
            flip(ib) = .true.
         end do

         if (normals) then

            if (qfile) then ! Deallocate and reuse its fields

               do ib = 1, nblocks
                  deallocate (grid(ib)%q)
               end do

            end if

            num_q = 3

            do ib = 1, nblocks

               ni = grid(ib)%ni;  grid(ib)%mi = ni
               nj = grid(ib)%nj;  grid(ib)%mj = nj
                                  grid(ib)%mk = 1
               allocate (grid(ib)%q(num_q,ni,nj,1))

               do j = 1, nj

                  if (j < nj) then  ! Must point to lower left of each cell
                     jc = j
                     q  = zero
                  else
                     jc = nj - 1
                     q  = one
                  end if

                  do i = 1, ni

                     if (i < ni) then
                        ic = i
                        p  = zero
                     else
                        ic = ni - 1
                        p  = one
                     end if

                     call surface_normal (ni, nj, &
                                          grid(ib)%x, grid(ib)%y, grid(ib)%z,  &
                                          ic, jc, p, q, un)
                     grid(ib)%q(:,i,j,1) = un(:)

                  end do

               end do

            end do

            open (lunq, file='unit_normals.f', status='unknown', iostat=ios)

            call q_write (lunq, .true., nblocks, num_q, grid, ios)

            do ib = 1, nblocks
               deallocate (grid(ib)%q)
            end do

            go to 999  ! Done

          end if

      end if

!!!   Display the block dimensions to reassure the user:

!!!   write (luncrt, '(a)')
!!!   write (luncrt, '(4(i6, 2X, 3i5))') &
!!!      (ib, grid(ib)%ni, grid(ib)%nj, grid(ib)%nk, ib = 1, nblocks)

!     We must check all blocks for right-handedness before we can write the
!     output header records.  If none need adjusting, there is no need to
!     open and write any file.

      nrectify = 0
      eps = epsilon (eps)

      do ib = 1, nblocks

         if (surface) then

            if (.not. flip(ib)) cycle

         else ! Calculate a representative cell volume, lower corner at (2,2,2)

            call cellvol (grid(ib)%ni, grid(ib)%nj, grid(ib)%nk,         &
                          2, 3, 2, 3, 2, 3,                              &
                          grid(ib)%x, grid(ib)%y, grid(ib)%z, volume, one)

            flip(ib) = volume(3,3,3) < -eps

         end if

         if (flip(ib)) then

            write (luncrt, '(a, i5)') ' Rectifying block', ib

            call rectify_block  ! Permute (i,j) in-place (internal procedure)

            nrectify = nrectify + 1

         end if

      end do

      if (nrectify > 0) then

         read (lunctl, *) filename  ! Output grid name
         read (lunctl, *) formatted

         i1 = 1;  if (formatted) i1 = 3

         open (lunxyz, file=filename, form=format(i1:11), status='unknown', &
               iostat=ios)

         if (ios /= 0) then
            write (luncrt, '(/, 2a)') &
               ' Unable to open output grid: ', filename(1:len_trim(filename))
            go to 999
         end if

         call xyz_write (lunxyz, formatted, nblocks, grid, ios)

         if (ios /= 0) go to 999

         if (qfile) then
            read (lunctl, *) filename
            read (lunctl, *) formatted

            open (lunq, file=filename, form=format(i1:11), status='unknown', &
                  iostat=ios)
            if (ios /= 0) then
               write (luncrt, '(/, 2a)') &
                  ' Unable to open output function file: ', &
                  filename(1:len_trim(filename))
                  go to 999
            end if

            call q_write (lunq, formatted, nblocks, num_q, grid, ios)

         end if

      end if

  999 continue

! *** stop ! Avoid system dependencies.

!     Internal procedure for program rectify_grid:

      contains

         subroutine rectify_block

!        Local variables:

         integer :: i, j, k

         type (grid_type) :: t

!        Execution:

         t%ni = grid(ib)%nj;  t%nj = grid(ib)%ni;  t%nk = grid(ib)%nk

         call xyz_allocate (t, ios)

         do k = 1, t%nk
            do j = 1, t%nj
               do i = 1, t%ni
                  t%x(i,j,k) = grid(ib)%x(j,i,k)
                  t%y(i,j,k) = grid(ib)%y(j,i,k)
                  t%z(i,j,k) = grid(ib)%z(j,i,k)
               end do
            end do
         end do

         deallocate (grid(ib)%x, grid(ib)%y, grid(ib)%z)

         grid(ib)%ni = t%ni;  grid(ib)%nj = t%nj

         call xyz_allocate (grid(ib), ios)

         do k = 1, t%nk
            do j = 1, t%nj
               do i = 1, t%ni
                  grid(ib)%x(i,j,k) = t%x(i,j,k)
                  grid(ib)%y(i,j,k) = t%y(i,j,k)
                  grid(ib)%z(i,j,k) = t%z(i,j,k)
               end do
            end do
         end do

         deallocate (t%x, t%y, t%z)

         if (num_q > 0) then

            t%mi = grid(ib)%mj;  t%mj = grid(ib)%mi;  t%mk = grid(ib)%nk

            call q_allocate (t, num_q, ios)

            do k = 1, t%nk
               do j = 1, t%nj
                  do i = 1, t%ni
                     t%q(:,i,j,k) = grid(ib)%q(:,j,i,k)
                  end do
               end do
            end do

            deallocate (grid(ib)%q)

            grid(ib)%mi = t%mi;  grid(ib)%mj = t%mj

            call q_allocate (grid(ib), num_q, ios)

            do k = 1, t%nk
               do j = 1, t%nj
                  do i = 1, t%ni
                     grid(ib)%q(:,i,j,k) = t%q(:,i,j,k)
                  end do
               end do
            end do

            deallocate (t%q)

         end if
         
         end subroutine rectify_block

      end program rectify_grid
