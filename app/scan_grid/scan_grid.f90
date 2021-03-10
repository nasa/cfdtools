!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program scan_grid
!
!  Description:
!
!        SCAN_GRID scans the specified block(s) of a grid to determine such
!     characteristics as the data range.  This version can also read a function
!     file as might be needed to locate (for instance) maximum residuals in a
!     flow solution.
!
!  Implementation Notes:
!
!        Simply scanning all blocks sequentially would allow storing of only
!     one block at a time.  However, providing an inner loop over index sub-
!     ranges of a specified block makes it awkward not to store the entire grid.
!
!     If a function file is present, the program does indeed store just one
!     block at a time and allows just one of the functions to be scanned per
!     run.  Use something like Tecplot to scrutinize particular points of
!     particular blocks if function values are needed as well as coordinates.
!
!  History:
!
!     01/11/00  D. Saunders  Initial implementation (PLOT3D grid data range),
!               ELORET Corp. using Mark Rimlinger's CFD_IO_PACKAGE.
!     08/27/03   "      "    Provide data range for each block in addition to
!                            the data range over the block range.
!                            This can help identify (say) which surface
!                            patches are which in a surface paneling.
!     10/14/03   "      "    Tabulate block face corner (x,y,z)s as a visual
!                            aid to determining connectivity.
!     12/03/03   "      "    Tabulate the shortest cell edges for all blocks
!                            in all index directions.
!     10/20/04   "      "    Show the local normal if a single surface point
!                            is specified when a subblock is scanned.
!     12/17/10  D. Saunders  Switch from cfd_io_package to xyzq_io package,
!               ERC, Inc.    and add handling of function files as needed to
!               NASA ARC     locate maximum residuals.
!     03/28/20  D. Saunders  Tabulate function minima as well as maxima,
!               AMA, Inc.    prompted by a mystery involving FLOW_INTERP output.
!
!  Author:
!
!     David Saunders, AMA, Inc. at NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_block_structure  ! Derived data type for one grid block
   use xyzq_io_module        ! PLOT3D I/O package

   implicit none

!  Constants:

   integer, parameter :: &
      lung   = 1, &
      lunf   = 2, &
      luntab = 3, &
      lunkbd = 5, &
      luncrt = 6, &
      lunmax = 7, &
      lunmin = 8

   real, parameter :: &
      one  = 1.0, &
      zero = 0.0

   logical,   parameter :: &
      false = .false., &
      true  = .true.

   character, parameter :: &
      formatted * 11 = 'unformatted'

!  Variables:

   integer :: &
      i, i1, i2, ic, ios, iq, j, j1, j2, jc, k1, k2, lun, n, n1, n2, nblock, &
      nf, ni, nj, nk, &
      imat, jmat, kmat, imax, jmax, kmax, imit, jmit, kmit, imin, jmin, kmin

   real :: &
      overall_length_scale, p, q, unit_normal(3), &
      xmin, xmax, ymin, ymax, zmin, zmax, &
      xmit, xmat, ymit, ymat, zmit, zmat

   logical :: &
      cell_centered, cr, eof, formatted_f, formatted_g, function_file_present

   type (grid_type), dimension (:), pointer :: &
      grid

   logical :: &
      first, show_surface_normal

   character :: &
      answer * 1, filename_f * 80, filename_g * 80

!  Execution:

!  Prompt for the grid file, determine its formatting, and open it:

   ios = 1

   do while (ios /= 0)

      call reads (luncrt, 'Grid file name: ', lunkbd, filename_g, cr, eof)
      if (cr .or. eof) go to 99  ! Single STOP philosophy

      call determine_grid_form (filename_g, lung, formatted_g, ios)

   end do

   i1 = 1;  if (formatted_g) i1 = 3

   open (lung, file=filename_g, form=formatted(i1:11), status='old')

!  Function file as well?

   function_file_present = false
   ios = 1

   do while (ios /= 0)

      filename_f = 'none'
      call reads (luncrt, 'Function file name? [<cr>=none]: ', &
                  lunkbd, filename_f, cr, eof)
      if (cr .or. eof) exit

      function_file_present = filename_f(1:4) /= 'none'

      if (.not. function_file_present) exit

      call determine_grid_form (filename_f, lunf, formatted_f, ios)
      
   end do

   if (function_file_present) then
      i1 = 1;  if (formatted_f) i1 = 3

      open (lunf, file=filename_f, form=formatted(i1:11), status='old')

      iq = 1
      call readi (luncrt, 'Function number to be scanned? [<cr>=1]: ', &
                  lunkbd, iq, cr, eof)

      call scan_function_data ()  ! Modularized below so the remainder can
                                  ! deal with the grid only as originally
      go to 99

   end if

!  Read the entire grid.  File attributes will be prompted for:

   call xyzq_read (lung, -lunf, formatted_g, nblock, nf, cell_centered, &
                   grid, ios)
   if (ios /= 0) go to 99

!  Display the grid block dimensions:

   if (nblock > 1) then
      write (luncrt, '(/, a)') ' Block dimensions:'
      write (luncrt, '(/, (i6, 2x, 3i5))') &
         (n, grid(n)%ni, grid(n)%nj, grid(n)%nk, n = 1, nblock)
   end if

!  Write the data range scan to the screen and to a file:

   open (luntab, file='data_range.txt', status='unknown')

!  Interactive scan of specified (sub)blocks:

   first = true

   do ! Until EOF at keyboard

      if (nblock > 1) then
         write (luncrt, '(/, a)', advance='no') &
            ' First and last block to scan: '
         read (lunkbd, *, iostat=ios) n1, n2
      else ! Avoid the redundant prompt for block number range
         n1 = 1;  n2 = 1
         if (first) then
            ios = 0
         else
50          write (luncrt, '(/, a)', advance='no') &
               ' More?  (y = yes; ^D = done) '
            read (lunkbd, *, iostat=ios) answer
            if (ios > 0) go to 50  ! Probably a bad input
            if (answer /= 'Y' .and. answer /= 'y') ios = -1
         end if
      end if

      if (ios > 0) then ! Error
         cycle
      else if (ios < 0) then ! EOF
         exit
      else if (n1 < 1 .or. n2 > nblock .or. n1 > n2) then
         write (luncrt, '(/, a)') ' Bad block range.  Try again.'
         cycle
      end if

      show_surface_normal = false

      if (n1 == n2) then ! Allow for searching a subblock

         n  = n1
         ni = grid(n)%ni;  nj = grid(n)%nj;  nk = grid(n)%nk
         write (luncrt, '(/, a, 3i6)') ' Dimensions: ', ni, nj, nk
         ios = 1

         do while (ios > 0)

            if (nk > 1) then
               write (luncrt, '(/, a)', advance='no') &
                  ' Subblock range (3 pairs): '
               read (lunkbd, *, iostat=ios) i1, i2, j1, j2, k1, k2
            else
               write (luncrt, '(/, a)', advance='no') &
                  ' Subblock range (2 pairs): '
               read (lunkbd, *, iostat=ios) i1, i2, j1, j2
               k1 = 1;  k2 = 1
            end if

            if (ios > 0) cycle
            if (ios < 0) exit

            if (i1 < 1 .or. i2 < 1 .or. i1 > ni .or. i2 > ni .or. &
                j1 < 1 .or. j2 < 1 .or. j1 > nj .or. j2 > nj .or. &
                k1 < 1 .or. k2 < 1 .or. k1 > nk .or. k2 > nk) then
               ios = 1
            end if

         end do

         call scan_grid_block (grid(n))

         xmit = xmin;  xmat = xmax
         ymit = ymin;  ymat = ymax
         zmit = zmin;  zmat = zmax

!        Show the local normal if this is one surface point:

         if (nk == 1) then ! Surface
            if (i1 == i2 .and. j1 == j2) show_surface_normal = true
         end if

      else ! More than one block specified

         xmit = huge (xmit) ! "t" = total (over all blocks from n1 to n2)
         ymit = xmit
         zmit = xmit
         xmat =-xmit
         ymat = xmat
         zmat = xmat

         do lun = luntab, luncrt, luncrt - luntab

            write (lun, '(/, 3a)') &
               ' Block            Xmin            Xmax', &
                   '              Ymin            Ymax', &
                   '              Zmin            Zmax'
         end do

         do n = n1, n2

            i1 = 1;  i2 = grid(n)%ni
            j1 = 1;  j2 = grid(n)%nj
            k1 = 1;  k2 = grid(n)%nk

            call scan_grid_block (grid(n))

            do lun = luntab, luncrt, luncrt - luntab
               write (lun, '(i4, 3(es18.7, es16.7))') &
                  n, xmin, xmax, ymin, ymax, zmin, zmax
            end do

            xmit = min (xmit, xmin);  xmat = max (xmat, xmax)
            ymit = min (ymit, ymin);  ymat = max (ymat, ymax)
            zmit = min (zmit, zmin);  zmat = max (zmat, zmax)

         end do

      end if

      overall_length_scale = &
         sqrt ((xmat - xmit)**2 + (ymat - ymit)**2 + (zmat - zmit)**2)

      do lun = luntab, luncrt, luncrt - luntab

         write (lun, '(/, (2x, a,  es16.7, es18.7))') &
            'Overall Xmin & Xmax:', xmit, xmat, &
            '        Ymin & Ymax:', ymit, ymat, &
            '        Zmin & Zmax:', zmit, zmat
         write (lun, '(/, a, es15.7)') &
            '  Overall length scale:', overall_length_scale
      end do

      if (show_surface_normal) then  ! Be sure to point to lower-left vertex

         if (i1 < ni) then
            p = zero;  ic = i1
         else
            p = one;   ic = ni - 1
         end if

         if (j1 < nj) then
            q = zero;  jc = j1
         else
            q = one;   jc = nj - 1
         end if

         n = n1

         call surface_normal (ni, nj, grid(n)%x, grid(n)%y, grid(n)%z, &
                              ic, jc, p, q, unit_normal)

         do lun = luntab, luncrt, luncrt - luntab
            write (lun, '(/, a, 3f11.7)') '  Local unit normal:', unit_normal
         end do

      end if

      first = false

   end do ! Next search range

   close (luntab)

!  Tabulate corner (x,y,z)s by block and face, to help with connectivity:

   call tabulate_corners ()  ! Gives trouble

!  Tabulate the shortest cell edges for all blocks in all directions:

   call tabulate_shortest_edges ()

99 continue


!  Internal procedures for SCAN_GRID:

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine scan_function_data ()

!     Processing one grid and function block at a time, tabulate the maximum
!     value of the indicated function for each block, and its coordinates.
!     Tabulate the minimum value similarly.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local variables:

      integer :: i, j, k, n, nf, nmat, nmit, npts
      real    :: f, fmat  ! Max. over all blocks (t = total)
      real    ::    fmit  ! Min. over all blocks (t = total)
      real    :: fmax     ! Max. within a block
      real    :: fmin     ! Min. within a block

!     Execution:

      open (lunmax, file='max_function_values.txt', status='unknown')
      open (lunmin, file='min_function_values.txt', status='unknown')

      write (lunmax, '(2a)') 'Block         Maximum    I    J    K',           &
                 '               X               Y               Z'
      write (lunmin, '(2a)') 'Block         Minimum    I    J    K',           &
                 '               X               Y               Z'

      call xyz_header_io (1, lung, formatted_g, nblock, grid, ios)
      if (ios /= 0) go to 99

      call q_header_io (1, lunf, formatted_f, nblock, nf, grid, ios)
      if (ios /= 0) go to 99

      fmat = -huge (fmat);  fmit = -fmat

      do n = 1, nblock

         call xyz_allocate (grid(n), ios)
         if (ios /= 0) go to 99

         ni = grid(n)%ni
         nj = grid(n)%nj
         nk = grid(n)%nk

         npts = ni*nj*nk

         call xyz_block_io (1, lung, formatted_g, npts, &
                            grid(n)%x, grid(n)%y, grid(n)%z, ios)
         if (ios /= 0) go to 99

         call q_allocate (grid(n), nf, ios)
         if (ios /= 0) go to 99

         call q_block_io (1, lunf, formatted_f, nf, ni, nj, nk, grid(n)%q, ios)
         if (ios /= 0) go to 99

         fmax = -huge (fmax);  fmin = -fmax
         do k = 1, nk
            do j = 1, nj
               do i = 1, ni
                  f = grid(n)%q(iq,i,j,k)
                  if (f > fmax) then
                      fmax = f;  imax = i;  jmax = j;  kmax = k
                      xmax = grid(n)%x(i,j,k)
                      ymax = grid(n)%y(i,j,k)
                      zmax = grid(n)%z(i,j,k)
                  else if (f < fmin) then
                      fmin = f;  imin = i;  jmin = j;  kmin = k
                      xmin = grid(n)%x(i,j,k)
                      ymin = grid(n)%y(i,j,k)
                      zmin = grid(n)%z(i,j,k)
                  end if
               end do
            end do
         end do

         write (lunmax, '(i5, es16.8, 3i5, 3es16.8)') &
            n, fmax, imax, jmax, kmax, xmax, ymax, zmax
         write (lunmin, '(i5, es16.8, 3i5, 3es16.8)') &
            n, fmin, imin, jmin, kmin, xmin, ymin, zmin

         if (fmax > fmat) then
             fmat = fmax;  imat = imax;  jmat = jmax;  kmat = kmax;  nmat = n
             xmat = xmax;  ymat = ymax;  zmat = zmax
         end if

         if (fmin < fmit) then
             fmit = fmin;  imit = imin;  jmit = jmin;  kmit = kmin;  nmit = n
             xmit = xmin;  ymit = ymin;  zmit = zmin
         end if

         deallocate (grid(n)%x, grid(n)%y, grid(n)%z, grid(n)%q)

      end do

      write (lunmax, '(/, i5, es16.8, 3i5, 3es16.8)') &
         nmat, fmat, imat, jmat, kmat, xmat, ymat, zmat
      write (lunmin, '(/, i5, es16.8, 3i5, 3es16.8)') &
         nmit, fmit, imit, jmit, kmit, xmit, ymit, zmit

      close (lunmax)
      close (lunmin)

      deallocate (grid)

   99 return

      end subroutine scan_function_data

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine scan_grid_block (block)

!     Determine the data range for the (sub)block defined by the argument block
!     and the current values of i1, i2, j1, j2, k1, k2.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      type (grid_type), intent (in) :: block

!     Local variables:

      integer :: i, j, k

!     Execution:

      xmin = huge (xmin)
      ymin = xmin
      zmin = xmin
      xmax =-xmin
      ymax = xmax
      zmax = xmax

      do k = k1, k2
         do j = j1, j2
            do i = i1, i2
               xmin = min (xmin, block%x(i,j,k))
               xmax = max (xmax, block%x(i,j,k))
               ymin = min (ymin, block%y(i,j,k))
               ymax = max (ymax, block%y(i,j,k))
               zmin = min (zmin, block%z(i,j,k))
               zmax = max (zmax, block%z(i,j,k))
            end do
         end do
      end do

      end subroutine scan_grid_block

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine tabulate_corners ()

!     Tabulate corner (x,y,z)s by block/face to help with connectivity.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: n

      open (luntab, file='corner_xyzs.txt', status='unknown')

      write (luntab, '(/, 2a)') &
         ' Corner coordinates for grid ', trim (filename_g)

      do n = 1, nblock

         call block_corners (n, grid(n))

      end do

      close (luntab)

      end subroutine tabulate_corners

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine block_corners (n, block)

!     Ancillary for TABULATE_CORNERS, (originally) allowing 3-D indexing.
!     Face numbers are 1-6 for i = 1, i = ni, ...., k = 1, k = nk respectively.
!     Corners of a given face are tabulated cyclically: (j,k), (k,i), (i,j)

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer,          intent (in) :: n      ! Block number
      type (grid_type), intent (in) :: block  ! Current grid block

      integer, parameter :: nchar = 51 ! Per group

      integer   :: ni, nj, nk
      integer   :: i, j, k, kinc, l, m1, m2
      character :: line * (4 * nchar)

!     Execution:

      ni = block%ni;  nj = block%nj;  nk = block%nk

      write (luntab, '(/, 4(a, i4))') &
         ' Block: ', n, '  Dimensions: ', ni, '  x ', nj, '  x ', nk

      if (nk > 1) then
         kinc = nk - 1
      else
         kinc = 1
      end if

      l = 2 * nchar

      if (nk > 1) then ! Else there is only one face at k = 1

!        Faces 1 (i1) and 2 (i2):

         do i = 1, ni, ni - 1
            m1 = 1; m2 = nchar
            do k = 1, nk, kinc
               do j = 1, nj, nj - 1
                  write (line(m1:m2), '(3i4, 3f13.6)') &
                     i, j, k, block%x(i,j,k), block%y(i,j,k), block%z(i,j,k)
                  m1 = m2 + 1; m2 = m2 + nchar
               end do
            end do
            write (luntab, '(/, (a))') line(1:l), line(l+1:2*l)
         end do

!        Faces 3 (j1) and 4 (j2):

         do j = 1, nj, nj - 1
            m1 = 1; m2 = nchar
            do i = 1, ni, ni - 1
               do k = 1, nk, kinc
                  write (line(m1:m2), '(3i4, 3f13.6)') &
                     i, j, k, block%x(i,j,k), block%y(i,j,k), block%z(i,j,k)
                  m1 = m2 + 1; m2 = m2 + nchar
               end do
            end do
            write (luntab, '(/, (a))') line(1:l), line(l+1:2*l)
         end do

      end if

!     Faces 5 (k1) and 6 (k2):

      do k = 1, nk, kinc
         m1 = 1; m2 = nchar
         do j = 1, nj, nj - 1
            do i = 1, ni, ni - 1
               write (line(m1:m2), '(3i4, 3f13.6)') &
                  i, j, k, block%x(i,j,k), block%y(i,j,k), block%z(i,j,k)
               m1 = m2 + 1; m2 = m2 + nchar
            end do
         end do
         write (luntab, '(/, (a))') line(1:l), line(l+1:2*l)
      end do

      end subroutine block_corners

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine tabulate_shortest_edges ()

!     Tabulate the smallest cell edges for all blocks.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, parameter :: big = 1.e+30

      integer :: i, j, k, ii, ji, ki, ij, jj, kj, ik, jk, kk, n
      real    :: edge, edgei, edgej, edgek

!     Execution:

      open (luntab, file='smallest_edges.txt', status='unknown')

      write (luntab, '(/, 2a)') &
         ' Smallest edges for grid ', trim (filename_g),   &
         ' BLOCK   NI  NJ  NK    I DIRECTION   I   J   K', &
         '    J DIRECTION   I   J   K    K DIRECTION   I   J   K'

      do n = 1, nblock

         ni = grid(n)%ni;  nj = grid(n)%nj;  nk = grid(n)%nk

!        I direction:

         edgei = big
         do k = 1, nk
            do j = 1, nj
               do i = 1, ni - 1
                  edge = (grid(n)%x(i+1,j,k) - grid(n)%x(i,j,k)) ** 2 + &
                         (grid(n)%y(i+1,j,k) - grid(n)%y(i,j,k)) ** 2 + &
                         (grid(n)%z(i+1,j,k) - grid(n)%z(i,j,k)) ** 2
                  if (edge < edgei) then
                     edgei = edge
                     ii = i; ji = j; ki = k
                  end if
               end do
            end do
         end do
         edgei = sqrt (edgei)

!        J direction:

         edgej = big
         do k = 1, nk
            do i = 1, ni
               do j = 1, nj - 1
                  edge = (grid(n)%x(i,j+1,k) - grid(n)%x(i,j,k)) ** 2 + &
                         (grid(n)%y(i,j+1,k) - grid(n)%y(i,j,k)) ** 2 + &
                         (grid(n)%z(i,j+1,k) - grid(n)%z(i,j,k)) ** 2
                  if (edge < edgej) then
                     edgej = edge
                     ij = i; jj = j; kj = k
                  end if
               end do
            end do
         end do
         edgej = sqrt (edgej)

!        K direction?

         if (nk == 1) then
            edgek = 0.
            ik = 0;  jk = 0;  kk = 0
         else
            edgek = big
            do j = 1, nj
               do i = 1, ni
                  do k = 1, nk - 1
                     edge = (grid(n)%x(i,j,k+1) - grid(n)%x(i,j,k)) ** 2 + &
                            (grid(n)%y(i,j,k+1) - grid(n)%y(i,j,k)) ** 2 + &
                            (grid(n)%z(i,j,k+1) - grid(n)%z(i,j,k)) ** 2
                     if (edge < edgek) then
                        edgek = edge
                        ik = i; jk = j; kk = k
                     end if
                  end do
               end do
            end do
            edgek = sqrt (edgek)
         end if

         write (luntab, '(i4, 3x, 3i4, 3(es15.6, 3i4))')      &
            n, ni, nj, nk, edgei, ii, ji, ki, edgej, ij, jj, kj, &
            edgek, ik, jk, kk

      end do

      close (luntab)

      end subroutine tabulate_shortest_edges

   end program scan_grid
