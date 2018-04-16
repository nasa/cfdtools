!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   module grid_block_utilities

!  This adaptation of the earlier surface_patch_utilities module offers some
!  similar functionality for structured multiblock volume grids and optional
!  volume data.  It was prompted by a need to scan a volume grid for its
!  bounding box (one block at a time in memory) to enable further processing
!  one block at time (at the cost of reading all blocks twice).
!
!  Where I/O is implied, files are assumed to be open ready for reading/writing
!  of the intended record(s).
!
!  Not all the surface patch functions have volume block analogues.  Conversely,
!  some of the volume block functions can be applied to surface grid patches.
!
!     block_data_range   (block)      assigns min & max values of x, y, z
!     centroids_volumes  (blockin, blockout) blockout%xyzf = centroids & volumes
!     clone_structured_grid  (nb, nf, g1, g2)      Allocates/copies header info.
!     get_bounding_boxes (lun, nb, grid, form, bbox)    b-boxes; read if lun > 0
!     reflect_block_xyz  (block, n)   reflects coordinate n about plane thru O
!     reverse_block_i    (block, nf)  reverses the i indices
!     reverse_block_j    (block, nf)  reverses the j indices
!     reverse_block_k    (block, nf)  reverses the k indices
!     scale_shift_block  (block, scale, dx, dy, dz, pout)  scales/shifts x, y, z
!     swap_block_xyz     (block, m, n) swaps coordinates m and n
!     update_block       (block1, block2, nf)    transfers block1 to block2
!     ????                            what else?
!
!  04/09/05  D. A. Saunders   Initial implementation of surface_patch_utilities.
!  07/23/10     "      "      Initial adaptation of grid_block_utilities.
!
!  Author:  David Saunders, ELORET/NASA Ames Research Center, Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use grid_block_structure ! From Tecplot_io module, or equivalent
   use xyzq_io_module       ! I/O package for PLOT3D files (for bbox option)

   implicit none

   public :: block_data_range
   public :: centroids_volumes
   public :: clone_structured_grid
   public :: get_bounding_boxes
   public :: reflect_block_xyz
   public :: reverse_block_i
   public :: reverse_block_j
   public :: reverse_block_k
   public :: scale_shift_block
   public :: swap_block_xyz
   public :: update_block

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine block_data_range (block)

!     Assign min & max values of x, y, z for the given volume block.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      type (grid_type), intent (inout) :: block

!     Local constants:

      real, parameter :: big = 1.e+32

!     Local variables:

      integer :: i, j, k

!     Execution:

      block%xmin = big;  block%xmax = -big
      block%ymin = big;  block%ymax = -big
      block%zmin = big;  block%zmax = -big

      do k = 1, block%nk
         do j = 1, block%nj
            do i = 1, block%ni
               block%xmin = min (block%x(i,j,k), block%xmin)
               block%xmax = max (block%x(i,j,k), block%xmax)
               block%ymin = min (block%y(i,j,k), block%ymin)
               block%ymax = max (block%y(i,j,k), block%ymax)
               block%zmin = min (block%z(i,j,k), block%zmin)
               block%zmax = max (block%z(i,j,k), block%zmax)
            end do
         end do
      end do

      end subroutine block_data_range

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine centroids_volumes (bin, bout)

!     For the given volume block and a companion "block" with none of its
!     fields allocated yet, assign companion dimensions to represent the number
!     of cells rather than vertices. Allocate new %x, %y, %z fields and a single
!     function field and return them with the cell centroids and volumes.
!
!     [Use of the Tecplot_io package with its grid_block_structure module does
!     not permit adding further fields to that derived data type, so the work-
!     around is to declare a second array of this data type at the higher level
!     and construct ancillary fields in the companion array elements.]
!
!     Centroids are simply the average coordinates of the vertices.
!
!     Cell volumes are calculated as the sums of six pyramids with four-sided
!     bases (adapted from the FLO87 solver of Antony Jameson).
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      type (grid_type), intent (in)    :: bin   ! %x,y,z and %ni,nj defined
      type (grid_type), intent (inout) :: bout  ! Fields are assumed to be
                                                ! unallocated on input;
                                                ! returned with cell-related
                                                ! quantities described above
!     Local constants:

      real, parameter :: eighth = 1.0/8.0, fourth = 1.0/4.0, &
                         hand   = 1.0,     factor = hand / 6.0

!     Local variables:

      integer :: i, j, k, l, m, n, ni, nj, nk
      real    :: vp1, vp2, vp3, vp4, vp5, vp6, xcg, ycg, zcg

!     Statement function (more efficient than an internal procedure):

      real :: xa, ya, za, xb, yb, zb, xc, yc, zc, xd, yd, zd, volpym

      volpym (xa, ya, za, xb, yb, zb, xc, yc, zc, xd, yd, zd) = &
             (xcg - fourth * (xa + xb + xc + xd))               &
           * ((ya - yc) * (zb - zd) - (za - zc) * (yb - yd))    &
           + (ycg - fourth * (ya + yb + yc + yd))               &
           * ((za - zc) * (xb - xd) - (xa - xc) * (zb - zd))    &
           + (zcg - fourth * (za + zb + zc + zd))               &
           * ((xa - xc) * (yb - yd) - (ya - yc) * (xb - xd))

!     volpym = volume of a pyramid with four-sided base, times 6 (Jameson/FLO87)

!     Execution:

      ni = bin%ni - 1;  bout%ni = ni;  bout%mi = ni  ! In case of I/O
      nj = bin%nj - 1;  bout%nj = nj;  bout%mj = nj
      nk = bin%nk - 1;  bout%nk = nk;  bout%mk = nk

      allocate (bout%x(ni,nj,nk), bout%y(ni,nj,nk), bout%z(ni,nj,nk), &
                bout%q(1,ni,nj,nk))

      do k = 1, nk
         n = k + 1
         do j = 1, nj
            m = j + 1
            do i = 1, ni
               l = i + 1

               xcg = (bin%x(i,j,k) + bin%x(l,j,k) + &
                      bin%x(i,m,k) + bin%x(l,m,k) + &
                      bin%x(i,j,n) + bin%x(l,j,n) + &
                      bin%x(i,m,n) + bin%x(l,m,n)) * eighth
               bout%x(i,j,k) = xcg

               ycg = (bin%y(i,j,k) + bin%y(l,j,k) + & 
                      bin%y(i,m,k) + bin%y(l,m,k) + & 
                      bin%y(i,j,n) + bin%y(l,j,n) + &
                      bin%y(i,m,n) + bin%y(l,m,n)) * eighth
               bout%y(i,j,k) = ycg

               zcg = (bin%z(i,j,k) + bin%z(l,j,k) + & 
                      bin%z(i,m,k) + bin%z(l,m,k) + & 
                      bin%z(i,j,n) + bin%z(l,j,n) + &
                      bin%z(i,m,n) + bin%z(l,m,n)) * eighth
               bout%z(i,j,k) = zcg

               vp1 = -volpym (bin%x(l,m,n), bin%y(l,m,n), bin%z(l,m,n), &
                              bin%x(l,j,n), bin%y(l,j,n), bin%z(l,j,n), &
                              bin%x(l,j,k), bin%y(l,j,k), bin%z(l,j,k), &
                              bin%x(l,m,k), bin%y(l,m,k), bin%z(l,m,k))
               vp2 =  volpym (bin%x(i,m,n), bin%y(i,m,n), bin%z(i,m,n), &
                              bin%x(i,j,n), bin%y(i,j,n), bin%z(i,j,n), &
                              bin%x(i,j,k), bin%y(i,j,k), bin%z(i,j,k), &
                              bin%x(i,m,k), bin%y(i,m,k), bin%z(i,m,k))
               vp3 = -volpym (bin%x(l,m,n), bin%y(l,m,n), bin%z(l,m,n), &
                              bin%x(l,m,k), bin%y(l,m,k), bin%z(l,m,k), &
                              bin%x(i,m,k), bin%y(i,m,k), bin%z(i,m,k), &
                              bin%x(i,m,n), bin%y(i,m,n), bin%z(i,m,n))
               vp4 =  volpym (bin%x(l,j,n), bin%y(l,j,n), bin%z(l,j,n), &
                              bin%x(l,j,k), bin%y(l,j,k), bin%z(l,j,k), &
                              bin%x(i,j,k), bin%y(i,j,k), bin%z(i,j,k), &
                              bin%x(i,j,n), bin%y(i,j,n), bin%z(i,j,n))
               vp5 = -volpym (bin%x(l,m,n), bin%y(l,m,n), bin%z(l,m,n), &
                              bin%x(i,m,n), bin%y(i,m,n), bin%z(i,m,n), &
                              bin%x(i,j,n), bin%y(i,j,n), bin%z(i,j,n), &
                              bin%x(l,j,n), bin%y(l,j,n), bin%z(l,j,n))
               vp6 =  volpym (bin%x(l,m,k), bin%y(l,m,k), bin%z(l,m,k), &
                              bin%x(i,m,k), bin%y(i,m,k), bin%z(i,m,k), &
                              bin%x(i,j,k), bin%y(i,j,k), bin%z(i,j,k), &
                              bin%x(l,j,k), bin%y(l,j,k), bin%z(l,j,k))

               bout%q(1,i,j,k) = (vp1 + vp2 + vp3 + vp4 + vp5 + vp6) * factor

            end do
         end do
      end do

      end subroutine centroids_volumes

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine clone_structured_grid (nb, nf, grid1, grid2)

!     Set up a copy of a multiblock grid by allocating the blocks and
!     transcribing the header information (only).

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in)          :: nb      ! # blocks
      integer, intent (in)          :: nf      ! nf > 0 means copy %mi, etc. too
      type (grid_type), intent (in) :: grid1(nb)
      type (grid_type), pointer     :: grid2(:)

!     Local variables:

      integer :: ib

!     Execution:

      allocate (grid2(nb))

      do ib = 1, nb
         grid2%ni = grid1%ni
         grid2%nj = grid1%nj
         grid2%nk = grid1%nk
         if (nf > 0) then
            grid2%mi = grid1%mi
            grid2%mj = grid1%mj
            grid2%mk = grid1%mk
         end if
      end do

      end subroutine clone_structured_grid

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_bounding_boxes (lun, nblocks, grid, formatted, bbox)

!     This high-level utility modularizes the chore of scanning an entire multi-
!     block grid for its data range, with the option to read and store one block
!     at a time at the application's expense of rereading blocks for later pro-
!     cessing.  It can also be used if the whole grid is already in memory.
!
!     If lun > 0, the file to be read is assumed to be open ready for reading
!     the first block; the grid dimensions, etc., have already been read and
!     the grid block array has already been allocated (nblocks is known).  The
!     coordinates of each block are allocated and deallocated here.  The file
!     should be rewound or closed by the application after calling this routine.
!
!     If lun < 0, the whole grid is assumed to be in memory, and it stays there.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in)  :: lun            ! See usage above
      integer, intent (in)  :: nblocks        ! Known number of grid blocks
      type (grid_type)      :: grid(nblocks)  ! Grid block array, dimensions set
      logical, intent (in)  :: formatted      ! T|F => formatted|unformatted
      real,    intent (out) :: bbox(6)        ! Overall grid bounding box;
                                              ! xmin, xmax = bbox(1:2), etc.
!     Local constants:

      real, parameter :: big = 1.e+32

!     Local variables:

      integer :: ib, ios, npts
      logical :: do_the_IO

!     Execution:

      do_the_IO = lun > 0

      do ib = 1, nblocks

         if (do_the_IO) then

            call xyz_allocate (grid(ib), ios)

            npts = grid(ib)%ni * grid(ib)%nj * grid(ib)%nk

            call xyz_block_io (1, lun, formatted, npts, grid(ib)%x,            &
                               grid(ib)%y, grid(ib)%z, ios)
         end if

         call block_data_range (grid(ib))  ! Sets grid(ib)%xmin, etc.

         if (do_the_IO) deallocate (grid(ib)%x, grid(ib)%y, grid(ib)%z)

      end do

      bbox(1:5:2) = big;  bbox(2:6:2) = -big

      do ib = 1, nblocks
         bbox(1) = min (bbox(1), grid(ib)%xmin)
         bbox(2) = max (bbox(2), grid(ib)%xmax)
         bbox(3) = min (bbox(3), grid(ib)%ymin)
         bbox(4) = max (bbox(4), grid(ib)%ymax)
         bbox(5) = min (bbox(5), grid(ib)%zmin)
         bbox(6) = max (bbox(6), grid(ib)%zmax)
      end do

      end subroutine get_bounding_boxes

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine reflect_block_xyz (block, n)

!     Reflect the indicated coordinate in the appropriate plane through (0,0,0).

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      type (grid_type), intent (inout) :: block
      integer, intent (in) :: n           ! n = 1, 2, 3  <->  reflect x, y, z

!     Execution:

      select case (n)
         case (1)
            block%x = -block%x
         case (2)
            block%y = -block%y
         case (3)
            block%z = -block%z
      end select

      end subroutine reflect_block_xyz

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine reverse_block_i (block, nf)

!     Reverse the i indices of a volume block, in-place.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      type (grid_type), intent (inout) :: block
      integer,          intent (in)    :: nf   ! # function values; nf >= 0

!     Local variables:

      integer :: i, ii, j, k, ni, niby2
      real    :: t, tq(nf)

!     Execution:

      ni = block%ni;  niby2 = (ni + 1) / 2

      do k = 1, block%nk
         do j = 1, block%nj
            ii = ni
            do i = 1, niby2
               t               = block%x(i, j,k)
               block%x(i, j,k) = block%x(ii,j,k)
               block%x(ii,j,k) = t
               t               = block%y(i, j,k)
               block%y(i, j,k) = block%y(ii,j,k)
               block%y(ii,j,k) = t
               t               = block%z(i, j,k)
               block%z(i, j,k) = block%z(ii,j,k)
               block%z(ii,j,k) = t
               ii = ii - 1 
            end do
         end do
      end do

      if (nf > 0) then

         ni = block%mi;  niby2 = (ni + 1) / 2

         do k = 1, block%mk
            do j = 1, block%mj
               ii = ni
               do i = 1, niby2
                  tq                = block%q(:,i, j,k)
                  block%q(:,i, j,k) = block%q(:,ii,j,k)
                  block%q(:,ii,j,k) = tq
                  ii = ii - 1  
               end do  
            end do
         end do

      end if

      end subroutine reverse_block_i

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine reverse_block_j (block, nf)

!     Reverse the j indices of a volume block, in-place.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      type (grid_type), intent (inout) :: block
      integer,          intent (in)    :: nf   ! # function values; nf >= 0

!     Local variables:

      integer :: i, j, jj, k, nj, njby2
      real    :: t, tq(nf)

!     Execution:

      nj = block%nj;  njby2 = (nj + 1) / 2

      do k = 1, block%nk
         do i = 1, block%ni
            jj = nj
            do j = 1, njby2
               t               = block%x(i,j, k)
               block%x(i,j, k) = block%x(i,jj,k)
               block%x(i,jj,k) = t
               t               = block%y(i,j, k)
               block%y(i,j, k) = block%y(i,jj,k)
               block%y(i,jj,k) = t
               t               = block%z(i,j, k)
               block%z(i,j, k) = block%z(i,jj,k)
               block%z(i,jj,k) = t
               jj = jj - 1
            end do
         end do
      end do

      if (nf > 0) then

         nj = block%mj;  njby2 = (nj + 1) / 2

         do k = 1, block%mk
            do i = 1, block%mi
               jj = nj
               do j = 1, njby2
                  tq                = block%q(:,i,j, k)
                  block%q(:,i,j, k) = block%q(:,i,jj,k)
                  block%q(:,i,jj,k) = tq
                  jj = jj - 1
               end do
            end do
         end do

      end if

      end subroutine reverse_block_j

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine reverse_block_k (block, nf)

!     Reverse the k indices of a volume block, in-place.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      type (grid_type), intent (inout) :: block
      integer,          intent (in)    :: nf   ! # function values; nf >= 0

!     Local variables:

      integer :: i, j, k, kk, nk, nkby2
      real    :: t, tq(nf)

!     Execution:

      nk = block%nk;  nkby2 = (nk + 1) / 2

      do j = 1, block%nj
         do i = 1, block%ni
            kk = nk
            do k = 1, nkby2
               t               = block%x(i,j,k )
               block%x(i,j,k ) = block%x(i,j,kk)
               block%x(i,j,kk) = t
               t               = block%y(i,j,k )
               block%y(i,j,k ) = block%y(i,j,kk)
               block%y(i,j,kk) = t
               t               = block%z(i,j,k )
               block%z(i,j,k ) = block%z(i,j,kk)
               block%z(i,j,kk) = t
               kk = kk - 1
            end do
         end do
      end do

      if (nf > 0) then

         nk = block%mk;  nkby2 = (nk + 1) / 2

         do j = 1, block%mj
            do i = 1, block%mi
               kk = nk
               do k = 1, nkby2
                  tq                = block%q(:,i,j,k )
                  block%q(:,i,j,k ) = block%q(:,i,j,kk)
                  block%q(:,i,j,kk) = tq
                  kk = kk - 1
               end do
            end do
         end do

      end if

      end subroutine reverse_block_k

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine scale_shift_block (block, scale, dx, dy, dz, block_out) 

!     Scale all of x, y, z and/or shift any of x, y, z, in-place or not:
!
!                  x_out <-- scale * x + dx          etc.
!
!     Unit scaling and zero shifts are not suppressed.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      type (grid_type), intent (inout) :: block
      real, intent (in)                :: scale, dx, dy, dz
      type (grid_type), intent (inout) :: block_out  ! May be the input block

!     Local variables:

      integer :: i, j, k

!     Execution:

      do k = 1, block%nk
         do j = 1, block%nj
            do i = 1, block%ni
               block_out%x(i,j,k) = scale * block%x(i,j,1) + dx
               block_out%y(i,j,k) = scale * block%y(i,j,1) + dy
               block_out%z(i,j,k) = scale * block%z(i,j,1) + dz
            end do
         end do
      end do

      end subroutine scale_shift_block

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine swap_block_xyz (block, m, n)

!     Exchange the indicated pair of coordinates.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      type (grid_type), intent (inout) :: block
      integer, intent (in) :: m, n        ! 1, 2, 3  <->  x, y, z

!     Local variables:

      integer :: ni, nj, nk

!     Execution:

      ni = block%ni;  nj = block%nj;  nk = block%nk

      select case (m)
         case (1)
            select case (n)
               case (2)
                  call swap (ni, nj, nk, block%x, block%y)
               case (3)
                  call swap (ni, nj, nk, block%x, block%z)
            end select
         case (2)
            select case (n)
               case (1)
                  call swap (ni, nj, nk, block%y, block%x)
               case (3)
                  call swap (ni, nj, nk, block%y, block%z)
            end select
         case (3)
            select case (n)
               case (1)
                  call swap (ni, nj, nk, block%z, block%x)
               case (2)
                  call swap (ni, nj, nk, block%z, block%y)
            end select
      end select

      contains

         subroutine swap (ni, nj, nk, cm, cn)

         integer, intent (in)    :: ni, nj, nk
         real,    intent (inout) :: cm(ni,nj,nk), cn(ni,nj,nk)
         integer :: i, j, k
         real    :: temp

!        Play safe by not assuming packed 2-D arrays.

         do k = 1, nk
            do j = 1, nj
               do i = 1, ni
                  temp      = cm(i,j,k)
                  cm(i,j,k) = cn(i,j,k)
                  cn(i,j,k) = temp
               end do
            end do
         end do

         end subroutine swap

      end subroutine swap_block_xyz

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine update_block (block_in, nf, release_in, release_out, block_out)

!     Transcribe one volume block to another.
!     The output block fields are (re)allocated here before the copy.
!     The input block fields are deallocated after the copy if indicated.
!     The "release" arguments are adopted because of observed misbehavior with
!     if (allocated(...)) ...
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      type (grid_type), intent (inout) :: block_in    ! Volume block to copy
      integer,          intent (in)    :: nf          ! # functions >= 0
      logical,          intent (in)    :: release_in  ! T => deallocate inputs
      logical,          intent (in)    :: release_out ! T => reallocate outputs
      type (grid_type), intent (inout) :: block_out   ! x

!     Execution:

      if (release_out) then
         deallocate (block_out%x, block_out%y, block_out%z)
         if (nf > 0) deallocate (block_out%q)
      end if

      block_out%ni = block_in%ni
      block_out%nj = block_in%nj
      block_out%nk = block_in%nk

      allocate (block_out%x(block_in%ni,block_in%nj,block_in%nk), &
                block_out%y(block_in%ni,block_in%nj,block_in%nk), &
                block_out%z(block_in%ni,block_in%nj,block_in%nk))

      block_out%x = block_in%x
      block_out%y = block_in%y
      block_out%z = block_in%z

      if (nf > 0) then
         block_out%mi = block_in%mi
         block_out%mj = block_in%mj
         block_out%mk = block_in%mk
         allocate (block_out%q(nf,block_in%mi,block_in%mj,block_in%mk))
         block_out%q  = block_in%q
      end if

      if (release_in) then
         deallocate (block_in%x, block_in%y, block_in%z)
         if (nf > 0) deallocate (block_in%q)
      end if

      end subroutine update_block

   end module grid_block_utilities
