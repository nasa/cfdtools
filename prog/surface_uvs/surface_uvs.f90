!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      program surface_uvs
!
!     Description:
!
!        Parameterize a multiblock surface grid by converting each patch to
!     normalized arc-length form.  For grid morphing purposes via MESHWARP, the
!     output is a 2-space iblanked PLOT3D file, where the integers appended to
!     the (u,v) coordinates for each block are all the block number.
!
!     Explanation:
!
!        MESHWARP is designed to work with a surface paneling distinct from the
!     surface grid.  UV_MAP maps the baseline grid to the baseline paneling, and
!     the iblank numbers point to the surface panel that each surface grid point
!     is found to lie in by UV_MAP.  SURFACE_UVS bypasses the UV_MAP step for
!     the case where the defining geometry is the multiblock surface grid.  The
!     (u,v) map in this case is trivial to produce, as is done here.
!
!     Further notes:
!
!        The input multiblock grid may be a surface grid or a volume grid.  In
!     the latter case, common for hypersonic applications, the surface is taken
!     to be at the k = 1 face for every block.  (A later version of this program
!     may need the block connectivity file to determine the surface patches.)
!
!     Method:
!
!        Blocks are processed one at a time.
!
!     Control file ('surface_uvs.inp'):
!
!        SURFACE_UVS control file
!        mygrid.g     Input grid
!        T            Formatted?
!        mygrid.uv    (u,v) map
!        T            Formatted?
!
!     Procedures:
!
!        XYZQ_IO package  I/O utilities for PLOT3D grid and function files
!
!     History:
!
!        03/28/05  D. Saunders  Adaptation of RECTIFY_GRID.
!
!     Author:  David Saunders, ELORET/NASA Ames Research Ctr, Moffett Field, CA
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
         lunxyz = 2,        &
         lunmap = 3

      real, parameter ::    &
         one  = 1.,         &
         zero = 0.

      character, parameter :: &
         format * 11 = 'unformatted'

!     Variables:

      integer :: &
         i, i1, ib, ios, j, nblocks, ni, nj, nk, npts

      logical :: &
         formatted, surface

      character :: &
         filename * 80

!     Composite data types:

      type (grid_type), pointer, dimension (:) :: &
         grid, uv

!     Execution:

      open (lunctl, file='surface_uvs.inp', status='old', iostat=ios)

      if (ios /= 0) then
         write (luncrt, '(/, a)') &
            ' Unable to open surface_uvs.inp control file.'
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

!     Allocate grid coordinate arrays and read the header records:

      call xyz_header_io (1, lunxyz, formatted, nblocks, grid, ios)
      if (ios /= 0) go to 999

      surface = grid(1)%nk == 1

      read (lunctl, *) filename  ! Output (u,v) map name
      read (lunctl, *) formatted

      i1 = 1;  if (formatted) i1 = 3

      open (lunmap, file=filename, form=format(i1:11), status='unknown', &
            iostat=ios)

      if (ios /= 0) then
         write (luncrt, '(/, 2a)') &
            ' Unable to open output map file: ', filename(1:len_trim(filename))
         go to 999
      end if

!     Set up the map file header:

      allocate (uv(nblocks))

      do ib = 1, nblocks
         uv(ib)%ni = grid(ib)%ni
         uv(ib)%nj = grid(ib)%nj
         uv(ib)%nk = 1
      end do

      if (formatted) then
         write (lunmap, '(i4)') nblocks
         write (lunmap, '(2i4)') (uv(ib)%ni, uv(ib)%nj, ib = 1, nblocks)
      else
         write (lunmap) nblocks
         write (lunmap) (uv(ib)%ni, uv(ib)%nj, ib = 1, nblocks)
      end if

!     Process one block at a time:

      do ib = 1, nblocks

         call xyz_allocate (grid(ib), ios)
         if (ios /= 0) go to 999

         ni = grid(ib)%ni;  nj = grid(ib)%nj;  nk = grid(ib)%nk
         npts = ni * nj * nk

         call xyz_block_io (1, lunxyz, formatted, npts,            &
                            grid(ib)%x, grid(ib)%y, grid(ib)%z, ios)
         if (ios /= 0) go to 999

!        Normalized arc lengths for the surface face:

         call xyz_allocate (uv(ib), ios)
         if (ios /= 0) go to 999

         call param2d (ni, nj, 1, ni, 1, nj, &
                       grid(ib)%x, grid(ib)%y, grid(ib)%z, &
                       uv(ib)%x, uv(ib)%y)

         deallocate (grid(ib)%x, grid(ib)%y, grid(ib)%z)

         if (formatted) then
            write (lunmap, '(1p, 6e18.11)') uv(ib)%x
            write (lunmap, '(1p, 6e18.11)') uv(ib)%y
            write (lunmap, '(20i4)') (ib, i = 1, ni * nj)
         else
            write (lunmap) uv(ib)%x, uv(ib)%y, (ib, i = 1, ni * nj)
         end if

         deallocate (uv(ib)%x, uv(ib)%y, uv(ib)%z)

      end do

  999 continue

! *** stop ! Avoid system dependencies.

      end program surface_uvs
