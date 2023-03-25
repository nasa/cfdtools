!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      program refine_grid
!
!     Description:
!
!        REFINE_GRID densifies all blocks of a multiblock grid the same way, at
!        least initially.  If different treatment for different blocks is ever
!        needed, we can cross that bridge when we come to it.  For now, use of
!        other available utilities may suffice to extract/refine/merge groups of
!        blocks appropriately.
!
!        This version handles 2D as well as 3D structured grids.
!
!        If grid densification were limited to integer multiples of the input
!        cell counts, treatment of function files could be avoided, since DPLR's
!        FCONVERT can perform grid sequencing on solution files.  However, since
!        fractional multiples of cell counts are permitted here, the option to
!        process function files as well as grid files is provided.
!
!        Any function file should have the same block dimensions and formatting
!        as the grid file.  One may choose to work with cell-centered function
!        data and corresponding cell-centered grid, but keep in mind that if
!        mere doubling of dimensions, say, is all that is required, then a DPLR
!        user should be able to employ FCONVERT's options and avoid processing
!        function files here.  Processing just a vertex-centered grid with no
!        function file is expected to be the most common usage.
!
!        The three index directions are treated independently (possibly with
!        different degrees of refinement), meaning it is up to the user to know
!        that different treatment will not produce mismatches at common block
!        faces.
!
!        Any sensible multiples of the input grid block dimensions (possibly
!        non-integer multiples greater or less than 1) may be specified to
!        define the output dimensions.  The same scaling applies to all blocks.
!
!        Nothing more than 1-D spline interpolation is provided, multi-dimen-
!        sional nonlinear interpolation being unreliable in the author's
!        experience.  Indices are treated in the i-j-k order, and each new
!        set of interpolations operates on the previous interpolated values.
!        This asymmetry may lead to slight mismatches at common faces even if
!        the point counts match.  [There are no easy answers here.]
!
!        [Originally:]
!        Arc-length-based monotonic Hermite cubic spline interpolation is
!        employed as in OUTBOUND's option to adjust the number of radial grid
!        points.  (Use OUTBOUND rather than REFINE_GRID if the grid topology
!        allows it, and the off-wall direction is the only one of interest.)
!        Note that 2nd-derivative continuity at the data points/spline knots is
!        NOT guaranteed by such splines, which is likely to introduce noisy
!        curvature at geometric surfaces.
!
!        [Currently:}
!        Conventional cubic splines are now used for the (x,y,z) coordinates.
!        Local cubic splines still suffice for flow solution interpolations.
!
!        For a completely general grid with arbitrary topology, treating all
!        index directions identically may be the only safe operation.  For
!        cases commonly treated in hypersonic blunt body analysis (single layer
!        of grid blocks, with k = 1 at the wall), an OUTBOUND-type option is
!        provided, namely preserving the input spacing off the wall rather than
!        preserving the existing relative spacing.  Use it with discretion!
!
!     Control File Format (on standard input):
!
!        [Note that the formatting of the input grid is determined by this
!         program.  For the output result, use *.g/gu and *.f/fuappropriately.]
!
!        REFINE_GRID control file. For output formatting, use *.g or *.gu names.
!        my_coarse_grid.gu   ! Input grid name
!        my_coarse_f.fu      ! Input function file name, or "none"
!        my_refined_grid.gu  ! Output grid name (*.g | *.gu)
!        my_refined_f.fu     ! Output function file name (*.f | *.fu | none)
!        2 2 1.5             ! si/sj/sk applied to cell counts; may be < 1.
!        1 1 2               ! mi/mj/mk: 1 retains relative spacing; 2 keeps ds1
!        1. 1. 0.2           ! ds2_frctns (mi=2): ds2i = 1-sided ds2 * ds2_frctn
!
!     Procedures:
!
!        XYQ_IO  module      2D analog of XYZQ_IO
!        XYZQ_IO package     I/O utilities for PLOT3D grid and function files
!        DENSIFY_GRID_BLOCK  Adjust the number of cells in i/j/k directions
!        CHORDS3D            Arc lengths along a 3-space curve
!        EXPDIS5             One-sided stretching (Vinokur/geometric hybrid)
!        BLGRID              Two-sided stretching (Vinokur)
!     !! PLSCRV3D            Arc-length-based spline interpolation, 3-space crv.
!        CSPLINE             Conventional cubic spline for the grid lines
!        LSCFIT              Monotonic Hermite (local) cubic spline for the flow
!
!     History:
!
!        11/19/09  D.A.Saunders  Adapted OUTBOUND's radial point count option.
!        11/21/09    "     "     Introduced the cspline_module and applied it
!                                to the grid coordinates (only) to improve
!                                curvature continuity on geometric surfaces.
!                                Block boundary effects are still likely.
!        11/23/09    "     "     If a scale factor is 2., preserve the odd/
!                                evenness of the input dimension.  I.e.,
!                                odd n --> 2n - 1, but even n --> 2n for the
!                                case of cell-centered data.  If it is 0.5,
!                                use of THIN_GRID & THIN_FLOW is recommended.
!                                Todd White suggested checking for blocks with
!                                negative cell volumes.
!        09/29/16    "     "     Aaron Brandis noticed that a scale factor of 4
!                                didn't preserve oddness as expected, the way 2
!                                does. Therefore, treat all integer multipliers
!                                the same way.  Remember, even point counts are
!                                still likely for cell-centered data.
!        09/23/20  Jeff Hill     Updated file handling to support 2D grids.
!
!     Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA
!              Now with AMA, Inc. at NASA ARC.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_block_structure  ! Derived data type for a grid block
   use xyzq_io_module        ! PLOT3D file I/O package
   use xyq_io_module         ! 2D grid I/O

   implicit none

!  Constants:

   integer, parameter :: &
      luning  = 1,       &  ! Input grid
      luninf  = 2,       &  ! Input function file, if any
      lunoutg = 3,       &  ! Output grid
      lunoutf = 4,       &  ! Output function file, if any
      lunctl  = 5,       &  ! Control file (on command line)
      luncrt  = 6           ! Screen diagnostics

!  Variables:

   integer :: &
      ib, ios, method_i, method_j, method_k, mi, mj, mk, &
      nblocks, ndigits, nf, ni, nj, nk, npad, npts, nspaces

   real :: &
      ds2_fraction_i, ds2_fraction_j, ds2_fraction_k, scale_i, scale_j, scale_k

   logical :: &
      threed, formatted_in, formatted_out

   character :: &
      b_format*41

!  Derived data types:

   type (grid_type), pointer, dimension (:) :: &
      gridin, gridout

!  Execution:
!  ----------

!  Read the controls from standard input:

   call read_controls ()  ! Local procedure below

   if (ios /= 0) go to 99

!  Check for thinning dimensions by an integer multiple (unnecessary round-off):

   call scale_check (scale_i)  ! Local procedure; may set ios /= 0 and
   call scale_check (scale_j)  ! urge use of THIN_GRID/THIN_FLOW
   call scale_check (scale_k)

   if  (ios /= 0) go to 99

!  The files are open.  Read the input grid header and allocate work-space:

   threed = .true.
   call xyz_header_io (1, luning, formatted_in, nblocks, gridin, ios)
   if (ios /= 0) then
      write(*,*) "Retry grid header read assuming 2D Plot3D format."
      rewind luning
      call xy_header_io (1, luning, formatted_in, nblocks, gridin, ios)
      if (ios /= 0) go to 99
      write(*,*) "The grid file appears to be two dimensional."
      gridin%nk = 1  ! Not initialized by xy_header_io
      threed = .false.
   endif

!  Read the input function file header, if present, and ensure dimensions match:

   if (nf > 0) then

      write(*,*) "Loading solution file"
      call q_header_io (1, luninf, formatted_in, nblocks, nf, gridin, ios)
      if (ios /= 0) then
         write(*,*) "Retry solution header read assuming 2D Plot3D format."
         rewind luninf
         call q_header_io_2d (1, luninf, formatted_in, nblocks, nf, gridin, ios)
         write(*,*) "The solution file appears to be two dimensional."
         if (ios /= 0) go to 99
         gridin%mk = 1
      end if

      do ib = 1, nblocks
         if (gridin(ib)%ni /= gridin(ib)%mi) then
             write(*,*) "ni mismatch for block ", ib, ": ", gridin(ib)%ni, gridin(ib)%mi
             ios = 1
         end if
         if (gridin(ib)%nj /= gridin(ib)%mj) then
             write(*,*) "nj mismatch for block ", ib, ": ", gridin(ib)%nj, gridin(ib)%mj
             ios = 1
         end if
         if (gridin(ib)%nk /= gridin(ib)%mk) then
             write(*,*) "nk mismatch for block ", ib, ": ", gridin(ib)%nk, gridin(ib)%mk
             ios = 1
         end if
      end do

      if (ios /= 0) then
         write (luncrt, '(a)') 'Function file dimensions must match the grid.'
         go to 99
      end if

   end if

!  Set up the output block dimensions:

   allocate (gridout(nblocks))

   do ib = 1, nblocks  ! See local procedure below
     call scale_dimensions (scale_i, gridin(ib)%ni, gridout(ib)%ni)
     call scale_dimensions (scale_j, gridin(ib)%nj, gridout(ib)%nj)
     call scale_dimensions (scale_k, gridin(ib)%nk, gridout(ib)%nk)
     gridout(ib)%mi = gridout(ib)%ni
     gridout(ib)%mj = gridout(ib)%nj  ! Flow dimensions must match the grid
     gridout(ib)%mk = gridout(ib)%nk
   end do

   if (threed) then
      call xyz_header_io (2, lunoutg, formatted_out, nblocks, gridout, ios)
      if (ios /= 0) go to 99
      if (nf > 0) then
         call q_header_io (2, lunoutf, formatted_out, nblocks, nf, gridout, ios)
         if (ios /= 0) go to 99
      end if
   else
      call xy_header_io (2, lunoutg, formatted_out, nblocks, gridout, ios)
      if (ios /= 0) go to 99
      if (nf > 0) then
         call q_header_io_2d (2, lunoutf, formatted_out, nblocks, nf, gridout, ios)
         if (ios /= 0) go to 99
      end if
   end if

   b_format = '(a, i?, ":", ?x, 3i4, a, 3i5, ",", i4, a)'
!              12345678901234567890123456789012345678901

   if (nblocks < 10) then  ! Avoid excessive spaces in the screen dialogue
      nspaces = 2
   else if (nblocks < 100) then
      nspaces = 3
   else if (nblocks < 1000) then
      nspaces = 4
   else
      nspaces = 5
   end if

   write (luncrt, '(a)')

!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Process the blocks one at a time:
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   do ib = 1, nblocks

      call xyz_allocate (gridin(ib), ios)
      if (ios /= 0) go to 99

      ni = gridin(ib)%ni;  nj = gridin(ib)%nj;  nk = gridin(ib)%nk
      npts = ni * nj * nk

      if (threed) then
         call xyz_block_io (1, luning, formatted_in, npts, &
                            gridin(ib)%x, gridin(ib)%y, gridin(ib)%z, ios)
      else
         call xy_block_io (1, luning, formatted_in, npts, &
                             gridin(ib)%x, gridin(ib)%y, ios)
         gridin(ib)%z = 0.0d0
      end if
      if (ios /= 0) go to 99

      call xyz_allocate (gridout(ib), ios)
      if (ios /= 0) go to 99

      if (nf > 0) then

         call q_allocate (gridin(ib), nf, ios)
         if (ios /= 0) go to 99

         call q_block_io (1, luninf, formatted_in, nf, ni, nj, nk,             &
                          gridin(ib)%q, ios)
         if (ios /= 0) go to 99

         call q_allocate (gridout(ib), nf, ios)
         if (ios /= 0) go to 99

      else  ! Avoid unallocated arguments in the following call
         allocate (gridin(ib)%q(1,1,1,1), gridout(ib)%q(1,1,1,1))
      end if

      mi = gridout(ib)%ni;  mj = gridout(ib)%nj;  mk = gridout(ib)%nk

      if (ib < 10) then
         ndigits = 2;  npad = nspaces
      else if (ib < 100) then
         ndigits = 3;  npad = nspaces - 1
      else if (ib < 1000) then
         ndigits = 4;  npad = nspaces - 2
      else
         ndigits = 5;  npad = nspaces - 3
      end if

      write (b_format(6:6),   '(i1)') ndigits
      write (b_format(14:14), '(i1)') npad

      write (luncrt, b_format) &
       'Processing block', ib, ni, nj, nk, '  -->', mi, mj, mk, nf, ' functions'

      call densify_grid_block (ib, nf, method_i, method_j, method_k,           &
                     ds2_fraction_i, ds2_fraction_j, ds2_fraction_k,           &
        ni, nj, nk, gridin(ib)%x,  gridin(ib)%y,  gridin(ib)%z,  gridin(ib)%q, &
        mi, mj, mk, gridout(ib)%x, gridout(ib)%y, gridout(ib)%z, gridout(ib)%q,&
        ios)

      if (ios /= 0) go to 99

      deallocate (gridin(ib)%x,  gridin(ib)%y,  gridin(ib)%z,  stat=ios)

      npts = mi * mj * mk

      if (threed) then
         call xyz_block_io (2, lunoutg, formatted_out, npts, &
                            gridout(ib)%x, gridout(ib)%y, gridout(ib)%z, ios)
      else
         call xy_block_io (2, lunoutg, formatted_out, npts, &
                            gridout(ib)%x, gridout(ib)%y, ios)
      end if
      if (ios /= 0) go to 99

!     Warn the user if any cells have negative volumes (local procedure below):

      if (nk /= 1) then
         call volume_check (mi, mj, mk, &
                            gridout(ib)%x, gridout(ib)%y, gridout(ib)%z)
      end if

      deallocate (gridout(ib)%x, gridout(ib)%y, gridout(ib)%z, stat=ios)

      if (nf > 0) then
         call q_block_io (2, lunoutf, formatted_out, nf, mi, mj, mk,           &
                          gridout(ib)%q, ios)
         if (ios /= 0) go to 99
      end if

      deallocate (gridin(ib)%q, gridout(ib)%q)

   end do

   deallocate (gridin, gridout)

   close (luning)
   close (lunoutg)

   if (nf > 0) then
      close (luninf)
      close (lunoutf)
   end if

99 continue

! *** stop ! Avoid system dependencies.

!  Internal procedures for program refine_grid, in the order they're used:

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine read_controls ()

!     Read REFINE_GRID controls from standard input.
!     Most variables are inherited from the main program.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local constants:

      real, parameter :: &
         one = 1.0

      character, parameter :: &
         format*11 = 'unformatted', &
         none*4    = 'none'

!     Local variables:

      integer   :: i1, l
      character :: filename * 128

!     Execution:

      read (lunctl, '(a)', iostat=ios)  ! Header

      if (ios /= 0) then
         write (luncrt, '(a)') 'Unable to read indicated control file.'
         go to 99  ! One-return philosophy
      end if

!     *** Input grid ***

      read (lunctl, *, iostat=ios) filename

      if (ios /= 0) then
         write (luncrt, '(a)') 'Trouble reading input grid file name.'
         go to 99
      end if

      call determine_grid_form (filename, luning, formatted_in, ios)

      if (ios /= 0) then
         write (luncrt, '(2a)') 'Unable to open input grid, ', &
            trim (filename)
         go to 99
      end if

      i1 = 1;  if (formatted_in) i1 = 3

      open (luning, file=filename, form=format(i1:11), status='old')

!     *** Input function file ***

      read (lunctl, *, iostat=ios) filename

      if (ios /= 0) then
         write (luncrt, '(a)') 'Trouble reading input function file name.'
         go to 99
      end if

      nf = 0
      if (trim (filename) /= none) then
         open (luninf, file=filename, form=format(i1:11), status='old', &
               iostat=ios)
         if (ios /= 0) then
            write (luncrt, '(2a)') 'Unable to open input function file, ', &
               trim (filename)
            go to 99
         end if
         nf = 1  ! Or more (later)
      end if

!     *** Output grid ***

      read (lunctl, *, iostat=ios) filename

      if (ios /= 0) then
         write (luncrt, '(a)') 'Trouble reading output grid file name.'
         go to 99
      end if

      l = len_trim (filename)
      formatted_out = filename(l-2:l) /= '.gu'
      i1 = 1;  if (formatted_out) i1 = 3

      open (lunoutg, file=filename(1:l), form=format(i1:11), status='unknown')

!     *** Output function file ***

      read (lunctl, *, iostat=ios) filename

      if (ios /= 0) then
         write (luncrt, '(a)') 'Trouble reading output function file name.'
         go to 99
      end if

      if (nf /= 0) then
         open (lunoutf, file=filename, form=format(i1:11), status='unknown', &
               iostat=ios)
         if (ios /= 0) then
            write (luncrt, '(2a)') 'Unable to open output function file, ', &
               trim (filename)
            go to 99
         end if
      end if

!     *** Block dimension scale factors: ***

      read (lunctl, *, iostat=ios) scale_i, scale_j, scale_k  ! May be fractions

      if (ios /= 0) then
         write (luncrt, '(a)') 'Trouble reading the cell count scale factors.'
         go to 99
      end if

!     *** Point redistribution methods in each index direction: ***

      read (lunctl, *, iostat=ios) method_i, method_j, method_k  ! 1 or 2

      if (ios /= 0) then
         write (luncrt, '(a)') 'Trouble reading the redistribution methods.'
         go to 99
      end if

      if (min (method_i, method_j, method_k) < 1 .or. &
          max (method_i, method_j, method_k) > 2) then
         ios = -1
         write (luncrt, '(a)') 'The redistribution methods should be 1 or 2.'
         go to 99
      end if

      if (method_i * method_j * method_k /= 1) then
         read (lunctl, *, iostat=ios) ds2_fraction_i, ds2_fraction_j, &
                                      ds2_fraction_k
         if (ios /= 0) then
            write (luncrt, '(a)') 'Trouble reading the ds2 fraction controls.'
            go to 99
         end if
      else  ! Don't leave them undefined
         ds2_fraction_i = one;  ds2_fraction_j = one;  ds2_fraction_k = one
      end if

  99  return

      end subroutine read_controls

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine scale_check (s)

!     Guard against thinning by integer multiples.  Set ios if necessary.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      real, intent (in) :: s     ! Scale factor to check

!     Local constants:

      real, parameter :: eps = 1.e-4  ! For comparison against a whole number
      real, parameter :: one = 1.

!     Local variables:

      integer :: i
      real    :: r

!     Execution:

      if (s < one) then
         r = one / s
         i = nint (r)
         if (abs (r - real (i)) < eps) then
            write (luncrt, '(/, a, f9.6, a, i2, a, /, a)') &
               '*** Scale factor', s, ' implies thinning/every nth point, n =',&
               i, '.', '    Use THIN_GRID/THIN_FLOW for that.'
            ios = 1
         end if
      end if

      end subroutine scale_check

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine scale_dimensions (s, n, nnew)

!     Scale the indicated dimension, with special treatment of integer scaling.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      real,    intent (in)  :: s     ! Nominal scale factor; s >=< 1.
      integer, intent (in)  :: n     ! Input dimension
      integer, intent (out) :: nnew  ! Scaled dimension

!     Local variables:

      integer :: int_s

!     Execution:

      nnew  = nint (s * real (n))
      int_s = int (s)

      if (real (int_s) == s) then
         if (mod (n, 2) == 1) nnew = int_s*(n-1) + 1
      end if

      end subroutine scale_dimensions

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine volume_check (ni, nj, nk, x, y, z)

!     Check the given grid block for negative cell volumes, and issue a warning
!     if any.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in)                       :: ni, nj, nk   ! Block dims.
      real,    intent (in), dimension (ni,nj,nk) :: x, y, z      ! Grid coords.

!     Local constants:

      real, parameter :: hand = 1.  ! The block is assumed to be right-handed
      real, parameter :: zero = 0.

!     Local variables:

      integer :: i, j, k, iv, jv, kv, nv
      real    :: v, vmin
      real, allocatable, dimension (:,:,:) :: vol

!     Execution:

      allocate (vol(2:ni,2:nj,2:nk))

      call cellvol (ni, nj, nk, 1, ni, 1, nj, 1, nk, x, y, z, vol, hand)

      nv = 0;  vmin = 1.E+20

      do k = 2, nk
         do j = 2, nj
            do i = 2, ni
               v = vol(i,j,k)
               if (v <= zero) then
                  nv = nv + 1
                  if (v <  vmin) then
                     vmin = v
                     iv = i;  jv = j;  kv = k
                  end if
               end if
            end do
         end do
      end do

      deallocate (vol)

      if (nv > 0) write (luncrt, '(/, a, i8, 2x, 3i6, /)') &
         '*** # cells with volumes < 0, and worst-case lower/left indices:',   &
         nv, iv-1, jv-1, kv-1

      end subroutine volume_check

   end program refine_grid
