!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   module grid_block_structure

!  Data structure adopted for compatibility with NASA Langley:

   public :: grid_type

   type grid_type
      real, dimension(:,:,:),   pointer :: x,y,z ! grid coordinates
      real, dimension(:,:,:,:), pointer :: q     ! dependent variables
      integer   :: ni          ! number of nodes in i direction
      integer   :: nj          ! number of nodes in j direction
      integer   :: nk          ! number of nodes in k direction
      integer   :: mdep        ! number of dependent variables
      integer   :: mi          ! number of dependent variables in i direction
      integer   :: mj          ! number of dependent variables in j direction
      integer   :: mk          ! number of dependent variables in k direction
   end type grid_type

   end module grid_block_structure

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   program grid_morph_1

!  This program performs computational grid morphing for a specialized situation.  It was prompted by the need to analyze the
!  effects of minor damage to the Space Shuttle thermal protection system, and to analyze possible repairs to such damage.
!
!  Intended application:
!
!     >  Minor dings (relatively smooth gouges, not cavities)
!     >  Repairs (relatively smooth protuberances)
!
!  Assumptions:
!
!     >  The baseline grid has been densified in the damage region (interactive capability at NASA LaRC)
!     >  The damage is defined by a structured surface grid (one block) covering some or all of the densified block surface faces
!     >  Blocks being morphed have k = 1 at the body surface; if not, each block is permuted in place here, then permuted back
!     >  All files are PLOT3D-type, multiblock
!
!  Algorithm:
!
!     >  For each block specified in the densified region of the baseline grid
!
!        >  For each point on the baseline surface
!
!           >  Locate the nearest quad. cell of the structured damaged surface
!
!           >  If it is within the damaged surface (p and q in [0, 1]) then
!                 The perturbed surface grid point is (under run-time control)  EITHER
!                    The foot of the projection to the damaged surface  OR
!                    The intersection of the baseline grid radial line with the damaged surface
!              else
!                 The surface grid point is not moved
!
!        >  Update the block interior:  EITHER
!              Treat each affected radial line independently  OR
!              Update just the radial faces (either one line at a time or via edge updates and WARPQ3D) then apply WARP3D
!
!        >  If specified, ensure orthogonality at the surface by adjusting the indicated number of points within the boundary layer
!
!  Control file format ('grid_morph_1.inp')
!
!     GRID_MORPH_1 controls for text case xxx
!     ------------ INPUT GRID ---------------
!     baseline.xyz            Input grid file name
!     T                       Formatted? [T|F]
!     4                       Number of blocks to be morphed
!     4 5 6 7                 Block numbers to morph
!     5 5 5 5                 Face numbers on the body (5 => k = 1, 6 => k = maximum k; 1 or 2 for i and 3 or 4 for j likewise)
!     ------------ PERTURBED SURFACE DEFINITION ------------
!     damage.xyz              Structured surface grid file name
!     T                       Formatted? [T|F]
!     ------------ MORPHED GRID --------------
!     morphed.xyz             Output grid file name
!     T                       Formatted? [T|F]
!     ------------ MISCELLANEOUS CONTROLS ------------------
!     3                       Surface point method [1 = foot of normal; 2 = line/surface intersection; 3 = 2 + B.L. adjustment]
!     1                       Volume update method [1 = 1-D update of radial lines; 3 = 3-space TFI-like WARP utilities]
!     30                      Point # approximating edge of boundary layer for surface point method = 3 option
!     3                       Symmetry plane definition:  1 = X, 2 = Y, 3 = Z [= 0 on symmetry plane]
!
!  Sponsor:
!
!     Dr. James Reuther, Reacting Flow Environments Branch, NASA Ames Research Center
!
!  History:
!
!     02/10/04  DAS  Initial design.
!     03/01/04   "   Miscellaneous controls form an integer array for easier extensions.
!     03/02/04   "   Orthogonality control in the boundary layer via MORPH_LINE_3D appears to work.
!     03/04/04   "   Permuting each block in place if necessary allows the key routine to assume the k = 1 plane is on the body.
!     03/08/04   "   Switched from packed grid block coordinates to the "grid_type" structure employed at NASA LaRC.
!
!  Author:
!
!     David Saunders, ELORET/NASA Ames Research Center, Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Data structure definition:
!  !!!!!!!!!!!!!!!!!!!!!!!!!!

   use grid_block_structure

   implicit none

!  Local constants:
!  !!!!!!!!!!!!!!!!

   integer, parameter ::   &
      lunctl = 1,          & ! Control file
      lunin  = 2,          & ! Input grid; avoid CFD_IO_PACKAGE for now; assume formatted I/O
      lunsrf = 3,          & ! Structured surface grid defining damage/repair
      lunout = 4,          & ! Output grid
      luncrt = 6             ! Screen diagnostics

   character, parameter :: &
      format * 11 = 'unformatted'

!  Local data structures:
!  !!!!!!!!!!!!!!!!!!!!!!

   type (grid_type) :: &
      surface_block

   type (grid_type), allocatable, dimension(:) :: &
      baseline_blocks

!  Local variables:
!  !!!!!!!!!!!!!!!!

   integer :: &
      i, i1, ib, ier, ios, n, nblocks_baseline, nblocks_to_morph, npts

   integer :: &
      misc_controls(20)      ! I.e., more than enough

   integer, allocatable, dimension (:) :: &
      iblocks_to_morph, num_face_on_body

   logical :: &
      formatted

   character :: &
      filename * 80

!  Execution:
!  !!!!!!!!!!

!  Modularize the I/O, or not?  Keep it simple for now by letting the main program handle it.

   open (lunctl, file='grid_morph_1.inp', status='old', iostat=ios)

   if (ios /= 0) then
      write (luncrt, '(/, a)') ' Unable to open grid_morph_1.inp.'
      go to 999
   end if

!  Baseline grid inputs:
!  ---------------------

   read (lunctl, *)                  ! Case description
   read (lunctl, *)                  ! Header
   read (lunctl, *) filename         ! Input grid name
   read (lunctl, *) formatted
   read (lunctl, *) nblocks_to_morph

   allocate (iblocks_to_morph(nblocks_to_morph), num_face_on_body(nblocks_to_morph))

   read (lunctl, *) iblocks_to_morph
   read (lunctl, *) num_face_on_body

   i1 = 1; if (formatted) i1 = 3

   open (lunin, file=filename, status='old', form=format(i1:11), iostat=ios)

   if (ios /= 0) then
      write (luncrt, '(/, 2a)') ' Unable to open input grid: ', filename(1:len_trim(filename))
      go to 999
   end if

   if (formatted) then
      read (lunin, *) nblocks_baseline
   else
      read (lunin)    nblocks_baseline
   end if

   allocate (baseline_blocks(nblocks_baseline))

   if (formatted) then
      read (lunin, *) (baseline_blocks(ib)%ni, baseline_blocks(ib)%nj, baseline_blocks(ib)%nk, ib = 1, nblocks_baseline)
   else
      read (lunin)    (baseline_blocks(ib)%ni, baseline_blocks(ib)%nj, baseline_blocks(ib)%nk, ib = 1, nblocks_baseline)
   end if

   do ib = 1, nblocks_baseline

      call allocate_block (baseline_blocks(ib))

      npts = baseline_blocks(ib)%ni * baseline_blocks(ib)%nj * baseline_blocks(ib)%nk

      call block_io (1, lunin, formatted, npts, baseline_blocks(ib)%x, baseline_blocks(ib)%y, baseline_blocks(ib)%z, ier)

      if (ier /= 0) then
         write (luncrt, '(/, a, i4, a, i4)') ' Trouble reading input grid.  Block #:', ib, '  ier:', ier
         go to 999
      end if

   end do

   close (lunin)

!  Perturbed surface definition:
!  -----------------------------

   read (lunctl, *)            ! Header
   read (lunctl, *) filename   ! Perturbed surface grid name
   read (lunctl, *) formatted
   i1 = 1; if (formatted) i1 = 3

   open (lunsrf, file=filename, status='old', form=format(i1:11), iostat=ios)

   if (ios /= 0) then
      write (luncrt, '(/, 2a)') ' Unable to open surface file: ', filename(1:len_trim(filename))
      go to 999
   end if

   if (formatted) then
      read (lunsrf, *) n       ! Should be 1 (multiblock form)
   else
      read (lunsrf)    n
   end if

   if (n /= 1) rewind (lunsrf) ! Try single-block form

   if (formatted) then
      read (lunsrf, *) surface_block%ni, surface_block%nj
   else
      read (lunsrf)    surface_block%ni, surface_block%nj
   end if

   surface_block%nk = 1

   call allocate_block (surface_block)

   npts = surface_block%ni * surface_block%nj

   call block_io (1, lunsrf, formatted, npts, surface_block%x, surface_block%y, surface_block%z, ier)

   if (ier /= 0) then
      write (luncrt, '(/, a, i4)') ' Trouble reading surface block.  ier:', ier
      go to 999
   end if

   close (lunsrf)

!  Morphed grid controls:
!  ----------------------

   read (lunctl, *)          ! Header
   read (lunctl, *) filename ! Output grid name
   read (lunctl, *) formatted
   i1 = 1; if (formatted) i1 = 3

   open (lunout, file=filename, status='unknown', form=format(i1:11), iostat=ios)

   if (ios /= 0) then
      write (luncrt, '(/, 2a)') ' Unable to open output grid: ', filename(1:len_trim(filename))
      go to 999
   end if

!  Miscellaneous controls:
!  -----------------------

   read (lunctl, *)          ! Header
   do i = 1, 999             ! Until EOF
      read (lunctl, *, iostat=ios) misc_controls(i)
      if (ios < 0) exit
      if (ios > 0) then
         write (luncrt, '(/, a, i2)') ' Trouble reading miscellaneous integer control #', i
         go to 999
      end if
   end do

   close (lunctl)

!  Perform the specialized type of morphing:
!  _________________________________________

   call morph_1 (nblocks_baseline, baseline_blocks, nblocks_to_morph, iblocks_to_morph, num_face_on_body, surface_block, &
                 misc_controls, ier)

   if (ier /= 0) go to 999

!  Save the morphed grid:
!  ----------------------

   if (formatted) then
      write (lunout, '(i3)') nblocks_baseline
   else
      write (lunout) nblocks_baseline
   end if

   if (formatted) then
      write (lunout, '(3i5)') (baseline_blocks(ib)%ni, baseline_blocks(ib)%nj, baseline_blocks(ib)%nk, ib = 1, nblocks_baseline)
   else
      write (lunout) (baseline_blocks(ib)%ni, baseline_blocks(ib)%nj, baseline_blocks(ib)%nk, ib = 1, nblocks_baseline)
   end if

   do ib = 1, nblocks_baseline

      npts = baseline_blocks(ib)%ni * baseline_blocks(ib)%nj * baseline_blocks(ib)%nk

      call block_io (2, lunout, formatted, npts, baseline_blocks(ib)%x, baseline_blocks(ib)%y, baseline_blocks(ib)%z, ier)

      if (ier /= 0) then
         write (luncrt, '(/, a, i4, a, i4)') ' Trouble saving morphed grid.  Block #:', ib, '  ier:', ier
         go to 999
      end if

   end do

   close (lunout)

   999 continue

   end program grid_morph_1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine allocate_block (block)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use grid_block_structure

   implicit none

!  Arguments:

   type (grid_type), intent(inout) :: block

!  Execution:

   allocate (block%x(block%ni, block%nj, block%nk), block%y(block%ni, block%nj, block%nk), block%z(block%ni, block%nj, block%nk))

   end subroutine allocate_block

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine block_io (mode, lun, formatted, npts, x, y, z, ier)
!
!  Read or write one 3-space PLOT3D-type grid block efficiently (all coordinates packed as one record of length 3 * npts).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in) ::                   mode       ! 1 = read, 2 = write
   integer, intent (in) ::                   lun        ! Logical unit number
   logical, intent (in) ::                   formatted  ! T|F
   integer, intent (in) ::                   npts       ! Number of (x,y,z) points in the block
   real, dimension (npts), intent (inout) :: x, y, z    ! Grid block coordinates
   integer, intent (out) ::                  ier        ! 0 = no error

!  Execution:

   select case (mode)
      case (1)
         if (formatted) then
            read (lun, *, iostat=ier) x, y, z
         else
            read (lun, iostat=ier) x, y, z
         end if
      case (2)
        if (formatted) then
            write (lun, '(1p, 4e19.11)', iostat=ier) x, y, z
         else
            write (lun, iostat=ier) x, y, z
         end if
      case default
         ier = 999
   end select

   end subroutine block_io

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine morph_1 (nblocks_baseline, baseline_blocks, nblocks_to_morph, iblocks_to_morph, num_face_on_body, surface_block, &
                       misc_controls, ier)
!
!  This routine implements the specialized type of grid morphing described in program grid_morph_1, q.v.
!  It is the argument-driven, reusable portion of that program.  Briefly, it morphs the specified blocks of a baseline grid
!  to match the moderately perturbed damage/repair defined by the given surface patch.
!
!  This version no longer assumes k = 1 is the geometry surface for all baseline blocks listed in iblocks_to_morph(*).
!
!  02/11/04  DAS  Initial implementation.
!  03/04/04   "   Permute a block if necessary to make k = 1 the face on the body, then permute back after morphing.
!
!  Author:  David Saunders, ELORET/NASA Ames Research Center, CA
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Data structure definition:
!  !!!!!!!!!!!!!!!!!!!!!!!!!!

   use grid_block_structure

   implicit none

!  Arguments:
!  !!!!!!!!!!

   integer, intent (in) :: &
      nblocks_baseline                      ! Number of blocks in baseline grid

   type (grid_type), intent (inout) :: &
      baseline_blocks(nblocks_baseline)     ! Baseline grid blocks in input; output with some blocks morphed

   integer, intent (in) :: &
      nblocks_to_morph,                   & ! Number of blocks to morph
      iblocks_to_morph(nblocks_to_morph), & ! List of blocks to morph
      num_face_on_body(nblocks_to_morph)    ! Face numbers 1-6 indicating body surface:  1 => i = 1 face, and so on

   type (grid_type), intent (in) :: &
      surface_block                         ! Structured surface patch defining damage/repair

   integer, intent (in) :: &
      misc_controls(*)                      ! Miscellaneous controls, renamed and documented in routine morph_block below

   integer, intent (out) :: &
      ier                                   ! 0 = no problem encountered, else fatal error

!  Local variables:
!  !!!!!!!!!!!!!!!!

   integer :: &
      ib, iblock, ix, mi, mj, mk, ni, nj, nk

!  Execution:
!  !!!!!!!!!!

!  For each block specified in the densified region of the baseline grid ...

   do iblock = 1, nblocks_to_morph

      ib = iblocks_to_morph(iblock)
      mi = baseline_blocks(ib)%ni
      mj = baseline_blocks(ib)%nj
      mk = baseline_blocks(ib)%nk

!     Push processing of one block down a level to simplify the mnemonics.

!     Ensure k = 1 is the face on the body, at least temporarily:

      call permute_block (1, mi, mj, mk, num_face_on_body(iblock), &
                          baseline_blocks(ib)%x, baseline_blocks(ib)%y, baseline_blocks(ib)%z, ni, nj, nk)

!     Perform the procedures common to all affected blocks:

      call morph_block (ni, nj, nk, baseline_blocks(ib)%x, baseline_blocks(ib)%y, baseline_blocks(ib)%z, &
                        surface_block%ni, surface_block%nj, surface_block%x, surface_block%y, surface_block%z, &
                        misc_controls, ier)
      if (ier /= 0) exit

!     Ensure the original (i,j,k) indexing is recovered:

      call permute_block (2, ni, nj, nk, num_face_on_body(iblock), &
                          baseline_blocks(ib)%x, baseline_blocks(ib)%y, baseline_blocks(ib)%z, mi, mj, mk)

   end do ! Next block to morph

   end subroutine morph_1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine morph_block (ni, nj, nk, x, y, z, ni_srf, nj_srf, surface_x, surface_y, surface_z, misc_controls, ier)

!  Morph one block by the method outlined in program morph_grid_1, q.v.
!  Briefly, we replace the geometry face of the given grid block with one interpolated to the given damaged/repaired surface grid,
!  and adjust the rest of the grid block accordingly, in-place.
!
!  At this level, it is assumed that k = 1 is the geometry surface for the block being morphed.
!
!  02/12/04  DAS  Initial implementation (1-D options involving radial lines).
!  02/17/04   "   Option to intersect radial lines with new surface.
!  02/18/04   "   Option to morph the block faces and interior via the 3-space methods of WARPQ3D and WARP3D.
!  03/01/04   "   Handle symmetry-plane radial lines more carefully - projections can introduce gaps.
!                 Option to control orthogonality within the boundary layer.
!  03/16/04   "   Had to handle the ier = 1 case for the initial calls to PROJECT4 and distinguish between p & q close to [0,1] vs.
!                 being converged but well outside the unit square.
!  03/18/04   "   Automatic arrays appeared to give some compilers trouble - avoid them where internal procedures inherit them.
!
!  Author:  David Saunders, ELORET/NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:
!  !!!!!!!!!!

   integer, intent (in) :: &
      ni, nj, nk                              ! Dimensions of the current block

   real, dimension (ni,nj,nk), intent (inout) :: &
      x, y, z                                 ! Input with baseline grid block; output with morphed grid block

   integer, intent (in) :: &
      ni_srf, nj_srf                          ! Dimensions of the damaged/repaired surface grid ...

   real, intent (in), dimension (ni_srf,nj_srf) :: &
      surface_x, surface_y, surface_z         ! ... and its coordinates

   integer, intent (in) :: &
      misc_controls(*)                        ! Miscellaneous controls, renamed and documented below

   integer, intent (out) :: &
      ier                                     ! 0 = no problem encountered, else fatal error

!  Local constants:
!  !!!!!!!!!!!!!!!!

   integer, parameter :: &
      len_line = 30,     & ! Reasonable limit on the number of radial line points treated during line/intersection calculations
      lunout = -6          ! INTSEC5 diagnostics go to unit 6

   real, parameter :: &
      big = 1.e+32, half = 0.5, one = 1., zero = 0.

   logical, parameter :: &
      false = .FALSE.

   character, parameter :: &
      method_line * 1 = 'L' ! Linear interpolation suffices here for line/surface intersections

!  Local variables:
!  !!!!!!!!!!!!!!!!

   integer :: &
      i, i1, i2, ibest, isrf, isym, j, j1, j2, jbest, jsrf, k, k_line, klimit, method_surface, method_volume, &
      nblayer, nblhalf

   integer, allocatable, dimension (:,:) :: &
      i_best, j_best

   real :: &
      data_range_sq, dp, dpq, dpqmin, dq, dsq, dsqmin, eps, p0, p0best, pc, q0, q0best, qc, &
      tint, tol, tol_pq, ttotal, uint, vint, xint, yint, zint

   real, dimension (3) :: &
      x0, x0best, x1, x2, x3, x4, xi, xt

   real, dimension (len_line) :: &
      line_t, line_x, line_y, line_z

   real, allocatable, dimension (:) :: &
      xold, yold, zold, xnew, ynew, znew

   real, allocatable, dimension (:,:) :: &
      surface_u, surface_v

   real, allocatable, dimension (:,:,:,:) :: &
      morph_xyz, s

   real, allocatable, dimension (:,:,:,:,:) :: &
      dfacei, dfacej, dfacek

   logical :: &
      intersect, on_symmetry_plane

   logical, allocatable, dimension (:,:) :: &
      new_surface_point

!  Execution:
!  !!!!!!!!!!

!  Transcribe the miscellaneous controls to mnemonic control variables:

   method_surface = misc_controls(1)    ! Surface point control:
                                        !    1 = foot of projection from original surface to new surface
                                        !    2 = intersection of radial line with new surface (using foot as starting guess)
                                        !    3 = 2 + slope control within the boundary layer; see nblayer
   method_volume  = misc_controls(2)    ! Volume point control:
                                        !    1 = 1-D method (independent update of radial lines)
                                        !    3 = 3-D method (3-space TFI-like methods of WARPQ3D and WARP3D)
   nblayer        = misc_controls(3)    ! Point number defining approximate edge of boundary layer; e.g., 30;
                                        !    used if method_surface = 3
   isym           = misc_controls(4)    ! Symmetry plane definition:
                                        !    1 = X, 2 = Y, 3 = Z [= zero on symmetry plane]

!  Determine a measure for when a grid point is clearly well away from the new surface:

   dsq =           (surface_x(1,1) - surface_x(ni_srf,nj_srf))**2 + (surface_y(1,1) - surface_y(ni_srf,nj_srf))**2 + &
                   (surface_z(1,1) - surface_z(ni_srf,nj_srf))**2
   data_range_sq = (surface_x(ni_srf,1) - surface_x(1,nj_srf))**2 + (surface_y(ni_srf,1) - surface_y(1,nj_srf))**2 + &
                   (surface_z(ni_srf,1) - surface_z(1,nj_srf))**2
   data_range_sq = max (dsq, data_range_sq)

!  Many of the morphed block points many not be affected, so initialize them all to begin with:

   allocate (morph_xyz(ni,nj,nk,3), new_surface_point(ni,nj))

   new_surface_point = .FALSE. ! 1:ni, 1:nj

   do k = 1, nk
      do j = 1, nj
         do i = 1, ni
            morph_xyz(i,j,k,1) = x(i,j,k)
            morph_xyz(i,j,k,2) = y(i,j,k)
            morph_xyz(i,j,k,3) = z(i,j,k)
         end do
      end do
   end do

!  Always set up for surface/line intersections.  Symmetry plane radial lines require them even if feet of normals are specified.
!  There is minor duplication here if multiple blocks are being processed at the higher level.

   allocate (i_best(ni,nj), j_best(ni,nj))  ! For the orthogonality option (adjust_slopes)

   allocate (xold(nk), yold(nk), zold(nk), xnew(nk), ynew(nk), znew(nk)) ! Automatic arrays gave trouble

   allocate (surface_u(ni_srf,nj_srf), surface_v(ni_srf,nj_srf)) ! Keep it local to shorten the argument list

   surface_u(1,1) = -1.                   ! Tells INTSEC5 to parameterize the new surface on its first call here (only)
   klimit = min (len_line, nk)            ! No need to work with entire radial lines
   eps = max (5.* epsilon (eps), 1.e-7)   ! For PBILINT/RIPPLE2D used by INTSEC5
   tol = eps * sqrt (data_range_sq)       ! For symmetry plane test
!!!tol_pq = 10.* sqrt (epsilon (eps))     ! For PROJECT4: p & q cannot be closer than within sqrt (machine eps) of correct minimum
   tol_pq = 1.E-2                         ! Why does this need to be this big???

!  For each point on the baseline surface ...

   do j = 1, nj

      do i = 1, ni

         xt(1) = x(i,j,1) ! Current target point
         xt(2) = y(i,j,1)
         xt(3) = z(i,j,1)

!        Locate the nearest new surface POINT, then translate that to the nearest CELL later:

         dsqmin = big
         do jsrf = 1, nj_srf
            do isrf = 1, ni_srf
               dsq = (xt(1) - surface_x(isrf,jsrf))**2 + (xt(2) - surface_y(isrf,jsrf))**2 + (xt(3) - surface_z(isrf,jsrf))**2
               if (dsq < dsqmin) then
                  dsqmin = dsq
                  ibest  = isrf
                  jbest  = jsrf
               end if
            end do
         end do

         if (dsqmin > data_range_sq) cycle ! The original point is OK and has already been inserted.

!        The nearest CELL of the new surface should be one of at most 4 cells defined by (ibest,jbest).
!        The target may lie outside this surface.

         i1 = max (1, ibest - 1);  i2 = min (ni_srf - 1, ibest)
         j1 = max (1, jbest - 1);  j2 = min (nj_srf - 1, jbest)
         dpqmin = big

         do jsrf = j1, j2

            do isrf = i1, i2

!              Set up the quad. cell defined by (isrf,jsrf) as vertices in counterclockwise order:

               x1(1) = surface_x(isrf,  jsrf  );  x1(2) = surface_y(isrf,  jsrf  );  x1(3) = surface_z(isrf,  jsrf  )
               x2(1) = surface_x(isrf+1,jsrf  );  x2(2) = surface_y(isrf+1,jsrf  );  x2(3) = surface_z(isrf+1,jsrf  )
               x3(1) = surface_x(isrf+1,jsrf+1);  x3(2) = surface_y(isrf+1,jsrf+1);  x3(3) = surface_z(isrf+1,jsrf+1)
               x4(1) = surface_x(isrf,  jsrf+1);  x4(2) = surface_y(isrf,  jsrf+1);  x4(3) = surface_z(isrf,  jsrf+1)

               call project4 (x1, x2, x3, x4, xt, x0, dsq, p0, q0, ier) ! Find foot of normal to ["plane" of] quad. cell

               if (ier == 0) then ! Converged with p0 and q0 both in [0, 1]; this is the cell we're looking for
                  morph_xyz(i,j,1,1:3) = x0(1:3)
                  exit
               else if (ier == 1) then ! Converged, but p0 and q0 aren't in [0, 1]; save the closest such cell
                  dp = abs (p0 - half);  dq = abs (q0 - half)  ! One of these must exceed one half
                  dpq = max (dp, dq)
                  if (dpq < dpqmin) then
                     dpqmin = dpq
                     ibest = isrf;  jbest = jsrf
                     x0best = x0
                  end if
               else
                  write (-lunout, '(a, 3i5)') ' morph_block trouble.  project4 ier, isrf, jsrf:', ier, isrf, jsrf
               end if

            end do

            if (ier == 0) exit

         end do

         if (ier > 0) then
            if (dpqmin - half < tol_pq) then ! p, q were very close to [0, 1]; how do we distinguish being slightly off new surface?
               isrf = ibest
               jsrf = jbest
               morph_xyz(i,j,1,1:3) = x0best(1:3)
            else ! Force leaving the target point alone
               ier = 2
            end if
!!!!!       write (-lunout, '(a, 2i5, 1pe15.6)') ' Possible projection trouble.  isrf, jsrf, dpqmin:', isrf, jsrf, dpqmin
            if (ier > 1) cycle ! Presumably the original surface point must be outside the new surface, and is already in place.
         end if

         new_surface_point(i,j) = .TRUE.
         i_best(i,j) = isrf
         j_best(i,j) = jsrf

!        Intersect the original radial line with the new surface?
!        --------------------------------------------------------

         on_symmetry_plane = abs (xt(isym)) < tol
         intersect = method_surface > 1 .or. on_symmetry_plane

         if (intersect) then

            do k = 1, klimit ! Transfer the relevant part of the base radial line
               line_x(k) = x(i,j,k)
               line_y(k) = y(i,j,k)
               line_z(k) = z(i,j,k)
            end do

            k_line = 1       ! Estimate of radial line interval containing the intersection
            ttotal = zero    ! Tells INTSEC5 to parameterize the line
            tint   = zero    ! Tells INTSEC5 to use line_t(k_line) for the starting guess
            uint   = zero    ! Use surface_u/v(isrf,jsrf) for the (u,v) starting guess
            vint   = zero

            call intsec5 (ni_srf, nj_srf, 1, ni_srf, 1, nj_srf, surface_x, surface_y, surface_z, surface_u, surface_v, &
                    isrf, jsrf, eps, 1, klimit, line_x, line_y, line_z, line_t, k_line, ttotal, method_line, false, &
                    tint, uint, vint, xint, yint, zint, lunout, ier)

            if (ier <= 2) then
               morph_xyz(i,j,1,1) = xint
               morph_xyz(i,j,1,2) = yint
               morph_xyz(i,j,1,3) = zint
            else
               write (-lunout, '(/, a, i3, a, /, 1p, 3e16.7)') ' INTSEC5 ier:', ier, '.  Current baseline surface point:', xt
            end if ! Continue with the foot of the normal rather than quit

            ier = 0 ! Avoid quitting at the higher level

         end if

!        We can't adjust the surface point for orthogonality until the interim volume points have been generated.

      end do ! Next i for this j of baseline block at nominal surface

   end do ! Next j for this baseline block at nominal surface

   deallocate (surface_u, surface_v)


!  Update the block interior:
!  --------------------------

   select case (method_volume)

   case (1) ! Treat each affected radial line independently 

      do j = 1, nj
         do i = 1, ni
            if (new_surface_point(i,j)) call update_radial_line  ! Local procedure below updates the block interior in place
         end do
      end do

   case (3) ! Update the entire volume based on the new faces and the original volume distribution

!     First, calculate new radial block edge lines by the 1-D method:

      do j = 1, nj, nj - 1
         do i = 1, ni, ni - 1
            if (new_surface_point(i,j)) call new_radial_line  ! Local procedure below updates the block edge in the copy
         end do
      end do

!     Calculate new block faces by the quasi-3D method of WARPQ3D.
!     First, we need the normalized arc lengths for the original block.

      allocate (s(ni,nj,nk,3), dfacei(3,nj,nk,2,4), dfacej(3,ni,nk,2,4), dfacek(3,ni,nj,2,4))

      call paramxyz (1, ni, 1, nj, 1, nk,  1, ni,  1, nj, 1, nk, x, y, z, s)

      call warpq3d  (1, ni, 1, nj, 1, nk,  1,  1,  1, nj, 1, nk, x, y, z, s, dfacei, dfacej, dfacek, &         ! i = 1  face
                     morph_xyz(1,1,1,1), morph_xyz(1,1,1,2), morph_xyz(1,1,1,3))

      call warpq3d  (1, ni, 1, nj, 1, nk, ni, ni,  1, nj, 1, nk, x, y, z, s, dfacei, dfacej, dfacek, &         ! i = ni face
                     morph_xyz(1,1,1,1), morph_xyz(1,1,1,2), morph_xyz(1,1,1,3))

      call warpq3d  (1, ni, 1, nj, 1, nk,  1, ni,  1,  1, 1, nk, x, y, z, s, dfacei, dfacej, dfacek, &         ! j = 1  face
                     morph_xyz(1,1,1,1), morph_xyz(1,1,1,2), morph_xyz(1,1,1,3))

      call warpq3d  (1, ni, 1, nj, 1, nk,  1, ni, nj, nj, 1, nk, x, y, z, s, dfacei, dfacej, dfacek, &         ! j = nj face
                     morph_xyz(1,1,1,1), morph_xyz(1,1,1,2), morph_xyz(1,1,1,3))

!     Finally, derive the new block interior from the new faces and the original interior spacings:

      call warp3d   (1, ni, 1, nj, 1, nk,  1, ni,  1, nj, 1, nk, x, y, z, s, dfacei, dfacej, dfacek, &
                     morph_xyz(1,1,1,1), morph_xyz(1,1,1,2), morph_xyz(1,1,1,3))
 
      deallocate (s, dfacei, dfacej, dfacek)

!     Replace the original block with the morphed block:

      do k = 1, nk
         do j = 1, nj
            do i = 1, ni
               x(i,j,k) = morph_xyz(i,j,k,1)
               y(i,j,k) = morph_xyz(i,j,k,2)
               z(i,j,k) = morph_xyz(i,j,k,3)
            end do
         end do
      end do

   case default

      write (*, '(/, a, i10)') ' Bad input for volume update method:', method_volume
      ier = 999
      go to 999

   end select ! Choice of volume grid methods


!  Impose orthogonality in the boundary layer?
!  -------------------------------------------

   if (method_surface == 3) call adjust_slopes ! Internal procedure below


!  Done with this block.
!  ---------------------

   deallocate (morph_xyz, new_surface_point, i_best, j_best, xold, yold, zold, xnew, ynew, znew)

   ier = 0

   999 continue

!  Local procedures for subroutine morph_block:
!  --------------------------------------------

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine adjust_slopes

!     Adjust radial lines for orthogonality in the boundary layer.
!
!     03/01/04  DAS  Initial implementation.
!     03/17/04   "   Wing leading edge wrap case revealed pitfalls of normal projections to the wrong wing surface.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local variables:

      integer ier_best, k
      real    dsqtol, dx, dy, dz, p, q, r, pm1, qm1, s1(3), s2(3), xfoot(3)

!     Execution:

!     Determine a representative index for the middle of the boundary layer:

      do k = 1, nblayer
         xold(k) = x(1,1,k)
         yold(k) = y(1,1,k)
         zold(k) = z(1,1,k)
      end do

      call chords3d (nblayer, xold, yold, zold, false, r, xnew)  ! Reuse r for total arc, xnew(*) for cumulative arcs

      dsqtol = r * r  ! Plausible limit on squared length of orthogonal projections, to avoid accepting a cell too far away
      r = r * 0.5
      nblhalf = (nblayer * 3) / 4

      call interval (nblayer, xnew, r, one, nblhalf) ! 1-D search utility

      do j = 1, nj

         do i = 1, ni

            if (.not. new_surface_point(i,j)) cycle

            xt(1) = x(i,j,nblhalf) ! Middle of boundary layer, roughly
            xt(2) = y(i,j,nblhalf)
            xt(3) = z(i,j,nblhalf)

            on_symmetry_plane = abs (xt(isym)) < tol

!           The original nearest surface cell will often be close to the cell associated with the new end of the radial line.

            i1 = max (1, i_best(i,j) - 2);  i2 = min (ni_srf - 1, i_best(i,j) + 2)
            j1 = max (1, j_best(i,j) - 2);  j2 = min (nj_srf - 1, j_best(i,j) + 2)

            do jsrf = j1, j2

               do isrf = i1, i2

!                 Set up the quad. cell defined by (isrf,jsrf) as vertices in counterclockwise order:

                  x1(1) = surface_x(isrf,  jsrf  );  x1(2) = surface_y(isrf,  jsrf  );  x1(3) = surface_z(isrf,  jsrf  )
                  x2(1) = surface_x(isrf+1,jsrf  );  x2(2) = surface_y(isrf+1,jsrf  );  x2(3) = surface_z(isrf+1,jsrf  )
                  x3(1) = surface_x(isrf+1,jsrf+1);  x3(2) = surface_y(isrf+1,jsrf+1);  x3(3) = surface_z(isrf+1,jsrf+1)
                  x4(1) = surface_x(isrf,  jsrf+1);  x4(2) = surface_y(isrf,  jsrf+1);  x4(3) = surface_z(isrf,  jsrf+1)

                  call project4 (x1, x2, x3, x4, xt, x0, dsq, p0, q0, ier)

                  if (ier == 0 .and. dsq < dsqtol) exit ! p0 and q0 are both in [0, 1]; use the foot of the nearby normal.

               end do

               if (ier == 0 .and. dsq < dsqtol) exit

            end do

            if (ier /= 0 .or. dsq >= dsqtol) then ! Widen the above search in case we missed

               dsqmin = big

               do jsrf = 1, nj_srf - 1

                  do isrf = 1, ni_srf - 1

!                    Quad. cell defined by (isrf,jsrf), in counterclockwise order:

                     x1(1) = surface_x(isrf,  jsrf  );  x1(2) = surface_y(isrf,  jsrf  );  x1(3) = surface_z(isrf,  jsrf  )
                     x2(1) = surface_x(isrf+1,jsrf  );  x2(2) = surface_y(isrf+1,jsrf  );  x2(3) = surface_z(isrf+1,jsrf  )
                     x3(1) = surface_x(isrf+1,jsrf+1);  x3(2) = surface_y(isrf+1,jsrf+1);  x3(3) = surface_z(isrf+1,jsrf+1)
                     x4(1) = surface_x(isrf,  jsrf+1);  x4(2) = surface_y(isrf,  jsrf+1);  x4(3) = surface_z(isrf,  jsrf+1)

                     call project4 (x1, x2, x3, x4, xt, x0, dsq, p0, q0, ier)

                     if (ier == 0 .and. dsq < dsqtol) exit ! p0 and q0 are both in [0, 1]; use the foot of the nearby normal.

!                    Use of dsq is not right - a distant cell can appear close but (p,q) may be well outside the unit square

                     p   = max (zero, min (one, p0));  pm1 = one - p  ! Evaluate (x,y,z) inside the cell or at an edge, not outside
                     q   = max (zero, min (one, q0));  qm1 = one - q

                     dx  = pm1 * (qm1 * x1(1) + q * x4(1)) + p * (qm1 * x2(1) + q * x3(1)) - xt(1)
                     dy  = pm1 * (qm1 * x1(2) + q * x4(2)) + p * (qm1 * x2(2) + q * x3(2)) - xt(2)
                     dz  = pm1 * (qm1 * x1(3) + q * x4(3)) + p * (qm1 * x2(3) + q * x3(3)) - xt(3)

                     dsq = dx**2 + dy**2 + dz**2

                     if (dsq < dsqmin) then ! Save the nearest adjusted projection for occasional cases, incl. on symmetry plane
                        dsqmin = dsq
                        xfoot  = x0  ! 1:3
                        ier_best    = ier
                        i_best(i,j) = isrf
                        j_best(i,j) = jsrf
                     end if

                  end do

                  if (ier == 0 .and. dsq < dsqtol) exit

               end do

               if (ier /= 0 .or. dsq >= dsqtol) then ! Pick the nearest cell

                  x0 = xfoot ! 1:3

                  if (ier_best > 1) then ! Convergence problem
                     write (-lunout, '(/, a, 2i4, i6)') ' Trouble projecting from boundary layer middle.  i,j, ier: ', i, j, ier
                     write (-lunout, '(a, 1p, 3e17.7)') ' Target point: ', xt
                     write (-lunout, '(a, 1p, 3e17.7)') ' Nearest foot: ', xfoot

                     isrf = i_best(i,j);   jsrf = j_best(i,j)
                     write (-lunout, '(a, 2i4)') ' Using surface patch cell: ', isrf, jsrf
                     write (-lunout, '(1p, 3e17.7, 3x, 3e17.7)') &
                        surface_x(isrf,  jsrf  ), surface_y(isrf,  jsrf  ), surface_z(isrf,  jsrf  ), &
                        surface_x(isrf+1,jsrf  ), surface_y(isrf+1,jsrf  ), surface_z(isrf+1,jsrf  ), &
                        surface_x(isrf+1,jsrf+1), surface_y(isrf+1,jsrf+1), surface_z(isrf+1,jsrf+1), &
                        surface_x(isrf,  jsrf+1), surface_y(isrf,  jsrf+1), surface_z(isrf,  jsrf+1)
                  end if

                  if (on_symmetry_plane) x0(isym) = zero ! Staying on the symmetry plane should lose only a little orthogonality

               end if

            end if

!           Redistribute the bulk of the boundary layer points:

            xnew(1) = x0(1)
            ynew(1) = x0(2)
            znew(1) = x0(3)

            do k = 1, nblayer
               xold(k) = x(i,j,k)
               yold(k) = y(i,j,k)
               zold(k) = z(i,j,k)
            end do

            xnew(nblayer) = xold(nblayer)
            ynew(nblayer) = yold(nblayer)
            znew(nblayer) = zold(nblayer)

            dx  = xt(1) - x0(1) ! Set up the unit normal to the surface ...
            dy  = xt(2) - x0(2)
            dz  = xt(3) - x0(3)
            dsq = dx**2 + dy**2 + dz**2
            r   = one / sqrt (dsq)
            s1(1) = r * dx
            s1(2) = r * dy
            s1(3) = r * dz

            dx  = xold(nblayer) - xold(nblayer - 1) ! ... and the unit radial vector at the (estimated) edge of the boundary layer
            dy  = yold(nblayer) - yold(nblayer - 1)
            dz  = zold(nblayer) - zold(nblayer - 1)
            dsq = dx**2 + dy**2 + dz**2
            r   = one / sqrt (dsq)
            s2(1) = r * dx
            s2(2) = r * dy
            s2(3) = r * dz

!           Enforce end-point slopes, with smooth variation in between:

            call morph_line_3D (1, nblayer, s1, s2, xold, yold, zold, xnew, ynew, znew)

            do k = 1, nblayer - 1
               x(i,j,k) = xnew(k)
               y(i,j,k) = ynew(k)
               z(i,j,k) = znew(k)
            end do

         end do ! Next surface grid i

      end do ! Next surface grid j

      end subroutine adjust_slopes

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine new_radial_line

!     Morph a radial grid line to match a new surface end point at (i,j,1) via the 1-D arc-length-based method of NULINE3D.
!     This version leaves the original block intact for the rest of the 3-D method; update_radial_line can update it directly.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local variables:

      integer k

!     Execution:

      do k = 1, nk
         xold(k) = x(i,j,k)
         yold(k) = y(i,j,k)
         zold(k) = z(i,j,k)
      end do

      xnew(nk) = xold(nk)
      ynew(nk) = yold(nk)
      znew(nk) = zold(nk)

      xnew(1) = morph_xyz(i,j,1,1)
      ynew(1) = morph_xyz(i,j,1,2)
      znew(1) = morph_xyz(i,j,1,3)

      call nuline3d (1, nk, xold, yold, zold, xnew, ynew, znew) ! No slope control

      do k = 2, nk - 1
         morph_xyz(i,j,k,1) = xnew(k)
         morph_xyz(i,j,k,2) = ynew(k)
         morph_xyz(i,j,k,3) = znew(k)
      end do

      end subroutine new_radial_line

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine update_radial_line

!     Update the interior points of the radial line through (i,j,1) via the 1-D arc-length-based method of NULINE3D.
!     This version updates the original block directly for the case of applying the 1-D method to the entire volume,
!     whereas new_radial_line updates the copy of the block for morphing by the 3-D method before updating the original block.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local variables:

      integer k

!     Execution:

      do k = 1, nk
         xold(k) = x(i,j,k)
         yold(k) = y(i,j,k)
         zold(k) = z(i,j,k)
      end do

      xnew(nk) = xold(nk)
      ynew(nk) = yold(nk)
      znew(nk) = zold(nk)

      xnew(1) = morph_xyz(i,j,1,1)
      ynew(1) = morph_xyz(i,j,1,2)
      znew(1) = morph_xyz(i,j,1,3)

      call nuline3d (1, nk, xold, yold, zold, xnew, ynew, znew)

      do k = 1, nk - 1
         x(i,j,k) = xnew(k)
         y(i,j,k) = ynew(k)
         z(i,j,k) = znew(k)
      end do

!     No need to update the (i,j,*) line in morph_xyz() as well.

      end subroutine update_radial_line

   end subroutine morph_block

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine permute_block (mode, mi, mj, mk, num_face_on_body, x, y, z, ni, nj, nk)
!
!  Limited form of index permutation for a grid block to make k = 1 the indicated face, and to reverse the permutation.
!  The permuted block is returned in place.
!
!  03/03/04  David Saunders  Initial implementation (keeping it straightforward).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in) :: &
      mode,                &  ! 1 = initial permutation; 2 = reverse permutation
      mi, mj, mk,          &  ! Input block dimensions
      num_face_on_body        ! Input face number to change to k = 1 if mode = 1, or to change from k = 1 if mode = 2;
                              ! 1 means the i = 1 plane, 2 means the i = i(max) plane, j means the j = 1 plane, and so on

   real, intent (inout), dimension (mi,mj,mk) :: &  ! Packed grid block coordinates, input in the order implied by mi, mj, mk;
      x, y, z                                       ! output with permuted indices in the order implied by ni, nj, nk

   integer, intent (out) :: &
      ni, nj, nk              ! Dimensions of permuted block; when mode = 2, ni/j/k swap roles with mi/j/k in the call

!  Local integers:

   integer &
      i, j, k, k1, kk, knc, npts

   real, allocatable, dimension (:,:,:) :: &
      xt, yt, zt

!  Execution:

   npts = mi * mj * mk

   if (mode == 1) then ! Permute indices to achieve the desired face at k = 1

      select case (num_face_on_body)

      case (1, 2) ! i planes become k planes, possibly reversed

         ni = mj;  nj = mk;  nk = mi

         if (num_face_on_body == 1) then
            k1 = 1;  knc = 1
         else
            k1 = nk; knc = -1
         end if

         allocate (xt(ni,nj,nk), yt(ni,nj,nk), zt(ni,nj,nk))

         kk = k1
         do k = 1, nk
            do j = 1, nj
               do i = 1, ni
                  xt(i,j,k) = x(kk,i,j);  yt(i,j,k) = y(kk,i,j);  zt(i,j,k) = z(kk,i,j)
               end do
            end do
            kk = kk + knc
         end do

      case (3, 4) ! j planes become k planes, possibly reversed

         ni = mk;  nj = mi;  nk = mj

         if (num_face_on_body == 3) then
            k1 = 1;  knc = 1
         else
            k1 = nk; knc = -1
         end if

         allocate (xt(ni,nj,nk), yt(ni,nj,nk), zt(ni,nj,nk))

         kk = k1
         do k = 1, nk
            do j = 1, nj
               do i = 1, ni
                  xt(i,j,k) = x(j,kk,i);  yt(i,j,k) = y(j,kk,i);  zt(i,j,k) = z(j,kk,i)
               end do
            end do
            kk = kk + knc
         end do

      case (5, 6) ! k planes remain k planes, possibly reversed

         ni = mi;  nj = mj;  nk = mk

         if (num_face_on_body == 5) then
!           Do nothing
         else
            k1 = mk; knc = -1

            allocate (xt(ni,nj,nk), yt(ni,nj,nk), zt(ni,nj,nk))

            kk = k1
            do k = 1, nk
               do j = 1, nj
                  do i = 1, ni
                     xt(i,j,k) = x(i,j,kk);  yt(i,j,k) = y(i,j,kk);  zt(i,j,k) = z(i,j,kk)
                  end do
               end do
               kk = kk + knc
            end do
         end if

      end select ! The overwrite in-place is common to both modes

   else ! Reverse mode = 2

     select case (num_face_on_body)

      case (1, 2) ! k planes become i planes

         ni = mk;  nj = mi;  nk = mj

         if (num_face_on_body == 1) then
            k1 = 1;  knc = 1
         else
            k1 = mk; knc = -1
         end if

         allocate (xt(ni,nj,nk), yt(ni,nj,nk), zt(ni,nj,nk))

         kk = k1
         do k = 1, mk
            do j = 1, mj
               do i = 1, mi
                  xt(kk,i,j) = x(i,j,k);  yt(kk,i,j) = y(i,j,k);  zt(kk,i,j) = z(i,j,k)
               end do
            end do
            kk = kk + knc
         end do

      case (3, 4) ! k planes become j planes

         ni = mj;  nj = mk;  nk = mi

         if (num_face_on_body == 3) then
            k1 = 1;  knc = 1
         else
            k1 = mk; knc = -1
         end if

         allocate (xt(ni,nj,nk), yt(ni,nj,nk), zt(ni,nj,nk))

         kk = k1
         do k = 1, mk
            do j = 1, mj
               do i = 1, mi
                  xt(j,kk,i) = x(i,j,k);  yt(j,kk,i) = y(i,j,k);  zt(j,kk,i) = z(i,j,k)
               end do
            end do
            kk = kk + knc
         end do

      case (5, 6) ! k planes remain k planes

         ni = mi;  nj = mj;  nk = mk

         if (num_face_on_body == 5) then
!           Do nothing
         else
            k1 = mk; knc = -1

            allocate (xt(ni,nj,nk), yt(ni,nj,nk), zt(ni,nj,nk))

            kk = k1
            do k = 1, mk
               do j = 1, mj
                  do i = 1, mi
                     xt(i,j,kk) = x(i,j,k);  yt(i,j,kk) = y(i,j,k);  zt(i,j,kk) = z(i,j,k)
                  end do
               end do
               kk = kk + knc
            end do
         end if

      end select ! The overwrite in-place is common to both modes

   end if

   if (num_face_on_body /= 5) then

      call transfer (npts, xt, yt, zt, x, y, z)  ! Local procedure below overwrites input block as packed vectors

      deallocate (xt, yt, zt)

   end if

!  Local procedure for subroutine permute_block:

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine transfer (npts, xt, yt, zt, x, y, z)
!
!     Copy a packed grid block.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in) :: npts

      real, intent (in),  dimension (npts) :: xt, yt, zt

      real, intent (out), dimension (npts) :: x,  y,  z

!     Local variable:

      integer i

!     Execution:

      do i = 1, npts
         x(i) = xt(i);  y(i) = yt(i);  z(i) = zt(i)
      end do

      end subroutine transfer

   end subroutine permute_block
