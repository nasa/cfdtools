!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program surface_patches
!
!  Description:
!
!     This is an adaptation of ADJUST_GRID for the case of surface datasets,
!     which lend themselves to further common requirements, such as extracting
!     subpatches.  PLOT3D multiblock grid files are treated, with automatic
!     detection of formatted/unformatted inputs.
!
!     SURFACE_PATCHES applies one or more of the ADJUST_GRID transformations
!     (shift, scale, rotate, etc.) to all blocks of the indicated surface grid.
!     It has additional options applied to one patch at a time, employing the
!     utilities provided by the earlier surface_patch_utilities.f90 module.
!
!     Some of these options may make sense for just one operation per run.
!     In particular, extracting a portion of a patch and reversing the indexing
!     is the requirement that first prompted a generalized driver for these
!     surface patch utilities.
!
!  History:
!
!     02/18/00  D.A.Saunders  Initial SCALE_GRID adaptation of RESHAPE3D, which
!                             operates on 3-column formatted datasets.
!     01/20/08    "      "    ADJUST_GRID adapted from SCALE_GRID to reduce the
!                             number of prompts.
!     05/06/09    "      "    SURFACE_PATCHES adapted from ADJUST_GRID as a
!                             driver for surface_patch_utilities.f90 and other
!                             likely operations on structured surfaces.
!     02/12/15    "      "    Added a "fissure" option to simulate surface
!                             cracks for asteroid studies.  See option 25.
!     02/13/15    "      "    Added an option for (index) diagonal fissures.
!     02/15/15    "      "    Added an option to edit point coordinates.
!                             This allows correcting the result of normals
!                             off symmetry plane points that aren't truly
!                             in the symmetry plane.
!     05/08/15    "      "    Added an option to apply fissures along all
!                             patch boundaries.
!     06/23/15    "      "    The option to apply fissures along patch edges
!                             serves to test new utility SURFACE_VECTORS, which
!                             attempts to produce identical results at common
!                             edge points.
!     11/30/15    "      "    Added the option for match patch edges following
!                             drawn-out extensions to surface_patch_utilities.
!     12/01/15    "      "    1.e-5 seems too tight for edge/corner matching.
!
!  Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA
!           Now with ERC, Inc. at NASA ARC (5 years).
!           Now with AMA, Inc. at NASA ARC.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_block_structure     ! Derived data type for one grid block
   use xyzq_io_module           ! PLOT3D multiblock file I/O package (3-space)
   use surface_patch_utilities  ! Utilities for structured surface grids

!  Constants:

   implicit none

   integer, parameter :: &
      lunin  = 1,        &
      lunout = 2,        &
      lunkbd = 5,        &
      luncrt = 6,        &
      mxmenu = 28 + 1,   &      ! + 1 allows for the 99 = done choice
      halfm  = (mxmenu + 4) / 2 ! 4 here allows for -2, -1, 0, and mxmenu + 1

   real, parameter ::    &
      inches_to_meters = 0.0254, &
      meters_to_inches = 1.0 / inches_to_meters

   character, parameter :: &
      format * 11 = 'unformatted'

!  Variables:

   integer :: &
      choice, i, i1, ib, ios, n, nblock, numf

   logical :: &
      cell_centered, cr, eof, formatted

   character :: &
      answer * 1, filename * 128, filename2 * 128, menu (-2 : mxmenu + 1) * 36

   type (grid_type), pointer, dimension (:) :: &
      block                                    ! Surface patch array

!  Storage:

   data menu / &
      '  -2: Start over',                      &
      '  -1: Undo last transformation',        &
      '   0: Review data',                     &
      '   1: Translate X',                     &
      '   2: Scale X',                         &
      '   3: Translate Y',                     &
      '   4: Scale Y',                         &
      '   5: Translate Z',                     &
      '   6: Scale Z',                         &
      '   7: Reflect Z about the XY-plane',    &
      '   8: Reflect X about the YZ-plane',    &
      '   9: Reflect Y about the ZX-plane',    &
      '  10: Reverse the order (1:N)',         &
      '  11: Switch X & Y',                    &
      '  12: Switch Y & Z',                    &
      '  13: Switch Z & X',                    &
      '  14: Scale X & Y & Z the same way',    &
      '  15: Rotate (X,Y) about (Xc,Yc)',      &
      '  16: Rotate (Y,Z) about (Yc,Zc)',      &
      '  17: Rotate (Z,X) about (Zc,Xc)',      &
      '  18: Convert inches to meters',        &
      '  19: Convert meters to inches',        &
      '  20: Rotate (X,Y,Z) about line PQ',    &
      '  21: Transpose a surface patch',       &
      '  22: Reverse i indices of a patch',    &
      '  23: Reverse j indices of a patch',    &
      '  24: Extract a surface patch subset',  &
      '  25: Impose a "fissure" on a patch',   &
      '  26: Impose fissures at all edges',    &
      '  27: Edit the coordinates of a pt.',   &
      '  28: Match patch edges; sym. y -> 0.', &
      '  99: Done',                            &
      '                                   '/   ! Last ' ' eases display
                                               ! of menu as two columns.
!  Execution:

!  Prompt for the input file and check for its existence and apparent form:

   ios = 1
   do while (ios /= 0)

      call reads (luncrt, 'Input grid file name (PLOT3D /mgrid): ', &
                  lunkbd, filename, cr, eof)
      if (eof) go to 999  ! Quit
      if (cr) cycle

      call determine_grid_form (filename, lunin, formatted, ios)

   end do

   i1 = 1;  if (formatted) i1 = 3

   open (lunin, file=filename, form=format(i1:11), status='OLD')

!  Read all blocks.  Trying to process one block at a time isn't practical
!  if more than one transformation per run is allowed for.

   call xyzq_read (lunin, -1, formatted, nblock, numf, cell_centered, &
                   block, ios)

   if (ios /= 0) then
      write (luncrt, '(/, a)') ' Trouble reading the grid.'
      go to 999
   end if

!  Display the block dimensions to reassure the user:

   write (luncrt, '(/, (i6, 2x, 3i5))') &
      (ib, block(ib)%ni, block(ib)%nj, block(ib)%nk, ib = 1, nblock)

!  Loop over possibly several transformations per run:
!  ---------------------------------------------------

200 continue

      write (luncrt, '(/, (a, 10x, a))') &           ! Basically, i=1:mxmenu/2
         (menu(i), menu(i + halfm), i = -2, halfm - 3)

      if (mod (mxmenu, 2) == 0) write (luncrt, '(a)')

210   call readi (luncrt, 'Make a choice. eof (^D) means no more. ', &
                  lunkbd, choice, cr, eof)

      if (choice == 99) go to 800
      if (eof) go to 800
      if (cr)  go to 210
      if (choice < -2 .or. choice >= mxmenu) go to 210

      if (choice > 0 .or. choice == -2) then   ! Save current values:

! ***    xlast = x ! Leave hooks in case the memory usage is worth it

      end if

      if (choice == -2) then      ! "Start over" from scratch:

! ***    x = xorig

         write (luncrt, '(/, a)') ' Starting over is not implemented.  Sorry.'

      else if (choice == -1) then ! "Undo" previous operation:

! ***    x = xlast
        
         write (luncrt, '(/, a)') ' "undo" is not implemented.  Sorry.'

      else if (choice == 0) then  ! "Review": Display the data.

         write (luncrt, '(/, a, /)') ' First few (x,y,z)s of block 1:'
         n = min (10, block(1)%ni)
         write (luncrt, '(1x, 1p, 3e16.7)') &
            (block(1)%x(i,1,1), block(1)%y(i,1,1), block(1)%z(i,1,1), i = 1, n)

      else if (choice == 28) then

!        Matching edges is odd in that there is just one call for all blocks:

         call match_patch_edges (nblock, block, 1.e-4)

      else

         if (choice <= 20) then  ! Original choices from ADJUST_GRID

!           Process all surface patches the same way, one at a time:

            do ib = 1, nblock

               n   = block(ib)%ni * block(ib)%nj

               cr  = .false. ! TRANSFORM may want to return to the menu
               eof = .false.

!              Internal procedure similar to RESHAPE3D, for all patches:

               call transform (ib, n, block(ib)%x, block(ib)%y, block(ib)%z) 

               if (cr .or. eof) exit

            end do ! Next block

         else ! choice >= 21 but not 28: treat patches selectively

            call adjust_patch ()

         end if

      end if

!!!go to 200
   go to 210  ! Repeating the menu is rather redundant

800 continue

!  Save transformed grid?

   filename2 = ' '
   call reads (luncrt, 'Output file name [^D = quit without saving]: ', &
               lunkbd, filename2, cr, eof)

   if (.not. eof) then

      if (cr) filename2 = filename
      if (filename2 == filename) then

         call readc (luncrt, 'Overwrite the input file?  [Y|N; <cr>=Yes]: ', &
                     lunkbd, answer, cr, eof)
         if (eof) go to 800

         cr = cr .or. answer == 'Y'

         if (.not. cr) go to 800

      end if

      formatted = .false.
      call ready (luncrt, 'Formatted (Y) or unformatted (N); <cr>=No: ', &
                  lunkbd, formatted, cr, eof)
      if (eof) go to 800

      i1 = 1;  if (formatted) i1 = 3

      open (lunout, file=filename2, form=format(i1:11), status='unknown')

      call xyz_write (lunout, formatted, nblock, block, ios)

      if (ios /= 0) then
         write (luncrt, '(/, a)') ' Trouble saving the grid.'
      end if

   end if

   do ib = 1, nblock
      deallocate (block(ib)%x, block(ib)%y, block(ib)%z)
   end do
 
999 continue

! *** STOP ! Avoid system dependencies.


!  SURFACE_PATCHES internal procedures:

   contains

!     -------------------------------------
      subroutine transform (ib, n, x, y, z)
!     -------------------------------------

!     The x, y, z arguments simplify reuse of code from RESHAPE3D.
!     Each call in this application operates on one grid block with each
!     coordinate contiguous in memory.
!     Coordinates are changed in place.
!     Suppress local prompts if ib is not 1.
!     Choice, cr, eof, and the luns are inherited from the calling program.

!     Arguments:

      integer, intent (in)                   :: ib, n
      real,    intent (inout), dimension (n) :: x, y, z

!     Local variables:

      integer    :: i, j
      real       :: temp, xt, yt, zt
      real, save :: angle, p, px, py, pz, q, qx, qy, qz, scale, shift
      logical    :: prompt

!     Execution:

      prompt = ib == 1

      select case (choice)

      case (1) ! "Translate X":

         if (prompt) then
            call readr (luncrt, '   Enter X shift (+ or -): ', &
                        lunkbd, shift, cr, eof)
            if (cr .or. eof) go to 999
         end if

         x(:) = x(:) + shift

      case (2) ! "Scale X":

         if (prompt) then
            call readr (luncrt, '   Enter X scale (+ or -): ', &
                        lunkbd, scale, cr, eof)
            if (cr .or. eof) go to 999
         end if

         x(:) = x(:) * scale

      case (3) ! "Translate Y":

         if (prompt) then
            call readr (luncrt, '   Enter Y shift (+ or -): ', &
                        lunkbd, shift, cr, eof)
            if (cr .or. eof) go to 999
         end if

         y(:) = y(:) + shift

      case (4) ! "Scale Y":

         if (prompt) then
            call readr (luncrt, '   Enter Y scale (+ or -): ', &
                        lunkbd, scale, cr, eof)
            if (cr .or. eof) go to 999
         end if

         y(:) = y(:) * scale

      case (5) ! "Translate Z"

         if (prompt) then
            call readr (luncrt, '   Enter Z shift (+ or -): ', &
                        lunkbd, shift, cr, eof)
            if (cr .or. eof) go to 999
         end if

         z(:) = z(:) + shift

      case (6) ! "Scale Z"

         if (prompt) then
            call readr (luncrt, '   Enter Z scale (+ or -): ', &
                        lunkbd, scale, cr, eof)
            if (cr .or. eof) go to 999
         end if

         z(:) = z(:) * scale

      case (7) ! "Reflect Z about the XY-plane":

         z(:) = -z(:)

      case (8) ! "Reflect X about the YZ-plane":

         x(:) = -x(:)

      case (9) ! "Reflect Y about the ZX-plane":

         y(:) = -y(:)

      case (10) ! "Reverse the order":

         do i = 1, (n + 1) / 2
            j = n + 1 - i
            temp = x(i)
            x(i) = x(j)
            x(j) = temp
            temp = y(i)
            y(i) = y(j)
            y(J) = temp
            temp = z(i)
            z(i) = z(j)
            z(j) = temp
         end do

      case (11) ! "Switch X and Y":

         do i = 1, n
            temp = x(i)
            x(i) = y(i)
            y(i) = temp
         end do

      case (12) ! "Switch Y and Z":

         do i = 1, n
            temp = y(i)
            y(i) = z(i)
            z(i) = temp
         end do

      case (13) ! "Switch Z and X":

         do i = 1, n
            temp = x(i)
            x(i) = z(i)
            z(i) = temp
         end do
       
      case (14) ! "Scale X & Y & Z"

         if (prompt) then
            call readr (luncrt, '   Enter scale (+ or -): ', &
                        lunkbd, scale, cr, eof)
            if (cr .or. eof) go to 999
         end if

         x(:) = x(:) * scale
         y(:) = y(:) * scale
         z(:) = z(:) * scale

      case (15) ! "Rotate (X,Y)":

         if (prompt) then
            call readr (luncrt, '   Degrees anticlockwise: ', &
                        lunkbd, angle, cr, eof)
            if (cr .or. eof) go to 999

150         write (luncrt, '(a)', advance='no')'  (Center (Xc, Yc): '
            read  (lunkbd, *, err=150) p, q
            write (luncrt, '(a)')
         end if

         call rotate2d (n, x, y, angle, p, q)

      case (16) ! "Rotate (Y,Z)":

         if (prompt) then
            call readr (luncrt, '   Degrees anticlockwise: ', &
                        lunkbd, angle, cr, eof)
            if (cr .or. eof) go to 999

160         write (luncrt, '(a)', advance='no')'  (Center (Yc, Zc): '
            read  (lunkbd, *, err=160) p, q
            write (luncrt, '(a)')
         end if

         call rotate2d (n, y, z, angle, p, q)

      case (17) ! "Rotate (Z,X)":

         if (prompt) then
            call readr (luncrt, '   Degrees anticlockwise: ', &
                        lunkbd, angle, cr, eof)
            if (cr .or. eof) go to 999

170         write (luncrt, '(a)', advance='no') '  (Center (Zc, Xc): '
            read  (lunkbd, *, err=170) p, q
            write (luncrt, '(a)')
         end if

         call rotate2d (n, z, x, angle, p, q)

      case (18) ! "Convert inches to meters"

         x(:) = x(:) * inches_to_meters
         y(:) = y(:) * inches_to_meters
         z(:) = z(:) * inches_to_meters

      case (19) ! "Convert meters to inches"

         x(:) = x(:) * meters_to_inches
         y(:) = y(:) * meters_to_inches
         z(:) = z(:) * meters_to_inches

      case (20) ! "Rotate (X,Y,Z) about line joining P and Q"

         if (prompt) then
            call readr (luncrt, '   Degrees (RH rule; thumb P -> Q): ', &
                        lunkbd, angle, cr, eof)
            if (cr .or. eof) go to 999

            write (luncrt, '(a)', advance='no') '  (Px, Py, Pz): '
            read  (lunkbd, *) px, py, pz
            write (luncrt, '(a)', advance='no') '  (Qx, Qy, Qz): '
            read  (lunkbd, *) qx, qy, qz
            write (luncrt, '(a)')
         end if

         call rotate3d (n, x, y, z, angle, px, py, pz, qx, qy, qz)

      case default

         ! Shouldn't be possible

      end select

999   return

      end subroutine transform

!     -------------------------------------------------------------
      subroutine adjust_patch ()
!     -------------------------------------------------------------

!     Apply the indicated operation to the selected patch(es) only.

      integer, parameter :: nf = 0             ! No processing of flow data
      real,    parameter :: unchanged   = 99.  ! May not edit all of x, y, & z
      logical, parameter :: release_in  = .true.
      logical, parameter :: release_out = .true.

      integer :: i, i1, i2, ib, ib1, ib2, j, j1, j2, ni, nj
      real    :: d, xt, yt, zt

      type (grid_type) :: tpatch

!     Execution:

      ib1 = 1;  ib2 = 1

      if (nblock > 1) then
         write (luncrt, '(a)', advance='no') '   1st & last patch to process: '
         read  (lunkbd, *) ib1, ib2
         ib1 = max (ib1, 1);    ib2 = min (ib2, nblock)
      end if
                   
      do ib = ib1, ib2

         ni = block(ib)%ni;  nj = block(ib)%nj
         tpatch%ni = ni;     tpatch%nj = nj;     tpatch%nk = 1

         if (choice >= 24) write (luncrt, '(a, i4, a, 2i6)') &
            '   Dimensions for patch', ib, ':', ni, nj

         select case (choice)

         case (21) ! "Transpose a surface patch"

            call transpose_patch (block(ib), nf)

         case (22) ! "Reverse i indices of a patch"

            call reverse_patch_i (block(ib), nf)

         case (23) ! "Reverse j indices of a patch"

            call reverse_patch_j (block(ib), nf)

         case (24) ! "Extract a surface patch subset"

            write (luncrt, '(a)', advance='no') '   i1, i2, j1, j2 to extract: '
            read  (lunkbd, *) i1, i2, j1, j2

            ni = i2 - i1 + 1;  tpatch%ni = ni
            nj = j2 - j1 + 1;  tpatch%nj = nj

            allocate (tpatch%x(ni,nj,1), tpatch%y(ni,nj,1), tpatch%z(ni,nj,1))

            tpatch%x(1:ni,1:nj,1) = block(ib)%x(i1:i2,j1:j2,1)
            tpatch%y(1:ni,1:nj,1) = block(ib)%y(i1:i2,j1:j2,1)
            tpatch%z(1:ni,1:nj,1) = block(ib)%z(i1:i2,j1:j2,1)

            call update_patch (tpatch, nf, release_in, release_out, block(ib))

         case (25) ! "Impose a fissure on a surface patch"

            if (ib == ib1) write (luncrt, '(a)') &
         '   A "fissure" is defined by i1:i2, j1:j2, probably one point wide.',&
         '   A square index range means make the fissure along the diagonal.'
            write (luncrt, '(a)', advance='no') '   Enter two index pairs: '
            read  (lunkbd, *) i1, i2, j1, j2
 25         write (luncrt, '(a)', advance='no') &
         '   Depth? [xyz units; d < 0. => |d| = multiple of local cell diag. '
            read  (lunkbd, *) d
            if (d <= 0.) then
               write (luncrt, '(a)') '   d < 0 option is not implemented yet.'
               go to 25
            end if

            allocate (tpatch%x(ni,nj,1), tpatch%y(ni,nj,1), tpatch%z(ni,nj,1))

            call surface_fissure (ni, nj, i1, i2, j1, j2, d, block(ib), tpatch)

            call update_patch (tpatch, nf, release_in, release_out, block(ib))

         case (26) ! "Impose fissures along all patch edges"

            if (ib == ib1) then
               first_unit_normal = .true.  ! See surface_patch_utilities

               write (luncrt, '(a)', advance='no') &
                  '   Depth? [xyz units] '
               read  (lunkbd, *) d
            else
               first_unit_normal = .false.
            end if

            allocate (tpatch%x(ni,nj,1), tpatch%y(ni,nj,1), tpatch%z(ni,nj,1))

            call edge_fissures (ib, ni, nj, d, block(ib), tpatch)

            call update_patch (tpatch, nf, release_in, release_out, block(ib))

         case (27) ! "Edit the coordinates of a point"

271         write (luncrt, '(a)', advance='no') &
               '   Enter i, j of the point to edit: '
            read  (lunkbd, *, err=271) i, j
            write (luncrt, '(a, 3es16.8)') '   Current coordinates: ', &
               block(ib)%x(i,j,1), block(ib)%y(i,j,1), block(ib)%z(i,j,1)
272         write (luncrt, '(a)', advance='no') &
               '   Enter desired x, y, z.  Use 99 to mean leave alone: '
            read  (lunkbd, *, err=272) xt, yt, zt
            if (xt /= unchanged) block(ib)%x(i,j,1) = xt
            if (yt /= unchanged) block(ib)%y(i,j,1) = yt
            if (zt /= unchanged) block(ib)%z(i,j,1) = zt

         case default

            ! Shouldn't be possible

         end select

      end do  ! Next patch to adjust

      end subroutine adjust_patch

!     --------------------------------------------------------------------------
      subroutine surface_fissure (ni, nj, i1, i2, j1, j2, d, patch_i, patch_o)
!
!     The input patch is deformed in the indicated index range, normal to the
!     undeformed surface by the indicated [small] amount.  An associated volume
!     grid can then be perturbed correspondingly using the RADIAL_INTERP scheme
!     available from this author, assuming it obeys the too-handy-to-ignore
!     assumption that k = 1 is the surface being perturbed.  This "fissure"
!     option was prompted by asteroid studies.
!
!     If i2 - i1 = j2 - j1, assume the fissure is intended to be diagonal.
!     --------------------------------------------------------------------------

!     Arguments:

      integer, intent (in) :: ni, nj             ! Patch dimensions
      integer, intent (in) :: i1, i2, j1, j2     ! Extent of "fissure"
      real,    intent (in) :: d                  ! Depth of fissure;
                                                 ! d < 0. means the depth is
                                                 ! |d| x local cell diagonal
      type (grid_type), intent (in)  :: patch_i  ! Undeformed patch
      type (grid_type), intent (out) :: patch_o  ! Deformed patch

!     Local constants:

      real, parameter :: p = 0., q = 0.  !  Fractional indices <--> pt. (i,j)

!     Local variables:

      integer :: i, j
      real    :: un(3)

!     Execution:

      patch_o%x(:,:,:) = patch_i%x(:,:,:)
      patch_o%y(:,:,:) = patch_i%y(:,:,:)
      patch_o%z(:,:,:) = patch_i%z(:,:,:)

      if (j2 - j1 /= i2 - i1) then
         do j = j1, j2
            do i = i1, i2
               call surface_normal (ni, nj, patch_i%x, patch_i%y, patch_i%z, &
                                    i, j, p, q, un)
               patch_o%x(i,j,1) = patch_o%x(i,j,1) - d*un(1)
               patch_o%y(i,j,1) = patch_o%y(i,j,1) - d*un(2)
               patch_o%z(i,j,1) = patch_o%z(i,j,1) - d*un(3)
            end do
         end do
      else  ! Go along the index diagonal
         j = j1
         do i = i1, i2
            call surface_normal (ni, nj, patch_i%x, patch_i%y, patch_i%z, &
                                 i, j, p, q, un)
            patch_o%x(i,j,1) = patch_o%x(i,j,1) - d*un(1)
            patch_o%y(i,j,1) = patch_o%y(i,j,1) - d*un(2)
            patch_o%z(i,j,1) = patch_o%z(i,j,1) - d*un(3)
            j = j + 1
         end do
      end if

      end subroutine surface_fissure

!     --------------------------------------------------------------------------
      subroutine edge_fissures (ib, ni, nj, d, patch_i, patch_o)
!
!     The input patch is deformed along all edges with the same kind of crude
!     "fissure" originally implemented for one patch at a time along an interior
!     grid line, so that the full asteroid surface grid for which this option is
!     intended ends up with fissures along all edges in both index directions.
!     An associated volume grid can then be perturbed correspondingly using the
!     RADIAL_INTERP scheme available from this author, assuming that it obeys
!     the too-handy-to-ignore assumption that k = 1 is the surface being
!     perturbed.  This "fissure" option was prompted by asteroid studies.
!
!     Note that the deformations along common edges of two patches will not be
!     identical because of the way the SURFACE_NORMAL does not use information
!     from neighboring patches.  The flow solver ancillary utility ZB_FIX can
!     force exact matches using the interface file that accompanies the volume
!     grid following use of RADIAL_INTERP on the fissured surface and the smooth
!     surface volume grid to get the slightly flawed fissured volume grid.
!
!     (Later:)  This version uses the later SURFACE_VECTORS, which DOES attempt
!     to use neighboring patch information to ensure identical results at common
!     edge points.  Initially, this driving routine serves to test the newer
!     surface normal utility (which is restricted to grid points, unlike the
!     SURFACE_NORMAL utility).
!     --------------------------------------------------------------------------

!     Arguments:

      integer, intent (in) :: ib                 ! Current patch number
      integer, intent (in) :: ni, nj             ! Patch dimensions
      real,    intent (in) :: d                  ! Depth of fissure
      type (grid_type), intent (in)  :: patch_i  ! Undeformed patch
      type (grid_type), intent (out) :: patch_o  ! Deformed patch

!     Local constants:

!!!   real, parameter :: p = 0., q = 0.  ! Lower left fractional indices;
!!!                                      ! SURFACE_NORMAL handles i = ni, j = nj

!     Local variables:

      integer :: i, j
      real    :: tani(3), tanj(3), un(3)

!     Execution:

      patch_o%x(:,:,:) = patch_i%x(:,:,:)
      patch_o%y(:,:,:) = patch_i%y(:,:,:)
      patch_o%z(:,:,:) = patch_i%z(:,:,:)

      do i = 1, ni
         do j = 1, nj, nj - 1

!!!         call surface_normal (ni, nj, patch_i%x, patch_i%y, patch_i%z, &
!!!                              i, j, p, q, un)

            call surface_vectors (nblock, block, ib, i, j, tani, tanj, un)

            patch_o%x(i,j,1) = patch_o%x(i,j,1) - d*un(1)
            patch_o%y(i,j,1) = patch_o%y(i,j,1) - d*un(2)
            patch_o%z(i,j,1) = patch_o%z(i,j,1) - d*un(3)
         end do
      end do

      do j = 2, nj - 1
         do i = 1, ni, ni - 1

!!!         call surface_normal (ni, nj, patch_i%x, patch_i%y, patch_i%z, &
!!!                              i, j, p, q, un)

            call surface_vectors (nblock, block, ib, i, j, tani, tanj, un)

            patch_o%x(i,j,1) = patch_o%x(i,j,1) - d*un(1)
            patch_o%y(i,j,1) = patch_o%y(i,j,1) - d*un(2)
            patch_o%z(i,j,1) = patch_o%z(i,j,1) - d*un(3)
         end do
      end do

      end subroutine edge_fissures

   end program surface_patches
