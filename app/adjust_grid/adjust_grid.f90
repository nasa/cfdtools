!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program adjust_grid
!
!  Description:
!
!     This is a variant of the earlier SCALE_GRID that uses the XYZQ_IO package
!     in place of the CFD_IO_PACKAGE.  It is restricted to 3-space multiblock
!     grids, and automatically determines whether the input grid is formatted
!     or not.  The output format must still be prompted for.
!
!     ADJUST_GRID applies one or more transformations (shift, scale, rotate,
!     etc.) to all blocks of the indicated grid.  This version also includes
!     some menu choices that affect only specified blocks.
!
!  Procedures:
!
!     XYZQ_IO    I/O utilities for PLOT3D-type multiblock grid files
!     READER     Prompting utility
!     ROTATE2D   2-space rotation utility
!
!  History:
!
!     02/18/00  D.A.Saunders  Initial SCALE_GRID adaptation of RESHAPE3D, which
!                             operates on 3-column formatted datasets.
!     09/09/06    "      "    Added explicit inches <-> meters options.
!     01/20/08    "      "    ADJUST_GRID adapted from SCALE_GRID to reduce the
!                             number of prompts.
!     07/15/08    "      "    Added general rotation about an axis defined by
!                             two points.
!     02/17/15    "      "    Provided for doing different things to different
!                             blocks as in SURFACE_PATCHES.  First two options:
!                             clean up a symmetry plane block face, and edit
!                             the coordinates of a point.
!
!  Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA
!           Now with ERC, Inc. at NASA ARC.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_block_structure   ! Derived data type for one grid block
   use xyzq_io_module         ! PLOT3D multiblock file I/O package (3-space)

!  Constants:

   implicit none

   integer, parameter :: &
      lunin  = 1,        &
      lunout = 2,        &
      lunkbd = 5,        &
      luncrt = 6,        &
      mxmenu = 22 + 1,   &      ! + 1 allows for the 99 = done choice
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

!  Composite data types:

   type (grid_type), pointer, dimension (:) :: &
      block                                    ! Grid block array

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
      '  21: Clean up a symmetry plane face',  &
      '  22: Edit the coordinates of a pt.',   &
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

      else

         if (choice <= 20) then  ! Same operation on all blocks, in place

!           Process grid blocks one at a time:

            do ib = 1, nblock

               n   = block(ib)%ni*block(ib)%nj*block(ib)%nk    ! # pts. in block
               cr  = .false. ! TRANSFORM may want to return to the menu
               eof = .false.

!              Internal procedure similar to RESHAPE3D:

               call transform (ib, n, block(ib)%x, block(ib)%y, block(ib)%z) 
               if (cr .or. eof) exit

            end do ! Next block

         else  ! Operate on just specified block(s)

            call adjust_block ()

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


!  ADJUST_GRID internal procedure:

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
      real       :: temp
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

      end select

999   return

      end subroutine transform

!     -------------------------------------------------------------
      subroutine adjust_block ()
!     -------------------------------------------------------------

!     Apply the indicated operation to the selected block(s) only.

      real, parameter :: unchanged = 99.  ! May not edit all of x, y, & z
      real, parameter :: zero = 0.

      integer :: i, ib, ib1, ib2, icoord, iface, j, k, ni, nj, nk
      integer :: i1, i2, j1, j2, k1, k2
      real    :: d, xt, yt, zt

!     Execution:

      ib1 = 1;  ib2 = 1

      if (nblock > 1) then
         write (luncrt, '(a)', advance='no') '   1st & last patch to process: '
         read  (lunkbd, *) ib1, ib2
         ib1 = max (ib1, 1);    ib2 = min (ib2, nblock)
      end if

      do ib = ib1, ib2

         ni = block(ib)%ni;  nj = block(ib)%nj;  nk = block(ib)%nk

         write (luncrt, '(a, i4, a, 3i6)') &
            '   Dimensions for block', ib, ':', ni, nj, nk

         select case (choice)

         case (21) ! "Clean up a symmetry plane"

 21         write (luncrt, '(a)', advance='no') &
               '   Enter block face (1-6) and coordinate (1-3) to zero out: '
            read  (lunkbd, *, err=21) iface, icoord

            select case (iface)
               case (1)
                  i1 = 1;  i2 = 1;  j1 = 1;  j2 = nj;  k1 = 1;  k2 = nk
               case (2)
                  i1 = ni; i2 = ni; j1 = 1;  j2 = nj;  k1 = 1;  k2 = nk
               case (3)
                  i1 = 1;  i2 = ni; j1 = 1;  j2 = 1;   k1 = 1;  k2 = nk
               case (4)
                  i1 = 1;  i2 = ni; j1 = nj; j2 = nj;  k1 = 1;  k2 = nk
               case (5)
                  i1 = 1;  i2 = ni; j1 = 1;  j2 = nj;  k1 = 1;  k2 = 1
               case (6)
                  i1 = 1;  i2 = ni; j1 = 1;  j2 = nj;  k1 = nk; k2 = nk
            end select

            select case (icoord)
               case (1)
                  block(ib)%x(i1:i2,j1:j2,k1:k2) = zero
               case (2)
                  block(ib)%y(i1:i2,j1:j2,k1:k2) = zero
               case (3)
                  block(ib)%z(i1:i2,j1:j2,k1:k2) = zero
            end select

         case (22) ! "Edit the coordinates of a point"

221         write (luncrt, '(a)', advance='no') &
               '   Enter i, j, k of the point to edit: '
            read  (lunkbd, *, err=221) i, j, k
            write (luncrt, '(a, 3es16.8)') '   Current coordinates: ', &
               block(ib)%x(i,j,k), block(ib)%y(i,j,k), block(ib)%z(i,j,k)
222         write (luncrt, '(a)', advance='no') &
               '   Enter desired x, y, z.  Use 99 to mean leave alone: '
            read  (lunkbd, *, err=222) xt, yt, zt
            if (xt /= unchanged) block(ib)%x(i,j,k) = xt
            if (yt /= unchanged) block(ib)%y(i,j,k) = yt
            if (zt /= unchanged) block(ib)%z(i,j,k) = zt

         case default

            ! Shouldn't be possible

         end select

      end do

      end subroutine adjust_block

   end program adjust_grid
