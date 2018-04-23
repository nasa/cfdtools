!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program adjust_flow
!
!  Description:
!
!     This adaptation of ADJUST_GRID offers just the transformations that would
!     affect the velocity components of a flow solution containing the state
!     variables, and applies them to both the flow data and the corresponding
!     grid.  Both input files should be either cell-centered or vertex-centered.
!
!     For convenience, the options to shift x, y, or z are retained.
!
!  Procedures:
!
!     XYZQ_IO    I/O utilities for PLOT3D-type multiblock grid files
!     READER     Prompting utility
!     ROTATE2D   2-space rotation utility
!
!  History:
!
!     07/15/08  D.A.Saunders  ADJUST_GRID version being adapted here.
!     07/19/12    "      "    ADJUST_FLOW adapted from ADJUST_GRID.
!
!  Author:  David Saunders, ERC, Inc. at NASA Ames Research Center, CA
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
      lunfi  = 3,        &
      lunfo  = 4,        &
      lunkbd = 5,        &
      luncrt = 6,        &
      mxmenu = 14,       &
      halfm  = (mxmenu + 4) / 2 ! 4 here allows for -2, -1, 0, and mxmenu + 1

   character, parameter :: &
      format * 11 = 'unformatted'

!  Variables:

   integer :: &
      choice, i, i1, ib, ios, n, nblock, nf, ni, nj, nk, ns

   logical :: &
      cell_centered, cr, eof, formatted

   character :: &
      answer*1, filenamef1*96, filenamef2*96, filenameg1*96, filenameg2*96, &
      menu (-2 : mxmenu + 1) * 35

!  Composite data types:

   type (grid_type), pointer, dimension (:) :: &
      block                                    ! Grid block array

!  Storage:

   data menu / &
      '  -2: Start over',                      &
      '  -1: Undo last transformation',        &
      '   0: Review data',                     &
      '   1: Translate X',                     &
      '   2: Translate Y',                     &
      '   3: Translate Z',                     &
      '   4: Reflect Z about the XY-plane',    &
      '   5: Reflect X about the YZ-plane',    &
      '   6: Reflect Y about the ZX-plane',    &
      '   7: Switch X & Y',                    &
      '   8: Switch Y & Z',                    &
      '   9: Switch Z & X',                    &
      '  10: Rotate (X,Y) about (Xc,Yc)',      &
      '  11: Rotate (Y,Z) about (Yc,Zc)',      &
      '  12: Rotate (Z,X) about (Zc,Xc)',      &
      '  13: Rotate (X,Y,Z) about line PQ',    &
      '  99: Done',                            &
      '                                   '/   ! Last ' ' eases display
                                               ! of menu as two columns.
!  Execution:

!  Prompt for the input grid and check for its existence and apparent form:

   ios = 1
   do while (ios /= 0)
      call reads (luncrt, 'Input grid file name (PLOT3D /mgrid): ', &
                  lunkbd, filenameg1, cr, eof)
      if (eof) go to 999  ! Quit
      if (cr) cycle

      call determine_grid_form (filenameg1, lunin, formatted, ios)
   end do

   i1 = 1;  if (formatted) i1 = 3

   open (lunin, file=filenameg1, form=format(i1:11), status='OLD')

   ios = 1
   do while (ios /= 0)
      call reads (luncrt, 'Input flow file name (PLOT3D /mgrid): ', &
                  lunkbd, filenamef1, cr, eof)
      if (eof) go to 999  ! Quit
      if (.not. cr) ios = 0
   end do  ! The same formatting is assumed for both files

   open (lunfi, file=filenamef1, form=format(i1:11), status='OLD')

!  Read all blocks.  Trying to process one block at a time isn't practical
!  if more than one transformation per run is allowed for.

   call xyzq_read (lunin, lunfi, formatted, nblock, nf, cell_centered, &
                   block, ios)

   if (ios /= 0) then
      write (luncrt, '(/, a)') ' Trouble reading the grid or flow solution.'
      go to 999
   end if

!  Display the block dimensions to reassure the user:

   write (luncrt, '(/, (i6, 2x, 3i5))') &
      (ib, block(ib)%ni, block(ib)%nj, block(ib)%nk, ib = 1, nblock)

   write (luncrt, '(/, a, i3)') ' Number of flow variables found:', nf
   ios = 1
   do while (ios /= 0)
      call readi (luncrt, 'To locate u/v/w, enter the number of species: ', &
                  lunkbd, ns, cr, eof)
      if (eof) go to 999
      if (.not. cr) ios = 0
      if (ns < 1 .or. ns > 30) ios = 1
   end do

!  Loop over possibly several transformations per run:
!  ---------------------------------------------------

200 continue

      write (luncrt, '(/, (a, 10x, a))') &           ! Basically, i=1:mxmenu/2
         (menu(i), menu(i + halfm), i = -2, halfm - 3)

      if (mod (mxmenu, 2) == 0) write (luncrt, '(a)')

210   call readi (luncrt, 'Make a choice. 99 or eof (^D) means no more. ', &
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
         write (luncrt, '(1x, 3es16.7)') &
            (block(1)%x(i,1,1), block(1)%y(i,1,1), block(1)%z(i,1,1), i = 1, n)
         write (luncrt, '(/, a, /)') ' First few flow variables of block 1:'
         do i = 1, n
            write (luncrt, '(1x, 20es16.7)') block(1)%q(:,i,1,1)
         end do  ! 20 here is meant to be more than enough

      else  ! Process blocks one at a time:

         do ib = 1, nblock

            ni  = block(ib)%ni
            nj  = block(ib)%nj
            nk  = block(ib)%nk
            n   = ni*nj*nk

            cr  = .false. ! TRANSFORM may want to return to the menu
            eof = .false.

!           Internal procedure efficient for x/y/z, less so for f:

            call transform (ib, ni, nj, nk, nf, n, &
                            block(ib)%x, block(ib)%y, block(ib)%z, block(ib)%q)

            if (cr .or. eof) exit

         end do ! Next block

      end if

!!!go to 200
   go to 210  ! Repeating the menu is rather redundant

800 continue

!  Save transformed files?

   filenameg2 = ' '
   call reads (luncrt, 'Output grid name [^D = quit without saving]: ', &
               lunkbd, filenameg2, cr, eof)

   if (.not. eof) then
      if (cr) filenameg2 = filenameg1
      if (filenameg2 == filenameg1) then
         call readc (luncrt, 'Overwrite the input grid?  [Y|N; <cr>=Yes]: ', &
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

      open (lunout, file=filenameg2, form=format(i1:11), status='unknown')

      call xyz_write (lunout, formatted, nblock, block, ios)

      if (ios /= 0) then
         write (luncrt, '(/, a)') ' Trouble saving the grid.'
         go to 999
      end if

810   continue
      filenamef2 = ' '
      call reads (luncrt, 'Output flow name [^D = quit without saving]: ', &
                  lunkbd, filenamef2, cr, eof)

     if (.not. eof) then
        if (cr) filenamef2 = filenamef1
        if (filenamef2 == filenamef1) then
           call readc (luncrt, 'Overwrite the input flow?  [Y|N; <cr>=Yes]: ', &
                       lunkbd, answer, cr, eof)
           if (eof) go to 810
           cr = cr .or. answer == 'Y'
           if (.not. cr) go to 810
        end if

!!!     formatted = .false.
!!!     call ready (luncrt, 'Formatted (Y) or unformatted (N); <cr>=No: ', &
!!!                 lunkbd, formatted, cr, eof)
!!!     if (eof) go to 800

!!!     i1 = 1;  if (formatted) i1 = 3

        open (lunfo, file=filenamef2, form=format(i1:11), status='unknown')

        call q_write (lunfo, formatted, nblock, nf, block, ios)

        if (ios /= 0) then
           write (luncrt, '(/, a)') ' Trouble saving the flow.'
           go to 999
        end if

     end if

   end if

   do ib = 1, nblock
      deallocate (block(ib)%x, block(ib)%y, block(ib)%z, block(ib)%q)
   end do
 
999 continue

! *** STOP ! Avoid system dependencies.


!  ADJUST_FLOW internal procedure:

   contains

!     --------------------------------------------------------
      subroutine transform (ib, ni, nj, nk, nf, n, x, y, z, f)
!     --------------------------------------------------------

!     Each call in this application operates on one grid block with each
!     coordinate contiguous in memory (but all the flow variables at a grid
!     point contiguous, in contrast to the PLOT3D function file format).
!     The flow data are ignored for the translate x/y/z options.
!
!     Suppress local prompts if ib is not 1.
!     Choice, cr, eof, ns, and the luns are inherited from the calling program.

!     Arguments:

      integer, intent (in)                   :: ib, ni, nj, nk, nf, n
      real,    intent (inout), dimension (n) :: x, y, z
      real,    intent (inout)                :: f(nf,ni,nj,nk)

!     Local variables:

      integer    :: i, j, k
      real       :: temp
      real, save :: angle, p, px, py, pz, q, qx, qy, qz, shift
      logical    :: prompt
      real, allocatable, dimension (:,:,:) :: vx, vy, vz

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

      case (2) ! "Translate Y":

         if (prompt) then
            call readr (luncrt, '   Enter Y shift (+ or -): ', &
                        lunkbd, shift, cr, eof)
            if (cr .or. eof) go to 999
         end if

         y(:) = y(:) + shift

      case (3) ! "Translate Z"

         if (prompt) then
            call readr (luncrt, '   Enter Z shift (+ or -): ', &
                        lunkbd, shift, cr, eof)
            if (cr .or. eof) go to 999
         end if

         z(:) = z(:) + shift

      case (4) ! "Reflect Z about the XY-plane":

         z(:) = -z(:)
         f(ns+3,:,:,:) = -f(ns+3,:,:,:)

      case (5) ! "Reflect X about the YZ-plane":

         x(:) = -x(:)
         f(ns+1,:,:,:) = -f(ns+1,:,:,:)

      case (6) ! "Reflect Y about the ZX-plane":

         y(:) = -y(:)
         f(ns+2,:,:,:) = -f(ns+2,:,:,:)

      case (7) ! "Switch X and Y":

         do i = 1, n
            temp = x(i)
            x(i) = y(i)
            y(i) = temp
         end do

         do k = 1, nk
            do j = 1, nj
               do i = 1, ni
                  temp = f(ns+1,i,j,k)
                  f(ns+1,i,j,k) = f(ns+2,i,j,k)
                  f(ns+2,i,j,k) = temp
               end do
            end do
         end do

      case (8) ! "Switch Y and Z":

         do i = 1, n
            temp = y(i)
            y(i) = z(i)
            z(i) = temp
         end do

         do k = 1, nk
            do j = 1, nj
               do i = 1, ni
                  temp = f(ns+2,i,j,k)
                  f(ns+2,i,j,k) = f(ns+3,i,j,k)
                  f(ns+3,i,j,k) = temp
               end do
            end do
         end do

      case (9) ! "Switch Z and X":

         do i = 1, n
            temp = x(i)
            x(i) = z(i)
            z(i) = temp
         end do

         do k = 1, nk
            do j = 1, nj 
               do i = 1, ni
                  temp = f(ns+1,i,j,k)
                  f(ns+1,i,j,k) = f(ns+3,i,j,k)
                  f(ns+3,i,j,k) = temp
               end do
            end do
         end do

      case (10) ! "Rotate (X,Y)":

         if (prompt) then
            call readr (luncrt, '   Degrees anticlockwise: ', &
                        lunkbd, angle, cr, eof)
            if (cr .or. eof) go to 999

150         write (luncrt, '(a)', advance='no')'  (Center (Xc, Yc): '
            read  (lunkbd, *, err=150) p, q
            write (luncrt, '(a)')
         end if

         call rotate2d (n, x, y, angle, p, q)

         allocate (vx(ni,nj,nk), vy(ni,nj,nk))
         vx(:,:,:) = f(ns+1,:,:,:)
         vy(:,:,:) = f(ns+2,:,:,:)
         call rotate2d (n, vx, vy, angle, p, q)
         f(ns+1,:,:,:) = vx(:,:,:)
         f(ns+2,:,:,:) = vy(:,:,:)
         deallocate (vx, vy)

      case (11) ! "Rotate (Y,Z)":

         if (prompt) then
            call readr (luncrt, '   Degrees anticlockwise: ', &
                        lunkbd, angle, cr, eof)
            if (cr .or. eof) go to 999

160         write (luncrt, '(a)', advance='no')'  (Center (Yc, Zc): '
            read  (lunkbd, *, err=160) p, q
            write (luncrt, '(a)')
         end if

         call rotate2d (n, y, z, angle, p, q)

         allocate (vy(ni,nj,nk), vz(ni,nj,nk))
         vy(:,:,:) = f(ns+2,:,:,:)
         vz(:,:,:) = f(ns+3,:,:,:)
         call rotate2d (n, vy, vz, angle, p, q)
         f(ns+2,:,:,:) = vy(:,:,:)
         f(ns+3,:,:,:) = vz(:,:,:)
         deallocate (vy, vz)

      case (12) ! "Rotate (Z,X)":

         if (prompt) then
            call readr (luncrt, '   Degrees anticlockwise: ', &
                        lunkbd, angle, cr, eof)
            if (cr .or. eof) go to 999

170         write (luncrt, '(a)', advance='no') '  (Center (Zc, Xc): '
            read  (lunkbd, *, err=170) p, q
            write (luncrt, '(a)')
         end if

         call rotate2d (n, z, x, angle, p, q)

         allocate (vx(ni,nj,nk), vz(ni,nj,nk))
         vx(:,:,:) = f(ns+1,:,:,:)
         vz(:,:,:) = f(ns+3,:,:,:)
         call rotate2d (n, vz, vx, angle, p, q)
         f(ns+1,:,:,:) = vx(:,:,:)
         f(ns+3,:,:,:) = vz(:,:,:)
         deallocate (vx, vz)

      case (13) ! "Rotate (X,Y,Z) about line joining P and Q"

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

         allocate (vx(ni,nj,nk), vy(ni,nj,nk), vz(ni,nj,nk))
         vx(:,:,:) = f(ns+1,:,:,:)
         vy(:,:,:) = f(ns+2,:,:,:)
         vz(:,:,:) = f(ns+3,:,:,:)
         call rotate3d (n, vx, vy, vz, angle, px, py, pz, qx, qy, qz)
         f(ns+1,:,:,:) = vx(:,:,:)
         f(ns+2,:,:,:) = vy(:,:,:)
         f(ns+3,:,:,:) = vz(:,:,:)
         deallocate (vx, vy, vz)

      end select

999   return

      end subroutine transform

   end program adjust_flow
