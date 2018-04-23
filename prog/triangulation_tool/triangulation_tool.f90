!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   program triangulation_tool

!  Description:
!
!     Read one unstructured dataset (surface or volume), perform the indicated
!  analysis or analyses, and (if indicated) save the result as another such
!  file.  Initially, this serves as a driver for testing utilities added to the
!  triangulation_io module for calculating geometric properties of surface
!  datasets.  Now, it also drives the analogous utilities for volume datasets.
!
!     The dataset may contain more than one zone (but probably doesn't).
!
!     No control file is needed - just some prompts.
!
!  Input Tecplot unstructured surface dataset format (vertex-centered):
!
!     VARIABLES = "X", "Y", "Z", "TEMP", "PRESS", "CP", "MACHE", "ASOUNDE", ...
!     ZONE N=96000, E=32000, F=FEPOINT, ET=TRIANGLE
!     0.000000  0.000000 0.000000 2838.386719 51330.925781 1.552663 0.609412 ...
!     0.000883 -0.007150 0.113643 2838.386719 51330.925781 1.552663 0.609412 ...
!     0.000883  0.000000 0.113868 2838.386719 51330.925781 1.552663 0.609412 ...
!     ::::::::::::::::
!     ::::::::::::::::
!     4.882953  0.000000 0.011285 950.867676 16.506409 -0.001166 5.062649 ...
!     1 2 3
!     4 5 6
!     7 8 9
!     10 11 12
!     ::::::::
!     ::::::::
!     95992 95993 95994
!     95995 95996 95997
!     95998 95999 96000
!
!  Alternative cell-centered Tecplot surface format (DPLR overset grid tools):
!
!     TITLE     = ""
!     VARIABLES = "x"
!     "y"
!     "z"
!     "p"
!     "Chm"
!     "h"
!     "qw"
!     "Re_c"
!     ZONE T="ZONE 001"
!      STRANDID=0, SOLUTIONTIME=0
!      Nodes=9367, Elements=18438, ZONETYPE=FETriangle
!      DATAPACKING=BLOCK
!      VARLOCATION=([4-8]=CELLCENTERED)
!      FACENEIGHBORCONNECTIONS=55014
!      FACENEIGHBORMODE=LOCALONETOONE
!      FEFACENEIGHBORSCOMPLETE=YES
!      DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE )
!      3.756540120E-01 3.601163924E-01 3.451932967E-01 3.309260905E-01  ...
!     ::::::::::::::::
!
!  Format for a Tecplot unstructured volume file:
!
!      variables = x y z
!      zone T="grid", n=      505286 , e=     2759284 ,et=brick,     f=fepoint
!       2.273000031709671E-002  1.339999958872795E-003  0.119850002229214
!       2.288999967277050E-002 -1.100000008591451E-004  0.119819998741150
!      ::::::::::::::::
!       7.481993472351248E-002 -1.206158433734931E-002 -1.534229517647054E-002
!       98307    98308    98309    98309    98310    98310    98310    98310
!       98308    98307    98309    98309    98311    98311    98311    98311
!      ::::::::::::::::
!      263981   238511   270276   270276   270277   270277   270277   270277
!      318887   318885   378503   378503   378505   378505   378505   378505
!
!
!      If, as in this example, each hex cell is actually formed by duplicating
!      vertex 3 of a tetrahedron as vertex 4, and vertex 4 as vertices 6,7,8,
!      then volumes are treated as single tetrahedra -- not as 5 or 6 of them.
!
!  History:
!
!     10/17/14  DAS  Initial adaptation of TRI_TO_TRI as a driver for applying
!                    geometric utilities that have been added to the surface
!                    triangulation module triangulation_io.f90:
!                    x/y/z/data ranges, surface area, and enclose volume.
!                    Further extensions are expected.
!     10/20/14   "   More geometric options along the lines of ADJUST_GRID.
!     10/21/14   "   Surface enter of mass and moments of inertia options.
!     11/17/14   "   Volume grid analogues of the triangulation utilities have
!                    been added to triangulation_io.f90 and are now driven here.
!     12/30/14   "   Added a smoothing option for surface triangulations (as
!                    needed for turning the Itokawa asteroid surface into a
!                    more reasonable case for volume gridding).
!     03/08/15   "   Added an option to write a NASTRAN file.
!
!  Author:  David Saunders, ERC, Inc./NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use tri_header_structure  ! Part of triangulation_io.f90
   use tri_zone_structure    ! Likewise
   use triangulation_io      ! Unstructured data file I/O

   implicit none

!  Constants:

   integer, parameter :: &
      lunin  = 1,        &   ! Input Tecplot triangulated surface dataset
      lunout = 2,        &   ! Output Tecplot multizone triangulation, if any
      lunkbd = 5,        &   ! For keyboard entries
      luncrt = 6,        &   ! For prompts and diagnostics
      mxmenu = 24 + 1,   &   ! The + 1 counts the 99/Done choice
      halfm  = (mxmenu + 4) / 2 ! 4 here allows for -2, -1, 0, and mxmenu + 1

   real, parameter ::    &
      inches_to_meters = 0.0254, &
      meters_to_inches = 1.0 / inches_to_meters

!  Variables:

   integer :: &
      choice, i, ios, iz, l, n

   logical :: &
      cr, eof, nastran, tri, yes

   character (1) :: &
      answer
   character (38) :: &
      menu(-2:mxmenu+1)
   character (80) :: &
      filename_out

!  Derived data types:

   type (tri_header_type) :: &
      header
   type (tri_type), pointer, dimension (:) :: &
      xyzf

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
      '  10: Switch X & Y',                    &
      '  11: Switch Y & Z',                    &
      '  12: Switch Z & X',                    &
      '  13: Scale X & Y & Z the same way',    &
      '  14: Rotate (X,Y) about (Xc,Yc)',      &
      '  15: Rotate (Y,Z) about (Yc,Zc)',      &
      '  16: Rotate (Z,X) about (Zc,Xc)',      &
      '  17: Convert inches to meters',        &
      '  18: Convert meters to inches',        &
      '  19: Rotate (X,Y,Z) about line PQ',    &
      '  20: Calculate surface area',          &
      '  21: [Surface area & enclosed] volume',&
      '  22: Calculate CM & inertia moments',  &
      '  23: Apply rotation matrix from (22)', &
      '  24: Smooth a 1-zone bumpy surface',   &
      '  99: Done',                            &
      '                                   '/   ! Last ' ' eases display
                                               ! of menu as two columns.
!  Execution:
!  !!!!!!!!!!

   call reads (luncrt, 'Input file name (unstructured surface or volume): ', &
               lunkbd, header%filename, cr, eof)
   if (cr .or. eof) go to 99

   header%formatted = .true.  ! Unformatted has not been handled

!  Determine the element type without a prompt:

   call vol_get_element_type (lunin, header, ios)
   if (ios /= 0) go to 99

   write (luncrt, '(a, i3)') &
      ' # vertices per element found for zone 1:', header%nvertices, &
      ' Any further zones are assumed to contain the same element type.'

   tri = header%nvertices == 3

   header%fileform = 1        ! Vertex-centered
   call readi (luncrt, &
               'Vertex-centered or cell-centered fn. data? [1|2; <cr> = 1] ', &
               lunkbd, header%fileform, cr, eof)
   if (eof) go to 99

!  (Re)open and read the input dataset:

   ios = 1  ! Verbose mode

   if (tri) then
      call tri_read (lunin, header, xyzf, ios)  ! Read entire dataset
   else
      call vol_read (lunin, header, xyzf, ios)
   end if

   if (ios /= 0) then
      write (luncrt, '(/, a)') &
         'Trouble reading input dataset: aborting.'
      go to 99
   end if

!  Determine the x/y/z data range.  This also sets a default for the enclosed
!  volume calculations on a triangulated surface, so do it regardless.

   if (tri) then
      call tri_data_range (header, xyzf)
   else
      call vol_data_range (header, xyzf)
   end if

   write (luncrt, '(/, 2a)') &
      '            xmin            xmax            ymin            ymax', &
      '            zmin            zmax'
   write (luncrt, '(6es16.8)') &
      header%xmin, header%xmax, header%ymin, header%ymax, &
      header%zmin, header%zmax
   write (luncrt, '(a, 3es16.8, /)') &
      ' Data range mid-point:          ', header%interior_point(:)

!  Loop over possibly several transformations per run:
!  ---------------------------------------------------

200 continue

      write (luncrt, '(/, (a, 10x, a))') &           ! Basically, i=1:mxmenu/2
         (menu(i), menu(i + halfm), i = -2, halfm - 3)

      if (mod (mxmenu, 2) == 0) write (luncrt, '(a)')

210   call readi (luncrt, 'Make a choice. 99 means no more. ', &
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

         write (luncrt, '(/, a, /)') ' First few (x,y,z)s of zone 1:'
         n = min (10, xyzf(1)%nnodes)
         write (luncrt, '(1x, 3es16.8)') xyzf(1)%xyz(:,1:n)

      else

!        Process grid zones one at a time:

         do iz = 1, header%nzones

            n   = xyzf(iz)%nnodes
            cr  = .false. ! TRANSFORM may want to return to the menu
            eof = .false.

!           Internal procedure similar to RESHAPE3D/ADJUST_GRID:

            call transform (iz, n, xyzf(iz)%xyz)

            if (cr .or. eof) exit

         end do ! Next zone

      end if

!!!go to 200
   go to 210  ! Repeating the menu is rather redundant

800 continue

!  Save transformed triangulation?

   filename_out = ' '
   write (luncrt, '(/, a)') &
      ' Output is Tecplotable unless its name ends in ''.nas'''
   call reads (luncrt, 'Output file name [^D = quit without saving]: ', &
               lunkbd, filename_out, cr, eof)

   if (.not. eof) then

      if (cr) filename_out = header%filename
      if (filename_out == header%filename) then

         call readc (luncrt, 'Overwrite the input file?  [Y|N; <cr>=Yes]: ', &
                     lunkbd, answer, cr, eof)
         if (eof) go to 800

         cr = cr .or. answer == 'Y'

         if (.not. cr) go to 800

      end if

      header%filename = filename_out
      l = len_trim (filename_out)
      nastran = filename_out(l-3:l) == '.nas'

      if (tri) then
         if (nastran) then
            call nas_sf_tri_write (lunout, header, xyzf, ios)
         else
            call tri_write (lunout, header, xyzf, ios)
         end if
      else
         call vol_write (lunout, header, xyzf, ios)
      end if

      if (ios /= 0) then
         write (luncrt, '(/, a)') ' Trouble saving results.'
      end if

   end if

   if (tri) then
      call deallocate_tri_zones (1, header%nzones, header%numf, xyzf, ios)
   else
      call deallocate_vol_zones (1, header%nzones, header%numf, xyzf, ios)
   end if

99 continue

   contains

!     Internal procedure for program TRIANGULATION_TOOL:

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine transform (iz, n, xyz)

!     Operate on the x/y/z coordinates of zone iz of a triangulation.
!     Passing xyzf(iz)%xyz as an argument simplifies nomenclature.
!     Any local prompts are suppressed if iz > 1.
!     Some options have already been implemented over all zones (e.g.,
!     surface area) in triangulation_io.f90.  In these cases, if iz > 1,
!     we just return early.
!     Variables choice, cr, eof, and the luns are inherited from the main
!     program.  So are header and xyzf for calls to utilities that are in
!     triangulation_io.f90.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in)    :: iz        ! Zone number (probably never > 1)
      integer, intent (in)    :: n         ! # nodes in this zone
      real,    intent (inout) :: xyz(3,n)  ! Packed node coordinates, this zone

!     Local variables:

      integer    :: i, niter
      real       :: data_range, fwhm, percent_data_range, temp
      real, save :: angle, p, px, py, pz, q, qx, qy, qz, scale, shift
      logical    :: prompt
      real, allocatable, dimension (:) :: x, y, z

!     Execution:

      prompt = iz == 1

      select case (choice)

      case (1) ! "Translate X":

         if (prompt) then
            call readr (luncrt, '   Enter X shift (+ or -): ', &
                        lunkbd, shift, cr, eof)
            if (cr .or. eof) go to 999
         end if

         xyz(1,:) = xyz(1,:) + shift

      case (2) ! "Scale X":

         if (prompt) then
            call readr (luncrt, '   Enter X scale (+ or -): ', &
                        lunkbd, scale, cr, eof)
            if (cr .or. eof) go to 999
         end if

         xyz(1,:) = xyz(1,:) * scale

      case (3) ! "Translate Y":

         if (prompt) then
            call readr (luncrt, '   Enter Y shift (+ or -): ', &
                        lunkbd, shift, cr, eof)
            if (cr .or. eof) go to 999
         end if

         xyz(2,:) = xyz(2,:) + shift

      case (4) ! "Scale Y":

         if (prompt) then
            call readr (luncrt, '   Enter Y scale (+ or -): ', &
                        lunkbd, scale, cr, eof)
            if (cr .or. eof) go to 999
         end if

         xyz(2,:) = xyz(2,:) * scale

      case (5) ! "Translate Z"

         if (prompt) then
            call readr (luncrt, '   Enter Z shift (+ or -): ', &
                        lunkbd, shift, cr, eof)
            if (cr .or. eof) go to 999
         end if

         xyz(3,:) = xyz(3,:) + shift

      case (6) ! "Scale Z"

         if (prompt) then
            call readr (luncrt, '   Enter Z scale (+ or -): ', &
                        lunkbd, scale, cr, eof)
            if (cr .or. eof) go to 999
         end if

         xyz(3,:) = xyz(3,:) * scale

      case (7) ! "Reflect Z about the XY-plane":

         xyz(3,:) = -xyz(3,:)

      case (8) ! "Reflect X about the YZ-plane":

         xyz(1,:) = -xyz(1,:)

      case (9) ! "Reflect Y about the ZX-plane":

         xyz(2,:) = -xyz(2,:)

      case (10) ! "Switch X and Y":

         do i = 1, n
            temp     = xyz(1,i)
            xyz(1,i) = xyz(2,i)
            xyz(2,i) = temp
         end do

      case (11) ! "Switch Y and Z":

         do i = 1, n
            temp     = xyz(2,i)
            xyz(2,i) = xyz(3,i)
            xyz(3,i) = temp
         end do

      case (12) ! "Switch Z and X":

         do i = 1, n
            temp     = xyz(1,i)
            xyz(1,i) = xyz(3,i)
            xyz(3,i) = temp
         end do

      case (13) ! "Scale X & Y & Z"

         if (prompt) then
            call readr (luncrt, '   Enter scale (+ or -): ', &
                        lunkbd, scale, cr, eof)
            if (cr .or. eof) go to 999
         end if

         xyz(:,:) = xyz(:,:) * scale

      case (14) ! "Rotate (X,Y)":

         if (prompt) then
            call readr (luncrt, '   Degrees anticlockwise: ', &
                        lunkbd, angle, cr, eof)
            if (cr .or. eof) go to 999

140         write (luncrt, '(a)', advance='no')'  (Center (Xc, Yc): '
            read  (lunkbd, *, err=140) p, q
            write (luncrt, '(a)')
         end if

         allocate (x(n), y(n))
         x(:) = xyz(1,:)
         y(:) = xyz(2,:)
         call rotate2d (n, x, y, angle, p, q)
         xyz(1,:) = x(:)
         xyz(2,:) = y(:)
         deallocate (x, y)

      case (15) ! "Rotate (Y,Z)":

         if (prompt) then
            call readr (luncrt, '   Degrees anticlockwise: ', &
                        lunkbd, angle, cr, eof)
            if (cr .or. eof) go to 999

150         write (luncrt, '(a)', advance='no')'  (Center (Yc, Zc): '
            read  (lunkbd, *, err=150) p, q
            write (luncrt, '(a)')
         end if

         allocate (y(n), z(n))
         y(:) = xyz(2,:)
         z(:) = xyz(3,:)
         call rotate2d (n, y, z, angle, p, q)
         xyz(2,:) = y(:)
         xyz(3,:) = z(:)
         deallocate (y, z)

      case (16) ! "Rotate (Z,X)":

         if (prompt) then
            call readr (luncrt, '   Degrees anticlockwise: ', &
                        lunkbd, angle, cr, eof)
            if (cr .or. eof) go to 999

160         write (luncrt, '(a)', advance='no') '  (Center (Zc, Xc): '
            read  (lunkbd, *, err=160) p, q
            write (luncrt, '(a)')
         end if

         allocate (z(n), x(n))
         z(:) = xyz(3,:)
         x(:) = xyz(1,:)
         call rotate2d (n, z, x, angle, p, q)
         xyz(3,:) = z(:)
         xyz(1,:) = z(:)
         deallocate (z, x)

      case (17) ! "Convert inches to meters"

         xyz(:,:) = xyz(:,:) * inches_to_meters

      case (18) ! "Convert meters to inches"

         xyz(:,:) = xyz(:,:) * meters_to_inches

      case (19) ! "Rotate (X,Y,Z) about line joining P and Q"

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

         allocate (x(n), y(n), z(n))
         x(:) = xyz(1,:)
         y(:) = xyz(2,:)
         z(:) = xyz(3,:)
         call rotate3d (n, x, y, z, angle, px, py, pz, qx, qy, qz)
         xyz(1,:) = x(:)
         xyz(2,:) = y(:)
         xyz(3,:) = z(:)
         deallocate (x, y, z)

      case (20)  ! "Surface area"

         if (prompt) then
            if (tri) then
               call tri_area (header, xyzf)
               write (luncrt, '(/, a, es16.8, a)') &
                   ' Total surface area:', header%surface_area, ' sq. units'
            else
               write (luncrt, '(/, a)') &
                  ' Surface area has not been implemented for volume datasets.'
            end if
         end if

      case (21)  ! "Surface area and enclosed volume", or "Solid volume"

         if (prompt) then
            if (tri) then
               yes = .false.
               call ready (luncrt, &
        'Interior pt. for volume calculation? [y|n; <cr> = data range avg.] ', &
                           lunkbd, yes, cr, eof)

               if (yes) then
                  write (luncrt, '(a)', advance='no') ' Enter interior x,y,z: '
                  read  (lunkbd, *) header%interior_point(:)
               end if

               call tri_volume (header, xyzf)

               write (luncrt, '(/, a, es16.8, a)') &
                  ' Surface area:   ', header%surface_area, ' sq. units', &
                  ' Enclosed volume:', header%enclosed_volume, ' cu. units'
            else

               call vol_volume (header, xyzf)

               write (luncrt, '(/, a, es16.8, a)') &
                  ' Solid volume:', header%solid_volume, ' cu. units'
            end if
         end if

      case (22)  ! "CM and Moments of Inertia"

         if (prompt) then
            if (tri) then
               call tri_moments_of_inertia (header, xyzf)
            else
               call vol_moments_of_inertia (header, xyzf)
            end if
         end if

      case (23)  ! "Apply rotation matrix R from tri_moments_of_inertia"

         if (prompt) then
            if (tri) then
               call tri_apply_rotation_R (header, xyzf)
            else
               call vol_apply_rotation_R (header, xyzf)
            end if
         end if

      case (24)  ! "Smooth a bumpy one-zone surface dataset"

         if (prompt) then
            if (tri) then
               data_range = max (header%xmax - header%xmin, &
                                 header%ymax - header%ymin, &
                                 header%zmax - header%zmin)
               percent_data_range = 1.
               fwhm = data_range * 0.01
               write (luncrt, '(4(/, a), es11.4)') &
                 ' The full width at half maximum of the Gaussian kernel', &
                 ' is ~2.35 sigma, in units of distance from the current pt.', &
                 ' Use small values for removing high frequency noise.', &
                 ' The default corresponds to 1% of the data range: ', fwhm

               call readr (luncrt, 'fwhm % [1]: ', &
                           lunkbd, percent_data_range, cr, eof)
               if (eof) go to 999
               fwhm = data_range * percent_data_range  * 0.01

               niter = 1
               call readi (luncrt, '# iterations [default is 1]: ', &
                           lunkbd, niter, cr, eof)
               if (eof) go to 999

               allocate (x(n), y(n), z(n))
               x(:) = xyz(1,:)
               y(:) = xyz(2,:)
               z(:) = xyz(3,:)
               call smoothxyz (n, x, y, z, fwhm, niter)
               xyz(1,:) = x(:)
               xyz(2,:) = y(:)
               xyz(3,:) = z(:)
               deallocate (x, y, z)

            else
               write (luncrt, '(/, a)') &
                  ' Smoothing of volume datasets is not an option.'
               go to 999
            end if
         end if

      end select

999   return

      end subroutine transform

   end program triangulation_tool
