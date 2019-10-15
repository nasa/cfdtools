!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   program triangulation_tool

!  Description:
!
!     Read one unstructured dataset (surface or volume), perform the indicated
!  analysis or analyses, and (if indicated) save the result as another such
!  file.  Initially, this serves as a driver for testing utilities added to the
!  triangulation_io module for calculating geometric properties of surface
!  datasets.  Now, it also drives the analogous utilities for volume datasets.
!  Most operations are on the x/y/z coordinates only, but one option for a
!  surface triangulation with cell-centered function values can interpolate
!  the functions to the vertices (area-weighted averaging).  Another option
!  merges multizone surfaces as one zone.  Both of these options are performed
!  within the tri_read routine of triangulation_io, so certain restrictions
!  apply (more than one operation in one run may be inappropriate).
!
!     The dataset may contain more than one zone (but often doesn't).
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
!  WARNING: the second set of face neighbor pointers in this format is NOT
!  handled and may be missing in any output file.
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
!     07/26/18   "   Made use of a new triangulation_io option to convert cell-
!                    centered function data to vertex-centered.  It also has an
!                    option to merge multi-zone surface triangulations into a
!                    single zone in new file header fields analogous to the
!                    relevant fields of a zone.  Any subsequent manipulations
!                    of x/y/z coordinates should be performed in a separate run.
!                    A third option to interpolate vertex function values to
!                    cell centroids has also been added.
!     07/27/18   "   Arranged to visualize the cell areas of a triangulation,
!                    which are needed for the CM and moments of inertia of menu
!                    option 22.  An additional surface dataset is written with a
!                    single function in favor of adding the areas to any
!                    existing functions (which may not be cell-centered).
!
!  Author:  David Saunders, ERC, Inc./NASA Ames Research Center, CA
!                Later with AMA, Inc. at NASA ARC.
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
      mxmenu = 27 + 1,   &   ! The + 1 counts the 99/Done choice
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
      xyzarea, xyzf

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
      '  25: Multizone surface --> one zone',  &
      '  26: Cell-centrd. surf. --> vertex-c.',&
      '  27: Vertex-centrd. surf --> centrds.',&
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
   header%combine_zones = tri  ! Do it during reading; may be all that is needed
   header%centroids_to_vertices = .false. ! Reread the file if this is requested

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

      else if (choice == 25) then  ! "Multizone surface --> single zone"

         go to 800 ! This has already been done during the read of all zones

      else if (choice == 26) then  ! "Cell-centered fns. to vertex-centered."

         if (.not. tri) then
            write (luncrt, '(/, a)') &
               ' This option has been implemented for triangulations only.'
            go to 99
         end if

!        Reread the file with this option turned on:

         header%centroids_to_vertices = .true.
         call tri_read (lunin, header, xyzf, ios)
         write (luncrt, '(a, i7)') ' # zones found:', header%nzones
         header%fileform = 1
         xyzf(:)%element_type = 'TRIANGLE'

      else if (choice == 27) then  ! "Vertex-centered fns. --> centroids"

         if (header%fileform == 2) then
            write (luncrt, '(/, a)') ' The functions are already cell-centered.'
            go to 210
         end if

         if (.not. tri) then
            write (luncrt, '(/, a)') &
               ' This option has been implemented for triangulations only.'
            go to 99
         end if

         do iz = 1, header%nzones
            call tri_zone_v2c (iz, header%numf, xyzf(iz))  ! Internal procedure
         end do
         header%fileform = 2

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
         else if (choice /= 25) then
            call tri_write (lunout, header, xyzf, ios)
         else  ! Choice = 25 is a special retrofitted case

!           The merged zones are not in a tri_type data structure.  Kludge it:

            call deallocate_tri_zones (1, header%nzones, header%numf, xyzf, ios)
            header%nzones     = 1
            xyzf(1)%nnodes    = header%nnodes
            xyzf(1)%nelements = header%nelements
            call tri_zone_allocate (header, xyzf(1), ios)
            if (ios /=  0) go to 99

            xyzf(1)%conn = header%conn
            xyzf(1)%xyz  = header%xyz
            xyzf(1)%f    = header%f
            
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

!     Internal procedures for program TRIANGULATION_TOOL:

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

      integer, intent (in)    :: iz        ! Zone number
      integer, intent (in)    :: n         ! # nodes in this zone
      real,    intent (inout) :: xyz(3,n)  ! Packed node coordinates, this zone

!     Local variables:

      integer    :: i, idp, iform, izone, l, niter, ne, nf, nn, nzones
      real       :: data_range, fwhm, percent_data_range, temp
      real, save :: angle, p, px, py, pz, q, qx, qy, qz, scale, shift
      logical    :: prompt
      real, allocatable, dimension (:) :: x, y, z
      character (32) :: f1name           ! For saving triangle areas
      character (80) :: header_filename  ! For temporary reuse of header

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

      case (22)  ! "CM and Moments of Inertia" + option to save triangle areas

         if (prompt) then
            if (tri) then
               call tri_moments_of_inertia (header, xyzf)  ! All zones are done

!              Option to visualize the triangle areas:

               yes = .true.
               call ready (luncrt, &
              'Save the triangle areas for visualization? [y|n; <cr> = yes] ', &
                           lunkbd, yes, cr, eof)

               if (yes) then  ! Least awkward to copy all zone info except f
                  nzones = header%nzones
                  allocate (xyzarea(nzones))
                  do izone = 1, nzones
                     xyzarea(izone)%zone_title   = xyzf(izone)%zone_title
                     xyzarea(izone)%element_type = xyzf(izone)%element_type
                     xyzarea(izone)%nelements    = xyzf(izone)%nelements
                     xyzarea(izone)%nnodes       = xyzf(izone)%nnodes
                     ne = xyzf(izone)%nelements
                     nn = xyzf(izone)%nnodes
                     allocate (xyzarea(izone)%conn(3,ne), &
                               xyzarea(izone)%xyz(3,nn),  &
                               xyzarea(izone)%f(1,ne))
                     xyzarea(izone)%conn(:,:) = xyzf(izone)%conn(:,:)
                     xyzarea(izone)%xyz(:,:)  = xyzf(izone)%xyz(:,:)
                     xyzarea(izone)%f(1,:)    = xyzf(izone)%area(:)
                  end do
                  header_filename = header%filename
                  l = len_trim (header_filename)
                  header%filename = header_filename(1:l-3) // 'areas.dat'
                  iform           = header%fileform
                  header%fileform = 2  ! Cell-centered
                  idp             = header%datapacking
               header%datapacking = 1  ! Block order
                  nf              = header%numf
                  header%numf     = 1
                  if (nf > 0) then
                     f1name = header%varname(4)
                  else
                     deallocate (header%varname);  allocate (header%varname(4))
                     header%varname(1) = 'x'
                     header%varname(2) = 'y'
                     header%varname(3) = 'z'
                  end if
                  header%varname(4) = 'Triangle area'

                  call tri_write (lunout, header, xyzarea, ios)
                  if (ios /= 0) write (luncrt, '(/, a)') &
                     ' Trouble saving triangulation with cell areas.'

                  header%filename = header_filename
                  header%fileform = iform
               header%datapacking = idp
                  header%numf     = nf
                  if (nf > 0) header%varname(4) = f1name
               end if
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

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine tri_zone_v2c (iz, nf, xyzf)  ! For one zone, f(vert.)-->f(cen.)

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer,         intent (in)    :: iz    ! Zone number
      integer,         intent (in)    :: nf    ! Number of functions
      type (tri_type), intent (inout) :: xyzf  ! Zone iz data structure

!     Local constants:

      real, parameter :: third = 1./3.

!     Local variables:

      integer :: i1, i2, i3, ic, nc
      real, allocatable :: fcentroid(:,:)

!     Execution:

      write (luncrt, '(a, i7)') ' Converting functions for zone', iz
      nc = xyzf%nelements
      allocate (fcentroid(nf,nc))

      do ic = 1, nc
         i1 = xyzf%conn(1,ic)
         i2 = xyzf%conn(2,ic)
         i3 = xyzf%conn(3,ic)
         fcentroid(:,ic) = (xyzf%f(:,i1) + xyzf%f(:,i2) + xyzf%f(:,i3))*third
      end do

      deallocate (xyzf%f);  allocate (xyzf%f(nf,nc))
      xyzf%f(:,:) = fcentroid(:,:)
      deallocate (fcentroid)

      end subroutine tri_zone_v2c

   end program triangulation_tool
