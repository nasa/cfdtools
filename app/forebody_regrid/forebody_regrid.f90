!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program forebody_regrid
!
!  Program HEAT_SHIELD produces a single-patch surface grid with spokes that are
!  uniformly distributed in the azimuthal direction and form a singular point at
!  the apex. This specialized application replaces the apex region with a quasi-
!  rectangular patch (or two) to avoid a singular axis in the associated volume
!  grid.  It reads the left half of the forebody as written by HEAT_SHIELD, but
!  results are written in the more usual form as the starboard (right) half with
!  Y positive and Z up.
!
!  Input Coordinate System (Right-Handed):
!
!     X points downstream for zero angle of attack
!     Y points up
!     Z points to the right looking from upstream
!
!  Output Coordinate System (Right-Handed):
!
!     X points downstream for zero angle of attack
!     Y points to the left looking from upstream
!     Z is up
!
!  This version has been extended for the specialized case of an umbrella-type
!  configuration, where some deflection needs to be applied to the conical
!  portion of the (assumed sphere-cone) geometry between the "spokes" to model
!  flexible thermal protection material.  This is done before regridding the
!  nose, and should be independent since the nose quadrant patches are expected
!  to be forward of the sphere-cone juncture.  The vertical centerline is
!  placed between umbrella spokes so that the symmetry plane shows the region
!  of maximum deflection.  The edges of the now-polygonal shoulder are morphed
!  to be straight lines in the azimuthal direction for any (semi)patch between
!  spokes, before the deflections are applied.  Because of the assumed symmetry,
!  the strategy is to regrid one deflected semipatch and transcribe it to all
!  the other locations.  The deflections are constructed from catenary curves.
!
!  Original Nose-Regridding Strategy:
!
!  >  Retain the existing outer portion of the grid, split into subpatches as
!     needed.
!  >  Some radial index defines the roughly circular perimeter of the new
!     patches in the nose region (not necessarily at constant X).  17 or 33
!     are the usual choices depending on grid dimensions.  These are heuristic
!     but are good choices for 101 or 201 points along the input grid spokes
!     and 65 or 121 spokes.  N.B.: the number of spokes minus 1 must be a
!     multiple of 4 for the 2-quadrant nose patch option.
!  >  An elliptically smoothed planar, normalized patch (quarter or half circle)
!     underlies the new patch(es) as a template.  (Thanks, Jim Brown.)
!  >  Its number of perimeter cells in the template is redistributed if
!     necessary to match the given surface grid.
!  >  This form is morphed to the desired perimeter (WARPQ3D, a TFI-like scheme
!     that employs the interior relative spacings).
!  >  The interior points are projected onto the original surface grid (ADT
!     scheme).
!
!     The implied form of surface interpolation has the potential to introduce
!     faceting, but that should be negligible for reasonable resolution in the
!     input surface grid.  Gridgen or equivalent appears the only alternative.
!     Densifying the input surface grid by a factor of 16 beforehand helps.
!
!  Control:
!
!     Answer the prompts.
!
!  History:
!
!     01/09/2008  D.A.Saunders   Initial design.
!     01/10/2008-   "    "       Initial implementation (quarter circle case).
!     01/14/2008
!     12/17/2010  DAS, ERC, Inc. Fudged the two poor corner points of the
!                                circular quadrants (only), after pondering
!                                squarer nose patches that aren't really
!                                doable (well) algebraically.
!     10/11/2011 -  "    "       Added the options to impose umbrella-type
!     10/14/2011                 shoulder edge faceting & catenary deflections.
!     10/18/2011    "    "       Added the option to round off the sharp ridges.
!     10/20/2011    "    "       Noticed method ='b' should be 'B' (loose fits).
!                                Added the option to resolve the rib ridges by
!                                restretching outboard of where the spacing is
!                                more than (half) the specified rib thickness.
!     11/03/2011    "    "       Actually, 'b' was as intended for ADJUSTN: it
!                                gives 'B' (loose fits) and uniform spacing.
!     03/18/2012    "    "       A capsule+sting grid, initiated with the full-
!                                body CAPSULE_GRID option then truncated, needed
!                                the revolved surface grid to look like a spoked
!                                surface from HEAT_SHIELD.  Therefore, if the
!                                input surface is not named 'heat_shield.xyz',
!                                swap y and z and reverse the j indices to make
!                                it match the original HEAT_SHIELD conventions.
!     02/26/2014    "    "       All ADT variants have been merged into a module
!                                with generic build & search calls.
!                                Having to rename a file as heat_shield.xyz to
!                                avoid changing its x/y/z convention is clunky.
!                                Therefore, prompt the user.
!     04/18/2014    "    "       Application to extreme ellipsoid shapes by
!                                John Theisinger prompted transcribing the
!                                CAPSULE_GRID changes that do not assume the
!                                template nose quadrant has uniform spacing,
!                                even though experience has shown that DPLR
!                                gives cleaner nose heat flux contours when
!                                the nose patching spacing is uniform.
!  Author:
!
!     David Saunders, ELORET Corporation/NASA Ames Research Center, CA
!                now: ERC Inc./NASA ARC
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use adt_utilities            ! All ADT build & search routines
   use grid_block_structure     ! Derived data type for one grid block
   use surface_patch_utilities  ! 3-space surface grid utilities
   use xyzq_io_module           ! PLOT3D file I/O package

   implicit none

!  Constants:

   integer, parameter :: &
      lunkbd = 5,        &      ! Keyboard inputs
      luncrt = 6,        &      ! Screen
      lunin1 = 1,        &      ! Input surface grid from HEAT_SHIELD
      lunin2 = 2,        &      ! Input normalized nose patch
      lunchk = 3,        &      ! For checking redistributed nose patch
      lunout = 4                ! Regridded surface for volume gridding

   logical, parameter :: &
      false  = .false.,  &
      true   = .true.

!  Variables:

   integer :: &
      ib, ios, nb, nblocks, nedges, ni_outer, ni_regrid, ni_template, &
      ni_spoke, nj_semi, nj_template, nj_spoke, nose_case, nj_quarter

   real :: &
      peak_deflection, rib_thickness, semiangle

   logical :: &
      cell_centered, resolve_the_ridges, round_the_ridges, umbrella_case

   character :: &
      filename * 80

   real, allocatable, dimension (:) :: &  ! For changing pt. counts along lines
      snorm, x1, y1, z1, x2, y2, z2

   type (grid_type), pointer, dimension (:) :: &
      template,   &             ! Normalized quarter or half circle, with x = 0.
      spoke_grid, &             ! Single-block, uniform azimuth, singular apex
      new_grid                  ! Desired regridded surface

!  Execution:
!  ::::::::::

   call control ()  ! Internal procedure below; prompt and set up in/out files

   call morph_to_umbrella ()    ! If requested in control

   call adjust_template ()      ! Match circumf. point count to initial grid

   call initial_regrid ()       ! Before projecting to desired nose surface

   call project_to_nose ()      ! ADT scheme

   call save_result ()          ! And clean up

   contains

!     Internal procedures for program forebody_regrid:
!     ::::::::::::::::::::::::::::::::::::::::::::::::

!     --------------------------------------------------------------------------
      subroutine control ()            ! Prompt and set up input & output files
!     --------------------------------------------------------------------------

      integer   :: numf
      character :: answer*1

!     Execution:

      ios = 1
      do while (ios /= 0)
         write (luncrt, '(/, a)', advance='no') &
            'Input single-block forebody spoked grid (formatted): '
         read  (lunkbd, '(a)') filename
         open  (lunin1, file=filename, status='old', iostat=ios)
      end do

      call xyzq_read (lunin1, -1, true, nb, numf, cell_centered, &
                      spoke_grid, ios)

      if (ios /= 0) then
         write (luncrt, '(/, a)') 'Trouble reading the spoke grid.'
         stop
      end if

      if (nb /= 1) then
         write (luncrt, '(/, a, i4)') 'Unexpected # spoked grid blocks:', nb
         stop
      end if

      ni_spoke   = spoke_grid(1)%ni    ! Radial, from apex
      nj_spoke   = spoke_grid(1)%nj    ! Azimuthal, from 6 o'clock position
      nj_quarter = (nj_spoke + 1) / 2  ! # spokes on 1/4th of (whole) forebody

      if (mod (nj_spoke - 1, 4) /= 0) then
         write (luncrt, '(/, a, i4)') &
            'The input number of spokes must be of the form 4m + 1: ', nj_spoke
         stop
      end if

!     The assumed xyz/ij convention is that of the HEAT_SHIELD program, which
!     outputs heat_shield.xyz.  If this is not the input file name, assume that
!     the forebody is the result of rotating a generatrix from CAPSULE_GRID.
!     Swapping y & z and reversing j will make it match HEAT_SHIELD conventions.
!     NO! Prompt the user rather than require the name to be heat_shield.xyz
!     to avoid the swapping/

      if (trim (filename) /= 'heat_shield.xyz') then
         write (luncrt, '(a)', advance='no') &
            'Does the input surface grid use HEAT_SHIELD''s x/y/z convention? '
         read  (lunkbd, '(a)') answer
         if (answer /= 'y' .and. answer /= 'Y') then
            write (luncrt, '(a)') &
            'Swapping y and z and reversing j to match HEAT_SHIELD conventions.'
            call swap_coordinates (spoke_grid(1), 2, 3)
            call reverse_patch_j  (spoke_grid(1), 0)
         end if
      end if

      ios = 1
      do while (ios /= 0)
         write (luncrt, '(a)', advance='no') &
            'Input normalized nose patch template (formatted):   '
         read  (lunkbd, '(a)') filename
         open  (lunin2, file=filename, status='old', iostat=ios)
      end do

      call xyzq_read (lunin2, -1, true, nb, numf, cell_centered, &
                      template, ios)

      if (ios /= 0) then
         write (luncrt, '(/, a)') 'Trouble reading normalized nose patch.'
         stop
      end if

      ni_template = template(1)%ni
      nj_template = template(1)%nj

      if (ni_template /= nj_template) then
         write (luncrt, '(/, a, 2i5)') &
           'Input nose patch dimensions don''t match:', ni_template, nj_template
         stop
      end if

      write (luncrt, '(/, (a))') &
         'Nose patch template options:', &
         '   1:  Semicircle, 2 "corners" on circle (inactive)', &
         '   2:  Quarter circle                      (active)', &
         '   3:  Semicircle, 1 "corner" on circle  (inactive)'

      write (luncrt, '(a)', advance='no') 'Case 1|2|3: '
      read  (lunkbd, *) nose_case

      write (luncrt, '(/, a)', advance='no') &
         'Radial spoke index i inside which nose is regridded: '
      read  (lunkbd, *) ni_regrid
      ni_outer = ni_spoke - ni_regrid + 1

      write (luncrt, '(/, a)', advance='no') &
         'Impose umbrella-type deflections? [y|n]: '
      read  (lunkbd, *) answer
      umbrella_case = answer == 'y' .or. answer == 'Y'

      if (umbrella_case) then
         write (luncrt, '(a)', advance='no') 'Number of full polygon edges: '
         read  (lunkbd, *) nedges
         nj_semi = 1 + (nj_spoke - 1) / nedges  ! # input grid spokes/semiangle
         if (nedges * (nj_semi - 1) + 1 /= nj_spoke) then
            write (luncrt, '(a, 2i4)') &
               '# grid spokes not suited to this many edges:', nj_spoke, nedges
            stop
         end if
         semiangle = 180. / real (nedges)      ! <-> half a polygon edge
         write (luncrt, '(a)', advance='no') 'Maximum deflection: '
         read  (lunkbd, *) peak_deflection
         write (luncrt, '(a)', advance='no') 'Resolve the rib ridges? '
         read  (lunkbd, *) answer
         resolve_the_ridges = answer == 'y' .or. answer == 'Y'
         if (resolve_the_ridges) then
            write (luncrt, '(a)', advance='no') 'Rib thickness to resolve: '
            read  (lunkbd, *) rib_thickness
         end if
         write (luncrt, '(a)', advance='no') 'Round off the sharp ridges? '
         read  (lunkbd, *) answer
         round_the_ridges = answer == 'y' .or. answer == 'Y'
      end if

      write (luncrt, '(/, a)', advance='no') &
         'Output regridded surface file name: '
      read  (lunkbd, '(a)') filename

      open  (lunout, file=filename, status='unknown')

      end subroutine control

!     -----------------------------------------------------------------------
      subroutine adjust_template ()  ! Match circumf. pt. count to input grid
!     -----------------------------------------------------------------------

      character, parameter :: method * 1 = 'B' ! Loose local cubic spline

      integer          :: i, j, new_n
      type (grid_type) :: new_patch(2)

!     Execution:

      select case (nose_case)

         case (1) ! Three edges along a semicircle

            write (luncrt, '(/, a)') 'This case has not been implemented.'
            stop

         case (2) ! Quarter circle with two edges on the circle

            if (2*ni_template - 1 /= nj_quarter) then  ! Redistribute i and j

               new_n = (nj_quarter + 1) / 2

               write (luncrt, '(/, a, i4, a, i4)') &
                  'Adjusting template size:', ni_template, ' -->', new_n

               new_patch(1)%ni = new_n
               new_patch(1)%nj = nj_template
               new_patch(1)%nk = 1

               call xyz_allocate (new_patch(1), ios)

               allocate (x1(ni_template), y1(ni_template), z1(ni_template),    &
                         x2(new_n),       y2(new_n),       z2(new_n))

               do j = 1, nj_template  ! Rows

                  do i = 1, ni_template
                     x1(i) = template(1)%x(i,j,1)
                     y1(i) = template(1)%y(i,j,1)
                     z1(i) = template(1)%z(i,j,1)
                  end do

                  call adjustn (ni_template, 1, ni_template, x1, y1, z1, &
                                1, new_n, x2, y2, z2, method)

                  do i = 1, new_n
                     new_patch(1)%x(i,j,1) = x2(i)
                     new_patch(1)%y(i,j,1) = y2(i)
                     new_patch(1)%z(i,j,1) = z2(i)
                  end do

               end do

               new_patch(2)%ni = new_n
               new_patch(2)%nj = new_n
               new_patch(2)%nk = 1

               call xyz_allocate (new_patch(2), ios)

               do i = 1, new_n  ! Columns of new rows

                  do j = 1, nj_template
                     x1(j) = new_patch(1)%x(i,j,1)
                     y1(j) = new_patch(1)%y(i,j,1)
                     z1(j) = new_patch(1)%z(i,j,1)
                  end do

                  call adjustn (nj_template, 1, nj_template, x1, y1, z1, &
                                1, new_n, x2, y2, z2, method)

                  do j = 1, new_n
                     new_patch(2)%x(i,j,1) = x2(j)
                     new_patch(2)%y(i,j,1) = y2(j)
                     new_patch(2)%z(i,j,1) = z2(j)
                  end do

               end do

               deallocate (template(1)%x, template(1)%y, template(1)%z)
               deallocate (x1, y1, z1, x2, y2, z2)

               template(1)%ni = new_n;  ni_template = new_n
               template(1)%nj = new_n;  nj_template = new_n

               call xyz_allocate (template(1), ios)

               template(1)%x = new_patch(2)%x
               template(1)%y = new_patch(2)%y
               template(1)%z = new_patch(2)%z

               do i = 1, 2
                  deallocate (new_patch(i)%x, new_patch(i)%y, new_patch(i)%z)
               end do

            end if

            open (lunchk, file='adjusted_template.g', status='unknown')
            call xyz_write (lunchk, true, nb, template, ios)
            close (lunchk)

         case (3) ! Two edges along a semicircle

            write (luncrt, '(/, a)') 'This case has not been implemented.'
            stop

      end select

      end subroutine adjust_template

!     --------------------------------------------------------------------------
      subroutine initial_regrid ()     ! Before projecting to true nose surface
!
!     Replace the nose portion of the spoked surface with patches to avoid the
!     singular point.  A planar template patch is morphed appropriately.
!     This version does not assume uniformity in the template patch.
!     See program morph_quadrant for generating nonuniform patches.
!     --------------------------------------------------------------------------

      real,      parameter :: half = 0.5
      character, parameter :: method * 1 = 'B' ! Loose fits in ADJUSTN2

      integer :: i, ii, j, jj
      real    :: total, xfudge, yfudge, zfudge, yscale, zscale, zshift
      real, allocatable, dimension (:,:,:,:) :: uvw
      type (grid_type) :: nose_patch           ! Scaling template helps WARPQ3D2

!     Execution:

!     A quarter circle with two edges on the circle is the only case handled.
!     Originally, the template was assumed to be uniform along all edges, but
!     only the circular edges must be so: allow for nonuniformity along the
!     straight edges in case tightening the spacing towards the apex helps.

!     HEAT_SHIELD has y pointing up and z pointing right from upstream.
!     Switch to the preferred convention of DPLR users:

      call swap_coordinates (spoke_grid(1), 2, 3)

!     Strategy (WARPQ3D):
!     Morph the template by setting up edges and filling the interior(s).
!     ni_template has already been adjusted to match the current point count
!     determined by nj_spoke.

      allocate (uvw(ni_template, nj_template, 1, 3))  ! Normalize arc lengths

      select case (nose_case)

         case (1) ! Three edges along a semicircle

            ! Complete if needed

         case (2) ! Quarter circle with two edges on the circle

            nblocks = 6  ! # output blocks

            allocate (new_grid(nblocks))  ! 1-2 = upper & lower quarters of nose
                                          ! 3-6 = split outer grid, low to high

            new_grid(1:2)%ni = ni_template; new_grid(3:nblocks)%ni = ni_outer
            new_grid(1:2)%nj = nj_template; new_grid(3:nblocks)%nj = nj_template
            new_grid(:)%nk   = 1

            do ib = 1, nblocks
               call xyz_allocate (new_grid(ib), ios)
            end do

!           Shifting and scaling the template helps WARPQ3D behavior.  The
!           template is assumed to be in quadrant 1 of the unit (y,z) circle.

            nose_patch%ni = template(1)%ni
            nose_patch%nj = template(1)%nj
            nose_patch%nk = 1

            call xyz_allocate (nose_patch, ios)

!           Upper quarter of nose:
!           !!!!!!!!!!!!!!!!!!!!!!

!           The lower quarter was coded first & is clearer because the spokes
!           are assumed to start at 6 o'clock & go clockwise (from upstream).

!           Transcribe the (roughly) circular edge from 9 o'clock to 10:30:

            jj = 1
            do j = nj_quarter, nj_quarter + nj_template - 1
               new_grid(1)%x(ni_template,jj,1) = spoke_grid(1)%x(ni_regrid,j,1)
               new_grid(1)%y(ni_template,jj,1) = spoke_grid(1)%y(ni_regrid,j,1)
               new_grid(1)%z(ni_template,jj,1) = spoke_grid(1)%z(ni_regrid,j,1)
               jj = jj + 1
            end do

!           Edge from 10:30 to 12 o'clock:

            jj = nj_template
            do j = nj_quarter + nj_template - 1, nj_spoke
               new_grid(1)%x(jj,nj_template,1) = spoke_grid(1)%x(ni_regrid,j,1)
               new_grid(1)%y(jj,nj_template,1) = spoke_grid(1)%y(ni_regrid,j,1)
               new_grid(1)%z(jj,nj_template,1) = spoke_grid(1)%z(ni_regrid,j,1)
               jj = jj - 1
            end do

!           The horizontal edge of new patch 1 is part of the middle spoke,
!           but we need to adjust the point count and relative spacing:

            allocate (x1(ni_regrid+1), y1(ni_regrid+1), z1(ni_regrid+1),  &
                      x2(ni_template), y2(ni_template), z2(ni_template),  &
                      snorm(ni_template))

            do i = 1, ni_regrid + 1
               x1(i) = spoke_grid(1)%x(i,nj_quarter,1)
               y1(i) = spoke_grid(1)%y(i,nj_quarter,1)
               z1(i) = spoke_grid(1)%z(i,nj_quarter,1)
            end do

            yscale = y1(ni_regrid)

!           Impose the (adjusted) template's relative horizontal edge spacing:

            do i = 1, ni_template
               x2(i) = template(1)%x(i,1,1)
               y2(i) = template(1)%y(i,1,1)
               z2(i) = template(1)%z(i,1,1)
            end do

            call chords3d (ni_template, x2, y2, z2, true, total, snorm)

!           Use data point ni_regrid + 1, but there's no need to go beyond that.
!           x2/y2/z2 are reused.

            call adjustn2 (ni_regrid + 1, 1, ni_regrid, x1, y1, z1, 1, &
                           ni_template, snorm, x2, y2, z2, method)

            do i = 1, ni_template  ! This edge is common to new patch 2
               new_grid(1)%x(i,1,1) = x2(i);  new_grid(2)%x(i,1,1) = x2(i)
               new_grid(1)%y(i,1,1) = y2(i);  new_grid(2)%y(i,1,1) = y2(i)
               new_grid(1)%z(i,1,1) = z2(i);  new_grid(2)%z(i,1,1) = z2(i)
            end do

!           The vertical edge of new patch 1 is part of the last spoke,
!           but we need to adjust the point count and relative spacing:

            do i = 1, ni_regrid + 1
               x1(i) = spoke_grid(1)%x(i,nj_spoke,1)
               y1(i) = spoke_grid(1)%y(i,nj_spoke,1)
               z1(i) = spoke_grid(1)%z(i,nj_spoke,1)
            end do

            zshift = z1(1)
            zscale = z1(ni_regrid) - zshift

!           Impose the (adjusted) template's relative vertical edge spacing,
!           assumed to be the same as for its horizontal edge:

            call adjustn2 (ni_regrid + 1, 1, ni_regrid, x1, y1, z1,  &
                           1, nj_template, snorm, x2, y2, z2, method)

            do j = 1, nj_template
               new_grid(1)%x(1,j,1) = x2(j)
               new_grid(1)%y(1,j,1) = y2(j)
               new_grid(1)%z(1,j,1) = z2(j)
            end do

            nose_patch%x = template(1)%x
            nose_patch%y = template(1)%y * yscale
            nose_patch%z = template(1)%z * zscale + zshift

!           Normalized arc lengths for the scaled (planar) template nose patch:

            call paramxyz (1, ni_template, 1, nj_template, 1, 1, &
                           1, ni_template, 1, nj_template, 1, 1, &
                           nose_patch%x, nose_patch%y, nose_patch%z, uvw)

!!!         write (10, *) 1
!!!         write (10, *) ni_template, nj_template, 1
!!!         write (10, *) uvw(:,:,1,1), uvw(:,:,1,2), uvw(:,:,1,3)

!           Morph the scaled template interior pts. to the interim new interior:

            call warpq3d2 (1, ni_template, 1, nj_template, 1, 1, &
                           1, ni_template, 1, nj_template, 1, 1, &
                           nose_patch%x, nose_patch%y, nose_patch%z, uvw, &
                           new_grid(1)%x, new_grid(1)%y, new_grid(1)%z)

            call transpose_patch (new_grid(1), 0)  ! Only patch 2 is rt.-handed

!           Lower quarter of nose:
!           !!!!!!!!!!!!!!!!!!!!!!

            do j = 1, ni_template  ! Circle edge from 6 o'clock to 7:30
               new_grid(2)%x(j,nj_template,1) = spoke_grid(1)%x(ni_regrid,j,1)
               new_grid(2)%y(j,nj_template,1) = spoke_grid(1)%y(ni_regrid,j,1)
               new_grid(2)%z(j,nj_template,1) = spoke_grid(1)%z(ni_regrid,j,1)
            end do

            jj = nj_template
            do j = ni_template, nj_quarter  ! Edge from 7:30 to 9 o'clock
               new_grid(2)%x(ni_template,jj,1) = spoke_grid(1)%x(ni_regrid,j,1)
               new_grid(2)%y(ni_template,jj,1) = spoke_grid(1)%y(ni_regrid,j,1)
               new_grid(2)%z(ni_template,jj,1) = spoke_grid(1)%z(ni_regrid,j,1)
               jj = jj - 1
            end do

!           The vertical edge of new patch 2 is part of the first spoke,
!           but we need to adjust the point count and spacing:

            do i = 1, ni_regrid + 1
               x1(i) = spoke_grid(1)%x(i,1,1)
               y1(i) = spoke_grid(1)%y(i,1,1)
               z1(i) = spoke_grid(1)%z(i,1,1)
            end do

            zscale = z1(ni_regrid) - zshift

            call adjustn2 (ni_regrid + 1, 1, ni_regrid, x1, y1, z1,  &
                           1, nj_template, snorm, x2, y2, z2, method)

            do j = 1, nj_template
               new_grid(2)%x(1,j,1) = x2(j)
               new_grid(2)%y(1,j,1) = y2(j)
               new_grid(2)%z(1,j,1) = z2(j)
            end do

!           The horizontal edge of patch 2 is that of patch 1, copied above.

!!!         nose_patch%x = template(1)%x                   ! Same as above
!!!         nose_patch%y = template(1)%y * yscale          !   "   "   "
            nose_patch%z = template(1)%z * zscale + zshift

            call paramxyz (1, ni_template, 1, nj_template, 1, 1, &
                           1, ni_template, 1, nj_template, 1, 1, &
                           nose_patch%x, nose_patch%y, nose_patch%z, uvw)

!!!         write (11, *) 1
!!!         write (11, *) ni_template, nj_template, 1
!!!         write (11, *) uvw(:,:,1,1), uvw(:,:,1,2), uvw(:,:,1,3)

!           Morph the template interior points to the interim new interior:

            call warpq3d2 (1, ni_template, 1, nj_template, 1, 1, &
                           1, ni_template, 1, nj_template, 1, 1, &
                           nose_patch%x, nose_patch%y, nose_patch%z, uvw, &
                           new_grid(2)%x, new_grid(2)%y, new_grid(2)%z)

            deallocate (x1, y1, z1, x2, y2, z2, snorm)

!           Transcribe the outer grid, split into new patches 3-6:
!           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            do j = 1, ni_template
               ii = 0
               do i = ni_regrid, ni_spoke
                  ii = ii + 1
                  new_grid(3)%x(ii,j,1) = spoke_grid(1)%x(i,j,1)
                  new_grid(3)%y(ii,j,1) = spoke_grid(1)%y(i,j,1)
                  new_grid(3)%z(ii,j,1) = spoke_grid(1)%z(i,j,1)
               end do
            end do

            jj = 0
            do j = ni_template, nj_quarter
               ii = 0
               jj = jj + 1
               do i = ni_regrid, ni_spoke
                  ii = ii + 1
                  new_grid(4)%x(ii,jj,1) = spoke_grid(1)%x(i,j,1)
                  new_grid(4)%y(ii,jj,1) = spoke_grid(1)%y(i,j,1)
                  new_grid(4)%z(ii,jj,1) = spoke_grid(1)%z(i,j,1)
               end do
            end do

            jj = 0
            do j = nj_quarter, nj_quarter + nj_template - 1
               ii = 0
               jj = jj + 1
               do i = ni_regrid, ni_spoke
                  ii = ii + 1
                  new_grid(5)%x(ii,jj,1) = spoke_grid(1)%x(i,j,1)
                  new_grid(5)%y(ii,jj,1) = spoke_grid(1)%y(i,j,1)
                  new_grid(5)%z(ii,jj,1) = spoke_grid(1)%z(i,j,1)
               end do
            end do

            jj = 0
            do j = nj_spoke - ni_template + 1, nj_spoke
               ii = 0
               jj = jj + 1
               do i = ni_regrid, ni_spoke
                  ii = ii + 1
                  new_grid(6)%x(ii,jj,1) = spoke_grid(1)%x(i,j,1)
                  new_grid(6)%y(ii,jj,1) = spoke_grid(1)%y(i,j,1)
                  new_grid(6)%z(ii,jj,1) = spoke_grid(1)%z(i,j,1)
               end do
            end do

            do ib = 3, 6
               call transpose_patch (new_grid(ib), 0)
            end do

!           Fudge the corner points a little to avoid the nearly 180-degree
!           angles by averaging the first and second points on the outer spokes:

            xfudge = (new_grid(6)%x(1,1,1) + new_grid(6)%x(1,2,1))*half
            yfudge = (new_grid(6)%y(1,1,1) + new_grid(6)%y(1,2,1))*half
            zfudge = (new_grid(6)%z(1,1,1) + new_grid(6)%z(1,2,1))*half

            new_grid(1)%x(ni_template,ni_template,1) = xfudge
            new_grid(1)%y(ni_template,ni_template,1) = yfudge
            new_grid(1)%z(ni_template,ni_template,1) = zfudge

            new_grid(5)%x(ni_template,1,1) = xfudge
            new_grid(5)%y(ni_template,1,1) = yfudge
            new_grid(5)%z(ni_template,1,1) = zfudge

            new_grid(6)%x(1,1,1) = xfudge
            new_grid(6)%y(1,1,1) = yfudge
            new_grid(6)%z(1,1,1) = zfudge

            xfudge = (new_grid(4)%x(1,1,1) + new_grid(4)%x(1,2,1))*half
            yfudge = (new_grid(4)%y(1,1,1) + new_grid(4)%y(1,2,1))*half
            zfudge = (new_grid(4)%z(1,1,1) + new_grid(4)%z(1,2,1))*half

            new_grid(2)%x(ni_template,ni_template,1) = xfudge
            new_grid(2)%y(ni_template,ni_template,1) = yfudge
            new_grid(2)%z(ni_template,ni_template,1) = zfudge

            new_grid(3)%x(ni_template,1,1) = xfudge
            new_grid(3)%y(ni_template,1,1) = yfudge
            new_grid(3)%z(ni_template,1,1) = zfudge

            new_grid(4)%x(1,1,1) = xfudge
            new_grid(4)%y(1,1,1) = yfudge
            new_grid(4)%z(1,1,1) = zfudge

         case (3) ! Two edges along a semicircle

            ! Complete if needed

      end select

      deallocate (uvw)

      do ib = 1, nb
         deallocate (template(ib)%x, template(ib)%y, template(ib)%z)
      end do

      deallocate (nose_patch%x, nose_patch%y, nose_patch%z)

!!!   call xyz_write (lunout, true, nblocks, new_grid, ios)
!!!   stop

      end subroutine initial_regrid

!     --------------------------------------------------------------------------
      subroutine morph_to_umbrella ()  ! Specialized regridding of cone/shoulder
!     --------------------------------------------------------------------------

!     Local constants:

      real,      parameter :: fraction = 0.50, half = 0.5, one = 1.0, zero = 0.0
      logical,   parameter :: false = .false., true = .true.
      character, parameter :: method*1 = 'B'  ! Loose local cubic interpolation

!     Local variables:

      integer :: i, ier, ieval, ir, it1, it2, itm, j, je1, jo1, je2, jo2, &
                 m, n, nicat, njsemim1, npt
      logical :: new
      real    :: a, angle, curvature, ds, peak_local_deflection, &
                 px, py, pz, qx, r, semispan, si, smid, stotal, sutotal, &
                 tol, xr, yr, zr, &
                 derivs(3), sridge(4), xridge(4), yridge(4), zridge(4), &
                 unit_normal(3)
      real, allocatable, dimension (:)   :: s, snew, su, xs, xss, ys, yss, kappa
      real, allocatable, dimension (:)   :: xcat, ycat, zcat, &
                                            xicat,yicat,zicat,&
                                            xu,   yu,   zu
      real, allocatable, dimension (:,:) :: xold, yold, zold, &
                                            xnew, ynew, znew, &
                                            xtmp, ytmp, ztmp
!     Execution:

      if (.not. umbrella_case) return

!     First main step of two:  Morph the shoulder from rounded to polygonal.
!     ----------------------------------------------------------------------

!     Set up a half-edge segment as before-and-after work-space.  This will
!     become the first (in j) portion starting at 6 o'clock on the port side.
!     The symmetry plane is deliberately down the middle of bottom and top
!     segments so as to show the maximum deflection slice in the CFD solutions.

      n = ni_spoke;  m = nj_semi;  njsemim1 = m - 1

      allocate (xold(n,m), yold(n,m), zold(n,m), &
                xnew(n,m), ynew(n,m), znew(n,m))

      xold(:,:) = spoke_grid(1)%x(:,1:m,1);  xnew = xold
      yold(:,:) = spoke_grid(1)%y(:,1:m,1);  ynew = yold
      zold(:,:) = spoke_grid(1)%z(:,1:m,1);  znew = zold

!     Identify the first and last spoke indices on the conical flank:

      allocate (s(n), xs(n), xss(n), ys(n), yss(n), kappa(n))

      call chords2d (n, xold(:,1), yold(:,1), false, stotal, s)

      do i = 2, n - 1
         call fdcntr (i, s, xold(:,1), xs(i), xss(i))
         call fdcntr (i, s, yold(:,1), ys(i), yss(i))
      end do

      call curv2d (n - 2, xs(2), xss(2), ys(2), yss(2), kappa(2))

      deallocate (xs, xss, ys, yss)

      m   = n / 2  ! Roughly mid-spoke
      it1 = 0
      it2 = 0
      curvature = kappa(4)  ! Nose value
      tol = abs (curvature) * fraction

      do i = m, 2, -1
         if (abs (kappa(i) - curvature) < tol) then
            it1 = i + 2
            exit
         end if
      end do

      do i = m, n  ! Find the shoulder curvature
         curvature = max (curvature, abs (kappa(i)))
      end do

      tol = curvature * fraction
      do i = m, n - 1
         if (abs (abs (kappa(i)) - curvature) < tol) then
            it2 = i - 2
            exit
         end if
      end do

      deallocate (kappa)

      write (luncrt, '(a, 2i4)') &
         'First & last spoke indices bounding deflections:', it1, it2
      if (it1 * it2 == 0) stop

!     Substitute straight lines from i = it2:ni_spoke and adjust the rest of
!     the spokes from i = it1:it2.  Just do it for one semiedge patch:

      do i = it2, ni_spoke
         xnew(i,1) = xold(i,nj_semi)
         ynew(i,1) = yold(i,nj_semi)
         znew(i,1) = zero
         do j = 2, nj_semi - 1
            r = real (j - 1) / real (nj_semi - 1)
            xnew(i,j) = xnew(i,1)
            ynew(i,j) = ynew(i,1)
            znew(i,j) = zold(i,nj_semi)*r
         end do
      end do

!     Adjust the spoke distribution between i = it1 and it2:

      do j = 1, nj_semi - 1
         call nuline3d (it1, it2, xold(:,j), yold(:,j), zold(:,j), &
                                  xnew(:,j), ynew(:,j), znew(:,j))
      end do

!!    write (48, '(i1)')  1
!!    write (48, '(3i4)') ni_spoke, nj_semi, 1
!!    write (48, '(1p, 8e15.7)') &
!!       xnew(:,1:nj_semi), ynew(:,1:nj_semi), znew(:,1:nj_semi)

      derivs(1) = -999.  ! Suppresses unwanted calculations

!     Impose the specified deflection on the working half-edge segment:
!     -----------------------------------------------------------------

      if (peak_deflection /= zero) then

         xold = xnew;  yold = ynew;  zold = znew  ! x/y/zold are undeflected

         call chords2d (n, xold(:,1), yold(:,1), false, stotal, s)

         smid = half * stotal    ! Locate the index nearest to mid-it1-it2 arc;
         itm  = (it1 + it2) / 2  ! arcs for the 6 o'clock spoke were done above

         call interval (ni_spoke, s, smid, one, itm)

         if (smid - s(itm) > s(itm+1) - smid) itm = itm + 1

!        Construct a uniform-in-spanwise-coordinate catenary (z,y) grid (right
!        half) for spanning the first half segment at i = itm:

         call point_A_to_point_B (ni_spoke, nj_semi, itm, 1, itm, nj_semi, &
                                  xold, yold, zold, semispan)

         allocate (xcat(nj_semi), ycat(nj_semi), zcat(nj_semi))

!        2-space form of the right-half catenary.  Since we only need this for
!        the peak deflections in the i direction, set the right end y to be 0.:

         call catenary_grid (semispan, peak_deflection, zero, 2, &
                             nj_semi, zcat, ycat, a, ier)
         if (ier /= 0) stop

         ycat(:) = -ycat(:)  ! Invert the droop
!!       write (9, '(3f12.7)') (zero, ycat(i), zcat(i), i = 1, nj_semi)

!        Uniform catenary spacing is appropriate here for the spanwise points
!        defining the maximum deflections along each new spoke of this patch.
!        Now we impose deflections along the i direction from it1 to it2.

         nicat = it2 - it1 + 1

         allocate (xicat(nicat), yicat(nicat), zicat(nicat))
         allocate (xu(nicat), yu(nicat), zu(nicat), su(nicat))

         do j = 1, nj_semi - 1

!           Construct a uniform catenary grid in (y,x) space (both halves)
!           spanning the distance between points (it1,j) and (it2,j):

            call point_A_to_point_B (ni_spoke, nj_semi, it1, j, it2, j, &
                                     xold, yold, zold, semispan)

            semispan = half * semispan
            peak_local_deflection = ycat(j)
            call catenary_grid (semispan, peak_local_deflection, xold(it1,j), &
                                3, nicat, yicat, xicat, a, ier)
            if (ier /= 0) stop

!!          write (9, '(a, i3)') 'j:', j
!!          write (9, '(3f12.7)') (xicat(i), yicat(i), zero, i = 1, nicat)

!           Flip it because we're working on the lower (6'oclock) half-segment:

            xicat(:) = 2.*xold(it1,j) - xicat(:)
            zicat(:) = zold(it1,j)

!!          write (9, '(a)') 'Flipped'
!!          write (9, '(3f12.7)') (xicat(i), yicat(i), zicat(i), i = 1, nicat)

!           Reverse the abscissas because they decrease with i on this patch:

            call rverse (nicat, xicat, xicat)
            call rverse (nicat, yicat, yicat)

!           Transform the 2-space catenary to 3-space:

            xu(1) = xold(it1,j)
            yu(1) = yold(it1,j)
            zu(1) = zold(it1,j)
            xu(nicat) = xold(it2,j)
            yu(nicat) = yold(it2,j)
            zu(nicat) = zold(it2,j)

            call rigid_transform (nicat, xicat, yicat, zicat, xu, yu, zu)

!!          write (9, '(a)') 'Transformed'
!!          write (9, '(3f12.7)') (xu(i), yu(i), zu(i), i = 1, nicat)

!           Replace the uniform spacing with the original relative spacing:

            call chords3d (nicat, xu, yu, zu, false, sutotal, su)
            call chords3d (n, xold(:,j), yold(:,j), zold(:,j), false, stotal, s)

            r = sutotal / (s(it2) - s(it1))
            new = true
            ieval = 1

            do i = it1 + 1, it2 - 1
               si = (s(i) - s(it1)) * r
               call plscrv3d (nicat, xu, yu, zu, su, method, new, false, si, &
                              ieval, px, py, pz, derivs)
               xnew(i,j) = px
               ynew(i,j) = py
               znew(i,j) = pz
               new = false
            end do

         end do  ! Next spoke-wise deflection

!!       write (51, '(i1)')  1
!!       write (51, '(3i4)') ni_spoke, nj_semi, 1
!!       write (51, '(1p, 8e15.7)') &
!!          xnew(:,1:nj_semi), ynew(:,1:nj_semi), znew(:,1:nj_semi)

         deallocate (xcat, ycat, zcat, xicat, yicat, zicat, xu, yu, zu, su)

      end if  ! Deflection option

!     Resolve the ridges?  If so, we find the first i station where the rib
!     spacing exceeds half the specified rib width, and restretch from there.

      if (resolve_the_ridges) then

         n = ni_spoke;  m = nj_semi;  i = ni_regrid;  ir = n + 1

         call point_A_to_point_B (n, m, i, m-1, i, m, xold, yold, zold, ds)

         rib_thickness = max (rib_thickness, ds + ds)  ! Avoid the new nose cap

         do i = ni_regrid + 1, ni_spoke

            call point_A_to_point_B (n, m, i, m-1, i, m, xold, yold, zold, ds)

            if (ds > half*rib_thickness) then
                ds = half*rib_thickness
                ir = i
                exit
            end if

         end do

         if (ir <= ni_spoke) &
            write (luncrt, '(a, i4)') 'Resolving rib thickness from i =', ir

         allocate (xu(m), yu(m), zu(m), su(m), snew(m))

         do i = ir, ni_spoke  ! Impose a Vinokur distrbn. at each i

            xu(:) = xnew(i,:);  yu(:) = ynew(i,:);  zu(:) = znew(i,:)

            call chords3d (m, xu, yu, zu, false, sutotal, su)

            call htdis4 (true, zero, sutotal, su(2), ds, m, snew, -luncrt, ier)

            new = true;  ieval = 1

            do j = 1, m
               call plscrv3d (m, xu, yu, zu, su, method, new, false, snew(j), &
                              ieval, px, py, pz, derivs)
               xnew(i,j) = px
               ynew(i,j) = py
               znew(i,j) = pz
               new = false
            end do

         end do

         deallocate (xu, yu, zu, su, snew)

      end if  ! Azimuthal-wise redistribution option to resolve ridges
      
!     Smooth out the sharp edges at the ridges?  If so, we rotate so that the
!     ridge is at 6 o'clock, and use symmetric 4-point cubics across the ridge
!     (omitting the ridge point) to interpolate at the ridge location, then
!     rotate back.

      px = xold(1,1);  qx = px - one  ! Rotation axis definition
      py = yold(1,1)
      pz = zold(1,1)

      if (round_the_ridges) then

!        Work with the last 3 radial lines of the half-edge patch:

         n = ni_spoke;  m = nj_semi;  npt = n*3
         allocate (xtmp(n,3), ytmp(n,3), ztmp(n,3))

         xtmp(:,:) = xnew(:,m-2:m)  ! Rotate the last 3 lines towards 6 o'clock
         ytmp(:,:) = ynew(:,m-2:m)
         ztmp(:,:) = znew(:,m-2:m);            angle = -semiangle

         call rotate3d (npt, xtmp, ytmp, ztmp, angle, px, py, pz, qx, py, pz)

         ieval = 2

         do i = it1 + 1, ni_spoke

            xridge(1:2) = xtmp(i,1:2);  xridge(3:4) =  xtmp(i,2:1:-1)
            yridge(1:2) = ytmp(i,1:2);  yridge(3:4) =  ytmp(i,2:1:-1)
            zridge(1:2) = ztmp(i,1:2);  zridge(3:4) = -ztmp(i,2:1:-1)

            call chords3d (4, xridge, yridge, zridge, false, stotal, sridge)

            si = half * stotal  ! Corresponds to ridge location

            call plscrv3d (4, xridge, yridge, zridge, sridge, method, true, &
                           false, si, ieval, xr, yr, zr, derivs)

            call rotate3d (1, xr, yr, zr, semiangle, px, py, pz, qx, py, pz)

            xnew(i,m) = xr
            ynew(i,m) = yr
            znew(i,m) = zr

         end do

         deallocate (xtmp, ytmp, ztmp)

      end if  ! Sharp edge smoothing option

!     Reflect this half-edge patch for transcribing, reusing x/y/zold(:,:):

      xold(it1+1:ni_spoke,nj_semi:1:-1) =  xnew(it1+1:ni_spoke,1:nj_semi)
      yold(it1+1:ni_spoke,nj_semi:1:-1) =  ynew(it1+1:ni_spoke,1:nj_semi)
      zold(it1+1:ni_spoke,nj_semi:1:-1) = -znew(it1+1:ni_spoke,1:nj_semi)

!     Rotate and transcribe in odd-even pairs:
!     ----------------------------------------

      npt = ni_spoke * nj_semi  ! # points in a half-edge patch
      jo1 = 1                   ! Indices for first odd half-edge patch
      jo2 = nj_semi 
      je1 = jo2                 ! Indices for first even half-edge patch
      je2 = je1 + njsemim1
      m   = nj_semi

      allocate (xtmp(n,m), ytmp(n,m), ztmp(n,m)) ! Because of in-place rotations

      angle = zero

      do n = 1, nedges

         if (mod (n, 2) == 1) then
            xtmp = xnew;  ytmp = ynew;  ztmp = znew
         else
            angle = semiangle * real (n)
            xtmp = xold;  ytmp = yold;  ztmp = zold
         end if

!!       write (6, '(a, i3, f6.1, 2i4,1x,2i4)') 'n,angle,jo1,jo2,je1,je2:', &
!!                                               n,angle,jo1,jo2,je1,je2

         if (angle /= zero) call rotate3d (npt, xtmp, ytmp, ztmp, angle, &
                                           px, py, pz, qx, py, pz)

         if (mod (n, 2) == 1) then
            spoke_grid(1)%x(it1+1:ni_spoke,jo1:jo2,1) = &
                       xtmp(it1+1:ni_spoke,1:nj_semi)
            spoke_grid(1)%y(it1+1:ni_spoke,jo1:jo2,1) = &
                       ytmp(it1+1:ni_spoke,1:nj_semi)
            spoke_grid(1)%z(it1+1:ni_spoke,jo1:jo2,1) = &
                       ztmp(it1+1:ni_spoke,1:nj_semi)
            jo1 = jo1 + 2*njsemim1
            jo2 = jo2 + 2*njsemim1
         else
            spoke_grid(1)%x(it1+1:ni_spoke,je1:je2,1) = &
                       xtmp(it1+1:ni_spoke,1:nj_semi)
            spoke_grid(1)%y(it1+1:ni_spoke,je1:je2,1) = &
                       ytmp(it1+1:ni_spoke,1:nj_semi)
            spoke_grid(1)%z(it1+1:ni_spoke,je1:je2,1) = &
                       ztmp(it1+1:ni_spoke,1:nj_semi)
            je1 = je1 + 2*njsemim1
            je2 = je2 + 2*njsemim1
         end if

      end do

!!    write (50, '(i1)')  1
!!    write (50, '(3i4)') ni_spoke, nj_spoke, 1
!!    write (50, '(1p, 8e15.7)') &
!!       spoke_grid(1)%x, spoke_grid(1)%y, spoke_grid(1)%z

      deallocate (xold, yold, zold, xnew, ynew, znew, xtmp, ytmp, ztmp)

      end subroutine morph_to_umbrella

!     --------------------------------------------------------------------------
      subroutine point_A_to_point_B (ni, nj, ia, ja, ib, jb, x, y, z, distance)
!
!     Ad hoc utility, named to be in alphabetical order here, for calculating
!     the distance between two points of the same surface grid patch.
!     --------------------------------------------------------------------------

      integer, intent (in)  :: ni, nj  ! Grid patch dimensions
      integer, intent (in)  :: ia, ja, ib, jb  ! Pair of surface points
      real,    intent (in)  :: x(ni,nj), y(ni,nj), z(ni,nj)  ! Grid coordinates
      real,    intent (out) :: distance  ! ... between the indicated points

      distance = sqrt ((x(ia,ja) - x(ib,jb))**2 + (y(ia,ja) - y(ib,jb))**2 + &
                       (z(ia,ja) - z(ib,jb))**2)

      end subroutine point_A_to_point_B

!     --------------------------------------------------------------------------
      subroutine project_to_nose ()    ! ADT scheme; densify nose spokes first
!     --------------------------------------------------------------------------

      integer,   parameter :: ntimes  = 16     ! Densification of nose spokes
      logical,   parameter :: densify = .true.
      real,      parameter :: tolerance = 1.e-9, zero = 0.
      character, parameter :: method * 1 = 'b' ! Loose fits/uniform spacing in
                                               ! ADJUSTN usage here

      integer :: i, ib1, ib2, iquad, j, nb, ni_dense, nj_data, nj_dense, npts, &
                 nquad
      real    :: dmax, dmean, dmin, dsqmin, p, q
      real    :: interp_xyz(3), target_xyz(3)
      type (grid_type) :: dense_spokes(1)  ! Generic interfaces expect an array

      integer, allocatable :: conn(:,:)    ! Search tree connectivity info.

!     Execution:

!     Densify the nose portion of the spokes to reduce faceting?

      if (densify) then

         ni_dense =  ni_regrid
         nj_dense = (nj_spoke - 1) * ntimes + 1
         nj_data  =  nj_spoke + 2  ! Spokes are periodic

         dense_spokes(1)%ni = ni_dense
         dense_spokes(1)%nj = nj_dense
         dense_spokes(1)%nk = 1

         call xyz_allocate (dense_spokes(1), ios)

         allocate (x1(nj_data),  y1(nj_data),  z1(nj_data),  &
                   x2(nj_dense), y2(nj_dense), z2(nj_dense))

         dense_spokes(1)%x(1,:,1) = spoke_grid(1)%x(1,1,1)
         dense_spokes(1)%y(1,:,1) = spoke_grid(1)%y(1,1,1)
         dense_spokes(1)%z(1,:,1) = spoke_grid(1)%z(1,1,1)

         do i = 2, ni_regrid
            x1(1)           =  spoke_grid(1)%x(i,2,1)         ! Periodic data
            x1(2:nj_data-1) =  spoke_grid(1)%x(i,:,1)
            x1(nj_data)     =  spoke_grid(1)%x(i,nj_spoke-1,1)
            y1(1)           = -spoke_grid(1)%y(i,2,1)
            y1(2:nj_data-1) =  spoke_grid(1)%y(i,:,1)
            y1(nj_data)     = -spoke_grid(1)%y(i,nj_spoke-1,1)
            z1(1)           =  spoke_grid(1)%z(i,2,1)
            z1(2:nj_data-1) =  spoke_grid(1)%z(i,:,1)
            z1(nj_data)     =  spoke_grid(1)%z(i,nj_spoke-1,1)

            call adjustn (nj_data, 2, nj_data - 1, x1, y1, z1, &
                          1, nj_dense, x2, y2, z2, method)

            do j = 1, nj_dense  ! Make sure of the symmetry plane
               if (abs (y2(j)) < tolerance) y2(j) = zero
            end do

            dense_spokes(1)%x(i,:,1) = x2(:)
            dense_spokes(1)%y(i,:,1) = y2(:)
            dense_spokes(1)%z(i,:,1) = z2(:)
         end do

         deallocate (x1, y1, z1, x2, y2, z2)

      else  ! Leave well alone

         ni_dense = ni_regrid
         nj_dense = nj_spoke

         dense_spokes(1)%ni = ni_dense
         dense_spokes(1)%nj = nj_dense
         dense_spokes(1)%nk = 1

         call xyz_allocate (dense_spokes(1), ios)

         dense_spokes(1)%x(:,:,1) = spoke_grid(1)%x(1:ni_regrid,:,1)
         dense_spokes(1)%y(:,:,1) = spoke_grid(1)%y(1:ni_regrid,:,1)
         dense_spokes(1)%z(:,:,1) = spoke_grid(1)%z(1:ni_regrid,:,1)

      end if

!     Construct a search tree from the densified surface definition.

      nquad = (ni_dense - 1) * (nj_dense - 1)

      allocate (conn(3,nquad))  ! For patch # and (i,j) (not really needed here)

      nb = 1  ! Generic interfaces can't use a constant for the # surface blocks

      call build_adt (nb, dense_spokes, nquad, conn)

      select case (nose_case)

         case (1) ! Three edges along a semicircle

            ! Complete if needed

         case (2) ! Quarter circle with two edges on the circle

            ib1 = 1;  ib2 = 2   ! Patches of new_grid(:) to treat

         case (3) ! Two edges along a semicircle

            ! Complete if needed

      end select

!     Project each morphed nose point by finding the foot of the normal to the
!     nearest spoked cell.  Avoid the points in common with the main surface.

      npts = 0;  dmax = zero;  dmean = zero

      do ib = ib1, ib2
         do j = 1, new_grid(ib)%nj - 1
            do i = 1, new_grid(ib)%ni - 1
               target_xyz(1) = new_grid(ib)%x(i,j,1)
               target_xyz(2) = new_grid(ib)%y(i,j,1)
               target_xyz(3) = new_grid(ib)%z(i,j,1)

               call search_adt (target_xyz, iquad, p, q, dsqmin, true, &
                                nb, dense_spokes, nquad, conn, interp_xyz)

               if (abs (interp_xyz(2)) < tolerance) interp_xyz(2) = zero

               dmin  = sqrt (dsqmin)
               dmax  = max (dmin, dmax)
               dmean = dmean + dmin
               npts  = npts + 1

               new_grid(ib)%x(i,j,1) = interp_xyz(1)
               new_grid(ib)%y(i,j,1) = interp_xyz(2)
               new_grid(ib)%z(i,j,1) = interp_xyz(3)
            end do
         end do
      end do

      deallocate (dense_spokes(1)%x, dense_spokes(1)%y, dense_spokes(1)%z, conn)

      dmean = dmean / real (npts)

      write (luncrt, '(/, a, 1p, 2e14.6)') &
         'Largest and average nose projection distance:', dmax, dmean

      end subroutine project_to_nose

!     --------------------------------------------------------------------------
      subroutine save_result ()        ! And clean up
!     --------------------------------------------------------------------------

      call xyz_write (lunout, true, nblocks, new_grid, ios)

      if (ios /= 0) then
         write (luncrt, '(/, a)') 'Trouble saving regridded surface.'
      end if

      do ib = 1, nblocks
         deallocate (new_grid(ib)%x, new_grid(ib)%y, new_grid(ib)%z)
      end do

      deallocate (spoke_grid(1)%x, spoke_grid(1)%y, spoke_grid(1)%z)

      end subroutine save_result

   end program forebody_regrid
