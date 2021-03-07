!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program radial_interp_2d
!
!  Description:
!
!  The original RADIAL_INTERP reads a single-layer 3-space volume grid and a
!  related surface grid to produce a new volume grid with that precise surface.
!  It takes advantage of the situation common to CFD for reentry vehicles where
!  all blocks have k = 1:nk from the wall to the outer boundary.  For each new
!  surface point, it performs a surface interpolation giving fractional coord-
!  ates (p,q) within some cell (i,j) of some surface patch n.  Then it applies
!  the same surface interpolation at every k plane for k = 2:nk at that (i,j)
!  of block n.  Along the way, any discrepancy between the new target surface
!  point and the best interpolated surface point is either decayed towards zero
!  at the outer boundary or retained at every k plane.
!
!  All of these interpolations are necessarily bilinear because of the rapid
!  structured surface search technique used.  This is not as bad as it sounds,
!  because RADIAL_INTERP always substitutes the desired new target surface point
!  for each interpolated surface point, and decays any mismatch (or not, as
!  already explained).
!
!  Note that the number of k planes remains unchanged in the output 3D grid.
!
!  2D Analogue:
!
!  RADIAL_INTERP_2D performs a similar, simpler function for the 2D case using
!  linear interpolations with the same sort of correction if needed at the
!  surface and decay options away from the surface.  It also includes part of
!  the functionality of REFINE_GRID by allowing the type of radial distribution
!  and/or the number of points in the off-wall direction (j) to be changed.
!
!  Expected DPLR2D-related Usage:
!
!  If a coarsened 2D grid needs reclustering in the off-wall direction during
!  the flow calculation, the result no longer matches every nth point of the
!  original fine grid.  Therefore, the requirement that prompted this utility is
!  to recluster that fine grid in the equivalent manner, which has typically in-
!  volved cell-Reynolds-number-based spacing at the wall.  Preserving the same
!  relative spacing in the radial direction while (normally) doubling the number
!  of cells is the right choice, along with recovering the original surface grid
!  exactly, which simpler use of REFINE_GRID would not guarantee.
!
!  For convenience, the desired surface grid may be the j = 1 line of a 2D
!  volume grid (1 or more blocks), or it may be a simple generatrix (one curve
!  in (x,y)-space, with any third coordinate ignored).
!
!  XYZ Conventions:
!
!  Since the DPLR postprocessor extracts only x and y for a 2D grid, the input
!  volume grid to be interpolated may be either 2D/xy or 3D/xyz with z all zero.
!  The DPLR preprocessor ignores any z coordinates if they are present for 2D
!  input grids, so the output grid here includes z = 0 for compatibility with
!  other utilities.  Thus y is "up" for input and output files here.
!
!  Control:
!
!  A handful of prompts suffice.
!
!  History:
!
!  Feb. 2-5, 2012  D.A.Saunders  Initial design and implementation.
!  Nov. 6,   2015    "     "     Introduced determine_grid_form and -dim.
!                                Writing 3D output is inconsistent with
!                                FLOW_INTERP_2D, so make it an option.
!  Apr. 19,  2016    "     "     The option for entering a simple curve as
!                                the new surface should not invoke automatic
!                                detection of format.
!
!  Author:  David Saunders, ERC, Inc. at NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_block_structure  ! Grid block derived data type
   use xyq_io_module         ! PLOT2D file I/O package
   use xyzq_io_module        ! PLOT3D file I/O package

   implicit none

!  Constants:

   integer, parameter :: &
      lunin1 = 1,        &   ! Input 2D volume grid
      lunin2 = 2,        &   ! Input surface line or 2D grid (j=1 used)
      lunout = 3,        &   ! Output grid
      lunkbd = 5,        &   ! Keyboard inputs
      luncrt = 6             ! Prompts/diagnostics

   real, parameter :: &
      zero = 0.

!  Variables:

   integer :: &
      i, ios, j, method_decay, method_radial, n, nblocks1, nblocks2, ndim, &
      njout, npts

   real :: &
      ds1, ds2, ds2fraction

   logical :: &
      formatted_out, twod_out

   type (grid_type), pointer, dimension (:) :: &
      oldvol, newvol, target

!  Execution:

   call control ()        ! Prompt for file names, etc.;
   if (ios /= 0) go to 99 ! read the input files and set up the output file

   call interp_volume ()  ! Interpolate new radial lines at the old j levels;
                          ! optionally redistribute in one of various ways

   call save_results ()   ! Write the new volume grid and wrap up

99 continue

!  Internal procedures for RADIAL_INTERP_2D:

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine control ()  ! Prompt for input/output files and other controls;
!                            ! read input files and set up output file.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local constants:

      character, parameter :: &
         form*11 = 'unformatted'

!     Local variables:

      logical :: &
         formatted_in1, formatted_in2, surface_curve, twod

      character :: &
         answer*1, filename*80

!     Execution:

!     File name inputs, starting with the input 2D volume grid:

      write (luncrt, '(a)', advance='no') ' 2D volume grid to be interpolated: '
      read  (lunkbd, *) filename
      call determine_grid_form (filename, lunin1, formatted_in1, ios)

      if (formatted_in1) then
         open (lunin1, file=filename, status='old', iostat=ios)
      else
         open (lunin1, file=filename, status='old', form='unformatted', &
               iostat=ios)
      end if
      if (ios /= 0) go to 90

      write (luncrt, '(a)', advance='no') &
         ' 2D grid or curve containing desired new surface: '
      read  (lunkbd, *) filename
      write (luncrt, '(a)', advance='no') ' Grid or curve? [g|c]: '
      read  (lunkbd, '(a)') answer
      surface_curve = answer == 'c'

      if (surface_curve) then
         formatted_in2 = .true.
      else
         call determine_grid_form (filename, lunin2, formatted_in2, ios)
      end if

      if (formatted_in2) then
         open (lunin2, file=filename, status='old', iostat=ios)
      else
         open (lunin2, file=filename, status='old', form='unformatted', &
               iostat=ios)
      end if
      if (ios /= 0) go to 99

      write (luncrt, '(a)', advance='no') 'Output grid file name: '
      read  (lunkbd, *) filename
      write (luncrt, '(a)', advance='no') 'Formatted (y) or unformatted (n):  '
      read  (lunkbd, *) answer
      formatted_out = answer == 'y' .or. answer == 'Y'

      if (formatted_out) then
         open (lunout, file=filename, status='unknown', iostat=ios)
      else
         open (lunout, file=filename, status='unknown', form='unformatted', &
               iostat=ios)
      end if
      if (ios /= 0) go to 90

      write (luncrt, '(a)', advance='no') '2-space (2) or 3-space (z=0) (3):  '
      read  (lunkbd, *) ios
      twod_out = ios == 2

!     Read the input files to help with other controls; start with old volume:

      call determine_grid_dim ('none', lunin1, formatted_in1, ndim, ios)
      twod = ndim == 2

      if (formatted_in1) then
         read (lunin1, *, iostat=ios) nblocks1
      else
         read (lunin1,    iostat=ios) nblocks1
      end if
      if (ios /= 0) then
         write (luncrt, '(2a)') 'Trouble reading # input volume grid blocks.'
         go to 99
      end if

      allocate (oldvol(nblocks1))

      if (formatted_in1) then
         if (twod) then
            read (lunin1, *, iostat=ios) &
               (oldvol(n)%ni, oldvol(n)%nj, n = 1, nblocks1)
         else
            read (lunin1, *, iostat=ios) &
               (oldvol(n)%ni, oldvol(n)%nj, oldvol(n)%nk, n = 1, nblocks1)
            do n = 1, nblocks1
               if (oldvol(n)%nk /= 1) ios = 1  ! Shouldn't be necessary
            end do
         end if
      else
         if (twod) then
            read (lunin1, iostat=ios) &
               (oldvol(n)%ni, oldvol(n)%nj, n = 1, nblocks1)
         else
            read (lunin1, iostat=ios) &
               (oldvol(n)%ni, oldvol(n)%nj, oldvol(n)%nk, n = 1, nblocks1)
            do n = 1, nblocks1
               if (oldvol(n)%nk /= 1) ios = 1  ! Shouldn't be necessary
            end do
         end if
      end if
      if (ios /= 0) then
         write (luncrt, '(2a)') &
            'Trouble reading input 2D grid block dimensions.'
         go to 99
      end if

      oldvol(:)%nk = 1
      write (luncrt, '(a, /, (3i5))') &
         'Input 2D volume grid block dimensions found:', &
         (oldvol(n)%ni, oldvol(n)%nj, 1, n = 1, nblocks1)

!     Read the input 2D volume grid coordinates:

      do n = 1, nblocks1
         call xyz_allocate (oldvol(n), ios)

         if (twod) then
            if (formatted_in1) then
               read (lunin1, *, iostat=ios) oldvol(n)%x, oldvol(n)%y
            else
               read (lunin1,    iostat=ios) oldvol(n)%x, oldvol(n)%y
            end if
            if (ios /= 0) then
               write (luncrt, '(a, i4)') &
                  'Trouble reading (x,y) from 2D input volume; block #:', n
               go to 99
            end if
            oldvol(n)%z(:,:,:) = zero
         else
            npts = oldvol(n)%ni * oldvol(n)%nj
            call xyz_block_io (1, lunin1, formatted_in1, npts, &
                               oldvol(n)%x, oldvol(n)%y, oldvol(n)%z, ios)
            if (ios /= 0) then
               write (luncrt, '(a, i4)') &
                  'Trouble reading (x,y,z) from 3D input volume; block #:', n
               go to 99
            end if
         end if
      end do

      close (lunin1)

!     Likewise for the curve or grid containing the desired new surface grid:

      if (surface_curve) then

         nblocks2 = 1
         allocate (target(1))
         npts = 0
         do ! Until EOF
            read (lunin2, *, iostat=ios)  ! Count the number of points
            if (ios < 0) exit
            npts = npts + 1
         end do
         rewind (lunin2)
         write (luncrt, '(a, i5)') '# points found in new surface curve:', npts
         target(1)%ni = npts;  target(1)%nj = 1;  target(1)%nk = 1
         call xyz_allocate (target(1), ios)

         do i = 1, npts  ! Avoid any third column
            read (lunin2, *, iostat=ios) target(1)%x(i,1,1), target(1)%y(i,1,1)
            if (ios /= 0) then
               write (luncrt, '(a, i6)') &
                  'Trouble reading target surface curve; line #:', i
               go to 99
            end if
            target(1)%z(i,1,1) = zero
         end do

      else  ! Read 2- or 3-space grid containing the target surface

         call determine_grid_dim ('none', lunin2, formatted_in2, ndim, ios)
         twod = ndim == 2

         if (formatted_in2) then
            read (lunin2, *, iostat=ios) nblocks2
         else
            read (lunin2,    iostat=ios) nblocks2
         end if
         if (ios /= 0) then
            write (luncrt, '(2a)') 'Trouble reading # input target grid blocks.'
            go to 99
         end if

         allocate (target(nblocks2))

         if (twod) then
            if (formatted_in2) then
               read (lunin2, *, iostat=ios) &
                  (target(n)%ni, target(n)%nj, n = 1, nblocks2)
            else
               read (lunin2,    iostat=ios) &
                  (target(n)%ni, target(n)%nj, n = 1, nblocks2)
            end if
         else
            if (formatted_in2) then
               read (lunin2, *, iostat=ios) &
                  (target(n)%ni, target(n)%nj, target(n)%nk, n = 1, nblocks2)
            else
               read (lunin2,    iostat=ios) &
                  (target(n)%ni, target(n)%nj, target(n)%nk, n = 1, nblocks2)
            end if
         end if
         if (ios /= 0) then
            write (luncrt, '(2a)') &
               'Trouble reading target grid block dimensions.'
            go to 99
         end if

         target(:)%nk = 1
         write (luncrt, '(a, /, (3i5))') &
            'Input target grid block dimensions found:', &
            (target(n)%ni, target(n)%nj, 1, n = 1, nblocks2)

!        Read the target grid coordinates:

         do n = 1, nblocks2
            call xyz_allocate (target(n), ios)

            if (twod) then
               if (formatted_in2) then
                  read (lunin2, *, iostat=ios) target(n)%x, target(n)%y
               else
                  read (lunin2,    iostat=ios) target(n)%x, target(n)%y
               end if
               if (ios /= 0) then
                  write (luncrt, '(a, i4)') &
                     'Trouble reading (x,y) from 2D target grid; block #:', n
                  go to 99
               end if
               target(n)%z(:,:,:) = zero
            else
               npts = target(n)%ni * target(n)%nj
               call xyz_block_io (1, lunin2, formatted_in2, npts, &
                                  target(n)%x, target(n)%y, target(n)%z, ios)
               if (ios /= 0) then
                  write (luncrt, '(a, i4)') &
                     'Trouble reading (x,y,z) from 3D target grid; block #:', n
                     go to 99
               end if
            end if
         end do

      end if

      close (lunin2)

      allocate (newvol(nblocks2)) ! Output grid has as many blocks as the target

!     Other control inputs:

      write (luncrt, '(a)') 'Radial grid spacing options:', &
         '  1 => same spacing and nj as in the input volume grid.', &
         '  2 => same relative spacing as input volume; nj from target grid.', &
         '  3 => specified new wall spacing ds1, outer ds2fraction, and nj.'
      write (luncrt, '(a)', advance='no') 'Radial spacing method?  [1|2|3] '
      read  (lunkbd, *) method_radial

      if (method_radial == 1) then
         njout = oldvol(1)%nj
      else if (method_radial == 2) then
         njout = target(1)%nj
      else if (method_radial == 3) then
         write (luncrt, '(a)', advance='no') 'New # radial pts. nj: '
         read  (lunkbd, *) njout
         write (luncrt, '(a)', advance='no') 'ds1 and ds2fraction:  '
         read  (lunkbd, *) ds1, ds2fraction
      end if

      write (luncrt, '(a)') 'Options for decay of surface mismatches:', &
         '  1 => decay any mismatches to zero at the outer boundary.', &
         '  2 => carry any mismatches all the way to the boundary.'
      write (luncrt, '(a)', advance='no') 'Mismatch decay method?  [1|2] '
      read  (lunkbd, *) method_decay

      go to 99

90    write (luncrt, '(2a)') 'Unable to open ', trim (filename)

99    continue

      end subroutine control

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine interp_volume ()

!     Construct the 2D volume grid that corresponds to the input volume grid and
!     a different surface grid (assumed to be at j = 1 for all blocks):
!
!     For each (i,1) point of [each block of] the new surface grid, locate the
!     cell ic of the piecewise linear curve formed by the j = 1 edges of the
!     input volume grid (oldvol), and apply the associated interpolation coef-
!     ficients at each j in oldvol(ic,2:nj) for that i.  Correct for any mis-
!     match at j = 1 as explained above.
!
!     Optionally, adjust the type of radial distribution and/or number of radial
!     points as we go.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local constants:

      integer,   parameter :: ngeometric = 2      ! No part-geometric option 
      real,      parameter :: growth_rate = 1.01  ! Not used
      real,      parameter :: one = 1.
      character, parameter :: lcs_method*1 = 'B'  ! "Bessel" (loose) fits

!     Local variables:

      integer :: nbest, ib, ibest, ic, ier, j, jp, ni, njold
      real    :: decay, dsq, dsqmin, dx, dy, p, q, pbest, qbest, r, rp, total, &
                 xb, yb, zb, xbest, ybest
      real, allocatable, dimension (:) :: si, xi, yi  ! Interim radial line
      real, allocatable, dimension (:) :: so, xo, yo  ! Output radial line
      real, allocatable, dimension (:) :: derivs      ! Unused

!     Execution:

      njold = oldvol(1)%nj

      allocate (si(njold), xi(njold), yi(njold))
      allocate (so(njout), xo(njout), yo(njout), derivs(njout))

      do ib = 1, nblocks2  ! For each target block

         newvol(ib)%ni = target(ib)%ni
         newvol(ib)%nj = njout
         newvol(ib)%nk = 1

         call xyz_allocate (newvol(ib), ios)

         do i = 1, target(ib)%ni  ! For each target i at j = 1

            dsqmin = 1.e+20

            do n = 1, nblocks1  ! For each input volume grid block

               call nearest_curve_point (oldvol(n)%ni,        &
                                         oldvol(n)%x(:,1,1),  &
                                         oldvol(n)%y(:,1,1),  &
                                         oldvol(n)%z(:,1,1),  &
                                         target(ib)%x(i,1,1), &
                                         target(ib)%y(i,1,1), &
                                         target(ib)%z(i,1,1), &
                                         ic, xb, yb, zb, p, q, dsq)
               if (dsq < dsqmin) then
                   dsqmin = dsq
                   nbest  = n
                   ibest  = ic
                   pbest  = p
                   qbest  = q
                   xbest  = xb
                   ybest  = yb
               end if

            end do  ! Next wall portion of input volume grid

            xi(1) = target(ib)%x(i,1,1);  dx = xi(1) - xbest
            yi(1) = target(ib)%y(i,1,1);  dy = yi(1) - ybest
            si(1) = zero

!           Apply the wall interpolation coefficients at each j off the wall:

            do j = 2, njold
               xi(j) = pbest * oldvol(nbest)%x(ibest,j,1) + &
                       qbest * oldvol(nbest)%x(ibest+1,j,1)
               yi(j) = pbest * oldvol(nbest)%y(ibest,j,1) + &
                       qbest * oldvol(nbest)%y(ibest+1,j,1)
               si(j) = sqrt ((xi(j) - xi(j-1))**2 + (yi(j) - yi(j-1))**2) + &
                       si(j-1)
            end do

!           Options to correct for any discrepancy at the surface:

            if (method_decay == 1) then  ! Decay any discrepancy to zero at nj

               total = one / si(njold)
               do j = 2, njold
                  decay = (si(njold) - si(j)) * total
                  xi(j) = xi(j) + decay * dx
                  yi(j) = yi(j) + decay * dy
               end do

            else  ! Apply the same surface disrepancy at each j

               do j = 2, njold
                  xi(j) = xi(j) + dx
                  yi(j) = yi(j) + dy
               end do

            end if

!           Update the radial arc lengths in case of decay:

            call chords2d (njold, xi, yi, .false., total, si)

!           Options to change the radial distribution and/or # points:

            if (method_radial == 1) then  ! Keep the old radial spacing

               xo(:) = xi(:)
               yo(:) = yi(:)

            else if (method_radial == 2) then  ! Target volume nj, same relative

               r = real (njold - 1) / real (njout - 1)
               so(1) = zero
               do j = 2, njout - 1
                  rp = one + r * real (j - 1)
                  jp = int (rp)
                  p  = rp - real (jp)
                  so(j) = (one - p) * si(jp) + p * si(jp+1)
               end do
               so(njout) = total

            else if (method_radial == 3) then  ! Specified ds1 and ds2fraction

               call expdis5 (1, zero, total, ds1, njout, so, -luncrt)  ! 1-sided

               ds2 = (total - so(njout-1)) * ds2fraction

               call blgrid (njout, ds1, ds2, ngeometric, growth_rate, so, &
                            luncrt, ier)                               ! 2-sided
            end if

            if (method_radial /= 1) then  ! Arc-length-based redistribution

               call lcsfit (njold, si, xi, .true., lcs_method, &
                            njout, so, xo, derivs)
               call lcsfit (njold, si, yi, .true., lcs_method, &
                            njout, so, yo, derivs)
            end if

            newvol(ib)%x(i,:,1) = xo(:)
            newvol(ib)%y(i,:,1) = yo(:)
            newvol(ib)%z(i,:,1) = zero

         end do  ! Next i for this target block

      end do  ! Next target block

      end subroutine interp_volume

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine save_results ()
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (twod_out) then
         call xy_write  (lunout, formatted_out, nblocks2, newvol, ios)
      else
         call xyz_write (lunout, formatted_out, nblocks2, newvol, ios)
      end if

      end subroutine save_results

   end program radial_interp_2d
