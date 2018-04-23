!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program compress2d
!
!  This specialized program adjusts the single-block hyperbolic 2D volume grid
!  associated with a generatrix from CAPSULE_GRID by compressing its forward
!  portion to improve initial flow calculations.  It prompts for a handful of
!  controls.  A heuristic new forward boundary is established as a 2-space local
!  parametric cubic spline (x and y vs. arc length) and the input radial grid
!  lines are intersected with it then redistributed.  The untouched radial lines
!  are also redistributed to ensure smooth blending of the radial spacings.
!
!  The number of points in the radial direction for the output grid can be
!  changed along with the wall spacing (some constant) and outermost spacings.
!
!  Input Grid Format (from Gridgen glf file, with z = 0 everywhere):
!
!           1
!         285      101        1
!     0.00000000000e+00  8.96179650000e-06  3.59308720000e-05  8.08756720000e-05
!     1.43791190000e-04  2.24670450000e-04  3.23504450000e-04  4.40282210000e-04
!      :                  :                  :                  :
!
!     The input grid is assumed to have i along the surface and j in the radial
!     direction.
!
!  Option to Derive a 3D Volume Grid From the Compressed 2D Grid:
!
!     "Umbrella" (faceted and open at the back) full-body surface grids have
!     been found difficult to grow hyperbolic volume grids from in Gridgen.
!     Therefore, an option is provided to read a spoked form of the faceted
!     surface grid and morph the plane 2D hyperbolic grid to each azimuthal
!     station of the spoked surface.
!
!     The expected surfaces are spoked_surface_fore.g and spoked_surface_aft.g
!     from CAPSULE_GRID.  The option is invoked automatically if the first of
!     these is found in the working directory, but a prompt allows skipping this
!     volume gridding if the 2D compression controls still need to be checked
!     graphically, if the surface is axisymmetric, or if there is no aft body.
!
!  History:
!
!     01/23/12  D.A.Saunders  Initial implementation, for fore+aft body cases.
!     01/24/12    "     "     Split the output grid into 2 blocks so that one
!                             can have free stream inflow BC at jmax and the
!                             other can have supersonic outflow BC.  Splitting
!                             is at an odd index at or aft of ymax by default,
!                             but allow an override.
!     03/16/12    "     "     Compressing forebody or forebody-like sting cases
!                             where the blending aft of the body is exactly at
!                             at the aft end of the outer boundary needs control
!                             over the slope of the revised boundary there, and
!                             it doesn't need the splitting into two blocks.
!     08/28/12    "     "     For non-axisymmetric "umbrella"-type (faceted)
!                             geometries, use WARPQ3D morphing to produce a
!                             (compressed) 3D volume grid as well as the 2D one
!                             corresponding to the underlying generatrix.
!     08/29/12    "     "     Dinesh noted that the forebody diameter is a
!                             more natural reference length than body length
!                             (~4 diameters from the nose for the downstream
!                             boundary that the hyperbolic gridding should reach
!                             and multiple of the diameter for the desired up-
!                             stream boundary, depending on Mach and Alpha).
!     09/04/12    "     "     I forgot about WARPQ3D2, which has fewer arguments
!                             because work-space is allocated locally.
!     09/05/12    "     "     A prompt now allows suppression of the output 3D
!                             volume grid, either because the surface is axi-
!                             symmetric, or no aft body is present, or the 2D
!                             compression controls still need to be iterated on.
!     09/07/12    "     "     A SIAD (inflatable) case with a sharp vertex in
!                             the forebody generatrix prompted an option to
!                             redistribute the radial lines with the same
!                             relative spacing as in the uncompressed input
!                             2D volume grid.  Lines that have been compressed
!                             conflict with the assumptions of ADJUSTN, but an
!                             initial shuffle works around this as best we can.
!     09/13/12    "     "     The same relative spacing isn't much help for
!                             sharp vertices unless compression is limited.
!                             Preserving the existing wall spacing proves even
!                             worse, but it's now an option.  Rounding sharp
!                             corners is preferable if they're not critical.
!     11/13/14    "     "     The test for forebody or not needed a tolerance.
!     04/19/16    "     "     The test for forebody or not was using x when it
!                             should use y!
!
!  Author:  David Saunders, ERC, Inc. at NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_block_structure  ! Employed by the xyzq_io package
   use xyzq_io_module        ! PLOT3D-type I/O package
   use trigd

   implicit none

!  Local constants:

   integer, parameter :: &
      lungrid = 1,       &   ! Input and output 2D grid
      lunsurf = 2,       &   ! Optional input spoked surface grid(s)
      lunvol  = 3,       &   ! Optional output morphed 3D volume grid
      lunkbd  = 5,       &   ! Keyboard inputs
      luncrt  = 6            ! Screen prompts

   real, parameter :: &
      eps = 1.e-9            ! For intersection tolerance (x body length)

   real, parameter :: &
      one   = 1.0,    &
      zero  = 0.0

   logical, parameter :: &
      false = .false.,   &
      true  = .true.

   character, parameter :: &
      lcs_method * 1 = 'B'   ! "Loose" local cubic spline fits

!  Local variables:

   integer :: &
      i, ier, ios, j, j1, j2, l, nblocks, nf, ni, ni1, ni2, nib, nipad, &
      nj, njnew, numj

   logical :: &
      calct1, calct2, cell_centered, formatted_in, formatted_out, &
      treat_as_forebody

   real :: &
      diameter, ds, dt, dt1, dt2, dt2fr, f1, f2, shock_angle, &
      te, tint1, tint2, tol, total, xf1, xf2, xint, xmin, yint, ymax, ymin

   real, allocatable, dimension (:) :: &
      derivs, tb, xb, yb, tbpad, xbpad, ybpad, xboundary, yboundary, &
      told, xold, yold, tnew, xnew, ynew, trel, xrel, yrel, z1, z2

   type (grid_type), pointer :: &
      grid_in(:), spoked_surface_fore(:), spoked_surface_aft(:)

   type (grid_type) :: &
      grid_out

   type (grid_type) :: &
      morphed_plane, morphed_volume, unmorphed_plane

!  Execution:

   call file_prompt (lungrid, 0, 'input 2D hyperbolic grid', 'old', true, &
                     formatted_in, ios)
   if (ios /= 0) go to 99

   call xyzq_read (lungrid, -lungrid, formatted_in, nblocks, nf, &
                   cell_centered, grid_in, ios)
   if (ios /= 0) go to 99

   ni = grid_in(1)%ni
   nj = grid_in(1)%nj

   write (luncrt, '(a, i4, a, i4)') ' Dimensions found:', ni, ' x', nj

   formatted_out = true
   call file_prompt (lungrid, 0, 'output grid', 'unknown', false, &
                     formatted_out, ios)

!  Transfer the input outer boundary, some of which will remain untouched:

   nib = 6  ! For now
   allocate (xboundary(ni), yboundary(ni), xb(nib), yb(nib), tb(nib))

   xboundary(:) = grid_in(1)%x(:,nj,1)
   yboundary(:) = grid_in(1)%y(:,nj,1)

   xmin     = minval (grid_in(1)%x(:,1,1))
   ymin     = minval (grid_in(1)%y(:,1,1))
   ymax     = maxval (grid_in(1)%y(:,1,1))
   diameter = (ymax - ymin)*2.
   tol      = diameter * eps
   treat_as_forebody = abs (grid_in(1)%y(ni,1,1) - grid_in(1)%y(1,1,1)) > tol

   write (luncrt, '(a)', advance='no') &
      ' # diameters forward of the nose defining new upstream boundary: '
   read  (lunkbd, *) f1

   xf1 = xmin - f1 * diameter

   if (treat_as_forebody) then
      f2 = zero
      write (luncrt, '(a)', advance='no') &
         ' Apparently a forebody? Estimate the limiting shock angle: '
      read (lunkbd, *) shock_angle
   else
      write (luncrt, '(a)', advance='no') &
         ' # diameters aft of the nose defining boundary blend point:      '
      read  (lunkbd, *) f2
   end if

   xf2 = xmin + f2 * diameter

   write (luncrt, '(a)', advance='no') ' Desired # radial points in output:   '
   read  (lunkbd, *) njnew

   write (luncrt, '(a)', advance='no') &
   ' Redistribution [-1 = same relative | 0 = existing ds1 | const. ds1 > 0.]: '
   read  (lunkbd, *) ds

   if (ds >= zero) then
      dt1 = ds
      write (luncrt, '(a)', advance='no') &
         ' Fraction of 1-sided ds2 for 2-sided: '
      read  (lunkbd, *) dt2fr
   end if

!  Locate the index nearest the blend point:

   if (treat_as_forebody) then
      l = ni - 2
   else
      l = ni / 2
      call interval (ni, xboundary, xf2, one, l)
   end if

!  Set up points defining the compressed boundary:

   xb(1) = xf1 + 1.e-5 * diameter;  yb(1) = -yboundary(3)
   xb(2) = xf1;                     yb(2) = zero
   xb(3) = xb(1);                   yb(3) = yboundary(3)
   xb(4) = xboundary(l);            yb(4) = yboundary(l)
   xb(5) = xboundary(l+1);          yb(5) = yboundary(l+1)
   xb(6) = xboundary(l+2);          yb(6) = yboundary(l+2)

   if (treat_as_forebody) then  ! Control the limiting boundary slope
      shock_angle = tand (shock_angle)
      do i = 4, 5
         yb(i) = yb(6) - (xb(6) - xb(i)) * shock_angle
      end do
      l = ni  ! Used below for suppressing redundant intersection calculations
   end if

!  Pad the intended boundary to help the intersection starting guesses:

   call chords2d (nib, xb, yb, false, total, tb)

!! write (50, '(i3, 1p, 3e16.6)') (i, xb(i), yb(i), tb(i), i = 1, nib)

   nipad = ni / 3
   dt = total / real (nipad - 1)
   allocate (tbpad(nipad), xbpad(nipad), ybpad(nipad), derivs(nipad))

!! open (9, file='compressed_boundary.dat', status='unknown')
   do i = 1, nipad
      tbpad(i) = dt * real (i - 1)
   end do
   call lcsfit (nib, tb, xb, true, lcs_method, nipad, tbpad, xbpad, derivs)
   call lcsfit (nib, tb, yb, true, lcs_method, nipad, tbpad, ybpad, derivs)
   deallocate (derivs)
!! write (9, '(1p, 2e16.8)') (xbpad(i), ybpad(i), i = 1, nipad)
!! close (9)

!  For each radial line affected by the revised boundary, find the intersection
!  and redistribute the points.  Radial lines unaffected are also redistributed.

   numj = max (nj, njnew)
   allocate (xold(nj), yold(nj), told(nj), &
             xnew(njnew), ynew(njnew), tnew(njnew), derivs(numj))
   allocate (grid_out%x(ni,njnew,1), grid_out%y(ni,njnew,1), &
             grid_out%z(ni,njnew,1))
   allocate (trel(nj), xrel(nj), yrel(nj), z1(nj), z2(njnew))

   z1(:) = zero

   ymax = zero

   do i = 1, ni

      xold(:) = grid_in(1)%x(i,:,1)
      yold(:) = grid_in(1)%y(i,:,1)

      call chords2d (nj, xold, yold, false, total, told)

      if (i == 1) then
         tint1 = xmin - xf1  ! Arc length along 1st radial line at new boundary
         j1 = nj / 2
         call interval (nj, xold, tint1, -one, j1)  ! j1 estimates radial index
         j2 = 2  ! New boundary intersection index estimate
      end if

      if (i < l) then  ! Perform the intersection
         call intsec2t (nj,    xold,  yold,  told,  j1, false, &
                        nipad, xbpad, ybpad, tbpad, j2, false, tol, &
                        xint, yint, tint1, tint2, -luncrt, ier)
         if (ier /= 0) then
            write (luncrt, '(a, 2i5, a)') &
               ' Intersection trouble; i, ier:', i, ier, '; proceeding.'
         end if

         if (i == 1) yint = zero  ! Precisely
      else
         tint1 = total
      end if

      if (ds < zero) then  ! Same relative may help cases with sharp corners

         if (i < l) then  ! Awkward situation: the line has been shortened

!           Impose the original relative spacing on the shortened arc length:

            trel(:) = (tint1/total) * told(:)
            numj = j1 + 1
            xold(numj) = xint
            yold(numj) = yint
            told(numj) = tint1

            call lcsfit (numj, told, xold, true, lcs_method, nj, trel, xrel, &
                         derivs)
            xold(:) = xrel(:)
            call lcsfit (numj, told, yold, true, lcs_method, nj, trel, yrel, &
                         derivs)
            yold(:) = yrel(:)
         end if

         call adjustn (nj, 1, nj, xold, yold, z1, 1, njnew, xnew, ynew, z2,  &
                       z2, lcs_method)  ! ADUSTN is a 3-space utility

      else  ! Existing, or new/constant, wall spacing with revised stretching:

         if (ds == zero) dt1 = told(2) - told(1)  ! Preserve existing wall sp.

!        Redistribute, first with 1-sided stretching:

         call expdis5 (1, zero, tint1, dt1, njnew, tnew, -luncrt)

         dt2 = (tint1 - tnew(njnew-1)) * dt2fr

!        Pure 2-sided Vinokur stretching:

         call blgrid (njnew, dt1, dt2, 2, 1.001, tnew, luncrt, ier)

!        Evaluate x, y as functions of redistributed arc lengths:

         call lcsfit (nj, told, xold, true, lcs_method, njnew, tnew, xnew, &
                      derivs)
         call lcsfit (nj, told, yold, true, lcs_method, njnew, tnew, ynew, &
                      derivs)
      end if

      grid_out%x(i,:,1) = xnew(:)
      grid_out%y(i,:,1) = ynew(:)
      grid_out%z(i,:,1) = zero

      if (ynew(njnew) > ymax) then
         ymax = ynew(njnew)
         ni1  = i
      end if

   end do  ! Next radial line

   deallocate (told, xold, yold, derivs, tnew, xnew, ynew, trel, xrel, yrel, &
               z1, z2)

!  Save the result.  For boundary condition purposes, split a full-body into 2:

   if (treat_as_forebody) then
      write (lungrid, '(i2)') 1
      write (lungrid, '(2i5, i2)') ni, njnew, 1
      write (lungrid, '(1p, 6e19.11)') grid_out%x, grid_out%y, grid_out%z
   else
      if (mod (ni1, 2) == 0) ni1 = ni1 + 1  ! Some guess near ymax
      write (luncrt, '(a, i4, a)', advance='no') &
         ' Suggested (ymax) i index for splitting as 2 blocks:', ni1, &
         '; preference: '
      read  (lunkbd, *) ni1
      ni2 = ni - ni1 + 1
      write (lungrid, '(i2)') 2
      write (lungrid, '(2i5, i2)') ni1, njnew, 1, ni2, njnew, 1
      write (lungrid, '(1p, 6e19.11)') &
         grid_out%x(1:ni1,:,1),  grid_out%y(1:ni1,:,1),  grid_out%z(1:ni1,:,1)
      write (lungrid, '(1p, 6e19.11)') &
         grid_out%x(ni1:ni,:,1), grid_out%y(ni1:ni,:,1), grid_out%z(ni1:ni,:,1)
   end if

!  Option to construct a 3D volume grid for a faceted surface if it exists:

   open (lunsurf, file='spoked_surface_fore.g', status='old', iostat=ios)

   if (ios == 0) call volume3D ()

   deallocate (grid_out%x, grid_out%y, grid_out%z)

   close (lungrid)

99 continue

!  Internal procedures for COMPRESS2D:

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine volume3D ()
!
!     Construct a 3D volume grid for a faceted (non-axisymmetric) surface by
!     morphing the 2D volume grid for each azimuthal station of the faceted
!     surface grid.  The option is provided to proceed even if the (full) body
!     is axisymmetric.  Handling only a forebody is not provided for (yet).
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local variables:

      integer :: j, jr, niaft, nifore, npts, numj, numk
      real    :: angle, delta_angle
      real, allocatable :: s0(:,:,:,:), vy(:,:), vz(:,:)
      character :: answer*1

!     Execution:

      call xyzq_read (lunsurf, -lunsurf, .true., nblocks, nf, &
                      cell_centered, spoked_surface_fore, ios)
      if (ios /= 0) go to 99

      open (lunsurf, file='spoked_surface_aft.g', status='old', iostat=ios)
      if (ios /= 0) go to 99

      call xyzq_read (lunsurf, -lunsurf, .true., nblocks, nf, &
                      cell_centered, spoked_surface_aft, ios)
      if (ios /= 0) go to 99

      write (luncrt, '(a)') &
         ' Constructing a spoked 3D volume grid by morphing the 2D volume is', &
         ' recommended only for non-axisymmetric (faceted?) surfaces.'
      write (luncrt, '(a)', advance='no') ' Proceed? [y|n] '
      read  (lunkbd, '(a)') answer
      if (answer == 'n' .or. answer == 'N') go to 99

      formatted_out = true
      call file_prompt (lunvol, 0, 'output 3D grid', 'unknown', false, &
                        formatted_out, ios)

!     For 2D, y is "up" but for 3D we want y spanwise and z up, so set up a
!     copy of the compressed 2D volume grid in that form (indices likewise):

      numk = njnew
      allocate (unmorphed_plane%x(ni,1,numk), &
                unmorphed_plane%y(ni,1,numk), &
                unmorphed_plane%z(ni,1,numk))

      unmorphed_plane%x(:,1,:) = grid_out%x(:,:,1)
      unmorphed_plane%y(:,1,:) = grid_out%z(1,1,1)
      unmorphed_plane%z(:,1,:) = grid_out%y(:,:,1)

!     Storage for a morphed plane:

      allocate (morphed_plane%x(ni,1,numk), &
                morphed_plane%y(ni,1,numk), &
                morphed_plane%z(ni,1,numk))

!     Work-space for the morphing:

      allocate (s0(ni,1,numk,3), vy(ni,numk), vz(ni,numk))

!     Normalized arc lengths of 2D volume grid being morphed:

      call paramxyz (1, ni, 1, 1, 1, numk, 1, ni, 1, 1, 1, numk, &
                     unmorphed_plane%x, unmorphed_plane%y, unmorphed_plane%z, &
                     s0)

!     Storage for the output volume grid:

      numj = spoked_surface_fore(1)%nj
      allocate (morphed_volume%x(ni,numj,numk), &
                morphed_volume%y(ni,numj,numk), &
                morphed_volume%z(ni,numj,numk))

!     Morph the compressed 2D grid to every station of the desired 3D surface.
!     The spoked surface starts with j = 1 at 6 o'clock, but we output the j
!     planes in reverse order to match REVOLVE_GRID, which starts at 12 o'clock.

      nifore      = spoked_surface_fore(1)%ni
      niaft       = spoked_surface_aft(1)%ni
      npts        = ni*numk
      delta_angle = 180./ real (numj - 1)

      do j = 1, numj  ! Spoked surface starts at 6 o'clock; REVOLVEd at 12
         angle = delta_angle * real (j-1)
         if (j == numj) angle = 180.          ! Exactly
         angle = 180. - angle

         morphed_plane%x = unmorphed_plane%x  ! Really need just the edges
         morphed_plane%y = unmorphed_plane%y
         morphed_plane%z = unmorphed_plane%z

!        Set up the desired new surface edge at this j:

         morphed_plane%x(1:nifore,1,1) = spoked_surface_fore(1)%x(:,j,1)
         morphed_plane%x(nifore:ni,1,1)= spoked_surface_aft(1)%x(niaft:1:-1,j,1)
         morphed_plane%y(1:nifore,1,1) = spoked_surface_fore(1)%y(:,j,1)
         morphed_plane%y(nifore:ni,1,1)= spoked_surface_aft(1)%y(niaft:1:-1,j,1)
         morphed_plane%z(1:nifore,1,1) = spoked_surface_fore(1)%z(:,j,1)
         morphed_plane%z(nifore:ni,1,1)= spoked_surface_aft(1)%z(niaft:1:-1,j,1)

!        Revolve the new surface edge to match the base 2D at 12 o'clock:

         if (j < numj) call rotate2d (ni, morphed_plane%y(:,1,1), &
                                      morphed_plane%z(:,1,1), angle, zero, zero)

!        Morph the base 2D volume grid for this (rotated) new surface edge:

         call warpq3d2 (1, ni, 1, 1, 1, numk, 1, ni, 1, 1, 1, numk, &
                       unmorphed_plane%x, unmorphed_plane%y, unmorphed_plane%z,&
                       s0, morphed_plane%x, morphed_plane%y, morphed_plane%z)

!        Revolve the full morphed 2D grid back to this surface j station:

         vy(:,:) = morphed_plane%y(:,1,:)  ! Not necessarily exactly zero
         vz(:,:) = morphed_plane%z(:,1,:)

         if (j < numj) call rotate2d (npts, vy, vz, -angle, zero, zero)

         if (j == numj) vy(:,:) = zero  ! Exactly

         jr = numj - j + 1
         morphed_volume%x(:,jr,:) = morphed_plane%x(:,1,:)
         morphed_volume%y(:,jr,:) = vy(:,:)
         morphed_volume%z(:,jr,:) = vz(:,:)
      end do

      write (lunvol, '(i2)') 1
      write (lunvol, '(3i5)') ni, numj, numk
      write (lunvol, '(1p, 6e19.11)') &
         morphed_volume%x, morphed_volume%y, morphed_volume%z
      close (lunvol)

      deallocate (unmorphed_plane%x, unmorphed_plane%y, unmorphed_plane%z)
      deallocate (  morphed_plane%x,   morphed_plane%y,   morphed_plane%z)
      deallocate ( morphed_volume%x,  morphed_volume%y,  morphed_volume%z)
      deallocate (s0, vy, vz)

   99 return

      end subroutine volume3D

   end program compress2d
