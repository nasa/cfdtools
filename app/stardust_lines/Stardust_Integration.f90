!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program Stardust_Integration
!
!  Description:
!
!     Stardust_Integration is a variant of the neqair_integration utility to
!  accompany the Stardust_Lines line-of-sight generation scheme for the step
!  that integrates the radiation solver calculations.  Instead of a triangulated
!  hemisphere defining lines of sight to a spacecraft body point, it works with
!  a rectangular array of pixels defining parallel lines orthogonal to the pixel
!  plane.  The same main options are provided: integration of a single radiance
!  at each pixel (over solid angle subtended at a remote viewing point, giving a
!  heat flux) or (more likely) a comparable integration for each wavelength in
!  the spectral radiance outputs from the solver.  In the case of NEQAIR, these
!  large outputs are intensity.out or intensity_scanned.out. In the simple case,
!  each radiance is read from the last line of neqair.out, or (as a check) from
!  the single-line file named radiance.out written during a preceding spectral
!  integration run of this utility.
!
!  Method:
!
!     o  Read the sensor location relative to the spacecraft position (details
!        TBD).  This relative location was also needed to determine the orienta-
!        tion of the pixel plane that defines the lines of sight processed here.
!     o  Determine from the /LINE-1 subdirectory the kind of outputs from
!        NEQAIR that are to be integrated.  This defines a column number as
!        0, 1, 2, or 3 to be read from.
!     o  Read the pixel array as pixels.g from Stardust_Lines (a PLOT3D surface
!        grid in 3-space).
!     o  Read the header of the line segments output by Stardust_Lines.  Any
!        line with nk = 2 is assigned zero for each output from NEQAIR, which
!        is not invoked for such lines that did not pass through the flow field.
!        The (i,j) corresponding to each meaningful line in /LINE-n for any n
!        can also be determined.
!
!     o  [Simple radiance case:] For each meaningful line, read a single
!        radiance for the indicated pixel, and enter zeros for the other lines.
!        Integrate over the rectangular pixel array w.r.t. solid angle subtended
!        at the remote sensor location.  We are done.
!
!     o  [Spectral case:] Read the wavelengths specified in the NEQAIR control
!        file, neqair.inp.
!     o  Scan the NEQAIR outputs from all meaningful lines and establish a
!        single set of wavelengths that are uniform in each region and as dense
!        as the densest line found in each region.
!     o  Suppress the duplicate wavelengths at the region boundaries that are
!        present in the NEQAIR outputs.
!     o  As many sums as there are uniformly distributed wave-lengths in all
!        the regions are accumulated one line at a time as in the earlier
!        neqair_integration scheme, as follows:
!        o  For each cell in the pixel array, calculate the solid angle domega
!           subtended at the sensor location.
!        o  For each quad. cell vertex, the "weight" corresponding to each
!           vertex of each triangle in the hemisphere case is domega/4 * cosine
!           of the angle between the vertex and the central pixel line to the
!           remote sensor, as needed to isolate the normal component at the
!           sensor.  (Even though these angles will be extremely small, include
!           them anyway.)
!        o  For each pixel, skip it if the line of sight didn't pass through the
!           flow field (nk = 2), else read the spectral radiances associated
!           with the line of sight and interpolate to the regularized wave-
!           lengths.
!        e  Increment the sums for all (regularized) wavelengths by the
!           vertex weight*(regularized) radiances.
!        o  Integrate the regularized radiances w.r.t. wavelength, and write
!           a one-line 'radiances.out' result in that line's subdirectory.
!        o  Save the solid angle integrations for each regularized wavelength
!           as 'integrations.dat'.  This file contains spectral irradiances.
!        o  Integrate these w.r.t. wavelength as a check on the heat-flux case.
!           This result is converted from angstroms to microns to give W/cm^2.
!
!  (Later): Dinesh and Jay believe the spacecraft is effectively a point source,
!           so we should just sum results at all pixels.  These are written as
!           spectral_sums.dat (one sum per regularized wavelength).
!           NO: For such a file to converge with increasing resolution, divide
!           by the number of active pixels, ASSUMING DENSER RESOLUTIONS RETAIN
!           THE SAME RELATIVE SPACING.  This means changing the d1 used for
!           one-sided stretching proportionately.  E.g.: go from  65  0.0005
!           to  129  0.00025 or to  97  0.00033333333333...
!           Call the file spectral_irradiance.dat and double the average to
!           account for a 360-degree flow field, not 180 deg.
!  (Later   This file was never looked at and the name spectral_irradiance.dat
!   still)  is really what integrations.dat should've been named.  Suppress it.
!
!  Input files:
!
!     neqair.inp                For the specified wavelength regions
!     pixels.g                  Rectangular array of coordinates defining all
!                               lines of sight, active & inactive
!     line_segments.g           Discretized lines of sight, ni x nj of them;
!                               inactive lines are identified by nk = 2
!     lines.inp                 For the viewing angle and distance on line 8,
!                               from the Stardust_Lines control file
!  Output files:
!
!     integrations.dat          Spectral irradiances (integrations of spectral
!                               radiance over the solid angle subtended by the
!                               pixel plane at the observation location, one
!                               for each regularized wavelength). The normal
!                               components are integrated, even though all the
!                               cosines will be extremely close to 1 normally.
!                               Each values is doubled here to account for the
!                               fact that the calculations are performed for a
!                               180-degree revolution of the underlying axi-
!                               symmetric flow field, not a full revolution.
!     spectral_irradiance.dat   Average integration over active pixels, one for
!                               each regularized wavelength, doubled to account
!                               for the missing half of the flow field.
!     center.line.cumulative.dat   For the pixel nearest the middle only:
!                                  cumulative integrals w.r.t. wavelength, the
!                                  final value being the total radiance
!     radiance.out              Total radiance at each pixel (zero for inactive
!                               pixels)
!     pixel-indices.dat         Enables finding the LOS number associated with
!                               each pixel (i,j) for diagnostic purposes
!                               (Later:  No. Write this to the standard log file
!                               in favor of less info. from the routine that
!                               reads radiances for the single-radiance case.)
!     pixels.f                  Total integrated radiance (normal component at
!                               the observation location, not including atmos-
!                               pheric absorption) for each LOS/pixel (zero for
!                               inactive lines), from the single-radiance case.
!  History:
!     02/26/2020  D.A.Saunders  Initial design and implementation.
!     03/24/2020     "    "     Many delays later, it compiles but we lack the
!                               details giving the relative position of the
!                               sensor in the CFD coordinate system.
!     04/20/2020     "    "     (After a hospitalization crisis): the simpler
!                               "heat-flux" mode works.  Introduced a distance
!                               alongside the viewing angle driving the
!                               Stardust_Lines step; assume 'lines.inp' is the
!                               control file name.  Varying the distance for a
!                               given viewing angle should allow a sanity check
!                               on the single-radiance-per-line case before we
!                               exercise the spectral case.
!     05/06/2020     "    "     Belated finding that the cosines needed for the
!                               normal components of radiance (at the sensor)
!                               were using (xm,ym,zm) instead of (px,py,pz).
!                               Write spectral_sums.dat as explained above under
!                               (Later). Cosine calcs. have since been revised.
!     05/11/2020     "    "     Double the spectral_sums.dat values.
!     05/19/2020     "    "     Fix nomenclature error (spectral irradiance
!                               after integrating spectral radiance w.r.t. solid
!                               angle).  Change spectral_sums.dat file name to
!                               (average) spectral_irradiance.dat because the
!                               sums are now divided by the number of active
!                               pixels in order for these quantities to converge
!                               as (proportional) pixel resolution increases.
!                               See details about stretching input d1 above.
!     05/21/2020     "    "     Brett found that f9.1 meant too few digits in
!                               integrations.dat for wavelengths and other similar
!                               outputs.  Go to f12.4.
!     08/30/2020     "    "     After running NEQAIR without the black body BC,
!                               added output of pixel-indices.dat to associate
!                               pixel indices with active LOS numbers.
!                               Added full input/out file descriptions.
!     09/03/2020-    "    "     Instead of defining distance to observer (slant
!     09/09/2020                range) as being from the nose (0,0,0) in the
!                               CFD coordinate system, define it as the distance
!                               to the middle of the pixel rectangle: negligibly
!                               different for large slant ranges, but it cleans
!                               up the cosine calculations for normal components
!                               of radiance at the observer location, and the
!                               calculation of that location, (px,py,pz).
!     09/26/2020     "    "     (New) negative integrations were due to a sign
!                               error in in the calculation of px.
!     10/07/2020     "    "     Belated realization that the integrations over
!                               solid angle for each wavelength in the file
!                               integrations.dat should be doubled to account
!                               for the missing half of the flow field.
!     11/23/2020     "    "     The units of the averaged results over active
!                               pixels, and the spectral_irradiance.dat choice
!                               of file name, are problematic.  This file has
!                               not been made use of, so suppress it.  Also,
!                               the cumulative integrations for the line of
!                               closest to the middle of the pixel plane (with
!                               y = 0) haven't been showing up as intended in
!                               center.line.cumulative.dat and nobody cared, so
!                               the reason remains a mystery.
!
!  Author:  David Saunders, AMA, Inc. at NASA Ames Research Center.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_block_structure    ! Employed by the xyzq_io package
   use xyzq_io_module          ! PLOT3D-type I/O package

   implicit none

!  Constants:

   integer, parameter :: &
      lunctl = 1,        &   ! Stardust_Lines control file (for angle+distance)
      lunpix = 2,        &   ! Input pixel PLOT3D surface grid (Stardust_Lines)
      lunLOS = 3,        &   ! All lines of sight, meaninful or (nk = 2) not
      lunrad = 4,        &   ! For NEQAIR control file and for output results
      lunkbd = 5,        &   ! Keyboard inputs
      luncrt = 6,        &   ! Screen prompts and diagnostics
      lunout = 7             ! For output pixel-indices.dat and pixels.f

   real, parameter :: &
      angstrom_to_micron = 1.e-4, &  ! Brett says wavelengths need to be microns
      half = 0.5, one = 1., zero = 0.

   logical, parameter :: &
      false = .false.,   &
      true  = .true.

!  Variables:
 
   integer :: &
      i, icol, imid, ios, j, jmid, line, niLOS, njLOS, nlines, npixelblocks, &
      npixels, nregion, numf, nuniform

   real :: &
      angle, biggest_cosine, cosine, D, distance, px, py, pz, qrad, xm, &
      xmax, xmin, ym, ymax, ymin, zm, zmax, zmin

   logical :: &
      formatted, unused

   character (10) :: &
      prefix         ! Path for saving wavelength-dependent results

   integer, allocatable, dimension (:) :: ip, jp       ! Nonzero pixel indices
   integer, allocatable, dimension (:) :: iuniform     ! Unif. region indices
   integer, allocatable, dimension (:,:) :: irange     ! Indices for all hemi
                                                       ! lines; only the data
                                                       ! lines are assigned
   real,    allocatable, dimension (:) :: cumulative   ! Cumulative integrals
                                                       ! w.r.t. wavelength, for
                                                       ! possible plotting,
                                                       ! center pixel line only
   real,    allocatable, dimension (:) :: lambda_sum   ! One integration sum for
                                                       ! each regularized wave-
                                                       ! length
   real,    allocatable, dimension (:) :: simple_sum   ! For simply summing all
                                                       ! pixels at each uniform
                                                       ! wavelength
   integer, allocatable, dimension (:) :: nrecords     ! # wavelengths found in
                                                       ! each input intensity
                                                       ! file (with possible
                                                       ! duplicates at interior
                                                       ! region boundaries)
   real,    allocatable, dimension (:) :: uni_lambda   ! Uniform wavelengths
                                                       ! common to all lines
   real,    allocatable, dimension (:) :: uni_radiance ! Interpolated spectral
                                                       ! radiances for one data
                                                       ! line
   real,    allocatable, dimension (:,:) :: vertex_weight  ! One weight for each
                                                           ! vertex, from one or
                                                           ! more quad. cells
   real,    pointer,     dimension (:) :: range        ! Region boundary wave-
                                                       ! lengths, dim. nregion+1
   type (grid_type), pointer, dimension (:) :: pixels  ! Rectangular plane
   type (grid_type), pointer, dimension (:) :: LOS     ! niLOS*njLOS line of
                                                       ! sight segments, some
                                                       ! not meaningful (nk = 2)
!  ----------
!  Execution:
!  ----------

!  The observing location relative to the spacecraft may be needed to determine
!  the viewing angle input to Stardust_Lines, but a relative distance at least
!  is needed here.  First approach:  read the viewing angle and an adjacent
!  distance from line 8 of the Stardust_Lines control file, hard-coded as
!  'lines.inp'.  Integrating for a range of distances should provide a sanity
!  check for the single-radiance case to start with.

   call get_sensor_details ()
   if (ios /= 0) go to 99

!  Read the pixels.g and line-segments files written by Stardust_Lines:

   call get_pixels_and_line_segments ()
   if (ios /= 0) go to 99

!  Determine from the /LINE-1 subdirectory the kind of outputs from NEQAIR that
!  are to be integrated.  This defines a column number as 1, 2, or 3 to be read,
!  or 0 for the simple radiance case [neqair.out].

   call get_neqair_case (lunrad, icol)

!  Retrieve pixel indices for meaningful lines and initialize the integrations:

   call initialize_integrations ()

!  ------------------------------------------------------------
!  Simple case (one NEQAIR result per meaningful line segment):
!  ------------------------------------------------------------

   if (icol <= 1) then  ! One radiance per neqair.out (0) or radiance.out (1)

      call single_radiance_case ()
      go to 99          ! Skip past the rest, which is for spectral cases only

   end if

!  ---------------
!  Spectral cases:
!  ---------------

!  Identify the wavelengths specified in the NEQAIR control file.
!  Handling NEQAIR v15+ (only) should suffice.

   open (lunrad, file='neqair.inp', status='old', iostat=ios)
   if (ios /= 0) then
      write (luncrt, '(a)') &
         ' Unable to open neqair.inp for the region definitions.'
      go to 99
   end if

   call get_region_data (lunrad, nregion, range, ios)
   if (ios /= 0) go to 99

!  Scan all data lines, store their region indices for rereading, and
!  set up uniform wavelengths that are as dense as the densest line found
!  in each region.

   call get_common_wavelengths ()
   if (ios /= 0) go to 99

!  Perform the integrations for all uniform wavelengths in all regions, by
!  processing one active line at a time and accumulating sums over solid
!  angle elements.  Also, integrate w.r.t. the uniform wavelengths and save a
!  'radiance.out' result for each active line for checking against NEQAIR's
!  own integrations.

   call perform_integrations ()
   if (ios /= 0) go to 99

99 continue

!  Internal procedures for program Stardust_Integration, in alphabetical order:

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine analyze_wavelengths ()
!
!     This adaptation of neqair_integration performs its original function for
!     each wavelength of each region analyzed by NEQAIR.  It is isolated to
!     leave the original single-radiance-per-linescheme essentially intact,
!     although that scheme can be reused following analysis of the
!     intensity_scanned.out or intensity.out files by processing the one-line
!     radiance.out files that are written at the end of the more elaborate
!     analysis by integrating w.r.t. wavelength for each line of sight.
!
!     To be more precise, a common set of wavelengths has to be imposed on each
!     line before any integrations.  This is done by scanning all lines for the
!     maximum number of (uniform?) wavelengths in each region and interpolating
!     the spectral radiances to uniform wavelengths that are dense enough not to
!     lose resolution in any region of any line.
!
!     Storing all of the functions (possibly millions of spectral radiances) for
!     all of the lines in order to do the integrations the way they are done in
!     the original scheme would require prohibitive amounts of memory. Instead,
!     one line can be processed at a time and a summation can be performed for
!     each (common) wavelength.  Thanks to Brett Cruden for pointing this out.
!     Consider the 1-D trapezoidal rule analog to see this more easily.
!
!     04/16/2018  David Saunders  Preliminary design.
!     03/20/2020    "      "      Adaptation for Stardust-type integrations.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Execution:

!     First, identify the wavelengths specified in the NEQAIR control file:

      open (lunrad, file='neqair.inp', status='old', iostat=ios)
      if (ios /= 0) then
         write (luncrt, '(a)') &
            ' Unable to open neqair.inp for the region definitions.'
         go to 99
      end if

      call get_region_data (lunrad, nregion, range, ios)
      if (ios /= 0) go to 99

!     Scan all data lines, store their region indices for rereading, and
!     set up uniform wavelengths that are as dense as the densest line found
!     in each region.  Store the cone angles for each hemi line as well.

      call get_common_wavelengths ()
      if (ios /= 0) go to 99

!     Perform the integrations for all uniform wavelengths in all regions, by
!     processing one active line at a time and accumulating sums over solid
!     angle elements.  Also, integrate w.r.t. the uniform wavelengths and save a
!     'radiance.out' result for each active line for checking against NEQAIR's
!     own integrations.

      call perform_integrations ()
      if (ios /= 0) go to 99

 99   return

      end subroutine analyze_wavelengths

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine get_common_wavelengths ()
!
!     For the type of output from NEQAIR indicated by icol, scan all relevant
!     line results for the wavelengths used in each region, and calculate a
!     single set of wavelengths that are uniform in each region and as dense
!     as the densest line found in each region.  The spectral radiances from
!     each line will later be interpolated to this common set of wavelengths.
!
!     The handling of radiometer cone angles in the NEQAIR_Integration version
!     has been deleted here.
!
!     04/16/2018  David Saunders  Initial implementation in NEQAIR_Integration.
!     05/18/2018    "      "      Suppress duplicate wavelengths at the region
!                                 boundaries.
!     02/26/2020    "      "      The Stardust_Integration version doesn't need
!                                 to deal with triangulated zones.  We have to
!                                 handle all the meaningful lines, and their
!                                 number has to be known prior to entry here.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local variables:  (Most variables are inherited from the main program.)

      integer :: iline, iregion, line, mode, nlam, nu
      integer, allocatable :: mxregion(:)
      logical :: dialog

!     Execution:

      allocate (mxregion(nregion))  ! For maximum region sizes found
      mxregion(:) = 0

      allocate (irange(nregion+1,nlines), & ! Lambda regions found in data lines
                iuniform(nregion+1),      & ! Regions imposed on all data lines
                nrecords(nlines))           ! For intensity data lines

      mode    = 1  ! Tells read_spectral_radiance to find region data (only)
      line    = 0  ! LOS count
      formatted = true      ! For now; how can this be automated?
      allocate (uni_lambda(1), uni_radiance(1))  ! Avoid passing unallocated

      do iline = 1, nlines  ! Process all active line segments

         dialog = mod (iline, 10) == 1
         if (dialog) write (luncrt, '(a, i6)') '   Scanning line', iline

         call read_spectral_radiance &
            (lunrad, formatted, iline, icol, mode, nregion, range, &
             nrecords(iline), irange(:,iline), iuniform, uni_lambda, &
             uni_radiance, ios)
         if (ios /= 0) go to 99

         if (dialog) write (luncrt, '(a, 7i9)') ' irange(:):', irange(:,iline)

         do iregion = 1, nregion
            nlam = irange(iregion+1,iline) - irange(iregion,iline) + 1
            if (nlam > mxregion(iregion)) then
               mxregion(iregion) = nlam
            end if
         end do

      end do  ! Next line of intensities

!     Establish uniform region wavelengths dense enough for all active lines:

      write (luncrt, '(/, a)') ' Max. wavelengths found in each region:'
      write (luncrt, '(i3, i9)') &
         (iregion, mxregion(iregion), iregion = 1, nregion)

      nu = sum (mxregion) - nregion + 1  ! Don't count end point overlaps
      write (luncrt, '(/, a, i9)') ' Total # uniform wavelengths:', nu

      deallocate (uni_lambda, uni_radiance)
      allocate (uni_lambda(nu), uni_radiance(nu))

      write (luncrt, '(/, a)') &
         '      Region         Uniform Indices    # Uniform Wavelengths'
      iuniform(1) = 1
      do iregion = 1, nregion
         nlam = mxregion(iregion)  ! Includes both end points
         iuniform(iregion+1) = iuniform(iregion) + nlam - 1
         call xgrid (nlam, 0, range(iregion), range(iregion+1), &
                     uni_lambda(iuniform(iregion)))
         write (luncrt, '(2f9.1, 2i9, i18)') &
            range(iregion:iregion+1), iuniform(iregion:iregion+1), nlam
      end do

      deallocate (mxregion)
      allocate (lambda_sum(nu))  ! An integration sum for each uniform lambda
      allocate (simple_sum(nu))  ! A sum over all pixels  "   "   "   "   "

 99   return

      end subroutine get_common_wavelengths

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine get_neqair_case (lun, icol)
!
!     Determine the kind of NEQAIR output to be integrated, based on the file
!     names found in /LINE-1 below the working directory.  Wavelengths appear
!     in column 1 of the large *.out files.
!
!     intensity.out           Largest output with spectral radiance in column 3
!     intensity_scanned.out   Reduced output with spectral radiance in column 2
!     radiance.out            One-line file written by the above two cases
!     neqair.out (only)       Original case, with radiance on the last line
!
!     04/16/2018  D.A.Saunders  Avoid a new prompt affecting original case.
!     04/20/2018    "      "    Added icol = 1 case, which means icol = 2 or 3
!                               has already been run, leaving behind one-line
!                               files containing the result of integrating each
!                               line w.r.t. wavelength, giving radiances that
!                               can be integrated to give a heat flux for
!                               comparison with original-type integration of
!                               radiances found in neqair.out.  The presence of
!                               'radiance.out' should override all others.
!
!     04/16/2018  David Saunders  Implemented to avoid a prompt.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

!     Arguments:

      integer, intent (in)  :: lun   ! Logical unit to use
      integer, intent (out) :: icol  ! 3 if intensity.out is found, else
                                     ! 2 if intensity_scanned.out is found, else
                                     ! 1 if radiance.out is found (from 3 or 2);
                                     ! 0 if none of the above
!     Local parameters:

      character (7), parameter :: path = 'LINE-1/'

!     Local variables:

      integer :: i, ios, l
      character (21) :: name

!     Execution:

      do i = 1, 3  ! To favor radiance.out
         select case (i)
            case (1)  ! radiance.out from a case 2 or 3 run
               icol = 1
               l = 12
               name(1:l) = 'radiance.out'
            case (2)  ! intensity_scanned.out
               icol = 2
               l = 21
               name(1:l) = 'intensity_scanned.out'
            case (3)  ! intensity.out'
               icol = 3
               l = 13
               name(1:l) = 'intensity.out'
         end select

         open  (lun, file=path//name(1:l), status='old', iostat=ios)
         close (lun)

         if (ios == 0) exit  ! Else ..
         icol = 0
      end do

      if (icol == 0) then
         write (*, '(/, 2a)') ' Standard case: processing radiances from the', &
                              ' last lines of neqair.out files.'
      else
          write (*, '(/, 2a)') ' Files to be processed: ', name(1:l)
      end if

      end subroutine get_neqair_case

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine get_pixels_and_line_segments ()
!
!     Read the pixels.g file written by Stardust_Lines.  It is a rectangular
!     plane [non]uniformly discretized in each direction, orthogonal to the
!     parallel lines of sight on which NEQAIR has been run.  Its header tells
!     which of those lines (with nk = 2) did not pass through the flow field.
!     The pixel spacings may be nonuniform, so all cells are not the same area.
!
!     The distance from the pixel plane to the sensor location is defined as
!     the given slant range (to the pixel array center coordinates xm, ym, zm).
!     This simplifies calculation of the observer location, px, py, pz.
!     Remember, ym = 0 because we have only half the effective pixel array.
!     Also, the pixels are probably nonuniformly spaced.
!
!     Also, read the lines of sight file, hard-coded here as line_segments.g.
!     Only the header is needed.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: i, j
      real    :: ds, h

      open (lunpix, file='pixels.g', status='old', iostat=ios)
      if (ios /= 0) then
         write (luncrt, '(a)') &
            ' *** Unable to open pixels.g from Stardust_Lines.'
         go to 99
      end if

      call xyzq_read (lunpix, -lunpix, true, npixelblocks, numf, unused, &
                      pixels, ios)
      if (ios /= 0) then
         write (luncrt, '(a)') ' *** Unable to read the pixel data, pixels.g.'
         go to 99
      end if

      niLOS = pixels(1)%ni;  njLOS = pixels(1)%nj
      npixels = niLOS*njLOS

!     By construction, the pixel plane is orthogonal to the line of sight
!     from the sensor location.  Calculate the mid-point, as needed for
!     cosine calculations.  ym = 0. though (on symmetry plane).

      xmax = maxval (pixels(1)%x);  xmin = minval (pixels(1)%x)
      ymax = maxval (pixels(1)%y);  ymin = minval (pixels(1)%y)
      zmax = maxval (pixels(1)%z);  zmin = minval (pixels(1)%z)

      xm = (xmax - xmin)*half;  ym = zero;  zm = (zmax - zmin)*half
      px = xm - D*cosd(angle);  py = zero;  pz = zm - D*sind(angle)

      write (luncrt, '(3x, a, f10.2)') &
         'Viewing angle, deg:    ', angle, &
         'Slant range, m:        ', D
      write (luncrt, '(3x, a, 3es15.7)') &
         'Middle pixel xm,ym,zm: ', xm, ym, zm, &
         'Sensor px,py,pz:       ', px, py, pz

!     Calculate the cone angles subtended at the sensor location by the pixel
!     locations.  Do the cosines (for normal components of radiance at the
!     sensor) here rather than possibly many times later (even though all angles
!     are bound to be extremely small).  Use pixels(1)%q to store them.

      allocate (pixels(1)%q(1,niLOS,njLOS,1))

      biggest_cosine = -one  ! We want the cosine nearest to 1.
                             ! (angle nearest to 0.)
      do j = 1, njLOS
         do i = 1, niLOS
            ds = sqrt ((pixels(1)%x(i,j,1) - xm)**2 + &  ! In pixel plane
                       (pixels(1)%y(i,j,1) - ym)**2 + &
                       (pixels(1)%z(i,j,1) - zm)**2)
            h = sqrt (ds**2 + D**2)     ! Hypotenuse
            pixels(1)%q(1,i,j,1) = D/h  ! Cosine for normal component at sensor
            biggest_cosine = max (biggest_cosine, pixels(1)%q(1,i,j,1))
         end do
      end do

      write (luncrt, '(3x, a, es23.16)') 'Biggest cosine:', biggest_cosine

!     Read the header (only) of the line of sight segments to determine which
!     pixels are inactive.

      open (lunLOS, file='line_segments.g', status='old', iostat=ios)
      if (ios /= 0) then
         write (luncrt, '(a)') &
            ' *** Unable to open line_segments.g from Stardust_Lines.'
         go to 99
      end if

      call xyz_header_io (1, lunLOS, true, npixels, LOS, ios)
      if (ios /= 0) then
         write (luncrt, '(a)') &
            ' *** Unable to read the header of LOS file line_segments.g.'
!!       go to 99
      end if

 99   close (lunLOS)
      return

      end subroutine get_pixels_and_line_segments

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine get_region_data (lun, nregion, range, ios)
!
!     This has been adapted from the NEQAIR_Integration version for V15+ only.
!
!     Determine the wavelength regions for which many lines of sight have been
!     processed by NEQAIR.  Each region of each line needs to be interpolated
!     to a common uniform wavelength grid that is as dense as that of the dens-
!     est line.
!
!     08/27/2018  D.A.Saunders  Date on the NEQAIR_Integration version that is
!                               adapted here.
!     02/26/2020    "     "     Handle the v15+ control file format only.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in)  :: lun         ! Logical unit on which the NEQAIR 
                                           ! control file is opened/closed here
      integer, intent (out) :: nregion     ! # wavelength regions found;
                                           ! nregion = 6 means there are 7
                                           ! region boundaries
      real,    intent (out), pointer :: range(:)
                                           ! Array of wavelengths defining the
                                           ! regions, dimension nregion+1
      integer, intent (out) :: ios         ! 0  means no error was detected;
                                           ! nonzero means failure; error
                                           ! messages are written to std out
!     Local constants:

      integer, parameter :: maxbuf = 80    ! Buffer length
      integer, parameter :: maxreg = 32    ! More than enough, for free form
      logical, parameter :: back = .true.  ! For isolating a last token

!     Local variables:

      integer :: i, lenbuf
      real    :: rregion
      logical :: free_format
      character (maxbuf) :: buffer

!     Procedures:

      logical, external :: number  ! String utility
      external :: upcase           !    "      "

!     Execution:

      open (lun, file='neqair.inp', status='old', iostat=ios)
      if (ios /= 0) then
         write (luncrt, '(a)') &
            ' Unable to open neqair.inp for the region definitions.'
         go to 99
      end if

      do ! Until a REGION line is found or an error is detected
         read (lun, '(a)', iostat=ios) buffer
         if (ios > 0) then
            write (*, '(a)') ' Trouble reading region data from neqair.inp.'
            exit
         else if (ios < 0) then
            write (*, '(a)') ' Region data not found in neqair.inp.'
            exit
         else
            lenbuf = len_trim (buffer)
            call upcase (buffer(1:lenbuf))
            if (index (buffer(1:lenbuf), 'REGION') > 0) exit
         end if
      end do  ! Next input line

!     A blank line ends a list of regions like this:

!     Line 7: Regions
!         855.50   2000.00 M    0.00133 R 600
!        2000.00   5800.00 M    0.00334 R 50
!        5800.00  16000.00 M    0.01135 R 50
!       16000.00  39600.00 M    0.03806 R 25
!       39600.00  60000.00 M    0.03806 R 25
!       60000.00 200000.00 M    0.10000 R 25

      allocate (range(maxreg))  ! More than enough

      nregion = 0
      do  ! Until an empty line is found
         read (lun, '(a)')  buffer
         lenbuf = len_trim (buffer)
         if (lenbuf == 0) exit

         nregion = nregion + 1
         read (buffer(1:lenbuf), *, iostat=ios) range(nregion), rregion
         if (ios /= 0) then
            write (*, '(a, i3)') 'Trouble reading wavelength region', nregion
            go to 99
         end if
      end do

      range(nregion+1) = rregion
      write (*, '(/, a)') ' Wavelength regions found:'
      write (*, '(2f10.2)') (range(i), range(i+1), i = 1, nregion - 1)
      close (lun)

 99   return

      end subroutine get_region_data

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_sensor_details ()  ! TBD
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer, parameter :: line_to_read = 8

!     Initial stab at something plausible:

!     The control file for Stardust_Lines is named lines.inp and
!     line 8 with the viewing angle of elevation above the horizontal
!     now should have a viewing distance in meters following the angle.

      open (lunctl, file='lines.inp', status='old', iostat=ios)
      if (ios /= 0) then
         write (luncrt, '(a)') '*** Trouble opening lines.inp control file.'
         go to 99
      end if

      do i = 1, line_to_read - 1
         read (lunctl, *)  ! Skip the first so many lines
      end do

      read (lunctl, *) angle, distance  ! Degrees and meters
      D = distance                      ! Shorter nomenclature

!     See get_pixels_and_line_segments for more about the sensor details,
!     after some reorganization.

 99   return

      end subroutine get_sensor_details

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine initialize_integrations ()
!
!     Determine the number of meaningful lines/pixels, and back out the (i,j)
!     for each one (may not be needed).
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: i, j, iline

!     Count the number of meaningful lines/pixels:

      nlines = 0
      do iline = 1, npixels
         if (LOS(iline)%nk > 2) nlines = nlines + 1
      end do

      write (luncrt, '(i6, a, i4, a, i4, a, i8)') &
         nlines, ' pixels active out of', niLOS, ' x', njLOS, ' =', npixels

      allocate (ip(nlines), jp(nlines))

      iline = 0;  nlines = 0
      do j = 1, njLOS
         do i = 1, niLOS
            iline = iline + 1
            if (LOS(iline)%nk > 2) then
               nlines = nlines + 1
               ip(nlines) = i
               jp(nlines) = j
            end if
         end do
      end do

      end subroutine initialize_integrations

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine perform_integrations ()
!
!     For each uniform wavelength derived from the wavelengths found in many
!     intensity[_scanned].out files produced by NEQAIR, integrate the spectral
!     radiances w.r.t. solid angle subtended at the indicated sensor location
!     by the underlying rectangular array of quadrilateral pixels cells.
!
!     The straightforward scheme used for single radiance outputs from NEQAIR
!     would require exorbitant amounts of memory to set up all the data ready
!     for integration. Instead, one line is processed at a time and a summation
!     is performed for each (common) wavelength.  Brett Cruden pointed out this
!     more efficient approach.  Consider the 1-D trapezoidal rule analog to
!     understand it more readily.  Outline:
!
!     o  For each cell in the pixel array, calculate the solid angle domega
!        subtended at the sensor location.
!     o  For each cell, the "weight" corresponding to each vertex of the cell
!        is domega/4 * cosine (angle between cell centroid and pixel plane
!        center).
!     o  For each pixel, skip it if the line of sight didn't pass through the
!        flow field (nk = 2), else read the spectral radiances associated
!        with the line of sight and interpolate to the regularized wavelengths.
!     o  Increment the sums for each wavelength by weight*radiance.
!     o  Integrate the regularized radiances w.r.t. wavelength, and write
!        a one-line 'radiances.out' result in that line's subdirectory.
!     o  Save the solid angle integrations for each regularized wavelength
!        as 'integrations.dat'.  These are spectral irradiance files.
!        The values are doubled to account for calculations on a 180-degree
!        revolved axisymmetric flow solution.
!     o  Integrate these w.r.t. wavelength as a check on the heat-flux case.
!        This result is converted from angstroms to microns to give W/cm^2.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local constants:

      real, parameter :: fourth = 0.25

!     Local variables:

      integer :: i, ic, iv, j, jc, mode
      real    :: darea, domega
      character (1), parameter :: method = 'L'  ! Linear quadrature
                                                ! = trapezoidal rule
!     Execution:

      mode = 2  ! Tells read_spectral_radiance to interp. to common lambdas

      lambda_sum(:) = zero  ! One integration sum for each uniform wavelength
      simple_sum(:) = zero  ! One simple sum       "    "    "   "   "   "

      allocate (vertex_weight(niLOS,njLOS));  vertex_weight(:,:) = zero

      do jc = 1, njLOS - 1
         do ic = 1, niLOS - 1 ! Process all pixel cells

            call solid_angle_quad (niLOS, njLOS, pixels(1)%x, pixels(1)%y, &
                                   pixels(1)%z, ic, jc, px, py, pz, domega)

!           One to four quad. cells can affect a vertex/LOS:

            do j = jc, jc+1
               do i = ic, ic+1
                  vertex_weight(i,j) = vertex_weight(i,j) + &
                     (fourth*domega)*pixels(1)%q(1,i,j,1) ! %q = cos cone angle
               end do
            end do

         end do
      end do

      iv   = 0  ! Pixel array vertex count
      line = 0  ! LOS count over all active pixel points/lines of sight
      do j = 1, njLOS
         do i = 1, niLOS
            iv = iv + 1
            if (LOS(iv)%nk > 2) then
               line = line + 1

               call read_spectral_radiance (lunrad, formatted, line, icol, &
                                            mode, nregion, range, &
                                            nrecords(line), &
                                            irange(:,line), iuniform, &
                                            uni_lambda, uni_radiance, ios)
               if (ios /= 0) go to 99

               call save_integrations (i, j)  ! Cumulative for center line only;
                                              ! 'radiance.out' via a further
            else                              ! integration w.r.t. wavelength
               uni_radiance(:) = zero
            end if

            lambda_sum(:) = lambda_sum(:)+vertex_weight(i,j)*uni_radiance(:)
            simple_sum(:) = simple_sum(:)        +           uni_radiance(:)
         end do
      end do

!     Save the solid angle integrations for each uniform wavelength.
!     These are spectral irradiance files, w/cm^2/micron.
!     This has been moved to last below.

!!!   open  (lunrad, file='integrations.dat', status='unknown', iostat=ios)
!!!   write (lunrad, '(f12.4, es12.5)') &
!!!      (uni_lambda(i), lambda_sum(i), i = 1, nuniform)
!!!   close (lunrad)

!     Save the simple summations (x 2), averaged by # active pixels:
!     NO: This file has not been made use of, so suppress it.

!!!   open  (lunrad, file='spectral_irradiance.dat', status='unknown', &
!!!          iostat=ios)
!!!   write (lunrad, '(f12.4, es12.5)') &
!!!      (uni_lambda(i), simple_sum(i)*(2./real (nlines)), i = 1, nuniform)
!!!   close (lunrad)

!     Integrate the spectral irradiances w.r.t. wavelength as a check:

      call lcsquad (nuniform, uni_lambda, lambda_sum, uni_lambda(1), &
                    uni_lambda(nuniform), method, qrad)

      qrad = qrad*angstrom_to_micron  ! Per Brett Cruden, 05/18/2018
      write (luncrt, '(/, (a, es12.5, a))') &
         ' Integration w.r.t. wavelength of solid angle integrals:', qrad, &
         ' W/cm^2', &
         ' Double for 360 degree axisymmetric volume grid:        ', qrad*2., &
         ' W/cm^2'

!     This has been moved to last so that the doubling doesn't affect other
!     results that apply to one line/pixel, not that and its mirror image.

!     Save the solid angle integrations for each uniform wavelength.
!     These are spectral irradiance files, w/cm^2/micron.
!     The values are doubled to account for a 360-degree solution.

      open  (lunrad, file='integrations.dat', status='unknown', iostat=ios)
      write (lunrad, '(f12.4, es12.5)') &
         (uni_lambda(i), 2.*lambda_sum(i), i = 1, nuniform)
      close (lunrad)

 99   return

      end subroutine perform_integrations

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine read_radiance (lun, icol, line, radiance, ios)
!
!     For the indicated line of sight number n = line, read the file named
!     /LINE-n/neqair.out or /LINE-n/radiance.out (depending on icol) and return
!     the radiance value from its last line.  The radiance.out option has been
!     added to enable a check on the spectral integrations where icol = 2 and 3
!     are used to indicate intensity_scanned.out and intensity.out respectively.
!     These cases write a one-line radiance.out file after integrating each line
!     of sight w.r.t. wavelength.  Radiance is the first token on that line.
!
!     NEQAIR v15 complication:  Now, the last line starts with Total, not Final,
!     and a different fixed-format read is needed to get the radiance without
!     parsing the last line.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments (could be inherited, but this is clearer):

      integer, intent (in)  :: lun      ! Logical unit to use
      integer, intent (in)  :: icol     ! 0 => neqair.out; 1 => radiance.out
      integer, intent (in)  :: line     ! LOS number across active pixels
      real,    intent (out) :: radiance ! Result for this LOS from NEQAIR
      integer, intent (out) :: ios      ! Nonzero means a read error

!     Local constants:

      character (11), parameter :: neqair_output_file = '/neqair.out'
      character (13), parameter :: integration_file = '/radiance.out'

!     Local variables:

      integer :: lline
      character (67) :: buffer
      character (32), save :: filename = 'LINE-n'

!     Execution:

      call int_to_char (line, filename(6:), lline)

      if (icol == 0) then
         filename(6+lline:) = neqair_output_file
      else  ! icol = 1
         filename(6+lline:) = integration_file
      end if

      open (lun, file=filename, status='old', iostat=ios)
      if (ios /= 0) then
         write (luncrt, '(a, i6, /, 2a)') &
            ' Trouble opening radiance file for line of sight #', line, &
            ' File name: ', filename
         go to 99
      end if

      if (icol == 0) then  ! Get the last line of neqair.out
         do  ! Until EOF
            read (lun, '(a)', iostat=ios) buffer
            if (ios /= 0) exit
         end do

!        Avoid parsing, etc., as follows.  NEQAIR v14 handling is suppressed.
!        The buffer should look like this for NEQAIR v14:
!           Final radiance value is 1.23456E-01
!        For NEQAIR v15:
!   Total radiance from    855.50 to  39600.00 angstroms =      1.23456 W/cm2-sr

!!!      if (buffer(1:5) == 'Final') then
!!!         read (buffer, '(23x, es12.5)', iostat=ios) radiance
!!!      else  ! buffer(1:5) == 'Total'
            read (buffer, '(54x, f13.5)', iostat=ios) radiance
!!!      end if
      else  ! icol = 1 => radiance is at the start of a one-line file
         read (lun, *, iostat=ios) radiance
      end if
      close (lun)

      if (ios /= 0) then
         write (luncrt, '(a, i6, /, 2a)') &
            ' Trouble reading radiance for line of sight #', line, &
            ' File name: ', filename
         go to 99
      end if

!!!   write (luncrt, '(i6, a, es12.5)') line, ' radiance:', radiance
!!!   Do more at the higher level now.

   99 return

      end subroutine read_radiance

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine read_spectral_radiance (lun, formatted, line, icol, mode, &
                                         nregion, range, nrecords, irange, &
                                         iuniform, uni_lambda, uni_radiance, &
                                         ios)
!
!     For the indicated output from NEQAIR, either scan the line for its
!     contents within each wavelength region (as needed for looking for the most
!     wavelengths over many such outputs) or linearly interpolate the spectral
!     radiances to the indicated number of uniformly spaced wavelengths within
!     each region.
!
!     04/13/2018  D.A.Saunders  Initial design and partial implementation.
!     04/14/2018    "      "    Rest of the initial implementation.
!     04/26/2018    "      "    The formatted intensity[_scanned].out files have
!                               two header lines to be skipped.
!     05/18/2018    "      "    Handle possible duplicate wavelengths at the
!                               interior region boundaries.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in)  :: lun     ! The spectral radiance file is opened
                                       ! and closed here
      logical, intent (in)  :: formatted          ! T/F => ASCII/unformatted
      integer, intent (in)  :: line    ! LOS number indicates the subdirectory
      integer, intent (in)  :: icol    ! Indicator of the file being read;
                                       ! col. 1 contains wavelengths in 2 cases:
                                       ! icol = 2 => intensity_scanned.out;
                                       !      = 3 => intensity.out
      integer, intent (in)  :: mode    ! mode = 1 => return irange(:) indices
                                       !          defining regions for the line;
                                       !      = 2 => interp. each region to the
                                       !          uniform wavelengths indicated
                                       !          by the iuniform(:) indices and
                                       !          generated once beforehand by
                                       !          the calling routine
      integer, intent (in)  :: nregion ! Wavelength boundary indices in the file
      real,    intent (in)  :: range(nregion+1)   ! Region boundary wavelengths
      integer, intent (inout) :: nrecords         ! # wavelengths found in this
                                                  ! intensity file (possibly
                                                  ! with duplicates at interior
                                                  ! region boundaries)
      integer, intent (out) :: irange(nregion+1)  ! Region boundary indices
                                                  ! found for this line when
                                                  ! mode = 1
      integer, intent (in)  :: iuniform(nregion+1)! Region boundary indices for
                                                  ! the desired uniform grid
      real,    intent (in)  :: uni_lambda(*)      ! Uniform wavelengths in each
                                                  ! region generated at the
                                                  ! higher level after the
                                                  ! mode = 1 call
      real,    intent (out) :: uni_radiance(*)    ! Spectral radiances interp-
                                                  ! olated to uniform region
                                                  ! wavelengths
      integer, intent (out) :: ios                ! ios = 0 means no error;
                                                  ! ios > 0 means a read error
!     Local variables:

      integer :: i, lline, ntotal
      real    :: unused
      real, allocatable :: lambda_data(:), radiance_data(:)
      character (32), save :: filename = 'LINE-n'

!     Execution:

      call int_to_char (line, filename(6:), lline)

      if (icol == 2) then
         filename(6+lline:) = '/intensity_scanned.out'
      else  ! icol = 3
         filename(6+lline:) = '/intensity.out'
      end if

      open (lun, file=filename, status='old', iostat=ios)
      if (ios /= 0) then
         write (luncrt, '(a, i6, /, 2a)') &
            ' Trouble opening spectral radiance file for line of sight #', &
            line, ' File name: ', filename
         go to 99
      end if

      select case (mode)

         case (1) ! Scan the line for its region dimensions (only)

!!          if (mod (line, 10) == 1) write (luncrt, '(a, i6)') &
            write (luncrt, '(a, i6)') &
             ' Scanning intensity[_scanned].out wavelengths. Data line #:', line

            call count_records (lun, formatted, nrecords, ios) ! Total input
            if (ios /= 0) then                                 ! wavelengths
               write (luncrt, '(a, i5)') &
                  ' Trouble counting spectral radiances; line:', line
               go to 99
            end if

            if (formatted) then
               nrecords = nrecords - 2  ! Header lines
               allocate (lambda_data(nrecords))
               read (lun, '(a)')
               read (lun, '(a)')
               do i = 1, nrecords
                  read (lun, *, iostat=ios) lambda_data(i)
                  if (ios /= 0) then
                     write (luncrt, '(a, i12)') &
                        ' Trouble reading spectral radiance; line:', i
                     go to 99
                  end if
               end do
            else
               allocate (lambda_data(nrecords))
               do i = 1, nrecords
                  read (lun, iostat=ios) lambda_data(i)
                  if (ios /= 0) then
                     write (luncrt, '(a, i12)') &
                        ' Trouble reading spectral radiance; line:', i
                     go to 99
                  end if
               end do
            end if
            close (lun)

            call suppress1 (nrecords, lambda_data, ntotal)  ! Remove duplicates

            irange(1) = 1
            do i = 2, nregion
               irange(i) = irange(i-1)  ! Starting guess
               call interval (ntotal, lambda_data, range(i), one, irange(i))
            end do
            irange(nregion+1) = ntotal
            deallocate (lambda_data)

         case (2) ! Read the line and interpolate to the indicated uniform
                  ! wavelengths

!!          if (mod (line, 10) == 1) write (luncrt, '(a, i6)') &
            write (luncrt, '(a, i6)') &
               ' Interpolating to common wavelengths. Data line #:', line

            prefix = filename(1:6+lline) ! Path for saving this line's results
            allocate (lambda_data(nrecords), radiance_data(nrecords))

            select case (icol)
               case (2)  ! Spectral radiance is in column 2 of
                         ! intensity_scanned.out
                  if (formatted) then
                     read (lun, '(a)')  ! Skip header lines
                     read (lun, '(a)')
                     read (lun, *, iostat=ios) &
                     (lambda_data(i), radiance_data(i), unused, i = 1, nrecords)
                  else
                     read (lun, iostat=ios) &
                     (lambda_data(i), radiance_data(i), unused, i = 1, nrecords)
                  end if
                  if (ios /= 0) then
                     write (luncrt, '(2a, i12, a, i5)') ' Trouble reading ', &
                        'intensity.scanned.out; expected # wavelengths:', &
                        nrecords, '; line:', line
                     go to 99
                  end if
               case (3)  ! Spectral radiance is in column 3 of intensity.out
                  if (formatted) then
                     read (lun, '(a)')  ! Skip header lines
                     read (lun, '(a)')
                     read (lun, *, iostat=ios) &
                       (lambda_data(i), unused, radiance_data(i), unused, &
                        i = 1, nrecords)
                  else
                     read (lun, iostat=ios) &
                       (lambda_data(i), unused, radiance_data(i), unused, &
                        i = 1, nrecords)
                  end if
                  if (ios /= 0) then
                     write (luncrt, '(2a, i12, a, i5)') ' Trouble reading', &
                        ' intensity.out; expected # wavelengths:', nrecords, &
                        '; line:', line
                     go to 99
                  end if
            end select
            close (lun)

!           Remove possible duplicates at interior region boundaries, in place:

            call suppress2 (nrecords, lambda_data, radiance_data, ntotal)

!           Interpolate each region to the uniform wavelengths found appropriate
!           for all lines analyzed in this run:

            nuniform = iuniform(nregion+1)

            call lcsfit (ntotal, lambda_data, radiance_data, true, 'L', &
                         nuniform, uni_lambda, uni_radiance, uni_radiance)

            deallocate (lambda_data, radiance_data)

         end select

 99   return

      end subroutine read_spectral_radiance

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine save_integrations (i, j)
!
!     Save the results of wavelength-dependent integrations.  For the center
!     pixel line (only), also save cumulative integrals w.r.t. uniform wave-
!     lengths.
!
!     Single line 'radiance.out' files also written here for each active pixel
!     can be processed by another run of this program in its single-radiance
!     integration mode.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in) :: i, j  ! Pixel array indices, for finding
                                    ! the pixel nearest to the center
!     Local variables:

      integer :: iu, l
      real    :: radiance
      character (1), parameter :: method = 'L'  ! Trapezoidal integrations

!     Execution:

      nuniform = iuniform(nregion+1)
      if (pixels(1)%q(1,i,j,1) == biggest_cosine) then 

         allocate (cumulative(nuniform))  ! Cumulative integrals, center only

         call lcsareas (nuniform, uni_lambda, uni_radiance, method, zero, &
                        cumulative)
         open  (lunrad, file='center.line.cumulative.dat', status='unknown', &
                iostat=ios)
         write (lunrad, '(f12.4, es12.5)') &
            (uni_lambda(iu), cumulative(iu), iu = 1, nuniform)
         close (lunrad)
         deallocate (cumulative)
      end if

      call lcsquad (nuniform, uni_lambda, uni_radiance, uni_lambda(1), &
                    uni_lambda(nuniform), method, radiance)

      l = len_trim (prefix)
      open  (lunrad, file=prefix(1:l) // 'radiance.out', status='unknown', &
             iostat=ios)
      radiance = radiance*angstrom_to_micron  ! Per Brett, 05/18/2018
      write (lunrad, '(es12.5, a)') radiance, ' W/cm^2/sr'
      close (lunrad)

      end subroutine save_integrations

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine single_radiance_case ()
!
!     This is the non-spectral case, where NEQAIR has written a single radiance
!     for each active line of sight segment.  It may also be invoked after a
!     spectral case here has saved a single line radiance.out file for each line
!     as a check on the integrations of spectral radiances w.r.t. wavelength.
!     The radiances are saved for plotting on the pixel array, and an integrated
!     result vs. solid angle subtended at the sensor location is printed.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local variables:

      integer :: i, iline, imax, ios, j, jmax, line, lmax, nblocks, nf
      real    :: quadrature(1), radiance, rmax, rnorm

!     Execution:

!     The pixel data structure has only one block and is one k plane.
!     First, gather one radiance from each active line of sight.
!     Strictly speaking, normal components (only) are used, and while the
!     angles involved are extremely small for remote viewing cases, do the
!     cosine calculations anyway.

!     The following enables association of an active line with a pixel (i,j):

!!!   open  (lunout, file='pixel-indices.dat', status='unknown')
!!!   Write this to the standard log file now.

      write (luncrt, '(/, a)') '#  i    j  line        radiance          cosine'

      line = 0;  iline = 0;  rmax = -one
      do j = 1, njLOS
         do i = 1, niLOS
            iline = iline + 1
            if (LOS(iline)%nk == 2) then
               radiance = zero
            else
               line = line + 1
               call read_radiance (lunrad, icol, line, radiance, ios)
               if (ios /= 0) go to 99
            end if

            cosine = pixels(1)%q(1,i,j,1)
            write (luncrt, '(i4, i5, i6, es16.8, es24.16)') &
               i, j, line, radiance, cosine
            rnorm = radiance*cosine
            pixels(1)%q(1,i,j,1) = rnorm  ! Normal component
            if (rnorm > rmax) then
                rmax = rnorm
                imax = i
                jmax = j
                lmax = line
            end if
         end do
      end do

      write (luncrt, '(/, a, es16.8, a, i4, i5, i6)') &
         '# Max. radiance (normal component):', rmax, &
         '  i, j, active line:', imax, jmax, lmax
!!!   close (lunout)

!     Write the pixel radiances (normal components) for plotting:

      open (lunout, file='pixels.f', status='unknown')
      nblocks = 1;  nf = 1  ! Because constants don't work for inout arguments
      pixels(1)%mi = niLOS;  pixels(1)%mj = njLOS;  pixels(1)%mk = 1

      call q_write (lunout, true, npixelblocks, nf, pixels, ios)

!     Integrate the radiances w.r.t. solid angle subtended at the sensor loc.:

      call solid_angle_quadrature_4 (niLOS, njLOS, 1, pixels(1)%x, &
                                     pixels(1)%y, pixels(1)%z, pixels(1)%q, &
                                     px, py, pz, quadrature)

!!!   write (luncrt, '(a, 3es15.7)') ' Sensor location:       ', px, py, pz
      write (luncrt, '(/, (a,  es15.7, a))') &
         ' Quadrature of radiance:', quadrature(1),    ' W/cm^2', &
         ' Double for 360 deg vol:', quadrature(1)*2., ' W/cm^2'
 99   return

      end subroutine single_radiance_case

   end program Stardust_Integration
