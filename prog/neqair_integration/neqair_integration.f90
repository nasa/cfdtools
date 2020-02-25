!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program neqair_integration
!
!  First Words:
!
!  Most of the documentation here applies to the original form that worked with
!  a single output from NEQAIR per line of sight, namely the radiance at a body
!  point integrated over all wavelengths.  This usage is expected to remain the
!  most common, so the original coding is retained as far as possible.  A major
!  extension here treats NEQAIR's more bulky outputs showing wave length regions
!  with spectral radiances computed for many wavelengths.  Either intensity.out
!  or intensity_scanned.out files are treated. Only one of these (not both) is
!  expected in each /LINE-* directory for which NEQAIR has been run, and their
!  presence or absence in /LINE-1 is determined automatically. If neither file
!  is present, the original single-radiance-per-line case is assumed.
!
!  Introduction (Original Form):
!
!  This somewhat specialized application is a companion to HEMISPHERES_OF_SIGHT
!  and NEQAIR_DATA, NEQAIR being a radiative heating solver long used at NASA
!  Ames Research Center.  Each run of NEQAIR operates on one line of sight's
!  worth of real gas flow data to compute the radiance impinging on some body
!  surface point of a hypersonic vehicle in the direction of that line.  Most
!  often, this line is normal to the body, and then a 1-D flow assumption is
!  made to compute a "tangent-slab" approximation to integrating the effect of
!  all possible lines of sight, giving an upper-bound estimate of the radiative
!  heat flux at that body point.  Approximating all possible lines of sight
!  reaching a body point requires much more work.
!
!  Typical atmospheric entry vehicles are relatively simple convex bodies, any
!  point of which can be "seen" by flow field radiation coming from all of the
!  directions defined by a hemisphere centered at the point, tangent to the
!  body.  HEMISPHERES_OF_SIGHT produces a discretization of such a hemisphere
!  at the specified body point(s) in the form of straight lines that have been
!  intersected with the outer boundary of the computational volume grid.  The
!  surface triangulation of points on that outer boundary associated with a
!  given body point, along with NEQAIR radiance results for each line of sight,
!  are gathered by this utility which then performs full angular integration
!  with respect to solid angle to give a more thorough estimate of the radiative
!  heat flux at the body point, the accuracy of which depends on the coarseness
!  of the underlying hemisphere triangulation.
!
!  Further Preliminary Details:
!
!  Consider a run of HEMISPHERES_OF_SIGHT for a pair of body points on the aft
!  parachute cone centerline, using "parachute" as the identifier.  Since both
!  points have y = 0, only half a hemisphere is needed.  This contains two of
!  the underlying discretized quadrants, suitably transformed.  Two sets of
!  three files are produced as follows:
!
!     parachute.1.hemi.dat      ! Transformed/unintersected/2-zone Tecplot file
!     parachute.1.boundary.dat  ! Intersected form of the above triangulation
!     parachute.1.lines.g       ! PLOT3D file of all 1x1xnk intersected lines
!
!     parachute.2.hemi.dat      ! Repeat for body point 2 ...
!     parachute.2.boundary.dat
!     parachute.2.lines.g
!
!  Also, the initial triangulated quadrant centered in (0,0,0) is saved as:
!
!     unit_hemi_quadrant_vertices.dat
!
!  Alternative nomenclature from USLOS or SLOS rather than HEMISPHERES_OF_SIGHT:
!
!  Only one set of hemisphere lines is produced per run of USLOS or SLOS:
!
!     transformed_unit_hemi.dat ! Transformed/unintersected/2-zone triangulation
!                               ! centered on the body point
!     intersected_boundary.dat  ! Intersected form of the above triangulation
!     los.g                     ! 1x1xnk lines of sight in PLOT3D format
!
!  Plottable radiance datasets:
!
!     radiance_normal_components.dat  ! On the intersected boundary surface
!     radiance_normal_hemi.dat        ! On the unintersected hemisphere surface
!
!  The ID.n.hemi.dat files are for visual inspection.  The ID.n.lines.g files
!  contain 2 x nnodes "blocks" of dimensions 1 x 1 x nk, where nnodes is the
!  number of vertices in one quadrant of the triangulated hemisphere.  Each
!  ID.n.lines.g is input to FLOW_INTERP along with temperature and species
!  number density volume data and the associated volume grid (half the full
!  volume in this case, but both halves if any body point is off the y = 0
!  centerline).  The ...boundary.dat files are retrieved here, but NEQAIR_DATA
!  and then NEQAIR have to be run first.
!
!  For body point 1, parachute.1.lines.g and parachute.1.lines.f are passed to
!  NEQAIR_DATA, which produces lines in the form expected by NEQAIR, named
!  LOS-1.dat, LOS-2.dat, ..., LOS-m.dat where m = 2*nnodes in this case, or
!  4*nnodes for an off-center body point.
!
!  NOTE:  Treating more than one body point at a time with HEMISPHERES_OF_SIGHT
!         and subsequent NEQAIR steps is not recommended, so the point number
!         in the nomenclature is highly likely to be 1 (only).  The later USLOS
!         and SLOS eliminate this unnecessary aspect of the nomenclature, but
!         the integration allows for both forms.  Their output lines of sight
!         are written to los.g.
!
!  A shell script can then run NEQAIR on all lines associated with one body
!  point, with its tangent-slab integration step turned off, producing sub-
!  directories /LINE-1, /LINE-2, ..., /LINE-m each containing a neqair.out
!  file with a radiance value on the last line.  (Common edge points are
!  retained in the quadrants to simplify the quadrature performed here, at the
!  cost of some repeated NEQAIR runs.)
!
!  If more than one body point is in the picture, each needs to be in its own
!  subdirectory because the LOS-*.dat files and /LINE-* subdirectory names do
!  not carry the ID along.  See the NOTE above.
!
!  Adjustments For Radiometer Cone Angle < 90 degrees:
!
!  If HEMISPHERES_OF_SIGHT has been used with a cone angle less than 90 deg.
!  entered as a 4th number (say 35.) on the body point x/y/z line, the lines
!  processed by NEQAIR_DATA should be those written to the additional file
!  cone.35.00.parachute.1.lines.g.  These are the only lines for which NEQAIR
!  needs to be run.  The angle calculations used to determine which hemisphere
!  lines to suppress are used here to read NEQAIR results for those lines only,
!  and to zero the radiance values for the remainder.  The quadrature can then
!  proceed over the whole hemispherical triangulation as before.
!
!  Provision is also made for examining a cone angle smaller than the cone angle
!  for which NEQAIR runs are in hand.  However, note that the smaller the cone
!  angle, the denser the underlying discretization can and should be for reason-
!  able resolution.
!
!  NEQAIR_INTEGRATION Procedure (Original Form; Cone Angle = 90 degrees):
!
!  For one body point, at the level above subdirectories /LINE-1, /LINE-2, ...
!  the user of this program enters the relevant ID and body point number (say
!  parachute and 1).  File parachute.1.boundary.dat is read in a loop over (2)
!  zones.  For each zone, nnodes lines are processed by extracting a radiance
!  value from the last line of the appropriate /LINE-*/neqair.out file into
!  the %f(1,:) elements of the zone/quadrant triangulation's function array
!  to match its %xyz(1:3,:) coordinates, followed with a call to low-level
!  subroutine solid_angle_quadrature_3 to produce a measure of radiative heat
!  flux at the body point for that quadrant/zone.  The sum over all quadrants
!  for one body point is the final heat flux for that body point.
!
!  Extensions For Distinguishing Wavelength Dependencies:
!
!  Wavelength-dependent integrations w.r.t. solid angle are saved as the file
!  integrations.dat, along with one radiance per line further described below.
!  The main complication is that storing all the data in intensity.out or
!  intensity_scanned.out could require vast amounts of memory. Ensuring a
!  common set of wavelengths for all regions of all lines, and storing only
!  one line's worth of data at a time is far more efficient (thanks to Brett
!  Cruden for this insight).  The added complexity is isolated as much as
!  possible so that the original scheme for which solid_angle_quadrature_3 was
!  written can be retained.  Both the original and the extended schemes allow
!  radiometer cone angles less than 90 (meaning NEQAIR output data are present
!  for only a subset of hemisphere lines) and possibly further analyses with
!  cone angles specified as smaller than the cone angle of the data.  The same
!  choice of zeroing function values for lines outside the indicated cone angle
!  simplifies the integrations in all cases. For each line, an integration
!  w.r.t. wavelength produces a radiance saved as /LINE-*/radiance.out, and
!  these can be integrated w.r.t. solid angle as a check that should match the
!  original integration scheme for single radiances produced by NEQAIR.
!  In /LINE-1 (only) the file line-1.cumulative.dat contains the partial
!  integrals w.r.t. wavelength for each common uniform wavelength.
!  Any radiance.out files in /LINE-* should be deleted before a run intended
!  to process intensity[_scanned].out files, because the program automates
!  the type of integration performed according to the files found in /LINE-1.
!
!  History:
!
!  03/28/2014  D.A.Saunders  Initial design.
!  03/31/2014    "     "     Initial testing.
!  04/02/2014    "     "     Brett Cruden pointed out that only the normal
!                            component of radiance is seen by the body point.
!  02/19/2015    "     "     The screen heat flux dialogue can be partly lost
!                            when there are many lines, so write it to a file,
!                            results.dat.
!  02/23/2015    "     "     Save the radiances with the intersected surface
!                            triangulation in contourable form, as the file
!                            radiance_normal_components.dat.
!  02/24/2015    "     "     Contours of radiance (normal component) on the
!                            intersected surface are hard for some people to
!                            interpret, so do likewise for the hemispherical
!                            surface as well.  This has been transformed to
!                            be centered on the body point, but it can be
!                            rotated during the plotting of the contours.
!                            The file name is radiance_normal_hemi.dat.
!                            (The need is to help design radiometers that
!                            might be flown on future missions.)
!  02/28/2015    "      "    Clarified some of the description above.
!  03/03/2015    "      "    The twiddling of %numf wasn't quite right.
!  02/16/2018    "      "    Extended to handle less-than-full-hemisphere data
!                            as needed for radiometer studies by Todd White.
!  02/21/2018    "      "    Allow for integrating over a smaller cone angle
!                            than that for which the line of sight radiances
!                            are in hand (but keep the data resolution in mind).
!  04/19/2018    "      "    Start of extensions to allow study of wavelength
!                            dependencies as further needed for radiometers.
!  04/24/2018    "      "    Lots of head-scratching later: ready for testing.
!  05/17/2018    "      "    (After a hiatus:) Several minor glitches later, it
!                            looks promising, although the results are four
!                            orders of magnitude too high.
!  05/18/2018    "      "    Integrate the wavelength quadratures over solid
!                            angle w.r.t. wavelength as a check.  Brett says
!                            that radiances computed vs. angstroms should be
!                            vs. microns, meaning scaling by 1.e-4.  Worse
!                            issue: the intensity[_scanned].out files from
!                            NEQAIR v14 have duplicate abscissas at the interior
!                            region boundaries.  For now, look for these and
!                            suppress them.
!  05/25/2018    "      "    NEQAIR v15 writes radiance results differently.
!                            Look for Final or Total on the last line of
!                            neqair.out for v14 and v15 respectively, and use
!                            two different formatted reads (to avoid parsing).
!  05/31/2018    "      "    A desired cone angle of 10 deg. didn't work when
!                            the data angle was 90 deg.
!  08/27/2018    "      "    The reading of wavelength regions needed to handle
!                            NEQAIR V15's free format.
!  08/30/2018    "      "    Add sqrt (eps) to the data cone angle to allow for
!                            noise in the angle_between_vectors results.
!  08/31/2018    "      "    Even 10*sqrt (eps) is not enough for a Mars 2020
!                            radiometer case.  Go to 0.1 deg.
!  09/05/2018    "      "    Print doubled W/cm^2 at the end of screen output
!                            for centerline/half-hemisphere cases (2 zones).
!  09/07/2018    "      "    Integration of radiance.out results with desired
!                            cone angle < data cone angle needed a fix.
!  11/09/2018    "      "    The advent of USLOS (unstructured surfaces) and
!                            SLOS (structured grid analog of USLOS) required
!                            handling of simplified nomenclature for the various
!                            files. If the simpler los.g file name is not found,
!                            the original nomenclature is handled instead.
!  02/24/2020    "      "    The vertex weights used in the perform_integrations
!                            procedure have not been initialized to zero until
!                            now, although it doesn't seem to have hurt so far.
!
!  Author:  David Saunders, ERC, Inc./NASA Ames Research Center, CA
!                  Now with AMA, Inc. at NASA ARC.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_block_structure  ! Employed by the XYZQ_IO package
   use xyzq_io_module        ! PLOT3D-type I/O package
   use string_justify        ! Character manipulation module
   use tri_header_structure  ! For multizone unstructured Tecplot files
   use tri_zone_structure    ! For one zone of an unstructured Tecplot file
   use triangulation_io      ! I/O package for multizone surface triangulations

   implicit none

!  Local constants:

   integer, parameter :: &
      lunlos = 1,        &   ! Lines of sight file, ID.n.lines.g, for BP xyz
      luntri = 2,        &   ! Input triangulated surface,  ID.n.boundary.dat
      luntr2 = 3,        &   ! Unintersected triangulation, ID.n.hemi.dat
      lunrad = 4,        &   ! Files containing radiance results from NEQAIR
      lunkbd = 5,        &   ! Keyboard inputs
      luncrt = 6,        &   ! Screen prompts and diagnostics
      lunout = 7,        &   ! For results.dat file summary (one radiance/line)
      luncon = 8,        &   ! For writing retrieved radiances in Tecplot form
      lunhem = 9             ! For writing radiances with the hemisphere srfc.

   real, parameter :: &
      angle_margin = 0.1,         &  ! Deg., to ensure being inside cone angle
      angstrom_to_micron = 1.e-4, &  ! Brett says wavelengths need to be microns
      ninety = 90.,               &
      zero   = 0.

   logical, parameter :: &
      false = .false.,   &
      true  = .true.

!  Local variables:

   integer :: &
      iactive, ibp, icol, icone, ios, izone, lbp, lid, line, lname, lun, nf, &
      nk, nlines, node, nnodes, ntris, nzones

   real :: &
      cone_angle, cos_theta, data_cone_angle, desired_angle, qrad, radiance, &
      total_heat_flux

   real, dimension (3) :: &
      bpxyz, v1, v2

   real, allocatable :: &
      zone_heat_flux(:)

   logical :: &
      cr, eof, formatted, simplified, solid_angle_test

   character (1) :: &
      answer

   character (10) :: &
      prefix        ! Path for saving wavelength-dependent results

   character (64) :: &
      filename

   character (32) :: &
      varname(4)     ! x,y,z + one function for output radiance components

   character (64) :: &
      identifier

   type (grid_type), pointer, dimension (:) :: &
      lines_of_sight

   type (tri_header_type) :: &
      tri_header

   type (tri_type), pointer, dimension (:) :: &
      intersected_hemi, unintersected_hemi

!  Variables introduced for analyzing wavelength effects.  See internal
!  procedure analyze_wavelengths, which unfortunately cannot have its own
!  internal procedures.

   integer :: nregion, nuniform
   integer, allocatable, dimension (:) :: iuniform     ! Unif. region indices
   integer, allocatable, dimension (:,:) :: irange     ! Indices for all hemi
                                                       ! lines; only the data
                                                       ! lines are assigned
   real,    allocatable, dimension (:) :: cumulative   ! Cumulative integrals
                                                       ! w.r.t. wavelength, for
                                                       ! possible plotting
                                                       ! (LINE-1 only)
   real,    allocatable, dimension (:) :: lambda_sum   ! One sum for each wave-
                                                       ! length of each uniform-
                                                       ! ly-spaced region
   integer, allocatable, dimension (:) :: nrecords     ! # wavelengths found in
                                                       ! each input intensity
                                                       ! file (with possible
                                                       ! duplicates at interior
                                                       ! region boundaries)
   real,    allocatable, dimension (:) :: uni_lambda   ! Suits all data lines
   real,    allocatable, dimension (:) :: uni_radiance ! Interpolated spectral
                                                       ! radiances for one data
                                                       ! line
   real,    allocatable, dimension (:) :: vertex_weight! One weight for each
                                                       ! vertex, from one or
                                                       ! more triangles
   real,    pointer,     dimension (:) :: range        ! Region boundary wave-
                                                       ! lengths, dim. nregion+1

!  Execution:
!  ----------

!  The advent of USLOS and SLOS was an opportunity to get away from the unlikely
!  possibility that allowed for more than one hemisphere per run of HEMISPHERES_
!  OF_SIGHT.  Look for simplified notation, but allow for the original.

   open  (lunlos, file='los.g', status='old', iostat=ios)
   simplified = ios == 0

   if (.not. simplified) then
      write (luncrt, '(/, a)', advance='no') ' Body point ID and number: '
      read  (lunkbd, *) identifier, ibp

      lid = len_trim (identifier) + 1
      identifier(lid:lid) = '.'
      tri_header%filename = identifier(1:lid)

      call int_to_char (ibp, tri_header%filename(lid+1:), lbp)

      tri_header%filename(lid+lbp+1:) = '.lines.g'

      open (lunlos, file=tri_header%filename, status='old', iostat=ios)
      if (ios /= 0) then
         write (luncrt, '(/, 2a)') ' Unable to open ', tri_header%filename
         go to 99
      end if
   end if

   call xyz_header_io (1, lunlos, true, nlines, lines_of_sight, ios)
   if (ios /= 0) go to 99

   call xyz_allocate (lines_of_sight(1), ios)
   if (ios /= 0) go to 99

   nk = lines_of_sight(1)%nk

   call xyz_block_io (1, lunlos, true, nk, lines_of_sight(1)%x, &
                      lines_of_sight(1)%y, lines_of_sight(1)%z, ios)
   if (ios /= 0) go to 99

   bpxyz(1) = lines_of_sight(1)%x(1,1,1)
   bpxyz(2) = lines_of_sight(1)%y(1,1,1)
   bpxyz(3) = lines_of_sight(1)%z(1,1,1)

   write (luncrt, '(a, 3f14.8)') ' Body pt. coordinates found:', bpxyz(:)

!  Make the body point normal identical to that in HEMISPHERES_OF_SIGHT:

   v1(1)  = lines_of_sight(1)%x(1,1,nk) - bpxyz(1)
   v1(2)  = lines_of_sight(1)%y(1,1,nk) - bpxyz(2)
   v1(3)  = lines_of_sight(1)%z(1,1,nk) - bpxyz(3)

   deallocate (lines_of_sight(1)%x, lines_of_sight(1)%y, lines_of_sight(1)%z)
   deallocate (lines_of_sight)
   close (lunlos)

!  Option for radiometer studies to work with less than a full hemisphere:

   data_cone_angle = ninety
   call readr (luncrt, &
               'Cone angle of the NEQAIR data? [<cr> = 90. = full hemi.]: ', &
               lunkbd, data_cone_angle, cr, eof)

   desired_angle = data_cone_angle
   do  ! Until a valid cone angle is indicated
      call readr (luncrt, &
         'Desired cone angle (equal or smaller); [<cr> = data angle]: ', &
         lunkbd, desired_angle, cr, eof)
      if (desired_angle <= data_cone_angle) exit
   end do

!  Handle likely noise in the angle-between-vectors calculations.
!  We don't want data lines that are there to seem outside the cone angle.

   data_cone_angle = data_cone_angle + angle_margin
   desired_angle   = desired_angle   + angle_margin

!  Read the header of the triangulation formed by the LOS intersections with
!  the underlying quadrants of a hemisphere (2 or 4 zones):

   if (simplified) then
      tri_header%filename = 'intersected_boundary.dat'
      lname = 24
   else
      tri_header%filename(lid+lbp+1:) = '.boundary.dat'
      lname = lid + lbp + 13
   end if

   tri_header%fileform  = 1        ! Vertex-centered
   tri_header%formatted = true
   tri_header%nvertices = 3        ! Triangles, not tets
   ios = 1                         ! Turn on verbose option

   call tri_header_read (luntri, tri_header, intersected_hemi, ios)
   if (ios /= 0) go to 99

   nzones = tri_header%nzones

   write (luncrt, '(3a, i2)') ' # zones found in ', &
      tri_header%filename(1:lname), ':', nzones

!  Determine whether we'll be reading neqair.out (only choice originally) or
!  intensity.out (icol = 3) or intensity_scanned.out (icol = 2) or radiance.out
!  (icol = 1, written if icol = 3 or 2 on a prior run).

   call get_neqair_case (lunrad, icol)

   if (icol > 1) then
      call analyze_wavelengths ()  ! Avoid obscuring the original scheme, but
      go to 99                     ! that scheme is reused if icol = 1
   end if

!  Original scheme for a single radiance on the last line of neqair.out, now
!  reused for checking the extended scheme by processing the radiance.out files
!  that the extended scheme leaves behind in each LINE-* directory:

   allocate (zone_heat_flux(nzones))

!  Allow for testing the solid angle utilities by making %f = constant (1.):

   write (luncrt, '(a, /, a)', advance='no') &
      ' Option to check solid angle calculation:', &
      ' Just set f = 1. and skip processing of NEQAIR results? [y|n]: '
   read  (lunkbd, *)  answer
   solid_angle_test = answer == 'y' .or. answer == 'Y'
   write (luncrt, '(a)')
   open  (lunout, file='results.dat', status='unknown', position='append')
   write (lunout, '(a)')

!  The number of functions read with the triangulations is assumed to be 0, but
!  we want one to store radiances in %f, so we allocate space for each quadrant/
!  zone (and deallocate) as we process one zone at a time:

   nf = 1  ! %numf = 0 from the triang. reads; we allocate output %f explicitly

   if (.not. solid_angle_test) then

!     We need the unintersected hemisphere surface as well:

      if (simplified) then
         tri_header%filename = 'transformed_unit_hemi.dat'
      else
         tri_header%filename(lid+lbp+1:) = '.hemi.dat'
      end if

      call tri_header_read (luntr2, tri_header, unintersected_hemi, ios)
      if (ios /= 0) go to 99

!     Save radiance results in two plottable forms:

      tri_header%numf = nf
      tri_header%filename = 'radiance_normal_components.dat'
      varname(1:3) = tri_header%varname(1:3)
      varname(4)   = 'Radiance_normal'
      deallocate (tri_header%varname)
      allocate (tri_header%varname(4))
      tri_header%varname(:) = varname(:)

      call tri_header_write (luncon, tri_header, intersected_hemi, ios)
      if (ios /= 0) go to 99

      tri_header%filename = 'radiance_normal_hemi.dat'
      call tri_header_write (lunhem, tri_header, unintersected_hemi, ios)
      if (ios /= 0) go to 99

   end if

   line    = 0  ! LOS count of all hemisphere zones/all nodes
   icone   = 0  ! LOS count of lines within the data cone angle
   iactive = 0  ! LOS count of lines within the desired cone angle

   do izone = 1, nzones  ! Process all vertices of each zone

      tri_header%numf = 0  ! Modified below, so zero it for each zone read
      call tri_zone_allocate (tri_header, intersected_hemi(izone), ios)
      if (ios /= 0) go to 99

      call tri_zone_read (luntri, tri_header, intersected_hemi(izone), ios)
      if (ios /= 0) go to 99

      nnodes = intersected_hemi(izone)%nnodes
      ntris  = intersected_hemi(izone)%nelements

      allocate (intersected_hemi(izone)%f(nf,nnodes))

      if (solid_angle_test) then
         intersected_hemi(izone)%f(:,:) = 1.
      else

         call tri_zone_allocate (tri_header, unintersected_hemi(izone), ios)
         if (ios /= 0) go to 99

         call tri_zone_read (luntr2, tri_header, unintersected_hemi(izone), ios)
         if (ios /= 0) go to 99

         allocate (unintersected_hemi(izone)%f(nf,nnodes))

         do node = 1, nnodes

            line = line + 1  ! Counts nodes of all zones
            v2(:) = intersected_hemi(izone)%xyz(:,node) - bpxyz(:)

            call angle_between_vectors (v1, v2, cone_angle)

            if (cone_angle <= data_cone_angle) then
               icone = icone + 1  ! Counts lines within the data cone angle

!              Allow for looking at a cone angle smaller than that of the
!              NEQAIR results:

               if (cone_angle <= desired_angle) then
                  iactive = iactive + 1  ! Counts lines within desired angle
                  call read_radiance (lunrad, icol, icone, radiance)
                  if (ios /= 0) go to 99

!                 The body point sees only the normal component.
!                 Happily, los 1 is the body-normal line in each quadrant:

                  cos_theta = cosd (cone_angle)
                  unintersected_hemi(izone)%f(1,node) = radiance*cos_theta
                    intersected_hemi(izone)%f(1,node) = radiance*cos_theta
               else
                  unintersected_hemi(izone)%f(1,node) = zero
                    intersected_hemi(izone)%f(1,node) = zero
               end if
            else  ! No data for this node
               unintersected_hemi(izone)%f(1,node) = zero
                 intersected_hemi(izone)%f(1,node) = zero
            end if

         end do  ! Next zone vertex

!        Save the normal components of radiance in two plottable forms:

         tri_header%numf = nf
         call tri_zone_write (luncon, tri_header, intersected_hemi(izone), ios)
         if (ios /= 0) go to 99
         call tri_zone_write (lunhem, tri_header, unintersected_hemi(izone),ios)
         if (ios /= 0) go to 99
      end if

      call solid_angle_quadrature_3 (nnodes, ntris, nf, &
                                     intersected_hemi(izone)%xyz , &
                                     intersected_hemi(izone)%f,    &
                                     intersected_hemi(izone)%conn, &
                                     false, bpxyz, zone_heat_flux(izone))

      deallocate (intersected_hemi(izone)%xyz, intersected_hemi(izone)%f, &
                  intersected_hemi(izone)%conn)
      if (.not. solid_angle_test) &
         deallocate (unintersected_hemi(izone)%xyz, &
                     unintersected_hemi(izone)%f,   &
                     unintersected_hemi(izone)%conn)

      do lun = lunout, luncrt, luncrt - lunout
         if (izone == 1) &
            write (lun, '(a, f24.2, a)') &
               ' Data cone angle:  ', data_cone_angle, ' degrees', &
               ' Integration angle:', desired_angle,   ' degrees'
         if (solid_angle_test) then
            write (lun, '(a, i2, a, f36.4, a)') &
               ' zone', izone, ':', zone_heat_flux(izone), ' sr'
         else
            write (lun, '(a, i2, a, f36.4, a)') &
               ' zone', izone, ':', zone_heat_flux(izone), ' W/cm^2'
         end if
      end do

   end do  ! Next triangulation zone

   if (.not. solid_angle_test) write (luncrt, '(/, a, i6)') &
      ' # lines found within the desired cone angle:', iactive

   total_heat_flux = sum (zone_heat_flux)

   do lun = lunout, luncrt, luncrt - lunout
      if (solid_angle_test) then
         write (lun, '(a, f25.4, a)') &
            ' Total solid angle:', total_heat_flux, ' sr'

         if (nzones == 2) then  ! Body point was on the y = 0. symmetry plane
            write (lun, '(a, f15.4, a)') &
               ' Double for y = 0 body point:', 2.*total_heat_flux, ' sr'
         end if
         go to 99  ! Avoid treating smaller cone angles
      else
         write (lun, '(a, f27.4, a)') &
            ' Total heat flux:', total_heat_flux, ' W/cm^2'

         if (nzones == 2) then  ! Body point was on the y = 0. symmetry plane
            write (lun, '(a, f15.4, a)') &
               ' Double for y = 0 body point:', 2.*total_heat_flux, ' W/cm^2'
         end if
      end if
   end do

99 continue

!  Internal procedures for program neqair_integration, in alphabetical order:

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine analyze_wavelengths ()
!
!     This major extension to neqair_integration performs its original function
!     for each wavelength of each region analyzed by NEQAIR.  It is isolated to
!     leave the original scheme essentially intact, although that scheme can be
!     reused following analysis of intensity_scanned.out or intensity.out files
!     by processing the one-line radiance.out files that are written at the end
!     of the more elaborate analysis by integrating w.r.t. wavelength for each
!     line of sight.
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine get_common_wavelengths ()
!
!     For the type of output from NEQAIR indicated by icol, scan all relevant
!     line results for the wavelengths used in each region, and calculate a
!     single set of wavelengths that are uniform in each region and as dense
!     as the densest line found in each region.  The spectral radiances from
!     each line will later be interpolated to this common set of wavelengths.
!
!     The high-level processing of all relevant lines (with radiometer cone
!     angles in the picture) is duplicated from the original radiance-only
!     scheme, and it has to be repeated following this preprocessing step.
!     Store vertex cone angles as a single vertex-centered function value.
!
!     04/16/2018  David Saunders  Initial implementation (internal procedure).
!     05/18/2018    "      "      Suppress duplicate wavelengths at the region
!                                 boundaries.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local variables:  (Most variables are inherited from the main program.)

      integer :: iregion, mode, nlam, nu
      integer, allocatable :: mxregion(:)

!     Execution:

      allocate (mxregion(nregion))  ! For maximum region sizes found
      mxregion(:) = 0

!     We need to know the total number of data lines (= # triangulation vertices
!     within the appropriate cone angle), so we have to read at least zone 1
!     before allocating wavelength range workspace.  We might as well read all
!     zones outside the loop over zones.

      ios = 1  ! Verbose option
      close (luntri)  ! The header has been read earlier.
      call tri_read (luntri, tri_header, intersected_hemi, ios)
      if (ios /= 0) go to 99

      nnodes = intersected_hemi(1)%nnodes  ! Same for all zones

      allocate (irange(nregion+1,nnodes*nzones), & ! Regions found in data lines
                iuniform(nregion+1),             & ! Regions imposed on all data
                nrecords(nnodes*nzones))           ! For data lines (possibly
                                                   ! containing duplicates at
                                                   ! interior region boundaries)

      mode    = 1  ! Tells read_spectral_radiance to find region data (only)
      line    = 0  ! LOS count over all hemisphere zones
      icone   = 0  ! LOS count over lines within the data cone angle
      iactive = 0  ! LOS count over lines within the active cone angle specified
      formatted = true      ! For now; how can this be automated?
      allocate (uni_lambda(1), uni_radiance(1))  ! Avoid passing unallocated

      do izone = 1, nzones  ! Process all vertices of each triangulated zone
!        nnodes = intersected_hemi(izone)%nnodes
         ntris  = intersected_hemi(izone)%nelements
         deallocate (intersected_hemi(izone)%f)     ! Because I/O pkg. has
         allocate (intersected_hemi(izone)%f(1,nnodes))  ! allocated (1,1)

         do node = 1, nnodes
            line = line + 1
!!          if (mod (line, 10) == 1) write (luncrt, '(2a, i6)') &

            v2(:) = intersected_hemi(izone)%xyz(:,node) - bpxyz(:)

            call angle_between_vectors (v1, v2, cone_angle)

            write (luncrt, '(/, a, i6, es25.15)') &
               ' Looking at boundary triangln. line #:', line, cone_angle

            intersected_hemi(izone)%f(1,node) = cone_angle

            if (cone_angle <= data_cone_angle) then
               icone = icone + 1  ! icone counts the lines of NEQAIR output
               write (luncrt, '(a, i6)') ' Checking data line #', icone
               if (cone_angle <= desired_angle) then
                  iactive = iactive + 1  ! Allow for analyzing smaller angles

                  write (luncrt, '(a, i6, 2es25.15)') &
                     ' Processing active line #', iactive, cone_angle, &
                          desired_angle

                  call read_spectral_radiance &
                     (lunrad, formatted, icone, icol, mode, nregion, &
                      range, nrecords(icone), irange(:,icone), iuniform, &
                      uni_lambda, uni_radiance, ios)
                  if (ios /= 0) go to 99

                  write (luncrt, '(a, 7i9)') ' irange(:):', irange(:,icone)

                  do iregion = 1, nregion
                     nlam = irange(iregion+1,icone) - irange(iregion,icone) + 1
                     if (nlam > mxregion(iregion)) then
                        mxregion(iregion) = nlam
                     end if
                  end do
               end if
            end if

         end do  ! Next zone node
      end do  ! Next triangulation zone

!     Establish uniform region wavelengths dense enough for all active lines:

      write (luncrt, '(/, a, i6)') &
         ' # lines found within the desired cone angle:', iactive
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine get_region_data (lun, nregion, range, ios)
!
!     As part of radiometer analysis, determine the wavelength regions for which
!     many lines of sight have been processed by NEQAIR.  Each region of each
!     line needs to be interpolated to a common uniform wavelength grid that is
!     as dense as that of the densest line.
!
!     04/12/2018  D.A.Saunders  A necessary first step, for neqair.inp in the
!                               old formatted form. Deal with free format later.
!     08/27/2018  D.A.Saunders  Handle both old and new input format by looking
!                               for the first occurrence of "REGION" (upcased).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

      read (lun, '(a)') buffer
      free_format = number (buffer)  ! Version number is on line 1 from V15 on.

!     Either way, look for the first occurrence of "REGION" after upcasing:

      do ! Until the line is found or an error is detected
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

      if (free_format) then  ! A blank line ends the list of regions

!        Line 7: Regions
!            855.50   2000.00 M    0.00133 R 600
!           2000.00   5800.00 M    0.00334 R 50
!           5800.00  16000.00 M    0.01135 R 50
!          16000.00  39600.00 M    0.03806 R 25
!          39600.00  60000.00 M    0.03806 R 25
!          60000.00 200000.00 M    0.10000 R 25

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
!!       range(nregion+1) = rregion  ! Common to both cases

      else                   ! Old fixed-format control file

!        Read the number of regions as the last token on a line like this:
!        REGION DATA:        Line8                     # of regions =  6
!                                                           iii
!           w1 [A]   w2 [A]   range   grid_type  delta_lambda   pointsPerLine
!           855.5   2000.0     600        0        0.00133          10
!          2000.0   5800.0      50        0        0.00334          10
!          5800.0  16000.0      50        0        0.01135          10
!         16000.0  39600.0      25        0        0.03806          10
!         39600.0  60000.0      25        0        0.03806          10
!         60000.0  200000.      25        0        0.10000          10

         i = index (buffer(1:lenbuf), ' ', back) + 1
         write (*, '(a, 2i4)') ' i, lenbuf:', i, lenbuf
         call decode_number (buffer(i:lenbuf), nregion, rregion, ios)
         if (ios /= 2) then  ! Neither integer nor real found
            ios = 1  ! Failure
            write (*, '(a)') 'Bad region count found in this line:', &
               buffer(1:lenbuf)
         end if

         allocate (range(nregion+1))

!        Scan further for the region boundaries:

         do i = 1, 2
            read (lun, '(a)')
         end do

         do i = 1, nregion
            read (lun, *, iostat=ios) range(i), rregion
            if (ios /= 0) then 
               write (*, '(a, i3)') 'Trouble reading region line', i
               go to 99
            end if
         end do

      end if  ! neqair.inp format choice

      range(nregion+1) = rregion
      write (*, '(/, a)') ' Wavelength regions found:'
      write (*, '(2f10.2)') (range(i), range(i+1), i = 1, nregion - 1)
      close (lun)

 99   return

      end subroutine get_region_data

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine perform_integrations ()
!
!     For each uniform wavelength derived from the wavelengths found in many
!     intensity[_scanned].out files produced by NEQAIR, integrate the spectral
!     radiances w.r.t. solid angle, probably over the hemisphere triangulation
!     elements contained within an active radiometer cone angle that may be
!     equal to or smaller than the data cone angle used in HEMISPHERES_OF_SIGHT
!     to determine the lines of sight for which NEQAIR is run.
!
!     The original scheme used for single radiance outputs from NEQAIR would
!     require exorbitant amounts of memory to set up all the function data ready
!     for integration. Instead, one line is processed at a time and a summation
!     is performed for each (common) wavelength.  Brett Cruden pointed out this
!     more efficient approach.  Consider the 1-D trapezoidal rule analog to
!     understand it more readily. The cell-centered function values are taken to
!     be the average of triangle vertex values.  Only the body-normal components
!     of radiance are integrated.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local constants:

      real, parameter :: third = 1./3.

!     Local variables:

      integer :: i, iv, itri, iz, mode
      real    :: darea, domega
      character (1), parameter :: method = 'L'  ! Linear quadrature

!     Execution:

      mode    = 2  ! Tells read_spectral_radiance to interp. to common lambdas
      line    = 0  ! LOS count over all hemisphere zones/all nodes
      icone   = 0  ! LOS count over lines within the data cone angle
      iactive = 0  ! LOS count over lines within the active cone angle specified

      lambda_sum(:) = zero  ! One for each uniform wavelength

!     The triangulated zones have been read and stored by get_common_wavelengths
!     and nnodes, ntris are common to all zones.  Cone angles have also been
!     stored in %f(1,node).

      allocate (vertex_weight(nnodes*nzones));  vertex_weight(:) = zero

      do izone = 1, nzones  ! Process all vertices of each triangulated zone

         iz = (izone - 1)*nnodes

         do itri = 1, ntris
            call solid_angle_tri (nnodes, ntris, intersected_hemi(izone)%xyz, &
                                  intersected_hemi(izone)%conn, itri, bpxyz,  &
                                  darea, domega)

!           One to six triangles can affect a vertex/node/LOS:

            do i = 1, 3
               iv = intersected_hemi(izone)%conn(i,itri)
               vertex_weight(iv+iz) = vertex_weight(iv+iz) + & ! %f = cone ang.
                  (third*domega) * cosd (intersected_hemi(izone)%f(1,iv))
            end do
         end do

         do node = 1, nnodes
            line = line + 1
            cone_angle = intersected_hemi(izone)%f(1,node) ! See get_common_wav.

            if (cone_angle <= data_cone_angle) then
               icone = icone + 1  ! icone counts the lines of NEQAIR output
               if (cone_angle <= desired_angle) then
                  iactive = iactive + 1  ! Allow for analyzing smaller angles

                  call read_spectral_radiance (lunrad, formatted, icone, icol, &
                                               mode, nregion, range, &
                                               nrecords(icone), &
                                               irange(:,icone), iuniform, &
                                               uni_lambda, uni_radiance, ios)
                  if (ios /= 0) go to 99

                  call save_integrations ()  ! Cumulative for line 1 only, else
                                             ! 'radiance.out' via a further
               else                          ! integration w.r.t. wavelength
                  uni_radiance(:) = zero
               end if
               lambda_sum(:) = lambda_sum(:)+vertex_weight(line)*uni_radiance(:)
            end if

         end do  ! Next zone node

      end do  ! Next triangulation zone

!     Save the solid angle integrations for each uniform wavelength:

      open  (lunrad, file='integrations.dat', status='unknown', iostat=ios)
      write (lunrad, '(f9.1, es12.5)') &
         (uni_lambda(i), lambda_sum(i), i = 1, nuniform)
      close (lunrad)

!     Integrate these results w.r.t. wavelength as a check:

      call lcsquad (nuniform, uni_lambda, lambda_sum, uni_lambda(1), &
                    uni_lambda(nuniform), method, qrad)

      qrad = qrad*angstrom_to_micron  ! Per Brett, 05/18/2018
      write (luncrt, '(/, a, es12.5, a)') &
         ' Integration w.r.t. wavelength of solid angle integrals:', qrad, &
         ' W/cm^2'
      if (nzones == 2) write (luncrt, '(a, es12.5, a)') &
         ' Double for full hemisphere:                            ', qrad*2., &
         ' W/cm^2'

 99   return

      end subroutine perform_integrations

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine read_radiance (lun, icol, line, radiance)
!
!     For the indicated line of sight number n = line, read the file named
!     /LINE-n/neqair.out or /LINE-n/radiance.out (depending on icol) and return
!     the radiance value from its last line.  The radiance.out option has been
!     added to enable a check on the extensions for analyzing wavelength effects
!     where icol = 2 and 3 are used to indicate intensity_scanned.out and
!     intensity.out respectively.  These cases write a one-line file after
!     integrating each line of sight w.r.t. wavelength.  Radiance is the first
!     token on that line.
!
!     NEQAIR v15 complication:  Now, the last line starts with Total, not Final,
!     and a different fixed-format read is needed to get the radiance without
!     parsing the last line.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments (could be inherited, but this is clearer):

      integer, intent (in)  :: lun      ! Logical unit to use
      integer, intent (in)  :: icol     ! 0 => neqair.out; 1 => radiance.out
      integer, intent (in)  :: line     ! LOS number across hemisphere quadrants
      real,    intent (out) :: radiance ! Result for this LOS from NEQAIR

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

      if (icol == 0) then
         do  ! Until EOF
            read (lun, '(a)', iostat=ios) buffer
            if (ios /= 0) exit
         end do

!        Avoid parsing, etc., as follows:
!        The buffer should look like this for NEQAIR v14:
!           Final radiance value is 1.23456E-01
!        NEQAIR v15:
!   Total radiance from    855.50 to  39600.00 angstroms =      1.23456 W/cm2-sr

         if (buffer(1:5) == 'Final') then
            read (buffer, '(23x, es12.5)', iostat=ios) radiance
         else  ! buffer(1:5) == 'Total'
            read (buffer, '(54x, f13.5)', iostat=ios) radiance
         end if
      else  ! icol = 1 => radiance is at the start of a one-line file
         read (lun, *, iostat=ios) radiance
      end if
      close (lun)

      if (ios /= 0) then
         write (luncrt, '(a, i6, /, 2a)') &
            ' Trouble reading radiance for line of sight #', line, &
            ' File name: ', filename
      end if

      write (luncrt, '(a, es12.5)') ' Radiance found: ', radiance

   99 return

      end subroutine read_radiance

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
!     Local constants:

      real, parameter :: one = 1.

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine save_integrations ()
!
!     Save the results of wavelength-dependent integrations. For the body-normal
!     line 1 (only), also save cumulative integrals w.r.t. uniform wavelengths.
!     Single line 'radiance.out' files can be processed by another run of this
!     program in its original radiance-integration mode.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: i, l
      real    :: radiance
      character (1), parameter :: method = 'L'  ! Linear quadratures

      nuniform = iuniform(nregion+1)
      if (iactive == 1) then  ! Body-normal LOS (only): cumulative integrals

         allocate (cumulative(nuniform))

         call lcsareas (nuniform, uni_lambda, uni_radiance, method, zero, &
                        cumulative)
         open  (lunrad, file='LINE-1/line-1.cumulative.dat', status='unknown', &
                iostat=ios)
         write (lunrad, '(f9.1, es12.5)') &
            (uni_lambda(i), cumulative(i), i = 1, nuniform)
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

   end program neqair_integration
