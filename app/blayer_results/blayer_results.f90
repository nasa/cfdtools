!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program blayer_results
!
!  Description:
!
!     BLAYER_RESULTS reads a multiblock real gas flow solution and derives a
!  one-layer set of results, some of which apply to the wall and some to the
!  boundary layer edge.  This version reads and writes PLOT3D-type files as
!  originally, but it can also read Tecplot ASCII files and write Tecplot ASCII
!  or Tecplot binary files (now the recommended usage).
!
!     SI units should be input.  SI units are output, but English units may also
!  be requested as an additional output (Tecplot format only).
!
!     Five gas species are assumed.  See also the later BLAYER2D and BLAYER3D,
!  which handle variable numbers of species, zero  or more nonstandard flow
!  quantities, and 1 or 2 temperatures.
!
!     For multidiscipline reasons, a secondary file of OML tile polygon data
!  may also be interpolated at the surface grid points via triangulation, and
!  a related file is searched for the TPS thickness at the relevant centroid.
!  As tile material can change with time and vehicle, just the tile ID is
!  recorded at each surface grid point along with the thickness.
!
!  Assumptions:
!
!     >  The surface was originally required to be at k = 1 for all grid blocks.
!        This is no longer required: the block face with the smallest average
!        initial increment off it is determined, and if necessary the block is
!        permuted to put it in the k = 1 position.
!        This generalization still does not handle cavity blocks well, and it
!        is recommended that they be suppressed, possibly via the optional
!        blayer_results.inp.2 file.
!        BELATED REALIZATION:  PLOT3D-type output makes swapping block faces
!        a challenge, because the header records with block dimensions must
!        be written before any block can be written.  Since storing the entire
!        grid or reading it twice are both highly undesirable for large grids,
!        the generalization is limited to Tecplot-type output (for now).
!
!     >  The radial lines are sufficiently normal to the surface for the
!        1-D boundary layer thickness method to make sense.
!
!     >  The flow quantities are given at the grid-points, not cell centers.
!
!     >  Lee-side results will not be used:  the boundary layer edge detection
!        method used is inappropriate for separated flows.
!
!  Input Flow Field File (PLOT3D function file format, SI units):
!
!     Fields 1-14 are expected to contain grid point values as follows:
!
!        1:   density                      9:  u velocity component
!        2:   pressure                    10:  v    "        "
!        3:   temperature                 11:  w    "        "
!        4:   N2 species density          12:  total enthalpy, H
!        5:   O2    "      "              13:  Mach number
!        6:   NO    "      "              14:  viscosity
!        7:   N     "      "
!        8:   O     "      "
!
!  Equivalent Tecplot ASCII Input Format (1):
!
!        TITLE     = ""
!        VARIABLES = "x, m"
!        "y, m"
!        "z, m"
!        "rho, kg/m^3"
!        "p, Pa"
!        "T, K"
!        "c_N_2"
!        "c_O_2"
!        "c_N_O"
!        "c_N"
!        "c_O"
!        "u, m/s"
!        "v, m/s"
!        "w, m/s"
!        "H0, J/kg"
!        "M"
!        "mu, Pa.s"
!        ZONE T="G1"
!         I=17, J=25, K=81, ZONETYPE=Ordered
!         DATAPACKING=BLOCK
!         DT=(SINGLE SINGLE ... (17 of them) ... SINGLE SINGLE )
!         6.647333145E+00 6.575638294E+00 6.489704609E+00 6.390774727E+00 6. ...
!          :               :               :               :               :
!
!  The Tecplot I/O package used can also read the simpler ASCII format employed
!  by DPLR's postprocessor:
!
!             [Optional title]
!             variables=x,y,z,rho,p,T,C_n2,C_o2,C_no,C_n,C_o,u,v,w,h,M,mu
!             zone t="flow2d" F=point, i= 161 j= 157 k=  81
!             6.64733315E+00 6.57563824E+00 6.48970469E+00 6.39077472E+00 6. ...
!              :              :              :              :              :
!
!  Input TPS Files (enter file name = none to suppress their processing):
!
!     Tile Centroid Thicknesses           OML Tile Polygon Vertices
!
!        190002003 1.121            190002003     1299.006  -97.581  264.372
!        190002004 1.100            190002003     1294.764 -101.823  264.824
!        190002005 1.094            190002003     1299.006 -106.066  265.182
!        190002006 1.063            190002003     1303.249 -101.823  264.722
!        190002019 1.004            190002004     1292.642  -95.459  264.262
!        190002061 1.088            190002004     1288.400  -99.702  264.716
!        190002062 1.042            190002004     1292.642 -103.945  265.058
!        190002065 1.054            190002004     1296.885  -99.702  264.595
!        190002069 1.049            190002005     1288.400  -99.702  264.716
!            :      :                   :             :        :        :
!
!  TPS Data Assumptions:
!
!     >  Tile numbers or IDs are in ascending order.  (Actually, the vertex data
!        just need to have all the vertices for a given tile contiguous.)
!     >  The units of the two files are consistent (SI or English).  The
!        control file should contain the scale factor needed to convert to
!        meters.
!     >  If an ID found by searching the tile triangulation is not found among
!        the thickness data, a warning is issued and the thickness of the
!        nearest preceding tile is used.
!     >  Tile IDs, though integers, are read and written as reals because they
!        are included among the other boundary layer results.
!
!  Output Results (PLOT3D Function File + Grid File, or Tecplot File):
!
!     Units are SI except heat flux is in W/cm^2 and film coefficient is in
!     kg/cm^2.s.
!
!     This version also writes results in English units unless the associated
!     file name is entered as 'none' (Tecplot format only).
!
!     The first 14 quantities (not counting x, y, z) apply to the k = 1 surface.
!     The next 19 quantities apply to the boundary layer edge.  Several further
!     quantities are appended:
!
!            Wall                           Boundary layer edge
!
!            x                         15:  density
!            y                         16:  pressure
!            z                         17:  temperature
!        1:  density                   18:  total enthalpy
!        2:  pressure                  19:  u
!        3:  temperature               20:  v
!        4:  enthalpy                  21:  w
!        5:  viscosity                 22:  Mach number
!        6:  N2 species density        23:  viscosity
!        7:  O2    "       "           24:  N2 species density
!        8:  NO    "       "           25:  O2    "       "
!        9:  N     "       "           26:  NO    "       "
!       10:  O     "       "           27:  N     "       "
!       11:  heat flux                 28:  O     "       "
!       12:  tau_x                     29:  boundary layer height
!       13:  tau_y                     30:  displacement thickness, delta-*
!       14:  tau_z                     31:  momentum thickness, theta
!                                      32:  unit Reynolds #, Re-ue (see Notes)
!                                      33:  film coefficient, CH
!     Miscellaneous quantities:
!
!       34:  Re-kk (see Notes)
!       35:  TPS tile number
!       36:  TPS tile thickness
!
!     Notes:
!
!        1:  Heat flux is not among the data inputs (problematic in the volume).
!            At the wall, derive it from the surface temperature and the single
!            emissivity input.
!        2:  Re-ue is the unit Reynolds # based on edge conditions.
!            To obtain Re-theta or Re-ke, multiply it by theta or k.
!        3:  Re-kk is Reynolds # based on roughness height k & conditions at k,
!            but using viscosity at the wall rather than at k.
!            An Re-kk profile is written to standard output for the surface
!            point indicated in the control file (following the enthalpy ratio
!            profile originally specified this way).
!
!            Reference for Re-kk:
!
!            "Review and Synthesis of Roughness-Dominated Transition
!             Correlations for Reentry Applications" by Daniel C. Reda,
!            Journal of Spacecraft and Rockets, Vol. 39, No. 2, Mar-Apr 2002
!
!  Control File ('blayer_results.inp'):
!
!     TITLE
!     INPUT PLOT3D GRID OR TECPLOT FILE
!     ***.g                 ! Grid file name (or Tecplot input file)
!     T                     ! T = formatted; F = unformatted
!     1                     ! 1 = DATAPACKING=POINT; 2 = DATAPACKING=BLOCK
!     INPUT PLOT3D FLOW FIELD
!     ***.f                 ! Function file name (none for Tecplot input)
!     T                     ! T = formatted; F = unformatted
!     INPUT OML TILE VERTICES and TILE CENTROID THICKNESSES
!     ***.dat  or  none     ! Tile IDs and centroid thicknesses
!     ***.dat               ! Tile IDs and OML polygon vertex (x,y,z)s
!     0.0254                ! Scale factor to convert tile data to meters
!     0.    0.    0.        ! x,y,z shifts to align tile data & grid, meters
!     F                     ! Save the triangulation for plotting?
!     F                     ! Print missing tile thickness diagnostics?
!     OUTPUT PLOT3D GRID FILE OR TECPLOT FILE (SI units)
!     ******.g              ! Derived results grid file (or Tecplot output)
!     T                     ! T = formatted; F = unformatted
!     1                     ! 1 = DATAPACKING=POINT; 2 = DATAPACKING=BLOCK
!     OUTPUT PLOT3D FUNCTION FILE
!     ******.f              ! Derived results file (or none for Tecplot output)
!     T                     ! T = formatted; F = unformatted
!     OUTPUT TECPLOT FILE (English units or suppressed)
!     ******.dat            ! Tecplot outputs, English units (or none)
!     T                     ! T = formatted; F = unformatted
!     2                     ! 1 = DATAPACKING=POINT; 2 = DATAPACKING=BLOCK
!     MISCELLANEOUS CONTROLS
!     0.89                  ! Emissivity (0.89 = RCG everywhere)
!     0.001524              ! Roughness height k (meters); 0.060" for Shuttle
!     2   33  65            ! Block # & (i,j) for sample H/Hinf & Re-kk profiles
!
!  Optional Ancillary Control File ('blayer_results.inp.2'):
!
!     This file can be used to suppress boundary layer processing of specified
!     grid blocks.  If present, it should contain a list of block numbers to
!     suppress, in any obvious convenient form, on a single line.  For example:
!
!     31:36 or 31-36   would both expand to 31, 32, 33, 34, 35, 36
!     30 32 34:40      or any other such intelligible list
!
!     If the file is not present, no blocks are suppressed.
!     At present, "suppressed" means the blocks are not processed, but dummy
!     results are written as output.  This avoids possible confusion resulting
!     from renumbering of blocks in the output.
!
!  Boundary Layer Edge Strategy (Each Radial Line):
!
!     See History and subroutine blayer_edge below.
!
!  Procedures:
!
!     Tecplot_io package  I/O utilities for Tecplot reads and writes
!     XYZQ_IO package     I/O utilities for PLOT3D grid and function files
!
!  History:
!
!     06/19/04  D. Saunders  Initial implementation.
!     07/06/04   "     "     Added TPS data processing.
!     07/20/04   "     "     Switched from two-layer to one-layer outputs.
!     07/27/04   "     "     Option to suppress saving the tile triangulation.
!     07/29/04   "     "     Option to read and write Tecplot files.
!     08/07/04   "     "     Output (H/Hinf, Wall distance) profile at specified
!                            surface point for diagnostic purposes.
!     08/08/04   "     "     The edge detection method now works with a
!                            moderated form of curvature to get more smoothly
!                            varying results across the surface.
!     08/09/04   "     "     Switched from parametric to non-parametric form in
!                            the curvature calculations used for edge detection.
!     08/19/04   "     "     Added calculation of roughness Reynolds numbers
!                            Re-ke and Re-kk.  Height k is a new input.
!     08/27/04   "     "     SEARCH_ADT has a third coefficient argument now.
!                "     "     (PROJECT3 usage is replaced by NEAREST_TRI_POINT.)
!     09/03/04   "     "     Replaced Re-theta and Re-ke with Re-ue, and saved
!                            the Re-kk profile for the same radial line as
!                            specified for the saved H/Hinf profile.  A single
!                            input k for Re-kk does not allow Re-kk for other
!                            values of k to be derived, but the consensus was to
!                            use the Shuttle distributed roughness figure of
!                            0.060" = 0.001524 m.
!     11/20/04   "     "     yline(1) used for surface_xyz(2) was no longer the
!                            surface y.  All the searches for tile ID and thick-
!                            ness were affected.  Ouch!
!     12/09/04   "     "     Added y' and y" to the curvature plot for the
!                            indicated profile, hoping to understand some noisy
!                            edge thickness calculations for a 2-D wedge.
!     12/10/04   "     "     These derivatives seem to have no useful features
!                            in the region of interest.  Therefore, stay with
!                            curvature, but introduce some explicit smoothing
!                            to help with cases observed where the curvature is
!                            essentially constant for several grid points.
!     12/13/04   "     "     Shuttle results suffer if every profile has its
!                            curvature smoothed, even with just one iteration.
!                            Therefore, confine the smoothing to profiles that
!                            have two or more curvatures within a few percent
!                            of each other at the relevant peak.
!     02/03/05   "     "     First step towards dealing with cavity blocks:
!                            transparent provision for blanking of such blocks
!                            via an optional 'blayer_results.inp.2' file.
!                            All output flow quantities for blanked blocks are
!                            set to 1.0.  Blocks are not renumbered.
!     02/11/05   "     "     Bug fix in subroutine blayer_edge, q.v.
!     05/14/05   "     "     Tecplot_io now returns variable names and data-
!                            packing from a read.  Leave input datapacking in
!                            the control file as for output datapacking, but
!                            the input datapacking is actually now immaterial.
!     08/23/05   "     "     Introduced PERMUTE_BLOCK to handle cases where
!                            k = 1 is not at the wall.  Rather than requiring
!                            input boundary condition data, the face with the
!                            smallest average off-face grid spacing is treated
!                            as the wall face.  However, this is not feasible
!                            for PLOT3D-type output except at a high performance
!                            cost, so for now it applies to Tecplot output only,
!                            although adjusted block dimensions are reported for
!                            the PLOT3D case and could (gulp) be edited in.
!     01/17/06   "     "     Retrofit of enhancements initially developed in
!                            BLAYER2D and -3D:
!                            Delta* and theta formulations now use tangential
!                            velocities, not total velocity magnitudes.
!                            Following experiments with curve fitting and
!                            iterating on the arc length scaling to overcome
!                            a known weakness in the edge calculation method
!                            (namely, profile curvatures are not independent of
!                            the arc length scaling), the best work-around so
!                            far is to adjust the curvature-based result by
!                            picking the location of 99.5% OF THE PROFILE PEAK
!                            IN THAT NEIGHBORHOOD.  For clean 2-D data, this
!                            matches the traditional method; for 3-D data, it
!                            handles the fact that the total enthalpy ratios
!                            can differ significantly from 1 in the region of
!                            peak curvature.
!     01/18/06   "     "     Optional additional Tec output in English units.
!     02/03/06   "     "     Momentum thickness had an extra factor of
!                            |vtangent|/|vtangent(edge)| in the integral. Thanks
!                            to Frank Greene for noticing it looked wrong.
!     07/07/06   "     "     The traditional edge method was not trapping non-
!                            monotonic enthalpy ratios properly.
!     07/21/06   "     "     Avoiding inverse interpolation in the traditional
!                            edge method was a bad idea: the 4-point spline
!                            method does not degrade as gracefully to 3 or 2
!                            points as was assumed.  Profiles with overshoots
!                            are best handled by either linear interpolation or
!                            by the nonlinear inverse interpolation now used.
!     12/04/06   "     "     Translated to use the Tecplot 360 version of the
!                            I/O (needed at LaRC).
!     12/05/06   "     "     The output # zone auxiliaries was set after the
!                            allocation instead of before.
!     01/08/08   "     "     Bill Wood noticed that the English units form of
!                            unit Reynolds number needs scaling to be per inch.
!     06/24/08   "     "     Discontinuities towards the Shuttle wing tip were
!                            traced to profiles that straighten up short of 1.0
!                            before achieving 1.0.  This means the heuristic
!                            size of the neighborhood of the peak curvature
!                            can include ~1.0 for one profile but not for a
!                            neighboring profile.  Stage 2 of the edge method
!                            then seeks 99.5% of quite different peaks in those
!                            neighborhoods.  Use of 99.5% means even small
!                            changes lead to large differences in edge thickness
!                            when the profile is so steep.  Mike Olsen suggested
!                            using 95% to reduce the effect greatly, but then
!                            all edge-related quantities would be significantly
!                            lower everywhere.  After much pondering, we stay
!                            with 99.5%, and accept that wing tip regions are of
!                            limited interest anyway, even on the wind side.
!
!  Author:  David Saunders, ELORET/NASA Ames Research Center, Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_header_structure ! See Tecplot_io.f90 for the first 3 of these
   use grid_block_structure
   use tecplot_io_module
   use xyzq_io_module
   use adt_utilities

   implicit none

!  Constants:

   integer, parameter :: &
      lunctl     = 1,    &
      luntps     = 2,    &
      lunkbd     = 5,    &
      luncrt     = 6,    &
      lunxyz_in  = 7,    &
      lunq_in    = 8,    &
      lunxyz_out = 9,    &
      lunq_out   = 10,   &
      lun_eu     = 11,   &   ! Optional English units Tecplot output
      name_limit = 32,   &   ! Must match the value in Tecplot_io_module
      ndim       = 3,    &
      numf_in    = 14,   &   ! Expected number of input functions per grid pt.
      numf_out   = 36        ! Number of fields in the derived results file,
                             ! not counting the wall x, y, z
   real, parameter ::    &
      stefan     = 5.66097E-12,  & ! Stefan-Boltzmann constant, w/(cm^2 K^4)
      one        = 1.0,  &
      zero       = 0.0

   character, parameter :: &
      blank  * 1  = ' ',   &
      format * 11 = 'unformatted', &
      none   * 4  = 'none'

   logical, parameter :: &
      false = .false.,   &
      new   = .true.,    &
      true  = .true.

!  Variables:

   integer :: &
      datapacking_eu, datapacking_in, datapacking_out, i, i1, ib, ios,         &
      iprofile, j, jprofile, mi, mj, mk, n, nblank, nblocks, ni, nj, nk,       &
      nnode, nprofile, npts, nthick, ntile, ntri, numq_in, numq_out

   integer, dimension (1) :: &
      iwall                  ! The intrinsic minloc needs an array for output

   integer, allocatable, dimension (:) :: &
      iblock

   integer, pointer, dimension (:,:) :: &
      conn                   ! Pointer, not allocatable, because we can't pass
                             ! unallocated arrays to triangulate_tiles
   real :: &
      emissivity, esigma, roughness_height, units_scale, xshift, yshift, zshift

   real, dimension (6) :: &
      average_ds

   real, allocatable, dimension (:) :: &
      tile_ids, tile_thicknesses

   real, pointer, dimension (:) :: &
      tri_ids       ! As for conn(:,:) above

   real, pointer, dimension (:,:) :: &
      tri_xyzs      ! As for conn(:,:) above

   logical :: &
      cell_centered, English, formatted_eu, formatted_in, formatted_out,       &
      save_triangulation, Tecplot_input, Tecplot_output, tile_diagnostics,     &
      tps_data

   logical, allocatable, dimension (:) :: &
      suppress

   character :: &
      filename_in * 80, filename2 * 80, filename_out * 80, filename_tps * 80,  &
      filename_eu * 80, title_in * 80, title_out * 80

   character, dimension (ndim + numf_out) :: &
      names_eu * (name_limit), names_out * (name_limit)

   character, pointer, dimension (:) :: &
      names_in * (name_limit)

!  Composite data types:

   type (grid_header) :: &
      header_in, header_out, header_eu

   type (grid_type), pointer, dimension (:) :: &
      xyzq_in, xyzq_out

!  Tecplot functions:

   integer :: TecEnd110, TecFil110

!  Explicit interface to allow array allocation in a called procedure:

   interface
      subroutine triangulate_tiles (lunrd, save_triangulation,          &
                                    ntile, ntri, tri_ids, tri_xyzs, conn)
      implicit none
      integer, intent (in)  :: lunrd
      logical, intent (in)  :: save_triangulation
      integer, intent (out) :: ntile, ntri
      real,    pointer      :: tri_ids(:), tri_xyzs(:,:)
      integer, pointer      :: conn(:,:)

      end subroutine triangulate_tiles
   end interface

!  Data:

   data names_eu          &  ! Avoid commas used as delimiters for TecIni
     /'xw (in)',          'yw (in)',          'zw (in)',          &
      'rhow (lbm/ft3)',   'pw (psf)',         'Tw (F)',           &
      'Hw (BTU/lbm)',     'muw (psf.s)',      'cN2w',             &
      'cO2w',             'cNOw',             'cNw',              &
      'cOw',              'qw (BTU/ft2.s)',   'tauwx (psf)',      &
      'tauwy (psf)',      'tauwz (psf)',      'rhoe (lbm/ft3)',   &
      'pe (psf)',         'Te (F)',           'He (BTU/lbm)',     &
      'ue (ft/s)',        've (ft/s)',        'we (ft/s)',        &
      'Me',               'mue (psf.s)',      'cN2e',             &
      'cO2e',             'cNOe',             'cNe',              &
      'cOe',              'delta (in)',       'deltastar (in)',   &
      'theta (in)',       'Re-ue',            'CH (lbm/ft2.s)',   &
      'Re-kk',            'tile number',      'tile thickness (in)'/

   data names_out         &
     /'xw (m)',           'yw (m)',           'zw (m)',           &
      'rhow (kg/m3)',     'pw (Pa)',          'Tw (K)',           &
      'Hw (J/kg)',        'muw (Pa.s)',       'cN2w',             &
      'cO2w',             'cNOw',             'cNw',              &
      'cOw',              'qw (W/cm2)',       'tauwx (Pa)',       &
      'tauwy (Pa)',       'tauwz (Pa)',       'rhoe (kg/m3)',     &
      'pe (Pa)',          'Te (K)',           'He (J/kg)',        &
      'ue (m/s)',         've (m/s)',         'we (m/s)',         &
      'Me',               'mue (Pa.s)',       'cN2e',             &
      'cO2e',             'cNOe',             'cNe',              &
      'cOe',              'delta (m)',        'deltastar (m)',    &
      'theta (m)',        'Re-ue',            'CH (kg/cm2.s)',    &
      'Re-kk',            'tile number',      'tile thickness (m)'/

!  Execution:
!  !!!!!!!!!!

   open (lunctl, file='blayer_results.inp', status='old', iostat=ios)

   if (ios /= 0) then
      write (luncrt, '(/, a)') &
         ' Unable to open blayer_results.inp control file.'
      go to 99
   end if

!  Echo it to the output log, using filename2 as a buffer:

   do
      read (lunctl, '(a)', iostat=ios) filename2
      if (ios < 0) exit

      write (luncrt, '(1x, a)') trim (filename2)
   end do

   rewind (lunctl)

   read (lunctl, *)                  ! Case description
   read (lunctl, *)                  ! Header
   read (lunctl, *) filename_in      ! Input grid or Tecplot file
   header_in%filename = filename_in
   read (lunctl, *) formatted_in
   header_in%formatted = formatted_in
   i1 = 1; if (formatted_in) i1 = 3
   read (lunctl, *) datapacking_in   ! 1 or 2 for POINT or BLOCK (Tecplot)
!  header_in%datapacking is an output from the read, not a control input
   header_in%ndim = ndim

   read (lunctl, *)                  ! Header
   read (lunctl, *) filename2        ! Input function file, or none if Tecplot
   Tecplot_input =  filename2(1:4) == none
   read (lunctl, *) !!! formatted_in ! Must be the same for XYZQ_IO usage

!  The XYZQ_IO package is more modular than the Tecplot I/O package, which
!  both opens and reads or writes in one call.

   if (.not. Tecplot_input) then ! Open the input PLOT3D files

      open (lunxyz_in, file=filename_in, status='old', form=format(i1:11), &
            iostat=ios)

      if (ios /= 0) then
         write (luncrt, '(/, 2a)') &
            ' Cannot open the input grid file: ', trim (filename_in)
         go to 99
      end if

      open (lunq_in, file=filename2, status='old', form=format(i1:11), &
            iostat=ios)

      if (ios /= 0) then
         write (luncrt, '(/, 2a)') &
            ' Cannot open the input flow file: ', trim (filename2)
         go to 99
      end if

   end if

!  Deal with the optional TPS data files:
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   read (lunctl, *)
   read (lunctl, *) filename_tps     ! Tile centroid thickness data file
   tps_data = filename_tps(1:4) /= none

   if (tps_data) then

      open (luntps, file=filename_tps, status='old', iostat=ios)
      if (ios /= 0) then
         write (luncrt, '(/, 2a)') &
            ' Cannot open the TPS thickness file: ', trim (filename_tps)
         go to 99
      end if

!     Count the number of lines of tile thickness data:

      nthick = 0
      do ! Until EOF
         read (luntps, *, iostat=ios)
         if (ios < 0) exit

         nthick = nthick + 1
      end do

      rewind (luntps)

      write (luncrt, '(/, a, i7, /)') ' # tile thicknesses found:', nthick

      allocate (tile_ids(nthick), tile_thicknesses(nthick))

      read (luntps, *, iostat=ios) &
         (tile_ids(n), tile_thicknesses(n), n = 1, nthick)

      if (ios /= 0) then
         write (luncrt, '(/, 2a, i8)') &
            ' Trouble reading tile thicknesses.  # tiles: ', nthick
         go to 99
      end if

      close (luntps)

      read (lunctl, *) filename_tps
      read (lunctl, *) units_scale            ! Tile units to meters
      read (lunctl, *) xshift, yshift, zshift ! Meters
      read (lunctl, *) save_triangulation     ! T|F
      read (lunctl, *) tile_diagnostics       ! T|F

      open (luntps, file=filename_tps, status='old', iostat=ios)

      if (ios /= 0) then
         write (luncrt, '(/, 2a)') ' Unable to open OML polygon file: ', &
            trim (filename_tps)
         go to 99
      end if

!     Read the OML tile polygon vertices and triangulate them:

      call triangulate_tiles (luntps, save_triangulation,                &
                              ntile, ntri, tri_ids, tri_xyzs, conn)
!     ----------------------

      nnode = ntri + ntile ! # (x,y,z)s counting the tile centroids

      do i = 1, nnode
         tri_xyzs(1,i) = units_scale * tri_xyzs(1,i) + xshift
         tri_xyzs(2,i) = units_scale * tri_xyzs(2,i) + yshift
         tri_xyzs(3,i) = units_scale * tri_xyzs(3,i) + zshift
      end do

      tile_thicknesses(:) = units_scale * tile_thicknesses(:)

!     Build a search tree for the triangulated tiles:

      call build_adt (nnode, ntri, conn, tri_xyzs)
!     --------------

   else ! TPS data processing has been suppressed

      read (lunctl, *) ! Skip the OML data file name
      read (lunctl, *) ! units_scale
      read (lunctl, *) ! xshift, yshift, zshift
      read (lunctl, *) ! save_triangulation
      read (lunctl, *) ! tile_diagnostics
      save_triangulation = false
      tile_diagnostics   = false

   end if

!  Deal with the output files:
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!

   read (lunctl, *)
   read (lunctl, *) filename_out     ! Output grid or Tecplot file
   header_out%filename = filename_out
   read (lunctl, *) formatted_out
   header_out%formatted = formatted_out
   i1 = 1; if (formatted_out) i1 = 3
   read (lunctl, *) datapacking_out  ! 1 or 2 for POINT or BLOCK (Tecplot)
   header_out%datapacking = datapacking_out - 1 ! Tecplot 360 usage
   header_out%ndim = ndim

   read (lunctl, *)
   read (lunctl, *) filename2        ! Output function file, or none for Tecplot
   Tecplot_output = filename2(1:4) == none
   read (lunctl, *) !!! formatted_out

   if (Tecplot_output) then
!     The initialization includes opening the file
   else

      open (lunxyz_out, file=filename_out, status='unknown', &
            form=format(i1:11), iostat=ios)

      if (ios /= 0) then
         write (luncrt, '(/, 2a)') &
            ' Cannot open the output grid file: ', trim (filename_out)
         go to 99
      end if

      open (lunq_out, file=filename2, status='unknown', form=format(i1:11), &
            iostat=ios)

      if (ios /= 0) then
         write (luncrt, '(/, 2a)') &
            ' Cannot open the output function file: ', trim (filename2)
         go to 99
      end if

   end if

   read (lunctl, *)
   read (lunctl, *) filename_eu      ! English units Tecplot output, or none
   header_eu%filename = filename_eu
   English = filename_eu(1:4) /= none
   read (lunctl, *) formatted_eu     ! Tec output initialization opens the file
   header_eu%formatted = formatted_eu
   read (lunctl, *) datapacking_eu   ! 1 or 2 for POINT or BLOCK (Tecplot)
   header_eu%datapacking = datapacking_eu - 1 ! Tecplot 360 usage
   header_eu%ndim = ndim

!  Miscellaneous controls:
!  !!!!!!!!!!!!!!!!!!!!!!!

   read (lunctl, *, iostat=ios)
   if (ios /= 0) then
      write (luncrt, '(/, a)') ' Miscellaneous controls missing?'
      go to 99
   end if
   read (lunctl, *) emissivity
   read (lunctl, *) roughness_height
   read (lunctl, *, iostat=ios) nprofile, iprofile, jprofile ! (n,i,j) to save

   if (ios /= 0) then
      nprofile = 1;  iprofile = 1;  jprofile = 1 ! Default to point 1 of block 1
   end if

   close (lunctl)

!  Read the input file headers:
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!

   if (Tecplot_input) then

      ios = 0  ! Non-verbose mode

      call Tec_header_read (lunxyz_in, header_in, xyzq_in, ios)

      if (ios /= 0) then
         write (luncrt, '(/, a)') ' Trouble reading the Tecplot file header.'
         go to 99
      end if

      numq_in = header_in%numq
      nblocks = header_in%nblocks

   else ! PLOT3D input files

      call xyz_header_io (1, lunxyz_in, formatted_in, nblocks, xyzq_in, ios)

      if (ios /= 0) then
         write (luncrt, '(/, a)') ' Trouble reading the grid file header.'
         go to 99
      end if

      call q_header_io (1, lunq_in, formatted_in, nblocks, numq_in, xyzq_in,   &
                        ios)
      if (ios /= 0) then
         write (luncrt, '(/, a)') ' Trouble reading the function file header.'
         go to 99
      end if

   end if

   if (numq_in /= numf_in) then
      write (luncrt, '(/, (a, i3))')                              &
         ' Unexpected number of input flow quantities:', numq_in, &
         ' Expected number: ', numf_in
      go to 99
   end if

!  Display the block dimensions to reassure the user:

   write (luncrt, '(/, a, /)') ' Input grid dimensions:'
   write (luncrt, '(4(i6, 2X, 3i5))') &
      (ib, xyzq_in(ib)%ni, xyzq_in(ib)%nj, xyzq_in(ib)%nk, ib = 1, nblocks)

!  Look for the optional control file allowing suppression of blocks:

   allocate (suppress(nblocks), iblock(nblocks))

   suppress(:) = false

   open (lunctl, file='blayer_results.inp.2', status='old', iostat=ios)

   if (ios == 0) then
      nblank = nblocks

      write (luncrt, '(a)')

      call rdlist (luncrt, &
         ' Reading blayer_results.inp.2 for block numbers to suppress. ', &
         lunctl, nblank, iblock)

      write (luncrt, '(a)')

      do i = 1, nblank
         ib = iblock(i)
         if (ib < 1 .or. ib > nblocks) then
            write (luncrt, '(a, i9)') &
               ' Bad block number to suppress:', ib
         else
            suppress(ib) = true
         end if
      end do

      close (lunctl)
   end if

   deallocate (iblock)

!  Set up the output files:
!  !!!!!!!!!!!!!!!!!!!!!!!!

   allocate (xyzq_out(nblocks), stat=ios)

   if (ios /= 0) then
      write (luncrt, '(/, a)') ' Trouble allocating the output blocks.'
      go to 99
   end if

   do ib = 1, nblocks
      xyzq_out(ib)%ni = xyzq_in(ib)%ni
      xyzq_out(ib)%nj = xyzq_in(ib)%nj
      xyzq_out(ib)%nk = 1
      xyzq_out(ib)%mi = xyzq_in(ib)%mi
      xyzq_out(ib)%mj = xyzq_in(ib)%mj
      xyzq_out(ib)%mk = 1
   end do

   if (tps_data) then
      numq_out = numf_out
   else
      numq_out = numf_out - 2
   end if

   header_out%numq = numq_out
   header_eu%numq  = numq_out

   if (Tecplot_output .or. English) then

      if (.not. Tecplot_input) then
         title_out = trim (filename_out)
      else if (len_trim (header_in%title) > 0) then
         title_out = trim (header_in%title)
      else
         title_out = trim (filename_in)
      end if

      header_out%title = title_out;  header_out%ndatasetaux = 0
      header_eu%title  = title_out;  header_eu%ndatasetaux  = 0
   end if

   if (Tecplot_output) then

      header_out%nblocks = nblocks

      allocate (header_out%varname(ndim + numq_out))

      header_out%varname(:) = names_out(1:ndim + numq_out)

      call Tec_header_write (lunxyz_out, header_out, ios)

      if (ios /= 0) then
         write (luncrt, '(/, a)') ' Trouble writing the SI Tecplot file header.'
         go to 99
      end if

   else ! PLOT3D-type output

      call xyz_header_io (2, lunxyz_out, formatted_out, nblocks, xyzq_out, ios)

      if (ios /= 0) then
         write (luncrt, '(/, a)') ' Trouble writing the grid file header.'
         go to 99
      end if

      call q_header_io (2, lunq_out, formatted_out, nblocks, numq_out,         &
                        xyzq_out, ios)
      if (ios /= 0) then
         write (luncrt, '(/, a)') ' Trouble writing the function file header.'
         go to 99
      end if

   end if

   if (English) then  ! Optional extra output in English units

      header_eu%nblocks = nblocks

      allocate (header_eu%varname(ndim + numq_out))

      header_eu%varname(:) = names_eu(1:ndim + numq_out)

      call Tec_header_write (lun_eu, header_eu, ios)

      if (ios /= 0) then
         write (luncrt,'(/, a)') ' Trouble initializing the English units file.'
         go to 99
      end if

      if (.not. formatted_eu) then
         if (.not. formatted_out) ios = TecFil110 (1) ! Switch back to SI file
      end if

   end if

   esigma = emissivity * stefan


!  Process one block at a time:
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!

   do ib = 1, nblocks

      ni = xyzq_in(ib)%ni
      nj = xyzq_in(ib)%nj
      nk = xyzq_in(ib)%nk

!     Read an input block:

      if (Tecplot_input) then

         call Tec_block_allocate (xyzq_in(ib), ndim, numq_in, ios)

         if (ios /= 0) then
            write (luncrt, '(/, a, i5)') &
               ' Allocation trouble for Tecplot file input.  Block #:', ib
            go to 99
         end if

         call Tec_block_read (lunxyz_in, header_in, xyzq_in(ib), ios)

         if (ios /= 0) then
            write (luncrt, '(/, a, i5)') &
               ' Trouble reading Tecplot file.  Block #:', ib
            go to 99
         end if

      else ! PLOT3D-type inputs

         call xyz_allocate (xyzq_in(ib), ios)

         if (ios /= 0) then
            write (luncrt, '(/, a, i5)') &
               ' Allocation trouble for PLOT3D grid input.  Block #:', ib
            go to 99
         end if

         npts = ni * nj * nk

         call xyz_block_io (1, lunxyz_in, formatted_in, npts,               &
                            xyzq_in(ib)%x, xyzq_in(ib)%y, xyzq_in(ib)%z, ios)
         if (ios /= 0) then
            write (luncrt, '(/, a, i5)') &
               ' Trouble reading PLOT3D grid.  Block #:', ib
            go to 99
         end if

         call q_allocate (xyzq_in(ib), numq_in, ios)

         if (ios /= 0) then
            write (luncrt, '(/, a, i5)') &
               ' Allocation trouble for PLOT3D function file.  Block #:', ib
            go to 99
         end if

         call q_block_io (1, lunq_in, formatted_in, numq_in,                &
                          xyzq_in(ib)%mi, xyzq_in(ib)%mj, xyzq_in(ib)%mk,   &
                          xyzq_in(ib)%q, ios)

         if (ios /= 0) then
            write (luncrt, '(/, a, i5)') &
               ' Trouble reading PLOT3D function file.  Block #:', ib
            go to 99
         end if

      end if

!     Identify the (likely) wall from the average off-face grid spacings:

      call average_increments (ni, nj, nk, xyzq_in(ib)%x, xyzq_in(ib)%y,       &
                               xyzq_in(ib)%z, average_ds)

      iwall = minloc (average_ds)

      if (iwall(1) /= 5) then ! Swap faces in-place:

         mi = ni;  mj = nj;  mk = nk

         call permute_block (1, numq_in, mi, mj, mk, iwall(1), xyzq_in(ib)%x,  &
                             xyzq_in(ib)%y, xyzq_in(ib)%z, xyzq_in(ib)%q,      &
                             ni, nj, nk)

         xyzq_in(ib)%ni  = ni;  xyzq_out(ib)%ni = ni;  xyzq_out(ib)%mi = ni
         xyzq_in(ib)%nj  = nj;  xyzq_out(ib)%nj = nj;  xyzq_out(ib)%mj = nj
         xyzq_in(ib)%nk  = nk

!!!      call Tecplot_write (60, 'permuted.block', true, 1, title_in,          &
!!!                          nblocks, numq_in, names_in, xyzq_in, ios)
!!!      if (ios /= 0) then
!!!         write (6, *) 'Trouble writing permuted block. ios:', ios
!!!         go to 99
!!!      end if

         if (.not. Tecplot_output) then ! Weak solution for now

            write (luncrt, '(3i4, a, i5, a)') ni, nj, 1, '  block', ib,        &
               '  dimensions need to be edited into the output file. **********'
         end if

         if (ib == nprofile) then
!           It's not obvious what to do about iprofile and jprofile.
!           The user needs to enter them as they apply to the OUTPUT block.
         end if

      end if

!     Allocate the corresponding output block:

      xyzq_out(ib)%nzoneaux = 0

      call Tec_block_allocate (xyzq_out(ib), ndim, numq_out, ios)

      if (ios /= 0) then
         write (luncrt, '(/, a, i5)') &
            ' Allocation trouble for output.  Block #:', ib
         go to 99
      end if

!     Process all the radial lines of this block, one at a time:
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     The input block may have been rearranged in-place, but the arrays were
!     not deallocated and reallocated.  Therefore, push things down a level.
!     With so many controls, it is most convenient to use a local procedure
!     rather than a bona fide subroutine.

      call process_block (ni, nj, nk, numq_in, xyzq_in(ib)%x, xyzq_in(ib)%y,   &
                          xyzq_in(ib)%z, xyzq_in(ib)%q, numq_out,              &
                          xyzq_out(ib)%x, xyzq_out(ib)%y, xyzq_out(ib)%z,      &
                          xyzq_out(ib)%q, ios)

!!!   if (ios /= 0) go to 99  ! Trouble, but keep going

      call deallocate_blocks (ib, ib, ndim, numq_in, xyzq_in, ios)

      if (ios /= 0) go to 99

!     Save results for this block:
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (Tecplot_output .or. English) then
         if (Tecplot_input) then
            xyzq_out(ib)%zone_title = xyzq_in(ib)%zone_title
         else
            xyzq_out(ib)%zone_title = blank
            write (xyzq_out(ib)%zone_title(1:8), '(a4, i4)') 'Zone', ib
         end if
      end if

      xyzq_out(ib)%solutiontime = -999.  ! Undefined

      if (Tecplot_output) then

         call Tec_block_write (lunxyz_out, header_out, xyzq_out(ib), ios)

         if (ios /= 0) then
            write (luncrt, '(/, a, i5)') &
               ' Trouble writing Tecplot output.  Block #:', ib
            go to 99
         end if

      else ! PLOT3D-type output

         npts = ni * nj

         call xyz_block_io (2, lunxyz_out, formatted_out, npts,                &
                            xyzq_out(ib)%x, xyzq_out(ib)%y, xyzq_out(ib)%z, ios)
         if (ios /= 0) then
            write (luncrt, '(/, a, i5)') &
               ' Trouble writing PLOT3D grid output.  Block #:', ib
            go to 99
         end if

         call q_block_io (2, lunq_out, formatted_out, numq_out,                &
                          xyzq_out(ib)%mi, xyzq_out(ib)%mj, xyzq_out(ib)%mk,   &
                          xyzq_out(ib)%q, ios)
         if (ios /= 0) then
            write (luncrt, '(/, a, i5)') &
               ' Trouble writing PLOT3D function output.  Block #:', ib
            go to 99
         end if

      end if

      if (English) then

         call SI_to_English (ni, nj, numq_out, tps_data, xyzq_out(ib)%x,       &
                             xyzq_out(ib)%y, xyzq_out(ib)%z, xyzq_out(ib)%q)

         if (.not. formatted_eu) then
            if (.not. formatted_out) ios = TecFil110 (2) ! Switch binary context
         end if

         call Tec_block_write (lun_eu, header_eu, xyzq_out(ib), ios)

         if (ios /= 0) then
            write (luncrt, '(/, a, i5)') &
               ' Trouble writing English units output.  Block #:', ib
            go to 99
         end if

         if (.not. formatted_eu) then
            if (.not. formatted_out) ios = TecFil110 (1) ! Revert to SI file
         end if

      end if

      call deallocate_blocks (ib, ib, ndim, numq_out, xyzq_out, ios)

      if (ios /= 0) then
         write (luncrt, '(a, i5)') &
            ' Trouble deallocating x,y,z,q output arrays.  Block #:', ib
         go to 99
      end if

   end do ! Next block to process

   if (Tecplot_output) then
      if (formatted_out) then
         close (lunxyz_out)
      else
         ios = TecEnd110 ()
      end if
   else
      close (lunxyz_out)
      close (lunq_out)
   end if

   if (English) then
      if (formatted_eu) then
         close (lun_eu)
      else
         ios = TecFil110 (2)
         ios = TecEnd110 ()
      end if
   end if

99 continue

! *** stop ! Avoid system dependencies.

!  Internal procedures for program blayer_results:

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine process_block (ni, nj, nk, nf, x, y, z, f, nfout,             &
                                xout, yout, zout, fout, ier)
!
!     Process all radial lines of the given block, which may have been shuffled
!     in-place to put k = 1 at the wall.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in) :: &
         ni, nj, nk,          &       ! Active dimensions of (re)packed blocks
         nf                           ! and the number of associated functions

      real, intent (in), dimension (ni,nj,nk) :: &
         x, y, z                      ! Making these arguments gets around the
                                      ! fact that the calling program's
                                      ! allocated dimensions are still in effect

      real, intent (in), dimension (nf,ni,nj,nk) :: &
         f

      integer, intent (in) :: &
         nfout                        ! # output functions, excluding x, y, z

      real, intent (out), dimension (ni,nj) :: &
         xout, yout, zout             ! Making these arguments cleans the code

      real, intent (out), dimension (nfout,ni,nj) :: &
         fout                         ! Output functions for this block

      integer, intent (out) :: &
         ier                          ! 0 => no problem detected

!     Local constants:

      real,      parameter :: &
         v_epsilon  = 1.e-10          ! Safeguards delta* & theta at a stag. pt.

      character, parameter :: &
         method * 1 = 'B'     ! "Bessel" option: plain 4-pt. local spline method

!     Local variables:

      integer :: &
         i, id, itri, j, k, k1, kedge, l, left, min, mout, nbad, nk_height_k,  &
         numk

      integer, dimension (5) :: &
         k_indices              ! Initialized in a data statement below

      real :: &
         delstar, dsq, emissivity, Hedge, p, q, r, rho_vtan_edge, redge,       &
         tedge, theta, total, Twall, V, vedge, vtangent, vtangent_edge,        &
         wall_mu, yeval, ypeval

      real, save :: &
         Hinf       ! Use free stream value from the first block for all blocks

      real, dimension (3) :: &
         interp_xyz, surface_xyz, unit_normal, v_tangent, v_total

      real, dimension (4) :: &
         height_k_values

      real, allocatable, dimension (:) :: &
         Hratio, rratio, vratio, fline, tline, xline, yline, zline

      logical :: &
         save_profile

      character :: &
         subtitle * 32

!     Data:

      data k_indices       &  ! Items to be interpolated at roughness height "k"
        /1, 9, 10, 11, 14/    ! (rho, u, v, w), along with wall viscosity

!     Execution:

      left = 1                   ! Initialize pointer for tile ID search
      nk_height_k = (2 * nk) / 3 ! Enough to contain roughness height k

      allocate (xline(nk), yline(nk), zline(nk), tline(nk), fline(nk),         &
                Hratio(nk), rratio(nk), vratio(nk))

      if (ib == 1) Hinf = one / f(12,1,1,nk) ! 1 / Free-stream total enthalpy

      nbad = 0  ! Suppress output of problematic profiles if too many

      do j = 1, nj

         do i = 1, ni

            if (suppress(ib)) then ! Fake it

               xout(i,j) = x(i,j,1)
               yout(i,j) = y(i,j,1)
               zout(i,j) = z(i,j,1)

               fout(:,i,j) = one
               cycle
            end if

            save_profile = ib == nprofile .and. &
                           i  == iprofile .and. j == jprofile

            if (save_profile) then
               write (subtitle, '(a6, i4, 2(5x, a2, i4))')                  &
                  'Block:', nprofile, 'i:', iprofile, 'j:', jprofile
            end if

!           Wall values:

            fout(1: 3,i,j) = f(1:3,i,j,1)             ! rho, p, T
            Twall          = f(  3,i,j,1)
            fout(   4,i,j) = f( 12,i,j,1)             ! enthalpy
            fout(   5,i,j) = f( 14,i,j,1)             ! viscosity
            fout(6:10,i,j) = f(4:8,i,j,1)             ! Species
            fout(  11,i,j) = esigma * Twall**4        ! Heat flux

!           Set up this radial line for estimating the boundary layer edge.
!           We need arc length data for the shear stress terms also.

            do k = 1, nk
               xline(k)  = x(i,j,k)
               yline(k)  = y(i,j,k)
               zline(k)  = z(i,j,k)
               Hratio(k) = f(12,i,j,k) * Hinf
            end do

!           Unnormalized arc lengths:

            call chords3d (nk, xline, yline, zline, false, total, tline)

!           Shear stress terms = normal derivatives of tangential velocity:

            unit_normal(1)  = (xline(2) - xline(1)) / tline(2)
            unit_normal(2)  = (yline(2) - yline(1)) / tline(2)
            unit_normal(3)  = (zline(2) - zline(1)) / tline(2)

            v_total(:)      = f(9:11,i,j,2)
            v_tangent(:)    = v_total(:) - dot_product (v_total, unit_normal)  &
                            * unit_normal(:)
            fout(12:14,i,j) = v_tangent(:) * &                     ! (vt2 - 0) *
                              (f(14,i,j,1) / tline(2))             ! (mu / dt)

!           Estimate the edge of the boundary layer, tedge:
!           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            call blayer_edge (nk, Hratio, tline, save_profile, subtitle,       &
                              kedge, Hedge, tedge, ier)

            if (ier /= 0) then
               nbad = nbad + 1
               if (nbad == 1) then
                  write (luncrt, '(/, a, 3i5)') &
                     ' Trouble detecting boundary layer edge. ib, i, j:', ib,i,j
                  write (luncrt, '(/, 2a)') &
                  '!                 x                  y                  z', &
                  '                  t     Enthalpy Ratio'
                  write (luncrt, '(1p, 5e19.11)') &
                  (xline(k), yline(k), zline(k), tline(k), Hratio(k), k = 1, nk)
               end if
            end if ! Keep going anyway


!           Interpolate the grid at the boundary layer edge.  [Not used now.]

!           The local spline method used requires just 4 points:

            k1 = kedge - 1 ! tedge is in the interval [t(kedge), t(kedge+1))

!!!         call lcsfit (4, tline(k1), xline(k1), new, method, 1, &
!!!                      tedge, xedge, ypeval)  ! etc. for y & z

            xout(i,j) = xline(1)
            yout(i,j) = yline(1)
            zout(i,j) = zline(1)


!           Interpolate output quantities 15:28 at the boundary layer edge.

            do mout = 15, 28

               select case (mout)        ! The order changed awkwardly

               case (15:17)              ! rho, p, T
                  min = mout - 14
               case (18)                 ! Total enthalpy
                  min = 12
               case (19:21)              ! u, v, w
                  min = mout - 10
               case (22:23)              ! Mach, mu
                  min = mout - 9
               case (24:28)              ! Species densities
                  min = mout - 20
               end select

               do k = k1, k1 + 3
                  yline(k) = f(min,i,j,k)
               end do

               call lcsfit (4, tline(k1), yline(k1), new, method, 1, tedge,    &
                            yeval, ypeval)
               fout(mout,i,j) = yeval

            end do

            fout(29,i,j) = tedge                ! Edge distance from wall, delta
            redge        =       fout(15,i,j)                   ! Edge density
            vedge        = sqrt (fout(19,i,j)**2 + &            ! Edge total |V|
                                 fout(20,i,j)**2 + &
                                 fout(21,i,j)**2)

            v_total(:)    = fout(19:21,i,j)
            v_tangent(:)  = v_total(:) - dot_product (v_total, unit_normal) *  &
                            unit_normal(:)
            vtangent_edge = sqrt (v_tangent(1)**2 + &           ! Edge |Vtangnt|
                                  v_tangent(2)**2 + v_tangent(3)**2)

!           Set up the ordinates to integrate for displacement thickness:
!           Safeguarding a stag. point requires forming (rho*v) / (rhoe*ve),
!           not (rho/rhoe) * (v/ve) as done originally:

            rho_vtan_edge = redge * vtangent_edge + v_epsilon
            vtangent_edge =         vtangent_edge + v_epsilon

            numk = kedge + 1

            do k = 1, numk
               v_total(:)   = f(9:11,i,j,k)
               v_tangent(:) = v_total(:) - dot_product (v_total, unit_normal)* &
                              unit_normal(:)
               vtangent     = sqrt (v_tangent(1)**2 + &
                                    v_tangent(2)**2 + v_tangent(3)**2)
               rratio(k)    = (f(1,i,j,k)*vtangent + v_epsilon) / &  ! (rho*v) /
                              rho_vtan_edge                     ! (edge rho*v)
               vratio(k)    = (vtangent + v_epsilon) / vtangent_edge
               fline(k)     = one - rratio(k)      ! 1 - (rho*v) / (edge rho*v)
            end do

            call lcsquad (numk, tline, fline, zero, tedge, method, delstar)

            fout(30,i,j) = delstar

            if (save_profile) then
               write (luncrt, '(a, f10.6, a)') &
                  'Profile Integrated for Delta-* = ', delstar, '  Meters'
               write (luncrt, '(a)') subtitle, &
                  '1 - (rho vtan) / (edge rho vtan)', 'Wall distance, meters', &
                  ' $options grid=1, symbol=3, $end'
               write (luncrt, '(1p, 2e19.11)') (fline(k), tline(k), k = 1, numk)
               write (luncrt, '(/, a, /)') ' $options line=''dash'', $end'
               write (luncrt, '(1p, 2e19.11)') zero, tedge, one, tedge
               write (luncrt, '(/, a, /)') 'END FRAME'
            end if

!           Likewise for momentum thickness:

            do k = 1, numk
               fline(k) = rratio(k) * (one - vratio(k))
            end do

            call lcsquad (numk, tline, fline, zero, tedge, method, theta)

            fout(31,i,j) = theta

            if (save_profile) then
               write (luncrt, '(a, f10.6, a)') &
                  'Profile Integrated for Theta = ', theta, '  Meters'
               write (luncrt, '(a)') subtitle, &
                  '(rho vtan / edge rho vtan) * (1 - (vtan / edge vtan))',     &
                  'Wall distance, meters', ' $options grid=1, symbol=3, $end'
               write (luncrt, '(1p, 2e19.11)') (fline(k), tline(k), k = 1, numk)
               write (luncrt, '(/, a, /)') ' $options line=''dash'', $end'
               write (luncrt, '(1p, 2e19.11)') zero, tedge, one, tedge
               write (luncrt, '(/, a, /)') 'END FRAME'
            end if

!           Re-ue (unit Reynolds number based on edge conditions, per meter):

            fout(32,i,j) = (redge * vedge) / fout(23,i,j)  ! rho V / edge mu

!           Film coefficient:

            fout(33,i,j) = fout(11,i,j) / (fout(18,i,j) - fout(4,i,j))
                           ! Wall Qdot  / (Hedge - Hwall)


!           Re-kk (Reynolds # based on roughness height k and conditions at k):
!           -------------------------------------------------------------------

!           Re-kk requires interpolation of density & velocity at k:

            do l = 1, 4 ! Four quantities to interpolate

               id = k_indices(l) ! 1, 2, 3, 4 <-> input rho, u, v, w

               yline(1:nk_height_k) = f(id,i,j,1:nk_height_k)

               call lcsfit (nk_height_k, tline, yline, new, method, 1,         &
                            roughness_height, height_k_values(l), ypeval)
            end do

            V = sqrt (height_k_values(2)**2 + height_k_values(3)**2 +          &
                      height_k_values(4)**2)

!           Re-kk:

            id      = k_indices(5)          ! Index of input viscosity
            wall_mu = f(id,i,j,1)
            fout(34,i,j) = (height_k_values(1) * V * roughness_height) / wall_mu

            if (save_profile) then

               do k = 1, nk_height_k   ! Just the inner radial line grid points
                  V        = sqrt (f( 9,i,j,k)**2 + f(10,i,j,k)**2 +           &
                                   f(11,i,j,k)**2)
                  yline(k) = f(1,i,j,k) * V * tline(k) / wall_mu
               end do

               write (luncrt, '(/, (a))') 'Re-kk Profile', subtitle,           &
                                          'k, meters', 'Re-kk'
               write (luncrt, '(1p, 2e19.11)')                                 &
                  (tline(k), yline(k), k = 1, nk_height_k)

            end if

!           TPS tile ID and thickness?
!           !!!!!!!!!!!!!!!!!!!!!!!!!!

            if (tps_data) then

!              Locate the surface grid point in the OML tile triangulation:

               surface_xyz(1) = xout(i,j)
               surface_xyz(2) = yout(i,j)
               surface_xyz(3) = zout(i,j)

               call search_adt (surface_xyz, itri, p, q, r, dsq, true, nnode,  &
                                ntri, conn, tri_xyzs, interp_xyz)
!              !!!!!!!!!!!!!!!

!              Translate from triangle number to tile number.
!              Use the previous result as the starting guess for the index:

               call interval (nthick, tile_ids, tri_ids(itri), one, left)

               fout(numf_out-1,i,j) = tile_ids(left)
               fout(numf_out  ,i,j) = tile_thicknesses(left)

!!!            write (luncrt, '(a,3i5,a,i7,a,1p,e12.4,a,0p,f13.0,f10.6)')      &
!!!              ' ib,i,j: ', ib, i, j, '  itri:', itri, '  dsq:', dsq,        &
!!!              '  tile id, thickness:', tile_ids(left), tile_thicknesses(left)

               if (tile_ids(left) /= tri_ids(itri)) then
                  fout(numf_out,i,j) = -tile_thicknesses(left)
                  if (tile_diagnostics) then
                     id = int (tri_ids(itri))
                     n  = int (tile_ids(left))
                     write (luncrt, '(a, i13, a, i13, a, 1p, 4e14.6, a, 3i6)') &
                        ' Tile # missing:', id, ';  using -thickness of #', n, &
                        ';  dsq,x,y,z: ', dsq, surface_xyz, '  ib,i,j: ', ib,i,j
                  end if
               end if

            end if

         end do ! Next i

      end do ! Next j

      deallocate (xline, yline, zline, tline, fline, Hratio, rratio, vratio)

99    continue

      end subroutine process_block

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine SI_to_English (ni, nj, nf, tps, x, y, z, f)

!     Convert one block of BLAYER_RESULTS output.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in) :: ni, nj  ! Surface block dimensions
      integer, intent (in) :: nf      ! # flow quantities, not counting x,y,z
      logical, intent (in) :: tps     ! T if tile thicknesses are present
      real, intent (inout) :: x(ni,nj), y(ni,nj), z(ni,nj)
      real, intent (inout) :: f(nf,ni,nj)

!     Local constants:

      real, parameter ::                &
         m_to_in        = 39.37007874,  &
         m_to_ft        = 3.280839895,  &
         K_to_F         = 1.8,          &
         K_to_F_offset  = 459.67,       &
         Pa_to_psf      = 0.020885472,  &
         Jcm2_to_BTUft2 = 0.8807645051, &
         Jkg_to_BTUlb   = 0.0004300269, &
         kgcm2_to_lbft2 = 2048.161436,  &
         kgm3_to_lbft3  = 0.06242796

!     Local variables:

      integer :: i, j

!     Execution:

      do j = 1, nj
         do i = 1, ni
               x(i,j) =    x(i,j) * m_to_in
               y(i,j) =    y(i,j) * m_to_in
               z(i,j) =    z(i,j) * m_to_in
            f( 1,i,j) = f( 1,i,j) * kgm3_to_lbft3          ! density (wall)
            f( 2,i,j) = f( 2,i,j) * Pa_to_psf              ! pressure
            f( 3,i,j) = f( 3,i,j) * K_to_F - K_to_F_offset ! temperature
            f( 4,i,j) = f( 4,i,j) * Jkg_to_BTUlb           ! enthalpy
            f( 5,i,j) = f( 5,i,j) * Pa_to_psf              ! viscosity
            f(11,i,j) = f(11,i,j) * Jcm2_to_BTUft2         ! heat flux
            f(12,i,j) = f(12,i,j) * Pa_to_psf              ! shear stresses
            f(13,i,j) = f(13,i,j) * Pa_to_psf
            f(14,i,j) = f(14,i,j) * Pa_to_psf
            f(15,i,j) = f(15,i,j) * kgm3_to_lbft3          ! density (edge)
            f(16,i,j) = f(16,i,j) * Pa_to_psf              ! pressure
            f(17,i,j) = f(17,i,j) * K_to_F - K_to_F_offset ! temperature
            f(18,i,j) = f(18,i,j) * Jkg_to_BTUlb           ! enthalpy
            f(19,i,j) = f(19,i,j) * m_to_ft                ! velocity components
            f(20,i,j) = f(20,i,j) * m_to_ft
            f(21,i,j) = f(21,i,j) * m_to_ft
            f(23,i,j) = f(23,i,j) * Pa_to_psf              ! viscosity
            f(29,i,j) = f(29,i,j) * m_to_in                ! delta
            f(30,i,j) = f(30,i,j) * m_to_in                ! delta-*
            f(31,i,j) = f(31,i,j) * m_to_in                ! theta
            f(32,i,j) = f(32,i,j) / m_to_in                ! unit Re # per inch
            f(33,i,j) = f(33,i,j) * kgcm2_to_lbft2         ! film coefficient
            if (tps) &
            f(36,i,j) = f(36,i,j) * m_to_in                ! tile thickness
         end do
      end do

      end subroutine SI_to_English

   end program blayer_results

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine average_increments (ni, nj, nk, x, y, z, average_ds)

!  Calculate the average off-face grid spacings for faces 1:6 of a grid block.
!
!  08/23/05  David Saunders  Initial implementation, for identifying wall faces.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)                       :: ni, nj, nk    ! Block dims.
   real,    intent (in), dimension (ni,nj,nk) :: x, y, z       ! Block coords.
   real,    intent (out)                      :: average_ds(6) ! (1) => i = 1,
                                                               ! (2) => i = ni,
                                                               ! and so on
!  Local constants:

   real, parameter :: zero = 0.

!  Local variables:

   integer :: i, j, k

!  Execution:

   average_ds(:) = zero

!  The i direction:

   do k = 1, nk
      do j = 1, nj
         average_ds(1) = average_ds(1) + sqrt ( &
            (x(2,j,k) - x(1,j,k))**2 + &
            (y(2,j,k) - y(1,j,k))**2 + &
            (z(2,j,k) - z(1,j,k))**2   )
         average_ds(2) = average_ds(2) + sqrt ( &
            (x(ni,j,k) - x(ni-1,j,k))**2 + &
            (y(ni,j,k) - y(ni-1,j,k))**2 + &
            (z(ni,j,k) - z(ni-1,j,k))**2   )
      end do
   end do

   average_ds(1:2) = average_ds(1:2) / real (nj * nk)

!  The j direction:

   do k = 1, nk
      do i = 1, ni
         average_ds(3) = average_ds(3) + sqrt ( &
            (x(i,2,k) - x(i,1,k))**2 + &
            (y(i,2,k) - y(i,1,k))**2 + &
            (z(i,2,k) - z(i,1,k))**2   )
         average_ds(4) = average_ds(4) + sqrt ( &
            (x(i,nj,k) - x(i,nj-1,k))**2 + &
            (y(i,nj,k) - y(i,nj-1,k))**2 + &
            (z(i,nj,k) - z(i,nj-1,k))**2   )
      end do
   end do

   average_ds(3:4) = average_ds(3:4) / real (nk * ni)

!  The k direction:

   do j = 1, nj
      do i = 1, ni
         average_ds(5) = average_ds(5) + sqrt ( &
            (x(i,j,2) - x(i,j,1))**2 + &
            (y(i,j,2) - y(i,j,1))**2 + &
            (z(i,j,2) - z(i,j,1))**2   )
         average_ds(6) = average_ds(6) + sqrt ( &
            (x(i,j,nk) - x(i,j,nk-1))**2 + &
            (y(i,j,nk) - y(i,j,nk-1))**2 + &
            (z(i,j,nk) - z(i,j,nk-1))**2   )
      end do
   end do

   average_ds(5:6) = average_ds(5:6) / real (ni * nj)

   end subroutine average_increments

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine blayer_edge (n, Hratio, Wdistance, save_profile, subtitle, &
                           kedge, Hedge, tedge, ier)
!
!  Estimate the boundary layer edge along a grid line assumed to be normal to
!  the wall in the region of interest.  A curvature-based method is used for the
!  given curve in (H/Hinf, distance) space.  The hope is to improve on the
!  standard approach of defining the edge as where the total enthalpy reaches
!  a value such as 99.5% of the free stream value.
!
!  After normalizing the wall distances, the first peak in curvature for which
!  enthalpy ratio is at least some cut-off value such as 0.95 is located via a
!  1-D minimization in a narrow range found by evaluating curvature initially at
!  the data points only.  Plain curvature proves too awkward to deal with, so
!  a moderated form is used.  Details:
!
!  Boundary Layer Edge Strategy (Each Radial Line):
!
!     >  Ignore the inner handful of the grid points; use points 8:kmax
!     >  Calculate arc lengths t(1:n) and normalize them by t(n)
!     >  Normalize enthalpy values H by the value Hinf at k = kmax
!     >  Consider the curve of H/Hinf ("y") vs. normalized wall distance t ("x")
!        >  Locate the first index i1 for which H/Hinf > 0.95 (say)
!        >  Calculate spline coefficients for "y" vs. "x" on 8:kmax
!        >  Parametric spline derivatives can be used for curvature as follows,
!           but we have now switched to the non-parametric form (below):
!              |k(s)| = sqrt ((d2x/ds2)**2 + (d2y/ds2)**2)    with sign from
!              d2y/ds2 * dx/ds  -  d2x/ds2 * dy/ds
!           Non-parametric form with the opposite sign of the second derivative:
!              k(x) = -d2y/dx2 / [(1 + (dy/dx)**2)**3/2]
!        >  Evaluate curvature initially at the data points only, and locate
!           the index of the first peak away from the wall
!        >  Moderate the curvature in the neighborhood of this peak:
!              f(k) = [1 + |k|]**power  where power = -0.1, say
!              if k < 0, preserve sign information via f(k) <-- 2 - f(k)
!           and perform a minimization of f(k) in [kpeak-1, kpeak+1] with
!           respect to "x" = normalized wall distance t
!        >  The denormalized wall distance tedge corresponding to this
!           minimum is the desired estimate of boundary layer thickness
!        >  [Later:  To avoid the dependence of curvature on radial grid line
!           arc length, adjust the above result by choosing the location of
!           99.5% of the profile peak IN THAT NEIGHBORHOOD.]
!
!  History:
!
!  06/19/04  David Saunders  Initial implementation in modular form.
!                            Who proposed the peak curvature idea?
!                            Probably Dinesh Prabhu and/or James Reuther.
!  06/24/04    "       "     The cut-off value is problematic.  Lowering it to
!                            give cleaner lee-side results can hurt wind-side
!                            cases which are the only cases where the idea is
!                            applicable.  Also, starting the spline at the
!                            cut-off data point means the calculated curvatures
!                            are not independent of the cut-off value.
!                            Therefore: spline the whole (normalized) profile,
!                            or at at least the bulk of it, but use just the
!                            data defined by the cut-off value.
!  08/08/04    "       "     The plain curvature method (and parametric splines)
!                            was really doing little else but picking the best
!                            grid point, in spite of the minimization.
!                            Therefore work with a moderated form of curvature
!                            to refine the initial location of the peak.
!                            The non-parametric form appears to give somewhat
!                            cleaner results as well, possibly because the
!                            abscissas (normalized distances) are varying
!                            slowly and steadily, while the arc-lengths along
!                            the normalized curve tend to increase then get
!                            much smaller right in the region of interest.
!                            The coarsening of the grid where the boundary
!                            layer is thicker doesn't help either.
!  08/11/04    "       "     Hignore = 0.975 allowed spurious minor peaks in
!                            curvature to be selected.  Raised it to 0.98,
!                            but will this now miss some legitimate peaks?
!                            The real problem is a flattening in the GASP
!                            profiles centered around 99%, and it's not clear
!                            which of two curvature peaks to pick.  One seems
!                            too near the surface; the other too far away.
!                            DPLR profiles do not behave this way.
!  08/12/04    "       "     Seek another peak if refined Hpeak < Hignore.
!                            Introduced tlimit = 0.5 to limit some extremes
!                            likely with lee-side data (ignored anyway).
!  12/09/04    "       "     Added y' and y" to the curvature plot for the
!                            indicated profile, hoping to understand some noisy
!                            edge thickness calculations for a 2-D wedge.
!  12/10/04    "       "     These derivatives seem to have no useful features
!                            in the region of interest.  Therefore, stay with
!                            curvature, but introduce some explicit smoothing
!                            to help with cases observed where the curvature is
!                            essentially constant for several grid points.
!                            The earlier moderation of (now smoothed) curvature
!                            is retained because without it, the peak is too
!                            inclined to be very close to a grid point.
!  12/13/04    "       "     Shuttle results suffer if every profile has its
!                            curvature smoothed, even with one iteration.
!                            Therefore, confine the smoothing to profiles that
!                            have two or more curvatures within a few percent
!                            of each other at the relevant peak.
!  02/11/05    "       "     Odd results on the side of a 3-D blunt wedge near
!                            the nose were traced to redefining kleft as kpeak
!                            before a retry due to Hedge < Hignore.  Thanks
!                            again go to Ryan McDaniel's wedge case.  The choice
!                            of Hignore is still more arbitrary than we would
!                            like, as some well-behaved profiles have their
!                            peak curvature slightly below 0.980, say.  Set it
!                            to 0.975 now instead of 0.980.
!  01/12/06    "       "     Further study of Shuttle solutions suggest lowering
!                            this limit to 0.95.  Curve fit experiments suggest
!                            normalizing so that delta is around 0.1 is a
!                            reasonable choice if we iterate to account for the
!                            dependence of curvature on scaling/normalization.
!                            However, normalizing by something like twice the
!                            initial edge height estimate is clearly not right.
!                            It changes the result drastically because of how
!                            the scaling affects curvature.
!  01/14/06    "       "     Iterating is dubious: any rescaling of arc lengths
!                            to standardize the (scaled) delta amounts to doing
!                            different scaling for each profile, thus distorting
!                            the edge distribution.  Leave the option in but use
!                            maxiter = 1.  Better work-around for the dependence
!                            of curvature on scaling: once the neighborhood of
!                            the edge height is located, pick the location of
!                            99.5% of the profile peak IN THAT NEIGHBORHOOD.
!                            This matches the traditional method for clean 2-D
!                            data, and appears the best compromise in 3-D, where
!                            total enthalpy ratio can differ significantly from
!                            1 in the likely region of peak curvature.
!                            Ironically, it means the refinements built into the
!                            curvature method are redundant because only kedge
!                            is really needed now.  Since the cost is small, and
!                            plots of specified profiles depend on them, leave
!                            the refinements in for now.
!  06/19/08-    "       "    Mike Olsen suggested using 95% for stage 2 instead
!  06/24/08                  of 99.5%.  This would be less sensitive to the
!                            heuristic definition of "neighborhood."  However,
!                            in the absence of a stage 3 that recovers the 99.5%
!                            location for well-behaved profiles, we leave well
!                            alone.  Only minor refinements have been made.
!
!  Author:  David Saunders, ELORET/NASA Ames Research Center, Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer,   intent (in)  :: n            ! # pts. in the radial line of data

   real,      intent (in)  :: Hratio(n)    ! Enthalpy ratio values

   real,      intent (in)  :: Wdistance(n) ! Wall distances

   logical,   intent (in)  :: save_profile ! T means write plottable data to
                                           ! standard output
   character, intent (in)  :: subtitle*32  ! Annotates the particular profile

   integer,   intent (out) :: kedge        ! [kedge, kedge+1) contains tedge

   real,      intent (out) :: Hedge        ! Enthalpy ratio value at the BL edge

   real,      intent (out) :: tedge        ! Boundary layer thickness estimate

   integer,   intent (out) :: ier          ! 0 means no error

!  Local constants:

   integer, parameter :: iendl   = 0, &    ! Spline end conditions
                         iendr   = 0
   integer, parameter :: maxfun  = 100     ! Many fewer function evaluations
                                           ! are likely for the minimization
   integer, parameter :: maxiter = 1       ! Number of times to repeat the
                                           ! calculation, seeking to put delta
                                           ! at normalized arc length of 0.1;
                                           ! iteration 1 is the plain method
   integer, parameter :: lunerr  = 6       ! For diagnostics

   real,    parameter :: one     = 1.0, &
                         two     = 2.0, &
                         zero    = 0.0
   real,    parameter :: Hignore = 0.950   ! Ignore curvature data below this
                                           ! value of H/Hinf
   real,    parameter :: power   = -0.1, & ! Moderates curvature; reciprocal
                         rpower  = one/power
   real,    parameter :: small_fraction  & ! Criterion for smoothing curvature
                                 = 0.10
   real,    parameter :: tlimit  = 0.5     ! Upper limit on normalized wall
                                           ! distance allowed for edge value
   real,    parameter :: tol     = 1.e-8   ! Acceptable uncertainty in the
                                           ! BL edge arc length (s, not t)

   logical, parameter :: false   = .false.

   character, parameter :: subname * 7 = 'bl_edge'

!  Local variables:

   integer :: i, istat, k, k1, k2, kedge_interim, klast, kleft, kpeak,         &
              kpeak_previous, kwindow, l1, l2, lunout, nfit, nfit2, niter, nt, &
              numfun

   real    :: t1, t2, tmin, flattened_min, fmin, Hmax, stotal, term, ttotal,   &
              Hedge_interim, tedge_interim

   real, dimension (5) :: &
      b3, c3, d3, df, dff, flattened_kappa

   real, dimension (:), allocatable :: &
      b1, c1, d1, dy, dyy, f, kappa, t

!  Execution:
!  !!!!!!!!!!

!  Ignore data near the wall with enthalpy ratio < Hignore?  NO:  The curvature
!  values used are then not independent of Hignore.  Therefore, err on the side
!  of more than enough data very near the wall.

!!!   kleft = n / 2
!!!   call interval (n, Hratio, Hignore, one, kleft)

   kleft = 8

!  Ignore data beyond some fraction of the maximum wall distance, assuming the
!  calling program hasn't already suppressed the outer portion of the data:

   klast = (3 * n) / 4
   term  = Wdistance(n) * tlimit

   call interval (n, Wdistance, term, one, klast)

   nfit = klast - kleft + 1 ! Enough for working data, but allow for showing all
   nt   = n     - kleft + 1 ! the normalized distances in the saved profile

   allocate (t(nt), b1(nfit), c1(nfit), d1(nfit), f(nfit), &
             dy(nfit), dyy(nfit), kappa(nfit))

   ttotal = Wdistance(n) ! For the first pass where edge height is unknown
   niter  = 1

10 continue  ! Top of loop over arc length normalizations (but maxiter may be 1)

!  Normalize the wall distances:

   t(1:nt) = Wdistance(kleft:n) / ttotal

!  Calculate spline coefficients for enthalpy ratio vs. t:

   call csfit (nfit, t, Hratio(kleft), iendl, zero, iendr, zero, b1, c1, d1,ier)

   if (ier /= 0) then
      write (lunerr, '(/, 2a, i4)') subname, &
         ':  Trouble splining enthalpy ratio vs. t; ier:', ier
      go to 99
   end if

!  Compute 1st and 2nd derivatives of enthalpy ratio at the data points:

   call csdval (nfit, t, Hratio(kleft), nfit, t, b1, c1, d1, f, dy, dyy)

!  Derive the curvature at the data points:

   do k = 1, nfit
      term = one + (min (abs (dy(k)), 1.e+10)) ** 2
      kappa(k) = -dyy(k) / (term * sqrt (term)) ! Sign is consistent
   end do                                       ! with parametric version

!  Explicit smoothing of these curvatures can help smooth variation of
!  results along the wall in regions where curvature is almost constant
!  for several points along the current profile.  However, smoothing for
!  all profiles is not good either (undesirable smearing?).  Therefore,
!  do it only if the peak curvature has a close neighbor.

!  Locate the first data point where enthalpy ratio is above the cut-off
!  value and curvature is at a local maximum:

   kpeak = nfit
   do k = 1, nfit - 1
      if (Hratio(k+kleft-1)  < Hignore) cycle
      if (kappa(k+1) < kappa(k)) then
         kpeak = k
         exit
      end if
   end do

   k1   = max (1, kpeak - 1)
   k2   = min (kpeak + 1, nfit)
   term = max (kappa(k1), kappa(k2))

   if (abs (kappa(kpeak) - term) < kappa(kpeak) * small_fraction) then

      call smooth1d (1, nfit, t, kappa)

!     Re-locate the relevant local maximum in the smoothed curvature:

      kpeak = nfit
      do k = 1, nfit - 1
         if (Hratio(k+kleft-1)  < Hignore) cycle
         if (kappa(k+1) < kappa(k)) then
            kpeak = k
            exit
         end if
      end do

   end if

20 continue ! Return here if Hedge < Hignore

   k1 = max (1, kpeak - 1);     t1 = t(k1)
   k2 = min (kpeak + 1, nfit);  t2 = t(k2)

   if (save_profile) then
      lunout = lunerr  ! For showing minimization iterations
      write (lunout, '(/, a)') '! For specified (H/Hinf, distance) profile:'
      write (lunout, '(a, i4, 1p, 2e19.11)') &
         '! kpeak, t1, t2:', kpeak+kleft-1, t1, t2
   else
      lunout = -lunerr
   end if

!  Refine the peak in this reduced interval by performing a minimization.
!  Plain curvature is too spiky to work with directly.  Moderate it for the
!  points in the neighborhood of the relevant peak.  Preserve the sign of the
!  curvature in the process.

   l1 = max (1, kpeak - 2)
   l2 = min (kpeak + 2, nfit)
   nfit2 = l2 - l1 + 1

   k = l1
   do i = 1, nfit2
      flattened_kappa(i) = (abs (kappa(k)) + one) ** power
      if (dyy(k) > zero) flattened_kappa(i) = two - flattened_kappa(i)
      k = k + 1
   end do

!  Now we work with a new spline (that of moderated curvature vs. norm. dist.),
!  and we want to allow some overshoot, so we DON'T use a monotonic spline:

   call csfit (nfit2, t(l1), flattened_kappa, iendl, zero, iendr, zero,        &
               b3, c3, d3, ier)
   if (ier /= 0) then
      write (lunerr, '(/, 2a, i4)') subname, &
         ':  Trouble splining moderated curvature vs. norm. distance; ier:', ier
      go to 99
   end if

!  Locate the minimum of the moderated, splined curvature in [kpeak-1, kpeak+1]:

   istat  = 2          ! Initialization flag
   numfun = maxfun

   do ! Until convergence

      call fminrc (t1, t2, tmin, fmin, tol, numfun, subname, lunout, istat)

      if (istat < -1) then ! Fatal error

         write (lunerr, '(/, 2a)') subname, &
            ': Trouble refining peak curvature location.'
         ier = -1
         exit

      else if (istat < 0) then ! Is MAXFUN too low?

         write (lunerr, '(/, 2a)') subname, &
            ': Unconverged minimization, but proceeding.'
         tedge = Wdistance(kpeak + kleft - 1)
         Hedge = Hratio(kpeak + kleft - 1)
         exit

      else if (istat > 0) then ! Evaluate moderated curvature at tmin

         call csdval (nfit2, t(l1), flattened_kappa, 1, tmin, b3, c3, d3,      &
                      fmin, df, dff)

      else ! istat = 0 means convergence to tmin; assign other edge values

         flattened_min = fmin
         fmin = flattened_min ** rpower - one    ! Equivalent smoothed kappa pk.

         call csdval (nfit, t, Hratio(kleft), 1, tmin, b1,c1,d1, Hedge, dy, dyy)

!!!      term = one + (min (abs (dy(1)), 1.e+10)) ** 2
!!!      fmin = -dyy(1) / (term * sqrt (term))   ! Equivalent plain kappa peak

         tedge = tmin * ttotal                   ! Denormalize edge distance
         kedge = kleft + kpeak                   ! Ensure kedge pnts. to the pt.
                                                 ! to the left of the b.l. edge
         call interval (n, Wdistance, tedge, one, kedge)
         exit

      end if

   end do ! End of minimization iterations

   if (Hedge < Hignore) then ! Must be a spurious early peak; seek another

      kpeak_previous = kpeak
      do k = kpeak + 1, nfit - 1
         if (kappa(k+1) < kappa(k)) then
            kpeak = k
            exit
         end if
      end do
      if (kpeak > kpeak_previous) then
!!!      write (lunerr, '(a, f8.5, a, f8.5, a)') &
!!!         ' Dubious edge found at (', Hedge, ',', tedge, '); retry.'
         go to 20
      end if

   end if

   if (save_profile) then
      tmin = tedge / ttotal
      write (lunout, '(a, i2, a, f10.6)') &
         ' Normalized delta at iteration', niter, ':', tmin
   end if

   if (niter < maxiter) then ! Iterate to achieve normalized edge height ~ .10
       niter = niter + 1
      ttotal = tedge / 0.10
      go to 10
   end if

!  Now that we have the most likely neighborhood of the edge height, adjust
!  for possible anomalies and for the arc-length dependency by seeking 99.5%
!  of the PEAK enthalpy ratio in this neighborhood.
!  Actually, since 99.5% puts us in the ill-conditioned region where small
!  differences in local peak can make large differences in edge height, it
!  would be helpful to err on the low but less sensitive side.  Mike Olsen
!  has proposed using 95% with some empirical way of shifting that result to
!  the likely edge.  In the absence of such an empirical shift, we leave it
!  at 99.5%, knowing the heuristic definition of "neighborhood" can lead to
!  discontinuities in the estimated boundary layer thickness distribution
!  where the flow field is more complex (e.g., Shuttle wing tip region).  We
!  don't want to compromise clean windside flow regions by trying to improve
!  behavior in awkward regions that are of less interest anyway.

   kwindow = nint (real (n) * 0.06)
   kwindow = max (5, min (15,    kwindow))
   k1      = max (   10, kedge - kwindow)
   k2      = min (klast, kedge + kwindow)
   Hmax    = zero

   do k = k1, k2
      Hmax = max (Hratio(k), Hmax)
   end do

   Hedge_interim = Hedge
   tedge_interim = tedge
   kedge_interim = kedge

   Hedge = Hmax * 0.995

!  Make sure that k1 isn't too high to include the Hedge:

   do k = kedge - 1, 10, -1
      if (Hratio(k) <= Hedge) then
         k1 = k
         exit
      end if
   end do

   call bl_edge_traditional (Hedge, n, k1, k2, Hratio, Wdistance, kedge,       &
                             tedge, ier)

   if (ier /= 0) then
      Hedge = Hedge_interim
      tedge = tedge_interim
      kedge = kedge_interim
   end if

   if (save_profile) call save_plottable_details () ! Local procedure below

   deallocate (b1, c1, d1, dy, dyy, f, kappa, t)

99 return

!  Internal procedure for routine blayer_edge:

   contains

!     --------------------------------------------------------------------------

      subroutine save_plottable_details ()

!     Save the details of the current profile in QPLOTable form.
!     QPLOT was developed at NASA Ames around the CA-DISSPLA package,
!     and is thus not easily ported, but this procedure could be
!     adapted to suit some other X-Y plot package.

!     --------------------------------------------------------------------------

      write (lunout, '(a)') ' ',                    &
         'Boundary Layer Profile',                  &
         subtitle,                                  &
         'Enthalpy ratio',                          &
         'Wall distance, m',                        &
         ' $options',                               &
         ' width=6.5, height=7,',                   &
         ' xmin=0., xmax=1.1, xstep=0.1, grid=1,',  &
         ' line=''solid'', symbol=16,',             &
         ' spline=''parametric'', fit=''loose'', ', &
         ' ident=T,',                               &
         ' $end'
      write (lunout, '(/, a)') '!      H/Hinf            Wall distance   k'
      write (lunout, '(1p, 2e19.11, i4)') (Hratio(k), Wdistance(k), k, k = 1, n)

      write (lunout, '(a)') &
         ' $options line=''symbol'', symbol=15,'
      write (lunout, '(a, f8.5, a, f8.5, a)')       &
         ' legend = ''Boundary layer edge is at (', Hedge, ',', tedge, ' )'',',&
         ' $end'
      write (lunout, '(/, a, 2i4, /, 1p, 2e18.11)') &
         '! kedge, kleft, Hedge, tedge:', kedge, kleft, Hedge, tedge
      write (lunout, '(/, a, /)') 'END FRAME'

      write (lunout, '(a)') &
         'Boundary Layer Profile',                  &
         subtitle,                                  &
         'Enthalpy ratio',                          &
         'Normalized wall distance',                &
         ' $options',                               &
         ' height=7, xstep=0.1, ystep=0.1, grid=1,',&
         ' line=''solid'', symbol=16,',             &
         ' spline=''parametric'', fit=''loose'', ', &
         ' plot=''scale'',',                        &
         ' ident=T,',                               &
         ' $end'
      write (lunout, '(/, a)') '!      H/Hinf    Normalized Wall Dist.   k'
      write (lunout, '(1p, 2e19.11, i4)')           &
         (Hratio(k), t(k-kleft+1), k, k = kleft, n)

      write (lunout, '(a)') &
         ' $options line=''symbol'', symbol=15,'
      write (lunout, '(a, f8.5, a, f8.5, a)')       &
         ' legend = ''Boundary layer edge is at (', Hedge, ',', tmin, ' )'',', &
         ' $end'
      write (lunout, '(/, a, 2i4, /, 1p, 2e18.11)') &
         '! kedge, kleft, Hedge, tmin:', kedge, kleft, Hedge, tmin
      write (lunout, '(/, a, /)') 'END FRAME'

      write (lunout, '(a)') &
         'Boundary Layer Profile',                  &
         subtitle,                                  &
         'Normalized wall distance',                &
         'Derivatives/curvature of H/Hinf vs. normalized wall distance',       &
         ' $options height=7, ymin=0, ymax=50,',    &
         ' line=''solid'', symbol=16, ident=T,',    &
         ' legend = ''Smoothed curvature'', $end',  &
         '!         t                  Curvature   k'
      write (lunout, '(1p, 2e19.11, i4)') &
         (t(k), kappa(k), k+kleft-1, k = 1, nfit)

      write (lunout, '(a)') &
         ' $options line=''dash'',',                &
         ' legend = ''First derivative'', $end',    &
         '!         t                      dy/dt   k'
      write (lunout, '(1p, 2e19.11, i4)') &
         (t(k), dy(k), k+kleft-1, k = 1, nfit)

      write (lunout, '(a)') &
         ' $options line=''longdash'',',            &
         ' legend = ''Second derivative'', $end',   &
         '!         t                    d2y/dt2   k'
      write (lunout, '(1p, 2e19.11, i4)') &
         (t(k), -dyy(k), k+kleft-1, k = 1, nfit)

      write (lunout, '(a)') &
         ' $options line=''symbol'', symbol=15,',   &
         ' legend = ''BL edge curvature, smoothed'', $end'
      write (lunout, '(1p, 2e19.11)') tmin, fmin
      write (lunout, '(/, a, /)') 'END FRAME'

      write (lunout, '(a)') &
         'Boundary Layer Profile',                  &
         subtitle,                                  &
         'Normalized wall distance',                &
         'Moderated curvature in region of peak',   &
         ' $options height=7, fit=''loose'',',      &
         '  spline=''standard'', ident=T, $end',    &
         '!         t        Moderated Curvature   k'
      write (lunout, '(1p, 2e19.11, i4)') &
         (t(l1+k-1), flattened_kappa(k), k+l1+kleft-2, k = 1, nfit2)

      write (lunout, '(a)') &
         ' $options line=''symbol'', symbol=15,',   &
         ' legend = ''Boundary layer edge'', $end'
      write (lunout, '(1p, 2e19.11)') tmin, flattened_min
      write (lunout, '(/, a, /)') 'END FRAME'

      end subroutine save_plottable_details

   end subroutine blayer_edge

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine bl_edge_traditional (Hedge, nk, k1, k2, Hratio, Wdistance,       &
                                   kedge, delta, ier)
!
!  Traditional approach to defining a boundary layer edge: determine the height
!  at which total enthalpy first achieves some percentage (usually 99.5%) of the
!  free-stream value.  Spline interpolation is used.  See also blayer_edge for a
!  curvature-based method which may be preferable.
!
!  12/01/05  D.A.Saunders  Initial implementation.
!  01/11/06    "     "     Limit the data range at the higher level via k1, k2.
!  07/07/06    "     "     The traditional edge method was not trapping non-
!                          monotonic enthalpy ratios properly.
!  07/21/06    "     "     Avoiding inverse interpolation was a bad idea: the
!                          4-point spline method does not degrade as gracefully
!                          to 3 or 2 points as was assumed.
!
!  Author:  David Saunders, ELORET/NASA Ames Research Center, Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   real,      intent (in)  :: Hedge         ! Usually 0.995 (target Hratio)
   integer,   intent (in)  :: nk            ! # pts. in the radial line of data
   integer,   intent (in)  :: k1, k2        ! Data range of interest
   real,      intent (in)  :: Hratio(nk)    ! Enthalpy ratio values
   real,      intent (in)  :: Wdistance(nk) ! Wall distances
   integer,   intent (out) :: kedge         ! [kedge, kedge+1) contains delta
   real,      intent (out) :: delta         ! Boundary layer thickness estimate
   integer,   intent (out) :: ier           ! 0 means no error

!  Local constants:

!! logical,   parameter :: new = .true.     ! Always new data for LCSFIT
!! character, parameter :: method * 1 = 'M' ! Monotonic local spline method

!  Local variables:

   integer :: k, kleft
!! integer :: k, klast, kleft, nfit
!! real    :: unused

!  Execution:

!  Locate the first interval containing the target enthalpy ratio by sequential
!  search because we can't be sure the ratios are monotonic.

   ier = 1

   do k = k1, k2 - 1
      if (Hratio(k) <= Hedge .and. Hedge < Hratio(k+1)) then
         kedge = k   ! kedge is to the left of the peak now
         ier   = 0
         exit
      end if
   end do

   if (ier /= 0) then
      kedge = k2
      delta = Wdistance(kedge) ! Probably garbage
      go to 99
   end if

!  Avoid inverse nonlinear interpolation but guard against nonmonotonic Hratios.
!  Instead of reducing LCSFIT's normal 4-point data range if the abscissas are
!  not increasing, expand the data range from 2 points, but only if we can.

!! kleft = kedge           ! We need no more than 4 points and have at least 2
!! if (Hratio(kleft-1) < Hratio(kleft)) kleft = kleft - 1

!! klast = kedge + 1
!! if (Hratio(klast) < Hratio(klast+1)) klast = klast + 1

!! nfit = klast - kleft + 1

!! call lcsfit (nfit, Hratio(kleft), Wdistance(kleft), new, method, 1, Hedge,  &
!!              delta, unused)

!  No!  Degrading from 4 to 3 or 2 points is not as graceful as was assumed,
!  especially for profiles that have overshoots.  Swap the roles of x and y.
!  Wall distances are known to increase steadily, so we apply an inverse method
!  in the interval now known to contain the target Hedge:

   kleft = kedge - 1  ! Points kleft, kedge, kedge+1, kedge + 2 (only) are used

   call lcsinverse (Wdistance(kleft), Hratio(kleft), Hedge, delta)

99 return

   end subroutine bl_edge_traditional

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE SMOOTH1D (I1, I2, S, V)
!
!     SMOOTH1D smooths a vector in-place by an explicit method with hard-coded
!     parameters.  End elements are unchanged.
!
!     06/23/97  D.Saunders  Gridgen-type explicit smoothing to cope with spikes
!                           in edge Phi, Psi.
!     12/10/04      "       5 iters. seem too many for boundary layer edge work.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IMPLICIT NONE

!     Arguments:

      INTEGER, INTENT (IN)    :: I1, I2   ! Index range to be smoothed
      REAL,    INTENT (IN)    :: S(I1:I2) ! Abscissas associated with V(I1:I2)
      REAL,    INTENT (INOUT) :: V(I1:I2) ! Vector being smoothed; elements
                                          ! I1 & I2 are left unchanged
!     Local constants:

      INTEGER, PARAMETER :: NITERS = 1    ! Originally 5
      REAL,    PARAMETER :: ALPHA = 0.5, ONE = 1.

!     Local variables:

      INTEGER  I, ITER
      REAL     AI, BI, CI, HL, HR, TERM
      REAL     W(I1:I2) ! Work-space for previous iterate (to vectorize)

!     Execution:

      DO ITER = 1, NITERS

         W(I1:I2) = V(I1:I2)

         DO I = I1 + 1, I2 - 1
            HL   = S(I) - S(I-1)
            HR   = S(I+1) - S(I)
            TERM = ALPHA * ((MIN (HL, HR) ** 2) / (HL + HR))
            AI   = TERM / HL
            CI   = TERM / HR
            BI   = ONE - (AI + CI)
            V(I) = AI * W(I-1) + BI * W(I) + CI * W(I+1)
         END DO

      END DO

      END SUBROUTINE SMOOTH1D
