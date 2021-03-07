!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program blayer
!
!  Description:
!
!     BLAYER reads a volume dataset from a 2-D or 3-D multiblock real gas flow
!  solution and derives a one-layer set of results, some of which apply to the
!  wall, some to the edge of the boundary layer, and (optionally) some to the
!  roughness height or (if that height is entered as -1. or -2.) to the momentum
!  thickness height.
!
!     BLAYER is a generalization of BLAYER_RESULTS without the handling of TPS
!  tile data.  It has been adapted from BLAYER2D and BLAYER3D, which were first
!  adapted from BLAYER_RESULTS.  A variable number of species is handled, along
!  with (possibly) more than one temperature and additional flow quantities
!  beyond those that are required as inputs.
!
!     This version reads Tecplot files (ASCII only), or PLOT3D grid and function
!  files (ASCII or binary), and writes Tecplot ASCII or binary files.  The
!  merger of BLAYER2D and -3D was prompted by an upgrade of the Tecplot_io
!  package to handle Tecplot 360 auxiliary data records.
!
!  Assumptions:
!
!     >  The surface is normally expected to be at k = 1 for all grid blocks,
!        or j = 1 for the 2-D case.  [Internally, 2-D arrays are treated as if
!        they are dimensioned (1:ni, 1, 1:nk).]  However, for each block, the
!        program does look for the face or edge with the smallest average
!        initial increment off it, and if necessary permutes the block to put
!        the wall at the k = 1 position.
!        [Overriding the automated wall detection is also possible - see below.]
!        If necessary, blocks can be suppressed via the optional blayer.inp.2
!        file.
!
!     >  The radial lines are sufficiently normal to the surface for the 1-D
!        boundary layer edge detection method to make sense.
!
!     >  The flow quantities are given in SI units, at the grid-points, not cell
!        centers, although the latter should still be OK.
!
!     >  Lee-side results will not be used:  any boundary layer edge detection
!        method used is likely to be inappropriate for separated flows.
!
!     LEE-SIDE NOTE:  In spite of the disclaimer about detecting the edge in
!                     wake regions, an option has been provided to save the
!                     calculated edge locations (and thicknesses) as a second
!                     surface dataset to help visualize results with dubiously
!                     large edge thicknesses.  See dual use of the output
!                     datapacking control below to activate this option.
!
!  Tecplot ASCII Input Format:
!
!            TITLE     = ""
!            VARIABLES = "x, m"
!            "y, m"
!           ["z, m"]          for the 3-D case
!        1:  "rho, kg/m^3"    density
!        2:  "p, Pa"          pressure
!        3:  "T, K"           translational temperature
!       [ :  "Tv, K"          vibrational temperature; the default is 1 temp.]
!         :  "c_N_2"          1 or more species densities
!         :  "c_O_2"
!         :  "c_N_O"
!         :  "c_N"
!         :  "c_O"
!       [ :   :::             more species if specified; the default is 5 sp.]
!         :  "u, m/s"         velocity components
!         :  "v, m/s"
!         : ["w, m/s"]        for the 3-D case
!         :  "H0, J/kg"       total enthalpy
!         :  "M"              Mach number
!         :  "mu, Pa.s"       viscosity
!       [ :  "kappa, W/m.K"   total thermal conductivity if # temperatures > 1]
!       [ :   :::             miscellaneous extras; the default is 0 extras]
!            ZONE T="Zone 1"
!             I=17, J=25, [K=81, ]ZONETYPE=Ordered
!             DATAPACKING=BLOCK
!             DT=(SINGLE SINGLE .................... SINGLE SINGLE )
!             6.64733315E+00 6.57563824E+00 6.48970469E+00 6.39077472E+00 6. ...
!              :              :              :              :              :
!
!        Also, the format written by DPLR's Postflow is handled:
!
!             [Optional title]
!             variables=x,y,[z,]rho,p,T,C_n2,C_o2,C_no,C_n,C_o,u,v,[w,]h,M,mu
!             zone t="flow2|3d" F=point, i= 161 j= 157 k=  1|81
!             6.64733315E+00 6.57563824E+00 6.48970469E+00 6.39077472E+00 6. ...
!              :              :              :              :              :
!
!  PLOT3D Input Format (ASCII or Binary):
!
!        File pairs xxx.g/xxx.f or (preferably) xxx.gu/xxx.fu are implied if the
!        control file contains (in place of a Tecplot xxx.dat input volume file)
!        a file named xxx.g or xxx.gu.
!
!        For DPLR users, the POSTFLOW file should use output format 3 and
!
!           ivarp =  0 100 110 120     1000 150 151 [152] 132 154 50    [extras]
!        or ivarp =  0 100 110 120 125 1000 150 151 [152] 132 154 50 52 [extras]
!
!        for single-temperature and two-temperature solutions, respectively,
!        analogous to the Tecplot format.
!
!  Output Results (Tecplot ASCII or Binary File, One Zone Per Grid Block):
!
!            Wall               Boundary layer edge     Roughness height (k > 0)
!
!         x                     density                 height k
!         y                     pressure                density
!         s | z                 temperature             |velocity|
!         density               total enthalpy          viscosity
!         pressure              u                       Re-kk
!         temperature           v
!       [ Tvw ]               [ w ]                     or (if k = -1. or -2.):
!         total enthalpy        Mach number
!         viscosity             viscosity               Theta height values
!         N2 species density    N2 species density
!         O2    "       "       O2    "       "         k (= theta)
!         NO    "       "       NO    "       "         density    at theta
!         N     "       "       N     "       "         |velocity| at theta
!         O     "       "       O     "       "         viscosity  at theta
!       [ ??    "       "  ]  [ ??    "       "  ]      Re-theta     (k = -1.) |
!         heat flux             delta                   Re-theta/Medge (= -2).
!         tau_x                 delta* or vel-thickness, depending on input k
!         tau_y                 theta                   (see more on k below)
!       [ tau_z  ]              Re-ue
!       [ kappaw ]              CH
!       [ extras ]            [ kappae ]
!                             [ extras ]
!
!  Optional Additional Output File (blayer_edge_surface.dat):
!
!         x      Calculated boundary layer edge coordinates ...
!         y
!         z
!         delta  ... and thicknesses
!
!         Enter DATAPACKING = 21 or 22 instead of 1 or 2 to invoke this option.
!
!     Notes:
!
!        1:  Input and output files are assumed to be binary if their names end
!            in 'plt', else they are assumed to be ASCII (normally with a 'dat'
!            extension).  However, binary input is presently not an option.
!
!        2:  The program distinguishes 2-D/3-D data by reading the input file
!            header and looking for 'z' or 'Z' as the third variable name.
!
!        3:  The 2-D case outputs running length in place of z, so the outputs
!            are written as though they were 3-D.  Note that running length is
!            along the wall from the upstream edge of the current block.
!            Obtaining cumulative running lengths across blocks is grid-specific
!            but could presumably be derived within Tecplot from these outputs.
!
!            NOTE CONCERNING ARC LENGTHS:
!            Program SORT_SURFACE_SLICE available from the present author can
!            determine cumulative arc lengths across block boundaries.  First,
!            use Tecplot to slice the dataset (typically at z = 0.001 m for a
!            centerline cut).  The sorting program can assemble the slice into
!            a contiguous curve (or curves) and optionally insert arc lengths.
!            These may also be treated as run lengths from the stagnation point
!            if requested under some circumstances suited to capsules.
!
!            NOTE CONCERNING SHOCK STAND-OFF DISTANCES:
!            CEV practice is to add shock stand-off distance as an extra
!            surface quantity appended to the BLAYER output via programs
!            SHOCK_STAND_OFF and MERGE_FILES available from the present author.
!
!        4:  Heat flux is not among the data inputs (problematic in the volume).
!            At the wall, derive it from the surface temperature and the single
!            emissivity input.  For cold wall applications, DPLR's POSTFLOW can
!            help.
!
!            NOTE:  Ryan McDaniel can explain how to use a Tecplot macro to
!            overwrite BLAYER heat fluxes with those from POSTFLOW.
!
!        5:  Re-ue is the unit Reynolds # based on edge conditions.
!            To obtain Re-theta or Re-ke, multiply it by theta or k (but see 7).
!
!        6:  Re-kk is Reynolds # based on roughness height k & conditions at k,
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
!        7.  Re-theta [/ edge Mach number] is required commonly enough that
!            provision has been made for including it among the outputs if the
!            roughness-height quantities aren't needed:  enter k = -1. to
!            give Re-theta plus other quantities at height theta, or k = -2.
!            if Re-theta / Medge is more appropriate (Medge > 1.).
!
!        8.  Velocity thickness can be specified in lieu of delta* (displacement
!            thickness) as follows [awkward, but avoids a new input]:
!
!            Input k = 20. + intended k >= 0.: same as for k >= 0. otherwise;
!            Input k = -21. or -22.:           same as for k = -1.|-2. otherwise
!
!        9.  Program EXTRACT_BLAYER_DATA is available from the present author
!            to help tabulate certain flow solution results from CFD points
!            along a trajectory in readiness for curve fitting and padding to
!            full time histories as needed for TPS sizing at some body point(s),
!            one body point per run.  Simply enter a control file like this:
!
!            1 163                  ! (iblock, i) for a 2-D case
!            5 26 6 58 32 7 27 28   ! pw, qconv, Tw, CH, Hedge, Hw, tauwx, tauwy
!            MHL-t52.0/blayer.dat   ! Any number of BLAYER output files, to EOF
!            MHL-t55.3/blayer.dat
!            MHL-t57.9/blayer.dat
!            MHL-t61.2/blayer.dat
!            MHL-t64.2/blayer.dat
!            MHL-t65.8/blayer.dat
!            MHL-t70.0/blayer.dat
!
!            This example would also append |tauw| as an additional column in
!            the output, which is designed for MERGE_TABLES.
!
!  Control File (standard input):
!
!     TITLE
!     INPUT VOLUME DATASET
!     ***.dat               ! Tecplot ASCII or ***.g|***.gu PLOT3D ASCII|binary
!     OUTPUT TECPLOT FILE
!     ******.plt            ! Wall and boundary layer edge results
!     2                     ! 1 = DATAPACKING=POINT; 2 = DATAPACKING=BLOCK
!     MISCELLANEOUS CONTROLS
!     0                     ! Edge method: three choices are explained below
!     0.89                  ! Emissivity (0.89 = RCG everywhere)
!     0.001524              ! Roughness ht. k, m (=0.060" for Shuttle) |0|-1|-2
!     1  1  1               ! Block #, i, j for one Hnorm & one Re-kk profile
!     5                     ! # species      [5 if omitted]
!     1                     ! # temperatures [1 if these numbers are omitted]
!     0                     ! # extra items  [0 "    "     "     "     "    ]
!     0.                    ! Optional input for Hinf, but see extensions below
!     0.                    ! Optional %Hnorm > 0. value for Hignore | 95%
!     0.                    ! Optional Hshift; use Hform(0K) = 8932880. for CO2
!                           ! Hshift is ignored if the total enthalpy ratio is
!                           ! (H - Hwall) / (Hinf - Hwall); see edge_method < 0.
!
!  Further Control Explanations:
!
!     Several controls now have multiple uses because of the need to remain
!     consistent with existing control files.  A different input scheme, such
!     as a namelist, would have been preferable.
!
!     Edge Method (EM) Usage:
!
!        EM   >= 0.  means use the enthalpy ratio (H + Hshift) / (Hinf + Hshift)
!        EM   <  0.  means use the enthalpy ratio (H - Hwall)  / (Hinf - Hwall)
!       |EM|  < 90.  means use the hybrid curvature-based method
!       !EM|  = 99.5 means use the traditional 99.5% method (or 99. or whatever)
!  Thus  EM   = 99.5 means hybrid with Hshift profile (H/Hinf for air/Hshift=0.)
!        EM   = -1.  means hybrid with Hwall profile
!        EM   = -99. means traditional 99% method with Hwall profile
!        
!     Hinf Usage:
!
!        Hinf >  0. =>  use that value to normalize all profiles
!        Hinf =  0. =>  use the value from block 1, point (1,1,nk)   [default]
!        Hinf = -1. =>  use the (i,j,nk) value for profile (i,j) of each block
!        Hinf = -2. =>  use the peak value along each profile
!
!     Roughness Height k Usage:
!
!        See Note 7 above for the dual use of the roughness height input, k.
!
!  Optional Ancillary Control File blayer.inp.2:
!
!                      Line 1:  Blocks to be suppressed
!
!     This file can be used to suppress boundary layer processing of specified
!     grid blocks.  If present, it should contain a list of block numbers to
!     suppress, in any obvious convenient form, on a single line.  Examples:
!
!        11:16 or 11-16   would both expand to 11, 12, 13, 14, 15, 16
!        10 12 14:20      or any other such intelligible list, on line 1
!
!     This first line can be empty, meaning no blocks are suppressed.
!     If the file is not present, no blocks are suppressed.
!     At present, "suppressed" means the blocks are not processed, but dummy
!     results are written as output.  This avoids possible confusion resulting
!     from renumbering of blocks in the output.
!
!            Lines 2-n:  Overrides for the automated wall detection scheme
!
!     A wing leading edge plug case encountered spacings finer off the faces
!     surrounding the plug than off the plug, so the following override scheme
!     has been provided.  Enter one block/face pair per line starting at line 2.
!     For example:
!
!        10:13   [or blank line 1 if no blocks are being suppressed]
!        2 5     [use the k = 1 face for blocks 2:5]
!        3 5
!        4 5
!        5 5
!
!  Optional Ancillary Variable Names File blayer.inp.3 For PLOT3D Input Cases:
!
!     If the input volume dataset is in PLOT3D form, the variable names for the
!     species densities and any extra variables may be supplied or defaulted.
!     To supply them, provide a file named blayer.inp.3 with one or two lines.
!     The first line should contain all the species names.  The second line
!     should contain the names of any extra variables.  All the other input
!     variables are in known locations and hence can be hard-coded in BLAYER.
!
!     Example blayer.inp.3 (11 species on line 1; pressure coef. and total
!                           mixture number density as extras on line 2):
!
!        n2 o2 no no+ n2+ o2+ n o n+ o+ e
!        C_p N_tot
!
!     If either line is empty, the variable names are defaulted like this:
!
!        sp_1 sp_2 sp_3 ..... sp_11
!        xtra_1 xtra_2
!
!  Boundary Layer Edge Strategy (Each Radial Line):
!
!     See subroutine blayer_edge below.  Briefly, a curvature-based scheme is
!     use to locate the most likely edge region, then the traditional method
!     (99.5% or 99%) is applied as a fraction of the PEAK IN THAT NEIGHBORHOOD.
!
!  Procedures:
!
!     Tecplot_io package     I/O module  for Tecplot files
!     xyq_io and xyzq_io     I/O modules for PLOT3D  files
!     <numerous numerics utilities>
!
!  History:
!
!     06/19/04  D. Saunders  Initial implementation of BLAYER_RESULTS.
!     08/23/05   "     "     Last enhancement to BLAYER_RESULTS before the
!                            2-D translation.  (Introduced PERMUTE_BLOCK to
!                            handle the wall at any block face.)
!     11/05/05   "     "     Initial 2-D translation (5 species, not TPS tile
!                            handling, and no PLOT3D file I/O option).
!     11/15/05   "     "     Generalized # species, temperatures, and extras.
!     11/30/05   "     "     BLAYER3D adapted from BLAYER2D.
!     12/01/05   "     "     Introduced "edge_method" control to allow easy
!                            comparison with traditional 99[.5]% method.
!     12/07/05   "     "     Heat flux units are now W/m^2.
!     12/08/05   "     "     Delta* and theta formulations now use tangential
!                            velocities, not total velocity magnitudes;
!                            careful safeguarding of a stag. point is required.
!     01/17/06   "     "     Following experiments with curve fitting and
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
!     02/03/06   "     "     Momentum thickness had an extra factor of
!                            |vtangent|/|vtangent(edge)| in the integral. Thanks
!                            to Frank Greene for noticing it looked wrong.
!     03/14/06   "     "     A wing leading edge plug case picked the wrong face
!                            of blocks above the plug.  Therefore, make use of
!                            the optional ancillary control file to override the
!                            automated smallest-average-increment scheme.
!     07/07/06   "     "     The traditional edge method was not trapping non-
!                            monotonic enthalpy ratios properly.
!     07/21/06   "     "     Avoiding inverse interpolation in the traditional
!                            edge method was a bad idea: the 4-point spline
!                            method does not degrade as gracefully to 3 or 2
!                            points as was assumed.  Profiles with overshoots
!                            are best handled by either linear interpolation or
!                            by the nonlinear inverse interpolation now used.
!     08/05/06   "     "     A simulation of the Shuttle in the LENS facility
!                            produced some bad profiles with Hwall = Hedge =
!                            Hinfinity, causing divides by zero for the film
!                            coefficient (now guarded against).
!     08/21/06-  "     "     Merged BLAYER2D/-3D as BLAYER as part of upgrading
!     09/08/06               to the Tecplot 360 version of Tecplot_io.f90.
!     04/22/07   "     "     High-density cases (256 points off the wall) showed
!                            the k window around the curvature-based peak needs
!                            to be larger for larger nk.  The root of this
!                            problem is spurious local curvature peaks, which
!                            are more likely with denser spacing.  A bigger
!                            window for the traditional method should find the
!                            proper peak enthalpy ratio for denser meshes.
!     04/25/07   "     "     The species density names have been off by one.
!                            For ns = 1 (air), we get wall and edge densities
!                            output twice (rhow and airw for DPLR input) but it
!                            is not worth the trouble suppressing airw & aire.
!     05/17/07   "     "     Switched to finite difference derivatives for the
!                            initial curvature calculations at the data points,
!                            because they are not affected by far away data.
!     06/26/07   "     "     Non-uniform in-flow (as from an arc-jet nozzle)
!                            means cutting off data below Hnorm = 0.95 is not
!                            viable.  Therefore, optional new inputs Hinf and
!                            Hignore are provided after the nx[tra] input.
!                            (Later: If Hignore > 0. is input, the second-step
!                            use of the traditional method - 99.5% of the peak
!                            in the edge neighborhood - is suppressed.  This
!                            helps deal with severe in-flow non-uniformity but
!                            overlooks the fact that curvature is not independ-
!                            ent of the normalization/data scaling.)
!     08/26/07   "     "     Permuted blocks were having the output array's nk
!                            set to the volume value instead of being left at 1.
!     01/22/08   "     "     Mars atmosphere cases with negative total enthalpy
!                            values prompted another optional Hshift input in
!                            the control file to avoid shifting via Tecplot.
!     02/01/08   "     "     Todd White needed to normalize each profile by its
!                            peak, so optional input Hinf now has 4 choices.
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
!     12/16/09   "     "     Arranged for Re-theta[/ Medge] output as an option
!                            via the roughness height control (-1. or -2.).
!     07/07/10   "     "     Dinesh noticed that Mars cases did not subtract the
!                            total enthalpy offset (needed for edge detection)
!                            from the interpolated enthalpies.  Sorry!
!     01/26/11   "     "     Jay Hyatt asked if the surface formed by the
!                            locations of the estimated boundary layer edge
!                            points along each radial line could be saved for
!                            visualization.  See the extended use of the
!                            output DATAPACKING control to invoke this option.
!     04/13/11   "     "     Trouble with a 2-D case from Jay Hyatt was traced
!                            to setting unit_normal(3), now set only if nd = 3.
!     09/10/12   "     "     Dinesh asked for velocity thickness as an option in
!                            place of displacement thickness.  Usage of the
!                            roughness height input k has been made even more
!                            involved as a result (but avoids a new input).
!                            See Note 8 above for the details.
!     02/12/13   "     "     This version replaces the Hshift fix (01/22/08)
!                            by working with (H - Hwall)/(Hinf - Hwall).
!                            Any Hshift input is ignored, to avoid affecting
!                            existing input control files.  But see Oct. 2017.
!     06/09/14   "     "     Belated merging of two changes: this version works
!                            with (H - Hw) / (Hinf - Hw) and the traditional
!                            method (if used) can save the indicated profile
!                            as originally implemented for the curvature-based
!                            method (only).
!                            Also, (x,y,z,f) output results for the body point
!                            specified in the control file are written to
!                            standard output.
!     06/11/14   "     "     Append output variable names to the list of results
!                            for the specified profile.  Also, the upper k limit
!                            for both edge methods is now 4*nk/5 rather than
!                            2*nk/3 (99.5% method) or 3*nk/4 (curvature method).
!                            Avoiding shock anomalies is the problem here.
!     08/28/14   "     "     CH (heat transfer coefficient) is now defined as
!                            qw/(Hinf - Hw), not qw/(Hedge - Hw), so as to be
!                            consistent with US3D developers.  The issue is
!                            that Hedge cannot be trusted for aft-body points.
!     03/18/15   "     "     Ross Chaudhry at U. of Minnesota questioned use of
!                            CSDVAL with spline coefficients b1, c1, d1 for
!                            which the associated CSFIT call had been commented
!                            out, forgetting this usage.  Somehow turning on
!                            all checks (-C with ifort) did not catch this,
!                            and other safeguards (including use of the 99.5%
!                            traditional method in the neighborhood of the
!                            curvature peak) masked the slip.
!                            Also:
!                "     "     Retrofitted the PLOT3D input file option that had
!                            originally been handled by BLAYER_RESULTS.
!                            Introduced blayer.inp.3 to allow entry of names
!                            for the species and any extras for PLOT3D cases.
!     10/05/17   "     "     Merged the two choices of total enthalpy profile
!                            via the edge_method control, q.v.
!     10/11/17   "     "     Dinesh didn't like the way delstar and theta are
!                            zero at a stagnation point (possibly only for
!                            axisymmetric solutions).  Therefore, trap zero
!                            delstar and, before results for that block are
!                            saved, extrapolate along the surface to the
!                            stag. point according to the recorded indices.
!                            Also: suppress duplicate profile output by the
!                            hybrid method (second stage was repeating what
!                            was already written by the first stage).
!     11/17/17  "       "    Raised tlimit in subroutine blayer_edge from 0.5
!                            to 0.7 on account of cases with nk = 101 (low) and
!                            thick boundary layers at high altitude. For the
!                            selected profile, a curvature-based preliminary
!                            Hedge was being updated by the traditional method
!                            before the save_plottable_details call.
!
!  Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA
!           Now with ERC, Inc. at NASA ARC (August 2010 through June 2015).
!           Now with AMA, Inc. at NASA ARC.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Modules:

   use grid_header_structure ! See Tecplot_io.f90 for all of these modules
   use grid_block_structure
   use tecplot_io_module
   use xyq_io_module         ! PLOT2D file I/O
   use xyzq_io_module        ! PLOT3D file I/O

   implicit none

!  Constants:

   integer, parameter :: &
      lun_in     = 1,    &   ! Input Tecplot dataset or PLOT3D volume grid
      lun_in2    = 2,    &   ! Input PLOT3D function file if present
      lun_out    = 3,    &   ! Normal output of surface + edge quantities
      lun_out2   = 4,    &   ! Optional 'blayer_edge_surface.dat' output
      lunctl     = 5,    &
      luncrt     = 6,    &
      name_limit = 32        ! Must match the value in Tecplot_io_module

   real, parameter ::    &
      stefan     = 5.66097E-8,  & ! Stefan-Boltzmann constant, w/(m^2 K^4)
      one        = 1.0,  &
      twenty     = 20.,  &
      zero       = 0.0

   character (1),  parameter :: &
      blank  = ' '
   character (11), parameter :: &
      format = 'unformatted'

   logical, parameter :: &
      false = .false.,   &
      new   = .true.,    &
      true  = .true.

!  Variables:

   integer :: &
      i, i1, ib, ibfirst, ios, iprofile, istag, j, jstag, jprofile, len, &
      mi, mj, mk, n, nblank, nblocks, nd, ni, nj, nk, normalize_choice, &
      nprofile, npts, nr, ns, nt, nx, numf_in, numf_out, numq_in

!  Indices set by set_bl_pointers local procedure (too many to be arguments):

   integer :: &
      id1, ip1, it1,it2, is1,is2, iu1,iu2, ih1, im1, iv1, ik1, ix1, ix2,       &
      jd1, jp1, jt1,jt2, jh1, jv1, js1,js2, jq1,jq2, ja1, ja2, jk1, jx1, jx2,  &
      kd1, kp1, kt1,kt2, kh1, ku1,ku2, km1, kv1, ks1,ks2, kdl, kds, kth, kue,  &
      kch, kk1, kx1, kx2, lhk, ldk, lvk, lmk, lkk

   integer, dimension (1) :: &
      iwall                  ! The intrinsic minloc needs an array for output

   integer, allocatable, dimension (:) :: &
      iblock, wall_override

   real :: &
      edge_method, emissivity, esigma, Hinf, Hignore, Hshift, Hwall, &
      roughness_height

   real, dimension (6) :: &
      average_ds

   logical :: &
      formatted_input, hybrid, k_zero, Re_kk, Re_theta, Re_theta_Medge, &
      save_edge_xyzs, shift_ratio, Tecplot_input, two_dim, &
      velocity_thickness, wall_ratio

   logical, allocatable, dimension (:) :: &
      suppress

   character (132) :: &
      buffer

!  Derived data types:

   type (grid_header) :: &
      header_in, header_out

   type (grid_type), pointer, dimension (:) :: &
      xyzq_in, xyzq_out

!  Tecplot function:

   integer :: TecEnd110

!  Execution:
!  !!!!!!!!!!

!  Read the controls from standard input:

   call read_controls ()  ! Local procedure below
   if (ios /= 0) go to 99 ! Single stop philosophy

!  Determine the number of dimensions and blocks, etc.:

   if (Tecplot_input) then

      header_in%ndim = -1  ! Tells the I/O package to determine if it's 2D or 3D
      ios            =  1  ! Verbose mode

      call Tec_header_read (lun_in, header_in, xyzq_in, ios)
      if (ios /= 0) then
         write (luncrt, '(a)') ' Trouble reading input file header.'
         go to 99
      end if

      nd      = header_in%ndim
      two_dim = nd == 2
      nblocks = header_in%nblocks
      numq_in = header_in%numq

   else  ! Open and read the headers of the input PLOT3D files, 2D or 3D

      len = len_trim (header_in%filename)
      i1  = 1;  if (formatted_input) i1 = 3

      open (lun_in, file=header_in%filename(1:len), status='old', &
            form=format(i1:), iostat=ios)
      if (ios /= 0) then
         write (luncrt, '(a)') ' Trouble opening the input grid file: ', &
            header_in%filename(1:len)
         go to 99
      end if

      call determine_grid_dim ('none', lun_in, formatted_input, nd, ios)

      two_dim = nd == 2  ! 'none' above means file is open, & rewound on return

      j = len - 1;  if (formatted_input) j = len
      header_in%filename(j:j) = 'f'

      open (lun_in2, file=header_in%filename(1:len), status='old', &
            form=format(i1:), iostat=ios)
      if (ios /= 0) then
         write (luncrt, '(a)') ' Trouble opening the input flow file: ', &
            header_in%filename(1:len)
         go to 99
      end if

!     Read the PLOT3D file headers:

      if (two_dim) then
         call xy_header_io (1, lun_in, formatted_input, nblocks, xyzq_in, ios)
         if (ios /= 0) then
            write (luncrt, '(/, a)') ' Trouble reading the input grid header.'
            go to 99
         end if
         call q_header_io_2d (1, lun_in2, formatted_input, nblocks, numq_in, &
                              xyzq_in, ios)
         if (ios /= 0) then
            write (luncrt, '(/, a)') ' Trouble reading the input flow header.'
            go to 99
         end if
         xyzq_in(:)%nk = 1
         xyzq_in(:)%mk = 1
      else
         call xyz_header_io (1, lun_in,  formatted_input, nblocks, xyzq_in, ios)
         if (ios /= 0) then
            write (luncrt, '(/, a)') ' Trouble reading the input grid header.'
            go to 99
         end if
         call q_header_io (1, lun_in2, formatted_input, nblocks, numq_in, &
                           xyzq_in, ios)
         if (ios /= 0) then
            write (luncrt, '(/, a)') ' Trouble reading the input flow header.'
            go to 99
         end if
      end if

      header_in%numq        = numq_in
      header_in%nblocks     = nblocks
      header_in%title       = blank
      header_in%ndatasetaux = 0
      xyzq_in(:)%nzoneaux   = 0
      xyzq_in(:)%zone_title = blank

      allocate (header_in%varname(nd+numq_in))

   end if  ! End of reading the volume dataset file header(s)

   header_out%ndim = 3  ! Because we insert 2D running length in the z position


!  Set up pointers into the input and output flow variables arrays based
!  on the apparent numbers of species, temperatures and extra quantities:
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   call set_bl_pointers (nd, ns, nt, nx, nr, numf_in, numf_out) ! Local proc.


!  Finish processing the input file header:
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   if (header_in%numq /= numf_in) then
      write (luncrt, '(/, (a, i3))')                                           &
         ' Unexpected number of input flow quantities:', header_in%numq,       &
         ' Expected number: ', numf_in
      go to 99
   end if

   nblocks = header_in%nblocks

   write (luncrt, '(/, a, /)') ' Input grid block dimensions:'
   if (two_dim) then
      write (luncrt, '(6(i6, 2X, 2i5))') &
         (ib, xyzq_in(ib)%ni, xyzq_in(ib)%nj, ib = 1, nblocks)
   else
      write (luncrt, '(4(i6, 2X, 3i5))') &
         (ib, xyzq_in(ib)%ni, xyzq_in(ib)%nj, xyzq_in(ib)%nk, ib = 1, nblocks)
   end if

   write (luncrt, '(/, a, i4)') &
      ' Number of input functions found:', header_in%numq


!  Look for the optional ancillary control file blayer.inp.2:
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   allocate (suppress(nblocks), iblock(nblocks), wall_override(nblocks))

   nblank = 0
   suppress(:) = false
   wall_override(:) = 0

   open (lunctl, file='blayer.inp.2', status='old', iostat=ios)

   if (ios == 0) then

!     Option to blank some blocks:

      nblank = nblocks  ! Upper limit for RDLIST

      write (luncrt, '(a)')

      call rdlist (luncrt, 'Reading blayer.inp.2 for blocks to suppress. ',    &
                   lunctl, nblank, iblock)

      do i = 1, nblank
         ib = iblock(i)
         if (ib < 1 .or. ib > nblocks) then
            write (luncrt, '(/, a, i9, a)') &
               ' Bad block number to suppress:', ib, '; proceeding.'
         else
            suppress(ib) = true
         end if
      end do

!     Option to override the automated detection of wall faces:

      do while (ios == 0) ! Until EOF
         read (lunctl, *, iostat=ios) ib, wall_override(ib)
      end do

      close (lunctl)
   end if

   if (nblank > 0) then
      write (luncrt, '(/, a, (20i4))') ' Blocks to blank:', iblock(1:nblank)
      do ib = 1, nblocks
         if (.not. suppress(ib)) then
            ibfirst = ib         ! First block not suppressed
            exit
         end if
      end do
   else
      write (luncrt, '(/, a)') ' No blocks are being blanked.'
      ibfirst = 1
   end if

   deallocate (iblock)

   write (luncrt, '(a)')
   buffer(1:4) = 'face';  if (two_dim) buffer(1:4) = 'edge'

   do ib = 1, nblocks
      i = wall_override(ib)
      if (i > 0) write (luncrt, '(a, i4, 3a, i2)') &
          ' For block', ib, ', wall ', buffer(1:4), ' number requested is', i
   end do

!  If PLOT3D input, blayer.inp.3 can be used to set the species + extras names:
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   if (.not. Tecplot_input) call set_input_variable_names ()

!  Set up the output file:
!  !!!!!!!!!!!!!!!!!!!!!!!

   header_out%numq    = numf_out
   header_out%nblocks = nblocks
   header_out%title   = header_in%title
   if (len_trim (header_out%title) == 0) header_out%title = header_in%filename

   allocate (header_out%varname(3 + numf_out))  ! 2-D run length replaces z

   call set_bl_names (nd, ns, nt, nx, nr, numf_in, numf_out, name_limit,       &
                      header_in%varname, header_out%varname)   ! Local procedure

   n = header_in%ndatasetaux
   header_out%ndatasetaux = n

   if (n > 0) then
      allocate (header_out%datasetauxname(n), header_out%datasetauxvalue(n))
      header_out%datasetauxname (:) = header_in%datasetauxname (:)
      header_out%datasetauxvalue(:) = header_in%datasetauxvalue(:)
   end if

   call Tec_header_write (lun_out, header_out, ios)
   if (ios /= 0) then
     write (luncrt, '(/, a)') ' Trouble writing the Tecplot output file header.'
     go to 99
   end if

!  Optional second output surface dataset (written one point at a time):

   if (save_edge_xyzs) then
      open  (lun_out2, file='blayer_edge_surface.dat', status='unknown')
      write (lun_out2, '(a)') &
         'TITLE = "Boundary layer edge surface"', &
         'VARIABLES =', '"x (m)"', '"y (m)"'
      if (.not. two_dim) write (lun_out2, '(a)') '"z (m)"'
      write (lun_out2, '(a)') '"delta (m)"'
   end if

   allocate (xyzq_out(nblocks))

!  The I/O package reads/writes 2-D data as x/y(i,j,1) data with z omitted.
!  We can treat 2-D inputs with (i,1,k) indexing and z set constant, but the
!  arrays proper don't need reordering.  Also, the output z is replaced by
!  running length for the 2-D case, which is therefore written as though it
!  were 3-D.  This differs from BLAYER2D, where running length was the first
!  function in the 2-D output.
!  We adjust the header info. after input and restore it before output.

   do ib = 1, nblocks
      xyzq_out(ib)%ni = xyzq_in(ib)%ni
      xyzq_out(ib)%nj = xyzq_in(ib)%nj
      if (two_dim) xyzq_out(ib)%nj = 1
      xyzq_out(ib)%nk = 1
      xyzq_out(ib)%mi = xyzq_in(ib)%mi
      xyzq_out(ib)%mj = xyzq_in(ib)%mj
      if (two_dim) xyzq_out(ib)%mj = 1
      xyzq_out(ib)%mk = 1
   end do

   esigma = emissivity * stefan


!  Process one block at a time:
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!

   do ib = 1, nblocks

      write (luncrt, '(a, i4)') ' Processing block', ib

      ni = xyzq_in(ib)%ni
      nj = xyzq_in(ib)%nj
      nk = xyzq_in(ib)%nk

!     Read an input block:

      if (Tecplot_input) then

         call Tec_block_allocate (xyzq_in(ib), nd, numf_in, ios)
         if (ios /= 0) then
            write (luncrt, '(/, a, i4)') 'Trouble allocating input block #', ib
            go to 99
         end if

         call Tec_block_read (lun_in, header_in, xyzq_in(ib), ios)
         if (ios /= 0) then
            write (luncrt, '(/, a, i4)') 'Trouble reading input block #', ib
            go to 99
         end if

      else ! PLOT3D-type inputs

         npts = ni*nj*nk

         if (two_dim) then

            call xy_allocate (xyzq_in(ib), ios)
            if (ios /= 0) then
               write (luncrt, '(/, a, i5)') &
                  ' Allocation trouble for PLOT2D grid input.  Block #:', ib
               go to 99
            end if

            call xy_block_io (1, lun_in, formatted_input, npts, &
                              xyzq_in(ib)%x, xyzq_in(ib)%y, ios)
            if (ios /= 0) then
               write (luncrt, '(/, a, i5)') &
                  ' Trouble reading PLOT2D grid.  Block #:', ib
               go to 99
            end if
            call q_allocate_2d (xyzq_in(ib), numq_in, ios)
            if (ios /= 0) then
               write (luncrt, '(/, a, i5)') &
                  ' Allocation trouble for PLOT2D function file.  Block #:', ib
               go to 99
            end if

            call q_block_io_2d (1, lun_in2, formatted_input, numq_in, &
                                xyzq_in(ib)%mi, xyzq_in(ib)%mj, &
                                xyzq_in(ib)%q, ios)
            if (ios /= 0) then
               write (luncrt, '(/, a, i5)') &
                  ' Trouble reading PLOT2D function file.  Block #:', ib
               go to 99
            end if

         else

            call xyz_allocate (xyzq_in(ib), ios)
            if (ios /= 0) then
               write (luncrt, '(/, a, i5)') &
                  ' Allocation trouble for PLOT3D grid input.  Block #:', ib
               go to 99
            end if

            call xyz_block_io (1, lun_in, formatted_input, npts,               &
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

            call q_block_io (1, lun_in2, formatted_input, numq_in,             &
                             xyzq_in(ib)%mi, xyzq_in(ib)%mj, xyzq_in(ib)%mk,   &
                             xyzq_in(ib)%q, ios)
            if (ios /= 0) then
               write (luncrt, '(/, a, i5)') &
                  ' Trouble reading PLOT3D function file.  Block #:', ib
               go to 99
            end if

         end if

      end if

!     Shuffle the indexing to make 2-space look like 3-space, and set z = 0?

      if (two_dim) then
         call switch_dimensions (1, xyzq_in(ib), numf_in)  ! 1 means 2D -> 3D
         nk = nj;  nj = 1
      end if

!     Option to deduce or to specify which block face is the wall:

      iwall(1) = wall_override(ib) ! May have been read from blayer.inp.2

      if (iwall(1) == 0) then

!        Identify the (likely) wall from the average off-face grid spacings:

         call average_increments (ni, nj, nk, xyzq_in(ib)%x, xyzq_in(ib)%y,    &
                                  xyzq_in(ib)%z, average_ds)

         iwall = minloc (average_ds)

      end if

      if (iwall(1) /= 5) then ! Swap faces [edges in 2D case] in-place:

         mi = ni;  mj = nj;  mk = nk

         call permute_block (1, numf_in, mi, mj, mk, iwall(1), xyzq_in(ib)%x,  &
                             xyzq_in(ib)%y, xyzq_in(ib)%z, xyzq_in(ib)%q,      &
                             ni, nj, nk)

         xyzq_in(ib)%ni = ni;  xyzq_out(ib)%ni = ni;  xyzq_out(ib)%mi = ni
         xyzq_in(ib)%nj = nj;  xyzq_out(ib)%nj = nj;  xyzq_out(ib)%mj = nj
         xyzq_in(ib)%nk = nk ! xyzq_out(ib)%nk = 1 ;  xyzq_out(ib)%mk = 1

!        if (ib == nprofile) then
!           It's not obvious what to do about iprofile and jprofile.
!           The user needs to enter them as they apply to the OUTPUT block.
!        end if

      end if

      xyzq_out(ib)%zone_title = xyzq_in(ib)%zone_title
      xyzq_out(ib)%nzoneaux   = xyzq_in(ib)%nzoneaux

!     Allocate the corresponding output block as 3-D even if it's 2-D input:

      call Tec_block_allocate (xyzq_out(ib), 3, numf_out, ios)
      if (ios /= 0) then
         write (luncrt, '(/, a, i4)') &
            ' Allocation trouble for output.  Block #:', ib
         go to 99
      end if

      if (xyzq_out(ib)%nzoneaux > 0) then
          xyzq_out(ib)%zoneauxname (:) = xyzq_in(ib)%zoneauxname (:)
          xyzq_out(ib)%zoneauxvalue(:) = xyzq_in(ib)%zoneauxvalue(:)
      end if


!     Process all the radial lines of this block, one at a time:
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     The input block may have been rearranged in-place, but the arrays were
!     not deallocated and reallocated.  The arrays are packed, so we can push
!     things down a level and employ the standard indexing.
!     With so many controls, it is most convenient to use a local procedure
!     rather than a bona fide subroutine.

      call process_block (ni, nj, nk, numf_in, xyzq_in(ib)%x, xyzq_in(ib)%y,   &
                          xyzq_in(ib)%z,  xyzq_in(ib)%q,  numf_out,            &
                          xyzq_out(ib)%x, xyzq_out(ib)%y, xyzq_out(ib)%z,      &
                          xyzq_out(ib)%q, ios)
!!!   if (ios /= 0) go to 99  ! Edge detection trouble, but keep going

      if (istag /= 0) call stag_pt_delstar_theta ()

      call deallocate_blocks (ib, ib, nd, numf_in, xyzq_in, ios)
      if (ios /= 0) go to 99

      if (two_dim) then
         deallocate (xyzq_in(ib)%z, stat=ios)
         if (ios /= 0) then
            write (luncrt, '(/, a, i4)') ' Trouble deallocating %z.  Block:', ib
            go to 99
         end if
      end if

!     Save results for this block:
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (ib == nprofile) then
         write (luncrt, '(/, a, 3i5, /)') &
            '#  Output results for ib, i, j =', ib, iprofile, jprofile
         write (luncrt, '(es20.8, 2x, a)') &
            xyzq_out(ib)%x(iprofile,jprofile,1), trim (header_out%varname(1)), &
            xyzq_out(ib)%y(iprofile,jprofile,1), trim (header_out%varname(2)), &
            xyzq_out(ib)%z(iprofile,jprofile,1), trim (header_out%varname(3))
         write (luncrt, '(i4, es16.8, 2x, a)') &
            (i, xyzq_out(ib)%q(i,iprofile,jprofile,1), &
             trim (header_out%varname(3+i)), i = 1, numf_out)
         write (luncrt, '(a)')
      end if

      call Tec_block_write (lun_out, header_out, xyzq_out(ib), ios)
      if (ios /= 0) then
         write (luncrt, '(/, a, i5)') &
            ' Trouble writing Tecplot output.  Block #:', ib
         go to 99
      end if

      call deallocate_blocks (ib, ib, 3, numf_out, xyzq_out, ios)
      if (ios /= 0) go to 99

   end do ! Next block to process


   if (header_out%formatted) then
      close (lun_out)
   else
      ios = TecEnd110 ()
   end if

   call deallocate_header (header_in)
   call deallocate_header (header_out)

   if (save_edge_xyzs) close (lun_out2)

99 continue

! *** stop ! Avoid system dependencies.


!  Internal procedures for program blayer, in the order they're used:

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine read_controls ()

!     Read BLAYER controls from standard input.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     All variables are inherited from the main program following modularization
!     of this piece of the code.

!     Execution:

!     The control file is echoed to standard output as we read it:

      ios = 0
      read  (lunctl, '(a)') buffer          ! Case description
      write (luncrt, '(/, 1x, a)') trim (buffer)
      read  (lunctl, '(a)') buffer          ! Text header for input file name
      write (luncrt, '(1x, a)') trim (buffer)
      read  (lunctl, *) buffer              ! Input volume dataset name
      len  =  len_trim (buffer)
      header_in%filename = buffer(1:len)
      write (luncrt, '(1x, a)') buffer(1:len)

      Tecplot_input = buffer(len-2:len) == 'dat'
      if (Tecplot_input) then
         formatted_input = true
      else
         Tecplot_input = buffer(len-2:len) == 'plt'  ! Not an option yet, though
         if (Tecplot_input) then
            formatted_input = false
            write (luncrt, '(/, a)') &
               ' Binary input Tecplot files (*.plt) are not yet an option.'
            ios = 1
            go to 99
         else  ! Must be a PLOT3D file, 2D or 3D
            formatted_input = buffer(len-1:len) /= 'gu'
         end if
      end if
      header_in%formatted = formatted_input

      read  (lunctl, '(a)') buffer          ! Text header for output file name
      write (luncrt, '(1x, a)') trim (buffer)
      read  (lunctl, *) header_out%filename
      len  =  len_trim (header_out%filename)
      write (luncrt, '(1x, a)') header_out%filename(1:len)
      header_out%formatted = header_out%filename(len-2:len) /= 'plt'

      read  (lunctl, *) n                   ! Datapacking
      save_edge_xyzs = n > 20               ! Kludge to avoid a new control
      if (save_edge_xyzs) n = n - 20
      write (luncrt, '(i2, 20x, a)') n,    '! 1 = POINT; 2 = BLOCK'
      header_out%datapacking = n - 1        ! Tecplot 360 usage

!     Optional miscellaneous controls:
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      edge_method = zero           ! Curvature-based method with shifted ratio
      emissivity = 0.89            ! Set defaults in case inputs are missing
      roughness_height = zero      ! Suppresses third set of outputs
      nprofile = 1                 ! Block and ...
      iprofile = 1;  jprofile = 1  ! ... i & j of sample profile output
      ns = 5;  nt = 1;  nx = 0     ! 5 species, 1 temperature, no extras
      nr = 0                       ! Unless roughness height > 0.
      normalize_choice = 0         ! Switch corresponding to Hinf
      Hinf    = zero               ! Use block 1 1,1,nk value unless input /= 0.
      Hignore = zero               ! Ignore Hnorm < 95% unless input > 0.
      Hshift  = zero               ! > 0. allows dealing with negative total
                                   ! enthalpies as for Mars atmosphere cases;
                                   ! this is ignored if the total enthalpy ratio
                                   ! choice is (H - Hwall) / (Hinf - Hwall)

      read  (lunctl, '(a)', iostat=ios) buffer  ! Text header for misc. controls

      if (ios /= 0) then
         write (luncrt, '(/, (a))') ' Miscellaneous controls missing?',        &
            ' Proceeding with 5 species, 1 temperature, 0 extras.'
      else
         write (luncrt, '(1x, a)') trim (buffer)
         read  (lunctl, *, iostat=ios) edge_method
         if (ios /= 0) then
            write (luncrt, '(/, (a))') ' Edge method, etc. missing?',          &
               ' Proceeding with hybrid curvature-based method.'
         else
            read  (lunctl, *, iostat=ios) emissivity
            if (ios /= 0) then
               write (luncrt, '(/, (a))') ' Emissivity, etc. missing?',        &
                  ' Proceeding with 0.89.'
            else
               read  (lunctl, *, iostat=ios) roughness_height
               if (ios /= 0) then
                  write (luncrt, '(/, (a))') ' Roughness height missing?',     &
                     ' Proceeding with k = 0.'
               else
                  read (lunctl, *, iostat=ios) nprofile, iprofile, jprofile
                  if (ios /= 0) then
                    write (luncrt, '(/, a)') ' Sample profile specs. missing?'
                    write (luncrt, '(a)') &
                       ' Proceeding with block 1, (i,j) = (1,1) profile.'
                  else
                     read (lunctl, *, iostat=ios) ns;  if (ios /= 0) ns = 5
                     read (lunctl, *, iostat=ios) nt;  if (ios /= 0) nt = 1
                     read (lunctl, *, iostat=ios) nx;  if (ios /= 0) nx = 0
                     if (ios /= 0) then
                        write (luncrt, '(/, a, /, a, 3i4)')                    &
                        ' # species, # temperatures and/or # extras missing?', &
                          ' Proceeding with these numbers:', ns, nt, nx
                     end if
                     read (lunctl, *, iostat=ios) Hinf
                     if (ios /= 0) Hinf = zero
                     read (lunctl, *, iostat=ios) Hignore
                     if (ios /= 0) Hignore = zero
                     read (lunctl, *, iostat=ios) Hshift
                     if (ios /= 0) Hshift = zero
                  end if
               end if
            end if
         end if
      end if

      close (lunctl);  ios = 0

!     Extended use of edge method input:

      shift_ratio = edge_method >= zero     ! Use (H + Hshift) / (Hinf + Hshift)
      wall_ratio  = edge_method <  zero     ! Use (H - Hwall)  / (Hinf - Hwall)
      hybrid      = abs (edge_method) < 90. ! Use curvature-based hybrid method

!     Multiple use of the roughness height input:

      velocity_thickness = false
      if (roughness_height  >= twenty) then
          roughness_height   = roughness_height - twenty
          velocity_thickness = true  ! In place of displacement thickness
      else if (roughness_height <= -twenty) then
          roughness_height   = roughness_height + twenty
          velocity_thickness = true
      end if

      k_zero         = roughness_height ==  zero;  if (.not. k_zero) nr = 5
      Re_kk          = roughness_height  >  zero
      Re_theta       = roughness_height ==  -one   ! Not used, but clearer
      Re_theta_Medge = roughness_height == -(one + one)

      write (luncrt, '(f7.3, 15x, a)') edge_method,        '! Edge method'
      write (luncrt, '(f7.3, 15x, a)') emissivity,         '! Emissivity'
      write (luncrt, '(f9.6, 13x, a)') roughness_height,   '! Roughness_height'
      write (luncrt, '(3i4,  10x, a)') nprofile, iprofile, jprofile, &
                                           '! Profile block & (i,j)'
      write (luncrt, '(i2,  20x,  a)') ns, '! # species'
      write (luncrt, '(i2,  20x,  a)') nt, '! # temperatures'
      write (luncrt, '(i2,  20x,  a)') nx, '! # extra items'
      write (luncrt, '(a)')

      if (Hinf > zero) then
         normalize_choice = 1
         write (luncrt, '(es16.7, 6x, a)') Hinf, &
            '! Input total enthalpy to treat as Hinfinity.'
      else if (Hinf == zero) then
         normalize_choice = 0
         write (luncrt, '(2a)') &
          ' ! Using block 1, (1,1,nk) value of total enthalpy for all profiles.'
      else if (Hinf == -one) then
         normalize_choice = -1
         write (luncrt, '(2a)') &
            ' ! Using outermost value of total enthalpy for each profile.'
      else if (Hinf == -2.) then
         normalize_choice = -2
         write (luncrt, '(2a)') &
            ' ! Using peak value of total enthalpy for each profile.'
      else
         write (luncrt, '(a, es16.7)') ' ! Invalid entry for Hinf:', Hinf
         ios = 1
      end if

      if (Hignore > zero) then
         if (Hignore < one) Hignore = Hignore * 100.  ! Fraction -> percentage
      else
         Hignore = 95.
      end if
      write (luncrt, '(f6.2, 16x, a, /)') Hignore, &
         '! Total enthalpy ratio cut-off: b.l. edge cannot be below this %'
      if (Hignore > one) Hignore = Hignore * 0.01     ! Percentage -> fraction

   99 return

      end subroutine read_controls

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine switch_dimensions (mode, block, nf)

!     For 3-space applications where i and j are the surface directions and k is
!     in the outward off-surface direction, it is convenient to treat a 2-space
!     example as though it were dimensioned (ni,1,nk) rather than (ni,nj,1).
!
!     If mode = 1, this utility adjusts one block after it is read, to allow
!     3-space-type indexing.  There is no need to move the data proper when the
!     arrays are packed, but %z hasn't been initialized by the read utility so
!     allocate it and set it constant here.
!
!     If mode = 2, the readjustment necessary for 2D output is made to the
!     block (which may not be the input block any more, so we can't deallocate
!     any %z here).
!
!     08/31/06  David Saunders  Initial implementation.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer,          intent (in)    :: mode   ! 1 or 2 as described above
      type (grid_type), intent (inout) :: block  ! One grid block, %q optional
      integer,          intent (in)    :: nf     ! # functions >= 0

!     Execution:

      select case (mode)

         case (1) ! 2D -> 3D

            block%nk = block%nj;  block%nj = 1

            if (nf > 0) then
               block%mk = block%mj;  block%mj = 1
            end if

            allocate (block%z(block%ni,1,block%nk))

            block%z(:,:,:) = 0.

         case (2) ! 3D -> 2D

            block%nj = block%nk;  block%nk = 1  ! nk isn't used for 2D I/O

            if (nf > 0) then
               block%mj = block%mk;  block%mk = 1
            end if

      end select

      end subroutine switch_dimensions

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine set_input_variable_names ()
!
!     For PLOT3D input mode (only), set the variable names for the species and
!     any extra variables by reading file blayer.inp.3.  E.g.:
!
!        n2 o2 no no+ n2+ o2+ n o n+ o+ e      [all species on line 1]
!        C_p N_tot                             [all  extras on line 2]
!
!     If either line is empty, the variable names are defaulted like this:
!
!        sp_1 sp_2 sp_3 ..... sp_11
!        xtra_1 xtra_2
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: i, l, nextras, nspecies
      logical :: default_species, default_extras

      default_extras = nx > 0

      open (lunctl, file='blayer.inp.3', status='old', iostat=ios)

      default_species = ios /= 0

      if (.not. default_species) then  ! Look for species names on line 1

         read (lunctl, '(a)', iostat=ios) buffer
         default_species = ios /= 0
         if (.not. default_species) then
            l = len_trim (buffer)
            default_species = l == 0
            if (.not. default_species) then  ! Tokenize without changing case
               nspecies = ns
               call token2l (buffer(1:l), ' ,', nspecies, &
                             header_in%varname(nd+is1:nd+is2))
               default_species = nspecies /= ns

!              Similarly for any extras names on line 2:

               if (nx > 0) then
                  read (lunctl, '(a)', iostat=ios) buffer
                  default_extras = ios /= 0
                  if (.not. default_extras) then
                     l = len_trim (buffer)
                     default_extras = l == 0
                     if (.not. default_extras) then
                        nextras = nx
                        call token2l (buffer(1:l), ' ,', nextras, &
                                      header_in%varname(nd+ix1:nd+ix2))
                        default_extras = nextras /= nx
                     end if
                  end if
               end if

            end if
         end if
      end if

      close (lunctl)

      if (default_species) then
         do i = 1, ns
            call numbered_name ('sp_', i, header_in%varname(nd+is1+i-1), l)
         end do
      end if

      if (default_extras) then
         do i = 1, nx
            call numbered_name ('xtra_', i, header_in%varname(nd+ix1+i-1), l)
         end do
      end if

      do i = 1, nblocks
         call numbered_name ('zone ', i, xyzq_in(i)%zone_title, l)
      end do

      end subroutine set_input_variable_names

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine set_bl_pointers (nd, ns, nt, nx, nr, nq_in, nq_out)
!
!     Local procedure for BLAYER, to set input and output variable pointers into
!     flow field %q arrays according to the numbers of dimensions, species,
!     temperatures, and miscellaneous extras.
!
!     09/01/06  DAS  When 2-D & 3-D were merged, it was better to treat the 2-D
!                    running length as a substitute for z than as an extra
!                    variable not present in the 3-D case.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in)  :: nd, ns, nt, nx  ! See description above
      integer, intent (in)  :: nr              ! # roughness height values, 0|5
      integer, intent (out) :: nq_in, nq_out   ! Implied numbers of flow field
                                               ! quantities excluding x, y, z
!     Execution:

!     Input quantity pointers, i**:

      id1 = 1           ! density
      ip1 = id1 + 1     ! pressure
      it1 = ip1 + 1     ! nt temperatures
      it2 = ip1 + nt
      is1 = it2 + 1     ! ns species densities
      is2 = it2 + ns
      iu1 = is2 + 1     ! nd velocity components
      iu2 = is2 + nd
      ih1 = iu2 + 1     ! total enthalpy
      im1 = ih1 + 1     ! Mach number
      iv1 = im1 + 1     ! viscosity

      if (nt == 1) then
         ik1 = 0        ! suppress kappa
         ix1 = iv1 + 1  ! extra variables
         ix2 = iv1 + nx
      else
         ik1 = iv1 + 1  ! kappa (thermal conductivity)
         ix1 = ik1 + 1  ! extras
         ix2 = ik1 + nx
      end if

      nq_in = 5 + nd + nt + ns + nx ! Should agree with input file

      if (nt > 1) nq_in = nq_in + 1 ! For thermal conductivity, kappa

!     Output pointers, at wall & at boundary layer edge, plus integrated values:

!     First, wall quantity pointers, j** (pointer into wall_names(:)):

      jd1 = 1           ! density
      jp1 = jd1 + 1     ! pressure
      jt1 = jp1 + 1     ! nt temperatures
      jt2 = jp1 + nt
      jh1 = jt2 + 1     ! total enthalpy
      jv1 = jh1 + 1     ! viscosity
      js1 = jv1 + 1     ! ns species densities
      js2 = jv1 + ns
      jq1 = js2 + 1     ! nt heat fluxes
      jq2 = js2 + nt
      ja1 = jq2 + 1     ! shear stresses
      ja2 = jq2 + nd

      if (nt == 1) then
         jk1 = 0        ! suppress thermal conductivity
         jx1 = ja2 + 1  ! extra variables
         jx2 = ja2 + nx
      else
         jk1 = ja2 + 1  ! kappa
         jx1 = jk1 + 1  ! extras
         jx2 = jk1 + nx
      end if

!     B.L. edge quantity pointers, k**:

      kd1 = jx2 + 1     ! density
      kp1 = kd1 + 1     ! pressure
      kt1 = kp1 + 1     ! nt temperatures
      kt2 = kp1 + nt
      kh1 = kt2 + 1     ! total enthalpy
      ku1 = kh1 + 1     ! nd velocity components
      ku2 = kh1 + nd
      km1 = ku2 + 1     ! Mach number
      kv1 = km1 + 1     ! viscosity
      ks1 = kv1 + 1     ! ns species densities
      ks2 = kv1 + ns
      kdl = ks2 + 1     ! delta (B.L. thickness)
      kds = kdl + 1     ! delta* (displacement thickness) or vel. thickness
      kth = kds + 1     ! theta   (momentum thickness)
      kue = kth + 1     ! unit Reynolds number
      kch = kue + 1     ! film coefficient, CH

      if (nt == 1) then
         kk1 = 0        ! suppress thermal conductivity
         kx1 = kch + 1  ! extra variables
         kx2 = kch + nx
      else
         kk1 = kch + 1  ! kappa
         kx1 = kk1 + 1  ! extras
         kx2 = kk1 + nx
      end if

!     Roughness height quantity pointers, l**, unless suppressed:

      if (nr > 0) then
         lhk = kx2 + 1  ! height k or theta
         ldk = lhk + 1  ! density   at height k or theta
         lvk = ldk + 1  ! velocity  at height k or theta
         lmk = lvk + 1  ! viscosity at height k or theta
         lkk = lmk + 1  ! Re-kk or Re-theta or Re-theta/Medge
      end if

!     What else???

!     Count the output quantities:

      nq_out = 14 + 2*nd + 2*ns + 3*nt + 2*nx

      if (nt > 1) nq_out = nq_out + 2       ! thermal conductivities
      if (nr > 0) nq_out = nq_out + 5       ! roughness height quantities

!!!   write (luncrt, '(/, a, /, a, /, a, /, a, /, (15i4))') &
!!! ' id1, ip1, it1,it2, is1,is2, iu1,iu2, ih1, im1, iv1, ik1, ix1, ix2,    ', &
!!! ' jd1, jp1, jt1,jt2, jh1, jv1, js1,js2, jq1,jq2, ja1, ja2, jk1, jx1, jx2', &
!!! ' kd1, kp1, kt1,kt2, kh1, ku1,ku2, km1, kv1, ks1,ks2, kdl, kds, kth, kue', &
!!! ' kch, kk1, kx1, kx2, lhk, ldk, lvk, lmk, lkk', &
!!!   id1, ip1, it1,it2, is1,is2, iu1,iu2, ih1, im1, iv1, ik1, ix1, ix2,       &
!!!   jd1, jp1, jt1,jt2, jh1, jv1, js1,js2, jq1,jq2, ja1, ja2, jk1, jx1, jx2,  &
!!!   kd1, kp1, kt1,kt2, kh1, ku1,ku2, km1, kv1, ks1,ks2, kdl, kds, kth, kue,  &
!!!   kch, kk1, kx1, kx2, lhk, ldk, lvk, lmk, lkk

      end subroutine set_bl_pointers

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine set_bl_names (nd, ns, nt, nx, nr, nq_in, nq_out,              &
                               name_limit, names_in, names_out)
!
!     Local procedure for BLAYER2D/3D, to set output variable names following
!     use of set_bl_pointers to determine their number, etc., from the given
!     number of dimensions, species, temperatures, and miscellaneous extras.
!
!     Names for species densities and miscellaneous extra variables must be
!     transferred from the input data, since their order is not fixed.  Other
!     variables should be in a definite order, and can be assigned from the data
!     statements established here.
!
!     Note that Tecplot's variable name arrays include x, y, z.  The (1-nd:nq)
!     declarations keep the pointers into the %q and name arrays the same here.
!     The calling program should just pass names_in and names_out.
!
!     09/01/06  DAS  When 2-D & 3-D were merged, it was better to treat the 2-D
!                    running length as a substitute for z than as an extra
!                    variable not present in the 3-D case.  Also, the input
!                    names are arranged to include 'z' for the 2-D case at the
!                    higher level (outside the I/O package).
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer,   intent (in) ::  nd, ns, nt, nx     ! See description above
      integer,   intent (in) ::  nr                 ! # roughness ht. vals., 0|5
      integer,   intent (in) ::  nq_in, nq_out      ! Derived from nd/ns/nt/nx
      integer,   intent (in) ::  name_limit         ! Maximum length of a name
      character (name_limit), intent (in),  dimension (1-nd:nq_in) :: names_in
      character (name_limit), intent (out), dimension (-2:nq_out)  :: names_out

!     Local constants:

      integer, parameter   :: local_limit  = 18     ! Enough for running length
      character (1), parameter :: wall_tag = 'w',   &
                                  edge_tag = 'e',   &
                                  ruff_tag = 'r'
      character (18), parameter :: RethetaMe = 'Re_theta/Medge    '

!     Local variables:

      integer   :: i, j
      character (local_limit) :: wall_names(15)
      character (local_limit) :: edge_names(16)
      character (local_limit) :: ruff_names(5)
      character (local_limit) :: theta_names(5)

!     Storage:

      data wall_names /        &
         'xw (m)            ', & ! 1
         'yw (m)            ', & ! 2
         'zw (m)            ', & ! 3  ! Replaced by running length for 2D
         'rhow (kg/m^3)     ', & ! 4
         'pw (Pa)           ', & ! 5
         'Tw (K)            ', & ! 6
         'Tvw (K)           ', & ! 7
         'Hw (J/kg)         ', & ! 8
         'muw (Pa.s)        ', & ! 9
         'qw (W/m^2)        ', & ! 10
         'qvw (W/m^2)       ', & ! 11
         'tauwx (Pa)        ', & ! 12
         'tauwy (Pa)        ', & ! 13
         'tauwz (Pa)        ', & ! 14
         'kappaw (W/m.K)    '/   ! 15

      data edge_names /        &
         'rhoe (kg/m^3)     ', & ! 1
         'pe (Pa)           ', & ! 2
         'Te (K)            ', & ! 3
         'Tve (K)           ', & ! 4
         'He (J/kg)         ', & ! 5
         'ue (m/s)          ', & ! 6
         've (m/s)          ', & ! 7
         'we (m/s)          ', & ! 8
         'Me                ', & ! 9
         'mue (Pa.s)        ', & ! 10
         'delta (m)         ', & ! 11
         'deltastar (m)     ', & ! 12
         'theta (m)         ', & ! 13
         'Re-ue             ', & ! 14
         'CH (kg/m^2.s)     ', & ! 15
         'kappae (W/m.K)    '/   ! 16

      data ruff_names /        &
         'roughness (m)     ', & ! 1
         'rhok (kg/m^3)     ', & ! 2
         'velk (m/s)        ', & ! 3
         'muk (Pa.s)        ', & ! 4
         'Re-kk             '/   ! 5

      data theta_names /       &
         'theta (t) (m)     ', & ! 1
         'rhot (kg/m^3)     ', & ! 2
         'velt (m/s)        ', & ! 3
         'mut (Pa.s)        ', & ! 4
         'Re-theta          '/   ! 5

!     Execution:

      names_out(-2:0)    = wall_names(1:3)        ! x, y, z
      if (nd == 2) names_out(0:0) = 'running length (m)'
      names_out(jd1)     = wall_names(4)          ! density
      names_out(jp1)     = wall_names(5)          ! pressure
      names_out(jt1:jt2) = wall_names(6:5+nt)     ! T[, Tv]
      names_out(jh1)     = wall_names(8)          ! total enthalpy
      names_out(jv1)     = wall_names(9)          ! viscosity
      names_out(js1:js2) = names_in(is1:is2)      ! species densities

      call append_tag (js1, js2, name_limit, names_out(1), wall_tag)

      names_out(jq1:jq2) = wall_names(10:9+nt)    ! heat flux[es]
      names_out(ja1:ja2) = wall_names(12:11+nd)   ! shear stresses

      if (nt > 1) names_out(jk1) = wall_names(15) ! kappa (thermal conductivity)

      if (nx > 0) then                            ! Optional extras
        names_out(jx1:jx2) = names_in(ix1:ix2)
        call append_tag (jx1, jx2, name_limit, names_out(1), wall_tag)
      end if

!     Boundary layer edge quantities (k pointers):

      names_out(kd1)     = edge_names(1)          ! density
      names_out(kp1)     = edge_names(2)          ! pressure
      names_out(kt1:kt2) = edge_names(3:2+nt)     ! T[, Tv]
      names_out(kh1)     = edge_names(5)          ! total enthalpy
      names_out(ku1:kh1+nd) = edge_names(6:5+nd)  ! u, v[, w]
      names_out(km1)     = edge_names(9)          ! Mach number
      names_out(kv1)     = edge_names(10)         ! viscosity
      names_out(ks1:ks2) = names_in(is1:is2)      ! species densities

      call append_tag (ks1, ks2, name_limit, names_out(1), edge_tag)

      names_out(kdl)     = edge_names(11)         ! delta
      names_out(kds)     = edge_names(12)         ! delta-*
      if (velocity_thickness) names_out(kds) = 'vel_thickness (m) '  ! See k use
      names_out(kth)     = edge_names(13)         ! theta
      names_out(kue)     = edge_names(14)         ! Re-ue
      names_out(kch)     = edge_names(15)         ! CH

      if (nt > 1) names_out(kk1) = edge_names(16) ! kappa

      if (nx > 0) then                            ! Optional extras
        names_out(kx1:kx2) = names_in(ix1:ix2)
        call append_tag (kx1, kx2, name_limit, names_out(1), edge_tag)
      end if

!     Roughness height quantities (l pointers) if any:

      if (nr > 0) then
         if (Re_kk) then
            names_out(lhk)  = ruff_names(1)       ! height k
            names_out(ldk)  = ruff_names(2)       ! density
            names_out(lvk)  = ruff_names(3)       ! velocity
            names_out(lmk)  = ruff_names(4)       ! viscosity
            names_out(lkk)  = ruff_names(5)       ! Re-kk
         else
            names_out(lhk)  = theta_names(1)      ! theta t
            names_out(ldk)  = theta_names(2)      ! density   at theta
            names_out(lvk)  = theta_names(3)      ! velocity  at theta
            names_out(lmk)  = theta_names(4)      ! viscosity at theta
            names_out(lkk)  = theta_names(5)      ! Re-theta
            if (Re_theta_Medge) names_out(lkk) = RethetaMe  ! Re-theta/Medge
         end if
      end if

      end subroutine set_bl_names

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine append_tag (i1, i2, name_limit, names, tag)
!
!     This is really local to set_bl_names, but local procedures cannot have
!     local procedures of their own.
!
!     Append 'w', 'e', or 'r' to the given name to keep all names unique.
!     If the input name includes units, it's not clear how to insert instead.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer,   intent (in)    :: i1, i2                   ! Names to adjust
      integer,   intent (in)    :: name_limit               ! Length limit
      character (name_limit), intent (inout) :: names (i2)  ! With room for tag
      character (1),          intent (in)    :: tag         ! See description

!     Local variables:

      integer :: i, l

!     Execution:

      do i = i1, i2
         l = min (len_trim (names(i)) + 1, name_limit)
         names(i)(l:l) = tag
      end do

      end subroutine append_tag

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

      real, intent (inout), dimension (nf,ni,nj,nk) :: &
         f                            ! Inout as total enthalpies may be shifted

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

      character (1), parameter :: &
         method = 'B'         ! "Bessel" option: plain 4-pt. local spline method

!     Local variables:

      integer :: &
         i, in, j, k, k1, ka, kb, kedge, l, n, nbad, nk_height_k, numk

      integer, dimension (2 + nd) :: &
         k_indices                    ! For density, velocity components, and mu

      real :: &
         delstar, dsq, emissivity, Hedge, height, Hinf_use, Hmin, p, q, r,     &
         rho_vtan_edge, redge, tedge, theta, total, Twall, V, vedge,           &
         vtangent, vtangent_edge, wall_mu, xedge, yedge, yeval, ypeval, zedge

      real, dimension (nd) :: &
         unit_normal, v_tangent, v_total

      real, dimension (5) :: &
         height_k_values

      real, allocatable, dimension (:) :: &
         Hratio, rratio, vratio, fline, tline, xline, yline, zline

      logical :: &
         save_profile

      character (32) :: &
         subtitle

!     Execution:

      Hedge = abs (edge_method) * 0.01  ! For the case when % Hinf is specified

!     Items to be interpolated at roughness height k along with wall viscosity:

      k_indices(1) = id1              ! Density
      do i = 1, nd                    ! Velocity components
         k_indices(1+i) = iu1 + i - 1
      end do
      k_indices(2+nd) = iv1           ! Wall viscosity

      nk_height_k = (2 * nk) / 3      ! Enough to contain roughness height k
      ka = 8                          ! Index range for traditional method
      kb = (4 * nk) / 5               ! Make it consistent with curvature method

      allocate (xline(nk), yline(nk), zline(nk), tline(nk), fline(nk),         &
                Hratio(nk), rratio(nk), vratio(nk))

!     Shifted total enthalpy ratio?

      if (.not. suppress(ib)) then
         if (shift_ratio) then ! Check that shifted enthalpies are all positive
            Hmin = 1.e+30
            do k = 1, nk
               do j = 1, nj
                  do i = 1, ni
                     Hmin =    min (f(ih1,i,j,k) + Hshift,  Hmin)
                  end do
               end do
            end do
            if (Hmin < zero) write (luncrt, '(/, a, es25.15, a, i5, /)') &
               ' *** Warning:  Minimum shifted total enthalpy:', Hmin,   &
               ';  block #:', ib
         end if
      end if

      if (normalize_choice == 0) then  ! Get Hinf from first block?
         if (Hinf == zero) then
            if (ib == ibfirst) then ! First unsuppressed block
               Hinf = f(ih1,1,1,nk) ! Free-strm. total enthalpy from 1st profile
               Hinf_use = Hinf
               write (luncrt, '(/, a, i3, a, es18.10, /)') &
               ' Free-stream total enthalpy to be used from block', ib, ':',   &
                  Hinf
            end if
         end if
      end if

      Hinf_use = Hinf  ! May be overidden if normalize_choice < 0
      height   = roughness_height  ! May be changed to theta
      nbad     = 0  ! Suppress output of problematic profiles if too many
      istag    = 0  ! A zero delstar will be replaced with an extrapolation;
      jstag    = 0  ! likewise for theta (probably exact zero for axisym. only)

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
               write (subtitle, '(a6, i4, 2(5x, a2, i4))')                     &
                  'Block:', nprofile, 'i:', iprofile, 'j:', jprofile
            end if

!           Wall values:

            fout(jd1,i,j)     = f(id1,i,j,1)           ! rho
            fout(jp1,i,j)     = f(ip1,i,j,1)           ! p
            fout(jt1:jt2,i,j) = f(it1:it2,i,j,1)       ! Temp(s)
            Twall             = f(it1,i,j,1)
            fout(jh1,i,j)     = f(ih1,i,j,1)           ! total enthalpy
            fout(jv1,i,j)     = f(iv1,i,j,1)           ! viscosity
            fout(js1:js2,i,j) = f(is1:is2,i,j,1)       ! Species
            fout(jq1,i,j)     = esigma * Twall**4      ! Total heat flux

            select case (normalize_choice)

               case (0, 1) ! First profile or user-supplied value, already set
!!!               Hinf_use = Hinf
               case (-1)   ! Outermost value from each profile
                  Hinf_use = f(ih1,i,j,nk)
               case (-2)   ! Peak value from each profile
                  Hinf_use = -1.e+30
                  do k = 1, nk
                      Hinf_use = max (Hinf_use, f(ih1,i,j,k))
                  end do

            end select

!           Set up this radial line for estimating the boundary layer edge.
!           We need arc length data for q(catalytic) & for shear terms also.

            do k = 1, nk
               xline(k)  = x(i,j,k)
               yline(k)  = y(i,j,k)
               zline(k)  = z(i,j,k)
            end do

            if (shift_ratio) then
               do k = 1, nk
                  Hratio(k) = (f(ih1,i,j,k) + Hshift) / (Hinf_use + Hshift)
               end do
            else
               Hwall = f(ih1,i,j,1)
               do k = 1, nk
                  Hratio(k) = (f(ih1,i,j,k) - Hwall) / (Hinf_use - Hwall)
               end do
            end if

!           Unnormalized arc lengths:

            call chords3d (nk, xline, yline, zline, false, total, tline)

!           Catalytic heat flux?

            if (nt > 1) then !  q(catalytic) = q(total) - kappa(wall) x dT/dt

               fout(jq2,i,j)  = fout(jq1,i,j) - f(ik1,i,j,1) *                 &
                                (f(it1,i,j,2) - Twall) / tline(2)
            end if

!           Shear stress terms = normal derivatives of tangential velocity:

            unit_normal(1)  = (xline(2) - xline(1)) / tline(2)
            unit_normal(2)  = (yline(2) - yline(1)) / tline(2)
            if (nd == 3) &
            unit_normal(3)  = (zline(2) - zline(1)) / tline(2)

            v_total(:)      = f(iu1:iu2,i,j,2)
            v_tangent(:)    = v_total(:) - dot_product (v_total, unit_normal)  &
                            * unit_normal(:)
            fout(ja1:ja2,i,j) = v_tangent(:) * &                   ! (vt2 - 0) *
                                (f(iv1,i,j,1) / tline(2))          ! (mu / dt)

            if (nt > 1) fout(jk1,i,j)     = f(ik1,i,j,1)  ! Thermal conductivity

            if (nx > 0) fout(jx1:jx2,i,j) = f(ix1:ix2,i,j,1)       ! Extras


!           Estimate the edge of the boundary layer, tedge:
!           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            if (hybrid) then ! Two-stage curvature-based method

               call blayer_edge (hybrid, shift_ratio, nk, Hratio, tline, &
                                 Hignore, save_profile, subtitle, kedge, &
                                 Hedge, tedge, ier)
            else

               call bl_edge_traditional (shift_ratio, Hedge, nk, ka, kb, &
                                         Hratio, tline, save_profile, kedge, &
                                         tedge, ier)
            end if

            if (ier /= 0) then
               nbad = nbad + 1
               if (nbad == 1) then
                  write (luncrt, '(/, a, 3i5)') &
                     ' Trouble detecting boundary layer edge. ib, i, j:', ib,i,j
                  write (luncrt, '(/, 2a)') &
                  '#                 x                  y                  z', &
                  '                  t     Enthalpy Ratio'
                  write (luncrt, '(5es19.11)') &
                  (xline(k), yline(k), zline(k), tline(k), Hratio(k), k = 1, nk)
               end if
            end if ! Keep going anyway

!           The local spline method used requires just 4 points:

            k1 = kedge - 1 ! tedge is in the interval [t(kedge), t(kedge+1))


            if (save_edge_xyzs) then  ! Option to visualize b.l. edge surface

               call lcsfit (4, tline(k1), xline(k1), new, method, 1,           &
                            tedge, xedge, ypeval)
               call lcsfit (4, tline(k1), yline(k1), new, method, 1,           &
                            tedge, yedge, ypeval)
               call lcsfit (4, tline(k1), zline(k1), new, method, 1,           &
                            tedge, zedge, ypeval)

               if (i + j == 2) then  ! Zone header
                  write (lun_out2, '(a)') 'ZONE T="bl_edge"'
                  if (two_dim) then
                     write (lun_out2, '(a, i4, a, i4, a)') &
                        ' I=', ni, ', J=', nj, ', ZONETYPE=Ordered'
                  else
                     write (lun_out2, '(a, i4, a, i4, a)') &
                        ' I=', ni, ', J=', nj, ', K= 1, ZONETYPE=Ordered'
                  end if
                  write (lun_out2, '(a)') ' DATAPACKING=POINT'
               end if

               if (two_dim) then
                  write (lun_out2, '(3es16.8)') xedge, yedge, tedge
               else
                  write (lun_out2, '(4es16.8)') xedge, yedge, zedge, tedge
               end if

            end if

            xout(i,j) = xline(1)
            yout(i,j) = yline(1)
            zout(i,j) = zline(1)


!           Interpolate output {rho, p, ..., species densities} at the
!           boundary layer edge.
!           We can't use variables in case (*:*), so it's awkward!

            do n = kd1, ks2

               if (n == kd1) in = id1                           ! rho
               if (n == kp1) in = ip1                           ! p
               if (n >= kt1 .and. n <= kt2) in = it1 + n - kt1  ! Temp(s)
               if (n == kh1) in = ih1                           ! Total enthalpy
               if (n >= ku1 .and. n <= ku2) in = iu1 + n - ku1  ! u or v
               if (n == km1) in = im1                           ! Mach
               if (n == kv1) in = iv1                           ! viscosity
               if (n >= ks1 .and. n <= ks2) in = is1 + n - ks1  ! spc. densities

               do k = k1, k1 + 3
                  yline(k) = f(in,i,j,k)
               end do

               call lcsfit (4, tline(k1), yline(k1), new, method, 1, tedge,    &
                            yeval, ypeval)

!!             if (n == kh1) yeval = yeval - Hshift ! Restore any enthalpy shift
                                                    ! No - no need now
               fout(n,i,j) = yeval

            end do

            fout(kdl,i,j) = tedge               ! Edge distance from wall, delta
            redge         = fout(kd1,i,j)                       ! Edge density
            v_total(:)    = fout(ku1:ku2,i,j)                   ! Edge u, v, [w]
            vedge         = sqrt (dot_product (v_total, v_total))
            v_tangent(:)  = v_total(:) - dot_product (v_total, unit_normal) *  &
                            unit_normal(:)
            vtangent_edge = sqrt (dot_product (v_tangent, v_tangent))

!           Set up the ordinates to integrate for displacement thickness.
!           Safeguarding a stag. point requires forming (rho*v) / (rhoe*ve),
!           not (rho/rhoe) * (v/ve) as done originally:

            rho_vtan_edge = redge * vtangent_edge + v_epsilon
            vtangent_edge =         vtangent_edge + v_epsilon

            numk = kedge + 1

            do k = 1, numk
               v_total(:)   = f(iu1:iu2,i,j,k)
               v_tangent(:) = v_total(:) - dot_product (v_total, unit_normal)* &
                              unit_normal(:)
               vtangent     = sqrt (dot_product (v_tangent, v_tangent))
               rratio(k)    = (f(id1,i,j,k)*vtangent + v_epsilon) / & ! (rho*v)/
                              rho_vtan_edge                      ! (edge rho*v)
               vratio(k)    = (vtangent + v_epsilon) / vtangent_edge
               fline(k)     = one - rratio(k)      ! 1 - (rho*v) / (edge rho*v)
            end do

!           Option to replace delstar with velocity thickness:

            if (velocity_thickness) then
               do k = 1, numk
                  fline(k) = one - vratio(k)
               end do
            end if

            call lcsquad (numk, tline, fline, zero, tedge, method, delstar)

            fout(kds,i,j) = delstar
            if (delstar == zero) then  ! 2D stag. pt. only?
                istag = i              ! We will extrapolate instead
                jstag = j
            end if

            if (save_profile) then
               if (velocity_thickness) then
                  write (luncrt, '(a, f10.6, a)') &
                     'Profile Integrated for Velocity Thickness = ', delstar,  &
                     '  Meters'
                  write (luncrt, '(a)') subtitle, &
                     '1 - (vtan / edge vtan)', 'Wall distance, meters', &
                     ' $options grid=1, symbol=3, $end'
               else
                  write (luncrt, '(a, f10.6, a)') &
                     'Profile Integrated for Delta-* = ', delstar, '  Meters'
                  write (luncrt, '(a)') subtitle, &
                     '1 - (rho vtan) / (edge rho vtan)', 'Wall distance, m',   &
                     ' $options grid=1, symbol=3, $end'
               end if
               write (luncrt, '(2es19.11)') (fline(k), tline(k), k = 1, numk)
               write (luncrt, '(/, a, /)') ' $options line=''dash'', $end'
               write (luncrt, '(2es19.11)') zero, tedge, one, tedge
               write (luncrt, '(/, a, /)') 'END FRAME'
            end if

!           Likewise for momentum thickness:

            do k = 1, numk
               fline(k) = rratio(k) * (one - vratio(k))
            end do

            call lcsquad (numk, tline, fline, zero, tedge, method, theta)

            fout(kth,i,j) = theta

            if (save_profile) then
               write (luncrt, '(a, f10.6, a)') &
                  'Profile Integrated for Theta = ', theta, '  Meters'
               write (luncrt, '(a)') subtitle, &
                  '(rho vtan / edge rho vtan) * (1 - (vtan / edge vtan))',     &
                  'Wall distance, meters', ' $options grid=1, symbol=3, $end'
               write (luncrt, '(2es19.11)') (fline(k), tline(k), k = 1, numk)
               write (luncrt, '(/, a, /)') ' $options line=''dash'', $end'
               write (luncrt, '(2es19.11)') zero, tedge, one, tedge
               write (luncrt, '(/, a, /)') 'END FRAME'
            end if

!           Re-ue (unit Reynolds number based on edge conditions):

            fout(kue,i,j) = (redge * vedge) / fout(kv1,i,j)     ! Edge viscosity

!           Film coefficient:
!           Bad profiles have been observed with Hwall = Hedge = Hinfinity

!!!         ypeval = max (fout(kh1,i,j) - fout(jh1,i,j), 1.e-6) ! Reuse ypeval
            ypeval = max (Hinf_use      - fout(jh1,i,j), 1.e-6)

            fout(kch,i,j) = fout(jq1,i,j) / ypeval
                            ! Wall Qdot   / (Hinf  - Hwall) now, not
                            ! Wall Qdot   / (Hedge - Hwall)

!           Interpolate kappa at the edge?

            if (nt > 1) then

               do k = k1, k1 + 3
                  yline(k) = f(ik1,i,j,k)
               end do

               call lcsfit (4, tline(k1), yline(k1), new, method, 1, tedge,    &
                            yeval, ypeval)
               fout(kk1,i,j) = yeval

            end if

!           Extras interpolated at the edge?

            if (nx > 0) then

               do n = kx1, kx2

                  in = ix1 + n - kx1

                  do k = k1, k1 + 3
                     yline(k) = f(in,i,j,k)
                  end do

                  call lcsfit (4, tline(k1), yline(k1), new, method, 1, tedge, &
                               yeval, ypeval)
                  fout(n,i,j) = yeval

               end do

            end if

!           Re-kk (Reynolds # based on roughness height k and conditions at k):
!           -------------------------------------------------------------------

            if (nr > 0) then

!              Re-kk requires interpolation of density & velocity at height k.
!              Interpolate viscosity as well.

               if (.not. Re_kk) height = max (theta, v_epsilon)  ! < 0. observed

               do l = 1, 2 + nd ! Quantities to interpolate

                  if (l == 1) in = id1   ! density
                  if (l == 2) in = iv1   ! viscosity
                  if (l == 3) in = iu1   ! u
                  if (l == 4) in = iu1+1 ! v
                  if (l == 5) in = iu2   ! w

                  yline(1:nk_height_k) = f(in,i,j,1:nk_height_k)

                  call lcsfit (nk_height_k, tline, yline, new, method, 1,      &
                               height, height_k_values(l), ypeval)
               end do

               v_total(:) = height_k_values(3:2+nd)
               V = sqrt (dot_product (v_total, v_total))

               fout(lhk,i,j) = height             !  k or theta
               fout(ldk,i,j) = height_k_values(1) !  density
               fout(lvk,i,j) = V                  ! |velocity|
               fout(lmk,i,j) = height_k_values(2) !  viscosity

!              Re-kk, Re-theta, or Re-theta / Medge:

               wall_mu = f(iv1,i,j,1)

               if (Re_kk) then
                  fout(lkk,i,j) = (height_k_values(1) * V * height) / wall_mu
               else
                  fout(lkk,i,j) = fout(kue,i,j) * height  ! Re-u * theta
                  if (Re_theta_Medge) &
                     fout(lkk,i,j) = fout(lkk,i,j) / max (fout(km1,i,j), &
                                                          v_epsilon)
               end if

               if (save_profile) then

                  do k = 1, nk_height_k   ! Just the inner radial line grid pts.
                     v_total(:) = f(iu1:iu2,i,j,k)
                     V          = sqrt (dot_product (v_total, v_total))
                     yline(k)   = f(id1,i,j,k) * V * tline(k) / wall_mu
                  end do

                  write (luncrt, '(/, (a))') 'Re-kk Profile', subtitle, &
                                             'k, meters', 'Re-kk'
                  write (luncrt, '(2es19.11)') &
                     (tline(k), yline(k), k = 1, nk_height_k)

               end if

            end if

         end do ! Next i for this block

      end do ! Next j for this block

      deallocate (xline, yline, zline, tline, fline, Hratio, rratio, vratio)

      if (nbad > 0) write (luncrt, '(/, a, i4, a, i7)') &
         ' *** Number of troublesome profiles for block', ib, ':', nbad

!     The 2-D case lends itself to use of running length along the surface.
!     We only return cumulative arc lengths for the current block.  Maybe
!     Tecplot can help concatenate these for multiblock cases.

      if (nd == 2) then  ! Avoid CHORDS2D in favor of CHORDS3D used elsewhere

         allocate (xline(ni), yline(ni), zline(ni), tline(ni))

         if (x(ni,1,1) < x(1,1,1)) then
            call rverse (ni, x, xline)
            call rverse (ni, y, yline)
            call rverse (ni, z, yline)
         else
            xline(:) = x(:,1,1)
            yline(:) = y(:,1,1)
            zline(:) = z(:,1,1)
         end if

         call chords3d (ni, xline, yline, zline, false, total, tline)

         zout(:,1) = tline(:)

         deallocate (xline, yline, zline, tline)

      end if

      end subroutine process_block

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine stag_pt_delstar_theta ()  ! Replace zero delstar and theta
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Local constants:

      character (1), parameter :: method = 'B'  ! Plain Hermite cubic

!     Local variables:

      integer :: i, i1, i2, il, inc, ndat
      real    :: fextrap, total, unused
      real, dimension (5) :: fline, tline, xline, yline, zline

!     Execution:

      if (istag == 1 .or. istag == ni) then  ! Treat as 2D/axisymmetric

!        Extrapolate from the nearest 4 grid points in the i direction:

         if (istag == 1) then
            i1  = 1
            i2  = min (5, ni)
            inc = 1
         else   ! istag = ni
            i1  = ni
            i2  = max (ni - 4, 1)
            inc = -1
         end if
         ndat = abs (i2 - i1)

         il = 1
         do i = i1, i2, inc
            xline(il) = xyzq_in(ib)%x(i,jstag,1)
            yline(il) = xyzq_in(ib)%y(i,jstag,1)
            zline(il) = xyzq_in(ib)%z(i,jstag,1)
            fline(il) = xyzq_out(ib)%q(kds,i,jstag,1)  ! Delta*
            il = il + 1
         end do

         call chords3d (ndat, xline, yline, zline, false, total, tline)
         call lcsfit (ndat, tline(2), fline(2), true, method, 1, tline(1), &
                      fextrap, unused)
         xyzq_out(ib)%q(kds,istag,jstag,1) = fextrap
         write (luncrt, '(a, es16.8, a, i4, a, 2i4)') &
            'Stag. point delta* by extrapolation:', fextrap, &
            ' for block', ib, ' at i, j =', istag, jstag

         il = 1
         do i = i1, i2, inc
            fline(il) = xyzq_out(ib)%q(kth,i,jstag,1)  ! Theta
            il = il + 1
         end do

         call lcsfit (ndat, tline(2), fline(2), true, method, 1, tline(1), &
                      fextrap, unused)
         xyzq_out(ib)%q(kth,istag,jstag,1) = fextrap
         write (luncrt, '(a, es16.8, a, i4, a, 2i4)') &
            'Stag. point theta  by extrapolation:', fextrap, &
            ' for block', ib, ' at i, j =', istag, jstag

      else  ! Stag. point at an interior grid point?  Very unlikely.
      end if

      end subroutine stag_pt_delstar_theta

   end program blayer

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

!  The j direction (unless 2D case):

   if (nj > 1) then
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
   else
      average_ds(3:4) = 99999.
   end if

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
   subroutine blayer_edge (hybrid, shift_ratio, n, Hratio, Wdistance, Hignore, &
                           save_profile, subtitle, kedge, Hedge, tedge, ier)
!
!  Estimate the boundary layer edge along a grid line assumed to be normal to
!  the wall in the region of interest.  A curvature-based method is used for the
!  given curve in (Hratio, distance) space.  The hope is to improve on the
!  standard approach of defining the edge as where the total enthalpy first
!  reaches a value such as 99.5% of the free stream value.
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
!     >  Normalize enthalpy values H by the value Hinf, normally at k = kmax
!     >  Actually, work with the total enthalphy profile H shifted by Hwall
!        to deal with cases like the Mars atmosphere where Hinf is negative
!     >  Consider Hratio ("y") vs. normalized wall distance t ("x")
!        >  Locate the first index i1 for which Hratio > 0.95 (say)
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
!     06/19/04 D.A. Saunders Initial implementation in modular form.
!                            Who proposed the peak curvature idea?
!                            Probably Dinesh Prabhu and/or James Reuther.
!     06/24/04   "     "     The cut-off value is problematic.  Lowering it to
!                            give cleaner lee-side results can hurt wind-side
!                            cases which are the only cases where the idea is
!                            applicable.  Also, starting the spline at the
!                            cut-off data point means the calculated curvatures
!                            are not independent of the cut-off value.
!                            Therefore: spline the whole (normalized) profile,
!                            or at at least the bulk of it, but use just the
!                            data defined by the cut-off value.
!     08/08/04   "     "     The plain curvature method (and parametric splines)
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
!     08/11/04   "     "     Hignore = 0.975 allowed spurious minor peaks in
!                            curvature to be selected.  Raised it to 0.98,
!                            but will this now miss some legitimate peaks?
!                            The real problem is a flattening in the GASP
!                            profiles centered around 99%, and it's not clear
!                            which of two curvature peaks to pick.  One seems
!                            too near the surface; the other too far away.
!                            DPLR profiles do not behave this way.
!     08/12/04   "     "     Seek another peak if refined Hpeak < Hignore.
!                            Introduced tlimit = 0.5 to limit some extremes
!                            likely with lee-side data (ignored anyway).
!     12/09/04   "     "     Added y' and y" to the curvature plot for the
!                            indicated profile, hoping to understand some noisy
!                            edge thickness calculations for a 2-D wedge.
!     12/10/04   "     "     These derivatives seem to have no useful features
!                            in the region of interest.  Therefore, stay with
!                            curvature, but introduce some explicit smoothing
!                            to help with cases observed where the curvature is
!                            essentially constant for several grid points.
!                            The earlier moderation of (now smoothed) curvature
!                            is retained because without it, the peak is too
!                            inclined to be very close to a grid point.
!     12/13/04   "     "     Shuttle results suffer if every profile has its
!                            curvature smoothed, even with one iteration.
!                            Therefore, confine the smoothing to profiles that
!                            have two or more curvatures within a few percent
!                            of each other at the relevant peak.
!     02/11/05   "     "     Odd results on the side of a 3-D blunt wedge near
!                            the nose were traced to redefining kleft as kpeak
!                            before a retry due to Hedge < Hignore.  Thanks
!                            again go to Ryan McDaniel's wedge case.  The choice
!                            of Hignore is still more arbitrary than we would
!                            like, as some well-behaved profiles have their
!                            peak curvature slightly below 0.980, say.  Set it
!                            to 0.975 now instead of 0.980.
!     01/12/06   "     "     Further study of Shuttle solutions suggest lowering
!                            this limit to 0.95.  Curve fit experiments suggest
!                            normalizing so that delta is around 0.1 is a
!                            reasonable choice if we iterate to account for the
!                            dependence of curvature on scaling/normalization.
!                            However, normalizing by something like twice the
!                            initial edge height estimate is clearly not right.
!                            It changes the result drastically because of how
!                            the scaling affects curvature.
!     01/14/06   "     "     Iterating is dubious: any rescaling of arc lengths
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
!     04/22/07   "     "     High-density cases (256 points off the wall) showed
!                            the k window around the curvature-based peak needs
!                            to be larger for larger nk.  The root of this
!                            problem is spurious local curvature peaks, which
!                            are more likely with denser spacing.  A bigger
!                            window for the traditional method should find the
!                            proper peak enthalpy ratio for denser meshes.
!     05/17/07   "     "     Switched to finite difference derivatives for the
!                            initial curvature calculations at the data points,
!                            because they are not affected by far away data.
!     06/26/07   "     "     Non-uniform in-flow (as from an arc-jet nozzle)
!                            means cutting off data below H/Hinf = 0.95 is not
!                            viable.  Therefore, if the outer part of the
!                            profile appears non-uniform, relax this cut-off
!                            to (say) 0.40 (now a user input).
!     06/19/08- "       "    Mike Olsen suggested using 95% for stage 2 instead
!     06/24/08               of 99.5%.  This would be less sensitive to the
!                            heuristic definition of "neighborhood."  However,
!                            in the absence of a stage 3 that recovers the 99.5%
!                            location for well-behaved profiles, we leave well
!                            alone.  Only minor refinements have been made.
!     02/12/13  "       "    Belated realization that we should work with the
!                            ratio (H - Hwall) / (Hinf - Hwall), not H / Hinf.
!                            This affects only the documentation here, since
!                            Hratio is evaluated at the higher level.  Results
!                            will differ, however, because the shape of the
!                            ratio profile is different.
!     06/11/14  "       "    Raised klast from 3*nk/4 to 4*nk/5.
!     10/11/17  "       "    Added hybrid and shift_ratio arguments to control
!                            saved profile details.
!     11/17/17  "       "    Raised tlimit from 0.5 to 0.7 on account of cases
!                            with nk = 101 (low) and thick boundary layers at
!                            high altitude.
!
!  Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA
!           Now with ERC, Inc. at NASA ARC (August 2010 through June 2015).
!           Now with AMA, Inc. at NASA ARC.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   logical,   intent (in)  :: hybrid       ! T: 2-stage curvature-based method
   logical,   intent (in)  :: shift_ratio  ! T: Hratio=(H+Hshift)/(Hinf+Hshift);
                                           ! F: Hratio=(H-Hwall) /(Hinf-Hwall)
   integer,   intent (in)  :: n            ! # pts. in the radial line of data
   real,      intent (in)  :: Hratio(n)    ! Total enthalpy ratio values
   real,      intent (in)  :: Wdistance(n) ! Wall distances
   real,      intent (in)  :: Hignore      ! Ignore curvature data below this
                                           ! value of Hratio
   logical,   intent (in)  :: save_profile ! T means write plottable data to
                                           ! standard output
   character (32), intent (in) :: subtitle ! Annotates the particular profile
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
   real,    parameter :: power   = -0.1, & ! Moderates curvature; reciprocal
                         rpower  = one/power
   real,    parameter :: small_fraction  & ! Criterion for smoothing curvature
                                 = 0.10
   real,    parameter :: tlimit  = 0.7     ! Upper limit on normalized wall
                                           ! distance allowed for edge value
   real,    parameter :: tol     = 1.e-8   ! Acceptable uncertainty in the
                                           ! BL edge arc length (s, not t)
   logical, parameter :: false   = .false.

   character (7), parameter :: subname = 'bl_edge'

!  Local variables:

   integer :: i, istat, k, k1, k2, kedge_interim, klast, kleft, kpeak,         &
              kpeak_previous, kwindow, l1, l2, lunout, nfit, nfit2, niter, nt, &
              numfun

   real    :: t1, t2, tmin, flattened_min, fmin, Hmax, redge, stotal, term, &
              ttotal, Hedge_interim, tedge_interim

   logical :: write_profile

   real, dimension (5) :: &
      b3, c3, d3, df, dff, flattened_kappa

   real, dimension (:), allocatable :: &
      b1, c1, d1, dy, dyy, f, kappa, t

!  Execution:
!  !!!!!!!!!!

   write_profile = save_profile .and. .not.hybrid  ! Avoid duplicate output

!  Ignore data near the wall with enthalpy ratio < Hignore?  NO:  The curvature
!  values used are then not independent of Hignore.  Therefore, err on the side
!  of more than enough data very near the wall.

!!!   kleft = n / 2
!!!   call interval (n, Hratio, Hignore, one, kleft)

   kleft = 8

!  Ignore data beyond some fraction of the maximum wall distance, assuming the
!  calling program hasn't already suppressed the outer portion of the data:

   klast = (4 * n) / 5
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

!  Calculate conventional spline coefficients for enthalpy ratio vs. t:

   call csfit (nfit, t, Hratio(kleft), iendl, zero, iendr, zero, b1, c1, d1,ier)

   if (ier /= 0) then
      write (lunerr, '(/, 2a, i4)') subname, &
         ':  Trouble splining enthalpy ratio vs. t; ier:', ier
      go to 99
   end if

!  We have switched to finite differece differences for calculating the
!  profile curvature distribution, but we still need the above spline
!  coefficients to interpolate some items below at the peak curvature location.

!! Compute 1st and 2nd derivatives of enthalpy ratio at the data points:

!! call csdval (nfit, t, Hratio(kleft), nfit, t, b1, c1, d1, f, dy, dyy)

!! Derive the curvature at the data points:

!! do k = 1, nfit
!!    term = one + (min (abs (dy(k)), 1.e+10)) ** 2
!!    kappa(k) = -dyy(k) / (term * sqrt (term)) ! Sign is consistent
!! end do                                       ! with parametric version

!  Finite difference derivatives should behave a bit better than spline
!  derivatives, and are not affected by data far from the edge region.

   call fd12k (nfit, t, Hratio(kleft), dy, dyy, kappa)

   kappa(:) = -kappa(:)  ! For consistency with earlier development work

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
      if (shift_ratio) then
         write (lunout, '(/, a)') &
            '# For specified ((H+Hshift)/(Hinf+Hshift), distance) profile:'
      else
         write (lunout, '(/, a)') &
            '# For specified ((H-Hwall)/(Hinf-Hwall), distance) profile:'
      end if
      write (lunout, '(a, i4, 2es19.11)') &
         '# kpeak, t1, t2:', kpeak+kleft-1, t1, t2
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
         redge = Hedge                           ! Interim ratio at peak kappa;
                                                 ! trad. method may change it
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
      tmin = tedge / ttotal  ! In case of fminrc failure?
      if (niter > 1) then
         write (lunout, '(a, i2, a, f10.6)') &
            '# Normalized peak-curvature delta at iteration', niter, ':', tmin
      end if
   end if

   if (niter < maxiter) then ! Iterate to achieve normalized edge height ~ .10
       niter = niter + 1     ! This was a bad idea; maxiter = 1 now
      ttotal = tedge / 0.10  ! ??
      go to 10
   end if

!  Now that we have the most likely neighborhood of the edge height, adjust
!  for possible anomalies and for the arc-length dependency by seeking 99.5%
!  of the PEAK enthalpy ratio in this neighborhood.

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

   if (Hignore >= 0.95) then  ! Do the second stage

      Hedge = Hmax * 0.995

!     Make sure that k1 isn't too high to include the Hedge:

      do k = kedge - 1, 10, -1
         if (Hratio(k) <= Hedge) then
            k1 = k
            exit
         end if
      end do

      call bl_edge_traditional (shift_ratio, Hedge, n, k1, k2, Hratio, &
                                Wdistance, write_profile, kedge, tedge, ier)
      if (ier /= 0) then
         Hedge = Hedge_interim
         tedge = tedge_interim
         kedge = kedge_interim
      end if

   else
      ! Live with the slight scaling dependency for non-uniform in-flow cases
      ! or for lee-side flows such as on the front of the Shuttle OMS pod,
      ! where Hignore = 0.5 or lower is necessary to find what is the more
      ! likely boundary layer thickness (much smaller than that found with
      ! the default of 0.95).
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
      if (shift_ratio) then
         write (lunout, '(/, a)') '# (H+Hshift)/(Hinf+Hshift)  Wall dist.   k'
      else
         write (lunout, '(/, a)') '#  (H-Hw)/(Hinf-Hw)      Wall distance   k'
      end if
      write (lunout, '(2es19.11, i4)') (Hratio(k), Wdistance(k), k, k = 1, n)

      write (lunout, '(a)') &
         ' $options line=''symbol'', symbol=15,'
      write (lunout, '(a, f8.5, a, f8.5, a)')       &
         ' legend = ''Boundary layer edge is at (', Hedge, ',', tedge, ' )'',',&
         ' $end'
      write (lunout, '(/, a, 2i4, /, 2es18.11)')    &
         '# kedge, kleft, Hedge, tedge:', kedge, kleft, Hedge, tedge
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
      if (shift_ratio) then
         write (lunout, '(/, a)') '# (H+Hshift)/(Hinf+Hshift) Norm. dist.   k'
      else
         write (lunout, '(/, a)') '#  (H-Hw)/(Hinf-Hw)   Norm. wall dist.   k'
      end if
      write (lunout, '(2es19.11,i4)') (Hratio(k), t(k-kleft+1), k, k = kleft, n)

      write (lunout, '(a)') &
         ' $options line=''symbol'', symbol=15,'
      write (lunout, '(a, f8.5, a, f8.5, a)')       &
         ' legend = ''Peak curvature edge is at (', redge, ',', tmin, ' )'',', &
         ' $end'
      write (lunout, '(/, a, 2i4, /, 2es18.11)') &
         '# kedge, kleft, Hedge, tmin:', kedge, kleft, Hedge, tmin
      write (lunout, '(/, a, /)') 'END FRAME'

      write (lunout, '(a)') &
         'Boundary Layer Profile',                  &
         subtitle,                                  &
         'Normalized wall distance',                &
         'Derivatives/curvature of Hratio vs. normalized wall distance',       &
         ' $options height=7, ymin=0, ymax=50,',    &
         ' line=''solid'', symbol=16, ident=T,',    &
         ' legend = ''Smoothed curvature'', $end',  &
         '#         t                  Curvature   k'
      write (lunout, '(2es19.11, i4)') &
         (t(k), kappa(k), k+kleft-1, k = 1, nfit)

      write (lunout, '(a)') &
         ' $options line=''dash'',',                &
         ' legend = ''First derivative'', $end',    &
         '#         t                      dy/dt   k'
      write (lunout, '(2es19.11, i4)') &
         (t(k), dy(k), k+kleft-1, k = 1, nfit)

      write (lunout, '(a)') &
         ' $options line=''longdash'',',            &
         ' legend = ''Second derivative'', $end',   &
         '#         t                    d2y/dt2   k'
      write (lunout, '(2es19.11, i4)') &
         (t(k), -dyy(k), k+kleft-1, k = 1, nfit)

      write (lunout, '(a)') &
         ' $options line=''symbol'', symbol=15,',   &
         ' legend = ''BL edge curvature, smoothed'', $end'
      write (lunout, '(2es19.11)') tmin, fmin
      write (lunout, '(/, a, /)') 'END FRAME'

      write (lunout, '(a)') &
         'Boundary Layer Profile',                  &
         subtitle,                                  &
         'Normalized wall distance',                &
         'Moderated curvature in region of peak',   &
         ' $options height=7, fit=''loose'',',      &
         '  spline=''standard'', ident=T, $end',    &
         '#         t        Moderated Curvature   k'
      write (lunout, '(2es19.11, i4)') &
         (t(l1+k-1), flattened_kappa(k), k+l1+kleft-2, k = 1, nfit2)

      write (lunout, '(a)') &
         ' $options line=''symbol'', symbol=15,',   &
         ' legend = ''Boundary layer edge'', $end'
      write (lunout, '(2es19.11)') tmin, flattened_min
      write (lunout, '(/, a, /)') 'END FRAME'

      end subroutine save_plottable_details

   end subroutine blayer_edge

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine bl_edge_traditional (shift_ratio, Hedge, nk, k1, k2, Hratio, &
                                   Wdistance, save_profile, kedge, delta, ier)
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
!  04/02/13    "     "     The specified enthalpy ratio profile is written to
!                          STDOUT if new argument save_profile = T.
!  10/11/17    "     "     Added shift_ratio argument to control profile output.
!
!  Author:  David Saunders, ELORET Corporation/NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   logical,   intent (in)  :: shift_ratio  ! T: Hratio=(H+Hshift)/(Hinf+Hshift);
                                           ! F: Hratio=(H-Hwall) /(Hinf-Hwall)
   real,      intent (in)  :: Hedge        ! Usually 0.995 (target Hratio)
   integer,   intent (in)  :: nk           ! # pts. in the radial line of data
   integer,   intent (in)  :: k1, k2       ! Data range of interest
   real,      intent (in)  :: Hratio(nk)   ! Enthalpy ratios, (H-Hw)/(Hinf-Hw)
   real,      intent (in)  :: Wdistance(nk)! Wall distances
   logical,   intent (in)  :: save_profile ! T => write profile data to STDOUT
   integer,   intent (out) :: kedge        ! [kedge, kedge+1) contains delta
   real,      intent (out) :: delta        ! Boundary layer thickness estimate
   integer,   intent (out) :: ier          ! 0 means no error

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

   if (save_profile) then
      write (*, '(/, a, f10.6)') '# Traditional method seeking Hedge =', Hedge
      if (shift_ratio) then
         write (*, '(/, a)') '# (H+Hshift)/(Hinf+Hshift)  Wall dist.   k'
      else
         write (*, '(/, a)') '#  (H-Hw)/(Hinf-Hw)      Wall distance   k'
      end if
      write (*, '(f19.15, es19.11, i4)') (Hratio(k), Wdistance(k), k, k = 1, nk)
      write (*, '(/, a, i4, es15.7)') '#  kedge, delta:', kedge, delta
   end if

99 return

   end subroutine bl_edge_traditional

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine smooth1d (i1, i2, s, v)
!
!     SMOOTH1D smooths a vector in-place by an explicit method with hard-coded
!     parameters.  End elements are unchanged.
!
!     06/23/97  D.Saunders  Gridgen-type explicit smoothing to cope with spikes
!                           in edge Phi, Psi.
!     12/10/04      "       5 iters. seem too many for boundary layer edge work.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

!     Arguments:

      integer, intent (in)    :: i1, i2   ! Index range to be smoothed
      real,    intent (in)    :: s(i1:i2) ! Abscissas associated with v(i1:i2)
      real,    intent (inout) :: v(i1:i2) ! Vector being smoothed; elements
                                          ! i1 & i2 are left unchanged
!     Local constants:

      integer, parameter :: niters = 1    ! Originally 5
      real,    parameter :: alpha = 0.5, one = 1.

!     Local variables:

      integer :: i, iter
      real    :: ai, bi, ci, hl, hr, term
      real    :: w(i1:i2) ! Work-space for previous iterate (to vectorize)

!     Execution:

      do iter = 1, niters

         w(i1:i2) = v(i1:i2)

         do i = i1 + 1, i2 - 1
            hl   = s(i) - s(i-1)
            hr   = s(i+1) - s(i)
            term = alpha * ((min (hl, hr) ** 2) / (hl + hr))
            ai   = term / hl
            ci   = term / hr
            bi   = one - (ai + ci)
            v(i) = ai * w(i-1) + bi * w(i) + ci * w(i+1)
         end do

      end do

      end subroutine smooth1d
