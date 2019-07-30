**Outline and Formalities**
---------------------------

On January 24, 2014, this project received final approval from NASA to
upload actual software.  This process was completed on February 5, 2014,
although a few items in the one-liners have been held back for now.
Provided are 30+ subroutine libraries (mostly lowercase library names, with
a total of about 1200 *.f and *.f90 subroutine files) and more than 100
application programs (uppercase names), almost all in Fortran [90]. Most
have been developed at NASA Ames Research Center over several decades of
work in the Aerodynamics Division and Space Technology Division.  Most are
computational fluid dynamics-related (CFD), and many apply to multiblock
structured data in PLOT2D/3D form (grids and function files).  A number of
the applications work with datasets in Tecplot format.  (Tecplot is a de
facto standard at NASA.)  The mixture of generalized and specialized
functionality offered is apparent from the one-liners below.  Subroutine
one-liner READMEs also appear with each library, while each application
folder contains a README adapted from the main program header.

**Use and Change Notice**

Please refer to **CFD_Utilities_License_NOSA.pdf** on this **Files** page
for details of the NASA Open Source Agreement under which the CFD Utilities
are released.

Most of the subroutines and application programs are single-purpose
utilities that are unlikely to warrant modification.  If any would-be user
chooses to make changes, the present authors urge that the History section
found in each file be updated appropriately, along with other in-line
documentation such as argument description or purpose outline that is
affected by such changes.  We do not require notification of any such
changes, but would welcome any associated feedback.

**Contacts**

David.A.Saunders@NASA.gov, Todd.R.White@NASA.gov, or Jeffrey.P.Hill@NASA.gov

**Downloading**

Guidelines appear below and in the file **Downloading.txt** on this page.

**Libraries/Modules**
---------------------

**aa_nurbs                ** | Supplements the DT_NURBS library, for B-spline airfoils
**adt                     ** | Rapid search/interpolation within computational meshes
**blaslib                 ** | Basic linear algebra subprogram utilities (vectors)
**blas2lib                ** | NAG versions of the Level 2 BLAS Matrix-vector routines
**C_utilities             ** | Handful of C translations from interplib, linsys, etc.
**dt_nurbs                ** | (Dated) David Taylor Naval Research Center NURBS library
**eigenlib                ** | Just the EISPACK routines used by CAVITY_MAP rotations
**forsythe                ** | Original Forsythe, Malcolm & Moler library (historical)
**geomlib                 ** | Areas, rotations, nearest points, arc lengths, etc.
**grid_block_utilities    ** | Volume grid analogues of the surface_patch_utilities
**gridlib                 ** | Mostly for structured 2/3-space grids; many 1-D distributions
**integrate               ** | Numerical quadrature, mostly 1-D, with a driving program
**interplib               ** | 1-D interpolation & smoothing utilities; see SMOOTH
**interp2d                ** | 2-D interpolation, smoothing, & intersection utilities; see SMOOTH2D
**interp3d                ** | 3-space interpolations, intersections, etc.
**intrinsics              ** | Contains only trig_functions_in_degrees: g95 cosd, etc.
**kdtree                  ** | Matthew Kennel's open source package, now an option in FLOW_INTERP
**linsys                  ** | Linear system solvers, mostly dense factorization methods
**numodules               ** | Miscellaneous utilities: finite differences, standard deviations, etc.
**obsolete                ** | Still linked to by INTEGRATE, MAXMIN, PROFILE, and SMOOTH
**Optimal_Interpolation   ** | Alexander Barth's Kriging-like package, with extensions
**optlib                  ** | 1-D zero/minimum finders; n-D unconstrained methods, etc.
**progtools               ** | String-manipulation & prompting utilities, incl. string_justify module
**searchlib               ** | 1-D search/sort routines; some for 2|3-D now deprecated
**special_functions       ** | A few utilities related to erf(x), normal distributions, etc.
**surface_patch_utilities ** | For structured surface grids: scale/shift/transpose, etc.
**table_io                ** | I/O utilities for plain text tables (headers + numeric rows)
**tecplot_io              ** | Tecplot I/O package for structured multizone data
**triangulation_io        ** | Tecplot I/O package for triangulated surface data
**ugridlib                ** | Triangle and tetrahedron utilities
**xyq_io                  ** | PLOT2D-type multiblock structured grid I/O utilities, in xyzq_io folder
**xyzq_io                 ** | PLOT3D-type multiblock structured grid I/O utilities; see also f_io.f90

**Applications**
----------------
    C = CFD, G = General Numerics, P = Programming Tool

**A2B                ** | G   Converts a (large) ASCII file to a binary stream file.
**ADJUST_FLOW        ** | C   Analogue of ADJUST_GRID for PLOT3D-format flow solutions
**ADJUST_GRID        ** | C   20-odd options: shift/scale/rotate/swap coordinates, etc.
**Aero_Coefs         ** | C   Aerodynamic force/moment coefs. for structured surface data containing p or Cp
**AEROSURF           ** | C   Wing/body/nacelle surface paneling; underlies shape optimization
**ANCHOR             ** | G   Applies Optimal_Interpolation package to aerodynamic or other unstructured data
**BLAYER             ** | C   Boundary layer analysis of a structured flow solution
**BSPROFILE          ** | C   B-spline curve analogue of PROFILE airfoil utility
**BUMP_FACTORS       ** | C   Compare related structured surface solutions as ratios
**CAPSULE_GRID       ** | C   Automated surface grid for sphere/cone; many aft body options
**CAPSULE_SPOKES     ** | C   Supplements POLAR_INTERP by setting up spoked body points for a capsule
**CAVITY_MAP         ** | C   Map local cavity flow solution to unit cube; slice it, etc.
**CH_GRID            ** | C   Stand-alone form of the automated C-H volume gridding in wing/body design code SYN87
**COLUMNEDIT         ** | P   Insert/delete/replace column(s) in a dataset
**COMBINE_BLOCKS     ** | C   Combine multiblock grids/solutions; cavity block options
**COMBINE_BLOCKS_TURB** | C   COMBINE_BLOCKS variant required for turbulent flows
**COMPARE_BLOCKS     ** | C   Compare related structured grids: identical? maximum differences?
**COMPARE_FLOWS      ** | C   Flow field analogue of COMPARE_BLOCKS (volume | surface)
**COMPARE_PATCHES    ** | C   Surface grid analogue of COMPARE_BLOCKS
**COMPRESS2D         ** | C   CAPSULE_GRID tool; adjust upstream boundary, wall spacing, # radial grid points, etc.
**CONES_OF_SIGHT     ** | C   Body-normal sets of 9 discretized lines forming cones for radiation calculations
**CROSS_SECTIONS     ** | C   Slice a surface dataset (structured | unstructured) at many stations
**CURVATURE          ** | C   Evaluate curvature along X-Y curve; redistribute points according to curvature
**CURVATURE3D        ** | C   X-Y-Z curve analogue of CURVATURE
**DECONSTRUCT        ** | C   Convert a multiblock structured surface or volume grid to unstructured form
**DISTRIBUTE         ** | C   Drives available 1-D grid-point distribution utilities
**EXTRACT_BLOCKS     ** | C   Extract indicated block(s) from a multiblock grid/flow
**EXTRACT_BLOCKS_2D  ** | C   Extract indicated block(s) from a multiblock 2D grid/flow
**EXTRACT_COLUMNS    ** | GP  Extract column(s); options to scale/shift/reformat
**EXTRACT_FUNCTIONS  ** ! C   Extract selected functions efficiently from a PLOT3D function file.
**EXTRACT_LINES      ** | GP  Extract every nth line starting from some line
**EXTRACT_PEAKS      ** | CGP For a list of files of column data, tabulate peak values
**FILTER_LINES       ** | GP  Remove non-numeric lines (excess headers?) from a table
**FLOW_INTERP        ** | C   Volume grid search & flow interpolation (ADT method and KDTREE method)
**FLOW_INTERP_2D     ** | C   2-space analogue of FLOW_INTERP
**GEN1D              ** | G   Generate (x,y) datasets from a choice of y = f(x) functions
**GRADIENT_BASED     ** | CG  Redistribute points along an (x,y) curve according to |df/ds| where s is arc length, for some function f
**GRID_FACES         ** | C   Extract grid block face(s) from a multiblock grid/flow
**GSMOOTH            ** | C   Elliptic volume grid smoothing; option to detect bad cells
**GU                 ** | C   Ryan McDaniel's multi-option tool for PLOT3D grids/flows
**HEAT_SHIELD        ** | C   Forebody shape optimization, axisymmetric option; Newtonian CL/CD
**HEMISPHERES_OF_SIGHT** | C  Sets of discretized lines forming hemispheres for radiation calculations; see also SLOS, USLOS
**HYPER_AERO         ** | C   Driver for available impact methods, including Modified Newtonian
**INSERT_LINES       ** | P   Insert lines from one file into another at some interval
**INTERP_1D          ** | G   Interactive linear/local spline interpolation, 1-D
**INTERP_COLUMNS     ** | GP  Adaptation of REGULARIZE; variable column to treat as "x"
**JOIN_GRID          ** | C   SPLIT_GRID companion for merging blocks
**LINE_GRID          ** | GP  Uniform grid between 2 points in 3-space; extrapolation OK
**LINES_OF_SIGHT     ** | C   Body-normal | shock-normal discretized lines for radiation calculations
**LINES_OF_SIGHT_2D  ** | C   2-space analogue of 3-space LINES_OF_SIGHT
**LOC                ** | P   Count lines of code & comments in a Fortran source file
**MAKELIST           ** | P   Make a list of integers for given i1, i2, increment
**MAXMIN             ** | G   Analogue of SMOOTH for looking at 1st & 2nd derivatives
**MERGE_BLOCKS       ** | C   Variant of JOIN_GRID allowing for optional function file
**MERGE_FILES        ** | C   Merge variables of two Tecplot files (same structured grid)
**MERGE_TABLES       ** | GP  Rearrange columns from any number of tables into a new table
**NBYN               ** | G   Interactive solution of square n x n systems; n <= 6
**NEQAIR_DATA        ** | C   Converts line-of-sight data from PLOT3D form to NEQAIR form
**NEQAIR_Integration ** | C   Integrates radiances from hemisphere lines of sight to estimate radiative heat flux
**NOZZLE_THROAT_CONDITIONS  ** | C   Axisymmetric variant of THROAT_CONDITIONS_3D
**NPOPT_DRIVER       ** | G   Framework for a [sequence of] constrained minimization[s]
**OUTBOUND           ** | C   Off-line grid/shock alignment; other redistribution options
**P3D2TEC            ** | C   Convert between PLOT3D grid/flow format and Tecplot form
**POLAR_INTERP       ** | C   Pad spoked forebody radiation data nonlinearly
**PRECISION          ** | G   Estimate objective function precision (difference table)
**PREPARE_LOCAL_ANALYSIS ** | C  Rapid set-up of DPLR-based damage/repair/feature calculations: cavity | protruding gap-filler
**PREPARE_NEQAIR_DATA** | C   Facilitates scripting of radiative heating calculations by NEQAIR
**PROFILE            ** | C   Airfoil geometry display/manipulation utility
**QNMDRIVER2         ** | G   Framework for one unconstrained optimization (QNMDIF2)
**QNMDRIVER3         ** | G   Framework for a sequence of unconstrained optimizations
**RADIAL_INTERP      ** | C   Rapid volume grid|flow interpolation for new surface grid (single layer of volume grid blocks)
**RADIAL_INTERP_2D   ** | C   2-space RADIAL_INTERP analogue; both one layer of grid blocks
**RECTIFY_GRID       ** | C   Ensure right-handedness for 3-space volume or surface grid
**REFINE_GRID        ** | C   Densify | thin multiblock grid/flow data, any multipliers
**REFLECT_BLOCKS     ** | C   Reflect 3-space grid/flow data; save reflected or both halves
**REFLECT_BLOCKS_2D  ** | C   Reflect 2-space grid/flow data; save reflected or both halves
**REORDER_BLOCKS     ** | C   Reorder some or all grid blocks; optional flow file
**REORDER_ROWS       ** | GP  Reorder tabular data so leading n columns vary differently
**REORDER_SEGMENTS   ** | GC  Reorder multiple 1-D Tecplot zones as 1 continuous zone
**RESHAPE            ** | G   For 2 of 2 or more columns, shift/scale/rotate/... (11 options)
**RESHAPE3D          ** | G   3-space analogue of RESHAPE
**REVOLVE_GRID       ** | C   2-space volume grid --> 3-space; may need RADIAL_INTERP next
**SCAN_GRID          ** | C   Tabulate data range, etc., by grid block; optional flow
**SHADOWGRAPH        ** | C   CFD density data --> shadowgraph/schlieren-like image
**SHIFTSCALE         ** | G   Interactive calculation of linear coefficients for transforming [a, b] -> [p, q]
**SHOCK_STAND_OFF    ** | C   For volume grid & flow field temperature, calculate shock stand-off distances
**SLOS               ** | C   USLOS-type merge of LINES_OF_SIGHT & HEMISPHERES_OF_SIGHT (structured grid)
**SMOOTH             ** | G   Drives available 1-D interpolation & smoothing utilities
**SMOOTH2D           ** | G   Bivariate analogue of SMOOTH
**SORT_ROWS          ** | GP  Sorts the output from grep (say) or a list of file names containing numeric substrings.
**SORT_SURFACE_SLICE ** | C   Sort Tecplot surface slice into 1 or more curves for line plots
**SPLIT_GRID         ** | C   Split multiblock grid/flow and/or permute indices
**SURFACE_CURVATURE  ** | C   Gaussian/mean/principal curvatures for surface/volume grid
**SURFACE_DIFFS      ** | C   Map Tecplot surface grid 1 to grid 2; save [%]differences
**SURFACE_INTERP     ** | C   Interpolate 3-space Tecplot surface data at target data point(s)
**SURFACE_INTERP_2D  ** | C   Interpolate 2-space Tecplot surface data at target data point(s)
**SURFACE_PAD        ** | C   Pad structured surface data nonlinearly: 1-D in i and|or j
**SURFACE_PATCHES    ** | C   ADJUST_GRID variant with more surface_patch_utilities options
**TEMPLATE           ** | C   Calculate block interface data and DPLR-type control file (structured volume grid)
**TET_INTERP         ** | C   Interpolate tetrahedral volume data to a structured multiblock grid
**THIN_FLOW          ** | C   THIN_GRID analogue for multiblock flow data
**THIN_GRID          ** | C   Extract every mth/nth[/kth] point in i/j[/k] directions
**THIN_GRID_2D       ** | C   THIN_GRID analogue for 2-space; grid only
**THROAT_CONDITIONS_3D   ** | C   Boundary conditions at arc-jet nozzle: rectangular|circular|semi-elliptic
**TRAIL              ** | GP  Truncate lines and|or strip trailing blanks
**Traj_Fit           ** | C   Nonlinear least squares fitting of f(rho(t),V(t) = C rho(t)^m V(t)^n via the QNMDRIVER2 framework
**TRAJ_OPT           ** | CG  Trajectory optimization: would need an NPOPT license from Stanford University
**TRI_TO_QUAD        ** | C   Impose a structured surface on a Tecplot triangulation
**TRI_TO_TRI         ** | C   Interpolate a single-zone Tecplot triangulation dataset to another 1+-zone target triangulation
**TRIANGULATION_TOOL ** | C   Drives scale/shift/rotate transformations and area/volume/CM/moments of inertia calculations for a Tecplot unstructured surface or volume dataset (1+ zones)
**UPDATE_GRID        ** | C   Replace one or more grid blocks with same-sized block(s) from other file or files
**UPSEQUENCE         ** | C   Upsequence coarse cavity|plu|gap filler solution; impose fine local boundary flow
**USLOS              ** | C   Merge of LINES_OF_SIGHT & HEMISPHERES_OF_SIGHT for unstructured surfaces
**V2C                ** | C   Convert grid [+ optional flow]: cell vertices to centers
**WINGSECTIONS       ** | C   B-spline sections + chord/thickness data --> B-spline wing
**XDECK              ** | P   Ancient precursor of TRAIL: remove trailing | leading characters
**XLINE              ** | P   Remove lines starting with a target string


**Downloading Guidelines**
--------------------------

Presently, it is assumed that would-be users will not want the entire
heterogeneous collection gathered here, but rather some application(s) and
the libraries required for linking.  Therefore, a *.tar.gz file is provided
for each library with a simple build script.  Each application folder
contains the main program, a sample control file if any, a simple build
script, and any other source files not in one of the libraries.

Use of a directory structure reflecting the folder structure of this web
page is recommended.  Suppose /xxx is the root of such a directory
structure, and program CURVATURE is the application of interest:

(1) From /xxx enter "mkdir CURVATURE" and download the files curvature.f,
build, and any others into /xxx/CURVATURE.  The simple build script shows
that the needed libraries are gridlib, geomlib, interplib, interp3d,
linsys, numodules, searchlib, and prmodules from the progtools folder.

(2) Therefore, directories /xxx/geomlib, /xxx/gridlib, and so on need to
be set up first if they are not already there from downloading some other
application.  E.g., place geomlib.tar.gz into /xxx/geomlib and extract the
geomlib files into that directory using "tar -xvozf geomlib.tar.gz".

(3) Edit the file "build" to suit your compiler.  **64-bit arithmetic is
strongly recommended, via the -r8 switch or equivalent.**  Note that all
the source code uses **real** declarations, but real*8 or double precision
is intended via the appropriate compiler switch.

(4) "./build" compiles all subroutines, places the .o files into an
object library such as f90geomlib.a, then removes the .o files.

(5) In the case of Fortran 90 modules, such as /xxx/adt/adt_utilities.f90,
.mod and .o files are **not** removed, but they can be once successful
compilation has been demonstrated.  **(Such modules are compiled into the
application directory by the build file provided for the application.)**

(6) Once all the needed libraries have been built, adjust the "build" in
/xxx/CURVATURE similarly, and "./build" should link the application.

(7) Optionally, a directory /xxx/bin might be added, placed in the path
of the user, and assigned symbolic links to each application, via commands
such as "ln -s /xxx/CURVATURE/curvature ."

**History**
-----------

V13.08

    Aug 30 2013  Initial preparation.
    Sep 20 2013  Initial summary uploads to SourceForge.
    Sep 30 2013  Still awaiting final NASA approval to upload software.
    Oct 17 2013  CAPSULE_GRID document uploaded following 2.5-week government
                 shutdown hiatus.
    Oct 22 2013  Folders for all libraries are now set up with subroutine one-
                 liner READMEs.
    Oct 30 2013  Folders for all applications are now set up with README files.
    Jan 24 2014  NASA approval for uploading the software proper.
    Jan 29 2014  All the library *.tar.gz files are now in place; applications
                 still to do.
    Feb 05 2014  All straightforward applications have been uploaded with simple
                 build scripts.
    Feb 10 2014  UPDATE_GRID has been added; it and GSMOOTH now use xyzq_io, not
                 the earlier cfd_io_package.
    Jun 30 2014  Someone needed the DECONSTRUCT utility.  A few updates since
                 Feb 10 include solid angle utilities in ugridlib, a handful of
                 functions in special_functions, and a reworked form of the
                 HEMISPHERES_OF_SIGHT application that defines lines of sight
                 at a point on a convex body via a triangulated quadrant of a
                 hemisphere rather than via latitude/longitude discretization.
    Jul 03 2014  TRI_TO_TRI has been added; a triangulation_io glitch has been
                 corrected.
    Aug 07 2014  EXTRACT_BLAYER_DATA has been added to the BLAYER folder.
    Aug 14 2014  Optimal_Interpolation now contains the package for treating
                 unstructured data as modified at NASA Ames from Alexander
                 Barth's original, following NASA approval of its added
                 BSD 2-Clause License.
    Sep 05 2014  FILTER_LINES and MERGE_TABLES have been added.
    Sep 24 2014  Traj_Fit has been added.  It curve-fits heat flux or pressure
                 pulse data using f(rho(t),V(t)) = C rho(t)**m V(t)**n within
                 the QNMDRIVER2 framework, as might be needed for TPS sizing
                 at an entry vehicle body point from selected CFD data points.
    Oct 03 2014  The table_io module used by Traj_Fit has been updated.
    Oct 22 2014  The triangulation_io package has been extended with options
                 for calculation of area, volume, CM, and moments of inertia.
                 See the new TRIANGULATION_TOOL for driving these options
                 along with shift/scale/rotate options.
    Nov 17 2014  Unstructured volume analogues of unstructured surface
                 utilities have been added to triangulation_io and are
                 driven by TRIANGULATION_TOOL.
    Mar 18 2015  An apparently inconsequential glitch in BLAYER has been
                 remedied.  This version can also read [unformatted] volume
                 datasets to speed processing of large-grid cases.
    Mar 13 2015  OUTBOUND now has an option to applied a CONSTANT margin
                 as a way of moving the forebody boundary a little while
                 barely affecting the wake boundary.
    Apr 15 2015  The ADT rapid search package (/adt) now has 2-space multi-
                 block curve analogues of the 3-space structured multiblock
                 surface grid utilities.  LINES_OF_SIGHT_2D uses these to
                 deal with reflected axisymmetric 2-space volume grid.
                 SURFACE_INTERP_2D has been added, using the same 2-space
                 ADT utilities.
    Aug 23 2016  TET_INTERP has been added to interpolate tetrahedral volume
                 data to a structured multiblock grid.
    Mar 17 2017  SORT_ROWS (added Feb 14) has been extended to allow more
                 than one numeric string in a text column, and to allow
                 output in descending order if desired.
                 Other recent refinements: COLUMNEDIT now allows extraction
                 of more than one column at a time; the table_io module has
                 an updated table_io_read_alpha utility prompted by SORT_ROWS
                 that handles the case of all lines/rows treated as text and
                 having the same number of columns (with no header records
                 assumed).  Previously, all rows were interpreted as header
                 records, meaning the table appeared empty.
    Jun 23 2017  NEQAIR_DATA, NEQAIR_Integration, and PREPARE_NEQAIR_DATA
                 have been added to assist radiation calculations with
                 NEQAIR.  NEQAIR itself should be obtained through the
                 Commercial Technology Office at NASA Ames Research Center.
    Jan 04 2018  Added A2B utility.  Since June, minor changes have been
                 made to GU, BLAYER, LINES_OF_SIGHT, SHOCK_STAND_OFF, and
                 OUTBOUND.
    Jun 06 2018  Added Aero_Coefs application (supplement to DPLR's Postflow).
                 Since January, updates have been made to geomlib, gridlib,
                 progtools, searchlib, CAPSULE_GRID, CONVERTQ, GU, HEMISPHERES_-
                 OF_SIGHT, NEQAIR_Integration, RADIAL_INTERP, and RESHAPE3D.
    Sep 19 2018  Added EXTRACT_BLOCKS_2D and EXTRACT_FUNCTIONS; xyq_io and f_io
                 were missing from the xyzq_io folder; since June, several
                 triangulation utilities have been extended to overcome any
                 single-zone assumptions, NEQAIR_Integration has been improved,
                 and POLAR_INTERP has had an option added to generate spoked
                 body points (but see also CAPSULE_SPOKES, which had been over-
                 looked here as well).
    Nov 07 2018  Added SLOS (structured grid form of USLOS for lines of sight).
