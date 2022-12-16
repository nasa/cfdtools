# CFDTOOLS

## Overview
The CFD Utility Software Library, a.k.a. CFDTOOLS, contains nearly 30 numerical
analysis libraries and close to 100 utility applications built upon those
libraries. These utilities were developed over a roughly fifty year span to
support aerospace research and development activities at NASA Ames Research
Center (ARC). They are mostly written in Fortran 90 or 77 and are designed with
potential reuse in mind. The library also includes C translations of roughly a
dozen numerics routines in the `C_utilities` module.

The library began as the ARC Aerodynamics Division Software Library, and is
currently maintained by ARC Entry Systems and Technology Division. David
Saunders and Robert Kennelly are the primary authors, with miscellaneous
contributions from many others over the years. CFDTOOLS is expected to grow
slowly as new functionality is required.


## Feature Summary
* Single-function tools for manipulating multiblock grids and flow solutions
* Rapid searching/interpolation for structured or unstructured grids (ADT methods)
* Automated gridding for axisymmetric capsules with curvature-based discretization
* Generic optimization frameworks with various applications
* General-purpose numerics subroutines and character-manipulation subroutines

A detailed summary of each library/application is provided below.


## Points of Contacts

* jeffrey.p.hill@nasa.gov
* david.a.saunders@nasa.gov


## Build and Install
The CFDTOOLS software archive is available for download on NASA's public
[Github](https://www.github.com/nasa/cfdtools). The software may be cloned
using `git` in the typical fashion:

    git clone https://www.github.com/nasa/cfdtools.git

Once downloaded, the software is compiled and installed using CMake. The
software has no dependencies, and is known to build successfully using
relatively recent compilers from Intel and GNU on Ubuntu and SUSE Linux. The
basic installation process is as follows:

    cd <cfdtools_source_dir>
    mkdir build && cd build
    cmake -DCMAKE_INSTALL_PREFIX=<cfdtools_install_dir> ..
    cmake --build .
    cmake --install .

Upon installation, all that is required to use the CFDTOOLS applications is to
add `<cfdtools_install_dir>/bin` to the PATH environment variable. If you are
building software that links against the CFDTOOLS library using CMake, simply
add `-DCMAKE_PREFIX_PATH=<cfdtools_install_dir>` to the initial CMake
configuration command for your project. Then, in your projects `CMakeLists.txt`,
add the following:

    find_package(cfdtools REQUIRED)

    add_executable(my_program ...)
    target_link_libraries(my_program cfdtools::<libname>)

CFDTOOLS also supports direct inclusion as a project submodule using CMake's
``add_subdirectory`` command.


## License
This software is released under the [NASA Open Source Agreement Version
1.3](LICENSE.pdf).

To contribute to this software package we require that non-NASA authors assign
their rights in the contributed code to NASA. Contributor License Agreements (CLAs)
are available in the [contrib](contrib) folder.


## Notices
Copyright 2021 United States Government as represented by the Administrator of
the National Aeronautics and Space Administration.  No copyright is claimed in
the United States under Title 17, U.S. Code. All Other Rights Reserved.


## Disclaimers
__No Warranty__: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF
ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED
TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY
IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR
FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL BE ERROR
FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE
SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER, CONSTITUTE AN
ENDORSEMENT BY GOVERNMENT AGENCY OR ANY PRIOR RECIPIENT OF ANY RESULTS,
RESULTING DESIGNS, HARDWARE, SOFTWARE PRODUCTS OR ANY OTHER APPLICATIONS
RESULTING FROM USE OF THE SUBJECT SOFTWARE.  FURTHER, GOVERNMENT AGENCY
DISCLAIMS ALL WARRANTIES AND LIABILITIES REGARDING THIRD-PARTY SOFTWARE, IF
PRESENT IN THE ORIGINAL SOFTWARE, AND DISTRIBUTES IT "AS IS."

__Waiver and Indemnity__:  RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS AGAINST THE
UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY
PRIOR RECIPIENT.  IF RECIPIENT'S USE OF THE SUBJECT SOFTWARE RESULTS IN ANY
LIABILITIES, DEMANDS, DAMAGES, EXPENSES OR LOSSES ARISING FROM SUCH USE,
INCLUDING ANY DAMAGES FROM PRODUCTS BASED ON, OR RESULTING FROM, RECIPIENT'S
USE OF THE SUBJECT SOFTWARE, RECIPIENT SHALL INDEMNIFY AND HOLD HARMLESS THE
UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY
PRIOR RECIPIENT, TO THE EXTENT PERMITTED BY LAW.  RECIPIENT'S SOLE REMEDY FOR
ANY SUCH MATTER SHALL BE THE IMMEDIATE, UNILATERAL TERMINATION OF THIS
AGREEMENT.
## Library Summary

| Component                     | Description |
|-------------------------------|-------------------------------------------------------------------------|
| **aa_nurbs**                  | Supplements the DT_NURBS library, for B-spline airfoils |
| **adt**                       | Rapid search/interpolation within computational meshes |
| **blaslib**                   | Basic linear algebra subprogram utilities (vectors) |
| **blas2lib**                  | NAG versions of the Level 2 BLAS Matrix-vector routines |
| **C_utilities**               | Handful of C translations from interplib, linsys, etc. |
| **dt_nurbs**                  | (Dated) David Taylor Naval Research Center NURBS library |
| **eigenlib**                  | Just the EISPACK routines used by CAVITY_MAP rotations |
| **forsythe**                  | Original Forsythe, Malcolm & Moler library (historical) |
| **geomlib**                   | Areas, rotations, nearest points, arc lengths, etc. |
| **grid_block_utilities**      | Volume grid analogues of the surface_patch_utilities |
| **gridlib**                   | Mostly for structured 2/3-space grids; many 1-D distributions |
| **integrate**                 | Numerical quadrature, mostly 1-D, with a driving program |
| **interplib**                 | 1-D interpolation & smoothing utilities; see SMOOTH |
| **interp2d**                  | 2-D interpolation, smoothing, & intersection utilities; see SMOOTH2D |
| **interp3d**                  | 3-space interpolations, intersections, etc. |
| **intrinsics**                | Contains only trig_functions_in_degrees: g95 cosd, etc. |
| **kdtree**                    | Matthew Kennel's open source package, now an option in FLOW_INTERP |
! **lapacksubset**              ! LAPACK routines needed by symevd, symevdsolve added to linsys |
| **linsys**                    | Linear system solvers, mostly dense factorization methods |
| **numodules**                 | Miscellaneous utilities: finite differences, standard deviations, etc. |
| **obsolete**                  | Still linked to by INTEGRATE, MAXMIN, PROFILE, and SMOOTH |
| **Optimal_Interpolation**     | Alexander Barth's Kriging-like package, with extensions |
| **optlib**                    | 1-D zero/minimum finders; n-D unconstrained methods, etc. |
| **progtools**                 | String-manipulation & prompting utilities, incl. string_justify module |
| **searchlib**                 | 1-D search/sort routines; some for 2/3-D now deprecated |
| **special_functions**         | A few utilities related to erf(x), normal distributions, etc. |
| **surface_patch_utilities**   | For structured surface grids: scale/shift/transpose, etc. |
| **table_io**                  | I/O utilities for plain text tables (headers + numeric rows) |
| **tecplot_io**                | Tecplot I/O package for structured multizone data |
| **triangulation_io**          | Tecplot I/O package for triangulated surface data |
| **ugridlib**                  | Triangle and other unstructured surface/volume data utilities |
| **xyq_io**                    | PLOT2D-type multiblock structured grid I/O utilities, in xyzq_io folder |
| **xyzq_io**                   | PLOT3D-type multiblock structured grid I/O utilities; see also f_io.f90 |

## Application Summary

C = CFD, G = General Numerics, P = Programming Tool

| Program                       | Type | Description |
|-------------------------------|------|----------------------------------------------------------|
| **a2b**                       | G    | Converts a (large) ASCII file to a binary stream file.   |
| **adjust_flow**               | C    | Analogue of ADJUST_GRID for PLOT3D-format flow solutions |
| **adjust_grid**               | C    | 20-odd options: shift/scale/rotate/swap coordinates, etc. |
| **aero_coefs**                | C    | Aerodynamic force/moment coefs. for structured surface data containing p or Cp |
| **aerosurf**                  | C    | Wing/body/nacelle surface paneling; underlies shape optimization |
| **anchor**                    | G    | Applies Optimal_Interpolation package to aerodynamic or other unstructured data |
| **blayer**                    | C    | Boundary layer analysis of a structured flow solution |
| **body_point_data**           | CG   | For one body point extract time histories of flow variables as a table |
| **bsprofile**                 | C    | B-spline curve analogue of PROFILE airfoil utility |
| **bump_factors**              | C    | Compare related structured surface solutions as ratios |
| **capsule_grid**              | C    | Automated surface grid for sphere/cone; many aft body options |
| **capsule_spokes**            | C    | Supplements POLAR_INTERP by setting up spoked body points for a capsule |
| **cavity_map**                | C    | Map local cavity flow solution to unit cube; slice it, etc. |
| **ch_grid**                   | C    | Stand-alone form of the automated C-H volume gridding in wing/body design code SYN87 |
| **columnedit**                | P    | Insert/delete/replace column(s) in a dataset |
| **combine_blocks**            | C    | Combine multiblock grids/solutions; cavity block options |
| **combine_blocks_turb**       | C    | COMBINE_BLOCKS variant required for turbulent flows |
| **compare_blocks**            | C    | Compare related structured grids: identical? maximum differences? |
| **compare_flows**             | C    | Flow field analogue of COMPARE_BLOCKS (volume / surface) |
| **compare_patches**           | C    | Surface grid analogue of COMPARE_BLOCKS |
| **compress2d**                | C    | CAPSULE_GRID tool; adjust upstream boundary, wall spacing, # radial grid points, etc. |
| **cones_of_sight**            | C    | Body-normal sets of 9 discretized lines forming cones for radiation calculations |
| **cross_sections**            | C    | Slice a surface dataset (structured / unstructured) at many stations |
| **curvature**                 | C    | Evaluate curvature along X-Y curve; redistribute points according to curvature |
| **curvature3d**               | C    | X-Y-Z curve analogue of CURVATURE |
| **deconstruct**               | C    | Convert a multiblock structured surface or volume grid to unstructured form |
| **distribute**                | C    | Drives available 1-D grid-point distribution utilities |
| **extract_blocks**            | C    | Extract indicated block(s) from a multiblock grid/flow |
| **extract_blocks_2d**         | C    | Extract indicated block(s) from a multiblock 2D grid/flow |
| **extract_columns**           | GP   | Extract column(s); options to scale/shift/reformat |
| **extract_functions**         | C    | Extract selected functions efficiently from a PLOT3D function file. |
| **extract_lines**             | GP   | Extract every nth line starting from some line |
| **extract_peaks**             | CGP  | For a list of files of column data, tabulate peak values |
| **filter_lines**              | GP   | Remove non-numeric lines (excess headers?) from a table |
| **flow_interp**               | C    | Volume grid search & flow interpolation (ADT method and KDTREE method) |
| **flow_interp_2d**            | C    | 2-space analogue of FLOW_INTERP |
| **gen1d**                     | G    | Generate (x,y) datasets from a choice of y = f(x) functions |
| **gradient_based**            | CG   | Redistribute points along a curve using abs(df/ds), where s is arc length, for some function f |
| **grid_faces**                | C    | Extract grid block face(s) from a multiblock grid/flow |
| **gsmooth**                   | C    | Elliptic volume grid smoothing; option to detect bad cells |
| **gu**                        | C    | Ryan McDaniel's multi-option tool for PLOT3D grids/flows |
| **heat_shield**               | C    | Forebody shape optimization, axisymmetric option; Newtonian CL/CD |
| **hemispheres_of_sight**      | C    |Sets of discretized lines forming hemispheres for radiation calculations; see also SLOS, USLOS |
| **hyper_aero**                | C    | Driver for available impact methods, including Modified Newtonian |
| **insert_lines**              | P    | Insert lines from one file into another at some interval |
| **interp_1d**                 | G    | Interactive linear/local spline interpolation, 1-D |
| **interp_columns**            | GP   | Adaptation of REGULARIZE; variable column to treat as "x" |
| **join_grid**                 | C    | SPLIT_GRID companion for merging blocks |
| **line_grid**                 | GP   | Uniform grid between 2 points in 3-space; extrapolation OK |
| **lines_of_sight**            | C    | Body-normal | shock-normal discretized lines for radiation calculations |
| **lines_of_sight_2d**         | C    | 2-space analogue of 3-space LINES_OF_SIGHT |
| **loc**                       | P    | Count lines of code & comments in a Fortran source file |
| **los_attenuation**           | C    | Radio attenuation calculation for one line of sight; multiple options |
| **makelist**                  | P    | Make a list of integers for given i1, i2, increment |
| **maxmin**                    | G    | Analogue of SMOOTH for looking at 1st & 2nd derivatives |
| **merge_blocks**              | C    | Variant of JOIN_GRID allowing for optional function file |
| **merge_files**               | C    | Merge variables of two Tecplot files (same structured grid) |
| **merge_tables**              | GP   | Rearrange columns from any number of tables into a new table |
| **nbyn**                      | G    | Interactive solution of square n x n systems; n <= 6 |
| **neqair_data**               | C    | Converts line-of-sight data from PLOT3D form to NEQAIR form |
| **neqair_integration**        | C    | Integrates radiances from hemisphere lines of sight to estimate radiative heat flux |
| **nozzle_throat_conditions**  | C    | Axisymmetric variant of THROAT_CONDITIONS_3D |
| **npopt_driver**              | G    | Framework for a [sequence of] constrained minimization[s] |
| **outbound**                  | C    | Off-line grid/shock alignment; other redistribution options |
| **p3d2tec**                   | C    | Convert between PLOT3D grid/flow format and Tecplot form |
| **polar_interp**              | C    | Pad spoked forebody radiation data nonlinearly |
| **precision**                 | G    | Estimate objective function precision (difference table) |
| **prepare_local_analysis**    | C    |Rapid set-up of DPLR-based damage/repair/feature calculations: cavity / protruding gap-filler |
| **prepare_neqair_data**       | C    | Facilitates scripting of radiative heating calculations by NEQAIR |
| **profile**                   | C    | Airfoil geometry display/manipulation utility |
| **qnmdriver2**                | G    | Framework for one unconstrained optimization (QNMDIF2) |
| **qnmdriver3**                | G    | Framework for a sequence of unconstrained optimizations |
| **radial_interp**             | C    | Rapid volume grid|flow interpolation for new surface grid (single layer of volume grid blocks) |
| **radial_interp_2d**          | C    | 2-space RADIAL_INTERP analogue; both one layer of grid blocks |
| **rectify_grid**              | C    | Ensure right-handedness for 3-space volume or surface grid |
| **redistribute_xy**           | C    | Redistribute points in a 2- or 3-space line segment (1- or 2-sided stretching) |
| **refine_grid**               | C    | Densify or thin multiblock grid/flow data, any multipliers |
| **reflect_blocks**            | C    | Reflect 3-space grid/flow data; save reflected or both halves |
| **reflect_blocks_2d**         | C    | Reflect 2-space grid/flow data; save reflected or both halves |
| **reorder_blocks**            | C    | Reorder some or all grid blocks; optional flow file |
| **reorder_rows**              | GP   | Reorder tabular data so leading n columns vary differently |
| **reorder_segments**          | GC   | Reorder multiple 1-D Tecplot zones as 1 continuous zone |
| **reshape**                   | G    | For 2 of 2 or more columns, shift/scale/rotate/... (11 options) |
| **reshape3d**                 | G    | 3-space analogue of RESHAPE |
| **revolve_grid**              | C    | 2-space volume grid --> 3-space; may need RADIAL_INTERP next |
| **scan_grid**                 | C    | Tabulate data range, etc., by grid block; optional flow |
| **shadowgraph**               | C    | CFD density data --> shadowgraph/schlieren-like image |
| **shiftscale**                | G    | Interactive calculation of linear coefficients for transforming [a, b] -> [p, q] |
| **shock_stand_off**           | C    | For volume grid & flow field temperature, calculate shock stand-off distances |
| **slos**                      | C    | USLOS-type merge of LINES_OF_SIGHT & HEMISPHERES_OF_SIGHT (structured grid) |
| **smooth**                    | G    | Drives available 1-D interpolation & smoothing utilities |
| **smooth2d**                  | G    | Bivariate analogue of SMOOTH |
| **sort_rows**                 | GP   | Sorts the output from grep (say) or a list of file names containing numeric substrings. |
| **sort_surface_slice**        | C    | Sort Tecplot surface slice into 1 or more curves for line plots |
| **split_grid**                | C    | Split multiblock grid/flow and/or permute indices |
| **Stardust_Lines**            | C    | Radiation lines of sight for airborne or ground observation of an entry capsule |
| **Stardust_Integration**      | C    | Companion to Stardust_Lines |
| **surface_curvature**         | C    | Gaussian/mean/principal curvatures for surface/volume grid |
| **surface_diffs**             | C    | Map Tecplot surface grid 1 to grid 2; save [%]differences |
| **surface_interp**            | C    | Interpolate 3-space Tecplot or Plot3D surface data at target surface or list of data point(s) |
| **surface_interp_2d**         | C    | Interpolate 2-space Tecplot surface data at target data point(s) |
| **surface_pad**               | C    | Pad structured surface data nonlinearly: 1-D in i and/or j |
| **surface_patches**           | C    | ADJUST_GRID variant with more surface_patch_utilities options |
| **table_arithmetic**          | G    | Manipulate one or two data tables by columns |
| **template**                  | C    | Calculate block interface data and DPLR-type control file (structured volume grid) |
| **tet_interp**                | C    | Interpolate tetrahedral volume data to a structured multiblock grid |
| **thin_flow**                 | C    | THIN_GRID analogue for multiblock flow data |
| **thin_grid**                 | C    | Extract every mth/nth[/kth] point in i/j[/k] directions |
| **thin_grid_2d**              | C    | THIN_GRID analogue for 2-space; grid only |
| **throat_conditions_3d**      | C    | Boundary conditions at arc-jet nozzle: rectangular, circular, semi-elliptic |
| **trail**                     | GP   | Truncate lines and/or strip trailing blanks |
| **traj_fit**                  | C    | Nonlinear least squares fitting of f(rho(t),V(t) = C rho(t)^m V(t)^n via the QNMDRIVER2 framework |
| **traj_opt**                  | CG   | Trajectory optimization: would need an NPOPT license from Stanford University |
| **tri_to_quad**               | C    | Impose a structured surface on a Tecplot triangulation |
| **tri_to_tri**                | C    | Interpolate a single-zone Tecplot triangulation dataset to another 1+-zone target triangulation |
| **triangulation_tool**        | C    | Drives scale/shift/rotate transformations and area/volume/CM/moments of inertia calculations for a Tecplot unstructured surface or volume dataset (1+ zones) |
| **update_grid**               | C    | Replace one or more grid blocks with same-sized block(s) from other file or files |
| **upsequence**                | C    | Upsequence coarse cavity/plug/gap filler solution; impose fine local boundary flow |
| **usflowinterp**              | C    | Interpolation within an unstructured 3D Tecplottable dataset (ADT search techniques) |
| **usinterp**                  | C    | Driver for the optimal interpolation package (scattered data interpolation in n dimensions) |
| **uslos**                     | C    | Merge of LINES_OF_SIGHT & HEMISPHERES_OF_SIGHT for unstructured surfaces |
| **v2c**                       | C    | Convert grid [+ optional flow]: cell vertices to centers |
! **wedge**                     ! C    ! Convert a structured 2D grid & optional solution to wedge form for US3D or LAURA
| **wingsections**              | C    | B-spline sections + chord/thickness data --> B-spline wing |
| **xdeck**                     | P    | Ancient precursor of TRAIL: remove trailing / leading characters |
| **xline**                     | P    | Remove lines starting with a target string |


## History

* **Aug 30 2013**: Initial preparation.
* **Sep 20 2013**: Initial summary uploads to SourceForge.
* **Sep 30 2013**: Still awaiting final NASA approval to upload software.
* **Oct 17 2013**: CAPSULE_GRID document uploaded following 2.5-week government
                   shutdown hiatus.
* **Oct 22 2013**: Folders for all libraries are now set up with subroutine one-
                   liner READMEs.
* **Oct 30 2013**: Folders for all applications are now set up with README files.
* **Jan 24 2014**: NASA approval for uploading the software proper.
* **Jan 29 2014**: All the library *.tar.gz files are now in place; applications
                   still to do.
* **Feb 05 2014**: All straightforward applications have been uploaded with simple
                   build scripts.
* **Feb 10 2014**: UPDATE_GRID has been added; it and GSMOOTH now use xyzq_io, not
                   the earlier cfd_io_package.
* **Jun 30 2014**: Someone needed the DECONSTRUCT utility.  A few updates since
                   Feb 10 include solid angle utilities in ugridlib, a handful of
                   functions in special_functions, and a reworked form of the
                   HEMISPHERES_OF_SIGHT application that defines lines of sight
                   at a point on a convex body via a triangulated quadrant of a
                   hemisphere rather than via latitude/longitude discretization.
* **Jul 03 2014**: TRI_TO_TRI has been added; a triangulation_io glitch has been
                   corrected.
* **Aug 07 2014**: EXTRACT_BLAYER_DATA has been added to the BLAYER folder.
* **Aug 14 2014**: Optimal_Interpolation now contains the package for treating
                   unstructured data as modified at NASA Ames from Alexander
                   Barth's original, following NASA approval of its added
                   BSD 2-Clause License.
* **Sep 05 2014**: FILTER_LINES and MERGE_TABLES have been added.
* **Sep 24 2014**: Traj_Fit has been added.  It curve-fits heat flux or pressure
                   pulse data using f(rho(t),V(t)) = C rho(t)**m V(t)**n within
                   the QNMDRIVER2 framework, as might be needed for TPS sizing
                   at an entry vehicle body point from selected CFD data points.
* **Oct 03 2014**: The table_io module used by Traj_Fit has been updated.
* **Oct 22 2014**: The triangulation_io package has been extended with options
                   for calculation of area, volume, CM, and moments of inertia.
                   See the new TRIANGULATION_TOOL for driving these options
                   along with shift/scale/rotate options.
* **Nov 17 2014**: Unstructured volume analogues of unstructured surface
                   utilities have been added to triangulation_io and are
                   driven by TRIANGULATION_TOOL.
* **Mar 18 2015**: An apparently inconsequential glitch in BLAYER has been
                   remedied.  This version can also read [unformatted] volume
                   datasets to speed processing of large-grid cases.
* **Mar 13 2015**: OUTBOUND now has an option to applied a CONSTANT margin
                   as a way of moving the forebody boundary a little while
                   barely affecting the wake boundary.
* **Apr 15 2015**: The ADT rapid search package (/adt) now has 2-space multi-
                   block curve analogues of the 3-space structured multiblock
                   surface grid utilities.  LINES_OF_SIGHT_2D uses these to
                   deal with reflected axisymmetric 2-space volume grid.
                   SURFACE_INTERP_2D has been added, using the same 2-space
                   ADT utilities.
* **Aug 23 2016**: TET_INTERP has been added to interpolate tetrahedral volume
                   data to a structured multiblock grid.
* **Mar 17 2017**: SORT_ROWS (added Feb 14) has been extended to allow more
                   than one numeric string in a text column, and to allow
                   output in descending order if desired.
                   Other recent refinements: COLUMNEDIT now allows extraction
                   of more than one column at a time; the table_io module has
                   an updated table_io_read_alpha utility prompted by SORT_ROWS
                   that handles the case of all lines/rows treated as text and
                   having the same number of columns (with no header records
                   assumed).  Previously, all rows were interpreted as header
                   records, meaning the table appeared empty.
* **Jun 23 2017**: NEQAIR_DATA, NEQAIR_Integration, and PREPARE_NEQAIR_DATA
                   have been added to assist radiation calculations with
                   NEQAIR.  NEQAIR itself should be obtained through the
                   Commercial Technology Office at NASA Ames Research Center.
* **Jan 04 2018**: Added A2B utility.  Since June, minor changes have been
                   made to GU, BLAYER, LINES_OF_SIGHT, SHOCK_STAND_OFF, and
                   OUTBOUND.
* **Jun 06 2018**: Added Aero_Coefs application (supplement to DPLR's Postflow).
                   Since January, updates have been made to geomlib, gridlib,
                   progtools, searchlib, CAPSULE_GRID, CONVERTQ, GU, HEMISPHERES_-
                   OF_SIGHT, NEQAIR_Integration, RADIAL_INTERP, and RESHAPE3D.
* **Sep 19 2018**: Added EXTRACT_BLOCKS_2D and EXTRACT_FUNCTIONS; xyq_io and f_io
                   were missing from the xyzq_io folder; since June, several
                   triangulation utilities have been extended to overcome any
                   single-zone assumptions, NEQAIR_Integration has been improved,
                   and POLAR_INTERP has had an option added to generate spoked
                   body points (but see also CAPSULE_SPOKES, which had been over-
                   looked here as well).
* **Oct 23 2022**: Added DECONSTRUCT (3D structured grid/flow --> unstructured form).
* **Nov 07 2018**: Added SLOS (structured grid form of USLOS for lines of sight).
* **Mar 07 2021**: Release updated CMake build system. Move repository to NASA's
                   public [GitHub](https://www.github.com/nasa/cfdtools). Various
                   small improvements to `radial_interp`, `refine`, `capsule_grid`,
                   and `prepare_neqair_data`.
* **Nov 07 2021**: Added BODY_POINT_DATA and miscellaneous updates to repository.
* **Nov 20 2021**: Re-reviewed and approved for open source release by Ames Research
                   center under the NASA Open Source License v1.3. Added updated
                   license file and contributor license agreements.
* **May 19 2022**: Added USINTERP (driver for the optimal interpolation package),
                   USFLOWINTERP and USREFLECT for operating on Tecplottable outputs
                   from US3D, multiple new routines in geomlib and ugridlib needed
                   by these applications, updates to USLOS, triangulation_io and
                   triangulation_tool, and minor updates to gsmooth, tri_to_tri,
                   and tri_to_quad.  An lapacksubset library has also been added
                   because of changes to lib/optinterp/optimal_interpolation.f90.
                   It is needed by the symmetric eigenvalue routines added to linsys.
* **Dec 05 2022**: Added WEDGE (2D structured --> 3D wedge grid & optional solution).
