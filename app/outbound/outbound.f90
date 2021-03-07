!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program outbound
!
!  Introduction:
!
!     OUTBOUND applies a choice of methods for tailoring a hypersonic flow grid.
!     Tailoring refers mainly to aligning the outer boundary with the shock.
!     Techniques demonstrated here are expected to be incorporated in the DPLR
!     flow solver for periodic adaptation of the grid as the solution converges.
!     Adjusting the radial grid point distributions is an implied requirement.
!     Again, several options are provided.
!
!     SAGE performs these functions along with further solution adaptations, and
!     the LAURA flow solver is able to tailor the grid during the flow solution.
!     The intent here is to implement ideas from SAGE and LAURA and elsewhere in
!     a form that lends itself to installation in DPLR.
!
!  One-Block-at-a-Time Strategy:
!
!     Initially, it is hoped that adaptations can be performed on grid blocks
!     in succession, independently of other grid blocks.  Even if grid blocks
!     are split at DPLR run time, DPLR possesses the mechanism needed to gather
!     each block, operate on it, and scatter the results back to the relevant
!     processors.
!
!     [LATER REALIZATION:  Actually, DPLR cannot afford to process more than
!     one split block on one processor (for memory reasons).  Thus, the only
!     parts of the tailoring scheme that must treat a whole block at a time are
!     the smoothing of the outer boundary changes in radial arc length, and the
!     smoothing of the cell-Reynolds-number-based wall increments.  Both of
!     these exceptions involve surface data only.
!
!     In this stand-alone driver for the tailoring scheme, block splitting is
!     not in the picture, so all steps operate on one whole block at a time.
!     The implementation in DPLR will perform everything in parallel on the
!     split blocks except for the smoothing of surface data.]
!
!     Thus, the grid block is considered the natural unit for processing here.
!     Common boundaries must therefore be treated exactly the same way.  This
!     is not possible if blocks don't necessarily have entire faces in common.
!     Therefore, it is highly recommended that OUTBOUND be applied only to grids
!     that are split anywhere a neighboring block has a corner (NO "SUBFACING").
!
!     Preserving smoothness in the outer boundary is one of the challenges.
!     The one-block-at-a-time strategy means block corner points must be omitted
!     from any smoothing:  interior edge points can be smoothed with edge data
!     (only), then interior surface points can be smoothed, but no guarantee of
!     smoothness across block boundaries is possible with this strategy.
!
!  Smoothing Strategy:
!
!     Rather than smoothing (x,y,z) directly, the tailored outer boundary arc
!     lengths were initially smoothed as functions of two variables - the arc
!     lengths in the two surface index directions of the untailored boundary.
!     Use of normalized arc lengths loses information from outer patches that
!     are far from square.  Instead, the bidirectional smoothing uses real
!     space arc-length-based weighting of the contributions from the two
!     directions, combined with decay of those weights towards the edges.
!
!     Smoothing radial arc lengths, though, is always affected by geometric
!     features of the inner boundary.  Therefore:
!
!     Better idea:  smooth the distances of the adjusted boundary from the
!     original boundary.  With repeated tailoring, this smoothing should
!     approach a steady state, assuming the initial boundary is reasonably
!     smooth and stays that way.
!
!     Note that smoothing full arc lengths (a la SAGE) would not leave a
!     perfectly tailored grid as is, whereas smoothing CHANGES in arc length
!     (ideally all zero at convergence) would.
!
!     Nevertheless, for situations where the initial outer boundary is not
!     smooth, this version of OUTBOUND has an option to perform such smoothing.
!     Rather than attempting to interpolate the solution and proceed with
!     tailoring, the program terminates after such smoothing.  [DPLR can now
!     perform either or BOTH types of smoothing.  Most recently, DPLR can now
!     smooth surface data across block boundaries, thanks to Matt Bartkowicz's
!     work that has not been incorporated in OUTBOUND.  (It would conflict with
!     the straightforward processing of one block at a time.)  Matt's smoothing
!     is index-based rather than surface-arc-length based, but this appears to
!     be a non-issue.  Note that after several alignments, the initial estimate
!     of the shock location before any smoothing should be good, with little
!     noise, so fewer smoothing iterations are recommended - e.g., 10-15 for the
!     final pass, depending on grid resolution.]
!
!  Other stipulations:
!
!     >  Adaptations are strictly along the radial grid lines of the original
!        grid except that any extrapolation via added_margin is strictly linear.
!     >  If the grid contains more than a single layer of blocks, only the outer
!        layer of blocks is normally adapted (i.e., the blocks with a freestream
!        boundary are tailored to the shock, or possibly left alone in places
!        determined to be in the wake).  However, the retrofitted option to fix
!        a block outer boundary and apply the basic tailoring method to k planes
!        1 : mk < nk does allow an outer block layer to remain untouched.  This
!        requires the retrofitted block-blanking option.
!     >  The DPLR boundary condition control data are employed to treat grids
!        with arbitrary indexing.  For instance, a face with a BC of 1 means it
!        is at the free-stream boundary, while its opposite face may or may not
!        be at a solid boundary.
!     >  All files are PLOT3D-type grid and functions files, multiblock, 3-D,
!        (or 2-D now), formatted or not.
!     >  The function file is vertex-centered with one or more functions.  The
!        function to be used for shock detection - say Mach number - may be in
!        any specified location, but if cell-Reynolds-number-based spacing is
!        specified, the first three functions should be DENSITY, SPEED OF SOUND,
!        and VISCOSITY in that order.
!     >  Output of an interpolated flow solution is not attempted here.  To be
!        of any use to DPLR, the function file would need to contain the state
!        variables at the cell centers, and these would not include the Mach
!        number.  [Michael Wright points out that an approximation to Mach
!        number could be made, but standard practice would not be to use an
!        interpolated solution.  Halo cell issues are best avoided for now.]
!     >  Initial coding employs conventional (i,j,k) indexing for compatibility
!        with the XYZQ_IO package.  DPLR works with (k,j,i) ordering.  This
!        conflict is handled by swapping conventions within OUTBOUND, not DPLR.
!     >  The argument-driven routine that tailors one grid block cannot change
!        the size of that block because it updates the block in-place for the
!        convenience of DPLR.  This does not preclude the stand-alone driver
!        from imposing new grid point counts in a routine that will not be used
!        by DPLR.
!
!  Option to Compress an Initial Hyperbolic Volume Grid:
!
!     Shock edge method = -1 (along with a fudge factor at the bottom of the
!     control file that's not required for an optional free stream value) is
!     available in this version as an attempt to reduce the need for double
!     smoothing when DPLR is started on a hyperbolic volume grid commonly used
!     for capsule forebodies.  (Cases with aft bodies may also be compressed,
!     although the wake region may grow in length as a consequence of allowing
!     a margin for the shock, and there is no way at present to shorten wake
!     regions.)
!
!     Ancillary program COMPRESSION_DATA can determine the transformation from
!     an initial grid to an existing aligned grid.  The transformation consists
!     of scale factors applied to the radial grid line arc lengths, saved as a
!     multiblock 3-space (r, theta, z - zu) dataset, where r is the distance of
!     the initial outer boundary point U from the point furthest upstream (not
!     its perpendicular distance to the line through that point parallel to the
!     X axis, which would require a steady increase of the initial grid's outer
!     radius).  An input fudge factor is applied to the scale factors.
!
!     Running OUTBOUND with fudge factor = 1.0 and the initial grid used when
!     the compression dataset was generated would reconstruct the associated
!     compressed grid.  With fudge factor > 1.0, the radial lines of the output
!     grid will be longer than for the original compressed grid. The compression
!     can be varied this way - use trial and error.  Preserving convexity in the
!     resulting outer boundary is not guaranteed, but the hope is that some
!     compression will be preferable to none when a new geometry is tackled.
!     Compression data from a related geometry should be applicable initially.
!     The geometry dimensions can change, as can the surface grid topology,
!     because the compression dataset (compression_data.g and .f) is adjusted by
!     OUTBOUND for the maximum radius of the input grid, and it is searched
!     efficiently for scale factors in (r, theta, dz) space using the structured
!     surface form of an Alternating Digital Tree search package.  As with DPLR,
!     the effects of non-zero angle of attack need to be kept in mind.  Best
!     results should be obtained using compression data from a similar angle of
!     attack case, particularly the same Alpha for similar/related geometries.
!
!     The standard ds1 >= 0 options for radial distribution method 2 are
!     provided (existing wall spacing or specified spacing, one- or two-sided
!     stretching).
!
!  Option to Scale the Volume Grid (SAGE-like, two variations):
!
!     Shock edge method = -2 (along with a scale factor at the bottom of the
!     control file that's not required for an optional free stream value) is
!     available for growing or shrinking the outer boundary in a way that
!     differs from all other methods:  it does not retain the existing radial
!     lines unless they happen to be perfectly straight.  Scale factors close
!     to 1 (smaller or larger) are recommended, although the added_margin
!     capability is likely to be a better choice.  Two variations are offered:
!
!     scale_factor > 0 (constant everywhere):
!
!             x(k) <-- x(1) + scale_factor * (x(k) - x(1))
!             y(k) <-- y(1) + scale_factor * (y(k) - y(1))
!             z(k) <-- z(1) + scale_factor * (z(k) - z(1))
!
!     scale_factor < 0 (adjusted scaling as follows):
!
!             dx        = x(k) - x(1)             [and likewise for dy, dz]
!             dstraight = sqrt (dx^2 + dy^2 + dz^2)
!             factor    = -scale_factor * (dstraight / sold(k))
!             xnew(k)   = x(1) + factor * dx      [and likewise for ynew, znew]
!
!     The second variation slightly improves the effect of curving radial lines,
!     but not enough to help much if |scale_factor| is large.
!
!  OUTBOUND Outline:
!
!     >  Open and read the control file (standard input - format shown below)
!
!     >  Open and read the DPLR-type file of control data including BCs
!
!     >  For each block of the input grid:
!
!        >  Read the grid block and the corresponding function file block
!
!        >  Determine outer boundary face indices; do nothing if no BC = 1
!
!        >  For each outer boundary point:
!
!           >  Locate the preliminary new outer boundary (shock detection)
!
!        >  Smooth all 4 edges of the new locations (1-D method, corners fixed)
!
!        >  Smooth interior points of the new outer boundary (2-D surf. method)
!
!        >  If cell Reynolds number is being used to define radial increments
!           at the wall, generate them all for this block and smooth them.
!
!        >  For each radial line:
!
!           >  Redistribute the radial line along the original arc
!
!        >  Update the original block in-place.
!
!        >  If specified, adjust the point counts [not applicable to DPLR]
!
!        >  Output the tailored grid block
!
!  Control File Format (Standard Input)
!
!     OUTBOUND controls for case xxx
!     ---------------------- INPUT VOLUME GRID ---------------------------------
!     baseline.g      Initial grid
!     T               Formatted? [T|F]
!     --------------------- INPUT FUNCTION FILE --------------------------------
!     baseline.f      Flow quantity(s), vertex-centered: [rho, a, mu,] M|p|rho
!     T               Formatted? [T|F]
!     1               Index of function to be used for tailoring
!     --------------- INPUT DPLR CONTROL FILE WITH BCS -------------------------
!     dplr.inputs     For ibc(6,nblocks)
!     --------------------- TAILORED VOLUME GRID -------------------------------
!     tailored.g      Output volume grid
!     T               Formatted? [T|F]
!
!     --------------- OUTER BOUNDARY TAILORING CONTROLS ------------------------
!     1               Shock edge method|0; 1|4: Mach; 2: p,rho,T; 3: flow grad.
!     0               # smoothing iters. applied to flow var. [method 3 only?]
!     0.95            Freestrm. multiplier or tolerance; depends on edge method
!     4               Multiple of local grid spacing for outer margin
!     1.50            Multiple of local spacing to limit smoothing changes
!     30              # smoothing iterations applied to outer boundary
!     2.              Additional margin (times outer spacing, possibly < 0)
!     0               nk_margin >= 0; plane nk - nk_margin is treated as for 0
!
!     --------------- RADIAL DISTRIBUTION CONTROLS -----------------------------
!     1               1: Re_cell; 2: given ds1; 3 (method_edge = 0): use |dT/ds|
!     0               nradial > 0 allows changing the no. of radial grid points
!
!     ............... Radial method 1 controls: Re-cell-based ds1 ..............
!     1.              Nominal cell Reynolds # to impose at the wall
!     0.000005        Minimum wall spacing (before smoothing); < 0 => % arc
!     0.000015        Maximum   "     "        "        "        "       "
!     10              ng where pts. 1:ng are geometric; ng <= 2 => pure Vinokur
!     1.05            Geometric growth rate if ng > 2
!     0.25            Factor applied to 1-sided ds2 to give 2-sided ds2, or 0.
!                     0. means preserve outer spacing.
!     ............... Radial method 2 controls: ds1 specified ..................
!     0.              ds1 > 0 => const.; 0 => existing; < 0 => % arc; 999. => du
!     10              ng where pts. 1:ng are geometric; ng <= 2 => pure Vinokur
!     1.05            Geometric growth rate if ng > 2
!     0.25            Factor applied to 1-sided ds2 to give 2-sided ds2
!
!     ............... Radial method 3: recluster outer pts; use edge method = 0)
!     0.5  0.1  0.5   Outer and blend fractions of nk, and shape fn. exponent
!
!     --------------- PLOTTABLE PROFILE CONTROLS + OPTIONAL INPUTS -------------
!     8               Block number to extract radial profiles from
!     1 17 1 1        Indices defining surface patch row(s)|column(s) to extract
!     17.88           Optional free-stream; < 0 => local peak; absent => block 1
!                     (also used as a scale factor if shock edge method = -1|-2,
!                     or as a constant outer margin (+ve or -ve) if present,
!                     nonzero, and shock edge method = 0)
!
!
!  Ancillary Control File ('outbound.inp.2'):
!
!     This optional file can be used to blank blocks (suppress all changes).
!     Enter any reasonable list of blocks on the first line.  E.g.,
!
!        11:16 or 11-16   would both expand to 11, 12, 13, 14, 15, 16
!        10 12 14:20      or any other such intelligible list, in any order
!
!     This first line can be empty, meaning no blocks are suppressed.
!     If the file is not present, no blocks are suppressed.
!
!
!  Further Notes:
!
!     (1)  Shock edge detection methods (method_edge):
!
!            -2 => Scaling options - see description above
!            -1 => Compression option - see description above
!             0 => a: if method_radial /= 3:
!                     smooth radial arcs directly (no further tailoring):
!                     if a function file is present, it may still be used with
!                     method_radial = 1 (cell Re # option); otherwise, only the
!                     file header will be read; use file name 'none' if no such
!                     file is available.
!                  b: if method_radial == 3:
!                     no edge detection -- only redistribution of the outer
!                     portion of each radial line so as to cluster towards the
!                     shock (and possibly other high-temperature-gradient flow
!                     features) without touching the boundary layer region;
!                     the function file is assumed to have translational T as
!                     the first (and probably only) function.
!             1 => Mach number: some fraction of freestream (or of local peak).
!             2 => LAURA-type significant increase in pressure, density or T.
!             3 => SAGE-type relative flow gradient change.
!             4 => As for 1 but search outwards, not inwards, & omit wake test.
!          Corresponding fs_scale:
!             1 => fraction < 1 such as 0.95 applied to the freestream value.
!             2 => multiplier applied to the freestream value as in LAURA.
!             3 => fraction applied to peak flow gradient as in SAGE.
!             4 => as for 1.
!
!     (2)  Safety margin (ds_local_multiple, added_margin, and constant_margin):
!
!             ds_local_multiple = multiple of local spacing at estimated edge;
!                                 e.g., 4.*local_ds, added to edge location
!                                 before smoothing.
!             smoothing_limiter * local_ds provides upper & lower bounds on how
!                                 much smoothing of delta arcs can change them,
!                                 given that the initial edge estimate should
!                                 be in the right neighborhood.
!             smoothing_limiter * outer ds bounds direct smoothing of the
!                                 outer boundary similarly.
!             added_margin      = multiple of outermost spacing (possibly < 0)
!             constant_margin   = constant absolute margin entered (with
!                                 shock edge method = 0) in place of the
!                                 optional input for free-stream Mach number
!                                 used when shock edge method /= 0.
!
!     (3)  ds1 input (method_radial == 2):
!
!             > 0.   => apply this constant spacing at the wall everywhere.
!             = 0.   => retain the existing wall spacing everywhere.
!             < 0.   => apply this percentage of total arc length at the wall.
!             = 999. => apply uniform spacing along all radial lines.
!
!     (4)  Gradient-based clustering towards the shock (method_radial = 3):
!
!             Use method_edge = 0 for this option.  Some outer portion (say
!             half the number of radial points) is redistributed to resolve the
!             shock better as is desirable for uncoupled radiative heating
!             calculations.  The outer boundary location is not changed; nor
!             is the boundary layer portion (innermost ~half of the points).
!             Inputs outer_fraction and blend_fraction are applied to nk to
!             control what fraction of radial points are changed (leaving the
!             boundary-layer points undisturbed) and what fraction either side
!             of that split is used to blend the spacing (2-sided Vinokur).
!             0.5 and 0.1 should be reasonable inputs.
!
!             No other redistribution options are permitted at the same time.
!             The outer margin should be somewhat larger than normal so that
!             there will still be at least 2 points at free-stream conditions
!             after the redistribution.
!
!     (5)  nk_margin > 0 and method_radial /= 3:
!
!             Any unblanked block will have its outer boundary left intact,
!             presumably for compatibility with another block layer such as
!             the nozzle block(s) in an arc-jet simulation.  The normal align-
!             ment algorithm will be applied to k planes 1 : nk - nk_margin.
!             The points between mk = nk - nk_margin and nk will be blended at
!             k = mk and reuse the original spacing for k = nk - 1 : nk.
!
!  Sponsor:
!
!     Reacting Flow Environments Branch, NASA Ames Research Center,
!     in particular, Michael Wright, main author of DPLR.
!
!  History:
!
!     05/27/05  DAS/MJW/SSY  Initial design (one whole block at a time).
!     06/06/05  DAS/SSY      Initial testing (Mach number profiles; retaining
!                            relative radial distributions suffices for now).
!     06/10/05   "   "       Pressure and density-based shock detection seems
!                            too difficult, at least for the Shuttle, in spite
!                            of heuristics to try and cope with unpredictable
!                            flow profiles.
!     06/24/05   "   "       (After a hiatus:) higher fractions of Mach number
!                            look promising for Shuttle Mach 9 case (0.90+).
!     06/28/05   "   "       0.90 x Minf + 0.20 x arc margin still tends to
!                            underestimate the lee-side boundary and over-
!                            estimate most of the wind-side boundary.
!                            Try backing up a specified multiple of the local
!                            grid spacing rather than of total arc length.
!                            Use a fraction of total length as a cap.
!     07/10/05   "   "       Added Vinokur-type radial redistribution options.
!                            Smooth distances to the old outer boundary, not
!                            radial distances from the OML.
!     07/14/05   "   "       Handle more than one flow quantity in the function
!                            file, to allow for cell-Reynolds-number-based
!                            spacing at the wall.
!     08/09/05   "   "       Reorganization to drive a remodularized form of
!                            the original TAILOR_BLOCK in which the smoothing
!                            steps are separated for use by DPLR on surface
!                            data for a whole block while the remaining steps
!                            are suited to parallel operation on split blocks.
!     08/17/05   "   "       Accommodated DPLR's (k,j,i) indexing but left the
!                            (i,j,k) convention alone at the higher level for
!                            compatibility with the I/O utilities.
!     11/23/05   "   "       Added "added_margin" at Mike Wright's request;
!                            added an option to smooth a wiggly outer boundary
!                            (only - no other tailoring until the solution has
!                            been reconverged) - use method_edge = 0.
!     11/28/05   "   "       Disallowing the cell Re # option if method_edge = 0
!                            was misguided.  Allow for that and also for not
!                            having any relevant function file.
!     12/22/05   DAS         Ensure linear extrapolation if the outer boundary
!                            expands - parametric cubics can't be trusted.
!     01/31/06    "          10 smoothing iterations for cell-Reynolds # based
!                            increments may be too few.  It should be an input,
!                            but bump it up to 30 for now.
!     02/04/06    "          The added_margin input should be more useful if it
!                            is applied as a multiple of the outermost spacing.
!     04/24/06    "          Added option to change the number of radial points.
!                            This is not appropriate for the DPLR installation,
!                            and is appropriate here for standard grids only
!                            (where the k direction is off the wall for all
!                            blocks).
!     05/05/06    "          Added option to impose % arc length increments
!                            at the wall (possible alternative to cell Reynolds
!                            number for difficult cases).
!     05/06/06    "          The limits on the cell-Re-based wall increments
!                            may now be relative to arc length (-% inputs).
!                            Both inputs should have the same sign.  Smoothing
!                            occurs AFTER the limits are imposed.
!     12/17/06    "          Option to constrain the smoothing, which has been
!                            observed (in DPLR) to cause the shock to be hit
!                            by the boundary in front of the wing leading edge.
!     12/18/06    "          Refinement to ds_local in shock location routine.
!                            Long-time bug: method 2 radial controls were
!                            clobbering the values read for method 1!
!     01/13/07    "          Chun Tang had trouble with a Mach 2.5 case: the
!                            outflow boundary Mach number can exceed the free
!                            stream value.  Therefore, treat boundary points not
!                            at (essentially) free-stream as for those already
!                            considered in the wake, where Minf * 0.97 (or
!                            whatever) cannot be found: leave the boundary
!                            location alone.
!                            Optional final control input allows for overriding
!                            the free-stream value to use if the first point
!                            from block 1 is not in the free stream.
!     02/22/07    "          No-outer-boundary blocks were not being transcribed
!                            to the output file.  Safeguarded EXPDIS4 from an
!                            initial "smallest" increment that is on the wrong
!                            side of the uniform increment.
!     09/05/07  Todd White   Edge method 4 introduced to search from the wall
!                            and use local peaks instead of a single f.s. Mach.
!                            The test for being in the wake is also by-passed.
!     10/04/07  DAS/TRW      Made Todd's option permanent, with the added option
!                            to enter a negative free-stream value to invoke use
!                            of the local peak for either edge method 1 or 4.
!                            Method 4 now differs from 1 both in its search
!                            direction and in its lack of test for wake blocks
!                            where the outermost flow value is at the specified
!                            free stream. [Unplanned-for application of OUTBOUND
!                            to nozzle flows is what prompted Todd's experiment-
!                            ation.  It may be that no iso-surface behaves well
!                            everywhere in such flows, but use of local peak
!                            Mach numbers is now an option anyway. Note that the
!                            distinction from edge method 2 is a bit fuzzy now.]
!     10/05/07  DAS          Added the block-blanking and nk_margin options with
!                            an optional outbound.inp.2 file.
!     12/17/07   "           Allow nsmoothg < 0, which means entire edges are
!                            not touched during smoothing (not just corners).
!     10/15/09   "           Replaced EXPDIS4 with EXPDIS5 (which switches to
!                            geometric distributions if trouble arises).  The
!                            safeguarding of EXPDIS4 is no longer needed.
!     12/24/09   "           For method_radial = 2, ds1 = 999. means apply
!                            uniform distributions along all radial lines, as
!                            may be helpful for elliptic smoothing purposes.
!     03/24/11- D. Saunders  Compression option (shock edge method = -1 along
!     03/29/11  ERC, Inc./   with fudge factor on the existling last line of
!               NASA ARC     the control file) implemented as an attempt to
!                            improve widely used initial hyperbolic volume
!                            grids.  See description above.
!                            Also:  SAGE-like scaling option (shock edge method
!                            = -2, two variations) - see description above.
!     04/11/11   "    "      Fudge factor = 1.0 now corresponds to reproducing
!                            the compressed grid from the initial grid used to
!                            generate the compression data.  Fudge factor > 1.0
!                            => longer radial lines/less compression.
!     04/12/11   "    "      The compression transformation is now stored as a
!                            3-space (r, theta, dz) dataset, not the original
!                            2-space (r, theta) form, which suffered from rare
!                            aft-body anomalies because of extreme skewness.
!     08/06/13   "    "      All ADT variants are now in a single module with
!                            generic build_adt and search_adt interfaces.
!     12/12/13   "    "      2-D grids are now an option, in anticipation of
!                            a T-gradient redistribution option intended for
!                            better shock resolution, with uncoupled radiative
!                            heating calculations by NEQAIR in mind.
!     12/18/13   "    "      Installed GRADDIS[3D]2 for the shock-clustering
!                            option (method_edge = 0 + method_radial = 3) as
!                            described above.
!     03/13/15   "    "      A case where the shock was too close to the
!                            forebody prompted a way of growing the forebody
!                            outer boundary without also extending the wake
!                            too much (as gmargin applied to local outer ds
!                            tends to): use method_edge = 0 along with a
!                            constant/absolute margin in place of the optional
!                            input for Mach (free) used when method_edge /= 0.
!     04/11/16   "    "      Dinesh needed a way of changing wall spacing for
!                            some inner blocks while preserving the outer
!                            spacing.  (The inner blocks had to be extracted
!                            to be operated on with OUTBOUND.)  Use radial
!                            method 2 with ds2_fraction = 0. and fixed d1.
!  Authors:
!
!     David Saunders, ELORET Corp./NASA Ames Research Center, Moffett Field, CA
!                     Now ERC Inc./NASA ARC
!     Seokkwan Yoon,  Applications Branch, NAS Division, NASA ARC.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Derived data type and PLOT3D-type I/O module:
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use grid_block_structure  ! Must be compatible with the xy[z]q_io packages
   use xyq_io_module         ! 2-space form of PLOT3D I/O package
   use xyzq_io_module        ! Compile these three first

   implicit none

!  Local constants:
!  !!!!!!!!!!!!!!!!

   integer, parameter ::   &
      lung   = 1,          & ! Input volume grid
      lunf   = 2,          & ! Input function file
      lunbcs = 3,          & ! Input DPLR control file with BCs
      lunctl = 5,          & ! Control file (standard input/command line)
      luncrt = 6,          & ! Screen diagnostics
      lungrd = 7             ! Output volume grid

   integer, parameter ::   &
      nsmooth_d1 = 30        ! For cell Reynolds # deltas

   real, parameter ::      &
      zero   = 0.

   logical, parameter ::   &
      false  = .false.,    &
      true   = .true.

   character, parameter :: &
      format * 11 = 'unformatted'

!  Local variables:
!  !!!!!!!!!!!!!!!!

   integer :: &
      i, i1, iadapt, ib, ib_profiles, ie_profiles, ifree, is_profiles, ios,    &
      iq, iwall, j, je_profiles, js_profiles, k, k1_face, kb, kpeak_min,       &
      method_edge, method_radial, mi, mj, mk, n, n_too_close, nblank, nblk,    &
      nblocks, ndim, nf, ngeometric, ni, nj, nk, nk_margin, nlines_tailored,   &
      nlines_to_average, npts, nradial, nskip, nsmoothf, nsmoothg

   integer :: &
      icontrols(20) ! More than enough

   integer, allocatable, dimension (:) :: &
      iblock

   integer, allocatable, dimension (:,:) :: &
      ibc

   real :: &
      added_margin, avg_of_fractions, blend_fraction, cell_Re, &
      constant_margin, d1, ds, ds_local_multiple, ds1, ds1_bl, ds1_bu, &
      ds1_max, ds1_min, ds2_fraction, fpeak_min, fs_scale, fs_value, &
      growth_rate, outer_fraction, power, smoothing_limiter, sum_of_fractions, &
      worst_fraction

   real :: &
      arc_added, arc_length, arc_limit, arc_unsmoothed, rcontrols(20)

   real :: &
      rmaxh, xyz_at_xmin(3)

   real, allocatable, dimension (:,:) :: &
      d1_cell_Re, original_arcs, radial_change, stored_deltas, u, v, xs, ys, zs

   logical :: &
      flow_needed, formatted, formatted_out, left_half, right_half

   logical, allocatable, dimension (:) :: &
      suppress

   character :: &
      filename * 80

   type (grid_type), pointer, dimension(:) :: &
      grid, new_grid

!  Explicit interface to allow array allocation in a called procedure:
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   interface
      subroutine swap_ijk_conventions (mode, nf, ni, nj, nk, x, y, z, f)

      implicit none
      integer, intent (in)               :: mode, nf, ni, nj, nk
      real, pointer, dimension   (:,:,:) :: x, y, z
      real, pointer, dimension (:,:,:,:) :: f

      end subroutine swap_ijk_conventions
   end interface

!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                   Execution:
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Echo the control file to the output log:

   do ! Until EOF
      read (lunctl, '(a)', iostat=ios) filename ! Handy buffer
      if (ios < 0) exit
      write (luncrt, '(1x, a)') trim (filename)
   end do

   rewind (lunctl)

!  Deal with the input grid:
!  -------------------------

   read (lunctl, *)                  ! Case description
   read (lunctl, *)                  ! Header
   read (lunctl, *) filename         ! Input grid name
   read (lunctl, *) formatted

   i1 = 1; if (formatted) i1 = 3

   open (lung, file=filename, status='old', form=format(i1:11), iostat=ios)

   if (ios /= 0) then
      write (luncrt, '(/, 2a)') ' Unable to open the input grid file: ',       &
         trim (filename)
      go to 99
   end if

!  2-D or 3-D?

   close (lung)  ! Because determine_grid_form hasn't been made use of

   call determine_grid_dim (filename, lung, formatted, ndim, ios)
   if (ios /= 0) go to 99

   open (lung, file=filename, status='old', form=format(i1:11), iostat=ios)

   if (ndim == 2) then
      call xy_header_io (1, lung, formatted, nblocks, grid, ios)
      grid(:)%nk = 1
   else
      call xyz_header_io (1, lung, formatted, nblocks, grid, ios)
   end if
   if (ios /= 0) go to 99

!  Deal with the input function file:
!  ----------------------------------

!  This won't be used if method_edge = 0 unless method_radial = 1 or 3.

   read (lunctl, *)                       ! Header
   read (lunctl, *) filename              ! Input function file name
   read (lunctl, *) ! formatted           ! I/O package assumes same as grid
   read (lunctl, *) iq                    ! Flow variable to use for tailoring

   flow_needed = filename(1:4) /= 'none'  ! This may be updated below

   if (flow_needed) then

      open (lunf, file=filename, status='old', form=format(i1:11), iostat=ios)

      if (ios /= 0) then
         write (luncrt, '(/, 2a)') ' Cannot open the input function file: ',   &
            trim (filename)
         go to 99
      end if

      if (ndim == 2) then
         call q_header_io_2d (1, lunf, formatted, nblocks, nf, grid, ios)
         grid(:)%mk = 1
      else
         call q_header_io    (1, lunf, formatted, nblocks, nf, grid, ios)
      end if
      if (ios /= 0) go to 99

      do ib = 1, nblocks ! The function file should be vertex-centered
         if (grid(ib)%mi /= grid(ib)%ni .or. &
             grid(ib)%mj /= grid(ib)%nj .or. &
             grid(ib)%mk /= grid(ib)%nk) ios = 1
      end do

      if (ios /= 0) then
         write (luncrt, '(/, a)') ' The function file must match grid vertices.'
         go to 99
      end if

   else  ! Avoid passing undefined arguments, etc.:

      nf = 0;  iq = 0
      do ib = 1, nblocks
         grid(ib)%mi = grid(ib)%ni
         grid(ib)%mj = grid(ib)%nj
         grid(ib)%mk = grid(ib)%nk
      end do

   end if

!  Deal with the DPLR control file:
!  --------------------------------

   read (lunctl, *)
   read (lunctl, *) filename

   open (lunbcs, file=filename, status='old', iostat=ios)

   if (ios /= 0) then
      write (luncrt, '(/, 2a)') ' Unable to open the DPLR control file: ',     &
         trim (filename)
      go to 99
   end if

!  Scan it for the number of grid blocks:

   do
      read (lunbcs, '(a)') filename
      if (index (filename, 'nblk') > 0) exit
   end do
   read (lunbcs, *) nblk
   if (nblk /= nblocks) then
      write (luncrt, '(/, a, 2i5)') &
         ' Mismatched block count in grid & BCs file:', nblocks, nblk
      go to 99
   end if

!  Scan it for the BCs:

   allocate (ibc(6,nblocks))

   do
      read (lunbcs, '(a)') filename(1:32)
      if (index (filename(1:32), 'Boundary condition type') > 0) exit
   end do

   nskip = 1
   do ib = 1, nblocks
      do i = 1, nskip
         read (lunbcs, *)
      end do
      read (lunbcs, *) ibc(:,ib)
      nskip = 22
   end do

   close (lunbcs)

!  Deal with opening the output volume grid:
!  -----------------------------------------

   read (lunctl, *)
   read (lunctl, *) filename
   read (lunctl, *) formatted_out

   i1 = 1; if (formatted_out) i1 = 3

   open (lungrd, file=filename, status='unknown', form=format(i1:11),          &
         iostat=ios)
   if (ios /= 0) then
      write (luncrt, '(/, 2a)') ' Cannot open the output volume grid file: ',  &
         trim (filename)
      go to 99
   end if

!  Read the outer boundary tailoring controls:
!  -------------------------------------------

   read (lunctl, *)
   read (lunctl, *)
   read (lunctl, *) method_edge     ! Shock leading edge detection method (or 0)
   read (lunctl, *) nsmoothf        ! # smoothing iters. applied to flow grad.
   read (lunctl, *) fs_scale        ! Scale factor applied to freestream value
   read (lunctl, *) ds_local_multiple ! Multiple of local ds for outer margin
   read (lunctl, *) smoothing_limiter ! For constraining smoothed results
   read (lunctl, *) nsmoothg        ! # smoothing iters. applied to outer bndry.
   read (lunctl, *) added_margin    ! Multiple of outer spacing , possibly < 0
   read (lunctl, *) nk_margin       ! Plane nk - nk_margin is treated as for 0

!  ... and the radial distribution controls:
!  -----------------------------------------

   read (lunctl, *)
   read (lunctl, *)
   read (lunctl, *) method_radial   ! Different methods need different controls
   read (lunctl, *) nradial         ! > 0 allows changing # radial points

   if (method_edge == 0) then
      if (method_radial == 2) flow_needed = false
   end if

   read (lunctl, *)
   read (lunctl, *)                 ! Radial method 1 controls: Re_cell-based
   read (lunctl, *) cell_Re         ! Nominal cell Reynolds # at wall
   read (lunctl, *) ds1_min         ! Minimum wall spacing before smoothing
   read (lunctl, *) ds1_max         ! Maximum wall spacing before smoothing
   read (lunctl, *) ngeometric      ! Used by BLGRID; <= 2 => pure Vinokur;
                                    ! -nbl < ngeometric < -2 => fix 1:nbl
   read (lunctl, *) growth_rate     ! Geometric growth rate for BLGRID
   read (lunctl, *) ds2_fraction    ! Applied to 1-sided ds2 to give 2-sided ds2

   read (lunctl, *)
   read (lunctl, *)                 ! Radial method 2: specified|existing ds1
   read (lunctl, *) ds1             ! ds1 to impose everywhere; 0. => existing
   read (lunctl, *) icontrols(1)    ! Used by BLGRID; <= 2 => pure Vinokur
   read (lunctl, *) rcontrols(1)    ! Geometric growth rate for BLGRID
   read (lunctl, *) rcontrols(2)    ! Applied to 1-sided ds2 to give 2-sided ds2

   if (method_radial == 2) then
      ngeometric   = icontrols(1)
      growth_rate  = rcontrols(1)
      ds2_fraction = rcontrols(2)
   end if

   read (lunctl, *)
   read (lunctl, *)                 ! Radial method 3 controls (|dT/ds|-based)
   read (lunctl, *) outer_fraction,&! Try 0.5
                    blend_fraction,&! Try 0.1
                    power           ! Try 0.5

   read (lunctl, *)
   read (lunctl, *)                 ! Plottable profile controls
   read (lunctl, *) ib_profiles     ! > 0 means print profiles from this block
   read (lunctl, *) is_profiles, &  ! (i,j) range of profiles printed
                    ie_profiles, &
                    js_profiles, &
                    je_profiles

!  Optional final input, with multiple uses:

   read (lunctl, *, iostat=ios) fs_value  ! Possibly for overriding block 1 use

   constant_margin = zero    ! Normally
   if (ios /= 0) then
      fs_value = zero
   else                      ! fs_value is used differently if method_edge = -1
      if (method_edge == 0) constant_margin = fs_value ! or if method_edge =  0
   end if

   close (lunctl)

!  Look for the optional ancillary control file:
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   allocate (suppress(nblocks), iblock(nblocks))

   nblank = 0
   suppress(:) = false

   open (lunctl, file='outbound.inp.2', status='old', iostat=ios)

   if (ios == 0) then

!     Option to blank some blocks:

      nblank = nblocks  ! Upper limit for RDLIST

      write (luncrt, '(a)')

      call rdlist (luncrt, 'Reading outbound.inp.2 for blocks to suppress. ',  &
                   lunctl, nblank, iblock)

      do i = 1, nblank
         ib = iblock(i)
         if (ib < 1 .or. ib > nblocks) then
            write (luncrt, '(/, a, i9, a)') &
               ' Bad block number to suppress:', ib, '; aborting.'
            go to 99
         else
            suppress(ib) = true
         end if
      end do

      close (lunctl)

   end if

   if (nblank > 0) then
      write (luncrt, '(/, a, (20i4))') ' Blocks to blank:', iblock(1:nblank)
   else
      write (luncrt, '(/, a)') ' No blocks are being blanked.'
   end if

   deallocate (iblock)

!  Some sanity checks:

   if (iq > nf) then
      write (luncrt, '(/, a, 2i4)') &
         ' Tailoring variable > # variables found:', iq, nf
      go to 99
   end if

   if (method_radial == 1) then
      if (nf < 3) then
         write (luncrt, '(/, a, 2i4)') &
            ' Cell_Re > 0. requires rho, a, mu in f(1:3,*,*,*); iq, nf:', iq, nf
         go to 99
      end if
   end if

!  Combine the integer and real controls for passing to the tailoring routine:

   icontrols(1)  = method_edge
   icontrols(2)  = method_radial
   icontrols(3)  = nsmoothf
   icontrols(4)  = nsmoothg
   icontrols(5)  = nradial
   icontrols(6)  = ngeometric
   icontrols(7)  = ib_profiles
   icontrols(8)  = is_profiles
   icontrols(9)  = ie_profiles
   icontrols(10) = js_profiles
   icontrols(11) = je_profiles
   icontrols(12) = nk_margin

   rcontrols(1)  = fs_scale
   rcontrols(2)  = ds_local_multiple
   rcontrols(3)  = smoothing_limiter
   rcontrols(4)  = cell_Re
   rcontrols(5)  = ds1_min
   rcontrols(6)  = ds1_max
   rcontrols(7)  = ds1
   rcontrols(8)  = growth_rate
   rcontrols(9)  = ds2_fraction
   rcontrols(10) = added_margin


!  Write the output grid header:
!  -----------------------------

   if (ndim == 2 .and. nradial == grid(1)%nj .or. &
       ndim == 3 .and. nradial == grid(1)%nk) then
       write (luncrt, '(/, 2a, i4)') &
          ' Ignoring request to change radial point count. ', &
          ' Block 1 already matches the request:', nradial
       nradial = 0
   end if

   if (nradial == 0) then ! The block dimensions will not be changed

      if (ndim == 2) then
         call xy_header_io  (2, lungrd, formatted_out, nblocks, grid, ios)
      else
         call xyz_header_io (2, lungrd, formatted_out, nblocks, grid, ios)
      end if

   else ! Option to change the point counts in the radial direction (only)

      allocate (new_grid(nblocks))

      if (ndim == 2) then

         nj = grid(1)%nj
         do ib = 1, nblocks
            new_grid(ib)%ni = grid(ib)%ni
            new_grid(ib)%nj = nradial
            if (grid(ib)%nj /= nj) nj = 0
         end do

         if (nj /= 0) then

            call xy_header_io (2, lungrd, formatted_out, nblocks, new_grid, ios)

         else ! We can't handle nonstandard block indexing and still process
              ! only 1 block at a time in memory because of this header output.

            write (luncrt, '(/, 2a, i5)') &
               ' Radial point counts nj are not the same for all blocks.', &
               ' Ignoring nradial input:', nradial
               nradial = 0
            deallocate (new_grid)

         end if

      else  ! 3-D

         nk = grid(1)%nk
         do ib = 1, nblocks
            new_grid(ib)%ni = grid(ib)%ni
            new_grid(ib)%nj = grid(ib)%nj
            new_grid(ib)%nk = nradial
            if (grid(ib)%nk /= nk) nk = 0
         end do

         if (nk /= 0) then

            call xyz_header_io (2, lungrd, formatted_out, nblocks, new_grid,ios)

         else

            write (luncrt, '(/, 2a, i5)') &
               ' Radial point counts nk are not the same for all blocks.', &
               ' Ignoring nradial input:', nradial
               nradial = 0
            deallocate (new_grid)

         end if

      end if

   end if

   if (ios /= 0) then
      write (luncrt, '(/, a)') ' Trouble writing the output grid header.'
      go to 99
   end if

!  Option to scale radial line offsets from each surface point:

   if (method_edge == -2) then

!     Avoid switching to and from DPLR indexing, etc.:

      call scale_grid (ndim, lung, formatted, nblocks, grid,   &
                       icontrols, rcontrols, fs_value, lungrd, &
                       formatted_out, ios)
      go to 90

   end if

!  Option to compress a hyperbolic volume grid for a better starting guess:

   if (method_edge == -1) then

      if (ndim == 2) then  ! More trouble than it's worth
         write (luncrt, '(/, a)') 'Grid compression is a 3-D option only.'
         go to 99
      end if

      call scan_hyperbolic_grid (lung, formatted, nblocks, grid, &
                                 xyz_at_xmin, rmaxh, right_half, left_half, ios)
      if (ios /= 0) then
         write (luncrt, '(/, a, i4)') ' Trouble scanning grid.  Status:', ios
         go to 99
      end if

!     The grid file has been rewound.  Reread its header:

      call xyz_header_io (1, lung, formatted, nblocks, grid, ios)
      if (ios /= 0) go to 99

!     Avoid switching to and from DPLR indexing, etc.:

      call compress_grid (lung, formatted, nblocks, grid,            &
                          xyz_at_xmin, rmaxh, right_half, left_half, &
                          icontrols, rcontrols, fs_value, lungrd,    &
                          formatted_out, ios)
      go to 90

   end if

   write (luncrt, '(a)')
   nlines_to_average = 0    ! Some blocks may not have an outer boundary
   avg_of_fractions  = zero ! For average fractional underestimate of new bndry.

!  Process the grid one block at a time:
!  -------------------------------------

   do ib = 1, nblocks

      ni = grid(ib)%ni
      nj = grid(ib)%nj
      nk = grid(ib)%nk

!     Read the current grid block:

      if (ndim == 2) then
         call xy_allocate  (grid(ib), ios)
      else
         call xyz_allocate (grid(ib), ios)
      end if

      if (ios /= 0) then
         write (luncrt, '(/, a, i4, a, i4)') &
            ' Trouble allocating input grid block #', ib, '.  iostat:', ios
         write (luncrt, '(a, 3i5)') ' Dimensions: ', ni, nj, nk
         go to 99
      end if

      if (ndim == 2) then
         npts = ni * nj
         call xy_block_io (1, lung, formatted, npts, grid(ib)%x, &
                           grid(ib)%y, ios)
         allocate (grid(ib)%z(ni,nj,1))
         grid(ib)%z(:,:,:) = zero
      else
         npts = ni * nj * nk
         call xyz_block_io (1, lung, formatted, npts, grid(ib)%x, &
                            grid(ib)%y, grid(ib)%z, ios)
      end if

      if (ios /= 0) then
         write (luncrt, '(/, a, i4, a, i4)') &
            ' Trouble reading grid block #',  ib, '.  iostat:', ios
         write (luncrt, '(a, 3i5)') ' Dimensions: ', ni, nj, nk
         go to 99
      end if

!     Allocate the input function file block:

      n = max (nf, 1)
      if (ndim == 2) then
         call q_allocate_2d (grid(ib), n, ios)
      else
         call q_allocate    (grid(ib), n, ios)
      end if

      if (ios /= 0) then
         write (luncrt, '(/, a, i4, a, i4)') &
            ' Trouble allocating fn. space for block ', ib, '  iostat:', ios
         write (luncrt, '(a, 4i5)') ' Dimensions: ', grid(ib)%mi, &
            grid(ib)%mj, grid(ib)%mk, n
         go to 99
      end if

      if (flow_needed) then

         if (ndim == 2) then
            call q_block_io_2d (1, lunf, formatted, nf, ni, nj, grid(ib)%q, ios)
         else
            call q_block_io (1, lunf, formatted, nf, ni, nj, nk, grid(ib)%q,ios)
         end if

         if (ios /= 0) then
            write (luncrt, '(/, a, i4, a, i4)') &
               ' Trouble reading function block #',  ib, '.  iostat:', ios
            write (luncrt, '(a, 3i5)') ' Dimensions: ', ni, nj, nk
            go to 99
         end if

      else ! Avoid undefined flow variables

         grid(ib)%q(:,:,:,:) = zero

      end if

!     Perform the indicated tailoring for this block, returning it in-place:
!     ----------------------------------------------------------------------

!!!   call tailor_block (ib, iq, nf, ni, nj, nk, ibc(:,ib),                    &
!!!                      icontrols, rcontrols,                                 &
!!!                      grid(ib)%x, grid(ib)%y, grid(ib)%z, grid(ib)%q,       &
!!!                      fpeak_min, kpeak_min, n_too_close, worst_fraction,    &
!!!                      sum_of_fractions, nlines_tailored)

!     The above module has been decomposed so that everything except the
!     smoothing of surface data can be performed in parallel on split blocks in
!     the DPLR implementation.  Here, all steps are performed on whole blocks.


!     Determine relevant faces and whether (sub)block indices need permuting:
!     -----------------------------------------------------------------------

      call identify_boundaries (ndim, ibc(:,ib), ifree, iwall, k1_face, iadapt)

      if (iadapt == 0 .or. &  ! No outer boundary
          suppress(ib)) then  ! Transcribe only - block is not touched

         if (ndim == 2) then
            call xy_block_io  (2, lungrd, formatted_out, npts, &
                               grid(ib)%x, grid(ib)%y, ios)
         else
            call xyz_block_io (2, lungrd, formatted_out, npts, &
                               grid(ib)%x, grid(ib)%y, grid(ib)%z, ios)
         end if

         if (ios /= 0) then
            write (luncrt, '(/, a, i4, a, i4)') &
               ' Trouble writing output grid block #', ib, '  I/O status:', ios
            write (luncrt, '(a, 3i5)') ' Dimensions: ', ni, nj, nk
            go to 99
         end if

         deallocate (grid(ib)%x, grid(ib)%y, grid(ib)%z, grid(ib)%q)

         write (luncrt, '(a, i4)') ' Transcribed block', ib

         cycle

      end if

!     Adapt to DPLR's (n,k,j,i) convention so that the key utilities are the
!     same in both codes:
!     ----------------------------------------------------------------------

      call swap_ijk_conventions (1, nf, ni, nj, nk,                            &
                                 grid(ib)%x, grid(ib)%y, grid(ib)%z, grid(ib)%q)

!     Ensure that the (sub)block indices are in the preferred order, namely
!     with k = 1 at the face opposite the outer boundary face.
!     The output mi, mj, mk become the active dimensions with k = mk at the
!     outer face.  The coordinates and function values are permuted in place.
!     -----------------------------------------------------------------------

      if (ndim == 2) then  ! Make it look like 3-D
         nk = nj;  mk = nk
         nj = 1;   mj = 1
                   mi = ni
         k1_face = 5  ! Assume no permutations
      else
         call permute_block (1, nf, ni, nj, nk, k1_face, &
                             grid(ib)%x, grid(ib)%y, grid(ib)%z, grid(ib)%q, &
                             mi, mj, mk)
      end if

      kb = mk - nk_margin   ! k plane to be treated as the outer boundary
                            ! for the standard shock alignment procedure

      if (flow_needed) then
         if (fs_value == zero) then
             if (ndim == 2) then
                fs_value = grid(ib)%q(iq,1,mk,1)  ! Once, whole grid
             else
                fs_value = grid(ib)%q(iq,mk,1,1)
             end if
             write (luncrt, '(/, a, es15.7)') &
                ' Using free stream flow value from block 1, point (1,1,nk):', &
                fs_value
         end if
      end if


!     Simple clean-up of existing boundary, or proper tailoring?
!     ----------------------------------------------------------

      if (method_radial /= 3) then       ! Not if we're reclustering outer grid
         allocate (radial_change(mj,mi)) ! This is what we smooth (whole block):
                                         ! delta arcs, unless method_edge = 0
         allocate (original_arcs(mj,mi)) ! Now used for smoothing delta arcs too
         allocate (stored_deltas(mj,mi)) ! To help constrain the smoothing
      end if

      if (method_edge == 0) then

         if (method_radial /= 3) then    ! Smooth the arc lengths directly
                                         ! (no tailoring):
            do i = 1, mi
               do j = 1, mj

!                 The following internal procedure gets around subscript-
!                 checking traps affecting the 2-D case:

                  call get_arc_length (mk, mj, mi, &
                                       grid(ib)%x, grid(ib)%y, grid(ib)%z, &
                                       j, i, ds, arc_length)

                  original_arcs(j,i) = arc_length
                  stored_deltas(j,i) = ds                         ! Outermost ds
                  radial_change(j,i) = arc_length + constant_margin &
                                                  + added_margin*ds  ! Turned to
               end do                                             ! deltas after
            end do                                                ! smoothing

         else  ! Radial method 3 (recluster outer radial pts. based on |dT/ds|)

            call recluster (mk, mj, mi, nf, &
                            grid(ib)%x, grid(ib)%y, grid(ib)%z, grid(ib)%q, &
                            outer_fraction, blend_fraction, power)
         end if

      else  ! Normal tailoring to the shock envelope:

!        Determine the preliminary (unsmoothed) new outer boundary location for
!        this (sub)block. Actually, we calculate unsmoothed CHANGES in location.
!        -----------------------------------------------------------------------


         call preliminary_boundary (ib, iq, nf, mi, mj, mk, icontrols,         &
                                    rcontrols, grid(ib)%x, grid(ib)%y,         &
                                    grid(ib)%z, grid(ib)%q, fs_value,          &
                                    radial_change, stored_deltas,              &
                                    fpeak_min, kpeak_min,                      &
                                    n_too_close, worst_fraction,               &
                                    sum_of_fractions)

         avg_of_fractions   = avg_of_fractions + sum_of_fractions
         nlines_tailored    = mi * mj
         nlines_to_average  = nlines_to_average + nlines_tailored
         original_arcs(:,:) = radial_change(:,:)              ! Before smoothing

      end if


!     Smooth the outer boundary (method_edge = 0) or boundary changes:
!     ----------------------------------------------------------------

      if (nsmoothg /= 0 .and. method_radial /= 3) then ! Smooth in place

!        We smooth with respect to arc lengths along outer boundary grid lines:

         allocate (xs(mj,mi), ys(mj,mi), zs(mj,mi), u(mj,mi), v(mj,mi))

         call get_k_boundary (mk, mj, mi, mk, &
                              grid(ib)%x, grid(ib)%y, grid(ib)%z, xs, ys, zs)

         u(1,1) = -999. ! Suppresses normalization of arc lengths

         call param2d (mj, mi, 1, mj, 1, mi, xs, ys, zs, u, v)

         call smooth2d (nsmoothg, mj, mi, u, v, radial_change)

         deallocate (xs, ys, zs, u, v)

      end if


      if (method_edge == 0) then  ! No tailoring

         if (method_radial /= 3) then  ! Not if reclustering outer portion

!           We're just smoothing full arcs.
!           Don't allow the smoothing to overdo things.
!           Recall that "radial_change" here is really the new full arc length.
!           Use a large smoothing_limiter if added_margin is nonzero.
!           Likewise if constant_margin is nonzero.

            do i = 1, mi
               do j = 1, mj
                  arc_limit          = stored_deltas(j,i)  * smoothing_limiter
                  arc_added          = stored_deltas(j,i)  * added_margin
                  arc_unsmoothed     = original_arcs(j,i)  + arc_added
                  radial_change(j,i) = min (arc_unsmoothed + arc_limit, &
                                       max (arc_unsmoothed - arc_limit, &
                                            radial_change(j,i)))
               end do
            end do

!           We avoid touching the redistribution routine for this retrofit by
!           deriving delta arc lengths from the smoothed arc lengths:

            radial_change(:,:) = original_arcs(:,:) - radial_change(:,:)

         end if

      else ! Arc length changes have been smoothed: constrain that smoothing

!        Recall that "original_arcs" here is really the unsmoothed arc changes.
!        Positive radial_change means shrink the outer boundary.

         do i = 1, mi
            do j = 1, mj
               arc_limit          = stored_deltas(j,i) * smoothing_limiter
               radial_change(j,i) = min (original_arcs(j,i) + arc_limit,       &
                                    max (original_arcs(j,i) - arc_limit,       &
                                         radial_change(j,i)))
            end do
         end do

      end if

!     Set up for radial redistribution, unless ...
!     --------------------------------------------

      if (method_radial /= 3) then
         deallocate (original_arcs, stored_deltas)
         allocate (d1_cell_Re(mj,mi))  ! Always, because it's an argument below
      end if

!     Do we impose a cell Reynolds number?
!     The set-up of unsmoothed wall spacings can be done on sub-blocks,
!     but the smoothing must be done for the whole block.

      if (method_radial == 1) then  ! Cell Reynolds number option

         if (ds1_max > 0) then ! Same absolute limits on ds1 everywhere

            do i = 1, mi
               do j = 1, mj
                  call wall_spacing (nf, mk, mj, mi, grid(ib)%q, j, i, &
                                     ds1_min, ds1_max, cell_Re, d1_cell_Re(j,i))
               end do
            end do

         else ! ds1_max < 0 implies limits are arc length percentages

            do i = 1, mi
               do j = 1, mj

                  call get_arc_length (mk, mj, mi, &
                                       grid(ib)%x, grid(ib)%y, grid(ib)%z, &
                                       j, i, ds, arc_length)

                  arc_length = arc_length - radial_change(j,i) ! Final length

                  ds1_bl = (-0.01 * ds1_min) * arc_length
                  ds1_bu = (-0.01 * ds1_max) * arc_length

                  call wall_spacing (nf, mk, mj, mi, grid(ib)%q, j, i, &
                                     ds1_bl, ds1_bu, cell_Re, d1_cell_Re(j,i))
               end do
            end do

         end if

!        Smooth the cell-Reynolds-number-based increments at the wall:

         allocate (xs(mj,mi), ys(mj,mi), zs(mj,mi), u(mj,mi), v(mj,mi))

         call get_k_boundary (mk, mj, mi, 1, &
                              grid(ib)%x, grid(ib)%y, grid(ib)%z, xs, ys, zs)

         u(1,1) = -999. ! Suppresses arc length normalization

         call param2d (mj, mi, 1, mj, 1, mi, xs, ys, zs, u, v)

         call smooth2d (nsmooth_d1, mj, mi, u, v, d1_cell_Re)

         deallocate (xs, ys, zs, u, v)

      end if


!     Redistribute the points for all radial lines of the current (sub)block:
!     -----------------------------------------------------------------------

      if (method_radial /= 3) then  ! Not if they're already reclustered

         call redistribute (method_radial, mi, mj, mk, kb, iwall, ngeometric,  &
                            growth_rate, radial_change, d1_cell_Re, ds1,       &
                            ds2_fraction, grid(ib)%x, grid(ib)%y, grid(ib)%z)

         deallocate (radial_change, d1_cell_Re)
      end if


!     Reverse any indexing permutation of this (sub)block:
!     ----------------------------------------------------

      if (ndim == 3) then
         call permute_block (2, nf, mi, mj, mk, k1_face, &
                             grid(ib)%x, grid(ib)%y, grid(ib)%z, grid(ib)%q, &
                             ni, nj, nk)
      end if


!     Switch back to conventional (i,j,k) order.  Dimensions match input x/y/z:
!     -------------------------------------------------------------------------

      if (ndim == 2) then  ! Undo making it look like 3-D
         nj = nk
         nk = 1
      end if

      call swap_ijk_conventions (2, nf, ni, nj, nk,                            &
                                 grid(ib)%x, grid(ib)%y, grid(ib)%z, grid(ib)%q)


!     Option to adjust the point counts (not applicable to DPLR):
!     -----------------------------------------------------------

      if (nradial > 0) then

         call xyz_allocate (new_grid(ib), ios)
         if (ios /= 0) go to 99

         if (ndim == 2) then  ! Make it look like 3-D
            nk = nj
            nj = 1
         end if

         call adjust_point_counts (ib, ni, nj, nk, nradial,                    &
                                   grid(ib)%x,     grid(ib)%y,     grid(ib)%z, &
                               new_grid(ib)%x, new_grid(ib)%y, new_grid(ib)%z, &
                                   icontrols, rcontrols)
         nk   = nradial
         npts = ni * nj * nk

         if (ndim == 2) then
            call xy_block_io  (2, lungrd, formatted_out, npts, &
                               new_grid(ib)%x, new_grid(ib)%y, ios)
         else
            call xyz_block_io (2, lungrd, formatted_out, npts, &
                            new_grid(ib)%x, new_grid(ib)%y, new_grid(ib)%z, ios)
         end if

         deallocate (new_grid(ib)%x, new_grid(ib)%y, new_grid(ib)%z)

      else

!        Just save the adjusted grid block:
!        ----------------------------------

         if (ndim == 2) then
            call xy_block_io  (2, lungrd, formatted_out, npts, &
                               grid(ib)%x, grid(ib)%y, ios)
         else
            call xyz_block_io (2, lungrd, formatted_out, npts, &
                               grid(ib)%x, grid(ib)%y, grid(ib)%z, ios)
         end if

      end if

      if (ios /= 0) then
         write (luncrt, '(/, a, i4, a, i4)') &
            ' Trouble writing output grid block #', ib, '  I/O status:', ios
         write (luncrt, '(a, 3i5)') ' Dimensions: ', ni, nj, nk
         go to 99
      end if

      deallocate (grid(ib)%x, grid(ib)%y, grid(ib)%z, grid(ib)%q)

      select case (method_edge)

         case (0) ! Just cleaning up the outer boundary

            write (luncrt, '(a, i4)') ' Completed block', ib

         case (1, 4)

            write (luncrt, '(a, i4, a, i6, a, f10.6)')                         &
               ' Block', ib, ':  # boundaries too close:',                     &
               n_too_close, ';  worst underestimate:', worst_fraction

         case default

            write (luncrt, '(a, i4, a, i6, a, f10.6, a, f10.6, a, i4)')        &
               ' Block', ib, ':  # boundaries too close:',                     &
               n_too_close, ';  worst underestimate:', worst_fraction,         &
               ';  min. shock peak:', fpeak_min, ' at k =', kpeak_min

      end select

   end do ! Next grid block


   if (method_edge > 0) then
      avg_of_fractions = avg_of_fractions / real (nlines_to_average)
      write (luncrt, '(/, a, f11.6)') &
         ' Average outer boundary underestimate as an arc length fraction:',   &
         avg_of_fractions
   end if

90 close (lung)
   close (lungrd)

!  Done.

99 continue

!  Internal procedure for program outbound:

   contains

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine get_arc_length (mk, mj, mi, x, y, z, jin, iin, ds, arc_length)
!
!     This routine overcomes subscript-checking trouble for the 2-D case.
!     It computes the total arc length and outermost spacing in the k direction
!     for the (j,i) surface point of the indicated grid block, using DPLR's
!     reversed i,j,k convention.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in)                       :: mk, mj, mi
      real,    intent (in), dimension (mk,mj,mi) :: x, y, z
      integer, intent (in)                       :: jin, iin
      real,    intent (out)                      :: ds          ! Last delta
      real,    intent (out)                      :: arc_length  ! Total arc

!     Local variables:

      integer :: i, j, k
      real    :: delta, sum

!     Execution:

      i = iin;  j = jin;  sum = zero

      do k = 2, mk
         delta = sqrt ((x(k,j,i) - x(k-1,j,i))**2 +  &
                       (y(k,j,i) - y(k-1,j,i))**2 +  &
                       (z(k,j,i) - z(k-1,j,i))**2)
         sum = sum + delta
      end do

      ds = delta
      arc_length = sum

      end subroutine get_arc_length

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine get_k_boundary (mk, mj, mi, kplane, x, y, z, xs, ys, zs)
!
!     This routine overcomes subscript-checking trouble for the 2-D case.
!     It transfers the indicated innermost or outermost k plane to 2-D arrays.
!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in)                        :: mk, mj, mi
      integer, intent (in)                        :: kplane   ! k plane to copy
      real,    intent (in),  dimension (mk,mj,mi) :: x, y, z
      real,    intent (out), dimension (mj,mi)    :: xs, ys, zs

!     Local variables:

      integer :: i, j, k

!     Execution:

      k = kplane

      do i = 1, mi
         do j = 1, mj
            xs(j,i) = x(k,j,i)
            ys(j,i) = y(k,j,i)
            zs(j,i) = z(k,j,i)
         end do
      end do

      end subroutine get_k_boundary

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine wall_spacing (nf, mk, mj, mi, q, j, i, d1min, d1max, &
                               cell_Re, d1_cell_Re)
!
!     This routine overcomes subscript-checking trouble for the 2-D case.
!     It computes cell-Reynold-number-based wall spacing and imposes upper and
!     lower bounds for the (j,i) surface point of the indicated grid block,
!     using DPLR's reversed i,j,k convention.

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      integer, intent (in)  :: nf
      integer, intent (in)  :: mk, mj, mi
      real,    intent (in)  :: q(nf,mk,mj,mi)
      integer, intent (in)  :: j, i
      real,    intent (in)  :: d1min, d1max  ! Bounds on the spacing
      real,    intent (in)  :: cell_Re       ! Requested cell Reynolds number
      real,    intent (out) :: d1_cell_Re    ! Corresp. constrained wall spacing

!     Local variables:

      real :: d1

!     Execution:

      d1 = cell_Re * grid(ib)%q(3,1,j,i) / &
                    (grid(ib)%q(1,1,j,i) * grid(ib)%q(2,1,j,i))
      d1_cell_Re = min (max (d1, d1min), d1max)

      end subroutine wall_spacing

   end program outbound

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine identify_boundaries (ndim, ibc, ifree, iwall, k1_face, iadapt)
!
!  This is the first step in a loop over (sub)blocks for tailoring the grid
!  outer boundary.
!
!  Determine the free-stream block face (if any), the solid wall face (if any),
!  and the face that (if it isn't already face 5 / k = 1) controls permuting of
!  the (sub)block indices into standard order so that k is the radial direction.
!
!  08/17/05  DAS  Mike asked for a more general "iadapt" return argument.
!  12/12/13   "   Handled the 2-D case.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)  :: ndim     ! 2 or 3
   integer, intent (in)  :: ibc(6)   ! DPLR boundary conditions:
                                     !    1 => free-stream;
                                     !    2:4 => supersonic exit;
                                     !    8 => wall, inviscid;
                                     !    9, 10, 25-30 => wall, viscous;
                                     !    20 => zonal boundary

   integer, intent (out) :: ifree    ! 0 => no outer boundary face found

   integer, intent (out) :: iwall    ! May affect radial redistribution;
                                     ! 0 => k = 1 is NOT a solid wall

   integer, intent (out) :: k1_face  ! Face number to change to k = 1 (face 5)
                                     ! if it isn't already
   integer, intent (out) :: iadapt   ! 0 => leave the block alone

!  Local variables:

   integer :: i

!  Execution:

   iadapt = 0;  ifree = 0;  iwall = 0

   do i = 1, 2*ndim
      select case (ibc(i))
         case (1)
            ifree = i
         case (8:10, 25:49, 101:999)
            iwall = i
         case default
            ! Not a wall or free-stream boundary
      end select
   end do

!  A supersonic exit face may need treating as for a free-stream face:

   if (ifree == 0) then

      do i = 1, 2*ndim
         select case (ibc(i))
            case (2:4)
               ifree = i
            case default
               ! Not an exit face
         end select
      end do

   end if

   if (ifree == 0) go to 99 ! No outer boundary face - leave the block unchanged


   iadapt = 1

!  Subroutine permute_block expects as input the face number to make k = 1,
!  so we have to find it for the case where there is no wall face.
!  Assume we won't permute a 2-D block.

   if (iwall == 0) then
      if (mod (ifree, 2) == 0) then
         k1_face = ifree - 1
      else
         k1_face = ifree + 1
      end if
   else
      k1_face = iwall
   end if

99 return

   end subroutine identify_boundaries

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine preliminary_boundary (ib, iq, nf, mi, mj, mk,                    &
                                    icontrols, rcontrols, x, y, z, f,          &
                                    fs_value, radial_change, stored_deltas,    &
                                    fpeak_min, kpeak_min, n_too_close,         &
                                    worst_fraction, sum_of_fractions)
!
!  This is the second main step in a loop over grid blocks for tailoring the
!  outer boundary to the shock envelope.
!  For one block or sub-block, determine the preliminary (unsmoothed) location
!  of the adjusted outer boundary, or rather the unsmoothed change in location.
!  The present method is ill-suited to allowing the outer boundary to grow
!  rather than shrink, although that may be managed at the higher level.
!
!  June, 2005  DS/SY  Initial OUTBOUND implementation (1-function tailoring).
!  07/14/05    DAS    Handle more than one flow variable for the cell-Re option.
!  08/09/05     "     Decompose original TAILOR_BLOCK routine for application to
!                     sub-blocks along with smoothing at the whole block level.
!  11/23/05     "     Mike asked for an added margin option, same everywhere.
!  12/17/06     "     Store the local edge spacings to help constrain smoothing.
!
!  Authors:  David Saunders, ELORET and Seokkwan Yoon, NASA Ames Research Center
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)  :: ib               ! Current block number
   integer, intent (in)  :: iq               ! Flow var. to use for tailoring
   integer, intent (in)  :: nf               ! # flow variables in f(1:nf,*,*,*)
   integer, intent (in)  :: mi, mj, mk       ! Block dimensions after permuting
   integer, intent (in)  :: icontrols(*)     ! Integer controls; see code below
   real,    intent (in)  :: rcontrols(*)     ! Real     "    "    "    "    "
   real, intent (in), dimension (mk,mj,mi) :: x, y, z ! Current (sub)block
   real,    intent (inout)::f(nf,mk,mj,mi)   ! Accompanying flow variables
   real,    intent (in)  :: fs_value         ! Ref. free-stream value in f(iq,*)
                                             ! normally > 0; if negative, use
                                             ! local peak along each k line
   real,    intent (out) :: radial_change(mj,mi) ! Define unsmoothed new locns.
   real,    intent (out) :: stored_deltas(mj,mi) ! Local spacings near the edge
   real,    intent (out) :: fpeak_min        ! Lowest flow pk found within shock
   integer, intent (out) :: kpeak_min        ! Corresponding radial index
   integer, intent (out) :: n_too_close      ! # lines for which the initial
                                             ! new boundary estimate was inside
                                             ! the fr.strm. flow (then adjusted)
   real,    intent (out) :: worst_fraction   ! Largest fraction of arc length by
                                             ! which an edge estimate is inside
                                             ! the freestream value
   real,    intent (out) :: sum_of_fractions ! For averaging the underestimates

!  Local constants:

   real,    parameter :: one = 1., zero = 0.
   logical, parameter :: false = .false.

!  Local variables:

   integer :: i, j, k, kb, kedge, kpeak, method_edge, nk_margin
   integer :: ib_profiles, ie_profiles, is_profiles, je_profiles, js_profiles

   real :: added_margin, ds_local, ds_local_multiple, fpeak, fraction,         &
           free_stream, fs_multiple, fs_scale, growth_rate, rncells,           &
           smoothing_limiter, s_kb, s_edge, s_shock, s_total_new, s_total_old

   logical   :: monitor_peaks

   character :: flow_variable * 32

   real, allocatable, dimension (:) :: sold, xold, yold, zold, gradient

!  Execution:

   method_edge = icontrols(1)  ! 1|4: Mach; 2: p,rho,T; 3: flow gradient
   ib_profiles = icontrols(7)  ! > 0 means output profile(s) from this block
   nk_margin   = icontrols(12) ! Option to leave outer boundary untouched
   kb = mk - nk_margin         ! k plane to treat as the outer boundary for the
                               ! normal shock alignment procedure
   if (ib_profiles > 0) then
      is_profiles = icontrols(8)
      ie_profiles = icontrols(9)
      js_profiles = icontrols(10)
      je_profiles = icontrols(11)
   end if

   fs_scale          = rcontrols(1)  ! Scale factor applied to free-stream val.
   ds_local_multiple = rcontrols(2)  ! Multiple of local ds for outer margin
!! smoothing_limiter = rcontrols(3)  ! Arc fraction to limit outer margin
   added_margin      = rcontrols(10) ! Multiple of the outermost 1-sided spacing

   free_stream = fs_value            ! May be reset to 1. below
   fs_multiple = fs_scale * free_stream

!  Option to output a series of QPLOTable profiles from a specified block:

   if (ib == ib_profiles) call save_profiles (1)  ! Plot header info.

!  Normalize pressure, density, or temperature:

   if (method_edge == 2) then  ! Normalize flow variable
      f = f * (one / free_stream)
      free_stream = one
   end if

!  Work-space for treating one radial line:

   allocate (xold(mk), yold(mk), zold(mk), sold(mk), gradient(mk))

!  Locate the preliminary new outer boundary for every radial line:

   fpeak_min        = 1.e+30
   kpeak_min        = 0
   n_too_close      = 0
   worst_fraction   = zero
   sum_of_fractions = zero
   monitor_peaks    = method_edge /= 1 .and. method_edge /= 4

   do i = 1, mi

      do j = 1, mj

         do k = 1, mk
            xold(k) = x(k,j,i)
            yold(k) = y(k,j,i)
            zold(k) = z(k,j,i)
         end do

         call chords3d (mk, xold, yold, zold, false, s_total_old, sold)  ! Arcs

         call shock_location (method_edge, fs_scale, free_stream,              &
                              iq, nf, mi, mj, mk, f, i, j, sold, gradient,     &
                              kedge, ds_local, kpeak, fpeak, s_shock)

!        s_shock = s_total_old means no edge was found - leave boundary alone.
!        This helps deal with blocks in the wake of a capsule.
!        How do we signal the need to expand the boundary?
!        Presumably, added_margin can help, but it's added to wake blocks too ..

         stored_deltas(j,i) = ds_local  ! For use after smoothing

         if (monitor_peaks) then
            if (fpeak < fpeak_min) then ! Track lowest peak in the shock
               fpeak_min = fpeak        ! (for p, rho, or T, but not Mach)
               kpeak_min = kpeak
            end if
         end if

!        Add the outer margin:

         s_kb   = sold(kb)

!!!      s_edge = s_shock + min (ds_local_multiple * ds_local, &
!!!                              smoothing_limiter * s_total_old)
         s_edge = s_shock + ds_local_multiple * ds_local

         s_edge = min (s_edge, s_kb) + added_margin * (s_kb - sold(kb-1))

         if (s_edge < sold(kedge)) then

            n_too_close = n_too_close + 1
            fraction    = (sold(kedge) - s_edge) / s_kb
            sum_of_fractions = sum_of_fractions + fraction
            worst_fraction = max (worst_fraction, fraction)

!!!!        if (fraction > 0.17) then
!!!!           write (6, '(a, f10.6, a, 3i4, a, i3)') &
!!!!              ' Underestimation fraction:', fraction, &
!!!!              '  ib, i, j:', ib, i, j, '  kedge:', kedge
!!!!        end if

!!!!        s_edge = sold(kedge)  ! Gives too bumpy results

         end if

         radial_change(j,i) = s_kb - s_edge      ! For smoothing

         if (ib == ib_profiles) call save_profiles (2)  ! If in (i,j) range

      end do ! Next j

   end do ! Next i

   deallocate (sold, xold, yold, zold, gradient)

!  Local procedure for subroutine preliminary boundary:

   contains

!     --------------------------------------------------------------------------
!
      subroutine save_profiles (icall)
!
!     Save one or more radial profiles in QPLOTable form.
!
!     --------------------------------------------------------------------------

      integer, parameter   :: l = 9 ! Logical unit for profiles
      integer, intent (in) :: icall ! 1 => plot header; 2 => one profile

      if (icall == 1) then ! Plot header

         if (free_stream < one) then
            flow_variable = 'Normalized density'
         else if (free_stream < 45.) then
            flow_variable = 'Mach number'
         else
            flow_variable = 'Normalized pressure'
         end if

         open (l, file='outbound.qplot', status='unknown')

         write (l, '(2a)') trim (flow_variable), ' profiles'
         write (l, '(a, i4, a, i6, i4, i6, i4)') &
            'Grid Block:', ib, '  (i,j) range:', &
            is_profiles, ie_profiles, js_profiles, je_profiles
         write (l, '(a)') trim (flow_variable)
         write (l, '(a)') 'Radial arc length'
         write (l, '(a)') ' $options grid=1, ident=T, xnumbers=''whole'', $end'

      else ! Save a profile?

         if (js_profiles <= j .and. j <= je_profiles) then
            if (is_profiles <= i .and. i <= ie_profiles) then
               if (i == is_profiles .and. j == js_profiles) then ! First profile
                  write (l, '(a)') ' $options line=''dash'','
                  if (free_stream > one) then
                     write (l, '(a, f8.4, a)') &
                        ' legend=''Target flow value =', fs_multiple, ''','
                  else
                     write (l, '(a, f5.2, a)') &
                        ' legend=''Freestream multiple =', fs_multiple, ''','
                  end if
                  write (l, '(a)') ' $end'
                  write (l, '(2es12.4)') fs_multiple, zero, &
                                         fs_multiple, s_total_old
               end if
               write (l, '(a)') ' $options line=''solid'', fit=''tight'','
               write (l, '(a)') ' $end'
               write (l, '(a, 2i4, a)') '!legend=''(i, j) =', i, j, ''','
               write (l, '(a, i4)') &
                  '! Flow variable     Arc length     # points:', mk
               write (l, '(2es15.7)') (f(iq,k,j,i), sold(k), k = mk, 1,-1)
               write (l, '(a)') ' $options line=''symbol'', symbol=16,'
               if (i == ie_profiles .and. j == je_profiles) then
                  write (l, '(a)') ' legend=''Unsmoothed new boundary'','
               end if
               write (l, '(a)') ' $end'
               write (l, '(2es15.7)') free_stream, s_edge
               write (l, '(a)') ' $options line=''symbol'', symbol=17,'
               if (i == ie_profiles .and. j == je_profiles) then
                  write (l, '(a)') ' legend=''Original outer boundary'','
               end if
               write (l, '(a)') ' xshift=1, $end'
               write (l, '(2es15.7)') free_stream, s_total_old
            end if
         end if

      end if

      end subroutine save_profiles

   end subroutine preliminary_boundary

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine recluster (mk, mj, mi, nf, x, y, z, q, &
                         outer_fraction, blend_fraction, power)

!  Redistribute the outer portions of all radial lines of one (sub)block in
!  DPLR standard form using the gradient w.r.t. arc length of temperature T,
!  assumed to be f(1,:,:,:).  Leave the inner boundary layer portion of the
!  points alone for each line, except for blending of spacing in the interval
!
!                        k  = km - kb : km + kb    where
!                        kb = blend_fraction * mk  (say 0.1 mk)
!                        nr = outer_fraction * mk  (say 0.5 mk)
!                        km = mk - nr + 1          (k = km:mk initially touched)
!
!  The (x,y,z) coordinates are updated in place.
!
!  12/18/2013  David Saunders  Adaptation of an option in NEQAIR_DATA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in) :: mk, mj, mi      ! (Sub)block dimensions
   integer, intent (in) :: nf              ! Number of functions in f(:,:,:,:)
   real,    intent (inout), dimension (mk,mj,mi) :: x, y, z
   real,    intent (in) :: q(nf,mk,mj,mi)  ! q(1,:,:,:) = T
   real,    intent (in) :: outer_fraction  ! E.g., .5 (of pts. redistributed)
   real,    intent (in) :: blend_fraction  ! E.g., .1 defining k = nk*(.5+|-.1)
                                           ! where a 2-sided Vinokur distribn.
                                           ! is imposed after reclustering the
                                           ! points k = .5*nk : nk.
   real,    intent (in) :: power           ! Moderates the |dT/ds|-based shape
                                           ! function; try 0.5 (or smaller if
                                           ! some lines don't converge)
!  Local constants:

   integer,       parameter :: ismooth = 1   ! Smooth |dT/ds|, not ensuing arcs
   integer,       parameter :: lunout  = -6  ! Diagnostics
   logical,       parameter :: normlz  = .false.  ! Don't normalize arcs here
   character (1), parameter :: method  = 'M' ! Monotonic local cubic interpltns.

!  Local variables:

   integer :: i, ier, j, k1, k2, kr, lunerr, nb, nblend, nredistrib
   real    :: d1, d2, exponent, total
   real, allocatable, dimension (:) :: f, s, snew, xnew, ynew, znew

!  Execution:

   lunerr     = abs (lunout)
   nredistrib = nint (outer_fraction * real (mk))
   nblend     = nint (blend_fraction * real (mk))
   exponent   = power
   kr         = mk - nredistrib + 1             ! Start of initial redistribn.
   k1         = kr - nblend - 1                 ! Start of blend
   k2         = kr + nblend + 1                 ! End of blend
   nb         = 2*nblend + 3                    ! # in blend region, w/ overlap

   allocate (f(mk), s(mk), snew(mk), xnew(mk), ynew(mk), znew(mk))

!  For each radial grid line ...

   do i = 1, mi
      do j = 1, mj

!        The inner part of the line is untouched, but keep the indexing simple:

         xnew(:) = x(:,j,i)
         ynew(:) = y(:,j,i)
         znew(:) = z(:,j,i)

         call chords3d (mk, xnew, ynew, znew, normlz, total, s)

         snew(:) = s(:)
         f(:)    = q(1,:,j,i)

!        Initial reclustering of outer portion:

         call graddis3d2 (nredistrib, x(kr:mk,j,i), y(kr:mk,j,i), &
                          z(kr:mk,j,i), s(kr:mk), f(kr:mk),       &
                          nredistrib, exponent, ismooth, lunout,  &
                          xnew(kr:mk), ynew(kr:mk), znew(kr:mk), snew(kr), ier)

!        Blend the arcs across the split at k = kr:

         d1 = snew(k1+1) - snew(k1)
         d2 = snew(k2) - snew(k2-1)

         call vinokur (k1, k2, d1, d2, snew, lunerr, ier)

         call lcsfit2 (mk, s, x(:,j,i), method, nb, snew(k1), xnew(k1))
         x(k1:mk,j,i) = xnew(k1:mk)

         call lcsfit2 (mk, s, y(:,j,i), method, nb, snew(k1), ynew(k1))
         y(k1:mk,j,i) = ynew(k1:mk)

         call lcsfit2 (mk, s, z(:,j,i), method, nb, snew(k1), znew(k1))
         z(k1:mk,j,i) = znew(k1:mk)

      end do
   end do

   deallocate (f, s, snew, xnew, ynew, znew)

   end subroutine recluster

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine redistribute (method_radial, mi, mj, mk, kb, iwall, ngeometric,  &
                            growth_rate, radial_change, d1_cell_Re, ds1,       &
                            ds2_fraction, x, y, z)
!
!  Redistribute the grid points for all radial lines of one (sub)block in DPLR
!  standard form, according to the indicated method.
!
!  Parametric local spline techniques are used to slide the existing points
!  along the original radial line.  The number of points, mk, is not changed.
!
!  June, 2005  D. Saunders  Initial implementation. Retaining the input relative
!              S. Yoon      spacing is expedient until good outer boundaries are
!                           in hand.
!  07/10/05    "  "         Geometric + Vinokur-type stretching (1- or 2-sided).
!  12/22/05    DAS          If a radial line is being lengthened, use linear
!                           extrapolation.  Parametric cubics can't be trusted.
!  05/05/06     "           ds1 < 0 means use % arc lengths for wall increments.
!  10/05/07     "           kb argument allows applying standard method to less
!                           than the full line and bridging the gap if kb < mk.
!  12/24/09     "           For method_radial = 2, ds1 = 999. means apply a
!                           uniform distribution, as may be desirable for
!                           elliptic smoothing purposes.
!  04/11/16     "           ds2_fraction = 0. means preserve outer spacings.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in) :: method_radial   ! 0 => same relative distribution;
                                           ! 1 => [geometric +] Vinokur;
                                           ! 2 => as for 1 (differs higher up);
                                           ! 3 => LAURA-type [inactive]
   integer, intent (in) :: mi, mj, mk      ! (Sub)block dimensions
   integer, intent (in) :: kb              ! kb <= mk is the k plane to treat as
                                           ! the outer boundary to be aligned
   integer, intent (in) :: iwall           ! 0 => k = 1 is NOT a solid wall
   integer, intent (in) :: ngeometric      ! > 1; 2 = pure Vinokur (method 1|2);
                                           ! -nbl < ngeometric < -2 => fix 1:nbl
   real,    intent (in) :: growth_rate     ! If ngeometric > 2 (method 1 or 2)
   real,    intent (in) :: radial_change(mj,mi) ! .. in outer bndary. arc length
   real,    intent (in) :: d1_cell_Re(mj,mi)    ! Initial increments (method 1)
   real,    intent (in) :: ds1             ! Initial increment (method 2);
                                           ! > 0.   => constant wall spacing;
                                           ! = 0.   => existing wall spacing;
                                           ! < 0.   => % arc length wall space;
                                           ! = 999. => uniform spacing
   real,    intent (in) :: ds2_fraction    ! Applied to 1-sided d1 to give
                                           ! the 2-sided d2; 1.0 => 1-sided
   real, intent (inout), dimension (mk,mj,mi) :: x, y, z  ! Current grid block

!  Local constants:

   integer,   parameter :: lunerr = 6
   real,      parameter :: zero   = 0.
   logical,   parameter :: false  = .false., true = .true.
   character, parameter :: method * 1 = 'M' ! Monotonic spline fits safeguard
                                            ! the symmetry plane
!  Local variables:

   integer :: i, ier, j, k, keval, kfix, kfixm1, npts
   real    :: d1, d2, derivs(3), ds01, r, s_total_new, s_total_old
   real    :: s_fix, s_kb, s_scale, s_shift
   logical :: new, uniform
   real, allocatable, dimension (:) :: sold, xold, yold, zold
   real, allocatable, dimension (:) :: snew, xnew, ynew, znew

!  Execution:

   uniform = ds1 == 999.  ! Kludge with elliptic volume grid smoothing in mind
   ds01    = 0.01 * ds1   ! Negative input means use % arc
   kfix    = -ngeometric  ! Negative input means fix points 1:kfix
   kfixm1  = kfix - 1
   npts    = mk - kfixm1 + 1

   allocate (xold(mk), yold(mk), zold(mk), sold(mk), &
             xnew(mk), ynew(mk), znew(mk), snew(mk))

   do i = 1, mi

      do j = 1, mj

         do k = 1, mk
            xold(k) = x(k,j,i)
            yold(k) = y(k,j,i)
            zold(k) = z(k,j,i)
         end do

         call chords3d (mk, xold, yold, zold, false, s_total_old, sold)

         s_kb = sold(kb)
         s_total_new = s_kb - radial_change(j,i)

         select case (method_radial) ! For new radial arc lengths

            case (0) ! Same relative distribution

               if (kfix <= 0) then

                  snew(1:kb) = (s_total_new / s_kb) * sold(1:kb)

               else  ! Don't touch points 1:kfix

                  snew(1:kfix) = sold(1:kfix)

                  s_fix   = sold(kfix)
                  s_scale = (s_total_new - s_fix) / (s_kb - s_fix)
                  s_shift = (s_kb - s_total_new) * s_fix / (s_kb - s_fix)

                  do k = kfix + 1, kb
                     snew(k) = s_scale * sold(k) + s_shift
                  end do

               end if

            case (1) ! Cell Reynolds #-based initial spacing; ignore ngeom < 0

               d1 = d1_cell_Re(j,i)

!              1-sided stretching first, then derive 2-sided d2 from it:

!              EXPDIS5 is equivalent to Vinokur's 1-sided tanh-based method.
!              MODE = 1 bunches towards the low end.
!              If increasing spacing (or close to it) is implied, it switches
!              to geometric spacing, with blending for d1 in [0.94 du, 0.97 du].

               call expdis5 (1, zero, s_total_new, d1, kb, snew, -lunerr)

               d2 = (s_total_new - snew(kb-1)) * ds2_fraction

!              Pure 2-sided Vinokur if ngeometric = 2 else geometric in b.layer:

               snew(1)  = sold(1)      ! Always zero?
               snew(kb) = s_total_new  ! BLGRID expects the end-points in place:

               call blgrid (kb, d1, d2, ngeometric, growth_rate, snew, lunerr, &
                            ier)

            case (2) ! Specified or existing d1 everywhere

               if (kfix > 2) then            ! Don't touch points 1:kfix
                  snew(1:kfix)    = sold(1:kfix)
                  d1 = sold(kfix) - sold(kfixm1)
               else                          ! As originally, before kfix option
                  kfix = 2;  kfixm1 = 1
                  npts = kb
                  if (ds1 == zero) then      ! Preserve existing initial spacing
                     d1 = sold(2) - sold(1)
                  else if (ds1 > zero) then  ! Constant initial spacing
                     d1 =  ds1
                  else                       !  ds1 < zero => % arc length
                     d1 = -ds01 * s_total_new
                  end if
               end if

               if (uniform) then

                  d1      = s_total_new / real (kb - 1)
                  snew(1) = zero

                  do k = 2, kb - 1
                     snew(k) = d1 * real (k - 1)
                  end do

                  snew(kb) = s_total_new

               else if (ds2_fraction <= zero) then  ! Preserve outer spacing

                  d2 = s_total_old - sold(mk-1)
                  snew(mk) = s_total_old
                  snew(1)  = sold(1)

                  call blgrid (mk, d1, d2, ngeometric, growth_rate, snew, &
                               lunerr, ier)
               else

!                 Preliminary 1-sided stretching as for case (1):

                  call expdis5 (1, sold(kfixm1), s_total_new, d1, npts,        &
                                snew(kfixm1), -lunerr)

                  d2 = (s_total_new - snew(kb-1)) * ds2_fraction

!                 2-sided [geometric + ] Vinokur:

                  snew(kfixm1) = sold(kfixm1)
                  snew(kb)     = s_total_new

                  call blgrid (npts, d1, d2, ngeometric, growth_rate,          &
                               snew(kfixm1), lunerr, ier)
               end if

            case default ! Trap a bad radial method at a higher level

         end select


!        Is the outer boundary really being held fixed?

         if (kb < mk) call bridge_to_fixed_boundary ()  ! Local procedure


!        Update the (x,y,z)s as functions of the new arc lengths (all methods):

         keval = 1;  new = true

         do k = 1, mk

            if (snew(k) < s_total_old) then

               call plscrv3d (mk, xold, yold, zold, sold, method, new, false,  &
                              snew(k), keval, xnew(k), ynew(k), znew(k), derivs)
               new = false

            else ! Only linear extrapolation is likely to be safe

               r = (snew(k) - s_total_old) / (s_total_old - sold(mk-1))
               xnew(k) = xold(mk) + r * (xold(mk) - xold(mk-1))
               ynew(k) = yold(mk) + r * (yold(mk) - yold(mk-1))
               znew(k) = zold(mk) + r * (zold(mk) - zold(mk-1))

            end if

         end do

         do k = 1, mk
            x(k,j,i) = xnew(k)
            y(k,j,i) = ynew(k)
            z(k,j,i) = znew(k)
         end do

      end do ! Next j

   end do ! Next i

   deallocate (xold, yold, zold, sold, xnew, ynew, znew, snew)

!  Local procedure for subroutine redistribute:

   contains

!     --------------------------------------------------------------------------

      subroutine bridge_to_fixed_boundary ()

!     Blend points kb + 1 : mk with shock-adapted points 1 : kb, using the
!     existing last increment at the outer boundary that is being held fixed.

!     --------------------------------------------------------------------------

      integer :: npts
      real    :: d1, d2

      npts     = mk - kb + 2       ! Overlap by one cell
      d1       = snew(kb) - snew(kb-1)
      d2       = sold(mk) - sold(mk-1)
      snew(mk) = sold(mk)

      call blgrid (npts, d1, d2, 2, growth_rate, snew(kb-1), lunerr, ier)

      end subroutine bridge_to_fixed_boundary

   end subroutine redistribute

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine shock_location (method_edge, fs_scale, fs_value_in,              &
                              iq, nf, mi, mj, mk, f, i, j, s, g,               &
                              kedge, ds_local, kpeak, fpeak, s_shock)
!
!  Detect the bow shock location for the indicated radial grid line of a grid
!  block, using the indicated method.  The arc-length distance returned should
!  be adjusted and smoothed at the higher level.  If no edge is detected, the
!  returned value of s_shock will be the existing location, s(mk).
!
!  06, 2005  DAS/S.Yoon Initial implementation (single-function tailoring).
!  07/19/05       "     Initial exercise of gradient-based method, using Mach #.
!  12/18/06      DAS    Refined the definition of ds_local (linear interp.).
!  01/13/07       "     Handle supersonic flow cases where the exit flow Mach
!                       number can exceed the free stream value: leave the
!                       boundary location alone if it's not (essentially) at
!                       the free-stream value.
!  09/05/07  Todd White Todd introduced a method 4 to use the local peak and
!                       to search outward, for a nozzle flow simulation.
!  10/04/07  DAS/TRW    Made method 4 permanent with fs_value < 0 meaning use
!                       the local peak for both methods 1 and 4.  Note that
!                       method 4's search from the wall also skips the test for
!                       blocks interpreted as being in the wake.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in)  :: method_edge     ! 1 => Mach number-based method;
                                            ! 2 => pressure|density|T method;
                                            ! 3 => relative flow grad. change;
                                            ! 4 => as for 1 but search outwards
                                            !      and skip the wake test;
                                            !      see fs_value
   real,    intent (in)  :: fs_scale        ! 1 => fraction of freestream Mach
                                            !      defining shock location;
                                            ! 2 => multiplier applied to the
                                            !      freestream value as in LAURA;
                                            ! 3 => fraction applied to peak
                                            !      flow gradient as in SAGE;
                                            ! 4 => as for 1 but search outwards
   real,    intent (in)  :: fs_value_in     ! Freestream value of the relevant
                                            ! flow quantity, normally > 0;
                                            ! if negative, use local peak along
                                            ! the k line (method_edge = 1 or 4)
   integer, intent (in)  :: iq              ! Flow variable index to use

   integer, intent (in)  :: nf              ! # flow variables in f(1:nf,*,*,*)

   integer, intent (in)  :: mi, mj, mk      ! Grid block dimensions

   real,    intent (in)  :: f(nf,mk,mj,mi)  ! Flow variable used for edge method

   integer, intent (in)  :: i, j            ! Current radial line

   real,    intent (in)  :: s(mk)           ! Arc lengths for this radial line

   real,    intent (out) :: g(mk)           ! For flow gradient along grid line

   integer, intent (out) :: kedge           ! Outermost grid point found not
                                            ! to be at the freestream flow value
   real,    intent (out) :: ds_local        ! Grid spacing local to s_shock
   integer, intent (out) :: kpeak           ! Index of flow peak in the interval
                                            ! [cutoff * mk, mk]
   real,    intent (out) :: fpeak           ! Corresponding flow peak; track
                                            ! this in order to choose fs_scale
   real,    intent (out) :: s_shock         ! Arc length satisfying the shock
                                            ! detection test (~mid-shock if the
                                            ! value of fs_scale is about 0.97
                                            ! (Mach) or 1.2 (p, rho, T);
                                            ! allow a safety margin at the
                                            ! higher level (normally), but
                                            ! s_shock = s_total on return
                                            ! signals "leave boundary alone"
                                            ! [or expand boundary?  TBD]
!  Local constants:

   real, parameter :: almost_one = 0.98, cutoff = 0.2, one = 1., zero = 0.

!  Local variables:

   integer :: k, kcutoff, kstart
   real    :: dg, ds1, ds2, fs98, fs_multiple, fs_value, gmax, gmax_fraction, r

!  Execution:

   kpeak = 0; fpeak = 1.e+30  ! Some methods don't look for a peak inside shock
   kedge = 1; s_shock = s(mk) ! In case inward searches fail
   ds_local = s_shock - s(mk-1)
   fs_value = fs_value_in

   if (fs_value < zero) then  ! Use the local peak Mach (methods 1 or 4)

      fs_value = maxval (f(iq,2:mk-1,j,i))

   end if

   fs_multiple = fs_scale * fs_value

   if (method_edge == 4) then  ! Skip the search for kedge

      ! Set kedge to chosen k below

   else  ! Locate the outermost point not at [almost] freestream, for possible
         ! use in safeguarding the shock location estimate:

      if (fs_scale < one) then
         fs98 = fs_value * almost_one
         do k = mk - 1, 2, -1
            if (f(iq,k,j,i) < fs98) then
               kedge = k + 1
               exit
            end if
         end do
      else
         fs98 = fs_value / almost_one
         do k = mk - 1, 2, -1
            if (f(iq,k,j,i) > fs98) then
               kedge = k + 1
               exit
            end if
         end do
      end if

      if (kedge == mk) go to 99  ! Must be an out-flow boundary; leave it alone

   end if

   select case (method_edge)

      case (1) ! Mach-based shock leading edge method (fs_scale ~0.97 (?))

!        After-thought prompted by a supersonic case where the out-flow Mach
!        number can exceed the free-stream Mach number:

         if (f(iq,mk,j,i) < fs_value * 0.999999 .or. &
             f(iq,mk,j,i) > fs_value * 1.000001) then ! Treat as in the wake
            kedge   =   mk
            s_shock = s(mk)
            go to 99
         end if

         do k = mk - 1, 2, -1
            if (f(iq,k,j,i) < fs_multiple) then
               r = (fs_multiple - f(iq,k,j,i)) / (f(iq,k+1,j,i) - f(iq,k,j,i))
               s_shock  = (one - r) * s(k) + r * s(k+1)
!!!            ds_local = s(k+1) - s(k)  ! Discontinuous; use linear interp.
               ds1 = s(k) - s(k-1)
               ds2 = s(k+1) - s(k)
               r = (s_shock - s(k)) / ds2
               ds_local = (one - r) * ds1 + r * ds2
               exit
            end if
         end do

      case (2) ! LAURA-like: significant change in flow (pressure, density or T)

!        Revised strategy:  locate the flow peak in [0.2 * s(mk), s(mk)] and
!        search outward from there for the indicated multiple of freestream.
!        This search from the true peak avoids spurious lesser peaks observed
!        too near the outer boundary in grid blocks at the far downstream.
!        Safeguard cases that never reach the indicated multiple, but track
!        the lowest peak over all (relevant) grid blocks so that a better
!        multiple can be specified for a rerun.

!        Avoid the region too near the wall by starting the search some
!        reasonable fractional distance above the wall.  (Use of kcutoff itself
!        here can encounter spurious peaks inside the desired shock peak.)
!        LOOSE END:  What if there is no wall?

         kcutoff = 2 * int (real (mk) * cutoff)
         kstart  = kcutoff

         call interval (mk, s, cutoff * s(mk), one, kstart)

         if (kstart >= kedge) kstart = max (kstart / 2, 1)

         fpeak = zero
         do k = kstart, kedge
            if (f(iq,k,j,i) > fpeak) then
               fpeak = f(iq,k,j,i)
               kpeak = k
            end if
         end do

         fs_multiple = fs_scale * fs_value
         if (fs_multiple > fpeak) fs_multiple = fpeak * 0.999

         do k = kpeak + 1, mk
            if (f(iq,k,j,i) <= fs_multiple) then
               r = (fs_multiple - f(iq,k,j,i)) / (f(iq,k-1,j,i) - f(iq,k,j,i))
               s_shock =  r * s(k-1) + (one - r) * s(k)
               exit
            else if (f(iq,k,j,i) > f(iq,k-1,j,i)) then ! No such multiple of fs.
               s_shock = s(k-1)
               exit
            end if
         end do

         k = min (k, mk)
         ds_local = s(k) - s(k-1)

      case (3) ! SAGE-like: relative change in flow gradient vs. peak gradient

         gmax = zero
         do k = 2, mk  ! Would central differencing behave better?
            g(k) = abs (f(iq,k,j,i) - f(iq,k-1,j,i)) / s(k)
            gmax = max (g(k), gmax)
         end do

         gmax_fraction = gmax * fs_scale

         do k = mk - 1, 2, -1
            dg = abs (g(k+1) - g(k))
            if (dg > gmax_fraction) then
               r = (gmax_fraction - g(k)) / dg
               r = max (zero, min (one, r)) ! Else wild extrapolation can occur
               s_shock  = (one - r) * s(k) + r * s(k+1)
               ds_local = s(k+1) - s(k)
               exit
            end if
         end do

      case (4) ! Modified Mach-based shock edge method (fs_scale ~0.97 (?)):
               ! skip the test for being in the wake, and search outwards.

         do k = 2, mk - 1
            if (f(iq,k+1,j,i) > fs_multiple) then
               r = (fs_multiple - f(iq,k,j,i)) / (f(iq,k+1,j,i) - f(iq,k,j,i))
               s_shock  = (one - r) * s(k) + r * s(k+1)
!!!            ds_local = s(k+1) - s(k)  ! Discontinuous; use linear interp.
               ds1 = s(k) - s(k-1)
               ds2 = s(k+1) - s(k)
               r = (s_shock - s(k)) / ds2
               ds_local = (one - r) * ds1 + r * ds2

               kedge = k  ! Avoid triggering a count at the higher level

               exit
            end if
         end do

      case default

         write (*, '(/, a, i3)') &
            ' Shock leading edge method not implemented yet:', method_edge
         stop

   end select

99 return ! Single return philosophy

   end subroutine shock_location

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine adjust_point_counts (ib, ni, nj, nk, nradial,                    &
                                   x, y, z, xout, yout, zout,                  &
                                   icontrols, rcontrols)
!
!  This is an optional additional step following normal OUTBOUND processing.
!
!  The given block is redistributed in the radial direction.  The new dimensions
!  must be known at the higher level for output header info. to be written to
!  disk.  Therefore, the radial direction must be assumed to be along k.
!
!  Cell-Reynolds-number-based initial increments, if specified at the higher
!  level, should already be present in the input grid.  Either the existing
!  initial increments are preserved or all spacings are adjusted according to
!  the change in the number of points, depending on the normal OUTBOUND control
!  for the type of distribution prior to changing the point count.
!
!  04/24/06  D. Saunders  Initial implementation (k direction only).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

!  Arguments:

   integer, intent (in) :: &
      ib,                  &                        ! Current block number
      ni, nj, nk,          &                        ! Input block dimensions
      nradial                                       ! Desired new k dimension

   real, intent (in), dimension (ni,nj,nk) :: &
      x, y, z                                       ! Grid block coordinates

   real, intent (out), dimension (ni,nj,nradial) :: &
      xout, yout, zout                              ! Redistributed coordinates

   integer, intent (in) :: &
      icontrols(*)                                  ! Integer controls;
                                                    ! see active mnemonics below
   real,    intent (in) :: &
      rcontrols(*)                                  ! Real controls;
                                                    ! see active mnemonics below
!  Local constants:

   integer,   parameter :: lunerr = 6
   real,      parameter :: one = 1., zero = 0.
   logical,   parameter :: false  = .false., true = .true.
   character, parameter :: method * 1 = 'M' ! Monotonic spline fits safeguard
                                            ! the symmetry plane
!  Local variables:

   integer :: i, ier, ip, j, k, keval, method_local, method_radial, ngeometric
   real    :: d1, d2, ds2_fraction, growth_rate, p, r, rip, s_total
   real    :: derivs(3)
   logical :: new
   real, allocatable, dimension (:) :: sold, xold, yold, zold
   real, allocatable, dimension (:) :: snew, xnew, ynew, znew

!  Execution:

   allocate (xold(nk), yold(nk), zold(nk), sold(nk), &
             xnew(nradial), ynew(nradial), znew(nradial), snew(nradial))

   method_radial = icontrols(2)  ! Normal input to OUTBOUND before changing nk
   ngeometric    = icontrols(6)
   growth_rate   = rcontrols(8)
   ds2_fraction  = rcontrols(9)

   if (method_radial == 0) then  ! Same relative spacing/different arc lengths
      method_local = 0           ! Same relative spacing/different # points
      r = real (nk - 1) / real (nradial - 1)
   else
      method_local = 2           ! Retain input d1s; adjust d2s appropriately
   end if

   do j = 1, nj

      do i = 1, ni

         do k = 1, nk
            xold(k) = x(i,j,k)
            yold(k) = y(i,j,k)
            zold(k) = z(i,j,k)
         end do

         call chords3d (nk, xold, yold, zold, false, s_total, sold)

         select case (method_local) ! For new radial arc lengths

            case (0) ! Same relative distribution

               snew(1) = zero
               do k = 2, nradial - 1
                  rip = one + r * real (k - 1)
                  ip  = int (rip)
                  p   = rip - real (ip)
                  snew(k) = (one - p) * sold(ip) + p * sold(ip+1)
               end do
               snew(nradial) = s_total

            case (2) ! Existing d1 everywhere

               d1 = sold(2) - sold(1)

!              One-sided Vinokur stretching:

               call expdis5 (1, zero, s_total, d1, nradial, snew, -lunerr)

               d2 = (s_total - snew(nradial-1)) * ds2_fraction

               call blgrid (nradial, d1, d2, ngeometric, growth_rate, snew,    &
                            lunerr, ier)

            case default ! Can't happen

         end select

!        Interpolate the (x,y,z)s as functions of the new arc lengths:

         keval = 1;  new = true

         do k = 1, nradial

            call plscrv3d (nk, xold, yold, zold, sold, method, new, false,     &
                           snew(k), keval, xnew(k), ynew(k), znew(k), derivs)
            new = false
            xout(i,j,k) = xnew(k)
            yout(i,j,k) = ynew(k)
            zout(i,j,k) = znew(k)

         end do

      end do ! Next i

   end do ! Next j

   deallocate (xold, yold, zold, sold, xnew, ynew, znew, snew)

   end subroutine adjust_point_counts

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine scan_hyperbolic_grid (lung, formatted, nblocks, grid, &
                                    xyz_at_xmin, data_range,        &
                                    right_half, left_half, ios)

!  Assume that the grid has been generated via hyperbolic growth from the
!  surface (k = 1) with k = nk at the outer boundary for all blocks, and x
!  points downstream.  Scan all blocks for the outer boundary point furthest
!  upstream, and for the maximum data range in y and z needed for applying the
!  same compression data to different-size cases.  Process one block at a time.
!
!  03/24/11  David Saunders  Initial implementation.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Derived data type and PLOT3D-type I/O module:
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use grid_block_structure  ! Must be compatible with the xyzq_io package
   use xyzq_io_module        ! Compile these two first

   implicit none

!  Arguments:

   integer, intent (in)  :: lung                  ! Grid ready to read block 1,
   logical, intent (in)  :: formatted             ! ... and rewound upon return
   integer, intent (in)  :: nblocks               ! Known number of grid blocks
   type (grid_type), intent (inout) :: grid(nblocks)     ! Array of grid blocks
   real,    intent (out) :: xyz_at_xmin(3)        ! Outer boundary pt. at min. x
   real,    intent (out) :: data_range            ! Maximum y or z data range of
                                                  ! the outer boundary
   logical, intent (out) :: right_half, left_half ! Both halves are left to luck
   integer, intent (out) :: ios                   ! 0 => no allocate/read error

!  Local constants:

   real, parameter :: large = 1.e+9, small = 0.01

!  Local variables:

   integer :: i, ib, j, ni, nj, nk, npts
   real    :: xmin, ymax, ymin, zmax, zmin, y_at_xmin, z_at_xmin

!  Execution:

   xmin = large
   ymin = large;  ymax = -large
   zmin = large;  zmax = -large

   do ib = 1, nblocks

      ni = grid(ib)%ni
      nj = grid(ib)%nj
      nk = grid(ib)%nk

!     Read the current grid block:

      call xyz_allocate (grid(ib), ios)

      if (ios /= 0) then
         write (*, '(/, a, i4, a, i4)') &
            ' Trouble allocating input grid block #', ib, '.  I/O status:', ios
         write (*, '(a, 3i5)') ' Dimensions: ', ni, nj, nk
         go to 99
      end if

      npts = ni * nj * nk

      call xyz_block_io (1, lung, formatted, npts, grid(ib)%x, grid(ib)%y, &
                         grid(ib)%z, ios)
      if (ios /= 0) then
         write (*, '(/, a, i4, a, i4)') ' Trouble reading grid block #',   &
            ib, '.  I/O status:', ios
         write (*, '(a, 3i5)') ' Dimensions: ', ni, nj, nk
         go to 99
      end if

      do j = 1, nj
         do i = 1, ni
            if (xmin      > grid(ib)%x(i,j,nk)) then
                xmin      = grid(ib)%x(i,j,nk)
                y_at_xmin = grid(ib)%y(i,j,nk)
                z_at_xmin = grid(ib)%z(i,j,nk)
            end if
         end do
      end do

      ymax = max (ymax, maxval (grid(ib)%y))
      ymin = min (ymin, minval (grid(ib)%y))
      zmax = max (zmax, maxval (grid(ib)%z))
      zmin = min (zmin, minval (grid(ib)%z))

      deallocate (grid(ib)%x, grid(ib)%y, grid(ib)%z)

   end do  ! Next grid block

   xyz_at_xmin(1) = xmin
   xyz_at_xmin(2) = y_at_xmin
   xyz_at_xmin(3) = z_at_xmin

   right_half = ymax >  small
   left_half  = ymin < -small
   data_range = max (ymax - ymin, zmax - zmin)

   write (*, '(a, 3f16.7)') 'Boundary point furthest upstream: ', xyz_at_xmin
   write (*, '(a,  f16.7)') 'Max. y/z data range:              ', data_range

   rewind (lung)

99 return

   end subroutine scan_hyperbolic_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine compress_grid (lung, formatted, nblocks, grid, xyz_at_xmin, &
                             data_range, right_half, left_half,           &
                             icontrols, rcontrols, fudge, lungrd,         &
                             formatted_out, ios)

!  Shorten radial lines of an initial hyperbolic grid under the assumptions
!  stated in subroutine scan_hyperbolic_grid.  Process one block at a time.
!  For each radial line, apply a controlled fraction to its length based on the
!  distance d of its initial outermost point from the initial point furthest
!  upstream identified by the initial scanning.  See the ancillary program
!  COMPRESSION_DATA for the construction of an empirical dataset from an initial
!  + aligned case.  The dataset is used here to control the compression for a
!  related case, possibly different shape and/or size and/or topology.
!
!  The hope is that the outer boundary will remain smooth enough for the flow
!  solver to proceed with a solution then iteratively move the outer boundary
!  towards the shock in the usual way, having had a better grid to start with.
!
!  The standard ds1 >= 0 options for radial distribution method 2 are provided:
!  existing wall spacing or specified spacing, one- or two-sided stretching.
!
!  03/25/11  David Saunders  Initial implementation.
!  04/12/11    "      "      Switched from a 2-space transformation to 3-space
!                            by introducing z - zu as the third coordinate.
!  08/06/13    "      "      All ADT variants are now in a single module with
!                            generic build_adt and search_adt interfaces.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Derived data type and PLOT3D-type I/O module:
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use grid_block_structure  ! Must be compatible with xyzq_io & adt packages
   use adt_utilities         ! All variants of the ADT build & search routines
   use xyzq_io_module        ! Compile these two first

   implicit none

!  Arguments:

   integer, intent (in)  :: lung                  ! Grid ready to read block 1
   logical, intent (in)  :: formatted             ! T/F
   integer, intent (in)  :: nblocks               ! Number of grid blocks
   type (grid_type), intent (inout) :: grid(nblocks)!Array of grid blocks
   real,    intent (in)  :: xyz_at_xmin(3)        ! Outer boundary pt. at min. x
   real,    intent (in)  :: data_range            ! Max. y or z range of grid
   logical, intent (in)  :: right_half, left_half ! Both halves are left to luck
   integer, intent (in)  :: icontrols(*)          ! Integer controls;
                                                  ! see active mnemonics below
   real,    intent (in)  :: rcontrols(*)          ! Real controls;
                                                  ! see active mnemonics below
   real,    intent (in)  :: fudge                 ! Applied to scale factors to
                                                  ! vary the compression; use
                                                  ! trial and error; fudge = 1.0
                                                  ! would recreate the case used
                                                  ! to generate the compr. data;
                                                  ! > 1.0 => less compression
   integer, intent (in)  :: lungrd                ! Unit for output grid ...
   logical, intent (in)  :: formatted_out         ! ... and its formatting
   integer, intent (out) :: ios                   ! 0 => no allocate/read error

!  Local constants:

   integer,   parameter :: luncrt = 6       ! For diagnostics
   integer,   parameter :: luncgrid = 10    ! For compression data grid ...
   integer,   parameter :: luncfun  = 11    ! ... and compression function
   real,      parameter :: eps = 1.e-7      ! To safeguard atan2
   real,      parameter :: one = 1.
   real,      parameter :: zero = 0.
   logical,   parameter :: false = .false.
   logical,   parameter :: true  = .true.
   character, parameter :: method * 1 = 'M' ! Monotonic spline fits safeguard
                                            ! the symmetry plane
!  Local variables:

   integer :: i, ib, ic, ier, iquad, j, jc, k, keval, n, nbc, ni, nj, nk
   integer :: ngeometric, ninside, noutside, np, npts, nquad, nshape, numf
   real    :: added_margin, d1, d2, dmax, dmean, ds1, ds2_fraction, dsqmin
   real    :: dtolsq, dx, dy, dz, fraction, fraction_max, fraction_min, fscale
   real    :: growth_rate, p, pi, piby2, pm1, q, qm1, r, rscale, s_total, theta
   real    :: xu, yu, zu
   real    :: derivs(3), interp_xyz(3), target_xyz(3)
   logical :: cell_centered, new, two_sided

   integer, allocatable :: conn(:,:)  ! For patch & (i,j) of "surface" quads

   real, allocatable, dimension (:) :: sold, xold, yold, zold
   real, allocatable, dimension (:) :: snew, xnew, ynew, znew

   type (grid_type), pointer :: blockc(:)  ! For compression surface dataset

!  Execution:

   open (luncgrid, file='compression_data.g', status='old', iostat=ios)
   if (ios /= 0) then
      write (luncrt, '(/, a)') ' Trouble opening compression_data.g.'
      go to 99
   end if

   open (luncfun, file='compression_data.f', status='old', iostat=ios)
   if (ios /= 0) then
      write (luncrt, '(/, a)') ' Trouble opening compression_data.f.'
      go to 99
   end if

   call xyzq_read (luncgrid, luncfun, true, nbc, numf, cell_centered, blockc,  &
                   ios)
   if (ios /= 0) go to 99

   fraction_max = zero;  fraction_min = one / eps
   nquad = 0;  np = 0

   do ib = 1, nbc
      ni = blockc(ib)%ni
      nj = blockc(ib)%nj
      np = ni * nj + np
      nquad = (ni - 1) * (nj - 1) + nquad
      fraction_max = max (maxval (blockc(ib)%q), fraction_max)
      fraction_min = min (minval (blockc(ib)%q), fraction_min)
   end do

   write (luncrt, '(/, a, 2f11.6)') &
      'Range of scale factors found in compression data:', &
      fraction_min, fraction_max

   fscale = fudge             ! Fudge = 1 would reproduce the aligned grid from
                              ! from the initial grid used for compression data
   rscale = one / data_range  ! Converts input grid values to compr. data units
   xu = xyz_at_xmin(1)
   yu = xyz_at_xmin(2)
   zu = xyz_at_xmin(3)

!  Construct a search tree from all patches of the compression data:

   allocate (conn(3,nquad)) ! For patch # and (i,j)

   call build_adt (nbc, blockc, nquad, conn)

   piby2    = atan2 (one, zero)
   pi       = piby2 + piby2
   dtolsq   = (0.001)**2  ! Tolerance for search diagnostics
   dmax     = zero
   dmean    = zero
   ninside  = 0
   noutside = 0

   ngeometric   = icontrols(6)
   ds1          = rcontrols(7)
   growth_rate  = rcontrols(8)
   ds2_fraction = rcontrols(9)
   added_margin = rcontrols(10)   ! Multiple of the outermost hyperbolic grid ds
   two_sided    = ds2_fraction /= one

   nk = grid(1)%nk

   allocate (xold(nk), yold(nk), zold(nk), sold(nk), &
             xnew(nk), ynew(nk), znew(nk), snew(nk))

   do ib = 1, nblocks

      ni = grid(ib)%ni
      nj = grid(ib)%nj

!     Read the current grid block:

      call xyz_allocate (grid(ib), ios)

      if (ios /= 0) then
         write (luncrt, '(/, a, i4, a, i4)') &
            ' Trouble allocating input grid block #', ib, '.  I/O status:', ios
         write (luncrt, '(a, 3i5)') ' Dimensions: ', ni, nj, nk
         go to 99
      end if

      npts = ni * nj * nk

      call xyz_block_io (1, lung, formatted, npts, grid(ib)%x,                 &
                         grid(ib)%y, grid(ib)%z, ios)
      if (ios /= 0) then
         write (luncrt, '(/, a, i4, a, i4)') ' Trouble reading grid block #',  &
            ib, '.  I/O status:', ios
         write (luncrt, '(a, 3i5)') ' Dimensions: ', ni, nj, nk
         go to 99
      end if

      do j = 1, nj

         do i = 1, ni

!           Radial arc lengths are needed regardless:

            do k = 1, nk
               xold(k) = grid(ib)%x(i,j,k)
               yold(k) = grid(ib)%y(i,j,k)
               zold(k) = grid(ib)%z(i,j,k)
            end do

            call chords3d (nk, xold, yold, zold, false, s_total, sold)

            d1 = sold(2) - sold(1)
            d2 = s_total - sold(nk-1)

            dx = xold(nk) - xu
            dy = yold(nk) - yu
            dz = zold(nk) - zu

            r  = sqrt (dx**2 + dy**2 + dz**2)

!           Clean symmetry planes (y = 0) are highly recommended - the author's
!           reflect_blocks utility can help.
!           Even then, atan2 needs help near the symmetry plane, but the
!           both-halves case is still problematic - just take what we get:

            if (r < eps) then
               theta = zero
            else
               theta = atan2 (dy, dz)
               if (abs (dy) < eps) then
                  if (dz > zero) then
                     theta = zero
                  else if (dz == zero) then
                     if (right_half) then
                        theta =  piby2
                     else if (left_half) then
                        theta = -piby2
                     end if
                  else  ! dz < zero and dy is tiny
                     if (right_half) then
                        theta =  pi
                     else if (left_half) then
                        theta = -pi
                     end if
                  end if
               end if
            end if

            target_xyz(1) = rscale * r
            target_xyz(2) = theta
            target_xyz(3) = rscale * dz

            call search_adt (target_xyz, iquad, p, q, dsqmin, true,            &
                             nbc, blockc, nquad, conn, interp_xyz)

            if (dsqmin < dtolsq) then ! The nearest quad was within tolerance
               ninside  = ninside + 1
            else
               noutside = noutside + 1
            end if

            dmax  = max (dmax, dsqmin)
            dmean = dmean + dsqmin
            n     = conn(1,iquad) ! Block #
            ic    = conn(2,iquad) ! Lower left quad. indices
            jc    = conn(3,iquad)
            pm1   = one - p
            qm1   = one - q

            fraction = &
               qm1 * (pm1 * blockc(n)%q(1,ic,  jc,  1)  + &
                        p * blockc(n)%q(1,ic+1,jc,  1)) + &
                 q * (pm1 * blockc(n)%q(1,ic,  jc+1,1)  + &
                        p * blockc(n)%q(1,ic+1,jc+1,1))

!!          if (ib==5 .and. j==49) then
!!             write (6, '(a,5i4,6es19.10)') &
!!                'i,j,n,ic,jc,r,theta,p,q,fraction:', i,j,n,ic,jc,&
!!                target_xyz,p,q,fraction*fscale
!!          end if

            s_total  = fscale * fraction * s_total + added_margin * d2

!           One-sided Vinokur stretching:

            if (ds1 > zero) d1 = ds1

            call expdis5 (1, zero, s_total, d1, nk, snew, -luncrt)

            if (two_sided) then

               d2 = (s_total - snew(nk-1)) * ds2_fraction

               call blgrid (nk, d1, d2, ngeometric, growth_rate, snew,         &
                            luncrt, ier)
            end if

!           Interpolate the (x,y,z)s as functions of the new arc lengths:

            keval = 1;  new = true

            do k = 1, nk

               call plscrv3d (nk, xold, yold, zold, sold, method, new, false,  &
                              snew(k), keval, xnew(k), ynew(k), znew(k), derivs)
               new = false
               grid(ib)%x(i,j,k) = xnew(k)
               grid(ib)%y(i,j,k) = ynew(k)
               grid(ib)%z(i,j,k) = znew(k)

            end do

         end do ! Next i

      end do ! Next j

      call xyz_block_io (2, lungrd, formatted_out, npts,                       &
                            grid(ib)%x, grid(ib)%y, grid(ib)%z, ios)
      if (ios /= 0) then
         write (luncrt, '(/, a, i4, a, i4)') &
            ' Trouble writing output grid block #', ib, '  I/O status:', ios
         write (luncrt, '(a, 3i5)') ' Dimensions: ', ni, nj, nk
         go to 99
      end if

      deallocate (grid(ib)%x, grid(ib)%y, grid(ib)%z)

      write (luncrt, '(a, i4)') ' Completed block', ib

   end do  ! Next grid block

   deallocate (xold, yold, zold, sold, xnew, ynew, znew, snew)

   do ib = 1, nbc
      deallocate (blockc(ib)%x, blockc(ib)%y, blockc(ib)%z, blockc(ib)%q)
   end do

   write (luncrt, '(a, 2i7, a, 2es12.5)') &
      ' # points in/outside tolerance:', ninside, noutside, &
      ';  max/mean devn.:', sqrt (dmax), sqrt (dmean / real (np))

99 return

   end subroutine compress_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine scale_grid (ndim, lung, formatted, nblocks, grid,        &
                          icontrols, rcontrols, scale_factor, lungrd,  &
                          formatted_out, ios)

!  Scale radial grid lines by the same factor everywhere as follows:
!
!             x(k) <-- x(1) + scale_factor * (x(k) - x(1))
!             y(k) <-- y(1) + scale_factor * (y(k) - y(1))
!            [z(k) <-- z(1) + scale_factor * (z(k) - z(1))]
!
!  If scale_factor < 0, the scaling at each k is adjusted as follows:
!
!             dx        = x(k) - x(1)             [and likewise for dy, dz]
!             dstraight = sqrt (dx^2 + dy^2 + dz^2)
!             factor    = -scale_factor * (dstraight / sold(k))
!             xnew(k)   = x(1) + factor * dx      [and likewise for ynew, znew]
!
!  Neither option preserves the initial radial lines unless they are perfectly
!  straight.
!
!  The resulting points are redistributed with the indicated option, such as
!  retaining the input initial spacings with 1- or 2-sided stretching.
!  For fairly large changes, this may be preferable to using the added margin
!  option (multiple of outer spacings, added or subtracted), but neither of
!  these is likely to be useful for large changes.
!
!  03/31/11  David Saunders  Initial implementation, after Ryan McDaniel
!                            pointed out that this is what SAGE must do.
!  12/12/13    "      "      Handled 2-D option.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Derived data type and PLOT3D-type I/O modules:
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use grid_block_structure
   use xyq_io_module
   use xyzq_io_module

   implicit none

!  Arguments:

   integer, intent (in)  :: ndim                  ! 2 or 3
   integer, intent (in)  :: lung                  ! Grid ready to read block 1
   logical, intent (in)  :: formatted             ! T/F
   integer, intent (in)  :: nblocks               ! Number of grid blocks
   type (grid_type), intent (inout) :: grid(nblocks)!Array of grid blocks
   integer, intent (in)  :: icontrols(*)          ! Integer controls;
                                                  ! see active mnemonics below
   real,    intent (in)  :: rcontrols(*)          ! Real controls;
                                                  ! see active mnemonics below
   real,    intent (in)  :: scale_factor          ! Fraction to be applied to
                                                  ! offsets from wall points;
                                                  ! if scale_factor < 0., see
                                                  ! the variation detailed above
   integer, intent (in)  :: lungrd                ! Unit for output grid ...
   logical, intent (in)  :: formatted_out         ! ... and its formatting
   integer, intent (out) :: ios                   ! 0 => no allocate/read error

!  Local constants:

   integer,   parameter :: luncrt = 6       ! For diagnostics
   real,      parameter :: one = 1.
   real,      parameter :: zero = 0.
   logical,   parameter :: false = .false.
   logical,   parameter :: true  = .true.
   character, parameter :: method * 1 = 'M' ! Monotonic spline fits safeguard
                                            ! the symmetry plane
!  Local variables:

   integer :: i, ib, ier, j, k, keval, ni, nj, nk, ngeometric, np, npts
   real    :: added_margin, d1, d2, dstraight, ds2_fraction, dx, dy, dz
   real    :: factor, growth_rate, s_total, x1, y1, z1, derivs(3)
   logical :: new, two_sided

   real, allocatable, dimension (:) :: sold, xold, yold, zold
   real, allocatable, dimension (:) :: snew, xnew, ynew, znew

!  Execution:

   ngeometric   = icontrols(6)
   growth_rate  = rcontrols(8)
   ds2_fraction = rcontrols(9)
   added_margin = rcontrols(10)   ! Applied to the outermost input grid ds
   two_sided    = ds2_fraction /= one

   if (ndim == 2) then
      np = grid(1)%nj
      nk = 1  ! For an error message
   else
      np = grid(1)%nk
      nk = np
   end if

   allocate (xold(np), yold(np), zold(np), sold(np), &
             xnew(np), ynew(np), znew(np), snew(np))

   do ib = 1, nblocks

      ni = grid(ib)%ni
      nj = grid(ib)%nj

!     Read the current grid block:

      if (ndim == 2) then
         call xy_allocate  (grid(ib), ios)
         npts = ni * nj
         call xy_block_io (1, lung, formatted, npts, grid(ib)%x, &
                           grid(ib)%y, ios)
         allocate (grid(ib)%z(ni,nj,1))
         grid(ib)%z(:,:,:) = zero
      else
         call xyz_allocate (grid(ib), ios)
         npts = ni * nj * nk
         call xyz_block_io (1, lung, formatted, npts, grid(ib)%x, &
                            grid(ib)%y, grid(ib)%z, ios)
      end if

      if (ios /= 0) then
         write (luncrt, '(/, a, i4, a, i4)') ' Trouble reading grid block #',  &
            ib, '.  I/O status:', ios
         write (luncrt, '(a, 3i5)') ' Dimensions: ', ni, nj, nk
         go to 99
      end if

      factor = scale_factor  ! Unless the input value was negative

      if (ndim == 2) then  ! Makng it look like 3-D is least messy
         nk = nj
         nj = 1
      end if

      do j = 1, nj

         do i = 1, ni

            if (ndim == 2) then
               xold(:) = grid(ib)%x(i,:,1)  ! Need to handle subscript checking
               yold(:) = grid(ib)%y(i,:,1)
               zold(:) = grid(ib)%z(i,:,1)
            else
               xold(:) = grid(ib)%x(i,j,:)
               yold(:) = grid(ib)%y(i,j,:)
               zold(:) = grid(ib)%z(i,j,:)
            end if

            call chords3d (nk, xold, yold, zold, false, s_total, sold)

            d1 = sold(2) - sold(1)
            d2 = s_total - sold(nk-1)

            x1 = xold(1)
            y1 = yold(1)
            z1 = zold(1)

!           Geometric scaling:

            do k = 2, nk
               dx = xold(k) - x1
               dy = yold(k) - y1
               dz = zold(k) - z1
               if (scale_factor < zero) then
                  dstraight = sqrt (dx**2 + dy**2 + dz**2)
!!                factor    = -scale_factor * (sold(k) / dstraight)
                  factor    = -scale_factor * (dstraight / sold(k))  ! Better
               end if
               xold(k)   = x1 + factor * dx
               yold(k)   = y1 + factor * dy
               zold(k)   = z1 + factor * dz
            end do

            call chords3d (nk, xold, yold, zold, false, s_total, sold)

            s_total = s_total + added_margin * d2

!           One-sided Vinokur stretching:

            call expdis5 (1, zero, s_total, d1, nk, snew, -luncrt)

            if (two_sided) then

               d2 = (s_total - snew(nk-1)) * ds2_fraction

               call blgrid (nk, d1, d2, ngeometric, growth_rate, snew,         &
                            luncrt, ier)
            end if

!           Redistribute the pts. with the usual arc-length-based interpolation:

            keval = 1;  new = true

            do k = 1, nk

               call plscrv3d (nk, xold, yold, zold, sold, method, new, false,  &
                              snew(k), keval, xnew(k), ynew(k), znew(k), derivs)
               new = false

               if (ndim == 2) then
                  grid(ib)%x(i,k,1) = xnew(k)
                  grid(ib)%y(i,k,1) = ynew(k)
               else
                  grid(ib)%x(i,j,k) = xnew(k)
                  grid(ib)%y(i,j,k) = ynew(k)
                  grid(ib)%z(i,j,k) = znew(k)
               end if

            end do

         end do ! Next i

      end do ! Next j

      if (ndim == 2) then
         call xy_block_io  (2, lungrd, formatted_out, npts, &
                            grid(ib)%x, grid(ib)%y, ios)
         nj = nk
         nk = 1
      else
         call xyz_block_io (2, lungrd, formatted_out, npts, &
                            grid(ib)%x, grid(ib)%y, grid(ib)%z, ios)
      end if

      if (ios /= 0) then
         write (luncrt, '(/, a, i4, a, i4)') &
            ' Trouble writing output grid block #', ib, '  I/O status:', ios
         write (luncrt, '(a, 3i5)') ' Dimensions: ', ni, nj, nk
         go to 99
      end if

      deallocate (grid(ib)%x, grid(ib)%y, grid(ib)%z)

      write (luncrt, '(a, i4)') ' Completed block', ib

   end do  ! Next grid block

   deallocate (xold, yold, zold, sold, xnew, ynew, znew, snew)

99 return

   end subroutine scale_grid
