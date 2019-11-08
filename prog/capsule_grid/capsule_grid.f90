!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   program capsule_grid
!
!  CAPSULE_GRID combines the functions of its precursors, HEAT_SHIELD+FOREBODY_-
!  REGRID and AFT_BODY, for axisymmetric capsules (only).  [AFT_BODY assumes
!  axisymmetry too; HEAT_SHIELD+FOREBODY_REGRID do not.]  CAPSULE_GRID can omit
!  the aft body if only a symmetric forebody surface grid is required.  If both
!  fore- and aft body are being treated, CAPSULE_GRID eliminates some of the
!  steps required when the surfaces are generated separately and turned into a
!  volume grid.  If only an aft body is to be treated, use AFT_BODY instead.
!
!  Note:  Both forebody-only and fore+aft-body cases now have an "umbrella"
!         option to model flexible TPS with folding ribs.  They both start,
!         of course, with axisymmetric spoked surface grids that are regridded
!         to model the faceting with optional deflections, rib-rounding, etc.
!
!  Input Coordinate System (2-space):
!
!     X points downstream for zero angle of attack
!     R points up
!
!  Output Coordinate System (Right-Handed):
!
!     X points downstream for zero angle of attack
!     Y points to the left looking from upstream
!     Z is up
!
!  As with AFT_BODY, the generatrix can be defined analytically by a handful of
!  inputs (then discretized), or it can be read from a file.  Either way, the
!  generatrix is used to make an interim body-of-revolution surface grid with
!  a singular point at the nose (and opposite end if the aft-body is present),
!  then any singular point region is replaced with mapped quadrant patches.
!
!  Analytic forebodies are limited to sphere-cones, biconics, or Apollo-type
!  spherical sections.  Since the discretization is close to uniform (while
!  capturing tangency points, etc., precisely), the forebody generatrix starts
!  with about twice as many points as specified for the output forebody grid,
!  then it is reduced to numi_forebody via curvature-based redistribution,
!  which should produce good resolution of a rounded shoulder.
!
!  A "ripple" option is provided to simulate skin flexing between underlying
!  toroids of inflatable concepts.  Limited to sphere-cone forebodies, this
!  option is specified by a number of equal toroids positioned between the
!  rigid nose cap and the shoulder toroid (which may be of different size), and
!  a peak deflection between the toroids.  A sinuisoidal approximation is made
!  to give the ripple effect, with only small deflections intended.  The input
!  number of toroids, ntoroids, does not include the shoulder toroid.
!
!  Analytic aft bodies are more general (any number of cone angles, with options
!  to round the vertices).  Since aft bodies are either rectilinear or have
!  vertices rounded with the indicated number of local points, the strategy is
!  a little different, but if any rounding is specified, curvature-based
!  redistribution to numi_aft_body points is performed.
!
!  Original Nose-Regridding Strategy (Reused Similarly On the Aft Body):
!
!  >  Retain the existing outer portion of the interim spoked grid, split into
!     subpatches as needed.
!  >  Some radial index defines the roughly circular perimeter of the new
!     patches in the nose region (not necessarily at constant X).  17 or 33
!     are the usual choices depending on grid dimensions.  These are heuristic
!     but are good choices for 101 or 201 points along the input grid spokes
!     and 65 or 121 spokes.  N.B.: the number of spokes minus 1 must be a
!     multiple of 4 for the 2-quadrant nose patch substitution.
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
!     generatrix.  Densifying the interim surface grid by a factor of 16 here
!     helps produce good results.
!
!  Control Inputs (namelist(s) via standard input; all lines start in column 2):
!
!     $FOREBODY_INPUTS
!     verbose = F,                 ! T turns on lots of pt. distribution output
!     aft_body_too = T,            ! Option to suppress the aft body
!     units_in_inches = F,         ! Option suited to arc-jets in the US of A
!     geometric_scale = 1.,        ! Option to scale the indicated geometry
!     input_generatrix = '...',    ! Full generatrix file; 'none' => analytic
!     spherecone = T,              ! T => half_cone_angle >= 0, else biconic
!     numi_forebody = 201,         ! # grid points on forebody generatrix
!     numj = 121,                  ! # spokes in body of revolution (4m + 1)
!     x_nose = 0.,                 ! Nose coordinate
!     r_nose = 0.,                 ! Nose coordinate; must be 0. if aft body on
!     radius_nose = 0.0889,        ! Radius of nose circular segment
!     radius_base = 0.1778,        ! Radius at maximum diameter
!     radius_shoulder = 0.00889,   ! Radius of shoulder circular segment
!     radius_vertex = 0.           ! Biconic option to round the cone juncture
!     half_cone_angle = 55.,       ! 0 => spherical section (Apollo-like)
!     half_cone_angle_fore = ...,  ! Forward cone angle, if biconic
!     half_cone_angle_aft = ...,   ! Aft cone angle, if biconic
!     radius_cone_juncture = ...,  ! Radius where cones meet, if biconic
!     skirt_angle = 35.,           ! Semi-angle if cone segment follows shoulder
!     skirt_length = 0.01,         ! Length in x, >=< 0., but in y if angle=90.
!     nose_patch_file_name = '...' ! Grid quadrant used to remove singular point
!     output_generatrix = '...',   ! Generatrix corresp. to the surface grid
!     surface_grid_file_name = '.. ! Output surface grid/no singular point(s)
!     flat_blend_fore = T,         ! T ==> extra blending between high/low curv.
!     power_fore = 0.5,            ! Initial power in curvature-based redistrbn.
!     ni_regrid_forebody = 33,     ! 1:ni_regrid, 1:numj nose region is replaced
!     ripple_case = ...,           ! T => model n toroids on cone with deflctns.
!     ntoroids = ...,              ! # inflatable toroids to simulate along cone
!     peak_ripple = ...,           ! Peak deflection between toroids
!     umbrella_case = ...,         ! T => polygonal forebody/optional deflectns.
!     nedges = ...,                ! # full polygon edges if umbrella_case
!     peak_deflection = ...,       ! > 0. => catenary-type fabric deflections
!     rib_deflection = ...,        ! Option to bend ribs (+/- => concave/convex)
!     frustum_radius = ...,        ! 0. => start fabric/rib deflns. @ tngcy. pt.
!     resolve_the_ridges = ...,    ! T => adjust spacing either side of ribs
!     rib_thickness = ...,         ! Off-rib spacing to impose (<= input value)
!     round_the_ridges = ...,      ! T => monotonic local cubic spline rounding
!     $END
!
!     $AFT_BODY_INPUTS
!     sting_case = F,              ! T means don't close the aft body
!     ncones = 2,                  ! ncones >= 1
!     r_aft = 99., 0.12, 0.06,     ! Radii from fore/aft juncture (1:ncones+1)
!     cone_angle_aft = 35., 65.,   ! Aft cone half angles (1:ncones); 0, 90 OK
!     rounding_mode = 2            ! 1 => tangent length; 2 => radius; 0 = none
!     ds_round = 0., 0.003, 0.003, ! Tangent lengths or rounding radii >= 0.
!     ni_round = 0, 17, 17,        ! # uniform rounding circle arc points >= 0
!     flat_blend_aft = T,          ! T ==> extra blending between high/low curv.
!     power_aft = 0.5,             ! Initial power in curvature-based redistrbn.
!     numi_aft_body = 201,         ! # grid points on aft body generatrix
!     sting_stretch_multiplier = 2.! Applied to interim end-of-sting spacing
!     ng_split = 0,                ! Split index for full-body input generatrix
!                                  ! 0 = max. radius index; < 0 1-piece/no split
!     nblend_fore_aft = 8,         ! # aft pts. used to match fore/aft spacing
!     ni_regrid_aft_body = 0,      ! Automated if analytic; enter desired index
!     i0_umbrella = 0,             ! For fat aft body cases (extra morphing)
!     i1_umbrella = 0,             ! Inboard index defining aft spoke morphs
!     i2_umbrella = 0,             ! Outboard index near shoulder likewise
!     $END                         ! if automation fails or for input generatrix
!
!  Clarification of Input Generatrix Options:
!
!     Since an input generatrix may or may not represent the desired point
!     distribution in the output surface grid, the options are as follows:
!
!     If the number of points found in the input generatrix matches the spec-
!     ified grid point count, numi_forebody [+ numi_aft_body - 1], then the
!     input points are used without change.  Otherwise, the curvature-based
!     redistribution(s) applied to the analytic case is/are applied likewise,
!     and the number of grid points can be adjusted to suit the intended flow
!     calculation.
!
!     For the full-body input generatrix case, use ng_split > 0 to apply
!     numi_forebody to geometry points 1:ng_split and numi_aft_body to geometry
!     points ng_split:ng.  If ng_split = 0 in this case, it defaults to the
!     index of the input geometry point with maximum radius.
!
!     See the next item for use of ng_split < 0.
!
!  Clarification of Curvature-Based Gridding Controls:
!
!     Fore and aft body generatrices are normally redistributed separately, but
!     there are two exceptions:
!
!        (1) A fore- or full-body input generatrix case where the output number
!     of points, numi_forebody [+ numi_aft_body - 1], matches the input number,
!     as mentioned above: no redistribution at all.
!
!        (2) A full-body input generatrix case where ng_split < 0 is entered.
!     This is used to force redistribution of the full generatrix as one piece
!     so as to avoid blending of fore- and aft-grid spacing across a split.
!     Note that with roughly double the number of points and data range (which
!     is normalized to [0, 1] before the curvature-based calculations), the
!     single-piece full-body generatrix redistribution problem is harder and
!     different from the two-part case, and there is not much control over
!     where the numi_forebody grid point lands other than by changing its
!     value through trial and error.  The gridded forebody portion will not
!     be the same as if it were redistributed separately, because of the
!     different scaling to [0, 1] that affects curvature values nonlinearly.
!
!     Curvature-based redistribution may be controlled with exponent inputs
!     power_fore and power_aft, which default to 0.5 and should be in the range
!     [0.1, 1.0], with 1.0 maximizing the curvature effect.  If the iterative
!     solution does not converge, the exponent is lowered by 0.1 until it does.
!     The extreme case of power = 0. (uniform) should never fail, although it is
!     unlikely to be useful except for plotting purposes.
!
!     For the full-body input generatrix case with ng_split < 0, the forebody
!     control power_fore is used and power_aft is ignored.
!
!  Full Hemisphere Nose:
!
!     This extreme case of the "Apollo" option (half_cone_angle = 0.) is
!     indicated if radius_base = radius_nose on input.  However, only if
!     radius_shoulder = 0. can radius_base actually equal radius_nose.
!     Otherwise, the maximum radius is determined by radius_nose and _shoulder.
!
!  "Umbrella" Forebodies:  Two Types of Deflection
!
!     Originally, faceted outer forebodies were constructed with straight-line,
!     undeflected ribs matching the cone portion of the underlying sphere-cone
!     generatrix, and any deflections were confined to the panels between the
!     ribs.  The "peak_deflection" input controls these catenary-type surface
!     deflections in both surface directions.  Now, a second rib_deflection
!     input allows the ribs to be bent in either direction (negative input for
!     a convex/outwards rib bend).  The initial discretization of the generatrix
!     is bent with a catenary-type deflection either inwards or outwards.
!     Blending at the circular shoulder is an awkwardness treated as follows.
!     The outer end slope of the nominal catenary is used to locate the shoulder
!     circle point with matching slope.  Then the nominal catenary is morphed
!     via an arc-length-based technique to move its outer end from the tangent
!     point of the straight cone line to that matching-slope circle point.  The
!     resulting blend is not perfectly tangential, but should suffice for small
!     deflections.
!
!  Definition of Frustum Radius:
!
!     Originally, "umbrella" cases assumed that the ribs started at the sphere-
!     cone tangency point.  Now, ribs may be shorter than that, to model an
!     effectively larger nose cap that includes part of the cone.  The control
!     input "frustum_radius" defines the start of the ribs where any deflections
!     may be imposed out to the shoulder tangency location.  Defaulting this
!     frustum radius to zero leads to use of the original definition, namely the
!     sphere-cone tangency point.  To be precise, any rib bending will start at
!     the first point of the eventual (redistributed) generatrix that is at or
!     nearest to (but downstream of) the location specified by frustum_radius.
!
!  Clarification of Aft Body ncones and Collapsed Aft Ends:
!
!     To avoid confusion, ncones is expected to be 1 for a single frustrum,
!     2 for a typical biconic aft body, and so on.  Normally, these cases have
!     a finite vertical end segment in the generatrix, where a half cone angle
!     of 90 degrees is implied and ncones+1 radii > 0 should also be entered.
!     Internally, r_aft(ncones+2) = 0 will be assigned for arc length calcs.
!
!     Conversely, if the aft end collapses to a point, r_aft(ncones+1) = 0.
!
!  Aft Body Vertex Rounding Options (Analytic Cases Only):
!
!     [Note that forebody rounding is limited to the one vertex of a biconic,
!     and is handled in the capsules.f90 forebody package via the radius_vertex
!     parameter.  Aft bodies are more general and any rounding of vertices is
!     handled in this program.]
!
!     rounding_mode = 0 turns off any rounding of [aft body] vertices;
!                   = 1 means input ds_round values > 0 are treated as tangent
!                       lengths from the vertex to the underlying circle;
!                   = 2 means input ds_round values > 0 are treated as radii
!                       of the rounding circle; ALSO: ds_round < 0 means a
!                       special case of rounding for SRC-type capsules with a
!                       large circle <= one quadrant for the aft-most rounding.
!
!     If ds_round > 0 is specified for a vertex, a circular arc with ni_round
!     uniform points from tangent point to tangent point is inserted.  After
!     all insertions, the new geometry definition is redistributed to the
!     original numi_aft_body points using curvature-based techniques.
!     A collapsed aft body may also be rounded similarly: it is treated just as
!     for other vertices, but only half the circular points are inserted.
!
!     Note that ds_round(1) and ni_round(1) are ignored (no rounding at the
!     forebody/aft body juncture), but zeros should be entered as place holders.
!
!  Special SRC (Sample Return Capsule)-type Option (Full Quadrant Aft Body):
!
!     Use of ds_round(ncones) < 0. is the flag for a special case where the aft
!     body ends with a full quadrant (or a little less) that leaves no room for
!     normal rounding at the second last vertex.  A vertical segment is also
!     expected to precede the last vertex to be rounded.  Where that meets the
!     large aft circle has to be rounded after the normal rounding, because it
!     involves the vertical segment and what is now a large circular segment,
!     not a straight segment that would remain after less severe rounding of
!     the last vertex.  (Hard to explain!)  If the SRC-type aft body is not a
!     full quadrant, some geometry and trigonometry is required to determine
!     the last two r_aft inputs and the last cone angle implied by rounding
!     of the last vertex with the large radius of the less-than-full quadrant:
!
!        Rb = ds_round(ncones+1)  ! Large aft body radius of near-quadrant
!        xe = x_aft(ncones+2)     ! End of body, presumably known a priori
!        xo = xe - Rb             ! Center of large aft circle
!        xv = x_aft(ncones)       ! Vertical segment abscissa
!        tc = arcsin ((xv-xo)/Rb) ! The angle needed for cone_angle_aft(ncones)
!        ri = Rb*cos(tc)          ! Intersection radius needed for r_aft(ncones)
!
!  Handling of Aft Body 0, 180 and 90 Degree Half Cone Angles:
!
!     Note that if cone_angle_aft(j) = 0 or 180, the analytic specification is
!     not well defined.  The work-around is to enter the length of the step (in
!     the x direction) via r_aft(j+1); this will be reset to r_aft(j) internally
!     by the program.  For the 180 degree case, make r_aft(j+1) negative.
!
!     An interior cone angle of 90 degrees is also permitted without special
!     handling other than avoiding the infinite tangent internally.
!
!  Skirt Angle = Aft Body Half Cone Angle > 90 With Skirt Length < 0 (SRC-type):
!
!     The SRC-type configuration also requires an aft body cone angle > 90, in
!     combination with an equal forebody skirt angle and a skirt length < 0.
!     Remember that skirt length is in the x direction, not along the arc,
!     except if the skirt angle is 90 degrees, when the opposite is true.
!
!  Aft Body Angle > 180 For Open SPRITE/ADEPT-Type With Flexible TPS:
!
!     The forebody/aft body juncture may need to be placed aft of the true
!     forebody shoulder via choice of skirt_length > 0 with skirt_angle = 15
!     say), combined with initial cone angles of 15 (say), 90, and 195 (say).
!     Looking at a working example is recommended.
!
!  Aft Body With "Umbrella" Forebody:
!
!     Default i1_umbrella and i2_umbrella = 0 (or umbrella_case = F) initially.
!     This suppresses the option for morphing an ADEPT-type aft body to an
!     umbrella-type forebody.
!     Pick indices on the resulting axisymmetric aft body spokes that make
!     sense for arc-length-based morphing (defining the essentially straight
!     line portions).  Outboard of that, the shoulder is simply shifted.
!
!     N.B.: These indices should be found from patch 7 of the interim surface
!     grid, NOT from the single-patch spoked_surface_aft.g file.
!     Also: for a fat aft body, if i1_umbrella on the aft body corresponds to
!     a radius that is greater than radius_frustum (where the forebody faceting
!     starts), use 0 < i0_umbrella < i1 to control additional morphing between
!     i0 and i1 because of the aft-body x-shifts from i1:i2 that maintain the
!     fabric thickness.  (Ideally these shifts go to zero at i1, but may not
!     for a fat payload housing.)  Default:  i0_umbrella = i1_umbrella - 20
!
!  Non-Analytic Full Body Case:
!
!     The inputs for numi_forebody and numi_aft_body should allow for a common
!     point where the two portions meet.  E.g., if 401 points are read in the
!     input generatrix, then 201 and 201 are plausible fore- and aft inputs.
!     See also "Clarification of Input Generatrix Option" above concerning the
!     sum of numi_forebody + numi_aft_body - 1.
!
!     If curvature-based redistribution is indicated, this is now done in two
!     parts.  The input geometry is split at the point of maximum radius, but
!     this can be overridden by entering ng_split > 0.  Keep in mind that the
!     split should be at or slightly aft of the shoulder max. diameter, for BC
!     and fore/aft block plotting purposes.
!
!  Sharp Vertices Not To Be Rounded:
!
!     Sharp vertices in a discretized generatrix produce curvature spikes that
!     are effectively invisible to the curvature-based grid point distribution
!     scheme, because they are just one data point wide.  This means the spacing
!     in the vicinity of an intentionally sharp vertex is essentially uniform
!     unless some work-around is devised, as one has been.  The CURVDIS utility
!     now has an option to detect such spikes and broaden them and moderate
!     their height, and this option is invoked by CAPSULE_GRID because it passes
!     ISMOOTH = -3 to CURVDIS[2].  The geometry is not touched -- only the curv-
!     ature, and no new controls are needed (so far).
!
!  History:
!
!     01/09/2008  D.A.Saunders  Initial design of FOREBODY_REGRID.
!     01/10/2008- ELORET Corp.  Initial implementation (quarter circle case
!     01/14/2008                for replacing the singular nose point).
!     12/17/2010  DAS, ERC Inc. Fudged the two poor corner points of the
!                               circular quadrants (only), after pondering
!                               squarer nose patches that aren't really
!                               doable (well) algebraically.
!     10/11/2011 -  "    "      Added the options to impose umbrella-type
!     10/14/2011                shoulder edge faceting & catenary deflections.
!     10/18/2011    "    "      Added the option to round off the sharp ridges.
!     10/20/2011    "    "      Noticed method ='b' should be 'B' (loose fits
!                               in subroutine initial_regrid).  But see 11/03.
!                               Added the option to resolve the rib ridges by
!                               restretching outboard of where the spacing is
!                               more than (half) the specified rib thickness.
!     10/28/2011-   "    "      Initial adaptation of FOREBODY_REGRID as
!     10/31/2011                program AFT_BODY, for the simple cases.
!     11/03/2011    "    "      Actually, 'b' was as intended for ADJUSTN: it
!                               gives 'B' (loose fits) and uniform spacing.
!     11/05/2011-   "    "      [In AFT_BODY:]  Added options to round the
!     11/08/2011                (analytic) vertices, and to adjust juncture
!                               spacing to match the intended forebody.
!     11/12/2011    "    "      If the last aft body segment is rounded
!                               [starting as collapsed] then the regridded
!                               portion should be confined to the rounded
!                               region.  Otherwise, the resulting uniform
!                               spacing is too large for good results in
!                               subsequent volume gridding.
!     11/14/2011-   "    "      Merged FOREBODY_REGRID and AFT_BODY as
!     12/31/2011                CAPSULE_GRID (many interruptions on the way).
!     01/04/2012    "    "      The Gridgen script needs a Z column in the
!                               generatrix (and # points at the top, but
!                               that's incompatible with Tecplot).
!     01/05/2012    "    "      The hard-coded power of 0.5 was found not to
!                               converge for CURVDIS on an unusual aft body, so
!                               iterate with smaller powers until it converges.
!     01/06/2012    "    "      Sample Return Capsule (SRC) cases required a
!                               few tweaks to deal with cone angles > 90 deg.
!     01/09/2012    "    "      Rounding the vertex left by an SRC case with
!                               no room for normal rounding required special
!                               handling after all normal rounding, indicated
!                               by ds_round(ncones) < 0 and other conditions.
!     01/31/2012    "    "      The "umbrella" option needed to swap y and z
!                               to match FOREBODY_REGRID conventions (then
!                               swap back before the nose-cap regridding).
!                               Also, don't allow ribs to be on the nose cap
!                               when their thickness is being resolved.
!     02/22/2012    "    "      Added the option to enter dimensions in inches.
!                               The aft-body regridding to remove the singular
!                               point should not affect the 2D generatrix.
!     02/23/2012    "    "      The non-analytic full-body case needed tweaking
!                               to match an analytic case by reading back its
!                               generatrix.  Keep numi_forebody + numi_aft_body
!                               equal to 1 + numi_generatrix from input.
!     02/27/2012    "    "      Scaling the CEV geometry to wind tunnel size
!                               prompted an input geometric scaling factor.
!     04/16/2012    "    "      If input units are in inches, show any aft-body
!                               vertices in inches as well as meters.  Use of
!                               (2es15.7) (say) is preferable to (1p, 2e15.7).
!     04/18/2012    "    "      Use of aft-body features to simulate an arc-jet
!                               wedge with recessed test article revealed a
!                               failure mode in the forebody spacing matching.
!     04/20/2012    "    "      "Ripple" option added to approximate inflatable
!                               toroid concepts.
!     06/30/2012    "    "      Automated the sting case (suppress closure of
!                               the aft body if sting_case = T).
!     07/01/2012    "    "      Redistribute the sting portion of a sting case
!                               generatrix via sting_stretch_multiplier input.
!     07/27/2012    "    "      Generalized the original SRC special case so
!                               that the large rounded aft-body circle sector
!                               may be less than or equal to a full quadrant.
!     07/31/2012    "    "      The raw aft generatrix data wasn't being saved
!                               if rounding_mode = 0.
!     08/01/2012    "    "      A cone angle of 180 degrees needs to be handled
!                               as for 0 degrees.  Also, ni_regrid_aft_body = 0
!                               may lose curvature-based spacing, so issue an
!                               appropriate warning.
!     08/23/2012    "    "      Arranged to include an aft body with the
!                               "umbrella" option.  New inputs i1_umbrella and
!                               i2_umbrella need to be nonzero to become active
!                               since defaulting is too awkward.  Aft spokes
!                               outboard of i2_umbrella (on the shoulder) are
!                               simply shifted to match the faceted forebody;
!                               between i1_ and i2_umbrella, they are morphed
!                               based on arc length, because that segment is
!                               expected to be (close to) a straight line.
!     09/07/2012    "    "      Suppressing the output generatrix when one was
!                               read was a mistake: the number of points may
!                               need adjusting, and we need what matches the
!                               output surface grid, for use with REVOLVE_GRID.
!                               Also: a proportionally-adjusted form of an
!                               input generatrix is saved for possible use as
!                               an alternative to the curvature-based form.
!                               Further:  WARPQ3D2 has replaced WARPQ3D to
!                               avoid the work-space arguments.
!     09/11/2012    "    "      Arranged for deflections to begin at any index,
!                               not necessarily the default one (slightly out-
!                               board of the sphere-cone tangency point).
!     10/24/2012    "    "      New rib_deflection input controls optional
!                               bending of the cone portion of a sphere-cone
!                               generatrix (concave or convex).
!     11/05/2012    "    "      The aft body option for umbrella cases is
!                               affected by any fabric or rib deflections
!                               or if the underlying forebody is not a sphere-
!                               cone but an arbitrary shape such as Config. Z3.
!                               Fabric/rib thickness is made ~constant via
!                               shifts in x (only) between i1- & i2_umbrella.
!                               Also: the rib rounding option now makes them
!                               flat (tight fit interpolations), not rounded.
!     11/08/2012    "    "      Introduced CURVDIS2, which normalizes the
!                               geometry before calling CURVDIS, in order to
!                               deal with the extreme curvature of tiny shock
!                               tunnel configurations.
!     11/15/2012    "    "      Fat payload housings for ADEPT configurations
!                               with the radius at i1_umbrella > radius_frustum
!                               require additional morphing inboard of i1 to
!                               deal with a non-zero x shift at i1 when the
!                               fabric thickness at i2 is imposed for i1:i2.
!                               This requires an i0_umbrella input (see above).
!     03/05/2012    "    "      The non-analytic full-body case had a glitch:
!                               numi <-- numi_aft_body was missing before the
!                               spoked_surface (2) call.
!     03/07/2012    "    "      Allowed for a nose quadrant patch that isn't
!                               necessarily uniform along the straight edges.
!                               Tightening the spacing slightly towards the
!                               apex may be slightly better.  See specialized
!                               utility MORPH_QUADRANT for adjusting the
!                               uniform-edge quadrant template.
!     04/15/2012    "    "      Actually, DPLR gives better results with uni-
!                               form spacing along the straight nose-patch
!                               edges, so nonuniform is not recommended.
!                               Arranged for the full-body input generatrix
!                               case to apply curvature-based redistribution
!                               separately to the fore and aft portions in the
!                               hope that the smaller nonlinear systems to be
!                               solved are more viable than was the case with
!                               Phoenix/InSight as a single system involving
!                               extremely small radii/high curvature at some
!                               aft body vertices.  New input ng_split helps
!                               control the two portions, with a default at
!                               the maximum radius point.
!                               Related inputs power_fore and power_aft are
!                               also new, to help control the curvature-based
!                               redistribution(s).
!     08/06/2013    "    "      All ADT variants have been merged into a module
!                               with generic build & search calls.
!     10/21/2013    "    "      A 90-degree skirt angle case required the skirt
!                               length (in x, not arc length) to be 0.  Now,
!                               this special case (only) interprets the skirt
!                               length to be in the y or r direction, and a
!                               skirt length > 0. improves the curvature-based
!                               redistribution around the shoulder because the
!                               the forebody ends beyond the skirt tangency pt.
!                               Also, if Rbase = Rnose is input for the Apollo
!                               case (half_cone_angle = 0.), then a semicircle/
!                               hemisphere is implied, and the maximum radius
!                               is only Rbase if Rshoulder = 0.  Rbase is
!                               updated accordingly in the printout.
!     10/24/2013    "    "      CURVDIS2 now iterates over lower powers until
!                               convergence is achieved (as well as normalizing
!                               and denormalizing the data), so the application
!                               program no longer needs to iterate on exponent.
!     11/01/2013    "    "      Automated clustering towards sharp vertices via
!                               ISMOOTH = -3 passed to modified CURVDIS[2].
!     11/04/2013    "    "      Chun Tang asked about growth rates in the
!                               grid spacing along the generatrix, so now they
!                               are tabulated as generatrix-growth-rates.dat.
!     11/05/2013    "    "      Blend non-analytic (redistributed input) gener-
!                               atrices aft of ng_split as for analytic cases.
!     12/04/2013    "    "      Avoid such blending by entering ng_split < 0,
!                               causing the full-body input generatrix to be
!                               redistributed as a single piece instead of two
!                               pieces.  This is harder, though, and may not be
!                               as good as two-part redistribution.
!     12/20/2013    "    "      Frustum_radius was not being converted from
!                               inches to meters.
!     12/31/2013    "    "      Confusion caused by entering dimensions in
!                               inches prompted printing of details in both
!                               meters and inches.
!     01/10/2013    "    "      Allow the specified rib thickness to increase a
!                               little to match the nearest grid spacing, thus
!                               avoiding some of the skewing that the option
!                               to round the ribs (actually flatten them) is
!                               then affected by.  Also, iterate once on the
!                               azimuthal redistribution ds to avoid more of
!                               the skewing.
!     02/07/2014    "    "      Still some cell skewing, so try iterating. It
!                               converges, but flattening the ridges the way
!                               we do still produces a slight kink along the
!                               rib where the desired thickness starts.  Also,
!                               the j = 2 spokewise grid line will contain
!                               some bend if the deflection is nonzero, so the
!                               only good way to resolve the rib semithickness
!                               with the j = 2 grid line and still get properly
!                               straight ridges is to leave them sharp,
!                               specially if SPLIT_GRID and UPDATE_GRID are
!                               going to be used to combine more than one
!                               deflection into a common surface grid.
!     05/29/2015    "    "      The Gridgen form of the desired generatrix
!                               was not being written for the nonanalytic case.
!     06/30/2016    "    "      A new option in CURVDIS to work with adjusted
!                               curvatures that are artificially broadened onto
!                               flat/low-curvature segments appears to be good
!                               for analytic geometries but unreliable for input
!                               generatrices of unknown quality/precision.
!                               New controls flat_blend_fore & flat_blend_aft
!                               allow this option to be suppressed.  It was
!                               prompted by a clear need to resolve heat flux
!                               spikes at the shoulder better by spreading the
!                               tight shoulder grid spacing forward onto the
!                               conical flank.  This should work for Apollo-
!                               type spherical sections too (large radius/low
!                               curvature) given that we work with normalized
!                               geometry during curvature-based redistribution.
!                               Namelist controls are now printed to standard
!                               output.
!     10/27/2016    "    "      Aft-end stretching towards the symmetry axis in
!                               the generatrix is undesirable.  3D surfaces
!                               replace this region with uniform quadrants.
!                               Therefore, do similar to the generatrix, using
!                               ni_regrid_aft_body and a heuristic number of
!                               uniform points at the end.
!     04/05/2018    "    "      An indexing error was found in the above uniform
!                               redistribution of the aft generatrix end.
!     05/11/2018    "    "      No change at this level, but some heuristics in
!                               detect_flatness and detect_vertices (geomlib)
!                               and verbose output to help interpret the
!                               curvature adjustments (see also handle_flatness)
!                               has been made permanent at the cost of bulky
!                               log files.  User control of the heuristics has
!                               not been provided (yet); change some parameter
!                               constants in the above if needed (sorry!).
!                               The original CURVDIS (gridlib) has also been
!                               translated to Fortran 90.  Example issue:
!                               A Mars Sample Return Lander turned out to have
!                               normalized curvatures on the spherical-section
!                               heat shield that oscillated between ~0.49 and
!                               ~0.51 (from a CAD generatrix). The "flatmax"
!                               constant in detect_flatness happened to be 0.5.
!                               The heat shield needs to be treated as for the
!                               zero-curvature aft-body conical sections, in
!                               order to force more grid points forward of the
!                               shoulder for better blending.  Similarly, the
!                               "spike" constant in detect_vertices needed to
!                               be lowered for the cusp at the shoulder/cone
!                               juncture to be treated as for a vertex.
!     05/12/2018    "    "      Added "verbose" to the namelist to allow
!                               suppression of lots of log file output from the
!                               curvature-based point distribution (an aid to
!                               understanding artificial broadening of spikes in
!                               the shape function at generatrix vertices and
!                               spreading onto flat segments).
!     08/23/2018    "    "      The analytic aft body case's forcing of uniform
!                               spacing away from the aft stag. point was not
!                               right because x/rgen(:) was being reused for the
!                               aft body (starting from 1, not numi_forebody).
!     09/25/2019    "    "      Rerunning on the capsule-on-a-sting analytic
!                               case gave troubling aft-body results, traced to
!                               interference of the vertex-detection scheme in
!                               CURVDIS with the (later) scheme for broadening
!                               the curvature-based shape function on to flat or
!                               low-curvature segments. Better blending on to
!                               such segments (and either side of sharp corners)
!                               is highly desirable, but fraught with potential
!                               difficulties due to the unavoidable heuristics.
!     10/30/2019    "    "      Introduced nblend_fore_aft control. The default
!                               of 8 can mean poor results going from tight
!                               spacing on the shoulder to larger spacing on
!                               the aft cone. Spreading the blending over more
!                               points should help.
!  Author:
!
!     David Saunders, ERC, Inc. at NASA Ames Research Center, Moffett Field, CA
!            Now with AMA, Inc. at NASA ARC.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Programmer Note:

!     All application-specific subroutines appear as internal procedures of the
!     short main program, from which they can inherit any variables desired.
!     This makes for minimal argument lists and eases adding new capability.

!  Modules:

      use adt_utilities            ! All ADT build & search routines
      use biconic_parameters       ! Derived data type for biconic forebody
      use grid_block_structure     ! Derived data type for one grid block
      use sphere_cone_parameters   ! Derived data type for sphere/cone forebody
      use sphere_cone              ! Analytic forebody generatrix utilities
      use surface_patch_utilities  ! 3-space surface grid utilities
      use xyzq_io_module           ! PLOT3D file I/O package

   implicit none

!  Constants:

   integer, parameter :: &
      lunstdin  = 5,     &         ! Namelist inputs on the command line
      luncrt    = 6,     &         ! Screen
      lunin1    = 1,     &         ! Input generatrix, if any
      lunin2    = 2,     &         ! Input normalized nose patch
!!!   lunchk    = 3,     &         ! For checking redistributed nose patch
      lunout    = 4,     &         ! Surface grid ready for 3D volume gridding
      lungen    = 7,     &         ! Output generatrix in Tecplottable form
      lungen0   = 8,     &         ! Generatrix before kappa-based redistribn.
      lungen1   = 9,     &         ! Output generatrix for 2D volume gridding
      lunspoke  = 10,    &         ! Interim spoked surface grid
      mxcones   = 12,    &         ! More than enough for any likely aft body
      nose_case = 2                ! Only active singular-pt. replacement option

   logical, parameter :: &
      false = .false.,   &
      true  = .true.

   real, parameter :: &
      half = 0.5, meters_per_inch = 0.0254, one = 1., zero = 0.

!  Variables:

   integer :: &
      i0_umbrella, i1_umbrella, i2_umbrella, ib, ier, imatch, ios, &
      ismooth_aft, ismooth_fore, lunprint, &
      nb, nblend_fore_aft, nblocks, ncones, nedges, ng_split, ni_outer, &
      ni_regrid, ni_regrid_aft_body, ni_regrid_forebody, ni_template, &
      ntoroids, ni_spoke, nj_quarter, nj_semi, nj_spoke, nj_template, &
      numi, numi_aft_body, numi_forebody, numi_generatrix, numi_internal, &
      numi_total, numf, numj, nvertices, rounding_mode, ivertex(mxcones), &
      ni_round(mxcones)

   real :: &
      chords(mxcones), cone_angle_aft(mxcones), &
      ds_forebody, ds_round(mxcones), frustum_radius, geometric_scale, &
      half_cone_angle, half_cone_angle_aft, half_cone_angle_fore, &
      peak_deflection, peak_ripple, power, power_aft, power_fore, &
      r_aft(mxcones), r_nose, radius_base, radius_cone_juncture, radius_nose, &
      radius_shoulder, radius_vertex, rib_deflection, rib_thickness, &
      semiangle, skirt_angle, skirt_length, sting_stretch_multiplier, &
      x_aft(mxcones), x_nose

   logical :: &
      aft_body_too, analytic, cell_centered, collapsed, &
      flat_blend_aft, flat_blend_fore, &
      resolve_the_ridges, ripple_case, round_the_ridges, spherecone, &
      sting_case, umbrella_case, units_in_inches, use_input_generatrix, &
      verbose

   character :: &
      input_generatrix*80, nose_patch_file_name*80, output_generatrix*80, &
      surface_grid_file_name*80

   real, allocatable, dimension (:) :: &
      xfore, rfore,                    &  ! For the interim discretization
      xgen, rgen, sgen,                &  ! For the output generatrix
      xgen_new, rgen_new, sgen_new        ! Moved here for non-analytic case

   real, allocatable, dimension (:) :: &  ! For changing pt. counts along lines
      snorm, x1, y1, z1, x2, y2, z2

   type (biconic_type) :: &
      biconic                   ! For biconic forebody parameters

   type (sphere_cone_type) :: &
      capsule                   ! For sphere/cone forebody parameters

   type (grid_type) :: &
      spoke_grid,      &        ! Single-block, uniform azimuth, singular apex
      spoke_fore,      &        ! Copies that allow same work-space fore & aft
      spoke_aft                 ! for replacing the singular point region

   type (grid_type) :: &
      new_fore(6),     &        ! Copies that allow same work-space fore & aft
      new_aft(6)                ! for replacing the singular point region

   type (grid_type), pointer, dimension (:) :: &
      new_grid, &               ! Desired regridded surface, 6 or 12 patches
      template                  ! Normalized quarter circle grid, with x = 0.

   namelist /FOREBODY_INPUTS/ &
      verbose, aft_body_too, units_in_inches, geometric_scale, &
      input_generatrix, spherecone, numi_forebody, numj, x_nose, r_nose, &
      radius_nose, radius_base, radius_shoulder, radius_vertex, &
      half_cone_angle, half_cone_angle_fore, half_cone_angle_aft, &
      radius_cone_juncture, skirt_angle, skirt_length, nose_patch_file_name, &
      output_generatrix, surface_grid_file_name, flat_blend_fore, power_fore, &
      ni_regrid_forebody, ripple_case, ntoroids, peak_ripple, umbrella_case, &
      nedges, peak_deflection, rib_deflection, frustum_radius, &
      resolve_the_ridges, rib_thickness, round_the_ridges

   namelist /AFT_BODY_INPUTS/ &
      sting_case, numi_aft_body, ncones, flat_blend_aft, power_aft, ng_split, &
      nblend_fore_aft, ni_regrid_aft_body, cone_angle_aft, r_aft, &
      rounding_mode, ds_round, ni_round, sting_stretch_multiplier, &
      i0_umbrella, i1_umbrella, i2_umbrella

!  Execution:
!  ::::::::::

   call control ()               ! Read controls and open files

   call set_up_generatrix ()     ! Read or construct a generatrix;
                                 ! derive a body of revolution spoked grid

   call morph_to_umbrella ()     ! Option to model flexible TPS

   call adjust_template ()       ! Match circumf. point count to initial grid

   call initial_regrid (true)    ! First stage of removing the singular point

   call project_to_nose ()       ! Morphed nose patches are projected precisely

   if (aft_body_too .and. .not. sting_case) then  ! Remove aft singular point

      call initial_regrid (false)

      call project_to_nose ()

   end if

   call save_result ()           ! And clean up

   contains

!     Internal procedures for program capsule_grid:
!     :::::::::::::::::::::::::::::::::::::::::::::

!     --------------------------------------------------------------------------
      subroutine control ()     ! Read controls; set up the output files
!     --------------------------------------------------------------------------

!     Local variables:

      real :: scale_factor

!     Execution:

!     Default the namelist inputs:

!     $FOREBODY_INPUTS
      verbose = false              ! T sheds light on some CURVDIS heuristics
      aft_body_too = false         ! Option to suppress the aft body
      units_in_inches = false      ! Option suited to arc-jets in the US of A
      geometric_scale = one        ! Option to scale the indicated geometry
      input_generatrix = 'none'    ! Full generatrix file; 'none' => analytic
      spherecone = true            ! T => half_cone_angle >= 0, else biconic
      numi_forebody = 201          ! # forebody grid points
      numj = 121                   ! # spokes in body of revolution (4m + 1)
      x_nose = zero                ! Nose coordinate
      r_nose = zero                ! Nose coordinate; must be 0. if aft body on
      radius_nose = 0.0889         ! Radius of nose circular segment
      radius_base = 0.1778         ! Radius at maximum diameter
      radius_shoulder = 0.00889    ! Radius of shoulder circular segment
      radius_vertex = zero         ! Biconic option to round the cone juncture
      half_cone_angle = 55.        ! 0 => spherical section (Apollo-like)
      half_cone_angle_fore = 70.   ! Forward cone angle, if biconic
      half_cone_angle_aft = 55.    ! Aft cone angle, if biconic
      radius_cone_juncture = 0.1   ! Radius where cones meet, if biconic
      skirt_angle = 35.            ! Semi-angle if cone segment follows shoulder
      skirt_length = 0.001         ! An x length, >=< 0. unless skirt_angle = 90
      nose_patch_file_name = 'NosePatch61B.p3da'
                                   ! Grid quadrant used to remove singular point
      output_generatrix = 'gen.dat'! Generatrix corresp. to the surface grid
      surface_grid_file_name = 'surface_grid.g'
                                   ! Output surface grid/no singular point(s)
      flat_blend_fore = true       ! Puts more points on flat/low curv. segments
      power_fore = 0.5             ! Initial exponent for curvature-based redis.
      ni_regrid_forebody = 33      ! i = 1:ni_regrid is regridded for all j
      ripple_case = false          ! Option to model inflatable toroid effects
      ntoroids = 0                 ! # underlying toroids along cone
      peak_ripple = 0.01           ! Peak deflection between toroids
      umbrella_case = false        ! T => polygonal forebody/optional deflectns.
      nedges = 8                   ! # full polygon edges if umbrella_case
      peak_deflection = 0.01       ! > 0. => catenary-type fabric deflections
      rib_deflection = 0.          ! Option to bend ribs (+/- => concave/convex)
      frustum_radius = 0.          ! 0. => start fabric/rib deflns. @ tngcy. pt.
      resolve_the_ridges = true    ! T => adjust spacing either side of ribs
      rib_thickness = 0.01         ! Off-rib spacing to impose (<= input value)
      round_the_ridges = true      ! T => tight local cubic spline rounding
!     $END

      read (lunstdin, nml=forebody_inputs)
      write (luncrt,  nml=forebody_inputs)

      lunprint      = -luncrt
      if (verbose) lunprint = luncrt
      analytic      = trim (input_generatrix) == 'none'
      ismooth_aft   = 3  ! Both types of CURVDIS smoothings; no added blending
      if (flat_blend_fore) ismooth_fore = -3  ! Allow added blending
      numi_total    = numi_forebody   ! Following curvature-based redistribution
      numi_internal = numi_forebody*2 - 1  ! To resolve shoulder adequately

!     Some aft body defaults help catch missing inputs:

!     $AFT_BODY_INPUTS
      sting_case = false           ! Option to suppress closure of aft body
      ncones = 2                   ! ncones >= 1
      r_aft(:) = -one              ! Radii from fore/aft juncture (1:ncones+1)
      cone_angle_aft(:) = -999.    ! Half cone half angles (1:ncones); 0, 90 OK
      rounding_mode = 2            ! 1 => tangent length; 2 => radius; 0 = none
      ds_round(:) = zero           ! Tangent lengths or rounding radii >=< 0.
      ni_round(:) = 0              ! # uniform rounding circle arc points >= 0
      numi_aft_body = 201          ! # grid points on aft body
      sting_stretch_multiplier = 2.! Applied to interim end-of-sting spacing
      flat_blend_aft = true        ! Puts more points on flat/low curv. segments
      power_aft = 0.5              ! Initial exponent for curvature-based redis.
      ng_split = 0                 ! Full-body input generatrix split index
      nblend_fore_aft = 8          ! # aft pts. used to match fore/aft spacing
      ni_regrid_aft_body = 33      ! Aft i = 1:ni_regrid is regridded for all j;
      i0_umbrella = 0              ! May be needed for fat ADEPT aft bodies
      i1_umbrella = 0              ! Inboard aft spoke index for morphing
      i2_umbrella = 0              ! Outboard index; 0s suppress aft morphing
!     $END                         ! ni_regrid is overridden for analytic cases

      if (aft_body_too) then
         read (lunstdin, nml=aft_body_inputs)
         write (luncrt,  nml=aft_body_inputs)
         ismooth_aft = 3  ! Both types of CURVDIs smoothings; no added blending
         if (flat_blend_aft) ismooth_aft = -3  ! Allow added blending
         numi_total = numi_total + numi_aft_body
         numi_generatrix = numi_total - 1      ! Avoid common point
      else
         sting_case = false
         numi_generatrix = numi_forebody
      end if

      scale_factor = geometric_scale
      if (units_in_inches) scale_factor = scale_factor * meters_per_inch

      x_nose               = scale_factor * x_nose
      r_nose               = scale_factor * r_nose
      radius_nose          = scale_factor * radius_nose
      radius_base          = scale_factor * radius_base
      radius_shoulder      = scale_factor * radius_shoulder
      radius_vertex        = scale_factor * radius_vertex
      radius_cone_juncture = scale_factor * radius_cone_juncture
      skirt_length         = scale_factor * skirt_length
      peak_ripple          = scale_factor * peak_ripple
      peak_deflection      = scale_factor * peak_deflection
      rib_deflection       = scale_factor * rib_deflection
      rib_thickness        = scale_factor * rib_thickness
      frustum_radius       = scale_factor * frustum_radius
      r_aft(:)             = scale_factor * r_aft(:)
      ds_round(:)          = scale_factor * ds_round(:)

      ni_spoke   = numi_forebody   ! To suit code in umbrella forebody option
      nj_spoke   = numj            ! To suit code from FOREBODY_REGRID
      nj_quarter = (numj + 1) / 2  ! # spokes on 1/4th of (whole) forebody

      if (mod (numj - 1, 4) /= 0) then
         write (luncrt, '(/, a, i4)') &
            'The input number of spokes must be of the form 4m + 1: ', numj
         go to 99  ! Single stop philosophy
      end if

!     The nose patch is the first quadrant of the (y, z) plane at x = 0.

      open (lunin2, file=nose_patch_file_name, status='old', iostat=ios)

      if (ios /= 0) then
         write (luncrt, '(/, 2a)') &
            'Trouble opening nose patch file: ', trim (nose_patch_file_name)
         go to 99
      end if

      call xyzq_read (lunin2, -1, true, nb, numf, cell_centered, template, ios)

      if (ios /= 0) then
         write (luncrt, '(/, a)') 'Trouble reading normalized nose patch.'
         go to 99
      end if

      ni_template = template(1)%ni
      nj_template = template(1)%nj

      if (ni_template /= nj_template) then
         write (luncrt, '(/, a, 2i5)') &
           'Input nose patch dimensions don''t match:', ni_template, nj_template
         go to 99
      end if

      open (lunout, file=surface_grid_file_name, status='unknown')

      return

 99   stop

      end subroutine control

!     --------------------------------------------------------------------------
      subroutine set_up_generatrix ()  ! Read or construct the generatrix
!
!     If an aft body is present, both portions of the generatrix are set up
!     before any 3-space surface gridding is done.  If a generatrix has been
!     input, it is used precisely unless the specified number of grid points
!     differs from the number of generatrix points read.  If the point counts
!     differ, curvature-based redistribution is performed as for analytic cases.
!     --------------------------------------------------------------------------

!     Local constants:

      real, parameter :: ninety = 90., iperm = one / meters_per_inch

!     Local variables:

      integer :: i, j, nisofar, nl
      real    :: angle, dl, dr, du, dx, L, segment, total

!     Execution:

      allocate (xgen(numi_total), rgen(numi_total))

      if (.not. analytic) then

         call non_analytic ()  ! Read a generatrix and impose numi_total points

      else  ! Discretize the specified forebody

         use_input_generatrix = false

!        Outputs are assumed to be in meters as expected by flow solver DPLR.
!        Salient details are printed by sphere|biconic_cone_control_points.
!        If inputs are initially in inches, it should help to avoid confusion
!        if the same output details are shown in inches as well (see below).

         write (luncrt, '(/, a)') 'Analytic defining points (meters):'

         allocate (xfore(numi_internal), rfore(numi_internal))

         if (spherecone) then  ! Could be Apollo-type too

            capsule%x_nose               = x_nose
            capsule%r_nose               = r_nose
            capsule%radius_nose          = radius_nose
            capsule%radius_base          = radius_base
            capsule%radius_shoulder      = radius_shoulder
            capsule%half_cone_angle      = half_cone_angle
            capsule%skirt_angle          = skirt_angle
            capsule%skirt_length         = skirt_length

            call sphere_cone_control_points (capsule)

            if (units_in_inches) then  ! See capsules.f90
               write (luncrt, '(/, a)') 'Analytic defining points (inches):'
               write (*, '(/, (a, 3x, 2f15.10))') &
                  'Nose                   x and r:', &
                     x_nose*iperm, r_nose*iperm, &
                  'Rnose, Rbase:                  ', &
                     radius_nose*iperm, radius_base*iperm, &
                  'Rshoulder, semi-cone angle:    ', &
                     radius_shoulder*iperm, half_cone_angle, &
                  'Skirt angle, skirt length:     ', &
                     skirt_angle, skirt_length*iperm
               if (half_cone_angle > zero) then  ! Not Apollo-type
                  L = sqrt ((capsule%xt_shoulder - capsule%xt_nose)**2 + &
                            (capsule%rt_shoulder - capsule%rt_nose)**2)
                  write (*, '(a, 3x, 2f15.10)') &
                     'Nose tangency          x and r:', &
                        capsule%xt_nose*iperm, capsule%rt_nose*iperm, &
                     'Conical flank arc length:      ', &
                        L*iperm
               end if
               write (*, '(a, 3x, 2f15.10)') &
                  'Shoulder tangency      x and r:',    &
                     capsule%xt_shoulder*iperm, capsule%rt_shoulder*iperm, &
                  'Skirt tangency         x and r:',    &
                     capsule%xt_skirt*iperm, capsule%rt_skirt*iperm, &
                  'Maximum-diameter       x and r:',    &
                     capsule%x_base*iperm, capsule%r_base*iperm, &
                  'Skirt cut-off          x and r:',    &
                     capsule%x_skirt*iperm, capsule%r_skirt*iperm, &
                  'Nose circle center     x and r:',    &
                     capsule%xcenter_nose*iperm, capsule%r_nose*iperm, &
                  'Shoulder circle center x and r:',    &
                     capsule%x_base*iperm, capsule%rcenter_shoulder*iperm, &
                  ' '
            end if

            call sphere_cone_discretization (capsule, numi_internal, &
                                             xfore, rfore)
            if (ripple_case) call toroid_simulation ()  ! Perturbs generatrix

            if (umbrella_case) call bend_rib ()  ! Option to bend cone portion

         else  ! Biconic option

            biconic%x_nose               = x_nose
            biconic%r_nose               = r_nose
            biconic%radius_nose          = radius_nose
            biconic%radius_cone_juncture = radius_cone_juncture
            biconic%radius_base          = radius_base
            biconic%radius_shoulder      = radius_shoulder
            biconic%radius_vertex        = radius_vertex
            biconic%half_cone_angle_fore = half_cone_angle_fore
            biconic%half_cone_angle_aft  = half_cone_angle_aft
            biconic%skirt_angle          = skirt_angle
            biconic%skirt_length         = skirt_length

            call biconic_control_points (biconic)

            if (units_in_inches) then  ! See capsules.f90
               write (luncrt, '(/, a)') 'Analytic defining points (inches):'
               write (*, '(/, (a, 3x, 2f15.10))') &
                  'Nose                   x and r:',    &
                     x_nose*iperm, r_nose*iperm, &
                  'Rnose, Rbase:                  ',    &
                     radius_nose*iperm, radius_base*iperm, &
                  'Rshoulder, Rcone_juncture:     ',    &
                     biconic%radius_shoulder*iperm,     &
                     biconic%radius_cone_juncture*iperm, &
                  'Fore and aft cone semi-angles: ',    &
                     half_cone_angle_fore, half_cone_angle_aft, &
                  'Fore and aft frustum lengths:  ',    &
                     biconic%cone_length_fore*iperm,    &
                     biconic%cone_length_aft*iperm,     &
                  'Skirt angle, skirt length:     ',    &
                     skirt_angle, skirt_length*iperm,   &
                  'Nose tangency          x and r:',    &
                     biconic%xt_nose*iperm, biconic%rt_nose*iperm, &
                  'Cone juncture          x and r:',    &
                     biconic%x_cone_juncture*iperm,     &
                     biconic%r_cone_juncture*iperm,     &
                  'Shoulder tangency      x and r:',    &
                     biconic%xt_shoulder*iperm, biconic%rt_shoulder*iperm, &
                  'Skirt tangency         x and r:',    &
                     biconic%xt_skirt*iperm, biconic%rt_skirt*iperm, &
                  'Maximum-diameter       x and r:',    &
                     biconic%x_base*iperm, biconic%r_base*iperm, &
                  'Skirt cut-off          x and r:',    &
                     biconic%x_skirt*iperm, biconic%r_skirt*iperm, &
                  'Nose circle center     x and r:',    &
                     biconic%xcenter_nose*iperm, biconic%r_nose*iperm, &
                  'Shoulder circle center x and r:',    &
                     biconic%x_base*iperm, biconic%rcenter_shoulder*iperm
               if (radius_vertex > zero) then
                  write (*, '(a, 3x, f15.10, /, (a, 3x, 2f15.10))') &
                     'Rvertex:                       ', radius_vertex*iperm, &
                     'Vertex circle center   x and r:', &
                         biconic%xcenter_vertex*iperm,  &
                         biconic%rcenter_vertex*iperm,  &
                     'Vertex tangency (fore) x and r:', &
                         biconic%xt_vertex_fore*iperm,  &
                         biconic%rt_vertex_fore*iperm,  &
                     'Vertex tangency (aft)  x and r:', &
                         biconic%xt_vertex_aft*iperm,   &
                         biconic%rt_vertex_aft*iperm
               end if
               write (luncrt, '(a)')
            end if

            call biconic_discretization (biconic, numi_internal, xfore, rfore)

         end if

!        In case the curvature-based redistribution has trouble, save its input:

         open  (lungen0, file='raw_generatrix_fore.dat', status='unknown')
         write (lungen0, '(2es20.12)') (xfore(i),rfore(i), i = 1, numi_internal)
         close (lungen0)

!        Impose curvature-based spacing for the requested # forebody points:

         power = power_fore  ! CURVDIS2 now tries lower powers until success

         call curvdis2 (numi_internal, xfore, rfore, numi_forebody, power, &
                        ismooth_fore, lunprint, false, xgen, rgen, ier)
         if (ier /= 0) then
            write (luncrt, '(a, i3, f5.1, /, a)') &
               'CURVDIS2 trouble (forebody).  IER, POWER:', ier, power, &
               'Aborting.'
            stop
         end if

         rgen(1) = r_nose  ! Exactly

!        Save the analytic forebody generatrix for possible use elsewhere.
!        If there's an aft body, its generatrix will be appended to this one.

         open  (lungen, file=output_generatrix, status='unknown')
         write (lungen, '(a1, i4)') '#', numi_generatrix
         write (lungen, '(2es15.7, f3.0)') &
            (xgen(i), rgen(i), zero, i = 1, numi_forebody)

!        That is Tecplotable.  Save the form expected by Gridgen scripts too.

         if (aft_body_too) then
            open (lungen1, file='full-body-generatrix.dat', status='unknown')
         else
            open (lungen1, file='forebody-generatrix.dat',  status='unknown')
         end if
         write (lungen1, '(i4)') numi_generatrix
         write (lungen1, '(2es15.7, f3.0)') &
            (xgen(i), rgen(i), zero, i = 1, numi_forebody)

      end if

!     In case slopes need to be matched at the forebody/aft body juncture:

      j = numi_forebody
      angle = atan2d (rgen(j-1) - rgen(j), xgen(j) - xgen(j-1))
      write (luncrt, '(/, a, f11.6, a)') &
         'Effective angle at forebody/aft body juncture:', angle, ' degrees'

!     Generate the interim spoked grid as a body of revolution, single block:

      numi = numi_forebody

      call spoked_surface (1)

!     Analytic aft body?  If so, we complete the generatrix before dealing with
!     the forebody "umbrella" option and/or regridding to avoid the singular pt.

      if (analytic .and. aft_body_too) then

         if (r_aft(ncones+1) == -one) then  ! Common user error to omit last r
            write (luncrt, '(/, a)') &
               '*** You must enter ncones radii > 0 plus r_aft(ncones+1) >= 0.'
            stop
         end if

         ds_forebody = distance (xgen(j-1), rgen(j-1), xgen(j), rgen(j))

!        Begin the analytic aft body construction:

         numi     = numi_aft_body
         x_aft(1) = xgen(j)
         r_aft(1) = rgen(j)
         xgen(1)  = x_aft(1)  ! Note reuse of forebody computational grid arrays
         rgen(1)  = r_aft(1)

!        Calculate the abscissas of the other vertices, fore to aft:

         do j = 1, ncones
            if (cone_angle_aft(j) == ninety) then
               x_aft(j+1) = x_aft(j)
            else if (cone_angle_aft(j) == zero) then  ! Kludge: use r_aft(j+1)
               x_aft(j+1) = x_aft(j) + r_aft(j+1)     ! for the length of the
               r_aft(j+1) = r_aft(j)                  ! implied cylindrical seg.
            else if (cone_angle_aft(j) == 180.) then  ! Similar kludge:
               x_aft(j+1) = x_aft(j) + r_aft(j+1)     ! the step length is -ve
               r_aft(j+1) = r_aft(j)
            else
               x_aft(j+1) = x_aft(j) + (r_aft(j) - r_aft(j+1)) / &
                                       tand (cone_angle_aft(j))
            end if
         end do

         collapsed = false  ! No rounding if last segment is vertical

         if (r_aft(ncones+1) > zero) then ! Close with implied vertical segment?
            if (sting_case) then
               nvertices = ncones + 1
            else
               r_aft(ncones+2) = zero     ! if cone_angle_aft(ncones) = 0., this
               x_aft(ncones+2) = x_aft(ncones+1) ! correctly implies no collapse
               nvertices = ncones + 2     ! because r_aft(ncones+1) = length > 0
            end if
         else  ! Last input radius is explicitly zero
            nvertices = ncones + 1
            if (cone_angle_aft(ncones) < ninety) collapsed = true
         end if

!        Display the vertices:

         if (units_in_inches) then
            write (luncrt, '(/, a, /, 2(17x, a1))') &
               'Aft body vertices (inches):', 'x', 'r'
            write (luncrt, '(2es18.7)') &
               (x_aft(j)/meters_per_inch, r_aft(j)/meters_per_inch, &
               j = 1, nvertices)
         end if

         write (luncrt, '(/, a, /, 2(17x, a1))') &
            'Aft body vertices (meters):', 'x', 'r'
         write (luncrt, '(2es18.7)') (x_aft(j), r_aft(j), j = 1, nvertices)
         write (luncrt, '(/, (a))') &
            'Initial discretization of cone segments:', &
            '   i1   nl     arc length     uniform ds'

         open  (lungen0, file='vertices.dat', status='unknown')
         write (lungen0, '(2es18.7)') (x_aft(j), r_aft(j), j = 1, nvertices)
         close (lungen0)

!        Discretize the specified conical frustum(s) by linear interpolation:

         call chords2d (nvertices, x_aft, r_aft, false, total, chords)

         du = total / real (numi - 1)
         nisofar = 1

!        Make each linear segment spacing as near as possible to this du while
!        capturing the vertices precisely:

         do j = 1, nvertices - 1
            segment = chords(j+1) - chords(j)
            nl = nint (segment / du)
            if (j == nvertices - 1) nl = numi - nisofar  ! Any slack in last seg
            if (j == 1) imatch = nl ! In case of no rounding; see match_forebody
            dx = (x_aft(j+1) - x_aft(j)) / real (nl)
            dr = (r_aft(j+1) - r_aft(j)) / real (nl)
            dl = sqrt (dx**2 + dr**2)
            ivertex(j) = nisofar
            write (luncrt, '(2i5, 2es15.7)') nisofar, nl, segment, dl
            if (j == nvertices - 1) write (luncrt, '(a)')
            do i = nisofar + 1, nisofar + nl
               xgen(i) = x_aft(j) + dx*real (i-nisofar)
               rgen(i) = r_aft(j) + dr*real (i-nisofar)
            end do
            nisofar = nisofar + nl
         end do

         ivertex(nvertices) = numi
         ni_regrid  = nl + 1            ! # points on last segment (-> 1st seg)
         xgen(numi) = x_aft(nvertices)  ! Precisely
         rgen(numi) = r_aft(nvertices)

!        Option to round the vertices then redistribute according to curvature:

         call round_vertices ()

!        Match spacing at the forebody/aftbody juncture:

         call match_forebody ()

!        Avoid stretching towards the symmetry axis at the aft end:

!!!      write (6, '(a, i6)') 'numi_generatrix:', numi_generatrix
!!!      write (6, '(i, 2es15.7)') (i, xgen(i), rgen(i), i = 1, numi)

         call uniform_aft_end (numi)

!        Append the aft-body generatrix to the forebody in the output files,
!        omitting the common point:

         rgen(numi) = r_aft(nvertices)  ! Exactly
         write (lungen, '(2es15.7, f3.0)') (xgen(i), rgen(i), zero, i = 2, numi)
         write (lungen1,'(2es15.7, f3.0)') (xgen(i), rgen(i), zero, i = 2, numi)

      end if  ! Analytic aft body choice

!     If an arbitrary full-body generatrix was read, arrange for the aft portion
!     to look the same as the analytic case from here on:

      if (.not. analytic .and. aft_body_too) then
         xgen(1:numi_aft_body) = xgen(numi_forebody:numi_generatrix)
         rgen(1:numi_aft_body) = rgen(numi_forebody:numi_generatrix)
         ni_regrid = ni_regrid_aft_body
      end if

      if (aft_body_too) then

         if (sting_case) then  ! No regridding of aft body

         else  ! Make i = 1 at the aftmost singular pt. to be removed:

            call rverse (numi_aft_body, xgen, xgen)
            call rverse (numi_aft_body, rgen, rgen)

         end if

!        Generate the interim spoked grid as a body of revolution, single block.
!        We keep copies around so that the nose/tail regridding can be done with
!        the same work-space and code.

         allocate (spoke_fore%x(numi_forebody,numj,1), &
                   spoke_fore%y(numi_forebody,numj,1), &
                   spoke_fore%z(numi_forebody,numj,1))

         spoke_fore%x(:,:,:) = spoke_grid%x(:,:,:)  ! From the above call to
         spoke_fore%y(:,:,:) = spoke_grid%y(:,:,:)  ! spoked_surface for the
         spoke_fore%z(:,:,:) = spoke_grid%z(:,:,:)  ! forebody

         deallocate (spoke_grid%x, spoke_grid%y, spoke_grid%z)

         numi = numi_aft_body

         call spoked_surface (2)  ! Single-patch interim aft body grid

         allocate (spoke_aft%x(numi_aft_body,numj,1), &
                   spoke_aft%y(numi_aft_body,numj,1), &
                   spoke_aft%z(numi_aft_body,numj,1))

         spoke_aft%x(:,:,:) = spoke_grid%x(:,:,:)
         spoke_aft%y(:,:,:) = spoke_grid%y(:,:,:)
         spoke_aft%z(:,:,:) = spoke_grid%z(:,:,:)

!        Set up for regridding the forebody as for the no-aft-body case:

         deallocate (spoke_grid%x, spoke_grid%y, spoke_grid%z)

         allocate (spoke_grid%x(numi_forebody,numj,1), &
                   spoke_grid%y(numi_forebody,numj,1), &
                   spoke_grid%z(numi_forebody,numj,1))

         spoke_grid%ni = numi_forebody  ! In case of the umbrella step
         spoke_grid%nj = numj
         spoke_grid%nk = 1
         spoke_grid%x(:,:,:) = spoke_fore%x(:,:,:)
         spoke_grid%y(:,:,:) = spoke_fore%y(:,:,:)
         spoke_grid%z(:,:,:) = spoke_fore%z(:,:,:)

         if (.not. umbrella_case) then
            deallocate (spoke_fore%x, spoke_fore%y, spoke_fore%z)
         else
            ! We replace the axisymmetric spokes with the faceted spokes later.
         end if  ! See the save_results procedure.

         ni_regrid_aft_body = ni_regrid  ! But any rounding can change it

      end if  ! Aft body choice

      if (analytic) then
         close (lungen)
         close (lungen1)
      end if

      deallocate (xgen, rgen)

!     Set up to regrid the nose of the forebody:

      numi      = numi_forebody
      ni_regrid = ni_regrid_forebody

      end subroutine set_up_generatrix

!     -----------------------------------------------------------------------
      subroutine non_analytic ()  ! Read a generatrix and adjust its # points
!     -----------------------------------------------------------------------

      character, parameter :: mono *1 = 'M'  ! Monotonic for proportional adjust
      character, parameter :: loose*1 = 'B'  ! Loose for blending fore & aft

      integer :: i, i1, i2, ios, ng, ng_aft, ng_fore, ni
      real    :: ds2, total
      real, allocatable, dimension (:) :: chords, snew, vnew, vnewp, xg, rg

!     Execution:

      open (lungen, file=input_generatrix, status='old', iostat=ios)

      if (ios /= 0) then
         write (luncrt, '(/, 2a)') 'Trouble opening input generatrix ', &
            trim (input_generatrix)
         stop
      end if

      ng = 0
      do  ! Until EOF
         read (lungen, *, iostat=ios)
         if (ios < 0) exit
         ng = ng + 1
      end do

      allocate (xg(ng), rg(ng))

      rewind (lungen)
      do i = 1, ng  ! Ignore any 3rd column
         read (lungen, *) xg(i), rg(i)
      end do
      close  (lungen)

      use_input_generatrix = ng == numi_generatrix

      if (use_input_generatrix) then  ! Don't touch the input definition

         xgen(1:ng) = xg(:)
         rgen(1:ng) = rg(:)

      else  ! We use the input generatrix but redistribute it two ways

!        In case it proves helpful, adjust the number of points proportionally
!        (3-space utility) and save the resulting generatrix (only):

         allocate (z1(ng), z2(numi_generatrix))

         z1(:) = zero

         call adjustn (ng, 1, ng, xg, rg, z1, 1, numi_generatrix, xgen, rgen, &
                       z2, mono)
         deallocate (z1, z2)

!        Save this form of the redistributed generatrix:

         open  (lungen, file='proportionally_adjusted_gen.dat',status='unknown')
         write (lungen, '(a1, i4)') '#', numi_generatrix
         write (lungen, '(2es15.7, f3.0)') &
            (xgen(i), rgen(i), zero, i = 1, numi_generatrix)
         close (lungen)

!        Impose curvature-based spacing on the input generatrix and proceed.
!        The full-body case may be done as a single curve, but very tight
!        aft-body spacing/small radii can be challenging.  Doing it in two
!        parts may help some cases.

         power = power_fore      ! May apply to the single-piece full body case

         if (aft_body_too) then  ! We normally split the input generatrix
                                 ! somewhere near/aft of the max. diameter
            ng_fore = ng_split

            if (ng_fore == 0) then  ! Place it at the maximum radius index
               do i = 2, ng
                  if (rg(i) < rg(i-1)) then
                     ng_fore = i - 1
                     exit
                  end if
               end do
               write (luncrt, '(a, i4)') &
                  'Splitting generatrix at index', ng_fore

            else if (ng_split < 0) then  ! Option to redistribute as one piece

               call curvdis2 (ng, xg, rg, numi_generatrix, power, &
                              ismooth_fore, lunprint, false, xgen, rgen, ier)
               if (ier /= 0) then
                  write (luncrt, '(a, i3, f5.1, /, a)') &
               'CURVDIS2 trouble (input generatrix, full-body).  IER, POWER:', &
                     ier, power, 'Aborting.'
                  stop
               end if

               go to 95  ! Rare case where go to is best not avoided
            end if
         else
            ng_fore = ng  ! No aft body
         end if

!        Forebody gridding:

         call curvdis2 (ng_fore, xg, rg, numi_forebody, power, &
                        ismooth_fore, lunprint, false, xgen, rgen, ier)
         if (ier /= 0) then
            write (luncrt, '(a, i3, f5.1, /, a)') &
               'CURVDIS2 trouble (input generatrix, forebody).  IER, POWER:', &
               ier, power, 'Aborting.'
            stop
         end if

         if (aft_body_too) then  ! Allow for tighter radii on the aft body

            ng_aft = ng - ng_fore + 1
            power = power_aft

            call curvdis2 (ng_aft, xg(ng_fore), rg(ng_fore), numi_aft_body, &
                           power, ismooth_aft, lunprint, false, &
                           xgen(numi_forebody), rgen(numi_forebody), ier)
            if (ier /= 0) then
               write (luncrt, '(a, i3, f5.1, /, a)') &
                'CURVDIS2 trouble (input generatrix, aft body).  IER, POWER:', &
                 ier, power, 'Aborting.'
               stop
            end if

!           Should we try to blend fore and aft?  The user can choose
!           numi_forebody and numi_aft_body appropriately to help.  See also
!           ng_split and power_aft.  But blend anyway over a moderate # points,
!           along the lines of match_forebody elsewhere for analytic cases:

            ni = numi_forebody  ! End of forebody = ...
            i1 = ni             ! ... start of aft body

            ds_forebody = sqrt ((xgen(ni) - xgen(ni-1))**2 + &
                                (rgen(ni) - rgen(ni-1))**2)

!!!         ni = 8   ! Anything reasonable, assuming ng_split is sensible;
!!!                  ! it's hard to make it scale safely with numi_aft_body
!           The choice of 8 has been found to be poor when tight spacing on the
!           shoulder is being blended with larger spacing on the aft cone.
!           Spacing can actually be largest in the middle of the blend interval.
!           Spreading the blending over a wider range seems to be all we can do,
!           so now it's a user input.

            ni = nblend_fore_aft
            i2 = i1 + ni - 1

            allocate (chords(ni), snew(ni), vnew(ni), vnewp(ni))

            ds2 = sqrt ((xgen(i2) - xgen(i2-1))**2 + (rgen(i2) - rgen(i2-1))**2)

            call chords2d (ni, xgen(i1), rgen(i1), false, total, chords)

!           Safeguard abnormal cases:

!!          write (luncrt, '(a, 3es15.7)') 'ds_fore, ds2, total:', &
!!                                          ds_forebody, ds2, total
            if (ds_forebody > (total - ds2)*half) go to 90

            snew(1) = zero;  snew(ni) = total

            call vinokur (1, ni, ds_forebody, ds2, snew, luncrt, ier)

            call lcsfit (ni, chords, xgen(i1), true, loose, ni, snew, vnew, &
                         vnewp)
            xgen(i1+1:i2-1) = vnew(2:ni-1)

            call lcsfit (ni, chords, rgen(i1), true, loose, ni, snew, vnew, &
                         vnewp)
            rgen(i1+1:i2-1) = vnew(2:ni-1)

 90         deallocate (chords, snew, vnew, vnewp)

         end if

 95      continue

!        Avoid stretching towards the symmetry axis at the aft end:

!!!      write (6, '(a, i6)') 'At 95, numi_generatrix:', numi_generatrix
!!!      write (6, '(i, 2es15.7)') (i, xgen(i), rgen(i), i = 1, numi_generatrix)

         if (aft_body_too) call uniform_aft_end (numi_generatrix)

!        Save the curvature-based generatrix in Tecplot form:

         open  (lungen, file=output_generatrix, status='unknown')
         write (lungen, '(a1, i5)') '#', numi_generatrix
         write (lungen, '(2es15.7, f3.0)') &
            (xgen(i), rgen(i), zero, i = 1, numi_generatrix)
         close (lungen)

      end if

!     Save the generatrix in Gridgen form (# points at the top):

      if (aft_body_too) then
         open (lungen1, file='full-body-generatrix.dat', status='unknown')
      else
         open (lungen1, file='forebody-generatrix.dat',  status='unknown')
      end if
      write (lungen1, '(i4)') numi_generatrix
      write (lungen1, '(2es15.7, f3.0)') &
         (xgen(i), rgen(i), zero, i = 1, numi_generatrix)
      close (lungen1)

      deallocate (xg, rg)

      end subroutine non_analytic

!     -----------------------------------------------------------------------
      real function distance (x1, y1, x2, y2)  ! Obvious purpose & arguments
!     -----------------------------------------------------------------------

      real, intent (in) :: x1, y1, x2, y2

      distance = sqrt ((x1 - x2)**2 + (y1 - y2)**2)

      end function distance

!     -----------------------------------------------------------------------
      subroutine cone_indices (i1, i2, L)  ! Retrieve cone indices from the
                                           ! initial discretization
!     -----------------------------------------------------------------------

      integer, intent (out) :: i1, i2      ! Indices defining cone portion
      real,    intent (out) :: L           ! Corresponding line length

!     Locate the interim generatrix indices defining the cone segment:

      i1 = numi_internal/5
      call interval (numi_internal, xfore, capsule%xt_nose, one, i1)
      if (abs (xfore(i1+1) - capsule%xt_nose) < &
          abs (xfore(i1)   - capsule%xt_nose)) i1 = i1 + 1  ! Remote possibility

      i2 = 4*numi_internal/5
      call interval (numi_internal, xfore, capsule%xt_shoulder, one, i2)
      if (abs (xfore(i2+1) - capsule%xt_shoulder) < &
          abs (xfore(i2)   - capsule%xt_shoulder)) i2 = i2 + 1

      L = distance (xfore(i1), rfore(i1), xfore(i2), rfore(i2))

      end subroutine cone_indices

!     -----------------------------------------------------------------------
      subroutine bend_rib ()  ! Option to bend cone portion of the generatrix
!     -----------------------------------------------------------------------

!     Local variables:

      integer :: i, i1, i2, ier, ikeep, nicat, nleft
      real    :: a, deflection, dx, dy, L, semispan, slope, theta
      real    :: x1, y1, xc, yc, xt, yt
      real, allocatable, dimension (:) :: xcat, ycat, xmorph, ymorph

!     Execution:

!     Note that this operates on the raw initial discretization; curvature-based
!     redistribution follows elsewhere.

!     Locate the interim generatrix indices defining the cone segment:

      call cone_indices (i1, i2, L);  semispan = L * half

      if (frustum_radius == zero) then
!        i1 defines it at the nose tangency pt.
      else
         frustum_radius = max (frustum_radius, capsule%rt_nose)
         frustum_radius = min (frustum_radius, capsule%rt_shoulder*0.9)
         call interval (i2, rfore, frustum_radius, one, i1)
         if (abs (rfore(i1+1) - frustum_radius) < &
             abs (rfore(i1)   - frustum_radius)) i1 = i1 + 1
      end if
      frustum_radius = rfore(i1)  ! Needed again in morph_to_umbrella

      write (luncrt, '(a, es15.7, a)') &
         'Active frustum radius:', frustum_radius, ' meters', &
         '                      ', frustum_radius/meters_per_inch, ' inches'
      write (luncrt, '(a, 2i4)') &
         'First & last raw forebody indices bounding any deflections:', i1, i2

      deflection = abs (rib_deflection)
      if (deflection == zero) go to 99

      nicat = i2 - i1 + 1

!     Construct both halves of a catenary curve "suspended" from the x-axis:

      allocate (xcat(nicat), ycat(nicat))

      call catenary_grid (semispan, deflection, zero, 3, nicat, xcat, ycat, &
                          a, ier)
      if (ier /= 0) stop

      if (rib_deflection < zero) ycat(:) = -ycat(:)

!     Align the cone segment with the x-axis where the catenary lies:

      x1 = xfore(i1);  y1 = rfore(i1)

      call rotate2d (nicat, xfore(i1), rfore(i1), -half_cone_angle, x1, y1)

      xfore(i1:i2) = xfore(i1:i2) - (x1 + semispan)

      do i = i1, i2
         rfore(i) = ycat(i-i1+1)
      end do

      if (capsule%radius_shoulder > zero) then  ! Blend with the shoulder circle

         xc = capsule%x_base
         yc = capsule%rcenter_shoulder
         call rotate2d (1, xc, yc, -half_cone_angle, x1, y1)
         xc = xc - (x1 + semispan)
         yc = yc -  y1

!        Find where the end slope of the catenary matches the shoulder slope:

         slope = (ycat(nicat) - ycat(nicat-1)) / (xcat(nicat) - xcat(nicat-1))
         call circle_slope (xc, yc, capsule%radius_shoulder, slope, theta, &
                            xt, yt)

!        Morph this portion of the cone to (nearly) match slopes at (xt, yt):

         allocate (xmorph(nicat), ymorph(nicat))

         xmorph(1) = xfore(i1)
         ymorph(1) = rfore(i1)
         xmorph(nicat) = xt
         ymorph(nicat) = yt

         call nuline2d (1, nicat, xfore(i1), rfore(i1), xmorph, ymorph)

         xfore(i1:i2) = xmorph(:)
         rfore(i1:i2) = ymorph(:)

         deallocate (xmorph, ymorph)

      end if

      deallocate (xcat, ycat)

!     Transform the bent rib segment back to real space:

      dx = x1 - xfore(i1)
      xfore(i1:i2) = xfore(i1:i2) + dx
      dy = y1 - rfore(i1)
      rfore(i1:i2) = rfore(i1:i2) + dy

      call rotate2d (nicat, xfore(i1), rfore(i1), half_cone_angle, x1, y1)

      if (rib_deflection < zero) then  ! There may be some overlap
         ikeep = i2
         do i = i2 + 1, numi_internal
            ikeep = ikeep + 1
            if (xfore(i2) < xfore(i)) exit
         end do
         if (ikeep /= i2 + 1) then
            nleft = numi_internal - ikeep + 1
            xfore(i2+1:i2+nleft) = xfore(ikeep:numi_internal)
            rfore(i2+1:i2+nleft) = rfore(ikeep:numi_internal)
            numi_internal = i2 + nleft
         end if
      end if

  99  return

      end subroutine bend_rib

!     -----------------------------------------------------------------------
      subroutine toroid_simulation ()  ! Impose a ripple effect on generatrix
!     -----------------------------------------------------------------------

!     Equal-diameter touching toroids are envisaged with centers spanning the
!     distance between the two tangency points defining the cone segment.
!     The shoulder toroid counts as an extra one here, with ntoroids segments
!     between the two tangency points to be filled with a repeating cosine.
!     If the common diameter implied by ntoroids differs from the shoulder
!     diameter, the last section to the shoulder is morphed to take up any
!     slack in the initial sequence of equal cosine curve segments.
!     The input raw generatrix may have its number of points (numi_internal)
!     reduced slightly upon return.

!     Local variables:

      integer :: i, i1, i2, ishift, n, ndiameter, nshoulder
      real    :: diameter, dx, L, r, xshift
      real, allocatable, dimension (:) :: xcosine, ycosine

!     Execution:

!     Locate the interim generatrix indices defining the cone segment:

      call cone_indices (i1, i2, L)

      nshoulder = numi_internal - i2
      diameter  = L / real (ntoroids)
      ndiameter = 1 + int ((i2 - i1) / ntoroids)  ! Won't exceed available room

!     Construct a full cosine period on the interval [0, diameter]:

      r = 360. / diameter
      dx = diameter / real (ndiameter - 1)

      allocate (xcosine(numi_internal), ycosine(numi_internal))  ! > enough

      do i = 1, ndiameter
         xcosine(i) = dx*real (i - 1)
         ycosine(i) = (half*peak_ripple)*(cosd (r*xcosine(i)) - one)
      end do
      ycosine(1) = zero;  ycosine(ndiameter) = zero  ! Make sure of it

!     Replicate it for the rest of the toroids:

      do n = 2, ntoroids
         xshift =   diameter * real (n - 1)
         ishift = (ndiameter - 1) * (n - 1)
         do i = 2, ndiameter
            xcosine(i + ishift) = xcosine(i) + xshift
            ycosine(i + ishift) = ycosine(i)
         end do
      end do

      n = 1 + ntoroids * (ndiameter - 1)

      call rotate2d (n, xcosine, ycosine, capsule%half_cone_angle, zero, zero)

      do i = 1, n
         xcosine(i) = xcosine(i) + capsule%xt_nose;  xfore(i+i1-1) = xcosine(i)
         ycosine(i) = ycosine(i) + capsule%rt_nose;  rfore(i+i1-1) = ycosine(i)
      end do

!     Adjust for any difference between shoulder radius and toroid radius.
!     NULINE2D expects the same indexing in the input and output arrays.

      xfore(n+i1-1) = capsule%xt_shoulder
      rfore(n+i1-1) = capsule%rt_shoulder

      call nuline2d (n-ndiameter+1, n, xcosine, ycosine, xfore(i1), rfore(i1))

      n = n + i1 - 1
      xfore(n:n+nshoulder) = xfore(i2:numi_internal)
      rfore(n:n+nshoulder) = rfore(i2:numi_internal)
      numi_internal = n + nshoulder

      deallocate (xcosine, ycosine)

      end subroutine toroid_simulation

!     -----------------------------------------------------------------------
      subroutine round_vertices ()  ! Round off vertices as specified
!     -----------------------------------------------------------------------

!     Local constants:

      real,      parameter :: tol = 1.e-7     ! Awkwardness in replacing points
      character, parameter :: method*1 = 'B'  ! Loose fits

!     Local variables:

      integer           :: i, i2, ier, ikeep, in1, in2, iv, j, left, &
                           nafter, new_n, ni, nround, nu
      real              :: ds1, ds2, eps, x1, r1, x2, r2, x3, r3, r, radius, &
                           slast, stangent, total
      real              :: Rs, Rb, rsum, costh, sinth, theta, xo, xt, rt, &
                           x0, r0, xv, rv    ! SRC special-handling variables
      real, allocatable :: xc(:), rc(:), sc(:), snew(:), vnew(:), vnewp(:)
      logical           :: SRC_case, radii

!     Execution:

      new_n = numi
      do j = 2, nvertices
         if (ds_round(j) == zero) ni_round(j) = 0  ! Make sure of it
         new_n = new_n + ni_round(j)               ! More than enough
      end do

      nu = (numj - 1) / 4      ! In case there's no rounding

      if (new_n == numi .or. rounding_mode == 0) then  ! No rounding to do

!        In case the curvature-based redistribution has trouble, save its input:
         open  (lungen0, file='raw_generatrix_aft.dat', status='unknown')
         write (lungen0, '(2es20.12)') (xgen(i), rgen(i), i = 1, numi)
         close (lungen0)
         SRC_case = false
         go to 90
      end if  ! Else this file is written after the rounding

      radii = rounding_mode == 2  ! Else 1 => ds_round(:) = tangent lengths

      allocate (sgen(numi))
      call chords2d (numi, xgen, rgen, false, total, sgen)

      if (verbose) then
         write (luncrt, '(/, a)') 'Raw generatrix before vertex rounding: x,r,s'
         write (luncrt, '(i4, 3es15.7)') &
            (i, xgen(i), rgen(i), sgen(i), i = 1, numi)
      end if

      allocate (xgen_new(new_n), rgen_new(new_n), sgen_new(new_n))
      xgen_new(1:numi) = xgen(:)
      rgen_new(1:numi) = rgen(:)
      sgen_new(1:numi) = sgen(:)
      new_n = numi

      left = 1

      do j = 2, nvertices

         if (ds_round(j) < zero) cycle ! Special SRC case has to be done later

         if (j == nvertices) then
            if (.not. collapsed) exit  ! Else we may do a half-round
            nround = ni_round(j)
            if (nround > 0) then
               if (mod (nround, 2) == 0) ni_round(j) = nround + 1
            end if
         end if

!!       write (luncrt, '(a, 2i5)') 'j, ni_round(j):', j, ni_round(j)

         nround = ni_round(j);  if (nround == 0) cycle

         allocate (xc(nround), rc(nround), sc(nround))

         iv = ivertex(j)
         x3 = xgen(iv)
         r3 = rgen(iv)

!        Even if rounding_mode = 2, treat ds_round(:) as a tangent length here,
!        since the arc routine still needs points equidistant from the vertex:

!        Locate tangent point (or similar) left of vertex:

         stangent = sgen(iv) - ds_round(j)
         if (radii) then  ! Guard against large radius as for SRC aft body
            stangent = sgen(iv) - (sgen(iv) - sgen(iv-1))*0.9
         end if

         call interval (new_n, sgen_new, stangent, one, left)

         r  = (stangent - sgen_new(left)) / (sgen_new(left+1) - sgen_new(left))
         x1 = (one - r)*xgen_new(left)  +  r*xgen_new(left+1)
         r1 = (one - r)*rgen_new(left)  +  r*rgen_new(left+1)

         in1 = left;  if (r > zero) in1 = left + 1  ! First point to replace
         if (j == 2) imatch = left  ! Used in match_forebody

!!       write (10, '(a, i4, 3es15.7)') &
!!          'left, x1, r1, stangent:', left, x1, r1, stangent

!        Repeat for tangent point right of vertex:

!!       stangent = sgen(iv) + ds_round(j)  ! In view of SRC safeguard
         stangent = sgen(iv)*2. - stangent

         call interval (new_n, sgen_new, stangent, one, left)

         r  = (stangent - sgen_new(left)) / (sgen_new(left+1) - sgen_new(left))
         x2 = (one - r)*xgen_new(left)  +  r*xgen_new(left+1)
         r2 = (one - r)*rgen_new(left)  +  r*rgen_new(left+1)

         in2 = left;  if (r == one) in2 = left + 1  ! Last point to replace

!!       write (10, '(a, i4, 3es15.7)') &
!!          'left, x2, r2, stangent:', left, x2, r2, stangent

         if (j == nvertices) then  ! Kludge
            x2  =  x1  ! Reflect point 1 about Ox
            r2  = -r1
            in2 = new_n  ! Last point to replace
         end if

!        Discretize the indicated circular arc:

         if (radii) then
            radius = ds_round(j)  ! It will be returned as the tangent length
         else
            radius = zero
         end if

!!       write (luncrt, '(a, es17.8)') &
!!          'x1     ', x1, &
!!          'r1     ', r1, &
!!          'x2     ', x2, &
!!          'r2     ', r2, &
!!          'x3     ', x3, &
!!          'r3     ', r3, &
!!          'radius ', radius

         call circular_arc (x1, r1, x2, r2, x3, r3, radius, nround, xc, rc)

         if (radii) then  ! Figure out the index range to replace as above
            write (luncrt, '(a, es15.7)') 'Implied tangent length:', radius
            stangent = sgen(iv) - radius  ! Radius has been updated accordingly
            call interval (new_n, sgen_new, stangent, one, left)
            in1 = left  ! First point replaced
            eps = tol * radius
            if (abs (stangent - sgen_new(left)) < eps) then  ! Handle round-off
               stangent = sgen_new(left)
            else if (abs (stangent - sgen_new(left+1)) < eps) then
               stangent = sgen_new(left+1);  in1 = left + 1  ! 1st to replace
            end if
            if (stangent > sgen_new(left)) in1 = left + 1    ! 1st to replace
!!          write (luncrt, '(a, 2i4, 3es20.12)') &
!!             'j, in1, stangent, sgen_new(left:left+1):', &
!!              j, in1, stangent, sgen_new(left:left+1)
!!          write (10, '(a, i4, es15.7)') 'left, stangent:', left, stangent

            if (j < nvertices) then  ! Likewise for aft of the vertex
               stangent = sgen(iv) + radius
               call interval (new_n, sgen_new, stangent, one, left)
               in2 = left  ! Last point replaced
               if (abs (stangent - sgen_new(left)) < eps) then ! Treat round-off
                  stangent = sgen_new(left)
               else if (abs (stangent - sgen_new(left+1)) < eps) then
                  stangent = sgen_new(left+1);  in2 = left + 1 ! Last replaced
               end if
!!             write (luncrt, '(a, 2i4, 3es20.12)') &
!!                'j, in2, stangent, sgen_new(left:left+1):', &
!!                 j, in2, stangent, sgen_new(left:left+1)
!!             write (10, '(a, i4, es15.7)') 'left, stangent:', left, stangent
            end if
         else
            write (luncrt, '(a, es15.7)') 'Implied vertex radius:', radius
         end if

!        Substitute the arc points for the indicated rectilinear points:

         if (j == nvertices) then  ! Collapsed case ends with a half-round
            nround = nround / 2 + 1
            in2 = new_n
         end if

!!       write (luncrt, '(a, 2i4)') 'j, nround:', j, nround
         call replace  (nround, xc, new_n, xgen_new, in1, in2, nafter)
         call replace  (nround, rc, new_n, rgen_new, in1, in2, nafter)

         if (j < nvertices) then  ! Prepare for a next vertex
            call chords2d (nround, xc, rc, false, total, sc)
            sc(:) = sc(:) + sgen_new(in1)
            slast = sc(nround)  ! End of rounding, for locating ni_regrid below
            call replace  (nround, sc, new_n, sgen_new, in1, in2, nafter)
         else ! j = nvertices: must be rounding the end of the collapsed body
            slast = sgen_new(in1)  ! Start of the last rounding
         end if

!!       write (luncrt, '(a, i1, 3i4, es15.7)') &
!!          'j, in1, in2, nafter, slast:', j, in1, in2, nafter, slast

         new_n = nafter

         deallocate (xc, rc, sc)

      end do  ! Next vertex to round

!     Special handling of a SRC (Sample Return Capsule) situation: the aft-most
!     rounded vertex leaves no room for standard rounding of the preceding
!     vertex.  Some rounding there (even if it's not in the true geometry) is
!     desirable for the sake of the hyperbolic volume gridding.  The plain
!     segment meeting the large rounded aft body is assumed to be vertical.
!     Doing the following in-line is less messy than moving it to a subroutine.

      if (ncones ==1) then
         SRC_case = false
      else
         SRC_case = ds_round(ncones) < zero          .and. &
                    ni_round(ncones+1) > 0           .and. &
                    radii                            .and. &
                    x_aft(ncones) == x_aft(ncones-1) .and. &
                    cone_angle_aft(ncones) >= zero   .and. &
                    cone_angle_aft(ncones) < 90.
         write (luncrt, '(/, a, l2)') 'SRC_case: ', SRC_case
      end if

      if (SRC_case .and. ni_round(ncones) > 0) then

         xv    = x_aft(ncones)        ! Vertical segment abscissa
         Rb    =  ds_round(ncones+1)  ! Big radius of rounded aft body vertex
         Rs    = -ds_round(ncones)    ! Small rounding radius that had no room
         rsum  = Rs + Rb              !    by the usual rounding method
         xo    = x_aft(ncones+2) - Rb ! At center of large aft circle
         costh = (Rs + xv - xo) / rsum
         theta = acosd (costh)
         sinth = sind (theta)
         xt    = Rb * costh + xo      ! Tangent point where circles touch
         rt    = Rb * sinth
         x0    = rsum * costh + xo    ! Rounding circle center (x0 not needed)
         r0    = rsum * sinth
         rv    = r0 - Rs * tand (half * theta) ! Vertex where tangents meet

         nround = ni_round(ncones)
         allocate (xc(nround), rc(nround))
         radius = Rs                  ! Overwritten by tangent length

!!       write (luncrt, '(a, es17.8)') &
!!          'xv     ', xv, &
!!          'r0     ', r0, &
!!          'xt     ', xt, &
!!          'rt     ', rt, &
!!          'xv     ', xv, &
!!          'rv     ', rv, &
!!          'radius ', radius

         call circular_arc (xv, r0, xt, rt, xv, rv, radius, nround, xc, rc)

         write (luncrt, '(a, es15.7)') 'Implied tangent length:', radius

!        Locate the index range to replace:

         do i = 1, new_n
            in1 = i
            if (rgen_new(i) <= r0) exit
         end do

         in2 = in1 + 1
         do i = in1, new_n
            if (xgen_new(i) > xt) exit
            in2 = i
         end do

!!       write (luncrt, '(a, i4, es15.7)') 'in1, r0:', in1, r0, &
!!                                         'in2, xt:', in2, xt

         call replace (nround, xc, new_n, xgen_new, in1, in2, nafter)
         call replace (nround, rc, new_n, rgen_new, in1, in2, nafter)
         new_n = nafter

         deallocate (xc, rc)

      end if  ! SRC kludge

!     In case the curvature-based redistribution has trouble, save its input:

      open  (lungen0, file='raw_generatrix_aft.dat', status='unknown')
      write (lungen0, '(2es20.12)') (xgen_new(i), rgen_new(i), i = 1, new_n)
      close (lungen0)

!     Impose curvature-based spacing for the original numi number of points:

      power = power_aft

      call curvdis2 (new_n, xgen_new, rgen_new, numi, power, ismooth_aft, &
                     lunprint, false, xgen, rgen, ier)
      if (ier /= 0) then
         write (luncrt, '(a, i3, f5.1, /, a)') &
            'CURVDIS2 trouble (aft body).  IER, POWER:', ier, power, 'Aborting.'
         stop
      end if

!     No easy way of knowing the index of the last [possibly rounded] off-axis
!     vertex now, so search for the point nearest to it:

!     For the non-rounded case, we'll regrid the whole last segment:

      if (.not. collapsed) then
         do i = 2, numi
            if (xgen(i) >= x_aft(nvertices-1)) then
               ni_regrid = numi - i + 1  ! The i indices are reversed before the
               ikeep = i                 ! regridding is done
               exit
            end if
         end do
      else
         do i = numi, 2, -1  ! In case the last input cone angle was ~90.
            if (xgen(i) <= x_aft(nvertices-1)) then
               ni_regrid = numi - i + 1
               ikeep = i
               exit
            end if
         end do
      end if

      write (luncrt, '(a, i5)') 'ni_regrid(1):', ni_regrid

!     In case of rounding at either or both of the last two vertices, determine
!     the end of the rounding at the second last vertex and/or the start of the
!     rounding of the last collapsed vertex, and adjust ni_regrid accordingly
!     so that 1:ni_regrid (after reversing later) covers essentially uniformly-
!     spaced points so far.  If this is true, simply adjusting numi can give
!     smooth spacing across ni_grid.  (Hard to explain.)

      call chords2d (numi, xgen, rgen, false, total, sgen)
      left = 1
      call interval (numi, sgen, slast, one, left)
      write (luncrt, '(a, i5)') &
        'Index of end of rounding at 2nd-last (or last collapsed) vertex:', left

!     We don't want to redistribute across a rounded vertex, while capturing a
!     non-rounded vertex with the regridded quadrants makes most sense.)

      if (ni_round(nvertices-1) > 0) then   ! Keep clear of any rounded vertex
         ni_regrid = (numi - left) * 2/3    ! Safety margin
         write (luncrt, '(a, i5)') 'ni_regrid(2):', ni_regrid
      end if

!     However, if the last, collapsed vertex is rounded, confine the regridding
!     to the rounded portion for better hyperbolic volume grid results:

      if (collapsed .and. ni_round(nvertices) > 0) then
         ni_regrid = numi - (left + 3)  ! Allow for the smoothing in CURVDIS
         write (luncrt, '(a, i5)') 'ni_regrid(3):', ni_regrid
      end if

      if (sting_case) then

         if (sting_stretch_multiplier /= zero) then

!           Restretch the points along the sting part of the generatrix:

            ds1 = sgen(left+1) - sgen(left)
            ds2 = sting_stretch_multiplier * (total - sgen(numi-1))
            ni  = numi - left + 1
            allocate (snew(ni), vnew(ni), vnewp(ni))
            snew(1)  = sgen(left)
            snew(ni) = total
            call vinokur (1, ni, ds1, ds2, snew, luncrt, ier)
            call lcsfit (ni, sgen(left), xgen(left), true, method, ni, snew, &
                         vnew, vnewp)
            xgen(left:numi) = vnew(:)
            call lcsfit (ni, sgen(left), rgen(left), true, method, ni, snew, &
                         vnew, vnewp)
            rgen(left:numi) = vnew(:)
            deallocate  (snew, vnew, vnewp)
         end if

      else

!        Since numj is what determines the number of points on the regridded aft
!        patches of the revolved generatrix, and they're essentially uniform,
!        there is no way of ensuring smooth spacing variation at ni_grid other
!        than by trial and error with numi.  Therefore, blend the spacing of the
!        portion between the second last vertex (allowing for any rounding) and
!        the ni_grid location (but only if there's room; the rounding of a
!        collapsed last vertex may take up the whole last segment):

         i2 = numi - ni_regrid + 1  ! Index <--> ni_regrid before i is reversed
         ni = i2 - ikeep + 1
!!       write (luncrt, '(a, 3i5)') 'ikeep, i2, ni:', ikeep, i2, ni

         if (ni > 3 .and. ni_regrid_aft_body == 0) then
            write (luncrt, '(/, (a))') &
               'Redistributing to blend at default aft-body regrid index.', &
               'Entering nonzero ni_regrid_aft_body may be preferable.'
            ds1 = sgen(ikeep+1) - sgen(ikeep)
            ds2 = (total - sgen(i2)) / real (nu)  ! ~du, aft seg.
            allocate (snew(ni), vnew(ni), vnewp(ni))
            snew(1)  = sgen(ikeep)
            snew(ni) = sgen(i2)
            call vinokur (1, ni, ds1, ds2, snew, luncrt, ier)
            call lcsfit (ni, sgen(ikeep), xgen(ikeep), true, method, ni, snew, &
                         vnew, vnewp)
            xgen(ikeep:i2) = vnew(:)
            call lcsfit (ni, sgen(ikeep), rgen(ikeep), true, method, ni, snew, &
                         vnew, vnewp)
            rgen(ikeep:i2) = vnew(:)
            deallocate  (snew, vnew, vnewp)
         end if

      end if

      deallocate (sgen, xgen_new, rgen_new, sgen_new)

 90   continue

!     Another awkward situation occurs when the rounding is the whole vertex.
!     But then, the right choice for ni_regrid is simply ni_template (adjusted):

      if (SRC_case) ni_regrid = nu + 1

      if (ni_regrid_aft_body > 0) ni_regrid = ni_regrid_aft_body  ! Last resort

      write (luncrt, '(a, i5)') 'Final ni_regrid:', ni_regrid

      end subroutine round_vertices

!     --------------------------------------------------------------------------
      subroutine replace (nreplace, xreplace, nbefore, x, i1, i2, nafter)
!
!     Replace values x(i1:i2) with values xreplace(1:nreplace) by shifting to
!     the right to make room first.  If duplicates are not intended, the calling
!     program should adjust what's being inserted, and where, accordingly.
!     --------------------------------------------------------------------------

!     Arguments:

      integer, intent (in)    :: nreplace           ! # real values to plug in
      real,    intent (in)    :: xreplace(nreplace) ! Values to insert
      integer, intent (in)    :: nbefore            ! Initial # values in x(:)
      real,    intent (inout) :: x(*)               ! Array being inserted into
      integer, intent (in)    :: i1, i2             ! x(i1:i2) are replaced
      integer, intent (out)   :: nafter             ! # x(:) values afterwards

!     Local variables:

      integer :: nextra

!     Execution:

      nextra = nreplace - (i2 - i1 + 1)
      nafter = nbefore + nextra
      x(i2+1+nextra:nafter) = x(i2+1:nbefore)
      x(i1:i1+nreplace-1)   = xreplace(:)

      end subroutine replace

!     -----------------------------------------------------------------------
      subroutine match_forebody ()

!     Blend the first aft-body generatrix segment spacing with the forebody.
!     -----------------------------------------------------------------------

!     Local constants:

      character, parameter :: method*1 = 'B'  ! Loose fits

!     Local variables:

      integer :: ier, ni
      real    :: ds2, total
      real, allocatable :: chords(:), snew(:), vnew(:), vnewp(:)

!     Execution:

      ni = imatch * 4 / 5   ! Anything reasonable, away from any rounding
      if (ni < 5) go to 99  ! The rounding must be too near the juncture

      allocate (chords(ni), snew(ni), vnew(ni), vnewp(ni))

      ds2 = sqrt ((xgen(ni) - xgen(ni-1))**2 + (rgen(ni) - rgen(ni-1))**2)

      call chords2d (ni, xgen, rgen, false, total, chords)

!     Safeguard abnormal cases:

      if (ds_forebody > (total - ds2)*half) go to 90

      snew(1) = zero;  snew(ni) = total

      call vinokur (1, ni, ds_forebody, ds2, snew, luncrt, ier)

      call lcsfit (ni, chords, xgen, true, method, ni, snew, vnew, vnewp)
      xgen(1:ni) = vnew(:)

      call lcsfit (ni, chords, rgen, true, method, ni, snew, vnew, vnewp)
      rgen(1:ni) = vnew(:)

 90   deallocate (chords, snew, vnew, vnewp)

 99   return

      end subroutine match_forebody

!     -----------------------------------------------------------------------
      subroutine uniform_aft_end (numi)

!     Avoid stretching towards the symmetry axis at the aft end of the
!     generatrix by imposing a heuristic number of uniformly spaced points
!     there and redistributing the remaining points back to ni_regrid_aft...
!     -----------------------------------------------------------------------

!     Argument:

      integer, intent (in) :: numi
                       ! Index of the last point in the aft body generatrix.
                       ! The analytic case reuses x/rgen(1:numi) for the aft
                       ! body discretization, while the nonanalytic case has
                       ! x/rgen(numi_generatrix) defined already.  This argu-
                       ! ment is passed accordingly in each of two calls.

!     Local constants:

      integer,   parameter :: nu = 6          ! # uniform pts. at the aft end
      character, parameter :: method*1 = 'B'  ! Loose fits

!     Local variables:

      integer :: i
      integer :: i1, i2, iu, ier, ni
      real    :: ds2, dsu, total
      real, allocatable :: chords(:), snew(:), vnew(:), vnewp(:)

!     Execution:

      ni = min (ni_regrid_aft_body, 3*nu)
      i1 = numi - ni + 1
      i2 = numi - nu + 1
      iu = ni   - nu + 1

!!!   write (6, '(a, (a, i5))') 'Subroutine uniform_aft_end:', &
!!!      'numi:', numi, &
!!!      'nu:  ', nu, &
!!!      'ni:  ', ni, &
!!!      'i1:  ', i1, &
!!!      'i2:  ', i2, &
!!!      'iu:  ', iu

      allocate (chords(ni), snew(ni), vnew(ni), vnewp(ni))

!!!   write (6, '(i5, 2es16.8)') (i, xgen(i), rgen(i), i = i1, numi)

      call chords2d (ni, xgen(i1), rgen(i1), false, total, chords)

!!!   write (6, '(i5, 3es16.8)') (i, xgen(i), rgen(i), chords(i-i1+1), i = i1, numi)

      call xgrid (nu, 0, chords(iu), total, snew(iu))  ! Uniform

!!!   write (6, '(i5, 2f12.7)') (i, snew(i-i1+1), rgen(i), i = i1, numi)

      dsu = total - snew(ni-1)
      snew(1) = zero;  snew(iu) = chords(iu)

!!!   write (6, '(a, es15.8)') 'chords(2):', chords(2), 'dsu:', dsu

      call vinokur (1, iu, chords(2), dsu, snew, luncrt, ier)

!!!   write (6, '(a, i5)') 'Vinokur ier:', ier
!!!   write (6, '(i5, f12.7)') (i, snew(i), i = 1, ni)

      call lcsfit (ni, chords, xgen(i1), true, method, ni, snew, vnew, vnewp)
      xgen(i1:numi-1) = vnew(1:ni-1)

      call lcsfit (ni, chords, rgen(i1), true, method, ni, snew, vnew, vnewp)
      rgen(i1:numi-1) = vnew(1:ni-1)

      deallocate (chords, snew, vnew, vnewp)

      end subroutine uniform_aft_end

!     ------------------------------------------------------------------------
      subroutine spoked_surface (icase)

!     Rotate the generatrix to produce a simple interim surface.
!     icase = 1 and 2 for the forebody and aft body respectively.
!     ------------------------------------------------------------------------

      integer, intent (in) :: icase   ! 1, 2 ==> forebody, aft body resp.
      integer :: i, j
      real    :: angle, dtheta, px, py, pz, qx, qy, qz

!     Derive the simple spoked grid by rotating the generatrix:

      spoke_grid%ni = numi
      spoke_grid%nj = numj
      spoke_grid%nk = 1

      allocate (spoke_grid%x(numi,numj,1), spoke_grid%y(numi,numj,1), &
                spoke_grid%z(numi,numj,1))

      dtheta = 2.*asind (one) / real (numj - 1)

      px = zero;  py = zero;  pz = zero  ! Rotation axis
      qx = one;   qy = zero;  qz = zero

      do j = 1, numj
         do i = 1, numi
            spoke_grid%x(i,j,1) =  xgen(i)
            spoke_grid%y(i,j,1) =  zero
            spoke_grid%z(i,j,1) = -rgen(i)  ! We start at the 6 o'clock position
         end do
      end do

      do j = 2, numj - 1
         angle = dtheta * real (j - 1)
         call rotate3d (numi, spoke_grid%x(:,j,1), spoke_grid%y(:,j,1), &
                        spoke_grid%z(:,j,1), angle, px, py, pz, qx, qy, qz)
      end do

      do i = 1, numi
         spoke_grid%z(i,numj,1) = rgen(i)  ! Precisely
      end do

      if (icase == 1) then
         open (lunspoke, file='spoked_surface_fore.g', status='unknown')
      else
         open (lunspoke, file='spoked_surface_aft.g',  status='unknown')
      end if

      write (lunspoke, '(i1)') 1
      write (lunspoke, '(2i4, i2)') numi, numj, 1
      write (lunspoke, '(8es15.7)') spoke_grid%x, spoke_grid%y, spoke_grid%z
      close (lunspoke)

      end subroutine spoked_surface

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

!!!         open (lunchk, file='adjusted_template.g', status='unknown')
!!!         call xyz_write (lunchk, true, nb, template, ios)
!!!         close (lunchk)

         case (3) ! Two edges along a semicircle

            write (luncrt, '(/, a)') 'This case has not been implemented.'
            stop

      end select

      end subroutine adjust_template

!     --------------------------------------------------------------------------
      subroutine initial_regrid (nose)  ! Before projecting to the true surface
!
!     Replace the nose portion of the spoked surface with patches to avoid the
!     singular point.  A planar template patch is morphed appropriately.
!     Only the patch transposes for right-handedness differ between the nose and
!     the nose-like reflected rounded aft body operated on.
!     --------------------------------------------------------------------------

!     Arguments:

      logical, intent (in) :: nose  ! T => nose; F => rounded aft body

!     Local constants:

      real,      parameter :: half = 0.5
      character, parameter :: method * 1 = 'B' ! Loose fits in ADJUSTN2

!     Local variables:

      integer :: i, ii, it, iu, j, jj, jt, ju, ni, nj
      real    :: total, xfudge, yfudge, zfudge, yscale, zscale, zshift
      real, allocatable, dimension (:,:,:,:) :: uvw
      type (grid_type) :: nose_patch           ! Scaling template helps WARPQ3D

!     Execution:

!     A quarter circle with two edges on the circle is the only case handled.
!     Originally, the template was assumed to be uniform along all edges, but
!     only the circular edges must be so: allow for nonuniformity along the
!     straight edges in case tightening the spacing towards the apex helps.

      nblocks = 6

      if (.not. nose) then

!        Move the forebody results aside and retrieve the aft body spoked grid:

         do ib = 1, nblocks
            ni = new_grid(ib)%ni;  new_fore(ib)%ni = ni
            nj = new_grid(ib)%nj;  new_fore(ib)%nj = nj
                                   new_fore(ib)%nk = 1
            allocate (new_fore(ib)%x(ni,nj,1), new_fore(ib)%y(ni,nj,1), &
                      new_fore(ib)%z(ni,nj,1))

            new_fore(ib)%x(:,:,:) = new_grid(ib)%x(:,:,:)
            deallocate (new_grid(ib)%x)
            new_fore(ib)%y(:,:,:) = new_grid(ib)%y(:,:,:)
            deallocate (new_grid(ib)%y)
            new_fore(ib)%z(:,:,:) = new_grid(ib)%z(:,:,:)
            deallocate (new_grid(ib)%z)
         end do

         deallocate (new_grid)
         deallocate (spoke_grid%x, spoke_grid%y, spoke_grid%z)

         numi      = numi_aft_body
         ni_regrid = ni_regrid_aft_body

         allocate (spoke_grid%x(numi,numj,1), spoke_grid%y(numi,numj,1), &
                   spoke_grid%z(numi,numj,1))

         spoke_grid%x = spoke_aft%x
         spoke_grid%y = spoke_aft%y
         spoke_grid%z = spoke_aft%z

         if (.not. umbrella_case) then  ! Else see the save_results procedure
            deallocate (spoke_aft%x, spoke_aft%y, spoke_aft%z)
         end if

      end if

      ni_outer = numi - ni_regrid + 1  ! Number of spoke points not touched

!     Morph the template by setting up edges and filling the interior(s).
!     ni_template has already been adjusted to match the current point count
!     determined by numj.

      allocate (uvw(ni_template, nj_template, 1, 3))  ! Normalized arc lengths

      allocate (new_grid(nblocks))  ! 1-2 = upper & lower quarters of nose
                                    ! 3-6 = split outer grid, low to high

      new_grid(1:2)%ni = ni_template; new_grid(3:nblocks)%ni = ni_outer
      new_grid(1:2)%nj = nj_template; new_grid(3:nblocks)%nj = nj_template
      new_grid(:)%nk   = 1

      do ib = 1, nblocks
         call xyz_allocate (new_grid(ib), ios)
      end do

!     Shifting and scaling the template helps WARPQ3D behavior.
!     The template is assumed to be in quadrant 1 of the unit (y,z) circle.

      nose_patch%ni = template(1)%ni
      nose_patch%nj = template(1)%nj
      nose_patch%nk = 1

      call xyz_allocate (nose_patch, ios)

!     Upper quarter of nose:
!     !!!!!!!!!!!!!!!!!!!!!!

!     The lower quarter was coded first and is clearer because the spokes
!     are assumed to start at 6 o'clock and go clockwise (from upstream).

!     Transcribe the (roughly) circular edge from 9 o'clock to 10:30:

      jj = 1
      do j = nj_quarter, nj_quarter + nj_template - 1
         new_grid(1)%x(ni_template,jj,1) = spoke_grid%x(ni_regrid,j,1)
         new_grid(1)%y(ni_template,jj,1) = spoke_grid%y(ni_regrid,j,1)
         new_grid(1)%z(ni_template,jj,1) = spoke_grid%z(ni_regrid,j,1)
         jj = jj + 1
      end do

!     Edge from 10:30 to 12 o'clock:

      jj = nj_template
      do j = nj_quarter + nj_template - 1, numj
         new_grid(1)%x(jj,nj_template,1) = spoke_grid%x(ni_regrid,j,1)
         new_grid(1)%y(jj,nj_template,1) = spoke_grid%y(ni_regrid,j,1)
         new_grid(1)%z(jj,nj_template,1) = spoke_grid%z(ni_regrid,j,1)
         jj = jj - 1
      end do

!     The horizontal edge of new patch 1 is part of the middle spoke,
!     but we need to adjust the point count and relative spacing:

      allocate (x1(ni_regrid+1), y1(ni_regrid+1), z1(ni_regrid+1),  &
                x2(ni_template), y2(ni_template), z2(ni_template),  &
                snorm(ni_template))

      do i = 1, ni_regrid + 1
         x1(i) = spoke_grid%x(i,nj_quarter,1)
         y1(i) = spoke_grid%y(i,nj_quarter,1)
         z1(i) = spoke_grid%z(i,nj_quarter,1)
      end do

      yscale = y1(ni_regrid)

!     Impose the (adjusted) template's relative horizontal edge spacing:

      do i = 1, ni_template
         x2(i) = template(1)%x(i,1,1)
         y2(i) = template(1)%y(i,1,1)
         z2(i) = template(1)%z(i,1,1)
      end do

      call chords3d (ni_template, x2, y2, z2, true, total, snorm)

      call adjustn2 (ni_regrid + 1, 1, ni_regrid, x1, y1, z1, 1, ni_template, &
                     snorm, x2, y2, z2, method)  ! x2/y2/z2 are reused

      do i = 1, ni_template  ! This edge is common to new patch 2
         new_grid(1)%x(i,1,1) = x2(i);  new_grid(2)%x(i,1,1) = x2(i)
         new_grid(1)%y(i,1,1) = y2(i);  new_grid(2)%y(i,1,1) = y2(i)
         new_grid(1)%z(i,1,1) = z2(i);  new_grid(2)%z(i,1,1) = z2(i)
      end do

!     The vertical edge of new patch 1 is part of the last spoke,
!     but we need to adjust the point count and relative spacing:

      do i = 1, ni_regrid + 1
         x1(i) = spoke_grid%x(i,numj,1)
         y1(i) = spoke_grid%y(i,numj,1)
         z1(i) = spoke_grid%z(i,numj,1)
      end do

      zshift = z1(1)
      zscale = z1(ni_regrid) - zshift

!     Impose the (adjusted) template's relative vertical edge spacing,
!     assumed to be the same as for its horizontal edge:

      call adjustn2 (ni_regrid + 1, 1, ni_regrid, x1, y1, z1, 1, nj_template, &
                     snorm, x2, y2, z2, method)

      do j = 1, nj_template
         new_grid(1)%x(1,j,1) = x2(j)
         new_grid(1)%y(1,j,1) = y2(j)
         new_grid(1)%z(1,j,1) = z2(j)
      end do

      nose_patch%x = template(1)%x
      nose_patch%y = template(1)%y * yscale
      nose_patch%z = template(1)%z * zscale + zshift

!     Normalized arc lengths for the scaled (planar) template nose patch:

      call paramxyz (1, ni_template, 1, nj_template, 1, 1, &
                     1, ni_template, 1, nj_template, 1, 1, &
                     nose_patch%x, nose_patch%y, nose_patch%z, uvw)

!!!   write (10, *) 1
!!!   write (10, *) ni_template, nj_template, 1
!!!   write (10, *) uvw(:,:,1,1), uvw(:,:,1,2), uvw(:,:,1,3)

!     Morph the scaled template interior points to the interim new interior:

      call warpq3d2 (1, ni_template, 1, nj_template, 1, 1, &
                     1, ni_template, 1, nj_template, 1, 1, &
                     nose_patch%x, nose_patch%y, nose_patch%z, uvw, &
                     new_grid(1)%x, new_grid(1)%y, new_grid(1)%z)

      call transpose_patch (new_grid(1), 0)  ! Only patch 2 is rt.-handed

!     Lower quarter of nose:
!     !!!!!!!!!!!!!!!!!!!!!!

      do j = 1, ni_template  ! Circle edge from 6 o'clock to 7:30
         new_grid(2)%x(j,nj_template,1) = spoke_grid%x(ni_regrid,j,1)
         new_grid(2)%y(j,nj_template,1) = spoke_grid%y(ni_regrid,j,1)
         new_grid(2)%z(j,nj_template,1) = spoke_grid%z(ni_regrid,j,1)
      end do

      jj = nj_template
      do j = ni_template, nj_quarter  ! Edge from 7:30 to 9 o'clock
         new_grid(2)%x(ni_template,jj,1) = spoke_grid%x(ni_regrid,j,1)
         new_grid(2)%y(ni_template,jj,1) = spoke_grid%y(ni_regrid,j,1)
         new_grid(2)%z(ni_template,jj,1) = spoke_grid%z(ni_regrid,j,1)
         jj = jj - 1
      end do

!     The vertical edge of new patch 2 is part of the first spoke,
!     but we need to adjust the point count and spacing:

      do i = 1, ni_regrid + 1
         x1(i) = spoke_grid%x(i,1,1)
         y1(i) = spoke_grid%y(i,1,1)
         z1(i) = spoke_grid%z(i,1,1)
      end do

      zscale = z1(ni_regrid) - zshift

      call adjustn2 (ni_regrid + 1, 1, ni_regrid, x1, y1, z1, 1, nj_template, &
                     snorm, x2, y2, z2, method)

      do j = 1, nj_template
         new_grid(2)%x(1,j,1) = x2(j)
         new_grid(2)%y(1,j,1) = y2(j)
         new_grid(2)%z(1,j,1) = z2(j)
      end do

!     The horizontal edge of patch 2 is that of patch 1, copied above.

!!!   nose_patch%x = template(1)%x                   ! Same as above
!!!   nose_patch%y = template(1)%y * yscale          !   "   "   "
      nose_patch%z = template(1)%z * zscale + zshift

      call paramxyz (1, ni_template, 1, nj_template, 1, 1, &
                     1, ni_template, 1, nj_template, 1, 1, &
                     nose_patch%x, nose_patch%y, nose_patch%z, uvw)

!     Morph the template interior points to the interim new interior:

      call warpq3d2 (1, ni_template, 1, nj_template, 1, 1, &
                     1, ni_template, 1, nj_template, 1, 1, &
                     nose_patch%x, nose_patch%y, nose_patch%z, uvw, &
                     new_grid(2)%x, new_grid(2)%y, new_grid(2)%z)

      deallocate (x1, y1, z1, x2, y2, z2, snorm)

!     Transcribe the outer grid, split into new patches 3-6:
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do j = 1, ni_template
         ii = 0
         do i = ni_regrid, numi
            ii = ii + 1
            new_grid(3)%x(ii,j,1) = spoke_grid%x(i,j,1)
            new_grid(3)%y(ii,j,1) = spoke_grid%y(i,j,1)
            new_grid(3)%z(ii,j,1) = spoke_grid%z(i,j,1)
         end do
      end do

      jj = 0
      do j = ni_template, nj_quarter
         ii = 0
         jj = jj + 1
         do i = ni_regrid, numi
            ii = ii + 1
            new_grid(4)%x(ii,jj,1) = spoke_grid%x(i,j,1)
            new_grid(4)%y(ii,jj,1) = spoke_grid%y(i,j,1)
            new_grid(4)%z(ii,jj,1) = spoke_grid%z(i,j,1)
         end do
      end do

      jj = 0
      do j = nj_quarter, nj_quarter + nj_template - 1
         ii = 0
         jj = jj + 1
         do i = ni_regrid, numi
            ii = ii + 1
            new_grid(5)%x(ii,jj,1) = spoke_grid%x(i,j,1)
            new_grid(5)%y(ii,jj,1) = spoke_grid%y(i,j,1)
            new_grid(5)%z(ii,jj,1) = spoke_grid%z(i,j,1)
         end do
      end do

      jj = 0
      do j = numj - ni_template + 1, numj
         ii = 0
         jj = jj + 1
         do i = ni_regrid, numi
            ii = ii + 1
            new_grid(6)%x(ii,jj,1) = spoke_grid%x(i,j,1)
            new_grid(6)%y(ii,jj,1) = spoke_grid%y(i,j,1)
            new_grid(6)%z(ii,jj,1) = spoke_grid%z(i,j,1)
         end do
      end do

      do ib = 1, 6
         if ((nose .and. ib > 2) .or. (.not. nose .and. ib <= 2)) then
            call transpose_patch (new_grid(ib), 0)
         end if
      end do

      deallocate (uvw)

      if (.not. nose) then
         do ib = 1, nb
            deallocate (template(ib)%x, template(ib)%y, template(ib)%z)
         end do
      end if

      deallocate (nose_patch%x, nose_patch%y, nose_patch%z)

!     Fudge the corner points a little to avoid the nearly 180-degree
!     angles by averaging the first and second points on the outer spokes:

      if (nose) then
         it = 1;  jt = 2
         iu = ni_template;  ju = 1
      else  ! Because of different transposes
         it = 2;  jt = 1
         iu = 1;  ju = nj_template
      end if

      xfudge = (new_grid(6)%x(1,1,1) + new_grid(6)%x(it,jt,1))*half
      yfudge = (new_grid(6)%y(1,1,1) + new_grid(6)%y(it,jt,1))*half
      zfudge = (new_grid(6)%z(1,1,1) + new_grid(6)%z(it,jt,1))*half

      new_grid(1)%x(ni_template,ni_template,1) = xfudge
      new_grid(1)%y(ni_template,ni_template,1) = yfudge
      new_grid(1)%z(ni_template,ni_template,1) = zfudge

      new_grid(5)%x(iu,ju,1) = xfudge
      new_grid(5)%y(iu,ju,1) = yfudge
      new_grid(5)%z(iu,ju,1) = zfudge

      new_grid(6)%x(1,1,1) = xfudge
      new_grid(6)%y(1,1,1) = yfudge
      new_grid(6)%z(1,1,1) = zfudge

      xfudge = (new_grid(4)%x(1,1,1) + new_grid(4)%x(it,jt,1))*half
      yfudge = (new_grid(4)%y(1,1,1) + new_grid(4)%y(it,jt,1))*half
      zfudge = (new_grid(4)%z(1,1,1) + new_grid(4)%z(it,jt,1))*half

      new_grid(2)%x(ni_template,ni_template,1) = xfudge
      new_grid(2)%y(ni_template,ni_template,1) = yfudge
      new_grid(2)%z(ni_template,ni_template,1) = zfudge

      new_grid(3)%x(iu,ju,1) = xfudge
      new_grid(3)%y(iu,ju,1) = yfudge
      new_grid(3)%z(iu,ju,1) = zfudge

      new_grid(4)%x(1,1,1) = xfudge
      new_grid(4)%y(1,1,1) = yfudge
      new_grid(4)%z(1,1,1) = zfudge

      end subroutine initial_regrid

!     --------------------------------------------------------------------------
      subroutine morph_to_umbrella ()  ! Specialized regridding of cone/shoulder
!     --------------------------------------------------------------------------

!     Local constants:

      real,      parameter :: fraction = 0.50, half = 0.5, one = 1.0, zero = 0.0
      logical,   parameter :: false = .false., true = .true.
      character, parameter :: method*1  = 'B'  ! Loose local cubic interpolation
      character, parameter :: methodr*1 = 'M'  ! Tight (monotonic) rib rounding

!     Local variables:

      integer :: i, ier, ieval, ir, it1, it2, iter, itm, &
                 j, je1, jo1, je2, jo2, m, n, nicat, njsemim1, npt
      logical :: new, invert
      real    :: a, angle, curvature, ds, dsredo, dz, dzobtained, dztarget, &
                 peak_local_deflection, px, py, pz, qx, r, semispan, &
                 si, smid, stotal, sutotal, tol, xr, yr, zr, &
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

!     The original scheme in FOREBODY_REGRID worked with the (x,y,z) convention
!     of HEAT_SHIELD's spoked grid :(.  It is most expedient to swap y and z to
!     reuse the tricky "umbrella" gyrations, then swap back:

      call swap_coordinates (spoke_grid, 2, 3)

!     Now we have y up and z > 0 on the port side (right half looking along Ox).

!     First main step of two:  Morph the shoulder from rounded to polygonal.
!     ----------------------------------------------------------------------

!     Set up a half-edge segment as before-and-after work-space.  This will
!     become the first (in j) portion starting at 6 o'clock on the port side.
!     The symmetry plane is deliberately down the middle of bottom and top
!     segments so as to show the maximum deflection slice in the CFD solutions.

      n = numi_forebody
      semiangle = 180. / real (nedges)   ! <-> half a polygon edge
      nj_semi = 1 + (numj - 1) / nedges  ! # grid spokes per semiangle
      m = nj_semi;  njsemim1 = m - 1

      if (nedges * njsemim1 + 1 /= numj) then
         write (luncrt, '(a, 2i4)') &
            '# grid spokes not suited to this many edges:', numj, nedges
         write (luncrt, '(a)') 'Constraints:', &
            '   numj - 1 = must be divisible by both 4 and nedges.', &
            '   Aborting.'
         stop
      end if

      allocate (xold(n,m), yold(n,m), zold(n,m), &
                xnew(n,m), ynew(n,m), znew(n,m))

      xold(:,:) = spoke_grid%x(:,1:m,1);  xnew = xold
      yold(:,:) = spoke_grid%y(:,1:m,1);  ynew = yold
      zold(:,:) = spoke_grid%z(:,1:m,1);  znew = zold

!     Identify the first and last spoke indices on the conical flank.
!     This was done prior to the rib-bending option and before frustum_radius
!     was introduced as an input, but do it anyway.

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
         'First & last indices on conical flank:    ', it1, it2
      if (it1 * it2 == 0) stop

!     Determine it1 from the [possibly adjusted] value of input frustum_radius.
!     Note that we're now working with the curvature-based spoke-wise grid pts.

      call interval (it2, yold(:,1), -frustum_radius, -one, it1)  ! 6 o'clock

      it1 = it1 + 1  ! Because of how interval works on descending data

      write (luncrt, '(a, 2i4)') &
         'First & last indices bounding deflections:', it1, it2

!     Substitute straight lines from i = it2:ni_spoke and adjust the rest of
!     the spokes from i = it1:it2.  Just do it for one semiedge patch:

      do i = it2, ni_spoke
         xnew(i,1) = xold(i,nj_semi)
         ynew(i,1) = yold(i,nj_semi)
         znew(i,1) = zero
         do j = 2, njsemim1
            r = real (j - 1) / real (njsemim1)
            xnew(i,j) = xnew(i,1)
            ynew(i,j) = ynew(i,1)
            znew(i,j) = zold(i,nj_semi)*r
         end do
      end do

!     Adjust the spoke distribution between i = it1 and it2:

      do j = 1, njsemim1
         call nuline3d (it1, it2, xold(:,j), yold(:,j), zold(:,j), &
                                  xnew(:,j), ynew(:,j), znew(:,j))
      end do

!!    write (48, '(i1)')  1
!!    write (48, '(3i4)') ni_spoke, nj_semi, 1
!!    write (48, '(8es15.7)') &
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

         do j = 1, njsemim1

!           Construct a uniform catenary grid in (y,x) space (both halves)
!           spanning the distance between points (it1,j) and (it2,j):

            call point_A_to_point_B (ni_spoke, nj_semi, it1, j, it2, j, &
                                     xold, yold, zold, semispan)

            semispan = half * semispan
            peak_local_deflection = ycat(j) + rib_deflection  ! Mid-pt. values
            invert = peak_local_deflection < zero
            peak_local_deflection = abs (peak_local_deflection)

            call catenary_grid (semispan, peak_local_deflection, xold(it1,j), &
                                3, nicat, yicat, xicat, a, ier)  ! Assumes droop
            if (ier /= 0) stop

            if (invert) xicat(:) = 2.*xold(it1,j) - xicat(:)

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
!!       write (51, '(8es15.7)') &
!!          xnew(:,1:nj_semi), ynew(:,1:nj_semi), znew(:,1:nj_semi)

         deallocate (xcat, ycat, zcat, xicat, yicat, zicat, xu, yu, zu, su)

      end if  ! Deflection option

!     Resolve the ridges?  If so, we find the first i station where the rib
!     spacing exceeds half the specified rib width, and restretch from there.

      if (resolve_the_ridges) then

         n = ni_spoke;  m = nj_semi;  ir = n + 1

         call point_A_to_point_B (n, m, it1, m-1, it1, m, xold, yold, zold, ds)

!        ds here includes any deflection, so use dz as measure of rib thickness.
!        At the rib, this ignores effect of the rib angle (facet semiangle), but
!        is more convenient.  Ask for a slightly thicker rib if it's important.
!        Remember that resolving rib thickness this is crude at best.

         dz = abs (zold(it1,m) - zold(it1,m-1))
         dz = dz + dz

         write (luncrt, '(a, i4, a, f9.6, a, f9.5, a)') &
            'Minimum resolvable rib thickness at i = it1 =', it1, ':', dz, &
            ' meters (', dz/meters_per_inch, '")'

         rib_thickness = max (rib_thickness, dz)

!        Locate where the first grid spacing exceeds the specified rib semi-th.:

         do i = ni_regrid + 1, ni_spoke  ! Avoid the new nose cap

!!!         call point_A_to_point_B (n, m, i, m-1, i, m, xold, yold, zold, ds)
            dz = abs (zold(i,m) - zold(i,m-1))

            if (dz >= half*rib_thickness) then
                ir  = i
                call point_A_to_point_B (n, m, i, m-1, i, m, xold, yold, zold, &
                                         ds)  ! Gives starting guess for ds
                exit
            end if

         end do

         write (luncrt, '(a, f9.6, a, f9.5, a)') &
            'Adjusting rib thickness slightly:', &
            2.*dz, ' meters (', 2.*dz/meters_per_inch,'")'
         write (luncrt, '(a, i4)') 'Resolving rib thickness from i =', ir

         allocate (xu(m), yu(m), zu(m), su(m), snew(m))

         dztarget = abs (zold(ir-1,m) - zold(ir-1,m-1))
!!!      write (6, '(a, es13.5)') '   dztarget:', dztarget

         do i = ir, ni_spoke  ! Impose a Vinokur distribution at each i

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

!!!         dsobtained = sqrt ((px-xnew(i,m-1))**2 + (py-ynew(i,m-1))**2 + &
!!!                            (pz-znew(i,m-1))**2)
            dzobtained = abs (pz-znew(i,m-1))

!           Try to eliminate skewing caused by redistribution of redistributed
!           pts. (effect of approximate arc length as the parametric variable):

            dsredo = ds
            do iter = 1, 3

               dsredo = ds*dztarget/dzobtained

!!!            write (6, '(a, i4, i2, 3es13.5)') &
!!!               '   i,iter,ds,dzobtained,dsredo:', &
!!!                   i,iter,ds,dzobtained,dsredo

               call htdis4 (true, zero, sutotal, su(2), dsredo, m, snew, &
                            -luncrt, ier)

               new = true;  ieval = 1

               do j = 1, m
                  call plscrv3d (m, xu,yu,zu, su, method, new, false, snew(j), &
                                 ieval, px, py, pz, derivs)
                  xnew(i,j) = px
                  ynew(i,j) = py
                  znew(i,j) = pz
                  new = false
               end do

               dzobtained = abs (pz-znew(i,m-1))
               ds = dsredo

            end do  ! Repeat the correction?

!!!         write (6, '(a, i4, i2, 3es13.5)') &
!!!            '   i,iter,ds,dzobtained,dsredo:', &
!!!                i,iter,ds,dzobtained,dsredo
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

            call plscrv3d (4, xridge, yridge, zridge, sridge, methodr, true, &
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
            spoke_grid%x(it1+1:ni_spoke,jo1:jo2,1) = &
                       xtmp(it1+1:ni_spoke,1:nj_semi)
            spoke_grid%y(it1+1:ni_spoke,jo1:jo2,1) = &
                       ytmp(it1+1:ni_spoke,1:nj_semi)
            spoke_grid%z(it1+1:ni_spoke,jo1:jo2,1) = &
                       ztmp(it1+1:ni_spoke,1:nj_semi)
            jo1 = jo1 + 2*njsemim1
            jo2 = jo2 + 2*njsemim1
         else
            spoke_grid%x(it1+1:ni_spoke,je1:je2,1) = &
                       xtmp(it1+1:ni_spoke,1:nj_semi)
            spoke_grid%y(it1+1:ni_spoke,je1:je2,1) = &
                       ytmp(it1+1:ni_spoke,1:nj_semi)
            spoke_grid%z(it1+1:ni_spoke,je1:je2,1) = &
                       ztmp(it1+1:ni_spoke,1:nj_semi)
            je1 = je1 + 2*njsemim1
            je2 = je2 + 2*njsemim1
         end if

      end do

!     Revert to DPLR's usual (x,y,z) convention:

      call swap_coordinates (spoke_grid, 2, 3)

!!    write (50, '(i1)')  1
!!    write (50, '(3i4)') ni_spoke, nj_spoke, 1
!!    write (50, '(8es15.7)') spoke_grid%x, spoke_grid%y, spoke_grid%z

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

         dense_spokes(1)%x(1,:,1) = spoke_grid%x(1,1,1)
         dense_spokes(1)%y(1,:,1) = spoke_grid%y(1,1,1)
         dense_spokes(1)%z(1,:,1) = spoke_grid%z(1,1,1)

         do i = 2, ni_regrid
            x1(1)           =  spoke_grid%x(i,2,1)         ! Periodic data
            x1(2:nj_data-1) =  spoke_grid%x(i,:,1)
            x1(nj_data)     =  spoke_grid%x(i,nj_spoke-1,1)
            y1(1)           = -spoke_grid%y(i,2,1)
            y1(2:nj_data-1) =  spoke_grid%y(i,:,1)
            y1(nj_data)     = -spoke_grid%y(i,nj_spoke-1,1)
            z1(1)           =  spoke_grid%z(i,2,1)
            z1(2:nj_data-1) =  spoke_grid%z(i,:,1)
            z1(nj_data)     =  spoke_grid%z(i,nj_spoke-1,1)

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

         dense_spokes(1)%x(:,:,1) = spoke_grid%x(1:ni_regrid,:,1)
         dense_spokes(1)%y(:,:,1) = spoke_grid%y(1:ni_regrid,:,1)
         dense_spokes(1)%z(:,:,1) = spoke_grid%z(1:ni_regrid,:,1)

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

      write (luncrt, '(/, a, 2es14.6)') &
         'Largest and mean nose/aft patch projection distance:', dmax, dmean

      end subroutine project_to_nose

!     --------------------------------------------------------------------------
      subroutine save_result ()        ! And clean up
!     --------------------------------------------------------------------------

!     Local variables:

      integer :: i, ios, j, j1, j2, ni, nj
      real    :: total, zg
      real, allocatable, dimension (:) :: xg, yg, sg, growth

!     Execution:

      if (sting_case) then  ! Shuffle to append single aft block split into 4

         do i = 1, nblocks  ! Need to reallocate new_grid(:)
            ni = new_grid(i)%ni;  new_fore(i)%ni = ni
            nj = new_grid(i)%nj;  new_fore(i)%nj = nj
                                  new_fore(i)%nk = 1
            allocate (new_fore(i)%x(ni,nj,1), new_fore(i)%y(ni,nj,1), &
                      new_fore(i)%z(ni,nj,1))

            new_fore(i)%x(:,:,:) = new_grid(i)%x(:,:,:)
            deallocate (new_grid(i)%x)
            new_fore(i)%y(:,:,:) = new_grid(i)%y(:,:,:)
            deallocate (new_grid(i)%y)
            new_fore(i)%z(:,:,:) = new_grid(i)%z(:,:,:)
            deallocate (new_grid(i)%z)
         end do

         deallocate (new_grid)

         allocate (new_grid(nblocks+4))

         do i = 1, nblocks
            new_grid(i)%ni = new_fore(i)%ni
            new_grid(i)%nj = new_fore(i)%nj
            new_grid(i)%nk = 1
            call xyz_allocate (new_grid(i), ios)
            new_grid(i)%x(:,:,:) = new_fore(i)%x(:,:,:)
            new_grid(i)%y(:,:,:) = new_fore(i)%y(:,:,:)
            new_grid(i)%z(:,:,:) = new_fore(i)%z(:,:,:)
            deallocate (new_fore(i)%x, new_fore(i)%y, new_fore(i)%z)
         end do

         j1 = 1
         do i = nblocks + 1, nblocks + 4
            new_grid(i)%ni = numi_aft_body;  nj = new_grid(i-4)%ni
            new_grid(i)%nj = nj
            new_grid(i)%nk = 1
            call xyz_allocate (new_grid(i), ios)
            j2 = j1 + nj - 1
            new_grid(i)%x(:,:,:) = spoke_aft%x(:,j1:j2,:)
            new_grid(i)%y(:,:,:) = spoke_aft%y(:,j1:j2,:)
            new_grid(i)%z(:,:,:) = spoke_aft%z(:,j1:j2,:)
            j1 = j2
         end do

         deallocate (spoke_aft%x, spoke_aft%y, spoke_aft%z)

         nblocks = nblocks + 4

      else if (aft_body_too) then  ! Shuffle to concatenate blocks for the I/O

         do i = 1, nblocks
            ni = new_grid(i)%ni;  new_aft(i)%ni = ni
            nj = new_grid(i)%nj;  new_aft(i)%nj = nj
                                  new_aft(i)%nk = 1
            allocate (new_aft(i)%x(ni,nj,1), new_aft(i)%y(ni,nj,1), &
                      new_aft(i)%z(ni,nj,1))

            new_aft(i)%x(:,:,:) = new_grid(i)%x(:,:,:)
            deallocate (new_grid(i)%x)
            new_aft(i)%y(:,:,:) = new_grid(i)%y(:,:,:)
            deallocate (new_grid(i)%y)
            new_aft(i)%z(:,:,:) = new_grid(i)%z(:,:,:)
            deallocate (new_grid(i)%z)
         end do

         deallocate (new_grid);  allocate (new_grid(2*nblocks))

         do i = 1, nblocks
            new_grid(i)%ni = new_fore(i)%ni
            new_grid(i)%nj = new_fore(i)%nj
            new_grid(i)%nk = 1
            call xyz_allocate (new_grid(i), ios)
            new_grid(i)%x(:,:,:) = new_fore(i)%x(:,:,:)
            new_grid(i)%y(:,:,:) = new_fore(i)%y(:,:,:)
            new_grid(i)%z(:,:,:) = new_fore(i)%z(:,:,:)
            deallocate (new_fore(i)%x, new_fore(i)%y, new_fore(i)%z)

            j = 2*nblocks - i + 1
            new_grid(j)%ni = new_aft(i)%ni
            new_grid(j)%nj = new_aft(i)%nj
            new_grid(j)%nk = 1
            call xyz_allocate (new_grid(j), ios)
            new_grid(j)%x(:,:,:) = new_aft(i)%x(:,:,:)
            new_grid(j)%y(:,:,:) = new_aft(i)%y(:,:,:)
            new_grid(j)%z(:,:,:) = new_aft(i)%z(:,:,:)
            deallocate (new_aft(i)%x, new_aft(i)%y, new_aft(i)%z)
         end do

         nblocks = 2*nblocks

         if (umbrella_case) call morph_aft_to_umbrella ()  ! See procedure below

      end if

      call xyz_write (lunout, true, nblocks, new_grid, ios)

      if (ios /= 0) then
         write (luncrt, '(/, a)') 'Trouble saving regridded surface.'
      end if

      do ib = 1, nblocks
         deallocate (new_grid(ib)%x, new_grid(ib)%y, new_grid(ib)%z)
      end do

      deallocate (spoke_grid%x, spoke_grid%y, spoke_grid%z)

!     Tabulate the generatrix grid spacing growth rate distribution.
!     It is pragmatic to recover it by reading it:

      if (aft_body_too) then
         open (lungen1, file='full-body-generatrix.dat', status='old')
      else
         open (lungen1, file='forebody-generatrix.dat', status='old')
      end if

      read  (lungen1, *) ni
      allocate (xg(ni), yg(ni), sg(ni), growth(ni))
      read  (lungen1, *) (xg(i), yg(i), zg, i = 1, ni)
      close (lungen1)

      open  (lungen1, file='generatrix-growth-rates.dat', status='unknown')
      call chords2d (ni, xg, yg, false, total, sg)

      growth(1) = one  ! So the printed indices match
      do i = 2, ni - 1
         growth(i) = (sg(i+1) - sg(i)) / (sg(i) - sg(i-1))
      end do
      growth(ni) = one

!     Include (x, y) to help identify trouble spots:

      write (lungen1, '(a)') &
         '#       s   growth              x              y    i'
      write (lungen1, '(2f9.5, 2es15.7, i5)') &
         (sg(i), growth(i), xg(i), yg(i), i, i = 1, ni)
      close (lungen1)
      deallocate (xg, yg, sg, growth)

      end subroutine save_result

!     --------------------------------------------------------------------------
      subroutine morph_aft_to_umbrella ()

!     Specialized adapting of concave aft body to faceted forebody.
!     The shoulder region of each spoke (likely to be very angular) is moved as
!     a rigid sub-curve; the remainder to the center-body is likely to be
!     part of the underlying cone (i.e., ~ a straight line) and hence lends
!     itself to simple arc-length-based morphing.  Applying the latter to the
!     whole outer spoke would distort the shoulder somewhat, so we avoid that.
!
!     We also reconstruct the spoked form of the faceted surfaces in the hope
!     of being able to morph a 2D hyperbolic volume grid into a faceted volume.
!     (See a related option in Compress2D.)
!
!     Also:  If the forebody includes deflections (or comes from an input
!     generatrix) we seek ~constant fabric thickness by shifting the aft body
!     points between iu1 and iu2 in the x direction (only), since the underlying
!     aft body generatrix consists of just straight lines and circles.
!     --------------------------------------------------------------------------

      character (1), parameter :: method = 'B'  ! Loose local spline fits

      integer :: i, j, m, mj, n, ni, nj, njfore, iu0, iu1, iu2
      logical :: fix_thickness, new
      real    :: dx, dy, dz, raft, xinterp, unused
      real, allocatable, dimension (:) :: xold, yold, zold, xfore, rfore

!     Execution:

      iu0 = i0_umbrella
      iu1 = i1_umbrella
      iu2 = i2_umbrella

      if (iu1*iu2 == 0) go to 90  ! Too hard to set reliable defaults

      if (iu0 == 0) iu0 = iu1 - 20  ! For fat payload bay cases

      m  = 6               ! Forebody patch number matching first aft-body patch
      mj = new_grid(m)%nj  ! Spokewise point count on forebody patches
      ni = new_grid(7)%ni  ! Spokewise point count on aft-body patches
      nj = new_grid(7)%nj  ! Azimuthal point count along the rim

      if (iu1 > ni .or. iu2 > ni) go to 90  ! Old inputs /= 0 forgotten?

!!    fix_thickness = peak_deflection /= zero .or. &
!!                    rib_deflection  /= zero
      fix_thickness = analytic  ! Because of the faceting

      if (fix_thickness) then

!        Be sure aft iu1 doesn't correspond to a radius/z on a nose quadrant:

         if (new_grid(7)%z(iu1,nj,1) < new_grid(m)%z(nj,1,1)) then
            write (luncrt, '(a, 2es13.5)') &
               '*** Warning: Raise i1_umbrella; it is inboard of nose patch:', &
               new_grid(7)%z(iu1,nj,1), new_grid(m)%z(nj,1,1)
            go to 90
         end if

!        Find forebody patch 6 index <-> aft-body patch 7, iu2, nj (roughly):

         do j = 2, mj  ! Search along 12 o'clock spoke to find monotonic z range
            njfore = j
            if (new_grid(m)%z(nj,j,1) <= new_grid(m)%z(nj,j-1,1)) exit
         end do
         njfore = njfore - 1 ;  write (luncrt, '(a, i4)') 'njfore:', njfore

         allocate (xfore(njfore), rfore(njfore))
      end if

      allocate (xold(ni), yold(ni), zold(ni))

      do n = 7, 10  ! Aft-body patches

         do j = 1, nj

            xold(:) = new_grid(n)%x(:,j,1)
            yold(:) = new_grid(n)%y(:,j,1)
            zold(:) = new_grid(n)%z(:,j,1)

!           Just shift the shoulder region:

            dx = new_grid(m)%x(j,mj,1) - new_grid(n)%x(ni,j,1)
            dy = new_grid(m)%y(j,mj,1) - new_grid(n)%y(ni,j,1)
            dz = new_grid(m)%z(j,mj,1) - new_grid(n)%z(ni,j,1)
            new_grid(n)%x(iu2:ni,j,1)  = new_grid(n)%x(iu2:ni,j,1) + dx
            new_grid(n)%y(iu2:ni,j,1)  = new_grid(n)%y(iu2:ni,j,1) + dy
            new_grid(n)%z(iu2:ni,j,1)  = new_grid(n)%z(iu2:ni,j,1) + dz

!           Morph the rest of the spoke not on the aft center body:

            call nuline3d (iu1, iu2, xold, yold, zold, &
                           new_grid(n)%x(:,j,1), new_grid(n)%y(:,j,1), &
                           new_grid(n)%z(:,j,1))

!           Option to seek ~uniform fabric/rib thickness:

            if (fix_thickness) then

               do i = 1, njfore  ! Set up forebody spoke for interpolation
                  xfore(i) = new_grid(m)%x(j,i,1)
                  rfore(i) = sqrt (new_grid(m)%y(j,i,1)**2 + &
                                   new_grid(m)%z(j,i,1)**2)
               end do

!              At iu1, a wide aft body may be outboard of frustum_radius.
               i = iu2  ! Calculate last dx for this spoke = fabric thickness
               raft = sqrt (new_grid(n)%y(i,j,1)**2 + new_grid(n)%z(i,j,1)**2)
               new  = true
               call lcsfit (njfore, rfore, xfore, new, method, 1, raft, &
                            xinterp, unused)
               dx = new_grid(n)%x(i,j,1) - xinterp

               do i = iu2, iu1, -1  ! So we can use raft at iu1
                  raft = sqrt (new_grid(n)%y(i,j,1)**2 + &
                               new_grid(n)%z(i,j,1)**2)
                  call lcsfit (njfore, rfore, xfore, new, method, 1, raft, &
                               xinterp, unused)
                  new = false
                  new_grid(n)%x(i,j,1) = xinterp + dx
               end do

!              If the forebody faceting starts inboard of a fat aft body (iu1),
!              we need to blend from some iu0 to iu1-1 to account for a nonzero
!              shift at iu1:

               if (raft > frustum_radius) then
                  call nuline3d (iu0, iu1, xold, yold, zold, &
                                 new_grid(n)%x(:,j,1), new_grid(n)%y(:,j,1), &
                                 new_grid(n)%z(:,j,1))
               end if

            end if

         end do

         m = m - 1
      end do

      deallocate (xold, yold, zold)
      if (fix_thickness) deallocate (xfore, rfore)

!     Recover the spoked form of the faceted surfaces for volume gridding:

      ni = new_grid(3)%ni  ! Regridded nose patch azimuthal count ...
      nj = new_grid(3)%nj  ! ... and spokewise count
      ni_regrid = ni_regrid_forebody

      j = 1  ! Spoke counter
      do n = 3, 6
         do i = 1, ni
            spoke_fore%x(ni_regrid:numi_forebody,j,1) = new_grid(n)%x(i,:,1)
            spoke_fore%y(ni_regrid:numi_forebody,j,1) = new_grid(n)%y(i,:,1)
            spoke_fore%z(ni_regrid:numi_forebody,j,1) = new_grid(n)%z(i,:,1)
            j = j + 1
         end do
         j = j - 1
      end do

      open  (lunspoke, file='spoked_surface_fore.g', status='unknown')
      write (lunspoke, '(i1)') 1
      write (lunspoke, '(2i4, i2)') numi_forebody, numj, 1
      write (lunspoke, '(8es15.7)') spoke_fore%x, spoke_fore%y, spoke_fore%z
      close (lunspoke)

      nj = new_grid(7)%nj  ! Aft patch azimuthal point count
      ni_regrid = ni_regrid_aft_body

      j = 1
      do n = 10, 7, -1
         do i = 1, nj
            spoke_aft%x(ni_regrid:numi_aft_body,j,1) = new_grid(n)%x(:,i,1)
            spoke_aft%y(ni_regrid:numi_aft_body,j,1) = new_grid(n)%y(:,i,1)
            spoke_aft%z(ni_regrid:numi_aft_body,j,1) = new_grid(n)%z(:,i,1)
            j = j + 1
         end do
         j = j - 1
      end do

      open  (lunspoke, file='spoked_surface_aft.g', status='unknown')
      write (lunspoke, '(i1)') 1
      write (lunspoke, '(2i4, i2)') numi_aft_body, numj, 1
      write (lunspoke, '(8es15.7)') spoke_aft%x, spoke_aft%y, spoke_aft%z
      close (lunspoke)
      go to 95

  90  write (luncrt, '(/, (a))') &
        'Suppressing the aft body morph to umbrella shape.', &
        'Plot the present surface grid and pick aft spoke indices that ', &
        'define the (essentially) straight portions inboard of the shoulder.', &
        ' '

  95  deallocate (spoke_fore%x, spoke_fore%y, spoke_fore%z)
      deallocate (spoke_aft%x,  spoke_aft%y,  spoke_aft%z)

      end subroutine morph_aft_to_umbrella

   end program capsule_grid
