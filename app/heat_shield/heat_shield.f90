
!!!!!!                    Program HEAT_SHIELD modules                     !!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      MODULE AERO_MOD                         ! Aerodynamic coefficients, etc. !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      REAL :: &
         Alpha, Beta, M_inf, Gamma, S_REF, REF_AREA, REF_LEN, XYZ_CG(3),       &
         XYZ_CG0(3), CD, CS, CL, ML, MM, MN, SWET, SBASE

      END MODULE AERO_MOD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      MODULE COEFFS_MOD            ! Initial & target aerodynamic coefficients !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      REAL :: &
         LD, DCMDA, HEAT_FLUX, TEMP_MAX, CM_TOL,                               &
         AL_TARG, CD_TARG, CM_TARG, CR_TARG, CRVG_TARG, CRVJ_TARG, TM_TARG,    &
         AL_INI, MA_INI, CL_INI, CD_INI, CM_INI, LD_INI,                       &
         TEMP_MAX_INI, HEAT_FLUX_INI, SECTIONAL_AREA_INI, SWET_INI,            &
         V_EFF_INI, VOLUME_INI

      END MODULE COEFFS_MOD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      MODULE CONST_MOD                                             ! Constants !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      INTEGER, PARAMETER :: &
         LWRIT   = 1,       & ! Printable outputs; NPOPT writes to standard out
         LNPOPT  = 2,       & ! Temporary file used by NPOPT
         LGEOM   = 3,       & ! Input geometry sections
         LGSAVE  = 4,       & ! Optimized geometry sections (same format)
         LREAD   = 5,       & ! Input control file; NPOPT uses unit 6
         LXYZ    = 7,       & ! Regularized geometry (PLOT3D surface, /mgrid)
         LKEEL   = 8,       & ! Keel-crown line for possible HB_GRID (2D) use
         LTRIANG = 9          ! Corresp. triangulation (Tecplot) with centroid
                              ! pressure coefficients
      LOGICAL, PARAMETER :: &
         FALSE   = .FALSE., &
         TRUE    = .TRUE.

      REAL :: &
         DEG2RAD, EPSMCH, PI, RAD2DEG, RPI, SMALL, BIGBND, TWOPI, TWOPISQ,     &
         SIXROOTPI

      END MODULE CONST_MOD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      MODULE DESIGN_VARS_MOD                         ! Design parameterization !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      INTEGER :: &
         IALPHA, NBUMPS, NPLANF

      INTEGER, ALLOCATABLE, DIMENSION (:) :: &
         ILCON

      REAL, ALLOCATABLE, DIMENSION (:) :: &
         V, VSCALE, IGROUP, TLCON,        &
         ZBUMP, ZBMIN, ZBMAX, ZBEXP,      &
         RBUMP, RBMIN, RBMAX, RBEXP

      LOGICAL :: &
         CALC_X_RESIDUAL

      CHARACTER, ALLOCATABLE, DIMENSION (:) :: &
         RBTYP * 6, ZBTYP * 4

      END MODULE DESIGN_VARS_MOD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      MODULE GEOM_MOD                ! Initial and current geometry definition !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      INTEGER :: &
         NI_SEC, NO_SEC, NT_RAD, NT_SEC, FULL, NON_DIM

      REAL :: &
         XICNTR, YICNTR, ZICNTR, XIRING, RIRING,  TOT_AZ, CIRCLE, &
         XWCNTR, YWCNTR, ZWCNTR, XWRING, RWRING, RTOT_AZ, CUSP

      REAL :: &  ! Sphere-cone parameters (inputs, not derived parameters)
         X_NOSE, R_NOSE, RADIUS_NOSE, RADIUS_BASE, RADIUS_SHOULDER, &
         HALF_CONE_ANGLE, SKIRT_ANGLE, SKIRT_LENGTH

      REAL :: &  ! Biconic parameters (additional inputs + some derived)
         HALF_CONE_ANGLE_FORE, HALF_CONE_ANGLE_AFT, RADIUS_CONE_JUNCTURE, &
         X_CONE_JUNCTURE, R_CONE_JUNCTURE, CONE_LENGTH_FORE, CONE_LENGTH_AFT

      LOGICAL :: &
         AXISYMMETRIC, SMOOTH_X, SPHERECONE, TWO_PARTS

      INTEGER, ALLOCATABLE, DIMENSION (:) :: &
         NI_INIT, NO_INIT

      REAL, ALLOCATABLE, DIMENSION (:) :: &
         ZI_INIT, CI_INIT, ZO_INIT, CO_INIT, VO_INIT, &
         ZI_WORK, CI_WORK, ZO_WORK, CO_WORK, VO_WORK, &
         RI1,     XI1,     RO1,     XO1,              &
         COS_ZI,  COS_ZO,  SIN_ZI,  SIN_ZO

      REAL, ALLOCATABLE, DIMENSION (:,:) ::  &
         RI_INIT, XI_INIT, RO_INIT, XO_INIT, &
         RI_WORK, XI_WORK, RO_WORK, XO_WORK, YO_PRTB, ZO_PRTB

      END MODULE GEOM_MOD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      MODULE GRID_MOD                                ! Surface grid definition !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      INTEGER :: &
         N_INNER, N_OUTER, N_RAD_MAX, N_AZM_MAX

      REAL :: &
         SECTIONAL_AREA, TVD_GC, VOLUME, V_EFF,  X_RESIDUAL, &
         RI_SPL, RI_SPQ, RI_SPS, RI_SPC, &
         RO_SPL, RO_SPQ, RO_SPS, RO_SPC, &
         CRVI_MAX, CRVI_TOT, &
         CRVJ_MAX, CRVJ_TOT, &
         CRVG_MAX, CRVG_TOT

      REAL, ALLOCATABLE, DIMENSION (:) :: &
         AZM, COSPHI, SINPHI

      REAL, ALLOCATABLE, DIMENSION (:,:) :: &
         RFIN, XFIN, YFIN, ZFIN, XLEN,      &
         GAUSS_CURVATURE, MEAN_CURVATURE

      END MODULE GRID_MOD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      MODULE OPT_MOD                                  ! Optimization variables !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      INTEGER :: &
         NITMAX,  NDITER,  NFTOTL,  NDV,     NCLIN,   NCNLN,    NITER,         &
         NROWA,   NROWJ,   NROWR,   LENIW,   LENW,    NBNDS,    INFORM,        &
         MAJITS,  MAJORPL, MINORPL, MINORIL, NPRFREQ, LEVELVER,                &
         ALPHA_MODE, GRAD_MODE

      INTEGER, ALLOCATABLE, DIMENSION (:) :: &
         ISTATE,  IWVEC,   INLCON,  JNLCON

      REAL :: &
         ETA,             STEPLIM,      EPSOBJ,      EPSM6TH,     OBJINI,      &
         OBJ_CONFUN,      PENPARAM,     TOLLIN,      TOLNLIN,     TOLOPT,      &
         RHO_AL,  RHO_CD, RHO_LD,       RHO_CD_INV,  RHO_CD_TARG, RHO_CM_MAX,  &
         RHO_CM_TARG,     RHO_CURVI,    RHO_CURVJ,   RHO_CVIMX,   RHO_CVJMX,   &
         RHO_SMTH_X,      RHO_TEMP_MAX, RHO_SWET,    RHO_VNORM,   RHO_AREA,    &
         RHO_VOLUME,      RHO_V_EFF,    RHO_CD_A,    RHO_CM_ALF,  RHO_TVD_GC,  &
         FRACTION_CRVIMX

      REAL(KIND=4) :: CPU0

      REAL, ALLOCATABLE, DIMENSION (:) :: &
         FFWD,    GRAD,    CVEC,    WVEC,    AITCH,   ONEOVERH, BL,      BU,   &
         C,       CJACBAK, CLAMDA,  CONLIN,  CPRINT,  BLLIN,    BULIN,         &
         XNLCON,  SNLCON

      REAL, ALLOCATABLE, DIMENSION (:,:) :: &
         AMAT, RMAT, CJAC

      LOGICAL :: &
         CONFUN_DOES_OBJ, FIXED_CG, NEED_DCMDA

      CHARACTER, ALLOCATABLE, DIMENSION (:) :: &
         LCTYPE*6, NLCTYPE*6

      END MODULE OPT_MOD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      MODULE TRIANG_MOD                                 ! Triangulated surface !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      INTEGER :: &
         N_Triangles, N_verts

      INTEGER, ALLOCATABLE, DIMENSION (:) ::   &
         ICOMP

      INTEGER, ALLOCATABLE, DIMENSION (:,:) :: &
         IPT_TRI

      REAL, ALLOCATABLE, DIMENSION (:) ::      &
         CP, FACE_AREA, SIGNS

      REAL, ALLOCATABLE, DIMENSION (:,:) ::    &
         XYZ, XYZ_CENTROID, XYZ_NORM

      END MODULE TRIANG_MOD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
      PROGRAM HEAT_SHIELD
!
!        HEAT_SHIELD performs shape optimization for space vehicle forebodies
!     using Newtonian-type aerodynamics, which are appropriate for hypersonic
!     flight only.
!
!        Provision is made for an optional annulus (cruise engine mount ring)
!     which must remain fixed in shape - the inner & outer shield are optimized
!     around the ring, which severely limits the range of motion and may be
!     impractical to treat well here in practice.  (Think of the curvature
!     discontinuities that could develop across the ring.)
!
!        If no ring is present (input TWO_PARTS = F), the defining sections
!     (spokes) are treated as "outer" sections, with no "inner" sections.
!
!        Either half the shield or the whole shield may be optimized, although
!     it is unlikely that the full shield case is ever the right choice: only
!     the "half" case ensures bilateral symmetry.
!
!        This version includes an option to retain axisymmetry throughout.  It
!     is invoked for the single-part case with a single defining section in
!     the input (and optimized) geometry file.
!
!     COORDINATE SYSTEM:
!
!        Right-handed, with X down-stream, Y up, and Z > 0 for port side.
!
!
!     CONTROL FILE:
!
!        Any file name, such as heat_shield.inp, specified on the command line.
!
!
!     PRIMARY CONTROL INPUTS:
!
!        NDV        Number of design variables
!        NCLIN      Number of linear constraints (NPOPT only)
!        NCNLN      Number of nonlinear constraints (NPOPT only)
!        NITMAX     Number of optimization iterations >= 0
!        GRAD_MODE  Gradient mode (in application, not in optimization pkg.):
!                      2 = 2-point forward differences with given h;
!                      3 = 3-point central differences with h * EPSOBJ**(-1/6)
!
!
!     OPTIMIZER INPUTS:
!
!        ETA        Controls the line search's acceptance of a
!                   sufficiently lower objective. 0 < ETA < 1.0;
!                   try 0.2 for expensive objectives
!
!        EPSOBJ     Minimum absolute value of a significant
!                   difference in OBJ
!
!        N_RAD_MAX  # regularized pts., radial (N_INNER + N_OUTER - 1)
!        N_AZM_MAX  # regularized pts. in azimuthal direction
!        N_INNER    # regularized pts., inner
!        N_OUTER    # regularized pts., outer
!        CUSP       Fraction (applied to L_ref) used to determine the
!                   points fudged at the apex to smooth out any cusps
!        TWO_PARTS  T means ring is present, else N_RAD_MAX = N_OUTER
!        SMOOTH_X   T means fit quartics to X vs. azimuth to smooth
!                   likely irregularities near phi = zero; this is
!                   on the perturbed geometry points, not the surface
!                   grid, so that smoothed geometry is saved; it means
!                   the spokes should all be defined with consistent
!                   point distributions.
!
!        RI_SPL, RI_SPQ, RI_SPS, RI_SPC  Inner and outer controls for
!        RO_SPL, RO_SPQ, RO_SPS, RO_SPC  FOILGRD on radial sections
!                   Long-time usage:     0.5 0 0 0.5
!                   December 2010 recommendations:
!                   (1) R*_SPL < 0 means reverse the distribution, as
!                       found helpful with [-]0.04, 0, 0.3, 0.66 inputs;
!                   (2) R*_SPL = 1 means use curvature-based spacing;
!                       R*_SPQ = 0.5 (curvature-effect power) is recommended;
!                       R*_SPS = 1 = R*_SPC is recommended for sphere-cones
!                       (smoothing of the curvature-based shape function and
!                       index-based smoothing of the resulting distribution)
!                       else enter 0 to suppress either type of smoothing.
!
!
!     OBJECTIVE FUNCTION INPUTS:
!
!        The "RHO"s are multipliers for the terms included in the objective
!        being minimized.  Most have nonlinear constraint analogues.  Normally,
!        constraints are preferable and picking one form or the other is
!        advisable (not both), but this is a gray area.  For instance, RHO_AL
!        appears to assist satisfaction of the ALPHA constraint.  (Treating
!        Alpha as a variable is better yet.)
!        Normally, RHO_* > 0. means optimization should improve the situation.
!        For instance, L/D is normally increased if RHO_LD > 0., although some
!        constraint(s) may actually require it to decrease.  In other cases,
!        a negative RHO may actually be appropriate (as in the case of trying
!        to INCREASE CD instead of decrease it).
!
!        RHO_CM_TARG Pitching moment multiplier
!        RHO_CM_MAX (1. - |Pitching moment|) ...
!        RHO_CD     Drag ...
!        RHO_LD     -L/D (better behaved than D/L)
!        RHO_TEMP_MAX Max. temperature ...
!        RHO_CD_INV Inverse drag ...
!        RHO_CD_TARG Target drag ...
!        RHO_AL     Square of MAX (AL_TARG - Alpha, 0.)
!        RHO_CURVI  Sum of MAX (radial CURV(i,j) - CR_TARG, 0.) ** 2
!        RHO_CURVJ  Sum of MAX (azimuthal CURV(i,j) - zero, 0.) ** 2 ! DON'T USE
!        RHO_CVIMX  Max. curvature in radial direction over all spokes
!        RHO_CVJMX  Max. curvature in azimuthal direction (probably worthless)
!        RHO_SMTH_X Sum of residuals**2  when quartics are fit to X vs. azimuth
!        RHO_SWET   > 0 means penalize SWET > SBASE;
!                   < 0 means penalize SWET < SBASE; for a violation,
!                   RHO_SWET is multiplied by (SWET - SBASE)/SREF
!        RHO_VNORM  Squared 2-norm of the sine bump variables; small multipliers
!                   should tend to regularize the solution by preferring the
!                   shortest-length set of sine bump multipliers; try 1.E-6
!        RHO_AREA   -cross-sectional area (rho > 0 maximizes the area)
!        RHO_VOLUME -volume: rho > 0 maximizes volume
!        RHO_V_EFF  -eff. vol. coef. [6 sqrt (pi) volume / Sarea^1.5]; rho > 0
!        RHO_CD_A   -(CD x cross-sectional area): rho > 0 maximizes CD x A
!        RHO_CM_ALF dCM/dAlpha
!        RHO_TVD_GC Total variation diminishing function of the Gaussian
!                   curvature distribution G(i,k) along the surface grid spokes;
!                   use some small multiplier (0.001, say) to help smooth the
!                   curvature (and hence the geometry) in the i index direction;
!                      obj. term = sum over k and i of (G(i+1,k) - G(i,k))**2
!                   Note that the variation in the k (azimuthal) direction is
!                   not included, because the grid lines tend not to follow the
!                   shoulder for asymmetric cases.
!
!     AERODYNAMIC INPUTS:
!
!        CM_TARG    Target pitching moment (needs ALPHA_MODE > 0)
!        CM_TOL     Tolerance for matching target pitching moment
!        CD_TARG    Target CD
!        AL_TARG    Target Alpha: angles less than AL_TARG are penalized
!        CR_TARG    Maximum radial curvature at any point;
!                   curvature >= 0. for a normal convex shape, but CR_TARG > 0.
!                   avoids flatness; used by obj. fn. and nonl. constraint forms
!        TM_TARG    Target maximum temperature
!
!        ALPHA_MODE 0 = Fixed Alpha;
!                   1 = Fixed pitching moment (CM_TARG);
!                   2 = Alpha is an optimization variable ('ALPHA ' type)
!                       in which case the starting guess and bounds are
!                       read along with the other variables and the following
!                       input Alpha is ignored
!        Alpha      [Initial] angle of attack, degrees (unless ALPHA_MODE = 2)
!        Beta       Yaw angle, deg
!        M_inf      Free stream Mach number
!        Gamma      Estimate of ratio of specific heats (e.g. 1.2)
!
!
!     DESIGN VARIABLE INPUTS:
!
!     The bulk of the design variables are Hicks-Henne-type "sine bumps."
!     More precisely, they are multipliers of such perturbing shape functions.
!     These are preceded by any "planform"-type variables (XCEN, YCEN, etc.),
!     while the optional ALPHA variable follows the sine bumps.  These sine
!     bumps (which have to be decayed in both surface directions) influence
!     the multi-column form of all of the variable inputs.  If variables are
!     entered out of order, the program will detect this and stop with a
!     diagnostic.  Descriptions for the variables and related inputs follow.
!
!     #        Ordinal number of design variable - not actually used,
!              so it doesn't really matter if they get out of order
!     RBTYP    6-character design variable type for functions that
!              operate on the planform or in the radial direction.
!              Available variable names:
!
!        XCEN     The X axis is taken to be the longitudinal axis of
!                 the vehicle.  XCEN changes the location of the
!                 apex, or forward-most point of the  heat shield.
!                 It is additive:  XCEN = 0. means no perturbation.
!        YCEN     Controls the "vertical" position of the apex.
!        ZCEN     Controls the "sideways" position of the apex (assumed
!                 fixed for now).
!        X_CG     Applied as perturbations to the CG coordinates read from
!        Y_CG     the input geometry file
!        XRNG     Controls the X position of the cruise engine mount ring
!        ROUT     Controls the radial distance to the outer edge of the
!                 shield; note that this may be unsymmetric in azimuth.
!        XOUT     Controls the X position of the outer edge of the shield;
!                 note that this may be unsymmetric in azimuth.
!        TOUT     Controls the position of the outer edge of the shield
!                 along the line joining it to the apex, as long as the
!                 input value of "RBTYPE" equals XOUT/ROUT = tangent of
!                 (90 - cone half angle) (but can this change during the
!                 optimization?).
!                 (The idea is to be able to grow a trim tab.)
!
!                 Most of the following Hicks-Henne-type perturbing
!                 functions are suited to airfoils, not heat shields.
!
!        SIN      Standard (modified) sine perturbing function
!        SINC     Cyclic sine function suited to lofting heat shield
!                 perturbations in the azimuthal direction
!        SINF     Flipped sine function
!        SIN1     Symmetric sine function
!        SIN2     Symmetric flipped sine function
!        SIN3     Ensures airfoil LE/TE symmetry; use
!                 [XA, XB] = [0, 1] & XC in (0, .5)
!        SIN4     Ensures airfoil LE/TE symmetry; use
!                 [XA, XB] = [0, .5] & XC in (0, 1)
!        COSL     Half [co]sine, peak at left  (LE)
!        COSR     Half [co]sine, peak at right (TE)
!        LCOS     Inverted form of COSL
!        RCOS     Inverted form of COSR
!        EXP      Standard Exponential function
!        LED      Leading edge droop function  (also LEAD)
!        TRL      Trailing edge droop function (also TRAIL)
!        WAG      Wagner shape function
!
!        ALPHA    If ALPHA_MODE = 2, this variable must follow the sine bumps.
!                 It allows the optimization to proceed such that the specified
!                 pitching moment is achieved at convergence but not necessarily
!                 before that.  This is more effective than the alternative of
!                 iterating on angle of attack for every geometry perturbation.
!                 Enter a value and scale factor such that their product is the
!                 desired starting guess for the angle of attack in degrees.
!                 Likewise for the differencing interval (~one-millionth of a
!                 degree) and the lower & upper bounds.  E.g., for -20 degrees:
!
!                 V = -2.00   VSCALE = 10.0   H = 1.E-5   BL = -2.2   BU = 0.0
!
!     ZBTYP    4-character design variable type for functions that operate
!              azimuthally; see SIN2, SINC, etc.
!
!     REXP     Exponent used by chordwise perturbing shape function
!     RBMP     Center location (as a fraction in [0, 1]) for shape function
!     RMIN     Normally 0.; may be > 0. to avoid affecting forward part of spoke
!     RMAX     Normally 1.; may be < 1. to avoid affecting aft part of spoke
!     ZEXP     Exponent used by azimuthal (lofting) shape function
!     ZBMP     Center location of lofting shape function
!     ZMIN     As for RMIN, RMAX but in aximuthal lofting direction
!     ZMAX
!     V        Design variable value as seen by the optimizer
!     VSCALE   Scale factor for the design variable to keep it ~O(1);
!              V * VSCALE is typically the shape function multiplier
!     AITCH    Step size used for estimating the gradient of the objective
!              w.r.t. the variable by forward differencing (GRAD_MODE = 2);
!              the actual perturbation used is AITCH * VSCALE;
!              AITCH is also used for forward derivatives of the nonlinear
!              constraints; see ALPHA_MODE = 3 above for central differencing
!     BL       Lower bound on the design variable as seen by the optimizer;
!              BL = -999. means the variable has no lower bound, else
!              BL scaling should match that of the variable - e.g., for a
!              variable with VSCALE = 0.1, BL = -10. means the effective
!              values are >= -10.*0.1 = -1., consistent with AITCH usage
!     BU       Upper bound on the design variable as seen by the optimizer;
!              BU = 999. means the variable has no upper bound, else BU
!              scaling should match that of the variable
!
!
!     LINEAR CONSTRAINTS:
!
!     If NCLIN = 0, just include two header lines.
!
!        #       Ordinal number of linear constraint - not used
!        LCTYPE  6-character linear constraint type:
!
!                EDGE_K  <Explain>
!                EDGE_Z  <Explain>
!
!        BL      Lower bound for linear constraint; -999. = -BIGBND
!        BU      Upper bound for linear constraint;  999. =  BIGBND
!
!
!     NONLINEAR CONSTRAINTS:
!
!     If NCNLN = 0, just include two header lines.
!
!        #          Ordinal number of nonlinear constraint - not used
!        NLCTYPE    6-character nonlinear constraint type:
!
!           H_FLUX  Inactive
!           W_TEMP  Inactive
!           CM      |Pitching moment|
!           CD      Drag coef.
!           LD      L/D
!           CD_A    CD x cross-sectional area of forebody
!           CM_ALF  Derivative of CM w.r.t. Alpha
!           CG_OFF  CG offset (TBD)
!           ALPHA   Angle of attack (don't use if ALPHA_MODE = 2)
!           AREA    Cross-sectional area at forebody base
!           VOLUME  Forebody volume
!           V_EFF   Effective volume coef. = 6 sqrt (pi) volume / Sarea^1.5
!           WAREA   (Wetted area - SBASE) / REF_AREA; use BU = 0.
!           CURVI   Sum of squares of MAX (curvature - CR_TARG, 0.)
!                   in radial direction at all geometry points in (r,x) space
!           CURVJ   Sum of squares of MAX (curvature - CRVJ_TARG, 0.)
!                   in azimuthal direction at all geometry pts. in (y,z) space;
!                   enter CRVJ_TARG >= 0. via XNLCON()
!           CURVG   Sum of squares of MAX (Gaussian curvature - CRVG_TARG, 0.)
!                   over the computational surface grid;
!                   enter CRVG_TARG >= 0. via XNLCON()
!           CRVIMX  Maximum (peak) curvature in the I direction, not counting
!                   the fraction beyond NI * (XNLCON() in [0., 1.])
!           CRVGMX  Maximum Gaussian curvature on the computational grid
!           SMTH_X  Constraint form of the SMOOTH_X option: penalize the
!                   sum of squared deviations from best-fit quartics for
!                   X vs. azimuth at each I of the defining sections;
!                   use BL 0., BU = +eps (TBD).
!
!        BL         Lower bound on nonlinear constraint value;
!                   -999. = -BIGBND
!        BU         Upper bound on nonlinear constraint value;
!                   +999. =  BIGBND
!                   N.B.:  these bounds should be in the same units
!                   as the quantities being bounded;
!                   SNLCON will be applied by this program to the
!                   given BL and BU; this is in contrast to the
!                   bounds on the variables themselves, which
!                   (as for their finite differencing intervals) refer
!                   to the scaled variables
!        XNLCON     Real quantity needed (or not) by the nonlinear
!                   constraint
!        INLCON     First  index needed (or not) by the nonlinear
!                   constraint
!        JNLCON     Second index needed (or not) by the nonlinear
!                   constraint
!        SNLCON     Scale factor used to provide the optimizer with
!                   nonlinear constraint gradients of order 1
!
!
!     OPTIONAL INPUTS FOR NPOPT:
!
!     NPOPT's optional inputs.  It should appear after the nonlinear
!     constraint inputs.   Anything beyond that is ignored.
!     See the NPSOL User Guide for the control parameters in the
!     following namelist.  (Some of them apply to SNOPT only.)
!
!     NAMELIST /NPOPTIONS/ &
!        LEVELVER, MAJORPL, MINORPL, MINORIL, NPRFREQ, &
!        PENPARAM, STEPLIM, TOLLIN, TOLNLIN, TOLOPT
!
!        STEPLIM    Limits stepsize of line search; default = 2.
!
!
!     INPUT/OUTPUT FILES:
!
!        See CONST_MOD module above and OPENs in main program and in SUMMARY.
!
!
!     GEOMETRY INPUT FORMAT (heat_shield.geo):
!
!        Title, e.g. Geometry for Mars 2005 lander
!             S_ref     L_ref   X_cg(1)   X_cg(2)   X_cg(3)   S_base
!         11.044662      3.75      0.81      0.09      0.00  12.2552
!              FULL    NI_SEC    NO_SEC   NON-DIM
!                 0        11        11         1
!            XICNTR    YICNTR    ZICNTR    XIRING    RIRING
!               0.0       0.0       0.0   0.24419   0.83179
!        OUTER RADIAL LINES
!        Section 1 Outer
!           NO_INIT   ZO_INIT   CO_INIT   VO_INIT
!                49       0.0   1.87498   0.69098
!           RO_INIT   XO_INIT
!           0.00000   0.00000
!           0.03980   0.00087
!            :         :
!           1.87433   0.67667
!           1.87498   0.69098
!        Section 2 Outer
!           NO_INIT   ZO_INIT   CO_INIT   VO_INIT
!                49       0.1   1.87498   0.69098
!           RO_INIT   XO_INIT
!           0.00000   0.00000
!           0.03980   0.00087
!            :         :
!
!     ANALYTIC OPTIONS FOR COMMON GEOMETRIES (heat_shield.analytic):
!
!        Two choices are provided via the same file, and the number of lines
!        in the file is used to distinguish them (10 and 12 lines respectively):
!
!        SPHERE-CONE (including Apollo-type spherical section) OPTION:
!
!        Analytic specification of a symmetric heat_shield geometry:
!        0.       ! x_nose
!        0.       ! r_nose
!        5.       ! radius_nose
!        10.      ! radius_base
!        0.25647  ! radius_shoulder
!        65.      ! half_cone_angle; 0 => spherical section (Apollo-like)
!        0.       ! skirt_angle
!        0.       ! skirt_length (from aft-shoulder tangency pt. x to cut-off x)
!
!        BICONIC OPTION:
!
!        Title
!        0.       ! x_nose
!        0.       ! r_nose
!        5.       ! radius_nose
!        6.       ! radius_cone_juncture
!        10.      ! radius_base
!        0.25647  ! radius_shoulder
!        70.      ! half_cone_angle_fore
!        55.      ! half_cone_angle_aft
!        35.      ! skirt_angle
!        0.01     ! skirt_length (from aft-shoulder tangency pt. x to cut-off x)
!
!     GEOMETRY NOTES:
!
!        A shield consists of either one or two parts.  A single part
!        is considered all "outer".  Two parts are separated by a
!        circular engine mount ring which must not change shape.
!
!        Defining sections are entered in cylindrical coordinates
!        (x, r, theta) where x points into the shield along its axis
!        of symmetry (if it is conical), r is the radial distance from
!        the x axis, and theta ("z" in the code) is entered as the
!        fraction (0. - 1.) of 180 or 360 degrees starting from the
!        12 o'clock position, depending on the FULL input.
!
!        The option to smooth X vs. azimuth at the geometry level
!        requires that all input sections have the same number of
!        points, distributed in consistent fashion.
!
!        S_ref      Reference area for the whole shield.
!                   If FULL = 0, half of S_ref is used internally.
!
!        S_base     Wetted area, used for "WAREA" constraint, q.v.
!
!        FULL       1 means the whole shield is being optimized;
!                   0 means the "right" half looking downstream.
!                   FULL = 1 is inadvisable during optimization.
!
!        NON_DIM    Controls normalization/shearing  of input sections
!                   upon reading them (RD_GEOM):  1 means normalize.
!                   If the shield is in two parts, the inner sections
!                   are normalized by the ring radius, while for outer
!                   sections, the difference between r and ring radius
!                   is normalized by the shield radius. (??)
!
!         XICNTR,   Coordinates of shield vertex (which may move).
!         YICNTR,
!         ZICNTR
!
!         XIRING,   (x,r) coordinates of points on the ring
!         RIRING    (not used for the 1-part case).
!
!         If heat_shield.analytic is present, a symmetric sphere-cone,
!         spherical-section or biconic forebody is generated, and S_ref and
!         S_base are derived from the alternative geometry inputs.  This is
!         intended to simplify surface gridding of standard geometries for
!         CFD purposes (with NITMAX = 0).
!
!         Note that an analytic generatrix is discretized almost uniformly
!         along the arc, with the key tangency points captured precisely,
!         then it is treated as though it had been read from heat_shield.geo.
!         The FOREBODY_REGRID program should still be applied to the resulting
!         heat_shield.xyz output from HEAT_SHIELD to remove the singular point,
!         then a Gridgen glf file can build an initial hyperbolic volume grid.
!
!     HISTORY:
!
!        Aug. 2000  J.Reuther   Initial NEW_AERO_OPT implementation,
!                               starting with Traj_opt/SYN87-SB, and
!                               applied to the Mars-05 shield.
!        Nov. 2000  D.Saunders  Plugged in argument-driven routines
!                               from program NEWTONIAN (let's keep
!                               them in common), and polished the
!                               source somewhat; 3-pt. curvatures.
!                               Renamed it HEAT_SHIELD.
!        Dec. 2000      "       Aerodynamics now use FACE_AREA(*),
!                               not PROJECTED_AREA(1:3,*). Geometry
!                               and Cp distribution are now saved.
!        May  2001  J.Reuther   Added wetted area constraint.
!        Nov. 2001  D.Saunders  Sections are denormalized/unsheared
!                               before regularization; NOCUSP option;
!                               minimize -L/D, not D/L; move the
!                               moment center with the apex; gridding
!                               in the azimuthal direction needs to be
!                               nonlinear if the apex can move away
!                               from the center...
!        01/28/02       "       ... or for any significant variations
!                               in the defining sections (3-D effects).
!                               Applying interior corrections to the
!                               linear result based on splines at the
!                               rim typically misses an X correction.
!                               Therefore, spline every grid I.
!        01/29/02       "       Added CURVJ constraint: r and x vs.
!                               azimuthal angle, however, doesn't do
!                               it.  Use y vs. x approximation.
!        01/31/02       "       Save the perturbed geometry in PLOT3D
!                               form, assuming the spokes have common
!                               point counts.  SUMMARY writes file
!                               'heat_shield.spokes' to unit LXYZ.
!                               Added simple-minded CREASE constraint.
!        02/04/02       "       Option to smooth X vs. azimuth:  do it
!                               on the geometry data, not the grid, so
!                               that the saved geometry is smoothed.
!        02/07/02       "       SMOOTH_X option is promising, but it can
!                               go unstable.  Retain it (possibly just
!                               for use with NITMAX = 0), but adapt it to
!                               penalize X deviations from quartics via
!                               RHO_SMTH_X > 0 or new SMTH_X constraint.
!        02/11/02       "       1: Apex-smoothing quartics are now tangent at
!                               a precise R from the apex, not at the index
!                               nearest and inboard of R = CUSP * Lref.  This
!                               should eliminate occasional bad gradients.
!                               2: ALPHA_MODE = 2 allows Alpha to be one of the
!                               optimization variables.  This should reduce
!                               the nonlinearity:  we just need CM = 0 at the
!                               solution, not at every evaluation along the way.
!        02/13/02       "       RHO_VNORM option encourages the shortest-length
!                               solution for the non-unique sine bump variables.
!                               Switched from local splines to conventional
!                               splines in NOCUSP, which was presumably the
!                               cause of occasional bad gradients w.r.t. YCEN.
!                               Removed simple-minded CREASE constraint.
!        03/27/02       "       Two-part case needed changes in RD_GEOM and
!                               PERTURB.
!        11/08/02       "       Save (x,y,z,Cp) file for structural analysis.
!        11/15/02       "       NOCUSP failed if r wasn't monotonic (as at the
!                               rim).  Therefore, spline inner halves only.
!        05/21/07       "       Cleaned out unconstrained option (QNMDIF2) in
!                               preparation for revival with conic section
!                               geometry options.
!        08/06/07       "       Converted to *.f90 format (80-column free form).
!        08/07/07       "       Added central differencing option (GRAD_MODE).
!        08/08/07       "       Added volume and cross-sectional area calcs.
!                               and constraints on these and on CD x A and on
!                               dCM/dAlpha; these quantities may also contribute
!                               to the objective.
!        08/09/07       "       Added effective volume coef. constraint and
!                               contribution to the objective.  Note that the
!                               sectional area is included in the surface area.
!        08/10/07       "       Added X_CG and Y_CG optimization variables,
!                               which behave as for XCEN and YCEN (added to the
!                               initial values found in the input geometry).
!        08/24/07       "       Since the surface grid is left-handed in the
!                               Jameson coordinate system, save it in right-
!                               handed form at the end of a run (but leave it
!                               left-handed during intermediate saves for
!                               I/O efficiency reasons).  Right-handedness is
!                               desirable for HYPER_AERO compatibility.
!                               RHO_CVIMX and -CVJMX options added to objective
!                               function.  Ignore shoulder region for _CVIMX.
!                               Confine _CVIMAX to lee side also.
!        11/25/07       "       Added an axisymmetric option (single defining
!                               section, single-part case).
!        11/27/07       "       Added CRVIMX constraint, analogous to RHO_CVIMX
!                               option for the objective function.
!        12/06/07       "       Save the keel-crown line for possible use with
!                               HB_GRID and DPLR2D.
!                               Use XNLCON() to provide greater control over
!                               RHO_CVIMX and the CRVIMX (peak crv.) constraint.
!        12/07/07       "       FRACTION_CRVIMX wasn't being initialized for the
!                               RHO_CVIMX > 0 case (no XNLCON() read if the
!                               analogous CRVIMX constraint is omitted).
!        12/11/07       "       The "NOCUSP" option now imposes 3rd derivative
!                               continuity to ensure smooth curvature.
!        12/12/07       "       Save the centerline curvature distribution for
!                               plotting comparisons (see SUMMARY).
!        01/07/08       "       PERTURB wasn't constraining all spokes for
!                               maximum curvature in the radial direction.
!        02/06/08       "       Added Gaussian curvature constraints (maximum
!                               and total violation of target minimum).
!        02/19/08       "       Omit i = 1 from Gaussian curvature violation
!                               calculation (collapsed cells); expand summary
!                               information so all controls are apparent.
!        02/21/08       "       Replace unused RHO_CG_OFFSET with RHO_TVD_GC
!                               as a way of smoothing Gaussian curvature
!                               on the surface grid (and hence of smoothing
!                               the geometry).
!        02/23/08       "       RHO_TVD_GC doesn't seem to be effective for
!                               asymmetric cases, possibly because the grid
!                               lines don't follow the varying shoulder.
!                               Therefore, suppress the contributions in the
!                               azimuthal direction.
!        09/29/08       "       A slight weakness for off-center apex cases has
!                               been the fact that opposing defining spokes do
!                               not necessarily meet smoothly except for the
!                               imposed zero first derivative.  Centerline
!                               curvature can therefore be discontinuous. What
!                               second derivative to impose there?  About the
!                               only choice is the average of the two, and this
!                               can only be done after the existing "no cusp"
!                               scheme has been applied to all spokes.  It needs
!                               new 6th-degree polynomial routine POLY6.
!        12/06/10       "       Improved the radial point distribution options
!                               for better resolution of forebody shoulders:
!                               (1) If R*_SPL < 0, reverse the distribution, as
!                               found helpful with 0.04, 0, 0.3, 0.66 inputs;
!                               (2) If R*_SPL = 1, use curvature-based spacing;
!                               see further details above.
!        12/14/10       "       Belated retrofit of common symmetric shapes:
!                               calculate the generatrix directly instead of
!                               reading it from heat_shield.geo if the file
!                               heat_shield.analytic is found.  See more above
!                               near GEOMETRY NOTES.
!        02/04/11       "       Setting CO_INIT(1) to RADIUS_BASE was wrong for
!                               the analytic case: it should be RO_INIT(NGEOM).
!        10/25/11       "       Added the biconic option to the analytic case,
!                               via the same heatshield.analytic geometry file.
!                               The number of lines in the file distinguishes
!                               it from the sphere-cone/spherical-section case.
!
!     ORIGINAL RESEARCHER:               ORIGINAL CONSULTANT:
!
!        James Reuther,                  Peter Gage,
!        ASA Branch,                     ELORET Corporation/NASA Ames,
!        NASA Ames Research Center,      Moffett Field, CA
!        Moffett Field, CA
!
!     ORIGINAL SOFTWARE AUTHORS:
!
!        James Reuther,   ASA Branch, NASA Ames Research Center
!        David Saunders,  ELORET Corporation/ASA Branch, NASA Ames
!
!     CURRENT EXTENSIONS (2007-2008):
!
!        Roman Jits and   ELORET Corporation/TSA Branch, NASA Ames
!        David Saunders   (now ERC, Inc./TSA/NASA Ames)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Global variables:

      USE CONST_MOD
      USE DESIGN_VARS_MOD
      USE OPT_MOD
      USE TRIANG_MOD

      IMPLICIT NONE

!     Local constants:

      REAL,      PARAMETER :: ONE = 1.
      CHARACTER, PARAMETER :: PROGNAME*11 = 'HEAT_SHIELD'

!     Local variables:

      INTEGER :: I, IOS, J, NCALL, MODE
      REAL    :: OBJPAR
      LOGICAL :: ANALYTIC, FAIL
      REAL(KIND=4) :: CPU2
!     Procedures:

      EXTERNAL   SOLVE, CONFUN, OBJFUN

!     Execution:
!     !!!!!!!!!!

      CALL SECOND (CPU0)

      EPSMCH    = EPSILON (EPSMCH)
      SMALL     = MAX (3.* EPSMCH, 1.E-13)
      PI        = 4.* ATAN (ONE)
      TWOPI     = 2.* PI
      TWOPISQ   = TWOPI * PI
      SIXROOTPI = 6.* SQRT (PI)  ! For full config., else divide it by sqrt (2)
      RPI       = ONE / PI
      DEG2RAD   = PI / 180.
      RAD2DEG   = 180./ PI
      BIGBND    = 1.E+10

!     Open all files:
!     !!!!!!!!!!!!!!!

!     Enter the control file (with any name) on the command line now.

!!!   OPEN (LREAD,  FILE='heat_shield.inp', STATUS='OLD') ! Standard input

!     Option to construct a generatrix for a sphere-cone or biconic forebody:

      OPEN (LGEOM,  FILE='heat_shield.analytic', STATUS='OLD', IOSTAT=IOS)

      ANALYTIC = IOS == 0
      IF (.NOT. ANALYTIC) THEN
         OPEN (LGEOM, FILE='heat_shield.geo', STATUS='OLD')
      END IF

      OPEN (LWRIT,  FILE='heat_shield.out', STATUS='UNKNOWN')
      OPEN (LGSAVE, FILE='heat_shield.opt', STATUS='UNKNOWN')
      OPEN (LXYZ,   FILE='heat_shield.xyz', STATUS='UNKNOWN')

      WRITE (LWRIT, '(/, (A))') &
         ' HEAT_SHIELD Version:  October 25, 2011', &
         ' Sponsor:              Aerothermodynamics Branch', &
         '                       NASA Ames Research Center'

!     Read optimization controls:
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!

      CALL RD_OPT

!     Read initial geometry:
!     !!!!!!!!!!!!!!!!!!!!!!

      CALL RD_GEOM (ANALYTIC)

!     Read design variables:
!     !!!!!!!!!!!!!!!!!!!!!!

      CALL RD_VARS

!     Read linear and nonlinear constraints:
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      CALL RD_CONS

!     Set up the linear constraint matrix and corresponding bounds.
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF (NCLIN > 0) CALL SETLCON

      NDITER = 0
      NCALL  = -3 ! Signals the first objective function call

!     Constrained optimization:
!     !!!!!!!!!!!!!!!!!!!!!!!!!

      CALL SETUP_NPOPT (MODE)

      IF (MODE /= 0) GO TO 999

      CALL NPOPT (NDV, NCLIN, NCNLN, NROWA, NROWJ, NROWR, AMAT, BL, BU, &
                  CONFUN, OBJFUN, INFORM, MAJITS, ISTATE,               &
                  CVEC, CJAC, CLAMDA, OBJPAR, GRAD, RMAT, V,            &
                  IWVEC, LENIW, WVEC, LENW)

!     Ensure we've got the optimized result?

      IF (NITMAX > 0) THEN

         NCALL = -4 ! Suppresses geometry output already done

         CALL SOLVE (NDV, V, OBJPAR, NCALL, FAIL)

         WRITE (LWRIT, '(/, A, I7)') &
            ' Total # objective calculations:', NFTOTL
      END IF

!     Save the surface triangulation and Cp distribution in Tecplot format.
!     Save the triangulated surface with a Cp at each vertex of each face,
!     rather than trying to use the existing connectivity:

      OPEN (UNIT=LTRIANG, FILE='heat_shield_tri.dat', STATUS='UNKNOWN')

      WRITE (LTRIANG, '(A)') 'VARIABLES = "X", "Y", "Z", "Cp"'
      WRITE (LTRIANG, '(A, I6, A, I6, A)') &
         'ZONE N=', 3*N_Triangles, ', E=', N_Triangles, &
         ', F=FEPOINT, ET=TRIANGLE'
      WRITE (LTRIANG, '(1P, 4E15.7)') &
          ((XYZ(:,IPT_TRI(I,J)), CP(J), I = 1, 3), J = 1, N_Triangles)
      WRITE (LTRIANG, '(3I8)') (I, I = 1, 3*N_Triangles)
      CLOSE (LTRIANG)

!     Save centroid (x,y,z,Cp)s for possible structural analysis.
!     Tecplot can also show this dataset as scattered data with symbols
!     colored by Cp, but it deals better with the above vertex-centered data.

!!!   OPEN (UNIT=LCPS,  FILE='heat_shield.xyzcp', STATUS='UNKNOWN')

!!!   WRITE (LCPS, '(I6, A)') N_Triangles, ' ! Heat shield (x,y,z,Cp)s'
!!!   WRITE (LCPS, '(1P, 3E15.7, 0P, F11.7)') &
!!!        (XYZ_CENTROID(1:3,J), CP(J), J = 1, N_Triangles)
!!!   CLOSE (LCPS)

      WRITE (LWRIT, '(/, A)') ' Normal termination.'


  999 CALL SECOND (CPU2)

      CALL CPUTIME (CPU0, CPU2, PROGNAME, 'for this run', LWRIT)

      END PROGRAM HEAT_SHIELD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE CHKLCON ()
!
!     Check the perturbed design against the linear constraints to
!     give a more readable table than is provided by the optimizer.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Global variables:

      USE CONST_MOD
      USE DESIGN_VARS_MOD
      USE OPT_MOD

      IMPLICIT NONE

!     Local variables:

      INTEGER   INTERP, M
      REAL      BOUNDL, BOUNDU, T, VALUE
      CHARACTER FLAG*3

!     Execution:

      WRITE (LWRIT, 1010)

!     Evaluate the linear constraints in the form seen by NPOPT:

      CONLIN = MATMUL (AMAT, V) ! CALL AX (NCLIN, NCLIN, NDV, AMAT, V, CONLIN)

      DO M = 1, NCLIN

         BOUNDL = BLLIN(M)
         BOUNDU = BULIN(M)
         T      = TLCON(M)
         INTERP = ILCON(M)
         VALUE  = 999.

         SELECT CASE (LCTYPE(M))

         END SELECT

         BOUNDL = BLLIN(M)
         BOUNDU = BULIN(M)

         IF (VALUE < BOUNDL .OR. VALUE > BOUNDU) THEN
            FLAG = '***'
         ELSE
            FLAG = '   '
         END IF

         WRITE (LWRIT, 1030) M, LCTYPE(M), ILCON(M), VALUE, &
                             BOUNDL, BOUNDU, CONLIN(M), FLAG
      END DO

      RETURN

 1010 FORMAT (/, ' CHKLCON: Linear constraint check:', //,  &
         '    #  LCTYPE   ILCON  CURRENT VALUE',            &
         '    LOWER BOUND    UPPER BOUND           AMAT * V')
 1030 FORMAT (1X, I4, 2X, A6, I6, 2X, 3F15.6, E19.10, 2X, A3)

      END SUBROUTINE CHKLCON

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE CONFUN (MODE, NCPAR, NDVPAR, NROWJPAR, NEEDC, VPAR, &
                         CPAR, CJACPAR, NSTATE)
!
!     CONFUN computes the nonlinear constraint functions C(*) and their
!     gradients.  Its calling sequence is that required by NPSOL/NPOPT.
!     Finite differencing uses the "H"s read with each design variable.
!
!     This version ignores NEEDC(*) -  it is more efficient to evaluate all
!     constraints in groups than to evaluate each one independently.
!
!     This version also allows CONFUN to evaluate the objective function and
!     its derivatives.  Otherwise, the optimization repeats the same calls to
!     SOLVE for the objective and for the nonlinear constraints.
!
!     Since OBJFUN is called following each CONFUN call, it appears safe for
!     CONFUN to carry around the most recent objective and its gradient.  (But
!     OBJFUN and SOLVE still have to handle the no-nonlinear-constraint case.)
!
!     08/07/07  D. Saunders  GRAD_MODE = 3 allows central differencing.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Global quantities:

      USE CONST_MOD
      USE OPT_MOD

      IMPLICIT NONE

!     Arguments:

      INTEGER    MODE     ! I MODE = 2 means the elements of C(*) indicated
                          !            by NEEDC(*) > 0 need to be set (and
                          !            similarly for the available elements
                          !            of the rows of CJAC(*,*));
                          !   MODE = 1 means the available elements of the
                          !            rows of CJAC(*,*) corresponding to
                          !            NEEDC(*) > 0 must be set; C(*) and
                          !            other rows of CJAC will be ignored;
                          !   MODE = 0 means the elements of C(*) indicated
                          !            by NEEDC(*) > 0 must be set; CJAC and
                          !            other elements of C(*) are ignored;
                          ! O MODE =-1 means terminate the solution of the
                          !            current problem.
      INTEGER    NCPAR    ! I The number of nonlinear constraints (NCNLN).
      INTEGER    NDVPAR   ! I The number of design variables (NDV).
      INTEGER    NROWJPAR ! I Maximum no. of nonlinear constraints (NROWJ).
      INTEGER  NEEDC(NCPAR) ! I Indices of the elements of C(*) and CJAC(*,*)
                          !   which must be returned by CONFUN (ignored).
      REAL   VPAR(NDVPAR) ! I Variables at which to evaluate constraints.
      REAL   CPAR(NCPAR)  ! O Constraint values, returned if MODE = 0 or 2
                          !   for constraints indicated by NEEDC(*) > 0.
      REAL   CJACPAR(NROWJPAR,NDVPAR)
                          ! O Constraint derivatives; see MODE
      INTEGER    NSTATE   ! I NSTATE = 1 means this is the first call (which
                          !              precedes the first call to OBJFUN).

!     Local constants:

      REAL,      PARAMETER :: ZERO = 0., HALF = 0.5
      LOGICAL,   PARAMETER :: NOPRINT = .FALSE.
      CHARACTER, PARAMETER :: SUBNAME * 6 = 'CONFUN'

!     Local variables:

      INTEGER    I, J, NCALL
      REAL       HCEN, OBJFOR, OBJBAK, VTEMP
      LOGICAL    FAIL
      REAL(KIND=4) :: TIME1, TIME2
!     Execution:

!!!   WRITE (LWRIT, '(/, A, 2I4)') ' CONFUN: MODE & NSTATE: ', MODE, NSTATE

      IF (NSTATE == 1) THEN ! First call to CONFUN (precedes OBJFUN)

         NCALL = -13

      ELSE ! Not the first call to CONFUN

         IF (MODE == 0) THEN ! Constraints only
            NCALL = -10
         ELSE                ! Constraints + derivatives
            NCALL = -12
         END IF

      END IF

!     Analyze the current heat shield configuration:

      CALL SOLVE (NDVPAR, VPAR, OBJFOR, NCALL, FAIL)

      NFTOTL = NFTOTL + 1

      IF (FAIL) THEN
         WRITE (LWRIT, '(/, 1X, 2A)') &
            SUBNAME, ':  Bad return from SOLVE.'
         GO TO 900
      END IF

      OBJ_CONFUN = OBJFOR  ! Avoid repeat: OBJFUN uses it

!     Evaluate all nonlinear constraints, scaled:

      CALL NLCON (NOPRINT, CPAR)

      IF (NITMAX == 0) THEN
         IF (CONFUN_DOES_OBJ) THEN
            WRITE (LWRIT, '(/, 1X, 2A)') &
               SUBNAME, ':  Forcing exit from NPOPT (NITMAX = 0).'
            GO TO 900
         ELSE ! No point in doing the constraint gradients
            CJACPAR = 999. ! Keep NPOPT happy
            GO TO 999
         END IF
      END IF

!     Constraint derivatives as well?

      IF (MODE > 0) THEN

         CALL SECOND (TIME1)

         IF (GRAD_MODE == 2) THEN ! Forward differencing

!           Perturb a given variable only once for all constraints:

            DO J = 1, NDV

               VTEMP   = VPAR(J)
               VPAR(J) = VTEMP + AITCH(J)

               CALL SOLVE (NDVPAR, VPAR, OBJFOR, J, FAIL) ! NCALL = J here

               VPAR(J) = VTEMP

               CALL NLCON (NOPRINT, CJACPAR(1,J))

               DO I = 1, NCNLN
                  CJACPAR(I,J) = (CJACPAR(I,J) - CPAR(I)) * ONEOVERH(J)
               END DO

               FFWD(J) = (OBJFOR - OBJ_CONFUN) * ONEOVERH(J)

            END DO

            NFTOTL = NFTOTL + NDV

         ELSE ! GRAD_MODE = 3: central differencing

!           Perturb a given variable twice for all constraints:

            DO J = 1, NDV

               VTEMP   = VPAR(J)
               HCEN    = AITCH(J) * EPSM6TH ! Since eps**(1/3) is ~optimal
               VPAR(J) = VTEMP + HCEN

               CALL SOLVE (NDVPAR, VPAR, OBJFOR, J, FAIL) ! NCALL = J here

               CALL NLCON (NOPRINT, CJACPAR(1,J))

               VPAR(J) = VTEMP - HCEN

               CALL SOLVE (NDVPAR, VPAR, OBJBAK, J, FAIL)

               NFTOTL  = NFTOTL + 2
               VPAR(J) = VTEMP

               CALL NLCON (NOPRINT, CJACBAK)

               HCEN = HALF / HCEN

               DO I = 1, NCNLN
                  CJACPAR(I,J) = (CJACPAR(I,J) - CJACBAK(I)) * HCEN
               END DO

               FFWD(J) = (OBJFOR - OBJBAK) * HCEN

            END DO

         END IF

         CALL SECOND (TIME2)

         CALL CPUTIME (TIME1, TIME2, SUBNAME, &
          'to evaluate derivatives of nonlinear constraints & objective', LWRIT)

         WRITE (LWRIT, '(/, A)' ) &
            ' Gradients of objective and constraints:', &
            ' IVAR       dOBJ/dVAR   dC1/dVAR   dC2/dVAR ...'

         I = MIN (NCNLN, 10)

         DO J = 1, NDV
            WRITE (LWRIT, '(I4, 1P, E17.8, 10E11.3)') J, FFWD(J), CJACPAR(1:I,J)
         END DO

      END IF

      GO TO 999

  900 MODE = -1 ! NPOPT should quit

  999 RETURN

      END SUBROUTINE CONFUN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE CPUTIME (CPU1, CPU2, CALLER, DESCRIPTION, LUN)
!
!     CPUTIME modularizes printing of CPU time between two measures of total
!     CPU time used so far as given by the SECOND utility.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IMPLICIT NONE

!     Arguments:

      REAL,      INTENT (IN) :: CPU1, CPU2
      CHARACTER, INTENT (IN) :: CALLER * (*), DESCRIPTION * (*)
      INTEGER,   INTENT (IN) :: LUN

!     Local variables:

      REAL DIFF

!     Execution:

      DIFF = CPU2 - CPU1
      IF (DIFF >= 0.01) THEN
         WRITE (LUN, '(/, 1X, 4A, F10.2)') &
            CALLER, ': CPU secs. ', DESCRIPTION, ':   ', DIFF
      END IF

      END SUBROUTINE CPUTIME

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE KAPPA (N, NCRV_MAX, R, X, S, CRV_TARG, CRV_MAX, CRV_TOT)
!
!     Calculate the local curvatures of a heat shield section in (r,x) space
!     and accumulate the sum of squares of the values less than the target.
!     This is done on the geometry rather than on the surface grid, to avoid
!     difficulties introduced by splines.  End points are ignored because the
!     3-point formulas just give repeats for the values next to the end-points.
!
!     01/29/02  The same utility works for the azimuthal direction, although
!               the nomenclature is off.
!     08/24/07  Ignoring shoulder region in radial direction requires another
!               argument.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IMPLICIT NONE

!     Arguments:

      INTEGER, INTENT (IN)    :: N,       & ! # pts. in the section
                                 NCRV_MAX   ! Last pt. to include in max. curv.
      REAL,    INTENT (IN)    :: R(N),    & ! (r,x) coords. and corresponding
                                 X(N),    & ! arcs needed outside for gridding
                                 S(N),    &
                                 CRV_TARG   ! Lower bound, probably zero
      REAL,    INTENT (INOUT) :: CRV_MAX, & ! Maximum curvature found so far
                                 CRV_TOT    ! Cumulative sum of violations**2

!     Local constants:

      REAL, PARAMETER :: ZERO = 0.

!     Local variables:

      INTEGER :: &
         I
      REAL, DIMENSION (N) :: &
         CURV, RS, RSS, XS, XSS

!     Execution:

!     Central differences at interior points:

      DO I = 2, N - 1
         CALL FDCNTR (I, S, R, RS(I), RSS(I))
         CALL FDCNTR (I, S, X, XS(I), XSS(I))
      END DO

!     Signed local curvatures (> 0 = convex):

      CALL CURV2D (N - 2, RS(2), RSS(2), XS(2), XSS(2), CURV(2))

!     Accumulate violations:

      DO I = 2, N - 1
         CRV_TOT = CRV_TOT + (MIN (CURV(I) - CRV_TARG, ZERO)) ** 2
      END DO

!     Ignore shoulder region for max. curvature region:

      DO I = 2, NCRV_MAX
         CRV_MAX = MAX (CURV(I), CRV_MAX)
      END DO

      END SUBROUTINE KAPPA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE LOFT (TYPE2, Z, ZCEN, ZMIN, ZMAX, PWR, CIRCLE, HEIGHT, LUNERR)
!
!     LOFT modularizes the choice between polynomial and Hicks-Henne-type
!     spanwise lofting of the effect of a wing design variable, for use by
!     both PERTURB and SETLCON.  HEIGHT is returned between 0. and 1.
!
!     08/28/00  DAS  SINC (cyclic) allows for applying symmetric sine bumps
!                    (with centers at 0.5 on [0,1]) at any center which behaves
!                    periodically across the 0/1 boundary.
!     01/23/02   "   The half-shield case needed a new CIRCLE argument.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IMPLICIT NONE

!     Arguments:

      CHARACTER, INTENT (IN)  :: TYPE2 * 4
      INTEGER,   INTENT (IN)  :: LUNERR
      REAL,      INTENT (IN)  :: Z, ZCEN, ZMIN, ZMAX, PWR, CIRCLE
      REAL,      INTENT (OUT) :: HEIGHT

!     Local constants:

      REAL, PARAMETER :: HALF = 0.5, ONE = 1., ZERO = 0.

!     Local variables:

      REAL      ZCNORM, ZCUSE, ZEVAL
      CHARACTER TYPEUSE * 4

!     Execution:

      IF (TYPE2 == 'POLY') THEN

!        Polynomial-type lofting (often linear or constant):
!        N.B.: The constant case requires special treatment at end-pts.

         CALL RIPPER (Z, ZCEN, ZMIN, ZMAX, PWR, HEIGHT)

      ELSE

!        Hicks-Henne-type lofting:

         HEIGHT = 0. ! SFEVAL adds
         ZCNORM = (ZCEN - ZMIN) / (ZMAX - ZMIN) ! Reqd. by SFEVAL

         TYPEUSE = TYPE2(1:4)

         IF (TYPEUSE == 'SINC') THEN ! Slide the 0.5 sine bump around the circle
            TYPEUSE = 'SIN '
            ZCUSE = HALF
            ZEVAL = HALF + (Z - ZCEN) * CIRCLE ! CIRCLE = 1. (full), 0.5 (half)
            IF (ZEVAL < ZERO) ZEVAL = ZEVAL + ONE
            IF (ZEVAL > ONE ) ZEVAL = ZEVAL - ONE
         ELSE
            ZEVAL = Z
            ZCUSE = ZCNORM
         END IF

         CALL SFEVAL (TYPEUSE, PWR, ZCUSE, ZMIN, ZMAX, &
                      ONE, 1, 1, ZEVAL, HEIGHT, LUNERR)

      END IF

      END SUBROUTINE LOFT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE NEWTONIAN_CLCDCM (NFACES, FACE_AREA, XYZ_NORM,    &
                                   XYZ_CENTROID, XYZ_CG,           &
                                   REF_LEN, REF_AREA, GAMMA,       &
                                   RMACH, ALPHA, BETA,             &
                                   CL, CD, CS, CM_P, CM_Y, CM_R, CP)
!
!     This is a Fortran 90 translation of the Newtonian aerodynamics utility
!     originally written in C by Alvin Ramsey.  Input is a triangulated surface;
!     output are lift/drag/moment coefficients determined by summing the effects
!     of the free stream on the forward-facing triangles.  Results have validity
!     for hypersonic Mach numbers only.
!
!     ??/??/00  A. Ramsey    Original implementation.
!     09/05/00  J. Reuther   Initial F90 translation.
!     10/18/00  D. Saunders  Argument-driven version.
!     10/29/00      "        Tried the more elaborate formula in Anderson's 1989
!                            book (Hypersonic & High Temperature Gas Dynamics).
!     11/28/00      "        Return the Cps calculated for each triangle.
!     12/13/00      "        Work with FACE_AREA(*), not PROJECTED_AREA(1:3,*).
!     ??/??/01  J. Reuther   Added wetted area calculation.
!
!     Sponsor:  NASA Ames Research Center, Mountain View, CA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IMPLICIT NONE

!     Arguments:

      INTEGER, INTENT (IN) :: &
         NFACES                  ! Number of triangles

      REAL, INTENT (IN) ::    &
         FACE_AREA(NFACES)       ! Triangle areas

      REAL, INTENT (IN), DIMENSION (3,NFACES) :: &
         XYZ_NORM,            &  ! Components of triangle outward unit normals
         XYZ_CENTROID            ! Triangle centroids

      REAL, INTENT (IN) ::    &
         XYZ_CG(3),           &  ! Moment center
         REF_LEN,             &  ! Reference length
         REF_AREA,            &  ! Reference area
         GAMMA,               &  ! Ratio of specific heats
         RMACH,               &  ! Free-stream Mach number
         ALPHA, BETA             ! Angle of attack and yaw angle (degrees)

      REAL, INTENT (OUT) ::   &
         CL, CD, CS,          &  ! Aerodynamic coefficients
         CM_P, CM_Y, CM_R,    &
         CP(NFACES)              ! Newtonian Cps at each triangle (>= 0.)

!     Local constants:

      REAL, PARAMETER :: &
         DEG2RAD = 0.0174532925199433, ONE = 1., TWO = 2., ZERO = 0.

!     Local variables:

      INTEGER :: &
         I

      REAL :: &
         AB, CA, CB, CBETA, CP_MAX, CX, CY, CZ, GMSQ, RMSQ, &
         PC, S1, S2, S3, RAD, SA, SB, T1, XM, YM, ZM

!     Execution:

      RAD = DEG2RAD * ALPHA
      SA  = SIN (RAD)
      CA  = COS (RAD)
      RAD = DEG2RAD * BETA
      SB  = SIN (RAD)
      CB  = COS (RAD)
      AB  = CA * CB

      CX  = ZERO
      CY  = ZERO
      CZ  = ZERO
      XM  = ZERO
      YM  = ZERO
      ZM  = ZERO

      RMSQ = RMACH * RMACH
!!!!!!CP_MAX = 4.* (ONE - ONE / RMSQ) / (GAMMA + ONE)

      IF (GAMMA == ONE) THEN
         CP_MAX = TWO * (RMSQ - ONE) / RMSQ
      ELSE
         GMSQ = GAMMA * RMSQ
         T1   = ((GAMMA + ONE) * RMACH) ** 2 /   &
                (TWO * (TWO * GMSQ - GAMMA + ONE))
         CP_MAX = (TWO / GMSQ) * (T1 ** (GAMMA / (GAMMA - ONE)) *   &
                  ((ONE - GAMMA + TWO * GMSQ) / (GAMMA + ONE)) - ONE)
      END IF

!     Loop over all triangles and integrate forces.

      DO I = 1, NFACES

!        Calculate the dot product of the outward pointing normal
!        of the triangle and the free-stream direction.

         CBETA = AB * XYZ_NORM(1,I) + &
                 SA * XYZ_NORM(2,I) + &
                 SB * XYZ_NORM(3,I)

!        If the triangle doesn't point into the free-stream, suppress it.

         IF (CBETA >= ZERO) THEN

            CP(I) = ZERO

         ELSE ! Accumulate forces and moments.

            CP(I) = CP_MAX * CBETA * CBETA
            PC    = -CP(I) * FACE_AREA(I)
            S1    = PC * XYZ_NORM(1,I)
            S2    = PC * XYZ_NORM(2,I)
            S3    = PC * XYZ_NORM(3,I)
            CX    = CX + S1
            CY    = CY + S2
            CZ    = CZ + S3
            XM    = XM - S2 * (XYZ_CENTROID(3,I) - XYZ_CG(3)) + &
                         S3 * (XYZ_CENTROID(2,I) - XYZ_CG(2))
            YM    = YM - S3 * (XYZ_CENTROID(1,I) - XYZ_CG(1)) + &
                         S1 * (XYZ_CENTROID(3,I) - XYZ_CG(3))
            ZM    = ZM + S2 * (XYZ_CENTROID(1,I) - XYZ_CG(1)) - &
                         S1 * (XYZ_CENTROID(2,I) - XYZ_CG(2))
         END IF

      END DO

!     Finish up force calculations.

      AB = ONE / REF_AREA
      CX = CX * AB
      CY = CY * AB
      CZ = CZ * AB

      CL = ((CY * CA) - (CX * SA))
      CD = ((CX * CA) + (CY * SA)) * CB + &
           ( CZ                  ) * SB
      CS = ((CX * CA) + (CY * SA)) * SB - &
           ( CZ                  ) * CB

      AB   = AB / REF_LEN
      CM_R = XM * AB
      CM_Y = YM * AB
      CM_P = ZM * AB

      END SUBROUTINE NEWTONIAN_CLCDCM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE NLCON (PRINT, CPAR)
!
!     NLCON evaluates all nonlinear constraints, CPAR(*), with the option to
!     tabulate UNscaled results.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Global variables:

      USE AERO_MOD
      USE COEFFS_MOD
      USE CONST_MOD
      USE GEOM_MOD
      USE GRID_MOD
      USE OPT_MOD

      IMPLICIT NONE

!     Arguments:

      LOGICAL, INTENT (IN)  :: PRINT   ! T means tabulate UNscaled constraints
      REAL,    INTENT (OUT) :: CPAR(*) ! Nonlinear constraints (= C(*) in OPT)

!     Local constants:

      REAL, PARAMETER :: ZERO = 0.

!     Local variables:

      INTEGER   I, K, M
      REAL      BOUNDL, BOUNDU, SCALE, TIME, VALUE
      CHARACTER FLAG*3

!     Execution:

      K = NDV + NCLIN ! Offset for the bounds

      IF (PRINT) WRITE (LWRIT, 1010)

      DO M = 1, NCNLN

         SELECT CASE (NLCTYPE(M))

            CASE ('H_FLUX')
               VALUE = HEAT_FLUX

            CASE ('W_TEMP')
               VALUE = TEMP_MAX

            CASE ('CM    ')
               VALUE = MN

            CASE ('CD    ')
               VALUE = CD

            CASE ('LD    ')
               VALUE = LD

            CASE ('CD_A  ')
               VALUE = CD * SECTIONAL_AREA

            CASE ('CM_ALF')
               VALUE = DCMDA  ! Derivative of CM w.r.t. Alpha

            CASE ('CG_OFF')
               VALUE = 999.   ! TBD (CG offset)

            CASE ('ALPHA ')
               VALUE = Alpha

            CASE ('AREA  ')
               VALUE = SECTIONAL_AREA

            CASE ('VOLUME')
               VALUE = VOLUME

            CASE ('V_EFF ')
               VALUE = V_EFF

            CASE ('WAREA ')
               VALUE = (SWET - SBASE) / REF_AREA

            CASE ('CURVI ')
               VALUE = CRVI_TOT ! Total of radial curvature violations (squared)

            CASE ('CURVJ ')
               VALUE = CRVJ_TOT ! Total of azimu. curvature violations (squared)

            CASE ('CURVG ')
               VALUE = CRVG_TOT ! Total of Gauss. curvature violations (squared)

            CASE ('CRVIMX')
               VALUE = CRVI_MAX ! Peak curvature in I direction, XNLCON() * NI

            CASE ('CRVGMX')
               VALUE = CRVG_MAX ! Peak Gaussian curvature over full grid

            CASE ('SMTH_X')
               VALUE = X_RESIDUAL

            CASE DEFAULT
               WRITE (LWRIT, '(/, 2A)') ' NLCON: Bad constraint type = ', &
                  NLCTYPE(M)
               STOP

         END SELECT

         SCALE   = SNLCON(M)
         CPAR(M) = VALUE * SCALE

         IF (PRINT) THEN
            K = K + 1
            BOUNDL = BL(K) / SCALE ! Since the main program applies SNLCON(M)
            BOUNDU = BU(K) / SCALE

            IF (VALUE < BOUNDL .OR. VALUE > BOUNDU) THEN
               FLAG = '***'
            ELSE
               FLAG = '   '
            END IF

            WRITE (LWRIT, 1030) M, NLCTYPE(M), VALUE, BOUNDL, BOUNDU, &
                                XNLCON(M), INLCON(M), FLAG
          END IF

      END DO

      RETURN

 1010 FORMAT (/, '  # NLCTYPE        VALUE  LOWER BOUND  UPPER BOUND', &
              '    XNLCON  INLCON  VIOLATION?')
 1030 FORMAT (1X, I2, 2X, A6, 3F13.6, F10.4, I8, 5X, A3)

      END SUBROUTINE NLCON

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE NOCUSP (RCUTOFF, N, R, X)
!
!     Replace points 2 : M - 1 with points on a quintic defined by (R1, X1, 0.)
!     & (R2, X2, D1, D2, D3), where R2 = R(1) + RCUTOFF, X2 is the corresponding
!     X, and M is the first point for which R(M) - R(1) > CUTOFF.  D1, D2 and D3
!     are the corresponding derivatives of X w.r.t. R at the cut-off pt.
!
!     11/02/01  D.Saunders  Initial implementation to deal with cusps
!                           at the apex of a heat shield (quintic).
!     11/06/01      "       Regularized geometry doesn't get saved.
!                           Therefore, apply this to the raw geometry
!                           and specify a distance, not an index.
!     11/14/01      "       Peter Gage says y" = 0 at the apex isn't
!                           what we meant - change the quintic to a quartic.
!     02/11/02      "       Occasional bad gradients arise from just picking
!                           the defining point nearest to the cut-off.
!                           Therefore, use the exact cut-off location.
!                           (NOTE: Specifying an index M and assuming all
!                           spoke geometries have consistent distributions
!                           did not work well, for unknown reasons.)
!     02/13/02      "       Lack of 2nd derivative continuity for local splines
!                           must explain occasional bad gradients - go to less
!                           convenient conventional splines.
!     11/15/02      "       Spline just the inner portions, to avoid non-
!                           monotonic r values at the rim.
!     12/11/07      "       Seek 3rd derivative continuity at the cut-off to
!                           ensure curvature smoothness.  No need for the
!                           conventional spline derivatives now, but match at
!                           a fixed index defined by the very first one chosen.
!                           The geometry points are assumed to be consistent
!                           from spoke to spoke.  The 02/11/02 note should
!                           not be an issue.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IMPLICIT NONE

!     Arguments:

      REAL,    INTENT (IN) :: RCUTOFF      ! A fraction of  L_ref; see above

      INTEGER, INTENT (IN) :: N            ! Number of points in R and X

      REAL, INTENT (INOUT) :: R(N), X(N)   ! Coordinates being fudged

!     Local constants:

      REAL, PARAMETER :: ZERO = 0.

!     Local variables:

      INTEGER :: I, IER

      INTEGER, SAVE :: M = 999  ! Keep using the index found on the first call

      REAL :: C(5), D1, D2, D3, H

!     Execution:

      IF (M == 999) THEN
         M = N
         DO I = 3, N
            IF (R(I) > RCUTOFF) THEN
               M = I
               EXIT
            END IF
         END DO
         WRITE (1, '(/, A, I4)') ' NOCUSP index M:', M
      END IF

!     We need continuous 1st, 2nd & 3rd derivatives at the cut-off point.

      CALL DERIVS_123 (N, R, X, M, D1, D2, D3)

!     Quintic defined by end-point quantities:

      CALL QUINTIC2 (R(1), X(1), ZERO, R(M), X(M), D1, D2, D3, C)

      DO I = 2, M - 1
        H = R(I) ! - 0.
        X(I) = X(1) + H*(C(1) + H*(C(2) + H*(C(3) + H*(C(4) + H*C(5)))))
      END DO

      END SUBROUTINE NOCUSP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE OBJECTIVE (OBJ, NCALL)
!
!     OBJECTIVE is called by SOLVE.  It keeps the objective function details
!     in one place.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Global variables:

      USE AERO_MOD
      USE COEFFS_MOD
      USE CONST_MOD
      USE DESIGN_VARS_MOD
      USE GRID_MOD
      USE OPT_MOD

      IMPLICIT NONE

!     Arguments:

      INTEGER :: NCALL
      REAL, INTENT (OUT) :: OBJ

!     Local constants:

      REAL, PARAMETER :: ONE = 1., ZERO = 0.

!     Local variables:

      INTEGER :: I, K

      REAL :: &
         PENA, PENB, PENC, PEND, PENE, PENF, PENG, PENH, PENI, PENJ, PENK,     &
         PENL, PENM, PENN, PENO, PENP, PENQ, PENR, PENS, PENT, PENU, SUM

!     Execution:

      PENA = RHO_CM_TARG
      PENB = RHO_CM_MAX
      PENC = RHO_CD
      PEND = RHO_LD
      PENE = RHO_TEMP_MAX
      PENF = RHO_CD_INV
      PENG = RHO_CD_TARG
      PENH = RHO_AL
      PENI = RHO_CURVI
      PENJ = RHO_CURVJ
      PENK = RHO_CVIMX
      PENL = RHO_CVJMX
      PENM = RHO_SMTH_X
      PENN = RHO_SWET
      PENO = RHO_VNORM
      PENP = RHO_AREA
      PENQ = RHO_VOLUME
      PENR = RHO_V_EFF
      PENS = RHO_CD_A
      PENT = RHO_CM_ALF
      PENU = RHO_TVD_GC

      IF (PENA > SMALL) PENA = PENA * (MN - CM_TARG) ** 2
      IF (PENB > SMALL) PENB = PENB * (ONE - ABS (MN))
      IF (ABS (PENC) > SMALL) PENC = PENC * CD
      IF (ABS (PEND) > SMALL) PEND = PEND * (-LD)
      IF (PENE > SMALL) PENE = PENE * MAX (TEMP_MAX - TM_TARG, ZERO)

      IF (PENF > SMALL) THEN
         IF (ABS (CD) > SMALL) THEN
            PENF = PENF / CD
         ELSE
            PENF = PENF * 1000.
         END IF
      END IF

      IF (PENG > SMALL) PENG = PENG * (CD - CD_TARG) ** 2

      IF (PENH > SMALL) THEN
         IF (ALPHA < AL_TARG) THEN
            PENH = PENH * (AL_TARG - ALPHA) ** 2
         ELSE
            PENH = ZERO
         END IF
      END IF

      IF (PENI > SMALL) PENI = PENI * CRVI_TOT
      IF (PENJ > SMALL) PENJ = PENJ * CRVJ_TOT
      IF (PENK > SMALL) PENK = PENK * CRVI_MAX
      IF (PENL > SMALL) PENL = PENL * CRVJ_MAX
      IF (PENM > SMALL) PENM = PENM * X_RESIDUAL

      IF (PENN > SMALL) THEN
         PENN = PENN * MAX (ZERO, (SWET - SBASE) / REF_AREA)
      ELSE IF (PENN < -SMALL) THEN
         PENN = PENN * MIN (ZERO, (SWET - SBASE) / REF_AREA)
      END IF

      IF (PENO > SMALL) THEN ! 2-norm-squared of sine bump variables
         SUM = ZERO
         DO I = NPLANF + 1, NPLANF + NBUMPS
            SUM = SUM + V(I) ** 2
         END DO
         PENO = PENO * SUM
      END IF

      IF (PENP > SMALL) PENP = PENP * (-SECTIONAL_AREA)
      IF (PENQ > SMALL) PENQ = PENQ * (-VOLUME)
      IF (PENR > SMALL) PENR = PENR * (-V_EFF)
      IF (PENS > SMALL) PENS = PENS * (-CD * SECTIONAL_AREA)
      IF (PENT > SMALL) PENT = PENT * DCMDA
      IF (PENU > SMALL) PENU = PENU * TVD_GC

      OBJ = PENA + PENB + PENC + PEND + PENE + PENF + PENG + PENH + PENI +    &
            PENJ + PENK + PENL + PENM + PENN + PENO + PENP + PENQ + PENR +    &
            PENS + PENT + PENU

      IF (NCALL < 0) THEN
         WRITE (LWRIT, '(/, A, //, (1X, A, 1P, E25.15, 0P, F11.5))') &
            ' Contributions to the objective and associated weights:', &
            'CM target  :', PENA, RHO_CM_TARG, &
            'CM max.    :', PENB, RHO_CM_MAX,  &
            'CD         :', PENC, RHO_CD,      &
            '-L/D       :', PEND, RHO_LD,      &
!!          'Temp. max. :', PENE, RHO_TEMP_MAX,&
            '1/CD       :', PENF, RHO_CD_INV,  &
            'CD target  :', PENG, RHO_CD_TARG, &
            'Alpha targ.:', PENH, RHO_AL,      &
            'CURVI      :', PENI, RHO_CURVI,   &
            'CURVJ      :', PENJ, RHO_CURVJ,   &
            'CVIMX      :', PENK, RHO_CVIMX,   &
            'CVJMX      :', PENL, RHO_CVJMX,   &
            'SMOOTH_X   :', PENM, RHO_SMTH_X,  &
            'Swet       :', PENN, RHO_SWET,    &
            'Vnorm^2    :', PENO, RHO_VNORM,   &
            '-Area      :', PENP, RHO_AREA,    &
            '-Volume    :', PENQ, RHO_VOLUME,  &
            '-Veffective:', PENR, RHO_V_EFF,   &
            '-CD x A    :', PENS, RHO_CD_A,    &
            'dCM/dAlpha :', PENT, RHO_CM_ALF,  &
            'TVD Gauss. :', PENU, RHO_TVD_GC
      END IF

      END SUBROUTINE OBJECTIVE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE OBJFUN (MODE, NDVPAR, VPAR, OBJPAR, GRDPAR, NSTATE)
!
!     OBJFUN acts as the interface between NPOPT and SOLVE.  See also CONFUN
!     if nonlinear constraints are present.  Its main task is to determine
!     whether NPOPT wants just the objective function or the objective function
!     plus gradient information, and then call SOLVE appropriately.
!
!     03/96     J. Reuther   Initial implementation (SYN87).
!     02/00     D. Saunders  Initial adaptations (Traj_opt).
!     08/07/07     "         GRAD_MODE = 3 allows for central differencing.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Global variables:

      USE CONST_MOD
      USE OPT_MOD

      IMPLICIT NONE

!     Arguments:

      INTEGER, INTENT (INOUT) :: MODE    ! 0 means return the objective;
                                         ! 1 means return the gradient;
                                         ! 2 means return both;
                                         !-1 output means the obj. calculation
                                         !   was aborted; a search step
                                         !   should be shortened & retried
      INTEGER, INTENT (IN)  :: NDVPAR    ! Number of design variables
      REAL,    INTENT (INOUT) :: VPAR(NDVPAR) ! Design variables
      REAL,    INTENT (OUT) :: OBJPAR    ! Objective function value
      REAL,    INTENT (OUT) :: GRDPAR(NDVPAR) ! Gradient vector
      INTEGER, INTENT (IN)  :: NSTATE    ! 1 means this is the first call

!     Local constants:

      CHARACTER, PARAMETER :: SUBNAME * 6 = 'OBJFUN'

!     Local variables:

      INTEGER :: I, NCALL, NHALF
      REAL    :: OBJFOR, OBJBAK, TEMPH, TEMPV
      REAL(KIND=4) :: TIME1, TIME2
      LOGICAL :: FAIL

!     Execution:

!     Treat the case in which nonlinear constraints have forced calls
!     to SOLVE, thus providing OBJPAR and GRDPAR.  See CONFUN.

      IF (CONFUN_DOES_OBJ) THEN

         OBJPAR = OBJ_CONFUN

         IF (MODE > 0) THEN

            GRDPAR = FFWD ! Elements 1 : NDV;  CONFUN now prints the gradient

!!!         WRITE (LWRIT, 1030) 'OBJFUN: Objective gradient via CONFUN:', ' ', &
!!!            '   I            V (I)            G (I)            H (I)'
!!!         WRITE (LWRIT, '(I5, 3E17.8)') &
!!!            (I, VPAR(I), GRDPAR(I), AITCH(I), I = 1, NDVPAR)
         END IF

      ELSE

!        Case of no nonlinear constraints:

!!!      WRITE (LWRIT, 1040) 'MODE and NSTATE: ', MODE, NSTATE

!        Calculate the objective function.

         IF (NSTATE == 1) THEN ! First OBJFUN call from NPOPT

            NCALL = -3

            CALL SOLVE (NDVPAR, VPAR, OBJPAR, NCALL, FAIL)

            NFTOTL = NFTOTL + 1

            IF (FAIL) THEN
               WRITE (LWRIT, 1040) 'Bad return from call to SOLVE.'
               GO TO 900
            END IF

            IF (NITMAX == 0) THEN
               WRITE (LWRIT, 1040) 'Exiting from NPOPT (NITMAX = 0).'
               GO TO 900
            END IF

         ELSE ! Not the first call to OBJFUN

            IF (MODE == 0) THEN ! Objective function only
               NCALL = 0        ! Line search function evaluation
            ELSE                ! Objective and gradient
               NCALL = -2       ! Recover best line search solution
            END IF

            CALL SOLVE (NDVPAR, VPAR, OBJPAR, NCALL, FAIL)

            NFTOTL = NFTOTL + 1

            IF (FAIL) THEN
               IF (MODE == 0) THEN
                  WRITE (LWRIT, 1040)'NPOPT may shorten step & retry.'
                  OBJPAR = 9999.
                  GO TO 999
               ELSE
                  WRITE (LWRIT, 1040) 'OBJFUN: Aborting.'
                  GO TO 900
               END IF
            END IF

         END IF

         IF (MODE > 0) THEN ! Calculate the gradient too

            CALL SECOND (TIME1)

            IF (GRAD_MODE == 2) THEN ! Forward differencing

               DO NCALL = 1, NDVPAR

                  TEMPV = VPAR(NCALL)
                  TEMPH = AITCH(NCALL)
                  NHALF = 0

   10             CONTINUE

                  VPAR(NCALL) = TEMPV + TEMPH

                  CALL SOLVE (NDVPAR, VPAR, OBJFOR, NCALL, FAIL)

                  NFTOTL = NFTOTL + 1

                  IF (FAIL) THEN
                     WRITE (LWRIT, 1040) &
                        'Function failure for gradient element ', NCALL

                     IF  (NHALF < 3) THEN
                        NHALF = NHALF + 1
                        TEMPH = TEMPH * 0.5
                        WRITE (LWRIT, 1030) 'Halve the stepsize for this var.'
                        GO TO 10
                     END IF

                     OBJFOR = OBJPAR
                     WRITE (LWRIT, 1030) &
                        'Three halvings of gradient step size failed.', &
                        'Proceeding with zero for this component.'
                  END IF

                  GRDPAR(NCALL) = (OBJFOR - OBJPAR) / TEMPH
                  VPAR(NCALL)   = TEMPV

               END DO

            ELSE ! GRAD_MODE = 3 (central differencing)

               DO NCALL = 1, NDVPAR

                  TEMPV = VPAR(NCALL)
                  TEMPH = AITCH(NCALL) * EPSM6TH ! Since eps ** (1/3) is ~optimal
                  VPAR(NCALL) = TEMPV + TEMPH

                  CALL SOLVE (NDVPAR, VPAR, OBJFOR, NCALL, FAIL)

                  VPAR(NCALL) = TEMPV - TEMPH

                  CALL SOLVE (NDVPAR, VPAR, OBJBAK, NCALL, FAIL)

                  GRDPAR(NCALL) = (OBJFOR - OBJBAK) / (TEMPH + TEMPH)
                  VPAR(NCALL)   = TEMPV
                  NFTOTL = NFTOTL + 2

               END DO

            END IF

            CALL SECOND  (TIME2)
            CALL CPUTIME (TIME1, TIME2, 'OBJFUN', &
               'to evaluate the objective gradient', LWRIT)

            WRITE (LWRIT, 1030) 'OBJFUN: Gradient vector:', ' ',    &
               '   I            V (I)            G (I)            H (I)'
            WRITE (LWRIT, '(I5, 3E17.8)') &
               (I, VPAR(I), GRDPAR(I), AITCH(I), I = 1, NDVPAR)

         END IF

      END IF

      GO TO 999

  900 MODE = -1 ! NPOPT should quit

  999 RETURN

 1030 FORMAT (/, (1X, A))
 1040 FORMAT (/, ' OBJFUN: ', A, 2I4)

      END SUBROUTINE OBJFUN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE PERTURB (VPAR, NCALL)
!
!     Variant for one-piece or two-piece heat shield.
!     This version applies the design variables to the defining sections,
!     and interpolates a regular surface mesh.
!     Geometric properties such as cross-sectional area and volume are computed
!     numerically from the interpolated mesh.
!
!     Note:  SOLVE passes MCALL >= -3, not NCALL >= -13.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      USE AERO_MOD
      USE COEFFS_MOD
      USE CONST_MOD
      USE DESIGN_VARS_MOD
      USE GEOM_MOD
      USE GRID_MOD
      USE OPT_MOD

      IMPLICIT NONE

!     Arguments:

      REAL,    INTENT (IN) :: &
         VPAR(*) ! Optimization variables

      INTEGER, INTENT (IN) :: &
         NCALL   ! -3 means first call

!     Local constants:

      REAL,      PARAMETER :: &
         ZERO = 0., ONE = 1., TWO = 2.

      LOGICAL,   PARAMETER :: &
         NEW = .TRUE., NORM = .FALSE.

      CHARACTER, PARAMETER :: &
         LOOSE * 1 = 'B', TIGHT * 1 = 'M'

!     Local variables:

      INTEGER :: &
         I, I1, IER, ISMOOTH, K, N, NCRV_MAX, NTOT

      INTEGER, SAVE :: &
         NI, NK

      REAL :: &
         CPH, DELTAV, DPHI, DR, DX, DY, GAUSSI, GAUSSIM1, GAUSSK, GAUSSKM1,    &
         R, R0, RK, RM1, RN, SMAX, SPH, STOTAL, VS, XR_SCALE, XR_SHIFT

      CHARACTER :: &
         TYP1 * 1, TYP3 * 3, TYP4 * 4

      REAL, ALLOCATABLE :: &
         STEMP(:)

      REAL, ALLOCATABLE, DIMENSION (:), SAVE :: &
         DERIVS, SGEOM, SGRID, SNORM1, SNORM2

      REAL, ALLOCATABLE, DIMENSION (:,:), SAVE :: &
         RRG, XRG

!     Execution:

!     Allocate local gridding work-space on the first call:

      IF (NCALL == -3) THEN

         NI = N_RAD_MAX
         NK = N_AZM_MAX

         ALLOCATE (AZM(NK), COSPHI(NK), SINPHI(NK), DERIVS(NI), &
                   SGEOM(NT_RAD), SGRID(NI), RRG(NI,NT_SEC), XRG(NI,NT_SEC))

!        NOTE:  By construction, all regularized inner defining sections
!               end at (RWRING, XWRING) then all regularized outer sections
!               start at  "       "     so R/XRG(*,*) do not need an extra
!               point in the I direction.

         DPHI = TOT_AZ / REAL (NK - 1)

         DO K = 1, NK - 1
            AZM(K)    = DPHI * REAL (K - 1)
            COSPHI(K) = COS (AZM(K))
            SINPHI(K) = SIN (AZM(K))
         END DO

         AZM(NK)    = TOT_AZ ! Exactly
         SINPHI(NK) = ZERO

         IF (FULL == 0) THEN
            COSPHI(NK) = -ONE
         ELSE
            COSPHI(NK) = +ONE
         END IF

         IF (TWO_PARTS) THEN

            ALLOCATE (SNORM1(N_INNER))

!           Normalized inner spoke point distribution options:

            IF (RI_SPL /= ONE) THEN  ! Airfoil-type distribution

               CALL FOILGRD (N_INNER, ZERO, ONE, ABS (RI_SPL), RI_SPQ, &
                             RI_SPS, RI_SPC, SNORM1)

               IF (RI_SPL < ZERO) THEN

                  ALLOCATE (STEMP(N_INNER))

                  STEMP(:) = SNORM1(:)
                  STOTAL = SNORM1(N_INNER)

                  DO I = 1, N_INNER
                     SNORM1(I) = STOTAL - STEMP(N_INNER - I + 1)
                  END DO

                  DEALLOCATE (STEMP)

               END IF

               IER = 0

            ELSE  ! Curvature-based distribution

               IF (RI_SPS == ZERO .AND. RI_SPC == ZERO) ISMOOTH = 0
               IF (RI_SPS /= ZERO .AND. RI_SPC == ZERO) ISMOOTH = 1
               IF (RI_SPS == ZERO .AND. RI_SPC /= ZERO) ISMOOTH = 2
               IF (RI_SPS /= ZERO .AND. RI_SPC /= ZERO) ISMOOTH = 3

               NTOT = NI_INIT(1)

!              CURVDIS recommends geometric scaling to the unit circle:

               XR_SHIFT = RI1(1)
               XR_SCALE = RI1(NTOT) - RI1(1)
               RI1(:)   = (RI1(:) - XR_SHIFT) / XR_SCALE
               XR_SHIFT = XI1(1)
               XI1(:)   = (XI1(:) - XR_SHIFT) / XR_SCALE

               CALL CURVDIS (NTOT, XI1, RI1, N_INNER, RI_SPQ, ISMOOTH, &
                             LWRIT, TRUE, SNORM1, SNORM1, IER)
               IF (IER /= 0) THEN
                  WRITE (LWRIT, '(/, A, I5)') &
                     '*** Bad CURVDIS return (inner); IER:', IER
                  GO TO 9  ! Keep the number of stops down
               END IF

               DEALLOCATE (XI1, RI1)
               STOTAL = SNORM1(N_INNER)
               SNORM1(:) = SNORM1(:) / STOTAL

            END IF

         END IF

         ALLOCATE (SNORM2(N_OUTER))

!        Normalized inner spoke point distribution options:

         IF (RO_SPL /= ONE) THEN  ! Airfoil-type distribution

            CALL FOILGRD (N_OUTER, ZERO, ONE, ABS (RO_SPL), RO_SPQ, &
                          RO_SPS, RO_SPC, SNORM2)

            IF (RO_SPL < ZERO) THEN

               ALLOCATE (STEMP(N_OUTER))

               STEMP(:) = SNORM2(:)
               STOTAL = SNORM2(N_OUTER)

               DO I = 1, N_OUTER
                  SNORM2(I) = STOTAL - STEMP(N_OUTER - I + 1)
               END DO

               DEALLOCATE (STEMP)

            END IF

            IER = 0

         ELSE  ! Curvature-based distribution

            IF (RO_SPS == ZERO .AND. RO_SPC == ZERO) ISMOOTH = 0
            IF (RO_SPS /= ZERO .AND. RO_SPC == ZERO) ISMOOTH = 1
            IF (RO_SPS == ZERO .AND. RO_SPC /= ZERO) ISMOOTH = 2
            IF (RO_SPS /= ZERO .AND. RO_SPC /= ZERO) ISMOOTH = 3

            NTOT = NO_INIT(1)

            XR_SHIFT = RO1(1)      ! Shape preserving normalization
            XR_SCALE = RO1(NTOT) - RO1(1)
            RO1(:)   = (RO1(:) - XR_SHIFT) / XR_SCALE
            XR_SHIFT = XO1(1)
            XO1(:)   = (XO1(:) - XR_SHIFT) / XR_SCALE

            CALL CURVDIS (NTOT, XO1, RO1, N_OUTER, RO_SPQ, ISMOOTH, &
                          LWRIT, TRUE, SNORM2, SNORM2, IER)
            IF (IER /= 0) THEN
               WRITE (LWRIT, '(/, A, I5)') &
                  '*** Bad CURVDIS return (outer); IER:', IER
            END IF

            DEALLOCATE (XO1, RO1)
            STOTAL = SNORM2(N_OUTER)
            SNORM2(:) = SNORM2(:) / STOTAL

         END IF

   9     IF (IER /= 0) STOP

      END IF

      IF (AXISYMMETRIC) THEN  ! Avoid a bunch of stuff

         CALL PERTURB_SYMMETRIC ()  ! Internal procedure below

         GO TO 90

      END IF

!     --------------------------------------
!     Transfer input geometry to work space:
!     --------------------------------------

      XWCNTR = XICNTR
      YWCNTR = YICNTR
      ZWCNTR = ZICNTR
      XWRING = XIRING
      RWRING = RIRING
      XYZ_CG = XYZ_CG0

      IF (TWO_PARTS) THEN
         DO K = 1, NI_SEC
            CI_WORK(K) = CI_INIT(K)
            ZI_WORK(K) = ZI_INIT(K)
            DO I = 1, NI_INIT(K)
               XI_WORK(I,K) = XI_INIT(I,K)
               RI_WORK(I,K) = RI_INIT(I,K)
            END DO
         END DO
      END IF

      DO K = 1, NO_SEC
         CO_WORK(K) = CO_INIT(K)
         VO_WORK(K) = VO_INIT(K)
         ZO_WORK(K) = ZO_INIT(K)
         DO I = 1, NO_INIT(K)
            XO_WORK(I,K) = XO_INIT(I,K)
            RO_WORK(I,K) = RO_INIT(I,K)
         END DO
      END DO

!     --------------------------
!     Planform design variables:
!     --------------------------

      DO N = 1, NPLANF

         TYP4 = RBTYP(N)(1:4)

         VS = VPAR(N) * VSCALE(N)
         IF (ABS (VS) < SMALL) CYCLE ! Next N

         TYP3 = RBTYP(N)(2:4)

         SELECT CASE (TYP4)

            CASE ('XCEN')
               XWCNTR = XWCNTR + VS
               XWRING = XWRING + VS  ! Keep mount ring at fixed dist. from apex

               IF (FIXED_CG) XYZ_CG(1) = XYZ_CG(1) + VS ! Likewise unless Xcg is
                                                        ! an optimization var.
            CASE ('YCEN')
               YWCNTR = YWCNTR + VS

!!!!           IF (FIXED_CG) XYZ_CG(2) = XYZ_CG(2) + VS ! Just X for now

!!!!        CASE ('ZCEN')  ! Not allowed for now
!!!!           ZWCNTR = ZWCNTR + VS
!!!!           XYZ_CG(3) = XYZ_CG(3) + VS

            CASE ('X_CG')
               XYZ_CG(1) = XYZ_CG(1) + VS

            CASE ('Y_CG')
               XYZ_CG(2) = XYZ_CG(2) + VS

            CASE ('XRNG')
               XWRING = XWRING + VS

            CASE DEFAULT

               IF (TYP3 == 'OUT' ) THEN

                  DO K = 1, NO_SEC

                     CALL LOFT (ZBTYP(N), ZO_WORK(K), ZBUMP(N), ZBMIN(N),      &
                                ZBMAX(N), ZBEXP(N), CIRCLE, DELTAV, LWRIT)

                     DELTAV = DELTAV * VS

                     IF (ABS (DELTAV) > SMALL) THEN

                        SELECT CASE (TYP4)
                           CASE ('ROUT')
                              CO_WORK(K) = CO_WORK(K) + DELTAV
                           CASE ('XOUT')
                              VO_WORK(K) = VO_WORK(K) + DELTAV
                           CASE ('TOUT')
                              CO_WORK(K) = CO_WORK(K) + DELTAV
                              VO_WORK(K) = VO_WORK(K) + DELTAV * RBEXP(N)
                        END SELECT

                     END IF

                  END DO ! Next azimuthal defining station

               END IF

         END SELECT

      END DO ! Next planform design variable

!     -------------------------------
!     Defining section perturbations:
!     -------------------------------

!     Apply any Hicks-Henne shape functions to the heat shield sections:

      DO N = NPLANF + 1, NPLANF + NBUMPS

         VS = VPAR(N) * VSCALE(N)
         IF (ABS (VS) < SMALL) CYCLE ! Next N

         TYP1 = RBTYP(N)(1:1)
         TYP4 = RBTYP(N)(2:5)

!        -------------------------------------------------------
!        Perturb the inner shield sections (normalized/sheared):
!        -------------------------------------------------------

         IF (TWO_PARTS) THEN

            IF (TYP1 == 'I') THEN

               DO K = 1, NI_SEC

                  CALL LOFT (ZBTYP(N), ZI_WORK(K), ZBUMP(N), ZBMIN(N), &
                             ZBMAX(N), ZBEXP(N), CIRCLE, DELTAV, LWRIT)

                  DELTAV = DELTAV * VS

                  IF (ABS (DELTAV) > SMALL) THEN

                     NTOT =  NI_INIT(K)

                     CALL SFEVAL (TYP4, RBEXP(N), RBUMP(N), &
                                  RBMIN(N), RBMAX(N), DELTAV, 1, NTOT, &
                                  RI_WORK(1,K), XI_WORK(1,K), LWRIT)
                  END IF

               END DO

            END IF

         ELSE

!           ----------------------------------
!           Perturb the outer shield sections:
!           ----------------------------------

            DO K = 1, NO_SEC

               CALL LOFT (ZBTYP(N), ZO_WORK(K), ZBUMP(N), ZBMIN(N), &
                          ZBMAX(N), ZBEXP(N), CIRCLE, DELTAV, LWRIT)

               DELTAV = DELTAV * VS

               IF (ABS (DELTAV) > SMALL) THEN

                  NTOT = NO_INIT(K)

                  CALL SFEVAL (TYP4, RBEXP(N), RBUMP(N), &
                               RBMIN(N), RBMAX(N), DELTAV, 1, NTOT, &
                               RO_WORK(1,K), XO_WORK(1,K), LWRIT)
               END IF

            END DO

         END IF

      END DO ! Next design variable

!     -------------------------------------------------------
!     Denormalize r, unshear x, and regularize in real space:
!     -------------------------------------------------------

      CRVI_MAX   = ZERO
      CRVJ_MAX   = ZERO
      CRVG_MAX   = ZERO
      CRVI_TOT   = ZERO
      CRVJ_TOT   = ZERO
      CRVG_TOT   = ZERO
      X_RESIDUAL = ZERO

      IF (TWO_PARTS) THEN

         DX = XWRING - XWCNTR
         DY = YWCNTR - XYZ_CG(2)

         DO K = 1, NI_SEC ! For each inner shield defining section ...

            ZI_WORK(K) = ZI_WORK(K) * TOT_AZ

!           Calculate the radius from the CG location to the ring:

            CPH = COS_ZI(K)
            RK  = MAX (DY**2 * (CPH**2 - ONE) + RWRING**2, ZERO)
            RK  = -DY * CPH + SQRT (RK) ! Appropriate quadratic root
            CI_WORK(K) = RK             ! Used by the outer spokes

            NTOT = NI_INIT(K)

            DO I = 1, NTOT
               RN = RI_WORK(I,K)
               RI_WORK(I,K) = RN * RK                         ! Denormalize
               XI_WORK(I,K) = XI_WORK(I,K) + RN * DX + XWCNTR ! Unshear
            END DO

!           Option to smooth out a possible cusp at the apex (in-place):

            IF (CUSP > SMALL) &
               CALL NOCUSP (CUSP, NTOT, RI_WORK(1,K), XI_WORK(1,K))

!           Arc lengths are needed both for curvature calculations and
!           for regularizing the sections:

            CALL CHORDS2D (NTOT, RI_WORK(1,K), XI_WORK(1,K), NORM, SMAX, SGEOM)

!           Accumulate the local curvature violations:

            NCRV_MAX = NTOT - 1  ! No longer a reason to avoid some spokes?

            CALL KAPPA (NTOT, NCRV_MAX, RI_WORK(1,K), XI_WORK(1,K), SGEOM,     &
                        CR_TARG, CRVI_MAX, CRVI_TOT)

!           Regularize the radial distribution:

            SGRID(1:N_INNER) = SMAX * SNORM1(1:N_INNER)

            CALL LCSFIT (NTOT, SGEOM, RI_WORK(1,K), NEW, TIGHT, &
                         N_INNER, SGRID, RRG(1,K), DERIVS)

            CALL LCSFIT (NTOT, SGEOM, XI_WORK(1,K), NEW, LOOSE, &
                         N_INNER, SGRID, XRG(1,K), DERIVS)

         END DO ! Next inner section

!        Option to smooth X vs. azimuthal angle in-place via best-fit quartics,
!        and/or to penalize the sum of the residuals over all the fits.
!        Careful: Azimuth arrays start at 0 for periodicity reasons.

         IF (SMOOTH_X .OR. CALC_X_RESIDUAL) THEN

            CALL SMOOTH_X_VS_PHI (NT_RAD, NI_SEC, ZI_WORK(1), XI_WORK, &
                                  SMOOTH_X, X_RESIDUAL)
         END IF

!        Loop over all outer defining sections:

         I1 = N_INNER

         DO K = 1, NO_SEC

            ZO_WORK(K) = ZO_WORK(K) * TOT_AZ

!           Determine the radius from the shifted center to the unshifted rim,
!           as defined by the (fixed) azimuthal common to inner & outer spokes:

            DY  = YWCNTR - YICNTR
            CPH = COS_ZO(K)
            RK  = CO_INIT(K)
            RK  = MAX (DY**2 * (CPH**2 - ONE) + RK**2, ZERO) ! If circle
            RK  = -DY * CPH + SQRT (RK) ! Appropriate quadratic root

            RK = RK * (CO_WORK(K) / CO_INIT(K)) ! Account for possible tab
            R0 = CI_WORK(K) ! How can this be avoided?
            DR = RK - R0
            CO_WORK(K) = RK ! Saved in the output geometry
            DX = VO_WORK(K) - XWRING

            DO I = 1, NO_INIT(K)
               RN = RO_WORK(I,K)
               RO_WORK(I,K) = R0 + RN * DR                    ! Denormalize
               XO_WORK(I,K) = XO_WORK(I,K) + RN * DX + XWRING ! Unshear
            END DO

         END DO ! Need to do all sections before possible smoothing of X

      ELSE ! One-part shield

         I1 = 1

!        Loop over all defining sections:

         DO K = 1, NO_SEC

            ZO_WORK(K) = ZO_WORK(K) * TOT_AZ

            NTOT = NO_INIT(K)

!           Allow for (only) vertical shifts of the apex.
!           Determine the radius from the shifted center to the unshifted rim,
!           as defined by the (fixed) azimuthal angle:

            DY = YWCNTR - YICNTR
            RK = CO_INIT(K)

            IF (ABS (DY) > SMALL) THEN ! Use the triangle cosine rule
               CPH = COS_ZO(K)
               RK  = MAX (DY**2 * (CPH**2 - ONE) + RK**2, ZERO) ! If circle
               RK  = -DY * CPH + SQRT (RK) ! Appropriate quadratic root
            END IF

            RK = RK * (CO_WORK(K) / CO_INIT(K)) ! Account for possible tab
            CO_WORK(K) = RK ! Saved in the output geometry
            DX = VO_WORK(K) - XWCNTR

            DO I = 1, NTOT
               RN = RO_WORK(I,K)
               RO_WORK(I,K) = RN * RK                         ! Denormalize
               XO_WORK(I,K) = XO_WORK(I,K) + RN * DX + XWCNTR ! Unshear
            END DO

!           Option to smooth a possible cusp at the apex (real space, in-place):

            IF (CUSP > SMALL) &
               CALL NOCUSP (CUSP, NTOT, RO_WORK(1,K), XO_WORK(1,K))

         END DO ! Need to do all sections before possible smoothing of X

!        Asymmetric cases do not necessarily have continuous curvature across
!        the apex, unless we control the second derivatives as well as the
!        slopes.  Even this doesn't guarantee SMOOTH curvature across the apex:

         IF (CUSP > SMALL) CALL CONTROL_APEX_CURVATURE ()  ! Local procedure

      END IF

!     Option to smooth X vs. azimuthal angle in-place via best-fit quartics,
!     and/or to penalize the sum of the residuals over all the fits.
!     Careful: Azimuth arrays start at 0 for periodicity reasons.

      IF (SMOOTH_X .OR. CALC_X_RESIDUAL) THEN

         CALL SMOOTH_X_VS_PHI (NT_RAD, NO_SEC, ZO_WORK(1), XO_WORK, &
                               SMOOTH_X, X_RESIDUAL)
      END IF


!     Regularize the radial distribution at the specified grid resolution:

      DO K = 1, NO_SEC

         NTOT = NO_INIT(K)

         CALL CHORDS2D (NTOT, RO_WORK(1,K), XO_WORK(1,K), NORM, SMAX, SGEOM)

!        Accumulate the local curvature violations and find peak curvature.
!        Avoid the shoulder region, and confine search to lee side.

         NCRV_MAX = MIN (NINT (REAL (NTOT) * FRACTION_CRVIMX), NTOT - 1)

         CALL KAPPA (NTOT, NCRV_MAX, RO_WORK(1,K), XO_WORK(1,K), SGEOM,        &
                     CR_TARG, CRVI_MAX, CRVI_TOT)

!        Redistribute in the radial direction:

         SGRID(1:N_OUTER) = SMAX * SNORM2(1:N_OUTER)

         CALL LCSFIT (NTOT, SGEOM, RO_WORK(1,K), NEW, TIGHT, &
                      N_OUTER, SGRID, RRG(I1,K), DERIVS)

         CALL LCSFIT (NTOT, SGEOM, XO_WORK(1,K), NEW, LOOSE, &
                      N_OUTER, SGRID, XRG(I1,K), DERIVS)

      END DO ! Next outer section

!     ------------------------------------------------------
!     Interpolate radial sections in the azimuthal direction
!     (denormalized polar coordinates).  Calculate curvature
!     violations while we're at it.
!     ------------------------------------------------------

!     Common apex point for all cases:

      DO K = 1, NK
         RFIN(1,K) = RRG(1,1)
         XFIN(1,K) = XRG(1,1)
      END DO

      IF (TWO_PARTS) THEN ! Inner shield first

         IF (NI_SEC == 1) THEN ! Special case of only one defining section.

            DO K = 1, NK
               DO I = 2, N_INNER
                  RFIN(I,K) = RRG(I,1)
                  XFIN(I,K) = XRG(I,1)  ! CRVJ_TOT remains zero
               END DO
            END DO

         ELSE

           CALL SURFACE_GRID (FULL, NI_SEC, NI, NK, 2, N_INNER,   &
                              RRG, XRG, ZI_WORK, RFIN, XFIN, AZM, &
                              YWCNTR, ZWCNTR, CRVJ_TARG, CRVJ_MAX, CRVJ_TOT)
         END IF

      END IF

!     Azimuthal redistribution for the outer shield:

      IF (NO_SEC == 1) THEN ! Special case of only one defining section.

         DO K = 1, NK
            DO I = I1, NI
               RFIN(I,K) = RRG(I,1)
               XFIN(I,K) = XRG(I,1)
            END DO
         END DO

      ELSE ! Interpolate in the azimuthal direction.

         CALL SURFACE_GRID (FULL, NO_SEC, NI, NK, I1 + 1, NI,   &
                            RRG, XRG, ZO_WORK, RFIN, XFIN, AZM, &
                            YWCNTR, ZWCNTR, CRVJ_TARG, CRVJ_MAX, CRVJ_TOT)

      END IF

!     ---------------------------------------------------------
!     Complete the shield construction by converting from real-
!     space cylindrical coordinates to Cartesian coordinates:
!     ---------------------------------------------------------

      DO K = 1, NK
         CPH = COSPHI(K)
         SPH = SINPHI(K)

         DO I = 1, NI
            R = RFIN(I,K) - RFIN(1,K)
            YFIN(I,K) = YWCNTR + R * CPH
            ZFIN(I,K) = ZWCNTR + R * SPH
         END DO

      END DO

 90   CONTINUE  ! Common to symmetric case

!     Calculate the cross-sectional area and volume. V = integral of X = X(Y,Z).
!     The following is legitimate as long as the spoke radii increase steadily.

      DO K = 1, NK
         DX = XFIN(NI,K)
         XLEN(:,K) = DX - XFIN(:,K)  ! X increases from the apex
      END DO

      CALL TRAPEZOIDAL_QUADRATURE_2D (NI, NK, 1, YFIN, ZFIN, XLEN, &
                                      SECTIONAL_AREA, VOLUME)

!     Calculate the Gaussian and mean curvature distributions for the structured
!     surface grid.  Note that the 2-space curvature constraints work with the
!     geometry data, and while there's no reason for the spokes not to form a
!     structured surface, they tend to be far apart (typically 11 of them for
!     half the forebody).  Ensuring these curvatures have the right sign for all
!     computational grid cells errs on the safe side.

      CALL GAUSSIAN_CURVATURE (NI, NK, XFIN, YFIN, ZFIN, &
                               GAUSS_CURVATURE, MEAN_CURVATURE)

!!!   write (20, *) 1
!!!   write (20, *) NI, NK, 1
!!!   write (20, *) XFIN, YFIN, ZFIN

!!!   write (21, *) 1
!!!   write (21, *) NI, NK, 1, 2
!!!   write (21, *) GAUSS_CURVATURE, MEAN_CURVATURE

!     Gaussian curvature is positive for convex right-handed surfaces.
!     Calculate the largest value and the sum of (squared) violations of the
!     minimum desirable curvature (zero by default or entered via XNLCON(:)):

      TVD_GC = ZERO

      DO K = 1, NK
         GAUSSIM1 = GAUSS_CURVATURE(1,K)
         DO I = 2, NI  ! Avoid collapsed cells at nose for CRVG_TOT
            GAUSSI   = GAUSS_CURVATURE(I,K)
            CRVG_MAX = MAX (CRVG_MAX, GAUSSI)
            CRVG_TOT = CRVG_TOT + (MIN (GAUSSI - CRVG_TARG, ZERO)) ** 2
            TVD_GC   = TVD_GC   + (GAUSSI - GAUSSIM1) ** 2
            GAUSSIM1 = GAUSSI
         END DO
      END DO

!     This contribution to TVD_GC is probably not a good idea for asymmetric
!     cases because the grid lines don't follow the varying shoulder.

!!!   DO I = 1, NI
!!!      GAUSSKM1 = GAUSS_CURVATURE(I,1)
!!!      DO K = 2, NK
!!!         GAUSSK   = GAUSS_CURVATURE(I,K)
!!!         TVD_GC   = TVD_GC   + (GAUSSK - GAUSSKM1) ** 2
!!!         GAUSSKM1 = GAUSSK
!!!      END DO
!!!   END DO

!     Internal procedure for subroutine PERTURB:

      CONTAINS

!        -----------------------------------------------------------------------

         SUBROUTINE CONTROL_APEX_CURVATURE ()

!        Ensure continuous (but still not necessarily smooth) curvature across
!        the apex for all pairs of opposing defining spokes (assumed to be odd
!        in number and uniformly distributed azimuthally).
!        Make this argument-driven if the TWO_PARTS option is ever revived.

!        -----------------------------------------------------------------------

!        Local variables:

         INTEGER       :: I, J, K, L, NSPOKE

         INTEGER, SAVE :: KL       ! One less than the middle spoke
         INTEGER, SAVE :: NI       ! # pts. on all defining spokes
         INTEGER, SAVE :: M = 999  ! As in NOCUSP

         REAL          :: H, XP1, XPP1, XPPAVG, XP2, XPP2, XPPP2
         REAL          :: C(6)

!        Execution:

         NSPOKE = NO_SEC

         IF (M == 999) THEN
            NI = NO_INIT(1)   ! Assume spokes have equal geometry point counts
            M  = NI
            DO I = 3, NI
               IF (RO_WORK(I,1) > CUSP) THEN
                  M = I
                  EXIT
               END IF
            END DO
            KL = (NSPOKE - 1) / 2
         END IF

         DO K = 1, KL

            CALL FD1SID (1, 1, RO_WORK(1,K), XO_WORK(1,K), XP1, XPP1)

            XPPAVG = XPP1

            L = NSPOKE - (K - 1)  ! Opposing spoke index

            CALL FD1SID (1, 1, RO_WORK(1,L), XO_WORK(1,L), XP1, XPP1)

            XPPAVG = (XPPAVG + XPP1) * 0.5

            DO J = K, L, L - K  ! For both opposing spokes, reinterpolate nose

               CALL DERIVS_123 (NI, RO_WORK(1,J), XO_WORK(1,J), M, &
                                XP2, XPP2, XPPP2)

               CALL POLY6 (RO_WORK(1,J), XO_WORK(1,J), ZERO, XPPAVG, &
                           RO_WORK(M,J), XO_WORK(M,J), XP2, XPP2, XPPP2, C)

               DO I = 2, M - 1
                  H = RO_WORK(I,J)  ! - zero
                  XO_WORK(I,J) = XO_WORK(1,J) + &
                   H*(C(1) + H*(C(2) + H*(C(3) + H*(C(4) + H*(C(5) + H*C(6))))))
               END DO

            END DO

         END DO

         END SUBROUTINE CONTROL_APEX_CURVATURE

!        -----------------------------------------------------------------------

         SUBROUTINE PERTURB_SYMMETRIC ()

!        Simplified version of PERTURB for the single-spoke axisymmetric case.
!        Use the host variables - no need to declare local ones.

!        -----------------------------------------------------------------------

         XWCNTR = XICNTR  ! YWCNTR and ZWCNTR can't move
         XYZ_CG = XYZ_CG0

         CO_WORK(1) = CO_INIT(1)
         VO_WORK(1) = VO_INIT(1)

         NTOT = NO_INIT(1)

         DO I = 1, NTOT
            XO_WORK(I,1) = XO_INIT(I,1)
            RO_WORK(I,1) = RO_INIT(I,1)
         END DO

!        --------------------------
!        Planform design variables:
!        --------------------------

         DO N = 1, NPLANF

            TYP4 = RBTYP(N)(1:4)

            VS = VPAR(N) * VSCALE(N)
            IF (ABS (VS) < SMALL) CYCLE ! Next N

            TYP3 = RBTYP(N)(2:4)

            SELECT CASE (TYP4)

               CASE ('XCEN')
                  XWCNTR = XWCNTR + VS

                  IF (FIXED_CG) XYZ_CG(1) = XYZ_CG(1) + VS

               CASE ('X_CG')
                  XYZ_CG(1) = XYZ_CG(1) + VS

               CASE ('Y_CG')
                  XYZ_CG(2) = XYZ_CG(2) + VS

            END SELECT

         END DO ! Next planform design variable

!        -------------------------------------------------------------
!        Defining section perturbations (Hicks-Henne shape functions):
!        -------------------------------------------------------------

         DO N = NPLANF + 1, NPLANF + NBUMPS

            VS = VPAR(N) * VSCALE(N)
            IF (ABS (VS) < SMALL) CYCLE ! Next N

            TYP1 = RBTYP(N)(1:1)
            TYP4 = RBTYP(N)(2:5)

            CALL SFEVAL (TYP4, RBEXP(N), RBUMP(N), &
                               RBMIN(N), RBMAX(N), VS, 1, NTOT, &
                               RO_WORK(1,1), XO_WORK(1,1), LWRIT)
         END DO

!        -------------------------------------------------------
!        Denormalize r, unshear x, and regularize in real space:
!        -------------------------------------------------------

         CRVI_MAX   = ZERO
         CRVJ_MAX   = ZERO
         CRVI_TOT   = ZERO
         CRVJ_TOT   = ZERO
         X_RESIDUAL = ZERO

         RK = CO_INIT(1)
         DX = VO_WORK(1) - XWCNTR

         DO I = 1, NTOT
            RN = RO_WORK(I,1)
            RO_WORK(I,1) = RN * RK                         ! Denormalize
            XO_WORK(I,1) = XO_WORK(I,1) + RN * DX + XWCNTR ! Unshear
         END DO

!        Option to smooth a possible cusp at the apex (real space, in-place):

         IF (CUSP > SMALL) CALL NOCUSP (CUSP, NTOT, RO_WORK(1,1), XO_WORK(1,1))

!        Impose the computational grid in the radial direction:

         CALL CHORDS2D (NTOT, RO_WORK(1,1), XO_WORK(1,1), NORM, SMAX, SGEOM)

!        Accumulate the local radial curvature violations & find peak curvature.
!        Avoid the shoulder region.

         NCRV_MAX = MIN (NINT (REAL (NTOT) * FRACTION_CRVIMX), NTOT - 1)

         CALL KAPPA (NTOT, NCRV_MAX, RO_WORK(1,1), XO_WORK(1,1), SGEOM,        &
                     CR_TARG, CRVI_MAX, CRVI_TOT)

         SGRID(1:N_OUTER) = SMAX * SNORM2(1:N_OUTER)

         CALL LCSFIT (NTOT, SGEOM, RO_WORK(1,1), NEW, TIGHT, &
                      N_OUTER, SGRID, RRG(1,1), DERIVS)

         CALL LCSFIT (NTOT, SGEOM, XO_WORK(1,1), NEW, LOOSE, &
                      N_OUTER, SGRID, XRG(1,1), DERIVS)

!!!      write (lwrit, '(a)') 'i, snorm2, sgrid, xrg, rrg'
!!!      write (lwrit, '(i3, 4f12.7)') &
!!!         (i, snorm2(i), sgrid(i), xrg(i,1), rrg(i,1), i = 1, n_outer)

!        Convert from cylindrical to Cartesian coordinates, axisymmetric:

         DO K = 1, NK
            CPH = COSPHI(K)
            SPH = SINPHI(K)

            DO I = 1, NI
               R = RRG(I,1) - RRG(1,1)
               YFIN(I,K) = YICNTR + R * CPH
               ZFIN(I,K) = ZICNTR + R * SPH
               XFIN(I,K) = XRG(I,1)
            END DO
         END DO

         END SUBROUTINE PERTURB_SYMMETRIC

      END SUBROUTINE PERTURB

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE POST_AERO
!
!     Post-process aerodynamics and prepare for objective function calculation.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Global variables:

      USE AERO_MOD
      USE COEFFS_MOD
      USE CONST_MOD

      IMPLICIT NONE

!     Execution:

      IF (ABS (CD) > SMALL) THEN
         LD = CL / CD
      ELSE
         LD = 1000.
      END IF

      END SUBROUTINE POST_AERO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE RD_CONS
!
!     Read the constraints.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      USE COEFFS_MOD
      USE CONST_MOD
      USE DESIGN_VARS_MOD
      USE GRID_MOD
      USE OPT_MOD

      IMPLICIT NONE

!     Local constants:

      REAL, PARAMETER :: ZERO = 0.

!     Local variables:

      INTEGER :: J, K, M, IOS
      REAL    :: BOUNDLO, BOUNDUP

!     Execution:
!     ----------

!     Linear constraints:

      READ (LREAD, *)
      READ (LREAD, *)
      READ (LREAD, *)

      DO M = 1, NCLIN

         READ (LREAD, *, IOSTAT=IOS) J, LCTYPE(M), BOUNDLO, BOUNDUP, &
                                     TLCON(M), ILCON(M)
         IF (IOS /= 0) THEN
            WRITE (LWRIT, *) 'Error reading linear constraints.'
            GO TO 999
         END IF

         IF (INT (BOUNDLO) == -999) BOUNDLO = -BIGBND
         IF (INT (BOUNDUP) ==  999) BOUNDUP =  BIGBND
         BLLIN(M) = BOUNDLO
         BULIN(M) = BOUNDUP

      END DO

!     Nonlinear constraints:

      CALC_X_RESIDUAL = RHO_SMTH_X > EPSILON (RHO_SMTH_X)
      CRVG_TARG       = ZERO   ! In case Gaussian curvature constraint is absent
      CRVJ_TARG       = ZERO   ! In case azimuth. curvature constraint is absent
      NEED_DCMDA      = FALSE

      READ (LREAD, *)
      READ (LREAD, *)
      READ (LREAD, *)

      K = NDV + NCLIN ! Offset for the bounds

      DO M = 1, NCNLN

         READ (LREAD, *, IOSTAT=IOS) J, NLCTYPE(M), BOUNDLO, BOUNDUP, &
            XNLCON(M), INLCON(M), JNLCON(M), SNLCON(M)

         IF (IOS /= 0) THEN
            WRITE (LWRIT, *) J, NLCTYPE(M), BOUNDLO, BOUNDUP, &
               XNLCON(M), INLCON(M), JNLCON(M), SNLCON(M)
            WRITE (LWRIT, *) 'IOS: ', IOS
            WRITE (LWRIT, *) 'Error reading nonlinear constraints.'
            GO TO 999
         END IF

!        Special action is needed for some constraints:

         SELECT CASE (NLCTYPE(M))

         CASE ('CM_ALF')
            NEED_DCMDA      = TRUE
         CASE ('CRVIMX')
            FRACTION_CRVIMX = XNLCON(M)
         CASE ('CURVG ')
            CRVG_TARG       = XNLCON(M)
         CASE ('CURVI ')
            XNLCON(M)       = CR_TARG   ! CR_TARG is reqd. for obj. fn. form too
         CASE ('CURVJ ')
            CRVJ_TARG       = XNLCON(M) ! Don't use objective function form
         CASE ('SMTH_X')
            CALC_X_RESIDUAL = TRUE
         CASE DEFAULT
            CONTINUE
         END SELECT

!        Scaling of the nonlinear constraints differs from scaling of
!        the variables and linear constraints:

         IF (INT (BOUNDLO) == -999) THEN
            BOUNDLO = -BIGBND
         ELSE
            BOUNDLO = BOUNDLO * SNLCON(M)
         END IF

         IF (INT (BOUNDUP) ==  999) THEN
            BOUNDUP = BIGBND
         ELSE
            BOUNDUP = BOUNDUP * SNLCON(M)
         END IF

         BU(M+K) = BOUNDUP
         BL(M+K) = BOUNDLO

      END DO

      RETURN

 999  CONTINUE

      STOP

      END SUBROUTINE RD_CONS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE RD_GEOM (ANALYTIC)
!
!     Read a heat shield geometry defined in cylindrical coordinates as
!     normalized or unnormalized sections.
!     NON_DIM == 1 means normalize/shear the sections here.
!     A single input spoke implies axisymmetry.  Assume this applies to the
!     single-part case only.
!     If heat_shield.analytic has been opened successfully, construct a
!     generatrix rather than read it.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      USE AERO_MOD
      USE COEFFS_MOD
      USE CONST_MOD
      USE GEOM_MOD
      USE BICONIC_PARAMETERS      ! For derived data type biconic_type
      USE SPHERE_CONE_PARAMETERS  ! For derived data type sphere_cone_type
      USE SPHERE_CONE             ! For constructing a generatrix

      IMPLICIT NONE

!     Arguments:

      LOGICAL, INTENT (IN) :: ANALYTIC  ! T => construct/don't read generatrix

!     Local constants:

      REAL, PARAMETER :: HALF = 0.5, ONE = 1., ZERO = 0.

!     Local variables:

      INTEGER :: I, IOS, K, NGEOM, NLINES, NUMI
      REAL    :: CPH, DR, DX, DY, H, RINVRS, RK, RNORM, RRIM

      TYPE (BICONIC_TYPE)     :: BICONIC
      TYPE (SPHERE_CONE_TYPE) :: CAPSULE

!     Execution:

      IF (ANALYTIC) THEN

!        Count the number of lines to distinguish biconic from sphere-cone:

         NLINES = 0
         DO  ! Until EOF
            READ (LGEOM, *, IOSTAT=IOS)
            IF (IOS < 0) EXIT
            NLINES = NLINES + 1
         END DO

         REWIND (LGEOM)

         SPHERECONE = NLINES == 10  ! Else 12 for the biconic

         IF (SPHERECONE) THEN
            READ (LGEOM, *)
            READ (LGEOM, *) X_NOSE;          CAPSULE%X_NOSE       = X_NOSE
            READ (LGEOM, *) R_NOSE;          CAPSULE%R_NOSE       = R_NOSE
            READ (LGEOM, *) RADIUS_NOSE;     CAPSULE%RADIUS_NOSE  = RADIUS_NOSE
            READ (LGEOM, *) RADIUS_BASE;     CAPSULE%RADIUS_BASE  = RADIUS_BASE
            READ (LGEOM, *) RADIUS_SHOULDER; CAPSULE%RADIUS_SHOULDER = &
                            RADIUS_SHOULDER
            READ (LGEOM, *) HALF_CONE_ANGLE; CAPSULE%HALF_CONE_ANGLE = &
                            HALF_CONE_ANGLE
            READ (LGEOM, *) SKIRT_ANGLE;     CAPSULE%SKIRT_ANGLE  = SKIRT_ANGLE
            READ (LGEOM, *) SKIRT_LENGTH;    CAPSULE%SKIRT_LENGTH = SKIRT_LENGTH
            READ (LGEOM, *) NGEOM
         ELSE
            READ (LGEOM, *)
            READ (LGEOM, *) X_NOSE;          BICONIC%X_NOSE       = X_NOSE
            READ (LGEOM, *) R_NOSE;          BICONIC%R_NOSE       = R_NOSE
            READ (LGEOM, *) RADIUS_NOSE;     BICONIC%RADIUS_NOSE  = RADIUS_NOSE
            READ (LGEOM, *) RADIUS_CONE_JUNCTURE
                            BICONIC%RADIUS_CONE_JUNCTURE  = RADIUS_CONE_JUNCTURE
            READ (LGEOM, *) RADIUS_BASE;     BICONIC%RADIUS_BASE  = RADIUS_BASE
            READ (LGEOM, *) RADIUS_SHOULDER; BICONIC%RADIUS_SHOULDER = &
                            RADIUS_SHOULDER
            READ (LGEOM, *) HALF_CONE_ANGLE_FORE
                            BICONIC%HALF_CONE_ANGLE_FORE = HALF_CONE_ANGLE_FORE
            READ (LGEOM, *) HALF_CONE_ANGLE_AFT
                            BICONIC%HALF_CONE_ANGLE_AFT  = HALF_CONE_ANGLE_AFT
            READ (LGEOM, *) SKIRT_ANGLE;     BICONIC%SKIRT_ANGLE  = SKIRT_ANGLE
            READ (LGEOM, *) SKIRT_LENGTH;    BICONIC%SKIRT_LENGTH = SKIRT_LENGTH
            READ (LGEOM, *) NGEOM
         END IF

         REF_AREA     = PI * RADIUS_BASE**2
         REF_LEN      = RADIUS_BASE
         AXISYMMETRIC = TRUE
         FULL         = 0
         NI_SEC       = 1
         NO_SEC       = 1
         NON_DIM      = 1
         XICNTR       = X_NOSE
         YICNTR       = R_NOSE
         ZICNTR       = ZERO

      ELSE

         READ (LGEOM, *)
         READ (LGEOM, *)
         READ (LGEOM, *) REF_AREA, REF_LEN, XYZ_CG0(1:3), SBASE

         CUSP = CUSP * REF_LEN

         READ (LGEOM, *)
         READ (LGEOM, *) FULL, NI_SEC, NO_SEC, NON_DIM

         IF (.NOT. TWO_PARTS) NI_SEC = 1  ! Let associated allocates be minimal

         AXISYMMETRIC = NO_SEC == 1

         READ (LGEOM, *)
         READ (LGEOM, *) XICNTR, YICNTR, ZICNTR, XIRING, RIRING

      END IF

!     Set up to calculate either a half or whole surface:

      IF (FULL == 1) THEN  ! Why would we ever use this?!
         S_REF  = REF_AREA
         TOT_AZ = TWOPI
         CIRCLE = ONE ! Used in LOFT
      ELSE ! Half config.
         SIXROOTPI = SIXROOTPI / SQRT (2.)  ! For effective volume coefficient
         S_REF  = REF_AREA * HALF
         TOT_AZ = PI
         CIRCLE = HALF
      END IF

      RTOT_AZ = ONE / TOT_AZ
      NT_SEC  = MAX (NI_SEC, NO_SEC)

      ALLOCATE (NI_INIT(NI_SEC), NO_INIT(NO_SEC),                              &
                ZI_INIT(NI_SEC), ZO_INIT(NO_SEC),                              &
                CI_INIT(NI_SEC), CO_INIT(NO_SEC), VO_INIT(NO_SEC),             &
                COS_ZI (NI_SEC), COS_ZO (NO_SEC),                              &
                SIN_ZI (NI_SEC), SIN_ZO (NO_SEC))

      IF (AXISYMMETRIC) THEN
         ALLOCATE (ZO_WORK(1), CO_WORK(1), VO_WORK(1))
      ELSE
         ALLOCATE (ZI_WORK(0:NI_SEC+1), ZO_WORK(0:NO_SEC+1),   & ! For periodic
                   CI_WORK(0:NI_SEC+1), CO_WORK(0:NO_SEC+1),   & ! interpolation
                   VO_WORK(0:NO_SEC+1))
      END IF

!     Either generate or read the defining section(s):

      IF (ANALYTIC) THEN

         ALLOCATE (RO_INIT(NGEOM,1), XO_INIT(NGEOM,1), &
                   RO_WORK(NGEOM,1), XO_WORK(NGEOM,1))

         IF (SPHERECONE) THEN
            CALL SPHERE_CONE_CONTROL_POINTS (CAPSULE)
            CALL SPHERE_CONE_DISCRETIZATION (CAPSULE, NGEOM, XO_INIT, RO_INIT)
         ELSE
            CALL BICONIC_CONTROL_POINTS (BICONIC)
            CALL BICONIC_DISCRETIZATION (BICONIC, NGEOM, XO_INIT, RO_INIT)

            X_CONE_JUNCTURE  = BICONIC%X_CONE_JUNCTURE  ! In case we need to
            R_CONE_JUNCTURE  = BICONIC%R_CONE_JUNCTURE  ! capture it in the
            CONE_LENGTH_FORE = BICONIC%CONE_LENGTH_FORE ! surface grid
            CONE_LENGTH_AFT  = BICONIC%CONE_LENGTH_AFT
         END IF

         NT_RAD     = NGEOM
         NO_INIT(1) = NGEOM
         ZO_INIT(1) = ZERO
         CO_INIT(1) = RO_INIT(NGEOM,1)
         VO_INIT(1) = XO_INIT(NGEOM,1)
         COS_ZO(1)  = ONE
         SIN_ZO(1)  = ZERO
         XYZ_CG0(1) = XO_INIT(1,1) + 0.25*(XO_INIT(NGEOM,1) - XO_INIT(1,1))
         XYZ_CG0(2) = RO_INIT(1,1) - 0.04*RADIUS_BASE
         XYZ_CG0(3) = ZERO                          ! The CG is a kludge

      ELSE

!        Scan the defining sections for the largest number of points,
!        allocate storage, then read the sections:

         NT_RAD = 0

         IF (TWO_PARTS) THEN

            READ (LGEOM, *)
            DO K = 1, NI_SEC
               READ (LGEOM, *)
               READ (LGEOM, *)
               READ (LGEOM, *) NI_INIT(K), ZI_INIT(K), CI_INIT(K)
               NT_RAD = MAX (NT_RAD, NI_INIT(K))
               READ (LGEOM, *)
               DO I = 1, NI_INIT(K)
                  READ (LGEOM, *)
               END DO
            END DO

            ALLOCATE (RI_INIT(NT_RAD,NI_SEC), XI_INIT(NT_RAD,NI_SEC), &
                      RI_WORK(NT_RAD,NI_SEC), XI_WORK(NT_RAD,NI_SEC))
         END IF

         READ (LGEOM, *)

         DO K = 1, NO_SEC
            READ (LGEOM, *)
            READ (LGEOM, *)
            READ (LGEOM, *) NO_INIT(K), ZO_INIT(K), CO_INIT(K), VO_INIT(K)
            NT_RAD = MAX (NT_RAD, NO_INIT(K))
            READ (LGEOM, *)
            DO I = 1, NO_INIT(K)
               READ (LGEOM, *)
            END DO
         END DO

         ALLOCATE (RO_INIT(NT_RAD,NO_SEC), XO_INIT(NT_RAD,NO_SEC), &
                   RO_WORK(NT_RAD,NO_SEC), XO_WORK(NT_RAD,NO_SEC))

         REWIND (LGEOM)

         DO I = 1, 8
            READ (LGEOM, *)
         END DO

         IF (TWO_PARTS) THEN

            DO K = 1, NI_SEC
               COS_ZI(K) = COS (ZI_INIT(K) * TOT_AZ)
               SIN_ZI(K) = SIN (ZI_INIT(K) * TOT_AZ)
               DO I = 1, 4
                  READ (LGEOM, *)
               END DO
               READ (LGEOM, *) (RI_INIT(I,K), XI_INIT(I,K), I = 1, NI_INIT(K))
            END DO

!           For the curvature-based spoke distribution option, save the first
!           defining section in real space before it is sheared:

            NUMI = NI_INIT(1)
            ALLOCATE (RI1(NUMI), XI1(NUMI))

            RI1(:) = RI_INIT(1:NUMI,1)
            XI1(:) = XI_INIT(1:NUMI,1)

            IF (NON_DIM == 1) THEN ! Normalize and shear the sections

!              Allow for a Ycg that has been changed from the value that
!              the input spokes correspond to.  The outer spokes need to
!              use the "chords" (max. radii) of the inner spokes.  We
!              update the chords, but normalize the r's with the input chords.

               DX = XIRING - XICNTR
               DY = YICNTR - XYZ_CG0(2) ! ? Not reconciled with X_CG/Y_CG vars.

               DO K = 1, NI_SEC
                  CPH = COS_ZI(K)
                  RK  = MAX (DY**2 * (CPH**2 - ONE) + RIRING**2, ZERO)
                  RK  = -DY * CPH + SQRT (RK) ! Appropriate quadratic root
                                              ! gives radius from apex to ring
                  RINVRS = ONE / CI_INIT(K)
                  CI_INIT(K) = RK             ! To be used for denormalizing
                                              ! and for "LE" of outer spoke
                  DO I = 1, NI_INIT(K)
                     RNORM        = RI_INIT(I,K) * RINVRS
                     RI_INIT(I,K) = RNORM
                     XI_INIT(I,K) = XI_INIT(I,K) - RNORM * DX - XICNTR
                  END DO
               END DO
            END IF

            DO K = 1, NI_SEC

               IF (NI_INIT(K) /= NUMI) GO TO 900

               IF (ABS (XI_INIT(1,K))          > SMALL .OR. &
                   ABS (XI_INIT(NI_INIT(K),K)) > SMALL) THEN
                  WRITE (LWRIT, *) &
                     ' Input section did not begin and end with X = 0.', &
                     ' Inner K: ', K, ' First & last normalized Xs: ',   &
                     XI_INIT(1,K), XI_INIT(NI_INIT(K),K)
                  STOP
               END IF

               IF (ABS (RI_INIT(1,K))                > SMALL .OR. &
                   ABS (RI_INIT(NI_INIT(K),K) - ONE) > SMALL) THEN
                  WRITE (LWRIT, *) &
                     ' Input section did not begin with R = 0',   &
                     ' and end with R = 1.',                      &
                     ' Inner K: ', K, ' First & last R(norm)s: ', &
                     RI_INIT(1,K), RI_INIT(NI_INIT(K),K)
                  STOP
               END IF

            END DO

            READ (LGEOM, *)
         END IF

         DO K = 1, NO_SEC
            COS_ZO(K) = COS (ZO_INIT(K) * TOT_AZ)
            SIN_ZO(K) = SIN (ZO_INIT(K) * TOT_AZ)
            DO I = 1, 4
               READ (LGEOM, *)
            END DO
            READ (LGEOM, *) (RO_INIT(I,K), XO_INIT(I,K), I = 1, NO_INIT(K))
         END DO

      END IF  ! End of non-analytic branch

      CLOSE (LGEOM)

      NUMI = NO_INIT(1)
      ALLOCATE (RO1(NUMI), XO1(NUMI))

      RO1(:) = RO_INIT(1:NUMI,1)  ! Save initial section for use by CURVDIS
      XO1(:) = XO_INIT(1:NUMI,1)

      IF (NON_DIM == 1) THEN

         DO K = 1, NO_SEC

            IF (TWO_PARTS) THEN

               IF (ZO_INIT(K) /= ZI_INIT(K)) THEN
                  WRITE (LWRIT, *) &
                     'Inner & outer sections must be at same azimuths.'
                  STOP
               END IF

               RK     = RO_INIT(1,K)
               RINVRS = ONE / (RO_INIT(NO_INIT(K),K) - RK)
               DX     = VO_INIT(K) - XIRING

               DO I = 1, NO_INIT(K)
                  RNORM        = (RO_INIT(I,K) - RK) * RINVRS
                  RO_INIT(I,K) = RNORM
                  XO_INIT(I,K) = XO_INIT(I,K) - RNORM * DX - XIRING
               END DO

            ELSE

               RINVRS = ONE / RO_INIT(NO_INIT(K),K)
               DX     = VO_INIT(K) - XICNTR
               DO I = 1, NO_INIT(K)
                  RNORM        = RO_INIT(I,K) * RINVRS
                  RO_INIT(I,K) = RNORM
                  XO_INIT(I,K) = XO_INIT(I,K) - RNORM * DX - XICNTR
               END DO

            END IF

         END DO

      END IF

      DO K = 1, NO_SEC

         IF (NO_INIT(K) /= NUMI) GO TO 900

         IF (ABS (XO_INIT(1,K))          > SMALL .OR. &
             ABS (XO_INIT(NO_INIT(K),K)) > SMALL) THEN
            WRITE (LWRIT, *) &
               ' Input section did not start & end with X(norm) = 0.', &
               ' Outer K: ', K, ' First & last X(norm)s: ', &
               XO_INIT(1,K), XO_INIT(NO_INIT(K),K)
            STOP
         END IF

         IF (ABS (RO_INIT(1,K))                > SMALL .OR. &
             ABS (RO_INIT(NO_INIT(K),K) - ONE) > SMALL) THEN
            WRITE (LWRIT, *) &
               ' Input section did not begin with R = 0. or', &
               ' end with R = 1.  Outer K: ', K, ' First & last R(norm)s: ', &
               RO_INIT(1,K), RO_INIT(NO_INIT(K),K)
            STOP
         END IF

      END DO

      GO TO 999

  900 WRITE (LWRIT, '(/, 2A)') &
         ' SMOOTH_X = T and ''*.spokes'' file output assume constant', &
         ' point counts for defining sections.'

  999 RETURN

      END SUBROUTINE RD_GEOM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE RD_OPT
!
!     Read the primary optimization controls.
!     This is specialized for optimizing heat shield shapes.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Global variables:

      USE AERO_MOD
      USE COEFFS_MOD
      USE CONST_MOD
      USE DESIGN_VARS_MOD
      USE GEOM_MOD
      USE GRID_MOD
      USE OPT_MOD
      USE TRIANG_MOD

      IMPLICIT NONE

!     Local constants:

      REAL, PARAMETER :: ONE = 1., ZERO = 0.

!     Local variables:

      INTEGER   :: I, IOS
      LOGICAL   :: NONZERO
      CHARACTER :: BUFFER*132

!     Execution:

!     Echo the control inputs to the output file:

      WRITE (LWRIT, '(/, A, /)') ' Optimization control inputs:'

      DO
         READ (LREAD, '(A)', IOSTAT=IOS) BUFFER
         IF (IOS < 0) EXIT

         I = LEN_TRIM (BUFFER)
         WRITE (LWRIT, '(1X, A)') BUFFER(1:I)
      END DO

      REWIND (LREAD)

!     Read the control inputs:

      READ (LREAD, *)
      READ (LREAD, *)
      READ (LREAD, *)
      READ (LREAD, *) NDV, NCLIN, NCNLN, NITMAX, GRAD_MODE, ETA, EPSOBJ

      EPSM6TH = EPSOBJ ** (-1./ 6.) ! Applied to diff. interval if GRAD_MODE = 3

      READ (LREAD, *)
      READ (LREAD, *)
      READ (LREAD, *)
      READ (LREAD, *) N_RAD_MAX, N_AZM_MAX, N_INNER, N_OUTER, CUSP, &
                      TWO_PARTS, SMOOTH_X  ! CUSP is scaled by L_ref in RD_GEOM

      IF (.NOT. TWO_PARTS) N_INNER = 1     ! Keep associated allocates minimal

      READ (LREAD, *)
      READ (LREAD, *) RI_SPL, RI_SPQ, RI_SPS, RI_SPC
      READ (LREAD, *)
      READ (LREAD, *) RO_SPL, RO_SPQ, RO_SPS, RO_SPC
      READ (LREAD, *)
      READ (LREAD, *)
      READ (LREAD, *)
      READ (LREAD, *) RHO_CM_TARG, RHO_CM_MAX, RHO_CD, RHO_LD, RHO_TEMP_MAX,   &
                      RHO_CD_INV, RHO_CD_TARG
      READ (LREAD, *)
      READ (LREAD, *) RHO_AL, RHO_CURVI, RHO_CURVJ, RHO_CVIMX, RHO_CVJMX,      &
                      RHO_SMTH_X, RHO_SWET
      FRACTION_CRVIMX = ONE  ! May also be set via XNLCON() if CRVIMX is present
      READ (LREAD, *)
      READ (LREAD, *) RHO_VNORM, RHO_AREA, RHO_VOLUME, RHO_V_EFF, RHO_CD_A,    &
                      RHO_CM_ALF, RHO_TVD_GC
      READ (LREAD, *)
      READ (LREAD, *) CM_TARG, CM_TOL, CD_TARG, AL_TARG, CR_TARG, TM_TARG
      READ (LREAD, *)
      READ (LREAD, *) ALPHA_MODE, Alpha, Beta, M_inf, Gamma

      NONZERO = FALSE
      IF (ABS (RHO_CM_TARG ) > SMALL) NONZERO = TRUE
      IF (ABS (RHO_CM_MAX  ) > SMALL) NONZERO = TRUE
      IF (ABS (RHO_CD      ) > SMALL) NONZERO = TRUE
      IF (ABS (RHO_LD      ) > SMALL) NONZERO = TRUE
      IF (ABS (RHO_TEMP_MAX) > SMALL) NONZERO = TRUE
      IF (ABS (RHO_CD_INV  ) > SMALL) NONZERO = TRUE
      IF (ABS (RHO_CD_TARG ) > SMALL) NONZERO = TRUE
      IF (ABS (RHO_AL      ) > SMALL) NONZERO = TRUE
      IF (ABS (RHO_CURVI   ) > SMALL) NONZERO = TRUE
      IF (ABS (RHO_CURVJ   ) > SMALL) NONZERO = TRUE
      IF (ABS (RHO_CVIMX   ) > SMALL) NONZERO = TRUE
      IF (ABS (RHO_CVJMX   ) > SMALL) NONZERO = TRUE
      IF (ABS (RHO_SMTH_X  ) > SMALL) NONZERO = TRUE
      IF (ABS (RHO_SWET    ) > SMALL) NONZERO = TRUE
      IF (ABS (RHO_VNORM   ) > SMALL) NONZERO = TRUE
      IF (ABS (RHO_AREA    ) > SMALL) NONZERO = TRUE
      IF (ABS (RHO_VOLUME  ) > SMALL) NONZERO = TRUE
      IF (ABS (RHO_V_EFF   ) > SMALL) NONZERO = TRUE
      IF (ABS (RHO_CD_A    ) > SMALL) NONZERO = TRUE
      IF (ABS (RHO_CM_ALF  ) > SMALL) NONZERO = TRUE
      IF (ABS (RHO_TVD_GC  ) > SMALL) NONZERO = TRUE

      IF (.NOT. NONZERO) THEN
         WRITE (LWRIT, *) 'All RHO_* inputs are zero - aborting.'
         STOP
      END IF

      IF (TWO_PARTS) THEN
         IF (N_INNER + N_OUTER - 1 /= N_RAD_MAX) THEN
            WRITE (LWRIT,*) 'N_RAD_MAX /= N_INNER + N_OUTER - 1'
            STOP
         END IF
      ELSE
         IF (N_OUTER /= N_RAD_MAX) THEN
            WRITE (LWRIT,*) 'N_RAD_MAX /= N_OUTER'
            STOP
         END IF
      END IF

!     ------------------------------------
!     Set up various sizes and dimensions:
!     ------------------------------------

      NBNDS = NDV + NCLIN + NCNLN      ! For variable + constraint bounds
      NROWA = MAX (NCLIN, 1)           ! Max. # linear constraints ...
      NROWJ = MAX (NCNLN, 1)           ! ... & nonlinear constraints; dims. > 0
      NROWR = MAX (NDV, 1)             ! Row dim. of Cholesky factor R
                                       ! of the Hessian of the Lagrangian
      LENIW = 3*NDV + NROWA + 2*NROWJ  ! Length of integer workspace
      LENW  = NDV*(2*NDV + NROWA +  &  ! Length of real workspace
              2*NROWJ + 20) + 11*NROWA + &
             21*NROWJ

      CONFUN_DOES_OBJ =  NCNLN > 0
      N_Verts         =  N_AZM_MAX * N_RAD_MAX
      N_Triangles     = (N_AZM_MAX - 1) * (2*(N_RAD_MAX - 1) - 1) ! Accounts for
                                                                  ! singular pt.
      ALLOCATE (XYZ (3,N_Verts),              &
                FACE_AREA (N_Triangles),      &
                SIGNS     (N_Triangles),      &
                ICOMP     (N_Triangles),      &
                CP        (N_Triangles),      &
                IPT_TRI      (3,N_Triangles), &
                XYZ_CENTROID (3,N_Triangles), &
                XYZ_NORM     (3,N_Triangles))

!     ---------------------------------------------------
!     Allocate (most of) the work-space that is dependent
!     on the number of design variables and constraints:
!     ---------------------------------------------------

      ALLOCATE (INLCON(NCNLN), JNLCON(NCNLN), ISTATE(NBNDS), IWVEC(LENIW))

      ALLOCATE &
        (FFWD(NDV),       GRAD(NDV),     CVEC(NROWJ),     WVEC(LENW),          &
         IGROUP(NDV),     AITCH(NDV),    ONEOVERH(NDV),                        &
         BL(NBNDS),       BU(NBNDS),                                           &
         CONLIN(NCLIN),   C(NCNLN),      CPRINT(NCNLN),                        &
         BLLIN(NCLIN),    BULIN(NCLIN),  CLAMDA(NBNDS),                        &
         XNLCON(NCNLN),   SNLCON(NCNLN), TLCON(NCLIN),    ILCON(NCLIN),        &
         AMAT(NCLIN,NDV), RMAT(NDV,NDV), CJAC(NROWJ,NDV), CJACBAK(NROWJ))

      ALLOCATE &
        (V(NDV), VSCALE(NDV), &
         ZBUMP(NDV), ZBMIN(NDV), ZBMAX(NDV), ZBEXP(NDV), &
         RBUMP(NDV), RBMIN(NDV), RBMAX(NDV), RBEXP(NDV))

      ALLOCATE (RBTYP(NDV), ZBTYP(NDV), LCTYPE(NCLIN), NLCTYPE(NCNLN))

      ALLOCATE (RFIN(N_RAD_MAX, N_AZM_MAX), XFIN(N_RAD_MAX, N_AZM_MAX),        &
                YFIN(N_RAD_MAX, N_AZM_MAX), ZFIN(N_RAD_MAX, N_AZM_MAX),        &
                XLEN(N_RAD_MAX, N_AZM_MAX),                                    &
                GAUSS_CURVATURE(N_RAD_MAX, N_AZM_MAX),                         &
                MEAN_CURVATURE (N_RAD_MAX, N_AZM_MAX))

      END SUBROUTINE RD_OPT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE RD_VARS
!
!     Read the design variables.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      USE CONST_MOD
      USE DESIGN_VARS_MOD
      USE GRID_MOD
      USE OPT_MOD

      IMPLICIT NONE

!     Local variables:

      INTEGER :: &
         I, N, NALPHA, NPLAN, NSINE, NCOS, NEXPS, NTRAIL, NLEAD, NWAG
      LOGICAL :: &
         ORDERED
      CHARACTER :: &
         TYPE * 6, TYPE2 * 4

!     Execution:

      NALPHA = 0
      NPLAN  = 0
      NSINE  = 0
      NCOS   = 0
      NEXPS  = 0
      NTRAIL = 0
      NLEAD  = 0
      NWAG   = 0
      FIXED_CG = .TRUE.

      READ (LREAD,*)
      READ (LREAD,*)
      READ (LREAD,*)

      DO I = 1, NDV

         READ (LREAD,*) N, TYPE, TYPE2, &
                        RBEXP(I), RBUMP(I), RBMIN(I), RBMAX(I), &
                        ZBEXP(I), ZBUMP(I), ZBMIN(I), ZBMAX(I), &
                        V(I), VSCALE(I), AITCH(I), BL(I), BU(I)
         IF (BL(I) == -999.) BL(I) = -BIGBND
         IF (BU(I) ==  999.) BU(I) =  BIGBND

         ONEOVERH(I) = 1./ AITCH(I)
         RBTYP(I) = TYPE
         ZBTYP(I) = TYPE2

         IF (TYPE2 == 'POLY') THEN     ! Just check for bad inputs
         ELSE IF (TYPE2 == 'SIN ') THEN
         ELSE IF (TYPE2 == 'SINF') THEN
         ELSE IF (TYPE2 == 'SIN1') THEN
         ELSE IF (TYPE2 == 'SIN2') THEN
         ELSE IF (TYPE2 == 'SINC') THEN
         ELSE IF (TYPE2 == 'COSL') THEN
         ELSE IF (TYPE2 == 'COSR') THEN
         ELSE IF (TYPE2 == 'LCOS') THEN
         ELSE IF (TYPE2 == 'RCOS') THEN
         ELSE IF (TYPE2 == 'LED ') THEN
         ELSE IF (TYPE2 == 'TRL ') THEN
         ELSE IF (TYPE2 == 'EXP ') THEN
         ELSE
            GO TO 800 ! Unknown variable type 2
         END IF

         IF      (TYPE == 'XCEN  ') THEN
            NPLAN = NPLAN + 1
            IGROUP(I) = 1
         ELSE IF (TYPE == 'YCEN  ') THEN
            NPLAN = NPLAN + 1
            IGROUP(I) = 1
         ELSE IF (TYPE == 'ZCEN  ') THEN
            NPLAN = NPLAN + 1
            IGROUP(I) = 1
         ELSE IF (TYPE == 'X_CG  ') THEN  ! Both X_CG & Y_CG should be present
            NPLAN = NPLAN + 1
            IGROUP(I) = 1
            FIXED_CG  = .FALSE.
         ELSE IF (TYPE == 'Y_CG  ') THEN
            NPLAN = NPLAN + 1
            IGROUP(I) = 1
            FIXED_CG  = .FALSE.
         ELSE IF (TYPE == 'XRNG  ') THEN
            NPLAN = NPLAN + 1
            IGROUP(I) = 1
         ELSE IF (TYPE == 'ROUT  ') THEN
            NPLAN = NPLAN + 1
            IGROUP(I) = 1
         ELSE IF (TYPE == 'XOUT  ') THEN
            NPLAN = NPLAN + 1
            IGROUP(I) = 1
         ELSE IF (TYPE == 'TOUT  ') THEN
            NPLAN = NPLAN + 1
            IGROUP(I) = 1
         ELSE IF (TYPE == 'ISIN  ') THEN
            NSINE  = NSINE  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'ISINF ') THEN
            NSINE  = NSINE  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'ISIN1 ') THEN
            NSINE  = NSINE  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'ISIN2 ') THEN
            NSINE  = NSINE  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'ICOSL ') THEN
            NCOS  = NCOS  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'ICOSR ') THEN
            NCOS  = NCOS  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'ILCOS ') THEN
            NCOS  = NCOS  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'IRCOS ') THEN
            NCOS  = NCOS  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'IEXP  ') THEN
            NEXPS  = NEXPS  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'ITRL  ') THEN
            NTRAIL = NTRAIL + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'ILED  ') THEN
            NLEAD  = NLEAD  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'IWAG  ') THEN
            NWAG = NWAG  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'OSIN  ') THEN
            NSINE  = NSINE  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'OSINF ') THEN
            NSINE  = NSINE  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'OSIN1 ') THEN
            NSINE  = NSINE  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'OSIN2 ') THEN
            NSINE  = NSINE  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'OCOSL ') THEN
            NCOS  = NCOS  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'OCOSR ') THEN
            NCOS  = NCOS  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'OLCOS ') THEN
            NCOS  = NCOS  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'ORCOS ') THEN
            NCOS  = NCOS  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'OEXP  ') THEN
            NEXPS  = NEXPS  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'OTRL  ') THEN
            NTRAIL = NTRAIL + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'OLED  ') THEN
            NLEAD  = NLEAD  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'OWAG  ') THEN
            NWAG = NWAG  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'ALPHA ') THEN
            NALPHA = NALPHA + 1
            IALPHA = I
            IGROUP(I) = 3
         ELSE
            GO TO 800  ! Unknown variable type
         END IF

      END DO

      NPLANF = NPLAN
      NBUMPS = NSINE + NCOS + NEXPS + NTRAIL + NLEAD + NWAG

      IF (NPLANF + NBUMPS + NALPHA /= NDV) THEN
         WRITE (LWRIT,*) 'Total # design variables did not match NDV.'
         GO TO 900
      END IF

      ORDERED = .TRUE.
      DO I = 2, NDV
         IF (IGROUP(I) < IGROUP(I-1)) ORDERED = .FALSE.
      END DO

      DEALLOCATE (IGROUP)

      IF (.NOT. ORDERED) THEN
         WRITE (LWRIT, '(/, 1X, A, A, /, (4X, A))')                   &
            'RDVARS:  The optimization variables must be grouped ',   &
            'in the following order:', 'Planform (XCEN, YCEN, etc.)', &
            'Spoke perturbing functions (ISIN1, OSIN1, etc.)',        &
            'Angle of attack (only if ALPHA_MODE = 2)'
         GO TO 900
      END IF

      GO TO 999

  800 WRITE (LWRIT, '(/, A, I5, 2X, A, 2X, A)') &
         ' RDVARS: Unknown design variable type(s): ', I, TYPE, TYPE2

  900 STOP

  999 RETURN

      END SUBROUTINE RD_VARS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE RIPPER (Z, ZCEN, ZMIN, ZMAX, PWR, HEIGHT)
!
!     The name is all that's left of James's bizarre original utility for
!     linear/polynomial spanwise lofting of the wing section perturbations.
!     The index-based scheme has been generalized to lofting between arbitrary
!     spanwise locations.  HEIGHT is returned between 0. and 1.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Arguments:

      REAL, INTENT (IN)  :: Z, ZCEN, ZMIN, ZMAX, PWR
      REAL, INTENT (OUT) :: HEIGHT

!     Constants:

      REAL, PARAMETER :: ONE = 1., ZERO = 0.

!     Execution:

      IF (Z == ZCEN) THEN
         HEIGHT = ONE
      ELSE IF (PWR == ZERO) THEN
         IF (Z < ZMIN .OR. Z > ZMAX) THEN
            HEIGHT = ZERO
         ELSE
            HEIGHT = ONE
         END IF
      ELSE IF (Z < ZCEN) THEN
         IF (Z <= ZMIN) THEN ! Handles ZMIN = ZCEN case
            HEIGHT = ZERO
         ELSE
            HEIGHT = ((Z - ZMIN) / (ZCEN - ZMIN)) ** PWR
         END IF
      ELSE
         IF (Z >= ZMAX) THEN
            HEIGHT = ZERO
         ELSE
            HEIGHT = ((ZMAX - Z) / (ZMAX - ZCEN)) ** PWR
         END IF
      END IF

      END SUBROUTINE RIPPER

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE SETLCON
!
!     SETLCON sets up the linear constraint matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Global variables:

      USE DESIGN_VARS_MOD
      USE CONST_MOD
      USE OPT_MOD
      USE GEOM_MOD

      IMPLICIT NONE

!     Local constants:

      REAL, PARAMETER :: ZERO = 0.

!     Local variables:

      REAL    :: BOUNDL, BOUNDU, DELTAV, T, ZW
      INTEGER :: KW, M, N

!     Execution:

!     Zero the linear constraint matrix - only some of the variables
!     are applicable.

      AMAT = ZERO ! (1:NCLIN,1:NDV)

!     Determine the nonzero matrix elements:

      DO M = 1, NCLIN

         KW     = ILCON(M)
         ZW     = TLCON(M)
         BOUNDL = BLLIN(M)
         BOUNDU = BULIN(M)

         BU(M+NDV) = BOUNDU
         BL(M+NDV) = BOUNDL

         IF (LCTYPE(M)(1:6) == 'EDGE_K') THEN

            DO N = 1, NPLANF

                DELTAV = ZERO

                IF (RBTYP(N)(1:4) == 'TOUT' .OR. &
                    RBTYP(N)(1:4) == 'ROUT') THEN

                   CALL LOFT (ZBTYP(N), ZO_INIT(KW), ZBUMP(N), ZBMIN(N), &
                              ZBMAX(N), ZBEXP(N), CIRCLE, DELTAV, LWRIT)

                   AMAT(M,N) = DELTAV
                END IF

            END DO

         ELSE IF (LCTYPE(M)(1:6) == 'EDGE_Z') THEN

            DO N = 1, NPLANF

                DELTAV = ZERO

                IF (RBTYP(N)(1:4) == 'TOUT' .OR. &
                    RBTYP(N)(1:4) == 'ROUT') THEN

                   CALL LOFT (ZBTYP(N), ZW, ZBUMP(N), ZBMIN(N), &
                              ZBMAX(N), ZBEXP(N), CIRCLE, DELTAV, LWRIT)

                   AMAT(M,N) = DELTAV
                END IF

            END DO

         END IF

      END DO ! Next linear constraint

      END SUBROUTINE SETLCON

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE SETUP_NPOPT (MODE)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      USE CONST_MOD
      USE OPT_MOD

      IMPLICIT NONE

!     Arguments:

      INTEGER, INTENT (OUT) :: MODE ! 0 means no error found

!     Namelist to override NPOPT's options file defaults:

      NAMELIST /NPOPTIONS/ &
         LEVELVER, MAJORPL, MINORPL, MINORIL, NPRFREQ, &
         PENPARAM, STEPLIM, TOLLIN, TOLNLIN, TOLOPT

!     Execution:

      STEPLIM  = 2.       ! 2 means allow a 200% change in the variables
      PENPARAM = 0.
      LEVELVER = -1
      TOLOPT   = 1.E-6
      TOLLIN   = TOLOPT
      TOLNLIN  = TOLOPT
      MAJORPL  = 10
      MINORPL  = 0
      MINORIL  = 1000
      NPRFREQ  = 100
      NFTOTL   = 0        ! Incremented in CONFUN and/or OBJFUN

!     Look for an optional namelist at the end of standard input:

      READ (LREAD, NPOPTIONS, ERR=700, END=700)

  700 CLOSE (LREAD)

!     -------------------------------------------------
!     Run-time set-up of the optional inputs for NPOPT:
!     -------------------------------------------------

      OPEN (UNIT=LNPOPT, FILE='heat_shield.nps', STATUS='UNKNOWN')

      WRITE (LNPOPT, '(A)') 'Begin'

      WRITE (LNPOPT, '(A)') &
         'Nonderivative line search', &
         'Hessian       Full memory'
      WRITE (LNPOPT, '(A, I6)') &
         'Major Iteration Limit          ', NITMAX,   &
         'Minor Iteration Limit          ', MINORIL,  &
         'Major Print Level              ', MAJORPL,  &
         'Minor Print Level              ', MINORPL,  &
         'Print frequency                ', NPRFREQ,  &
         'Verify level                   ', LEVELVER, &
         'Derivative Level               ', 3         ! Why was it ever 0?

      WRITE (LNPOPT, '(A, E12.3)') &
!!!!!    'Crash Tolerance                ', 0.0,      &
         'Function Precision             ', EPSOBJ,   &
         'Infinite Bound                 ', BIGBND,   &
         'Step Limit                     ', STEPLIM,  &
         'Linear Feasibility Tolerance   ', TOLLIN,   &
         'Nonlinear Feasibility Tolerance', TOLNLIN,  &
         'Optimality Tolerance           ', TOLOPT,   &
         'Linesearch Tolerance           ', ETA,      &
         'Penalty parameter              ', PENPARAM, &
         'End'

      REWIND (LNPOPT)

      WRITE (LWRIT, '(//, A)') ! *** ' Options file for NPOPT:'

      CALL NPPRMS (LNPOPT, INFORM)

      CLOSE (LNPOPT)

      IF (INFORM /= 0) THEN
         WRITE (LWRIT, 1010) &
            'An error has been detected in the NPOPT options file.'
         GO TO 900
      END IF

      MODE = 0

      RETURN

  900 CONTINUE

      MODE = 1

      RETURN

 1010 FORMAT (/, ' SETUP_NPOPT: ', A)

      END SUBROUTINE SETUP_NPOPT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE SFEVAL (NAME, EX, XC, XA, XB, YMULT, I1, I2, X, Y, LUN)
!
!        SFEVAL evaluates a specified shape function at X (I1:I2) and
!     adds the result to the corresponding Y.  X need not be normalized
!     because the option to apply the shape function to just part of the
!     data range, [XA, XB], forces normalization locally.  XC, however,
!     refers to the normalized sub-range, so it must be normalized.
!
!        These are the Hicks-Henne shape functions used for aerodynamic
!     design.  SFEVAL is an adaptation of BEVAL from program PROFILE in
!     the form most convenient for design code SYN87, where a modular
!     form is needed for precise set-up of linear constraints as well
!     as for the geometry perturbations.
!
!        Distinguishing wing and body applications should be done at
!     the calling level to avoid too many similar cases here.
!
!     06/21/96    DAS    Adaptation of BEVAL and PERTURB routines.
!     08/22/96     "     Added SINF. Apart from the SIN* functions,
!                        testing just 3 characters avoids problems.
!     08/24/96     "     LED & LEA now both work for LEADing;
!                        TRL & TRA  "    "    "   "  TRAILing.
!     12/06/96  DAS/JJR  Added COSL & COSR for body keel & crown,
!                        and inverted forms LCOS & RCOS.
!     06/03/97    DAS    Converting flap & slap angles from degrees to
!                        radians was using 180./PI instead of PI/180.
!     08/09/97     "     FLAP/SLAT now expect hinge X in same units as
!                        X(*) and Y(*) (not fraction of XB - XA).
!     08/10/99     "     Added SIN3 (fore/aft symmetry via SIN1 on [0,1])
!                        and   SIN4 (  "   "   "   "   "   "   " each half).
!
!     Author:  David Saunders, Sterling Software/NASA Ames, Moffett Field, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IMPLICIT NONE

!     Arguments:

      CHARACTER  NAME * (*) ! I   Shape function name (upper case;
                            !     only characters 1:4 are checked here,
                            !     or 1:3 where that clearly suffices
      REAL       EX         ! I   Exponent, usually affecting "width";
                            !     for 'WAGN'[er] fn. N pass REAL (N)
      REAL       XC         ! I   "Center" X in [0, 1] where fn. peaks
      REAL       XA, XB     ! I   X sub-range in the same units as X (*)
      REAL       YMULT      ! I   Shape function multiplier
      INTEGER    I1, I2     ! I   Index range of X & Y to be treated
      REAL       X (*)      ! I   Abscissas, not necessarily in [0, 1]
      REAL       Y (*)      ! I/O Ordinates being perturbed
      INTEGER    LUN        ! I   Logical unit for error messages

!     Local constants:

      REAL, PARAMETER :: HALF = 0.5, ONE = 1., TWO = 2., ZERO = 0.

!     Local variables:

      INTEGER   :: I
      REAL      :: AEXP, BEXP, CENTER, CENTER2, DEGRAD, EPSMCH, PT5LOG,    &
                   ONEMC, PI, PIBY2, POWER, POWERL, POWERR, RN, RNM1, RPI, &
                   RRANGE, SCALE, TANGNT, THETA, WIDTH, X0, XI, XSCALE
      LOGICAL   :: FIRST
      CHARACTER :: KEY4 * 4, KEY3 * 3

!     Storage:

      DATA       FIRST /.TRUE./
      SAVE       FIRST, DEGRAD, PI, PIBY2, RPI, PT5LOG

!     Statement functions:

      REAL :: &
         EBUMP, SBUMP, PWR, WDTH, XNRM

      EBUMP (WDTH, PWR, XNRM) = &
         XNRM ** PWR * (ONE - XNRM) * EXP (-WDTH * XNRM)

      SBUMP (WDTH, PWR, XNRM) = (SIN (PI * (XNRM ** PWR))) ** WDTH

!     Execution:

      IF (FIRST) THEN

         FIRST = .FALSE.

!        Protect sine bump exponentiation from slightly negative arguments
!        by calculating a value for PI such that SIN (PI) > 0 (barely).
!        This allows the limiting of X to [0, 1] that is necessary for
!        the sub-range capability to cover the sine protection too.

         EPSMCH = EPSILON (EPSMCH)

         PI = ASIN (ONE) * TWO

         DO I = 1, 999 ! 6 iterations observed for DEC AXP  (32 bits)
                       ! 5     "         "      "  CRAY C90 (64 bits)
!!!!        write (6,*) 'I   PI  SIN (PI): ', i, pi, sin (pi)
            IF (SIN (PI) > ZERO) GO TO 10

            PI = PI - EPSMCH * REAL (I)
         END DO

         WRITE (LUN, '(A)') ' *** SFEVAL: Safeguard iteration failed.'
         GO TO 900

   10    CONTINUE

         PIBY2  = PI * HALF
         DEGRAD = PI / 180.
         RPI    = ONE / PI
         PT5LOG = LOG (HALF)
      END IF

!     Minimize argument references in loops:

      WIDTH  = EX
      CENTER = XC ! Normalized
      X0     = XA
      RRANGE = ONE / (XB - X0)
      SCALE  = YMULT
      KEY4   = NAME (1:4) ! Avoid comparison of different-length strings
      KEY3   = KEY4 (1:3)

!     Treat the cases in ~ the most likely order they'll be encountered:

      IF (KEY4 == 'SIN ') THEN ! Modified "SINE" (unsymmetric)

         POWER = PT5LOG / LOG (CENTER)
         DO I = I1, I2
            XI = MAX (ZERO, MIN ((X (I) - X0) * RRANGE, ONE))
            Y (I) = SCALE * SBUMP (WIDTH, POWER, XI) + Y (I)
         END DO

      ELSE IF (KEY4 == 'SIN1') THEN ! Modified "SINE" (symmetric)

         IF (CENTER <= HALF) THEN
            POWER = PT5LOG / LOG (CENTER)
            DO I = I1, I2
               XI = MAX (ZERO, MIN ((X (I) - X0) * RRANGE, ONE))
               Y (I) = SCALE * SBUMP (WIDTH, POWER, XI) + Y (I)
            END DO
         ELSE
            POWER = PT5LOG / LOG (ONE - CENTER)
            DO I = I1, I2
               XI = MAX (ZERO, MIN ((X (I) - X0) * RRANGE, ONE))
               Y (I) = SCALE * SBUMP (WIDTH, POWER, ONE - XI) + Y (I)
            END DO
         END IF

      ELSE IF (KEY4 == 'SIN2') THEN ! Modified "SINE" (symmetric)

         IF (CENTER <= HALF) THEN
            POWER = PT5LOG / LOG (ONE - CENTER)
            DO I = I1, I2
               XI = MAX (ZERO, MIN ((X (I) - X0) * RRANGE, ONE))
               Y (I) = SCALE * SBUMP (WIDTH, POWER, ONE - XI) + Y (I)
            END DO
         ELSE
            POWER = PT5LOG / LOG (CENTER)
            DO I = I1, I2
               XI = MAX (ZERO, MIN ((X (I) - X0) * RRANGE, ONE))
               Y (I) = SCALE * SBUMP (WIDTH, POWER, XI) + Y (I)
            END DO
         END IF

      ELSE IF (KEY4 == 'SIN3') THEN ! Fore/aft symmetry via SIN1 on [0, 1]

!        Assume [XA, XB] = [0, 1] and CENTER in (0, 0.5]:

         POWER = PT5LOG / LOG (CENTER)

         DO I = I1, I2
            XI = MAX (ZERO, MIN (X (I), ONE))
            IF (XI > HALF) XI = ONE - XI
            Y (I) = SCALE * SBUMP (WIDTH, POWER, XI) + Y (I)
         END DO

      ELSE IF (KEY4 == 'SIN4') THEN ! Fore/aft symmetry via SIN1 on each half

!        Assume [XA, XB] = [0, 0.5] and CENTER in (0, 1).  Left half:

         POWERL  = PT5LOG / LOG (CENTER)
         CENTER2 = ONE - CENTER
         POWERR  = PT5LOG / LOG (CENTER2)

         IF (CENTER <= HALF) THEN
            DO I = I1, I2
               XI = MAX (ZERO, MIN (X (I) * RRANGE, ONE))
               Y (I) = SCALE * SBUMP (WIDTH, POWERL, XI) + Y (I)
            END DO
         ELSE
            DO I = I1, I2
               XI = MAX (ZERO, MIN (X (I) * RRANGE, ONE))
               Y (I) = SCALE * SBUMP (WIDTH, POWERR, ONE - XI) + Y (I)
            END DO
         END IF

!        Now the [0.5, 0.1] half with center <-- 1 - center

         IF (CENTER2 <= HALF) THEN
            DO I = I1, I2
               XI = MAX (ZERO, MIN ((X (I) - HALF) * RRANGE, ONE))
               Y (I) = SCALE * SBUMP (WIDTH, POWERR, XI) + Y (I)
            END DO
         ELSE
            DO I = I1, I2
               XI = MAX (ZERO, MIN ((X (I) - HALF) * RRANGE, ONE))
               Y (I) = SCALE * SBUMP (WIDTH, POWERL, ONE - XI) + Y (I)
            END DO
         END IF

      ELSE IF (KEY4 == 'SINF') THEN ! Flipped "SIN " (unsymmetric)

         POWER = PT5LOG / LOG (ONE - CENTER)
         DO I = I1, I2
            XI = MAX (ZERO, MIN ((X (I) - X0) * RRANGE, ONE))
            Y (I) = SCALE * SBUMP (WIDTH, POWER, ONE - XI) + Y (I)
         END DO

      ELSE IF (KEY3 == 'EXP') THEN ! "EXPONENTIAL" (peak height = 1.)

         ONEMC = ONE - CENTER
         AEXP  = CENTER * (ONE + ONEMC * WIDTH) / ONEMC
         BEXP  = SCALE / EBUMP (WIDTH, AEXP, CENTER)

         DO I = I1, I2
            XI = MAX (ZERO, MIN ((X (I) - X0) * RRANGE, ONE))
            Y (I) = BEXP * EBUMP (WIDTH, AEXP, XI)  +  Y (I)
         END DO

      ELSE IF (KEY3(1:2) == 'LE') THEN ! "LEADING"-edge (LED or LEA)

         DO I = I1, I2
            XI = MAX (ZERO, MIN ((X (I) - X0) * RRANGE, ONE))
            Y (I) = ((ONE - XI) ** WIDTH) * SCALE + Y (I)
         END DO

      ELSE IF (KEY3(1:2) == 'TR') THEN ! "TRAILING"-edge (TRL or TRA)

         DO I = I1, I2
            XI = MAX (ZERO, MIN ((X (I) - X0) * RRANGE, ONE))
            Y (I) = (XI ** WIDTH) * SCALE + Y (I)
         END DO

      ELSE IF (KEY4 == 'COSL') THEN ! Half a cosine (or sine), peak at left

!        SIN is used instead of COS because SIN (PI) is protected above.

         XSCALE = PIBY2 * RRANGE
         X0 = TWO * X0 - XB

         DO I = I1, I2
            XI = MAX (PIBY2, MIN ((X (I) - X0) * XSCALE, PI))
            Y (I) = SCALE * (SIN (XI)) ** WIDTH + Y (I)
         END DO

      ELSE IF (KEY4 == 'COSR') THEN ! Half a cosine (or sine), peak at right

!        SIN is used instead of COS for consistency with COSL.

         XSCALE = PIBY2 * RRANGE

         DO I = I1, I2
            XI = MAX (ZERO, MIN ((X (I) - X0) * XSCALE, PIBY2))
            Y (I) = SCALE * (SIN (XI)) ** WIDTH + Y (I)
         END DO

      ELSE IF (KEY4 == 'LCOS') THEN ! Inverted half (co)sine, peak at left

         XSCALE = PIBY2 * RRANGE

         DO I = I1, I2
            XI = MAX (ZERO, MIN ((X (I) - X0) * XSCALE, PIBY2))
            Y (I) = SCALE * (ONE - (SIN (XI)) ** WIDTH) + Y (I)
         END DO

      ELSE IF (KEY4 == 'RCOS') THEN ! Inverted half (co)sine, peak at right

         XSCALE = PIBY2 * RRANGE
         X0 = TWO * X0 - XB

         DO I = I1, I2
            XI = MAX (PIBY2, MIN ((X (I) - X0) * XSCALE, PI))
            Y (I) = SCALE * (ONE - (SIN (XI)) ** WIDTH) + Y (I)
         END DO

      ELSE IF (KEY3 == 'WAG') THEN ! Wagner function # N

!        Reference:  Ramamoorthy, P. and Padmavathi, K.  "Airfoil Design
!        by Optimization" in J. Aircraft, Vol. 14 (Feb. 1977), 219-221.

         IF (WIDTH == ONE) THEN ! N = 1
            DO I = I1, I2
               XI = MAX (ZERO, MIN ((X (I) - X0) * RRANGE, ONE))
               THETA = TWO * ASIN (SQRT (XI))
               Y (I) = ((THETA + SIN (THETA)) * RPI - &
                        (SIN (HALF * THETA)) ** 2) * SCALE + Y (I)
            END DO
         ELSE
            RN   = ONE / WIDTH ! 1 / N
            RNM1 = WIDTH - ONE ! N - 1
            DO I = I1, I2
               XI = MAX (ZERO, MIN ((X (I) - X0) * RRANGE, ONE))
               THETA = TWO * ASIN (SQRT (XI))
               Y (I) = ((SIN (WIDTH * THETA) * RN  + &
                        SIN (RNM1 * THETA)) * RPI) * SCALE + Y (I)
            END DO
         END IF

      ELSE IF (KEY3 == 'FLA') THEN

!        "FLAP"-like function (shearing transformation only):
!        XC is the hinge point in the same frame as X & Y;
!        YMULT is the flap angle in DEGREEs (positive is "down").

         TANGNT = TAN (DEGRAD * SCALE)
         DO I = I1, I2
            XI = MAX (ZERO, X (I) - CENTER)
            Y (I) = Y (I) - XI * TANGNT
         END DO

      ELSE IF (KEY3 == 'SLA') THEN

!        "SLAT"-like function (shearing transformation only):
!        XC is the hinge point in the same frame as X & Y;
!        YMULT is the slat angle in DEGREEs (positive is "down").

         TANGNT = TAN (DEGRAD * SCALE)
         DO I = I1, I2
            XI = MAX (ZERO, CENTER - X (I))
            Y (I) = Y (I) - XI * TANGNT
         END DO

      ELSE IF (KEY3 == 'DRO') THEN ! "DROOP":

         DO I = I1, I2
            XI = MAX (ZERO, MIN ((X (I) - X0) * RRANGE, ONE))
            Y (I) = ((ONE - XI) * EXP (-WIDTH * XI)) * SCALE + Y (I)
         END DO

      ELSE IF (KEY3 == 'RAM') THEN ! "RAMP":  Y = YMULT * X

!        Commonly used in conjunction with Wagner functions.

         DO I = I1, I2
            XI = MAX (ZERO, MIN ((X (I) - X0) * RRANGE, ONE))
            Y (I) = SCALE * XI + Y (I)
         END DO

      ELSE IF (KEY3 == 'SCA') THEN

!        Simple scaling - X(I) are probably ordinates.  Note that the
!        effect is arranged to be additive, not multiplicative, so that
!        scaling can be treated just as for the other functions when it
!        is included in a group for perturbing purposes.

         DO I = I1, I2
            Y (I) = X (I) * (SCALE - ONE) + Y (I)
         END DO

      ELSE
         WRITE (LUN, '(A, A)') ' *** SFEVAL: Unknown shape fn.: ', NAME
         GO TO 900
      END IF

      GO TO 999

  900 STOP

  999 RETURN

      END SUBROUTINE SFEVAL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE SHIELD_TRI
!
!     Triangulate the heat shield, using the routine from program HYPER_AERO.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Global variables:

      USE AERO_MOD
      USE CONST_MOD
      USE GRID_MOD
      USE TRIANG_MOD

      IMPLICIT NONE

!     Local variables: arrays of length 1 patch avoid compiler warnings.

      INTEGER :: &
         I, K, NFACES, NPTS, ICOMPONENT(1), IP(1), NIJ(2)
      REAL :: &
         SWITCH_NORM(1)

!     Execution:

      NPTS = 0
      DO K = 1, N_AZM_MAX
         DO I = 1, N_RAD_MAX
            NPTS = NPTS + 1
            XYZ(1,NPTS) = XFIN(I,K)
            XYZ(2,NPTS) = YFIN(I,K)
            XYZ(3,NPTS) = ZFIN(I,K)
         END DO
      END DO

      NIJ(1) = N_RAD_MAX
      NIJ(2) = N_AZM_MAX
      IP(1)  = 1
      ICOMPONENT(1)  =  1
      SWITCH_NORM(1) = -1. ! k should point into the wind with (i,j,k) RH rule

      CALL TRIANGULATE_PATCHES (1, NIJ, XYZ, IP, ICOMPONENT, &
                                SWITCH_NORM, REF_LEN,        &
                                N_Verts, NFACES, IPT_TRI, ICOMP, SIGNS)

      IF (NFACES /= N_Triangles) THEN
         WRITE (LWRIT, '(/, A, /, (A, I10))') &
            ' Number of triangles error.',    &
            ' Allocated:', N_Triangles,       &
            ' Assigned :', NFACES
         CALL SUMMARY (999., 0)
         STOP
      END IF

!     Calculate surface geometric quantities independent of flight conditions:

      CALL TRIANGLE_PROPERTIES (N_Verts, NFACES, IPT_TRI, XYZ, SIGNS,  &
                                FACE_AREA, XYZ_NORM, XYZ_CENTROID, SWET)

      END SUBROUTINE SHIELD_TRI

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE SMOOTH_X_VS_PHI (NI, NSEC, AZ, X, SMOOTH_X, X_RESIDUAL)
!
!     Smooth X vs. Azimuth via best-fit quartics.  Only the highest coefficient
!     is unknown if the quartic is forced to have a maximum at the AZ = pi pt.
!
!     02/04/02  David Saunders  Initial implementation.
!     02/07/02    "      "      Option to suppress the update but calculate
!                               the deviations for forcing towards zero.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Global variables:

      USE CONST_MOD

      IMPLICIT NONE

!     Arguments:

      INTEGER, INTENT (IN) :: &
         NI,                  & ! First dimension of X(*,NSEC)
         NSEC                   ! Number of sections to smooth across

      REAL, INTENT (IN)    :: &
         AZ(NSEC)               ! Defining section azimuthal angles (radians)

      REAL, INTENT (INOUT) :: &
         X(NI,NSEC)             ! Defining section X coords., smoothed in-place

      LOGICAL, INTENT (IN) :: &
         SMOOTH_X               ! T means update the Xs & return X_RESIDUAL;
                                ! F means just return X_RESIDUAL

      REAL, INTENT (INOUT) :: &
         X_RESIDUAL             ! Sum of squares of deviations from the best fit
                                ! quartics over all fits

!     Local constants:

      REAL, PARAMETER :: &
         ZERO = 0.

!     Local variables:

      INTEGER :: &
         I, J

      REAL :: &
         A, VJ, VTR, VTV, WSQ, X_AT_PI, V_LHS(NSEC)

!     Execution:

!     If vj = wj**2(wj**2 - 2pi**2), wj = phij - pi, rj = xj - x(pi), then
!     the coefficients for each quartic in (phi - pi) are as follows:
!
!        a = (v'r) / (v'v)   ! Normal equations formula for one coefficient
!        b = 0
!        c = -2 pi**2 a
!        d = 0
!        e = x(pi)

      DO I = 2, NI - 1

         VTR = ZERO
         VTV = ZERO
         X_AT_PI = X(I,NSEC) ! Half shield case

         DO J = 1, NSEC
            WSQ = (AZ(J) - PI) ** 2
            VJ  = WSQ * (WSQ - TWOPISQ)
            VTR = VTR + VJ * (X(I,J) - X_AT_PI)
            VTV = VTV + VJ * VJ
            V_LHS(J) = VJ
         END DO

         A = VTR / VTV

!        Evaluate the quartic at each azimuth:

         DO J = 1, NSEC - 1 ! Unchanged at phi = pi
            V_LHS(J) = A * V_LHS(J) + X_AT_PI ! Overwrite with smoothed X
            X_RESIDUAL = X_RESIDUAL + (V_LHS(J) - X(I,J)) ** 2
         END DO

         IF (SMOOTH_X) THEN ! Do the update of X with the smoothed Xs
            DO J = 1, NSEC - 1 ! Unchanged at phi = pi
               X(I,J) = V_LHS(J)
            END DO
         END IF

      END DO

      END SUBROUTINE SMOOTH_X_VS_PHI

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE SOLVE (NDVPAR, VPAR, OBJPAR, NCALL, FAIL)
!
!     SOLVE is the objective function routine expected by the unconstrained
!     optimizer QNMDIF2 (no longer included).  It is also called by the OBJFUN
!     and CONFUN routines expected by the constrained optimizer NPOPT.
!
!     NCALL usage:
!
!     -13 = first call to SOLVE from CONFUN
!           (precedes first call to OBJFUN)
!     -12 = constraint evaluations following a line search
!     -10 = constraint evaluations within a line search
!      -3 = first call to SOLVE from OBJFUN
!      -2 = objective evaluation after a line search
!      -1 = local search evaluation (QNMDIF2 mode only, inactive)
!       0 = objective evaluation within a line search
!      >0 = gradient element NCALL evaluation
!           (constraints, objective function, or both)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Global variables:

      USE AERO_MOD
      USE COEFFS_MOD
      USE CONST_MOD
      USE GRID_MOD
      USE DESIGN_VARS_MOD
      USE OPT_MOD
      USE TRIANG_MOD

      IMPLICIT NONE

!     Arguments:

      INTEGER, INTENT (IN)    :: NDVPAR  ! Number of design variables
      INTEGER, INTENT (INOUT) :: NCALL   ! See usage above
      REAL,    INTENT (OUT)   :: OBJPAR  ! Objective function value
      REAL,    INTENT (INOUT) :: VPAR(*) ! * because NDVPAR may be zero.
      LOGICAL, INTENT (OUT)   :: FAIL    ! T means the function value is bad

!     Local constants:

      INTEGER,   PARAMETER :: ICALL = 1
      REAL,      PARAMETER :: ZERO  = 0.
      LOGICAL,   PARAMETER :: PRINT = .TRUE.
      CHARACTER, PARAMETER :: SUBNAME*5 = 'SOLVE'

!     Local variables:

      INTEGER    :: I, MCALL, N
      REAL       :: Alpha_TMP, DA, CM_TMP
      REAL(KIND=4), SAVE :: CPU1, CPU2, CPUGRAD1, DAA, TIME, TIME1, TIME2
      LOGICAL    :: FIRST, OUTPUTS

!     Execution:

!     -------------------------------------------------------------------
!     Clarifications:
!
!     (1) CONFUN_DOES_OBJ = T means OBJFUN does not call SOLVE at all.
!     (2) CONFUN_DOES_OBJ = T means CONFUN calls must do the outputs that
!         would otherwise be done by calls from OBJFUN.
!     (3) NCALL = -13, -12, -10  correspond to CONFUN calls equivalent to
!         NCALL =  -3,  -2,   0  calls from OBJFUN.
!     -------------------------------------------------------------------

      FAIL  = .FALSE.
      MCALL = NCALL

      IF (MCALL < 0) MCALL = MOD (MCALL, 10)

      FIRST = MCALL == -3
      IF (FIRST) DAA = SQRT (EPSILON (DAA))

      OUTPUTS = MCALL < 0

      CALL SECOND (TIME1)

      IF (MCALL == 1) THEN ! Start of a gradient calculation
         CPUGRAD1 = TIME1
      ELSE IF (MCALL == -2) THEN
         NDITER = NDITER + 1
      END IF

      IF (OUTPUTS) THEN
         WRITE (LWRIT, '(/, 1X, 2A, I4, A, I5)') SUBNAME, &
            ': Design iteration:', NDITER, '      NCALL:', NCALL
         WRITE (LWRIT, '(/, A, /)') ' DESIGN VARIABLES:'
         WRITE (LWRIT, '(5(I5, E21.13))') (I, VPAR(I), I = 1, NDVPAR)
      END IF

!     ---------------------------
!     Apply any design variables:
!     ---------------------------

      CALL PERTURB (VPAR, MCALL)  ! Volume and sectional area are computed here

!     -------------------------
!     Triangulate the geometry:
!     -------------------------

      CALL SHIELD_TRI             ! Wetted area is calculated here

!     Effective volume coefficient:
!     SIXROOTPI has been divided by sqrt (2) if FULL = 0 (half config.).
!     The sectional area at the base may be omitted in the literature, but
!     values greater than 1 (highest possible, for a closed sphere) would be
!     typical without it.

      DA = SWET + SECTIONAL_AREA
      V_EFF = SIXROOTPI * VOLUME / (DA * SQRT (DA))

!     ---------------------------
!     Calculate the aerodynamics:
!     ---------------------------

      IF (MCALL <= 0) CALL SECOND (CPU1)

      IF (ALPHA_MODE /= 1) THEN ! 0 = Alpha fixed; 2 = optimization variable

         IF (ALPHA_MODE == 2) Alpha = V(IALPHA) * VSCALE(IALPHA)

         IF (NEED_DCMDA) THEN  ! Do forward step first (need normal CL, CD ..)

            Alpha_TMP = Alpha
            Alpha     = Alpha + DAA

            CALL NEWTONIAN_CLCDCM (N_Triangles, FACE_AREA,         &
                                   XYZ_NORM, XYZ_CENTROID, XYZ_CG, &
                                   REF_LEN, S_REF, Gamma,          &
                                   M_inf, Alpha, Beta,             &
                                   CL, CD, CS, MN, MM, ML, CP)

            CM_TMP = MN
            Alpha  = Alpha_TMP

         END IF

         CALL NEWTONIAN_CLCDCM (N_Triangles, FACE_AREA,            &
                                XYZ_NORM, XYZ_CENTROID, XYZ_CG,    &
                                REF_LEN, S_REF, Gamma,             &
                                M_inf, Alpha, Beta,                &
                                CL, CD, CS, MN, MM, ML, CP)

         IF (NEED_DCMDA) DCMDA = (CM_TMP - MN) / DAA

      ELSE                      ! Forced trim mode

         DO N = 1, 21

            CALL NEWTONIAN_CLCDCM (N_Triangles, FACE_AREA,         &
                                   XYZ_NORM, XYZ_CENTROID, XYZ_CG, &
                                   REF_LEN, S_REF, Gamma,          &
                                   M_inf, Alpha, Beta,             &
                                   CL, CD, CS, MN, MM, ML, CP)

            IF (ABS (MN - CM_TARG) < CM_TOL) EXIT

            IF (N == 21) THEN
               WRITE (LWRIT, *) ' WARNING: No CM convergence.'
!              STOP
            END IF

            CM_TMP = MN
            Alpha  = Alpha + DAA

            CALL NEWTONIAN_CLCDCM (N_Triangles, FACE_AREA,         &
                                   XYZ_NORM, XYZ_CENTROID, XYZ_CG, &
                                   REF_LEN, S_REF, Gamma,          &
                                   M_inf, Alpha, Beta,             &
                                   CL, CD, CS, MN, MM, ML, CP)

            Alpha = Alpha - DAA
            DCMDA = (MN - CM_TMP) / DAA

            IF (ABS (DCMDA) < SMALL) THEN
               WRITE (LWRIT, *) ' Derivative is too small: DCMDA = ', DCMDA
               STOP
            END IF

            DA    = -(MN - CM_TARG) / DCMDA
            DA    = MIN (MAX (-2.0, DA), 2.0)
            Alpha = Alpha + DA
            Alpha = MIN (MAX (-30., Alpha), 30.) ! Alpha = +176. observed!

         END DO

      END IF

!     Post-process the aerodynamics, ready for the objective & constraints:

      CALL POST_AERO

      IF (OUTPUTS) THEN
         CALL SECOND (CPU2)
         CALL CPUTIME(CPU1, CPU2, SUBNAME, 'to compute the aerodynamics', LWRIT)
      END IF

      TEMP_MAX  = ZERO ! For now
      HEAT_FLUX = ZERO

!     -------------------
!     Objective function:
!     -------------------

      CALL OBJECTIVE (OBJPAR, MCALL)

      IF (FIRST) THEN
         OBJINI             = OBJPAR
         AL_INI             = Alpha
         MA_INI             = M_inf
         CL_INI             = CL
         CD_INI             = CD
         LD_INI             = LD
         CM_INI             = MN
         TEMP_MAX_INI       = TEMP_MAX
         HEAT_FLUX_INI      = HEAT_FLUX
         SECTIONAL_AREA_INI = SECTIONAL_AREA
         SWET_INI           = SWET
         VOLUME_INI         = VOLUME
         V_EFF_INI          = V_EFF
      END IF

!     ---------------
!     Output results?
!     ---------------

      IF (OUTPUTS) THEN

         IF (NCLIN > 0) CALL CHKLCON ()      ! Review the linear constraints

         IF (NCNLN > 0) CALL NLCON (PRINT, CPRINT) ! & nonlinear constraints

         CALL SUMMARY (OBJPAR, MCALL) ! Tabulate; save two forms of geometry

      END IF

      IF      (MCALL == -3 ) THEN
         WRITE (LWRIT, 1080) 'Initial'
      ELSE IF (MCALL == -2 ) THEN
         WRITE (LWRIT, 1080) 'End-of-line-search'
      ELSE IF (MCALL == -1 ) THEN
         WRITE (LWRIT, 1080) 'Local search'
      ELSE IF (MCALL ==  0 ) THEN
         WRITE (LWRIT, 1080) 'Line search'
      END IF


      IF (MCALL == NDVPAR) THEN
         IF (.NOT. CONFUN_DOES_OBJ) THEN
            CALL SECOND (TIME2)
            CALL CPUTIME (CPUGRAD1, TIME2, SUBNAME, 'to calculate gradient',   &
                          LWRIT)
         END IF
      END IF

      IF (MCALL < 0) THEN
         CALL SECOND (TIME2)
         CALL CPUTIME (CPU0, TIME2, SUBNAME, 'for this run so far', LWRIT)
      END IF

      GO TO 999

!     Termination:

  900 STOP

  999 CALL FLUSH (LWRIT) ! Flush output buffer

      RETURN

!     Formats:

 1080 FORMAT (/, ' SOLVE: ', A, ' evaluation.')
 1300 FORMAT (/, ' SOLVE: ', (1X, A))

      END SUBROUTINE SOLVE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE SUMMARY (OBJ, MCALL)
!
!     Summarize the initial and optimized result and save geometry.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Global variables:

      USE AERO_MOD
      USE COEFFS_MOD
      USE CONST_MOD
      USE GEOM_MOD
      USE GRID_MOD
      USE OPT_MOD
      USE TRIANG_MOD

      IMPLICIT NONE

!     Arguments:

      REAL    :: OBJ
      INTEGER :: MCALL ! -4 means SOLVE is being called at the end of a run

!     Local variables:

      INTEGER :: I, K, N, NI1
      REAL    :: CPH, R, SPH, STOTAL, Y0, Z0

      REAL, ALLOCATABLE, DIMENSION (:) :: S, RS, RSS, XS, XSS, CURVATURE

!     Execution:

      WRITE (LWRIT, 1200) NDITER,            &
         OBJINI, OBJ,                        &
         AL_INI, Alpha,                      &
         MA_INI, M_inf,                      &
         CL_INI, CL,                         &
         CD_INI, CD,                         &
         CM_INI, MN,                         &
         LD_INI, LD,                         &
!!       TEMP_MAX_INI, TEMP_MAX,             &
!!       HEAT_FLUX_INI, HEAT_FLUX,           &
         SECTIONAL_AREA_INI, SECTIONAL_AREA, &
         SWET_INI, SWET,                     &
         VOLUME_INI, VOLUME,                 &
         V_EFF_INI, V_EFF

      WRITE (LWRIT, '(/, 10X, 2A)') &
         '   Xcenter   Ycenter   Zcenter     Xring     Rring', &
         '      X_cg      Y_cg      Z_cg'
      WRITE (LWRIT, '(1X, A, 8F10.5)') &
         'Initial: ', XICNTR, YICNTR, ZICNTR, XIRING, RIRING, XYZ_CG0, &
         'Current: ', XWCNTR, YWCNTR, ZWCNTR, XWRING, RWRING, XYZ_CG

!     Save the current geometry in surface mesh form:

      REWIND (LXYZ)
      WRITE  (LXYZ, '(I8)') 1
      WRITE  (LXYZ, '(3I8)') N_RAD_MAX, N_AZM_MAX, 1
      WRITE  (LXYZ, '(1P, 8E15.7)') XFIN, YFIN, ZFIN  ! Left-handed surface grid

!     Save the current defining sections in input geometry form:

      REWIND (LGSAVE)

      CALL WR_GEOM

!     Save the geometry in plottable form, and the grid in right-handed form?

      IF (MCALL == -4 .OR. NITMAX == 0) THEN

         REWIND (LXYZ)
         WRITE  (LXYZ, '(I8)') 1
         WRITE  (LXYZ, '(3I8)') N_RAD_MAX, N_AZM_MAX, 1
         WRITE  (LXYZ, '(1P, 8E15.7)') XFIN(:,N_AZM_MAX:1:-1), &  ! Right-handed
                                       YFIN(:,N_AZM_MAX:1:-1), &
                                       ZFIN(:,N_AZM_MAX:1:-1)
         CLOSE  (LXYZ)

!        For possible 2-D analysis of the centerline profile (HB_GRID gridding),
!        save the keel-crown line in wrap-around form:

         OPEN   (LKEEL,  FILE='centerline.dat', STATUS='UNKNOWN')

         WRITE  (LKEEL, '(A)') 'Keel-crown line for HB_GRID'
         WRITE  (LKEEL, '(I4)') 2*N_RAD_MAX - 1
         WRITE  (LKEEL, '(1P, 2E15.7)') &
            (XFIN(I,N_AZM_MAX), YFIN(I,N_AZM_MAX), I =    N_RAD_MAX, 1, -1),   &
            (XFIN(I,1),         YFIN(I,1),         I = 2, N_RAD_MAX)
         CLOSE  (LKEEL)

!        Plottable form of spokes, assuming common point counts:

         CLOSE  (LGSAVE)
         OPEN   (LGSAVE, FILE='heat_shield.spokes', STATUS='UNKNOWN')

         ALLOCATE (YO_PRTB(NT_RAD,NT_SEC), ZO_PRTB(NT_RAD,NT_SEC))

         IF (TWO_PARTS) THEN

            NI1 = NI_INIT(1)
            WRITE (LGSAVE, '(I1)') 2
            IF (FULL == 1) THEN
               WRITE (LGSAVE, '(3I6)') NI1, NI_SEC, 1
               WRITE (LGSAVE, '(3I6)') NO_INIT(1), NO_SEC, 1
            ELSE
               WRITE (LGSAVE, '(3I6)') NI1, 2*NI_SEC-1, 1
               WRITE (LGSAVE, '(3I6)') NO_INIT(1), 2*NO_SEC-1, 1
            END IF

            DO K = 1, NI_SEC
               CPH = COS_ZI(K)
               SPH = SIN_ZI(K)
               DO I = 1, NI_INIT(K)
                  R = RI_WORK(I,K) - RI_WORK(1,K)
                  YO_PRTB(I,K) = YWCNTR + R * CPH
                  ZO_PRTB(I,K) = ZWCNTR + R * SPH
               END DO
            END DO

            IF (FULL == 1) THEN
               WRITE (LGSAVE, '(1P, 8E15.7)') &
                  XI_WORK(1:NI1,1:NI_SEC),    &
                  YO_PRTB(1:NI1,1:NI_SEC),    &
                  ZO_PRTB(1:NI1,1:NI_SEC)
            ELSE
               WRITE (LGSAVE, '(1P, 8E15.7)') &
                  XI_WORK(1:NI1,1:NI_SEC), XI_WORK(1:NI1,NI_SEC-1:1:-1), &
                  YO_PRTB(1:NI1,1:NI_SEC), YO_PRTB(1:NI1,NI_SEC-1:1:-1), &
                  ZO_PRTB(1:NI1,1:NI_SEC),-ZO_PRTB(1:NI1,NI_SEC-1:1:-1)
            END IF

            DO K = 1, NO_SEC
               CPH = COS_ZO(K)
               SPH = SIN_ZO(K)
               Y0  = YO_PRTB(NI1,K)
               Z0  = ZO_PRTB(NI1,K)
               DO I = 1, NO_INIT(K)
                  R = RO_WORK(I,K) - RO_WORK(1,K)
                  YO_PRTB(I,K) = Y0 + R * CPH
                  ZO_PRTB(I,K) = Z0 + R * SPH
               END DO
            END DO

         ELSE

            NI1 = NO_INIT(1)
            WRITE (LGSAVE, '(I1)') 1
            IF (FULL == 1) THEN
               WRITE (LGSAVE, '(3I6)') NI1, NO_SEC, 1
            ELSE
               WRITE (LGSAVE, '(3I6)') NI1, 2*NO_SEC-1, 1
            END IF

            DO K = 1, NO_SEC
               CPH = COS_ZO(K)
               SPH = SIN_ZO(K)
               DO I = 1, NO_INIT(K)
                  R = RO_WORK(I,K) - RO_WORK(1,K)
                  YO_PRTB(I,K) = YWCNTR + R * CPH
                  ZO_PRTB(I,K) = ZWCNTR + R * SPH
               END DO
            END DO

!           Curvature distribution along center line (geometry data, not grid):

            N = NI1
            ALLOCATE (S(N), RS(N), RSS(N), XS(N), XSS(N), CURVATURE(N))

            CALL CHORDS2D (N, RO_WORK(1,1), XO_WORK(1,1), FALSE, STOTAL, S)

            DO I = 2, N - 1
               CALL FDCNTR (I, S, RO_WORK(1,1), RS(I), RSS(I))
               CALL FDCNTR (I, S, XO_WORK(1,1), XS(I), XSS(I))
            END DO

            CALL CURV2D (N - 2, RS(2), RSS(2), XS(2), XSS(2), CURVATURE(2))

!           Write the curvature in Tecplotable form, top to bottom.
!           Normalize arc lengths from nose to rim to facilitate comparisons.

            S(:) = S(:) / STOTAL

            OPEN   (LKEEL,  FILE='center_curvature.dat', STATUS='UNKNOWN')
            WRITE  (LKEEL, '(1P, 2E15.7)') (-S(I), CURVATURE(I), I = N-1, 2, -1)

            IF (.NOT. AXISYMMETRIC) THEN
               K = NO_SEC
               N = NO_INIT(K)

               CALL CHORDS2D (N, RO_WORK(1,K), XO_WORK(1,K), FALSE, STOTAL, S)

               DO I = 2, N - 1
                  CALL FDCNTR (I, S, RO_WORK(1,K), RS(I), RSS(I))
                  CALL FDCNTR (I, S, XO_WORK(1,K), XS(I), XSS(I))
               END DO

               CALL CURV2D (N - 2, RS(2), RSS(2), XS(2), XSS(2), CURVATURE(2))

               S(:) = S(:) / STOTAL
            END IF

            WRITE  (LKEEL, '(1P, 2E15.7)') (S(I), CURVATURE(I), I = 2, N-1)
            CLOSE  (LKEEL)

            DEALLOCATE (S, RS, RSS, XS, XSS, CURVATURE)

         END IF

         NI1 = NO_INIT(1)
         IF (FULL == 1) THEN
            WRITE (LGSAVE, '(1P, 8E15.7)') &
               XO_WORK(1:NI1,1:NO_SEC),    &
               YO_PRTB(1:NI1,1:NO_SEC),    &
               ZO_PRTB(1:NI1,1:NO_SEC)
         ELSE
            WRITE (LGSAVE, '(1P, 8E15.7)') &
               XO_WORK(1:NI1,1:NO_SEC), XO_WORK(1:NI1,NO_SEC-1:1:-1), &
               YO_PRTB(1:NI1,1:NO_SEC), YO_PRTB(1:NI1,NO_SEC-1:1:-1), &
               ZO_PRTB(1:NI1,1:NO_SEC),-ZO_PRTB(1:NI1,NO_SEC-1:1:-1)
         END IF

         CLOSE (LGSAVE)

      END IF

      RETURN

 1200 FORMAT (/, &
         ' PERFORMANCE PARAMETERS AT ITERATION', I4, /, 42X, &
         ' Initial', 12X, 'Current', 1P, /,             &
         ' Objective function            ', 2E19.10, /, &
         ' Alpha                   (deg.)', 2E19.10, /, &
         ' Mach number                   ', 2E19.10, /, &
         ' CL                            ', 2E19.10, /, &
         ' CD                            ', 2E19.10, /, &
         ' CM                            ', 2E19.10, /, &
         ' L/D                           ', 2E19.10, /, &
!!       ' Peak temperature      (deg. K)', 2E19.10, /, &
!!       ' Peak heat flux        (W/cm^2)', 2E19.10, /, &
         ' Cross-sectional area     (m^2)', 2E19.10, /, &
         ' Wetted area              (m^2)', 2E19.10, /, &
         ' Volume                   (m^3)', 2E19.10, /, &
         ' Effective volume coefficient  ', 2E19.10)

      END SUBROUTINE SUMMARY

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE SURFACE_GRID (FULL, NKG, NI, NK, I1, I2, RRG, XRG, ZG, R, X,  &
                               AZM, YCNTR, ZCNTR, CRV_TARG, CRV_MAX, CRV_TOT)
!
!     Spline interpolation in the azimuthal direction of heat shield sections
!     which have already been regularized.  All cases of inner and/or outer
!     shields and half/full geometries are handled.  The calling routine skips
!     the singular point at the apex and the repeated points at the ring if
!     the shield is in two parts.
!
!     NOTE:  The radial interpolation remains done at the higher level because
!            the arc lengths used are also needed for radial curvatures.
!
!     01/28/02  DAS  Initial argument-driven implementation, for two calls.
!     01/29/02   "   Added azimuthal curvature calculations.  (r,x,phi) do NOT
!                    work - had to evaluate y and use (y,z,s).  Nomenclature
!                    becomes confusing, with Z for both azimuth and z.
!     08/23/07   "   KAPPA now returns a maximum curvature across calls.
!                    It may not be usable in the azimuthal direction, since
!                    those curvatures are greatest near the apex (smallest r).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IMPLICIT NONE

!     Arguments:

      INTEGER, INTENT (IN) :: &
         FULL,                & ! 1 = whole shield; 0 = half shield
         NKG,                 & ! # geometry defining sections
         NI, NK,              & ! Radial & azimuthal dimensions of grid arrays
         I1, I2                 ! Output grid index range in radial direction

      REAL, INTENT (IN), DIMENSION (NI, NKG) :: &
         RRG, XRG               ! Regularized sections (cylindrical coordinates)

      REAL, INTENT (INOUT) :: &
         ZG(0:NKG+1)            ! Azimuthal angles of defining sections

      REAL, INTENT (INOUT), DIMENSION (NI, NK) :: &
         R, X                   ! Surface grid patch (cylindrical coordinates)

      REAL, INTENT (IN) ::    & ! Desired azimuthal grid stations
         AZM(NK),             &
         YCNTR, ZCNTR,        & ! Apex coordinates needed for curvature calcs.
         CRV_TARG               ! Lower bound, probably zero

      REAL, INTENT (INOUT) :: &
         CRV_MAX,             & ! Maximum curvature found so far in k direction
         CRV_TOT                ! Cumulative sum of (azimuthal crv. violatns.)^2

!     Local constants:

      LOGICAL, PARAMETER ::   &
         NEW = .TRUE., NORM = .FALSE.

      CHARACTER, PARAMETER :: &
         LOOSE * 1 = 'B'

!     Local variables:

      INTEGER :: &
         I, K, NPTS

      REAL :: &
         T

      REAL, DIMENSION (0:NKG+1) :: &
         RG, SG, XG, YG, Z      ! For gathering geometry points at a given i ...

      REAL, DIMENSION (NK) :: &
         RM, XM, DERIVS         ! ... and the resulting interpolations

!     Execution:

      NPTS = NKG + 2

      DO I = I1, I2

         DO K = 1, NKG
            T     = RRG(I,K) - RRG(1,1)
            YG(K) = YCNTR + T * COS (ZG(K))
            Z (K) = ZCNTR + T * SIN (ZG(K))
         END DO

         IF (FULL == 1) THEN

            RG(0) = RRG(I,NKG-1)
            XG(0) = XRG(I,NKG-1)
            YG(0) = YG(NKG-1)
            ZG(0) = -(ZG(NKG) - ZG(NKG-1))
            Z (0) = -Z(2)

            DO K = 1, NKG
               RG(K) = RRG(I,K)
               XG(K) = XRG(I,K)
            END DO

            RG(NKG+1) = RRG(I,2)
            XG(NKG+1) = XRG(I,2)
            YG(NKG+1) = YG(2)
            ZG(NKG+1) = ZG(NKG) + ZG(2)
            Z (NKG+1) = Z(2)

         ELSE ! Half shield

            RG(0) = RRG(I,2)
            XG(0) = XRG(I,2)
            YG(0) =  YG(2)
            ZG(0) = -ZG(2)
            Z (0) = -Z (2)

            DO K = 1, NKG
               RG(K) = RRG(I,K)
               XG(K) = XRG(I,K)
            END DO

            RG(NKG+1) = RRG(I,NKG-1)
            XG(NKG+1) = XRG(I,NKG-1)
            YG(NKG+1) = YG(NKG-1)
            ZG(NKG+1) = ZG(NKG) * 2. - ZG(NKG-1)
            Z (NKG+1) = -Z(NKG-1)
         END IF

         CALL LCSFIT (NPTS, ZG, RG, NEW, LOOSE, NK, AZM, RM, DERIVS)
         CALL LCSFIT (NPTS, ZG, XG, NEW, LOOSE, NK, AZM, XM, DERIVS)

         DO K = 1, NK
            R(I,K) = RM(K)
            X(I,K) = XM(K)
         END DO

         CALL CHORDS2D (NPTS, YG, Z, NORM, T, SG)
         CALL KAPPA    (NPTS, NPTS - 1, YG, Z, SG, CRV_TARG, CRV_MAX, CRV_TOT)

      END DO

      END SUBROUTINE SURFACE_GRID

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE TRIANGLE_PROPERTIES (NPTS, NFACES, IV, XYZ, SIGNS, &
                                      FACE_AREA, XYZ_NORM, XYZ_CENTROID, SWET)
!
!     Calculate areas, components of unit normals, and centroids for all
!     faces of a triangulated surface, as needed for calculating Newtonian
!     aerodynamic coefficients.
!
!     08/25/00  J. Reuther   Original implementation.
!     10/18/00  D. Saunders  Argument-driven version.
!     12/13/00     "   "     Return FACE_AREA(*), not PROJECTED_AREA(1:3,*).
!     10/23/01     "   "     Moved wetted area here (James had it in NEWTONIAN).
!
!     Sponsor:  NASA Ames Research Center, Mountain View, CA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IMPLICIT NONE

!     Arguments:

      INTEGER, INTENT (IN) :: &
         NPTS,                &  ! Total mumber of points in XYZ(3,*)
         NFACES,              &  ! Number of triangles
         IV(3,*)                 ! Pointers to the three vertices of
                                 ! each triangle in XYZ(1:3,*)
      REAL, INTENT (IN) ::    &
         XYZ(3,NPTS),         &  ! Surface point coordinates
         SIGNS(NFACES)           ! +1. for outward normal, else -1.

      REAL, INTENT (OUT) ::   &
         FACE_AREA(NFACES)       ! Triangle areas >= 0.

      REAL, INTENT (OUT), DIMENSION (3,NFACES) :: &
         XYZ_NORM,            &  ! Components of triangle outward unit normals
         XYZ_CENTROID            ! Triangle centroids

      REAL, INTENT (OUT) ::   &
         SWET                    ! Wetted area

!     Local constants:

      REAL, PARAMETER :: &
         HALF = 0.5, ONE = 1., SMALL = 1.E-30, THIRD = 1./ 3., ZERO = 0.

!     Local variables:

      INTEGER :: I, I1, I2, I3
      REAL    :: AREAX, AREAY, AREAZ, RAREA, SWITCH

!     Execution:

      SWET = ZERO

      DO I = 1, NFACES

         I1 = IV(1,I)
         I2 = IV(2,I)
         I3 = IV(3,I)
         SWITCH = HALF * SIGNS(I)

         AREAX = ((XYZ(2,I2) - XYZ(2,I1))  * &
                  (XYZ(3,I3) - XYZ(3,I1))  - &
                  (XYZ(3,I2) - XYZ(3,I1))  * &
                  (XYZ(2,I3) - XYZ(2,I1))) * SWITCH
         AREAY = ((XYZ(3,I2) - XYZ(3,I1))  * &
                  (XYZ(1,I3) - XYZ(1,I1))  - &
                  (XYZ(1,I2) - XYZ(1,I1))  * &
                  (XYZ(3,I3) - XYZ(3,I1))) * SWITCH
         AREAZ = ((XYZ(1,I2) - XYZ(1,I1))  * &
                  (XYZ(2,I3) - XYZ(2,I1))  - &
                  (XYZ(2,I2) - XYZ(2,I1))  * &
                  (XYZ(1,I3) - XYZ(1,I1))) * SWITCH

         FACE_AREA(I)  = SQRT (AREAX ** 2 + AREAY ** 2 + AREAZ ** 2)
         SWET          = SWET + FACE_AREA(I)
         RAREA         = ONE / MAX (FACE_AREA(I), SMALL)

         XYZ_NORM(1,I) = AREAX * RAREA ! Outward unit normal components
         XYZ_NORM(2,I) = AREAY * RAREA
         XYZ_NORM(3,I) = AREAZ * RAREA

         XYZ_CENTROID(1,I) = (XYZ(1,I1) + XYZ(1,I2) + XYZ(1,I3)) * THIRD
         XYZ_CENTROID(2,I) = (XYZ(2,I1) + XYZ(2,I2) + XYZ(2,I3)) * THIRD
         XYZ_CENTROID(3,I) = (XYZ(3,I1) + XYZ(3,I2) + XYZ(3,I3)) * THIRD

      END DO

      END SUBROUTINE TRIANGLE_PROPERTIES

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE TRIANGULATE_PATCHES (NPATCH, NIJ, XYZ, IP, ICOMPONENT, &
                                      SWITCH_NORM, REF_LEN,             &
                                      NPTS, NFACES, IV, ICOMP, SIGNS)
!
!     One-liner:  Triangulate structured surface patches.
!
!     TRIANGULATE_PATCHES converts a surface paneling (one or more structured
!     surface patches) to unstructured form.  The shorter diagonal of each
!     quadrilateral cell is used to define each pair of triangles.  Account
!     is taken of collapsed cell edges: all output triangles have nonzero
!     area.  However, no check is made for non-convex cells.
!
!     The patch points are expected to be in triplet form (possibly via
!     appropriate reading of an AEROSURF/PLOT3D form).  Thus there is no
!     need to copy (x,y,z)s.
!
!     10/12/00  DAS  Initial adaptation of Mark Rimlinger's Struct2Unstruct
!                    routine from UV_MAP.
!
!     Author:  David Saunders, ELORET/NASA Ames Research Center, Mtn. View, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IMPLICIT NONE

!     Arguments:

      INTEGER, INTENT (IN) ::  &
         NPATCH                  ! Number of patches

      INTEGER, INTENT (IN) ::  &
         NIJ(2,NPATCH)           ! Patch dimensions

      REAL, INTENT (IN) ::     &
         XYZ(3,*)                ! Surface point coordinates (packed)

      INTEGER, INTENT (IN) ::  &
         IP(NPATCH),           & ! Indices in XYZ(1:3,*) of the start of
         ICOMPONENT(NPATCH)      ! each patch, and geometry component #s

      REAL, INTENT (IN) ::     &
         SWITCH_NORM(NPATCH),  & ! +/-1. handles or right/left-handed patches
         REF_LEN                 ! Reference length: cell edges less than
                                 ! REF_LEN * 1.E-15 are considered to be
                                 ! collapsed
      INTEGER, INTENT (OUT) :: &
         NPTS                    ! Total number of points in XYZ(3,*):
                                 ! all must be included if the triangulation
                                 ! is written to disk in FAST format
      INTEGER, INTENT (OUT) :: &
         NFACES                  ! Number of triangles (tet. faces) found
                                 ! (not counting collapsed triangles)
      INTEGER, INTENT (OUT) :: &
         IV(3,*)                 ! Pointers to the three vertices of
                                 ! each triangle in XYZ(1:3,*)
      INTEGER, INTENT (OUT) :: &
         ICOMP(*)                ! Flag for each triangle, set to the
                                 ! appropriate ICOMPONENT element
      REAL, INTENT (OUT) :: &
         SIGNS(*)                ! Flag for each triangle derived from
                                 ! the SWITCH_NORM flag for its patch
!     Local constants:

      REAL, PARAMETER :: &
         ZERO = 0.

!     Local variables:

      INTEGER :: &
         I, I0, I1, I2, I3, I4, IC, J, J0, L, N, NI, NJ, NTRI

      REAL :: &
         DIAG1, DIAG2, EDGE1, EDGE2, EDGE3, EDGE4, SN, TOL

!     Execution:

      TOL  = 1.E-30 * (REF_LEN ** 2)
      N    = NPATCH
      NPTS = IP(N) + NIJ(1,N) * NIJ(2,N) - 1 ! Possibly IP(NPATCH+1)
      NTRI = 0

      DO N = 1, NPATCH

         I0 = IP(N) - 1 ! Offset in XYZ(1:3,*) for (i,j) elements of patch N
         IC = ICOMPONENT(N)
         NI = NIJ(1,N)
         NJ = NIJ(2,N)
         SN = SWITCH_NORM(N)

         DO J = 1, NJ - 1

            J0 = I0 + (J - 1) * NI

            DO I = 1, NI - 1 ! For each quad. cell ...

               I1 = J0 + I  ! (i,j)
               I2 = I1 + 1  ! (i+1,j)
               I3 = I2 + NI ! (i+1,j+1)
               I4 = I3 - 1  ! (i,j+1)

               DIAG1 = ZERO
               DIAG2 = ZERO
               EDGE1 = ZERO
               EDGE2 = ZERO
               EDGE3 = ZERO
               EDGE4 = ZERO

               DO L = 1, 3
                  DIAG1 = DIAG1 + (XYZ(L,I1) - XYZ(L,I3)) ** 2
                  DIAG2 = DIAG2 + (XYZ(L,I2) - XYZ(L,I4)) ** 2
                  EDGE1 = EDGE1 + (XYZ(L,I1) - XYZ(L,I2)) ** 2
                  EDGE2 = EDGE2 + (XYZ(L,I2) - XYZ(L,I3)) ** 2
                  EDGE3 = EDGE3 + (XYZ(L,I3) - XYZ(L,I4)) ** 2
                  EDGE4 = EDGE4 + (XYZ(L,I4) - XYZ(L,I1)) ** 2
               END DO

               IF (DIAG1 <= DIAG2) THEN

                  IF (EDGE1 > TOL .AND. EDGE2 > TOL) THEN
                     NTRI = NTRI + 1
                     IV(1,NTRI)  = I1
                     IV(2,NTRI)  = I2
                     IV(3,NTRI)  = I3
                     ICOMP(NTRI) = IC
                     SIGNS(NTRI) = SN
                  END IF

                  IF (EDGE3 > TOL .AND. EDGE4 > TOL) THEN
                     NTRI = NTRI + 1
                     IV(1,NTRI)  = I3
                     IV(2,NTRI)  = I4
                     IV(3,NTRI)  = I1
                     ICOMP(NTRI) = IC
                     SIGNS(NTRI) = SN
                  END IF

               ELSE ! DIAG1 > DIAG2

                  IF (EDGE2 > TOL .AND. EDGE3 > TOL) THEN
                     NTRI = NTRI + 1
                     IV(1,NTRI)  = I2
                     IV(2,NTRI)  = I3
                     IV(3,NTRI)  = I4
                     ICOMP(NTRI) = IC
                     SIGNS(NTRI) = SN
                  END IF

                  IF (EDGE4 > TOL .AND. EDGE1 > TOL) THEN
                     NTRI = NTRI + 1
                     IV(1,NTRI)  = I4
                     IV(2,NTRI)  = I1
                     IV(3,NTRI)  = I2
                     ICOMP(NTRI) = IC
                     SIGNS(NTRI) = SN
                  END IF

               END IF

            END DO ! Next I

         END DO ! Next J

      END DO ! Next patch

      NFACES = NTRI

      END SUBROUTINE TRIANGULATE_PATCHES

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE WR_GEOM
!
!     Write the geometry of a general heat shield in cylindrical coordinates
!     as normalized sections.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      USE AERO_MOD
      USE COEFFS_MOD
      USE CONST_MOD
      USE GEOM_MOD

      IMPLICIT NONE

!     Local constants:

      REAL, PARAMETER :: ONE = 1.

!     Local variables:

      INTEGER :: I, K
      REAL    :: AZ_FRACTION, RSTAT, XTEMP(NT_RAD), RTEMP(NT_RAD)

!     Execution:

      WRITE (LGSAVE, '(A)') &
         'HEAT_SHIELD output defining sections', &
         '     S_ref     L_ref      X_cg      Y_cg      Z_cg    S_Base'
      WRITE (LGSAVE, '(6F10.5)') REF_AREA, REF_LEN, XYZ_CG, SBASE
      WRITE (LGSAVE, '(A)') '      FULL    NI_SEC    NO_SEC   NON-DIM'
      WRITE (LGSAVE, '(4I10)') FULL, NI_SEC, NO_SEC, 1
      WRITE (LGSAVE, '(A)') &
         '     XWCNTR     YWCNTR     ZWCNTR     XIRING     RIRING'
      WRITE (LGSAVE, '(5F11.6)') XWCNTR, YWCNTR, ZWCNTR, XWRING, RWRING

      IF (TWO_PARTS) THEN

         WRITE (LGSAVE, '(A)') 'INNER RADIAL SECTIONS'

         DO K = 1, NI_SEC
            WRITE (LGSAVE, '(A, I4)') 'Inner section', K
            WRITE (LGSAVE, '(A)') '   NI_INIT    ZI_INIT      CI_INIT'
            AZ_FRACTION = ZI_WORK(K) * RTOT_AZ
            WRITE (LGSAVE, '(I10, F11.4, F13.8)') &
               NI_INIT(K), AZ_FRACTION, CI_WORK(K)
            WRITE (LGSAVE, '(A)') '     RI_INIT      XI_INIT'
            WRITE (LGSAVE, '(2F13.8)') &
               (RI_WORK(I,K), XI_WORK(I,K), I = 1, NI_INIT(K))
         END DO

      END IF

      WRITE (LGSAVE, '(A)') 'OUTER RADIAL SECTIONS'

      DO K = 1, NO_SEC

         WRITE (LGSAVE, '(A, I4)') 'Outer section', K
         WRITE (LGSAVE, '(A)') &
            '   NO_INIT    ZO_INIT      CO_INIT      VO_INIT'
         AZ_FRACTION = ZO_WORK(K) * RTOT_AZ
         WRITE (LGSAVE, '(I10, F11.4, 2F13.8)') &
            NO_INIT(K), AZ_FRACTION, CO_WORK(K), VO_WORK(K)
         WRITE (LGSAVE, '(A)') '     RO_INIT      XO_INIT'
         WRITE (LGSAVE, '(2F13.8)') &
            (RO_WORK(I,K), XO_WORK(I,K), I = 1, NO_INIT(K))

      END DO

      END SUBROUTINE WR_GEOM
