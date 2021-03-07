
      MODULE CONSTS       ! *** Start of Traj_opt modules *** !

      INTEGER, PARAMETER ::
     >   LREAD    = 1,    ! traj_opt.inp
     >   LWRIT    = 2,    ! traj_opt.out
     >   LSAVE    = 3,    ! CL/CD sensitivities (if MODE_VAR = 0)
     >   LRESULTS = 4,    ! traj_opt.qplot (alpha and bank results)
     >   LNPOUT   = 6,    ! NPOPT output = standard output
     >   LNPOPT   = 7,    ! traj_opt.nps (NPOPT control input summary)
     >   LDAT     = 8,    ! APC dataset, ascent trajectory, etc.
     >   N_ENTRY  = 6,    ! # trajectory entry conditions provided for
     >   N_MAX    = 11    ! # maximum values provided for

      REAL
     >   DEG2RAD, EPSMCH, PI, RAD2DEG, RPI, RSTEFAN

      END MODULE CONSTS


      MODULE OPT

      INTEGER, PARAMETER ::
     >   MAXARGS  = 8,    ! Limit on the number of command line arguments
     >   MAXCHARS = 24    ! Limit on the length of a command line argument

******** N.B.:  MAXARGS and MAXCHARS are hard-coded in trajectory.c ************
*                   (Also in the traj2.c and traj3.c variants.)                *

      INTEGER
     >   NNOPT,  NNDIFF, NNGRAD, NNTIME,  NNSPLINE, 
     >   NALLOW, NARGS,  NDITER, NFTOTL,  NITMAX,
     >   NDV,    NCLIN,  NCNLN,  NBNDS,
     >   NROWA,  NROWJ,  NROWR,  LENCW,   LENIW,   LENW,
     >   LINCON,     NLINCON_RD, MODE_VAR,
     >   NJOURNAL,   NSTEPS,     LAST_STEP, MAX_STEP,
     >   NDV_ALPHA,  NDV_BANK,   NDV_BETA,  NDV_ANGLES,
     >   NDV_CL,     NDV_CD,     NDV_CY,    NDV_COEFFS,
     >   NK_ALPHA,   NK_BANK,    NK_BETA,   KNOT1_A,   KNOT1_B,
     >   NK_CL,      NK_CD,      NK_CY,
     >   N_MACH,     N_ALPHA,    N_BETA,    N_RE,      N_COEFS,
     >   N_APC,      N_ASCENT,   N_UPPER_C, N_TARGET,
     >   N_MACH_DH,  N_ALPHA_DH, N_BETA_DH, N_QBAR_DH,
     >   N_SURF_POINTS, N_SURF_TYPES,
     >   LABORT,     LEN_VEHICLE

*     The following equivalences avoid new control inputs that will seldom be
*     used.  Only one set of these 4 variables can be used in a given run.

      EQUIVALENCE
     >  (NDV_COEFFS, NDV_ANGLES),
     >  (NDV_CL, NDV_ALPHA),  (NDV_CD, NDV_BANK),  (NDV_CY, NDV_BETA)

      INTEGER, ALLOCATABLE, DIMENSION (:) ::
     >   IERROR, IGROUP, ILCON,  INLCON, JNLCON, ISTATE, NFDEL, IWVEC

      LOGICAL
     >   CONFUN_DOES_OBJ, DISTRIBUTED_HEATING, NORMALIZE,
     >   OPTIMIZE_TABORT,
     >   SHOW_APOAPS, SHOW_ECCENT, SHOW_INCLIN, SHOW_SHORTFALL

      REAL
     >   APOAPSIS,      APOAPSIS_INI,
     >   CROSS_RANGE,   CROSS_RANGE_INI,
     >   DOWN_RANGE,    DOWN_RANGE_INI,
     >   END_ALPHA,     END_ALPHA_INI,
     >   END_ALTITUDE,  END_ALTITUDE_INI,
     >   END_BANK,      END_BANK_INI,
     >   END_FPA,       END_FPA_INI,
     >   END_HEADING,   END_HEADING_INI,
     >   END_LATITUDE,  END_LATITUDE_INI,
     >   END_LONGITUDE, END_LONGITUDE_INI,
     >   END_VELOCITY,  END_VELOCITY_INI,
     >   HEAT_LOAD,     HEAT_LOAD_INI,
     >   HEATLD,        HEATLD_INI,
     >   ORBITAL_ECC,   ORBITAL_ECC_INI,
     >   ORBITAL_INC,   ORBITAL_INC_INI,
     >   SHORTFALL,     SHORTFALL_INI,
     >   TABORT,        TABORT_INI,
     >   TOTAL_TIME,    TOTAL_TIME_INI

      REAL
     >   CPU0, ENTRY_DV, EPSOBJ, EPSM6TH, GTP, GNM,
     >   OBJBND, OBJ_CONFUN, OBJ_INI, PERCENT_LOVERD,
     >   PITCH_RATE,       ROLL_RATE,        YAW_RATE,
     >   RHO_ALPHA_TVD,    RHO_BANK_TVD,     RHO_BETA_TVD,
     >   RHO_CL_TVD,       RHO_CD_TVD,       RHO_CY_TVD,
     >   RHO_CROSS_RANGE,  RHO_DOWN_RANGE,   RHO_DURATION,
     >   RHO_HEAT_LOAD,    RHO_HEATLD,       RHO_LST_SQRS,
     >   RHO_MAX_ACCEL,    RHO_MAX_W_TEMP,   RHO_MAX_DYN_PR,
     >   RHO_MAX_H_FLUX,   RHO_TABORT,       RHO_END_ALTITUDE,
     >   RHO_ECCENTRICITY, RHO_TARGET_POINT, RHO_TRIM,
     >   TARGET_LATITUDE,  TARGET_LONGITUDE, TARGET_TRIM

      EQUIVALENCE
     >  (RHO_CL_TVD, RHO_ALPHA_TVD), (RHO_CD_TVD, RHO_BANK_TVD),
     >  (RHO_CY_TVD, RHO_BETA_TVD)

      REAL
     >   ORBITAL(0:10)

      REAL, ALLOCATABLE, DIMENSION (:) ::
     >   ASCENT_TIME,    ASCENT_VELOCITY,
     >   ASCENT_GAMMA,   ASCENT_HEADING,
     >   ASCENT_ALTITUDE,ASCENT_LATITUDE, ASCENT_LONGITUDE,
     >   FFWD,   GRAD,   PWORK,  SH,     URH,    V,
     >   YWORK,  ZWORK,  CVEC,   WVEC,
     >   VSCALE, AITCH,  BL,     BU,     CONLIN, TLCON,  C, CPRINT,
     >   BLLIN,  BULIN,  XNLCON, SNLCON, OLDG,   PSRCH,
     >   DEE,    ELL,    CLAMDA,
     >   ALPHA_STEPS,    BANK_STEPS,     DURATION,
     >   ALPHA_K,        BANK_K,         BETA_K,
     >   TK_ALPHA,       TK_BANK,        TK_BETA,
     >   CL_K,           CD_K,           CY_K,
     >   TK_CL,          TK_CD,          TK_CY,
     >   V_APC,          H_APC,          V_TRJ,         H_TRJ,
     >   V_UPPER_C,      H_UPPER_C,
     >   MAX_VALUES,     MAX_VALUES_INI,
     >   MAX_TIMES,      MAX_TIMES_INI,  ENTRY_STATE,   FINAL_STATE,
     >   MACH, ALPHA, BETA, RE, CL, CD, CJACBAK,
     >   MACH_DH, ALPHA_DH, BETA_DH, QBAR_DH,
     >   SURF_HEAT_LOADS, SURF_HOT_AREAS

      REAL, ALLOCATABLE, DIMENSION (:,:) ::
     >   AMAT, RMAT, CJAC, SURF_BOUNDS, SURF_PEAKS, SURF_VALUES,
     >   TRAJ_STEPS, TRAJ_TARGET

      REAL, ALLOCATABLE, DIMENSION (:,:,:) ::
     >   SURF_QUANTITIES, SURF_VIOLATIONS   

      REAL, ALLOCATABLE, DIMENSION (:,:,:,:,:) ::
     >   ATDB

      LOGICAL, ALLOCATABLE, DIMENSION (:) ::
     >   LFAIL

      CHARACTER
     >   AERO_FILENAME*48, ARGS(0:MAXARGS)*(MAXCHARS), VEHICLE*32

      CHARACTER, ALLOCATABLE, DIMENSION (:) ::
     >   CW*8, VTYPE*6, LCTYPE*6, NLCTYPE*6

      END MODULE OPT


********************************************************************************
*
      PROGRAM TRAJ_OPT
*
*        Traj_opt performs trajectory optimization with respect to variables
*     such as bank angle and angle of attack, using the Traj analysis package
*     and the NPOPT nonlinear programming package.  The quantity being optimized
*     may be some measure of cross range, down range, peak acceleration, etc.,
*     or some weighted combination of these.  Constraints may be imposed on
*     similar quantities.
*
*        This program does not treat trajectory optimization with respect to the
*     shape of the vehicle.  However, an option to calculate sensitivities with
*     respect to the aerodynamic database of lift and drag coefficients is
*     included.  This can facilitate the choice of flight conditions at which
*     multipoint shape optimization might be performed for a given trajectory.
*
*        Nonzero sideslip and side force coefficient are treated by this
*     version of Traj/Traj_opt in response to the loss of the Columbia orbiter.
*
*     APPROACH:  (Fixed-geometry/unpowered case/normal usage)
*
*        For a vehicle with fixed geometry, the variables for adjusting an
*     unpowered trajectory are normally the bank angle and angle of attack,
*     as functions of time.  Originally, Traj_opt had to manipulate these angles
*     via the Traj maneuver script mechanism.  Now, Traj has the option to
*     interpolate bank and alpha at every time step.  Traj also reads the
*     control points (spline knots) from traj.alp and traj.bnk (originally done
*     in Traj_opt).  Traj_opt merely updates these control points using their
*     ordinates as the main optimization variables.
*
*        Any of the Traj command line arguments may be used as arguments to
*     Traj_opt (at least under SGI IRIX; IARGC and GETARG system utilities may
*     be missing from other operating systems).
*
*        Since the total trajectory time is unknown a priori, the tail ends of
*     the bank and alpha schedules may not be well optimized.  However, the time
*     ranges of the variables should normally exceed the trajectory duration;
*     the last few variables for each schedule are "wasted" in this case.
*
*        If the optimized time interval t is short of a trajectory time T, the
*     relevant angle is held constant from time t:T.  Normally, starting with
*     t comfortably greater than T(final) is recommended.
*
*        The initial part of the trajectory may be frozen by specifying the
*     indices of the first variable control points for alpha and bank to be
*     greater than 1 - see KNOT1_A/B.  Keep in mind that the "local" spline
*     technique employed within Traj uses 4-point formulas.
*
*        A new trajectory is calculated for each change of variable.  Quantities
*     used for the objective function and constraints may be those from the last
*     time step of the trajectory, or from some fixed time step.  (See NNTIME.)
*     [Later: Traj now interpolates between time steps to capture the terminal
*     condition precisely, making gradients much more reliable.]
*
*        Gradient information is calculated by finite differencing.  [An adjoint
*     method for much-improved efficiency is very unlikely to be incorporated.]
*
*        Constraints imposed on the alpha and/or bank angles are simple bounds -
*     the same BL and BU for all angle of attack control points and a second
*     pair of bounds for all bank control points.  Other constraints such as on
*     maximum deceleration are nonlinear.  [Later:  Crude control of roll rate
*     and pitch rate has been implemented as linear constraints on the alpha and
*     bank control points.]
*
*        Given that NPOPT always calls CONFUN before OBJFUN, duplicate function
*     evaluations are avoided here by having CONFUN also calculate the objective
*     function.  See NNGRAD = 2 or 3.
*
*        Use of an APC (aerothermal performance constraint) is appropriate for
*     controlling heating at one leading edge point (presumably a sharp leading
*     edge).  Traj also has the option to perform stagnation point heating
*     calculations for blunt noses.  Normally, other key points on the surface
*     should also have their temperatures constrained via the distributed
*     heating capability within Traj_opt, more on which below.  Constraints on
*     dynamic pressure and G load are generally recommended also.
*
*        Maximum down-range or minimum heat-load calculations should benefit
*     from the "UPPERC" constraint (originally "CONSTQ", meaning a curve in
*     velocity/altitude space corresponding to a low dynamic pressure such as
*     10 psf).  This defines a meaningful upper corridor boundary with at least
*     some control surface effectiveness.  The "OSCILL" constraint may also
*     help avoid skipping or other altitude oscillations, although oscillations
*     aren't necessarily bad.  Excessive roll reversals likely for certain
*     types of trajectory can be inhibited via a TVD (total-variation-diminish-
*     ing) term available for the objective.  A further TVD term can also
*     encourage smoothness in the solution for the angle of attack schedule.
*     [Likewise for Beta (MODE_VAR = 5) and CL, CD, CY (MODE_VAR = 6).]
*
*        Heat load and heat flux calculations by Traj are strictly for the
*     stagnation point, suitable for blunt bodies.  More thorough treatment
*     of surface heating is provided by the constraints SURF1[R], SURF2, ...
*     and an aerothermal database containing m quantities at each of n key
*     surface points as functions of Mach, angle of attack, and dynamic
*     pressure.
*
*        The time history of these three independent variables as stored in
*     Traj's "journal" block of output time steps is used to interpolate the
*     aerothermal database at the Traj_opt level (as opposed to within Traj).
*     From these interpolations, constraint values are determined either via
*     peaks over time (as for peak temperature), or via integrations (as for
*     heat load derived from the heat flux history at each point).  Actually,
*     this version integrates elements of area above the bound for each point
*     rather than just the peak violations.
*
*        The order of the aerothermal database quantities (which may grow) is:
*
*     1: temperature  2: heat flux ...  m: ???
*
*        The order of the constraints is:
*
*     SURF1: distrib. temperature  SURF2: distrib. heat flux  ...  SURFm: ???
*
*        Rather than a distinct nonlinear constraint at each surface point,
*     a total violation of the constraints over all surface points is used,
*     for each quantity implied by SURF1, SURF2, ...  The upper bound (on
*     total violation) should therefore be zero.  Typically, only one or two
*     of the points touch their limits, so summing over points is not as
*     crude as it sounds, and certainly keeps the constraint nomenclature
*     more manageable.
*
*        A meaningful integrated heat load over the surface requires an area-
*     based weight associated with each surface point.  Associated constraint
*     or objective function calculations are invoked via HEATLD or RHO_HEATLD.
*     <This still needs thought.>
*
*        These calculations require possibly different temperature/heat flux
*     limits at each surface point along with a point weight.  Therefore, if a
*     SURFj or HEATLD constraint is present, or RHO_HEATLD is nonzero, a file
*     of point weights and upper limits should be supplied as follows:
*
*     "traj_opt.surface_bounds" format:
*
*        Title
*        m   n                  ! # quantities at each point; # points
*        W  BU1  BU2  ...  BUm  ! Wt. + upper bounds for m quantities at point 1
*        W  BU1  BU2  ...  BUm  ! .......................................point 2
*        .....................  ! ..............................................
*        .....................  ! .......................................point n
*
*     Note:  It has been belatedly realized that temperature and heat flux, as
*            typically calculated for these tables, are not independent.
*            Therefore, constraining both of them is redundant.
*            Moreover, trilinear interpolation of temperature in the presence
*            of insufficient data points (particularly in the dynamic pressure
*            direction - HAVOC is limited to 5 per run) is not as accurate as
*            deriving temperature from interpolated heat flux.  The retrofitted
*            SURF1R constraint (using heat-Rate-derived temperatures) therefore
*            expects surface emissivities in the BU2 column above.  This avoids
*            a separate input.
*
*        The aerothermal database is rectangular with format as follows:
*
*     "aerothermal.database" format:  (Units should be mks, although
*                                      Qbar is shown in psi.)
*        Title
*        N_MACH_DH             ! # Mach numbers
*        N_ALPHA_DH            ! # angles of attack
*        N_QBAR_DH             ! # dynamic pressures
*        N_SURF_TYPES  (m)     ! # surface data types at each pt in the database
*        N_SURF_POINTS (n)     ! # surface points in the database
*        SURFACE POINT 1
*        Mach      Alpha     Qbar     Temp   [Qdot ....]    ! Column headers
*        2.        0.        100.     xxx     xxx           ! Mach varies first,
*        3.        0.        100.     xxx     xxx           ! then Alpha, then
*        4.        0.        100.     xxx     xxx           ! dynamic pressure
*        :         :         :        :       :
*        25.       45.       500.     xxx     xxx
*        SURFACE POINT 2
*        Mach      Alpha     Qbar     Temp   [Qdot ....]    ! Column headers
*        2.        0.        100.     xxx     xxx
*        3.        0.        100.     xxx     xxx
*        4.        0.        100.     xxx     xxx
*        :         :         :        :       :
*        25.       45.       500.     xxx     xxx
*        :         :         :        :       :
*        :         :         :        :       :
*
*
*     ASCENT ABORTS:
*
*        Calculation of ascent abort trajectories prompted introduction of a
*     target landing site (latitude, longitude) and Traj input script options
*        Run_to   Impact
*        Floor    10.
*     which allow reaching the site at the specified altitude (10 km here).
*
*        Optimizing such abort trajectories (as opposed to use of trial-and-
*     error starting points) requires reading booster launch trajectory data
*     and treating abort time as an optimization variable.  See the TABORT
*     variable description below for the ascent data format.
*
*     AEROCAPTURE:
*
*        Aero-capture applications prompted the transfer of orbital parameters
*     from Traj and use of the "run to apoapsis" option, initially using the
*     eccentricity as an objective with constraints on apoapsis & inclination.
*
*        Aero-capture also prompted the step-function option for Alpha and bank
*     (MODE_VAR = 3 or 4).  The format for entering the angle control points in
*     this form via "traj_opt.steps" is as follows:
*
*        Alpha & bank steps
*        4  !    Alpha          Bank  Duration (s)
*        28.000000E+00 -5.286404E+01  4.292902E+01
*        28.000000E+00 -4.645406E+01  3.505116E+01
*        28.000000E+00  1.314490E+02  7.209257E+01
*        28.000000E+00  8.619661E+01  5.000000E+01
*
*        Optimized steps are output as "traj_opt.steps_opt" in this case.
*
*     COMMAND LINE OPTIONS:
*
*        See the descriptions in Traj - the same mechanism is employed here.
*     However, Linux and Windows-based systems appear to lack IARGC and GETARG
*     system utilities, so all Traj options may have to be entered in "traj.in".
*
*     CONTROL INPUT DESCRIPTION (traj_opt.inp):
*
*        NDV        Number of "design" (optimization) variables (total);
*                   alpha variables normally precede bank variables;
*                   others may follow (see MODE_VAR)
*        NDV_ALPHA  Number of variables affecting angle of attack vs. time;
*        (NDV_CL)   needs "Schedul Fnct_alpha_..." (or Fnct_cl_...) in traj.in;
*                   if MODE_VAR = 3 or 4 (step functions), enter NSTEPS;
*                   if MODE_VAR = 6, enter NDV_CL here
*        NDV_BANK   Number of variables affecting bank angle vs. time;
*        (NDV_CD)   needs "Schedul Fnct_bank_..."  (or Fnct_cd_...) in traj.in;
*                   if MODE_VAR = 3 or 4 (step functions), enter NSTEPS;
*                   if MODE_VAR = 6, enter NDV_CL here
*        NDV_BETA   Number of variables affecting sideslip angle vs. time;
*        (NDV_CY)   needs "Schedul Fnct_beta_..."  (or Fnct_cy_...) in traj.in;
*                   MODE_VAR must be 5 for nonzero Beta;
*                   Beta variables follow Alpha and bank variables in V(*);
*                   if MODE_VAR = 6, enter NDV_CY here
*        MODE_VAR   Controls optimization variable choice:
*                      0 = CLs, CDs of aero. database, for sensitivities only;
*                          NDV = 2 * N_MACH * N_ALPHA * N_RE; see also H_ALPHA;
*                      1 = Alpha and/or bank only (fixed control point times);
*                          NDV = NDV_ALPHA + NDV_BANK;
*                      2 = As for 1 but provides for variable entry conditions,
*                          which are required to be inertial;
*                          NDV = NDV_(ALPHA + BANK) allows use of ENTRY_DV > 0;
*                          NDV = NDV_(ALPHA + BANK) + 1 means an ascent abort;
*                      3 = Alpha & bank are step functions (fixed durations):
*                          alpha & bank are assumed to step at the same times;
*                          NDV = 2 * NSTEPS in "traj_opt.steps";
*                      4 = Alpha & bank are step functions (variable durations);
*                          NDV = 3 * NSTEPS in "traj_opt.steps" - 1
*                      5 = Alpha + bank + Beta control points at fixed times;
*                          nonzero Beta is intended to allow calculation of a
*                          trajectory which minimizes heating at some damaged
*                          point while still constraining heating elsewhere
*                      6 = CL + CD + CY control points at fixed times;
*                          intended to be used with known schedules for Alpha,
*                          bank and sideslip (as from Columbia tracking data)
*                          to match the known trajectory as well as possible
*                          to estimate the side forces experienced and hence
*                          glean some information about the likely damage
*
*                   The aero. database dimensions should be outputs from Traj
*                   initialization, but are instead inputs to Traj_opt:
*
*        N_MACH     Dimensions of ...
*        N_ALPHA    ... the aerodynamic database ...
*        N_RE       ... which is a rectangular table in (M,AoA,log(Re)) space
*        N_BETA     ... or (M,AoA,log(Re),Beta) space if MODE_VAR = 5
*        PERCENT_LOVERD  Mechanism for checking sensitivities to L/D; enter a
*                        positive or negative value p in order for the CLs in
*                        the aero. database to be scaled in a way that changes
*                        L/D by p%.  E.g., -5.0 lowers L/D (and CL) by 5%.
*                        Enter 0. to suppress any scaling of lift coefficients.
*
*        LINCON     Control for specifying linear constraint TYPES (NPOPT only);
*                      0 means no linear constraints (NCLIN = 0);
*                      1 means pitch rate is limited for all alpha control pnts;
*                        NCLIN is set to NDV_ALPHA - 1;
*                      2 means roll rate is limited for all bank control points;
*                        NCLIN is set to NDV_BANK - 1;
*                      3 means both pitch and roll rate are constrained;
*                        NCLIN is set to NDV_ALPHA + NDV_BANK - 2
*                   [Analogous treatment of Beta is not worth the trouble.]
*
*        NCNLN      Number of nonlinear constraints (NPOPT only)
*        NITMAX     Number of optimization steps or major iterations;
*                   NITMAX = 0 means one analysis only;
*                   NITMAX = 1 means initial analysis + initial gradient as
*                            probably desired for a study of gradients
*        KNOT1_A    First alpha knot to optimize (if NDV_ALPHA > 0)
*        KNOT1_B    First bank   "   "   "   "   (if NDV_BANK  > 0)
*                   [The first Beta knot is assumed to be 1 if MODE_VAR = 5;
*                   likewise for the first CL, CD, CY knots if MODE_VAR = 6.]
*
*        NNOPT      Optimizer switch:
*                      1 = QNMDIF2 (unconstrained [disabled now]);
*                      2 = NPOPT (constrained)
*        NNSPLINE   Controls splining of control points for alpha and bank;
*                   Traj interpolates the angle control points it finds in
*                   alpha.dat and/or bank.dat at every time step;
*                   NNSPLINE > 0 in traj_opt.inp should match use of
*                   "Schedul Fnct_alpha_mono" etc. in the Traj control file,
*                   in order for similar interpolations to be indicated by
*                   Traj_opt in traj_opt.qplot:
*                      1 = piecewise linear interpolation;
*                      2 = monotonic cubic spline;
*                      3 = non-monotonic ("Bessel") cubic spline
*        NNGRAD     Controls gradient scheme:
*                      0 = conventional finite differencing;
*                      1 = adjoint method;
*                      2 = finite differencing, but CONFUN does OBJFUN F & g;
*                      3 = adjoint method, with CONFUN used as for NNGRAD = 2
*        NNDIFF     Controls finite differencing for objective & constraints:
*                      2 = 2-point forward differencing with given h;
*                      3 = 3-point central differencing with h * EPSOBJ**(-1/6)
*        NNTIME     Controls time step to use from trajectory results:
*                      0 = last time step;
*                      N > 0 is the Traj journal block time step to use;
*                      N will be bounded by the actual last step
*        NALLOW     Upper limit on number of trajectory calculations
*        MAX_STEP   This should be at least as large as Step_mx in the Traj
*                   control file if APC or SURF* constraints are present;
*                   it is a limit on the number of uniform steps in the time
*                   history maintained by Traj, some items of which are copied
*                   in trajectory.c for use by Traj_opt.
*
*                   OPTIMIZER INPUTS:
*
*        OBJBND     Lower bound on objective function
*        ETA        Controls line search's acceptance of a sufficiently lower
*                   objective. 0 < ETA < 1.0; try 0.2 for expensive objectives
*        OPTOL      Minimum stepsize of QNMDIF2 line search;
*                   see reference to optional NPOPT/SNOPT inputs below
*        STEPMX     Maximum stepsize of line search; units are obscure; use a
*                   smaller value to help avoid a big initial step (but too
*                   small can easily make NPOPT think it can't move)
*        EPSOBJ     Minimum absolute value of a significant difference in OBJ
*        ZETA       Reduces perturbation in CENDIF2 if a function eval. fails
*        ENTRY_DV   Velocity increment added to entry velocity found in traj.in
*                   and also to velocities calculated from ascent data if a
*                   TABORT design variable is present
*
*        H_ALPHA    2-pt. differencing interval for alpha variables (degrees);
*        (H_CL)     adjusted larger if central differencing specified (NNDIFF);
*                   also used for coefficient perturbations, forward & back-
*                   ward) when estimating aero. sensitivities (MODE_VAR = 0);
*        BL_ALPHA,  Bounds on alpha control points (degrees); note that
*        BU_ALPHA   alpha and bank variables are scaled for NPOPT purposes;
*        (BL_CL,    if MODE_VAR = 6, enter H_CL, BL_CL, BU_CL on this line
*         BU_CL)
*        PITCH_RATE Finite rate used when expanding alpha steps to spline form
*
*        H_BANK     Differencing interval for bank variables (degrees)
*        (H_CD      if MODE_VAR = 6)
*        BL_BANK,   Bounds on bank angle control points (degrees)
*        BU_BANK
*        (BL_CD, BU_CD if MODE_VAR = 6)
*        ROLL_RATE  Finite rate used when expanding bank steps to spline form
*
*        H_BETA     Differencing interval for beta variables (degrees)
*        (H_CY      if MODE_VAR = 6)
*        BL_BETA,   Bounds on sideslip angle controls (degrees)
*        BU_BETA
*        (BL_CY, BU_CY if MODE_VAR = 6)
*        YAW_RATE   <Inactive>
*
*                   OBJECTIVE FUNCTION INPUTS:
*
*        RHO_xx is a multiplier for the corresponding contribution of xx
*        to the (possibly composite) objective function being minimized.
*        Multipliers are also used to scale the objective function to be
*        O(1) at the minimum.
*
*        Be sure to have at least one RHO* input nonzero.
*
*        RHO_CROSS_RANGE   [latitude is optimized for now, for equatorial orbit]
*                          N.B.:  Bank > 0 rolls the vehicle into the southern
*                          hemisphere from the equator.  Use either bank < 0 or
*                          RHO_CROSS_RANGE = -1. (but not both) to maximize the
*                          northern hemisphere latitude.  RHO * (-CROSS_RANGE)
*                          is always the quantity minimized.
*        RHO_DOWN_RANGE    [-(surface arc) is the quantity minimized]
*        RHO_HEAT_LOAD     Stagnation point heat load (joules/cm^2)
*                          [-(heat load) is minimized; see also RHO_HEATLD]
*        RHO_MAX_ACCEL     Max. acceleration magnitude in Gs
*        RHO_MAX_W_TEMP    Max. temperature            deg K  (stagnation point)
*        RHO_MAX_DYN_PR    Max. dynamic pressure       pascals
*        RHO_MAX_H_FLUX    Max. total heat flux        watts/cm^2  (stagn. pt.)
*
*        RHO_TABORT        If a TABORT variable is present after variable #
*                          NDV_ANGLES, RHO_TABORT > 0 minimizes the ascent abort
*                          time for which a target landing site can be reached;
*                          RHO_TABORT < 0 maximizes the abort time;
*                          use in conjunction with ENDLAT & ENDLON or ENDIST
*                          equality constraints, or with RHO_TARGET_POINT > 0
*
*        RHO_END_ALTITUDE  Use this multiplier in conjunction with terminating
*                          trajectories at (say) Mach 2, in order to maximize
*                          altitude over a target landing site (see TARGET_*).
*        RHO_ECCENTRICITY  Minimizing the orbit eccentricity means minimizing
*                          the velocity increment needed to circularize an
*                          orbit (later).  Use this for aero-capture problems
*                          (run to apoapsis), with constraints on apoapsis and
*                          orbital inclination.
*        RHO_HEATLD        Integrated heat load using an aerothermal database
*                          for multiple surface points with area-based weighting
*        RHO_TRIM          Time-integrated squared deviation of pitching moment
*                          (or flap deflection) from TARGET_TRIM.
*
*        RHO_TARGET_POINT  This multiple of the (squared) distance from target
*                          end point (TARGET_LATITUDE, TARGET_LONGITUDE) is
*                          minimized; the units are degrees ** 2;
*                          the appropriate ending altitude should appear in
*                          the traj.in control file via the Floor keyword
*        TARGET_LATITUDE   Target end-of-trajectory pt. (degrees), ...
*        TARGET_LONGITUDE  ... if RHO_TARGET_POINT > 0
*        TARGET_TRIM       Desirable bound on |Cm| or |flap deflection|.  Use
*                          the value 0. if RHO_TRIM > 0; use some small positive
*                          value if the corresponding constraint is being used
*                          instead (because zero everywhere is unlikely to be
*                          feasible).
*
*        RHO_ALPHA_TVD     Total-variation-diminishing contribution to the
*        (RHO_CL_TVD)      objective; adds a (small) multiple of the sum of
*                          (Alpha(i+1) - Alpha(i))**2 / TK_ALPHA(NK_ALPHA),
*                          where i is a control pt., not a time history value;
*                          0.01 is a typical input for RHO_ALPHA_TVD, which
*                          serves for smoothing of CL vs. T if MODE_VAR = 6;
*                          try 1.0 in this case
*        RHO_BANK_TVD      Total-variation-diminishing contribution to the
*        (RHO_CD_TVD)      objective, as for RHO_ALPHA_TVD but for the bank
*                          control pts.
*                          Rule of thumb:
*                          Use RHO_*_TVD = 1.E-m where m is the decimal place
*                          of the total objective function intended to be
*                          affected by this contribution.
*        RHO_BETA_TVD      Analogous contribution for sideslip if MODE_VAR = 5
*        (RHO_CY_TVD)      (or for CY if MODE_VAR = 6)
*        RHO_DURATION      Duration (total time) of the trajectory - probably
*                          in conjunction with reasonable constraints on the
*                          terminal velocity and flight path angle
*        RHO_LST_SQRS      Minimize the sum of squared differences between
*                          the current trajectory and a target trajectory,
*                          in (latitude, longitude, altitude) space; enter
*                          a negative multiple to indicate use of time
*                          ranges normalized to [0, 1]; otherwise, the current
*                          trajectory is compared only at times within the
*                          range of the target trajectory times; the target
*                          trajectory should look like this:
*
*                             ! Time, s  Lat, deg  Long, deg  Alt, km
*                                 0.500    12.345    -23.567   20.000
*                                 1.000    12.345    -23.578   20.000
*                                 1.500    12.346    -23.591   19.998
*                                  :         :          :        :
*
*
*                   SUPPLEMENTARY OPTIMIZATION VARIABLE INPUTS
*
*                   If NDV = NDV_ANGLES = NDV_ALPHA + NDV_BANK [+ NDV_BETA],
*                   meaning no additional variables, just include two header
*                   lines in traj_opt.inp.  Likewise for MODE_VAR = 3 or 4.
*
*        #          Ordinal number of design variable added to NDV_ANGLES
*
*        VTYPE      6-character design variable type or name:
*
*                   TABORT   Time at which to initiate an ascent abort;
*                            if present, a launch trajectory dataset is
*                            read from "traj_opt.ascent" in this format,
*                            using INERTIAL coordinates:
*
*                               Title of ascent trajectory dataset
*                               N_ASCENT (# pts. defining the ascent trajectory)
*                               T (sec) V (kps) GAMMA HEADING ALT (km) LAT LONG
*                               T       V       GAMMA HEADING ALT      LAT LONG
*                               T       V       GAMMA HEADING ALT      LAT LONG
*                               :       :       :     :       :        :   :
*
*                   [Allow entry conditions as variables here some day?]
*
*        V          Optimization variable value as seen by the optimizer;
*                   see VSCALE
*        VSCALE     Scale factor for the optimization variable to keep it ~O(1);
*                   V * VSCALE is the actual quantity used by Traj_opt & Traj
*        AITCH      Step size used for estimating the gradient of the objective
*                   w.r.t. the variable by forward differencing (NNGRAD = 0, 2);
*                   the actual perturbation used is AITCH * VSCALE;
*                   AITCH is also used for forward derivatives of any nonlinear
*                   constraints; see NNDIFF above for central differencing
*        BL         Lower bound on the design variable as seen by the optimizer;
*                   BL = -999. means the variable has no lower bound, else
*                   BL scaling should match that of the variable - e.g., for a
*                   variable with VSCALE = 0.1, BL = -10. means the effective
*                   values are >= -10.*0.1 = -1., consistent with AITCH usage
*        BU         Upper bound on the design variable as seen by the optimizer;
*                   BU = 999. means the variable has no upper bound, else
*                   BU scaling should match that of the variable
*
*                   LINEAR CONSTRAINTS:
*
*        If LINCON = 0, just include two header lines.
*        If LINCON > 0, further linear constraint description lines should be
*        entered as indicated by the LINCON description above.
*        Now that use of perturbing shape functions for the initial alpha and
*        bank schedules has been eliminated, the only linear constraints are
*        on the alpha knots (to limit pitch rate) and/or the bank knots (to
*        limit roll rate).  Rather than requiring a constraint input for each
*        pair of adjacent control points, a single input line here activates
*        each TYPE of linear constraint.
*        If NDV_ALPHA = 0 but LINCON = 1 or 3, the PITCH inputs here will be
*        ignored; likewise for the NDV_BANK = 0 case.
*
*        #          Ordinal number of linear constraint type - not used
*        LCTYPE     6-character linear constraint type:
*
*                   PITCH   Pitch rate applied to alpha control point pairs;
*                           both BL and BU are used (see below)
*                   ROLL    Roll rate applied to bank control point pairs;
*                           BL is not used; see BU below
*                  [Yaw rate constraints are not worth the bother.]
*
*        BL         Lower bound for linear constraint [type]; -999. = -BIGBND;
*                   if LCTYPE = PITCH, BL and BU both apply (degrees per second)
*        BU         Upper bound for linear constraint [type];  999. =  BIGBND;
*                   if LCTYPE = ROLL, BU is the upper bound on the roll rate
*                   magnitude in degrees per second for all adjacent pairs of
*                   bank control points; BL is ignored
*        TLCON      Time at which linear constraint applies (seconds);
*                   if LCTYPE = PITCH or ROLL, TLCON is ignored
*        ILCON      Integer control for the linear constraint;
*                   if LCTYPE = PITCH or ROLL, ILCON is ignored
*
*                   NONLINEAR CONSTRAINTS:
*
*        If NCNLN = 0, just include two header lines.
*
*        #          Ordinal number of nonlinear constraint - not used
*        NLCTYPE    6-character nonlinear constraint type of name:
*
*                   ACCEL   Peak acceleration magnitude; use XNLCON to specify
*                           the desired upper bound (Gs);
*                           BU should be zero; BL is either zero (for exact
*                           peak G) or -XNLCON value if lower is OK, e.g., -3.
*                   CRANGE  Cross range <degrees latitude for now>
*                   DRANGE  Down  range (surface arc, km)
*                   DYN_PR  Dynamic pressure; use XNLCON to specify the desired
*                           upper bound (pascals, e.g., 23940. for 500 psf);
*                           BU should be zero; BL is either zero (for exact
*                           peak Q) or -XNLCON value if lower is OK
*                   ENDALP  Ending alpha (deg)
*                   ENDALT  Ending altitude (km)
*                   ENDBNK  Ending bank angle (deg)
*                   ENDFPA  Ending flight path angle (gamma, degrees)
*                   ENDHED  Ending heading angle (psi, degrees; 0 = N, 90 = E)
*                   ENDLAT  Ending latitude (deg) ! ENDLAT & ENDLON seem to do
*                   ENDLON  Ending longitude  "   ! better than ENDIST
*                   ENDVEL  Ending velocity (km/sec, relative)
*                   ENDIST  Ending distance from (TARGET_LATITUDE, -_LONGITUDE)
*                           degrees (not deg ** 2); use zero for both bounds
*                   HEATLD  Integrated heat load over multiple surface points
*                   H_LOAD  Stagnation point heat load (joules/cm^2)
*                   H_FLUX  Peak total heat flux at stagnation pt. (watts/cm^2)
*                   QALPHA  Dynamic pressure * alpha;
*                           use XNLCON (psf * degrees) for the upper bound;
*                           15,000 is a likely value; set BU = 0. & BL = -XNLCON
*                   STG_PR  Peak stagnation pressure (pascals)
*                   W_TEMP  Peak wall temperature at stagnation pt. (deg K)
*
*                   DESCNT  Descent rate at end of trajectory, probably for
*                           ascent abort trajectories (km/sec)
*                   PK_ALT  Peak altitude (km); use XNLCON for the upper bound
*                           (km); if this is lower than the entry altitude,
*                           start measuring from where it first drops below the
*                           bound; set BU = 0. and BL = -XNLCON
*
*                   APOAPS  Target apoapsis for aero-capture (km from planet CG)
*                   ECCENT  Ending orbital eccentricity, probably for aero-
*                           capture
*                   INCLIN  Ending orbital inclination (degrees), probably for
*                           aero-capture:  0 = equatorial due East; 90 = polar
*
*                   APC     Aerothermal performance constraint:
*                           file "apc.dat" should contain points as follows:
*
*                              Title
*                              n      ! # points
*                              V1  H1 ! High to low, km/sec and km;
*                              V2  H2 ! Traj_opt reverses the order if necessary
*                              ::  ::
*                              Vn  Hn
*
*                           The constraint tends to force the minimum of
*                           "current altitude minus interpolated APC altitude"
*                           to be >= BL = 0;
*                           set BU = max. likely altitude difference
*
*                   UPPERC  Upper corridor constraint analogous to the APC,
*                           forming an upper boundary in (velocity, altitude)
*                           space (probably corresponding to Q = 10 psf).
*                           File "upper_corridor.dat" should match the APC file
*                           format.  Program CONSTANT_Q is available for
*                           constructing the initial curve, which may need
*                           padding via spline interpolation.  The suggested
*                           choice is 478.8 pascals (10 psf).  The abscissas
*                           should be redistributed ~ uniformly to 50 - 100
*                           points in the range ~ 0 - 12 km/sec;
*                           set BU = zero and BL ~ -50.
*
*                   OSCILL  An upper bound of 0. means any increase in altitude
*                           during reentry is penalized; oscillations may thus
*                           be damped; severe skipping may also be overcome.
*                           (All altitude increases at each step in the journal
*                           history are summed; try scale = 1 & BL = -10.)
*
*                   SURF1   Distributed heating constraint type 1 (temperature):
*                           the total violation of the temperature limits over
*                           all surface points is constrained to the upper bound
*                           of 0.  Use BL ~ -1000. (largest likely undershoot in
*                           degrees K).  See further description above.
*                   SURF1R  Variant of SURF1 using temperatures derived from
*                           interpolated heat RATEs.  This requires surface
*                           emissivities in place of the upper bounds on the
*                           heat rates in 'traj_opt.surface_bounds'.  See
*                           the NOTE above where file formats are described.
*                   SURF2   Distributed heating constraint type 2 (heat flux):
*                           the total violation of heat rate limits over all
*                           surface points is constrained to BU = 0.
*
*                   TRIM    The time integral of the square of the instantaneous
*                           pitching moment (or possibly flap deflection) dev-
*                           iation from TARGET_TRIM.  Choose BL = 0., BU > 0.
*                           according to the estimated value of the integral
*                           at convergence for the indicated TARGET_TRIM.
*                           Use of a RHO_TRIM contribution to the objective is
*                           probably better than using this constraint.
*                           
*
*        BL         Lower bound on nonlinear constraint value; -999. = -BIGBND
*        BU         Upper bound on nonlinear constraint value;  999. =  BIGBND
*                   N.B.:  these bounds should be in the same units as the
*                   quantities being bounded, which is sometimes an "area" or
*                   "distance", not the intended constraint quantity, for which
*                   (as in the case of ACCEL) XNLCON is used to enter the upper
*                   limit; the bounds and XNLCON (if used) are entered in real
*                   space units; SNLCON will be applied by Traj_opt to the given
*                   BL and BU; this is in contrast to the bounds on the optimiz-
*                   ation variables themselves, which (as for their finite diff-
*                   erencing intervals) refer to the scaled variables
*        XNLCON     Real quantity needed (or not) by the nonlinear constraint;
*                   used for entering the simple bound imposed on a quantity
*                   such as acceleration for which all violations are integrated
*                   w.r.t. time and forced to the lower bound of zero
*        INLCON     First  index needed (or not) by the nonlinear constraint;
*                   [Unused so far.]
*        JNLCON     Second index needed (or not) by the nonlinear constraint;
*                   [Unused so far.]
*        SNLCON     Scale factor used to provide the optimizer with nonlinear
*                   constraint derivatives comparable to those of the objective
*                   and of other constraints - preferably O(1)
*
*                   OPTIONAL INPUTS FOR NPOPT/SNOPT:
*
*        Traj_opt generates a run-time "specs" file for NPOPT.  It also
*        provides a namelist, $NPOPTIONS, for entering further options via
*        traj_opt.inp.  Remember to start each line in column 2.  Sample:
*
*         $NPOPTIONS
*         LEVELVER=-1, MAJORPL=1, MINORPL=0, MINORIL=2000, NPRFREQ=100,
*         TOLLIN=1.E-6, TOLNLIN=1.E-2, TOLOPT=1.E-2,
*         PENPARAM=0., STEPLIM=0.01, TIGHTEN=0.1,
*         $END
*
*        Note:  If NNDIFF /= 3 and 0. < TIGHTEN < 1., Traj_opt will reset
*               NNDIFF to 3, scale TOLNLIN, TOLOPT, and STEPLIM by TIGHTEN
*               upon return from NPOPT, and call it again for a warm start.
*               Central differencing offers some hope of satisfying the tighter
*               tolerances where forward differencing normally does not.
*               Starting with looser tolerances and 2-point differencing allows
*               switching to more accurate gradients once the problem has been
*               largely solved, and may reduce excessive numbers of major
*               iterations.  However, there is no easy way of terminating when
*               the solution is clearly "good enough" for practical purposes.
*               TIGHTEN = 0. (the default) suppresses the warm restart option.
*
*        See the NPSOL or SNOPT User Guides for descriptions of the other
*        namelist parameters.
*
*        Also:  File npopt.specs may be used to override further defaults
*               not addressed by the namelist.  This file is optional
*               (need not be present).
*
*
*     HISTORY:
*
*        02/11/00  DAS  Initial framework adapted from SYN87-SB.
*        02/28/00   "   Faked alpha/beta data; minimized "arc" length (in time).
*        04/04/00   "   Installed further quantities compatible with the new
*                       trajectory.c interface to Traj.
*        05/02/00   "   Implemented nonlinear constraints on peak quantities.
*        05/08/00   "   Trajectory.c interface now determines the numbers of
*                       bank and alpha maneuvers during initialization.
*        05/09/00   "   Implemented linear constraints on alpha and bank angles.
*        05/15/00   "   NNGRAD = 2 or 3 allows CONFUN to provide OBJFUN results.
*        05/18/00   "   Arranged to pass Traj_opt command-line arguments into
*                       the Traj package (possibly SGI/IRIX-specific).
*        05/24/00   "   Added NNTIME option, and MAX_TIMES.
*        06/01/00   "   NNVAR option allows use of spline knots as variables,
*                       but kludging is required until Traj can use them too.
*        06/12/00   "   Incorporated an APC constraint option.
*        06/19/00   "   Variable scaling for spline knots case: leave STEPMX,
*                       H_ANGLES, and BU_ALPHA in degrees.
*        06/20/00   "   Introduced KNOT1_A and -B.
*        07/11/00   "   NNDIFF = 3 allows central differencing.
*        07/21/00   "   Traj now interpolates alpha and bank at every time step,
*                       so we no longer interract with its maneuvering script.
*        08/25/00   "   Introduced target (latitude, longitude) option, and
*                       maximum ending-altitude option.
*        08/30/00   "   Eliminated sine bump-type perturbing variables;
*                       introduced TABORT optimization variable and nonlinear
*                       constraints for ending latitude and longitude.
*        09/07/00   "   Added ENDIST constraint option.
*        04/19/01   "   Accessed Traj's orbital parameters for possible use
*                       in aero-capture calculations; inclination may be
*                       constrained.
*        04/24/01   "   Orbital eccentricity may be minimized for aero-capture;
*                       Roman Jits pointed out the need for a target apoapsis.
*        04/27/01   "   Bounds on the alpha and bank variables are now inputs.
*        05/08/01   "   ECCENT constraint allows aero-capture calculations to
*                       seek a low but not necessarily minimum eccentricity.
*        06/02/01   "   Added Step Limit to optional NPOPT inputs.
*        07/17/01   "   Added a PK_ALT constraint, but it may not stop maximum
*                       down-range trajectories from skipping.
*        07/20/01   "   Started implementing distributed heating constraints.
*        08/02/01   "   Implemented the UPPERC constraint for max. down-range.
*        08/22/01   "   Finished implementing distributed heating constraints.
*        08/29/01   "   Dump the surface point histories at the end of the run.
*        09/07/01   "   Added descent-rate constraint (DESCNT) because some
*                       abort trajectories were plummeting at the end.
*        09/10/01 DS/JR Distributed surface constraints now integrate areas
*                       above the bounds for each point, rather than summing
*                       the peak violations.
*        09/19/01  DAS  Upgrade of Traj, which now captures end-point conditions
*                       precisely.  Traj's number of maximum values grew from
*                       10 to 11, but the 11th is currently ignored.
*        09/21/01 DS/JR Implemented pitch & roll rate constraint options.  James
*                       suggested the linear constraint approach as opposed
*                       to the nonlinear integrated-violations alternative.
*        10/14/01  DAS  Switched from (dense) NPOPT to SNOPT-based NPOPT, only
*                       in case the latest niceties help somewhat.
*        10/18/01   "   An optional "npopt.specs" file avoids hard-coding new
*                       optional control parameters used by NPOPT/SNOPT.
*        11/08/01   "   Larger numbers of surface heating points forced
*                       separate files for each data type ('.temperatures',
*                       '.fluxes', ...).  Summing the temperature undershoots
*                       for a constraint value was bad - pick the smallest
*                       undershoot instead.
*        12/09/01   "   Introduced CURVE_VIOLATIONS in anticipation of using
*                       time histories rather than just peaks for G load and
*                       other possible constraints; used it in place of the
*                       original APC and upper corridor implementations, and
*                       for the DYN_PR and (new) QALPHA constraints.
*        03/20/02   "   ACCEL constraint now uses the time history, which
*                       had to be added to what trajectory.c returns.
*        04/08/02   "   RHO_TRIM or TRIM constraint, along with TARGET_TRIM,
*                       provide for minimizing the flap deflections over time
*                       (or possibly the pitching moment coefficient).  Used
*                       for a Crew Escape Module application, but beware:
*                       we don't want to push CM towards second zeros of the
*                       CM vs. Alpha curves for given Mach numbers.
*        04/23/02   "   Added RHO_ALPHA_TVD and RHO_BANK_TVD after observing
*                       that oscillations in both Alpha and bank can be smoothed
*                       (by hand) with negligible effect on abort trajectories,
*                       an indication of flat objective functions with poorly-
*                       defined minima.  These might help define the optimum
*                       better.
*        04/25/02   "   Switched from integral of |Alpha(i+1) - Alpha(i)| w.r.t.
*                       time to (Sum dAlpha**2) / time(n_alpha) for TVD idea.
*        05/02/02   "   Check for bad aero. database dimensions (count lines).
*                       Traj set-up really should return these dimensions.
*        05/07/02   "   PERCENT_LOVERD input allows for scaling drag and hence
*                       estimate sensitivities to L/D; termination constraints
*                       ENDVEL and ENDFPA have been added belatedly.
*        06/12/02   "   Renamed "CONSTQ" constraint as "UPPERC[orridor]", etc.
*        06/19/02   "   Avoided integrating w.r.t. index number for APC and
*                       UPPERC by passing time steps as well as vel. & alt.;
*                       suppressed trajectory calcs. for gradient elements
*                       corresponding to Alpha & bank controls "off the end";
*                       printed objective gradient in parallel with Jacobian.
*        06/21/02   "   Added option to minimize the trajectory duration.
*        06/25/02   "   Added ENDALP, ENDBNK, and ENDHED constraints.
*        07/11/02   "   Added step-function mode (MODE_VAR = 3, 4), requiring
*                       set-up and reading of "traj.alp" and "traj.bnk" from
*                       "traj_opt.steps".
*        07/23/02   "   Arranged for switching to central differencing via a
*                       warm restart.  See discussion above.
*        07/26/02   "   Removed QNMDIF2 option.
*        08/16/02   "   Warm restarts can be problematic.  Introduced the
*                       TIGHTEN namelist parameter for more control, and
*                       tightened STEPLIM as well when it is used.
*                       CPU time for gradients shouldn't be printed from SOLVE.
*        09/10/02   "   SURF1R constraint introduced to reduce trilinear
*                       interpolation errors in temperature by deriving
*                       temperature from interpolated heat flux.
*        01/22/03   "   Exercised aero. database sensitivity mode and made the
*                       outputs more readable.  Left CENDIF2's inappropriate
*                       printout alone for this unusual application.
*        02/07/03   "   Safeguarded the case of a trajectory not overlapping
*                       an APC at all: no good action, so stop with a message
*                       to remove the APC.  See CURVE_VIOLATIONS.
*        03/05/03   "   Option to match a target trajectory forced changes in
*                       trajectory.c as well as introducing RHO_LST_SQRS.
*        05/15/03   "   Added sideslip capability (MODE_VAR = 5) and option
*                       to optimize CL, CD, CY w.r.t. time (MODE_VAR = 6).
*                       Cleaned out unconstrained option (QNMDIF2), but CENDIF2
*                       is still employed by MODE_VAR = 0 case.
*        05/16/03   "   CL, CD, CY need smoothing options, so RHO_BETA_TVD has
*                       been added; the same 3 inputs work for MODE_VAR = 5 & 6
*        04/14/07   "   A capsule case required 180 - alpha data, making alpha
*                       descending in the aerodynamic and aerothermal tables.
*                       Only the aerothermal table interpolation needed a fix.
*        10/14/10   "   The LCTYPE array was too short for MODE_VAR = 3 or 4.
*        10/20/10   "   PERCENT_LOVERD is now used to scale CL instead of CD.
*
*     SPONSOR:
*
*        Aerothermodynamics Branch (TSA), especially James Reuther
*        NASA Ames Research Center
*        Moffett Field, CA
*
*     SOFTWARE AUTHORS:
*
*        David Saunders,  ELORET Corporation/ERC, Inc., NASA Ames (Traj_opt)
*        Gary Allen, Jr., ELORET Corporation/ERC, Inc., NASA Ames (Traj)
*
********************************************************************************

*     Global variables:

      USE CONSTS
      USE OPT

      IMPLICIT NONE

*     Local constants:

      REAL, PARAMETER ::
     >   BIGBND = 1.E+10, ONE = 1., ZERO = 0.

      CHARACTER, PARAMETER ::
     >   PROGNAME*8 = 'Traj_opt'

*     Local variables:

      INTEGER
     >   I, IER, ILINCON, INFORM, IOS, IRETURN, ISTART,
     >   J, K, L, LEVELVER, M, MAJITS, MAJORPL, MINORPL, MINORIL,
     >   N, NCALL, NFUNCT, NITER, NPRFREQ

      REAL
     >   H_ALPHA, BL_ALPHA, BU_ALPHA,  H_CL, BL_CL, BU_CL,
     >   H_BANK,  BL_BANK,  BU_BANK,   H_CD, BL_CD, BU_CD,
     >   H_BETA,  BL_BETA,  BU_BETA,   H_CY, BL_CY, BU_CY,
     >   BOUNDLO, BOUNDUP,  CPU1, CPU2, DT, ETA, OBJPAR, OPTOL,
     >   PENPARAM, STEPLIM, STEPMX, TIGHTEN, TLINCON, TOLLIN, TOLNLIN,
     >   TOLOPT, ZETA

*     The following equivalences avoid additional rarely-used inputs:

      EQUIVALENCE
     >  (H_CL, H_ALPHA), (BL_CL, BL_ALPHA), (BU_CL, BU_ALPHA),
     >  (H_CD, H_BANK),  (BL_CD, BL_BANK),  (BU_CD, BU_BANK),
     >  (H_CY, H_BETA),  (BL_CY, BL_BETA),  (BU_CY, BU_BETA)

      LOGICAL
     >   FAIL

      CHARACTER
     >   BUFFER*132, LINCON_TYPE*6, NULL*1

*     Namelist to override NPOPT's options file defaults:

      NAMELIST /NPOPTIONS/
     >   LEVELVER, MAJORPL, MINORPL, MINORIL, NPRFREQ,
     >   PENPARAM, STEPLIM, TIGHTEN, TOLLIN, TOLNLIN, TOLOPT

*     Procedures:

      INTEGER
     >   IARGC        ! SGI/IRIX intrinsic (command line argument count)

      EXTERNAL
     >   CENDIF2,     ! Utility for estimating gradient, diagonal of Hessian,
                      ! and (optionally) good finite differencing intervals
     >   CHECK_INPUTS,! Catch some of the possible input errors
     >   CONFUN,      ! Nonlinear constraint function in NPOPT form
     >   DYNAMIC,     ! Dynamic allocation of Fortran 90 work-space
     >   GETARG,      ! SGI/IRIX subroutine (extracts one argument per call)
     >   IARGC,       ! SGI/IRIX function (see above)
     >   NPINIT,      ! SNOPT utility for initialization
     >   NPSPEC,      ! SNOPT utility for overriding some defaults
     >   NPOPT,       ! Constrained optimization package
     >   OBJFUN,      ! Objective function routine in NPOPT form
     >   OBJFUN_AERO  ! Objective function in CENDIF2 form, for estimating
      EXTERNAL        ! sensitivities to the aerodynamic database entries
     >   PERTURB,     ! Applies shape functions to initial bank/alpha schedules
**** >   QNMDIF2,     ! Unconstrained optimization package
     >   RDVARS,      ! Reads design variables/shape functions
     >   RVERSE,      ! Reverses a vector, possibly in-place
     >   SETLCON,     ! Sets up the linear constraint matrix
     >   SOLVE,       ! Objective function routine in QNMDIF2 form, also used
                      ! by OBJFUN and CONFUN
     >   TRAJECTORY,  ! C interface to the Traj package
     >   USER         ! Application-specific routine called by QNMDIF2 to
     >                ! reevaluate the best line search function & gradient or
     >                ! to summarize the MODE_VAR = 0 results.

*     Execution:

      CALL SECOND (CPU0)

      OPEN (LREAD,    FILE='traj_opt.inp',   STATUS='OLD')
      OPEN (LWRIT,    FILE='traj_opt.out',   STATUS='NEW')
      OPEN (LRESULTS, FILE='traj_opt.qplot', STATUS='NEW')

      WRITE (LWRIT, 1)
    1 FORMAT (
     >   /, ' Traj_opt Version:     October 14, 2010',
     >   /, ' Sponsor:              Aerothermodynamics Branch,',
     >   /, '                       (TSA), NASA Ames Research Center',
     >   /, ' Traj_opt Development: David Saunders,  ELORET',
     >   /, ' Traj Package:         Gary Allen, Jr., ELORET',
     >   /, ' Optimization:         SNOPT (NPOPT interface), Stanford',
     >   /, ' Inspiration:          James Reuther, ASA Branch')

      NULL    = CHAR (0)
      EPSMCH  = EPSILON (EPSMCH)
      PI      = 4.* ATAN (ONE)
      RPI     = ONE / PI
      DEG2RAD = PI / 180.
      RAD2DEG = 180./ PI
      RSTEFAN = ONE / 5.66097E-12 ! Stefan-Boltzmann constant, watts/(cm^2 K^4)

*     Gather any command line arguments to be passed to Traj.  We cannot
*     quite set up the char *argv[] structure that a C program does, but
*     the TRAJECTORY interface routine can reformat the ARGS character array:

      NARGS = MIN (IARGC (), MAXARGS)     ! Unix intrinsic omits program name
      ARGS(0) = PROGNAME                  ! but C expects it as element 0
      L = MIN (LEN_TRIM (PROGNAME) + 1, MAXCHARS)
      ARGS(0)(L:L) = NULL

      DO I = 1, NARGS ! May be zero

         CALL GETARG (I, ARGS(I))         ! Unix intrinsic

         L = MIN (LEN_TRIM (ARGS(I)) + 1, MAXCHARS)
         ARGS(I)(L:L) = NULL
      END DO

      NARGS = NARGS + 1

*     Echo the optimization inputs to the output file:

      WRITE (LWRIT, '(/, A, /)') ' Optimization control inputs:'

      DO
         READ (LREAD, '(A)', IOSTAT=IOS) BUFFER
         IF (IOS < 0) EXIT

         I = LEN_TRIM (BUFFER)
         WRITE (LWRIT, '(1X, A)') BUFFER(1:I)
      END DO

      REWIND (LREAD)

      READ (LREAD, *) ! Major design controls
      READ (LREAD, *)
      READ (LREAD, *) NDV, NDV_ALPHA, NDV_BANK, NDV_BETA, MODE_VAR,
     >                VEHICLE

      LEN_VEHICLE = LEN_TRIM (VEHICLE)

      READ (LREAD, *)
      READ (LREAD, *) N_MACH, N_ALPHA, N_RE, N_BETA, PERCENT_LOVERD
      READ (LREAD, *)
      READ (LREAD, *) LINCON, NCNLN, NITMAX, KNOT1_A, KNOT1_B
      READ (LREAD, *)
      READ (LREAD, *) NNOPT, NNSPLINE, NNGRAD, NNDIFF, NNTIME, NALLOW,
     >                MAX_STEP

      READ (LREAD, *) ! Optimizer inputs
      READ (LREAD, *)
      READ (LREAD, *) OBJBND, ETA, OPTOL, STEPMX, EPSOBJ, ZETA, ENTRY_DV
      READ (LREAD, *)
      READ (LREAD, *) H_ALPHA, BL_ALPHA, BU_ALPHA, PITCH_RATE
      READ (LREAD, *)
      READ (LREAD, *) H_BANK,  BL_BANK,  BU_BANK,  ROLL_RATE
      READ (LREAD, *)
      READ (LREAD, *) H_BETA,  BL_BETA,  BU_BETA,  YAW_RATE
      READ (LREAD, *)

      READ (LREAD, *) RHO_CROSS_RANGE, RHO_DOWN_RANGE, RHO_HEAT_LOAD,
     >                RHO_MAX_ACCEL,   RHO_MAX_W_TEMP, RHO_MAX_DYN_PR,
     >                RHO_MAX_H_FLUX
      READ (LREAD, *)
      READ (LREAD, *) RHO_TABORT, RHO_END_ALTITUDE, RHO_ECCENTRICITY,
     >                RHO_HEATLD, RHO_TRIM
      READ (LREAD, *)
      READ (LREAD, *) RHO_TARGET_POINT, TARGET_LATITUDE,
     >                TARGET_LONGITUDE, TARGET_TRIM
      READ (LREAD, *)
      READ (LREAD, *) RHO_ALPHA_TVD, RHO_BANK_TVD, RHO_BETA_TVD,
     >                RHO_DURATION, RHO_LST_SQRS


*     Catch some of the possible bad inputs and do some initialization:

      CALL CHECK_INPUTS (IER)  ! NK_ALPHA, NK_BANK, etc. are read here

      IF (IER /= 0) GO TO 999


*     Set up the linear constraints proper:

      L = MAX (NDV_ALPHA - 1, 0)
      M = MAX (NDV_BANK  - 1, 0)

      IF (MODE_VAR == 3 .OR. MODE_VAR == 4) THEN  ! Effective knot counts have
         L = NK_ALPHA - 1                         ! been doubled
         M = NK_BANK  - 1
      END IF

      SELECT CASE (LINCON)
         CASE (0)
            NCLIN = 0
            NLINCON_RD = 0
         CASE (1)
            NCLIN = L
            NLINCON_RD = 1
         CASE (2)
            NCLIN = M
            NLINCON_RD = 1
         CASE (3)
            NCLIN = L + M
            NLINCON_RD = 2
      END SELECT

      NBNDS = NDV + NCLIN + NCNLN      ! For variable + constraint bounds
      NROWA = MAX (NCLIN, 1)           ! Max. # linear constraints ...
      NROWJ = MAX (NCNLN, 1)           ! ... & nonlinear constraints; dims. > 0
      NROWR = MAX (NDV, 1)             ! Row dim. of Cholesky factor R
                                       ! of the Hessian of the Lagrangian
      LENIW = 3*NDV + NROWA + 2*NROWJ  ! Length of integer workspace
      LENW  = NDV*(2*NDV + NROWA +     ! Length of real workspace
     >        2*NROWJ + 20) + 11*NROWA +
     >        21*NROWJ

      LENIW = LENIW + 100000  ! Above are for dense NPOPT; sparse techniques
      LENW  = LENW  * 5       ! can't predict precisely, so pad these for SNOPT
      LENW  = MAX (LENW, 900) ! Minimum for SNOPT is more than 600
      LENCW = 900             ! 500 = min. for SNOPT; used by SNSET for restart

      IF (MODE_VAR == 0) THEN ! Special case: just aero. coef. sensitivities
         LENCW = 1
         LENIW = 1
         LENW  = 1
      END IF

      CONFUN_DOES_OBJ = NNGRAD >= 2 .AND. NCNLN > 0
      NNGRAD = MOD (NNGRAD, 2)


*     ****** Initialize the trajectory package ******

*     First, allocate space for copies of its outputs.

      ALLOCATE (ALPHA_K (NK_ALPHA),  BANK_K(NK_BANK),  BETA_K(NK_BETA),
     >          TK_ALPHA(NK_ALPHA), TK_BANK(NK_BANK), TK_BETA(NK_BETA),
     >          MAX_VALUES(0:N_MAX-1),    MAX_VALUES_INI(0:N_MAX-1),
     >          MAX_TIMES (0:N_MAX-1),    MAX_TIMES_INI (0:N_MAX-1),
     >          ENTRY_STATE(0:N_ENTRY-1), FINAL_STATE(0:N_ENTRY+1), !Alpha+bank
     >          V(NDV),       VSCALE(NDV),
     >          MACH(N_MACH), ALPHA(N_ALPHA), BETA(N_BETA), RE(N_RE),
     >          CL(N_COEFS), CD(N_COEFS), TRAJ_STEPS(1,1))

      IF (MODE_VAR == 3 .OR. MODE_VAR == 4) THEN ! Step function mode

*        Transcribe compact form of "traj_opt.steps" to the equivalent
*        "traj.alp" and "traj.bnk" expected by Traj initialization:

         ALLOCATE (ALPHA_STEPS(NSTEPS), BANK_STEPS(NSTEPS),
     >             DURATION(NSTEPS))

         CALL ANGLE_STEPS (0, MODE_VAR, NSTEPS, ALPHA_STEPS, BANK_STEPS,
     >                     DURATION, PITCH_RATE, ROLL_RATE,
     >                     ALPHA_K, BANK_K, TK_ALPHA, TK_BANK,
     >                     NDV, V, VSCALE, LDAT, LWRIT)

      ELSE IF (MODE_VAR == 6) THEN ! CL, CD, CY are functions of time

         ALLOCATE (CL_K (NK_CL), CD_K (NK_CD), CY_K (NK_CY),
     >             TK_CL(NK_CL), TK_CD(NK_CD), TK_CY(NK_CY))

      END IF

      NJOURNAL = 0 ! Suppress copying parts of Traj time history unless it is
                   ! needed for APC or for other time-history-based constraints.
                   ! (We haven't read the constraints yet.)

      IF (MODE_VAR <= 4) THEN ! All zero-side-force cases

         CALL TRAJECTORY (0, NARGS, ARGS, MODE_VAR,
     >                    NNTIME,   N_ENTRY,    ENTRY_STATE,
     >                    N_MAX,    MAX_VALUES, MAX_TIMES,
     >                    NK_BANK,  TK_BANK,    BANK_K,
     >                    NK_ALPHA, TK_ALPHA,   ALPHA_K,
     >                    N_MACH, MACH, N_ALPHA, ALPHA, N_RE, RE,
     >                    CL, CD,   PERCENT_LOVERD,
     >                    FINAL_STATE, TOTAL_TIME,
     >                    CROSS_RANGE, DOWN_RANGE, HEAT_LOAD, ORBITAL,
     >                    NJOURNAL, MAX_STEP, LAST_STEP, TRAJ_STEPS,
     >                    IRETURN)

      ELSE IF (MODE_VAR == 5) THEN ! Alpha, bank, and beta are variables

         CALL TRAJ2 (0, NARGS, ARGS, ENTRY_STATE,
     >               MAX_VALUES, MAX_TIMES,
     >               NK_BANK,  TK_BANK,  BANK_K,
     >               NK_ALPHA, TK_ALPHA, ALPHA_K,
     >               NK_BETA,  TK_BETA,  BETA_K,
     >               FINAL_STATE, TOTAL_TIME,
     >               CROSS_RANGE, DOWN_RANGE, HEAT_LOAD,
     >               MAX_STEP, LAST_STEP, TRAJ_STEPS,
     >               IRETURN)

      ELSE IF (MODE_VAR == 6) THEN ! CL, CD, and CY are variables

         CALL TRAJ3 (0, NARGS, ARGS, ENTRY_STATE,
     >               MAX_VALUES, MAX_TIMES,
     >               NK_CL,  TK_CL,  CL_K,
     >               NK_CD,  TK_CD,  CD_K,
     >               NK_CY,  TK_CY,  CY_K,
     >               FINAL_STATE, TOTAL_TIME,
     >               CROSS_RANGE, DOWN_RANGE, HEAT_LOAD,
     >               MAX_STEP, LAST_STEP, TRAJ_STEPS,
     >               IRETURN)

      END IF

      IF (IRETURN == 0) THEN
         WRITE (LWRIT, '(/, A)')
     >      ' Traj initialization failed. Check NDV, MAX_STEP, etc.'
         GO TO 999
      END IF

      ENTRY_STATE(1) = ENTRY_STATE(1) + ENTRY_DV


*     Other work-space affected by the numbers of variables & constraints:

      CALL DYNAMIC ()


*     *********************************************************************
*     If it's the special case of calculating sensitivities with respect to
*     aerodynamic coefs., do it & avoid the rest of the normal procedures:
*     *********************************************************************

      IF (MODE_VAR == 0) THEN

         OPEN (UNIT=LSAVE, STATUS='NEW',
     >         FILE='traj_opt.CL_CD_gradients')

         V(1:N_COEFS)       = CL(1:N_COEFS) ! Does it make sense to
         V(1+N_COEFS : NDV) = CD(1:N_COEFS) ! check both CL & CD?

         NCALL = -3 ! Signals the first objective function call

         CALL OBJFUN_AERO (NDV, V, OBJPAR, NCALL, FAIL)

         NFUNCT = 0

*        Estimate gradients, and optimal delta Vs if input AITCH() < 0.:

         AITCH = H_ALPHA

         CALL SECOND (CPU1)

         CALL CENDIF2 (NDV, V, OBJPAR, EPSOBJ, AITCH, GRAD, DEE,
     >                 NFUNCT, NALLOW, ZETA, LWRIT, OBJFUN_AERO, URH,
     >                 NFDEL, IERROR)

         CALL SECOND (CPU2)

         NFTOTL = 1 + NFUNCT

         WRITE (LWRIT, '(/, A, I6)')
     >      ' Total # trajectories calculated by CENDIF2: ', NFTOTL

         IF (NFTOTL >= NALLOW) THEN
            WRITE (LWRIT, 1010) 'Maximum # trajectories reached.'
         END IF

         CALL CPUTIME (CPU1, CPU2, PROGNAME,
     >      'for calculating aero database sensitivities', LWRIT)

*        Save results:

         CALL USER (3, NDV, 1, V, OBJPAR, GRAD, AITCH, ELL,
     >              DEE, OLDG, PSRCH, GTP, GNM, NFTOTL, NALLOW,
     >              NCALL, NDITER, 1, .FALSE., .TRUE., LWRIT)

         GO TO 999

      END IF


*     ***********************************************************
*     Normal optimization of a trajectory for fixed vehicle shape
*     ***********************************************************

*     Initialize the alpha & bank [& beta, or CL, CD & CY] variables,
*     including scaling and finite difference intervals:

      CALL SET_VARIABLES () ! Internal procedure at end of main program


      READ (LREAD, *) ! Supplementary design variable header lines
      READ (LREAD, *) ! ------------------------------------------

*     Read further variables NDV_ANGLES + 1, ... ?

      IF (NDV > NDV_ANGLES .AND. MODE_VAR <= 2) THEN

         CALL RDVARS (LREAD, LWRIT, LDAT) ! E.g., TABORT, in which case we
                                          ! also read ascent trajectory data
      END IF

      DO I = 1, NDV
         IF (BL(I) == -999.) BL(I) = -BIGBND
         IF (BU(I) ==  999.) BU(I) =  BIGBND
      END DO

*     Otherwise, BL & BU scaling should match that of the variables.


*     **********************************
*     Read any linear constraint inputs.
*     **********************************

      READ (LREAD, *) ! Linear constraints
      READ (LREAD, *) ! Headers

      M = 0
      DO K = 1, NLINCON_RD

         READ (LREAD, *, IOSTAT=IOS) J, LINCON_TYPE, BOUNDLO, BOUNDUP,
     >                               TLINCON, ILINCON
         IF (IOS /= 0) THEN
            WRITE (LWRIT, 1010) 'Error reading linear constraints.'
            GO TO 999
         END IF

         IF (BOUNDLO == -999.) BOUNDLO = -BIGBND
         IF (BOUNDUP ==  999.) BOUNDUP =  BIGBND

*        Warning:  # pitch constraints = NDV_ALPHA - 1; likewise for roll:

         IF (LINCON_TYPE == 'PITCH ' .AND. NDV_ALPHA > 0) THEN
            IF (M /= 0) THEN
               WRITE (LWRIT, '(/, A)')
     >            ' PITCH constraint must precede ROLL.'
               GO TO 999
            END IF
            DO J = KNOT1_A, NK_ALPHA - 1
               M = M + 1
               LCTYPE(M) = LINCON_TYPE
               TLCON(M)  = TLINCON
               ILCON(M)  = ILINCON
               BLLIN(M)  = BOUNDLO / VSCALE(J)
               BULIN(M)  = BOUNDUP / VSCALE(J)
            END DO

         ELSE IF (LINCON_TYPE == 'ROLL  ' .AND. NDV_BANK > 0) THEN
            DO J = KNOT1_B, NK_BANK - 1
               M = M + 1
               LCTYPE(M) = LINCON_TYPE
               TLCON(M)  = TLINCON
               ILCON(M)  = ILINCON
               BULIN(M)  = BOUNDUP / VSCALE(J+NDV_ALPHA)
               BLLIN(M)  = -BULIN(M)
            END DO
         END IF

      END DO


*     *************************************
*     Read any nonlinear constraint inputs:
*     *************************************

      READ (LREAD, *) ! Nonlinear constraints
      READ (LREAD, *) ! Headers

      K = NDV + NCLIN ! Offset for the bounds

      SHOW_SHORTFALL = RHO_TABORT /= ZERO .OR. RHO_TARGET_POINT /= ZERO
      SHOW_APOAPS    = .FALSE.
      SHOW_ECCENT    = RHO_ECCENTRICITY /= ZERO
      SHOW_INCLIN    = .FALSE.

      DISTRIBUTED_HEATING = RHO_HEATLD /= ZERO

      DO M = 1, NCNLN

         READ (LREAD, *, IOSTAT=IOS) J, NLCTYPE(M), BOUNDLO, BOUNDUP,
     >      XNLCON(M), INLCON(M), JNLCON(M), SNLCON(M)

         IF (IOS /= 0) THEN
            WRITE (LWRIT, 1010) 'Error reading nonlinear constraints.'
            GO TO 999
         END IF

*        Scaling of the nonlinear constraints differs from scaling of
*        the variables and linear constraints:

         IF (BOUNDLO == -999.) THEN
            BOUNDLO = -BIGBND
         ELSE
            BOUNDLO = BOUNDLO * SNLCON(M)
         END IF
         BL(M+K) = BOUNDLO

         IF (BOUNDUP ==  999.) THEN
            BOUNDUP = BIGBND
         ELSE
            BOUNDUP = BOUNDUP * SNLCON(M)
         END IF
         BU(M+K) = BOUNDUP

         IF (NLCTYPE(M) == 'ENDLON' .OR.
     >       NLCTYPE(M) == 'ENDIST') THEN
            SHOW_SHORTFALL = .TRUE.
         ELSE IF (NLCTYPE(M) == 'APOAPS') THEN
            SHOW_APOAPS = .TRUE.
         ELSE IF (NLCTYPE(M) == 'ECCENT') THEN
            SHOW_ECCENT = .TRUE.
         ELSE IF (NLCTYPE(M) == 'INCLIN') THEN
            SHOW_INCLIN = .TRUE.
         END IF

*        APC constraint?

         IF (NLCTYPE(M)(1:3) == 'APC') THEN
            NJOURNAL = MAX (NJOURNAL, 8)  ! 0, 8 or 10 are the only cases now

            OPEN (LDAT, FILE='apc.dat', STATUS='OLD', IOSTAT=IOS)

            IF (IOS /= 0) THEN
               WRITE (LWRIT, 1010) 'Unable to open apc.dat.'
               GO TO 999
            END IF

            READ (LDAT, '(A)') ! Title
            READ (LDAT, *) N_APC

            ALLOCATE (V_APC(N_APC), H_APC(N_APC))

            READ (LDAT, *) (V_APC(I), H_APC(I), I = 1, N_APC)
            CLOSE (LDAT)

            IF (V_APC(1) < V_APC(N_APC)) THEN   ! Ensure vel. & alt. decrease
               CALL RVERSE (N_APC, V_APC, V_APC)
               CALL RVERSE (N_APC, H_APC, H_APC)
            END IF
         END IF

*        Upper corridor constraint?

         IF (NLCTYPE(M) == 'UPPERC') THEN
            NJOURNAL = MAX (NJOURNAL, 8)

            OPEN (LDAT, FILE='upper_corridor.dat', STATUS='OLD')

            IF (IOS /= 0) THEN
               WRITE (LWRIT, 1010) 'Unable to open upper_corridor.dat.'
               GO TO 999
            END IF

            READ (LDAT, '(A)') ! Title
            READ (LDAT, *) N_UPPER_C

            ALLOCATE (V_UPPER_C(N_UPPER_C), H_UPPER_C(N_UPPER_C))

            READ (LDAT, *) (V_UPPER_C(I), H_UPPER_C(I), I = 1,N_UPPER_C)
            CLOSE (LDAT)

            IF (V_UPPER_C(1) < V_UPPER_C(N_UPPER_C)) THEN ! Want V descending
               CALL RVERSE (N_UPPER_C, V_UPPER_C, V_UPPER_C)
               CALL RVERSE (N_UPPER_C, H_UPPER_C, H_UPPER_C)
            END IF
         END IF

*        Descent rate constraint?

         IF (NLCTYPE(M) == 'DESCNT') THEN
            NJOURNAL = MAX (NJOURNAL, 8) ! We need time as well as APC data
         END IF                          ! (Abandon original NJOURNAL = 3 case.)

*        Dynamic pressure constraint?

         IF (NLCTYPE(M) == 'DYN_PR') THEN
            NJOURNAL = MAX (NJOURNAL, 8) ! Include time, Mach, Alpha, Q, G, CM
         END IF

*        Q * Alpha constraint?

         IF (NLCTYPE(M) == 'QALPHA') THEN
            NJOURNAL = MAX (NJOURNAL, 8)
         END IF

*        Distributed heating calculations?

         IF (NLCTYPE(M)(1:4) == 'SURF' .OR.
     >       NLCTYPE(M)      == 'HEATLD') THEN
            DISTRIBUTED_HEATING = .TRUE.
            NJOURNAL = MAX (NJOURNAL, 8)
         END IF

*        Acceleration constraint?

         IF (NLCTYPE(M) == 'ACCEL ') THEN
            NJOURNAL = MAX (NJOURNAL, 8)
         END IF

*        Trim constraint?

         IF (NLCTYPE(M) == 'TRIM  ') THEN
            NJOURNAL = MAX (NJOURNAL, 8)
         END IF

      END DO

*     Alternatively, ...

      IF (RHO_TRIM      > ZERO) NJOURNAL = MAX (NJOURNAL, 8)
      IF (RHO_LST_SQRS /= ZERO) NJOURNAL = 10  ! Need lat. & long.
      IF (MODE_VAR == 5 .OR. MODE_VAR == 6) NJOURNAL = 10 ! Keep it simple

      IF (NJOURNAL /= 0) THEN
         DEALLOCATE (TRAJ_STEPS)
         ALLOCATE (TRAJ_STEPS(MAX_STEP,NJOURNAL+1)) ! +1 for derived quantities
      END IF


*     ******************
*     Target trajectory?
*     ******************

      IF (RHO_LST_SQRS /= ZERO) THEN

         CALL GET_TARGET ()     ! Internal procedure below

         IF (IER /= 0) GO TO 999
      END IF


*     *******************************
*     Set up an aerothermal database?
*     *******************************

      IF (DISTRIBUTED_HEATING) THEN

         CALL GET_AEROTHERMAL_DATA ()  ! Internal procedure

         IF (IOS /= 0) GO TO 999
      END IF


*     *************************************************************
*     Set up the linear constraint matrix and corresponding bounds.
*     *************************************************************

      IF (NCLIN > 0) CALL SETLCON (BIGBND)

      NDITER = 0
      NCALL  = -3 ! Signals the first objective function call

      CALL SET_NPOPT ()  ! Internal set-up procedure below

      IF (INFORM /= 0) GO TO 999


*     **********************************************************
*     Constrained optimization, with possible warm restart loop:
*     **********************************************************

  800 CONTINUE

         CALL NPOPT (NDV, NCLIN, NCNLN, NROWA, NROWJ, NROWR,
     >               AMAT, BL, BU,
     >               CONFUN, OBJFUN,
     >               INFORM, MAJITS, ISTATE,
     >               CVEC, CJAC, CLAMDA, OBJPAR, GRAD, RMAT, V,
     >               IWVEC, LENIW, WVEC, LENW)

*        Switch to central differencing to improve convergence?

         IF (NITMAX > 0 .AND. NNDIFF /= 3 .AND.
     >       TIGHTEN > ZERO .AND. TIGHTEN < ONE) THEN

            NNDIFF  = 3
            INFORM  = 0
            TOLNLIN = TOLNLIN * TIGHTEN
            TOLOPT  = TOLOPT  * SQRT (TIGHTEN)
            STEPLIM = MAX (STEPLIM * TIGHTEN, 0.001)

            WRITE (LWRIT, '(/, A, /, A, A, 1P, E18.10)')
     >         ' Perform warm restart with tighter tolerances.',
     >         ' EPSOBJ ** (-1./ 6.) is applied to Hs',
     >         ' for central differencing:', EPSM6TH ! See CHECK_INPUTS

            WRITE (LNPOUT, '(/, A, /)') ' Tighten tolerances:'

            CALL SNSETR ('Major feasibility tolerance', TOLNLIN,
     >                   LNPOUT, LNPOUT, INFORM, CW, LENCW,
     >                   IWVEC, LENIW, WVEC, LENW)
            CALL SNSETR ('Optimality tolerance', TOLOPT,
     >                   LNPOUT, LNPOUT, INFORM, CW, LENCW,
     >                   IWVEC, LENIW, WVEC, LENW)
            CALL SNSETR ('Major step limit', STEPLIM,
     >                   LNPOUT, LNPOUT, INFORM, CW, LENCW,
     >                   IWVEC, LENIW, WVEC, LENW)
            CALL SNSET  ('Warm start',
     >                   LNPOUT, LNPOUT, INFORM, CW, LENCW,
     >                   IWVEC, LENIW, WVEC, LENW)
            GO TO 800

         END IF


      CONTINUE ! Recover the best trajectory results + plots

      IF (NITMAX > 0) THEN

         NCALL = -4

         CALL SOLVE (NDV, V, OBJPAR, NCALL, FAIL)

      END IF

      IF (DISTRIBUTED_HEATING) THEN ! Save the time histories of the surface pts

         CALL PUT_AEROTHERMAL_DATA ()

      END IF

      WRITE (LWRIT, '(/, A, I7, //, A)')
     >   ' Total # trajectory calculations: ', NFTOTL,
     >   ' Normal termination.'


  999 CALL SECOND (CPU2)

      CALL CPUTIME (CPU0, CPU2, PROGNAME, 'for this run', LWRIT)

*     Formats:

 1010 FORMAT (/, ' Traj_opt: ', A)


*     **************************************************************************
*
*                  Internal procedures for Traj_opt main program,
*                              in alphabetical order:
*
*                               GET_AEROTHERMAL_DATA
*                               GET_TARGET
*                               PUT_AEROTHERMAL_DATA
*                               SET_NPOPT
*                               SET_VARIABLES
*
*     **************************************************************************

      CONTAINS

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         SUBROUTINE GET_AEROTHERMAL_DATA () ! Read aerothermal database, etc.

         INTEGER :: LINE

         OPEN (LDAT, FILE='aerothermal.database', STATUS='OLD',
     >         IOSTAT=IOS)

         IF (IOS /= 0) THEN
            WRITE (LWRIT, '(/, 2A)')
     >         ' GET_AEROTHERMAL_DATA:',
     >         ' Unable to open aerothermal.database.'
            RETURN
         END IF

         READ (LDAT, '(A)') ! Title
         READ (LDAT, *) N_MACH_DH
         READ (LDAT, *) N_ALPHA_DH
         READ (LDAT, *) N_QBAR_DH
         READ (LDAT, *) N_SURF_TYPES
         READ (LDAT, *) N_SURF_POINTS

         LINE = 6

         ALLOCATE
     >     (MACH_DH(N_MACH_DH), ALPHA_DH(N_ALPHA_DH),
     >      QBAR_DH(N_QBAR_DH), ATDB
     >      (N_SURF_TYPES,N_SURF_POINTS,N_MACH_DH,N_ALPHA_DH,N_QBAR_DH))

         DO J = 1, N_SURF_POINTS
            READ (LDAT, *) ! SURFACE POINT #
            READ (LDAT, *) ! Column headers
            LINE = LINE + 2
            DO M = 1, N_QBAR_DH
               DO L = 1, N_ALPHA_DH
                  DO K = 1, N_MACH_DH

                     LINE = LINE + 1
                     READ (LDAT, *, IOSTAT=IOS)
     >                  MACH_DH(K), ALPHA_DH(L), QBAR_DH(M),
     >                  ATDB(1:N_SURF_TYPES,J,K,L,M)

                     IF (IOS /= 0) THEN
                        WRITE (LWRIT, '(/,2A,2I6, /,2A,4I6, /,A,I7)')
     >                  ' Trouble reading aerothermal.database.',
     >                  ' # surface data types, # surface points: ',
     >                  N_SURF_TYPES, N_SURF_POINTS,
     >                  ' Surface point index, Mach index,',
     >                  ' alpha index, Qbar index: ', J, K, L, M,
     >                  ' Database line number:', LINE
                        RETURN
                     END IF

                  END DO
               END DO
            END DO
         END DO

         CLOSE (LDAT)

         ! Read upper bounds for all quantities at all surface points:

         OPEN (LDAT, FILE='traj_opt.surface_bounds', STATUS='OLD',
     >         IOSTAT=IOS)

         IF (IOS /= 0) THEN
            WRITE (LWRIT, '(/, 2A)')
     >         ' GET_AEROTHERMAL_DATA:',
     >         ' Unable to open traj_opt.surface_bounds.'
            RETURN
         END IF

         READ (LDAT, '(A)') ! Title
         READ (LDAT, *) I, J

         IF (I /= N_SURF_TYPES .OR. J /= N_SURF_POINTS) THEN
            WRITE (LWRIT, '(/, 2A, /, (A, 2I6))')
     >         ' Mismatch between surface bound counts and',
     >         ' aerothermal database: ',
     >         ' traj_opt.surface_bounds: ', I, J,
     >         ' aerothermal.database:    ', N_SURF_TYPES, N_SURF_POINTS
            IOS = -1
            RETURN
         END IF

         ALLOCATE (SURF_BOUNDS (0:N_SURF_TYPES, N_SURF_POINTS))

         READ (LDAT, *, IOSTAT=IOS) SURF_BOUNDS ! Element 0 contains point wts.

         IF (IOS /= 0) THEN
            WRITE (LWRIT, '(/, 2A, 2I6, /, 2A, 4I6)')
     >         ' Trouble reading surface heating controls in',
     >         ' traj_opt.surface_bounds: ', ' # types, # pts.: ', I, J
            RETURN
         END IF

         CLOSE (LDAT)

         ALLOCATE
     >     (SURF_QUANTITIES(MAX_STEP,N_SURF_TYPES,N_SURF_POINTS),
     >      SURF_VIOLATIONS(MAX_STEP,N_SURF_TYPES,N_SURF_POINTS),
     >      SURF_PEAKS              (N_SURF_TYPES,N_SURF_POINTS),
     >      SURF_VALUES             (N_SURF_TYPES,N_SURF_POINTS),
     >      SURF_HOT_AREAS                       (N_SURF_POINTS),
     >      SURF_HEAT_LOADS                      (N_SURF_POINTS))

         END SUBROUTINE GET_AEROTHERMAL_DATA

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         SUBROUTINE GET_TARGET () ! Read target trajectory; isolate other crud

         REAL R

         OPEN (LDAT, FILE='traj_opt.target_trajectory', STATUS='OLD',
     >         IOSTAT=IOS)
         IF (IOS /= 0) THEN
            WRITE (LWRIT, '(/, A)')
     >         ' GET_TARGET: Unable to open traj_opt.target_trajectory.'
            IER = -1
            RETURN
         END IF

         ALLOCATE (TRAJ_TARGET(MAX_STEP,4))

         READ (LDAT, '(A)')  ! Skip headers
         N_TARGET = MAX_STEP ! In case the limit is reached

         DO I = 1, MAX_STEP  ! Read T, lat, long, altitude
            READ (LDAT, *, IOSTAT=IOS) TRAJ_TARGET(I, 1:4)
            IF (IOS /= 0) THEN
               N_TARGET = I - 1
               EXIT
            END IF
         END DO

         CLOSE (LDAT)

         WRITE (LWRIT, '(/, A, I6)')
     >      ' GET_TARGET:  Points in target trajectory: ', N_TARGET

         R = ONE / TRAJ_TARGET(N_TARGET,1) ! Last target time
         NORMALIZE = RHO_LST_SQRS < ZERO

         IF (NORMALIZE) THEN  ! A negative multiplier was entered
            RHO_LST_SQRS = ABS (RHO_LST_SQRS)
            TRAJ_TARGET(1:N_TARGET,1) = R * TRAJ_TARGET(1:N_TARGET,1)
         END IF

         IER = 0

         END SUBROUTINE GET_TARGET

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         SUBROUTINE PUT_AEROTHERMAL_DATA () ! Save heating time history

         OPEN (LDAT, FILE='traj_opt.temperatures', STATUS='UNKNOWN')

         DO K = 1, LAST_STEP
            WRITE (LDAT, '(F6.1, 24F9.2)') TRAJ_STEPS(K,3),
     >         SURF_QUANTITIES(K,1,1:N_SURF_POINTS)
         END DO

         CLOSE (LDAT)

         IF (N_SURF_TYPES > 1) THEN

            OPEN (LDAT, FILE='traj_opt.fluxes', STATUS='UNKNOWN')

            DO K = 1, LAST_STEP
               WRITE (LDAT, '(F6.1, 24F8.3)') TRAJ_STEPS(K,3),
     >            SURF_QUANTITIES(K,2,1:N_SURF_POINTS)
            END DO

            CLOSE (LDAT)
         END IF

         END SUBROUTINE PUT_AEROTHERMAL_DATA

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         SUBROUTINE SET_NPOPT () ! Removes NPOPT set-up clutter from main prog.

         STEPLIM = ONE      ! 1. means allow a 100% change in the variables
         TOLOPT  = 1.E-4
         TOLLIN  = TOLOPT
         TOLNLIN = TOLOPT
         TIGHTEN = ZERO
         PENPARAM= ZERO
         LEVELVER= -1
         MAJORPL = 10
         MINORPL = 0
         NPRFREQ = 100
         MINORIL = 2000
         NFTOTL  = 0        ! Incremented in CONFUN and/or OBJFUN

         ! Look for an optional namelist at the end of standard input:

         READ (LREAD, NPOPTIONS, ERR=700, END=700)

  700    CLOSE (LREAD)

         ! Run-time set-up of the optional inputs for NPOPT:

         OPEN (LNPOPT, FILE='traj_opt.nps', STATUS='UNKNOWN')

         WRITE (LNPOPT, '(A)')
     >      'Begin',
     >      'Cold Start',
     >      'Nonderivative line search',
     >      'Hessian       Full memory'
         WRITE (LNPOPT, '(A, I6)')
     >      'Major iteration limit          ', NITMAX,
     >      'Minor iteration limit          ', MINORIL,
     >      'Major print Level              ', MAJORPL,
     >      'Minor print Level              ', MINORPL,
     >      'Print frequency                ', NPRFREQ,
     >      'Verify level                   ', LEVELVER,
     >      'Derivative level               ', 3
         WRITE (LNPOPT, '(A, E12.3)')
     >      'Crash tolerance                ', ZERO,
     >      'Function precision             ', EPSOBJ,
     >      'Infinite bound                 ', BIGBND,
     >      'Infinite step                  ', STEPMX,
     >      'Major step limit               ', STEPLIM,
     >      'Minor feasibility tolerance    ', TOLLIN,
     >      'Major feasibility tolerance    ', TOLNLIN,
     >      'Optimality tolerance           ', TOLOPT,
     >      'Linesearch tolerance           ', ETA,
     >      'Penalty parameter              ', PENPARAM,
     >      'End'

         REWIND (LNPOPT)

         ! NPOPT writes to standard output (unit 6 = LNPOUT).
         ! Suppress optional screen output:

         CALL NPINIT (LNPOUT,      0, IWVEC, LENIW, WVEC, LENW)
         CALL NPSPEC (LNPOPT, INFORM, IWVEC, LENIW, WVEC, LENW)

         CLOSE (LNPOPT)

         WRITE (LWRIT, '(//, A)')

         IF (INFORM /= 0) THEN
            WRITE (LWRIT, '(/, 2A)') ' SET_NPOPT: ',
     >         'Trouble with the run-time traj_opt.nps file.'
            RETURN
         END IF

*        Check for a second options file that can be edited:

         OPEN (LNPOPT, FILE='npopt.specs', STATUS='OLD',
     >         IOSTAT=IOS)

         IF (IOS == 0) THEN

            CALL NPSPEC (LNPOPT, INFORM, IWVEC, LENIW, WVEC, LENW)

            CLOSE (LNPOPT)

            IF (INFORM /= 0) THEN
               WRITE (LWRIT, '(/, 2A)') ' SET_NPOPT: ',
     >            'Trouble with the optional npopt.specs file.'
            END IF

         END IF

         END SUBROUTINE SET_NPOPT

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         SUBROUTINE SET_VARIABLES () ! Set up variables, scaling, bounds, etc.

         IF (MODE_VAR <= 2 .OR.  ! Continuously variable alpha and/or bank
     >       MODE_VAR == 5) THEN ! Beta variables handled below

            I = 1

            IF (NDV_ALPHA > 0) THEN
               DO J = KNOT1_A, NK_ALPHA
                  VTYPE(I)  = 'A_KNOT'
                  VSCALE(I) = 20.
                  V(I)      = ALPHA_K(J)/ VSCALE(I)
                  AITCH(I)  = H_ALPHA   / VSCALE(I)
                  BL(I)     = BL_ALPHA  / VSCALE(I)
                  BU(I)     = BU_ALPHA  / VSCALE(I)
                  I = I + 1
               END DO
            END IF

            IF (NDV_BANK > 0) THEN
               DO J = KNOT1_B, NK_BANK
                  VTYPE(I)  = 'B_KNOT'
                  VSCALE(I) = 100.
                  V(I)      = BANK_K(J)/ VSCALE(I)
                  AITCH(I)  = H_BANK   / VSCALE(I)
                  BL(I)     = BL_BANK  / VSCALE(I)
                  BU(I)     = BU_BANK  / VSCALE(I)
                  I = I + 1
               END DO
            END IF

            IF (MODE_VAR == 5) THEN ! Sideslip is handled
               DO J = 1, NK_BETA
                  VTYPE(I)  = 'S_KNOT'
                  VSCALE(I) = 1.
                  V(I)      = BETA_K(J)/ VSCALE(I)
                  AITCH(I)  = H_BETA   / VSCALE(I)
                  BL(I)     = BL_BETA  / VSCALE(I)
                  BU(I)     = BU_BETA  / VSCALE(I)
                  I = I + 1
               END DO
            END IF

         ELSE IF (MODE_VAR <= 4) THEN ! 3 or 4 = step function mode

*           V(*) and VSCALE(*) have been set up by CALL ANGLE_STEPS (0, ...)

            DO I = 1, NSTEPS
               J = I + NSTEPS
               VTYPE(I) = 'A_STEP'
               VTYPE(J) = 'B_STEP'
               AITCH(I) = H_ALPHA  / VSCALE(I)
               AITCH(J) = H_BANK   / VSCALE(J)
               BL(I)    = BL_ALPHA / VSCALE(I)
               BL(J)    = BL_BANK  / VSCALE(J)
               BU(I)    = BU_ALPHA / VSCALE(I)
               BU(J)    = BU_BANK  / VSCALE(J)
            END DO

            IF (MODE_VAR == 4) THEN

               DT = MAX ((BU_ALPHA - BL_ALPHA) / PITCH_RATE,
     >                   (BU_BANK  - BL_BANK)  /  ROLL_RATE) + ONE

               DO I = 2 * NSTEPS + 1, NDV
                  VTYPE(I) = 'T_STEP'
                  AITCH(I) = SQRT (EPSOBJ)
                  BL(I)    = DT / VSCALE(I)
                  BU(I)    = REAL (NSTEPS - 1) ! Total duration / Vscale
               END DO

            END IF

         ELSE IF (MODE_VAR == 6) THEN ! CL, CD, CY are functions of time

            I = 1

            DO J = 1, NK_CL
               VTYPE(I)  = 'L_KNOT'
               VSCALE(I) = 1.
               V(I)      = CL_K(J)/ VSCALE(I)
               AITCH(I)  = H_CL   / VSCALE(I)
               BL(I)     = BL_CL  / VSCALE(I)
               BU(I)     = BU_CL  / VSCALE(I)
               I = I + 1
            END DO

            DO J = 1, NK_CD
               VTYPE(I)  = 'D_KNOT'
               VSCALE(I) = 1.
               V(I)      = CD_K(J)/ VSCALE(I)
               AITCH(I)  = H_CD   / VSCALE(I)
               BL(I)     = BL_CD  / VSCALE(I)
               BU(I)     = BU_CD  / VSCALE(I)
               I = I + 1
            END DO

            DO J = 1, NK_CY
               VTYPE(I)  = 'Y_KNOT'
               VSCALE(I) = 1.
               V(I)      = CY_K(J)/ VSCALE(I)
               AITCH(I)  = H_CY   / VSCALE(I)
               BL(I)     = BL_CY  / VSCALE(I)
               BU(I)     = BU_CY  / VSCALE(I)
               I = I + 1
            END DO

         END IF

*        The details of STEPMX usage are obscure - best to keep it very big.

         STEPMX = (STEPMX + STEPMX) / (VSCALE(1) + VSCALE(NDV_ANGLES))

         END SUBROUTINE SET_VARIABLES

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      END PROGRAM TRAJ_OPT

********************************************************************************

      SUBROUTINE ANGLE_STEPS (MODE, MODE_VAR, NSTEPS,
     >                        ALPHA, BANK, DURATION,
     >                        PITCH_RATE, ROLL_RATE,
     >                        ALPHA_K, BANK_K, TK_ALPHA, TK_BANK,
     >                        NVAR, VAR, VSCALE, LUNIO, LUNERR)

*     This routine manages step-function-type variation of angle of attack and
*     roll angle as needed for closed-loop real-time trajectory optimization.
*     It converts between step-function format and standard spline knot format
*     for alpha and bank as functions of time.  Alpha and bank are assumed to be
*     adjusted simultaneously with a small time lag suited to monotonic splines.
*
*     "traj_opt.steps" is read and scaled variables are output   (MODE = 0);
*     "traj.alp" and "traj.bnk" are also written for Traj        (MODE = 0);
*     alpha and bank arrays in Traj spline form are set up   (MODE = 0 : 2);
*     "traj_opt.steps_opt" is written from the input variables   (MODE = 2).
*
*     07/05/02  David Saunders  Initial implementation.
*
********************************************************************************

      IMPLICIT NONE

*     Arguments:

      INTEGER, INTENT (IN) ::
     >   MODE                   ! MODE = 0 means initialize: read, unpack, scale
                                !      = 1 unscale variables as expanded angles
                                !      = 2 as for 1 + write packed alpha & bank
      INTEGER, INTENT (IN) ::
     >   MODE_VAR               ! MODE_VAR = 3 means durations are fixed;
                                !          = 4 means durations are further
                                !            variables following alpha & bank

      INTEGER, INTENT (INOUT) ::
     >   NSTEPS                 ! Number of steps in time variation

      REAL, INTENT (OUT) ::
     >   ALPHA(*), BANK(*),     ! Angles for alpha and bank steps, degrees
     >   DURATION(*)            ! Corresponding durations, seconds

      REAL, INTENT (IN) ::
     >   PITCH_RATE,            ! Average rates used to calculate finite time
     >   ROLL_RATE              ! intervals between adjacent step values;
                                ! units are degrees per second
      REAL, INTENT (OUT) ::
     >   ALPHA_K(*),            ! Alpha & bank control pts. in spline form
     >   BANK_K(*),
     >   TK_ALPHA(*),           ! Times for alpha/bank controls in spline form
     >   TK_BANK(*)

      INTEGER, INTENT (INOUT) ::
     >   NVAR                   ! Number of optimization variables

      REAL, INTENT (INOUT) ::
     >   VAR(*),                ! Variables seen by the optimizer
     >   VSCALE(*)              ! Corresponding scale factors

      INTEGER, INTENT (IN) ::
     >   LUNIO,                 ! Logical unit for packed step representation
     >   LUNERR                 ! Logical unit for error messages

*     Local constants:

      REAL, PARAMETER ::
     >   DT_MIN      = 0.01,    ! In case adjacent steps are equal
     >   SCALE_ALPHA = 20.,     ! As for normal Traj_opt scaling
     >   SCALE_BANK  = 100.,
     >   ZERO        = 0.

*     Local variables:

      INTEGER
     >   I, IOS, J, K

      INTEGER, SAVE ::
     >   NSPLINE

      REAL
     >   DT

      REAL, SAVE ::
     >   SCALE_DURATION

*     Execution:

      SELECT CASE (MODE)

                  ! -------------------------------------------------
         CASE (0) ! Initialize optimization variables & spline forms.
                  ! -------------------------------------------------

            OPEN (LUNIO, FILE='traj_opt.steps', STATUS='OLD',
     >            IOSTAT=IOS)
            IF (IOS /= 0) GO TO 800

            READ (LUNIO, '(A)', IOSTAT=IOS) ! Title
            IF (IOS /= 0) GO TO 810

            READ (LUNIO, *, IOSTAT=IOS) NSTEPS
            IF (IOS /= 0) GO TO 810

            READ (LUNIO, *, IOSTAT=IOS)
     >         (ALPHA(I), BANK(I), DURATION(I), I = 1, NSTEPS)
            IF (IOS /= 0) GO TO 810

            CLOSE (LUNIO)

*           Expand the angle steps into the form suited to monotonic splines:

            NSPLINE = NSTEPS * 2
            NVAR    = NSPLINE

            VSCALE(1)        = SCALE_ALPHA
            VAR(1)           = ALPHA(1) / SCALE_ALPHA
            ALPHA_K(1:2)     = ALPHA(1)
            TK_ALPHA(1)      = ZERO
            TK_ALPHA(2)      = DURATION(1)

            VSCALE(NSTEPS+1) = SCALE_BANK
            VAR(NSTEPS+1)    = BANK(1)  / SCALE_BANK
            BANK_K(1:2)      = BANK(1)
            TK_BANK(1)       = ZERO
            TK_BANK(2)       = DURATION(1)

            SCALE_DURATION   = ZERO

            DO I = 2, NSTEPS
               J = I + I
               VSCALE(I)        = SCALE_ALPHA
               VAR(I)           = ALPHA(I) / SCALE_ALPHA
               ALPHA_K(J-1)     = ALPHA(I)
               ALPHA_K(J)       = ALPHA(I)
               DT               = ABS (ALPHA(I) - ALPHA(I-1)) /
     >                            PITCH_RATE
               TK_ALPHA(J-1)    = TK_ALPHA(J-2) + MAX (DT, DT_MIN)
               TK_ALPHA(J)      = TK_ALPHA(J-2) + DURATION(I)

               VSCALE(I+NSTEPS) = SCALE_BANK
               VAR(I+NSTEPS)    = BANK(I) / SCALE_BANK
               BANK_K(J-1)      = BANK(I)
               BANK_K(J)        = BANK(I)
               DT               = ABS (BANK(I)  - BANK(I-1)) /
     >                            ROLL_RATE
               TK_BANK(J-1)     = TK_BANK(J-2)  + MAX (DT, DT_MIN)
               TK_BANK(J)       = TK_ALPHA(J)

               SCALE_DURATION   = DURATION(I-1) + SCALE_DURATION
            END DO

*           Write "traj.alp" and "traj.bnk" for Traj initialization:

            OPEN  (LUNIO, FILE='traj.alp', STATUS='UNKNOWN')
            WRITE (LUNIO, '(A, /, I3)') 'Alpha knots', NSPLINE
            WRITE (LUNIO, '(1P, 2E14.6)')
     >         (TK_ALPHA(I), ALPHA_K(I), I = 1, NSPLINE)
            CLOSE (LUNIO)

            OPEN  (LUNIO, FILE='traj.bnk', STATUS='UNKNOWN')
            WRITE (LUNIO, '(A, /, I3)') 'Bank knots', NSPLINE
            WRITE (LUNIO, '(1P, 2E14.6)')
     >         (TK_BANK(I), BANK_K(I), I = 1, NSPLINE)
            CLOSE (LUNIO)

            IF (MODE_VAR == 4) THEN ! Durations are variables too
               J = NVAR
               NVAR = J + NSTEPS - 1
               SCALE_DURATION   = SCALE_DURATION / REAL (NSTEPS - 1)
               VAR(J+1:NVAR)    = DURATION(1:NSTEPS-1) / SCALE_DURATION
               VSCALE(J+1:NVAR) = SCALE_DURATION
               WRITE (LUNERR, '(/, A, 3F14.8)')
     >            ' Scale factors for alpha, bank, duration: ',
     >            SCALE_ALPHA, SCALE_BANK, SCALE_DURATION
            END IF

                     ! ---------------------------------------------
         CASE (1, 2) ! Set up the alpha & bank steps in spline form.
                     ! ---------------------------------------------

            IF (MODE_VAR == 4) THEN ! Control point times are varying

               J = 2 * NSTEPS
               DO I = 1, NSTEPS - 1
                  DURATION(I) = VAR(I+J) * VSCALE(I+J)
               END DO

            END IF

            ALPHA_K(1)  = VAR(1) * SCALE_ALPHA
            ALPHA_K(2)  = ALPHA_K(1)
***         TK_ALPHA(1) = ZERO
            TK_ALPHA(2) = DURATION(1)

            BANK_K(1)   = VAR(NSTEPS+1) * SCALE_BANK
            BANK_K(2)   = BANK_K(1)
***         TK_BANK(1)  = ZERO
            TK_BANK(2)  = DURATION(1)

            DO I = 2, NSTEPS
               J = I + I
               ALPHA_K(J-1)    = VAR(I) * SCALE_ALPHA
               ALPHA_K(J)      = ALPHA_K(J-1)
               DT              = ABS (ALPHA_K(J) - ALPHA_K(J-2)) /
     >                           PITCH_RATE
               TK_ALPHA(J-1)   = TK_ALPHA(J-2) + MAX (DT, DT_MIN)
               TK_ALPHA(J)     = TK_ALPHA(J-2) + DURATION(I)

               BANK_K(J-1)     = VAR(I+NSTEPS) * SCALE_BANK
               BANK_K(J)       = BANK_K(J-1)
               DT              = ABS (BANK_K(J) - BANK_K(J-2)) /
     >                           ROLL_RATE
               TK_BANK(J-1)    = TK_BANK(J-2)  + MAX (DT, DT_MIN)
               TK_BANK(J)      = TK_ALPHA(J)
            END DO

            IF (MODE == 2) THEN ! Save angles in compact form

               OPEN (LUNIO, FILE='traj_opt.steps_opt', STATUS='UNKNOWN')

               WRITE (LUNIO, '(A)')  'Optimized alpha & bank steps'
               WRITE (LUNIO, '(I3)') NSTEPS

*              May want to update ALPHA(*) and BANK(*) (= step arrays), but ...

               WRITE (LUNIO, '(1P, 3E14.6)')
     >           (ALPHA_K(I+I), BANK_K(I+I), DURATION(I), I = 1, NSTEPS)

               CLOSE (LUNIO)
            END IF

      END SELECT

      RETURN

  800 WRITE (LUNERR, '(/, A)')
     >   ' Trouble opening "traj_opt.steps".  Aborting.'
      GO TO 900

  810 WRITE (LUNERR, '(/, A)')
     >   ' Trouble reading "traj_opt.steps".  Aborting.'

  900 STOP
      
      END SUBROUTINE ANGLE_STEPS

********************************************************************************
*
      SUBROUTINE CHECK_INPUTS (IER)
*
*     Look for some of the possible input errors, and perform some other
*     initialization.  Ideally, all counts such as the number of alpha and bank
*     control points would be determined from the relevant ancillary files.
*     However, it's not that simple.  For instance, bank variation may be
*     suppressed by entering NDV_BANK = 0 (yet Traj still reads traj.bnk).
*     See also the MODE_VAR description.
*
********************************************************************************

*     Global variables:

      USE CONSTS
      USE OPT

      IMPLICIT NONE

*     Arguments:

      INTEGER, INTENT (OUT) :: IER ! 0 = No error, else quit

*     Local constants:

      INTEGER, PARAMETER :: NO_GOOD = -1
      REAL,    PARAMETER :: ZERO = 0.

*     Local variables:

      INTEGER
     >   IOS, L, N

      CHARACTER
     >   KEYWORD * 4, KEY_VALUE * 48

*     Execution:

      IER = 0

      IF (NDV_ALPHA == 0) RHO_ALPHA_TVD = ZERO
      IF (NDV_BANK  == 0) RHO_BANK_TVD  = ZERO
      IF (NDV_BETA  == 0) RHO_BETA_TVD  = ZERO

      IF (RHO_CROSS_RANGE  == ZERO .AND.
     >    RHO_DOWN_RANGE   == ZERO .AND.
     >    RHO_HEAT_LOAD    == ZERO .AND.
     >    RHO_MAX_ACCEL    == ZERO .AND.
     >    RHO_MAX_W_TEMP   == ZERO .AND.
     >    RHO_MAX_DYN_PR   == ZERO .AND.
     >    RHO_MAX_H_FLUX   == ZERO .AND.
     >    RHO_TABORT       == ZERO .AND.
     >    RHO_END_ALTITUDE == ZERO .AND.
     >    RHO_ECCENTRICITY == ZERO .AND.
     >    RHO_HEATLD       == ZERO .AND.
     >    RHO_TRIM         == ZERO .AND.
     >    RHO_TARGET_POINT == ZERO .AND.
     >    RHO_ALPHA_TVD    == ZERO .AND.
     >    RHO_BANK_TVD     == ZERO .AND.
     >    RHO_BETA_TVD     == ZERO .AND.
     >    RHO_DURATION     == ZERO .AND.
     >    RHO_LST_SQRS     == ZERO) THEN
         WRITE (LWRIT, '(/, A)') ' All RHO_* inputs are zero.'
         IER = NO_GOOD
      END IF

      EPSM6TH = EPSOBJ ** (-1./ 6.) ! NNDIFF /= 3 for the warm start option

      IF (NNDIFF == 3) THEN
         WRITE (LWRIT, '(/, A, 1P, E18.10)')
     >   ' EPSOBJ ** (-1./ 6.) applied to Hs for central differencing:',
     >     EPSM6TH
      END IF

*     Traj doesn't return the dimensions of the aero. database, so they
*     are input to Traj_opt.  (3 integers at the top of the file would
*     have saved a lot of mucking around in both packages.)
*     Check for a probable error, first by identifying the file name,
*     then by counting the number of lines in that file.  We assume that
*     Traj's control file is "traj.in", so its -s <file> switch should
*     not be used at the Traj_opt command line.

      N_COEFS = N_MACH * N_ALPHA * N_RE
      IF (MODE_VAR == 5) N_COEFS = N_COEFS * N_BETA

      OPEN (LDAT, FILE='traj.in', STATUS='OLD', IOSTAT=IOS)

      IF (IOS /= 0) THEN
         WRITE (LWRIT, '(/, A)')
     >      ' Unable to open traj.in to find aero. database file name.'
         IER = NO_GOOD
      ELSE
         DO L = 1, 999
            READ (LDAT, *, IOSTAT=IOS) KEYWORD, KEY_VALUE

            IF (IOS <  0) THEN ! EOF; no aero. database - forget the check
               CLOSE (LDAT)
               EXIT
            END IF
               
            IF (IOS /= 0) THEN
               WRITE (LWRIT, '(/, A, I6)')
     >            ' Trouble reading traj.in; line: ', L
               IER = NO_GOOD
               EXIT
            END IF

            IF (KEYWORD == 'File') THEN

               CLOSE (LDAT)
               OPEN (LDAT, FILE=KEY_VALUE, STATUS='OLD', IOSTAT=IOS)

               IF (IOS /= 0) THEN
                  WRITE (LWRIT, '(/, A, A)')
     >               ' Unable to open aero. database: ', KEY_VALUE
                  IER = NO_GOOD
                  EXIT
               END IF

               AERO_FILENAME = KEY_VALUE ! See kludge in USER, MODE_VAR = 0 case)

               N = 0
               DO WHILE (IOS == 0)
                  N = N + 1
                  READ (LDAT, '(A)', IOSTAT=IOS) ! Skip a line
               END DO
               N = N - 1

               IF (N /= N_COEFS) THEN
                  WRITE (LWRIT, '(/, A, /, A, 2I8)')
     >               ' Aero. database dimension inputs appear wrong.',
     >               ' Product of dimensions and line count: ',
     >               N_COEFS, N
                  IER = NO_GOOD
               END IF

               CLOSE (LDAT)
               EXIT

            END IF

         END DO ! Next traj.in line         

      END IF

*     Checks involving the number of optimization variables:

      NK_ALPHA   = 0 ! In case some variables are suppressed
      NK_BANK    = 0
      NK_BETA    = 0
      NDV_ANGLES = NDV_ALPHA + NDV_BANK + NDV_BETA

      SELECT CASE (MODE_VAR)

         CASE (0) ! Sensitivities to CLs and CDs in aero. database

            N = 2 * N_COEFS
            IF (NDV /= N) THEN
               WRITE (LWRIT, '(/, A, I5)')
     >            ' MODE_VAR = 0 requires NDV = 2 * # lift coefs. =', N
               IER = NO_GOOD
            END IF

         CASE (1, 2) ! Alpha and/or bank (fixed times, continuous variation)

            IF (NDV_ALPHA > 0) THEN

               CALL READ_COUNT ('traj.alp', NK_ALPHA) ! Internal procedure

               IF (NK_ALPHA - KNOT1_A + 1 /= NDV_ALPHA) THEN
                  WRITE (LWRIT, '(/, A, 3I6)')
     >              ' Mismatched NDV_ALPHA, KNOT1_A, & N in traj.alp: ',
     >              NDV_ALPHA, KNOT1_A, NK_ALPHA
                  IER = NO_GOOD
               END IF

            END IF

            IF (NDV_BANK > 0) THEN

               CALL READ_COUNT ('traj.bnk', NK_BANK)

               IF (NK_BANK - KNOT1_B + 1 /= NDV_BANK) THEN
                  WRITE (LWRIT, '(/, A, 3I6)')
     >               ' Mismatched NDV_BANK, KNOT1_B, & N in traj.bnk: ',
     >               NDV_BANK, KNOT1_B, NK_BANK
                  IER = NO_GOOD
               END IF

            END IF

            IF (MODE_VAR == 1) THEN
               IF (NDV /= NDV_ANGLES) THEN
                  WRITE (LWRIT, '(/, A, I5)')
     >              ' MODE_VAR = 1 requires NDV = NDV_(ALPHA + BANK) =',
     >               NDV_ANGLES
                  IER = NO_GOOD
               END IF
            ELSE ! MODE_VAR == 2
               N = NDV_ANGLES + 1
               IF (NDV == NDV_ANGLES) THEN
                  IF (ENTRY_DV == ZERO) THEN
                     WRITE (LWRIT, '(/, 2A)')
     >                  ' WARNING:  MODE_VAR = 2 yet NDV =',
     >                  ' NDV_(ALPHA + BANK) and ENTRY_DV = 0.'
                     IER = NO_GOOD
                  END IF
               ELSE IF (NDV /= N) THEN
                  WRITE (LWRIT, '(/, A, I5)')
     >          ' MODE_VAR = 2 requires NDV = NDV_(ALPHA + BANK) + 1 =',
     >               N
                  IER = NO_GOOD
               END IF
            END IF

         CASE (3, 4) ! Alpha and bank are step functions

            CALL READ_COUNT ('traj_opt.steps', NSTEPS)

            IF (NDV_ALPHA /= NSTEPS .AND. NDV_BANK /= NSTEPS) THEN
               WRITE (LWRIT, '(/, 2A, 3I6)')
     >            ' MODE_VAR = 3 or 4 requires NDV_ALPHA = NDV_BANK',
     >            ' = NSTEPS: ', NDV_ALPHA, NDV_BANK, NSTEPS
               IER = NO_GOOD
            END IF

            IF (KNOT1_A /= 1 .OR. KNOT1_B /= 1) THEN
               WRITE (LWRIT, '(/, 2A, 2I6)')
     >            ' MODE_VAR = 3 or 4 requires all steps active.',
     >            ' KNOT1_A, KNOT1_B: ', KNOT1_A, KNOT1_B
               IER = NO_GOOD
            END IF

            N = 2 * NSTEPS
            NK_ALPHA = N   ! After expansion of steps to standard Traj_opt form
            NK_BANK  = N

            IF (MODE_VAR == 3) THEN
               IF (NDV /= N) THEN
                  WRITE (LWRIT, '(/, A, I5)')
     >             ' MODE_VAR = 3 requires NDV = 2 * NSTEPS =', N
                  IER = NO_GOOD
               END IF
            ELSE ! MODE_VAR == 4
               N = N + NSTEPS - 1
               IF (NDV /= N) THEN
                  WRITE (LWRIT, '(/, A, I5)')
     >               ' MODE_VAR = 2 requires NDV = 3 * NSTEPS - 1 =', N
                  IER = NO_GOOD
               END IF
            END IF

         CASE (5) ! Alpha, bank, and beta are functions of time

            CALL READ_COUNT ('traj.alp', NK_ALPHA)

            IF (NK_ALPHA /= NDV_ALPHA) THEN
               WRITE (LWRIT, '(/, A, 2I6)')
     >      ' Mismatched NDV_ALPHA & N in traj.alp; use KNOT1_A = 1: ',
     >            NDV_ALPHA, NK_ALPHA
               IER = NO_GOOD
            END IF

            CALL READ_COUNT ('traj.bnk', NK_BANK)

            IF (NK_BANK /= NDV_BANK) THEN
               WRITE (LWRIT, '(/, A, 2I6)')
     >      ' Mismatched NDV_BANK & N in traj.bnk; use KNOT1_B = 1: ',
     >            NDV_BANK, NK_BANK
               IER = NO_GOOD
            END IF

            CALL READ_COUNT ('traj.bet', NK_BETA)

            IF (NK_BETA /= NDV_BETA) THEN
               WRITE (LWRIT, '(/, A, 2I6)')
     >            ' Mismatched NDV_BETA & N in traj.bet: ',
     >            NDV_BETA, NK_BETA
               IER = NO_GOOD
            END IF

         CASE (6) ! CL, CD, & CY are functions of time

            CALL READ_COUNT ('traj.lft', NK_CL)

            IF (NK_CL /= NDV_CL) THEN
               WRITE (LWRIT, '(/, A, 2I6)')
     >            ' Mismatched NDV_CL & N in traj.lft: ', NDV_CL, NK_CL
               IER = NO_GOOD
            END IF

            CALL READ_COUNT ('traj.drg', NK_CD)

            IF (NK_CD /= NDV_CD) THEN
               WRITE (LWRIT, '(/, A, 2I6)')
     >            ' Mismatched NDV_CD & N in traj.drg: ', NDV_CD, NK_CD
               IER = NO_GOOD
            END IF

            CALL READ_COUNT ('traj.sid', NK_CY)

            IF (NK_CY /= NDV_CY) THEN
               WRITE (LWRIT, '(/, A, 2I6)')
     >            ' Mismatched NDV_CY & N in traj.sid: ', NDV_CY, NK_CY
               IER = NO_GOOD
            END IF

            NDV_COEFFS = NDV_CL + NDV_CD + NDV_CY ! Equivalenced to NDV_ANGLES

         END SELECT

      IF (IER /= 0) THEN
         WRITE (LWRIT, '(/, A)') ' CHECK_INPUTS:  Aborting.' ! In main program
      END IF

*     Internal procedure for CHECK_INPUTS:

      CONTAINS

         SUBROUTINE READ_COUNT (FILENAME, NPOINTS)

!        Find the data point count on the 2nd line of the indicated file.

!        Arguments:

         CHARACTER, INTENT (IN)  :: FILENAME * (*)
         INTEGER,   INTENT (OUT) :: NPOINTS

!        Execution:

         OPEN (LDAT, FILE=FILENAME, STATUS='OLD', IOSTAT=IOS)

         IF (IOS /= 0) THEN
            WRITE (LWRIT, '(/, 2A)') ' Unable to open file ', FILENAME
            NPOINTS = 0
         ELSE
            READ (LDAT, '(A)') ! Skip the title
            READ (LDAT, *, IOSTAT=IOS) NPOINTS
            IF (IOS /= 0) THEN
               WRITE (LWRIT, '(/, 2A)')
     >            ' Error reading # points from ', FILENAME
               NPOINTS = 0
            END IF
         END IF

         CLOSE (LDAT)

         END SUBROUTINE READ_COUNT

      END SUBROUTINE CHECK_INPUTS

********************************************************************************
*
      SUBROUTINE CHKLCON
*
*     Check the optimized trajectory against the linear constraints to give
*     a more readable table than is provided by the optimizer.
*
*     06/25/96  DAS  Adaptation of SETLCON.
*     02/28/00   "   Adaptation for trajectory optimization (no linear c's.).
*     05/09/00   "   Implemented linear constraints on alpha & bank angle.
*     09/27/00   "   Implemented pitch and roll rate constraints.
*
********************************************************************************

*     Global variables:

      USE CONSTS
      USE OPT

      IMPLICIT NONE

*     Local variables:

      INTEGER   K, M, N_ALPHA_OFFSET
      REAL      BOUNDL, BOUNDU, S, VALUE
      CHARACTER FLAG*3

*     Execution:

      WRITE (LWRIT, 1010)

*     Evaluate the linear constraints in the form seen by NPOPT:

      CALL AX (NCLIN, NCLIN, NDV, AMAT, V, CONLIN)

      IF (LINCON == 3) THEN
         N_ALPHA_OFFSET = NDV_ALPHA
      ELSE
         N_ALPHA_OFFSET = 1
      END IF

      DO M = 1, NCLIN

         SELECT CASE (LCTYPE(M))

         CASE ('PITCH ')

*           Internal procedure treats pitch or roll rate, returning VALUE:

            K = KNOT1_A + M - 1

            CALL ANGLE_RATE (NK_ALPHA, TK_ALPHA, ALPHA_K)

            S = VSCALE(1)
            BOUNDL = BLLIN(M) * S
            BOUNDU = BULIN(M) * S

         CASE ('ROLL ')

            K = KNOT1_B + M - N_ALPHA_OFFSET
            CALL ANGLE_RATE (NK_BANK,  TK_BANK,  BANK_K)

            S = VSCALE(NDV_ANGLES)
            BOUNDL = BLLIN(M) * S
            BOUNDU = BULIN(M) * S

         END SELECT

         IF (VALUE < BOUNDL .OR. VALUE > BOUNDU) THEN
            FLAG = '***'
         ELSE
            FLAG = '   '
         END IF

         WRITE (LWRIT, 1030) M, LCTYPE(M), ILCON(M), VALUE,
     >                       BOUNDL, BOUNDU, CONLIN(M) * S, FLAG
      END DO

      RETURN

 1010 FORMAT (/, ' CHKLCON: Linear constraint check:', //,
     >   '    #  LCTYPE   ILCON  CURRENT VALUE',
     >   '    LOWER BOUND    UPPER BOUND           AMAT * V')
 1030 FORMAT (1X, I4, 2X, A6, I6, 2X, 3F15.6, E19.10, 2X, A3)

*     Internal procedures for CHKLCON:

      CONTAINS

         SUBROUTINE ANGLE_RATE (NPTS, TIMES, ANGLES)

!        Evaluates pitch or roll rate at Kth angle control pt., returning VALUE.

!        Arguments:

         INTEGER NPTS        ! Number of angle vs. time data points
         REAL    TIMES(NPTS)
         REAL    ANGLES(NPTS)

!        Execution:

         VALUE = (ANGLES(K+1) - ANGLES(K)) /
     >           ( TIMES(K+1) -  TIMES(K))

         END SUBROUTINE ANGLE_RATE

         subroutine ax (idim, m, n, a, x, b)

!        Computes b = Ax for m x n matrix A.

!        Arguments:

         integer, intent (in)  :: idim, m, n
         real,    intent (in)  :: a(idim,n), x(n)
         real,    intent (out) :: b(m)

!        Local variables:

         integer  i, j
         real     sum

!        Execution:

         do i = 1, m
            sum = 0.
            do j = 1, n
               sum = a(i,j) * x(j) + sum
            end do
            b(i) = sum
         end do

         end subroutine ax

      END SUBROUTINE CHKLCON

********************************************************************************
*
      SUBROUTINE COMPARE_TRAJECTORIES (NC, TC, XC, YC, ZC, NORM,
     >                                 NT, TT, XT, YT, ZT, AREA)
*
*     Calculate a measure of how well a 3-space trajectory matches a target
*     trajectory.  (X, Y, Z) may be (Lat, Long, Altitude) in any reasonable
*     units, as long as Lat and Long don't contain discontinuities around
*     0, 90, 180, or 360 degrees.
*
*     If normalized times are indicated, the entire histories of the target and
*     current trajectories are compared.  Otherwise, only the portions which
*     overlap in real time are compared.  It remains to be seen which approach
*     behaves better during optimization with respect to the control schedules.
*     Either way, the current data are evaluated at the times of the target time
*     history via local cubic spline interpolation.
*
*     Each history is assumed to imply a start at time zero, even if the first
*     data point is slightly later.  (Traj time histories omit the entry point.)
*
*     History:
*
*     03/04/03   DAS   Initial implementation (plain sum of squared deviations).
*     03/06/03    "    Normalized case works but unnormalized case bogs down;
*                      therefore, integrate the |deviations| w.r.t. time.
*
*     Author:  David Saunders, ELORET/NASA Ames Research Center.
*
********************************************************************************

      IMPLICIT NONE

*     Arguments:

      INTEGER, INTENT (IN) ::
     >   NC                 ! Number of points in the current time history

      REAL, INTENT (IN), DIMENSION (NC) ::
     >   TC, XC, YC, ZC     ! Unnormalized times & spacial coordinates for
                            ! the current trajectory
      LOGICAL, INTENT (IN) ::
     >   NORM               ! T means normalize the current data to match the
                            ! target data; TT(*) is assumed normalized already
      INTEGER, INTENT (IN) ::
     >   NT                 ! Number of point in the target time history

      REAL, INTENT (IN), DIMENSION (NT) ::
     >   TT, XT, YT, ZT     ! Unnormalized times & spacial coordinates for
                            ! the current trajectory
      REAL, INTENT (OUT) ::
     >   AREA               ! Measure of the trajectory match (or lack of it)
                            ! as described above
*     Local constants:

      REAL,      PARAMETER :: ONE = 1., ZERO = 0.
      LOGICAL,   PARAMETER :: NEW = .TRUE.
      CHARACTER, PARAMETER :: LOOSE * 1 = 'B', TIGHT * 1 = 'M'

*     Local variables:

      INTEGER ::
     >   I, NI

      REAL, ALLOCATABLE, DIMENSION (:) ::
     >   DEV, TCNORM, XI, YI, ZI

*     Procedures:

      EXTERNAL
     >   INTERVAL,          ! Search utility
     >   LCSAREAS,          ! Spline quadrature utility
     >   LCSFIT             ! Local cubic spline utility

*     Execution:

      IF (NORM) THEN ! Interpolate at all target points

         NI = NT
         ALLOCATE (TCNORM(NC), DEV(NI), XI(NI), YI(NI), ZI(NI))

         TCNORM = TC * (ONE / TC(NC))  ! Assume a start near time = 0.
                                       ! DEV serves for unused derivatives

         CALL LCSFIT (NC, TCNORM, XC, NEW, LOOSE, NI, TT, XI, DEV)
         CALL LCSFIT (NC, TCNORM, YC, NEW, LOOSE, NI, TT, YI, DEV)
         CALL LCSFIT (NC, TCNORM, ZC, NEW, LOOSE, NI, TT, ZI, DEV)

         DEALLOCATE (TCNORM)

      ELSE

         IF (TC(NC) <= TT(NT)) THEN ! Locate last target time to match

            NI = NC - 1
            CALL INTERVAL (NT, TT, TC(NC), ONE, NI)  ! NI < NT

            IF (TC(NC) == TT(NI + 1)) NI = NI + 1

         ELSE
            NI = NT
         END IF

         ALLOCATE (DEV(NI), XI(NI), YI(NI), ZI(NI))
 
         CALL LCSFIT (NC, TC, XC, NEW, LOOSE, NI, TT, XI, DEV)
         CALL LCSFIT (NC, TC, YC, NEW, LOOSE, NI, TT, YI, DEV)
         CALL LCSFIT (NC, TC, ZC, NEW, LOOSE, NI, TT, ZI, DEV)

      END IF

      DO I = 1, NI
         DEV(I) = SQRT ((XI(I) - XT(I)) ** 2 + (YI(I) - YT(I)) ** 2 +
     >                  (ZI(I) - ZT(I)) ** 2)
      END DO

*     LCSAREAS is slightly more efficient than LCSQUAD.
*     Reuse XI(*) for the partial integrals.

      CALL LCSAREAS (NI, TT, DEV, TIGHT, ZERO, XI)

      AREA = XI(NI)

      DEALLOCATE (DEV, XI, YI, ZI)

      END SUBROUTINE COMPARE_TRAJECTORIES

********************************************************************************
*
      SUBROUTINE CONFUN (MODE, NCPAR, NDVPAR, NROWJPAR, NEEDC, VPAR,
     >                   CPAR, CJACPAR, NSTATE)
*
*     CONFUN computes the nonlinear constraint functions C(*) and their
*     gradients.  Its calling sequence is that required by NPSOL/NPOPT.
*     Finite differencing uses the "H"s read with each design variable.
*
*     This version ignores NEEDC(*) -  it is more efficient to evaluate
*     all constraints in groups than to evaluate each one independently.
*
*     This version also allows CONFUN to evaluate the objective function
*     and its derivatives.  Otherwise, pure trajectory optimization repeats
*     the same trajectories for the objective as for the nonlinear constraints.
*     Since OBJFUN is called following each CONFUN call, it appears safe for
*     CONFUN to carry around the most recent objective and its gradient.
*
*     06/10/96  DAS/JJR  Initial implementation (fuel volume,
*                        floor angles, and cabin radii).
*     12/11/96    DAS    Avoid recalculating zero derivatives after the
*                        first time through.
*     02/28/00     "     Initial adaptation for trajectory optimization.
*     05/03/00     "     Constraints aren't geometric now - we need to
*                        calculate a trajectory.
*     05/15/00     "     NNGRAD = 2 or 3 option saves OBJFUN from repeating
*                        trajectories just calculated by CONFUN.
*                        Communication is via OBJ_CONFUN and the FFWD(*)
*                        array already available from the QNMDIF2 option.
*                        NCALL <= 0 is set to the OBJFUN value minus 10
*                        as a way of telling SOLVE this is CONFUN, not OBJFUN.
*     07/11/00     "     NNDIFF = 3 allows central differencing.
*     06/19/02     "     Suppressed gradient calculations for variables
*                        suitably beyond the end of the current trajectory.
*     05/20/03     "     The suppression was suppressing one too many.
*
********************************************************************************

*     Global quantities:

      USE CONSTS
      USE OPT             ! Needed for AITCH(*), etc.

      IMPLICIT NONE

*     Arguments:

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

*     Local constants:

      REAL,      PARAMETER :: HALF = 0.5, ONE = 1., ZERO = 0.
      LOGICAL,   PARAMETER :: NOPRINT = .FALSE.
      CHARACTER, PARAMETER :: SUBNAME * 6 = 'CONFUN'

*     Local variables:

      INTEGER    I, J, LASTA, LASTB, LASTY, NCALL
      REAL       HCEN, OBJFOR, OBJBAK, TIME1, TIME2, VTEMP
      LOGICAL    FAIL

*     Execution:

***      WRITE (LWRIT, '(/, A, 2I4)')
***     >   ' CONFUN: MODE and NSTATE: ', MODE, NSTATE

      IF (NSTATE == 1) THEN ! First call to CONFUN (precedes OBJFUN)

         NCALL = -13

      ELSE ! Not the first call to CONFUN

         IF (MODE == 0) THEN ! Constraints only
            NCALL = -10
         ELSE                ! Constraints + derivatives
            NCALL = -12
         END IF

      END IF

*     Calculate a trajectory:

      CALL SOLVE (NDVPAR, VPAR, OBJFOR, NCALL, FAIL)

      NFTOTL = NFTOTL + 1

      IF (FAIL) THEN
         WRITE (LWRIT, '(/, 1X, 2A)')
     >      SUBNAME, ':  Bad return from SOLVE.'
         GO TO 900
      END IF

      OBJ_CONFUN = OBJFOR ! Avoid repeat: OBJFUN uses it if CONFUN_DOES_OBJ

      CALL NLCON (NOPRINT, CPAR) ! Evaluate all nonlinear constraints, scaled

      IF (NITMAX == 0) THEN
         IF (CONFUN_DOES_OBJ) THEN
            WRITE (LWRIT, '(/, 1X, 2A)')
     >         SUBNAME, ':  Forcing exit from NPOPT (NITMAX = 0).'
            GO TO 900
         ELSE ! No point doing the constraint gradients
            CJACPAR = 999. ! Keep NPOPT happy
            GO TO 999
         END IF
      END IF

*     Constraint derivatives as well?

      IF (MODE > 0) THEN

         IF (MODE_VAR <= 5) THEN

*           Locate the last significant Alpha & bank [& Beta] control points:

            IF (NDV_ALPHA > 0) THEN
               LASTA = NDV_ALPHA - 1
               CALL INTERVAL (NK_ALPHA,TK_ALPHA, TOTAL_TIME, ONE, LASTA)
            END IF

            IF (NDV_BANK > 0) THEN
               LASTB = NDV_BANK - 1
               CALL INTERVAL (NK_BANK, TK_BANK,  TOTAL_TIME, ONE, LASTB)
            END IF

            IF (MODE_VAR == 5) THEN
               LASTY = NDV_BETA - 1
               CALL INTERVAL (NK_BETA, TK_BETA,  TOTAL_TIME, ONE, LASTY)
            END IF

         ELSE ! MODE_VAR = 6:  CL, CD, CY are the variables

            LASTA = NDV_CL - 1
            CALL INTERVAL (NK_CL, TK_CL, TOTAL_TIME, ONE, LASTA)

            LASTB = NDV_CD - 1
            CALL INTERVAL (NK_CD, TK_CD, TOTAL_TIME, ONE, LASTB)

            LASTY = NDV_CY - 1
            CALL INTERVAL (NK_CY, TK_CY, TOTAL_TIME, ONE, LASTY)

         END IF

         CALL SECOND (TIME1)

         IF (NNDIFF == 2) THEN ! Forward differencing

*           Perturb a given variable only once for all constraints:

            DO J = 1, NDV

               IF (OFF_THE_END (J)) CYCLE ! Avoid a wasted trajectory

               VTEMP   = VPAR(J)
               VPAR(J) = VTEMP + AITCH(J)

               CALL SOLVE (NDVPAR, VPAR, OBJFOR, J, FAIL) ! NCALL = J here

               NFTOTL  = NFTOTL + 1
               VPAR(J) = VTEMP

               CALL NLCON (NOPRINT, CJACPAR(1,J))

               DO I = 1, NCNLN
                  CJACPAR(I,J) = (CJACPAR(I,J) - CPAR(I)) / AITCH(J)
               END DO

               IF (CONFUN_DOES_OBJ)
     >                    FFWD(J) = (OBJFOR - OBJ_CONFUN) / AITCH(J)
            END DO

         ELSE ! NNDIFF = 3: central differencing

*           Perturb a given variable twice for all constraints:

            DO J = 1, NDV

               IF (OFF_THE_END (J)) CYCLE ! Avoid a wasted trajectory

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

               IF (CONFUN_DOES_OBJ)
     >                    FFWD(J) = (OBJFOR - OBJBAK) * HCEN

            END DO

         END IF

         CALL SECOND (TIME2)

         I = MIN (NCNLN, 10)

         IF (CONFUN_DOES_OBJ) THEN

            CALL CPUTIME (TIME1, TIME2, SUBNAME,
     >   'to evaluate derivatives of nonlinear constraints & objective',
     >         LWRIT)

            WRITE (LWRIT, '(/, A)' )
     >         ' Gradients of objective and constraints:',
     >         ' IVAR       dOBJ/dVAR   dC1/dVAR   dC2/dVAR ...'

            DO J = 1, NDV
               WRITE (LWRIT, '(I4, 1P, E17.8, 10E11.3)')
     >           J, FFWD(J), CJACPAR(1:I,J)
            END DO

         ELSE

            CALL CPUTIME (TIME1, TIME2, SUBNAME,
     >         'to evaluate nonlinear constraint derivatives', LWRIT)

            WRITE (LWRIT, '(/, A)' ) ' Jacobian matrix transpose:',
     >         ' IVAR    dC1/dVAR    dC2/dVAR ...'
            DO J = 1, NDV
               WRITE (LWRIT, '(I5, 1P, 10E12.3)') J, CJACPAR(1:I,J)
            END DO

         END IF

      END IF

      GO TO 999

  900 MODE = -2 ! NPOPT should quit

  999 RETURN

*     Internal procedure for CONFUN:

      CONTAINS

         LOGICAL FUNCTION OFF_THE_END (J)

!        Check for zero derivatives w.r.t. an angle or coef. control point
!        that is sufficiently beyond the duration of the current trajectory.

!        Argument:

         INTEGER, INTENT (IN) ::
     >      J                 ! Optimization variable number

!        Local variables:

         INTEGER
     >      JT
         LOGICAL
     >      ZEROS

!        Execution:

         ZEROS = .FALSE.

         IF (J <= NDV_ALPHA) THEN ! NDV_CL is equivalenced to NDV_ALPHA, etc.

            JT = J + KNOT1_A - 1
            IF (JT > LASTA + 2) THEN
               FFWD(J) = ZERO
               CJACPAR(1:NCNLN,J) = ZERO
               ZEROS = .TRUE.
            END IF

         ELSE IF (J <= NDV_ALPHA + NDV_BANK) THEN

            JT = J + KNOT1_B - 1
            IF (JT - NDV_ALPHA > LASTB + 2) THEN
               FFWD(J) = ZERO
               CJACPAR(1:NCNLN,J) = ZERO
               ZEROS = .TRUE.
            END IF

         ELSE IF (MODE_VAR >= 5) THEN

            IF (J - NDV_ALPHA - NDV_BANK > LASTY + 2) THEN
               FFWD(J) = ZERO
               CJACPAR(1:NCNLN,J) = ZERO
               ZEROS = .TRUE.
            END IF

         END IF

         OFF_THE_END = ZEROS

         END FUNCTION OFF_THE_END

      END SUBROUTINE CONFUN

************************************************************************
*
      SUBROUTINE CPUTIME (CPU1, CPU2, CALLER, DESCRIPTION, LUN)
*
*     CPUTIME modularizes printing of CPU time between two measures of
*     total CPU time used so far as given by the SECOND utility.
*
************************************************************************

      IMPLICIT NONE

*     Arguments:

      REAL,      INTENT (IN) :: CPU1, CPU2
      CHARACTER, INTENT (IN) :: CALLER * (*), DESCRIPTION * (*)
      INTEGER,   INTENT (IN) :: LUN

*     Local variables:

      REAL DIFF

*     Execution:

      DIFF = CPU2 - CPU1
      WRITE (LUN, '(/, 1X, 4A, F10.2)')
     >   CALLER, ': CPU secs. ', DESCRIPTION, ':   ', DIFF

      END SUBROUTINE CPUTIME

********************************************************************************
*
      SUBROUTINE CURVE_VIOLATIONS (ABOVE, N, X, Y, T, NC, XC, YC,
     >                             VALUE, XPEAK)
*
*     Calculate a measure of how much a curve defined by quantities X and Y from
*     a time history violates (or doesn't violate) a constraint curve defined in
*     the (X,Y) space of the two quantities.
*
*     For trajectory applications, X and Y might be Velocity and Altitude.  A
*     violation may be defined as being above the constraint or below it - both
*     cases are handled.
*
*     If a violation exists, it is calculated as a sum of all the elemental
*     "areas" rather than as a peak Y "distance".  Otherwise, a minimum delta Y
*     is returned rather than zero.  This change in units on either side of a
*     zero violation is considered preferable to returning zero for the no-
*     violation case because it should provide usable gradient information.
*     Experience shows that such gradient discontinuities are tolerable.
*
*     Note that while PEAK violation information on its own may well suffice,
*     it is not as satisfactory as a summing of all violations because the peak
*     location may jump around and give confusing gradients.  This can cause
*     "hiccups" near convergence.  The extra expense here should be compensated
*     for by fewer optimization iterations.
*
*     The simpler case of YC = constant (Y might be a temperature or G load) is
*     also treated - more efficiently than the general case.
*
*     Assumptions:
*
*     >  X may not be strictly monotonic.  Therefore, any area integrations are
*        performed with respect to time T, not X.
*
*     >  If NC > 1, X(*) and XC(*) are assumed to be both ordered in the same
*        general direction (even though X(*) may not be strictly monotonic).
*
*     >  For the "above" case with NC > 1, violations are not counted until
*        after the data Y drops below the constraint.  It is assumed that this
*        must happen.  This exception handles an atmospheric reentry trajectory,
*        which probably violates any upper corridor boundary at the start of the
*        data.
*
*     >  For the simpler case (NC = 1, YC = constant), "above" is assumed.
*
*     >  If NC > 1 and the constraint curve is not defined at X(1), as in the
*        case of a trajectory to the "left" of an APC, there is no good way
*        to define the constraint value that is sure to be continuous.
*        The constraint should be turned off in this case.
*
*     >  However, XC(NC) need not extend to X(N), so a search is performed in
*        the direction indicated by X(N) - X(1).  This does not apply to the
*        simple case, where X is assumed to be increasing (being Time).
*
*     >  The goal is to eliminate constraint violations, so the appropriate
*        (upper or lower) bound used by the optimizer should be zero:
*
*                                        BL          BU
*        "Above" = violation case:      -BIG         0.
*        "Below" = violation case:       0.         +BIG
*
*        Here, BIG should suit the likely delta Y in units of Y.
*
*     History:
*
*     June, 2000  DAS  Initial APC implementation ("below", but integrations
*                      were with respect to Velocity, not time step).
*     Aug., 2001   "   Upper corridor case ("above", and integrating with
*                      respect to time step).
*     Dec.  2001   "   Generalization to improve the APC implementation and
*                      handle the "simple" case properly rather than with
*                      just peak data.
*     June, 2002   "   Added T argument to avoid use of integer abscissas for
*                      APC, UPPERC.
*     Feb., 2003   "   Trapped the case where the curves being compared don't 
*                      overlap at all (or barely overlap, as in the APC case,
*                      where V may increase initially): no good solution.
*
*     Author:  David Saunders, ELORET/NASA Ames Research Center.
*
********************************************************************************

      IMPLICIT NONE

*     Arguments:

      LOGICAL, INTENT (IN) ::
     >   ABOVE              ! T means violations are above the constraint;
                            ! F means below, as in the trajectory APC case
      INTEGER, INTENT (IN) ::
     >   N                  ! Length of time histories

      REAL,    INTENT (IN) ::
     >   X(N), Y(N)         ! Data curve being constrained

      REAL,    INTENT (IN) ::
     >   T(N)               ! Corresponding time steps for cases where
                            ! X(*) /= T(*), as in the case of APC and UPPERC
      INTEGER, INTENT (IN) ::
     >   NC                 ! Length of constraint curve;
                            ! NC = 1 means a simple upper bound, YC
      REAL,    INTENT (IN) ::
     >   XC(NC), YC(NC)     ! Constraint curve

      REAL,    INTENT (OUT) ::
     >   VALUE              ! Measure of constraint violation (or lack of it)
                            ! as described above
      REAL,    INTENT (OUT) ::
     >   XPEAK              ! Abscissa of peak data point violation or least
                            ! non-violation (informational only - no attempt
                            ! is made to interpolate)
*     Local constants:

      INTEGER,   PARAMETER :: LUNERR = 6 ! Avoid yet another argument
      REAL,      PARAMETER :: BIG = 1.E+20, ONE = 1., ZERO = 0.
      LOGICAL,   PARAMETER :: NEW = .TRUE.
      CHARACTER, PARAMETER :: LOOSE * 1 = 'B', MONO * 1 = 'M'

*     Local variables:

      INTEGER    I, I1, I2, J1, J2, NPTS
      REAL       ARROW, DISTANCE_MAX, DISTANCE_MIN, HEIGHT, XCNC
      REAL       DERIVATIVES(N), YDIFF(N)

*     Procedures:

      EXTERNAL   LCSFIT, ! Local cubic spline utility
     >           LCSQUAD ! ... and its quadrature form

*     Execution:

      I2 = N

      IF (NC > 1) THEN

*        Trap the case of (essentially) no overlap at all:

         XCNC  = XC(NC)
         ARROW = SIGN (ONE, XCNC - XC(1))

         IF (ARROW * X(1) > ARROW * XCNC) THEN
            WRITE (LUNERR, '(/, A, 1P, 2E15.5)')
     >         ' CURVE_VIOLATIONS: No overlap of data and constraint: ',
     >         X(1), XCNC, ' Turn this constraint off.  Aborting.'
            STOP
         END IF

*        Find the last data point where a curve comparison is possible:

         IF (ARROW * X(N) > ARROW * XCNC) THEN

            I2 = I2 - 1 ! A smart search should be faster than brute force;
                        ! still a tiny chance of trouble if X(*) isn't monotonic

            CALL INTERVAL (N, X, XCNC, ARROW, I2) ! I2 <= N - 1

            IF (ABS (X(I2+1) - XCNC) < TINY (XCNC)) I2 = I2 + 1

         END IF

*        Interpolate the constraint Y at each X in-range:

         CALL LCSFIT (NC, XC, YC, NEW, LOOSE, I2, X, YDIFF,
     >                DERIVATIVES)

*        Set up a difference curve, Ydata - Yconstraint:

         DO I = 1, I2
            YDIFF(I) = Y(I) - YDIFF(I)
         END DO

      ELSE ! Simple upper bound case: set up differences

         DO I = 1, I2
            YDIFF(I) = Y(I) - YC(1)
         END DO

      END IF

      I1 = 1

      IF (ABOVE) THEN

         IF (NC > 1) THEN

*           An initial portion above the constraint doesn't count until the
*           data curve drops below the constraint OR until the data curve
*           starts increasing while still above the constraint.
*           However, we first eliminate data for which Y (altitude) is
*           increasing at the very beginning (possible for an ascent abort):

            DO I = 1, I2 - 1
               IF (Y(I+1) < Y(I)) THEN
                  I1 = I
                  EXIT
               END IF
            END DO

            IF (YDIFF(I1) > ZERO) THEN
               DO I = I1 + 1, I2
                  IF (YDIFF(I) <= ZERO .OR. Y(I) > Y(I-1)) THEN
                     I1 = I
                     EXIT
                  END IF
               END DO
            END IF

         END IF

*        Locate the worst violation (or the smallest non-violation):

         DISTANCE_MAX = -BIG
         DO I = I1, I2
            HEIGHT = YDIFF(I)
            IF (HEIGHT > DISTANCE_MAX) THEN
               DISTANCE_MAX = HEIGHT
               XPEAK = X(I)
            END IF
         END DO

*        Integrate area patches above the constraint curve w.r.t. T, not X:

         IF (DISTANCE_MAX > ZERO) THEN

            DO I = I1, I2
               YDIFF(I) = MAX (YDIFF(I), ZERO)
            END DO

*           Avoid the bulk of the zeros:

            J1 = I1
            DO I = I1, I2
               IF (YDIFF(I) > ZERO) THEN
                  J1 = MAX (I - 1, I1)
                  EXIT
               END IF
            END DO

            J2 = I2
            DO I = I2, J1, -1
               IF (YDIFF(I) > ZERO) THEN
                  J2 = MIN (I + 1, I2)
                  EXIT
               END IF
            END DO

            IF (J2 <= J1 + 1) THEN
               J1 = MAX (J1 - 1, I1)
               J2 = MIN (J2 + 1, I2)
            END IF

            NPTS = J2 - J1 + 1

            CALL LCSQUAD (NPTS, T(J1), YDIFF(J1), T(J1), T(J2), MONO,
     >                    VALUE)

         ELSE ! No violation

            VALUE = DISTANCE_MAX ! <= 0., to provide some gradient info.

         END IF

      ELSE ! "Below" case

*        Locate the worst violation (or the smallest non-violation):

         DISTANCE_MIN = BIG
         DO I = 1, I2
            HEIGHT = YDIFF(I)
            IF (HEIGHT < DISTANCE_MIN) THEN
               DISTANCE_MIN = HEIGHT
               XPEAK = X(I)
            END IF
         END DO

*        Integrate area patches below the constraint curve w.r.t. T, not X:

         IF (DISTANCE_MIN < ZERO) THEN

            DO I = 1, I2
               YDIFF(I) = MIN (YDIFF(I), ZERO)
            END DO

*           Avoid the bulk of the zeros:

            J1 = I1
            DO I = I1, I2
               IF (YDIFF(I) < ZERO) THEN
                  J1 = MAX (I - 1, I1)
                  EXIT
               END IF
            END DO

            J2 = I2
            DO I = I2, J1, -1
               IF (YDIFF(I) < ZERO) THEN
                  J2 = MIN (I + 1, I2)
                  EXIT
               END IF
            END DO

            IF (J2 <= J1 + 1) THEN
               J1 = MAX (J1 - 1, I1)
               J2 = MIN (J2 + 1, I2)
            END IF

            NPTS = J2 - J1 + 1

            CALL LCSQUAD (NPTS, T(J1), YDIFF(J1), T(J1), T(J2), MONO,
     >                    VALUE)

         ELSE ! No violation

            VALUE = DISTANCE_MIN ! >= 0., to provide some gradient info.

         END IF

      END IF ! End "below" case

      END SUBROUTINE CURVE_VIOLATIONS

********************************************************************************
*
      SUBROUTINE DYNAMIC ()
*
*     Allocate (most of) the workspace that is dependent upon the number of
*     optimization variables and constraints.
*
********************************************************************************

*     Global variables:

      USE OPT

*     Execution:

      ALLOCATE
     >  (IERROR(NDV),   IGROUP(NDV),   ILCON(NCLIN),  INLCON(NCNLN),
     >   JNLCON(NCNLN), ISTATE(NBNDS), NFDEL(NDV),    IWVEC(LENIW),
     >   CW(LENCW))

      ALLOCATE
     >  (FFWD(NDV),     GRAD(NDV),     PWORK(NDV),    SH(NDV),
     >   URH(NDV),      YWORK(NDV),    ZWORK(NDV),    AITCH(NDV),
     >   BL(NBNDS),     BU(NBNDS),     CVEC(NROWJ),   WVEC(LENW),
     >   BLLIN(NCLIN),  BULIN(NCLIN),  CONLIN(NCLIN), TLCON(NCLIN),
     >   C(NCNLN),      CPRINT(NCNLN), XNLCON(NCNLN), SNLCON(NCNLN),
     >   CLAMDA(NBNDS),
     >   OLDG(NDV),     PSRCH(NDV),    DEE(NDV),      ELL(1),
     >   AMAT(NCLIN,NDV), RMAT(NDV,NDV),
     >   CJAC(NROWJ,NDV), CJACBAK(NROWJ))

      ALLOCATE
     >  (LFAIL(NDV))
 
      ALLOCATE
     >  (VTYPE(NDV), LCTYPE(NCLIN), NLCTYPE(NCNLN))

      END SUBROUTINE DYNAMIC

************************************************************************
*
      SUBROUTINE NLCON (PRINT, CPAR)
*
*     NLCON evaluates all nonlinear constraints, CPAR(*), with the option
*     to tabulate UNscaled results.
*
*     5/96  DAS  Initial implementation for monitoring cabin radii.
*     6/96   "   Adapted for use by CONFUN.
*     5/00   "   Adaptation for trajectory optimization.
*     6/00   "   Added APC constraint.
*     8/00   "   Added ending latitude, longitude, & altitude constraints.
*     9/00   "   Added distance from target (latitude, longitude) constraint.
*     4/01   "   Added apoapsis and orbit inclination, for aero-capture.
*     5/01   "   Added apoapsis eccentricity, for aero-capture.
*     7/01   "   Added peak altitude, for max. down-range difficulties;
*                added distributed heating constraints assuming the
*                aerothermal database is processed in Traj (wrong).
*     8/01   "   Added UPPERC and OSCILL constraints; revised the
*                distributed heating constraints.
*     9/01   "   Added DESCNT constraint.
*    12/01   "   Introduced CURVE_VIOLATIONS for APC, upper corridor,
*                DYN_PR, and (new) Q * Alpha constraints.
*     4/02   "   Added TRIM constraint.
*     5/02   "   Added ENDFPA and ENDVEL constraints.
*     6/02   "   Added ENDALP, ENDBNK, and ENDHED constraint.
*
************************************************************************

*     Global variables:

      USE CONSTS
      USE OPT

      IMPLICIT NONE

*     Arguments:

      LOGICAL, INTENT (IN)  :: PRINT   ! T means tabulate UNscaled constraints
      REAL,    INTENT (OUT) :: CPAR(*) ! Nonlinear constraints (= C(*) in OPT)

*     Local constants:

      INTEGER, PARAMETER ::
     >   IV = 1, IH = 2, IT = 3,       ! Column numbers of TRAJ_STEPS
     >   IM = 4, IA = 5, IQ = 6,
     >   IG = 7, IP = 8

      LOGICAL, PARAMETER ::
     >   FALSE = .FALSE., TRUE = .TRUE.

      REAL, PARAMETER ::
     >   PA2PSF = 0.020885472, ! psf per pascal
     >   ZERO   = 0.

*     Local variables:

      INTEGER   I, IWK, J, K, M
      REAL      BOUNDL, BOUNDU, HEIGHT, SCALE, TIME, VALUE
      REAL      T_SIMPLE_BOUND(1), Y_SIMPLE_BOUND(1)
      CHARACTER FLAG*3

*     Execution:

      IF (PRINT) WRITE (LWRIT, 1010)

      K = NDV + NCLIN ! Offset for the bounds

      DO M = 1, NCNLN

         K = K + 1

         SELECT CASE (NLCTYPE(M))

            CASE ('CRANGE')          ! Cross range

               VALUE = CROSS_RANGE
               TIME  = TOTAL_TIME

            CASE ('DRANGE')          ! Down range

               VALUE = DOWN_RANGE
               TIME  = TOTAL_TIME

            CASE ('H_LOAD')          ! Stagnation pt. heat load

               VALUE = HEAT_LOAD
               TIME  = TOTAL_TIME

            CASE ('H_FLUX')          ! Peak stagnation pt. heat flux (total)

               VALUE = MAX_VALUES(2)
               TIME  = MAX_TIMES(2)

            CASE ('W_TEMP')          ! Peak stagnation pt. wall temperature

               VALUE = MAX_VALUES(6)
               TIME  = MAX_TIMES(6)

            CASE ('ACCEL ')          ! Acceleration magnitude

!!!               VALUE = MAX_VALUES(7)
!!!               TIME  = MAX_TIMES(7)

               T_SIMPLE_BOUND(1) = ZERO ! Array avoids compiler warnings
               Y_SIMPLE_BOUND(1) = XNLCON(M) ! Gs

               CALL CURVE_VIOLATIONS (TRUE, LAST_STEP,
     >                               TRAJ_STEPS(1,IT), TRAJ_STEPS(1,IG),
     >                               TRAJ_STEPS(1,IT),
     >                               1, T_SIMPLE_BOUND, Y_SIMPLE_BOUND,
     >                               VALUE, TIME)

            CASE ('DYN_PR')          ! Dynamic pressure

!!!               VALUE = MAX_VALUES(8)
!!!               TIME  = MAX_TIMES(8)

               T_SIMPLE_BOUND(1) = ZERO ! Array avoids compiler warnings
               Y_SIMPLE_BOUND(1) = XNLCON(M) ! Pascals

               CALL CURVE_VIOLATIONS (TRUE, LAST_STEP,
     >                               TRAJ_STEPS(1,IT), TRAJ_STEPS(1,IQ),
     >                               TRAJ_STEPS(1,IT),
     >                               1, T_SIMPLE_BOUND, Y_SIMPLE_BOUND,
     >                               VALUE, TIME)

            CASE ('QALPHA')          ! Q * Alpha (psf-degrees)

               T_SIMPLE_BOUND(1) = ZERO
               Y_SIMPLE_BOUND(1) = XNLCON(M) ! psf-degrees
               IWK = NJOURNAL + 1

               DO I = 1, LAST_STEP
                  TRAJ_STEPS(I,IWK) = TRAJ_STEPS(I,IQ) *
     >                                TRAJ_STEPS(I,IA) * PA2PSF
               END DO

               CALL CURVE_VIOLATIONS (TRUE, LAST_STEP,
     >                               TRAJ_STEPS(1,IT),TRAJ_STEPS(1,IWK),
     >                               TRAJ_STEPS(1,IT),
     >                               1, T_SIMPLE_BOUND, Y_SIMPLE_BOUND,
     >                               VALUE, TIME)

            CASE ('STG_PR')          ! Peak static pressure

               VALUE = MAX_VALUES(9)
               TIME  = MAX_TIMES(9)

            CASE ('PK_ALT')          ! Peak altitude

               T_SIMPLE_BOUND(1) = ZERO ! Array avoids compiler warnings
               Y_SIMPLE_BOUND(1) = XNLCON(M) ! km

               CALL CURVE_VIOLATIONS (TRUE, LAST_STEP,
     >                               TRAJ_STEPS(1,IT), TRAJ_STEPS(1,IH),
     >                               TRAJ_STEPS(1,IT),
     >                               1, T_SIMPLE_BOUND, Y_SIMPLE_BOUND,
     >                               VALUE, TIME)

            CASE ('APC   ')          ! Aerothermal performance (lower corridor)

               CALL CURVE_VIOLATIONS (FALSE, LAST_STEP,
     >                               TRAJ_STEPS(1,IV), TRAJ_STEPS(1,IH),
     >                               TRAJ_STEPS(1,IT),
     >                               N_APC, V_APC, H_APC, VALUE, TIME)

            CASE ('UPPERC')          ! Upper corridor boundary

               CALL CURVE_VIOLATIONS (TRUE, LAST_STEP,
     >                               TRAJ_STEPS(1,IV), TRAJ_STEPS(1,IH),
     >                               TRAJ_STEPS(1,IT),
     >                               N_UPPER_C, V_UPPER_C, H_UPPER_C,
     >                               VALUE, TIME)

            CASE ('OSCILL')          ! Penalize altitude oscillations

               CALL OSCILLATIONS (VALUE, TIME)

            CASE ('ENDALP')          ! Ending angle of attack (deg)

               VALUE = END_ALPHA
               TIME  = TOTAL_TIME

            CASE ('ENDALT')          ! Ending altitude

               VALUE = END_ALTITUDE
               TIME  = TOTAL_TIME

            CASE ('ENDBNK')          ! Ending bank angle (deg)

               VALUE = END_BANK
               TIME  = TOTAL_TIME

            CASE ('ENDFPA')          ! Ending flight path angle, gamma (deg)

               VALUE = END_FPA
               TIME  = TOTAL_TIME

            CASE ('ENDHED')          ! Ending heading, psi (deg)

               VALUE = END_HEADING
               TIME  = TOTAL_TIME

            CASE ('ENDLAT')          ! Ending latitude

               VALUE = END_LATITUDE
               TIME  = TOTAL_TIME

            CASE ('ENDLON')          ! Ending longitude

               VALUE = END_LONGITUDE
               TIME  = TOTAL_TIME

            CASE ('ENDIST')          ! Ending distance from target (deg)

               VALUE = SQRT ((END_LATITUDE  - TARGET_LATITUDE ) ** 2 +
     >                       (END_LONGITUDE - TARGET_LONGITUDE) ** 2)
               ! The SQRT avoids distance ** 4 in the merit function
               TIME  = TOTAL_TIME

            CASE ('ENDVEL')          ! Ending velocity, relative (km/s)

               VALUE = END_VELOCITY
               TIME  = TOTAL_TIME

            CASE ('APOAPS')          ! Aero-capture apoapsis (km from CG)

               ! Apoapsis = semi-latus rectum / (1 - eccentricity)
               VALUE = ORBITAL(2) / (1000.* (1.- ORBITAL(4)))
               TIME  = TOTAL_TIME

            CASE ('INCLIN')          ! Ending orbit inclination (deg)

               VALUE = ORBITAL(5) * RAD2DEG
               TIME  = TOTAL_TIME

            CASE ('ECCENT')          ! Ending orbit eccentricity

               VALUE = ORBITAL(4)
               TIME  = TOTAL_TIME

            CASE ('SURF1 ')          ! Surface heating (1): temperature

               CALL SURFACE_HEATING (1, VALUE)

               TIME = 999.           ! Not defined

            CASE ('SURF1R')          ! Surface heating (1R): T from heat RATE

               CALL SURFACE_HEATING (-1, VALUE) ! Overwrite temperatures
               CALL SURFACE_HEATING ( 1, VALUE) ! As for SURF1

               TIME = 999.           ! Not defined

            CASE ('SURF2 ')          ! Surface heating (2): Peak heat flux

               CALL SURFACE_HEATING (2, VALUE)

               TIME = 999.           ! Not defined

            CASE ('HEATLD')          ! Integrated surface heat load

               VALUE = HEATLD
               TIME  = 999.          ! Not defined

            CASE ('DESCNT')          ! Final descent rate (km/sec)

               I = LAST_STEP - 3     ! The very last two steps may be too close

               CALL LCSFIT (I, TRAJ_STEPS(I,IT), TRAJ_STEPS(I,IH),
     >                      TRUE,  'M',  1,  TRAJ_STEPS(I+2,IT),
     >                      HEIGHT, VALUE)
               VALUE = -VALUE        ! Derivative = ascent rate
               TIME  = TOTAL_TIME    ! Not defined

            CASE ('TRIM  ')          ! Integrated square of Cm or flap angle

               T_SIMPLE_BOUND(1) = ZERO
               Y_SIMPLE_BOUND(1) = TARGET_TRIM**2 ! Tolerable |Cm| or |defln.|
               IWK = NJOURNAL + 1

               DO I = 1, LAST_STEP
                  TRAJ_STEPS(I,IWK) = TRAJ_STEPS(I,IP) ** 2
               END DO

               CALL CURVE_VIOLATIONS (TRUE, LAST_STEP,
     >                               TRAJ_STEPS(1,IT),TRAJ_STEPS(1,IWK),
     >                               TRAJ_STEPS(1,IT),
     >                               1, T_SIMPLE_BOUND, Y_SIMPLE_BOUND,
     >                               VALUE, TIME)

            CASE DEFAULT

               WRITE (LWRIT, '(/, 2A)')
     >            ' NLCON: Bad constraint type = ', NLCTYPE(M)
               STOP

         END SELECT

         SCALE   = SNLCON(M)
         CPAR(M) = VALUE * SCALE

         IF (PRINT) THEN

            BOUNDL = BL(K) / SCALE ! Since the main program applies SNLCON(M)
            BOUNDU = BU(K) / SCALE

            IF (VALUE < BOUNDL .OR. VALUE > BOUNDU) THEN
               FLAG = '***'
            ELSE
               FLAG = '   '
            END IF

            WRITE (LWRIT, 1030) M, NLCTYPE(M), INLCON(M), VALUE,
     >                          BOUNDL, BOUNDU, TIME, FLAG

          END IF

      END DO

      IF (DISTRIBUTED_HEATING) THEN
         IF (PRINT) THEN
            WRITE (LWRIT, 1050)
            DO J = 1, N_SURF_POINTS
               WRITE (LWRIT,
     >            '(I5, 1P, E17.8, E12.4, 0P, 2(2X, F15.7, F12.3))') J,
     >            SURF_HEAT_LOADS(J), SURF_BOUNDS(0,J),
     >           (SURF_PEAKS(I,J),    SURF_BOUNDS(I,J),
     >            I = 1, N_SURF_TYPES)
            END DO
         END IF
      END IF

      RETURN

 1010 FORMAT (/, ' NLCON: Nonlinear constraint check:', //,
     >   '    # NLCTYPE  INLCON  CURRENT VALUE',
     >   '    LOWER BOUND    UPPER BOUND  TIME or VELOCITY')
 1030 FORMAT (1X, I4, 2X, A6, I6, 2X, 3F15.6, F15.3, 2X, A3)
 1050 FORMAT (/, ' Surface heating summary:', //, ' PT #',
     >        '  HT LOAD, J/CM^2      WEIGHT',
     >        '     PEAK TEMP, K       BOUND',
     >        '  PK QDOT, W/CM^2       BOUND')

      END SUBROUTINE NLCON

************************************************************************
*
      SUBROUTINE OBJECTIVE (OBJ)
*
*     OBJECTIVE is called by both OBJFUN_AERO and SOLVE.  It keeps the
*     trajectory-related objective function details in one place.
*
*     04/00  D.A.Saunders  Initial separation from SOLVE.
*     08/00      "         Added target latitude, longitude, and
*                          altitude optimization; added TABORT.
*     08/01      "         Added distributed heating calculations.
*     03/03      "         Added target trajectory option.
*
************************************************************************

*     Global variables:

      USE OPT

      IMPLICIT NONE

*     Arguments:

      REAL, INTENT (OUT) :: OBJ

*     Local constants:

      INTEGER, PARAMETER ::
     >   IAL = 2,  ! Altitude index in time history from Traj
     >   IT  = 3,  ! Time index
     >   IP  = 8,  ! Pitch (Cm/flap angle) index
     >   ILA = 9,  ! Latitude index
     >   ILO = 10  ! Longitude

      REAL,    PARAMETER ::
     >   ONE = 1., NINETIETH = ONE / 90., ZERO = 0.

*     Local variables:

      INTEGER
     >   I, IWK
      REAL
     >   PENA, PENB, PENC, PEND, PENE, PENF, PENG, PENH, PENI, PENJ,
     >   PENK, PENL, PENM, PENN, PENO, PENP, PENQ, PENR
      REAL
     >   TIME, VALUE
      REAL
     >   T_SIMPLE_BOUND(1), Y_SIMPLE_BOUND(1)

*     Execution:

      END_ALTITUDE  = FINAL_STATE(0)
      END_VELOCITY  = FINAL_STATE(1)
      END_LONGITUDE = FINAL_STATE(2)
      END_FPA       = FINAL_STATE(3)
      END_LATITUDE  = FINAL_STATE(4)
      END_HEADING   = FINAL_STATE(5)
      END_ALPHA     = FINAL_STATE(6)
      END_BANK      = FINAL_STATE(7)

      PENA = RHO_CROSS_RANGE
      IF (PENA /= ZERO)
     >    PENA = PENA * (ONE - NINETIETH * ABS (CROSS_RANGE))
******IF (PENA /= ZERO) PENA = PENA * (-CROSS_RANGE)

      PENB = RHO_DOWN_RANGE
      IF (PENB /= ZERO) PENB = PENB * (-DOWN_RANGE)

      PENC = RHO_HEAT_LOAD
      IF (PENC /= ZERO) PENC = PENC * HEAT_LOAD

*     Penalizing only values above a specified bound would make more sense
*     for some of these options (but we don't have inputs for the bounds):

      PEND = RHO_MAX_H_FLUX
      IF (PEND /= ZERO) PEND = PEND * MAX_VALUES(2)

      PENE = RHO_MAX_W_TEMP
      IF (PENE /= ZERO) PENE = PENE * MAX_VALUES(6)

      PENF = RHO_MAX_ACCEL
      IF (PENF /= ZERO) PENF = PENF * MAX_VALUES(7)

      PENG = RHO_MAX_DYN_PR
      IF (PENG /= ZERO) PENG = PENG * MAX_VALUES(8)

      PENH = RHO_TARGET_POINT
      IF (PENH /= ZERO) PENH = PENH *
     >   ((END_LATITUDE  - TARGET_LATITUDE ) ** 2 +
     >    (END_LONGITUDE - TARGET_LONGITUDE) ** 2)

      PENI = RHO_END_ALTITUDE
      IF (PENI /= ZERO) PENI = PENI * (ONE - 0.1 * END_ALTITUDE)
                               ! Goes to zero at 10 km; negative is OK
      PENJ = RHO_ECCENTRICITY
      IF (PENJ /= ZERO) PENJ = PENJ * ORBITAL(4)

      IF (DISTRIBUTED_HEATING) THEN

*        Initialize by interpolating the aerothermal database for all surface
*        quantities at every step in the time history.

         CALL SURFACE_HEATING (0, HEATLD) ! 0 = initialize; HEATLD is unused
         CALL SURFACE_HEATING (9, HEATLD) ! 9 = integrated heat load

      END IF

      PENK = RHO_HEATLD
      IF (PENK /= ZERO) PENK = PENK * HEATLD

      PENL = RHO_TABORT
      IF (PENL /= ZERO) PENL = PENL * TABORT

      PENM = RHO_TRIM

      IF (PENM /= ZERO) THEN

         T_SIMPLE_BOUND(1) = ZERO
         Y_SIMPLE_BOUND(1) = TARGET_TRIM ** 2 ! Tolerable |Cm| or |deflection|
         IWK = NJOURNAL + 1

         DO I = 1, LAST_STEP
            TRAJ_STEPS(I,IWK) = TRAJ_STEPS(I,IP) ** 2
         END DO

         CALL CURVE_VIOLATIONS (.TRUE., LAST_STEP,
     >                          TRAJ_STEPS(1,IT),TRAJ_STEPS(1,IWK),
     >                          TRAJ_STEPS(1,IT),
     >                          1, T_SIMPLE_BOUND, Y_SIMPLE_BOUND,
     >                          VALUE, TIME)
         PENM = PENM * VALUE

      END IF

      PENN = RHO_ALPHA_TVD ! Equivalenced to RHO_CL_TVD

      IF (PENN /= ZERO) THEN ! Contribution tending to smooth Alpha vs. time

*        For a gradient estimate, and f(xi) = (xi+1 - xi) ** 2, we have
*        g(xi) = (f(xi + h) - f(xi)) / h = 2(xi+1 - xi) + O(h), so the
*        gradient for PENN is RHO * 2dX / T(nAlpha).  Furthermore, if
*        T(nAlpha) = 1000 and the average change in Alpha is 5 degrees
*        and nAlpha = 50, then we have
*        PENN = RHO * 50 * 25 / 1000 = RHO * 1.25, suggesting RHO should
*        be about 1.E-m where m is the decimal place intended to be affected.

         VALUE = ZERO

         IF (MODE_VAR /= 6) THEN

            DO I = 1, NK_ALPHA - 1
               VALUE = (ALPHA_K(I+1) - ALPHA_K(I)) ** 2 + VALUE
            END DO

            PENN = PENN * (VALUE / TK_ALPHA(NK_ALPHA))

         ELSE ! MODE_VAR = 6

            DO I = 1, NK_CL - 1
               VALUE = (CL_K(I+1) - CL_K(I)) ** 2 + VALUE
            END DO

            PENN = PENN * (VALUE / TK_CL(NK_CL))

         END IF

      END IF

      PENO = RHO_BANK_TVD ! or RHO_CD_TVD

      IF (PENO /= ZERO) THEN

         IF (MODE_VAR /= 6) THEN

            VALUE = ZERO
            DO I = 1, NK_BANK - 1
               VALUE = (BANK_K(I+1) - BANK_K(I)) ** 2 + VALUE
            END DO

            PENO = PENO * (VALUE / TK_BANK(NK_BANK))

         ELSE ! MODE_VAR = 6

            DO I = 1, NK_CD - 1
               VALUE = (CD_K(I+1) - CD_K(I)) ** 2 + VALUE
            END DO

            PENO = PENO * (VALUE / TK_CD(NK_CD))

         END IF

      END IF

      PENP = RHO_BETA_TVD ! or RHO_CY_TVD

      IF (PENP /= ZERO) THEN

         IF (MODE_VAR /= 6) THEN

            VALUE = ZERO
            DO I = 1, NK_BETA - 1
               VALUE = (BETA_K(I+1) - BETA_K(I)) ** 2 + VALUE
            END DO

            PENP = PENP * (VALUE / TK_BETA(NK_BETA))

         ELSE ! MODE_VAR = 6

            DO I = 1, NK_CY - 1
               VALUE = (CY_K(I+1) - CY_K(I)) ** 2 + VALUE
            END DO

            PENP = PENP * (VALUE / TK_CY(NK_CY))

         END IF

      END IF

      PENQ = RHO_DURATION

      IF (PENQ /= ZERO) PENQ = PENQ * TOTAL_TIME

      PENR = RHO_LST_SQRS

      IF (PENR /= ZERO) THEN ! Target trajectory option

         CALL COMPARE_TRAJECTORIES (LAST_STEP,
     >                TRAJ_STEPS(1,IT),  TRAJ_STEPS(1,ILA),
     >                TRAJ_STEPS(1,ILO), TRAJ_STEPS(1,IAL),
     >                NORMALIZE, N_TARGET,
     >                TRAJ_TARGET(1,1), TRAJ_TARGET(1,2),
     >                TRAJ_TARGET(1,3), TRAJ_TARGET(1,4),
     >                VALUE)
         PENR = PENR * VALUE

      END IF

      OBJ = PENA + PENB + PENC + PEND + PENE + PENF + PENG + PENH +
     >      PENI + PENJ + PENK + PENL + PENM + PENN + PENO + PENP +
     >      PENQ + PENR

      END SUBROUTINE OBJECTIVE

************************************************************************
*
      SUBROUTINE OBJFUN (MODE, NDVPAR, VPAR, OBJPAR, GRDPAR, NSTATE)
*
*     OBJFUN acts as the interface between NPOPT and SOLVE. Its main task
*     is to determine whether NPOPT wants just the objective function or
*     the objective function plus gradient information, and then call
*     SOLVE appropriately.
*
*        03/96  J. Reuther  Initial implementation (SYN87).
*        02/00  D. Saunders Initial adaptation for trajectory optimization.
*     05/15/00     "        NNGRAD = 2 or 3 allows CONFUN to do the work.
*     07/11/00     "        NNDIFF = 3 allows for central differencing.
*
************************************************************************

*     Global variables:

      USE CONSTS
      USE OPT

      IMPLICIT NONE

*     Arguments:

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

*     Local constants:

      CHARACTER, PARAMETER :: SUBNAME * 6 = 'OBJFUN'

*     Local variables:

      INTEGER I, NCALL, NHALF
      REAL    OBJFOR, OBJBAK, TEMPH, TEMPV, TIME1, TIME2
      LOGICAL FAIL

*     Execution:

      IF (CONFUN_DOES_OBJ) THEN ! Don't call SOLVE at all

         OBJPAR = OBJ_CONFUN

         IF (MODE > 0) THEN

            GRDPAR = FFWD ! Elements 1 : NDV;  CONFUN now prints the gradient

**          WRITE (LWRIT, 1030)
**    >         'OBJFUN: Objective gradient via CONFUN:', ' ',
**    >         '   I            V (I)            G (I)            H (I)'
**           WRITE (LWRIT, '(I5, 1P, 3E17.8)')
**    >         (I, VPAR(I), GRDPAR(I), AITCH(I), I = 1, NDVPAR)
         END IF

         GO TO 999

      END IF

*     Conventional case: no nonlinear constraints, or above option suppressed.

***      WRITE (LWRIT, 1040) 'MODE and NSTATE: ', MODE, NSTATE

      IF (NSTATE == 1) THEN ! This is the first OBJFUN call from NPOPT

         NCALL = -3

         CALL SOLVE (NDVPAR, VPAR, OBJPAR, NCALL, FAIL)

         NFTOTL = NFTOTL + 1

         IF (FAIL) THEN
            WRITE (LWRIT, 1040) 'Bad return from first call to SOLVE.'
            GO TO 900
         END IF

         IF (NITMAX == 0) THEN
            WRITE (LWRIT, 1040) 'Forcing exit from NPOPT (NITMAX = 0).'
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
               WRITE (LWRIT, 1040)'NPOPT may shorten the step & retry.'
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

         IF (NNDIFF == 2) THEN ! Forward differencing

            DO NCALL = 1, NDVPAR

               TEMPV = VPAR(NCALL)
               TEMPH = AITCH(NCALL)
               NHALF = 0

   10          CONTINUE

               VPAR(NCALL) = TEMPV + TEMPH

               CALL SOLVE (NDVPAR, VPAR, OBJFOR, NCALL, FAIL)

               NFTOTL = NFTOTL + 1

               IF (FAIL) THEN
                  WRITE (LWRIT, 1040)
     >               'Function failure for gradient element ', NCALL

                  IF  (NHALF < 3) THEN
                     NHALF = NHALF + 1
                     TEMPH = TEMPH * 0.5
                     WRITE (LWRIT, 1030)
     >                  'Halve the stepsize for this element.'
                     GO TO 10
                  END IF

                  OBJFOR = OBJPAR
                  WRITE (LWRIT, 1030)
     >               'Three halvings of gradient step size failed.',
     >               'Proceeding with zero for this component.'
               END IF

               GRDPAR(NCALL) = (OBJFOR - OBJPAR) / TEMPH
               VPAR(NCALL)   = TEMPV

            END DO

         ELSE ! NNDIFF == 3: central differencing

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
         CALL CPUTIME (TIME1, TIME2, 'OBJFUN',
     >      'to evaluate the objective gradient', LWRIT)

         WRITE (LWRIT, 1030) 'OBJFUN: Gradient vector:', ' ',
     >      '   I            V (I)            G (I)            H (I)'
         WRITE (LWRIT, '(I5, 1P, 3E17.8)')
     >      (I, VPAR(I), GRDPAR(I), AITCH(I), I = 1, NDVPAR)

*        Save restart information:   ! (Seems redundant for trajectories)

***      CALL USER (2, NDV, 1, VPAR, OBJPAR, GRDPAR, AITCH,
***  >              ELL, DEE, OLDG, PSRCH, GTP, GNM, NFTOTL, NALLOW,
***  >              NCALL, NDITER, NTYPE, .FALSE., .TRUE., LWRIT)
      END IF

      GO TO 999

  900 MODE = -1 ! NPOPT should quit

  999 RETURN

 1030 FORMAT (/, (1X, A))
 1040 FORMAT (/, ' OBJFUN: ', A, 2I4)

      END SUBROUTINE OBJFUN

************************************************************************
*
      SUBROUTINE OBJFUN_AERO (NDVPAR, VPAR, OBJPAR, NCALL, FAIL)
*
*     OBJFUN_AERO is a specialized objective function routine in CENDIF2/
*     QNMDIF2 form, for estimating trajectory sensitivities to the aero-
*     dynamic coefficients (CL & CD vs. Mach, Alpha, Re) in the Traj aero
*     database.
*
*     NCALL usage:
*
*        -3 = initial function evaluation
*        >0 = gradient element NCALL evaluation
*
*     4/00  D.A.Saunders  Stripped-down form of SOLVE for the aero
*                         database sensitivities case.
*
************************************************************************

*     Global variables:

      USE CONSTS
      USE OPT

      IMPLICIT NONE

*     Arguments:

      INTEGER, INTENT (IN)    :: NDVPAR       ! Number of design variables
      INTEGER, INTENT (INOUT) :: NCALL        ! See usage above
      REAL,    INTENT (OUT)   :: OBJPAR       ! Objective function value
      REAL,    INTENT (INOUT) :: VPAR(NDVPAR) ! Variables
      LOGICAL, INTENT (OUT)   :: FAIL         ! T means the function is bad

*     Local constants:

      INTEGER, PARAMETER ::
     >   ICALL = 1        ! Normal trajectory calculation

      CHARACTER, PARAMETER ::
     >   SUBNAME*11 = 'OBJFUN_AERO'

*     Local variables:

      INTEGER
     >   I, IRETURN

      REAL, SAVE ::
     >   CD_ORIG, CL_ORIG, CPU1, CPU2

      LOGICAL
     >   FIRST

*     Procedures:

      EXTERNAL
     >   CPUTIME, OBJECTIVE, SECOND, TRAJECTORY

*     Execution:
*     ----------

      FIRST = NCALL == -3

      IF (FIRST) THEN

         WRITE (LWRIT, '(//, 1X, 2A)') SUBNAME,
     >      ': Trajectory sensitivities with respect to aero database.'
         WRITE (LWRIT, '(/, 3(A, I4))')
     >      ' N_Mach:', N_MACH, '  N_Alpha:', N_ALPHA, '  N_Re:', N_RE
         CALL SECOND (CPU1)

      ELSE ! Update an aerodynamic coefficient

         IF (NCALL <= N_COEFS) THEN
            CL_ORIG = CL(NCALL)
            CL(NCALL) = VPAR(NCALL)
         ELSE
            I = NCALL - N_COEFS
            CD_ORIG = CD(I)
            CD(I) = VPAR(NCALL)
         END IF

      END IF

*     Calculate a trajectory:
*     -----------------------

      CALL TRAJECTORY (ICALL, NARGS, ARGS, MODE_VAR,
     >                 NNTIME,   N_ENTRY,    ENTRY_STATE,
     >                 N_MAX,    MAX_VALUES, MAX_TIMES,
     >                 NK_BANK,  TK_BANK,    BANK_K,
     >                 NK_ALPHA, TK_ALPHA,   ALPHA_K,
     >                 N_MACH, MACH, N_ALPHA, ALPHA, N_RE, RE,
     >                 CL, CD,   PERCENT_LOVERD,
     >                 FINAL_STATE, TOTAL_TIME,
     >                 CROSS_RANGE, DOWN_RANGE, HEAT_LOAD, ORBITAL,
     >                 NJOURNAL, MAX_STEP, LAST_STEP, TRAJ_STEPS,
     >                 IRETURN)

      FAIL = IRETURN == 0

      IF (FAIL) THEN
         WRITE (LWRIT, '(/, 1X, 2A, I6)')
     >      SUBNAME, ': Traj return code:', IRETURN
         GO TO 999
      END IF

      IF (FIRST) THEN

         CALL SECOND (CPU2)
         CALL CPUTIME (CPU1, CPU2, SUBNAME,
     >                 'to compute the trajectory', LWRIT)

      ELSE ! Restore the perturbed aero. coefficient

         IF (NCALL <= N_COEFS) THEN
            CL(NCALL) = CL_ORIG 
         ELSE
            CD(NCALL - N_COEFS) = CD_ORIG
         END IF

      END IF

*     Objective function:
*     -------------------

      CALL OBJECTIVE (OBJPAR)

      WRITE (LWRIT, '(/, A, I6, A, 1P, E25.15)')
     >   ' NCALL:', NCALL, '    Objective:', OBJPAR

      RETURN

  999 STOP

      END SUBROUTINE OBJFUN_AERO

************************************************************************
*
      SUBROUTINE OSCILLATIONS (VALUE, TIME)
*
*     Calculate a measure of how much the current altitude history
*     oscillates.  The lower bound for this nonlinear constraint should
*     be zero to enforce a steadily descending trajectory.
*
*     08/03/01  DAS  Initial implementation, for down-range maximization.
*
************************************************************************

*     Global variables:

      USE CONSTS
      USE OPT

      IMPLICIT NONE

*     Arguments:

      REAL, INTENT (OUT) :: VALUE  ! Sum of altitude increases with time
      REAL, INTENT (OUT) :: TIME   ! Time step of pt. most (nearly) in violation

*     Local constants:

      INTEGER, PARAMETER :: IH = 2 ! Column number of altitude in TRAJ_STEPS
      REAL,    PARAMETER :: ZERO = 0.

*     Local variables:

      INTEGER  I, I_MAX
      REAL     DHEIGHT, DHEIGHT_MAX, SUM

*     Execution:

      SUM = ZERO
      DHEIGHT_MAX = -1.E10

      DO I = 2, LAST_STEP
         DHEIGHT = TRAJ_STEPS(I,IH) - TRAJ_STEPS(I-1,IH)
         IF (DHEIGHT > DHEIGHT_MAX) THEN
            DHEIGHT_MAX = DHEIGHT
            I_MAX = I
         END IF
         IF (DHEIGHT > ZERO) SUM = DHEIGHT + SUM
      END DO

      IF (SUM > ZERO) THEN
         VALUE = SUM
      ELSE
         VALUE = DHEIGHT_MAX
      END IF

      TIME = REAL (I_MAX)

      END SUBROUTINE OSCILLATIONS

********************************************************************************
*
      SUBROUTINE OUTPUT (OBJPAR, NCALL)
*
*     Save plottable results at each step (supplements Traj outputs).
*
********************************************************************************

*     Global variables:

      USE CONSTS
      USE OPT

      IMPLICIT NONE

*     Arguments:

      REAL,    INTENT (IN) :: OBJPAR
      INTEGER, INTENT (IN) :: NCALL

*     Local variables:

      INTEGER  I

*     Execution:

      REWIND (LRESULTS)
      WRITE (LRESULTS, '(3A)') 'Traj_opt Trajectory Optimization  (',
     >   VEHICLE(1:LEN_VEHICLE), ')'

      IF (MODE_VAR <= 4) THEN
         WRITE (LRESULTS, '(A, I4, 3X, A, I4)')
     >      'Alpha variables:', NDV_ALPHA, 'Bank variables:', NDV_BANK
         WRITE (LRESULTS, '(A)')
     >      'Time (seconds)', 'Angle (degrees)', '""'
      ELSE IF (MODE_VAR == 5) THEN
         WRITE (LRESULTS, '(A, I4, 3X, A, I4, 3X, A, I4)')
     >      'Alpha variables:', NDV_ALPHA, 'Bank variables:', NDV_BANK,
     >      'Beta variables:', NDV_BETA
         WRITE (LRESULTS, '(A)')
     >      'Time (seconds)', 'Angle (degrees)', '""'
      ELSE !   MODE_VAR == 6
         WRITE (LRESULTS, '(A, I4, 3X, A, I4, 3X, A, I4)')
     >      'CL variables:', NDV_CL, 'CD variables:', NDV_CD,
     >      'CY variables:', NDV_CY
         WRITE (LRESULTS, '(A)') 'Time (seconds)', 'Coefficient', '""'
      END IF

      WRITE (LRESULTS, '(A, I4, A, 1P, E15.6)')
     >   'At optimization iteration', NDITER,
     >   ',  objective function  =', OBJPAR, '""'

      WRITE (LRESULTS, '(2A)')
     >   ' $options xnumbers=''whole'', ynumbers=''whole'',',
     >   ' grid=1, $end'

      IF (MODE_VAR <= 5) THEN
         IF (NDV_ALPHA > 0) THEN
            WRITE (LRESULTS, '(2A)')
     >         ' $options legend=''Alpha knots'', fit=''tight'',',
     >         ' spline=''standard'', symbol=0, $end'
            WRITE (LRESULTS, '(A, /, A, I4)')
     >         '!Alpha knots', '!', NK_ALPHA
            WRITE (LRESULTS, '(1P, 2E16.6)')
     >         (TK_ALPHA(I), ALPHA_K(I), I = 1, NK_ALPHA)
         END IF

         IF (NDV_BANK > 0) THEN
            WRITE (LRESULTS, '(2A)')
     >         ' $options legend=''Bank knots'', fit=''tight'',',
     >         ' spline=''standard'', symbol=1, $end'
            WRITE (LRESULTS, '(A, /, A, I4)')
     >         '!Bank knots', '!', NK_BANK
            WRITE (LRESULTS, '(1P, 2E16.6)')
     >         (TK_BANK(I), BANK_K(I), I = 1, NK_BANK)
         END IF
      END IF

      IF (MODE_VAR == 5) THEN
         WRITE (LRESULTS, '(2A)')
     >      ' $options legend=''Beta knots'', fit=''tight'',',
     >      ' spline=''standard'', symbol=2, $end'
         WRITE (LRESULTS, '(A, /, A, I4)')
     >      '!Beta knots', '!', NK_BETA
         WRITE (LRESULTS, '(1P, 2E16.6)')
     >      (TK_BETA(I), BETA_K(I), I = 1, NK_BETA)
      ELSE IF (MODE_VAR == 6) THEN
         WRITE (LRESULTS, '(2A)')
     >      ' $options legend=''CL knots'', fit=''tight'',',
     >      ' spline=''standard'', symbol=0, $end'
         WRITE (LRESULTS, '(A, /, A, I4)') '!CL knots', '!', NK_CL
         WRITE (LRESULTS, '(1P, 2E16.6)')
     >      (TK_CL(I), CL_K(I), I = 1, NK_CL)
         WRITE (LRESULTS, '(2A)')
     >      ' $options legend=''CD knots'', fit=''tight'',',
     >      ' spline=''standard'', symbol=1, $end'
         WRITE (LRESULTS, '(A, /, A, I4)') '!CD knots', '!', NK_CD
         WRITE (LRESULTS, '(1P, 2E16.6)')
     >      (TK_CD(I), CD_K(I), I = 1, NK_CD)
         WRITE (LRESULTS, '(2A)')
     >      ' $options legend=''CY knots'', fit=''tight'',',
     >      ' spline=''standard'', symbol=2, $end'
         WRITE (LRESULTS, '(A, /, A, I4)') '!CY knots', '!', NK_CY
         WRITE (LRESULTS, '(1P, 2E16.6)')
     >      (TK_CY(I), CY_K(I), I = 1, NK_CY)
      END IF

      END SUBROUTINE OUTPUT

********************************************************************************
*
      SUBROUTINE PERTURB (VPAR, ICALL)
*
*     This is all that is left of the aero. shape optimization routine for
*     perturbing geometry.  Here, it perturbs trajectory alpha and/or bank
*     schedules and possibly further quantities such as TABORT.
*     See the main program for descriptions of the inputs.
*
********************************************************************************

*     Global variables:

      USE CONSTS
      USE OPT

      IMPLICIT NONE

*     Arguments:

      REAL,    INTENT (IN) ::
     >   VPAR(*)            ! Current optimization variables

      INTEGER, INTENT (IN) ::
     >   ICALL              ! Used by step function mode:
                            ! 1 = set up alpha and bank in spline form;
                            ! 2 = as for 1 plus write them in compact form
*     Local constants:

      REAL,      PARAMETER :: ONE = 1.
      LOGICAL,   PARAMETER :: NEW = .TRUE.
      CHARACTER, PARAMETER :: METHOD * 1 = 'M'

*     Local variables:

      INTEGER    L, M, N
      REAL       D

*     Execution:
*     ----------

      IF (MODE_VAR <= 2) THEN ! Continuously-variable alpha and/or bank

         N = NDV_ANGLES

         IF (NDV_ALPHA > 0) THEN

            ALPHA_K (KNOT1_A : NK_ALPHA) = VPAR (1 : NDV_ALPHA) *
     >                                    VSCALE(1)
         END IF

         IF (NDV_BANK > 0) THEN

            BANK_K (KNOT1_B : NK_BANK) = VPAR (1 + NDV_ALPHA : N) *
     >                                  VSCALE(1 + NDV_ALPHA)
         END IF

         DO N = NDV_ANGLES + 1, NDV

            SELECT CASE (VTYPE(N))

               CASE ('TABORT') ! Time of ascent abort

                  TABORT = VPAR(N) * VSCALE(N)
               
                  CALL INTERVAL (N_ASCENT, ASCENT_TIME, TABORT, ONE,
     >                           LABORT) ! LABORT is initialized in RDVARS

                  L = MAX (LABORT - 1, 1) ! Narrow the data range for LCSFIT
                  M = MIN (4, N_ASCENT - L + 1) ! # pts. for 4-pt. method

                  CALL LCSFIT (M, ASCENT_TIME(L), ASCENT_ALTITUDE(L),
     >                        NEW, METHOD, 1, TABORT, ENTRY_STATE(0), D)

                  CALL LCSFIT (M, ASCENT_TIME(L), ASCENT_VELOCITY(L),
     >                        NEW, METHOD, 1, TABORT, ENTRY_STATE(1), D)
                  ENTRY_STATE(1) = ENTRY_DV + ENTRY_STATE(1)

                  CALL LCSFIT (M, ASCENT_TIME(L), ASCENT_LONGITUDE(L),
     >                        NEW, METHOD, 1, TABORT, ENTRY_STATE(2), D)

                  CALL LCSFIT (M, ASCENT_TIME(L), ASCENT_GAMMA(L),
     >                        NEW, METHOD, 1, TABORT, ENTRY_STATE(3), D)

                  CALL LCSFIT (M, ASCENT_TIME(L), ASCENT_LATITUDE(L),
     >                        NEW, METHOD, 1, TABORT, ENTRY_STATE(4), D)

                  CALL LCSFIT (M, ASCENT_TIME(L), ASCENT_HEADING(L),
     >                        NEW, METHOD, 1, TABORT, ENTRY_STATE(5), D)

               CASE DEFAULT

                  WRITE (LWRIT, '(/, 2A)')
     >               ' PERTURB: Bad variable type = ', VTYPE(N)
                  STOP

            END SELECT

         END DO

      ELSE IF (MODE_VAR <= 4) THEN ! Step function mode

         CALL ANGLE_STEPS (ICALL, MODE_VAR, NSTEPS,
     >                     ALPHA_STEPS, BANK_STEPS,
     >                     DURATION, PITCH_RATE, ROLL_RATE,
     >                     ALPHA_K, BANK_K, TK_ALPHA, TK_BANK,
     >                     NDV, VPAR, VSCALE, LDAT, LWRIT)

      ELSE IF (MODE_VAR == 5) THEN ! Alpha, bank, Beta are variables

         ALPHA_K(1:NK_ALPHA) = VPAR(1:NDV_ALPHA) * VSCALE(1)

         BANK_K (1:NK_BANK)  = VPAR(1+NDV_ALPHA : NDV_ALPHA + NDV_BANK)
     >                                          * VSCALE(1 + NDV_ALPHA)
         BETA_K (1:NK_BETA)  = VPAR(1+NDV_ALPHA + NDV_BANK : NDV_ANGLES)
     >                                              * VSCALE(NDV_ANGLES)
      ELSE ! MODE_VAR == 6; CL, CD, CY are variables

         CL_K (1 : NK_CL) = VPAR (1 : NDV_CL) * VSCALE(1)

         CD_K (1 : NK_CD) = VPAR (1 + NDV_CL : NDV_CL + NDV_CD)
     >                                     * VSCALE(1 + NDV_CL)
         CY_K (1 : NK_CY) = VPAR (1 + NDV_CL + NDV_CD : NDV_COEFFS)
     >                                         * VSCALE(NDV_COEFFS)
      END IF

      END SUBROUTINE PERTURB

************************************************************************
*
      SUBROUTINE RDVARS (LUNRD, LUNERR, LDAT)
*
*     File-driven control scheme for trajectory optimization variables
*     other than the alpha and bank schedules.  The variables read here
*     become V(NDV_ANGLES + 1), ...., V(NDV).
*
*     See the main program for descriptions.  Variable order is important.
*     For instance, if ascent abort time is being optimized, then this
*     is assumed to be V(NDV_ANGLES + 1).
*
*     08/30/00  David Saunders  Eliminated sine bumps, etc.; provided
*                               for trajectory variables supplementary
*                               to the alpha and bank control points.
*
************************************************************************

*     Global variables:

      USE OPT

      IMPLICIT NONE

*     Arguments:

      INTEGER, INTENT (IN) :: LUNRD, LUNERR, LDAT

*     Local variables:

      INTEGER  I, N
      LOGICAL  ORDERED

*     Execution:

      OPTIMIZE_TABORT = .FALSE.

      DO I = NDV_ANGLES + 1, NDV

         READ (LUNRD, *)
     >      N, VTYPE(I), V(I), VSCALE(I), AITCH(I), BL(I), BU(I)

         IF      (VTYPE(I) == 'TABORT') THEN
            IGROUP(I) = 1
            OPTIMIZE_TABORT = .TRUE.

!!!      ELSE IF (VTYPE(I) == 'VENTRY') THEN  ! Entry conditions some day
!!!         IGROUP(I) = 2
         ELSE
            WRITE (LUNERR, '(/, A, I5, 2X, A)')
     >         ' RDVARS: Unknown design variable:', I, VTYPE(I)
            GO TO 900
         END IF

      END DO

      ORDERED = .TRUE.
      DO I = NDV_ANGLES + 2, NDV
         IF (IGROUP(I) < IGROUP(I-1)) ORDERED = .FALSE.
      END DO

      IF (.NOT. ORDERED) THEN
         WRITE (LUNERR, '(/, 1X, A, A, /, (4X, A))')
     >      'RDVARS:  Supplementary variables must be in this order:',
     >      'TABORT', 'VENTRY, etc.'
         GO TO 900
      END IF

      IF (OPTIMIZE_TABORT) THEN ! Read the ascent trajectory data

         OPEN (LDAT, FILE='traj_opt.ascent', STATUS='OLD')

         READ (LDAT, *) ! Title
         READ (LDAT, *) N

         N_ASCENT = N
         LABORT   = N / 2 ! Pointer into ASCENT_TIME for TABORT; see PERTURB

         ALLOCATE
     >     (ASCENT_TIME(N),     ASCENT_VELOCITY(N),
     >      ASCENT_GAMMA(N),    ASCENT_HEADING(N),
     >      ASCENT_ALTITUDE(N), ASCENT_LATITUDE(N),
     >      ASCENT_LONGITUDE(N))

         READ (LDAT, *)
     >     (ASCENT_TIME(I),     ASCENT_VELOCITY(I),
     >      ASCENT_GAMMA(I),    ASCENT_HEADING(I),
     >      ASCENT_ALTITUDE(I), ASCENT_LATITUDE(I),
     >      ASCENT_LONGITUDE(I), I = 1, N)

         CLOSE (LDAT)

      END IF

      RETURN

  900 STOP

      END SUBROUTINE RDVARS

C+---------------------------------------------------------------------
C
      SUBROUTINE RVERSE ( NX, X1, X2 )
C
C  PURPOSE:  RVERSE reverses the order of the elements of X1(*), re-
C            turning the result in X2(*). In-place reordering is OK.
C
C  ARGUMENTS:
C  ARG  DIM   TYPE I/O/S DESCRIPTION
C  NX   -       I    I   No. of elements in given array (odd or even).
C  X1   NX      R    I   Array whose ordering is to be reversed.
C  X2   NX      R    O   Reordered version of X1(*).  The calling pro-
C                        gram may pass the same array for X1 & X2.
C
C  ENVIRONMENT:  FORTRAN IV with slight updates on 08/02/01.
C
C  AUTHOR: David Saunders, Informatics, Palo Alto, CA.  (07/30/82)
C
C----------------------------------------------------------------------

      REAL, INTENT (INOUT) :: X1(NX), X2(NX)

      INTEGER I, II, NXBY2
      REAL    XI

      NXBY2 = (NX + 1) / 2

      DO I = 1, NXBY2
         XI     = X1(I)
         II     = (NX + 1) - I
         X2(I)  = X1(II)
         X2(II) = XI
      END DO

      END SUBROUTINE RVERSE

********************************************************************************
*
      SUBROUTINE SETLCON (BIGBND)
*
*     SETLCON sets up the linear constraint matrix.
*
*     3/96  J.Reuther  Initial SYN87 implementation for T/C constraints
*                      (lower bounds only).
*     6/96     "       Upper bounds provided for. WARNING: The linear constraint
*                      set-up must be done with XINIT/YINIT (the geometry being
*                      perturbed). Applying any initial nonzero variables then
*                      setting up with XPTRB/YPTRB is wrong.
*     8/96  D.Saunders SFEVAL avoids interpolation of unit perturbations at
*                      constraint Xs.
*     02/00    "       Adaptation for trajectory optimization: no constraints.
*  05/09/00    "       Implemented linear constraints on bank and alpha.
*  08/30/00    "       Constraints on bank and alpha perturbations are no longer
*                      relevant - took them out.
*  09/24/01    "       Pitch rate and roll rate constraints implemented.
*
********************************************************************************

*     Global variables:

      USE CONSTS
      USE OPT

      IMPLICIT NONE

*     Arguments:

      REAL, INTENT (IN) :: BIGBND ! Large number for no-constraint cases

*     Local constants:

      REAL, PARAMETER :: ZERO = 0.

*     Local variables:

      INTEGER L, M, N

*     Execution:

*     Zero the linear constraint matrix - only some of the variables
*     are applicable.

      AMAT = ZERO ! (1:NCLIN,1:NDV)

*     Determine the nonzero matrix elements.
*     An internal procedure treats pitch and bank similarly:

      L = NDV ! Offset in BL(*) and BU(*) for linear constraint bounds
      M = 0   ! Row
      N = 0   ! Column

      IF (LINCON == 1 .OR. LINCON == 3) THEN
         IF (NDV_ALPHA > 0) THEN
            CALL RATE_CONSTRAINTS (KNOT1_A, NK_ALPHA, TK_ALPHA)
         END IF
      END IF

      N = NDV_ALPHA ! Column

      IF (LINCON >= 2) THEN
         IF (NDV_BANK > 0) THEN
            CALL RATE_CONSTRAINTS (KNOT1_B, NK_BANK, TK_BANK)
         END IF
      END IF

*     Internal procedure for SETLCON:

      CONTAINS

         SUBROUTINE RATE_CONSTRAINTS (KNOT1, NKNOTS, TIMES)

!        Set up linear constraints on either pitch rate or roll rate

!        Arguments:

         INTEGER, INTENT (IN) ::
     >      KNOT1,             ! First alpha or bank control point varying
     >      NKNOTS             ! Number of control points
         REAL,    INTENT (IN) ::
     >      TIMES(NKNOTS)      ! Associated times

!        Local constants:

         REAL, PARAMETER ::
     >      ONE = 1.

!        Local variables:

         INTEGER
     >      K
         REAL
     >      DT

!        Execution:

         DO K = KNOT1, NKNOTS - 1
            DT = ONE / (TIMES(K+1) - TIMES(K))

            M = M + 1
            N = N + 1
            AMAT(M,N)   = -DT
            AMAT(M,N+1) =  DT

            L = L + 1
            BL(L) = BLLIN(M)
            BU(L) = BULIN(M)
         END DO

         END SUBROUTINE RATE_CONSTRAINTS

      END SUBROUTINE SETLCON

********************************************************************************
*
      SUBROUTINE SOLVE (NDVPAR, VPAR, OBJPAR, NCALL, FAIL)
*
*     SOLVE is the trajectory-based objective function routine expected
*     by the unconstrained optimizer QNMDIF2.  It is also called by the
*     OBJFUN and CONFUN routines expected by the constrained optimizer NPOPT.
*
*     NCALL usage:
*
*       -13 = first call to SOLVE from CONFUN (precedes first call to OBJFUN)
*       -12 = constraint evaluations following a line search
*       -10 = constraint evaluations within a line search
*        -4 = final call to SOLVE to recover best trajectory & plot results
*        -3 = first call to SOLVE from OBJFUN (or Traj_opt if QNMDIF2 mode)
*        -2 = objective evaluation after a line search
*        -1 = local search evaluation (QNMDIF2 mode only)
*         0 = objective evaluation within a line search
*        >0 = gradient element NCALL evaluation (constraints, obj.fn., or both)
*
*     1986-  Robert Kennelly     FLO6QNM version.
*     1990+  James Reuther       OPT67 version.
*     8/94+  David Saunders/JJR  Overhaul of common variables, etc.
*     4/95+    "      "      "   SYN87 adaptation.
*     9/98+  DAS                 Fortran 90 translation.
*     5/00    "                  Trajectory optimization version of CONFUN
*                                must call SOLVE, not just PERTURB.
* 05/15/00    "                  NNGRAD = 2 or 3 case suppresses trajectory
*                                & objective function calculations that will
*                                have been done by CONFUN:  SOLVE is not
*                                even called in these cases now if NCNLN > 0.
*                                NCALL <= 0 case for CONFUN now has 10
*                                subtracted from the OBJFUN equivalent.
*                                NCALL usage has become very awkward!
* 08/22/01    "                  Added distributed heating calculations.
*
********************************************************************************

*     Global variables:

      USE CONSTS
      USE OPT

      IMPLICIT NONE

*     Arguments:

      INTEGER, INTENT (IN)    :: NDVPAR  ! Number of design variables
      INTEGER, INTENT (INOUT) :: NCALL   ! See usage above
      REAL,    INTENT (OUT)   :: OBJPAR  ! Objective function value
      REAL,    INTENT (INOUT) :: VPAR(*) ! * because NDVPAR may be zero.
      LOGICAL, INTENT (OUT)   :: FAIL    ! T means the function value is bad

*     Local constants:

      REAL, PARAMETER ::
     >   HALF = 0.5, ONE = 1., ZERO = 0.

      LOGICAL, PARAMETER ::
     >   PRINT = .TRUE.

      CHARACTER, PARAMETER ::
     >   PROGNAME*8 = 'TRAJ_OPT', SUBNAME*5 = 'SOLVE'

*     Local variables:

      INTEGER
     >   I, ICALL, IRETURN, MCALL
      REAL, SAVE ::
     >   CPU1, CPU2, TIME1, TIME2
      LOGICAL
     >   FIRST, OUTPUTS

*     Procedures:

      EXTERNAL
     >   CHKLCON, CPUTIME, FLUSH, NLCON, OUTPUT, OBJECTIVE, PERTURB,
     >   SECOND, TRAJECTORY

************************************************************************
*
*     Clarifications:
*
*     (1) CONFUN_DOES_OBJ = T means OBJFUN does not call SOLVE at all.
*     (2) CONFUN_DOES_OBJ = T means CONFUN calls must do the outputs that
*         would otherwise be done by calls from OBJFUN (or QNMDIF2/main).
*     (3) NCALL = -13, -12, -10  correspond to CONFUN calls equivalent to
*         NCALL =  -3,  -2,   0  calls from OBJFUN.

************************************************************************

*     Execution:

      CALL SECOND (TIME1)

      IF (CONFUN_DOES_OBJ) THEN ! For the first objective evaluation ...
         FIRST = NCALL == -13   ! CONFUN calls to SOLVE serve for OBJFUN
      ELSE
         FIRST = NCALL == -3
      END IF

*     MCALL allows either OBJFUN or CONFUN to trigger outputs.

      MCALL = NCALL
      IF (MCALL < 0) MCALL = MOD (MCALL, 10)
      OUTPUTS = MCALL <= 0    ! If both CONFUN & OBJFUN, CONFUN shouldn't output
      IF (.NOT. CONFUN_DOES_OBJ) OUTPUTS = OUTPUTS .AND. MCALL == NCALL

      IF (MCALL == -2. AND. OUTPUTS) THEN ! A line search has completed
         NDITER = NDITER + 1
      END IF

      IF (OUTPUTS) THEN
         WRITE (LWRIT, '(/, 1X, 2A, I4, A, I5)') SUBNAME,
     >      ': Optimization iteration:', NDITER, '      NCALL:', NCALL
         IF (MCALL <= -2) THEN
            WRITE (LWRIT, '(/, A, /)') ' OPTIMIZATION VARIABLES:'
            WRITE (LWRIT, '(5(I5, 1P, E21.13))')
     >         (I, VPAR(I), I = 1, NDVPAR)
         END IF
      END IF

      IF (NCALL == -4) THEN      ! Trajectory + plot & clean-up
         ICALL = 2
      ELSE IF (NITMAX == 0) THEN ! Activate the plots here too
         ICALL = 2
      ELSE                       ! Trajectory only
         ICALL = 1
      END IF

*     Apply the optimization variables:
*     ---------------------------------

      CALL PERTURB (VPAR, ICALL) ! ICALL is used in step function mode

*     Calculate a trajectory:
*     -----------------------

      IF (MCALL <= -2) CALL SECOND (CPU1)

      IF (MODE_VAR <= 4) THEN ! All zero-side-force cases

         CALL TRAJECTORY (ICALL, NARGS, ARGS, MODE_VAR,
     >                    NNTIME,   N_ENTRY,    ENTRY_STATE,
     >                    N_MAX,    MAX_VALUES, MAX_TIMES,
     >                    NK_BANK,  TK_BANK,    BANK_K,
     >                    NK_ALPHA, TK_ALPHA,   ALPHA_K,
     >                    N_MACH, MACH, N_ALPHA, ALPHA, N_RE, RE,
     >                    CL, CD,   PERCENT_LOVERD,
     >                    FINAL_STATE, TOTAL_TIME,
     >                    CROSS_RANGE, DOWN_RANGE, HEAT_LOAD, ORBITAL,
     >                    NJOURNAL, MAX_STEP, LAST_STEP, TRAJ_STEPS,
     >                    IRETURN)

      ELSE IF (MODE_VAR == 5) THEN ! Alpha, bank, and beta are variables

         CALL TRAJ2 (ICALL, NARGS, ARGS, ENTRY_STATE,
     >               MAX_VALUES, MAX_TIMES,
     >               NK_BANK,  TK_BANK,  BANK_K,
     >               NK_ALPHA, TK_ALPHA, ALPHA_K,
     >               NK_BETA,  TK_BETA,  BETA_K,
     >               FINAL_STATE, TOTAL_TIME,
     >               CROSS_RANGE, DOWN_RANGE, HEAT_LOAD,
     >               MAX_STEP, LAST_STEP, TRAJ_STEPS,
     >               IRETURN)

      ELSE IF (MODE_VAR == 6) THEN ! CL, CD, and CY are variables

         CALL TRAJ3 (ICALL, NARGS, ARGS, ENTRY_STATE,
     >               MAX_VALUES, MAX_TIMES,
     >               NK_CL,  TK_CL,  CL_K,
     >               NK_CD,  TK_CD,  CD_K,
     >               NK_CY,  TK_CY,  CY_K,
     >               FINAL_STATE, TOTAL_TIME,
     >               CROSS_RANGE, DOWN_RANGE, HEAT_LOAD,
     >               MAX_STEP, LAST_STEP, TRAJ_STEPS,
     >               IRETURN)

      END IF

      FAIL = IRETURN == 0

      IF (FAIL) THEN
         WRITE (LWRIT, '(/, 1X, A, (A, I4))')
     >      SUBNAME, ': Traj return code:', IRETURN,
     >      ' ICALL input to TRAJECTORY:', ICALL,
     >      ' NCALL input to SOLVE:', NCALL
         GO TO 999
      END IF

      IF (MCALL <= -2) THEN

         CALL SECOND (CPU2)
         CALL CPUTIME (CPU1, CPU2, SUBNAME,
     >                 'to compute the trajectory', LWRIT)
      END IF

*     Objective function:
*     -------------------

      IF (MCALL <= 0 .OR. NNGRAD == 0) THEN ! Main prog. does MOD (NNGRAD, 2)

*        Line search, or finite-difference mode:

         CALL OBJECTIVE (OBJPAR)

         IF (OUTPUTS) THEN

            IF (SHOW_SHORTFALL) THEN
               SHORTFALL = 60. *
     >            SQRT ((END_LATITUDE  - TARGET_LATITUDE ) ** 2 +
     >                  (END_LONGITUDE - TARGET_LONGITUDE) ** 2)
            END IF
            IF (SHOW_APOAPS) APOAPSIS    = ORBITAL(2) /
     >                       (1000.* (1. - ORBITAL(4)))
            IF (SHOW_ECCENT) ORBITAL_ECC = ORBITAL(4)
            IF (SHOW_INCLIN) ORBITAL_INC = ORBITAL(5) * RAD2DEG

            IF (FIRST) THEN
               OBJ_INI           = OBJPAR
               CROSS_RANGE_INI   = CROSS_RANGE
               DOWN_RANGE_INI    = DOWN_RANGE
               END_ALPHA_INI     = END_ALPHA
               END_ALTITUDE_INI  = END_ALTITUDE
               END_BANK_INI      = END_BANK
               END_FPA_INI       = END_FPA
               END_HEADING_INI   = END_HEADING
               END_LATITUDE_INI  = END_LATITUDE
               END_LONGITUDE_INI = END_LONGITUDE
               END_VELOCITY_INI  = END_VELOCITY
               HEAT_LOAD_INI     = HEAT_LOAD
               TOTAL_TIME_INI    = TOTAL_TIME
               MAX_VALUES_INI    = MAX_VALUES ! 0 : N_MAX - 1
               MAX_TIMES_INI     = MAX_TIMES  ! "     "     "
               IF (OPTIMIZE_TABORT) TABORT_INI      = TABORT
               IF (SHOW_SHORTFALL)  SHORTFALL_INI   = SHORTFALL
               IF (SHOW_APOAPS)     APOAPSIS_INI    = APOAPSIS
               IF (SHOW_ECCENT)     ORBITAL_ECC_INI = ORBITAL_ECC
               IF (SHOW_INCLIN)     ORBITAL_INC_INI = ORBITAL_INC
               IF (DISTRIBUTED_HEATING)  HEATLD_INI = HEATLD
            END IF

            WRITE (LWRIT, '(/, 1P, A, E22.13, A, E22.13)')
     >         ' OBJECTIVE:', OBJPAR, '   TIME:', TOTAL_TIME

         ELSE IF (MCALL > 0) THEN

            IF (NDITER <= 1) THEN
               IF (MCALL == 1) WRITE (LWRIT, '(A)')
               WRITE (LWRIT,
     >         '(1X, 2A, I4, 3A, 1P, E19.11, A, E21.13, A, 0P, F11.4)')
     >            SUBNAME, ': Variable perturbed:', MCALL,
     >            '   Type: ', VTYPE(MCALL), '  Value:', VPAR(MCALL),
     >            '   Obj:', OBJPAR, '  Time:', TOTAL_TIME
            ELSE
               IF (NDITER == 2) THEN
                  IF (MCALL == 1) WRITE (LWRIT, '(/, A)')
     >               ' Suppressing perturbed variable info.'
               END IF
            END IF

         END IF

      ELSE ! NNGRAD = 1

*        Adjoint mode:

***      CALL ADJOINT ! Unlikely to be implemented

*        These terms represent the variation in the (composite) objective fn.
*        needed for gradient information.  The optimizer still differences
*        function values, so:

***      OBJPAR = OBJMIN + <various terms>

***      WRITE (LWRIT, '(/, 1X, A, 1P, E21.12)') 'OBJPAR:', OBJPAR

***      IF (NDVPAR <= 20 .OR. NCALL <= 2)
***      WRITE (LWRIT, '(1X, A, 1P, E21.12)')
***  >      'OBJMIN:', OBJMIN

      END IF

*     -----------------------------------------------
*     Save/display results at each optimization step:
*     -----------------------------------------------

      IF (OUTPUTS .AND. MCALL <= -2) THEN     ! 1st call or end of line search

         CALL OUTPUT (OBJPAR, NCALL)          ! Save control schedules

         IF (NCLIN > 0) CALL CHKLCON          ! Review the linear constraints

         IF (NCNLN > 0) CALL NLCON (PRINT, CPRINT)  ! & nonlinear constraints

         CALL SUMMARY (OBJPAR)

      END IF

      IF (MCALL == NDVPAR) THEN

         IF (NDITER <= 1) THEN
            CALL SECOND  (TIME2)
            CALL CPUTIME (TIME1, TIME2, SUBNAME,
     >                    'for this call to SOLVE', LWRIT)
         END IF

      ELSE IF (MCALL <= 0) THEN

         CALL SECOND (TIME2)

         IF       (MCALL == -3) THEN
            WRITE (LWRIT, 1080) 'Initial trajectory.'
         ELSE IF  (MCALL == -2) THEN
            WRITE (LWRIT, 1080) 'End-of-line-search evaluation.'
         ELSE IF  (MCALL == -1) THEN
            WRITE (LWRIT, 1080) 'Local search evaluation.'
         ELSE IF  (MCALL ==  0) THEN ! May be -4 at end of run
            WRITE (LWRIT, 1080) 'Line search evaluation.'
         END IF

         IF (NDITER <= 1) THEN
            CALL CPUTIME (TIME1, TIME2, SUBNAME,
     >                    'for this call to SOLVE', LWRIT)
         END IF

         IF (MCALL <= -2) THEN
            CALL CPUTIME (CPU0, TIME2, SUBNAME, 
     >                    'for this run so far', LWRIT)
         END IF

      END IF


*     Check for avoiding a wasteful gradient calculation (specially
*     if finite differencing is requested, or if results are intended
*     for comparing adjoint-based gradients with finite differences):

      IF (MCALL == -2) THEN ! End of line search; NPOPT makes it
                            ! hard to handle the NITMAX > 1 case
         IF (NITMAX == 1) THEN
            WRITE (LWRIT, 1300)
     >         'NITMAX = 1; skipping the second gradient.'
            GO TO 900
         END IF

      END IF

      GO TO 999


*     Termination:

  900 STOP

  999 IF (MCALL <= 0) THEN
         CALL FLUSH (LWRIT) ! Traj_opt output
         CALL FLUSH (6)     ! NPOPT output
      END IF

      RETURN

*     Formats:

 1080 FORMAT (/, ' SOLVE: ', A)
 1300 FORMAT (/, ' SOLVE: ', (1X, A))

      END SUBROUTINE SOLVE

********************************************************************************
*
      SUBROUTINE SUMMARY (OBJ)
*
*     Summarize the initial trajectory and optimized result.
*
********************************************************************************

*     Global variables:

      USE CONSTS
      USE OPT

      IMPLICIT NONE

*     Argument:

      REAL    OBJ ! True current objective; may not be the same as OBJMIN

*     Local variables:

      INTEGER I

*     Execution:

      WRITE (LWRIT, 1200)   NDITER, 
     >   OBJ_INI,           OBJ,
     >   TOTAL_TIME_INI,    TOTAL_TIME,
     >   CROSS_RANGE_INI,   CROSS_RANGE,
     >   DOWN_RANGE_INI,    DOWN_RANGE,
     >   END_ALTITUDE_INI,  END_ALTITUDE,
     >   END_VELOCITY_INI,  END_VELOCITY,
     >   END_FPA_INI,       END_FPA,
     >   END_HEADING_INI,   END_HEADING,
     >   HEAT_LOAD_INI,     HEAT_LOAD,
     >  (MAX_VALUES_INI(I), MAX_VALUES(I),
     >   MAX_TIMES_INI(I),  MAX_TIMES(I), I = 0, N_MAX - 2)

      IF (SHOW_APOAPS)      WRITE (LWRIT, 1220) APOAPSIS_INI, APOAPSIS
      IF (SHOW_ECCENT)      WRITE (LWRIT, 1240) ORBITAL_ECC_INI,
     >                                          ORBITAL_ECC
      IF (SHOW_INCLIN)      WRITE (LWRIT, 1260) ORBITAL_INC_INI,
     >                                          ORBITAL_INC
      IF (SHOW_SHORTFALL)   WRITE (LWRIT, 1300) SHORTFALL_INI, SHORTFALL
      IF (OPTIMIZE_TABORT)  WRITE (LWRIT, 1400) TABORT_INI, TABORT
      IF (DISTRIBUTED_HEATING)
     >                      WRITE (LWRIT, 1500) HEATLD_INI, HEATLD

      WRITE (LWRIT, 1600) END_LATITUDE, END_LONGITUDE

 1200 FORMAT
     >  (//, ' PERFORMANCE PARAMETERS AT OPTIMIZATION ITERATION', I4,
     >   //, 43X, 'Initial', 12X, 'Current',
     >   '     Corresponding Times', 1P, /,
     >   ' Objective function            ', 2E19.10, /,
     >   ' Total time               (sec)', 2E19.10, /,
     >   ' Cross-range (degrees latitude)', 2E19.10, /,
     >   ' Down-range   (surface arc, km)', 2E19.10, /,
     >   ' Terminal altitude         (km)', 2E19.10, /,
     >   '     "    velocity     (km/sec)', 2E19.10, /,
     >   '     "  flight path angle (deg)', 2E19.10, /,
     >   '     "  heading angle     (deg)', 2E19.10, /,
     >   ' Fay-Riddell heat load (J/cm^2)', 2E19.10, /,
     >   ' Convective heat flux    (peak)', 2E19.10, 0P, 2F12.3, 1P, /,
     >   ' Radiative  heat flux  (W/cm^2)', 2E19.10, 0P, 2F12.3, 1P, /,
     >   ' Total      heat flux (stag.pt)', 2E19.10, 0P, 2F12.3, 1P, /,
     >   ' Convective/ablative       "   ', 2E19.10, 0P, 2F12.3, 1P, /,
     >   ' Radiative /ablative       "   ', 2E19.10, 0P, 2F12.3, 1P, /,
     >   ' Total     /ablative       "   ', 2E19.10, 0P, 2F12.3, 1P, /,
     >   ' Peak temperature   (degrees K)', 2E19.10, 0P, 2F12.3, 1P, /,
     >   '   "  acceleration          (G)', 2E19.10, 0P, 2F12.3, 1P, /,
     >   '   "  dynamic pressure (pascal)', 2E19.10, 0P, 2F12.3, 1P, /,
     >   '   "  stagnation pr.       "   ', 2E19.10, 0P, 2F12.3)

 1220 FORMAT  (' Orbit apoapsis    (km from CG)', 2F19.7)
 1240 FORMAT  (' Terminal orbit eccentricity   ', 2F19.7)
 1260 FORMAT  (' Terminal orbit inclin.   (deg)', 2F19.7)
 1300 FORMAT  (' Shortfall from target   (n mi)', 2F19.7)
 1400 FORMAT  (' Ascent abort time        (sec)', 2F19.6)
 1500 FORMAT  (' Integrated heat load  (J/cm^2)', 2F19.2)
 1600 FORMAT (/' Terminal (latitude, longitude)', 2F19.8)

      END SUBROUTINE SUMMARY

********************************************************************************
*
      SUBROUTINE SURFACE_HEATING (ITEM, VALUE)
*
*     Perform distributed surface heating calculations, including initial
*     evaluations for the time history, peak detections, and integrated
*     heat load over the surface.
*
*     The originally-planned peak constraints are implemented as follows:
*
*     Calculate the total peak violation of surface heating limits over all
*     surface points, for the quantity indicated by ITEM.  The surface
*     points may have different upper bounds, but the constraint upper
*     bound should be 0.  If all points are under their limits, the
*     constraint value returned is the sum of the undershoots, in order
*     to provide useful gradient information.
*
*     This version integrates area elements above the bounds, if any,
*     and sums them rather than summing the peak violations.
*
*     07/20/01  D.Saunders  Initial implementation (database in Traj).
*     08/22/01      "       There was no need for Traj to manage the
*                           aerothermal database; Traj_opt does it.
*     09/10/01    DAS/JJR   James suggested areas rather than peaks.
*     09/10/02      DAS     ITEM = -1 case (SURF1R constraint) derives
*                           temperatures from interpolated heat RATEs.
*
********************************************************************************

*     Global variables:

      USE CONSTS
      USE OPT

      IMPLICIT NONE

*     Arguments:

      INTEGER, INTENT (IN)  :: ITEM   ! 0 = initialize the surface heating
                                      !     calculations by interpolating the
                                      !     aerothermal database at each step
                                      !     of the current time history;
                                      ! 1 = evaluate surface temperature
                                      !     constraint;
                                      !-1 = overwrite interpolated Ts with
                                      !     heat-flux-derived values in
                                      !     readiness for a second +1 call;
                                      ! 2 = evaluate surface heat flux constr.;
                                      ! ................. ??? constraint;
                                      ! 9 = evaluate integrated heat load
      REAL,    INTENT (OUT) :: VALUE  ! > 0. means average over-the-limit value
                                      ! < 0. means average under-the-limit value
                                      !      (for 0 < ITEM < 9 in both cases)

*     Local constants:

      INTEGER, PARAMETER ::
     >   IT = 3, IM = 4, IA = 5, IQ = 6 ! Column numbers in TRAJ_STEPS
*        time    Mach    alpha   Qbar     (1 & 2 are for vel. & alt.)

      REAL, PARAMETER ::
     >   ONE = 1., ZERO = 0.

      CHARACTER, PARAMETER ::
     >   METHOD * 1 = 'M'

*     Local variables:

      INTEGER 
     >   I, J, K, I_M, J_A, K_Q
      REAL
     >   ALFA, ARROW, CLOSEST, MNUM, QBAR, P, Q, R, PM1, QM1, RM1,
     >   RSIGEPS, SUM1, VIOLATION

*     Execution:

      I = ITEM

      IF (I == 0) THEN ! Initialize by interpolating the aerothermal database
                       ! for all quantities at all points at each time step

         SURF_PEAKS = ZERO ! All types at all points

         I_M = 1
         J_A = 1
         K_Q = 1

         DO K = 1, LAST_STEP

            MNUM  = TRAJ_STEPS(K,IM)
            ALFA  = TRAJ_STEPS(K,IA)
            QBAR  = TRAJ_STEPS(K,IQ)
            ARROW = SIGN (ONE, ALPHA_DH(2) - ALPHA_DH(1))

*           Locate the relevant cell in the rectangular 3-variable table:

            CALL INTERVAL ( N_MACH_DH,  MACH_DH, MNUM, ONE,   I_M)
            CALL INTERVAL (N_ALPHA_DH, ALPHA_DH, ALFA, ARROW, J_A)
            CALL INTERVAL ( N_QBAR_DH,  QBAR_DH, QBAR, ONE,   K_Q)

            P = MIN (ONE, MAX (ZERO, (MNUM     -  MACH_DH(I_M)) /
     >                        ( MACH_DH(I_M+1) -  MACH_DH(I_M))))
            Q = MIN (ONE, MAX (ZERO, (ALFA     - ALPHA_DH(J_A)) /
     >                        (ALPHA_DH(J_A+1) - ALPHA_DH(J_A))))
            R = MIN (ONE, MAX (ZERO, (QBAR     -  QBAR_DH(K_Q)) /
     >                        ( QBAR_DH(K_Q+1) -  QBAR_DH(K_Q))))
            PM1 = ONE - P
            QM1 = ONE - Q
            RM1 = ONE - R

            DO J = 1, N_SURF_POINTS
               DO I = 1, N_SURF_TYPES
                  SURF_QUANTITIES(K,I,J) =
     >               RM1 * (QM1 * (PM1 * ATDB(I,J,I_M,  J_A,  K_Q)    +
     >                             P   * ATDB(I,J,I_M+1,J_A,  K_Q))   +
     >                      Q   * (PM1 * ATDB(I,J,I_M,  J_A+1,K_Q)    +
     >                             P   * ATDB(I,J,I_M+1,J_A+1,K_Q)))  +
     >               R   * (QM1 * (PM1 * ATDB(I,J,I_M,  J_A,  K_Q+1)  +
     >                             P   * ATDB(I,J,I_M+1,J_A,  K_Q+1)) +
     >                      Q   * (PM1 * ATDB(I,J,I_M,  J_A+1,K_Q+1)  +
     >                             P   * ATDB(I,J,I_M+1,J_A+1,K_Q+1)))
                   SURF_PEAKS(I,J) =
     >                MAX (SURF_QUANTITIES(K,I,J),  SURF_PEAKS(I,J))
                   SURF_VIOLATIONS(K,I,J) =
     >                MAX (SURF_QUANTITIES(K,I,J) - SURF_BOUNDS(I,J),
     >                     ZERO)
               END DO
            END DO
         END DO

      ELSE IF (I == -1) THEN ! Overwrite Ts with Ts from heat flux values;
                             ! emissivities should replace heat flux bounds
         DO J = 1, N_SURF_POINTS
            RSIGEPS = RSTEFAN / SURF_BOUNDS(2,J)
            SURF_PEAKS(1,J) = ZERO
            DO K = 1, LAST_STEP
               SURF_QUANTITIES(K,1,J) =
     >            SQRT (SQRT (RSIGEPS * SURF_QUANTITIES(K,2,J)))
               SURF_PEAKS(1,J) =
     >                MAX (SURF_QUANTITIES(K,1,J),  SURF_PEAKS(1,J))
                   SURF_VIOLATIONS(K,1,J) =
     >                MAX (SURF_QUANTITIES(K,1,J) - SURF_BOUNDS(1,J),
     >                     ZERO)
            END DO
         END DO

!        Now this routine should be called again with ITEM = +1

      ELSE IF (I < 9) THEN ! 1 = distrib. temperature; 2 = distrib. heat rate

         SUM1    = ZERO
         CLOSEST = -1.E+30

         DO J = 1, N_SURF_POINTS
            VIOLATION = SURF_PEAKS(I,J) - SURF_BOUNDS(I,J)

            IF (VIOLATION > ZERO) THEN ! Integrate area elements above bound

               CALL LCSQUAD (LAST_STEP, TRAJ_STEPS(1,IT),
     >                       SURF_VIOLATIONS(1,I,J),
     >                       TRAJ_STEPS(1,IT), TRAJ_STEPS(LAST_STEP,IT),
     >                       METHOD, SURF_HOT_AREAS(J))

               SUM1 = SURF_HOT_AREAS(J) + SUM1
            ELSE
               SURF_HOT_AREAS(J) = ZERO
               CLOSEST = MAX (VIOLATION, CLOSEST) ! Distance under the bound
            END IF

         END DO

         IF (SUM1 > ZERO) THEN
            VALUE = SUM1
         ELSE ! Provide some information if no surface point's limit is exceeded
            VALUE = CLOSEST
         END IF

      ELSE ! ITEM = 9 (area-weighted integrated heat load)

         VALUE = ZERO

         DO J = 1, N_SURF_POINTS ! Integrate the interpolated heat fluxes

            CALL LCSQUAD (LAST_STEP, TRAJ_STEPS(1,IT),
     >                    SURF_QUANTITIES(1,2,J),
     >                    TRAJ_STEPS(1,IT), TRAJ_STEPS(LAST_STEP,IT),
     >                    METHOD, SURF_HEAT_LOADS(J))

            VALUE = VALUE + SURF_BOUNDS(0,J) * SURF_HEAT_LOADS(J)

         END DO

      END IF

      END SUBROUTINE SURFACE_HEATING

************************************************************************
*
      SUBROUTINE USER (IENTRY, NDVPAR, NLPAR, VPAR, OBJPAR, GPAR, HPAR,
     >   ELPAR, DPAR, OLDGPAR, PPAR, GTPPAR, GNMPAR, NFPAR, NFMAX,
     >   NCALL, NITPAR, NTYPAR, CONTINUE, PRINT, LUNLOG)
*
*     Description:
*
*           User-supplied routine called by QNMDIF2, intended to give the
*        user control over printing and/or plotting as the optimization
*        proceeds, to allow the flow field to be updated before the next
*        gradient is calculated, and to handle writing a restart file if
*        requested.  Different applications may require modification to
*        this routine, or to the points in QNMDIF2 where it is called.
*
*           The current best solution is recalculated with full printout
*        after each successful line search (IENTRY = 1), while full restart
*        data is written to disk only after the gradient and Hessian have
*        been brought up to date (IENTRY = 2).
*
*     History:
*
*     03/15/96 D.Saunders All arguments now differ from variables added
*                         to /OPT/ for NPOPT purposes.  This routine is
*                         called when either NPOPT or QNMDIF2 is used.
*     02/26/00   "    "   Initial version for trajectory optimization.
*     01/21/03   "    "   Sensitivities to aero. database needed to be
*                         more easily interpreted.
*
************************************************************************

*     Global variables:

      USE CONSTS
      USE OPT

      IMPLICIT NONE

*     Arguments:

      INTEGER
     >   IENTRY, NDVPAR, NLPAR, NFPAR, NFMAX, NCALL, NITPAR,
     >   NTYPAR, LUNLOG

      REAL
     >   VPAR(NDVPAR), OBJPAR, GPAR(NDVPAR), HPAR(NDVPAR),
     >   ELPAR(NLPAR), DPAR(NDVPAR), OLDGPAR(NDVPAR),
     >   PPAR(NDVPAR), GTPPAR, GNMPAR

      LOGICAL
     >   CONTINUE, PRINT

*     Local constants:

      REAL, PARAMETER :: ZERO = 0.

*     Local variables:

      INTEGER
     >   I, IOS, J, K, M, N, N_MA, NSORT

      INTEGER, ALLOCATABLE ::
     >   NONZERO (:)

      REAL, ALLOCATABLE ::
     >   SENSITIVITY (:)

      LOGICAL
     >   FAIL

*     Execution:

      IF (IENTRY == 1) THEN ! Redo best line search solution

         NCALL = -2

         CALL SOLVE (NDVPAR, VPAR, OBJPAR, NCALL, FAIL)

         IF (FAIL) THEN
            WRITE (LWRIT, '(/,A)')
     >         ' USER: End-of-line-search evaluation failed.'
            STOP
         END IF

         NFPAR = NFPAR + 1
         WRITE (LWRIT, '(/,A,I7)')
     >      ' USER: Total # function evaluations so far:', NFPAR
         IF (NFPAR >= NFMAX) THEN
            WRITE (LWRIT, *) 'Max. # function evaluations reached.'
            STOP
         END IF

      ELSE IF (IENTRY == 2) THEN ! Save optimization info.

***      IF (CONTINUE) THEN
***         REWIND (LSAVE)

***         IF (NNOPT == 1) THEN ! QNMDIF2 restart quantities:

***            WRITE (LSAVE) OBJPAR, NDITER, NTYPAR, GTPPAR, GNMPAR,
***  >                       VPAR, GPAR, ELPAR, DPAR, OLDGPAR, PPAR
***            WRITE (LSAVE) HPAR ! May be skipped when reread

***         ELSE         ! NNOPT = 2; NPOPT restart quantities:

***            WRITE (LSAVE) OBJPAR, NDITER, VPAR, GPAR, OLDGPAR, PPAR
***            WRITE (LSAVE) HPAR
***            WRITE (LSAVE) (BL(I), I = NDV + 1, NDV + NCLIN),
***  >                       (BU(I), I = NDV + 1, NDV + NCLIN),
***  >                       ((AMAT(I,J), I = 1, NCLIN), J = NCHCKS + 1,
***  >                                                  NCHCKS + NCHCK1)
***            WRITE (LSAVE) ISTATE, CLAMDA, RMAT
***         END IF

***         WRITE (LWRIT, '(/,A,I3)')
***  >         ' USER: Optimization restart data saved on unit', LSAVE
***      END IF

      ELSE IF (IENTRY == 3) THEN ! Save aero coef. database sensitivities

         WRITE (LSAVE, '(3I5, A)') N_MACH, N_ALPHA, N_RE,
     >                        ' !  N_MACH, N_ALPHA, N_RE'
         WRITE (LSAVE, '(/, A, 1P, E15.5)')
     >      ' Coefficient perturbation used for central differences:',
     >      AITCH(1)

*        Kludge needed because Traj does not return the Mach, Alpha, and
*        log (Re) values from the aerodynamic database.  Read them here.

         OPEN (LDAT, FILE=AERO_FILENAME, STATUS='OLD', IOSTAT=IOS)

         IF (IOS /= 0) THEN
            WRITE (LWRIT, '(/, 2A)')
     >         ' USER: Unable to open aero. database: ', AERO_FILENAME,
     >         ' Proceeding.'
         ELSE
            DO K = 1, N_RE
               DO J = 1, N_ALPHA
                  DO I = 1, N_MACH
                     READ (LDAT, *) MACH(I), ALPHA(J), RE(K)
                  END DO
               END DO
            END DO
            CLOSE (LDAT)
         END IF

         N_MA = N_MACH * N_ALPHA
         ALLOCATE (NONZERO (N_COEFS), SENSITIVITY (N_COEFS))

         WRITE (LSAVE, '(/, A, /)')
     >      ' Non-zero CL sensitivities, unsorted:'

         NSORT = 0
         DO N = 1, N_COEFS
            IF (GPAR(N) /= ZERO) THEN
               NSORT = NSORT + 1
               SENSITIVITY(NSORT) = GPAR(N)
               NONZERO(NSORT)     = N
               K =  (N - 1) / N_MA + 1
               J = ((N - 1) - (K - 1) * N_MA) / N_MACH + 1
               I =   N - (J - 1) * N_MACH - (K - 1) * N_MA

               WRITE (LSAVE, '(I6, 3I4, 1P, E15.5)')
     >            N, I, J, K, GPAR(N)
            END IF
         END DO

         CALL HSORTRI (SENSITIVITY, NONZERO, NSORT) ! In-place heap sort

         WRITE (LSAVE, '(/, A, //, 2A, /)')
     >      ' Non-zero CL sensitivities, sorted:',
     >      ' ICOEF   I   J   K   d(OBJ)/d(CL)   ',
     >      ' MACH   ALPHA  LOG_RE    CL'

         DO M = 1, NSORT
            N = NONZERO(M)
            K =  (N - 1) / N_MA + 1
            J = ((N - 1) - (K - 1) * N_MA) / N_MACH + 1
            I =   N - (J - 1) * N_MACH - (K - 1) * N_MA

            WRITE (LSAVE, '(I6, 3I4, 1P, E15.5, 0P, 3F8.2, F10.6)')
     >         N, I, J, K, SENSITIVITY(M),
     >         MACH(I), ALPHA(J), RE(K), CL(N)
         END DO

*        Repeat for CDs:

         WRITE (LSAVE, '(/, A, /)')
     >      ' Non-zero CD sensitivities, unsorted:'

         NSORT = 0
         DO N = 1, N_COEFS
            IF (GPAR(N + N_COEFS) /= ZERO) THEN
               NSORT = NSORT + 1
               SENSITIVITY(NSORT) = GPAR(N + N_COEFS)
               NONZERO(NSORT)     = N
               K =  (N - 1) / N_MA + 1
               J = ((N - 1) - (K - 1) * N_MA) / N_MACH + 1
               I =   N - (J - 1) * N_MACH - (K - 1) * N_MA

               WRITE (LSAVE, '(I6, 3I4, 1P, E15.5)')
     >            N, I, J, K, GPAR(N + N_COEFS)
            END IF
         END DO

         CALL HSORTRI (SENSITIVITY, NONZERO, NSORT) ! In-place heap sort

         WRITE (LSAVE, '(/, A, //, 2A, /)')
     >      ' Non-zero CD sensitivities, sorted:',
     >      ' ICOEF   I   J   K   d(OBJ)/d(CD)   ',
     >      ' MACH   ALPHA  LOG_RE    CD'

         DO M = 1, NSORT
            N = NONZERO(M)
            K =  (N - 1) / N_MA + 1
            J = ((N - 1) - (K - 1) * N_MA) / N_MACH + 1
            I =   N - (J - 1) * N_MACH - (K - 1) * N_MA

            WRITE (LSAVE, '(I6, 3I4, 1P, E15.5, 0P, 3F8.2, F10.6)')
     >         N, I, J, K, SENSITIVITY(M),
     >         MACH(I), ALPHA(J), RE(K), CD(N)
         END DO

         DEALLOCATE (NONZERO, SENSITIVITY)

      END IF

      END SUBROUTINE USER
