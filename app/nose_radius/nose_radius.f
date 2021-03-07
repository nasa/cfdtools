****                Program NOSE_RADIUS modules                 ******

************************************************************************
      MODULE CONST_MOD                                     ! Constants !
************************************************************************

      INTEGER, PARAMETER ::
     >   LREAD   = 1,  ! Input control file
     >   LWRIT   = 6,  ! Printable outputs
     >   LNPOPT  = 2,  ! Temporary file used by NPOPT
     >   LNPOUT  = 6,  ! NPOPT output = standard output
     >   LGEOM   = 3,  ! Input airfoil geometry (wrapped clockwise)
     >   LGEOPT  = 4,  ! Optimized airfoil (same format)
     >   LTARGET = 7,  ! Target [radius of] curvature distribution
     >   LRADOPT = 8,  ! Optimized radius of curvature distribution
     >   LCRVOPT = 9   ! Optimized curvature distribution

      REAL
     >   BIGBND, EPSMCH, SMALL

      END MODULE CONST_MOD

************************************************************************
      MODULE DESIGN_VARS_MOD                        ! Design variables !
************************************************************************

      INTEGER
     >   NBUMPS

      INTEGER, ALLOCATABLE, DIMENSION (:) ::
     >   IGROUP, ILCON

      REAL, ALLOCATABLE, DIMENSION (:) ::
     >   V, VSCALE, DOLO, DOUP, POWER, TLCON, XBUMP, XBMIN, XBMAX

      CHARACTER, ALLOCATABLE, DIMENSION (:) ::
     >   XBTYP * 5

      END MODULE DESIGN_VARS_MOD

************************************************************************
      MODULE GEOM_MOD        ! Initial and current geometry definition !
************************************************************************

      INTEGER
     >   NL, NU, NW, NLTARG, NTARGET, NUTARG

      CHARACTER * 5
     >   WBTYP (2)

      REAL, DIMENSION (2) ::
     >   WPOWER, WXCNTR, WXMIN, WXMAX, WPEAK, ADDTO

      REAL, ALLOCATABLE, DIMENSION (:) ::
     >   X_INIT, Y_INIT, S_INIT, CRV_INIT, R_INIT,
     >   X_PRTB, Y_PRTB, S_PRTB, CRV_PRTB, R_PRTB,
     >   CRV_TARGET, R_TARGET, S_TARGET, WT_TARGET,
     >   XOC_TARGET, YOC_TARGET

      END MODULE GEOM_MOD

************************************************************************
      MODULE OPT_MOD                          ! Optimization variables !
************************************************************************

      INTEGER, PARAMETER ::
     >   MXOBJS  = 3   ! Max. number of contributions to the objective

      INTEGER
     >   MODUS_OPERANDI,   NITMAX,  NDITER,  NALLOW,  NFTOTL,  NDV,
     >   NCLIN,   NCNLN,   NROWA,   NROWJ,   NROWR,   LENIW,   LENW,
     >   NBNDS,   INFORM,  MAJITS,  MAJORPL, MINORPL, MINORIL, NPRFREQ,
     >   LEVELDER, LEVELVER

      INTEGER, DIMENSION (MXOBJS) ::
     >   I1_OBJ,  I2_OBJ

      INTEGER, ALLOCATABLE, DIMENSION (:) ::
     >   ISTATE,  IWVEC,   INLCON,  JNLCON

      REAL
     >   ETA,     EPSOBJ,  OBJINI,  OBJMIN,  OBJ_CONFUN,
     >   STEPMX,  STEPLIM, TOLLIN,  TOLNLIN, TOLOPT,  PENPARAM,
     >   DIFFERENCE_INTERVAL

      REAL(KIND=4) :: CPU0

      REAL, DIMENSION (MXOBJS) ::
     >   RHO,     X1_OBJ,  X2_OBJ

      REAL, ALLOCATABLE, DIMENSION (:) ::
     >   FFWD,    GRAD,    CVEC,    WVEC,    AITCH,   ONEOVERH,
     >   BL,      BU,      CONLIN,  C,       CPRINT,
     >   BLLIN,   BULIN,   XNLCON,  SNLCON,  CLAMDA

      REAL, ALLOCATABLE, DIMENSION (:,:) ::
     >   AMAT,    RMAT,    CJAC

      LOGICAL
     >   SOLVE_CALLED

      CHARACTER * 10
     >   RHOTYPE(MXOBJS)

      CHARACTER, ALLOCATABLE, DIMENSION (:) ::
     >   LCTYPE*6, NLCTYPE*6

      END MODULE OPT_MOD

********************************************************************************
*
      PROGRAM NOSE_RADIUS
*
*        NOSE_RADIUS performs shape optimization for an airfoil using a target
*     radius of curvature distribution or possibly other geometric objective
*     functions and constraints.  It was prompted by a need to reduce excessive
*     heating near the wing leading edge of a small reusable launch vehicle
*     during reentry from orbit.
*
*        Control inputs are specified in 'nose_radius.inp' as follows:
*
*     PRIMARY CONTROL INPUTS:
*
*        MODE       Operating mode:  1 = target radius of curvature distribution
*                                    2 = target curvature distribution
*                                    3 = target airfoil
*        NDV        Number of design variables
*        NCLIN      Number of linear constraints
*        NCNLN      Number of nonlinear constraints
*        NITMAX     Number of optimization iterations >= 0
*        NALLOW     Upper limit on number of function calculations
*
*     OPTIMIZER INPUTS:
*
*        ETA        Controls line search's acceptance of a sufficiently
*                   lower objective. 0 < ETA < 1.0; try 0.1 or higher
*        STEPMX     Infinite step
*        EPSOBJ     Minimum absolute value of a significant difference
*                   in the objective function
*
*     DESIGN VARIABLE INPUTS:
*
*        The bulk of the design variables are Hicks-Henne-type "bumps."
*     Normally, they are multipliers of such perturbing shape functions,
*     although variable exponents may be an option some day.  Sample:
*
*        #  VTYPE   UP LW WIDTH XCNTR  XMIN  XMAX      V  VSCALE   H    BL  BU
*        1 'LEAD '  1. 1.  25.   999.  0.00  1.00  5.000   0.001 1.E-6 -10. 10.
*        2 'SIN2 '  0. 1.   2.  0.005  0.00  1.00  5.000   0.001 1.E-6 -10. 10.
*        3 'SIN2 '  0. 1.   2.  0.025  0.00  1.00  5.000   0.001 1.E-6 -10. 10.
*
*     #        Ordinal number of design variable - not actually used,
*              so it doesn't really matter if they get out of order
*     VTYPE    5-character design variable type as follows:
*
*        SCALE    Simple Y scaling
*        SIN      Standard (modified) sine perturbing function
*        SINC     Cyclic sine function suited to lofting heat shield
*                 perturbations in the azimuthal direction
*        SINF     Flipped sine function
*        SIN1     Symmetric sine function
*        SIN2     Symmetric flipped sine function
*        SIN3     Ensures airfoil LE/TE symmetry; use
*                 [XA, XB] = [0, 1] & XC in (0, .5)
*        SIN4     Ensures airfoil LE/TE symmetry; use
*                 [XA, XB] = [0, .5] & XC in (0, 1)
*        COSL     Half [co]sine, peak in [XA, XC], 0 at XB
*        COSR     Half [co]sine, peak in [XC, XB], 0 at XA
*        LCOS     Inverted form of COSL, peak in [XA, XC], 0 at XB
*        RCOS     Inverted form of COSR, peak in [XC, XB], 0 at XA
*        EXP      Exponential function (variable center)
*        DROOP    Simplified EXP with peak at X/C = 0;
*                 affects both surfaces the same way
*        LEAD     Leading  edge droop function (also LED; see also DROOP);
*                 affects both surfaces the same way
*        TRAIL    Trailing edge droop function (also TRL);
*                 affects both surfaces the same way
*        WAG      Wagner shape function; N is indicated via "POWER"
*
*     UP       Toggles for upper & lower surface application of shape function;
*     LO       DROOP, LEAD, TRAIL assume both surfaces are perturbed
*     WIDTH    Exponent or power:  higher means narrower bump;
*              use WIDTH = N for specifying Wagner function N
*     XCNTR    Location of peak
*     XMIN     Normally 0 and 1; may be less than full chord
*     XMAX
*     V        Design variable value as seen by the optimizer
*     VSCALE   Scale factor for the design variable to keep it ~O(1);
*              V * VSCALE is typically the shape function multiplier
*     H        Step size used for estimating the gradient of the objective
*              w.r.t. the variable by finite differencing;
*              the actual perturbation used is H * VSCALE;
*              H is also used for the nonlinear constraint derivatives;
*              H is NOT USED if Derivative Level = 0 because then NPOPT
*              determines an interval for each level UNLESS a positive
*              Difference Interval is entered via the NAMELIST (confusing!)
*     BL       Lower bound on the design variable as seen by the optimizer;
*              BL = -999. means the variable has no lower bound, else
*              BL scaling should match that of the variable - e.g., for a
*              variable with VSCALE = 0.1, BL = -10. means the effective
*              values are >= -10.*0.1 = -1., consistent with AITCH usage
*     BU       Upper bound on the design variable as seen by the optimizer;
*              BU = 999. means the variable has no upper bound, else BU
*              scaling should match that of the variable
*
*     OBJECTIVE FUNCTION INPUTS:
*
*        The "RHO"s refer to multipliers on terms added to the objective.
*     Some may have nonlinear constraint analogues.  Normally, picking one
*     form or the other is advisable (not both), but this is a gray area.
*
*        These contributions to the objective must all be present, but the
*     order is arbitrary.  Use a zero multiplier to suppress a contribution.
*     The accompanying pair of integers and pair of reals are available if
*     necessary to help define the objective contribution.
*
*     Format, where use of 999 is suggested to mean "not applicable":
*
*        TYPE       MULTIPLIER  I1_OBJ  I2_OBJ  X1_OBJ  X2_OBJ
*        RHO_***      0.001       999     999     999.    999.
*
*     Choice of contributions:
*
*        RHO_TARGET Target distribution (sum of squared differences)
*        RHO_VNORM  Squared 2-norm of the sine bump variables; small multipliers
*                   should tend to regularize the solution by preferring the
*                   shortest-length set of sine bump multipliers; try 1.E-4
*        RHO_R_TVD  Total squared variation of curvature radii for X/C in
*                   [X1_OBJ, X2_OBJ], to help smooth the result;
*                   presently applies to one surface at a time only: use
*                   I1_OBJ = 1 to mean upper surface; 2 to mean lower surface
*
*     NONLINEAR WEIGHTING OF TARGET DISTRIBUTIONS:
*
*        The idea here is to promote matching of the target near some X/C
*     while deemphasizing less critical targets.  Scaling by 1. leaves the
*     weighting of all target points at 1.  The "ADDTO" input is an offset
*     that counters the decay to zero of sinusoidal weighting at 0 and 1.
*
*              SRF  WTYPE  POWER  CENTER  XMIN  XMAX  PEAK  ADDTO
*              U   'SCALE'  999.    0.5     0.    1.    1.     0.
*              L   'LCOS '    3.    999.    0.   .15    10.    1.
*
*        The upper surface specs. should precede the lower surface specs.
*
*     LINEAR CONSTRAINTS:
*
*        If NCLIN = 0, just include two header lines.
*
*        #       Ordinal number of linear constraint - not used
*        LCTYPE  6-character linear constraint type:
*
*           XXXXXX  <None implemented yet.>
*
*        BL      Lower bound for linear constraint; -999. = -BIGBND
*        BU      Upper bound for linear constraint;  999. =  BIGBND
*
*     NONLINEAR CONSTRAINTS:
*
*        If NCNLN = 0, just include two header lines.
*
*        #          Ordinal number of nonlinear constraint - not used
*        NLCTYPE    6-character nonlinear constraint type:
*
*           RADIUS  Local radius of curvature at X/C = XNLCON in [0, 1];
*                   use INLCON = 1 for upper surface, 2 for lower surface
*           KAPPA   Local curvature alternative to RADIUS constraint, used
*                   analogously; note that kappa is mostly negative for a
*                   typical airfoil when evaluated vs. arc-length measured
*                   clockwise from lower to upper TE
*           MINRAD  Minimum local radius of curvature at any X/C in [0, X]
*                   where X <= 1.0 is entered via XNLCON;
*                   use INLCON = 1 for upper surface, 2 for lower surface
*           MINCRV  Minimum local curvature, intended to prevent undesired
*                   zero crossing of the signed 2-D curvature distribution;
*                   usage matches that of MINRAD, q.v.; see the note on the
*                   sign of 2-space curvature under KAPPA constraint above
*
*        BL         Lower bound on nonlinear constraint value;
*                   -999. = -BIGBND
*        BU         Upper bound on nonlinear constraint value;
*                   +999. =  BIGBND
*                   N.B.:  these bounds should be in the same units
*                   as the quantities being bounded;
*                   SNLCON will be applied by this program to the
*                   given BL and BU; this is in contrast to the
*                   bounds on the variables themselves, which
*                   (as for their finite differencing intervals) refer
*                   to the scaled variables
*        XNLCON     Real quantity needed (or not) by the nonlinear
*                   constraint
*        INLCON     First  index needed (or not) by the nonlinear
*                   constraint
*        JNLCON     Second index needed (or not) by the nonlinear
*                   constraint
*        SNLCON     Scale factor used to provide the optimizer with
*                   nonlinear constraint gradients of order 1
*
*     OPTIONAL INPUTS FOR NPOPT:
*
*        Certain of NPOPT's optional inputs may be entered via a namelist,
*     which should appear after the nonlinear constraint inputs.
*     Anything beyond that is ignored.
*
*        See the NPSOL or SNOPT User Guide for the control parameters in the
*     following namelist.  (Some of them apply to SNOPT only.)
*
*     NAMELIST /NPOPTIONS/
*    >   LEVELVER, MAJORPL, MINORPL, MINORIL, NPRFREQ,
*    >   LEVELDER, DIFFERENCE_INTERVAL, PENPARAM, STEPLIM,
*    >   TOLLIN, TOLNLIN, TOLOPT
*
*        STEPLIM    Limits stepsize of line search; default = 2.
*
*     INPUT/OUTPUT FILES:
*
*        See CONST_MOD module above and OPENs below.
*
*     GEOMETRY INPUT FORMAT ('nose_radius.dat'):
*
*        The airfoil should be normalized and wrap-around from lower TE
*     to upper TE.  Program PROFILE is handy for such reformatting.
*     Sample:
*
*        X-37 wing section at 64" station
*         225 ! Coordinates Clockwise
*           1.0000000      -0.1750951E-01
*           0.9772109      -0.1769160E-01
*           0.9490933      -0.1834753E-01
*           0.9141979      -0.1958385E-01
*           0.8983076      -0.2022068E-01
*           0.8474295      -0.2371158E-01
*            :               :
*
*     TARGET RADIUS DISTRIBUTION FORMAT (MODE = 1, 'nose_radius.target'):
*
*        Initially, treating radius of curvature as a function of X/C
*     appears most convenient.  The following format produced by program
*     CURVATURE was initially chosen for the target distribution:
*
*        !     S             K             R             X/C
*          0.000000E+00  6.013019E-01  1.663058E+00  9.995447E-01
*          2.233454E-02  6.013019E-01  1.663058E+00  9.772109E-01
*          5.045979E-02  3.835683E-01  2.607097E+00  9.490933E-01
*          8.537709E-02  1.826320E-01  5.475490E+00  9.141979E-01
*          1.012801E-01  8.576715E-01  1.165948E+00  8.983076E-01
*          1.513309E-01  3.508954E-01  2.849853E+00  8.483737E-01
*           :             :             :             :
*
*        However, the extra columns are inconvenient when R is being edited
*     as a target, so the format for both the input 'nose_radius.target' and
*     the output 'nose_radius.radius' files is simply this:
*
*        !     X/C           R
*          9.995447E-01  1.663058E+00
*          9.772109E-01  1.663058E+00
*          9.490933E-01  2.607097E+00
*          9.141979E-01  5.475490E+00
*          8.983076E-01  1.165948E+00
*          8.483737E-01  2.849853E+00
*           :             :
*
*        The radius of curvature values should apply to a normalized airfoil.
*
*     TARGET CURVATURE DISTRIBUTION FORMAT (MODE = 2, 'nose_radius.target'):
*
*        This option was added for completeness, and is analogous to the
*     target radius distribution case.
*
*     TARGET AIRFOIL FORMAT (MODE = 3, 'nose_radius.target'):
*
*        Target airfoil coordinates should be normalized and in the same format
*     as the input airfoil (see 'nose_radius.dat' above) except that one of the
*     surfaces may be omitted if it is not being modified.
*
*     SPONSOR:
*
*        James Reuther, Reacting Flow Environments Branch, NASA ARC.
*
*     HISTORY:
*
*        10/24/03  D.A.Saunders  Initial adaptatation of HEAT_SHIELD,
*                                which had the right high level design.
*        11/03/03    "      "    Added target airfoil option.
*        11/06/03    "      "    Save output curvatures as well as radii.
*        11/10/03    "      "    Added MINRAD constraint.
*        11/12/03    "      "    Switched to SNOPT (via its NPOPT interface),
*                                improving convergence behavior as hoped;
*                                added RHO_R_TVD contribution to the objective,
*                                requiring accompanying integer & real controls.
*        11/13/03    "      "    Added MINCRV constraint.
*        11/14/03    "      "    Added KAPPA constraint and target curvature
*                                distribution option for completeness.
*
*     AUTHOR:
*
*        David Saunders,  ELORET/NASA Ames Research Center, Moffett Field, CA
*
********************************************************************************

*     Global variables:

      USE CONST_MOD
      USE DESIGN_VARS_MOD
      USE OPT_MOD

      IMPLICIT NONE

*     Local constants:

      CHARACTER, PARAMETER ::
     >   PROGNAME*11 = 'NOSE_RADIUS'

*     Local variables:

      INTEGER ::
     >   I, NCALL
      REAL ::
     >   OBJPAR
      REAL(KIND=4) ::
     >   CPU2
      LOGICAL ::
     >   FAIL

*     Procedures:

      EXTERNAL
     >   NPOPT,  ! Constrained optimization, Stanford University
     >   SOLVE,  ! Function evaluation in QNMDIF form, used by ...
     >   CONFUN, ! ... and ...
     >   OBJFUN

*     **********
*     Execution:
*     **********

      CALL SECOND (CPU0)

      EPSMCH = EPSILON (EPSMCH)
      SMALL  = MAX (3.* EPSMCH, 1.E-13)
      BIGBND = 1.E+10

*     ***************
*     Open all files:
*     ***************

      OPEN (UNIT=LREAD,   FILE='nose_radius.inp',    STATUS='OLD')
      OPEN (UNIT=LGEOM,   FILE='nose_radius.dat',    STATUS='OLD')
      OPEN (UNIT=LTARGET, FILE='nose_radius.target', STATUS='OLD')
      OPEN (UNIT=LWRIT,   FILE='nose_radius.out',    STATUS='UNKNOWN')
      OPEN (UNIT=LGEOPT,  FILE='nose_radius.opt',    STATUS='UNKNOWN')
      OPEN (UNIT=LRADOPT, FILE='nose_radius.radius', STATUS='UNKNOWN')
      OPEN (UNIT=LCRVOPT, FILE='nose_radius.crv',    STATUS='UNKNOWN')

      WRITE (LWRIT, '(/, (A))')
     >   ' NOSE_RADIUS Version:  November 14, 2003',
     >   ' Sponsor:              Reacting Flow Environments Branch,',
     >   '                       NASA Ames Research Center'

*     ***************************
*     Read optimization controls:
*     ***************************

      CALL RD_OPT

*     **********************
*     Read initial geometry:
*     **********************

      CALL RD_GEOM

*     **********************
*     Read design variables:
*     **********************

      CALL RD_VARS

*     ******************************
*     Read objective function specs:
*     ******************************

      CALL RD_OBJECTIVE

*     **************************************
*     Read linear and nonlinear constraints:
*     **************************************

      CALL RD_CONS

*     *************************************************************
*     Set up the linear constraint matrix and corresponding bounds:
*     *************************************************************

      IF (NCLIN > 0) CALL SETLCON

*     *****************************
*     Read the target distribution:
*     *****************************

      CALL RD_TARGET (MODUS_OPERANDI)

*     *************************
*     Constrained optimization:
*     *************************

      NDITER = 0
      OBJMIN = BIGBND
      NCALL  = -3 ! Signals the first objective function call

      CALL SETUP_NPOPT

      IF (INFORM /= 0) GO TO 999

      CALL NPOPT (NDV, NCLIN, NCNLN, NROWA, NROWJ, NROWR,
     >            AMAT, BL, BU, CONFUN, OBJFUN,
     >            INFORM, MAJITS, ISTATE,
     >            CVEC, CJAC, CLAMDA, OBJPAR, GRAD, RMAT, V,
     >            IWVEC, LENIW, WVEC, LENW)

*     **************************************
*     Ensure we've got the optimized result?
*     **************************************

      IF (NITMAX > 0) THEN

         NCALL = -4 ! Force saving of optimized geometry
                    ! (not worth doing it every optimization iteration)

         CALL SOLVE (NDV, V, OBJPAR, NCALL, FAIL)

         WRITE (LWRIT, '(/, A, I7)')
     >      ' Total # objective calculations:', NFTOTL

      END IF

*     *************************************************
*     Save the optimized [radius of] curvature results:
*     *************************************************

      CALL SAVE_CRV

*     **********************************************
      WRITE (LWRIT, '(/, A)') ' Normal termination.'
*     **********************************************

  999 CALL SECOND (CPU2)

      CALL CPUTIME (CPU0, CPU2, PROGNAME, 'for this run', LWRIT)

      END PROGRAM NOSE_RADIUS

************************************************************************
*                                                                      *
      SUBROUTINE CHKLCON ()
*                                                                      *
*     Check the perturbed design against the linear constraints to     *
*     give a more readable table than is provided by the optimizer.    *
*                                                                      *
************************************************************************

*     Global variables:

      USE CONST_MOD
      USE DESIGN_VARS_MOD
      USE OPT_MOD

      IMPLICIT NONE

*     Local variables:

      INTEGER   INTERP, M
      REAL      BOUNDL, BOUNDU, T, VALUE
      CHARACTER FLAG*3

*     Execution:

      WRITE (LWRIT, 1010)

*     Evaluate the linear constraints in the form seen by NPOPT:

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

         WRITE (LWRIT, 1030) M, LCTYPE(M), ILCON(M), VALUE,
     >                       BOUNDL, BOUNDU, CONLIN(M), FLAG
      END DO

      RETURN

 1010 FORMAT (/, ' CHKLCON: Linear constraint check:', //,
     >   '    #  LCTYPE   ILCON  CURRENT VALUE',
     >   '    LOWER BOUND    UPPER BOUND           AMAT * V')
 1030 FORMAT (1X, I4, 2X, A6, I6, 2X, 3F15.6, E19.10, 2X, A3)

      END SUBROUTINE CHKLCON

************************************************************************
*                                                                      *
      SUBROUTINE CONFUN (MODE, NCPAR, NDVPAR, NROWJPAR, NEEDC, VPAR,
     >                   CPAR, CJACPAR, NSTATE)
*                                                                      *
*     CONFUN computes the nonlinear constraint functions C(*).         *
*     This version assumes the optimizer determines its own first      *
*     derivatives via finite differencing.                             *
*                                                                      *
*     This version ignores NEEDC(*) -  it is more efficient to         *
*     evaluate all constraints in groups than to evaluate each one     *
*     independently.                                                   *
*                                                                      *
************************************************************************

*     Global quantities:

      USE CONST_MOD
      USE OPT_MOD

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

      REAL,      PARAMETER :: ZERO = 0.
      LOGICAL,   PARAMETER :: NOPRINT = .FALSE.
      CHARACTER, PARAMETER :: SUBNAME * 6 = 'CONFUN'

*     Local variables:

      INTEGER    I, J, NCALL
      REAL       OBJPAR, VTEMP
      REAL(KIND=4) :: TIME1, TIME2
      LOGICAL    FAIL

*     Execution:

****  WRITE (LWRIT, '(/, A, 2I4)')
**** >   ' CONFUN: MODE and NSTATE: ', MODE, NSTATE

      IF (NSTATE == 1) THEN ! First call to CONFUN (precedes OBJFUN)

         NCALL = -13

      ELSE ! Not the first call to CONFUN

         IF (MODE == 0) THEN ! Constraints only
            NCALL = -10
         ELSE                ! Constraints + derivatives
            NCALL = -12
         END IF

      END IF

*     Analyze the current geometry:

      CALL SOLVE (NDVPAR, VPAR, OBJPAR, NCALL, FAIL)

      OBJ_CONFUN = OBJPAR

      NFTOTL = NFTOTL + 1

      IF (FAIL) THEN
         WRITE (LWRIT, '(/, 1X, 2A)')
     >      SUBNAME, ':  Bad return from SOLVE.'
         GO TO 900
      END IF

*     Evaluate all nonlinear constraints, scaled:

      CALL NLCON (NOPRINT, CPAR)

      IF (NITMAX == 0) THEN
         IF (SOLVE_CALLED) THEN
            NCALL = -4 ! Force saving of perturbed geometry

            CALL SUMMARY (OBJPAR, NCALL)

            WRITE (LWRIT, '(/, 1X, 2A)')
     >         SUBNAME, ':  Forcing exit from NPOPT (NITMAX = 0).'
            GO TO 900
         END IF
      END IF

*     Constraint derivatives:

      IF (MODE > 0 .AND. LEVELDER /= 0) THEN

         CALL SECOND (TIME1)

*        Perturb a given variable only once for all constraints:

         DO J = 1, NDV

            VTEMP   = VPAR(J)
            VPAR(J) = VTEMP + AITCH(J)

            CALL SOLVE (NDVPAR, VPAR, OBJPAR, J, FAIL) ! NCALL = J here

            VPAR(J) = VTEMP

            CALL NLCON (NOPRINT, CJACPAR(1,J))

            DO I = 1, NCNLN
               CJACPAR(I,J) = (CJACPAR(I,J) - CPAR(I)) * ONEOVERH(J)
            END DO

            FFWD(J) = (OBJPAR - OBJ_CONFUN) * ONEOVERH(J)

         END DO

         NFTOTL = NFTOTL + NDV

         CALL SECOND (TIME2)

***      IF (SOLVE_CALLED) THEN
***         CALL CPUTIME (TIME1, TIME2, SUBNAME,
***  >   'to evaluate derivatives of nonlinear constraints & objective',
***  >         LWRIT)
***      ELSE
***         CALL CPUTIME (TIME1, TIME2, SUBNAME,
***  >         'to evaluate nonlinear constraint derivatives', LWRIT)
***      END IF

         WRITE (LWRIT, '(/, A)' ) ' CONFUN: Jacobian matrix transpose:',
     >      ' IVAR       dC1/dVAR       dC2/dVAR ...'

         I = MIN (NCNLN, 8)

         DO J = 1, NDV
            WRITE (LWRIT, '(I5, 1P, 8E15.6)') J, CJACPAR(1:I,J)
         END DO

      END IF

      GO TO 999

  900 MODE = -1 ! NPOPT should quit

  999 RETURN

      END SUBROUTINE CONFUN

************************************************************************
*                                                                      *
      SUBROUTINE CPUTIME (CPU1, CPU2, CALLER, DESCRIPTION, LUN)
*                                                                      *
*     CPUTIME modularizes printing of CPU time between two measures of *
*     total CPU time used so far as given by the SECOND utility.       *
*                                                                      *
************************************************************************

      IMPLICIT NONE

*     Arguments:

      REAL(KIND=4), INTENT (IN) :: CPU1, CPU2
      CHARACTER,    INTENT (IN) :: CALLER * (*), DESCRIPTION * (*)
      INTEGER,      INTENT (IN) :: LUN

*     Local variables:

      REAL DIFF

*     Execution:

      DIFF = CPU2 - CPU1
      IF (DIFF >= 0.01) THEN
         WRITE (LUN, '(/, 1X, 4A, F10.2)')
     >      CALLER, ': CPU secs. ', DESCRIPTION, ':   ', DIFF
      END IF

      END SUBROUTINE CPUTIME

************************************************************************
*                                                                      *
      SUBROUTINE KAPPA (N, X, Y, S, CRV, RADIUS)
*                                                                      *
*     Calculate an (x,y) curve curvature distribution via arc-length-  *
*     based 3-point finite difference derivatives.  At an end point,   *
*     these formulas give the same value as at the neighboring point.  *
*     The radius of curvature is also returned - safeguarded against   *
*     zero curvature.                                                  *
*                                                                      *
************************************************************************

      IMPLICIT NONE

*     Arguments:

      INTEGER, INTENT (IN)  :: N         ! # pts. in the curve
      REAL,    INTENT (IN)  :: X(N),     ! (x,y) coordinates
     >                         Y(N)
      REAL,    INTENT (OUT) :: S(N),     ! Arc lengths and
     >                         CRV(N),   ! curvatures and
     >                         RADIUS(N) ! radii
*     Local constants:

      REAL,    PARAMETER :: CURV_MIN  = 1.E-10, ONE = 1.
      LOGICAL, PARAMETER :: NORMALIZE = .FALSE.

*     Local variables:

      INTEGER  I
      REAL     STOTAL
      REAL, DIMENSION (N) :: XS, XSS, YS, YSS

*     Execution:

*     Unnormalized arc lengths:

      CALL CHORDS2D (N, X, Y, NORMALIZE, STOTAL, S)

*     Central differences at interior points:

      DO I = 2, N - 1
         CALL FDCNTR (I, S, X, XS(I), XSS(I))
         CALL FDCNTR (I, S, Y, YS(I), YSS(I))
      END DO

*     Signed local curvatures:

      CALL CURV2D (N - 2, XS(2), XSS(2), YS(2), YSS(2), CRV(2))

*     Unsigned radii:

      DO I = 2, N - 1
         RADIUS(I) = ONE / MAX (ABS (CRV(I)), CURV_MIN)
      END DO

*     End points:

      CRV(1)    = CRV(2);     CRV(N)    = CRV(N-1)
      RADIUS(1) = RADIUS(2);  RADIUS(N) = RADIUS(N-1)

      END SUBROUTINE KAPPA

************************************************************************
*                                                                      *
      SUBROUTINE NLCON (PRINT, CPAR)
*                                                                      *
*     NLCON evaluates all nonlinear constraints, CPAR(*), with the     *
*     option to tabulate UNscaled results.                             *
*                                                                      *
************************************************************************

*     Global variables:

      USE CONST_MOD
      USE GEOM_MOD
      USE OPT_MOD

      IMPLICIT NONE

*     Arguments:

      LOGICAL, INTENT (IN)  :: PRINT   ! T means tabulate UNscaled constraints
      REAL,    INTENT (OUT) :: CPAR(*) ! Nonlinear constraints (= C(*) in OPT)

*     Local constants:

      REAL, PARAMETER ::
     >   ZERO = 0.

      LOGICAL, PARAMETER ::
     >   NEW = .TRUE.

      CHARACTER, PARAMETER ::
     >   TIGHT * 1 = 'M'

*     Local variables:

      INTEGER   I, I1, I2, K, M, N
      REAL      BOUNDL, BOUNDU, DERIV, R_MIN, SCALE, TIME, VALUE, XOC
      CHARACTER FLAG*3

*     Execution:

      K = NDV + NCLIN ! Offset for the bounds

      IF (PRINT) WRITE (LWRIT, 1010)

      DO M = 1, NCNLN

         SELECT CASE (NLCTYPE(M))

            CASE ('RADIUS') ! Radius of curvature at specified X

*              Interpolate the current radius distribution at the indicated X/C:

               IF (INLCON(M) == 1) THEN ! Upper surface
                  I = NL; N = NU
               ELSE                     ! Lower surface
                  I = 1;  N = NL
               END IF

               CALL LCSFIT (N, X_PRTB(I), R_PRTB(I), NEW, TIGHT,
     >                      1, XNLCON(M), VALUE, DERIV)

***            write (lwrit, '(a, i1, 2f12.6)')
***  >            ' Surface, X/C, radius: ', INLCON(M), XNLCON(M), VALUE

            CASE ('KAPPA ') ! Curvature at specified X

*              Interpolate the current kappa distribution at the indicated X/C:

               IF (INLCON(M) == 1) THEN ! Upper surface
                  I = NL; N = NU
               ELSE                     ! Lower surface
                  I = 1;  N = NL
               END IF

               CALL LCSFIT (N, X_PRTB(I), CRV_PRTB(I), NEW, TIGHT,
     >                      1, XNLCON(M), VALUE, DERIV)

***            write (lwrit, '(a, i1, 2f12.6)')
***  >            ' Surface, X/C, kappa: ', INLCON(M), XNLCON(M), VALUE

            CASE ('MINRAD') ! Minimum radius of curvature in X/C = [0, XNLCON]

*              Interpolate the current radius distribution at the indicated X/C:

               IF (INLCON(M) == 1) THEN ! Upper surface
                  I1 = NL; I2 = NU
               ELSE                     ! Lower surface
                  I1 = 1;  I2 = NL
               END IF

               R_MIN = 1.E+30
               DO I = I1, I2
                  IF (X_PRTB(I) <= XNLCON(M)) THEN
                     IF (R_PRTB(I) < R_MIN) THEN
                        R_MIN = R_PRTB(I)
                        XOC   = X_PRTB(I)
                     END IF
                  END IF
               END DO
               VALUE = R_MIN

***            write (lwrit, '(a, i1, 2f12.6)')
***  >            ' Surface, X/C, min. R: ', INLCON(M), XOC, VALUE

            CASE ('MINCRV') ! Minimum |curvature| in X/C = [0, XNLCON]

*              Interpolate the current kappa distribution at the indicated X/C:

               IF (INLCON(M) == 1) THEN ! Upper surface
                  I1 = NL; I2 = NU
               ELSE                     ! Lower surface
                  I1 = 1;  I2 = NL
               END IF

               R_MIN = -1.E+30 ! We want kappa to stay below zero
               DO I = I1, I2
                  IF (X_PRTB(I) <= XNLCON(M)) THEN
                     IF (CRV_PRTB(I) > R_MIN) THEN
                        R_MIN = CRV_PRTB(I)
                        XOC   = X_PRTB(I)
                     END IF
                  END IF
               END DO
               VALUE = R_MIN

***            write (lwrit, '(a, i1, 2f12.6)')
***  >            ' Surface, X/C, min. kappa: ', INLCON(M), XOC, VALUE

            CASE DEFAULT

               WRITE (LWRIT, '(/, 2A)')
     >            ' NLCON: Bad constraint type = ', NLCTYPE(M)
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

            WRITE (LWRIT, 1030) M, NLCTYPE(M), INLCON(M), XNLCON(M),
     >                          VALUE, BOUNDL, BOUNDU, FLAG
          END IF

      END DO

      RETURN

 1010 FORMAT (/, ' NLCON: Nonlinear constraint check:', //,
     >   '    # NLCTYPE  INLCON  XNLCON  CURRENT VALUE',
     >   '    LOWER BOUND    UPPER BOUND   VIOLATION?')
 1030 FORMAT (1X, I4, 2X, A6, I6, F10.4, 3F15.6, 10X, A3)

      END SUBROUTINE NLCON

************************************************************************
*                                                                      *
      SUBROUTINE OBJECTIVE (NDVPAR, VPAR, OBJ, NCALL)
*                                                                      *
*     OBJECTIVE is called by SOLVE.  It keeps the objective function   *
*     details in one place.                                            *
*                                                                      *
************************************************************************

*     Global variables:

      USE CONST_MOD
      USE DESIGN_VARS_MOD
      USE GEOM_MOD
      USE OPT_MOD

      IMPLICIT NONE

*     Arguments:

      INTEGER, INTENT (IN)  :: NDVPAR         ! Dimension of VPAR(*)
      REAL,    INTENT (IN)  :: VPAR(NDVPAR)   ! Current variables,
                                              ! not necessarily V(*)
      REAL,    INTENT (OUT) :: OBJ            ! Current objective
      INTEGER, INTENT (IN)  :: NCALL          ! Controls printing here

*     Local constants:

      REAL, PARAMETER :: ZERO = 0.
      CHARACTER, PARAMETER :: NAME * 12 = ' OBJECTIVE: '

*     Local variables:

      INTEGER I, I1, I2, N
      REAL    CONTRIBUTION(MXOBJS), SUMSQ

*     Execution:

*     Evaluate the current curvature distribution.
*     Do it for all modes in case some local RADIUS constraints are present.

      CALL KAPPA (NW, X_PRTB, Y_PRTB, S_PRTB, CRV_PRTB, R_PRTB)

      OBJ = ZERO

      DO N = 1, MXOBJS ! For each possible contribution to the objective ...

         IF (RHO(N) == ZERO) THEN
            CONTRIBUTION(N) = ZERO
            CYCLE
         END IF

         SELECT CASE (RHOTYPE(N))

         CASE ('RHO_TARGET') ! Sum of squares difference from target distribn.

*           Unwrap the lower surface to simplify interpolation at target X/Cs:

            IF (NLTARG /= 1) X_PRTB(1:NL) = -X_PRTB(1:NL)

            SELECT CASE (MODUS_OPERANDI)

            CASE (1) ! Target radius distribution

               CALL COMPARE_DISTRIBUTIONS
     >            (NW, X_PRTB, R_PRTB, NTARGET, XOC_TARGET, R_TARGET)

            CASE (2) ! Target curvature distribution

               CALL COMPARE_DISTRIBUTIONS
     >            (NW, X_PRTB, CRV_PRTB, NTARGET, XOC_TARGET,CRV_TARGET)

            CASE (3) ! Target (x,y)s

               CALL COMPARE_DISTRIBUTIONS
     >            (NW, X_PRTB, Y_PRTB, NTARGET, XOC_TARGET, YOC_TARGET)

            CASE DEFAULT

               WRITE (LWRIT, '(/, 2A, I10)')
     >            NAME, 'Bad operating mode: ', MODUS_OPERANDI
            STOP

            END SELECT

            IF (NLTARG /= 1) X_PRTB(1:NL) = -X_PRTB(1:NL) ! Restore lower X/Cs

            CONTRIBUTION(N) = RHO(N) * SUMSQ

         CASE ('RHO_VNORM ') ! 2-norm-squared of sine bump variables

            SUMSQ = ZERO
            DO I = 1, NBUMPS
               SUMSQ = SUMSQ + VPAR(I) ** 2
            END DO
            CONTRIBUTION(N) = RHO(N) * SUMSQ

         CASE ('RHO_R_TVD ') ! Total squared variation of radii in [X1, X2]

            SELECT CASE (I1_OBJ(N))

            CASE (1) ! Upper surface

               I1 = NL + 1; I2 = NW

            CASE (2) ! Lower surface

               I1 = 2;  I2 = NL

            CASE DEFAULT

               WRITE (LWRIT, '(/, 2A, 2I10)')
     >            NAME, 'Bad I1_OBJ: ', I1_OBJ(N), N
               STOP

            END SELECT

            SUMSQ = ZERO
            DO I = I1, I2
               IF (X_PRTB(I) >= X1_OBJ(N) .AND.
     >             X_PRTB(I) <= X2_OBJ(N)) THEN
                  SUMSQ = SUMSQ + (R_PRTB(I) - R_PRTB(I-1)) ** 2
               END IF
            END DO
            CONTRIBUTION(N) = RHO(N) * SUMSQ

         CASE DEFAULT

            WRITE (LWRIT, '(/, 3A, I10)')
     >         NAME, 'Bad RHO type: ', RHOTYPE(N), N
            STOP

         END SELECT

         OBJ = OBJ + CONTRIBUTION(N)

      END DO ! Next contribution to the objective

      IF (NCALL < 0) THEN
         WRITE (LWRIT, '(/, A, //, (1X, A, 1P, E25.15))')
     >      ' Contributions to the objective:',
     >      (RHOTYPE(N), CONTRIBUTION(N), N = 1, MXOBJS)
      END IF

      RETURN

*     Local procedure for subroutine OBJECTIVE:

      CONTAINS

         SUBROUTINE COMPARE_DISTRIBUTIONS (N, X, Y,
     >                                     NTARGET, XTARGET, YTARGET)
!        Arguments:

         INTEGER, INTENT (IN) ::
     >      N, NTARGET        ! # points in current and target distributions
         REAL, INTENT (IN) ::
     >      X(N), Y(N),       ! Current distribution
     >      XTARGET(NTARGET), ! Target distribution
     >      YTARGET(NTARGET)

!        Automatic arrays:

         REAL, DIMENSION (NTARGET) ::
     >      YINTERP, DERIVS   ! Current Y interpolated at XTARGET(*),
                              ! and unused derivatives
!        Local constants:

         LOGICAL,   PARAMETER :: NEW = .TRUE.
         CHARACTER, PARAMETER :: TIGHT * 1 = 'M'

!        Execution:

!        Interpolate the current Y distribution at the target Xs:

         CALL LCSFIT (N, X, Y, NEW, TIGHT,
     >                NTARGET, XTARGET, YINTERP, DERIVS)

!        Weighted sum of squared differences:

         SUMSQ = ZERO
         DO I = 1, NTARGET
            SUMSQ = SUMSQ +  WT_TARGET(I) *
     >              (YINTERP(I) - YTARGET(I)) ** 2
         END DO

         END SUBROUTINE COMPARE_DISTRIBUTIONS

      END SUBROUTINE OBJECTIVE

************************************************************************
*                                                                      *
      SUBROUTINE OBJFUN (MODE, NDVPAR, VPAR, OBJPAR, GRDPAR, NSTATE)
*                                                                      *
*     OBJFUN acts as the interface between NPOPT and SOLVE.            *
*     See also CONFUN if nonlinear constraints are present.            *
*     Its main task is to determine whether NPOPT wants just           *
*     the objective function or the objective function plus            *
*     gradient information, and then call SOLVE appropriately.         *
*                                                                      *
*     3/96      J. Reuther  Initial implementation (SYN87).            *
*     2/00      D. Saunders Initial adaptations (Traj_opt).            *
*                                                                      *
************************************************************************

*     Global variables:

      USE CONST_MOD
      USE OPT_MOD

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

      INTEGER :: I, IENTRY, NCALL, NHALF
      REAL    :: OBJNEW, TEMPH, TEMPV
      LOGICAL :: FAIL

*     Execution:

*     ---------------------------------------------------------------
*     Treat the case in which nonlinear constraints have forced calls
*     to SOLVE, thus providing OBJPAR and GRDPAR.  See CONFUN.
*     ---------------------------------------------------------------

      IF (SOLVE_CALLED) THEN

         OBJPAR = OBJ_CONFUN

         IF (MODE > 0 .AND. LEVELDER /= 0) THEN

            GRDPAR = FFWD ! Elements 1 : NDV

            WRITE (LWRIT, 1030)
     >         'OBJFUN: Objective gradient via CONFUN:', ' ',
     >         '   I            V (I)            G (I)            H (I)'
            WRITE (LWRIT, '(I5, 3E17.8)')
     >         (I, VPAR(I), GRDPAR(I), AITCH(I), I = 1, NDVPAR)
         END IF

      ELSE

*        -------------------------------------------------------------
*        Case of no nonlinear constraints, or above option suppressed:
*        -------------------------------------------------------------

***      WRITE (LWRIT, 1040) 'MODE and NSTATE: ', MODE, NSTATE

*        Calculate the objective function.

         IF (NSTATE == 1) THEN

            NCALL = -3

            CALL SOLVE (NDVPAR, VPAR, OBJPAR, NCALL, FAIL)

            NFTOTL = NFTOTL + 1

            IF (FAIL) THEN
               WRITE (LWRIT, 1040) 'Bad return from call to SOLVE.'
               GO TO 900
            END IF

            IF (NITMAX == 0) THEN
               NCALL = -4 ! Force saving of perturbed geometry

               CALL SUMMARY (OBJPAR, NCALL)

               WRITE (LWRIT, 1040) 'Exiting from NPOPT (NITMAX = 0).'
               GO TO 900
            END IF

         ELSE

            IF (MODE == 0) THEN ! Objective function only
               NCALL = 0        ! Line search function evaluation
            ELSE                ! Objective and gradient
               NCALL = -2       ! "Central" solution after a line search
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

*        Calculate the gradient if it is requested.

         IF (MODE > 0 .AND. LEVELDER /= 0) THEN

            WRITE (LWRIT, 1040) 'Calculate the objective gradient.'

            DO NCALL = 1, NDVPAR

               TEMPV = VPAR(NCALL)
               TEMPH = AITCH(NCALL)
               NHALF = 0

   10          CONTINUE

               VPAR(NCALL) = TEMPV + TEMPH

               CALL SOLVE (NDVPAR, VPAR, OBJNEW, NCALL, FAIL)

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

                  OBJNEW = OBJPAR
                  WRITE (LWRIT, 1030)
     >               'Three halvings of gradient step size failed.',
     >               'Proceeding with zero for this component.'
               END IF

               GRDPAR(NCALL) = (OBJNEW - OBJPAR) / TEMPH
               VPAR(NCALL)   = TEMPV

            END DO

            WRITE (LWRIT, 1030) 'OBJFUN: Gradient vector:', ' ',
     >      '   I            V (I)            G (I)            H (I)'
            WRITE (LWRIT, '(I5, 3E17.8)')
     >         (I, VPAR(I), GRDPAR(I), AITCH(I), I = 1, NDVPAR)

         END IF

      END IF

      GO TO 999

  900 MODE = -1 ! NPOPT should quit

  999 RETURN

 1030 FORMAT (/, (1X, A))
 1040 FORMAT (/, ' OBJFUN: ', A, 2I4)

      END SUBROUTINE OBJFUN

************************************************************************
*                                                                      *
      SUBROUTINE PERTURB (VPAR, NCALL)
*                                                                      *
*     Variant for a normalized airfoil.                                *
*                                                                      *
*     Note:  SOLVE passes MCALL >= -3, not NCALL >= -13.               *
*                                                                      *
************************************************************************

      USE CONST_MOD
      USE DESIGN_VARS_MOD
      USE GEOM_MOD

      IMPLICIT NONE

*     Arguments:

      REAL,    INTENT (IN) ::
     >   VPAR(*) ! Optimization variables

      INTEGER, INTENT (IN) ::
     >   NCALL   ! -3 means first call

*     Local constants:

      REAL, PARAMETER ::
     >   ZERO = 0., ONE = 1.

*     Local variables:

      INTEGER ::
     >   I, N

      REAL ::
     >   DL, DU, VS

      CHARACTER ::
     >   TYP4 * 4

*     Execution:

*     --------------------------------------
*     Initialize as the unperturbed airfoil:
*     --------------------------------------

      DO I = 1, NW
         X_PRTB(I) = X_INIT(I)
         Y_PRTB(I) = Y_INIT(I)
      END DO

*     -----------------------------
*     Apply perturbations in-place:
*     -----------------------------

      DO N = 1, NBUMPS

         VS = VPAR(N) * VSCALE(N)
         IF (ABS (VS) < SMALL) CYCLE ! Next N

         TYP4 = XBTYP(N)(1:4)

C        Some variables will perturb the leading edge twice if
C        we don't perturb both surfaces in one call, so:

         IF (TYP4(1:2) == 'LE' .OR. TYP4(1:2) == 'DR') THEN

            CALL SFEVAL (TYP4, POWER(N), XBUMP(N), XBMIN(N), XBMAX(N),
     >                   VS, 1, NW, X_PRTB, Y_PRTB, LWRIT)
         ELSE

            IF (DOLO(N) /= ZERO) THEN

               DL = VS * SIGN (ONE, DOLO(N)) ! Not clear why

               CALL SFEVAL (TYP4, POWER(N), XBUMP(N), XBMIN(N),XBMAX(N),
     >                      DL, 1, NL, X_PRTB, Y_PRTB, LWRIT)
            END IF

            IF (DOUP(N) /= ZERO) THEN

               DU = VS * SIGN (ONE, DOUP(N))

               CALL SFEVAL (TYP4, POWER(N), XBUMP(N), XBMIN(N),XBMAX(N),
     >                      DU, NL, NW, X_PRTB, Y_PRTB, LWRIT)
            END IF

         END IF

      END DO ! Next design variable

      END SUBROUTINE PERTURB

************************************************************************
*                                                                      *
      SUBROUTINE RD_CONS
*                                                                      *
*     Read the constraints.                                            *
*                                                                      *
************************************************************************

      USE CONST_MOD
      USE DESIGN_VARS_MOD
      USE OPT_MOD

      IMPLICIT NONE

*     Local variables:

      INTEGER :: J, K, M, IOS
      REAL    :: BOUNDLO, BOUNDUP

*     Execution:
*     ----------

*     Linear constraints:

      READ (LREAD, *)
      READ (LREAD, *)
      READ (LREAD, *)

      DO M = 1, NCLIN

         READ (LREAD, *, IOSTAT=IOS) J, LCTYPE(M), BOUNDLO, BOUNDUP,
     >                               TLCON(M), ILCON(M)
         IF (IOS /= 0) THEN
            WRITE (LWRIT, *) 'Error reading linear constraints.'
            GO TO 999
         END IF

         IF (INT (BOUNDLO) == -999) BOUNDLO = -BIGBND
         IF (INT (BOUNDUP) ==  999) BOUNDUP =  BIGBND
         BLLIN(M) = BOUNDLO
         BULIN(M) = BOUNDUP

      END DO

*     Nonlinear constraints:

      READ (LREAD, *)
      READ (LREAD, *)
      READ (LREAD, *)

      K = NDV + NCLIN ! Offset for the bounds

      DO M = 1, NCNLN

         READ (LREAD, *, IOSTAT=IOS) J, NLCTYPE(M), BOUNDLO, BOUNDUP,
     >      XNLCON(M), INLCON(M), JNLCON(M), SNLCON(M)

         IF (IOS /= 0) THEN
            WRITE (LWRIT, *) J, NLCTYPE(M), BOUNDLO, BOUNDUP,
     >         XNLCON(M), INLCON(M), JNLCON(M), SNLCON(M)
            WRITE (LWRIT, *) 'IOS: ', IOS
            WRITE (LWRIT, *) 'Error reading nonlinear constraints.'
            GO TO 999
         END IF

*        Scaling of the nonlinear constraints differs from scaling of
*        the variables and linear constraints:

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

************************************************************************
*                                                                      *
      SUBROUTINE RD_GEOM
*                                                                      *
*     Read an airfoil in wrap-around form as output by PROFILE.        *
*     NOSE_RADIUS assumes it is normalized, as required for applying   *
*     shape functions.  Sample:                                        *
*                                                                      *
*        X-37 wing section at 64" station                              *
*         225 ! Coordinates Clockwise                                  *
*           1.0000000      -0.1750951E-01                              *
*           0.9772109      -0.1769160E-01                              *
*           0.9490933      -0.1834753E-01                              *
*            :               :                                         *
*                                                                      *
************************************************************************

      USE CONST_MOD
      USE GEOM_MOD

      IMPLICIT NONE

*     Local constants:

      REAL,      PARAMETER :: ONE = 1., ZERO = 0.
      CHARACTER, PARAMETER :: FILENAME * 15 = 'nose_radius.dat'

*     Local variables:

      INTEGER I, IOS
      REAL    XLE

*     Execution:

      OPEN (LGEOM, FILE=FILENAME, STATUS='OLD', IOSTAT=IOS)

      IF (IOS /= 0) THEN
         WRITE (LWRIT, '(/, 2A)') ' Trouble opening ', filename
         GO TO 999
      END IF

      READ (LGEOM, '(A)') ! Title
      READ (LGEOM, *, IOSTAT=IOS) NW

      IF (IOS /= 0) THEN
         WRITE (LWRIT, '(/, 2A)') ' Trouble reading NW from ', filename
         GO TO 999
      END IF

      ALLOCATE
     >  (X_INIT(NW), Y_INIT(NW), S_INIT(NW), CRV_INIT(NW), R_INIT(NW),
     >   X_PRTB(NW), Y_PRTB(NW), S_PRTB(NW), CRV_PRTB(NW), R_PRTB(NW))

      XLE = 1.E+10

      DO I = 1, NW
         READ (LGEOM, *, IOSTAT=IOS) X_INIT(I), Y_INIT(I)
         IF (IOS /= 0) THEN
            WRITE (LWRIT, '(/, A, I4)')
     >         ' Trouble reading airfoil point ', I
            GO TO 999
         END IF
         IF (X_INIT(I) < XLE) THEN
            XLE = X_INIT(I)
            NL = I
         END IF
      END DO

      NU = NW - NL + 1

      CLOSE (LGEOM)

      IF (XLE /= ZERO .OR. MAX (X_INIT(1), X_INIT(NW)) /= ONE) THEN
         WRITE (LWRIT, '(/, A)')
     >      ' Airfoil is not normalized.  Aborting.'
         GO TO 999
      END IF

      IF (X_INIT(NL-1) == XLE .OR. X_INIT(NL+1) == XLE) THEN
         WRITE (LWRIT, '(/, A)')  ! Breaks splining, etc.
     >      ' Nonunique leading edge point?  Aborting.'
         GO TO 999
      END IF

      IF (Y_INIT(NW) < Y_INIT(1)) THEN
         WRITE (LWRIT, '(/, 2A)') ' Airfoil does NOT appear wrapped',
     >      ' from lower TE to upper TE.  Aborting.'
         GO TO 999
      END IF

*     Evaluate the initial airfoil's curvature distribution:

      CALL KAPPA (NW, X_INIT, Y_INIT, S_INIT, CRV_INIT, R_INIT)

      RETURN

  999 STOP

      END SUBROUTINE RD_GEOM

************************************************************************
*                                                                      *
      SUBROUTINE RD_OBJECTIVE
*                                                                      *
*     Read objective function specifications.                          *
*                                                                      *
************************************************************************

      USE CONST_MOD
      USE GEOM_MOD
      USE OPT_MOD

      IMPLICIT NONE

*     Local variables:

      INTEGER
     >   I
      CHARACTER * 1
     >   SRF(2)

*     Execution:

      READ (LREAD, *) ! Contributions to the objective
      READ (LREAD, *)
      READ (LREAD, *)

      DO I = 1, MXOBJS
         READ (LREAD, *) RHOTYPE(I), RHO(I),
     >      I1_OBJ(I), I2_OBJ(I), X1_OBJ(I), X2_OBJ(I)
      END DO

      READ (LREAD, *) ! Nonlinear weighting of target distributions
      READ (LREAD, *)
      READ (LREAD, *)

      DO I = 1, 2
         READ (LREAD, *) SRF(I), WBTYP(I),
     >      WPOWER(I), WXCNTR(I), WXMIN(I), WXMAX(I), WPEAK(I), ADDTO(I)
      END DO

      IF (SRF(1) /= 'U' .AND. SRF(1) /= 'u') THEN
         WRITE (LWRIT, '(/, A)')
     >      ' Enter upper surface target weighting spec. first.'
         STOP
      END IF

      END SUBROUTINE RD_OBJECTIVE

************************************************************************
*                                                                      *
      SUBROUTINE RD_OPT
*                                                                      *
*     Read the primary optimization controls.                          *
*                                                                      *
************************************************************************

*     Global variables:

      USE CONST_MOD
      USE DESIGN_VARS_MOD
      USE GEOM_MOD
      USE OPT_MOD

      IMPLICIT NONE

*     Local variables:

      INTEGER   :: I, IOS
      CHARACTER :: BUFFER*132

*     Execution:

*     Echo the control inputs to the output file:

      WRITE (LWRIT, '(/, A, /)') ' Optimization control inputs:'

      DO
         READ (LREAD, '(A)', IOSTAT=IOS) BUFFER
         IF (IOS < 0) EXIT

         I = LEN_TRIM (BUFFER)
         WRITE (LWRIT, '(1X, A)') BUFFER(1:I)
      END DO

      REWIND (LREAD)

*     Read the control inputs:

      READ (LREAD, *)
      READ (LREAD, *)
      READ (LREAD, *)
      READ (LREAD, *) MODUS_OPERANDI, NDV, NCLIN, NCNLN, NITMAX, NALLOW
      READ (LREAD, *)
      READ (LREAD, *) ETA, STEPMX, EPSOBJ

*     ------------------------------------
*     Set up various sizes and dimensions:
*     ------------------------------------

      NBNDS = NDV + NCLIN + NCNLN      ! For variable + constraint bounds
      NROWA = MAX (NCLIN, 1)           ! Max. # linear constraints ...
      NROWJ = MAX (NCNLN, 1)           ! ... & nonlinear constraints; dims. > 0
      NROWR = MAX (NDV, 1)             ! Row dim. of Cholesky factor R
                                       ! of the Hessian of the Lagrangian
      LENIW = 3*NDV + NROWA + 2*NROWJ  ! Length of integer workspace
      LENW  = NDV*(2*NDV + NROWA +     ! Length of real workspace
     >        2*NROWJ + 20) + 11*NROWA +
     >        21*NROWJ

*     The above are for the dense-matrix NPOPT.  Sparse techniques can't predict
*     workspace precisely, so pad these for SNOPT:

      LENIW = LENIW + 50000
      LENW  = LENW  * 5

      SOLVE_CALLED = NCNLN > 0 ! CONFUN will calculate the objective & gradient

*     ---------------------------------------------------
*     Allocate (most of) the work-space that is dependent
*     on the number of design variables and constraints:
*     ---------------------------------------------------

      ALLOCATE (ILCON(NCLIN),  INLCON(NCNLN), JNLCON(NCNLN),
     >          ISTATE(NBNDS), IWVEC(LENIW))

      ALLOCATE (FFWD(NDV),     GRAD(NDV),    CVEC(NROWJ),   WVEC(LENW),
     >          IGROUP(NDV),   AITCH(NDV),   ONEOVERH(NDV),
     >          BL(NBNDS),     BU(NBNDS),    CLAMDA(NBNDS),
     >          BLLIN(NCLIN),  BULIN(NCLIN),
     >          CONLIN(NCLIN), TLCON(NCLIN),
     >          C(NCNLN),      CPRINT(NCNLN),
     >          XNLCON(NCNLN), SNLCON(NCNLN),
     >          AMAT(NCLIN,NDV), RMAT(NDV,NDV), CJAC(NROWJ,NDV))

      ALLOCATE (V(NDV), VSCALE(NDV), DOLO(NDV), DOUP(NDV),
     >          XBUMP(NDV), XBMIN(NDV), XBMAX(NDV), POWER(NDV),
     >          XBTYP(NDV), LCTYPE(NCLIN), NLCTYPE(NCNLN))

      RETURN

  999 STOP

      END SUBROUTINE RD_OPT

************************************************************************
*                                                                      *
      SUBROUTINE RD_TARGET (MODE)
*                                                                      *
*     Read a target distribution for optimizing an airfoil.            *
*     If MODE = 1, the targets are curvature radii;                    *
*     if MODE = 2, the targets are signed curvatures;                  *
*     if MODE = 3, the targets are airfoil coordinates.                *
*     Clockwise wrap-around order is expected in both cases, and one   *
*     or both surfaces may be present.  X/C = 0. must be present if    *
*     both surfaces are present.                                       *
*                                                                      *
*     Target radii or kappas:                                          *
*                                                                      *
*        Lines with '!' in column 1 are ignored.                       *
*                                                                      *
*        !     X/C       R [or KAPPA]                                  *
*          9.995447E-01  1.663058E+00                                  *
*          9.772109E-01  1.663058E+00                                  *
*          9.490933E-01  2.607097E+00                                  *
*          9.141979E-01  5.475490E+00                                  *
*           :             :                                            *
*                                                                      *
*     Target airfoil (compatible with PROFILE):                        *
*                                                                      *
*        Title                                                         *
*        N                                                             *
*        1.000000 -0.017500                                            *
*        0.999950 -0.017700                                            *
*        0.999900 -0.018010                                            *
*         :         :                                                  *
*                                                                      *
************************************************************************

      USE CONST_MOD
      USE GEOM_MOD

      IMPLICIT NONE

*     Argument:

      INTEGER, INTENT (IN) :: MODE  ! See above

*     Local constants:

      REAL, PARAMETER :: ZERO = 0.
      CHARACTER, PARAMETER :: FILENAME * 18 = 'nose_radius.target'

*     Local variables:

      INTEGER   I, IOS
      CHARACTER BUFFER * 80, TARGET_NAME * 5

*     Execution:

      OPEN (LTARGET, FILE=FILENAME, STATUS='OLD', IOSTAT=IOS)

      IF (IOS /= 0) THEN
         WRITE (LWRIT, '(/, 2A)') ' Trouble opening ', filename
         GO TO 999
      END IF

      SELECT CASE (MODE)

      CASE (1, 2) ! Target radii or kappas

*        Count the lines not starting with '!':

         NTARGET = 0

         DO ! Until EOF
            READ (LTARGET, '(A)', IOSTAT=IOS) BUFFER(1:1)
            IF (IOS < 0) EXIT

            IF (BUFFER(1:1) /= '!') NTARGET = NTARGET + 1
         END DO

         REWIND (LTARGET)

         ALLOCATE (S_TARGET(NTARGET), CRV_TARGET(NTARGET),
     >             R_TARGET(NTARGET), XOC_TARGET(NTARGET))

         NTARGET = 0

         DO ! Until EOF
            READ (LTARGET, '(A)', IOSTAT=IOS) BUFFER
            IF (IOS < 0) EXIT

            IF (BUFFER(1:1) /= '!') THEN
               NTARGET = NTARGET + 1
               READ (BUFFER, *, IOSTAT=IOS)
     >            XOC_TARGET(NTARGET), R_TARGET(NTARGET)
               IF (IOS /= 0) THEN
                  WRITE (LWRIT, '(/, A, I4, /, 2A)')
     >               ' Trouble reading target point:', NTARGET,
     >               ' Offending line :', BUFFER
                  GO TO 999
               END IF
            END IF
         END DO

         IF (MODE == 1) THEN
            TARGET_NAME = '    R'
         ELSE ! MODE = 2
            TARGET_NAME = 'KAPPA'
            CRV_TARGET  = R_TARGET ! 1:NTARGET; keep nomenclature clear
         END IF

      CASE (3) ! Target airfoil has slightly different format

         TARGET_NAME = '  Y/C'

         READ (LTARGET, '(A)')
         READ (LTARGET, *) NTARGET

         ALLOCATE (XOC_TARGET(NTARGET), YOC_TARGET(NTARGET))

         READ (LTARGET, *, IOSTAT=IOS)
     >      (XOC_TARGET(I), YOC_TARGET(I), I = 1, NTARGET)
         IF (IOS /= 0) THEN
            WRITE (LWRIT, '(/, (A))')
     >         ' Trouble reading the target airfoil.',
     >         ' The format should match that of nose_radius.dat,',
     >         ' though one surface may be omitted.'
            GO TO 999
         END IF

      END SELECT

      CLOSE (LTARGET)

*     Determine whether targets are on both surfaces or just one:

      NLTARG = 1

      DO I = 2, NTARGET
         IF (XOC_TARGET(I) < XOC_TARGET(I-1)) NLTARG = I
      END DO

*     NLTARG = 1       means targets are on the upper surface only
*     NLTARG = NTARGET means targets are on the lower surface only

      IF (NLTARG /= 1 .AND. NLTARG /= NTARGET) THEN
         IF (XOC_TARGET(NLTARG) /= ZERO) THEN ! Ambiguous situation
            WRITE (LWRIT, '(/, 2A)') ' Target data must include ',
     >         'X/C = 0 if both surfaces are present.'
            GO TO 999
         END IF
      END IF

      NUTARG = NTARGET - NLTARG + 1


*     Set up nonlinear weights for the target distribution.
*     Uh, oh!  The leading edge needs two weights!  Live with it for now.

      ALLOCATE (WT_TARGET(NTARGET))

      IF (NUTARG > 1) THEN  ! Upper surface

         WT_TARGET(NLTARG : NTARGET) = ADDTO(1) ! Offset

         IF (WBTYP(1)(1:4) == 'SCAL') THEN
            WT_TARGET(NLTARG : NTARGET) =
     >      WT_TARGET(NLTARG : NTARGET) * WPEAK(1)
         ELSE

            CALL SFEVAL (WBTYP(1)(1:4), WPOWER(1), WXCNTR(1),
     >                   WXMIN(1), WXMAX(1), WPEAK(1),
     >                   NLTARG, NTARGET, XOC_TARGET, WT_TARGET, LWRIT)
         END IF

      END IF

      IF (NLTARG > 1) THEN  ! Lower surface

         WT_TARGET(1 : NLTARG) = ADDTO(2) ! Offset

         IF (WBTYP(2)(1:4) == 'SCAL') THEN
            WT_TARGET(1 : NLTARG) =
     >      WT_TARGET(1 : NLTARG) * WPEAK(2)
         ELSE

            CALL SFEVAL (WBTYP(2)(1:4), WPOWER(2), WXCNTR(2),
     >                   WXMIN(2), WXMAX(2), WPEAK(2),
     >                   1, NLTARG, XOC_TARGET, WT_TARGET, LWRIT)
         END IF

      END IF

      WRITE (LWRIT, '(/, A, //, 3A, /)') ' Target distribution:',
     >   '   I            X/C          ', TARGET_NAME,
     >   '         WEIGHT'

      IF (MODE == 1) THEN ! Target radii
         WRITE (LWRIT, '(I4, 1P, 3E15.7)')
     >     (I, XOC_TARGET(I), R_TARGET(I), WT_TARGET(I), I = 1, NTARGET)
      ELSE IF (MODE == 2) THEN ! Target curvatures
         WRITE (LWRIT, '(I4, 1P, 3E15.7)')
     >     (I, XOC_TARGET(I),CRV_TARGET(I),WT_TARGET(I), I = 1, NTARGET)
      ELSE ! MODE = 3: Target airfoil
         WRITE (LWRIT, '(I4, 1P, 3E15.7)')
     >   (I, XOC_TARGET(I), YOC_TARGET(I), WT_TARGET(I), I = 1, NTARGET)
      END IF

*     Trick to simplify treatment of both surfaces together:

      IF (NLTARG /= 1) THEN  ! Some lower surface targets are present
         XOC_TARGET(1:NLTARG) = -XOC_TARGET(1:NLTARG)
      END IF

      RETURN

  999 STOP

      END SUBROUTINE RD_TARGET

************************************************************************
*                                                                      *
      SUBROUTINE RD_VARS
*                                                                      *
*     Read the design variables.                                       *
*                                                                      *
************************************************************************

      USE CONST_MOD
      USE DESIGN_VARS_MOD
      USE OPT_MOD

      IMPLICIT NONE

*     Local variables:

      INTEGER
     >   I, N,
     >   NSINE, NCOS, NEXPS, NFLAP, NSCALE, NSLAT, NTRAIL, NLEAD, NWAG
      LOGICAL
     >   ORDERED
      CHARACTER
     >   TYPE * 5

*     Execution:

      NSINE  = 0
      NCOS   = 0
      NEXPS  = 0
      NSCALE = 0
      NTRAIL = 0
      NLEAD  = 0
      NWAG   = 0
      NFLAP  = 0
      NSLAT  = 0

      READ (LREAD,*)
      READ (LREAD,*)
      READ (LREAD,*)

      DO I = 1, NDV

         READ (LREAD,*) N, TYPE,  DOUP(I), DOLO(I),
     >                  POWER(I), XBUMP(I), XBMIN(I), XBMAX(I),
     >                  V(I), VSCALE(I), AITCH(I), BL(I), BU(I)

         IF (XBUMP(I) < XBMIN(I) .OR. XBUMP(I) > XBMAX(I)) THEN
            N = I
            GO TO 850
         END IF

         IF (BL(I) == -999.) BL(I) = -BIGBND
         IF (BU(I) ==  999.) BU(I) =  BIGBND

         ONEOVERH(I) = 1./ AITCH(I)
         XBTYP(I)    = TYPE

         IF      (TYPE == 'SCALE') THEN
            NSCALE = NSCALE + 1
            IGROUP(I) = 1
         ELSE IF (TYPE == 'SIN  ') THEN
            NSINE  = NSINE  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'SINF ') THEN
            NSINE  = NSINE  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'SIN1 ') THEN
            NSINE  = NSINE  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'SIN2 ') THEN
            NSINE  = NSINE  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'SIN3 ') THEN
            NSINE  = NSINE  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'SIN4 ') THEN
            NSINE  = NSINE  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'COSL ') THEN
            NCOS  = NCOS  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'COSR ') THEN
            NCOS  = NCOS  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'LCOS ') THEN
            NCOS  = NCOS  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'RCOS ') THEN
            NCOS  = NCOS  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'EXP  ') THEN
            NEXPS  = NEXPS  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'TRL  ') THEN
            NTRAIL = NTRAIL + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'LED  ') THEN
            NLEAD  = NLEAD  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'LEAD ') THEN
            NLEAD  = NLEAD  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE(1:3) == 'WAG') THEN
            TYPE = 'WAG  ' ! Make certain for PERTURB
            NWAG = NWAG  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'FLAP ') THEN
            NFLAP  = NFLAP  + 1
            IGROUP(I) = 3
         ELSE IF (TYPE == 'SLAT ') THEN
            NSLAT  = NSLAT  + 1
            IGROUP(I) = 3
         ELSE
            GO TO 800  ! Unknown variable type
         END IF

      END DO

      NBUMPS = NSCALE + NSINE + NCOS + NEXPS + NTRAIL + NLEAD + NWAG +
     >         NFLAP + NSLAT

      IF (NBUMPS /= NDV) THEN
         WRITE (LWRIT,*) 'Total # design variables did not match NDV.'
         GO TO 900
      END IF

      ORDERED = .TRUE.
      DO I = 2, NDV
         IF (IGROUP(I) < IGROUP(I-1)) ORDERED = .FALSE.
      END DO

      DEALLOCATE (IGROUP)

      IF (.NOT. ORDERED) THEN
         WRITE (LWRIT, '(/, 1X, A, A, /, (4X, A))')
     >      'RD_VARS:  The optimization variables must be grouped ',
     >      'in the following order:',
     >      'Simple Y scaling',
     >      'Hicks-Henne-type sine bumps and Wagner functions',
     >      'Flap and slat variables'
         GO TO 900
      END IF

      GO TO 999

  800 WRITE (LWRIT, '(/, A, I5, 2X, A)')
     >   ' RD_VARS: Unknown design variable type: ', I, TYPE
      GO TO 900

  850 WRITE (LWRIT, '(/, A, I5)')
     >   ' RD_VARS: X center is outside [min, max]. Variable:', N

  900 STOP

  999 RETURN

      END SUBROUTINE RD_VARS

************************************************************************
*                                                                      *
      SUBROUTINE SAVE_CRV
*                                                                      *
*     Save optimized [radius of] curvature results.                    *
*                                                                      *
************************************************************************

*     Global variables:

      USE CONST_MOD
      USE GEOM_MOD

      IMPLICIT NONE

      INTEGER I

      WRITE (LRADOPT, '(A)') '!     X/C           R'
      WRITE (LRADOPT, '(1P, 2E14.6)') (X_PRTB(I),  R_PRTB(I), I = 1, NW)

      WRITE (LCRVOPT, '(A)') '!     X/C         KAPPA'
      WRITE (LCRVOPT, '(1P, 2E14.6)') (X_PRTB(I),CRV_PRTB(I), I = 1, NW)

      END SUBROUTINE SAVE_CRV

************************************************************************
*                                                                      *
      SUBROUTINE SETLCON
*                                                                      *
*     SETLCON sets up the linear constraint matrix.                    *
*                                                                      *
************************************************************************

*     Global variables:

      USE DESIGN_VARS_MOD
      USE CONST_MOD
      USE OPT_MOD
      USE GEOM_MOD

      IMPLICIT NONE

*     Local constants:

      REAL, PARAMETER :: ZERO = 0.

*     Local variables:

      REAL    :: BOUNDL, BOUNDU, DELTAV, T, ZW
      INTEGER :: KW, M, N

*     Execution:

*     Zero the linear constraint matrix - only some of the variables
*     are applicable.

      AMAT = ZERO ! (1:NCLIN,1:NDV)

*     Determine the nonzero matrix elements:

      DO M = 1, NCLIN

         KW     = ILCON(M)
         ZW     = TLCON(M)
         BOUNDL = BLLIN(M)
         BOUNDU = BULIN(M)

!        < No linear constraints yet >

         BU(M+NDV) = BOUNDU
         BL(M+NDV) = BOUNDL

      END DO ! Next linear constraint

      END SUBROUTINE SETLCON

************************************************************************
*                                                                      *
      SUBROUTINE SETUP_NPOPT
*                                                                      *
************************************************************************

      USE CONST_MOD
      USE OPT_MOD

      IMPLICIT NONE

*     Namelist to override NPOPT's options file defaults if desired:

      NAMELIST /NPOPTIONS/
     >   LEVELVER, MAJORPL, MINORPL, MINORIL, NPRFREQ,
     >   LEVELDER, DIFFERENCE_INTERVAL, PENPARAM, STEPLIM,
     >   TOLLIN, TOLNLIN, TOLOPT

*     Execution:

      DIFFERENCE_INTERVAL = -1.  ! If > 0, else NPOPT calculates one
                                 ! for each variable
      STEPLIM  = 2.              ! 2 means allow a 200% change in the variables
      PENPARAM = 0.
      LEVELVER = -1
      LEVELDER = 0
      TOLLIN   = 1.E-6
      TOLNLIN  = TOLOPT
      TOLOPT   = 1.E-5
      MAJORPL  = 10
      MINORPL  = 0
      MINORIL  = 1000
      NPRFREQ  = 100
      NFTOTL   = 0               ! Incremented in CONFUN and/or OBJFUN

*     Look for an optional namelist at the end of standard input:

      READ (LREAD, NPOPTIONS, ERR=100, END=100)

  100 CLOSE (LREAD)

*     -------------------------------------------------
*     Run-time set-up of the optional inputs for NPOPT:
*     -------------------------------------------------

      OPEN (UNIT=LNPOPT, FILE='nose_radius.nps', STATUS='UNKNOWN')

      WRITE (LNPOPT, '(A)') 'Begin'

      WRITE (LNPOPT, '(A)')
     >   'Cold Start',
     >   'Nonderivative line search',
     >   'Hessian       Full memory'
      WRITE (LNPOPT, '(A, I6)')
     >   'Major Iteration Limit          ', NITMAX,
     >   'Minor Iteration Limit          ', MINORIL,
     >   'Major Print Level              ', MAJORPL,
     >   'Minor Print Level              ', MINORPL,
     >   'Print frequency                ', NPRFREQ,
     >   'Verify level                   ', LEVELVER,
     >   'Derivative Level               ', LEVELDER

      IF (DIFFERENCE_INTERVAL > 0.)       ! Else NPOPT finds one for each var.
     >    WRITE (LNPOPT, '(A, E12.3)')    ! if Derivative Level = 0
     >   'Difference Interval            ', DIFFERENCE_INTERVAL

      WRITE (LNPOPT, '(A, E12.3)')
     >   'Function Precision             ', EPSOBJ,
     >   'Infinite Bound                 ', BIGBND,
     >   'Infinite Step                  ', STEPMX,
     >   'Major Step Limit               ', STEPLIM,
     >   'Minor Feasibility Tolerance    ', TOLLIN,
     >   'Major Feasibility Tolerance    ', TOLNLIN,
     >   'Optimality Tolerance           ', TOLOPT,
     >   'Linesearch Tolerance           ', ETA,
     >   'Penalty parameter              ', PENPARAM,
     >   'End'

      REWIND (LNPOPT)

*     NPOPT writes to standard output (unit 6 = LNPOUT).
*     Suppress optional screen output:

      CALL NPINIT (LNPOUT,      0, IWVEC, LENIW, WVEC, LENW)
      CALL NPSPEC (LNPOPT, INFORM, IWVEC, LENIW, WVEC, LENW)

      CLOSE (LNPOPT)

      WRITE (LWRIT, '(/, A)')

      IF (INFORM /= 0) THEN
         WRITE (LWRIT, '(/, 2A)') ' SET_NPOPT: ',
     >      'Trouble with the run-time nose_radius.nps options file.'
      END IF

      RETURN

      END SUBROUTINE SETUP_NPOPT

********************************************************************************
*
      SUBROUTINE SFEVAL (NAME, EX, XC, XA, XB, YMULT, I1, I2, X, Y, LUN)
*
*        SFEVAL evaluates a specified shape function at X (I1:I2) and
*     adds the result to the corresponding Y.  X need not be normalized
*     because the option to apply the shape function to just part of the
*     data range, [XA, XB], forces normalization locally.  XC, however,
*     refers to the normalized sub-range, so it must be normalized.
*
*        These are the Hicks-Henne shape functions used for aerodynamic
*     design.  SFEVAL is an adaptation of BEVAL from program PROFILE in
*     the form most convenient for design code SYN87, where a modular
*     form is needed for precise set-up of linear constraints as well
*     as for the geometry perturbations.
*
*        Distinguishing wing and body applications should be done at
*     the calling level to avoid too many similar cases here.
*
*     06/21/96    DAS    Adaptation of BEVAL and PERTURB routines.
*     08/22/96     "     Added SINF. Apart from the SIN* functions,
*                        testing just 3 characters avoids problems.
*     08/24/96     "     LED & LEA now both work for LEADing;
*                        TRL & TRA  "    "    "   "  TRAILing.
*     12/06/96  DAS/JJR  Added COSL & COSR for body keel & crown,
*                        and inverted forms LCOS & RCOS.
*     06/03/97    DAS    Converting flap & slap angles from degrees to
*                        radians was using 180./PI instead of PI/180.
*     08/09/97     "     FLAP/SLAT now expect hinge X in same units as
*                        X(*) and Y(*) (not fraction of XB - XA).
*     08/10/99     "     Added SIN3 (fore/aft symmetry via SIN1 on [0,1])
*                        and   SIN4 (  "   "   "   "   "   "   " each half).
*     11/03/03     "     Generalized the COS* and *COS functions to extend
*                        the peak value any distance.  E.g., for COSL,
*                        XC = XA gives the original form, else COSL = 1.0
*                        for X in [XA, XC], going to 0.0 at XB.
*     11/06/03     "     XC units should match those of XA, XB for the
*                        generalized COS* and *COS functions.
*
*     Author:  David Saunders, Sterling Software/NASA Ames, Moffett Field, CA
*
********************************************************************************

      IMPLICIT NONE

*     Arguments:

      CHARACTER, INTENT (IN) ::
     >   NAME * (*)          ! Shape function name (upper case;
                             ! only characters 1:4 are checked here,
                             ! or 1:3 where that clearly suffices
      REAL, INTENT (IN) ::
     >   EX,                 ! Exponent, usually affecting "width";
                             ! for 'WAG'[ner] fn. N pass REAL (N)
     >   XC,                 ! "Center" X in [0, 1] where fn. peaks;
                             ! XC is in [XA, XB] for *COS and COS*
     >   XA, XB,             ! X sub-range in the same units as X (*)
     >   YMULT               ! Shape function multiplier

      INTEGER, INTENT (IN) ::
     >   I1, I2              ! Index range of X & Y to be treated

      REAL, INTENT (IN) ::
     >   X (*)               ! Abscissas, not necessarily in [0, 1]

      REAL, INTENT (INOUT) ::
     >   Y (*)               ! Ordinates being perturbed

      INTEGER, INTENT (IN) ::
     >   LUN                 ! Logical unit for error messages

*     Local constants:

      REAL, PARAMETER ::
     >   HALF = 0.5, ONE = 1., TWO = 2., ZERO = 0.

*     Local variables:

      INTEGER
     >   I

      REAL
     >   AEXP, BEXP, CENTER, CENTER2, DEGRAD, EPSMCH, PT5LOG,
     >   ONEMC, PI, PIBY2, POWER, POWERL, POWERR, RN, RNM1, RPI,
     >   RRANGE, SCALE, TANGNT, THETA, WIDTH, X0, XI, XSCALE

      CHARACTER
     >   KEY4 * 4, KEY3 * 3

      LOGICAL ::
     >   FIRST = .TRUE.

*     Storage:

      SAVE
     >   FIRST, DEGRAD, PI, PIBY2, RPI, PT5LOG

*     Statement functions:

      REAL
     >   EBUMP, SBUMP, PWR, WDTH, XNRM

      EBUMP (WDTH, PWR, XNRM) =
     >   XNRM ** PWR * (ONE - XNRM) * EXP (-WDTH * XNRM)

      SBUMP (WDTH, PWR, XNRM) = (SIN (PI * (XNRM ** PWR))) ** WDTH

*     Execution:

      IF (FIRST) THEN

         FIRST = .FALSE.

*        Protect sine bump exponentiation from slightly negative arguments
*        by calculating a value for PI such that SIN (PI) > 0 (barely).
*        This allows the limiting of X to [0, 1] that is necessary for
*        the sub-range capability to cover the sine protection too.

         EPSMCH = EPSILON (EPSMCH)

         PI = ASIN (ONE) * TWO

         DO I = 1, 999 ! 6 iterations observed for DEC AXP    (32 bits)
                       ! 5     "         "      "  CRAY C90   (64 bits)
                       ! 1     "         "      "  SGI R10000 ( "   " )
***         write (6,*) 'I   PI  SIN (PI): ', i, pi, sin (pi)
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

*     Minimize argument references in loops:

      WIDTH  = EX
      CENTER = XC ! Normalized
      X0     = XA
      RRANGE = ONE / (XB - X0)
      SCALE  = YMULT
      KEY4   = NAME (1:4) ! Avoid comparison of different-length strings
      KEY3   = KEY4 (1:3)

*     Treat the cases in ~ the most likely order they'll be encountered:

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

*        Assume [XA, XB] = [0, 1] and CENTER in (0, 0.5]:

         POWER = PT5LOG / LOG (CENTER)

         DO I = I1, I2
            XI = MAX (ZERO, MIN (X (I), ONE))
            IF (XI > HALF) XI = ONE - XI
            Y (I) = SCALE * SBUMP (WIDTH, POWER, XI) + Y (I)
         END DO

      ELSE IF (KEY4 == 'SIN4') THEN ! Fore/aft symmetry via SIN1 on each half

*        Assume [XA, XB] = [0, 0.5] and CENTER in (0, 1).  Left half:

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

*        Now the [0.5, 0.1] half with center <-- 1 - center

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

*        SIN is used instead of COS because SIN (PI) is protected above.
*        This version gives the peak value for X (I) in [XA, XC].

         XSCALE = PIBY2 / (XB - XC)
         X0 = TWO * XC - XB

         DO I = I1, I2
            XI = MAX (PIBY2, MIN ((X (I) - X0) * XSCALE, PI))
            Y (I) = SCALE * (SIN (XI)) ** WIDTH + Y (I)
         END DO

      ELSE IF (KEY4 == 'COSR') THEN ! Half a cosine (or sine), peak at right

*        SIN is used instead of COS for consistency with COSL.
*        This version gives the peak value for X (I) in [XC, XB].

         XSCALE = PIBY2 / (XC - XA)

         DO I = I1, I2
            XI = MAX (ZERO, MIN ((X (I) - XA) * XSCALE, PIBY2))
            Y (I) = SCALE * (SIN (XI)) ** WIDTH + Y (I)
         END DO

      ELSE IF (KEY4 == 'LCOS') THEN ! Inverted half (co)sine, peak at left

         XSCALE = PIBY2 / (XB - XC)

         DO I = I1, I2
            XI = MAX (ZERO, MIN ((X (I) - XC) * XSCALE, PIBY2))
            Y (I) = SCALE * (ONE - (SIN (XI)) ** WIDTH) + Y (I)
         END DO

      ELSE IF (KEY4 == 'RCOS') THEN ! Inverted half (co)sine, peak at right

         XSCALE = PIBY2 / (XC - XA)
         X0 = TWO * XA - XC

         DO I = I1, I2
            XI = MAX (PIBY2, MIN ((X (I) - X0) * XSCALE, PI))
            Y (I) = SCALE * (ONE - (SIN (XI)) ** WIDTH) + Y (I)
         END DO

      ELSE IF (KEY3 == 'WAG') THEN ! Wagner function # N

*        Reference:  Ramamoorthy, P. and Padmavathi, K.  "Airfoil Design
*        by Optimization" in J. Aircraft, Vol. 14 (Feb. 1977), 219-221.

         IF (WIDTH == ONE) THEN ! N = 1
            DO I = I1, I2
               XI = MAX (ZERO, MIN ((X (I) - X0) * RRANGE, ONE))
               THETA = TWO * ASIN (SQRT (XI))
               Y (I) = ((THETA + SIN (THETA)) * RPI -
     >                  (SIN (HALF * THETA)) ** 2) * SCALE + Y (I)
            END DO
         ELSE
            RN   = ONE / WIDTH ! 1 / N
            RNM1 = WIDTH - ONE ! N - 1
            DO I = I1, I2
               XI = MAX (ZERO, MIN ((X (I) - X0) * RRANGE, ONE))
               THETA = TWO * ASIN (SQRT (XI))
               Y (I) = ((SIN (WIDTH * THETA) * RN  +
     >                  SIN (RNM1 * THETA)) * RPI) * SCALE + Y (I)
            END DO
         END IF

      ELSE IF (KEY3 == 'FLA') THEN

*        "FLAP"-like function (shearing transformation only):
*        XC is the hinge point in the same frame as X & Y;
*        YMULT is the flap angle in DEGREEs (positive is "down").

         TANGNT = TAN (DEGRAD * SCALE)
         DO I = I1, I2
            XI = MAX (ZERO, X (I) - CENTER)
            Y (I) = Y (I) - XI * TANGNT
         END DO

      ELSE IF (KEY3 == 'SLA') THEN

*        "SLAT"-like function (shearing transformation only):
*        XC is the hinge point in the same frame as X & Y;
*        YMULT is the slat angle in DEGREEs (positive is "down").

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

*        Commonly used in conjunction with Wagner functions.

         DO I = I1, I2
            XI = MAX (ZERO, MIN ((X (I) - X0) * RRANGE, ONE))
            Y (I) = SCALE * XI + Y (I)
         END DO

      ELSE IF (KEY3 == 'SCA') THEN

*        Simple scaling - X(I) are probably ordinates.  Note that the
*        effect is arranged to be additive, not multiplicative, so that
*        scaling can be treated just as for the other functions when it
*        is included in a group for perturbing purposes.

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

********************************************************************************
*
      SUBROUTINE SOLVE (NDVPAR, VPAR, OBJPAR, NCALL, FAIL)
*
*     SOLVE is the objective function routine expected by the unconstrained
*     optimizer QNMDIF2.  It is also called by the OBJFUN and CONFUN routines
*     expected by the constrained optimizer NPOPT.
*
*     NCALL usage:
*
*     -13 = first call to SOLVE from CONFUN
*           (precedes first call to OBJFUN)
*     -12 = constraint evaluations following a line search
*     -10 = constraint evaluations within a line search
*      -4 = final call to SOLVE from main program
*      -3 = first call to SOLVE from OBJFUN
*      -2 = objective evaluation after a line search
*       0 = objective evaluation within a line search
*      >0 = gradient element NCALL evaluation
*           (constraints, objective function, or both)
*
********************************************************************************

*     Global variables:

      USE CONST_MOD
      USE DESIGN_VARS_MOD
      USE OPT_MOD

      IMPLICIT NONE

*     Arguments:

      INTEGER, INTENT (IN)    :: NDVPAR  ! Number of design variables
      INTEGER, INTENT (INOUT) :: NCALL   ! See usage above
      REAL,    INTENT (OUT)   :: OBJPAR  ! Objective function value
      REAL,    INTENT (INOUT) :: VPAR(*) ! * because NDVPAR may be zero.
      LOGICAL, INTENT (OUT)   :: FAIL    ! T means the function value is bad

*     Local constants:

      INTEGER,   PARAMETER :: ICALL = 1
      REAL,      PARAMETER :: ZERO  = 0.
      LOGICAL,   PARAMETER :: PRINT = .TRUE.
      CHARACTER, PARAMETER :: SUBNAME*5 = 'SOLVE'

*     Local variables:

      INTEGER    :: I, MCALL, N
      REAL(KIND=4), SAVE :: CPU1, CPU2, CPUGRAD1, TIME, TIME1, TIME2
      LOGICAL    :: FIRST, OUTPUTS

*     Execution:

*     -------------------------------------------------------------------
*     Clarifications:
*
*     (1) SOLVE_CALLED = T means OBJFUN does not call SOLVE at all.
*     (2) SOLVE_CALLED = T means CONFUN calls must do the outputs that
*         would otherwise be done by calls from OBJFUN (or main).
*     (3) NCALL = -13, -12, -10  correspond to CONFUN calls equivalent to
*         NCALL =  -3,  -2,   0  calls from OBJFUN.
*     -------------------------------------------------------------------

      FAIL  = .FALSE.
      MCALL = NCALL

      IF (MCALL < 0) MCALL = MOD (MCALL, 10)

      FIRST   = MCALL == -3
      OUTPUTS = MCALL < 0

      CALL SECOND (TIME1)

      IF (MCALL == 1) THEN ! Start of a gradient calculation
         CPUGRAD1 = TIME1
      ELSE IF (MCALL == -2) THEN
         NDITER = NDITER + 1
      END IF

      IF (OUTPUTS) THEN
         WRITE (LWRIT, '(/, A, I4, /)')
     >      ' DESIGN VARIABLES AT DESIGN ITERATION', NDITER
         WRITE (LWRIT, '(1P, 5(I5, E21.13))')
     >      (I, VPAR(I), I = 1, NDVPAR)
      END IF

*     ---------------------------
*     Apply any design variables:
*     ---------------------------

      CALL PERTURB (VPAR, MCALL)

*     -------------------
*     Objective function:
*     -------------------

      CALL OBJECTIVE (NDVPAR, VPAR, OBJPAR, MCALL)

      IF (FIRST) OBJINI = OBJPAR

*     ---------------
*     Output results?
*     ---------------

      IF (OUTPUTS) THEN

         IF (NCLIN > 0) CALL CHKLCON ()      ! Review the linear constraints

         IF (NCNLN > 0) CALL NLCON (PRINT, CPRINT) ! & nonlinear constraints

         CALL SUMMARY (OBJPAR, MCALL)        ! Tabulate

      END IF

***   IF      (MCALL == -3 ) THEN
***      WRITE (LWRIT, 1080) 'Initial'
***   ELSE IF (MCALL == -2 ) THEN
***      WRITE (LWRIT, 1080) 'End-of-line-search'
***   ELSE IF (MCALL == -1 ) THEN
***      WRITE (LWRIT, 1080) 'Local search'
***   ELSE IF (MCALL ==  0 ) THEN
***      WRITE (LWRIT, 1080) 'Line search'
***   END IF


***   IF (MCALL == NDVPAR) THEN
***      IF (.NOT. SOLVE_CALLED) THEN ! SOLVE_CALLED means CONFUN does it
***         CALL SECOND (TIME2)
***         CALL CPUTIME (CPUGRAD1, TIME2, SUBNAME,
***  >                    'to calculate gradient', LWRIT)
***      END IF
***   END IF

***   IF (MCALL < 0) THEN
***      CALL SECOND (TIME2)
***      CALL CPUTIME (CPU0, TIME2, SUBNAME,
***  >                 'for this run so far', LWRIT)
***   END IF

      GO TO 999

*     Termination:

  900 STOP

  999 CALL FLUSH (LWRIT) ! Flush output buffer

      RETURN

*     Formats:

 1080 FORMAT (/, ' SOLVE: ', A, ' evaluation.')
 1300 FORMAT (/, ' SOLVE: ', (1X, A))

      END SUBROUTINE SOLVE

************************************************************************
*                                                                      *
      SUBROUTINE SUMMARY (OBJ, MCALL)
*                                                                      *
*     Summarize the initial and optimized results.  Save geometry?     *
*                                                                      *
************************************************************************

*     Global variables:

      USE CONST_MOD
      USE GEOM_MOD
      USE OPT_MOD

      IMPLICIT NONE

*     Arguments:

      REAL    OBJ
      INTEGER MCALL ! -4 means SOLVE is being called at the end of a run

*     Local variables:

      INTEGER I

*     Execution:

      WRITE (LWRIT, 1200)
     >   NDITER,
     >   OBJINI, OBJ,
     >   OBJINI, OBJ   ! Put something else here later

*     Save best result at the end of a design run (only):

      IF (MCALL == -4) THEN

         WRITE (LGEOPT, '(A)') 'Optimized airfoil'
         WRITE (LGEOPT, '(I6)') NW
         WRITE (LGEOPT, '(1P, 2E15.7)')
     >      (X_PRTB(I), Y_PRTB(I), I = 1, NW)

         CLOSE (LGEOPT)

      END IF

      RETURN

 1200 FORMAT (/,
     >   ' AIRFOIL PARAMETERS AT DESIGN ITERATION', I4, //, 42X,
     >   ' Initial', 12X, 'Current', 1P, /,
     >   ' Objective function            ', 2E19.10, /,
     >   ' What else here?               ', 2E19.10)

      END SUBROUTINE SUMMARY
