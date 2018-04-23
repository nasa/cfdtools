!     CH_GRID Modules:

      MODULE DIMENS

         INTEGER, PARAMETER ::
     >      MXIWING = 265, ! Input geometry points per wing section wraparound
     >      MXKWING = 80,  ! Input wing sections + ZCRANK insertions
     >      MXIBODY = 210, ! Input geometry fuselage stations, nose to tail
     >      MXJBODY = 151, ! Input geometry points per (half) fuselage station
     >      MXIMAX  = 265, ! MAX (MXIWING, MXIBODY)
     >      MXKMAX  = 151, ! MAX (MXKWING, MXJBODY)
     >      MXFLOOR = 5,   ! Max. # points defining the cabin floor segments
     >      MXRADII = 6,   ! Max. # radial pts. defining a cabin section polygon
     >      MXCRANK = 6,   ! The code needs to append the tip station
     >      NJBODY  = 151, ! # uniform points for each regularized body section
     >      NREFN   = 101, ! May be used by Z0PLANE: <= IL assumed for ELLIP2D
     >      NCASE2D = 15,  ! Elliptic smoothing cases
     >      NCASE3D = 5

      END MODULE DIMENS

      MODULE CONSTS

         REAL DEG2RAD, EPSMCH, LOGPT5, PI, RAD2DEG, RPI

      END MODULE CONSTS

      MODULE DMAX

         REAL XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, RADIUS

      END MODULE DMAX

      MODULE ELLIPTIC

         USE DIMENS

         INTEGER
     >      ITFLOAT, ITFREEZE, JLAYER,
     >      ITMAX2D(NCASE2D), ITMAX3D(NCASE3D)

         REAL
     >      DMAX2D, DMAX3D,   FGLIM2D, FGLIM3D, URFG2D, URFG3D, URFLOAT,
     >      CONV2D(NCASE2D),  CONV3D(NCASE3D),  CONVMIN(NCASE2D),
     >      EXPI12D(NCASE2D), EXPI22D(NCASE2D), EXPJ12D(NCASE2D),
     >      EXPJ22D(NCASE2D),
     >      EXPI13D(NCASE3D), EXPI23D(NCASE3D), EXPJ13D(NCASE3D),
     >      EXPJ23D(NCASE3D), EXPK13D(NCASE3D), EXPK23D(NCASE3D),
     >      OMG2D(NCASE2D),   OMG3D(NCASE3D),
     >      POWERI(NCASE2D),  POWERJ(NCASE2D)
         LOGICAL
     >      PRINT2D, PRINT3D
         CHARACTER
     >      BGMODE2D(NCASE2D)*4,
     >      BGPHI(NCASE3D)*4,    BGPSI(NCASE3D)*4,    BGOMG(NCASE3D)*4,
     >      FGMODE2D(NCASE2D)*4, FGMODE3D(NCASE3D)*6, SPMODE(NCASE2D)*4

      END MODULE ELLIPTIC

      MODULE FACES

         REAL, ALLOCATABLE ::
     >      DFACEI(:,:,:,:,:), DFACEJ(:,:,:,:,:), DFACEK(:,:,:,:,:)

      END MODULE FACES

      MODULE FUEL

         USE DIMENS
         INTEGER
     >      NFUEL, KFUEL1(MXKWING), KFUEL2(MXKWING)
         REAL
     >      XAFUEL1(MXKWING), XBFUEL1(MXKWING),
     >      XAFUEL2(MXKWING), XBFUEL2(MXKWING)

      END MODULE FUEL

      MODULE FUSEL

         USE DIMENS
         INTEGER
     >      NF(MXIBODY), NFSTN, NENVL, NENVU
         REAL
     >      FNF, XF(MXIBODY), YF(NJBODY, MXIBODY), ZF(NJBODY, MXIBODY)

      END MODULE FUSEL

      MODULE FUSEL2

         USE DIMENS

         INTEGER
     >      NRADII(MXFLOOR-1), NFLOOR
         REAL
     >      YFINIT(MXJBODY,MXIBODY),  ZFINIT(MXJBODY,MXIBODY),
     >      YFPRTB(MXJBODY,MXIBODY),  ZFPRTB(MXJBODY,MXIBODY),
     >      YFUDGE(MXJBODY,MXIBODY),  ZFUDGE(MXJBODY,MXIBODY),
     >      RBODY(MXRADII,MXFLOOR-1), THBODY(MXRADII,MXFLOOR-1),
     >      XFLOORI(MXFLOOR), YFLOORI(MXFLOOR),
     >      XFLOORF(MXFLOOR), YFLOORF(MXFLOOR),
     >      TCABIN(NJBODY), YCABIN(NJBODY), ZCABIN(NJBODY)

      END MODULE FUSEL2

      MODULE INDEX

         INTEGER
     >      NBODY, NWING, NHALF, KWING1, KWING2, IBOD1, IBOD2, IDEGBOD

      END MODULE INDEX

      MODULE LIM

         INTEGER IL, JL, KL, ITL, IBL, ILE, ITU, IBU, JCR, KTIP,
     >           NKREGRID, IFACE(5), KFACE(3), KREGRID(9)

      END MODULE LIM

      MODULE LOGIC

         LOGICAL
     >      WINGONLY, GRIDONLY, READGRID, PERTURBONLY, PARABODY,
     >      CHECKWARP, REDISTRIBUTE, RESMOOTH, OPTIMIZING

      END MODULE LOGIC

      MODULE LUNS

         INTEGER, PARAMETER ::
     >      LGEOM  = 1,
     >      IREAD  = 5,
     >      IWRIT  = 6,
     >      LWING  = 2,
     >      LUNBOD = 3,
     >      LUNXYZ = 4,
     >      LUNXY0 = 7

      END MODULE LUNS

      MODULE OPT

         USE DIMENS

         INTEGER
     >      NDV,    NBOUND, NCLIN,  NCNLN,  NBNDS,  IYSMOO,
     >      NTWING, NCHCKS, NCHCK1, NCHCK2, NTWIST, NTHICK,
     >      NXLEAD, NYLEAD, NZLEAD, NCHORD, NIROOT,
     >      NCMBRF, NSPANF, NAREAF, NXFLR,  NYFLR,  NTKINK,
     >      NSINE,  NCOS,   NEXPS,  NTRAIL, NLEAD,  NWAG,   NSINVC,
     >      NSINEP, NEXPSP, NTRAIP, NLEADP, NFLAP,  NSLAT,  NTOTFS,
     >      NSINFC, NCOSFC, NEXPFC, NLEDFC, NTRLFC,
     >      NSINFS, NCOSFS, NEXPFS, NLEDFS, NTRLFS,
     >      NSINFA, NCOSFA, NEXPFA, NLEDFA, NTRLFA, NTOTFU,
     >      NUPPER(MXKWING),  NLOWER(MXKWING),  NTOTAL(MXKWING),
     >      NUPPER0(MXKWING), NLOWER0(MXKWING), NTOTAL0(MXKWING)

         INTEGER, ALLOCATABLE ::
     >      KMIN(:),  KCEN(:), KMAX(:), IGROUP(:),
     >      KLCON(:), INLCON(:), JNLCON(:)

         REAL
     >      XINIT(MXIWING,MXKWING), XPRTB(MXIWING,MXKWING),
     >      YINIT(MXIWING,MXKWING), YPRTB(MXIWING,MXKWING),
     >      ZINIT(MXIWING,MXKWING), ZPRTB(MXIWING,MXKWING),
     >      XLEADI(MXKWING), YLEADI(MXKWING), ZLEADI(MXKWING),
     >      XLEADW(MXKWING), YLEADW(MXKWING), ZLEADW(MXKWING),
     >      AINIT(MXKWING),  CINIT(MXKWING),  TINIT(MXKWING),
     >      AWORK(MXKWING),  CWORK(MXKWING),  TWORK(MXKWING),
     >      THMAX(MXKWING),  XTHMAX(MXKWING), FAREA(MXKWING),
     >      THBUMP,  EPSY,   EXPY,   WINGVL,  FUELVL,
     >      ALPHAW,  ALPHAI, ALMAXF, ALPHAF,  TMAXF

         REAL, ALLOCATABLE ::
     >      V(:),
     >      WBUMP(:), XBUMP(:), XWMIN(:), XWMAX(:),
     >      PBUMP(:), UBUMP(:), UBMIN(:), UBMAX(:),
     >      DOUP(:),  DOLO(:),  VSCALE(:), AITCH(:), BL(:), BU(:),
     >      XLCON(:), BLLIN(:), BULIN(:), CONLIN(:), C(:),
     >      XNLCON(:), SNLCON(:), AMAT(:,:)

         CHARACTER, ALLOCATABLE ::
     >      VTYPE(:)*6, UTYPE(:)*4, LCTYPE(:)*5, NLCTYPE(:)*6

      END MODULE OPT

      MODULE OPTIONS

         INTEGER
     >      NNACEL, NNFLOOR, NNGRID, NNROOT, NNTIP, WINGIN, WINGOUT

      END MODULE OPTIONS

      MODULE OUTR

         REAL, ALLOCATABLE ::
     >      XLBACK(:,:), YLBACK(:,:), ZLBACK(:,:), XOUT(:),
     >      XUBACK(:,:), YUBACK(:,:), ZUBACK(:,:), YOUT(:),
     >      XCROWN(:),   YCROWN(:)

      END MODULE OUTR

      MODULE SIZES

         REAL ZSPAN, ASPECT, SREF, SREF2, CREF, XREF, YREF, ZREF

      END MODULE SIZES

      MODULE SPACING

         USE DIMENS

         INTEGER
     >      KTOUTER, NBLAYER, NBLAYEREU, NBLAYERNS
         REAL
     >      CRDSPL,    CRDSPQ,   CRDSPS,   CRDSPC,   SPNSPC,   SPNSPS,
     >      OSPNOSE,   OSPTAIL,  OSPSYM,   OSPTIP,
     >      SWEEP,     DIHED,    BSWEEP,   TANSWEEP,
     >      TEROOT,    TETIP,    TEKMAX,   RBLAYER,
     >      D1CROWN,   D1NOSE,   D1TAIL,   DWTAIL,   DCTAIL,   D2JL,
     >      D1CROWNEU, D1NOSEEU, D1TAILEU, DWTAILEU, DCTAILEU, D2JLEU,
     >      D1ROOTEU,  D2ROOTEU, D3ROOTEU, D1TIPEU,  D2TIPEU,  D3TIPEU,
     >      D1KMAXEU,  D2KMAXEU, D3KMAXEU, RBLAYEREU,
     >      D1CROWNNS, D1NOSENS, D1TAILNS, DWTAILNS, DCTAILNS, D2JLNS,
     >      D1ROOTNS,  D2ROOTNS, D3ROOTNS, D1TIPNS,  D2TIPNS,  D3TIPNS,
     >      D1KMAXNS,  D2KMAXNS, D3KMAXNS, RBLAYERNS
         REAL, ALLOCATABLE ::
     >      DRADIALEU(:,:), DRADIALNS(:,:), ARCOUTER(:),
     >      DRADIAL(:,:), TRADIAL(:), XRADIAL(:), YRADIAL(:), ZRADIAL(:)

      END MODULE SPACING

      MODULE TFI

         REAL, ALLOCATABLE :: TFIWORK(:)

      END MODULE TFI

      MODULE TIT

         CHARACTER TITLE * 80

      END MODULE TIT

      MODULE WALL

         REAL, ALLOCATABLE ::
     >      UWALL(:,:), VWALL(:,:), XWALL(:,:), YWALL(:,:), ZWALL(:,:)

      END MODULE WALL

      MODULE WNGSRF

         USE DIMENS

         INTEGER
     >      NLG(MXKWING), NW(MXKWING), NWSEC,
     >      KCRANK(MXCRANK), MCRANK(MXCRANK), NCRANK

         LOGICAL
     >      LCRANK(MXCRANK)

         REAL
     >      XWG(MXIWING,MXKWING), YWG(MXIWING,MXKWING),
     >      ZWG(MXIWING,MXKWING), CHORDG(MXKWING), ALG(MXKWING),
     >      ROUNDED(MXKWING), ZCRANK(MXCRANK),
     >      YLESHIFT, YSHELF, ZSHELF, FNWSEC

         INTEGER, ALLOCATABLE ::
     >      IINTR(:), JINTR(:), KINTR(:)

         REAL, ALLOCATABLE ::
     >      TINTR(:), UINTR(:), VINTR(:),
     >      XRG(:,:), YRG(:,:), ZRG(:,:), CHORDM(:),
     >      XW(:,:),  YW(:,:),  ZW(:,:),  SNORM(:), TBASE(:),
     >      XWLE(:),  YWLE(:),  ZWLE(:),  XWTE(:),  YWTE(:)

      END MODULE WNGSRF

      MODULE XYZ

         REAL HAND
         REAL, ALLOCATABLE :: X(:,:,:,:)

      END MODULE XYZ

      MODULE XYZ0

         REAL, ALLOCATABLE :: X0(:,:,:,:), S0(:,:,:,:)
         LOGICAL NEEDX0

      END MODULE XYZ0

C*******************************************************************************
C
      PROGRAM CH_GRID
C
C
C     CH_GRID is a stand-alone form of the C-H grid generation in SYN87-SB,
C     a single-block wing/body design code.  It also provides the geometry
C     perturbation and constraint-checking capabilities of SYN87-SB.
C
C     The grid generator began as the WBGRID package from Lockheed, but it has
C     been almost completely rewritten.  The present scheme employs 8 subblocks,
C     with artificial boundaries forward of the leading edge in the plane of the
C     wing (~horizontal), at the wing tip (~vertical), and at the trailing edge
C     (~vertical/spanwise).  Thomas-Middlecoff control of spacing and Sorenson-
C     type control of orthogonality are provided as options at all boundaries.
C
C     The recommended body surface grid option uses parametric techniques,
C     although a WBGRID-like nonparametric scheme is still supported.  The
C     nose and tail regions are inherently the weakest aspects of this single-
C     block C-H topology because of the stretching required.  See the main
C     GRIDWB module for further grid generation details.
C
C
C     AUTHORS:    James Reuther              David Saunders
C     -------     RIACS                      Sterling Software/Raytheon STX
C                 NASA Ames Research Center, MS 227-6
C                 Moffett Field, CA 94035-1000
C                 (650) 604-1516             (650) 604-1480
C
C                 Also: Stephen Edwards, U.C. Davis, contributed significantly
C                 to the improved TFI and elliptic smoothing utilities.
C
C
C     SPONSORS:   High Speed Aerodynamics Branch, NASA Ames Research Center
C     --------
C
C     COORDINATE SYSTEM:
C     -----------------
C
C     The coordinate system is right-handed as used by the flow solvers of
C     Antony Jameson such as FLO87 and AIRPLANE:
C
C                 I  &  X  Streamwise
C                 J  &  Y  Vertical
C                 K  &  Z  Spanwise
C
C                 I = 1 is at or beyond the lower trailing edge.
C                 J = 1 is at the wing surface.
C                 K = 1 is at the center plane or fuselage surface.
C
C     Downstream, up, and towards the (left) tip are all positive.
C
C
C     APPLICATION-INDEPENDENT NUMERICAL UTILITIES USED:   numutil.f
C     ------------------------------------------------
C
C
C     FILES USED:
C     ----------
C
C        File names are hard-coded as 'ch_grid.EXT' for various EXTs shown.
C        Logical unit numbers appear in the 'LUNS' module.
C
C        LUN    EXT I/O TYPE
C
C        LGEOM  geo  I  ASCII        Geometry inputs as for SYN87-SB
C        IREAD  inp  I  ASCII        Standard inputs except geometry
C        IWRIT  out  O  ASCII        Standard output file
C        LWING  wng  O  ASCII        Perturbed geometry in syn87.geo format
C        LUNBOD body O  Unformatted  Regularized body geometry in PLOT3D form
C        LUNXYZ xyz  O  Unformatted  PLOT3D grid   (/mgrid, but just 1 block)
C
C
C     CONTROL INPUTS:
C     ---------------
C
C        NDV       Number of design variables (perturbing shape functions) >= 0
C        NBOUND    Number of design variable bound inputs overriding the
C                  defaults of [-BIGBND, BIGBND]:
C                    >0 = read this many bounds ahead of the linear constraints;
C                     0 = no such bound inputs except for the two heading lines;
C                    -1 = read bounds from the last two columns of the standard
C                         variable inputs; include the two heading lines below
C                         for ease of switching to the NBOUND > 0 mode
C        NCLIN     Number of linear constraints
C        NCNLN     Number of nonlinear constraints
C        NNFLOOR   Controls monitoring of cabin floor angles and section radii:
C                     0 = no cabin section monitoring;
C                     1 = read cabin floor and radius data following body data;
C                     2 = as for 1 but adjust the data if body camber changes;
C                         the data in syn87.geo has the following format:
C                         # CABIN FLOOR BREAK POINTS (= # FLOOR SEGMENTS + 1)
C                         nfloor
C                         XFLOOR    YFLOOR
C                         xfloor    yfloor
C                         :         :
C                         # CABIN SECTION RADII FOR FLOOR SEGMENT 1
C                         nradii
C                         RBODY     THBODY   (ORIGIN AT FLOOR ON SYMMETRY PLANE)
C                         rbody     thbody
C                         :         :
C                         <Repeat for nfloor-1 cabin sections>
C
C                         Initially, tables of floor angles and cabin radii
C                         are printed at the end of each design iteration.
C        WINGIN    Indicates whether wing sections are planar or not;
C                  the first section may be an exception if NNROOT = 1;
C                     2 = 2-column input (and output) sections;
C                     3 = 3-column input (and output) sections;
C                     4 = 2-column input sections and 3-column output sections
C        WINGOUT   Indicates desired form of saved wing geometry:
C                     0 = sections are untwisted and normalized;
C                     1 =    "      "    twisted and normalized;
C                     2 =    "      "    twisted and denormalized
C        ALPHA     Angle of attack (needed with any floor angle constraints)
C
C
C                  GRID GENERATION CONTROLS:
C
C        WINGONLY     T means no fuselage (wing-on-wall case)
C        PERTURBONLY  T means perturb the geometry & save it, but skip the
C                       grid generation; fuel volume printed if NFUEL > 0
C        READGRID     T means read a grid file (*.xy0); if RESMOOTH = T, it will
C                       be smoothed further; if REDISTRIBUTE = T, its radial
C                       spacing will be changed; cannot apply to design mode
C        PARABODY     T means generate the body surface grid via parametric
C                       bicubic interpolation; NNROOT > 0 assumed;
C                       otherwise, use WBGRID-type linear techniques
C        CHECKWARP    T means after generating an initial grid, perturb the
C                       FIRST variable by the input AITCH(1) then immediately
C                       perturb the geometry, regrid it in "warp" mode, save
C                       the grid for checking, and stop;
C                       VSCALE(1) = 0. can force zero surface perturbation for
C                       detecting spurious grid perturbations
C        RESMOOTH     T means a partial grid (*.xy0 with all boundaries but no
C                       or only partial interior smoothing) is to be read for
C                       further smoothing; see READGRID
C        REDISTRIBUTE T means apply D1ROOTNS, etc., to change the initial
C                       stretching via radial redistribution; see READGRID
C
C        NNROOT    Controls handling of the wing/body intersection:
C                     0 = planar mode - a wing section is lofted at Z = ZPLANAR
C                         *** disabled: use NNROOT = 3 now ***
C                     1 = read the wing/body intersection
C                         *** disabled: use NNROOT = 3 now ***
C                     2 = compute the intersection using sections KWING1:KWING2
C                     3 = as for NNROOT = 2 but check if the wing is below (or
C                         nearly below) the body; if necessary modify sections
C                         IBOD1:IBOD2 using the lower surfaces of wing sections
C                         1:KWING2 so that the intersection calculation may
C                         proceed; see also YSHELF and ZSHELF
C
C        NNTIP     Controls the grid generation at and beyond the wing tip:
C                     0 = squared tip, C-H grid outboard;
C                     1 = squared tip, C-O grid outboard;
C                     2 = rounded tip, C-O grid outboard;    ** INCOMPLETE **
C                     3 = input tip cap geometry between the wing & body data,
C                         with C-O grid outboard
C
C        IL        Number of grid points in the streamwise direction
C        JL        Number of grid points in the radial direction
C        KL        Number of grid points in the spanwise direction
C        KTIP      Spanwise grid index of the (first) tip station
C        JCR       Radial grid index of the body crown/keel line
C        NWING     Number of grid points on wing sections (wrap-around)
C        NBODY     Number of grid points along the body crown/keel line
C
C        KWING1,   Wing sections to use for the wing/body intersection
C         KWING2   calculation if NNROOT >= 2. E.g., 1 & 2. Local cubic splines
C                  connect corresponding chordwise points (= straight lines if
C                  KWING2 = KWING1 + 1). KWING1 need not be 1.
C        IBOD1,    Body sections to use for the NNROOT = 3 option to protect
C         IBOD2    the wing/body intersection calculation
C        IDEGBOD   1 means linear interpolation of body geometry stations to
C                    a common number of points for all stations, and all other
C                    body surface grid interpolations are bilinear, not bicubic;
C                  3 means cubic spline regularization of body sections ("loose"
C                    if NNROOT < 3; monotonic if NNROOT = 3), and bicubic body
C                    surface grid interpolations
C        NBLAYER   Grid spacing for J = 1:NBLAYER varies geometrically off the
C        RBLAYER   wing along subblock edges at I = ILE and I = ITL/ITU;
C                  RBLAYER is the multiplier (e.g., 4 and 1.1 for Euler grids;
C                  20 and 1.05 for Navier-Stokes); Vinokur distributions are
C                  used for J = NBLAYER-1:JL, with D2JL controlling the outer-
C                  most increment
C
C        CRDSPL,   Controls for chordwise surface grid spacing via FOILGRD;
C        CRDSPQ,   try 0.04, 0., 0.3, 0.66 for the linear, quadratic, sine, and
C        CRDSPS,   cosine terms respectively; the sum should be 1.0
C        CRDSPC
C        SPNSPC,   Weights on the cosine and sine distributions combined with
C        SPNSPS    a uniform distribution for the spanwise surface grid spacing:
C                  cosine weighting increases density towards the root and tip;
C                  sine weighting increases density towards the tip
C
C        D1ROOT,   Controls for the initial increments of radial lines off the
C        D2ROOT,   wing root, at the leading edge and trailing edge and the far
C        D3ROOT    downstream boundary respectively, all relative to root chord:
C                     D1ROOT ~ 0.003 for Euler, or 5.E-7 for Navier-Stokes;
C                     D2ROOT ~ 0.004  "    "    "  1.E-6  "    "    "    "
C                     D3ROOT ~ 0.005  "    "    "  0.005  "    "    "    "
C                  If REDISTRIBUTE = T, the D*NS inputs below are imposed
C                  following initial use of these inputs.
C        TEROOT    Fudge factor controlling the slope of the grid line off the
C                  root trailing edge:
C                     TEROOT = 1. gives the traiing edge mean-line slope;
C                     TEROOT > 1. amplifies the mean-line slope (same sign);
C                     TEROOT < 1. reduces the slope (same sign);
C                     TEROOT = 0. zeros the slope;
C                     TEROOT < 0. reverses the slope
C
C        D1TIP,    Controls analogous to D1ROOT, etc., applied at the wing tip
C        D2TIP,    and relative to the tip chord; the root, tip, and Kmax values
C        D3TIP,    are interpolated linearly for other spanwise grid stations
C        TETIP
C
C        D1KMAX,   Controls analogous to D1ROOT, etc., applied at the far
C        D2KMAX,   boundary beyong the wing tip (K = Kmax)
C        D3KMAX,
C        TEKMAX
C
C        D1CROWN   Controls the off-nose increment for an artificial boundary
C                  forward of the nose, relative to the root chord
C        D1NOSE    Controls the axial increments at the nose, relative to the
C                  root chord
C        D1TAIL    Controls the axial increments at the end of the body,
C                  relative to the root chord
C        DWTAIL    Controls the radial increments at the end of the body,
C                  at the wake/water line, relative to the root chord;
C                  should not be much bigger than D2ROOT for N-S cases,
C                  else TFI3D overshoots in the wake between TE and end of body
C        DCTAIL    Controls the radial increments at the end of the body,
C                  at the crown/keel lines, relative to the root chord
C
C        OSPNOSE,  Factors applied to the uniform grid interval for the outer
C        OSPTAIL   C boundary spacings opposite the nose and aft of the tail;
C                  try OSPNOSE = 1.0 and OSPTAIL = 3.0
C        OSPSYM,   For the C-H + C-O case, similar control of the back plane
C        OSPTIP    outer boundary spacing at the symmetry and wing planes;
C                  try OSPSYM = 0.20 and OSPTIP = 2.0
C        D2JL      Increment at outer boundary used with D1CROWN, D1ROOT, etc.,
C                  for edge distributions at I = 1, ITL, ILE, ITU, and IL;
C                  D2JL is scaled by the root chord before use
C
C        SWEEP,    Sweep and dihedral of grid lines beyond the tip; enter
C        DIHED     999. to use angles derived from the last two wing sections
C        BSWEEP    Sweep of outer boundary grid planes; 999. defaults to SWEEP
C        YLESHIFT  Vertical shift applied to wing geometry as it is read, for
C                  fudging extreme cases of high- or low-mounted wings
C        YSHELF    Safety margin (e.g. 1") used with NNROOT = 3 option to help
C                  ensure a findable wing/body intersection; if body section Ys
C                  are replaced by interpolated wing Ys, YSHELF is included
C        ZSHELF    For the NNROOT = 3 option, this forces a minimum width of any
C                  "shelf" built on to body sections to provide room for the
C                  body grid between the keel line & the wing/body intersection;
C                  ZSHELF is applied as a fraction of the local body width;
C                  the input value should be a compromise between the value
C                  associated with the trailing edge (e.g. 0.6) and the higher
C                  value that would increase room for the body grid under the
C                  rest of the wing; the trailing edge region is assisted by
C                  monitoring the upper wing surface to avoid a too-wide shelf
C
C        NCRANK    Number of span stations at which the wing surface mesh should
C                  include grid stations at ZCRANK(1:NCRANK), either to capture
C                  a crank in the leading (or trailing) edge or to confine non-
C                  planar grid sections to the inboard wing panel.  If leading
C                  edge (slat?) deflections are likely, include a crank inboard
C                  of the deflections because the I = ILE sheet will capture the
C                  slope for planar wing grid sections only. (Nonplanar sections
C                  have ill-defined leading edge slopes.)  NCRANK >= 0
C        ZCRANK    Span stations to capture in the grid (see NCRANK); input wing
C                  section Zs are used if they match to within semispan * 0.001,
C                  else a new section must be inserted to allow arc-length-based
C                  spanwise gridding
C
C        MCRANK    Mesh Ks corresponding to crank Zs.  These are normally
C                  calculated by the wing surface gridding, and so the inputs
C                  are just place holders, but the READGRID case needs correct
C                  values because that case bypasses the wing surface gridding.
C
C        LCRANK    T or F flags enabling activation or suppression of MCRANK(*)
C                  K stations as additional boundaries during the REGRID option
C
C        XMIN      Outer boundary center-line coordinates, streamwise
C        XMAX
C        YMIN      Outer boundary coordinates, vertical
C        YMAX
C        ZMIN      Center-line and beyond-tip boundary coordinates, spanwise
C        ZMAX
C        RADIUS    Radius of upstream outer-boundary circular quadrants
C
C                  RADIAL REDISTRIBUTION CONTROLS (used if REDISTRIBUTE = T)
C
C        NBLAYERNS Try 20 and 1.05 for Navier-Stokes grids
C        RBLAYERNS
C
C        D1ROOTNS  Try 0.5E-6 for N-S
C        D2ROOTNS
C        D3ROOTNS  Try the Euler value
C
C        D1TIPNS
C        D2TIPNS
C        D3TIPNS
C
C        D1KMAXNS
C        D2KMAXNS
C        D3KMAXNS
C
C        D1CROWNNS See D1CROWN, etc.
C        DWTAILNS
C        DCTAILNS
C        D2JLNS
C
C                  2D ELLIPTIC SMOOTHING CONTROLS (14 cases annotated in *.inp):
C
C        ITMAX2D   # smoothing iterations
C        OMG2D     Over-relaxation factor for SLOR smoothing iterations
C        CONV2D    # orders of magnitude reduction in max. dX or dY
C        CONVMIN   >= 0.; > 0. suppresses foreground (orthogonality) control
C                  until there is some degree of background-only convergence
C        EXPI12D   Exponential decay factors for foreground terms:
C        EXPI22D   FGMODE2D(n:n) = 'Y': scheme is index-based - try 0.4-0.5;
C        EXPJ12D   FGMODE2D(n:n) = 'A': input is the fractional arc length off
C        EXPJ22D                        the boundary at which the foreground
C                                       coef. should decay from 1.0 to 0.1;
C                  FGMODE2D(n:n) = 'N': orthogonality control is suppressed
C        POWERI    1. means linear combination of opposite edge Phi & Psi for
C        POWERJ    interior background terms: > 1. emphasizes lower edge;
C                  < 1. deemphasizes lower edge; 1. and are 0.5 suggested
C        BGMODE2D  'YYYY' turns on background control at I1, I2, J1, J2 edges;
C                  use 'N' to turn an edge off
C        FGMODE2D  'NNNN' turns off foreground control at the edges likewise;
C                  see EXPI12D, etc., for two ways of controlling orthogonality
C        SPMODE    'YYYY' means spacing off interior edges is interpolated
C                  from corner values; use 'N' to retain starting guess
C                  increments along an entire edge
C
C                  3D ELLIPTIC SMOOTHING CONTROLS (5 cases):
C
C        ITMAX3D   Analogous to 2D controls
C        OMG3D
C        CONV3D
C        EXPI13D   Inputs depend on FGMODE3D as for the 2D case
C        EXPI23D
C        EXPJ13D
C        EXPJ23D
C        EXPK13D
C        EXPK23D
C        BGPHI     'YYYY' turns on background control of interior spacing from
C                  the I faces (Phi effects) at edges J1, J2, K1, K2 resp.
C        BGPSI     'YYYY' turns on background control of interior spacing from
C                  the J faces (Psi effects) at edges I1, I2, K1, K2 resp.
C        BGOMG     'YYYY' turns on background control of interior spacing from
C                  the K faces (Omega effects) at edges I1, I2, J1, J2, resp.
C        FGMODE3D  'NNNNNN' turns off foreground effects at all 6 faces;
C                  'YYYYYY' turns on index-based orthogonality control;
C                  'AAAAAA' turns on arc-based orthogonality control;
C                  any mixture of controls is permitted (e.g., 'NYAANN')
C
C                  MISCELLANEOUS ELLIPTIC SMOOTHING CONTROLS:
C
C        PRINT2D,  .FALSE. suppresses elliptic smoothing iteration print-out
C        PRINT3D
C        FGLIM2D   Foreground term growth limiters;
C        FGLIM3D   1.0 normally suffices
C        URFG2D    Under-relaxation factors for changes in foreground terms;
C        URFG3D    0.1 normally suffices
C
C                  SURFACE SMOOTHING CONTROLS:
C
C        THBUMP    Exponent of sine bump used to recover maximum thickness;
C                  currently inactive
C        IYSMOO    Controls smoothing of Y coordinates of airfoil sections;
C                  may be helpful to mitigate wiggles at the trailing edge:
C                     0 = no smoothing;
C                     1 = smoothing on
C        EPSY      Implicit smoothing coefficient, such as 0.01
C        EXPY      Exponent of smoothing: higher exponents such as 5. will
C                  concentrate the smoothing towards the trailing edge
C
C                  DESIGN VARIABLE (SHAPE FUNCTION) INPUTS:
C
C                  *** NOTE:  The different groups of design variables  ***
C                      should be entered in the order indicated below.
C
C        #         Ordinal number of design variable - not actually used, so
C                  it doesn't really matter if they get out of order
C
C        VTYPE     6-character STREAMWISE design variable type or name;
C                  quotes are now optional ('SIN' and SIN are equivalent):
C
C                  TWT      Twist function                    (planform)
C                  THK      Thickness function                (planform)
C                  XLD      X Leading edge position function  (planform)
C                  YLD      Y Leading edge position function  (planform)
C                  ZLD      Z Leading edge position function  (planform)
C                  CRD      Chord length function             (planform)
C
C                  SIN      Standard (Modified) Sine function (section Y)
C                  SINF     Flipped Sine function             (section Y)
C                  SIN1     Symmetric Sine function           (section Y)
C                  SIN2     Symmetric Flipped Sine function   (section Y)
C                  SIN3     Ensures LE/TE symmetry; use       (section Y)
C                           [XA, XB] = [0, 1] & XC in (0, .5)
C                  SIN4     Ensures LE/TE symmetry; use       (section Y)
C                           [XA, XB] = [0, .5] & XC in (0, 1)
C                  COSL     Half [co]sine, peak at left  (LE) (section Y)
C                  COSR     Half [co]sine, peak at right (TE) (section Y)
C                  LCOS     Inverted form of COSL             (section Y)
C                  RCOS     Inverted form of COSR             (section Y)
C                  EXP      Standard Exponential function     (section Y)
C                  LED      Leading edge droop function       (section Y)
C                  TRL      Trailing edge droop function      (section Y)
C                  WAG      Wagner shape function             (section Y)
C
C                  DSIN     Y" sine function                  (wing Y")
C                  DEXP     Y" exponential function           (wing Y")
C                  DLED     Y" leading edge function          (wing Y")
C                  DTRL     Y" trailing edge function         (wing Y")
C
C                  SINEVC   Variable Center Sine function    (flap/slat/shock Y)
C                  SINFVC   and its variations;              (flap/slat/shock Y)
C                  SIN1VC   used for perturbing flaps/slats  (flap/slat/shock Y)
C                  SIN2VC   or for matching shock angles     (flap/slat/shock Y)
C
C                  FLAP     Trailing edge flap function       (rotates X and Y)
C                  FLAPSH       "     "     "     "           (shears Y only)
C                  SLAT     Leading edge slat function        (rotates X and Y)
C                  SLATSH       "     "     "     "           (shears Y only)
C
C                  BCSIN    Standard sine function            (body camber)
C                  BCCOS    Standard cosine function          (body camber)
C                  BCEXP    Standard exponential function     (body camber)
C                  BCLED    Standard leading edge function    (body camber)
C                  BCTRL    Standard trailing edge function   (body camber)
C
C                  BSSIN    Standard sine function            (body span)
C                  BSCOS    Standard cosine function          (body span)
C                  BSEXP    Standard exponential function     (body span)
C                  BSLED    Standard leading edge function    (body span)
C                  BSTRL    Standard trailing edge function   (body span)
C
C                  BRSIN    Standard sine function            (body radius)
C                  BRCOS    Standard cosine function          (body radius)
C                  BREXP    Standard exponential function     (body radius)
C                  BRLED    Standard leading edge function    (body radius)
C                  BRTRL    Standard trailing edge function   (body radius)
C
C                  BASIN    Standard sine function            (body area)
C                  BACOS    Standard cosine function          (body area)
C                  BAEXP    Standard exponential function     (body area)
C                  BALED    Standard leading edge function    (body area)
C                  BATRL    Standard trailing edge function   (body area)
C                  BASV1    Symmetric volume conserving sine  (body area)
C                  BASV2    Alternative sym. "   "   "   "    (body area)
C
C                  XFLR     X location of cabin floor         (kink point)
C                  YFLR     Y location of cabin floor         (kink point)
C                           (perturbations to initial coordinates;
C                            enter kink point # via KCEN, aka UBUMP)
C
C        UTYPE     4-character name of the SPANWISE or CIRCUMFERENTIAL shape fn.
C                  used to "taper off" the shape function VTYPE (whose multiple
C                  is a design variable) in the transverse direction;
C                  for a wing VTYPE, UTYPE = 'POLY' indicates the original
C                  polynomial (usually linear) spanwise lofting of the
C                  perturbing function, using PBUMP, KCEN, KMIN, & KMAX below;
C                  UTYPE = 'POLY' and PBUMP = 0. are required for flaps/slats;
C                  otherwise, UTYPE may indicate sinusoidal lofting for wing or
C                  body via one of the appropriate VTYPEs such as SIN2 or COSL
C
C        DOUP      Control for application of wing variables to UPPER surface;
C                  DOUP in conjunction with DOLO allows the shape function
C                  to modify camber, thickness, or a single surface:
C                     +1. =  positive influence of design variable;
C                     -1. =  negative influence of design variable;
C                     +/-2.  forces max. thickness conservation   ** INACTIVE **
C                     +/-3.  controls perturbing of Y" rather than Y  ** ???? **
C        DOLO      As for DOUP, applied to LOWER wing surface
C
C                  STREAMWISE (Wing or Body) Shape Function Controls
C                  -------------------------------------------------
C
C        WBUMP     "Width" (exponent) (Sine and Exponential types) or
C                  order of term N    (WAG[ner] function type)
C        XBUMP     "Center" in [0, 1] (Sine and Exponential types) or
C                  center of rotation (Twist type); see also XWMIN/XWMAX
C        XWMIN,    Normalized X (wing or body) to which the shape function is
C        XWMAX     applied - normally 0. and 1., but may be, say, 0.8 and 1.0
C                  to confine it to a flap with hinge at constant X/C; XBUMP
C                  is input in [0, 1] and transformed to [XWMIN, XWMAX]
C
C                  Special FLAP/SLAT Usages
C                  ------------------------
C
C        UTYPE     must be 'POLY' to allow use of lofting input PBUMP = 0.
C        XWMIN,    define the inboard and outboard flap/slat hinge X/Cs resp.
C        XWMAX     (not necessarily the same; hinge is assumed to be straight)
C        PBUMP     (power for spanwise lofting) must be 0. = no tapering off of
C                  the effects of the angle-of-rotation design variable
C
C                  Special Variable Center Function Usages (SIN*VC)
C                  ------------------------------------------------
C
C                  (a) for flaps or slats:
C
C        XBUMP     X/C in (0., 1.) is the shape function "center" referred to
C                  the local chord of the flap or slat at each section
C        XWMIN,    define the inboard and outboard flap/slat hinge X/Cs as for
C        XWMAX     the FLAP/SLAT functions;
C                  XWMAX < 0.5 indicates a slat; XWMAX >= 0.5 indicates a flap
C
C                  (b) for matching shock angles on wing surface:
C
C        XBUMP     must be 0. to distinguish this from the flap/slat case;
C        XWMIN,    define the inboard and outboard wing X/C defining the line
C        XWMAX     parallel to the shock on which to center the shape function
C
C                  SPANWISE (Wing) or CIRCUMFERENTIAL (Body) Shape Fn. Controls
C                  ------------------------------------------------------------
C
C                  There are TWO options for spanwise lofting of perturbations.
C                  CH_GRID uses the PRESENCE or ABSENCE of DECIMAL POINTS to
C                  distinguish the cases.  Input groups KCEN, KMIN, KMAX and
C                  UBUMP, UBMIN, UBMAX are ALTERNATIVES to each other.
C
C                  [The "UB" terminology refers to parametric body coordinates,
C                  with 0. at the keel; the same REAL inputs may also be used
C                  for spanwise lofting of the wing variables, although the
C                  original index-based spanwise lofting is still supported.]
C
C              (1) The original style based on wing defining section K indices:
C                  KCEN, KMIN, & KMAX are input as INTEGERs (no decimals pts.).
C
C        PBUMP     Power applied to the lofting shape function between KCEN and
C                  KMIN, and between KCEN and KMAX.
C                  If UTYPE = 'POLY', the lofting is polynomial:
C                     PBUMP = 1. means the shape function goes to 0. linearly
C                                at stations KMIN & KMAX from its peak at KCEN;
C                                if KMIN = KCEN, the peak applies at KMIN;
C                                if KMAX = KCEN, the peak applies at KMAX.
C                     PBUMP = 0. means the perturbation is constant for K = KMIN
C                                through KMAX inclusive, and 0. elsewhere (as
C                                required for FLAP/SLAT functions).
C                  If UTYPE = 'SIN', etc., the lofting is Hicks-Henne-type.
C        KCEN      Center span station of influence for wing design variable;
C                  also used for cabin floor kink point # (XFLR/YFLR variables)
C        KMIN      First (possible) geometry span station influenced by the wing
C                  design variable
C        KMAX      Last (possible) station of influence for wing design variable
C
C              (2) Lofting in the transverse direction is controlled by spanwise
C                  coords. (wing) or by normalized circumferential body coords.
C                  The 3 inputs after PBUMP must be REALs with DECIMAL POINTS.
C
C                  [Use of normalized spanwise quantities for the wing, treated
C                  as a parametric surface in a fashion similar to the body
C                  treatment, was planned as a third option but this has been
C                  deferred as probably unnecessary, and is not consistent with
C                  the automated means of distinguishing between INTEGER and
C                  REAL inputs for options (1) and (2), since REALs in the range
C                  [0, 1] could be either normalized or not.]
C
C        PBUMP     Power applied to the UTYPE shape function on [UBMIN, UBMAX]
C        UBUMP     Center of the UTYPE shape function
C        UBMIN,    Lower and upper limits to which UTYPE is confined:
C        UBMAX     for body VTYPEs, 0. is at the keel and 1. is at the crown;
C                  for wing VTYPEs, use unnormalized spanwise (Z) coordinates
C
C                  Design variable inputs most related to the optimizer
C                  ----------------------------------------------------
C
C        V         Design variable value as seen by the optimizer (see VSCALE)
C        VSCALE    Scale factor for the design variable to keep it ~O(1);
C                  V * VSCALE is typically the shape function multiplier;
C                  retained for compatibility with SYN87-SB input scheme
C        AITCH     Step size used for finite differencing of geometry changes;
C                  retained for compatibility with SYN87-SB input scheme
C
C                  See NBOUND description for these optional inputs:
C
C        BL        Lower bound on the design variable as seen by the optimizer;
C                  BL = -999. means the variable has no lower bound, else
C                  BL scaling should match that of the variable - e.g., for a
C                  'TWT' variable with VSCALE = 0.1, BL = -10. means
C                  twist >= -10.*0.1 = -1., which is consistent with AITCH usage
C        BU        Upper bound on the design variable as seen by the optimizer;
C                  BU = 999. means the variable has no upper bound, else
C                  BU scaling should match that of the variable
C
C                  VARIABLE BOUND INPUTS (Optional - see NBOUND above):
C
C        #         Ordinal number of variable bound input - not actually used,
C                  so it doesn't really matter if they get out of order
C        IVAR      Design variable number for which [-BIGBND, BIGBND] bounds
C                  are being overridden with this NBOUND option
C        BL        Lower and upper bounds to be applied to variable IVAR;
C        BU        VSCALE applies to these inputs as described for BL, BU above
C
C                  LINEAR CONSTRAINT INPUTS:
C
C                  Note:    Some quantities such as wing section T/C at a
C                           specified X/C can be constrained efficiently by
C                           exploiting the fact that CHANGES in the quantity
C                           will be linear combinations of perturbing variables;
C                           the constraints are actually imposed on any
C                           PERTURBATIONS in the specified quantity, but the
C                           inputs refer to the quantity itself.
C
C        #         Ordinal number of linear constraint - not actually used
C        LCTYPE    5-character linear constraint type or name:
C
C                  THK      Wing section thickness constraint: T/C bounds given
C                  THK0     Wing section thickness constraint: T/C (min) read is
C                           overwritten by initial T/C found, which is tabulated
C                           for possible reuse; THK0 requires care - it may
C                           impose T/C larger than is really meant;
C                           ignored if ISTART > 0
C                  YUPR,    Linear constraint on upper or lower Y/C
C                  YLWR
C                  YUPR0,   Variants of YUPR/YLWR constraints for determining
C                  YLWR0    initial values of Y/C at the specified Y/C and K;
C                           values are tabulated and used to overwrite both
C                           bound inputs; ignored if ISTART > 0
C                  CAB,     Linear forms of cabin radius constraints, viable if
C                  CABI0    the body camber line is held fixed (else see the
C                           nonlinear constraint description); the CABI0 variant
C                           allows the initial radius to serve as a lower bound;
C                           a cabin floor definition and cabin radii constraints
C                           must be appended to the fuselage geometry in the
C                           syn87.geo file, and NNFLOOR > 0; enter cabin class
C                           and polar angle number as a compound integer:
C
C                           KLCON() = 100 * floor segment + radial angle number
C
C        KLCON     Wing geometry section at which the linear constraint applies
C        XLCON     The X/C value          "    "    "    "    "    "    "    "
C        BL        The lower bound being imposed on T/C or Y/C; -999. = no bound
C        BU        The upper bound being imposed on T/C or Y/C;  999. = no bound
C        [Z        Span station may be included here as a visual aid - not used]
C
C                  NONLINEAR CONSTRAINT INPUTS (see NCNLN above):
C
C        #         Ordinal number of nonlinear constraint - not used
C        NLCTYPE   6-character nonlinear constraint type or name:
C
C                  FUEL     Fuel volume constraint; requires NFUEL > 0 for the
C                           fuel volume definition; should not be used if there
C                           are no appropriate design variables to vary volume
C                  CABIN    Cabin radius-type nonlinear constraint; requires the
C                           cabin floor definition and cabin radii constraints
C                           appended to the fuselage geometry in syn87.geo, and
C                           NNFLOOR > 0
C                  FLOOR    Cabin floor angle constraint
C                  DFLOOR   Constraint on change in cabin floor angle where
C                           two floor segments meet
C                  TCMAX    Wing section max. T/C constraint (at any X/C);
C                           use of 'THK' linear constraints at specified X/Cs
C                           may be more efficient
C                  TEANG    Trailing edge angle constraint; use of 'THK' linear
C                           constraints at specified X/Cs may be more efficient
C                  CLEAR    Clearance between cabin floor and upper wing surface
C                           at specified X and specified wing station (INLCON)
C
C        BL        Lower bound on nonlinear constraint value; -999. = -BIGBND
C        BU        Upper bound on nonlinear constraint value;  999. =  BIGBND;
C                  N.B.:  these bounds should be in the same units as the
C                  quantities being bounded (e.g., fuel volume); SNLCON will
C                  be applied by SYN87 to the given BL and BU; this is in
C                  contrast to the bounds on the variables themselves, which
C                  (as for their finite differencing intervals) refer to the
C                  scaled variables
C        XNLCON    Coordinate needed (or not) by the nonlinear constraint;
C                  CABIN type needs a body X station
C        INLCON    First index needed (or not) by the nonlinear constraint;
C                  CABIN type uses INLCON = floor segment # (~cabin class);
C                  FLOOR   "    "    "    "    "    "    "    "    "    "
C                  DFLOOR  "    "    "    " floor kink point #
C                  TCMAX   "    "    "    " wing geometry section #
C                  TEANG   "    "    "    "    "    "    "    "
C                  CLEAR   "    "    "    "    "    "    "    "
C        JNLCON    Second index needed (or not) by the nonlinear constraint;
C                  CABIN type uses JNLCON = radial angle associated with INLCON
C        SNLCON    Scale factor used to provide the optimizer with nonlinear
C                  constraint values of order 1.
C
C                  FUEL VOLUME CALCULATION INPUTS:
C
C        NFUEL     Number of fuel sub-volume definitions to follow.
C                  0 <= NFUEL <= MXKWING;  NFUEL must be present; suppress
C                  the following inputs if NFUEL = 0
C        KFUEL1    Indices of the inboard wing geometry sections bounding the
C                  sub-volumes of fuel
C        ZFUEL1    Corresponding Zs (to aid in preparing inputs - not used)
C        XAFUEL1,  Values of X/C bounding fore & aft fuel sub-volume at KFUEL1
C        XBFUEL1
C        KFUEL2    Indices and Z coordinates of the outboard sections bounding
C        ZFUEL2    the sub-volumes of fuel
C        XAFUEL2,  Values of X/C bounding fore & aft fuel sub-volume at KFUEL2
C        XBFUEL2
C
C
C     GEOMETRY INPUTS:
C     ---------------
C
C
C     *** NOTE:    Wing sections inside the body are permissible, and may be
C                  desirable for nonlinear fuselage midsections, given that
C                  all input sections (except the first if NNROOT = 1) are
C                  presently 2-dimensional.
C
C                  SAMPLE WING GEOMETRY:
C
C     Bizjet wing <title>
C            ZSPAN      SREF      CREF     XREF    YREF    ZREF
C           263.50  15272.93     62.50     268.     56.     0.
C              FNC    ALPHAW
C         6.000000  0.000000
C                Z       XLE       YLE    CHORD   THICK   TWIST   NEWSEC  ROUND
C               0.  237.4450   55.9363 123.9940   1.000   0.000   1.0     1.0
C             YSYM        NU        NL
C             0.00    100.00    100.00
C    UPPER SURFACE
C         0.000000  0.000000  [0.000000]
C         0.000252  0.004334  [ ...... ]
C         0.001007  0.007771  [ ...... ]
C         0.002264  0.011163  [ ...... ]
C          :   :   :   :   :   :   :   :
C
C                  Description:
C
C        TITLE     Case title, used on the pressure plots
C        ZSPAN     Reference (semi-)span
C        SREF      Reference (semi-)area
C        CREF      Reference chord
C        XREF,     Moment center (e.g., center of gravity)
C        YREF,
C        ZREF
C
C        FNC       Number of input wing sections
C        ALPHAW    Wing angle of attack (in degrees, added to local twists)
C
C        ZLE, XLE, YLE, CHORD, THICK, TWIST, NEWSEC, ROUND
C                  are mostly self-explanatory:
C                  THICK is a scale factor applied to the Y coordinates only;
C                  NEWSEC = 0. means reuse the previous section;
C                  ROUND = 0. or 1. if the leading edge is sharp or blunt resp.;
C                  if WINGIN  = 3, the third (Z) column is expected on input;
C                  if WINGOUT = 1, CHORD refers to the chord BEFORE twisting,
C                                  and TWIST is set to 0.;
C                  if WINGOUT = 2, CHORD = 1., and ZLE = XLE = YLE = TWIST = 0.
C
C        YSYM      Must be 0.; symmetric sections are not handled specially
C        FNU       Number of input points on the upper & lower surfaces
C        FNL       If NNROOT = 2, sections KWING1:KWING2 must have the same
C                  FNU and FNL; if NNROOT = 3, this applies to 1:KWING2
C
C        XU, YU[, ZU]
C                  Wing section coordinates; for the first section, ZU = 0.
C                  is assumed unless NNROOT = 1 is being used; if NNROOT is
C                  not 1, any third column for the first section is ignored
C
C                  SAMPLE FUSELAGE GEOMETRY (APPENDED TO WING GEOMETRY):
C                  [See also NNFLOOR for optional further inputs.]
C
C              FNF      FLIP
C           174.00      0.00
C               FN          XF        FNEW
C         1.000000   34.000000    1.000000
C               ZF          YF
C         0.000000   79.500000
C               FN          XF        FNEW
C        70.000000   34.479900    1.000000
C               ZF          YF
C         0.000000   77.971000
C         0.081900   77.972100
C         0.163700   77.975500
C         0.245400   77.979900
C          :   :   :   :   :
C
C                  Description:
C
C        FNF       Number of input fuselage sections
C        FLIP      Allows for flipping Y and Z:
C                     0. = Z (spanwise) and Y (vertical);
C                     1. = Y followed by Z; also, FN & XF are flipped
C
C        FN        Number of points defining the section (1. at the nose)
C        XF        X coordinate of the section
C        FNEW      Use 0. to repeat a section, else 1.
C
C        ZF        Section coordinates (possibly in the other order - see FLIP)
C        YF
C
C                  NOTE:  Body sections may be entered from keel to crown or
C                         crown to keel; the program ensures they are from
C                         keel to crown by reversing the order if necessary.
C
C                  [See also NNFLOOR for optional further inputs.]
C
C
C     HISTORY:
C     -------
C
C     The history of the grid generation package has been moved to GRIDWB.
C     The following covers other aspects such as design variables & constraints.
C
C------------------------------ OPT67/SYN87 Origins ----------------------------
C
C     ??/92-  J.Reuther  Original OPT67, using WBGRID wing/body grid generation
C     ??/93              from Lockheed.
C     08/94-  D.Saunders OPT67 overhauled to improve modularity and workspace
C     03/95              use.  The grid generation was largely rewritten.
C     04/95        "     SYN87 updated with OPT67 improvements.
C     07/95-       "     Geometry & grid generation now use the flow solver
C     12/95              convention for Y and Z throughout.  Geometry inputs
C                        are now separated from other standard inputs.
C     02/96        "     Added flap, slat, and Wagner shape functions.
C     03/96     JJR/DAS  Initial incorporation of constrained optimization
C                        via NPOPT, with QNMDIF2 option retained.
C     04/96        "     Stage 2 of 3-stage WARP3D scheme is improved;
C                        Hicks-Henne functions apply to [XWMIN, XWMAX].
C     05/96       DAS    NNFLOOR introduced along with initial monitoring
C                        of cabin radii;
C                        KWING1, KWING2 inputs replace IDEGWING; added
C                        'THK0 ', 'YUPR ', 'YUPR0', ... constraints;
C                        added fuel volume option (NFUEL); added max. %T/C
C                        and TE angle to wing section summary table.
C     06/96        "     Arranged for PERTURB to be called first from one of
C                        several places (initialization issue);
C                        provided for bounding the variables (NBOUND);
C                        CONFUN allows NPOPT to constrain cabin radii,
C                        fuel volume, floor angles, T/Cmax, or trailing
C                        edge angles; NLCON called from CONFUN or PERTURB
C                        either evaluates nonlinear constraints (scaled)
C                        or prints (some of) them (unscaled) respectively;
C                        perturbations of cabin floor kink point (X,Y)s
C                        may be used as design variables (to help meet
C                        constraints on floor angles or cabin radii);
C                        linear constraints now have upper bounds and are
C                        summarized in T/C or Y/C form at start and end.
C     07/03/96     "     Print T/C, not %T/C (consistent with linear and
C                        nonlinear constraints);
C                        cabin radius printout shows Ys & Zs too.
C     08/07/96     "     Constrained floor/wing clearance ('CLEAR').
C     08/23/96     "     SFEVAL introduced in SETLCON for precise linear
C                        constraint set-up, and PERTURB (wing Y, Y" only).
C     09/30/96     "     Body design variables have been extended (radial and
C                        vertical in addition to spanwise and area); variable
C                        input scheme is correspondingly affected (for the wing
C                        as well); SFEVAL vectorizes the body variables;
C                        0/1 flag eliminated from body section inputs.
C     10/07/96     "     Introduced NNROOT = 3 option, with IBOD1 & -2 to
C                        protect the wing/body intersection calculation if
C                        the wing drops below the body.
C     10/09/96     "     YSHELF & ZSHELF are additional inputs (NNROOT=3).
C     10/16/96     "     Refined WBFUDGE; ZSHELF is relative to Zmax now.
C     11/02/96     "     BILINT used by 2D search now has LUSOLVE in-lined
C                        (but cf77 in-lining of numutil.f routines turns
C                        out to be counter-productive).
C     11/13/96  JJR/DAS  Body design variables force 6-char. design variable
C                        names; unit body perturbations introduced.
C     11/14/96     "     Body fudging in-place meant body sections in *.wng
C                        were no longer oval ... Y/ZFUDGE(*,*) introduced.
C     11/15/96     "     XBMIN, XBMAX eliminated - body variables now use
C                        XWMIN(N), XWMAX(N).
C     11/16/96     "     Cabin radius constraints can be LINEAR if the body
C                        camber line is fixed.
C     12/06/96  DAS/JJR  Shape functions COSL, COSR, LCOS, RCOS added to
C                        SFEVAL for body crown/keel (& maybe wing LE/TE).
C     01/17/97    DAS    Added T as well as T/C to constraint table.
C     01/31/97     "     RDVLINE provides for original index-based lofting
C                        of wing perturbations or ZCEN/ZMIN/ZMAX option,
C                        in anticipation of use of Hicks-Henne functions.
C     02/26/97     "     Completed application of Hicks-Henne functions to
C                        spanwise lofting.  LOFT combines revised RIPPER and
C                        SFEVAL for use in PERTURB and SETLCON.
C     03/10/97     "     No-floor-data case & wing-only case revealed glitches.
C
C------------------------------ CH_GRID Adaptation -----------------------------
C
C     03/11/97    DAS    Initial extraction of CH_GRID from SYN87-SB.
C     03/27/97  DAS/JJR  WINGIN, WINGOUT controls added as part of a more
C                        careful treatment of local twist, chord, etc.
C     Aug. 1997   DAS  o FLAP/FLAPSH & SLAT/SLATSH (rotate & shear pairs)
C                        provide for arbitrary hinge lines, replacing the
C                        original FLAP/SLAT functions (which were shears
C                        with constant X/C for the hinge assumed).
C                        The originals are still achievable with LED/TRL
C                        and WBUMP=1. (linear) & PBUMP=0. (no taper-off).
C                      o SIN*VC design variables (Variable Center family)
C                        provide for perturbing slats or flaps or for
C                        centering bumps along an arbitrary line such as a
C                        line parallel to a shock off a nacelle/diverter.
C     Apr. 1998    "     Some improvements enabled by Fortran 90, and
C                        refinements consistent with AEROSURF.  E.g.,
C                        generalized WBINTR is now used by WFINTR.
C     May  1998    "   o Made CH_GRID executable independent of the grid
C                        size and the number of variables and constraints.
C                      o Fuel volume can be computed in PERTURBONLY mode.
C     June 1998    "     Symmetric wing+body revealed NF(*)/NJBODY mix-up.
C     July 1998    "     A radial redistribution option to improve N-S grids
C                        leads to reinstating the design code's "warping"
C                        option and a top-level reorganization so that the
C                        GRIDWB module is common to CH_GRID and SYN87-SB.
C     Sep. 1998    "     GRIDWB needed READGRID and NNACEL usage to stay
C                        identical to what is needed in the design code.
C     Nov. 1998    "     Minor change in WSURF to match a SYN87-SB change
C                        involving the tip as a crank for ghost nacelle cases.
C                        Also: WSURF redistributes the tip mean line now.
C     Dec. 1998    "     WBFUDGE now captures the shelf edge only if the
C                        previous J point was fudged.  Z0PLANE was copying
C                        the water-line from 1:ITU instead of 1:IBU.
C                        MXJBODY and NJBODY raised from 125 & 101 to 151.
C     Jan. 1999    "     Eliminated UFINIT/VFINIT/XFINIT(*,*) in PERTURB.
C     Aug. 1999    "     Revised OFFTEDGE handles rounded trailing edges;
C                        SIN3 and SIN4 functions ensure fore/aft symmetry.
C
C*******************************************************************************

C     Global variables:

      USE DIMENS
      USE CONSTS
      USE FUEL
      USE FUSEL
      USE FUSEL2
      USE INDEX
      USE LIM
      USE LUNS
      USE LOGIC
      USE OPT
      USE OPTIONS
      USE SPACING
      USE SIZES
      USE TIT
      USE WNGSRF
      USE XYZ
      USE XYZ0

      IMPLICIT NONE

C     Local constants:

      REAL, PARAMETER ::
     >   BIGBND = 1.E+10, HALF = 0.5, ONE = 1., ZERO = 0.
      CHARACTER, PARAMETER ::
     >   PROGNAME * 7 = 'CH_GRID'

C     Local variables:

      INTEGER
     >   I, IOS, J, K, K1, K2, M, N, NCALL, NFUELIN, NL, NU

      REAL
     >   BOUNDLO, BOUNDUP, COSTWT, EFFTWT, R1, R2,
     >   SINTWT, TEANGL, XA1, XA2, XB1, XB2,
     >   XTEDGE, XTRAIL, YTEDGE, YTRAIL, ZK

      REAL(KIND=4) CPU1, CPU2, TIME1, TIME2

      LOGICAL
     >   FAIL, NEWGEOM

      CHARACTER
     >   BUFFER*132, LTYP*5

C     Procedures:

      EXTERNAL
     >   CHKLCON, CPUTIME, DYNAMIC, FOILGRD, FUELVOL, FUSL, GEOMWR,
     >   GETANGLE, GRIDINP, GRID_ELLIP3D, GRIDWB, GRIDIO, NLCON,
     >   PERTURB, RADIAL_INIT, RADIAL_SET, RDVARS, REGRID_FACES,
     >   REGRID_VOLS, SECOND, SETLCON, WING, WREGUL

C-----------------------------------------------------------------------

C     Execution:

      CALL SECOND (TIME1)

      I = MXJBODY ! Kludge to avoid "no path" diagnostic.
      IF (I < NJBODY) THEN
         WRITE (IWRIT,1040) 'MXJBODY < NJBODY'
         GO TO 999
      END IF

      I = MXIMAX
      IF (I < MAX (MXIWING, MXIBODY)) THEN
         WRITE (IWRIT,1040) 'MXIMAX < MAX (MXIWING, MXIBODY)'
         GO TO 999
      END IF

      I = MXKMAX
      IF (I < MAX (MXKWING, MXJBODY)) THEN
         WRITE (IWRIT,1040) 'MXKMAX < MAX (MXKWING, MXJBODY)'
         GO TO 999
      END IF

      OPEN (UNIT=IREAD, FILE='ch_grid.inp', STATUS='OLD')
      OPEN (UNIT=LGEOM, FILE='ch_grid.geo', STATUS='OLD')
      OPEN (UNIT=IWRIT, FILE='ch_grid.out', STATUS='NEW')

      WRITE (IWRIT,1)
    1 FORMAT (' CH_GRID Version 11-AUG-1999', /,
     >        ' James Reuther, David Saunders, Stephen Edwards', /,
     >        ' High Speed Aerodynamics Branch', /,
     >        ' NASA Ames Research Center, MS 227-6', /,
     >        ' Moffett Field, CA  94035-1000', /)

C     Handy constants:

      EPSMCH  = EPSILON (EPSMCH)
      PI      = 4.* ATAN (ONE)
      RPI     = ONE / PI
      DEG2RAD = PI / 180.
      RAD2DEG = 180./ PI
      LOGPT5  = LOG (HALF)

C     Echo the inputs to the output file:

      WRITE (IWRIT,'(/,A,/)') ' Control inputs:'

      DO
         READ (IREAD,1001,IOSTAT=IOS) BUFFER
         IF (IOS < 0) EXIT

         I = LEN_TRIM (BUFFER)
         WRITE (IWRIT,1010) BUFFER(1:I)
      END DO

      REWIND (IREAD)

      READ (IREAD,*)
      READ (IREAD,*)
      READ (IREAD,*) NDV, NBOUND, NCLIN, NCNLN, NNFLOOR,
     >               WINGIN, WINGOUT, ALPHAI

      NNACEL = 0      ! WSURF needs this
      ALPHAF = ALPHAI ! Alpha is needed for floor slopes
      NBNDS  = NDV + NCLIN + NCNLN

      IF (NBOUND > NDV) THEN
         WRITE (IWRIT,1040) 'More variable bounds than variables.'
         GO TO 999
      END IF

      READ (IREAD,*)
      READ (IREAD,*) WINGONLY, PERTURBONLY, READGRID, PARABODY,
     >               CHECKWARP, RESMOOTH, REDISTRIBUTE
      GRIDONLY   = .TRUE.  ! Design code option
      OPTIMIZING = .FALSE. ! "  "

      READ (IREAD,*)
      READ (IREAD,*) NNROOT, NNTIP

      NNGRID = 2 ! Only case relevant to CH_GRID (see CHECKWARP option below)
                 ! The following is from the SYN87-SB input description:

C     NNGRID    Controls the grid generation strategy (see REDISTRIBUTE too):
C                  0 = full grid from scratch
C                      (for all grids in the design case)
C                  1 = full grid from scratch
C                      (for initial grid in design case - perturb thereafter)
C                  2 = full grid at the end of each design step (NCALL <= -2)
C                  3 = full grid in the line search too (NCALL <= 0)

      READ (IREAD,*)
      READ (IREAD,*) IL, JL, KL, KTIP, JCR, NWING, NBODY

C     Allocate dynamic storage for the grid and the variables & constraints:

      CALL DYNAMIC ()

C     Read remaining grid generation inputs:

      CALL GRIDINP ()


C     The normalized chordwise distribution is needed for fuel volume
C     calculations as well as for the grid generation, so set it up once:

      CALL FOILGRD (NHALF, ZERO, ONE, CRDSPL, CRDSPQ, CRDSPS, CRDSPC,
     >              SNORM)


      IF (.NOT. PERTURBONLY) THEN

         HAND = ONE ! Right-handed xyz convention for cell volume checker

         OPEN (UNIT=LUNXYZ, FILE='ch_grid.xyz', STATUS='NEW',
     >         FORM='UNFORMATTED')
      END IF


C     Option to read a partially-prepared grid for possible further smoothing
C     and/or radial redistribution.  Skip all the geometry processing.

      IF (READGRID .AND. .NOT. PERTURBONLY) THEN

         OPEN (UNIT=LUNXY0, FILE='ch_grid.xy0', STATUS='OLD',
     >         FORM='UNFORMATTED')

         CALL GRIDIO (0, IL, JL, KL, X, LUNXY0, IWRIT)

         CLOSE (LUNXY0)

         NCALL = -3 ! First call flag

         IF (RESMOOTH .OR. REDISTRIBUTE) ALLOCATE (S0(IL,JL,KL,3))

         IF (RESMOOTH) THEN

            CALL SECOND (CPU1)

            CALL GRID_ELLIP3D (NCALL)

            CALL SECOND (CPU2)

            CALL CPUTIME (CPU1, CPU2, PROGNAME,
     >                    'to smooth the grid', IWRIT)
         END IF

         IF (REDISTRIBUTE) THEN

            CALL SECOND (CPU1)

            CALL RADIAL_INIT (0) ! 0 means get CHORDM & ZWLE from X(*,*,*,*)

            CALL RADIAL_SET (.FALSE.) ! EULER = F means activate N-S controls

            CALL REGRID_EDGES (NCALL)

            CALL REGRID_FACES (NCALL)

            CALL REGRID_VOLS (NCALL)

            CALL SECOND (CPU2)

            CALL CPUTIME (CPU1, CPU2, PROGNAME,
     >                    'to redistribute the grid', IWRIT)
         END IF

         IF (RESMOOTH .OR. REDISTRIBUTE) THEN

            DEALLOCATE (S0)

C           Save the grid (LUNERR = 1 is unused).

            CALL GRIDIO (1, IL, JL, KL, X, LUNXYZ, 1)

         END IF

C        Check the cell volumes:

         CALL CHECK_VOLS (PROGNAME, IL, JL, KL, HAND, X, IWRIT)

         GO TO 999

      END IF


C     Read various surface smoothing controls:

      READ (IREAD,*)
      READ (IREAD,*)
      READ (IREAD,*) THBUMP, IYSMOO, EPSY, EXPY

      NEWGEOM = IYSMOO > 0 .OR. NDV > 0


C     Read the design variable details:

      IF (NDV > 0) THEN

         BL(1:NDV) = -BIGBND ! Default the variable bounds
         BU(1:NDV) =  BIGBND

         CALL RDVARS (IREAD)

      END IF


C     Read the wing geometry from unit LGEOM:

      CALL WING ()


C     Read the fuselage geometry if any, from LGEOM:

      NFLOOR = 0

      IF (.NOT. WINGONLY) THEN

         CALL FUSL ()

C        Regularized body geometry is saved in PLOT3D form by PERTURB.

         OPEN (LUNBOD, FILE='ch_grid.body', STATUS='NEW',
     >         FORM='UNFORMATTED')

C        Option to read cabin floor and radius data:
C        -------------------------------------------

         IF (NNFLOOR > 0) THEN
            READ (LGEOM,*)
            READ (LGEOM,*) NFLOOR
            IF (NFLOOR > MXFLOOR) THEN
               WRITE (IWRIT,1040) 'Too many floor segment defining pts.'
               GO TO 999
            END IF

            READ (LGEOM,*)
            DO I = 1, NFLOOR
               READ (LGEOM,*) XFLOORI(I), YFLOORI(I)
            END DO

            DO I = 1, NFLOOR - 1
               READ (LGEOM,*)
               READ (LGEOM,*) NRADII(I)

               IF (NRADII(I) > MXRADII) THEN
                  WRITE (IWRIT,1040) 'Too many cabin section radial pts'
                  GO TO 999
               END IF

               READ (LGEOM,*)
               READ (LGEOM,*) (RBODY(J,I), THBODY(J,I), J =1, NRADII(I))
            END DO

C           Echo the data to standard output:

            WRITE (IWRIT,1001)
            WRITE (IWRIT,1001)
     >         '# CABIN FLOOR BREAK POINTS (= 1 + # FLOOR SEGMENTS)'
            WRITE (IWRIT,'(I3)') NFLOOR

            WRITE (IWRIT,1001) '   XFLOOR   YFLOOR'
            WRITE (IWRIT,'(2F9.2)')
     >         (XFLOORI(I), YFLOORI(I), I = 1, NFLOOR)

            DO I = 1, NFLOOR - 1
               WRITE (IWRIT,'(A,I2)')
     >            '# CABIN RADII FOR FLOOR SEGMENT', I
               WRITE (IWRIT,'(I3)') NRADII(I)

               WRITE (IWRIT,1001) '    RBODY   THBODY'
               WRITE (IWRIT,'(2F9.4)')
     >            (RBODY(J,I), THBODY(J,I), J = 1, NRADII(I))
            END DO

         END IF

         WRITE (IWRIT,'(A)')
         WRITE (IWRIT,'(A,2I6)') ' End-of-body   grid indices:', IBL,IBU

      ELSE ! Wing only
         WRITE (IWRIT,'(A)')
      END IF

      WRITE (IWRIT,'(A,2I6)')    ' Trailing edge grid indices:', ITL,ITU

      CLOSE (LGEOM)


      READ (IREAD,*) ! Variable bounds (optional)
      READ (IREAD,*) ! --------------------------

      DO M = 1, NBOUND ! NBOUND may be 0 or -1
         READ (IREAD,*,IOSTAT=IOS) J, I, BL(I), BU(I)
         IF (IOS /= 0) THEN
            WRITE (IWRIT,1040) 'Error reading variable bounds.'
            GO TO 999
         END IF
      END DO

      DO I = 1, NDV
         IF (BL(I) == -999.) BL(I) = -BIGBND
         IF (BU(I) ==  999.) BU(I) =  BIGBND
      END DO

C     Otherwise, BL & BU scaling should match that of the variables.


      READ (IREAD,*) ! Linear constraints
      READ (IREAD,*) ! ------------------

C     Separate reads allows trailing text such as unused Zs:

      DO M = 1, NCLIN
         READ (IREAD,*,IOSTAT=IOS) J, LCTYPE(M), KLCON(M), XLCON(M),
     >                             BOUNDLO, BOUNDUP
         IF (IOS /= 0) THEN
            WRITE (IWRIT,1040) 'Error reading linear constraints.'
            GO TO 999
         END IF

         IF (BOUNDLO == -999.) BOUNDLO = -BIGBND
         BLLIN(M) = BOUNDLO

         IF (BOUNDUP ==  999.) BOUNDUP =  BIGBND
         BULIN(M) = BOUNDUP
      END DO


      READ (IREAD,*) ! Nonlinear constraints
      READ (IREAD,*) ! ---------------------

      K = NDV + NCLIN ! Offset for the bounds
      DO M = 1, NCNLN

         READ (IREAD,*,IOSTAT=IOS) J, NLCTYPE(M), BOUNDLO, BOUNDUP,
     >      XNLCON(M), INLCON(M), JNLCON(M), SNLCON(M)

         IF (IOS /= 0) THEN
            WRITE (IWRIT,1040) 'Error reading nonlinear constraints.'
            GO TO 999
         END IF

C        Scaling of the nonlinear constraints differs from scaling of
C        the variables and linear constraints:

         IF (BOUNDLO == -999.) THEN
            BOUNDLO = -BIGBND
         ELSE
            BOUNDLO = BOUNDLO * SNLCON(M)
         END IF
         BL(M+K) = BOUNDLO

         IF (BOUNDUP ==  999.) THEN
            BOUNDUP =  BIGBND
         ELSE
            BOUNDUP = BOUNDUP * SNLCON(M)
         END IF
         BU(M+K) = BOUNDUP
      END DO


      READ (IREAD,*) ! Fuel volume set-up; Zs are there to help prepare inputs
      READ (IREAD,*) ! ------------------

      NFUEL = 0
      READ (IREAD,*) NFUELIN
      READ (IREAD,*)

      DO M = 1, NFUELIN ! Careful: NFUEL >= NFUELIN

         READ (IREAD,*,IOSTAT=IOS) K1, ZK, XA1, XB1, K2, ZK, XA2, XB2

         IF (IOS /= 0) THEN
            WRITE (IWRIT,1040) 'Error reading fuel volume data.'
            GO TO 999
         END IF

         IF (K1 + 1 == K2) THEN
            NFUEL = NFUEL + 1
            KFUEL1(NFUEL)  = K1
            XAFUEL1(NFUEL) = XA1
            XBFUEL1(NFUEL) = XB1
            KFUEL2(NFUEL)  = K2
            XAFUEL2(NFUEL) = XA2
            XBFUEL2(NFUEL) = XB2
         ELSE
            ZK = ZLEADI(K2) - ZLEADI(K1)
            DO K = K1, K2 - 1
               NFUEL = NFUEL + 1
               KFUEL1(NFUEL) = K
               KFUEL2(NFUEL) = K + 1
               R1 = (ZLEADI(K)   - ZLEADI(K1)) / ZK
               R2 = (ZLEADI(K+1) - ZLEADI(K1)) / ZK
               XAFUEL1(NFUEL) = (ONE - R1) * XA1 + R1 * XA2
               XBFUEL1(NFUEL) = (ONE - R1) * XB1 + R1 * XB2
               XAFUEL2(NFUEL) = (ONE - R2) * XA1 + R2 * XA2
               XBFUEL2(NFUEL) = (ONE - R2) * XB1 + R2 * XB2
            END DO
         END IF

      END DO

      CLOSE (IREAD)
      WRITE (IWRIT,1040) 'All inputs read successfully.'

      IF (NFUEL > 0)
     >   WRITE (IWRIT,'(//,A,//,A,/,(2(1X,I6,2F10.6)))')
     >      ' Expanded fuel volume specification:',
     >      ' KFUEL1   XAFUEL1   XBFUEL1 KFUEL2   XAFUEL2   XBFUEL2',
     >      (KFUEL1(M), XAFUEL1(M), XBFUEL1(M),
     >       KFUEL2(M), XAFUEL2(M), XBFUEL2(M), M = 1, NFUEL)


C     Set up any design variable perturbations to be stored in unit form:
C     -------------------------------------------------------------------

      IF (NDV > 0) THEN

         CALL SECOND (CPU1)
         NCALL = -4

         CALL PERTURB (NCALL)

         CALL SECOND (CPU2)

         CALL CPUTIME (CPU1, CPU2, PROGNAME,
     >                 'to establish unit perturbations', IWRIT)

      END IF


C     Set up the linear constraint matrix and corresponding bounds?
C     -------------------------------------------------------------

      IF (NCLIN > 0 .AND. NDV > 0) THEN

         CALL SETLCON (BIGBND)

      END IF

      NCALL = -3 ! First call

  100 CONTINUE ! Come back here for CHECKWARP case

C     Apply any design variables and/or smooth the wing sections.
C     -----------------------------------------------------------

C     PERTURB is needed to denormalize the wing even if NDV = 0:

      CALL SECOND (CPU1)

      CALL PERTURB (NCALL)

      CALL SECOND (CPU2)

      IF (NEWGEOM) THEN

         CALL CPUTIME (CPU1, CPU2, PROGNAME,
     >                 'to perturb the geometry', IWRIT)

C        Save the perturbed (or just smoothed) geometry:

         OPEN (UNIT=LWING, FILE='ch_grid.wng', STATUS='UNKNOWN')

         CALL GEOMWR (LWING, PROGNAME)

         CLOSE (LWING)

      END IF


C     Summarize the (possibly perturbed) wing sections.
C     -------------------------------------------------

      WRITE (IWRIT,1130)
      TMAXF  = ZERO
      ALMAXF = ZERO

      DO K = 1, NWSEC
         NL     = NLOWER(K)
         NU     = NUPPER(K)
         N      = NTOTAL(K)
         TMAXF  = MAX (THMAX(K), TMAXF)
         ALMAXF = MAX (AWORK(K), ALMAXF)
         COSTWT = COS (ALG(K) * DEG2RAD)
         SINTWT = SIN (ALG(K) * DEG2RAD)
         XTEDGE = (XPRTB(1,K) + XPRTB(N,K)) * HALF
         YTEDGE = (YPRTB(1,K) + YPRTB(N,K)) * HALF
         XTRAIL = XLEADW(K) + CWORK(K) *
     >            ((XTEDGE - XPRTB(NL,K)) * COSTWT +
     >             (YTEDGE - YPRTB(NL,K)) * SINTWT * TWORK(K))
         YTRAIL = YLEADW(K) + CWORK(K) *
     >            ((YTEDGE - YPRTB(NL,K)) * COSTWT * TWORK(K) -
     >             (XTEDGE - XPRTB(NL,K)) * SINTWT)
         EFFTWT = ATAN ((YLEADW(K) - YTRAIL) / (XTRAIL - XLEADW(K)))
     >            * RAD2DEG

         CALL GETANGLE (NL, NU, XPRTB(1,K), YPRTB(1,K), TWORK(K),
     >                  TEANGL)

         WRITE (IWRIT,1140) K, XLEADW(K), YLEADW(K), ZLEADW(K),
     >      CWORK(K), XTRAIL, YTRAIL, TWORK(K), THMAX(K), XTHMAX(K),
     >      AWORK(K), ALG(K), EFFTWT, TEANGL

      END DO


C     Check any linear constraints:
C     -----------------------------

      IF (NCLIN > 0) THEN

         CALL SECOND (CPU1)

         CALL CHKLCON ()

         CALL SECOND (CPU2)

         CALL CPUTIME (CPU1, CPU2, PROGNAME,
     >                 'to check the linear constraints', IWRIT)
      END IF


C     Evaluate any nonlinear constraints:
C     -----------------------------------

      IF (NCNLN > 0) THEN

         CALL SECOND (CPU1)

         CALL NLCON (.TRUE.)

         CALL SECOND (CPU2)

         CALL CPUTIME (CPU1, CPU2, PROGNAME,
     >                 'to evaluate nonlinear constraints', IWRIT)
      END IF


C     Generate the grid?
C     ------------------

      IF (.NOT. PERTURBONLY) THEN

         CALL SECOND (CPU1)

         CALL GRIDWB (NCALL, FAIL)

         IF (FAIL) GO TO 999

         CALL SECOND (CPU2)

         CALL CPUTIME (CPU1, CPU2, PROGNAME,
     >                 'to generate the grid', IWRIT)

C        Option to check "warp" mode by perturbing the geometry with V(1)
C        (or keeping it the same if VSCALE(1) = 0) and redoing the grid:

         IF (CHECKWARP .AND. NCALL == -3) THEN

            V(1) = V(1) + AITCH(1)

            IF (VSCALE(1) == ZERO) THEN
               NCALL = 1 ! Gradient-type - suppress flap/slat angle stuff
            ELSE
               NCALL = 0 ! Allow capturing of flap/slat angles after WARPing
            END IF

            GO TO 100 ! We have NNGRID = 2 above to allow warping

         END IF

C        Save the grid (LUNERR = 1 is unused).

         CALL GRIDIO (1, IL, JL, KL, X, LUNXYZ, 1)

      END IF


C     The fuel volume calculation may be affected by artificial crank stations.
C     Therefore, reregularize the defining sections - don't use those from the
C     grid generation.

      IF (NFUEL > 0) THEN

C        Regularize relevant wing sections with a common distribution:

         CALL WREGUL (KFUEL1(1), KFUEL2(NFUEL), MXIWING, MXKWING,
     >                XWG, YWG, ZWG, ROUNDED, NLG, NW, IL, ILE,
     >                NHALF, SNORM, XRG, YRG, ZRG)

         CALL FUELVOL (FUELVL)

      ELSE

         FUELVL = ZERO

      END IF


      WRITE (IWRIT,1200) THMAX(1), THMAX(NWSEC), TMAXF, AWORK(1),
     >   AWORK(NWSEC), ALMAXF, WINGVL, FUELVL

  999 CONTINUE

      CALL SECOND (TIME2)

      CALL CPUTIME (TIME1, TIME2, PROGNAME, 'for this run', IWRIT)

      STOP

C-----------------------------------------------------------------------
C     Formats.
C-----------------------------------------------------------------------

 1001 FORMAT (A)
 1010 FORMAT (1X, A)
 1040 FORMAT (/, ' CH_GRID: ', A)
 1130 FORMAT (/, ' Wing section summary (including YLESHIFT):', //,
     >        ' Section    Xle        Yle        Zle      ',
     >        'Chord        Xte        Yte    Thick  T/C Max',
     >        '   at X/C Lo.Twist To.Twist Ef.Twist TE.Ang.', /)
 1140 FORMAT (I4, 6F11.4, F9.4, F9.6, 4F9.4, F8.3)
 1200 FORMAT (//, ' GEOMETRIC PARAMETERS', //, 1P,
     >   ' Thickness root     ', E17.8, /,
     >   ' Thickness tip      ', E17.8, /,
     >   ' Thickness maximum  ', E17.8, /,
     >   ' Twist root         ', E17.8, /,
     >   ' Twist tip          ', E17.8, /,
     >   ' Twist maximum      ', E17.8, /,
     >   ' Wing volume        ', E17.8, /,
     >   ' Fuel volume        ', E17.8, //,
     >   ' CH_GRID: Normal termination.')

      END PROGRAM CH_GRID

C***********************************************************************
C
      SUBROUTINE CABINR (X, XB4, MXJBODY, NJBODY, NFSTN, XF, YF, ZF,
     >                   NFLOOR, XFLOOR, YFLOOR, RPOLY, THPOLY,
     >                   TCABIN, YCABIN, ZCABIN, TINT, RMOD, RADIUS,
     >                   IWRIT, ICLASS, IER)
C
C     CABINR calculates the actual cabin RADIUS corresponding to the
C     given "polygon" point (RPOLY, THPOLY) with origin on the center
C     line (cabin floor if RMOD = 0. else camber line if RMOD = 1.).
C
C     It handles interpolating a body station at the specified X if
C     X /= XB4, and the 2-space line/line intersection calculation.
C     IWRIT < 0 suppresses tabulation of results.
C
C     11/15/96  D.Saunders  Modularization of nonlinear constraint
C                           code for use in setting up and checking
C                           linear constraint application.
C
C***********************************************************************

C     Global constants:

      USE CONSTS

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   MXJBODY, NJBODY, NFSTN, NFLOOR, IWRIT, ICLASS, IER
      REAL
     >   X, XB4, XF(NFSTN), YF(MXJBODY,NFSTN), ZF(MXJBODY,NFSTN),
     >   XFLOOR(NFLOOR), YFLOOR(NFLOOR), RPOLY, THPOLY, TCABIN(NJBODY),
     >   YCABIN(NJBODY), ZCABIN(NJBODY), TINT, RMOD, RADIUS

C     Local constants:

      REAL, PARAMETER :: HALF = 0.5, ONE = 1., TOLER = 1.E-6, ZERO = 0.

C     Local variables:

      INTEGER
     >   J, JCABIN, JSPOKE, LEFT, LUNERR
      REAL
     >   TSPOKE(2), YSPOKE(2), ZSPOKE(2), DELY, DELZ, R, THETA,
     >   THSCALE, THSHIFT, YINT, ZINT

C     Execution:

C     Since the regularized body sections have uniform points along the arc,
C     we can transform angles in [-90, 90] to the interval [1, NJBODY] to
C     get reasonable starting indices for the body section curves:

      THSCALE = REAL (NJBODY - 1) / 180.
      THSHIFT = REAL (NJBODY + 1) * HALF

      IF (X /= XB4) THEN

C        Interpolate a cabin section at the specified station:

         XB4 = X
         IF (IWRIT > 0) WRITE (IWRIT,1001)

         LEFT = NFSTN / 4

         CALL INTERVAL (NFSTN, XF, X, ONE, LEFT)

         R = (X - XF(LEFT)) / (XF(LEFT+1) - XF(LEFT))

         DO J = 1, NJBODY
            YCABIN(J) = (ONE - R) * YF(J,LEFT) + R * YF(J,LEFT+1)
            ZCABIN(J) = (ONE - R) * ZF(J,LEFT) + R * ZF(J,LEFT+1)
         END DO

C        Parameterize it:

         CALL CHORDS2D (NJBODY, YCABIN, ZCABIN, .FALSE., TCABIN(NJBODY),
     >                  TCABIN)
      END IF

C     Two-point "spoke" from cabin floor to constraint polygon vertex:

      R = (X - XFLOOR(ICLASS)) / (XFLOOR(ICLASS+1) - XFLOOR(ICLASS))
      YSPOKE(1) = (ONE - R) * YFLOOR(ICLASS) + R * YFLOOR(ICLASS+1)
      ZSPOKE(1) = ZERO
      R         = RPOLY ! For printing
      THETA     = THPOLY
      YSPOKE(2) = R * SIN (THETA * DEG2RAD) + YSPOKE(1)
      ZSPOKE(2) = R * COS (THETA * DEG2RAD)

C     Relocate polar coordinates origin from floor to camber line?

      IF (RMOD == ONE) THEN ! Linear constraint case
         YSPOKE(1) = (YCABIN(1) + YCABIN(NJBODY)) * HALF
         DELY      = YSPOKE(2) - YSPOKE(1)
         DELZ      = ZSPOKE(2)
         THETA     = ATAN2 (DELY, DELZ) * RAD2DEG
         R         = SQRT (DELY**2 + DELZ**2)
         RMOD      = R
      END IF

C     Find where the "spoke" intersects the interpolated cabin section.

      JCABIN = NINT (THETA * THSCALE + THSHIFT)
      JSPOKE = 2
      LUNERR = ABS (IWRIT)

      CALL INTSEC2 (NJBODY, YCABIN, ZCABIN, TCABIN, JCABIN,
     >              .FALSE., 2, YSPOKE, ZSPOKE, TSPOKE, JSPOKE,
     >              .TRUE., TOLER, YINT, ZINT, -LUNERR, IER)

      IF (IER /= 0) THEN
         WRITE (LUNERR,'(/,A,I2,A,3F12.4)')
     >      ' CABINR: INTSEC2 IER =', IER, '.  X, R, THETA:',
     >      X, R, THETA
      END IF

      RADIUS = SQRT ((YINT - YSPOKE(1)) ** 2 + ZINT ** 2)

      IF (IWRIT > 0) THEN
         WRITE (IWRIT,1002) X, YSPOKE(1), ICLASS, THETA, R, RADIUS,
     >                  RADIUS - R, YSPOKE(2), ZSPOKE(2), YINT, ZINT
      END IF

C     Normalized circumferential location of intersection needed by SETLCON:

      TINT = (YINT - YCABIN(JCABIN))**2 + (ZINT - ZCABIN(JCABIN))**2
      TINT = (TCABIN(JCABIN) + SQRT (TINT)) / TCABIN(NJBODY)

      RETURN

 1001 FORMAT (/, '    XCABIN    YFLOOR  CLASS     THETA      RMIN     ',
     >        'RBODY  R - RMIN      YMIN      ZMIN     YBODY     ZBODY')
 1002 FORMAT (F10.2, F10.4, I7, 8F10.4)

      END SUBROUTINE CABINR

C***********************************************************************
C
      SUBROUTINE CHKLCON ()
C
C     Check the perturbed geometry against the linear constraints to
C     give a more readable table than is provided by the optimizer.
C
C     06/25/96  DAS  Adaptation of SETLCON.
C     11/16/96   "   Cabin radius constraint type added.
C     11/27/96   "   Arguments are now in /OPT/ for calling from SOLVE.
C     12/09/96   "   Added b = Ax column to table.
C     01/17/96   "   Added T column so spar thickness can be plotted.
C
C***********************************************************************

C     Global quantities:

      USE DIMENS
      USE FUSEL
      USE FUSEL2
      USE LUNS
      USE OPT

C     Local constants.

      LOGICAL,   PARAMETER :: NEW = .TRUE.  ! New data for each call to LCSFIT
      CHARACTER, PARAMETER :: TIGHT*1 = 'M' ! Monotonic fits

C     Local variables:

      CHARACTER FLAG*3

C     Execution:

      WRITE (IWRIT,1010)
      XB4 = -999.

C     Evaluate the linear constraints in the form seen by NPOPT:

      CALL AX (NCLIN, NCLIN, NDV, AMAT, V, CONLIN)

      DO M = 1, NCLIN
         X = XLCON(M)
         K = KLCON(M)

         IF (LCTYPE(M)(1:3) /= 'CAB') THEN
            NUPPK = NUPPER(K)
            NLOWK = NLOWER(K)
            Z     = ZLEADW(K)
            CHORD = CWORK(K)
         END IF

         IF (LCTYPE(M)(1:3) == 'THK') THEN ! Thickness constraint

C           Find the T/C at this X/C:

            CALL LCSFIT (NLOWK, XPRTB(1,K), YPRTB(1,K),
     >                   NEW, TIGHT, 1, X, YLO, YLO)

            CALL LCSFIT (NUPPK, XPRTB(NLOWK,K), YPRTB(NLOWK,K),
     >                   NEW, TIGHT, 1, X, YUP, YUP)

            VALUE = YUP - YLO

         ELSE IF (LCTYPE(M)(1:3) == 'YUP') THEN ! Upper Y/C lower bound

C           Find the Y/C at this X/C:

            CALL LCSFIT (NUPPK, XPRTB(NLOWK,K), YPRTB(NLOWK,K),
     >                   NEW, TIGHT, 1, X, YUP, YUP)

            VALUE = YUP

         ELSE IF (LCTYPE(M)(1:3) == 'YLW') THEN ! Lower Y/C upper bound

            CALL LCSFIT (NLOWK, XPRTB(1,K), YPRTB(1,K),
     >                   NEW, TIGHT, 1, X, YLO, YLO)

            VALUE = YLO

         ELSE IF (LCTYPE(M)(1:3) == 'CAB') THEN ! Cabin radius constraint

            ISEG = K / 100 ! K = KLCON(M) = 100*floor segment + radial angle #
            JANG = K - 100 * ISEG

C           Find the initial perturbed cabin radius at this X and Theta:

            RPOLYMOD = 1.  ! Kludge to signal recentering without another arg.

            CALL CABINR (X, XB4, NJBODY, NJBODY, NFSTN, XF, YF, ZF,
     >                   NFLOOR, XFLOORF, YFLOORF,
     >                   RBODY(JANG,ISEG), THBODY(JANG,ISEG),
     >                   TCABIN, YCABIN, ZCABIN, TINT, RPOLYMOD, RADIUS,
     >                   -IWRIT, ISEG, IER)

            VALUE = RADIUS
         END IF

         BOUNDL = BLLIN(M)
         BOUNDU = BULIN(M)

         IF (VALUE < BOUNDL .OR. VALUE > BOUNDU) THEN
            FLAG = '***'
         ELSE
            FLAG = '   '
         END IF

         IF (LCTYPE(M)(1:3) /= 'CAB') THEN
            CVALUE = CHORD * VALUE
            WRITE (IWRIT,1020) M, LCTYPE(M), K, Z, X, VALUE, CVALUE,
     >                         BOUNDL, BOUNDU, CONLIN(M), FLAG
         ELSE
            WRITE (IWRIT,1030) M, LCTYPE(M), K, 999., X, VALUE, VALUE,
     >                         BOUNDL, BOUNDU, CONLIN(M), FLAG
         END IF

      END DO

      RETURN

 1010 FORMAT (/, ' CHKLCON: Linear constraint check:', //,
     >   '    # LCTYPE   K         Z  X or X/C RADIUS or T/C or Y/C',
     >   '    R or T or Y       LOWER BOUND       UPPER BOUND',
     >   '          AMAT * V')
 1020 FORMAT (1X, I4, 2X, A5, I4, F10.4, F10.6, E21.13, F15.6, 3E18.10,
     >        2X, A3)
 1030 FORMAT (1X, I4, 2X, A5, I4, F10.0, F10.2, F21.4,  F15.6, 3E18.10,
     >        2X, A3)

C     Internal procedure for CHKLCON:

      CONTAINS

      subroutine ax (idim, m, n, a, x, b)

!     Computes b = Ax for m x n matrix A.

!     Arguments:

      integer, intent (in)  :: idim, m, n
      real,    intent (in)  :: a(idim,n), x(n)
      real,    intent (out) :: b(m)

!     Local variables:

      integer  i, j
      real     sum

!     Execution:

      do i = 1, m
         sum = 0.
         do j = 1, n
            sum = a(i,j) * x(j) + sum
         end do
         b(i) = sum
      end do

      end subroutine ax

      END SUBROUTINE CHKLCON

C***********************************************************************
C
      SUBROUTINE CPUTIME (CPU1, CPU2, CALLER, DESCRIPTION, LUN)
C
C     CPUTIME modularizes printing of CPU time between two measures of
C     total CPU time used so far as given by the SECOND utility.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      REAL(KIND=4), INTENT (IN) :: CPU1, CPU2
      CHARACTER,    INTENT (IN) :: CALLER * (*), DESCRIPTION * (*)
      INTEGER,      INTENT (IN) :: LUN

C     Local variables:

      REAL DIFF

C     Execution:

      DIFF = CPU2 - CPU1
      WRITE (LUN, '(/, 1X, 4A, F10.2)')
     >   CALLER, ': CPU secs. ', DESCRIPTION, ':   ', DIFF

      END SUBROUTINE CPUTIME

C***********************************************************************
C
      SUBROUTINE DEFLECT ()
C
C     Perform flap/slat rotations in REAL space.
C
C     08/29/97  DAS  Extracted from PERTURB, which is getting too big.
C     04/21/99   "   A SINA > ZERO test should have been ANGLE > ZERO.
C
C***********************************************************************

C     Global variables:

      USE CONSTS
      USE DIMENS
      USE LUNS
      USE OPT
      USE WNGSRF

C     Local constants:

      REAL, PARAMETER :: HALF = 0.5, ONE = 1.0, TINY = 0.001, ZERO = 0.0

C     Local variables:

      INTEGER  INDICES(MXIWING)
      LOGICAL  FLAP, SHEAR

C     Execution:

      DO N = NCHCKS + NCHCK1 + NCHCK2 + NSINVC + 1, NTWING

         KI = KMIN(N)           ! Inboard end station of slat/flap
         KO = KMAX(N)           ! Outboard ...
         XHI = XWG(NLOWER(KI),KI) + XWMIN(N) * CHORDG(KI) ! Inboard hinge X
         ZHI = ZWG(NLOWER(KI),KI)
         XHO = XWG(NLOWER(KO),KO) + XWMAX(N) * CHORDG(KO) ! Outboard ...
         ZHO = ZWG(NLOWER(KO),KO)
         COSSWEEP = (ZHO - ZHI) / SQRT ((XHO-XHI)**2 + (ZHO-ZHI)**2) ! Hinge

         DY = V(N) * VSCALE(N)  ! Degrees of rotation about hinge (+ve down)
         TANGNT = TAN (DEG2RAD * DY) * COSSWEEP ! Equivalent chordwise angle

         FLAP  = VTYPE(N)(1:4) == 'FLAP'
         SHEAR = VTYPE(N)(5:6) == 'SH'

         IF (.NOT. SHEAR) THEN
            ANGLE = ATAN (TANGNT)
            COSA  = COS (ANGLE)
            SINA  = SIN (ANGLE)
            IF (FLAP) THEN ! Positive is clockwise; in ROTATE2D it isn't
               SINA = -SINA
            END IF
         END IF

C        Apply the flap/slat shear/rotation to all affected wing sections:

         DO K = KMIN(N), KMAX(N)

            NTOTK = NTOTAL(K)
            NUPPK = NUPPER(K)
            NLOWK = NLOWER(K)
            Z     = ZWG(NLOWK,K)

C           Hinge X at this span station:

            XHK = XHI + (Z - ZHI) * (XHO - XHI) / (ZHO - ZHI)

            IF (SHEAR) THEN

               IF (FLAP) THEN

                  DO I = 1, NTOTK
                     XI = MAX (ZERO, XWG(I,K) - XHK)
                     YWG(I,K) = YWG(I,K) - XI * TANGNT
                  END DO

               ELSE ! Sheared slat

                  DO I = 1, NTOTK
                     XI = MAX (ZERO, XHK - XWG(I,K))
                     YWG(I,K) = YWG(I,K) - XI * TANGNT
                  END DO

               END IF

            ELSE

C              Rotate about hinge (X, Ymean) through the adjusted angle

               ILO = NLOWK / 2
               CALL INTERVAL (NLOWK, XWG(1,K), XHK, -ONE, ILO)
               R = (XHK - XWG(ILO+1,K)) / (XWG(ILO,K) - XWG(ILO+1,K))
               YHLO = (ONE - R) * YWG(ILO+1,K) + R * YWG(ILO,K)

               IUP = NUPPK / 2
               CALL INTERVAL (NUPPK, XWG(NLOWK,K), XHK, ONE, IUP)
               IUP = IUP + NLOWK - 1
               R = (XHK - XWG(IUP,K)) / (XWG(IUP+1,K) - XWG(IUP,K))
               YHUP = (ONE - R) * YWG(IUP,K) + R * YWG(IUP+1,K)

               YHK = (YHUP + YHLO) * HALF ! Hinge Y

               IF (FLAP) THEN

                  DO I = 1, ILO
                     XP = XWG(I,K) - XHK
                     YP = YWG(I,K) - YHK
                     XWG(I,K) = XP * COSA - YP * SINA + XHK
                     YWG(I,K) = XP * SINA + YP * COSA + YHK
                  END DO

                  DO I = IUP + 1, NTOTK
                     XP = XWG(I,K) - XHK
                     YP = YWG(I,K) - YHK
                     XWG(I,K) = XP * COSA - YP * SINA + XHK
                     YWG(I,K) = XP * SINA + YP * COSA + YHK
                  END DO

               ELSE ! Slat

                  DO I = ILO + 1, IUP
                     XP = XWG(I,K) - XHK
                     YP = YWG(I,K) - YHK
                     XWG(I,K) = XP * COSA - YP * SINA + XHK
                     YWG(I,K) = XP * SINA + YP * COSA + YHK
                  END DO

               END IF

C              Suppress any points now inside the wing section:

               IF (ANGLE > TINY) THEN ! Lower surface affects upper

                  CALL REORDER (1, NLOWK, XWG(1,K), INDICES, NNEW)

                  IF (NNEW /= NLOWK) THEN

                     NLOWER(K) = NNEW
                     NLG(K)    = NNEW
                     NTOTAL(K) = NNEW + NUPPK - 1
                     NW(K)     = NTOTAL(K)

                     DO I = 1, NNEW
                        J = INDICES(I)
                        XWG(I,K) = XWG(J,K)
                        YWG(I,K) = YWG(J,K)
                        ZWG(I,K) = ZWG(J,K)
                     END DO

                     DO I = 1, NUPPK - 1
                        XWG(I+NNEW,K) = XWG(I+NLOWK,K)
                        YWG(I+NNEW,K) = YWG(I+NLOWK,K)
                     END DO

                  END IF

               ELSE IF (ANGLE < -TINY) THEN ! Upper surface only is affected

                  CALL REORDER (NLOWK, NTOTK, XWG(1,K), INDICES, NNEW)

                  IF (NNEW /= NUPPK) THEN

                     NUPPER(K) = NNEW
                     NTOTAL(K) = NNEW + NLOWK - 1
                     NW(K)     = NTOTAL(K)

                     DO I = 1, NNEW
                        J = INDICES(I)
                        XWG(I,K) = XWG(J,K)
                        YWG(I,K) = YWG(J,K)
                        ZWG(I,K) = ZWG(J,K)
                     END DO

                  END IF

               END IF

            END IF

         END DO ! Next section

      END DO ! Next slat/flap design variable

      END SUBROUTINE DEFLECT

C*******************************************************************************
C
      SUBROUTINE DYNAMIC ()
C
C     Allocate (most of) the workspace that is dependent upon the grid size or
C     the number of design variables and constraints.
C
C     04/30/98  DAS  Initial implementation (ignoring geometry variables).
C
C*******************************************************************************

C     Global variables:

      USE DIMENS
      USE FACES
      USE LIM
      USE OPT
      USE OUTR
      USE SPACING
      USE TFI
      USE WALL
      USE WNGSRF
      USE XYZ

C     Local variables:

      INTEGER LDM, MAXBND

C     Execution:

C     Module FACES

      ALLOCATE
     >   (DFACEI(3,JL,KL,2,4), DFACEJ(3,IL,KL,2,4), DFACEK(3,IL,JL,2,4))

C     Module OPT

      ALLOCATE (KMIN(NDV), KCEN(NDV), KMAX(NDV), IGROUP(NDV),
     >          KLCON(NCLIN), INLCON(NCNLN), JNLCON(NCNLN))

      MAXBND = NDV + NCLIN + NCNLN

      ALLOCATE (V(NDV),
     >          WBUMP(NDV),    XBUMP(NDV),    XWMIN(NDV),  XWMAX(NDV),
     >          PBUMP(NDV),    UBUMP(NDV),    UBMIN(NDV),  UBMAX(NDV),
     >          DOUP(NDV),     DOLO(NDV),     VSCALE(NDV), AITCH(NDV),
     >          BL(MAXBND),    BU(MAXBND),
     >          XLCON(NCLIN),  BLLIN(NCLIN),  BULIN(NCLIN),
     >          CONLIN(NCLIN), C(NCNLN), XNLCON(NCNLN), SNLCON(NCNLN),
     >          AMAT(NCLIN,NDV))

      ALLOCATE (VTYPE(NDV), UTYPE(NDV), LCTYPE(NCLIN), NLCTYPE(NCNLN))

C     Module OUTR

      ALLOCATE (XLBACK(JL,KL), YLBACK(JL,KL), ZLBACK(JL,KL), XOUT(IL),
     >          XUBACK(JL,KL), YUBACK(JL,KL), ZUBACK(JL,KL), YOUT(IL),
     >          XCROWN(IL),    YCROWN(IL))

C     Module SPACING

      ALLOCATE (ARCOUTER(IL), DRADIAL(KL,3), TRADIAL(JL),
     >          XRADIAL(JL),  YRADIAL(JL),   ZRADIAL(JL),
     >          DRADIALEU(KL,3), DRADIALNS(KL,3))

C     Module TFI

      LDM = MAX (JL, KL)

      ALLOCATE (TFIWORK(2*(IL + LDM))) ! Assumes IL >> LDM

C     Module WALL

      ALLOCATE (UWALL(IL,JL), VWALL(IL,JL),
     >          XWALL(IL,JL), YWALL(IL,JL), ZWALL(IL,JL))

C     Module WNGSRF

      ALLOCATE (IINTR(IL), JINTR(IL), KINTR(IL),
     >          TINTR(IL), UINTR(IL), VINTR(IL),
     >          XRG(IL,MXKWING), YRG(IL,MXKWING), ZRG(IL,MXKWING),
     >          CHORDM(KL), XW(IL,KL), YW(IL,KL), ZW(IL,KL), SNORM(IL),
     >          TBASE(KL), XWLE(KL), YWLE(KL), ZWLE(KL),
     >          XWTE(KL), YWTE(KL))

C     Module XYZ

      ALLOCATE (X(IL,JL,KL,3))

      END SUBROUTINE DYNAMIC

C***********************************************************************
C
      SUBROUTINE FUELVOL (VOLUME)
C
C     05/24/96 DAS/JJR Initial form of fuel volume calculation, using
C                      the raw, normalized geometry.
C     05/30/96    "    Trapezoidal rule is too crude - go to Simpson's
C                      rule.  Use of the regularized sections helps, but
C                      could be invalid if section 1 is involved - it
C                      may now be the wing/body intersection.
C
C***********************************************************************

C     Global variables:

      USE FUEL
      USE INDEX
      USE LIM
      USE LUNS
      USE WNGSRF

      IMPLICIT NONE

C     Arguments:

      REAL, INTENT (OUT) :: VOLUME

C     Local constants:

      REAL,      PARAMETER :: HALF = 0.5, SIXTH = 1./ 6., ZERO = 0.
      CHARACTER, PARAMETER :: METHOD * 1 = 'B' ! Bessel (loose) fit

C     Local variables:

      INTEGER    I, KLEFT, KRGHT, L
      REAL       ALOWER, AUPPER,
     >           ARLEFT, ARMIDL, ARRGHT, CHLEFT, CHMIDL, CHRGHT,
     >           XALEFT, XAMIDL, XARGHT, XBLEFT, XBMIDL, XBRGHT,
     >           XLLEFT, XLMIDL, XLRGHT
      REAL       XMIDL(ITL:ITU), YMIDL(ITL:ITU)

C     Execution:

      VOLUME = ZERO
      DO L = 1, NFUEL

C        Partial area of left section segment:

         KLEFT  = KFUEL1(L)
         XLLEFT = XRG(ILE,KLEFT)
         CHLEFT = XRG(ITL,KLEFT) - XLLEFT
         XALEFT = XAFUEL1(L) * CHLEFT + XLLEFT
         XBLEFT = XBFUEL1(L) * CHLEFT + XLLEFT

         CALL LCSQUAD (NHALF, XRG(ILE,KLEFT), YRG(ILE,KLEFT),
     >                 XALEFT, XBLEFT, METHOD, AUPPER)

         CALL LCSQUAD (NHALF, XRG(ITL,KLEFT), YRG(ITL,KLEFT),
     >                 XALEFT, XBLEFT, METHOD, ALOWER)

         ARLEFT = AUPPER - ALOWER

C        Repeat for right section segment:

         KRGHT  = KFUEL2(L)
         XLRGHT = XRG(ILE,KRGHT)
         CHRGHT = XRG(ITL,KRGHT) - XLRGHT
         XARGHT = XAFUEL2(L) * CHRGHT + XLRGHT
         XBRGHT = XBFUEL2(L) * CHRGHT + XLRGHT

         CALL LCSQUAD (NHALF, XRG(ILE,KRGHT), YRG(ILE,KRGHT),
     >                 XARGHT, XBRGHT, METHOD, AUPPER)

         CALL LCSQUAD (NHALF, XRG(ITL,KRGHT), YRG(ITL,KRGHT),
     >                 XARGHT, XBRGHT, METHOD, ALOWER)

         ARRGHT = AUPPER - ALOWER

C        Loft a section half-way between:

         DO I = ITL, ITU
            XMIDL(I) = (XRG(I,KLEFT) + XRG(I,KRGHT)) * HALF
            YMIDL(I) = (YRG(I,KLEFT) + YRG(I,KRGHT)) * HALF
         END DO

         XAMIDL = (XALEFT + XARGHT) * HALF
         XBMIDL = (XBLEFT + XBRGHT) * HALF

         CALL LCSQUAD (NHALF, XMIDL(ILE), YMIDL(ILE),
     >                 XAMIDL, XBMIDL, METHOD, AUPPER)

         CALL LCSQUAD (NHALF, XMIDL(ITL), YMIDL(ITL),
     >                 XAMIDL, XBMIDL, METHOD, ALOWER)

         ARMIDL = AUPPER - ALOWER

         VOLUME = SIXTH * (ZRG(ILE,KRGHT) - ZRG(ILE,KLEFT)) *
     >            (ARRGHT + 4.* ARMIDL + ARLEFT) + VOLUME
      END DO

      END SUBROUTINE FUELVOL

C***********************************************************************
C
      SUBROUTINE FUSL ()
C
C     FUSL reads the fuselage geometry.
C
C***********************************************************************

C     Global variables:

      USE DIMENS
      USE FUSEL
      USE FUSEL2
      USE LUNS
      USE OPTIONS

C     Local variables:

      REAL
     >   YBUF(MXJBODY), ZBUF(MXJBODY), YP(MXJBODY), ZP(MXJBODY),
     >   YCRN(MXIBODY), YKEL(MXIBODY), ZWAT(MXIBODY)

C     Execution:

      READ (LGEOM,*)
      READ (LGEOM,*) FNF, FLIP
      NFSTN = FNF

      IF (NFSTN > MXIBODY) THEN
         WRITE (IWRIT,*) 'FUSL:  Too many stations: ', NFSTN, MXIBODY
         STOP
      END IF

      DO I = 1, NFSTN

         YCRN(I) = -99999.
         YKEL(I) = -YCRN(I)
         ZWAT(I) =  YCRN(I)

         READ (LGEOM,*)
         IF (FLIP == 0.) THEN
            READ (LGEOM,*) FN, XF(I), FSEC
         ELSE
            READ (LGEOM,*) XF(I), FN, FSEC
         END IF

         N = FN
         NF(I) = N

         IF (FSEC > 0.) THEN ! New fuselage section

            IF (N > MXJBODY) THEN
               WRITE (IWRIT,*) 'Too many fuselage points at station ', I
               WRITE (IWRIT,*) 'Input: ', N, '   Limit: ', MXJBODY
               STOP
            END IF

            READ (LGEOM,*)

            IF (FLIP == 0.) THEN
               DO J = 1, N
                  READ (LGEOM,*) ZP(J), YP(J)
               END DO
            ELSE
               DO J = 1, N
                  READ (LGEOM,*) YP(J), ZP(J)
               END DO
            END IF

C           Ensure that fuselage points go from keel to crown.

            IF (N > 1) THEN
               IF (YP(1) < YP(N)) THEN
                  YBUF(1:N) = YP(1:N)
                  ZBUF(1:N) = ZP(1:N)
               ELSE
                  YBUF(1:N) = YP(N:1:-1)
                  ZBUF(1:N) = ZP(N:1:-1)
               END IF
            ELSE
               YBUF(1) = YP(1)
               ZBUF(1) = ZP(1)
            END IF

         END IF

         DO J = 1, N
            YFINIT(J,I) = YBUF(J)
            ZFINIT(J,I) = ZBUF(J)
            YCRN(I) = MAX (YBUF(J), YCRN(I))
            YKEL(I) = MIN (YBUF(J), YKEL(I))
            ZWAT(I) = MAX (ZBUF(J), ZWAT(I))
         END DO

      END DO

      WRITE (IWRIT,'(//,A,//,A,/)')
     >   ' Input fuselage station summary:',
     >   ' Station        X     Y Crown      Y Keel       Z Max'
      WRITE (IWRIT,'(I4,1X,4F12.5)')
     >   (I, XF(I), YCRN(I), YKEL(I), ZWAT(I), I = 1, NFSTN)

      END SUBROUTINE FUSL

C***********************************************************************
C
      SUBROUTINE GEOMWR (LUN, PROGNAME)
C
C     GEOMWR writes the wing/body geometry to the indicated logical unit.
C     Any floor angle and cabin radii data are also appended.
C
C***********************************************************************

C     Global variables:

      USE CONSTS
      USE DIMENS
      USE FUSEL
      USE FUSEL2
      USE LOGIC
      USE LUNS
      USE OPT
      USE OPTIONS
      USE SIZES
      USE WNGSRF

      IMPLICIT NONE

C     Arguments:

      INTEGER,   INTENT (IN) :: LUN
      CHARACTER, INTENT (IN) :: PROGNAME * (*)

C     Local constants:

      REAL, PARAMETER :: ONE = 1., ZERO = 0.

C     Local variables:

      INTEGER  I, J, K, NLOWK, NTOTK
      REAL     CA, SA, FNFP, FNU, FNL, YLETEMP,
     >         XTEMP(MXIWING), YTEMP(MXIWING), ZTEMP(MXIWING)
      LOGICAL  THREECOLS

C     Execution:

C     Write the current geometry in single-block design code format:

      WRITE (LUN,1002) 'Geometry from ', PROGNAME,
     >   '     ZSPAN      SREF      CREF      XREF      YREF      ZREF'
      WRITE (LUN,1006) ZSPAN, SREF, CREF, XREF, YREF, ZREF
      WRITE (LUN,1001) '       FNC    ALPHAW'
      WRITE (LUN,1030) FNWSEC, ALPHAW

      DO K = 1, NWSEC

         THREECOLS = WINGIN >= 3 .OR. (K == 1 .AND. NNROOT == 1)
         NTOTK = NTOTAL(K)
         NLOWK = NLOWER(K)
         FNL   = NLOWK
         FNU   = NUPPER(K)

         WRITE (LUN,1002) '         Z       XLE       YLE',
     >      '     CHORD     THICK     TWIST       NEW     ROUND'

         IF (WINGOUT == 0) THEN ! Sections are untwisted and normalized

            YLETEMP = YLEADW(K) - YLESHIFT

            WRITE (LUN,1030) ZLEADW(K), XLEADW(K), YLETEMP, CWORK(K),
     >                       TWORK(K), AWORK(K), ONE, ROUNDED(K)
            WRITE (LUN,1010) ZERO, FNU, FNL

            IF (THREECOLS) THEN
               WRITE (LUN,1005) (XPRTB(I,K), YPRTB(I,K), ZPRTB(I,K),
     >                           I = NLOWK, NTOTK)
               WRITE (LUN,1011)
               WRITE (LUN,1005) (XPRTB(I,K), YPRTB(I,K), ZPRTB(I,K),
     >                           I = NLOWK, 1, -1)
            ELSE
               WRITE (LUN,1004) (XPRTB(I,K), YPRTB(I,K),
     >                           I = NLOWK, NTOTK)
               WRITE (LUN,1011)
               WRITE (LUN,1004) (XPRTB(I,K), YPRTB(I,K),
     >                           I = NLOWK, 1, -1)
            END IF

         ELSE IF (WINGOUT == 1) THEN ! Twisted & normalized

            CA = AWORK(K) * DEG2RAD
            SA = SIN (CA)
            CA = COS (CA)

            DO I = 1, NTOTK ! Twist about [0.,0.]
               XTEMP(I) = XPRTB(I,K) * CA + YPRTB(I,K) * SA
               YTEMP(I) = YPRTB(I,K) * CA - XPRTB(I,K) * SA
            END DO

            YLETEMP = YLEADW(K) - YLESHIFT

            WRITE (LUN,1030) ZLEADW(K), XLEADW(K), YLETEMP, CWORK(K),
     >                       TWORK(K), ZERO, ONE, ROUNDED(K)
            WRITE (LUN,1010) ZERO, FNU, FNL

            IF (THREECOLS) THEN
               WRITE (LUN,1005) (XTEMP(I), YTEMP(I), ZPRTB(I,K),
     >                           I = NLOWK, NTOTK)
               WRITE (LUN,1011)
               WRITE (LUN,1005) (XTEMP(I), YTEMP(I), ZPRTB(I,K),
     >                           I = NLOWK, 1, -1)
            ELSE
               WRITE (LUN,1004) (XTEMP(I), YTEMP(I), I = NLOWK, NTOTK)
               WRITE (LUN,1011)
               WRITE (LUN,1004) (XTEMP(I), YTEMP(I), I = NLOWK, 1, -1)
            END IF

         ELSE !  WINGOUT = 2  (Fully denormalized)

            DO I = 1, NTOTK
               YTEMP(I) = YWG(I,K) - YLESHIFT
            END DO

            WRITE (LUN,1030) ZLEADW(K), ZERO, ZERO, ONE, ONE, ZERO, ONE,
     >                       ROUNDED(K)
            WRITE (LUN,1010) ZERO, FNU, FNL

            IF (THREECOLS) THEN
               WRITE (LUN,1008) (XWG(I,K), YTEMP(I), ZWG(I,K),
     >                           I = NLOWK, NTOTK)
               WRITE (LUN,1011)
               WRITE (LUN,1008) (XWG(I,K), YTEMP(I), ZWG(I,K),
     >                           I = NLOWK, 1, -1)
            ELSE
               WRITE (LUN,1007) (XWG(I,K), YTEMP(I), I = NLOWK, NTOTK)
               WRITE (LUN,1011)
               WRITE (LUN,1007) (XWG(I,K), YTEMP(I), I = NLOWK, 1, -1)
            END IF

         END IF
      END DO

      IF (.NOT. WINGONLY) THEN
         WRITE (LUN,1001) '       FNF      FLIP'
         WRITE (LUN,1013) FNF, ZERO

         DO I = 1, NFSTN
            WRITE (LUN,1001) '        FNFP          XF        FSEC'
            FNFP = NF(I)
            WRITE (LUN,1012) FNFP, XF(I), ONE
            WRITE (LUN,1001) '        Z(J)        Y(J)'
            WRITE (LUN,1014) (ZFPRTB(J,I), YFPRTB(J,I), J = 1, NF(I))
         END DO

C        Option to transmit cabin floor and radius data:
C        -----------------------------------------------

         IF (NNFLOOR > 0) THEN
            WRITE (LUN,1001)
     >         '# CABIN FLOOR BREAK POINTS (= 1 + # FLOOR SEGMENTS)'
            WRITE (LUN,'(I3)') NFLOOR

            WRITE (LUN,1001) '   XFLOOR   YFLOOR'
            WRITE (LUN,'(2F9.2)') (XFLOORF(I), YFLOORF(I), I = 1,NFLOOR)

            DO I = 1, NFLOOR - 1
               WRITE (LUN,'(A,I2)') '# CABIN RADII FOR FLOOR SEGMENT', I
               WRITE (LUN,'(I3)') NRADII(I)

               WRITE (LUN,1001) '    RBODY   THBODY'
               WRITE (LUN,'(2F9.4)')
     >            (RBODY(J,I), THBODY(J,I), J = 1, NRADII(I))
            END DO
         END IF
      END IF

      RETURN

C     Formats:

 1001 FORMAT (A)
 1002 FORMAT (A, A)
 1004 FORMAT (2F12.8)
 1005 FORMAT (3F12.8)
 1006 FORMAT (F10.3,F10.1,4F10.3)
 1007 FORMAT (2F15.5)
 1008 FORMAT (3F15.5)
 1009 FORMAT (F10.0, 2I10)
 1010 FORMAT ('      YSYM        NU        NL', /,
     >        3F10.2, /, ' Upper surface')
 1011 FORMAT (' Lower surface')
 1012 FORMAT (3F12.6)
 1013 FORMAT (2F10.2)
 1014 FORMAT (2F12.6)
 1030 FORMAT (8F10.4)

      END SUBROUTINE GEOMWR

C***********************************************************************
C
      SUBROUTINE GETANGLE (NL, NU, X, Y, THICK, TEANGL)
C
C     Trailing edge angle for wrap-around airfoil (in degrees).
C
C***********************************************************************

      USE TRIGD

C     Arguments:

      INTEGER, INTENT (IN)  :: NL, NU
      REAL,    INTENT (IN)  :: X(*), Y(*), THICK
      REAL,    INTENT (OUT) :: TEANGL

C     Constants:

      REAL, PARAMETER :: ARROW = 1., FRACTION = 0.99

C     Execution:

C     Pick points at ~99% X/C:

      IUP = NU - 5
      CALL INTERVAL (NU, X(NL), FRACTION, ARROW, IUP)
      IUP = IUP + NL - 1

      ILO = 3
      CALL INTERVAL (NL, X(1), FRACTION, -ARROW, ILO)
      ILO = ILO + 1

      N = NL + NU - 1
      TEANGL = ATAND (THICK*(Y(1) - Y(ILO)) / (X(1) - X(ILO))) +
     >         ATAND (THICK*(Y(IUP) - Y(N)) / (X(N) - X(IUP)))

      END SUBROUTINE GETANGLE

C***********************************************************************
C
      SUBROUTINE GETTHICK (NL, NU, X, Y, YEVAL, THMAX, XTHMAX)
C
C     Wing section maximum thickness calculation for wrap-around data:
C     interpolate the lower surface at the upper surface Xs, then pick
C     the point of maximum thickness (returned in the input units).
C     All dimensions should be at least NU.
C
C     05/28/96  D.Saunders  Adaptation of familiar stuff.
C     11/27/96      "       Switched to monotonic fits to match SETLCON.
C     ??? Some day:         Do it more thoroughly for gradients.
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN)  :: NL, NU
      REAL,    INTENT (IN)  :: X(*), Y(*)
      REAL,    INTENT (OUT) :: YEVAL(*), THMAX, XTHMAX

C     Local variables:

      INTEGER  I
      REAL     THICK

C     Execution:

C     Interpolate the lower surface at the upper surface Xs:

      CALL LCSFIT (NL, X, Y, .TRUE., 'M', NU, X(NL), YEVAL(NL),
     >             YEVAL(NL))

C     Pick the point with maximum thickness:

      THMAX = 0.
      XTHMAX = X(NL)

      DO I = NL, NL + NU - 1
         THICK = Y(I) - YEVAL(I)
         IF (THICK > THMAX) THEN
            THMAX  = THICK
            XTHMAX = X(I)
         END IF
      END DO

      END SUBROUTINE GETTHICK

C***********************************************************************
C
      SUBROUTINE GRIDINP ()
C
C     GRIDINP reads most of the grid generation inputs.
C
C***********************************************************************

C     Global variables:

      USE DMAX
      USE ELLIPTIC
      USE INDEX
      USE LIM
      USE LUNS
      USE OPTIONS
      USE SIZES
      USE SPACING
      USE WNGSRF

C     Execution:

      READ (IREAD,*)
      READ (IREAD,*)
      READ (IREAD,*) KWING1, KWING2, IBOD1, IBOD2, IDEGBOD,
     >               NBLAYEREU, RBLAYEREU

      ILE   = (IL + 1) / 2    ! I at leading edge
      NHALF = (NWING + 1) / 2 ! # pts on one surface of each wing grid section
      I     = (IL - NWING) / 2
      ITU   = IL - I          ! I at trailing edge, upper & lower
      ITL   = 1  + I
      I     = (IL - NBODY) / 2
      IBU   = IL - I          ! I at body tail, upper & lower
      IBL   = 1  + I

      IFACE(1) = 1
      IFACE(2) = ITL
      IFACE(3) = ILE
      IFACE(4) = ITU
      IFACE(5) = IL

      KFACE(1) = 1
      KFACE(2) = KTIP
      KFACE(3) = KL

      READ (IREAD,*)
      READ (IREAD,*) CRDSPL, CRDSPQ, CRDSPS, CRDSPC, SPNSPC, SPNSPS
      READ (IREAD,*)
      READ (IREAD,*) D1ROOTEU, D2ROOTEU, D3ROOTEU, TEROOT
      READ (IREAD,*)
      READ (IREAD,*) D1TIPEU,  D2TIPEU,  D3TIPEU,  TETIP
      READ (IREAD,*)
      READ (IREAD,*) D1KMAXEU, D2KMAXEU, D3KMAXEU, TEKMAX
      READ (IREAD,*)
      READ (IREAD,*) D1CROWNEU, D1NOSEEU, D1TAILEU, DWTAILEU, DCTAILEU
      READ (IREAD,*)
      READ (IREAD,*) OSPNOSE, OSPTAIL, OSPSYM, OSPTIP, D2JLEU
      READ (IREAD,*)
      READ (IREAD,*) SWEEP, DIHED, BSWEEP, YLESHIFT, YSHELF, ZSHELF
      READ (IREAD,*)
      READ (IREAD,*) NCRANK, (ZCRANK(K), K = 1, MIN (MXCRANK, NCRANK))
      READ (IREAD,*)

      K = NCRANK
      IF (K >= MXCRANK) THEN  ! Have to allow for adding one at the tip
         NCRANK = MXCRANK - 1
         WRITE (IWRIT,'(/,(1X,A,I4))')
     >      'Too many wing cranks specified:', K,
     >      'Limit:', NCRANK, 'Proceeding ...'
      END IF

      IF (NCRANK > 0) THEN ! These are needed for the READGRID option
         READ (IREAD,*) (MCRANK(K), K = 1, NCRANK)
         READ (IREAD,*) (LCRANK(K), K = 1, NCRANK)
      ELSE
         READ (IREAD,*)
         READ (IREAD,*)
      END IF

      READ (IREAD,*)
      READ (IREAD,*) XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, RADIUS

      READ (IREAD,*) ! Radial redistribution controls
      READ (IREAD,*)
      READ (IREAD,*) NBLAYERNS, RBLAYERNS
      READ (IREAD,*)
      READ (IREAD,*) D1ROOTNS, D2ROOTNS, D3ROOTNS
      READ (IREAD,*)
      READ (IREAD,*) D1TIPNS,  D2TIPNS,  D3TIPNS
      READ (IREAD,*)
      READ (IREAD,*) D1KMAXNS, D2KMAXNS, D3KMAXNS
      READ (IREAD,*)
      READ (IREAD,*) D1CROWNNS, DWTAILNS, DCTAILNS, D2JLNS

      D1NOSENS = D1NOSEEU ! Axial distributions (nose crown/keel) ...
      D1TAILNS = D1TAILEU ! ... (and tail water line) must not change

      READ (IREAD,*) ! Elliptic smoothing controls

      ITFLOAT  = 0
      ITFREEZE = 1
      JLAYER   = 4
      URFLOAT  = 0.1
      DMAX2D   = D1TIPEU ! (?)
      DMAX3D   = D1TIPEU

      READ (IREAD,*)
      DO I = 1, NCASE2D
         READ (IREAD,*)
     >      ITMAX2D(I),  OMG2D(I),    CONV2D(I),   CONVMIN(I),
     >      EXPI12D(I),  EXPI22D(I),  EXPJ12D(I),  EXPJ22D(I),
     >      POWERI(I),   POWERJ(I),   BGMODE2D(I), FGMODE2D(I),
     >      SPMODE(I)
      END DO
      READ (IREAD,*)
      READ (IREAD,*)
      DO I = 1, NCASE3D
         READ (IREAD,*)
     >      ITMAX3D(I), OMG3D(I),   CONV3D(I),
     >      EXPI13D(I), EXPI23D(I), EXPJ13D(I),
     >      EXPJ23D(I), EXPK13D(I), EXPK23D(I),
     >      BGPHI(I),   BGPSI(I),   BGOMG(I),   FGMODE3D(I)
      END DO
      READ (IREAD,*) ! Miscellaneous elliptic smoothing controls
      READ (IREAD,*)
      READ (IREAD,*) PRINT2D, PRINT3D, FGLIM2D, FGLIM3D, URFG2D, URFG3D

      IF (NNTIP /= 0) THEN
         IF (ABS (YMIN + YMAX) > YMAX * 0.000001) THEN
            WRITE (IWRIT,'(/,A)')
     >         ' GRIDINP: Warning - setting YMIN = -YMAX for C-O tip.'
            YMIN = -YMAX
         END IF
      END IF

      END SUBROUTINE GRIDINP

C***********************************************************************
C
      SUBROUTINE GRIDIO (INOUT, IL, JL, KL, X, LUNXYZ, LUNERR)
C
C     GRIDIO reads or writes a single-block grid in PLOT3D /mgrid format.
C
C***********************************************************************

C     Arguments:

      INTEGER, INTENT (IN)    :: INOUT
      INTEGER, INTENT (INOUT) :: IL, JL, KL
      INTEGER, INTENT (IN)    :: LUNXYZ, LUNERR
      REAL,    INTENT (INOUT) :: X(IL,JL,KL,3)

C     Execution:

      IF (INOUT == 0) THEN ! Read

         READ (LUNXYZ) ! 1 block
         READ (LUNXYZ) IL, JL, KL
         READ (LUNXYZ) X
         CLOSE (LUNXYZ)

      ELSE ! INOUT = 1 (Write)

         REWIND (LUNXYZ)
         WRITE (LUNXYZ) 1
         WRITE (LUNXYZ) IL, JL, KL
         WRITE (LUNXYZ) X

      END IF

      END SUBROUTINE GRIDIO

C***********************************************************************
C
      SUBROUTINE LOFT (TYPE2, Z, ZCEN, ZMIN, ZMAX, PWR, HEIGHT, LUNERR)
C
C     LOFT modularizes the choice between polynomial and Hicks-Henne-type
C     spanwise lofting of the effect of a wing design variable, for use
C     by both PERTURB and SETLCON.  HEIGHT is returned between 0. and 1.
C
C***********************************************************************

C     Arguments:

      CHARACTER, INTENT (IN)  :: TYPE2 * 4
      INTEGER,   INTENT (IN)  :: LUNERR
      REAL,      INTENT (IN)  :: Z, ZCEN, ZMIN, ZMAX, PWR
      REAL,      INTENT (OUT) :: HEIGHT

C     Constants:

      REAL, PARAMETER :: ONE = 1.

C     Execution:

      IF (TYPE2 == 'POLY') THEN

C        Polynomial-type lofting (often linear or constant):
C        N.B.: The constant case requires special treatment at end-pts.

         CALL RIPPER (Z, ZCEN, ZMIN, ZMAX, PWR, HEIGHT)

      ELSE

C        Hicks-Henne-type lofting:

         HEIGHT = 0. ! SFEVAL adds
         ZCNORM = (ZCEN - ZMIN) / (ZMAX - ZMIN) ! Reqd. by SFEVAL

         CALL SFEVAL (TYPE2, PWR, ZCNORM, ZMIN, ZMAX,
     >                ONE, 1, 1, Z, HEIGHT, LUNERR)

      END IF

      END SUBROUTINE LOFT

C***********************************************************************
C
      SUBROUTINE NLCON (PRINT)
C
C     NLCON evaluates all nonlinear constraints, C(*), with the option to
C     tabulate results.  It may be called from CONFUN with PRINT = .FALSE.
C     or from PERTURB at the start of a run and at the end of each design
C     step to give a more readable summary with scaling suppressed.
C
C     5/96  DAS  Initial implementation for monitoring cabin radii.
C     6/96   "   Adapted for use by CONFUN.
C     8/96   "   Added 'CLEAR' constraint (centerline only).
C     11/96  "   CABINR replaces in-line code; 'CLEAR' applies to a
C                specified wing section (but not wing/body intersection
C                because we don't have that prior to grid generation).
C     9/97   "   WREGUL is now fully argument-driven, and does not
C                generate its own chordwise distribution.
C
C***********************************************************************

C     Global quantities:

      USE CONSTS
      USE DIMENS
      USE FUEL
      USE FUSEL
      USE FUSEL2
      USE INDEX
      USE LIM
      USE LUNS
      USE OPT
      USE SPACING
      USE TRIGD
      USE WNGSRF

      IMPLICIT  NONE

C     Arguments:

      LOGICAL, INTENT (IN) :: PRINT ! I .TRUE. suppresses scaling, and tabulates

C     Local constants:

      REAL,    PARAMETER :: ONE = 1., ZERO = 0.
      LOGICAL, PARAMETER :: NEW = .TRUE.

C     Local variables:

      INTEGER
     >   I, IER, ISEG, J, JANG, K, LUNWRIT, M, NL, NU
      REAL
     >   DSLOPEF(MXFLOOR), DSLOPEI(MXFLOOR),
     >   SLOPEF(MXFLOOR),  SLOPEI(MXFLOOR),
     >   CLEAR, FHEIGHT, RADIUS, RPOLYMOD, SCALE, TEANGL, THETA, TINT,
     >   WHEIGHT, X, XB4
      LOGICAL
     >   CLEAR1, FLOOR1
      CHARACTER
     >   TYPE * 6

C     Execution:

      XB4     = -999.  ! First cabin radius constraint
      CLEAR1  = .TRUE. ! First 'CLEAR' constraint?
      FLOOR1  = .TRUE. ! First 'FLOOR' or 'DFLOOR'?

      DO M = 1, NCNLN
         X     = XNLCON(M)
         ISEG  = INLCON(M)
         JANG  = JNLCON(M)
         SCALE = SNLCON(M)
         TYPE  = NLCTYPE(M)

         IF (TYPE == 'CABIN ') THEN ! Cabin radii

            LUNWRIT = -IWRIT
            IF (PRINT) LUNWRIT = IWRIT

            RPOLYMOD = 0. ! Kludge to suppress recentering without another arg.

            CALL CABINR (X, XB4, NJBODY, NJBODY, NFSTN, XF, YF, ZF,
     >                   NFLOOR, XFLOORF, YFLOORF,
     >                   RBODY(JANG,ISEG), THBODY(JANG,ISEG),
     >                   TCABIN, YCABIN, ZCABIN, TINT, RPOLYMOD, RADIUS,
     >                   LUNWRIT, ISEG, IER)

            IF (.NOT. PRINT) C(M) = RADIUS * SCALE ! Scaled constraint value

         ELSE IF (TYPE == 'FUEL  ') THEN ! Fuel volume

C           Regularize relevant wing sections:

            CALL WREGUL (KFUEL1(1), KFUEL2(NFUEL), MXIWING, MXKWING,
     >                   XWG, YWG, ZWG, ROUNDED, NLG, NW, IL, ILE,
     >                   NHALF, SNORM, XRG, YRG, ZRG)

            CALL FUELVOL (FUELVL)

            IF (PRINT) THEN
               WRITE (IWRIT,'(//,A,1P,E16.8)') ' Fuel volume:', FUELVL
            ELSE
               C(M) = FUELVL * SCALE
            END IF

         ELSE IF (TYPE == 'FLOOR ' .OR.  ! Floor angles & angle changes
     >            TYPE == 'DFLOOR') THEN

            IF (FLOOR1) THEN  ! Set up all of the floor angles if any
               FLOOR1 = .FALSE.
               IF (PRINT) WRITE (IWRIT,1003)

               DO I = 1, NFLOOR
                  IF (I < NFLOOR) THEN
                     SLOPEI(I) = ATAN2D (YFLOORI(I) - YFLOORI(I+1),
     >                           XFLOORI(I+1) - XFLOORI(I)) + ALPHAI
                     SLOPEF(I) = ATAN2D (YFLOORF(I) - YFLOORF(I+1),
     >                           XFLOORF(I+1) - XFLOORF(I)) + ALPHAF
                     IF (I > 1) THEN
                        DSLOPEI(I) = SLOPEI(I) - SLOPEI(I-1)
                        DSLOPEF(I) = SLOPEF(I) - SLOPEF(I-1)
                     END IF
                  END IF

                  IF (PRINT) THEN
                     IF (I == 1 .OR. I == NFLOOR) THEN
                        WRITE (IWRIT,'(2F10.2,20X,2F10.2)')
     >                    XFLOORI(I), YFLOORI(I), XFLOORF(I), YFLOORF(I)
                     ELSE
                        WRITE (IWRIT,'(2F10.2,F15.2,5X,2F10.2,F15.2)')
     >                    XFLOORI(I), YFLOORI(I), DSLOPEI(I),
     >                    XFLOORF(I), YFLOORF(I), DSLOPEF(I)
                     END IF

                     IF (I < NFLOOR)
     >                  WRITE (IWRIT,'(20X,F7.2,33X,F7.2)') SLOPEI(I),
     >                    SLOPEF(I)
                  END IF
               END DO
            END IF

            IF (.NOT. PRINT) THEN ! Constraint value (scaled)
               IF (TYPE == 'FLOOR ') THEN
                  C(M) = SLOPEF(ISEG) * SCALE
               ELSE
                  C(M) = DSLOPEF(ISEG) * SCALE
               END IF
            END IF

         ELSE IF (TYPE == 'TCMAX ') THEN ! (Untwisted) wing section T/C max

            C(M) = THMAX(ISEG) * SCALE   ! Calculated in PERTURB

         ELSE IF (TYPE == 'TEANG ') THEN ! Wing section trailing edge angle

            K = ISEG
            CALL GETANGLE (NLOWER(K), NUPPER(K), XPRTB(1,K), YPRTB(1,K),
     >                     TWORK(K), TEANGL)

            C(M) = TEANGL * SCALE

         ELSE IF (TYPE == 'CLEAR ') THEN ! Floor/wing clearance

            CALL LCSFIT (NFLOOR, XFLOORF, YFLOORF, NEW, 'L', 1, X,
     >                   FHEIGHT, FHEIGHT)

            K  = ISEG
            NL = NLG(K)
            NU = NW(K) - NL

            CALL LCSFIT (NU, XWG(NL,K), YWG(NL,K), NEW, 'M', 1, X,
     >                   WHEIGHT, WHEIGHT)

            CLEAR = FHEIGHT - WHEIGHT

            IF (PRINT) THEN
               IF (CLEAR1) THEN
                  WRITE (IWRIT,1004)
                  CLEAR1 = .FALSE.
               END IF
               WRITE (IWRIT,'(2F10.4,3F14.5)')
     >            X, ZWG(1,K), FHEIGHT, WHEIGHT, CLEAR
            ELSE
               C(M) = CLEAR * SCALE
            END IF

         END IF

      END DO

  999 RETURN

 1003 FORMAT (/, 13X, 'INITIAL FLOOR GEOMETRY', 15X,
     >        'CURRENT (INCLUDING ALPHA)', //,
     >        '    XFLOOR    YFLOOR  SLOPE  DSLOPE', 5X,
     >        '    XFLOOR    YFLOOR  SLOPE  DSLOPE', /)
 1004 FORMAT (/, '         X         Z  FLOOR HEIGHT   WING HEIGHT',
     >        '     CLEARANCE', /)

      END SUBROUTINE NLCON

C*******************************************************************************
C
      SUBROUTINE PERTURB (NCALL)
C
C     Application of shape functions/design variables to wing and/or fuselage.
C     See main program for descriptions of the functions and other inputs.
C
C     c. 1992  J.Reuther  Original file-driven scheme.
C     Aug. 95  D.Saunders Vectorized version.
C     Mar. 96   JJR/DAS   Unit perturbations are now stored to facilitate
C                         application of constraints; Hicks-Henne initially.
C     Apr. 96      "      Hicks-Henne fns. apply to [XWMIN, XWMAX] now.
C     May  96      "      Fuselage break points handled (NNFLOOR).
C     June 96     DAS     PERTURB may be called from several places now; had to
C                         do any initializations only once; left NCALL = -4 case
C                         alone (unit bumps) so main prog. can measure CPU time.
C     Aug. 96      "      Introduced SFEVAL (Y & Y" only).
C     Sep. 96   DAS/JJR   Body variables vectorized with SFEVAL; body vertical &
C                         radial variables added to body span variables, and the
C                         variation is two-directional (affecting wing inputs).
C     Jan. 97     DAS     Spanwise lofting of wing perturbations may now be
C                         controlled with Hicks-Henne shape functions, which may
C                         in turn be controlled by Z inputs (REALs). But the
C                         original K-index-based control is still supported.
C     Aug. 97      "      FLAP/FLAPSH and SLAT/SLATSH (rotate and shear pairs)
C                         introduced for arbitrary hinge lines, replacing the
C                         original FLAP/SLAT functions (shear/const. X/C hinge).
C                         SIN*VC design variables (Variable Center family)
C                         introduced for perturbing slats or flaps or for
C                         centering bumps along an arbitrary line such as a line
C                         parallel to a shock off a nacelle/diverter.
C*******************************************************************************

C     Global variables:

      USE CONSTS
      USE DIMENS
      USE FUSEL
      USE FUSEL2
      USE INDEX
      USE LOGIC
      USE LUNS
      USE OPT
      USE OPTIONS
      USE WNGSRF

C     Arguments:

      INTEGER, INTENT (IN) :: NCALL

C     Local constants:

      INTEGER
     >   IFLAG
      REAL
     >   HALF, ONE, TOLER, TWO, ZERO
      LOGICAL
     >   NEW

      PARAMETER
     >  (IFLAG = -999, ! Assigned to KCEN/KMIN/KMAX in RDVLINE if REALs found
     >   HALF  = 0.5,
     >   ONE   = 1.,
     >   TOLER = 1.E-6,
     >   TWO   = 2.,
     >   ZERO  = 0.,
     >   NEW   = .TRUE.) ! New data for each call to LCSFIT

C     Local variables:

      REAL
     >   A(MXIWING), B(MXIWING), CC(MXIWING), DD(MXIWING),
     >   COSTHETA(MXJBODY), SINTHETA(MXJBODY),
     >   RBUMP(MXJBODY), SFUS(MXJBODY), SUNIF(NJBODY),
     >   YENVL2(MXIBODY), YENVU2(MXIBODY), YENVINT(MXIBODY),
     >   YNORM(MXJBODY,MXIBODY), ZNORM(MXJBODY,MXIBODY)
      REAL, ALLOCATABLE, SAVE ::
     >   YUNIT(:,:,:)
      LOGICAL
     >   BODYFIXD, BOTH, INITIZE, NEWBODY, YPP
      LOGICAL, SAVE ::
     >   FIRST = .TRUE.
      CHARACTER
     >   METHOD*1, TYPE*6, TYPE2*4, TYP1*1, TYP4*4

C     Storage:

C     /SCRATCH/ is not used anywhere else - it just ensures plenty of
C     scratch space for passing C to NLCON.

      COMMON /SCRATCH/ CC, A, B, DD, SFUS, SUNIF

C     Shape function statement functions:

C     VBUMP1 & 2 are of the form (AX + B) * <Hicks-Henne SINE function>
C     such that they have zero net area under them on [0, 1].  Both
C     forms are symmetric about XD = 0.5.

      REAL
     >   EM, EN, EA, XD, VBUMP1, VBUMP2

      VBUMP1 (EM, EN, EA, XD) = (EA * (XD - EN) + ONE) *
     >   SIGN (ABS (SIN (TWO * PI * XD **
     >   (LOGPT5 / LOG (EN)))) ** EM,
     >   SIN (TWO * PI * XD ** (LOGPT5 / LOG (EN))))

      VBUMP2 (EM, EN, EA, XD) = (EA * (XD - EN) + ONE) *
     >   SIGN (ABS (SIN (TWO * PI * (ONE - XD) **
     >   (LOGPT5 / LOG (ONE - EN)))) ** EM,
     >   SIN (TWO * PI * XD ** (LOGPT5 / LOG (ONE - EN))))


C     Execution:
C     ----------

      IF (NCALL == -4 .AND. NDV > 0) THEN

C        Allocate as little unit perturbation storage as practical:

         IF (WINGONLY .OR. NDV == NTWING) THEN
            NK = NWSEC
         ELSE
            NK = MAX (NWSEC, MXJBODY)
         END IF

         ALLOCATE (YUNIT(MXIMAX,NK,NDV))

C        Zero out all unit perturbations:

         YUNIT(1:MXIMAX,1:NK,1:NDV) = ZERO

C        Store unit perturbations for the wing design variables.
C        -------------------------------------------------------

         NWVAR = NCHCK1 + NCHCK2 ! Only one can be nonzero
         YPP   = NCHCK2 > 0      ! TRUE means Y" is being perturbed, not Y

C        First, translate INTEGER lofting variables to wing section Zs:

         DO N = 1, NTWING                 ! Includes planform variables & flaps
            IF (KCEN(N) /= IFLAG) THEN    ! KCEN = -999 if REALs were read
               UBUMP(N) = ZLEADI(KCEN(N)) ! Z at shape function peak,
               UBMIN(N) = ZLEADI(KMIN(N)) ! etc.
               UBMAX(N) = ZLEADI(KMAX(N))
            END IF
         END DO

C        Outer loop over wing sections, not variables, is barely more efficient:

         DO K = 1, NWSEC

            NTOTK = NTOTAL(K)
            NUPPK = NUPPER(K)
            NLOWK = NLOWER(K)
            Z     = ZLEADI(K)

C           Determine a multiplier in [0.,1.] at this K for each wing design
C           variable (spanwise lofting), and where appropriate (as for sine
C           bumps) store the corresponding chordwise variation as well.

            DO N = 1, NTWING

               CALL LOFT (UTYPE(N), Z, UBUMP(N), UBMIN(N), UBMAX(N),
     >                    PBUMP(N), DYMAX, IWRIT)

               IF (DYMAX /= ZERO) THEN

                  IF (N <= NCHCKS) THEN ! Planform variable

                     YUNIT(1,K,N) = DYMAX ! Multiplier in (0.,1.]

                  ELSE IF (N <= NCHCKS + NWVAR) THEN ! Wing section Y or Y"
                                                     ! (constant center X/Cs)

C                    Some variables will perturb the leading edge twice if
C                    we don't perturb both surfaces in one call, so:

                     IF (.NOT. YPP) THEN
                        TYPE = VTYPE(N)
                        BOTH = TYPE(1:2) == 'LE' .OR.  ! LED or LEAD
     >                         TYPE(1:2) == 'DR'       ! DROOP
                     ELSE  ! Y" case:  VTYPE = 'Dxxxxx'
                        TYPE = VTYPE(N)(2:)
                        BOTH = TYPE(1:2) == 'LE'       ! LED or LEAD
                     END IF

                     IF (BOTH) THEN

                        CALL SFEVAL (TYPE, WBUMP(N), XBUMP(N),
     >                               XWMIN(N), XWMAX(N), DYMAX, 1,NTOTK,
     >                               XINIT(1,K), YUNIT(1,K,N), IWRIT)
                     ELSE

                        IF (DOLO(N) /= ZERO) THEN

                           DL = DYMAX * SIGN (ONE, DOLO(N))

                           CALL SFEVAL (TYPE, WBUMP(N), XBUMP(N),
     >                                  XWMIN(N), XWMAX(N), DL, 1,NLOWK,
     >                                  XINIT(1,K), YUNIT(1,K,N), IWRIT)
                        END IF

                        IF (DOUP(N) /= ZERO) THEN

                           DU = DYMAX * SIGN (ONE, DOUP(N))

                           CALL SFEVAL (TYPE, WBUMP(N), XBUMP(N),
     >                                 XWMIN(N),XWMAX(N),DU,NLOWK,NTOTK,
     >                                 XINIT(1,K), YUNIT(1,K,N), IWRIT)
                        END IF

                     END IF

                  ELSE IF (N <= NCHCKS + NWVAR + NSINVC) THEN ! SIN*VC bump:
                                                              ! Varying Center
                     TYP4 = VTYPE(N)(1:4)
                     IF (TYP4(4:4) == 'E') TYP4(4:4) = ' '    ! See SFEVAL

C                    Determine the hinge X (or X of line parallel to shock)
C                    at this span station:

                     KI = KMIN(N)       ! Inboard end station of slat/flap
                     XHI = XLEADI(KI) + XWMIN(N) * CINIT(KI)  ! X at hinge,
                                                              ! inboard end
                     KO = KMAX(N)       ! Outboard ...
                     XHO = XLEADI(KO) + XWMAX(N) * CINIT(KO)

                     XHK = XHI + (Z - ZLEADI(KI)) * (XHO - XHI) /
     >                     (ZLEADI(KO) - ZLEADI(KI))          ! Hinge X
                     XCHK = (XHK - XLEADI(K)) / CINIT(K)      ! Corresp. X/C

                     IF (XBUMP(N) > ZERO) THEN

                        XCENTER = XBUMP(N) ! Refers to slat or bump chord

                        IF (XWMAX(N) < HALF) THEN ! Slat bump
                           XA = ZERO
                           XB = XCHK
                        ELSE ! XWMAX(N) >= HALF means flap bump
                           XA = XCHK
                           XB = ONE
                        END IF

                     ELSE ! XBUMP(N) = 0. means shock-following line
                         XCENTER = XCHK
                         XA = ZERO
                         XB = ONE
                     END IF

                     IF (DOLO(N) /= ZERO) THEN

                        DL = DYMAX * SIGN (ONE, DOLO(N))

                        CALL SFEVAL (TYP4, WBUMP(N), XCENTER,
     >                               XA, XB, DL, 1, NLOWK,
     >                               XINIT(1,K), YUNIT(1,K,N), IWRIT)
                     END IF

                     IF (DOUP(N) /= ZERO) THEN

                        DU = DYMAX * SIGN (ONE, DOUP(N))

                        CALL SFEVAL (TYP4, WBUMP(N), XCENTER,
     >                               XA, XB, DU, NLOWK, NTOTK,
     >                               XINIT(1,K), YUNIT(1,K,N), IWRIT)
                     END IF

                  ELSE ! Flap or slat rotation or shear (may not be used)

                     YUNIT(1,K,N) = ONE ! For this K, else default to 0.

                  END IF

               END IF

            END DO ! Next wing design variable

         END DO ! Next wing section


         IF (.NOT. WINGONLY) THEN

C           Store unit perturbations for certain body design variables:

            XBLONG  = XF(NFSTN) - XF(1)
            RXBLONG = ONE / XBLONG

            DO N = NTWING + 1, NTWING + NCMBRF ! Body camber variables
               TYPE = VTYPE(N)  ! 'BCxxxx'
               TYP4 = TYPE(3:6) ! 'xxxx'; SFEVAL expects NAME*4
               RXBRANGE = ONE / (XWMAX(N) - XWMIN(N)) ! XWMIN/MAX are in [0,1]
               DO I = 1, NFSTN
                  XNORM = ((XF(I) - XF(1))*RXBLONG - XWMIN(N))*RXBRANGE
                  IF (XNORM < ZERO) CYCLE ! Next station
                  IF (XNORM > ONE)  EXIT  ! No more stations to do

                  DELTAY = ZERO ! SFEVAL adds

                  CALL SFEVAL (TYP4, WBUMP(N), XBUMP(N), ZERO, ONE,
     >                         XBLONG, 1, 1, XNORM, DELTAY, IWRIT)

                  YUNIT(I,1:NF(I),N) = DELTAY

               END DO ! Next body section
            END DO ! Next body camber design variable


            IF (NSPANF > 0) THEN ! Body span variables

C              Determine unit vector components at the initial body points:

               DO I = 2, NFSTN ! Nose point has zero span/radius - avoid it

C                 Use SFUS(*) for what had been UFINIT(*,I) - no need
C                 for XFINIT, UFINIT, VFINIT for just NCALL = -4 usage.

                  NFP = NF(I)

                  CALL CHORDS2D (NFP, YFINIT(1,I), ZFINIT(1,I), .TRUE.,
     >                           STOTAL, SFUS)

C                 Originally, unit normals were chosen, but treatment of
C                 radial perturbations with linear constraints requires
C                 a common origin for the radial vectors, namely the
C                 body camber line:

C****             CALL UNITNORM ('B', NFP, ZFINIT(1,I), YFINIT(1,I),
C****>                           SFUS, YFUDGE, ZNORM(1,I), YNORM(1,I))

                  YCAMBER = (YFINIT(1,I) + YFINIT(NFP,I)) * HALF

                  YNORM(1,I) = -ONE ! Crown/keel pts. don't leave symmetry plane
                  ZNORM(1,I) = ZERO

                  DO J = 2, NFP - 1
                     ANGLE = ATAN2 (YFINIT(J,I) - YCAMBER, ZFINIT(J,I))
                     YNORM(J,I) = SIN (ANGLE)
                     ZNORM(J,I) = COS (ANGLE)
                  END DO

                  YNORM(NFP,I) = ONE
                  ZNORM(NFP,I) = ZERO
               END DO
            END IF

            DO N = NTWING + NCMBRF + 1, NTWING + NCMBRF + NSPANF
               TYPE2 = UTYPE(N)  ! Circumferential = U in (U,V) parameterization
               TYPE  = VTYPE(N)  ! 'BSxxxx', 'BVxxx', or 'BRxxx'; V = axial
               TYP1  = TYPE(2:2) ! 'S', 'V', or 'R'
               TYP4  = TYPE(3:6) ! 'xxxx'; SFEVAL expects NAME*4
               WIDTH = PBUMP(N)
               CENTR = UBUMP(N)
               UMIN  = UBMIN(N)
               UMAX  = UBMAX(N)
               RXBRANGE = ONE / (XWMAX(N) - XWMIN(N))

               DO I = 2, NFSTN

C                 Axial direction:

                  XDD = ((XF(I) - XF(1))*RXBLONG - XWMIN(N)) * RXBRANGE
                  IF (XDD < ZERO) CYCLE
                  IF (XDD > ONE ) EXIT

                  DZMAX = ZERO ! Max. dY, dZ, or dR at this X; SFEVAL adds

                  CALL SFEVAL (TYP4, WBUMP(N), XBUMP(N), ZERO, ONE,
     >                         XBLONG, 1, 1, XDD, DZMAX, IWRIT)

C                 Circumferential direction:

                  NFP = NF(I)
                  RBUMP(1:NFP) = ZERO ! SFEVAL adds

                  CALL SFEVAL (TYPE2, WIDTH, CENTR, UMIN, UMAX, DZMAX,
     >                         1, NFP, SFUS, RBUMP, IWRIT)

                  YUNIT(I,1:NFP,N) = RBUMP(1:NFP)

               END DO ! Next body section

            END DO ! Next body span design variable

         END IF ! Not-wing-only block

         GO TO 999

      END IF ! NCALL = -4 case


      INITIZE = NCALL == -3 .AND. FIRST

C     Initialize airfoil geometry for (optional) perturbations:

      DO K = 1, NWSEC
         AWORK(K)  = AINIT(K)
         TWORK(K)  = TINIT(K)
         CWORK(K)  = CINIT(K)
         XLEADW(K) = XLEADI(K)
         YLEADW(K) = YLEADI(K)
         ZLEADW(K) = ZLEADI(K)
         NUPPER(K) = NUPPER0(K) ! Flaps or slats may change these
         NLOWER(K) = NLOWER0(K)
         NTOTAL(K) = NTOTAL0(K)
      END DO

      DO K = 1, NWSEC
         DO I = 1, NTOTAL(K)
            XPRTB(I,K) = XINIT(I,K)
            YPRTB(I,K) = YINIT(I,K)
            ZPRTB(I,K) = ZINIT(I,K)
         END DO
      END DO

      IF (NDV == 0) GO TO 200


      IF (INITIZE) THEN
         FIRST = .FALSE.

         WRITE (IWRIT,1050) ' ',
     >                      '----------------------------',
     >                      '   OPTIMIZATION VARIABLES',
     >                      '----------------------------'
         WRITE (IWRIT,1100) '# WING TWIST      :', NTWIST,
     >                      '# WING THICKNESS  :', NTHICK,
     >                      '# WING X LEADING  :', NXLEAD,
     >                      '# WING Y LEADING  :', NYLEAD,
     >                      '# WING Z LEADING  :', NZLEAD,
     >                      '# CHORD           :', NCHORD,
     >                      '----------------------------'
         WRITE (IWRIT,1100) 'TOTAL # PLANFORM  :', NCHCKS,
     >                      '----------------------------'
         WRITE (IWRIT,1100) '# WING Y  SINE    :', NSINE,
     >                      '# WING Y  COS     :', NCOS,
     >                      '# WING Y  EXP     :', NEXPS,
     >                      '# WING Y  TRAILING:', NTRAIL,
     >                      '# WING Y  LEADING :', NLEAD,
     >                      '# WING Y  WAGNER  :', NWAG,
     >                      '----------------------------'
         WRITE (IWRIT,1100) 'TOTAL # WING Ys   :', NCHCK1,
     >                      '----------------------------'
         WRITE (IWRIT,1100) '# WING Y" SINE    :', NSINEP,
     >                      '# WING Y" EXP     :', NEXPSP,
     >                      '# WING Y" TRAILING:', NTRAIP,
     >                      '# WING Y" LEADING :', NLEADP,
     >                      '----------------------------'
         WRITE (IWRIT,1100) 'TOTAL # WING Y"s  :', NCHCK2,
     >                      '----------------------------'
         WRITE (IWRIT,1100) 'TOTAL VARIABLE CTR:', NSINVC,
     >                      '----------------------------'
         WRITE (IWRIT,1100) '# WING Y  FLAP    :', NFLAP,
     >                      '# WING Y  SLAT    :', NSLAT,
     >                      '----------------------------'
         WRITE (IWRIT,1100) 'TOTAL # FLAP/SLAT :', NTOTFS,
     >                      '----------------------------'
         WRITE (IWRIT,1100) 'TOTAL # WING      :', NTWING,
     >                      '----------------------------'
         WRITE (IWRIT,1100) '# BODY CAMBER SINE:', NSINFC,
     >                      '# BODY CAMBER COS :', NCOSFC,
     >                      '# BODY CAMBER EXP :', NEXPFC,
     >                      '# BODY CAMBER LEAD:', NLEDFC,
     >                      '# BODY CAMBER TRL :', NTRLFC,
     >                      '# BODY SPAN   SINE:', NSINFS,
     >                      '# BODY SPAN   COS :', NCOSFS,
     >                      '# BODY SPAN   EXP :', NEXPFS,
     >                      '# BODY SPAN   LEAD:', NLEDFS,
     >                      '# BODY SPAN   TRL :', NTRLFS,
     >                      '# BODY AREA   SINE:', NSINFA,
     >                      '# BODY AREA   COS :', NCOSFA,
     >                      '# BODY AREA   EXP :', NEXPFA,
     >                      '# BODY AREA   LEAD:', NLEDFA,
     >                      '# BODY AREA   TRL :', NTRLFA,
     >                      '----------------------------'
         WRITE (IWRIT,1100) 'TOTAL # BODY      :', NTOTFU,
     >                      '----------------------------'
         WRITE (IWRIT,1100) '# FLOOR KINK  X   :', NXFLR,
     >                      '# FLOOR KINK  Y   :', NYFLR,
     >                      '----------------------------'
         WRITE (IWRIT,1100) 'TOTAL # FLOOR KINK:', NTKINK,
     >                      '----------------------------'
         WRITE (IWRIT,1100) 'TOTAL # VARIABLES :', NDV,
     >                      '----------------------------'
      END IF


C     Perturb the wing?

      DO 100 K = 1, NWSEC

         NTOTK = NTOTAL(K)
         NUPPK = NUPPER(K)
         NLOWK = NLOWER(K)

         IF (ABS (XPRTB(NLOWK,K)) > 0.0000001) THEN
            WRITE (IWRIT,*) 'Leading edge not at X = 0.  K: ', K
            WRITE (IWRIT,*) 'Leading edge X: ', XPRTB(NLOWK,K)
C***        GO TO 990
         ELSE IF (ABS (MAX (XPRTB(1,K), XPRTB(NTOTK,K)) - ONE)
     >          >  0.0000001) THEN
            WRITE (IWRIT,*) 'Trailing edge not at X = 1. K: ', K
            WRITE (IWRIT,*) 'Trailing edge Xs: ',
     >         XPRTB(1,K), XPRTB(NTOTK,K)
C***        GO TO 990
         END IF

C        Planform variables:
C        -------------------

         DO N = 1, NCHCKS

C           Determine this design variable's contribution at this
C           span station:   0. <= YUNIT(1,K,N) <= 1.

            DELTAV = YUNIT(1,K,N) * V(N) * VSCALE(N)

            IF (DELTAV /= ZERO) THEN
               TYPE = VTYPE(N)

               IF      (TYPE == 'TWT   ') THEN ! Twist about XBUMP in [0, 1]
                  ALPHA1    = AWORK(K)*DEG2RAD
                  AWORK(K)  = AWORK(K)  + DELTAV
                  ALPHA2    = AWORK(K)*DEG2RAD
                  XLEADW(K) = XLEADW(K) + XBUMP(N)*CWORK(K)*
     >                        (COS (ALPHA1) - COS (ALPHA2))
                  YLEADW(K) = YLEADW(K) + XBUMP(N)*CWORK(K)*
     >                        (SIN (ALPHA2) - SIN (ALPHA1))
               ELSE IF (TYPE == 'THK   ') THEN
                  TWORK(K)  = TWORK(K)  + DELTAV
               ELSE IF (TYPE == 'XLD   ') THEN
                  XLEADW(K) = XLEADW(K) + DELTAV
               ELSE IF (TYPE == 'YLD   ') THEN
                  YLEADW(K) = YLEADW(K) + DELTAV
               ELSE IF (TYPE == 'ZLD   ') THEN
                  ZLEADW(K) = ZLEADW(K) + DELTAV
               ELSE IF (TYPE == 'CRD   ') THEN
                  CWORK(K)  = CWORK(K)  + DELTAV
               ELSE
                  GO TO 800
               END IF
            END IF
         END DO


C        Wing section Y perturbations (including flap/slat bumps)
C        --------------------------------------------------------

C        Apply a multiple of the stored unit perturbations:

         DO N = NCHCKS + 1, NCHCKS + NCHCK1 + NSINVC
            DY = V(N) * VSCALE(N)

            IF (DY /= ZERO) THEN
               DO I = 1, NTOTK
                  YPRTB(I,K) = YPRTB(I,K) + DY * YUNIT(I,K,N)
               END DO
            END IF
         END DO


C        Perturbations applied to section Y"s?
C        -------------------------------------

         IF (NCHCK2 == 0) GO TO 50

         DO I = 2, NTOTK - 1
            A(I)  =  ONE / (XPRTB(I,K) - XPRTB(I-1,K))
            B(I)  = -ONE / (XPRTB(I+1,K) - XPRTB(I,K))
     >              -ONE / (XPRTB(I,K) - XPRTB(I-1,K))
            CC(I) =  ONE / (XPRTB(I+1,K) - XPRTB(I,K))
            DD(I) = (YPRTB(I+1,K) - YPRTB(I,K)) * CC(I) -
     >              (YPRTB(I,K) - YPRTB(I-1,K)) *  A(I)
         END DO

         DO N = NCHCKS + 1, NCHCKS + NCHCK2 ! NCHCK1 is 0
            DY = V(N) * VSCALE(N)

            IF (DY /= ZERO) THEN
               DO I = 2, NTOTK - 1
                  DD(I) = DD(I) + DY * YUNIT(I,K,N)
               END DO
            END IF
         END DO ! Next design variable

         A(1)  = ZERO
         B(1)  = ONE
         CC(1) = ZERO
         DD(1) = YPRTB(1,K)
         A(NLOWK)  = ZERO
         B(NLOWK)  = ONE
         CC(NLOWK) = ZERO
         DD(NLOWK) = YPRTB(NLOWK,K)
         A(NTOTK)  = ZERO
         B(NTOTK)  = ONE
         CC(NTOTK) = ZERO
         DD(NTOTK) = YPRTB(NTOTK,K)

C        Recover the Ys from the Y"s by solving a tridiagonal system.

         CALL TRDIAG (A, B, CC, DD, DD, NTOTK)

         N = NCHCKS + 1 ! Assume same upper/lower control for all Y" variables

         IF (ABS (DOUP(N)) >= 3. .OR. ABS (DOLO(N)) >= 3.) THEN
            YPRTB(1:NTOTK,K) = DD(1:NTOTK)
         ELSE
            DO I = NLOWK, NTOTK
               DU = YPRTB(I,K)
               YPRTB(I,K) = DD(I)
               J = 2*NLOWK - I
               YPRTB(J,K) = YPRTB(J,K) + (DD(I) - DU)
            END DO
         END IF

   50    CONTINUE

         IF (IYSMOO /= 0) THEN ! Smooth the airfoil Ys in place.

            CALL YSMOOP (NTOTK, XPRTB(1,K), YPRTB(1,K), EPSY, EXPY,
     >                   A, B, CC, DD)
         END IF


C        (Untwisted) section maximum thickness:

         CALL GETTHICK (NLOWK, NUPPK, XPRTB(1,K), YPRTB(1,K), A,
     >                  THMAX(K), XTHMAX(K))

C        THICK scaling is done below.


C        Apply section twists and denormalize:
C        -------------------------------------

         ALPHA     = AWORK(K) + ALPHAW
         ALG(K)    = ALPHA
         THICK     = TWORK(K)
         THMAX(K)  = THMAX(K) * THICK
         SCALE     = CWORK(K)
         CHORDG(K) = SCALE
         NW(K)     = NTOTK
         NLG(K)    = NLOWK

         XXW = XPRTB(NLOWK,K)
         YYW = YPRTB(NLOWK,K)
         XLW = XLEADW(K) + XXW * SCALE
         YLW = YLEADW(K) + YYW * SCALE * THICK
         ZLW = ZLEADW(K)
         CA  = COS (ALPHA*DEG2RAD)
         SA  = SIN (ALPHA*DEG2RAD)

         DO I = 1, NTOTK
            XWG(I,K) = XLW + SCALE * ((XPRTB(I,K) - XXW)*CA +
     >                                (YPRTB(I,K) - YYW)*SA * THICK)
            YWG(I,K) = YLW + SCALE * ((YPRTB(I,K) - YYW)*CA * THICK -
     >                                (XPRTB(I,K) - XXW)*SA)
            ZWG(I,K) = ZLW + SCALE *   ZPRTB(I,K)
         END DO

         FAREA(K) = ZERO
         DO I = 2, NLOWK
            J = NTOTK + 1 - I
            DXDX1 = XWG(J,K) - XWG(I-1,K)
            DXDX2 = XWG(I,K) - XWG(J+1,K)
            DYDY1 = YWG(J,K) - YWG(I-1,K)
            DYDY2 = YWG(I,K) - YWG(J+1,K)
            DAREA = HALF * ABS (DXDX1*DYDY2 - DXDX2*DYDY1)
            FAREA(K) = FAREA(K) + DAREA
         END DO

  100 CONTINUE ! Next wing section


      WINGVL = ZERO
      DO K = 2, NWSEC
         WINGVL = (ZWG(1,K) - ZWG(1,K-1))*(FAREA(K) + FAREA(K-1)) +
     >            WINGVL
      END DO
      WINGVL = WINGVL * HALF ! Before any flap/slat deflections


C     Flap/slat rotations in REAL space:
C     ----------------------------------

      IF (NTOTFS > 0) CALL DEFLECT ()


  200 IF (WINGONLY) GO TO 999


C-----------------------------------------------------------------------
C                        Fuselage Modifications
C-----------------------------------------------------------------------

      IF (INITIZE) THEN
         FIRST = .FALSE.
      END IF

      IF (NNFLOOR /= 0) THEN

C        Set up cabin floor break points for the case of no body camber changes:

         XFLOORF(1:NFLOOR) = XFLOORI(1:NFLOOR)
         YFLOORF(1:NFLOOR) = YFLOORI(1:NFLOOR)
      END IF

C     Initialize the current body geometry:

      DO I = 1, NFSTN
         DO J = 1, NF(I)
            YFPRTB(J,I) = YFINIT(J,I)
            ZFPRTB(J,I) = ZFINIT(J,I)
         END DO
      END DO

      BODYFIXD = NCMBRF == 0 .AND. NSPANF == 0 .AND. NAREAF == 0

      IF (BODYFIXD) GO TO 500


C-----------------------------------------------------------------------
C     Fuselage camber modifications.
C-----------------------------------------------------------------------

      XBLONG  = XF(NFSTN) - XF(1)
      RXBLONG = ONE / XBLONG

      DO N = NTWING + 1, NTWING + NCMBRF
         DELTAY = V(N) * VSCALE(N)
         DO I = 1, NFSTN
            DO J = 1, NF(I)
               YFPRTB(J,I) = YFPRTB(J,I) + DELTAY * YUNIT(I,J,N)
            END DO
         END DO
      END DO

C     Check for consequences of body cambering:

      IF (NCMBRF > 0) THEN

         IF (NNFLOOR == 2) THEN

C           Apply body camber to the cabin floor break points:

            DO I = 1, NFLOOR
               DELTAY = ZERO

               DO N = NTWING + 1, NTWING + NCMBRF
                  RXBRANGE = ONE / (XWMAX(N) - XWMIN(N))
                  XDD = ((XFLOORF(I)-XF(1))*RXBLONG-XWMIN(N))*RXBRANGE
                  IF (XDD < ZERO) CYCLE
                  IF (XDD > ONE)  CYCLE

                  TYPE  = VTYPE(N)  ! 'BCxxxx'
                  TYP4  = TYPE(3:6) ! 'xxxx'; SFEVAL expects NAME*4
                  YMULT = V(N) * VSCALE(N)

                  CALL SFEVAL (TYP4, WBUMP(N), XBUMP(N), ZERO, ONE,
     >                         YMULT, 1, 1, XDD, DELTAY, IWRIT)
               END DO

               YFLOORF(I) = YFLOORF(I) + DELTAY * XBLONG

            END DO ! Next floor break point

         END IF

      END IF

C-----------------------------------------------------------------------
C     Fuselage section modifications (spanwise, vertical, or radial).
C-----------------------------------------------------------------------

      DO N = NTWING + NCMBRF + 1, NTWING + NCMBRF + NSPANF
         DELTAY = V(N) * VSCALE(N)
         TYP1   = VTYPE(N)(2:2)

         IF (TYP1 == 'S') THEN ! Spanwise
            DO I = 1, NFSTN
               DO J = 1, NF(I)
                  ZFPRTB(J,I) = ZFPRTB(J,I) + DELTAY * YUNIT(I,J,N)
               END DO
            END DO
         ELSE IF (TYP1 == 'V') THEN ! Vertical
            DO I = 1, NFSTN
               DO J = 1, NF(I)
                  YFPRTB(J,I) = YFPRTB(J,I) + DELTAY * YUNIT(I,J,N)
               END DO
            END DO
         ELSE ! TYP1 = 'R' (Radial)
            DO I = 1, NFSTN
               DO J = 1, NF(I)
                  ZFPRTB(J,I) = (DELTAY * YUNIT(I,J,N)) * ZNORM(J,I) +
     >                          ZFPRTB(J,I)
                  YFPRTB(J,I) = (DELTAY * YUNIT(I,J,N)) * YNORM(J,I) +
     >                          YFPRTB(J,I)
               END DO
            END DO
         END IF
      END DO

C-----------------------------------------------------------------------
C     Fuselage area modifications (sections scaled/same shape).
C-----------------------------------------------------------------------

      DO N = NTWING + NCMBRF + NSPANF + 1,
     >           NTWING + NCMBRF + NSPANF + NAREAF

         TYPE  = VTYPE(N)  ! 'BAxxxx'
         TYP4  = TYPE(3:6) ! 'xxxx'; SFEVAL expects NAME*4
         WIDTH = WBUMP(N)
         CENTR = XBUMP(N)
         SLOPE = PBUMP(N)  ! Calculated in RDVARS (not an input)
         AMULT = TWO * V(N) * VSCALE(N) * XBLONG ** 2
         RXBRANGE = ONE / (XWMAX(N) - XWMIN(N))

         DO I = 2, NFSTN - 1

            XDD = ((XF(I) - XF(1)) * RXBLONG - XWMIN(N)) * RXBRANGE
            IF (XDD < ZERO) CYCLE  ! Next body station
            IF (XDD > ONE)  EXIT   ! Next body design variable

C           Calculate the unperturbed area, and set up angular coordinates
C           (or rather sines and cosines of them) for the section points:

            NFP   = NF(I)
            YKEEL = YFPRTB(1,I)
            YCAM  = (YKEEL + YFPRTB(NFP,I)) * HALF
            CCX   = YCAM - YKEEL
            TOTA  = ZERO
            SINTHETA(1) = ZERO
            COSTHETA(1) = ONE

            DO J = 2, NFP - 1
               COF1 = YFPRTB(J,I)   - YCAM
               COF2 = ZFPRTB(J,I)
               COF3 = YFPRTB(J-1,I) - YCAM
               COF4 = ZFPRTB(J-1,I)
               TOTA = TOTA + (COF1*COF4 - COF3*COF2) ! Full section, not half

               AAX = SQRT ((YFPRTB(J,I) - YKEEL)**2 + ZFPRTB(J,I)**2)
               BBX = SQRT (                 COF1**2 + ZFPRTB(J,I)**2)

               COSTHETA(J) = (AAX**2 - BBX**2 - CCX**2) / (-2.*BBX*CCX)
               SINTHETA(J) = SQRT (MAX (ZERO, ONE - COSTHETA(J)**2))
            END DO
            SINTHETA(NFP) = ZERO
            COSTHETA(NFP) = -ONE

            TOTR = SQRT (TOTA / PI) ! Average radius of unperturbed section

C           Determine desired area:

            IF (TYP4 == 'SV1 ') THEN
               IF (XDD <= HALF) THEN
                  BUMP1 = VBUMP1 (WIDTH, CENTR, SLOPE, XDD)
               ELSE
                  BUMP1 = VBUMP2 (WIDTH, CENTR, SLOPE, XDD)
               END IF
            ELSE IF (TYP4 == 'SV2 ') THEN
               IF (XDD <= HALF) THEN
                  BUMP1 = VBUMP2 (WIDTH, CENTR, SLOPE, XDD)
               ELSE
                  BUMP1 = VBUMP1 (WIDTH, CENTR, SLOPE, XDD)
               END IF
            ELSE
               BUMP1 = ZERO ! SFEVAL adds

               CALL SFEVAL (TYP4, WIDTH, CENTR, ZERO, ONE, ONE,
     >                      1, 1, XDD, BUMP1, IWRIT)
            END IF

            TAGA = AMULT * BUMP1 + TOTA ! Desired area
            TAGR = SQRT (TAGA / PI)     ! Initial estimate of corresp. radius
            DIFR = TAGR - TOTR
            TAGB = TAGA

C           Iterate on averaged radius until the desired area is achieved:

            AREATOL = MAX (100.* EPSMCH, TOLER)

CCC            write (6,'(/,a,i4,2f14.6,1p,2e12.3)')
CCC     >         ' i, tota, taga, amult, bump1: ',
CCC     >           i, tota, taga, amult, bump1

            DO ITER = 1, 30

               DO J = 1, NFP
                  YFPRTB(J,I) = YFPRTB(J,I) - DIFR * COSTHETA(J)
                  ZFPRTB(J,I) = ZFPRTB(J,I) + DIFR * SINTHETA(J)
               END DO

               TOTB = ZERO
               DO J = 2, NFP
                  COF1 = YFPRTB(J,I) - YCAM
                  COF2 = ZFPRTB(J,I)
                  COF3 = YFPRTB(J-1,I) - YCAM
                  COF4 = ZFPRTB(J-1,I)
                  TOTB = TOTB + (COF1*COF4 - COF3*COF2)
               END DO

               DIFA = TOTB - TAGA

CCC               write (6,'(1x,a,i3,2f14.6,2f12.6)')
CCC     >            'iter, totb, taga, totr, tagr: ',
CCC     >             iter, totb, taga, totr, tagr

               IF (ABS (DIFA) < TOTA * AREATOL) THEN
                  EXIT
               ELSE
                  TOTR = TAGR
                  TAGB = TAGB - DIFA
                  TAGR = SQRT (TAGB / PI)
                  DIFR = TAGR - TOTR
               END IF

               IF (ITER == 30) THEN
                  WRITE (IWRIT,'(A,/,A,I5)')
     >          ' No convergence on target fuselage area distribution.',
     >          ' Aborting at fuselage station', I
                  GO TO 990
               END IF

            END DO ! Next iteration

         END DO ! Next body station

      END DO ! Next body area design variable


C     Cabin floor kink point locations:

      DO N = NTWING + NCMBRF + NSPANF + NAREAF + 1,
     >       NTWING + NCMBRF + NSPANF + NAREAF + NTKINK

         TYPE = VTYPE(N)
         I = KCEN(N) ! Kink point #

         IF (TYPE == 'XFLR ') THEN
            XFLOORF(I) = XFLOORI(I) + V(N) * VSCALE(N)
         ELSE IF (TYPE == 'YFLR ') THEN
            YFLOORF(I) = YFLOORI(I) + V(N) * VSCALE(N)
         END IF

      END DO

  500 CONTINUE

  600 CONTINUE

C     --------------------------------------------------------------------
C     Option to fudge the body sections if the wing drops below the body:
C     --------------------------------------------------------------------

      DO I = 1, NFSTN
         NFP = NF(I)
         YFUDGE(1:NFP,I) = YFPRTB(1:NFP,I)
         ZFUDGE(1:NFP,I) = ZFPRTB(1:NFP,I)
      END DO

      NEWBODY = .FALSE. ! WBFUDGE indicates if it fudged the body

      IF (NNROOT == 3) THEN

         CALL WBFUDGE (NCALL, IBOD1, IBOD2, MXJBODY, MXIBODY, NF,
     >                 XF, YFUDGE, ZFUDGE, KWING2,
     >                 MXIWING, MXKWING, XWG, YWG, ZWG, NW, NLG,
     >                 YSHELF, ZSHELF, IWRIT, NEWBODY)
      END IF

C     --------------------------------------------------------------------
C     Interpolate the possibly-perturbed fuselage stns. to a regular grid?
C     --------------------------------------------------------------------

      IF (.NOT. INITIZE) THEN
         IF (BODYFIXD) THEN
            IF (.NOT. NEWBODY) GO TO 999
         END IF
      END IF

      IF (IDEGBOD == 3) THEN
         IF (NNROOT < 3) THEN
            METHOD = 'B' ! Loose local cubic spline fits
         ELSE
            METHOD = 'M' ! Monotonic ...
         END IF
      ELSE
         METHOD = 'L'    ! Piecewise linear
      END IF

C     Singular point at the nose:

      YF(1:NJBODY,1) = YFUDGE(1,1)
      ZF(1:NJBODY,1) = ZFUDGE(1,1)

      SFUS(1)  = ZERO
      SUNIF(1) = ZERO

      DO I = 2, NFSTN

         NFP = NF(I)

         CALL CHORDS2D (NFP, YFUDGE(1,I), ZFUDGE(1,I),.FALSE., DS, SFUS)

         SUNIF(NJBODY) = DS
         DS = DS / REAL (NJBODY-1)

         DO J = 2, NJBODY - 1
            SUNIF(J) = DS * REAL (J-1)
         END DO

         CALL LCSFIT (NFP, SFUS, YFUDGE(1,I), NEW, METHOD, NJBODY,
     >                SUNIF, YF(1,I), YF(1,I))

         CALL LCSFIT (NFP, SFUS, ZFUDGE(1,I), NEW, METHOD, NJBODY,
     >                SUNIF, ZF(1,I), ZF(1,I))
      END DO

      IF (NCALL == -3 .AND. .NOT. CHECKWARP) THEN
         WRITE (LUNBOD) 1
         WRITE (LUNBOD) NJBODY, NFSTN, 1
         WRITE (LUNBOD) ((XF(I),   J = 1, NJBODY), I = 1, NFSTN),
     >                  ((YF(J,I), J = 1, NJBODY), I = 1, NFSTN),
     >                  ((ZF(J,I), J = 1, NJBODY), I = 1, NFSTN)
         CLOSE (LUNBOD)
      END IF

      GO TO 999


  800 WRITE (IWRIT, 1010) TYPE

  990 STOP

  999 RETURN

C     Formats:

 1000 FORMAT (A)
 1010 FORMAT (/, ' PERTURB: Bad design variable order or type: ', A)
 1050 FORMAT (1X, A)
 1100 FORMAT (1X, A, I9)

      END SUBROUTINE PERTURB

C***********************************************************************
C
      SUBROUTINE RDVARS (LUNRD)
C
C     File-driven control scheme for wing/body design variables.
C     See PERTURB for the descriptions.  Variable order is important.
C
C***********************************************************************

C     Global quantities:

      USE LUNS
      USE OPT

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) :: LUNRD

C     Local constants:

      REAL, PARAMETER :: ONE = 1., TWO = 2.

C     Local variables:

      INTEGER    I, N, IER
      LOGICAL    ORDERED
      CHARACTER  TYPE * 6, TYPE2 * 4

C     Execution:

      NTWIST = 0
      NTHICK = 0
      NXLEAD = 0
      NYLEAD = 0
      NZLEAD = 0
      NCHORD = 0
      NSINE  = 0
      NCOS   = 0
      NEXPS  = 0
      NTRAIL = 0
      NLEAD  = 0
      NWAG   = 0
      NSINVC = 0
      NSINEP = 0
      NEXPSP = 0
      NTRAIP = 0
      NLEADP = 0
      NFLAP  = 0
      NSLAT  = 0
      NSINFC = 0
      NCOSFC = 0
      NEXPFC = 0
      NLEDFC = 0
      NTRLFC = 0
      NSINFS = 0
      NCOSFS = 0
      NEXPFS = 0
      NLEDFS = 0
      NTRLFS = 0
      NSINFA = 0
      NCOSFA = 0
      NEXPFA = 0
      NLEDFA = 0
      NTRLFA = 0
      NXFLR  = 0
      NYFLR  = 0

      READ (LUNRD,*)
      READ (LUNRD,*)

      DO I = 1, NDV

         CALL RDVLINE (I, VTYPE, UTYPE, DOUP, DOLO,
     >                 WBUMP, XBUMP, XWMIN, XWMAX,
     >                 PBUMP, UBUMP, UBMIN, UBMAX,
     >                 KCEN, KMIN, KMAX, V, VSCALE, AITCH, BL, BU,
     >                 NBOUND, LUNRD, IWRIT)

         TYPE  = VTYPE(I)
         TYPE2 = UTYPE(I)

         IF      (TYPE2 == 'POLY') THEN ! Just check for bad inputs
         ELSE IF (TYPE2 == 'SIN ') THEN
         ELSE IF (TYPE2 == 'SINF') THEN
         ELSE IF (TYPE2 == 'SIN1') THEN
         ELSE IF (TYPE2 == 'SIN2') THEN
         ELSE IF (TYPE2 == 'SIN3') THEN
         ELSE IF (TYPE2 == 'SIN4') THEN
         ELSE IF (TYPE2 == 'COSL') THEN
         ELSE IF (TYPE2 == 'COSR') THEN
         ELSE IF (TYPE2 == 'LCOS') THEN
         ELSE IF (TYPE2 == 'RCOS') THEN
         ELSE IF (TYPE2 == 'LED ') THEN
         ELSE IF (TYPE2 == 'TRL ') THEN
         ELSE IF (TYPE2 == 'EXP ') THEN
         ELSE
            GO TO 800 ! Unknown transverse lofting variable name
         END IF

         IF      (TYPE == 'TWT   ') THEN
            NTWIST = NTWIST + 1
            IGROUP(I) = 1
         ELSE IF (TYPE == 'THK   ') THEN
            NTHICK = NTHICK + 1
            IGROUP(I) = 1
         ELSE IF (TYPE == 'XLD   ') THEN
            NXLEAD = NXLEAD + 1
            IGROUP(I) = 1
         ELSE IF (TYPE == 'YLD   ') THEN
            NYLEAD = NYLEAD + 1
            IGROUP(I) = 1
         ELSE IF (TYPE == 'ZLD   ') THEN
            NZLEAD = NZLEAD + 1
            IGROUP(I) = 1
         ELSE IF (TYPE == 'CRD   ') THEN
            NCHORD = NCHORD + 1
            IGROUP(I) = 1
         ELSE IF (TYPE == 'SIN   ' .OR. TYPE == 'SINE  ') THEN
            NSINE  = NSINE  + 1
            TYPE = 'SIN   ' ! Make certain for PERTURB
            IGROUP(I) = 2
         ELSE IF (TYPE == 'SINF  ') THEN
            NSINE  = NSINE  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'SIN1  ') THEN
            NSINE  = NSINE  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'SIN2  ') THEN
            NSINE  = NSINE  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'SIN3  ') THEN
            NSINE  = NSINE  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'SIN4  ') THEN
            NSINE  = NSINE  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'COSL  ') THEN
            NCOS  = NCOS  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'COSR  ') THEN
            NCOS  = NCOS  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'LCOS  ') THEN
            NCOS  = NCOS  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'RCOS  ') THEN
            NCOS  = NCOS  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'EXP   ') THEN
            NEXPS  = NEXPS  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'TRL   ') THEN
            NTRAIL = NTRAIL + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'LED   ') THEN
            NLEAD  = NLEAD  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE(1:3) == 'WAG') THEN
            TYPE = 'WAG   ' ! Make certain for PERTURB
            NWAG = NWAG  + 1
            IGROUP(I) = 2
         ELSE IF (TYPE == 'DSIN  ') THEN
            NSINEP = NSINEP + 1
            IGROUP(I) = 3
         ELSE IF (TYPE == 'DEXP  ') THEN
            NEXPSP = NEXPSP + 1
            IGROUP(I) = 3
         ELSE IF (TYPE == 'DTRL  ') THEN
            NTRAIP = NTRAIP + 1
            IGROUP(I) = 3
         ELSE IF (TYPE == 'DLED  ') THEN
            NLEADP = NLEADP + 1
            IGROUP(I) = 3
         ELSE IF (TYPE(5:6) == 'VC') THEN ! SIN*VC (variable center)
            NSINVC = NSINVC  + 1
            IGROUP(I) = 4
         ELSE IF (TYPE == 'FLAP  ') THEN
            NFLAP  = NFLAP  + 1
            IGROUP(I) = 5
         ELSE IF (TYPE == 'FLAPSH') THEN
            NFLAP  = NFLAP  + 1
            IGROUP(I) = 5
         ELSE IF (TYPE == 'SLAT  ') THEN
            NSLAT  = NSLAT  + 1
            IGROUP(I) = 5
         ELSE IF (TYPE == 'SLATSH') THEN
            NSLAT  = NSLAT  + 1
            IGROUP(I) = 5
         ELSE IF (TYPE == 'BCSIN ') THEN
            NSINFC = NSINFC + 1
            IGROUP(I) = 6
         ELSE IF (TYPE == 'BCSIN1') THEN
            NSINFC = NSINFC + 1
            IGROUP(I) = 6
         ELSE IF (TYPE == 'BCSIN2') THEN
            NSINFC = NSINFC + 1
            IGROUP(I) = 6
         ELSE IF (TYPE == 'BCSINF') THEN
            NSINFC = NSINFC + 1
            IGROUP(I) = 6
         ELSE IF (TYPE == 'BCCOS ') THEN
            NCOSFC = NCOSFC + 1
            IGROUP(I) = 6
         ELSE IF (TYPE == 'BCEXP ') THEN
            NEXPFC = NEXPFC + 1
            IGROUP(I) = 6
         ELSE IF (TYPE == 'BCLED ') THEN
            NLEDFC = NLEDFC + 1
            IGROUP(I) = 6
         ELSE IF (TYPE == 'BCTRL ') THEN
            NTRLFC = NTRLFC + 1
            IGROUP(I) = 6
         ELSE IF (TYPE == 'BSSIN ') THEN
            NSINFS = NSINFS + 1
            IGROUP(I) = 7
         ELSE IF (TYPE == 'BSSIN1') THEN
            NSINFS = NSINFS + 1
            IGROUP(I) = 7
         ELSE IF (TYPE == 'BSSIN2') THEN
            NSINFS = NSINFS + 1
            IGROUP(I) = 7
         ELSE IF (TYPE == 'BSSINF') THEN
            NSINFS = NSINFS + 1
            IGROUP(I) = 7
         ELSE IF (TYPE == 'BSCOS ') THEN
            NCOSFS = NCOSFS + 1
            IGROUP(I) = 7
         ELSE IF (TYPE == 'BSEXP ') THEN
            NEXPFS = NEXPFS + 1
            IGROUP(I) = 7
         ELSE IF (TYPE == 'BSLED ') THEN
            NLEDFS = NLEDFS + 1
            IGROUP(I) = 7
         ELSE IF (TYPE == 'BSTRL ') THEN
            NTRLFS = NTRLFS + 1
            IGROUP(I) = 7
         ELSE IF (TYPE == 'BVSIN ') THEN
            NSINFS = NSINFS + 1
            IGROUP(I) = 7
         ELSE IF (TYPE == 'BVSIN1') THEN
            NSINFS = NSINFS + 1
            IGROUP(I) = 7
         ELSE IF (TYPE == 'BVSIN2') THEN
            NSINFS = NSINFS + 1
            IGROUP(I) = 7
         ELSE IF (TYPE == 'BVSINF') THEN
            NSINFS = NSINFS + 1
            IGROUP(I) = 7
         ELSE IF (TYPE == 'BVCOS ') THEN
            NCOSFS = NCOSFS + 1
            IGROUP(I) = 7
         ELSE IF (TYPE == 'BVEXP ') THEN
            NEXPFS = NEXPFS + 1
            IGROUP(I) = 7
         ELSE IF (TYPE == 'BVLED ') THEN
            NLEDFS = NLEDFS + 1
            IGROUP(I) = 7
         ELSE IF (TYPE == 'BVTRL ') THEN
            NTRLFS = NTRLFS + 1
            IGROUP(I) = 7
         ELSE IF (TYPE == 'BRSIN ') THEN
            NSINFS = NSINFS + 1
            IGROUP(I) = 7
         ELSE IF (TYPE == 'BRSIN1') THEN
            NSINFS = NSINFS + 1
            IGROUP(I) = 7
         ELSE IF (TYPE == 'BRSIN2') THEN
            NSINFS = NSINFS + 1
            IGROUP(I) = 7
         ELSE IF (TYPE == 'BRSINF') THEN
            NSINFS = NSINFS + 1
            IGROUP(I) = 7
         ELSE IF (TYPE == 'BRCOS ') THEN
            NCOSFS = NCOSFS + 1
            IGROUP(I) = 7
         ELSE IF (TYPE == 'BREXP ') THEN
            NEXPFS = NEXPFS + 1
            IGROUP(I) = 7
         ELSE IF (TYPE == 'BRLED ') THEN
            NLEDFS = NLEDFS + 1
            IGROUP(I) = 7
         ELSE IF (TYPE == 'BRTRL ') THEN
            NTRLFS = NTRLFS + 1
            IGROUP(I) = 7
         ELSE IF (TYPE == 'BASIN ') THEN
            NSINFA = NSINFA + 1
            IGROUP(I) = 8
         ELSE IF (TYPE == 'BASIN1') THEN
            NSINFA = NSINFA + 1
            IGROUP(I) = 8
         ELSE IF (TYPE == 'BASIN2') THEN
            NSINFA = NSINFA + 1
            IGROUP(I) = 8
         ELSE IF (TYPE == 'BASINF') THEN
            NSINFA = NSINFA + 1
            IGROUP(I) = 8
         ELSE IF (TYPE == 'BACOS ') THEN
            NCOSFA = NCOSFA + 1
            IGROUP(I) = 8
         ELSE IF (TYPE == 'BAEXP ') THEN
            NEXPFA = NEXPFA + 1
            IGROUP(I) = 8
         ELSE IF (TYPE == 'BALED ') THEN
            NLEDFA = NLEDFA + 1
            IGROUP(I) = 8
         ELSE IF (TYPE == 'BATRL ') THEN
            NTRLFA = NTRLFA + 1
            IGROUP(I) = 8
         ELSE IF (TYPE == 'BASV1 ') THEN
            NSINFA = NSINFA + 1
            IGROUP(I) = 8

C           Determine slope of line giving zero area for composite function:

            CALL ZEROAREA (XBUMP(I), WBUMP(I), PBUMP(I), ONE,-IWRIT,IER)

            IF (IER /= 0) THEN
               WRITE (IWRIT,*)
     >            'Error initializing design variable ', I, ': ', TYPE
               GO TO 900
            END IF

         ELSE IF (TYPE == 'BASV2 ') THEN
            NSINFA = NSINFA + 1
            IGROUP(I) = 8

            CALL ZEROAREA (XBUMP(I), WBUMP(I), PBUMP(I), TWO,-IWRIT,IER)

            IF (IER /= 0) THEN
               WRITE (IWRIT,*)
     >            'Error initializing design variable ', I, ': ', TYPE
               GO TO 900
            END IF

         ELSE IF (TYPE == 'XFLR  ') THEN ! Kink point # (INTEGER) is
            NXFLR = NXFLR + 1            ! read into KCEN(I) by RDVLINE
            IGROUP(I) = 9

         ELSE IF (TYPE == 'YFLR  ') THEN !  "   "
            NYFLR = NYFLR + 1
            IGROUP(I) = 9

         ELSE
            GO TO 800 ! Unknown design variable name
         END IF

      END DO

      NCHCKS = NTWIST + NTHICK + NXLEAD + NZLEAD + NYLEAD + NCHORD
      NCHCK1 = NSINE  + NCOS   + NEXPS  + NTRAIL + NLEAD  + NWAG
      NCHCK2 = NSINEP + NEXPSP + NTRAIP + NLEADP
      NTOTFS = NFLAP  + NSLAT
      NTWING = NCHCKS + NCHCK1 + NCHCK2 + NSINVC + NTOTFS
      NCMBRF = NSINFC + NCOSFC + NEXPFC + NLEDFC + NTRLFC
      NSPANF = NSINFS + NCOSFS + NEXPFS + NLEDFS + NTRLFS
      NAREAF = NSINFA + NCOSFA + NEXPFA + NLEDFA + NTRLFA
      NTOTFU = NCMBRF + NSPANF + NAREAF
      NTKINK = NXFLR  + NYFLR

      IF (NCHCK1 /= 0 .AND. NCHCK2 /= 0) THEN
         WRITE (IWRIT,*) 'Y & Y" may not be modified simultaneously.'
         GO TO 900
      END IF

      IF (NTWING + NTOTFU + NTKINK /= NDV) THEN
         WRITE (IWRIT,*) 'Total # design variables did not match NDV.'
         GO TO 900
      END IF

      ORDERED = .TRUE.
      DO I = 2, NDV
         IF (IGROUP(I) < IGROUP(I-1)) ORDERED = .FALSE.
      END DO

      IF (ORDERED) GO TO 999

      WRITE (IWRIT, '(/, 1X, A, A, /, (4X, A))')
     >   'RDVARS:  The design variables must be grouped ',
     >   'in the following order:',
     >   'Planform', 'Wing Y', 'Wing Y"', 'Wing Y/Variable Center',
     >   'Flap/Slat Rotation or Shear', 'Body Camber',
     >   'Body Spanwise/Vertical/Radial', 'Body Area',
     >   'Floor Kink Point'
      GO TO 900

  800 WRITE (IWRIT, '(/, A, I5, 2X, A, 2X, A)')
     >   ' RDVARS: Unknown design variable or lofting variable: ',
     >   I, TYPE, TYPE2

  900 STOP

  999 RETURN

      END SUBROUTINE RDVARS

C+------------------------------------------------------------------------------
C
      SUBROUTINE RDVLINE (IVAR, VTYPE, UTYPE, DOUP, DOLO,
     >                    WBUMP, XBUMP, XWMIN, XWMAX,
     >                    PBUMP, UBUMP, UBMIN, UBMAX,
     >                    KCEN, KMIN, KMAX, V, VSCALE, AITCH, BL, BU,
     >                    NBOUND, LUNRD, LUNWR)
C
C     RDVLINE is a specialized utility for reading one line of design variable
C     inputs.  It was prompted by the several options desirable for controlling
C     spanwise lofting.  Distinguishing between INTEGER and REAL for certain
C     inputs was preferred to adding another input column or making an awkward
C     extension to the use of UTYPE, but it required parsing the tokens on each
C     line and forgoing the convenience of list-directed input.  A side-benefit
C     is that the variable names can now be input as (say) SIN or 'SIN' both.
C
C     Most of the arguments deliberately match SYN87 nomenclature.  Use of the
C     IVAR argument was considered preferable to passing VTYPE(IVAR), etc., in
C     the call.
C
C     01/30/97  D.Saunders  Introduced as part of generalizing spanwise lofting.
C
C     AUTHOR:  David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   IVAR, KCEN(*), KMIN(*), KMAX(*), NBOUND, LUNRD, LUNWR
      REAL
     >   DOUP(*), DOLO(*),
     >   WBUMP(*), XBUMP(*), XWMIN(*), XWMAX(*),
     >   PBUMP(*), UBUMP(*), UBMIN(*), UBMAX(*),
     >   V(*), VSCALE(*), AITCH(*), BL(*), BU(*)
      CHARACTER
     >   VTYPE(*)*(*), UTYPE(*)*(*)

C     Local constants:

      INTEGER
     >   IFLAG
      CHARACTER
     >   BLANK * 1, SEPS * 2

      PARAMETER
     >  (IFLAG = -999,  ! Flag assigned to KCEN/KMIN/KMAX if REALs were found
     >   BLANK = ' ',
     >   SEPS  = ' ''') ! I.e., blank or ' as needed for design variable names

C     Local variables:

      INTEGER
     >   FIRST, IER, IV, LAST, MARK, MARK0, N
      CHARACTER
     >   FORM*4,
     >   LINE*131 ! 131 allows a carriage control character in error messages

C     Execution:

      READ (LUNRD, '(A)') LINE

C     Trying to parse the items in a loop over all items would be awkward,
C     so we look for each in series, treating certain items specially.

C     The first item is not required to be the correct variable number.
C     Parse it but don't use it:

      FIRST = 1
      LAST = LEN (LINE)

      CALL SCAN2 (LINE, BLANK, FIRST, LAST, MARK)

      IV = IVAR ! As opposed to what was on the line

C     Chordwise/axial variable name - ignore previously-needed apostrophes:

      FIRST = MARK + 2

      CALL SCAN2 (LINE, SEPS, FIRST, LAST, MARK)

      READ (LINE (FIRST : MARK), '(A)', ERR=900) VTYPE(IV)

C     Spanwise/circumferential variable name similarly:

      FIRST = MARK + 2

      CALL SCAN2 (LINE, SEPS, FIRST, LAST, MARK)

      READ (LINE (FIRST : MARK), '(A)', ERR=900) UTYPE(IV)

C     Upper & lower surface switches (0., 1., or other possible REALs):

      CALL RDREAL (LINE, FIRST, LAST, MARK, DOUP(IV), IER)
      IF (IER /= 0) GO TO 900

      CALL RDREAL (LINE, FIRST, LAST, MARK, DOLO(IV), IER)
      IF (IER /= 0) GO TO 900

C     Streamwise width/power/exponent, center, min, and max:

      CALL RDREAL (LINE, FIRST, LAST, MARK, WBUMP(IV), IER)
      IF (IER /= 0) GO TO 900

      CALL RDREAL (LINE, FIRST, LAST, MARK, XBUMP(IV), IER)
      IF (IER /= 0) GO TO 900

      CALL RDREAL (LINE, FIRST, LAST, MARK, XWMIN(IV), IER)
      IF (IER /= 0) GO TO 900

      CALL RDREAL (LINE, FIRST, LAST, MARK, XWMAX(IV), IER)
      IF (IER /= 0) GO TO 900

C     Spanwise or circumferential width/power/exponent:

      CALL RDREAL (LINE, FIRST, LAST, MARK, PBUMP(IV), IER)
      IF (IER /= 0) GO TO 900

C     The center, min, and max for spanwise WING variables require
C     distinguishing between INTEGER and REAL inputs.  These items
C     will be REAL for BODY variables.  The XFLR/YFLR variables also
C     make use of KCEN for entering floor kink point #, which should
C     therefore be INTEGER.

      MARK0 = MARK ! So we can reparse once we know if it's REAL or not
      FIRST = MARK + 2

      CALL SCAN2 (LINE, BLANK, FIRST, LAST, MARK)

      IF (INDEX (LINE (FIRST : MARK), '.') > 0) THEN

C        A decimal point was found - disable the INTEGER possibilities
C        and treat all three items as REALs:

         KCEN(IV) = IFLAG
         KMIN(IV) = IFLAG
         KMAX(IV) = IFLAG

         MARK = MARK0

         CALL RDREAL (LINE, FIRST, LAST, MARK, UBUMP(IV), IER)
         IF (IER /= 0) GO TO 900

         CALL RDREAL (LINE, FIRST, LAST, MARK, UBMIN(IV), IER)
         IF (IER /= 0) GO TO 900

         CALL RDREAL (LINE, FIRST, LAST, MARK, UBMAX(IV), IER)
         IF (IER /= 0) GO TO 900

      ELSE ! An INTEGER was found in the KCEN position; assume 3 integers

         N = MARK - FIRST + 1
         FORM = '(IN)'
         WRITE (FORM (3 : 3), '(I1)') N ! 0 < N < 10 assumed
         READ (LINE (FIRST : MARK), FORM, ERR=900) KCEN(IV)

         FIRST = MARK + 2
         CALL SCAN2 (LINE, BLANK, FIRST, LAST, MARK)
         N = MARK - FIRST + 1
         FORM = '(IN)'
         WRITE (FORM (3 : 3), '(I1)') N ! 0 < N < 10 assumed
         READ (LINE (FIRST : MARK), FORM, ERR=900) KMIN(IV)

         FIRST = MARK + 2
         CALL SCAN2 (LINE, BLANK, FIRST, LAST, MARK)
         N = MARK - FIRST + 1
         FORM = '(IN)'
         WRITE (FORM (3 : 3), '(I1)') N ! 0 < N < 10 assumed
         READ (LINE (FIRST : MARK), FORM, ERR=900) KMAX(IV)

      END IF

C     Design variable value, scale factor, and differencing interval:

      CALL RDREAL (LINE, FIRST, LAST, MARK, V(IV), IER)
      IF (IER /= 0) GO TO 900

      CALL RDREAL (LINE, FIRST, LAST, MARK, VSCALE(IV), IER)
      IF (IER /= 0) GO TO 900

      CALL RDREAL (LINE, FIRST, LAST, MARK, AITCH(IV), IER)
      IF (IER /= 0) GO TO 900

C     Variable bounds?

      IF (NBOUND == -1) THEN

         CALL RDREAL (LINE, FIRST, LAST, MARK, BL(IV), IER)
         IF (IER /= 0) GO TO 900

         CALL RDREAL (LINE, FIRST, LAST, MARK, BU(IV), IER)
         IF (IER /= 0) GO TO 900

      END IF

      GO TO 999

  900 WRITE (LUNWR, '(/, 1X, A, I4, A)')
     >   'RDVLINE: Error reading design variable inputs.  Variable #:',
     >   IVAR, ';  offending line:', LINE

      STOP

  999 RETURN

C     Internal procedures for RDVLINE:

      CONTAINS

C+----------------------------------------------------------------------
C
      SUBROUTINE RDREAL (LINE, FIRST, LAST, MARK, VALUE, IER)
C
C     RDREAL is a somewhat-specialized utility needed to remove much
C     repetition from RDVLINE.  It is intended for parsing a REAL value
C     from the next item on a line, which may or may not contain a
C     decimal point.  MARK is assumed to be input as pointing to the
C     end of the previous item.
C
C     01/30/97  D.Saunders  Introduced with RDVLINE.
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      CHARACTER * (*) LINE
      INTEGER   FIRST, LAST, MARK, IER
      REAL      VALUE

C     Local variables:

      INTEGER   N
      CHARACTER FORM * 7

C     Execution:

      FIRST = MARK + 2

      CALL SCAN2 (LINE, ' ', FIRST, LAST, MARK)

      N = MARK - FIRST + 1

      IF (N < 10) THEN
         FORM = '(FN.0)'  ! Handles 1, 1., or 1.0, etc.
         WRITE (FORM (3 : 3), '(I1)') N
      ELSE
         FORM = '(FNN.0)' ! Needed for -9.999E+03, say
         WRITE (FORM (3 : 4), '(I2)') N
      END IF

      READ (LINE (FIRST : MARK), FORM, ERR=900) VALUE
      IER = 0
      GO TO 999

  900 IER = 1

  999 RETURN

      END SUBROUTINE RDREAL

C+----------------------------------------------------------------------
C
      SUBROUTINE SCAN2 (STRING, SEPS, FIRST, LAST, MARK)
C
C *** SCAN2 ought to be separated as an unchanging library-type utility,
C     comparable to the many numerics utilities used by this program.
C     However, it is not numeric and is the only one of its kind needed
C     so far, so it is kept with RDVLINE and RDREAL (above) for now.
C
C     Description and usage:
C
C           Looks for non-blank fields ("tokens") in a string, where the
C        fields are of arbitrary length and separated by any of a set of
C        user-specified separators (e.g., blanks or commas).  The position
C        of the end of the first token is also returned so that this
C        routine may be conveniently used within a loop to process an
C        entire line of text.
C
C           The procedure examines a substring, STRING (FIRST : LAST), which
C        may of course be the entire string (in which case just call SCAN2
C        with FIRST = 1 and LAST = LEN (STRING) ).  The indices returned
C        are relative to STRING itself, not the substring.
C
C
C     Arguments:
C
C        Name    Dimension  Type  I/O/S  Description
C        STRING    *         C    I      Text string containing data to be
C                                        scanned.
C        SEPS      *         C    I      String containing the separators.
C                                        Each character in SEPS counts as a
C                                        token delimiter.
C        FIRST               I    I/O    Index of first character of interest
C                                        in STRING.  FIRST >= 1 is assumed.
C                                        Output is index of the beginning of
C                                        the first non-separator, or 0 if no
C                                        token was found.
C        LAST                I    I/O    Index of last character of interest
C                                        in STRING.  LAST <= LEN (STRING) is
C                                        assumed.  Output is index of the end
C                                        of the last non-separator, or 0 if no
C                                        token was found.
C        MARK                I      O    Points to the end of the first token
C                                        found in the specified portion of
C                                        STRING.  Use FIRST = MARK + 2 for
C                                        further searches in the same string.
C                                        MARK is set to 0 if no token was found.
C
C
C     Notes:
C
C        (1)  FIRST >= 1 and LAST <= LEN (STRING) are assumed upon entry,
C             because these are the obvious bounds needed by the calling
C             program on the first scan of a given string - no need to
C             evaluate LEN (STRING) here with every call, while use of
C             MIN and MAX tends to disguise programming errors.
C
C        (2)  Entering with FIRST > LAST is an eventual normal consequence
C             of setting FIRST = MARK + 2 for further scans of the same string.
C             The zero DO loop iteration count feature of the language is
C             taken advantage of here - the result is simply to drop
C             out the bottom with FIRST = LAST = MARK = 0.
C
C        (3)  This version adjusts LAST using a BACKward search, which is
C             degenerate for all tokens after the first.  (Doing it in
C             the one forward loop means unnecessary repeated tokenizing
C             to find the end.)
C
C        (4)  Zeroing all of FIRST, LAST, and MARK (rather than just MARK)
C             when no token is found may, in retrospect, be undesirable in
C             that information is lost (esp. LAST), but there are too many
C             applications to change this now.
C
C
C     Author:  Robert Kennelly, Informatics General Corporation.
C
C
C     Development history:
C
C         4 Mar 1986    RAK    Variation of SCANNR, which is hard-coded
C                              (for historical reasons and a modest speed
C                              advantage) for separators BLANK, TAB,
C                              COMMA, COLON, and EQUAL.
C
C         5 May  1988   DAS    Reverse search used to find LAST; MAX, MIN,
C                              and LEN also eliminated (see Notes).
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   FIRST, LAST, MARK
      CHARACTER
     >   SEPS * (*), STRING * (*)

C     Local variables.

      INTEGER
     >   HEAD, I, J, NSEPS, TAIL
      LOGICAL
     >   FOUND

C     Execution.
C     ----------

      NSEPS = LEN (SEPS)
      HEAD = FIRST
      TAIL = LAST

      FIRST = 0
      LAST = 0
      MARK = 0
      FOUND = .FALSE.

C     Look at each character in STRING (HEAD : TAIL) from left to right
C     until a token is found.  (Drop through if LAST > FIRST.)

      DO 30, I = HEAD, TAIL

C        Is the character a separator?

         DO 10, J = 1, NSEPS
            IF (STRING (I : I) == SEPS (J : J)) GO TO 20
   10    CONTINUE

C           Not a separator.  Check for the beginning of the FIRST token.

            IF (.NOT. FOUND) THEN
               FIRST = I
               FOUND = .TRUE.
            END IF
            GO TO 30

   20    CONTINUE

C           We found a separator.

            IF (FOUND) THEN

C              We just passed the "trailing edge" of the first token.

               MARK = I - 1
               GO TO 40
            END IF
   30 CONTINUE

C     We reached the last character while still seeking the first token.
C     Either this is part of the first token, or MARK and LAST are still zero.

      IF (FOUND) THEN
         MARK = TAIL
         LAST = TAIL
      END IF

      GO TO 99

   40 CONTINUE

C     We found the first token but haven't reached the end yet to adjust LAST.

      DO 60, I = TAIL, MARK + 1, -1

C        Keep searching as long as we have a separator.

         DO 50, J = 1, NSEPS
            IF (STRING (I : I) == SEPS (J : J)) GO TO 60
   50    CONTINUE

C           Not a separator, so we've found LAST.

            LAST = I
            GO TO 99

   60 CONTINUE

C     Dropped through, so there were no more tokens.

      LAST = MARK

   99 RETURN

      END SUBROUTINE SCAN2

      END SUBROUTINE RDVLINE

C***********************************************************************
C
      SUBROUTINE REORDER (I1, I2, X, INDICES, NNEW)
C
C     Utility which searches X(I1:I2) for nonmonotonicity, and returns
C     INDICES(1:NNEW) pointing to strictly monotonic Xs.  The order of
C     X(I1) and X(I1+1) is assumed to be correct.
C
C     08/09/97  D.Saunders  Initial implementation, needed to suppress
C                           flap/slat points inside a wing section.
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN)  :: I1, I2
      REAL,    INTENT (IN)  :: X(1:I2)
      INTEGER, INTENT (OUT) :: INDICES(*), NNEW

C     Local variables:

      INTEGER  I, J, L
      REAL     ARROW

C     Execution:

      L = I1
      ARROW = SIGN (1., X(L+1) - X(L))
      INDICES(1) = L
      J = 1

      DO I = I1 + 1, I2
         IF (X(I) * ARROW > X(L) * ARROW) THEN
            L = L + 1
            J = J + 1
            INDICES(J) = I
         END IF
      END DO

      NNEW = J

      END SUBROUTINE REORDER

C***********************************************************************
C
      SUBROUTINE RIPPER (Z, ZCEN, ZMIN, ZMAX, PWR, HEIGHT)
C
C     The name is all that's left of James's bizarre original utility for
C     linear/polynomial spanwise lofting of the wing section perturbations.
C     The index-based scheme has been generalized to lofting between
C     arbitrary spanwise locations.  HEIGHT is returned between 0. and 1.
C
C***********************************************************************

C     Arguments:

      REAL, INTENT (IN)  :: Z, ZCEN, ZMIN, ZMAX, PWR
      REAL, INTENT (OUT) :: HEIGHT

C     Constants:

      REAL, PARAMETER :: ONE = 1., ZERO = 0.

C     Execution:

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

C*******************************************************************************
C
      SUBROUTINE SETLCON (BIGBND)
C
C     SETLCON sets up the linear constraint matrix for appropriate geometric
C     constraints, and adjusts the linear constraint bounds.
C     The '0' forms of the constraints enable initial (unperturbed) T/Cs or
C     Y/Cs to be determined for future use via the table printed here.
C     (However, the linear constraint set-up is read from disk if the
C     optimizer is being restarted.  SETLCON is not called in this case.)
C
C     3/96  J.Reuther  Initial implementation for T/C constraints (lower bounds
C                      only).
C     5/96  DAS/RAK    Linear interpolation wasn't good enough for l.e. radius-
C                      type thickness constraints.
C      "    DAS/JJR    'THK0 ' case added;
C      "       "       'YUPR ', 'YLWR ', 'YUPR0', 'YLWR0' added.
C     6/96     "       Upper bounds provided for. WARNING: The linear constraint
C                      set-up must be done with XINIT/YINIT (the geometry being
C                      perturbed). Applying any initial nonzero variables then
C                      setting up with XPTRB/YPTRB is wrong.
C     7/96     "       SETLCON isn't called for an NPOPT restart. BL, BU, and
C                      AMAT had to be moved to /OPT/.
C     8/96    DAS      SFEVAL avoids interpolation of unit perturbations at
C                      constraint Xs.
C    11/96     "       James realized that body radius constraints can be
C                      linear if the body camber is not changing.
C*******************************************************************************

C     Global quantities:

      USE CONSTS
      USE FUSEL
      USE FUSEL2
      USE LOGIC
      USE LUNS
      USE OPT
      USE OPTIONS

C     Arguments:

      REAL, INTENT (IN) :: BIGBND ! Large number for no-constraint cases

C     Local constants:

      REAL      ONE, ZERO
      LOGICAL   NEW
      CHARACTER TIGHT*1

      PARAMETER
     >  (ONE   = 1.,
     >   ZERO  = 0.,
     >   TIGHT = 'M',    ! Monotonic fits
     >   NEW   = .TRUE.) ! New data for each call to LCSFIT

C     Local variables:

      CHARACTER TYPE*6, TYP1*1, TYP4*4
      LOGICAL   FIRST

C     Execution:

C     Zero the linear constraint matrix - only some of the variables
C     are applicable.

      AMAT(1:NCLIN,1:NDV) = ZERO

C     Determine the nonzero matrix elements:

      IF (.NOT. WINGONLY) THEN
         NJ     = NF(2) ! CABINR assumes regular geometry
         XB4    = -999.
         XBLONG = XF(NFSTN) - XF(1)
      END IF

      NBAD  = 0
      FIRST = .TRUE.

      DO M = 1, NCLIN

         X = XLCON(M)
         K = KLCON(M)
         BOUNDL = BLLIN(M)
         BOUNDU = BULIN(M)

         IF (LCTYPE(M)(1:3) /= 'CAB') THEN
            NUPPK = NUPPER(K)
            NLOWK = NLOWER(K)
            Z     = ZLEADI(K)
         END IF

         IF (LCTYPE(M)(1:3) == 'THK') THEN ! Thickness constraint

C           Find the unperturbed section T/C at this X/C:

            CALL LCSFIT (NLOWK, XINIT(1,K), YINIT(1,K),
     >                   NEW, TIGHT, 1, X, YLO, YLO)

            CALL LCSFIT (NUPPK, XINIT(NLOWK,K), YINIT(NLOWK,K),
     >                   NEW, TIGHT, 1, X, YUP, YUP)

            THINIT = YUP - YLO

            IF (LCTYPE(M) == 'THK0 ') THEN ! Impose initial T/C as a minimum
               BOUNDL   = THINIT
               BLLIN(M) = THINIT
               IF (FIRST) THEN
                  FIRST = .FALSE.
                  WRITE (IWRIT,1020)
               END IF
               WRITE (IWRIT,1030) M, ' ''THK0 ''', K, X, THINIT, 999., Z
            END IF

C           Adjust for the fact that CHANGES in T/C are being constrained:

            IF (BOUNDL /= -BIGBND) BOUNDL = BOUNDL - THINIT
            IF (BOUNDU /=  BIGBND) BOUNDU = BOUNDU - THINIT
            BU(M+NDV) = BOUNDU
            BL(M+NDV) = BOUNDL

C           Evaluate the shape functions at the constraint X/C:

            Z = ZLEADI(K)

            DO N = NCHCKS + 1, NCHCKS + NCHCK1

C              Determine a unit design variable's contribution at this
C              span station:   0. <= |DYMAX| <= 1.

               CALL LOFT (UTYPE(N), Z, UBUMP(N), UBMIN(N), UBMAX(N),
     >                    PBUMP(N), DYMAX, IWRIT)

               YUP = ZERO
               YLO = ZERO

               IF (DYMAX /= ZERO) THEN
                  IF (DOLO(N) /= ZERO) THEN
                     DL = DYMAX * SIGN (ONE, DOLO(N))

                     CALL SFEVAL (VTYPE(N), WBUMP(N), XBUMP(N),XWMIN(N),
     >                            XWMAX(N), DL, 1, 1, X, YLO, IWRIT)
                  END IF

                  IF (DOUP(N) /= ZERO) THEN
                     DU = DYMAX * SIGN (ONE, DOUP(N))

                     CALL SFEVAL (VTYPE(N), WBUMP(N), XBUMP(N),XWMIN(N),
     >                            XWMAX(N), DU, 1, 1, X, YUP, IWRIT)
                  END IF
               END IF

               AMAT(M,N) = YUP - YLO
            END DO

         ELSE IF (LCTYPE(M)(1:3) == 'YUP') THEN ! Upper Y/C constraint

C           Find the initial section Y/C at this X/C:

            CALL LCSFIT (NUPPK, XINIT(NLOWK,K), YINIT(NLOWK,K),
     >                   NEW, TIGHT, 1, X, YUP, YUP)

            IF (LCTYPE(M) == 'YUPR0') THEN ! Impose initial Y/C as a minimum
               BOUNDL   = YUP
               BLLIN(M) = YUP
               BOUNDU   = BIGBND
               IF (FIRST) THEN
                  FIRST = .FALSE.
                  WRITE (IWRIT,1020)
               END IF
               WRITE (IWRIT,1030) M, ' ''YUPR0''', K, X, YUP, 999., Z
            END IF

C           Adjust for the fact that CHANGES in Y/C are being constrained:

            IF (BOUNDL /= -BIGBND) BOUNDL = BOUNDL - YUP
            IF (BOUNDU /=  BIGBND) BOUNDU = BOUNDU - YUP
            BU(M+NDV) = BOUNDU
            BL(M+NDV) = BOUNDL

C           Evaluate each shape function at the constraint X/C:

            DO N = NCHCKS + 1, NCHCKS + NCHCK1

               CALL LOFT (UTYPE(N), Z, UBUMP(N), UBMIN(N), UBMAX(N),
     >                    PBUMP(N), DYMAX, IWRIT)

               YUP = ZERO

               IF (DYMAX /= ZERO) THEN
                  IF (DOUP(N) /= ZERO) THEN
                     DU = DYMAX * SIGN (ONE, DOUP(N))

                     CALL SFEVAL (VTYPE(N), WBUMP(N), XBUMP(N),XWMIN(N),
     >                            XWMAX(N), DU, 1, 1, X, YUP, IWRIT)
                  END IF
               END IF

               AMAT(M,N) = YUP
            END DO

         ELSE IF (LCTYPE(M)(1:3) == 'YLW') THEN ! Lower Y/C upper bound

C           Find the initial section Y/C at this X/C:

            CALL LCSFIT (NLOWK, XINIT(1,K), YINIT(1,K),
     >                   NEW, TIGHT, 1, X, YLO, YLO)

            IF (LCTYPE(M) == 'YLWR0') THEN ! Impose initial Y/C as a maximum
               BOUNDL   = -BIGBND
               BOUNDU   = YLO
               BULIN(M) = YLO
               IF (FIRST) THEN
                  FIRST = .FALSE.
                  WRITE (IWRIT,1020)
               END IF
               WRITE (IWRIT,1030) M, ' ''YLWR0''', K, X, -999., YLO, Z
            END IF

            IF (BOUNDL /= -BIGBND) BOUNDL = BOUNDL - YLO
            IF (BOUNDU /=  BIGBND) BOUNDU = BOUNDU - YLO
            BU(M+NDV) = BOUNDU
            BL(M+NDV) = BOUNDL

C           Evaluate each shape function at the constraint X/C:

            DO N = NCHCKS + 1, NCHCKS + NCHCK1

               CALL LOFT (UTYPE(N), Z, UBUMP(N), UBMIN(N), UBMAX(N),
     >                    PBUMP(N), DYMAX, IWRIT)

               YLO = ZERO

               IF (DYMAX /= ZERO) THEN
                  IF (DOLO(N) /= ZERO) THEN
                     DL = DYMAX * SIGN (ONE, DOLO(N))

                     CALL SFEVAL (VTYPE(N), WBUMP(N), XBUMP(N),XWMIN(N),
     >                            XWMAX(N), DL, 1, 1, X, YLO, IWRIT)
                  END IF
               END IF

               AMAT(M,N) = YLO
            END DO

         ELSE IF (LCTYPE(M)(1:3) == 'CAB') THEN ! Cabin radius constraint

            ISEG = K / 100 ! K = KLCON(M) = 100*floor segment + radial angle #
            JANG = K - 100 * ISEG

C           Note: We could search XFLOORI(*) to find ISEG, but inputting it
C           is more consistent with the nonlinear form, which does not
C           search either (too wasteful in the case of many derivatives).

C           Find the initial cabin radius corresp. to this polygon point:

            RPOLYMOD = ONE ! Kludge to signal recentering without another arg.

            CALL CABINR (X, XB4, MXJBODY, NJ, NFSTN, XF, YFINIT, ZFINIT,
     >                   NFLOOR, XFLOORI, YFLOORI,
     >                   RBODY(JANG,ISEG), THBODY(JANG,ISEG),
     >                   TCABIN, YCABIN, ZCABIN, TINT, RPOLYMOD, RADIUS,
     >                   IWRIT, ISEG, IER)

            IF (LCTYPE(M) == 'CABI0') THEN ! Initial radius is the bound
               BOUNDL   = RADIUS
               BLLIN(M) = RADIUS
               BOUNDU   = BIGBND
               IF (FIRST) THEN
                  FIRST = .FALSE.
                  WRITE (IWRIT,1020)
               END IF

               WRITE (IWRIT,1040) M, ' ''CABI0''', K, X,RADIUS,999.,999.

            ELSE                           ! Recentered radius is the bound
               BOUNDL   = RPOLYMOD
               BLLIN(M) = BOUNDL
            END IF

C           Adjust for the fact that CHANGES in radius are being constrained:

            IF (BOUNDL /= -BIGBND) BOUNDL = BOUNDL - RADIUS
            IF (BOUNDU /=  BIGBND) BOUNDU = BOUNDU - RADIUS
            BU(M+NDV) = BOUNDU
            BL(M+NDV) = BOUNDL

C           Evaluate each shape function at the constraint X and Theta:

            DO N = NTWING + NCMBRF + 1, NTWING + NCMBRF + NSPANF
               TYPE = VTYPE(N)  ! 'BSxxxx', 'BVxxx', or 'BRxxx'; V = axial
               TYP1 = TYPE(2:2) ! 'S', 'V', or 'R'
               IF (TYP1 /= 'R') CYCLE

C              Axial direction:

               XDD = ((X-XF(1))/XBLONG-XWMIN(N)) / (XWMAX(N)-XWMIN(N))
               IF (XDD < ZERO .OR. XDD > ONE) CYCLE

               TYP4  = TYPE(3:6) ! 'xxxx'; SFEVAL expects NAME*4
               DRMAX = ZERO      ! Max. dR at this X; SFEVAL adds

               CALL SFEVAL (TYP4, WBUMP(N), XBUMP(N), ZERO, ONE, XBLONG,
     >                      1, 1, XDD, DRMAX, IWRIT)

C              Circumferential direction:

               UMIN  = UBMIN(N)
               UMAX  = UBMAX(N)
               UINT  = (TINT - UMIN) / (UMAX - UMIN)
               IF (UINT < ZERO .OR. UINT > ONE) CYCLE

               RBUMP = ZERO ! SFEVAL adds

               CALL SFEVAL (UTYPE(N), PBUMP(N), UBUMP(N), UMIN, UMAX,
     >                      DRMAX, 1, 1, UINT, RBUMP, IWRIT)

               AMAT(M,N) = RBUMP

            END DO ! Next body radius design variable

         ELSE
            NBAD = NBAD + 1
            WRITE (IWRIT,1010) 'Bad linear constraint:', M, LCTYPE(M)
         END IF
      END DO

      IF (NBAD > 0) THEN
         WRITE (IWRIT,1010) 'Aborting.'
         STOP
      END IF

      IF (.NOT. FIRST) THEN
         NBAD = 0
         DO N = NCHCKS + 1, NCHCKS + NCHCK1
            IF (V(N) /= ZERO) NBAD = NBAD + 1
         END DO
         IF (NBAD > 0) THEN
            WRITE (IWRIT,1010)
     >        'Misuse of ''0'' forms with nonzero variables - aborting.'
            STOP
         END IF
      END IF

      RETURN

 1010 FORMAT (/, ' SETLCON: ', A, I5, 3X, A)
 1020 FORMAT (/, '   #  LCTYPE KLCON      XLCON         BL         BU',
     >           '         Z')
 1030 FORMAT (1X, I3, A8, I6, 3F11.6, F10.4)
 1040 FORMAT (1X, I3, A8, I6, F11.2, F11.4, 2F10.0)

      END SUBROUTINE SETLCON

C+----------------------------------------------------------------------
C
      SUBROUTINE SFEVAL (NAME, EX, XC, XA, XB, YMULT, I1, I2, X, Y, LUN)
C
C        SFEVAL evaluates a specified shape function at X (I1:I2) and
C     adds the result to the corresponding Y.  X need not be normalized
C     because the option to apply the shape function to just part of the
C     data range, [XA, XB], forces normalization locally.  XC, however,
C     refers to the normalized sub-range, so it must be normalized.
C
C        These are the Hicks-Henne shape functions used for aerodynamic
C     design.  SFEVAL is an adaptation of BEVAL from program PROFILE in
C     the form most convenient for design code SYN87, where a modular
C     form is needed for precise set-up of linear constraints as well
C     as for the geometry perturbations.
C
C        Distinguishing wing and body applications should be done at
C     the calling level to avoid too many similar cases here.
C
C     06/21/96    DAS    Adaptation of BEVAL and PERTURB routines.
C     08/22/96     "     Added SINF. Apart from the SIN* functions,
C                        testing just 3 characters avoids problems.
C     08/24/96     "     LED & LEA now both work for LEADing;
C                        TRL & TRA  "    "    "   "  TRAILing.
C     12/06/96  DAS/JJR  Added COSL & COSR for body keel & crown,
C                        and inverted forms LCOS & RCOS.
C     06/03/97    DAS    Converting flap & slap angles from degrees to
C                        radians was using 180./PI instead of PI/180.
C     08/09/97     "     FLAP/SLAT now expect hinge X in same units as
C                        X(*) and Y(*) (not fraction of XB - XA).
C     08/10/99     "     Added SIN3 (fore/aft symmetry via SIN1 on [0,1])
C                        and   SIN4 (  "   "   "   "   "   "   " each half).
C
C     Author:  David Saunders, Sterling Software/NASA Ames, Moffett Field, CA
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

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

C     Local constants:

      REAL, PARAMETER :: HALF = 0.5, ONE = 1., TWO = 2., ZERO = 0.

C     Local variables:

      INTEGER    I
      REAL       AEXP, BEXP, CENTER, CENTER2, DEGRAD, EPSMCH, PT5LOG,
     >           ONEMC, PI, PIBY2, POWER, POWERL, POWERR, RN, RNM1, RPI,
     >           RRANGE, SCALE, TANGNT, THETA, WIDTH, X0, XI, XSCALE
      LOGICAL    FIRST
      CHARACTER  KEY4 * 4, KEY3 * 3

C     Storage:

      DATA       FIRST /.TRUE./
      SAVE       FIRST, DEGRAD, PI, PIBY2, RPI, PT5LOG

C     Statement functions:

      REAL
     >   EBUMP, SBUMP, PWR, WDTH, XNRM

      EBUMP (WDTH, PWR, XNRM) =
     >   XNRM ** PWR * (ONE - XNRM) * EXP (-WDTH * XNRM)

      SBUMP (WDTH, PWR, XNRM) = (SIN (PI * (XNRM ** PWR))) ** WDTH

C     Execution:

      IF (FIRST) THEN

         FIRST = .FALSE.

C        Protect sine bump exponentiation from slightly negative arguments
C        by calculating a value for PI such that SIN (PI) > 0 (barely).
C        This allows the limiting of X to [0, 1] that is necessary for
C        the sub-range capability to cover the sine protection too.

         EPSMCH = EPSILON (EPSMCH)

         PI = ASIN (ONE) * TWO

         DO I = 1, 999 ! 6 iterations observed for DEC AXP  (32 bits)
                       ! 5     "         "      "  CRAY C90 (64 bits)
C***        write (6,*) 'I   PI  SIN (PI): ', i, pi, sin (pi)
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

C     Minimize argument references in loops:

      WIDTH  = EX
      CENTER = XC ! Normalized
      X0     = XA
      RRANGE = ONE / (XB - X0)
      SCALE  = YMULT
      KEY4   = NAME (1:4) ! Avoid comparison of different-length strings
      KEY3   = KEY4 (1:3)

C     Treat the cases in ~ the most likely order they'll be encountered:

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

C        Assume [XA, XB] = [0, 1] and CENTER in (0, 0.5]:

         POWER = PT5LOG / LOG (CENTER)

         DO I = I1, I2
            XI = MAX (ZERO, MIN (X (I), ONE))
            IF (XI > HALF) XI = ONE - XI
            Y (I) = SCALE * SBUMP (WIDTH, POWER, XI) + Y (I)
         END DO

      ELSE IF (KEY4 == 'SIN4') THEN ! Fore/aft symmetry via SIN1 on each half

C        Assume [XA, XB] = [0, 0.5] and CENTER in (0, 1).  Left half:

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

C        Now the [0.5, 0.1] half with center <-- 1 - center

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

C        SIN is used instead of COS because SIN (PI) is protected above.

         XSCALE = PIBY2 * RRANGE
         X0 = TWO * X0 - XB

         DO I = I1, I2
            XI = MAX (PIBY2, MIN ((X (I) - X0) * XSCALE, PI))
            Y (I) = SCALE * (SIN (XI)) ** WIDTH + Y (I)
         END DO

      ELSE IF (KEY4 == 'COSR') THEN ! Half a cosine (or sine), peak at right

C        SIN is used instead of COS for consistency with COSL.

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

C        Reference:  Ramamoorthy, P. and Padmavathi, K.  "Airfoil Design
C        by Optimization" in J. Aircraft, Vol. 14 (Feb. 1977), 219-221.

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

C        "FLAP"-like function (shearing transformation only):
C        XC is the hinge point in the same frame as X & Y;
C        YMULT is the flap angle in DEGREEs (positive is "down").

         TANGNT = TAN (DEGRAD * SCALE)
         DO I = I1, I2
            XI = MAX (ZERO, X (I) - CENTER)
            Y (I) = Y (I) - XI * TANGNT
         END DO

      ELSE IF (KEY3 == 'SLA') THEN

C        "SLAT"-like function (shearing transformation only):
C        XC is the hinge point in the same frame as X & Y;
C        YMULT is the slat angle in DEGREEs (positive is "down").

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

C        Commonly used in conjunction with Wagner functions.

         DO I = I1, I2
            XI = MAX (ZERO, MIN ((X (I) - X0) * RRANGE, ONE))
            Y (I) = SCALE * XI + Y (I)
         END DO

      ELSE IF (KEY3 == 'SCA') THEN

C        Simple scaling - X(I) are probably ordinates.  Note that the
C        effect is arranged to be additive, not multiplicative, so that
C        scaling can be treated just as for the other functions when it
C        is included in a group for perturbing purposes.

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

C***********************************************************************
C
      SUBROUTINE WBFUDGE (NCALL, IBOD1, IBOD2, MXJBODY, MXIBODY, NFJ,
     >                    XF, YBODY, ZBODY, KWING2, MXIWING, MXKWING,
     >                    XWG, YWG, ZWG, NW, NLG, YSHELF, ZSHELF, IWRIT,
     >                    NEWBODY)
C
C     WBFUDGE protects wing/body intersection calculations by fudging
C     any body sections with interpolated wing geometry if the wing has
C     dropped below the body.  A safety margin (YSHELF) is included.
C
C     10/07/96  DAS  Initial implementation of a long-talked-about idea.
C     10/14/97   "   Fully-argument-driven version for multiblock code.
C     04/30/98   "   IER = 1 from first RIPPLE2D call led to undefined
C                    cell indices input for the next call.
C     12/04/98   "   Don't capture the shelf edge if prior J not fudged.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   NCALL,                  ! -3 means first call for this case
     >   IBOD1, IBOD2,           ! Range of body sections to treat
     >   MXJBODY,                ! Circumferential/axial dimensions of
     >   MXIBODY,                ! fuselage geometry arrays
     >   NFJ(MXIBODY)            ! Numbers of points for body sections
      REAL, INTENT (IN) ::
     >   XF(MXIBODY)             ! Body stations
      REAL, INTENT (INOUT) ::
     >   YBODY(MXJBODY,MXIBODY), ! Body section coordinates, possibly fudged
     >   ZBODY(MXJBODY,MXIBODY)  ! on output
      INTEGER, INTENT (IN) ::
     >   KWING2,                 ! Wing sections 1:KWING2 are treated
     >   MXIWING, MXKWING        ! I/K dimensions of wing section arrays
      REAL, INTENT (IN) ::
     >   XWG(MXIWING,MXKWING),   ! Wing geometry: sections 1:KWING2 must be
     >   YWG(MXIWING,MXKWING),   ! regular for 2-D interpolation
     >   ZWG(MXIWING,MXKWING)
      INTEGER, INTENT (IN) ::
     >   NW(MXKWING),            ! Numbers of wing section points (wrap)
     >   NLG(MXKWING)            ! Numbers on lower surfaces
      REAL, INTENT (IN) ::
     >   YSHELF                  ! Safety margin (e.g. 1") included in
                                 ! any fudged body section Y coordinate
      REAL, INTENT (IN) ::
     >   ZSHELF                  ! Fraction (e.g. 0.6) applied to the local
                                 ! body width to give a minimum shelf width
      INTEGER, INTENT (IN) ::
     >   IWRIT                   ! For diagnostics
      LOGICAL, INTENT (OUT) ::
     >   NEWBODY                 ! T means some body point(s) changed

C     Local constants:

      REAL, PARAMETER :: TOLER = 1.E-7

C     System functions:

      REAL, INTRINSIC :: EPSILON

C     Local variables:

      INTEGER
     >   I, IER, ILASTLO, ILASTUP, IWLO, IWUP, J, K, KLAST, KW, NFP,
     >   NFUDGE, NL, NTOT
      REAL
     >   EPS, P, Q, XLE, XPT, XTE, YJI, YLO, YUP, ZBODMAX, ZLAST, ZMIN,
     >   ZPT
      LOGICAL
     >   PREVIOUSJ

C     Procedures:

      EXTERNAL
     >   RIPPLE2D

C     Storage:

      SAVE
     >   EPS, NL, NTOT

C     Execution:

C     Check for feasibility of interpolating the wing lower surface?

      IF (NCALL == -3) THEN

         NL = NLG(1)
         NTOT = NW(1)

         DO K = 2, KWING2
            IF (NLG(K) /= NL)  GO TO 900
            IF (NW(K) /= NTOT) GO TO 900
         END DO

         EPS = MAX (5.* EPSILON (TOLER), TOLER) ! For p/q in [0,1] in RIPPLE2D
      END IF

      IWLO    = NL
      IWUP    = NL
      ILASTLO = IWLO ! In case first search fails
      KLAST   = 1

      NFUDGE = 0

      XLE = MIN (XWG(NL,1), XWG(NL,KWING2))
      XTE = MAX (XWG(1,1),   XWG(1,KWING2))

      DO I = IBOD1, IBOD2

         XPT   = XF(I)
         IF (XPT < XLE .OR. XPT > XTE) CYCLE ! Saves some wasted searching

         KW    = 1
         ZLAST = -1.
         NFP   = NFJ(I)

         ZBODMAX = ZBODY(1,I)
         DO J = 2, NFP
            ZBODMAX = MAX (ZBODY(J,I), ZBODMAX)
         END DO
         ZMIN = ZBODMAX * ZSHELF ! ZSHELF is now a fraction of ZMAX
         PREVIOUSJ = .FALSE.

         DO J = 1, NFP

            ZPT = ZBODY(J,I)

            IF (ZPT > ZMIN) THEN ! Not true for J = 1 or 2

               IF (PREVIOUSJ) THEN ! Previous J was fudged

C                 Try to capture the edge of the shelf more precisely:

                  P = (ZMIN - ZBODY(J-2,I)) /
     >                MAX (1.E-6, ZBODY(J-1,I) - ZBODY(J-2,I))
                  YJI = (1.- P) * YBODY(J-2,I) + P * YBODY(J-1,I)
                  ZBODY(J,I) = ZMIN
                  YBODY(J,I) = YJI
                  NFUDGE = NFUDGE + 1
               END IF

               EXIT ! Next I

            END IF

            IF (ZPT <= ZLAST) EXIT ! Next I ('cause we're above the water line)

            YJI = YBODY(J,I)

            CALL RIPPLE2D (MXIWING, MXKWING, 1, NL, 1, KWING2, XWG, ! Lower
     >                     ZWG, XPT, ZPT, IWLO, KW, EPS, P, Q, IER)

            IF (IER == 0) THEN
               YLO =
     >            (1.- Q)*((1.- P)*YWG(IWLO,KW)   + P*YWG(IWLO+1,KW))  +
     >                 Q *((1.- P)*YWG(IWLO,KW+1) + P*YWG(IWLO+1,KW+1))

               IF (YJI > YLO - YSHELF) THEN ! Quit if we're above upper wing too

                  IWUP = MIN (NTOT - 1, 2*NL - IWLO)

                  CALL RIPPLE2D (MXIWING, MXKWING, NL, NTOT, 1, KWING2, ! Upper
     >                           XWG, ZWG, XPT, ZPT, IWUP, KW, EPS,
     >                           P, Q, IER)
                  YUP =
     >            (1.- Q)*((1.- P)*YWG(IWUP,KW)   + P*YWG(IWUP+1,KW))  +
     >                 Q *((1.- P)*YWG(IWUP,KW+1) + P*YWG(IWUP+1,KW+1))

                  IF (YJI >= YUP) EXIT ! Next I; avoids TE trouble for large ZSHELF

                  YBODY(J,I) = YLO - YSHELF
                  NFUDGE = NFUDGE + 1
                  PREVIOUSJ = .TRUE.
               ELSE
                  PREVIOUSJ = .FALSE.
               END IF

            ELSE
               IWLO = ILASTLO ! Imagine a swept-back trailing edge crossing
               KW   = KLAST   ! a body section; hard to avoid other waste
               PREVIOUSJ = .FALSE.
            END IF

            ZLAST   = ZPT
            ILASTLO = IWLO
            KLAST   = KW

         END DO ! Next J

      END DO ! Next I

      NEWBODY = NFUDGE > 0

      IF (NCALL <= 0) THEN
         IF (IWRIT > 0) WRITE (IWRIT, '(/,A,I6)')
     >    ' WBFUDGE: # Body points fudged to be below the wing:', NFUDGE
      END IF

      RETURN

 900  WRITE (IWRIT, '(/,A)')
     >   ' WBFUDGE: Sections 1:KWING2 must have the same nos. of pts.'

      STOP

      END SUBROUTINE WBFUDGE

C***********************************************************************
C
      SUBROUTINE WING ()
C
C     WING reads the wing geometry and establishes the sections in
C     untwisted, normalized, wrap-around form.
C
C***********************************************************************

C     Global variables:

      USE CONSTS
      USE DIMENS
      USE LOGIC
      USE LUNS
      USE OPT
      USE OPTIONS
      USE SIZES
      USE TIT
      USE WNGSRF

C     Local constants:

      REAL, PARAMETER :: ZERO = 0.

C     Local variables:

      REAL      XP(MXIWING), YP(MXIWING), ZP(MXIWING), YTHICK(MXIWING)
      LOGICAL   THREECOLS
      CHARACTER BUFFER*80, ROUND*1

C     Execution:

      WRITE (IWRIT,'(/,A)') ' Geometry inputs:', ' '

      IF (YLESHIFT /= ZERO) THEN
         WRITE (IWRIT,'(/,A,/)')
     >      ' *** N.B.: Wing geometry will be shifted vertically. ***'
      END IF

      DO K = 1, 50 ! Transcribe the first portion to the printable output.
         READ  (LGEOM,   '(A)') BUFFER
         I = LEN_TRIM (BUFFER)
         WRITE (IWRIT,'(1X,A)') BUFFER(1:I)
      END DO

      WRITE (IWRIT,'(/,A)') ' ... et cetera ...'
      REWIND (LGEOM)

      READ (LGEOM,'(A)') TITLE
      READ (LGEOM,*)
      READ (LGEOM,*) ZSPAN, SREF, CREF, XREF, YREF, ZREF
      ASPECT = (2.* ZSPAN) ** 2 / (2.* SREF) ! Needed for plotting
      READ (LGEOM,*)
      READ (LGEOM,*) FNWSEC, ALPHAW
      NWSEC = FNWSEC

      IF (NWSEC < 2 .OR. NWSEC > MXKWING) THEN
         WRITE (IWRIT,'(/,1X,A,I4)')
     >      'WING:  Bad input no. of sections. Limits:  2  &', MXKWING
         STOP
      END IF

C     WINGIN = 2, 3, 4 means 2-, 3-, 2-column inputs

      THREECOLS = WINGIN == 3
      WRITE (IWRIT,1130)

      DO K = 1, NWSEC

         READ (LGEOM,*)
         READ (LGEOM,*) ZL, XL, YL, XYZSCALE, THKSCALE, AL, FSEC,
     >                  ROUNDED(K)

         XLEADI(K) = XL
         YL        = YL + YLESHIFT ! Eases fudging low or high-mounted wings
         YLEADI(K) = YL
         ZLEADI(K) = ZL
         AINIT(K)  = AL
         ALPHA     = AL + ALPHAW
         ALG(K)    = ALPHA
         CINIT(K)  = XYZSCALE
         TINIT(K)  = THKSCALE

         IF (FSEC /= ZERO) THEN
            READ (LGEOM,*)
            READ (LGEOM,*) YSYM, FNU, FNL
            READ (LGEOM,*)
            IF (YSYM > ZERO) THEN
               WRITE (IWRIT,*) 'No YSYM option supported - aborting.'
               STOP
            END IF

            NU = FNU
            NL = FNL
            N  = NU + NL - 1

            IF (N > MXIWING) THEN
               WRITE (IWRIT,'(/,1X,A,I4,A,I4)')
     >            'Too many pts. at wing section', K,'. Limit:', MXIWING
               STOP
            END IF
         END IF

         NUPPER(K)  = NU
         NUPPER0(K) = NU
         NLOWER(K)  = NL
         NLOWER0(K) = NL
         NTOTAL(K)  = N
         NTOTAL0(K) = N

         IF (FSEC /= ZERO) THEN

            IF (THREECOLS) THEN
               DO I = NL, N
                  READ (LGEOM,*) XP(I), YP(I), ZP(I)
                  XINIT(I,K) = XP(I)
                  YINIT(I,K) = YP(I)
                  ZINIT(I,K) = ZP(I)
               END DO
            ELSE
               DO I = NL, N
                  READ (LGEOM,*) XP(I), YP(I)
                  XINIT(I,K) = XP(I)
                  YINIT(I,K) = YP(I)
                  ZINIT(I,K) = ZL
               END DO
            END IF

            L = NL + 1
            READ (LGEOM,*)

            IF (THREECOLS) THEN
               DO I = 1, NL
                  J = L - I
                  READ (LGEOM,*) XP(J), YP(J), ZP(J)
                  XINIT(J,K) = XP(J)
                  YINIT(J,K) = YP(J)
                  ZINIT(J,K) = ZP(J)
               END DO
            ELSE
               DO I = 1, NL
                  J = L - I
                  READ (LGEOM,*) XP(J), YP(J)
                  XINIT(J,K) = XP(J)
                  YINIT(J,K) = YP(J)
                  ZINIT(J,K) = ZL
               END DO
            END IF

         ELSE

            DO I = 1, N
               XINIT(I,K) = XINIT(I,K-1)
               YINIT(I,K) = YINIT(I,K-1)
               ZINIT(I,K) = ZINIT(I,K-1)
            END DO

         END IF

C        Denormalize and apply any twist; find the true leading edge;
C        untwist about that, and renormalize precisely to ensure safe
C        application of shape functions on [0., 1.].  The normalized
C        and untwisted section may also be reused for the next K.

         CALL WINGSEC (K)


C        Tabulate initial wing section characteristics.
C        Untwisted section maximum thickness:

         NL = NLOWER(K) ! May have changed in WINGSEC
         NU = NUPPER(K)

         CALL GETTHICK (NL, NU, XINIT(1,K), YINIT(1,K), YTHICK,
     >                  THMAX(K), XTHMAX(K))

         THMAX(K) = THMAX(K) * THKSCALE


C        Trailing edge angle:

         CALL GETANGLE (NL, NU, XINIT(1,K), YINIT(1,K), THKSCALE,
     >                  TEANGL)


C        Total effective twist:

         COSTWT = COS (ALG(K) * DEG2RAD)
         SINTWT = SIN (ALG(K) * DEG2RAD)
         XTEDGE = (XINIT(1,K) + XINIT(N,K)) * 0.5
         YTEDGE = (YINIT(1,K) + YINIT(N,K)) * 0.5
         XTRAIL = XLEADI(K) + CINIT(K) *
     >            ((XTEDGE - XINIT(NL,K)) * COSTWT +
     >             (YTEDGE - YINIT(NL,K)) * SINTWT * THKSCALE)
         YTRAIL = YLEADI(K) + CINIT(K) *
     >            ((YTEDGE - YINIT(NL,K)) * COSTWT * THKSCALE -
     >             (XTEDGE - XINIT(NL,K)) * SINTWT)
         EFFTWT = ATAN2 (YLEADI(K) - YTRAIL, XTRAIL - XLEADI(K))*RAD2DEG

         IF (ROUNDED(K) == ZERO) THEN
            ROUND = 'S' ! Sharp leading edge
         ELSE
            ROUND = 'B' ! Blunt ...
         END IF

         WRITE (IWRIT,1140) K, ROUND, XLEADI(K), YLEADI(K), ZLEADI(K),
     >      CINIT(K), XTRAIL, YTRAIL, THKSCALE, THMAX(K), XTHMAX(K),
     >      AL, ALG(K), EFFTWT, TEANGL

      END DO ! Next wing section

      RETURN

 1130 FORMAT
     >  (/,' Initial wing summary (unperturbed but including YLESHIFT;',
     >   ' B = Blunt leading edge, S = Sharp):',
     >   //, ' Section     Xle       Yle        Zle   ',
     >   '   Chord        Xte        Yte    Thick  T/C Max',
     >   '   at X/C Lo.Twist To.Twist Ef.Twist TE.Ang.', /)
 1140 FORMAT (I4, A1, F11.4, F10.4, 4F11.4, F9.4, F9.6, 4F9.4, F8.3)

C     Internal procedure (subroutine WING):

      CONTAINS

!***********************************************************************
!
      SUBROUTINE WINGSEC (KSEC)
!
!     WINGSEC organizes an input wing section into a strictly untwisted,
!     normalized section, guaranteed to have X on [0, 1] by adjusting
!     the leading edge & chord if necessary.
!
!     03/25/97  D.Saunders  Introduced to simplify WING.
!     07/19/98      "       Internal procedure form.
!
!***********************************************************************

!     Argument:

      INTEGER, INTENT (IN) :: KSEC

!     Local variables:

      INTEGER
     >   I, ILW, K, NTOTK, NLOWK
      REAL
     >   ALPHA, CA, DSQ, DSQMAX, SA, SCALE, THICK, XLW, YLW, ZLW,
     >   XTRAIL, YTRAIL, XXW, YYW, ZZW

!     Execution:

      K = KSEC

!     Apply section twists about the GIVEN leading edge,
!     and denormalize (as in PERTURB):

      ALPHA = AINIT(K)
      THICK = TINIT(K)
      SCALE = CINIT(K)
      NTOTK = NTOTAL(K)
      NLOWK = NLOWER(K)
      XXW = XINIT(NLOWK,K)
      YYW = YINIT(NLOWK,K)
      ZZW = ZINIT(NLOWK,K)
      XLW = XLEADI(K) + XXW * SCALE
      YLW = YLEADI(K) + YYW * SCALE * THICK
      ZLW = ZLEADI(K)
      CA  = COS (ALPHA*DEG2RAD)
      SA  = SIN (ALPHA*DEG2RAD)

      DO I = 1, NTOTK
         XWG(I,K) = XLW + SCALE * ((XINIT(I,K) - XXW)*CA +
     >                             (YINIT(I,K) - YYW)*SA * THICK)
         YWG(I,K) = YLW + SCALE * ((YINIT(I,K) - YYW)*CA * THICK -
     >                             (XINIT(I,K) - XXW)*SA)
         ZWG(I,K) = ZLW + SCALE *  (ZINIT(I,K) - ZZW)
      END DO

!     Find the "true" leading edge as the furthest point from the
!     trailing edge in (X,Y) space:

      XTRAIL = (XWG(1,K) + XWG(NTOTK,K)) * 0.5
      YTRAIL = (YWG(1,K) + YWG(NTOTK,K)) * 0.5
      DSQMAX = 0.

      DO I = 2, NTOTK - 1
         DSQ = (XWG(I,K) - XTRAIL)**2 + (YWG(I,K) - YTRAIL)**2
         IF (DSQ > DSQMAX) THEN
            DSQMAX = DSQ
            ILW = I
         END IF
      END DO

!     Adjust the leading edge, chord, and local twist:

      NLOWER(K) = ILW
      NUPPER(K) = NTOTK - ILW + 1

      XLW = XWG(ILW,K)
      YLW = YWG(ILW,K)
      ZLW = ZWG(ILW,K)
      XLEADI(K) = XLW
      YLEADI(K) = YLW
      ZLEADI(K) = ZLW

      ALPHA    = ATAN2 (YLW - YTRAIL, XTRAIL - XLW)
      AINIT(K) = ALPHA * RAD2DEG
      ALG(K)   = AINIT(K) + ALPHAW

!     Untwist about the true leading edge, and remove the thickness
!     scaling (since the nominal section may be reused, and PERTURB
!     applies the scaling as a design variable):

      CA = COS (-ALPHA)
      SA = SIN (-ALPHA)
      THICK = 1./ THICK

      DO I = 1, NTOTK
         XINIT(I,K) =  XLW + (XWG(I,K) - XLW)*CA + (YWG(I,K) - YLW)*SA
         YINIT(I,K) = (YLW + (YWG(I,K) - YLW)*CA - (XWG(I,K) - XLW)*SA)
     >                * THICK
         ZINIT(I,K) = ZWG(I,K)
      END DO

      CINIT(K) = MAX (XINIT(1,K), XINIT(NTOTK,K)) - XLW

!     Normalize:

      SCALE = 1./ CINIT(K)
      YLW = YINIT(ILW,K)

      DO I = 1, NTOTK
         XINIT(I,K) = (XINIT(I,K) - XLW) * SCALE
         YINIT(I,K) = (YINIT(I,K) - YLW) * SCALE
         ZINIT(I,K) = (ZINIT(I,K) - ZLW) * SCALE
      END DO

      END SUBROUTINE WINGSEC

      END SUBROUTINE WING

C***********************************************************************
C
      SUBROUTINE YSMOOP (N, X, Y, EPSY, EXPY, A, B, C, S)
C
C     YSMOOP smooths airfoil coordinates by solving an equation of the
C     form
C                Ynew - eps d(dYnew/dS)/dS  =  Yold
C     where eps is a function of X arranged to help more towards the
C     troublesome trailing edge.
C
C     ??/??/??  J.Reuther   Original implementation.
C     09/15/94  D.Saunders  Properly remodularized; reworked to use
C                           TRDIAG instead of TRIAD, and less storage.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN)    :: N
      REAL,    INTENT (IN)    :: X(N)
      REAL,    INTENT (INOUT) :: Y(N)
      REAL,    INTENT (IN)    :: EPSY, EXPY
      REAL,    INTENT (OUT)   :: A(N), B(N), C(N), S(N)

C     Local variables:

      INTEGER  I
      REAL     TERM

C     Execution:

      S(1) = 0.
      DO I = 2, N
         S(I) = S(I-1) + SQRT ((X(I) - X(I-1))**2 + (Y(I) - Y(I-1))**2)
      END DO

      A(1) = 0.
      B(1) = 1.
      C(1) = 0.

      DO I = 2, N - 1
         TERM = -2.* EPSY * X(I) ** EXPY / (S(I+1) - S(I-1))
         A(I) = TERM / (S(I) - S(I-1))
         C(I) = TERM / (S(I+1) - S(I))
         B(I) = 1.- A(I) - C(I)
      END DO

      A(N) = 0.
      B(N) = 1.
      C(N) = 0.

      CALL TRDIAG (A, B, C, Y, Y, N)

      END SUBROUTINE YSMOOP

C***********************************************************************
C
      SUBROUTINE ZEROAREA (P1, P2, P3, TYPE, LUNERR, IER)
C
C PURPOSE:
C
C     ZEROAREA initializes evaluation of a zero-area, sine-type shape
C     function by determining the slope of the straight line which
C     combines with a full-period modified sine function to produce
C     zero area.  Two types of such a composite function arise from
C     the fact that the left and right (positive and negative) portions
C     of the modified sine function are not the same shape, under a
C     rotation about (P1, 0), except when P1 = 0.5.
C
C     For instance, if P1 = 0.3, the left/positive portion of the curve
C     is not the same shape as the right/negative portion for P1 = 0.7.
C     Since it is not clear which family of shapes is preferable, the
C     curves are grouped into two families symmetric about P1 = 0.5.
C     (Evidence suggests, however, that combinations of the first N
C     type 1. functions with power 1. provide the best least-squares
C     approximations to a typical airfoil.)
C
C     This initialization portion has been separated from the original
C     BEVAL2 routine.
C
C ARGUMENTS:
C     ARG    DIM  TYPE  I/O/S DESCRIPTION
C     P1            R     I   Zero-crossing value of X in [0, 1].
C     P2            R     I   Power of the sine function (real, > 0.).
C     P3            R     O   Slope of the linear weight fn. giving zero area.
C     TYPE          R     I   1. uses P1 <= 0.5 (or reflections if P1 > 0.5);
C                             2. uses P2 >= 0.5 (or reflections if P1 < 0.5).
C     LUNERR  -     I     I   LUNERR < 0 suppresses iteration printout;
C                             |LUNERR| is used for error messages.
C     IER     -     I     O   IER = 0 means successful initialization;
C                             IER < 0 means a failure in the zero-finding
C                                   or in the quadrature - probably the
C                                   no. of fn. evals. limit was exceeded.
C PROCEDURES:
C     QUANC8RC  Reverse-communication adaptive quadrature
C     ZERORC       "      "      "    1-D zero-finder
C
C HISTORY:
C   02/04/94  D.A.Saunders  Initial adaptation of BEVAL as BEVAL2 for
C                           some new functions proposed by James.
C   Pre-8/94  JJR           Separated BEVAL2 into ZEROAREA and two
C                           statement functions.
C   08/18/94  DAS           Restored the lost description, etc.
C
C AUTHORS: James Reuther, David Saunders, NASA Ames, Mt. View, CA
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      REAL,    INTENT (IN)  :: P1, P2
      REAL,    INTENT (OUT) :: P3
      REAL,    INTENT (IN)  :: TYPE
      INTEGER, INTENT (IN)  :: LUNERR
      INTEGER, INTENT (OUT) :: IER

C     Local constants:

      INTEGER
     >   MAXFUNZ
      REAL
     >   HALF, ONE, ZERO, LOGPT5, TWOPI, ABSERR, P3A, P3B, RELERR, TOL
      CHARACTER
     >   SUBNAME * 8
      PARAMETER
     >  (SUBNAME = 'ZEROAREA',
     >   MAXFUNZ = 50,
     >   HALF    =  5.E-1,
     >   ONE     =  1.E+0,
     >   ZERO    =  0.E+0,
     >   LOGPT5  = -0.693147180559945E+0,
     >   TWOPI   =  6.283185307179586E+0,
     >   ABSERR  = 1.E-6, ! Quadrature's absolute error tolerance
     >   RELERR  = 1.E-5, ! Tolerance relative to the integral size
     >   TOL     = ZERO,  ! Zero-finder's tolerance
     >   P3A     = -100., ! Search interval for P (3)
     >   P3B     = +100.)

C     Local variables:

      INTEGER
     >   I, I2, ISTATQ, ISTATZ, NFUNQ, NFUNZ
      REAL
     >   ERRESTQ, FLAGQ, FUNQ, FUNZ, HOLDZ(13), P4, POWER, SINE1, XQ
      LOGICAL
     >   TYPE1

C     Procedures:

      EXTERNAL
     >   QUANC8RC, ZERORC

C     Execution:

      IER = 0

C     Two variations of area-conserving modified sine function.
C     Each produces a family which is symmetric about P1 = 0.5.
C     Translate family type to sine-expression type:

      TYPE1 = TYPE == 1. .AND. P1 <= HALF .OR.
     >        TYPE == 2. .AND. P1 > HALF

      IF (TYPE1) THEN
         POWER = LOGPT5 / LOG (P1)
      ELSE
         POWER = LOGPT5 / LOG (ONE - P1)
      END IF

      IF (P1 == HALF) THEN
         P3 = ZERO
         GO TO 999
      END IF

C     We have a zero-finding iteration wrapped around an
C     adaptive quadrature iteration.  Use reverse-communication
C     to avoid passing data via common blocks.

      ISTATZ = 2 ! Initializes the zero-finder
      NFUNZ  = MAXFUNZ

   20 CONTINUE

         CALL ZERORC (P3A, P3B, P3, FUNZ, TOL, NFUNZ, SUBNAME,
     >                LUNERR, HOLDZ, ISTATZ)

         IF (ISTATZ < 0) THEN ! Probable application error
            WRITE (ABS (LUNERR), 1000) ISTATZ
            IER = ISTATZ
            GO TO 999

         ELSE IF (ISTATZ > 0) THEN
            ISTATQ = 0 ! Initialize the quadrature for this P3
            XQ = ZERO  ! Left end of interval

   30       CONTINUE

C              Evaluate the shape function (one of two types).
C              Arrange for higher powers to take the sign of sin ** 1.

               IF (TYPE1) THEN
                  SINE1 = SIN (TWOPI * XQ ** POWER)
                  FUNQ  = (P3 * (XQ - P1) + ONE) *
     >                    SIGN (ABS (SINE1) ** P2, SINE1)
               ELSE
                  SINE1 = SIN (TWOPI * (ONE - XQ) ** POWER)
                  FUNQ  = -(P3 * (XQ - P1) + ONE) *
     >                    SIGN (ABS (SINE1) ** P2, SINE1)
               END IF

               CALL QUANC8RC (XQ, FUNQ, ZERO, ONE, ABSERR, RELERR,
     >                        FUNZ, ERRESTQ, NFUNQ, FLAGQ, ISTATQ)

               IF (ISTATQ > 0)
     >      GO TO 30

            IF (ISTATQ /= 0) THEN
               WRITE (ABS (LUNERR), 1001) FLAGQ, NFUNQ, FUNZ
               IER = ISTATQ
               GO TO 999
            END IF

C           Else we have evaluated the integral successfully.

            GO TO 20 ! Continue the search for a zero integral

C        ELSE        ! ISTATZ = 0; P3 giving a zero integral has been found
         END IF

  999 RETURN

C     Formats.

 1000 FORMAT (/, ' *** ZEROAREA trouble (ZERORC):  status code = ', I2)
 1001 FORMAT (/, ' *** ZEROAREA trouble (QUANC8RC):', /, 5X,
     >        'FLAG = ', F7.3, '   NFUN = ', I5,  '   Integral = ', 1P,
     >        E12.3)

      END SUBROUTINE ZEROAREA

C*******************************************************************************
C
              SUBROUTINE GRIDWB (NCALL, FAIL)
C
C                        --------------------------------
C                        WING/BODY GRID GENERATION MODULE
C                        --------------------------------
C
C     The GRIDWB wing/body grid package originated as WBGRID (Lockheed), which
C     generated a C-H grid around a wing/body configuration using algebraic
C     starting guesses in combination with 2- and 3-D elliptic smoothers.  The
C     initial adaptation at NASA Ames was for the OPT67 wing/body design code.
C     Incorporation in the adjoint-based SYN87 followed, then the single-block
C     grid generation was adapted in stand-alone form as program CH_GRID, with
C     the grid-perturbation scheme disabled.  This was later resurrected so that
C     the GRIDWB module is now common to both CH_GRID and the design code.
C
C     Major enhancements have amounted to a 99% rewrite as an 8-subblock scheme.
C     An option for a C-O-type grid beyond the wing tip is included but has not
C     been exercised and is probably buggy.  The option to derive an N-S grid
C     from an initial Euler grid of the same dimensions forced packaging of
C     existing and new routines into 5 high-level modules for the design case.
C     For maintenance reasons, this is also the form used in stand-alone mode.
C
C     Only in CHECKWARP mode does the "warping" option add significantly to the
C     stand-alone memory requirements, but such a checking mode (especially for
C     the zero-perturbation case) is a valuable part of ensuring noise-free
C     gradient calculations in the design code, and is more easily tested in
C     stand-alone mode.  It could enable rapid derivation of a new grid from
C     an existing one if necessary (not implemented yet).
C
C
C     AUTHORS:    James Reuther              David Saunders
C     -------     RIACS/MCAT                 Sterling Software/Raytheon
C                 NASA Ames Research Center, MS 227-6
C                 Moffett Field, CA 94035-1000
C                 (650) 604-1516             (650) 604-1480
C
C                 Also: Stephen Edwards, U.C. Davis, contributed significantly
C                 to the improved TFI and elliptic smoothing utilities.
C
C
C     SPONSORS:   High Speed Aerodynamics Branch, NASA Ames Research Center
C     --------
C
C     COORDINATE SYSTEM:
C     -----------------
C
C                 The coordinate system is right-handed as used by the flow
C                 solvers of Antony Jameson such as FLO87 and AIRPLANE:
C
C                    I  &  X  Streamwise
C                    J  &  Y  Vertical
C                    K  &  Z  Spanwise
C
C                    I = 1 is at or beyond the lower trailing edge.
C                    J = 1 is at the wing surface.
C                    K = 1 is at the center plane or fuselage surface.
C
C                 Downstream, up, and towards the (left) tip are all positive.
C
C
C     SUMMARY OF GRID GENERATION ROUTINES:
C     -----------------------------------
C
C                      Top level grid generation routine
C                      ---------------------------------
C
C        GRIDWB:       READGRID = F is assumed by this top level; if
C                      READGRID = T, the application calls the appropriate
C                      second-level routines directly
C
C                      Second-level grid routines
C                      --------------------------
C
C        GRID_FACES:   Generates all subblock faces; some may be perturbed
C        GRID_TFI3D:   Initializes subblock volumes from scratch
C        GRID_ELLIP3D: Smooths the subblock volumes elliptically
C
C        REGRID_EDGES: Redistributes radial edges of subblock faces
C        REGRID_FACES: Redistributes radial lines of subblock faces
C        REGRID_VOLS:  Redistributes radial lines of subblock volumes
C
C        WARP_VOLS:    Perturbs subblock volumes given perturbed faces
C
C                      Lower-level grid routines
C                      -------------------------
C
C        BILINCUT:      Parametric bilinear form of PLBICUT
C        BLGRID:        Composite geometric + Vinokur distribution
C        CELLVOL:       Calculates cell volumes
C        CHECK_VOLS:    Calculates cell volumes; counts & reports volumes < 0.
C        CHORDS2D:      Arc lengths along a 2-space curve
C        CHORDS3D:      Arc lengths along a 3-space curve
C        CHORDSRF:      Arc lengths along a row of a 3-space surface
C        CHORDS_FACES:  Unnormalized radial arcs for subblock faces
C        CHORDS_J1:     Unnormalized arcs for the wing surface (J = 1 plane)
C        CHORDS_VOLUME: Unnormalized radial arcs for entire grid
C        CROWN:    Determines the crown line of fuselage
C        ELLIP2D:  2-D elliptic grid smoother
C        ELLIP3D:  3-D elliptic grid smoother
C        ELLIPQ3D: Elliptic smoother for 3-space grid block faces
C        EXCHANGE: Swaps C-grid smoothing controls for I1 & I2 boundaries
C        EXPDIS4:  Exponential stretching of points along a line
C        FIXBGRID: Adjusts the radial lines on the parametric body grid
C        FIXLEDGE: Adjusts lines forward of the leading edge, in case of slats
C        FIXWAKE:  Adjusts lines in the wake, in case of flaps
C        FMINRC:   1-D minimizer
C        FOILGRD:  Linear/quadratic/sine/cosine hybrid 1-D distribution
C        HTDIS4:   Vinokur distribution of points along a line
C        LCSFIT:   Local 1-D cubic spline interpolation
C        NULINE2D: Warps a 2-space line given new end points
C        OFFCROWN: Radial point distribution off crown/keel in symmetry plane
C        OFFLEDGE: Radial point distribution off leading edge ...
C        OFFTEDGE: ... and off trailing edge
C        OUTER:    Determines outer grid boundary
C        PARAM2D:  Parameterizes a regular surface grid as needed by PLBICUBE
C        PARAMXY:  Variation of PARAM2D used by ELLIP2D
C        PARAMXYZ: Parameterizes a 3-D regular volume grid as needed by WARP3D
C        PBILINT:  Parametric bilinear interpolation
C        PLBICUBE: Parametric local bicubic spline interpolation
C        PLBICUT:  Find Z (and U,V) for given X,Y on a surface, etc.
C        PLSCRV3D: Parametric interpolation for a 3-space line
C        PLXCUT:   Find T for given X on a parametric 2-space curve
C        RADIAL_BLEND:  Redistributes arcs given old arcs & two new edge arcs
C        RADIAL_INIT:   Sets up two sets of radial distribution controls
C        RADIAL_SET:    Activates the indicated set of radial controls
C        RADIAL_TFI2D:  Redistributes radial lines of a subblock volume
C        RECTIFY:  Adjusts wing surface grid sections (ILE at true LE, etc.)
C        REGRID_IPLANE: Redistributes radial lines of an I face interior
C        REGRID_KPLANE: Redistributes radial lines of a  K face interior
C        REGRID_LINE:   Redistributes a radial edge line
C        SMOOTH1D: Explicit smoothing of a vector
C        SMOOTHX:  Adjusts abscissas to include given one(s) smoothly
C        TFI2D:    Transfinite interpolation for a 2-D plane
C        TFI3D:    3-D transfinite interpolation with optimal Soni blending
C        TFINT3F:  3-function TFI utility used by TFI3D
C        TFIQ3D:   Transfinite interpolation for a 3-D surface
C        VINOKUR:  Simplifies use of HTDIS4 on descending data
C        VSHEET:   Deals with vortex sheet leaving wing trailing edge in C-mesh
C        WARP3D:   Grid perturbation utility for subblock volumes
C        WARPQ3D:  Grid perturbation utility for subblock faces
C        WARP_VOLS: Perturb all subblock volumes
C        WBINTR:   Reusable "wing/body" intersection utility
C        WFINTR:   Calculates the wing/fuselage intersection
C        WINSERT:  Checks input crank stations, possibly inserting sections
C        WREGUL:   Regularizes a range of wing sections
C        WSURF:    Wing surface grid generation
C        Z0PLANE:  Generates symmetry plane grid & body surface grid if present
C
C
C     HISTORY OF THE GRID GENERATION PACKAGE:
C     --------------------------------------
C
C                        -------------------
C                        OPT67/SYN87 Origins
C                        -------------------
C
C     ??/92-  J.Reuther  Original OPT67, using WBGRID wing/body grid generation
C     ??/93              from Lockheed.
C     08/94-  D.Saunders OPT67 overhauled to improve modularity and workspace
C     03/95              use.  The grid generation was essentially rewritten,
C                        including a body surface grid option using 2D para-
C                        ametric interpolation and smoothing, and the grid
C                        perturbation scheme is via the new WARPQ3D and WARP3D.
C     04/95       "      SYN87 updated with OPT67 improvements.
C     05/95-      "      Numerous grid refinements, including rounded-tip
C     06/95              options, and four-block C-H structure with artificial
C                        boundaries forward of the leading edge and at the tip
C                        to overcome negative cell volumes near the nose.
C     07/95-      "      Geometry & grid generation now use the flow solver
C     12/95              convention for Y and Z throughout.
C     02/96       "      The first N (e.g. 2) leading edge cranks are now
C                        captured in the surface grid.
C     04/96    JJR/DAS   Stage 2 of 3-stage WARP3D scheme is improved.
C     08/10-  S.Edwards/ Replaced TTM2D with ELLIP2D, and
C     08/16/96   DAS     replaced TTM3D with ELLIP3D.
C     08/20/96   DAS     Independent control of calls to ELLIP2D/3D.
C     08/28-      "      Various grid refinements inspired by James:
C     08/30/96           D1ROOT/TIP/KMAX replace D1WING; CRDSPL/Q/S/C are
C                        inputs for FOILGRD vs. CRDSPC for FOILGRID;
C                        the K = 1 edge of the back plane uses D1CROWN,
C                        and TFIQ3D is used on the whole back plane;
C                        efforts to use ELLIPQ3D on the sheet forward of
C                        the leading edge unsuccessful: irregular results.
C     09/11/96  SJE/DAS  TFI3D replaces TFIQ3D for initial volume grids.
C     09/13/96  DAS/JJR  Wing/body intersection calculations start both
C                        surfaces at the leading edge now.
C     09/19/96     "     D2ROOT, D2TIP introduced; 6-block scheme now.
C     10/02/96     "     Wing sections now need a ROUND flag to control
C                        gridding as two surfaces or as one wrap-around.
C
C                        ------------------------------
C                        CH_GRID Stand-alone Adaptation
C                        ------------------------------
C
C     03/11/97    DAS    Initial extraction of CH_GRID from SYN87-SB.
C     03/13/97  DAS/JJR  Introduced D3ROOT/TIP/KMAX and TEROOT/TIP/KMAX,
C                        although D3TIP and D3KMAX cannot be used without
C                        going to an 8-subblock scheme.
C     03/18/97    DAS    The inboard panel of the 2-panel wing grid option
C                        now has uniform spanwise distributions, after the
C                        Vinokur distributions were found liable to fail.
C     03/19/97     "     Added READGRID option after James found he needs
C                        to smooth large grids in pieces to stay within the
C                        debug queue CPU limit.
C     04/22/97  DAS/JJR  Introduced D1NOSE and D1TAIL to help N-S cases.
C     04/24/97     "     Introduced DWTAIL and DCTAIL.
C     05/02/97     "     JCR was being used in WINGONLY mode in OUTER and
C                        VOLGRID.  FIXBGRID needed smarter variation of the
C                        last radial increments from nose to tail.
C     05/09/97     "     Went from 6 subblocks to 8 (new OFFTEDGE routine).
C     05/30/97    DAS    The wing-only case needed X/Y/ZWALL copied to X()
C                        for OFFTEDGE purposes in Z0PLANE.
C     06/14/97     "     FIXTGRID appears redundant in the presence of the
C                        improved elliptic smoothing (in (u,v) space).
C     06/21/97  SJE/DAS  ELLIPQ3D improvement allows use on I = ILE plane;
C                        ELLIP3D work-space reduced; SINTERP eliminated.
C     06/29/97    DAS    Introduced TIPLE control, and adjusted interior
C                        increment DM of outer boundary for monotonicity.
C     07/11/97  SJE/DAS  ELLIP2D foreground control may be arc-based now.
C                        Distributions at I = 1, ILE, IL are 2-sided.
C                        D2JL introduced; TIPLE (and IC3) dispensed with.
C     07/27/97     "     (1 - "s")**power terms applied to all foreground
C                        decay terms in ELLIP2D, ELLIPQ3D, and ELLIP3D.
C     Aug. 1997   DAS    NBLAYER & RBLAYER introduced to make composite
C                        geometric + Vinokur edge distributions forward
C                        of leading edge and above & below trailing edge.
C     Sep. 1997    "   o Spanwise wing surface gridding is now arc-based,
C                        as would be needed by winglets.  This prompted
C                        regularizing the wing sections BEFORE the wing/
C                        body intersection, at the price of rectifying
C                        the leading edge index after the initial grid.
C                      o FIXWAKE and FIXLEDGE refinements account for body
C                        curvature better.  The planform crank stations
C                        are read instead of being detected automatically,
C                        because slat deflections can interfere.  ZPLANAR
C                        and NNROOT = 0 or 1 options have been eliminated:
C                        use NCRANK, ZCRANK(*) or NNROOT = 3.
C     May  1998    "     Made CH_GRID executable independent of the grid
C                        size and the number of variables and constraints.
C     June 1998    "   o Poor TFI3D behavior led to an I-plane boundary at
C                        I = IBL and IBU (but not for ELLIP3D).  Behavior
C                        is still poor, but a straight boundary in the
C                        symmetry plane off the crown/keel opposite the
C                        trailing edge does not help. Also, for N-S grids,
C                        DWTAIL must not be much bigger than D2ROOT, else
C                        TFI3D overshoots at the wake sheet.
C                      o RECTIFY no longer picks the smallest-X point as
C                        the leading edge because the discretization can
C                        produce wrinkles in the I=ILE sheet.  Locate the
C                        leading edge by minimizing X w.r.t. arc length.
C                        WREGUL no longer uses monotonic fits for X vs. T.
C                      o FIXLEDGE fiddling is reduced for small LE angles
C                        and avoided altogether for non-planar K sections.
C
C                        ----------------------------------------------
C                        GRIDWB Package Common to CH_GRID & Design Code
C                        ----------------------------------------------
C
C     July 1998   DAS    Top-level reorganization to incorporate a radial
C                        redistribution option and restore "warping" mode
C                        to the CH_GRID form so that the GRIDWB package is
C                        common to CH_GRID and the wing/body design code.
C                        The READGRID case is handled outside GRIDWB to help
C                        simplify the logic here and because of different
C                        interpretations by the two applications.
C
C     Sep. 1998    "     The design code case of reading a grid ready for a flow
C                        solution required introducing READGRID and NNACEL here.
C
C     Aug. 1999    "     Revised OFFTEDGE to allow rounded trailing edges;
C                        2-D interpolation utilities now return (p,q).
C
C*******************************************************************************

C     Global variables:

      USE LIM
      USE LOGIC
      USE LUNS
      USE OPTIONS
      USE SPACING
      USE XYZ
      USE XYZ0

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN)  :: NCALL
      LOGICAL, INTENT (OUT) :: FAIL

C     Local variables:

      LOGICAL FULLGRID

C     Execution:


      IF (READGRID) THEN ! Applies only to the design code

         IF (NNACEL .EQ. 0) THEN ! No more to do

            FAIL = .FALSE.

         ELSE ! NACSHIFT needs regularized wing sections and the surface grid

            CALL WSURF (NCALL, FAIL)

         END IF

         GO TO 999

      END IF


      IF (NCALL == -3) THEN        ! First call
         ALLOCATE (S0(IL,JL,KL,3)) ! For TFI3D, ELLIP3D, WARP[Q]3D, & REGRID_*
         NEEDX0 = .TRUE.           ! X0 is allocated only if we're perturbing
      END IF

      FULLGRID = (NCALL  == -3) .OR.                   ! First call
     >           (NNGRID ==  0) .OR.                   ! Full grid always
     >           (NNGRID ==  2 .AND. NCALL <= -2) .OR. ! End of line search
     >           (NNGRID ==  3 .AND. NCALL <=  0)      ! In line search too

C     Set up the subblock face grids.
C     If we're warping (FULLGRID = F), all faces except the wing and the body
C     may be done by perturbation given perturbed edges, which likewise are
C     perturbed given perturbed end points as far as possible.
C     If we're in indirect mode (REDISTRIBUTE = T) and warping (FULLGRID = F),
C     the final radial distributions apply throughout except for the body
C     surface grid, which initially uses the standard controls but employs
C     the *NS controls in the FIXBGRID call.

      CALL GRID_FACES (NCALL, FULLGRID, FAIL)

      IF (FAIL) GO TO 999


C     Set up the subblock volume grids:

      IF (FULLGRID) THEN

C        Initialize the subblock grid volumes:

         CALL GRID_TFI3D (NCALL)

C        Smooth the subblock grid volumes:

         CALL GRID_ELLIP3D (NCALL)


C        Redistribute interim radial lines?

         IF (REDISTRIBUTE) THEN

C           Activate the radial redistribution controls:

            CALL RADIAL_SET (.FALSE.) ! T means initial ("Euler") controls

C           Redistribution of all radial edges of all subblocks:

            CALL REGRID_EDGES (NCALL)

C           Radial redistribution of the subblock faces:

            CALL REGRID_FACES (NCALL)

C           Radial redistribution of the subblock volumes:

            CALL REGRID_VOLS (NCALL)

         END IF

C        Save a copy of the grid for ensuing perturbations?

         IF ((CHECKWARP) .OR.
     >       (.NOT. GRIDONLY .AND. NNGRID > 0 .AND. OPTIMIZING)) THEN

            IF (NEEDX0) THEN ! The first full grid is ready

               NEEDX0 = .FALSE.
               ALLOCATE (X0(IL,JL,KL,3))

            END IF

            X0 = X ! ... and parameterize it:

            CALL PARAMXYZ (1, IL, 1, JL, 1, KL, 1, IL, 1, JL, 1, KL,
     >                     X0(1,1,1,1), X0(1,1,1,2), X0(1,1,1,3), S0)
         END IF


      ELSE ! The grid is being perturbed:

C        Perturb the subblock volumes:

         CALL WARP_VOLS (NCALL)

      END IF

C     Check for negative cell volumes unless we're doing gradients:

      IF (NCALL <= 0 .OR. CHECKWARP) THEN

         CALL CHECK_VOLS ('GRIDWB', IL, JL, KL, HAND, X, IWRIT)

      END IF

  999 RETURN

      END SUBROUTINE GRIDWB

C*******************************************************************************
C
      SUBROUTINE BLGRID (N, D1, D2, NBLAYER, RBLAYER, X, LUNERR, IER)
C
C        BLGRID generates a specialized 1-D grid intended for radial lines
C     starting in a boundary layer and ending in the far field.  The first
C     and last X values are assumed to be in X(1) and X(N).  Points 1
C     through NBLAYER have spacing that varies geometrically according to
C     the ratio RBLAYER (or uniformly if RBLAYER is 1.0).  Spacing is
C     Vinokur-type beyond that.  X is probably arc-length (zero to STOTAL),
C     but could be decreasing if desired.  D1 & D2 are both positive.
C
C     Use NBLAYER <= 2 to suppress the special treatment.
C     See HTDIS4 for further details, including the LUNERR and IER arguments.
C
C     07/08/97  DAS  Effort to emulate the single-piece C grid result for
C                    the two-halves 2-D case forward of the leading edge.
C     08/06/97   "   Generalized to allow application off the trailing edge.
C     08/22/98   "   In-lined the VINOKUR call because it was easy.
C
C*******************************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN)    :: N, NBLAYER, LUNERR
      REAL,    INTENT (IN)    :: D1, D2, RBLAYER
      REAL,    INTENT (INOUT) :: X(N)
      INTEGER, INTENT (OUT)   :: IER

C     Local constants:

      REAL, PARAMETER :: ONE = 1., ZERO = 0.

C     Local variables:

      INTEGER  I, I1, NB
      REAL     DX, RANGE, XI1, XI2

C     Execution:

      DX = SIGN (D1, X(N) - X(1))
      NB = MAX (NBLAYER, 2) ! Force a single pass if special treatment
                            ! is suppressed
      DO I = 2, NB
         X(I) = X(I-1) + DX
         DX   = DX * RBLAYER
      END DO

      DX    = ABS (DX / RBLAYER)
      I1    = NB - 1
      XI1   = X(I1)
      XI2   = X(N)
      RANGE = XI2 - XI1

C     Overlap the two distributions by one cell:

C *** CALL VINOKUR (I1, N, DX, D2, X, LUNERR, IER) ! In-lined below

C     Normalized Vinokur distribution in ascending order:

      CALL HTDIS4 (.TRUE., ZERO, ONE, ABS (DX/RANGE), ABS (D2/RANGE),
     >             N - I1 + 1, X(I1), -LUNERR, IER)

C     Denormalize:

      DO I = I1, N - 1
         X(I) = XI1 + X(I) * RANGE
      END DO

      X(N) = XI2 ! Exactly

      END SUBROUTINE BLGRID

C**********************************************************************
C
      SUBROUTINE CELLVOL (IL, JL, KL, X, VOL, HAND)
C
C     Argument-driven portion of FLO87's METRIC routine for calculating
C     cell volumes, as needed for a grid quality check. HAND = +/-1 for
C     RH/LH xyz resp. VP1, VP3, & VP5 negative signs added empirically.
C
C**********************************************************************

      IMPLICIT   NONE

C     Arguments:

      INTEGER, INTENT (IN)  :: IL, JL, KL
      REAL,    INTENT (IN)  :: X(IL,JL,KL,3), HAND
      REAL,    INTENT (OUT) :: VOL(2:IL,2:JL,2:KL)

C     Local constants:

      REAL, PARAMETER :: FOURTH = 1./ 4., SIXTH = 1./ 6., R8TH = 1./ 8.

C     Local variables:

      INTEGER
     >   I, J, K, L, M, N
      REAL
     >   FACTOR, XP, YP, ZP, VP1, VP2, VP3, VP4, VP5, VP6,
     >   VOLPYM, XA, YA, ZA, XB, YB, ZB, XC, YC, ZC, XD, YD, ZD

C     Statement function (more efficient than internal procedure):

      VOLPYM (XA, YA, ZA, XB, YB, ZB, XC, YC, ZC, XD, YD, ZD) =
     >      (XP - FOURTH * (XA + XB + XC + XD))
     >   * ((YA - YC) * (ZB - ZD) - (ZA - ZC) * (YB - YD))
     >   +  (YP - FOURTH * (YA + YB + YC + YD))
     >   * ((ZA - ZC) * (XB - XD) - (XA - XC) * (ZB - ZD))
     >   +  (ZP - FOURTH * (ZA + ZB + ZC + ZD))
     >   * ((XA - XC) * (YB - YD) - (YA - YC) * (XB - XD))

C     VOLPYM = Volume of a pyramid with four-sided base, times 6
C     (adapted from Jameson code in FLO87).

C     Execution:

      FACTOR = HAND * SIXTH

      DO K = 2, KL
         N = K - 1
         DO J = 2, JL
            M = J - 1
            DO I = 2, IL
               L = I - 1
               XP =(X(I,J,K,1) +X(I,M,K,1) +X(I,M,N,1) +X(I,J,N,1)
     >             +X(L,J,K,1) +X(L,M,K,1) +X(L,M,N,1) +X(L,J,N,1))*R8TH
               YP =(X(I,J,K,2) +X(I,M,K,2) +X(I,M,N,2) +X(I,J,N,2)
     >             +X(L,J,K,2) +X(L,M,K,2) +X(L,M,N,2) +X(L,J,N,2))*R8TH
               ZP =(X(I,J,K,3) +X(I,M,K,3) +X(I,M,N,3) +X(I,J,N,3)
     >             +X(L,J,K,3) +X(L,M,K,3) +X(L,M,N,3) +X(L,J,N,3))*R8TH
               VP1 =-VOLPYM (X(I,J,K,1),X(I,J,K,2),X(I,J,K,3),
     >                       X(I,M,K,1),X(I,M,K,2),X(I,M,K,3),
     >                       X(I,M,N,1),X(I,M,N,2),X(I,M,N,3),
     >                       X(I,J,N,1),X(I,J,N,2),X(I,J,N,3))
               VP2 = VOLPYM (X(L,J,K,1),X(L,J,K,2),X(L,J,K,3),
     >                       X(L,M,K,1),X(L,M,K,2),X(L,M,K,3),
     >                       X(L,M,N,1),X(L,M,N,2),X(L,M,N,3),
     >                       X(L,J,N,1),X(L,J,N,2),X(L,J,N,3))
               VP3 =-VOLPYM (X(I,J,K,1),X(I,J,K,2),X(I,J,K,3),
     >                       X(I,J,N,1),X(I,J,N,2),X(I,J,N,3),
     >                       X(L,J,N,1),X(L,J,N,2),X(L,J,N,3),
     >                       X(L,J,K,1),X(L,J,K,2),X(L,J,K,3))
               VP4 = VOLPYM (X(I,M,K,1),X(I,M,K,2),X(I,M,K,3),
     >                       X(I,M,N,1),X(I,M,N,2),X(I,M,N,3),
     >                       X(L,M,N,1),X(L,M,N,2),X(L,M,N,3),
     >                       X(L,M,K,1),X(L,M,K,2),X(L,M,K,3))
               VP5 =-VOLPYM (X(I,J,K,1),X(I,J,K,2),X(I,J,K,3),
     >                       X(L,J,K,1),X(L,J,K,2),X(L,J,K,3),
     >                       X(L,M,K,1),X(L,M,K,2),X(L,M,K,3),
     >                       X(I,M,K,1),X(I,M,K,2),X(I,M,K,3))
               VP6 = VOLPYM (X(I,J,N,1),X(I,J,N,2),X(I,J,N,3),
     >                       X(L,J,N,1),X(L,J,N,2),X(L,J,N,3),
     >                       X(L,M,N,1),X(L,M,N,2),X(L,M,N,3),
     >                       X(I,M,N,1),X(I,M,N,2),X(I,M,N,3))
               VOL(I,J,K) = (VP1 + VP2 + VP3 + VP4 + VP5 + VP6) * FACTOR
            END DO
         END DO
      END DO

      END SUBROUTINE CELLVOL

C**********************************************************************
C
      SUBROUTINE CHECK_VOLS (NAME, IL, JL, KL, HAND, X, IWRIT)
C
C     Calculate cell volumes; count volumes < 0. & print with least volume.
C
C**********************************************************************

      IMPLICIT NONE

C     Arguments:

      CHARACTER, INTENT (IN) :: NAME * (*)
      INTEGER,   INTENT (IN) :: IL, JL, KL, IWRIT
      REAL,      INTENT (IN) :: HAND, X(IL,JL,KL,3)

C     Local constants:

      REAL, PARAMETER :: ZERO = 0.

C     Local variables:

      INTEGER I, J, K, L, IVMIN, JVMIN, KVMIN
      REAL    VOL(2:IL,2:JL,2:KL), VOLMIN
      REAL(KIND=4) CPU1, CPU2
C     Execution:

      CALL SECOND (CPU1)

      CALL CELLVOL (IL, JL, KL, X, VOL, HAND)

      VOLMIN = VOL(2,2,2)
      IVMIN = 2
      JVMIN = 2
      KVMIN = 2
      L = 0

      DO K = 2, KL
         DO J = 2, JL
            DO I = 2, IL
               IF (VOL(I,J,K) < ZERO) L = L + 1
               IF (VOL(I,J,K) < VOLMIN) THEN
                  VOLMIN = VOL(I,J,K)
                  IVMIN = I
                  JVMIN = J
                  KVMIN = K
               END IF
            END DO
         END DO
      END DO

      CALL SECOND (CPU2)
      CALL CPUTIME (CPU1, CPU2, NAME,
     >              'to check cell volumes', IWRIT)

      WRITE (IWRIT,'(/,A,3I5,A,1P,E15.4,/,A,I10)') ' Cell ',
     >   IVMIN, JVMIN, KVMIN, '  has least volume:', VOLMIN,
     >   ' Number of cells with negative volume:', L

      END SUBROUTINE CHECK_VOLS

C***********************************************************************
C
      SUBROUTINE CHORDS_FACES ()
C
C     Calculate unnormalized arc lengths for the radial lines of all
C     subblock faces, as needed for radial redistribution of faces.
C
C***********************************************************************

C     Global variables:

      USE LIM
      USE XYZ
      USE XYZ0

      IMPLICIT NONE

C     Local constants:

      REAL, PARAMETER :: ZERO = 0.

C     Local variables:

      INTEGER I, J, K, M

C     Local function:

      REAL    DELJ

      DELJ (I, J, K) = SQRT ((X(I,J,K,1) - X(I,J-1,K,1)) ** 2 +
     >                       (X(I,J,K,2) - X(I,J-1,K,2)) ** 2 +
     >                       (X(I,J,K,3) - X(I,J-1,K,3)) ** 2)

C     Execution:

      DO M = 1, NKREGRID
         K = KREGRID(M) ! 1, [MCRANK(*),] KTIP, KL
         DO I = 1, IL
            S0(I,1,K,2) = ZERO
            DO J = 2, JL
               S0(I,J,K,2) = S0(I,J-1,K,2) + DELJ (I, J, K)
            END DO
         END DO
      END DO

      DO M = 1, 5
         I = IFACE(M) ! 1, ITL, ILE, ITU, IL
         DO K = 2, KL - 1
            S0(I,1,K,2) = ZERO
            DO J = 2, JL
               S0(I,J,K,2) = S0(I,J-1,K,2) + DELJ (I, J, K)
            END DO
         END DO
      END DO

      END SUBROUTINE CHORDS_FACES

C***********************************************************************
C
      SUBROUTINE CHORDS_J1 ()
C
C     Calculate unnormalized arc lengths for the non-radial lines of the
C     J = 1 surface, as needed for blending radial face distributions.
C
C***********************************************************************

C     Global variables:

      USE LIM
      USE XYZ
      USE XYZ0

      IMPLICIT NONE

C     Local constants:

      REAL, PARAMETER :: ZERO = 0.

C     Local variables:

      INTEGER I, K

C     Execution:

      DO K = 1, KL
         S0(1,1,K,1) = ZERO
         DO I = 2, IL
            S0(I,1,K,1) = S0(I-1,1,K,1) + SQRT (
     >                      (X(I,1,K,1) - X(I-1,1,K,1)) ** 2 +
     >                      (X(I,1,K,2) - X(I-1,1,K,2)) ** 2 +
     >                      (X(I,1,K,3) - X(I-1,1,K,3)) ** 2)
         END DO
      END DO

      DO I = 1, IL
         S0(I,1,1,3) = ZERO
         DO K = 2, KL
            S0(I,1,K,3) = S0(I,1,K-1,3) + SQRT (
     >                      (X(I,1,K,1) - X(I,1,K-1,1)) ** 2 +
     >                      (X(I,1,K,2) - X(I,1,K-1,2)) ** 2 +
     >                      (X(I,1,K,3) - X(I,1,K-1,3)) ** 2)
         END DO
      END DO

      END SUBROUTINE CHORDS_J1

C***********************************************************************
C
      SUBROUTINE CHORDS_VOLUME ()
C
C     Calculate unnormalized arc lengths for the radial lines of the
C     whole grid, as needed for radial redistribution of subblock volumes.
C
C***********************************************************************

C     Global variables:

      USE LIM
      USE XYZ
      USE XYZ0

      IMPLICIT NONE

C     Local constants:

      REAL, PARAMETER :: ZERO = 0.

C     Local variables:

      INTEGER I, J, K

C     Execution:

      DO K = 1, KL
         DO I = 1, IL
            S0(I,1,K,2) = ZERO
            DO J = 2, JL
               S0(I,J,K,2) = S0(I,J-1,K,2) + SQRT (
     >                         (X(I,J,K,1) - X(I,J-1,K,1)) ** 2 +
     >                         (X(I,J,K,2) - X(I,J-1,K,2)) ** 2 +
     >                         (X(I,J,K,3) - X(I,J-1,K,3)) ** 2)
            END DO
         END DO
      END DO

      END SUBROUTINE CHORDS_VOLUME

C***********************************************************************
C
      SUBROUTINE CROWN ()
C
C     CROWN sets up the computational points on the fuselage crown/keel
C     line.  These are no longer influenced by the outer boundary pts.
C     This version uses D1NOSE (2-sided instead of 1-sided from the nose
C     to opposite the trailing edge).
C
C***********************************************************************

C     Global variables:

      USE CONSTS
      USE FUSEL
      USE INDEX
      USE LIM
      USE LUNS
      USE OUTR
      USE SPACING
      USE WNGSRF

      IMPLICIT NONE

C     Local constants:

      REAL,      PARAMETER :: ZERO = 0., TOLER = 1.E-7 ! PLXCUT tolerance
      LOGICAL,   PARAMETER :: NEW = .TRUE. ! New data for each call to LCSFIT
      CHARACTER, PARAMETER :: LOOSE*1 = 'B'

C     Local variables:

      INTEGER I, I1, I2, IER, NFORE
      REAL    TGRID(IL), XGRID(IL), YGRID(IL),
     >        TGEOM(MXIBODY), YGEOM(MXIBODY),
     >        BETA, D1, D2, R, TEDGE, TOL, TOTGEOM, XEDGE
      LOGICAL FIRST

C     Execution:

      YCROWN(1)   = YLBACK(JCR,1)
      YCROWN(IBL) = YF(1,NFSTN)
      YCROWN(IBU) = YF(NJBODY,NFSTN)
      YCROWN(IL)  = YUBACK(JCR,1)

C     Reuse VSHEET's X distribution aft of the trailing edge:

      DO I = ITU, IL
         XCROWN(I)          = XW(I,1)
         XCROWN(IL + 1 - I) = XW(I,1)
      END DO

C     Linear interpolation of Y aft of tail:

      R = (YCROWN(1) - YCROWN(IBL)) / (XCROWN(1) - XCROWN(IBL))
      DO I = 2, IBL - 1
         YCROWN(I) = YCROWN(IBL) + (XCROWN(I) - XCROWN(IBL)) * R
      END DO

      R = (YCROWN(IL) - YCROWN(IBU)) / (XCROWN(IL) - XCROWN(IBU))
      DO I = IBU + 1, IL - 1
         YCROWN(I) = YCROWN(IBU) + (XCROWN(I) - XCROWN(IBU)) * R
      END DO

C     Conventional spline interpolation of Y between t.e. and tail:

      YGEOM(1:NFSTN) = YF(NJBODY,1:NFSTN)

      CALL LCSFIT (NFSTN, XF, YGEOM, NEW, LOOSE, IBU - ITU,
     >             XCROWN(ITU), YCROWN(ITU), YCROWN(ITU))

C     Arc-based distribution forward of the trailing edge:
C     We need T corresponding to the trailing edge X = XCROWN(ITU/L).

      CALL CHORDS2D (NFSTN, XF, YGEOM, .FALSE., TOTGEOM, TGEOM)

      XEDGE = XCROWN(ITU)
      I1    = NFSTN / 2
      TOL   = MAX (5.* EPSMCH, TOLER)

      CALL PLXCUT (NFSTN, XF, YGEOM, TGEOM, I1, .FALSE., XEDGE, TOL,
     >             TEDGE, -IWRIT, IER)

      IF (IER > 1) THEN
         WRITE (IWRIT,*) 'CROWN 1: I1, X, TEDGE = ', I1, XEDGE, TEDGE
         GO TO 900
      END IF

C     Vinokur T distribution from the nose to opposite the trailing edge:

      D1 = D1NOSE
      D2 = XCROWN(ITU+1) - XEDGE ! T ~ X here
      NFORE = ITU - ILE + 2      ! Overlap by 1

      CALL HTDIS4 (.TRUE., ZERO, TEDGE + D2, D1, D2, NFORE, TGRID,
     >             -IWRIT, IER)

      IF (IER /= 0 .AND. IER /= 3) THEN
         WRITE (IWRIT,*) 'CROWN (a): NFORE, TEDGE, D1, D2 = ',
     >                               NFORE, TEDGE, D1, D2
         GO TO 900
      END IF

C     Interpolate X, Y at these Ts:

      CALL LCSFIT (NFSTN, TGEOM, XF, NEW, LOOSE, NFORE-2, TGRID,
     >             XCROWN(ILE), XCROWN(ILE))

      CALL LCSFIT (NFSTN, TGEOM, YGEOM, NEW, LOOSE, NFORE-2, TGRID,
     >             YCROWN(ILE), YCROWN(ILE))

C     Repeat for the keel line:

      DO I = 1, NFSTN
         YGEOM(I) = YF(1,I)
      END DO

      DO I = ITL, IBL + 1, -1
         CALL LCSFIT (NFSTN, XF, YGEOM, NEW, LOOSE, 1, XCROWN(I),
     >                YCROWN(I), YCROWN(I))
      END DO

      CALL CHORDS2D (NFSTN, XF, YGEOM, .FALSE., TOTGEOM, TGEOM)

      CALL PLXCUT (NFSTN, XF, YGEOM, TGEOM, I1, .FALSE., XEDGE, TOL,
     >             TEDGE, -IWRIT, IER)

      IF (IER > 1) THEN
         WRITE (IWRIT,*) 'CROWN 2: I1, X, TEDGE = ', I1, XEDGE, TEDGE
         GO TO 900
      END IF

      CALL HTDIS4 (.TRUE., ZERO, TEDGE + D2, D1, D2, NFORE, TGRID,
     >             -IWRIT, IER)

      IF (IER /= 0 .AND. IER /= 3) THEN
         WRITE (IWRIT,*) 'CROWN (b): NFORE, TEDGE, D1, D2 = ',
     >                               NFORE, TEDGE, D1, D2
         GO TO 900
      END IF

      CALL LCSFIT (NFSTN, TGEOM, XF, NEW, LOOSE, NFORE-2, TGRID,
     >             XGRID, XGRID)

      CALL LCSFIT (NFSTN, TGEOM, YGEOM, NEW, LOOSE, NFORE-2, TGRID,
     >             YGRID, YGRID)

      XCROWN(ITL+1:ILE-1) = XGRID(NFORE-2:2:-1)
      YCROWN(ITL+1:ILE-1) = YGRID(NFORE-2:2:-1)

      RETURN

  900 STOP

      END SUBROUTINE CROWN

C***********************************************************************
C
      SUBROUTINE EXCHANGE (MODE1, MODE2)
C
C     Swap I1 & I2 C-boundary controls for 2- or 3-D elliptic smoothing.
C
C***********************************************************************

      CHARACTER MODE1*(*), MODE2*(*)

      MODE2      = MODE1
      MODE2(1:1) = MODE1(2:2)
      MODE2(2:2) = MODE1(1:1)

      END SUBROUTINE EXCHANGE

C***********************************************************************
C
      SUBROUTINE FIXBGRID (IDIM, JDIM, IBL, ITL, ILE, ITU, IBU, JCR,
     >                     DELTA1, DELTA2, DNOSE, DW, DC, NBLAYER,
     >                     RBLAYER, X, Y, Z, LUNERR)
C
C     FIXBGRID redistributes the radial lines of the body surface grid
C     using Vinokur distributions along the arcs.
C
C     DELTA1 = D1 at the leading edge;
C     DELTA2 = D1 at the trailing edge;
C     DNOSE  = D2 at the nose;
C     DW     = D1 at the tail ("DWake")
C     DC     = D2 at the tail ("DCrown/keel")
C
C     09/19/96  DAS  Introduced DELTA2 & DI(*) - now D1I(*).
C     04/24/97   "   Added DNOSE, DW, DC, & D2I(*).
C     05/02/97   "   D2I(*) needs limiting near low or high wings.
C     07/01/98   "   Introduced NBLAYER & RBLAYER and BLGRID.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   IDIM, JDIM, IBL, ITL, ILE, ITU, IBU, JCR, NBLAYER, LUNERR
      REAL, INTENT (IN) ::
     >   DELTA1, DELTA2, DNOSE, DW, DC, RBLAYER
      REAL, INTENT (INOUT) ::
     >   X(IDIM,JDIM), Y(IDIM,JDIM), Z(IDIM,JDIM)

C     Local constants:

      REAL,      PARAMETER :: ONE = 1., ZERO = 0., SUPPRESS = -999.
      LOGICAL,   PARAMETER :: CLOSED = .FALSE.
      CHARACTER, PARAMETER :: METHOD * 1 = 'B'

C     Local variables:

      INTEGER
     >   I, IER, J, JEVAL, NJ
      REAL
     >   BETA, D1, D2, DENOM, DIFF, DUNIFORM, POWER, R, TTOTAL
      REAL, DIMENSION (IDIM) ::
     >   D1I, D2I, TJ
      REAL, DIMENSION (JCR) ::
     >   TGRID, XJ, YJ, ZJ
      LOGICAL
     >   NEW

C     Execution:

C     Vary initial increment smoothly from leading edge to end of body.
C     Reuse existing J-direction variables; TJ(*) has to be bigger though.

      NJ = IBU - IBL + 1

C     Normalized arc lengths along the J = 1 line:

      CALL CHORDS3D (NJ, X(IBL,1), Y(IBL,1), Z(IBL,1), .TRUE., TTOTAL,
     >               TJ(IBL))

C     Set up the initial radial increments as functions of I.
C     Navier-Stokes cases force nonlinear variation from TE to tail.

      IF (DW <= DELTA2 * 2.) THEN
         POWER = 1.
      ELSE IF (DW <= DELTA2 * 5.) THEN
         POWER = 2.
      ELSE
         POWER = 3.
      END IF

      DENOM = ONE / (TJ(IBL) - TJ(ITL))
      DIFF  = DW - DELTA2

      DO I = IBL, ITL - 1
         R = (TJ(I) - TJ(ITL)) * DENOM
         D1I(I) = DELTA2 + DIFF * R ** POWER
      END DO

      DENOM = ONE / (TJ(ILE) - TJ(ITL))
      DO I = ITL, ILE ! Linear
         R = (TJ(I) - TJ(ITL)) * DENOM
         D1I(I) = (ONE - R) * DELTA2 + R * DELTA1
      END DO

      DENOM = ONE / (TJ(ITU) - TJ(ILE))
      DO I = ILE, ITU ! Linear
         R = (TJ(I) - TJ(ILE)) * DENOM
         D1I(I) = (ONE - R) * DELTA1 + R * DELTA2
      END DO

      DENOM = ONE / (TJ(IBU) - TJ(ITU))
      DO I = ITU + 1, IBU ! Polynomial
         R = (TJ(I) - TJ(ITU)) * DENOM
         D1I(I) = DELTA2 + DIFF * R ** POWER
      END DO

C     Variation of last radial increments from DNOSE to DC(rown) along the
C     keel is nominally linear, but may be adjusted below:

      DENOM = ONE / (TJ(ILE) - TJ(IBL))
      DO I = IBL, ILE
         R = (TJ(I) - TJ(IBL)) * DENOM
         D2I(I) = (ONE - R) * DC + R * DNOSE
      END DO

C     Likewise for the last radial increments along the crown:

      DENOM = ONE / (TJ(IBU) - TJ(ILE))
      DO I = ILE, IBU
         R = (TJ(I) - TJ(ILE)) * DENOM
         D2I(I) = (ONE - R) * DNOSE + R * DC
      END DO

C     For each radial line, impose an arc-length-based Vinokur distribution:

      NJ = JCR
      DUNIFORM = 0.98 / REAL (NJ - 1) ! I.e. a bit less to protect EXPDIS4

      DO I = IBL, IBU

C        Gather the current J line as needed for parametric interpolation:

         XJ(1:NJ) = X(I,1:NJ)
         YJ(1:NJ) = Y(I,1:NJ)
         ZJ(1:NJ) = Z(I,1:NJ)

C        Unnormalized arc length increments and total length:

         CALL CHORDS3D (NJ, XJ, YJ, ZJ, .FALSE., TTOTAL, TJ)

C        Check D2 against the one-sided stretching value:

         D1 = MIN (D1I(I), DUNIFORM * TTOTAL)

         CALL EXPDIS4 (3, ZERO, TTOTAL, D1, NJ, TGRID, BETA, -LUNERR)

         D2 = MIN (D2I(I), TTOTAL - TGRID(NJ-1))

C        Unnormalized geometric + Vinokur distribution:

         CALL BLGRID (NJ, D1, D2, NBLAYER, RBLAYER, TGRID, LUNERR, IER)

         IF (IER /= 0 .AND. IER /= 3) THEN
            WRITE (LUNERR,*) 'FIXBGRID: BLGRID trouble. IER: ', IER
            WRITE (LUNERR,*) 'I, NJ, D1, D2: ', I, NJ, D1, D2
            STOP
         END IF

C        Interpolate X, Y, Z at these arc lengths:

         NEW = .TRUE.
         JEVAL = 1

         DO J = 2, NJ - 1

            CALL PLSCRV3D (NJ, XJ, YJ, ZJ, TJ, METHOD, NEW, CLOSED,
     >                     TGRID(J), JEVAL, X(I,J), Y(I,J), Z(I,J),
     >                     SUPPRESS) ! Derivatives
            NEW = .FALSE.
         END DO

      END DO

      END SUBROUTINE FIXBGRID

C*******************************************************************************
C
      SUBROUTINE FIXLEDGE (IL, JL, KL, ILE, JLAST, K1, K2, PLANEILE,
     >                     T, XYZ)
C
C     FIXLEDGE (fix leading edge) modularizes redistributing some chordwise
C     points forward of the wing following use of TFIQ3D on the plane I = ILE,
C     in a way that accounts for non-trivial leading edge deflections.
C
C     08/26/97  D.Saunders  Initial implementation, prompted by slats, following
C                           no success with the HERMQ3D elaboration of TFIQ3D.
C     09/30/97      "       Vary last J from JLAST/2 to JLAST with K to follow
C                           a body nose better; also, better to leave Z alone.
C     06/29/98      "       Low-wing HSCT gave trouble near the body because of
C                           3-D effects inboard of artificial crank station.
C                           Therefore, take trivial LE slopes into account.
C     06/30/98      "       Wrinkles still appear.  Therefore, don't touch Ks
C                           for which the wing grid section is not vertical
C                           because the leading edge angle is ill-defined.
C                           This means specifying a crank inboard of any slat.
C*******************************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN)    :: IL, JL, KL, ILE, JLAST, K1, K2
      REAL,    INTENT (INOUT) :: PLANEILE(JL,KL,3)
      REAL,    INTENT (OUT)   :: T(JLAST)
      REAL,    INTENT (IN)    :: XYZ(IL,JL,KL,3)

C     Local constants:

      REAL,      PARAMETER :: HALF = 0.5, ONE = 1.,
     >                        ANGLESIG = 0.35 ! radians, ~20 degrees
      INTEGER,   PARAMETER :: NLINE = 8
      LOGICAL,   PARAMETER :: NEW = .TRUE.
      CHARACTER, PARAMETER :: LOOSE * 1 = 'B'

C     Local variables:

      INTEGER    I, IER, J, JM, JU, K
      REAL       ANGLE, DJ, DX2, DY2, EPS, R, RJ, TOTAL1, TOTAL2, YJ,
     >           XLINE(NLINE), YLINE(NLINE), ZLINE(NLINE), TLINE(NLINE)

C     Execution:

      JM = JLAST / 2 ! Near a body, follow it more closely
      DJ = REAL (JLAST - JM) / REAL (K2 - K1)
      RJ = REAL (JM)
      EPS = 10.* EPSILON (ONE)

      DO K = K1, K2

C        Set up an 8-point curve with the right slope at the leading edge:

         DO I = 1, 3
            J = 3 - I
            XLINE(I) = (XYZ(ILE+J,1,K,1) + XYZ(ILE-J,1,K,1)) * HALF
            YLINE(I) = (XYZ(ILE+J,1,K,2) + XYZ(ILE-J,1,K,2)) * HALF
            ZLINE(I) = (XYZ(ILE+J,1,K,3) + XYZ(ILE-J,1,K,3)) * HALF
         END DO

         DX2 = XLINE(1) - XLINE(3)
         DY2 = YLINE(1) - YLINE(3)
         XLINE(4) = XLINE(3) - DX2
         YLINE(4) = YLINE(3) - DY2
         ZLINE(4) = ZLINE(3)

         RJ = RJ + DJ
         JM = MIN (NINT (RJ), JLAST)

C        Wing grid sections affected by a wing/body intersection do not have
C        well-defined leading edge angles.  Leave such Ks alone:

         J = ILE / 4
         IF (ABS (XYZ(ILE+J,1,K,3) - XYZ(ILE-J,1,K,3)) > EPS) CYCLE

C        Reduce the extent of the adjustment for smaller leading edge slopes,
C        because the TFI/ELLIPQ3D result is normally preferable.
C        The upper J index, JU, is reduced quadratically from the nominal value,
C        JM, to as low as NLINE as the slope varies from "significant"
C        (ANGLESIG) to insignificant (zero).

         ANGLE = ABS (ATAN2 (DY2, DX2))

         IF (ANGLE >= ANGLESIG) THEN
            JU = JM
         ELSE
            YJ = REAL (NLINE-JM) * (ANGLE / ANGLESIG-ONE)**2 + REAL (JM)
            JU = NINT (YJ)
         END IF

         IF (JU <= NLINE) CYCLE

         J = JU - 4

         DO I = 5, NLINE
            J = J + 1
            XLINE(I) = PLANEILE(J,K,1)
            YLINE(I) = PLANEILE(J,K,2)
            ZLINE(I) = PLANEILE(J,K,3)
         END DO

         CALL CHORDS3D (NLINE, XLINE, YLINE, ZLINE, .FALSE., TOTAL1,
     >                  TLINE)

C        Unnormalized arc-length distribution of the interim line:

         CALL CHORDS3D (JU, PLANEILE(1,K,1), PLANEILE(1,K,2),
     >                  PLANEILE(1,K,3), .FALSE., TOTAL2, T)

C        Adjust the distribution for the revised total length:

         R = (TOTAL1 - TLINE(3)) / TOTAL2
         DO J = 1, JU - 1
            T(J) = TLINE(3) + R * T(J)
         END DO
         T(JU) = TOTAL1

C        Interpolate X & Y as functions of the redistributed T:

         CALL LCSFIT (NLINE, TLINE, XLINE, NEW, LOOSE, JU, T,
     >                PLANEILE(1,K,1), PLANEILE(1,K,1))
         CALL LCSFIT (NLINE, TLINE, YLINE, NEW, LOOSE, JU, T,
     >                PLANEILE(1,K,2), PLANEILE(1,K,2))

C        Leave Z alone to reduce possible trouble near the nose.

      END DO

      END SUBROUTINE FIXLEDGE

C***********************************************************************
C
      SUBROUTINE FIXWAKE (IL, KL, ITL, ITU, ILAST, K1, K2, T,
     >                    XW, YW, ZW)
C
C     FIXWAKE adjusts streamwise lines in the wake to account for flap
C     deflections.
C
C     09/17/97  DAS  Modularization of VSHEET code for multiple use.
C     09/30/97   "   Vary last I from ILAST/2 to ILAST with K to follow
C                    aft body better; also, better to leave Z alone.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) :: IL, KL, ITL, ITU, ILAST, K1, K2
      REAL, INTENT (OUT)   :: T(ILAST)
      REAL, INTENT (INOUT) :: XW(IL,KL), YW(IL,KL)
      REAL, INTENT (IN)    :: ZW(IL,KL)

C     Local constants:

      REAL,      PARAMETER :: HALF = 0.5
      INTEGER,   PARAMETER :: NLINE = 7
      LOGICAL,   PARAMETER :: NEW = .TRUE.
      CHARACTER, PARAMETER :: LOOSE*1 = 'B'

C     Local variables:

      INTEGER  I, IM, J, K, NI
      REAL     DI, DX2, DY2, DZ2, R, RI, TOTAL1, TOTAL2,
     >         TLINE(NLINE), XLINE(NLINE), YLINE(NLINE), ZLINE(NLINE)

C     Execution:

      IM = (ITU + ILAST) / 2
      DI = REAL (ILAST - IM) / REAL (K2 - K1)
      RI = REAL (IM)

      DO K = K1, K2

         DO I = 1, 3
            J = 3 - I
            XLINE(I) = (XW(ITL+J,K) + XW(ITU-J,K)) * HALF
            YLINE(I) = (YW(ITL+J,K) + YW(ITU-J,K)) * HALF
            ZLINE(I) = (ZW(ITL+J,K) + ZW(ITU-J,K)) * HALF
         END DO

         DX2 = XLINE(3) - XLINE(1)
         DY2 = YLINE(3) - YLINE(1)
         DZ2 = ZLINE(3) - ZLINE(1)

         XLINE(4) = XLINE(3) + DX2
         YLINE(4) = YLINE(3) + DY2
         ZLINE(4) = ZLINE(3) + DZ2

         RI = RI + DI
         IM = MIN (NINT (RI), ILAST)
         NI = IM - ITU + 1
         J  = IM - 3

         DO I = 5, NLINE
            J = J + 1
            XLINE(I) = XW(J,K)
            YLINE(I) = YW(J,K)
            ZLINE(I) = ZW(J,K)
         END DO

         CALL CHORDS3D (NLINE, XLINE, YLINE, ZLINE, .FALSE., TOTAL1,
     >                  TLINE)

C        Unnormalized arc lengths of the interim line:

         CALL CHORDS3D (NI, XW(ITU,K), YW(ITU,K), ZW(ITU,K), .FALSE.,
     >                  TOTAL2, T)

C        Adjust the interim distribution for the revised total length:

         R = (TOTAL1 - TLINE(3)) / TOTAL2
         DO I = 1, NI - 1
            T(I) = TLINE(3) + R * T(I)
         END DO
         T(NI) = TOTAL1

C        Interpolate X & Y as functions of the redistributed T:

         CALL LCSFIT (NLINE, TLINE, XLINE, NEW, LOOSE, NI, T,
     >                XW(ITU,K), XW(ITU,K))
         CALL LCSFIT (NLINE, TLINE, YLINE, NEW, LOOSE, NI, T,
     >                YW(ITU,K), YW(ITU,K))

C        Leave Z alone to reduce possible trouble near the aft body.

      END DO

      END SUBROUTINE FIXWAKE

C***********************************************************************
C
      SUBROUTINE GRID_ELLIP3D (NCALL)
C
C     Smooth all subblock interiors via 3D elliptic smoothing.
C
C***********************************************************************

C     Global variables:

      USE ELLIPTIC
      USE LIM
      USE LOGIC
      USE LUNS
      USE OPTIONS
      USE XYZ
      USE XYZ0

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) :: NCALL

C     Local constants:

      CHARACTER, PARAMETER :: E3D*9 = ' ELLIP3D '

C     Local variables:

      INTEGER   L
      REAL(KIND=4) CPU1, CPU2
      CHARACTER BGMODE*4, BGPH*4, BGPS*4, BGOM*4, FG3D*6, FGMODE*4

C     Execution:

      IF (NCALL == -3) CALL SECOND (CPU1)

      IF (NNTIP == 0) THEN

         L = 1

         IF (PRINT3D)
     >      WRITE (IWRIT,'(/,2A)') E3D, 'below & aft of wing:'

         CALL ELLIP3D (IL, JL, KL, 1, ITL, 1, JL, 1, KTIP,
     >                 X, S0, PRINT3D,
     >                 ITMAX3D(L), CONV3D(L), DMAX3D, OMG3D(L),
     >                 BGPHI(L), BGPSI(L), BGOMG(L), FGMODE3D(L),
     >                 FGLIM3D, URFG3D, EXPI13D(L), EXPI23D(L),
     >                 EXPJ13D(L), EXPJ23D(L), EXPK13D(L), EXPK23D(L),
     >                 NNTIP, ILE, KTIP)

         L = 2

         IF (PRINT3D) WRITE (IWRIT,'(/,2A)') E3D, 'below wing:'

         CALL ELLIP3D (IL, JL, KL, ITL, ILE, 1, JL, 1, KTIP,
     >                 X, S0, PRINT3D,
     >                 ITMAX3D(L), CONV3D(L), DMAX3D, OMG3D(L),
     >                 BGPHI(L), BGPSI(L), BGOMG(L), FGMODE3D(L),
     >                 FGLIM3D, URFG3D, EXPI13D(L), EXPI23D(L),
     >                 EXPJ13D(L), EXPJ23D(L), EXPK13D(L), EXPK23D(L),
     >                 NNTIP, ILE, KTIP)

         CALL EXCHANGE (BGPHI(L), BGPH)
         CALL EXCHANGE (BGPSI(L), BGPS)
         CALL EXCHANGE (BGOMG(L), BGOM)
         CALL EXCHANGE (FGMODE3D(L), FG3D)

         IF (PRINT3D) WRITE (IWRIT,'(/,2A)') E3D, 'above wing:'

         CALL ELLIP3D (IL, JL, KL, ILE, ITU, 1, JL, 1, KTIP,
     >                 X, S0, PRINT3D,
     >                 ITMAX3D(L), CONV3D(L), DMAX3D, OMG3D(L),
     >                 BGPH, BGPS, BGOM, FG3D,
     >                 FGLIM3D, URFG3D, EXPI23D(L), EXPI13D(L),
     >                 EXPJ13D(L), EXPJ23D(L), EXPK13D(L), EXPK23D(L),
     >                 NNTIP, ILE, KTIP)

         L = 1
         CALL EXCHANGE (BGPHI(L), BGPH)
         CALL EXCHANGE (BGPSI(L), BGPS)
         CALL EXCHANGE (BGOMG(L), BGOM)
         CALL EXCHANGE (FGMODE3D(L), FG3D)

         IF (PRINT3D)
     >      WRITE (IWRIT,'(/,2A)') E3D, 'above & aft of wing:'

         CALL ELLIP3D (IL, JL, KL, ITU, IL, 1, JL, 1, KTIP,
     >                 X, S0, PRINT3D,
     >                 ITMAX3D(L), CONV3D(L), DMAX3D, OMG3D(L),
     >                 BGPH, BGPS, BGOM, FG3D,
     >                 FGLIM3D, URFG3D, EXPI23D(L), EXPI13D(L),
     >                 EXPJ13D(L), EXPJ23D(L), EXPK13D(L), EXPK23D(L),
     >                 NNTIP, ILE, KTIP)

C        Beyond the tip:

         L = 3

         IF (PRINT3D)
     >      WRITE (IWRIT,'(/,2A)') E3D, 'beyond tip, lower/aft:'

         CALL ELLIP3D (IL, JL, KL, 1, ITL, 1, JL, KTIP, KL,
     >                 X, S0, PRINT3D,
     >                 ITMAX3D(L), CONV3D(L), DMAX3D, OMG3D(L),
     >                 BGPHI(L), BGPSI(L), BGOMG(L), FGMODE3D(L),
     >                 FGLIM3D, URFG3D, EXPI13D(L), EXPI23D(L),
     >                 EXPJ13D(L), EXPJ23D(L), EXPK13D(L), EXPK23D(L),
     >                 NNTIP, ILE, KTIP)

         L = 4

         IF (PRINT3D)
     >      WRITE (IWRIT,'(/,2A)') E3D, 'beyond tip, lower/fore:'

         CALL ELLIP3D (IL, JL, KL, ITL, ILE, 1, JL, KTIP, KL,
     >                 X, S0, PRINT3D,
     >                 ITMAX3D(L), CONV3D(L), DMAX3D, OMG3D(L),
     >                 BGPHI(L), BGPSI(L), BGOMG(L), FGMODE3D(L),
     >                 FGLIM3D, URFG3D, EXPI13D(L), EXPI23D(L),
     >                 EXPJ13D(L), EXPJ23D(L), EXPK13D(L), EXPK23D(L),
     >                 NNTIP, ILE, KTIP)

         CALL EXCHANGE (BGPHI(L), BGPH)
         CALL EXCHANGE (BGPSI(L), BGPS)
         CALL EXCHANGE (BGOMG(L), BGOM)
         CALL EXCHANGE (FGMODE3D(L), FG3D)

         IF (PRINT3D)
     >      WRITE (IWRIT,'(/,2A)') E3D, 'beyond tip, upper/fore:'

         CALL ELLIP3D (IL, JL, KL, ILE, ITU, 1, JL, KTIP, KL,
     >                 X, S0, PRINT3D,
     >                 ITMAX3D(L), CONV3D(L), DMAX3D, OMG3D(L),
     >                 BGPH, BGPS, BGOM, FG3D,
     >                 FGLIM3D, URFG3D, EXPI23D(L), EXPI13D(L),
     >                 EXPJ13D(L), EXPJ23D(L), EXPK13D(L), EXPK23D(L),
     >                 NNTIP, ILE, KTIP)

         L = 3
         CALL EXCHANGE (BGPHI(L), BGPH)
         CALL EXCHANGE (BGPSI(L), BGPS)
         CALL EXCHANGE (BGOMG(L), BGOM)
         CALL EXCHANGE (FGMODE3D(L), FG3D)

         IF (PRINT3D)
     >      WRITE (IWRIT,'(/,2A)') E3D, 'beyond tip, upper/aft:'

         CALL ELLIP3D (IL, JL, KL, ITU, IL, 1, JL, KTIP, KL,
     >                 X, S0, PRINT3D,
     >                 ITMAX3D(L), CONV3D(L), DMAX3D, OMG3D(L),
     >                 BGPH, BGPS, BGOM, FG3D,
     >                 FGLIM3D, URFG3D, EXPI23D(L), EXPI13D(L),
     >                 EXPJ13D(L), EXPJ23D(L), EXPK13D(L), EXPK23D(L),
     >                 NNTIP, ILE, KTIP)

      ELSE ! Smooth the C-H + C-O volume grid:

         L = 5

         IF (PRINT3D)
     >      WRITE (IWRIT,'(/,2A)') E3D, 'below wing:'

         CALL ELLIP3D (IL, JL, KL, 1, ILE, 1, JL, 1, KL,
     >                 X, S0, PRINT3D,
     >                 ITMAX3D(L), CONV3D(L), DMAX3D, OMG3D(L),
     >                 BGPHI(L), BGPSI(L), BGOMG(L), FGMODE3D(L),
     >                 FGLIM3D, URFG3D, EXPI13D(L), EXPI23D(L),
     >                 EXPJ13D(L), EXPJ23D(L), EXPK13D(L), EXPK23D(L),
     >                 NNTIP, ILE, KTIP)

         CALL EXCHANGE (BGPHI(L), BGPH)
         CALL EXCHANGE (BGPSI(L), BGPS)
         CALL EXCHANGE (BGOMG(L), BGOM)
         CALL EXCHANGE (FGMODE3D(L), FG3D)

         IF (PRINT3D)
     >      WRITE (IWRIT,'(/,2A)') E3D, 'above wing:'

         CALL ELLIP3D (IL, JL, KL, ILE, IL, 1, JL, 1, KL,
     >                 X, S0, PRINT3D,
     >                 ITMAX3D(L), CONV3D(L), DMAX3D, OMG3D(L),
     >                 BGPH, BGPS, BGOM, FG3D,
     >                 FGLIM3D, URFG3D, EXPI23D(L), EXPI13D(L),
     >                 EXPJ13D(L), EXPJ23D(L), EXPK13D(L), EXPK23D(L),
     >                 NNTIP, ILE, KTIP)
      END IF

      IF (NCALL == -3) THEN
         CALL SECOND (CPU2)
         CALL CPUTIME (CPU1, CPU2, 'GRID_ELLIP3D',
     >                 'to smooth the volume grid', IWRIT)
      END IF

      END SUBROUTINE GRID_ELLIP3D

C***********************************************************************
C
      SUBROUTINE GRID_FACES (NCALL, FULLGRID, FAIL)
C
C     Establish all subblock face grids, either from scratch or by warping.
C     If FULLGRID = F and REDISTRIBUTE = T, care must be taken to match the
C     initial full/indirect grid exactly for the case of no surface change,
C     by starting with Euler spacing first on the body surface and some
C     edges and regridding them.
C
C***********************************************************************

C     Global variables:

      USE DMAX
      USE ELLIPTIC
      USE FACES
      USE LIM
      USE LOGIC
      USE LUNS
      USE OPTIONS
      USE OUTR
      USE SPACING
      USE TFI
      USE WALL
      USE WNGSRF
      USE XYZ
      USE XYZ0

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN)  :: NCALL
      LOGICAL, INTENT (IN)  :: FULLGRID
      LOGICAL, INTENT (OUT) :: FAIL

C     Local constants:

      REAL,      PARAMETER :: HALF = 0.5, ONE = 1., ZERO = 0.
      CHARACTER, PARAMETER :: E2D*9 = ' ELLIP2D '

C     Local variables:

      INTEGER
     >   I, I1, I2, ILAST, ITERS, J, JLAST, K, K1, K2, L, M, N

      REAL
     >   PLANEI(JL,KL,3), YTEU(KL), XJ(JL), YJ(JL), ZJ(JL),
     >   COSTH, D1KTIP, R, SINTH, THETA, XSHIFT, ZWK

      REAL(KIND=4) CPU1, CPU2

      LOGICAL
     >   EULER, WARPNS

      CHARACTER
     >   BGMODE*4, FGMODE*4

C     Execution:

      IF (NCALL < 0) CALL SECOND (CPU1)


C     Calculate the wing/body intersection (if any); grid the wing surface:

      CALL WSURF (NCALL, FAIL)

      IF (FAIL) GO TO 999


C     Establish the radial grid controls, which can vary if REDISTRIBUTE = T.
C     Both sets are initialized once using CHORDM(*) and ZWLE(*) from WSURF:

      IF (NCALL == -3) THEN

         CALL RADIAL_INIT (1) ! 0 means use X(*,*,*,) instead (READGRID case)

      ELSE IF (.NOT. CHECKWARP) THEN ! Suppress excessive output

         PRINT2D = .FALSE.
         PRINT3D = .FALSE.

      END IF


C     Activate the correct set of radial controls:

      EULER  = FULLGRID .OR. .NOT. REDISTRIBUTE
      WARPNS = .NOT. EULER ! Surprisingly


C     ***************************************************************
C            Clarification of the radial control usage, which
C            is complicated by the grid perturbation option:
C
C     EULER = F means use the *NS set of controls AND we are warping,
C     so GRID_FACES should set up the edges of warped faces in a way
C     that matches how REGRID_FACES sets them up initially.  Also:
C     EULER = T can still mean N-S-type grids if REDISTRIBUTE = F and
C     the initial radial (*EU) and smoothing controls are appropriate.
C     ***************************************************************


      CALL RADIAL_SET (EULER)


C     Grid the outer boundaries:

      CALL OUTER (FULLGRID)


C     Set up (most of) the wake boundaries - along the body is done in Z0PLANE:

      CALL VSHEET ()


C     Fuselage crown/keel line grid?

      IF (.NOT. WINGONLY) CALL CROWN ()


C     Generate the grid at the symmetry plane (Z = 0), and on the fuselage
C     surface, if it exists.  Also, complete the wake boundary.

      CALL Z0PLANE (FULLGRID, FAIL)

      IF (FAIL) GO TO 999


C     The rest of GRID_FACES used to be the first half of VOLGRID.

      X(ITL:ITU,1,1:KL,1) = XW(ITL:ITU,1:KL) ! Wing surface (J = 1)
      X(ITL:ITU,1,1:KL,2) = YW(ITL:ITU,1:KL)
      X(ITL:ITU,1,1:KL,3) = ZW(ITL:ITU,1:KL)

      X(1:IL,1:JL,1,1)  = XWALL ! Symmetry plane, including any body
      X(1:IL,1:JL,1,2)  = YWALL
      X(1:IL,1:JL,1,3)  = ZWALL

      X(1, 1:JL,1:KL,1) = XLBACK ! Downstream face boundaries
      X(IL,1:JL,1:KL,1) = XUBACK
      X(1, 1:JL,1:KL,2) = YLBACK
      X(IL,1:JL,1:KL,2) = YUBACK
      X(1, 1:JL,1:KL,3) = ZLBACK
      X(IL,1:JL,1:KL,3) = ZUBACK

C     Upper part of the wake sheet:

      IF (WINGONLY) THEN
         ILAST = ITU + 2 * (IL - ITU) / 3 ! Extent of FIXWAKE redistribution
      ELSE
         ILAST = IBU
      END IF

      IF (FULLGRID) THEN

C        Preserve any trailing edge thickness:

         DO K = 1, KTIP
            YTEU(K)   =  YW(ITU,K)
            YW(ITU,K) = (YW(ITL,K) + YW(ITU,K)) * HALF
         END DO

         IF (NNTIP == 0) THEN
            CALL TFIQ3D (IL, ITU, IL, 1,  KTIP, XW, YW, ZW, TFIWORK)
            CALL TFIQ3D (IL, ITU, IL, KTIP, KL, XW, YW, ZW, TFIWORK)
         ELSE
            CALL TFIQ3D (IL, ITU, IL, 1,    KL, XW, YW, ZW, TFIWORK)
         END IF

C        Recover any trailing edge thickness:

         YW(ITU,1:KTIP) = YTEU(1:KTIP)

C        Account for any trailing edge deflections:

         CALL FIXWAKE (IL, KL, ITL, ITU, ILAST, 2, KTIP - 1, TFIWORK,
     >                 XW, YW, ZW)

C        Smooth the wake sheet:

         DO K = 1, KTIP
            YW(ITU,K) = (YW(ITL,K) + YW(ITU,K)) * HALF
         END DO

         L = 3

         IF (PRINT2D)
     >      WRITE (IWRIT,'(/,A)') ' ELLIPQ3D on the wake sheet:'

         CALL ELLIPQ3D (IL, KL, ITU, IL, 1, KTIP, XW, ZW, YW,
     >                  DFACEJ,       PRINT2D,
     >                  BGMODE2D(L),  FGMODE2D(L),  SPMODE(L),
     >                  ITMAX2D(L),   ITFLOAT,      ITFREEZE, JLAYER,
     >                  CONV2D(L),    CONVMIN(L),   DMAX2D, OMG2D(L),
     >                  POWERI(L),    POWERJ(L),    EXPI12D(L),
     >                  EXPI22D(L),   EXPJ12D(L),   EXPJ22D(L),
     >                  FGLIM2D,      URFG2D,       URFLOAT)

         YW(ITU,1:KTIP) = YTEU(1:KTIP)

      ELSE ! Perturb the upper wake sheet

C        Transfer remaining upper wake boundaries:

         I = ITU + 1
         X(I:IL,1,KTIP,1) = XW(I:IL,KTIP)
         X(I:IL,1,  KL,1) = XW(I:IL,KL)
         X(I:IL,1,KTIP,2) = YW(I:IL,KTIP)
         X(I:IL,1,  KL,2) = YW(I:IL,KL)
         X(I:IL,1,KTIP,3) = ZW(I:IL,KTIP)
         X(I:IL,1,  KL,3) = ZW(I:IL,KL)

         IF (NNTIP == 0) THEN

            DO M = 1, 2

               K1 = KFACE(M)
               K2 = KFACE(M+1)

               CALL WARPQ3D (1, IL, 1, JL, 1, KL, ITU, IL, 1, 1, K1, K2,
     >                       X0(1,1,1,1), X0(1,1,1,2), X0(1,1,1,3), S0,
     >                       DFACEI,     DFACEJ,     DFACEK,
     >                       X(1,1,1,1), X(1,1,1,2), X(1,1,1,3))
            END DO

         ELSE

            CALL WARPQ3D (1, IL, 1, JL, 1, KL, ITU, IL, 1, 1, 1, KL,
     >                    X0(1,1,1,1), X0(1,1,1,2), X0(1,1,1,3), S0,
     >                    DFACEI,     DFACEJ,     DFACEK,
     >                    X(1,1,1,1), X(1,1,1,2), X(1,1,1,3))
         END IF

         DO K = 2, KL
            XW(I:IL,K) = X(I:IL,1,K,1)
            YW(I:IL,K) = X(I:IL,1,K,2)
            ZW(I:IL,K) = X(I:IL,1,K,3)
         END DO

      END IF

C     Reaccount for trailing edge slopes except for gradients (for now):

      IF (NCALL <= 0) THEN

         CALL FIXWAKE (IL, KL, ITL, ITU, ILAST, 2, KTIP - 1, TFIWORK,
     >                 XW, YW, ZW)
      END IF

C     Transfer the wake sheet:

      DO K = 2, KL
         DO I = ITU + 1, IL - 1
            J = IL + 1 - I
            X(I,1,K,1) = XW(I,K)
            X(J,1,K,1) = XW(I,K)
            X(I,1,K,2) = YW(I,K)
            X(J,1,K,2) = YW(I,K)
            X(I,1,K,3) = ZW(I,K)
            X(J,1,K,3) = ZW(I,K)
         END DO
      END DO


C     The KL and KTIP faces depend on the topology.

      D1KTIP = DRADIAL(KTIP,1)

      IF (NNTIP == 0) THEN ! C-H everywhere

         DO I = 1, IL      ! Missing J = JL edge
            X(I,JL,KL,1) = XOUT(I) + ZMAX * TANSWEEP
            X(I,JL,KL,2) = YOUT(I)
         END DO

         ZWK = ZW(1,KL)
         X(1:IL,1:JL,KL,3) = ZWK

C        Lines off the (extended) trailing edge allow the radial spacing
C        to be independent of the spacing at I = 1, ILE, & IL:

         DO I = ITL, ITU, ITU - ITL

C           We can never warp it if we're capturing trailing edge normals,
C           so indirect mode must go through the Euler spacing always.

            CALL OFFTEDGE (I, KL, DRADIALEU(KL,2), D2JLEU,
     >                     NBLAYEREU, RBLAYEREU, TRADIAL, IWRIT)

            IF (WARPNS) THEN

               CALL REGRID_LINE (IL, JL, KL, I, KL, X, DRADIAL(KL,2),
     >                           D2JL, NBLAYER, RBLAYER, TRADIAL,
     >                           XRADIAL, YRADIAL, ZRADIAL, IWRIT)
            END IF

         END DO

C        Radial line forward of the (extended) leading edge at K = KL:
C        We can never warp it if we're capturing leading edge angles,
C        so indirect mode must go through the Euler spacing always.

         CALL OFFLEDGE (IL, JL, KL, ILE, KL, X, DRADIALEU(KL,1),
     >                  D2JLEU, NBLAYEREU, RBLAYEREU, TRADIAL,
     >                  XRADIAL, YRADIAL, ZRADIAL, IWRIT)

         IF (WARPNS) THEN

            CALL REGRID_LINE (IL, JL, KL, ILE, KL, X, DRADIAL(KL,1),
     >                        D2JL, NBLAYER, RBLAYER, TRADIAL,
     >                        XRADIAL, YRADIAL, ZRADIAL, IWRIT)
         END IF


C        4-part TFI and smoothing:

         IF (FULLGRID) THEN

            DO N = 1, 4

               CALL TFI2D (IL, IFACE(N), IFACE(N+1), 1, JL,
     >                     X(1,1,KL,1), X(1,1,KL,2), TFIWORK)
            END DO

            IF (PRINT2D)
     >         WRITE (IWRIT,'(/,2A,I4)') E2D, 'on plane KL =', KL

            L = 2 ! K-plane regions 1 and 4

            CALL ELLIP2D (IL, JL, 1, ITL, 1, JL,
     >                    X(1,1,KL,1), X(1,1,KL,2), DFACEK,    PRINT2D,
     >                    BGMODE2D(L), FGMODE2D(L), SPMODE(L),
     >                    ITMAX2D(L),  ITFLOAT,     ITFREEZE,  JLAYER,
     >                    CONV2D(L),   CONVMIN(L),  DMAX2D,    OMG2D(L),
     >                    POWERI(L),   POWERJ(L),   EXPI12D(L),
     >                    EXPI22D(L),  EXPJ12D(L),  EXPJ22D(L),
     >                    FGLIM2D,     URFG2D,      URFLOAT)

            L = 8 ! KL plane regions 2 and 3

            CALL ELLIP2D (IL, JL, ITL, ILE, 1, JL,
     >                    X(1,1,KL,1), X(1,1,KL,2), DFACEK,    PRINT2D,
     >                    BGMODE2D(L), FGMODE2D(L), SPMODE(L),
     >                    ITMAX2D(L),  ITFLOAT,     ITFREEZE,  JLAYER,
     >                    CONV2D(L),   CONVMIN(L),  DMAX2D,    OMG2D(L),
     >                    POWERI(L),   POWERJ(L),   EXPI12D(L),
     >                    EXPI22D(L),  EXPJ12D(L),  EXPJ22D(L),
     >                    FGLIM2D,     URFG2D,      URFLOAT)

            CALL EXCHANGE (BGMODE2D(L), BGMODE) ! Exchange I1 & I2 controls
            CALL EXCHANGE (FGMODE2D(L), FGMODE)

            CALL ELLIP2D (IL, JL, ILE, ITU, 1, JL,
     >                    X(1,1,KL,1), X(1,1,KL,2), DFACEK,    PRINT2D,
     >                    BGMODE,      FGMODE,      SPMODE(L),
     >                    ITMAX2D(L),  ITFLOAT,     ITFREEZE,  JLAYER,
     >                    CONV2D(L),   CONVMIN(L),  DMAX2D,    OMG2D(L),
     >                   -POWERI(L),   POWERJ(L),   EXPI22D(L),
     >                    EXPI12D(L),  EXPJ12D(L),  EXPJ22D(L),
     >                    FGLIM2D,     URFG2D,      URFLOAT)

            L = 2
            CALL EXCHANGE (BGMODE2D(L), BGMODE)
            CALL EXCHANGE (FGMODE2D(L), FGMODE)

            CALL ELLIP2D (IL, JL, ITU, IL, 1, JL,
     >                    X(1,1,KL,1), X(1,1,KL,2), DFACEK,    PRINT2D,
     >                    BGMODE,      FGMODE,      SPMODE(L),
     >                    ITMAX2D(L),  ITFLOAT,     ITFREEZE,  JLAYER,
     >                    CONV2D(L),   CONVMIN(L),  DMAX2D,    OMG2D(L),
     >                   -POWERI(L),   POWERJ(L),   EXPI22D(L),
     >                    EXPI12D(L),  EXPJ12D(L),  EXPJ22D(L),
     >                    FGLIM2D,     URFG2D,      URFLOAT)

         ELSE ! Perturb the KL plane interior

            DO N = 1, 4 ! Quadrants

               I1 = IFACE(N)
               I2 = IFACE(N+1)

               CALL WARPQ3D (1, IL, 1, JL, 1, KL, I1, I2, 1, JL, KL, KL,
     >                       X0(1,1,1,1), X0(1,1,1,2), X0(1,1,1,3), S0,
     >                       DFACEI,     DFACEJ,     DFACEK,
     >                       X(1,1,1,1), X(1,1,1,2), X(1,1,1,3))
            END DO

         END IF

C        Remaining C outer boundaries:

         DO K = 2, KL - 1
            XSHIFT = ZUBACK(JL,K) * TANSWEEP
            DO I = 2, IL - 1
               X(I,JL,K,1) = XOUT(I) + XSHIFT
               X(I,JL,K,2) = YOUT(I)
               X(I,JL,K,3) = ZUBACK(JL,K) ! See OUTER
            END DO
         END DO


C        A C grid at the tip plane is an artificial boundary.
C        Radial edges off the trailing edge:

         K = KTIP

         DO I = ITL, ITU, ITU - ITL

            CALL OFFTEDGE (I, K, DRADIALEU(K,2), D2JLEU,
     >                     NBLAYEREU, RBLAYEREU, TRADIAL, IWRIT)

            IF (WARPNS) THEN

               CALL REGRID_LINE (IL, JL, KL, I, K, X, DRADIAL(K,2),
     >                           D2JL, NBLAYER, RBLAYER, TRADIAL,
     >                           XRADIAL, YRADIAL, ZRADIAL, IWRIT)
            END IF

         END DO

C        Radial line forward of the tip leading edge:

         CALL OFFLEDGE (IL, JL, KL, ILE, K, X, DRADIALEU(K,1),
     >                  D2JLEU, NBLAYEREU, RBLAYEREU, TRADIAL,
     >                  XRADIAL, YRADIAL, ZRADIAL, IWRIT)

         IF (WARPNS) THEN

            CALL REGRID_LINE (IL, JL, KL, ILE, K, X, DRADIAL(K,1),
     >                        D2JL, NBLAYER, RBLAYER, TRADIAL,
     >                        XRADIAL, YRADIAL, ZRADIAL, IWRIT)
         END IF

         IF (FULLGRID) THEN

C           Quasi-3D TFI of tip plane and 2D smoothing of X and Y, in 4 parts:

            DO N = 1, 4

               CALL TFIQ3D (IL, IFACE(N), IFACE(N+1), 1, JL, X(1,1,K,1),
     >                      X(1,1,K,2), X(1,1,K,3), TFIWORK)
            END DO

            IF (PRINT2D)
     >         WRITE (IWRIT,'(/,2A,I4)') E2D, 'on plane KTIP =', K

            L = 2 ! Regions 1 and 4

            CALL ELLIP2D (IL, JL, 1, ITL, 1, JL,
     >                    X(1,1,K,1),  X(1,1,K,2),  DFACEK,    PRINT2D,
     >                    BGMODE2D(L), FGMODE2D(L), SPMODE(L),
     >                    ITMAX2D(L),  ITFLOAT,     ITFREEZE,  JLAYER,
     >                    CONV2D(L),   CONVMIN(L),  DMAX2D,    OMG2D(L),
     >                    POWERI(L),   POWERJ(L),   EXPI12D(L),
     >                    EXPI22D(L),  EXPJ12D(L),  EXPJ22D(L),
     >                    FGLIM2D,     URFG2D,      URFLOAT)

            L = 9 ! Regions 2 and 3

            CALL ELLIP2D (IL, JL, ITL, ILE, 1, JL,
     >                    X(1,1,K,1),  X(1,1,K,2),  DFACEK,    PRINT2D,
     >                    BGMODE2D(L), FGMODE2D(L), SPMODE(L),
     >                    ITMAX2D(L),  ITFLOAT,     ITFREEZE,  JLAYER,
     >                    CONV2D(L),   CONVMIN(L),  DMAX2D,    OMG2D(L),
     >                    POWERI(L),   POWERJ(L),   EXPI12D(L),
     >                    EXPI22D(L),  EXPJ12D(L),  EXPJ22D(L),
     >                    FGLIM2D,     URFG2D,      URFLOAT)

            CALL EXCHANGE (BGMODE2D(L), BGMODE)
            CALL EXCHANGE (FGMODE2D(L), FGMODE)

            CALL ELLIP2D (IL, JL, ILE, ITU, 1, JL,
     >                    X(1,1,K,1),  X(1,1,K,2),  DFACEK,    PRINT2D,
     >                    BGMODE,      FGMODE,      SPMODE(L),
     >                    ITMAX2D(L),  ITFLOAT,     ITFREEZE,  JLAYER,
     >                    CONV2D(L),   CONVMIN(L),  DMAX2D,    OMG2D(L),
     >                   -POWERI(L),   POWERJ(L),   EXPI22D(L),
     >                    EXPI12D(L),  EXPJ12D(L),  EXPJ22D(L),
     >                    FGLIM2D,     URFG2D,      URFLOAT)

            L = 2
            CALL EXCHANGE (BGMODE2D(L), BGMODE)
            CALL EXCHANGE (FGMODE2D(L), FGMODE)

            CALL ELLIP2D (IL, JL, ITU, IL, 1, JL,
     >                    X(1,1,K,1),  X(1,1,K,2),  DFACEK,    PRINT2D,
     >                    BGMODE,      FGMODE,      SPMODE(L),
     >                    ITMAX2D(L),  ITFLOAT,     ITFREEZE,  JLAYER,
     >                    CONV2D(L),   CONVMIN(L),  DMAX2D,    OMG2D(L),
     >                   -POWERI(L),   POWERJ(L),   EXPI22D(L),
     >                    EXPI12D(L),  EXPJ12D(L),  EXPJ22D(L),
     >                    FGLIM2D,     URFG2D,      URFLOAT)

         ELSE ! Perturb the tip K plane

            DO N = 1, 4

               I1 = IFACE(N)
               I2 = IFACE(N+1)

               CALL WARPQ3D (1, IL, 1, JL, 1, KL, I1, I2, 1, JL, K, K,
     >                       X0(1,1,1,1), X0(1,1,1,2), X0(1,1,1,3), S0,
     >                       DFACEI,     DFACEJ,     DFACEK,
     >                       X(1,1,1,1), X(1,1,1,2), X(1,1,1,3))
            END DO

         END IF


C        Set up edges of an artificial plane boundary forward of the wing:

         DO J = 1, JL
            PLANEI(J,1,1)    = XWALL(ILE,J)    ! Water-line edge
            PLANEI(J,1,2)    = YWALL(ILE,J)
            PLANEI(J,1,3)    = ZWALL(ILE,J)
            PLANEI(J,KTIP,1) = X(ILE,J,KTIP,1) ! Edge forward of K = KTIP
            PLANEI(J,KTIP,2) = X(ILE,J,KTIP,2)
            PLANEI(J,KTIP,3) = X(ILE,J,KTIP,3)
            PLANEI(J,KL,1)   = X(ILE,J,KL,1)   ! Edge forward of K = KL
            PLANEI(J,KL,2)   = X(ILE,J,KL,2)
            PLANEI(J,KL,3)   = X(ILE,J,KL,3)
         END DO

         DO K = 1, KL - 1
            PLANEI(1, K,1) = XW(ILE,K)         ! Wing leading edge
            PLANEI(1, K,2) = YW(ILE,K)
            PLANEI(1, K,3) = ZW(ILE,K)
            PLANEI(JL,K,1) = X(ILE,JL,K,1)     ! Upstream edge
            PLANEI(JL,K,2) = X(ILE,JL,K,2)
            PLANEI(JL,K,3) = X(ILE,JL,K,3)
         END DO

C        Fill the artificial I = ILE sheet as two parts:

         IF (WINGONLY) THEN
            JLAST = 2 * JL / 3
         ELSE
            JLAST = JCR
         END IF

         IF (FULLGRID) THEN

            DO M = 1, 2

               K1 = KFACE(M)
               K2 = KFACE(M+1)

               CALL TFIQ3D (JL, 1, JL, K1, K2, PLANEI(1,1,1),
     >                      PLANEI(1,1,2), PLANEI(1,1,3), TFIWORK)
            END DO

C           Account for any leading edge droop:

            CALL FIXLEDGE (IL, JL, KL, ILE, JLAST, 2, KTIP-1, PLANEI,
     >                     TRADIAL, X)

C           Smooth it?

            L = 10
            ITERS = ITMAX2D(L) / 2

            IF (PRINT2D) WRITE (IWRIT,'(/,2A,I4)')
     >         E2D, '& -Q3D on plane ILE =', ILE

            CALL ELLIP2D  (JL, KL, 1, JL, 1, KTIP,     PLANEI(1,1,3),
     >                     PLANEI(1,1,1), DFACEI,      PRINT2D,
     >                     BGMODE2D(L),  FGMODE2D(L),  SPMODE(L),
     >                     ITERS,        ITFLOAT,      ITFREEZE, JLAYER,
     >                     CONV2D(L),    CONVMIN(L),   DMAX2D, OMG2D(L),
     >                     POWERI(L),    POWERJ(L),    EXPI12D(L),
     >                     EXPI22D(L),   EXPJ12D(L),   EXPJ22D(L),
     >                     FGLIM2D,      URFG2D,       URFLOAT)

            CALL ELLIPQ3D (JL, KL, 1, JL, 1, KTIP,     PLANEI(1,1,3),
     >                     PLANEI(1,1,1),PLANEI(1,1,2),DFACEI, PRINT2D,
     >                     BGMODE2D(L),  FGMODE2D(L),  SPMODE(L),
     >                     ITERS,        ITFLOAT,      ITFREEZE, JLAYER,
     >                     CONV2D(L),    CONVMIN(L),   DMAX2D, OMG2D(L),
     >                     POWERI(L),    POWERJ(L),    EXPI12D(L),
     >                     EXPI22D(L),   EXPJ12D(L),   EXPJ22D(L),
     >                     FGLIM2D,      URFG2D,       URFLOAT)

C           Re-account for any leading edge droop:

            CALL FIXLEDGE (IL, JL, KL, ILE, JLAST, 2, KTIP-1, PLANEI,
     >                     TRADIAL, X)

C           And beyond the tip:

            CALL ELLIPQ3D (JL, KL, 1, JL, KTIP, KL,    PLANEI(1,1,3),
     >                     PLANEI(1,1,1),PLANEI(1,1,2),DFACEI, PRINT2D,
     >                     BGMODE2D(L),  FGMODE2D(L),  SPMODE(L),
     >                     ITMAX2D(L)/2, ITFLOAT,      ITFREEZE, JLAYER,
     >                     CONV2D(L),    CONVMIN(L),   DMAX2D, OMG2D(L),
     >                     POWERI(L),    POWERJ(L),    EXPI12D(L),
     >                     EXPI22D(L),   EXPJ12D(L),   EXPJ22D(L),
     >                     FGLIM2D,      URFG2D,       URFLOAT)

            J = JL - 1
            DO K = 2, KL - 1 ! Transfer the artificial boundary
               X(ILE,2:J,K,1) = PLANEI(2:J,K,1)
               X(ILE,2:J,K,2) = PLANEI(2:J,K,2)
               X(ILE,2:J,K,3) = PLANEI(2:J,K,3)
            END DO

         ELSE ! Perturb the sheet forward of the leading edge in two parts

            CALL WARPQ3D (1, IL, 1, JL, 1, KL, ILE, ILE, 1, JL, 1, KTIP,
     >                    X0(1,1,1,1), X0(1,1,1,2), X0(1,1,1,3), S0,
     >                    DFACEI,     DFACEJ,     DFACEK,
     >                    X(1,1,1,1), X(1,1,1,2), X(1,1,1,3))

C           Account for any leading edge deflection:

            J = JL - 1
            DO K = 2, KTIP - 1 ! Transfer the ILE plane in front of the wing
               PLANEI(2:J,K,1) = X(ILE,2:J,K,1)
               PLANEI(2:J,K,2) = X(ILE,2:J,K,2)
               PLANEI(2:J,K,3) = X(ILE,2:J,K,3)
            END DO

            IF (NCALL <= 0) THEN

               CALL FIXLEDGE (IL, JL, KL, ILE, JLAST, 2, KTIP-1,
     >                        PLANEI, TRADIAL, X)
            END IF

            DO K = 2, KTIP - 1
               X(ILE,2:J,K,1) = PLANEI(2:J,K,1)
               X(ILE,2:J,K,2) = PLANEI(2:J,K,2)
               X(ILE,2:J,K,3) = PLANEI(2:J,K,3)
            END DO

            CALL WARPQ3D (1, IL, 1, JL, 1, KL, ILE, ILE, 1,JL, KTIP, KL,
     >                    X0(1,1,1,1), X0(1,1,1,2), X0(1,1,1,3), S0,
     >                    DFACEI,     DFACEJ,     DFACEK,
     >                    X(1,1,1,1), X(1,1,1,2), X(1,1,1,3))
         END IF


C        Break the blocks above and below the wing into two with
C        spanwise/vertical planes at the trailing edge:

C        Transfer I = ITL plane edges as needed for filling/smoothing:

         IF (FULLGRID) THEN

            DO J = 1, JL, JL - 1
               PLANEI(J,1:KL,1) = X(ITL,J,1:KL,1)
               PLANEI(J,1:KL,2) = X(ITL,J,1:KL,2)
               PLANEI(J,1:KL,3) = X(ITL,J,1:KL,3)
            END DO

            DO K = 1, KTIP, KTIP - 1
               PLANEI(1:JL,K,1) = X(ITL,1:JL,K,1)
               PLANEI(1:JL,K,2) = X(ITL,1:JL,K,2)
               PLANEI(1:JL,K,3) = X(ITL,1:JL,K,3)
            END DO

            CALL TFIQ3D (JL, 1, JL, 1, KTIP, PLANEI(1,1,2),
     >                   PLANEI(1,1,3), PLANEI(1,1,1), TFIWORK)

            L = 15

            IF (PRINT2D)
     >         WRITE (IWRIT,'(/,2A,I4)') E2D, 'on plane ITL =', ITL

CCC         CALL ELLIPQ3D (JL, KL, 1, JL, 1, KTIP,     PLANEI(1,1,2),
CCC     >                  PLANEI(1,1,3),PLANEI(1,1,1),DFACEI, PRINT2D,
CCC     >                  BGMODE2D(L),  FGMODE2D(L),  SPMODE(L),
CCC     >                  ITMAX2D(L),   ITFLOAT,      ITFREEZE, JLAYER,
CCC     >                  CONV2D(L),    CONVMIN(L),   DMAX2D, OMG2D(L),
CCC     >                  POWERI(L),    POWERJ(L),    EXPI12D(L),
CCC     >                  EXPI22D(L),   EXPJ12D(L),   EXPJ22D(L),
CCC     >                  FGLIM2D,      URFG2D,       URFLOAT)

C           Until ELLIPQ3D works, smooth just Y & Zs:

            CALL ELLIP2D  (JL, KL, 1, JL, 1, KTIP,     PLANEI(1,1,3),
     >                     PLANEI(1,1,2), DFACEI,      PRINT2D,
     >                     BGMODE2D(L),  FGMODE2D(L),  SPMODE(L),
     >                     ITMAX2D(L),   ITFLOAT,      ITFREEZE, JLAYER,
     >                     CONV2D(L),    CONVMIN(L),   DMAX2D, OMG2D(L),
     >                     POWERI(L),    POWERJ(L),    EXPI12D(L),
     >                     EXPI22D(L),   EXPJ12D(L),   EXPJ22D(L),
     >                     FGLIM2D,      URFG2D,       URFLOAT)

            J = JL - 1
            DO K = 2, KTIP - 1
               X(ITL,2:J,K,1) = PLANEI(2:J,K,1)
               X(ITL,2:J,K,2) = PLANEI(2:J,K,2)
               X(ITL,2:J,K,3) = PLANEI(2:J,K,3)
            END DO

C           ... and beyond the tip, lower half:

            PLANEI(1:JL,KL,1) = X(ITL,1:JL,KL,1)
            PLANEI(1:JL,KL,2) = X(ITL,1:JL,KL,2)
            PLANEI(1:JL,KL,3) = X(ITL,1:JL,KL,3)

            CALL TFIQ3D (JL, 1, JL, KTIP, KL, PLANEI(1,1,2),
     >                   PLANEI(1,1,3), PLANEI(1,1,1), TFIWORK)

            CALL ELLIP2D  (JL, KL, 1, JL, KTIP, KL,    PLANEI(1,1,3),
     >                     PLANEI(1,1,2), DFACEI,      PRINT2D,
     >                     BGMODE2D(L),  FGMODE2D(L),  SPMODE(L),
     >                     ITMAX2D(L),   ITFLOAT,      ITFREEZE, JLAYER,
     >                     CONV2D(L),    CONVMIN(L),   DMAX2D, OMG2D(L),
     >                     POWERI(L),    POWERJ(L),    EXPI12D(L),
     >                     EXPI22D(L),   EXPJ12D(L),   EXPJ22D(L),
     >                     FGLIM2D,      URFG2D,       URFLOAT)

            DO K = KTIP + 1, KL - 1
               X(ITL,2:J,K,1) = PLANEI(2:J,K,1)
               X(ITL,2:J,K,2) = PLANEI(2:J,K,2)
               X(ITL,2:J,K,3) = PLANEI(2:J,K,3)
            END DO


C           Repeat for the plane above the trailing edge:

            DO J = 1, JL, JL - 1
               PLANEI(J,1:KL,1) = X(ITU,J,1:KL,1)
               PLANEI(J,1:KL,2) = X(ITU,J,1:KL,2)
               PLANEI(J,1:KL,3) = X(ITU,J,1:KL,3)
            END DO

            DO K = 1, KTIP, KTIP - 1
               PLANEI(1:JL,K,1) = X(ITU,1:JL,K,1)
               PLANEI(1:JL,K,2) = X(ITU,1:JL,K,2)
               PLANEI(1:JL,K,3) = X(ITU,1:JL,K,3)
            END DO

            CALL TFIQ3D (JL, 1, JL, 1, KTIP, PLANEI(1,1,2),
     >                   PLANEI(1,1,3), PLANEI(1,1,1), TFIWORK)

            CALL EXCHANGE (BGMODE2D(L), BGMODE)
            CALL EXCHANGE (FGMODE2D(L), FGMODE)

            IF (PRINT2D)
     >         WRITE (IWRIT,'(/,2A,I4)') E2D, 'on plane ITU =', ITU

CCC         CALL ELLIPQ3D  (JL, KL, 1, JL, 1, KTIP,    PLANEI(1,1,2),
CCC     >                  PLANEI(1,1,3), PLANEI,      DFACEI, PRINT2D,
CCC     >                  BGMODE,       FGMODE,       SPMODE(L),
CCC     >                  ITMAX2D(L),   ITFLOAT,      ITFREEZE, JLAYER,
CCC     >                  CONV2D(L),    CONVMIN(L),   DMAX2D, OMG2D(L),
CCC     >                  POWERI(L),    POWERJ(L),    EXPI22D(L),
CCC     >                  EXPI12D(L),   EXPJ12D(L),   EXPJ22D(L),
CCC     >                  FGLIM2D,      URFG2D,       URFLOAT)

            CALL ELLIP2D  (JL, KL, 1, JL, 1, KTIP,     PLANEI(1,1,2),
     >                     PLANEI(1,1,3), DFACEI,      PRINT2D,
     >                     BGMODE,       FGMODE,       SPMODE(L),
     >                     ITMAX2D(L),   ITFLOAT,      ITFREEZE, JLAYER,
     >                     CONV2D(L),    CONVMIN(L),   DMAX2D, OMG2D(L),
     >                     POWERI(L),    POWERJ(L),    EXPI22D(L),
     >                     EXPI12D(L),   EXPJ12D(L),   EXPJ22D(L),
     >                     FGLIM2D,      URFG2D,       URFLOAT)

            J = JL - 1
            DO K = 2, KTIP - 1
               X(ITU,2:J,K,1) = PLANEI(2:J,K,1)
               X(ITU,2:J,K,2) = PLANEI(2:J,K,2)
               X(ITU,2:J,K,3) = PLANEI(2:J,K,3)
            END DO

            PLANEI(1:JL,KL,1) = X(ITU,1:JL,KL,1)
            PLANEI(1:JL,KL,2) = X(ITU,1:JL,KL,2)
            PLANEI(1:JL,KL,3) = X(ITU,1:JL,KL,3)

            CALL TFIQ3D (JL, 1, JL, KTIP, KL, PLANEI(1,1,2),
     >                   PLANEI(1,1,3), PLANEI(1,1,1), TFIWORK)

            CALL ELLIP2D  (JL, KL, 1, JL, KTIP, KL,    PLANEI(1,1,2),
     >                     PLANEI(1,1,3), DFACEI,      PRINT2D,
     >                     BGMODE,       FGMODE,       SPMODE(L),
     >                     ITMAX2D(L),   ITFLOAT,      ITFREEZE, JLAYER,
     >                     CONV2D(L),    CONVMIN(L),   DMAX2D, OMG2D(L),
     >                     POWERI(L),    POWERJ(L),    EXPI22D(L),
     >                     EXPI12D(L),   EXPJ12D(L),   EXPJ22D(L),
     >                     FGLIM2D,      URFG2D,       URFLOAT)

            DO K = KTIP + 1, KL - 1
               X(ITU,2:J,K,1) = PLANEI(2:J,K,1)
               X(ITU,2:J,K,2) = PLANEI(2:J,K,2)
               X(ITU,2:J,K,3) = PLANEI(2:J,K,3)
            END DO

         ELSE ! Perturb the I = ITL and ITU planes in 4 parts

            DO N = 2, 4, 2 ! Quadrants

               I = IFACE(N)

               DO M = 1, 2

                  K1 = KFACE(M)
                  K2 = KFACE(M+1)

                  CALL WARPQ3D (1, IL, 1, JL, 1, KL, I,I, 1, JL, K1, K2,
     >                          X0(1,1,1,1), X0(1,1,1,2), X0(1,1,1,3),
     >                          S0, DFACEI, DFACEJ, DFACEK,
     >                          X(1,1,1,1), X(1,1,1,2), X(1,1,1,3))
               END DO

            END DO

         END IF


      ELSE ! C-H + C-O case is more awkward

C ***    Boundaries at t.e. not implemented yet ***
C ***    Capturing slat & flap angles off the wing not handled yet ***
C ***    Warping not handled yet ***

C        Easier inboard part of C-H portion:

         DO K = 2, KTOUTER ! KTOUTER is last full-size outer boundary
            DO I = 2, IL - 1
               X(I,JL,K,2) = YOUT(I)
               R = ABS (YOUT(I) / YOUT(IL))
               X(I,JL,K,3) = R * ZUBACK(JL,K) + (ONE - R) * ZW(1,K)
               X(I,JL,K,1) = X(I,JL,K,3) * TANSWEEP + XOUT(I)
            END DO
         END DO

C        Rest of outer C boundary by rotation of K = KTOUTER C about wing tip,
C        plus shearing to follow the wing sweep.
C        The radius of the tip quadrant is YOUT(IL).

         ZWK = ZW(ILE,KTIP) ! Circle centers are at (ZWK, 0.)
         DO K = KTOUTER + 1, KL
            THETA = (ARCOUTER(KL) - ARCOUTER(K)) / YOUT(IL)
            COSTH = COS (THETA)
            SINTH = SIN (THETA)
            DO I = ILE, IL - 1
               J = IL + 1 - I
               X(I,JL,K,2) = YOUT(I) * SINTH
               X(J,JL,K,2) =-X(I,JL,K,2)
               R = ABS (X(I,JL,K,2) / YOUT(IL))
               X(I,JL,K,3) = YOUT(I) * COSTH + R*ZWK + (ONE - R)*ZW(1,K)
               X(J,JL,K,3) = X(I,JL,K,3)
               X(I,JL,K,1) = X(I,JL,K,3) * TANSWEEP + XOUT(I)
               X(J,JL,K,1) = X(I,JL,K,1)
            END DO
         END DO

C        Introduce an artificial boundary plane forward of the tip leading edge.
C        The outboard edge is needed as part of the K = KL boundary, and the
C        interior lines are needed for K planes to be filled by halves:

C ****** Fix this as above if it's ever used

         CALL OFFLEDGE (IL, JL, KL, ILE, KTIP, X, D1KTIP,
     >                  D2JL, NBLAYER, RBLAYER, TRADIAL,
     >                  XRADIAL, YRADIAL, ZRADIAL, IWRIT)

         PLANEI(1:JL,KTIP,1) = XRADIAL
         PLANEI(1:JL,KTIP,2) = YRADIAL
         PLANEI(1:JL,KTIP,3) = ZWK

C        Singular line forward of the tip for K >= KTIP:

         DO K = KTIP + 1, KL
            X(ILE,1:JL,K,1) = XRADIAL
            X(ILE,1:JL,K,2) = YRADIAL
            X(ILE,1:JL,K,3) = ZWK
         END DO

C        Set up remaining edges:

         DO K = 1, KTIP - 1
            PLANEI(1, K,1) = XW(ILE,K)   ! Wing leading edge
            PLANEI(1, K,2) = YW(ILE,K)
            PLANEI(1, K,3) = ZW(ILE,K)
            PLANEI(JL,K,2) = ZERO        ! Upstream edge
            PLANEI(JL,K,3) = ZUBACK(1,K) ! May want to use relative l.e. distrb.
            PLANEI(JL,K,1) = ZUBACK(1,K) * TANSWEEP + XOUT(ILE)
         END DO

         DO J = 2, JL - 1
            PLANEI(J,1,1) = XWALL(ILE,J) ! Water-line edge
            PLANEI(J,1,2) = YWALL(ILE,J)
            PLANEI(J,1,3) = ZWALL(ILE,J)
         END DO

         IF (FULLGRID) THEN

C           Fill the I = ILE plane:

            CALL TFIQ3D (JL, 1, JL, 1, KTIP, PLANEI(1,1,1),
     >                   PLANEI(1,1,2), PLANEI(1,1,3), TFIWORK)

C           Smooth it?  (Remember Y is not constant.)

            L = 10

            IF (PRINT2D) WRITE (IWRIT,'(/,2A)')
     >         E2D, 'on sheet forward of leading edge:'

            CALL ELLIP2D (JL, KL, 1, JL, 1, KTIP,    PLANEI(1,1,1),
     >                    PLANEI(1,1,3), DFACEI,     PRINT2D,
     >                    BGMODE2D(L),  FGMODE2D(L), SPMODE(L),
     >                    ITMAX2D(L),   ITFLOAT,     ITFREEZE, JLAYER,
     >                    CONV2D(L),    CONVMIN(L),  DMAX2D,   OMG2D(L),
     >                    POWERI(L),    POWERJ(L),   EXPI12D(L),
     >                    EXPI22D(L),   EXPJ12D(L),  EXPJ22D(L),
     >                    FGLIM2D,      URFG2D,      URFLOAT)

            DO K = 1, KTIP - 1
               X(ILE,2:JL,K,1) = PLANEI(2:JL,K,1)
               X(ILE,2:JL,K,2) = PLANEI(2:JL,K,2)
               X(ILE,2:JL,K,3) = PLANEI(2:JL,K,3)
            END DO

         ELSE ! Perturb it if we ever use this option

         END IF


C        Fill the KL half planes via quasi-3D TFI:

         DO N = 1, 3, 2 ! Quadrants

            CALL TFIQ3D (IL, IFACE(N), IFACE(N+2), 1, JL,
     >                   X(1,1,KL,1), X(1,1,KL,2), X(1,1,KL,3), TFIWORK)
         END DO

      END IF

      IF (NCALL < 0) THEN
         CALL SECOND (CPU2)
         CALL CPUTIME (CPU1, CPU2, 'GRID_FACES',
     >                 'to grid the subblock faces', IWRIT)
      END IF

  999 RETURN

      END SUBROUTINE GRID_FACES

C***********************************************************************
C
      SUBROUTINE GRID_TFI3D (NCALL)
C
C     Fill all subblock grids via 3-D TFI given the subblock faces.
C
C***********************************************************************

C     Global variables:

      USE ELLIPTIC
      USE FACES
      USE LIM
      USE LOGIC
      USE LUNS
      USE OPTIONS
      USE SPACING
      USE TFI
      USE XYZ
      USE XYZ0

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) :: NCALL

C     Local constants:

      INTEGER, PARAMETER ::
     >   NITSONI = 10  ! Max. # iters. for Soni-type blending fns. in TFI3D

      REAL, PARAMETER ::
     >   TOLTFI3D = 1. ! Tolerance for terminating this iteration, relative to
                       ! an estimate of the smallest normalized grid spacing
C     Local variables:

      INTEGER
     >   I1, I2, J, K, L, N

      REAL, ALLOCATABLE ::
     >   BF(:,:,:,:)

      REAL
     >   PLANEI(JL,KL,3), TIPCHORD, TOLSONI, YRANGE

      REAL(KIND=4) CPU1, CPU2
C     Execution:

      IF (NCALL == -3) CALL SECOND (CPU1)

C     TFI3D's iteration tolerance is on normalized arc lengths.
C     All subblocks are affected by the tip trailing edge radial increment.
C     Use it as a fraction of the (vertical) data range:

      YRANGE   = X(1,1,1,2)      - X(1,JL,1,2)
      TIPCHORD = X(ITL,1,KTIP,1) - X(ILE,1,KTIP,1)
      TOLSONI  = TOLTFI3D * (DRADIAL(KTIP,2) * TIPCHORD) / YRANGE
      TOLSONI  = MIN (TOLSONI, 0.0001)

      ALLOCATE (BF(IL,JL,KL,3))

      IF (NNTIP == 0) THEN ! C-H everywhere

C        Initialize lower and upper halves via 3-D TFI.

C        Poor TFI3D results above the wake for the wing/body case are
C        improved by introducing an I plane at the end of body:

         IF (.NOT. WINGONLY) THEN

            DO J = 1, JL, JL - 1
               PLANEI(J,1:KL,1) = X(IBL,J,1:KL,1)
               PLANEI(J,1:KL,2) = X(IBL,J,1:KL,2)
               PLANEI(J,1:KL,3) = X(IBL,J,1:KL,3)
            END DO

            DO K = 1, KTIP, KTIP - 1
               PLANEI(1:JL,K,1) = X(IBL,1:JL,K,1)
               PLANEI(1:JL,K,2) = X(IBL,1:JL,K,2)
               PLANEI(1:JL,K,3) = X(IBL,1:JL,K,3)
            END DO

            CALL TFIQ3D (JL, 1, JL, 1, KTIP, PLANEI(1,1,2),
     >                   PLANEI(1,1,3), PLANEI(1,1,1), TFIWORK)

            L = 15

            IF (PRINT2D) WRITE (IWRIT,'(/,A,I4)')
     >         ' ELLIP2D on plane IBL =', IBL

            CALL ELLIP2D (JL, KL, 1, JL, 1, KTIP,   PLANEI(1,1,3),
     >                    PLANEI(1,1,2), DFACEI,    PRINT2D,
     >                    BGMODE2D(L),  FGMODE2D(L),SPMODE(L),
     >                    ITMAX2D(L),   ITFLOAT,    ITFREEZE, JLAYER,
     >                    CONV2D(L),    CONVMIN(L), DMAX2D, OMG2D(L),
     >                    POWERI(L),    POWERJ(L),  EXPI12D(L),
     >                    EXPI22D(L),   EXPJ12D(L), EXPJ22D(L),
     >                    FGLIM2D,      URFG2D,     URFLOAT)

            J = JL - 1
            DO K = 2, KTIP - 1
               X(IBL,2:J,K,1) = PLANEI(2:J,K,1)
               X(IBL,2:J,K,2) = PLANEI(2:J,K,2)
               X(IBL,2:J,K,3) = PLANEI(2:J,K,3)
            END DO

            I1 = 1
            I2 = IBL

            DO N = 1, 2

               CALL TFI3D (IL, JL, KL, I1, I2, 1, JL, 1, KTIP,
     >                     X(1,1,1,1), X(1,1,1,2), X(1,1,1,3),
     >                     S0, BF, PRINT3D, NITSONI, TOLSONI)
               I1 = IBL
               I2 = ITL
            END DO

C           Likewise for above the wake:

            DO J = 1, JL, JL - 1
               PLANEI(J,1:KL,1) = X(IBU,J,1:KL,1)
               PLANEI(J,1:KL,2) = X(IBU,J,1:KL,2)
               PLANEI(J,1:KL,3) = X(IBU,J,1:KL,3)
            END DO

            DO K = 1, KTIP, KTIP - 1
               PLANEI(1:JL,K,1) = X(IBU,1:JL,K,1)
               PLANEI(1:JL,K,2) = X(IBU,1:JL,K,2)
               PLANEI(1:JL,K,3) = X(IBU,1:JL,K,3)
            END DO

            CALL TFIQ3D (JL, 1, JL, 1, KTIP, PLANEI(1,1,2),
     >                   PLANEI(1,1,3), PLANEI(1,1,1), TFIWORK)

            IF (PRINT2D) WRITE (IWRIT,'(/,A,I4)')
     >         ' ELLIP2D on plane IBU =', IBU

            CALL ELLIP2D (JL, KL, 1, JL, 1, KTIP,   PLANEI(1,1,2),
     >                    PLANEI(1,1,3), DFACEI,    PRINT2D,
     >                    BGMODE2D(L),  FGMODE2D(L),SPMODE(L),
     >                    ITMAX2D(L),   ITFLOAT,    ITFREEZE, JLAYER,
     >                    CONV2D(L),    CONVMIN(L), DMAX2D, OMG2D(L),
     >                    POWERI(L),    POWERJ(L),  EXPI22D(L),
     >                    EXPI12D(L),   EXPJ12D(L), EXPJ22D(L),
     >                    FGLIM2D,      URFG2D,     URFLOAT)

            J = JL - 1
            DO K = 2, KTIP - 1
               X(IBU,2:J,K,1) = PLANEI(2:J,K,1)
               X(IBU,2:J,K,2) = PLANEI(2:J,K,2)
               X(IBU,2:J,K,3) = PLANEI(2:J,K,3)
            END DO

            I1 = ITU
            I2 = IBU

            DO N = 1, 2

               CALL TFI3D (IL, JL, KL, I1, I2, 1, JL, 1, KTIP,
     >                     X(1,1,1,1), X(1,1,1,2), X(1,1,1,3),
     >                     S0, BF, PRINT3D, NITSONI, TOLSONI)
               I1 = IBU
               I2 = IL
            END DO

         ELSE ! Wing-only case

            CALL TFI3D (IL, JL, KL, 1,   ITL, 1, JL, 1, KTIP,
     >                  X(1,1,1,1), X(1,1,1,2), X(1,1,1,3),
     >                  S0, BF, PRINT3D, NITSONI, TOLSONI)
         END IF

C        Below and above the wing:

         DO N = 2, 3 ! Quadrants

            I1 = IFACE(N)
            I2 = IFACE(N+1)

            CALL TFI3D (IL, JL, KL, I1, I2, 1, JL, 1, KTIP,
     >                  X(1,1,1,1), X(1,1,1,2), X(1,1,1,3),
     >                  S0, BF, PRINT3D, NITSONI, TOLSONI)
         END DO

         IF (WINGONLY) THEN ! Above the wake in one part

            CALL TFI3D (IL, JL, KL, ITU,  IL, 1, JL, 1, KTIP,
     >                  X(1,1,1,1), X(1,1,1,2), X(1,1,1,3),
     >                  S0, BF, PRINT3D, NITSONI, TOLSONI)
         END IF

C        Beyond the wing tip:

         DO N = 1, 4

            I1 = IFACE(N)
            I2 = IFACE(N+1)

            CALL TFI3D (IL, JL, KL, I1, I2, 1, JL, KTIP, KL,
     >                  X(1,1,1,1), X(1,1,1,2), X(1,1,1,3),
     >                  S0, BF, PRINT3D, NITSONI, TOLSONI)
         END DO

      ELSE ! C-H + C-O case

C        Fill the volume by halves:

         DO N = 1, 3, 2

            I1 = IFACE(N)
            I2 = IFACE(N+2)

            CALL TFI3D (IL, JL, KL, I1, I2, 1, JL, 1, KL,
     >                  X(1,1,1,1), X(1,1,1,2), X(1,1,1,3),
     >                  S0, BF, PRINT3D, NITSONI, TOLSONI)
         END DO

      END IF

      DEALLOCATE (BF)

      IF (NCALL == -3) THEN
         CALL SECOND (CPU2)
         CALL CPUTIME (CPU1, CPU2, 'GRID_TFI3D',
     >                 'to initialize the volume grid', IWRIT)
      END IF

      END SUBROUTINE GRID_TFI3D

C***********************************************************************
C
      SUBROUTINE OFFCROWN (IEDGE, IL, JL, JCR, D1, D2, TRADIAL,
     >                     XWALL, YWALL, LUNERR, FAIL)
C
C     Set up a symmetry plane edge boundary between the crown/keel at
C     the specified I = IEDGE and the outer boundary, as a straight line
C     to improve the TFI3D result.
C
C     06/10/98  DAS  Introduced for I = ITL and ITU to help TFI3D (?);
C                    usable for I = ILE as well.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN)    :: IEDGE, IL, JL, JCR, LUNERR
      REAL,    INTENT (IN)    :: D1, D2
      REAL,    INTENT (OUT)   :: TRADIAL(JL)
      REAL,    INTENT (INOUT) :: XWALL(IL,JL), YWALL(IL,JL)
      LOGICAL, INTENT (OUT)   :: FAIL

C     Local constants:

      REAL, PARAMETER :: ONE = 1., ZERO = 0.

C     Local variables:

      INTEGER  IER, J
      REAL     D1NORM, D2NORM, DRANGE, XRANGE, YRANGE

C     Execution:

      FAIL   = .FALSE.
      XRANGE = XWALL(IEDGE,JCR) - XWALL(IEDGE,JL)
      YRANGE = YWALL(IEDGE,JCR) - YWALL(IEDGE,JL)
      DRANGE = SQRT (XRANGE**2 + YRANGE**2)
      D1NORM = D1 / DRANGE
      D2NORM = D2 / DRANGE

      CALL HTDIS4 (.TRUE., ZERO, ONE, D1NORM, D2NORM, JL - JCR + 1,
     >             TRADIAL(JCR), -LUNERR, IER)

      IF (IER /= 0 .AND. IER /= 3) THEN
         WRITE (LUNERR, '(/,A,I2,A,I4,/,A,2E14.6)')
     >      ' OFFCROWN: IER from HTDIS4 =', IER, '  IEDGE:', IEDGE,
     >      ' D1NORM, D2NORM: ', D1NORM, D2NORM
         FAIL = .TRUE.
         GO TO 990
      END IF

      DO J = JCR + 1, JL - 1
         XWALL(IEDGE,J) = XWALL(IEDGE,JCR) - XRANGE * TRADIAL(J)
         YWALL(IEDGE,J) = YWALL(IEDGE,JCR) - YRANGE * TRADIAL(J)
      END DO

  990 RETURN

      END SUBROUTINE OFFCROWN

C***********************************************************************
C
      SUBROUTINE OFFLEDGE (IL, JL, KL, ILE, KP, XYZ, D1, D2,
     >                     NBLAYER, RBLAYER, T, X, Y, Z, LUNERR)
C
C     OFFLEDGE (off leading edge) modularizes setting up a composite
C     point distribution forward of a wing section grid along a line
C     which accounts for any leading edge droop. See BLGRID for more.
C
C     08/14/97  D.Saunders  Initial implementation, prompted by slats.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) :: IL, JL, KL, ILE, KP, NBLAYER, LUNERR
      REAL, INTENT (INOUT) :: XYZ(IL,JL,KL,3)
      REAL, INTENT (IN)    :: D1, D2, RBLAYER
      REAL, INTENT (OUT)   :: T(JL), X(JL), Y(JL), Z(JL)

C     Local constants:

      REAL,      PARAMETER :: HALF = 0.5
      LOGICAL,   PARAMETER :: NEW = .TRUE.
      CHARACTER, PARAMETER :: LOOSE * 1 = 'B' ! New data for each LCSFIT fit

C     Local variables:

      INTEGER   I, IER, J, K
      REAL      DX2, DY2, TTOTAL, TLINE(6), XLINE(6), YLINE(6), ZLINE(6)

C     Execution:

      K = KP
      DO I = 1, 3
         J = 3 - I
         XLINE(I) = (XYZ(ILE+J,1,K,1) + XYZ(ILE-J,1,K,1)) * HALF
         YLINE(I) = (XYZ(ILE+J,1,K,2) + XYZ(ILE-J,1,K,2)) * HALF
         ZLINE(I) =  XYZ(ILE+J,1,K,3)
      END DO

      DX2 = XLINE(1) - XLINE(3)
      DY2 = YLINE(1) - YLINE(3)
      XLINE(4) = XLINE(3) - DX2
      YLINE(4) = YLINE(3) - DY2
      ZLINE(4) = ZLINE(3)
      XLINE(6) = XYZ(ILE,JL,K,1)
      XLINE(5) = XLINE(6) + D2
      YLINE(5) = XYZ(ILE,JL,K,2)
      YLINE(6) = YLINE(5)
      ZLINE(5) = XYZ(ILE,JL,K,3)
      ZLINE(6) = ZLINE(5)

      CALL CHORDS3D (6, XLINE, YLINE, ZLINE, .FALSE., TTOTAL, TLINE)

C     Composite geometric + Vinokur arc-length distribution:

      T(1) = TLINE(3)
      T(JL) = TTOTAL

      CALL BLGRID (JL, D1, D2, NBLAYER, RBLAYER, T, LUNERR, IER)

      IF (IER /= 0 .AND. IER /= 3) THEN
         WRITE (LUNERR, '(/,A,I4,/,A,2I6,1P,3E12.3)')
     >      ' OFFLEDGE: K =', K, ' IER, JL, D1, D2, T1, T2 =',
     >      IER, JL, D1, D2, TLINE(3), TTOTAL
         STOP
      END IF

C     Interpolate X, Y, Z as functions of the distributed T:

      CALL LCSFIT (6, TLINE, XLINE, NEW, LOOSE, JL, T, X, X)
      CALL LCSFIT (6, TLINE, YLINE, NEW, LOOSE, JL, T, Y, Y)
      CALL LCSFIT (6, TLINE, ZLINE, NEW, LOOSE, JL, T, Z, Z)

      J = JL - 1
      XYZ(ILE,2:J,K,1) = X(2:J)
      XYZ(ILE,2:J,K,2) = Y(2:J)
      XYZ(ILE,2:J,K,3) = Z(2:J)

      END SUBROUTINE OFFLEDGE

C*******************************************************************************
C
      SUBROUTINE OFFTEDGE (ITEDGE, KPLANE, D1TRAIL, D2JL,
     >                     NBLAYER, RBLAYER, TGRID, LUNERR)
C
C     OFFTEDGE sets up an artificial grid boundary line from the trailing edge
C     to the outer C boundary in the plane K = 1, Ktip, or Kmax, as part of
C     using an 8-subblock CH scheme rather than 4 or 6 subblocks in order to
C     provide radial spacing control at the trailing edge that is independent
C     of the spacings at the leading edge and far downstream boundary.
C
C     Normally, a 5-point spline is established orthogonal to the trailing edge
C     and to the outer boundary, and this is evaluated at points forming a
C     geometric distribution along the arc out to point NBLAYER, followed by a
C     Vinokur distribution beyond that.  For rounded trailing edges (NNTRAIL=1)
C     the 5-pt. spline bisects the angle between the trailing edge and the wake.
C
C     OFFTEDGE serves for either above or below the wing (if awkwardly). Making
C     it completely argument driven would have been cumbersome.
C
C     05/09/96  D.Saunders  Initial implementation.
C     08/12/97      "       Adapted to use BLGRID as in 2-D application, CGRID.
C     08/09/99      "       Revised to allow rounded trailing edges.
C
C*******************************************************************************

C     Global variables:

      USE LIM
      USE LUNS
      USE XYZ

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) :: ITEDGE   ! ITL or ITU
      INTEGER, INTENT (IN) :: KPLANE   ! 1, KTIP, or KL
      REAL,    INTENT (IN) :: D1TRAIL  ! Desired radial inc. at trailing edge
      REAL,    INTENT (IN) :: D2JL     ! ... and at the outer boundary
      INTEGER, INTENT (IN) :: NBLAYER  ! Points 1:NBLAYER vary geometrically
      REAL,    INTENT (IN) :: RBLAYER  ! (or uniformly if RBLAYER = 1.0)
      REAL,    INTENT (OUT):: TGRID(*) ! Work-space for a JL-sized 1-D distribn.
      INTEGER, INTENT (IN) :: LUNERR   ! For error messages

C     Local constants:

      REAL, PARAMETER ::
     >   FRACTION = 0.05, ! Of chord, for normal distance off TE
     >   HALF = 0.5, ONE = 1., ZERO = 0., SUPPRESS = -999. ! Derivatives
      LOGICAL, PARAMETER ::
     >   CLOSED = .FALSE.
      CHARACTER, PARAMETER ::
     >   METHOD = 'B'

C     Local variables:

      INTEGER
     >   IER, I1, I2, I3, I4, J, JEVAL, K
      REAL
     >   ARROW, CFRACTION, D2, D3, DX, DY, R, TTOTAL,
     >   X2, X3, X4, X5, XT, Y2, Y3, Y4, Y5, YT,
     >   TLINE(5), XLINE(5), YLINE(5), ZLINE(5)
      LOGICAL
     >   NEW

C     Execution:

      I1 = ITEDGE       ! I2 = neighboring pt.
      I3 = 2 * ILE - I1 ! I3 = other TE pt.; I4 = its neighbor

      IF (I1 < ILE) THEN ! Below wing
         I2 = I1 + 1
         I4 = I3 - 1
         ARROW = -ONE
      ELSE               ! Above
         I2 = I1 - 1
         I4 = I3 + 1
         ARROW = ONE
      END IF

      K = KPLANE

C     Construct a 5-point curve to be splined between TE and outer C boundary:

      XLINE(1) = X(I1,1,K,1) ! Current TE pt.
      YLINE(1) = X(I1,1,K,2)
      ZLINE(1) = X(I1,1,K,3)

      XT = (XLINE(1) + X(I3,1,K,1)) * HALF    ! Middle of TE
      YT = (YLINE(1) + X(I3,1,K,2)) * HALF
      DX = XT - X(ILE,1,K,1)
      DY = YT - X(ILE,1,K,2)
      CFRACTION = SQRT (DX*DX + DY*DY) * FRACTION ! Reasonable dist. off surface

C     Construct an isosceles triangle to enable bisecting the angle between the
C     wake (extrapolated meanline) and the current surface at the trailing edge:

      X2 = X(I2,1,K,1)                        ! One pt. upstream of TE
      Y2 = X(I2,1,K,2)
      D2 = SQRT ((XT - X2)**2 + (YT - Y2)**2) ! Dist. to extrap. along meanline

      X3 = (X2 + X(I4,1,K,1)) * HALF          ! On meanline inside TE
      Y3 = (Y2 + X(I4,1,K,2)) * HALF
      D3 = SQRT ((XT - X3)**2 + (YT - Y3)**2) ! Segment to extrap. into wake

      R  = D2 / D3
      X4 = XT + R * (XT - X3) ! Pt. along wake line
      Y4 = YT + R * (YT - Y3)

      X5 = (X2 + X4) * HALF   ! Pt. on bisector of exterior wake/surface angle
      Y5 = (Y2 + Y4) * HALF

      D3 = SQRT ((XT - X5)**2 + (YT - Y5)**2) ! Seg. to extrap. along bisector

      IF (D3 > CFRACTION * 1.E-5) THEN
         R  = (CFRACTION - D3) / D3
         XLINE(2) = X5 + R * (X5 - XT)      ! Pt. along bisector
         YLINE(2) = Y5 + R * (Y5 - YT)
      ELSE                                  ! Very thin, or degenerate at KL
         DY = CFRACTION * ARROW             ! Orthogonal by original method
         R  = (YT - Y2) / (XT - X2)         ! using similar right triangles
         DX = -DY * R
         XLINE(2) = XT + DX
         YLINE(2) = YT + DY
      END IF

      XLINE(2) = XLINE(2) + (XLINE(1) - XT) ! Shift from the middle of the TE to
      YLINE(2) = YLINE(2) + (YLINE(1) - YT) ! the actual TE pt.
      ZLINE(2) = ZLINE(1)

      ZLINE(3) = ZLINE(1)
      XLINE(3) = XLINE(2) + (XLINE(2) - XLINE(1)) * HALF
      YLINE(3) = YLINE(2) + (YLINE(2) - YLINE(1))

      XLINE(5) = X(I1,JL,K,1)
      YLINE(5) = X(I1,JL,K,2)
      ZLINE(5) = X(I1,JL,K,3)

      XLINE(4) = XLINE(5)
      YLINE(4) = YLINE(5) - D2JL * ARROW
      ZLINE(4) = ZLINE(5)

C     Unnormalized arc lengths:

      CALL CHORDS3D (5, XLINE, YLINE, ZLINE, .FALSE., TTOTAL, TLINE)

C     Unnormalized geometric + Vinokur distribution:

      TGRID(1)  = ZERO ! 1 & N are reqd. in place by BLGRID
      TGRID(JL) = TTOTAL

      CALL BLGRID (JL, D1TRAIL, D2JL, NBLAYER, RBLAYER, TGRID, LUNERR,
     >             IER)
      IF (IER /= 0 .AND. IER /= 3) THEN
         WRITE (LUNERR,*) 'OFFTEDGE: IER, N, D1, D2, TTOTAL = ',
     >                    IER, JL, D1TRAIL, D2JL, TTOTAL
         STOP
      END IF

C     Evaluate X, Y, Z along the arc:

      NEW = .TRUE.
      JEVAL = 1

      DO J = 2, JL - 1

         CALL PLSCRV3D (5, XLINE, YLINE, ZLINE, TLINE, METHOD, NEW,
     >                  CLOSED, TGRID(J), JEVAL, X(I1,J,K,1),
     >                  X(I1,J,K,2), X(I1,J,K,3), SUPPRESS)
         NEW = .FALSE.
      END DO

      END SUBROUTINE OFFTEDGE

C***********************************************************************
C
      SUBROUTINE OUTER (FULLGRID)
C
C     OUTER determines (X,Y) (chordwise & vertical) grid distributions
C     on the top, forward, bottom, and downstream boundaries of a
C     typical C-mesh in the range of the wing/body.  The outer boundary
C     "C" is rectangular apart from circular upstream corners.
C     It may be skewed to match the wing sweep.
C     This version allows for a C-O mesh off the tip rather than C-H by
C     setting up the back planes for I = 1 & IL as well in both cases.
C
C***********************************************************************

C     Global variables:

      USE CONSTS
      USE DIMENS
      USE DMAX
      USE ELLIPTIC
      USE FACES
      USE FUSEL
      USE LIM
      USE LOGIC
      USE LUNS
      USE OPTIONS
      USE OUTR
      USE SPACING
      USE TFI
      USE WNGSRF
      USE XYZ
      USE XYZ0

      IMPLICIT NONE

C     Arguments:

      LOGICAL, INTENT (IN) :: FULLGRID

C     Local constants:

      REAL,      PARAMETER :: HALF = 0.5, ONE = 1., ZERO = 0.
      LOGICAL,   PARAMETER :: NEW = .TRUE.
      CHARACTER, PARAMETER :: TIGHT * 1 = 'M'

C     Local variables:

      INTEGER
     >   I, IEDGE, IER, J, K, K1, K2, KEY, L, M, M1
      REAL
     >   SA(5), SOUT(IL), SREVISE(JL),
     >   ANGLE, CS, D1, D2, DJ, DM, DU, PHIB, QCIRC, R, RJCRM1, S, SMIN,
     >   STOTAL, XLPRC, YBPRC, YTMRC

C     Execution:

      IF (.NOT. FULLGRID) GO TO 100 ! Avoid the recalculation and indenting

C     Generate the outer boundary shape and distribute points along it.

      QCIRC  = HALF * PI * RADIUS
      STOTAL = 2.* (XMAX - XMIN + QCIRC) + (YMAX - YMIN) - 4.* RADIUS
      XLPRC  = XMIN + RADIUS
      YTMRC  = YMAX - RADIUS
      YBPRC  = YMIN + RADIUS

C     Original distribution: uniform in arc length.
C     Set up arc lengths at ends of each segment of the outer boundary,
C     with SA(0) = 0. understood for the bottom right corner.

      SA(1) = XMAX  - XLPRC
      SA(2) = SA(1) + QCIRC
      SA(3) = SA(2) + YTMRC - YBPRC
      SA(4) = SA(3) + QCIRC
      SA(5) = STOTAL

      XOUT(1) = XMAX
      YOUT(1) = YMIN
      KEY = 1
      DU  = STOTAL / REAL (IL - 1)

      DO I = 2, IL - 1

         S = DU * REAL (I - 1)

         IF (S > SA(KEY)) KEY = KEY + 1

         IF (KEY == 1) THEN
            XOUT(I) = XMAX - S
            YOUT(I) = YMIN

         ELSE
            CS = S - SA(KEY-1)

            IF (KEY == 2) THEN
               PHIB    = -HALF * PI - CS / RADIUS
               XOUT(I) = XLPRC + RADIUS * COS(PHIB)
               YOUT(I) = YBPRC + RADIUS * SIN(PHIB)

            ELSE IF (KEY == 3) THEN
               XOUT(I) = XMIN
               YOUT(I) = YBPRC + CS

            ELSE IF (KEY == 4) THEN
               PHIB = PI - CS / RADIUS
               XOUT(I) = XLPRC + RADIUS * COS(PHIB)
               YOUT(I) = YTMRC + RADIUS * SIN(PHIB)

            ELSE ! KEY = 5
               XOUT(I) = XLPRC + CS
               YOUT(I) = YMAX

            END IF

         END IF

      END DO

      XOUT(IL) = XMAX
      YOUT(IL) = YMAX

C     Redistribute the uniform symmetry-plane outer boundary?

      IF (OSPNOSE /= ONE .OR. OSPTAIL /= ONE) THEN

         CALL CHORDS2D (IL, XOUT, YOUT, .FALSE., STOTAL, SOUT)

         D1 = DU * OSPTAIL ! "Outer spacing, tail boundary"; try 3.
         D2 = DU * OSPNOSE ! "Outer spacing, nose boundary"; try 1.
         DM = D1 * 0.05    ! Heuristic interior increment, adjusted below

C        Set up an artificial boundary point opposite the root trailing edge:

         IEDGE = ITL
         S = XMAX - XW(ITL,1)
         CALL INTERVAL (ILE / 2, SOUT, S, ONE, IEDGE)

C        Vinokur distributions along the lower half:

         CALL HTDIS4 (.TRUE., ZERO, SOUT(IEDGE), D1, DM, ITL,
     >                ARCOUTER, -IWRIT, IER)

         IF (IER /= 0 .AND. IER /= 3) THEN
            WRITE (IWRIT,*) 'OUTER: HTDIS4 trouble (1a). IER: ', IER
            WRITE (IWRIT,*) 'RANGE, D1, DM: ', SOUT(IEDGE), D1, DM
            GO TO 900
         END IF

C        Adjust DM to avoid non-monotonic variation between TE and XMAX:

         SMIN = 1.E+9
         DO I = 2, ITL
            SMIN = MIN (ARCOUTER(I) - ARCOUTER(I-1), SMIN)
         END DO

         IF (SMIN < DM) THEN
            DM = SMIN

            CALL HTDIS4 (.TRUE., ZERO, SOUT(IEDGE), D1, DM, ITL,
     >                   ARCOUTER, -IWRIT, IER)

            IF (IER /= 0 .AND. IER /= 3) THEN
               WRITE (IWRIT,*) 'OUTER: HTDIS4 trouble (1b). IER: ', IER
               WRITE (IWRIT,*) 'RANGE, D1, DM: ', SOUT(IEDGE), D1, DM
               GO TO 900
            END IF
         END IF

         S = ARCOUTER(ITL-1)

         CALL HTDIS4 (.TRUE., S, SOUT(ILE), DM, D2, ILE - ITL + 2,
     >                ARCOUTER(ITL-1), -IWRIT, IER)

         IF (IER /= 0 .AND. IER /= 3) THEN
            WRITE (IWRIT,*) 'OUTER: HTDIS4 trouble (1c). IER: ', IER
            WRITE (IWRIT,*) 'A, B, DM, D2: ', S, SOUT(ILE), DM, D2
            GO TO 900
         END IF

         TFIWORK(1:ILE)       = XOUT(1:ILE) ! Can't redistribute in place
         TFIWORK(IL+1:IL+ILE) = YOUT(1:ILE)

         CALL LCSFIT (ILE, SOUT, TFIWORK(1),    NEW, TIGHT, ILE,
     >                ARCOUTER, XOUT, XOUT)

         CALL LCSFIT (ILE, SOUT, TFIWORK(IL+1), NEW, TIGHT, ILE,
     >                ARCOUTER, YOUT, YOUT)

C        Likewise for the upper half:

         IEDGE = IL - IEDGE + 1

         CALL HTDIS4 (.TRUE., SOUT(IEDGE), SOUT(IL), DM, D1, ITL,
     >                ARCOUTER(ITU), -IWRIT, IER)

         IF (IER /= 0 .AND. IER /= 3) THEN
            WRITE (IWRIT,*) 'OUTER: HTDIS4 trouble (2a). IER: ', IER
            WRITE (IWRIT,*) 'S1,S2,DM,D2: ', SOUT(IEDGE),SOUT(IL),DM,D1
            GO TO 900
         END IF

         S = ARCOUTER(ITU+1)

         CALL HTDIS4 (.TRUE., SOUT(ILE), S, D2, DM, ITU - ILE + 2,
     >                ARCOUTER(ILE), -IWRIT, IER)

         IF (IER /= 0 .AND. IER /= 3) THEN
            WRITE (IWRIT,*) 'OUTER: HTDIS4 trouble (2b). IER: ', IER
            WRITE (IWRIT,*) 'S1,S2,D2,DM: ', SOUT(ILE),S,D2,DM
            GO TO 900
         END IF

         TFIWORK(ILE:IL)       = XOUT(ILE:IL)
         TFIWORK(IL+ILE:IL+IL) = YOUT(ILE:IL)

         CALL LCSFIT (ILE, SOUT(ILE), TFIWORK(ILE),    NEW, TIGHT, ILE,
     >                ARCOUTER(ILE), XOUT(ILE), XOUT(ILE))

         CALL LCSFIT (ILE, SOUT(ILE), TFIWORK(IL+ILE), NEW, TIGHT, ILE,
     >                ARCOUTER(ILE), YOUT(ILE), YOUT(ILE))
      END IF

  100 CONTINUE

C     Downstream boundary on the symmetry plane:

      IF (NNTIP == 0 .AND. .NOT. WINGONLY) THEN
         YUBACK(1,1) = HALF * (YF(1,NFSTN) + YF(NJBODY,NFSTN))
      ELSE
         YUBACK(1,1) = ZERO ! Simplifies C-O grid at tip
      END IF

      YLBACK(1,1) = YUBACK(1,1)

      IF (FULLGRID) THEN
         XUBACK(1:JL,1) = XMAX ! I = 1 and IL symmetry plane lines
         XLBACK(1:JL,1) = XMAX
         YLBACK(JL,1)   = YOUT(1)
         YUBACK(JL,1)   = YOUT(IL)
         ZUBACK(1:JL,1) = ZERO
         ZLBACK(1:JL,1) = ZERO
      ELSE ! Preserve the size of the full-grid fan
         DJ = X0(1,1,1,2) - X0(1,JCR,1,2)
         YLBACK(JCR,1) = YLBACK(1,1) - DJ
         DJ = X0(IL,JCR,1,2) - X0(IL,1,1,2)
         YUBACK(JCR,1) = YUBACK(1,1) + DJ
      END IF

C     The C-H + C-O option requires the full back plane to be set up here,
C     so do it for the C-H case as well.

C     Transfer the relative spanwise grid from the trailing edge to downstream
C     (originally done in VSHEET, but needed here for back plane to be usable):

      R = (ZW(ITL,KL) - ZERO) / (ZW(ITL,KL) - ZW(ITL,1))
      DO K = 2, KL
         ZUBACK(1,K) = ZERO + R * (ZW(ITL,K) - ZW(ITL,1))
         ZLBACK(1,K) = ZUBACK(1,K)
         XUBACK(1,K) = ZUBACK(1,K) * TANSWEEP + XMAX
         XLBACK(1,K) = XUBACK(1,K)
         YUBACK(1,K) = YUBACK(1,1)
         YLBACK(1,K) = YUBACK(1,1)
      END DO

      IF (NNTIP == 0) THEN ! C-H everywhere

C        J = JL edge of back plane:

         XLBACK(JL,1:KL) = XUBACK(1,1:KL)
         XUBACK(JL,1:KL) = XUBACK(1,1:KL)
         YLBACK(JL,1:KL) = YOUT(1)
         YUBACK(JL,1:KL) = YOUT(IL)
         ZLBACK(JL,1:KL) = ZUBACK(1,1:KL)
         ZUBACK(JL,1:KL) = ZUBACK(1,1:KL)

C        Vertical edges of back plane:

         IF (WINGONLY) THEN

            M1 = 1 ! First full K plane to process

         ELSE ! K = 1: J = JCR is a corner point

            M1 = 2

            IF (FULLGRID) THEN

C              Vary dYs for 1:JCR arithmetically from D1 (root) to D2 (crown):

               D1 = DRADIAL(1,3)
               D2 = D1CROWN
               RJCRM1 = ONE / REAL (JCR - 1)

               DO J = 2, JCR
                  R = REAL (J - 1) * RJCRM1
                  DJ = (ONE - R) * D1 + R * D2
                  YUBACK(J,1) = YUBACK(J-1,1) + DJ
                  YLBACK(J,1) = YLBACK(J-1,1) - DJ
               END DO

C              Vinokur distributions from JCR to JL, overlapped for JCR-1:JCR.

               J = JCR - 1
               CALL VINOKUR (J, JL, D2, D2JL, YUBACK(1,1), IWRIT, IER)
               CALL VINOKUR (J, JL, D2, D2JL, YLBACK(1,1), IWRIT, IER)

            ELSE

               YRADIAL = X0(1,1:JL,1,2)
               ZRADIAL = ZERO

               CALL NULINE2D (1,  JCR, YRADIAL, ZRADIAL,
     >                        YLBACK(1,1), ZLBACK(1,1))
               CALL NULINE2D (JCR, JL, YRADIAL, ZRADIAL,
     >                        YLBACK(1,1), ZLBACK(1,1))

               YRADIAL = X0(IL,1:JL,1,2)

               CALL NULINE2D (1,  JCR, YRADIAL, ZRADIAL,
     >                        YUBACK(1,1), ZUBACK(1,1))
               CALL NULINE2D (JCR, JL, YRADIAL, ZRADIAL,
     >                        YUBACK(1,1), ZUBACK(1,1))

            END IF

         END IF

         DO M = M1, 3 ! K = [1,] KTIP, & KL

            K = KFACE(M)

            XUBACK(2:JL,K) = XUBACK(1,K)
            XLBACK(2:JL,K) = XLBACK(1,K)

            IF (FULLGRID) THEN

               D1 = DRADIAL(K,3)

               CALL VINOKUR (1, JL, D1, D2JL, YUBACK(1,K), IWRIT, IER)
               CALL VINOKUR (1, JL, D1, D2JL, YLBACK(1,K), IWRIT, IER)

               ZUBACK(2:JL,K) = ZUBACK(1,K)
               ZLBACK(2:JL,K) = ZLBACK(1,K)

            ELSE

               YRADIAL = X0(1,1:JL,K,2)
               ZRADIAL = X0(1,1:JL,K,3)

               CALL NULINE2D (1, JL, YRADIAL, ZRADIAL,
     >                        YLBACK(1,K), ZLBACK(1,K))

               YRADIAL = X0(IL,1:JL,K,2)
               ZRADIAL = X0(IL,1:JL,K,3)

               CALL NULINE2D (1, JL, YRADIAL, ZRADIAL,
     >                        YUBACK(1,K), ZUBACK(1,K))

            END IF

         END DO


         IF (FULLGRID) THEN

C           Fill the interior of the back plane and (maybe) smooth it:

            DO M = 1, 2
               K1 = KFACE(M)
               K2 = KFACE(M+1)

               CALL TFIQ3D (JL, 1, JL, K1, K2, XUBACK, YUBACK, ZUBACK,
     >                      TFIWORK)
            END DO

            L = 11

            IF (ITMAX2D(L) > 0) THEN ! Split it at KTIP if we ever smooth

               CALL ELLIPQ3D (JL, KL, 1, JL, 1, KL,    XUBACK,
     >                        YUBACK,      ZUBACK,     DFACEI, PRINT2D,
     >                        BGMODE2D(L), FGMODE2D(L),SPMODE(L),
     >                        ITMAX2D(L),  ITFLOAT,    ITFREEZE, JLAYER,
     >                        CONV2D(L),   CONVMIN(L), DMAX2D, OMG2D(L),
     >                        POWERI(L),   POWERJ(L),  EXPI12D(L),
     >                        EXPI22D(L),  EXPJ12D(L), EXPJ22D(L),
     >                        FGLIM2D,     URFG2D,     URFLOAT)
            END IF

            DO M = 1, 2
               K1 = KFACE(M)
               K2 = KFACE(M+1)

               CALL TFIQ3D (JL, 1, JL, K1, K2, XLBACK, YLBACK, ZLBACK,
     >                   TFIWORK)
            END DO

            IF (ITMAX2D(L) > 0) THEN

               CALL ELLIPQ3D (JL, KL, 1, JL, 1, KL,    XLBACK,
     >                        ZLBACK,      YLBACK,     DFACEI,  PRINT2D,
     >                        BGMODE2D(L), FGMODE2D(L),SPMODE(L),
     >                        ITMAX2D(L),  ITFLOAT,    ITFREEZE,JLAYER,
     >                        CONV2D(L),   CONVMIN(L), DMAX2D, OMG2D(L),
     >                        POWERI(L),   POWERJ(L),  EXPI12D(L),
     >                        EXPI22D(L),  EXPJ12D(L), EXPJ22D(L),
     >                        FGLIM2D,     URFG2D,     URFLOAT)
            END IF

         ELSE ! Warping has to be done in the big grid because of S0:

            X(1,1:JL,1:KL,1) = XLBACK ! Whole face saves fussing with edges
            X(1,1:JL,1:KL,2) = YLBACK
            X(1,1:JL,1:KL,3) = ZLBACK

            DO M = 1, 2
               K1 = KFACE(M)
               K2 = KFACE(M+1)

               CALL WARPQ3D (1, IL, 1, JL, 1, KL, 1, 1, 1, JL, K1, K2,
     >                       X0(1,1,1,1), X0(1,1,1,2), X0(1,1,1,3), S0,
     >                       DFACEI,     DFACEJ,     DFACEK,
     >                       X(1,1,1,1), X(1,1,1,2), X(1,1,1,3))
            END DO

C           Copy back to be compatible with the full-grid & C-O cases:

            XLBACK = X(1,1:JL,1:KL,1)
            YLBACK = X(1,1:JL,1:KL,2)
            ZLBACK = X(1,1:JL,1:KL,3)

            X(IL,1:JL,1:KL,1) = XUBACK
            X(IL,1:JL,1:KL,2) = YUBACK
            X(IL,1:JL,1:KL,3) = ZUBACK

            DO M = 1, 2
               K1 = KFACE(M)
               K2 = KFACE(M+1)

               CALL WARPQ3D (1, IL, 1, JL, 1, KL, IL, IL, 1, JL, K1, K2,
     >                       X0(1,1,1,1), X0(1,1,1,2), X0(1,1,1,3), S0,
     >                       DFACEI,     DFACEJ,     DFACEK,
     >                       X(1,1,1,1), X(1,1,1,2), X(1,1,1,3))
            END DO

            XUBACK = X(IL,1:JL,1:KL,1)
            YUBACK = X(IL,1:JL,1:KL,2)
            ZUBACK = X(IL,1:JL,1:KL,3)

         END IF

      ELSE ! C-H + C-O case: back plane is more involved

C ***    Do the edges for the warp case more carefully if this ever gets used.

         D1 = DRADIAL(1,3)

         CALL VINOKUR (1, JL, D1, D2JL, YLBACK(1,1), IWRIT, IER)
         CALL VINOKUR (1, JL, D1, D2JL, YUBACK(1,1), IWRIT, IER)

         S = ZW(ILE,KTIP) ! Semispan
         R = YOUT(IL)     ! Radius of O boundary beyond tip, center (S, 0.)
         STOTAL = S + HALF * PI * R
         DU = STOTAL / REAL (KL - 1)
         D1 = DU * OSPSYM ! "Outer spacing, symmetry plane"
         D2 = DU * OSPTIP ! "Outer spacing, beyond tip"
         IF (D1 > HALF*S) THEN
            D1 = HALF*S
            WRITE (IWRIT,*) 'OUTER WARNING: D1 reduced to half semispan'
         END IF

         CALL HTDIS4 (.TRUE., ZERO, STOTAL, D1, D2, KL, ARCOUTER,
     >                -IWRIT, IER)

         IF (IER /= 0 .AND. IER /= 3) THEN
            WRITE (IWRIT,*) 'OUTER: HTDIS4 trouble (3). IER: ', IER
            WRITE (IWRIT,*) 'STOTAL, D1, D2: ', STOTAL, D1, D2
            GO TO 900
         END IF

C        ARCOUTER(*) & outermost K within the tip span are reqd. by GRID_FACES:

         KTOUTER = 2
         CALL INTERVAL (KL, ARCOUTER, S, ONE, KTOUTER)

         DO K = 2, KTOUTER
            ZUBACK(JL,K) = ARCOUTER(K)
            XUBACK(JL,K) = ARCOUTER(K) * TANSWEEP + XMAX
            YUBACK(JL,K) = YUBACK(JL,1)
         END DO

         DO K = KTOUTER + 1, KL - 1
            ANGLE = (STOTAL - ARCOUTER(K)) / R
            ZUBACK(JL,K) = R * COS (ANGLE) + S
            XUBACK(JL,K) = ZUBACK(JL,K) * TANSWEEP + XMAX
            YUBACK(JL,K) = R * SIN (ANGLE)
         END DO

C        K = KL edge is now horizontal:

         DO J = 2, JL
            ZUBACK(J,KL) = S + YUBACK(J,1)
            XUBACK(J,KL) = ZUBACK(J,KL) * TANSWEEP + XMAX
            YUBACK(J,KL) = ZERO
         END DO

C        Fill the interior of the upper back plane:

         CALL TFI2D (JL, 1, JL, 1, KL, YUBACK, ZUBACK, TFIWORK)

C        May need to smooth it too?

C        Enforce initial increments by radial redistribution:

         DO K = 2, KL

C           Arc length increments: reuse SOUT(*).

            CALL CHORDS2D (JL, YUBACK(1,K), ZUBACK(1,K), .FALSE.,
     >                     STOTAL, SOUT)

            D1 = DRADIAL(K,3)
            D2 = STOTAL - SOUT(JL-1)

C           Vinokur distribution:

            CALL HTDIS4 (.TRUE., ZERO, STOTAL, D1, D2, JL, SREVISE,
     >                   -IWRIT, IER)

            IF (IER /= 0 .AND. IER /= 3) THEN
               WRITE (IWRIT,*) 'OUTER: HTDIS4 trouble. IER: ', IER
               WRITE (IWRIT,*) 'Back upper K: ', K
               WRITE (IWRIT,*) 'D1, D2: ', D1, D2
               GO TO 900
            END IF

C           Interpolate Y, Z at these arc-lengths.  Copy needed to overwrite.

            TFIWORK(1:JL)       = YUBACK(1:JL,K)
            TFIWORK(JL+1:JL+JL) = ZUBACK(1:JL,K)

            CALL LCSFIT (JL, SOUT, TFIWORK(1),    NEW, 'B', JL,
     >                   SREVISE, YUBACK(1,K), YUBACK(1,K))

            CALL LCSFIT (JL, SOUT, TFIWORK(JL+1), NEW, 'B', JL,
     >                   SREVISE, ZUBACK(1,K), ZUBACK(1,K))

         END DO

C        Interior Xs:

         DO K = 2, KL - 1
            DO J = 2, JL - 1
               XUBACK(J,K) = ZUBACK(J,K) * TANSWEEP + XMAX
            END DO
         END DO

C        Reflect to the lower back plane:

         DO K = 1, KL
            DO J = 1, JL
               XLBACK(J,K) = XUBACK(J,K)
               YLBACK(J,K) =-YUBACK(J,K)
               ZLBACK(J,K) = ZUBACK(J,K)
            END DO
         END DO
      END IF

      GO TO 999

  900 STOP

  999 RETURN

      END SUBROUTINE OUTER

C*******************************************************************************
C
      SUBROUTINE RADIAL_BLEND (N, SA, SB, SC, WEIGHT, D1, NBL, RBL, TC)
C
C     RADIAL_BLEND performs 1-D blending of two arc-length distributions of the
C     type from BLGRID, where the first NBL points are precisely controlled in
C     real space.  Here, given point distributions at A and B and a weight for
C     blending them at some in-between location C, the first NBL points are
C     constructed as in BLGRID using D1 and RBL.  Points NBL + 1 : N are then
C     derived as a combination of the corresponding points from locations A & B.
C     The interim distribution SC(*) at C is also input for its total length;
C     it is needed externally for X/Y/Z redistribution at output points TC(*).
C
C     The intended application is to radial redistribution of block faces in a
C     manner consistent with TFI2D-type redistribution of block volumes.
C     See RADIAL_TFI2D for the case of redistributing block volume lines.
C
C     07/22/98  D.A.Saunders  Initial implementation.
C
C*******************************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN)  :: N       ! # points along each radial line
      REAL,    INTENT (IN)  :: SA(N),  ! Real-space arcs along face edges A & B
     >                         SB(N),
     >                         SC(N),  ! ... and interim arcs at location C
     >                         WEIGHT, ! Fractional distance of C between A & B
     >                         D1      ! Initial increment applied at C
      INTEGER, INTENT (IN)  :: NBL     ! # geometrically varying points at C > 0
      REAL,    INTENT (IN)  :: RBL     ! Geometric ratio applied at C
      REAL,    INTENT (OUT) :: TC(N)   ! Redistributed arcs along radial line C

C     Local constants:

      REAL, PARAMETER :: ONE = 1.

C     Local variables:

      INTEGER  I, NB
      REAL     DS, RANGEA, RANGEB, RANGEC, RI

C     Execution:

      DS = D1
      TC(1) = SC(1) ! Probably 0, but may be nonzero off a body crown line

      DO I = 2, NBL
         TC(I) = TC(I-1) + DS
         DS = DS * RBL
      END DO

      RANGEA = (ONE - WEIGHT) / (SA(N) - SA(NBL))
      RANGEB =        WEIGHT  / (SB(N) - SB(NBL))
      RANGEC =                   SC(N) - TC(NBL)

      DO I = NBL + 1, N - 1
         RI = RANGEA * (SA(I) - SA(NBL)) + RANGEB * (SB(I) - SB(NBL))
         TC(I) = TC(NBL) + RI * RANGEC
      END DO

      TC(N) = SC(N)

      END SUBROUTINE RADIAL_BLEND

C***********************************************************************
C
      SUBROUTINE RADIAL_INIT (MODE)
C
C     Initialize both sets of radial controls.  The READGRID case
C     requires some set-up normally done by WSURF.
C
C***********************************************************************

C     Global variables:

      USE LIM
      USE OPTIONS
      USE SPACING
      USE WNGSRF
      USE XYZ

      IMPLICIT NONE

C     Argument:

      INTEGER, INTENT (IN) :: MODE ! 0 means WSURF not called - use X(*,*,*,)

C     Local constants:

      REAL, PARAMETER :: ONE = 1.

C     Local variables:

      INTEGER K, N
      REAL    C, DZ, R

C     Execution:

      IF (MODE .EQ. 0) THEN ! READGRID case, else WSURF sets these

         DO K = 1, KL
            CHORDM(K) = X(ITL,1,K,1) - X(ILE,1,K,1)
            ZWLE(K)   = X(ILE,1,K,3)
         END DO

      END IF

      C = CHORDM(1)
      D1CROWNEU = C * D1CROWNEU
      D1NOSEEU  = C * D1NOSEEU
      D1TAILEU  = C * D1TAILEU
      D2JLEU    = C * D2JLEU
      DWTAILEU  = C * DWTAILEU
      DCTAILEU  = C * DCTAILEU

      D1CROWNNS = C * D1CROWNNS
      D1NOSENS  = C * D1NOSENS
      D1TAILNS  = C * D1TAILNS
      D2JLNS    = C * D2JLNS
      DWTAILNS  = C * DWTAILNS
      DCTAILNS  = C * DCTAILNS

      DZ = ONE / (ZWLE(KTIP) - ZWLE(1))

      DO K = 1, KTIP
         C = CHORDM(K)
         R = (ZWLE(K) - ZWLE(1)) * DZ

         DRADIALEU(K,1) = ((ONE - R) * D1ROOTEU + R * D1TIPEU) * C
         DRADIALEU(K,2) = ((ONE - R) * D2ROOTEU + R * D2TIPEU) * C
         DRADIALEU(K,3) = ((ONE - R) * D3ROOTEU + R * D3TIPEU) * C

         DRADIALNS(K,1) = ((ONE - R) * D1ROOTNS + R * D1TIPNS) * C
         DRADIALNS(K,2) = ((ONE - R) * D2ROOTNS + R * D2TIPNS) * C
         DRADIALNS(K,3) = ((ONE - R) * D3ROOTNS + R * D3TIPNS) * C
      END DO

      IF (NNTIP == 0) THEN ! 1-sided stretching beyond the tip

         DZ = ONE / (ZWLE(KL) - ZWLE(KTIP))

         DO K = KTIP + 1, KL
            R = (ZWLE(K) - ZWLE(KTIP)) * DZ
            DRADIALEU(K,1) = ((ONE - R) * D1TIPEU + R * D1KMAXEU) * C
            DRADIALEU(K,2) = ((ONE - R) * D2TIPEU + R * D2KMAXEU) * C
            DRADIALEU(K,3) = ((ONE - R) * D3TIPEU + R * D3KMAXEU) * C

            DRADIALNS(K,1) = ((ONE - R) * D1TIPNS + R * D1KMAXNS) * C
            DRADIALNS(K,2) = ((ONE - R) * D2TIPNS + R * D2KMAXNS) * C
            DRADIALNS(K,3) = ((ONE - R) * D3TIPNS + R * D3KMAXNS) * C
         END DO

      ELSE ! C-O grid beyond tip

         DO N = 1, 3
            DRADIALEU(KTIP:KL,N) = DRADIALEU(KTIP,N)
            DRADIALNS(KTIP:KL,N) = DRADIALNS(KTIP,N)
         END DO

      END IF

      END SUBROUTINE RADIAL_INIT

C***********************************************************************
C
      SUBROUTINE RADIAL_SET (EULER)
C
C     Activate the indicated radial grid controls.
C
C***********************************************************************

C     Global variables:

      USE SPACING

      IMPLICIT NONE

C     Argument:

      LOGICAL, INTENT (IN) :: EULER

C     Execution:

      IF (EULER) THEN ! Initial controls (possibly N-S-type if REDISTRIBUTE = F)

         D1CROWN = D1CROWNEU
         D1NOSE  = D1NOSEEU
         D1TAIL  = D1TAILEU
         DWTAIL  = DWTAILEU
         DCTAIL  = DCTAILEU
         D2JL    = D2JLEU
         NBLAYER = NBLAYEREU
         RBLAYER = RBLAYEREU

         DRADIAL = DRADIALEU ! (1:KL,1:3) scaled by local chord now

      ELSE ! REDISTRIBUTE = T

         D1CROWN = D1CROWNNS
         D1NOSE  = D1NOSENS
         D1TAIL  = D1TAILNS
         DWTAIL  = DWTAILNS
         DCTAIL  = DCTAILNS
         D2JL    = D2JLNS
         NBLAYER = NBLAYERNS
         RBLAYER = RBLAYERNS

         DRADIAL = DRADIALNS ! (1:KL,1:3) scaled by local chord

      END IF

      END SUBROUTINE RADIAL_SET

C*******************************************************************************
C
      SUBROUTINE RADIAL_TFI2D (IQUAD, I1, I2, K1, K2)
C
C     RADIAL_TFI2D applies the TFI2D algorithm to radial grid line arc-length
C     distributions of one subblock.  First, the (relative) boundary radial
C     distributions and the normalized edge distributions for the J = 1 plane
C     are established.  Then, for each interior (I,K), a relative distribution
C     in the radial direction is derived from these (normalized) distributions
C     edge and corner distributions and applied to the actual arc length.
C     Finally, new (X,Y,Z)s are interpolated at the redistributed arcs.
C
C     Provision is made for precise BLGRID-type control of points 1:NBL in real
C     space.  In fact, the normalized edge radial distributions must allow for
C     NBL.  See RADIAL_BLEND for the analogous 1-D case.
C
C     Unnormalized radial arcs are assumed to be input in S0, along with
C     arcs in the J = 1 plane as needed at the J = 1 face edges for the TFI.
C     They must be renormalized for subfaces here.
C
C     07/27/98  D. Saunders  Adaptation of TFI2D.
C
C*******************************************************************************

C     Global variables:

      USE LIM
      USE SPACING
      USE XYZ
      USE XYZ0

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   IQUAD,         ! Quadrant # needed to set up D1(*) as fns. of I
     >   I1, I2, K1, K2 ! Subblock limits in the J = 1 plane

C     Local constants:

      REAL,      PARAMETER :: ONE = 1., SUPPRESS = -999., ZERO = 0.
      LOGICAL,   PARAMETER :: CLOSED = .FALSE.
      CHARACTER, PARAMETER :: METHOD = 'B'

C     Local variables:

      REAL, DIMENSION (I1:I2) ::
     >   D1I,                ! D1 at I for given K
     >   EDGEK1, EDGEK2      ! Normalized arcs along J = 1 edges

      REAL, DIMENSION (K1:K2) ::
     >   EDGEI1, EDGEI2      !   "   "   "

      REAL, DIMENSION (JL) ::
     >   SIK,                ! Unnormalized interim distribution at (I,K)
     >   TIK                 ! Blended radial distribution at (I,K)

      REAL, DIMENSION (JL,I1:I2) ::
     >   SK1, SK2            ! Normalized radial distributions at K faces
      REAL, DIMENSION (JL,K1:K2) ::
     >   SI1, SI2            ! ... and I faces (allowing for NBLAYER)

C     Local variables:

      INTEGER I, J, JEVAL, K, N1, N2, NBL
      REAL    DI, DK, DT, DV, P, P1, Q, Q1, R
      LOGICAL NEW

C     Execution:

C     Normalized arcs along the J = 1 edges (Cf. TFI2D's edge distributions):

      DT = ONE / (S0(I2,1,K1,1) - S0(I1,1,K1,1))
      DV = ONE / (S0(I2,1,K2,1) - S0(I1,1,K2,1))

      DO I = I1, I2
         EDGEK1(I) = (S0(I,1,K1,1) - S0(I1,1,K1,1)) * DT
         EDGEK2(I) = (S0(I,1,K2,1) - S0(I1,1,K2,1)) * DV
      END DO

      DT = ONE / (S0(I1,1,K2,3) - S0(I1,1,K1,3))
      DV = ONE / (S0(I2,1,K2,3) - S0(I2,1,K1,3))

      DO K = K1, K2
         EDGEI1(K) = (S0(I1,1,K,3) - S0(I1,1,K1,3)) * DT
         EDGEI2(K) = (S0(I2,1,K,3) - S0(I2,1,K1,3)) * DV
      END DO

C     Normalize the face radial distributions only once per block.
C     For each J, these are analogous to the edge Xs and Ys in TFI2D -
C     they are what are being interpolated into the interior of the J = 1 face.

      IF (IQUAD == 1 .OR. IQUAD == 4) THEN
         NBL = 1 ! The back plane doesn't have a boundary-layer-type grid
      ELSE
         NBL = NBLAYER
      END IF

      DO I = I1, I2
         DT = ONE / (S0(I,JL,K1,2) - S0(I,NBL,K1,2))
         DV = ONE / (S0(I,JL,K2,2) - S0(I,NBL,K2,2))
         DO J = NBL + 1, JL - 1
            SK1(J,I) = (S0(I,J,K1,2) - S0(I,NBL,K1,2)) * DT
            SK2(J,I) = (S0(I,J,K2,2) - S0(I,NBL,K2,2)) * DV
         END DO
      END DO

      DO K = K1, K2
         DT = ONE / (S0(I1,JL,K,2) - S0(I1,NBL,K,2))
         DV = ONE / (S0(I2,JL,K,2) - S0(I2,NBL,K,2))
         DO J = NBL + 1, JL - 1
            SI1(J,K) = (S0(I1,J,K,2) - S0(I1,NBL,K,2)) * DT
            SI2(J,K) = (S0(I2,J,K,2) - S0(I2,NBL,K,2)) * DV
         END DO
      END DO

      N1 = ABS (IQUAD - 3) + 1 ! Radial increment pointer for lower I face
      N2 = ABS (IQUAD - 2) + 1 ! ... and upper I face

      DO K = K1 + 1, K2 - 1

C        Set up initial increments as functions of I for this K:

         DT = ONE / (S0(I2,1,K,1) - S0(I1,1,K,1))

         DO I = I1, I2
            R = (S0(I,1,K,1) - S0(I1,1,K,1)) * DT
            D1I(I) = (ONE - R) * DRADIAL(K,N1) + R * DRADIAL(K,N2)
         END DO

C        Process interior radial lines:

         DO I = I1 + 1, I2 - 1

            SIK = S0(I,1:JL,K,2) ! Interim Euler-type arcs

C           Option for geometric spacing in the boundary layer:

            DT = D1I(I)
            TIK(1) = SIK(1) ! Probably zero

            DO J = 2, NBL
               TIK(J) = TIK(J-1) + DT
               DT = DT * RBLAYER
            END DO

            DT = SIK(JL) - TIK(NBL) ! Length beyond b.layer to be redistrbd.

C           Soni blending functions are derived as the (P,Q) intersection
C           coordinates of the straight lines connecting opposite I and K
C           points in the normalized edge distributions:

            DI = EDGEK2(I) - EDGEK1(I)
            DK = EDGEI2(K) - EDGEI1(K)
            DV = ONE / (ONE  - DI * DK)
            P  = (EDGEK1(I) + EDGEI1(K) * DI) * DV
            Q  = (EDGEI1(K) + EDGEK1(I) * DK) * DV
            P1 = ONE - P
            Q1 = ONE - Q

C           2-D TFI interpolation of relative radial distributions from the
C           edges into the interior, for points beyond the boundary layer:

            DO J = NBL + 1, JL - 1
               TIK(J) = TIK(NBL) + DT *
     >                  (P1 * SI1(J,K) + P * SI2(J,K)
     >                 + Q1 * SK1(J,I) + Q * SK2(J,I)
     >                 - Q1 * (P1 * SK1(J,I1) + P * SK1(J,I2))
     >                 - Q  * (P1 * SK2(J,I1) + P * SK2(J,I2)))
            END DO

C****       TIK(JL) = SIK(JL) ! Not needed

C           Calculate (X,Y,Z)s at the revised Ts:

            XRADIAL = X(I,1:JL,K,1)
            YRADIAL = X(I,1:JL,K,2)
            ZRADIAL = X(I,1:JL,K,3)

            NEW = .TRUE.
            JEVAL = 1

            DO J = 2, JL - 1

               CALL PLSCRV3D (JL, XRADIAL, YRADIAL, ZRADIAL, SIK,
     >                        METHOD, NEW, CLOSED, TIK(J), JEVAL,
     >                        X(I,J,K,1), X(I,J,K,2), X(I,J,K,3),
     >                        SUPPRESS)
               NEW = .FALSE.

            END DO

         END DO ! Next I

      END DO ! Next K

      END SUBROUTINE RADIAL_TFI2D

C***********************************************************************
C
      SUBROUTINE RECTIFY (K1, K2, KL, IL, ILE, NHALF, SNORM,
     >                    XW, YW, ZW)
C
C     RECTIFY forces sections of an arc-based wing surface grid to have
C     the ILE index at their true leading edges.  Rounded leading edges
C     are assumed where it matters (near the body), but the wrap-around
C     redistribution is suppressed if the leading edge angle is less
C     than 45 degrees.  In fact, the redistribution is continued as long
C     as the leading edge does not appear sharp, for better blending.
C
C     09/29/97  DAS  Adaptation of WREGUL for patching very 3-dimensional
C                    grid sections near a fuselage.
C     06/25/98   "   Discrete choice of leading edge can add irregularity
C                    in the Y direction.  Use continuous techniques.
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) :: K1, K2       ! Subset of wing sections to treat
      INTEGER, INTENT (IN) :: KL, IL       ! Dimensions of regularized arrays
      INTEGER, INTENT (IN) :: ILE, NHALF   ! Middle index & number either side
      REAL,    INTENT (IN) :: SNORM(NHALF) ! Normalized chordwise distribution
      REAL, INTENT (INOUT) :: XW(IL,KL),   ! Grid sections
     >                        YW(IL,KL), ZW(IL,KL)

C     Local constants:

      INTEGER,   PARAMETER :: LUNERR = 6, NFUNMAX = 15
      REAL,      PARAMETER :: COSMIN = 1.414, ! ~2 * COS (PI/4)
     >                        ONE = 1.
      LOGICAL,   PARAMETER :: NEW = .TRUE., NORM = .FALSE.
      CHARACTER, PARAMETER :: LOOSE * 1 = 'B', SUBNAME * 7 = 'RECTIFY'

C     Local variables:

      INTEGER    I, ILEAD, ILEFT, ISTAT, ITL, ITU, K, NFUN, NLEAD, NWING
      REAL       DERIVS(2*NHALF), S(IL), SEVAL(2*NHALF), XYZ(IL),
     >           A, B, C, EPS, SLOWER, SMAX, SUPPER, SLE, TOL, XLE
      LOGICAL    DOIT, NEWLE

C     Execution:

      NWING = 2*NHALF - 1
      NLEAD = 5 ! 4-pt. method in one of two adjacent intervals
      ITL   = ILE - NHALF + 1
      ITU   = ILE + NHALF - 1
      EPS   = MAX (3.* EPSILON (ONE), 1.E-7)

      DO K = K1, K2

C        Locate input leading edge index:

         DO I = ILE - 12, ILE + 12
            IF (XW(I+1,K) > XW(I,K)) EXIT
         END DO
         ILEAD = I

         DOIT = ILEAD /= ILE

         IF (.NOT. DOIT) THEN ! Do it anyway unless it's a sharp section
            A = (XW(I-1,K) - XW(I,K)) ** 2 + (YW(I-1,K) - YW(I,K)) ** 2
            B = (XW(I+1,K) - XW(I,K)) ** 2 + (YW(I+1,K) - YW(I,K)) ** 2
            C = (XW(I-1,K) - XW(I+1,K))**2 + (YW(I-1,K) - YW(I+1,K))**2
            DOIT = A + B - C < COSMIN * SQRT (A * B)
         END IF

         IF (DOIT) THEN

C           Apply the rounded leading edge case of WREGUL, but make it vary
C           smoothly across Ks by not just picking the nearest data point:

            CALL CHORDS3D (NWING, XW(ITL,K), YW(ITL,K), ZW(ITL,K), NORM,
     >                     SMAX, S(ITL))

C           Initialize an iteration for minimizing X w.r.t. S:

            ISTAT = 2
            ILEFT = ILEAD - 2
            A     = S(ILEFT)
            B     = S(ILEAD + 2)
            TOL   = B * EPS
            NEWLE = .TRUE.
            NFUN  = NFUNMAX

            DO I = 1, NFUNMAX ! Or until convergence

               CALL FMINRC (A, B, SLE, XLE, TOL, NFUN, SUBNAME, -LUNERR,
     >                      ISTAT)

               IF (ISTAT < -1) THEN

                  WRITE (LUNERR, '(/, A, 2I4)')
     >               ' RECTIFY: Bad FMINRC return. ISTAT, K =', ISTAT, K
                  STOP

               ELSE IF (ISTAT < 0) THEN

                  WRITE (LUNERR, '(/, A, I4)')
     >               ' RECTIFY: FMINRC hit NFUNMAX. K =', K

               ELSE IF (ISTAT > 0) THEN ! Evaluate X at S = SLE

                  CALL LCSFIT (NLEAD, S(ILEFT), XW(ILEFT,K), NEW, LOOSE,
     >                         1, SLE, XLE, DERIVS)

                  NEWLE = .FALSE.

               ELSE ! ISTAT = 0 - success

                  EXIT

               END IF

            END DO

C           Redistribute the section in place:

            SLOWER = SLE
            SUPPER = SMAX - SLE

            DO I = 1, NHALF
               SEVAL(I) = SLOWER * (ONE - SNORM(NHALF + 1 - I))
               SEVAL(I + NHALF - 1) = SLOWER + SUPPER * SNORM(I)
            END DO

            CALL LCSFIT (NWING, S(ITL), XW(ITL,K), NEW, LOOSE, NWING,
     >                   SEVAL, XYZ(ITL), DERIVS)

            XW(ITL:ITU,K) = XYZ(ITL:ITU)

            CALL LCSFIT (NWING, S(ITL), YW(ITL,K), NEW, LOOSE, NWING,
     >                   SEVAL, XYZ(ITL), DERIVS)

            YW(ITL:ITU,K) = XYZ(ITL:ITU)

            CALL LCSFIT (NWING, S(ITL), ZW(ITL,K), NEW, LOOSE, NWING,
     >                   SEVAL, XYZ(ITL), DERIVS)

            ZW(ITL:ITU,K) = XYZ(ITL:ITU)

        END IF

      END DO

      END SUBROUTINE RECTIFY

C*******************************************************************************
C
      SUBROUTINE REGRID_EDGES (NCALL)
C
C     Redistribution of radial edges in the subblock faces.
C     This version includes wing cranks as additional block faces.
C
C*******************************************************************************

C     Global variables:

      USE LIM
      USE LOGIC
      USE LUNS
      USE SPACING
      USE WNGSRF
      USE XYZ
      USE XYZ0

      IMPLICIT NONE

C     Argument:

      INTEGER, INTENT (IN) :: NCALL

C     Local constants:

      REAL,      PARAMETER :: ZERO = 0., SUPPRESS = -999.
      LOGICAL,   PARAMETER :: CLOSED = .FALSE.
      CHARACTER, PARAMETER :: METHOD * 1 = 'B'

C     Local variables:

      INTEGER I, IER, J, JEVAL, K, L, M, M1, N, NBL, NJ
      REAL    D1, TOTAL
      LOGICAL NEW
      REAL, DIMENSION (JL) :: TNEW

C     Execution:

      IF (NCALL == -3) THEN

C        Set up a list of K faces which includes flagged wing crank stations:

         KREGRID(1) = 1
         NKREGRID = 2

         DO K = 1, NCRANK
            IF (LCRANK(K)) THEN
               KREGRID(NKREGRID) = MCRANK(K)
               NKREGRID = NKREGRID + 1
            END IF
         END DO

         KREGRID(NKREGRID) = KTIP
         NKREGRID = NKREGRID + 1
         KREGRID(NKREGRID) = KL

      END IF


      IF (WINGONLY) THEN

         M1 = 1 ! First full K plane to process

      ELSE

         M1 = 2

C        Treat radial edges not on the body in the K = 1 plane:

         D1 = DRADIAL(1,3) ! For the end of the aft "fan"

         DO N = 1, 5, 4

            I = IFACE(N) ! 1 & IL for K = 1

            YRADIAL = X(I,1:JL,1,2)

            CALL VINOKUR (1,    JCR, D1,   D1CROWN, YRADIAL, IWRIT, IER)
            CALL VINOKUR (JCR-1, JL, D1CROWN, D2JL, YRADIAL, IWRIT, IER)

            X(I,1:JL,1,2) = YRADIAL

         END DO

         NJ = JL - JCR + 1

         ZRADIAL(JCR:JL) = ZERO
         TNEW(JCR)       = ZERO

         DO N = 2, 4

            I = IFACE(N) ! ITL, ILE, & ITU for K = 1

            XRADIAL(JCR:JL) = X(I,JCR:JL,1,1) ! Get initial Euler distribn. ...
            YRADIAL(JCR:JL) = X(I,JCR:JL,1,2)

            CALL CHORDS2D (NJ, XRADIAL(JCR), YRADIAL(JCR), .FALSE.,
     >                     TOTAL, TRADIAL(JCR))

            TNEW(JL) = TOTAL ! ... and redistribute it

            CALL VINOKUR (JCR, JL, D1CROWN, D2JL, TNEW, IWRIT, IER)

            NEW = .TRUE.
            JEVAL = 1 ! Not JCR

            DO J = JCR + 1, JL - 1

               CALL PLSCRV3D (NJ, XRADIAL(JCR), YRADIAL(JCR),
     >                        ZRADIAL(JCR), TRADIAL(JCR),
     >                        METHOD, NEW, CLOSED, TNEW(J), JEVAL,
     >                        X(I,J,1,1), X(I,J,1,2), X(I,J,1,3),
     >                        SUPPRESS) ! Derivatives
               NEW = .FALSE.

            END DO

         END DO ! Next I edge off body

      END IF

C     Edges of all quadrants in full K planes:

      TNEW(1) = ZERO

      DO M = M1, NKREGRID

         K = KREGRID(M) ! [1,] [MCRANK(*),] KTIP, & KL

         DO N = 1, 5 ! 1, ITL, ILE, ITU, & IL

            I = IFACE(N)

            XRADIAL = X(I,1:JL,K,1) ! Get the initial Euler distribution ...
            YRADIAL = X(I,1:JL,K,2)
            ZRADIAL = X(I,1:JL,K,3)

            CALL CHORDS3D (JL, XRADIAL, YRADIAL, ZRADIAL, .FALSE.,
     >                     TOTAL, TRADIAL)

            TNEW(JL) = TOTAL ! ... and redistribute it

            L  = ABS (N - 3) + 1
            D1 = DRADIAL(K,L)

            IF (L == 3) THEN ! Match backplane edges in subroutine OUTER
               NBL = 1
            ELSE
               NBL = NBLAYER
            END IF

            CALL BLGRID (JL, D1, D2JL, NBL, RBLAYER, TNEW, IWRIT, IER)

            NEW = .TRUE.
            JEVAL = 1

            DO J = 2, JL - 1

               CALL PLSCRV3D (JL, XRADIAL, YRADIAL, ZRADIAL, TRADIAL,
     >                        METHOD, NEW, CLOSED, TNEW(J), JEVAL,
     >                        X(I,J,K,1), X(I,J,K,2), X(I,J,K,3),
     >                        SUPPRESS)
               NEW = .FALSE.

            END DO

         END DO ! Next I edge in plane K

      END DO ! Next K subblock face

      END SUBROUTINE REGRID_EDGES

C*******************************************************************************
C
      SUBROUTINE REGRID_FACES (NCALL)
C
C     Redistribution of all radial lines in the subblock faces given new edges.
C
C*******************************************************************************

C     Global variables:

      USE LIM
      USE LOGIC
      USE LUNS
      USE SPACING
      USE XYZ
      USE XYZ0

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) :: NCALL

C     Local constants:

      REAL,      PARAMETER :: ONE = 1., SUPPRESS = -999.
      LOGICAL,   PARAMETER :: CLOSED = .FALSE.
      CHARACTER, PARAMETER :: METHOD * 1 = 'B'

C     Local variables:

      INTEGER
     >   I, I1, I2, IR, IY, J, K, K1, K2, L, M, M1, N, NBL
      REAL
     >   D1, D1I(IL), D2, DENOM, R, TOTAL, XRANGE, YRANGE
      REAL(KIND=4) CPU1, CPU2

C     Execution:

      IF (NCALL < 0) CALL SECOND (CPU1)


C     Non-radial arcs in the J = 1 plane serve as weights for blending
C     face edge distributions into interior face radial lines:
C     PARAMXYZ is a menace setting degenerate directions to 1, and it
C     normalizes unnecessarily here, so avoid it.

      CALL CHORDS_J1 ()


      IF (WINGONLY) THEN

         M1 = 1 ! First full K plane to do

C        Unnormalized radial arcs for new face edges and interim interior
C        radial lines of all subblock faces:

         CALL CHORDS_FACES ()

      ELSE ! K = 1 face

         M1 = 2

C        Body surface:

         CALL FIXBGRID (IL, JL, IBL, ITL, ILE, ITU, IBU, JCR,
     >                  DRADIAL(1,1), DRADIAL(1,2), D1NOSE,
     >                  DWTAIL, DCTAIL, NBLAYER, RBLAYER,
     >                  X(1,1,1,1), X(1,1,1,2), X(1,1,1,3), IWRIT)

C        Unnormalized radial arcs for new face edges and interim interior
C        radial lines of all subblock faces:

         CALL CHORDS_FACES () ! Must be done after FIXBGRID


C        Aft "fan": Awkward!

         I1 = 1
         I2 = IBL - 1 ! For off-body part of fan
         IY = IBL     ! End of body
         IR = I2      ! One right of end of body
         DENOM = ONE / (X(1,1,1,1) - X(IR,1,1,1))

         DO N = 1, 2

C           Transfer the last body section distribution as in Z0PLANE:

            CALL CHORDSRF (IL, JL, IY, 1, JCR, X(1,1,1,1), X(1,1,1,2),
     >                     X(1,1,1,3), .TRUE., TOTAL, TRADIAL)

            XRANGE = X(IR,1,1,1) - X(IR,JCR,1,1)
            YRANGE = X(IR,1,1,2) - X(IR,JCR,1,2)

            DO J = 2, JCR - 1
               X(IR,J,1,1) = X(IR,1,1,1) - XRANGE * TRADIAL(J)
               X(IR,J,1,2) = X(IR,1,1,2) - YRANGE * TRADIAL(J)
            END DO

C           Update S0(IR,2:JCR-1,1,2) for REGRID_KPLANE:

            CALL CHORDSRF (IL, JL, IR, 1, JCR, X(1,1,1,1), X(1,1,1,2),
     >                     X(1,1,1,3), .FALSE., TOTAL, TRADIAL)

            S0(IR,2:JCR-1,1,2) = TRADIAL(2:JCR-1)

            IF (N == 1) THEN
               D1 = DRADIAL(1,3)
               D2 = X(IR,1,1,2) - X(IR,2,1,2)
            ELSE
               D1 = X(IR,2,1,2) - X(IR,1,1,2)
               D2 = DRADIAL(1,3)
            END IF

C           Set up D1I(*) = D1 as a fn. of I for I = I1:I2.

            DO I = I1 + 1, I2 - 1
               R = ABS (X(I,1,1,1) - X(I1,1,1,1)) * DENOM
               D1I(I) = (ONE - R) * D1 + R * D2
            END DO

C           TFI2D should do here, but in case the fan had been smoothed,
C           use REGRID_KPLANE with NBLAYER = 1 to blend the I1/I2 distribns.

            CALL REGRID_KPLANE (1, I1, I2, 1, JCR, IL, JL, KL,
     >                          1, RBLAYER, D1I, S0, X)

            I1 = IBU + 1
            I2 = IL
            IY = IBU
            IR = I1

         END DO

C        Off-body part of the symmetry plane:

C        Arcs along the JCR line are needed for weights (J = 1 elsewhere):

         CALL CHORDS2D (IL, X(1,JCR,1,1), X(1,JCR,1,2), .FALSE., TOTAL,
     >                  S0(1,JCR,1,1))

         D1I = D1CROWN ! (1:IL) trivial case - no need for IVARIATION procedures

         DO N = 1, 4 ! Quadrants

            CALL REGRID_KPLANE (1, IFACE(N), IFACE(N+1), JCR, JL,
     >                          IL, JL, KL, 1, RBLAYER, D1I, S0, X)
         END DO

      END IF


C     Redistribute the K-plane boundaries:

      DO M = M1, NKREGRID ! [1,] [MCRANK(*),] KTIP, & KL

         K = KREGRID(M)

C        Set up initial increments D1I(1:IL) for this K:

         CALL IVARIATION (K) ! Internal procedure

         DO N = 1, 4 ! Quadrants

            IF (N == 1 .OR. N == 4) THEN
               NBL = 1 ! Simple blend, since back plane doesn't have a b.layer
            ELSE
               NBL = NBLAYER
            END IF

            CALL REGRID_KPLANE (K, IFACE(N), IFACE(N+1), 1, JL,
     >                          IL, JL, KL, NBL, RBLAYER, D1I, S0, X)
         END DO

      END DO


C     Regrid the I = 1, ITL, ILE, ITU, & IL planes in two parts by blending
C     the radial edge distributions from the bounding K planes:

      DO N = 1, 5

         I = IFACE(N)
         L = ABS (N - 3) + 1

         IF (L == 3) THEN
            NBL = 1 ! Back plane doesn't have a boundary-layer-type distribution
         ELSE
            NBL = NBLAYER
         END IF

         DO M = 2, NKREGRID

            K1 = KREGRID(M-1)
            K2 = KREGRID(M)

            CALL REGRID_IPLANE (I, K1, K2, IL, JL, KL,
     >                          NBL, RBLAYER, DRADIAL(1,L), S0, X)
         END DO

      END DO

      IF (NCALL < 0) THEN
         CALL SECOND (CPU2)
         CALL CPUTIME (CPU1, CPU2, 'REGRID_FACES',
     >                 'to redistribute subblock faces', IWRIT)
      END IF

      RETURN

      CONTAINS ! Internal procedure for REGRID_FACES

         SUBROUTINE IVARIATION (KPLANE)

!        Set up initial radial increments as functions of I for this K.

!        Argument:

         INTEGER, INTENT (IN) :: KPLANE

!        Local constants:

         REAL, PARAMETER :: ONE = 1.

!        Local variables:

         INTEGER I, I1, I2, K, N, L1, L2
         REAL    DENOM, R

!        Execution:

         K = KPLANE

!        Normalized arc lengths along the J = 1 line of this K plane:

         DO N = 1, 4 ! Intervals between I face boundaries

            L1 = ABS (N - 3) + 1 ! Radial increment pointer for lower I face
            L2 = ABS (N - 2) + 1 ! ... and upper I face
            I1 = IFACE(N)
            I2 = IFACE(N+1)      ! S0(*,1,*,1) are normalized
            DENOM = ONE / (S0(I2,1,K,1) - S0(I1,1,K,1))

            DO I = I1, I2
               R = (S0(I,1,K,1) - S0(I1,1,K,1)) * DENOM
               D1I(I) = (ONE - R) * DRADIAL(K,L1) + R * DRADIAL(K,L2)
            END DO

         END DO

         END SUBROUTINE IVARIATION

      END SUBROUTINE REGRID_FACES

C***********************************************************************
C
      SUBROUTINE REGRID_IPLANE (I, K1, K2, IL, JL, KL,
     >                          NBLAYER, RBLAYER, D1K, S0, X)
C
C     Redistribute the interior radial lines of the given I sub-plane of
C     a C grid using geometric + blended edge distributions.
C
C     D1K(*) = D1 at each K for this I.
C
C     07/08/98  DAS  Adaptation of RADIAL_KPLANE.
C     07/22/98   "   Introduced RADIAL_BLEND; ignored JCR boundary.
C     07/24/98   "   S0 input avoids arc calculations here that will be
C                    needed during volume regridding.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) :: I, K1, K2, IL, JL, KL, NBLAYER
      REAL,    INTENT (IN) :: RBLAYER, D1K(K2), S0(IL,JL,KL,3)
      REAL, INTENT (INOUT) :: X(IL,JL,KL,3)

C     Local constants:

      REAL,      PARAMETER :: ONE = 1., SUPPRESS = -999.
      LOGICAL,   PARAMETER :: CLOSED = .FALSE.
      CHARACTER, PARAMETER :: METHOD * 1 = 'B'

C     Local variables:

      INTEGER    J, JEVAL, K
      LOGICAL    NEW
      REAL       DENOM, WEIGHT
      REAL, DIMENSION (1:JL) :: TK1, TK2, TNEW, TJ, XJ, YJ, ZJ

C     Execution:

C     Unnormalized K1 and K2 edge distributions for blending:

      TK1 = S0(I,1:JL,K1,2)
      TK2 = S0(I,1:JL,K2,2)

      DENOM = ONE / (S0(I,1,K2,3) - S0(I,1,K1,3))

C     For each interior radial line ...

      DO K = K1 + 1, K2 - 1

         TJ = S0(I,1:JL,K,2) ! Unnormalized radial arcs

C        Redistribute the arc lengths via a combination of those at K1 & K2.

         WEIGHT = (S0(I,1,K,3) - S0(I,1,K1,3)) * DENOM

         CALL RADIAL_BLEND (JL, TK1, TK2, TJ, WEIGHT, D1K(K),
     >                      NBLAYER, RBLAYER, TNEW)

C        Redistribute X, Y, Z:

         XJ = X(I,1:JL,K,1)
         YJ = X(I,1:JL,K,2)
         ZJ = X(I,1:JL,K,3)

         NEW = .TRUE.
         JEVAL = 1

         DO J = 2, JL - 1

            CALL PLSCRV3D (JL, XJ, YJ, ZJ, TJ, METHOD, NEW, CLOSED,
     >                     TNEW(J), JEVAL, X(I,J,K,1), X(I,J,K,2),
     >                     X(I,J,K,3), SUPPRESS) ! Derivatives
            NEW = .FALSE.

         END DO

      END DO ! Next K

      END SUBROUTINE REGRID_IPLANE

C***********************************************************************
C
      SUBROUTINE REGRID_KPLANE (K, I1, I2, J1, J2, IL, JL, KL,
     >                          NBLAYER, RBLAYER, D1I, S0, X)
C
C     Redistribute the interior radial lines of the given K sub-plane of
C     a C grid using geometric + blended edge distributions.
C
C     D1I(*) = D1 at each I.
C
C     07/22/98  DAS  Analogue of revised REGRID_IPLANE.
C     07/24/98   "   S0 input avoids arc calculations here that will be
C                    needed during volume regridding.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) :: K, I1, I2, J1, J2, IL, JL, KL, NBLAYER
      REAL,    INTENT (IN) :: RBLAYER, D1I(I2), S0(IL,JL,KL,3)
      REAL, INTENT (INOUT) :: X(IL,JL,KL,3)

C     Local constants:

      REAL,      PARAMETER :: ONE = 1., SUPPRESS = -999.
      LOGICAL,   PARAMETER :: CLOSED = .FALSE.
      CHARACTER, PARAMETER :: METHOD * 1 = 'B'

C     Local variables:

      INTEGER    I, J, JEVAL, NJ
      REAL       DENOM, WEIGHT
      LOGICAL    NEW
      REAL, DIMENSION (J1:J2) :: TI1, TI2, TNEW, TJ, XJ, YJ, ZJ

C     Execution:

C     Unnormalized I1 and I2 edge distributions for blending:

      TI1 = S0(I1,J1:J2,K,2)
      TI2 = S0(I2,J1:J2,K,2)

      NJ = J2 - J1 + 1
      DENOM = ONE / (S0(I2,J1,K,1) - S0(I1,J1,K,1))

C     For each interior radial line ...

      DO I = I1 + 1, I2 - 1

         TJ = S0(I,J1:J2,K,2) ! Unnormalized radial arcs

C        Redistribute the arc lengths via a combination of those at I1 & I2:

         WEIGHT = (S0(I,J1,K,1) - S0(I1,J1,K,1)) * DENOM

         CALL RADIAL_BLEND (NJ, TI1, TI2, TJ, WEIGHT, D1I(I),
     >                      NBLAYER, RBLAYER, TNEW)

C        Redistribute X, Y, Z:

         XJ = X(I,J1:J2,K,1)
         YJ = X(I,J1:J2,K,2)
         ZJ = X(I,J1:J2,K,3)

         NEW = .TRUE.
         JEVAL = J1

         DO J = J1 + 1, J2 - 1

            CALL PLSCRV3D (NJ, XJ, YJ, ZJ, TJ, METHOD, NEW, CLOSED,
     >                     TNEW(J), JEVAL, X(I,J,K,1), X(I,J,K,2),
     >                     X(I,J,K,3), SUPPRESS) ! Derivatives
            NEW = .FALSE.

         END DO

      END DO ! Next I

      END SUBROUTINE REGRID_KPLANE

C*******************************************************************************
C
      SUBROUTINE REGRID_LINE (IL, JL, KL, I, K, XYZ, D1, D2,
     >                        NBLAYER, RBLAYER, T, X, Y, Z, LUNERR)
C
C     Redistribution of a radial edge line as needed for correct grid warping
C     in indirect mode, where the initial distribution must use Euler spacing.
C
C     08/19/98  DAS  Adaptation of OFFLEDGE.
C
C*******************************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) :: IL, JL, KL, I, K, NBLAYER, LUNERR
      REAL, INTENT (INOUT) :: XYZ(IL,JL,KL,3)
      REAL, INTENT (IN)    :: D1, D2, RBLAYER
      REAL, INTENT (OUT)   :: T(JL), X(JL), Y(JL), Z(JL)

C     Local constants:

      REAL,      PARAMETER :: SUPPRESS = -999., ZERO = 0.
      LOGICAL,   PARAMETER :: CLOSED = .FALSE.
      CHARACTER, PARAMETER :: METHOD * 1 = 'B'

C     Local variables:

      INTEGER    IER, J, JEVAL
      REAL       TOTAL, TNEW(JL)
      LOGICAL    NEW

C     Execution:

      X = XYZ(I,1:JL,K,1) ! Get the initial Euler distribution ...
      Y = XYZ(I,1:JL,K,2)
      Z = XYZ(I,1:JL,K,3)

      CALL CHORDS3D (JL, X, Y, Z, .FALSE., TOTAL, T)

      TNEW(1)  = ZERO ! ... and redistribute it
      TNEW(JL) = TOTAL

      CALL BLGRID (JL, D1, D2, NBLAYER, RBLAYER, TNEW, LUNERR, IER)

      NEW = .TRUE.
      JEVAL = 1

      DO J = 2, JL - 1

         CALL PLSCRV3D (JL, X, Y, Z, T, METHOD, NEW, CLOSED, TNEW(J),
     >                  JEVAL, XYZ(I,J,K,1), XYZ(I,J,K,2), XYZ(I,J,K,3),
     >                  SUPPRESS)
         NEW = .FALSE.

      END DO

      END SUBROUTINE REGRID_LINE

C*******************************************************************************
C
      SUBROUTINE REGRID_VOLS (NCALL)
C
C     Redistribution of all radial lines in the subblock interiors via
C     2-D TFI applied to the radial distributions off the J = 1 face edges.
C
C     08/08/98  DAS  Introduced (flagged) cranks as boundaries on the wing.
C
C*******************************************************************************

C     Global variables:

      USE LIM
      USE LUNS

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) :: NCALL

C     Local variables:

      INTEGER   I1, I2, K1, K2, M, N
      REAL(KIND=4) CPU1, CPU2

C     Execution:

      IF (NCALL < 0) CALL SECOND (CPU1)

C     Unnormalized radial arcs for new faces and interim interior
C     radial lines of all subblock volumes (entire grid):

      CALL CHORDS_VOLUME ()

C     Process all subblocks.  At each (i,k) of the j = 1 face, apply 2-D TFI
C     to the radial distributions of subblock faces (for j > NBLAYER, that is).

      DO M = 2, NKREGRID

         K1 = KREGRID(M-1)
         K2 = KREGRID(M)

         DO N = 1, 4 ! Quadrants

            I1 = IFACE(N)
            I2 = IFACE(N+1)

            CALL RADIAL_TFI2D (N, I1, I2, K1, K2)

         END DO

      END DO

      IF (NCALL < 0) THEN
         CALL SECOND (CPU2)
         CALL CPUTIME (CPU1, CPU2, 'REGRID_VOLS',
     >                 'to redistribute subblock volumes', IWRIT)
         CPU1 = CPU2
      END IF

      END SUBROUTINE REGRID_VOLS

C***********************************************************************
C
      SUBROUTINE VSHEET ()
C
C     VSHEET generates the edges of the vortex sheet grid. Filling the
C     interior of the wake sheet is now done in GRID_FACES.
C     Z0PLANE does the body boundary.  This version forces inclusion of
C     the last fuselage station.
C
C***********************************************************************

C     Global variables:

      USE DIMENS
      USE FUSEL
      USE INDEX
      USE LIM
      USE LOGIC
      USE LUNS
      USE OPTIONS
      USE OUTR
      USE SPACING
      USE WNGSRF

      IMPLICIT NONE

C     Local constants:

      REAL,      PARAMETER :: HALF = 0.5, ONE = 1., ZERO = 0.
      LOGICAL,   PARAMETER :: NEW = .TRUE.
      CHARACTER, PARAMETER :: LOOSE * 1 = 'B'

C     Local variables:

      INTEGER
     >   I, IEDGE, IER, II, J, JCUT, K, NUMI
      REAL
     >   T(IL/2), TLINE(6), XLINE(6), YLINE(6), ZLINE(6), YTEU(KL),
     >   BETA, DX1, DX2, DY1, DY2, DZ2, DYPCT, DYPCT1, R, XEDGE, XTAIL,
     >   YBOT, YEDGE, YTOP
      LOGICAL
     >   LINEAR

C     Execution:

C     Downstream boundaries of the upper/lower vortex sheet:

      XW(1,1:KL)  = XLBACK(1,1:KL)
      YW(1,1:KL)  = YLBACK(1,1:KL)
      ZW(1,1:KL)  = ZLBACK(1,1:KL)
      XW(IL,1:KL) = XUBACK(1,1:KL)
      YW(IL,1:KL) = YUBACK(1,1:KL)
      ZW(IL,1:KL) = ZUBACK(1,1:KL)

C     Outer boundary at or beyond the tip (upper half of vortex sheet).
C     Smooth it as for other lines aft of the trailing edge.

      NUMI = IL - ITU + 1
      DX1  = XW(ITU-1,KL) - XW(ITU-2,KL)

      CALL EXPDIS4 (3, XW(ITU,KL), XW(1,KL), DX1, NUMI,
     >              XW(ITU,KL), BETA, -IWRIT)

      DO I = 1, 3
         J = 3 - I
         XLINE(I) = XW(ITL+J,KL)
         YLINE(I) = YW(ITL+J,KL)
      END DO

      DX2 = XLINE(3) - XLINE(1)
      DY2 = YLINE(3) - YLINE(1)
      XLINE(4) = XLINE(3) + DX2
      YLINE(4) = YLINE(3) + DY2 * TEKMAX
      XLINE(5) = XW(1,KL) - DX2
      YLINE(5) = YW(1,KL)
      XLINE(6) = XW(1,KL)
      YLINE(6) = YW(1,KL)

      CALL LCSFIT (6, XLINE, YLINE, NEW, LOOSE, NUMI,
     >             XW(ITU,KL), YW(ITU,KL), YW(ITU,KL))

      DO I = ITU + 1, IL
         ZW(I,1)  = ZERO ! Helps copying to lower edge; Z0PLANE redoes body
         ZW(I,KL) = ZW(1,KL)
      END DO

C     Similarly for the wake edge aft of the tip, but allow for TE thickness:

      DX1 = XW(ITU-1,KTIP) - XW(ITU-2,KTIP)

      CALL EXPDIS4 (3, XW(ITU,KTIP), XW(1,KTIP), DX1, NUMI,
     >              XW(ITU,KTIP), BETA, -IWRIT)

      DO I = 1, 3
         J = 3 - I
         XLINE(I) = (XW(ITL+J,KTIP) + XW(ITU-J,KTIP)) * HALF
         YLINE(I) = (YW(ITL+J,KTIP) + YW(ITU-J,KTIP)) * HALF
         ZLINE(I) =  ZW(ITL+J,KTIP)
      END DO
      DX2 = XLINE(3) - XLINE(1)
      DY2 = YLINE(3) - YLINE(1)
      DZ2 = ZLINE(3) - ZLINE(1)
      XLINE(4) = XLINE(3) + DX2
      YLINE(4) = YLINE(3) + DY2 * TETIP
      ZLINE(4) = ZLINE(3) + DZ2
      XLINE(5) = XW(1,KTIP) - DX2
      YLINE(5) = YW(1,KTIP)
      ZLINE(5) = ZW(1,KTIP)
      XLINE(6) = XW(1,KTIP)
      YLINE(6) = YW(1,KTIP)
      ZLINE(6) = ZW(1,KTIP)

C     Preserve any trailing edge thickness by leaving point ITU alone:

      CALL LCSFIT (6, XLINE, YLINE, NEW, LOOSE, NUMI - 1,
     >             XW(ITU+1,KTIP), YW(ITU+1,KTIP), YW(ITU+1,KTIP))

      CALL LCSFIT (6, XLINE, ZLINE, NEW, LOOSE, NUMI,
     >             XW(ITU,KTIP), ZW(ITU,KTIP), ZW(ITU,KTIP))


C     Aft of the root trailing edge, in the plane of the wing, start with a
C     one-sided stretching of X to the downstream boundary. If a fuselage
C     is present, the (interior) increments at the end-of-body will be
C     imposed later, but for TE-to-tail and tail-to-back-plane, the initial
C     one-sided distribution gives one of the two Vinokur increments needed.
C
C     Upper half of inboard vortex edge, TE-to-back-plane:

      DX1 = XW(ITU-1,1) - XW(ITU-2,1)

      CALL EXPDIS4 (3, XW(ITU,1), XW(1,1), DX1, NUMI,
     >              XW(ITU,1), BETA, -IWRIT)

      IF (WINGONLY) THEN

         DO I = 1, 3
            J = 3 - I
            XLINE(I) = (XW(ITL+J,1) + XW(ITU-J,1)) * HALF
            YLINE(I) = (YW(ITL+J,1) + YW(ITU-J,1)) * HALF
         END DO
         DX2 = XLINE(3) - XLINE(1)
         DY2 = YLINE(3) - YLINE(1)
         XLINE(4) = XLINE(3) + DX2
         YLINE(4) = YLINE(3) + DY2 * TEROOT
         XLINE(5) = XW(1,1)  - DX2
         YLINE(5) = YW(1,1)
         XLINE(6) = XW(1,1)
         YLINE(6) = YW(1,1)

         CALL LCSFIT (6, XLINE, YLINE, NEW, LOOSE, NUMI - 1,
     >                XW(ITU+1,1), YW(ITU+1,1), YW(ITU+1,1))

      ELSE

C        Force inclusion of the last fuselage station:

         NUMI  = IBU - ITU + 1
         YEDGE = YW(ITU,1)
         XEDGE = XW(ITU,1)
         DX1   = XW(ITU+1,1) - XEDGE
         DX2   = D1TAIL
         XTAIL = XF(NFSTN)

         CALL HTDIS4 (.TRUE., XEDGE, XTAIL, DX1, DX2, NUMI, XW(ITU,1),
     >                -IWRIT, IER)

         IF (IER /= 0 .AND. IER /= 3) THEN
            WRITE (IWRIT,*) 'VSHEET: IER, N, XEDGE, XTAIL, D1, D2 = ',
     >                      IER, NUMI, XEDGE, XTAIL, DX1, DX2
            STOP
         END IF

         NUMI = IL - IBU + 2
         DX1  = XTAIL - XW(IBU-1,1)
         DX2  = XW(1,1) - XW(IL-1,1)

         CALL HTDIS4 (.TRUE., XW(IBU-1,1), XW(1,1), DX1, DX2, NUMI,
     >                XW(IBU-1,1), -IWRIT, IER)

         IF (IER /= 0 .AND. IER /= 3) THEN
            WRITE (IWRIT,*) 'VSHEET: IER, N, XTAIL, XEND, D1, D2 = ',
     >                      IER, NUMI, XTAIL, XW(1,1), DX1, DX2
            STOP
         END IF

         XW(IBU,1) = XTAIL ! Restore it precisely


C        Interpolation of Y between trailing edge and tail.
C        The parametric body grid is assisted if this boundary passes through
C        a geometry point on the last station.

         YLINE(6) = (YF(1,NFSTN) + YF(NJBODY,NFSTN)) * HALF
         JCUT = NJBODY / 2

         CALL INTERVAL (NJBODY, YF(1,NFSTN), YLINE(6), ONE, JCUT)
         YLINE(6) = YF(JCUT,NFSTN)

         LINEAR = .FALSE. ! *** Revert to linear if odd geometry forces it

         IF (LINEAR) THEN
            IEDGE = NFSTN / 2

            CALL INTERVAL (NFSTN, XF, XEDGE, ONE, IEDGE)

            R     = (XEDGE - XF(IEDGE)) / (XF(IEDGE+1) - XF(IEDGE))
            YBOT  = (ONE - R)*YF(1,IEDGE)      + R*YF(1,IEDGE+1)
            YTOP  = (ONE - R)*YF(NJBODY,IEDGE) + R*YF(NJBODY,IEDGE+1)
            DYPCT = (YEDGE - YBOT) / (YTOP - YBOT)

            DO I = ITU + 1, IBU - 1
               DO II = IEDGE, NFSTN - 1
                  IF (XW(I,1) >= XF(II) .AND. XW(I,1) <= XF(II+1))
     >               EXIT
               END DO
               R    = (XW(I,1) - XF(II)) / (XF(II+1) - XF(II))
               YBOT = (ONE - R) * YF(1,II)      + R * YF(1,II+1)
               YTOP = (ONE - R) * YF(NJBODY,II) + R * YF(NJBODY,II+1)
               R    = (XW(I,1) - XEDGE) / (XTAIL - XEDGE)
               DYPCT1  = DYPCT + (HALF - DYPCT) * R
               YW(I,1) = YBOT + DYPCT1 * (YTOP - YBOT)
            END DO
            YW(IBU,1) = YLINE(6)

         ELSE ! Bisect the trailing edge angle with a spline:

            DO I = 1, 3
               J = 3 - I
               XLINE(I) = (XW(ITL+J,1) + XW(ITU-J,1)) * HALF
               YLINE(I) = (YW(ITL+J,1) + YW(ITU-J,1)) * HALF
            END DO
            DX1 = XLINE(3) - XLINE(2)
            DY1 = YLINE(3) - YLINE(2)
            XLINE(4) = XLINE(3) + DX1
            YLINE(4) = YLINE(3) + DY1 * TEROOT
            XLINE(5) = XTAIL - DX1
            XLINE(6) = XTAIL
            YLINE(5) = YLINE(6)

            CALL LCSFIT (6, XLINE, YLINE, NEW, LOOSE, IBU - ITU,
     >                   XW(ITU+1,1), YW(ITU+1,1), YW(ITU+1,1))
         END IF

C        Aft of the tail:

         R = (YW(1,1) - YW(IBU,1)) / (XW(1,1) - XTAIL)
         DO I = IBU + 1, IL
            YW(I,1) = YW(IBU,1) + (XW(I,1) - XTAIL) * R
         END DO
      END IF


C     Lower edges of the vortex sheet (but the body Y is redone in Z0PLANE):

      DO I = ITU + 1, IL
         II = IL + 1 - I
         XW(II,1)    = XW(I,1)
         YW(II,1)    = YW(I,1)
         ZW(II,1)    = ZW(I,1)
         XW(II,KTIP) = XW(I,KTIP) ! See above
         YW(II,KTIP) = YW(I,KTIP)
         ZW(II,KTIP) = ZW(I,KTIP)
         XW(II,KL)   = XW(I,KL)
         YW(II,KL)   = YW(I,KL)
         ZW(II,KL)   = ZW(I,KL)
      END DO

      END SUBROUTINE VSHEET

C*******************************************************************************
C
      SUBROUTINE WARP_VOLS (NCALL)
C
C     Perturb all subblock interiors given the subblock faces.
C
C*******************************************************************************

C     Global variables:

      USE FACES
      USE LIM
      USE LUNS
      USE OPTIONS
      USE XYZ
      USE XYZ0

      IMPLICIT NONE

C     Argument:

      INTEGER, INTENT (IN) :: NCALL

C     Local variables:

      INTEGER I1, I2, K1, K2, M, N
      REAL(KIND=4) CPU1, CPU2

C     Execution:

      IF (NCALL < 0) CALL SECOND (CPU1)

      IF (NNTIP == 0) THEN ! 8-block grid perturbation

         DO M = 1, 2 ! Blocks in the K direction

            K1 = KFACE(M)
            K2 = KFACE(M+1)

            DO N = 1, 4 ! Quadrants

               I1 = IFACE(N)
               I2 = IFACE(N+1)

               CALL WARP3D (1, IL, 1, JL, 1, KL, I1, I2, 1, JL, K1, K2,
     >                      X0(1,1,1,1), X0(1,1,1,2), X0(1,1,1,3), S0,
     >                      DFACEI,     DFACEJ,     DFACEK,
     >                      X(1,1,1,1), X(1,1,1,2), X(1,1,1,3))
            END DO

         END DO

      ELSE ! Rounded tip; perturb as 4 blocks

         DO M = 1, 2 ! Blocks in the K direction

            K1 = KFACE(M)
            K2 = KFACE(M+1)

            DO N = 1, 3, 2 ! I.e., 1:ILE and ILE:IL for each K range

               I1 = IFACE(N)
               I2 = IFACE(N+2)

               CALL WARP3D (1, IL, 1, JL, 1, KL, I1, I2, 1, JL, K1, K2,
     >                      X0(1,1,1,1), X0(1,1,1,2), X0(1,1,1,3), S0,
     >                      DFACEI,     DFACEJ,     DFACEK,
     >                      X(1,1,1,1), X(1,1,1,2), X(1,1,1,3))
            END DO

         END DO

      END IF

      IF (NCALL < 0) THEN
         CALL SECOND (CPU2)
         CALL CPUTIME (CPU1, CPU2, 'WARP_VOLS',
     >                 'to perturb the volume grid', IWRIT)
      END IF

      END SUBROUTINE WARP_VOLS

C***********************************************************************
C
      SUBROUTINE WFINTR (NCALL, KWING1, KWING2, MXKWING, IL, ILE, ITL,
     >                   ITU, XRG, YRG, ZRG, NJF, NIF, MXJBODY, MXIBODY,
     >                   XF, YF, ZF, IDEGBOD, XINTR, YINTR, ZINTR,
     >                   TINTR, UINTR, VINTR, KINTR, JINTR, IINTR,
     >                   IWRIT, FAIL)
C
C     WFINTR calculates a wing/fuselage intersection given a regular
C     geometry mesh for each.  This version assumes the entire wing is
C     regularized already.  It returns the "t" values associated with
C     each spanwise line at the intersection, as needed for generating
C     an arc-based spanwise wing surface grid.  Intersection (x,y,z)s
C     can be evaluated by the caller with the rest of the surface mesh.
C
C     09/26/97  DAS  Adaptation of the original SYN87 routine for use
C                    AFTER regularizing the wing sections, not before.
C     04/27/98   "   Retrofitted WBINTR along AEROSURF lines.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER NCALL               ! I   -3 means first call for this case
      INTEGER KWING1, KWING2      ! I   Subset of wing sections to treat
      INTEGER MXKWING             ! I   K and I dimensions of wing section
      INTEGER IL                  ! I   arrays
      INTEGER ILE, ITL, ITU       ! I   Leading & trailing edge wing indices
      REAL    XRG(IL,MXKWING),    ! I   Regularized wing sections
     >        YRG(IL,MXKWING),
     >        ZRG(IL,MXKWING)
      INTEGER NJF                 ! I   Active circumferential and axial
      INTEGER NIF                 !     indices of regular fuselage sections
      INTEGER MXJBODY             ! I   J and I dimensions of fuselage
      INTEGER MXIBODY             !     section arrays
      REAL    XF(MXIBODY),        ! I   Regularized fuselage sections
     >        YF(MXJBODY,MXIBODY),
     >        ZF(MXJBODY,MXIBODY)
      INTEGER IDEGBOD             ! I   3 means bicubics (INTSEC4) else
                                  !     1 means bilinear (INTSEC5)
      REAL    XINTR(IL),          !   O Intersection coordinates for
     >        YINTR(IL),          !     I = ITL : ITU
     >        ZINTR(IL)
      REAL    TINTR(IL)           !   O ITL:ITU elements are returned
                                  !     with "t" values at the intersection
      REAL    UINTR(IL),          !   O Corresponding body surface "u" and
     >        VINTR(IL)           !     "v" values (circumferential & axial)
      INTEGER KINTR(IL),          !   O Corresponding wing k & body cell (j,i)
     >        JINTR(IL),          !     indices
     >        IINTR(IL)
      INTEGER IWRIT               ! I   For diagnostics
      LOGICAL FAIL                !   O Success flag

C     Local constants:

      REAL, PARAMETER :: ONE = 1.

C     Local variables:

      INTEGER I, I1, I2, IFAIL, ISTN1, ISTN2, J
      LOGICAL FIRST
      REAL    YFMN, YFMX
      REAL, DIMENSION (MXJBODY,MXIBODY) :: UBODY, VBODY, XBODY

C     Storage:

      SAVE    ISTN1, ISTN2

C     Execution:

      FAIL  = .FALSE. ! Allow retrying in WBINTR
      FIRST = NCALL == -3

      IF (FIRST) THEN

C        Determine a subset of body sections containing the intersection:

         I1 = NIF / 5
         CALL INTERVAL (NIF, XF, XRG(ILE,KWING1), ONE, I1)

         I2 = I1
         CALL INTERVAL (NIF, XF, XRG(ILE,KWING1+1), ONE, I2)

         IINTR(ILE) = (I1 + I2) / 2
         ISTN1 = MAX (1, MIN (I1, I2) - 2)

         I1 = 4 * NIF / 5
         CALL INTERVAL (NIF, XF, XRG(ITL,KWING1), ONE, I1)

         I2 = I1
         CALL INTERVAL (NIF, XF, XRG(ITL,KWING1+1), ONE, I2)

         ISTN2 = MIN (NIF, MAX (I1, I2) + 3) ! INTERVAL always points to left

C        Set up for starting from the leading edge:

         I1 = ILE
         I2 = IINTR(I1)

         YFMN = 1.E+10
         YFMX = -YFMN

         DO J = 1, NJF
            YFMN = MIN (YF(J,I2), YFMN)
            YFMX = MAX (YF(J,I2), YFMX)
         END DO

         UINTR(I1) = (YRG(I1,KWING1) - YFMN) / (YFMX - YFMN)
         UINTR(I1) = MAX (0.3, MIN (UINTR(I1), 0.7)) ! Since Y is poor for LE u
         VINTR(I1) = (XF(I2) - XF(ISTN1)) / (XF(ISTN2) - XF(ISTN1))
         TINTR(I1) = MIN ((ZRG(I1,KWING1+1) - ZRG(I1,KWING1)) /
     >                    (ZRG(I1,KWING2)   - ZRG(I1,KWING1)), 0.9)
         JINTR(I1) = NJF * UINTR(I1)
         KINTR(I1) = (KWING1 + KWING2) / 2

      END IF

      DO I = ISTN1, ISTN2
         XBODY(1:NJF,I) = XF(I)
      END DO

      UBODY(1,ISTN1) = -ONE ! Tells INTSEC4/5 to parameterize the body

      I1 = ILE ! Upper surface first
      I2 = ITU

      DO ! Two-pass loop

         CALL WBINTR (IL, MXKWING, I1, I2, KWING1, KWING2,
     >                XRG, YRG, ZRG, MXJBODY, MXIBODY, 1, NJF,
     >                ISTN1, ISTN2, XBODY, YF, ZF, 1, IDEGBOD,
     >                UBODY, VBODY,
     >                XINTR, YINTR, ZINTR, UINTR, VINTR, TINTR,
     >                JINTR, IINTR, KINTR, FIRST, IWRIT, IFAIL, FAIL)

         IF (FAIL) THEN
            WRITE (IWRIT,'(/,A,I5)')
     >         ' WFINTR: Bad return from WBINTR. IFAIL:', IFAIL
            GO TO 999
         END IF

         IF (I1 < I2) THEN ! Do the lower surface from the leading edge
            I2 = ITL
         ELSE
            EXIT
         END IF

      END DO ! Two-pass loop

  999 RETURN

      END SUBROUTINE WFINTR

C***********************************************************************
C
      SUBROUTINE WINSERT (K1, K2, MXKWING, IL, ILE, ITL, ITU,
     >                    XRG, YRG, ZRG, NCRANK, ZCRANK, KCRANK, IWRIT)
C
C     WINSERT checks regularized sections of [one panel of] a wing-type
C     surface against the input crank stations, inserting sections if
C     necessary to allow arc-based spanwise gridding to capture the
C     specified ZCRANK(*) stations exactly.  Forcing a mesh station at
C     some inboard Z can avoid unwanted spanwise shifts in nacelle
C     locations caused by fuselage diameter changes during optimization.
C
C     10/02/97  DAS  Replaces automated detection of leading edge cranks,
C                    which are now inputs; absorbs the original "ZPLANAR"
C                    option as part of gridding w.r.t. arc length, not Z.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER K1, K2            ! I/O Subset of wing sections to treat;
                                !     K2 may be incremented on output
      INTEGER MXKWING           ! I   Max. # input wing sections
      INTEGER IL                ! I   Declared dimen. of regularized arrays
      INTEGER ILE, ITL, ITU     ! I   Leading & trailing edge mesh indices
      REAL    XRG(IL,MXKWING),  ! I/O Input regularized sections, possibly
     >        YRG(IL,MXKWING),  !     returned with new section(s) inserted
     >        ZRG(IL,MXKWING)   !     at some ZCRANK locations
      INTEGER NCRANK            ! I   No. of Z stations to capture in grid
      REAL    ZCRANK(NCRANK)    ! I/O Span stations to be captured by the
                                !     eventual grid: replaced by ZRG(ILE,*)
                                !     if within 0.001 * semispan, else a
                                !     new section is inserted & K2 increased
      INTEGER KCRANK(NCRANK)    !   O Indices of cranks in output sections
      INTEGER IWRIT             ! I   Unit # for printing crank locations

C     Local constants:

      REAL, PARAMETER :: ONE = 1., ZERO = 0.

C     Local variables:

      INTEGER I, K, KB, KC, KI, NMATCH
      REAL    R, TOLER, ZC

C     Execution:

C     Distinguish intended matches from insertions:

      TOLER = (ZRG(ILE,K2) - ZRG(ILE,K1)) * 0.001
      KB = K1
      NMATCH = 0

      DO KC = 1, NCRANK
         DO K = KB, K2
            IF (ABS (ZCRANK(KC) - ZRG(ILE,K)) < TOLER) THEN
               ZCRANK(KC) = ZRG(ILE,K)
               KCRANK(KC) = K
               KB = K + 1
               NMATCH = NMATCH + 1
               EXIT
            END IF
         END DO
      END DO

C     Insert section(s)?

      IF (NMATCH < NCRANK) THEN

         DO KC = 1, NCRANK
            K = K1 ! May be only 2 sections
            ZC = ZCRANK(KC)

   20       IF (ZC >= ZRG(ILE,K+1)) THEN
               K = K + 1
               IF (K < K2) GO TO 20
            END IF

            R = (ZC - ZRG(ILE,K)) / (ZRG(ILE,K+1) - ZRG(ILE,K))

            IF (R == ZERO) THEN

               KCRANK(KC) = K

            ELSE ! Insert a section

               KI = K + 1
               KCRANK(KC) = KI

               IF (K2 == MXKWING) THEN
                  WRITE (IWRIT, '(/,A,F11.4)')
     >               ' WINSERT: No room to insert a section at Z =', ZC
                  GO TO 900
               END IF

               DO K = K2, KI, -1 ! Make room
                  DO I = ITL, ITU
                     XRG(I,K+1) = XRG(I,K)
                     YRG(I,K+1) = YRG(I,K)
                     ZRG(I,K+1) = ZRG(I,K)
                  END DO
               END DO

               K2 = K2 + 1
               K  = KI
               DO I = ITL, ITU
                  XRG(I,K) = (ONE - R) * XRG(I,K-1) + R * XRG(I,K)
                  YRG(I,K) = (ONE - R) * YRG(I,K-1) + R * YRG(I,K)
                  ZRG(I,K) = (ONE - R) * ZRG(I,K-1) + R * ZRG(I,K)
               END DO
               ZRG(ILE,K) = ZC

            END IF
         END DO
      END IF

      GO TO 999

  900 STOP

  999 RETURN

      END SUBROUTINE WINSERT

C***********************************************************************
C
      SUBROUTINE WREGUL (K1, K2, MXIWING, MXKWING, XWG, YWG, ZWG,
     >                   ROUNDED, NLG, NW, IL, ILE, NHALF, SNORM,
     >                   XRG, YRG, ZRG)
C
C     WREGUL regularizes a range of wing geometry sections by imposing
C     the indicated chordwise mesh along the upper and lower surfaces.
C
C     06/14/96  DAS  Extracted from SYN87'S WSURF for reuse by NLCON.
C     10/02/96   "   Rounded leading edges are handled wrap-around now.
C     09/25/97   "   Fully argument-driven, Fortran 90 version with
C                    application to multiple panels in mind.
C     06/27/98   "   Abandoned the monotonic fit for X vs. T in view of
C                    the non-discrete method in RECTIFY now.
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER K1, K2                ! I   Subset of wing sections to treat
      INTEGER MXIWING, MXKWING      ! I   Declared dimens. of wing sections
      REAL    XWG(MXIWING,MXKWING), ! I   Defining section coordinates
     >        YWG(MXIWING,MXKWING),
     >        ZWG(MXIWING,MXKWING)
      REAL    ROUNDED(MXKWING)      ! I   0. means sharp leading edge
      INTEGER NLG(MXKWING),         ! I   #s of geometry pts. on lower surfaces
     >        NW(MXKWING)           !     and on whole sections (wrap-around)
      INTEGER IL                    ! I   Declared dimen. of regularized arrays
      INTEGER ILE                   ! I   Leading edge index of regular mesh
      INTEGER NHALF                 ! I   # mesh pts. for lower & upper halves
      REAL    SNORM(NHALF)          ! I   Normalized chordwise distribution
      REAL    XRG(IL,MXKWING),      !   O Regularized sections
     >        YRG(IL,MXKWING),
     >        ZRG(IL,MXKWING)

C     Local constants:
C     Use "loose" fits except monotonic X vs. T is best at leading edge.

      REAL,      PARAMETER :: ONE = 1., ZERO = 0.
      LOGICAL,   PARAMETER :: NEW = .TRUE., NORM = .FALSE.
      CHARACTER, PARAMETER :: LOOSE * 1 = 'B'

C     Local variables:

      INTEGER    I, ITL, K, N, NL, NU, NWING
      REAL       SLOWER, SMAX, SUPPER,
     >           SEVAL(2*NHALF), SWG(MXIWING), DERIVS(2*NHALF)

C     Execution:

      NWING = 2*NHALF - 1
      ITL   = ILE - NHALF + 1

      DO K = K1, K2
         N  = NW(K)
         NL = NLG(K)
         NU = N - NL + 1

         IF (ROUNDED(K) == ZERO) THEN ! Not rounded - do surfaces separately

C           Unnormalized geometry arc lengths along the lower surface:

            CALL CHORDS3D (NL, XWG(1,K), YWG(1,K), ZWG(1,K), NORM,
     >                     SMAX, SWG)

C           Denormalized chordwise mesh distribution:

            DO I = 1, NHALF
               SEVAL(I) = SMAX * (ONE - SNORM(NHALF + 1 - I))
            END DO

C           Lower surface interpolation along the arc:

            CALL LCSFIT (NL, SWG, XWG(1,K), NEW, LOOSE, NHALF, SEVAL,
     >                   XRG(ITL,K), DERIVS)
            CALL LCSFIT (NL, SWG, YWG(1,K), NEW, LOOSE, NHALF, SEVAL,
     >                   YRG(ITL,K), DERIVS)
            CALL LCSFIT (NL, SWG, ZWG(1,K), NEW, LOOSE, NHALF, SEVAL,
     >                   ZRG(ITL,K), DERIVS)

C           Likewise for the upper surface:

            CALL CHORDS3D (NU, XWG(NL,K), YWG(NL,K), ZWG(NL,K), NORM,
     >                     SMAX, SWG)

            DO I = 1, NHALF
               SEVAL(I) = SMAX * SNORM(I)
            END DO

            CALL LCSFIT (NU, SWG, XWG(NL,K), NEW, LOOSE, NHALF, SEVAL,
     >                   XRG(ILE,K), DERIVS)
            CALL LCSFIT (NU, SWG, YWG(NL,K), NEW, LOOSE, NHALF, SEVAL,
     >                   YRG(ILE,K), DERIVS)
            CALL LCSFIT (NU, SWG, ZWG(NL,K), NEW, LOOSE, NHALF, SEVAL,
     >                   ZRG(ILE,K), DERIVS)

        ELSE ! Rounded leading edge case:

            CALL CHORDS3D (N, XWG(1,K), YWG(1,K), ZWG(1,K), NORM,
     >                     SMAX, SWG)

            SLOWER = SWG(NL)
            SUPPER = SMAX - SLOWER

            DO I = 1, NHALF
               SEVAL(I) = SLOWER * (ONE - SNORM(NHALF + 1 - I))
               SEVAL(I + NHALF - 1) = SLOWER + SUPPER * SNORM(I)
            END DO

            CALL LCSFIT (N, SWG, XWG(1,K), NEW, LOOSE, NWING, SEVAL,
     >                   XRG(ITL,K), DERIVS)
            CALL LCSFIT (N, SWG, YWG(1,K), NEW, LOOSE, NWING, SEVAL,
     >                   YRG(ITL,K), DERIVS)
            CALL LCSFIT (N, SWG, ZWG(1,K), NEW, LOOSE, NWING, SEVAL,
     >                   ZRG(ITL,K), DERIVS)
        END IF

      END DO

      END SUBROUTINE WREGUL

C***********************************************************************
C
      SUBROUTINE WSURF (NCALL, FAIL)
C
C     WSURF generates a structured mesh on a wing-type surface.
C     First, it regularizes the wing defining sections, then it intersects
C     them with the fuselage if present.  Any leading edge cranks are
C     identified, and an option to force a mesh line at the given ZPLANAR
C     (as needed for SYN87's pseudo-nacelle option) is handled by inserting
C     a section at ZPLANAR then including it among the crank stations.
C
C     This version bases the spanwise grid on arc length, not Z, in
C     anticipation of treating winglets.  Performing the intersection
C     calculation after regularizing rather than before is part of the
C     arc-based strategy.
C
C     09/09/94  DAS  Replaced TAINT degree 3 interpolation with LCSFIT.
C     10/04/94   "   Eliminated overwriting of XW/YW/ZW geometry arrays
C                    with surface mesh values (bad dimension-wise) by
C                    introducing XWG/YWG/ZWG from SOLVE, not XW/YW/ZW,
C                    and outputting XW/YW/ZW here.
C     03/30/95   "   Incorporated FOILGRID for chordwise distributions.
C     04/26/95   "   Moved some vortex sheet edge stuff to VSHEET and
C                    arranged for C-O mesh option beyond tip (NNTIP = 1).
C     05/29/95   "   Installed round-tip option (NNTIP = 2).
C     10/26/95   "   Planar mode mods: Overwrite first regularized
C                    section with a section at Z = ZPLANAR.
C     02/18/96   "   Automated handling of leading edge crank(s).
C     03/01/96   "   Sections inside the body are OK now.
C     11/25/96 DS/JR ZPLANAR option added to avoid spanwise shifts of
C                    delta Cp effects.
C     03/17/97  DAS  Vinokur spanwise distribution on inner panel can
C                    fail - went to uniform instead.
C     09/27/97   "   Switched from Z-based to arc-based spanwise grid,
C                    along the lines needed for SYN87-MB paneling;
C                    abandoned "planar" option (NNROOT = 0); ZPLANAR
C                    option now covered by input ZCRANKs.
C     11/23/98   "   Tip mean line is now redistributed to reflect SNORM.
C
C***********************************************************************

C     Global variables:

      USE CONSTS
      USE DIMENS
      USE DMAX
      USE FUSEL
      USE INDEX
      USE LIM
      USE LOGIC
      USE LUNS
      USE OPTIONS
      USE OUTR
      USE SIZES
      USE SPACING
      USE TFI
      USE WNGSRF

      IMPLICIT  NONE

C     Arguments:

      INTEGER, INTENT (IN)  :: NCALL
      LOGICAL, INTENT (OUT) :: FAIL

C     Local constants:

      INTEGER,   PARAMETER :: NPTS = 7 ! For each I of tip rounding (NNTIP=2)
      REAL,      PARAMETER :: HALF = 0.5, ONE = 1., ZERO = 0.,
     >                        SUPPRESS = -999.
      LOGICAL,   PARAMETER :: NEW = .TRUE., NORM = .FALSE.
      CHARACTER, PARAMETER :: LINEAR * 1 = 'L', LOOSE * 1 = 'B'

C     Local variables:

      INTEGER
     >   I, IER, II, J, K, K1, K2, KTIPM1, L, NREGSEC, NTIP
      REAL
     >   TCRANK(MXCRANK), TEVAL(KTIP), TNORM(KTIP), TRG(MXKWING),
     >   TTIPM(NHALF), XINTR(IL), YINTR(IL), ZINTR(IL),
     >   TITIP(NPTS), XITIP(NPTS), YITIP(NPTS), ZITIP(NPTS),
     >   XTIPM(IL), YTIPM(IL), XYZG(MXKWING), XYZM(KTIP),
     >   ANGLE, BETA, C, CHDCL, CTIP, DTKTIPM1, DT, DY, DZ, HALFT,
     >   R, RNK, T1, T2, TDIHED, TK, TSCALE, TSHIFT, TSWEEP, TTOTAL,
     >   XIK, XMKTIP, XMKTIPM1, YIK, YMKTIP, YMKTIPM1, ZIK, ZMKTIP,
     >   ZMEAN, ZNEWC, ZNEWS, ZNEWT, ZUNIF
      LOGICAL
     >   NEWDATA

C     Execution:

      IF (NCALL == -3) THEN

C        Set up the normalized spanwise grid just once:

         TNORM(1) = ZERO
         DO K = 2, KTIP - 1
            ZUNIF    = REAL (K - 1) / REAL (KTIP - 1)
            ANGLE    = PI * ZUNIF
            ZNEWS    = SIN (ANGLE * HALF)
            ZNEWC    = (ONE - COS (ANGLE)) * HALF
            ZNEWT    = (ONE - SPNSPC) * ZUNIF + SPNSPC * ZNEWC
            TNORM(K) = (ONE - SPNSPS) * ZNEWT + SPNSPS * ZNEWS
         END DO
         TNORM(KTIP) = ONE

         IF (WINGONLY) THEN
            KWING1 = 1 ! Saves resetting it when suppressing the body
            TINTR(ITL:ITU) = ZERO
         END IF

      END IF


C     Update the planform reference area (not used at present):

      SREF2 = ZERO
      DO K = 2, NWSEC
         SREF2 = SREF2 + (ZWG(1,K)-ZWG(1,K-1))*(CHORDG(K)+CHORDG(K-1))
      END DO
      SREF2 = SREF2 * HALF


C     Regularize the wing sections to suit the eventual chordwise grid:

      CALL WREGUL (1, NWSEC, MXIWING, MXKWING, XWG, YWG, ZWG, ROUNDED,
     >             NLG, NW, IL, ILE, NHALF, SNORM, XRG, YRG, ZRG)


C     Check input leading edge cranks, inserting sections if necessary:

      NREGSEC = NWSEC ! May be incremented upon return

      IF (NCRANK > 0) THEN

         CALL WINSERT (1, NREGSEC, MXKWING, IL, ILE, ITL, ITU,
     >                 XRG, YRG, ZRG, NCRANK, ZCRANK, KCRANK, IWRIT)
      END IF


C     Wing/fuselage intersection?

      IF (.NOT. WINGONLY) THEN ! Else TINTR (*) are zeroed above

         CALL WFINTR (NCALL, KWING1,KWING2, MXKWING, IL, ILE, ITL, ITU,
     >                XRG, YRG, ZRG, NJBODY, NFSTN, NJBODY, MXIBODY,
     >                XF, YF, ZF, IDEGBOD, XINTR, YINTR, ZINTR,
     >                TINTR, UINTR, VINTR, KINTR, JINTR, IINTR,
     >                IWRIT, FAIL)

         IF (FAIL) GO TO 999

      ELSE
         FAIL = .FALSE.
      END IF


C     Arc-length-based wing surface gridding considerata:

C     Blending the crank Ts into the nominal spanwise distribution can
C     lead to different K indices at the cranks for different I indices,
C     so use a single base distribution containing any cranks.
C     Shape changes could also lead to different Ks for different shapes,
C     so use the base distribution from the initial shape.
C     The leading edge distribution may be muddied by slat deflections,
C     so pick a representative one around quarter-chord.
C     For the ghost nacelle case, restarts with the nacelle set-up must
C     also use the base distribution to avoid unintended spanwise shifts
C     of the nacelle effects caused by different wing grids, so the base
C     distribution and crank information are saved as part of the set-up.

      IF (NCALL == -3 .AND. ABS (NNACEL) /= 2) THEN

         I = ILE - NHALF / 3

C        Unnormalized arc lengths between wing sections for this I:

         CALL CHORDSRF (IL, MXKWING, I, 1, NREGSEC, XRG, YRG, ZRG,
     >                  NORM, TTOTAL, TRG)

C        Denormalized nominal base distribution:

         TBASE(1) = TRG(KWING1) + TINTR(I) *
     >             (TRG(KWING2) - TRG(KWING1)) ! Intersection T
         TBASE(KTIP) = TTOTAL
         TTOTAL = TTOTAL - TBASE(1)

         DO K = 2, KTIP - 1
            TBASE(K) = TBASE(1) + TTOTAL * TNORM(K)
         END DO

         IF (NCRANK > 0) THEN

C           Blend crank T(s) into the base T distribution:

            DO K = 1, NCRANK
               TCRANK(K) = TRG(KCRANK(K))
            END DO

            CALL SMOOTHX (TBASE, KTIP, TCRANK, NCRANK, 3, TFIWORK, IER)

C           Match base grid Ks with geometry Ks at cranks:

            K1 = 1
            DO J = 1, NCRANK
               DO K = K1, KTIP
                  IF (TCRANK(J) == TBASE(K)) THEN
                     MCRANK(J) = K
                     K1 = K1 + 1
                     EXIT
                  END IF
               END DO
            END DO

            WRITE (IWRIT, '(/,A,6F11.4)')
     >         ' WSURF: Mesh stations to be forced at these Zs:',
     >         ZCRANK(1:NCRANK)
            WRITE (IWRIT, '(/,A,6I4)')  ' Corresponding mesh Ks:',
     >         MCRANK(1:NCRANK)

         END IF

C        Appending the tip as a final crank helps paneling between cranks
C        even if the input NCRANK is 0.

         KCRANK(NCRANK+1) = NREGSEC
         MCRANK(NCRANK+1) = KTIP

      END IF

      IF (NCALL == -3 .AND. ABS (NNACEL) == 2) THEN
         KCRANK(NCRANK+1) = NREGSEC ! These are missing from the .rso file
         MCRANK(NCRANK+1) = KTIP
      END IF

C     Arc-length-based spanwise gridding:

      DO I = ITL, ITU

C        Unnormalized spanwise arcs for this I:

         CALL CHORDSRF (IL, MXKWING, I, 1, NREGSEC, XRG, YRG, ZRG,
     >                  NORM, TTOTAL, TRG)

C        Transform the base distribution by panels:

         T1 = TRG(KWING1) + TINTR(I) * (TRG(KWING2) - TRG(KWING1))
         TEVAL(1) = T1
         K1 = 1

         DO J = 1, NCRANK + 1
            K2 = MCRANK(J)
            T2 = TRG(KCRANK(J))
            DT = TBASE(K2) - TBASE(K1)
            TSCALE = (T2 - T1) / DT
            TSHIFT = (TBASE(K2) * T1 - TBASE(K1) * T2) / DT

            DO K = K1 + 1, K2
               TEVAL(K) = TBASE(K) * TSCALE + TSHIFT
            END DO

            K1 = K2
            T1 = T2
         END DO
         TEVAL(KTIP) = TTOTAL

C        Piecewise linear interpolation for spanwise grid:

         XYZG(1:NREGSEC) = XRG(I,1:NREGSEC)

         CALL LCSFIT (NREGSEC, TRG, XYZG, NEW, LINEAR, KTIP, TEVAL,
     >                XYZM, TFIWORK)

         XW(I,1:KTIP) = XYZM(1:KTIP)

         XYZG(1:NREGSEC) = YRG(I,1:NREGSEC)

         CALL LCSFIT (NREGSEC, TRG, XYZG, NEW, LINEAR, KTIP, TEVAL,
     >                XYZM, TFIWORK)

         YW(I,1:KTIP) = XYZM(1:KTIP)

         XYZG(1:NREGSEC) = ZRG(I,1:NREGSEC)

         CALL LCSFIT (NREGSEC, TRG, XYZG, NEW, LINEAR, KTIP, TEVAL,
     >                XYZM, TFIWORK)

         ZW(I,1:KTIP) = XYZM(1:KTIP)

      END DO


C     Make sure all sections have index ILE at the true leading edges.
C     They are redistributed anyway unless they are sharp, to blend better.

      CALL RECTIFY (1, KTIP, KL, IL, ILE, NHALF, SNORM, XW, YW, ZW)

C     Planform quantities:

      DO K = 1, KTIP
         XWLE(K) = XW(ILE,K)
         YWLE(K) = YW(ILE,K)
         ZWLE(K) = ZW(ILE,K) ! The design code needs these in NACSHIFT
         XWTE(K) = XW(ITL,K)
         YWTE(K) = YW(ITL,K)
         CHORDM(K) = XWTE(K) - XWLE(K)
      END DO

      CTIP = CHORDM(KTIP)
      CHORDM(KTIP:KL) = CTIP

C     Beyond the tip for all calls:

      DZ = ZWLE(KTIP) - ZWLE(KTIP-1)

      IF (SWEEP == 999.) THEN
         TSWEEP = ((XWLE(KTIP)   + XWTE(KTIP)) -
     >             (XWLE(KTIP-1) + XWTE(KTIP-1))) / (DZ + DZ)
         SWEEP  = RAD2DEG * ATAN (TSWEEP)
      ELSE
         TSWEEP = TAN (SWEEP*DEG2RAD)
      END IF

      IF (BSWEEP == 999.) BSWEEP = SWEEP ! Default; allow 0. too
      TANSWEEP = TAN (BSWEEP*DEG2RAD)

      IF (DIHED == 999.) THEN
         TDIHED = ((YWLE(KTIP)   + YWTE(KTIP)) -
     >             (YWLE(KTIP-1) + YWTE(KTIP-1))) / (DZ + DZ)
         DIHED  = RAD2DEG * ATAN (TDIHED)
      ELSE
         TDIHED = TAN (DIHED*DEG2RAD)
      END IF


      IF (NNTIP == 0) THEN ! C-H grid all the way to the boundary

         CALL  EXPDIS4 (3, ZWLE(KTIP), ZMAX, DZ, KL - KTIP + 1,
     >                  ZWLE(KTIP), BETA, -IWRIT)

C        Tip section mean-line:

         K = KTIP + 1 ! Temporary usage before redistributing as X/YTIPM(*)

         L = ILE
         DO I = ILE, ITU
            XW(I,K) = (XW(I,KTIP) + XW(L,KTIP)) * HALF
            YW(I,K) = (YW(I,KTIP) + YW(L,KTIP)) * HALF
            ZW(I,K) =  ZWLE(KTIP)
            L = L - 1
         END DO

C        Impose the wing chordwise arc length distribution on the mean line:

         NTIP = ITU - ILE + 1

         CALL CHORDS3D (NTIP, XW(ILE,K), YW(ILE,K), ZW(ILE,K), NORM, DT,
     >                  TTIPM)

         NEWDATA = .TRUE.
         J  = 1
         L  = ILE
         II = ILE

         DO I = ILE, ITU

            T1 = SNORM(J) * DT

            CALL PLSCRV3D (NTIP, XW(ILE,K), YW(ILE,K), ZW(ILE,K), TTIPM,
     >                     LOOSE, NEWDATA, .FALSE., T1, II,
     >                     XTIPM(I), YTIPM(I), ZNEWT, SUPPRESS)

            NEWDATA = .FALSE.
            XTIPM(L) = XTIPM(I)
            YTIPM(L) = YTIPM(I)
            L = L - 1
            J = J + 1

         END DO

C        Degenerate section at K = KL:

         K = KL
         DZ = ZWLE(K) - ZWLE(KTIP)
         XWLE(K) = XWLE(KTIP) + DZ * TSWEEP
         XWTE(K) = XWLE(K)    + CTIP
         YWLE(K) = YWLE(KTIP) + DZ * TDIHED
         YWTE(K) = YWLE(K) ! Untwist at KL

         L = ILE
         J = 1
         DO I = ILE, ITU
            XW(I,K) = XWLE(K) + SNORM(J) * CTIP
            XW(L,K) = XW(I,K)
            YW(I,K) = YWLE(K)
            YW(L,K) = YWLE(K)
            ZW(I,K) = ZWLE(K)
            ZW(L,K) = ZWLE(K)
            L = L - 1
            J = J + 1
         END DO

C        Linear interpolation between mean tip line and K = KL section:

         DO K = KTIP + 1, KL - 1
            R = (ZWLE(K) - ZWLE(KTIP)) / DZ
            XWLE(K) = (ONE - R) * XWLE(KTIP) + R * XWLE(KL)
            XWTE(K) = XWLE(K) + CTIP
            YWLE(K) = (ONE - R) * YWLE(KTIP) + R * YWLE(KL)
            YWTE(K) = (ONE - R) * YWTE(KTIP) + R * YWLE(KL)

            DO I = ITL, ITU
               XW(I,K) = (ONE - R) * XTIPM(I) + R * XW(I,KL)
               YW(I,K) = (ONE - R) * YTIPM(I) + R * YW(I,KL)
               ZW(I,K) = ZWLE(K)
            END DO

         END DO

      ELSE ! Grid will be C-O beyond the tip

         XWLE(KTIP:KL) = XWLE(KTIP)
         YWLE(KTIP:KL) = YWLE(KTIP)
         XWTE(KTIP:KL) = XWTE(KTIP)
         YWTE(KTIP:KL) = YWTE(KTIP)

C        NNTIP = 1 case:  Grid the squared-off tip.
C        Do it if NNTIP = 2 too to cover the lead-/trailing edges and meanline:

         RNK = ONE / REAL (KL - KTIP)

         DO K = KTIP + 1, KL
            R = REAL (K - KTIP) * RNK
            DO I = ITL, ITU
               ZW(I,K) = ZW(I,KTIP)
               II = ITL + ITU - I
               XW(I,K) = R * (XW(I,KTIP) + XW(II,KTIP)) * HALF +
     >                   (ONE - R) * XW(I,KTIP)
               YW(I,K) = R * (YW(I,KTIP) + YW(II,KTIP)) * HALF +
     >                   (ONE - R) * YW(I,KTIP)
            END DO
         END DO

         IF (NNTIP == 2) THEN ! Generate a rounded tip cap grid

C           Parametric spline techniques on 7-pt. curves should suffice:

            KTIPM1 = KTIP - 1
            DO I = ITL + 1, ILE - 1
               XITIP(1) = XW(I,KTIPM1)
               YITIP(1) = YW(I,KTIPM1)
               ZITIP(1) = ZW(I,KTIPM1)
               XITIP(3) = XW(I,KTIP)
               YITIP(3) = YW(I,KTIP)
               ZITIP(3) = ZW(I,KTIP)
               II = ITL + ITU - I
               XITIP(5) = XW(II,KTIP)
               YITIP(5) = YW(II,KTIP)
               ZITIP(5) = ZW(II,KTIP)
               XITIP(7) = XW(II,KTIPM1)
               YITIP(7) = YW(II,KTIPM1)
               ZITIP(7) = ZW(II,KTIPM1)
               R = 0.99 ! Control end slopes of tip-cap curve at I
               XITIP(2) = (ONE - R) * XITIP(1) + R * XITIP(3)
               YITIP(2) = (ONE - R) * YITIP(1) + R * YITIP(3)
               ZITIP(2) = (ONE - R) * ZITIP(1) + R * ZITIP(3)
               XITIP(6) = (ONE - R) * XITIP(7) + R * XITIP(5)
               YITIP(6) = (ONE - R) * YITIP(7) + R * YITIP(5)
               ZITIP(6) = (ONE - R) * ZITIP(7) + R * ZITIP(5)

C              Extrapolate the mean-dihedral line half the tip thickness at I:

               HALFT = SQRT ((XITIP(5) - XITIP(3)) ** 2 +
     >                       (YITIP(5) - YITIP(3)) ** 2) * HALF
               XMKTIPM1 = (XITIP(1) + XITIP(7)) * HALF
               YMKTIPM1 = (YITIP(1) + YITIP(7)) * HALF
               XMKTIP   = XW(I,KL) ! From NNTIP = 1 case above
               YMKTIP   = YW(I,KL)
               ZMKTIP   = ZW(I,KL)
               DTKTIPM1 = SQRT ((XMKTIP - XMKTIPM1) ** 2 +
     >                          (YMKTIP - YMKTIPM1) ** 2 +
     >                          (ZMKTIP - ZITIP(1)) ** 2)
               R = (DTKTIPM1 + HALFT) / DTKTIPM1
               XITIP(4) = (ONE - R) * XMKTIPM1 + R * XMKTIP
               YITIP(4) = (ONE - R) * YMKTIPM1 + R * YMKTIP
               ZITIP(4) = (ONE - R) * ZITIP(1) + R * ZMKTIP

C              Remove any dihedral to enable a monotonic fit for Y:

               ANGLE = ATAN2 (YMKTIP   - YMKTIPM1,
     >                        ZITIP(3) - ZITIP(1)) * RAD2DEG ! X is unaffected

               CALL ROTATE2D (NPTS, ZITIP, YITIP,-ANGLE, ZMKTIP, YMKTIP)

               CALL CHORDS3D (NPTS, XITIP, YITIP, ZITIP, NORM, TTOTAL,
     >                        TITIP)

C              Lower tip cap stations:

               DT = (TITIP(4) - TITIP(3)) * RNK
               DO K = KTIP + 1, KL
                  TK = TITIP(3) + DT * REAL (K - KTIP)
                  CALL LCSFIT (NPTS,TITIP,XITIP,NEW,'B',1,TK,XIK,XIK)
                  CALL LCSFIT (NPTS,TITIP,YITIP,NEW,'M',1,TK,YIK,YIK)
                  CALL LCSFIT (NPTS,TITIP,ZITIP,NEW,'B',1,TK,ZIK,ZIK)
                  CALL ROTATE2D (1, ZIK, YIK, ANGLE, ZMKTIP, YMKTIP)
                  XW(I,K) = XIK
                  YW(I,K) = YIK
                  ZW(I,K) = ZIK
               END DO

C              Upper tip cap:

               DT = (TITIP(5) - TITIP(4)) * RNK
               DO K = KTIP + 1, KL
                  TK = TITIP(5) - DT * REAL (K - KTIP)
                  CALL LCSFIT (NPTS,TITIP,XITIP,NEW,'B',1,TK,XIK,XIK)
                  CALL LCSFIT (NPTS,TITIP,YITIP,NEW,'M',1,TK,YIK,YIK)
                  CALL LCSFIT (NPTS,TITIP,ZITIP,NEW,'B',1,TK,ZIK,ZIK)
                  CALL ROTATE2D (1, ZIK, YIK, ANGLE, ZMKTIP, YMKTIP)
                  XW(II,K) = XIK
                  YW(II,K) = YIK
                  ZW(II,K) = ZIK
               END DO

            END DO ! Next I

         ELSE ! NNTIP = 3: Grid an input tip cap geometry (incomplete)

         END IF

      END IF

  999 RETURN

      END SUBROUTINE WSURF

C***********************************************************************
C
      SUBROUTINE Z0PLANE (FULLGRID, FAIL)
C
C     Z0PLANE fills the Z = 0 symmetry plane and the generates the
C     fuselage surface mesh if there is one.
C
C     NREFN is the (arbitrary) number of points used here by the
C     original nonparametric method to obtain initial grid points
C     on the body (which are then redistributed).
C
C     Mar 1995  DAS  This version has a parametric body grid option.
C     May 1995   "   Removed all scaling by tip chord.
C     Oct 1995   "   Revived and refined the bilinear body option.
C     Feb 1996   "   Apply WARPQ3D off the body if possible.
C     May 1996   "   WARPQ3D must come after aft-body redistribution.
C     08/07/96   "   Adjusted axial redistribution for open tr. edge.
C     08/10/96   "   Installed ELLIP2D.
C     Jun 1998   "   A straight-line edge boundary off the crown/keel
C                    at I = ITL & ITU doesn't help much - abandoned it.
C
C***********************************************************************

C     Global variables:

      USE CONSTS
      USE DIMENS
      USE ELLIPTIC
      USE FACES
      USE FUSEL
      USE INDEX
      USE LIM
      USE LOGIC
      USE LUNS
      USE OUTR
      USE SPACING
      USE TFI
      USE WALL
      USE WNGSRF
      USE XYZ
      USE XYZ0

      IMPLICIT NONE

C     Arguments:

      LOGICAL, INTENT (IN)  :: FULLGRID
      LOGICAL, INTENT (OUT) :: FAIL

C     Local constants:

      REAL,    PARAMETER :: HALF = 0.5, ONE = 1., TOLER = 1.E-7,
     >                      ZERO = 0., SUPPRESS = -999.
      LOGICAL, PARAMETER :: NEWDAT = .TRUE.

C     Local variables:

      INTEGER
     >   I, I1, I2, ICROWN, IER, IKEEL, IS, ISMAT, ISTN, J, JBLEND,
     >   JCUT, JHALF, JJ, JJ1, JS, JSMAT, K, L, LF, M, N, NBAD,
     >   NCROWN, NREFNM1
      REAL
     >   BETA, D1, D2, D2MAX, DRANGE, DUL, DUNIFORM, DVL, DVU, DUU,
     >   DX, DY, DY1, EPS, P, PWRI, Q, R, RATIO, REFNM1, RJCRM1, RJM1,
     >   RLC, RLW, RUC, RUW, SFTOT, SINT, UCUT, USHIFT, UTOTAL, V1, V2,
     >   VCROWNI, VCUT, VKEELI, VTOTAL, XG, XHALF, XJ, XRANGE, XWALLJ,
     >   YCR, YG, YKL, YMID, YNOSE, YR, YRANGE, YU, YWALLJ, ZCUT, ZG,
     >   SMATRIX(4,4,3), UVRANGE(4), VCROWN(MXIBODY), VKEEL(MXIBODY),
     >   XLINE(JL), YLINE(JL), XWAT(5), YWAT(5), XYZCUT(3), YP(NJBODY),
     >   ZP(NJBODY)
      REAL, DIMENSION (NJBODY,MXIBODY) ::
     >   UBODY, VBODY, XBODY
      REAL, ALLOCATABLE ::
     >   XTMP(:,:), YTMP(:,:), ZTMP(:,:), STMP(:), XWATER(:)
      LOGICAL
     >   NEW
      CHARACTER
     >   FGMODE*4, BGMODE*4

C     Execution:

      FAIL  = .FALSE.
      ZWALL = ZERO ! (IL x JL); used for scratch elsewhere, so always reset it

      IF (WINGONLY) GO TO 800 ! Avoid a massive indent

      NBAD   = 0
      RJCRM1 = ONE / REAL (JCR - 1)

C     Set up the symmetry plane ("wall") as two regions above and below
C     the body, suitable for TFI and elliptic smoothing.

      XWALL(1:IL,1)   = XW(1:IL,1)   ! Wing/body intersection
      XWALL(1:IL,JCR) = XCROWN(1:IL) ! ... and C parts of ...
      XWALL(1:IL,JL)  = XOUT(1:IL)   ! ... off-body boundaries
      YWALL(1:IL,1)   = YW(1:IL,1)
      YWALL(1:IL,JCR) = YCROWN(1:IL)
      YWALL(1:IL,JL)  = YOUT(1:IL)

      XWALL( 1,1:JL) = XLBACK(1:JL,1) ! Far boundaries, incl. aft "fan"
      XWALL(IL,1:JL) = XUBACK(1:JL,1)
      YWALL( 1,1:JL) = YLBACK(1:JL,1)
      YWALL(IL,1:JL) = YUBACK(1:JL,1)

C     Body surface grid.
C     ------------------

      IF (.NOT. PARABODY) THEN ! Original body grid by linear techniques

         ALLOCATE (XTMP(NREFN,IL), YTMP(NREFN,IL), ZTMP(NREFN,IL),
     >             STMP(NREFN), XWATER(NREFN))

C        Determine the X & Y coordinates of "spokes" between the wing/body
C        intersection line and the corresponding points on the crown line.
C        A dense regular mesh in the plane is generated initially, smoothed,
C        and used to interpolate the body surface linearly.  The radial lines
C        are then redistributed linearly.

         NREFNM1 = NREFN - 1
         REFNM1  = ONE / REAL (NREFNM1)

C        Set up an artificial boundary between the nose and the projection
C        of the root leading edge on the symmetry plane.

         XRANGE = XW(ILE,1) - XCROWN(ILE)
         DUNIFORM = XRANGE * REFNM1
         D1 = DUNIFORM * 1.8 ! D1 & D2 are relatively the same as for the
         D2 = DRADIAL(1,1) * REAL (JCR) / REAL (NREFN) ! final distribution

         CALL HTDIS4 (.TRUE., XCROWN(ILE), XW(ILE,1), D1, D2, NREFN,
     >                XWATER, -IWRIT, IER)

         IF (IER /= 0 .AND. IER /= 3) THEN
            WRITE (IWRIT,*) 'Z0PLANE: HTDIS4 trouble. IER: ', IER
            WRITE (IWRIT,*) 'D1, D2, XRANGE: ', D1, D2, XRANGE
            GO TO 990
         END IF

C        Handling a low-or high-mounted wing is the problem here.
C        We want to split the crown-keel line by ~ half at the nose, yet
C        meet the root leading edge.  The following blends between the two.

         IKEEL  = ILE / 2   ! Pointer estimates
         ICROWN = ILE + IKEEL
         JBLEND = NREFN / 5 ! Say (allowing for bunching towards l.e.)
         JHALF  = 1
         XHALF  = XCROWN(ILE) + XRANGE * HALF

         DO J = 1, NREFNM1
            XJ = XWATER(NREFN + 1 - J)
            XTMP(J,ILE) = XJ
            IF (XJ > XHALF) JHALF = J

            CALL INTERVAL (ILE, XCROWN(ILE), XJ, ONE, ICROWN)
            I = ICROWN + ILE - 1
            R = (XJ - XCROWN(I)) / (XCROWN(I+1) - XCROWN(I))
            YCR = (ONE - R) * YCROWN(I) + R * YCROWN(I+1)

            CALL INTERVAL (ILE, XCROWN(1), XJ, -ONE, IKEEL)
            I = IKEEL + 1
            R = (XJ - XCROWN(I)) / (XCROWN(I-1) - XCROWN(I))
            YKL = (ONE - R) * YCROWN(I) + R * YCROWN(I-1)

            IF (J == 1) RATIO = (YCR - YW(ILE,1)) / (YCR - YKL)

            YR = (ONE - RATIO) * YCR + RATIO * YKL
            IF (J <= JBLEND) THEN
               YTMP(J,ILE) = YR
            ELSE                    ! Blend to half way at the nose
               YU = (YCR + YKL) * HALF
               R  = REAL (J - JBLEND) / REAL (NREFN - JBLEND)
               YTMP(J,ILE) = (ONE - R) * YR + R * YU
            END IF
         END DO

C        Use a monotonic spline to overcome discontinuities in YTMP(*,ILE):

         XWAT(1) = XTMP(1,ILE)
         YWAT(1) = YTMP(1,ILE)
         XWAT(2) = XTMP(2,ILE)
         YWAT(2) = YTMP(2,ILE)
         XWAT(3) = XTMP(JHALF,ILE)
         YWAT(3) = YTMP(JHALF,ILE)
         XWAT(4) = XTMP(NREFN-1,ILE)
         YWAT(4) = YTMP(NREFN-1,ILE)
         XWAT(5) = XTMP(NREFN,ILE)
         YWAT(5) = YTMP(NREFN,ILE)

         CALL LCSFIT (5, XWAT, YWAT, NEWDAT, 'M', NREFN, XTMP(1,ILE),
     >                YTMP(1,ILE), YTMP(1,ILE))

C        Set up the other boundaries:

         DO I = 1, IL
            XTMP(1,I) = XW(I,1)
            YTMP(1,I) = YW(I,1)
            XTMP(NREFN,I) = XCROWN(I)
            YTMP(NREFN,I) = YCROWN(I)
            ZTMP(NREFN,I) = ZERO
         END DO

         DO I = 1, IL, IL - 1 ! Uniform at the far downstream boundary
            DX = (XCROWN(I) - XW(I,1)) * REFNM1
            DY = (YCROWN(I) - YW(I,1)) * REFNM1
            DO J = 1, NREFN
               XTMP(J,I) = XW(I,1) + DX * REAL (J-1)
               YTMP(J,I) = YW(I,1) + DY * REAL (J-1)
            END DO
         END DO

C        Force capturing of the last body station:

         DO I = IBL, IBU, IBU - IBL
            DY = (YCROWN(I) - YW(I,1)) * REFNM1
            DO J = 1, NREFN
               XTMP(J,I) = XF(NFSTN)
               YTMP(J,I) = YW(I,1) + DY * REAL (J-1)
            END DO
         END DO

C        Fill and smooth the interior in four pieces.  NOTE:
C        This use of TFI2D and ELLIP2D differs from all others in the
C        dimensions NREFN x IL rather than IL x JL.

         CALL TFI2D (NREFN, 1, NREFN,   1, IBL, XTMP, YTMP, DFACEI)
         CALL TFI2D (NREFN, 1, NREFN, IBL, ILE, XTMP, YTMP, DFACEI)
         CALL TFI2D (NREFN, 1, NREFN, ILE, IBU, XTMP, YTMP, DFACEI)
         CALL TFI2D (NREFN, 1, NREFN, IBU,  IL, XTMP, YTMP, DFACEI)

         L = 12

         IF (PRINT2D) WRITE (IWRIT,'(/,A)')
     >      ' ELLIP2D on nonparametric body grid:'

         CALL ELLIP2D (NREFN, IL, 1, NREFN, 1, IBL,
     >                 XTMP,        YTMP,        DFACEK,   PRINT2D,
     >                 BGMODE2D(L), FGMODE2D(L), SPMODE(L),
     >                 ITMAX2D(L),  ITFLOAT,     ITFREEZE, JLAYER,
     >                 CONV2D(L),   CONVMIN(L),  DMAX2D,   OMG2D(L),
     >                 POWERI(L),   POWERJ(L),   EXPI12D(L),
     >                 EXPI22D(L),  EXPJ12D(L),  EXPJ22D(L),
     >                 FGLIM2D,     URFG2D,      URFLOAT)

         L = 13

         CALL ELLIP2D (NREFN, IL, 1, NREFN, IBL, ILE,
     >                 XTMP,        YTMP,        DFACEK,   PRINT2D,
     >                 BGMODE2D(L), FGMODE2D(L), SPMODE(L),
     >                 ITMAX2D(L),  ITFLOAT,     ITFREEZE, JLAYER,
     >                 CONV2D(L),   CONVMIN(L),  DMAX2D,   OMG2D(L),
     >                 POWERI(L),   POWERJ(L),   EXPI12D(L),
     >                 EXPI22D(L),  EXPJ12D(L),  EXPJ22D(L),
     >                 FGLIM2D,     URFG2D,      URFLOAT)

         CALL EXCHANGE (BGMODE2D(L), BGMODE)
         CALL EXCHANGE (FGMODE2D(L), FGMODE)

         CALL ELLIP2D (NREFN, IL, 1, NREFN, ILE, IBU,
     >                 XTMP,        YTMP,        DFACEK,   PRINT2D,
     >                 BGMODE,      FGMODE,      SPMODE(L),
     >                 ITMAX2D(L),  ITFLOAT,     ITFREEZE, JLAYER,
     >                 CONV2D(L),   CONVMIN(L),  DMAX2D,   OMG2D(L),
     >                -POWERI(L),   POWERJ(L),   EXPI22D(L),
     >                 EXPI12D(L),  EXPJ12D(L),  EXPJ22D(L),
     >                 FGLIM2D,     URFG2D,      URFLOAT)

         L = 12

         CALL EXCHANGE (BGMODE2D(L), BGMODE)
         CALL EXCHANGE (FGMODE2D(L), FGMODE)

         CALL ELLIP2D (NREFN, IL, 1, NREFN, IBU, IL,
     >                 XTMP,        YTMP,        DFACEK,   PRINT2D,
     >                 BGMODE,      FGMODE,      SPMODE(L),
     >                 ITMAX2D(L),  ITFLOAT,     ITFREEZE, JLAYER,
     >                 CONV2D(L),   CONVMIN(L),  DMAX2D,   OMG2D(L),
     >                -POWERI(L),   POWERJ(L),   EXPI22D(L),
     >                 EXPI12D(L),  EXPJ12D(L),  EXPJ22D(L),
     >                 FGLIM2D,     URFG2D,      URFLOAT)

C        Some lines may be too tight for low- or high-mounted wings,
C        so revert to Laplace smoothing forward of the trailing edge.
C
C        *** This could conceivably cause objective discontinuities ***
C            if a design has RATIO near 1/3 or 2/3 and is varying it.

         L = 14

         IF (RATIO > 0.6667) THEN ! Low-mounted wing

            IF (PRINT2D) WRITE (IWRIT,'(/,A)')
     >         ' Laplace smoothing of body (X,Y) grid below wing.'

            CALL ELLIP2D (NREFN, IL, 1, NREFN, ITL, ILE,
     >                    XTMP,        YTMP,        DFACEK,   PRINT2D,
     >                    BGMODE2D(L), FGMODE2D(L), SPMODE(L),
     >                    ITMAX2D(L),  ITFLOAT,     ITFREEZE, JLAYER,
     >                    CONV2D(L),   CONVMIN(L),  DMAX2D,   OMG2D(L),
     >                    POWERI(L),   POWERJ(L),   EXPI12D(L),
     >                    EXPI22D(L),  EXPJ12D(L),  EXPJ22D(L),
     >                    FGLIM2D,     URFG2D,      URFLOAT)

         ELSE IF (RATIO < 0.3333) THEN ! High-mounted

            CALL EXCHANGE (BGMODE2D(L), BGMODE)
            CALL EXCHANGE (FGMODE2D(L), FGMODE)

            IF (PRINT2D) WRITE (IWRIT,'(/,A)')
     >         ' Laplace smoothing of body (X,Y) grid above wing:'

            CALL ELLIP2D (NREFN, IL, 1, NREFN, ILE, ITU,
     >                    XTMP,        YTMP,        DFACEK,   PRINT2D,
     >                    BGMODE,      FGMODE,      SPMODE(L),
     >                    ITMAX2D(L),  ITFLOAT,     ITFREEZE, JLAYER,
     >                    CONV2D(L),   CONVMIN(L),  DMAX2D,   OMG2D(L),
     >                   -POWERI(L),   POWERJ(L),   EXPI22D(L),
     >                    EXPI12D(L),  EXPJ12D(L),  EXPJ22D(L),
     >                    FGLIM2D,     URFG2D,      URFLOAT)
         END IF


C        Bilinear interpolation  of the body surface at the dense mesh points:
C        --------------------------------------------------------------------

         LF = NFSTN - 1
         DO I = 1, IL
            DO J = 1, NREFNM1
               ZTMP(J,I) = ZERO
               XWALLJ = XTMP(J,I)

               IF (XWALLJ <= XF(NFSTN)) THEN ! Bracket this X by body stns.

                  CALL INTERVAL (NFSTN, XF, XWALLJ, ONE, LF)

C                 Generate an interpolated fuselage station at this X.

                  DO M = 1, NJBODY
                     R = (XWALLJ - XF(LF)) / (XF(LF+1) - XF(LF))
                     YP(M) = (ONE - R) * YF(M,LF) + R * YF(M,LF+1)
                     ZP(M) = (ONE - R) * ZF(M,LF) + R * ZF(M,LF+1)
                  END DO

C                 Interpolate the body at YWALLJ on this interpolated stn.

                  YWALLJ = YTMP(J,I)
                  DO M = 2, NJBODY
                     IF (YWALLJ >= YP(M-1) .AND. YWALLJ <= YP(M)) THEN
                        DY1 = YP(M) - YP(M-1)
                        IF (ABS (DY1) > 0.001) THEN
                           ZTMP(J,I) = ZP(M-1) + (ZP(M)  - ZP(M-1)) *
     >                                 (YWALLJ - YP(M-1)) / DY1
                        ELSE IF (J /= 1) THEN
                           ZTMP(J,I) = ZTMP(J-1,I) * REAL (J-NREFN-2) /
     >                                               REAL (J-NREFN-3)
                        END IF ! Else 0. above
                        EXIT
                     END IF
                  END DO

               END IF

               IF (J == 1) THEN
                  IF (I >= IBL .AND. I <= IBU) THEN
                     IF (I >= ITL .AND. I <= ITU) THEN
                        ZTMP(1,I) = MAX (ZTMP(1,I), ZW(I,1)) ! For planar mode
                     END IF
                     ZW(I,1) = ZTMP(1,I) ! Includes missing line aft of tr. edge
                  END IF
               END IF

            END DO ! Next J
         END DO ! Next I


C        Generate the computational mesh on the body by redistributing the arcs.
C        A FIXBGRID-like scheme is retrofitted here to avoid separate inputs.

         DO I = 1, IL

            CALL CHORDS3D (NREFN, XTMP(1,I), YTMP(1,I), ZTMP(1,I),
     >                     .FALSE., SFTOT, STMP)

C           Unnormalized one-sided stretching in TRADIAL(1:JCR):

            IF (I <= ITL) THEN
               R  = REAL (I - 1) / REAL (ITL - 1)
               D1 = (ONE - R) * DRADIAL(1,3) + R * DRADIAL(1,2)
            ELSE IF (I < ITU) THEN
               R  = ABS (REAL (I - ILE) / REAL (ITU - ILE))
               D1 = (ONE - R) * DRADIAL(1,1) + R * DRADIAL(1,2)
            ELSE ! I >= ITU
               R  = REAL (I - ITU) / REAL (IL - ITU)
               D1 = (ONE - R) * DRADIAL(1,2) + R * DRADIAL(1,3)
            END IF

            DUNIFORM = 0.98 * RJCRM1 * SFTOT
            D1 = MIN (D1, DUNIFORM) ! Safeguards EXPDIS4

            CALL EXPDIS4 (3, ZERO, SFTOT, D1, JCR, TRADIAL, BETA,-IWRIT)

C           Revert to two-sided if the outer increment is too big:

            D2 = SFTOT - TRADIAL(JCR-1)
            D2MAX = DUNIFORM * 1.8
            IF (D2 > D2MAX) THEN

               CALL HTDIS4 (.TRUE., ZERO, SFTOT, D1, D2MAX, JCR,
     >                      TRADIAL, -IWRIT, IER)

               IF (IER /= 0 .AND. IER /= 3) THEN
                  WRITE (IWRIT,*) 'Z0PLANE: HTDIS4 trouble. IER: ', IER
                  WRITE (IWRIT,*) 'I, D1, D2, SFTOT: ', I,D1,D2MAX,SFTOT
                  GO TO 990
               END IF
            END IF

C           Linear interpolation for desired distribution along the arc:

            ZWALL(I,1) = ZTMP(1,I)
            JJ1 = 2
            DO J = 2, JCR - 1
               SINT = TRADIAL(J)
               DO JJ = JJ1, NREFN
                  L = JJ
                  IF (SINT <= STMP(JJ)) EXIT
               END DO
               JJ1 = L
               R = (SINT - STMP(L-1)) / (STMP(L) - STMP(L-1))
               XWALL(I,J) = (ONE - R) * XTMP(L-1,I) + R * XTMP(L,I)
               YWALL(I,J) = (ONE - R) * YTMP(L-1,I) + R * YTMP(L,I)
               ZWALL(I,J) = (ONE - R) * ZTMP(L-1,I) + R * ZTMP(L,I)
            END DO
         END DO

         DEALLOCATE (XTMP, YTMP, ZTMP, STMP, XWATER)

      ELSE

C        Parametric method for the body grid:
C        -----------------------------------

         DO I = 1, NFSTN
            XBODY(1:NJBODY,I) = XF(I)
         END DO

C        Parameterize the body surface, fudging the singular nose point:

         YMID  = (YF(1,2) + YF(NJBODY,2)) * HALF
         YNOSE = YF(1,1)
         DO J = 1, NJBODY
            YF(J,1) = 0.01 * (YF(J,2) - YMID) + YNOSE
            ZF(J,1) = 0.01 *  ZF(J,2)
         END DO

         UBODY(1,1) = -999.   ! Kludge to avoid normalization

         CALL PARAM2D (NJBODY, MXIBODY, 1, NJBODY, 1, NFSTN,
     >                 XBODY, YF, ZF, UBODY, VBODY)

         UVRANGE(1) =  1.E+32 ! Umin
         UVRANGE(2) = -1.E+32 ! Umax
         UVRANGE(3) =  ZERO   ! Vmin
         UVRANGE(4) = -1.E+32 ! Vmax

         DO J = 1, NJBODY
            UVRANGE(4) = MAX (UVRANGE(4), VBODY(J,NFSTN))
         END DO

C        Equilibrate the V lines to avoid a ragged edge at I = NFSTN:

         DO J = 1, NJBODY
            R = UVRANGE(4) / VBODY(J,NFSTN)
            DO I = 1, NFSTN - 1
               VBODY(J,I) = VBODY(J,I) * R
            END DO
            VBODY(J,NFSTN) = UVRANGE(4)
         END DO

C        Center the arc lengths about U = 0, determine the U range,
C        and set up the crown/keel line parameterization for searching.

         DO I = 1, NFSTN
            USHIFT = UBODY (NJBODY,I) * HALF
            DO J = 1, NJBODY
               UBODY(J,I) = UBODY(J,I) - USHIFT
            END DO
            UVRANGE(1) = MIN (UVRANGE(1), UBODY(1,I))
            UVRANGE(2) = MAX (UVRANGE(2), UBODY(NJBODY,I))
            VCROWN(I)  = VBODY(NJBODY,I)
            VKEEL(I)   = VBODY(1,I)
         END DO

C***         write (50) NJBODY, NFSTN
C***         write (50) ((VBODY(J,I), J = 1, NJBODY), I = 1, NFSTN),
C***     >              ((UBODY(J,I), J = 1, NJBODY), I = 1, NFSTN)

C        Set up the crown-line grid boundary (U,V)s:

         UWALL(ILE,JCR) = UBODY(NJBODY/2,1)
         VWALL(ILE,JCR) = ZERO

         NCROWN = IBU - ILE + 1
         CALL CHORDS2D (NCROWN, XCROWN(ILE), YCROWN(ILE), .FALSE.,
     >                  VTOTAL, TFIWORK)

         R = UVRANGE(4) / VTOTAL ! In view of equilibration above
         J = 1
         ISTN = 1

         DO I = ILE + 1, IBU - 1
            J = J + 1
            VCROWNI = TFIWORK(J) * R
            VWALL(I,JCR) = VCROWNI

            CALL INTERVAL (NFSTN, VCROWN, VCROWNI, ONE, ISTN)

            UWALL(I,JCR) = UBODY(NJBODY,ISTN) +
     >         (VCROWNI - VCROWN(ISTN)) *
     >         (UBODY(NJBODY,ISTN+1) - UBODY(NJBODY,ISTN)) /
     >         (VCROWN(ISTN+1) - VCROWN(ISTN))
         END DO

         UWALL(IBU,JCR) = UBODY(NJBODY,NFSTN)
         VWALL(IBU,JCR) = VBODY(NJBODY,NFSTN)

C        Likewise for the keel-line grid boundary (U,V)s:

         CALL CHORDS2D (NCROWN, XCROWN(IBL), YCROWN(IBL), .FALSE.,
     >                  VTOTAL, TFIWORK)

         R = UVRANGE(4) / VTOTAL

         DO I = 1, NCROWN
            TFIWORK(I) = VTOTAL - TFIWORK(I) ! Arcs increase in wrong direction
         END DO

         J = NCROWN
         ISTN = 1

         DO I = ILE - 1, IBL + 1, -1
            J = J - 1
            VKEELI = TFIWORK(J) * R
            VWALL(I,JCR) = VKEELI

            CALL INTERVAL (NFSTN, VKEEL, VKEELI, ONE, ISTN)

            UWALL(I,JCR) = UBODY(1,ISTN) + (VKEELI - VKEEL(ISTN)) *
     >         (UBODY(1,ISTN+1) - UBODY(1,ISTN)) /
     >         (VKEEL(ISTN+1) - VKEEL(ISTN))
         END DO

         UWALL(IBL,JCR) = UBODY(1,NFSTN)
         VWALL(IBL,JCR) = VBODY(1,NFSTN)

C        (U,V)s of the inner boundary of the C-mesh on the body are harder.
C        We don't have Zs between the trailing edge and the tail either.
C        The following finds (U,V) as a by-product of finding Z for given (X,Y).

         IS = NFSTN - 1

         CALL INTERVAL (NFSTN, XF, XW(IBL+1,1), ONE, IS)

         JS = NJBODY / 2
         UCUT = UBODY(JS,IS)
         VCUT = VBODY(JS,IS)
         DRANGE = XW(ITL,1) - XW(ILE,1)
         EPS = MAX (5.* EPSMCH, TOLER) ! Scaled by DRANGE
         ISMAT = -1
         JSMAT = -1

         DO I = IBL + 1, IBU - 1
            XYZCUT(1) = XW(I,1)
            XYZCUT(2) = YW(I,1)

            IF (IDEGBOD == 3) THEN

               CALL PLBICUT (1, NJBODY, MXIBODY, 1, NJBODY, 1, NFSTN,
     >                       XBODY, YF, ZF, UBODY, VBODY, JS, IS,
     >                       JSMAT, ISMAT, SMATRIX, EPS, DRANGE, XYZCUT,
     >                       UVRANGE, UCUT, VCUT, P, Q, -IWRIT, IER)
            ELSE

               CALL BILINCUT (1, NJBODY, MXIBODY, 1, NJBODY, 1, NFSTN,
     >                       XBODY, YF, ZF, UBODY, VBODY, JS, IS,
     >                       EPS, DRANGE, XYZCUT,
     >                       UVRANGE, UCUT, VCUT, P, Q, -IWRIT, IER)
            END IF

            IF (IER > 2) THEN
               WRITE (IWRIT, '(/,A,I2,/,A,I5,3F10.4,/,A,2I5,2F13.6)')
     >            ' Z0PLANE: (u,v) boundary error at root =', IER,
     >            ' I, X, Y, Z:', I, XYZCUT,
     >            ' IS, JS, U, V:', IS, JS, UCUT, VCUT
            END IF

            UWALL(I,1) = UCUT
            VWALL(I,1) = VCUT
            ZWALL(I,1) = XYZCUT(3)
         END DO

C        Ensure common points aft of the trailing edge:

         J = ITL
         DO I = ITU, IBU - 1
            UCUT = (UWALL(J,1) + UWALL(I,1)) * HALF
            UWALL(J,1) = UCUT
            UWALL(I,1) = UCUT
            VCUT = (VWALL(J,1) + VWALL(I,1)) * HALF
            VWALL(J,1) = VCUT
            VWALL(I,1) = VCUT
            ZCUT = (ZWALL(J,1) + ZWALL(I,1)) * HALF
            ZWALL(J,1) = ZCUT
            ZWALL(I,1) = ZCUT
            ZW(J,1)    = ZCUT
            ZW(I,1)    = ZCUT
            J = J - 1
         END DO

C        Special treatment of the (U,V) boundary grid at the last body station.
C        First, locate the value of U at the vortex sheet cut more reliably
C        than by PLBICUT above.  Output UCUT <-> YW(IBL,1).

         JCUT = NJBODY / 2

         CALL PLXCUT (NJBODY, YF(1,NFSTN), ZF(1,NFSTN), UBODY(1,NFSTN),
     >                JCUT, .FALSE., YW(IBL,1), EPS, UCUT, -IWRIT, IER)

         IF (IER > 1) THEN
            WRITE (IWRIT,'(A,E15.6)') ' Z0PLANE: Target Y =', YW(IBL,1)
            STOP
         END IF

         UWALL(IBL,1) = UCUT
         UWALL(IBU,1) = UCUT

C        We need V at the ends of the inner C boundary, but it will not
C        be used by PLBICUBE - just by the TFI and smoothing.
C        (If V lines are equilibrated above, VCUT will just be Vmax.)

         CALL LCSFIT (4, UBODY(JCUT-1,NFSTN), VBODY(JCUT-1,NFSTN),
     >                NEWDAT, 'B', 1, UCUT, VCUT, VCUT)

         VWALL(IBL,1) = VCUT
         VWALL(IBU,1) = VCUT

C        YWALL here is done below via PLSCRV3D.

C        Uniform distribution in U (circumferential):

         DUL = (UWALL(IBL,JCR) - UCUT) * RJCRM1
         DUU = (UWALL(IBU,JCR) - UCUT) * RJCRM1

C        Corresp. tail boundary V distribution for the interim grid (uniform):

         DO J = 2, JCR - 1
            RJM1 = REAL (J - 1)
            UWALL(IBL,J) = UWALL(IBL,1) + RJM1 * DUL

            CALL LCSFIT (NJBODY, UBODY(1,NFSTN), VBODY(1,NFSTN), NEWDAT,
     >                'B', 1, UWALL(IBL,J), VWALL(IBL,J), VWALL(IBL,J))

            UWALL(IBU,J) = UWALL(IBU,1) + RJM1 * DUU

            CALL LCSFIT (NJBODY, UBODY(1,NFSTN), VBODY(1,NFSTN), NEWDAT,
     >                'B', 1, UWALL(IBU,J), VWALL(IBU,J), VWALL(IBU,J))
         END DO


C        Initialize interior (U,V)s via transfinite interpolation.
C        Introduce a boundary between the nose and the root leading edge.
C        Always use the initial radial increment - FIXBGRID imposes the final
C        radial distribution:

         D1 = D1NOSEEU
         D2 = MIN (D1 * HALF, 0.005 * CHORDM(1))
         V1 = VWALL(ILE,JCR) ! Actually 0.
         V2 = VWALL(ILE,1)

         CALL HTDIS4 (.TRUE., V1, V2, D1, D2, JCR, TFIWORK, -IWRIT,
     >                IER)

         IF (IER /= 0 .AND. IER /= 3) THEN
            WRITE (IWRIT, '(/,A,I2,/,A,4F13.6)')
     >         ' Z0PLANE: IER from HTDIS4 =', IER,
     >         ' V1, V2, D1, D2: ', V1, V2, D1, D2
            GO TO 990
         END IF

         RATIO = (UWALL(ILE,JCR) - UWALL(ILE,1)) / (V1 - V2)
         DO J = 2, JCR - 1
            VWALL(ILE,J) = TFIWORK(JCR + 1 - J)
            UWALL(ILE,J) = UWALL(ILE,1) + RATIO * (VWALL(ILE,J) - V2)
         END DO

C        Introduce 2 more boundaries at the trailing edge (uniform in u):

         DUL = (UWALL(ITL,JCR) - UWALL(ITL,1)) * RJCRM1
         DUU = (UWALL(ITU,JCR) - UWALL(ITU,1)) * RJCRM1
!!!      DVL = (VWALL(ITL,JCR) - VWALL(ITL,1)) * RJCRM1
!!!      DVU = (VWALL(ITU,JCR) - VWALL(ITU,1)) * RJCRM1

         RLW = VWALL(ITL,1) / VCUT ! V ratio at lower water line
         RUW = VWALL(ITU,1) / VCUT ! .. and upper (= RLW for closed t.e.)
         RLC = VWALL(ITL,JCR) / VWALL(IBL,JCR) ! V ratio at lower crown pt.
         RUC = VWALL(ITU,JCR) / VWALL(IBU,JCR) ! ... and at upper ...

         DO J = 2, JCR - 1
            RJM1 = REAL (J - 1)
            UWALL(ITL,J) = UWALL(ITL,1) + RJM1 * DUL
            UWALL(ITU,J) = UWALL(ITU,1) + RJM1 * DUU
!!!         VWALL(ITL,J) = VWALL(ITL,1) + RJM1 * DVL ! Can give irregular
!!!         VWALL(ITU,J) = VWALL(ITU,1) + RJM1 * DVU ! line in XYZ space
            R = RJM1 * RJCRM1
            VWALL(ITL,J) = VWALL(IBL,J) * (R*RLC + (ONE - R)*RLW)
            VWALL(ITU,J) = VWALL(IBU,J) * (R*RUC + (ONE - R)*RUW)
         END DO

         CALL TFI2D (IL, IBL, ITL, 1, JCR, UWALL, VWALL, TFIWORK)
         CALL TFI2D (IL, ITL, ILE, 1, JCR, UWALL, VWALL, TFIWORK)
         CALL TFI2D (IL, ILE, ITU, 1, JCR, UWALL, VWALL, TFIWORK)
         CALL TFI2D (IL, ITU, IBU, 1, JCR, UWALL, VWALL, TFIWORK)

C        Smooth the (U,V)s.  Note that spacing control from J1 bunches
C        the radial lines too much between the root LE and the nose.

         L = 4

         IF (PRINT2D)
     >      WRITE (IWRIT,'(/,A)') ' ELLIP2D on body (u,v)s:'

         CALL ELLIP2D (IL, JL, IBL, ITL, 1, JCR,
     >                 VWALL,       UWALL,       DFACEK,   PRINT2D,
     >                 BGMODE2D(L), FGMODE2D(L), SPMODE(L),
     >                 ITMAX2D(L),  ITFLOAT,     ITFREEZE, JLAYER,
     >                 CONV2D(L),   CONVMIN(L),  DMAX2D,   OMG2D(L),
     >                 POWERI(L),   POWERJ(L),   EXPI12D(L),
     >                 EXPI22D(L),  EXPJ12D(L),  EXPJ22D(L),
     >                 FGLIM2D,     URFG2D,      URFLOAT)

         L = 5

         CALL ELLIP2D (IL, JL, ITL, ILE, 1, JCR,
     >                 VWALL,       UWALL,       DFACEK,   PRINT2D,
     >                 BGMODE2D(L), FGMODE2D(L), SPMODE(L),
     >                 ITMAX2D(L),  ITFLOAT,     ITFREEZE, JLAYER,
     >                 CONV2D(L),   CONVMIN(L),  DMAX2D,   OMG2D(L),
     >                 POWERI(L),   POWERJ(L),   EXPI12D(L),
     >                 EXPI22D(L),  EXPJ12D(L),  EXPJ22D(L),
     >                 FGLIM2D,     URFG2D,      URFLOAT)

         CALL EXCHANGE (BGMODE2D(L), BGMODE)
         CALL EXCHANGE (FGMODE2D(L), FGMODE)

         CALL ELLIP2D (IL, JL, ILE, ITU, 1, JCR,
     >                 VWALL,       UWALL,       DFACEK,   PRINT2D,
     >                 BGMODE,      FGMODE,      SPMODE(L),
     >                 ITMAX2D(L),  ITFLOAT,     ITFREEZE, JLAYER,
     >                 CONV2D(L),   CONVMIN(L),  DMAX2D,   OMG2D(L),
     >                -POWERI(L),   POWERJ(L),   EXPI22D(L),
     >                 EXPI12D(L),  EXPJ12D(L),  EXPJ22D(L),
     >                 FGLIM2D,     URFG2D,      URFLOAT)

         L = 4

         CALL EXCHANGE (BGMODE2D(L), BGMODE)
         CALL EXCHANGE (FGMODE2D(L), FGMODE)

         CALL ELLIP2D (IL, JL, ITU, IBU, 1, JCR,
     >                 VWALL,       UWALL,       DFACEK,   PRINT2D,
     >                 BGMODE,      FGMODE,      SPMODE(L),
     >                 ITMAX2D(L),  ITFLOAT,     ITFREEZE, JLAYER,
     >                 CONV2D(L),   CONVMIN(L),  DMAX2D,   OMG2D(L),
     >                -POWERI(L),   POWERJ(L),   EXPI22D(L),
     >                 EXPI12D(L),  EXPJ12D(L),  EXPJ22D(L),
     >                 FGLIM2D,     URFG2D,      URFLOAT)

C***     write (51) IBU - IBL + 1, JCR
C***     write (51) ((VWALL(I,J), I = IBL, IBU), J = 1, JCR),
C*** >              ((UWALL(I,J), I = IBL, IBU), J = 1, JCR)

C        Interpolate the body surface at the smoothed (u,v) grid points:

         ISMAT = -1
         JSMAT = -1
         JS = NJBODY / 2

         DO J = 1, JCR - 1
            IS = NFSTN - 1
            JS = NJBODY - JS

            DO I = IBL + 1, IBU - 1 ! IBL, IBU are handled above

               IF (J == 1) THEN ! (U,V) was averaged above along cut
                  IF (I <= ITL .OR. I >= ITU) CYCLE ! Next I
               END IF

               IF (IDEGBOD == 3) THEN

                 CALL PLBICUBE (0, NJBODY, MXIBODY, 1, NJBODY, 1, NFSTN,
     >                          XBODY, YF, ZF, UBODY, VBODY,
     >                          UWALL(I,J), VWALL(I,J), JS, IS,
     >                          JSMAT, ISMAT, SMATRIX, EPS, P, Q,
     >                          XG, YG, ZG, XYZCUT, XYZCUT, IER)
               ELSE

                  CALL PBILINT (0, NJBODY, MXIBODY, 1, NJBODY, 1, NFSTN,
     >                          XBODY, YF, ZF, UBODY, VBODY,
     >                          UWALL(I,J), VWALL(I,J), JS, IS, EPS,
     >                          P, Q, XG, YG, ZG, XYZCUT, XYZCUT, IER)
               END IF

               IF (IER /= 0) THEN
                  WRITE (IWRIT, '(/,A,A,I2,A,4I5,/,A,2F13.5,A,3F13.5)')
     >               ' Z0PLANE: Marginal PLBICUBE or PBILINT result. ',
     >               ' IER:', IER, '; I, J, IS, JS:', I, J, IS, JS,
     >               ' U, V:', UWALL(I,J), VWALL(I,J),
     >               '    X, Y, Z used:', XG, YG, ZG
                  NBAD = NBAD + 1
                  IF (NBAD > JCR * 2) THEN
                     WRITE (IWRIT,'(/,A)')
     >                  ' *** Z0PLANE: Aborting body surface grid.'
                     GO TO 990
                  END IF
               END IF

               XWALL(I,J) = XG
               YWALL(I,J) = YG
               ZWALL(I,J) = ZG

            END DO ! Next I
         END DO ! Next J

C        Special treatment of the XYZ grid on the last body station.
C        UBODY(*,NFSTN) is no longer needed, so normalize it in-place.

         USHIFT = UBODY(NJBODY,NFSTN)
         UTOTAL = ONE / (USHIFT + USHIFT)
         UBODY(1,NFSTN) = ZERO
         DO J = 2, NJBODY - 1
            UBODY(J,NFSTN) = (UBODY(J,NFSTN) + USHIFT) * UTOTAL
         END DO
         UBODY(NJBODY,NFSTN) = ONE

         NEW = .TRUE.
         JCUT = NJBODY / 2

         DO J = 1, JCR - 1
            UCUT = (UWALL(IBL,J) + USHIFT) * UTOTAL

            CALL PLSCRV3D (NJBODY, XBODY(1,NFSTN), YF(1,NFSTN),
     >                     ZF(1,NFSTN), UBODY(1,NFSTN), 'M', NEW,
     >                     .FALSE., UCUT, JCUT, XWALL(IBL,J),
     >                     YWALL(IBL,J), ZWALL(IBL,J), SUPPRESS)
            NEW = .FALSE.
         END DO

C        Likewise for the upper aft boundary of the body C mesh:

         JCUT = NJBODY / 2
         DO J = 1, JCR - 1
            UCUT = (UWALL(IBU,J) + USHIFT) * UTOTAL

            CALL PLSCRV3D (NJBODY, XBODY(1,NFSTN), YF(1,NFSTN),
     >                     ZF(1,NFSTN), UBODY(1,NFSTN), 'M', NEW,
     >                     .FALSE., UCUT, JCUT, XWALL(IBU,J),
     >                     YWALL(IBU,J), ZWALL(IBU,J), SUPPRESS)
         END DO

C        Update the wing/body intersection to match X/Y/ZWALL:

         XW(1:IBU,1) = XWALL(1:IBU,1) ! The upper wake sheet edge is copied
         YW(1:IBU,1) = YWALL(1:IBU,1) ! in VSHEET
         ZW(1:IBU,1) = ZWALL(1:IBU,1)


C        Redistribute the radial lines on the body to ensure precise increments:
C        -----------------------------------------------------------------------

         CALL FIXBGRID (IL, JL, IBL, ITL, ILE, ITU, IBU, JCR,
     >                  DRADIALEU(1,1), DRADIALEU(1,2), D1NOSEEU,
     >                  DWTAILEU, DCTAILEU, NBLAYEREU, RBLAYEREU,
     >                  XWALL, YWALL, ZWALL, IWRIT)

         IF (.NOT. FULLGRID .AND. REDISTRIBUTE) THEN

C           We're warping a grid in indirect mode.  The FIXBGRID call above
C           using *EU controls followed by *NS controls here gives the
C           equivalent of what REGRID_FACES does initially for a full grid.

            CALL FIXBGRID (IL, JL, IBL, ITL, ILE, ITU, IBU, JCR,
     >                     DRADIAL(1,1), DRADIAL(1,2), D1NOSE,
     >                     DWTAIL, DCTAIL, NBLAYER, RBLAYER,
     >                     XWALL, YWALL, ZWALL, IWRIT)
         END IF


C        Retrieve the singular point at the nose:

         YF(1:NJBODY,1) = YNOSE
         ZF(1:NJBODY,1) = ZERO

      END IF

C **      write (52) IBU - IBL + 1, JCR, 1
C **      write (52) ((XWALL(I,J), I = IBL, IBU), J = 1, JCR),
C **     >           ((YWALL(I,J), I = IBL, IBU), J = 1, JCR),
C **     >           ((ZWALL(I,J), I = IBL, IBU), J = 1, JCR)
C **      if (IBL > 1) stop 'Surface grid on unit 52.'


C     Off-body part of the symmetry plane.
C     ------------------------------------

C     Aft-of-body fan:  Transfer the last body station distribution:

      I1 = IBL
      I2 = IBL - 1

      DO N = 1, 2

         CALL CHORDSRF (IL, JL, I1, 1, JCR, XWALL, YWALL, ZWALL,
     >                  .TRUE., SFTOT, TRADIAL)

         XRANGE = XWALL(I2,1) - XWALL(I2,JCR)
         YRANGE = YWALL(I2,1) - YWALL(I2,JCR)

         DO J = 2, JCR - 1
            XWALL(I2,J) = XWALL(I2,1) - XRANGE * TRADIAL(J)
            YWALL(I2,J) = YWALL(I2,1) - YRANGE * TRADIAL(J)
         END DO

         I1 = IBU
         I2 = IBU + 1

      END DO


      IF (FULLGRID) THEN

         I1 = 1
         I2 = IBL - 1
         L  = 6

         DO N = 1, 2 ! Fan halves

            CALL TFI2D (IL, I1, I2, 1, JCR, XWALL, YWALL, TFIWORK)

            IF (ITMAX2D(L) > 0) THEN

               IF (N == 1) THEN
                  BGMODE = BGMODE2D(L)
                  FGMODE = FGMODE2D(L)
                  PWRI   = POWERI(L)
                  IF (PRINT2D) WRITE (IWRIT,'(/,A)')
     >               ' ELLIP2D on aft-of-body "fan":'
               ELSE
                  CALL EXCHANGE (BGMODE2D(L), BGMODE)
                  CALL EXCHANGE (FGMODE2D(L), FGMODE)
                  PWRI = -PWRI
               END IF

               CALL ELLIP2D (IL, JL, I1, I2, 1, JCR,
     >                       XWALL, YWALL, DFACEK, PRINT2D,
     >                       BGMODE,      FGMODE,     SPMODE(L),
     >                       ITMAX2D(L),  ITFLOAT,    ITFREEZE, JLAYER,
     >                       CONV2D(L),   CONVMIN(L), DMAX2D, OMG2D(L),
     >                       PWRI,        POWERJ(L),  EXPI12D(L),
     >                       EXPI22D(L),  EXPJ12D(L), EXPJ22D(L),
     >                       FGLIM2D,     URFG2D,     URFLOAT)
            END IF

            I1 = IBU + 1
            I2 = IL

         END DO

C        Line forward of the nose:

         CALL OFFCROWN (ILE, IL, JL, JCR, D1CROWN, D2JL, TRADIAL,
     >                  XWALL, YWALL, IWRIT, FAIL)
         IF (FAIL) GO TO 990

         L = 1 ! Off-body symmetry plane

         DO N = 1, 3, 2 ! ... in two halves, not quadrants

            I1 = IFACE(N)
            I2 = IFACE(N+2)

            CALL TFI2D (IL, I1, I2, JCR, JL, XWALL, YWALL, TFIWORK)

            IF (N == 1) THEN
               BGMODE = BGMODE2D(L)
               FGMODE = FGMODE2D(L)
               PWRI   = POWERI(L)
               IF (PRINT2D) WRITE (IWRIT,'(/,A)')
     >            ' ELLIP2D on off-body symmetry plane:'
            ELSE
               CALL EXCHANGE (BGMODE2D(L), BGMODE)
               CALL EXCHANGE (FGMODE2D(L), FGMODE)
               PWRI = -PWRI
            END IF

            CALL ELLIP2D (IL, JL, I1, I2, JCR, JL,
     >                    XWALL,      YWALL,      DFACEK,   PRINT2D,
     >                    BGMODE,     FGMODE,     SPMODE(L),
     >                    ITMAX2D(L), ITFLOAT,    ITFREEZE, JLAYER,
     >                    CONV2D(L),  CONVMIN(L), DMAX2D,   OMG2D(L),
     >                    PWRI,       POWERJ(L),  EXPI12D(L),
     >                    EXPI22D(L), EXPJ12D(L), EXPJ22D(L),
     >                    FGLIM2D,    URFG2D,     URFLOAT)
         END DO

      ELSE ! Perturb the off-body symmetry plane grid

         I1 = 1     ! Aft fan, lower part first
         I2 = IBL - 1

         DO N = 1, 2

            CALL WARPQ3D (1, IL, 1, JL, 1, KL, I1, I2, 1, JCR, 1, 1,
     >                    X0(1,1,1,1), X0(1,1,1,2), X0(1,1,1,3), S0,
     >                    DFACEI, DFACEJ, DFACEK, XWALL, YWALL, ZWALL)
            I1 = IBU + 1
            I2 = IL

         END DO

C        Outside the crown/keel line in two halves:

C        Warp the line forward of the nose:

         XRADIAL(JCR:JL) = X0(ILE,JCR:JL,1,1)
         YRADIAL(JCR:JL) = X0(ILE,JCR:JL,1,2)
         J = JL - JCR
         XLINE(JCR:JL:J) = XWALL(ILE,JCR:JL:J)
         YLINE(JCR:JL:J) = YWALL(ILE,JCR:JL:J)

         CALL NULINE2D (JCR, JL, XRADIAL, YRADIAL, XLINE, YLINE)

         XWALL(ILE,JCR:JL) = XLINE(JCR:JL)
         YWALL(ILE,JCR:JL) = YLINE(JCR:JL)

         DO N = 1, 3, 2 ! Quadrants

            I1 = IFACE(N)
            I2 = IFACE(N+2)

            CALL WARPQ3D (1, IL, 1, JL, 1, KL, I1, I2, JCR, JL, 1, 1,
     >                    X0(1,1,1,1), X0(1,1,1,2), X0(1,1,1,3), S0,
     >                    DFACEI, DFACEJ, DFACEK, XWALL, YWALL, ZWALL)
         END DO

      END IF

      GO TO 999


  800 CONTINUE

C**********************************************************************
C     For wing-alone configurations the following generates the 2-D
C     grid at the symmetry plane.
C**********************************************************************

      XWALL(1, 1:JL) = XLBACK(1:JL,1) ! I edges
      XWALL(IL,1:JL) = XUBACK(1:JL,1)
      YWALL(1 ,1:JL) = YLBACK(1:JL,1)
      YWALL(IL,1:JL) = YUBACK(1:JL,1)

      XWALL(1:IL,1)  = XW(1:IL,1)     ! J edges
      XWALL(1:IL,JL) = XOUT
      YWALL(1:IL,1)  = YW(1:IL,1)
      YWALL(1:IL,JL) = YOUT

C     Edges off TE: OFFTEDGE works with X(*,*,*,*), not X/Y/ZWALL, so ...

      J = JL - 1
      X(1:IL,1:JL:J,1,1) = XWALL(1:IL,1:JL:J)
      X(1:IL,1:JL:J,1,2) = YWALL(1:IL,1:JL:J)
      X(1:IL,1:JL:J,1,3) = ZERO

      X(1, 2:J,1,1) = XWALL(1, 2:J)
      X(IL,2:J,1,1) = XWALL(IL,2:J)
      X(1, 2:J,1,2) = YWALL(1, 2:J)
      X(IL,2:J,1,2) = YWALL(IL,2:J)
      X(1, 2:J,1,3) = ZERO
      X(IL,2:J,1,3) = ZERO

      DO I = ITL, ITU, ITU - ITL

         CALL OFFTEDGE (I, 1, DRADIAL(1,2), D2JL, NBLAYER, RBLAYER,
     >                  TRADIAL, IWRIT)
      END DO

      XWALL(ITL,2:J) = X(ITL,2:J,1,1)
      XWALL(ITU,2:J) = X(ITU,2:J,1,1)
      YWALL(ITL,2:J) = X(ITL,2:J,1,2)
      YWALL(ITU,2:J) = X(ITU,2:J,1,2)

C     Line forward of the root leading edge (K = 1):

      CALL OFFLEDGE (IL, JL, KL, ILE, 1, X, DRADIALEU(1,1),
     >               D2JLEU, NBLAYEREU, RBLAYEREU, TRADIAL,
     >               XRADIAL, YRADIAL, ZRADIAL, IWRIT)

      IF (.NOT. FULLGRID .AND. REDISTRIBUTE) THEN

         CALL REGRID_LINE (IL, JL, KL, ILE, 1, X, DRADIAL(1,1),
     >                     D2JL, NBLAYER, RBLAYER, TRADIAL,
     >                     XRADIAL, YRADIAL, ZRADIAL, IWRIT)
      END IF

      XWALL(ILE,1:JL) = XRADIAL
      YWALL(ILE,1:JL) = YRADIAL
      ZWALL(ILE,1:JL) = ZERO


      IF (FULLGRID) THEN

C        TFI and smoothing of the symmetry plane in four parts:

         IF (PRINT2D)
     >      WRITE (IWRIT,'(/,A)') ' ELLIP2D on symmetry plane:'

         CALL TFI2D (IL, 1, ITL,  1, JL, XWALL, YWALL, TFIWORK)

         L = 2 ! Regions 1 and 4

         CALL ELLIP2D (IL, JL, 1, ITL, 1, JL,
     >                 XWALL,       YWALL,       DFACEK,    PRINT2D,
     >                 BGMODE2D(L), FGMODE2D(L), SPMODE(L),
     >                 ITMAX2D(L),  ITFLOAT,     ITFREEZE,  JLAYER,
     >                 CONV2D(L),   CONVMIN(L),  DMAX2D,    OMG2D(L),
     >                 POWERI(L),   POWERJ(L),   EXPI12D(L),
     >                 EXPI22D(L),  EXPJ12D(L),  EXPJ22D(L),
     >                 FGLIM2D,     URFG2D,      URFLOAT)

         CALL TFI2D (IL, ITL, ILE, 1, JL, XWALL, YWALL, TFIWORK)

         L = 7 ! Regions 2 and 3

         CALL ELLIP2D (IL, JL, ITL, ILE, 1, JL,
     >                 XWALL,       YWALL,       DFACEK,    PRINT2D,
     >                 BGMODE2D(L), FGMODE2D(L), SPMODE(L),
     >                 ITMAX2D(L),  ITFLOAT,     ITFREEZE,  JLAYER,
     >                 CONV2D(L),   CONVMIN(L),  DMAX2D,    OMG2D(L),
     >                 POWERI(L),   POWERJ(L),   EXPI12D(L),
     >                 EXPI22D(L),  EXPJ12D(L),  EXPJ22D(L),
     >                 FGLIM2D,     URFG2D,      URFLOAT)

         CALL TFI2D (IL, ILE, ITU, 1, JL, XWALL, YWALL, TFIWORK)

         CALL EXCHANGE (BGMODE2D(L), BGMODE)
         CALL EXCHANGE (FGMODE2D(L), FGMODE)

         CALL ELLIP2D (IL, JL, ILE, ITU,  1, JL,
     >                 XWALL,       YWALL,       DFACEK,    PRINT2D,
     >                 BGMODE,      FGMODE,      SPMODE(L),
     >                 ITMAX2D(L),  ITFLOAT,     ITFREEZE,  JLAYER,
     >                 CONV2D(L),   CONVMIN(L),  DMAX2D,    OMG2D(L),
     >                -POWERI(L),   POWERJ(L),   EXPI22D(L),
     >                 EXPI12D(L),  EXPJ12D(L),  EXPJ22D(L),
     >                 FGLIM2D,     URFG2D,      URFLOAT)

         CALL TFI2D (IL, ITU, IL, 1, JL, XWALL, YWALL, TFIWORK)

         L = 2
         CALL EXCHANGE (BGMODE2D(L), BGMODE)
         CALL EXCHANGE (FGMODE2D(L), FGMODE)

         CALL ELLIP2D (IL, JL, ITU, IL, 1, JL,
     >                 XWALL,       YWALL,       DFACEK,    PRINT2D,
     >                 BGMODE,      FGMODE,      SPMODE(L),
     >                 ITMAX2D(L),  ITFLOAT,     ITFREEZE,  JLAYER,
     >                 CONV2D(L),   CONVMIN(L),  DMAX2D,    OMG2D(L),
     >                -POWERI(L),   POWERJ(L),   EXPI22D(L),
     >                 EXPI12D(L),  EXPJ12D(L),  EXPJ22D(L),
     >                 FGLIM2D,     URFG2D,      URFLOAT)

      ELSE ! Perturb the wing-only symmetry plane

         DO N = 1, 4 ! Quadrants

            I1 = IFACE(N)
            I2 = IFACE(N+1)

            CALL WARPQ3D (1, IL, 1, JL, 1, KL, I1, I2, 1, JL, 1, 1,
     >                    X0(1,1,1,1), X0(1,1,1,2), X0(1,1,1,3), S0,
     >                    DFACEI, DFACEJ, DFACEK, XWALL, YWALL, ZWALL)
         END DO

      END IF

      GO TO 999


C     Fatal error handling:

  990 CONTINUE

      FAIL = .TRUE. ! Allow for retrying with adjusted design variables
                    ! (design code only)
  999 RETURN

      END SUBROUTINE Z0PLANE
