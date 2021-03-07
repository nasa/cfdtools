!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      MODULE GEOM0
!
!     AEROSURF module GEOM0 is used to match geometry-related work-space to the
!     problem size.  In order to avoid reading the geometry more than once, the
!     strategy is to allocate interim arrays (GEOM0) that should be more than
!     large enough, read the data, allocate the more precise work-space (GEOM1),
!     copy the data, and release the interim arrays.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     The following constants should err on the large side.  They enable interim
!     work-space to be allocated for the input geometry then (soon) released.

      INTEGER, PARAMETER ::
     >   MXWSURF = 8,   ! Max. # "wing" surfaces (wing/pylon/nac./canard/tail)
     >   MXIWING = 600, ! Max. # geometry points per wing-type defining section
     >   MXKWING = 180, !  "   "   "   "  defining sections per "wing" surface
     >   MXIBODY = 251, !  "   "   "   "   "   "   "   "    for the fuselage
     >   MXJBODY = 261, !  "   "   "   "  points per fuselage defining section
     >   MXCRANK = 6    !  "   " crank sections per wing-type component (+ tip)

      INTEGER, ALLOCATABLE ::
     >   IQUAD0(:,:), KCRANK0(:,:), NFJ0(:)

      INTEGER, ALLOCATABLE, DIMENSION (:) ::
     >   ILWING0, KLWING0, IDEGREE0, IDIHED0, IFLAP0, IROLL0, NWSEC0,
     >   KWING10, KWING20, LCRANK0, NCRANK0,
     >   NFIXA0,  NFIXB0,  NFIXL0,  NFIXR0, NNROOT0, NNTAIL0, NNTIP0

      INTEGER, ALLOCATABLE, DIMENSION (:,:) ::
     >   NLG0, NUG0, NWG0

      REAL, ALLOCATABLE, DIMENSION (:) ::
     >   CRDSPL0, CRDSPQ0, CRDSPS0, CRDSPC0, SPNSPC0, SPNSPS0,
     >   XC10, YC10, ZC10, XC20, YC20, ZC20, CLOCK0, XTRAP0,
     >   ZINNER0, ZOUTER0

      REAL, ALLOCATABLE, DIMENSION (:,:) ::
     >   AINIT0, ALG0, CINIT0, DFLAP0, DIHED0, DSLAT0, ROLL0, ROUNDED0,
     >   SWFLAP0, SWSLAT0, TINIT0, XFLAP0, XSLAT0, XLEADI0, YLEADI0,
     >   ZLEADI0, ZCRANK0

      REAL, ALLOCATABLE, DIMENSION (:) ::
     >   XF0

      REAL, ALLOCATABLE, DIMENSION (:,:) ::
     >   YFINIT0, ZFINIT0

      REAL, ALLOCATABLE, DIMENSION (:,:,:) ::
     >   XINIT0, YINIT0, ZINIT0, XWG0, YWG0, ZWG0

      REAL, ALLOCATABLE, DIMENSION (:,:,:,:) ::
     >   DIAXIS0

      CHARACTER * 32, ALLOCATABLE, DIMENSION (:) ::
     >   PIECE0

      END MODULE GEOM0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      MODULE GEOM1
!
!     Module GEOM1 arrays should parallel those in module GEOM0.  They will be
!     sized more precisely to the problem at run time.
!
!     EXCEPTION:  NLGI, NUGI, NWGI (no need to have these copies in GEOM0).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      INTEGER, ALLOCATABLE ::
     >   IQUAD(:,:), KCRANK(:,:), NFJ(:)

      INTEGER, ALLOCATABLE, DIMENSION (:) ::
     >   ILWING, KLWING, IDEGREE, IDIHED, IFLAP, IROLL, NWSEC,
     >   KWING1, KWING2, LCRANK, NCRANK,
     >   NFIXA,  NFIXB,  NFIXL,  NFIXR, NNROOT, NNTAIL, NNTIP

      INTEGER, ALLOCATABLE, DIMENSION (:,:) ::
     >   NLG,  NUG,  NWG,
     >   NLGI, NUGI, NWGI

      REAL, ALLOCATABLE, DIMENSION (:) ::
     >   CRDSPL, CRDSPQ, CRDSPS, CRDSPC, SPNSPC, SPNSPS,
     >   XC1, YC1, ZC1, XC2, YC2, ZC2, CLOCK, XTRAP, ZINNER, ZOUTER

      REAL, ALLOCATABLE, DIMENSION (:,:) ::
     >   AINIT, ALG, CINIT, DFLAP, DIHED, DSLAT, ROLL, ROUNDED, SWFLAP,
     >   SWSLAT, TINIT, XFLAP, XSLAT, XLEADI, YLEADI, ZLEADI, ZCRANK

      REAL, ALLOCATABLE, DIMENSION (:) ::
     >   XF

      REAL, ALLOCATABLE, DIMENSION (:,:) ::
     >   YFINIT, ZFINIT

      REAL, ALLOCATABLE, DIMENSION (:,:,:) ::
     >   XINIT, YINIT, ZINIT, XWG, YWG, ZWG

      REAL, ALLOCATABLE, DIMENSION (:,:,:,:) ::
     >   DIAXIS

      CHARACTER * 32, ALLOCATABLE, DIMENSION (:) ::
     >   PIECE

      END MODULE GEOM1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      MODULE GRIDS
!
!     Module GRIDS allows precise sizing of arrays for gridded/patched surface.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      INTEGER, ALLOCATABLE, DIMENSION (:) ::
     >   IMXCOM, JMXCOM, IMXCOMO, JMXCOMO, IDIMENS, JDIMENS, IPIECE,
     >   IXB, NSUBI, NSUBJ

      INTEGER, ALLOCATABLE, DIMENSION (:,:) ::
     >   ISUB, JSUB, ISUBO, JSUBO, IINTR, JINTR, KINTR, MCRANK

      REAL, ALLOCATABLE, DIMENSION (:) ::
     >   XBASE, YBASE, ZBASE, XRG, YRG, ZRG

      REAL, ALLOCATABLE, DIMENSION (:,:) ::
     >   SNORM, TINTR, UINTR, VINTR, TBASE, TNORM, YF, ZF

      END MODULE GRIDS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      MODULE SHAPE_FUNCTIONS
!
!     Module SHAPE_FUNCTIONS contains shape function-related variables.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      INTEGER
     >   NDV,
     >   NTWING, NCHCKS, NCHCK1, NCHCK2, NTWIST, NTHICK,
     >   NXLEAD, NYLEAD, NZLEAD, NCHORD, NIROOT, NTOTFU,
     >   NCMBRF, NSPANF, NAREAF, NXFLR,  NYFLR,  NTKINK,
     >   NSINE,  NCOS,   NEXPS,  NTRAIL, NLEAD,  NWAG,
     >   NSINEP, NEXPSP, NTRAIP, NLEADP, NTOTFS,
     >   NSINVC, NFLAP,  NSLAT,  NDIHED,
     >   NSINFC, NCOSFC, NEXPFC, NLEDFC, NTRLFC,
     >   NSINFS, NCOSFS, NEXPFS, NLEDFS, NTRLFS,
     >   NSINFA, NCOSFA, NEXPFA, NLEDFA, NTRLFA

      INTEGER, ALLOCATABLE, DIMENSION (:) ::
     >   KCEN, KMIN, KMAX, IDWING

      REAL, ALLOCATABLE, DIMENSION (:) ::
     >   V, VSCALE, AITCH, DOUP, DOLO,
     >   WBUMP, XBUMP, XWMIN, XWMAX,
     >   PBUMP, UBUMP, UBMIN, UBMAX

      CHARACTER, ALLOCATABLE, DIMENSION (:) ::
     >   UTYPE * 4, VTYPE * 6

C     The following belong with the geometry modules; ignore for now.

C****      INTEGER, ALLOCATABLE, DIMENSION (:) ::
C****     >   IDIVLD, KDIVLD, IDIVTR, KDIVTR,
C****     >   IDIV1,  IDIV2,  KDIV1,  KDIV2

C****      REAL, ALLOCATABLE, DIMENSION (:) ::
C****     >   XFLOORI, YFLOORI,
C****     >   XTVCI1, YTVCI1, ZTVCI1,
C****     >   XTVCI2, YTVCI2, ZTVCI2,
C****     >   XDIVLD, YDIVLD, ZDIVLD,
C****     >   XDIVTR, YDIVTR, ZDIVTR

      END MODULE SHAPE_FUNCTIONS

C*******************************************************************************
C
      PROGRAM AEROSURF
C
C
C        AEROSURF is the stand-alone driver of the geometry paneling modules
C     used by the multiblock design code, SYN107-MB.  It is needed to provide
C     the initial paneling used with an initial surface grid by program UV to
C     parameterize that grid for efficient mesh perturbations during the design
C     process.  It treats geometric entities such as wings, fuselages, pylons,
C     nacelles, etc., by calculating various intersections and imposing regular
C     panels on the intersected surfaces representing a complete or partial
C     aerodynamic configuration.
C
C        Main routine SURFER and the routines it calls should be identical to
C     the corresponding routines in the design software.  AEROSURF is a follow
C     on to WBSURF, which treated just a wing and body.  It saves results in
C     PLOT3D /mgrid form.
C
C
C     AEROSURF Control Inputs (aero.inp file)
C     ---------------------------------------
C
C     Main Controls:
C
C        NNCONFIG   Indicates vehicle configuration:
C                   0 means not a particularly special case - any case
C                     handled prior to CTV (which has overlapping wing + tail);
C                   1 means CTV case (Crew Transfer Vehicle - odd body patches);
C                   2 means more specialized CTV case (affects wing patches);
C                   3 means center-line fin variant of case 2 (off-center fin);
C                   4 means <some future possibility>
C
C        NNWING     Number of wing-type surfaces to be read as defining sections
C        NNNAC      ......... nacelle ..........................................
C        NNFUSE   >=1 means read fuselage sections following wings/nacelles;
C                   0 means do not read a fuselage
C                   2 means write body sections from WBFUDGE, possibly with any
C                     open base (almost) filled in, to aero.body and stop;
C                     format: as for AEROSURF fuselage inputs;
C                     NNFUDGE should be set to the appropriate IW
C                   3 means as for 2 but use PLOT3D format, assuming that the
C                     point counts match; file name: fuselage.plot3d, formatted
C        NNSPLIT    1 means split ALL fuselage panels at all wing leading edges;
C                   2 means split along lines connecting TE and LE points;
C                     NNSPLIT = 1 can give degenerate or varying numbers of
C                     patches if the wing leading edges are closely aligned;
C                     NNSPLIT = 2 gives fewer patches (but more pts.) and should
C                     be better except for extremely misaligned and close wings
C        NNFORMAT   0 means write paneled surface patches unformatted (aero.xyz)
C                   1 ................................... formatted
C                   2 means as for NNFORMAT = 0 but also write the regularized
C                     surfaces to aero.components (unformatted)
C                   3 means as for NNFORMAT = 2, but formatted
C        NNREPEAT   0 means panel the geometry once (normal NCALL = -3 case)
C                   1 means repeat the paneling with NCALL > -3 as a check
C        NNINSERT     permits interpolation of additional defining sections as
C                     may be needed for more effective use of design variables
C                     or (in the case of the fuselage) to assist intersection
C                     calculations:
C                   0 means don't insert new defining sections
C                  +n means insert new sections in wing-type surface n, and stop
C                  -1 means insert new sections in the fuselage, and stop
C        NNKAPPA      permits output of curvature in the I and J directions for
C                     the regularized form of the specified component (or for
C                     radius of curvature if specified at the prompt):
C                   0 means no curvature calculations;
C                  +n means write kappa or rho to aero.ikappa & .jkappa, and the
C                     corresponding  (x,y,z)s  to aero.regular  in either
C                     unformatted or formatted single-function FAST form
C                     according to NNFORMAT, and stop;
C                  -1 means likewise for the regularized fuselage
C
C     Control of subpatch sizes and of packed panel work-space size:
C
C        NNPANEL      controls SURFER's option to adjust patch point counts:
C                   1 means accept whatever point counts are initially produced;
C                   2 means AEROSURF reads subpatch info. from aero.xy0 and
C                     SURFER forces its panels to match those point counts
C
C        NNEXPAND   0 means don't expand subpatches;
C                   1 means expand them after returning from SURFER
C
C        LENXB      Estimated total number of points in the output patches;
C                   if a fuselage is present, the NNSPLIT = 2 option requires
C                   an extra patch temporarily which should be allowed for in
C                   LENXB; any input that is too small will cause termination
C                   with an indication (in aero.out) of how large LENXB needs
C                   to be; err on the large side to begin with, then lower it
C                   to the value reported in aero.out for repeated use.
C
C     Wing-type Surface Controls (including Nacelle Surfaces IW):
C
C        PIECE      Descriptor for this component (up to 32 characters)
C        NNROOT     Surface intersection control:
C        NNROOT(IW) = 0 suppresses an intersection calculation for surface IW
C                   = +N >= 1 means surface IW intersects wing surface N
C                   = -1 means a wing(IW)/body intersection     (if NNFUSE = 1)
C                   = -2 means a vert.fin(IW)/body intersection ( "    "    " )
C                   = -N <= -3 means (tail) surface IW intersects v.fin surf. N
C        NNTAIL     Vertical/horizontal tail association control, meaningful
C                   only if NNROOT(IW) = -2:
C        NNTAIL(IW) = 0 means no horizontal tail (or wing) intersects the body
C                       in a way that overlaps the body intersection of fin IW;
C                   = +N means horizontal tail or wing  N  <-> dorsal  fin IW;
C                   = -N means horizontal tail or wing |N| <-> ventral fin IW
C        ILWING     Streamwise dimension of (plain) paneled surface
C        KLWING     Transverse ............................
C        IROLL      Wing/nacelle section roll control:
C                   1 means geometry sections contain roll angles to be applied
C                   0 means no roll angles for this surface - all are 0.
C        IFLAP      Wing section flap/slat control (not needed for nacelles):
C                   1 means geometry sections contain D/X/SWSLAT & D/X/SWFLAP
C                   0 means no deflections for this wing surface - 0 is assumed
C        IDIHED     Dihedral control:
C                   1 means geometry sections contain dihedral angles/rotn. axes
C                   0 means no dihedral rotations for this surface
C
C        KWING1,    Range of sections to use for the lines during the
C        KWING2     surface/surface intersection calculations
C        IDEGREE    1 or 3 for bilinear or bicubic interpolation when this
C                   surface is the "body" in an intersection calculation;
C                   N.B.: WINTERP also uses this to control the 1-D spanwise
C                   interpolations, and IDEGREE = 1 is normally appropriate
C        XTRAP      X distance beyond each surface edge at which an extra
C                   station is extrapolated to protect the intersection
C                   calculation:
C                   for diverter IW, XTRAP(IW) is used for extrapolating
C                                    the WING trailing edge;
C                   for nacelle IW,  XTRAP(IW)  > 0. means extend the nacelle
C                                    aft, while < 0. means fudge the associated
C                                    DIVERTER from XTE - |XTRAP| to XTE,
C                                    closing it at XTE
C        ZINNER     For XTRAP(IW) < 0. cases, these are the spanwise locations
C        ZOUTER     of the wing trailing edge where the diverter associated with
C                   nacelle IW is forced to close
C
C        NFIXL      For rounded leading edges of wing-type intersections with
C        NFIXR      the fuselage, NFIXL is the number of patched body sections
C                   left of the leading edge to redistribute by temporarily
C                   merging them with NFIXR sections to the right.
C        NFIXA      NFIXA and -B (above and below the leading edge) allow
C        NFIXB      similar control in the vertical direction.
C                   These inputs are inactive for wing/nacelle-type surfaces
C                   not intersecting the fuselage.  Suggested usage:
C
C                   NFIXL, -R, -A, -B = 2, 4, 10, 10 (HSCT wing)
C                   NFIXL, -R, -A, -B = 4, 12, 4, 4  (bizjet wing or aft pylon)
C
C        NCRANK     # of cranks specified for this surface (LE and/or TE)
C        LCRANK     Crank # at which to stop rectifying paneled leading edge;
C                   0 <= LCRANK <= NCRANK
C        ZCRANK     0 or more spanwise (Z) stations at which cranks are to be
C                   contained precisely in the paneled surface;
C                   enter Y stations for vertical fins
C
C     Streamwise and Spanwise Paneling Distribution Controls:
C
C        CRDSPL     Streamwise linear coefficient, e.g. 0.04  |
C        CRDSPQ     .......... quadratic .............. 0.0   | Sum is 1.0
C        CRDSPS     .......... sine ................... 0.3   |
C        CRDSPC     .......... cosine ................. 0.66  |
C        SPNSPC     Spanwise cosine coefficient,   e.g. 0.35
C        SPNSPS     ........ sine ..................... 0.0  (0., 0. = uniform)
C
C     Fuselage Controls:
C
C        PIECE      Descriptor (up to 32 characters)
C        JLBODY     Number of paneled points to impose circumferentially
C        IDEGFUS    1 or 3 means piecewise [bi]linear or [bi]cubic interpolation
C        NNFUDGE    IW > 0 allows fudging of fuselage sections if wing IW drops
C                          below it;
C                   0 suppresses the fudging, which applies to just one wing
C        IBOD1,     Fuselage section #s to which possible fudging is confined
C        IBOD2
C        IBOD3      IBOD3 = 0     suppresses the body shelf after the wing UPPER
C                                 surface is found below the body in WBFUDGE;
C                   IBOD3 > IBOD2 forces the shelf to continue all the way to
C                                 IBOD2 regardless; it is then decayed smoothly
C                                 to zero at station IBOD3 unless ...
C                   IBOD3 = 999   carries the shelf to the end of body without
C                                 decay
C        NSHELF     Number of data points defining YSHELF & ZSHELF (next)
C
C        XSHELF     Abscissas in [0,1] transformed internally to the X range
C                   of fuselage stations IBOD1 and IBOD2 then used with Y/ZSHELF
C        YSHELF     Data points defining a safety margin added to fudged body
C                   section Ys ensure a wing/body intersection (units of body Y)
C        ZSHELF     Data points defining the fraction of local fuselage width
C                   used to force a minimum "shelf" width built onto body
C                   sections, providing room for paneling between the keel and
C                   the wing/body intersection (with monitoring of the upper
C                   wing surface to prevent the shelf from being unnecessarily
C                   wide); see also IBOD3
C
C     Shape Functions:
C
C        NDV        Number of design variables, >= 0
C        WINGOUT    Indicates desired form of saved wing geometry if NDV > 0:
C                      0 = sections are untwisted and normalized;
C                      1 =    "      "    twisted and normalized;
C                      2 =    "      "    twisted and denormalized
C
C        VTYPE      6-character shape function typically applied to 2-D sections
C        UTYPE      4-character shape function controlling spanwise lofting
C
C                   Basic Hicks-Henne-type perturbing function descriptions:
C
C          SCALE    Simple Y scaling;
C          SIN      Standard (modified) sine perturbing function
C          SINC     Cyclic sine function suited to lofting heat shield
C                   perturbations in the azimuthal direction
C          SINF     Flipped sine function
C          SIN1     Symmetric sine function
C          SIN2     Symmetric flipped sine function
C          SIN3     Ensures airfoil LE/TE symmetry; use
C                   [XA, XB] = [0, 1] & XC in (0, .5)
C          SIN4     Ensures airfoil LE/TE symmetry; use
C                   [XA, XB] = [0, .5] & XC in (0, 1)
C          COSL     Half [co]sine, peak in [XA, XC], 0 at XB
C          COSR     Half [co]sine, peak in [XC, XB], 0 at XA
C          LCOS     Inverted form of COSL, peak in [XA, XC], 0 at XB
C          RCOS     Inverted form of COSR, peak in [XC, XB], 0 at XA
C          EXP      Exponential function (variable center)
C          DROOP    Simplified EXP with peak at X/C = 0;
C                   affects both surfaces the same way
C          LEAD     Leading  edge droop function (also LED; see also DROOP);
C                   affects both surfaces the same way
C          TRAIL    Trailing edge droop function (also TRL);
C                   affects both surfaces the same way
C          WAG      Wagner shape function; N is indicated via "POWER"
C
C        See SYN107-MB for more complete shape function descriptions.
C        Here's a new one:
C
C        'DIHED '   Applies dihedral to indicated wing sections:
C                   KMIN & KMAX = UBMIN & UBMAX = 1st & last sections rotated;
C                   KCEN = UBUMP = section whose real-space LE & TE (presumably
C                          with no flap or slat deflection) define the axis of
C                          rotation in a constant-Z plane; KCEN = KMIN usually;
C                   V * VSCALE = rotation angle about the axis (in degrees);
C                              < 0 = anhedral
C
C
C     AEROSURF Geometry Inputs (aero.geo file)
C     ----------------------------------------
C
C     Wing-type Sections:
C
C        Wing/canard/horizontal tail surfaces must be in forward-to-aft order.
C        They should be followed by any vertical tail/under-wing pylon/diverter
C        surfaces.  Any nacelle surfaces appear last.
C
C        TITLE  Description of wing    [or nacelle]
C
C        FSEC                          [FSEC, FQUAD, FINTR, CLOCK]
C
C        # sections                    [FQUAD & FINTR refer to the quadrant # of
C        (ALPHAW is no longer used)     the first nacelle section and that of
C                                       the intersection; enter 3, 6, 9 or 12
C                                       o'clock as viewed from upstream;
C                                       further sections are specified clockwise
C                                       at fractions of a turn from 12 o'clock;
C                                       CLOCK is added to the azimuthal angles
C                                       input with nacelle sections (also input
C                                       as a fraction from 12 o'clock); CLOCK
C                                       represents nacelle roll about its axis.]
C
C                                      [XC, YC, ZC (inlet), XC, YC, ZC (outlet)]
C
C                                      [XT1, YT1, ZT1, XT2, YT2, ZT2
C                                       thrust center and pt. forming a normal
C                                       vector at the outlet plane]
C
C        ZLE, XLE, YLE, CHD, THK, TWIST, FSEC, ROUND [, ROLL] (if IROLL(IW) = 1)
C       [XLE, RLE, TLE, .... for nacelles ..............; TLE=frctn. of 360 deg]
C
C                                       ROLL>0 is clockwise from upstream (deg.)
C
C        YSYM, FNU, FNL                [Symmetry option is not supported.]
C        0. & # pts. upper & lower
C
C       [DSLAT, XSLAT, SWSLAT,         [If IFLAP(IW) = 1:
C           DFLAP, XFLAP, SWFLAP]       DSLAT = slat deflection about X = XSLAT,
C                                       DFLAP = flap    "   "   "   " X = XFLAP,
C                                       where +ve in degrees is down for both;
C                                       XSLAT & XFLAP are in real space;
C                                       DSLAT & DFLAP are rotations about the
C                                       hinge line; SWSLAT & SWFLAP are the
C                                       sweeps of the hinge lines, needed to
C                                       determine each point rotation in the
C                                       plane of the defining section]
C
C       [DIHED, DIAXIS(1:3,1:2)]       [If IDIHED(IW) = 1:
C                                       DIHED = dihedral angle, degrees;
C                                       DIAXIS(*,1) = upstrm. (x,y,z) of axis
C                                       DIAXIS(*,2) = downst. (x,y,z) of axis]
C
C        X, Y                           (x,y) coordinates for upper surface,
C        X, Y                           followed by lower surface
C        :  :
C                                      [or (r,theta) coordinates for nacelles]
C
C     Fuselage Sections (following all wing & nacelle surfaces):
C
C        FNFSEC, FLIP ! Note the lack of title, for compatibility with SYN87-SB
C
C        FN, XF, FSEC  (FLIP = 0.)  or  XF, FN, FSEC  (FLIP = 1.)
C
C        Z, Y                       or  Y, Z
C        Z, Y                           Y, Z
C        :  :                           :  :
C
C     The program ensures body section points go from keel to crown.
C
C
C     Other quantities are documented below and in the SURFER header.
C
C
C     Environment:  Fortran 90
C     ------------
C
C     Funding:  Aerodynamics Division, NASA Ames Research Center, CA
C     --------
C
C     History:
C     --------
C
C     Sep. 1997  James Reuther   Initial implementation, adding nacelles and
C                                pylons to WBSURF capabilities.
C     Oct. 1997  David Saunders  General overhaul; incorporated missing spanwise
C                                gridding of wing surfaces and the wing-below-
C                                the-body protection scheme from SYN87-SB.
C     11/22/97        "          All intersection calculations now use WBINTR.
C     11/26/97        "          Introduced XTRAP to protect aft pylon intersec.
C     11/28/97     DAS/JJR       Introduced IROLL/ROLL and IQUAD controls.
C     12/01/97       DAS         Coded WDINTR, extending lower wing surface if
C                                XTRAP>0. to protect each diverter intersection.
C     12/27/97        "          Introduced KSTN1/2 during WDSURF development.
C     01/02/98        "          Underwing nacelles need (pylon) roll angles to
C                                be clockwise (from upstream) and input pylon
C                                section Ys such that "upper" <-> outboard.
C     01/29/98        "          KSTN1/2 are now determined by SURFER.
C     02/02/98        "          Added NFIXL, -R, -A, -B inputs for FSURF's
C                                redistributions via FIXROOT.
C     02/26/98        "          Clock angles added to nacelle geometry inputs.
C     03/04/98        "          Pylon/diverter centerlines need not be strictly
C                                streamwise; corresponding crank station inputs
C                                need be nominal only.
C     03/18/98        "          ZINNER and ZOUTER inputs needed with the new
C                                XTRAP < 0. option for closing diverters at the
C                                wing trailing edge and paneling aft portions of
C                                HSCT-type nacelles.
C     04/28/98        "          NNFUDGE = IW overcomes original assumption of 1
C                                which broke for canard cases.
C     08/17/98        "          NNFUSE = 2 kludge allows saving of aero.body
C                                file (sections needed for original AEROSURF).
C     08/28/98        "          New NNSPLIT input and subroutine FREPANEL give
C                                an alternative to the original FSPLIT for cases
C                                where nearly-aligned wing leading edges could
C                                lead to variable numbers of fuselage patches.
C     09/23/98        "          James found he needs a structured wing surface
C                                in PERTURB to control diverter height.  Thus
C                                regularizing the wing(s) (and the body), though
C                                rightfully part of SURFER's paneling, is now
C                                moved up to the PERTURB level.  X/Y/ZBASE is
C                                used now to pass SURFER the wing sections, and
C                                fuselage control inputs are moved to be last.
C     10/16/98        "          Introduced some dynamic allocation in AEROSURF
C                                to match the design code's use of SURFER more
C                                closely.  In particular, the regularized fuse-
C                                lage J dimension is no longer the same as that
C                                of the patch dimensions.
C     10/22/98     DAS/JJR       The nacelle transformations are now done in
C                                PERTURB for consistency; nacelle geometry may
C                                need fudging (WNFUDGE) near the wing trailing
C                                edge to protect the scheme of closing diverters
C                                there at ZINNER/ZOUTER.
C     10/26/98       DAS         NNFORMAT = 2 or 3 adds saving of regularized
C                                surfaces input to SURFER, as file aero.regular.
C     12/04/98        "          WBFUDGE shouldn't capture the edge of the shelf
C                                if the previous J point was not fudged.
C     12/09/98        "          NNFUSE = 3 allows access to the unregularized
C                                body sections as fuselage.plot3d, assuming the
C                                pt. counts allow that; changed aero.regular
C                                name to aero.components (see NNFORMAT usage).
C     12/12/98        "          NNINSERT permits inserting wing sections ...
C     12/15/98        "          ... and inserting nacelle & fuselage sections;
C                                tagged patches with component # at the end of
C                                aero.xyz, plus descriptive name in aero.out.
C     12/18/98        "          NNKAPPA permits evaluation of curvature at each
C                                point of the indicated regularized components,
C                                in files aero.ikappa & .jkappa, then stopping.
C     12/21/98        "          Prompt for curvature or radius of curvature.
C     01/20/99        "          Combined curvature outputs in index directions
C                                into one file for FAST, and wrote corresponding
C                                (x,y,z)s to aero.regular instead of relying on
C                                aero.components.
C     01/21/99        "          Reverted to aero.ikappa & .jkappa because FAST
C                                treats more than 1 function as vector data.
C     01/28/99        "          Dynamic allocation of (unpacked) geometry and
C                                surface patch arrays.
C     01/30/99        "          Packed the regularized wing components.
C     02/11/99        "          Packed the surface patches.
C     03/15/99        "          Introduced IFLAP(*) and DFLAP/XFLAP & D/XSLAT.
C     03/19/99        "          SWFLAP/SWSLAT were also needed for proper
C                                recovery of deflections stored in aero.geo.
C     03/29/99        "          Overcame the wing root LE patch splitting
C                                problem by rectifying the plain paneled wing
C                                for K = 1:MCRANK(LCRANK(IW),IW).
C     04/23/99        "          Subpatch version.
C     06/02/99        "          Introduced NNTAIL to deal with vertical fins.
C     08/18/99        "          Retrofitted all the shape functions in PERTURB.
C     10/13/99        "          NNPANEL input permits AEROSURF to force current
C                                paneling to match that of aero.xy0 size-wise.
C     01/04/00        "          NNEXPAND input controls use of EXPAND_PATCHES
C                                as needed by the multiblock warping scheme.
C     06/23/00        "          DIHED shape function simplifies flap/slat
C                                deflections for surfaces with large dihedrals.
C     07/19/00        "          NNCONFIG introduced to handle different body
C                                paneling for CTV-class vehicles with over-
C                                lapping wing & tail surfaces.
C     07/27/00        "          NNCONFIG = 2 invokes FCTV2, not FCTV in SURFER.
C     08/07/00        "          IBOD3 > IBOD2 overrides possible suppression of
C                                the shelf in WBFUDGE.
C     08/08/00        "          NNCONFIG = 3 invokes FCTV3 (center-line fin).
C     08/09/00        "          Generalized the controls for WBFUDGE from the
C                                same YSHELF & ZSHELF everywhere to variation
C                                with respect to X via splining data entered
C                                in aero.inp.
C     09/12/00        "          Added surface area calculations for possible
C                                use with viscous drag estimates.
C     10/03/00        "          REORDER needed protecting in FLAPSLAT: ensure
C                                that the true leading edge index is used.
C     04/11/01        "          Fuselage base area may be (almost) filled with
C                                additional defining sections (see NNFUSE = 2).
C     04/13/01        "          Better idea: an extra patch is constructed for
C                                an open fuselage base.
C     04/18/01        "          WINTERP shouldn't zero out the twists of the
C                                original sections.
C     07/16/03        "          If NNPANEL = 2 (for forcing patch dimensions
C                                to match those of aero.xy0), the check for the
C                                anticipated number of patches to match the
C                                number found in aero.xy0 was causing premature
C                                termination.  Some cases need an extra patch
C                                for temporary work-space.  Therefore, use the
C                                larger number of patches when allocating the
C                                arrays for the aero.xy0 patch dimensions.
C     08/20/03        "          Trapped more geometry reading errors.
C     10/17/03        "          Spanwise sinusoidal lofting wasn't clipping at
C                                ZMIN and ZMAX.  Trapping of ZCEN < ZMIN can
C                                save some head-scratching.
C     11/03/03        "          Generalized COS* and *COS shape functions to
C                                ease some kinds of spanwise lofting.
C     11/06/03        "          LOFT needed a work-around for the generalized
C                                COS* and *COS options.  SFEVAL's description of
C                                the "center" argument is affected.
C
C*******************************************************************************

C     Global variables:

      USE GEOM0
      USE GEOM1
      USE GRIDS
      USE SHAPE_FUNCTIONS

      IMPLICIT NONE

C     Local constants:

      INTEGER, PARAMETER ::
     >   IREAD   = 1, ! Control inputs
     >   IWRIT   = 2, ! Printable outputs
     >   IKBD    = 5, ! Keyboard inputs (NNINSERT option)
     >   ICRT    = 6, ! Screen outputs
     >   IGEOM   = 7, ! Geometry inputs (& possibly outputs)
     >   ICOMP   = 8, ! Regularized components output from PERTURB
     >   ISAVE   = 9, ! Patched surface outputs (& curvature outputs)
     >   MXSHELF = 10 ! Limit on NSHELF

      REAL, PARAMETER ::
     >   HALF = 0.5, ONE = 1., ZERO = 0.

C     Local variables:

C     NOTE: MXJGRID became redundant once patches were packed.  MXKGRID
C           serves for wings and nacelles; JLBODY serves for the fuselage.

      INTEGER
     >   I, I1, IBOD1,  IBOD2,  IBOD3, IDEGFUS, ILE, IOS, IW, J, JLBODY,
     >   K, L, LASTXB, LENRG, LENXB,
     >   M, MOSTXB, MXIGRID, MXKGRID, MXINTR, MXPATCH, MXSUBI, MXSUBJ,
     >   N, NBASE, NCALL, NEWPTS, NFSTN, NI, NIWING, NJ, NJFUSE, NKWING,
     >   NL, NNCONFIG, NNEXPAND, NNFORMAT, NNFUDGE, NNFUSE, NNINSERT,
     >   NNKAPPA, NNNAC, NNPANEL, NNSPLIT, NNREPEAT, NNSURF, NNWING,
     >   NPATCH, NSHELF, NU, NWCRANK, NWSURF, WINGOUT

      REAL
     >   AREA_PATCH, AREA_TOTAL, CHORD, DX, EFFTWT, FN,
     >   RAD2DEG, SCALE, TCMAX, TEANGL, XFI, XLAST, XLE, XLENGTH, XMIN,
     >   XTCMAX, YFI, YLE, ZFI, ZLE, XTRAIL, YTRAIL, YCRN, YKEL, YM,
     >   ZWAT

      REAL(KIND=4) :: CPU1, CPU2

      REAL, DIMENSION (MXSHELF) ::
     >   XSHELF, YSHELF, ZSHELF

      LOGICAL
     >   FAIL, PERTURBED

      CHARACTER
     >   TITLE * 132

C     Miscellaneous work-space allocated after geometry dimensions are known.
C     Declare them here, not in module GEOM1, because GEOM1 parallels GEOM0.

      INTEGER, ALLOCATABLE, DIMENSION (:) ::
     >   IDIM, IRG, KDIM, NUMK


C     Execution:
C     ----------

      RAD2DEG = 45./ ATAN (ONE)

C     Open input and diagnostic output files:

      OPEN (UNIT=IREAD, FILE='aero.inp', STATUS='OLD')
      OPEN (UNIT=IWRIT, FILE='aero.out', STATUS='UNKNOWN')
      OPEN (UNIT=IGEOM, FILE='aero.geo', STATUS='OLD')

      WRITE (IWRIT, '(A)') ' AEROSURF version of November 6, 2003'

C     Echo the inputs to the output file:

      WRITE (IWRIT, '(/, A, /)') ' Control inputs:'

      DO
         READ  (IREAD, '(A)', IOSTAT=IOS) TITLE
         IF (IOS < 0) EXIT

         L = LEN_TRIM (TITLE)
         WRITE (IWRIT, '(1X, A)') TITLE(1:L)
      END DO

      WRITE (IWRIT, '(A)')
      REWIND (IREAD)


C     Allocate temporary geometry arrays, to avoid reading the geometry
C     more than once:

      CALL DYNAMIC_GEOM0 (1)


C     Read the control inputs other than the optional design variables:

      CALL CONTROL (IREAD,   IWRIT,   MXWSURF, MXKWING, MXCRANK,
     >              LENXB,   PIECE0,
     >              NNCONFIG,NNWING,  NNNAC,   NNFUSE,  NNSPLIT,
     >              NNFORMAT,NNFUDGE, NNREPEAT,NNINSERT,NNKAPPA,
     >              NNPANEL, NNEXPAND,NNROOT0, NNTAIL0,
     >              ILWING0, JLBODY,  KLWING0, KWING10, KWING20,
     >              IROLL0,  IFLAP0,  IDIHED0, IDEGREE0,IDEGFUS,
     >              IBOD1,   IBOD2,   IBOD3,
     >              MXSHELF, NSHELF,  XSHELF,  YSHELF,  ZSHELF,
     >              XTRAP0,  ZINNER0, ZOUTER0,
     >              CRDSPL0, CRDSPQ0, CRDSPS0, CRDSPC0,
     >              SPNSPC0, SPNSPS0, NFIXL0,  NFIXR0,
     >              NFIXA0,  NFIXB0,  NCRANK0, LCRANK0, ZCRANK0)


      NNSURF = NNWING + NNNAC
      FAIL = .FALSE.

C     Some cases have restrictions:

      I = NNCONFIG
      IF (I == 3) I = 2 ! FCTV3 has same requirements as FCTV2

      SELECT CASE (I)

         CASE (0) ! All cases prior to CTV

         CASE (1) ! FCTV body paneling

            IF (ILWING0(1) /= ILWING0(2)) THEN
               WRITE (ICRT, '(/, A)')
     >            ' Wing & tail ILWINGs must match in FCTV routine.'
               FAIL = .TRUE.
            END IF

            IF (MOD (JLBODY, 2) /= 1) THEN
               WRITE (ICRT, '(/, A)')
     >            ' JLBODY must be odd for CTV cases.'
               FAIL = .TRUE.
            END IF

         CASE (2) ! FCTV2 or -3 body paneling

            IF (ILWING0(1) < ILWING0(2) + 4) THEN
               WRITE (ICRT, '(/, A)')
     >            ' Wing ILWING >= tail ILWING + 4 for FCTV2/3 method.'
               FAIL = .TRUE.
            END IF

            IF (MOD (JLBODY, 2) /= 1) THEN
               WRITE (ICRT, '(/, A)')
     >            ' JLBODY must be odd for CTV cases.'
               FAIL = .TRUE.
            END IF

      END SELECT

      IF (FAIL) GO TO 999


      IF (NNINSERT > 0 .AND. NNFUSE <= 1) THEN
         NNFUSE = 0
         NNSURF = NNINSERT ! No need to read beyond the target wing component
         IF (NNINSERT <= NNWING) THEN
            NNWING = NNINSERT
            NNNAC = 0
         ELSE
            NNNAC = NNINSERT - NNWING
         END IF
      END IF

      IF (NNKAPPA > 0) THEN ! Likewise ...
         NNFUSE = 0
         NNSURF = NNKAPPA
         IF (NNKAPPA <= NNWING) THEN
            NNWING = NNKAPPA
            NNNAC = 0
         ELSE
            NNNAC = NNKAPPA - NNWING
         END IF
      END IF

      NWSURF = MAX (NNSURF, 1) ! Avoid passing 0 as a dimension


C     Read wing, pylon, and nacelle geometries:

      DO IW = 1, NNSURF

         CALL WINGREAD (IGEOM, IW, NNWING, MXWSURF, MXIWING, MXKWING,
     >                  IROLL0, IFLAP0, IDIHED0,
     >                  NWSEC0, NLG0, NUG0, NWG0, ALG0,
     >                  XLEADI0, YLEADI0, ZLEADI0, ROUNDED0, ROLL0,
     >                  DFLAP0, DSLAT0, XFLAP0, XSLAT0, SWFLAP0,SWSLAT0,
     >                  AINIT0, CINIT0, TINIT0, XINIT0, YINIT0, ZINIT0,
     >                  XC10, YC10, ZC10, XC20, YC20, ZC20,
     >                  DIHED0, DIAXIS0, IQUAD0, CLOCK0, TITLE, IWRIT)
      END DO


C     Read fuselage geometry:

      IF (NNFUSE /= 0) THEN

         CALL FUSEREAD (IGEOM, MXIBODY, MXJBODY, NFSTN, NFJ0,
     >                  XF0, YFINIT0, ZFINIT0, IWRIT)
      END IF

      CLOSE (IGEOM)


C     Determine the largest number of sections across components
C     and the largest number of points across sections for wings.
C     Determine the size of the packed, regularized wing sections.

      NIWING  = 1 ! Avoid 0 dimensions in case NNSURF = 0
      NKWING  = 1
      NWCRANK = 1
      MXIGRID = 1
      MXKGRID = 1

      DO IW = 1, NNSURF

         DO K = 1, NWSEC0(IW)
            NIWING  = MAX (NIWING, NWG0(K,IW))
         END DO

         NKWING  = MAX (NKWING,   NWSEC0(IW))
         NWCRANK = MAX (NWCRANK, NCRANK0(IW))
         MXIGRID = MAX (MXIGRID, ILWING0(IW))
         MXKGRID = MAX (MXKGRID, KLWING0(IW))

      END DO

      NWCRANK = NWCRANK + 1 ! Allow for adding the tip as a crank
      NJFUSE  = 1
      MXINTR  = MXIGRID ! Unless there's a fuselage + wing-mounted pylons

      IF (NNFUSE /= 0) THEN

         DO I = 1, NNWING ! For wing-mounted pylon/diverter cases
            IF (NNROOT0(I) > 0) MXINTR = (3 * MXIGRID) / 2
         END DO

         DO I = 1, NFSTN
            NJFUSE = MAX (NJFUSE, NFJ0(I))
         END DO
      ELSE
         NFSTN  = 1 ! Avoid passing 0 as a dimension
         JLBODY = 1
      END IF


C     Allocate more precise geometry-related arrays:

      CALL DYNAMIC_GEOM1 (NWSURF, NIWING, NKWING, NWCRANK,
     >                    NFSTN,  NJFUSE)


C     Transfer the geometry-related arrays ...

      CALL COPY_GEOM (NWSURF, NNFUSE, NIWING, NKWING, NWCRANK,
     >                NFSTN,  NJFUSE)

C     ... and release the interim arrays:

      CALL DYNAMIC_GEOM0 (2)


C     Look for optional shape functions/design variables:

      DO I = 1, 3
         READ (IREAD, *)
      END DO

      READ (IREAD, *) NDV, WINGOUT

      CALL RDVARS (IREAD, IWRIT)

      CLOSE (IREAD)


C     Set up helpful regularized geometry dimensions and pointers:

      ALLOCATE (IDIM(NWSURF), KDIM(NWSURF), NUMK(NWSURF), IRG(NWSURF))

C     IDIM(IW) = Allocated array I dimen. for regularized wing-type surface IW
C     KDIM(IW) = ............... K ...........................................
C     NUMK(IW) = Active number of regularized sections ...
C     IRG(IW)  = First word of regularized sections X/Y/ZRG(*) for surface IW

      DO IW = 1, NNWING
         IDIM(IW) = ILWING(IW) + 1         ! Allow for extrapolating at TE.
         NUMK(IW) = NWSEC(IW)              ! SURFER inserts crank stations,
         KDIM(IW) = NWSEC(IW) + NCRANK(IW) ! not PERTURB
      END DO

      DO IW = NNWING + 1, NNSURF ! SURFER inserts a pt. at pylon intersection
         IDIM(IW) = ILWING(IW) + 1         ! Allow for extrapolating at TE
         NUMK(IW) = KLWING(IW) - 1         ! NCYLINDER overwrites its inputs
         KDIM(IW) = MAX (NWSEC(IW), NUMK(IW)) + NCRANK(IW)
      END DO


C     Determine the total size of the packed regularized wing geometry arrays,
C     & the pointers IRG(*) to the start of each wing surface in X/Y/ZRG(*):

      IRG(1) = 1
      LENRG  = IDIM(1) * KDIM(1)

      DO IW = 2, NNSURF
         IRG(IW) = LENRG + 1
         LENRG   = IDIM(IW) * KDIM(IW) + LENRG
      END DO


C     Determine the number of paneled surface patches to allocate:

      IF (NNFUSE >= 2 .OR. NNINSERT /= 0 .OR. NNKAPPA /= 0) THEN

         NPATCH  = NNSURF + MIN (1, NNFUSE) ! No SURFER; I/JDIMENS need + NNFUSE
         MXSUBI  = 1
         MXSUBJ  = 1
         NNPANEL = 0 ! Suppress possible processing of aero.xy0

      ELSE

         CALL COUNT_PATCHES (NNCONFIG, NNWING,  NNNAC,  NNFUSE,
     >                       NNSPLIT,  NNROOT,  NNTAIL, NNEXPAND,
     >                       NWSURF,   NWCRANK, ILWING, XTRAP,
     >                       MXSUBI,   MXSUBJ,  NPATCH, IWRIT)
      END IF

      MXPATCH = NPATCH
      WRITE (ICRT, '(/, A, I3)') ' Patches allocated:', MXPATCH

C     Allocate and read size arrays from aero.xy0 if NNPANEL = 2.
C     Even if NNPANEL = 1, SURFER needs I/JMXCOMO(*) to report possibly
C     changing patch sizes.

      CALL READ_SIZES (NNPANEL, MXPATCH, MXSUBI, MXSUBJ,
     >                 ISAVE, ICRT, NNFORMAT)


C     Allocate the normal surface-patch-related arrays:

      CALL DYNAMIC_GRIDS (NWSURF, MXPATCH, NFSTN, JLBODY, LENRG,
     >                    NWCRANK, MXIGRID, MXKGRID, MXINTR,
     >                    MXSUBI, MXSUBJ, LENXB)


C     Proceed with work arrays as precisely sized as possible.
C     The MODULEs aren't used beyond here - everything is argument-driven.
C
C     Correction:  the shape functions forced use of a new module by PERTURB,
C     but this is the only hidden communication below here.  SURFER remains
C     fully argument-driven.


      NCALL = -3 ! First call flag

  500 CONTINUE ! Come back here for the repeat check


C     Apply any shape functions, denormalize, possibly adjust fuselage sections
C     to protect intersection calculations, and regularize the sections:

      CALL PERTURB (NCALL, NWSURF, NIWING, NKWING, NFSTN, NJFUSE,
     >              MXIGRID, NNCONFIG, NNWING, NNNAC, NNFUSE, LENRG,
     >              NWSEC, IDIM, KDIM, IRG, CLOCK, ALG, ROLL, ROUNDED,
     >              XLEADI, YLEADI, ZLEADI, DIHED, DIAXIS,
     >              DFLAP, DSLAT, XFLAP, XSLAT, SWFLAP, SWSLAT,
     >              AINIT, CINIT, TINIT, XINIT, YINIT, ZINIT,
     >              NLGI, NUGI, NWGI, NLG, NUG, NWG,
     >              ILWING, KLWING, XWG, YWG, ZWG,
     >              XC1, YC1, ZC1, XC2, YC2, ZC2,
     >              SNORM, CRDSPL, CRDSPQ, CRDSPS, CRDSPC,
     >              XRG, YRG, ZRG, NFSTN, NFJ, JLBODY, IDEGFUS,
     >              XF, YFINIT, ZFINIT, YF, ZF, XTRAP, ZINNER, ZOUTER,
     >              NNROOT, NNFUDGE, NNREPEAT,
     >              IBOD1, IBOD2, IBOD3, KWING2,
     >              NSHELF, XSHELF, YSHELF, ZSHELF, IWRIT)


C     Save perturbed geometry sections?

      IF (NCALL == -3 .AND. NDV > 0) THEN

         OPEN (UNIT=IGEOM, FILE='aero.geo2', STATUS='NEW')

         DO IW = 1, NNSURF

            PERTURBED = .FALSE.

            DO N = 1, NDV
               IF (IDWING(N) == IW) PERTURBED = .TRUE.
            END DO

            IF (PERTURBED) CALL WINGWRITE (IGEOM, IW, WINGOUT)

         END DO

         IF (NNFUSE /= 0) THEN

            IF (NTOTFU + NTKINK > 0) CALL FUSEWRITE (IGEOM, NFSTN)

         END IF

         CLOSE (IGEOM)

      END IF


C     Regularized sections are returned from PERTURB in packed form.

      IF (NCALL == -3 .AND. NNFORMAT >= 2) THEN

C        Option to save the regularized surfaces that are input to SURFER.
C        Set up the active dimensions:

         DO N = 1, NNSURF
            IDIMENS(N) = ILWING(N)
            JDIMENS(N) = NUMK(N)
         END DO

         N = NNSURF + MIN (1, NNFUSE)

         IF (NNFUSE > 0) THEN
            IF (NNCONFIG == 0) THEN
               M = JLBODY - 1 ! Allows for inserting intersection pts. in FSURF
            ELSE
               M = JLBODY     ! FCTV[2|3] cases
            END IF
            IDIMENS(N) = M
            JDIMENS(N) = NFSTN
         END IF

         IF (NNFORMAT == 2) THEN ! Unformatted

            OPEN (UNIT=ICOMP, FILE='aero.components', STATUS='UNKNOWN',
     >            FORM='UNFORMATTED')

            WRITE (ICOMP) N
            WRITE (ICOMP) (IDIMENS(I), JDIMENS(I), 1, I = 1, N)

         ELSE ! Formatted

            OPEN (UNIT=ICOMP, FILE='aero.components', STATUS='UNKNOWN')

            WRITE (ICOMP, '(I4)') N
            WRITE (ICOMP, '(8 (2I4, I2))')
     >         (IDIMENS(I), JDIMENS(I), 1, I = 1, N)

         END IF

         DO N = 1, NNSURF

            I1 = IRG(N)

            CALL XYZ_WRITE (NNFORMAT, ICOMP, IDIM(N), KDIM(N),
     >                      IDIMENS(N), JDIMENS(N),
     >                      XRG(I1), YRG(I1), ZRG(I1))
         END DO

         IF (NNFUSE > 0) THEN

            IF (NNFORMAT == 2) THEN

               WRITE (ICOMP)
     >            ((XF(I),   J = 1, M), I = 1, NFSTN),
     >            ((YF(J,I), J = 1, M), I = 1, NFSTN),
     >            ((ZF(J,I), J = 1, M), I = 1, NFSTN)

            ELSE ! NNFORMAT == 3

               WRITE (ICOMP, '(6F12.6)')
     >            ((XF(I),   J = 1, M), I = 1, NFSTN),
     >            ((YF(J,I), J = 1, M), I = 1, NFSTN),
     >            ((ZF(J,I), J = 1, M), I = 1, NFSTN)

            END IF

         END IF

         CLOSE (ICOMP)

      END IF


C     Option to evaluate I & J [radius of] curvature on a regularized component:

      IF (NNKAPPA /= 0) THEN

         IW = MAX (NNKAPPA, 1) ! May be -1 for fuselage
         I1 = IRG(IW)

         CALL KAPPA (NNKAPPA, NNFORMAT, ICRT, IKBD, ISAVE,
     >               IDIM(IW), KDIM(IW), ILWING(IW), NUMK(IW),
     >               NFSTN, JLBODY, XRG(I1), YRG(I1), ZRG(I1),
     >               XF, YF, ZF)

         GO TO 999

      END IF


C     Option to insert defining sections via interpolation.
C     The raw geometry is required to be regular.

      IF (NNINSERT >= 1 .AND. NNSURF >= 1 .AND. NNINSERT <= NNSURF) THEN

         CALL WINTERP (NNINSERT, NNWING, NNNAC,
     >                 NWSURF, NIWING, NKWING, NWSEC, IDEGREE,
     >                 CLOCK, ROLL, XLEADI, YLEADI, ZLEADI,
     >                 XINIT, YINIT, ZINIT, AINIT, CINIT, TINIT, ALG,
     >                 NLG, NUG, NWG, XWG, YWG, ZWG, ROUNDED,
     >                 IKBD, ICRT, ISAVE)
         GO TO 999

      ELSE IF (NNINSERT == -1 .AND. NNFUSE /= 0) THEN

         CALL FINTERP (NFSTN, NJFUSE, NFSTN, NFJ, IDEGFUS,
     >                 XF, YFINIT, ZFINIT, IKBD, ICRT, ISAVE)
         GO TO 999

      END IF


C     Summarize the denormalized defining sections, and save other results:

      IF (NCALL == -3) THEN

         DO IW = 1, NNSURF

            L = LEN_TRIM (PIECE(IW))

            WRITE (IWRIT, '(//, A, I2, A, 2X, A, //, A, A, A, /)')
     >         ' Defining section summary for "wing" surface', IW, ':',
     >         PIECE(IW)(1:L),
     >         ' Section    Xle        Yle        Zle      ',
     >         'Chord        Xte        Yte    Thick  T/C Max',
     >         '   at X/C Lo.Twist To.Twist Ef.Twist TE.Ang.'

            DO K = 1, NWSEC(IW)

               N   = NWG(K,IW)
               NU  = NUG(K,IW)
               NL  = NLG(K,IW)
               ILE = NL
               XLE = XWG(ILE,K,IW)

               DO I = NL - NL/4, NL + NL/4
                  IF (XWG(I,K,IW) < XLE) THEN
                     XLE = XWG(I,K,IW)
                     ILE = I
                  END IF
               END DO
               YLE = YWG(ILE,K,IW)
               ZLE = ZWG(ILE,K,IW)

C              Thick sections with large flap deflections (for instance) can
C              lead to nonmonotonic Xs, so determine T/Cmax from the normalized
C              untwisted sections, not the denormalized sections:

               CALL GETTHICK (NL, NU, XINIT(1,K,IW), YINIT(1,K,IW),
     >                        TCMAX, XTCMAX)

               XTRAIL = (XWG(1,K,IW) + XWG(N,K,IW)) * HALF
               YTRAIL = (YWG(1,K,IW) + YWG(N,K,IW)) * HALF
               EFFTWT = ATAN ((YLE - YTRAIL) / (XTRAIL - XLE)) * RAD2DEG

               CALL GETANGLE (NL, NU, XWG(1,K,IW), YWG(1,K,IW), ONE,
     >                        TEANGL)
               TEANGL = TEANGL * RAD2DEG

               WRITE (IWRIT, '(I4,6F11.4,F9.4,F9.6,4F9.4,F8.3)') K, XLE,
     >            YLE, ZLE, CINIT(K,IW), XTRAIL, YTRAIL, TINIT(K,IW),
     >            TCMAX, XTCMAX, AINIT(K,IW), ALG(K,IW), EFFTWT, TEANGL

             END DO

         END DO

         IF (NNFUSE /= 0) THEN

            WRITE (IWRIT, '(//, A, //, A, /)')
     >         ' Input fuselage station summary:',
     >         ' Station        X     Y Crown      Y Keel       Z Max'

            DO I = 1, NFSTN

               YCRN = YFINIT(1,I)
               YKEL = YCRN
               ZWAT = ZFINIT(1,I)

               DO J = 1, NFJ(I)
                  YCRN = MAX (YFINIT(J,I), YCRN)
                  YKEL = MIN (YFINIT(J,I), YKEL)
                  ZWAT = MAX (ZFINIT(J,I), ZWAT)
               END DO

               WRITE (IWRIT, '(I4, 1X, 4F12.5)')
     >            I, XF(I), YCRN, YKEL, ZWAT

            END DO

         END IF


         IF (NNFUSE == 2) THEN ! Kludge to save possibly-fudged fuselage
                               ! and/or close an open base

            I  = NFSTN
            J  = NFJ(I)
            YM = (YFINIT(1,I) + YFINIT(J,I)) * HALF

            IF (MAX (ZFINIT(J/2,I), ABS (YM - YFINIT(1,I))) >
     >          (XF(I) - XF(1)) * 1.E-5) THEN

               WRITE (ICRT, '(/, A, /, A, 2F10.4)')
     >            ' Option to close fuselage base:',
     >            ' Nose and current aft body Xs: ', XF(1), XF(I)
               WRITE (ICRT, '(2A)', ADVANCE='NO')
     >            ' # body stns. to add, and X of last one',
     >            ' (N and X; 0 and X means none): '
               READ (IKBD, *) NBASE, XLAST
               IF (NBASE > 0) DX = (XLAST - XF(I)) / REAL (NBASE)
            ELSE
               NBASE = 0
            END IF

            OPEN (UNIT=IGEOM, FILE='aero.body', STATUS='NEW')

            WRITE (IGEOM, '(A)') '    FNF    FLIP'
            FN = NFSTN + NBASE
            WRITE (IGEOM, '(2F7.0)') FN, ZERO

            DO I = 1, NFSTN
               WRITE (IGEOM, '(A)')
     >            '      FNFP        XF      FSEC'
               FN = NFJ(I)
               WRITE (IGEOM, '(3F10.4)') FN, XF(I), ONE
               WRITE (IGEOM, '(A)') '    Z         Y'
               WRITE (IGEOM, '(2F10.5)')
     >            (ZFINIT(J,I), YFINIT(J,I), J = 1, NFJ(I))
            END DO

            DO I = 1, NBASE
               WRITE (IGEOM, '(A)')
     >            '      FNFP        XF      FSEC'
               XFI = XF(NFSTN) + DX * REAL (I)
               WRITE (IGEOM, '(3F10.4)') FN, XFI, ONE
               WRITE (IGEOM, '(A)') '    Z         Y'
               SCALE = MAX (REAL (NBASE - I) / REAL (NBASE), 0.01)
               DO J = 1, NFJ(NFSTN)
                  YFI = SCALE * (YFINIT(J,NFSTN) - YM) + YM
                  ZFI = SCALE *  ZFINIT(J,NFSTN)
                  WRITE (IGEOM, '(2F10.5)') ZFI, YFI
               END DO
            END DO

            CLOSE (IGEOM)
            WRITE (ICRT, '(/, 2A)')
     >         ' Fuselage (possibly fudged and/or closed)',
     >         ' written to aero.body in aero.geo format.'
            GO TO 999

         END IF

         IF (NNFUSE == 3) THEN ! Kludge to save raw fuselage points

            J = NFJ(1)
            IF (J == 1) THEN
               J = NFJ(2)
               NFJ(1) = J
               YFINIT(1:J,1) = YFINIT(1,1)
               ZFINIT(1:J,1) = ZFINIT(1,1)
            END IF

            DO I = 1, NFSTN

               IF (NFJ(I) /= J) THEN
                  WRITE (IWRIT, '(/, (A, I5))')
     >               ' Body point count mismatch at section', I,
     >               ' NFJ(I):', NFJ(I), ' NFJ(I-1):', J
                  GO TO 999
               END IF

            END DO

            OPEN (UNIT=IGEOM, FILE='fuselage.plot3d', STATUS='NEW')

            WRITE (IGEOM, '(I1)') 1
            WRITE (IGEOM, '(3I5)') J, NFSTN, 1
            WRITE (IGEOM, '(6F12.6)')
     >         ((XF(I),       J = 1, NFJ(1)), I = 1, NFSTN),
     >         ((YFINIT(J,I), J = 1, NFJ(1)), I = 1, NFSTN),
     >         ((ZFINIT(J,I), J = 1, NFJ(1)), I = 1, NFSTN)

            CLOSE (IGEOM)
            WRITE (IWRIT, '(/, A)')
     >         ' Raw fuselage sections written to fuselage.plot3d.'
            GO TO 999

         END IF

      END IF


C     Perform the intersections and paneling:

C *** IWRIT = 6                 ! Not needed here; in the design code,
C *** IF (IPROC > 1) IWRIT = -6 ! this avoids COMMON /PRC/ IPROC, NPROC

      CALL SECOND (CPU1)

      CALL SURFER
     >
     >  (NCALL, NWSURF, NWCRANK, MXIGRID, MXKGRID, MXINTR, MXPATCH,
     >   NNCONFIG, NNPANEL, NNWING, NNNAC, NNFUSE, NNSPLIT, NNROOT,
     >   NNTIP, NNTAIL, NNSURF, NWSEC, IQUAD, XTRAP, ZINNER, ZOUTER,
     >   LENRG, IDIM, KDIM, IRG, XRG, YRG, ZRG,
     >   NFSTN, JLBODY, XF, YF, ZF, IDEGFUS, IDEGREE,
     >   ILWING, KLWING, KWING1, KWING2,
     >   NCRANK, LCRANK, KCRANK, MCRANK, ZCRANK,
     >   SPNSPC, SPNSPS, TNORM, TBASE,
     >   TINTR, UINTR, VINTR, IINTR, JINTR, KINTR,
     >   NFIXL, NFIXR, NFIXA, NFIXB,
     >   IMXCOM, JMXCOM, IMXCOMO, JMXCOMO,
     >   MXSUBI, MXSUBJ, NSUBI, NSUBJ, ISUB, JSUB, ISUBO, JSUBO,
     >   LENXB, XBASE, YBASE, ZBASE, NEWPTS, NPATCH, LASTXB, MOSTXB,
     >   IXB, IPIECE, IWRIT, FAIL)

      CALL SECOND (CPU2)
      CPU2 = CPU2 - CPU1

      WRITE (IWRIT, '(/, A, F10.2)')
     >   ' CPU seconds to panel geometry:', CPU2


      IF (NNEXPAND /= 0) THEN ! Expand subpatches as in SYN107-MB

         CALL SECOND (CPU1)

         CALL EXPAND_PATCHES (MXPATCH, NPATCH, IMXCOM, JMXCOM,
     >                        NSUBI, NSUBJ, MXSUBI, MXSUBJ,
     >                        ISUB, JSUB, LENXB, LASTXB, MOSTXB, IXB,
     >                        IPIECE, XBASE, YBASE, ZBASE, IWRIT)

         CALL SECOND (CPU2)
         CPU2 = CPU2 - CPU1

         WRITE (IWRIT, '(/, A, F10.2)')
     >      ' CPU seconds to expand subpatches:', CPU2

      END IF


      WRITE (ICRT, '(A, I3)') ' Patches generated:', NPATCH

      WRITE (IWRIT, '(/, A, 30X, 2A)')
     >   ' Patch #   Dimensions    Component', 'Mean Xlength & Xmin',
     >   '       Area        Subpatch Indices (I & J)'

      AREA_TOTAL = ZERO
      J = 0

      DO N = 1, NPATCH

C        Patch areas and mean lengths may help viscous drag estimates:

         I = IXB(N)

         CALL PATCH_AREA (IMXCOM(N), JMXCOM(N), XBASE(I), YBASE(I),
     >                    ZBASE(I), XMIN, XLENGTH, AREA_PATCH)

         AREA_TOTAL = AREA_TOTAL + AREA_PATCH

C        Assemble the subpatch indices as a character string:

         TITLE = ' '
         M = 0
         DO I = 0, NSUBI(N)
            M = M + 4
            WRITE (TITLE(M-3:M), '(I4)') ISUB(I,N)
         END DO

         M = M + 4
         DO I = 0, NSUBJ(N)
            M = M + 4
            WRITE (TITLE(M-3:M), '(I4)') JSUB(I,N)
         END DO

         I = IPIECE(N)

         IF (I /= J) THEN
            J = I
            L = LEN_TRIM (PIECE(I))
            WRITE (IWRIT,
     >         '(/, I6, I8, I7, I5, 7X, A, T66, 2F9.3, 1P, E16.7, A)')
     >         N, IMXCOM(N), JMXCOM(N), I, PIECE(I)(1:L), XLENGTH, XMIN,
     >         AREA_PATCH, TITLE(1:M)
         ELSE
            WRITE (IWRIT, '(I6, I8, I7, I5, T66, 2F9.3, 1P, E16.7, A)')
     >         N, IMXCOM(N), JMXCOM(N), I, XLENGTH, XMIN, AREA_PATCH,
     >         TITLE(1:M)
         END IF

      END DO

      WRITE (IWRIT, '(/, A, 1P, E16.7)')
     >   ' Total wetted area:    ', AREA_TOTAL
      WRITE (IWRIT, '(/, (A, I10))')
     >   ' Total number of points:     ', LASTXB,
     >   ' Room needed in X/Y/ZBASE(*):', MOSTXB


C     Repeat the paneling as a check?

      IF (NNREPEAT /= 0) THEN

         NNREPEAT = 0
         NCALL = -2
         GO TO 500

      END IF


C     Save results:

      IF (MOD (NNFORMAT, 2) == 0) THEN ! Unformatted patches

         OPEN (UNIT=ISAVE, FILE='aero.xyz', STATUS='UNKNOWN',
     >         FORM='UNFORMATTED')

         WRITE (ISAVE) NPATCH
         WRITE (ISAVE) (IMXCOM(N), JMXCOM(N), 1, N = 1, NPATCH)

      ELSE ! Formatted patches

         OPEN (UNIT=ISAVE, FILE='aero.xyz', STATUS='UNKNOWN')

         WRITE (ISAVE, '(I4)') NPATCH
         WRITE (ISAVE, '(8 (2I4, I2))')
     >      (IMXCOM(N), JMXCOM(N), 1, N = 1, NPATCH)

      END IF

      DO N = 1, NPATCH

         I1 = IXB(N)

         CALL XYZ_WRITE (NNFORMAT, ISAVE, IMXCOM(N), JMXCOM(N),
     >                   IMXCOM(N), JMXCOM(N),
     >                   XBASE(I1), YBASE(I1), ZBASE(I1))
      END DO

C     Append component numbers, etc., to help CFD surface grid mapping:

      N = NPATCH ! May be MXPATCH - 1
      L = ISAVE

      IF (MOD (NNFORMAT, 2) == 0) THEN
         WRITE (L) IPIECE(1:N)
         WRITE (L) MXSUBI, MXSUBJ
         WRITE (L) NSUBI(1:N)
         WRITE (L) NSUBJ(1:N)
         WRITE (L) ISUB(0:MXSUBI,1:N)
         WRITE (L) JSUB(0:MXSUBJ,1:N)
      ELSE
         WRITE (L, '(A, /, (26I3))') 'Component numbers:', IPIECE(1:N)
         WRITE (L, '(A, /,   2I6)' ) 'MXSUBI  MXSUBJ:', MXSUBI, MXSUBJ
         WRITE (L, '(A, /, (26I3))') 'NSUBI(1:NPATCH):', NSUBI(1:N)
         WRITE (L, '(A, /, (26I3))') 'NSUBJ(1:NPATCH):', NSUBJ(1:N)
         WRITE (L, '(A, /, (13I6))') 'Subpatch Is:', ISUB(0:MXSUBI,1:N)
         WRITE (L, '(A, /, (13I6))') 'Subpatch Js:', JSUB(0:MXSUBJ,1:N)
      END IF

      CLOSE (L)

  999 CONTINUE

      END PROGRAM AEROSURF

C***********************************************************************
C
      SUBROUTINE CONTROL (IREAD,   IWRIT,   MXWSURF, MXKWING, MXCRANK,
     >                    LENXB,   PIECE,
     >                    NNCONFIG,NNWING,  NNNAC,   NNFUSE,  NNSPLIT,
     >                    NNFORMAT,NNFUDGE, NNREPEAT,NNINSERT,NNKAPPA,
     >                    NNPANEL, NNEXPAND,NNROOT,  NNTAIL,
     >                    ILWING,  JLBODY,  KLWING,  KWING1,  KWING2,
     >                    IROLL,   IFLAP,   IDIHED,  IDEGREE, IDEGFUS,
     >                    IBOD1,   IBOD2,   IBOD3,
     >                    MXSHELF, NSHELF,  XSHELF,  YSHELF,  ZSHELF,
     >                    XTRAP,   ZINNER,  ZOUTER,
     >                    CRDSPL,  CRDSPQ,  CRDSPS,  CRDSPC,
     >                    SPNSPC,  SPNSPS,  NFIXL,   NFIXR,
     >                    NFIXA,   NFIXB,   NCRANK,  LCRANK,  ZCRANK)
C
C     Read AEROSURF control inputs.  See descriptions in the main program.
C
C     01/29/99  DAS  Most error checking removed for dynamic version.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   IREAD, IWRIT, MXWSURF, MXKWING, MXCRANK

      INTEGER, INTENT (OUT) ::
     >   LENXB

      CHARACTER, INTENT (OUT) ::
     >   PIECE(MXWSURF+1)*(*)

      INTEGER, INTENT (OUT) ::
     >   NNCONFIG, NNWING, NNNAC, NNSPLIT, NNFUSE, NNFORMAT, NNFUDGE,
     >   NNREPEAT, NNINSERT, NNKAPPA, NNPANEL, NNEXPAND,
     >   NNROOT(MXWSURF), NNTAIL(MXWSURF), ILWING(MXWSURF), JLBODY

      INTEGER, INTENT (OUT), DIMENSION (MXWSURF) ::
     >   KLWING, KWING1, KWING2, IROLL, IFLAP, IDIHED, IDEGREE

      INTEGER, INTENT (OUT) ::
     >   IDEGFUS, IBOD1, IBOD2, IBOD3

      INTEGER, INTENT (IN) ::
     >   MXSHELF

      INTEGER, INTENT (OUT) ::
     >   NSHELF

      REAL, INTENT (OUT), DIMENSION (MXSHELF) ::
     >   XSHELF, YSHELF, ZSHELF

      REAL, INTENT (OUT), DIMENSION (MXWSURF) ::
     >   XTRAP, ZINNER, ZOUTER,
     >   CRDSPL, CRDSPQ, CRDSPS, CRDSPC, SPNSPC, SPNSPS

      INTEGER, INTENT (OUT), DIMENSION (MXWSURF) ::
     >   NFIXL, NFIXR, NFIXA, NFIXB, NCRANK, LCRANK

      REAL, INTENT (OUT) ::
     >   ZCRANK(MXCRANK,MXWSURF)

C     Local variables:

      INTEGER
     >   K, LUN, N, NC

      CHARACTER
     >   TEXT * 10

C     Execution:

      LUN = IREAD
      READ (LUN,*)
      READ (LUN,*)
      READ (LUN,*) NNCONFIG
      READ (LUN,*)
      READ (LUN,*)
      READ (LUN,*) NNWING, NNNAC, NNFUSE, NNSPLIT, NNFORMAT,
     >             NNREPEAT, NNINSERT, NNKAPPA

      IF (NNWING + NNNAC > MXWSURF) THEN
         WRITE (IWRIT, '(/, A)') ' NNWING + NNNAC > MXWSURF.'
         GO TO 900
      END IF

      READ (LUN,*)
      READ (LUN,*) NNPANEL, NNEXPAND, LENXB

      DO N = 1, NNWING + NNNAC
         READ (LUN,*)
         READ (LUN, '(A)') PIECE(N)
         READ (LUN,*)
         IF (N <= NNWING) THEN
            READ (LUN,*) NNROOT(N), NNTAIL(N), ILWING(N), KLWING(N),
     >                   IROLL(N), IFLAP(N), IDIHED(N)
         ELSE
            READ (LUN,*) NNROOT(N), ILWING(N), KLWING(N), IROLL(N)
         END IF
         READ (LUN,*)
         READ (LUN,*) KWING1(N), KWING2(N), IDEGREE(N), XTRAP(N),
     >                ZINNER(N), ZOUTER(N)
         READ (LUN,*)
         READ (LUN,*) NFIXL(N), NFIXR(N), NFIXA(N), NFIXB(N)
         READ (LUN,*)
         READ (LUN,*) NC, LCRANK(N),
     >                (ZCRANK(K,N), K = 1, MIN (MXCRANK, NC))

         IF (NC > MXCRANK) THEN
            WRITE (IWRIT, '(/, (1X, A, I4))')
     >         'Too many cranks specified for wing surface', N,
     >         'Input:', NC, 'Limit:', MXCRANK, 'Proceeding ...'
            NC = MXCRANK
         END IF

         NCRANK(N) = NC
         LCRANK(N) = MIN (NC + 1, LCRANK(N)) ! Wing tip becomes a crank

         READ (LUN,*)
         READ (LUN,*) CRDSPL(N), CRDSPQ(N), CRDSPS(N), CRDSPC(N),
     >                SPNSPC(N), SPNSPS(N)
      END DO

C     Read fuselage data whether there is a fuselage or not, because
C     shape function inputs now follow the fuselage inputs.

      N = NNWING + NNNAC + 1
      READ (LUN,*)
      READ (LUN,'(A)') PIECE(N)
      READ (LUN,*)
      READ (LUN,*) JLBODY, IDEGFUS, NNFUDGE, IBOD1,IBOD2, IBOD3, NSHELF

      IF (NSHELF < 2 .OR. NSHELF > MXSHELF) THEN
         WRITE (IWRIT, '(/, A)') ' Must have 2 <= NSHELF <= MXSHELF.'
         GO TO 900
      END IF

      READ (LUN,*) TEXT, XSHELF(1:NSHELF)
      READ (LUN,*) TEXT, YSHELF(1:NSHELF)
      READ (LUN,*) TEXT, ZSHELF(1:NSHELF)

      RETURN

  900 STOP

      END SUBROUTINE CONTROL

C***********************************************************************
C
      SUBROUTINE COPY_GEOM (NWSURF, NNFUSE, NIWING, NKWING, NWCRANK,
     >                      NFSTN,  NJFUSE)
C
C     Transfer geometry data from interim GEOM0 arrays to GEOM1 arrays
C     (part of avoiding reading all the geometry components twice).
C
C***********************************************************************

C     Global variables:

      USE GEOM0
      USE GEOM1

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   NWSURF, NNFUSE, NIWING, NKWING, NWCRANK, NFSTN, NJFUSE

C     Execution:

      IQUAD   =   IQUAD0(1:2,1:NWSURF)

      KCRANK  =  KCRANK0(1:NWCRANK,1:NWSURF)
      ZCRANK  =  ZCRANK0(1:NWCRANK,1:NWSURF)

      ILWING  =  ILWING0(1:NWSURF)
      KLWING  =  KLWING0(1:NWSURF)
      IDEGREE = IDEGREE0(1:NWSURF)
      IROLL   =   IROLL0(1:NWSURF)
      IFLAP   =   IFLAP0(1:NWSURF)
      IDIHED  =  IDIHED0(1:NWSURF)
      NWSEC   =   NWSEC0(1:NWSURF)
      KWING1  =  KWING10(1:NWSURF)
      KWING2  =  KWING20(1:NWSURF)
      NCRANK  =  NCRANK0(1:NWSURF)
      LCRANK  =  LCRANK0(1:NWSURF)
      NFIXA   =   NFIXA0(1:NWSURF)
      NFIXB   =   NFIXB0(1:NWSURF)
      NFIXL   =   NFIXL0(1:NWSURF)
      NFIXR   =   NFIXR0(1:NWSURF)
      NNROOT  =  NNROOT0(1:NWSURF)
      NNTAIL  =  NNTAIL0(1:NWSURF)
      NNTIP   =   NNTIP0(1:NWSURF)

      NLG     =  NLG0(1:NKWING,1:NWSURF)
      NLGI    =  NLG                    ! Flap/slat deflections forced this
      NUG     =  NUG0(1:NKWING,1:NWSURF)
      NUGI    =  NUG
      NWG     =  NWG0(1:NKWING,1:NWSURF)
      NWGI    =  NWG

      CRDSPL  =  CRDSPL0(1:NWSURF)
      CRDSPQ  =  CRDSPQ0(1:NWSURF)
      CRDSPS  =  CRDSPS0(1:NWSURF)
      CRDSPC  =  CRDSPC0(1:NWSURF)
      SPNSPC  =  SPNSPC0(1:NWSURF)
      SPNSPS  =  SPNSPS0(1:NWSURF)
      XC1     =     XC10(1:NWSURF)
      YC1     =     YC10(1:NWSURF)
      ZC1     =     ZC10(1:NWSURF)
      XC2     =     XC20(1:NWSURF)
      YC2     =     YC20(1:NWSURF)
      ZC2     =     ZC20(1:NWSURF)
      CLOCK   =   CLOCK0(1:NWSURF)
      XTRAP   =   XTRAP0(1:NWSURF)
      ZINNER  =  ZINNER0(1:NWSURF)
      ZOUTER  =  ZOUTER0(1:NWSURF)
      PIECE   =   PIECE0(1:NWSURF+1)

      AINIT   =   AINIT0(1:NKWING,1:NWSURF)
      ALG     =     ALG0(1:NKWING,1:NWSURF)
      CINIT   =   CINIT0(1:NKWING,1:NWSURF)
      DFLAP   =   DFLAP0(1:NKWING,1:NWSURF)
      DSLAT   =   DSLAT0(1:NKWING,1:NWSURF)
      DIHED   =   DIHED0(1:NKWING,1:NWSURF)
      DIAXIS  =  DIAXIS0(1:3,1:2,1:NKWING,1:NWSURF)
      ROLL    =    ROLL0(1:NKWING,1:NWSURF)
      ROUNDED = ROUNDED0(1:NKWING,1:NWSURF)
      SWFLAP  =  SWFLAP0(1:NKWING,1:NWSURF)
      SWSLAT  =  SWSLAT0(1:NKWING,1:NWSURF)
      TINIT   =   TINIT0(1:NKWING,1:NWSURF)
      XFLAP   =   XFLAP0(1:NKWING,1:NWSURF)
      XSLAT   =   XSLAT0(1:NKWING,1:NWSURF)
      XLEADI  =  XLEADI0(1:NKWING,1:NWSURF)
      YLEADI  =  YLEADI0(1:NKWING,1:NWSURF)
      ZLEADI  =  ZLEADI0(1:NKWING,1:NWSURF)

      XINIT   =   XINIT0(1:NIWING,1:NKWING,1:NWSURF)
      YINIT   =   YINIT0(1:NIWING,1:NKWING,1:NWSURF)
      ZINIT   =   ZINIT0(1:NIWING,1:NKWING,1:NWSURF)
      XWG     =     XWG0(1:NIWING,1:NKWING,1:NWSURF)
      YWG     =     YWG0(1:NIWING,1:NKWING,1:NWSURF)
      ZWG     =     ZWG0(1:NIWING,1:NKWING,1:NWSURF)

      IF (NNFUSE /= 0) THEN
         NFJ      =    NFJ0(1:NFSTN)
         XF       =     XF0(1:NFSTN)
         YFINIT   = YFINIT0(1:NJFUSE,1:NFSTN)
         ZFINIT   = ZFINIT0(1:NJFUSE,1:NFSTN)
      END IF

      END SUBROUTINE COPY_GEOM

C***********************************************************************
C
      SUBROUTINE COUNT_PATCHES (NNCONFIG, NNWING,  NNNAC,  NNFUSE,
     >                          NNSPLIT,  NNROOT,  NNTAIL, NNEXPAND,
     >                          NWSURF,   NWCRANK, ILWING, XTRAP,
     >                          MXSUBI,   MXSUBJ,  NPATCH, IWRIT)
C
C     Count the number of surface patches which SURFER will produce.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   NWSURF,               ! MAX (NNWING + NNNAC, 1)
     >   NNCONFIG,             ! 0 except for CTV (1, 2, or 3)
     >   NNWING,               ! # wing-type surface, >= 0
     >   NNNAC,                ! # nacelles, >= 0
     >   NNFUSE,               ! 0 means no fuselage is present
     >   NNSPLIT,              ! Type of fuselage panel splitting
     >   NNROOT(NWSURF),       ! See lengthy definition in AEROSURF or SURFER
     >   NNTAIL(NWSURF),       !  "   "  ...
     >   NNEXPAND,             ! 1 means allow for expanding the subpatches
     >   NWCRANK,              ! Max. # wing cranks + 1 for the tip
     >   ILWING(NWSURF)        ! Needed for vertical fin = tail pt. count check

      REAL, INTENT (IN) ::
     >   XTRAP(NWSURF)         ! For a diverter, XTRAP > 0. extends wing TE;
                               ! For a nacelle,  XTRAP < 0. indicates it is
                               ! attached to an under-wing diverter, not a pylon
      INTEGER, INTENT (OUT) ::
     >   MXSUBI, MXSUBJ,       ! Max. # subpatches expected over all patches
     >   NPATCH                ! Total # patches SURFER will produce

      INTEGER, INTENT (IN) ::
     >   IWRIT

C     Local variables:

      INTEGER
     >   IW, M, N, NFIN0, NFIN1, NFINT, NWBINT, NWPINT(NWSURF)
      LOGICAL
     >   DIVERTERS, PYLONS

C     Execution:

      NWBINT = 0 ! # wing/body intersections
      NWPINT = 0 ! NWPINT(IW) = # wing/pylon intersections for wing IW
      NFIN0  = 0 ! # fin/body intersections with no associated wing/body int.
      NFIN1  = 0 ! # fin/body intersections with associated wing/body int.
      NFINT  = 0 ! # fin/tail intersections
      PYLONS = .FALSE.
      DIVERTERS = .FALSE.

C     Count the wing/body & wing/pylon intersections:

      DO IW = 1, NNWING

         N = NNROOT(IW)

         IF (N <= -3) THEN ! Tail IW intersects vertical fin |N|
            NFINT = NFINT + 1
         ELSE IF (N == -2) THEN ! Vertical fin IW intersects the body
            IF (NNFUSE /= 0) THEN
               M = NNTAIL(IW)
               IF (M == 0) THEN ! 0 = no associated tail on the body
                  NFIN0 = NFIN0 + 1
               ELSE IF (NNCONFIG == 0) THEN ! Point counts are required to match
                  M = ABS (M)
                  IF (ILWING(M) /= ILWING(IW)) THEN
                     IF (IWRIT > 0) THEN
                        WRITE (IWRIT, '(/, A, /, (A, I3, A, I4))')
     >                     ' Vertical fin & associated tail mismatch:',
     >                     ' Fin  surface:', IW, '  # pts:', ILWING(IW),
     >                     ' Tail surface:',  M, '  # pts:', ILWING(M)
                     END IF
                     CALL SYNCH
                     STOP
                  END IF
                  NFIN1 = NFIN1 + 1
               END IF
            END IF
         ELSE IF (N == -1) THEN
            IF (NNFUSE /= 0) NWBINT = NWBINT + 1
         ELSE IF (N > 0) THEN ! It's an under-wing pylon or diverter
            NWPINT(N) = NWPINT(N) + 1
            IF (XTRAP(IW) > 0.) DIVERTERS = .TRUE. ! > 0. means extend wing TE
         END IF

      END DO

C     Count the patches and subpatches:

      MXSUBI = 1
      MXSUBJ = 1

      IF (NNFUSE == 0) THEN
         NPATCH = 0
      ELSE IF (NWBINT == 0) THEN
         NPATCH = 2 + NFIN0 + NFINT ! Allow for possible base patch
         MXSUBI = 2 * NFIN0 + 1     ! For the body
      ELSE
         IF (NNCONFIG == 0) THEN ! All cases prior to CTV
            NPATCH = 3 * NWBINT + 1 + NFIN0 + NFIN1 + NFINT
            MXSUBI = 2 * NFIN0 + 1 ! At worst, on the body

            IF (NNSPLIT == 1) THEN ! Split full length for every intersection
               MXSUBJ = NWBINT + 1
            ELSE ! Just connect forward root TE to following root LE
               MXSUBJ = 2
            END IF
         ELSE IF (NNCONFIG == 1) THEN ! 1st CTV case of special body paneling
            NPATCH = 5
            MXSUBJ = MAX (MXSUBJ, 3)
         ELSE IF (NNCONFIG == 2) THEN ! 2nd CTV case of special body paneling
            NPATCH = 6
            MXSUBI = MAX (MXSUBI, 2)
            MXSUBJ = MAX (MXSUBJ, 3)
         ELSE IF (NNCONFIG == 3) THEN ! 2nd CTV case of special body paneling
            NPATCH = 6                ! Includes fin panel not counted above
            MXSUBI = MAX (MXSUBI, 2)
            MXSUBJ = MAX (MXSUBJ, 2)
         END IF
         NPATCH = NPATCH + 1 ! For possible fuselage base patch
      END IF

C     The fin-mounted tail case is uncertain:

      IF (NFINT > 0) MXSUBI = MAX (MXSUBI, 3) ! ?

C     Add any wing/pylon/diverter patches:

      IF (NNWING > 0) MXSUBJ = MAX (MXSUBJ, NWCRANK) ! Assume REPACK is used

      DO IW = 1, NNWING

         IF (NNROOT(IW) == -2) CYCLE ! Fin half-surface is counted above

         M = NWPINT(IW)
         NPATCH = NPATCH + 2 * (M + 1) ! Explicit split at LE hard to avoid

         IF (M > 0) THEN ! Not a plain surface
            IF (DIVERTERS) THEN
               MXSUBI = MAX (MXSUBI, 2)
               NPATCH = NPATCH + M ! Upper surface patches <-> diverter cutouts
            ELSE ! Wing-mounted pylon
               MXSUBI = MAX (MXSUBI, 3)
            END IF
         END IF

      END DO

C     Add any nacelle patches:

      DO IW = NNWING + 1, NNWING + NNNAC

         MXSUBI = MAX (MXSUBI, 2)

         IF (NNROOT(IW) == 0) THEN ! Plain nacelle; explicit split at LE
            NPATCH = NPATCH + 2
         ELSE ! Pylon or diverter attached
            NPATCH = NPATCH + 5         ! Interior + 4 exterior, normally
            IF (XTRAP(IW) >= 0.) THEN   ! Pylon attached, not diverter
               PYLONS = .TRUE.          ! NPATCH may be 1 bigger than necessary
               MXSUBJ = MAX (MXSUBJ, 2) ! Implicit nacelle split at pylon LE
            ELSE
               MXSUBJ = MAX (MXSUBJ, 3) ! Upper aft nacelle cap adds 1
            END IF
         END IF

      END DO

      IF (NNEXPAND /= 0) THEN ! Allow for expanding the subpatches

         NPATCH = NPATCH * (MXSUBI * MXSUBJ) ! Upper bound

      END IF

      END SUBROUTINE COUNT_PATCHES

C***********************************************************************
C
      SUBROUTINE DEFLECT (IWRIT, IW, N, VTYPE, DEGREES,
     >                    KMIN, KMAX, XWMIN, XWMAX,
     >                    CHORD, XLEAD, ZLEAD,
     >                    DFLAP, XFLAP, SWFLAP,
     >                    DSLAT, XSLAT, SWSLAT)
C
C     Translate a design variable flap or slat rotation to the equivalent
C     definition carried along with the corresponding wing sections.
C
C     08/23/99  DAS  Initial reconciliation of design variable definition
C                    with geometry file specification of deflections.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   IWRIT,                ! For error messages, if > 0
     >   IW,                   ! Wing component number,  "   "
     >   N                     ! Design variable number, for diagnostics

      CHARACTER, INTENT (IN) ::
     >   VTYPE * (*)           ! Character 1 is 'F' or 'S'

      REAL, INTENT (IN) ::
     >   DEGREES               ! Deflection angle (> 0 is down)

      INTEGER, INTENT (IN) ::
     >   KMIN, KMAX            ! First and last control surface station

      REAL, INTENT (IN) ::
     >   XWMIN, XWMAX          ! X/Cs of hinge end points

      REAL, INTENT (IN), DIMENSION (KMAX) ::
     >   CHORD,                ! Section chords
     >   XLEAD, ZLEAD          ! Leading edge (X,Z)s

      REAL, INTENT (INOUT), DIMENSION (KMAX) ::
     >   DFLAP, XFLAP, SWFLAP, ! Equivalent chordwise rotation, ...
     >   DSLAT, XSLAT, SWSLAT  ! ... X of hinge in X/ZLEAD units, ...
                               ! ... and sweep of hinge line in degrees

C     Local constants:

      REAL, PARAMETER ::
     >   RAD2DEG = 57.29577951, ZERO = 0.

C     Local variables:

      INTEGER
     >   K

      REAL
     >   COSSWEEP, DINIT, DX, DZ, R, SINIT, SWEEP, TOL,
     >   XHI, XHK, XHO, XINIT, ZHI, ZHO

      LOGICAL
     >   FLAP

C     Execution:

      TOL = (CHORD(KMIN) + CHORD(KMAX)) * 0.00005
      XHI = XLEAD(KMIN) + XWMIN * CHORD(KMIN)
      XHO = XLEAD(KMAX) + XWMAX * CHORD(KMAX)
      DX  = XHO - XHI
      ZHI = ZLEAD(KMIN)
      ZHO = ZLEAD(KMAX)
      DZ  = ZHO - ZHI
      R   = DX / DZ

      COSSWEEP = DZ / SQRT (DX ** 2 + DZ ** 2)
      SWEEP    = RAD2DEG * ACOS (COSSWEEP)
      FLAP     = VTYPE(1:1) == 'F'

      DO K = KMIN, KMAX

         XHK = XHI + R * (ZLEAD(K)  - ZHI)

         IF (FLAP) THEN
            DINIT = DFLAP(K)
            XINIT = XFLAP(K)
            SINIT = SWFLAP(K)
         ELSE
            DINIT = DSLAT(K)
            XINIT = XSLAT(K)
            SINIT = SWSLAT(K)
         END IF

         IF (DINIT /= ZERO .AND. IWRIT > 0) THEN
            IF (ABS (XHK - XINIT) > TOL) THEN
               WRITE (IWRIT, '(/, A, 3(A, I4), /, (A, 3F12.6))')
     >            ' DEFLECT warning:',
     >            ' Component #', IW, '   Section #', K,
     >            '   Design variable #', N,
     >            ' Initial angle, hinge X, sweep: ', DINIT,XINIT,SINIT,
     >            ' Current angle, hinge X, sweep: ', DEGREES,XHK,SWEEP
            END IF
         END IF

         IF (FLAP) THEN
            DFLAP(K)  = DEGREES
            XFLAP(K)  = XHK
            SWFLAP(K) = SWEEP
         ELSE
            DSLAT(K)  = DEGREES
            XSLAT(K)  = XHK
            SWSLAT(K) = SWEEP
         END IF

      END DO

      END SUBROUTINE DEFLECT

C***********************************************************************
C
      SUBROUTINE DYNAMIC_GEOM0 (MODE)
C
C     Allocate or deallocate interim work-space which avoids reading the
C     geometry data twice.
C
C***********************************************************************

C     Global variables:

      USE GEOM0

C     Argument:

      INTEGER, INTENT (IN) :: MODE ! 1 = allocate; 2 = deallocate

C     Execution:

      IF (MODE == 1) THEN

         ALLOCATE (IQUAD0(2,MXWSURF), KCRANK0(MXCRANK,MXWSURF),
     >             NFJ0(MXIBODY))

         ALLOCATE (ILWING0(MXWSURF), KLWING0(MXWSURF),IDEGREE0(MXWSURF),
     >             IROLL0(MXWSURF),  IFLAP0(MXWSURF),  IDIHED0(MXWSURF),
     >             NWSEC0(MXWSURF),  KWING10(MXWSURF), KWING20(MXWSURF),
     >             NCRANK0(MXWSURF), LCRANK0(MXWSURF),
     >             NFIXA0(MXWSURF),  NFIXB0(MXWSURF),
     >             NFIXL0(MXWSURF),  NFIXR0(MXWSURF),
     >             NNROOT0(MXWSURF), NNTAIL0(MXWSURF), NNTIP0(MXWSURF))

         ALLOCATE (NLG0(MXKWING,MXWSURF), NUG0(MXKWING,MXWSURF),
     >             NWG0(MXKWING,MXWSURF))

         ALLOCATE (CRDSPL0(MXWSURF), CRDSPQ0(MXWSURF),
     >             CRDSPS0(MXWSURF), CRDSPC0(MXWSURF),
     >             SPNSPC0(MXWSURF), SPNSPS0(MXWSURF),
     >             XC10(MXWSURF),    YC10(MXWSURF),    ZC10(MXWSURF),
     >             XC20(MXWSURF),    YC20(MXWSURF),    ZC20(MXWSURF),
     >             CLOCK0(MXWSURF),  XTRAP0(MXWSURF),
     >             ZINNER0(MXWSURF), ZOUTER0(MXWSURF),
     >             PIECE0(MXWSURF+1)) ! Allow for fuselage

         ALLOCATE (AINIT0(MXKWING,MXWSURF),   ALG0(MXKWING,MXWSURF),
     >             CINIT0(MXKWING,MXWSURF),   DFLAP0(MXKWING,MXWSURF),
     >             DIHED0(MXKWING,MXWSURF),
     >             DIAXIS0(3,2,MXKWING,MXWSURF),
     >             DSLAT0(MXKWING,MXWSURF),   ROLL0(MXKWING,MXWSURF),
     >             ROUNDED0(MXKWING,MXWSURF), SWFLAP0(MXKWING,MXWSURF),
     >             SWSLAT0(MXKWING,MXWSURF),  TINIT0(MXKWING,MXWSURF),
     >             XFLAP0(MXKWING,MXWSURF),   XSLAT0(MXKWING,MXWSURF),
     >             XLEADI0(MXKWING,MXWSURF),  YLEADI0(MXKWING,MXWSURF),
     >             ZLEADI0(MXKWING,MXWSURF))

         ALLOCATE (ZCRANK0(MXCRANK,MXWSURF))

         ALLOCATE (XF0(MXIBODY),
     >             YFINIT0(MXJBODY,MXIBODY), ZFINIT0(MXJBODY,MXIBODY))

         ALLOCATE (XINIT0(MXIWING,MXKWING,MXWSURF),
     >             YINIT0(MXIWING,MXKWING,MXWSURF),
     >             ZINIT0(MXIWING,MXKWING,MXWSURF),
     >             XWG0(MXIWING,MXKWING,MXWSURF),
     >             YWG0(MXIWING,MXKWING,MXWSURF),
     >             ZWG0(MXIWING,MXKWING,MXWSURF))

      ELSE ! MODE = 2

         DEALLOCATE (IQUAD0, KCRANK0, NFJ0, ILWING0, KLWING0, IDEGREE0,
     >               IROLL0, IFLAP0, IDIHED0, NWSEC0, KWING10, KWING20,
     >               NCRANK0, LCRANK0, NFIXA0, NFIXB0, NFIXL0, NFIXR0,
     >               NNROOT0, NNTAIL0, NNTIP0, NLG0, NUG0, NWG0,
     >               CRDSPL0, CRDSPQ0, CRDSPS0, CRDSPC0,
     >               SPNSPC0, SPNSPS0,
     >               XC10, YC10, ZC10, XC20, YC20, ZC20,
     >               CLOCK0, XTRAP0, ZINNER0, ZOUTER0,
     >               AINIT0, ALG0, CINIT0, DFLAP0, DIHED0, DIAXIS0,
     >               DSLAT0, ROLL0, ROUNDED0, SWFLAP0, SWSLAT0, TINIT0,
     >               XFLAP0, XSLAT0, XLEADI0, YLEADI0, ZLEADI0,
     >               ZCRANK0, XF0, YFINIT0, ZFINIT0,
     >               XINIT0, YINIT0, ZINIT0, XWG0, YWG0, ZWG0)
      END IF

      END SUBROUTINE DYNAMIC_GEOM0

C***********************************************************************
C
      SUBROUTINE DYNAMIC_GEOM1 (NWSURF, MXIWING, MXKWING, NWCRANK,
     >                          NFSTN, MXJBODY)
C
C     Allocate more precise work-space for the geometry read temporarily
C     into module GEOM0 arrays.  The "MX"s in the arguments are still
C     appropriate because we're not packing the data, but they are now
C     no bigger than necessary for the geometry encountered.
C
C***********************************************************************

C     Global variables:

      USE GEOM1

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   NWSURF,  ! MAX (Actual # "wing" surfaces, 1) to avoid 0 dimensions
     >   MXIWING, ! Max. # geometry points per wing-type defining section
     >   MXKWING, ! Max. # defining sections per "wing" surface
     >   NFSTN,   ! MAX (Actual # fuselage sections, 1) to avoid 0 dimensions
     >   MXJBODY, ! Max. # points per fuselage defining section
     >   NWCRANK  ! Max. # crank sections read for any wing-type component;
                  ! NWCRANK elements are allocated to allow adding the tip

C     Execution:

      ALLOCATE (IQUAD(2,NWSURF), KCRANK(NWCRANK+1,NWSURF), NFJ(NFSTN))

      ALLOCATE (ILWING(NWSURF), KLWING(NWSURF), IDEGREE(NWSURF),
     >          IROLL(NWSURF),  IFLAP(NWSURF),  IDIHED(NWSURF),
     >          NWSEC(NWSURF),  KWING1(NWSURF), KWING2(NWSURF),
     >          NCRANK(NWSURF), LCRANK(NWSURF),
     >          NFIXA(NWSURF),  NFIXB(NWSURF),
     >          NFIXL(NWSURF),  NFIXR(NWSURF),
     >          NNROOT(NWSURF), NNTAIL(NWSURF), NNTIP(NWSURF))

      ALLOCATE (NLG(MXKWING,NWSURF), NUG(MXKWING,NWSURF),
     >          NWG(MXKWING,NWSURF))

      ALLOCATE (NLGI(MXKWING,NWSURF), NUGI(MXKWING,NWSURF),
     >          NWGI(MXKWING,NWSURF))

      ALLOCATE (CRDSPL(NWSURF), CRDSPQ(NWSURF),
     >          CRDSPS(NWSURF), CRDSPC(NWSURF),
     >          SPNSPC(NWSURF), SPNSPS(NWSURF),
     >          XC1(NWSURF),    YC1(NWSURF),    ZC1(NWSURF),
     >          XC2(NWSURF),    YC2(NWSURF),    ZC2(NWSURF),
     >          CLOCK(NWSURF),  XTRAP(NWSURF),
     >          ZINNER(NWSURF), ZOUTER(NWSURF),
     >          PIECE(NWSURF+1)) ! Allow for fuselage

      ALLOCATE (AINIT(MXKWING,NWSURF),   ALG(MXKWING,NWSURF),
     >          CINIT(MXKWING,NWSURF),   DFLAP(MXKWING,NWSURF),
     >          DIHED(MXKWING,NWSURF),   DIAXIS(3,2,MXKWING,NWSURF),
     >          DSLAT(MXKWING,NWSURF),   ROLL(MXKWING,NWSURF),
     >          ROUNDED(MXKWING,NWSURF), SWFLAP(MXKWING,NWSURF),
     >          SWSLAT(MXKWING,NWSURF),  TINIT(MXKWING,NWSURF),
     >          XFLAP(MXKWING,NWSURF),   XSLAT(MXKWING,NWSURF),
     >          XLEADI(MXKWING,NWSURF),  YLEADI(MXKWING,NWSURF),
     >          ZLEADI(MXKWING,NWSURF))

      ALLOCATE (ZCRANK(NWCRANK,NWSURF))

      ALLOCATE (XF(NFSTN),
     >          YFINIT(MXJBODY,NFSTN), ZFINIT(MXJBODY,NFSTN))

      ALLOCATE (XINIT(MXIWING,MXKWING,NWSURF),
     >          YINIT(MXIWING,MXKWING,NWSURF),
     >          ZINIT(MXIWING,MXKWING,NWSURF),
     >          XWG(MXIWING,MXKWING,NWSURF),
     >          YWG(MXIWING,MXKWING,NWSURF),
     >          ZWG(MXIWING,MXKWING,NWSURF))

      END SUBROUTINE DYNAMIC_GEOM1

C***********************************************************************
C
      SUBROUTINE DYNAMIC_GRIDS (NWSURF, NPATCH, NFSTN, JLBODY, LENRG,
     >                          NWCRANK, MXIGRID, MXKGRID, MXINTR,
     >                          MXSUBI, MXSUBJ, LENXB)
C
C     Allocate surface grid/patch-related arrays with precise dimensions.
C
C***********************************************************************

C     Global variables:

      USE GRIDS

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   NWSURF, NPATCH, NFSTN, JLBODY, LENRG, NWCRANK, MXIGRID,
     >   MXKGRID, MXINTR, MXSUBI, MXSUBJ, LENXB

C     Execution:

      ALLOCATE (IMXCOM(NPATCH),  JMXCOM(NPATCH),
     >          IDIMENS(NPATCH), JDIMENS(NPATCH),
     >          IPIECE(NPATCH),  IXB(NPATCH),
     >          NSUBI(NPATCH),   NSUBJ(NPATCH),
     >          ISUB(0:MXSUBI,NPATCH),  JSUB(0:MXSUBJ,NPATCH),
     >          IINTR(MXINTR,NWSURF),  JINTR(MXINTR,NWSURF),
     >          KINTR(MXINTR,NWSURF),  TINTR(MXINTR,NWSURF),
     >          UINTR(MXINTR,NWSURF),  VINTR(MXINTR,NWSURF),
     >          SNORM(MXIGRID,NWSURF), TBASE(MXKGRID,NWSURF),
     >          TNORM(MXKGRID,NWSURF), MCRANK(NWCRANK,NWSURF))

      ALLOCATE (YF(JLBODY,NFSTN), ZF(JLBODY,NFSTN))

      ALLOCATE (XRG(LENRG),   YRG(LENRG),   ZRG(LENRG),
     >          XBASE(LENXB), YBASE(LENXB), ZBASE(LENXB))

      END SUBROUTINE DYNAMIC_GRIDS

C*******************************************************************************
C
      SUBROUTINE EXPAND_PATCHES (MXPATCH, NPATCH, IMXCOM, JMXCOM,
     >                           NSUBI, NSUBJ, MXSUBI, MXSUBJ,
     >                           ISUB, JSUB, LENXB, LASTXB, MOSTXB, IXB,
     >                           IPIECE, XBASE, YBASE, ZBASE, IWRIT)
C
C     Convert AEROSURF paneling from subpatch form to original no-subpatch form
C     as needed by the multiblock mesh morphing scheme.  We repack in-place in
C     reverse patch order, although within one patch, a temporary copy cannot
C     be avoided.  (For a counter-example, consider repacking as two panels a
C     single 4 x 7 panel with just one implicit split at I = 3.  Reverse order
C     does not avoid overwriting elements before they have been copied.  Am I
C     missing something?)
C
C     12/31/99  DAS  Initial approach attempted to modularize MESHWARP scheme.
C     01/04/00   "   Abandoned MESHWARP approach - have to do it in-place.
C
C     Author:  David Saunders, Raytheon/NASA Ames Research Center.
C
C*******************************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   MXPATCH                          ! Max. # patches allowed for

      INTEGER, INTENT (INOUT) ::
     >   NPATCH,                          ! # patches before & after
     >   IMXCOM(MXPATCH), JMXCOM(MXPATCH) ! Patch dimensions

      INTEGER, INTENT (INOUT) ::
     >   NSUBI(MXPATCH), NSUBJ(MXPATCH)   ! # subpatches in each direction

      INTEGER, INTENT (IN) ::
     >   MXSUBI, MXSUBJ                   ! Dimensions of ISUB & JSUB

      INTEGER, INTENT (INOUT) ::
     >   ISUB(0:MXSUBI,MXPATCH),          ! Subpatch dimensions;
     >   JSUB(0:MXSUBJ,MXPATCH)           ! ISUB(0,*) = JSUB(0,*) = 1

      INTEGER, INTENT (IN) ::
     >   LENXB                            ! (Max.) length of X/Y/ZBASE(*)

      INTEGER, INTENT (INOUT) ::
     >   LASTXB,                          ! Current space used in X/Y/ZBASE
     >   MOSTXB,                          ! Most space reqd. so far (XBCHECK)
     >   IXB(MXPATCH),                    ! 1st indices of patches in X/Y/ZBASE
     >   IPIECE(MXPATCH)                  ! Component #s (for XBCHECK)

      REAL, INTENT (INOUT), DIMENSION (LENXB) ::
     >   XBASE, YBASE, ZBASE              ! Packed patches in subpatch form

      INTEGER, INTENT (IN) ::
     >   IWRIT                            ! For XBCHECK error message

C     Local variables:

      INTEGER
     >   IPC, IS, JS, K, L, NEEDED, NPATCHNEW

      INTEGER, ALLOCATABLE, DIMENSION (:) ::
     >   IDIM, IPCT, IXBT, JDIM

      REAL, ALLOCATABLE ::
     >   PANEL(:,:,:)

C     Execution:

C     Count the expanded number of patches, L, and quit if L > MXPATCH:

      L = 0

      DO K = 1, NPATCH
         L = L + NSUBI(K) * NSUBJ(K)
      END DO

      IF (L == NPATCH) GO TO 999 ! Nothing to do


      CALL XBCHECK (LENXB, MXPATCH, L, LASTXB, 0, 0, IPIECE(NPATCH),
     >              IWRIT, MOSTXB)

C     Establish the new patch dimensions, initial pointers, etc.

      ALLOCATE (IDIM(L), JDIM(L), IPCT(L), IXBT(L))

      NEEDED = 0
      L = 0

      DO K = 1, NPATCH
         IPC = IPIECE(K)
         DO JS = 1, NSUBJ(K)
            DO IS = 1, NSUBI(K)
               L = L + 1
               IPCT(L) = IPC
               IDIM(L) = ISUB(IS,K) - ISUB(IS-1,K) + 1
               JDIM(L) = JSUB(JS,K) - JSUB(JS-1,K) + 1
               IXBT(L) = NEEDED + 1
               NEEDED  = NEEDED + IDIM(L) * JDIM(L)
            END DO
         END DO
      END DO

      NEEDED = NEEDED - LASTXB ! Now it's the extra length needed

      CALL XBCHECK (LENXB, MXPATCH, L, LASTXB, NEEDED, 1, IPC,
     >              IWRIT, MOSTXB)

      LASTXB = LASTXB + NEEDED
      NPATCHNEW = L

C     The patches can be processed in reverse order to avoid clobbering data,
C     but the subpatches within a patch cannot, so we keep a copy of one
C     patch at a time.  Pushing it down a level allows (i,j) subscripting.

      DO K = NPATCH, 1, -1

         IS = IMXCOM(K)
         JS = JMXCOM(K)

         ALLOCATE (PANEL(IS,JS,3))

C        L points to the highest patch derived from patch K.

         CALL EXPAND_PANEL (NSUBI(K), NSUBJ(K), ISUB(0,K), JSUB(0,K),
     >                      IXB(K))

C        L is decremented by NSUBI(K) * NSUBJ(K) upon return.

         DEALLOCATE (PANEL)

      END DO

      NPATCH = NPATCHNEW

      DO L = 1, NPATCHNEW
         NSUBI(L)  = 1
         NSUBJ(L)  = 1
         IMXCOM(L) = IDIM(L)
         ISUB(1,L) = IDIM(L)
         JMXCOM(L) = JDIM(L)
         JSUB(1,L) = JDIM(L)
         IXB(L)    = IXBT(L)
         IPIECE(L) = IPCT(L)
      END DO

      DEALLOCATE (IDIM, JDIM, IPCT, IXBT)

  999 RETURN


C     EXPAND_PATCHES internal procedure:

      CONTAINS

         SUBROUTINE EXPAND_PANEL (NSUBI, NSUBJ, ISUB, JSUB, IP) ! Treat 1 panel

!        Arguments:

         INTEGER, INTENT (IN) ::
     >      NSUBI, NSUBJ,       ! Numbers of implicit splits for this patch
     >      ISUB(0:NSUBI),      ! Split indices for this patch
     >      JSUB(0:NSUBJ),
     >      IP                  ! Index in X/Y/ZBASE of start of this patch

!        Local variables:

         INTEGER
     >      I, J, IT, JT, LT, M

!        Execution:


         M = IP ! Make a temporary copy of the current old panel

         DO J = 1, JS
            DO I = 1, IS
               PANEL(I,J,1) = XBASE(M)
               PANEL(I,J,2) = YBASE(M)
               PANEL(I,J,3) = ZBASE(M)
               M = M + 1
            END DO
         END DO

!        Copy subpatches out of PANEL() back to X/Y/ZBASE:

         L  = L - NSUBI * NSUBJ
         LT = L

         DO JT = 1, NSUBJ
            DO IT = 1, NSUBI
               LT = LT + 1    ! Current new panel
               M  = IXBT(LT)
               DO J = JSUB(JT-1), JSUB(JT)
                  DO I = ISUB(IT-1), ISUB(IT)
                     XBASE(M) = PANEL(I,J,1)
                     YBASE(M) = PANEL(I,J,2)
                     ZBASE(M) = PANEL(I,J,3)
                     M = M + 1
                  END DO
               END DO
            END DO
         END DO

         END SUBROUTINE EXPAND_PANEL

      END SUBROUTINE EXPAND_PATCHES

C***********************************************************************
C
      SUBROUTINE FINTERP (MXIBODY, MXJBODY, NFSTN, NFJ, IDEGFUS,
     >                    XF, YFINIT, ZFINIT, IKBD, ICRT, ISAVE)
C
C     Insert fuselage defining sections interactively, saving them
C     interleaved with input sections as we go.  The raw geometry must
C     be regular and suited to streamwise lofting.
C
C     12/14/98  David Saunders  Adaptation of program Interp_Body.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   MXIBODY, MXJBODY, NFSTN, NFJ(MXIBODY), IDEGFUS

      REAL, INTENT (IN) ::
     >   XF(MXIBODY)

      REAL, DIMENSION (MXJBODY,MXIBODY), INTENT (IN) ::
     >   YFINIT, ZFINIT

      INTEGER, INTENT (IN) ::
     >   IKBD, ICRT, ISAVE

C     Local constants:

      REAL, PARAMETER ::
     >   ONE = 1., XFLAG = 9.E7

      LOGICAL, PARAMETER ::
     >   NEW = .TRUE.

      CHARACTER, PARAMETER ::
     >   METHOD * 1 = 'B'

C     Local variables:

      REAL, ALLOCATABLE, DIMENSION (:) ::
     >   YNEW, ZNEW

      REAL, ALLOCATABLE, DIMENSION (:,:) ::
     >   Y, Z

      INTEGER
     >   I, IOS, J, L, M, NFP, NNEW

      REAL
     >   DERIV, FNFP, XINT

C     Execution:

      NFP  = NFJ(1)
      FNFP = NFP

C     Check for feasibility:

      DO I = 2, NFSTN
         IF (NFJ(I) /= NFP) THEN
            WRITE (ICRT, '(/, A, 2I5)')
     >         ' FINSERT: Mismatched point count.  I, NFJ(I):',
     >         I, NFJ(I)
            GO TO 99
         END IF
      END DO

      WRITE (ICRT, '(/, A, /, (6F12.4))') ' Current X stations:',
     >   XF(1:NFSTN)
      WRITE (ICRT, '(A)')

      OPEN (UNIT=ISAVE, FILE='aero.sections', STATUS='UNKNOWN')

      ALLOCATE (Y(NFSTN,NFP), Z(NFSTN,NFP), YNEW(NFP), ZNEW(NFP))

      DO I = 1, NFSTN
         Y(I,1:NFP) = YFINIT(1:NFP,I)
         Z(I,1:NFP) = ZFINIT(1:NFP,I)
      END DO

      L = 1 ! Input section pointer
      NNEW = 0

      DO ! Until EOF at the prompt

         WRITE(ICRT, '(A)', ADVANCE='NO') ' X at which to interpolate? '
         READ (IKBD, *, IOSTAT=IOS) XINT

         IF (IOS < 0) XINT = XFLAG ! Force transfer of remaining sections

C        Transfer any untransferred original sections preceding this one:

         DO I = L, NFSTN

            IF (XF(I) < XINT) THEN
               M = I
               WRITE (ISAVE, 100)
               WRITE (ISAVE, 101) FNFP, XF(I), ONE
               WRITE (ISAVE, 102)
               WRITE (ISAVE, 103) (ZFINIT(J,I), YFINIT(J,I), J = 1, NFP)
            ELSE
               L = M + 1
               EXIT
            END IF

         END DO

         IF (XINT == XFLAG) EXIT ! Done

         NNEW = NNEW + 1

         DO J = 1, NFP

            CALL LCSFIT (NFSTN, XF, Y(1,J), NEW, METHOD,
     >                   1, XINT, YNEW(J), DERIV)

            CALL LCSFIT (NFSTN, XF, Z(1,J), NEW, METHOD,
     >                   1, XINT, ZNEW(J), DERIV)
         END DO

         WRITE (ISAVE, 100)
         WRITE (ISAVE, 101) FNFP, XINT, ONE
         WRITE (ISAVE, 102)
         WRITE (ISAVE, 103) (ZNEW(J), YNEW(J), J = 1, NFP)

      END DO ! Next new X

      DEALLOCATE (Y, Z, YNEW, ZNEW)
      CLOSE (UNIT=ISAVE)

      NNEW = NFSTN + NNEW
      WRITE (ICRT, '(/, I4, A)') NNEW,
     >   ' sections written to aero.sections.'

   99 CONTINUE

  100 FORMAT (' FNFP      XF   FSEC')
  101 FORMAT (F5.0, F10.4, F5.0)
  102 FORMAT ('      Z            Y')
  103 FORMAT (2F12.6)

      END SUBROUTINE FINTERP

C***********************************************************************
C
      SUBROUTINE FLAPSLAT (CASE, NL, NU, NTOT, XWG, YWG, ZWG,
     >                     ANGLE, XHINGE, SWEEP)
C
C     Perform a flap or slat rotation on one section.  The numbers of
C     defining points may be changed as a result.
C
C     03/16/99  DAS  Argument-driven portion of DEFLECT from SYN87-SB,
C                    adapted for applying deflections read with geometry.
C     10/03/00   "   Don't assume NL is the leading edge: REORDER can be
C                    fooled if it is not.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      CHARACTER, INTENT (IN) ::
     >   CASE * 1             ! 'F' = flap, 'S' = slat

      INTEGER, INTENT (INOUT) ::
     >   NL, NU, NTOT         ! The usual numbers of airfoil points

      REAL, INTENT (INOUT), DIMENSION (NTOT) ::
     >   XWG, YWG, ZWG        ! Airfoil coordinates

      REAL, INTENT (IN) ::
     >   ANGLE,               ! Deflection about hinge, in degrees; +ve is down
     >   XHINGE,              ! Hinge location in units of XWG
     >   SWEEP                ! Hinge line sweep, in degrees

C     Local constants:

      REAL, PARAMETER ::
     >   DEG2RAD = 0.0174532925, HALF = 0.5, ONE = 1.0, TINY = 0.001

C     Local variables:

      INTEGER
     >   I, J, L, ILE, ILO, IUP, NNEW, INDICES(NTOT)

      REAL
     >   ANGLE_CHORDWISE, COSA, R, SINA, TANGENT, XP, YHINGE,
     >   YLO, YP, YUP

C     Execution:

      TANGENT = TAN (DEG2RAD * ANGLE) * COS (DEG2RAD * SWEEP)
      ANGLE_CHORDWISE = ATAN (TANGENT)
      COSA = COS (ANGLE_CHORDWISE)
      SINA = SIN (ANGLE_CHORDWISE) ! -SINA used for flap

C     Rotate about (XHINGE, Ymean) through the adjusted angle:

      ILO = NL / 2
      CALL INTERVAL (NL, XWG, XHINGE, -ONE, ILO)

      R   = (XHINGE - XWG(ILO+1)) / (XWG(ILO) - XWG(ILO+1))
      YLO = (ONE - R) * YWG(ILO+1) + R * YWG(ILO)

      IUP = NU / 2
      CALL INTERVAL (NU, XWG(NL), XHINGE, ONE, IUP)

      IUP = IUP + NL - 1
      R   = (XHINGE - XWG(IUP)) / (XWG(IUP+1) - XWG(IUP))
      YUP = (ONE - R) * YWG(IUP) + R * YWG(IUP+1)

      YHINGE = (YUP + YLO) * HALF

      IF (CASE == 'F') THEN ! Flap: SINA <-- -SINA makes it clockwise

         DO I = 1, ILO
            XP = XWG(I) - XHINGE
            YP = YWG(I) - YHINGE
            XWG(I) =  XP * COSA + YP * SINA + XHINGE
            YWG(I) = -XP * SINA + YP * COSA + YHINGE
         END DO

         DO I = IUP + 1, NTOT
            XP = XWG(I) - XHINGE
            YP = YWG(I) - YHINGE
            XWG(I) =  XP * COSA + YP * SINA + XHINGE
            YWG(I) = -XP * SINA + YP * COSA + YHINGE
         END DO

      ELSE ! Slat; +ve is anticlockwise as in ROTATE2D

         DO I = ILO + 1, IUP
            XP = XWG(I) - XHINGE
            YP = YWG(I) - YHINGE
            XWG(I) = XP * COSA - YP * SINA + XHINGE
            YWG(I) = XP * SINA + YP * COSA + YHINGE
         END DO

      END IF

C     Suppress any points now inside the wing section.
C     First find the true leading edge, to protect REORDER:

      DO I = NL - 3, NL + 2
         IF (XWG(I+1) > XWG(I)) THEN
            ILE = I
            EXIT
         END IF
      END DO

      IF (ANGLE > TINY) THEN ! Lower surface changes affect upper surface

         CALL REORDER (1, ILE, XWG, INDICES, NNEW)

         IF (NNEW /= ILE) THEN

            DO I = 1, NNEW
               J = INDICES(I)
               XWG(I) = XWG(J)
               YWG(I) = YWG(J)
               ZWG(I) = ZWG(J)
            END DO

            NL   = NNEW
            NU   = NTOT - ILE + 1
            NTOT = NNEW + NU - 1

            DO I = 1, NU - 1
               XWG(I+NNEW) = XWG(I+ILE)
               YWG(I+NNEW) = YWG(I+ILE)
            END DO

         END IF

      ELSE IF (ANGLE < -TINY) THEN ! Negative angle affects upper surface only

         CALL REORDER (ILE, NTOT, XWG, INDICES, NNEW)

         IF (NNEW /= NTOT - ILE + 1) THEN

            NL   = ILE
            NU   = NNEW
            NTOT = NNEW + ILE - 1
            L    = ILE

            DO I = 2, NNEW
               L = L + 1
               J = INDICES(I)
               XWG(L) = XWG(J)
               YWG(L) = YWG(J)
               ZWG(L) = ZWG(J)
            END DO

         END IF

      END IF

      END SUBROUTINE FLAPSLAT

C***********************************************************************
C
      SUBROUTINE FUSEREAD (IGEOM, MXIBODY, MXJBODY, NFSTN, NFJ,
     >                     XF, YFINIT, ZFINIT, IWRIT)
C
C     FUSEREAD reads a fuselage geometry.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   IGEOM, MXIBODY, MXJBODY

      INTEGER, INTENT (OUT) ::
     >   NFSTN, NFJ(MXIBODY)

      REAL, INTENT (OUT) ::
     >   XF(MXIBODY), YFINIT(MXJBODY,MXIBODY), ZFINIT(MXJBODY,MXIBODY)

      INTEGER, INTENT (IN) ::
     >   IWRIT

C     Local constants:

      REAL,      PARAMETER :: ZERO = 0.
      CHARACTER, PARAMETER :: NAME * 10 = 'FUSEREAD: '

C     Local variables:

      INTEGER
     >   I, IOS, J, JJ, N
      REAL
     >   FLIP, FN, FNF, FSEC
      REAL, DIMENSION (MXJBODY) ::
     >   YBUF, ZBUF, YP, ZP

C     Execution.

      READ (IGEOM,*)
      READ (IGEOM,*,IOSTAT=IOS) FNF, FLIP

      IF (IOS /= 0) THEN
         WRITE (IWRIT, '(/, 2A)') NAME, 'Trouble reading FNF, FLIP.'
         GO TO 900
      END IF
      NFSTN = FNF

      IF (NFSTN > MXIBODY) THEN
         WRITE (IWRIT, '(/, 2A, 2I5)') NAME,
     >      'Too many body sections:', NFSTN, MXIBODY
         GO TO 900
      END IF

      DO I = 1, NFSTN

         READ (IGEOM,*)
         IF (FLIP == ZERO) THEN
            READ (IGEOM,*,IOSTAT=IOS) FN, XF(I), FSEC
         ELSE
            READ (IGEOM,*,IOSTAT=IOS) XF(I), FN, FSEC
         END IF

         IF (IOS /= 0) THEN
            WRITE (IWRIT, '(/, 2A, I4)') NAME,
     >         'Trouble reading section #', I
            GO TO 900
         END IF

         N = FN
         IF (FSEC > ZERO) THEN ! New fuselage section

            IF (N > MXJBODY) THEN
               WRITE (IWRIT, '(/, A, I4, A, 2I5)')
     >            ' FUSEREAD: Too many fuselage pts. at station', I,
     >            '.  Input & limit:', N, MXJBODY
               GO TO 900
            END IF

            NFJ(I) = N
            READ (IGEOM,*)

            IF (FLIP == ZERO) THEN
               DO J = 1, N
                  READ (IGEOM,*,IOSTAT=IOS) ZP(J), YP(J)
               END DO
            ELSE
               DO J = 1, N
                  READ (IGEOM,*,IOSTAT=IOS) YP(J), ZP(J)
               END DO
            END IF

            IF (IOS /= 0) THEN
               WRITE (IWRIT, '(/, 2A, I4)') NAME,
     >            'Trouble reading section #', I
               GO TO 900
            END IF


C           Ensure that fuselage points go from keel to crown:

            IF (N > 1) THEN
               IF (YP(1) < YP(N)) THEN
                  YBUF(1:N) = YP(1:N)
                  ZBUF(1:N) = ZP(1:N)
               ELSE
                  DO J = 1, N
                     JJ = N - J + 1
                     YBUF(J) = YP(JJ)
                     ZBUF(J) = ZP(JJ)
                  END DO
               END IF
            ELSE
               YBUF(1) = YP(1)
               ZBUF(1) = ZP(1)
            END IF

         ELSE
            NFJ(I) = NFJ(I-1)
         END IF

         YFINIT(1:N,I) = YBUF(1:N)
         ZFINIT(1:N,I) = ZBUF(1:N)

      END DO

      RETURN

  900 WRITE (IWRIT, '(/, A)') ' Aborting ...'

      STOP

      END SUBROUTINE FUSEREAD

C***********************************************************************
C
      SUBROUTINE FUSEWRITE (LUN, NFSTN)
C
C***********************************************************************

C     Global variables:

      USE GEOM1

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) :: LUN, NFSTN

C     Local constants:

      REAL, PARAMETER :: ONE = 1., ZERO = 0.

C     Local variables:

      INTEGER I, J

      REAL    FNFP

C     Execution:

      WRITE (LUN, 1001) '       FNF      FLIP'
      FNFP = NFSTN
      WRITE (LUN, 1013) FNFP, ZERO

      DO I = 1, NFSTN
         WRITE (LUN, 1001) '        FNFP          XF        FSEC'
         FNFP = NFJ(I)
         WRITE (LUN, 1012) FNFP, XF(I), ONE
         WRITE (LUN, 1001) '        Z(J)        Y(J)'
         WRITE (LUN, 1014) (ZFINIT(J,I), YFINIT(J,I), J = 1, NFJ(I))
      END DO

      RETURN

 1001 FORMAT (A)
 1012 FORMAT (3F12.6)
 1013 FORMAT (2F10.2)
 1014 FORMAT (2F12.6)

      END SUBROUTINE FUSEWRITE

C***********************************************************************
C
      SUBROUTINE GETANGLE (NL, NU, X, Y, THICK, TEANGL)
C
C     Trailing edge angle for wrap-around airfoil (in radians).
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN)  :: NL, NU
      REAL,    INTENT (IN)  :: X(*), Y(*), THICK
      REAL,    INTENT (OUT) :: TEANGL

C     Local constants:

      REAL, PARAMETER :: ARROW = 1., FRACTION = 0.99

C     Local variables:

      INTEGER ILO, IUP, N
      REAL    XTARG

C     Execution:

C     Pick points at ~99% X/C:

      XTARG = X(NL) + (X(1) - X(NL)) * FRACTION
      IUP = NU - 5
      CALL INTERVAL (NU, X(NL), XTARG, ARROW, IUP)
      IUP = IUP + NL - 1

      ILO = 3
      CALL INTERVAL (NL, X(1), XTARG, -ARROW, ILO)
      ILO = ILO + 1

      N = NL + NU - 1
      TEANGL = ATAN (THICK*(Y(1) - Y(ILO)) / (X(1) - X(ILO))) +
     >         ATAN (THICK*(Y(IUP) - Y(N)) / (X(N) - X(IUP)))

      END SUBROUTINE GETANGLE

C***********************************************************************
C
      SUBROUTINE GETTHICK (NL, NU, X, Y, THMAX, XTHMAX)
C
C     Wing section maximum thickness calculation for wrap-around data:
C     interpolate the lower surface at the upper surface Xs, then pick
C     the point of maximum thickness (returned in the input units).
C
C     05/28/96  D.Saunders  Adaptation of familiar stuff.
C     11/27/96      "       Switched to monotonic fits to match SETLCON.
C     ??? Some day:         Do it more thoroughly for gradients.
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN)  :: NL, NU
      REAL,    INTENT (IN)  :: X(*), Y(*)
      REAL,    INTENT (OUT) :: THMAX, XTHMAX

C     Local variables:

      INTEGER I
      REAL    THICK, YEVAL(NL+NU)

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
      SUBROUTINE KAPPA (NNKAPPA, NNFORMAT, ICRT, IKBD, ISAVE,
     >                  IDIM, KDIM, ILWING, NUMK, NFSTN, JLBODY,
     >                  XRG, YRG, ZRG, XF, YF, ZF)
C
C     Evaluation of curvature or radius of curvature in the index directions
C     for the indicated regularized geometry component.
C
C     12/18/98  D. Saunders  Initial implementation (curvature).
C     12/21/98     "         Allowed for curvature or radius.
C     01/20/99     "         Combined index directions in one FAST file.
C     01/21/99     "         Went back to two scalar files; introduced
C                            X/Y/ZCURV because X/Y/ZBASE have wrong
C                            dimensions for the fuselage.
C     01/30/99     "         Regularized wings/nacelles are packed now.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   NNKAPPA,     ! +n means process wing surface n; -1 means fuselage
     >   NNFORMAT,    ! 0 or 2 means unformatted outputs; 1 or 3 means formatted
     >   ICRT,        ! Logical unit for screen output ...
     >   IKBD,        ! ... and for keyboard input ...
     >   ISAVE,       ! ... and for output curvature and corresp. (x,y,z)s
     >   IDIM,        ! Array sizes for THIS regularized wing-type surface
     >   KDIM,
     >   ILWING,      ! Active dimensions
     >   NUMK         ! Active dimensions

      INTEGER, INTENT (IN) ::
     >   NFSTN,       ! Array and active dimensions for regularized fuselage
     >   JLBODY

      REAL, INTENT (INOUT), DIMENSION (IDIM, KDIM) ::
     >   XRG,         ! Regularized wing-type surfaces
     >   YRG, ZRG

      REAL, INTENT (IN), DIMENSION (NFSTN) ::
     >   XF           ! Fuselage stations Xs

      REAL, INTENT (IN), DIMENSION (JLBODY, NFSTN) ::
     >   YF, ZF       ! Regularized body sections

C     Local constants:

      REAL, PARAMETER ::
     >   HUGE = 1.E7, ONE = 1., TINY = 1.E-7, ZERO = 0.
      LOGICAL, PARAMETER ::
     >   NORM = .FALSE.

C     Local variables:

      INTEGER
     >   I, J, J1, NCUTS, NPTS
      REAL
     >   TOTAL
      REAL, ALLOCATABLE, DIMENSION (:) ::
     >   ARC, X, XS, XSS, Y, YS, YSS, Z, ZS, ZSS
      REAL, ALLOCATABLE, DIMENSION (:,:) ::
     >   CURV, XCURV, YCURV, ZCURV
      LOGICAL
     >   RADIUS, WING
      CHARACTER
     >   RESPONSE * 1

C     Execution:

      WRITE (ICRT, '(/, A)', ADVANCE='NO')
     >   ' Curvature (c) or Radius of Curvature (r)? '
      READ (IKBD, '(A)') RESPONSE
      RADIUS = RESPONSE == 'r' .OR. RESPONSE == 'R'
      WING = NNKAPPA > 0

      IF (WING) THEN
         J1 = 1
         NPTS  = ILWING
         NCUTS = NUMK
      ELSE ! NNKAPPA = -1 for the fuselage
         J1 = 2 ! Suppress the singular nose section
         NPTS  = JLBODY - 1 ! Pitfall - regularization allows for intersection
         NCUTS = NFSTN
      END IF

C     Curvature along 2-space cuts can have a useful sign attached:

      ALLOCATE (ARC(NPTS), CURV(NPTS,NCUTS), XCURV(NPTS,NCUTS),
     >          YCURV(NPTS,NCUTS), ZCURV(NPTS,NCUTS),
     >          XS(NPTS), XSS(NPTS), YS(NPTS), YSS(NPTS))

      IF (WING) THEN

         XCURV = XRG(1:NPTS,1:NCUTS)
         YCURV = YRG(1:NPTS,1:NCUTS)
         ZCURV = ZRG(1:NPTS,1:NCUTS)

      ELSE ! Order body (x,y,z)s to match those for 2D wing section curvatures

         DO J = 1, NCUTS
            XCURV(1:NPTS,J) = YF(1:NPTS,J)
            YCURV(1:NPTS,J) = ZF(1:NPTS,J)
            ZCURV(1:NPTS,J) = XF(J)
         END DO

         IF (RADIUS) THEN ! For the nose point
            CURV(1:NPTS,1) = HUGE
         ELSE
            CURV(1:NPTS,1) = ZERO
         END IF

      END IF

C     Evaluate curvature along all the cuts:

      DO J = J1, NCUTS

C        Arc lengths along this cut:

         CALL CHORDS2D (NPTS, XCURV(1,J), YCURV(1,J), NORM, TOTAL, ARC)

C        3-point forward difference dX/dS and dY/dS at point 1:

         CALL FD1SID (1, 1, ARC, XCURV(1,J), XS(1), XSS(1))
         CALL FD1SID (1, 1, ARC, YCURV(1,J), YS(1), YSS(1))

         DO I = 2, NPTS - 1 ! Central differences at interior points

            CALL FDCNTR (I, ARC, XCURV(1,J), XS(I), XSS(I))
            CALL FDCNTR (I, ARC, YCURV(1,J), YS(I), YSS(I))

         END DO

C        3-point backward difference dX/dS and dY/dS at last point:

         I = NPTS
         CALL FD1SID (I, -1, ARC, XCURV(1,J), XS(I), XSS(I))
         CALL FD1SID (I, -1, ARC, YCURV(1,J), YS(I), YSS(I))

C        Curvature along the cut:

         CALL CURV2D (NPTS, XS, XSS, YS, YSS, CURV(1,J))

      END DO ! Next cut

      IF (RADIUS) THEN

         DO J = 1, NCUTS
            DO I = 1, NPTS
               CURV(I,J) = ONE / MAX (ABS (CURV(I,J)), TINY)
            END DO
         END DO

      END IF

C     Save results:

      IF (.NOT. WING) THEN ! Recover the fuselage (x,y,z) order

         DO J = 1, NCUTS
            XCURV(1:NPTS,J) = XF(J)
            YCURV(1:NPTS,J) = YF(1:NPTS,J)
            ZCURV(1:NPTS,J) = ZF(1:NPTS,J)
         END DO

      END IF

      IF (MOD (NNFORMAT, 2) == 0) THEN ! Unformatted

         OPEN (UNIT=ISAVE, FILE='aero.ikappa', STATUS='UNKNOWN',
     >         FORM='UNFORMATTED')

         WRITE (ISAVE) 1
         WRITE (ISAVE) NPTS, NCUTS, 1, 1
         WRITE (ISAVE) CURV

         CLOSE (ISAVE)

         OPEN (UNIT=ISAVE, FILE='aero.regular', STATUS='UNKNOWN',
     >         FORM='UNFORMATTED')

         WRITE (ISAVE) 1
         WRITE (ISAVE) NPTS, NCUTS, 1
         WRITE (ISAVE) XCURV, YCURV, ZCURV

      ELSE

         OPEN (UNIT=ISAVE, FILE='aero.ikappa', STATUS='UNKNOWN')

         WRITE (ISAVE, '(I4)') 1
         WRITE (ISAVE, '(4I4)') NPTS, NCUTS, 1, 1
         WRITE (ISAVE, '(1P, 10E13.5)') CURV

         CLOSE (ISAVE)

         OPEN (UNIT=ISAVE, FILE='aero.regular', STATUS='UNKNOWN')

         WRITE (ISAVE, '(I4)') 1
         WRITE (ISAVE, '(3I4)') NPTS, NCUTS, 1
         WRITE (ISAVE, '(6F12.6)') XCURV, YCURV, ZCURV

      END IF

      CLOSE (ISAVE)

      DEALLOCATE (ARC, CURV, XS, XSS, YS, YSS)


      IF (NCUTS <= 2) THEN
         WRITE (ICRT, '(//, A)')
     >      ' Spanwise calculations require at least 3 sections.'
         GO TO 800
      END IF

C     Curvature along spanwise lines (wing) or streamwise lines (fuselage):

      ALLOCATE (ARC(NCUTS), CURV(NCUTS,NPTS),
     >          X(NCUTS), XS(NCUTS), XSS(NCUTS),
     >          Y(NCUTS), YS(NCUTS), YSS(NCUTS),
     >          Z(NCUTS), ZS(NCUTS), ZSS(NCUTS))

      DO I = 1, NPTS

         X = XCURV(I,1:NCUTS)
         Y = YCURV(I,1:NCUTS)
         Z = ZCURV(I,1:NCUTS)

C        Arc lengths across cuts for this point:

         CALL CHORDS3D (NCUTS, X, Y, Z, NORM, TOTAL, ARC)

C        3-point forward differences at point 1:

         CALL FD1SID (1, 1, ARC, X, XS(1), XSS(1))
         CALL FD1SID (1, 1, ARC, Y, YS(1), YSS(1))
         CALL FD1SID (1, 1, ARC, Z, ZS(1), ZSS(1))

         DO J = 2, NCUTS - 1 ! Central differences at interior points

            CALL FDCNTR (J, ARC, X, XS(J), XSS(J))
            CALL FDCNTR (J, ARC, Y, YS(J), YSS(J))
            CALL FDCNTR (J, ARC, Z, ZS(J), ZSS(J))

         END DO

C        3-point backward differences at last point:

         J = NCUTS
         CALL FD1SID (J, -1, ARC, X, XS(J), XSS(J))
         CALL FD1SID (J, -1, ARC, Y, YS(J), YSS(J))
         CALL FD1SID (J, -1, ARC, Z, ZS(J), ZSS(J))

C        Curvature magnitudes along the line:

         CALL CURV3D (NCUTS, XS, XSS, YS, YSS, ZS, ZSS, CURV(1,I))

      END DO ! Next line across cuts

      IF (RADIUS) THEN

         DO I = 1, NPTS
            DO J = 1, NCUTS
               CURV(J,I) = ONE / MAX (ABS (CURV(J,I)), TINY)
            END DO
         END DO

      END IF

C     Save results:

      IF (MOD (NNFORMAT, 2) == 0) THEN

         OPEN (UNIT=ISAVE, FILE='aero.jkappa', STATUS='UNKNOWN',
     >         FORM='UNFORMATTED')

         WRITE (ISAVE) 1
         WRITE (ISAVE) NPTS, NCUTS, 1, 1
         WRITE (ISAVE) ((CURV(J,I), I = 1, NPTS), J = 1, NCUTS)

      ELSE

         OPEN (UNIT=ISAVE, FILE='aero.jkappa', STATUS='UNKNOWN')

         WRITE (ISAVE, '(I4)') 1
         WRITE (ISAVE, '(4I4)') NPTS, NCUTS, 1, 1
         WRITE (ISAVE, '(1P, 10E13.5)')
     >         ((CURV(J,I), I = 1, NPTS), J = 1, NCUTS)
      END IF

      CLOSE (ISAVE)

      DEALLOCATE (ARC, CURV, XCURV, YCURV, ZCURV, X, XS, XSS,
     >            Y, YS, YSS, Z, ZS, ZSS)

  800 IF (RADIUS) THEN
         WRITE (ICRT, '(/, (A))')
     >      ' Surface radius of curvature evaluated as follows:',
     >      ' aero.ikappa  = rho along each wing or body section',
     >      ' aero.jkappa  = rho along spanwise or body axial lines'
      ELSE
         WRITE (ICRT, '(/, (A))')
     >      ' Surface curvature evaluated as follows:',
     >      ' aero.ikappa  = +/-kappa along each wing or body section',
     >      ' aero.jkappa  = |kappa| along spanwise or body axial lines'
      END IF

      WRITE (ICRT, '(A)')
     >   ' aero.regular = corresponding regular wing or body (x,y,z)s'

      END SUBROUTINE KAPPA

C***********************************************************************
C
      SUBROUTINE LOFT (TYPE2, Z, ZCEN, ZMIN, ZMAX, PWR, HEIGHT, LUNERR)
C
C     LOFT modularizes the choice between polynomial and Hicks-Henne-type
C     spanwise lofting of the effect of a wing design variable, for use
C     by both PERTURB and SETLCON.  HEIGHT is returned between 0. and 1.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      CHARACTER, INTENT (IN)  :: TYPE2 * 4
      INTEGER,   INTENT (IN)  :: LUNERR
      REAL,      INTENT (IN)  :: Z, ZCEN, ZMIN, ZMAX, PWR
      REAL,      INTENT (OUT) :: HEIGHT

C     Constants:

      REAL, PARAMETER :: ONE = 1., ZERO = 0.

C     Local variables:

      REAL :: ZC

C     Execution:

      IF (TYPE2 == 'POLY') THEN

C        Polynomial-type lofting (often linear or constant):
C        The constant PWR = 0. case requires special treatment at end-pts.

         CALL RIPPER (Z, ZCEN, ZMIN, ZMAX, PWR, HEIGHT)

      ELSE ! Hicks-Henne-type lofting:

         HEIGHT = ZERO ! SFEVAL adds

         IF (Z >= ZMIN .AND. Z <= ZMAX) THEN ! Else clip

            IF (INDEX (TYPE2, 'COS') > 0) THEN ! Generalized COS* and *COS
               ZC =  ZCEN
            ELSE
               ZC = (ZCEN - ZMIN) / (ZMAX - ZMIN) ! Reqd. by SFEVAL
            END IF

            CALL SFEVAL (TYPE2, PWR, ZC, ZMIN, ZMAX,
     >                   ONE, 1, 1, Z, HEIGHT, LUNERR)
         END IF

      END IF

      END SUBROUTINE LOFT

C***********************************************************************
C
      SUBROUTINE NCYLINDER (MXIGRID, MXKGRID, NWSEC, IL, KL,
     >                      XRG, YRG, ZRG)
C
C     NCYLINDER converts nacelle sections from cylindrical coordinates to
C     (X,Y,Z)s in the nacelle axis system, and interpolates them to a
C     uniform distribution in the azimuthal direction.
C     ZRG(*,*) should include any clock angle (nacelle roll angle about
C     its axis) on input.  X/Y/ZRG are overwritten on output.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   MXIGRID, MXKGRID, NWSEC, IL, KL

      REAL, INTENT (INOUT), DIMENSION (MXIGRID,MXKGRID) ::
     >   XRG, YRG, ZRG ! Regularized Xs, Rs, Thetas on input;
                       ! redistributed (X,Y,Z)s on output (still axis-centered)

C     Local constants:

      REAL, PARAMETER :: ONE = 1.

C     Local variables:

      INTEGER
     >   I, K, KF, KK
      REAL
     >   CP, DPHI, PHIK, PHIL, R, RM1, RNEW, SP, TWOPI
      REAL, DIMENSION (IL,KL) ::
     >   XNAC, YNAC, ZNAC

C     Execution:

      TWOPI = 4.* ASIN (ONE)
      DPHI  = TWOPI / REAL (KL - 1)

      IF (NWSEC == 1) THEN

         DO K = 1, KL
            PHIK = DPHI * REAL (K - 1) + ZRG(1,1)
            CP   = COS (PHIK)
            SP   = SIN (PHIK)
            DO I = 1, IL
               XNAC(I,K) = XRG(I,1)
               YNAC(I,K) = YRG(I,1) * CP
               ZNAC(I,K) = YRG(I,1) * SP
            END DO
         END DO

      ELSE

         PHIL = ZRG(1,NWSEC)
         KF   = 1

         DO K = 1, KL
            PHIK = DPHI * REAL (K - 1) + ZRG(1,1)
            CP   = COS (PHIK)
            SP   = SIN (PHIK)

            IF (PHIK > PHIL) THEN

               R   = (PHIK - PHIL) / (TWOPI - PHIL)
               RM1 = ONE - R

               DO I = 1, IL
                  RNEW      = RM1 * YRG(I,NWSEC) + R * YRG(I,1)
                  XNAC(I,K) = RM1 * XRG(I,NWSEC) + R * XRG(I,1)
                  YNAC(I,K) = RNEW * CP
                  ZNAC(I,K) = RNEW * SP
               END DO

            ELSE

               DO KK = KF, NWSEC - 1

                  IF (PHIK >= ZRG(1,KK) .AND.
     >                PHIK <= ZRG(1,KK+1)) THEN
                     KF  = KK
                     R   = (PHIK - ZRG(1,KK)) / (ZRG(1,KK+1)-ZRG(1,KK))
                     RM1 = ONE - R
                     DO I = 1, IL
                        RNEW      = RM1 * YRG(I,KK) + R * YRG(I,KK+1)
                        XNAC(I,K) = RM1 * XRG(I,KK) + R * XRG(I,KK+1)
                        YNAC(I,K) = RNEW * CP
                        ZNAC(I,K) = RNEW * SP
                     END DO
                     EXIT
                  END IF

               END DO

            END IF

         END DO ! Next uniform phi

      END IF

      XRG(1:IL,1:KL) = XNAC ! Overwrite the inputs
      YRG(1:IL,1:KL) = YNAC
      ZRG(1:IL,1:KL) = ZNAC

      END SUBROUTINE NCYLINDER

C***********************************************************************
C
      SUBROUTINE NPITCHYAW (MXIGRID, MXKGRID, IL, KL, XC1, YC1, ZC1,
     >                      XC2, YC2, ZC2, XNAC, YNAC, ZNAC)
C
C     NPITCHYAW transforms one regularized, padded nacelle surface from
C     axes centered on the nacelle axis to airframe coordinates, in-place.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   MXIGRID, MXKGRID, IL, KL

      REAL, INTENT (IN) ::
     >   XC1, YC1, ZC1, XC2, YC2, ZC2

      REAL, INTENT (INOUT), DIMENSION (MXIGRID,MXKGRID) ::
     >   XNAC, YNAC, ZNAC

C     Local variables:

      INTEGER
     >   I, K
      REAL
     >   CPCH, CYAW, PTCH, SPCH, SYAW, YAW, XT, YT, ZT

C     Execution:


C     Transform coordinates:

      YAW  = ATAN ((ZC2 - ZC1) / (XC2 - XC1))
      SYAW = SIN (YAW)
      CYAW = COS (YAW)
      PTCH =       (YC2 - YC1) / (XC2 - XC1)
      PTCH = ATAN (PTCH * CYAW)
      SPCH = SIN (PTCH)
      CPCH = COS (PTCH)

      DO K = 1, KL
         DO I = 1, IL
            XT = XNAC(I,K)
            YT = YNAC(I,K)
            ZT = ZNAC(I,K)
            XNAC(I,K) = XC1 + CYAW*(CPCH*XT - SPCH*YT) - SYAW*ZT
            YNAC(I,K) = YC1 +       SPCH*XT + CPCH*YT
            ZNAC(I,K) = ZC1 + SYAW*(CPCH*XT - SPCH*YT) + CYAW*ZT
         END DO
      END DO

      END SUBROUTINE NPITCHYAW

C***********************************************************************
C
      SUBROUTINE PATCH_AREA (NI, NJ, X, Y, Z, XMIN, XLENGTH, AREA_PATCH)
C
C     Area and streamwise length for a packed surface patch.
C
C     09/12/00  David Saunders  Introduced for viscous drag estimates.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   NI, NJ                              ! Packed patch dimensions
      REAL, DIMENSION (NI, NJ), INTENT (IN) ::
     >   X, Y, Z                             ! Coordinates
      REAL, INTENT (OUT) ::
     >   XMIN, XLENGTH,                      ! Mean Xmin and Xlength over
     >   AREA_PATCH                          ! all Js, and surface area

C     Local constants:

      REAL, PARAMETER ::
     >   ZERO = 0.

C     Local variables:

      INTEGER
     >   I, J
      REAL
     >   XLEFT, XLMEAN, XRIGHT, XRMEAN

C     Procedure:

      REAL, EXTERNAL ::
     >   AREA4        ! Area of a quadrilateral

C     Execution:

      XLMEAN     = ZERO
      XRMEAN     = ZERO
      AREA_PATCH = ZERO

      DO J = 1, NJ

         XLEFT  = X(1,J)
         XRIGHT = XLEFT

         DO I = 2, NI
            XLEFT  = MIN (XLEFT,  X(I,J))
            XRIGHT = MAX (XRIGHT, X(I,J))
         END DO

         XLMEAN = XLMEAN + XLEFT
         XRMEAN = XRMEAN + XRIGHT

         IF (J == 1) CYCLE

         DO I = 2, NI
            AREA_PATCH = AREA_PATCH + AREA4
     >        (X(I-1,J-1), Y(I-1,J-1), Z(I-1,J-1),
     >         X(I,J-1),   Y(I,J-1),   Z(I,J-1),
     >         X(I,J),     Y(I,J),     Z(I,J),
     >         X(I-1,J),   Y(I-1,J),   Z(I-1,J))
         END DO

      END DO

      XRMEAN  = XRMEAN / REAL (NJ)
      XLMEAN  = XLMEAN / REAL (NJ)
      XMIN    = XLMEAN
      XLENGTH = XRMEAN - XLMEAN

      END SUBROUTINE PATCH_AREA

C***********************************************************************
C
      SUBROUTINE PERTURB
     >  (NCALL, NWSURF, MXIWING, MXKWING, MXIBODY, MXJBODY, MXIGRID,
     >   NNCONFIG, NNWING, NNNAC, NNFUSE, LENRG, NWSEC, IDIM, KDIM, IRG,
     >   CLOCK, ALG, ROLL, ROUNDED, XLEADI, YLEADI, ZLEADI,
     >   DIHED, DIAXIS, DFLAP, DSLAT, XFLAP, XSLAT, SWFLAP, SWSLAT,
     >   AINIT, CINIT, TINIT, XINIT, YINIT, ZINIT,
     >   NLGI, NUGI, NWGI, NLG, NUG, NWG, ILWING, KLWING, XWG, YWG, ZWG,
     >   XC1, YC1, ZC1, XC2, YC2, ZC2,
     >   SNORM, CRDSPL, CRDSPQ, CRDSPS, CRDSPC, XRG, YRG, ZRG,
     >   NFSTN, NFJ, JLBODY, IDEGFUS, XF, YFINIT, ZFINIT, YF, ZF,
     >   XTRAP, ZINNER, ZOUTER, NNROOT, NNFUDGE, NNREPEAT,
     >   IBOD1, IBOD2, IBOD3, KWING2,
     >   NSHELF, XSHELF, YSHELF, ZSHELF, IWRIT)
C
C     PERTURB applies any shape functions to the input geometry, denormalizes,
C     and perhaps adjusts fuselage sections to protect the wing/body
C     intersection calculation.  It returns sections in regularized form.
C
C     Oct 1997     DAS  Initial pared-down PERTURB (no shape functions).
C     09/23/98      "   Wing (and body) regularization is done in PERTURB
C                       now as part of controlling diverter height.
C     10/22/98      "   Nacelle transformations moved here from SURFER.
C     10/26/98      "   Introduced WNFUDGE, requiring NNROOT.
C     01/30/99      "   Regularized wings are packed in X/Y/ZRG(*).
C     8/18-8/26/99  "   Retrofitted all the shape functions.
C     06/23/00      "   Introduced DIHED shape function.
C     07/19/00      "   Introduced CONFIG for CTV case (NJBODY = JLBODY).
C
C***********************************************************************

      USE SHAPE_FUNCTIONS

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   NCALL, NWSURF, MXIWING, MXKWING, MXIBODY, MXJBODY,
     >   MXIGRID, NNCONFIG, NNWING, NNNAC, NNFUSE, LENRG

      INTEGER, INTENT (IN), DIMENSION (NWSURF) ::
     >   NWSEC, IDIM, KDIM, IRG

      REAL, INTENT (IN) ::
     >   CLOCK(NWSURF)

      REAL, INTENT (IN), DIMENSION (MXKWING,NWSURF) ::
     >   ROLL, ROUNDED

      REAL, INTENT (INOUT), DIMENSION (MXKWING,NWSURF) ::
     >   ALG, XLEADI, YLEADI, ZLEADI, DIHED

      REAL, INTENT (INOUT), DIMENSION (3,2,MXKWING,NWSURF) ::
     >   DIAXIS

      REAL, INTENT (INOUT), DIMENSION (MXKWING,NWSURF) ::
     >   DFLAP, DSLAT, XFLAP, XSLAT, SWFLAP, SWSLAT,
     >   AINIT, CINIT, TINIT

      REAL, INTENT (INOUT), DIMENSION (MXIWING,MXKWING,NWSURF) ::
     >   XINIT, YINIT, ZINIT

      INTEGER, INTENT (INOUT), DIMENSION (MXKWING,NWSURF) ::
     >   NLGI, NUGI, NWGI, NLG, NUG, NWG

      INTEGER, INTENT (IN), DIMENSION (NWSURF) ::
     >   ILWING, KLWING

      REAL, INTENT (OUT), DIMENSION (MXIWING,MXKWING,NWSURF) ::
     >   XWG, YWG, ZWG

      REAL, INTENT (IN), DIMENSION (NWSURF) ::
     >   XC1, YC1, ZC1, XC2, YC2, ZC2

      REAL, INTENT (OUT), DIMENSION (MXIGRID,NWSURF) ::
     >   SNORM

      REAL, INTENT (IN), DIMENSION (NWSURF) ::
     >   CRDSPL, CRDSPQ, CRDSPS, CRDSPC

      REAL, INTENT (OUT), DIMENSION (LENRG) ::
     >   XRG, YRG, ZRG

      INTEGER, INTENT (IN) ::
     >   NFSTN, NFJ(MXIBODY), JLBODY, IDEGFUS

      REAL, INTENT (IN) ::
     >   XF(MXIBODY)

      REAL, INTENT (INOUT), DIMENSION (MXJBODY,MXIBODY) ::
     >   YFINIT, ZFINIT

      REAL, INTENT (OUT), DIMENSION (JLBODY,NFSTN) ::
     >   YF, ZF ! Regularized body sections

      REAL, INTENT (IN), DIMENSION (NWSURF) ::
     >   XTRAP, ZINNER, ZOUTER

      INTEGER, INTENT (IN) ::
     >   NNROOT(NWSURF), NNFUDGE, NNREPEAT, IBOD1, IBOD2, IBOD3,
     >   KWING2(NWSURF)

      INTEGER, INTENT (IN) ::
     >   NSHELF

      REAL, INTENT (INOUT), DIMENSION (NSHELF) ::
     >   XSHELF                              ! Denormalized in WBFUDGE

      REAL, INTENT (IN), DIMENSION (NSHELF) ::
     >   YSHELF, ZSHELF

      INTEGER, INTENT (IN) ::
     >   IWRIT

C     Local constants:

      REAL, PARAMETER ::
     >   HALF = 0.5, ONE = 1., TOLER = 1.E-6, TWO = 2., ZERO = 0.

      LOGICAL, PARAMETER ::
     >   TRUE = .TRUE.

C     Local variables:

      INTEGER
     >   I, I1, I1W, ID, IF1, IL, ITER, IWING, J, IW, K, KC, KI, KL, KO,
     >   N, NFP, NHALF, NI, NJ, NJBODY, NK, NLOWK, NTOTK, NUPPK, NWVAR

      REAL
     >   AAX, ALPHA, ALPHA1, ALPHA2, AMULT, ANGLE, AREATOL, BBX, BUMP1,
     >   CAB, CCX, CENTR, COF1, COF2, COF3, COF4, DELTAV, DIFA, DIFR,
     >   DL, DS, DU, DZMAX, RLW, RXBLONG, RXBRANGE, SAB, SCALE,
     >   TAGA, TAGB, TAGR, THICK, TLW, TOTA, TOTB, TOTR, UMAX, UMIN, VS,
     >   WIDTH, XA, XB, XBLONG, XBRANGE, XCENTER, XDD, XHK, XLW, XNORM,
     >   XXW, YCAM, YKEEL, YLW, YYW, ZLW

      REAL, SAVE ::
     >   DEG2RAD, PI, TWOPI

      REAL, DIMENSION (JLBODY) ::
     >   DERIVS, SUNIF

      REAL, DIMENSION (MXJBODY) ::
     >   COSTHETA, RBUMP, SFUS, SINTHETA, YNORM, ZNORM

      LOGICAL
     >   CONSTANT_XOVERC, LE_ONCE, NEWBODY

      CHARACTER
     >   TYPE * 6, TYPE2 * 4, TYP1 * 1, TYP3 * 3, TYP4 * 4

      CHARACTER, SAVE ::
     >   METHOD * 1

C     Procedures:

      EXTERNAL
     >   CHORDS2D, DEFLECT, FLAPSLAT, FOILGRD, LCSFIT, LOFT, NCYLINDER,
     >   NPITCHYAW, ROTATE2D, ROTATE3D, SFEVAL, WBFUDGE, WNFUDGE, WREGUL


C     Execution:
C     ----------

      IF (NCALL == -3) THEN ! First call

         PI      = 4.* ATAN (ONE)
         TWOPI   = TWO * PI
         DEG2RAD = PI / 180.

C        Translate any INTEGER lofting variables to wing section Zs:

         DO N = 1, NTWING                    ! Includes planform vars. & flaps
            IF (KCEN(N) /= 0) THEN
               IW = IDWING(N)
               UBUMP(N) = ZLEADI(KCEN(N),IW) ! Z at shape function peak,
               UBMIN(N) = ZLEADI(KMIN(N),IW) ! etc.
               UBMAX(N) = ZLEADI(KMAX(N),IW)
            END IF
         END DO

      ELSE ! Flap deflections may have changed pt. counts

         DO IW = 1, NWSURF
            DO K = 1, NWSEC(IW)
               NLG(K,IW) = NLGI(K,IW)
               NUG(K,IW) = NUGI(K,IW)
               NWG(K,IW) = NWGI(K,IW)
            END DO
         END DO

      END IF

C     All geometry variables are changed in-place, since no more than
C     one set of perturbations is intended in the AEROSURF application.

C     --------------------------
C     Planform design variables:
C     --------------------------

      DO N = 1, NCHCKS

         VS   = V(N) * VSCALE(N)
         IF (VS == ZERO) CYCLE ! Next N

         TYP3 = VTYPE(N)(1:3)
         XHK  = XBUMP(N) ! Fixed for 'TWT' but not for 'TWT2'
         IW   = IDWING(N)

         DO K = 1, NWSEC(IW)

            CALL LOFT (UTYPE(N), ZLEADI(K,IW), UBUMP(N), UBMIN(N),
     >                 UBMAX(N), PBUMP(N), DELTAV, IWRIT) ! DELTAV is in [0.,1.]

CCC         write (iwrit, '(a,i3,1x,a,1x,6f12.6)')
CCC  >         ' K,UT,Z,UB,UMN,UMX,PB,DV:',
CCC  >         k, UTYPE(N), ZLEADI(K,IW), UBUMP(N), UBMIN(N),
CCC  >         UBMAX(N), PBUMP(N), DELTAV

            DELTAV = DELTAV * VS

            IF (DELTAV /= ZERO) THEN

               IF (TYP3 == 'TWT') THEN

                  IF (VTYPE(N)(4:4) == '2') THEN ! 'TWT2' has arbitrary axis

                     CALL XHINGE_AT_THIS_K () ! Internal procedure sets XHK X/C

                  END IF ! Else 'TWT' = twist about constant-X/C axis

                  ALPHA1       = AINIT(K,IW)  * DEG2RAD
                  AINIT(K,IW)  = AINIT(K,IW)  + DELTAV
                  ALPHA2       = AINIT(K,IW)  * DEG2RAD

                  XLEADI(K,IW) = XLEADI(K,IW) + XHK * CINIT(K,IW) *
     >                          (COS (ALPHA1) - COS (ALPHA2))

                  YLEADI(K,IW) = YLEADI(K,IW) + XHK * CINIT(K,IW) *
     >                          (SIN (ALPHA2) - SIN (ALPHA1))

               ELSE IF (TYP3 == 'THK') THEN
                  TINIT(K,IW)  = TINIT(K,IW)  + DELTAV
               ELSE IF (TYP3 == 'XLD') THEN
                  XLEADI(K,IW) = XLEADI(K,IW) + DELTAV
               ELSE IF (TYP3 == 'YLD') THEN
                  YLEADI(K,IW) = YLEADI(K,IW) + DELTAV
               ELSE IF (TYP3 == 'ZLD') THEN
                  ZLEADI(K,IW) = ZLEADI(K,IW) + DELTAV
               ELSE IF (TYP3 == 'CRD') THEN
                  CINIT(K,IW)  = CINIT(K,IW)  + DELTAV
               END IF

            END IF

         END DO ! Next wing or nacelle section

      END DO ! Next planform design variable


C     --------------------------------------
C     Wing or nacelle section perturbations:
C     --------------------------------------

C     Check for improperly normalized sections.  (Is this really needed?)

      DO IW = 1, NNWING + NNNAC

         DO K = 1, NWSEC(IW)

            NLOWK = NLG(K,IW)

            IF (XINIT(NLOWK,K,IW) < ZERO) THEN
               WRITE (IWRIT, '(/, A, 1P, E13.5, A, I4)')
     >            ' Leading edge X < 0:', XINIT(NLOWK,K,IW), '  K:', K
               GO TO 990
            END IF

            NTOTK = NWG(K,IW)

            IF (XINIT(1,K,IW) > ONE .OR.
     >          XINIT(NTOTK,K,IW) > ONE) THEN
               WRITE (IWRIT, '(/, A, 1P, 2E13.5, A, I4)')
     >            ' Lower or upper TE X > 1:',
     >            XINIT(1,K,IW), XINIT(NTOTK,K,IW), '  K:', K
               GO TO 990
            END IF

         END DO

      END DO


C     Apply any standard Hicks-Henne shape functions to wing/nacelle sections:

      NWVAR = NCHCK1 + NCHCK2 ! Only 1 can be nonzero; NSINVC = 0 if NCHCK2 /= 0
                              ! Perturbing Y" is not fully implemented here

      DO N = NCHCKS + 1, NCHCKS + NWVAR + NSINVC

         VS = V(N) * VSCALE(N)
         IF (VS == ZERO) CYCLE ! Next N

         CONSTANT_XOVERC = N <= NCHCKS + NWVAR

         TYP4 = VTYPE(N)(1:4)
         LE_ONCE = TYP4(1:2) == 'LE' .OR. TYP4(1:2) == 'DR'

         KI = KMIN(N) ! Inboard lofting station used by varying-X/C functions
         KO = KMAX(N) ! Outboard ...
         IW = IDWING(N)

         DO K = 1, NWSEC(IW)

            CALL LOFT (UTYPE(N), ZLEADI(K,IW), UBUMP(N), UBMIN(N),
     >                 UBMAX(N), PBUMP(N), DELTAV, IWRIT) ! DELTAV is in [0.,1.]

CCC         write (iwrit, '(a,i3,1x,a,1x,6f12.6)')
CCC  >         ' K,UT,Z,UB,UMN,UMX,PB,DV:',
CCC  >         k, UTYPE(N), ZLEADI(K,IW), UBUMP(N), UBMIN(N),
CCC  >         UBMAX(N), PBUMP(N), DELTAV

            DELTAV = DELTAV * VS

            IF (DELTAV /= ZERO) THEN

               NTOTK = NWG(K,IW)
               NLOWK = NLG(K,IW)
               NUPPK = NUG(K,IW)

               IF (CONSTANT_XOVERC) THEN

C                 Some variables will perturb the leading edge twice if
C                 we don't perturb both surfaces in one call, so:

                  IF (LE_ONCE) THEN

                     CALL SFEVAL (TYP4, WBUMP(N), XBUMP(N),
     >                            XWMIN(N), XWMAX(N), DELTAV, 1, NTOTK,
     >                            XINIT(1,K,IW), YINIT(1,K,IW), IWRIT)
                  ELSE

                     IF (DOLO(N) /= ZERO) THEN

                        DL = DELTAV * SIGN (ONE, DOLO(N))

                        CALL SFEVAL (TYP4, WBUMP(N), XBUMP(N),
     >                               XWMIN(N), XWMAX(N), DL, 1, NLOWK,
     >                               XINIT(1,K,IW), YINIT(1,K,IW),IWRIT)
                     END IF

                     IF (DOUP(N) /= ZERO) THEN

                        DU = DELTAV * SIGN (ONE, DOUP(N))

                        CALL SFEVAL (TYP4, WBUMP(N), XBUMP(N),
     >                               XWMIN(N), XWMAX(N), DU,NLOWK,NTOTK,
     >                               XINIT(1,K,IW), YINIT(1,K,IW),IWRIT)
                     END IF

                  END IF

               ELSE ! SIN*VC bump (varying center)

C                 Determine the hinge X/C (or X/C of line parallel to shock)
C                 at this span station:

                  CALL XHINGE_AT_THIS_K () ! Internal procedure sets XHK

                  IF (XBUMP(N) > ZERO) THEN

                     XCENTER = XBUMP(N) ! Refers to chord of slat/flap/bump

                     IF (XWMAX(N) < HALF) THEN ! Slat bump
                        XA = ZERO
                        XB = XHK
                     ELSE ! XWMAX(N) >= HALF means flap bump
                        XA = XHK
                        XB = ONE
                     END IF

                  ELSE ! XBUMP(N) = 0. means shock-following line
                     XCENTER = XHK
                     XA = ZERO
                     XB = ONE
                  END IF

                  IF (DOLO(N) /= ZERO) THEN

                     DL = DELTAV * SIGN (ONE, DOLO(N))

                     CALL SFEVAL (TYP4, WBUMP(N), XCENTER,
     >                            XA, XB, DL, 1, NLOWK,
     >                            XINIT(1,K,IW), YINIT(1,K,IW), IWRIT)
                  END IF

                  IF (DOUP(N) /= ZERO) THEN

                     DU = DELTAV * SIGN (ONE, DOUP(N))

                     CALL SFEVAL (TYP4, WBUMP(N), XCENTER,
     >                            XA, XB, DU, NLOWK, NTOTK,
     >                            XINIT(1,K,IW), YINIT(1,K,IW), IWRIT)
                  END IF

               END IF ! Fixed vs. variable center or hinge X/C

            END IF ! Non-zero design variable contribution

         END DO ! Next wing or nacelle K station

      END DO ! Next design variable


C     --------------------------
C     Denormalize wing sections:
C     --------------------------

      DO IW = 1, NNWING

         DO K = 1, NWSEC(IW)

            NTOTK = NWG(K,IW)
            NLOWK = NLG(K,IW)
            NUPPK = NUG(K,IW)

C           Apply section twists and denormalize:

            ALPHA     = AINIT(K,IW)
            ALG(K,IW) = ALPHA
            THICK     = TINIT(K,IW)
            SCALE     = CINIT(K,IW)

            XXW = XINIT(NLOWK,K,IW)
            YYW = YINIT(NLOWK,K,IW)
            XLW = XLEADI(K,IW) + XXW * SCALE
            YLW = YLEADI(K,IW) + YYW * SCALE * THICK
            ZLW = ZLEADI(K,IW)
            CAB = COS (ALPHA*DEG2RAD)
            SAB = SIN (ALPHA*DEG2RAD)

            DO I = 1, NTOTK
               XWG(I,K,IW) = SCALE*((XINIT(I,K,IW) - XXW)*CAB +
     >                              (YINIT(I,K,IW) - YYW)*SAB * THICK)
     >                                                          + XLW
               YWG(I,K,IW) = SCALE*((YINIT(I,K,IW) - YYW)*CAB * THICK -
     >                              (XINIT(I,K,IW) - XXW)*SAB)  + YLW
               ZWG(I,K,IW) = SCALE*  ZINIT(I,K,IW)              + ZLW
            END DO


C           Roll the section about its leading edge (as for winglets)?

            IF (ROLL(K,IW) /= ZERO) THEN ! +ve is clockwise from upstream

               CALL ROTATE2D (NTOTK, ZWG(1,K,IW), YWG(1,K,IW),
     >                       -ROLL(K,IW), ZLW, YLW)
            END IF

         END DO ! Next section

      END DO ! Next wing-type component


C     ------------------------------------------------------------
C     Wing flap or slat deflections must be applied in real space:
C     ------------------------------------------------------------

      IF (NTOTFS > 0) THEN

         DO N = NCHCKS + NWVAR + NSINVC + 1, NTWING - NDIHED

            VS = V(N) * VSCALE(N)
            IW = IDWING(N)

C           Translate the deflection to the DFLAP/XFLAP/etc. quantities
C           of the affected K stations.  These quantities may have been
C           specified in the input geometry file, so actually applying the
C           the deflection later covers both possibilities:

            CALL DEFLECT (IWRIT, IW, N, VTYPE(N), VS,
     >                    KMIN(N), KMAX(N), XWMIN(N), XWMAX(N),
     >                    CINIT(1,IW), XLEADI(1,IW), ZLEADI(1,IW),
     >                    DFLAP(1,IW), XFLAP(1,IW), SWFLAP(1,IW),
     >                    DSLAT(1,IW), XSLAT(1,IW), SWSLAT(1,IW))
         END DO

      END IF

C     --------------------------------------------
C     Apply any wing flap and/or slat deflections:
C     --------------------------------------------

C     Careful:  the point counts may change, so don't use temporaries.

      DO IW = 1, NNWING

         DO K = 1, NWSEC(IW)

            IF (DFLAP(K,IW) /= ZERO) THEN ! Flap deflection; +ve is down

               CALL FLAPSLAT ('F', NLG(K,IW), NUG(K,IW), NWG(K,IW),
     >                        XWG(1,K,IW), YWG(1,K,IW), ZWG(1,K,IW),
     >                        DFLAP(K,IW), XFLAP(K,IW), SWFLAP(K,IW))
            END IF

            IF (DSLAT(K,IW) /= ZERO) THEN ! Slat deflection; +ve is down

               CALL FLAPSLAT ('S', NLG(K,IW), NUG(K,IW), NWG(K,IW),
     >                        XWG(1,K,IW), YWG(1,K,IW), ZWG(1,K,IW),
     >                        DSLAT(K,IW), XSLAT(K,IW), SWSLAT(K,IW))
            END IF

         END DO ! Next section

      END DO ! Next wing-type component

C     -------------------------------------------------------------
C     Applying dihedral last simplifies the above flap deflections:
C     -------------------------------------------------------------

      IF (NDIHED > 0) THEN

C        Translate the dihedral and the rotation axis to the variables
C        DIHED() and DIAXIS() of the affected K stations.  These quantities
C        may have been specified in the input geometry file, so actually
C        applying the the dihedral later covers both possibilities:

         DO N = NTWING - NDIHED + 1, NTWING

            VS = V(N) * VSCALE(N)
            IW = IDWING(N)
            KC = KCEN(N)
            I  = NLG(KC,IW)
            J  = NWG(KC,IW)

            DO K = KMIN(N), KMAX(N)

               DIHED(K,IW) = VS

               DIAXIS(1,1,K,IW) = XWG(I,KC,IW)
               DIAXIS(2,1,K,IW) = YWG(I,KC,IW)
               DIAXIS(3,1,K,IW) = ZWG(I,KC,IW)

               DIAXIS(1,2,K,IW) = (XWG(1,KC,IW) + XWG(J,KC,IW)) * HALF
               DIAXIS(2,2,K,IW) = (YWG(1,KC,IW) + YWG(J,KC,IW)) * HALF
               DIAXIS(3,2,K,IW) = (ZWG(1,KC,IW) + ZWG(J,KC,IW)) * HALF

            END DO

         END DO

      END IF

      DO IW = 1, NNWING
         DO K = 1, NWSEC(IW)
            VS = -DIHED(K,IW)    ! RH rule applies in ROTATE3D
            IF (VS /= ZERO) THEN ! Apply dihedral

               CALL ROTATE3D
     >           (NWG(K,IW), XWG(1,K,IW), YWG(1,K,IW), ZWG(1,K,IW), VS,
     >            DIAXIS(1,1,K,IW), DIAXIS(2,1,K,IW), DIAXIS(3,1,K,IW),
     >            DIAXIS(1,2,K,IW), DIAXIS(2,2,K,IW), DIAXIS(3,2,K,IW))

            END IF
         END DO
      END DO


C     --------------------------------------------------------------------------
C     Denormalize nacelle sections (still in cylindrical coords. at this stage):
C     --------------------------------------------------------------------------

      DO IW = NNWING + 1, NNWING + NNNAC

         DO K = 1, NWSEC(IW)

            NTOTK = NWG(K,IW)
            NLOWK = NLG(K,IW)
            NUPPK = NUG(K,IW)

C           Apply nacelle roll and section twists and denormalize as
C           cylindrical coordinates still in the nacelle axis system:

            ALPHA     = AINIT(K,IW)
            ALG(K,IW) = ALPHA
            THICK     = TINIT(K,IW)
            SCALE     = CINIT(K,IW)

            XLW =  XLEADI(K,IW)
            RLW =  YLEADI(K,IW)
            TLW = (ZLEADI(K,IW) + CLOCK(IW)) * TWOPI

            CAB = COS (ALPHA*DEG2RAD)
            SAB = SIN (ALPHA*DEG2RAD)

            DO I = 1, NTOTK
               XWG(I,K,IW) = SCALE*(XINIT(I,K,IW)*CAB +
     >                              YINIT(I,K,IW)*SAB*THICK) + XLW
               YWG(I,K,IW) = SCALE*(YINIT(I,K,IW)*CAB*THICK -
     >                              XINIT(I,K,IW)*SAB)       + RLW
               ZWG(I,K,IW) = TLW
            END DO

C           Roll the section about its leading edge?

            IF (ROLL(K,IW) /= ZERO) THEN ! +ve is clockwise from upstream

               CALL ROTATE2D (NTOTK, YWG(1,K,IW), ZWG(1,K,IW),
     >                       -ROLL(K,IW), RLW, TLW)
            END IF

         END DO ! Next nacelle section

      END DO ! Next nacelle


C     ---------------------------------------------------------
C     Regularize any wing/pylon/diverter/tail/nacelle sections:
C     ---------------------------------------------------------

      DO IW = 1, NNWING + NNNAC

         IL = ILWING(IW)
         NHALF = (IL + 1) / 2

         IF (NCALL == -3) THEN ! Normalized streamwise distribution

            CALL FOILGRD (NHALF, ZERO, ONE, CRDSPL(IW), CRDSPQ(IW),
     >                    CRDSPS(IW), CRDSPC(IW), SNORM(1,IW))
         END IF

         NI = IDIM(IW)
         NK = KDIM(IW)
         I1 = IRG(IW) ! Start of this surface in X/Y/ZRG(*)

         CALL WREGUL (1, NWSEC(IW), MXIWING, MXKWING,
     >                XWG(1,1,IW), YWG(1,1,IW), ZWG(1,1,IW),
     >                ROUNDED(1,IW), NLG(1,IW), NUG(1,IW), NWG(1,IW),
     >                NHALF, NHALF, SNORM(1,IW), NI, NK,
     >                XRG(I1), YRG(I1), ZRG(I1))

         IF (IW > NNWING) THEN

            KL = KLWING(IW) - 1 ! SURFER inserts pt. at pylon intersection

C           Convert regularized sections from cylindrical coordinates
C           to XYZ space and interpolate them to sections at uniform
C           azimuthal angles, still centered on the nacelle axis:

            CALL NCYLINDER (NI, NK, NWSEC(IW), IL, KL,
     >                      XRG(I1), YRG(I1), ZRG(I1))

C           Pitch/yaw/shift the redistributed nacelle surface, in-place:

            CALL NPITCHYAW (NI, NK, IL, KL,
     >                      XC1(IW), YC1(IW), ZC1(IW),
     >                      XC2(IW), YC2(IW), ZC2(IW),
     >                      XRG(I1), YRG(I1), ZRG(I1))

C           Fudge the nacelle to help close the diverter at the trailing edge?

            ID = NNROOT(IW) ! Diverter/pylon surface corresp. to nacelle IW

            IF (ID /= 0) THEN

               IWING = NNROOT(ID) ! Index of wing surface (or -1 for fuselage)

               IF (XTRAP(IW) < ZERO .AND. IWING > 0) THEN

                  I1W = IRG(IWING)

                  CALL WNFUDGE (NCALL, NWSURF, IW, IWING,
     >                          ILWING, KLWING, NWSEC,
     >                          NI, NK, IDIM(IWING), KDIM(IWING),
     >                          XRG(I1),  YRG(I1),  ZRG(I1),
     >                          XRG(I1W), YRG(I1W), ZRG(I1W),
     >                          XTRAP, ZINNER, ZOUTER)
               END IF

            END IF

         END IF

      END DO


C     -----------------------
C     Fuselage perturbations:
C     -----------------------

      IF (NNFUSE /= 0) THEN
         IF (ZFINIT(2,1) == ZFINIT(1,1)) THEN ! Avoid singular nose station
            IF1 = 2
         ELSE
            IF1 = 1
         END IF
      END IF

      IF (NNFUSE /= 0 .AND. NTOTFU /= 0) THEN

         XBLONG  = XF(NFSTN) - XF(1)
         RXBLONG = ONE / XBLONG

C        Body camber variables:
C        ----------------------

         DO N = NTWING + 1, NTWING + NCMBRF

            VS = V(N) * VSCALE(N)
            IF (VS == ZERO) CYCLE ! Next N

            TYPE = VTYPE(N)  ! 'BCxxxx'
            TYP4 = TYPE(3:6) ! 'xxxx'; SFEVAL expects NAME*4
            RXBRANGE = ONE / (XWMAX(N) - XWMIN(N)) ! XWMIN/MAX are in [0,1]

            DO I = 1, NFSTN
               XNORM = ((XF(I) - XF(1))*RXBLONG - XWMIN(N))*RXBRANGE
               IF (XNORM < ZERO) CYCLE ! Next station
               IF (XNORM > ONE)  EXIT  ! No more stations to do

               DELTAV = ZERO ! SFEVAL adds

               CALL SFEVAL (TYP4, WBUMP(N), XBUMP(N), ZERO, ONE,
     >                      XBLONG, 1, 1, XNORM, DELTAV, IWRIT)

               DELTAV = DELTAV * VS

               DO J = 1, NFJ(I)
                  YFINIT(J,I) = YFINIT(J,I) + DELTAV
               END DO

            END DO ! Next body section

         END DO ! Next body camber design variable


C        Body span variables:
C        --------------------

         DO N = NTWING + NCMBRF + 1, NTWING + NCMBRF + NSPANF

            VS = V(N) * VSCALE(N)
            IF (VS == ZERO) CYCLE ! Next N

            TYPE2 = UTYPE(N)  ! Circumferential = U in (U,V) parameterization
            TYPE  = VTYPE(N)  ! 'BSxxxx', 'BVxxx', or 'BRxxx'; V = axial
            TYP1  = TYPE(2:2) ! 'S', 'V', or 'R'
            TYP4  = TYPE(3:6) ! 'xxxx'; SFEVAL expects NAME*4
            WIDTH = PBUMP(N)
            CENTR = UBUMP(N)
            UMIN  = UBMIN(N)
            UMAX  = UBMAX(N)
            RXBRANGE = ONE / (XWMAX(N) - XWMIN(N))

            DO I = IF1, NFSTN

               NFP = NFJ(I)

C              Determine the bump height in the axial direction ...

               XDD = ((XF(I) - XF(1))*RXBLONG - XWMIN(N)) * RXBRANGE
               IF (XDD < ZERO) CYCLE ! Next I
               IF (XDD > ONE ) EXIT  ! Next N

               DZMAX = ZERO ! Max. dY, dZ, or dR at this X; SFEVAL adds

               CALL SFEVAL (TYP4, WBUMP(N), XBUMP(N), ZERO, ONE,
     >                      XBLONG, 1, 1, XDD, DZMAX, IWRIT)

C              ... and in the circumferential direction:

               CALL CHORDS2D (NFP, YFINIT(1,I), ZFINIT(1,I), TRUE,
     >                        DS, SFUS)

               RBUMP = ZERO ! (1:NFP); SFEVAL adds

               CALL SFEVAL (TYPE2, WIDTH, CENTR, UMIN, UMAX, DZMAX,
     >                      1, NFP, SFUS, RBUMP, IWRIT)


               IF      (TYP1 == 'S') THEN ! Spanwise
                                          ! --------
                  DO J = 1, NFP
                     ZFINIT(J,I) = ZFINIT(J,I) + VS * RBUMP(J)
                  END DO

               ELSE IF (TYP1 == 'V') THEN ! Vertical
                                          ! --------
                  DO J = 1, NFP
                     YFINIT(J,I) = YFINIT(J,I) + VS * RBUMP(J)
                  END DO

               ELSE   ! TYP1 = 'R'        ! (Radial)
                                          ! --------

C                 Originally, unit normals were chosen, but treatment of
C                 radial perturbations with linear constraints requires
C                 a common origin for the radial vectors, namely the
C                 body camber line:

C****             CALL UNITNORM ('B', NFP, ZFINIT(1,I), YFINIT(1,I),
C****>                           SFUS, YFUDGE, ZNORM, YNORM)

                  YCAM = (YFINIT(1,I) + YFINIT(NFP,I)) * HALF

                  YNORM(1) = -RBUMP(1) ! Crown/keel stays on sym. plane
                  ZNORM(1) = ZERO

                  DO J = 2, NFP - 1
                     ANGLE = ATAN2 (YFINIT(J,I) - YCAM, ZFINIT(J,I))
                     YNORM(J) = RBUMP(J) * SIN (ANGLE)
                     ZNORM(J) = RBUMP(J) * COS (ANGLE)
                  END DO

                  YNORM(NFP) = RBUMP(NFP)
                  ZNORM(NFP) = ZERO

                  DO J = 1, NFP
                     YFINIT(J,I) = VS * YNORM(J) + YFINIT(J,I)
                     ZFINIT(J,I) = VS * ZNORM(J) + ZFINIT(J,I)
                  END DO

               END IF

            END DO ! Next body section

         END DO ! Next body span design variable


C        Fuselage area modifications (sections scaled/same shape):
C        ---------------------------------------------------------

         DO N = NTWING + NCMBRF + NSPANF + 1,
     >          NTWING + NCMBRF + NSPANF + NAREAF

            VS = V(N) * VSCALE(N)
            IF (VS == ZERO) CYCLE ! Next N

            TYPE  = VTYPE(N)  ! 'BAxxxx'
            TYP4  = TYPE(3:6) ! 'xxxx'; SFEVAL expects NAME*4
            WIDTH = WBUMP(N)
            CENTR = XBUMP(N)
            AMULT = TWO * VS * XBLONG ** 2
            RXBRANGE = ONE / (XWMAX(N) - XWMIN(N))
            AREATOL  = MAX (100.* EPSILON (TOLER), TOLER)

            DO I = 2, NFSTN

               XDD = ((XF(I) - XF(1)) * RXBLONG - XWMIN(N)) * RXBRANGE
               IF (XDD < ZERO) CYCLE  ! Next body station
               IF (XDD > ONE)  EXIT   ! Next N

C              Calculate the unperturbed area, and set up angular coordinates
C              (or rather sines and cosines of them) for the section points:

               NFP   = NFJ(I)
               YKEEL = YFINIT(1,I)
               YCAM  = (YKEEL + YFINIT(NFP,I)) * HALF
               CCX   = YCAM - YKEEL
               TOTA  = ZERO
               SINTHETA(1) = ZERO
               COSTHETA(1) = ONE

               DO J = 2, NFP - 1
                  COF1 = YFINIT(J,I)   - YCAM
                  COF2 = ZFINIT(J,I)
                  COF3 = YFINIT(J-1,I) - YCAM
                  COF4 = ZFINIT(J-1,I)
                  TOTA = TOTA + (COF1*COF4 - COF3*COF2) ! Full section, not half

                  AAX  = (YFINIT(J,I)  - YKEEL) ** 2 + COF2 ** 2 ! AAX ** 2
                  BBX  = SQRT (            COF1 ** 2 + COF2 ** 2)

                  COSTHETA(J) = (AAX - BBX ** 2 - CCX ** 2) /
     >                          (-TWO * BBX * CCX)
                  SINTHETA(J) = SQRT (MAX (ZERO, ONE - COSTHETA(J) **2))
               END DO

               SINTHETA(NFP) = ZERO
               COSTHETA(NFP) = -ONE

               TOTR = SQRT (TOTA / PI) ! Average radius of unperturbed section

C              Determine desired area:

               BUMP1 = ZERO ! SFEVAL adds

               CALL SFEVAL (TYP4, WIDTH, CENTR, ZERO, ONE, ONE,
     >                      1, 1, XDD, BUMP1, IWRIT)

               TAGA = AMULT * BUMP1 + TOTA ! Desired area
               TAGR = SQRT (TAGA / PI)     ! Initial estimate of corresp. radius
               DIFR = TAGR - TOTR
               TAGB = TAGA

C              Iterate on averaged radius until the desired area is achieved:

CCC            write (iwrit,'(/,a,i4,2f14.6,1p,2e12.3)')
CCC     >         ' i, tota, taga, amult, bump1: ',
CCC     >           i, tota, taga, amult, bump1

               DO ITER = 1, 30

                  DO J = 1, NFP
                     YFINIT(J,I) = YFINIT(J,I) - DIFR * COSTHETA(J)
                     ZFINIT(J,I) = ZFINIT(J,I) + DIFR * SINTHETA(J)
                  END DO

                  TOTB = ZERO
                  DO J = 2, NFP
                     COF1 = YFINIT(J,I)   - YCAM
                     COF2 = ZFINIT(J,I)
                     COF3 = YFINIT(J-1,I) - YCAM
                     COF4 = ZFINIT(J-1,I)
                     TOTB = TOTB + (COF1*COF4 - COF3*COF2)
                  END DO

                  DIFA = TOTB - TAGA

CCC               write (iwrit,'(1x,a,i3,2f14.6,2f12.6)')
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
                     WRITE (IWRIT,'(2A, /, (A, I5))')
     >                  ' No convergence on target fuselage',
     >                  ' area distribution.',
     >                  ' Design variable:', N,
     >                  ' Aborting at fuselage station', I
                     GO TO 990
                  END IF

               END DO ! Next body area iteration

            END DO ! Next body station

         END DO ! Next body area design variable

      END IF ! Fuselage perturbation block


      NEWBODY = .FALSE.

      IF (NNFUSE /= 0 .AND. NCALL == -3) THEN ! Avoid for NNREPEAT = 1 case

C        ---------------------------------------------
C        Option to protect the wing/body intersection:
C        ---------------------------------------------

C        SINGLE FUSELAGE/SINGLE WING COMPONENT ASSUMED:

         IF (NNFUDGE /= 0) THEN ! Adjust some body sections in-place

            IW = NNFUDGE

            CALL WBFUDGE (NCALL, IBOD1, IBOD2, IBOD3, MXJBODY, MXIBODY,
     >                    NFSTN, NFJ, XF, YFINIT, ZFINIT, KWING2(IW),
     >                    MXIWING, MXKWING, XWG(1,1,IW), YWG(1,1,IW),
     >                    ZWG(1,1,IW), NWG(1,IW), NLG(1,IW),
     >                    NSHELF, XSHELF, YSHELF, ZSHELF,
     >                    IWRIT, NEWBODY)
         END IF

C        Regularize the fuselage sections:

         IF (NCALL == -3) THEN

            IF (IDEGFUS == 3) THEN
               METHOD = 'B'  ! Loose local cubic spline fits
            ELSE IF (IDEGFUS == 2) THEN  ! Kludge
               METHOD = 'M'  ! Monotonic ...
            ELSE   ! IDEGFUS = 1
               METHOD = 'L'  ! Piecewise linear
            END IF

         END IF

         IF (NNCONFIG >= 1) THEN ! CTV cases with probable singular line at nose
            NJBODY = JLBODY
         ELSE
            NJBODY = JLBODY - 1 ! Allow for inserting a point <-> root LE
         END IF

C        If first call, or body has changed ...

         IF (NCALL == -3 .OR. NEWBODY) THEN

            IF (IF1 /= 1) THEN ! Treat singular nose point
               YF(1:NJBODY,1) = YFINIT(1,1)
               ZF(1:NJBODY,1) = ZFINIT(1,1)
            END IF

            SFUS(1)  = ZERO
            SUNIF(1) = ZERO

            DO I = IF1, NFSTN

               NJ = NFJ(I)

               CALL CHORDS2D (NJ, YFINIT(1,I), ZFINIT(1,I), .FALSE.,
     >                        DS, SFUS)

               SUNIF(NJBODY) = DS ! Total arc length
               DS = DS / REAL (NJBODY - 1)

               DO J = 2, NJBODY - 1
                  SUNIF(J) = DS * REAL (J - 1)
               END DO

               CALL LCSFIT (NJ, SFUS, YFINIT(1,I), TRUE, METHOD, NJBODY,
     >                      SUNIF, YF(1,I), DERIVS)

               CALL LCSFIT (NJ, SFUS, ZFINIT(1,I), TRUE, METHOD, NJBODY,
     >                      SUNIF, ZF(1,I), DERIVS)

               ZF(NJBODY,I) = ZERO ! Exactly

            END DO

         END IF

      END IF

      RETURN

  990 WRITE (IWRIT, '(/, A, I3)') ' Aborting on surface IW =', IW

      STOP


C     PERTURB internal procedure:
C     ---------------------------

      CONTAINS

!        -----------------------------------------------------------------------
         SUBROUTINE XHINGE_AT_THIS_K () ! Sets XHK = Hinge line X/C at station K
!        -----------------------------------------------------------------------

!        Before the first call for a given N, set KI & KO = KMIN(N) & KMAX(N).

!        Local variables:

         REAL, SAVE :: RH, XHI, XHO, ZLI

!        Execution:

         IF (K == KI) THEN ! Some values are lost after section KI is treated

            XHI = XLEADI(KI,IW) + XWMIN(N) * CINIT(KI,IW)
            XHO = XLEADI(KO,IW) + XWMAX(N) * CINIT(KO,IW)
            ZLI = ZLEADI(KI,IW)
            RH  = (XHO - XHI) / (ZLEADI(KO,IW) - ZLI)

         END IF

         XHK = XHI + (ZLEADI(K,IW)  - ZLI) * RH
         XHK = (XHK - XLEADI(K,IW)) / CINIT(K,IW) ! Corresp. X/C

         END SUBROUTINE XHINGE_AT_THIS_K

      END SUBROUTINE PERTURB

C***********************************************************************
C
      SUBROUTINE RDVARS (IREAD, IWRIT)
C
C     File-driven control scheme for wing/body design variables.
C     See PERTURB or SYN107-MB for the descriptions.  Variable order
C     is important.
C
C***********************************************************************

      USE SHAPE_FUNCTIONS

      IMPLICIT NONE

C     Argument:

      INTEGER, INTENT (IN) ::
     >   IREAD, IWRIT

C     Local variables:

      INTEGER, ALLOCATABLE ::
     >   IGROUP(:)

      INTEGER
     >   I, N, IER

      REAL
     >   RIDWING

      LOGICAL
     >   ORDERED

      CHARACTER
     >   TYPE * 6, TYPE2 * 4

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
      NSINEP = 0
      NEXPSP = 0
      NTRAIP = 0
      NLEADP = 0
      NSINVC = 0
      NFLAP  = 0
      NSLAT  = 0
      NDIHED = 0
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

C     Allocate associated arrays:

      IF (NDV > 0) THEN

         ALLOCATE (KCEN(NDV), KMIN(NDV), KMAX(NDV), IDWING(NDV), V(NDV),
     >             VSCALE(NDV), AITCH(NDV), DOUP(NDV), DOLO(NDV),
     >             WBUMP(NDV), XBUMP(NDV), XWMIN(NDV), XWMAX(NDV),
     >             PBUMP(NDV), UBUMP(NDV), UBMIN(NDV), UBMAX(NDV),
     >             UTYPE(NDV), VTYPE(NDV), IGROUP(NDV))

         READ (IREAD, *)

      END IF

      DO I = 1, NDV

         READ (IREAD, *) N, TYPE, TYPE2, DOUP(I), DOLO(I),
     >                   WBUMP(I), XBUMP(I), XWMIN(I), XWMAX(I),
     >                   PBUMP(I), UBUMP(I), UBMIN(I), UBMAX(I),
     >                   V(I), VSCALE(I), AITCH(I), RIDWING

         IDWING(I) = RIDWING
         VTYPE(I)  = TYPE
         UTYPE(I)  = TYPE2
         KCEN(I)   = 0 ! PERTURB won't plug Zs into UBUMP, etc., unless ...

         IF      (TYPE2 == 'STAT') THEN ! James's choice for "station #, not Z"
            UTYPE(I) = 'POLY'
            KCEN(I)  = INT (UBUMP(I))
            KMIN(I)  = INT (UBMIN(I))
            KMAX(I)  = INT (UBMAX(I))
         ELSE IF (TYPE2 == 'POLY') THEN ! Just check for bad inputs
         ELSE IF (TYPE2 == 'SIN ') THEN
         ELSE IF (TYPE2 == 'SINF') THEN
         ELSE IF (TYPE2 == 'SIN1') THEN
         ELSE IF (TYPE2 == 'SIN2') THEN
         ELSE IF (TYPE2 == 'COSL') THEN
         ELSE IF (TYPE2 == 'COSR') THEN
         ELSE IF (TYPE2 == 'LCOS') THEN
         ELSE IF (TYPE2 == 'RCOS') THEN
         ELSE IF (TYPE2 == 'LED ') THEN
         ELSE IF (TYPE2 == 'TRL ') THEN
         ELSE IF (TYPE2 == 'EXP ') THEN
         ELSE
            GO TO 800 ! Unknown variable UTYPE
         END IF

         IF (XBUMP(I) < XWMIN(I) .OR. XBUMP(I) > XWMAX(I) .OR.
     >       UBUMP(I) < UBMIN(I) .OR. UBUMP(I) > UBMAX(I)) THEN
            N = I
            GO TO 850
         END IF

         IF      (TYPE == 'TWT   ') THEN
            NTWIST = NTWIST + 1
            IGROUP(I) = 1
         ELSE IF (TYPE == 'TWT2  ') THEN
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
         ELSE IF (TYPE == 'SIN   ') THEN
            NSINE  = NSINE  + 1
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
         ELSE IF (TYPE == 'DIHED ') THEN
            NDIHED = NDIHED + 1
            IGROUP(I) = 6
         ELSE IF (TYPE == 'BCSIN ') THEN
            NSINFC = NSINFC + 1
            IGROUP(I) = 7
         ELSE IF (TYPE == 'BCSIN1') THEN
            NSINFC = NSINFC + 1
            IGROUP(I) = 7
         ELSE IF (TYPE == 'BCSIN2') THEN
            NSINFC = NSINFC + 1
            IGROUP(I) = 7
         ELSE IF (TYPE == 'BCSINF') THEN
            NSINFC = NSINFC + 1
            IGROUP(I) = 7
         ELSE IF (TYPE == 'BCCOS ') THEN
            NCOSFC = NCOSFC + 1
            IGROUP(I) = 7
         ELSE IF (TYPE == 'BCEXP ') THEN
            NEXPFC = NEXPFC + 1
            IGROUP(I) = 7
         ELSE IF (TYPE == 'BCLED ') THEN
            NLEDFC = NLEDFC + 1
            IGROUP(I) = 7
         ELSE IF (TYPE == 'BCTRL ') THEN
            NTRLFC = NTRLFC + 1
            IGROUP(I) = 7
         ELSE IF (TYPE == 'BSSIN ') THEN
            NSINFS = NSINFS + 1
            IGROUP(I) = 8
         ELSE IF (TYPE == 'BSSIN1') THEN
            NSINFS = NSINFS + 1
            IGROUP(I) = 8
         ELSE IF (TYPE == 'BSSIN2') THEN
            NSINFS = NSINFS + 1
            IGROUP(I) = 8
         ELSE IF (TYPE == 'BSSINF') THEN
            NSINFS = NSINFS + 1
            IGROUP(I) = 8
         ELSE IF (TYPE == 'BSCOS ') THEN
            NCOSFS = NCOSFS + 1
            IGROUP(I) = 8
         ELSE IF (TYPE == 'BSEXP ') THEN
            NEXPFS = NEXPFS + 1
            IGROUP(I) = 8
         ELSE IF (TYPE == 'BSLED ') THEN
            NLEDFS = NLEDFS + 1
            IGROUP(I) = 8
         ELSE IF (TYPE == 'BSTRL ') THEN
            NTRLFS = NTRLFS + 1
            IGROUP(I) = 8
         ELSE IF (TYPE == 'BVSIN ') THEN
            NSINFS = NSINFS + 1
            IGROUP(I) = 8
         ELSE IF (TYPE == 'BVSIN1') THEN
            NSINFS = NSINFS + 1
            IGROUP(I) = 8
         ELSE IF (TYPE == 'BVSIN2') THEN
            NSINFS = NSINFS + 1
            IGROUP(I) = 8
         ELSE IF (TYPE == 'BVSINF') THEN
            NSINFS = NSINFS + 1
            IGROUP(I) = 8
         ELSE IF (TYPE == 'BVCOS ') THEN
            NCOSFS = NCOSFS + 1
            IGROUP(I) = 8
         ELSE IF (TYPE == 'BVEXP ') THEN
            NEXPFS = NEXPFS + 1
            IGROUP(I) = 8
         ELSE IF (TYPE == 'BVLED ') THEN
            NLEDFS = NLEDFS + 1
            IGROUP(I) = 8
         ELSE IF (TYPE == 'BVTRL ') THEN
            NTRLFS = NTRLFS + 1
            IGROUP(I) = 8
         ELSE IF (TYPE == 'BRSIN ') THEN
            NSINFS = NSINFS + 1
            IGROUP(I) = 8
         ELSE IF (TYPE == 'BRSIN1') THEN
            NSINFS = NSINFS + 1
            IGROUP(I) = 8
         ELSE IF (TYPE == 'BRSIN2') THEN
            NSINFS = NSINFS + 1
            IGROUP(I) = 8
         ELSE IF (TYPE == 'BRSINF') THEN
            NSINFS = NSINFS + 1
            IGROUP(I) = 8
         ELSE IF (TYPE == 'BRCOS ') THEN
            NCOSFS = NCOSFS + 1
            IGROUP(I) = 8
         ELSE IF (TYPE == 'BREXP ') THEN
            NEXPFS = NEXPFS + 1
            IGROUP(I) = 8
         ELSE IF (TYPE == 'BRLED ') THEN
            NLEDFS = NLEDFS + 1
            IGROUP(I) = 8
         ELSE IF (TYPE == 'BRTRL ') THEN
            NTRLFS = NTRLFS + 1
            IGROUP(I) = 8
         ELSE IF (TYPE == 'BASIN ') THEN
            NSINFA = NSINFA + 1
            IGROUP(I) = 9
         ELSE IF (TYPE == 'BASIN1') THEN
            NSINFA = NSINFA + 1
            IGROUP(I) = 9
         ELSE IF (TYPE == 'BASIN2') THEN
            NSINFA = NSINFA + 1
            IGROUP(I) = 9
         ELSE IF (TYPE == 'BASINF') THEN
            NSINFA = NSINFA + 1
            IGROUP(I) = 9
         ELSE IF (TYPE == 'BACOS ') THEN
            NCOSFA = NCOSFA + 1
            IGROUP(I) = 9
         ELSE IF (TYPE == 'BAEXP ') THEN
            NEXPFA = NEXPFA + 1
            IGROUP(I) = 9
         ELSE IF (TYPE == 'BALED ') THEN
            NLEDFA = NLEDFA + 1
            IGROUP(I) = 9
         ELSE IF (TYPE == 'BATRL ') THEN
            NTRLFA = NTRLFA + 1
            IGROUP(I) = 9
         ELSE IF (TYPE == 'XFLR  ') THEN
            NXFLR = NXFLR + 1
            IGROUP(I) = 10
         ELSE IF (TYPE == 'YFLR  ') THEN
            NYFLR = NYFLR + 1
            IGROUP(I) = 10
         ELSE
            GO TO 800 ! Unknown variable type
         END IF

      END DO

      NCHCKS = NTWIST + NTHICK + NXLEAD + NZLEAD + NYLEAD + NCHORD
      NCHCK1 = NSINE  + NCOS   + NEXPS  + NTRAIL + NLEAD  + NWAG
      NCHCK2 = NSINEP + NEXPSP + NTRAIP + NLEADP
      NTOTFS = NFLAP  + NSLAT
      NTWING = NCHCKS + NCHCK1 + NCHCK2 + NSINVC + NTOTFS + NDIHED
      NCMBRF = NSINFC + NCOSFC + NEXPFC + NLEDFC + NTRLFC
      NSPANF = NSINFS + NCOSFS + NEXPFS + NLEDFS + NTRLFS
      NAREAF = NSINFA + NCOSFA + NEXPFA + NLEDFA + NTRLFA
      NTOTFU = NCMBRF + NSPANF + NAREAF
      NTKINK = NXFLR  + NYFLR

      IF (NCHCK1 /= 0 .AND. NCHCK2 /= 0) THEN
         WRITE (IWRIT, '(/, A)')
     >      ' RDVARS:  Y & Y" may not be modified simultaneously.'
         GO TO 900
      END IF

      N = NTWING + NTOTFU + NTKINK
      IF (N /= NDV) THEN
         WRITE (IWRIT, '(/, A, I6)')
     >      ' RDVARS:  Total # design variables did not match NDV:', N
         GO TO 900
      END IF

      ORDERED = .TRUE.
      DO I = 2, NDV
         IF (IGROUP(I) < IGROUP(I-1)) ORDERED = .FALSE.
      END DO

      IF (NDV > 0) DEALLOCATE (IGROUP)

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
     >   ' RDVARS: Unknown design variable type(s): ', I, TYPE, TYPE2
      GO TO 900

  850 WRITE (IWRIT, '(/, A, I5)')
     >   ' RDVARS: X or Z center is outside [min, max]. Variable:', N

  900 STOP

  999 RETURN

      END SUBROUTINE RDVARS

C***********************************************************************
C
      SUBROUTINE READ_SIZES (NNPANEL, MXPATCH, MXSUBI, MXSUBJ,
     >                       LUN, ICRT, IFORMAT)
C
C     Read patch/subpatch sizes from aero.xy0 if NNPANEL = 2, else
C     just allocate the I/JMXCOM(*) arrays still needed by SURFER.
C
C     10/13/99  D.Saunders  Initial implementation.
C     07/16/03      "       Possible differences between maximum and
C                           final counts forced working harder here.
C
C***********************************************************************

      USE GRIDS

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   NNPANEL,            ! So we can do the unused allocates here
     >   MXPATCH,            ! Maximums set aside for the initial paneling
     >   MXSUBI, MXSUBJ,     ! and used here for the aero.xy0 data arrays
     >   LUN,                ! Logical unit for temporary use
     >   ICRT,               ! For diagnostics
     >   IFORMAT             ! Even = unformatted; odd = formatted

C     Local constants:

      CHARACTER, PARAMETER ::
     >   FORM * 11 = 'UNFORMATTED'

C     Local variables:

      INTEGER
     >   I, J, IOS, MXSUBIO, MXSUBJO, N, NPATCHO

      REAL
     >   X

      LOGICAL
     >   FORMATTED

C     Execution:

      ALLOCATE (IMXCOMO(MXPATCH), JMXCOMO(MXPATCH),
     >          ISUBO(0:MXSUBI,MXPATCH), JSUBO(0:MXSUBJ,MXPATCH))

      IF (NNPANEL /= 2) THEN ! Accept whatever patch sizes SURFER gets

         IMXCOMO = 1 ! Surfer looks for changed dimensions
         JMXCOMO = 1
         ISUBO   = 1
         JSUBO   = 1

         GO TO 700 ! Done; avoid the indent for NNPANEL = 2

      END IF


      FORMATTED = MOD (IFORMAT, 2) /= 0

      IF (FORMATTED) THEN
         N = 3
      ELSE
         N = 1
      END IF

      OPEN (LUN, FILE='aero.xy0', FORM=FORM(N:11), STATUS='OLD',
     >      IOSTAT=IOS)

      IF (IOS /= 0) THEN
         WRITE (ICRT, '(/, A, I6)')
     >      ' READ_SIZES: Trouble opening aero.xy0.  IOSTAT:', IOS
         GO TO 900
      END IF

      IF (FORMATTED) THEN
         READ (LUN, *, IOSTAT=IOS) NPATCHO
      ELSE
         READ (LUN,    IOSTAT=IOS) NPATCHO
      END IF

      IF (IOS /= 0) THEN
         WRITE (ICRT, '(/, A, I6)')
     >   ' READ_SIZES: Trouble reading NPATCH from aero.xy0.  IOSTAT:',
     >   IOS
         GO TO 900
      END IF

      IF (NPATCHO /= MXPATCH) THEN
         WRITE (ICRT, '(/, (A, I10))')
     >      ' # patches found in aero.xy0:', NPATCHO,
     >      ' # estimated by AEROSURF:    ', MXPATCH
         IF (MXPATCH < NPATCHO) GO TO 900
         WRITE (ICRT, '(/, A)') ' Proceeding with the larger figure.'
      END IF

      IF (FORMATTED) THEN
         READ (LUN, *, IOSTAT=IOS)
     >      (IMXCOMO(N), JMXCOMO(N), I, N = 1, NPATCHO)
      ELSE
         READ (LUN,    IOSTAT=IOS)
     >      (IMXCOMO(N), JMXCOMO(N), I, N = 1, NPATCHO)
      END IF

      IF (IOS /= 0) GO TO 800

C     Skip the patches:

      DO N = 1, NPATCHO

         IF (FORMATTED) THEN
            READ (LUN, *, IOSTAT=IOS)
     >         (X, I = 1, 3 * IMXCOMO(N) * JMXCOMO(N)) ! Can't just skip
         ELSE
            READ (LUN, IOSTAT=IOS)
         END IF

         IF (IOS /= 0) GO TO 800

      END DO

      IF (FORMATTED) THEN
         READ (LUN, '(A)',    IOSTAT=IOS) ! 'Component numbers:'
         IF (IOS /= 0) GO TO 800
         READ (LUN, '(26I3)', IOSTAT=IOS) (J, I = 1, NPATCHO)
         IF (IOS /= 0) GO TO 800
         READ (LUN, '(A)',    IOSTAT=IOS) ! 'MXSUBI   MXSUBJ:'
         IF (IOS /= 0) GO TO 800
         READ (LUN, *,        IOSTAT=IOS)    MXSUBIO, MXSUBJO
         IF (IOS /= 0) GO TO 800
      ELSE
         READ (LUN, IOSTAT=IOS) ! IPIECE(1:NPATCHO)
         IF (IOS /= 0) GO TO 800
         READ (LUN, IOSTAT=IOS) MXSUBIO, MXSUBJO
         IF (IOS /= 0) GO TO 800
      END IF

      IF (MXSUBIO /= MXSUBI .OR. MXSUBJO /= MXSUBJ) THEN
         WRITE (ICRT, '(/, (A, 2I10))')
     >      ' Max. # subpatches found in aero.xy0:', MXSUBIO, MXSUBJO,
     >      ' Max. #s estimated by AEROSURF:      ', MXSUBI,  MXSUBJ
         GO TO 900
      END IF

      IF (FORMATTED) THEN
         READ (LUN, '(A)',    IOSTAT=IOS) ! 'NSUBIO(1:NPATCHO):'
         IF (IOS /= 0) GO TO 800
         READ (LUN, '(26I3)', IOSTAT=IOS)  (J, I = 1, NPATCHO)
         IF (IOS /= 0) GO TO 800
         READ (LUN, '(A)',    IOSTAT=IOS) ! 'NSUBJO(1:NPATCHO):'
         IF (IOS /= 0) GO TO 800
         READ (LUN, '(26I3)', IOSTAT=IOS)  (J, I = 1, NPATCHO)
         IF (IOS /= 0) GO TO 800
         READ (LUN, '(A)',    IOSTAT=IOS) ! 'Subpatch Is:'
         IF (IOS /= 0) GO TO 800
         READ (LUN, '(13I6)', IOSTAT=IOS) ISUBO(0:MXSUBI, 1:NPATCHO)
         IF (IOS /= 0) GO TO 800
         READ (LUN, '(A)',    IOSTAT=IOS) ! 'Subpatch Js:'
         IF (IOS /= 0) GO TO 800
         READ (LUN, '(13I6)', IOSTAT=IOS) JSUBO(0:MXSUBJ, 1:NPATCHO)
         IF (IOS /= 0) GO TO 800
      ELSE
         READ (LUN, IOSTAT=IOS) (J, I = 1, NPATCHO) ! NSUBIO(1:NPATCHO)
         IF (IOS /= 0) GO TO 800
         READ (LUN, IOSTAT=IOS) (J, I = 1, NPATCHO) ! NSUBJO(1:NPATCHO)
         IF (IOS /= 0) GO TO 800
         READ (LUN, IOSTAT=IOS) ISUBO(0:MXSUBI, 1:NPATCHO)
         IF (IOS /= 0) GO TO 800
         READ (LUN, IOSTAT=IOS) JSUBO(0:MXSUBJ, 1:NPATCHO)
         IF (IOS /= 0) GO TO 800
      END IF

      CLOSE (LUN)

  700 RETURN

  800 WRITE (ICRT, '(/, A, I6)')
     >   ' READ_SIZES: Trouble reading aero.xy0.  IOSTAT:', IOS
  900 STOP

      END SUBROUTINE READ_SIZES

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
C
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
C     11/03/03     "     Generalized the COS* and *COS functions to extend
C                        the peak value any distance.  E.g., for COSL,
C                        XC = XA gives the original form, else COSL = 1.0
C                        for X in [XA, XC], going to 0.0 at XB.
C     11/06/03     "     XC units should match those of XA, XB for the
C                        generalized COS* and *COS functions.
C
C     Author:  David Saunders, Sterling Software/NASA Ames, Moffett Field, CA
C
C*******************************************************************************

      IMPLICIT NONE

C     Arguments:

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

C     Local constants:

      REAL, PARAMETER ::
     >   HALF = 0.5, ONE = 1., TWO = 2., ZERO = 0.

C     Local variables:

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

C     Storage:

      SAVE
     >   FIRST, DEGRAD, PI, PIBY2, RPI, PT5LOG

C     Statement functions:

      REAL
     >   EBUMP, SBUMP, PWR, WDTH, XNRM

      EBUMP (WDTH, PWR, XNRM) =
     >   XNRM ** PWR * (ONE - XNRM) * EXP (-WDTH * XNRM)

      SBUMP (WDTH, PWR, XNRM) = (SIN (PI * (XNRM ** PWR))) ** WDTH

C     Execution:
C     ----------

      IF (FIRST) THEN

         FIRST = .FALSE.

C        Protect sine bump exponentiation from slightly negative arguments
C        by calculating a value for PI such that SIN (PI) > 0 (barely).
C        This allows the limiting of X to [0, 1] that is necessary for
C        the sub-range capability to cover the sine protection too.

         EPSMCH = EPSILON (EPSMCH)

         PI = ASIN (ONE) * TWO

         DO I = 1, 999 ! 6 iterations observed for DEC AXP    (32 bits)
                       ! 5     "         "      "  CRAY C90   (64 bits)
                       ! 1     "         "      "  SGI R10000 ( "   " )

C***        write (6,*) 'I   PI  SIN (PI): ', i, pi, sin (pi)

            IF (SIN (PI) > ZERO) EXIT

            PI = PI - EPSMCH * REAL (I)

         END DO

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
*        This version gives the peak value for X (I) in [XA, XC].

         XSCALE = PIBY2 / (XB - XC)
         X0 = TWO * XC - XB

         DO I = I1, I2
            XI = MAX (PIBY2, MIN ((X (I) - X0) * XSCALE, PI))
            Y (I) = SCALE * (SIN (XI)) ** WIDTH + Y (I)
         END DO

      ELSE IF (KEY4 == 'COSR') THEN ! Half a cosine (or sine), peak at right

C        SIN is used instead of COS for consistency with COSL.
C        This version gives the peak value for X (I) in [XC, XB].

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
      SUBROUTINE WBFUDGE (NCALL, IBOD1, IBOD2, IBOD3, MXJBODY, MXIBODY,
     >                    NFSTN, NFJ, XF, YBODY, ZBODY, KWING2,
     >                    MXIWING, MXKWING, XWG, YWG, ZWG, NWG, NLG,
     >                    NSHELF, XSHELF, YSHELF, ZSHELF,
     >                    IWRIT, NEWBODY)
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
C     08/07/00   "   Introduced IBOD3 option.
C     08/09/00   "   Generalized YSHELF/ZSHELF as functions of X.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   NCALL,                  ! -3 means first call for this case
     >   IBOD1, IBOD2,           ! Range of body sections to treat
     >   IBOD3,                  ! IBOD3 > IBOD2 overrides possible suppression
                                 ! of the shelf towards the trailing edge; it
                                 ! it is decayed smoothly to zero at stn. IBOD3
                                 ! except IBOD3 = 999 suppresses the decay
     >   MXJBODY,                ! Circumferential/axial dimensions of
     >   MXIBODY,                ! fuselage geometry arrays
     >   NFSTN,                  ! Active number of body sections
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
     >   NWG(MXKWING),           ! Numbers of wing section points (wrap)
     >   NLG(MXKWING)            ! Numbers on lower surfaces

      INTEGER, INTENT (IN) ::
     >   NSHELF

      REAL, INTENT (INOUT), DIMENSION (NSHELF) ::
     >   XSHELF                  ! Denormalized here on first call

      REAL, INTENT (IN), DIMENSION (NSHELF) ::
     >   YSHELF, ZSHELF          ! YSHELF(*) = safety margin (e.g. 1")
                                 ! included in any fudged body Y coordinate;
                                 ! ZSHELF(*) = fraction (e.g. 0.8) applied to
                                 ! local body width to give a minimum shelf
                                 ! width;  YSHELF & ZSHELF are splined wrt X
      INTEGER, INTENT (IN) ::
     >   IWRIT                   ! For diagnostics

      LOGICAL, INTENT (OUT) ::
     >   NEWBODY                 ! T means some body point(s) changed

C     Local constants:

      REAL,      PARAMETER :: ONE = 1., TOLER = 1.E-7
      LOGICAL,   PARAMETER :: NEW = .TRUE.
      CHARACTER, PARAMETER :: METHOD * 1 = 'M'

C     Local variables:

      INTEGER
     >   I, IEND, IER, ILAST, ILASTLO, ILASTUP, IWLO, IWUP, J, K,
     >   KLAST, KW, NFP, NFUDGE, NL, NTOT

      REAL
     >   EPS, P, Q, XLE, XPT, XTE, YEXTRA, YJI, YLO, YUP,
     >   ZBODMAX, ZFRACTION, ZLAST, ZMIN, ZPT

      REAL, ALLOCATABLE, SAVE, DIMENSION (:) ::
     >   YDIFF, ZDIFF

      LOGICAL
     >   NO_DECAY, NO_SHELF_ABOVE_WING, PREVIOUSJ

C     Procedures:

      EXTERNAL
     >   LCSFIT, RIPPLE2D

C     Storage:

      SAVE
     >   EPS, NL, NO_DECAY, NO_SHELF_ABOVE_WING, NTOT

C     Execution:

C     Check for feasibility of interpolating the wing lower surface, etc.

      IF (NCALL == -3) THEN

         EPS = MAX (5.* EPSILON (TOLER), TOLER) ! For p/q in [0,1] in RIPPLE2D
         NO_SHELF_ABOVE_WING = IBOD3 == 0
         NO_DECAY = IBOD3 == 999

         NL   = NLG(1)
         NTOT = NWG(1)

         DO K = 2, KWING2
            IF (NLG(K) /= NL)   GO TO 900
            IF (NWG(K) /= NTOT) GO TO 900
         END DO

C        Denormalize the abscissas used for defining the "shelf":

         P = XF(IBOD2) - XF(IBOD1)
         DO I = 1, NSHELF
            XSHELF(I) = XF(IBOD1) + P * XSHELF(I)
         END DO

      END IF

      IF (.NOT. NO_SHELF_ABOVE_WING) THEN
         ALLOCATE (YDIFF(MXJBODY), ZDIFF(MXJBODY)) ! For decaying the shelf
      END IF

      IWLO    = NL
      IWUP    = NL
      ILASTLO = IWLO ! In case first search fails
      KLAST   = 1
      NFUDGE  = 0

      XLE = MIN (XWG(NL,1), XWG(NL,KWING2))
      XTE = MAX (XWG (1,1), XWG (1,KWING2))

      DO I = IBOD1, IBOD2

         XPT = XF(I)
         IF (XPT < XLE) CYCLE ! Saves some wasted searching

         IF (XPT > XTE) EXIT

         NFP = NFJ(I)

         IF (.NOT. NO_SHELF_ABOVE_WING) THEN ! Save the last unfudged section
            DO J = 1, NFP
               YDIFF(J) = YBODY(J,I)
               ZDIFF(J) = ZBODY(J,I)
            END DO
         END IF

         KW = 1
         ZLAST = -ONE

         ZBODMAX = ZBODY(1,I)

         DO J = 2, NFP
            ZBODMAX = MAX (ZBODY(J,I), ZBODMAX)
         END DO

C        Allow the fractional width and the Y safety margin to vary with X:

         CALL LCSFIT (NSHELF, XSHELF, YSHELF, NEW, METHOD, 1, XPT,
     >                YEXTRA, Q)

         CALL LCSFIT (NSHELF, XSHELF, ZSHELF, NEW, METHOD, 1, XPT,
     >                ZFRACTION, Q)

         ZMIN = ZBODMAX * ZFRACTION

         PREVIOUSJ = .FALSE.

         DO J = 1, NFP

            ZPT = ZBODY(J,I)

            IF (ZPT > ZMIN) THEN ! Not true for J = 1 or 2

               IF (PREVIOUSJ) THEN ! Previous J was fudged

C                 Try to capture the edge of the shelf more precisely:

                  P = (ZMIN - ZBODY(J-2,I)) /
     >                MAX (1.E-6, ZBODY(J-1,I) - ZBODY(J-2,I))
                  YJI = (ONE - P) * YBODY(J-2,I) + P * YBODY(J-1,I)
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
     >            (ONE-Q)*((ONE-P)*YWG(IWLO,KW)   + P*YWG(IWLO+1,KW))  +
     >                 Q *((ONE-P)*YWG(IWLO,KW+1) + P*YWG(IWLO+1,KW+1))

               IF (YJI > YLO - YEXTRA) THEN ! Quit if above upper wing too?

                  IF (NO_SHELF_ABOVE_WING) THEN

                     IWUP = MIN (NTOT - 1, 2*NL - IWLO)

                     CALL RIPPLE2D (MXIWING, MXKWING, NL, NTOT, 1, ! Upper
     >                              KWING2, XWG, ZWG, XPT, ZPT,
     >                              IWUP, KW, EPS, P, Q, IER)

                     YUP = (ONE-Q)*((ONE-P)*YWG(IWUP,KW)    +
     >                                   P *YWG(IWUP+1,KW)) +
     >                          Q *((ONE-P)*YWG(IWUP,KW+1)  +
     >                                   P *YWG(IWUP+1,KW+1))

                     IF (YJI >= YUP) EXIT ! Avoid TE trouble for large ZSHELF

                  END IF

                  YBODY(J,I) = YLO - YEXTRA
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

         ILAST = I

      END DO ! Next I

      NEWBODY = NFUDGE > 0

      IF (NCALL <= 0) THEN
         IF (IWRIT > 0) WRITE (IWRIT, '(/, A, I6)')
     >    ' WBFUDGE: # Body points fudged to be below the wing:', NFUDGE
      END IF

      IF (.NOT. NO_SHELF_ABOVE_WING) THEN

         IF (NEWBODY) THEN

            NFP = NFJ(ILAST)

            DO J = 1, NFP
               YDIFF(J) = YBODY(J,ILAST) - YDIFF(J)
               ZDIFF(J) = ZBODY(J,ILAST) - ZDIFF(J)
            END DO

            IF (NO_DECAY) THEN
               IEND = NFSTN
            ELSE
               IEND = IBOD3 - 1
            END IF

            P = ONE

            DO I = ILAST + 1, IEND

               IF (NFJ(I) /= NFP) GO TO 910

               IF (.NOT. NO_DECAY) THEN
                  P = (XF(IBOD3) - XF(I)) / (XF(IBOD3) - XF(ILAST))
               END IF

               DO J = 1, NFP
                  YBODY(J,I) = YBODY(J,I) + P * YDIFF(J)
                  ZBODY(J,I) = ZBODY(J,I) + P * ZDIFF(J)
               END DO

            END DO

         END IF

         DEALLOCATE (YDIFF, ZDIFF)

      END IF


      RETURN

  900 WRITE (IWRIT, '(/, A)')
     >   ' WBFUDGE: Sections 1:KWING2 must have equal pt. counts.'
      GO TO 999

  910 WRITE (IWRIT, '(/, A)')
     >   ' WBFUDGE: Sections ILAST:IBOD3 must have equal pt. counts.'

  999 STOP

      END SUBROUTINE WBFUDGE

C***********************************************************************
C
      SUBROUTINE WINGREAD (IGEOM, IW, NNWING, NWSURF, MXIWING, MXKWING,
     >                     IROLL, IFLAP, IDIHED,
     >                     NWSEC, NLG, NUG, NWG, ALG,
     >                     XLEADI, YLEADI, ZLEADI, ROUNDED, ROLL,
     >                     DFLAP, DSLAT, XFLAP, XSLAT, SWFLAP, SWSLAT,
     >                     AINIT, CINIT, TINIT, XINIT, YINIT, ZINIT,
     >                     XC1, YC1, ZC1, XC2, YC2, ZC2,
     >                     DIHED, DIAXIS, IQUAD, CLOCK, TITLE, IWRIT)
C
C     WINGREAD reads the defining sections for one wing-type geometry,
C     including a pylon, tail, or canard surface, or possibly a nacelle,
C     which is handled differently according to component number IW.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   IGEOM, IW, NNWING, NWSURF, MXIWING, MXKWING

      INTEGER, INTENT (IN), DIMENSION (NWSURF) ::
     >   IROLL, IFLAP, IDIHED

      INTEGER, INTENT (OUT) ::
     >   NWSEC(NWSURF)

      INTEGER, INTENT (OUT), DIMENSION (MXKWING,NWSURF) ::
     >   NLG, NUG, NWG

      REAL, INTENT (OUT), DIMENSION (MXKWING,NWSURF) ::
     >   ALG, XLEADI, YLEADI, ZLEADI, ROUNDED, ROLL,
     >   DFLAP, DSLAT, XFLAP, XSLAT, SWFLAP, SWSLAT,
     >   AINIT, CINIT, TINIT

      REAL, INTENT (OUT), DIMENSION (MXIWING,MXKWING,NWSURF) ::
     >   XINIT, YINIT, ZINIT

      REAL, INTENT (OUT), DIMENSION (NWSURF) ::
     >   XC1, YC1, ZC1, XC2, YC2, ZC2

      REAL, INTENT (OUT), DIMENSION (MXKWING,NWSURF) ::
     >   DIHED

      REAL, INTENT (OUT), DIMENSION (3,2,MXKWING,NWSURF) ::
     >   DIAXIS

      INTEGER, INTENT (OUT) ::
     >   IQUAD(2,NWSURF)

      REAL, INTENT (OUT) ::
     >   CLOCK(NWSURF)

      CHARACTER, INTENT (OUT) ::
     >   TITLE * 80

      INTEGER, INTENT (IN) ::
     >   IWRIT

C     Local constants:

      REAL, PARAMETER :: DELTA = 0.0000001, ONE = 1., ZERO = 0.

C     Local variables:

      INTEGER
     >   I, J, K, L, N, NL, NU

      REAL
     >   AL, FINTR, FNL, FNU, FNWSEC, FQUAD, FSEC, TCHORD, THKSCALE,
     >   XLD, YLD, ZLD, XYZSCALE, XL, YL, ZL, YSYM

      REAL, DIMENSION (MXIWING) ::
     >   XP, YP, ZP, XWG, YWG, ZWG

      LOGICAL
     >   DEFLECTIONS, DIHEDRAL, NOTNACELLE

C     Execution:

      READ (IGEOM, '(A)') TITLE
      READ (IGEOM,*)

      IF (IW <= NNWING) THEN

         READ (IGEOM,*) FNWSEC ! ALPHAW is no longer used

      ELSE ! Nacelle defining section

         READ (IGEOM,*) FNWSEC, FQUAD, FINTR, CLOCK(IW)
         IF (FQUAD == FINTR) THEN
            WRITE (IWRIT, '(/, A, A)') ' WINGREAD: Nacelle sections ',
     >         'must not start/end in the intersection neighborhood.'
            GO TO 900
         END IF
         IQUAD(1,IW) = FQUAD
         IQUAD(2,IW) = FINTR
         READ (IGEOM,*)
         READ (IGEOM,*) XC1(IW), YC1(IW), ZC1(IW),
     >                  XC2(IW), YC2(IW), ZC2(IW)
         READ (IGEOM,*)
         READ (IGEOM,*) ! Thrust vector coordinates, not used stand-alone
      END IF

      NWSEC(IW) = FNWSEC

      NOTNACELLE = IW <= NNWING ! Avoid testing undefineds for nacelles

      IF (NOTNACELLE) THEN
         DEFLECTIONS = IFLAP(IW) /= 0
         DIHEDRAL    = IDIHED(IW) /= 0
      ELSE
         DEFLECTIONS = .FALSE.
         DIHEDRAL    = .FALSE.
      END IF

      IF (NWSEC(IW) > MXKWING .OR.
     >   (NWSEC(IW) < 2 .AND. NOTNACELLE) .OR.
     >   (NWSEC(IW) < 1 .AND. .NOT. NOTNACELLE)) THEN
         WRITE (IWRIT, '(/, A, I4, A, I4)')
     >      ' WINGREAD:  Bad input no. of sections:', NWSEC(IW),
     >      '   Limit:', MXKWING
         GO TO 900
      END IF

      DO K = 1, NWSEC(IW)

         READ (IGEOM,*)
         IF (IROLL(IW) == 0) THEN
            READ (IGEOM,*) ZL, XL, YL, XYZSCALE, THKSCALE, AL, FSEC,
     >                     ROUNDED(K,IW)
            ROLL(K,IW) = ZERO
         ELSE
            READ (IGEOM,*) ZL, XL, YL, XYZSCALE, THKSCALE, AL, FSEC,
     >                     ROUNDED(K,IW), ROLL(K,IW)
         END IF

C        *** Temporary kludge to test pylon non-intersection (Beech) ***

         IF (NWSEC(IW) == 2) THEN
C ***       IF (ABS (XYZSCALE - 74.74) < 0.001) XYZSCALE = 80.
            IF (ABS (XYZSCALE - 74.74) < 0.001) XYZSCALE = 70.
         END IF

         IF (IW <= NNWING) THEN
            XLEADI(K,IW) = XL
            YLEADI(K,IW) = YL
            ZLEADI(K,IW) = ZL
         ELSE
            XLEADI(K,IW) = ZL ! XLE (cylindrical coordinates)
            YLEADI(K,IW) = XL ! RLE
            ZLEADI(K,IW) = YL ! TLE (fraction of 360 deg.)
         END IF

         AINIT(K,IW) = AL
         ALG(K,IW)   = AL
         CINIT(K,IW) = XYZSCALE
         TINIT(K,IW) = THKSCALE

         IF (FSEC /= ZERO) THEN
            READ (IGEOM,*)
            READ (IGEOM,*) YSYM, FNU, FNL

            IF (YSYM /= ZERO) THEN
               WRITE (IWRIT, '(/, A)')
     >            ' WINGREAD: YSYM option is not supported.'
               GO TO 900
            END IF

            IF (.NOT. DEFLECTIONS) THEN
               IF (NOTNACELLE) THEN
                  DSLAT(K,IW)  = ZERO
                  DFLAP(K,IW)  = ZERO
                  SWSLAT(K,IW) = ZERO
                  SWSLAT(K,IW) = ZERO
               END IF
            ELSE
               READ (IGEOM,*)
               READ (IGEOM,*) DSLAT(K,IW), XSLAT(K,IW), SWSLAT(K,IW),
     >                        DFLAP(K,IW), XFLAP(K,IW), SWFLAP(K,IW)
            END IF

            IF (DIHEDRAL) THEN
               READ (IGEOM,*)
               READ (IGEOM,*) DIHED(K,IW), DIAXIS(1:3,1:2,K,IW)
            ELSE
               DIHED(K,IW) = ZERO
               DIAXIS(1:3,1:2,K,IW) = ZERO
            END IF

            NU = FNU
            NL = FNL
            N  = NU + NL - 1

            IF (N > MXIWING) THEN
               WRITE (IWRIT, '(/, A, I4, A, I4)')
     >          ' WINGREAD: Too many pts. at section', K,'. Limit:',
     >             MXIWING
               GO TO 900
            END IF

            READ (IGEOM,*)
         END IF

         NLG(K,IW) = NL
         NUG(K,IW) = NU
         NWG(K,IW) = N

         IF (FSEC /= ZERO) THEN

            DO I = NL, N
               READ (IGEOM,*) XP(I), YP(I)
               XINIT(I,K,IW) = XP(I)
               YINIT(I,K,IW) = YP(I)
               ZINIT(I,K,IW) = ZERO
            END DO

            L = NL + 1
            READ (IGEOM,*)

            DO I = 1, NL
               J = L - I
               READ (IGEOM,*) XP(J), YP(J)
               XINIT(J,K,IW) = XP(J)
               YINIT(J,K,IW) = YP(J)
               ZINIT(J,K,IW) = ZERO
            END DO

         ELSE

            XINIT(1:N,K,IW) = XINIT(1:N,K-1,IW)
            YINIT(1:N,K,IW) = YINIT(1:N,K-1,IW)
            ZINIT(1:N,K,IW) = ZINIT(1:N,K-1,IW)

         END IF

C        Denormalize and apply any twist; find the true leading edge;
C        untwist about that, and renormalize precisely to ensure safe
C        application of shape functions on [0., 1.].  The normalized
C        and untwisted section may also be reused for the next K.

         CALL WINGSEC () ! Internal procedure

      END DO

      RETURN

  900 WRITE (IWRIT, '(/, A, I4)') ' Aborting at wing/nac. surface #', IW

      STOP

C     Internal procedure for WINGREAD:

      CONTAINS

!***********************************************************************
!
      SUBROUTINE WINGSEC ()
!
!     WINGSEC organizes an input wing section into a strictly untwisted,
!     normalized section, guaranteed to have X on [0, 1] by adjusting
!     the leading edge & chord if necessary.
!
!     03/25/97  D.Saunders  Introduced to simplify WING in CH_GRID.
!     07/19/98      "       Internal procedure form.
!     08/25/99      "       Added to AEROSURF for proper shape fn. treatment.
!
!***********************************************************************

!     Local constants:

      REAL, PARAMETER ::
     >   DEG2RAD = 0.0174532925, RAD2DEG = 57.29577951,
     >   HALF = 0.5, ONE = 1., ZERO = 0.

!     Local variables:

      INTEGER
     >   I, ILW, NTOTK, NLOWK
      REAL
     >   ALPHA, CA, DSQ, DSQMAX, SA, SCALE, THICK, XLW, YLW, ZLW,
     >   XTRAIL, YTRAIL, XXW, YYW, ZZW

!     Execution:

!     Apply section twists about the GIVEN leading edge,
!     and denormalize (as in PERTURB):

      ALPHA = AINIT(K,IW)
      THICK = TINIT(K,IW)
      SCALE = CINIT(K,IW)
      NTOTK = NWG(K,IW)
      NLOWK = NLG(K,IW)
      XXW = XINIT(NLOWK,K,IW)
      YYW = YINIT(NLOWK,K,IW)
      ZZW = ZINIT(NLOWK,K,IW)
      XLW = XLEADI(K,IW) + XXW * SCALE
      YLW = YLEADI(K,IW) + YYW * SCALE * THICK
      ZLW = ZLEADI(K,IW)
      CA  = COS (ALPHA*DEG2RAD)
      SA  = SIN (ALPHA*DEG2RAD)

      DO I = 1, NTOTK
         XWG(I) = XLW + SCALE * ((XINIT(I,K,IW) - XXW)*CA +
     >                           (YINIT(I,K,IW) - YYW)*SA * THICK)
         YWG(I) = YLW + SCALE * ((YINIT(I,K,IW) - YYW)*CA * THICK -
     >                           (XINIT(I,K,IW) - XXW)*SA)
         ZWG(I) = ZLW + SCALE *  (ZINIT(I,K,IW) - ZZW)
      END DO

!     Find the "true" leading edge as the furthest point from the
!     trailing edge in (X,Y) space:

      XTRAIL = (XWG(1) + XWG(NTOTK)) * HALF
      YTRAIL = (YWG(1) + YWG(NTOTK)) * HALF
      DSQMAX = ZERO

      DO I = 2, NTOTK - 1
         DSQ = (XWG(I) - XTRAIL)**2 + (YWG(I) - YTRAIL)**2
         IF (DSQ > DSQMAX) THEN
            DSQMAX = DSQ
            ILW = I
         END IF
      END DO

!     Adjust the leading edge, chord, and local twist:

      NLG(K,IW) = ILW
      NUG(K,IW) = NTOTK - ILW + 1

      XLW = XWG(ILW)
      YLW = YWG(ILW)
      ZLW = ZWG(ILW)
      XLEADI(K,IW) = XLW
      YLEADI(K,IW) = YLW
      ZLEADI(K,IW) = ZLW

      ALPHA    = ATAN2 (YLW - YTRAIL, XTRAIL - XLW)
      AINIT(K,IW) = ALPHA * RAD2DEG
      ALG(K,IW)   = AINIT(K,IW)

!     Untwist about the true leading edge, and remove the thickness
!     scaling (since the nominal section may be reused, and PERTURB
!     applies the scaling as a design variable):

      CA = COS (-ALPHA)
      SA = SIN (-ALPHA)
      THICK = ONE / THICK

      DO I = 1, NTOTK
         XINIT(I,K,IW) =  XLW + (XWG(I) - XLW)*CA + (YWG(I) - YLW)*SA
         YINIT(I,K,IW) = (YLW + (YWG(I) - YLW)*CA - (XWG(I) - XLW)*SA)
     >                   * THICK
         ZINIT(I,K,IW) = ZWG(I)
      END DO

      CINIT(K,IW) = MAX (XINIT(1,K,IW), XINIT(NTOTK,K,IW)) - XLW

!     Normalize:

      SCALE = ONE / CINIT(K,IW)
      YLW = YINIT(ILW,K,IW)

      DO I = 1, NTOTK
         XINIT(I,K,IW) = (XINIT(I,K,IW) - XLW) * SCALE
         YINIT(I,K,IW) = (YINIT(I,K,IW) - YLW) * SCALE
         ZINIT(I,K,IW) = (ZINIT(I,K,IW) - ZLW) * SCALE
      END DO

      END SUBROUTINE WINGSEC

      END SUBROUTINE WINGREAD

C***********************************************************************
C
      SUBROUTINE WINGWRITE (LUN, IW, WINGOUT)
C
C     Write the indicated wing-type component in AEROSURF geometry format.
C
C     08/24/99  DAS  Part of introducing shape functions.
C     08/26/99   "   Introduced NLGI/NUGI/NWGI in case of flap deflections.
C     06/23/00   "   Introduced dihedral.  Somehow, ROLL was missing.
C
C***********************************************************************

C     Global variables:

      USE GEOM1

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   LUN,               ! Unit number to write to
     >   IW,                ! Wing/nacelle component number
     >   WINGOUT            ! 0 = sections are untwisted and normalized;
                            ! 1 =    "      "    twisted and normalized;
                            ! 2 =    "      "    twisted and denormalized
C     Local constants:

      REAL, PARAMETER ::
     >   DEG2RAD = 0.0174532925, ONE = 1., ZERO = 0.

C     Local variables:

      INTEGER
     >   I, J, K, NLOWK, NSEC, NTOTK

      REAL
     >   CA, SA, FNFP, FNU, FNL, SCALE, THICK, XI, XLW, YLW, XXW, YYW

      REAL, ALLOCATABLE, DIMENSION (:) ::
     >   XTEMP, YTEMP, ZTEMP

      LOGICAL, SAVE ::
     >   FIRST = .TRUE.

C     Execution:

      NSEC  = NWSEC(IW)
      NTOTK = MAXVAL (NWGI(1:NSEC,IW))

      ALLOCATE (XTEMP(NTOTK), YTEMP(NTOTK))

      IF (FIRST) THEN
         FIRST = .FALSE.
         WRITE (LUN, '(A)') 'Geometry from AEROSURF'
      END IF

      WRITE (LUN, '(A, /, 2I4)') 'NSEC  IW', NSEC, IW

      DO K = 1, NSEC

         FNU   = NUGI(K,IW) ! Section counts prior to any flaps/slats
         NLOWK = NLGI(K,IW)
         FNL   = NLOWK
         NTOTK = NWGI(K,IW)

         DO I = 1, NTOTK
            XTEMP(I) = XINIT(I,K,IW)
            YTEMP(I) = YINIT(I,K,IW)
         END DO

         WRITE (LUN, '(2A)') '         Z       XLE       YLE',
     >      '     CHORD     THICK     TWIST       NEW RND  ROLL'

         IF (WINGOUT == 0) THEN ! Sections are untwisted and normalized

            WRITE (LUN, 1030) ZLEADI(K,IW), XLEADI(K,IW), YLEADI(K,IW),
     >                        CINIT(K,IW),  TINIT(K,IW),  AINIT(K,IW),
     >                        ONE, ROUNDED(K,IW), ROLL(K,IW)

         ELSE IF (WINGOUT == 1) THEN ! Twisted & normalized

            IF (AINIT(K,IW) /= ZERO) THEN

               SA = AINIT(K,IW) * DEG2RAD
               CA = COS (SA)
               SA = SIN (SA)

               DO I = 1, NTOTK ! Twist about [0.,0.]
                  XI = XTEMP(I)
                  XTEMP(I) = XI * CA + YTEMP(I) * SA
                  YTEMP(I) = YTEMP(I) * CA - XI * SA
               END DO

            END IF

            WRITE (LUN, 1030) ZLEADI(K,IW), XLEADI(K,IW), YLEADI(K,IW),
     >                        CINIT(K,IW),  TINIT(K,IW),  ZERO, ONE,
     >                        ROUNDED(K,IW), ROLL(K,IW)

         ELSE ! WINGOUT = 2  (Fully denormalized except for flap deflections
                             !and dihedral)
            THICK = TINIT(K,IW)
            SCALE = CINIT(K,IW)
            SA    = AINIT(K,IW) * DEG2RAD
            CA    = COS (SA)
            SA    = SIN (SA)
            XXW   = XTEMP(NLOWK)
            YYW   = YTEMP(NLOWK)
            XLW   = XLEADI(K,IW) + XXW * SCALE
            YLW   = YLEADI(K,IW) + YYW * SCALE * THICK

            DO I = 1, NTOTK
               XI = XTEMP(I)
               XTEMP(I) = XLW + SCALE*((XI - XXW)*CA +
     >                                 (YTEMP(I) - YYW)*SA * THICK)
               YTEMP(I) = YLW + SCALE*((YTEMP(I) - YYW)*CA * THICK -
     >                                 (XI - XXW)*SA)
            END DO

            WRITE (LUN, 1030) ZLEADI(K,IW), ZERO, ZERO, ONE, ONE, ZERO,
     >                        ONE, ROUNDED(K,IW), ROLL(K,IW)
         END IF

         WRITE (LUN, 1008) ZERO, FNU, FNL
         WRITE (LUN, 1009) DSLAT(K,IW), XSLAT(K,IW), SWSLAT(K,IW),
     >                     DFLAP(K,IW), XFLAP(K,IW), SWFLAP(K,IW)
         WRITE (LUN, 1010) DIHED(K,IW), DIAXIS(1:3,1:2,K,IW)
         WRITE (LUN, 1011)
         WRITE (LUN, 1004) (XTEMP(I), YTEMP(I), I = NLOWK, NTOTK)
         WRITE (LUN, 1012)
         WRITE (LUN, 1004) (XTEMP(I), YTEMP(I), I = NLOWK, 1, -1)

      END DO

      DEALLOCATE (XTEMP, YTEMP)

      RETURN

 1004 FORMAT (2F12.8)
 1007 FORMAT (2F15.5)
 1008 FORMAT (' YSYM   NU   NL ', /, 3F5.0)
 1009 FORMAT ('     DSLAT       XSLAT    SWSLAT',
     >        '     DFLAP       XFLAP    SWFLAP', /,
     >        2(F10.5, F12.5, F10.5))
 1010 FORMAT ('     DIHED',
     >        '          PX        PY        PZ',
     >        '          QX        QY        QZ', /,
     >        F10.5, 2(F12.5, 2F10.5))
 1011 FORMAT ('Upper surface')
 1012 FORMAT ('Lower surface')
 1030 FORMAT (7F10.4, F3.0, F8.4)

      END SUBROUTINE WINGWRITE

C***********************************************************************
C
      SUBROUTINE WINTERP (NNINSERT, NNWING, NNNAC,
     >                    NWSURF, MXIWING, MXKWING, NWSEC, IDEGREE,
     >                    CLOCK, ROLL, XLEADI, YLEADI, ZLEADI,
     >                    XINIT, YINIT, ZINIT, AINIT, CINIT, TINIT, ALG,
     >                    NLG, NUG, NWG, XWG, YWG, ZWG, ROUNDED,
     >                    IKBD, ICRT, ISAVE)
C
C     Insert defining sections interactively, saving them interleaved
C     with input sections as we go.  The raw geometry must be regular
C     and suited to spanwise lofting.  Nacelles are also treated.
C     The inputs are denormalized, real-space outputs from PERTURB;
C     the outputs are normalized.
C
C     12/12/98  David Saunders  Initial implementation for wings.
C     12/14/98    "      "      Handled nacelles as well.
C     01/18/01    "      "      Original stations were having their
C                               local twists zeroed out explicitly by
C                               mistake.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   NNINSERT, NNWING, NNNAC, NWSURF, MXIWING, MXKWING

      INTEGER, INTENT (IN), DIMENSION (NWSURF) ::
     >   NWSEC, IDEGREE

      REAL, INTENT (IN) ::
     >   CLOCK(NWSURF)

      REAL, INTENT (IN), DIMENSION (MXKWING,NWSURF) ::
     >   ROLL, XLEADI, YLEADI

      REAL, INTENT (INOUT), DIMENSION (MXKWING,NWSURF) ::
     >   ZLEADI ! Nacelle case modifies ZLEADI

      REAL, INTENT (INOUT), DIMENSION (MXIWING,MXKWING,NWSURF) ::
     >   XINIT, YINIT

      REAL, INTENT (IN), DIMENSION (MXIWING,MXKWING,NWSURF) ::
     >   ZINIT

      REAL, INTENT (IN), DIMENSION (MXKWING,NWSURF) ::
     >   AINIT, CINIT, TINIT, ALG

      INTEGER, INTENT (IN), DIMENSION (MXKWING,NWSURF) ::
     >   NLG, NUG, NWG

      REAL, INTENT (IN), DIMENSION (MXIWING,MXKWING,NWSURF) ::
     >   XWG, YWG, ZWG

      REAL, INTENT (IN), DIMENSION (MXKWING,NWSURF) ::
     >   ROUNDED

      INTEGER, INTENT (IN) ::
     >   IKBD, ICRT, ISAVE

C     Local constants:

      REAL,    PARAMETER ::
     >   HALF = 0.5, ONE = 1.0, RAD2DEG = 57.29577951, TWOPI = 360.,
     >   ZERO = 0.0, ZFLAG = 9.E7
      LOGICAL, PARAMETER ::
     >   NEW = .TRUE.

C     Local variables:

      REAL, DIMENSION (:), ALLOCATABLE ::
     >   XI, YI, ZI, XNEW, YNEW, ZNEW
      INTEGER
     >   I, IOS, IW, K, L, M, NL, NNEW, NSEC, NTOT
      REAL
     >   ANEW, CNEW, DERIV, FNL, FNU, ROLLNEW, ROUNDNEW,
     >   XLENEW, YLENEW, XTRAIL, YTRAIL, ZINT, ZLEAD
      LOGICAL
     >   NACELLE
      CHARACTER
     >   METHOD * 1

C     Execution:

      IW   = NNINSERT
      NSEC = NWSEC(IW)
      NL   = NLG(1,IW)
      FNL  = NL
      FNU  = NUG(1,IW)
      NTOT = NWG(1,IW)
      NACELLE = IW > NNWING

C     Check for feasibility:

      DO K = 2, NSEC
         IF (NLG(K,IW) /= NL .OR. NWG(K,IW) /= NTOT) THEN
            WRITE (ICRT, '(/, A, 3I5)')
     >         ' WINTERP: Mismatched point count.  K, NL, NTOT:',
     >         K, NLG(K,IW), NWG(K,IW)
            GO TO 99
         END IF
      END DO

      IF (.NOT. NACELLE) THEN
         WRITE (ICRT, '(/, A, /, (6F12.4))') ' Current Z stations:',
     >      ZLEADI(1:NSEC,IW)
      ELSE
         ZLEADI(1:NSEC,IW) = ZLEADI(1:NSEC,IW) * TWOPI

         WRITE (ICRT, '(/, A, /, (6F12.4))')
     >      ' Current nacelle section angular coordinates (degrees):',
     >      ZLEADI(1:NSEC,IW)
      END IF

      WRITE (ICRT, '(A)')

      IF (IDEGREE(IW) == 1) THEN
         METHOD = 'L'
      ELSE
         METHOD = 'B'
      END IF

      OPEN (UNIT=ISAVE, FILE='aero.sections', STATUS='UNKNOWN')

      ALLOCATE (XI(NSEC), YI(NSEC), ZI(NSEC),
     >          XNEW(NTOT), YNEW(NTOT), ZNEW(NTOT))

      L = 1 ! Input section pointer
      M = 0 ! Output   "    "
      NNEW = 0

      DO ! Until EOF at the prompt

         IF (.NOT. NACELLE) THEN
            WRITE(ICRT, '(A)', ADVANCE='NO')
     >         ' Z at which to interpolate? '
         ELSE
            WRITE(ICRT, '(A)', ADVANCE='NO')
     >         ' Angle at which to interpolate? '
         END IF
         READ (IKBD, *, IOSTAT=IOS) ZINT

         IF (IOS < 0) ZINT = ZFLAG ! Force transfer of remaining sections

C        Transfer any untransferred original sections preceding this one:


         DO K = L, NSEC
            ZLEAD = ZLEADI(K,IW)

            IF (ZLEAD >= ZINT) EXIT

            M = K

            IF (.NOT. NACELLE) THEN
               WRITE (ISAVE, 100)
               WRITE (ISAVE, 101)
     >            ZLEAD, XLEADI(K,IW), YLEADI(K,IW), CINIT(K,IW),
     >            ONE, AINIT(K,IW), ONE, ROUNDED(K,IW), ROLL(K,IW)
            ELSE
               ZLEAD = ZLEAD / TWOPI
               WRITE (ISAVE, 105)
               WRITE (ISAVE, 101)
     >            XLEADI(K,IW), YLEADI(K,IW), ZLEAD, CINIT(K,IW),
     >            ONE, AINIT(K,IW), ONE, ROUNDED(K,IW)
            END IF

            WRITE (ISAVE, 102) ZERO, FNU, FNL
            WRITE (ISAVE, 103) (XINIT(I,K,IW), YINIT(I,K,IW),
     >         I = NL, NTOT)
            WRITE (ISAVE, 104)
            WRITE (ISAVE, 103) (XINIT(I,K,IW), YINIT(I,K,IW),
     >         I = NL, 1, -1)

         END DO

         IF (ZINT == ZFLAG) EXIT ! Done

         L = M + 1
         NNEW = NNEW + 1

         DO I = 1, NTOT

            XI(1:NSEC) = XWG(I,1:NSEC,IW)
            YI(1:NSEC) = YWG(I,1:NSEC,IW)
            ZI(1:NSEC) = ZWG(I,1:NSEC,IW)

            CALL LCSFIT (NSEC, ZLEADI(1,IW), XI, NEW, METHOD,
     >                   1, ZINT, XNEW(I), DERIV)

            CALL LCSFIT (NSEC, ZLEADI(1,IW), YI, NEW, METHOD,
     >                   1, ZINT, YNEW(I), DERIV)

            CALL LCSFIT (NSEC, ZLEADI(1,IW), ZI, NEW, METHOD,
     >                   1, ZINT, ZNEW(I), DERIV)
         END DO

         CALL LCSFIT (NSEC, ZLEADI(1,IW), XLEADI(1,IW), NEW, METHOD,
     >                1, ZINT, XLENEW, DERIV)

         CALL LCSFIT (NSEC, ZLEADI(1,IW), YLEADI(1,IW), NEW, METHOD,
     >                1, ZINT, YLENEW, DERIV)

         CALL LCSFIT (NSEC, ZLEADI(1,IW),  CINIT(1,IW), NEW, METHOD,
     >                1, ZINT, CNEW, DERIV)

         CALL LCSFIT (NSEC, ZLEADI(1,IW),   ROLL(1,IW), NEW, METHOD,
     >                1, ZINT, ROLLNEW, DERIV)

         IF (ROLLNEW /= ZERO) THEN ! Undo any roll put in by PERTURB

            CALL ROTATE2D (NTOT, ZNEW, YNEW, ROLLNEW, ZINT, YLENEW)

         END IF

         IF (.NOT. NACELLE) THEN

C           Local twist is a property of the new section - not a simple
C           function of neighboring twists:

            XTRAIL = (XNEW(1) + XNEW(NTOT)) * HALF
            YTRAIL = (YNEW(1) + YNEW(NTOT)) * HALF
            ANEW   = ATAN2 (YLENEW - YTRAIL, XTRAIL - XLENEW)

C           Untwist about the leading edge:

            IF (ANEW /= ZERO) THEN
               ANEW   = ANEW * RAD2DEG

               CALL ROTATE2D (NTOT, XNEW, YNEW, ANEW, XLENEW, YLENEW)

               XTRAIL = (XNEW(1) + XNEW(NTOT)) * HALF
               CNEW   = XTRAIL - XLENEW
            END IF

         ELSE
            ANEW = ZERO
         END IF

         DO I = 1, NTOT ! Renormalize
            XNEW(I) = (XNEW(I) - XLENEW) / CNEW
            YNEW(I) = (YNEW(I) - YLENEW) / CNEW
         END DO

         IF (.NOT. NACELLE) THEN ! Else can't be sure ...
            XNEW(1)    = ONE
            XNEW(NTOT) = ONE
         END IF

         IF (.NOT. NACELLE) THEN
            WRITE (ISAVE, 100)
            WRITE (ISAVE, 101) ZINT, XLENEW, YLENEW, CNEW,
     >         ONE, ANEW, ONE, ROUNDED(K,IW), ROLLNEW
         ELSE
            ZINT = ZINT / TWOPI
            WRITE (ISAVE, 105)
            WRITE (ISAVE, 101) XLENEW, YLENEW, ZINT, CNEW,
     >         ONE, ZERO, ONE, ROUNDED(K,IW)
         END IF
         WRITE (ISAVE, 102) ZERO, FNU, FNL
         WRITE (ISAVE, 103) (XNEW(I), YNEW(I), I = NL, NTOT)
         WRITE (ISAVE, 104)
         WRITE (ISAVE, 103) (XNEW(I), YNEW(I), I = NL, 1, -1)

      END DO ! Next new Z or angle

      DEALLOCATE (XI, YI, XNEW, YNEW)
      CLOSE (UNIT=ISAVE)

      NNEW = NSEC + NNEW
      WRITE (ICRT, '(/, I4, A)') NNEW,
     >   ' sections written to aero.sections.'

   99 CONTINUE

  100 FORMAT ('         Z       XLE       YLE     CHORD     THICK',
     >        '     TWIST NEW ROUND      ROLL')
  101 FORMAT (6F10.4, F4.0, F6.0, F10.4)
  102 FORMAT ('  YSYM    NU    NL', /, 3F6.0, /, ' Upper surface')
  103 FORMAT (2F10.6)
  104 FORMAT (' Lower surface')
  105 FORMAT ('       XLE       RLE       TLE     CHORD     THICK',
     >        '     TWIST NEW ROUND')

      END SUBROUTINE WINTERP

C***********************************************************************
C
      SUBROUTINE WNFUDGE (NCALL, NWSURF, IN, IW, ILWING, KLWING, NWSEC,
     >                    NIN, NKN, NIW, NKW, XNAC, YNAC, ZNAC,
     >                    XWNG, YWNG, ZWNG, XTRAP, ZINNER, ZOUTER)
C
C     WNFUDGE ensures that an HSCT nacelle lends itself to the scheme for
C     closing the diverter at the wing trailing edge at specified points.
C     It constructs a two-piece surface fore and aft of the wing trailing
C     edge between ZINNER and ZOUTER and interpolates for Y = Y(X,Z) at the
C     relevant nacelle section points.
C
C     10/26/98  DAS/JJR  Simpler approach after abandoning smooth warping.
C     01/30/99  DAS      X/Y/ZRG are packed now, so pass nacelle & wing.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   NCALL,                ! -3 means first call for this case
     >   NWSURF,               ! # "wing" + "nacelle" type surfaces
     >   IN,                   ! Surface # of the nacelle
     >   IW,                   ! Surface # of the wing
     >   ILWING(NWSURF),       ! Paneling dimensions specified for each
     >   KLWING(NWSURF),       ! wing-type or nacelle surface
     >   NWSEC(NWSURF),        ! Pre-spanwise-gridding number of wing sections
     >   NIN, NKN,             ! Dimensions of packed regularized nacelle data
     >   NIW, NKW              ! ... and packed wing

      REAL, INTENT (INOUT), DIMENSION (NIN,NKN) ::
     >   XNAC, YNAC, ZNAC      ! Regularized nacelle sections; some parts are
                               ! fudged upon return

      REAL, INTENT (IN), DIMENSION (NIW,NKW) ::
     >   XWNG, YWNG, ZWNG      ! Regularized wing sections

      REAL, INTENT (IN) ::
     >   XTRAP(NWSURF),        ! Controls used by SURFER's scheme for closing
     >   ZINNER(NWSURF),       ! the diverters at the wing trailing edge
     >   ZOUTER(NWSURF)

C     Local constants:

      INTEGER, PARAMETER ::
     >   NJPATCH = 33 ! Enough points in the artificial surface to capture the
                      ! wing trailing edge and fore/aft edges on the nacelle
      REAL, PARAMETER ::
     >   HALF = 0.5, ONE = 1., TOLER = 1.E-7

C     Local variables:

      INTEGER
     >   I, I1, I2, IER, ILE, IP, J, JP, K, K1, K2, K12, L, NI
      REAL
     >   DZ, P, Q, R, XLEFT, XMIDDLE, XRIGHT, ZJ,
     >   XPATCH(3,NJPATCH), YPATCH(3,NJPATCH), ZPATCH(3,NJPATCH),
     >   SPATCH(6+NJPATCH*2)

C     System functions:

      REAL, INTRINSIC :: EPSILON

C     Procedures:

      EXTERNAL RIPPLE2D

C     Storage:

      REAL, SAVE :: EPS

C     Execution:

      IF (NCALL == -3) EPS = MAX (5.* EPSILON (EPS), TOLER) ! For p/q in [0,1]

C     Locate the wing trailing edge points where the diverter is to be closed,
C     and set up the interior patch line along the trailing edge:

      DZ = (ZOUTER(IN) - ZINNER(IN)) / REAL (NJPATCH - 1)
      K1 = 1

      DO J = 1, NJPATCH
         ZJ = ZINNER(IN) + DZ * REAL (J - 1)

         DO K = K1, NWSEC(IW)
            IF (ZJ < ZWNG(1,K)) EXIT
         END DO

         K1 = K - 1
         R  = (ZJ - ZWNG(1,K1)) / (ZWNG(1,K1+1) - ZWNG(1,K1))
         ZPATCH(2,J) = ZJ
         YPATCH(2,J) = (ONE - R) * YWNG(1,K1) + R * YWNG(1,K1+1)
         XPATCH(2,J) = (ONE - R) * XWNG(1,K1) + R * XWNG(1,K1+1)
      END DO

C     I range of nacelle sections to fudge:

      XMIDDLE = (XPATCH(2,1) + XPATCH(2,NJPATCH)) * HALF
      XLEFT   = XMIDDLE + XTRAP(IN) ! XTRAP < 0
      XRIGHT  = XMIDDLE - XTRAP(IN)

      ILE = (ILWING(IN) + 1) / 2 + 2 ! Must be monotonic X from ILE on
      NI  = ILWING(IN) - ILE
      I1  = ILE / 2

      K12 = KLWING(IN) / 2 ! Roughly the 12 o'clock nacelle section number,
                           ! assuming nac. sections start at 6 o'c. & go clockw.

      CALL INTERVAL (NI, XNAC(ILE,K12), XLEFT, ONE, I1)

      I1  = ILE + I1 - 1 ! Nacelle section at the upstream edge of the fudge
      I2  = ILE + I1 / 2

      CALL INTERVAL (NI, XNAC(ILE,K12), XRIGHT, ONE, I2)

      I2  = ILE + I2     ! Nacelle section at the downstream edge ...

C     Upstream and downstream edges of the nacelle patch:

      L = 1
      DO I = I1, I2, I2 - I1

         DO K = K12, 1, -1
            IF (ZINNER(IN) > ZNAC(I,K)) EXIT
         END DO

         K1 = K

         DO K = K12, KLWING(IN) - 1 ! -1 because NCYLINDER sets up this many
            IF (ZOUTER(IN) < ZNAC(I,K)) EXIT
         END DO

         K2 = K

         DO J = 1, NJPATCH
            ZJ = ZINNER(IN) + DZ * REAL (J - 1)

            DO K = K1, K2
              IF (ZJ < ZNAC(I,K)) EXIT
            END DO

            K1 = K - 1
            R  = (ZJ - ZNAC(I,K1)) / (ZNAC(I,K1+1) - ZNAC(I,K1))
            ZPATCH(L,J) = ZJ
            YPATCH(L,J) = (ONE - R) * YNAC(I,K1) + R * YNAC(I,K1+1)
            XPATCH(L,J) = (ONE - R) * XNAC(I,K1) + R * XNAC(I,K1+1)
         END DO

         L = 3 ! Go to the downstream edge of the patch

      END DO

      ZPATCH(1:3,NJPATCH) = ZOUTER(IN) ! Exactly

C     Fill the patch interior in two pieces:

      CALL TFIQ3D (3, 1, 2, 1, NJPATCH, XPATCH, YPATCH, ZPATCH, SPATCH)
      CALL TFIQ3D (3, 2, 3, 1, NJPATCH, XPATCH, YPATCH, ZPATCH, SPATCH)

C     For each nacelle section to be fudged ...

      IP = 1

      DO I = I1 + 1, I2 - 1

         JP = 1

         DO K = K12, 1, -1
            IF (ZINNER(IN) > ZNAC(I,K)) EXIT
         END DO

         K1 = K + 1

         DO K = K12, KLWING(IN) - 1
            IF (ZOUTER(IN) < ZNAC(I,K)) EXIT
         END DO

         K2 = K - 1

         DO K = K1, K2

            CALL RIPPLE2D (3, NJPATCH, 1, 3, 1, NJPATCH, XPATCH, ZPATCH,
     >                     XNAC(I,K), ZNAC(I,K), IP, JP, EPS, P, Q,
     >                     IER)
            YNAC(I,K) =
     >       (ONE - Q)*((ONE - P)*YPATCH(IP,JP)   + P*YPATCH(IP+1,JP)) +
     >              Q *((ONE - P)*YPATCH(IP,JP+1) + P*YPATCH(IP+1,JP+1))

         END DO

      END DO

      END SUBROUTINE WNFUDGE

C***********************************************************************
C
      SUBROUTINE WREGUL (K1, K2, MXIWING, MXKWING, XWG, YWG, ZWG,
     >                   ROUNDED, NLG, NUG, NWG, ILE, NHALF, SNORM,
     >                   MXIGRID, MXKGRID, XRG, YRG, ZRG)
C
C     WREGUL regularizes a range of wing geometry sections by imposing
C     the indicated chordwise mesh along the upper and lower surfaces.
C
C     06/14/96  DAS  Extracted from SYN87'S WSURF for reuse by NLCON.
C     10/02/96   "   Rounded leading edges are handled wrap-around now.
C     09/25/97   "   Fully argument-driven, Fortran 90 version with
C                    application to multiple panels in mind.
C     09/23/98   "   SURFER now expects the regularized sections input
C                    in X/Y/ZBASE, forcing new dimensions for X/Y/ZRG.
C     01/30/99   "   This is no longer true, but the same code works.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   K1, K2,                    ! Subset of wing sections to treat
     >   MXIWING, MXKWING           ! Declared dimensions of wing sections

      REAL, INTENT (IN), DIMENSION (MXIWING,MXKWING) ::
     >   XWG, YWG, ZWG              ! Defining section coordinates

      REAL, INTENT (IN) ::
     >   ROUNDED(MXKWING)           ! 0. means sharp leading edge

      INTEGER, INTENT (IN), DIMENSION (MXKWING) ::
     >   NLG, NUG, NWG              ! Numbers of geometry pts. on lower & upper
                                    ! surface and on whole sections, wrap-around
      INTEGER, INTENT (IN) ::
     >   ILE,                       ! Leading edge index of regular mesh
     >   NHALF                      ! # mesh pts. for lower & upper halves

      REAL, INTENT (IN) ::
     >   SNORM(NHALF)               ! Normalized chordwise distribution

      INTEGER, INTENT (IN) ::
     >   MXIGRID, MXKGRID           ! Declared dimensions of X/Y/ZRG(*,*)

      REAL, INTENT (OUT), DIMENSION (MXIGRID,MXKGRID) ::
     >   XRG, YRG, ZRG              ! Regularized sections

C     Local constants:

      REAL,      PARAMETER :: ONE = 1., ZERO = 0.
      LOGICAL,   PARAMETER :: NEW = .TRUE., NORM = .FALSE.
      CHARACTER, PARAMETER :: LOOSE * 1 = 'B', TIGHT * 1 = 'M'

C     Use "loose" fits except monotonic X vs. T is best at leading edge.

C     Local variables:

      INTEGER    I, ITL, K, N, NL, NU, NWING
      REAL       SLOWER, SMAX, SUPPER,
     >           DERIVS(2*NHALF), SEVAL(2*NHALF), SWG(MXIWING)

C     Execution:

      NWING = 2*NHALF - 1
      ITL   = ILE - NHALF + 1

      DO K = K1, K2

         NL = NLG(K)
         NU = NUG(K)
         N  = NWG(K)

         IF (ROUNDED(K) == ZERO) THEN  ! Not rounded - do surfaces separately

C           Unnormalized geometry arc lengths along the lower surface:

            CALL CHORDS3D (NL, XWG(1,K), YWG(1,K), ZWG(1,K), NORM,
     >                     SMAX, SWG)

C           Denormalized chordwise mesh distribution:

            DO I = 1, NHALF
               SEVAL(I) = SMAX * (ONE - SNORM(NHALF + 1 - I))
            END DO

C           Lower surface interpolation along the arc:

            CALL LCSFIT (NL, SWG, XWG(1,K), NEW, TIGHT, NHALF, SEVAL,
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

            CALL LCSFIT (NU, SWG, XWG(NL,K), NEW, TIGHT, NHALF, SEVAL,
     >                   XRG(ILE,K), DERIVS)
            CALL LCSFIT (NU, SWG, YWG(NL,K), NEW, LOOSE, NHALF, SEVAL,
     >                   YRG(ILE,K), DERIVS)
            CALL LCSFIT (NU, SWG, ZWG(NL,K), NEW, LOOSE, NHALF, SEVAL,
     >                   ZRG(ILE,K), DERIVS)

        ELSE  ! Rounded leading edge case:

            CALL CHORDS3D (N, XWG(1,K), YWG(1,K), ZWG(1,K), NORM,
     >                     SMAX, SWG)

            SLOWER = SWG(NL)
            SUPPER = SMAX - SLOWER

            DO I = 1, NHALF
               SEVAL(I) = SLOWER * (ONE - SNORM(NHALF + 1 - I))
               SEVAL(I + NHALF - 1) = SLOWER + SUPPER * SNORM(I)
            END DO

            CALL LCSFIT (N, SWG, XWG(1,K), NEW, TIGHT, NWING, SEVAL,
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
      SUBROUTINE XYZ_WRITE (IFORMAT, LUN, IDIM, JDIM, NI, NJ, X, Y, Z)
C
C     Write one block of a PLOT3D /mgrid-type (X,Y,Z) file.
C
C     01/30/99  DAS  Needed to deal with packed data.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   IFORMAT,             ! Even = unformatted; odd = formatted
     >   LUN,                 ! Logical unit
     >   IDIM, JDIM,          ! Block dimensions
     >   NI, NJ               ! Active portions

      REAL, INTENT (IN), DIMENSION(IDIM,JDIM) ::
     >   X, Y, Z

C     Local variables:

      INTEGER
     >   I, J

      LOGICAL
     >   FULLYPACKED

C     Execution:

      FULLYPACKED = NI == IDIM .AND. NJ == JDIM

      IF (MOD (IFORMAT, 2) == 0) THEN ! Unformatted

         IF (FULLYPACKED) THEN
            WRITE (LUN) X, Y, Z
         ELSE
            WRITE (LUN)
     >         ((X(I,J), I = 1, NI), J = 1, NJ),
     >         ((Y(I,J), I = 1, NI), J = 1, NJ),
     >         ((Z(I,J), I = 1, NI), J = 1, NJ)
         END IF

      ELSE ! Formatted

         IF (FULLYPACKED) THEN
            WRITE (LUN, '(6F12.6)') X, Y, Z
         ELSE
            WRITE (LUN, '(6F12.6)')
     >         ((X(I,J), I = 1, NI), J = 1, NJ),
     >         ((Y(I,J), I = 1, NI), J = 1, NJ),
     >         ((Z(I,J), I = 1, NI), J = 1, NJ)
         END IF
      END IF

      END SUBROUTINE XYZ_WRITE
