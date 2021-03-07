********************************************************************************
*
      PROGRAM HB_GRID
*
*        HB_GRID generates a 2-D structured grid for a hypersonic blunt body
*     defined by a discretized curve, not necessarily symmetric.  The outer
*     boundary is a function of the intended Mach number, angle of attack, and
*     surface radius at the stagnation point where the free stream is normal to
*     the surface.  (The effective body deflection angle and gamma also enter
*     into the outer boundary shape.  The hope is that the result will be an
*     improvement over automated hyperbolic gridding, which produces an outer
*     boundary much too far upstream than the desired boundary that should be
*     aligned with the shock.)
*
*        The downstream boundary is placed at a specified fraction of the body
*     length from the nose, as this may vary with Mach number.
*
*        The single-block topology is constructed in lower and upper halves.
*     The initial grid has Euler-type spacing.  This may optionally be
*     redistributed in the radial direction for Navier-Stokes purposes.
*
*        For axisymmetric bodies and ballistic entry (Alpha = 0), the upper
*     half of the two-halves grid is written as an additional output,
*     axisymmetric.gu, as suited to the DPLR2D flow solver.
*
*        Originally, the geometry had to be defined by a full centerline curve
*     starting at the lower trailing edge.  For this non-analytic definition,
*     the geometry format is airfoil-like, defining reference length:
*
*        Title ! Description
*        n     ! Number of points
*        x  y  ! Lower trailing edge point ...
*        x  y  ! ... proceeding forward then back in wrap-around fashion
*        :  :
*        :  :
*        x  y  ! Point n at the upper trailing edge
*
*        For the standard sphere-cone and spherical-section geometries, the
*     analytic definition can be specified and the generatrix is calculated
*     internally instead of being read.  Enter 'none' for geometry file to
*     invoke this case.  Then a namelist is looked for at the end of the
*     control file to define the sphere-cone more simply - see below.
*     A biconic forebody can also be specified similarly now.
*
*     Control inputs (standard input):
*
*                       ! HB_GRID control file
*        'abc.dat'      ! Input geometry file name; 'none' => analytic
*        'xyz.grid'     ! Output grid name: single-block PLOT3D /mgrid format
*        FORMATTED?     ! T = formatted output; F = unformatted output
*        TWO_D?         ! T = (x,y) only; F = (x,y,z) for Gridgen access
*        Mach           ! Intended Mach number
*        Alpha          ! Intended angle of attack
*        Gamma          ! Ratio of specific heats, affecting shock estimate
*        RFUDGE         ! For fudging body radius at stagnation point (1.5?)
*        BFUDGE         ! For fudging limiting shock angle, beta (2.0?)
*        AFUDGE         ! For reducing the Alpha-related rotation
*        THETA          ! Deflection angle represented by the body (degrees)
*        XC_DB          ! Cut-off location of downstream boundary (% length)
*        XC_SB          ! Distance of sonic bubble from stag. pt. ( "    " )
*        NI_EU          ! Dimensions of half of initial Euler-type grid; ...
*        NJ_EU          ! .. J is the radial direction
*        D1BODY         ! Surface spacing at stag. pt. & aft body (* arc l.) ...
*        D2BODY         ! D1BODY = 0. => curvature-based; D2BODY < 0. => 1-sided
*        D1SHOCK        ! Corresp. spacings on outer bndry. (fraction of ...
*        D2SHOCK        ! ... active lower bndry.); D1SHOCK < 0. => use body D1
*        NBLAYER        ! Euler value - see below - try 4
*        RBLAYER        !   "     "      "     "     "  1.1
*        D1NOSE         ! Euler-spacing of inner radial increment at nose ...
*        D1TAIL         ! ... and at downstream bndry. (fractn. of radial lgth.)
*        D2NOSE         ! Euler-spacing of outer radial increment at nose ...
*        D2TAIL         ! ... and at downstream bndry. (fractn. of radial lgth.)
*        SMOOTH?        ! T = perform elliptic smoothing of Euler-type grid
*        ITERATIONS?    ! T = show iterations for Newton solns. & smoothing
*        REGRID?        ! T = redistribute in radial direction using NJ_NS
*        NJ_NS          ! Radial dimension of N-S-type grid
*        NBLAYERNS      ! Try  20  for N-S grids
*        RBLAYERNS      !  "  1.05  "   "   "
*        D1NOSENS       ! N-S-spacing of inner radial increment at nose ...
*        D1TAILNS       ! ... and at downstrm. boundary (frn. of full body lgth)
*        D2NOSENS       ! N-S-spacing of outer radial increment at nose ...
*        D2TAILNS       ! ... and at downstrm. boundary (frn. of full body lgth)
*
*        In all cases, if grid spacing input D2* is positive, Vinokur-type
*        distributions are used.  However, the chore of estimating D2* may be
*        avoided by entering negative values, in which case one-sided stretching
*        is employed.  Curvature-based gridding along the wall (only) is also an
*        option via D1BODY = 0.
*
*        The D1* and D2* inputs are fractions of radial grid line lengths,
*        not body length, because of the unusually short lengths involved.
*
*        For the cases of D2* > 0:
*        Grid spacing on the body is 2-sided Vinokur, starting at the stag. pt.
*        Grid spacing for J = 1 : NBLAYER varies geometrically off the body
*        using RBLAYER as the multiplier; Vinokur distributions are used
*        for J = NBLAYER-1 : Jmax, with D2Jmax controlling the outer increment.
*
*     Optional elliptic smoothing inputs at the end of hb_grid.inp:
*
*        The defaults should make these inputs redundant.
*        Note that orthogonality is forced at the stagnation I boundary;
*        FGMODE = 'AAAN' applies to all other I boundaries by default.
*        See ellip2d.f for descriptions.
*        Here are the defaults (all lines starting in column 2):
*
*        $ELLIPTIC
*        PRINT2D = .TRUE.,  CONV2D  = 10.,    CONVMIN  = 0.,   DMAX2D = 0.,
*        ITMAX   = 2000,    ITFLOAT = 0,      ITFREEZE = 100,  JLAYER = 4,
*        FGLIMIT = 1.,      URFG    = 0.1,    URFLOAT  = 0.1,
*        OMG2D   = 1.1,     POWERI  = 1.0,    POWERJ   = 0.5,
*        BGMODE  = 'YYYY',  FGMODE  = 'AAAN', SPMODE   = 'YYYY',
*        EXPI1   = 0.20,    EXPI2   = 0.45,
*        EXPJ1   = 0.45,    EXPJ2   = 0.45,
*        $END
*
*     Optional analytic sphere-cone specification following $ELLIPTIC $END:
*
*        This namelist is read if the geometry file name is entered as 'none'.
*        A biconic forebody is distinguished from a sphere-cone or spherical
*        section if RADIUS_CONE_JUNCTURE is entered with a positive value, in
*        which case HALF_CONE_ANGLE_FORE and HALF_CONE_ANGLE_AFT should also
*        be entered.
*
*        Sample for a sphere-cone (each line starting in column 2):
*
*        $SPHERE_CONE_INPUTS
*        x_nose = 0.,            ! Nose (x,y) coordinates
*        r_nose = 0.,
*        radius_nose = 5.,
*        radius_base = 10.,
*        radius_shoulder = 0.25647,
*        half_cone_angle = 65.   ! 0 => spherical section (Apollo-like)
*        skirt_angle = 25.,
*        skirt_length = 0.1,     ! from aft-shoulder tangency pt. x to cut-off x
*        $END
*
*     Reference:  Hypersonic & High Temperature Gas Dynamics, John D. Anderson
*
*     10/16/02  DAS  Initial design, for airfoil leading edge studies.
*     10/23/02   "   Elliptic smoothing is ineffective near the nose;
*                    introduced further artificial boundaries for the TFI.
*     10/31/02   "   Greater NI_EU and RFUDGE help the smoothing;
*                    one-sided stretching options via D2* < 0.
*     11/13/02   "   Introduced D1/2SHOCK; switched D1/2BODY to be fractions of
*                    lower arc length, not active X length.
*     09/20/07   "   Application to MSL forebody prompted entry of deflection
*                    angle THETA, and calculation of BETA from the relation
*                    2.16 of Anderson (THETA-BETA-Mach relation), which also
*                    requires a Gamma input.  It doesn't seems to work, so
*     09/21/07       fall back on the THETA input.  Perhaps the theory applies
*                    more to slender blunt bodies than to very blunt capsules?
*     12/30/10   "   EXPDIS5 is preferable to EXPDIS4.
*     01/04/11   "   The outflow boundaries are no longer simple straight lines.
*                    They're now orthogonal to the body and the estimated shock.
*     01/05/11   "   Curvature-based surface grid option via D1BODY = 0.
*     01/06/11   "   Analytic sphere-cone option (namelist $SPHERE_CONE_INPUTS).
*     01/07/11   "   If Alpha = 0, save the upper portion of the grid as an
*                    additional output for axisymmetric flow calculations.
*     01/14/11   "   Also, ensure y = 0 exactly along the stagnation point line,
*                    prior to the elliptic smoothing (not after).
*                    ITMAX = 200 was too low: 1000 is better for the default.
*     03/07/11   "   Changed Eqs. 5.37, 5.38 (Anderson, p. 189) from the
*                    cylinder-wedge variants to the sphere-cone variants;
*                    raised default ITMAX to 2000 for dense grid cases.
*     10/26/11   "   Added the analytic biconic option; defaulted ITMAX to 5000.
*     11/30/11   "   The biconic option now has a way of rounding the vertex
*                    where the two cones meet (radius_vertex > 0.).
*
*     David Saunders, ELORET Corp./NASA Ames Research Center, Moffett Field, CA
*                    Now ERC, Inc./NASA ARC
*
********************************************************************************

      USE BICONIC_PARAMETERS      ! For derived data type biconic_type
      USE SPHERE_CONE_PARAMETERS  ! For derived data type sphere_cone_type
      USE SPHERE_CONE             ! For constructing a generatrix

      IMPLICIT NONE

*     Constants:

      INTEGER, PARAMETER ::
     >   LUNGEOM  = 1,     ! Input geometry
     >   LUNSHOCK = 2,     ! Bow shock estimated shape for possible inspection
     >   LUNGRID  = 3,     ! Output grid
     >   LUNINP   = 5,     ! Control file
     >   LUNCRT   = 6      ! Screen

      REAL, PARAMETER ::
     >   ZERO = 0.

*     Local variables except ELLIP2D controls:

      INTEGER ::
     >   I, IER, IOS, J, LINE,
     >   NBLAYER, NBLAYERNS, NGEOM, NI, NI_EU, NJ_EU, NJ_NS

      REAL ::
     >   AFUDGE, ALPHA, BFUDGE, BODY_LENGTH, CHORD,
     >   D1BODY, D1NOSE, D1NOSENS, D1SHOCK, D1TAIL, D1TAILNS,
     >   D2BODY, D2NOSE, D2NOSENS, D2SHOCK, D2TAIL, D2TAILNS,
     >   FMACH, GAMMA, R, RADIUS, RBLAYER, RBLAYERNS, RFUDGE, THETA,
     >   XC_DB, XC_SB

      REAL, ALLOCATABLE, DIMENSION (:) ::
     >   XGEOM, YGEOM

      REAL, ALLOCATABLE, DIMENSION (:,:) ::
     >   X, Y, Z, XNS, YNS, ZNS

      LOGICAL ::
     >   ANALYTIC, FORMATTED, ITERATIONS, REGRID, SPHERECONE, SMOOTH,
     >   TWO_D

      CHARACTER ::
     >   FILENAME * 80

      TYPE (BICONIC_TYPE) ::
     >   BICONIC

      TYPE (SPHERE_CONE_TYPE) ::
     >   CAPSULE

*     (Most of the) arguments used by ELLIP2D:

      INTEGER ::
     >   ITMAX, ITFLOAT, ITFREEZE, JLAYER

      REAL ::
     >   CONV2D, CONVMIN, DMAX2D, OMG2D,
     >   POWERI, POWERJ, EXPI1, EXPI2, EXPJ1, EXPJ2,
     >   FGLIMIT, URFG, URFLOAT

      CHARACTER ::
     >   BGMODE*4,  FGMODE*4, SPMODE*4

      NAMELIST /ELLIPTIC/
     >   ITMAX, ITFLOAT, ITFREEZE, JLAYER,
     >   CONV2D, CONVMIN, DMAX2D, OMG2D,
     >   BGMODE, FGMODE, SPMODE,
     >   POWERI, POWERJ, EXPI1, EXPI2, EXPJ1, EXPJ2,
     >   FGLIMIT, URFG, URFLOAT

      REAL ::
     >   X_NOSE, R_NOSE, RADIUS_NOSE, RADIUS_CONE_JUNCTURE, RADIUS_BASE,
     >   RADIUS_SHOULDER, RADIUS_VERTEX, HALF_CONE_ANGLE,
     >   HALF_CONE_ANGLE_FORE, HALF_CONE_ANGLE_AFT, SKIRT_ANGLE,
     >   SKIRT_LENGTH, X_CONE_JUNCTURE, R_CONE_JUNCTURE,
     >   CONE_LENGTH_FORE, CONE_LENGTH_AFT

      INTEGER ::
     >   NGEOM_SPOKE

      NAMELIST /SPHERE_CONE_INPUTS/
     >   X_NOSE, R_NOSE, RADIUS_NOSE, RADIUS_CONE_JUNCTURE, RADIUS_BASE,
     >   RADIUS_SHOULDER, RADIUS_VERTEX, HALF_CONE_ANGLE,
     >   HALF_CONE_ANGLE_FORE, HALF_CONE_ANGLE_AFT, SKIRT_ANGLE,
     >   SKIRT_LENGTH, NGEOM_SPOKE

*     Execution:
*     **********

      WRITE (LUNCRT, '(A)')

*     Optional namelist for controlling elliptic smoothing.
*     Note that orthogonality is forced at the internal stag. pt. bndry. too.

      CONV2D  = 10.0;   CONVMIN = 0.0;    DMAX2D   = 0.0
      ITMAX   = 5000;   ITFLOAT = 0;      ITFREEZE = 100;  JLAYER = 4
      FGLIMIT = 1.0;    URFG    = 0.1;    URFLOAT  = 0.1
      OMG2D   = 1.1;    POWERI  = 1.0;    POWERJ   = 0.5
      BGMODE  = 'YYYY'; FGMODE  = 'AAAN'; SPMODE   = 'YYYY'
      EXPI1   = 0.45;   EXPI2   = 0.45
      EXPJ1   = 0.45;   EXPJ2   = 0.45

*     Read the required control inputs (any file on standard input):

***   OPEN (LUNINP, FILE='hb_grid.inp', STATUS='OLD', IOSTAT=IOS)

      LINE = 1
      READ (LUNINP, *, IOSTAT=IOS) ! Skip title
      IF (IOS /= 0) GO TO 800
      LINE = LINE + 1
      READ (LUNINP, *, IOSTAT=IOS) FILENAME  ! Input geometry file or 'none'
      IF (IOS /= 0) GO TO 800

      ANALYTIC = FILENAME(1:4) == 'none' .OR. FILENAME(1:4) == 'NONE'

      IF (.NOT. ANALYTIC) THEN

         OPEN (LUNGEOM, FILE=FILENAME, STATUS='OLD', IOSTAT=IOS)
         IF (IOS /= 0) THEN
            WRITE (LUNCRT, '(/, 2A)') ' Trouble opening ', FILENAME
            GO TO 999
         END IF

*        Read the geometry (wrap-around format):

         READ (LUNGEOM, *) ! Title
         READ (LUNGEOM, *, IOSTAT=IOS) NGEOM
         IF (IOS /= 0) THEN
            WRITE (LUNCRT, '(/, A)') ' Trouble reading NGEOM.'
            GO TO 999
         END IF

         ALLOCATE (XGEOM(NGEOM), YGEOM(NGEOM))

         READ (LUNGEOM, *, IOSTAT=IOS)
     >      (XGEOM(I), YGEOM(I), I = 1, NGEOM)
         IF (IOS /= 0) THEN
            WRITE (LUNCRT, '(/, A)') ' Trouble reading geometry.'
            GO TO 999
         END IF

         CLOSE (LUNGEOM)

      END IF

*     An internal procedure below moves most of the reading clutter out:

      CALL READ_CONTROLS ()     ! ****** Bulk of the control inputs ******

      IF (IOS /= 0) GO TO 800

*     Any overrides for the elliptic smoothing controls?

      READ (LUNINP, NML=ELLIPTIC, IOSTAT=IOS)

      IF (IOS > 0) THEN
         WRITE (LUNCRT, '(/, 2A)')
     >      ' *** Trouble reading namelist ELLIPTIC; proceeding. ***'
      END IF

*     Construct a generatrix instead of reading the geometry?

      IF (ANALYTIC) THEN

         RADIUS_CONE_JUNCTURE = ZERO  ! Default means not a biconic
         RADIUS_VERTEX        = ZERO  ! Default means leave biconic vertex sharp

         READ (LUNINP, NML=SPHERE_CONE_INPUTS, IOSTAT=IOS)

         IF (IOS > 0) THEN
            WRITE (LUNCRT, '(/, 2A)')
     >         ' *** Trouble reading namelist SPHERE_CONE_INPUTS. ***'
            GO TO 999
         END IF

         SPHERECONE = RADIUS_CONE_JUNCTURE <= ZERO

         IF (SPHERECONE) THEN
            CAPSULE%X_NOSE          = X_NOSE
            CAPSULE%R_NOSE          = R_NOSE
            CAPSULE%RADIUS_NOSE     = RADIUS_NOSE
            CAPSULE%RADIUS_BASE     = RADIUS_BASE
            CAPSULE%RADIUS_SHOULDER = RADIUS_SHOULDER
            CAPSULE%HALF_CONE_ANGLE = HALF_CONE_ANGLE
            CAPSULE%SKIRT_ANGLE     = SKIRT_ANGLE
            CAPSULE%SKIRT_LENGTH    = SKIRT_LENGTH
         ELSE
            BICONIC%X_NOSE               = X_NOSE
            BICONIC%R_NOSE               = R_NOSE
            BICONIC%RADIUS_NOSE          = RADIUS_NOSE
            BICONIC%RADIUS_CONE_JUNCTURE = RADIUS_CONE_JUNCTURE
            BICONIC%RADIUS_BASE          = RADIUS_BASE
            BICONIC%RADIUS_SHOULDER      = RADIUS_SHOULDER
            BICONIC%RADIUS_VERTEX        = RADIUS_VERTEX
            BICONIC%HALF_CONE_ANGLE_FORE = HALF_CONE_ANGLE_FORE
            BICONIC%HALF_CONE_ANGLE_AFT  = HALF_CONE_ANGLE_AFT
            BICONIC%SKIRT_ANGLE          = SKIRT_ANGLE
            BICONIC%SKIRT_LENGTH         = SKIRT_LENGTH
         END IF

         NGEOM = 2*NGEOM_SPOKE - 1

         ALLOCATE (XGEOM(NGEOM), YGEOM(NGEOM))

         IF (SPHERECONE) THEN
            CALL SPHERE_CONE_CONTROL_POINTS (CAPSULE)
            CALL SPHERE_CONE_DISCRETIZATION (CAPSULE, NGEOM_SPOKE,
     >                                       XGEOM(NGEOM_SPOKE),
     >                                       YGEOM(NGEOM_SPOKE))
         ELSE
            CALL BICONIC_CONTROL_POINTS (BICONIC)
            CALL BICONIC_DISCRETIZATION (BICONIC, NGEOM_SPOKE,
     >                                   XGEOM(NGEOM_SPOKE),
     >                                   YGEOM(NGEOM_SPOKE))
            X_CONE_JUNCTURE  = BICONIC%X_CONE_JUNCTURE
            R_CONE_JUNCTURE  = BICONIC%R_CONE_JUNCTURE
            CONE_LENGTH_FORE = BICONIC%CONE_LENGTH_FORE
            CONE_LENGTH_AFT  = BICONIC%CONE_LENGTH_AFT
         END IF

         J = NGEOM_SPOKE
         DO I = NGEOM_SPOKE + 1, NGEOM
            J = J - 1
            XGEOM(J) =  XGEOM(I)
            YGEOM(J) = R_NOSE - (YGEOM(I) - R_NOSE)
         END DO

         OPEN  (LUNGEOM, FILE='generatrix.dat', status='unknown')
         WRITE (LUNGEOM, '(A)')
     >      '#                      X                       Y'
         WRITE (LUNGEOM, '(1P, 2E24.15)')
     >      (XGEOM(I), YGEOM(I), I = NGEOM_SPOKE, NGEOM)
         CLOSE (LUNGEOM)

      END IF

      CLOSE (LUNINP)

      CALL SANITY_CHECKS ()     ! ****** Another internal procedure ******

      IF (IOS /= 0) GO TO 999


*     ******************** Initial Euler-type gridding *******************

      NI = 2 * NI_EU - 1  ! Full wrap-around dimension

      ALLOCATE (X(NI,NJ_EU), Y(NI,NJ_EU))

*     Push the rest down an argument-driven level in case it can be reused
*     for a wing on a wall or perhaps a 3-D body.

      CALL HB_EULER (LUNCRT, LUNSHOCK, NGEOM, XGEOM, YGEOM, BODY_LENGTH,
     >               FMACH, ALPHA, GAMMA, RFUDGE, BFUDGE, AFUDGE, THETA,
     >               XC_DB, XC_SB, NI, NI_EU, NJ_EU, NBLAYER, RBLAYER,
     >               D1BODY, D2BODY, D1SHOCK, D2SHOCK,
     >               D1NOSE, D2NOSE, D1TAIL,  D2TAIL,
     >               SMOOTH, ITMAX, ITFLOAT, ITFREEZE, JLAYER,
     >               CONV2D, CONVMIN, DMAX2D, OMG2D,
     >               POWERI, POWERJ, EXPI1, EXPI2, EXPJ1, EXPJ2,
     >               FGLIMIT, URFG, URFLOAT, ITERATIONS,
     >               BGMODE, FGMODE, SPMODE,
     >               X, Y, IER)

      IF (FORMATTED) THEN
         OPEN (LUNGRID, FILE=FILENAME, STATUS='UNKNOWN')
      ELSE
         OPEN (LUNGRID, FILE=FILENAME, STATUS='UNKNOWN',
     >         FORM='UNFORMATTED')
      END IF

      IF (.NOT. REGRID) THEN

         IF (TWO_D) THEN
            ALLOCATE (Z(1,1))
         ELSE
            ALLOCATE (Z(NI, NJ_EU))
         END IF

         Z = ZERO

         CALL HB_SAVE (LUNGRID, FORMATTED, TWO_D, ALPHA, NI, NJ_EU,
     >                 X, Y, Z)

      ELSE ! Redistribute in the radial direction for N-S spacing:

         ALLOCATE (XNS(NI, NJ_NS), YNS(NI, NJ_NS))

         CALL HB_NS (LUNCRT, NI, NI_EU, NJ_EU, NJ_NS,
     >               NBLAYERNS, RBLAYERNS,
     >               BODY_LENGTH, D1NOSENS, D1TAILNS,
     >               D2NOSENS, D2TAILNS, ITERATIONS, X, Y, XNS, YNS)

         IF (TWO_D) THEN
            ALLOCATE (ZNS(1,1))
         ELSE
            ALLOCATE (ZNS(NI, NJ_NS))
         END IF

         ZNS = ZERO

         CALL HB_SAVE (LUNGRID, FORMATTED, TWO_D, ALPHA, NI, NJ_NS,
     >                 XNS, YNS, ZNS)

      END IF

      CLOSE (LUNGRID)

      GO TO 999


*     Error handling:
*     ***************

  800 WRITE (LUNCRT, '(/, A, I2)')
     >   ' Trouble reading control inputs.  Line number: ', LINE

  999 WRITE (LUNCRT, '(A)')


*     ************* Internal procedures for program HB_GRID **************

      CONTAINS

         SUBROUTINE READ_CONTROLS ()

!        Internal procedure for most of the control inputs, to remove clutter
!        from the main program.

         LINE = LINE + 1
         READ (LUNINP, *, IOSTAT=IOS) FILENAME ! Output grid file name
         IF (IOS /= 0) GO TO 800
         LINE = LINE + 1
         READ (LUNINP, *, IOSTAT=IOS) FORMATTED
         IF (IOS /= 0) GO TO 800
         LINE = LINE + 1
         READ (LUNINP, *, IOSTAT=IOS) TWO_D
         IF (IOS /= 0) GO TO 800
         LINE = LINE + 1
         READ (LUNINP, *, IOSTAT=IOS) FMACH
         IF (IOS /= 0) GO TO 800
         LINE = LINE + 1
         READ (LUNINP, *, IOSTAT=IOS) ALPHA
         IF (IOS /= 0) GO TO 800
         LINE = LINE + 1
         READ (LUNINP, *, IOSTAT=IOS) GAMMA
         IF (IOS /= 0) GO TO 800
         LINE = LINE + 1
         READ (LUNINP, *, IOSTAT=IOS) RFUDGE
         IF (IOS /= 0) GO TO 800
         LINE = LINE + 1
         READ (LUNINP, *, IOSTAT=IOS) BFUDGE
         IF (IOS /= 0) GO TO 800
         LINE = LINE + 1
         READ (LUNINP, *, IOSTAT=IOS) AFUDGE
         IF (IOS /= 0) GO TO 800
         LINE = LINE + 1
         READ (LUNINP, *, IOSTAT=IOS) THETA
         IF (IOS /= 0) GO TO 800
         LINE = LINE + 1
         READ (LUNINP, *, IOSTAT=IOS) XC_DB
         IF (IOS /= 0) GO TO 800
         LINE = LINE + 1
         READ (LUNINP, *, IOSTAT=IOS) XC_SB
         IF (IOS /= 0) GO TO 800
         LINE = LINE + 1
         READ (LUNINP, *, IOSTAT=IOS) NI_EU
         IF (IOS /= 0) GO TO 800
         LINE = LINE + 1
         READ (LUNINP, *, IOSTAT=IOS) NJ_EU
         IF (IOS /= 0) GO TO 800
         LINE = LINE + 1
         READ (LUNINP, *, IOSTAT=IOS) D1BODY
         IF (IOS /= 0) GO TO 800
         LINE = LINE + 1
         READ (LUNINP, *, IOSTAT=IOS) D2BODY
         IF (IOS /= 0) GO TO 800
         LINE = LINE + 1
         READ (LUNINP, *, IOSTAT=IOS) D1SHOCK
         IF (IOS /= 0) GO TO 800
         LINE = LINE + 1
         READ (LUNINP, *, IOSTAT=IOS) D2SHOCK
         IF (IOS /= 0) GO TO 800
         LINE = LINE + 1
         READ (LUNINP, *, IOSTAT=IOS) NBLAYER
         IF (IOS /= 0) GO TO 800
         LINE = LINE + 1
         READ (LUNINP, *, IOSTAT=IOS) RBLAYER
         IF (IOS /= 0) GO TO 800
         LINE = LINE + 1
         READ (LUNINP, *, IOSTAT=IOS) D1NOSE
         IF (IOS /= 0) GO TO 800
         LINE = LINE + 1
         READ (LUNINP, *, IOSTAT=IOS) D1TAIL
         IF (IOS /= 0) GO TO 800
         LINE = LINE + 1
         READ (LUNINP, *, IOSTAT=IOS) D2NOSE
         IF (IOS /= 0) GO TO 800
         LINE = LINE + 1
         READ (LUNINP, *, IOSTAT=IOS) D2TAIL
         IF (IOS /= 0) GO TO 800
         LINE = LINE + 1
         READ (LUNINP, *, IOSTAT=IOS) SMOOTH
         IF (IOS /= 0) GO TO 800
         LINE = LINE + 1
         READ (LUNINP, *, IOSTAT=IOS) ITERATIONS
         IF (IOS /= 0) GO TO 800
         LINE = LINE + 1
         READ (LUNINP, *, IOSTAT=IOS) REGRID
         IF (IOS /= 0) GO TO 800
         LINE = LINE + 1
         READ (LUNINP, *, IOSTAT=IOS) NJ_NS
         IF (IOS /= 0) GO TO 800
         LINE = LINE + 1
         READ (LUNINP, *, IOSTAT=IOS) NBLAYERNS
         IF (IOS /= 0) GO TO 800
         LINE = LINE + 1
         READ (LUNINP, *, IOSTAT=IOS) RBLAYERNS
         IF (IOS /= 0) GO TO 800
         LINE = LINE + 1
         READ (LUNINP, *, IOSTAT=IOS) D1NOSENS
         IF (IOS /= 0) GO TO 800
         LINE = LINE + 1
         READ (LUNINP, *, IOSTAT=IOS) D1TAILNS
         IF (IOS /= 0) GO TO 800
         LINE = LINE + 1
         READ (LUNINP, *, IOSTAT=IOS) D2NOSENS
         IF (IOS /= 0) GO TO 800
         LINE = LINE + 1
         READ (LUNINP, *, IOSTAT=IOS) D2TAILNS
         IF (IOS /= 0) GO TO 800

  800    RETURN

         END SUBROUTINE READ_CONTROLS

         SUBROUTINE SANITY_CHECKS ()

!        Internal procedure for HB_GRID to do some control input checking.

         ! <To be completed>

         END SUBROUTINE SANITY_CHECKS

      END PROGRAM HB_GRID

C*******************************************************************************
*
      SUBROUTINE HB_EULER (LUNCRT, LUNSHOCK, NGEOM, XGEOM, YGEOM,
     >                     BODY_LENGTH, FMACH, ALPHA, GAMMA,
     >                     RFUDGE, BFUDGE, AFUDGE,
     >                     THETA, XC_DB, XC_SB, NI, NI_EU, NJ_EU,
     >                     NBLAYER, RBLAYER,
     >                     D1BODY, D2BODY, D1SHOCK, D2SHOCK,
     >                     D1NOSE, D2NOSE, D1TAIL,  D2TAIL,
     >                     SMOOTH, ITMAX, ITFLOAT, ITFREEZE, JLAYER,
     >                     CONV2D, CONVMIN, DMAX2D, OMG2D,
     >                     POWERI, POWERJ, EXPI1, EXPI2, EXPJ1, EXPJ2,
     >                     FGLIMIT, URFG, URFLOAT, ITERATIONS,
     >                     BGMODE, FGMODE, SPMODE,
     >                     X, Y, IER)
*
*        Generate a 2-D structured grid for a hypersonic body defined by a
*     discretized curve, not necessarily symmetric.  The outer boundary is
*     a function of the intended Mach number, angle of attack, and surface
*     radius at the stagnation point where the free stream is normal to the
*     surface.  The downstream boundary is placed at a specified fraction of
*     the body length from the nose, as this may vary with Mach.
*
*        The single-block topology is constructed in lower and upper halves
*     via transfinite interpolation and elliptic smoothing with Euler-type
*     spacing.  A separate routine can redistribute to N-S-type spacing.
*
********************************************************************************

      USE TRIGD
      IMPLICIT NONE

*     Arguments:

      INTEGER, INTENT (IN) ::
     >   LUNCRT,            ! Screen
     >   LUNSHOCK,          ! Estimated shock curve for plotting
     >   NGEOM              ! Length of wrap-around body curve

      REAL, INTENT (IN) ::
     >   XGEOM(NGEOM),      ! Body coordinates, lower-to-upper
     >   YGEOM(NGEOM)

      REAL, INTENT (OUT) ::
     >   BODY_LENGTH        ! Full body length found, for use by HB_NS

      REAL, INTENT (IN) ::
     >   FMACH, ALPHA,      ! Intended Mach & angle of attack (degrees)
     >   GAMMA,             ! Ratio of specific heats
     >   RFUDGE,            ! For fudging the nose radius (e.g., 1.5)
     >   BFUDGE,            ! For fudging limiting shock angle, beta (~2.?)
     >   AFUDGE,            ! For reducing the Alpha-related rotation
     >   THETA,             ! Effective body deflection angle, degrees
     >   XC_DB,             ! % body length at which to place downstream bndry.
     >   XC_SB              ! "   "   "   " from stag. pt. (sonic bubble edges)

      INTEGER, INTENT (IN) ::
     >   NI,                ! Desired full wrap-around grid dimension
     >   NI_EU, NJ_EU,      ! Desired half-grid dimensions (Euler-type spacing)
     >   NBLAYER            ! See main program or BLGRID

      REAL, INTENT (IN) ::
     >   RBLAYER,           !  "    "
     >   D1BODY, D2BODY,    ! Body surface grid spacings (fraction of lower arc
                            ! length) at stag. point and downstream boundary;
                            ! D1BODY = 0. => curvature-based distribution, else
                            ! D2BODY > 0. => 2-sided Vinokur, else
                            ! D2BODY < 0. => one-sided stretching from stag. pt.
     >   D1SHOCK, D2SHOCK,  ! Outer boundary spacings (fraction of lower bound-
                            ! ary arc length) as for D1BODY and D2BODY
     >   D1NOSE, D2NOSE,    ! Radial increments at surface (linear in between)
     >   D1TAIL, D2TAIL     ! ... and at the upstream boundary, as fractions of
                            ! the radial grid line lengths; d2* < 0 are not
                            ! used, in favor of one-sided stretching
      LOGICAL, INTENT (IN) ::
     >   SMOOTH             ! F suppresses elliptic smoothing

      INTEGER, INTENT (IN) ::             ! See ELLIP2D descriptions
     >   ITMAX, ITFLOAT, ITFREEZE, JLAYER

      REAL, INTENT (IN) ::                ! See ELLIP2D descriptions
     >   CONV2D, CONVMIN, DMAX2D, OMG2D,
     >   POWERI, POWERJ, EXPI1, EXPI2, EXPJ1, EXPJ2,
     >   FGLIMIT, URFG, URFLOAT

      LOGICAL, INTENT (IN) ::
     >   ITERATIONS         ! T turns on elliptic smoothing iteration printing

      CHARACTER * 4, INTENT (IN) ::       ! See ELLIP2D descriptions
     >   BGMODE, FGMODE, SPMODE

      REAL, INTENT (OUT), DIMENSION (NI, NJ_EU) ::
     >   X, Y               ! Desired grid (single-block)

      INTEGER, INTENT (OUT) ::
     >   IER                ! 0 means no error

*     Local constants:

      REAL, PARAMETER ::
     >   EIGHTH = 0.125, HALF = 0.5, ONE = 1.0, TWO = 2.0, ZERO = 0.0

      LOGICAL, PARAMETER ::
     >   FALSE = .FALSE., NEW = .TRUE., TRUE = .TRUE.

      CHARACTER, PARAMETER ::
     >   LINEAR = 'L', LOOSE * 1 = 'B', TIGHT='M'

*     Automatic arrays:

      REAL
     >   ARCS(NI,NJ_EU,2), DERIVSI(NI_EU), DERIVSJ(NJ_EU),
     >   TEVAL(NI),       TGEOM(NGEOM),    TGRID(NI),    TOUTER(NI),
     >   TSHOCK(NI),      XSHOCK(NI),      YSHOCK(NI),
     >   TOUTFLOW(NJ_EU), XOUTFLOW(NJ_EU), YOUTFLOW(NJ_EU),
     >   TRADIAL(NJ_EU),  XRADIAL(NJ_EU),  YRADIAL(NJ_EU)

*     Local variables:

      INTEGER
     >   I, I1, I2, ILE, ISH, ISL, ISU, IT, ITL, ITU, J, LUNITER,
     >   NGNOSE, NGNOSE2

      REAL
     >   ANGLE, B, BETA, BLENGTH, D1SHOCKL, D1SHOCKU, D1SONIC, D2SONIC,
     >   D2XDY2, DELTA, DS, DT, DX, DXDT, DXDY, DY, DYDT,
     >   FMSQ, PI, R, RC, RLENGTH, SLENGTH, SLOPE,
     >   T1, T2, TANBETA, TLENGTH, TOL, TOPT, TSONIC, TSTAG, TTOTAL,
     >   X_SB, XINT, XMAX, XMIN, XTE, XVERTEX, YINT, YMAX, YVERTEX,
     >   TLINE(4), XLINE(4), YLINE(4)

      CHARACTER * 4
     >   BGSWAP, FGSWAP

*     Execution:
*     **********

      PI  = 2.0 * ASIN (ONE)
      TOL = MAX (10.0 * EPSILON (TOL), 1.E-14)

      IF (ITERATIONS) THEN ! Both Newton-type and elliptic smoothing
         LUNITER =  LUNCRT
      ELSE
         LUNITER = -LUNCRT
      END IF

*     Determine the leading edge, body length, and trailing edge of grid:

      XMIN = XGEOM(1);  ILE = 1

      DO I = 2, NGEOM
         IF (XGEOM(I) < XMIN) THEN
            XMIN = XGEOM(I)
            ILE  = I
         END IF
      END DO

      BODY_LENGTH = HALF * (XGEOM(1) + XGEOM(NGEOM)) - XMIN  ! Full, for HB_NS
      X_SB    = BODY_LENGTH * XC_SB * 0.01     ! From stag.pt. to edge of s.b.
      BLENGTH = BODY_LENGTH * XC_DB * 0.01     ! Active length of body
      XTE     = XGEOM(ILE) + BLENGTH

      ITL = ILE - 1
      CALL INTERVAL (ILE, XGEOM, XTE, -ONE, ITL)
      IF (XTE - XGEOM(ITL + 1) < XGEOM(ITL) - XTE) ITL = ITL + 1

      ITU = 2
      CALL INTERVAL (NGEOM - ILE + 1, XGEOM(ILE), XTE, ONE, ITU)
      ITU = ITU + ILE - 1
      IF (XGEOM(ITU+1) - XTE < XTE - XGEOM(ITU)) ITU = ITU + 1

      NGNOSE = ITU - ITL + 1  ! Geometry index range to be gridded

*     Arc lengths for this part of the body (originally; all of it now):

**    CALL CHORDS2D (NGNOSE, XGEOM(ITL), YGEOM(ITL), FALSE, TTOTAL,
**   >               TGEOM(ITL))
      CALL CHORDS2D (NGEOM, XGEOM, YGEOM, FALSE, TTOTAL, TGEOM)

      TTOTAL = TGEOM(ITU) - TGEOM(ITL)


*     Locate the stagnation point for this Alpha,
*     where the appropriate tangent touches the body:
*     ***********************************************

      IF (ALPHA == ZERO) THEN
         TSTAG = TGEOM(ILE)
         IT = ILE
      ELSE
         SLOPE = -180.0 / (ALPHA * PI)
         I1 = 1
         I2 = ILE - ITL + 1 ! I.e., the leading edge point
         IF (ALPHA < ZERO) THEN
            I1 = I2
            I2 = NGNOSE
         END IF

         CALL PLSTANGENT (NGNOSE, XGEOM(ITL), YGEOM(ITL), TGEOM(ITL),
     >                    I1, I2, FALSE, SLOPE, TOL, TSTAG, IT,
     >                    LUNITER, IER)
         IF (IER /= 0) THEN
            WRITE (LUNCRT, '(/, A, I3)')
     >         ' Trouble locating stag. point. PLSTANGENT code:', IER
            IT = ILE - 1 ! No need to offset it - not used elsewhere
            TSTAG = TGEOM(IT)
         END IF

      END IF

*     Computational point distribution along the body, lower surface first:
*     *********************************************************************

      TEVAL(1)     = TGEOM(ITL)
      TEVAL(NI_EU) = TSTAG
      TLENGTH      = TSTAG - TEVAL(1)
      NGNOSE2      = IT - ITL + 1

      IF (D1BODY == ZERO) THEN  ! Curvature-based surface grid (arcs only)

         CALL CURVDIS (NGNOSE2, XGEOM(ITL), YGEOM(ITL), NI_EU, HALF,
     >                 3, LUNITER, TRUE, TGRID, TGRID, IER)

         IF (IER /= 0) THEN
            WRITE (LUNCRT, '(/, A, I3)')
     >         ' CURVDIS trouble with lower surface grid.  IER:', IER
            GO TO 999
         END IF

         TGRID(1:NI_EU) = TGRID(1:NI_EU) + TGEOM(ITL)

      ELSE IF (D2BODY > ZERO) THEN  ! 2-sided Vinokur, no geometric portion

         CALL BLGRID (NI_EU, D1BODY*TLENGTH, D2BODY*TLENGTH,
     >                2, ONE, TEVAL, LUNCRT, IER)

         IF (IER /= 0) THEN
            WRITE (LUNCRT, '(/, A, I3)')
     >         ' BLGRID trouble with lower surface grid.  IER:', IER
            GO TO 999
         END IF

      ELSE ! One-sided stretching for the body:

         CALL EXPDIS5 (1, TGEOM(ITL), TSTAG, D1BODY*TLENGTH, NI_EU,
     >                 TEVAL, LUNITER)
      END IF

*     Impose the computational grid on the lower body:

      IF (D1BODY /= ZERO) THEN
         DO I = 1, NI_EU ! Reverse the clustering
            TGRID(I) = TSTAG - TEVAL(NI_EU + 1 - I)
         END DO
      END IF

      D1SHOCKL = TGRID(NI_EU) - TGRID(NI_EU-1)  ! For possible use at shock

      CALL LCSFIT (NGNOSE, TGEOM(ITL), XGEOM(ITL), NEW, TIGHT,
     >             NI_EU, TGRID, X(1,1), DERIVSI)

      CALL LCSFIT (NGNOSE, TGEOM(ITL), YGEOM(ITL), NEW, LOOSE,
     >             NI_EU, TGRID, Y(1,1), DERIVSI)

*     Upper body surface:

      TGRID(NI_EU) = TSTAG
      TGRID(NI)    = TGEOM(ITU)
      TLENGTH      = TGEOM(ITU) - TSTAG
      NGNOSE2      = ITU - IT + 1

      IF (D1BODY == ZERO) THEN  ! Curvature-based surface grid (arcs only)

         CALL CURVDIS (NGNOSE2, XGEOM(IT), YGEOM(IT), NI_EU, HALF,
     >                3, LUNITER, TRUE, TGRID(NI_EU), TGRID(NI_EU), IER)

         IF (IER /= 0) THEN
            WRITE (LUNCRT, '(/, A, I3)')
     >         ' CURVDIS trouble with upper surface grid.  IER:', IER
            GO TO 999
         END IF

         TGRID(NI_EU:NI) = TGRID(NI_EU:NI) + TGEOM(IT)

      ELSE IF (D2BODY > ZERO) THEN  ! 2-sided Vinokur

         CALL BLGRID (NI_EU, D1BODY*TLENGTH, D2BODY*TLENGTH,
     >                2, ONE, TGRID(NI_EU), LUNCRT, IER)

         IF (IER /= 0) THEN
            WRITE (LUNCRT, '(/, A, I3)')
     >         ' BLGRID trouble with upper surface grid.  IER:', IER
            GO TO 999
         END IF

      ELSE ! One-sided stretching for the body:

         CALL EXPDIS5 (1, TSTAG, TGEOM(ITU), D1BODY*TLENGTH, NI_EU,
     >                 TGRID(NI_EU), LUNITER)
      END IF

      D1SHOCKU = TGRID(NI_EU+1) - TGRID(NI_EU)  ! For possible use at shock

      CALL LCSFIT (NGNOSE, TGEOM(ITL), XGEOM(ITL), NEW, TIGHT,
     >             NI_EU, TGRID(NI_EU), X(NI_EU,1), DERIVSI)

      CALL LCSFIT (NGNOSE, TGEOM(ITL), YGEOM(ITL), NEW, LOOSE,
     >             NI_EU, TGRID(NI_EU), Y(NI_EU,1), DERIVSI)


*     Radius of curvature at the stagnation point (by switching X and Y):
*     *******************************************************************

      CALL FDCNTR (NI_EU, Y(1,1), X(1,1), DXDY, D2XDY2)
      CALL FDCURV (DXDY, D2XDY2, R)  ! R = curvature here; could be ~zero

      IF (ABS (R) < BODY_LENGTH * 0.01) R = BODY_LENGTH * 0.01
      R = RFUDGE / R

*     Shock stand-off distance with a pad of 10% (p. 189, Eq. 5.37):

      FMSQ  = FMACH * FMACH
***   DELTA = 0.386 * R * EXP (4.67 / FMSQ)  ! Cylinder-wedge
      DELTA = 0.143 * R * EXP (3.24 / FMSQ)  ! Sphere-cone

*     Shock radius of curvature at vertex (p. 189, Eq. 5.38):

***   RC = 1.386 * R * EXP (1.80 / (FMACH - ONE) ** 0.75)  ! Cylinder-wedge
      RC = 1.143 * R * EXP (0.54 / (FMACH - ONE) ** 1.20)  ! Sphere-cone

*     Estimate the limiting shock angle BETA (p. 36, Eq. 2.16):

      CALL SHOCK_ANGLE ()  ! Local procedure below

      BETA = BFUDGE * BETA;  TANBETA = TAN (BETA)

*     Shock Y downstream (zero Alpha for now; shock Y = 0 at X = R + delta):

      XMAX = XGEOM(ITL) + R + TWO*DELTA  ! Add shift applied below
      IF (ALPHA /= ZERO) XMAX = XMAX + BLENGTH  ! Allow for rotation
      YMAX = (RC / TANBETA) *
     >   SQRT ((ONE + TANBETA**2 * (XMAX - R - DELTA) / RC)**2 - ONE)

      WRITE (LUNCRT, '(1x, a, i15)') 'ITL:                       ', ITL
      WRITE (LUNCRT, '(1x, a, f15.10)')
     >    'BLENGTH:                   ', BLENGTH,
     >    'XGEOM(ITL):                ', XGEOM(ITL),
     >    'XMAX:                      ', XMAX,
     >    'YMAX:                      ', YMAX,
     >    'Shock stand-off, DELTA:    ', DELTA,
     >    'Limiting angle, BETA, deg: ', BETA * (180./PI),
     >    'Stag. pt. body radius, R:  ', R,
     >    '... which includes RFUDGE: ', RFUDGE,
     >    'Corresp. shock radius, RC: ', RC

*     Generate the shock curve shape for interpolating later (p. 189, Eq. 5.36):
*     **************************************************************************

      DY = YMAX / REAL (NI_EU - 1);  J = NI_EU

      DO I = NI_EU, NI
         YSHOCK(I) = DY * REAL (I - NI_EU)
         YSHOCK(J) = -YSHOCK(I)
         XSHOCK(I) = (R + DELTA) + (RC / TANBETA**2) *
     >               (SQRT (ONE  + (YSHOCK(I) * TANBETA / RC)**2) - ONE)
         XSHOCK(J) = XSHOCK(I)
         J = J - 1
      END DO

*     Translate and rotate the shock boundary to match Alpha [*AFUDGE]:

      XVERTEX = X(NI_EU,1) - DELTA * COSD (ALPHA)
      YVERTEX = Y(NI_EU,1) - DELTA * SIND (ALPHA)

      DX = XSHOCK(NI_EU) - XVERTEX;  XSHOCK = XSHOCK - DX
      DY = YSHOCK(NI_EU) - YVERTEX;  YSHOCK = YSHOCK - DY

      IF (ALPHA /= ZERO) THEN

         CALL ROTATE2D (NI, XSHOCK, YSHOCK, ALPHA*AFUDGE, XVERTEX,
     >                  YVERTEX)

      ELSE  ! Enforce exactness forward of the stag. pt.:

         YSHOCK(NI_EU) = Y(NI_EU,1)  ! Normally zero

      END IF

      CALL CHORDS2D (NI, XSHOCK, YSHOCK, FALSE, TTOTAL, TSHOCK)

      OPEN (LUNSHOCK, FILE='shock-x+y+t.dat', status='unknown')
      WRITE (LUNSHOCK, '(A)')
     >   '#        XSHOCK         YSHOCK         TSHOCK'
      WRITE (LUNSHOCK, '(1P, 3E15.6)')
     >   (XSHOCK(I), YSHOCK(I), TSHOCK(I), I = 1, NI)
      CLOSE (LUNSHOCK)

*     Interior boundary at the stagnation point (straight line):
*     **********************************************************

*     It is crucial to get this exact for the axisymmetric case.

      XLINE(1) = X(NI_EU,1);  XLINE(2) = XSHOCK(NI_EU)
      YLINE(1) = Y(NI_EU,1);  YLINE(2) = YSHOCK(NI_EU)

      CALL CHORDS2D (2, XLINE, YLINE, FALSE, RLENGTH, TLINE)

      TRADIAL(1)     = ZERO
      TRADIAL(NJ_EU) = RLENGTH

      IF (D2NOSE > ZERO) THEN ! Geometric + Vinokur

         CALL BLGRID (NJ_EU, D1NOSE*RLENGTH, D2NOSE*RLENGTH,
     >                NBLAYER, RBLAYER, TRADIAL, LUNCRT, IER)

         IF (IER /= 0) THEN
            WRITE (LUNCRT, '(/, A, I3)')
     >         ' BLGRID trouble at stagnation point.  IER:', IER
            GO TO 999
         END IF

      ELSE ! One-sided stretching

         CALL EXPDIS5 (1, ZERO, RLENGTH, D1NOSE*RLENGTH, NJ_EU,
     >                 TRADIAL, LUNITER)
      END IF

      CALL LCSFIT (2, TLINE, XLINE, NEW, LINEAR, NJ_EU, TRADIAL,
     >             XRADIAL, DERIVSJ)
      CALL LCSFIT (2, TLINE, YLINE, NEW, LINEAR, NJ_EU, TRADIAL,
     >             YRADIAL, DERIVSJ)

      IF (ALPHA == ZERO) YRADIAL(1:NJ_EU) = Y(NI_EU,1)

      DO J = 1, NJ_EU
         X(NI_EU,J) = XRADIAL(J)
         Y(NI_EU,J) = YRADIAL(J)
      END DO


*     Downstream outflow boundary, normal to the end-of-body and shock bndry.:
*     ************************************************************************

*     Lower outflow edge first.
**    Original scheme:  find where a straight-line outflow boundary intersects
**    the shock boundary.  Set up a 2-point line:

**    ANGLE = ATAN ((YGEOM(ITL) - YGEOM(ITL+1)) /
**   >              (XGEOM(ITL) - XGEOM(ITL+1))) ! Slope of lower TE

**    XLINE(1) = XGEOM(ITL)
**    YLINE(1) = YGEOM(ITL)
**    DY       = YLINE(1) - YSHOCK(1)
**    XLINE(2) = XLINE(1) + DY * SIN (ANGLE)
**    YLINE(2) = YLINE(1) - DY * COS (ANGLE)

**    I1 = 2;  I2 = 1  ! Starting guesses

**    CALL INTSEC2T (2, XLINE, YLINE, TLINE, I1, TRUE,
**   >               NI_EU, XSHOCK, YSHOCK, TSHOCK, I2, FALSE, TOL,
**   >               XINT, YINT, T1, T2, LUNITER, IER)

**    IF (IER /= 0) THEN
**       WRITE (LUNCRT, '(/, A, I3)')
**   >      ' INTSEC2T trouble finding lower corner.  IER: ', IER
**       GO TO 999
**    END IF

**    RLENGTH = SQRT ((XINT - XGEOM(ITL))**2 + (YINT - YGEOM(ITL))**2)
**    TRADIAL(NJ_EU) = RLENGTH

*     Improved outflow edge:  slide its outer end along the shock curve so as
*     to minimize its peak curvature while forcing orthogonality at both ends
*     (a 1-D optimization problem):

*     The lower outer boundary point must be below the end of the geometry:

      DO I = 1, NI_EU
         IF (YSHOCK(I) > YGEOM(ITL)) THEN
            ISH = I - 1
            EXIT
         END IF
      END DO

      IER = 1  ! 0 suppresses iteration printout

      CALL MIN_MAX_CURV_EDGE_2D
     >   (NGEOM, XGEOM, YGEOM, TGEOM, ITL, TRUE,
     >    NI, XSHOCK, YSHOCK, TSHOCK, 1, ISH, FALSE, TOPT,
     >    NJ_EU, XOUTFLOW, YOUTFLOW, TOUTFLOW, IER)

      IF (IER /= 0) THEN
         WRITE (LUNCRT, '(/, A, I3)')
     >    ' MIN_MAX_CURV_EDGE_2D trouble, lower outflow boundary; IER:',
     >    IER
         GO TO 999
      END IF

      T2 = TOPT  ! To match straight-line scheme

*     Update the outflow boundary arc lengths (currently from a 6-point form):

      CALL CHORDS2D (NJ_EU, XOUTFLOW, YOUTFLOW, FALSE, RLENGTH,
     >               TOUTFLOW)

      IF (D2TAIL > ZERO) THEN ! Two-sided stretching

         CALL BLGRID (NJ_EU, D1TAIL*RLENGTH, D2TAIL*RLENGTH,
     >                NBLAYER, RBLAYER, TRADIAL, LUNCRT, IER)

         IF (IER /= 0) THEN
            WRITE (LUNCRT, '(/, A, I3)')
     >         ' BLGRID trouble at lower downstream boundary.  IER: ',
     >         IER
            GO TO 999
         END IF

      ELSE ! One-sided

         CALL EXPDIS5 (1, ZERO, RLENGTH, D1TAIL*RLENGTH, NJ_EU,
     >                 TRADIAL, LUNITER)
      END IF

**    CALL LCSFIT (2, TLINE, XLINE, NEW, LINEAR, NJ_EU, TRADIAL,
**   >             XRADIAL, DERIVSJ)
**    CALL LCSFIT (2, TLINE, YLINE, NEW, LINEAR, NJ_EU, TRADIAL,
**   >             YRADIAL, DERIVSJ)

      CALL LCSFIT (NJ_EU, TOUTFLOW, XOUTFLOW, NEW, LOOSE, NJ_EU,
     >             TRADIAL, XRADIAL, DERIVSJ)
      CALL LCSFIT (NJ_EU, TOUTFLOW, YOUTFLOW, NEW, LOOSE, NJ_EU,
     >             TRADIAL, YRADIAL, DERIVSJ)

      DO J = 1, NJ_EU
         X(1,J) = XRADIAL(J)
         Y(1,J) = YRADIAL(J)
      END DO

*     Lower computational outer (shock) boundary:

      TEVAL(1)     = T2
      TEVAL(NI_EU) = TSHOCK(NI_EU)
      SLENGTH      = TSHOCK(NI_EU) - T2

      IF (D1SHOCK < ZERO) THEN ! Use the body spacing at the stag. pt.
         DT = D1SHOCKL
      ELSE
         DT = D1SHOCK*SLENGTH
      END IF

      IF (D2SHOCK > ZERO) THEN ! 2-sided Vinokur (no geometric growth portion)

         CALL BLGRID (NI_EU, DT, D2SHOCK*SLENGTH, 2, ONE, TEVAL, LUNCRT,
     >                IER)

         IF (IER /= 0) THEN
            WRITE (LUNCRT, '(/, A, I3, /, A, 1P, 2E16.8)')
     >         ' BLGRID trouble with lower outer boundary.  IER:', IER,
     >         ' D1, D2:', DT, D2SHOCK*SLENGTH
            GO TO 999
         END IF

      ELSE ! One-sided stretching for the body:

         CALL EXPDIS5 (1, T2, TSHOCK(NI_EU), DT, NI_EU, TEVAL, LUNITER)

      END IF

!!!   write (6, '(/, a)') 'i, teval'
!!!   write (6, '(i4, 1p, e15.6)') (i, teval(i), i = 1, ni_eu)

*     Reverse the clustering for the lower shock boundary, and interpolate:

      DO I = 1, NI_EU
         TGRID(I) = T2 + (TEVAL(NI_EU) - TEVAL(NI_EU + 1 - I))
      END DO
      TGRID(NI_EU) = TSHOCK(NI_EU)

      CALL LCSFIT (NI, TSHOCK, XSHOCK, NEW, TIGHT,
     >             NI_EU, TGRID, X(1,NJ_EU), DERIVSI)

      CALL LCSFIT (NI, TSHOCK, YSHOCK, NEW, LOOSE,
     >             NI_EU, TGRID, Y(1,NJ_EU), DERIVSI)

!!!   write (6, '(/, a)') 'i, tgrid, xgrid, ygrid'
!!!   write (6, '(i4, 1p, 3e15.6)')
!!!  >   (i, tgrid(i), x(i,nj_eu), y(i,nj_eu), i = 1, ni_eu)

*     Repeat for upper boundaries:
*     ****************************

**    Original scheme for upper outflow boundary (straight line):

**    ANGLE = ATAN ((YGEOM(ITU) - YGEOM(ITU-1)) /
**   >              (XGEOM(ITU) - XGEOM(ITU-1))) ! Slope of upper TE

**    XLINE(1) = XGEOM(ITU)
**    YLINE(1) = YGEOM(ITU)
**    DY = YSHOCK(NI) - YLINE(1)
**    XLINE(2) = XLINE(1) - DY * SIN (ANGLE)
**    YLINE(2) = YLINE(1) + DY * COS (ANGLE)

**    I1 = 2;  I2 = NI

**    CALL INTSEC2T (2, XLINE, YLINE, TLINE, I1, TRUE,
**   >               NI, XSHOCK, YSHOCK, TSHOCK, I2, FALSE, TOL,
**   >               XINT, YINT, T1, T2, LUNITER, IER)

**    IF (IER /= 0) THEN
**       WRITE (LUNCRT, '(/, A, I3)')
**   >      ' INTSEC2T trouble finding upper corner.  IER:', IER
**       GO TO 999
**    END IF

**    RLENGTH = SQRT ((XINT - XGEOM(ITU))**2 + (YINT - YGEOM(ITU))**2)
**    TRADIAL(NJ_EU) = RLENGTH

*     Improved outflow edge:  slide its outer end along the shock curve so as
*     to minimize its peak curvature while forcing orthogonality at both ends.

*     The upper outer boundary point must be above the end of the geometry:

      DO I = NI_EU, NI
         IF (YSHOCK(I) > YGEOM(ITU)) THEN
            ISH = I
            EXIT
         END IF
      END DO

      IER = 1  ! 0 suppresses iteration printout

      CALL MIN_MAX_CURV_EDGE_2D
     >   (NGEOM, XGEOM, YGEOM, TGEOM, ITU, TRUE,
     >    NI, XSHOCK, YSHOCK, TSHOCK, ISH, NI, FALSE, TOPT,
     >    NJ_EU, XOUTFLOW, YOUTFLOW, TOUTFLOW, IER)

      IF (IER /= 0) THEN
         WRITE (LUNCRT, '(/, A, I3)')
     >    ' MIN_MAX_CURV_EDGE_2D trouble, upper outflow boundary; IER:',
     >    IER
         GO TO 999
      END IF

      T2 = TOPT  ! To match straight-line scheme

      CALL INTERVAL (NI, TSHOCK, T2, ONE, I2)

*     Update the outflow boundary arc lengths (currently from a 6-point form):

      CALL CHORDS2D (NJ_EU, XOUTFLOW, YOUTFLOW, FALSE, RLENGTH,
     >               TOUTFLOW)

      IF (D2TAIL > ZERO) THEN

         CALL BLGRID (NJ_EU, D1TAIL*RLENGTH, D2TAIL*RLENGTH,
     >                NBLAYER, RBLAYER, TRADIAL, LUNCRT, IER)

         IF (IER /= 0) THEN
            WRITE (LUNCRT, '(/, A, I3)')
     >         ' BLGRID trouble at upper downstream boundary.  IER: ',
     >         IER
            GO TO 999
         END IF

      ELSE ! One-sided

         CALL EXPDIS5 (1, ZERO, RLENGTH, D1TAIL*RLENGTH, NJ_EU,
     >                 TRADIAL, LUNITER)
      END IF

**    CALL LCSFIT (2, TLINE, XLINE, NEW, LINEAR, NJ_EU, TRADIAL,
**   >             XRADIAL, DERIVSJ)
**    CALL LCSFIT (2, TLINE, YLINE, NEW, LINEAR, NJ_EU, TRADIAL,
**   >             YRADIAL, DERIVSJ)

      CALL LCSFIT (NJ_EU, TOUTFLOW, XOUTFLOW, NEW, LOOSE, NJ_EU,
     >             TRADIAL, XRADIAL, DERIVSJ)
      CALL LCSFIT (NJ_EU, TOUTFLOW, YOUTFLOW, NEW, LOOSE, NJ_EU,
     >             TRADIAL, YRADIAL, DERIVSJ)

      DO J = 1, NJ_EU
         X(NI,J) = XRADIAL(J)
         Y(NI,J) = YRADIAL(J)
      END DO

*     Upper computational shock boundary:

      TGRID(NI) = T2

      IF (D1SHOCK < ZERO) THEN ! Use the body spacing at the stag. pt.
         DT = D1SHOCKU
      ELSE
         DT = D1SHOCK*SLENGTH
      END IF

      IF (D2SHOCK > ZERO) THEN ! 2-sided Vinokur, no geometric growth portion

         CALL BLGRID (NI_EU, DT, D2SHOCK*(T2 - TSHOCK(NI_EU)), 2, ONE,
     >                TGRID(NI_EU), LUNCRT, IER)

         IF (IER /= 0) THEN
            WRITE (LUNCRT, '(/, A, I3)')
     >         ' BLGRID trouble with upper outer boundary.  IER:', IER
            GO TO 999
         END IF

      ELSE ! One-sided stretching for the shock:

         CALL EXPDIS5 (1, TSHOCK(NI_EU), T2, DT, NI_EU, TGRID(NI_EU),
     >                 LUNITER)
      END IF

*     Interpolate:

      CALL LCSFIT (NI, TSHOCK, XSHOCK, NEW, TIGHT,
     >             NI_EU, TGRID(NI_EU), X(NI_EU,NJ_EU), DERIVSI)

      CALL LCSFIT (NI, TSHOCK, YSHOCK, NEW, LOOSE,
     >             NI_EU, TGRID(NI_EU), Y(NI_EU,NJ_EU), DERIVSI)

!!!   write (6, '(/, a)') 'i, tgrid, xgrid, ygrid'
!!!   write (6, '(i4, 1p, 3e15.6)')
!!!  >   (i, tgrid(i), x(i,nj_eu), y(i,nj_eu), i = ni_eu, ni)

*     Interior boundaries near the edges of the sonic bubble:
*     *******************************************************

*     Grid point (not geometry-based) arc-lengths in the right order:

      CALL CHORDS2D (NI, X(1,1), Y(1,1), FALSE, TTOTAL, TGRID)
      CALL CHORDS2D (NI, X(1,NJ_EU), Y(1,NJ_EU), FALSE, TTOTAL, TOUTER)

*     Lower sonic bubble edge location:

      TSONIC = TGRID(NI_EU) - X_SB;  ISL = NI_EU - NI_EU / 10

      CALL INTERVAL (NI_EU, TGRID, TSONIC, ONE, ISL)

      ISL = MIN (ISL, NI_EU - NI_EU / 20)
      TSONIC = TGRID(ISL)

*     Normal at lower sonic bubble edge: find the slope.

      CALL LCSFIT (NI_EU, TGRID, X, NEW, LOOSE, 1, TSONIC,
     >             XLINE(1), DXDT)
      CALL LCSFIT (NI_EU, TGRID, Y, NEW, LOOSE, 1, TSONIC,
     >             YLINE(1), DYDT)

      ANGLE = ATAN (DYDT / DXDT)
      DS = EIGHTH * DELTA
      XLINE(2) = XLINE(1) - DS * ABS (SIN (ANGLE))
      YLINE(2) = YLINE(1) + DS * SIGN (ONE, ANGLE) * COS (ANGLE)

*     Finding where the normal intersects the interim outer boundary then
*     redistributing the interim boundary leads to slight irregularities.
*     Therefore, just connect the boundaries at index ISL.

      XLINE(3) = XLINE(2) * TWO - XLINE(1);  XLINE(4) = X(ISL,NJ_EU)
      YLINE(3) = YLINE(2) * TWO - YLINE(1);  YLINE(4) = Y(ISL,NJ_EU)

      CALL CHORDS2D (4, XLINE, YLINE, FALSE, RLENGTH, TLINE)

*     For grid points along the lower sonic bubble edge, derive initial
*     and final increments from the those at the nose and tail:

      TTOTAL  = TGRID(NI_EU)
      D1SONIC = (TGRID(ISL) * D1NOSE + (TTOTAL - TGRID(ISL)) * D1TAIL) /
     >          TTOTAL
      D1SONIC = D1SONIC * RLENGTH

      TTOTAL  = TOUTER(NI_EU)
      D2SONIC = (TOUTER(I2) * D2NOSE + (TTOTAL - TOUTER(I2)) * D2TAIL) /
     >          TTOTAL
      D2SONIC = D2SONIC * RLENGTH

      TRADIAL(NJ_EU) = RLENGTH

      IF (D2SONIC > ZERO) THEN ! Two-sided

         CALL BLGRID (NJ_EU, D1SONIC, D2SONIC,
     >                NBLAYER, RBLAYER, TRADIAL, LUNCRT, IER)

         IF (IER /= 0) THEN
            WRITE (LUNCRT, '(/, A, I3)')
     >         ' BLGRID trouble at lower sonic bubble edge.  IER: ', IER
            GO TO 999
         END IF

      ELSE ! One-sided

         CALL EXPDIS5 (1, ZERO, RLENGTH, D1SONIC, NJ_EU, TRADIAL,
     >                 LUNITER)
      END IF

      CALL LCSFIT (4, TLINE, XLINE, NEW, TIGHT, NJ_EU, TRADIAL,
     >             XRADIAL, DERIVSJ)
      CALL LCSFIT (4, TLINE, YLINE, NEW, LOOSE, NJ_EU, TRADIAL,
     >             YRADIAL, DERIVSJ)

      DO J = 2, NJ_EU - 1
         X(ISL,J) = XRADIAL(J)
         Y(ISL,J) = YRADIAL(J)
      END DO


*     **** Repeat for upper outer boundary ****

*     Upper sonic bubble edge location:

      TSONIC = TGRID(NI_EU) + X_SB;  ISU = NI_EU + NI_EU / 10

      CALL INTERVAL (NI, TGRID, TSONIC, ONE, ISU)

      ISU = MAX (ISU, NI_EU + NI_EU / 20)
      TSONIC = TGRID(ISU)

      WRITE (LUNCRT, '(/, A, 2I5)')
     >   ' Sonic bubble indices: ', ISL, ISU

*     Normal at upper sonic bubble edge: find the slope.

      CALL LCSFIT (NI_EU, TGRID(NI_EU), X(NI_EU,1), NEW, LOOSE,
     >             1, TSONIC, XLINE(1), DXDT)
      CALL LCSFIT (NI_EU, TGRID(NI_EU), Y(NI_EU,1), NEW, LOOSE,
     >             1, TSONIC, YLINE(1), DYDT)

      ANGLE = ATAN (DYDT / DXDT)

      XLINE(2) = XLINE(1) - DS * ABS (SIN (ANGLE))
      YLINE(2) = YLINE(1) + DS * SIGN (ONE, ANGLE) * COS (ANGLE)
      XLINE(3) = XLINE(2) *  TWO  - XLINE(1);  XLINE(4) = X(ISU,NJ_EU)
      YLINE(3) = YLINE(2) *  TWO  - YLINE(1);  YLINE(4) = Y(ISU,NJ_EU)

      CALL CHORDS2D (4, XLINE, YLINE, FALSE, RLENGTH, TLINE)

      D2SONIC = ((TOUTER(NI) - TOUTER(I2))    * D2NOSE  +
     >           (TOUTER(I2) - TOUTER(NI_EU)) * D2TAIL) /
     >           (TOUTER(NI) - TOUTER(NI_EU))
      D2SONIC = D2SONIC * RLENGTH

      TRADIAL(NJ_EU) = RLENGTH

      IF (D2SONIC > ZERO) THEN ! Two-sided

         CALL BLGRID (NJ_EU, D1SONIC, D2SONIC,
     >                NBLAYER, RBLAYER, TRADIAL, LUNCRT, IER)

         IF (IER /= 0) THEN
            WRITE (LUNCRT, '(/, A, I3)')
     >         ' BLGRID trouble at upper sonic bubble edge.  IER: ', IER
            GO TO 999
         END IF

      ELSE ! One-sided

         CALL EXPDIS5 (1, ZERO, RLENGTH, D1SONIC, NJ_EU, TRADIAL,
     >                 LUNITER)
      END IF

      CALL LCSFIT (4, TLINE, XLINE, NEW, TIGHT, NJ_EU, TRADIAL,
     >             XRADIAL, DERIVSJ)
      CALL LCSFIT (4, TLINE, YLINE, NEW, LOOSE, NJ_EU, TRADIAL,
     >             YRADIAL, DERIVSJ)

      DO J = 2, NJ_EU - 1
         X(ISU,J) = XRADIAL(J)
         Y(ISU,J) = YRADIAL(J)
      END DO


*     Starting guess for the interior grid in 4 zones:
*     ***************************************************

      CALL TFI2D (NI, 1,     ISL,   1, NJ_EU, X, Y, ARCS)
      CALL TFI2D (NI, ISL,   NI_EU, 1, NJ_EU, X, Y, ARCS)
      CALL TFI2D (NI, NI_EU, ISU,   1, NJ_EU, X, Y, ARCS)
      CALL TFI2D (NI, ISU,   NI,    1, NJ_EU, X, Y, ARCS)

      IF (SMOOTH) THEN

         CALL ELLIP2D (NI, NJ_EU, 1,   ISL,    1, NJ_EU, X, Y, ARCS,
     >                 ITERATIONS, BGMODE, FGMODE, SPMODE,
     >                 ITMAX, ITFLOAT, ITFREEZE, JLAYER,
     >                 CONV2D, CONVMIN, DMAX2D, OMG2D,
     >                 POWERI, POWERJ, EXPI1, EXPI2, EXPJ1, EXPJ2,
     >                 FGLIMIT, URFG, URFLOAT)

         FGSWAP = FGMODE
         FGSWAP(2:2) = 'Y' ! Force integer-based orthogonality at stag. pt. I

         CALL ELLIP2D (NI, NJ_EU, ISL, NI_EU,  1, NJ_EU, X, Y, ARCS,
     >                 ITERATIONS, BGMODE, FGSWAP, SPMODE,
     >                 ITMAX, ITFLOAT, ITFREEZE, JLAYER,
     >                 CONV2D, CONVMIN, DMAX2D, OMG2D,
     >                 POWERI, POWERJ, EXPI1, EXPI2, EXPJ1, EXPJ2,
     >                 FGLIMIT, URFG, URFLOAT)

         BGSWAP      = BGMODE
         BGSWAP(1:1) = BGMODE(2:2)
         BGSWAP(2:2) = BGMODE(1:1)

         FGSWAP      = FGMODE
         FGSWAP(1:1) = 'Y'  ! Force integer-based orthogonality at stag. pt. I
         FGSWAP(2:2) = FGMODE(1:1)

         CALL ELLIP2D (NI, NJ_EU, NI_EU, ISU, 1, NJ_EU, X, Y, ARCS,
     >                 ITERATIONS, BGSWAP, FGSWAP, SPMODE,
     >                 ITMAX, ITFLOAT, ITFREEZE, JLAYER,
     >                 CONV2D, CONVMIN, DMAX2D, OMG2D,
     >                 -POWERI, POWERJ, EXPI2, EXPI1, EXPJ1, EXPJ2,
     >                 FGLIMIT, URFG, URFLOAT)

         FGSWAP(1:1) = FGMODE(2:2)

         CALL ELLIP2D (NI, NJ_EU, ISU,   NI,  1, NJ_EU, X, Y, ARCS,
     >                 ITERATIONS, BGSWAP, FGSWAP, SPMODE,
     >                 ITMAX, ITFLOAT, ITFREEZE, JLAYER,
     >                 CONV2D, CONVMIN, DMAX2D, OMG2D,
     >                 -POWERI, POWERJ, EXPI2, EXPI1, EXPJ1, EXPJ2,
     >                 FGLIMIT, URFG, URFLOAT)
      END IF

  999 RETURN

*     Internal procedure for subroutine HB_EULER:

      CONTAINS

*        -----------------------------------------------------------------------

         SUBROUTINE SHOCK_ANGLE ()

*        Solve eq. 2.16, p. 36 for limiting shock angle BETA.

*        -----------------------------------------------------------------------

         USE TRIGD
         IMPLICIT NONE

*        Local constants:

         INTEGER,   PARAMETER :: MAXFUN = 30    ! Max. # zero-finder fn. evals.
         REAL,      PARAMETER :: TOL = 1.E-7    ! No need for full precision
         REAL,      PARAMETER :: BETAMIN = 0.1  ! Search range for zero finder
         REAL,      PARAMETER :: BETAMAX = 1.5  ! (radians)
         CHARACTER, PARAMETER :: SUBNAME * 11 = 'SHOCK_ANGLE'

*        Local variables:

         INTEGER :: ISTAT, NUMFUN
         REAL    :: F, HOLD(13), TAN_THETA

*        Execution:

         TAN_THETA = TAND (THETA)

         ISTAT = 2          ! Initialize the iteration
         NUMFUN = MAXFUN

 10      CONTINUE

            CALL ZERORC (BETAMIN, BETAMAX, BETA, F, TOL, NUMFUN,
     >                   SUBNAME, LUNITER, HOLD, ISTAT)

            IF (ISTAT .LT. -1) THEN      ! Fatal error

               WRITE (LUNCRT, '(/, 2A)')
     >            ' Shock angle calculation failed:',
     >            ' proceeding with THETA*(gamma + 1)/2.'
               BETA = THETA * (GAMMA + ONE) * HALF * (PI / 180.)

            ELSE IF (ISTAT < 0) THEN  ! Iteration limit reached

               WRITE (LUNCRT, '(/, A)')
     >         ' Shock angle iteration limit reached; proceeding.'

            ELSE IF (ISTAT > 0) THEN  ! Evaluate the function

               F = TWO * TAN_THETA * TAN (BETA) -
     >             (FMSQ * (SIN (BETA))**2 - ONE) /
     >             (FMSQ * (COS (BETA*TWO) + GAMMA) + ONE)
               GO TO 10

            ! ELSE ISTAT = 0 - a zero has been found.

            END IF

         END SUBROUTINE SHOCK_ANGLE

      END SUBROUTINE HB_EULER

********************************************************************************
*
      SUBROUTINE HB_NS (LUNCRT, NI, NI_EU, NJ_EU, NJ_NS,
     >                  NBLAYERNS, RBLAYERNS,
     >                  BODY_LENGTH, D1NOSENS, D1TAILNS,
     >                  D2NOSENS, D2TAILNS, ITERATIONS, X, Y, XNS, YNS)
*
*        Redistribute the radial lines of a 2-space Euler-type hypersonic body
*     grid to achieve Navier-Stokes-type clustering.
*
*     10/31/02  DAS  Initial implementation.
*
*     David Saunders, ELORET/NASA Ames Research Center, Moffett Field, CA.
*
********************************************************************************

      IMPLICIT NONE

*     Arguments:

      INTEGER, INTENT (IN) ::
     >   LUNCRT,             ! Unit for iteration printing; < 0 suppresses it
     >   NI,                 ! Full wrap-around dimension
     >   NI_EU,              ! Streamwise dimension of upper & lower halves
     >   NJ_EU,              ! Radial dimension of input grid
     >   NJ_NS,              ! Output dimensions are NI x NJ_NS
     >   NBLAYERNS           ! Desired # points in boundary layer (if D2* > 0.)

      REAL, INTENT (IN) ::
     >   BODY_LENGTH,        ! Body length used to scale D1*
     >   RBLAYERNS,          ! Geometric growth rate in boundary layer
     >   D1NOSENS, D1TAILNS, ! Desired increments at the surface ...
     >   D2NOSENS, D2TAILNS  ! ... and at the upstream boundary, as in D2_EULER;
                             ! D2* < 0. gives 1-sided stretching, else D2* is a
                             ! fraction of the radial length, with linear
                             ! variation between the nose and tail
      LOGICAL, INTENT (IN) ::
     >   ITERATIONS          ! T turns on Vinokur-type iteration printing

      REAL, INTENT (IN) ::
     >   X(NI, NJ_EU),       ! Input Euler-type grid
     >   Y(NI, NJ_EU)

      REAL, INTENT (OUT), DIMENSION (NI, NJ_NS) ::
     >   XNS, YNS            ! Redistributed grid

*     Local constants:

      REAL, PARAMETER ::
     >   ONE = 1., ZERO = 0.

      LOGICAL, PARAMETER ::
     >   FALSE = .FALSE., TRUE = .TRUE.

      CHARACTER, PARAMETER ::
     >   LOOSE * 1 = 'B'

*     Automatic arrays:

      REAL
     >   D1(NI), D2(NI), DERIVS(NJ_NS),
     >   TBLEND(NJ_NS), TBODY(NI_EU),  TEULER(NJ_EU), TKEEP(NJ_NS),
     >   TLOWER(NI_EU), TSPOKE(NJ_NS), TSTAG(NJ_NS),
     >   XBLEND(NJ_NS), XEULER(NJ_EU), XSPOKE(NJ_NS),
     >   YBLEND(NJ_NS), YEULER(NJ_EU), YSPOKE(NJ_NS)

*     Local variables:

      INTEGER
     >   I, IER, J, K, LUNITER, NBLEND

      REAL
     >   B, PIBY2, R, RLENGTH, TOTAL, TOTSTAG

*     Execution:

      IF (ITERATIONS) THEN
         LUNITER =  LUNCRT
      ELSE
         LUNITER = -LUNCRT
      END IF

*     Establish the specified initial increments at the lower surface:

      CALL CHORDS2D (NI_EU, X, Y, FALSE, TOTAL, TLOWER)

      DO I = 1, NI_EU
         TBODY(I) = TOTAL - TLOWER(I)
      END DO

*     Parabolic variation of D1 from nose to tail (D1 ** 2 ~ a t + b):

      CALL INVERSE_PARABOLA (NI_EU, TBODY, ZERO, D1NOSENS * BODY_LENGTH,
     >                       TOTAL, D1TAILNS * BODY_LENGTH, D1)

*     Initial increments at the upper surface:

      CALL CHORDS2D (NI_EU, X(NI_EU,1), Y(NI_EU,1), FALSE, TOTAL, TBODY)

      CALL INVERSE_PARABOLA (NI_EU, TBODY, ZERO, D1NOSENS * BODY_LENGTH,
     >                       TOTAL, D1TAILNS * BODY_LENGTH, D1(NI_EU))


*     Interim restretching (but it probably produces cusps off the stag. pt.):

      IF (D2NOSENS > ZERO) THEN ! Two-sided stretching (geometric + Vinokur)

*        Establish the specified final increments at the lower outer boundary:

         CALL CHORDS2D (NI_EU, X(1,NJ_EU), Y(1,NJ_EU), TRUE, TOTAL,
     >                  TBODY) ! Reuse TBODY(*) for normalized outer arc lengths

*        Linear variation of D2 between tail and nose:

         DO I = 1, NI_EU
            RLENGTH = SQRT ((X(I,1) - X(I,NJ_EU)) ** 2 +  ! Near enough
     >                      (Y(I,1) - Y(I,NJ_EU)) ** 2)
            R = TBODY(I)
            D2(I) = RLENGTH * ((ONE - R) * D2TAILNS + R * D2NOSENS)
         END DO

*        Final increments at the upper outer boundary:

         CALL CHORDS2D (NI_EU, X(NI_EU,NJ_EU), Y(NI_EU,NJ_EU), TRUE,
     >                  TOTAL, TBODY) ! Reuse TBODY(*) for normalized outer arcs

         J = 1
         DO I = NI_EU, NI
            RLENGTH = SQRT ((X(I,1) - X(I,NJ_EU)) ** 2 +  ! Near enough
     >                      (Y(I,1) - Y(I,NJ_EU)) ** 2)
            R = TBODY(J)
            D2(I) = RLENGTH * ((ONE - R) * D2NOSENS + R * D2TAILNS)
            J = J + 1
         END DO

         TSPOKE(1) = ZERO

         DO I = 1, NI

            DO J = 1, NJ_EU
               XEULER(J) = X(I,J)
               YEULER(J) = Y(I,J)
            END DO

            CALL CHORDS2D (NJ_EU, XEULER, YEULER, FALSE, TOTAL, TEULER)

            TSPOKE(NJ_NS) = TOTAL

            CALL BLGRID (NJ_NS, D1(I), D2(I),
     >                   NBLAYERNS, RBLAYERNS, TSPOKE, LUNCRT, IER)

            IF (IER /= 0) THEN
               WRITE (LUNCRT, '(/, (A, I3))')
     >            ' BLGRID trouble restretching at I = ', I,
     >            '.  IER: ', IER, '.  Revert to 1-sided.'

               CALL EXPDIS5 (1, ZERO, TOTAL, D1(I), NJ_NS, TSPOKE,
     >                       LUNITER)

            END IF

            CALL LCSFIT (NJ_EU, TEULER, XEULER, TRUE, LOOSE, NJ_NS,
     >                   TSPOKE, XSPOKE, DERIVS)

            CALL LCSFIT (NJ_EU, TEULER, YEULER, TRUE, LOOSE, NJ_NS,
     >                   TSPOKE, YSPOKE, DERIVS)

            DO J = 1, NJ_NS
               XNS(I,J) = XSPOKE(J)
               YNS(I,J) = YSPOKE(J)
            END DO

         END DO

      ELSE ! One-sided stretching

         DO I = 1, NI

            DO J = 1, NJ_EU
               XEULER(J) = X(I,J)
               YEULER(J) = Y(I,J)
            END DO

            CALL CHORDS2D (NJ_EU, XEULER, YEULER, FALSE, TOTAL, TEULER)

            CALL EXPDIS5 (1, ZERO, TOTAL, D1(I), NJ_NS, TSPOKE, LUNITER)

            CALL LCSFIT (NJ_EU, TEULER, XEULER, TRUE, LOOSE, NJ_NS,
     >                   TSPOKE, XSPOKE, DERIVS)

            CALL LCSFIT (NJ_EU, TEULER, YEULER, TRUE, LOOSE, NJ_NS,
     >                   TSPOKE, YSPOKE, DERIVS)

            DO J = 1, NJ_NS
               XNS(I,J) = XSPOKE(J)
               YNS(I,J) = YSPOKE(J)
            END DO

         END DO

      END IF

*     Smooth out the inevitable cusp at the stag. pt. radial line by
*     blending its relative distribution with the interim distribution
*     in a nonlinear way for a reasonable range of I either side:

*     Normalized arcs at stag. point:

      DO J = 1, NJ_NS
         XSPOKE(J) = XNS(NI_EU,J)
         YSPOKE(J) = YNS(NI_EU,J)
      END DO

      CALL CHORDS2D (NJ_NS, XSPOKE, YSPOKE, TRUE, TOTSTAG, TSTAG)

*     Normalized arcs at lower station that will remain untouched:

***   NBLEND = MAX (NI_EU / 3, 40)  ! Make it an input?
***   NBLEND = MIN (NI_EU / 2, NBLEND)

      NBLEND = NI_EU - 1

      I = NI_EU - NBLEND

      DO J = 1, NJ_NS
         XSPOKE(J) = XNS(I,J)
         YSPOKE(J) = YNS(I,J)
      END DO

      CALL CHORDS2D (NJ_NS, XSPOKE, YSPOKE, TRUE, TOTAL, TKEEP)

*     Normalized arcs along body where we're blending:

      CALL CHORDS2D (NBLEND + 1, XNS(I,1), YNS(I,1), TRUE, TOTAL, TBODY)

      PIBY2 = ACOS (ZERO)

      K = 2
      DO I = NI_EU - NBLEND + 1, NI_EU - 1

         DO J = 1, NJ_NS
            XSPOKE(J) = XNS(I,J)
            YSPOKE(J) = YNS(I,J)
         END DO

         CALL CHORDS2D (NJ_NS, XSPOKE, YSPOKE, FALSE, TOTAL, TSPOKE)

         R = TBODY(K)
         R = (SIN ((R + ONE) * PIBY2)) ** 5 ! Favors stag. pt. region smoothly
         K = K + 1                          ! See RCOS function in SFEVAL

         DO J = 2, NJ_NS
            TBLEND(J) = (R * TKEEP(J) + (ONE - R) * TSTAG(J)) * TOTAL
         END DO

         CALL LCSFIT (NJ_NS, TSPOKE, XSPOKE, TRUE, LOOSE, NJ_NS - 1,
     >                TBLEND(2), XBLEND(2), DERIVS)

         CALL LCSFIT (NJ_NS, TSPOKE, YSPOKE, TRUE, LOOSE, NJ_NS - 1,
     >                TBLEND(2), YBLEND(2), DERIVS)

         DO J = 2, NJ_NS
            XNS(I,J) = XBLEND(J)
            YNS(I,J) = YBLEND(J)
         END DO

      END DO

*     Repeat for above the stag. point:

      I = NI_EU + NBLEND

      DO J = 1, NJ_NS
         XSPOKE(J) = XNS(I,J)
         YSPOKE(J) = YNS(I,J)
      END DO

      CALL CHORDS2D (NJ_NS, XSPOKE, YSPOKE, TRUE, TOTAL, TKEEP)

*     Normalized arcs along body where we're blending:

      I = NI_EU
      CALL CHORDS2D (NBLEND + 1, XNS(I,1), YNS(I,1), TRUE, TOTAL, TBODY)

      K = 2
      DO I = NI_EU + 1, NI_EU + NBLEND - 1

         DO J = 1, NJ_NS
            XSPOKE(J) = XNS(I,J)
            YSPOKE(J) = YNS(I,J)
         END DO

         CALL CHORDS2D (NJ_NS, XSPOKE, YSPOKE, FALSE, TOTAL, TSPOKE)

         R = TBODY(K)
         R = (SIN (R * PIBY2)) ** 5 ! Favors stag. pt. region smoothly
         K = K + 1                  ! See LCOS function in SFEVAL

         DO J = 2, NJ_NS
            TBLEND(J) = ((ONE - R) * TSTAG(J) + R * TKEEP(J)) * TOTAL
         END DO

         CALL LCSFIT (NJ_NS, TSPOKE, XSPOKE, TRUE, LOOSE, NJ_NS - 1,
     >                TBLEND(2), XBLEND(2), DERIVS)

         CALL LCSFIT (NJ_NS, TSPOKE, YSPOKE, TRUE, LOOSE, NJ_NS - 1,
     >                TBLEND(2), YBLEND(2), DERIVS)

         DO J = 2, NJ_NS
            XNS(I,J) = XBLEND(J)
            YNS(I,J) = YBLEND(J)
         END DO

      END DO

      END SUBROUTINE HB_NS

********************************************************************************
*
      SUBROUTINE HB_SAVE (LUN, FORMATTED, TWO_D, ALPHA, NI, NJ, X, Y, Z)
*
********************************************************************************

      IMPLICIT NONE

*     Arguments:

      INTEGER, INTENT (IN) :: LUN        ! Logical unit (assumed open)
      LOGICAL, INTENT (IN) :: FORMATTED, ! T or F
     >                        TWO_D      ! T = (x,y) only; F = (x,y,z)
      REAL,    INTENT (IN) :: ALPHA      ! 0 => add axisymmetric.gu for DPLR
      INTEGER, INTENT (IN) :: NI, NJ     ! Grid dimensions
      REAL,    INTENT (IN) :: X(NI, NJ)  ! Grid coordinates
      REAL, INTENT (INOUT) :: Y(NI, NJ)  ! Allow for forcing Y = Ystag (AoA = 0)
      REAL,    INTENT (IN) :: Z(NI, NJ)  ! Grid coordinates

*     Local variables:

      INTEGER :: I, NUMI
***   REAL    :: Y_STAG_PT

*     Execution:

      NUMI = (NI + 1) / 2

*     Ensure y = 0 exactly in front of stag. pt. for ballistic entry cases,
*     assuming geometric symmetry.  (NO:  Do it before the elliptic smoothing.)

***   IF (ALPHA == 0.) THEN
***      Y_STAG_PT = Y(NUMI,1)
***      IF (ABS (Y_STAG_PT) < 1.E-16) Y_STAG_PT = 0.
***      Y(NUMI,:) = Y_STAG_PT
***   END IF

      IF (FORMATTED) THEN
         IF (TWO_D) THEN
            WRITE (LUN, '(I1)') 1
            WRITE (LUN, '(I3, I4)') NI, NJ
            WRITE (LUN, '(1P, 5E23.15)') X, Y
         ELSE ! For possible Gridgen access
            WRITE (LUN, '(I1)') 1
            WRITE (LUN, '(I3, I4, I3)') NI, NJ, 1
            WRITE (LUN, '(1P, 5E23.15)') X, Y
            WRITE (LUN, '(26F3.0)') Z ! Assume all zeros
         END IF
      ELSE
         IF (TWO_D) THEN
            WRITE (LUN) 1
            WRITE (LUN) NI, NJ
            WRITE (LUN) X, Y
         ELSE
            WRITE (LUN) 1
            WRITE (LUN) NI, NJ, 1
            WRITE (LUN) X, Y, Z
         END IF
      END IF

      IF (ALPHA == 0.) THEN
         CLOSE (LUN)
         OPEN  (LUN, FILE='axisymmetric.gu', FORM='UNFORMATTED',
     >          STATUS='UNKNOWN')
         IF (TWO_D) THEN
            WRITE (LUN) 1
            WRITE (LUN) NUMI, NJ
            WRITE (LUN) X(NUMI:NI,:), Y(NUMI:NI,:)
         ELSE
            WRITE (LUN) 1
            WRITE (LUN) NUMI, NJ, 1
            WRITE (LUN) X(NUMI:NI,:), Y(NUMI:NI,:), Z(NUMI:NI,:)
         END IF
      END IF

      END SUBROUTINE HB_SAVE
