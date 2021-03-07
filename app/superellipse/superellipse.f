C+------------------------------------------------------------------------------
C
      PROGRAM SUPERELLIPSE
C
C     PURPOSE:
C
C        SUPERELLIPSE generates planform and thickness data for Oblique All-
C     Wing transport (OAW) designs, with MODIFIED elliptic variation of the
C     chord and thickness outboard of the specified span stations.
C
C        Parameter POW is the exponent in the sub/super-ellipse formula, here
C     assumed to be the same for both axes. Values less than 2.0 are "sharper"
C     than an ellipse; greater than 2.0 "blunter." Note that the exact formula
C     of the original has been replaced by numerical approximations to the
C     figures' area and volume.
C
C        This version of the program also generates a QPLOT file for graphic
C     output.
C
C        Program R22OPT uses the output from this program to produce a
C     complete geometry definition (with lofted wing sections and possibly
C     some perturbations applied by its application-specific PERTRB routine).
C
C     PROCEDURES:
C
C        ARBDIS    Arbitrary-shape 1-D grid distributions
C        INTERVAL  Search utility (also used at a lower level)
C        LCSQUAD   Numerical integration (area, volume)
C        QINTRP    Quadratic interpolation (2 points, 1 slope)
C        READER    Prompting utility
C        SMOOTHX   Ensures inclusion of specified grid stations
C
C     HISTORY:
C
C        11/28/92   D.Saunders   Initial elliptic adaptation of parabolic/
C                                linear taper program TAPER.
C        01/11/93   R.Kennelly   Modified from ELLIPTIC-OAW. Try super-
C                                (actually sub-) ellipse functional form
C                                with plot file as output. Integrate
C                                numerically.
C        01/12/93   DAS/RAK      Ensure that forward tip wake clears the
C                                trailing edge.  Linear tip trailing edge
C                                (larger chords) should also help R22OPT's
C                                poor behavior with highly curved tips.
C        01/13/93      "         R22OPT still has trouble at the tips.
C                                Make the tip chord ~0.3 x root chord.
C        01/14/93      "         "THICK" was being tapered where "TMAX" was
C                                probably intended.  However, tapering TMAX
C                                as for CHORD gives constant T/C.  Therefore,
C                                scale TMAX beyond cabin such that T/C is
C                                as specified (e.g. half the root T/C) at
C                                the near-tip station (tip-ratio*tip chord
C                                inboard).
C        01/15/93      "         A linearly-varying scale factor on T to
C                                constrain near-tip thickness is not good
C                                (T/C is rectilinear).  Provide independent
C                                superellipse exponents, and allow these to
C                                vary linearly with span station.
C        01/19/93      "         1. Add elliptic tip cap code from TAPER.
C                                Use of complete composite superellipses
C                                for the non-squared tip case is eschewed
C                                for fear of producing sections that are
C                                too thin towards the tip unless the volume
C                                grows considerably.  Instead, we try to
C                                match the squared-tip case (analyzed in
C                                R22OPT) as closely as possible.
C                                2. Provide a smoother spanwise distribution
C                                for the round-tip case (elliptic shape function
C                                for ARBDIS) while still including end of cabin
C                                stations, etc., via SMOOTHX.
C        01/31/93      "         Added an option to vary the exponents
C                                quadratically (zero slope at the tip).
C        02/03/93      "         Added the alternative quadratic option
C                                (zero slope at Z1).
C        02/08/93      "         More precise application of taper ratio and
C                                tip cap ratio (at the Z2 station, not ZT).
C                                Prompt for tip ratio.
C        02/15/93      "         Allow for larger spans by prompting.
C        02/18/93      "         Added option for semi-elliptic trailing tip
C                                following Ilan Kroo's suggestion to carry
C                                the straight leading edge all the way out
C                                to the tip for maximum swept span.  (What
C                                to do at the forward tip is more complex.)
C        03/01/93      "         R22OPT's more precise thickness calculation
C                                means root T/C is slightly greater than
C                                that indicated by PROFILE.
C        03/08/93      "         Start the taper from 550" rather than 700".
C        03/19/93      "         Taper the tip thickness parabolically now
C                                since elliptic is too blunt.
C-------------------------------------------------------------------------------

      USE TRIGD
      IMPLICIT NONE

C     Constants defining the planform:
C     NOTE: Some of these assume planform OAW-1 as baseline.

      INTEGER
     >   MXZ, NCAP
      REAL
     >   B1, B2, CROOT, S, SWEEP, XC, Z1
      PARAMETER
     >  (MXZ   = 101,       ! Maximum number of span stations allowed for
     >   NCAP  = 2,         ! Number of stations on the square tip beyond Z2
     >   CROOT = 55.8,      ! Root chord (ft., as for other dimensions)
     >   SWEEP = 68.0,      ! Wing sweep (degrees)
     >   S     = 0.32,      ! X/C line which is the major axis of the ellipses
     >   XC    = S * CROOT, ! X line which is at constant X/C for all sections
     >   Z1    = 550./12.,  ! Station at which chord & thickness start tapering
     >   B1    = S * CROOT, ! Minor axis of leading edge ellipses
     >   B2    = CROOT - B1)! ............. trailing ............

C     Parameters defining the thickness distribution (YLE = 0. assumed).
C     NOTE: This set assumes airfoil OAW 60-085-16 in center section.
C     Also: the same airfoil is assumed to apply at all stations,
C           with Y coordinates scaled by the THICK parameter.

      REAL
     >   AREARATIO, T1, TCROOT
      PARAMETER
     >  (TCROOT = 0.160028,           ! Max. thickness/chord for cabin section
     >   AREARATIO = 0.1087446/TCROOT,! Area frctn. of section wrt bounding rect
     >   T1     = TCROOT * CROOT)     ! Thickness of cabin section = minor axis

C     Other constants:

      INTEGER
     >   LUNCRT, LUNCTL, LUNKBD, LUNOUT, LUNRND, LUNSQR
      REAL
     >   ONE, ZERO
      CHARACTER
     >   CTLFILE * 16, OUTFILE * 16, RNDFILE * 13, SQRFILE * 14
      PARAMETER
     >  (LUNCRT  = 6,
     >   LUNCTL  = 1,
     >   LUNKBD  = 5,
     >   LUNOUT  = 2,
     >   LUNRND  = 3,
     >   LUNSQR  = 4,
     >   CTLFILE = 'superellipse.ctl',  ! QPLOT control file
     >   OUTFILE = 'superellipse.out',  ! R22OPT geometry data
     >   RNDFILE = 'round-tip.plt',     ! QPLOTable data, multi-column
     >   SQRFILE = 'square-tip.plt',    ! ............................
     >   ZERO    = 0.0E+0,
     >   ONE     = 1.0E+0)

C     Variables:

      INTEGER
     >   IER, IPLAN, IPLAN1, IPLAN2, IVARY, IVOL, IVOL1, IVOL2,
     >   K, K1, K2, KL, KMID, KS, NZ
      REAL
     >   A, AX, AREA (MXZ), BX, C2, CHORD (MXZ), CHORDCAP, CHSLOPE,
     >   DZ, LESLOPE, PLAN, POWC, POWC1, POWC2, POWT, POWT1, POWT2,
     >   QSLOPE, TEANGLE, TESLOPE, THICK (MXZ), THSLOPE, TAPERCAP,
     >   TAPERRATIO, TCCAP, TIPRATIO, TMAX (MXZ), TMSLOPE, TOVERC (MXZ),
     >   TSCALE, VOL, WAKESLOPE, WORK ((1 + MXZ/2) * 5), XLE (MXZ),
     >   XLE1, XLECAP, XT, XTECAP, XTEDGE, YCSLOPE, YCTL (MXZ), Z (MXZ),
     >   Z2, ZC, ZCTL (MXZ), ZE, ZELLIPSE, ZFIXED (3), ZSTRAIGHT, ZT
      CHARACTER
     >   AREATXT * 7, DAY * 9, VARYCASE (3) * 32, VOLUMETXT * 7
      LOGICAL
     >   CR, EOF, SEMI_SQUARE_LEFT_TIP

C     Procedures:

      EXTERNAL
     >   ARBDIS, INTERVAL, LCSQUAD, QINTRP, READI, READR, READY,
     >   SMOOTHX

C     Storage:

      DATA
     >   VARYCASE / 'linear', 'quadratic (zero slope outboard)',
     >              'quadratic (zero slope inboard)'/

C     Statement functions:

      REAL
     >   LINEAR,           ! Linear interpolation/extrapolation
     >   ABSCISSA, ORDINATE, SLOPE, TARGET,
     >   ELLIPSE,          ! First-quadrant superellipse
     >   AXIS, BXIS, POW, X

      LINEAR (ABSCISSA, ORDINATE, SLOPE, TARGET) =
     >   ORDINATE + SLOPE * (TARGET - ABSCISSA)

      ELLIPSE (AXIS, BXIS, POW, X) =
     >   BXIS * (ONE - (X / AXIS) ** POW) ** (ONE / POW)


C     Execution:
C     ----------

      ZT = 406.8
      CALL READR (LUNCRT, 'Wing span?  [406.8 feet]:  ',
     >   LUNKBD, ZT, CR, EOF)
      ZT = ZT * 0.5
      A = ZT - Z1        ! Major axis of all ellipses except tip cap

   90 NZ = 101
      CALL READI (LUNCRT,
     >   'Number of span stations (an odd number)?  [101]: ',
     >   LUNKBD, NZ, CR, EOF)
      IF (EOF) GO TO 999
      IF (NZ .GT. MXZ) GO TO 90

C     R22OPT suffers with too small a tip chord.

      TAPERRATIO = 0.3
      CALL READR (LUNCRT,
     >   'Taper ratio (tip-cap chord / root chord)?  [0.3]: ',
     >   LUNKBD, TAPERRATIO, CR, EOF)

      TIPRATIO = 0.5
      CALL READR (LUNCRT,
     >   'Tip ratio (tip-cap span / tip-cap chord)?  [0.5]: ',
     >   LUNKBD, TIPRATIO, CR, EOF)

      C2 = CROOT * TAPERRATIO
      Z2 = ZT - C2 * TIPRATIO

      WRITE (LUNCRT, '(A, F7.3)')
     >   ' Corresponding tip cap station where T/C is monitored: ',
     >   Z2

      SEMI_SQUARE_LEFT_TIP = .TRUE.
      CALL READY (LUNCRT,
     >   'Semi-elliptic left tip cap?  ' //
     >   '(Y/N; No = fully elliptic)  [Y]: ',
     >   LUNKBD, SEMI_SQUARE_LEFT_TIP, CR, EOF)

      WRITE (LUNCRT, '(A)')
     >   '0Superellipse exponent variation options:', ' '
      WRITE (LUNCRT, '(I6, '':  '', A)') (K, VARYCASE (K), K = 1, 3)
      WRITE (LUNCRT, '(A)')

      IVARY = 2
      CALL READI (LUNCRT, 'Pick one.  [2]: ',
     >   LUNKBD, IVARY, CR, EOF)

      OPEN (UNIT=LUNOUT, FILE=OUTFILE, STATUS='NEW')
      OPEN (UNIT=LUNCTL, FILE=CTLFILE, STATUS='NEW')
      OPEN (UNIT=LUNRND, FILE=RNDFILE, STATUS='NEW')
      OPEN (UNIT=LUNSQR, FILE=SQRFILE, STATUS='NEW')

      WRITE (LUNCRT, '(A)')
     >   '0Enter ^Z to terminate the following trial-and-error loop.',
     >   ' '

  100 CONTINUE

C     Start of squared-tip retry loop for checking volume, tip-cap T/C, etc.
C     ----------------------------------------------------------------------

      WRITE (LUNCRT, '(A)')
     >   '$Superellipse exponents for chord (inboard, outboard): '
      READ (LUNKBD, *, ERR=100, END=600) POWC1, POWC2
      WRITE (LUNCRT, '(A)')
     >   '$Superellipse exponents for Tmax  (inboard, outboard): '
      READ (LUNKBD, *, ERR=100, END=600) POWT1, POWT2

C     Determine a reasonable index for end-of-cabin station, Z1, as will be
C     used by R22OPT with squared-off tips.  (The rounded-tip case has a
C     much smoother spanwise distribution once the break station has been
C     established.)  The 0.5 fudge for K1 gives fewer stations for the straight
C     cabin (better resolution of the nonlinear panel) than a uniform spanwise
C     distribution would.  Note that the outer linear panel's chord taper
C     gradient is continuous only in the sense of a simple linear extrapolation
C     from the last two elliptic stations.  The slopes are small enough that
C     this is surely acceptable, although the result for "ZSTRAIGHT" varies
C     slightly with NZ as a consequence.  The more stations between Z1 and Z2
C     for the square-tip case the better.

      KMID = 1 + NZ / 2
      K1 = KMID + (0.5 * Z1 / ZT) * KMID
      K1 = MAX (K1, KMID + 2)

C     Right cabin panel:

      DZ = Z1 / (K1 - KMID)

      DO K = KMID, K1
         Z (K) = (K - KMID) * DZ
         XLE (K) = ZERO
         TMAX (K) = T1
         AREA (K) = T1 * CROOT * AREARATIO
         CHORD (K) = CROOT
         THICK (K) = ONE
         TOVERC (K) = TCROOT
      END DO

C     Right outer panel, Z1 to Z2: leave NCAP stations on the tip cap.

      DZ = (Z2 - Z1) / ((NZ - NCAP) - K1)
      THSLOPE = (POWT2 - POWT1) / (ZT - Z1) ! Slope of linear variation of POWT
      CHSLOPE = (POWC2 - POWC1) / (ZT - Z1) ! ............................ POWC

      DO K = K1 + 1, NZ - NCAP
         Z (K) = (K - K1) * DZ + Z1

         IF (IVARY .EQ. 1) THEN  ! Linear variation of exponents

            POWC = LINEAR  (Z1, POWC1, CHSLOPE, Z (K))
            POWT = LINEAR  (Z1, POWT1, THSLOPE, Z (K))

         ELSE IF (IVARY .EQ. 2) THEN  ! Quadratic (zero slope outboard)

            CALL QINTRP (1, Z (K), ZT, POWC2, ZERO, Z1, POWC1,
     >                   QSLOPE, POWC)
            CALL QINTRP (1, Z (K), ZT, POWT2, ZERO, Z1, POWT1,
     >                   QSLOPE, POWT)

         ELSE IF (IVARY .EQ. 3) THEN  ! Quadratic (zero slope inboard)

            CALL QINTRP (1, Z (K), Z1, POWC1, ZERO, ZT, POWC2,
     >                   QSLOPE, POWC)
            CALL QINTRP (1, Z (K), Z1, POWT1, ZERO, ZT, POWT2,
     >                   QSLOPE, POWT)

         END IF

         ZELLIPSE  = Z (K) - Z1
         XLE (K)   = B1 - ELLIPSE (A, B1, POWC, ZELLIPSE)
         CHORD (K) = B1 + ELLIPSE (A, B2, POWC, ZELLIPSE) - XLE (K)
         TMAX (K)  =      ELLIPSE (A, T1, POWT, ZELLIPSE)
      END DO

C     Right outer tip, Z2 to ZT (same as above, but DZ differs):

      DZ = (ZT - Z2) / NCAP
      DO K = NZ - NCAP + 1, NZ
         Z (K) = (K - (NZ - NCAP)) * DZ + Z2

         IF (IVARY .EQ. 1) THEN
            POWC = LINEAR  (Z1, POWC1, CHSLOPE, Z (K))
            POWT = LINEAR  (Z1, POWT1, THSLOPE, Z (K))

         ELSE IF (IVARY .EQ. 2) THEN
            CALL QINTRP (1, Z (K), ZT, POWC2, ZERO, Z1, POWC1,
     >                   QSLOPE, POWC)
            CALL QINTRP (1, Z (K), ZT, POWT2, ZERO, Z1, POWT1,
     >                   QSLOPE, POWT)

         ELSE IF (IVARY .EQ. 3) THEN
            CALL QINTRP (1, Z (K), Z1, POWC1, ZERO, ZT, POWC2,
     >                   QSLOPE, POWC)
            CALL QINTRP (1, Z (K), Z1, POWT1, ZERO, ZT, POWT2,
     >                   QSLOPE, POWT)
         END IF

         ZELLIPSE  = Z (K) - Z1
         XLE (K)   = B1 - ELLIPSE (A, B1, POWC, ZELLIPSE)
         CHORD (K) = B1 + ELLIPSE (A, B2, POWC, ZELLIPSE) - XLE (K)
         TMAX (K)  =      ELLIPSE (A, T1, POWT, ZELLIPSE)
      END DO


C     Fudge the tip trailing edge to avoid forward tip wake problems.
C     Find the span station pair to extrapolate linearly from.
C     Treat the specified tip-cap chord as lower bound.

      WAKESLOPE = TAND (SWEEP - 90.)   ! More negative is bad for TESLOPE

      DO K = NZ - 1, KMID + 1, -1   ! Search from the tip backwards
         KS = K
         DZ = Z (K + 1) - Z (K)
         XTEDGE  = XLE (K) + CHORD (K)
         TESLOPE = ((XLE (K + 1) + CHORD (K + 1)) - XTEDGE) / DZ
         XTECAP  = LINEAR (Z (K), XTEDGE, TESLOPE, Z2)

         LESLOPE = (XLE (K + 1) - XLE (K)) / DZ
         XLECAP  = LINEAR (Z (K), XLE (K), LESLOPE, Z2)

         CHORDCAP = XTECAP - XLECAP

         IF (TESLOPE .GT. WAKESLOPE .AND. CHORDCAP .GE. C2)
     >      GO TO 500
      END DO

      STOP 'Tip fudge failed.'


  500 CONTINUE

      ZSTRAIGHT = Z (KS + 1)
      TEANGLE = ATAND (TESLOPE)
      TAPERCAP = CHORDCAP / CROOT

      WRITE (LUNCRT, '(/, A, I3, A, F7.3, A, F6.2, A)')
     >   ' Trailing edge slope adjusted beyond K = ', KS + 1,
     >   ' (Z = ', ZSTRAIGHT, ') to ', TEANGLE, ' degrees'

C     Adjust the outer trailing edge; linearize the leading edge
C     and thickness taper likewise.

      TMSLOPE = (TMAX (KS + 1) - TMAX (KS)) / DZ

      DO K = KS + 2, NZ
         XLE (K)  = LINEAR (Z (KS), XLE (KS),  LESLOPE, Z (K))
         CHORD (K)= LINEAR (Z (KS), XTEDGE,    TESLOPE, Z (K)) - XLE (K)
         TMAX (K) = LINEAR (Z (KS), TMAX (KS), TMSLOPE, Z (K))
      END DO


C     Remainder of T/C, "THICK" factor, and area distributions:

      DO K = K1 + 1, NZ
         TOVERC (K) = TMAX (K) / CHORD (K)
         THICK (K) = TOVERC (K) / TCROOT
         AREA (K) = TMAX (K) * CHORD (K) * AREARATIO
      END DO

      TCCAP = TOVERC (NZ - NCAP) * 100.

C     Left half of wing:

      DO K = KMID + 1, NZ
         KL = KMID - (K - KMID)
         Z (KL) = -Z (K)
         XLE (KL) = XLE (K)
         TMAX (KL) = TMAX (K)
         AREA (KL) = AREA (K)
         CHORD (KL) = CHORD (K)
         THICK (KL) = THICK (K)
         TOVERC (KL) = TOVERC (K)
      END DO

      XLE1 = XLE (1)       ! Needed below if left tip is to be semi-square

C     Wing area and volume:

      CALL LCSQUAD (NZ, Z, CHORD, -ZT, ZT, 'B', PLAN)
      CALL LCSQUAD (NZ, Z, AREA,  -ZT, ZT, 'B', VOL)

C     Print numbers like  "12,345"  rather than  12345.

      IVOL = NINT (VOL)
      IVOL1 = IVOL / 1000
      IVOL2 = IVOL - IVOL1 * 1000
      IPLAN = NINT (PLAN)
      IPLAN1 = IPLAN / 1000
      IPLAN2 = IPLAN - IPLAN1 * 1000
      AREATXT = '   ,   '
      VOLUMETXT = AREATXT
      WRITE (AREATXT (1:3), '(I3)') IPLAN1
      WRITE (AREATXT (5:7), '(I3.3)') IPLAN2
      WRITE (VOLUMETXT (1:3), '(I3)') IVOL1
      WRITE (VOLUMETXT (5:7), '(I3.3)') IVOL2

      WRITE (LUNCRT, '(/, (6X, 2A))')
     >   'Area:         ', AREATXT,
     >   'Volume:       ', VOLUMETXT
      WRITE (LUNCRT, '(6X, A, F8.4)')
     >   'Taper ratio: ', TAPERCAP,
     >   'Tip cap %T/C:', TCCAP

      GO TO 100  ! Repeat until ^Z


  600 CONTINUE

C     Log the results:

      CALL DATE_AND_TIME(DAY)

      WRITE (LUNOUT, '(A, 20X, A)')
     >  'Modified Elliptic OAW Planform Design', DAY
      WRITE (LUNOUT, '(//, A, /)') 'Defining parameters:'
      WRITE (LUNOUT, '(6X, F10.5, 6X, A)')
     >   SWEEP, 'Wing sweep (degrees)',
     >   CROOT, 'Root chord',
     >   ZT,    'Semi-span',
     >   Z1,    'End of cabin minus 70" of shaved chord and thickness',
     >   TCROOT,'Thickness/chord ratio for cabin section',
     >   T1,    'Maximum thickness of cabin section',
     >   TIPRATIO, 'Tip ratio',
     >   Z2,    'Tip-cap station where T/C is monitored',
     >   TAPERRATIO, 'Taper ratio (at tip-cap station)',
     >   S,     'X/C line which is held straight across the span',
     >   A,     'Length of major axis of superellipses except tip cap',
     >   B1,    'Length of minor axis of leading edge superellipses',
     >   B2,    'Length of minor axis of trailing edge superellipses',
     >   POWC1, 'Superellipse exponent for chord (inboard)',
     >   POWC2, 'Superellipse exponent for chord (outboard)',
     >   POWT1, 'Superellipse exponent for thickness (inboard)',
     >   POWT2, 'Superellipse exponent for thickness (outboard)'
      WRITE (LUNOUT, '(/, 2A)')
     >   'Superellipse exponent variation is ', VARYCASE (IVARY)

      IF (KS .LT. NZ - 1) THEN
         WRITE (LUNOUT, '(/, A, I3, A, F7.3, A, F6.2, A)')
     >      'Trailing edge slope adjusted beyond K = ', KS + 1,
     >      ' (Z = ', ZSTRAIGHT, ') to ', TEANGLE, ' degrees'
      END IF

      WRITE (LUNOUT, '(/, A)')
     >   'Planform defining stations as expected by R22OPT:'
      WRITE (LUNOUT, '(/, 2A)')
     >   ' KG          Z       XLE       YLE     CHORD     THICK',
     >   '     TWIST NEWSEC CRANK'
      WRITE (LUNOUT, '(I3, F11.5, 3F10.5, 2F10.6, 2I5)')
     >   (K, Z (K), XLE (K), ZERO, CHORD (K), THICK (K), ZERO, 0, 0,
     >    K = 1, NZ)

      WRITE (LUNOUT, '(/, (6X, 2A))')
     >   'Area:         ', AREATXT,
     >   'Volume:       ', VOLUMETXT
      WRITE (LUNOUT, '(6X, A, F8.4)')
     >   'Taper ratio: ', TAPERCAP,
     >   'Tip cap %T/C:', TCCAP

      IF (KS .EQ. NZ - 1) THEN
         WRITE (LUNOUT, '(/, A)') 'Linear extrapolation to the tip:'
         WRITE (LUNOUT, '(/, I3, F11.5, 3F10.5, 2F10.6, 2I5)')
     >    NZ, Z (NZ), XLE (NZ), ZERO, CHORD (NZ), THICK (NZ), ZERO, 0, 0
      END IF

C     Save results in QPLOTable form:

      WRITE (LUNCTL, '(A)')
     >   'OAW-1 Cabin with Superelliptic/linear Outer Panels'
      WRITE (LUNCTL, '(5(A, F6.4))')
     >   'Superellipse exponents:  ',
     >   POWC1, ' - ', POWC2, ' (chord)  &  ',
     >   POWT1, ' - ', POWT2, ' (thickness)'
      WRITE (LUNCTL, '(A)') 'Span station (ft.)',
     >  'Chord (ft.)  or  Section Area (sq. ft.)', '""',
     >  'Squared-tip case', '""'
      WRITE (LUNCTL, '(2A)')
     >   'Superellipse exponent variation is ', VARYCASE (IVARY), '""'
      WRITE (LUNCTL, '(3(A, F8.4))')
     >   'Taper ratio: ', TAPERCAP, '   Tip ratio: ', TIPRATIO,
     >   '   Tip cap %T/C: ', TCCAP, '""'
      WRITE (LUNCTL, '(6A)')
     >   'Wing area:  ', AREATXT, ' sq. ft.    ',
     >   'Volume:  ', VOLUMETXT, ' cu. ft.'
      WRITE (LUNCTL, '(A)') '""', DAY, ' '
      WRITE (LUNCTL, '(A)')
     >   ' $OPTIONS',
     >   ' xmin=0., xmax=250., ymin=0., ystep=50., ymax=450.,',
     >   ' y2label=''Maximum thickness (ft.)  or  T/C (%)'',',
     >   ' y2scale=0.04, grid=1,',
     >   ' $END'
      WRITE (LUNCTL, '(/, A)')
     >   ' $OPTIONS legend=''Section area'', ycol=5, line=''thi'' $END',
     >   '@square-tip.plt',
     >   ' $OPTIONS legend=''Chord'', ycol=3, line=''solid'' $END',
     >   '@square-tip.plt'
      WRITE (LUNCTL, '(/, 2A)')
     >   ' $OPTIONS legend=''% T/C'', ycol=4, line=''longdash'',',
     >   ' yscale=25. $END',
     >   '@square-tip.plt'
      WRITE (LUNCTL, '(/, 2A)')
     >   ' $OPTIONS legend=''Thickness'', ycol=2, line=''dash'',',
     >   ' yscale=25. $END',
     >   '@square-tip.plt'

      WRITE (LUNSQR, '(A)')
     >'!       Z            T            CHORD        %T/C         AREA'
      WRITE (LUNSQR, '(5(2X, F11.5))')
     >   (Z (K), TMAX (K), CHORD (K), 100. * TOVERC (K), AREA (K),
     >   K = KMID, NZ)

C     Rounded-tip form for codes other than R22OPT.  Use the same NZ.
C     ---------------------------------------------------------------

C     First, set up an elliptic set of control points used to shape
C     the spanwise variation of Z spacings.  (Other standard distributions
C     tend to give too large an increment at the center if they capture
C     the elliptic tip region properly.)

      DZ = ZT / (KMID - 1)
      DO K = KMID, NZ
         ZELLIPSE = ZT - (K - KMID) * DZ        ! Uniform
         ZCTL (K) = ELLIPSE (ZT, ZT, 2.E+0, ZELLIPSE)  ! Bunched toward ZT
         YCTL (K) = ELLIPSE (ZT, CROOT, 2.E+0, ZCTL (K))    ! ~ Chord variation
      END DO

C     Ensure non-zero control pt. at the tip:

      YCSLOPE = (YCTL (NZ - 1) - YCTL (NZ - 2)) /
     >          (ZCTL (NZ - 1) - ZCTL (NZ - 2))
      YCTL (NZ) = LINEAR (ZCTL (NZ - 1), YCTL (NZ - 1), YCSLOPE, ZT)


C     Spanwise Z distribution with dZ varying in the manner of YCTL (*):

      CALL ARBDIS (KMID, ZERO, ZT, KMID, ZCTL (KMID), YCTL (KMID),
     >             'AUTO', LUNCRT, WORK, Z (KMID), IER)
      IF (IER .NE. 0) GO TO 999


C     Include the end of cabin and other significant stations exactly:

      ZFIXED (1) = Z1
      ZFIXED (2) = ZSTRAIGHT
      ZFIXED (3) = Z2

      CALL SMOOTHX (Z (KMID), KMID, ZFIXED, 3, 3, WORK, IER)

      IF (IER .NE. 0) STOP 'SMOOTHX failed.'


C     We now have a smoothly varying distribution of span stations
C     suited to a rounded tip and including 3 defining stations exactly.
C     But we need to determine the subscript ranges defined by these
C     fixed stations:

      CALL INTERVAL (NZ, Z, Z1, ONE, K1)
      KS = K1
      CALL INTERVAL (NZ, Z, ZSTRAIGHT, ONE, KS)
      K2 = NZ
      CALL INTERVAL (NZ, Z, Z2, ONE, K2)

C     Evaluate the chord, etc., in the now-defined ranges:

C     Right cabin panel:

      DO K = KMID, K1
         XLE (K) = ZERO
         TMAX (K) = T1
         AREA (K) = T1 * CROOT * AREARATIO
         CHORD (K) = CROOT
         THICK (K) = ONE
         TOVERC (K) = TCROOT
      END DO

C     Right outer panel, Z1 to ZSTRAIGHT:

      DO K = K1 + 1, KS

         IF (IVARY .EQ. 1) THEN
            POWC = LINEAR  (Z1, POWC1, CHSLOPE, Z (K))
            POWT = LINEAR  (Z1, POWT1, THSLOPE, Z (K))

         ELSE IF (IVARY .EQ. 2) THEN
            CALL QINTRP (1, Z (K), ZT, POWC2, ZERO, Z1, POWC1,
     >                   QSLOPE, POWC)
            CALL QINTRP (1, Z (K), ZT, POWT2, ZERO, Z1, POWT1,
     >                   QSLOPE, POWT)

         ELSE IF (IVARY .EQ. 3) THEN
            CALL QINTRP (1, Z (K), Z1, POWC1, ZERO, ZT, POWC2,
     >                   QSLOPE, POWC)
            CALL QINTRP (1, Z (K), Z1, POWT1, ZERO, ZT, POWT2,
     >                   QSLOPE, POWT)
         END IF

         ZELLIPSE  = Z (K) - Z1
         XLE (K)   = B1 - ELLIPSE (A, B1, POWC, ZELLIPSE)
         CHORD (K) = B1 + ELLIPSE (A, B2, POWC, ZELLIPSE) - XLE (K)
         TMAX (K)  =      ELLIPSE (A, T1, POWT, ZELLIPSE)
      END DO

C     Fudge the tip trailing edge as above, to match R22OPT version as
C     much as possible and avoid too-thin sections (or too large volumes
C     if tip T/C is held by using higher superellipse exponents).
C
C     WARNING:  KS now corresponds to ZSTRAIGHT (where it was KS + 1 above).

      WRITE (LUNCRT, '(/, A, I3, A, F7.3, A, F6.2, A)')
     >   ' Trailing edge slope adjusted beyond K = ', KS,
     >   ' (Z = ', ZSTRAIGHT, ') to ', TEANGLE, ' degrees'

C     Linearize in the range ZSTRAIGHT to Z2:

      XTEDGE = XLE (KS) + CHORD (KS)
      DO K = KS + 1, K2
         XLE (K)  = LINEAR (ZSTRAIGHT, XLE (KS),  LESLOPE, Z (K))
         CHORD (K)= LINEAR (ZSTRAIGHT, XTEDGE, TESLOPE, Z (K)) - XLE (K)
         TMAX (K) = LINEAR (ZSTRAIGHT, TMAX (KS), TMSLOPE, Z (K))
      END DO


C     Round-tip ellipse for the leading edge, Z2 to ZT, tangent at (Z2, XT).
C     The origin is best located at Z2 along the constant X/C straight line.

      LESLOPE = -LESLOPE       ! X axis is pointing forward for this purpose
      ZE = ZT - Z2             ! End coordinate of major axis
      XT = S * CHORD (K2)      ! X coordinate of point of tangency
      ZC = (LESLOPE * ZE * ZE) / (XT + 2.E+0 * LESLOPE * ZE)   ! Center
      AX = ZE - ZC                             ! Major axis
      BX = AX * SQRT (LESLOPE * XT / ZC)       ! Minor axis

      DO K = K2 + 1, NZ - 1
         ZELLIPSE = (Z (K) - Z2) - ZC
         XLE (K) = XC - ELLIPSE (AX, BX, 2.E+0, ZELLIPSE)
      END DO

C     Repeat for the trailing edge:

      XT = CHORD (K2) - XT
      ZC = (TESLOPE * ZE * ZE) / (XT + 2.E+0 * TESLOPE * ZE) ! Center
      AX = ZE - ZC                             ! Major axis
      BX = AX * SQRT (TESLOPE * XT / ZC)       ! Minor axis

      DO K = K2 + 1, NZ - 1
         ZELLIPSE = (Z (K) - Z2) - ZC
         CHORD (K) = XC + ELLIPSE (AX, BX, 2.E+0, ZELLIPSE) - XLE (K)
      END DO

C     ... and for the maximum thickness (parabolic now, not elliptic).

      CALL QINTRP (NZ - K2 - 1, Z (K2 + 1), Z (K2), TMAX (K2), TMSLOPE,
     >             ZT, ZERO, YCSLOPE, TMAX (K2 + 1))

!      XT = TMAX (K2)
!      ZC = (TMSLOPE * ZE * ZE) / (XT + 2.E+0 * TMSLOPE * ZE) ! Center
!      AX = ZE - ZC                             ! Major axis
!      BX = AX * SQRT (TMSLOPE * XT / ZC)       ! Minor axis

!      DO K = K2 + 1, NZ - 1
!         ZELLIPSE = (Z (K) - Z2) - ZC
!         TMAX (K) = ELLIPSE (AX, BX, 2.E+0, ZELLIPSE)
!      END DO


C     Avoid zero chord and thickness at the very tip:

      DZ = Z (NZ - 1) - Z (NZ - 2)
      YCSLOPE = (XLE (NZ - 1) - XLE (NZ - 2)) / DZ
      XLE (NZ) = LINEAR (Z (NZ - 1), XLE (NZ - 1), YCSLOPE, ZT)

      YCSLOPE = (TMAX (NZ - 1) - TMAX (NZ - 2)) / DZ
      TMAX (NZ) = LINEAR (Z (NZ - 1), TMAX (NZ - 1), YCSLOPE, ZT)

      YCSLOPE = (CHORD (NZ - 1) - CHORD (NZ - 2)) / DZ
      CHORD (NZ) = LINEAR (Z (NZ - 1), CHORD (NZ - 1), YCSLOPE, ZT)


C     Remainder of T/C, "THICK" factor, and area distributions:

      DO K = K1 + 1, NZ
         TOVERC (K) = TMAX (K) / CHORD (K)
         THICK (K) = TOVERC (K) / TCROOT
         AREA (K) = TMAX (K) * CHORD (K) * AREARATIO
      END DO

      TCCAP = TOVERC (K2) * 100.

C     Left half of wing:

      DO K = KMID + 1, NZ
         KL = KMID - (K - KMID)
         Z (KL) = -Z (K)
         XLE (KL) = XLE (K)
         TMAX (KL) = TMAX (K)
         AREA (KL) = AREA (K)
         CHORD (KL) = CHORD (K)
         THICK (KL) = THICK (K)
         TOVERC (KL) = TOVERC (K)
      END DO


C     Semi-elliptic left tip instead of fully-rounded?

      IF (SEMI_SQUARE_LEFT_TIP) THEN

C        The origin is best located at (Z = Z-tip-cap, X = XLE1):

         TESLOPE = -TESLOPE       ! Left wing now
         K2 = NZ - K2 + 1         ! Left wing tip-cap station
         XT = XLE (K2) + CHORD (K2) - XLE1   ! Where tangent to ellipse is known
         ZC = (TESLOPE * ZE * ZE) / (XT + 2.E+0 * TESLOPE * ZE)   ! Center
         AX = ZE - ZC                             ! Major axis
         BX = AX * SQRT (TESLOPE * XT / ZC)       ! Minor axis

         DO K = 2, K2 - 1
            XLE (K) = LINEAR (Z (1), XLE1, LESLOPE, Z (K))  ! Recovers originals
            ZELLIPSE = -(Z (K) - Z (K2)) - ZC
            CHORD (K) = XLE1 + ELLIPSE (AX, BX, 2.E+0, ZELLIPSE) -
     >                  XLE (K)
         END DO

C        Avoid zero chord at the very tip.  (Thickness has not changed.)

         DZ = Z (3) - Z (2)
         YCSLOPE = (XLE (3) - XLE (2)) / DZ
         XLE (1) = LINEAR (Z (2), XLE (2), YCSLOPE, Z (1))

         YCSLOPE = (CHORD (3) - CHORD (2)) / DZ
         CHORD (1) = LINEAR (Z (2), CHORD (2), YCSLOPE, Z (1))

         DO K = 1, K2 - 1
            TOVERC (K) = TMAX (K) / CHORD (K)
            THICK (K) = TOVERC (K) / TCROOT
            AREA (K) = TMAX (K) * CHORD (K) * AREARATIO
         END DO

      END IF


C     Wing area and volume:

      CALL LCSQUAD (NZ, Z, CHORD, -ZT, ZT, 'B', PLAN)
      CALL LCSQUAD (NZ, Z, AREA,  -ZT, ZT, 'B', VOL)

      IVOL = NINT (VOL)
      IVOL1 = IVOL / 1000
      IVOL2 = IVOL - IVOL1 * 1000
      IPLAN = NINT (PLAN)
      IPLAN1 = IPLAN / 1000
      IPLAN2 = IPLAN - IPLAN1 * 1000
      WRITE (AREATXT (1:3), '(I3)') IPLAN1
      WRITE (AREATXT (5:7), '(I3.3)') IPLAN2
      WRITE (VOLUMETXT (1:3), '(I3)') IVOL1
      WRITE (VOLUMETXT (5:7), '(I3.3)') IVOL2

      WRITE (LUNCRT, '(/, (6X, 2A))')
     >   'Round-tip planform area:  ', AREATXT,
     >   'Volume:                   ', VOLUMETXT
      WRITE (LUNCRT, '(6X, A, F8.4)')
     >   'Tip cap %T/C:            ', TCCAP


C     Planform definition for R22OPT (still used to insert wing sections):

      WRITE (LUNOUT, '(/, A)')
     >  'Round-tip planform defining stations in R22OPT format:'
      WRITE (LUNOUT, '(/, 2A)')
     >   ' KG          Z       XLE       YLE     CHORD     THICK',
     >   '     TWIST NEWSEC CRANK'
      WRITE (LUNOUT, '(I3, F11.5, 3F10.5, 2F10.6, 2I5)')
     >   (K, Z (K), XLE (K), ZERO, CHORD (K), THICK (K), ZERO, 0, 0,
     >    K = 1, NZ)

      WRITE (LUNOUT, '(/, (6X, 2A))')
     >   'Area:         ', AREATXT,
     >   'Volume:       ', VOLUMETXT
      WRITE (LUNOUT, '(6X, A, F8.4)')
     >   'Tip cap %T/C:', TCCAP


C     Save round-tip planform results in QPLOTable form:

      WRITE (LUNCTL, '(A)') 'END FRAME',
     >   'OAW-1 Cabin with Superelliptic/linear Outer Panels'
      WRITE (LUNCTL, '(5(A, F6.4))')
     >   'Superellipse exponents:  ',
     >   POWC1, ' - ', POWC2, ' (chord)  &  ',
     >   POWT1, ' - ', POWT2, ' (thickness)'
      WRITE (LUNCTL, '(A)') 'Span station (ft.)',
     >   'Chord (ft.)  or  Section Area (sq. ft.)', '""',
     >   'Round-tip case', '""'
      WRITE (LUNCTL, '(2A)')
     >   'Superellipse exponent variation is ', VARYCASE (IVARY), '""'
      WRITE (LUNCTL, '(3(A, F9.4))')
     >   'Taper ratio: ', TAPERCAP, '   Tip ratio: ', TIPRATIO,
     >   '   Tip cap %T/C: ', TCCAP, '""'
      WRITE (LUNCTL, '(6A)')
     >   'Wing area:  ', AREATXT, ' sq. ft.    ',
     >   'Volume:  ', VOLUMETXT, ' cu. ft.'
      WRITE (LUNCTL, '(A)') '""', DAY, ' '
      WRITE (LUNCTL, '(A)')
     >   ' $OPTIONS',
     >   ' xmin=0., xmax=250., ymin=0., ystep=50., ymax=450.,',
     >   ' y2label=''Maximum thickness (ft.)  or  T/C (%)'',',
     >   ' y2scale=0.04, grid=1,',
     >   ' $END'
      WRITE (LUNCTL, '(/, A)')
     >   ' $OPTIONS legend=''Section area'', ycol=5, line=''thi'' $END',
     >   '@round-tip.plt',
     >   ' $OPTIONS legend=''Chord'', ycol=3, line=''solid'' $END',
     >   '@round-tip.plt'
      WRITE (LUNCTL, '(/, 2A)')
     >   ' $OPTIONS legend=''% T/C'', ycol=4, line=''longdash'',',
     >   ' yscale=25. $END',
     >   '@round-tip.plt'
      WRITE (LUNCTL, '(/, 2A)')
     >   ' $OPTIONS legend=''Thickness'', ycol=2, line=''dash'',',
     >   ' yscale=25. $END',
     >   '@round-tip.plt'

      WRITE (LUNRND, '(A)')
     >'!       Z            T            CHORD        %T/C         AREA'
      WRITE (LUNRND, '(5(2X, F11.5))')
     >   (Z (K), TMAX (K), CHORD (K), 100. * TOVERC (K), AREA (K),
     >   K = KMID, NZ)



C     Termination:
C     ------------

      WRITE (LUNCRT, '(/, 1X, A)')
     >  'See "superellipse.out" and "superellipse.ctl" for the results.'

  999 STOP ' '
      END
