C***********************************************************************
C
      PROGRAM CGRID
C
C     CGRID sets up a C-grid boundary for a specified airfoil,
C     initializes the interior in one of two ways, and smooths it.
C
C     06/12/95  DAS  Adaptation of SYN87's OUTER routine
C                    for testing TTM2D orthogonality mods.
C     09/11/95   "   Variant to drive ELLIP2D.
C     09/17/95   "   ELLIP2D is now fully argument-driven.
C     08/08/96   "   ELLIP2D argument count reduced.
C     06/30/97   "   Swap I1 & I2 controls for upper half.
C     07/07/97   "   Went to 2-sided distributions at I=1/ILE/IL,
C                    with partially-geometric option for I=ILE at
C                    Steve Edwards's suggestion.
C     07/28/97   "   Read arbitrary geometry rather than a gridded
C                    airfoil.
C     08/02/97   "   4-quadrant option.
C     08/05/97 DS/SE Partially-geometric option off trailing edge too.
C     05/07/98  DAS  Fortran 90 upgrade.  ELLIP2D is fully portable now.
C     05/15/98   "   Save 32-bit form of grid (only).
C     10/24/07   "   Removed the 32-bit form.
C     10/25/07   "   The "1" was missing from the output PLOT3D grid,
C                    (which is now 'cgrid.xy', not 'cgrid.xyz') for
C                    compatibility with COMPARE_PATCHES.
C
C     David Saunders, Sterling Software/NASA Ames, Mt. View, CA
C
C***********************************************************************

      IMPLICIT NONE

C     Constants:

      INTEGER
     >   IDM, JDM, LUNCRT, LUNFOIL, LUNKBD, LUNNML, LUNPLOT, MXGEOM
      REAL
     >   C1, C2, C3, C4, HALF, ONE, ZERO
      LOGICAL
     >   NEW
      CHARACTER
     >   LOOSE * 1, TIGHT * 1

      PARAMETER
     >  (IDM     = 297,
     >   JDM     = 97,
     >   LUNCRT  = 6,
     >   LUNFOIL = 1,
     >   LUNKBD  = 5,
     >   LUNNML  = 2,
     >   LUNPLOT = 3,
     >   MXGEOM  = 265,
     >   C1      = 0.04,
     >   C2      = 0.,
     >   C3      = 0.3,
     >   C4      = 0.66,
     >   HALF    = 0.5,
     >   ONE     = 1.,
     >   ZERO    = 0.,
     >   NEW     = .TRUE.,  ! New data for all calls to LCSFIT
     >   LOOSE   = 'B',     ! Bessel = "loose" fit
     >   TIGHT   = 'M')     ! Monotonic fit for X at the leading edge

C     (Most of the) arguments used by ELLIP2D:

      INTEGER
     >   IT2MAX, IT2MAX4, ITFLOAT, ITFREEZE, JLAYER
      LOGICAL
     >   PRINT2D
      REAL
     >   CONV2D, CONVMIN, DMAX2D, OMG2D, 
     >   POWERI, POWERJ, EXPI1, EXPI2, EXPJ1, EXPJ2,
     >   POWERI4, POWERJ4, EXPI14, EXPI24, EXPJ14, EXPJ24, OMG2D4,
     >   FGLIMIT, FGLIMI1, FGLIMI2, FGLIMJ1, FGLIMJ2,
     >   URFG, URI1, URI2, URJ1, URJ2, URFLOAT
      CHARACTER
     >   BGMODE*4,  FGMODE*4, SPMODE*4, BGMODE4*4, FGMODE4*4,
     >   BGSWAP*4,  FGSWAP*4

C     Local variables:

      INTEGER
     >   I, IEDGE, IER, IL, ILE, ITL, ITU, J, JL, KEY,
     >   NBLAYER, NBLOCKS, NGEOM, NFULL, NHALF, NL, NU
      REAL
     >   ARCOUTER(IDM), ARCS(IDM,JDM,2), ARCTOT, BETA, CHORD, CS,
     >   D1, D1WING, D2, D2JL, D2WING, D3WING, DA, DM,
     >   DU, DX1, OSPNOSE, OSPTAIL, PHIB, PI, QCIRC, R, RADIUS,
     >   RBLAYER, S, SA(5), SEVAL(IDM), SGEOM(MXGEOM), SLOWER,
     >   SNORM(IDM), SOUT(IDM), STOTAL, SUPPER, TFIWORK(2*(IDM + JDM)),
     >   XLE, XLPRC, XMAX, XMIN, XTE,
     >   YBPRC, YLE, YTE, YTMRC, YMAX, YMIN,
     >   X(IDM,JDM), XGEOM(MXGEOM), XOUT(IDM), XILE(JDM),
     >   Y(IDM,JDM), YGEOM(MXGEOM), YLBACK(JDM), YUBACK(JDM), YOUT(IDM)
      LOGICAL
     >   EROUND, ESHARP, NROUND, ROUNDED
      CHARACTER
     >   ANSWER*1

C     Control:

      NAMELIST /INPUTS/
     >   IL, ITL, ITU, JL, NBLOCKS, NBLAYER, RBLAYER,
     >   CHORD, D1WING, D2JL, D2WING, D3WING, OSPNOSE, OSPTAIL, RADIUS,
     >   XLE, XMAX, XMIN, YLE, YMAX, YMIN,
     >   IT2MAX, ITFLOAT, ITFREEZE, JLAYER, PRINT2D,
     >   CONV2D, CONVMIN, DMAX2D, OMG2D, 
     >   BGMODE, FGMODE, SPMODE, BGMODE4, FGMODE4, OMG2D4, IT2MAX4,
     >   POWERI, POWERJ, EXPI1, EXPI2, EXPJ1, EXPJ2,
     >   POWERI4, POWERJ4, EXPI14, EXPI24, EXPJ14, EXPJ24,
     >   FGLIMIT, URFG, URFLOAT

C****>   FGLIMI1, FGLIMI2, FGLIMJ1, FGLIMJ2,
C****>   URI1, URI2, URJ1, URJ2, URFLOAT,

C     Execution:


      PI = 2.0 * ASIN (ONE)

      WRITE (LUNCRT, '(A)', ADVANCE='NO')
     >   ' Euler round/sharp or N-S round/sharp? (a, b, c, or d): '

      READ  (LUNKBD, '(A)') ANSWER
      EROUND  = ANSWER == 'a' .OR. ANSWER == 'A'
      ESHARP  = ANSWER == 'b' .OR. ANSWER == 'B'
      NROUND  = ANSWER == 'c' .OR. ANSWER == 'C'
      ROUNDED = EROUND .OR. NROUND

      IF (EROUND) THEN
         OPEN (UNIT=LUNNML,  FILE='cgrid.inp-euler-round', STATUS='OLD')
         OPEN (UNIT=LUNFOIL, FILE='cgrid.geo-round',       STATUS='OLD')
      ELSE IF (ESHARP) THEN
         OPEN (UNIT=LUNNML,  FILE='cgrid.inp-euler-sharp', STATUS='OLD')
         OPEN (UNIT=LUNFOIL, FILE='cgrid.geo-sharp',       STATUS='OLD')
      ELSE  IF (NROUND) THEN
         OPEN (UNIT=LUNNML,  FILE='cgrid.inp-ns-round', STATUS='OLD')
         OPEN (UNIT=LUNFOIL, FILE='cgrid.geo-round',    STATUS='OLD')
      ELSE
         OPEN (UNIT=LUNNML,  FILE='cgrid.inp-ns-sharp', STATUS='OLD')
         OPEN (UNIT=LUNFOIL, FILE='cgrid.geo-sharp',    STATUS='OLD')
      END IF

      OPEN (UNIT=LUNPLOT, FILE='cgrid.xy', STATUS='UNKNOWN',
     >      FORM='UNFORMATTED')

      READ (LUNNML, NML=INPUTS)

C     Read the airfoil geometry (normalized, standard PROFILE format):

      READ (LUNFOIL,*) ! Title
      READ (LUNFOIL,*) NU
      READ (LUNFOIL,*) (XGEOM(I), YGEOM(I), I = 1, NU)
      READ (LUNFOIL,*) NL
      NGEOM = NU + NL - 1

      IF (NGEOM > MXGEOM) THEN
         WRITE (LUNCRT, '(A)') ' Too many airfoil coordinates.'
         GO TO 999
      END IF

      XGEOM(NL:NGEOM) = XGEOM(1:NU)
      YGEOM(NL:NGEOM) = YGEOM(1:NU)

      READ (LUNFOIL,*) (XGEOM(I), YGEOM(I), I = NL, 1, -1)

      CLOSE (UNIT=LUNFOIL)

      DO I = 1, NGEOM
         XGEOM(I) = XGEOM(I) * CHORD + XLE
         YGEOM(I) = YGEOM(I) * CHORD + YLE
      END DO

C     Impose the computational grid on the airfoil:

      ILE = (ITL + ITU) / 2
      NHALF = ILE - ITL + 1
      NFULL = ITU - ITL + 1

      CALL FOILGRD (NHALF, ZERO, ONE, C1, C2, C3, C4, SNORM)

      IF (ROUNDED) THEN

C        Arc lengths along both surfaces:

         CALL CHORDS2D (NGEOM, XGEOM, YGEOM, .FALSE., STOTAL, SGEOM)

         SLOWER = SGEOM(NL)
         SUPPER = STOTAL - SLOWER

         DO I = 1, NHALF
            SEVAL(I) = SLOWER * (ONE - SNORM(NHALF + 1 - I))
            SEVAL(I + NHALF - 1) = SLOWER + SUPPER * SNORM(I)
         END DO

C        Interpolation along the arc for both surfaces:

         CALL LCSFIT (NGEOM, SGEOM, XGEOM, NEW, TIGHT, NFULL, SEVAL,
     >                X(ITL,1), X(ITL,1))

         CALL LCSFIT (NGEOM, SGEOM, YGEOM, NEW, LOOSE, NFULL, SEVAL,
     >                Y(ITL,1), Y(ITL,1))

      ELSE  ! Sharp leading edge case:

C        Arc lengths along the upper surface:

         CALL CHORDS2D (NU, XGEOM(NL), YGEOM(NL), .FALSE., SUPPER,
     >                  SGEOM)

         DO I = 1 , NHALF
            SEVAL(I) = SUPPER * SNORM(I)
         END DO

C        Upper surface interpolation along the arc.

         CALL LCSFIT (NU, SGEOM, XGEOM(NL), NEW, TIGHT, NHALF, SEVAL,
     >                X(ILE,1), X(ILE,1))

         CALL LCSFIT (NU, SGEOM, YGEOM(NL), NEW, LOOSE, NHALF, SEVAL,
     >                Y(ILE,1), Y(ILE,1))

C        Repeat for the lower surface.

         CALL CHORDS2D (NL, XGEOM, YGEOM, .FALSE., SLOWER, SGEOM)

         DO I = 1, NHALF
            SEVAL(I) = SLOWER * (ONE - SNORM(NHALF + 1 - I))
         END DO

         CALL LCSFIT (NL, SGEOM, XGEOM, NEW, TIGHT, NHALF, SEVAL,
     >                X(ITL,1), X(ITL,1))

         CALL LCSFIT (NL, SGEOM, YGEOM, NEW, LOOSE, NHALF, SEVAL,
     >                Y(ITL,1), Y(ITL,1))

      END IF

      XTE = (X(ITL,1) + X(ITU,1)) * HALF
      YTE = (Y(ITL,1) + Y(ITU,1)) * HALF

C     C-grid boundary beyond airfoil: leave the trailing edge smoothly.

      XOUT(1) = (X(ITL+2,1) + X(ITU-2,1)) * HALF
      YOUT(1) = (Y(ITL+2,1) + Y(ITU-2,1)) * HALF
      XOUT(2) = (X(ITL+1,1) + X(ITU-1,1)) * HALF
      YOUT(2) = (Y(ITL+1,1) + Y(ITU-1,1)) * HALF
      XOUT(3) = XTE
      YOUT(3) = YTE
      DX1     = XTE - XOUT(1)
      XOUT(4) = XTE + DX1
      YOUT(4) = YOUT(3) * 2. - YOUT(1)
      XOUT(5) = XMAX - DX1
      YOUT(5) = ZERO
      XOUT(6) = XMAX
      YOUT(6) = ZERO

      DX1 = XOUT(2) - XOUT(1)

      CALL EXPDIS4 (3, XTE, XMAX, DX1, IL - ITU + 1, X(ITU,1),
     >              BETA, -LUNCRT)

      CALL LCSFIT (6, XOUT, YOUT, .TRUE., 'B', IL - ITU + 1, X(ITU,1),
     >             Y(ITU,1), Y(ITU,1))

      DO I = ITU, IL
         J = IL + 1 - I
         Y(J,1) = Y(I,1)
         X(J,1) = X(I,1)
      END DO


C     Generate the outer boundary shape and distribute points along it.

      QCIRC  = HALF * PI * RADIUS
      STOTAL = 2.* (XMAX - XMIN + QCIRC) + (YMAX - YMIN) - 4.* RADIUS
      DA = STOTAL / REAL (IL - 1)
      XLPRC = XMIN + RADIUS
      YTMRC = YMAX - RADIUS
      YBPRC = YMIN + RADIUS

C     Original distribution: uniform in arc length.
C     Set up arc lengths at ends of each segment of the outer boundary,
C     with SA(0) = 0. understood for bottom right corner.

      SA(1) = XMAX  - XLPRC
      SA(2) = SA(1) + QCIRC
      SA(3) = SA(2) + YTMRC - YBPRC
      SA(4) = SA(3) + QCIRC
      SA(5) = STOTAL

      XOUT(1) = XMAX
      YOUT(1) = YMIN
      KEY = 1

      DO I = 2, IL - 1

         S = DA * REAL (I - 1)

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

C     Kludge to ensure true symmetry if needed for test purposes:

      IF (ABS (YMAX + YMIN) < YMAX * 0.00001) THEN
         YOUT(ILE) = ZERO
         J = IL + 1
         DO I = 1, ILE
            J = J - 1
            XOUT(J) = XOUT(I)
            YOUT(J) = -YOUT(I)
         END DO
      END IF

C     Redistribute the uniform symmetry-plane boundary?

      IF (OSPNOSE /= ONE .OR. OSPTAIL /= ONE) THEN

         CALL CHORDS2D (IL, XOUT, YOUT, .FALSE., ARCTOT, SOUT)

         DU = ARCTOT / REAL (IL - 1)
         D1 = DU * OSPTAIL     ! "Outer spacing, tail boundary"
         D2 = DU * OSPNOSE     ! "Outer spacing, nose boundary"
         DM = D1 * 0.05        ! Heuristic choice for interior increment

C        Set up an artificial boundary point opposite the root trailing edge:

         IEDGE = ITL
         S = XMAX - X(ITL,1)
         CALL INTERVAL (ILE/2, SOUT, S, ONE, IEDGE)

C        Vinokur distributions along the lower half:

         CALL HTDIS4 (.TRUE., ZERO, SOUT(IEDGE), D1, DM, ITL,
     >                ARCOUTER, -LUNCRT, IER)

         IF (IER /= 0 .AND. IER /= 3) THEN
            WRITE (LUNCRT,*) 'HTDIS4 trouble (1a). IER: ', IER
            WRITE (LUNCRT,*) 'RANGE, D1, DM: ', SOUT(IEDGE), D1, DM
            GO TO 999
         END IF

         S = ARCOUTER(ITL-1)
         CALL HTDIS4 (.TRUE., S, SOUT(ILE), DM, D2, ILE - ITL + 2,
     >                ARCOUTER(ITL-1), -LUNCRT, IER)

         IF (IER /= 0 .AND. IER /= 3) THEN
            WRITE (LUNCRT,*) 'HTDIS4 trouble (1b). IER: ', IER
            WRITE (LUNCRT,*) 'A, B, DM, D2: ', S, SOUT(ILE), DM, D2
            GO TO 999
         END IF

C        Can't redistribute in place, so ...

         TFIWORK(1:IL)       = XOUT(1:IL)
         TFIWORK(IL+1:IL+IL) = YOUT(1:IL)

         CALL LCSFIT (ILE, SOUT, TFIWORK(1),    .TRUE., 'M', ILE,
     >                ARCOUTER, XOUT, XOUT)

         CALL LCSFIT (ILE, SOUT, TFIWORK(IL+1), .TRUE., 'M', ILE,
     >                ARCOUTER, YOUT, YOUT)

C        Likewise for the upper half:

         IEDGE = IL - IEDGE + 1

         CALL HTDIS4 (.TRUE., SOUT(IEDGE), SOUT(IL), DM, D1, ITL,
     >                ARCOUTER(ITU), -LUNCRT, IER)

         IF (IER /= 0 .AND. IER /= 3) THEN
            WRITE (LUNCRT,*) 'HTDIS4 trouble (2a). IER: ', IER
            WRITE (LUNCRT,*) 'S1,S2,DM,D2: ', SOUT(IEDGE),SOUT(IL),DM,D1
            GO TO 999
         END IF

         S = ARCOUTER(ITU+1)
         CALL HTDIS4 (.TRUE., SOUT(ILE), S, D2, DM, ITU - ILE + 2,
     >                ARCOUTER(ILE), -LUNCRT, IER)

         IF (IER /= 0 .AND. IER /= 3) THEN
            WRITE (LUNCRT,*) 'HTDIS4 trouble (2b). IER: ', IER
            WRITE (LUNCRT,*) 'S1,S2,D2,DM: ', SOUT(ILE),S,D2,DM
            GO TO 999
         END IF

         TFIWORK(1:IL)       = XOUT(1:IL)
         TFIWORK(IL+1:IL+IL) = YOUT(1:IL)

         CALL LCSFIT (ILE, SOUT(ILE), TFIWORK(ILE),    .TRUE., 'M', ILE,
     >                ARCOUTER(ILE), XOUT(ILE), XOUT(ILE))

         CALL LCSFIT (ILE, SOUT(ILE), TFIWORK(IL+ILE), .TRUE., 'M', ILE,
     >                ARCOUTER(ILE), YOUT(ILE), YOUT(ILE))
      END IF

      X(1:IL,JL) = XOUT(1:IL)
      Y(1:IL,JL) = YOUT(1:IL)

C     Downstream boundary:

      IF (NBLOCKS < 4) THEN

C        The 1- or 2-block scheme cannot control radial spacing at the
C        trailing edge directly.  Adjust the D2 input to reflect control
C        at the far boundary:

         R = (XMAX - X(ILE,1)) / CHORD
         D2WING = D1WING + R * (D2WING - D1WING)
         D1 = CHORD * D2WING   ! Careful - not D1WING here
      ELSE
         D1 = CHORD * D3WING
      END IF

      D2 = CHORD * D2JL
      YUBACK(1)  = ZERO
      YUBACK(JL) = YOUT(IL)

      CALL VINOKUR (1, JL, D1, D2, YUBACK, LUNCRT, IER)

      YLBACK(1)  = YUBACK(1)
      YLBACK(JL) = YOUT(1)

      CALL VINOKUR (1, JL, D1, D2, YLBACK, LUNCRT, IER)

      X(1,1:JL)  = XMAX
      Y(1,1:JL)  = YLBACK(1:JL)
      X(IL,1:JL) = XMAX
      Y(IL,1:JL) = YUBACK(1:JL)

C     Artificial boundary forward of leading edge?
 
      IF (NBLOCKS > 1) THEN

         XILE(1) = X(ILE,1)
         XILE(JL) = XOUT(ILE)
         D1 = CHORD * D1WING

C        NBLAYER <= 2 suppresses adjustment of basic Vinokur distribution:

         CALL BLGRID (JL, D1, D2, NBLAYER, RBLAYER, XILE, LUNCRT, IER)

         R = (YOUT(ILE) - Y(ILE,1)) / (XOUT(ILE) - XILE(1))
         DO J = 1, JL
            X(ILE,J) = XILE(J)
            Y(ILE,J) = Y(ILE,1) + (XILE(J) - XILE(1)) * R
         END DO

      END IF

C     Artificial boundaries off the trailing edge?

      IF (NBLOCKS == 4) THEN

         D1 = CHORD * D2WING

         CALL TEBOUND (IDM, JDM, IL, JL, ILE, ITL, ITU,  1, X, Y, D1,
     >                 NBLAYER, RBLAYER, TFIWORK, LUNCRT)

         CALL TEBOUND (IDM, JDM, IL, JL, ILE, ITL, ITU, IL, X, Y, D1,
     >                 NBLAYER, RBLAYER, TFIWORK, LUNCRT)
      END IF

C     Initialize and smooth the interior:

      IF (NBLOCKS == 1) THEN

         CALL TFI2D (IDM, 1,   IL, 1, JL, X, Y, TFIWORK)

         CALL ELLIP2D (IDM, JDM, 1,   IL, 1, JL, X, Y, ARCS, PRINT2D,
     >                 BGMODE, FGMODE, SPMODE,
     >                 IT2MAX, ITFLOAT, ITFREEZE, JLAYER,
     >                 CONV2D, CONVMIN, DMAX2D, OMG2D,
     >                 POWERI, POWERJ, EXPI1, EXPI2, EXPJ1, EXPJ2, 
     >                 FGLIMIT, URFG, URFLOAT)

      ELSE IF (NBLOCKS == 2) THEN

         CALL TFI2D (IDM, 1,  ILE, 1, JL, X, Y, TFIWORK)

         CALL ELLIP2D (IDM, JDM, 1, ILE,  1, JL, X, Y, ARCS, PRINT2D,
     >                 BGMODE, FGMODE, SPMODE,
     >                 IT2MAX, ITFLOAT, ITFREEZE, JLAYER,
     >                 CONV2D, CONVMIN, DMAX2D, OMG2D,
     >                 POWERI, POWERJ, EXPI1, EXPI2, EXPJ1, EXPJ2, 
     >                 FGLIMIT, URFG, URFLOAT)

         CALL TFI2D (IDM, ILE, IL, 1, JL, X, Y, TFIWORK)

         BGSWAP      = BGMODE
         BGSWAP(1:1) = BGMODE(2:2)
         BGSWAP(2:2) = BGMODE(1:1)

         FGSWAP      = FGMODE
         FGSWAP(1:1) = FGMODE(2:2)
         FGSWAP(2:2) = FGMODE(1:1)

         CALL ELLIP2D (IDM, JDM, ILE, IL, 1, JL, X, Y, ARCS, PRINT2D,
     >                 BGSWAP, FGSWAP, SPMODE,
     >                 IT2MAX, ITFLOAT, ITFREEZE, JLAYER,
     >                 CONV2D, CONVMIN, DMAX2D, OMG2D,
     >                 -POWERI, POWERJ, EXPI2, EXPI1, EXPJ1, EXPJ2, 
     >                 FGLIMIT, URFG, URFLOAT)

      ELSE  ! 4 quadrants

         CALL TFI2D (IDM, 1,  ITL, 1, JL, X, Y, TFIWORK)

         CALL ELLIP2D (IDM, JDM, 1, ITL,  1, JL, X, Y, ARCS, PRINT2D,
     >                 BGMODE4, FGMODE4, SPMODE,
     >                 IT2MAX4, ITFLOAT, ITFREEZE, JLAYER,
     >                 CONV2D, CONVMIN, DMAX2D, OMG2D4,
     >                 POWERI4, POWERJ4, EXPI14, EXPI24, EXPJ14, EXPJ24, 
     >                 FGLIMIT, URFG, URFLOAT)


         CALL TFI2D (IDM, ITL,ILE, 1, JL, X, Y, TFIWORK)

         CALL ELLIP2D (IDM, JDM, ITL, ILE, 1, JL, X, Y, ARCS, PRINT2D,
     >                 BGMODE, FGMODE, SPMODE,
     >                 IT2MAX, ITFLOAT, ITFREEZE, JLAYER,
     >                 CONV2D, CONVMIN, DMAX2D, OMG2D,
     >                 POWERI, POWERJ, EXPI1, EXPI2, EXPJ1, EXPJ2, 
     >                 FGLIMIT, URFG, URFLOAT)

         CALL TFI2D (IDM, ILE,ITU, 1, JL, X, Y, TFIWORK)

         BGSWAP      = BGMODE
         BGSWAP(1:1) = BGMODE(2:2)
         BGSWAP(2:2) = BGMODE(1:1)

         FGSWAP      = FGMODE
         FGSWAP(1:1) = FGMODE(2:2)
         FGSWAP(2:2) = FGMODE(1:1)

         CALL ELLIP2D (IDM, JDM, ILE, ITU, 1, JL, X, Y, ARCS, PRINT2D,
     >                 BGSWAP, FGSWAP, SPMODE,
     >                 IT2MAX, ITFLOAT, ITFREEZE, JLAYER,
     >                 CONV2D, CONVMIN, DMAX2D, OMG2D,
     >                 -POWERI, POWERJ, EXPI2, EXPI1, EXPJ1, EXPJ2, 
     >                 FGLIMIT, URFG, URFLOAT)

         CALL TFI2D (IDM, ITU, IL, 1, JL, X, Y, TFIWORK)

         BGSWAP      = BGMODE4
         BGSWAP(1:1) = BGMODE4(2:2)
         BGSWAP(2:2) = BGMODE4(1:1)

         FGSWAP      = FGMODE4
         FGSWAP(1:1) = FGMODE4(2:2)
         FGSWAP(2:2) = FGMODE4(1:1)

         CALL ELLIP2D (IDM, JDM, ITU, IL, 1, JL, X, Y, ARCS, PRINT2D,
     >                 BGSWAP, FGSWAP, SPMODE,
     >                 IT2MAX4, ITFLOAT, ITFREEZE, JLAYER,
     >                 CONV2D, CONVMIN, DMAX2D, OMG2D4,
     >                 -POWERI4,POWERJ4, EXPI24, EXPI14, EXPJ14, EXPJ24, 
     >                 FGLIMIT, URFG, URFLOAT)
      END IF

      WRITE (LUNPLOT) 1
      WRITE (LUNPLOT) IL, JL
      WRITE (LUNPLOT) ((X(I,J), I = 1, IL), J = 1, JL),
     >                ((Y(I,J), I = 1, IL), J = 1, JL)

  999 WRITE (LUNCRT, '(A)')

      END PROGRAM CGRID

C*******************************************************************************
C
      SUBROUTINE BLGRID (N, D1, D2, NBLAYER, RBLAYER, X, LUNERR, IER)
C
C        BLGRID generates a specialized 1-D grid intended for radial lines
C     starting in a boundary layer and ending in the far field.  The first
C     and last X values are assumed to be in X(1) and X(N).  Points 1
C     through NBLAYER have spacing that varies geometrically according to
C     the ratio RBLAYER (or uniformly if RBLAYER is 1.0).  Spacing is
C     Vinokur-type beyond that.  X is probably arc-length (O. to STOTAL),
C     but could be decreasing if desired.  D1 & D2 are both positive.
C
C     Use NBLAYER <= 2 to suppress the special treatment.
C     See VINOKUR and HTDIS4 for further details, including the LUNERR
C     and IER arguments.
C
C     07/08/97  DAS  Effort to emulate the single-piece C grid result for
C                    the two-halves case forward of the leading edge.
C     08/06/97   "   Generalized to allow application off the trailing edge.
C
C*******************************************************************************

      IMPLICIT NONE

C     Arguments.

      INTEGER  N, NBLAYER, LUNERR, IER
      REAL     D1, D2, RBLAYER, X(N)

C     Local variables.

      INTEGER  I, NF
      REAL     DX

C     Execution.

      DX = SIGN (D1, X(N) - X(1))
      NF = MAX (NBLAYER, 2)  ! Force a single pass if special treatment
                             ! is suppressed
      DO I = 2, NF
         X(I) = X(I-1) + DX
         DX = DX * RBLAYER
      END DO

      DX = ABS (DX / RBLAYER)

C     Overlap the two distributions by one cell:

      CALL VINOKUR (NF - 1, N, DX, D2, X, LUNERR, IER)

      RETURN
      END


C*******************************************************************************
C
      SUBROUTINE TEBOUND (IDM, JDM, IL, JL, ILE, ITL, ITU, IEND,
     >                    X, Y, D1TRAIL, NBLAYER, RBLAYER, T, LUNERR)
C
C     TEBOUND sets up an artificial grid boundary line from the trailing edge
C     to the outer C boundary as part of using a 4-subblock CH scheme in order
C     to provide radial spacing control at the trailing edge that is independent
C     of the spacings at the leading edge and far downstream boundary.
C
C     A 6-point spline is established orthogonal to the trailing edge and to
C     the outer boundary, and this is evaluated at points forming a hybrid
C     geometric + Vinokur distribution along the arc.
C
C     TEBOUND serves for either above or below the wing (if awkwardly).
C
C     05/09/96  D.Saunders  Initial implementation for CH_GRID.
C     08/02/97      "       Belated adaptation for 2-D application (CGRID).
C
C*******************************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER    IDM, JDM   ! Grid limits
      INTEGER    IL, JL     ! Active limits
      INTEGER    ILE, ITL, ITU  ! The usual ...
      INTEGER    IEND       ! 1 for below the wing, IL for above it
      REAL       X(IDM,JDM),
     >           Y(IDM,JDM)
      REAL       D1TRAIL    ! Desired radial increment at the trailing edge
      INTEGER    NBLAYER    ! Number of points in the boundary layer;
                            ! NBLAYER <=2 gives pure Vinokur throughout
      REAL       RBLAYER    ! Geometric boundary layer grid factor (e.g. 1.01)
      REAL       T(JL)      ! Space for a 1-D distribution
      INTEGER    LUNERR     ! For error messages

C     Local constants:

      REAL, PARAMETER ::
     >   FRACTION = 0.05,   ! Of chord, for normal distance off TE
     >   ONE      = 1.,
     >   ZERO     = 0.

C     Local variables:

      INTEGER    IER, IM, I1, I2, J, JEVAL, K
      REAL       ARROW, CFRACTION, DX, DY, TTOTAL,
     >           XLINE(6), YLINE(6), TLINE(6)
      LOGICAL    NEW

C     Execution:

      IM = IEND
      IF (IM == 1) THEN  ! Below wing
         I1 = ITL
         I2 = ITL + 1
         ARROW = -ONE
      ELSE               ! IEND = IL (above wing)
         I1 = ITU
         I2 = ITU - 1
         ARROW = ONE
      END IF

      DX = X(I1,1) - X(ILE,1)  ! Basically chord, but allow for twist
      DY = Y(I1,1) - Y(ILE,1)
      CFRACTION = SQRT (DX*DX + DY*DY) * FRACTION

C     Construct a 6-point curve to be splined between TE and outer C boundary:

      XLINE(1) = X(I1,1)
      YLINE(1) = Y(I1,1)
      DY = CFRACTION * ARROW  ! ... and from similar right triangles:
      DX = -DY * ((YLINE(1) - Y(I2,1)) / (XLINE(1) - X(I2,1)))
      XLINE(2) = XLINE(1) + DX
      YLINE(2) = YLINE(1) + DY
      XLINE(3) = XLINE(2) + DX
      YLINE(3) = YLINE(2) + DY
      XLINE(4) = X(I1,JL)
      YLINE(4) = Y(IM,JL-2)
      XLINE(5) = XLINE(4)
      YLINE(5) = Y(IM,JL-1)
      XLINE(6) = XLINE(4)
      YLINE(6) = Y(IM,JL)

C     Unnormalized arc lengths:

      CALL CHORDS2D (6, XLINE, YLINE, .FALSE., TTOTAL, TLINE)

C     Composite geometric + Vinokur distribution:

      DX = D1TRAIL  ! Not dX, dY, but avoids new variables
      DY = TTOTAL - TLINE(5)
      T(1)  = ZERO  ! 1 & N reqd. in place by BLGRID
      T(JL) = TTOTAL

      CALL BLGRID (JL, DX, DY, NBLAYER, RBLAYER, T, LUNERR, IER)

      IF (IER /= 0 .AND. IER /= 3) THEN
         WRITE (LUNERR,*) 'TEBOUND: IER, N, D1, D2, TTOTAL = ',
     >                    IER, JL, DX, DY, TTOTAL
         STOP
      END IF

C     Distribute X and Y along the arc.  One point at a time avoids
C     a work-space vector.

      NEW = .TRUE.
      DO J = 2, JL - 1
         CALL LCSFIT (6, TLINE, XLINE, NEW, 'B', 1, T(J), X(I1,J), DX)
         NEW = .FALSE.
      END DO

      NEW = .TRUE.
      DO J = 2, JL - 1
         CALL LCSFIT (6, TLINE, YLINE, NEW, 'B', 1, T(J), Y(I1,J), DY)
         NEW = .FALSE.
      END DO

      RETURN
      END
