C*******************************************************************************
C
      SUBROUTINE SURFER
     >
     >  (NCALL, NWSURF, MXCRANK, MXIGRID, MXKGRID, MXINTR, MXPATCH,
     >   NNCONFIG, NNPANEL, NNWING, NNNAC, NNFUSE, NNSPLIT, NNROOT,
     >   NNTIP, NNTAIL, NNSURF, NWSEC, IQUAD, XTRAP, ZINNER, ZOUTER,
     >   LENRG, IDIM, KDIM, IRG, XRG, YRG, ZRG,
     >   NFSTN, JLBODY, XF, YF, ZF, IDEGFUS, IDEGREE,
     >   ILWING, KLWING, KWING1, KWING2,
     >   NCRANK, LCRANK, KCRANK, MCRANK, ZCRANK,
     >   SPNSPC, SPNSPS, TNORM, TBASE,
     >   TINTR, UINTR, VINTR, IINTR, JINTR, KINTR,
     >   NFIXL, NFIXR, NFIXA, NFIXB, IMXCOM, JMXCOM, IMXCOMO, JMXCOMO,
     >   MXSUBI, MXSUBJ, NSUBI, NSUBJ, ISUB, JSUB, ISUBO, JSUBO,
     >   LENXB, XBASE, YBASE, ZBASE, NEWPTS, NPATCH, LASTXB, MOSTXB,
     >   IXB, IPIECE, IWRIT, FAIL)
C
C     SURFER is the main routine for intersecting aircraft geometry components
C     and converting the outer mold line into regular surface patches suitable
C     for (u,v) parameterization.  Such a paneling of some initial intersected
C     components by the stand-alone AEROSURF, along with an initial multiblock
C     computational grid from GRIDGEN (say) and a corresponding (u,v) map from
C     program UVMAP, enables the computational grid to be "warped" efficiently
C     given another SURFER paneling of the perturbed surfaces.
C
C     SURFER is intended to be the portion of the surface patching scheme that
C     is common to AEROSURF and the aerodynamic design code.
C
C     Environment:  Fortran 90
C
C     Funding:  High Speed Research Program, Aerodynamics Division, NASA Ames
C
C     Sep. 1997  J. Reuther   Initial implementation, adding nacelles and
C                             pylons to WBSURF capabilities, but still no
C                             spanwise gridding of wing surfaces.
C     Oct. 1997: D. Saunders  General overhaul: separated geometry and grid
C     Nov. 1997               work arrays; allowed for less than max. grid
C                             sizes; generalized the connectivity scheme;
C                             installed arc-based spanwise gridding from
C                             SYN87-SB, allowing for end-of-pylon nacelle;
C                             separated out a general-purpose WBINTR for use
C                             in all intersection calculations.
C     Dec. 1997     DAS       Initial handling of under-wing nacelle/diverters.
C     Jan. 1998      "        Allowed for diverters beyond the trailing edge,
C                             and for wing patches aft of pylons (WDSURF).
C                             FSURF fixes up the fuselage near wing roots.
C                             Fortran 90 array syntax applied where reasonable.
C     Feb. 1998      "        Split interim patches as needed for continuous
C                             surface grid perturbations across patch edges
C                             (i.e., 4 edges at every corner).  Exceptions:
C                             > The edge corresp. to diverter LE is not carried
C                               into the fuselage paneling for now.
C                             > Splitting the wing anywhere but at the leading
C                               edge makes no sense for sharp sections, so we
C                               live with possible consequences at the root.
C                             > Nacelle patches along nacelle/pylon/diverter
C                               intersections don't meet at first/last sections.
C     Mar. 1998      "        NSURF broken up as part of paneling aft portions
C                             of HSCT-type nacelle/diverters.
C     Aug. 1998      "        NNSPLIT option provides a second way of splitting
C                             the interim body patches at wing leading edges.
C     Sep. 1998      "        FIXROOT confines (x,y,z) redistributions of the
C                             evaluations at the (u,v)s to the NFIXA/B region.
C                             Nacelles are now split at their leading edges.
C                             Regularized sections are input now, because they
C                             are needed by PERTURB to maintain diverter height.
C     Oct. 1998      "        MXJBODY argument was no longer needed; eliminated
C                             MXIBODY also by assuming the caller allocates
C                             regular fuselage arrays using NFSTN & NJBODY;
C                             the nacelle transformations to real space are now
C                             done in PERTURB for consistency.
C     11/18/98       "        FINSERT's search for ILE now allows for Xs that
C                             are not monotonic (high sweep + low resolution).
C     12/15/98       "        IPIECE(*) argument allows tagging patches.
C     01/31/99       "        X/Y/ZRG(*) are packed regularized inputs now.
C     02/09/99 +     "        X/Y/ZBASE(*) are packed output patches now.
C     03/10/99       "        MXINTR argument needed in addition to MXIGRID,
C                             because pylons/diverters can raise # intersn. pts.
C     04/03/99       "        Use of RECTIFY on LE points introduced via LCRANK.
C     04/15/99       "        BODYSRCH improved (1-variable projection method).
C     04/23/99       "        Started subpatch version; finished 05/15/99.
C     05/21/99 -     "        Vertical fin(s) handled, requiring NNTAIL(*).
C     06/15/99
C     07/13-20/00    "        FCTV case (body with 2 overlapping wing surfaces).
C     07/27-31/00    "        FCTV2 case (more specialized; splits the wing).
C     04/13/01       "        FBASE adds a fuselage base patch.
C     08/23/03       "        Slight glitch in FIXROOT has been corrected.
C
C*******************************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   NCALL,               ! -3 means first call
     >   NWSURF,              ! # "wing" + "nacelle" type surfaces (or 1 if 0)
     >   MXCRANK,             ! Max. # crank stations/wing (+ room for tip)
     >   MXIGRID,             ! Largest ILWING(*) value
     >   MXKGRID,             ! Max. # spanwise/azim. grid stns. (wing/nacel.)
     >   MXINTR,              ! MXINTR >= MXIGRID; WDSURF can need more
     >   MXPATCH,             ! Max. # gridded patches for all components
     >   NNCONFIG,            ! 0 = not a particularly special case: any case
     >                        !     handled prior to CTV (2 overlapping wings);
     >                        ! 1 = Crew Transfer Vehicle case ( "    "    " );
     >                        ! 2 = More specialized CTV case (splits the wing);
     >                        ! 3 = As for 2 but fin is on the center line;
     >                        ! 4 = <some future possibility>
     >   NNPANEL              ! Controls the paneling scheme:
                              ! 1 = parametric patches (variable point counts),
                              !     or nonparametric scheme but first call
                              !     (I/JSUBO(*,*) not used);
                              ! 2 = nonparametric scheme and not the first call:
                              !     point counts are forced to match those input
                              !     via I/JSUBO(*,*)
      INTEGER, INTENT (IN) ::
     >   NNWING,              ! # wing/pylon/canard/tail/diverter surfaces
     >   NNNAC,               ! # nacelle surfaces following NNWING surfaces
     >   NNFUSE,              ! 0 = no fuselage; 1 = fuselage present
     >   NNSPLIT,             ! 1 = orig. FSPLIT (all wing LEs); 2 = HSCT-type
     >   NNROOT(NWSURF)       ! See tricky usage below

!     Clarification of NNROOT(*) usage for controlling intersections:
!     ---------------------------------------------------------------
!
!        IW = 1 : NNWING refers to non-nacelle wing-type surfaces;
!        IW = NNWING + 1 : NNWING + NNNAC points to nacelle surfaces, if any.
!
!        NNROOT(IW) =  0 suppresses an intersection calculation for surface IW
!        NNROOT(IW) = +N >= +1 means surface IW intersects wing surface N
!        NNROOT(IW) = -1 means a wing(IW)/body intersection     (if NNFUSE = 1)
!        NNROOT(IW) = -2 means a vert.fin(IW)/body intersection ( "    "    " )
!        NNROOT(IW) = -N <= -3 means (tail) surface IW intersects v.fin surf. N

      INTEGER, INTENT (INOUT) ::
     >   NNTIP(NWSURF)        ! Nacelle indices derived from NNROOT on 1st call

      INTEGER, INTENT (IN) ::
     >   NNTAIL(NWSURF)       ! NNTAIL(IW) = +n means tail +n <->  dorsal fin IW
                              !    "   "   = -n   "    "   +n <-> ventral fin IW
                              ! where NNROOT(n) = -1 (body-mounted tail);
                              ! otherwise, NNTAIL(IW) = 0
      INTEGER, INTENT (IN) ::
     >   NNSURF,              ! Actual # "wing" + "nacelle" type surfaces >= 0
     >   NWSEC(NWSURF)        ! Actual # input "wing" defining sections; no
                              ! longer used for nacelles (KLWING() - 1 instead)
      INTEGER, INTENT (IN) ::
     >   IQUAD(2,NWSURF)      ! (1,*) = quadrant of 1st nacelle section;
                              ! (2,*) = quadrant nearest the intersection;
                              ! enter 3, 6, 9 or 12 o'clock, looking downstream

      REAL, INTENT (IN), DIMENSION (NWSURF) ::
     >   XTRAP                ! X distance beyond each surface edge at
                              ! which an extra station is extrapolated
                              ! to protect the intersection calculation;
                              ! for diverter IW, XTRAP(IW) is used for
                              ! extrapolating the WING trailing edge;
                              ! for nacelle IW, XTRAP(IW) > 0. means
                              ! extend the nacelle aft, while < 0. means
                              ! fudge the associated DIVERTER from
                              ! XTE - |XTRAP| to XTE, closing it at XTE
      REAL, INTENT (IN) ::
     >   ZINNER(NWSURF),      ! For XTRAP(IW) < 0. cases, these are the
     >   ZOUTER(NWSURF)       ! spanwise locations of the wing trailing
                              ! edge where the diverter associated with
                              ! nacelle IW is forced to close
      INTEGER, INTENT (IN) ::
     >   LENRG,               ! Length of X/Y/ZRG(*) arrays
     >   IDIM(NWSURF),        ! Effective dimensions of each wing-type surface
     >   KDIM(NWSURF),        ! packed into X/Y/ZRG(*)
     >   IRG(NWSURF)          ! IRG(IW) = start of regularized wing surface IW

      REAL, INTENT (INOUT), DIMENSION (LENRG) ::
     >   XRG, YRG, ZRG        ! Input with regularized wing & nacelle surfaces;
                              ! crank sections are inserted in place by SURFER
      INTEGER, INTENT (IN) ::
     >   NFSTN,               ! # body sections input; 1 if NNFUSE = 0
     >   JLBODY               ! # gridded pts. on body sections; "   "

      REAL, INTENT (IN) ::
     >   XF(NFSTN),           ! Current body section coordinates
     >   YF(JLBODY,NFSTN),    ! Regularized body sections
     >   ZF(JLBODY,NFSTN)

      INTEGER, INTENT (IN) ::
     >   IDEGFUS,             ! 1 or 3 for [bi]linear or [bi]cubic fus. splines
     >   IDEGREE(NWSURF),     ! 1 or 3 for bilinear or bicubic splines when
                              ! each surface is the "body" for an intersection
     >   ILWING(NWSURF),      ! # (wrap-around) chordwise gridded pts. per wing
     >   KLWING(NWSURF)       ! # gridded pts. per wing/nacelle spanwise/azim.

      INTEGER, INTENT (INOUT) ::
     >   KWING1(NWSURF)       ! Inboard geom. sections for w/b intersections,
                              ! forced to 0 here if NNFUSE or NNROOT(IW) = 0

      INTEGER, INTENT (IN), DIMENSION(NWSURF) ::
     >   KWING2,              ! Outboard geom. sections for w/b intersections
     >   NCRANK,              ! # wing surface cranks to include in paneling
     >   LCRANK               ! Crank # K >= 0 such that plain paneled wing LE
                              ! is rectified for grid stations 1:MCRANK(K);
                              ! LCRANK = NCRANK + 1 means all the way to the tip

      INTEGER, INTENT (OUT), DIMENSION (MXCRANK,NWSURF) ::
     >   KCRANK,              ! Geometry section indices of input cranks
     >   MCRANK               ! Mesh indices of crank stations

      REAL, INTENT (IN), DIMENSION (MXCRANK,NWSURF) ::
     >   ZCRANK               ! Zs of input crank stations (or
                              ! Ys if the surface is a vertical fin)

      REAL, INTENT (IN), DIMENSION (NWSURF) ::
     >   SPNSPC, SPNSPS       ! Spanwise gridding controls

      REAL, INTENT (INOUT), DIMENSION (MXKGRID,NWSURF) ::
     >   TNORM, TBASE         ! Spanwise distributions, set on first call

      REAL, INTENT (INOUT), DIMENSION (MXINTR,NWSURF) ::
     >   TINTR, UINTR, VINTR  ! Intersection restart data

      INTEGER, INTENT (INOUT), DIMENSION (MXINTR,NWSURF) ::
     >   IINTR, JINTR, KINTR  ! Intersection restart indices

      INTEGER, INTENT (IN), DIMENSION (NWSURF) ::
     >   NFIXL,               ! Controls for size of repatched fuselage
     >   NFIXR,               ! regions left & right of intersection LEs
     >   NFIXA,               ! and associated controls above and below
     >   NFIXB                ! intersection leading edges; see AEROSURF

      INTEGER, INTENT (INOUT), DIMENSION (MXPATCH) ::
     >   IMXCOM,  JMXCOM      ! In/out with old/new patch sizes

      INTEGER, INTENT (INOUT), DIMENSION (MXPATCH) ::
     >   IMXCOMO, JMXCOMO     ! Output with copy of old sizes if NNPANEL = 1;
                              ! input with specified sizes if NNPANEL = 2
      INTEGER, INTENT (IN) ::
     >   MXSUBI, MXSUBJ       ! Max. # subpatches per patch over all patches

      INTEGER, INTENT (INOUT), DIMENSION (MXPATCH) ::
     >   NSUBI, NSUBJ         ! Numbers of subpatches in each direction;
                              ! output if NCALL = -3 or NNPANEL = 1;
                              ! input if NNPANEL = 2; constant for a given
                              ! configuration either way

      INTEGER, INTENT (INOUT), DIMENSION (0:MXSUBI, MXPATCH) ::
     >   ISUB,                ! Interim subpatch dimensions = ISUB(0:NSUBI(*),*)
     >   ISUBO                ! Desired dimensions; ISUB(0,*) = ISUBO(0,*) = 1;
                              ! ISUB(*,*) can vary if NNPANEL = 1, in which case
                              ! ISUBO(*,*) is not used

      INTEGER, INTENT (INOUT), DIMENSION (0:MXSUBJ, MXPATCH) ::
     >   JSUB,                ! Analogous subpatch definitions in J direction;
     >   JSUBO                ! JSUB(0,*) = JSUBO(0,*) = 1

      INTEGER, INTENT (IN) ::
     >   LENXB                ! Length of X/Y/ZBASE(*) arrays

      REAL, INTENT (INOUT), DIMENSION (LENXB) ::
     >   XBASE, YBASE, ZBASE  ! Output with packed surface patches

      INTEGER, INTENT (OUT) ::
     >   NEWPTS,              ! 1 if NCALL <= 0 or I/JMXCOMO(*) /= I/JMXCOM(*)
     >   NPATCH,              ! Number of surface patches generated
     >   LASTXB,              ! Last element of last patch in X/Y/ZBASE(*)
     >   MOSTXB               ! Most room needed in X/Y/ZBASE(*) (>= LASTXB)

      INTEGER, INTENT (OUT), DIMENSION (MXPATCH) ::
     >   IXB,                 ! IXB(L) = start of patch L in X/Y/ZBASE(*)
     >   IPIECE               ! Component numbers for each patch

      INTEGER, INTENT (IN) ::
     >   IWRIT                ! For diagnostics

      LOGICAL, INTENT (OUT) ::
     >   FAIL                 ! The optimizer may be able to back up & try again

C     Local constants:

      INTEGER,   PARAMETER :: MXNNAC  = 2,  ! Local max. assumed for NNNAC
     >                        MXNWING = 8,  ! ... and for NNWING
     >                        MXWSURF = 10  ! ... and for NWSURF
      CHARACTER, PARAMETER :: NAME * 9 = ' SURFER: '
      REAL,      PARAMETER :: HALF = 0.5, ONE = 1., ZERO = 0.,
     >                        EPSY = 1.E-6  ! For diverter fudge on nacelle top

C     Local variables:

      INTEGER
     >   I, I1, I1R, I1X, I2, IBSPLIT(3), ID, ID2, IDIV(NNWING),
     >   IFIN, IFUS, IL, ILD, IN, IW, IWING, J, J1, J2, K, KL,
     >   KPATCH, KPATCH0, L, LWING, M, N, NBSPLIT, NDIV, NFP, NFUSE,
     >   NH, NI, NJBODY, NK, NL, NP, NU

      INTEGER, DIMENSION (NWSURF) ::
     >   NHALF, NINTR, NREGSEC

      REAL
     >   ANGLE, DS, YFUDGE, ZMIN, ZMAX, ZNEWC, ZNEWS, ZNEWT, ZUNIF,
     >   SNORM(MXINTR)

      REAL, DIMENSION (MXINTR,NWSURF) ::
     >   XINTR, YINTR, ZINTR

      LOGICAL
     >   DFUDGE, DORSAL, EXTRAP

C     Arrays below are short enough that they can easily be kept out of the
C     argument list for SURFER via hard-coding of dimensions certain to suffice.

      INTEGER
     >   ITL(MXWSURF),   ! Start/end of pylon/diverter intersections
     >   ITU(MXWSURF)

      INTEGER, SAVE ::
     >   ISTN1(MXNWING), ! Fuselage ranges for wing or vert. fin intersections
     >   ISTN2(MXNWING),
     >   KSTN1(MXWSURF), ! Wing section ranges for diverter/pylon intersections
     >   KSTN2(MXWSURF)

      INTEGER, DIMENSION (MXNNAC) ::
     >   IPDIVLE, IPDIVTE     ! Subpatch scheme needs these to find patch sizes

      REAL
     >   CLOSEDIV(3,2,MXNNAC) ! HSCT diverters are forced to close at wing TE

      LOGICAL
     >   BLUNT

      LOGICAL, SAVE ::
     >   DIVCLOSED(MXNWING),
     >   LOSTINTR(MXNWING) ! Wing pylons/diverters or rounded wing roots prevent
                           ! reuse of initial wing/body intersection indices
      REAL, SAVE ::
     >   EPS, PI

      CHARACTER, SAVE ::
     >   METHOD * 1

C     Execution:
C     ----------

      FAIL = .FALSE.           ! Needed for the wing-alone case

      IF (NNPANEL == 1) THEN   ! Probably the parametric scheme
         IF (NCALL /= -3) THEN ! Not the first call - parametric for sure
            IMXCOMO = IMXCOM   ! Copy the patch sizes from the previous call
            JMXCOMO = JMXCOM
         END  IF
      ELSE ! NNPANEL = 2: Input I/JMXCOMO and I/JSUBO are imposed before return
      END IF

      KPATCH = 1        ! Next (or current) patch number
      LASTXB = 0        ! Last index used in X/Y/ZBASE(*)
      MOSTXB = 0        ! Most space needed normally exceeds final LASTXB value

      DO N = 1, MXPATCH ! Lowest subpatch indices are always 1
         NSUBI(N)  = 1  ! In many cases - default it
         ISUB(0,N) = 1
         ISUB(1:MXSUBI,N) = 0 ! Avoid undefined indices appended to aero.xyz
         NSUBJ(N)  = 1  !
         JSUB(0,N) = 1
         JSUB(1:MXSUBJ,N) = 0
      END DO

C     First-call-only stuff:
C     ----------------------

      IF (NCALL == -3) THEN

         PI  = 4.* ATAN (ONE)
         EPS = MAX (5.* EPSILON (ONE), 1.E-7) ! For RIPPLE2D

         IF (NWSURF > MXWSURF .OR. NNWING > MXNWING) THEN
            IF (IWRIT > 0) THEN
               WRITE (IWRIT,'(/,A,A,4I5)') NAME,
     >            'NWSURF > MXWSURF or NNWING > MXNWING:',
     >            NWSURF, MXWSURF, NNWING, MXNWING
            END IF
            GO TO 990
         END IF

         DO IW = 1, NNWING

C           Normalized spanwise distributions:

            TNORM(1,IW) = ZERO
            KL = KLWING(IW)

            DO K = 2, KL - 1
               ZUNIF = REAL (K-1) / REAL (KL-1)
               ANGLE = PI * ZUNIF
               ZNEWS = SIN (ANGLE * HALF)
               ZNEWC = HALF * (ONE - COS (ANGLE))
               ZNEWT =
     >            (ONE - SPNSPC(IW)) * ZUNIF + SPNSPC(IW) * ZNEWC
               TNORM(K,IW) =
     >            (ONE - SPNSPS(IW)) * ZNEWT + SPNSPS(IW) * ZNEWS
            END DO

            TNORM(KL,IW) = ONE

            IF (NNFUSE == 0 .OR. NNROOT(IW) == 0) THEN
               KWING1(IW) = 1 !  This is needed for wing gridding
               TINTR(1:ILWING(IW),IW) = ZERO
            END IF

C           Any pylons or diverters attached to this wing surface?

            CALL GET_DIV_LIST (IW) ! Internal procedure

            LOSTINTR(IW) = NDIV > 0 ! Because WDSURF will fiddle with I/JINTR(*)

C           Also, if the root LE may be rectified, we can't reuse I/JINTR(*):

            IF (NNROOT(IW) == -1) THEN
               IF (NNFUSE > 0) THEN
                  IF (LCRANK(IW) > 0) LOSTINTR(IW) = .TRUE.
               END IF
            END IF

C           Any nacelle associated with this wing/pylon/diverter?

            NNTIP(IW) = 0
            DO IN = NNWING + 1, NNSURF
               IF (NNROOT(IN) == IW) THEN
                  NNTIP(IW) = IN
                  EXIT
               END IF
            END DO

            DIVCLOSED(IW) = .FALSE. ! Except for HSCT diverters (below)
         END DO

         IF (NNFUSE > 0) THEN ! Fuselage section interpolation method
            IF (IDEGFUS == 3) THEN
               METHOD = 'B' ! Loose local cubic spline fits
            ELSE IF (IDEGFUS == 2) THEN ! Kludge
               METHOD = 'M' ! Monotonic ...
            ELSE   ! IDEGFUS = 1
               METHOD = 'L' ! Piecewise linear
            END IF
         END IF

      END IF ! End of first-call stuff


C     Check for wing crank stations, inserting sections if necessary:
C     ---------------------------------------------------------------

      DO IW = 1, NNSURF

         IL = ILWING(IW)
         NHALF(IW) = (IL + 1) / 2

         IF (IW <= NNWING) THEN
            NK = NWSEC(IW)
         ELSE ! For nacelle sections after PERTURB transformations
            NK = KLWING(IW) - 1
         END IF

         IF (NCRANK(IW) > 0) THEN

            I1 = IRG(IW)
            I2 = I1 + IDIM(IW) ! Element (1,2) of ZRG(I,K)

            IF (ZRG(I2) /= ZRG(I1)) THEN ! It's not a vertical fin

              CALL WINSERT (NCALL, IW, 1, NK, IDIM(IW), KDIM(IW),
     >                      NHALF(IW), 1, IL, XRG(I1), YRG(I1), ZRG(I1),
     >                      NCRANK(IW), ZCRANK(1,IW), KCRANK(1,IW),
     >                      IWRIT)

            ELSE ! Switch Y and Z to treat cranks in a vertical fin

              CALL WINSERT (NCALL, IW, 1, NK, IDIM(IW), KDIM(IW),
     >                      NHALF(IW), 1, IL, XRG(I1), ZRG(I1), YRG(I1),
     >                      NCRANK(IW), ZCRANK(1,IW), KCRANK(1,IW),
     >                      IWRIT)
            END IF

         END IF

         NREGSEC(IW) = NK ! May be incremented upon return

      END DO


C     Cover the no-fuselage case: only under-wing diverters differ later.

      DO IW = 1, NNWING ! Needed for splitting wrapped patches after WSURF
         ITL(IW) = 1
         ITU(IW) = ILWING(IW)
      END DO

C     Calculate intersections with the fuselage:
C     ------------------------------------------

      IF (NNFUSE > 0) THEN

         IF (NNCONFIG == 0) THEN
            NJBODY = JLBODY - 1 ! Allow for inserting a point <-> root LE
         ELSE
            NJBODY = JLBODY     ! For FCTV, FCTV2, or FCTV3
         END IF

C        Wing/body and/or vertical-fin/body intersections:
C        -------------------------------------------------

         DO IW = 1, NNWING

            IL = ILWING(IW)
            I1 = IRG(IW)
            IWING = NNROOT(IW)

            IF (IWING == -2) THEN ! Vertical fin/fuselage intersection

               CALL VFINTR (NCALL, NNCONFIG, IW, LOSTINTR(IW),
     >                      KWING1(IW), KWING2(IW),
     >                      IDIM(IW), KDIM(IW), MXINTR, NHALF(IW),
     >                      IL, XRG(I1), YRG(I1), ZRG(I1),
     >                      NJBODY, NFSTN, JLBODY, NFSTN, XF, YF, ZF,
     >                      ISTN1(IW), ISTN2(IW), IDEGFUS, NINTR(IW),
     >                      XINTR(1,IW), YINTR(1,IW), ZINTR(1,IW),
     >                      UINTR(1,IW), VINTR(1,IW), TINTR(1,IW),
     >                      IINTR(1,IW), JINTR(1,IW), KINTR(1,IW),
     >                      IWRIT, FAIL)

               IF (FAIL) THEN
                  IF (IWRIT > 0) WRITE (IWRIT,'(/,A,A,I2)') NAME,
     >               'Body intersection failed for vertical fin', IW
                  GO TO 990
               END IF

            ELSE IF (IWING == -1) THEN ! Wing/fuselage

               CALL WFINTR (NCALL, NNCONFIG, IW, LOSTINTR(IW),
     >                      KWING1(IW), KWING2(IW),
     >                      IDIM(IW), KDIM(IW), MXINTR, NHALF(IW),
     >                      IL, XRG(I1), YRG(I1), ZRG(I1),
     >                      NJBODY, NFSTN, JLBODY, NFSTN, XF, YF, ZF,
     >                      ISTN1(IW), ISTN2(IW), IDEGFUS, NINTR(IW),
     >                      XINTR(1,IW), YINTR(1,IW), ZINTR(1,IW),
     >                      UINTR(1,IW), VINTR(1,IW), TINTR(1,IW),
     >                      IINTR(1,IW), JINTR(1,IW), KINTR(1,IW),
     >                      IWRIT, FAIL)

               IF (FAIL) THEN
                  IF (IWRIT > 0) WRITE (IWRIT,'(/,A,A,I2)') NAME,
     >               'Fuselage intersection failed for surface', IW
                  GO TO 990
               END IF

            END IF

         END DO

      END IF


C     Wing/diverter intersections and/or vertical fin/horiz.tail intersections:
C     -------------------------------------------------------------------------

      ID = 0 ! Diverter #

      DO IW = 1, NNWING

         IWING = NNROOT(IW) ! IWING = wing index if +ve, else -(vert. fin index)

         IF (IWING > 0) THEN ! IW = diverter/pylon index

            IL = ILWING(IW)
            I1 = IRG(IW)

            IF (NCALL == -3) THEN

C              Determine the wing section K range for this intersection.
C              Note that this applies after possible insertions at cranks.

               K = KWING1(IW)              ! 1st relevant diverter/pylon section
               L = (K - 1) * IDIM(IW) + I1 ! 1st pt. of section K
               J = L + IL - 1              ! Last  "   "   "

               ZMIN = ZRG(L)
               ZMAX = ZRG(J)

               DO I = 2, NHALF(IW)
                  L = L + 1
                  ZMIN = MIN (ZRG(L), ZMIN)
                  J = J - 1
                  ZMAX = MAX (ZRG(J), ZMAX)
               END DO

               L = IRG(IWING) ! First point of first section of wing IWING
               J = IDIM(IWING)

               DO K = 1, NREGSEC(IWING)
                  IF (ZMIN <= ZRG(L)) EXIT
                  L = L + J
               END DO
               KSTN1(IW) = K - 1

               L = IRG(IWING)

               DO K = 1, NREGSEC(IWING)
                  IF (ZMAX <= ZRG(L)) EXIT
                  L = L + J
               END DO
               KSTN2(IW) = K

            END IF

            L = IRG(IWING) ! First pt. of wing IWING

            CALL WDINTR (NCALL, IDIM(IW), KDIM(IW), MXINTR, NHALF(IW),
     >                   IL, XRG(I1), YRG(I1), ZRG(I1),
     >                   IDIM(IWING), KDIM(IWING), 1, NHALF(IWING),
     >                   KSTN1(IW), KSTN2(IW), KWING1(IW), KWING2(IW),
     >                   XRG(L), YRG(L), ZRG(L), XTRAP(IW),
     >                   IDEGREE(IWING), ITL(IW), ITU(IW), NINTR(IW),
     >                   XINTR(1,IW), YINTR(1,IW), ZINTR(1,IW),
     >                   UINTR(1,IW), VINTR(1,IW), TINTR(1,IW),
     >                   IINTR(1,IW), JINTR(1,IW), KINTR(1,IW),
     >                   IWRIT, FAIL)

            IF (FAIL) THEN
               IF (IWRIT > 0) WRITE (IWRIT,'(/,A,A,I2)') NAME,
     >            'Wing/diverter intersection failed for surface', IW
               GO TO 990
            END IF

C           Preserve diverter LE/TE indices (before MATCHDIV) so that the
C           subpatch scheme can determine the whole patch sizes in WDSURF.
C           These are just reasonably close if RECTIFY is used below.

            ID = ID + 1
            IPDIVTE(ID) = IINTR(ITL(IW),IW)
            IPDIVLE(ID) = IINTR(NHALF(IW),IW)

         ELSE IF (IWING == -3) THEN ! Vertical fin (-IWING) intersects
                                    ! horizontal tail (IW)
            IFIN = -IWING

            CALL VHINTR (NCALL, IFIN, IW, IWRIT, FAIL) ! Etc.

            IF (FAIL) THEN
               IF (IWRIT > 0)
     >            WRITE (IWRIT,'(/,A,A,I2,A,I2)') NAME,
     >               'Bad intersection for vertical fin', IFIN,
     >               ' & horizontal tail', IW
               GO TO 990
            END IF

         END IF

      END DO


C     Match on-wing point counts of opposing diverter intersections:
C     --------------------------------------------------------------

C     MATCHDIV redistributes half of some regularized diverter sections, so
C     this needs to be done before the diverter/nacelle intersections.

      DO IW = 1, NNWING

C        Find any diverters attached to this wing surface:

         CALL GET_DIV_LIST (IW)

         DO J = 2, NDIV

            ID  = IDIV(J-1)
            ID2 = IDIV(J)
            K   = IRG(IW)
            L   = IRG(ID2)

            CALL MATCHDIV (IW, ID, ID2, NWSURF, MXINTR,
     >                     IDIM(IW),  KDIM(IW),  XRG(K), YRG(K), ZRG(K),
     >                     IDIM(ID2), KDIM(ID2), XRG(L), YRG(L), ZRG(L),
     >                     XTRAP, EPS, IDEGREE, ILWING, NHALF, ITL, ITU,
     >                     KSTN1, KSTN2, KWING1, KWING2,
     >                     NINTR, XINTR, YINTR, ZINTR,
     >                     TINTR, UINTR, VINTR, IINTR, JINTR, KINTR,
     >                     IWRIT, FAIL)

            IF (FAIL) THEN
               IF (IWRIT > 0) WRITE (IWRIT,'(/,A,A,I2)') NAME,
     >            'Wing/diverter reintersection failed for div.', J
               GO TO 990
            END IF

         END DO

      END DO


C     Nacelle intersections:
C     ----------------------

      DO J = 1, NNNAC

         IN = NNWING + J ! Index of nacelle J sections
         ID = NNROOT(IN) !   "   "  corresp. diverter or pylon "wing" sections
         IL = ILWING(IN)
         KL = NREGSEC(IN)
         NH = NHALF(IN)

         IF (ID == 0) CYCLE ! Not much to do except NACSURF splits it at the LE


         IW = NNROOT(ID) ! Index of wing surface (or -1 for fuselage)
         EXTRAP = XTRAP(IN) > ZERO
         DFUDGE = XTRAP(IN) < ZERO .AND. IW > 0
         DIVCLOSED(ID) = IW > 0

         IF (EXTRAP) THEN ! Append nacelle hoop to protect intersection

            I1 = IRG(IN)
            I2 = I1 + IL
            IL = IL + 1

            DO K = 1, KL
               XRG(I2) = XRG(I2-1) + XTRAP(IN)
               YRG(I2) = YRG(I2-1) ! More precisely some day?
               ZRG(I2) = ZRG(I2-1)
               I2 = I2 + IDIM(IN)
            END DO

         END IF


C        Kludge for closing diverters at the wing trailing edge:
C        -------------------------------------------------------

         ILD = ILWING(ID)

         IF (DFUDGE) THEN

C           Locate where the wing trailing edge meets the nacelle.
C           (Internal procedure returns CLOSEDIV(1:3,1:2,J).)

            I1R = IRG(IW)

            CALL TE_POINTS (J, IDIM(IW), KDIM(IW), NREGSEC(IW),
     >                      XRG(I1R), YRG(I1R), ZRG(I1R)) 

C           Warp two parts of the diverter and wing/diverter intersection:

            I1R = IRG(ID)

            CALL WARP_DIVERTER (J, IDIM(ID), KDIM(ID), NREGSEC(ID),
     >                          XRG(I1R), YRG(I1R), ZRG(I1R)) 

         END IF


C        Intersect the nacelle with a diverter or pylon:

         K = IRG(IN)
         L = IRG(ID)

         CALL WNINTR (NCALL, KWING1(IN), KWING2(IN), IQUAD(1,IN),
     >                IDIM(ID), KDIM(ID), MXINTR, NHALF(ID), ILD,
     >                XRG(L), YRG(L), ZRG(L),
     >                IL, KL, IDIM(IN), KDIM(IN),
     >                XRG(K), YRG(K), ZRG(K),
     >                IDEGREE(IN), ITL(IN), ITU(IN),
     >                XINTR(1,IN), YINTR(1,IN), ZINTR(1,IN),
     >                UINTR(1,IN), VINTR(1,IN), TINTR(1,IN),
     >                IINTR(1,IN), JINTR(1,IN), KINTR(1,IN),
     >                IWRIT, FAIL)

         IF (FAIL) THEN
            IF (IWRIT > 0) WRITE (IWRIT,'(/,A,A,I2)')
     >         NAME, 'Intersection failed for nacelle', J
            GO TO 990
         END IF


C        Fudge two intersection pts. to include the last nacelle section?

         IF (EXTRAP) THEN

            IL = IL - 1 ! Suppress extrapolated nacelle section

            CALL CAPTURE_NACELLE_TE (J, IDIM(IN), KDIM(IN), XRG(K))

         END IF


C        Fudge the nacelle/diverter intersection off the wing trailing edge?

         IF (DFUDGE) THEN ! Else WSURF will fail for the diverter

            I1 = IRG(IN)
            I2 = I1 + IL / 2
            YFUDGE = EPSY * (XRG(I1) - XRG(I2))

            DO I = 1, ITL(ID)
               YINTR(I,IN) = YINTR(I,ID) - YFUDGE
            END DO

            DO I = ITU(ID), ILD
               YINTR(I,IN) = YINTR(I,ID) - YFUDGE
            END DO

         END IF

      END DO ! Next nacelle


C     Wing-type surface gridding:
C     ---------------------------

      IBSPLIT(1) = 0 ! Unless there are under-wing pylons or diverters

      DO IW = 1, NNWING

C        Any pylons/diverters attached to this wing surface?

         CALL GET_DIV_LIST (IW) ! Returns list IDIV(1:NDIV), where NDIV >= 0

         IL = ILWING(IW)
         KL = KLWING(IW)
         NH = NHALF(IW)
         IN = NNTIP(IW) ! Tip nacelle index or 0

         IF (IN > 0) THEN

            IF (DIVCLOSED(IW)) THEN ! IW is a diverter

C              Force the on-nacelle portion of the nacelle/diverter intersection
C              to end at the wing TE and undo the above fudge of YINTR(<te>,IN):

               ITL(IN) = ITL(IW)
               YINTR(ITL(IN),IN) = YINTR(ITL(IW),IW)
               ITU(IN) = ITU(IW)
               YINTR(ITU(IN),IN) = YINTR(ITU(IW),IW)
               J = IN - NNWING ! For CLOSEDIV(*,*,J) after call to WSURF
             END IF

         END IF

         I1R = IRG(IW)
         I2  = I1R + IDIM(IW) ! First pt. of second section

         IF (ZRG(I2) /= ZRG(I1R)) THEN ! It's not a vertical fin
            NI = IL
         ELSE
            NI = NH
         END IF

         L = KPATCH

         CALL XBCHECK (LENXB, MXPATCH, L, LASTXB, NI, KL, IW, IWRIT,
     >                 MOSTXB)

         I1X    = LASTXB + 1
         IXB(L) = I1X              ! Start of next patch
         LASTXB = LASTXB + NI * KL ! End of patch when done

         IF (NI == IL) THEN ! It's not a vertical fin

            CALL WSURF (NCALL, IW, IN, NNSURF, IDIM(IW), KDIM(IW),
     >                  MXCRANK, MXINTR, NREGSEC(IW), NH, IL, KL,
     >                  XRG(I1R), YRG(I1R), ZRG(I1R),
     >                  XINTR, YINTR, ZINTR, TINTR, KWING1, KWING2,
     >                  NCRANK(IW), KCRANK(1,IW), MCRANK(1,IW),
     >                  TNORM(1,IW), TBASE(1,IW),
     >                  XBASE(I1X), YBASE(I1X), ZBASE(I1X), IWRIT)

         ELSE ! NI = NHALF; panel just the left surface of a fin

            CALL VSURF (NCALL, IW, IN, NNSURF, IDIM(IW), KDIM(IW),
     >                  MXCRANK, MXINTR, NREGSEC(IW), NH, IL, KL,
     >                  XRG(I1R), YRG(I1R), ZRG(I1R),
     >                  XINTR, YINTR, ZINTR, TINTR, KWING1, KWING2,
     >                  NCRANK(IW), KCRANK(1,IW), MCRANK(1,IW),
     >                  TNORM(1,IW), TBASE(1,IW),
     >                  XBASE(I1X), YBASE(I1X), ZBASE(I1X), IWRIT)

            IMXCOM(L) = NH
            JMXCOM(L) = KL
            ISUB(1,L) = NH
            NSUBJ(L)  = NCRANK(IW) + 1

            DO K = 1, NSUBJ(L) ! Works for NCRANK = 0; + 1 is at the tip
               JSUB(K,L) = MCRANK(K,IW) ! Mesh crank stations
            END DO

            KPATCH = L + 1 ! Next available

         END IF

C        Force plain wing paneling to have leading edge at middle index?
C        ---------------------------------------------------------------

         IF (NNROOT(IW) == -1 .AND. NNFUSE > 0) THEN

            IF (LCRANK(IW) > 0) THEN

               K = MCRANK(LCRANK(IW),IW) ! Paneled stations 1:K are rectified

               I1 = I1R + NH - 1 ! Recover a normalized arc distribution

               CALL CHORDS2D (NH, XRG(I1), YRG(I1), .TRUE., DS, SNORM)

               CALL RECTIFY (1, K, KL, IL, NH, NH, SNORM,
     >                       XBASE(I1X), YBASE(I1X), ZBASE(I1X))

               I1 = I1X
               DO I = 1, IL
                  XINTR(I,IW) = XBASE(I1) ! I/JINTR(*) are lost
                  YINTR(I,IW) = YBASE(I1)
                  ZINTR(I,IW) = ZBASE(I1)
                  I1 = I1 + 1
               END DO

            END IF

         END IF

C        Split the wrap-around surface at the leading edge?
C        --------------------------------------------------

         IF (NDIV == 0 .AND. NI == IL) THEN ! It's a plain wing-type surface

C           Implicit splitting at the leading edge would be much simpler than
C           explicit splitting, but unfortunately both have drawbacks if the
C           point counts are adjusted later for the nonparametric case:
C           implicit is bad at the leading edge for a sharp wing surface;
C           explicit is less than ideal near a rounded leading edge.
C           Since some surfaces have both round and sharp leading edges, even
C           an input control could not help.  Trying to do it implicitly for the
C           parametric case would conflict with the strategy of fixing the point
C           counts as one extra step at the end for the nonparametric case.
C           Therefore the decision has to be in favor of explicit splitting.

C           The following splits either a plain wing or a diverter which
C           needs to be suppressed beyond the wing trailing edge.

            M  = L + 1            ! For upper patch after splitting
            NL = NH - ITL(IW) + 1 ! NL = NH for a plain wing
            NU = ITU(IW) - NH + 1 ! NU = NH  "  "  "

C           (1) Move upper portions to the next patch space, but pack them:

            CALL XBCHECK (LENXB, MXPATCH, M, LASTXB, NU, KL, IW, IWRIT,
     >                    MOSTXB)

            I1 = IXB(L) + NH - 1  ! Mid-pt. (LE) of first section
            I2 = LASTXB + 1       ! 1st element beyond interim WSURF patch

            CALL XBCOPY (LENXB, XBASE, YBASE, ZBASE,
     >                   I1, NU, KL, IL, I2, NU)

C           (2) Pack lower portions of all sections in place of whole patch:

            I2 = IXB(L)           ! Start of lower half after repacking
            I1 = I2 + ITL(IW) - 1 ! ... and before repacking

            CALL XBCOPY (LENXB, XBASE, YBASE, ZBASE,
     >                   I1, NL, KL, IL, I2, NL)

C           (3) Move the upper surface patch to the end of the lower patch:

            NP = NU * KL          ! Move all upper sections as one chunk
            I2 = I2 + NL * KL     ! Start of packed upper surface patch
            IXB(M) = I2
            I1 = LASTXB + 1       ! Start of interim  "     "     "
            LASTXB = I2 + NP - 1  ! When we're done

            CALL XBCOPY (LENXB, XBASE, YBASE, ZBASE,
     >                   I1, NP, 1, NP, I2, NP)

            IMXCOM(L) = NL
            JMXCOM(L) = KL
            ISUB(1,L) = NL

            IMXCOM(M) = NU
            JMXCOM(M) = KL
            ISUB(1,M) = NU

C           The following implicit splits are needed if REPACK is to be used.
C           Even though WSURF always captures planform cranks precisely, they
C           would be blurred if REPACK redistributed across cranks.

            NSUBJ(L) = NCRANK(IW) + 1 ! NSUBI/J are defaulted at 1 at the start
            NSUBJ(M) = NSUBJ(L)

            DO K = 1, NSUBJ(L) ! Works for NCRANK = 0; + 1 is at the tip
               JSUB(K,L) = MCRANK(K,IW) ! Mesh crank stations
               JSUB(K,M) = MCRANK(K,IW)
            END DO

            KPATCH = M + 1 ! Next available

            IF (DIVCLOSED(IW)) THEN ! Ensure a singular pt. at the trailing edge

               I1 = IXB(L)
               I2 = I2 + NU - 1

               DO K = 1, KL
                  XBASE(I1) = CLOSEDIV(1,1,J)
                  YBASE(I1) = CLOSEDIV(2,1,J)
                  ZBASE(I1) = CLOSEDIV(3,1,J)
                  I1 = I1 + NL
                  XBASE(I2) = CLOSEDIV(1,2,J)
                  YBASE(I2) = CLOSEDIV(2,2,J)
                  ZBASE(I2) = CLOSEDIV(3,2,J)
                  I2 = I2 + NU
               END DO

            END IF

         ELSE IF (NI == IL) THEN ! NDIV > 0: Repanel around the pylons/diverters

            CALL WDSURF (IW, NDIV, IDIV, NWSURF, KPATCH, MXPATCH,
     >                   MXINTR, NCRANK(IW), MCRANK(1,IW), IL, KL,
     >                   NNFUSE, EPS, IPDIVLE, IPDIVTE, LENXB,
     >                   XBASE, YBASE, ZBASE, LASTXB, MOSTXB, IXB,
     >                   MXSUBI, MXSUBJ, NSUBI, NSUBJ, ISUB, JSUB, 
     >                   NINTR, XINTR, YINTR, ZINTR, IINTR, JINTR,
     >                   NHALF, ITL, ITU, NBSPLIT, IBSPLIT,
     >                   IMXCOM, JMXCOM, LWING, IWRIT, FAIL)

            IF (FAIL) THEN
               IF (IWRIT > 0) WRITE (IWRIT,'(/,A,A,A,I2)')
     >            NAME, 'Pylon/diverter repaneling failed for ',
     >            'wing surface', IW
               GO TO 990
            END IF

         END IF

         IPIECE(L:KPATCH-1) = IW ! Geometry component number

      END DO


C     Nacelle paneling:
C     -----------------

      DO J = 1, NNNAC

         L = KPATCH  ! Next available
         KPATCH0 = L

         IN = NNWING + J ! Index of nacelle J sections
         IL = ILWING(IN)
         ID = NNROOT(IN) !   "   "  corresp. diverter or pylon sections

         IF (ID /= 0) THEN
            ILD = ILWING(ID)
            KL  = KLWING(IN) - 1 ! See nacelle loop above - allows inserting
            DFUDGE = DIVCLOSED(ID)
         ELSE                    ! Plain nacelle
            ILD = 1
            KL  = NREGSEC(IN)
            DFUDGE = .FALSE.
         END IF

C        Original nacelle paneling, possibly with aft portion suppressed:

         I1R = IRG(IN)

         CALL NACSURF (IN, ID, L, MXPATCH, IDIM(IN), KDIM(IN), IL, KL,
     >                 DFUDGE, XRG(I1R), YRG(I1R), ZRG(I1R),
     >                 NHALF(IN), ITL(IN), ITU(IN), ILD,
     >                 XINTR(1,IN), YINTR(1,IN), ZINTR(1,IN),
     >                 IINTR(1,IN), JINTR(1,IN),
     >                 LENXB, XBASE, YBASE, ZBASE, LASTXB, MOSTXB, IXB,
     >                 MXSUBI, MXSUBJ, NSUBJ, ISUB, JSUB,
     >                 IMXCOM, JMXCOM, IWRIT)

C        Special handling of HSCT-type aft nacelle:

         IF (DFUDGE) THEN

            CALL NACSURF2 (J, L, MXPATCH, IDIM(IN), KDIM(IN),
     >                     XRG(I1R), YRG(I1R), ZRG(I1R),
     >                     IL, KL, IL/2, LWING, CLOSEDIV(1,1,J), EPS,
     >                     LENXB, XBASE, YBASE, ZBASE, LASTXB, MOSTXB,
     >                     IXB, MXSUBI, MXSUBJ, NSUBJ, ISUB, JSUB,
     >                     IMXCOM, JMXCOM, IWRIT, FAIL)

            IF (FAIL) THEN
               IF (IWRIT > 0) WRITE (IWRIT,'(/,A,A,I2)')
     >            NAME, 'NACSURF2 failed for nacelle', J
               GO TO 990
            END IF

            LWING = LWING + 1 ! Upper wing patch # for next nacelle

         END IF

         KPATCH = L
         IPIECE(KPATCH0:L-1) = IN

      END DO


C     Fuselage paneling:
C     ------------------

      IF (NNFUSE > 0) THEN

         NFUSE = KPATCH ! First fuselage patch, needed by REPACK
         IFUS  = NNSURF + 1 ! Component number

         SELECT CASE (NNCONFIG) 

         CASE (0) ! All cases prior to CTV

            CALL FSURF (NWSURF, MXINTR, MXPATCH, KPATCH, NFSTN, NNWING,
     >                  JLBODY, NNROOT, NNSPLIT, NBSPLIT, IBSPLIT,
     >                  NNTAIL, IFUS, METHOD, EPS, XF, YF, ZF, MXNWING,
     >               LOSTINTR, NINTR, XINTR, YINTR, ZINTR, IINTR, JINTR,
     >                  ISTN1, ISTN2, NFIXL, NFIXR, NFIXA, NFIXB,
     >                  LENXB, XBASE, YBASE, ZBASE, LASTXB, MOSTXB, IXB,
     >                  MXSUBI, MXSUBJ, NSUBI, NSUBJ, ISUB, JSUB,
     >                  IMXCOM, JMXCOM, IWRIT, FAIL)

         CASE (1) !  CTV case with two overlapping wings (general; may fail)

C           Note switch of U/VINTR (not I/JINTR) to match WFINTR usage.

            CALL FCTV (NWSURF, MXINTR, MXPATCH, KPATCH, NFSTN,
     >                 JLBODY, IFUS, IDEGFUS, EPS, XF, YF, ZF,
     >                 NFIXL, NFIXR, NFIXA, NFIXB, NINTR,
     >                 XINTR, YINTR, ZINTR, VINTR, UINTR, IINTR, JINTR,
     >                 LOSTINTR, LENXB, XBASE, YBASE, ZBASE,
     >                 LASTXB, MOSTXB, IXB, MXSUBI, MXSUBJ,
     >                 NSUBJ, ISUB, JSUB, IMXCOM, JMXCOM, IWRIT, FAIL)

         CASE (2) !  CTV case with fin much smaller than wing (specialized)

C           Note switch of U/VINTR (not I/JINTR) to match WFINTR usage.

            CALL FCTV2 (NWSURF, MXINTR, MXPATCH, KPATCH, NFSTN,
     >                  JLBODY, IFUS, IDEGFUS, EPS, XF, YF, ZF,
     >                  NFIXL, NFIXR, NFIXA, NFIXB, NINTR,
     >                  XINTR, YINTR, ZINTR, VINTR, UINTR, IINTR, JINTR,
     >                  LOSTINTR, LENXB, XBASE, YBASE, ZBASE,
     >                  LASTXB, MOSTXB, IXB, MXSUBI, MXSUBJ,
     >                  NSUBI, NSUBJ, ISUB, JSUB, IMXCOM, JMXCOM,
     >                  IWRIT, FAIL)

         CASE (3) !  As for (2) but with fin on the center line

C           Note switch of U/VINTR (not I/JINTR) to match WFINTR usage.

            CALL FCTV3 (NWSURF, MXINTR, MXPATCH, KPATCH, NFSTN,
     >                  JLBODY, IFUS, IDEGFUS, EPS, XF, YF, ZF,
     >                  NFIXL, NFIXR, NFIXA, NFIXB, NINTR,
     >                  XINTR, YINTR, ZINTR, VINTR, UINTR, IINTR, JINTR,
     >                  LOSTINTR, LENXB, XBASE, YBASE, ZBASE,
     >                  LASTXB, MOSTXB, IXB, MXSUBI, MXSUBJ,
     >                  NSUBI, NSUBJ, ISUB, JSUB, IMXCOM, JMXCOM,
     >                  IWRIT, FAIL)

         END SELECT

         IF (FAIL) THEN
            IF (IWRIT > 0)
     >         WRITE (IWRIT,'(/,A,A)') NAME, 'Fuselage paneling failed.'
            GO TO 990
         END IF

C        Add a fuselage base patch (suppressed if singular point):

         CALL FBASE (MXPATCH, KPATCH, NFSTN, JLBODY, IFUS,
     >               XF, YF, ZF, LENXB, XBASE, YBASE, ZBASE,
     >               LASTXB, MOSTXB, IXB, MXSUBI, MXSUBJ,
     >               NSUBI, NSUBJ, ISUB, JSUB, IMXCOM, JMXCOM, IWRIT)

         IPIECE(NFUSE:KPATCH-1) = NNSURF + 1

      ELSE
         NFUSE = 999 ! No first fuselage patch
      END IF


      NPATCH = KPATCH - 1
      NEWPTS = 0

      IF (NNPANEL == 1) THEN ! Probably the parametric scheme

         IF (NCALL <= 0) THEN
            NEWPTS = 1
         ELSE ! NCALL > 0 means geometry was perturbed for a gradient calc.
            DO N = 1, NPATCH
               IF (IMXCOMO(N) /= IMXCOM(N)) NEWPTS = 1
               IF (JMXCOMO(N) /= JMXCOM(N)) NEWPTS = 1
            END DO
         END IF

      ELSE ! NNPANEL == 2: Nonparametric scheme: fix the subpatch dimensions

         CALL REPACK (LENXB, XBASE, YBASE, ZBASE, LASTXB, MOSTXB,
     >                NPATCH, NFUSE, MXSUBI, MXSUBJ, NSUBI, NSUBJ,
     >                ISUB, JSUB, ISUBO, JSUBO, IXB, IWRIT)

         DO N = 1, NPATCH
            IMXCOM(N) = IMXCOMO(N)
            JMXCOM(N) = JMXCOMO(N)
            I = NSUBI(N)
            J = NSUBJ(N)
            ISUB(1:I,N) = ISUBO(1:I,N)
            JSUB(1:J,N) = JSUBO(1:J,N)
         END DO

      END IF

      GO TO 999

 990  IF (NCALL == -3 .OR. NCALL > 0) THEN
         CALL SYNCH
         STOP
      END IF

 999  RETURN


C     SURFER internal procedures:

      CONTAINS

!        --------------------------------------------------
         SUBROUTINE CAPTURE_NACELLE_TE (J, IDIM, KDIM, XRG)
!        --------------------------------------------------

!        Fudge nacelle/pylon intersection points to capture the nacelle
!        trailing edge.  The nearest upper and lower intersection points
!        are moved to avoid changing the point count.

!        Arguments:

         INTEGER, INTENT (IN) ::
     >      J,                  ! Nacelle number starting at 1
     >      IDIM, KDIM          ! Nacelle surface array dimensions

         REAL, INTENT (IN), DIMENSION (IDIM,KDIM) ::
     >      XRG                 ! Regularized nacelle sections

!        Local variables:

         INTEGER I, I1, I2, INC, JLEFT
         REAL    R, RM1

!        Execution:

!        If either upper or lower, do both:

         IF (IINTR(1,IN) == IL .OR. IINTR(ILD,IN) == IL) THEN

            I1  = 1         ! Lower trailing edge
            I2  = ILD / 2
            INC = 1

            DO ! Two passes

               DO I = I1, I2, INC
                  IF (IINTR(I,IN) < IL) EXIT
               END DO

               JLEFT = JINTR(I,IN)
               R = (XRG(IL,JLEFT) - XINTR(I+INC,IN)) / ! Extrapolate in
     >               (XINTR(I,IN) - XINTR(I+INC,IN))   ! case I = I1
               RM1 = ONE - R
               XINTR(I,IN) = RM1 * XINTR(I+INC,IN) + R * XINTR(I,IN)
               YINTR(I,IN) = RM1 * YINTR(I+INC,IN) + R * YINTR(I,IN)
               ZINTR(I,IN) = RM1 * ZINTR(I+INC,IN) + R * ZINTR(I,IN)
               TINTR(I,IN) = RM1 * TINTR(I+INC,IN) + R * TINTR(I,IN)

               IF (I1 == 1) THEN ! Likewise for the upper trailing edge
                  ITL(IN) = I
                  I1  = ILD
                  INC = -1
               ELSE
                  ITU(IN) = I
                  EXIT
               END IF

            END DO ! The two-pass loop

         END IF

         END SUBROUTINE CAPTURE_NACELLE_TE

!        ----------------------------
         SUBROUTINE GET_DIV_LIST (IW)
!        ----------------------------

!        Locate pylon(s)/diverter(s) attached to wing surface IW.

!        Argument:

         INTEGER IW ! Wing surface index in [1, NNWING]

!        Local variables:

         INTEGER I

!        Execution:

         NDIV = 0
         DO I = 1, NNWING
            IF (NNROOT(I) == IW) THEN
               NDIV = NDIV + 1
               IDIV(NDIV) = I
            END IF
         END DO

         END SUBROUTINE GET_DIV_LIST

!        -------------------------------------------------------------
         SUBROUTINE INTERP_WING (IDIM, KDIM, I1, I2, K1, K2, IP, KP,
     >                           XRG, YRG, ZRG, XINT, ZINT, YINT, IER)
!        -------------------------------------------------------------

!        Project point (XINT, ZINT) onto the lower surface of wing IW.

!        Arguments:

         INTEGER, INTENT (IN) ::
     >      IDIM, KDIM,           ! Regularized wing array dimensions
     >      I1, I2, K1, K2        ! Subsurface of the wing to search
         INTEGER, INTENT (INOUT) ::
     >      IP, KP                ! Input cell guesses, output for reuse
         REAL, INTENT (IN), DIMENSION (IDIM,KDIM) ::
     >      XRG, YRG, ZRG         ! Regularized wing
         REAL, INTENT (IN) ::
     >      XINT, ZINT            ! (X,Z) coords. of point to be projected
         REAL, INTENT (OUT) ::
     >      YINT                  ! Wing surface Y (bilinearly interpolated)
         INTEGER, INTENT (OUT) ::
     >      IER                   ! Probably fatal if not 0

!        Local variables:

         REAL
     >      P, PM1, Q, QM1

!        Execution:

         CALL RIPPLE2D (IDIM, KDIM, I1, I2, K1, K2, XRG, ZRG,
     >                  XINT, ZINT, IP, KP, EPS, P, Q, IER)

         IF (IER /= 0) THEN
            IF (IWRIT > 0) THEN
               WRITE (IWRIT, '(/,A,I3,/,A,2F11.4,/,A,6I5)')
     >            ' INTERP_WING: Search failure for wing surface', IW,
     >            ' Target X and Z:', XINT, ZINT,
     >            ' I1, I2, K1, K2, IP, KP:', I1, I2, K1, K2, IP, KP,
     >            ' Proceeding.'
            END IF
         END IF

         PM1  = ONE - P
         QM1  = ONE - Q
         YINT = QM1 * (PM1*YRG(IP,KP)   + P*YRG(IP+1,KP))  +
     >          Q   * (PM1*YRG(IP,KP+1) + P*YRG(IP+1,KP+1))

         END SUBROUTINE INTERP_WING

!        ----------------------------------------------------------
         SUBROUTINE TE_POINTS (J, IDIM, KDIM, NWSEC, XRG, YRG, ZRG)
!        ----------------------------------------------------------

!        Locate the wing trailing edge points where diverter J is to be
!        artificially terminated with zero height.
!        Ideally, a wing TE/nacelle intersection would serve if the nacelle
!        isn't clipped off above the wing.  Instead, this version depends on
!        the (ugly) ZINNER and ZOUTER inputs.

!        Arguments:

         INTEGER, INTENT (IN) ::
     >      J,                  ! Nacelle/diverter number starting at 1
     >      IDIM, KDIM,         ! Wing surface array dimensions
     >      NWSEC               ! # regularized wing sections

         REAL, INTENT (IN), DIMENSION (IDIM,KDIM) ::
     >      XRG, YRG, ZRG       ! Regularized wing sections

!        Local variables:

         INTEGER I, K
         REAL    R

!        Execution:

         I = ITL(ID)
         DO K = JINTR(I,ID), NWSEC
            IF (ZINNER(IN) < ZRG(1,K)) EXIT
         END DO

         K = K - 1
         R = (ZINNER(IN) - ZRG(1,K)) / (ZRG(1,K+1) - ZRG(1,K))
         CLOSEDIV(3,1,J) = ZINNER(IN)
         CLOSEDIV(2,1,J) = (ONE - R) * YRG(1,K) + R * YRG(1,K+1)
         CLOSEDIV(1,1,J) = (ONE - R) * XRG(1,K) + R * XRG(1,K+1)

         DO K = JINTR(I,ID), NWSEC
            IF (ZOUTER(IN) < ZRG(1,K)) EXIT
         END DO

         K = K - 1
         R = (ZOUTER(IN) - ZRG(1,K)) / (ZRG(1,K+1) - ZRG(1,K))
         CLOSEDIV(3,2,J) = ZOUTER(IN)
         CLOSEDIV(2,2,J) = (ONE - R) * YRG(1,K) + R * YRG(1,K+1)
         CLOSEDIV(1,2,J) = (ONE - R) * XRG(1,K) + R * XRG(1,K+1)

         END SUBROUTINE TE_POINTS

!        --------------------------------------------------------------
         SUBROUTINE WARP_DIVERTER (J, IDIMD, KDIMD, NDSEC,
     >                             XDIV, YDIV, ZDIV)
!        --------------------------------------------------------------

!        Artificial closing of diverter at wing trailing edge.

!        Arguments:

         INTEGER, INTENT (IN) ::
     >      J,                  ! Nacelle/diverter number starting at 1
     >      IDIMD, KDIMD,       ! Diverter surface array dimensions
     >      NDSEC               ! # regularized diverter sections

         REAL, INTENT (INOUT), DIMENSION (IDIMD,KDIMD) ::
     >      XDIV, YDIV, ZDIV    ! Regularized diverter sections

!        Local constants:

         REAL,      PARAMETER :: DERIVS = -999. ! Suppress them
         LOGICAL,   PARAMETER :: CLOSED = .FALSE., NORM = .FALSE.
         CHARACTER, PARAMETER :: METHOD * 1 = 'B'

!        Local variables:

         INTEGER
     >      I, I1W, I2, IER, IL, ILEFT, IP, IT, IU, K, K1, K2, KP,
     >      MI, NI, MK, NK
         REAL
     >      XLEFT
         REAL, DIMENSION (IDIMD) ::
     >      XTMP, YTMP, ZTMP
         LOGICAL
     >      NEW

!        Execution:

         NI = NHALF(ID)
         NK = NDSEC
         IL = ITL(ID)

!        Inboard diverter surface:
!        -------------------------

!        Warp each diverter section to go through the nacelle/wing TE intersecn.

!        Determine the on-wing index at which to start the fudge:

         XLEFT = XINTR(IL,ID) + XTRAP(IN) ! XTRAP(IN) < 0.
         ILEFT = IL + 2

         CALL INTERVAL (NI, XINTR(1,ID), XLEFT, -ONE, ILEFT)

         ILEFT = MAX (ILEFT, IL + 2)

         XTMP(IL) = CLOSEDIV(1,1,J)
         ZTMP(IL) = CLOSEDIV(3,1,J)

         DO K = 1, NK
            XTMP(1)     = XDIV(1,K)  ! Copy end points
            YTMP(1)     = YDIV(1,K)
            ZTMP(1)     = ZDIV(1,K)
            YTMP(IL)    = YDIV(IL,K) ! ... and rest of TE point
            XTMP(ILEFT) = XDIV(ILEFT,K)
            YTMP(ILEFT) = YDIV(ILEFT,K)
            ZTMP(ILEFT) = ZDIV(ILEFT,K)

            CALL NULINE3D (1, IL,
     >                     XDIV(1,K), YDIV(1,K), ZDIV(1,K),
     >                     XTMP, YTMP, ZTMP)

            CALL NULINE3D (IL, ILEFT,
     >                     XDIV(1,K),  YDIV(1,K),  ZDIV(1,K),
     >                     XTMP, YTMP, ZTMP)

            DO I = 1, ILEFT ! Warped section
               XDIV(I,K) = XTMP(I)
               YDIV(I,K) = YTMP(I)
               ZDIV(I,K) = ZTMP(I)
            END DO

         END DO

!        Reintersect this part of the diverter with the wing:

!        TOO AWKWARD! (without repeated reparameterization of the wing, etc.)
!        For now, apply the above method to the intersection:

         XTMP(1)     = XINTR(1,ID)
         YTMP(1)     = YINTR(1,ID)
         ZTMP(1)     = ZINTR(1,ID)
         YTMP(IL)    = CLOSEDIV(2,1,J)
         XTMP(ILEFT) = XINTR(ILEFT,ID)
         YTMP(ILEFT) = YINTR(ILEFT,ID)
         ZTMP(ILEFT) = ZINTR(ILEFT,ID)

         CALL NULINE3D (1, IL,
     >                  XINTR(1,ID), YINTR(1,ID), ZINTR(1,ID),
     >                  XTMP, YTMP, ZTMP)

         CALL NULINE3D (IL, ILEFT,
     >                  XINTR(1,ID), YINTR(1,ID), ZINTR(1,ID),
     >                  XTMP, YTMP, ZTMP)

         DO I = 1, ILEFT
            XINTR(I,ID) = XTMP(I)
            YINTR(I,ID) = YTMP(I)
            ZINTR(I,ID) = ZTMP(I)
         END DO

!        Project the warped intersection portion onto the wing:

         I2 = NHALF(IW)
         K1 = KSTN1(ID)
         K2 = KSTN2(ID)
         KP = (K1 + K2) / 2
         IP = 1

         I1W = IRG(IW)

         DO I = IL + 1, ILEFT - 1

            CALL INTERP_WING (IDIM(IW), KDIM(IW), 1, I2, K1, K2, IP, KP,
     >                        XRG(I1W), YRG(I1W), ZRG(I1W), XINTR(I,ID),
     >                        ZINTR(I,ID), YINTR(I,ID), IER)
         END DO


!        Outboard diverter surface:
!        --------------------------

!        Determine the on-wing index at which to start the fudge:

         IU = ITU(ID)
         XLEFT = XINTR(IU,ID) + XTRAP(IN) ! XTRAP(IN) < 0.
         ILEFT = IU - 2 - NI

         CALL INTERVAL (NI, XINTR(NI,ID), XLEFT, ONE, ILEFT)

         ILEFT = MIN (ILEFT + NI, IU - 2)

         XTMP(IU) = CLOSEDIV(1,2,J)
         ZTMP(IU) = CLOSEDIV(3,2,J)
         IT = ILWING(ID)

         DO K = 1, NK
            XTMP(ILEFT) = XDIV(ILEFT,K)
            YTMP(ILEFT) = YDIV(ILEFT,K)
            ZTMP(ILEFT) = ZDIV(ILEFT,K)
            YTMP(IU)    = YDIV(IU,K)
            XTMP(IT)    = XDIV(IT,K)
            YTMP(IT)    = YDIV(IT,K)
            ZTMP(IT)    = ZDIV(IT,K)

            CALL NULINE3D (ILEFT, IU,
     >                     XDIV(1,K),  YDIV(1,K),  ZDIV(1,K),
     >                     XTMP, YTMP, ZTMP)

            CALL NULINE3D (IU, IT,
     >                     XDIV(1,K),  YDIV(1,K),  ZDIV(1,K),
     >                     XTMP, YTMP, ZTMP)

            DO I = ILEFT, IT
               XDIV(I,K) = XTMP(I)
               YDIV(I,K) = YTMP(I)
               ZDIV(I,K) = ZTMP(I)
            END DO

         END DO

!        For now, apply the above method to the intersection:

         XTMP(ILEFT) = XINTR(ILEFT,ID)
         YTMP(ILEFT) = YINTR(ILEFT,ID)
         ZTMP(ILEFT) = ZINTR(ILEFT,ID)
         YTMP(IU)    = CLOSEDIV(2,2,J)
         XTMP(IT)    = XINTR(IT,ID)
         YTMP(IT)    = YINTR(IT,ID)
         ZTMP(IT)    = ZINTR(IT,ID)

         CALL NULINE3D (ILEFT, IU,
     >                  XINTR(1,ID), YINTR(1,ID), ZINTR(1,ID),
     >                  XTMP, YTMP, ZTMP)

         CALL NULINE3D (IU, IT,
     >                  XINTR(1,ID), YINTR(1,ID), ZINTR(1,ID),
     >                  XTMP, YTMP, ZTMP)

         DO I = ILEFT, IT
            XINTR(I,ID) = XTMP(I)
            YINTR(I,ID) = YTMP(I)
            ZINTR(I,ID) = ZTMP(I)
         END DO

!        Project the warped intersection portion onto the wing:

         KP = K2 - 1

         DO I = ILEFT + 1, IU - 1

            CALL INTERP_WING (IDIM(IW), KDIM(IW), 1, I2, K1, K2, IP, KP,
     >                        XRG(I1W), YRG(I1W), ZRG(I1W), XINTR(I,ID),
     >                        ZINTR(I,ID), YINTR(I,ID), IER)
         END DO

         END SUBROUTINE WARP_DIVERTER

      END SUBROUTINE SURFER

C*******************************************************************************
C
      SUBROUTINE BODYSRCH (IDIM, JDIM, I1B, I2B, J1B, J2B, XB, YB, ZB,
     >                     XT, YT, ZT, IT, JT, IER)
C
C     BODYSRCH performs a specialized 2-D search within a fuselage-type
C     surface mesh consisting of sections at constant X stations, as
C     needed to determine surface cell indices corresponding to a target
C     point (XT,YT,ZT) on or near the surface.  It was prompted by a need
C     to redistribute points along part of a wing/body intersection curve.
C
C     This version effectively constructs a body section linearly at X = XT,
C     and determines its J interval giving the shortest perpendicular distance
C     from (YT,ZT).  This is consistent with FINSERT's construction of a new
C     section at the root leading edge.  (The nonlinear method of projecting
C     onto quad cells may give trouble on tapered bodies.)
C
C     01/08/98  DAS  Effort to avoid the more proper but less convenient
C                    parametric surface/line intersection calculation.
C     04/15/99   "   Rectification of root leading edges seemed to produce
C                    LE points further off the bilinear fuselage.  Switched
C                    from Y and Z comparisons to projection method in 2-space.
C
C     Author: David Saunders, Sterling Software/NASA Ames, Moffett Field, CA
C
C*******************************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER IDIM, JDIM         ! I   Dimensions of "body" surface arrays
      INTEGER I1B, I2B, J1B, J2B ! I   Active region to search
      REAL    XB(IDIM),          ! I   Regularized "body" sections
     >        YB(JDIM,IDIM),
     >        ZB(JDIM,IDIM)
      REAL    XT, YT, ZT         ! I   Target surface point
      INTEGER IT, JT             ! I/O Target cell lower left indices (input
                                 !     with starting guesses assumed in-range)
      INTEGER IER                !   O 0 means the cell was found
                                 !     1 means the J-direction search failed
C     Local constants:

      REAL, PARAMETER :: ONE = 1., ZERO = 0.

C     Local variables:

      INTEGER    I, IC, ISHIFT, J, JBEST, JC, JINC, NI, NJ
      REAL       DSQ, DSQBEST, R

C     Procedures:

      EXTERNAL   INTERVAL
     
C     Execution:

C     Locate the two body sections to use for interpolating one at X = XT:

      ISHIFT = I1B - 1
      NI = I2B - ISHIFT
      IC = IT  - ISHIFT

      CALL INTERVAL (NI, XB(I1B), XT, ONE, IC)

      IT = IC + ISHIFT
      I  = IT
      R  = (XT - XB(I)) / (XB(I+1) - XB(I))

C     Interpolate new body points at X at each J as we need them, rather than
C     storing them all.  Normally, the starting guess should be close.
C     An oscillatory J search can use the input JT reasonably efficiently.

      J = JT
      IER = 0
      DSQBEST = 9999999.

      IF (CHECK ()) GO TO 999

      JINC = 1

      DO JC = J1B, J2B ! There are probably better ways to allow checking all Js

         J = JT + JINC
         IF (J < J2B) THEN
            IF (CHECK ()) THEN
               JT = J
               GO TO 999
            END IF
         END IF

         J = JT - JINC
         IF (J >= J1B) THEN
            IF (CHECK ()) THEN
               JT = J
               GO TO 999
            END IF
         END IF

         JINC = JINC + 1
      END DO

C     We dropped through without apparent success.  Use the nearest J
C     unless the nearest distance is excessive.

      IF (DSQBEST < ((XB(I2B) - XB(I1B)) * 0.001) ** 2) THEN
         JT = JBEST
      ELSE
         IER = 1
      END IF

C**** write (6, '(a, 2i5, a, 1p, e13.5, /)')
C****>       ' BODYSRCH: IT, JT =', it, jt, '  DSQBEST =', dsqbest

  999 RETURN

C     BODYSRCH internal procedure:

      CONTAINS

!        ------------------------------------------------------------
         LOGICAL FUNCTION CHECK () ! Check the [J,J+1] interval at XT
!        ------------------------------------------------------------

         REAL    A1, A2, ATA, ATB, B1, B2, T, Y1, Y2, Z1, Z2

!        Interpolate body section coordinates at XT:

         Y1 = (ONE - R) * YB(J,I)   + R * YB(J,I+1)
         Z1 = (ONE - R) * ZB(J,I)   + R * ZB(J,I+1)
         Y2 = (ONE - R) * YB(J+1,I) + R * YB(J+1,I+1)
         Z2 = (ONE - R) * ZB(J+1,I) + R * ZB(J+1,I+1)

         A1 = Y2 - Y1
         A2 = Z2 - Z1
         B1 = YT - Y1
         B2 = ZT - Z1

!        The foot of the perpendicular from (YT, ZT) to the line segment
!        is defined by P1 + T * P2 where T = A'B / A'A:

         ATA = A1 * A1 + A2 * A2
         ATB = A1 * B1 + A2 * B2
         T   = ATB / ATA

         CHECK = ZERO <= T .AND. T <= ONE

         IF (.NOT. CHECK) THEN ! Track the shortest distance just in case

            DSQ = B1 * B1 + B2 * B2 - T * ATB

!****       write (6,'(a, 3f11.5, a, 2i5, f11.7, 1p, e13.5)')
!****>         ' x/y/zt: ', xt,yt,zt, '  i,j,t,dsq:', i,j,t,dsq

            IF (DSQ < DSQBEST) THEN
               DSQBEST = DSQ
               JBEST = J
            END IF
!****    else
!****       write (6,'(a, 3f11.5, a, 2i5, f11.7)')
!****>         ' x/y/zt: ', xt,yt,zt, '  i,j,t    :', i,j,t
         END IF

         END FUNCTION CHECK

      END SUBROUTINE BODYSRCH

C***********************************************************************
C
      SUBROUTINE FBARREL (I, MXIGRID, MXJGRID, NJ, METHOD, XF, YF, ZF,
     >                    JL, SNORM, XTEMP, XBASE, YBASE, ZBASE)
C
C     FBARREL imposes a normalized distribution on the given portion of a
C     body section, returning results as the Ith row of the given patch.
C
C     I  = Patch I index where the interpolated barrel is placed
C     NJ = # points used as data for the redistribution
C     JL = # points interpolated at relative locations SNORM(1:JL)
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   I, MXIGRID, MXJGRID, NJ, JL
      CHARACTER, INTENT (IN) ::
     >   METHOD * 1
      REAL, INTENT (IN) ::
     >   XF, YF(NJ), ZF(NJ), SNORM(JL)
      REAL, INTENT (OUT) ::
     >   XTEMP(NJ)
      REAL, INTENT (OUT), DIMENSION (MXIGRID,MXJGRID) ::
     >   XBASE, YBASE, ZBASE

C     Local constants:

      LOGICAL, PARAMETER :: CLOSED = .FALSE., NORMALIZE = .TRUE.

C     Local variables:

      INTEGER
     >   J, JEVAL
      REAL
     >   STOTAL, STEMP(NJ)
      LOGICAL
     >   NEW

C     Storage:

      REAL, PARAMETER :: DERIVS = -999. ! Suppress them

C     Execution:

      CALL CHORDS2D (NJ, YF, ZF, NORMALIZE, STOTAL, STEMP)

      XTEMP = XF ! 1:NJ
      JEVAL = 1
      NEW   = .TRUE.

      DO J = 1, JL - 1

         CALL PLSCRV3D (NJ, XTEMP, YF, ZF, STEMP,
     >                  METHOD, NEW, CLOSED, SNORM(J), JEVAL,
     >                  XBASE(I,J), YBASE(I,J), ZBASE(I,J), DERIVS)
         NEW = .FALSE.
      END DO

      XBASE(I,JL) = XF   ! Exactly
      YBASE(I,JL) = YF(NJ)
      ZBASE(I,JL) = ZF(NJ)

      END SUBROUTINE FBARREL

C***********************************************************************
C
      SUBROUTINE FBASE (MXPATCH, KPATCH, NFSTN, JLBODY, IFUS,
     >                  XF, YF, ZF, LENXB, XBASE, YBASE, ZBASE,
     >                  LASTXB, MOSTXB, IXB, MXSUBI, MXSUBJ,
     >                  NSUBI, NSUBJ, ISUB, JSUB, IMXCOM, JMXCOM, IWRIT)
C
C     FBASE constructs a fuselage base patch as one ~semicircle with a
C     singular point at the mid-point of the last body section at the
C     symmetry plane.  Transfinite interpolation suffices.
C
C     04/12/01  David Saunders  Prompted by the large base area of the
C                               Crew Transfer Vehicle.
C     10/23/02    "     "       Work-space for TFIQ3D was inadequate for
C                               HSCT case; allocate it explicitly.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   MXPATCH               ! Max. # patches expected

      INTEGER, INTENT (INOUT) ::
     >   KPATCH                ! In & out with the next patch # to fill;
                               ! KPATCH - 1 is assumed to be the aft body
      INTEGER, INTENT (IN) ::
     >   NFSTN,                ! # body defining sections
     >   JLBODY,               ! J dimension of regularized body arrays
     >   IFUS                  ! Component # for diagnostics

      REAL, INTENT (IN) ::
     >   XF(NFSTN),            ! Regularized body sections;
     >   YF(JLBODY,NFSTN),     ! YF/ZF(JLBODY-1,NFSTN) is used to handle
     >   ZF(JLBODY,NFSTN)      ! both regular and CTV cases

      INTEGER, INTENT (IN) ::
     >   LENXB                 ! Packed patch space limit

      REAL, INTENT (INOUT), DIMENSION (LENXB) ::
     >   XBASE, YBASE, ZBASE   ! Patch coordinates

      INTEGER, INTENT (INOUT) ::
     >   LASTXB, MOSTXB,       ! Packed patch indices (last, biggest, first)
     >   IXB(MXPATCH)

      INTEGER, INTENT (IN) ::  ! Subpatch info.
     >   MXSUBI,
     >   MXSUBJ

      INTEGER, INTENT (OUT) ::
     >   NSUBI(MXPATCH),
     >   NSUBJ(MXPATCH),
     >   ISUB(0:MXSUBI,MXPATCH),
     >   JSUB(0:MXSUBJ,MXPATCH)

      INTEGER, INTENT (INOUT) ::
     >   IMXCOM(MXPATCH), JMXCOM(MXPATCH) ! Patch dimensions

      INTEGER, INTENT (IN) ::
     >   IWRIT

C     Local constants:

      REAL, PARAMETER :: HALF = 0.5, ZERO = 0.

C     Local variables:

      INTEGER
     >   I, IA, IB, IE, IDIM, J, JDIM, KLAST, NI

      REAL
     >   D1, D2, DY, TOL, YM, YS

      REAL, ALLOCATABLE ::
     >   ARCS(:) ! Room for TFIQ3D is awkward

C     Execution:
C     ----------

      KLAST = KPATCH - 1
      IDIM  = IMXCOM(KLAST)
      JDIM  = JMXCOM(KLAST)
      TOL   = (XF(NFSTN) - XF(1)) * 1.E-5
      IE    = IXB(KLAST) - 1
      IB    = IE + (JDIM / 2 - 1) * IDIM ! ~Water line pt. of last body stn.

C     Avoid a degenerate base patch:

      IE =  IE + IDIM                    ! End of keel point index
      YS =  YBASE(IE)
      YM = (YBASE(LASTXB) + YS) * HALF
      DY =  YBASE(LASTXB) - YS

      IF (ZBASE(IB) < TOL .AND. DY < TOL) GO TO 999

C     Is there enough room for the base patch?

      I  = NFSTN
      D1 = SQRT ((YF(3,I) - YF(2,I))**2 + (ZF(3,I) - ZF(2,I))**2)
      J  = JLBODY - 1  ! Avoid NJBODY/JLBODY confusion
      D2 = SQRT ((YF(J,I) - YF(J-1,I))**2 + (ZF(J,I) - ZF(J-1,I))**2)
      D1 = (D1 + D2) * HALF
      DY = DY * HALF
      NI = NINT (DY / D1)

      CALL XBCHECK (LENXB, MXPATCH, KPATCH, LASTXB, NI, JDIM,
     >              IFUS, IWRIT, MOSTXB)

      IMXCOM(KPATCH) = NI
      JMXCOM(KPATCH) = JDIM
      NSUBI (KPATCH) = 1
      ISUB(1,KPATCH) = NI
      NSUBJ (KPATCH) = NSUBJ(KLAST)
      JSUB(1:NSUBJ(KLAST),KPATCH) = JSUB(1:NSUBJ(KLAST),KLAST)
      DY = DY / REAL (NI - 1)
      IB = LASTXB + 1
      IXB(KPATCH) = IB

C     Set up the base patch edges in-place:
 
      LASTXB = LASTXB + NI * JDIM ! When done
      IA     = LASTXB - NI        ! For base edge above the Z = 0. mid-pt.
 
      DO I = 1, NI
         XBASE(IB) = XF(NFSTN)
         YBASE(IB) = YM - DY * REAL (I - 1)
         ZBASE(IB) = ZERO
         IB        = IB + 1
         IA        = IA + 1
         XBASE(IA) = XF(NFSTN)
         YBASE(IA) = YM + DY * REAL (I - 1)
         ZBASE(IA) = ZERO
      END DO

      IB = IXB(KPATCH) ! Start of base I = 1  edge (singular pt.)
      IA = IB + NI - 1 ! Start of base I = NI edge (base perimeter)

      DO J = 1, JDIM
         XBASE(IB) = XF(NFSTN)
         YBASE(IB) = YM
         ZBASE(IB) = ZERO
         IB        = IB + NI
         XBASE(IA) = XF(NFSTN)
         YBASE(IA) = YBASE(IE)
         ZBASE(IA) = ZBASE(IE)
         IA        = IA + NI
         IE        = IE + IDIM
      END DO

C     Fill the base interior.  TFI2D would be better, but it's not used
C     by other users of the numutil.f collection ...

      IB = IXB(KPATCH)

      ALLOCATE (ARCS(2*(NI + JDIM)))

      CALL TFIQ3D (NI, 1, NI, 1, JDIM, XBASE(IB), YBASE(IB), ZBASE(IB),
     >             ARCS)

      DEALLOCATE (ARCS)

      XBASE(IB:LASTXB) = XF(NFSTN) ! Make sure of it

      KPATCH = KPATCH + 1 ! For the next patch

  999 RETURN

      END SUBROUTINE FBASE

C***********************************************************************
C
      SUBROUTINE FCTV (NWSURF, MXINTR, MXPATCH, KPATCH, NFSTN,
     >                 JLBODY, IFUS, IDEGFUS, EPS, XF, YF, ZF,
     >                 NFIXL, NFIXR, NFIXA, NFIXB, NINTR,
     >                 XINTR, YINTR, ZINTR, UINTR, VINTR, IINTR, JINTR,
     >                 LOSTINTR, LENXB, XBASE, YBASE, ZBASE,
     >                 LASTXB, MOSTXB, IXB, MXSUBI, MXSUBJ,
     >                 NSUBJ, ISUB, JSUB, IMXCOM, JMXCOM, IWRIT, FAIL)
C
C     FCTV imposes regular surface patches on the fuselage of a CTV (Crew
C     Transfer Vehicle)-type configuration with overlapping wing and tail
C     surfaces.
C
C     Simplifying assumptions:
C
C     >  Just two wing-type surfaces, with wing surface 1 lower than
C        tail surface 2, and each having the same number of points per
C        regularized wing section.
C     >  Both intersections employ the whole body (ISTN1/2 = 1:NFSTN)
C        because the whole body parameterization is used here.
C     >  If the nose is a singular line, NJBODY should be odd, but so
C        should JLBODY, so make them the same in all cases (unlike the
C        standard case where JLBODY = NJBODY + 1).
C
C     The driving program should check that these assumptions are met.
C
C     07/13/00  David Saunders  Initial adaptation of FSURF, which is too
C                               long and too general to add this more
C                               specialized case to.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   NWSURF,               ! # wing/nacelle surfaces
     >   MXINTR,               ! Max. # pts. per intersection
     >   MXPATCH               ! Max. # patches expected

      INTEGER, INTENT (INOUT) ::
     >   KPATCH                ! In & out with the next patch # to fill

      INTEGER, INTENT (IN) ::
     >   NFSTN,                ! # body defining sections
     >   JLBODY,               ! # output points per section of (unsplit) nose
     >                         ! panel; JLBODY = NJBODY = # J pts. input per
     >                         ! regularized geometry section (not NJBODY + 1
     >                         ! as in FSURF); JLBODY should be ODD
     >   IFUS,                 ! Component # for diagnostics
     >   IDEGFUS               ! 1 or 3 for [bi]linear or [bi]cubic fus. splines

      REAL, INTENT (IN) ::
     >   EPS,                  ! Tolerance for 2-D searches/interpolations
     >   XF(NFSTN),            ! Regularized body sections;
     >   YF(JLBODY,NFSTN),     ! Y/ZF(1:JLBODY,*) are input (not 1:JLBODY-1)
     >   ZF(JLBODY,NFSTN)

      INTEGER, INTENT (IN), DIMENSION (NWSURF) ::
     >   NFIXL, NFIXR,         ! FIXROOT controls, as for FSURF (wing only)
     >   NFIXA, NFIXB

      INTEGER, INTENT (IN) ::
     >   NINTR(NWSURF)         ! # pts. for each intersection

      REAL, INTENT (IN), DIMENSION (MXINTR,NWSURF) ::
     >   XINTR, YINTR, ZINTR   ! Intersection coordinates

      REAL, INTENT (INOUT), DIMENSION (MXINTR,NWSURF) ::
     >   UINTR, VINTR          ! Corresp. (u,v)s referred to whole body

      INTEGER, INTENT (INOUT), DIMENSION (MXINTR,NWSURF) ::
     >   IINTR, JINTR          ! Body cell indices for intersection points

      LOGICAL, INTENT (IN) ::
     >   LOSTINTR(NWSURF)      ! Controls whether intersection data need to
                               ! be updated here
      INTEGER, INTENT (IN) ::
     >   LENXB                 ! Packed patch space limit

      REAL, INTENT (INOUT), DIMENSION (LENXB) ::
     >   XBASE, YBASE, ZBASE   ! Patch coordinates

      INTEGER, INTENT (INOUT) ::
     >   LASTXB, MOSTXB,       ! Packed patch indices (last, biggest, first)
     >   IXB(MXPATCH)

      INTEGER, INTENT (IN) ::  ! Subpatch info.
     >   MXSUBI, MXSUBJ

      INTEGER, INTENT (OUT) ::
     >   NSUBJ(MXPATCH),
     >   ISUB(0:MXSUBI,MXPATCH),
     >   JSUB(0:MXSUBJ,MXPATCH)

      INTEGER, INTENT (INOUT) ::
     >   IMXCOM(MXPATCH), JMXCOM(MXPATCH) ! Patch dimensions

      INTEGER, INTENT (IN) ::
     >   IWRIT

      LOGICAL, INTENT (OUT) ::
     >   FAIL

C     Local constants:

      REAL,      PARAMETER :: HALF = 0.5, ONE = 1., ZERO = 0.
      CHARACTER, PARAMETER :: METHOD = 'B' ! For CHANGEN/PLSCRV3D

C     Local variables:

      INTEGER
     >   I, I1, I2, IB, IE, IER, ILE, IMAT, IX, IXSHIFT,
     >   J, J1, J2, JB, JE, JM, JMAT,
     >   K, K1, K2, K3, K4, K5, KASE, L, LUNERR, M,
     >   N, NFIX, NI, NIN, NIT, NITOT, NJA, NJB, NJM

      REAL
     >   DRANGE, P, Q, R, SMATRIX(4,4,3), UCUT, VCUT, XYZCUT(3),
     >   XLINE(NFSTN),  YLINE(NFSTN),  ZLINE(NFSTN),
     >   XLIN1(JLBODY), YLIN1(JLBODY), ZLIN1(JLBODY),
     >   XLIN2(JLBODY), YLIN2(JLBODY), ZLIN2(JLBODY)

      REAL, DIMENSION (JLBODY,NFSTN) ::
     >   UBODY, VBODY, XBODY

      REAL, ALLOCATABLE, DIMENSION (:,:) ::
     >   UEVAL, VEVAL

      LOGICAL
     >   BICUBIC

C     Storage:

      REAL, DIMENSION(4), SAVE ::
     >   UVRANGE

      DATA
     >   UVRANGE /ZERO, ONE, ZERO, ONE/ ! See PLBICUT

C     Execution:
C     ----------

      K1 = KPATCH ! First patch to fill
      K2 = K1 + 1
      K3 = K1 + 2
      K4 = K1 + 3
      K5 = K1 + 4

C     Parameterize the body geometry:

      UBODY(1,1) = ZERO ! -999. suppresses normalization

      DO I = 1, NFSTN
         XBODY(1:JLBODY,I) = XF(I)
      END DO

      CALL PARAM2D (JLBODY, NFSTN, 1, JLBODY, 1, NFSTN, ! Note U/V order
     >              XBODY, YF, ZF, VBODY, UBODY)        ! U is streamwise

      BICUBIC = IDEGFUS == 3
      LUNERR  = -ABS (IWRIT)
      DRANGE  = XF(NFSTN) - XF(1)
      KASE    = 1 ! Find Zcut given X- & YCUT
      J = JLBODY
      N = NINTR(1) ! (= NINTR(2))

      DO L = 1, 2

         IF (LOSTINTR(L)) THEN

            UCUT = UINTR(2,L)
            VCUT = VINTR(2,L)
            IE   = IINTR(2,L)
            JE   = JINTR(2,L)
            IMAT = 0
            JMAT = 0

            DO I = 2, N - 1

               XYZCUT(1) = XINTR(I,L)
               XYZCUT(2) = YINTR(I,L)
               XYZCUT(3) = ZINTR(I,L)

               IF (BICUBIC) THEN

                  CALL PLBICUT (KASE, J, NFSTN, 1, J, 1, NFSTN,
     >                          XBODY, YF, ZF, VBODY, UBODY, JE, IE,
     >                          JMAT, IMAT, SMATRIX,
     >                          EPS, DRANGE, XYZCUT, UVRANGE,
     >                          VCUT, UCUT, Q, P, LUNERR, IER)
               ELSE

                  CALL BILINCUT (KASE, J, NFSTN, 1, J, 1, NFSTN,
     >                           XBODY, YF, ZF, VBODY, UBODY, JE, IE,
     >                           EPS, DRANGE, XYZCUT, UVRANGE,
     >                           VCUT, UCUT, Q, P, LUNERR, IER)
               END IF

               IF (IER > 2) THEN
                  IF (IWRIT > 0) THEN
                     WRITE (IWRIT,
     >                  '(2(/, A, I2), /, A, I4, /, A, 1P, 3E15.6)')
     >                  ' FCTV: PLBICUT/BILINCUT error =', IER,
     >                  ' Wing/fin surface =', L,
     >                  ' Intersection pt. =', I,
     >                  ' Corresp.(x,y,z)  =', XYZCUT
                  END IF

                  CALL SYNCH
                  STOP

               END IF

               UINTR(I,L) = UCUT
               VINTR(I,L) = VCUT
               IINTR(I,L) = IE
               JINTR(I,L) = JE

            END DO

         END IF

         KASE = 3 ! Find YCUT given ZCUT and XCUT (for the fin, if redone)

      END DO

C     Establish the dimensions of all subpatches:

      ILE = (N + 1) / 2
      NIN = MAX (IINTR(ILE,1), IINTR(ILE,2)) + 1 ! Length of nose panel

      I1  = MIN (IINTR(1,1), IINTR(N,1))
      I2  = MIN (IINTR(1,2), IINTR(N,2))
      NIT = NFSTN - MIN (I1, I2) + 1             ! ... and tail panel

      IMXCOM(K1)    = NIN ! Nose
      IMXCOM(K2:K4) = ILE ! Panels below/between/above wing & tail surfaces
      IMXCOM(K5)    = NIT ! Aft body

      ISUB(1,K1)    = NIN ! No subpatches in the I direction
      ISUB(1,K2:K4) = ILE
      ISUB(1,K5)    = NIT
      NITOT         = NIN + ILE + NIT - 2

      NJB = JINTR(ILE,1) + 1      ! # js in patch below wing
      NJA = JLBODY - JINTR(ILE,2) ! # js in patch above tail
      J2  = JLBODY - NJA + 1      ! Nose patch j index at tail LE
      NJM = J2     - NJB + 1      ! # js in middle patch

      JMXCOM(K1)    = JLBODY ! Nose
      JMXCOM(K2)    = NJB    ! Below wing
      JMXCOM(K3)    = NJM    ! Between wings
      JMXCOM(K4)    = NJA    ! Above tail
      JMXCOM(K5)    = JLBODY ! Aft body

      NSUBJ(K1)     = 3
      NSUBJ(K5)     = 3

      JSUB(1,K1)    = NJB
      JSUB(2,K1)    = J2
      JSUB(3,K1)    = JLBODY
      JSUB(1,K2:K4) = JMXCOM(K2:K4)
      JSUB(1:3, K5) = JSUB(1:3, K1)

C     Check that there is room for all fuselage panels:

      IX = LASTXB + 1
      I1 = IX
      M  = 0

      DO L = K1, K5
         NI = IMXCOM(L) * JMXCOM(L)
         M  = NI + M
         IXB(L) = I1
         I1 = I1 + NI
      END DO

      CALL XBCHECK (LENXB, MXPATCH, K5, LASTXB, M, 1,
     >              IFUS, IWRIT, MOSTXB)

      LASTXB = LASTXB + M ! Upon return from FCTV

C     Strategy for the subpatches:
C
C     >  Work from keel to crown, forward to aft.
C     >  Exception: for subpatches between wing & tail, reduce likely
C        trouble with edge (u,v)s by using the TFI result for the full
C        length of the body.
C     >  Establish the corner points in (u,v) space.
C     >  Generate appropriate edges  in (u,v) space.
C     >  Generate interior (u,v)s via TFI.
C     >  Interpolate interior and edge (x,y,z)s (bilinear or bicubic).

C     Allocate (u,v)s for interpolating all patches (edges overlapping):

      ALLOCATE (UEVAL(NITOT,JLBODY), VEVAL(NITOT,JLBODY))

      DO J = 1, JLBODY
         UEVAL(1,J)      = ZERO
         UEVAL(NITOT,J)  = ONE
      END DO

      DO I = 1, NITOT
         VEVAL(I,1)      = ZERO
         VEVAL(I,JLBODY) = ONE
      END DO

C     Deal with the nose point or line for three subpatches:
C     ------------------------------------------------------

      I = IX

      IF (ZF(2,1) == ZF(1,1)) THEN ! Point

         DO J = 1, JLBODY
            XBASE(I) = XF(1)
            YBASE(I) = YF(1,1)
            ZBASE(I) = ZF(1,1)
            I = I + NIN
            VEVAL(1,J) = VBODY(J,1)
         END DO

      ELSE ! Singular line needs redistribution

         JM = (JLBODY + 1) / 2 ! JLBODY is odd to capture the water-line

         DO J = 1, JLBODY
            XLIN1(J) = XF(1)
            YLIN1(J) = YF(J,1)
            ZLIN1(J) = ZF(J,1)
         END DO

         XLIN2(1)   = XLIN1(1)
         YLIN2(1)   = YLIN1(1)
         ZLIN2(1)   = ZLIN1(1)
         XLIN2(NJB) = XLIN1(JM)
         YLIN2(NJB) = YLIN1(JM)
         ZLIN2(NJB) = ZLIN1(JM)

         XLIN2(J2)     = XLIN1(JM)
         YLIN2(J2)     = YLIN1(JM)
         ZLIN2(J2)     = ZLIN1(JM)
         XLIN2(JLBODY) = XLIN1(JLBODY)
         YLIN2(JLBODY) = YLIN1(JLBODY)
         ZLIN2(JLBODY) = ZLIN1(JLBODY)

         CALL CHANGEN (1,  JM,     XLIN1, YLIN1, ZLIN1,
     >                 1,  NJB,    XLIN2, YLIN2, ZLIN2, METHOD)

         CALL CHANGEN (JM, JLBODY, XLIN1, YLIN1, ZLIN1,
     >                 J2, JLBODY, XLIN2, YLIN2, ZLIN2, METHOD)

         R = VBODY(JM,1) / REAL (NJB - 1)
         DO J = 1, NJB
            XBASE(I) = XLIN2(J)
            YBASE(I) = YLIN2(J)
            ZBASE(I) = ZLIN2(J)
            I = I + NIN
            VEVAL(1,J) = R * REAL (J - 1)
         END DO

         DO J = NJB + 1, J2
            XBASE(I) = XLIN2(NJB)
            YBASE(I) = YLIN2(NJB)
            ZBASE(I) = ZLIN2(NJB)
            I = I + NIN
            VEVAL(1,J) = VBODY(JM,1)
         END DO

         R = (ONE - VBODY(JM,1)) / REAL (NJA - 1)
         DO J = J2 + 1, JLBODY
            XBASE(I) = XLIN2(J)
            YBASE(I) = YLIN2(J)
            ZBASE(I) = ZLIN2(J)
            I = I + NIN
            VEVAL(1,J) = VBODY(JM,1) + R * REAL (J - J2)
         END DO

      END IF

C     First subpatch (on the nose at the keel):
C     -----------------------------------------

C     Upper right corner (wing LE):

      UEVAL(NIN,NJB) = UINTR(ILE,1)
      VEVAL(NIN,NJB) = VINTR(ILE,1)

C     Right-hand edge (u,v)s:

      NI =  IINTR(ILE,1) + 1
      J  =  JINTR(ILE,1)
      R  =  VINTR(ILE,1) / REAL (NJB - 1)
      P  = (UINTR(ILE,1) - UBODY(J,NI-1)) /
     >     (UBODY(J,NI)  - UBODY(J,NI-1))

      DO J = 1, NJB - 1
         UEVAL(NIN,J) = (ONE - P) * UBODY(J,NI-1) + P * UBODY(J,NI)
         VEVAL(NIN,J) = R * REAL (J - 1)
      END DO

C     Keel-line (u,v)s:

      ZLINE = ZERO ! For all I, to be safe

      IF (NIN == NI) THEN

         DO I = 2, NIN - 1
            UEVAL(I,1) = UBODY(1,I)
         END DO

      ELSE ! 3-space CHANGEN can be used for a 2-space (u,v) line

         DO I = 1, NIN - 1
            XLINE(I) = UBODY(1,I)
            YLINE(I) = ZERO
         END DO

         XLINE(NI) = UEVAL(NIN,1)
         YLINE(NI) = ZERO

         CALL CHANGEN (1, NI, XLINE, YLINE, ZLINE, 1, NIN,
     >                 UEVAL, VEVAL, ZLINE, METHOD)
      END IF

C     Upper edge (u,v)s of lower nose subpatch:

      R = ONE / (XINTR(ILE,1) - XF(1)) ! A straight (u,v) line suffices

      IF (NIN == NI) THEN

         DO I = 2, NIN - 1
            P = (XF(I) - XF(1)) * R
            UEVAL(I,NJB) = P * UEVAL(NIN,NJB)
            VEVAL(I,NJB) = P * VEVAL(NIN,NJB) + (ONE - P) * VEVAL(1,NJB)
         END DO

      ELSE

         DO I = 1, NI - 1
            P = (XF(I) - XF(1)) * R
            XLINE(I) = P * UEVAL(NIN,NJB)
            YLINE(I) = P * VEVAL(NIN,NJB) + (ONE - P) * VEVAL(1,NJB)
         END DO

         XLINE(NI) = UEVAL(NIN,NJB)
         YLINE(NI) = VEVAL(NIN,NJB)

         CALL CHANGEN (1, NI, XLINE, YLINE, ZLINE, 1, NIN,
     >                 UEVAL(1,NJB), VEVAL(1,NJB), ZLINE, METHOD)
      END IF


C     Second subpatch on the nose:
C     ----------------------------

      UEVAL(NIN,J2) = UINTR(ILE,2) ! Upper right corner (tail LE)
      VEVAL(NIN,J2) = VINTR(ILE,2)

C     Upper edge (u,v)s:

      NI = IINTR(ILE,2) + 1
      R  = ONE / (XINTR(ILE,2) - XF(1))

      IF (NIN == NI) THEN

         DO I = 2, NIN - 1
            P = (XF(I) - XF(1)) * R
            UEVAL(I,J2) = P * UEVAL(NIN,J2) 
            VEVAL(I,J2) = P * VEVAL(NIN,J2) + (ONE - P) * VEVAL(1,J2)
         END DO

      ELSE

         DO I = 1, NI - 1
            P = (XF(I) - XF(1)) * R
            XLINE(I) = P * UEVAL(NIN,J2) 
            YLINE(I) = P * VEVAL(NIN,J2) + (ONE - P) * VEVAL(1,J2)
         END DO

         XLINE(NI) = UEVAL(NIN,J2)
         YLINE(NI) = VEVAL(NIN,J2)

         CALL CHANGEN (1, NI, XLINE, YLINE, ZLINE, 1, NIN,
     >                 UEVAL(1,J2), VEVAL(1,J2), ZLINE, METHOD)
      END IF

C     Linear (u,v)s connecting leading edges may give trouble at wing LE,
C     so fill the (u,v)s for the full body length later via a single TFI.

C     Upper nose subpatch:
C     --------------------

C     Right-hand edge (u,v)s:

      NI =  IINTR(ILE,2) + 1
      J  =  JINTR(ILE,2)
      P  = (UINTR(ILE,2) - UBODY(J,NI-1)) /
     >     (UBODY(J,NI)  - UBODY(J,NI-1))
      R  = (ONE - VINTR(ILE,2)) / REAL (NJA - 1)

      DO J = J2 + 1, JLBODY
         UEVAL(NIN,J) = (ONE - P) * UBODY(J,NI-1) + P * UBODY(J,NI)
         VEVAL(NIN,J) = VINTR(ILE,2) + R * REAL (J - J2)
      END DO

C     Crown-line (u,v)s:

      IF (NIN == NI) THEN

         DO I = 2, NIN - 1
            UEVAL(I,JLBODY) = UBODY(JLBODY,I)
         END DO

      ELSE

         J = JLBODY
         DO I = 1, NIN - 1
            XLINE(I) = UBODY(J,I)
            YLINE(I) = ONE
         END DO

         XLINE(NI) = UEVAL(NIN,J)
         YLINE(NI) = ONE

         CALL CHANGEN (1, NI, XLINE, YLINE, ZLINE, 1, NIN,
     >                 UEVAL(1,J), VEVAL(1,J), ZLINE, METHOD)
      END IF


C     Panel under the wing:
C     ---------------------

C     Close the trailing edge in (u,v) space:

      I2 = NIN + ILE - 1
      UEVAL(I2,NJB) = (UINTR(1,1) + UINTR(N,1)) * HALF
      VEVAL(I2,NJB) = (VINTR(1,1) + VINTR(N,1)) * HALF

C     Right edge (u,v)s:

      I =  IINTR(1,1)
      J =  JINTR(1,1)
      R =  VEVAL(I2,NJB) / REAL (NJB - 1)
      P = (UEVAL(I2,NJB) - UBODY(J,I)) /
     >    (UBODY(J,I+1)  - UBODY(J,I))

      DO J = 1, NJB - 1
         UEVAL(I2,J) = (ONE - P) * UBODY(J,I) + P * UBODY(J,I+1)
         VEVAL(I2,J) = R * REAL (J - 1)
      END DO

C     Keel (u,v)s + upper patch edge (= lower half of intersection).
C     Use "arc length" along the intersection (u,v)s in case U(LE) is not
C     at the same index as X(LE).

      IX = IXB(K2) ! Start of ample work-space for arc lengths

      CALL CHORDS2D (ILE, UINTR(1,1), VINTR(1,1), .TRUE., P, XBASE(IX))

      P  = UEVAL(NIN,1)
      R  = UEVAL(I2,1) - P
      I1 = IX + ILE - 1
      J  = ILE

      DO I = NIN + 1, I2 - 1
         I1 = I1 - 1
         UEVAL(I,1) = P + R * (ONE - XBASE(I1))
         J  = J - 1
         UEVAL(I,NJB) = UINTR(J,1)
         VEVAL(I,NJB) = VINTR(J,1)
      END DO

C     Fill the interior of patch K2 now so we can reuse its upper edge,
C     U/VEVAL(*,NJB), for the lower edge of patch K3.

      CALL TFI2D (NITOT, NIN, I2, 1, NJB, UEVAL, VEVAL, XBASE(IX))


C     Do the panel between the wing & tail later.


C     Panel above the tail:
C     ---------------------

C     Close the trailing edge:

      UEVAL(I2,J2) = (UINTR(1,2) + UINTR(N,2)) * HALF
      VEVAL(I2,J2) = (VINTR(1,2) + VINTR(N,2)) * HALF

C     Right edge (uv)s:

      I =  IINTR(1,2)
      J =  JINTR(1,2)
      R = (ONE - VEVAL(I2,J2)) / REAL (NJA - 1)
      P = (UEVAL(I2,J2) - UBODY(J,I)) /
     >    (UBODY(J,I+1) - UBODY(J,I))

      DO J = J2 + 1, JLBODY
         UEVAL(I2,J) = (ONE - P) * UBODY(J,I) + P * UBODY(J,I+1)
         VEVAL(I2,J) = VEVAL(I2,J2) + R * REAL (J - J2)
      END DO

C     Crown-line (u,v)s, + lower patch edge (= upper half of intersection).

      CALL CHORDS2D (ILE, UINTR(ILE,2), VINTR(ILE,2), .TRUE., P,
     >               XBASE(IX))

      P  = UINTR(ILE,2)
      Q  = UEVAL(NIN,JLBODY)
      R  = UEVAL(I2,JLBODY) - Q
      I1 = IX
      J  = ILE

      DO I = NIN + 1, I2 - 1
         J = J + 1
         UEVAL(I,J2) = UINTR(J,2)
         VEVAL(I,J2) = VINTR(J,2)
         I1 = I1 + 1
         UEVAL(I,JLBODY) = Q + R * XBASE(I1)
      END DO

C     Interior of K4:

      CALL TFI2D (NITOT, NIN, I2, J2, JLBODY, UEVAL, VEVAL, XBASE(IX))


C     Aft body panel in three subpatches:
C     -----------------------------------

C     Upper right corner:

      I1 = IINTR(1,1)
      J1 = JINTR(1,1)
      Q  = (VEVAL(I2,NJB)  - VBODY(J1,I1)) /
     >     (VBODY(J1+1,I1) - VBODY(J1,I1))
      P  = (ONE - Q) * VBODY(J1,  NFSTN) +   ! P = VEVAL(NITOT,NJB)
     >            Q  * VBODY(J1+1,NFSTN)

C     Right-hand edge (u,v)s:

      R = P / REAL (NJB - 1)

      DO J = 2, NJB
         VEVAL(NITOT,J) = R * REAL (J - 1)
      END DO

C     Keel-line edge (u,v)s, and edge aft of wing TE:

      IF (NIT == NFSTN - I1 + 1) THEN

         DO I = I2 + 1, NITOT - 1
            I1 = I1 + 1
            UEVAL(I,1)   = UBODY(1,I1)
            UEVAL(I,NJB) = (ONE - Q) * UBODY(J1,I1) + Q * UBODY(J1+1,I1)
            VEVAL(I,NJB) = (ONE - Q) * VBODY(J1,I1) + Q * VBODY(J1+1,I1)
         END DO

      ELSE ! Adjust the number of points

         XLINE(I1) = UEVAL(I2,1) ! Keel line first
         YLINE(I1) = ZERO

         DO I = I1 + 1, NFSTN
            XLINE(I) = UBODY(1,I)
            YLINE(I) = ZERO
         END DO

         CALL CHANGEN (I1, NFSTN, XLINE, YLINE, ZLINE,
     >                 I2, NITOT, UEVAL, VEVAL, ZLINE, METHOD)

         XLINE(I1) = UEVAL(I2,NJB) ! Now for wing TE --> end of body
         YLINE(I1) = VEVAL(I2,NJB)

         DO I = I1 + 1, NFSTN
            XLINE(I) = (ONE - Q) * UBODY(J1,I) + Q * UBODY(J1+1,I)
            YLINE(I) = (ONE - Q) * VBODY(J1,I) + Q * VBODY(J1+1,I)
         END DO

         CALL CHANGEN (I1, NFSTN, XLINE, YLINE, ZLINE,
     >                 I2, NITOT, UEVAL(1,NJB), VEVAL(1,NJB), ZLINE,
     >                 METHOD)
      END IF

C     Likewise for the upper aft-body subpatch edges:

C     Lower right corner:

      I1 = IINTR(1,2)
      J1 = JINTR(1,2)
      Q  = (VEVAL(I2,J2)   - VBODY(J1,I1)) /
     >     (VBODY(J1+1,I1) - VBODY(J1,I1))
      P  = (ONE - Q) * VBODY(J1,  NFSTN) +   ! P = VEVAL(NITOT,J2)
     >            Q  * VBODY(J1+1,NFSTN)

C     Right-hand edge (u,v)s:

      R = (ONE - P) / REAL (NJA - 1)

      DO J = J2, JLBODY - 1
         VEVAL(NITOT,J) = P + R * REAL (J - J2)
      END DO

C     Crown-line edge (u,v)s, and edge aft of tail TE:

      IF (NIT == NFSTN - I1 + 1) THEN

         DO I = I2 + 1, NITOT - 1
            I1 = I1 + 1
            UEVAL(I,J2) = (ONE - Q) * UBODY(J1,I1) + Q * UBODY(J1+1,I1)
            VEVAL(I,J2) = (ONE - Q) * VBODY(J1,I1) + Q * VBODY(J1+1,I1)
            UEVAL(I,JLBODY) = UBODY(JLBODY,I1)
         END DO

      ELSE ! Adjust the number of points

         XLINE(I1) = UEVAL(I2,J2) ! Tail TE --> end of body
         YLINE(I1) = VEVAL(I2,J2)

         DO I = I1 + 1, NFSTN
            XLINE(I) = (ONE - Q) * UBODY(J1,I) + Q * UBODY(J1+1,I)
            YLINE(I) = (ONE - Q) * VBODY(J1,I) + Q * VBODY(J1+1,I)
         END DO

         CALL CHANGEN (I1, NFSTN, XLINE, YLINE, ZLINE,
     >                 I2, NITOT, UEVAL(1,J2), VEVAL(1,J2), ZLINE,
     >                 METHOD)

         XLINE(I1) = UEVAL(I2,JLBODY) ! Now for the crown line
         YLINE(I1) = ONE

         DO I = I1 + 1, NFSTN
            XLINE(I) = UBODY(JLBODY,I)
            YLINE(I) = ONE
         END DO

         CALL CHANGEN (I1, NFSTN, XLINE, YLINE, ZLINE,
     >                 I2, NITOT, UEVAL(1,JLBODY), VEVAL(1,JLBODY),
     >                 ZLINE, METHOD)
      END IF

C     Aft edge of middle aft-body subpatch:

      P =  VEVAL(NITOT,NJB)
      R = (VEVAL(NITOT,J2) - P) / REAL (NJM - 1)

      DO J = NJB + 1, J2 - 1
         VEVAL(NITOT,J) = P + R * REAL (J - NJB)
      END DO


C     Fill the interior (u,v)s of the nose & aft-body subpatches at keel/crown:
C     -------------------------------------------------------------------------

      DO KASE = 1, 5, 4

         SELECT CASE (KASE)

            CASE (1)
               K  = K1
               I1 = 1
               L  = NIN
            CASE (5)
               K  = K5
               I1 = I2
               L  = NITOT

         END SELECT

         DO J = 1, 3, 2

            J1 = JSUB(J-1,K)
            M  = JSUB(J,K)

            CALL TFI2D (NITOT, I1, L, J1, M, UEVAL, VEVAL, XBASE(IX))

         END DO

      END DO

C     Middle region is filled as one large subpatch.  Transcribe the
C     upper wing and lower tail intersection curves as edges:

      DO I = 1, ILE - 2
         UEVAL(NIN + I,NJB) = UINTR(ILE + I,1)
         VEVAL(NIN + I,NJB) = VINTR(ILE + I,1)
         UEVAL(I2 - I, J2)  = UINTR(I + 1,  2)
         VEVAL(I2 - I, J2)  = VINTR(I + 1,  2)
      END DO

      CALL TFI2D (NITOT, 1, NITOT, NJB, J2, UEVAL, VEVAL, XBASE(IX))

!!!      write (60) 3
!!!      write (60) nitot, jlbody, n, 1, n, 1
!!!      write (60) ueval, veval
!!!      write (60) uintr(1:n,1), vintr(1:n,1)
!!!      write (60) uintr(1:n,2), vintr(1:n,2)

C     Parametric bilinear or bicubic interpolation of the body panel (x,y,z)s:
C     ------------------------------------------------------------------------

      K = K1
      IMAT = 0
      JMAT = 0
      IXSHIFT = 0

      DO KASE = 1, 5

         SELECT CASE (KASE)

            CASE (1)
               I1 = 2
               L  = NIN
               J1 = 1
               M  = JLBODY
               IB = 2
               JB = 1
            CASE (2)
               I1 = NIN
               L  = I2
!              J1 = 1 (still)
               M  = NJB - 1
               IB = IINTR(ILE,1)
               JB = 1
            CASE (3)
               J1 = NJB + 1
               M  = J2  - 1
               IB = IINTR(ILE,1)
               JB = JINTR(ILE,1)
               IXSHIFT = ILE   ! Since J is not starting at the edge
            CASE (4)
               J1 = J2  + 1
               M  = JLBODY
               IB = IINTR(ILE,2)
               JB = JINTR(ILE,2)
            CASE (5)
               I1 = I2
               L  = NITOT
               J1 = 1
               IB = IINTR(1,1)
               JB = 1
               IXSHIFT = 0

         END SELECT

         IX = IXB(K) + IXSHIFT

         DO J = J1, M

            IF (K == K1) IX = IX + 1 ! Awkwardness from omitting nose pts.

            IE = IB
            JE = JB

            DO I = I1, L

               IF (BICUBIC) THEN

                  CALL PLBICUBE (0, JLBODY, NFSTN, 1, JLBODY, 1, NFSTN,
     >                           XBODY, YF, ZF, VBODY, UBODY,
     >                           VEVAL(I,J), UEVAL(I,J), JE, IE, 
     >                           JMAT, IMAT, SMATRIX, EPS, Q, P,
     >                           XBASE(IX), YBASE(IX), ZBASE(IX),
     >                           XLINE, YLINE, IER)
               ELSE

                  CALL PBILINT  (0, JLBODY, NFSTN, 1, JLBODY, 1, NFSTN,
     >                           XBODY, YF, ZF, VBODY, UBODY,
     >                           VEVAL(I,J), UEVAL(I,J), JE, IE,
     >                           EPS, Q, P,
     >                           XBASE(IX), YBASE(IX), ZBASE(IX),
     >                           XLINE, YLINE, IER)
               END IF

               IF (IER /= 0) THEN
                  IF (IWRIT > 0) THEN
                    WRITE (IWRIT,'(/,A,A,I2,A,4I5,/,A,2F13.5,A,3F13.5)')
     >               ' FCTV: Marginal 2D interpolation. ',
     >               ' IER:', IER, '; I, J, IE, JE:', I, J, IE, JE,
     >               ' U, V:', UEVAL(I,J), VEVAL(I,J),
     >               '   X, Y, Z used:', XBASE(IX), YBASE(IX), ZBASE(IX)
                  END IF
               END IF

               IF (I == I1) THEN
                  IB = IE ! Ready for the next J
                  JB = JE
               END IF

               IX = IX + 1

            END DO

         END DO

         K = K + 1

      END DO ! Next panel


      DEALLOCATE (UEVAL, VEVAL)


C     Transcribe intersections (x,y,z)s and ensure common edge points:

      IX = IXB(K1) + NIN * NJB - 1   ! Wing LE pt. of nose panel
      XBASE(IX) = XINTR(ILE,1)
      YBASE(IX) = YINTR(ILE,1)
      ZBASE(IX) = ZINTR(ILE,1)

      IX = IXB(K1) + NIN * J2  - 1   ! Tail LE pt. of nose panel
      XBASE(IX) = XINTR(ILE,2)
      YBASE(IX) = YINTR(ILE,2)
      ZBASE(IX) = ZINTR(ILE,2)

      IX = IXB(K5) + NIT * (NJB - 1) ! Wing TE pt. from aft-body panel
      I1 = IXB(K3) - 1               ! Corresp. pt. in under-wing panel
      XBASE(I1) = XBASE(IX)
      YBASE(I1) = YBASE(IX)
      ZBASE(I1) = ZBASE(IX)

      I1 = IXB(K3) + ILE - 1         ! Corresp. pt. of panel between wing & tail
      XBASE(I1) = XBASE(IX)
      YBASE(I1) = YBASE(IX)
      ZBASE(I1) = ZBASE(IX)

      IX = IXB(K5) + NIT * (J2 - 1)  ! Tail TE pt. from aft-body panel
      I1 = IXB(K4) - 1               ! Corresp. pt. in panel between wing & tail
      XBASE(I1) = XBASE(IX)
      YBASE(I1) = YBASE(IX)
      ZBASE(I1) = ZBASE(IX)

      I1 = IXB(K4) + ILE - 1         ! Corresp. pt. of panel above tail
      XBASE(I1) = XBASE(IX)
      YBASE(I1) = YBASE(IX)
      ZBASE(I1) = ZBASE(IX)

      IX = IXB(K3) - 1               ! Lower wing intersection (except TE)
      DO I = 2, ILE
         IX = IX - 1
         XBASE(IX) = XINTR(I,1)
         YBASE(IX) = YINTR(I,1)
         ZBASE(IX) = ZINTR(I,1)
      END DO

      IX = IXB(K3)                   ! Upper wing intersection (except TE)
      DO I = ILE, N - 1
         XBASE(IX) = XINTR(I,1)
         YBASE(IX) = YINTR(I,1)
         ZBASE(IX) = ZINTR(I,1)
         IX = IX + 1
      END DO

      IX = IXB(K4) - 1               ! Lower tail intersection (except TE)
      DO I = 2, ILE
         IX = IX - 1
         XBASE(IX) = XINTR(I,2)
         YBASE(IX) = YINTR(I,2)
         ZBASE(IX) = ZINTR(I,2)
      END DO

      IX = IXB(K4)                   ! Upper tail intersection (except TE)
      DO I = ILE, N - 1
         XBASE(IX) = XINTR(I,2)
         YBASE(IX) = YINTR(I,2)
         ZBASE(IX) = ZINTR(I,2)
         IX = IX + 1
      END DO

      
C     Fix up extreme skewing at blunt wing leading edge?

      NFIX = NFIXL(1) + NFIXR(1) + 1

      IF (NFIX > 2) THEN ! Patches K1, K1+1, K1+2 are partially revised

         CALL FIXROOT (MXPATCH, LENXB, IXB, IMXCOM, JMXCOM, K1,
     >                 JLBODY, NFIXL(1), NFIXR(1), NFIX, NFIXA(1), 
     >                 NFIXB(1), EPS, XBASE, YBASE, ZBASE)
      END IF


      KPATCH = K5 + 1 ! Next available patch
      FAIL   = .FALSE.

      END SUBROUTINE FCTV

C***********************************************************************
C
      SUBROUTINE FCTV2 (NWSURF, MXINTR, MXPATCH, KPATCH, NFSTN,
     >                  JLBODY, IFUS, IDEGFUS, EPS, XF, YF, ZF,
     >                  NFIXL, NFIXR, NFIXA, NFIXB, NINTR,
     >                  XINTR, YINTR, ZINTR, UINTR, VINTR, IINTR, JINTR,
     >                  LOSTINTR, LENXB, XBASE, YBASE, ZBASE,
     >                  LASTXB, MOSTXB, IXB, MXSUBI, MXSUBJ,
     >                  NSUBI, NSUBJ, ISUB, JSUB, IMXCOM, JMXCOM,
     >                  IWRIT, FAIL)
C
C     FCTV2 imposes regular surface patches on the fuselage of a CTV (Crew
C     Transfer Vehicle)-type configuration with overlapping wing and tail
C     surfaces.  This version is more specialized than FCTV, but overcomes
C     weaknesses near the wing leading edge if the tail is small/well aft.
C
C     Simplifying assumptions:
C
C     >  Just two wing-type surfaces, with wing surface 1 lower than
C        tail surface 2.
C        ILWING(1) >> ILWING(2).  E.g., double.  (ILWING = NINTR here.)
C     >  Both intersections employ the whole body (ISTN1/2 = 1:NFSTN)
C        because the whole body parameterization is used here.
C     >  If the nose is a singular line, NJBODY should be odd, but so
C        should JLBODY, so make them the same in all cases (unlike the
C        standard case where JLBODY = NJBODY + 1).
C
C     The driving program should check that these assumptions are met.
C
C     07/13/00  David Saunders  Initial FCTV adaptation of FSURF.
C     07/27/00    "      "      FCTV2 is more specialized (requiring
C                               splitting of the wing patches also).
C     08/02/00    "      "      Very blunt leading edges force careful
C                               capture of the true intersection leading
C                               edges, in conjunction with higher level
C                               redistribution of wing sections via
C                               input LCRANK > 0, which in turn means
C                               having to update the intersection data.
C     08/03/00    "      "      Installed FIXROOT option.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   NWSURF,               ! # wing/nacelle surfaces
     >   MXINTR,               ! Max. # pts. per intersection
     >   MXPATCH               ! Max. # patches expected

      INTEGER, INTENT (INOUT) ::
     >   KPATCH                ! In & out with the next patch # to fill

      INTEGER, INTENT (IN) ::
     >   NFSTN,                ! # body defining sections
     >   JLBODY,               ! # output points per section of (unsplit) nose
     >                         ! panel; JLBODY = NJBODY = # J pts. input per
     >                         ! regularized geometry section (not NJBODY + 1
     >                         ! as in FSURF); JLBODY should be ODD
     >   IFUS,                 ! Component # for diagnostics
     >   IDEGFUS               ! 1 or 3 for [bi]linear or [bi]cubic fus. splines

      REAL, INTENT (IN) ::
     >   EPS,                  ! Tolerance for 2-D searches/interpolations
     >   XF(NFSTN),            ! Regularized body sections;
     >   YF(JLBODY,NFSTN),     ! Y/ZF(1:JLBODY,*) are input (not 1:JLBODY-1)
     >   ZF(JLBODY,NFSTN)

      INTEGER, INTENT (IN), DIMENSION (NWSURF) ::
     >   NFIXL, NFIXR,         ! FIXROOT controls, as for FSURF (wing only)
     >   NFIXA, NFIXB

      INTEGER, INTENT (IN) ::
     >   NINTR(NWSURF)         ! # pts. for each intersection; (1) >> (2)

      REAL, INTENT (IN), DIMENSION (MXINTR,NWSURF) ::
     >   XINTR, YINTR, ZINTR   ! Intersection coordinates

      REAL, INTENT (INOUT), DIMENSION (MXINTR,NWSURF) ::
     >   UINTR, VINTR          ! Corresp. (u,v)s referred to whole body;
                               ! modified here if LOSTINTR() = T

      INTEGER, INTENT (INOUT), DIMENSION (MXINTR,NWSURF) ::
     >   IINTR, JINTR          ! Body cell indices for intersection points

      LOGICAL, INTENT (IN) ::
     >   LOSTINTR(NWSURF)      ! Controls whether intersection data need to
                               ! be updated here
      INTEGER, INTENT (IN) ::
     >   LENXB                 ! Packed patch space limit

      REAL, INTENT (INOUT), DIMENSION (LENXB) ::
     >   XBASE, YBASE, ZBASE   ! Patch coordinates

      INTEGER, INTENT (INOUT) ::
     >   LASTXB, MOSTXB,       ! Packed patch indices (last, biggest, first)
     >   IXB(MXPATCH)

      INTEGER, INTENT (IN) ::  ! Subpatch info.
     >   MXSUBI, MXSUBJ

      INTEGER, INTENT (OUT) ::
     >   NSUBI(MXPATCH),
     >   NSUBJ(MXPATCH),
     >   ISUB(0:MXSUBI,MXPATCH),
     >   JSUB(0:MXSUBJ,MXPATCH)

      INTEGER, INTENT (INOUT) ::
     >   IMXCOM(MXPATCH), JMXCOM(MXPATCH) ! Patch dimensions

      INTEGER, INTENT (IN) ::
     >   IWRIT

      LOGICAL, INTENT (OUT) ::
     >   FAIL

C     Local constants:

      REAL,      PARAMETER :: HALF = 0.5, ONE = 1., ZERO = 0.
      CHARACTER, PARAMETER :: METHOD = 'B' ! For CHANGEN/PLSCRV3D

C     Local variables:

      INTEGER
     >   I, I1, I2, I2S, IB, IE, IER, ILE1, ILE2, IMAT, IU, IX, IXSHIFT,
     >   J, J1, J2, JB, JE, JM, JMAT,
     >   K, K1, K2, K3, K4, K5, K6, KASE, L, LUNERR, M,
     >   N, N1, N2, NFIX, NI3, NIN, NIT, NITOT, NJ3, NJA, NJB, NJM

      REAL
     >   DRANGE, P, Q, R, SMATRIX(4,4,3), UCUT, VCUT, XYZCUT(3),
     >   XLINE(NFSTN),  YLINE(NFSTN),  ZLINE(NFSTN),
     >   XLIN1(JLBODY), YLIN1(JLBODY), ZLIN1(JLBODY),
     >   XLIN2(JLBODY), YLIN2(JLBODY), ZLIN2(JLBODY)

      REAL, DIMENSION (JLBODY,NFSTN) ::
     >   UBODY, VBODY, XBODY

      REAL, ALLOCATABLE, DIMENSION (:,:) ::
     >   UEVAL, VEVAL

      LOGICAL
     >   BICUBIC

C     Storage:

      REAL, DIMENSION(4), SAVE ::
     >   UVRANGE

      DATA
     >   UVRANGE /ZERO, ONE, ZERO, ONE/ ! See PLBICUT

C     Execution:
C     ----------

      K1 = KPATCH ! First patch to fill
      K2 = K1 + 1
      K3 = K1 + 2
      K4 = K1 + 3
      K5 = K1 + 4
      K6 = K1 + 5

C     Intersection leading edges must be at the middle index for proper
C     grid perturbation.  Use of LCRANK > 0 ensures this.  But if the
C     intersection (x,y,z)s have been redistributed, we need to update the
C     (u,v)s and (i,j)s here.

C     Parameterize the body geometry:

      UBODY(1,1) = ZERO ! -999. suppresses normalization

      DO I = 1, NFSTN
         XBODY(1:JLBODY,I) = XF(I)
      END DO

      CALL PARAM2D (JLBODY, NFSTN, 1, JLBODY, 1, NFSTN, ! Note U/V order
     >              XBODY, YF, ZF, VBODY, UBODY)        ! U is streamwise

      BICUBIC = IDEGFUS == 3
      LUNERR  = -ABS (IWRIT)
      DRANGE  = XF(NFSTN) - XF(1)
      KASE    = 1 ! Find Zcut given X- & YCUT
      J       = JLBODY

      DO L = 1, 2

         IF (LOSTINTR(L)) THEN

            UCUT = UINTR(2,L)
            VCUT = VINTR(2,L)
            IE   = IINTR(2,L)
            JE   = JINTR(2,L)
            IMAT = 0
            JMAT = 0

            DO I = 2, NINTR(L) - 1

               XYZCUT(1) = XINTR(I,L)
               XYZCUT(2) = YINTR(I,L)
               XYZCUT(3) = ZINTR(I,L)

               IF (BICUBIC) THEN

                  CALL PLBICUT (KASE, J, NFSTN, 1, J, 1, NFSTN,
     >                          XBODY, YF, ZF, VBODY, UBODY, JE, IE,
     >                          JMAT, IMAT, SMATRIX,
     >                          EPS, DRANGE, XYZCUT, UVRANGE,
     >                          VCUT, UCUT, Q, P, LUNERR, IER)
               ELSE

                  CALL BILINCUT (KASE, J, NFSTN, 1, J, 1, NFSTN,
     >                           XBODY, YF, ZF, VBODY, UBODY, JE, IE,
     >                           EPS, DRANGE, XYZCUT, UVRANGE,
     >                           VCUT, UCUT, Q, P, LUNERR, IER)
               END IF

               IF (IER > 2) THEN
                  IF (IWRIT > 0) THEN
                     WRITE (IWRIT,
     >                  '(2(/, A, I2), /, A, I4, /, A, 1P, 3E15.6)')
     >                  ' FCTV2: PLBICUT/BILINCUT error =', IER,
     >                  ' Wing/fin surface =', L,
     >                  ' Intersection pt. =', I,
     >                  ' Corresp.(x,y,z)  =', XYZCUT
                  END IF

                  CALL SYNCH
                  STOP

               END IF

               UINTR(I,L) = UCUT
               VINTR(I,L) = VCUT
               IINTR(I,L) = IE
               JINTR(I,L) = JE

            END DO

         END IF

         KASE = 3 ! Find YCUT given ZCUT and XCUT (for the fin, if redone)

      END DO

      N1 = NINTR(1)
      N2 = NINTR(2)

      ILE1 = (N1 + 1) / 2
      ILE2 = (N2 + 1) / 2

C     Establish the dimensions of all subpatches:

      IU  = N1 - ILE2 + 1               ! Index of wing intersection <-> tail LE
      NI3 = IU - ILE1 + 1               ! Length of patch K3 forward of tail
      NIN = IINTR(ILE1,1) + 1           ! Length of nose panel

      I1  = MIN (IINTR(1,1), IINTR(N1,1))
      I2  = MIN (IINTR(1,2), IINTR(N2,2))
      NIT = NFSTN - MIN (I1, I2) + 1    ! Length of tail panel

      IMXCOM(K1) = NIN  ! Nose
      IMXCOM(K2) = ILE1 ! Below wing
      IMXCOM(K3) = NI3  ! Above forward part of wing
      IMXCOM(K4) = ILE2 ! Below tail
      IMXCOM(K5) = ILE2 ! Above tail
      IMXCOM(K6) = NIT  ! Aft of tail

      NSUBI(K2)  = 2

      ISUB(1,K1) = NIN  ! Subpatch indices
      ISUB(1,K2) = NI3
      ISUB(1,K3) = NI3
      ISUB(2,K2) = ILE1
      ISUB(1,K4) = ILE2
      ISUB(1,K5) = ILE2
      ISUB(1,K6) = NIT
      I2S        = NIN + NI3  - 1  ! I split index in UEVAL(*,J) for patch K2
      I2         = NIN + ILE1 - 1  ! I index in UEVAL(*,J) at end of patch K2
      NITOT      = I2  + NIT  - 1

      NJB = JINTR(ILE1,1) + 1      ! # js in patch below wing
      NJA = JLBODY - JINTR(ILE2,2) ! # js in patch above tail
      J2  = JLBODY - NJA + 1       ! Nose patch j index at tail LE
      NJM = J2     - NJB + 1       ! # js in middle (sub)patches
      NJ3 = NJM + NJA - 1          ! # js in patch K3

      JMXCOM(K1) = JLBODY          ! Nose
      JMXCOM(K2) = NJB             ! Below wing
      JMXCOM(K3) = NJ3             ! Above forward part of wing
      JMXCOM(K4) = NJM             ! Between wing and tail
      JMXCOM(K5) = NJA             ! Above tail
      JMXCOM(K6) = JLBODY          ! Aft body

      NSUBJ(K1)  = 3
      NSUBJ(K3)  = 2
      NSUBJ(K6)  = 3

      JSUB(1,K1) = NJB
      JSUB(2,K1) = J2
      JSUB(3,K1) = JLBODY
      JSUB(1,K2) = NJB
      JSUB(1,K3) = NJM
      JSUB(2,K3) = NJ3
      JSUB(1,K4) = NJM
      JSUB(1,K5) = NJA
      JSUB(1,K6) = NJB
      JSUB(2,K6) = J2
      JSUB(3,K6) = JLBODY

C     Check that there is room for all fuselage panels:

      IX = LASTXB + 1
      I1 = IX
      M  = 0

      DO L = K1, K6
         N = IMXCOM(L) * JMXCOM(L)
         M = N + M
         IXB(L) = I1
         I1 = I1 + N
      END DO

      CALL XBCHECK (LENXB, MXPATCH, K6, LASTXB, M, 1,
     >              IFUS, IWRIT, MOSTXB)

      LASTXB = LASTXB + M ! Upon return from FCTV2

C     Strategy for the subpatches:
C
C     >  Work from keel to crown, forward to aft.
C     >  Establish the corner points in (u,v) space.
C     >  Generate appropriate edges  in (u,v) space.
C     >  Generate interior (u,v)s via TFI.
C     >  Interpolate interior and edge (x,y,z)s (bilinear or bicubic).

C     Allocate (u,v)s for interpolating all patches (edges overlapping):

      ALLOCATE (UEVAL(NITOT,JLBODY), VEVAL(NITOT,JLBODY))

      DO J = 1, JLBODY
         UEVAL(1,J)      = ZERO
         UEVAL(NITOT,J)  = ONE
      END DO

      DO I = 1, NITOT
         VEVAL(I,1)      = ZERO
         VEVAL(I,JLBODY) = ONE
      END DO

C     Deal with the nose point or line for three subpatches:
C     ------------------------------------------------------

      I = IX

      IF (ZF(2,1) == ZF(1,1)) THEN ! Point

         DO J = 1, JLBODY
            XBASE(I) = XF(1)
            YBASE(I) = YF(1,1)
            ZBASE(I) = ZF(1,1)
            I = I + NIN
            VEVAL(1,J) = VBODY(J,1)
         END DO

      ELSE ! Singular line needs redistribution

         JM = (JLBODY + 1) / 2 ! JLBODY is odd to capture the water-line

         DO J = 1, JLBODY
            XLIN1(J) = XF(1)
            YLIN1(J) = YF(J,1)
            ZLIN1(J) = ZF(J,1)
         END DO

         XLIN2(1)   = XLIN1(1)
         YLIN2(1)   = YLIN1(1)
         ZLIN2(1)   = ZLIN1(1)
         XLIN2(NJB) = XLIN1(JM)
         YLIN2(NJB) = YLIN1(JM)
         ZLIN2(NJB) = ZLIN1(JM)

         XLIN2(J2)     = XLIN1(JM)
         YLIN2(J2)     = YLIN1(JM)
         ZLIN2(J2)     = ZLIN1(JM)
         XLIN2(JLBODY) = XLIN1(JLBODY)
         YLIN2(JLBODY) = YLIN1(JLBODY)
         ZLIN2(JLBODY) = ZLIN1(JLBODY)

         CALL CHANGEN (1,  JM,     XLIN1, YLIN1, ZLIN1,
     >                 1,  NJB,    XLIN2, YLIN2, ZLIN2, METHOD)

         CALL CHANGEN (JM, JLBODY, XLIN1, YLIN1, ZLIN1,
     >                 J2, JLBODY, XLIN2, YLIN2, ZLIN2, METHOD)

         R = VBODY(JM,1) / REAL (NJB - 1)
         DO J = 1, NJB
            XBASE(I) = XLIN2(J)
            YBASE(I) = YLIN2(J)
            ZBASE(I) = ZLIN2(J)
            I = I + NIN
            VEVAL(1,J) = R * REAL (J - 1)
         END DO

         DO J = NJB + 1, J2
            XBASE(I) = XLIN2(NJB)
            YBASE(I) = YLIN2(NJB)
            ZBASE(I) = ZLIN2(NJB)
            I = I + NIN
            VEVAL(1,J) = VBODY(JM,1)
         END DO

         R = (ONE - VBODY(JM,1)) / REAL (NJA - 1)
         DO J = J2 + 1, JLBODY
            XBASE(I) = XLIN2(J)
            YBASE(I) = YLIN2(J)
            ZBASE(I) = ZLIN2(J)
            I = I + NIN
            VEVAL(1,J) = VBODY(JM,1) + R * REAL (J - J2)
         END DO

      END IF

C     First subpatch (on the nose at the keel):
C     -----------------------------------------

C     Upper right corner (wing LE):

      UEVAL(NIN,NJB) = UINTR(ILE1,1)
      VEVAL(NIN,NJB) = VINTR(ILE1,1)

C     Right-hand edge (u,v)s:

      I =  NIN - 1
      J =  JINTR(ILE1,1)
      R =  VINTR(ILE1,1) / REAL (NJB - 1)
      P = (UINTR(ILE1,1) - UBODY(J,I)) / (UBODY(J,NIN)  - UBODY(J,I))

      DO J = 1, NJB - 1
         UEVAL(NIN,J) = (ONE - P) * UBODY(J,I) + P * UBODY(J,NIN)
         VEVAL(NIN,J) = R * REAL (J - 1)
      END DO

C     Do the rest of the right-hand edge while we're at it:

      Q = VEVAL(NIN,NJB)
      R = (ONE - Q) / REAL (NJ3 - 1)

      DO J = NJB + 1, JLBODY
         UEVAL(NIN,J) = (ONE - P) * UBODY(J,I) + P * UBODY(J,NIN)
         VEVAL(NIN,J) = Q + R * REAL (J - NJB)
      END DO

C     Streamwise edge (u,v)s of nose subpatches:

      R = ONE / (XINTR(ILE1,1) - XF(1)) ! Straight (u,v) lines suffice

      DO I = 2, NIN - 1
         P = (XF(I) - XF(1)) * R
         UEVAL(I,1)   = UBODY(1,I)                      ! Keel
         UEVAL(I,NJB) = P * UEVAL(NIN,NJB)              ! First J split
         VEVAL(I,NJB) = P * Q + (ONE - P) * VEVAL(1,NJB)
         UEVAL(I,J2)  = P * UEVAL(NIN,J2)               ! Second J split
         VEVAL(I,J2)  = P * VEVAL(NIN,J2) + (ONE - P) * VEVAL(1,J2)
         UEVAL(I,JLBODY) = UBODY(JLBODY,I)              ! Crown
      END DO

C     Panel K2 under the wing:
C     ------------------------

C     Close the trailing edge in (u,v) space:

      UEVAL(I2,NJB) = (UINTR(1,1) + UINTR(N1,1)) * HALF
      VEVAL(I2,NJB) = (VINTR(1,1) + VINTR(N1,1)) * HALF

C     Right edge (u,v)s:

      I =  IINTR(1,1)
      J =  JINTR(1,1)
      P = (UEVAL(I2,NJB) - UBODY(J,I)) / (UBODY(J,I+1)  - UBODY(J,I))
      R =  VEVAL(I2,NJB) / REAL (NJB - 1)

      DO J = 1, NJB - 1
         UEVAL(I2,J) = (ONE - P) * UBODY(J,I) + P * UBODY(J,I+1)
         VEVAL(I2,J) = R * REAL (J - 1)
      END DO

C     Keel (u,v)s + upper patch edge (= lower half of intersection).
C     Use "arc length" along the intersection (u,v)s in case U(LE) is not
C     at the same index as X(LE).

      IX = IXB(K2) ! Start of ample work-space for arc lengths

      CALL CHORDS2D (ILE1, UINTR(1,1), VINTR(1,1), .TRUE., P, XBASE(IX))

      P  = UEVAL(NIN,1)
      R  = UEVAL(I2,1) - P
      I1 = IX + ILE1 - 1
      J  = ILE1

      DO I = NIN + 1, I2 - 1
         I1 = I1 - 1
         UEVAL(I,1)   = P + R * (ONE - XBASE(I1))
         J  = J - 1
         UEVAL(I,NJB) = UINTR(J,1)
         VEVAL(I,NJB) = VINTR(J,1)
      END DO


C     Fill the interior of patch K2 now so we can reuse its upper edge,
C     U/VEVAL(*,NJB), for the lower edge of patch K3.  Do the nose too:
C     -----------------------------------------------------------------

      DO J = 1, 3

         CALL TFI2D (NITOT, 1, NIN, JSUB(J-1,K1), JSUB(J,K1),
     >               UEVAL, VEVAL, XBASE(IX))
      END DO

      CALL TFI2D (NITOT, NIN, I2, 1, NJB, UEVAL, VEVAL, XBASE(IX)) ! K2 interior


C     Patch K3 above forward part of wing:
C     ------------------------------------

C     Lower edge of K3 & K4 is the upper part of the wing/body intersection:

      DO I = 1, ILE1 - 2 ! Don't clobber the closed TE
         UEVAL(NIN + I,NJB) = UINTR(ILE1 + I,1)
         VEVAL(NIN + I,NJB) = VINTR(ILE1 + I,1)
      END DO

C     Right-hand edge (u,v)s of lower subpatch of K3 (straight/uniform):

      UEVAL(I2S,J2) = UINTR(ILE2,2) ! Tail LE
      VEVAL(I2S,J2) = VINTR(ILE2,2)
      R = ONE / REAL (NJM - 1)
      P = (UEVAL(I2S,J2) - UEVAL(I2S,NJB)) * R
      Q = (VEVAL(I2S,J2) - VEVAL(I2S,NJB)) * R

      DO J = NJB + 1, J2 - 1
         R = REAL (J - NJB)
         UEVAL(I2S,J) = UEVAL(I2S,NJB) + R * P
         VEVAL(I2S,J) = VEVAL(I2S,NJB) + R * Q
      END DO

C     Upper edge of lower subpatch of K3 comes from TFI.

C     Right-hand edge of upper subpatch of K3:

      I =  IINTR(ILE2,2)
      J =  JINTR(ILE2,2)
      P = (UINTR(ILE2,2) - UBODY(J,I)) / (UBODY(J,I+1) - UBODY(J,I))
      Q =  VINTR(ILE2,2)
      R = (ONE - Q) / REAL (NJA - 1)

      DO J = J2 + 1, JLBODY
         UEVAL(I2S,J) = (ONE - P) * UBODY(J,I) + P * UBODY(J,I+1)
         VEVAL(I2S,J) = Q + R * REAL (J - J2)
      END DO

C     Crown line (u,v)s match the forward wing intersection arc distribution:

      CALL CHORDS2D (NI3, UINTR(ILE1,1), VINTR(ILE1,1), .TRUE., P,
     >               XBASE(IX))
      P  = UINTR(ILE1,1)
      Q  = UEVAL(NIN,JLBODY)
      R  = UEVAL(I2S,JLBODY) - Q
      I1 = IX

      DO I = NIN + 1, I2S - 1
         I1 = I1 + 1
         UEVAL(I,JLBODY) = Q + R * XBASE(I1)
      END DO

C     Interiors of both K3 subpatches:

      CALL TFI2D (NITOT, NIN, I2S, NJB, JLBODY, UEVAL, VEVAL, XBASE(IX))


C     Panel K4 below the tail/above the aft wing:
C     -------------------------------------------

C     Close the trailing edge:

      UEVAL(I2,J2) = (UINTR(1,2) + UINTR(N2,2)) * HALF
      VEVAL(I2,J2) = (VINTR(1,2) + VINTR(N2,2)) * HALF

C     Right edge (linear/uniform):

      R = ONE / REAL (NJM - 1)
      P = (UEVAL(I2,J2) - UEVAL(I2,NJB)) * R
      Q = (VEVAL(I2,J2) - VEVAL(I2,NJB)) * R

      DO J = NJB + 1, J2 - 1
         R = REAL (J - NJB)
         UEVAL(I2,J) = UEVAL(I2,NJB) + R * P
         VEVAL(I2,J) = VEVAL(I2,NJB) + R * Q
      END DO

C     Upper edge is the lower edge of the tail/body intersection:

      DO I = 1, ILE2 - 2
         UEVAL(I2-I,J2) = UINTR(I+1,2)
         VEVAL(I2-I,J2) = VINTR(I+1,2)
      END DO

C     Interior of panel K4:

      CALL TFI2D (NITOT, I2S, I2, NJB, J2, UEVAL, VEVAL, XBASE(IX))


C     Panel above the tail:
C     ---------------------

C     Right edge (u,v)s:

      I =  IINTR(N2,2)
      J =  JINTR(N2,2)
      Q =  VEVAL(I2,J2)
      R = (ONE - Q) / REAL (NJA - 1)
      P = (UEVAL(I2,J2) - UBODY(J,I)) / (UBODY(J,I+1) - UBODY(J,I))

      DO J = J2 + 1, JLBODY
         UEVAL(I2,J) = (ONE - P) * UBODY(J,I) + P * UBODY(J,I+1)
         VEVAL(I2,J) = Q + R * REAL (J - J2)
      END DO

C     Crown line + lower patch edge (= upper half of intersection):

      J  =  JLBODY
      P  =  UEVAL(I2S,J)
      R  = (UEVAL(I2,J) - P) / (UEVAL(I2,J2) - UEVAL(I2S,J2))
      I1 =  ILE2

      DO I = I2S + 1, I2 - 1
         I1 = I1 + 1
         UEVAL(I,J2) = UINTR(I1,2)
         VEVAL(I,J2) = VINTR(I1,2)
         UEVAL(I,J)  = P + R * (UINTR(I1,2) - UINTR(ILE2,2))
      END DO

C     Interior of panel K5:

      CALL TFI2D (NITOT, I2S, I2, J2, J, UEVAL, VEVAL, XBASE(IX))


C     Aft body panel in three subpatches:
C     -----------------------------------

C     Upper right corner:

      I1 = IINTR(1,1)
      J1 = JINTR(1,1)
      Q  = (VEVAL(I2,NJB)  - VBODY(J1,I1)) /
     >     (VBODY(J1+1,I1) - VBODY(J1,I1))
      P  = (ONE - Q) * VBODY(J1,  NFSTN) +   ! P = VEVAL(NITOT,NJB)
     >            Q  * VBODY(J1+1,NFSTN)

C     Right-hand edge:

      R = P / REAL (NJB - 1)

      DO J = 2, NJB
         VEVAL(NITOT,J) = R * REAL (J - 1)
      END DO

C     Keel line and edge aft of wing TE:

      IF (NIT == NFSTN - I1 + 1) THEN

         DO I = I2 + 1, NITOT - 1
            I1 = I1 + 1
            UEVAL(I,1)   = UBODY(1,I1)
            UEVAL(I,NJB) = (ONE - Q) * UBODY(J1,I1) + Q * UBODY(J1+1,I1)
            VEVAL(I,NJB) = (ONE - Q) * VBODY(J1,I1) + Q * VBODY(J1+1,I1)
         END DO

      ELSE ! Adjust the number of points

         XLINE(I1) = UEVAL(I2,1) ! Keel line first
         YLINE(I1) = ZERO

         DO I = I1 + 1, NFSTN
            XLINE(I) = UBODY(1,I)
            YLINE(I) = ZERO
         END DO

         ZLINE = ZERO ! For all I, to be safe

         CALL CHANGEN (I1, NFSTN, XLINE, YLINE, ZLINE,
     >                 I2, NITOT, UEVAL, VEVAL, ZLINE, METHOD)

         XLINE(I1) = UEVAL(I2,NJB) ! Now for wing TE --> end of body
         YLINE(I1) = VEVAL(I2,NJB)

         DO I = I1 + 1, NFSTN
            XLINE(I) = (ONE - Q) * UBODY(J1,I) + Q * UBODY(J1+1,I)
            YLINE(I) = (ONE - Q) * VBODY(J1,I) + Q * VBODY(J1+1,I)
         END DO

         CALL CHANGEN (I1, NFSTN, XLINE, YLINE, ZLINE,
     >                 I2, NITOT, UEVAL(1,NJB), VEVAL(1,NJB), ZLINE,
     >                 METHOD)
      END IF

C     Likewise for the upper aft-body subpatch edges:

C     Lower right corner:

      I1 = IINTR(1,2)
      J1 = JINTR(1,2)
      Q  = (VEVAL(I2,J2)   - VBODY(J1,I1)) /
     >     (VBODY(J1+1,I1) - VBODY(J1,I1))
      P  = (ONE - Q) * VBODY(J1,  NFSTN) +   ! P = VEVAL(NITOT,J2)
     >            Q  * VBODY(J1+1,NFSTN)

C     Right-hand edge:

      R = (ONE - P) / REAL (NJA - 1)

      DO J = J2, JLBODY - 1
         VEVAL(NITOT,J) = P + R * REAL (J - J2)
      END DO

C     Crown line and edge aft of tail TE:

      IF (NIT == NFSTN - I1 + 1) THEN

         DO I = I2 + 1, NITOT - 1
            I1 = I1 + 1
            UEVAL(I,J2) = (ONE - Q) * UBODY(J1,I1) + Q * UBODY(J1+1,I1)
            VEVAL(I,J2) = (ONE - Q) * VBODY(J1,I1) + Q * VBODY(J1+1,I1)
            UEVAL(I,JLBODY) = UBODY(JLBODY,I1)
         END DO

      ELSE ! Adjust the number of points

         XLINE(I1) = UEVAL(I2,J2) ! Tail TE --> end of body
         YLINE(I1) = VEVAL(I2,J2)

         DO I = I1 + 1, NFSTN
            XLINE(I) = (ONE - Q) * UBODY(J1,I) + Q * UBODY(J1+1,I)
            YLINE(I) = (ONE - Q) * VBODY(J1,I) + Q * VBODY(J1+1,I)
         END DO

         CALL CHANGEN (I1, NFSTN, XLINE, YLINE, ZLINE,
     >                 I2, NITOT, UEVAL(1,J2), VEVAL(1,J2), ZLINE,
     >                 METHOD)

         XLINE(I1) = UEVAL(I2,JLBODY) ! Now for the crown line
         YLINE(I1) = ONE

         DO I = I1 + 1, NFSTN
            XLINE(I) = UBODY(JLBODY,I)
            YLINE(I) = ONE
         END DO

         CALL CHANGEN (I1, NFSTN, XLINE, YLINE, ZLINE,
     >                 I2, NITOT, UEVAL(1,JLBODY), VEVAL(1,JLBODY),
     >                 ZLINE, METHOD)
      END IF

C     Aft edge of middle aft-body subpatch:

      P =  VEVAL(NITOT,NJB)
      R = (VEVAL(NITOT,J2) - P) / REAL (NJM - 1)

      DO J = NJB + 1, J2 - 1
         VEVAL(NITOT,J) = P + R * REAL (J - NJB)
      END DO

C     Interior (u,v)s of the aft-body subpatches:

      DO J = 1, 3

         CALL TFI2D (NITOT, I2, NITOT, JSUB(J-1,K6), JSUB(J,K6),
     >               UEVAL, VEVAL, XBASE(IX))
      END DO

!!!      write (60) 3
!!!      write (60) nitot, jlbody, n1, 1, n2, 1
!!!      write (60) ueval, veval
!!!      write (60) uintr(1:n1,1), vintr(1:n1,1)
!!!      write (60) uintr(1:n2,2), vintr(1:n2,2)

C     Parametric bilinear/bicubic interpolation of (x,y,z)s for all panels:
C     ---------------------------------------------------------------------

      K = K1
      IMAT = 0
      JMAT = 0
      IXSHIFT = 0

      DO KASE = 1, 6

         SELECT CASE (KASE)

            CASE (1)
               I1 = 2
               L  = NIN
               J1 = 1
               M  = JLBODY
               IB = 2
               JB = 1
            CASE (2)
               I1 = NIN
               L  = I2
!              J1 = 1 (still)
               M  = NJB - 1
               IB = IINTR(ILE1,1)
               JB = 1
            CASE (3)
               L  = I2S
               J1 = NJB + 1
               M  = JLBODY
               IB = IINTR(ILE1,1)
               JB = JINTR(ILE1,1)
               IXSHIFT = NI3 ! Since J is not starting at the edge
            CASE (4)
               I1 = I2S
               L  = I2
               M  = J2 - 1
               IB = IINTR(IU,1)
               JB = JINTR(IU,1)
               IXSHIFT = ILE2
            CASE (5)
               J1 = J2  + 1
               M  = JLBODY
               IB = IINTR(ILE2,2)
               JB = JINTR(ILE2,2)
            CASE (6)
               I1 = I2
               L  = NITOT
               J1 = 1
               IB = IINTR(1,1)
               JB = 1
               IXSHIFT = 0

         END SELECT

         IX = IXB(K) + IXSHIFT

         DO J = J1, M

            IF (K == K1) IX = IX + 1 ! Awkwardness from omitting nose pts.

            IE = IB
            JE = JB

            DO I = I1, L

               IF (BICUBIC) THEN

                  CALL PLBICUBE (0, JLBODY, NFSTN, 1, JLBODY, 1, NFSTN,
     >                           XBODY, YF, ZF, VBODY, UBODY,
     >                           VEVAL(I,J), UEVAL(I,J), JE, IE, 
     >                           JMAT, IMAT, SMATRIX, EPS, Q, P,
     >                           XBASE(IX), YBASE(IX), ZBASE(IX),
     >                           XLINE, YLINE, IER)
               ELSE

                  CALL PBILINT  (0, JLBODY, NFSTN, 1, JLBODY, 1, NFSTN,
     >                           XBODY, YF, ZF, VBODY, UBODY,
     >                           VEVAL(I,J), UEVAL(I,J), JE, IE,
     >                           EPS, Q, P,
     >                           XBASE(IX), YBASE(IX), ZBASE(IX),
     >                           XLINE, YLINE, IER)
               END IF

               IF (IER /= 0) THEN
                  IF (IWRIT > 0) THEN
                    WRITE (IWRIT,'(/,A,A,I2,A,4I5,/,A,2F13.5,A,3F13.5)')
     >               ' FCTV2: Marginal 2D interpolation. ',
     >               ' IER:', IER, '; I, J, IE, JE:', I, J, IE, JE,
     >               ' U, V:', UEVAL(I,J), VEVAL(I,J),
     >               '   X, Y, Z used:', XBASE(IX), YBASE(IX), ZBASE(IX)
                  END IF
               END IF

               IF (I == I1) THEN
                  IB = IE ! Ready for the next J
                  JB = JE
               END IF

               IX = IX + 1

            END DO

         END DO

         K = K + 1

      END DO ! Next panel


      DEALLOCATE (UEVAL, VEVAL)


C     Transcribe intersections (x,y,z)s and ensure common edge points:

      IX = IXB(K1) + NIN * NJB - 1   ! Wing LE pt. of nose panel
      XBASE(IX) = XINTR(ILE1,1)
      YBASE(IX) = YINTR(ILE1,1)
      ZBASE(IX) = ZINTR(ILE1,1)

      IX = IXB(K3) + NI3 * NJM  - 1  ! Tail LE pt. of K3 panel
      XBASE(IX) = XINTR(ILE2,2)
      YBASE(IX) = YINTR(ILE2,2)
      ZBASE(IX) = ZINTR(ILE2,2)

      IX = IXB(K6) + NIT * (NJB - 1) ! Wing TE pt. from aft-body panel
      I1 = IXB(K3) - 1               ! Corresp. pt. in under-wing panel
      XBASE(I1) = XBASE(IX)
      YBASE(I1) = YBASE(IX)
      ZBASE(I1) = ZBASE(IX)

      I1 = IXB(K4) + ILE2 - 1        ! Corresp. pt. of panel between wing & tail
      XBASE(I1) = XBASE(IX)
      YBASE(I1) = YBASE(IX)
      ZBASE(I1) = ZBASE(IX)

      IX = IXB(K6) + NIT * (J2 - 1)  ! Tail TE pt. from aft-body panel
      I1 = IXB(K5) - 1               ! Corresp. pt. in panel between wing & tail
      XBASE(I1) = XBASE(IX)
      YBASE(I1) = YBASE(IX)
      ZBASE(I1) = ZBASE(IX)

      I1 = IXB(K5) + ILE2 - 1        ! Corresp. pt. of panel above tail
      XBASE(I1) = XBASE(IX)
      YBASE(I1) = YBASE(IX)
      ZBASE(I1) = ZBASE(IX)

      IX = IXB(K3) - 1               ! Lower wing intersection (except TE)
      DO I = 2, ILE1
         IX = IX - 1
         XBASE(IX) = XINTR(I,1)
         YBASE(IX) = YINTR(I,1)
         ZBASE(IX) = ZINTR(I,1)
      END DO

      IX = IXB(K3)                   ! Upper wing intersection, forward
      DO I = ILE1, IU
         XBASE(IX) = XINTR(I,1)
         YBASE(IX) = YINTR(I,1)
         ZBASE(IX) = ZINTR(I,1)
         IX = IX + 1
      END DO

      IX = IXB(K4)                   ! Upper wing int., aft (except TE)
      DO I = IU, N1 - 1
         XBASE(IX) = XINTR(I,1)
         YBASE(IX) = YINTR(I,1)
         ZBASE(IX) = ZINTR(I,1)
         IX = IX + 1
      END DO

      IX = IXB(K5) - 1               ! Lower tail intersection (except TE)
      DO I = 2, ILE2
         IX = IX - 1
         XBASE(IX) = XINTR(I,2)
         YBASE(IX) = YINTR(I,2)
         ZBASE(IX) = ZINTR(I,2)
      END DO

      IX = IXB(K5)                   ! Upper tail intersection (except TE)
      DO I = ILE2, N2 - 1
         XBASE(IX) = XINTR(I,2)
         YBASE(IX) = YINTR(I,2)
         ZBASE(IX) = ZINTR(I,2)
         IX = IX + 1
      END DO

      
C     Fix up extreme skewing at blunt wing leading edge?

      NFIX = NFIXL(1) + NFIXR(1) + 1

      IF (NFIX > 2) THEN ! Patches K1, K1+1, K1+2 are partially revised

         CALL FIXROOT (MXPATCH, LENXB, IXB, IMXCOM, JMXCOM, K1,
     >                 JLBODY, NFIXL(1), NFIXR(1), NFIX, NFIXA(1), 
     >                 NFIXB(1), EPS, XBASE, YBASE, ZBASE)
      END IF


C     Finally:  split the wing patches implicitly at the point corresp. to I2S:

      NSUBI(1:2) = 2
      ISUB(1,1)  = ILE2
      ISUB(2,1)  = ILE1
      ISUB(1,2)  = NI3
      ISUB(2,2)  = ILE1

      KPATCH = K6 + 1 ! Next available patch
      FAIL   = .FALSE.

      END SUBROUTINE FCTV2

C***********************************************************************
C
      SUBROUTINE FCTV3 (NWSURF, MXINTR, MXPATCH, KPATCH, NFSTN,
     >                  JLBODY, IFUS, IDEGFUS, EPS, XF, YF, ZF,
     >                  NFIXL, NFIXR, NFIXA, NFIXB, NINTR,
     >                  XINTR, YINTR, ZINTR, UINTR, VINTR, IINTR, JINTR,
     >                  LOSTINTR, LENXB, XBASE, YBASE, ZBASE,
     >                  LASTXB, MOSTXB, IXB, MXSUBI, MXSUBJ,
     >                  NSUBI, NSUBJ, ISUB, JSUB, IMXCOM, JMXCOM,
     >                  IWRIT, FAIL)
C
C     FCTV3 imposes regular surface patches on the fuselage of a CTV (Crew
C     Transfer Vehicle)-type configuration with overlapping wing and tail
C     surfaces.  FCTV3 is a variant of FCTV2 for the case where the vertical
C     tail is at the center-line.  Standard paneling (NNCONFIG = 0) of this
C     case suffers from "grid shock" opposite the wing leading edge after the
C     FSURF results are redistributed to accommodate the dorsal fin.
C
C     Simplifying assumptions:
C
C     >  Just two wing-type surfaces, with fin surface 2 at the crown.
C     >  ILWING(1) >> ILWING(2).  E.g., double.  (ILWING = NINTR here.)
C     >  Both intersections employ the whole body (ISTN1/2 = 1:NFSTN)
C        because the whole body parameterization is used here.
C     >  If the nose is a singular line, NJBODY should be odd, but so
C        should JLBODY, so make them the same in all cases (unlike the
C        standard case where JLBODY = NJBODY + 1).
C
C     The driving program should check that these assumptions are met.
C
C     08/08/00  David Saunders  FCTV3 adapted from FCTV2.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   NWSURF,               ! # wing/nacelle surfaces
     >   MXINTR,               ! Max. # pts. per intersection
     >   MXPATCH               ! Max. # patches expected

      INTEGER, INTENT (INOUT) ::
     >   KPATCH                ! In & out with the next patch # to fill

      INTEGER, INTENT (IN) ::
     >   NFSTN,                ! # body defining sections
     >   JLBODY,               ! # output points per section of (unsplit) nose
     >                         ! panel; JLBODY = NJBODY = # J pts. input per
     >                         ! regularized geometry section (not NJBODY + 1
     >                         ! as in FSURF); JLBODY should be ODD
     >   IFUS,                 ! Component # for diagnostics
     >   IDEGFUS               ! 1 or 3 for [bi]linear or [bi]cubic fus. splines

      REAL, INTENT (IN) ::
     >   EPS,                  ! Tolerance for 2-D searches/interpolations
     >   XF(NFSTN),            ! Regularized body sections;
     >   YF(JLBODY,NFSTN),     ! Y/ZF(1:JLBODY,*) are input (not 1:JLBODY-1)
     >   ZF(JLBODY,NFSTN)

      INTEGER, INTENT (IN), DIMENSION (NWSURF) ::
     >   NFIXL, NFIXR,         ! FIXROOT controls, as for FSURF (wing only)
     >   NFIXA, NFIXB

      INTEGER, INTENT (IN) ::
     >   NINTR(NWSURF)         ! # pts. for each intersection; (1) >> (2)

      REAL, INTENT (IN), DIMENSION (MXINTR,NWSURF) ::
     >   XINTR, YINTR, ZINTR   ! Intersection coordinates

      REAL, INTENT (INOUT), DIMENSION (MXINTR,NWSURF) ::
     >   UINTR, VINTR          ! Corresp. (u,v)s referred to whole body;
                               ! modified here if LOSTINTR() = T

      INTEGER, INTENT (INOUT), DIMENSION (MXINTR,NWSURF) ::
     >   IINTR, JINTR          ! Body cell indices for intersection points;
                               ! left half of fin uses NHALF(2) : ILWING(2)
      LOGICAL, INTENT (IN) ::
     >   LOSTINTR(NWSURF)      ! Controls whether intersection data need to
                               ! be updated here
      INTEGER, INTENT (IN) ::
     >   LENXB                 ! Packed patch space limit

      REAL, INTENT (INOUT), DIMENSION (LENXB) ::
     >   XBASE, YBASE, ZBASE   ! Patch coordinates

      INTEGER, INTENT (INOUT) ::
     >   LASTXB, MOSTXB,       ! Packed patch indices (last, biggest, first)
     >   IXB(MXPATCH)

      INTEGER, INTENT (IN) ::  ! Subpatch info.
     >   MXSUBI, MXSUBJ

      INTEGER, INTENT (OUT) ::
     >   NSUBI(MXPATCH),
     >   NSUBJ(MXPATCH),
     >   ISUB(0:MXSUBI,MXPATCH),
     >   JSUB(0:MXSUBJ,MXPATCH)

      INTEGER, INTENT (INOUT) ::
     >   IMXCOM(MXPATCH), JMXCOM(MXPATCH) ! Patch dimensions

      INTEGER, INTENT (IN) ::
     >   IWRIT

      LOGICAL, INTENT (OUT) ::
     >   FAIL

C     Local constants:

      REAL,      PARAMETER :: HALF = 0.5, ONE = 1., ZERO = 0.
      CHARACTER, PARAMETER :: METHOD = 'B' ! For CHANGEN/PLSCRV3D

C     Local variables:

      INTEGER
     >   I, I1, I2, I2S, IB, IE, IER, ILE1, ILE2, IMAT, IU, IX, IXSHIFT,
     >   J, J1, J2, JB, JE, JMAT,
     >   K, K1, K2, K3, K4, K5, KASE, L, LUNERR, M,
     >   N, N1, N2, NFIX, NI3, NIN, NIT, NITOT, NJ3, NJB

      REAL
     >   DRANGE, P, Q, R, SMATRIX(4,4,3), UCUT, VCUT, XYZCUT(3),
     >   XLINE(NFSTN),  YLINE(NFSTN),  ZLINE(NFSTN),
     >   XLIN1(JLBODY), YLIN1(JLBODY), ZLIN1(JLBODY),
     >   XLIN2(JLBODY), YLIN2(JLBODY), ZLIN2(JLBODY)

      REAL, DIMENSION (JLBODY,NFSTN) ::
     >   UBODY, VBODY, XBODY

      REAL, ALLOCATABLE, DIMENSION (:,:) ::
     >   UEVAL, VEVAL

      LOGICAL
     >   BICUBIC

C     Storage:

      REAL, DIMENSION(4), SAVE ::
     >   UVRANGE

      DATA
     >   UVRANGE /ZERO, ONE, ZERO, ONE/ ! See PLBICUT

C     Execution:
C     ----------

      K1 = KPATCH ! First patch to fill
      K2 = K1 + 1
      K3 = K1 + 2
      K4 = K1 + 3
      K5 = K1 + 4

C     Intersection leading edges must be at the middle index for proper
C     grid perturbation.  Use of LCRANK > 0 ensures this.  But if the
C     intersection (x,y,z)s have been redistributed, we need to update the
C     (u,v)s and (i,j)s here.

C     Parameterize the body geometry:

      UBODY(1,1) = ZERO ! -999. suppresses normalization

      DO I = 1, NFSTN
         XBODY(1:JLBODY,I) = XF(I)
      END DO

      CALL PARAM2D (JLBODY, NFSTN, 1, JLBODY, 1, NFSTN, ! Note U/V order
     >              XBODY, YF, ZF, VBODY, UBODY)        ! U is streamwise

      BICUBIC = IDEGFUS == 3
      LUNERR  = -ABS (IWRIT)
      DRANGE  = XF(NFSTN) - XF(1)
      KASE    = 1 ! Find Zcut given X- & YCUT

      N1   = NINTR(1)
      N2   = NINTR(2)
      ILE1 = (N1 + 1) / 2
      ILE2 = (N2 + 1) / 2
      I1   = 2
      I2   = N1 - 1
      J    = JLBODY

      DO L = 1, 2

         IF (LOSTINTR(L)) THEN

            UCUT = UINTR(I1,L)
            VCUT = VINTR(I1,L)
            IE   = IINTR(I1,L)
            JE   = JINTR(I1,L)
            IMAT = 0
            JMAT = 0

            DO I = I1, I2

               XYZCUT(1) = XINTR(I,L)
               XYZCUT(2) = YINTR(I,L)
               XYZCUT(3) = ZINTR(I,L)

               IF (BICUBIC) THEN

                  CALL PLBICUT (KASE, J, NFSTN, 1, J, 1, NFSTN,
     >                          XBODY, YF, ZF, VBODY, UBODY, JE, IE,
     >                          JMAT, IMAT, SMATRIX,
     >                          EPS, DRANGE, XYZCUT, UVRANGE,
     >                          VCUT, UCUT, Q, P, LUNERR, IER)
               ELSE

                  CALL BILINCUT (KASE, J, NFSTN, 1, J, 1, NFSTN,
     >                           XBODY, YF, ZF, VBODY, UBODY, JE, IE,
     >                           EPS, DRANGE, XYZCUT, UVRANGE,
     >                           VCUT, UCUT, Q, P, LUNERR, IER)
               END IF

               IF (IER > 2) THEN
                  IF (IWRIT > 0) THEN
                     WRITE (IWRIT,
     >                  '(2(/, A, I2), /, A, I4, /, A, 1P, 3E15.6)')
     >                  ' FCTV3: PLBICUT/BILINCUT error =', IER,
     >                  ' Wing/fin surface =', L,
     >                  ' Intersection pt. =', I,
     >                  ' Corresp.(x,y,z)  =', XYZCUT
                  END IF

                  CALL SYNCH
                  STOP

               END IF

               UINTR(I,L) = UCUT
               VINTR(I,L) = VCUT
               IINTR(I,L) = IE
               JINTR(I,L) = JE

            END DO

         END IF

         KASE = 3 ! Find YCUT given ZCUT and XCUT (for the fin, if redone)
         I1   = ILE2 + 1
         I2   = N2 - 1

      END DO

C     Establish the dimensions of all subpatches:

      IU  = N1 - ILE2 + 1               ! Index of wing intersection <-> fin LE
      NI3 = IU - ILE1 + 1               ! Length of patch K3 forward of fin
      NIN = IINTR(ILE1,1) + 1           ! Length of nose panel

      I1  = MIN (IINTR(1,1), IINTR(N1,1))
      I2  = IINTR(N2,2)
      NIT = NFSTN - MIN (I1, I2) + 1    ! Length of aft-body panel

      IMXCOM(K1) = NIN  ! Nose
      IMXCOM(K2) = ILE1 ! Below wing
      IMXCOM(K3) = NI3  ! Above forward part of wing
      IMXCOM(K4) = ILE2 ! Below fin
      IMXCOM(K5) = NIT  ! Aft of fin

      NSUBI(K2)  = 2

      ISUB(1,K1) = NIN  ! Subpatch indices
      ISUB(1,K2) = NI3
      ISUB(1,K3) = NI3
      ISUB(2,K2) = ILE1
      ISUB(1,K4) = ILE2
      ISUB(1,K5) = NIT
      I2S        = NIN + NI3  - 1  ! I split index in UEVAL(*,J) for patch K2
      I2         = NIN + ILE1 - 1  ! I index in UEVAL(*,J) at end of patch K2
      NITOT      = I2  + NIT  - 1

      NJB = JINTR(ILE1,1) + 1      ! # js in patch below wing
      NJ3 = JLBODY - NJB  + 1      ! # js in patches K3 : K5

      JMXCOM(K1) = JLBODY          ! Nose
      JMXCOM(K2) = NJB             ! Below wing
      JMXCOM(K3) = NJ3             ! Above forward part of wing
      JMXCOM(K4) = NJ3             ! Between wing and fin
      JMXCOM(K5) = JLBODY          ! Aft body

      NSUBJ(K1)  = 2
      NSUBJ(K5)  = 2

      JSUB(1,K1) = NJB
      JSUB(2,K1) = JLBODY
      JSUB(1,K2) = NJB
      JSUB(1,K3) = NJ3
      JSUB(1,K4) = NJ3
      JSUB(1,K5) = NJB
      JSUB(2,K5) = JLBODY

C     Check that there is room for all fuselage panels:

      IX = LASTXB + 1
      I1 = IX
      M  = 0

      DO L = K1, K5
         N = IMXCOM(L) * JMXCOM(L)
         M = N + M
         IXB(L) = I1
         I1 = I1 + N
      END DO

      CALL XBCHECK (LENXB, MXPATCH, K5, LASTXB, M, 1,
     >              IFUS, IWRIT, MOSTXB)

      LASTXB = LASTXB + M ! Upon return from FCTV3

C     Strategy for the subpatches:
C
C     >  Work from keel to crown, forward to aft.
C     >  Establish the corner points in (u,v) space.
C     >  Generate appropriate edges  in (u,v) space.
C     >  Generate interior (u,v)s via TFI.
C     >  Interpolate interior and edge (x,y,z)s (bilinear or bicubic).

C     Allocate (u,v)s for interpolating all patches (edges overlapping):

      ALLOCATE (UEVAL(NITOT,JLBODY), VEVAL(NITOT,JLBODY))

      DO J = 1, JLBODY
         UEVAL(1,J)      = ZERO
         UEVAL(NITOT,J)  = ONE
      END DO

      DO I = 1, NITOT
         VEVAL(I,1)      = ZERO
         VEVAL(I,JLBODY) = ONE
      END DO

C     Deal with the nose point or line for two subpatches:
C     ----------------------------------------------------

      I = IX

      IF (ZF(2,1) == ZF(1,1)) THEN ! Point

         DO J = 1, JLBODY
            XBASE(I) = XF(1)
            YBASE(I) = YF(1,1)
            ZBASE(I) = ZF(1,1)
            I = I + NIN
            VEVAL(1,J) = VBODY(J,1)
         END DO

      ELSE ! Singular line does NOT need the redistribution of FCTV2

         DO J = 1, JLBODY
            XBASE(I) = XF(1)
            YBASE(I) = YF(J,1)
            ZBASE(I) = ZF(J,1)
            I = I + NIN
            VEVAL(1,J) = VBODY(J,1)
         END DO

      END IF

C     First subpatch (on the nose at the keel):
C     -----------------------------------------

C     Upper right corner (wing LE):

      UEVAL(NIN,NJB) = UINTR(ILE1,1)
      VEVAL(NIN,NJB) = VINTR(ILE1,1)

C     Right-hand edge (u,v)s:

      I =  NIN - 1
      J =  JINTR(ILE1,1)
      R =  VINTR(ILE1,1) / REAL (NJB - 1)
      P = (UINTR(ILE1,1) - UBODY(J,I)) / (UBODY(J,NIN)  - UBODY(J,I))

      DO J = 1, NJB - 1
         UEVAL(NIN,J) = (ONE - P) * UBODY(J,I) + P * UBODY(J,NIN)
         VEVAL(NIN,J) = R * REAL (J - 1)
      END DO

C     Do the rest of the right-hand edge while we're at it:

      Q = VEVAL(NIN,NJB)
      R = (ONE - Q) / REAL (NJ3 - 1)

      DO J = NJB + 1, JLBODY
         UEVAL(NIN,J) = (ONE - P) * UBODY(J,I) + P * UBODY(J,NIN)
         VEVAL(NIN,J) = Q + R * REAL (J - NJB)
      END DO

C     Streamwise edge (u,v)s of nose subpatches:

      R = ONE / (XINTR(ILE1,1) - XF(1)) ! Straight (u,v) lines suffice

      DO I = 2, NIN - 1
         P = (XF(I) - XF(1)) * R
         UEVAL(I,1)   = UBODY(1,I)                      ! Keel
         UEVAL(I,NJB) = P * UEVAL(NIN,NJB)              ! J split at wing LE
         VEVAL(I,NJB) = P * Q + (ONE - P) * VEVAL(1,NJB)
         UEVAL(I,JLBODY) = UBODY(JLBODY,I)              ! Crown
      END DO

C     Panel K2 under the wing:
C     ------------------------

C     Close the trailing edge in (u,v) space:

      UEVAL(I2,NJB) = (UINTR(1,1) + UINTR(N1,1)) * HALF
      VEVAL(I2,NJB) = (VINTR(1,1) + VINTR(N1,1)) * HALF

C     Right edge (u,v)s:

      I =  IINTR(1,1)
      J =  JINTR(1,1)
      P = (UEVAL(I2,NJB) - UBODY(J,I)) / (UBODY(J,I+1)  - UBODY(J,I))
      R =  VEVAL(I2,NJB) / REAL (NJB - 1)

      DO J = 1, NJB - 1
         UEVAL(I2,J) = (ONE - P) * UBODY(J,I) + P * UBODY(J,I+1)
         VEVAL(I2,J) = R * REAL (J - 1)
      END DO

C     Keel (u,v)s + upper patch edge (= lower half of intersection).
C     Use "arc length" along the intersection (u,v)s in case U(LE) is not
C     at the same index as X(LE).

      IX = IXB(K2) ! Start of ample work-space for arc lengths

      CALL CHORDS2D (ILE1, UINTR(1,1), VINTR(1,1), .TRUE., P, XBASE(IX))

      P  = UEVAL(NIN,1)
      R  = UEVAL(I2,1) - P
      I1 = IX + ILE1 - 1
      J  = ILE1

      DO I = NIN + 1, I2 - 1
         I1 = I1 - 1
         UEVAL(I,1) = P + R * (ONE - XBASE(I1))
         J  = J - 1
         UEVAL(I,NJB) = UINTR(J,1)
         VEVAL(I,NJB) = VINTR(J,1)
      END DO


C     Fill the interior of patch K2 now so we can reuse its upper edge,
C     U/VEVAL(*,NJB), for the lower edge of patch K3.  Do the nose too:
C     -----------------------------------------------------------------

      DO J = 1, 2

         CALL TFI2D (NITOT, 1, NIN, JSUB(J-1,K1), JSUB(J,K1),
     >               UEVAL, VEVAL, XBASE(IX))
      END DO

      CALL TFI2D (NITOT, NIN, I2, 1, NJB, UEVAL, VEVAL, XBASE(IX)) ! K2 interior


C     Patch K3 above forward part of wing:
C     ------------------------------------

C     Lower edge of K3 & K4 is the upper part of the wing/body intersection:

      DO I = 1, ILE1 - 2 ! Don't clobber the closed TE
         UEVAL(NIN + I,NJB) = UINTR(ILE1 + I,1)
         VEVAL(NIN + I,NJB) = VINTR(ILE1 + I,1)
      END DO

C     Right-hand edge (u,v)s of lower subpatch of K3 (straight/uniform):

      UEVAL(I2S,JLBODY) = UINTR(ILE2,2) ! Fin LE
      R = ONE / REAL (NJ3 - 1)
      P = (UEVAL(I2S,JLBODY) - UEVAL(I2S,NJB)) * R
      Q = (ONE               - VEVAL(I2S,NJB)) * R

      DO J = NJB + 1, JLBODY - 1
         R = REAL (J - NJB)
         UEVAL(I2S,J) = UEVAL(I2S,NJB) + R * P
         VEVAL(I2S,J) = VEVAL(I2S,NJB) + R * Q
      END DO

C     Crown line (u,v)s match the forward wing intersection arc distribution:

      CALL CHORDS2D (NI3, UINTR(ILE1,1), VINTR(ILE1,1), .TRUE., P,
     >               XBASE(IX))
      P  = UINTR(ILE1,1)
      Q  = UEVAL(NIN,JLBODY)
      R  = UEVAL(I2S,JLBODY) - Q
      I1 = IX

      DO I = NIN + 1, I2S - 1
         I1 = I1 + 1
         UEVAL(I,JLBODY) = Q + R * XBASE(I1)
      END DO

C     Interior of K3:

      CALL TFI2D (NITOT, NIN, I2S, NJB, JLBODY, UEVAL, VEVAL, XBASE(IX))


C     Panel K4 below the fin/above the aft wing:
C     ------------------------------------------

C     Close the trailing edge:

      UEVAL(I2,JLBODY) = UINTR(N2,2) ! VEVAL(I2,JLBODY) = 1.0 already

C     Right edge (linear/uniform):

      R = ONE / REAL (NJ3 - 1)
      P = (UEVAL(I2,JLBODY) - UEVAL(I2,NJB)) * R
      Q = (ONE              - VEVAL(I2,NJB)) * R

      DO J = NJB + 1, JLBODY - 1
         R = REAL (J - NJB)
         UEVAL(I2,J) = UEVAL(I2,NJB) + R * P
         VEVAL(I2,J) = VEVAL(I2,NJB) + R * Q
      END DO

C     Upper edge is the upper edge of the fin/body intersection:

      DO I = 1, ILE2 - 2
         UEVAL(I2S+I,JLBODY) = UINTR(ILE2+I,2)
         VEVAL(I2S+I,JLBODY) = VINTR(ILE2+I,2)
      END DO

C     Interior of panel K4:

      CALL TFI2D (NITOT, I2S, I2, NJB, JLBODY, UEVAL, VEVAL, XBASE(IX))


C     Aft body panel in three subpatches:
C     -----------------------------------

C     Upper right corner:

      I1 = IINTR(1,1)
      J1 = JINTR(1,1)
      Q  = (VEVAL(I2,NJB)  - VBODY(J1,I1)) /
     >     (VBODY(J1+1,I1) - VBODY(J1,I1))
      P  = (ONE - Q) * VBODY(J1,  NFSTN) +   ! P = VEVAL(NITOT,NJB)
     >            Q  * VBODY(J1+1,NFSTN)

C     Right-hand edge:

      R = P / REAL (NJB - 1)

      DO J = 2, NJB
         VEVAL(NITOT,J) = R * REAL (J - 1)
      END DO

C     Keel line and edge aft of wing TE:

      IF (NIT == NFSTN - I1 + 1) THEN

         DO I = I2 + 1, NITOT - 1
            I1 = I1 + 1
            UEVAL(I,1)   = UBODY(1,I1)
            UEVAL(I,NJB) = (ONE - Q) * UBODY(J1,I1) + Q * UBODY(J1+1,I1)
            VEVAL(I,NJB) = (ONE - Q) * VBODY(J1,I1) + Q * VBODY(J1+1,I1)
         END DO

      ELSE ! Adjust the number of points

         XLINE(I1) = UEVAL(I2,1) ! Keel line first
         YLINE(I1) = ZERO

         DO I = I1 + 1, NFSTN
            XLINE(I) = UBODY(1,I)
            YLINE(I) = ZERO
         END DO

         ZLINE = ZERO ! For all I, to be safe

         CALL CHANGEN (I1, NFSTN, XLINE, YLINE, ZLINE,
     >                 I2, NITOT, UEVAL, VEVAL, ZLINE, METHOD)

         XLINE(I1) = UEVAL(I2,NJB) ! Now for wing TE --> end of body
         YLINE(I1) = VEVAL(I2,NJB)

         DO I = I1 + 1, NFSTN
            XLINE(I) = (ONE - Q) * UBODY(J1,I) + Q * UBODY(J1+1,I)
            YLINE(I) = (ONE - Q) * VBODY(J1,I) + Q * VBODY(J1+1,I)
         END DO

         CALL CHANGEN (I1, NFSTN, XLINE, YLINE, ZLINE,
     >                 I2, NITOT, UEVAL(1,NJB), VEVAL(1,NJB), ZLINE,
     >                 METHOD)
      END IF

C     Likewise for the upper aft-body subpatch edges:

C     Right-hand edge:

      R = (ONE - P) / REAL (NJ3 - 1)

      DO J = NJB + 1, JLBODY - 1
         VEVAL(NITOT,J) = P + R * REAL (J - NJB)
      END DO

C     Crown line:

      IF (NIT == NFSTN - I1 + 1) THEN

         DO I = I2 + 1, NITOT - 1
            I1 = I1 + 1
            UEVAL(I,JLBODY) = UBODY(JLBODY,I1)
         END DO

      ELSE ! Adjust the number of points

         XLINE(I1) = UEVAL(I2,JLBODY)
         YLINE(I1) = ONE

         DO I = I1 + 1, NFSTN
            XLINE(I) = UBODY(JLBODY,I)
            YLINE(I) = ONE
         END DO

         CALL CHANGEN (I1, NFSTN, XLINE, YLINE, ZLINE,
     >                 I2, NITOT, UEVAL(1,JLBODY), VEVAL(1,JLBODY),
     >                 ZLINE, METHOD)
      END IF

C     Interior (u,v)s of the aft-body subpatches:

      DO J = 1, 2

         CALL TFI2D (NITOT, I2, NITOT, JSUB(J-1,K5), JSUB(J,K5),
     >               UEVAL, VEVAL, XBASE(IX))
      END DO

!!!      write (60) 3
!!!      write (60) nitot, jlbody, n1, 1, ile2, 1
!!!      write (60) ueval, veval
!!!      write (60) uintr(1:n1,1), vintr(1:n1,1)
!!!      write (60) uintr(ile2:n2,2), vintr(ile2:n2,2)

C     Parametric bilinear/bicubic interpolation of (x,y,z)s for all panels:
C     ---------------------------------------------------------------------

      K = K1
      IMAT = 0
      JMAT = 0
      IXSHIFT = 0

      DO KASE = 1, 5

         SELECT CASE (KASE)

            CASE (1)
               I1 = 2
               L  = NIN
               J1 = 1
               M  = JLBODY
               IB = 2
               JB = 1
            CASE (2)
               I1 = NIN
               L  = I2
!              J1 = 1 (still)
               M  = NJB - 1
               IB = IINTR(ILE1,1)
               JB = 1
            CASE (3)
               L  = I2S
               J1 = NJB + 1
               M  = JLBODY
               IB = IINTR(ILE1,1)
               JB = JINTR(ILE1,1)
               IXSHIFT = NI3 ! Since J is not starting at the edge
            CASE (4)
               I1 = I2S
               L  = I2
               IB = IINTR(IU,1)
               JB = JINTR(IU,1)
               IXSHIFT = ILE2
            CASE (5)
               I1 = I2
               L  = NITOT
               J1 = 1
               IB = IINTR(1,1)
               JB = 1
               IXSHIFT = 0

         END SELECT

         IX = IXB(K) + IXSHIFT

         DO J = J1, M

            IF (K == K1) IX = IX + 1 ! Awkwardness from omitting nose pts.

            IE = IB
            JE = JB

            DO I = I1, L

               IF (BICUBIC) THEN

                  CALL PLBICUBE (0, JLBODY, NFSTN, 1, JLBODY, 1, NFSTN,
     >                           XBODY, YF, ZF, VBODY, UBODY,
     >                           VEVAL(I,J), UEVAL(I,J), JE, IE, 
     >                           JMAT, IMAT, SMATRIX, EPS, Q, P,
     >                           XBASE(IX), YBASE(IX), ZBASE(IX),
     >                           XLINE, YLINE, IER)
               ELSE

                  CALL PBILINT  (0, JLBODY, NFSTN, 1, JLBODY, 1, NFSTN,
     >                           XBODY, YF, ZF, VBODY, UBODY,
     >                           VEVAL(I,J), UEVAL(I,J), JE, IE,
     >                           EPS, Q, P,
     >                           XBASE(IX), YBASE(IX), ZBASE(IX),
     >                           XLINE, YLINE, IER)
               END IF

               IF (IER /= 0) THEN
                  IF (IWRIT > 0) THEN
                    WRITE (IWRIT,'(/,A,A,I2,A,4I5,/,A,2F13.5,A,3F13.5)')
     >               ' FCTV3: Marginal 2D interpolation. ',
     >               ' IER:', IER, '; I, J, IE, JE:', I, J, IE, JE,
     >               ' U, V:', UEVAL(I,J), VEVAL(I,J),
     >               '   X, Y, Z used:', XBASE(IX), YBASE(IX), ZBASE(IX)
                  END IF
               END IF

               IF (I == I1) THEN
                  IB = IE ! Ready for the next J
                  JB = JE
               END IF

               IX = IX + 1

            END DO

         END DO

         K = K + 1

      END DO ! Next panel


      DEALLOCATE (UEVAL, VEVAL)


C     Transcribe intersections (x,y,z)s and ensure common edge points:

      IX = IXB(K1) + NIN * NJB - 1   ! Wing LE pt. of nose panel
      XBASE(IX) = XINTR(ILE1,1)
      YBASE(IX) = YINTR(ILE1,1)
      ZBASE(IX) = ZINTR(ILE1,1)

      IX = IXB(K4) - 1               ! Fin LE pt. of K3 panel
      XBASE(IX) = XINTR(ILE2,2)
      YBASE(IX) = YINTR(ILE2,2)
      ZBASE(IX) = ZINTR(ILE2,2)

      IX = IXB(K5) + NIT * (NJB - 1) ! Wing TE pt. from aft-body panel
      I1 = IXB(K3) - 1               ! Corresp. pt. in under-wing panel
      XBASE(I1) = XBASE(IX)
      YBASE(I1) = YBASE(IX)
      ZBASE(I1) = ZBASE(IX)

      I1 = IXB(K4) + ILE2 - 1        ! Corresp. pt. of panel between wing & fin
      XBASE(I1) = XBASE(IX)
      YBASE(I1) = YBASE(IX)
      ZBASE(I1) = ZBASE(IX)

      IX = IXB(K5) + NIT * (JLBODY - 1)  ! Fin TE pt. from aft-body panel
      I1 = IXB(K5) - 1               ! Corresp. pt. in panel between wing & fin
      XBASE(I1) = XBASE(IX)
      YBASE(I1) = YBASE(IX)
      ZBASE(I1) = ZBASE(IX)

      IX = IXB(K3) - 1               ! Lower wing intersection (except TE)
      DO I = 2, ILE1
         IX = IX - 1
         XBASE(IX) = XINTR(I,1)
         YBASE(IX) = YINTR(I,1)
         ZBASE(IX) = ZINTR(I,1)
      END DO

      IX = IXB(K3)                   ! Upper wing intersection, forward
      DO I = ILE1, IU
         XBASE(IX) = XINTR(I,1)
         YBASE(IX) = YINTR(I,1)
         ZBASE(IX) = ZINTR(I,1)
         IX = IX + 1
      END DO

      IX = IXB(K4)                   ! Upper wing int., aft (except TE)
      DO I = IU, N1 - 1
         XBASE(IX) = XINTR(I,1)
         YBASE(IX) = YINTR(I,1)
         ZBASE(IX) = ZINTR(I,1)
         IX = IX + 1
      END DO

      IX = IXB(K5) - ILE2            ! Lower fin intersection (except TE)
      DO I = ILE2, N2 - 1
         XBASE(IX) = XINTR(I,2)
         YBASE(IX) = YINTR(I,2)
         ZBASE(IX) = ZINTR(I,2)
         IX = IX + 1
      END DO

      
C     Fix up extreme skewing at blunt wing leading edge?

      NFIX = NFIXL(1) + NFIXR(1) + 1

      IF (NFIX > 2) THEN ! Patches K1, K1+1, K1+2 are partially revised

         CALL FIXROOT (MXPATCH, LENXB, IXB, IMXCOM, JMXCOM, K1,
     >                 JLBODY, NFIXL(1), NFIXR(1), NFIX, NFIXA(1), 
     >                 NFIXB(1), EPS, XBASE, YBASE, ZBASE)
      END IF


C     Finally:  split the wing patches implicitly at the point corresp. to I2S:

      NSUBI(1:2) = 2
      ISUB(1,1)  = ILE2
      ISUB(2,1)  = ILE1
      ISUB(1,2)  = NI3
      ISUB(2,2)  = ILE1

      KPATCH = K5 + 1 ! Next available patch
      FAIL   = .FALSE.

      END SUBROUTINE FCTV3

C***********************************************************************
C
      SUBROUTINE FDORSAL (L1, L2, L3, IDIM, JDIM, METHOD, X1, Y1, Z1,
     >                    X2, Y2, Z2, X3, Y3, Z3, NFIN,
     >                    XFIN, YFIN, ZFIN, IWRIT)
C
C     FDORSAL redistributes the three standard fuselage (sub)patches
C     associated with a horizontal tail in order to accommodate an
C     accompanying dorsal fin (represented by its fin/body intersection).
C     The point counts of the fin and tail are required to match.
C     Since the point counts do not change, the input patches are
C     overwritten with the desired results.
C
C     06/11/99  David Saunders  Initial implementation.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   L1, L2, L3,           ! Patches forward, above, and aft of the tail
     >   IDIM(L3), JDIM(L3)    ! Patch dimensions

      CHARACTER, INTENT (IN) ::
     >   METHOD * 1            ! 1-D method for fuselage barrels also controls
                               ! use of PLBICUT/PLBICUBE or BILINCUT/PBILINT

      REAL, INTENT (INOUT), DIMENSION (IDIM(L1), JDIM(L1)) ::
     >   X1, Y1, Z1            ! Patch left of tail

      REAL, INTENT (INOUT), DIMENSION (IDIM(L2), JDIM(L2)) ::
     >   X2, Y2, Z2            ! Patch above tail

      REAL, INTENT (INOUT), DIMENSION (IDIM(L3), JDIM(L3)) ::
     >   X3, Y3, Z3            ! Patch aft of tail

      INTEGER, INTENT (IN) ::
     >   NFIN                  ! NFIN = IDIM(L2) is assumed

      REAL, INTENT (IN), DIMENSION (NFIN) ::
     >   XFIN, YFIN, ZFIN      ! Fin/body intersection pts.

      INTEGER, INTENT (IN) ::
     >   IWRIT                 ! For diagnostics

C     Local constants:

      REAL,      PARAMETER :: ONE = 1., ZERO = 0.
      LOGICAL,   PARAMETER :: NEW = .TRUE.
      CHARACTER, PARAMETER :: LOOSE = 'B'

C     Local variables:

      INTEGER
     >   I, I0, I2, IER, IP, ISIZE, IS, ISMAT, J, J0, J1, JP, JS,
     >   JSIZE, JSMAT, LI, LUNERR, MI, NI
      REAL
     >   DRANGE, DU, EPS, P, Q, R, SMATRIX(4,4,3), UCUT, VCUT,
     >   UVRANGE(4), XINT, YINT, ZINT, XYZCUT(3)
      REAL, ALLOCATABLE, DIMENSION (:) ::
     >   EDGES
      REAL, ALLOCATABLE, DIMENSION (:,:) ::
     >   X, Y, Z, U, V, UNEW, VNEW
      LOGICAL
     >   BICUBIC

C     Storage:

      DATA UVRANGE /ZERO, ONE, ZERO, ONE/ ! See PLBICUT

C     Execution:

      EPS = MAX (5.* EPSILON (EPS), 1.E-7) ! For PLBICUT/RIPPLE2D
      BICUBIC = METHOD == 'B'

C     We need to gather the three (sub)patches as one in order to redistribute.
C     Set up the workspace:

      NI    = IDIM(L1)
      MI    = NI + NFIN - 1
      LI    = IDIM(L3)
      ISIZE = MI + LI - 1
      JSIZE = JDIM(L2)

      ALLOCATE (X(ISIZE,JSIZE), Y(ISIZE,JSIZE), Z(ISIZE,JSIZE))

      J1 = JDIM(L1) - JSIZE
      JP = J1

      DO J = 1, JSIZE
         JP = JP + 1
         DO I = 1, NI
            X(I,J) = X1(I,JP)
            Y(I,J) = Y1(I,JP)
            Z(I,J) = Z1(I,JP)
         END DO
      END DO

      DO J = 1, JSIZE
         IP = NI
         DO I = 2, NFIN
            IP = IP + 1
            X(IP,J) = X2(I,J)
            Y(IP,J) = Y2(I,J)
            Z(IP,J) = Z2(I,J)
         END DO
      END DO

      JP = J1
      DO J = 1, JSIZE
         JP = JP + 1
         IP = MI
         DO I = 2, LI
            IP = IP + 1
            X(IP,J) = X3(I,JP)
            Y(IP,J) = Y3(I,JP)
            Z(IP,J) = Z3(I,JP)
         END DO
      END DO

      ALLOCATE (U(ISIZE,JSIZE), UNEW(ISIZE,JSIZE),
     >          V(ISIZE,JSIZE), VNEW(ISIZE,JSIZE))

      U(1,1) = ZERO ! -999. suppresses normalization

      CALL PARAM2D (ISIZE, JSIZE, 1, ISIZE, 1, JSIZE, X, Y, Z, U, V)

C     Three of the edges in (u,v) space stay the same:

      DO I = 1, ISIZE
         UNEW(I,1) = U(I,1)
         VNEW(I,1) = ZERO
      END DO

      DO J = 2, JSIZE
         UNEW(1,J) = ZERO
         VNEW(1,J) = V(1,J)
         UNEW(ISIZE,J) = ONE
         VNEW(ISIZE,J) = V(ISIZE,J)
      END DO

C     Locate the fin/body intersection points in (u,v) space:

      DRANGE = X(ISIZE,1) - X(1,1)
      JP = JSIZE
      JS = JSIZE - 1
      IP = NI
      IS = NI

      CALL INTERVAL (ISIZE, X(1,JP), XFIN(1), ONE, IS) ! Fin LE pt.

      R = (XFIN(1) - X(IS,JP)) / (X(IS+1,JP) - X(IS,JP))
      UCUT = (ONE - R) * U(IS,JP) + R * U(IS+1,JP)
      UNEW(IP,JP) = UCUT
      VNEW(IP,JP) = ONE
      VCUT = ONE
      J = JSIZE / 2 ! Should be plenty big enough
      ISMAT = 0
      JSMAT = 0
      LUNERR = -ABS (IWRIT)

      DO I = 1, NFIN ! Highly nonuniform spacing affects PLBICUBE, so redo I=1

         XYZCUT(1) = XFIN(I)
         XYZCUT(3) = ZFIN(I)

         IF (BICUBIC) THEN

            CALL PLBICUT (3, ISIZE, JSIZE, 1, ISIZE, J, JSIZE,
     >                    X, Y, Z, U, V, IS, JS, ISMAT, JSMAT,
     >                    SMATRIX, EPS, DRANGE, XYZCUT, UVRANGE,
     >                    UCUT, VCUT, P, Q, LUNERR, IER)
         ELSE

            CALL BILINCUT (3, ISIZE, JSIZE, 1, ISIZE, J, JSIZE,
     >                     X, Y, Z, U, V, IS, JS,
     >                     EPS, DRANGE, XYZCUT, UVRANGE,
     >                     UCUT, VCUT, P, Q, LUNERR, IER)
         END IF

         IF (IER > 2) THEN
            IF (IWRIT > 0) THEN
               WRITE (IWRIT, '(/, A, I2, /, A, I4, /, A, 1P, 3E15.6)')
     >            ' FDORSAL: PLBICUT/BILINCUT error =', IER,
     >            ' Fin/body intersection pt. =', I,
     >            ' Intersection pt. (x,y,z)  =', XYZCUT
            END IF
            CALL SYNCH
            STOP
         END IF

         UNEW(IP,JP) = UCUT
         VNEW(IP,JP) = VCUT
         IP = IP + 1

      END DO

C     Match the relative water-line distribution fore & aft of the fin:

      R = UNEW(NI,JSIZE) / UNEW(NI,1)

      DO I = 2, NI - 1
         UNEW(I,JSIZE) = R * UNEW(I,1)
         VNEW(I,JSIZE) = ONE
      END DO

      R = (ONE - UNEW(MI,JSIZE)) / (ONE - UNEW(MI,1))

      DO I = 2, LI - 1
         UNEW(IP,JP) = UNEW(MI,JP) + R * (UNEW(IP,1) - UNEW(MI,1))
         VNEW(IP,JP) = ONE
         IP = IP + 1
      END DO

C     New interior (u,v)s:

      ALLOCATE (EDGES(2*(ISIZE+JSIZE)))

      CALL TFI2D (ISIZE, 1, ISIZE, 1, JSIZE, UNEW, VNEW, EDGES)

      DEALLOCATE (EDGES)

C     Redistribute the interior points:

      I2 = 1
      I0 = 0
      J0 = J1

      DO I = 2, ISIZE - 1

         IF (I == NI) THEN
            I0 = NI - 1
            J0 = 0
         ELSE IF (I == MI) THEN
            I0 = MI - 1
            J0 = J1
         END IF

         IP = I - I0
         IS = I2
         JS = 1

         DO J = 2, JSIZE - 1

            IF (BICUBIC) THEN

               CALL PLBICUBE (0, ISIZE, JSIZE, 1, ISIZE, 1, JSIZE,
     >                        X, Y, Z, U, V, UNEW(I,J), VNEW(I,J),
     >                        IS, JS, ISMAT, JSMAT, SMATRIX, EPS, P, Q,
     >                        XINT, YINT, ZINT, XYZCUT, XYZCUT, IER)
            ELSE

               CALL PBILINT (0, ISIZE, JSIZE, 1, ISIZE, 1, JSIZE,
     >                       X, Y, Z, U, V, UNEW(I,J), VNEW(I,J),
     >                       IS, JS, EPS, P, Q, XINT, YINT, ZINT,
     >                       XYZCUT, XYZCUT, IER)
            END IF

            IF (IER /= 0) THEN
               IF (IWRIT > 0)
     >            WRITE (IWRIT, '(/,A,A,I2,A,4I5,/,A,2F13.5,A,3F13.5)')
     >               ' FDORSAL: Marginal 2D interpolation. ',
     >               ' IER:', IER, '; I, J, IS, JS:', I, J, IS, JS,
     >               ' U, V:', UNEW(I,J), VNEW(I,J),
     >               '    X, Y, Z used:', XINT, YINT, ZINT
            END IF

            IF (J == 2) I2 = IS ! Crown pts. can be offset from water-line pts.

            JP = J + J0

            IF (I < NI) THEN
               X1(IP,JP) = XINT
               Y1(IP,JP) = YINT
               Z1(IP,JP) = ZINT
            ELSE IF (I < MI) THEN
               X2(IP,JP) = XINT
               Y2(IP,JP) = YINT
               Z2(IP,JP) = ZINT
            ELSE
               X3(IP,JP) = XINT
               Y3(IP,JP) = YINT
               Z3(IP,JP) = ZINT
            END IF

         END DO

      END DO

C     New crown points:

      JP = JDIM(L1)

      CALL LCSFIT (ISIZE, U(1,JSIZE), X(1,JSIZE), NEW, LOOSE,
     >             NI - 2, UNEW(2,JSIZE), X1(2,JP), Z(1,JSIZE))

      CALL LCSFIT (ISIZE, U(1,JSIZE), Y(1,JSIZE), NEW, LOOSE,
     >             NI - 2, UNEW(2,JSIZE), Y1(2,JP), Z(1,JSIZE))

      X1(NI,JP) = XFIN(1)
      Y1(NI,JP) = YFIN(1)
      Z1(NI,JP) = ZFIN(1)

      DO I = 1, NFIN
         X2(I,JSIZE) = XFIN(I)
         Y2(I,JSIZE) = YFIN(I)
         Z2(I,JSIZE) = ZFIN(I)
      END DO

      X3(1,JP) = XFIN(NFIN)
      Y3(1,JP) = YFIN(NFIN)
      Z3(1,JP) = ZFIN(NFIN)

      CALL LCSFIT (ISIZE, U(1,JSIZE), X(1,JSIZE), NEW, LOOSE,
     >             LI - 2, UNEW(MI+1,JSIZE), X3(2,JP), Z(1,JSIZE))

      CALL LCSFIT (ISIZE, U(1,JSIZE), Y(1,JSIZE), NEW, LOOSE,
     >             LI - 2, UNEW(MI+1,JSIZE), Y3(2,JP), Z(1,JSIZE))

      DEALLOCATE (X, Y, Z, U, V, UNEW, VNEW)

C     Duplicate two edges:

      IP = NFIN
      JP = J1 + 2

      DO J = 2, JSIZE - 1
         X1(NI,JP) = X2(1,J)
         Y1(NI,JP) = Y2(1,J)
         Z1(NI,JP) = Z2(1,J)
         X2(IP,J)  = X3(1,JP)
         Y2(IP,J)  = Y3(1,JP)
         Z2(IP,J)  = Z3(1,JP)
         JP = JP + 1
      END DO

      END SUBROUTINE FDORSAL

C***********************************************************************
C
      SUBROUTINE FINSERT (IW, MXPATCH, KPATCH, NFSTN, JLBODY,
     >                    NMIN, IFUS, METHOD, XF, YF, ZF, LOSTINTR,
     >                    NINTR, XINTR, YINTR, ZINTR, IINTR, JINTR,
     >                    ILE, ISTN1, ISTN2, LENXB,
     >                    XBASE, YBASE, ZBASE, IXB,
     >                    XTEMP, YTEMP, ZTEMP, SFINL,
     >                    IMXCOM, JMXCOM, IWRIT, FAIL)
C
C     FINSERT panels the two body regions above and below a wing/body
C     intersection.  It inserts a new body section at the leading and
C     trailing edges of the intersection.
C
C     This version expects all patch dimensions to be input in I/JMXCOM,
C     and the intersection trailing edge to be (temporarily) closed.
C     It also transcribes the LE & TE barrels to the neighboring patches.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   IW, MXPATCH, KPATCH
      INTEGER, INTENT (IN) ::
     >   NFSTN, JLBODY, NMIN, IFUS
      CHARACTER, INTENT (IN) ::
     >   METHOD * 1
      REAL, INTENT (IN) ::
     >   XF(NFSTN), YF(JLBODY,NFSTN), ZF(JLBODY,NFSTN)
      LOGICAL, INTENT (IN) ::
     >   LOSTINTR
      INTEGER, INTENT (IN) ::
     >   NINTR
      REAL, INTENT (IN) ::
     >   XINTR(NINTR), YINTR(NINTR), ZINTR(NINTR)
      INTEGER, INTENT (IN) ::
     >   IINTR(NINTR), JINTR(NINTR)
      INTEGER, INTENT (IN) ::
     >   ILE, ISTN1, ISTN2, LENXB
      REAL, INTENT (INOUT), DIMENSION (LENXB) ::
     >   XBASE, YBASE, ZBASE
      INTEGER, INTENT (INOUT) ::
     >   IXB (MXPATCH)
      REAL, INTENT (OUT), DIMENSION (JLBODY) ::
     >   XTEMP, YTEMP, ZTEMP, SFINL
      INTEGER, INTENT (INOUT), DIMENSION (MXPATCH) ::
     >   IMXCOM, JMXCOM
      INTEGER, INTENT (IN) ::
     >   IWRIT
      LOGICAL, INTENT (OUT) ::
     >   FAIL

C     Local constants:

      REAL,    PARAMETER :: ONE = 1.0
      LOGICAL, PARAMETER :: CLOSED = .FALSE.

C     Local variables:

      INTEGER
     >   I, I1, I2, IER, II, ILEFT, J, JLEFT,
     >   K1, K2, K3, K4, NJBODY, NJU, NUMI, NUMJ
      REAL
     >   DS, R, STOTAL, TOL, XKEY, YKEY, ZKEY
      LOGICAL
     >   NEW

C     Execution:

      K1 = KPATCH ! Forward of intersection
      K2 = K1 + 1 ! Below       "   "
      K3 = K1 + 2 ! Above       "   "
      K4 = K1 + 3 ! Aft of      "   "
      NJBODY = JLBODY - 1

C     Construct the patch below the wing:
C     -----------------------------------

      I1 = IXB(K2)
      NUMJ = JMXCOM(K2)

      DS = ONE / REAL (NUMJ - 1)
      DO J = 1, NUMJ - 1       ! Normalized uniform distribution
         SFINL(J) = DS * REAL (J - 1)
      END DO
      SFINL(NUMJ) = ONE

      ILEFT = IINTR(ILE) ! Body cell indices for true LE
      JLEFT = JINTR(ILE)

      DO I = 1, ILE

         II   = ILE - I + 1
         XKEY = XINTR(II)
         YKEY = YINTR(II)
         ZKEY = ZINTR(II)

         IF (LOSTINTR) THEN ! Search, starting with previous pointers

            CALL BODYSRCH (NFSTN, JLBODY, ISTN1, ISTN2, 1, NJBODY,
     >                     XF, YF, ZF, XKEY, YKEY, ZKEY, ILEFT, JLEFT,
     >                     IER)
            IF (IER /= 0) THEN
               IF (IWRIT > 0) THEN
                  WRITE (IWRIT, '(/,A,I3,/,A,I4,1P,3X,A,3E15.6)')
     >               ' FINSERT: Trouble for lower wing surface', IW,
     >               ' II:', II, 'X/Y/ZKEY:', XKEY, YKEY, ZKEY
                  FAIL = .TRUE.
                  GO TO 999
               END IF
            END IF

         ELSE

            ILEFT = IINTR(II)
            JLEFT = JINTR(II)

         END IF

C        Construct a (lower) body section at the intersection point:

         R = (XKEY - XF(ILEFT)) / (XF(ILEFT+1) - XF(ILEFT))

         DO J = 1, JLEFT
            YTEMP(J) = (ONE - R) * YF(J,ILEFT) + R * YF(J,ILEFT+1)
            ZTEMP(J) = (ONE - R) * ZF(J,ILEFT) + R * ZF(J,ILEFT+1)
         END DO

         J = JLEFT + 1
         YTEMP(J) = YKEY
         ZTEMP(J) = ZKEY

C        Respline the lower barrel for this I:

         CALL FBARREL (I, ILE, NUMJ, J, METHOD,
     >                 XKEY, YTEMP, ZTEMP, NUMJ, SFINL, XTEMP,
     >                 XBASE(I1), YBASE(I1), ZBASE(I1))
      END DO

C     Construct the patch above the wing:
C     -----------------------------------

      I1   = IXB(K3)
      NUMI = IMXCOM(K3)
      NUMJ = JMXCOM(K3)

      DS = ONE / REAL (NUMJ - 1)
      DO J = 1, NUMJ
         SFINL(J) = DS * REAL (J - 1)
      END DO
      SFINL(NUMJ) = ONE

      ILEFT = IINTR(ILE)
      JLEFT = JINTR(ILE)

      DO I = 1, NUMI

         II   = ILE + I - 1
         XKEY = XINTR(II)
         YKEY = YINTR(II)
         ZKEY = ZINTR(II)

         IF (LOSTINTR) THEN ! Search, starting with previous pointers

            CALL BODYSRCH (NFSTN, JLBODY, ISTN1, ISTN2, 1, NJBODY,
     >                     XF, YF, ZF, XKEY, YKEY, ZKEY, ILEFT, JLEFT,
     >                     IER)
            IF (IER /= 0) THEN
               IF (IWRIT > 0) THEN
                  WRITE (IWRIT, '(/,A,I3,/,A,I4,1P,3X,A,3E15.6)')
     >               ' FINSERT: Trouble for upper wing surface', IW,
     >               ' II:', II, 'X/Y/ZKEY:', XKEY, YKEY, ZKEY
                  FAIL = .TRUE.
                  GO TO 999
               END IF
            END IF

         ELSE

            ILEFT = IINTR(II)
            JLEFT = JINTR(II)

         END IF

         R = (XKEY - XF(ILEFT)) / (XF(ILEFT+1) - XF(ILEFT))

         YTEMP(JLEFT) = YKEY
         ZTEMP(JLEFT) = ZKEY

         DO J = JLEFT + 1, NJBODY
            YTEMP(J) = (ONE - R) * YF(J,ILEFT) + R * YF(J,ILEFT+1)
            ZTEMP(J) = (ONE - R) * ZF(J,ILEFT) + R * ZF(J,ILEFT+1)
         END DO

         J = JLEFT
         NJU = NJBODY - J + 1

C        Respline the upper barrel:

         CALL FBARREL (I, NUMI, NUMJ, NJU, METHOD,
     >                 XKEY, YTEMP(J), ZTEMP(J), NUMJ, SFINL, XTEMP,
     >                 XBASE(I1), YBASE(I1), ZBASE(I1))
      END DO

C     Load the end of the forward patch:

      J  = JMXCOM(K2)
      I2 = IMXCOM(K1)

      CALL XBCOPY2 (LENXB, XBASE, YBASE, ZBASE, IXB, IMXCOM,
     >              K2, 1, 1, 1, J, K1, I2, 1)

      CALL XBCOPY2 (LENXB, XBASE, YBASE, ZBASE, IXB, IMXCOM,
     >              K3, 1, 1, 2, JMXCOM(K3), K1, I2, J + 1)

C     Load the start of the aft patch:

      I2 = IMXCOM(K2)

      CALL XBCOPY2 (LENXB, XBASE, YBASE, ZBASE, IXB, IMXCOM,
     >              K2, I2, I2, 1, J, K4, 1, 1)

      I2 = IMXCOM(K3)

      CALL XBCOPY2 (LENXB, XBASE, YBASE, ZBASE, IXB, IMXCOM,
     >              K3, I2, I2, 2, JMXCOM(K3), K4, 1, J + 1)

      FAIL = .FALSE.

  999 RETURN

      END SUBROUTINE FINSERT

C***********************************************************************
C
      SUBROUTINE FIXROOT (MXPATCH, LENXB, IXB, IMXCOM, JMXCOM, L,
     >                    JLBODY, NFIXL, NFIXR, NFIX, NFIXA, NFIXB,
     >                    EPS, XBASE, YBASE, ZBASE)
C
C     Redistribute parts of fuselage patches around the leading edge of a
C     wing-type intersection to overcome cell skewing caused by bluntness.
C     Patches L, L+1, L+2 are partially modified.
C
C     NFIXL & NFIXR lines left & right of LE, including LE, are revised.
C     NFIXA & -B lines above & below LE, inclusive, have v splined (axial);
C     beyond that, the vs are uniform between LE - NFIXL and LE + NFIXR.
C     NFIX is NFIXL + NFIXR + 1.
C
C     01/30/98  DAS  In-line implementation in FSURF.
C     02/02/98   "   Modularized so NFIXL/R can be user inputs.  Added
C                    NFIXA/B for more control in the vertical direction.
C     08/28/98   "   Redistributed the circumferential lines uniformly
C                    to smooth out some irregularities from normalized
C                    (u,v) manipulations.
C     09/14/98   "   Redistribution is confined to NFIXA/B region now.
C     08/10/00   "   Uniform V along line forward of the leading edge.
C     08/23/03   "   Upper portion forward of root was being affected
C                    by missing update of DV2 and VFIX(NJ,I).
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   MXPATCH, LENXB
      INTEGER, INTENT (IN), DIMENSION (MXPATCH) ::
     >   IXB, IMXCOM, JMXCOM
      INTEGER, INTENT (IN) ::
     >   L, JLBODY, NFIXL, NFIXR, NFIX, NFIXA, NFIXB
      REAL, INTENT (IN) ::
     >   EPS
      REAL, INTENT (INOUT), DIMENSION (LENXB) ::
     >   XBASE, YBASE, ZBASE

C     Local constants:

      INTEGER,   PARAMETER :: LINUV = 2 ! Linear in the v direction
      REAL,      PARAMETER :: DERIVS = -999., ONE = 1., ZERO = 0.
      LOGICAL,   PARAMETER :: CLOSED = .FALSE., NORM = .FALSE.
      CHARACTER, PARAMETER :: LOOSE * 1 = 'B', TIGHT * 1 = 'M'

C     Local variables:

      INTEGER
     >   I, I1, I2, IB, IC, IDIM, IDIM1, IDIM2, IER, IP, J, JE, JEVAL,
     >   JP, MJ, NEVAL, NJ
      REAL
     >   DU, DV, DV2, RN, RNJ, RMJ, TEVAL, TOTAL, UDAT(4), VDAT(4)
      REAL, DIMENSION (JLBODY) ::
     >   T, X, Y, Z
      REAL, DIMENSION (JLBODY,NFIX) ::
     >   UBOD, VBOD, XBOD, YBOD, ZBOD, UFIX, VFIX, XFIX, YFIX, ZFIX
      LOGICAL
     >   NEW

C     Execution:

C     Gather the indicated parts of patches L, L+1 (below wing), & L+2 (above).
C     Address (B(i,j)) = Address (B(1,1)) + (j - 1)idim + i - 1

      IDIM  = IMXCOM(L)
      IDIM1 = IMXCOM(L+1)
      IDIM2 = IMXCOM(L+2)
      I     = IXB(L) + (IDIM - NFIXL) - 1
      NJ    = JMXCOM(L+1)
      MJ    = JMXCOM(L+2)
      RNJ   = ONE / REAL (NFIXB)
      RMJ   = ONE / REAL (NFIXA)

      DO IP = 1, NFIXL
         IB = I ! (i,j) = (IDIM - NFIXL, 1) on 1st pass
         DO J = 1, JLBODY
            XBOD(J,IP) = XBASE(IB)
            YBOD(J,IP) = YBASE(IB)
            ZBOD(J,IP) = ZBASE(IB)
            IB = IB + IDIM
         END DO
         I  = I + 1
      END DO

      IP = NFIXL
      J  = NJ + 1
      JP = NJ + MJ - 1
      IB = IXB(L+1)
      IC = IXB(L+2)

      DO I = 1, NFIXR + 1
         IP = IP + 1
         I1 = IB
         DO J = 1, NJ
            XBOD(J,IP) = XBASE(I1) ! (I,J,L+1)
            YBOD(J,IP) = YBASE(I1)
            ZBOD(J,IP) = ZBASE(I1)
            I1 = I1 + IDIM1
         END DO
         IB = IB + 1

         I2 = IC
         DO J = NJ + 1, JP ! Skip J = 1 of L2 till after lower portion is done
            I2 = I2 + IDIM2
            XBOD(J,IP) = XBASE(I2) ! (I,2,L+2) on 1st pass
            YBOD(J,IP) = YBASE(I2)
            ZBOD(J,IP) = ZBASE(I2)
         END DO
         IC = IC + 1

      END DO

      UBOD(1,1) = ZERO ! -999. suppresses normalization

      CALL PARAM2D (JLBODY, NFIX, 1, NJ, 1, NFIX, ! Lower portion
     >              XBOD, YBOD, ZBOD, UBOD, VBOD)

!!!   write (61) 2
!!!   write (61) nj, nfix, mj, nfix
!!!   write (61) vbod(1:nj,1:nfix), ubod(1:nj,1:nfix)

C     Copy the edges.  Avoiding the interior is not worth the trouble.

      DO I = 1, NFIX
         DO J = 1, NJ
            UFIX(J,I) = UBOD(J,I)
            VFIX(J,I) = VBOD(J,I)
         END DO
      END DO

C     Make the v (axial) distribution uniform for most Js:

      JP = NJ - NFIXB

      DV = ONE / REAL (NFIX - 1)
      DO I = 2, NFIX - 1
         VFIX(1:JP,I) = DV * REAL (I - 1)
      END DO

C     Complete the partial lines via 4-pt. splines (v as a function of u):

      NEVAL = NFIXB - 1
      DV2   = VFIX(NJ,NFIXL + 1) / REAL (NFIXL)
      NEW   = .TRUE.

      DO I = 2, NFIX - 1

         UDAT(1:3) = UFIX(JP-2:JP,I)
         VDAT(1:3) = VFIX(JP-2:JP,I)
         UDAT(4)   = ONE

         IF (I <= NFIXL) THEN ! In case body sections are oddly spaced
            VDAT(4) = DV2 * REAL (I - 1)
            VFIX(NJ,I) = VDAT(4)
         ELSE
            VDAT(4) = VFIX(NJ,I)
         END IF

         DU = (ONE - UDAT(3)) * RNJ

         DO J = JP + 1, NJ - 1
            UFIX(J,I) = UFIX(JP,I) + REAL (J - JP) * DU
         END DO

         CALL LCSFIT (4, UDAT, VDAT, NEW, TIGHT, NEVAL, UFIX(JP+1,I),
     >                VFIX(JP+1,I), T)
      END DO

!!!   write (62) 2
!!!   write (62) nj, nfix, mj, nfix
!!!   write (62) vfix(1:nj,1:nfix), ufix(1:nj,1:nfix)

C     Interpolate (x,y,z) at each new interior (u,v) + forward of LE:

      JE = NJ - NFIXB
      IP = 1

      DO I = 2, NFIX - 1

         JP = 1
         DO J = 1, NJ

            CALL PLINCUB (JLBODY, NFIX, 1, NJ, 1, NFIX,
     >                    XBOD, YBOD, ZBOD, UBOD, VBOD,
     >                    UFIX(J,I), VFIX(J,I), LINUV, TIGHT, CLOSED,
     >                    JP, IP, EPS, X(J), Y(J), Z(J), IER)
         END DO

         IF (I > NFIXL) THEN
            X(NJ) = XBOD(NJ,I)
            Y(NJ) = YBOD(NJ,I)
            Z(NJ) = ZBOD(NJ,I)
         END IF

C        Redistribute uniformly in the circumferential direction:

         DO J = 1, JE
            XFIX(J,I) = X(J)
            YFIX(J,I) = Y(J)
            ZFIX(J,I) = Z(J)
         END DO

         CALL CHORDS3D (NJ, X, Y, Z, NORM, TOTAL, T)

         JEVAL = JE
         NEW = .TRUE.
         DU  = (TOTAL - T(JE)) * RNJ

         DO J = JE + 1, NJ - 1
            TEVAL = T(JE) + DU * REAL (J - JE)

            CALL PLSCRV3D (NJ, X, Y, Z, T, LOOSE, NEW, CLOSED, TEVAL,
     >                     JEVAL, XFIX(J,I), YFIX(J,I), ZFIX(J,I),
     >                     DERIVS)
            NEW = .FALSE.
         END DO

         XFIX(NJ,I) = X(NJ)
         YFIX(NJ,I) = Y(NJ)
         ZFIX(NJ,I) = Z(NJ)

      END DO ! Next I

      I = IXB(L) + IMXCOM(L) - NFIXL
      DO IP = 2, NFIXL + 1
         IB = I
         DO J = 1, NJ
            XBASE(IB) = XFIX(J,IP) ! (i,j) = (IDIM-NFIXL+1, 1) on 1st pass
            YBASE(IB) = YFIX(J,IP)
            ZBASE(IB) = ZFIX(J,IP)
            IB = IB + IDIM
         END DO
         I = I + 1
      END DO

      IB = IXB(L+1)
      IP = NFIXL
      DO I = 1, NFIXR
         IP = IP + 1
         I1 = IB
         DO J = 1, NJ - 1
            XBASE(I1) = XFIX(J,IP) ! xbase(i,j,l+1)
            YBASE(I1) = YFIX(J,IP)
            ZBASE(I1) = ZFIX(J,IP)
            I1 = I1 + IDIM1
         END DO
         IB = IB + 1
      END DO

C     Above the leading edge:
C     -----------------------

      I2 = IXB(L+2)
      IP = NFIXL + 1
      DO I = 1, NFIXR ! Finish the copy done mostly above
         I2 = I2 + 1
         IP = IP + 1
         XBOD(NJ,IP) = XBASE(I2)
         YBOD(NJ,IP) = YBASE(I2)
         ZBOD(NJ,IP) = ZBASE(I2)
      END DO

      UBOD(NJ,1) = ZERO

      CALL PARAM2D (JLBODY, NFIX, NJ, JLBODY, 1, NFIX,
     >              XBOD, YBOD, ZBOD, UBOD, VBOD)

!!!   write (61) vbod(nj:jlbody,1:nfix), ubod(nj:jlbody,1:nfix)

      DO I = 1, NFIX
         DO J = NJ, JLBODY
            UFIX(J,I) = UBOD(J,I) ! Edges needed
            VFIX(J,I) = VBOD(J,I)
         END DO
      END DO

      JP = NJ + NFIXA
      DO I = 2, NFIX - 1
         VFIX(JP:JLBODY,I) = DV * REAL (I-1)
      END DO

      NEVAL = NFIXA - 1
      DV2   = VFIX(NJ,NFIXL + 1) / REAL (NFIXL)
      NEW   = .TRUE.

      DO I = 2, NFIX - 1

         UDAT(1) = ZERO

         IF (I <= NFIXL) THEN ! In case body sections are oddly spaced
            VDAT(1) = DV2 * REAL (I - 1)
            VFIX(NJ,I) = VDAT(1)
         ELSE
            VDAT(1) = VFIX(NJ,I)
         END IF

         UDAT(2:4) = UFIX(JP:JP+2,I)
         VDAT(2:4) = VFIX(JP:JP+2,I)

         DU = UDAT(2) * RMJ

         DO J = NJ + 1, JP - 1
            UFIX(J,I) = REAL (J - NJ) * DU
         END DO

         CALL LCSFIT (4, UDAT, VDAT, NEW, TIGHT, NEVAL, UFIX(NJ+1,I),
     >                VFIX(NJ+1,I), T)
      END DO

!!!   write (62) vfix(nj:jlbody,1:nfix), ufix(nj:jlbody,1:nfix)

      JE = NJ + NFIXA
      IP = 1

      DO I = 2, NFIX - 1

         JP = NJ
         DO J = NJ, JLBODY

            CALL PLINCUB (JLBODY, NFIX, NJ, JLBODY, 1, NFIX,
     >                    XBOD, YBOD, ZBOD, UBOD, VBOD,
     >                    UFIX(J,I), VFIX(J,I), LINUV, TIGHT, CLOSED,
     >                    JP, IP, EPS, X(J), Y(J), Z(J), IER)
         END DO

         Z(JLBODY) = ZERO

         IF (I > NFIXL) THEN
            X(NJ) = XBOD(NJ,I)
            Y(NJ) = YBOD(NJ,I)
            Z(NJ) = ZBOD(NJ,I)
         END IF

C        Redistribute uniformly in the circumferential direction:

         CALL CHORDS3D (MJ, X(NJ), Y(NJ), Z(NJ), NORM, TOTAL, T(NJ))

         JEVAL = 1
         NEW = .TRUE.
         DU  = T(JE) * RMJ

         DO J = NJ + 1, JE - 1
            TEVAL = DU * REAL (J - NJ)

            CALL PLSCRV3D (MJ, X(NJ), Y(NJ), Z(NJ), T(NJ),
     >                     LOOSE, NEW, CLOSED, TEVAL, JEVAL,
     >                     XFIX(J,I), YFIX(J,I), ZFIX(J,I), DERIVS)
            NEW = .FALSE.
         END DO

         DO J = JE, JLBODY
            XFIX(J,I) = X(J)
            YFIX(J,I) = Y(J)
            ZFIX(J,I) = Z(J)
         END DO

      END DO ! Next I

C     Address (B(i,j)) = Address (B(1,1)) + (j - 1)idim + i - 1

      I = IXB(L) + NJ * IDIM + (IDIM - NFIXL)
      DO IP = 2, NFIXL + 1
         IB = I
         DO J = NJ + 1, JLBODY
            XBASE(IB) = XFIX(J,IP) ! (i,j) = (IDIM-NFIXL+1, NJ+1) on 1st pass
            YBASE(IB) = YFIX(J,IP)
            ZBASE(IB) = ZFIX(J,IP)
            IB = IB + IDIM
         END DO
         I = I + 1
      END DO

      IP = NFIXL
      JP = JLBODY - NJ + 1

      IC = IXB(L+2) + IDIM2
      DO I = 1, NFIXR
         IP = IP + 1
         I2 = IC
         DO J = NJ + 1, JLBODY
            XBASE(I2) = XFIX(J,IP) ! (i,j) = (1, 2) on 1st pass
            YBASE(I2) = YFIX(J,IP)
            ZBASE(I2) = ZFIX(J,IP)
            I2 = I2 + IDIM2
         END DO
         IC = IC + 1
      END DO

      END SUBROUTINE FIXROOT

C***********************************************************************
C
      SUBROUTINE FSPLIT (MXPATCH, K1, KPATCH, NWNGINT,
     >                   JMXCOM, MXSUBJ, NSUBJ, JSUB, IWRIT)
C
C     FSPLIT implicitly splits the initial fuselage patches along axial
C     lines through all "wing" root leading edges.  This is needed because
C     the design code's surface grid perturbation scheme must have patch
C     u/v = 0 or 1 at all critical intersection points so that grid points
C     near patch edges are perturbed continuously across the edges.
C     See also the NNSPLIT = 2 option.
C
C     02/11/98  DAS  Adaptation of James's original SPLITPAN.
C     02/18/99   "   Panels are packed now.
C     04/24/99   "   Subpatch version - doesn't need X/Y/ZBASE(*), etc.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   MXPATCH,             ! Most patches allowed for
     >   K1, KPATCH,          ! Patch number of nose & aft-body patches
     >   NWNGINT,             ! Number of fuselage/wing intersections
     >   JMXCOM(MXPATCH),     ! Active patch J dimensions
     >   MXSUBJ

      INTEGER, INTENT (OUT) ::
     >   NSUBJ(MXPATCH),      ! Subpatch info.
     >   JSUB(0:MXSUBJ,MXPATCH)

      INTEGER, INTENT (IN) ::
     >   IWRIT

C     Local variables:

      INTEGER
     >   I, ICUT, J, JDIM, K, NCUT, JCUT(NWNGINT+1), JFCUT(NWNGINT)

C     Execution:

      NCUT = NWNGINT ! Unless some J split indices are equal

      K = K1 + 1
      DO I = 1, NCUT ! Split indices = J dims. of under-wing interim patches
         JFCUT(I) = JMXCOM(K)
         K = K + 3
      END DO

C     Bubble sort of the split indices into ascending order:

      DO I = 2, NCUT
         DO J = NCUT, I, -1
            IF (JFCUT(J-1) > JFCUT(J)) THEN
               K          = JFCUT(J-1)
               JFCUT(J-1) = JFCUT(J)
               JFCUT(J)   = K
            END IF
         END DO
      END DO

C     Eliminate duplicates.

      JCUT(1) = JFCUT(1)
      J = 1
      DO I = 2, NCUT
         IF (JFCUT(I) /= JCUT(J)) THEN
            J = J + 1
            JCUT(J) = JFCUT(I)
         END IF
      END DO
      NCUT = J + 1 ! Normally NWNGINT + 1
      JCUT(NCUT) = JMXCOM(K1) ! Crown J of nose

C     Split the patches implicitly:

      DO K = K1, KPATCH ! For each interim patch

         IF (MOD (K - K1, 3) < 2) THEN ! Not above a wing

            JDIM = JMXCOM(K)
            DO I = 1, NCUT
               IF (JCUT(I) <= JDIM) THEN
                  JSUB(I,K) = JCUT(I)
                  NSUBJ(K)  = I
               END IF
            END DO

         ELSE ! Above a wing

            I = 0
            DO ICUT = 2, NCUT
               J = JCUT(ICUT) - JDIM
               IF (J > 0) THEN
                  I = I + 1
                  JSUB(I,K) = J + 1
               END IF
               NSUBJ(K) = I
            END DO

         END IF

      END DO

      END SUBROUTINE FSPLIT

C***********************************************************************
C
      SUBROUTINE FSURF (NWSURF, MXINTR, MXPATCH, KPATCH,
     >                  NFSTN, NNWING, JLBODY, NNROOT, NNSPLIT,
     >                  NBSPLIT, IBSPLIT, NNTAIL, IFUS, METHOD, EPS,
     >                  XF, YF, ZF, MXNWING, LOSTINTR,
     >                  NINTR, XINTR, YINTR, ZINTR, IINTR, JINTR,
     >                  ISTN1, ISTN2, NFIXL, NFIXR, NFIXA, NFIXB,
     >                  LENXB, XBASE, YBASE, ZBASE, LASTXB, MOSTXB, IXB,
     >                  MXSUBI, MXSUBJ, NSUBI, NSUBJ, ISUB, JSUB,
     >                  IMXCOM, JMXCOM, IWRIT, FAIL)
C
C     FSURF imposes regular patches on a regularized fuselage surface,
C     which may be intersected by one or more wing-type surfaces.
C
C     Pre-09/97  JJR  Original FSURF/FINSERT/FBARREL scheme of WBSURF.
C     Oct. 1997  DAS  Cleaned it up.
C     Feb. 1999   "   Packed-patch version.
C     04/26/99    "   Avoided redistributions of the NNSPLIT = 2 option
C                     by using the correct numbers of J points directly.
C     05/12/99    "   Body patch(es) split implicitly <-> pylon/diverter.
C     05/29/99    "   Started handling the presence of vertical fins.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   NWSURF,              ! # wing/nacelle surfaces (or 1 if 0)
     >   MXINTR,              ! Max. # pts. per intersection
     >   MXPATCH              ! Max. # patches expected

      INTEGER, INTENT (INOUT) ::
     >   KPATCH               ! In & out with the next patch # to fill

      INTEGER, INTENT (IN) ::
     >   NFSTN,               ! # body defining sections
     >   NNWING,              ! # wing surfaces (>= 0)
     >   JLBODY,              ! JLBODY - 1 = NJBODY = # J pts. input per section
     >   NNROOT(NWSURF),      ! NNROOT(IW) = -1 if surf. IW intersects the body
     >   NNSPLIT              ! 1 = original splitting at all intersection LEs;
                              ! 2 = just 1 split for off-wing patches

      INTEGER, INTENT (IN) ::
     >   NBSPLIT              ! NBSPLIT = 1: up/lo body splits <-> diverter LE;
                              !         = 2: lower body splits <-> pylon LE/TE;
      INTEGER, INTENT (IN) ::
     >   IBSPLIT(3)           ! IBSPLIT(1)   = wing surface <-> splitting, or 0
                              ! IBSPLIT(2:3) = implicit body split indices
      INTEGER, INTENT (IN) ::
     >   NNTAIL(*)            ! NNTAIL(IW) = +n means tail +n <->  dorsal fin IW
                              !    "   "   = -n   "    "   +n <-> ventral fin IW
                              ! where NNROOT(n) = -1 (body-mounted tail);
                              ! otherwise NNTAIL(IW) = 0
      INTEGER, INTENT (IN) ::
     >   IFUS                 ! Component # for diagnostics

      CHARACTER, INTENT (IN) ::
     >   METHOD * 1           ! 1-D spline method for fuselage barrels

      REAL, INTENT (IN) ::
     >   EPS,                 ! Tolerance for PLINCUB/RIPPLE2D
     >   XF(NFSTN),           ! Regularized body sections
     >   YF(JLBODY,NFSTN),
     >   ZF(JLBODY,NFSTN)

      INTEGER, INTENT (IN) ::
     >   MXNWING              ! Max. # wing surfaces expected

      LOGICAL, INTENT (IN) ::
     >   LOSTINTR(MXNWING)    ! T means intersection indices need recovering

      INTEGER, INTENT (IN) ::
     >   NINTR(NWSURF)        ! # pts. for each intersection

      REAL, INTENT (INOUT), DIMENSION (MXINTR,NWSURF) ::
     >   XINTR, YINTR, ZINTR  ! Intersection coordinates

      INTEGER, INTENT (INOUT), DIMENSION (MXINTR,NWSURF) ::
     >   IINTR, JINTR         ! Body cell indices, possibly recovered here

      INTEGER, INTENT (IN), DIMENSION (NWSURF) ::
     >   ISTN1, ISTN2,        ! Body section ranges for each wing surface
     >   NFIXL, NFIXR,        ! FIXROOT controls
     >   NFIXA, NFIXB

      INTEGER, INTENT (IN) ::
     >   LENXB                ! Packed patch space limit

      REAL, INTENT (INOUT), DIMENSION (LENXB) ::
     >   XBASE, YBASE, ZBASE  ! Patch coordinates

      INTEGER, INTENT (INOUT) ::
     >   LASTXB, MOSTXB,      ! Packed patch indices (last, biggest, first)
     >   IXB(MXPATCH)

      INTEGER, INTENT (IN) :: ! Subpatch info.
     >   MXSUBI, MXSUBJ

      INTEGER, INTENT (OUT) ::
     >   NSUBI(MXPATCH),
     >   NSUBJ(MXPATCH),
     >   ISUB(0:MXSUBI,MXPATCH),
     >   JSUB(0:MXSUBJ,MXPATCH)

      INTEGER, INTENT (INOUT) ::
     >   IMXCOM(MXPATCH), JMXCOM(MXPATCH) ! Patch dimensions

      INTEGER, INTENT (IN) ::
     >   IWRIT

      LOGICAL, INTENT (OUT) ::
     >   FAIL

C     Local constants:

      INTEGER, PARAMETER :: NMIN = 10 ! Min. # body sections between wings
      REAL,    PARAMETER :: HALF = 0.5, ONE = 1.
      LOGICAL, PARAMETER :: NORMALIZE = .TRUE.

C     Local variables:

      INTEGER
     >   I, I1, I2, IB1, IB2, IDIM, IER, IF, II, ILAST, ILE, ILEAD,
     >   ILST, IFST, ISUB1, IW, IW1, IX, IXB1, J, J1, J2, J3, JLEAD,
     >   K1, K2, K3, K4, L, L1, L2, L3, N, NFIN0, NFIX, NH, NI, NJ,
     >   NJA, NJB, NJBODY, NJT, NWNGINT

      INTEGER, DIMENSION (NWSURF) ::
     >   IFIN, ILEDGE, IINTTL, IINTTU, IWING, JINTTL, JINTTU

      REAL
     >   DX, R1, R2, RANGE, STOTAL, XLST, XFST, XKEY, XYZKEY(3),
     >   XYZITL(3,NWSURF), XYZITU(3,NWSURF), YKEY, ZKEY

      REAL, ALLOCATABLE, DIMENSION (:) ::
     >   SFINL, SFINL1, SFINL2, XTEMP, YTEMP, ZTEMP

      LOGICAL
     >   DORSAL

C     Execution:
C     ----------

      K1      = KPATCH     ! First (& possibly only) patch to fill
      NJBODY  = JLBODY - 1 ! Regularization allows for inserting a point
      NWNGINT = 0          ! # wing/body intersections
      NFIN0   = 0          ! # vertical fins NOT opposite a wing/tail/canard

      DO IW = 1, NNWING
         N = NNROOT(IW)
         IF (N == -1) THEN       ! Wing/body intersection
            NWNGINT = NWNGINT + 1
            IWING(NWNGINT) = IW
         ELSE IF (N == -2) THEN  ! Vertical fin/body intersection
            I = NNTAIL(IW)
            IF (I == 0) THEN     ! Fin is not opposite a wing/tail/canard
               NFIN0 = NFIN0 + 1
            END IF
         END IF
      END DO

C     Treat the no-wing case (complicated by possible fins):
C     ------------------------------------------------------

      IF (NWNGINT == 0) THEN ! No wing/body intersections

         JMXCOM(K1) = NJBODY
         JSUB(1,K1) = NJBODY

         IF (NFIN0 == 0) THEN ! No fin/body intersections either

C           Check for room for a plain body paneling:

            CALL XBCHECK (LENXB, MXPATCH, K1, LASTXB, NFSTN, NJBODY,
     >                    IFUS, IWRIT, MOSTXB)

            IMXCOM(K1) = NFSTN
            ISUB(1,K1) = NFSTN

            I1 = LASTXB + 1
            IXB(K1) = I1
            LASTXB  = LASTXB + NFSTN * NJBODY ! When done

            DO J = 1, NJBODY
               DO I = 1, NFSTN ! Fuselage patch (i,j) indices keep i streamwise
                  XBASE(I1) = XF(I)
                  YBASE(I1) = YF(J,I)
                  ZBASE(I1) = ZF(J,I)
                  I1 = I1 + 1
               END DO
            END DO

         ELSE ! At least one fin intersects the plain body

C           Set up the I subpatch boundaries:

            II = 1 ! IFIN(*) index
            I  = 1 ! I subpatch counter

            DO IW = 1, NNWING    ! Assume fins are ordered by X as for wings ...

               N = NNROOT(IW)    ! ... (though fins & wings may be interleaved)

               IF (N == -2) THEN ! Surface IW is a fin intersecting the body

                  IFIN(II) = IW
                  II = II + 1
                  NI = NINTR(IW)
                  NH = (NI + 1) / 2
                  I1 = IINTR(NH,IW) ! Body I at fin LE

                  IF (I == 1) THEN
                     ISUB(I,K1) = I1 + 1
                  ELSE
                     ISUB(I,K1) = ISUB(I-1,K1) + (I1 - I2 + 1)
                  END IF

                  ISUB(I+1,K1) = ISUB(I,K1) + NH - 1
                  I2 = IINTR(NI,IW) ! Body I at fin TE, for next fin or wrap-up
                  I  = I + 2

               END IF

            END DO ! Next surface to check for being a fin

            NSUBI(K1) = I
            NI = ISUB(I-1,K1) + (NFSTN - I2)
            ISUB(I,K1) = NI
            IMXCOM(K1) = NI
            IDIM = NI
            NJT  = NJBODY

C           Check for room for a plain body + fin(s):

            CALL XBCHECK (LENXB, MXPATCH, K1, LASTXB, NI, NJT,
     >                    IFUS, IWRIT, MOSTXB)

            I1 = LASTXB + 1
            LASTXB  = LASTXB + NI * NJT ! When done
            IXB1    = I1
            IXB(K1) = I1

C           Fill each subpatch in the I direction:

            ALLOCATE (SFINL(NJT), XTEMP(NJT), YTEMP(NJT), ZTEMP(NJT))

            R1 = ONE / REAL (NJT - 1)

            DO J = 1, NJT - 1 ! Normalized uniform distribution
               SFINL(J) = R1 * REAL (J - 1)
            END DO

            SFINL(NJT) = ONE

            II  = 1 ! IFIN(*) counter
            IB1 = 1 ! Body section counter with IB2

            DO L = 1, NSUBI(K1)

               ISUB1 = ISUB(L-1,K1) ! First patch I for subpatch L

               IF (MOD (L, 2) == 1) THEN ! Odd subpatch number - not at a fin

                  IF (L < NSUBI(K1)) THEN
                     IF  = IFIN(II)
                     NI  = NINTR(IF)
                     NH  = (NI + 1) / 2
                     IB2 = IINTR(NH,IF)
                  ELSE
                     IB2 = NFSTN
                  END IF

                  IF (L > 1) THEN ! Don't include fin TE station
                     I1 = IXB1 + ISUB1 ! First element of subpatch to fill
                  END IF

                  DO J = 1, NJT
                     IX = I1
                     I1 = I1 + IDIM
                     DO I = IB1, IB2
                        XBASE(IX) = XF(I)
                        YBASE(IX) = YF(J,I)
                        ZBASE(IX) = ZF(J,I)
                        IX = IX + 1
                     END DO
                  END DO

               ELSE ! Even subpatch number - at a fin

                  DORSAL = YINTR(NH,IF) > YF(NJT/2,IB2)

                  CALL FINSERTV () ! Internal procedure;

                  IB1 = IINTR(NI,IF) + 1 ! Body stn. for next plain subpatch
                  II  = II + 1 ! Next fin pointer

               END IF

            END DO ! Next subpatch

            DEALLOCATE (SFINL, XTEMP, YTEMP, ZTEMP)

         END IF ! Finished plain body + fin(s) case

         GO TO 990 ! Avoid indenting the bulk of the code; increment KPATCH

      END IF ! End of no-wing-intersection case.


C     There is at least one wing/body intersection.
C     ---------------------------------------------

C     The NNSPLIT = 2 case can avoid redistribution and awkward repacking
C     if the numbers of J points above and below each wing are determined
C     ahead of time rather than indirectly from the loop calling FINSERT.
C     We might as well establish the numbers of I points while we're at it.

      NJA = 999 ! See NNSPLIT = 2 case below
      NJB = 0

      DO II = 1, NWNGINT

         IW = IWING(II)

C        Find the leading edge of the intersection, which may not be the
C        middle index in spite of any rectifying after WSURF, because WDSURF
C        may have had to fiddle with the intersection.

         N  = NINTR(IW)
         I1 = N / 3
         XYZKEY(1) = XINTR(I1,IW) ! Nonmonotonic Xs possible for high sweeps

         DO I = I1, I1 + I1
            IF (XINTR(I,IW) < XYZKEY(1)) THEN
               XYZKEY(1) = XINTR(I,IW)
               ILE = I
            END IF
         END DO

         ILEDGE(II) = ILE
         XYZKEY(2)  = YINTR(ILE,IW)
         XYZKEY(3)  = ZINTR(ILE,IW)

C        Establish the body cell indices associated with this true LE:

         ILEAD = IINTR(ILE,IW) ! Unless we have to recover them
         JLEAD = JINTR(ILE,IW)

         IF (LOSTINTR(IW)) THEN

            I = (N + 1) / 2 ! Middle indices from WFINTR are best guess
            ILEAD = IINTR(I,IW)
            JLEAD = JINTR(I,IW)

            CALL BODYSRCH (NFSTN, JLBODY, ISTN1(IW), ISTN2(IW),
     >                     1, NJBODY, XF, YF, ZF, XYZKEY(1),
     >                     XYZKEY(2), XYZKEY(3), ILEAD, JLEAD, IER)
 
            IF (IER /= 0) THEN
               IF (IWRIT > 0) THEN
                  WRITE (IWRIT, '(/,A,I3,/,A,I4,1P,3X,A,3E15.6)')
     >               ' FSURF: LE trouble for wing surface', IW,
     >               ' ILE:', ILE, 'X/Y/ZLE:', XYZKEY
                  FAIL = .TRUE.
                  GO TO 999
               END IF
            END IF

            IINTR(ILE,IW) = ILEAD
            JINTR(ILE,IW) = JLEAD

         END IF

         NJB = MAX (NJB, JLEAD)
         NJA = MIN (NJA, JLEAD)


C        Treat the intersection trailing edge similarly so that patch sizes
C        can be established.  Temporarily close it to simplify FINSERT:

         XYZITL(1,II) = XINTR(1,IW) ! Save these for restoring below
         XYZITL(2,II) = YINTR(1,IW)
         XYZITL(3,II) = ZINTR(1,IW)
         XYZITU(1,II) = XINTR(N,IW)
         XYZITU(2,II) = YINTR(N,IW)
         XYZITU(3,II) = ZINTR(N,IW)
         DO I = 1, 3
            XYZKEY(I) = (XYZITL(I,II) + XYZITU(I,II)) * HALF
         END DO
         XINTR(1,IW)  = XYZKEY(1)
         YINTR(1,IW)  = XYZKEY(2)
         ZINTR(1,IW)  = XYZKEY(3)
         XINTR(N,IW)  = XYZKEY(1)
         YINTR(N,IW)  = XYZKEY(2)
         ZINTR(N,IW)  = XYZKEY(3)

         I = IINTR(1,IW)
         J = JINTR(1,IW)
         IINTTL(II) = I
         JINTTL(II) = J
         IINTTU(II) = IINTR(N,IW)
         JINTTU(II) = JINTR(N,IW)

C        Establish the body cell indices associated with the closed TE:

         CALL BODYSRCH (NFSTN, JLBODY, ISTN1(IW), ISTN2(IW),
     >                  1, NJBODY, XF, YF, ZF, XYZKEY(1),
     >                  XYZKEY(2), XYZKEY(3), I, J, IER)
 
         IF (IER /= 0) THEN
            IF (IWRIT > 0) THEN
               WRITE (IWRIT, '(/,A,I3,/,A,I4,1P,3X,A,3E15.6)')
     >            ' FSURF: TE trouble for wing surface', IW,
     >            ' N:', N, 'X/Y/ZTEclosed:', XYZKEY
               FAIL = .TRUE.
               GO TO 999
            END IF
         END IF

         IINTR(1,IW) = I ! Temporarily, for the LOSTINTR(IW) = F case
         IINTR(N,IW) = I
         JINTR(1,IW) = J
         JINTR(N,IW) = J

      END DO ! Next intersection


      IF (NNSPLIT == 2) THEN
         NJA = JLBODY - NJA  ! "Above" LE (inclusive) for all patches
         NJB = NJB + 1       ! "Below"; allow for insertion at the LE
         NJT = NJA + NJB - 1 ! For nose, aft-body, and between intersections
      ELSE
         NJT = JLBODY
      END IF


C     We are now in a position to set up all (sub)patch dimensions:
C     -------------------------------------------------------------

C     Nose patch:

      L = K1
      JMXCOM(L) = NJT
      IW1 = IWING(1) ! Surface for first wing/body intersection
      ILE = ILEDGE(1)
      I2  = IINTR(ILE,IW1) ! Body I corresp. to first intersection LE

      CALL FINDFIN (0, IW1, IF) ! Internal procedure locates a possible fin

      IF (IF == 0) THEN
         IMXCOM(L) = I2 + 1
         ISUB(1,L) = IMXCOM(L)
      ELSE
         NI = NINTR(IF)
         NH = (NI + 1) / 2
         ISUB(1,L) = IINTR(NH,IF) + 1
         ISUB(2,L) = ISUB(1,L) + NH - 1
         I1 = IINTR(NI,IF)
         ISUB(3,L) = ISUB(2,L) + (I2 - I1 + 1)
         IMXCOM(L) = ISUB(3,L)
         NSUBI(L)  = 3
      END IF

      IF (NNSPLIT == 2) THEN
         NSUBJ(L)  = 2 ! Left of intersection; NSUBI defaults to 1
         JSUB(1,L) = NJB
         JSUB(2,L) = NJT
      ELSE
         JSUB(1,L) = NJT
      END IF

C     Patches below/above/between/behind wing intersection(s):

      DO II = 1, NWNGINT
         IW  = IWING(II)
         ILE = ILEDGE(II)
         N   = NINTR(IW)
         L1  = L + 1
         L2  = L + 2
         IMXCOM(L1) = ILE         ! Below intersection
         ISUB(1,L1) = ILE
         IMXCOM(L2) = N - ILE + 1 ! Above
         ISUB(1,L2) = IMXCOM(L2)

         IF (IW == IBSPLIT(1)) THEN

C           Implicit I splits corresponding to wing pylon/diverter(s):

            IF (NBSPLIT == 1) THEN ! Under-wing diverter: LE only
               NSUBI(L1)  = 2
               ISUB(1,L1) = IBSPLIT(2)
               ISUB(2,L1) = ILE
               NSUBI(L2)  = 2
               ISUB(1,L2) = IBSPLIT(3)
               ISUB(2,L2) = IMXCOM(L2)
            ELSE ! NBSPLIT = 2     ! Under-wing pylon: LE & TE
               NSUBI(L1)  = 3
               ISUB(1,L1) = IBSPLIT(2)
               ISUB(2,L1) = IBSPLIT(3)
               ISUB(3,L1) = ILE
            END IF
         END IF

         IF (NNSPLIT == 2) THEN
            JMXCOM(L1) = NJB
            JMXCOM(L2) = NJA
         ELSE
            JMXCOM(L1) = JINTR(ILE,IW) + 1
            JMXCOM(L2) = NJT - JMXCOM(L1) + 1
         END IF

         JSUB(1,L1) = JMXCOM(L1)
         JSUB(1,L2) = JMXCOM(L2)

         L = L + 3        ! Right of intersection
         JMXCOM(L) = NJT
         I1 = IINTR(1,IW) ! Body cell I for closed TE

         IF (II < NWNGINT) THEN ! Between intersections II & II + 1

            IW1 = IWING(II+1)
            ILE = ILEDGE(II+1)
            I2  = IINTR(ILE,IW1)

            CALL FINDFIN (IW, IW1, IF) ! Locate a possible fin

            IF (IF == 0) THEN
               IMXCOM(L) = MAX (I2 - I1 + 2, NMIN) ! Includes a TE & LE pt.
               ISUB(1,L) = IMXCOM(L)
            ELSE
               NI = NINTR(IF)
               NH = (NI + 1) / 2
               ISUB(1,L) = IINTR(NH,IF) - I1 + 2
               ISUB(2,L) = ISUB(1,L) + NH - 1
               I1 = IINTR(NI,IF)
               ISUB(3,L) = ISUB(2,L) + (I2 - I1 + 1)
               IMXCOM(L) = ISUB(3,L)
               NSUBI(L)  = 3
            END IF

         ELSE ! Aft body

            CALL FINDFIN (IW, 99, IF) ! Locate a possible fin

            IF (IF == 0) THEN
               IMXCOM(L) = NFSTN - I1 + 1 ! Aft body
               ISUB(1,L) = IMXCOM(L)
            ELSE
               NI = NINTR(IF)
               NH = (NI + 1) / 2
               ISUB(1,L) = IINTR(NH,IF) - I1 + 2
               ISUB(2,L) = ISUB(1,L) + NH - 1
               ISUB(3,L) = ISUB(2,L) + NFSTN - IINTR(NI,IF)
               IMXCOM(L) = ISUB(3,L)
               NSUBI(L)  = 3
            END IF

         END IF

         IF (NNSPLIT == 2) THEN
            NSUBJ(L)  = 2
            JSUB(1,L) = NJB
            JSUB(2,L) = NJT
         ELSE
            JSUB(1,L) = NJT
         END IF

      END DO
         

C     Check that there is room for all fuselage patches:
C     --------------------------------------------------

      I1 = LASTXB + 1
      K4 = L ! Aft body patch
      N  = 0

      DO L = K1, K4
         NI = IMXCOM(L) * JMXCOM(L)
         N  = NI + N
         IXB(L) = I1
         I1 = I1 + NI
      END DO

      CALL XBCHECK (LENXB, MXPATCH, K4, LASTXB, N, 1,
     >              IFUS, IWRIT, MOSTXB)

      LASTXB = LASTXB + N ! Upon return from FSURF

      ALLOCATE (SFINL(NJT), SFINL1(NJT), SFINL2(NJT),
     >          XTEMP(NJT),  YTEMP(NJT),  ZTEMP(NJT))


C     Panels above and below the intersections are treated first,
C     because the barrel distributions at the leading and trailing
C     edges provide the relative distributions for neighboring patches.
C     -----------------------------------------------------------------

      L = K1 ! Patch forward of each intersection

      DO II = 1, NWNGINT

         IW = IWING(II)
         N  = NINTR(IW)

         CALL FINSERT (IW, MXPATCH, L, NFSTN, JLBODY,
     >                 NMIN, IFUS, METHOD, XF, YF, ZF, LOSTINTR(IW),
     >                 N, XINTR(1,IW), YINTR(1,IW), ZINTR(1,IW),
     >                 IINTR(1,IW), JINTR(1,IW),
     >                 ILEDGE(II), ISTN1(IW), ISTN2(IW),
     >                 LENXB, XBASE, YBASE, ZBASE, IXB,
     >                 XTEMP, YTEMP, ZTEMP, SFINL,
     >                 IMXCOM, JMXCOM, IWRIT, FAIL)

         IF (FAIL) THEN
            IF (IWRIT > 0) THEN
               WRITE (IWRIT, '(/,A,I3)')
     >            ' FSURF: FINSERT trouble for wing surface', IW
               GO TO 999
            END IF
         END IF

         L = L + 3 ! Forward of next intersection, or aft body patch

      END DO


C     Panel the forward fuselage:
C     ---------------------------

      NI = IMXCOM(K1)
      NJ = NJT
      IXB1 = IXB(K1)
      IX = IXB1

      DO J = 1, NJ ! Singular nose point
         XBASE(IX) = XF(1)
         YBASE(IX) = YF(1,1)
         ZBASE(IX) = ZF(1,1)
         IX = IX + NI
      END DO

C     Normalized arcs for the end of the nose patch (inserted by FINSERT):

      CALL CHORDSRF (NI, NJ, NI, 1, NJ, XBASE(IXB1), YBASE(IXB1),
     >               ZBASE(IXB1), NORMALIZE, STOTAL, SFINL)

C     Impose these on the interior sections of the nose patch:

      CALL FINDFIN (0, IWING(1), IF) ! Locate a possible fin

      IF (IF == 0) THEN ! No nose fin

         DO I = 2, NI - 1

            CALL FBARREL (I, NI, NJ, NJBODY, METHOD,
     >                    XF(I), YF(1,I), ZF(1,I), NJ, SFINL, XTEMP,
     >                    XBASE(IXB1), YBASE(IXB1), ZBASE(IXB1))
         END DO

      ELSE ! Handle one nose fin (only, so far)

         ISUB1 = ISUB(1,K1)

         DO I = 2, ISUB1 - 1 ! Forward of the fin

            CALL FBARREL (I, NI, NJ, NJBODY, METHOD,
     >                    XF(I), YF(1,I), ZF(1,I), NJ, SFINL, XTEMP,
     >                    XBASE(IXB1), YBASE(IXB1), ZBASE(IXB1))
         END DO

         IDIM = NI         ! FINSERTV needs IDIM, NI, NH
         NI = NINTR(IF)
         NH = (NI + 1) / 2
         II = IINTR(NI,IF)
         DORSAL = YINTR(NI,IF) > YF(NJT/2,II)

         CALL FINSERTV () ! At the fin

         DO I = ISUB(2,K1) + 1, IDIM - 1 ! Aft of the fin
            II = II + 1

            CALL FBARREL (I, IDIM, NJ, NJBODY, METHOD,
     >                    XF(II), YF(1,II), ZF(1,II), NJ, SFINL, XTEMP,
     >                    XBASE(IXB1), YBASE(IXB1), ZBASE(IXB1))
         END DO

      END IF ! End of the nose paneling

      K4 = K1 + 3 ! Panel after the first intersection
 

C     Handle possible mid-fuselage panels between intersections:
C     ----------------------------------------------------------

      DO L = 2, NWNGINT ! Panel the portion between intersections L-1 & L

         IW1   = IWING(L-1)
         IW    = IWING(L)
         IFST  = IINTR(1,IW1)
         XFST  = XINTR(1,IW1)
         ILE   = ILEDGE(L)
         ILST  = IINTR(ILE,IW)
         XLST  = XINTR(ILE,IW)
         RANGE = XLST - XFST

         IF (RANGE < 10.* EPS) THEN
            IF (IWRIT > 0) WRITE (IWRIT,'(/,1X,A,2I5)')
     >         'FSURF: Overlapping wing roots: ', IW1, IW
            FAIL = .TRUE.
            GO TO 999
         END IF

C        K4 points to the panel between intersections L-1, L.
C        The first & last barrels have been loaded by FINSERT.
C        Normalized arcs of the first barrel between intersections:

         I1 = IXB(K4)
         NI = IMXCOM(K4)
         NJ = JMXCOM(K4)

         CALL CHORDSRF (NI, NJ, 1, 1, NJ,
     >                  XBASE(I1), YBASE(I1), ZBASE(I1),
     >                  NORMALIZE, STOTAL, SFINL1)

C        Normalized arcs of the last barrel:

         CALL CHORDSRF (NI, NJ, NI, 1, NJ,
     >                  XBASE(I1), YBASE(I1), ZBASE(I1),
     >                  NORMALIZE, STOTAL, SFINL2)

         ILAST = ILST - IFST + 2 ! Includes 2 pts. at TE and LE

         IF (ILAST < NMIN) THEN ! Force at least NMIN sections

            DX = RANGE / REAL (NMIN - 1)
            RANGE = ONE / RANGE

            DO II = 2, NMIN - 1

               XKEY = XFST + DX * REAL (II - 1)

               DO I = IFST, ILST
                  IF (XKEY >= XF(I) .AND. XKEY <= XF(I+1)) EXIT
               END DO

               R1 = (XKEY - XF(I)) / (XF(I+1) - XF(I))
               R2 = (XKEY - XFST) * RANGE

               DO J = 1, NJ
                  YTEMP(J) = (ONE - R1)*YF(J,I)   + R1*YF(J,I+1)
                  ZTEMP(J) = (ONE - R1)*ZF(J,I)   + R1*ZF(J,I+1)
                  SFINL(J) = (ONE - R2)*SFINL1(J) + R2*SFINL2(J)
               END DO

               CALL FBARREL (II, NI, NJ, NJBODY, METHOD,
     >                       XKEY, YTEMP, ZTEMP, NJ, SFINL, XTEMP,
     >                       XBASE(I1), YBASE(I1), ZBASE(I1))
            END DO

         ELSE ! At least NMIN input sections - use their Xs

            RANGE = ONE / RANGE

            CALL FINDFIN (IW1, IW, IF) ! Locate a possible fin

            IF (IF == 0) THEN ! No fin

               I = IFST
               DO II = 2, ILAST - 1
                  I = I + 1
                  XKEY = XF(I)
                  R2 = (XKEY - XFST) * RANGE

                  DO J = 2, NJ
                     SFINL(J) = (ONE - R2)*SFINL1(J) + R2*SFINL2(J)
                  END DO

                  CALL FBARREL (II, NI, NJ, NJBODY, METHOD, XKEY,
     >                          YF(1,I), ZF(1,I), NJ, SFINL, XTEMP,
     >                          XBASE(I1), YBASE(I1), ZBASE(I1))
               END DO

            ELSE ! Handle one fin between intersections

               ISUB1 = ISUB(1,K4)
               II = IFST

               DO I = 2, ISUB1 - 1 ! Forward of the fin
                  II = II + 1
                  XKEY = XF(II)
                  R2 = (XKEY - XFST) * RANGE

                  DO J = 2, NJ
                     SFINL(J) = (ONE - R2)*SFINL1(J) + R2*SFINL2(J)
                  END DO

                  CALL FBARREL (I, NI, NJ, NJBODY, METHOD,
     >                          XKEY, YF(1,II), ZF(1,II),
     >                          NJ, SFINL, XTEMP, XBASE(I1),
     >                          YBASE(I1), ZBASE(I1))
               END DO

C              We have to kludge NI and NH to use FINSERTV at each fin point:

               IXB1 = I1 ! FINSERTV needs IXB1, ISUB1, IDIM, NI, NH
               IDIM = NI
               NI = NINTR(IF)
               NH = (NI + 1) / 2
               I1 = NH
               I2 = NI
               II = IINTR(NI,IF)
               DORSAL = YINTR(NI,IF) > YF(NJT/2,II)

               DO I = I1, I2 ! For each fin point

                  NH = I ! FINSERTV does just one pass through its loop
                  NI = I
                  R2 = (XINTR(I,IF) - XFST) * RANGE

                  DO J = 2, NJ
                     SFINL(J) = (ONE - R2)*SFINL1(J) + R2*SFINL2(J)
                  END DO

                  CALL FINSERTV ()

                  ISUB1 = ISUB1 + 1 ! More of the kludge

               END DO

               DO I = ISUB(2,K4) + 1, IDIM - 1 ! Aft of the fin
                  II = II + 1
                  XKEY = XF(II)
                  R2 = (XKEY - XFST) * RANGE

                  DO J = 2, NJ
                     SFINL(J) = (ONE - R2)*SFINL1(J) + R2*SFINL2(J)
                  END DO

                  CALL FBARREL (I, IDIM, NJ, NJBODY, METHOD,
     >                          XKEY, YF(1,II), ZF(1,II),
     >                          NJ, SFINL, XTEMP, XBASE(IXB1),
     >                          YBASE(IXB1), ZBASE(IXB1))
               END DO

            END IF

         END IF

         K4 = K4 + 3 ! Next mid-fuselage patch or aft body

      END DO ! Next wing/body intersection


C     Aft body patch:
C     ---------------

      KPATCH = K4
      NI = IMXCOM(K4)
      NJ = JMXCOM(K4)
      IXB1 = IXB(K4)

C     Normalized arcs of the first barrel of the aft body:

      CALL CHORDSRF (NI, NJ, 1, 1, NJ, XBASE(IXB1), YBASE(IXB1),
     >               ZBASE(IXB1), NORMALIZE, STOTAL, SFINL)

      CALL FINDFIN (IW, 99, IF) ! Locate a possible fin

      IFST = IINTR(1,IW) + 1

      IF (IF == 0) THEN

         II = 1
         DO I = IFST, NFSTN
            II = II + 1

            CALL FBARREL (II, NI, NJ, NJBODY, METHOD,
     >                    XF(I), YF(1,I), ZF(1,I), NJ, SFINL, XTEMP,
     >                    XBASE(IXB1), YBASE(IXB1), ZBASE(IXB1))
         END DO

      ELSE ! Handle one tail fin (only, so far)

         ISUB1 = ISUB(1,K4)
         II = IFST

         DO I = 2, ISUB1 - 1 ! Forward of the fin

            CALL FBARREL (I, NI, NJ, NJBODY, METHOD,
     >                    XF(II), YF(1,II), ZF(1,II), NJ, SFINL, XTEMP,
     >                    XBASE(IXB1), YBASE(IXB1), ZBASE(IXB1))
            II = II + 1
         END DO

         IDIM = NI    ! FINSERTV uses IDIM, NI, NH
         NI = NINTR(IF)
         NH = (NI + 1) / 2
         II = IINTR(NI,IF)
         DORSAL = YINTR(NI,IF) > YF(NJT/2,II)

         CALL FINSERTV () ! At the fin

         DO I = ISUB(2,K4) + 1, IDIM ! Aft of the fin
            II = II + 1

            CALL FBARREL (I, IDIM, NJ, NJBODY, METHOD,
     >                    XF(II), YF(1,II), ZF(1,II), NJ, SFINL, XTEMP,
     >                    XBASE(IXB1), YBASE(IXB1), ZBASE(IXB1))
         END DO

      END IF ! End of the aft-body paneling

      DEALLOCATE (SFINL, SFINL1, SFINL2, XTEMP, YTEMP, ZTEMP)


C     Redistribute subpatches near any dorsal fins associated with tails:
C     -------------------------------------------------------------------

      DO IF = 1, NNWING

         IF (NNROOT(IF) == -2) THEN ! Surface IF is a vertical fin

            IW = NNTAIL(IF)

            IF (IW > 0) THEN ! Dorsal fin IF is alongside tail surface IW

               DO II = 1, NWNGINT
                  IF (IWING(II) == IW) EXIT
               END DO

               L1 = K1 + 3 * II - 3 ! Patch left of the tail leading edge
               L2 = L1 + 2
               L3 = L1 + 3
               J1 = IXB(L1)
               J2 = IXB(L2)
               J3 = IXB(L3)
               NH = (NINTR(IF) + 1) / 2
            
               CALL FDORSAL (L1, L2, L3, IMXCOM, JMXCOM, METHOD,
     >                       XBASE(J1), YBASE(J1), ZBASE(J1),
     >                       XBASE(J2), YBASE(J2), ZBASE(J2),
     >                       XBASE(J3), YBASE(J3), ZBASE(J3), NH,
     >                       XINTR(NH,IF), YINTR(NH,IF), ZINTR(NH,IF),
     >                       IWRIT)

            ELSE IF (IW < 0) THEN ! Ventral fin IF is alongside tail surface -IW

               DO II = 1, NWNGINT
                  IF (IWING(II) == -IW) EXIT
               END DO

               L1 = K1 + 3 * II - 3 ! Patch left of the tail leading edge
               L2 = L1 + 1
               L3 = L1 + 3
               J1 = IXB(L1)
               J2 = IXB(L2)
               J3 = IXB(L3)
               NH = (NINTR(IF) + 1) / 2
            
               CALL FVENTRAL (L1, L2, L3, IMXCOM, JMXCOM, METHOD,
     >                        XBASE(J1), YBASE(J1), ZBASE(J1),
     >                        XBASE(J2), YBASE(J2), ZBASE(J2),
     >                        XBASE(J3), YBASE(J3), ZBASE(J3), NH,
     >                        XINTR(NH,IF), YINTR(NH,IF), ZINTR(NH,IF),
     >                        IWRIT)
            END IF

         END IF

      END DO


C     Fix skewed cells near rounded leading edge regions?
C     ---------------------------------------------------

      DO II = 1, NWNGINT

         IW = IWING(II)
         NFIX = NFIXL(IW) + NFIXR(IW) + 1

         IF (NFIX > 2) THEN

            L = K1 + 3 * II - 3 ! Patch left of the root leading edge

            CALL FIXROOT (MXPATCH, LENXB, IXB, IMXCOM, JMXCOM, L,
     >                    NJT, NFIXL(IW), NFIXR(IW), NFIX, NFIXA(IW), 
     >                    NFIXB(IW), EPS, XBASE, YBASE, ZBASE)
         END IF

C        Restore the trailing edge data:

         XINTR(1,IW) = XYZITL(1,II)
         YINTR(1,IW) = XYZITL(2,II)
         ZINTR(1,IW) = XYZITL(3,II)
         IINTR(1,IW) = IINTTL(II)
         JINTR(1,IW) = JINTTL(II)
         XINTR(N,IW) = XYZITU(1,II)
         YINTR(N,IW) = XYZITU(2,II)
         ZINTR(N,IW) = XYZITU(3,II)
         IINTR(N,IW) = IINTTU(II)
         JINTR(N,IW) = JINTTU(II)

      END DO


C     Define subpatches along axial lines through each wing leading edge?
C     -------------------------------------------------------------------

      IF (NNSPLIT == 1) THEN

         CALL FSPLIT (MXPATCH, K1, K4, NWNGINT, JMXCOM,
     >                MXSUBJ, NSUBJ, JSUB, IWRIT)

      END IF ! Else if NNSPLIT = 2 the subpatches have been defined above


  990 KPATCH = KPATCH + 1 ! Next available patch
      FAIL = .FALSE.

  999 RETURN

C     FSURF internal procedures:

      CONTAINS

!        -------------------------------------------------------------------
         SUBROUTINE FINDFIN (ILEFT, IRIGHT, IFIN) ! Find a fin between wings
!        -------------------------------------------------------------------

!        Arguments:

         INTEGER, INTENT (IN) ::
     >      ILEFT,              ! 0  = nose point,  else wing IW  left of fin;
     >      IRIGHT              ! 99 = tail end of body, else IW right of fin
         INTEGER, INTENT (OUT) ::
     >      IFIN                ! 0  means there no fin was found, else surface
                                ! IW > 0 is a fin between the indicated points
!        Local variables:

         INTEGER IW, N, NH, NI
         REAL    XLEFT, XRIGHT

!        Execution:

         IFIN = 0

         IF (ILEFT == 0) THEN
            XLEFT = XF(1)
         ELSE
            XLEFT = XINTR(1,ILEFT) ! Wing/body intersection trailing edge pt.
         END IF

         IF (IRIGHT == 99) THEN
            XRIGHT = XF(NFSTN)
         ELSE
            NI = NINTR(IRIGHT)
            XRIGHT = XINTR((NI+1)/2,IRIGHT) ! ~ Wing/body intersection LE pt.
         END IF

         DO IW = 1, NNWING
            N = NNROOT(IW)
            IF (N == -2) THEN ! It's a fin
               N = NNTAIL(IW)
               IF (N == 0) THEN ! Fin is not opposite a wing-type intersection
                  NI = NINTR(IW)
                  NH = (NI + 1) / 2
                  IF (XLEFT < XINTR(NH,IW) .AND. XINTR(NI,IW) < XRIGHT)
     >               IFIN = IW
               END IF
            END IF
         END DO

         END SUBROUTINE FINDFIN

!        -----------------------------------------------------------------
         SUBROUTINE FINSERTV () ! Generate a plain-body subpatch at a fin.
!        -----------------------------------------------------------------

!        Local variables:

         INTEGER I, ILEFT, IS, J, JLEFT
         REAL    R

!        Execution:

         IS = ISUB1 ! First patch I of this subpatch to fill

         IF (DORSAL) THEN

            DO I = NH, NI

               XKEY  = XINTR(I,IF)
               YKEY  = YINTR(I,IF)
               ZKEY  = ZINTR(I,IF)
               ILEFT = IINTR(I,IF)
               JLEFT = JINTR(I,IF)

!              Construct a body section at the intersection point:

               R = (XKEY - XF(ILEFT)) / (XF(ILEFT+1) - XF(ILEFT))

               DO J = 1, JLEFT
                  YTEMP(J) = (ONE - R) * YF(J,ILEFT) + R * YF(J,ILEFT+1)
                  ZTEMP(J) = (ONE - R) * ZF(J,ILEFT) + R * ZF(J,ILEFT+1)
               END DO

               J = JLEFT + 1
               YTEMP(J) = YKEY
               ZTEMP(J) = ZKEY

!              Respline the interpolated barrel:

               CALL FBARREL (IS, IDIM, NJT, J, METHOD,
     >                       XKEY, YTEMP, ZTEMP, NJT, SFINL, XTEMP,
     >                       XBASE(IXB1), YBASE(IXB1), ZBASE(IXB1))
               IS = IS + 1 

            END DO

         ELSE ! Ventral fin

            DO I = NH, NI

               XKEY  = XINTR(I,IF)
               YKEY  = YINTR(I,IF)
               ZKEY  = ZINTR(I,IF)
               ILEFT = IINTR(I,IF)
               JLEFT = JINTR(I,IF)

!              Construct a body section at the intersection point:

               R = (XKEY - XF(ILEFT)) / (XF(ILEFT+1) - XF(ILEFT))

               YTEMP(JLEFT) = YKEY
               ZTEMP(JLEFT) = ZKEY

               DO J = JLEFT + 1, NJBODY
                  YTEMP(J) = (ONE - R) * YF(J,ILEFT) + R * YF(J,ILEFT+1)
                  ZTEMP(J) = (ONE - R) * ZF(J,ILEFT) + R * ZF(J,ILEFT+1)
               END DO

               J = JLEFT

!              Respline the interpolated barrel:

               CALL FBARREL (IS, IDIM, NJT, NJBODY - J + 1, METHOD,
     >                       XKEY, YTEMP(J), ZTEMP(J), NJT, SFINL,
     >                       XTEMP, XBASE(IXB1),YBASE(IXB1),ZBASE(IXB1))
               IS = IS + 1 

            END DO

         END IF

         END SUBROUTINE FINSERTV

      END SUBROUTINE FSURF

C***********************************************************************
C
      SUBROUTINE FVENTRAL (L1, L2, L3, IDIM, JDIM, METHOD, X1, Y1, Z1,
     >                     X2, Y2, Z2, X3, Y3, Z3, NFIN,
     >                     XFIN, YFIN, ZFIN, IWRIT)
C
C     FVENTRAL redistributes the three standard fuselage (sub)patches
C     associated with a horizontal tail in order to accommodate an
C     accompanying ventral fin (represented by its fin/body intersection).
C     The point counts of the fin and tail are required to match.
C     Since the point counts do not change, the input patches are
C     overwritten with the desired results.
C
C     06/15/99  David Saunders  Analogue of FDORSAL.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   L1, L2, L3,           ! Patches forward, below, and aft of the tail
     >   IDIM(L3), JDIM(L3)    ! Patch dimensions

      CHARACTER, INTENT (IN) ::
     >   METHOD * 1            ! 1-D method for fuselage barrels also controls
                               ! use of PLBICUT/PLBICUBE or BILINCUT/PBILINT

      REAL, INTENT (INOUT), DIMENSION (IDIM(L1), JDIM(L1)) ::
     >   X1, Y1, Z1            ! Patch left of tail

      REAL, INTENT (INOUT), DIMENSION (IDIM(L2), JDIM(L2)) ::
     >   X2, Y2, Z2            ! Patch below tail

      REAL, INTENT (INOUT), DIMENSION (IDIM(L3), JDIM(L3)) ::
     >   X3, Y3, Z3            ! Patch aft of tail

      INTEGER, INTENT (IN) ::
     >   NFIN                  ! NFIN = IDIM(L2) is assumed

      REAL, INTENT (IN), DIMENSION (NFIN) ::
     >   XFIN, YFIN, ZFIN      ! Fin/body intersection pts.

      INTEGER, INTENT (IN) ::
     >   IWRIT                 ! For diagnostics

C     Local constants:

      REAL,      PARAMETER :: ONE = 1., ZERO = 0.
      LOGICAL,   PARAMETER :: NEW = .TRUE.
      CHARACTER, PARAMETER :: LOOSE = 'B'

C     Local variables:

      INTEGER
     >   I, I0, I2, IER, IP, ISIZE, IS, ISMAT, J, JS, JSIZE, JSMAT,
     >   LI, LUNERR, MI, NI
      REAL
     >   DRANGE, DU, EPS, P, Q, R, SMATRIX(4,4,3), UCUT, VCUT,
     >   UVRANGE(4), XINT, YINT, ZINT, XYZCUT(3)
      REAL, ALLOCATABLE, DIMENSION (:) ::
     >   EDGES
      REAL, ALLOCATABLE, DIMENSION (:,:) ::
     >   X, Y, Z, U, V, UNEW, VNEW
      LOGICAL
     >   BICUBIC

C     Storage:

      DATA UVRANGE /ZERO, ONE, ZERO, ONE/ ! See PLBICUT

C     Execution:

      EPS = MAX (5.* EPSILON (EPS), 1.E-7) ! For PLBICUT/RIPPLE2D
      BICUBIC = METHOD == 'B'

C     We need to gather the three (sub)patches as one in order to redistribute.
C     Set up the workspace:

      NI    = IDIM(L1)
      MI    = NI + NFIN - 1
      LI    = IDIM(L3)
      ISIZE = MI + LI - 1
      JSIZE = JDIM(L2)

      ALLOCATE (X(ISIZE,JSIZE), Y(ISIZE,JSIZE), Z(ISIZE,JSIZE))

      DO J = 1, JSIZE
         DO I = 1, NI
            X(I,J) = X1(I,J)
            Y(I,J) = Y1(I,J)
            Z(I,J) = Z1(I,J)
         END DO
      END DO

      DO J = 1, JSIZE
         IP = NI
         DO I = 2, NFIN
            IP = IP + 1
            X(IP,J) = X2(I,J)
            Y(IP,J) = Y2(I,J)
            Z(IP,J) = Z2(I,J)
         END DO
      END DO

      DO J = 1, JSIZE
         IP = MI
         DO I = 2, LI
            IP = IP + 1
            X(IP,J) = X3(I,J)
            Y(IP,J) = Y3(I,J)
            Z(IP,J) = Z3(I,J)
         END DO
      END DO

      ALLOCATE (U(ISIZE,JSIZE), UNEW(ISIZE,JSIZE),
     >          V(ISIZE,JSIZE), VNEW(ISIZE,JSIZE))

      U(1,1) = ZERO ! -999. suppresses normalization

      CALL PARAM2D (ISIZE, JSIZE, 1, ISIZE, 1, JSIZE, X, Y, Z, U, V)

C     Three of the edges in (u,v) space stay the same:

      DO I = 1, ISIZE
         UNEW(I,JSIZE) = U(I,JSIZE)
         VNEW(I,JSIZE) = ONE
      END DO

      DO J = 1, JSIZE - 1
         UNEW(1,J) = ZERO
         VNEW(1,J) = V(1,J)
         UNEW(ISIZE,J) = ONE
         VNEW(ISIZE,J) = V(ISIZE,J)
      END DO

C     Locate the fin/body intersection points in (u,v) space:

      DRANGE = X(ISIZE,1) - X(1,1)
      IP = NI
      IS = NI
      JS = 1

      CALL INTERVAL (ISIZE, X(1,1), XFIN(1), ONE, IS) ! Fin LE pt.

      R = (XFIN(1) - X(IS,1)) / (X(IS+1,1) - X(IS,1))
      UCUT = (ONE - R) * U(IS,1) + R * U(IS+1,1)
      UNEW(IP,1) = UCUT
      VNEW(IP,1) = ZERO
      VCUT = ZERO
      J = JSIZE / 2 ! Should be plenty
      ISMAT = 0
      JSMAT = 0
      LUNERR = -ABS (IWRIT)

      DO I = 1, NFIN ! Highly nonuniform spacing affects PLBICUBE, so redo I=1

         XYZCUT(1) = XFIN(I)
         XYZCUT(3) = ZFIN(I)

         IF (BICUBIC) THEN

            CALL PLBICUT (3, ISIZE, JSIZE, 1, ISIZE, 1, J,
     >                    X, Y, Z, U, V, IS, JS, ISMAT, JSMAT,
     >                    SMATRIX, EPS, DRANGE, XYZCUT, UVRANGE,
     >                    UCUT, VCUT, P, Q, LUNERR, IER)
         ELSE

            CALL BILINCUT (3, ISIZE, JSIZE, 1, ISIZE, 1, J,
     >                     X, Y, Z, U, V, IS, JS,
     >                     EPS, DRANGE, XYZCUT, UVRANGE,
     >                     UCUT, VCUT, P, Q, LUNERR, IER)
         END IF

         IF (IER > 3) THEN
            IF (IWRIT > 0) THEN
               WRITE (IWRIT, '(/, A, I2, /, A, I4, /, A, 1P, 3E15.6)')
     >            ' FVENTRAL: PLBICUT/BILINCUT error =', IER,
     >            ' Fin/body intersection pt. =', I,
     >            ' Intersection pt. (x,y,z)  =', XYZCUT
            END IF
            CALL SYNCH
            STOP
         END IF

         UNEW(IP,1) = UCUT
         VNEW(IP,1) = VCUT
         IP = IP + 1

      END DO

C     Match the relative water-line distribution fore & aft of the fin:

      R = UNEW(NI,1) / UNEW(NI,JSIZE)

      DO I = 2, NI - 1
         UNEW(I,1) = R * UNEW(I,JSIZE)
         VNEW(I,1) = ZERO
      END DO

      R = (ONE - UNEW(MI,1)) / (ONE - UNEW(MI,JSIZE))

      DO I = 2, LI - 1
         UNEW(IP,1) = UNEW(MI,1) + R * (UNEW(IP,JSIZE) - UNEW(MI,JSIZE))
         VNEW(IP,1) = ZERO
         IP = IP + 1
      END DO

C     New interior (u,v)s:

      ALLOCATE (EDGES(2*(ISIZE+JSIZE)))

      CALL TFI2D (ISIZE, 1, ISIZE, 1, JSIZE, UNEW, VNEW, EDGES)

      DEALLOCATE (EDGES)

C     Redistribute the interior points:

      I2 = 1
      I0 = 0

      DO I = 2, ISIZE - 1

         IF (I == NI) THEN
            I0 = NI - 1
         ELSE IF (I == MI) THEN
            I0 = MI - 1
         END IF

         IP = I - I0
         IS = I2
         JS = 1

         DO J = 2, JSIZE - 1

            IF (BICUBIC) THEN

               CALL PLBICUBE (0, ISIZE, JSIZE, 1, ISIZE, 1, JSIZE,
     >                        X, Y, Z, U, V, UNEW(I,J), VNEW(I,J),
     >                        IS, JS, ISMAT, JSMAT, SMATRIX, EPS, P, Q,
     >                        XINT, YINT, ZINT, XYZCUT, XYZCUT, IER)
            ELSE

               CALL PBILINT (0, ISIZE, JSIZE, 1, ISIZE, 1, JSIZE,
     >                       X, Y, Z, U, V, UNEW(I,J), VNEW(I,J),
     >                       IS, JS, EPS, P, Q, XINT, YINT, ZINT,
     >                       XYZCUT, XYZCUT, IER)
            END IF

            IF (IER /= 0) THEN
               IF (IWRIT > 0)
     >            WRITE (IWRIT, '(/,A,A,I2,A,4I5,/,A,2F13.5,A,3F13.5)')
     >               ' FVENTRAL: Marginal 2D interpolation. ',
     >               ' IER:', IER, '; I, J, IS, JS:', I, J, IS, JS,
     >               ' U, V:', UNEW(I,J), VNEW(I,J),
     >               '    X, Y, Z used:', XINT, YINT, ZINT
            END IF

            IF (J == 2) I2 = IS ! Crown pts. can be offset from water-line pts.

            IF (I < NI) THEN
               X1(IP,J) = XINT
               Y1(IP,J) = YINT
               Z1(IP,J) = ZINT
            ELSE IF (I < MI) THEN
               X2(IP,J) = XINT
               Y2(IP,J) = YINT
               Z2(IP,J) = ZINT
            ELSE
               X3(IP,J) = XINT
               Y3(IP,J) = YINT
               Z3(IP,J) = ZINT
            END IF

         END DO

      END DO

C     New keel points:

      CALL LCSFIT (ISIZE, U(1,1), X(1,1), NEW, LOOSE,
     >             NI - 2, UNEW(2,1), X1(2,1), Z(1,1))

      CALL LCSFIT (ISIZE, U(1,1), Y(1,1), NEW, LOOSE,
     >             NI - 2, UNEW(2,1), Y1(2,1), Z(1,1))

      X1(NI,1) = XFIN(1)
      Y1(NI,1) = YFIN(1)
      Z1(NI,1) = ZFIN(1)

      DO I = 1, NFIN
         X2(I,1) = XFIN(I)
         Y2(I,1) = YFIN(I)
         Z2(I,1) = ZFIN(I)
      END DO

      X3(1,1) = XFIN(NFIN)
      Y3(1,1) = YFIN(NFIN)
      Z3(1,1) = ZFIN(NFIN)

      CALL LCSFIT (ISIZE, U(1,1), X(1,1), NEW, LOOSE,
     >             LI - 2, UNEW(MI+1,1), X3(2,1), Z(1,1))

      CALL LCSFIT (ISIZE, U(1,1), Y(1,1), NEW, LOOSE,
     >             LI - 2, UNEW(MI+1,1), Y3(2,1), Z(1,1))

      DEALLOCATE (X, Y, Z, U, V, UNEW, VNEW)

C     Duplicate two edges:

      IP = NFIN

      DO J = 2, JSIZE - 1
         X1(NI,J) = X2(1,J)
         Y1(NI,J) = Y2(1,J)
         Z1(NI,J) = Z2(1,J)
         X2(IP,J) = X3(1,J)
         Y2(IP,J) = Y3(1,J)
         Z2(IP,J) = Z3(1,J)
      END DO

      END SUBROUTINE FVENTRAL

C*******************************************************************************
C
      SUBROUTINE MATCHDIV (IWING, ID, ID2, NWSURF, MXINTR,
     >                     IDIMW, KDIMW, XWING, YWING, ZWING,
     >                     IDIMD, KDIMD, XDIV2, YDIV2, ZDIV2,
     >                     XTRAP, EPS, IDEGREE,
     >                     ILWING,NHALF, ITL, ITU,
     >                     KSTN1, KSTN2, KWING1, KWING2,
     >                     NINTR, XINTR, YINTR, ZINTR,
     >                     TINTR, UINTR, VINTR, IINTR, JINTR, KINTR,
     >                     IWRIT, FAIL)

C     MATCHDIV compares, for adjacent diverters, the inboard edge of the
C     wing/diverter intersection for diverter ID2 with the outboard edge for
C     diverter ID, and adjusts the sections of diverter ID2 if necessary to
C     give a point match for the on-wing portions of the intersections as
C     needed by WDSURF paneling around the diverters.  The redistribution
C     does not change the number of points, NHALF(ID2).  Recomputing
C     the intersection appears unavoidable for determining the associated
C     quantities (t, u, v, i, j, k).
C
C     01/12/98  D. Saunders  Initial implementation.  Mark Rimlinger
C                            suggested use of PLBICUT on the diverter
C                            surface to determine the trailing edge t values.
C     01/31/99     "   "     Moved the loop over diverters up a level, now
C                            that X/Y/ZRG(*) are packed.
C
C*******************************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   IWING,               ! Index of wing associated with these diverters
     >   ID,                  ! Inboard diverter surface index
     >   ID2,                 ! Adjacent outboard diverter surface index
     >   NWSURF,              ! # "wing" + "nacelle" type surfaces
     >   MXINTR,              ! Intersection array dimension
     >   IDIMW,               ! Wing surface array dimensions
     >   KDIMW,
     >   IDIMD,               ! Diverter ID2 array dimensions
     >   KDIMD

      REAL, INTENT (IN), DIMENSION (IDIMW,KDIMW) ::
     >   XWING, YWING, ZWING  ! Regularized wing sections

      REAL, INTENT (INOUT), DIMENSION (IDIMD,KDIMD) ::
     >   XDIV2, YDIV2, ZDIV2  ! Regularized diverter ID2 sections

      REAL, INTENT (IN) ::
     >   EPS,                 ! For RIPPLE2D
     >   XTRAP(NWSURF)        ! X distances to extrapolate wing trailing edge

      INTEGER, DIMENSION (NWSURF) ::
     >   IDEGREE,             ! 3 or 1 = INTSEC4 or -4 for the wing
     >   ILWING,              ! # (wrap-around) chordwise gridded pts. per wing
     >   NHALF,               ! ILWING(*) / 2 + 1
     >   ITL,                 ! First/last diverter intersection pts. on wing
     >   ITU,       
     >   KSTN1,               ! Wing section ranges for diverter intersections
     >   KSTN2,       
     >   KWING1,              ! Diverter K ranges for intersections
     >   KWING2,       
     >   NINTR                ! Intersection point counts

      REAL, DIMENSION (MXINTR,NWSURF) ::
     >   XINTR, YINTR, ZINTR, TINTR, UINTR, VINTR

      INTEGER, DIMENSION (MXINTR,NWSURF) ::
     >   IINTR, JINTR, KINTR

      INTEGER, INTENT (IN) ::
     >   IWRIT

      LOGICAL, INTENT (OUT) ::
     >   FAIL

C     Local constants:

      INTEGER, PARAMETER :: KASE = 1
      LOGICAL, PARAMETER :: NEW  = .TRUE.
      REAL,    PARAMETER :: HALF = 0.5, ONE = 1., ZERO = 0.

C     Local variables:

      INTEGER
     >   I, I1, I2, IER, IS, ISMAT, J, K, K1, K2, KS, KSMAT,
     >   N, NHALF2, NLEFT
      REAL
     >   DRANGE, P, Q, SMATRIX(4,4,3), UCUT, UVRANGE(4), VCUT, XYZCUT(3)
      REAL, DIMENSION (IDIMD,KDIMD) ::
     >   UDIV, VDIV
      REAL, DIMENSION (NHALF(ID2)) ::
     >   XT, YT, ZT
      DATA
     >   UVRANGE /ZERO, ONE, ZERO, ONE/

C     Execution:

      FAIL   = .FALSE.
      N      = ITU(ID) - NHALF(ID) + 1 ! # points on inboard diverter upper edge
      NHALF2 = NHALF(ID2)

      IF (NHALF2 - ITL(ID2) + 1 /= N) THEN ! Adjust outboard diverter lower edge

         XYZCUT(1) = XINTR(ITL(ID2),ID2)
         XYZCUT(2) = YINTR(ITL(ID2),ID2)
         I1 = 1
         I2 = NHALF2
         K1 = KWING1(ID2)
         K2 = KWING2(ID2)
         UDIV(I1,K1) = ZERO ! -999. would not normalize diverter (u,v)s

         CALL PARAM2D (IDIMD, KDIMD, I1, I2, K1, K2,
     >                 XDIV2, YDIV2, ZDIV2, UDIV, VDIV)

C        Determine the diverter (u,v) at the wing TE as a biproduct of
C        calculating Z on the diverter surface corresp. to the (X,Y) of
C        the intersection trailing edge:

         IS     = ITL(ID2)
         KS     = K1
         ISMAT  = 0
         KSMAT  = 0
         DRANGE = XDIV2(I1,K1) - XDIV2(I2,K1) ! ~ Match u,v precision
         UCUT   = UDIV(IS,KS)
         VCUT   = VDIV(IS,KS)

         CALL PLBICUT (KASE, IDIMD, KDIMD, I1, I2, K1, K2,
     >                 XDIV2, YDIV2, ZDIV2, UDIV, VDIV,
     >                 IS, KS, ISMAT, KSMAT, SMATRIX, EPS,
     >                 DRANGE, XYZCUT, UVRANGE, UCUT, VCUT,
     >                 P, Q, -ABS (IWRIT), IER)

         IF (IER > 3) THEN
            IF (IWRIT > 0) THEN
               WRITE (IWRIT, '(/, A, I4, A, I3)')
     >            ' MATCHDIV: IER from PLBICUT =', IER,
     >            ' Diverter surface # =', ID2
            END IF
            FAIL = .TRUE.
            GO TO 999
         END IF

C        Impose N points on the "lower" surfaces of diverter ID2 sections
C        between t = UCUT and t = 1.

         NLEFT = NHALF2 - N + 1
         ITL(ID2) = NLEFT

         IF (NLEFT <= 1) THEN
            IF (IWRIT > 0) THEN
               WRITE (IWRIT, '(/, A, I3, A, 2I4)')
     >            ' MATCHDIV: Too few pts. on diverter surface', ID2,
     >            '   NHALF2 & N =', NHALF2, N
            END IF
            FAIL = .TRUE.
            GO TO 999
         END IF

         DO K = K1, K2

            DO I = 1, NHALF2
               XT(I) = XDIV2(I,K)
               YT(I) = YDIV2(I,K)
               ZT(I) = ZDIV2(I,K)
            END DO

C           Insert the point corresponding to UCUT at the wing trailing edge,
C           and match the number of points on the wing.  NLEFT is such that
C           the total remains NHALF2 (but MATCHN is more general).
C           UDIV(1:NHALF2,K) serves for all K because X/Y/ZRG are regular.

            CALL MATCHN (NHALF2, NLEFT, N, UCUT, XT, YT, ZT,
     >                   UDIV(1,K), XDIV2(1,K), YDIV2(1,K), ZDIV2(1,K))

         END DO

C        Recompute the inboard half of the wing/diverter J intersection.
C        NCALL = -4 tells WDINTR to do just the first half.
C        The starting guesses should still be close enough.

         CALL WDINTR (-4, IDIMD, KDIMD, MXINTR, NHALF2, ILWING(ID2),
     >                XDIV2, YDIV2, ZDIV2, IDIMW, KDIMW,
     >                1, NHALF(IWING), KSTN1(ID2), KSTN2(ID2),
     >                KWING1(ID2), KWING2(ID2), XWING, YWING, ZWING,
     >                XTRAP(ID2), IDEGREE(IWING),
     >                ITL(ID2), ITU(ID2), NINTR(ID2),
     >                XINTR(1,ID2), YINTR(1,ID2), ZINTR(1,ID2),
     >                UINTR(1,ID2), VINTR(1,ID2), TINTR(1,ID2),
     >                IINTR(1,ID2), JINTR(1,ID2), KINTR(1,ID2),
     >                IWRIT, FAIL)

         IF (FAIL) THEN
            IF (IWRIT > 0)
     >         WRITE (IWRIT,'(/,A,A,I3)') ' MATCHDIV: ',
     >            'Wing (re)intersection failed for diverter', ID2
            GO TO 999
         END IF

      END IF

  999 RETURN

      END SUBROUTINE MATCHDIV

C*******************************************************************************
C
      SUBROUTINE MATCHN (NDATA, NLEFT, NRIGHT, TINSERT, X, Y, Z, T,
     >                   XOUT, YOUT, ZOUT)
C
C        MATCHN is a somewhat specialized utility for inserting a point in
C     a 3-space curve and adjusting the number of points either side of it.
C     Preserving the character of the data point distribution appears to be
C     intractable, so the output points are uniform on either side of the
C     inserted point.
C
C     01/13/98  DAS  Initial implementation to deal with under-wing diverters.
C
C     Author:  David Saunders, Sterling Software/NASA Ames, Moffett Field, CA
C
C*******************************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER    NDATA      ! I   Input number of data points
      INTEGER    NLEFT,     ! I   Desired numbers of points either side
     >           NRIGHT     !     of (and including) the inserted point
      REAL       TINSERT    ! I   Parametric location of   "      "
      REAL       X(NDATA),  ! I   Data points
     >           Y(NDATA),
     >           Z(NDATA)
      REAL       T(NDATA)   ! I   Data parameterization
      REAL       XOUT(*),   !   O Redistributed coordinates:
     >           YOUT(*),   !     NLEFT + NRIGHT - 1 of them, with the
     >           ZOUT(*)    !     inserted point at I = NLEFT

      EXTERNAL   INTERVAL,  !     1-D search utility
     >           PLSCRV3D   !     3-space local cubic spline utility

C     Local constants:

      REAL,          PARAMETER :: DERIVS = -999., ONE = 1.
      LOGICAL,       PARAMETER :: CLOSED = .FALSE.
      CHARACTER * 1, PARAMETER :: METHOD = 'B' ! "Loose" fit

C     Local variables:

      INTEGER    I, IEVAL, NEVAL
      REAL       DT, R, TEVAL(NLEFT + NRIGHT)
      LOGICAL    NEW

C     Execution:

      DT = (TINSERT - T(1)) / REAL (NLEFT - 1)
      DO I = 1, NLEFT - 1
         TEVAL(I) = T(1) + REAL (I-1) * DT
      END DO
      TEVAL(NLEFT) = TINSERT

      DT = (T(NDATA) - TINSERT) / REAL (NRIGHT - 1)
      NEVAL = NLEFT + NRIGHT - 1

      DO I = NLEFT + 1, NEVAL - 1
         TEVAL(I) = TINSERT + REAL (I-NLEFT) * DT
      END DO

      IEVAL = 1
      NEW = .TRUE.
      DO I = 1, NEVAL - 1

         CALL PLSCRV3D (NDATA, X, Y, Z, T, METHOD, NEW, CLOSED,
     >                  TEVAL(I), IEVAL,
     >                  XOUT(I), YOUT(I), ZOUT(I), DERIVS)
         NEW = .FALSE.
      END DO

      XOUT(NEVAL) = X(NDATA)
      YOUT(NEVAL) = Y(NDATA)
      ZOUT(NEVAL) = Z(NDATA)

      END SUBROUTINE MATCHN

C*******************************************************************************
C
      SUBROUTINE NACSURF (IN, ID, KPATCH, MXPATCH, IDIM, KDIM, ILNAC,
     >                    KLNAC, HSCT, XNAC, YNAC, ZNAC, NHALF, ITL,ITU,
     >                    NINTR, XINTR, YINTR, ZINTR, IINTR, JINTR,
     >                    LENXB, XBASE, YBASE, ZBASE, LASTXB, MOSTXB,
     >                    IXB, MXSUBI, MXSUBJ, NSUBJ, ISUB, JSUB,
     >                    IMXCOM, JMXCOM, IWRIT)
C
C     NACSURF panels a nacelle surface following any pylon or diverter
C     intersection calculation.  The regularized, padded nacelle surface
C     is input as one surface and output as 2, 6, or 8 patches including
C     splits needed at the intersection leading & trailing edges and at
C     the nacelle leading edge.
C
C     Sept. 97  JJR  Original nacelle paneling as part of NSURF.
C     10-11/97  DAS  Streamlined, but still linear or 1-D parametric methods.
C     March 98   "   HSCT = T forces special treatment of the aft portion
C                    of an HSCT-type nacelle at the higher level.
C     09/14/98   "   Nacelles must be split at the leading edge.
C     02/06/99   "   Packed-regular-geometry version.
C     02/25/99   "   Packed-output-patch version: merged NINSERT into NACSURF.
C     04/29/99   "   Subpatch version.
C
C*******************************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   IN, ID                ! Nacelle & pylon/diverter component #s;
                               ! ID = 0 means plain nacelle - no intersection
      INTEGER, INTENT (INOUT) ::
     >   KPATCH                ! Next patch to use (in & out)

      INTEGER, INTENT (IN) ::
     >   MXPATCH,              ! Patch count limit
     >   IDIM, KDIM,           ! Plain nacelle array dimensions
     >   ILNAC, KLNAC          ! Active plain nacelle dimensions

      LOGICAL, INTENT (IN) ::
     >   HSCT                  ! TRUE means aft nacelle is left to NACSURF2

      REAL, INTENT (IN), DIMENSION (IDIM,KDIM) ::
     >   XNAC, YNAC, ZNAC      ! Full interior + exterior in (1:ILNAC,1:KLNAC)

      INTEGER, INTENT (IN) ::
     >   NHALF,                ! Nacelle leading edge index
     >   ITL, ITU,             ! Pylon/div. trailing edge indices on nacelle
     >   NINTR                 ! # nacelle/pylon intersection points

      REAL, INTENT (IN), DIMENSION (NINTR) ::
     >   XINTR, YINTR, ZINTR   ! Nacelle/pylon intersection coordinates ...

      INTEGER, INTENT (IN), DIMENSION (NINTR) ::
     >   IINTR, JINTR          ! ... and associated surface cell indices

      INTEGER, INTENT (IN) ::
     >   LENXB                 ! Length of X/Y/ZBASE(*) arrays

      REAL, INTENT (INOUT), DIMENSION (LENXB) ::
     >   XBASE, YBASE, ZBASE

      INTEGER, INTENT (INOUT) ::
     >   LASTXB,               ! Last element used so far by finished patches
     >   MOSTXB                !   "    "     needed for interim workspace

      INTEGER, INTENT (INOUT), DIMENSION (MXPATCH) ::
     >   IXB                   ! IXB(L) = 1st element of patch L in X/Y/ZBASE(*)

      INTEGER, INTENT (IN) ::  ! Subpatch info.
     >   MXSUBI, MXSUBJ

      INTEGER, INTENT (OUT) ::
     >   NSUBJ(MXPATCH),
     >   ISUB(0:MXSUBI,MXPATCH),
     >   JSUB(0:MXSUBJ,MXPATCH)

      INTEGER, INTENT (INOUT) ::
     >   IMXCOM(MXPATCH), JMXCOM(MXPATCH) ! Patch dimensions

      INTEGER, INTENT (IN) ::
     >   IWRIT

C     Local constants:

      REAL,    PARAMETER :: HALF = 0.5, ONE = 1., ZERO = 0.
      LOGICAL, PARAMETER :: NORMALIZE = .TRUE.

C     Local variables:

      INTEGER
     >   I, IEND, IENDL, IENDU, IER, II, ILE, ILEAD, ILEFT,
     >   IX, K, KL, KLEAD, KLEFT, L, L1, L2, L3, L4, L5,
     >   NK, NKU, NUMI, NUMK
      REAL
     >   DS, R, STOTAL, SFINL(KLNAC+1), XKEY, YKEY, ZKEY

      REAL, DIMENSION (KLNAC) ::
     >   SN, XN, YN, ZN

C     Execution:

      NK = KLNAC  ! Before insertion of intersection (decremented in SURFER)
      KL = NK + 1 ! After ...

      IF (ID == 0) THEN ! Plain nacelle

         L = KPATCH

C        Check for room for both halves:

         CALL XBCHECK (LENXB, MXPATCH, L, LASTXB, 2*NHALF, NK, IN,
     >                 IWRIT, MOSTXB)

         DO K = L, L + 1
            IMXCOM(K) = NHALF
            ISUB(1,K) = NHALF
            JMXCOM(K) = NK
            JSUB(1,K) = NK
         END DO

         II = LASTXB + 1
         IXB(L) = II   ! Interior

         DO K = 1, NK
            DO I = 1, NHALF
               XBASE(II) = XNAC(I,K)
               YBASE(II) = YNAC(I,K)
               ZBASE(II) = ZNAC(I,K)
               II = II + 1
            END DO
         END DO

         IXB(L+1) = II ! Exterior

         DO K = 1, NK
            DO I = NHALF, ILNAC
               XBASE(II) = XNAC(I,K)
               YBASE(II) = YNAC(I,K)
               ZBASE(II) = ZNAC(I,K)
               II = II + 1
            END DO
         END DO

         LASTXB = II - 1
         KPATCH = L + 2 ! Next available

         GO TO 990 ! Rather than indent the main stuff

      END IF

C     Pylon/diverter intersection cases:
C     ----------------------------------

C     Insert nacelle sections at the intersection leading & trailing edges,
C     and form the panels along each side of the intersection:

      L1 = KPATCH ! Nacelle interior
      L2 = L1 + 1 ! Exterior, forward of intersection
      L3 = L2 + 1 ! Exterior, below/left of intersection
      L4 = L3 + 1 ! Exterior, above/right of intersection
      L5 = L4 + 1 ! Aft of intersection

      DO ILE = NINTR / 3, NINTR         ! Make sure of the true LE
         IF (XINTR(ILE+1) > XINTR(ILE)) EXIT
      END DO

      ILEAD = IINTR(ILE) + 1 ! # Points "left" of intersection (incl. interior)
      KLEAD = JINTR(ILE) + 1 ! "    "   in patch "below" the intersection
      NUMK  = KL - KLEAD + 1 ! "    "    "    "  "above"

      IF (ITL == 1 .AND. ITU == NINTR) THEN ! TE must be closed in the paneling:
         IEND = MAX (IINTR(1), IINTR(NINTR)) + 1
      ELSE
         IEND = ILNAC
      END IF

      IMXCOM(L1) = NHALF
      IMXCOM(L2) = ILEAD - NHALF + 1
      IMXCOM(L3) = ILE - ITL + 1
      IMXCOM(L4) = ITU - ILE + 1

      DO L = L1, L2
         ISUB(1,L) = IMXCOM(L)
         NSUBJ(L)  = 2
         JSUB(1,L) = KLEAD
         JSUB(2,L) = KL
         JMXCOM(L) = KL
      END DO

      ISUB(1,L3) = IMXCOM(L3)
      JMXCOM(L3) = KLEAD
      JSUB(1,L3) = KLEAD
      ISUB(1,L4) = IMXCOM(L4)
      JMXCOM(L4) = NUMK
      JSUB(1,L4) = NUMK

      IX = LASTXB
      DO L = L1, L4
         IXB(L) = IX + 1
         NUMI   = IMXCOM(L) * JMXCOM(L)
         IX     = IX + NUMI
      END DO
      NUMI = IX - LASTXB

      CALL XBCHECK (LENXB, MXPATCH, L4, LASTXB, NUMI, 1, IN,
     >              IWRIT, MOSTXB)

      LASTXB = IX

C     Construct the patch below (left of) the pylon:
C     ----------------------------------------------

      DS = ONE / REAL (KLEAD - 1)

      DO K = 1, KLEAD - 1
         SFINL(K) = DS * REAL (K - 1)
      END DO

      SFINL(KLEAD) = ONE
      NUMI = IMXCOM(L3)
      IX = IXB(L3)
      II = 1

      DO I = ILE, ITL, -1

         XKEY  = XINTR(I)
         YKEY  = YINTR(I)
         ZKEY  = ZINTR(I)
         ILEFT = IINTR(I)
         KLEFT = JINTR(I)

C        Close the trailing edge for (u,v) definition on the nacelle:

         IF (I == 1) THEN
            IF (ITU == NINTR) THEN
               XKEY = (XKEY + XINTR(NINTR)) * HALF
               YKEY = (YKEY + YINTR(NINTR)) * HALF
               ZKEY = (ZKEY + ZINTR(NINTR)) * HALF
            END IF
         END IF

C        Interpolate a new barrel at the intersection point:

         R = (XKEY - XNAC(ILEFT,KLEFT)) /
     >       (XNAC(ILEFT+1,KLEFT) - XNAC(ILEFT,KLEFT))

         DO K = 1, KLEFT
            XN(K) = (ONE - R) * XNAC(ILEFT,K) + R * XNAC(ILEFT+1,K)
            YN(K) = (ONE - R) * YNAC(ILEFT,K) + R * YNAC(ILEFT+1,K)
            ZN(K) = (ONE - R) * ZNAC(ILEFT,K) + R * ZNAC(ILEFT+1,K)
         END DO

         KLEFT = KLEFT + 1
         XN(KLEFT) = XKEY
         YN(KLEFT) = YKEY
         ZN(KLEFT) = ZKEY

C        Respline the lower barrel:

         CALL NBARREL (II, NUMI, KLEAD, KLEFT,
     >                 XN, YN, ZN, SN, KLEAD, SFINL,
     >                 XBASE(IX), YBASE(IX), ZBASE(IX))
         II = II + 1

      END DO ! Next I

      IENDL = ILEFT

C     Construct the patch above (right of) the pylon:
C     -----------------------------------------------

      DS = ONE / REAL (NUMK - 1)
      DO K = 1, NUMK - 1
         SFINL(K) = DS * REAL (K - 1)
      END DO

      SFINL(NUMK) = ONE
      NUMI = IMXCOM(L4)
      IX = IXB(L4)
      II = 1

      DO I = ILE, ITU

         XKEY  = XINTR(I)
         YKEY  = YINTR(I)
         ZKEY  = ZINTR(I)
         ILEFT = IINTR(I)
         KLEFT = JINTR(I)

         IF (I == NINTR) THEN
            IF (ITL == 1) THEN
               XKEY = (XKEY + XINTR(1)) * HALF
               YKEY = (YKEY + YINTR(1)) * HALF
               ZKEY = (ZKEY + ZINTR(1)) * HALF
            END IF
         END IF

         R = (XKEY - XNAC(ILEFT,KLEFT)) /
     >       (XNAC(ILEFT+1,KLEFT) - XNAC(ILEFT,KLEFT))

         XN(KLEFT) = XKEY
         YN(KLEFT) = YKEY
         ZN(KLEFT) = ZKEY

         DO K = KLEFT + 1, NK
            XN(K) = (ONE - R) * XNAC(ILEFT,K) + R * XNAC(ILEFT+1,K)
            YN(K) = (ONE - R) * YNAC(ILEFT,K) + R * YNAC(ILEFT+1,K)
            ZN(K) = (ONE - R) * ZNAC(ILEFT,K) + R * ZNAC(ILEFT+1,K)
         END DO

         NKU = NK - KLEFT + 1

C        Respline the upper barrel:

         CALL NBARREL (II, NUMI, NUMK, NKU,
     >                 XN(KLEFT), YN(KLEFT), ZN(KLEFT), SN(KLEFT),
     >                 NUMK, SFINL, XBASE(IX), YBASE(IX), ZBASE(IX))
         II = II + 1

      END DO ! Next I

      IENDU = ILEFT
      IEND = MAX (IENDL, IENDU) + 1


C     Load the aft edge of the forward exterior nacelle patch:
C     --------------------------------------------------------

      II = IMXCOM(L2)

      CALL XBCOPY2 (LENXB, XBASE, YBASE, ZBASE, IXB, IMXCOM,
     >              L3, 1, 1, 1, KLEAD, L2, II, 1)

      CALL XBCOPY2 (LENXB, XBASE, YBASE, ZBASE, IXB, IMXCOM,
     >              L4, 1, 1, 2, NUMK, L2, II, KLEAD + 1)

C     Normalized arcs for this barrel through the intersection LE:

      IX = IXB(L2)

      CALL CHORDSRF (II, KL, II, 1, KL, XBASE(IX), YBASE(IX),
     >               ZBASE(IX), NORMALIZE, STOTAL, SFINL)

C     Panel the interior:
C     -------------------

      IX = IXB(L1)

      DO I = 1, NHALF

         DO K = 1, NK
            XN(K) = XNAC(I,K)
            YN(K) = YNAC(I,K)
            ZN(K) = ZNAC(I,K)
         END DO

         CALL NBARREL (I, NHALF, KL, NK, XN, YN, ZN, SN,
     >                 KL, SFINL, XBASE(IX), YBASE(IX), ZBASE(IX))
      END DO

C     Likewise for the forward exterior:

      NUMI = IMXCOM(L2)
      IX   = IXB(L2)
      II   = 1

      DO I = NHALF, ILEAD - 1

         DO K = 1, NK
            XN(K) = XNAC(I,K)
            YN(K) = YNAC(I,K)
            ZN(K) = ZNAC(I,K)
         END DO

         CALL NBARREL (II, NUMI, KL, NK, XN, YN, ZN, SN,
     >                 KL, SFINL, XBASE(IX), YBASE(IX), ZBASE(IX))
         II = II + 1

      END DO

      KPATCH = L5 ! Next empty patch, unless there's an aft panel


C     Panel the aft nacelle section, if any:
C     --------------------------------------

      NUMI = ILNAC - IEND + 2 ! Careful; includes pt. inserted at pylon TE

      IF (NUMI > 1 .AND. .NOT. HSCT) THEN

         II = NUMI * KL

         CALL XBCHECK (LENXB, MXPATCH, L5, LASTXB, II, 1, IN,
     >                 IWRIT, MOSTXB)

         IMXCOM(L5) = NUMI
         ISUB(1,L5) = NUMI
         NSUBJ(L5)  = 2
         JSUB(1,L5) = KLEAD
         JSUB(2,L5) = KL
         JMXCOM(L5) = KL

         IX      = LASTXB + 1
         IXB(L5) = IX
         LASTXB  = LASTXB + II
         KPATCH  = L5 + 1

C        First barrel of the aft nacelle patch at the intersection TE:

         II = IMXCOM(L3)

         CALL XBCOPY2 (LENXB, XBASE, YBASE, ZBASE, IXB, IMXCOM,
     >                 L3, II, II, 1, KLEAD, L5, 1, 1)

         II = IMXCOM(L4)

         CALL XBCOPY2 (LENXB, XBASE, YBASE, ZBASE, IXB, IMXCOM,
     >                 L4, II, II, 2, NUMK, L5, 1, KLEAD + 1)

C        Normalized arcs for this barrel through the intersection TE:

         CALL CHORDSRF (NUMI, KL, 1, 1, KL, XBASE(IX), YBASE(IX),
     >                  ZBASE(IX), NORMALIZE, STOTAL, SFINL)

         II = IEND

         DO I = 2, NUMI

            DO K = 1, NK
               XN(K) = XNAC(II,K)
               YN(K) = YNAC(II,K)
               ZN(K) = ZNAC(II,K)
            END DO

            II = II + 1

            CALL NBARREL (I, NUMI, KL, NK, XN, YN, ZN, SN,
     >                    KL, SFINL, XBASE(IX), YBASE(IX), ZBASE(IX))
         END DO

      END IF

  990 RETURN

      END SUBROUTINE NACSURF

C***********************************************************************
C
      SUBROUTINE NACSURF2 (J, KPATCH, MXPATCH, IDIM, KDIM,
     >                     XNAC, YNAC, ZNAC,
     >                     IL, KL, NHALF, LWING, CLOSEDIV, EPS,
     >                     LENXB, XBASE, YBASE, ZBASE, LASTXB, MOSTXB,
     >                     IXB, MXSUBI, MXSUBJ, NSUBJ, ISUB, JSUB,
     >                     IMXCOM, JMXCOM, IWRIT, FAIL)
C
C     NACSURF2 deals with the aft portion of an HSCT-type nacelle which
C     extends beyond the wing trailing edge.  It uses 2-D parametric
C     methods, unlike the bulk of the nacelle paneling.
C
C     03/25/98  DAS  Initial implementation.
C     11/01/98   "   Aft cap is one piece, not two, now that WDSURF
C                    adjusts upper aft wing patches to match lower.
C     03/10/99   "   Packed patch version.
C     05/14/99   "   Subpatch version.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   J                     ! Nacelle number
      INTEGER, INTENT (INOUT) ::
     >   KPATCH                ! In & out with the next available patch
      INTEGER, INTENT (IN) ::
     >   MXPATCH,
     >   IDIM, KDIM           ! Dimensions of input nacelle arrays
      REAL, INTENT (IN), DIMENSION (IDIM,KDIM) ::
     >   XNAC, YNAC, ZNAC      ! Nacelle coordinates in real space
      INTEGER, INTENT (IN) ::
     >   IL, KL,               ! Active (plain) nacelle dimensions
     >   NHALF,                ! UEVAL/VEVAL needed on exterior only
     >   LWING                 ! Upper wing patch to match at wing TE
      REAL, INTENT (IN) ::
     >   CLOSEDIV(3,2),        ! (x,y,z)s where diverter meets wing TE
     >   EPS
      INTEGER, INTENT (IN) ::
     >   LENXB
      REAL, INTENT (INOUT), DIMENSION (LENXB) ::
     >   XBASE, YBASE, ZBASE   ! Packed patches
      INTEGER, INTENT (INOUT) ::
     >   LASTXB,               ! Last element used so far by finished patches
     >   MOSTXB                !   "    "     needed for interim workspace
      INTEGER, INTENT (INOUT) ::
     >   IXB(MXPATCH)          ! IXB(L) = 1st element of patch L in X/Y/ZBASE(*)
      INTEGER, INTENT (IN) ::  ! Subpatch info.
     >   MXSUBI, MXSUBJ
      INTEGER, INTENT (OUT) ::
     >   NSUBJ(MXPATCH),
     >   ISUB(0:MXSUBI,MXPATCH),
     >   JSUB(0:MXSUBJ,MXPATCH)
      INTEGER, INTENT (INOUT), DIMENSION (MXPATCH) ::
     >   IMXCOM, JMXCOM        ! Patch dimensions
      INTEGER, INTENT (IN) ::
     >   IWRIT
      LOGICAL, INTENT (OUT) ::
     >   FAIL

C     Local constants:

      INTEGER,   PARAMETER :: KASE1 = 1, KASE3 = 3, ! See PLBICUT
     >                        LINUV = 1  ! See PLINCUB (linear in the U dirn.)
      REAL,      PARAMETER :: ONE = 1., ZERO = 0.
      LOGICAL,   PARAMETER :: CLOSED = .FALSE.
      CHARACTER, PARAMETER :: LOOSE * 1 = 'B' ! Loose cubics

C     Local variables:

      INTEGER
     >   I, I1, IER, IP, ISMAT, ITE, IX, IX1, IX2, IX3, JW, JCALL,
     >   K, K1, K2, K3, K9, KASE, KCROWN, KP, KSMAT, KTE,
     >   L, LI, LO, LW, LUNERR, NI, NJ, MI, NUMI, NUMK
      REAL
     >   D1, D2, D3, DRANGE, P, PM1, Q, QM1, RN, SMATRIX(4,4,3),
     >   UCUT, VCUT, UVRANGE(4), XYZCUT(3)
      REAL, DIMENSION (IDIM,KDIM) ::
     >   UNAC, VNAC
      REAL, DIMENSION (NHALF,KL) ::
     >   UEVAL, VEVAL

      SAVE UVRANGE
      DATA UVRANGE /ZERO, ONE, ZERO, ONE/ ! See PLBICUT

C     Execution:

      L  = KPATCH ! Next available patch
      LO = L - 1  ! Outboard under-wing nacelle patch from NACSURF
      LI = L - 2  ! Inboard ...
      LW = LWING  ! LW is the relevant upper wing surface patch

C     Determine the I range of the aft nacelle to work with:

      NI = IMXCOM(LW)   ! Indices of wing patch TE at nacelle centerline
      NJ = JMXCOM(LW)
      IX = IXB(LW) + (NJ / 2) * NI - 1 ! (NI,NJ/2,LW)
      I1 = IL / 4       ! Estimate of wing TE index relative to nacelle LE
      KCROWN = KL / 2   ! Nacelle crown index
      NUMI = IL / 2 - 2 ! Monotonic number of nacelle pts. to search
      I = IL - NUMI + 1 ! Corresp. initial index

      CALL INTERVAL (NUMI, XNAC(I,KCROWN), XBASE(IX), ONE, I1)

      I1 = I1 + I - 1   ! Wing TE index in fully wrapped nacelle surface
      I1 = I1 - 5       ! Safety margin to include diverter pts. at TE
      K9 = KL / 4       ! Azimuthal range of upper aft nacelle (~9 & 3 o'clock)
      K3 = 3 * K9

C     Parameterize the aft nacelle:

      UNAC(I1,1) = ZERO ! -999. would suppress normalization

      CALL PARAM2D (IDIM, KDIM, I1, IL, 1, KL, XNAC, YNAC, ZNAC,
     >              UNAC, VNAC)

      DRANGE = XNAC(IL,1) - XNAC(IL/2,1) ! Data range
      DRANGE = MAX (DRANGE, XNAC(IL/2,1), XNAC(IL,1)) ! Similar to INTSEC4 mod.
      LUNERR = -ABS (IWRIT)
      ISMAT  = 0
      KSMAT  = 0


C     Inboard lower aft nacelle subpatch:
C     -----------------------------------

C     Locate the wing/diverter trailing edge pt. in the nacelle surface mesh:

      XYZCUT(1:3) = CLOSEDIV(1:3,1)

      IP = I1 + 5
      KP = (K9 + KCROWN) / 2

      UCUT = UNAC(IP,KP)
      VCUT = VNAC(IP,KP)

      CALL PLBICUT (KASE3, IDIM, KDIM, I1, IL, K9, KCROWN,
     >              XNAC, YNAC, ZNAC, UNAC, VNAC, IP, KP,
     >              ISMAT, KSMAT, SMATRIX, EPS, DRANGE,
     >              XYZCUT, UVRANGE, UCUT, VCUT, P, Q, LUNERR, IER)

      IF (IER > 3) THEN
         JW = 99
         JCALL = 1
         GO TO 800
      END IF

      PM1  = ONE - P
      QM1  = ONE - Q
      ITE  = IP
      KTE  = KP
      NUMI = IL - IP + 2 ! Any reasonable number of aft nacelle pts. streamwise
      NUMK = JMXCOM(LI)  ! From inboard nacelle patch under the wing
      MI   = IMXCOM(LI)

C     Check for exceeding LENXB:

      CALL XBCHECK (LENXB, MXPATCH, L, LASTXB, NUMI, NUMK, J, IWRIT,
     >              MOSTXB)

      IMXCOM(L) = NUMI
      ISUB(1,L) = NUMI
      JSUB(1,L) = NUMK
      IXB(L)    = LASTXB + 1
      LASTXB    = LASTXB + NUMI * NUMK
      IX        = IXB(L) + (NUMK - 1) * NUMI ! (1,NUMK,L)
      XBASE(IX) = CLOSEDIV(1,1)
      YBASE(IX) = CLOSEDIV(2,1)
      ZBASE(IX) = CLOSEDIV(3,1)

      UEVAL(1,NUMK)    = UCUT
      VEVAL(1,NUMK)    = VCUT
      UEVAL(NUMI,NUMK) = ONE
      VEVAL(NUMI,NUMK) = QM1 * VNAC(IL,KP) + Q * VNAC(IL,KP+1)

C     (u,v)s along the edge below the wing trailing edge and the nacelle exit:

      IX   = IXB(L)
      IX1  = IXB(LI) + MI - 1 ! (MI,1,LI)
      D3   = VEVAL(NUMI,NUMK) / REAL (NUMK - 1)
      KP   = 1
      K1   = 1
      K2   = K9
      UCUT = UNAC(IP,1)
      VCUT = ZERO

      DO K = 1, NUMK - 1

         XYZCUT(1) = XBASE(IX1) ! (MI,K,LI)
         XYZCUT(2) = YBASE(IX1)
         XYZCUT(3) = ZBASE(IX1)

         XBASE(IX) = XYZCUT(1) ! (1,K,L)
         YBASE(IX) = XYZCUT(2)
         ZBASE(IX) = XYZCUT(3)
         IX = IX + NUMI

         IX2 = IX1 + MI ! Next K
         D1  = YBASE(IX1) - YBASE(IX2)
         D2  = ZBASE(IX1) - ZBASE(IX2)
         IX1 = IX2

         IF (ABS (D1) > ABS (D2)) THEN
            KASE = KASE1
         ELSE
            KASE = KASE3
         END IF

         IF (KP + 3 >= K9) THEN ! Try to avoid top/bottom confusion
            K1 = KP - 2
            K2 = KCROWN
         END IF

         CALL PLBICUT (KASE, IDIM, KDIM, I1, IL, K1, K2,
     >                 XNAC, YNAC, ZNAC, UNAC, VNAC, IP, KP,
     >                 ISMAT, KSMAT, SMATRIX, EPS, DRANGE,
     >                 XYZCUT, UVRANGE, UCUT, VCUT, P, Q, LUNERR, IER)

         IF (IER > 3) THEN
            JW = K
            JCALL = 2
            GO TO 800
         END IF

         UEVAL(1,K) = UCUT
         VEVAL(1,K) = VCUT

         UEVAL(NUMI,K) = ONE
         VEVAL(NUMI,K) = D3 * REAL (K - 1)

      END DO

C     Axial edges:

      RN = ONE / REAL (NUMI - 1)
      D1 = RN * (ONE - UEVAL(1,1))
      D2 = RN * (ONE - UEVAL(1,NUMK))
      D3 = RN * (VEVAL(NUMI,NUMK) - VEVAL(1,NUMK))

      DO I = 2, NUMI - 1
         VEVAL(I,1)    = ZERO
         UEVAL(I,1)    = UEVAL(1,1)    + D1 * REAL (I - 1)
         UEVAL(I,NUMK) = UEVAL(1,NUMK) + D2 * REAL (I - 1)
         VEVAL(I,NUMK) = VEVAL(1,NUMK) + D3 * REAL (I - 1)
      END DO

C     UEVAL(*,KL-5) serves as TFI workspace:

      K = KL - 5

      CALL TFI2D (NHALF, 1, NUMI, 1, NUMK, UEVAL, VEVAL, UEVAL(1,K))

C     Interpolate interior (x,y,z)s and some edges:

      IX1 = IXB(L) + 1 ! (2,1,L)
      KP  = 1
      K1  = 1
      K2  = K9

      DO K = 1, NUMK

         IX  = IX1
         IX1 = IX1 + NUMI
         IP  = ITE

         IF (KP + 3 >= K9) THEN
            K1 = KP - 2
            K2 = KCROWN
         END IF

         DO I = 2, NUMI

            CALL PLINCUB (IDIM, KDIM, I1, IL, K1, K2,
     >                    XNAC, YNAC, ZNAC, UNAC, VNAC,
     >                    UEVAL(I,K), VEVAL(I,K), LINUV, LOOSE,
     >                    CLOSED, IP, KP, EPS, XBASE(IX), 
     >                    YBASE(IX), ZBASE(IX), IER)
            IX = IX + 1
         END DO
      END DO

C     Aft nacelle cap:
C     ----------------

      CALL XBCHECK (LENXB, MXPATCH, L, LASTXB, NUMI, NJ, J, IWRIT,
     >              MOSTXB)

      JSUB(2,L) = NUMK + NJ - 1
      IX        = LASTXB + 1  ! (1,2) element of cap subpatch before it's lost
      IX1       = IX + 1      ! (2,2)   "    "    "
      LASTXB    = LASTXB + NUMI * (NJ - 1)

C     Copy the base of the upper surface triangle patch:

      CALL XBCOPY2 (LENXB, XBASE, YBASE, ZBASE, IXB, IMXCOM,
     >              LW, NI, NI, 1, NJ, L, 1, NUMK)

      IP = ITE
      KP = KTE
      UCUT = UNAC(IP,KP)
      VCUT = VNAC(IP,KP)

C     Set up (u,v)s for the cap subpatch:

      DO K = 2, NJ

         XYZCUT(1) = XBASE(IX) ! Trailing edge pt.
         XYZCUT(3) = ZBASE(IX)
         IX = IX + NUMI

         CALL PLBICUT (KASE3, IDIM, KDIM, I1, IL, K9, K3,
     >                 XNAC, YNAC, ZNAC, UNAC, VNAC, IP, KP,
     >                 ISMAT, KSMAT, SMATRIX, EPS, DRANGE,
     >                 XYZCUT, UVRANGE, UCUT, VCUT, P, Q, LUNERR, IER)

         IF (IER > 3) THEN
            JW = K
            JCALL = 3
            GO TO 800
         END IF

         QM1 = ONE - Q
         UEVAL(1,K) = UCUT
         VEVAL(1,K) = VCUT
         UEVAL(NUMI,K) = ONE
         VEVAL(NUMI,K) = QM1 * VNAC(IL,KP) + Q * VNAC(IL,KP+1)

         D1 = RN * (ONE - UCUT)
         D2 = RN * (VEVAL(NUMI,K) - VCUT)

         DO I = 2, NUMI - 1
            UEVAL(I,K) = UCUT + D1 * REAL (I - 1)
            VEVAL(I,K) = VCUT + D2 * REAL (I - 1)
         END DO

      END DO

      ITE = IP
      KP  = KTE

      DO K = 2, NJ ! Evaluate remaining (x,y,z)s

         IX  = IX1
         IX1 = IX1 + NUMI
         IP  = ITE

         DO I = 2, NUMI

            CALL PLINCUB (IDIM, KDIM, I1, IL, K9, K3,
     >                    XNAC, YNAC, ZNAC, UNAC, VNAC,
     >                    UEVAL(I,K), VEVAL(I,K), LINUV, LOOSE,
     >                    CLOSED, IP, KP, EPS, XBASE(IX), 
     >                    YBASE(IX), ZBASE(IX), IER)
            IX = IX + 1
         END DO
      END DO

C     Outboard aft subpatch below the level of the wing TE:
C     -----------------------------------------------------

      MI   = IMXCOM(LO)
      NUMK = JMXCOM(LO)

      CALL XBCHECK (LENXB, MXPATCH, L, LASTXB, NUMI, NUMK, J, IWRIT,
     >              MOSTXB)

      IX        = LASTXB + 1  ! (1,2) element of subpatch
      IX3       = IX          ! Needed again later
      LASTXB    = LASTXB + NUMI * (NUMK - 1)
      JMXCOM(L) = JSUB(2,L) + NUMK - 1
      JSUB(3,L) = JMXCOM(L)
      NSUBJ(L)  = 3

      DO I = 1, NUMI ! Subpatch edge
         UEVAL(I,1) = UEVAL(I,NJ)
         VEVAL(I,1) = VEVAL(I,NJ)
      END DO

C     (u,v)s along the edge below the wing trailing edge and the nacelle exit:

      D3 = (ONE - VEVAL(NUMI,1)) / REAL (NUMK - 1)
      IP = ITE
      K1 = KCROWN
      K2 = K3
      UCUT = UNAC(IP,KP)
      VCUT = VNAC(IP,KP)
      IX1  = IXB(LO) + MI + MI - 1 ! (MI,2,LO)

      DO K = 2, NUMK

         XYZCUT(1) = XBASE(IX1) ! (MI,K,LO)
         XYZCUT(2) = YBASE(IX1)
         XYZCUT(3) = ZBASE(IX1)

         XBASE(IX) = XYZCUT(1) ! (1,K) of subpatch
         YBASE(IX) = XYZCUT(2)
         ZBASE(IX) = XYZCUT(3)
         IX        = IX + NUMI

         IX2 = IX1 - MI ! (MI,K-1,LO)
         D1  = YBASE(IX1) - YBASE(IX2)
         D2  = ZBASE(IX1) - ZBASE(IX2)
         IX1 = IX1 + MI

         IF (ABS (D1) > ABS (D2)) THEN
            KASE = KASE1
         ELSE
            KASE = KASE3
         END IF

         IF (KP + 3 >= K2) THEN
            K1 = KP - 2
            K2 = KL
         END IF

         CALL PLBICUT (KASE, IDIM, KDIM, I1, IL, K1, K2,
     >                 XNAC, YNAC, ZNAC, UNAC, VNAC, IP, KP,
     >                 ISMAT, KSMAT, SMATRIX, EPS, DRANGE,
     >                 XYZCUT, UVRANGE, UCUT, VCUT, P, Q, LUNERR, IER)

         IF (IER > 3) THEN
            JW = K
            JCALL = 5
            GO TO 800
         END IF

         UEVAL(1,K) = UCUT
         VEVAL(1,K) = VCUT

         UEVAL(NUMI,K) = ONE
         VEVAL(NUMI,K) = VEVAL(NUMI,1) + D3 * REAL (K - 1)

      END DO

C     Keel edge (u,v)s:

      D1 = RN * (ONE - UEVAL(1,NUMK))

      DO I = 2, NUMI - 1
         VEVAL(I,NUMK) = ONE
         UEVAL(I,NUMK) = UEVAL(1,NUMK) + D1 * REAL (I - 1)
      END DO


      K = KL - 5 ! UEVAL(*,KL-5) serves as TFI workspace:

      CALL TFI2D (NHALF, 1, NUMI, 1, NUMK, UEVAL, VEVAL, UEVAL(1,K))

C     Interpolate interior (x,y,z)s and some edges:

      IX1 = IX3 + 1 ! (2,2) of subpatch
      KP  = KTE
      K1  = KCROWN
      K2  = K3

      DO K = 2, NUMK

         IX  = IX1
         IX1 = IX1 + NUMI
         IP  = ITE

         IF (KP + 3 >= K2) THEN
            K1 = KP - 2
            K2 = KL
         END IF

         DO I = 2, NUMI

            CALL PLINCUB (IDIM, KDIM, I1, IL, K1, K2,
     >                    XNAC, YNAC, ZNAC, UNAC, VNAC,
     >                    UEVAL(I,K), VEVAL(I,K), LINUV, LOOSE,
     >                    CLOSED, IP, KP, EPS, XBASE(IX), 
     >                    YBASE(IX), ZBASE(IX), IER)
            IX = IX + 1
         END DO
      END DO

      L = L + 1


      KPATCH = L ! Next empty patch
      FAIL = .FALSE.
      GO TO 999

  800 IF (IWRIT > 0) THEN
         WRITE (IWRIT, '(3(/,A,I2),/,A,3I4,/,A,1P,3E12.3)')
     >      ' NACSURF2: IER from PLBICUT =', IER,
     >      ' Nacelle/diverter # =', J,
     >      ' Call # to PLBICUT  =', JCALL,
     >      ' Target J/K and nacelle IP, KP =', JW, IP, KP,
     >      ' XYZCUT =', XYZCUT
      END IF

  900 FAIL = .TRUE.

  999 RETURN

      END SUBROUTINE NACSURF2

C***********************************************************************
C
      SUBROUTINE NBARREL (I, MXIGRID, MXKGRID, NK, XN, YN, ZN, SN,
     >                    KL, SNORM, XBASE, YBASE, ZBASE)
C
C     NBARREL imposes a normalized distribution on the given portion of a
C     nacelle section, returning results as the Ith row of the given patch.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   I, MXIGRID, MXKGRID, NK, KL
      REAL, INTENT (IN) ::
     >   XN(NK), YN(NK), ZN(NK)
      REAL, INTENT (OUT) ::
     >   SN(NK)
      REAL, INTENT (IN) ::
     >   SNORM(KL)
      REAL, INTENT(OUT), DIMENSION (MXIGRID,MXKGRID) ::
     >   XBASE, YBASE, ZBASE

C     Local constants:

      REAL,      PARAMETER :: DERIVS = -999.   ! Suppress them
      LOGICAL,   PARAMETER :: CLOSED = .FALSE., NORMALIZE = .TRUE.
      CHARACTER, PARAMETER :: METHOD * 1 = 'B' ! Loose cubics

C     Local variables:

      INTEGER
     >   K, KEVAL
      REAL
     >   STOTAL
      LOGICAL
     >   NEW

C     Execution:

      CALL CHORDS3D (NK, XN, YN, ZN, NORMALIZE, STOTAL, SN)

      KEVAL = 1
      NEW   = .TRUE.

      DO K = 1, KL - 1

         CALL PLSCRV3D (NK, XN, YN, ZN, SN, METHOD, NEW, CLOSED,
     >                  SNORM(K), KEVAL, XBASE(I,K), YBASE(I,K),
     >                  ZBASE(I,K), DERIVS)
         NEW = .FALSE.
      END DO

      XBASE(I,KL) = XN(NK) ! Exactly
      YBASE(I,KL) = YN(NK)
      ZBASE(I,KL) = ZN(NK)

      END SUBROUTINE NBARREL

C*******************************************************************************
C
      SUBROUTINE NEWPTS1D (NSUB, ISUB1, X1, Y1, Z1, ISUB2, X2, Y2, Z2,
     >                     METHOD)
C
C        NEWPTS1D is a generalization of CHANGEN for adjusting the point counts
C     within multiple subcurves of a curve defined by points in 3-space.  It
C     avoids unnecessary boundary condition effects which use of CHANGEN on each
C     subcurve would introduce.  All degenerate cases of no change are handled.
C
C        The diagram illustrates a case with 3 subcurves, where the "critical"
C     points marked o are preserved precisely, and the numbers of + points are
C     modified via interpolation within the entire input curve:
C
C     Input:     o-+-+-+-+-+-o-+-+-+-+-+-+-+-+-o-+-+-o
C
C     Output:    o--+--+--+--o--+--+--+--+--+--o+++++o
C
C     04/14/99  DAS  Initial adaptation of CHANGEN.
C
C     Author:  David Saunders, Raytheon/NASA Ames, Moffett Field, CA
C
C*******************************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   NSUB,                ! Number of subcurves >= 1
     >   ISUB1(0:NSUB)        ! Indices of critical points on input curve,
                              ! where ISUB1(0) = 1

      REAL, INTENT (IN), DIMENSION (ISUB1(NSUB)) ::
     >   X1, Y1, Z1           ! Input data points

      INTEGER, INTENT (IN) ::
     >   ISUB2(0:NSUB)        ! Critical point indices for output curve,
                              ! where ISUB2(0) = 1

      REAL, INTENT (OUT), DIMENSION (ISUB2(NSUB)) ::
     >   X2, Y2, Z2           ! Output data points

      CHARACTER, INTENT (IN) ::
     >   METHOD*1             ! Type of fit to be used by PLSCRV3D, q.v.

C     Procedures:

      EXTERNAL
     >   CHORDS3D,            ! Arc-length utility
     >   PLSCRV3D             ! 3-space local cubic spline utility

C     Local constants:

      REAL,    PARAMETER :: DERIVS = -999.,   ! Suppresses them
     >                      ONE    = 1.
      LOGICAL, PARAMETER :: NORM   = .FALSE., ! No need to normalize
     >                      CLOSED = .FALSE.

C     Local variables:

      INTEGER
     >   I, I1, I2, IA, IB, IEVAL, IP, N, NDATA
      REAL
     >   T(ISUB1(NSUB)), P, R, RI1, RIP, TEVAL, TOTAL
      LOGICAL
     >   NEW

C     Execution:

      NEW   = .TRUE.
      NDATA = ISUB1(NSUB)

      CALL CHORDS3D (NDATA, X1, Y1, Z1, NORM, TOTAL, T)

      DO N = 1, NSUB

         I1 = ISUB1(N-1)
         I2 = ISUB1(N)
         IA = ISUB2(N-1)
         IB = ISUB2(N)

         IF (IB - IA == I2 - I1) THEN ! No change, but we have to transcribe

            IP = I1
            DO I = IA, IB
               X2(I) = X1(IP)
               Y2(I) = Y1(IP)
               Z2(I) = Z1(IP)
               IP = IP + 1
            END DO

         ELSE ! Preserve the nature of the input point distribution

            R     = REAL (I2 - I1) / REAL (IB - IA)
            RI1   = REAL (I1)
            IEVAL = I1

            X2(IA) = X1(I1)
            Y2(IA) = Y1(I1)
            Z2(IA) = Z1(I1)

            DO I = IA + 1, IB - 1

               RIP   = RI1 + R * REAL (I - IA)
               IP    = INT (RIP)
               P     = RIP - REAL (IP)
               TEVAL = (ONE - P) * T(IP) + P * T(IP+1)

               CALL PLSCRV3D (NDATA, X1, Y1, Z1, T,
     >                        METHOD, NEW, CLOSED, TEVAL, IEVAL,
     >                        X2(I), Y2(I), Z2(I), DERIVS)
               NEW = .FALSE.

            END DO

            X2(IB) = X1(I2)
            Y2(IB) = Y1(I2)
            Z2(IB) = Z1(I2)

         END IF

      END DO ! Next subcurve

      END SUBROUTINE NEWPTS1D

C*******************************************************************************
C
      SUBROUTINE NEWPTS2D (NSUBI, NSUBJ, ISUB1, JSUB1, ISUB2, JSUB2,
     >                     IMETHOD, JMETHOD, X1, Y1, Z1, X2, Y2, Z2)
C
C        NEWPTS2D redistributes the rows and/or columns of a packed surface
C     patch containing one or more subpatches, imposing specified point counts
C     within each subpatch.  The character of the input point distribution
C     along each subpatch row or column is preserved.  The degenerate case
C     where all point counts already match is handled here, but avoiding the
C     call at the higher level might be preferable.
C
C     Two-dimensional interpolation methods are avoided, partly because the
C     changes in point counts are expected to be small or zero, and partly
C     because the application may well require two different 1-D methods, as
C     in the case of patches derived from planar sections for wing surfaces,
C     where linear is appropriate spanwise and cubic chordwise.  (Bicubic
C     would be inappropriate, and we don't have hybrid 2-D methods.)
C     Efficiency is also a consideration.
C
C     03/23/99  DAS  Surface subpatch utility wrapped around NEWPTS1D,
C                    adapted from CHANGEMBYN/CHANGEN (whole patch analogue).
C     10/13/99   "   Singular nose point needs special treatment.
C
C     Author:  David Saunders, Raytheon/NASA Ames, Moffett Field, CA
C
C*******************************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   NSUBI, NSUBJ,      ! Number of subpatches in each direction
     >   ISUB1(0:NSUBI),    ! Indices of subpatch boundaries in the
     >   JSUB1(0:NSUBJ)     ! input patch; ISUB1(0) = JSUB1(0) = 1

      INTEGER, INTENT (IN) ::
     >   ISUB2(0:NSUBI),    ! Indices of subpatch boundaries desired in
     >   JSUB2(0:NSUBJ)     ! the output patch; ISUB2(0) = JSUB2(0) = 1

      CHARACTER * 1, INTENT (IN) ::
     >   IMETHOD, JMETHOD   ! Types of fit to be used by PLSCRV3D, q.v.

      REAL, INTENT (IN), DIMENSION (ISUB1(NSUBI), JSUB1(NSUBJ)) ::
     >   X1, Y1, Z1         ! Input patch coordinates (packed)

      REAL, INTENT (OUT), DIMENSION (ISUB2(NSUBI), JSUB2(NSUBJ)) ::
     >   X2, Y2, Z2         ! Output patch coordinates (not overlapping X/Y/Z1)

C     Procedures:

      EXTERNAL NEWPTS1D     ! Subcurve utility wrapped around PLSCRV3D

C     Local variables:

      INTEGER
     >   I, J, M1, M2, N1, N2

      REAL, DIMENSION (JSUB1(NSUBJ)) ::
     >   XA, YA, ZA         ! For gathering one input row at a time

      REAL, DIMENSION (JSUB2(NSUBJ)) ::
     >   XB, YB, ZB         ! For one redistributed row at a time

      REAL, ALLOCATABLE, DIMENSION (:, :) ::
     >   XT, YT, ZT         ! For an M2 x N1 intermediate patch

      LOGICAL
     >   SAMEI, SAMEJ, SINGULAR

C     Execution:

      SAMEI = .TRUE.
      SAMEJ = .TRUE.

      DO I = 1, NSUBI
         IF (ISUB2(I) - ISUB2(I-1) /=
     >       ISUB1(I) - ISUB1(I-1)) SAMEI = .FALSE.
      END DO

      DO J = 1, NSUBJ
         IF (JSUB2(J) - JSUB2(J-1) /=
     >       JSUB1(J) - JSUB1(J-1)) SAMEJ = .FALSE.
      END DO

      M1 = ISUB1(NSUBI)
      M2 = ISUB2(NSUBI)
      N1 = JSUB1(NSUBJ)
      N2 = JSUB2(NSUBJ)

      IF (.NOT. SAMEI) THEN

         IF (.NOT. SAMEJ) THEN ! Revised columns will be changed again row-wise

            ALLOCATE (XT(M2,N1), YT(M2,N1), ZT(M2,N1))

            DO J = 1, N1 ! Redistribute in the I direction into interim M2 x N1

               CALL NEWPTS1D (NSUBI, ISUB1,   X1(1,J), Y1(1,J), Z1(1,J),
     >                        ISUB2, XT(1,J), YT(1,J), ZT(1,J), IMETHOD)
            END DO

            DO I = 1, M2 ! Redistribute each interim line in the J direction

               DO J = 1, N1 ! Gather an interim row
                  XA(J) = XT(I,J)
                  YA(J) = YT(I,J)
                  ZA(J) = ZT(I,J)
               END DO

               SINGULAR = XA(1) == XA(2) .AND. YA(1) == YA(2) .AND.
     >                    ZA(1) == ZA(2)

               IF (SINGULAR) THEN
                  DO J = 1, N2
                     XB(J) = XA(1)
                     YB(J) = YA(1)
                     ZB(J) = ZA(1)
                  END DO
               ELSE

                  CALL NEWPTS1D (NSUBJ, JSUB1,  XA, YA, ZA,
     >                           JSUB2, XB, YB, ZB, JMETHOD)
               END IF

               DO J = 1, N2
                  X2(I,J) = XB(J)
                  Y2(I,J) = YB(J)
                  Z2(I,J) = ZB(J)
               END DO

            END DO

            DEALLOCATE (XT, YT, ZT)

         ELSE ! SAMEJ but .NOT. SAMEI; avoid the interim patch

            DO J = 1, N1 ! The new columns can go into place directly

               CALL NEWPTS1D (NSUBI, ISUB1,   X1(1,J), Y1(1,J), Z1(1,J),
     >                        ISUB2, X2(1,J), Y2(1,J), Z2(1,J), IMETHOD)
            END DO

         END IF

      ELSE IF (.NOT. SAMEJ) THEN ! SAMEI; redistribute in J direction directly

         DO I = 1, M1

            DO J = 1, N1 ! Gather an interim row
               XA(J) = X1(I,J)
               YA(J) = Y1(I,J)
               ZA(J) = Z1(I,J)
            END DO

            SINGULAR = XA(1) == XA(2) .AND. YA(1) == YA(2) .AND.
     >                 ZA(1) == ZA(2)

            IF (SINGULAR) THEN
               DO J = 1, N2
                  XB(J) = XA(1)
                  YB(J) = YA(1)
                  ZB(J) = ZA(1)
               END DO
            ELSE

               CALL NEWPTS1D (NSUBJ, JSUB1,  XA, YA, ZA,
     >                        JSUB2, XB, YB, ZB, JMETHOD)
            END IF

            DO J = 1, N2
               X2(I,J) = XB(J)
               Y2(I,J) = YB(J)
               Z2(I,J) = ZB(J)
            END DO

         END DO

      ELSE ! SAMEI and SAMEJ

         DO J = 1, N1
            DO I = 1, M1
               X2(I,J) = X1(I,J)
               Y2(I,J) = Y1(I,J)
               Z2(I,J) = Z1(I,J)
            END DO
         END DO

      END IF

      END SUBROUTINE NEWPTS2D

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
C                    grid sections near a fuselage (CH_GRID/SYN87-SB).
C     06/25/98   "   Discrete choice of leading edge can add irregularity
C                    in the Y direction.  Use continuous techniques.
C     04/03/99   "   Reused in SURFER, which assumes round LE from K1 to K2,
C                    so there's no need to test for sharpness.
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) :: K1, K2       ! Subset of grid stations to treat
      INTEGER, INTENT (IN) :: KL, IL       ! Dimensions of regularized arrays
      INTEGER, INTENT (IN) :: ILE, NHALF   ! Middle index & number either side
      REAL,    INTENT (IN) :: SNORM(NHALF) ! Normalized chordwise distribution
      REAL, INTENT (INOUT) :: XW(IL,KL),   ! Grid stations
     >                        YW(IL,KL), ZW(IL,KL)

C     Local constants:

      INTEGER,   PARAMETER :: LUNERR = 6, NFUNMAX = 20
      REAL,      PARAMETER :: ONE = 1.
      LOGICAL,   PARAMETER :: NEW = .TRUE., NORM = .FALSE.
      CHARACTER, PARAMETER :: LOOSE * 1 = 'B', SUBNAME * 7 = 'RECTIFY'

C     Local variables:

      INTEGER
     >   I, II, ILEAD, ILEFT, ISTAT, ITL, ITU, K, NFUN, NLEAD, NWING
      REAL
     >   DERIVS(2*NHALF), S(IL), SEVAL(2*NHALF), XYZ(IL),
     >   A, B, EPS, SLOWER, SMAX, SUPPER, SLE, TOL, XLE

C     Execution:

      NWING = 2*NHALF - 1
      II    = NWING / 3
      NLEAD = 5 ! 4-pt. method in one of two adjacent intervals
      ITL   = ILE - NHALF + 1
      ITU   = ILE + NHALF - 1
      EPS   = MAX (3.* EPSILON (ONE), 1.E-7)

      DO K = K1, K2

C        Locate input leading edge index, allowing for nonmonotonicity:

         A = XW(II,K)
         DO I = II + 1, II + II
            IF (XW(I,K) < A) THEN
               A = XW(I,K)
               ILEAD = I
            END IF
         END DO

C        Apply the rounded leading edge case of WREGUL, but make it vary
C        smoothly across Ks by not just picking the nearest data point:

         CALL CHORDS3D (NWING, XW(ITL,K), YW(ITL,K), ZW(ITL,K), NORM,
     >                  SMAX, S(ITL))

C        Initialize an iteration for minimizing X w.r.t. S:

         ISTAT = 2
         ILEFT = ILEAD - 2
         A     = S(ILEFT)
         B     = S(ILEAD + 2)
         TOL   = B * EPS
         NFUN  = NFUNMAX

         DO I = 1, NFUNMAX ! Or until convergence

            CALL FMINRC (A, B, SLE, XLE, TOL, NFUN, SUBNAME, -LUNERR,
     >                   ISTAT)

            IF (ISTAT < -1) THEN

               WRITE (LUNERR, '(/, A, 2I4)')
     >            ' RECTIFY: Bad FMINRC return. ISTAT, K =', ISTAT, K
               STOP

            ELSE IF (ISTAT < 0) THEN

               WRITE (LUNERR, '(/, A, I4)')
     >            ' RECTIFY: FMINRC hit NFUNMAX. K =', K

            ELSE IF (ISTAT > 0) THEN ! Evaluate X at S = SLE

               CALL LCSFIT (NLEAD, S(ILEFT), XW(ILEFT,K), NEW, LOOSE,
     >                      1, SLE, XLE, DERIVS)

            ELSE ! ISTAT = 0 - success

               EXIT

            END IF

         END DO

C        Redistribute the section in place:

         SLOWER = SLE
         SUPPER = SMAX - SLE

         DO I = 1, NHALF
            SEVAL(I) = SLOWER * (ONE - SNORM(NHALF + 1 - I))
            SEVAL(I + NHALF - 1) = SLOWER + SUPPER * SNORM(I)
         END DO

         CALL LCSFIT (NWING, S(ITL), XW(ITL,K), NEW, LOOSE, NWING,
     >                SEVAL, XYZ(ITL), DERIVS)

         XW(ITL:ITU,K) = XYZ(ITL:ITU)

         CALL LCSFIT (NWING, S(ITL), YW(ITL,K), NEW, LOOSE, NWING,
     >                SEVAL, XYZ(ITL), DERIVS)

         YW(ITL:ITU,K) = XYZ(ITL:ITU)

         CALL LCSFIT (NWING, S(ITL), ZW(ITL,K), NEW, LOOSE, NWING,
     >                SEVAL, XYZ(ITL), DERIVS)

         ZW(ITL:ITU,K) = XYZ(ITL:ITU)

      END DO

      END SUBROUTINE RECTIFY

C***********************************************************************
C
      SUBROUTINE REPACK (LENXB, XBASE, YBASE, ZBASE, LASTXB, MOSTXB,
     >                   NPATCH, NFUSE, MXSUBI, MXSUBJ, NSUBI, NSUBJ,
     >                   ISUB, JSUB, ISUBO, JSUBO, IXB, IWRIT)
C
C     REPACK adjusts the dimensions of a set of packed surface patches
C     to match some required set of dimensions.  The dimension changes
C     for each patch are expected to be zero most often/small otherwise.
C     The 1-D interpolations are via cubic splines chordwise/piecewise
C     linear spanwise (I and J resp.) for wing/nacelle surfaces.  For
C     fuselage patches NFUSE:NPATCH, they are cubic in both directions.
C
C     This version adjusts point counts for subpatches within patches.
C
C     03/23/99  DAS  Initial implementation (no subpatches).
C     04/14/99   "   Subpatch version.
C
C     Author:  David Saunders, Raytheon/NASA Ames, Moffett Field, CA.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   LENXB                 ! Available space in X/Y/ZBASE(*)
      REAL, INTENT (INOUT), DIMENSION (LENXB) ::
     >   XBASE, YBASE, ZBASE   ! Packed patch coordinates
      INTEGER, INTENT (INOUT) ::
     >   LASTXB,               ! Last element used so far by finished patches
     >   MOSTXB                !   "    "     needed for interim workspace
      INTEGER, INTENT (IN) ::
     >   NPATCH,               ! Number of active patches
     >   NFUSE,                ! Patches NFUSE:NPATCH are fuselage patches
     >   MXSUBI, MXSUBJ        ! Max. # subpatches per patch over all patches
      INTEGER, INTENT (IN), DIMENSION (NPATCH) ::
     >   NSUBI, NSUBJ          ! Numbers of subpatches in each direction
      INTEGER, INTENT (IN), DIMENSION (0:MXSUBI, NPATCH) ::
     >   ISUB,                 ! Input subpatch dimensions = ISUB(0:NSUBI(*),*)
     >   ISUBO                 ! Desired dimensions; ISUB(0,*) = ISUBO(0,*) = 1
      INTEGER, INTENT (IN), DIMENSION (0:MXSUBJ, NPATCH) ::
     >   JSUB,                 ! Analogous subpatch definitions in J direction;
     >   JSUBO                 ! JSUB(0,*) = JSUBO(0,*) = 1
      INTEGER, INTENT (INOUT) ::
     >   IXB(NPATCH)           ! IXB(L) = 1st element of patch L in X/Y/ZBASE(*)
      INTEGER, INTENT (IN) ::
     >   IWRIT                 ! For diagnostics

C     Local constants:

      CHARACTER * 1, PARAMETER ::
     >   LINEAR = 'L', CUBIC * 1 = 'B' ! "Bessel"

C     Local variables:

      INTEGER
     >   I, I1, I2, J, L, M, NADD, NIJ, NSHIFT, NSI, NSJ

      REAL, ALLOCATABLE, DIMENSION (:) ::
     >   XT, YT, ZT

      LOGICAL
     >   SAMEI, SAMEJ

      CHARACTER
     >   JMETHOD * 1

C     Execution:

      JMETHOD = LINEAR ! Until the fuselage patches

      DO L = 1, NPATCH

         NSI = NSUBI(L)
         SAMEI = .TRUE.

         DO I = 1, NSI
            IF (ISUB (I,L) - ISUB (I-1,L) /=
     >          ISUBO(I,L) - ISUBO(I-1,L)) SAMEI = .FALSE.
         END DO

         NSJ = NSUBJ(L)
         SAMEJ = .TRUE.

         DO J = 1, NSJ
            IF (JSUB (J,L) - JSUB (J-1,L) /=
     >          JSUBO(J,L) - JSUBO(J-1,L)) SAMEJ = .FALSE.
         END DO

         IF (SAMEI) THEN
            IF (SAMEJ) CYCLE
         END IF

         NIJ  = ISUBO(NSI,L) * JSUBO(NSJ,L)
         NADD = NIJ - ISUB(NSI,L) * JSUB(NSJ,L)

         IF (NADD > 0) THEN ! Check for sufficient room in X/Y/ZBASE(*)

            CALL XBCHECK (LENXB, NPATCH, L, LASTXB, NADD, 1, L, IWRIT,
     >                    MOSTXB)

         END IF

C        Allocate a temporary patch of the new size:

         ALLOCATE (XT(NIJ), YT(NIJ), ZT(NIJ))

C        Redistribute the wrong-sized patch into the temporary patch:

         I1 = IXB(L)
         IF (L >= NFUSE) JMETHOD = CUBIC

         CALL NEWPTS2D (NSI, NSJ, ISUB(0,L), JSUB(0,L),
     >                  ISUBO(0,L), JSUBO(0,L), CUBIC, JMETHOD,
     >                  XBASE(I1), YBASE(I1), ZBASE(I1), XT, YT, ZT)

C        Shift patches below the target patch to make a right-sized space:

         IF (L < NPATCH) THEN

            NSHIFT = LASTXB - IXB(L+1) + 1

            I2 = I1 + NIJ ! Start of shuffled patches below current patch

            CALL XBCOPY (LENXB, XBASE, YBASE, ZBASE,
     >                   IXB(L+1), NSHIFT, 1, NSHIFT, I2, NSHIFT)

            DO M = L + 1, NPATCH
               IXB(M) = IXB(M) + NADD
            END DO

         END IF

C        Copy the temporary patch into place:

         I = I1
         DO J = 1, NIJ
            XBASE(I) = XT(J)
            YBASE(I) = YT(J)
            ZBASE(I) = ZT(J)
            I = I + 1
         END DO

         DEALLOCATE (XT, YT, ZT)

         LASTXB = LASTXB + NADD
         MOSTXB = MAX (LASTXB, MOSTXB)

      END DO ! Next patch

      END SUBROUTINE REPACK

C*******************************************************************************
C
      SUBROUTINE VFINTR (NCALL, NNCONFIG, IFIN, LOSTINTR,
     >                   KFIN1, KFIN2, MXIFIN, MXKFIN, MXINTR,
     >                   NHALF, ILFIN, XRG, YRG, ZRG,
     >                   NJF, NIF, MXJBODY, MXIBODY, XF, YF, ZF,
     >                   ISTN1, ISTN2, IDEGFUS, NINTR,
     >                   XINTR, YINTR, ZINTR, UINTR, VINTR, TINTR,
     >                   IINTR, JINTR, KINTR, IWRIT, FAIL)
C
C     VFINTR calculates a vertical fin/fuselage intersection given regular
C     geometry sections for each.  The fin may be above or below the fuselage.
C     The "upper" surfaces of its sections (after rotation through 90 degrees)
C     are intersected with the left half of the body, at the crown or keel.
C
C     05/18/99  David Saunders  Adaptation of WFINTR.
C     08/08/00    "      "      FCTV3 requires use of whole body.
C
C*******************************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   NCALL,                     ! -3 means first call for this case
     >   NNCONFIG,                  ! 3 means use the whole body for FCTV3
     >   IFIN                       ! Fin surface # (for LOSTINTR case)
      LOGICAL, INTENT (IN) ::
     >   LOSTINTR                   ! T means [re]store intersection LE data
      INTEGER, INTENT (IN) ::
     >   KFIN1, KFIN2,              ! Subset of wing sections to treat
     >   MXIFIN, MXKFIN,            ! Dimensions of fin section arrays
     >   MXINTR,                    ! Dimension of intersection arrays
     >   NHALF, ILFIN               ! Leading & upper trailing edge fin indices
      REAL, INTENT (IN), DIMENSION (MXIFIN,MXKFIN) ::
     >   XRG, YRG, ZRG              ! Regularized fin sections
      INTEGER, INTENT (IN) ::
     >   NJF,                       ! Active circumferential and axial
     >   NIF,                       ! indices of regular fuselage sections
     >   MXJBODY, MXIBODY           ! J and I dimensions of body section arrays
      REAL, INTENT (IN) ::
     >   XF(MXIBODY),               ! Regularized fuselage sections
     >   YF(MXJBODY,MXIBODY),
     >   ZF(MXJBODY,MXIBODY)
      INTEGER, INTENT (INOUT) ::
     >   ISTN1,                     ! Fuselage station range for this fin,
     >   ISTN2                      ! output on first call, input thereafter
      INTEGER, INTENT (IN) ::
     >   IDEGFUS                    ! 3 means bicubics (INTSEC4) else
                                    ! 1 means bilinear (INTSEC5)
      INTEGER, INTENT (OUT) ::
     >   NINTR                      ! No. of pts. in the intersection
                                    ! (ILFIN, but may change externally)
      REAL, INTENT (OUT), DIMENSION (MXINTR) ::
     >   XINTR, YINTR, ZINTR,       ! Intersection coordinates in NHALF:ILFIN
     >   UINTR,                     ! Corresponding body surface "u" and
     >   VINTR,                     ! "v" values (circumferential & axial)
     >   TINTR                      ! "t" values at the intersection
      INTEGER, INTENT (INOUT), DIMENSION (MXINTR) ::
     >   IINTR, JINTR, KINTR        ! Corresponding body & fin cell indices
      INTEGER, INTENT (IN) ::
     >   IWRIT                      ! For diagnostics if > 0
      LOGICAL, INTENT (OUT) ::
     >   FAIL                       ! Success flag

C     Local constants:

      REAL, PARAMETER :: ONE = 1., ZERO = 0.

C     Local variables:

      INTEGER
     >   I, I1, I2, IFAIL, J, JSTN1, JSTN2
      REAL
     >   R, XT
      REAL, DIMENSION (MXJBODY,MXIBODY) ::
     >   UBODY, VBODY, XBODY ! Dimensions must all match for PARAM2D
      LOGICAL
     >   CROWN, FIRST

C     Local leading edge data saved for LOSTINTR case.  10 >= NNWING

      INTEGER, DIMENSION (10), SAVE :: IINTLE, JINTLE, KINTLE
      REAL,    DIMENSION (10), SAVE :: TINTLE, UINTLE, VINTLE

C     Execution:

      FAIL  = .FALSE. ! Allow retrying in WBINTR
      NINTR = ILFIN   ! But allow this to change in the LOSTINTR case
      FIRST = NCALL == -3

C     Determine the subset of body sections containing the intersection:

      IF (FIRST) THEN
         R  = (XRG(NHALF,KFIN1) - XF(1)) / (XF(NIF) - XF(1))
         I1 = NINT (R * REAL (NIF))
         I2 = I1 + 1
      ELSE
         I1 = ISTN1
         I2 = ISTN2 - 1
      END IF

      XT = MIN (XRG(NHALF,KFIN1), XRG(NHALF,KFIN2))

      CALL INTERVAL (NIF, XF, XT, ONE, I1)

      ISTN1 = I1

      XT = MAX (XRG(ILFIN,KFIN1), XRG(ILFIN,KFIN2))

      CALL INTERVAL (NIF, XF, XT, ONE, I2)

      ISTN2 = MIN (I2 + 1, NIF)

      CROWN = YRG(NHALF,KFIN2) > YF(NJF,I1) ! Fin tip above body crown?

      I1 = NHALF
      I2 = ILFIN

      IF (FIRST) THEN

         IINTR(I1) = ISTN1 + 1

         IF (CROWN) THEN
            JINTR(I1) = NJF - 1
            UINTR(I1) = ONE
         ELSE
            JINTR(I1) = 1
            UINTR(I1) = ZERO
         END IF

         IF (NNCONFIG /= 3) THEN
            VINTR(I1) = (XF(ISTN1+1) - XF(ISTN1)) /
     >                  (XF(ISTN2)   - XF(ISTN1))
         ELSE
            VINTR(I1) = (XF(ISTN1+1) - XF(1)) /
     >                  (XF(NIF)     - XF(1))
         END IF

         KINTR(I1) = (KFIN1 + KFIN2) / 2
         TINTR(I1) = 0.5

      ELSE IF (LOSTINTR) THEN ! Restore leading edge results

         IINTR(I1) = IINTLE(IFIN) 
         JINTR(I1) = JINTLE(IFIN) 
         KINTR(I1) = KINTLE(IFIN) 
         TINTR(I1) = TINTLE(IFIN) 
         UINTR(I1) = UINTLE(IFIN) 
         VINTR(I1) = VINTLE(IFIN) 
         FIRST = .TRUE.

      END IF

      IF (NNCONFIG /= 3) THEN

         IF (CROWN) THEN
            JSTN1 = (4 * NJF) / 5
            JSTN2 = NJF
         ELSE
            JSTN1 = 1
            JSTN2 = NJF / 5
         END IF

      ELSE       ! NNCONFIG == 3: use the whole body
         ISTN1 = 1
         ISTN2 = NIF
         JSTN1 = 1
         JSTN2 = NJF
      END IF

      DO I = ISTN1, ISTN2
         XBODY(JSTN1:JSTN2,I) = XF(I)
      END DO

      UBODY(JSTN1,ISTN1) = -ONE ! Tells INTSEC4/5 to parameterize the body

      CALL WBINTR (MXIFIN, MXKFIN, I1, I2, KFIN1, KFIN2,
     >             XRG, YRG, ZRG, MXJBODY, MXIBODY, JSTN1, JSTN2,
     >             ISTN1, ISTN2, XBODY, YF, ZF, 1, IDEGFUS, UBODY,VBODY,
     >             XINTR, YINTR, ZINTR, UINTR, VINTR, TINTR,
     >             JINTR, IINTR, KINTR, FIRST, IWRIT, IFAIL, FAIL)

      IF (.NOT. FAIL) THEN

         IF (LOSTINTR) THEN ! Save leading edge results for later reuse
            IINTLE(IFIN) = IINTR(I1)
            JINTLE(IFIN) = JINTR(I1)
            KINTLE(IFIN) = KINTR(I1)
            TINTLE(IFIN) = TINTR(I1)
            UINTLE(IFIN) = UINTR(I1)
            VINTLE(IFIN) = VINTR(I1)
         END IF

      END IF

      END SUBROUTINE VFINTR

C***********************************************************************
C
      SUBROUTINE VHINTR (NCALL, IFIN, IW, IWRIT, FAIL) ! etc.
C
C     VHINTR calculates a vertical fin/horizontal tail intersection
C     given a regular geometry mesh for each.
C
C***********************************************************************

      IF (IWRIT > 0) WRITE (IWRIT,*) '??? Dummy VHINTR called.'
      CALL SYNCH
      STOP

      END SUBROUTINE VHINTR

C***********************************************************************
C
      SUBROUTINE VSURF (NCALL, IW, IN, NNSURF, MXIWING, MXKWING,
     >                  MXCRANK, MXINTR, NREGSEC, NHALF, ILWING,KLWING,
     >                  XRG, YRG, ZRG, XINTR, YINTR, ZINTR, TINTR,
     >                  KWING1, KWING2, NCRANK, KCRANK, MCRANK,
     >                  TNORM, TBASE, XBASE, YBASE, ZBASE, IWRIT)
C
C     VSURF imposes a structured mesh on the left surface of a vertical fin
C     given both surfaces of the defining sections in regularized form.
C     Any intersection with a fuselage is input.  Tailplanes aren't handled.
C     Any input crank stations are captured by the mesh.
C
C     05/21/99  D. Saunders  Initial adaptation of WSURF.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   NCALL,              ! -3 means first call
     >   IW,                 ! Vertical fin surface #
     >   IN,                 ! Tip nacelle surface #, or 0 if none (inactive)
     >   NNSURF,             ! Number of wings + fins + nacelles
     >   MXIWING, MXKWING,   ! Dimensions of regularized section arrays
     >   MXCRANK,            ! Max. # cranks handled
     >   MXINTR,             ! Max. # intersection points handled
     >   NREGSEC,            ! Number of regularized fin sections
     >   NHALF,              ! Leading edge index of input stations
     >   ILWING, KLWING      ! Desired paneling dimensions (as if it were
                             ! a wrap-around wing-type panel)

      REAL, INTENT (IN), DIMENSION (MXIWING,MXKWING) ::
     >   XRG, YRG, ZRG       ! Regularized fin defining sections (both surfaces)

      REAL, INTENT (IN), DIMENSION (MXINTR,NNSURF) ::
     >   XINTR, YINTR,       ! Intersection data for all surfaces
     >   ZINTR, TINTR

      INTEGER, INTENT (IN) ::
     >   KWING1(NNSURF),     ! Ranges of stations used for all of the
     >   KWING2(NNSURF),     ! intersections calculations
     >   NCRANK              ! Number of crank stations for this fin

      INTEGER, INTENT (INOUT), DIMENSION (MXCRANK) ::
     >   KCRANK, MCRANK      ! Input and output indices of crank stations

      REAL, INTENT (IN) ::
     >   TNORM(KLWING)       ! Normalized distribution, root to tip

      REAL, INTENT (INOUT) ::
     >   TBASE(KLWING)       ! Denormalized nominal distribution, root to tip

      REAL, INTENT (OUT), DIMENSION (NHALF,KLWING) ::
     >   XBASE, YBASE, ZBASE ! Surface patch (left surface only)

      INTEGER, INTENT (IN) ::
     >   IWRIT               ! For error messages

C     Local constants:

      REAL,      PARAMETER :: ZERO = 0.
      LOGICAL,   PARAMETER :: NEW = .TRUE.
      CHARACTER, PARAMETER :: LINEAR * 1 = 'L'

C     Local variables:

      INTEGER
     >   I, IER, IX, J, K, K1, K2
      REAL
     >   TCRANK(MXCRANK), TEVAL(KLWING), TRG(NREGSEC),
     >   WORK(KLWING * 2), XYZG(NREGSEC,3), XYZM(KLWING,3),
     >   DT, T1, T2, TSCALE, TSHIFT, TTIP, TTOTAL

C     Execution:

C     See WSURF for comments regarding cranks.
C     Set up a nominal root-to-tip distribution on the first call:

      IF (NCALL == -3) THEN

         I = 4 * NHALF / 3

C        Unnormalized arc lengths between fin sections for this I:

         CALL CHORDSRF (MXIWING, MXKWING, I, 1, NREGSEC, XRG, YRG, ZRG,
     >                  .FALSE., TTOTAL, TRG)

C        Set up a denormalized nominal base distribution using this I.
C        A non-intersection at the low end is anticipated in SURFER's set-up.

         TBASE(1) = TRG(KWING1(IW)) + TINTR(I,IW) *
     >             (TRG(KWING2(IW)) - TRG(KWING1(IW))) ! Left intersection T

         IF (IN == 0) THEN ! Plain right tip
            TBASE(KLWING) = TTOTAL
         ELSE
            TBASE(KLWING) = TRG(KWING1(IN)) + TINTR(I,IN) *
     >             (TRG(KWING2(IN)) - TRG(KWING1(IN))) ! Right intersection T
         END IF

         TTOTAL = TBASE(KLWING) - TBASE(1)

         DO K = 2, KLWING - 1
            TBASE(K) = TBASE(1) + TTOTAL * TNORM(K)
         END DO

         IF (NCRANK > 0) THEN

C           Blend crank T(s) into the base T distribution:

            DO K = 1, NCRANK
               TCRANK(K) = TRG(KCRANK(K))
            END DO

            CALL SMOOTHX (TBASE, KLWING, TCRANK, NCRANK, 3, WORK, IER)

C           Match base grid Ks with geometry Ks at cranks:

            K1 = 1
            DO J = 1, NCRANK
               DO K = K1, KLWING
                  IF (TCRANK(J) == TBASE(K)) THEN
                     MCRANK(J) = K
                     K1 = K1 + 1
                     EXIT
                  END IF
               END DO
            END DO

         END IF

C        Appending the tip as a final crank helps paneling between cranks
C        even if the input NCRANK is 0.

         KCRANK(NCRANK+1) = NREGSEC
         MCRANK(NCRANK+1) = KLWING

      END IF

C     Arc-length-based spanwise gridding (left half only):

      IX = 1

      DO I = NHALF, ILWING

C        Root-to-tip arcs for this I:

         CALL CHORDSRF (MXIWING, MXKWING, I, 1, NREGSEC, XRG, YRG, ZRG,
     >                  .FALSE., TTOTAL, TRG)

C        Transform the base distribution by panels:

         T1 = TRG(KWING1(IW)) + TINTR(I,IW) *
     >       (TRG(KWING2(IW)) - TRG(KWING1(IW)))
         TEVAL(1) = T1
         K1 = 1

         IF (IN == 0) THEN
            TTIP = TTOTAL
         ELSE
            TTIP = TRG(KWING1(IN)) + TINTR(I,IN) *
     >            (TRG(KWING2(IN)) - TRG(KWING1(IN)))
         END IF

         DO J = 1, NCRANK + 1

            IF (J <= NCRANK) THEN
               T2 = TRG(KCRANK(J))
            ELSE
               T2 = TTIP
            END IF

            K2 = MCRANK(J)
            DT = TBASE(K2) - TBASE(K1)
            TSCALE = (T2 - T1) / DT
            TSHIFT = (TBASE(K2) * T1 - TBASE(K1) * T2) / DT

            DO K = K1 + 1, K2
               TEVAL(K) = TBASE(K) * TSCALE + TSHIFT
            END DO

            K1 = K2
            T1 = T2
         END DO
         TEVAL(KLWING) = T2

C        Piecewise linear interpolation for spanwise grid:

         DO K = 1, NREGSEC
            XYZG(K,1) = XRG(I,K)
            XYZG(K,2) = YRG(I,K)
            XYZG(K,3) = ZRG(I,K)
         END DO

         CALL LCSFIT (NREGSEC, TRG, XYZG(1,1), NEW, LINEAR, KLWING,
     >                TEVAL, XYZM(1,1), WORK)

         CALL LCSFIT (NREGSEC, TRG, XYZG(1,2), NEW, LINEAR, KLWING,
     >                TEVAL, XYZM(1,2), WORK)

         CALL LCSFIT (NREGSEC, TRG, XYZG(1,3), NEW, LINEAR, KLWING,
     >                TEVAL, XYZM(1,3), WORK)

         DO K = 1, KLWING
            XBASE(IX,K) = XYZM(K,1)
            YBASE(IX,K) = XYZM(K,2)
            ZBASE(IX,K) = XYZM(K,3)
         END DO

         IX = IX + 1

      END DO

C     If there's an intersection at either end, ensure exactness:

      IF (TINTR(1,IW) /= ZERO) THEN ! See set-up in SURFER
         IX = 1
         DO I = NHALF, ILWING
            XBASE(IX,1) = XINTR(I,IW)
            YBASE(IX,1) = YINTR(I,IW)
            ZBASE(IX,1) = ZINTR(I,IW)
            IX = IX + 1
         END DO
      END IF

      IF (IN /= 0) THEN
         IX = 1
         DO I = NHALF, ILWING
            XBASE(IX,KLWING) = XINTR(I,IN)
            YBASE(IX,KLWING) = YINTR(I,IN)
            ZBASE(IX,KLWING) = ZINTR(I,IN)
            IX = IX + 1
         END DO
      END IF

      END SUBROUTINE VSURF

C***********************************************************************
C
      SUBROUTINE WDINTR (NCALL, MXIDIV, MXKDIV, MXINTR, ILE, IL,
     >                   XRG, YRG, ZRG, MXIBODY, MXJBODY,
     >                   ISTN1, ISTN2, KSTN1, KSTN2, KWING1, KWING2,
     >                   XF, YF, ZF, XTRAP, IDEGREE, ITL, ITU, NINTR,
     >                   XINTR, YINTR, ZINTR, UINTR, VINTR, TINTR,
     >                   IINTR, JINTR, KINTR, IWRIT, FAIL)
C
C     WDINTR calculates a wing/diverter intersection given a regular
C     geometry mesh for each.  Here, "body" and J refer to the wing and K,
C     while the wing in WFINTR is now the diverter.
C
C     12/01/97  DAS  Adaptation of WFINTR.
C     01/14/98   "   NCALL = -4 means (re)do just inboard half of intersection.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   NCALL,                  ! -3 = first call for this case;
     >                           ! -4 = a MATCHDIV call for 1:ILE only
     >   MXIDIV, MXKDIV,         ! Dimensions for diverter sections
     >   MXINTR,                 ! Dimension of intersection arrays
     >   ILE, IL                 ! Diverter leading edge & last indices
      REAL, INTENT (IN), DIMENSION (MXIDIV,MXKDIV) ::
     >   XRG, YRG, ZRG           ! Regularized diverter sections
      INTEGER, INTENT (IN) ::
     >   MXIBODY, MXJBODY,       ! I & J dimens. for wing sections X/Y/ZF
     >   ISTN1, ISTN2,           ! Wing section I range for intersection
     >   KSTN1, KSTN2,           ! ............ K ......................
     >   KWING1, KWING2          ! Diverter K range for intersection
      REAL, INTENT (IN), DIMENSION (MXIBODY,MXJBODY) ::
     >   XF, YF, ZF              ! Regularized wing sections
      REAL, INTENT (IN) ::
     >   XTRAP                   ! X distance beyond lower trailing edge 
                                 ! at which wing pts. are extrapolated to
                                 ! protect the intersection calculation
      INTEGER, INTENT (IN) ::
     >   IDEGREE                 ! 3 = bicubics (INTSEC4) for wing;
                                 ! 1 = bilinear (INTSEC5) ........
      INTEGER, INTENT (OUT) ::
     >   ITL, ITU,               ! Indices of the first/last intersection
     >                           ! points actually on the lower wing
     >   NINTR                   ! Number of intersection points (IL)
      REAL, INTENT (OUT), DIMENSION (MXINTR) ::
     >   XINTR, YINTR, ZINTR,    ! Intersection coordinates, and corresp.
     >   UINTR, VINTR,           ! wing "u"s & "v"s (strmwise & spanwise)
     >   TINTR                   ! and "t"s in elements ITL:ITU
      INTEGER, INTENT (OUT), DIMENSION (MXINTR) ::
     >   IINTR, JINTR, KINTR     ! Corresp. wing & diverter cell indices
      INTEGER, INTENT (IN) ::
     >   IWRIT                   ! For diagnostics if > 0
      LOGICAL, INTENT (OUT) ::
     >   FAIL                    ! Success flag

C     Local constants:

      REAL, PARAMETER :: ONE = 1., ZERO = 0.

C     Local variables:

      INTEGER
     >   I, I1, I2, IDIM, IFAIL, ILAST, INC, ITE, K, KLEFT
      REAL
     >   CD, CW, DX, R, SD, SW
      REAL, DIMENSION (ISTN2+1,KSTN2) ::
     >   XBODY, YBODY, ZBODY, UBODY, VBODY
      LOGICAL
     >   FIRST

C     Execution:

      IDIM = ISTN2 + 1

C     To protect the wing/diverter intersection, we extend the lower wing
C     trailing edge.  Since there is no room in X/Y/ZF, make local copies:

      IF (XTRAP > ZERO) THEN

         DO K = KSTN1, KSTN2
            R = ONE + XTRAP / (XF(1,K) - XF(2,K))
            XBODY(1,K) = R * XF(1,K) + (ONE - R) * XF(2,K)
            YBODY(1,K) = R * YF(1,K) + (ONE - R) * YF(2,K)
            ZBODY(1,K) = R * ZF(1,K) + (ONE - R) * ZF(2,K)

            I = IDIM
            XBODY(2:I,K) = XF(1:ISTN2,K) ! All really wing geometry
            YBODY(2:I,K) = YF(1:ISTN2,K)
            ZBODY(2:I,K) = ZF(1:ISTN2,K)
         END DO

         ILAST = IDIM

      ELSE ! Wing geometry

         ILAST = ISTN2
         DO K = KSTN1, KSTN2
            DO I = ISTN1, ILAST
               XBODY(I,K) = XF(I,K)
               YBODY(I,K) = YF(I,K)
               ZBODY(I,K) = ZF(I,K)
            END DO
         END DO

      END IF

      FIRST = NCALL == -3

      IF (FIRST) THEN

         I1 = ILE ! for diverter
         K  = (KSTN1 + KSTN2) / 2
         UINTR(I1) = (XF(ISTN1,K) - XRG(ILE,1)) / 
     >               (XF(ISTN1,K) - XF(ISTN2,K))
         IINTR(I1) = NINT (REAL (ISTN2 - ISTN1) * UINTR(I1))
         JINTR(I1) = K
         KINTR(I1) = 1
         VINTR(I1) = 0.5
         TINTR(I1) = 0.1
         IF (YRG(ILE,1) < YRG(ILE,2)) TINTR(I1) = 0.9

      END IF

      FAIL = .FALSE. ! Allows retries in WBINTR
      UBODY(ISTN1,KSTN1) = -ONE ! Tells INTSEC4/5 to parameterize the wing
      I1 = ILE ! Do the diverter intersection lower half (in index) first
      I2 = 1

      DO ! Two-pass loop

         CALL WBINTR (MXIDIV, MXKDIV, I1, I2, KWING1, KWING2,
     >                XRG, YRG, ZRG,
     >                IDIM, KSTN2, ISTN1, ILAST, KSTN1, KSTN2,
     >                XBODY, YBODY, ZBODY, 1, IDEGREE, UBODY, VBODY,
     >                XINTR, YINTR, ZINTR, UINTR, VINTR, TINTR,
     >                IINTR, JINTR, KINTR, FIRST, IWRIT, IFAIL, FAIL)

         IF (FAIL) GO TO 999

         IF (NCALL == -4) GO TO 999 ! No need to determine ITL, ITU

         IF (I1 > I2) THEN ! Do the 2nd half diverter from the leading edge
            ITL = 1
            I2  = IL
         ELSE
            EXIT
         END IF

      END DO ! Two-pass loop

      ITU = IL

C     Fudge one or two intersection pts. to be at the lower wing trailing edge?

      IF (XTRAP > ZERO) THEN

         I1  = 1 ! Lower-index end of diverter intersection
         I2  = ILE
         INC = 1

         DO ! Two-pass loop

            IF (IINTR(I1) == 1) THEN

               DO I = I1, I2, INC
                  IF (IINTR(I) > 1) EXIT
               END DO

               KLEFT = JINTR(I)
               DX = XBODY(2,KLEFT+1) - XBODY(2,KLEFT)

               IF (DX == ZERO) THEN
                  R = (XBODY(2,KLEFT) - XINTR(I)) /
     >                  (XINTR(I-INC) - XINTR(I))
               ELSE ! Intersection of Wing/Div. straight lines, Z = S X + C
                  SW = (ZBODY(2,KLEFT+1) - ZBODY(2,KLEFT)) / DX
                  CW =  ZBODY(2,KLEFT)   - XBODY(2,KLEFT)  * SW
                  DX =  XINTR(I-INC) - XINTR(I)
                  SD = (ZINTR(I-INC) - ZINTR(I)) / DX
                  CD =  ZINTR(I)     - XINTR(I)  * SD
                  R  = (CD - CW) / (SW - SD) ! X at intersection
                  R  = (R - XINTR(I)) / DX
               END IF

               ITE = I
               IF (R > 0.5) ITE = I - INC ! Fudge the nearer pt. of the two

               XINTR(ITE) = (ONE - R) * XINTR(I) + R * XINTR(I-INC)
               YINTR(ITE) = (ONE - R) * YINTR(I) + R * YINTR(I-INC)
               ZINTR(ITE) = (ONE - R) * ZINTR(I) + R * ZINTR(I-INC)
               TINTR(ITE) = (ONE - R) * TINTR(I) + R * TINTR(I-INC)

               IF (I1 == 1) THEN ! Adjust other end of the intersection 
                  ITL = ITE
                  I1  = IL
                  INC = -1
               ELSE
                  ITU = ITE
                  EXIT
               END IF

            END IF

         END DO ! Two-pass loop

      END IF

      NINTR = IL

  999 RETURN

      END SUBROUTINE WDINTR

C*******************************************************************************
C
      SUBROUTINE WDSURF (IW, NDIV, IDIV, NWSURF, KPATCH, MXPATCH,
     >                   MXINTR, NCRANK, MCRANK, ILWING, KLWING,
     >                   NNFUSE, EPS, IPDIVLE, IPDIVTE, LENXB,
     >                   XBASE, YBASE, ZBASE, LASTXB, MOSTXB, IXB,
     >                   MXSUBI, MXSUBJ, NSUBI, NSUBJ, ISUB, JSUB, 
     >                   NINTR, XINTR, YINTR, ZINTR, IINTR, JINTR,
     >                   NHALF, ITL, ITU, NBSPLIT, IBSPLIT,
     >                   IMXCOM, JMXCOM, LWING, IWRIT, FAIL)
C
C     WDSURF (re)panels a wing surface in the presence of one or more under-wing
C     diverters or pylons.  The initial paneling of the plain wing is input as
C     KPATCH; this already captures fixed spanwise features defined as crank
C     stations, and serves for much of the wing unaffected by the diverter(s).
C
C     Under-wing diverter/pylon convention:
C     -------------------------------------
C
C     INBOARD:   I = 1 : leading edge
C     OUTBOARD:  I = leading edge : NINTR
C
C     12/18/97- DAS  Initial implementation (but some interpolations are non-
C     01/02/98   "   parametric).
C     01/12/98   "   Handled diverters going off the trailing edge (MATCHDIV).
C     01/20/98   "   Handled patches aft of pylons.
C     01/27/98   "   All surface interpolations are parametric now.
C     02/27/98   "   Allowed for yawed pylons/diverters.  The associated crank
C                    stations in the plain wing grid are nominal only.
C     10/30/98   "   Upper wing of diverter cases now reflect lower wing panels
C                    at the trailing edge (forced by perturbed grid problems).
C     03/06/99   "   Packed patch version.
C     05/01/99   "   Subpatch version started; finished 05/12/99 (awkward!).
C
C*******************************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   IW,                   ! Wing surface #
     >   NDIV, IDIV(NDIV),     ! # diverters and their surface #s
     >   NWSURF                ! # wing-type surfaces

      INTEGER, INTENT (INOUT) ::
     >   KPATCH                ! Input as patch containing plain wrapped wing;
                               ! output as next available patch
      INTEGER, INTENT (IN) ::
     >   MXPATCH,
     >   MXINTR,
     >   NCRANK,               ! These are all wing quantities (not diverter)
     >   MCRANK(NCRANK+1),
     >   ILWING, KLWING,
     >   NNFUSE                ! 0 means no fuselage

      REAL, INTENT (IN) ::     ! Search tolerance
     >   EPS

      INTEGER, INTENT (INOUT), DIMENSION (NDIV) ::
     >   IPDIVLE, IPDIVTE      ! Input with diverter LE/TE intersection index
                               ! estimates; adjusted in case RECTIFY was used
      INTEGER, INTENT (IN) ::
     >   LENXB                 ! Available space in X/Y/ZBASE(*)

      REAL, DIMENSION (LENXB) ::
     >   XBASE, YBASE, ZBASE   ! Packed patches

      INTEGER, INTENT (INOUT) ::
     >   LASTXB,               ! Last active element of X/Y/ZBASE(*) ...
     >   MOSTXB,               ! ... and highest (temporarily?) used element
     >   IXB(MXPATCH)          ! IXB(L) = start of patch L

      INTEGER, INTENT (IN) ::  ! Subpatch info.
     >   MXSUBI, MXSUBJ

      INTEGER, INTENT (OUT) ::
     >   NSUBI(MXPATCH),
     >   NSUBJ(MXPATCH),
     >   ISUB(0:MXSUBI,MXPATCH),
     >   JSUB(0:MXSUBJ,MXPATCH)

      INTEGER, INTENT (INOUT) ::
     >   NINTR(NWSURF)         ! Intersection info

      REAL, INTENT (INOUT), DIMENSION (MXINTR,NWSURF) ::
     >   XINTR, YINTR, ZINTR   ! Needed for diverter/wing and wing/fuselage

      INTEGER, INTENT (INOUT), DIMENSION (MXINTR,NWSURF) ::
     >   IINTR, JINTR          ! Needed for wing/fuselage (only)

      INTEGER, INTENT (IN), DIMENSION (NWSURF) ::
     >   NHALF,                ! Saves searching for diverter leading edges
     >   ITL, ITU              ! Wing/div. intersections may go off the planform

      INTEGER, INTENT (OUT) ::
     >   NBSPLIT               ! NBSPLIT = 1: up/lo body splits <-> diverter LE;
                               !         = 2: lower body splits <-> pylon LE/TE;
      INTEGER, INTENT (OUT) ::
     >   IBSPLIT(3)            ! IBSPLIT(1)   = wing index IW
                               ! IBSPLIT(2:3) = implicit split indices for FSURF

      INTEGER, INTENT (INOUT), DIMENSION (MXPATCH) ::
     >   IMXCOM, JMXCOM        ! Patch dimensions

      INTEGER, INTENT (OUT) ::
     >   LWING                 ! For NACSURF2 to match with 1st aft nacelle cap

      INTEGER, INTENT (IN) ::
     >   IWRIT

      LOGICAL, INTENT (OUT) ::
     >   FAIL

C     Local constants:

      INTEGER, PARAMETER ::
     >   KASE  = 3, ! See PLBICUT
     >   LINUV = 2, ! See PLINCUB
     >   NKTRI = 20 ! # pts. on upper triangle patches

      REAL, PARAMETER ::
     >   DERIVS = -999., HALF = 0.5, ONE = 1., ZERO = 0.

      LOGICAL, PARAMETER ::
     >   CLOSED = .FALSE.

      CHARACTER, PARAMETER ::
     >   LINEAR * 1 = 'L', LOOSE * 1 = 'B'

C     Local variables:

      INTEGER, DIMENSION (NDIV) ::
     >   KPDIVLE, KPDIVTE

      INTEGER
     >   KCNUMBER(0:NDIV), KDCRANK(NDIV+1)

      INTEGER
     >   I, I1, I2, IB, IC, ID, ID2, IER, IL, ILE, ILED, ILIMIT, IM, IP,
     >   IS, ISIZE, ISMAT, ISTART, IX, IX1, J, JD, K, K1, K2, KBREAK,
     >   KC, KD, KD2, KLEFT, KP, KSMAT, KT, L, L1, L2, LCAP, LUNERR,
     >   M, MI, MI2, MK, NADD, NAFT, NI, NK, NMONOX, NUMI, NUMK

      REAL, DIMENSION (NDIV) ::
     >   UDIVLE, VDIVLE

      REAL, DIMENSION (KLWING) ::
     >   ULINE, VLINE, ULINE0, VLINE0, XLINE, YLINE, ZLINE

      REAL, DIMENSION (ILWING,KLWING) ::
     >   UEVAL, VEVAL, UPLAIN, VPLAIN, XPLAIN, YPLAIN, ZPLAIN

      REAL
     >   SMATRIX(4,4,3), TFIWORK(ILWING+2*KLWING), UVRANGE(4),
     >   XYZCUT(3),
     >   DRANGE, DV, DZ, DZMIN, P, PM1, Q, R, RU, RV, RVRANGE,
     >   UCUT, VCUT, UE, VE, V1, V2, VSCALE, VSHIFT, XT, YT, ZT,
     >   XTE, YTE, ZLE

      LOGICAL
     >   PYLON

C     Storage:

      SAVE UVRANGE
      DATA UVRANGE /ZERO, ONE, ZERO, ONE/ ! See PLBICUT

C     Execution:
C     ----------

C     Copy the plain-wing paneling so its patch location can be reused:

      L1 = KPATCH     ! Remains the first wing patch
      IX = IXB(L1)    ! Start of input wing patch
      LASTXB = IX - 1 ! After releasing the plain wing patch

      DO K = 1, KLWING
         DO I = 1, ILWING
            XPLAIN(I,K) = XBASE(IX)
            YPLAIN(I,K) = YBASE(IX)
            ZPLAIN(I,K) = ZBASE(IX)
            IX = IX + 1
         END DO
      END DO

      UPLAIN(1,1) = ZERO ! -999. would suppress normalization

      CALL PARAM2D (ILWING, KLWING, 1, ILWING, 1, KLWING,
     >              XPLAIN, YPLAIN, ZPLAIN, UPLAIN, VPLAIN)

C     There is supposed to be a (nominal) crank station, and hence a plain-wing
C     mesh K, associated with each diverter/pylon. Locate it for each diverter:

      DO J = 1, NDIV
         ID    = IDIV(J)
         ILED  = NHALF(ID)
         ZLE   = ZINTR(ILED,ID)
         DZMIN = 1.E+10

         DO KC = 1, NCRANK
            KD = MCRANK(KC)
            DZ = ABS (ZPLAIN(1,KD) - ZLE)
            IF (DZ < DZMIN) THEN
               DZMIN = DZ
               K = KC
            END IF
         END DO

         KCNUMBER(J) = K
         KDCRANK(J)  = MCRANK(K)
      END DO

      KCNUMBER(0)     = 0 ! Helps the following loop
      KDCRANK(NDIV+1) = MCRANK(NCRANK+1) ! At the tip; helps the following loop

C     Cranks at pylon/diverter(s) define patch boundaries.
C     Establish subpatch K boundaries at all cranks:
C     ----------------------------------------------

      J  = 1
      K1 = 1

      DO L = L1, L1 + NDIV
         I = 0

         DO KC = KCNUMBER(J-1) + 1, NCRANK + 1
            IF (MCRANK(KC) <= KDCRANK(J)) THEN
               I = I + 1
               JSUB(I,L) = MCRANK(KC) - K1 + 1
            ELSE
               K1 = KDCRANK(J) ! Mesh K of previous diverter
               EXIT
            END IF
         END DO

         NSUBJ(L)  = I
         JMXCOM(L) = JSUB(I,L)
         J = J + 1 ! Next diverter, or tip
      END DO ! Next patch

      IBSPLIT(1) = IW           ! Wing surface <-> body patch splits @ div. TE
      ILE    = NHALF(IW)
      NMONOX = ILE - 6          ! To protect some uses of INTERVAL on X
      ILIMIT = ILE + 1          ! Assuming RECTIFY is employed appropriately
                                ! + 1 preserves 4-pt. formula at the wing LE
C*****ILIMIT = (3 * ILWING) / 5 ! Limits the searching for most of the
                                ! under-wing paneling; above-wing pylons
                                ! will need another version of WDSURF

      XTE    = (XINTR(1,ID) + XINTR(NINTR(ID),ID)) * HALF
      DRANGE = XTE - XPLAIN(ILE,1) ! ~ Match (u,v) precision in PLBICUT
      DRANGE = MAX (ABS (XTE), ABS (XPLAIN(ILE,1)), DRANGE)
      ISMAT  = 0
      KSMAT  = 0
      LUNERR = -ABS (IWRIT)

C     -----------------------------------------------------------------
C     Awkward consequence of subpatching: we have to know the full size
C     of each (larger) patch before we can store any of its subpatches.
C     Make sure of the wing cell indices for each pylon/diverter TE/LE,
C     because the values saved after WSURF may be affected by RECTIFY.
C     -----------------------------------------------------------------

      PYLON = .FALSE.
      DO J = 1, NDIV

         KD  = KDCRANK(J)
         ID  = IDIV(J)
         XTE = (XINTR(1,ID) + XINTR(NINTR(ID),ID)) * HALF

         IF (XTE < XPLAIN(1,KD)) THEN ! Pylon, not diverter

            PYLON = .TRUE. ! For all J is assumed
            XYZCUT(1) = XTE
            XYZCUT(3) = (ZINTR(1,ID) + ZINTR(NINTR(ID),ID)) * HALF
            IP = IPDIVTE(J) ! Now saved after WSURF call

            CALL INTERVAL (NMONOX, XPLAIN(1,KD), XTE, -ONE, IP)

            KP = KD
            UCUT = UPLAIN(IP,KP)
            VCUT = VPLAIN(IP,KP)

            CALL PLBICUT (KASE, ILWING, KLWING, 1, ILE, 1, KLWING,
     >                    XPLAIN, YPLAIN, ZPLAIN, UPLAIN, VPLAIN,
     >                    IP, KP, ISMAT, KSMAT, SMATRIX, EPS, DRANGE,
     >                    XYZCUT, UVRANGE, UCUT, VCUT, P, Q, LUNERR,
     >                    IER)

            IF (IER > 3) THEN
               JD = 1
               I  = 1
               GO TO 800
            END IF

            IPDIVTE(J) = IP
            KPDIVTE(J) = KP

         ELSE ! Diverter case

            IPDIVTE(J) = 1
            KPDIVTE(J) = KD - 1 ! Still not precise

         END IF

C        Likewise for the pylon/diverter LE:

         ILED = NHALF(ID)
         XYZCUT(1) = XINTR(ILED,ID)
         XYZCUT(3) = ZINTR(ILED,ID)
         IP = IPDIVLE(J)
         KP = KD
         UCUT = UPLAIN(IP,KP)
         VCUT = VPLAIN(IP,KP)

         CALL PLBICUT (KASE, ILWING, KLWING, 1, ILE, 1, KLWING,
     >                 XPLAIN, YPLAIN, ZPLAIN, UPLAIN, VPLAIN,
     >                 IP, KP, ISMAT, KSMAT, SMATRIX, EPS, DRANGE,
     >                 XYZCUT, UVRANGE, UCUT, VCUT, P, Q, LUNERR, IER)

         IF (IER > 3) THEN
            JD = 1
            GO TO 800
         END IF
            
         IPDIVLE(J) = IP
         KPDIVLE(J) = KP

      END DO ! Next pylon/diverter

C     Establish subpatch I boundaries:
C     --------------------------------

      ID   = IDIV(1)
      ILED = NHALF(ID)

C     Lower surface patches:

      IF (PYLON) THEN

         NI = IPDIVTE(1) + 1 ! For all aft subpatches; assume constant NHALF(ID)

         DO L = L1, L1 + NDIV
            NSUBI(L)  = 3
            ISUB(1,L) = NI
            ISUB(2,L) = NI + ILED - 1
            ISIZE     = ISUB(2,L) + ILE - IPDIVLE(1) ! Full patch size
            ISUB(3,L) = ISIZE
            IMXCOM(L) = ISIZE
            IXB(L)    = LASTXB + 1
            LASTXB    = LASTXB + ISIZE * JMXCOM(L)
         END DO

         NBSPLIT = 2 ! # body patch splits below wing (implicit)
         IBSPLIT(2) = ISIZE - ISUB(2,L1) + 1 ! Counting from LE for FSURF
         IBSPLIT(3) = ISIZE - NI + 1

      ELSE ! Diverter case

         L = L1 ! Beside the body
         NSUBI(L)  = 2
         NI        = ILED - ITL(ID) + 1
         ISUB(1,L) = NI
         ISIZE     = NI + ILE - IPDIVLE(1)
         ISUB(2,L) = ISIZE
         IMXCOM(L) = ISIZE
         IXB(L)    = LASTXB + 1
         LASTXB    = LASTXB + ISIZE * JMXCOM(L)

         NBSPLIT = 1 ! # body patch splits above/below wing (implicit)
         IBSPLIT(2) = ISIZE - NI + 1 !

         J = 1
         DO L = L1 + 1, L1 + NDIV
            NUMI = ITU(ID) - ILED + 1 ! # pts. along outboard edge of div. J
            J  = J + 1
            IF  (J <= NDIV) THEN ! For next L
               ID   = IDIV(J)
               ILED = NHALF(ID)
            END IF
            NSUBI(L)  = 2
            ISUB(1,L) = NUMI
            ISIZE     = IMXCOM(L-1) + (NUMI - NI)
            NI        = NUMI
            ISUB(2,L) = ISIZE
            IMXCOM(L) = ISIZE
            IXB(L)    = LASTXB + 1
            LASTXB    = LASTXB + ISIZE * JMXCOM(L)
         END DO

      END IF

C     Upper surface patch boundaries:
C     -------------------------------

      IL = ILWING - MAX (10, ILWING / 20) ! Index for upper triangle apexes
      NI = IL - ILE + 1
      M  = L1
      L2 = L1 + 2 * NDIV + 1

      DO L = L1 + NDIV + 1, L2

         IF (PYLON) THEN
            ISUB(1,L)  = ILE
         ELSE
            NSUBI(L)   = 2
            ISUB(1,L)  = NI
            ISUB(2,L)  = ILE
            IBSPLIT(3) = NI
         END IF
         IMXCOM(L) = ILE

         J = NSUBJ(M) ! Match the lower surface spanwise subpatching
         NSUBJ(L)    = J
         JSUB(1:J,L) = JSUB(1:J,M)
         JMXCOM(L)   = JMXCOM(M)
         IXB(L)      = LASTXB + 1
         LASTXB      = LASTXB + ILE * JMXCOM(L)

         M = M + 1
      END DO

      L = L2

C     Triangular patches above diverter cutouts?

      IF (.NOT. PYLON) THEN
         NI = ILWING - IL + 1
         DO J = 1, NDIV
            L = L + 1
            IMXCOM(L) = NI
            ISUB(1,L) = NI
            JMXCOM(L) = NKTRI
            JSUB(1,L) = NKTRI
            IXB(L)    = LASTXB + 1
            LASTXB    = LASTXB + NI * NKTRI
         END DO
      END IF

C     Exceeded packed patch space?  Checker assumes patch not counted yet.

      I = LASTXB - IMXCOM(L) * JMXCOM(L)

      CALL XBCHECK (LENXB, MXPATCH, L, I, IMXCOM(L), JMXCOM(L), IW,
     >              IWRIT, MOSTXB)

C     Prepare for subpatch 1 of lower wing patch 1:

      L      = L1
      ISIZE  = IMXCOM(L)
      IX     = IXB(L)
      ID     = IDIV(1)
      ILED   = NHALF(ID)
      KD     = KDCRANK(1)
      I      = NSUBJ(L)
      KBREAK = JSUB(I-1,L) ! K at crank nearest to/inboard of pylon # 1, or 1

C     Patches aft of pylon(s)?
C     ------------------------

      IF (PYLON) THEN  

         XYZCUT(1) = (XINTR(1,ID) + XINTR(NINTR(ID),ID)) * HALF
         XYZCUT(3) = (ZINTR(1,ID) + ZINTR(NINTR(ID),ID)) * HALF
         IP = IPDIVTE(1)
         KP = KPDIVTE(1)
         UCUT = UPLAIN(IP,KP)
         VCUT = VPLAIN(IP,KP)

         CALL PLBICUT (KASE, ILWING, KLWING, 1, ILE, 1, KLWING,
     >                 XPLAIN, YPLAIN, ZPLAIN, UPLAIN, VPLAIN,
     >                 IP, KP, ISMAT, KSMAT, SMATRIX, EPS, DRANGE,
     >                 XYZCUT, UVRANGE, UCUT, VCUT, P, Q, LUNERR, IER)

         PM1 = ONE - P

C        Subpatch next to the fuselage, inboard and aft of pylon # 1:
C        ------------------------------------------------------------

         NI = IP + 1 ! For all aftmost subpatches

         DO K = 1, KD ! Wing trailing edge
            XBASE(IX) = XPLAIN(1,K)
            YBASE(IX) = YPLAIN(1,K)
            ZBASE(IX) = ZPLAIN(1,K)
            IX = IX + ISIZE
         END DO

C        Inboard of the break is straightforward:

         IX1 = IXB(L)
         DO K = 1, KBREAK
            IX  = IX1
            IX1 = IX1 + ISIZE
            DO I = 2, IP ! Skip (1,K) and NI,K)
               IX = IX + 1
               XBASE(IX) = XPLAIN(I,K)
               YBASE(IX) = YPLAIN(I,K)
               ZBASE(IX) = ZPLAIN(I,K)
            END DO
         END DO

         IX = IXB(L) + IP ! (NI,1)
         I  = IP

         DO K = 1, KBREAK ! Spanwise line corresp. to pylon TE

            UE = PM1 * UPLAIN(IP,K) + P * UPLAIN(IP+1,K)
            VE = PM1 * VPLAIN(IP,K) + P * VPLAIN(IP+1,K)
            IP = I
            KP = K

            CALL PLINCUB (ILWING, KLWING, 1, ILIMIT, 1, KLWING,
     >                    XPLAIN, YPLAIN, ZPLAIN, UPLAIN, VPLAIN,
     >                    UE, VE, LINUV, LOOSE, CLOSED, IP, KP, EPS,
     >                    XBASE(IX), YBASE(IX), ZBASE(IX), IER)
            IX = IX + ISIZE
         END DO

C        Use TFI for subpatch (1:NI, KBREAK:KD):

         IP = I
         UEVAL(1:IP,KBREAK) = UPLAIN(1:IP,KBREAK) ! K = KBREAK edge
         VEVAL(1:IP,KBREAK) = VPLAIN(1:IP,KBREAK)

         UEVAL(NI,KBREAK)   = UE ! Subpatch corner from above loop
         VEVAL(NI,KBREAK)   = VE

         UEVAL(1,KBREAK:KD) = UPLAIN(1,KBREAK:KD) ! Wing TE
         VEVAL(1,KBREAK:KD) = VPLAIN(1,KBREAK:KD)

C        Spanwise edge through pylon TE (internal procedure below):

         RVRANGE = ONE / (VPLAIN(1,KD) - VPLAIN(1,KBREAK))

         CALL UV_LINE_INBOARD (KBREAK, KD, IPDIVTE(1), KPDIVTE(1), NI)

C        Warp the nominal edge along K = KD to a line exactly through the
C        pylon trailing edge, in (u,v) space, according to end-pt. motion:

         UEVAL(1,KD) = UPLAIN(1,KD) ! UCUT, VCUT end is done by UV_LINE_INBOARD
         VEVAL(1,KD) = VPLAIN(1,KD)

         UE = UPLAIN(NI,KD) ! Temporarily replace point NI in the plain mesh
         VE = VPLAIN(NI,KD)
         UPLAIN(NI,KD) = PM1 * UPLAIN(IP,KD) + P * UPLAIN(NI,KD)
         VPLAIN(NI,KD) = PM1 * VPLAIN(IP,KD) + P * VPLAIN(NI,KD)

         CALL NULINE2D (1, NI, UPLAIN(1,KD), VPLAIN(1,KD),
     >                  UEVAL(1,KD), VEVAL(1,KD))

         UPLAIN(NI,KD) = UE ! Recover original values
         VPLAIN(NI,KD) = VE

C        Interior (u,v)s:

         CALL TFI2D (ILWING, 1, NI, KBREAK, KD, UEVAL, VEVAL, TFIWORK)

C        Interpolate interior (x,y,z)s & some edges:

C        Address (B(i,k)) = Address (B(1,1)) + (k - 1)idim + i - 1

         IX1 = IXB(L) + KBREAK * ISIZE + 1 ! (2,KBREAK+1)
         KP  = KBREAK + 1

         DO K = KBREAK + 1, KD
            IX  = IX1
            IX1 = IX1 + ISIZE
            IP  = 2
            DO I = 2, NI

               CALL PLINCUB (ILWING, KLWING, 1, ILIMIT, 1, KLWING,
     >                       XPLAIN, YPLAIN, ZPLAIN, UPLAIN, VPLAIN,
     >                       UEVAL(I,K), VEVAL(I,K), LINUV, LOOSE,
     >                       CLOSED, IP, KP, EPS, 
     >                       XBASE(IX), YBASE(IX), ZBASE(IX), IER)
               IX = IX + 1
            END DO
         END DO

         IX = IX - 1 ! (NI,KD) again
         XBASE(IX) = XYZCUT(1)
         YBASE(IX) = XYZCUT(2)
         ZBASE(IX) = XYZCUT(3)

         MK = KD ! # Ks for J - 1 patch


C        Patch(es) aft of and between pylons:
C        ------------------------------------

         DO J = 2, NDIV

            ID2 = IDIV(J)
            KD2 = KDCRANK(J)
            NK  = KD2 - KD + 1
            L   = L + 1

            CALL XBCOPY2 (LENXB, XBASE, YBASE, ZBASE, IXB, IMXCOM,
     >                    L - 1, 1, NI, MK, MK, L, 1, 1)

            DO I = 1, NI
               UEVAL(I,1) = UEVAL(I,MK) ! Common chordwise edge
               VEVAL(I,1) = VEVAL(I,MK)
            END DO

C           Locate the current pylon trailing edge in the plain wing mesh:

            XYZCUT(1) = (XINTR(1,ID2) + XINTR(NINTR(ID2),ID2)) * HALF
            XYZCUT(3) = (ZINTR(1,ID2) + ZINTR(NINTR(ID2),ID2)) * HALF
            IP = IPDIVTE(J)
            KP = KPDIVTE(J)
            UCUT = UPLAIN(IP,KP)
            VCUT = VPLAIN(IP,KP)

            CALL PLBICUT (KASE, ILWING, KLWING, 1, ILE, 1, KLWING,
     >                    XPLAIN, YPLAIN, ZPLAIN, UPLAIN, VPLAIN,
     >                    IP, KP, ISMAT, KSMAT, SMATRIX, EPS, DRANGE,
     >                    XYZCUT, UVRANGE, UCUT, VCUT, P, Q, LUNERR,
     >                    IER)

            PM1 = ONE - P

C           Warp the nominal points aft of the current pylon:

            UEVAL(1,NK)  = UPLAIN(1,KD2)
            VEVAL(1,NK)  = VPLAIN(1,KD2)
            UEVAL(NI,NK) = UCUT
            VEVAL(NI,NK) = VCUT

            UE = UPLAIN(NI,KD2) ! Temporarily replace ...
            VE = VPLAIN(NI,KD2)
            UPLAIN(NI,KD2) = PM1 * UPLAIN(IP,KD2) + P * UPLAIN(NI,KD2)
            VPLAIN(NI,KD2) = PM1 * VPLAIN(IP,KD2) + P * VPLAIN(NI,KD2)

            CALL NULINE2D (1, NI, UPLAIN(1,KD2), VPLAIN(1,KD2),
     >                     UEVAL(1,NK), VEVAL(1,NK))

            UPLAIN(NI,KD2) = UE ! Recover ...
            VPLAIN(NI,KD2) = VE

C           Spanwise edge (u,v)s:

            RVRANGE = ONE / (VPLAIN(1,KD2) - VPLAIN(1,KD))
            KP = KD
            IX = IXB(L)

            DO K = 2, NK
               IX = IX + ISIZE ! (1,K,L)
               KP = KP + 1
               UEVAL(1,K) = ZERO
               VEVAL(1,K) = VPLAIN(1,KP)
               XBASE(IX)  = XPLAIN(1,KP)
               YBASE(IX)  = YPLAIN(1,KP)
               ZBASE(IX)  = ZPLAIN(1,KP)
               R  = (VPLAIN(1,KP) - VPLAIN(1,KD)) * RVRANGE
               UEVAL(NI,K) = (ONE - R) * UEVAL(NI,1) + R * UEVAL(NI,NK)
               VEVAL(NI,K) = (ONE - R) * VEVAL(NI,1) + R * VEVAL(NI,NK)
            END DO

            UEVAL(NI,NK) = UCUT ! Again, precisely
            VEVAL(NI,NK) = VCUT

C           Interior (u,v)s:

            CALL TFI2D (ILWING, 1, NI, 1, NK, UEVAL, VEVAL, TFIWORK)

C           Interpolate interior (x,y,z)s & some edges:

            IX1 = IXB(L) + 1 ! (2,1)
            KP  = KD + 1
            DO K = 2, NK
               IX1 = IX1 + ISIZE ! (2,K)
               IX  = IX1
               IP  = 2
               DO I = 2, NI

                  CALL PLINCUB (ILWING, KLWING, 1, ILIMIT, 1, KLWING,
     >                          XPLAIN, YPLAIN, ZPLAIN, UPLAIN, VPLAIN,
     >                          UEVAL(I,K), VEVAL(I,K), LINUV, LOOSE,
     >                          CLOSED, IP, KP, EPS,
     >                          XBASE(IX), YBASE(IX), ZBASE(IX), IER)
                  IX = IX + 1
               END DO
            END DO

            IX = IX - 1 ! (NI,NK,L) again
            XBASE(IX) = XYZCUT(1)
            YBASE(IX) = XYZCUT(2)
            ZBASE(IX) = XYZCUT(3)

            KD = KD2
            MK = NK

         END DO ! Next pylon


C        Patch aft and outboard of the outermost pylon:
C        ----------------------------------------------

         L  = L + 1
         NK = KLWING - KD + 1
         K1 = KD
         KC = KCNUMBER(NDIV) + 1
         K2 = MCRANK(KC) ! Possibly the tip

C        Common chordwise edge:

         CALL XBCOPY2 (LENXB, XBASE, YBASE, ZBASE, IXB, IMXCOM,
     >                 L - 1, 1, NI, MK, MK, L, 1, 1)

         DO I = 1, NI
            UEVAL(I,1) = UEVAL(I,MK)
            VEVAL(I,1) = VEVAL(I,MK)
         END DO

         KP = K1 + NK - 1
         KT = 1
         IX = IXB(L)

         DO K = K1 + 1, KP  ! Trailing edge out to the tip
            IX = IX + ISIZE ! (1,2,L) first time
            XBASE(IX) = XPLAIN(1,K)
            YBASE(IX) = YPLAIN(1,K)
            ZBASE(IX) = ZPLAIN(1,K)
            KT = KT + 1
            UEVAL(1,KT) = ZERO
            VEVAL(1,KT) = VPLAIN(1,K)
         END DO

C        Spanwise edge through outboard pylon TE out to next crank:

         RVRANGE = ONE / (VPLAIN(1,K2) - VPLAIN(1,K1))

         CALL UV_LINE_OUTBOARD (K1,K2, IPDIVTE(NDIV), KPDIVTE(NDIV), NI)

C        Impose the relative (u,v) distribution aft of the outermost
C        pylon on lines K = crank to tip in the plain wing mesh.  Do it now
C        because we need the K2 result to do TFI between K1 and K2.

         NUMI = IPDIVTE(NDIV) + 2
         NUMK = K2 - K1 + 1 ! For TFI
         IP = 1
         KP = K1

         DO I = 2, NI

            CALL RIPPLE2D (ILWING, KLWING, 1, NUMI, 1, K2,
     >                     UPLAIN, VPLAIN, UEVAL(I,1), VEVAL(I,1),
     >                     IP, KP, EPS, P, Q, IER)
            PM1 = ONE - P

            J = NUMK
            DO K = K2, KLWING
               UEVAL(I,J) = PM1 * UPLAIN(IP,K) + P * UPLAIN(IP+1,K)
               VEVAL(I,J) = PM1 * VPLAIN(IP,K) + P * VPLAIN(IP+1,K)
               J = J + 1
            END DO

         END DO

C        Interior (u,v)s aft of outer pylon, between pylon and crank (or tip):

         CALL TFI2D (ILWING, 1, NI, 1, NUMK, UEVAL, VEVAL, TFIWORK)

C        Interpolate (x,y,z)s between outer pylon K and tip (aft of pylon):

         J   = 2
         KP  = K1
         IX1 = IXB(L)

         DO K = K1 + 1, KLWING
            IX1 = IX1 + ISIZE ! (1,2,L) first time
            IX  = IX1
            IP  = 1
            DO I = 2, NI
               IX = IX + 1

               CALL PLINCUB (ILWING, KLWING, 1, ILE, 1, KLWING,
     >                       XPLAIN, YPLAIN, ZPLAIN, UPLAIN, VPLAIN,
     >                       UEVAL(I,J), VEVAL(I,J), LINUV, LOOSE,
     >                       CLOSED, IP, KP, EPS,
     >                       XBASE(IX), YBASE(IX), ZBASE(IX), IER)
            END DO
            J = J + 1
         END DO

         L  = L - NDIV    ! Patch next to body
         KD = KDCRANK(1)  ! Recover inboard pylon cell index
         I1 = NI          ! First wing/fuselage intersection pt. to update
         IB = IXB(L) + NI ! Index in X/Y/ZBASE(*) of point (ISTART,1,L)
         NAFT = NI - 1    ! Versus 0 for diverter case
         NI = ILED        ! Alongside pylon
         ISTART = 2       ! Avoid redoing the subpatch edges

      ELSE ! Aft diverter case

         I1 = 1
         IB = IXB(L)
         NAFT = 0
         NI = ISUB(1,L)
         ISTART = 1

      END IF

      IS = ISTART ! 1 or 2 is also the next I subpatch number


C     Patch between the fuselage and pylon/diverter # 1:
C     --------------------------------------------------

      IP = IPDIVTE(1)
      KP = KPDIVTE(1)
      UCUT = UPLAIN(IP,KP)
      VCUT = VPLAIN(IP,KP)
      RVRANGE = ONE / (VPLAIN(1,KD) - VPLAIN(1,KBREAK)) ! See UV_LINE_INBOARD

      DO I = ITL(ID) + ISTART - 1, ILED

C        Back out (u,v) for this inboard intersection pt. in the plain mesh:

         XYZCUT(1) = XINTR(I,ID)
         XYZCUT(3) = ZINTR(I,ID)

         CALL PLBICUT (KASE, ILWING, KLWING, 1, ILE, 1, KLWING,
     >                 XPLAIN, YPLAIN, ZPLAIN, UPLAIN, VPLAIN,
     >                 IP, KP, ISMAT, KSMAT, SMATRIX, EPS, DRANGE,
     >                 XYZCUT, UVRANGE, UCUT, VCUT, P, Q, LUNERR, IER)

         IF (IER > 3) THEN
            JD = 1
            GO TO 800
         END IF
            
C        Interpolate a corresponding spanwise line in the (u,v) mesh:

         PM1 = ONE - P
         DO K = 1, KBREAK
            ULINE(K)   = PM1 * UPLAIN(IP,K) + P * UPLAIN(IP+1,K)
            UEVAL(1,K) = ULINE(K)
            VLINE(K)   = PM1 * VPLAIN(IP,K) + P * VPLAIN(IP+1,K)
            VEVAL(1,K) = VLINE(K)
         END DO

C        Impose the trailing edge relative v distribution on the panel next
C        to the inboard pylon in a way that handles non-straight (u,v) lines:

         CALL UV_LINE_INBOARD (KBREAK, KD, IP, KP, 1)

C        Evaluate (x,y,z) along the interpolated (u,v) line:

         IX = IB
         KP = 1

         DO K = 1, KD ! Not KD - 1, so we get the leading edge IP and (u,v)

            CALL PLINCUB (ILWING, KLWING, 1, ILE, 1, KLWING,
     >                    XPLAIN, YPLAIN, ZPLAIN, UPLAIN, VPLAIN,
     >                    UEVAL(1,K), VEVAL(1,K), LINUV, LOOSE, CLOSED,
     >                    IP, KP, EPS,
     >                    XBASE(IX), YBASE(IX), ZBASE(IX), IER)
            IX = IX + ISIZE
         END DO

         IX = IX - ISIZE ! (IB,KD,L) again
         IB = IB + 1
         XBASE(IX) = XINTR(I,ID)
         YBASE(IX) = YINTR(I,ID)
         ZBASE(IX) = ZINTR(I,ID)

      END DO ! Next I along inboard diverter edge

      UDIVLE(1) = UCUT ! For the subpatch upstream of this one
      VDIVLE(1) = VCUT


C     Update the corresponding part of the wing/fuselage intersection?
C     ----------------------------------------------------------------

      IF (NNFUSE /= 0) THEN

         I2   = IP ! Last wing/fuselage intersection pt. to replace
         NADD = NI - (I2 - I1 + 1) ! Number of pts. to add/subtract
         NUMI = NINTR(IW) + NADD   ! New intersection count

         IF (NADD > 0) THEN ! Make more room

            IF (NUMI > MXINTR) THEN
               IF (IWRIT > 0) THEN
                  WRITE (IWRIT, '(/, A, 2I5)')
     >               ' WDSURF:  NUMI > MXINTR:', NUMI, MXINTR
                  GO TO 900
               END IF
            END IF

            J = NUMI
            DO I = NINTR(IW), I2 + 1, -1
               XINTR(J,IW) = XINTR(I,IW)
               YINTR(J,IW) = YINTR(I,IW)
               ZINTR(J,IW) = ZINTR(I,IW)
               IINTR(J,IW) = IINTR(I,IW)
               JINTR(J,IW) = JINTR(I,IW)
               J = J - 1
            END DO

         END IF

         I = IXB(L) + NAFT ! Pylon/div. TE pt. location
         J = I1

         DO IX = I, I + NI - 1
            XINTR(J,IW) = XBASE(IX)
            YINTR(J,IW) = YBASE(IX)
            ZINTR(J,IW) = ZBASE(IX)
            IINTR(J,IW) = 0 ! Tells fuselage paneling to recover via BODYSRCH
            J  = J  + 1
         END DO

         IF (NADD < 0) THEN ! Shift down

            DO I = I2 + 1, NINTR(IW)
               XINTR(J,IW) = XINTR(I,IW)
               YINTR(J,IW) = YINTR(I,IW)
               ZINTR(J,IW) = ZINTR(I,IW)
               IINTR(J,IW) = IINTR(I,IW)
               JINTR(J,IW) = JINTR(I,IW)
               J = J + 1
            END DO

         END IF

         NINTR(IW) = NUMI

      END IF

C     Handle patches between diverters:
C     ---------------------------------

      DO J = 2, NDIV

         L   = L + 1
         KC  = KCNUMBER(J)
         KD2 = MCRANK(KC)
         NK  = KD2 - KD + 1
         NI  = ITU(ID) - ILED + 1
         ISIZE = IMXCOM(L)

C        Determine (u,v)s along the outboard half of pylon/diverter J - 1:

         IB = NI                ! Two pointers corresponding to diverter LE
         IX = IXB(L) + NAFT + IB
         IP = IPDIVLE(J-1)
         KP = KPDIVLE(J-1)
         UCUT = UPLAIN(IP,KP)
         VCUT = VPLAIN(IP,KP)

         DO I = ILED, ITU(ID)

            IX = IX - 1
            XBASE(IX) = XINTR(I,ID)
            YBASE(IX) = YINTR(I,ID)
            ZBASE(IX) = ZINTR(I,ID)
            XYZCUT(1) = XINTR(I,ID)
            XYZCUT(3) = ZINTR(I,ID)

            CALL PLBICUT (KASE, ILWING, KLWING, 1, ILE, 1, KLWING,
     >                    XPLAIN, YPLAIN, ZPLAIN, UPLAIN, VPLAIN,
     >                    IP, KP, ISMAT, KSMAT, SMATRIX, EPS, DRANGE,
     >                    XYZCUT, UVRANGE, UCUT, VCUT, P, Q, LUNERR,
     >                    IER)

            IF (IER > 3) THEN
               JD = J - 1
               IF (I /= ITU(ID)) THEN ! Div. pt. at wing TE is problematic
                  GO TO 800           ! since wing surface is not extended now
               ELSE IF (IER > 4) THEN
                  GO TO 800
               END IF
            END IF

            UEVAL(IB,1) = UCUT
            VEVAL(IB,1) = VCUT
            IB = IB - 1
         END DO

C        Find (u,v)s for inboard half of pylon J intersection, starting at TE:

         IP   = IPDIVTE(J)
         KP   = KPDIVTE(J)
         UCUT = UPLAIN(IP,KP)
         VCUT = VPLAIN(IP,KP)
         ID2  = IDIV(J)
         ILED = NHALF(ID2)
         IB   = 1
         IX   = IXB(L) + NAFT + (NK - 1) * ISIZE

         DO I = ITL(ID2), ILED

            XBASE(IX) = XINTR(I,ID2)
            YBASE(IX) = YINTR(I,ID2)
            ZBASE(IX) = ZINTR(I,ID2)
            IX = IX + 1
            XYZCUT(1) = XINTR(I,ID2)
            XYZCUT(3) = ZINTR(I,ID2)

            CALL PLBICUT (KASE, ILWING, KLWING, 1, ILE, 1, KLWING,
     >                    XPLAIN, YPLAIN, ZPLAIN, UPLAIN, VPLAIN,
     >                    IP, KP, ISMAT, KSMAT, SMATRIX, EPS, DRANGE,
     >                    XYZCUT, UVRANGE, UCUT, VCUT, P, Q, LUNERR,
     >                    IER)

            IF (IER > 3) THEN
               JD = J
               IF (I /= ITL(ID2)) THEN ! Div. pt. at wing TE is problematic
                  GO TO 800            ! since wing surf. is not extended now
               ELSE IF (IER > 4) THEN
                  GO TO 800
               END IF
            END IF

            UEVAL(IB,NK) = UCUT
            VEVAL(IB,NK) = VCUT
            IB = IB + 1
         END DO

         UDIVLE(J) = UCUT ! For patches upstream
         VDIVLE(J) = VCUT

C        Spanwise edge (u,v)s between the diverter trailing edges, which
C        should match the subpatch aft done above if this is the pylon case.
C        Cranks between pylons aren't handled properly yet.

         IF (PYLON) THEN

            RVRANGE = ONE / (VPLAIN(1,KD2) - VPLAIN(1,KD))
            KP = KD
            DO K = 2, NK - 1
               KP = KP + 1
               R = (VPLAIN(1,KP) - VPLAIN(1,KD)) * RVRANGE
               UEVAL(1,K) = (ONE - R) * UEVAL(1,1) + R * UEVAL(1,NK)
               VEVAL(1,K) = (ONE - R) * VEVAL(1,1) + R * VEVAL(1,NK)
            END DO

         ELSE ! Diverter case: Hard to do anything but uniform:

            RV = (VEVAL(1,NK) - VEVAL(1,1)) / REAL (NK - 1)
            DO K = 2, NK - 1
               UEVAL(1,K) = ZERO
               VEVAL(1,K) = VEVAL(1,1) + RV * REAL (K - 1)
            END DO

         END IF

C        (u,v)s between diverter leading edges:

         R  = ONE / REAL (NK - 1)
         RU = R * (UEVAL(NI,NK) - UEVAL(NI,1))
         RV = R * (VEVAL(NI,NK) - VEVAL(NI,1))

         DO K = 2, NK - 1
            R = REAL (K - 1)
            UEVAL(NI,K) = UEVAL(NI,1) + R * RU
            VEVAL(NI,K) = VEVAL(NI,1) + R * RV
         END DO

C        Interior (u,v)s:

         CALL TFI2D (ILWING, 1, NI, 1, NK, UEVAL, VEVAL, TFIWORK)

C        Interpolate the interior (x,y,z)s:

         IX1 = IXB(L) + NAFT + ISTART - 1
         KP = KD

         DO K = 2, NK - 1

            R = REAL (K - 1) / REAL (NK - 1)
            IP = 1
            IF (PYLON) IP = NINT ((ONE - R) * REAL (IPDIVTE(J-1)) +
     >                                   R  * REAL (IPDIVTE(J)))
            IX1 = IX1 + ISIZE
            IX  = IX1

            DO I = ISTART, NI

               CALL PLINCUB (ILWING, KLWING, 1, ILIMIT, 1, KLWING,
     >                       XPLAIN, YPLAIN, ZPLAIN, UPLAIN, VPLAIN,
     >                       UEVAL(I,K), VEVAL(I,K), LINUV, LOOSE,
     >                       CLOSED, IP, KP, EPS,
     >                       XBASE(IX), YBASE(IX), ZBASE(IX), IER)
               IX = IX + 1
            END DO
         END DO

         KD = KD2
         ID = ID2

      END DO ! Next diverter J


C     Patch outboard of the outermost diverter:
C     -----------------------------------------

      L  = L + 1
      NI = ITU(ID) - ILED + 1
      NK = KLWING - KD + 1
      IB = IXB(L) + NAFT + ISTART - 1
      IP = IPDIVTE(NDIV)
      KP = KD
      IF (PYLON) KP = KPDIVTE(NDIV)
      K1 = KD
      KC = KCNUMBER(NDIV) + 1 ! # of first crank beyond last pylon
      K2 = MCRANK(KC)         ! Possibly the tip
      RVRANGE = ONE / (VPLAIN(1,K2) - VPLAIN(1,K1))
      UCUT    = UPLAIN(IP,KP)
      VCUT    = VPLAIN(IP,KP)
      ISIZE   = IMXCOM(L)

      DO I = ITU(ID) - ISTART + 1, ILED, -1

C        Determine (u,v) for this intersection point in the plain wing mesh:

         XBASE(IB) = XINTR(I,ID) ! (IB,1,L) where IB starts at ISTART
         YBASE(IB) = YINTR(I,ID)
         ZBASE(IB) = ZINTR(I,ID)
         XYZCUT(1) = XINTR(I,ID)
         XYZCUT(3) = ZINTR(I,ID)

         CALL PLBICUT (KASE, ILWING, KLWING, 1, ILE, 1, KLWING,
     >                 XPLAIN, YPLAIN, ZPLAIN, UPLAIN, VPLAIN,
     >                 IP, KP, ISMAT, KSMAT, SMATRIX, EPS, DRANGE,
     >                 XYZCUT, UVRANGE, UCUT, VCUT, P, Q, LUNERR, IER)

         IF (IER > 3) THEN
            JD = NDIV
            IF (.NOT. PYLON) THEN ! Diverter case, problematic at TE since the
               IF (I /= ITU(ID)) GO TO 800 ! wing surface is not extended now
            ELSE IF (IER > 4) THEN
               GO TO 800
            END IF
         END IF

C        Interpolate a corresponding spanwise line in the plain (u,v) mesh:

         PM1 = ONE - P

C        Transform the trailing edge spanwise distribn. out to the next crank.

         CALL UV_LINE_OUTBOARD (K1, K2, IP, KP, 1)

         J = K2 - K1 + 1
         DO K = K2, KLWING ! Any further panels are simpler
            UEVAL(1,J) = PM1 * UPLAIN(IP,K) + P * UPLAIN(IP+1,K)
            VEVAL(1,J) = PM1 * VPLAIN(IP,K) + P * VPLAIN(IP+1,K)
            J = J + 1
         END DO

         IX = IB

         DO K = 2, NK

            IX = IX + ISIZE ! (IB,K,L) where IB starts at ISTART

            CALL PLINCUB (ILWING, KLWING, 1, ILIMIT, 1, KLWING,
     >                    XPLAIN, YPLAIN, ZPLAIN, UPLAIN, VPLAIN,
     >                    UEVAL(1,K), VEVAL(1,K), LINUV, LOOSE, CLOSED,
     >                    IP, KP, EPS, XBASE(IX), YBASE(IX), ZBASE(IX),
     >                    IER)
         END DO

         IB = IB + 1

      END DO


C     Patches forward of diverter(s)?
C     -------------------------------

      I1 = MAX (IPDIVLE(1), IPDIVLE(NDIV))

      IF (I1 < ILE - 1) THEN ! ? Weak test

C        Mid-chord patch near the fuselage, in two pieces:

         L  = L - NDIV ! Back to the first under-wing patch
         I2 = (I1 + ILE + 1) / 2 ! Artificial mid-chord I boundary
         I1 = IPDIVLE(1) + 1
         NI = I2 - I1 + 2        ! For all mid-chord sub-subpatches
         MI = ISUB(IS,L)         ! I index of diverter LE
         NUMI  = ILE - I1 + 2    ! For full subpatches up to leading edge
         ISIZE = IMXCOM(L)
         NK    = JMXCOM(L)
         KD    = NK

C        Chordwise lines out to the crank inboard of pylon 1 (else K = 1 only):

         IX1 = IXB(L) + MI - 1
         DO K = 1, KBREAK
            IX = IX1
            IX1 = IX1 + ISIZE
            DO I = I1, ILE
               IX = IX + 1
               XBASE(IX) = XPLAIN(I,K)
               YBASE(IX) = YPLAIN(I,K)
               ZBASE(IX) = ZPLAIN(I,K)
            END DO
         END DO

C        Forward sub-subpatch for K = KBREAK+1:KD:

         IX1 = IXB(L) + KBREAK * ISIZE + MI + NI - 1 ! (MI+NI,KBRK+1,L)

         DO K = KBREAK + 1, KD
            IX = IX1
            IX1 = IX1 + ISIZE
            DO I = I2 + 1, ILE
               XBASE(IX) = XPLAIN(I,K)
               YBASE(IX) = YPLAIN(I,K)
               ZBASE(IX) = ZPLAIN(I,K)
               IX = IX + 1
            END DO
         END DO

C        Recover the diverter LE (p,q):

         UCUT = UDIVLE(1)
         VCUT = VDIVLE(1)
         IP   = IPDIVLE(1)
         KP   = KPDIVLE(1)

         CALL BILINT (ILWING, KLWING, UPLAIN, VPLAIN, UCUT, VCUT,
     >                IP, KP, EPS, P, Q, IER)

C        Set up the (u,v) edge through the diverter LE for K = KBREAK:KD.

         PM1 = ONE - P
         RVRANGE = ONE / (VPLAIN(1,KD) - VPLAIN(1,KBREAK))

         CALL UV_LINE_INBOARD (KBREAK, KD, IP, KP, 1)

         UEVAL(2:NI,KBREAK)  = UPLAIN(I1:I2,KBREAK) ! Chordwise edge
         VEVAL(2:NI,KBREAK)  = VPLAIN(I1:I2,KBREAK)
         UEVAL(NI,KBREAK:KD) = UPLAIN(I2,KBREAK:KD) ! Spanwise edge
         VEVAL(NI,KBREAK:KD) = VPLAIN(I2,KBREAK:KD)

C        Warp the nominal line in front of the diverter:

         I  = I1 - 1
         UE = UPLAIN(I,KD) ! Temporarily replace point I1-1 in the plain mesh
         VE = VPLAIN(I,KD)
         UPLAIN(I,KD) = PM1 * UPLAIN(I,KD) + P * UPLAIN(I1,KD)
         VPLAIN(I,KD) = PM1 * VPLAIN(I,KD) + P * VPLAIN(I1,KD)

C        NULINE2D assumes the same I range for input & output, so use KD+1:

         K = KD + 1
         UEVAL(I,K)  = UCUT
         VEVAL(I,K)  = VCUT
         UEVAL(I2,K) = UPLAIN(I2,KD)
         VEVAL(I2,K) = VPLAIN(I2,KD)

         CALL NULINE2D (I, I2, UPLAIN(1,KD), VPLAIN(1,KD),
     >                  UEVAL(1,K), VEVAL(1,K))

         UEVAL(1:NI,KD) = UEVAL(I:I2,K)
         VEVAL(1:NI,KD) = VEVAL(I:I2,K)
         UPLAIN(I,KD)   = UE ! Recover original values
         VPLAIN(I,KD)   = VE

C        Fill the interior (u,v)s between KBREAK & KD:

         CALL TFI2D (ILWING, 1, NI, KBREAK, KD, UEVAL, VEVAL, TFIWORK)

C        Evaluate the interior sub-patch & some edges:

         IX1 = IXB(L) + KBREAK * ISIZE + MI ! (MI+1,KBREAK+1,L)
         KP  = KBREAK + 1

         DO K = KBREAK + 1, KD
            IX  = IX1
            IX1 = IX1 + ISIZE
            IP  = IPDIVLE(1)
            DO I = 2, NI

               CALL PLINCUB (ILWING, KLWING, 1, ILIMIT, 1, KLWING,
     >                       XPLAIN, YPLAIN, ZPLAIN, UPLAIN, VPLAIN,
     >                       UEVAL(I,K), VEVAL(I,K), LINUV, LOOSE,
     >                       CLOSED, IP, KP, EPS,
     >                       XBASE(IX), YBASE(IX), ZBASE(IX), IER)
               IX = IX + 1
            END DO
         END DO


C        Subpatch(es) forward of diverters, between them:
C        ------------------------------------------------

         MK = NK  ! # Ks in preceding patch left

         DO J = 2, NDIV

            L     = L + 1
            ISIZE = IMXCOM(L)
            NK    = JMXCOM(L)
            MI2   = ISUB(IS,L) ! I index of diverter J LE

C           Common chordwise edge of patches L, L - 1:

            CALL XBCOPY2 (LENXB, XBASE, YBASE, ZBASE, IXB, IMXCOM,
     >                    L - 1, MI, IMXCOM(L-1), MK, MK, L, MI2, 1)

            DO I = 1, NI
               UEVAL(I,1) = UEVAL(I,MK)
               VEVAL(I,1) = VEVAL(I,MK)
            END DO

C           Spanwise edge of buffer sub-subpatch corresp. to plain I = I2 ...

            IX = IXB(L) + MI2 + NI - 2
            KP = KD

            DO K = 2, NK
               IX = IX + ISIZE
               KP = KP + 1
               XBASE(IX)   = XPLAIN(I2,KP)
               YBASE(IX)   = YPLAIN(I2,KP)
               ZBASE(IX)   = ZPLAIN(I2,KP)
               UEVAL(NI,K) = UPLAIN(I2,KP)
               VEVAL(NI,K) = VPLAIN(I2,KP)
            END DO

C           Adjust the number of nominal pts. in front of the current diverter:

            IX  = IXB(L) + (NK - 1) * ISIZE + MI2 - 1 ! (MI2,NK,L)
            IX1 = IX ! Needed for restoring LE pt.
            I1  = IPDIVLE(J) + 1
            KD2 = KD + NK - 1

            CALL CHANGEN
     >           (I1, I2, XPLAIN(1,KD2), YPLAIN(1,KD2), ZPLAIN(1,KD2),
     >            2,  NI, XBASE(IX), YBASE(IX), ZBASE(IX), LOOSE)

C           Corresponding (u,v)s, in U/VEVAL(*,NK+1) as nominal values:

            XT = XBASE(IX)
            YT = YBASE(IX)
            ZT = ZBASE(IX)
            IP = I1 - 1
            KP = KD2
            XBASE(IX) = XPLAIN(IP,KD2) ! Temporarily - grubby one!
            YBASE(IX) = YPLAIN(IP,KD2)
            ZBASE(IX) = ZPLAIN(IP,KD2)
            UCUT = UDIVLE(J)
            VCUT = VDIVLE(J)
            K = NK + 1

            DO I = 1, NI - 1
               XYZCUT(1) = XBASE(IX)
               XYZCUT(3) = ZBASE(IX)
               IX = IX + 1

               CALL PLBICUT (KASE, ILWING, KLWING, 1, ILE, 1, KLWING,
     >                       XPLAIN, YPLAIN, ZPLAIN, UPLAIN, VPLAIN,
     >                       IP, KP, ISMAT, KSMAT, SMATRIX, EPS, DRANGE,
     >                       XYZCUT, UVRANGE, UCUT, VCUT, P, Q, LUNERR,
     >                       IER)

               IF (IER > 3) THEN
                  JD = J
                  GO TO 800
               END IF

               UEVAL(I,K) = UCUT
               VEVAL(I,K) = VCUT
            END DO

            UEVAL(NI,K) = UEVAL(NI,NK)
            VEVAL(NI,K) = VEVAL(NI,NK)

            XBASE(IX1) = XT ! Restore the diverter LE pt.
            YBASE(IX1) = YT
            ZBASE(IX1) = ZT

C           Warp the nominal (u,v)s to pass through the true diverter LE (u,v):

            UEVAL(1,NK) = UDIVLE(J)
            VEVAL(1,NK) = VDIVLE(J)

            CALL NULINE2D (1, NI, UEVAL(1,K), VEVAL(1,K),
     >                     UEVAL(1,NK), VEVAL(1,NK))

C           (u,v) edge between diverter leading edges as done above:

            R  = ONE / REAL (NK - 1)
            RU = R * (UEVAL(1,NK) - UEVAL(1,1))
            RV = R * (VEVAL(1,NK) - VEVAL(1,1))

            DO K = 2, NK - 1
               R = REAL (K - 1)
               UEVAL(1,K) = UEVAL(1,1) + R * RU
               VEVAL(1,K) = VEVAL(1,1) + R * RV
            END DO

C           Fill the interior (u,v)s:

            CALL TFI2D (ILWING, 1, NI, 1, NK, UEVAL, VEVAL, TFIWORK)

C           Interpolate the interior (x,y,z)s:

            IX1 = IXB(L) + ISIZE + MI2 ! (MI2+1,2,L)
            KP  = KD
            Q   = ONE / REAL (NK - 1)

            DO K = 2, NK

               IX  = IX1
               IX1 = IX1 + ISIZE
               R   = REAL (K - 1) * Q
               IP  = NINT ((ONE - R) * REAL (IPDIVLE(J-1)) +
     >                            R  * REAL (IPDIVLE(J)))
               DO I = 2, NI - 1

                  CALL PLINCUB (ILWING, KLWING, 1, ILIMIT, 1, KLWING,
     >                          XPLAIN, YPLAIN, ZPLAIN, UPLAIN, VPLAIN,
     >                          UEVAL(I,K), VEVAL(I,K), LINUV, LOOSE,
     >                          CLOSED, IP, KP, EPS, XBASE(IX),
     >                          YBASE(IX), ZBASE(IX), IER)
                  IX = IX + 1
               END DO
            END DO

C           Rest of this patch forward of the "buffer" subpatch:

            IX1 = IXB(L) + MI2 + NI - 2 ! (MI2 + NI - 1,1,L)
            KP  = KD
            DO K = 2, NK
               IX1 = IX1 + ISIZE
               IX  = IX1
               KP  = KP + 1
               DO I = I2 + 1, ILE
                  IX = IX + 1
                  XBASE(IX) = XPLAIN(I,KP)
                  YBASE(IX) = YPLAIN(I,KP)
                  ZBASE(IX) = ZPLAIN(I,KP)
               END DO
            END DO

            MI = MI2
            MK = NK
            KD = KD2

         END DO ! Next diverter J


C        Patch outboard/forward of the outermost diverter:
C        -------------------------------------------------

         L     = L + 1
         ISIZE = IMXCOM(L)
         MI2   = ISUB(IS,L)
         NK    = JMXCOM(L)
         K1    = KD
         KC    = KCNUMBER(NDIV) + 1
         K2    = MCRANK(KC) ! Possibly the tip
         NUMK  = K2 - K1 + 1 ! For TFI region

C        Common chordwise edge:

         CALL XBCOPY2 (LENXB, XBASE, YBASE, ZBASE, IXB, IMXCOM,
     >                 L - 1, MI, IMXCOM(L-1), MK, MK, L, MI2, 1)

         DO I = 1, NI
            UEVAL(I,1) = UEVAL(I,MK)
            VEVAL(I,1) = VEVAL(I,MK)
         END DO

         IX1 = IXB(L) + ISIZE + MI2 + NI - 2 ! (MI2 + NI - 1,2,L)
         KT  = 1

         DO K = K1 + 1, KLWING
            IX = IX1
            IX1 = IX1 + ISIZE
            DO I = I2, ILE
               XBASE(IX) = XPLAIN(I,K) ! Subpatch forward of
               YBASE(IX) = YPLAIN(I,K) ! the "buffer" line I2
               ZBASE(IX) = ZPLAIN(I,K)
               IX = IX + 1
            END DO
            KT = KT + 1
            UEVAL(NI,KT) = UPLAIN(I2,K)
            VEVAL(NI,KT) = VPLAIN(I2,K)
         END DO

C        Spanwise edge through outboard pylon LE out to next crank:

         UCUT = UDIVLE(NDIV)
         VCUT = VDIVLE(NDIV)
         IP   = IPDIVLE(NDIV)
         KP   = KPDIVLE(NDIV)

         CALL BILINT (ILWING, KLWING, UPLAIN, VPLAIN, UCUT, VCUT, ! Recover p
     >                IP, KP, EPS, P, Q, IER)

         PM1 = ONE - P
         RVRANGE = ONE / (VPLAIN(1,K2) - VPLAIN(1,K1))

         CALL UV_LINE_OUTBOARD (K1, K2, IP, KP, 1)

C        Adjust the number of points along the crank/tip in the buffer region:

         IX = IXB(L) + (NUMK - 1) * ISIZE + MI2 - 1 ! (MI2,NUMK,L)

         CALL CHANGEN (I1, I2, XPLAIN(1,K2), YPLAIN(1,K2), ZPLAIN(1,K2),
     >                 2,  NI, XBASE(IX), YBASE(IX), ZBASE(IX), LOOSE)

C        We need the corresponding (u,v)s along the crank/tip.
C        Impose this relative (u,v) distribution on lines K = crank to tip
C        in the plain wing mesh.  Do it now because we need the K2 result
C        to do TFI between K1 and K2.

         IP = I1
         KP = K2
         UCUT = UPLAIN(I1,K2)
         VCUT = VPLAIN(I1,K2)

         DO I = 2, NI - 1

C           Deflected slats seem to have caused PLBICUT trouble - use 1D method:

            IX = IX + 1
            XYZCUT(1) = XBASE(IX)

            CALL PLXCUT (ILE, XPLAIN(1,K2), YPLAIN(1,K2), UPLAIN(1,K2),
     >                   IP, .FALSE., XYZCUT(1), EPS, UCUT, LUNERR, IER)

            IF (IER > 2) THEN
               JD = NDIV + 1 ! Reuse the diverter error message
               GO TO 800
            END IF

            P = (UCUT            - UPLAIN(IP,K2)) /
     >          (UPLAIN(IP+1,K2) - UPLAIN(IP,K2))
            PM1 = ONE - P
            VCUT = PM1 * VPLAIN(IP,K2) + P * VPLAIN(IP+1,K2)

            UEVAL(I,NUMK) = UCUT
            VEVAL(I,NUMK) = VCUT

            J = NUMK + 1
            DO K = K2 + 1, KLWING
               UEVAL(I,J) = PM1 * UPLAIN(IP,K) + P * UPLAIN(IP+1,K)
               VEVAL(I,J) = PM1 * VPLAIN(IP,K) + P * VPLAIN(IP+1,K)
               J = J + 1
            END DO

         END DO

C        Interior (u,v)s forward of outer pylon, between pylon and crank/tip:

         CALL TFI2D (ILWING, 1, NI, 1, NUMK, UEVAL, VEVAL, TFIWORK)

C        Interpolate (x,y,z)s between outer diverter and tip (forward of div.):

         IX1 = IXB(L) + ISIZE + MI2 ! (MI2 + 1,2,L)
         J   = 2
         KP  = K1

         DO K = K1 + 1, KLWING
            IX  = IX1
            IX1 = IX1 + ISIZE
            IP  = I1
            DO I = 2, NI - 1

               CALL PLINCUB (ILWING, KLWING, 1, ILE, 1, KLWING,
     >                       XPLAIN, YPLAIN, ZPLAIN, UPLAIN, VPLAIN,
     >                       UEVAL(I,J), VEVAL(I,J), LINUV, LOOSE,
     >                       CLOSED, IP, KP, EPS, XBASE(IX),
     >                       YBASE(IX), ZBASE(IX), IER)
               IX = IX + 1
            END DO
            J = J + 1
         END DO

      END IF


C     Upper wing, split spanwise to match diverter/pylon splits:
C     ----------------------------------------------------------

      L  = L + 1

      MI = ISUB(1,L)
      IL = ILE + MI - 1 ! Plain wing index for upper triangle apexes, or ILWING
      K1 = 1

      DO J = 1, NDIV + 1

         NK  = JMXCOM(L)
         K2  = K1 + NK - 1
         IX1 = IXB(L)

         DO K = K1, K2
            IX = IX1
            IX1 = IX1 + ILE
            DO I = ILE, IL
               XBASE(IX) = XPLAIN(I,K)
               YBASE(IX) = YPLAIN(I,K)
               ZBASE(IX) = ZPLAIN(I,K)
               IX = IX + 1
            END DO
         END DO

         K1 = K2
         L  = L + 1
      END DO


C     More awkward paneling along the upper wing trailing edge?
C     ---------------------------------------------------------
C     The diverter case needs to match lower wing patch corners at the TE.
C     The patches along the upper TE are done in order from inboard to tip.

      IF (.NOT. PYLON) THEN

         LWING = L             ! For NACSURF2 to match with the 1st nacelle cap
         LCAP  = L
         L     = L1 + NDIV + 1 ! Back to the first upper wing patch
         M     = L1            ! For matching upper TE pts. with lower TE pts.
         IX1   = IXB(L) + MI - 1
         NI    = ILWING - IL + 1
         IM    = IL - 1        ! To retain cubics in the chordwise dirn. near IL
         NK    = JMXCOM(L)

         DO K = 1, KBREAK ! 1 <= KBREAK < NK (nearest crank to diverter if any)
            IX  = IX1
            IX1 = IX1 + ILE
            IP  = 1
            DO I = IL, ILWING
               XBASE(IX)   = XPLAIN(I,K)
               YBASE(IX)   = YPLAIN(I,K)
               ZBASE(IX)   = ZPLAIN(I,K)
               UEVAL(IP,K) = UPLAIN(I,K)
               VEVAL(IP,K) = VPLAIN(I,K)
               IP = IP + 1
               IX = IX + 1
            END DO
         END DO

         IX = IXB(L) + KBREAK * ILE + MI - 1 ! (MI,KBREAK+1,L)

         DO K = KBREAK + 1, NK ! Rest of upstream edge
            XBASE(IX)  = XPLAIN(IL,K)
            YBASE(IX)  = YPLAIN(IL,K)
            ZBASE(IX)  = ZPLAIN(IL,K)
            UEVAL(1,K) = UPLAIN(IL,K)
            VEVAL(1,K) = VPLAIN(IL,K)
            IX = IX + ILE
         END DO

C        Interpolate the rest of upper wing TE at the lower wing patch TE pts.:

         DO K = 1, KLWING ! Upper wing TE
            XLINE(K) = XPLAIN(ILWING,K)
            YLINE(K) = YPLAIN(ILWING,K)
            ZLINE(K) = ZPLAIN(ILWING,K)
         END DO

         KP = KBREAK

         CALL LOWER_TO_UPPER_TE (KP + 1, NK) ! From lower patch M to upper L

C        Construct an edge in the upper wing (u,v) space which roughly
C        shadows the inboard side of the first diverter on the lower wing:

         CALL UV_EDGE_UPPER (NK, NK)

C        Fill the interior (u,v)s:

         CALL TFI2D (ILWING, 1, NI, KBREAK, NK, UEVAL, VEVAL, TFIWORK)

C        Evaluate the (x,y,z)s:

         IX1 = IXB(L) + KBREAK * ILE + MI ! (MI+1,KBREAK+1,L)
         KP  = KBREAK

         DO K = KBREAK + 1, NK
            IX  = IX1
            IX1 = IX1 + ILE
            IP  = IL
            DO I = 2, NI - 1

               CALL PLINCUB (ILWING, KLWING, IM, ILWING, 1, KLWING,
     >                       XPLAIN, YPLAIN, ZPLAIN, UPLAIN, VPLAIN,
     >                       UEVAL(I,K), VEVAL(I,K), LINUV, LOOSE,
     >                       CLOSED, IP, KP, EPS, XBASE(IX),
     >                       YBASE(IX), ZBASE(IX), IER)
               IX = IX + 1
            END DO
         END DO

C        Handle patches corresponding to diverters, and in-between patches:
C        ------------------------------------------------------------------

         KD = NK
         KC = KCNUMBER(1)

         DO J = 1, NDIV

C           The trailing edge of patch L + 1 gives a point of the triangle:

            L  = L + 1 ! Upper patch outboard of diverter J
            M  = M + 1 ! Lower wing patch outboard of diverter J
            KP = KD    ! All these are needed by LOWER_TO_UPPER_TE

            CALL LOWER_TO_UPPER_TE (1, 1) ! From lower patch M to upper patch L

            CALL UV_EDGE_UPPER (KD, 1)    ! Inb. edge of L = outer triangle edge

C           (NI,NKTRI,LCAP) <-- (ILE,1,L)  ! TE pt. of outer edge of triangle

            CALL XBCOPY2 (LENXB, XBASE, YBASE, ZBASE, IXB, IMXCOM,
     >                    L, ILE, ILE, 1, 1, LCAP, NI, NKTRI)

            DO I = 2, NI ! Outboard edge of triangle
               UEVAL(I,NKTRI) = UEVAL(I,1)
               VEVAL(I,NKTRI) = VEVAL(I,1)
            END DO

C           (1:NI,1,LCAP) <-- (MI:ILE,NK,L-1) ! Inboard edge of the triangle

            CALL XBCOPY2 (LENXB, XBASE, YBASE, ZBASE, IXB, IMXCOM,
     >                    L - 1, MI, ILE, NK, NK, LCAP, 1, 1)

            DO I = 1, NI
               UEVAL(I,1) = UEVAL(I,NK)
               VEVAL(I,1) = VEVAL(I,NK)
            END DO

            IX1 = IXB(LCAP)
            IX  = IX1
            DO K = 2, NKTRI ! Singular point of the triangle
               IX = IX + NI
               XBASE(IX)  = XBASE(IX1)
               YBASE(IX)  = YBASE(IX1)
               ZBASE(IX)  = ZBASE(IX1)
               UEVAL(1,K) = UEVAL(1,1)
               VEVAL(1,K) = VEVAL(1,1)
            END DO

C           Interpolate the upper trailing edge along the base of the triangle:

            DV = (VEVAL(NI,NKTRI) - VEVAL(NI,1)) / REAL (NKTRI - 1)

            DO K = 2, NKTRI - 1
               UEVAL(NI,K) = ONE
               VEVAL(NI,K) = VEVAL(NI,1) + DV * REAL (K - 1)
            END DO

C           Fill the interior (u,v)s of the triangle:

            CALL TFI2D (ILWING, 1, NI, 1, NKTRI, UEVAL, VEVAL, TFIWORK)

C           Evaluate the unassigned (x,y,z)s of the triangular patch:

            IX1 = IX1 + NI + 1 ! (2,2,LCAP)
            KP  = KD

            DO K = 2, NKTRI
               IX  = IX1
               IX1 = IX1 + NI
               IP  = IL
               DO I = 2, NI

                  CALL PLINCUB (ILWING, KLWING, IM, ILWING, 1, KLWING,
     >                          XPLAIN, YPLAIN, ZPLAIN, UPLAIN, VPLAIN,
     >                          UEVAL(I,K), VEVAL(I,K), LINUV, LOOSE,
     >                          CLOSED, IP, KP, EPS, XBASE(IX),
     >                          YBASE(IX), ZBASE(IX), IER)
                  IX = IX + 1
               END DO
            END DO

C           Subpatch between triangular patches?

            IF (J < NDIV) THEN

               KC  = KCNUMBER(J+1)
               KD2 = MCRANK(KC)
               NK  = KD2 - KD + 1

C              (MI:ILE,1,L) <-- (1:NI,NKTRI,LCAP) ! Edge common to triangle

               CALL XBCOPY2 (LENXB, XBASE, YBASE, ZBASE, IXB, IMXCOM,
     >                       LCAP, 1, NI, NKTRI, NKTRI, L, MI, 1)

               DO I = 1, NI
                  UEVAL(I,1) = UEVAL(I,NKTRI)
                  VEVAL(I,1) = VEVAL(I,NKTRI)
               END DO

               IX = IXB(L) + MI - 1
               KT = 1

               DO K = KD, KD2 ! Upstream edge
                  XBASE(IX)   = XPLAIN(IL,K)
                  YBASE(IX)   = YPLAIN(IL,K)
                  ZBASE(IX)   = ZPLAIN(IL,K)
                  UEVAL(1,KT) = UPLAIN(IL,K)
                  VEVAL(1,KT) = VPLAIN(IL,K)
                  KT = KT + 1
                  IX = IX + ILE
               END DO

               CALL LOWER_TO_UPPER_TE (2, NK) ! From lower M to upper L

               CALL UV_EDGE_UPPER (KD2, NK) ! Outer edge of in-between patch

               CALL TFI2D (ILWING, 1, NI, 1, NK, UEVAL, VEVAL, TFIWORK)

C              Evaluate the unassigned (x,y,z)s of the in-between patch:

               IX1 = IXB(L) + MI ! (MI+1,1,L)
               KP  = KD

               DO K = 2, NK
                  IX1 = IX1 + ILE
                  IX  = IX1
                  IP  = IL
                  DO I = 2, NI - 1

                     CALL PLINCUB (ILWING, KLWING, IM, ILWING, 1,KLWING,
     >                             XPLAIN,YPLAIN,ZPLAIN, UPLAIN, VPLAIN,
     >                             UEVAL(I,K), VEVAL(I,K), LINUV, LOOSE,
     >                             CLOSED, IP, KP, EPS, XBASE(IX),
     >                             YBASE(IX), ZBASE(IX), IER)
                     IX = IX + 1
                  END DO
               END DO

               KD = KD2

            END IF

            LCAP = LCAP + 1

         END DO ! Next triangular patch

C        Subpatch outboard of outermost triangular patch:

         NK = KLWING - KD + 1
         K1 = KD
         KP = KD ! For LOWER_TO_UPPER_TE
         KC = KC + 1
         K2 = MCRANK(KC) ! Possibly the tip
         NUMK = K2 - K1 + 1 ! For TFI region out to nearest crank
         LCAP = LCAP - 1

C        (MI:ILE,1,L) <-- (1:NI,NKTRI,LCAP) ! Edge common to triangle

         CALL XBCOPY2 (LENXB, XBASE, YBASE, ZBASE, IXB, IMXCOM,
     >                 LCAP, 1, NI, NKTRI, NKTRI, L, MI, 1)

         DO I = 1, NI
            UEVAL(I,1) = UEVAL(I,NKTRI)
            VEVAL(I,1) = VEVAL(I,NKTRI)
         END DO

         KT = 1
         DO K = K1, K2 ! (Part of) upstream spanwise edge
            UEVAL(1,KT) = UPLAIN(IL,K)
            VEVAL(1,KT) = VPLAIN(IL,K)
            KT = KT + 1
         END DO

         CALL LOWER_TO_UPPER_TE (2, NUMK - 1) ! Adapt TE from lower M to upper L

C        Crank:

         IP = 1
         DO I = IL, ILWING
            UEVAL(IP,NUMK) = UPLAIN(I,K2)
            VEVAL(IP,NUMK) = VPLAIN(I,K2)
            IP = IP + 1
         END DO

         CALL TFI2D (ILWING, 1, NI, 1, NUMK, UEVAL, VEVAL, TFIWORK)

         IX1 = IXB(L) + MI ! (MI+1,1,L)
         KP  = KD

         DO K = 2, NUMK - 1
            IX1 = IX1 + ILE
            IX  = IX1
            IP  = IL
            DO I = 2, NI - 1

               CALL PLINCUB (ILWING, KLWING, IM, ILWING, 1, KLWING,
     >                       XPLAIN,YPLAIN,ZPLAIN, UPLAIN, VPLAIN,
     >                       UEVAL(I,K), VEVAL(I,K), LINUV, LOOSE,
     >                       CLOSED, IP, KP, EPS, XBASE(IX),
     >                       YBASE(IX), ZBASE(IX), IER)
               IX = IX + 1
            END DO
         END DO

C        Crank to tip outboard of the last triangle:

         IX1 = IXB(L) + (NUMK - 1) * ILE + MI ! (MI + 1,NUMK,L)

         DO K = K2, KLWING
            IX  = IX1
            IX1 = IX1 + ILE
            DO I = IL + 1, ILWING
               XBASE(IX) = XPLAIN(I,K)
               YBASE(IX) = YPLAIN(I,K)
               ZBASE(IX) = ZPLAIN(I,K)
               IX = IX + 1
            END DO
         END DO

         L = LCAP + 1

      END IF ! Diverter case, matching lower surface patches at wing TE

      KPATCH = L ! Next empty patch
      FAIL = .FALSE.
      GO TO 999


  800 IF (IWRIT > 0) THEN
         WRITE (IWRIT, '(/,A,I2,/,A,I2,/,A,I4,/,A,1P,3E15.6)')
     >      ' WDSURF: IER from PLBICUT =', IER,
     >      ' Pylon/diverter # J =', JD,
     >      ' Intersection pt. I =', I,
     >      ' Intersection X,Y,Z =', XYZCUT
      END IF

  900 FAIL = .TRUE.

  999 RETURN


C     WDSURF internal procedures:

      CONTAINS

!        -------------------------------------
         SUBROUTINE LOWER_TO_UPPER_TE (K1, K2)
!        -------------------------------------

!        Interpolate upper wing trailing edge points K1:K2 for patch L at the
!        lower wing trailing edge points from patch M (both (x,y,z) and (u,v)).

!        Arguments:

         INTEGER, INTENT (IN) :: K1, K2

!        Local variables:

         INTEGER IXL, IXM, K
         REAL    ZTE

!        Execution:

!        Control KP at the higher level.

         IXM = IXB(M) + (K1 - 1) * IMXCOM(M) ! (1,K1,M)
         IXL = IXB(L) + K1 * IMXCOM(L) - 1   ! (NI,K1,L)

         DO K = K1, K2

C           Using the lower surface XTE is not good enough for the upper XTE.

            ZTE = ZBASE(IXM) ! On lower surface
            IXM = IXM + IMXCOM(M)

            CALL INTERVAL (KLWING, ZLINE, ZTE, ONE, KP)

            R = (ZTE - ZLINE(KP)) / (ZLINE(KP+1) - ZLINE(KP))
            XBASE(IXL)  = (ONE - R) * XLINE(KP) + R * XLINE(KP+1)
            YBASE(IXL)  = (ONE - R) * YLINE(KP) + R * YLINE(KP+1)
            ZBASE(IXL)  =  ZTE
            IXL = IXL + IMXCOM(L)
            UEVAL(NI,K) =  ONE
            VEVAL(NI,K) = (ONE - R) * VPLAIN(ILWING,KP) +
     >                           R  * VPLAIN(ILWING,KP+1)
         END DO

         END SUBROUTINE LOWER_TO_UPPER_TE

!        -----------------------------------------
         SUBROUTINE UV_EDGE_UPPER (KPLAIN, KPANEL)
!        -----------------------------------------

!        Establish a plain wing (u,v) edge for an upper wing panel which shadows
!        a lower wing diverter region centered at K = KPLAIN.  KPANEL = 1 or NK.
!        End-points 1 & NI corresp. to IL & ILWING are assumed to be in place.

!        Arguments:

         INTEGER, INTENT (IN) :: KPLAIN, KPANEL

!        Local variables:

         INTEGER I, IP
         REAL    DENOM, ULEFT

!        Execution:

         ULEFT = UPLAIN(IL,KPLAIN)
         DENOM = ONE / (ONE - ULEFT)
         I = 1

         DO IP = IL + 1, ILWING - 1 ! Any reasonable line will do
            I = I + 1
            UEVAL(I,KPANEL) = UPLAIN(IP,KPLAIN)
            R = (UEVAL(I,KPANEL) - ULEFT) * DENOM
            VEVAL(I,KPANEL) = (ONE - R) * VPLAIN(IL,KPLAIN) +
     >                        R * VEVAL(NI,KPANEL)
         END DO

         END SUBROUTINE UV_EDGE_UPPER

!        -----------------------------------------------
         SUBROUTINE UV_LINE_INBOARD (K1, K2, IP, KP, IE)
!        -----------------------------------------------

!        Impose the trailing edge V distribution between K1 & K2 on a
!        spanwise line interpolated between plain mesh lines IP & IP + 1
!        according to the values of IP, KP, P, Q, UCUT, VCUT.
!        Results are returned in U/VEVAL(IE,K1:K2).
!        This is more thorough than working with just the end points as
!        NULINE2D does.  (Remember, U = 0 all along the trailing edge.)

!        Arguments:

         INTEGER, INTENT (IN) :: K1, K2, IP, KP, IE

!        Local variables:

         INTEGER K, KL, KLEFT, NLINE
         REAL    R

!        Establish a plain mesh (u,v) data line to interpolate within:

         DO K = K1, KP
            ULINE(K) = PM1 * UPLAIN(IP,K) + P * UPLAIN(IP+1,K)
            VLINE(K) = PM1 * VPLAIN(IP,K) + P * VPLAIN(IP+1,K)
         END DO

         NLINE = KP - K1 + 1
         IF (ABS (Q) > EPS) THEN
            NLINE = NLINE + 1
            ULINE(KP+1) = UCUT
            VLINE(KP+1) = VCUT
         END IF

         VSCALE = RVRANGE * (VCUT - VLINE(K1))
         VSHIFT = RVRANGE * (VPLAIN(1,K2) * VLINE(K1) -
     >                       VPLAIN(1,K1) * VCUT)
         KLEFT = 1

         DO K = K1, K2 - 1
            VEVAL(IE,K) = VPLAIN(1,K) * VSCALE + VSHIFT

            CALL INTERVAL (NLINE, VLINE(K1), VEVAL(IE,K), ONE, KLEFT)

            KL = KLEFT + K1 - 1
            R  = (VEVAL(IE,K) - VLINE(KL)) / (VLINE(KL+1) - VLINE(KL))
            UEVAL(IE,K) = (ONE - R) * ULINE(KL) + R * ULINE(KL+1)
         END DO

         UEVAL(IE,K2) = UCUT
         VEVAL(IE,K2) = VCUT

         END SUBROUTINE UV_LINE_INBOARD

!        ------------------------------------------------
         SUBROUTINE UV_LINE_OUTBOARD (K1, K2, IP, KP, IE)
!        ------------------------------------------------

!        Analogue of UV_LINE_INBOARD for lines beyond outer pylon/diverter.
!        Results are returned in U/VEVAL(IE,1:K2-K1+1).

!        Arguments:

         INTEGER, INTENT (IN) :: K1, K2, IP, KP, IE

!        Local variables:

         INTEGER J, K, KL, KLEFT, NLINE
         REAL    R

!        Establish a plain mesh (u,v) data line to interpolate within:

         ULINE(KP) = UCUT
         VLINE(KP) = VCUT

         DO K = KP + 1, K2
            ULINE(K) = PM1 * UPLAIN(IP,K) + P * UPLAIN(IP+1,K)
            VLINE(K) = PM1 * VPLAIN(IP,K) + P * VPLAIN(IP+1,K)
         END DO

         VSCALE = RVRANGE * (VLINE(K2) - VCUT)
         VSHIFT = RVRANGE * (VPLAIN(1,K2) * VCUT -
     >                       VPLAIN(1,K1) * VLINE(K2))

         UEVAL(IE,1) = UCUT
         VEVAL(IE,1) = VCUT

         KLEFT = 1
         NLINE = K2 - KP + 1
         J = 1

         DO K = K1 + 1, K2
            J = J + 1
            VEVAL(IE,J) = VPLAIN(1,K) * VSCALE + VSHIFT

            CALL INTERVAL (NLINE, VLINE(KP), VEVAL(IE,J), ONE, KLEFT)

            KL = KLEFT + KP - 1
            R  = (VEVAL(IE,J) - VLINE(KL)) / (VLINE(KL+1) - VLINE(KL))
            UEVAL(IE,J) = (ONE - R) * ULINE(KL) + R * ULINE(KL+1)
         END DO

         END SUBROUTINE UV_LINE_OUTBOARD

      END SUBROUTINE WDSURF

C*******************************************************************************
C
      SUBROUTINE WFINTR (NCALL, NNCONFIG, IW, LOSTINTR,
     >                   KWING1, KWING2, MXIWING, MXKWING, MXINTR,
     >                   NHALF, ILWING, XRG, YRG, ZRG,
     >                   NJF, NIF, MXJBODY, MXIBODY, XF, YF, ZF,
     >                   ISTN1, ISTN2, IDEGFUS, NINTR,
     >                   XINTR, YINTR, ZINTR, UINTR, VINTR, TINTR,
     >                   IINTR, JINTR, KINTR, IWRIT, FAIL)
C
C     WFINTR calculates a wing/fuselage intersection given a regular
C     geometry mesh for each.  The version from SYN87-SB just returned
C     the intersection T values, producing the intersection coordinates
C     as part of its arc-based mesh generation.  But SURFER panels the
C     body before the wings, so intersection (x,y,z)s are now returned too.
C
C     09/26/97  DAS  SYN87 adaptation of the earlier SYN87 routine for use
C                    AFTER regularizing the wing sections, not before.
C     10/22/97   "   SURFER needs X/Y/ZINTR (not just TINTR).
C     11/22/97   "   Made use of generalized WBINTR; avoided some of the
C                    initialization after the first call.
C     01/21/98   "   LOSTINTR argument needed in case part of the intersection
C                    is updated by WDINTR, preventing reuse of initial results.
C     07/19/00   "   NNCONFIG needed for CTV case (full body for intersections).
C
C*******************************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   NCALL,                     ! -3 means first call for this case
     >   NNCONFIG,                  ! > 0 means CTV cases (whole body)
     >   IW                         ! Wing surface # (for LOSTINTR case)
      LOGICAL, INTENT (IN) ::
     >   LOSTINTR                   ! T means [re]store intersection LE data
      INTEGER, INTENT (IN) ::
     >   KWING1, KWING2,            ! Subset of wing sections to treat
     >   MXIWING, MXKWING,          ! Dimensions of wing section arrays
     >   MXINTR,                    ! Dimension of intersection arrays
     >   NHALF, ILWING              ! Leading & upper trailing edge indices
      REAL, INTENT (IN), DIMENSION (MXIWING,MXKWING) ::
     >   XRG, YRG, ZRG              ! Regularized wing sections
      INTEGER, INTENT (IN) ::
     >   NJF,                       ! Active circumferential and axial
     >   NIF,                       ! indices of regular fuselage sections
     >   MXJBODY,                   ! J and I dimensions of fuselage
     >   MXIBODY                    ! section arrays
      REAL, INTENT (IN) ::
     >   XF(MXIBODY),               ! Regularized fuselage sections
     >   YF(MXJBODY,MXIBODY),
     >   ZF(MXJBODY,MXIBODY)
      INTEGER, INTENT (INOUT) ::
     >   ISTN1,                     ! Fuselage station range for this wing,
     >   ISTN2                      ! output on first call, input thereafter
      INTEGER, INTENT (IN) ::
     >   IDEGFUS                    ! 3 means bicubics (INTSEC4) else
                                    ! 1 means bilinear (INTSEC5)
      INTEGER, INTENT (OUT) ::
     >   NINTR                      ! No. of pts. in the intersection
                                    ! (ILWING, but may change in WDSURF)
      REAL, INTENT (OUT), DIMENSION (MXINTR) ::
     >   XINTR, YINTR, ZINTR,       ! Intersection coordinates in 1:ILWING
     >   UINTR,                     ! Corresponding body surface "u" and
     >   VINTR,                     ! "v" values (circumferential & axial)
     >   TINTR                      ! "t" values at the intersection
      INTEGER, INTENT (INOUT), DIMENSION (MXINTR) ::
     >   IINTR, JINTR, KINTR        ! Corresponding body & wing cell indices
      INTEGER, INTENT (IN) ::
     >   IWRIT                      ! For diagnostics if > 0
      LOGICAL, INTENT (OUT) ::
     >   FAIL                       ! Success flag

C     Local constants:

      REAL, PARAMETER :: ONE = 1.

C     Local variables:

      INTEGER
     >   I, I1, I2, IFAIL, J
      REAL
     >   YFMN, YFMX, ZFMX
      REAL, DIMENSION (MXJBODY,NIF) ::
     >   UBODY, VBODY, XBODY
      LOGICAL
     >   FIRST

C     Local leading edge data saved for LOSTINTR case.  10 >= NNWING

      INTEGER, DIMENSION (10), SAVE :: IINTLE, JINTLE, KINTLE
      REAL,    DIMENSION (10), SAVE :: TINTLE, UINTLE, VINTLE

C     Execution:

      FAIL  = .FALSE. ! Allow retrying in WBINTR
      NINTR = ILWING  ! But WDSURF may update this in the LOSTINTR case
      FIRST = NCALL == -3

      IF (FIRST) THEN

C        Determine a subset of body sections containing the intersection:

         I1 = NIF / 5
         CALL INTERVAL (NIF, XF, XRG(NHALF,KWING1), ONE, I1) ! LE of KWING1 stn.

         I2 = I1
         CALL INTERVAL (NIF, XF, XRG(NHALF,KWING1+1), ONE, I2) ! LE of next stn.

         ISTN1 = MAX (1, MIN (I1, I2) - 2)

         IINTR(NHALF) = (I1 + I2) / 2 ! Estimate of body stn. for LE of intersn.

         I1 = 4 * NIF / 5
         CALL INTERVAL (NIF, XF, XRG(1,KWING1), ONE, I1)     ! TE of KWING1 stn.

         I2 = I1
         CALL INTERVAL (NIF, XF, XRG(1,KWING1+1), ONE, I2)   ! TE of next stn.

         ISTN2 = MIN (NIF, MAX (I1, I2) + 3)   ! INTERVAL always points to left

         IINTR(1) = (I1 + I2) / 2     ! Estimate of body stn. for TE of intersn.

         IF (NNCONFIG == 0) THEN ! Pre-CTV cases start at LE
            I2 = IINTR(NHALF)
         ELSE ! For starting at the TE
            I2 = IINTR(1)
         END IF

         YFMN = 1.E+10 ! For body stn. ~ initial intersection pt.
         YFMX = -YFMN
         ZFMX =  YFMX

         DO J = 1, NJF
            YFMN = MIN (YF(J,I2), YFMN)
            YFMX = MAX (YF(J,I2), YFMX)
            ZFMX = MAX (ZF(J,I2), ZFMX)
         END DO

         IF (NNCONFIG == 0) THEN ! Pre-CTV cases

C           Set up for starting from the leading edge:

            I1 = NHALF

            UINTR(I1) = (YRG(I1,KWING1) - YFMN) / (YFMX - YFMN)
            UINTR(I1) = MAX (0.3, MIN (UINTR(I1), 0.7)) ! Y is poor for LE u
            VINTR(I1) = (XF(I2) - XF(ISTN1)) / (XF(ISTN2) - XF(ISTN1))
            TINTR(I1) = (ZFMX           - ZRG(I1,KWING1)) /
     >                  (ZRG(I1,KWING2) - ZRG(I1,KWING1))
            TINTR(I1) = MAX (0.1, MIN (TINTR(I1), 0.9))
            JINTR(I1) = NJF * UINTR(I1)
            KINTR(I1) = (KWING1 + KWING2) / 2

         ELSE ! NNCONFIG >= 1: CTV cases benefit from using the whole body

            ISTN1 = 1
            ISTN2 = NIF

C           Set up for starting from the trailing edge:

            IF (ZRG(1,KWING2) - ZRG(1,KWING1) >
     >          YRG(1,KWING2) - YRG(1,KWING1)) THEN ! Wing, not fin

               UINTR(1) = (YRG(1,KWING1) - YFMN) / (YFMX - YFMN)
               UINTR(1) = MAX (0.2, MIN (UINTR(1), 0.5)) ! Y is poor for LE u
               VINTR(1) = (XF(I2) - XF(1)) / (XF(NIF) - XF(1))
               TINTR(1) = (ZFMX          - ZRG(1,KWING1)) /
     >                    (ZRG(1,KWING2) - ZRG(1,KWING1))
               TINTR(1) = MAX (0.1, MIN (TINTR(1), 0.9))
               JINTR(1) = NJF * UINTR(1)
               KINTR(1) = (KWING1 + KWING2) / 2

            ELSE ! Assume tail fin is ~ vertical, on top

               UINTR(1) = ONE - 0.5 * (ZRG(1,KWING1) / ZFMX)
               UINTR(1) = MAX (0.6, MIN (UINTR(1), 0.9)) ! Y is poor for TE u
               VINTR(1) = (XF(I2) - XF(1)) / (XF(NIF) - XF(1))
               TINTR(1) = (YFMX          - YRG(1,KWING1)) /
     >                    (YRG(1,KWING2) - YRG(1,KWING1))
               TINTR(1) = MAX (0.1, MIN (TINTR(1), 0.5))
               JINTR(1) = NJF * UINTR(1)
               KINTR(1) = 1

            END IF

         END IF

      ELSE IF (LOSTINTR) THEN ! Restore leading edge results

         I1 = NHALF
         IINTR(I1) = IINTLE(IW) 
         JINTR(I1) = JINTLE(IW) 
         KINTR(I1) = KINTLE(IW) 
         TINTR(I1) = TINTLE(IW) 
         UINTR(I1) = UINTLE(IW) 
         VINTR(I1) = VINTLE(IW) 
         FIRST = .TRUE.

      END IF


      DO I = ISTN1, ISTN2
         XBODY(1:NJF,I) = XF(I)
      END DO

      UBODY(1,ISTN1) = -ONE ! Tells INTSEC4/5 to parameterize the body

      IF (NNCONFIG == 0) THEN ! Pre-CTV cases start from the LE
         I1 = NHALF ! Do the intersection upper surface first
         I2 = ILWING
      ELSE ! CTV cases start at TE
         I1 = 1
         I2 = ILWING
      END IF

      DO ! One- or two-pass loop

         CALL WBINTR (MXIWING, MXKWING, I1, I2, KWING1, KWING2,
     >                XRG, YRG, ZRG, MXJBODY, MXIBODY, 1, NJF,
     >                ISTN1, ISTN2, XBODY, YF, ZF, 1, IDEGFUS,
     >                UBODY, VBODY,
     >                XINTR, YINTR, ZINTR, UINTR, VINTR, TINTR,
     >                JINTR, IINTR, KINTR, FIRST, IWRIT, IFAIL, FAIL)

         IF (FAIL) GO TO 999

         IF (I1 /= 1 .AND. I2 /= 1) THEN ! Do the lower surface from the LE
            I2 = 1
         ELSE
            EXIT
         END IF

      END DO ! Two-pass loop

      IF (LOSTINTR) THEN ! Save leading edge results for later reuse
         IINTLE(IW) = IINTR(I1)
         JINTLE(IW) = JINTR(I1)
         KINTLE(IW) = KINTR(I1)
         TINTLE(IW) = TINTR(I1)
         UINTLE(IW) = UINTR(I1)
         VINTLE(IW) = VINTR(I1)
      END IF

  999 RETURN

      END SUBROUTINE WFINTR

C***********************************************************************
C
      SUBROUTINE WINSERT (NCALL, IW, K1, K2, IDIM, KDIM, ILE, ITL, ITU,
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
C                    which are now inputs.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER NCALL                 ! I   -3 means write on first call
      INTEGER IW                    ! I   Surface #, for printing only
      INTEGER K1, K2                ! I/O Subset of wing sections to treat;
                                    !     K2 may be incremented on output
      INTEGER IDIM, KDIM            ! I   Dimensions of regularized arrays
      INTEGER ILE, ITL, ITU         ! I   Leading & trailing edge mesh indices
      REAL    XRG(IDIM,KDIM),       ! I/O Input regularized sections, possibly
     >        YRG(IDIM,KDIM),       !     returned with new section(s) inserted
     >        ZRG(IDIM,KDIM)        !     at some ZCRANK locations
      INTEGER NCRANK                ! I   No. of Z stations to capture in grid
      REAL    ZCRANK(NCRANK)        ! I/O Span stations to be captured by the
                                    !     eventual grid: replaced by ZRG(ILE,*)
                                    !     if within 0.001 * semispan, else a
                                    !     new section is inserted & K2 increased
      INTEGER KCRANK(NCRANK)        !   O Indices of cranks in output sections
      INTEGER IWRIT                 ! I   Unit # for printing crank locations

C     Local constants:

      REAL,   PARAMETER :: ONE = 1., ZERO = 0.

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

            ZC = ZCRANK(KC)

            DO K = K1, K2 - 1
               IF (ZC < ZRG(ILE,K+1)) EXIT
            END DO

            R = (ZC - ZRG(ILE,K)) / (ZRG(ILE,K+1) - ZRG(ILE,K))

            IF (R == ZERO) THEN

               KCRANK(KC) = K

            ELSE ! Insert a section

               KI = K + 1
               KCRANK(KC) = KI

               IF (K2 == KDIM) THEN
                  IF (IWRIT > 0) THEN
                     WRITE (IWRIT, '(/,1X,A,F11.4)')
     >                 'WINSERT: No room to insert a section at Z =', ZC
                  END IF
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

      IF (NCALL == -3) THEN
         IF (IWRIT > 0) THEN
            WRITE (IWRIT, '(/,A,I2,A,6F11.4)') ' Surface', IW,
     >         ' paneling will capture station(s)', ZCRANK
            WRITE (IWRIT,'(A,6I4)')
     >         ' Corresponding defining section(s):', KCRANK
         END IF
      END IF

      GO TO 999

  900 CALL SYNCH
      STOP

  999 RETURN

      END SUBROUTINE WINSERT

C***********************************************************************
C
      SUBROUTINE WNINTR (NCALL, KWING1, KWING2, IQUAD, MXIWING, MXKWING,
     >                   MXINTR, ILE, IL, XRG, YRG, ZRG, NIF, NJF,
     >                   MXIBODY, MXJBODY, XF, YF, ZF, IDEGREE, ITL,ITU,
     >                   XINTR, YINTR, ZINTR, UINTR, VINTR, TINTR,
     >                   IINTR, JINTR, KINTR, IWRIT, FAIL)
C
C     WNINTR calculates a pylon/nacelle intersection given a regular
C     geometry mesh for each.  It has virtually the same arguments as
C     WFINTR, except XF (for the regularized nacelle) has 2 subscripts,
C     not 1, and I, J now mean I, K for the nacelle surface.
C     WNINTR serves also for the under-wing nacelle/diverter case.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   NCALL,                ! -3 means first call for this case
     >   KWING1, KWING2,       ! Subset of pylon sections to treat
     >   IQUAD(2),             ! (1) = quadrant of 1st nacelle section;
     >                         ! (2) = quadrant near the intersection;
     >                         !       3, 6, 9 or 12 o'clock, looking aft
     >   MXIWING, MXKWING,     ! Dimensions of pylon section arrays
     >   MXINTR,               ! Dimensions of intersection arrays
     >   ILE, IL               ! Middle/last indices of pylon sections
      REAL, INTENT (IN), DIMENSION (MXIWING,MXKWING) ::
     >   XRG, YRG, ZRG         ! Regularized pylon sections
      INTEGER, INTENT (IN) ::
     >   NIF,                  ! Active axial & circumferential
     >   NJF,                  ! indices of regular nacelle sections
     >   MXIBODY,              ! I and J dimensions of nacelle
     >   MXJBODY               ! section arrays
      REAL, INTENT (IN), DIMENSION (MXIBODY,MXJBODY) ::
     >   XF, YF, ZF            ! Regularized nacelle sections
      INTEGER, INTENT (IN) ::
     >   IDEGREE               ! 3 means bicubics (INTSEC4);
                               ! 1 means bilinear (INTSEC5) (nacelle)
      INTEGER, INTENT (OUT) ::
     >   ITL, ITU              ! First and last indices of the inter-
                               ! section (which may hit end of nacelle)
      REAL, INTENT (OUT), DIMENSION (MXINTR) ::
     >   XINTR, YINTR, ZINTR,  ! Intersection coordinates
     >   UINTR,                ! Corresponding nacelle surface "u" and
     >   VINTR,                ! "v" values (circumferential & axial)
     >   TINTR                 ! ITL:ITU elements are returned
                               ! with "t" values at the intersection
      INTEGER, INTENT (OUT), DIMENSION (MXINTR) ::
     >   IINTR, JINTR, KINTR   ! Corresp. nacelle & pylon cell indices
      INTEGER, INTENT (IN) ::
     >   IWRIT                 ! For diagnostics
      LOGICAL, INTENT (OUT) ::
     >   FAIL                  ! Success flag

C     Local constants:

      REAL, PARAMETER :: ONE = 1.

C     Local variables:

      INTEGER
     >   I, I1, I2, IFAIL, ISTN1, ISTN2, JSTN1, JSTN2
      REAL, DIMENSION (MXIBODY,MXJBODY) ::
     >   UBODY, VBODY
      LOGICAL
     >   FIRST

C     Execution:

      FIRST = NCALL == -3
      ITL   = 1 ! If intersection doesn't go off the end ...
      ITU   = IL

C     Determine the appropriate half of the outer nacelle surface.
C     The case of IQUAD(1) = IQUAD(2) should be trapped at a higher level,
C     so there are only three other possibilities:

      I = MOD (4 + (IQUAD(1) - IQUAD(2)) / 3, 4)

      IF (I == 1) THEN
         JSTN1 = NJF / 2 + 2
         JSTN2 = NJF - 2
      ELSE IF (I == 2) THEN
         JSTN1 = NJF / 4 + 2
         JSTN2 = NJF - JSTN1
      ELSE ! I = 3
         JSTN1 = 3
         JSTN2 = NJF / 2 - 2
      END IF

      ISTN1 = NIF / 2 + 2
      ISTN2 = NIF
      UBODY(ISTN1,JSTN1) = -ONE ! Tells INTSEC4/5 to parameterize the nacelle

C     Set up for starting from the pylon leading edge, upper surface first:

      I1 = ILE
      I2 = IL
      IINTR(I1) = ISTN1 + 1
      JINTR(I1) = (JSTN1 + JSTN2) / 2
      KINTR(I1) = KWING2 - 1
      UINTR(I1) = 0.1
      VINTR(I1) = 0.5
      TINTR(I1) = 0.9

      DO ! Two-pass loop

         FAIL = .FALSE. ! Disable detection of an off-the-end intersection

         CALL WBINTR (MXIWING, MXKWING, I1, I2, KWING1, KWING2,
     >                XRG, YRG, ZRG, MXIBODY, MXJBODY, ISTN1, ISTN2,
     >                JSTN1, JSTN2, XF, YF, ZF, 1, IDEGREE,
     >                UBODY, VBODY,
     >                XINTR, YINTR, ZINTR, UINTR, VINTR, TINTR,
     >                IINTR, JINTR, KINTR, FIRST, IWRIT, IFAIL, FAIL)

C ***    The commented-out code illustrates use of WBINTR's input FAIL=T option,
C        abandoned here - SURFER may extend the nacelle XTRAP units aft instead.

         IF (FAIL) THEN
C ***       IF (IFAIL == I1) THEN
               GO TO 999 ! Fatal
C ***       ELSE IF (IFAIL < I1) THEN ! Presume end-of-nacelle encountered
C ***          ITL = IFAIL + 1
C ***          DO I = 1, IFAIL
C ***             TINTR(I) = ONE ! All that WSURF needs along with ITL, ITU
C ***          END DO
C ***       ELSE
C ***          ITU = IFAIL - 1
C ***          DO I = IFAIL, IL
C ***             TINTR(I) = ONE
C ***          END DO
C ***       END IF
         END IF

         IF (I1 < I2) THEN ! Do the lower surface from the leading edge
            I2 = 1
         ELSE
            EXIT
         END IF

      END DO ! Two-pass loop

  999 RETURN

      END SUBROUTINE WNINTR

C***********************************************************************
C
      SUBROUTINE WSURF (NCALL, IW, IN, NNSURF, MXIWING, MXKWING,
     >                  MXCRANK, MXINTR, NREGSEC, NHALF, ILWING,KLWING,
     >                  XRG, YRG, ZRG, XINTR, YINTR, ZINTR, TINTR,
     >                  KWING1, KWING2, NCRANK, KCRANK, MCRANK,
     >                  TNORM, TBASE, XBASE, YBASE, ZBASE, IWRIT)
C
C     WSURF generates a structured mesh on a wing-type surface, which
C     is in regularized form.  Any intersections with other surfaces
C     (fuselage and/or nacelle, and wing if this is a diverter) are input.
C     Any input crank stations are captured by the mesh.
C
C     09/27/97  DAS  Initial arc-length-based approach from SYN87-SB.
C     10/08/97   "   Generalized in SURFER to allow tip-mounted nacelles.
C     07/31/98   "   Avoided incrementing NCRANK to help the paneling.
C                    MXCRANK is needed in case NCRANK = 0.
C     02/09/99   "   Moved splitting at LE up a level as part of packing
C                    X/Y/ZBASE(*).  ITL(*) & ITU(*) aren't needed.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   NCALL,              ! -3 means first call
     >   IW,                 ! Wing surface #
     >   IN,                 ! Tip nacelle surface #, or 0 if none
     >   NNSURF,             ! Number of wings + nacelles
     >   MXIWING, MXKWING,   ! Dimensions of regularized section arrays
     >   MXCRANK,            ! Max. # cranks handled
     >   MXINTR,             ! Max. # intersection points handled
     >   NREGSEC,            ! Number of regularized sections
     >   NHALF,              ! Leading edge index of input/output stations
     >   ILWING, KLWING      ! Desired paneling dimensions

      REAL, INTENT (IN), DIMENSION (MXIWING,MXKWING) ::
     >   XRG, YRG, ZRG       ! Regularized defining sections

      REAL, INTENT (IN), DIMENSION (MXINTR,NNSURF) ::
     >   XINTR, YINTR,       ! Intersection data for all surfaces
     >   ZINTR, TINTR

      INTEGER, INTENT (IN) ::
     >   KWING1(NNSURF),     ! Ranges of stations used for all of the
     >   KWING2(NNSURF),     ! intersections calculations
     >   NCRANK              ! Number of crank stations for this surface

      INTEGER, INTENT (INOUT), DIMENSION (MXCRANK) ::
     >   KCRANK, MCRANK      ! Input and output indices of crank stations

      REAL, INTENT (IN) ::
     >   TNORM(KLWING)       ! Normalized spanwise distribution

      REAL, INTENT (INOUT) ::
     >   TBASE(KLWING)       ! Denormalized nominal spanwise distribution

      REAL, INTENT (OUT), DIMENSION (ILWING,KLWING) ::
     >   XBASE, YBASE, ZBASE ! Surface patch (both surfaces, wrap-around)

      INTEGER, INTENT (IN) ::
     >   IWRIT               ! For error messages

C     Local constants:

      REAL,      PARAMETER :: ZERO = 0.
      LOGICAL,   PARAMETER :: NEW = .TRUE.
      CHARACTER, PARAMETER :: LINEAR * 1 = 'L'

C     Local variables:

      INTEGER
     >   I, IER, J, K, K1, K2
      REAL
     >   TCRANK(MXCRANK), TEVAL(KLWING), TRG(NREGSEC),
     >   WORK(KLWING * 2), XYZG(NREGSEC,3), XYZM(KLWING,3),
     >   DT, T1, T2, TSCALE, TSHIFT, TTIP, TTOTAL

C     Execution:

C     Considerata:
C     Blending the crank Ts into the nominal spanwise distribution can
C     lead to different K indices at the cranks for different I indices,
C     so use a single base distribution containing any cranks.
C     Shape changes could also lead to different Ks for different shapes,
C     so use the base distribution from the initial shape.
C     The leading edge distribution may be muddied by slat deflections,
C     so pick a representative one around 1/4 chord.

      IF (NCALL == -3) THEN

         I = 4 * NHALF / 3

C        Unnormalized arc lengths between wing sections for this I:

         CALL CHORDSRF (MXIWING, MXKWING, I, 1, NREGSEC, XRG, YRG, ZRG,
     >                  .FALSE., TTOTAL, TRG)

C        Set up a denormalized nominal base distribution using this I.
C        A non-intersection at the low end is anticipated in SURFER's set-up.

         TBASE(1) = TRG(KWING1(IW)) + TINTR(I,IW) *
     >             (TRG(KWING2(IW)) - TRG(KWING1(IW))) ! Left intersection T

         IF (IN == 0) THEN ! Plain right tip
            TBASE(KLWING) = TTOTAL
         ELSE
            TBASE(KLWING) = TRG(KWING1(IN)) + TINTR(I,IN) *
     >             (TRG(KWING2(IN)) - TRG(KWING1(IN))) ! Right intersection T
         END IF

         TTOTAL = TBASE(KLWING) - TBASE(1)

         DO K = 2, KLWING - 1
            TBASE(K) = TBASE(1) + TTOTAL * TNORM(K)
         END DO

         IF (NCRANK > 0) THEN

C           Blend crank T(s) into the base T distribution:

            DO K = 1, NCRANK
               TCRANK(K) = TRG(KCRANK(K))
            END DO

            CALL SMOOTHX (TBASE, KLWING, TCRANK, NCRANK, 3, WORK, IER)

C           Match base grid Ks with geometry Ks at cranks:

            K1 = 1
            DO J = 1, NCRANK
               DO K = K1, KLWING
                  IF (TCRANK(J) == TBASE(K)) THEN
                     MCRANK(J) = K
                     K1 = K1 + 1
                     EXIT
                  END IF
               END DO
            END DO

         END IF

C        Appending the tip as a final crank helps paneling between cranks
C        even if the input NCRANK is 0.

         KCRANK(NCRANK+1) = NREGSEC
         MCRANK(NCRANK+1) = KLWING

      END IF

C     Arc-length-based spanwise gridding:

      DO I = 1, ILWING

C        Spanwise arcs for this I:

         CALL CHORDSRF (MXIWING, MXKWING, I, 1, NREGSEC, XRG, YRG, ZRG,
     >                  .FALSE., TTOTAL, TRG)

C        Transform the base distribution by panels:

         T1 = TRG(KWING1(IW)) + TINTR(I,IW) *
     >       (TRG(KWING2(IW)) - TRG(KWING1(IW)))
         TEVAL(1) = T1
         K1 = 1

         IF (IN == 0) THEN
            TTIP = TTOTAL
         ELSE
            TTIP = TRG(KWING1(IN)) + TINTR(I,IN) *
     >            (TRG(KWING2(IN)) - TRG(KWING1(IN)))
         END IF

         DO J = 1, NCRANK + 1

            IF (J <= NCRANK) THEN
               T2 = TRG(KCRANK(J))
            ELSE
               T2 = TTIP
            END IF

            K2 = MCRANK(J)
            DT = TBASE(K2) - TBASE(K1)
            TSCALE = (T2 - T1) / DT
            TSHIFT = (TBASE(K2) * T1 - TBASE(K1) * T2) / DT

            DO K = K1 + 1, K2
               TEVAL(K) = TBASE(K) * TSCALE + TSHIFT
            END DO

            K1 = K2
            T1 = T2
         END DO
         TEVAL(KLWING) = T2

C        Piecewise linear interpolation for spanwise grid:

         DO K = 1, NREGSEC
            XYZG(K,1) = XRG(I,K)
            XYZG(K,2) = YRG(I,K)
            XYZG(K,3) = ZRG(I,K)
         END DO

         CALL LCSFIT (NREGSEC, TRG, XYZG(1,1), NEW, LINEAR, KLWING,
     >                TEVAL, XYZM(1,1), WORK)

         CALL LCSFIT (NREGSEC, TRG, XYZG(1,2), NEW, LINEAR, KLWING,
     >                TEVAL, XYZM(1,2), WORK)

         CALL LCSFIT (NREGSEC, TRG, XYZG(1,3), NEW, LINEAR, KLWING,
     >                TEVAL, XYZM(1,3), WORK)

         DO K = 1, KLWING
            XBASE(I,K) = XYZM(K,1)
            YBASE(I,K) = XYZM(K,2)
            ZBASE(I,K) = XYZM(K,3)
         END DO

      END DO

C     If there's an intersection at either end, ensure exactness:

      IF (TINTR(1,IW) /= ZERO) THEN ! See set-up in SURFER
         DO I = 1, ILWING
            XBASE(I,1) = XINTR(I,IW)
            YBASE(I,1) = YINTR(I,IW)
            ZBASE(I,1) = ZINTR(I,IW)
         END DO
      END IF

      IF (IN /= 0) THEN
         DO I = 1, ILWING
            XBASE(I,KLWING) = XINTR(I,IN)
            YBASE(I,KLWING) = YINTR(I,IN)
            ZBASE(I,KLWING) = ZINTR(I,IN)
         END DO
      END IF

C     Split a plain wrap-around surface at the leading edge upon return.

      END SUBROUTINE WSURF

C***********************************************************************
C
      SUBROUTINE XBCHECK (LENXB, MXPATCH, IPATCH, LASTXB, NI, NJ,
     >                    IPIECE, IWRIT, MOSTXB)
C
C     Check for exceeding the room available in surface patch arrays.
C
C     02/09/99  DAS  Initial implementation.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   LENXB,              ! Length of X/Y/ZBASE(*)
     >   MXPATCH,            ! Maximum # patches expected
     >   IPATCH,             ! Number of the prospective patch
     >   LASTXB,             ! Last index of X/Y/ZBASE(*) used so far
     >   NI, NJ,             ! Dimensions of prospective patch
     >   IPIECE,             ! Component number
     >   IWRIT               ! IWRIT < 0 suppresses any diagnostic

      INTEGER, INTENT (INOUT) ::
     >   MOSTXB              ! Most space needed so far

C     Local variables:

      INTEGER
     >   NEEDED

C     Execution:

      NEEDED = LASTXB + NI * NJ
      MOSTXB = MAX (MOSTXB, NEEDED)

      IF (IPATCH > MXPATCH .OR. NEEDED > LENXB) THEN
         IF (IWRIT > 0) WRITE (IWRIT, '(/, A, /, (1X, A, I9))')
     >      'Patch limits exceeded:',
     >      'MXPATCH =', MXPATCH,
     >      'IPATCH  =', IPATCH,
     >      'LASTXB  =', LASTXB,
     >      'LENXB   =', LENXB,
     >      'NEEDED  =', NEEDED,
     >      'NI      =', NI,
     >      'NJ      =', NJ,
     >      'IPIECE  =', IPIECE
         CALL SYNCH
         STOP
      END IF

      END SUBROUTINE XBCHECK

C***********************************************************************
C
      SUBROUTINE XBCOPY (LENXB, XBASE, YBASE, ZBASE,
     >                   I1, NI, NJ, INC, I2, INC2)
C
C     Copy packed X/Y/ZBASE(*) data from one location to another.
C     Use XBCHECK first where appropriate to make sure there is room.
C
C     02/09/99  DAS  Initial implementation.
C     02/22/99   "   Reverse order can be required to avoid clobbering.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   LENXB               ! Length of X/Y/ZBASE(*)

      REAL, INTENT (INOUT), DIMENSION (LENXB) ::
     >   XBASE, YBASE, ZBASE

      INTEGER, INTENT (IN) ::
     >   I1,                 ! First element to copy
     >   NI,                 ! # elements to copy for each J
     >   NJ,                 ! # sets of NI elements to copy
     >   INC,                ! Distance between starts of each set to copy
     >   I2,                 ! Start of copied data
     >   INC2                ! Distance between copied sets of NI elements

C     Local variables:

      INTEGER
     >   I, IA, IB, IY, IZ, J, N1

C     Execution:

      N1 = NI - 1

      IF (I2 <= I1) THEN ! Copying ascending indices should be OK

         IA = I1
         IZ = I2

         DO J = 1, NJ
            IB = IA + N1
            IY = IZ
            DO I = IA, IB ! ~ third of code for XBASE(IA:IB) = XBASE(IY:IZ) style
               XBASE(IY) = XBASE(I)
               YBASE(IY) = YBASE(I)
               ZBASE(IY) = ZBASE(I)
               IY = IY + 1
            END DO
            IA = IA + INC
            IZ = IZ + INC2
         END DO

      ELSE ! Copy in reverse order

         IA = I1 + NJ * INC
         IZ = I2 + NJ * INC2 + N1

         DO J = 1, NJ
            IZ = IZ - INC2
            IA = IA - INC
            IB = IA + N1
            IY = IZ
            DO I = IB, IA, -1
               XBASE(IY) = XBASE(I)
               YBASE(IY) = YBASE(I)
               ZBASE(IY) = ZBASE(I)
               IY = IY - 1
            END DO
         END DO

      END IF

      END SUBROUTINE XBCOPY

C***********************************************************************
C
      SUBROUTINE XBCOPY2 (LENXB, XBASE, YBASE, ZBASE, IXB, IDIM,
     >                    L, I1, I2, K1, K2, M, IA, KA)
C
C     Copy (I1:I2,K1:K2) of patch L to (IA:IB,KA:KB) of patch M with no
C     overlap assumed.  IB and KB would be redundant arguments.
C
C     03/06/99  DAS  Adaptation of XBCOPY for simpler use in WDSURF.
C     03/08/99   "   Edges seem to be the only usage, so tailor it.
C
C***********************************************************************

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   LENXB               ! Length of X/Y/ZBASE(*)

      REAL, INTENT (INOUT), DIMENSION (LENXB) ::
     >   XBASE, YBASE, ZBASE

      INTEGER, INTENT (IN) ::
     >   IXB(*),             ! IXB(L) = First element of patch L
     >   IDIM(*),            ! I dimension of each patch
     >   L, I1, I2, K1, K2,
     >   M, IA, KA

C     Local variables:

      INTEGER
     >   I, K, IL, IM, IXL, IXM

C     Execution:

C     Address (B(i,j)) = Address (B(1,1)) + (j - 1)idim + i - 1

      IL = IXB(L) + (K1 - 1) * IDIM(L) + I1 - 1
      IM = IXB(M) + (KA - 1) * IDIM(M) + IA - 1

      IF (K1 == K2) THEN ! An edge in the I direction

         DO I = IL, IL + I2 - I1
            XBASE(IM) = XBASE(I)
            YBASE(IM) = YBASE(I)
            ZBASE(IM) = ZBASE(I)
            IM = IM + 1
         END DO

      ELSE IF (I1 == I2) THEN ! K-direction edge

         DO K = K1, K2
            XBASE(IM) = XBASE(IL)
            YBASE(IM) = YBASE(IL)
            ZBASE(IM) = ZBASE(IL)
            IL = IL + IDIM(L)
            IM = IM + IDIM(M)
         END DO

      ELSE ! More than just an edge; no overlap assumed

         DO K = K1, K2
            IXL = IL
            IXM = IM
            DO I = I1, I2
               XBASE(IXM) = XBASE(IXL)
               YBASE(IXM) = YBASE(IXL)
               ZBASE(IXM) = ZBASE(IXL)
               IXL = IXL + 1
               IXM = IXM + 1
            END DO
            IL = IL + IDIM(L)
            IM = IM + IDIM(M)
         END DO

      END IF

      END SUBROUTINE XBCOPY2
