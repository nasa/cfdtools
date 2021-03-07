C*******************************************************************************
C
      PROGRAM SURFACE_CUTS
C
C     Description:
C     ------------
C
C     SURFACE_CUTS is intended to turn a triangulated surface into a
C     structured surface.  It slices the triangles at the indicated
C     stations for constant X, Y, or Z (one cut direction only per run),
C     and regularizes the slices via spline interpolation.
C
C     The input surface should represent a single geometry component
C     (wing, fuselage, or tail surface, presumably).  Otherwise,
C     regularizing the cuts would normally not make a lot of sense.
C     Handling of a nacelle is not attempted here yet.  Maybe a kludge
C     for handling cylindrical coordinates could work some day.
C
C     A right HALF body is presently treated correctly, but not a left half.
C     In this case, an option to reflect the half as a whole is provided.
C
C     The special case of TWO line segments for a cut is handled on the
C     assumption that an intersecting component causes a gap that needs
C     to be bridged for the present purposes.  E.g., body cuts may
C     encounter a gap at a wing root (or two gaps if a whole body is
C     being cut - still two line segments).
C
C     The special case of THREE line segments is also handled (as for a
C     half-fuselage where an off-center tail overlaps a wing).
C
C     This version has the option to output results in HAVOC form.
C     Since HAVOC expects the defining sections for winglets and fins to
C     have the same orientation as for a wing (with a rotation applied
C     by HAVOC), an option to rotate the geometry here before cutting it
C     is provided.  The direction of the cuts should be chosen accordingly.
C     For instance, a winglet should use ROTATION = 90 degrees (or less)
C     before 'Y' cuts are performed (not 'Z' cuts).  The rotation and
C     shifts printed for use by HAVOC assume that HAVOC rotates anti-
C     clockwise (looking downstream) before it shifts.
C
C     XYZ convention (right-handed):
C     ------------------------------
C
C        X increases downstream (nose to tail/leading to trailing edge);
C        Y increases towards the tip of a right wing;
C        Z increases upwards.
C
C     Input triangulation formats:
C     ----------------------------
C
C        The input surface should be one of the following:
C
C        Multi-zone Tecplot surface triangulation
C        ----------------------------------------
C
C        VARIABLES = "X", "Y", "Z"[, "TEMP", "PRESS", ...]
C        ZONE N=96000, E=32000, F=FEPOINT, ET=TRIANGLE
C        0.000000  0.000000 0.000000 [2838.386719 51330.925781 ...]
C        0.000883 -0.007150 0.113643 [2838.386719 51330.925781 ...]
C        0.000883  0.000000 0.113868 [2838.386719 51330.925781 ...]
C        ::::::::::::::::
C        ::::::::::::::::
C        4.882953  0.000000 0.011285 [ 950.867676    16.506409 ...]
C        1 2 3
C        4 5 6
C        7 8 9
C        ::::::::
C        ::::::::
C        95995 95996 95997
C        95998 95999 96000
C        [ZONE ....
C        :::::::::: Optional further zone(s)
C        ::::::::::
C        13579 24688 24689]
C
C        Multi-zone format of program FAST (formatted or unformatted)
C        ------------------------------------------------------------
C
C           NZONES                 ! May be omitted from single-zone surfaces
C           NPOINTS1, NTRIANGLES1, 0, ..., NPOINTSNZ, NTRIANGLESNZ, 0 tets
C           X1,X2, ..., XNP,Y1, ..., YNP,Z1, ..., ZNP    ! Coords. for zone 1
C           IV11,IV12,IV13, ..., IVNT1,IVNT2,IVNT3       ! Vertex indices
C           IC1,IC2, ..., ICNT                           ! Triangle flags
C           [Repeat for next zone(s)]
C
C        AIRPLANE format (formatted only; fixed files names)
C        ---------------------------------------------------
C
C           'press.d' (x,y,z,p) 4-tuples         'conn.d' pointers into press.d
C
C           Mach number                          # triangles
C           # points                             i1 i2 i3  ! Vertices of tri. 1
C           x  y  z  p  ! The ps are ignored     :  :  :   !    "     "    "  2
C           :  :  :  :                           :  :  :   !    :     :    :  :
C
C        NASTRAN format (formatted only)
C        -------------------------------
C
C        Any number of surface patches defined by m (x,y,z) triples, 1 per line,
C        and n (i1,i2,i3) triples, 1 per line defining the triangles.  E.g.,:
C
C              10 !  Point ID  x            y           z
C           BODY1          1   88.837000    0.000000   -7.484000
C           BODY1          2   88.837000    0.000000   -7.637000
C           BODY1          3   88.837000    0.812000   -7.637000
C           BODY1          4   88.837000    1.625000   -7.637000
C           BODY1          5   88.837000    2.437000   -7.637000
C           BODY1          6   88.837000    3.250000   -7.637000
C           BODY1          7   88.837000    3.250000   -7.484000
C           BODY1          8   88.837000    2.437000   -7.484000
C           BODY1          9   88.837000    1.625000   -7.484000
C           BODY1         10   88.837000    0.812000   -7.484000
C               8 ! Triangle ID Subcode?   3 vertex pointers
C           TRIA1          1       1       1       2      10
C           TRIA1          2       1       6       7       5
C           TRIA1          3       1       4       5       8
C           TRIA1          4       1       8       5       7
C           TRIA1          5       1       8       9       4
C           TRIA1          6       1       4       9       3
C           TRIA1          7       1       9      10       3
C           TRIA1          8       1      10       2       3
C           <Repeat (x,y,z) and (i1,i2,i3) pairs of groups till EOF;
C            each pair of groups represents a patch or zone>
C
C           WARNING:  The point and triangle numbers start at 1 for each zone,
C                     and are no more than ordinal numbers within each zone.
C
C        STL format (formatted only, single zone only)
C        ---------------------------------------------
C
C        solid ABC
C           facet normal 1.  0.  0.
C             outer loop
C               vertex 50. 30. 10.
C               vertex 50. 30. -10.
C               vertex 50. 31. -5.
C             endloop
C           endfacet
C           facet normal ...
C             ::::::::::::
C           endfacet
C           ::::::::::::::
C           ::::::::::::::
C        endsolid ABC
C
C
C     Output format:
C     --------------
C
C        The sliced surface is output as a single patch in PLOT3D /mgrid form,
C        formatted or unformatted.  If a half body is reflected, it remains
C        one patch:
C
C           NGRID = 1
C           NI, NJ  1 = # pts/section, # sections, 1
C           All Xs for surface patch
C           All Ys  "     "      "
C           All Zs  "     "      "
C
C        Such an output is not feasible if any of the initial cuts finds
C        more than one line segment.  Thus the input is expected to be one
C        plain component of an aerospace vehicle.  For diagnostic purposes,
C        the points along each slice, prior to regularization, are written to
C        the file 'line.segments':
C
C           Slice #
C           Segment # (should not exceed 1)       # points on this segment
C           x   y   z
C           :   :   :
C
C        A third output, suitable for HAVOC input, is optional.  Use 'none'
C        as the file name in order to suppress it.  Otherwise a formatted
C        file is produced as follows.
C
C           HAVOC body definition:
C
C           50 ! # cuts
C           51 ! # points per cut
C           0.   0.   0.   ! X, Y, Z at right water line of first section
C           0.   0.   0.   ! Proceed clockwise looking downstream
C            :    :    :
C           1.   .04  .03  ! X, Y, Z at water line of last section
C           1.   .04  .035 !    :  ! The body length is scaled to 1.
C            :    :    :
C           -999           ! End of component
C           -99999         ! End of data
C
C           HAVOC right wing definition:
C
C          [10 ! # cuts]   ! Omitted in standard practice
C           40 ! # points per cut
C           0.   0.   0.   ! X, Y, Z at root leading edge
C           .003 0. .0006  ! Proceed along upper surface, back along lower surf.
C            :    :    :
C            :    :    :
C           .2   .04  .03
C           -999
C           -99999
C
C           HAVOC fin/right-winglet definition:
C
C          [10 ! # cuts]   ! Omitted in standard practice
C           40 ! # points per cut
C           0.   0.   0.   ! X, Y, Z at root leading edge
C           .002 0. .0002  ! Proceed along left surface, back along right surf.
C            :    :    :
C            :    :    :
C           .12  .01  .02
C           -999
C           -99999
C
C        When a body is cut for HAVOC purposes, the NEGATIVES of the shifts
C        used to move the nose to (0, 0, 0) are displayed.  These apply to the
C        outputs that have been normalized by REF_LENGTH:
C
C           HAVOC x/y/z = (full scale x/y/z) / REF_LENGTH + shift
C
C        This same REF_LENGTH should be entered when other related components
C        are cut.  The NEGATIVES of the shifts needed to move a root leading
C        edge to (0, 0, 0) are displayed in these cases.
C
C     Input control file:
C     -------------------
C
C        Control inputs are entered on the command line.        Sample:
C
C     SURFACE_CUTS control file for ...        ! Title
C     'Geometry.nas'          ! Input surface triangulation file name in quotes
C     T                       ! Formatted triangulation? T = yes; F = unform.
C     T                       ! T = right half body input; F = whole body input
C     T                       ! T = reflect if right half body -> whole body out
C     NASTRAN                 ! Triang. type: TECPLOT|FAST|AIRPLANE|NASTRAN|STL
C     2                       ! Number of active zones specified on next line(s)
C     1  2
C     T                       ! Formatted output? T = yes; F = unformatted
C     'wing-no-te.xyz'        ! Plottable output file name (PLOT3D /mgrid)
C     'none'                  ! HAVOC-type output file name ('none' = suppress)
C     1.                      ! REF_LENGTH, applied to HAVOC-type output (only)
C     0.                      ! Rotation ab. OX (deg. anticl. looking downstrm.)
C     Y                       ! Type of cuts:  X|Y|Z for constant X|Y|Z
C     T                       ! T = clockwise direction of PLOT3D pts along cuts
C     51                      ! NU (# pts. output on upper wing or above chine)
C     FOILGRD                 ! UNIFORM | FOILGRD | VINOKUR
C     0.04  0.  0.3  0.66     ! FOILGRD weights for linear/quad/sin/cos terms
C     0.007  0.090            ! VINOKUR first and last increments on [0., 1.]
C     51                      ! NL (# pts. output on lower wing or below chine)
C     FOILGRD                 ! UNIFORM | FOILGRD | VINOKUR
C     0.04  0.  0.3  0.66     ! FOILGRD weights for linear/quad/sin/cos terms
C     0.007  0.090            ! VINOKUR first and last increments on [0., 1.]
C     0.    0.    0.          ! X/Y/Znose; Xcut(1) = Xnose gives a singular pt.
C     999.  999.  999.        ! X/Y/Ztail; Xcut(Ncuts) = Xtail " " " " " " " "
C     0                       ! # refined cuts calculated here (0 or >= 10)
C     27                      ! # cuts defined on following lines (1 per line)
C     1  4.3                  ! cut # and station (before scaling)
C     2  5.                   !  "  "  "  "
C     3  6.
C     4  7.
C     :   :
C     :   :
C
C     Further usage notes:
C     --------------------
C
C        > Cut stations (for fuselages in particular) are tedious to prepare.
C          It is suggested that the preliminary cuts be used to produce a
C          refined set of size specified on the line above the current number
C          of cuts.  The rate of change of sectional area is used to compute
C          the refined set of cuts.
C
C        > REF_LENGTH applies to HAVOC-type outputs only.  The cutting is done
C          in the units of the input triangulation, and the PLOT3D-type struc-
C          tured surface is in the original coordinate system and units.
C
C        > For a body being normalized for HAVOC purposes, the nose is always
C          moved to (0, 0, 0).  The nose point here should be the singular
C          point specified in the input control file, with its X matching
C          the first X cut station.
C
C        > Any singular point coordinates for a nose or tail should be in the
C          system of the input geometry (regardless of REF_LENGTH and whether
C          or not HAVOC-type results are requested).
C
C        > The CLOCKWISE input is ignored if HAVOC-type outputs are requested,
C          because these outputs are written assuming the internal order of
C          anticlockwise for a body and clockwise for other components.
C
C        > For wings or tails, the total number of points per output section
C          is NINTERP = NU + NL - 1.  However, the optional HAVOC-type output
C          has NU + NL points because the trailing edge is assumed to be finite.
C
C        > For fuselages, the total number of points per output section is
C          2 * (NU + NL) - 3, even if uniform distributions are specified,
C          except for the case of a right half body where NINTERP = NU + NL - 1.
C
C        > Three types of point distribution are provided:
C             UNIFORM is often appropriate for fuselages;
C             FOILGRD combines specified multiples of linear, quadratic, sine,
C                     and cosine terms;
C             VINOKUR iterates to achieve precise initial and final increments.
C
C        > Choice of nonuniform point distributions depends on the application.
C          Possible FOILGRD coefficients for normalized distributions:
C
C          WTLIN     WTQUAD    WTSIN     WTCOS     Grid Type
C
C          0.04      0.0       0.3       0.66      fine, bunch at both ends
C          0.0783536 0.4210620 0.5005844 0.0000000 coarse, steady increase in dX
C          0.5       0.0       0.0       0.5       for body with chine?
C          1.0       0.0       0.0       0.0       uniform
C
C        > Vinokur distributions may be preferable for wing/tail cuts.  The
C          specified first and last increment are relative to the local chord.
C
C        > If UNIFORM distributions are specified for both upper and lower
C          portions of the cuts, then the entire cut is redistributed as one
C          segment.  This should be preferable for bodies with flat sides or
C          under-surfaces, where identification of the water line and keel
C          indices is awkward and prone to erratic results.
C
C        > [For treating solid bodies:]  If the HAVOC output file name
C          contains 'volume', output instead an interior volume grid formed
C          by discretizing cross-section radii formed by joining each
C          regularized cut point to the nose-tail point for that cut.
C          Enter # refined cuts = -33 to get 33 points along each radius.
C
C     History:
C     --------
C
C        12/22/01  D.Saunders  Initial adaptation of MULTICUT, which treats
C                              a structured surface by turning each quad.
C                              into two triangles and using the unstructured
C                              techniques of Scott Thomas.
C        12/26/01       "      Allowed for FAST or AIRPLANE triangulations.
C        01/07/02       "      Switched from reading all zones and possibly
C                              blanking some to reading a list of zones to cut,
C                              and reading only those, skipping the rest.
C        01/08/02       "      Handled NASTRAN inputs via scheme of 01/07/02.
C        01/09/02       "      Provided likely nonuniform point distributions
C                              in pieces, retaining one-piece uniform option.
C        01/10/02       "      Handled right half bodies (only); whole
C                              components are assumed otherwise; option to
C                              reflect a half as a whole; option for singular
C                              point at nose and/or tail.
C        01/14/02       "      Provided HAVOC-type output options, and the
C                              REF_LENGTH input.
C        01/16/02       "      Handled two-segment cuts specially, as needed
C                              for gaps in bodies at a wing root.
C        01/22/02       "      Efforts to optimize FOILGRD coefficients were
C                              only partially successful.  (The initial growth
C                              rates of the increments tend to be too big.)
C                              Therefore, Vinokur distributions are provided
C                              as an alternative.
C        02/14/02       "      A volume estimate is calculated from the current
C                              cuts along with a suggested set of refined cut
C                              stations based on the rate of change of the
C                              sectional area (later, on local curvature, but
C                              it is still unreliable).
C        02/21/02       "      Winglets for HAVOC need to be rotated flat before
C                              cutting.  This affects the shifts displayed.
C        02/23/02       "      Extended the 2-line-segment case to include the
C                              3-line-segment case, but full generalization
C                              appears difficult.
C        07/19/02       "      Handled single-zone STL format (formatted only).
C        10/23/14       "      Revived to deal with asteroids, requiring reading
C                              of Tecplot surface triangulations (1+ zones).
C        11/04/14       "      Option to generate an interior volume grid from
C                              the regularized cuts as described in the last
C                              bullet above (-33 -> 33 points along a radius).
C                              The volume grid is intended for CM and moments
C                              of inertia calculations.
C        11/05/14       "      Testing on a unit-radius sphere produce some
C                              1-point extra line segments at 2 of 129 cut
C                              stations that can't be explained by EPS.
C                              Workaround:  if VOLUME_GRID, force NLINES = 1.
C
C     Author:  David Saunders, ELORET/NASA Ames, Moffett Field, CA
C              Now with ERC, Inc. at NASA/ARC.
C
C*******************************************************************************

!     Modules:

      USE TRI_HEADER_STRUCTURE  ! Part of triangulation_io.f90
      USE TRI_ZONE_STRUCTURE    ! Likewise
      USE TRIANGULATION_IO      ! Unstructured data file I/O
      USE TRIGD

      IMPLICIT NONE

C     Constants:

      INTEGER, PARAMETER ::
     >   NADDIM  = 100,     ! 1 + max. # line strings per cut over all zones
     >   LUNDAT  = 1,       ! Input triangulation in one of several formats
     >   LUNCUT  = 2,       ! Output file of cuts in PLOT3D /mgrid format
     >   LUNSEG  = 3,       ! Output file 'line.segments' before regularization
     >   LUNHVC  = 4,       ! HAVOC-type output file | interior vol.grid | none
     >   LUNCTL  = 5,       ! Control inputs (on standard input)
     >   LUNCRT  = 6,       ! Screen
     >   MAXSEG  = 3000,    ! Max. # 2-point line segments handled per cut
     >   MAXNODE = 2*MAXSEG ! Max. # end-points of line segments handled per cut

      REAL, PARAMETER ::
     >   EPS     = 1.E-6,   ! Fraction of data range for CUT_TRIANGLES
     >   ONE     = 1.,
     >   ZERO    = 0.

C     Variables:

      INTEGER ::
     >   I, IDUMMY, IER, IKL, ILE, IOS, IP1, IP2, IT1, IT2, IWL, IWR,
     >   IXYZ, IZ, J, K, L, LAST_ZONE, N, NCUTS, NINTERP, NL, NPSKIP,
     >   NPTOTAL, NPTS, NREFINED, NTETS, NTSKIP, NTTOTAL, NU, NZONES,
     >   NZONE_ACTIVE, NK,
     >   NAD(NADDIM),
     >   NDGCUT(2,MAXSEG)  ! Pointers used by CUT_TRIANGLES & CUTORDER

      INTEGER, ALLOCATABLE ::
     >   IC(:), IP(:), IT(:), IV(:,:), IZONE_ACTIVE(:),
     >   NPOINTS_FOUND(:), NTRIANGLES_FOUND(:),
     >   NPOINTS(:), NTRIANGLES(:)

      REAL ::
     >   COSTHETA, DELTA1, DELTA2, DUMMY, EPSXYZ, GSCALE,
     >   REF_LENGTH, ROTATION, SINTHETA, VOLUME,
     >   WTLIN, WTQUAD, WTSIN, WTCOS,
     >   XMAX, XMIN, XNOSE, XTAIL,
     >   YMAX, YMIN, YNOSE, YTAIL, YTEMP,
     >   ZMAX, ZMIN, ZNOSE, ZTAIL, ZTEMP,
     >   P(3), S(3), SEG(4,MAXNODE), FCRV(MAXSEG), SCRV(MAXSEG),
     >   XCRV(MAXSEG), YCRV(MAXSEG), ZCRV(MAXSEG)

      REAL, ALLOCATABLE, DIMENSION (:) ::
     >   AREAS, CUTS, DERIVS, F, REFINED_CUTS,
     >   SNORM, SNORML, SNORML_REVERSED, SNORMU, SNORMU_REVERSED

      REAL, ALLOCATABLE, DIMENSION (:,:) ::
     >   X, Y, Z, XYZ

      REAL, ALLOCATABLE, DIMENSION (:,:,:) ::
     >   XVOL, YVOL, ZVOL

      LOGICAL ::
     >   CLOCKWISE, CONSTANT, FORMATTED, HALFBODY, HAVOC,
     >   REFLECT, UNIFORM, VOLUME_GRID

      CHARACTER ::
     >   DISTRIB*7, FILENAME*64, FORM*11, HAVOC_FILENAME*64, LABEL*6,
     >   TEXT*1, TYPE*8, XYZCASE*1

C     Derived data types:

      TYPE (TRI_HEADER_TYPE) ::
     >   TRI_HEADER
      TYPE (TRI_TYPE), POINTER, DIMENSION (:) ::
     >   TRI_XYZF

C     Execution:

CCC   OPEN (LUNCTL, FILE='surface_cuts.inp', STATUS='OLD')  ! Standard input now

C     ******************** Start of reading control inputs *********************

C     Read controls as needed, as opposed to reading them all ahead of time.

      READ (LUNCTL, '(A)')           ! Title - not used
      READ (LUNCTL, *) FILENAME      ! Surface triangulation name
      READ (LUNCTL, *) FORMATTED     ! T or F
      READ (LUNCTL, *) HALFBODY      ! T means RIGHT half
      READ (LUNCTL, *) REFLECT       ! T means reflect right half to give whole
      READ (LUNCTL, '(A)') TYPE      ! FAST|AIRPLANE|NASTRAN|STL (upper case)
      READ (LUNCTL, *) NZONE_ACTIVE  ! # active zones to read on next line(s)

      ALLOCATE (IZONE_ACTIVE(NZONE_ACTIVE))

      READ (LUNCTL, *) IZONE_ACTIVE

      LAST_ZONE = IZONE_ACTIVE(NZONE_ACTIVE)

      ALLOCATE (NPOINTS(NZONE_ACTIVE), NTRIANGLES(NZONE_ACTIVE),
     >          IP(NZONE_ACTIVE), IT(NZONE_ACTIVE))

      IF (FORMATTED) THEN
         FORM = 'FORMATTED'
      ELSE
         FORM = 'UNFORMATTED'
      END IF

C     **************** Start of reading a surface triangulation  ***************

      SELECT CASE (TYPE)

      CASE ('TECPLOT ') !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         TRI_HEADER%FILENAME  = FILENAME
         TRI_HEADER%FILEFORM  = 1       ! Vertex-centered fn. data for now
         TRI_HEADER%FORMATTED = .TRUE.  ! Can't read binary Tecplot files
         TRI_HEADER%NVERTICES = 3       ! Triangles, not tetrahedra

         IOS = 1                        ! Verbose mode

         CALL TRI_READ (LUNDAT, TRI_HEADER, TRI_XYZF, IOS)

         IF (IOS /= 0) THEN
            WRITE (LUNCRT, '(/, A)')
     >         'Trouble reading input triangulated dataset: aborting.'
            GO TO 999
         END IF

C        For this retrofit, follow the FAST format approach by transferring
C        the active zones to the work arrays originally employed, and
C        deallocating the derived data type arrays.

         NZONES = TRI_HEADER%NZONES

         ALLOCATE (NPOINTS_FOUND(NZONES), NTRIANGLES_FOUND(NZONES))

         NPOINTS_FOUND(:)    = TRI_XYZF(:)%NNODES
         NTRIANGLES_FOUND(:) = TRI_XYZF(:)%NELEMENTS

         NPTOTAL = 0
         NTTOTAL = 0

         DO L = 1, NZONE_ACTIVE
            N = IZONE_ACTIVE(L)
            NPOINTS(L)    = NPOINTS_FOUND(N)
            NTRIANGLES(L) = NTRIANGLES_FOUND(N)
            IP(L)   = NPTOTAL + 1 ! 1st element in XYZ(1:3,*), F(*) for zone L
            NPTOTAL = NPTOTAL + NPOINTS(L)
            IT(L)   = NTTOTAL + 1 ! 1st element in IV(1:3,*), IC(*) for zone L
            NTTOTAL = NTTOTAL + NTRIANGLES(L)

            WRITE (LUNCRT, '(I7, I10, I11, I14)')
     >         N, L, NPOINTS(L), NTRIANGLES(L)
         END DO

         ALLOCATE (XYZ(3,NPTOTAL), F(NPTOTAL),
     >             IV(3,NTTOTAL), IC(NTTOTAL))

C        Transfer the active coordinates and pointers:

         DO L = 1, NZONE_ACTIVE
            N = IZONE_ACTIVE(L)
            IP1 = IP(L)
            IP2 = IP1 + NPOINTS(L) - 1
            IT1 = IT(L)
            IT2 = IT1 + NTRIANGLES(L) - 1
            XYZ(:,IP1:IP2) = TRI_XYZF(N)%XYZ( :,:)
            IV( :,IT1:IT2) = TRI_XYZF(N)%CONN(:,:)
         END DO
         F(:) = ZERO  ! CUT_TRIANGLES looks at a scalar fn., so it must be set

         CALL DEALLOCATE_TRI_ZONES (1, NZONES, TRI_HEADER%NUMF, TRI_XYZF, IOS)

         IF (IOS /= 0) GO TO 999

      CASE ('FAST    ') !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         OPEN (LUNDAT, FILE=FILENAME, STATUS='OLD', FORM=FORM)

C        The single-zone case may not have NZONES present.  Try for 3 integers:

         IF (FORMATTED) THEN
            READ (LUNDAT, *, IOSTAT=IOS) IP1, IP2, NTETS
         ELSE
            READ (LUNDAT,    IOSTAT=IOS) IP1, IP2, NTETS
         END IF

         REWIND (LUNDAT)

         IF (IOS == 0) THEN ! Single zone
            NZONES = 1
         ELSE
            IF (FORMATTED) THEN
               READ (LUNDAT, *) NZONES
            ELSE
               READ (LUNDAT) NZONES
            END IF
         END IF

         IF (NZONES < LAST_ZONE) GO TO 900

         ALLOCATE (NPOINTS_FOUND(NZONES), NTRIANGLES_FOUND(NZONES))

         IF (FORMATTED) THEN
            READ (LUNDAT, *)
     >         (NPOINTS_FOUND(N), NTRIANGLES_FOUND(N), NTETS,
     >          N = 1, NZONES)
         ELSE
            READ (LUNDAT)
     >         (NPOINTS_FOUND(N), NTRIANGLES_FOUND(N), NTETS,
     >          N = 1, NZONES)
         END IF

         WRITE (LUNCRT, '(/, A)')
     >      ' Zone #  Active #   # points   # triangles'

         NPTOTAL = 0
         NTTOTAL = 0

         DO L = 1, NZONE_ACTIVE
            N = IZONE_ACTIVE(L)
            NPOINTS(L)    = NPOINTS_FOUND(N)
            NTRIANGLES(L) = NTRIANGLES_FOUND(N)
            IP(L)   = NPTOTAL + 1 ! 1st element in XYZ(1:3,*), F(*) for zone L
            NPTOTAL = NPTOTAL + NPOINTS(L)
            IT(L)   = NTTOTAL + 1 ! 1st element in IV(1:3,*), IC(*) for zone L
            NTTOTAL = NTTOTAL + NTRIANGLES(L)

            WRITE (LUNCRT, '(I7, I10, I11, I14)')
     >         N, L, NPOINTS(L), NTRIANGLES(L)
         END DO

         ALLOCATE (XYZ(3,NPTOTAL), F(NPTOTAL),
     >             IV(3,NTTOTAL), IC(NTTOTAL))

         L = 1    ! Read only the active zones:

         DO N = 1, NZONES

            IF (N < IZONE_ACTIVE(L)) THEN ! Skip this zone
               NPSKIP = 3 * NPOINTS_FOUND(N)
               NTSKIP = 4 * NTRIANGLES_FOUND(N)

               IF (FORMATTED) THEN
                  READ (LUNDAT, *) (DUMMY,  I = 1, NPSKIP),
     >                             (IDUMMY, I = 1, NTSKIP)
               ELSE
                  READ (LUNDAT) ! Skip a record
               END IF
            ELSE
               IP1 = IP(L)
               IP2 = IP1 + NPOINTS(L) - 1
               IT1 = IT(L)
               IT2 = IT1 + NTRIANGLES(L) - 1
               L   = L + 1

               IF (FORMATTED) THEN
                  READ (LUNDAT, *)
     >               ((XYZ(J,I), I = IP1, IP2), J = 1, 3),
     >               (IV(1:3,I), I = IT1, IT2), IC(IT1:IT2)
               ELSE
                  READ (LUNDAT)
     >               ((XYZ(J,I), I = IP1, IP2), J = 1, 3),
     >               (IV(1:3,I), I = IT1, IT2), IC(IT1:IT2)
               END IF
            END IF

         END DO

      CASE ('AIRPLANE') ! Single zone !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         OPEN (LUNDAT, FILE='press.d', STATUS='OLD', FORM=FORM)

         READ (LUNDAT, *) ! Mach number
         READ (LUNDAT, *) NPTOTAL

         NPOINTS(1) = NPTOTAL
         IP(1) = 1 !  1st element in XYZ(1:3,*) & F(*) for zone 1
         IT(1) = 1 !  1st element in IV(1:3,*) & IC(*) for zone 1

         ALLOCATE (XYZ(3,NPTOTAL), F(NPTOTAL))

         DO I = 1, NPTOTAL
            READ (LUNDAT, *) XYZ(1:3,I), F(I)
         END DO

         CLOSE (LUNDAT)

         OPEN (LUNDAT, FILE='conn.d', STATUS='OLD', FORM=FORM)

         READ (LUNDAT, *) NTTOTAL
         NTRIANGLES(1)  = NTTOTAL

         ALLOCATE (IV(3,NTTOTAL), IC(NTTOTAL))

         IT1 = IT(1)
         IC  = 1    ! Identifier for all triangles (not used)

         READ (LUNDAT, *) (IV(1:3,I), I = 1, NTTOTAL)

      CASE ('NASTRAN ') !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         OPEN (LUNDAT, FILE=FILENAME, STATUS='OLD', FORM=FORM)

C        Scan the data to determine the storage needed:

         NPTOTAL = 0
         NTTOTAL = 0
         L = 1
         N = 0

         WRITE (LUNCRT, '(/, A)')
     >      ' Zone #  Active #   # points   # triangles'

         READ (LUNDAT, *, IOSTAT=IOS) NPSKIP

         DO WHILE (IOS == 0 .AND. N < LAST_ZONE)

            N = N + 1
            DO I = 1, NPSKIP
               READ (LUNDAT, *) ! label, point #, x, y, z
            END DO

            READ (LUNDAT, *) NTSKIP
            DO I = 1, NTSKIP
               READ (LUNDAT, *) ! label, triangle #, subcode, i1, i2, i3
            END DO

            IF (N == IZONE_ACTIVE(L)) THEN
               NPOINTS(L)    = NPSKIP
               NTRIANGLES(L) = NTSKIP
               IP(L)   = NPTOTAL + 1 ! 1st elem. in XYZ(1:3,*), F(*) for zone L
               NPTOTAL = NPTOTAL + NPSKIP
               IT(L)   = NTTOTAL + 1 ! 1st elem. in IV(1:3,*), IC(*) for zone L
               NTTOTAL = NTTOTAL + NTSKIP

               WRITE (LUNCRT, '(I7, I10, I11, I14)')
     >            N, L, NPSKIP, NTSKIP

               L = L + 1
            END IF

            READ (LUNDAT, *, IOSTAT=IOS) NPSKIP
         END DO

         NZONES = N

         IF (NZONES /= LAST_ZONE) GO TO 900

         ALLOCATE (XYZ(3,NPTOTAL), F(NPTOTAL),
     >             IV(3,NTTOTAL), IC(NTTOTAL))

C        Read the specified zones (only):

         REWIND (LUNDAT)
         L = 1

         DO N = 1, NZONES

            READ (LUNDAT, *) NPSKIP

            IF (N < IZONE_ACTIVE(L)) THEN ! Skip this zone

               DO I = 1, NPSKIP
                  READ (LUNDAT, *)
               END DO

               READ (LUNDAT, *) NTSKIP
               DO I = 1, NTSKIP
                  READ (LUNDAT, *)
               END DO

            ELSE

               K = IP(L)
               DO I = 1, NPSKIP
                  READ (LUNDAT, *) LABEL, J, XYZ(1:3,K)
                  K = K + 1
               END DO

               READ (LUNDAT, *) NTSKIP
               K = IT(L)
               DO I = 1, NTSKIP
                  READ (LUNDAT, *) LABEL, J, J, IV(1:3,K)
                  IC(K) = N
                  K = K + 1
               END DO

               L = L + 1
            END IF

         END DO ! Next zone in input file

      CASE ('STL     ') ! Single zone !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         OPEN (LUNDAT, FILE=FILENAME, STATUS='OLD', FORM=FORM)
         READ (LUNDAT, *) ! Title

         NTTOTAL = 0

         DO ! Until EOF, count # points and # triangles

            READ (LUNDAT, *) ! 'facet' or 'endsolid'
            READ (LUNDAT, *, IOSTAT=IOS) ! 'outer' or EOF
            IF (IOS < 0) EXIT

            NTTOTAL = NTTOTAL + 1
            READ (LUNDAT, *)
            READ (LUNDAT, *)
            READ (LUNDAT, *)
            READ (LUNDAT, *)
            READ (LUNDAT, *)

         END DO

         NPTOTAL       = NTTOTAL * 3
         NPOINTS(1)    = NPTOTAL
         NTRIANGLES(1) = NTTOTAL

         IP(1) = 1 !  1st element in XYZ(1:3,*) & F(*) for zone 1
         IT(1) = 1 !  1st element in IV(1:3,*) & IC(*) for zone 1

         WRITE (LUNCRT, '(/, A, I7, A, I7)')
     >      ' # points:', NPTOTAL, '    # triangles:', NTTOTAL

         ALLOCATE (XYZ(3,NPTOTAL), F(NPTOTAL),
     >             IV(3,NTTOTAL), IC(NTTOTAL))

         REWIND (LUNDAT)
         READ   (LUNDAT, *) ! Title
         NPTS = -2

         DO N = 1, NTTOTAL
            NPTS    = NPTS + 3
            IV(1,N) = NPTS
            IV(2,N) = NPTS + 1
            IV(3,N) = NPTS + 2
            IC(N)   = 1

            READ (LUNDAT, *) ! 'facet'; ignore unit normal information
            READ (LUNDAT, *) ! 'outer'
            READ (LUNDAT, *) TEXT, XYZ(1:3,NPTS)
            READ (LUNDAT, *) TEXT, XYZ(1:3,NPTS+1)
            READ (LUNDAT, *) TEXT, XYZ(1:3,NPTS+2)
            READ (LUNDAT, *) ! 'endloop'
            READ (LUNDAT, *) ! 'endfacet'
         END DO

C        Save the triangulation in FAST format for visual inspection,
C        unless it already exists:

         CLOSE (LUNDAT)
         OPEN  (LUNDAT, FILE='FAST.xyz', STATUS='NEW', IOSTAT=IOS)

         IF (IOS == 0) THEN

            WRITE (LUNDAT, '(3I10)') NPTOTAL, NTTOTAL, 0
            DO I = 1, 3
               WRITE (LUNDAT, '(8ES15.7)')
     >            (XYZ(I,IV(1,J)), XYZ(I,IV(2,J)), XYZ(I,IV(3,J)),
     >            J = 1, NTTOTAL)
            END DO
            WRITE (LUNDAT, '(15I8)') IV
            WRITE (LUNDAT, '(60I2)') IC

            WRITE (LUNCRT, '(/, A)')
     >         ' Triangulation saved as "FAST.xyz".'
         END IF

      CASE DEFAULT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         WRITE (LUNCRT, '(/, 2A)')
     >      ' Bad input for triangulation type: ', TYPE
         GO TO 999

      END SELECT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      CLOSE (LUNDAT)

C     ***************** End of reading a surface triangulation  ****************

      READ (LUNCTL, *) FORMATTED        ! Output file type
      IF (FORMATTED) THEN
         FORM = 'FORMATTED'
      ELSE
         FORM = 'UNFORMATTED'
      END IF

      READ (LUNCTL, *) FILENAME         ! Output structured surface name

      OPEN (LUNCUT, FILE=FILENAME, STATUS='UNKNOWN', FORM=FORM)
      OPEN (LUNSEG, FILE='line.segments', STATUS='UNKNOWN')

      READ (LUNCTL, *) HAVOC_FILENAME   ! HAVOC-type output file name

      VOLUME_GRID = INDEX (HAVOC_FILENAME, 'volume') > 0
      HAVOC = HAVOC_FILENAME(1:4) /= 'none' .AND. .NOT. VOLUME_GRID

      READ (LUNCTL, *) REF_LENGTH       ! Scale factor applied to HAVOC outputs
      GSCALE = ONE / REF_LENGTH

      READ (LUNCTL, *) ROTATION         ! Degrees about X axis,
                                        ! clockwise looking downstream
      READ (LUNCTL, *)    XYZCASE
      IF (XYZCASE == 'x') XYZCASE = 'X' ! Type of cut; ensure upper case
      IF (XYZCASE == 'y') XYZCASE = 'Y'
      IF (XYZCASE == 'z') XYZCASE = 'Z'

      READ (LUNCTL, *) CLOCKWISE ! Order of PLOT3D points looking along cut axis

      READ (LUNCTL, *) NU        ! # output pts. on upper wing/above water line
      READ (LUNCTL, '(A)') DISTRIB ! Upper case please
      READ (LUNCTL, *) WTLIN, WTQUAD, WTSIN, WTCOS
      READ (LUNCTL, *) DELTA1, DELTA2

      UNIFORM  = DISTRIB == 'UNIFORM' .OR. WTLIN == ONE
      CONSTANT = UNIFORM

      IF (UNIFORM) THEN ! Make sure of it
         WTLIN  = ONE
         WTQUAD = ZERO
         WTSIN  = ZERO
         WTCOS  = ZERO
      END IF

      ALLOCATE (SNORMU(NU), SNORMU_REVERSED(NU))

      IF (UNIFORM .OR. DISTRIB == 'FOILGRD') THEN

         CALL FOILGRD (NU, ZERO, ONE, WTLIN, WTQUAD, WTSIN, WTCOS,
     >                 SNORMU)

      ELSE IF (DISTRIB == 'VINOKUR') THEN

         SNORMU(1)  = ZERO
         SNORMU(NU) = ONE

         CALL VINOKUR (1, NU, DELTA1, DELTA2, SNORMU, LUNCRT, IER)
         IF (IER /= 0) GO TO 910

      ELSE
         GO TO 920
      END IF

      READ (LUNCTL, *) NL        ! # output pts. on lower wing/below water line
      READ (LUNCTL, '(A)') DISTRIB
      READ (LUNCTL, *) WTLIN, WTQUAD, WTSIN, WTCOS
      READ (LUNCTL, *) DELTA1, DELTA2

      UNIFORM  = DISTRIB == 'UNIFORM' .OR. WTLIN == ONE
      CONSTANT = CONSTANT .AND. UNIFORM

      IF (UNIFORM) THEN ! Make sure of it
         WTLIN  = ONE
         WTQUAD = ZERO
         WTSIN  = ZERO
         WTCOS  = ZERO
      END IF

      ALLOCATE (SNORML(NL), SNORML_REVERSED(NL))

      IF (UNIFORM .OR. DISTRIB == 'FOILGRD') THEN

         CALL FOILGRD (NL, ZERO, ONE, WTLIN, WTQUAD, WTSIN, WTCOS,
     >                 SNORML)

      ELSE IF (DISTRIB == 'VINOKUR') THEN

         SNORML(1)  = ZERO
         SNORML(NL) = ONE

         CALL VINOKUR (1, NL, DELTA1, DELTA2, SNORML, LUNCRT, IER)
         IF (IER /= 0) GO TO 910

      ELSE
         GO TO 920
      END IF

      READ (LUNCTL, *) XNOSE, YNOSE, ZNOSE ! XNOSE = CUTS(1) gives singular pt.
      READ (LUNCTL, *) XTAIL, YTAIL, ZTAIL ! XTAIL = CUTS(NCUTS) "   "   "   "

      READ (LUNCTL, *) NREFINED ! # cuts desired in a suggested refined set
      IF (NREFINED < 0) NK = -NREFINED  ! For the volume grid case

      READ (LUNCTL, *) NCUTS

      ALLOCATE (CUTS(NCUTS))

      DO N = 1, NCUTS
         READ (LUNCTL, *) I, CUTS(N)
      END DO

      CLOSE (LUNCTL)

C     ********************* End of control inputs **********************


      P = ZERO ! P will be a point on each cut
      S = ZERO ! S will be a vector defining the cutting plane

      NINTERP = NU + NL - 1 ! # regularized pts. per cut, except for full bodies

      SELECT CASE (XYZCASE)

      CASE ('X')

         IXYZ = 1
         IF (.NOT. HALFBODY) NINTERP = 2 * NINTERP - 1

      CASE ('Y')

         IXYZ = 2

      CASE ('Z')

         IXYZ = 3

      CASE DEFAULT

         WRITE (LUNCRT, '(/, 2A)')
     >      ' Bad input for direction of cuts: ', XYZCASE
         GO TO 999

      END SELECT


C     Winglet or fin to be laid flat before cutting?

      IF (ROTATION /= ZERO) THEN ! Rotate (Y,Z) in-place

         COSTHETA = COSD (ROTATION)
         SINTHETA = SIND (ROTATION)

         DO N = 1, NZONE_ACTIVE
            DO I = IP(N), IP(N) + NPOINTS(N) - 1
               YTEMP = XYZ(2,I)
               ZTEMP = XYZ(3,I)
               XYZ(2,I) =  COSTHETA * YTEMP + SINTHETA * ZTEMP
               XYZ(3,I) = -SINTHETA * YTEMP + COSTHETA * ZTEMP
            END DO
         END DO

      END IF

C     Determine the data ranges to help with the input cut stations and
C     to provide a tolerance for the cutting utilities:

      XMAX = XYZ(1,1)
      XMIN = XMAX
      YMAX = XYZ(2,1)
      YMIN = YMAX
      ZMAX = XYZ(3,1)
      ZMIN = ZMAX

      DO N = 1, NZONE_ACTIVE
         DO I = IP(N), IP(N) + NPOINTS(N) - 1
            XMIN = MIN (XYZ(1,I), XMIN)
            XMAX = MAX (XYZ(1,I), XMAX)
            YMIN = MIN (XYZ(2,I), YMIN)
            YMAX = MAX (XYZ(2,I), YMAX)
            ZMIN = MIN (XYZ(3,I), ZMIN)
            ZMAX = MAX (XYZ(3,I), ZMAX)
         END DO
      END DO

      WRITE (LUNCRT, '(/, (A, 2F14.6))')
     >   ' Data range in X direction: ', XMIN, XMAX,
     >   ' ............. Y .........: ', YMIN, YMAX,
     >   ' ............. Z .........: ', ZMIN, ZMAX,
     >   ' '

      SELECT CASE (IXYZ)
      CASE (1)
         EPSXYZ = EPS * (XMAX - XMIN)
      CASE (2)
         EPSXYZ = EPS * (YMAX - YMIN)
      CASE (3)
         EPSXYZ = EPS * (ZMAX - ZMIN)
      END SELECT

C     Redistribute each cut as one piece if fully uniform has been specified:

      IF (CONSTANT) THEN
         ALLOCATE (SNORM(NINTERP))

         CALL FOILGRD (NINTERP, ZERO, ONE, ONE, ZERO, ZERO, ZERO, SNORM)
      ELSE
         SNORMU_REVERSED(1) = ZERO
         J = NU
         DO I = 2, NU - 1
            SNORMU_REVERSED(I) = SNORMU_REVERSED(I-1) +
     >                           SNORMU(J) - SNORMU(J-1)
            J = J - 1
         END DO
         SNORMU_REVERSED(NU) = ONE

         SNORML_REVERSED(1) = ZERO
         J = NL
         DO I = 2, NL - 1
            SNORML_REVERSED(I) = SNORML_REVERSED(I-1) +
     >                           SNORML(J) - SNORML(J-1)
            J = J - 1
         END DO
         SNORML_REVERSED(NL) = ONE
      END IF

      F = ZERO ! Not used, but avoid NANs in the cutting utilities

      ALLOCATE (X(NINTERP,NCUTS), Y(NINTERP,NCUTS), Z(NINTERP,NCUTS),
     >          DERIVS(NINTERP), AREAS(NCUTS))

      IF (NREFINED > 0) ALLOCATE (REFINED_CUTS(NREFINED))
     >


C     ********************* Cuts at specified stations ***********************

      S(IXYZ) = ONE

      DO J = 1, NCUTS

         IF (ABS (CUTS(J) - XNOSE) < EPSXYZ) THEN ! Singular point

            DO I = 1, NINTERP
               X(I,J) = XNOSE
               Y(I,J) = YNOSE
               Z(I,J) = ZNOSE
            END DO

         ELSE IF (ABS (CUTS(J) - XTAIL) < EPSXYZ) THEN

            DO I = 1, NINTERP
               X(I,J) = XTAIL
               Y(I,J) = YTAIL
               Z(I,J) = ZTAIL
            END DO

         ELSE

            P(IXYZ) = CUTS(J)

C           Cut all zones at this station and regularize the cut:

            CALL CUT (J, NINTERP, X(1,J), Y(1,J), Z(1,J))

         END IF

      END DO


C     Save the structured surface definition in PLOT3D /mgrid form:

      REFLECT = REFLECT .AND. HALFBODY .AND. IXYZ == 1

      N = NINTERP
      IF (REFLECT) N = 2 * N - 1

      IF (FORMATTED) THEN
         WRITE (LUNCUT, '(I2)') 1
         WRITE (LUNCUT, '(3I5)') N, NCUTS, 1
         IF (REFLECT) THEN ! No need to store the reflection in memory
            WRITE (LUNCUT, 1001)
     >         (X(1:NINTERP,J),  X(NINTERP-1:1:-1,J), J = 1, NCUTS),
     >         (Y(1:NINTERP,J), -Y(NINTERP-1:1:-1,J), J = 1, NCUTS),
     >         (Z(1:NINTERP,J),  Z(NINTERP-1:1:-1,J), J = 1, NCUTS)
         ELSE
            WRITE (LUNCUT, 1001) X, Y, Z
         END IF
      ELSE
         WRITE (LUNCUT) 1
         WRITE (LUNCUT) N, NCUTS, 1
         IF (REFLECT) THEN
            WRITE (LUNCUT)
     >         (X(1:NINTERP,J),  X(NINTERP+1:1:-1,J), J = 1, NCUTS),
     >         (Y(1:NINTERP,J), -Y(NINTERP+1:1:-1,J), J = 1, NCUTS),
     >         (Z(1:NINTERP,J),  Z(NINTERP+1:1:-1,J), J = 1, NCUTS)
         ELSE
            WRITE (LUNCUT) X, Y, Z
         END IF
      END IF

      CLOSE (LUNCUT)
      CLOSE (LUNSEG)

      WRITE (LUNCRT, '(/, (2A))')
     >   ' Plottable cuts: ', FILENAME,
     >   ' Raw cut points: ', 'line.segments'

      IF (VOLUME_GRID) THEN

         ALLOCATE (XVOL(NINTERP,NCUTS,NK), YVOL(NINTERP,NCUTS,NK),
     >             ZVOL(NINTERP,NCUTS,NK))

         CALL INTERIOR_VOLUME_GRID ()

         OPEN (LUNCUT, FILE=HAVOC_FILENAME, STATUS='UNKNOWN',
     >         FORM='UNFORMATTED')
         WRITE (LUNCUT) 1
         WRITE (LUNCUT) NINTERP, NCUTS, NK
         WRITE (LUNCUT) XVOL, YVOL, ZVOL
         CLOSE (LUNCUT)
      END IF

C     Estimate the volume from the sectional areas:

      CALL LCSQUAD (NCUTS, CUTS, AREAS, CUTS(1), CUTS(NCUTS),
     >              'M', VOLUME)

      WRITE (LUNCRT, '(/, A, ES16.8)') ' Volume via quadrature:', VOLUME


C     Save HAVOC-type outputs?

      IF (HAVOC) THEN
         OPEN (LUNHVC, FILE=HAVOC_FILENAME, STATUS='NEW', FORM=FORM)

         CALL HAVOC_OUTPUTS () ! Internal procedure

         CLOSE (LUNHVC)
         WRITE (LUNCRT, '(/, 2A)') ' HAVOC inputs:   ', HAVOC_FILENAME
      END IF


C     Generate a refined set of cuts based on curvature?

      IF (NREFINED > 0) THEN

         CALL REDISTRIBUTE_CUTS (LUNCRT, NCUTS, CUTS, NINTERP, X, Y, Z,
     >                           NREFINED, REFINED_CUTS)
      END IF

      GO TO 999


C     Common error handling:

  900 WRITE (LUNCRT, '(/, (A, I6))')
     >   ' Last active zone:   ', LAST_ZONE,
     >   ' No. of zones found: ', NZONES
      GO TO 999

  910 WRITE (LUNCRT, '(/, A, I6, /, A, 2F10.6)')
     >   ' Vinokur distribution failed.     IER: ', IER,
     >   ' Initial & final increments specified: ', DELTA1, DELTA2
      GO TO 999

  920 WRITE (LUNCRT, '(/, 2A)')
     >   ' Bad distribution type (must be upper case): ', DISTRIB

  999 CONTINUE

C     Formats:

 1001 FORMAT (5ES16.8)


C     Internal procedure for SURFACE_CUTS:

      CONTAINS

!        ---------------------------------------------------------------

         SUBROUTINE CUT (ICUT, NINTERP, X, Y, Z)

!        CUT makes one cut of all surface zones and regularizes the cut.
!        ---------------------------------------------------------------

!        Arguments:

         INTEGER, INTENT (IN) ::
     >      ICUT,              ! Cut number
     >      NINTERP            ! Making it an argument eases zeroing X, Y, Z

         REAL, INTENT (OUT), DIMENSION (NINTERP) ::
     >      X, Y, Z            ! Regularized cut

!        Local constants:

         LOGICAL, PARAMETER ::
     >      TRUE = .TRUE.

         CHARACTER, PARAMETER ::
     >      CURVE * 14 = 'open.  closed.',
     >      LOOSE * 1  = 'B',
     >      TIGHT * 1  = 'M',
     >      XZ    * 2  = 'XZ'

!        Local variables:

         INTEGER
     >      I1(1), I14, I34, ICV, IER, IXZ, L, L1, L2, L3, L4,
     >      LINE, N, NLINES, NNODE, NPTS, NSEG, NSEGN

         REAL
     >      DIFF, DYL, DYMIN, R, STOTAL, XCUT, YCUT, ZCUT

         LOGICAL
     >      CLOSED

!        Execution:

         NSEG  = 0
         NNODE = 0

         DO N = 1, NZONE_ACTIVE

!           Add the 2-pt. line segments from this zone for this cut:

            NSEGN = NSEG

            CALL CUT_TRIANGLES (NPOINTS(N), NTRIANGLES(N),
     >                          IV(1,IT(N)), XYZ(1,IP(N)), F(IP(N)),
     >                          SEG, MAXNODE, MAXSEG, NNODE, NDGCUT,
     >                          NSEG, P, S, EPSXYZ, EPSXYZ, IER)

            IF (IER /= 0) THEN
               WRITE (LUNCRT, 1001) 'CUT_TRIANGLES', IER, N,
     >            XYZCASE, P(IXYZ), NSEG, 'Proceeding.'
            END IF

            IF (NSEG > NSEGN) THEN

!              Convert the 2-pt. line segments to contiguous curves:

               CALL CUTORDER (NAD, NADDIM, XCRV, YCRV, ZCRV, FCRV,
     >                        MAXSEG, SEG, NNODE, MAXNODE, NDGCUT,
     >                        NSEG, MAXSEG, EPSXYZ, EPSXYZ, IER)

               IF (IER /= 0) THEN
                  WRITE (LUNCRT, 1001) 'CUTORDER', IER, N, XYZCASE,
     >               P(IXYZ), NSEG, 'Aborting this cut.'
                  EXIT  ! Skip remaining zones
               END IF

            END IF

         END DO  ! Next zone

!        We expect only one line segment per cut, but save all line segments
!        for diagnostic purposes and possibly alternative redistributions.
!        The 2-|3-segment case is now handled specially by connecting them up.

         NLINES = NAD(1)

         IF (NLINES > 0) THEN

            WRITE (LUNSEG, '(I4, A)') ICUT, ' ! Cut number'

            NAD(1) = 1  ! So the following works for line segment 1

            IF (VOLUME_GRID) NLINES = 1  ! Unexplained 1-point "segments" on a
                                         ! test sphere need to be suppressed
            DO LINE = 1, NLINES

               L1 = NAD(LINE)
               L2 = NAD(LINE+1) - 1

               WRITE (LUNSEG, '(I4, I6, A)')
     >            LINE, L2 - L1 + 1, ' ! Segment # and # points'

               WRITE (LUNSEG, '(3ES17.9)')
     >            (XCRV(L), YCRV(L), ZCRV(L), L = L1, L2)

               IF (NLINES == 2 .OR. NLINES == 3) THEN

C                 Bridge what is probably a body gap or gaps at wing/tail roots.
C                 Complete the loop, but exit below after one pass, not NLINES.

                  L3 = L2 + 1
                  L4 = NAD(3) - 1

                  WRITE (LUNSEG, '(I4, I6, A)')
     >               2, L4 - L3 + 1, ' ! Segment # and # points'

                  WRITE (LUNSEG, '(3ES17.9)')
     >               (XCRV(L), YCRV(L), ZCRV(L), L = L3, L4)

                  IF (NLINES == 3) THEN
                     L3 = L4 + 1
                     L4 = NAD(4) - 1

                     WRITE (LUNSEG, '(I4, I6, A)')
     >                  3, L4 - L3 + 1, ' ! Segment # and # points'

                     WRITE (LUNSEG, '(3ES17.9)')
     >                  (XCRV(L), YCRV(L), ZCRV(L), L = L3, L4)
                  END IF

                  WRITE (LUNCRT, '(A, I2, A, I5, A, ES16.8, /, A)')
     >               ' Bridging gap(s) between', NLINES,
     >               ' segments for cut #',
     >               ICUT, '   at', P(IXYZ), ' See line.segments file.'

                  CALL JOIN_LINES (NLINES, L2) ! Internal procedure


               ELSE IF (NLINES > 3) THEN ! No attempt to join these yet

                  IF (LINE == 1) THEN
                     WRITE (LUNCRT, '(2(A, I5), A, ES16.8, /, A)')
     >                  ' WARNING:', NLINES,
     >                  ' segments for cut #', ICUT, '   at', P(IXYZ),
     >                  ' See line.segments file.'
                     X = ZERO
                     Y = ZERO
                     Z = ZERO
                  END IF

                  CYCLE ! Avoid indenting the main case

               END IF

!              Ensure all cuts start at equivalent points and the
!              points are in the same direction.  Looking down the
!              cut axis, for X cuts we use ANTICLOCKWISE order, for
!              compatibility with the right-half-body case, where
!              crown-to-keel order is desirable because the crown
!              is less likely to be flat than the under-body.
!              For wing or fin cuts, the CLOCKWISE order is more
!              conventional, trailing edge to trailing edge,
!              because it handles sharp and blunt trailing edges.

!              First, a duplicate end point may have to be suppressed:

               CLOSED = ABS (XCRV(1) - XCRV(L2)) < EPSXYZ .AND.
     >                  ABS (YCRV(1) - YCRV(L2)) < EPSXYZ .AND.
     >                  ABS (ZCRV(1) - ZCRV(L2)) < EPSXYZ

               NPTS = L2

               IF (CLOSED) THEN
                  L2   = L2 - 1
                  ICV  = 8
               ELSE
                  ICV = 1
               END IF

               IF (IXYZ == 1) THEN ! X cuts should start at the crown
                  I1  = MAXLOC (ZCRV(1:L2)) ! Index of peak Z
                  IXZ = 2
!
!                 This doesn't suffice for an asteroid where the crown is
!                 ill-defined.  Instead, pick the index of the slice point
!                 that is closest to the NOSE-TAIL line and above it in Z:

                  IF (YNOSE /= YTAIL) THEN
                     XCUT = CUTS(ICUT)
                     R    = (XCUT - XNOSE) / (XTAIL - XNOSE)
                     YCUT = (ONE - R)*YNOSE + R*YTAIL
                     ZCUT = (ONE - R)*ZNOSE + R*ZTAIL

                     DYMIN = 1.E+6
                     DO L = 1, L2
                        IF (ZCRV(L) < ZCUT) CYCLE
                        DYL = ABS (YCRV(L) - YCUT)
                        IF (DYL < DYMIN) THEN
                            DYMIN = DYL
                            I1(1) = L
                        END IF
                     END DO
                  END IF

                  WRITE (LUNCRT, '(A, I4, A, ES16.8, 3A, I5, 3A, I5)')
     >               ' Cut #', ICUT, ' at', P(IXYZ), ' is ',
     >               CURVE(ICV:ICV+6), '  # points:', NPTS,
     >               '  Index taken as crown ', XZ(IXZ:IXZ), ':', I1(1)

               ELSE                ! Y & Z cuts start at the trailing edge
                  I1  = MAXLOC (XCRV(1:L2)) ! Index of max. X
                  IXZ = 1

                  WRITE (LUNCRT, '(A, I4, A, ES16.8, 3A, I5, 3A, I5)')
     >               ' Cut #', ICUT, ' at', P(IXYZ), ' is ',
     >               CURVE(ICV:ICV+6), '  # points:', NPTS,
     >               '  Index of max. ', XZ(IXZ:IXZ), ':', I1(1)
               END IF

               IF (I1(1) /= 1) THEN ! Make the curve start at the peak
                  IF (CLOSED) THEN
                     CALL CIRC_SHIFT (L2, I1(1), 1, XCRV)
                     CALL CIRC_SHIFT (L2, I1(1), 1, YCRV)
                     CALL CIRC_SHIFT (L2, I1(1), 1, ZCRV)
                  ELSE ! Circular shifting can't reverse the order!
                     CALL RVERSE (L2, XCRV, XCRV)
                     CALL RVERSE (L2, YCRV, YCRV)
                     CALL RVERSE (L2, ZCRV, ZCRV)
                  END IF
               END IF

!              Reclose the curve?

               IF (CLOSED) THEN
                  L2 = L2 + 1
                  XCRV(L2) = XCRV(1)
                  YCRV(L2) = YCRV(1)
                  ZCRV(L2) = ZCRV(1)
               END IF

!              Ensure known order prior to regularizing:

               I14 = L2 / 4  ! Roughly 1/4 along the curve
               I34 = 3 * I14 ! Roughly 3/4 along

               IF (IXYZ == 1) THEN                ! Body (anticlockwise)
                  IF (HALFBODY) THEN
                     DIFF = ZCRV(I14) - ZCRV(I34) ! Should already be +ve
                  ELSE                            ! Full body section
                     DIFF = YCRV(I14) - YCRV(I34)
                  END IF

                  CALL AREAXY (L2, YCRV, ZCRV, AREAS(ICUT)) ! Ignores order

               ELSE IF (IXYZ == 2) THEN           ! Wing (clockwise)
                  DIFF = ZCRV(I34) - ZCRV(I14)

                  CALL AREAXY (L2, ZCRV, XCRV, AREAS(ICUT))

               ELSE
                  DIFF = YCRV(I34) - YCRV(I14)    ! Fin  (clockwise)

                  CALL AREAXY (L2, XCRV, YCRV, AREAS(ICUT))

               END IF

               IF (DIFF < ZERO) THEN
                  CALL RVERSE (L2, XCRV, XCRV)
                  CALL RVERSE (L2, YCRV, YCRV)
                  CALL RVERSE (L2, ZCRV, ZCRV)
               END IF

!!!            WRITE (50, '(I4, 3ES14.6)')
!!!  >            (L, XCRV(L), YCRV(L), ZCRV(L), L = 1, L2)

!              Regularize in one piece if uniform:

               IF (CONSTANT) THEN

                  CALL CHORDS3D (L2, XCRV, YCRV, ZCRV, TRUE,
     >                           STOTAL, SCRV)
                  CALL LCSFIT   (L2, SCRV, XCRV, TRUE, LOOSE,
     >                           NINTERP, SNORM, X, DERIVS)
                  CALL LCSFIT   (L2, SCRV, YCRV, TRUE, LOOSE,
     >                           NINTERP, SNORM, Y, DERIVS)
                  CALL LCSFIT   (L2, SCRV, ZCRV, TRUE, LOOSE,
     >                           NINTERP, SNORM, Z, DERIVS)
               ELSE

!                 Regularize in pieces:

                  SELECT CASE (XYZCASE)

                  CASE ('X') ! Body cuts

!                    Crown to right water line point:

                     I1  = MAXLOC (YCRV(1:L2))
                     IWR = I1(1)

                     CALL CHORDS3D (IWR, XCRV, YCRV, ZCRV, TRUE,
     >                              STOTAL, SCRV)
                     CALL LCSFIT   (IWR, SCRV, XCRV, TRUE, LOOSE,
     >                              NU, SNORMU, X, DERIVS)
                     CALL LCSFIT   (IWR, SCRV, YCRV, TRUE, TIGHT,
     >                              NU, SNORMU, Y, DERIVS)
                     CALL LCSFIT   (IWR, SCRV, ZCRV, TRUE, LOOSE,
     >                              NU, SNORMU, Z, DERIVS)

!                    Right water line to keel:

                     I1  = MINLOC (ZCRV(1:L2))
                     IKL = I1(1)
                     N   = IKL - IWR + 1

                     CALL CHORDS3D (N, XCRV(IWR), YCRV(IWR), ZCRV(IWR),
     >                              TRUE, STOTAL, SCRV)
                     CALL LCSFIT   (N, SCRV, XCRV(IWR), TRUE, LOOSE,
     >                              NL, SNORML_REVERSED, X(NU), DERIVS)
                     CALL LCSFIT   (N, SCRV, YCRV(IWR), TRUE, TIGHT,
     >                              NL, SNORML_REVERSED, Y(NU), DERIVS)
                     CALL LCSFIT   (N, SCRV, ZCRV(IWR), TRUE, LOOSE,
     >                              NL, SNORML_REVERSED, Z(NU), DERIVS)

                     IF (.NOT. HALFBODY) THEN

!                       Keel to left water line:

                        I1  = MINLOC (YCRV(1:L2)) ! Left water line
                        IWL = I1(1)
                        N   = IWL - IKL + 1
                        L1  = NU + NL - 1

                        CALL CHORDS3D (N, XCRV(IKL),YCRV(IKL),ZCRV(IKL),
     >                                 TRUE, STOTAL, SCRV)
                        CALL LCSFIT   (N, SCRV, XCRV(IKL), TRUE, LOOSE,
     >                                 NL, SNORML, X(L1), DERIVS)
                        CALL LCSFIT   (N, SCRV, YCRV(IKL), TRUE, TIGHT,
     >                                 NL, SNORML, Y(L1), DERIVS)
                        CALL LCSFIT   (N, SCRV, ZCRV(IKL), TRUE, LOOSE,
     >                                 NL, SNORML, Z(L1), DERIVS)

!                       Left water line to crown:

                        N  = L2 - IWL + 1
                        L1 = NU + 2 * (NL - 1)

                        CALL CHORDS3D (N, XCRV(IWL),YCRV(IWL),ZCRV(IWL),
     >                                 TRUE, STOTAL, SCRV)
                        CALL LCSFIT   (N, SCRV, XCRV(IWL), TRUE, LOOSE,
     >                                 NU, SNORMU_REVERSED,X(L1),DERIVS)
                        CALL LCSFIT   (N, SCRV, YCRV(IWL), TRUE, TIGHT,
     >                                 NU, SNORMU_REVERSED,Y(L1),DERIVS)
                        CALL LCSFIT   (N, SCRV, ZCRV(IWL), TRUE, LOOSE,
     >                                 NU, SNORMU_REVERSED,Z(L1),DERIVS)
                     END IF

                  CASE ('Y', 'Z') ! Wing or tail cuts

!                    Lower surface:

                     I1  = MINLOC (XCRV(1:L2)) ! Leading edge index
                     ILE = I1(1)

                     CALL CHORDS3D (ILE, XCRV, YCRV, ZCRV, TRUE,
     >                              STOTAL, SCRV)
                     CALL LCSFIT   (ILE, SCRV, XCRV, TRUE, TIGHT,
     >                              NL, SNORML_REVERSED, X, DERIVS)
                     CALL LCSFIT   (ILE, SCRV, YCRV, TRUE, LOOSE,
     >                              NL, SNORML_REVERSED, Y, DERIVS)
                     CALL LCSFIT   (ILE, SCRV, ZCRV, TRUE, LOOSE,
     >                              NL, SNORML_REVERSED, Z, DERIVS)

!                    Upper surface:

                     N = L2 - ILE + 1

                     CALL CHORDS3D (N, XCRV(ILE), YCRV(ILE), ZCRV(ILE),
     >                              TRUE, STOTAL, SCRV)
                     CALL LCSFIT   (N, SCRV, XCRV(ILE), TRUE, TIGHT,
     >                              NU, SNORMU, X(NL), DERIVS)
                     CALL LCSFIT   (N, SCRV, YCRV(ILE), TRUE, LOOSE,
     >                              NU, SNORMU, Y(NL), DERIVS)
                     CALL LCSFIT   (N, SCRV, ZCRV(ILE), TRUE, LOOSE,
     >                              NU, SNORMU, Z(NL), DERIVS)
                  END SELECT

               END IF

               IF (HAVOC) THEN ! Ignore the "clockwise" control input
                  IF (NLINES /= 2 .AND. NLINES /= 3) CYCLE
               ELSE
                  IF ((IXYZ == 1 .AND. CLOCKWISE) .OR.
     >                (IXYZ /= 1 .AND. .NOT. CLOCKWISE)) THEN
                     CALL RVERSE (NINTERP, X, X)
                     CALL RVERSE (NINTERP, Y, Y)
                     CALL RVERSE (NINTERP, Z, Z)
                  END IF
               END IF

               IF (NLINES == 2 .OR. NLINES == 3) EXIT ! Gap(s) now bridged

            END DO  ! Next line segment (normally no more)

         END IF

  999    RETURN

!        Formats:

 1001    FORMAT (/, ' *** ', A, ' error:', I3, '   Zone:', I3, 3X, A1,
     >           ' cut:', F10.3, '   NSEG:', I5, /, 1X, A)

         END SUBROUTINE CUT

!        --------------------------------------------------------------------

         REAL FUNCTION DIST_SQUARED (L1, L2)

!        Squared distance between two points in X/Y/ZCRV(*).  See JOIN_LINES.
!        --------------------------------------------------------------------

         INTEGER L1, L2

         DIST_SQUARED =
     >      (XCRV(L1) - XCRV(L2)) ** 2 +
     >      (YCRV(L1) - YCRV(L2)) ** 2 +
     >      (ZCRV(L1) - ZCRV(L2)) ** 2

         END FUNCTION DIST_SQUARED

!        --------------------------------------------------------------------

         SUBROUTINE HAVOC_OUTPUTS ()

!        Output the cuts in normalized/shifted form for HAVOC.
!        --------------------------------------------------------------------

!        Local variables:

         INTEGER
     >      I, IT, J
         REAL
     >      XOFFSET, YOFFSET, ZOFFSET, YMAX, YMIN

!        Execution:

         IF (REF_LENGTH /= ONE) THEN
            X = GSCALE * X
            Y = GSCALE * Y
            Z = GSCALE * Z
         END IF

         IF (IXYZ == 1) THEN
            I = 1  ! With J = 1, presumably a body nose
         ELSE
            I = NL ! Current leading edge index
         END IF

         XOFFSET = -X(I,1)
         X = X + XOFFSET
         YOFFSET = -Y(I,1)
         Y = Y + YOFFSET
         ZOFFSET = -Z(I,1)
         Z = Z + ZOFFSET

         WRITE (LUNCRT, '(/, A, F16.6)')
     >      ' Body length used for normalization:   ', REF_LENGTH

         IF (ROTATION /= ZERO) THEN
            WRITE (LUNCRT, '(/, A, F12.6)')
     >         ' HAVOC rotation about X axis (degrees):', -ROTATION
            YTEMP = -YOFFSET
            ZTEMP = -ZOFFSET
            YOFFSET = COSTHETA * YTEMP - SINTHETA * ZTEMP
            ZOFFSET = SINTHETA * YTEMP + COSTHETA * ZTEMP
         END IF

         WRITE (LUNCRT, '(/, A, 3ES16.8)')
     >      ' X, Y, Z offsets for HAVOC:', -XOFFSET, YOFFSET, ZOFFSET

         SELECT CASE (XYZCASE)

         CASE ('X') ! Body assumed

            IF (CONSTANT) THEN ! Pick a right water line index from section 3
               YMAX = -HUGE (YMAX)
               DO I = 1, NU + NL - 1
                  IF (Y(I,3) > YMAX) THEN
                     YMAX = Y(I,3)
                     IWR = I
                  END IF
               END DO
            ELSE
               IWR = NU ! By construction
            END IF

            WRITE (LUNHVC, '(I4)') NCUTS

            IF (HALFBODY) THEN ! Right half assumed

               N = 2 * NINTERP - 1
               DO J = 1, NCUTS ! HAVOC proceeds clockwise from right water line
                  WRITE (LUNHVC, '(I4)') N
                  WRITE (LUNHVC, 1002)
     >               (X(I,J),  Y(I,J), Z(I,J), I = IWR, 1, -1),
     >               (X(I,J), -Y(I,J), Z(I,J), I = 2, NINTERP),
     >               (X(I,J),  Y(I,J), Z(I,J), I = NINTERP - 1, IWR, -1)
               END DO

            ELSE ! Whole body is present, anticlockwise from crown

               DO J = 1, NCUTS
                  WRITE (LUNHVC, '(I4)') NINTERP
                  WRITE (LUNHVC, 1002)
     >               (X(I,J), Y(I,J), Z(I,J), I = IWR, 1, -1),
     >               (X(I,J), Y(I,J), Z(I,J), I = NINTERP - 1, IWR, -1)
               END DO

            END IF

         CASE ('Y', 'Z') ! Assume a finite thickness trailing edge

            DO J = 1, NCUTS
               WRITE (LUNHVC, '(I4)') NINTERP + 1 ! Duplicate LE pt.
               WRITE (LUNHVC, 1002)
     >            (X(I,J), Y(I,J), Z(I,J), I = NL, NINTERP),
     >            (X(I,J), Y(I,J), Z(I,J), I = 1, NL)
            END DO

         END SELECT

         WRITE (LUNHVC, '(I4)') -999
         WRITE (LUNHVC, '(I6)') -99999

 1002    FORMAT (3ES16.8)

         END SUBROUTINE HAVOC_OUTPUTS

!        --------------------------------------------------------------------

         SUBROUTINE JOIN_LINES (NLINES, L2)

!        Join n = NLINES line segments as a single line segment by finding
!        the combination which minimizes the total length of the gaps between
!        segments.  The general case of n segments has 2^n * (n-1)! possible
!        combinations to consider, a number which grows rapidly:
!
!             n        |    2    3    4    5    6   ....
!        -----------------------------------------------
!        2^n & (n-1)!  |    4   16   96  768  7680  ....
!
!        The input packed line segments X/Y/ZCRV(*) are overwritten with the
!        reordered line segments, and pointers NAD(*) are updated accordingly.
!
!        Only the cases n = 2 and 3 are handled here, by explicit itemization.
!        n = 2 bridges a single wing root gap in a fuselage (half body or
!              whole body).
!        n = 3 treats an overlapping pair of wing and tail root gaps in a
!              half body, but a whole body implies n = 4.
!
!        Is anyone up to solving the general case?
!
!        02/23/02  David Saunders  Hurried effort to handle the n = 3 case;
!                                  n = 2 case reimplemented in the process.
!
!        --------------------------------------------------------------------

!        Arguments:

         INTEGER, INTENT (IN) ::
     >      NLINES             ! # line segments, n (only n = 2 or 3 handled)

         INTEGER, INTENT (OUT) ::
     >      L2                 ! Last point of joined line segment

!        Local variables:

         INTEGER
     >      H1, H2, H3, T1, T2, T3,
     >      I1, ICASE(1), NI, NP(3), NT

         REAL
     >      DSQ(16)

         REAL, ALLOCATABLE, DIMENSION (:) ::
     >      X, Y, Z

!        Execution:
!        ----------


!        The n = 2 case can be handled in-place:

         SELECT CASE (NLINES)

         CASE (2) ! n = 2 segments:  4 possible pairs of end-points to connect
!        --------------------------

!                          H1            T1
!                         /  \          /  \
!                        H2   T2       H2   T2

            H1 = NAD(1)
            H2 = NAD(2)
            T1 = H2 - 1
            T2 = NAD(3) - 1

            DSQ(1) = DIST_SQUARED (H1, H2)
            DSQ(2) = DIST_SQUARED (H1, T2)
            DSQ(3) = DIST_SQUARED (T1, H2)
            DSQ(4) = DIST_SQUARED (T1, T2)

            ICASE = MINLOC (DSQ(1:4))

!           The final ordering here doesn't matter - the joined segment is
!           circular-shifted or reversed at the higher level so that all cuts
!           start at a consistent place and go in a consistent direction.

!           |--------| |--------|
!           H1      T1 H2      T2
!
!           Arrange for the closest pair to be at locations T1 & H2, by
!           reversing in-place as necessary:

            SELECT CASE (ICASE(1))
            CASE (1, 2)
               CALL RVERSE (T1, XCRV, XCRV)
               CALL RVERSE (T1, YCRV, YCRV)
               CALL RVERSE (T1, ZCRV, ZCRV)
            CASE DEFAULT
            END SELECT

            SELECT CASE (ICASE(1))
            CASE (2, 4)
               NI = T2 - H2 + 1
               CALL RVERSE (NI, XCRV(H2), XCRV(H2))
               CALL RVERSE (NI, YCRV(H2), YCRV(H2))
               CALL RVERSE (NI, ZCRV(H2), ZCRV(H2))
            CASE DEFAULT
            END SELECT

            L2 = T2 ! New tail of joined segment

            IF (.NOT. HALFBODY) THEN ! Assume it's a full body cut with 2 gaps
               L2 = L2 + 1           ! Close the curve
               XCRV(L2) = XCRV(1)
               YCRV(L2) = YCRV(1)
               ZCRV(L2) = ZCRV(1)
            END IF


         CASE (3) ! n = 3 segments:  16 possible combinations of points to join
!        --------------------------

            NT = NAD(NLINES+1) - 1 ! Total packed length

            ALLOCATE (X(NT), Y(NT), Z(NT)) ! Can't repack in-place for n = 3

            H1 = NAD(1)
            H2 = NAD(2)
            H3 = NAD(3)
            T1 = H2 - 1
            T2 = H3 - 1
            T3 = NT

            NP(1) = T1
            NP(2) = T2 - H2 + 1
            NP(3) = T3 - H3 + 1

            DSQ(1)  = DIST_SQUARED (H1, H2) + DIST_SQUARED (T2, H3)
            DSQ(2)  = DIST_SQUARED (H1, H2) + DIST_SQUARED (T2, T3)
            DSQ(3)  = DIST_SQUARED (H1, T2) + DIST_SQUARED (H2, H3)
            DSQ(4)  = DIST_SQUARED (H1, T2) + DIST_SQUARED (H2, T3)
            DSQ(5)  = DIST_SQUARED (H1, H3) + DIST_SQUARED (T3, H2)
            DSQ(6)  = DIST_SQUARED (H1, H3) + DIST_SQUARED (T3, T2)
            DSQ(7)  = DIST_SQUARED (H1, T3) + DIST_SQUARED (H3, H2)
            DSQ(8)  = DIST_SQUARED (H1, T3) + DIST_SQUARED (H3, T2)

            DSQ(9)  = DIST_SQUARED (T1, H2) + DIST_SQUARED (T2, H3)
            DSQ(10) = DIST_SQUARED (T1, H2) + DIST_SQUARED (T2, T3)
            DSQ(11) = DIST_SQUARED (T1, T2) + DIST_SQUARED (H2, H3)
            DSQ(12) = DIST_SQUARED (T1, T2) + DIST_SQUARED (H2, T3)
            DSQ(13) = DIST_SQUARED (T1, H3) + DIST_SQUARED (T3, H2)
            DSQ(14) = DIST_SQUARED (T1, H3) + DIST_SQUARED (T3, T2)
            DSQ(15) = DIST_SQUARED (T1, T3) + DIST_SQUARED (H3, H2)
            DSQ(16) = DIST_SQUARED (T1, T3) + DIST_SQUARED (H3, T2)

            ICASE = MINLOC (DSQ)

!           Repack the points in X/Y/Z:
!           The order in X/Y/Z/CRV is    H1----T1 H2----T2 H3----T3

!           Avoid too much repetition for the first line segment:

            SELECT CASE (ICASE(1))
            CASE (1:8)
               CALL RVERSE (T1, XCRV, X)
               CALL RVERSE (T1, YCRV, Y)
               CALL RVERSE (T1, ZCRV, Z)
            CASE (9:16)
               X(1:T1) = XCRV(1:T1)
               Y(1:T1) = YCRV(1:T1)
               Z(1:T1) = ZCRV(1:T1)
            CASE DEFAULT
            END SELECT

!           Second and third line segments:

            I1 = NP(1) + NP(3) + 1
            NI = NP(2) + NP(3)

            SELECT CASE (ICASE(1))
            CASE (1)               ! T1----H1 H2----T2 H3----T3
               CALL COPY   (NI,    XCRV(H2), X(H2))
               CALL COPY   (NI,    YCRV(H2), Y(H2))
               CALL COPY   (NI,    ZCRV(H2), Z(H2))
            CASE (2)               ! T1----H1 H2----T2 T3----H3
               CALL COPY   (NP(2), XCRV(H2), X(H2))
               CALL COPY   (NP(2), YCRV(H2), Y(H2))
               CALL COPY   (NP(2), ZCRV(H2), Z(H2))
               CALL RVERSE (NP(3), XCRV(H3), X(H3))
               CALL RVERSE (NP(3), YCRV(H3), Y(H3))
               CALL RVERSE (NP(3), ZCRV(H3), Z(H3))
            CASE (3)               ! T1----H1 T2----H2 H3----T3
               CALL RVERSE (NP(2), XCRV(H2), X(H2))
               CALL RVERSE (NP(2), YCRV(H2), Y(H2))
               CALL RVERSE (NP(2), ZCRV(H2), Z(H2))
               CALL COPY   (NP(3), XCRV(H3), X(H3))
               CALL COPY   (NP(3), YCRV(H3), Y(H3))
               CALL COPY   (NP(3), ZCRV(H3), Z(H3))
            CASE (4)               ! T1----H1 T2----H2 T3----H3
               CALL RVERSE (NP(2), XCRV(H2), X(H2))
               CALL RVERSE (NP(2), YCRV(H2), Y(H2))
               CALL RVERSE (NP(2), ZCRV(H2), Z(H2))
               CALL RVERSE (NP(3), XCRV(H3), X(H3))
               CALL RVERSE (NP(3), YCRV(H3), Y(H3))
               CALL RVERSE (NP(3), ZCRV(H3), Z(H3))
            CASE (5)               ! T1----H1 H3----T3 H2----T2
               CALL COPY   (NP(3), XCRV(H3), X(H2))
               CALL COPY   (NP(3), YCRV(H3), Y(H2))
               CALL COPY   (NP(3), ZCRV(H3), Z(H2))
               CALL COPY   (NP(2), XCRV(H2), X(I1))
               CALL COPY   (NP(2), YCRV(H2), Y(I1))
               CALL COPY   (NP(2), ZCRV(H2), Z(I1))
            CASE (6)               ! T1----H1 H3----T3 T2----H2
               CALL COPY   (NP(3), XCRV(H3), X(H2))
               CALL COPY   (NP(3), YCRV(H3), Y(H2))
               CALL COPY   (NP(3), ZCRV(H3), Z(H2))
               CALL RVERSE (NP(2), XCRV(H2), X(I1))
               CALL RVERSE (NP(2), YCRV(H2), Y(I1))
               CALL RVERSE (NP(2), ZCRV(H2), Z(I1))
            CASE (7)               ! T1----H1 T3----H3 H2----T2
               CALL RVERSE (NP(3), XCRV(H3), X(H2))
               CALL RVERSE (NP(3), YCRV(H3), Y(H2))
               CALL RVERSE (NP(3), ZCRV(H3), Z(H2))
               CALL COPY   (NP(2), XCRV(H2), X(I1))
               CALL COPY   (NP(2), YCRV(H2), Y(I1))
               CALL COPY   (NP(2), ZCRV(H2), Z(I1))
            CASE (8)               ! T1----H1 T3----H3 T2----H2
               CALL RVERSE (NP(3), XCRV(H3), X(H2))
               CALL RVERSE (NP(3), YCRV(H3), Y(H2))
               CALL RVERSE (NP(3), ZCRV(H3), Z(H2))
               CALL RVERSE (NP(2), XCRV(H2), X(I1))
               CALL RVERSE (NP(2), XCRV(H2), X(I1))
               CALL RVERSE (NP(2), XCRV(H2), X(I1))
            CASE (9)               ! H1----T1 H2----T2 H3----T3
               CALL COPY   (NI,    XCRV(H2), X(H2))
               CALL COPY   (NI,    YCRV(H2), Y(H2))
               CALL COPY   (NI,    ZCRV(H2), Z(H2))
            CASE (10)              ! H1----T1 H2----T2 T3----H3
               CALL COPY   (NP(2), XCRV(H2), X(H2))
               CALL COPY   (NP(2), YCRV(H2), Y(H2))
               CALL COPY   (NP(2), ZCRV(H2), Z(H2))
               CALL RVERSE (NP(3), XCRV(H3), X(H3))
               CALL RVERSE (NP(3), YCRV(H3), Y(H3))
               CALL RVERSE (NP(3), ZCRV(H3), Z(H3))
            CASE (11)              ! H1----T1 T2----H2 H3----T3
               CALL RVERSE (NP(2), XCRV(H2), X(H2))
               CALL RVERSE (NP(2), YCRV(H2), Y(H2))
               CALL RVERSE (NP(2), ZCRV(H2), Z(H2))
               CALL COPY   (NP(3), XCRV(H3), X(H3))
               CALL COPY   (NP(3), YCRV(H3), Y(H3))
               CALL COPY   (NP(3), ZCRV(H3), Z(H3))
            CASE (12)              ! H1----T1 T2----H2 T3----H3
               CALL RVERSE (NP(2), XCRV(H2), X(H2))
               CALL RVERSE (NP(2), YCRV(H2), Y(H2))
               CALL RVERSE (NP(2), ZCRV(H2), Z(H2))
               CALL RVERSE (NP(3), XCRV(H3), X(H3))
               CALL RVERSE (NP(3), YCRV(H3), Y(H3))
               CALL RVERSE (NP(3), ZCRV(H3), Z(H3))
            CASE (13)              ! H1----T1 H3----T3 H2----T2
               CALL COPY   (NP(3), XCRV(H3), X(H2))
               CALL COPY   (NP(3), YCRV(H3), Y(H2))
               CALL COPY   (NP(3), ZCRV(H3), Z(H2))
               CALL COPY   (NP(2), XCRV(H2), X(I1))
               CALL COPY   (NP(2), YCRV(H2), Y(I1))
               CALL COPY   (NP(2), ZCRV(H2), Z(I1))
            CASE (14)              ! H1----T1 H3----T3 T2----H2
               CALL COPY   (NP(3), XCRV(H3), X(H2))
               CALL COPY   (NP(3), YCRV(H3), Y(H2))
               CALL COPY   (NP(3), ZCRV(H3), Z(H2))
               CALL RVERSE (NP(2), XCRV(H2), X(I1))
               CALL RVERSE (NP(2), YCRV(H2), Y(I1))
               CALL RVERSE (NP(2), ZCRV(H2), Z(I1))
            CASE (15)              ! H1----T1 T3----H3 H2----T2
               CALL RVERSE (NP(3), XCRV(H3), X(H2))
               CALL RVERSE (NP(3), YCRV(H3), Y(H2))
               CALL RVERSE (NP(3), ZCRV(H3), Z(H2))
               CALL COPY   (NP(2), XCRV(H2), X(I1))
               CALL COPY   (NP(2), YCRV(H2), Y(I1))
               CALL COPY   (NP(2), ZCRV(H2), Z(I1))
            CASE (16)              ! H1----T1 T3----H3 T2----H2
               CALL RVERSE (NP(3), XCRV(H3), X(H2))
               CALL RVERSE (NP(3), YCRV(H3), Y(H2))
               CALL RVERSE (NP(3), ZCRV(H3), Z(H2))
               CALL RVERSE (NP(2), XCRV(H2), X(I1))
               CALL RVERSE (NP(2), YCRV(H2), Y(I1))
               CALL RVERSE (NP(2), ZCRV(H2), Z(I1))
            CASE DEFAULT
            END SELECT

!           Overwrite the input points with the joined segment:

            CALL COPY (NT, X, XCRV)
            CALL COPY (NT, Y, YCRV)
            CALL COPY (NT, Z, ZCRV)

            DEALLOCATE (X, Y, Z)

            L2 = NT ! New tail of joined segment

         CASE DEFAULT

            WRITE (LUNCRT, '(/, A, I5)')
     >         ' Cannot connect this many line segments: ', NLINES

         END SELECT

         END SUBROUTINE JOIN_LINES

!        --------------------------------------------------------------------
         SUBROUTINE INTERIOR_VOLUME_GRID ()
!
!        Turn the regularized slices of a presumably closed body into an
!        interior volume grid for CM and moments of inertia calculations.
!        The surface slices are dimensions NINTERP x NCUTS x 1, and the
!        volume grid is NINTERP x NCUTS x NK.  For each slice, it is formed
!        by joining the point at X = XCUT on the nose-tail line to each
!        regularized surface cut point, and discretizing the resulting radius
!        with NK uniformly-space points.  The intent is to perform CM and
!        moments of inertia calculations on the structured volume grid.
!        --------------------------------------------------------------------

!        Local variables:

         INTEGER :: I, J, K, NBAD
         REAL :: CVOL, DY, DZ, R, TOTAL_VOLUME, VMIN, XCUT, YCUT, ZCUT
         REAL, ALLOCATABLE :: CELL_VOLUMES(:,:,:)

!        Execution:

         DO J = 1, NCUTS

            XCUT = CUTS(J)

            IF (ABS (XCUT - XNOSE) < EPSXYZ) THEN ! Singular point

               XVOL(:,1,:) = XNOSE
               YVOL(:,1,:) = YNOSE
               ZVOL(:,1,:) = ZNOSE

            ELSE IF (ABS (XCUT - XTAIL) < EPSXYZ) THEN

               XVOL(:,NCUTS,:) = XTAIL
               YVOL(:,NCUTS,:) = YTAIL
               ZVOL(:,NCUTS,:) = ZTAIL

            ELSE

               XVOL(:,J,:) = XCUT
               R    = (XCUT - XNOSE) / (XTAIL - XNOSE)
               YCUT = (ONE - R)*YNOSE + R*YTAIL
               ZCUT = (ONE - R)*ZNOSE + R*ZTAIL

               DO I = 1, NINTERP
                  DY = (Y(I,J) - YCUT) / REAL (NK - 1)
                  DZ = (Z(I,J) - ZCUT) / REAL (NK - 1)
                  DO K = 1, NK
                     YVOL(I,J,K) = YCUT + DY * REAL (K - 1)
                     ZVOL(I,J,K) = ZCUT + DZ * REAL (K - 1)
                  END DO
               END DO

            END IF

         END DO

!        Volume check: do the collapsed cells cause trouble?

         ALLOCATE (CELL_VOLUMES(2:NINTERP,2:NCUTS, 2:NK))

         CALL CELLVOL (NINTERP, NCUTS, NK, 1, NINTERP, 1, NCUTS, 1, NK,
     >                 XVOL, YVOL, ZVOL, CELL_VOLUMES, ONE)

         NBAD = 0
         TOTAL_VOLUME = ZERO
         VMIN = 1.E+10

         DO K = 2, NK
            DO J = 2, NCUTS
               DO I = 2, NINTERP
                  CVOL = CELL_VOLUMES(I,J,K)
                  IF (CVOL <= ZERO) NBAD = NBAD + 1
                  IF (CVOL <  VMIN) VMIN = CVOL
                  TOTAL_VOLUME = TOTAL_VOLUME + CVOL
               END DO
            END DO
         END DO

         WRITE (LUNCRT, '(A, I8, (/, A, ES16.8))')
     >      ' # interior volume cells with non-positive volumes:', NBAD,
     >      ' Minimum  cell  volume:', VMIN,
     >      ' Total interior volume:', TOTAL_VOLUME

         DEALLOCATE (CELL_VOLUMES)

         END SUBROUTINE INTERIOR_VOLUME_GRID

      END PROGRAM SURFACE_CUTS

********************************************************************************

      SUBROUTINE REDISTRIBUTE_CUTS (LUNCRT, NCUTS, CUTS, NINTERP,
     >                              X, Y, Z, NREFINED, REFINED_CUTS)

*     Generate a refined set of geometry stations from the curvatures of the
*     profiles associated with a preliminary set.  The desired distribution
*     is adapted smoothly to the variations in the given distribution, which
*     may be relatively coarse.
*
*     02/19/02  D.A.Saunders  Initial implementation.
*     02/20/02       "        Area derivatives are no good at the nose -
*                             switched to (radius of) curvature.
*
********************************************************************************

      IMPLICIT NONE

*     Arguments:

      INTEGER, INTENT (IN) ::
     >   LUNCRT,            ! For desired cuts and any error message
     >   NCUTS              ! Number of preliminary cuts

      REAL, INTENT (IN) ::
     >   CUTS(NCUTS)        ! Preliminary cuts in "X" direction

      INTEGER, INTENT (IN) ::
     >   NINTERP            ! Number of pts. along regularized cuts

      REAL, DIMENSION (NINTERP,NCUTS), INTENT (IN) ::
     >   X, Y, Z            ! Regularized sections

      INTEGER, INTENT (IN) ::
     >   NREFINED           ! Number of refined cuts to compute

      REAL, INTENT (OUT) ::
     >   REFINED_CUTS(NREFINED) ! Redistributed cuts in "X" direction

*     Local constants:

      REAL, PARAMETER ::
     >   ADD = 0.1, EXPONENT = -2.0, ! Ad hoc, body only
     >   ONE = 1.0, ZERO = 0.0

*     Local variables:

      INTEGER
     >   I, IC, IER

      REAL
     >   D1, D2, FRACTION, SCALE, SHIFT,
     >   CONTROLS(NCUTS), CURVATURES(NCUTS), RINDICES(NCUTS),
     >   F(NCUTS), FP(NCUTS), FPP(NCUTS), WORK(5*NCUTS)

*     Execution:

      IF (NREFINED > 0) THEN
         WRITE (LUNCRT, '(/, A)')
     >      ' Curvature-based cut stations are not reliable yet.'
         RETURN
      END IF

*     The control points for ARBDIS should be proportional to the spacing
*     of the revised distribution, and must be positive.
*     Calculate the curvature magnitudes along crown/keel/water lines:

      F(1:NCUTS) = Z(1,1:NCUTS)

      CALL FD12K (NCUTS, CUTS, F, FP, FPP, CURVATURES)

      CONTROLS = ABS (CURVATURES)

      WRITE (LUNCRT, '(/, A, /, (I4, F12.6))')
     >   ' Curvatures:', (I, CURVATURES(I), I = 1, NCUTS)

      IC = (NINTERP + 1) / 2

      F(1:NCUTS) = Y(IC,1:NCUTS)

      CALL FD12K (NCUTS, CUTS, F, FP, FPP, CURVATURES)

      CONTROLS = CONTROLS + ABS (CURVATURES)

      WRITE (LUNCRT, '(/, A, /, (I4, F12.6))')
     >   ' Curvatures:', (I, CURVATURES(I), I = 1, NCUTS)

      F(1:NCUTS) = Z(NINTERP,1:NCUTS)

      CALL FD12K (NCUTS, CUTS, F, FP, FPP, CURVATURES)

      WRITE (LUNCRT, '(/, A, /, (I4, F12.6))')
     >   ' Curvatures:', (I, CURVATURES(I), I = 1, NCUTS)

*     Working with CUTS(*) as the abscissas gives too much bunching
*     at the nose.  Therefore, work with indices, but they need to
*     be transformed to the desired data range.
*     However, ARBDIS works better on [0,1]:

      CALL GETXFORM (CUTS(1), CUTS(NCUTS), ZERO, ONE,
     >               SCALE, SHIFT)

      D1 = ONE / REAL (NCUTS - 1)
      DO I = 1, NCUTS
         RINDICES(I) = REAL (I-1) * D1
         CONTROLS(I) =
     >      (CONTROLS(I) + ABS (CURVATURES(I)) + ADD) ** EXPONENT
      END DO

      CONTROLS(3) = CONTROLS(4) * 0.7
      CONTROLS(2) = CONTROLS(3) * 0.6
      CONTROLS(1) = CONTROLS(2) * 0.5

      WRITE (LUNCRT, '(/, A, /, (I4, 2F12.6))')
     >   ' Curvature-based controls:',
     >   (I, RINDICES(I), CONTROLS(I), I = 1, NCUTS)

      CALL ARBDIS (NREFINED, ZERO, ONE, NCUTS, RINDICES,
     >             CONTROLS, 'A', LUNCRT, WORK, REFINED_CUTS, IER)

      IF (IER /= 0) THEN
         WRITE (LUNCRT, '(/, A, I2)')
     >      ' Bad return from ARBDIS.  IER:', IER, ' Try Vinokur guess.'

         REFINED_CUTS(1) = ZERO
         REFINED_CUTS(NREFINED) = ONE
         FRACTION = (REAL (NCUTS) / REAL (NREFINED)) /
     >              (CUTS(NCUTS) - CUTS(1))
         D1 = FRACTION * (CUTS(2) - CUTS(1)) * 2.
         D2 = FRACTION * (CUTS(NCUTS) - CUTS(NCUTS-1))

         CALL VINOKUR (1, NREFINED, D1, D2, REFINED_CUTS, LUNCRT, IER)

         CALL ARBDIS (NREFINED, ZERO, ONE, NCUTS, RINDICES,
     >                CONTROLS, 'I', LUNCRT, WORK, REFINED_CUTS, IER)

         IF (IER /= 0) WRITE (LUNCRT, '(/, A, I2)')
     >      ' Bad return from ARBDIS.  IER:', IER, ' Proceeding.'

      END IF

      CALL GETXFORM (ZERO, ONE, CUTS(1), CUTS(NCUTS),
     >               SCALE, SHIFT)

      DO I = 1, NREFINED
         REFINED_CUTS(I) = SCALE * REFINED_CUTS(I) + SHIFT
      END DO

      REFINED_CUTS(1) = CUTS(1)
      REFINED_CUTS(NREFINED) = CUTS(NCUTS)

      WRITE (LUNCRT, '(/, A, /, (I4, F12.6))')
     >   ' Refined cut stations:',
     >   (I, REFINED_CUTS(I), I = 1, NREFINED)

      END SUBROUTINE REDISTRIBUTE_CUTS
