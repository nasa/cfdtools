********************************************************************************
*
      PROGRAM HYPER_AERO
*
*     Description:
*     ------------
*
*     Calculate Newtonian-type aerodynamic lift/drag/moment coefficients for
*     a hypersonic vehicle at a range of Mach numbers and angles of attack.
*     A choice of impact pressure methods and shadow methods is provided.
*
*     The geometry OML (outer mold line) is input either as a triangulated
*     surface in FAST or STL formats, or as a surface paneling of the form
*     produced by AEROSURF (structured multiblock surface grid in PLOT3D/FAST
*     format).  If the OML is in surface patch form, it is triangulated here
*     prior to the aerodynamic calculations.
*
*     Coefficients for one or more flight conditions are calculated and
*     tabulated as standard output.  A surface Cp distribution is an optional
*     output, corresponding to the last case calculated.  The associated
*     triangulation may be suppressed (since it may have been produced by
*     an earlier run), and does NOT match the initial triangulation.  This
*     is because the Cp values are associated with the triangle centroids,
*     not their vertices.  Rather than interpolating to the vertices, each
*     Cp is simply written three times and the corresponding triangulation
*     contains all the implied repetitions of points with trivial connectivity.
*
*     Special case for Cp studies:
*     ---------------------------
*
*        If a single triangle is being analyzed (two faces), and the Cp file
*     name is not 'none', that file is used for saving (M, Alpha, Cp, Cp-shadow)
*     at all of the indicated angles of attack (= impact angles), where the
*     range 0 : 180 degrees is suggested.
*
*        Here is a sample such single triangle with two faces:
*
*           3  2  0                               ! points, faces, tetrahedra
*           -.5773502693 -.5773502693 1.154700539 ! X coordinates
*           0. 0. 0.                              ! Y coordinates
*           -1. 1. 0.                             ! Z coordinates
*           1 3 2 1 2 3                           ! Connectivity
*           1 1                                   ! Component numbers
*
*     Sample control file (standard input):
*     -------------------------------------
*
*        INPUT SURFACE      QUADS/TRIS/STL?    NORMAL_FLAGS?    FORMATTED?
*        sharp-v5.xyz           quads          .TRUE.           .TRUE.
*
*        OUTPUT TABULATION      OUTPUT TRIANGULATION         Cp file (last case)
*        sharp-v5.coefs         none                         mach20.cps
*
*        REF. LENGTH   REF. AREA   X_CG   Y_CG   Z_CG    GAMMA
*        39.7534       173.95     20.862    0.     0.    1.2
*
*        MACH1  MACH2  INCREMENT
*        20.    20.    1.
*
*        ALPHA1 ALPHA2 INCREMENT
*        0.     10.    2.
*
*        BETA1  BETA2  INCREMENT
*        0.     0.     1.
*
*        PITCH/YAW/ROLL MOMENT SWITCHES (to overcome convention uncertainties)
*        -1.   1.   1.
*
*        X/Y/Z CONVENTION SWITCHES (for conversion to Jameson convention)
*        1  0  0
*        0  1  0
*        0  0  1
*
*        IMPACT METHOD  SHADOW METHOD
*        1              1
*
*     Notes:
*     ------
*
*     1. The internal X/Y/Z convention is that of the Jameson flow solvers:
*        X points downstream; Y points up; Z points outboard (left wing).
*        The permutation inputs allow for switching the sign of X and for
*        swapping and/or switching the signs of Y and Z (only).
*        Output coordinates are transformed back to the input convention
*        if necessary.
*
*     2. NORMAL_FLAGS applies to structured surface patches only.
*        If NORMAL_FLAGS is TRUE, the surface paneling file should end with the
*        following data (as from the AEROSURF program):
*
*        Components:
*          1  1  1  1  2  2  2  ...  (as many integers as surface patches)
*        Signs to apply to normals:
*         -1 -1 -1 -1  1  1  1  ...  (means for switching the sign of normals)
*
*        Otherwise, all component numbers are set to 1, and all signs to +1,
*        meaning the input paneling should have outward normals implied by the
*        right-hand (i,j,k) rule: k should be out of the surface.
*
*     3. The (final case) Cp distribution will not be written if the file
*        name is 'none' or 'NONE'.  The file is formatted (ASCII text).
*        It contains 3 * NFACES values - the same Cp at each vertex of each
*        triangle, for all triangles.  See above for the special case of a
*        single triangle, where instead Mach, Alpha, Cp, Cp-shadow are saved.
*
*     4. Likewise, the corresponding triangulation will not be saved on disk if
*        the specified file name is 'none' or .NONE'.  Otherwise, the file is
*        formatted (not binary), for FAST.  It contains all vertices of every
*        triangle (meaning replicated points) in order to match the Cp file.
*        A Tecplottable form of this file is also written now.
*
*     5. STL triangulation format:
*
*        solid ABC
*           facet normal 1.  0.  0.
*             outer loop
*               vertex 50. 30. 10.
*               vertex 50. 30. -10.
*               vertex 50. 31. -5.
*             endloop
*           endfacet
*           facet normal ...
*             ::::::::::::
*           endfact
*           ::::::::::::::
*           ::::::::::::::
*        endsolid ABC
*
*     6. Asteroid studies with an ellipsoid promped the option to work with
*        projected area and length for REF_AREA and REF_LENGTH.  To use
*        this option, enter negative values for these quantities.
*        Projected length is taken to be the projected Y data range
*        as a function of Alpha.  Beta = 0 is assumed.
*
*     7. Available impact pressure methods:
*
*        1  Modified Newtonian
*        2  Improved Tangent-Cone
*        3  Tangent-Wedge (Empirical)
*        4  Tangent-Wedge (Oblique shock relationships of NACA TR 1135)
*        5  Newtonian + Prandtl-Meyer expansion method (slope matched)
*        6  Inclined Circular Cone
*        7  Dahlem-Buck method
*        8  Free Molecular Flow
*
*     8. Available shadow methods:
*
*        1. Cp = 0.
*        2. Prandtl-Meyer expansion from free stream
*        3. Newtonian + Prandtl-Meyer
*        4. Base pressure relationship
*
*     History:
*     --------
*
*     10/18/00  D.A.Saunders  Program NEWTONIAN is an initial completion of
*                             James Reuthers beginnings for AEROSURF-type
*                             paneled surfaces, with potentially reusable
*                             utilities in argument-driven form.
*     10/20/00      "         The surface normals are problematic.  Ideally,
*                             the input paneling should be such that they are
*                             outward normals.  If not, a switch can now be
*                             imposed by this program (see NORMAL_FLAGS).
*     10/24/00      "         Added L/D output and second tabulation for
*                             more convenient plotting of polars.
*     11/29/00      "         Store the surface Cp distribution, and write it
*                             and/or the triangulation to disk if specified.
*                             Thanks to Scott Thomas (Raytheon/Ames) for the
*                             suggestion that avoids interpolation.
*     12/13/00      "         FACE_AREA(*), not PROJECTED_AREA(1:3,*), is now
*                             stored at the cost of one multiply per triangle.
*     01/03/01      "         Added option to read a triangulated surface.
*     02/02/01      "         Provided X/Y/Z convention switches.
*                             Arranged for choice of impact pressure method
*                             (all components and all flight conditions for
*                             now, but this needs generalizing), prompted by
*                             installation of a free molecular flow method.
*     02/08/01      "         Installed improved tangent-cone method from HABP.
*     02/12/01      "         Installed tangent-wedge empirical method (HABP).
*     02/14/01      "         Installed tangent-wedge method COMPR from HABP.
*     03/01/01      "         Installed Newtonian + Prandtl-Meyer method of
*                             NEWTPM from HABP.  (COMPR never really needs it.)
*     03/02/01      "         Installed INCLINED_CONE method (ACONE in HABP).
*     03/05/01      "         Installed Dahlem-Buck method from in HABP.
*                             The ACONE method is not right yet.
*     03/16/01      "         Introduced shadow method options in TANGENT_WEDGE_
*                             EMPIRICAL (only).
*     03/23/01      "         Put the impact and shadow methods inside the
*                             loop over surface elements, as opposed to
*                             replicating the loop for each impact method.
*     04/06/01      "         Introduced STL triangulated geometry input option.
*     04/23/01   DAS/JJR      NEWTONIAN program name has changed to HYPER_AERO.
*     07/05/01     DAS        Arranged for saving (M, Alpha, Cp, Cp-shadow) for
*                             the special case of a single triangle.
*     06/07/02      "         Handled zero increments for ranges of Mach and
*                             Alpha; trapped likely misuse of NORMAL_FLAGS.
*     07/19/02      "         STL input had a bug.
*     04/01/08      "         The outputs in FAST format are no longer useful,
*                             but writing an addition file in Tecplot format is
*                             easier than eliminating them.  The file name is
*                             fixed at Tecplotable_Cp.dat.
*     09/15/10      "         The control file is now read from standard input
*                             instead of being hardcoded as hyper_aero.inp.
*                             Any coordinate transformation is applied to the
*                             CG coordinates as well, and is undone before any
*                             surface Cp distributions are written.
*     10/29/14      "         Asteroid studies involving an ellipsoid suggested
*                             the reference area to use for CL, CD, CM should be
*                             the projected area that is not shadowed.  Enter
*                             REF_AREA < 0 to invoke this option.  Likewise, the
*                             moments that depend on REF_LENGTH can be adjusted
*                             to work with projected vertical data range (of Y
*                             in the working coordinate system) by entering
*                             REF_LENGTH < 0.  Only the effect of Alpha is
*                             accounted for, so this may not be right unless
*                             Beta = 0.
*
*     Authors:  David Saunders (ELORET Corporation), James Reuther (TSA Branch).
*               DAS is now with ERC, Inc. at NASA ARC.
*
*     Sponsor:  Reacting Flow Environments Branch, NASA Ames Research Center, CA
*               Now Aerothermodynamics Branch, NASA ARC.
*
********************************************************************************

      IMPLICIT NONE

*     Constants:
*     ----------

      INTEGER, PARAMETER ::
     >   INPUT_QUADS = 1, ! Mnemonics for input geometry type (1 = AEROSURF)
     >   INPUT_TRIS  = 2, ! FAST format
     >   INPUT_STL   = 3, ! STL (CAD) format
     >   LUNXYZ      = 2, ! AEROSURF-type surface paneling, or triangulated OML
     >   LUNOUT      = 3, ! Tabulated results and optional unstructured grid
     >   LUNCPS      = 4, ! For unstructured Cps or 1-triangle Alpha, Cp, Cpsh
     >   LUNINP      = 5, ! Control file on the command line
     >   LUNERR      = 6  ! For file I/O error messages

      REAL, PARAMETER ::
     >   ONE    = 1.,
     >   ZERO   = 0.

*     Variables:
*     ----------

      INTEGER ::
     >   I, I1, I2, IER, IMPACT_METHOD, INPUT_GEOMETRY, IOS, J, K, N,
     >   NALPHA, NBETA, NMACH, NPATCH, NPTS, NFACES, NPOINTS,
     >   NTRIANGLES, SHADOW_METHOD

      INTEGER, ALLOCATABLE, DIMENSION (:) ::
     >   ICOMP, ICOMPONENT, IP

      INTEGER, ALLOCATABLE, DIMENSION (:,:) ::
     >   IV, NIJ

      REAL ::
     >   ALPHA1, ALPHA2, ALPHA_INC, ALPHAIJK,
     >   BETA1,  BETA2,  BETA_INC,  BETAIJK,
     >   RMACH1, RMACH2, RMACH_INC, RMACHIJK,
     >   CLCD,   GAMMA,  REF_AREA,  REF_LEN,
     >   SIGN2,  SIGN3,  SIGNP,  SIGNR,  SIGNY,  SWAP,
     >   XYZ_CG(3), XYZ_SWITCH(3,3)

      REAL(KIND=4) :: CPU1, CPU2

      REAL, ALLOCATABLE, DIMENSION (:) ::
     >   ALPHA, BETA, CP, FACE_AREA, RMACH, SWITCH_NORM, SIGNS

      REAL, ALLOCATABLE, DIMENSION (:,:) ::
     >   XYZ, XYZ_NORM, XYZ_CENTROID

      REAL, ALLOCATABLE, DIMENSION (:,:,:) ::
     >   CL, CD, CS, CM_P, CM_Y, CM_R

      LOGICAL
     >   FORMATTED, NORMAL_FLAGS, SINGLE_TRIANGLE

      CHARACTER
     >   FILEIN * 64, FILEOUT * 64, FILETRI * 64, FILECPS * 64,
     >   GEOMETRY * 1, METHOD (8) * 29, SHADOW (4) * 29, TEXT * 1

*     Storage:

      DATA
     >   METHOD
     >     /'           Modified Newtonian',
     >      '        Improved Tangent-Cone',
     >      '    Tangent-Wedge (Empirical)',
     >      'Tangent-Wedge (Oblique Shock)',
     >      '    Newtonian + Prandtl-Meyer',
     >      '       Inclined Circular Cone',
     >      '     Dahlem-Buck Relationship',
     >      '          Free Molecular Flow'/

      DATA
     >   SHADOW
     >     /'                      Cp = 0.',
     >      'Prandtl-Meyer from free strm.',
     >      '    Newtonian + Prandtl-Meyer',
     >      '   Base pressure relationship'/

*     ----------
*     Execution:
*     ----------

***   OPEN (LUNINP, FILE='hyper_aero.inp', STATUS='OLD', IOSTAT=IOS)

***   IF (IOS /= 0) THEN
***      WRITE (LUNERR, '(/, A)') ' Cannot open hyper_aero.inp.'
***      GO TO 999
***   END IF

***   Read from standard input:

      READ (LUNINP, *)
      READ (LUNINP, *) FILEIN, GEOMETRY, NORMAL_FLAGS, FORMATTED
      READ (LUNINP, *)
      READ (LUNINP, *)
      READ (LUNINP, *) FILEOUT, FILETRI, FILECPS
      READ (LUNINP, *)
      READ (LUNINP, *)
      READ (LUNINP, *) REF_LEN, REF_AREA, XYZ_CG, GAMMA
      READ (LUNINP, *)
      READ (LUNINP, *)
      READ (LUNINP, *) RMACH1, RMACH2, RMACH_INC
      READ (LUNINP, *)
      READ (LUNINP, *)
      READ (LUNINP, *) ALPHA1, ALPHA2, ALPHA_INC
      READ (LUNINP, *)
      READ (LUNINP, *)
      READ (LUNINP, *) BETA1,  BETA2,  BETA_INC
      READ (LUNINP, *)
      READ (LUNINP, *)
      READ (LUNINP, *) SIGNP,  SIGNY,  SIGNR
      READ (LUNINP, *)
      READ (LUNINP, *)
      READ (LUNINP, *) (XYZ_SWITCH(I,1:3), I = 1, 3)
      READ (LUNINP, *)
      READ (LUNINP, *)
      READ (LUNINP, *) IMPACT_METHOD, SHADOW_METHOD

      CLOSE (LUNINP)

*     Open the output tabulation file now in case it already exists:

      OPEN (LUNOUT, FILE=FILEOUT, STATUS='NEW', IOSTAT=IOS)

      IF (IOS /= 0) THEN
         WRITE (LUNERR, '(/, 2A)') ' Tabulation file already exists: ',
     >      FILEOUT
         GO TO 999
      END IF

*     Open the surface geometry file:

      IF (FORMATTED) THEN
         OPEN (LUNXYZ, FILE=FILEIN, STATUS='OLD', IOSTAT=IOS)
      ELSE
         OPEN (LUNXYZ, FILE=FILEIN, STATUS='OLD', FORM='UNFORMATTED',
     >         IOSTAT=IOS)
      END IF

      IF (IOS /= 0) THEN
         WRITE (LUNERR, '(/, 2A)') ' Cannot open geometry file: ',
     >      FILEIN
         GO TO 999
      END IF

*     Decode the geometry type:

      IF      (GEOMETRY == 'q' .OR.
     >         GEOMETRY == 'Q') THEN ! Quads (AEROSURF format)
         INPUT_GEOMETRY = INPUT_QUADS
      ELSE IF (GEOMETRY == 't' .OR.
     >         GEOMETRY == 'T') THEN ! Tris (FAST format)
         INPUT_GEOMETRY = INPUT_TRIS
      ELSE IF (GEOMETRY == 's' .OR.
     >         GEOMETRY == 'S') THEN ! STL (CAD format)
         INPUT_GEOMETRY = INPUT_STL
      END IF

      SINGLE_TRIANGLE = .FALSE.

*     Read the header info (and trap read errors some day?):

      IF (FORMATTED) THEN

         SELECT CASE (INPUT_GEOMETRY)

            CASE (INPUT_QUADS)

               READ (LUNXYZ, *) NPATCH

            CASE (INPUT_TRIS)

               READ (LUNXYZ, *) NPOINTS, NTRIANGLES

            CASE (INPUT_STL) ! Scan the file

               READ (LUNXYZ, *) ! Title

               NTRIANGLES = 0

               DO ! Until EOF, count NPOINTS and NTRIANGLES

                  READ (LUNXYZ, *) ! 'facet' or 'endsolid'
                  READ (LUNXYZ, *, IOSTAT=IOS) ! 'outer' or EOF
                  IF (IOS < 0) EXIT

                  NTRIANGLES = NTRIANGLES + 1
                  READ (LUNXYZ, *)
                  READ (LUNXYZ, *)
                  READ (LUNXYZ, *)
                  READ (LUNXYZ, *)
                  READ (LUNXYZ, *)

               END DO

               NPOINTS = NTRIANGLES * 3

               REWIND (LUNXYZ)

         END SELECT

      ELSE

         SELECT CASE (INPUT_GEOMETRY)

            CASE (INPUT_QUADS)
               READ (LUNXYZ) NPATCH
            CASE (INPUT_TRIS)
               READ (LUNXYZ) NPOINTS, NTRIANGLES
            CASE (INPUT_STL)
               WRITE (*, '(/, A)')
     >            ' STL geometry cannot be unformatted.'
               GO TO 999

         END SELECT

      END IF

*     Allocate storage, and complete the geometry reading:

      SELECT CASE (INPUT_GEOMETRY)

         CASE (INPUT_QUADS) ! AEROSURF-type paneling

            ALLOCATE (ICOMPONENT(NPATCH), IP(NPATCH + 1), NIJ(2,NPATCH),
     >                SWITCH_NORM(NPATCH))

            IF (FORMATTED) THEN
               READ (LUNXYZ, *) (NIJ(1:2,N), I, N = 1, NPATCH) ! 3rd dim. is 1
            ELSE
               READ (LUNXYZ)    (NIJ(1:2,N), I, N = 1, NPATCH)
            END IF

            IP(1)  = 1 ! Pointers to start of each patch in XYZ(1:3,*)
            NFACES = 0

            DO N = 1, NPATCH
               IP(N+1) = IP(N)  +  NIJ(1,N) * NIJ(2,N)
               NFACES  = NFACES + (NIJ(1,N) - 1) * (NIJ(2,N) - 1)
            END DO

            NFACES = NFACES * 2 ! Upper bound; some quad cells may be degenerate
            NPTS = IP(NPATCH+1) - 1

            ALLOCATE (XYZ(3,NPTS))

*           Reorder the coordinates during the read - no need to store 2 forms:

            DO N = 1, NPATCH

               I1 = IP(N)
               I2 = IP(N+1) - 1

               IF (FORMATTED) THEN
                  READ (LUNXYZ, *) ((XYZ(I,J), J = I1, I2), I = 1, 3)
               ELSE
                  READ (LUNXYZ)    ((XYZ(I,J), J = I1, I2), I = 1, 3)
               END IF

            END DO

            IER = 0
            IF (NORMAL_FLAGS) THEN
               IF (FORMATTED) THEN
                  READ (LUNXYZ, *, IOSTAT=IOS)
                  IF (IOS /= 0) THEN
                     IER = 1
                     GO TO 900
                  END IF
                  READ (LUNXYZ, *, IOSTAT=IOS) ICOMPONENT
                  IF (IOS /= 0) THEN
                     IER = 2
                     GO TO 900
                  END IF
                  READ (LUNXYZ, *, IOSTAT=IOS)
                  IF (IOS /= 0) THEN
                     IER = 3
                     GO TO 900
                  END IF
                  READ (LUNXYZ, *, IOSTAT=IOS) SWITCH_NORM
                  IF (IOS /= 0) THEN
                     IER = 4
                     GO TO 900
                  END IF
               ELSE
                  READ (LUNXYZ, IOSTAT=IOS) ICOMPONENT
                  IF (IOS /= 0) THEN
                     IER = 2
                     GO TO 900
                  END IF
                  READ (LUNXYZ, IOSTAT=IOS) SWITCH_NORM
                  IF (IOS /= 0) THEN
                     IER = 4
                     GO TO 900
                  END IF
               END IF
            ELSE
               ICOMPONENT  = 1
               SWITCH_NORM = ONE
            END IF

            CLOSE (LUNXYZ)

*           Triangulate the paneled surface:
*           --------------------------------

            ALLOCATE (IV(3,NFACES), ICOMP(NFACES), ! NFACES = upper bound here
     >                SIGNS(NFACES))

            CALL TRIANGULATE_PATCHES (NPATCH, NIJ, XYZ, IP, ICOMPONENT,
     >                                SWITCH_NORM, REF_LEN, NPTS,
     >                                NFACES, IV, ICOMP, SIGNS)


         CASE (INPUT_TRIS) ! The input geometry is a FAST-type triangulation.

            NPTS   = NPOINTS
            NFACES = NTRIANGLES

*           Special case for Cp studies:

            SINGLE_TRIANGLE = NFACES == 2 .AND.
     >         FILECPS(1:4) /= 'none' .AND. FILECPS(1:4) /= 'NONE'

            IF (SINGLE_TRIANGLE) THEN
               OPEN  (LUNCPS, FILE=FILECPS, STATUS='UNKNOWN')
               WRITE (LUNCPS, '(2A)')
     >            '! Impact pressure method:  ', METHOD(IMPACT_METHOD),
     >            '! Shadow method:           ', SHADOW(SHADOW_METHOD),
     >            '!  Mach    Alpha        Cp Cp-shadow'
            END IF


            ALLOCATE (XYZ(3,NPTS), IV(3,NFACES), ICOMP(NFACES),
     >                SIGNS(NFACES))

            IF (FORMATTED) THEN
               READ (LUNXYZ, *)
     >            ((XYZ(I,N), N = 1, NPTS), I = 1, 3), IV, ICOMP
            ELSE
               READ (LUNXYZ)
     >            ((XYZ(I,N), N = 1, NPTS), I = 1, 3), IV, ICOMP
            END IF

            CLOSE (LUNXYZ)

            SIGNS = ONE

         CASE (INPUT_STL) ! STL (CAD) format

            NPTS   = NPOINTS
            NFACES = NTRIANGLES

            ALLOCATE (XYZ(3,NPTS), IV(3,NFACES), ICOMP(NFACES),
     >                SIGNS(NFACES))

            READ (LUNXYZ, *) ! Title

            NPTS = -2

            DO N = 1, NFACES
               NPTS     = NPTS + 3
               IV(1,N)  = NPTS
               IV(2,N)  = NPTS + 1
               IV(3,N)  = NPTS + 2
               ICOMP(N) = 1

               READ (LUNXYZ, *) ! 'facet'; ignore unit normal information
               READ (LUNXYZ, *) ! 'outer'
               READ (LUNXYZ, *) TEXT, XYZ(1:3,NPTS)
               READ (LUNXYZ, *) TEXT, XYZ(1:3,NPTS+1)
               READ (LUNXYZ, *) TEXT, XYZ(1:3,NPTS+2)
               READ (LUNXYZ, *) ! 'endloop'
               READ (LUNXYZ, *) ! 'endfacet'
            END DO

            CLOSE (LUNXYZ)

            SIGNS = ONE

      END SELECT


*     Allow for foreign X/Y/Z conventions:
*     ------------------------------------

      DO I = 1, 3
         IF (XYZ_SWITCH(I,I) == -ONE) THEN
            DO N = 1, NPTS
               XYZ(I,N) = -XYZ(I,N)
            END DO
            XYZ_CG(I) = -XYZ_CG(I)
         END IF
      END DO

      IF (XYZ_SWITCH(2,2) == ZERO) THEN ! Assume swap of Y/Z & possibly signs
         SIGN2 = XYZ_SWITCH(2,3)
         SIGN3 = XYZ_SWITCH(3,2)
         DO N = 1, NPTS
            SWAP     = XYZ(2,N) * SIGN3
            XYZ(2,N) = XYZ(3,N) * SIGN2
            XYZ(3,N) = SWAP
         END DO
         SWAP      = XYZ_CG(2) * SIGN3
         XYZ_CG(2) = XYZ_CG(3) * SIGN2
         XYZ_CG(3) = SWAP
      END IF

*     Calculate triangulation quantities independent of flight conditions:
*     --------------------------------------------------------------------

      ALLOCATE (FACE_AREA(NFACES),   ! NFACES now excludes collapsed triangles
     >          XYZ_NORM(3,NFACES),
     >          XYZ_CENTROID(3,NFACES))

      CALL TRIANGLE_PROPERTIES (NPTS, NFACES, IV, XYZ, SIGNS,
     >                          FACE_AREA, XYZ_NORM, XYZ_CENTROID)

*     Calculate aerodynamic coefficients and Cps for the indicated cases:
*     -------------------------------------------------------------------

      IF (RMACH1 == RMACH2) RMACH_INC = ONE ! Guard against possible zeros
      IF (ALPHA1 == ALPHA2) ALPHA_INC = ONE
      IF (BETA1  == BETA2 ) BETA_INC  = ONE

      NMACH  = 1 + NINT ((RMACH2 - RMACH1) / RMACH_INC)
      NALPHA = 1 + NINT ((ALPHA2 - ALPHA1) / ALPHA_INC)
      NBETA  = 1 + NINT ((BETA2  - BETA1 ) / BETA_INC )

      ALLOCATE (  CL(NMACH, NALPHA, NBETA),   CD(NMACH, NALPHA, NBETA),
     >            CS(NMACH, NALPHA, NBETA),
     >          CM_P(NMACH, NALPHA, NBETA), CM_Y(NMACH, NALPHA, NBETA),
     >          CM_R(NMACH, NALPHA, NBETA),
     >            CP(NFACES))

      CALL SECOND (CPU1)

      DO K = 1, NBETA

         BETAIJK = BETA1 + BETA_INC * REAL (K - 1)

         DO J = 1, NALPHA

            ALPHAIJK = ALPHA1 + ALPHA_INC * REAL (J - 1)

            DO I = 1, NMACH

               RMACHIJK = RMACH1 + RMACH_INC * REAL (I - 1)

               CALL HIGH_SPEED_AERO
     >             (NFACES, FACE_AREA, XYZ_NORM, XYZ_CENTROID,
     >              XYZ_CG, REF_LEN, REF_AREA, GAMMA,
     >              IMPACT_METHOD, SHADOW_METHOD,
     >              RMACHIJK, ALPHAIJK, BETAIJK,
     >              CL(I,J,K),   CD(I,J,K),   CS(I,J,K),
     >              CM_P(I,J,K), CM_Y(I,J,K), CM_R(I,J,K), CP)

               IF (SINGLE_TRIANGLE)
     >            WRITE (LUNCPS, '(F7.3, F9.4, 2F10.6)')
     >               RMACHIJK, ALPHAIJK, CP ! Cp(1:2)

            END DO
         END DO
      END DO

      CALL SECOND (CPU2)

      CPU2 = CPU2 - CPU1

      WRITE (*, '(/, A, F8.2)')
     >   ' CPU seconds to calculate coefficients:', CPU2

      IF (SINGLE_TRIANGLE) THEN
         CLOSE (LUNCPS)
         WRITE (*, '(2A)')
     >      ' M, Alpha, Cp, Cp-shadow: ', TRIM (FILECPS)
      END IF

      WRITE (*, '(2A)')
     >   ' Impact pressure method:  ', METHOD(IMPACT_METHOD),
     >   ' Shadow method:           ', SHADOW(SHADOW_METHOD)
      WRITE (*, '(A, F8.3)')
     >   ' Gamma:                                ', GAMMA
      N = NMACH * NALPHA * NBETA
      WRITE (*, '(A, I8)')
     >   ' Number of flight conditions:          ', N,
     >   ' Number of surface triangles:          ', NFACES,
     >   ' Number of surface points:             ', NPTS
      WRITE (*, '(2A)')
     >   ' Output coefficients:     ', TRIM (FILEOUT),
     >   ' '

*     Tabulate results:
*     -----------------

      IF (SIGNP /= ONE) CM_P = -CM_P
      IF (SIGNY /= ONE) CM_Y = -CM_Y
      IF (SIGNR /= ONE) CM_R = -CM_R

      WRITE (LUNOUT, '(A, A)')
     >   '!  Mach    Alpha    Beta        CL        CD        CS',
     >   '      CM_P      CM_Y      CM_R       L/D'

      DO K = 1, NBETA ! This table roughly matches Traj inputs
         BETAIJK = BETA1 + BETA_INC * REAL (K - 1)
         DO J = 1, NALPHA
            ALPHAIJK = ALPHA1 + ALPHA_INC * REAL (J - 1)
            DO I = 1, NMACH
               RMACHIJK = RMACH1 + RMACH_INC * REAL (I - 1)
               CLCD = CL(I,J,K) / MAX (CD(I,J,K), 1.E-30)
               CLCD = MIN (MAX (CLCD, -99.), 99.) ! Avoid *** in output

               WRITE (LUNOUT, '(F7.3, F9.4, F8.3, 7F10.5)')
     >            RMACHIJK, ALPHAIJK, BETAIJK,
     >            CL(I,J,K),   CD(I,J,K),   CS(I,J,K),
     >            CM_P(I,J,K), CM_Y(I,J,K), CM_R(I,J,K), CLCD
            END DO
         END DO
      END DO

*     A second tabulation is better suited to standard plots:

      IF (NMACH > 1) THEN

         WRITE (LUNOUT, '(/, A, A)')
     >      '!  Mach    Alpha    Beta        CL        CD        CS',
     >      '      CM_P      CM_Y      CM_R       L/D'

         DO K = 1, NBETA
            BETAIJK = BETA1 + BETA_INC * REAL (K - 1)
            DO I = 1, NMACH
               RMACHIJK = RMACH1 + RMACH_INC * REAL (I - 1)
               DO J = 1, NALPHA
                  ALPHAIJK = ALPHA1 + ALPHA_INC * REAL (J - 1)
                  CLCD = CL(I,J,K) / MAX (CD(I,J,K), 1.E-30)
                  CLCD = MIN (MAX (CLCD, -99.), 99.) ! Avoid *** in plot file

                  WRITE (LUNOUT, '(F7.3, F9.4, F8.3, 7F10.5)')
     >               RMACHIJK, ALPHAIJK, BETAIJK,
     >               CL(I,J,K),   CD(I,J,K),   CS(I,J,K),
     >               CM_P(I,J,K), CM_Y(I,J,K), CM_R(I,J,K), CLCD
               END DO
            END DO
         END DO

      END IF

      CLOSE (LUNOUT)

*     Undo any coordinate transformation?

      IF ((FILETRI(1:4) /= 'none' .AND. FILETRI(1:4) /= 'NONE') .OR.
     >    (FILECPS(1:4) /= 'none' .AND. FILECPS(1:4) /= 'NONE')) THEN

         DO I = 1, 3
            IF (XYZ_SWITCH(I,I) == -ONE) THEN
               DO N = 1, NPTS
                  XYZ(I,N) = -XYZ(I,N)
               END DO
            END IF
         END DO

         IF (XYZ_SWITCH(2,2) == ZERO) THEN ! Assume swap of Y/Z & possibly signs
            SIGN2 = XYZ_SWITCH(2,3)
            SIGN3 = XYZ_SWITCH(3,2)
            DO N = 1, NPTS
               SWAP     = XYZ(2,N) * SIGN2  ! Inverse; really dividing by +/-1
               XYZ(2,N) = XYZ(3,N) * SIGN3
               XYZ(3,N) = SWAP
            END DO
         END IF

      END IF

*     Save a surface triangulation to go with the Cps?

      IF (FILETRI(1:4) /= 'none' .AND. FILETRI(1:4) /= 'NONE') THEN

         OPEN  (LUNOUT, FILE=FILETRI, STATUS='UNKNOWN')

         WRITE (LUNOUT, '(3I10)') 3*NFACES, NFACES, 0
         DO I = 1, 3
            WRITE (LUNOUT, '(8ES15.7)')
     >         (XYZ(I,IV(1,J)), XYZ(I,IV(2,J)), XYZ(I,IV(3,J)),
     >          J = 1, NFACES)
         END DO
         WRITE (LUNOUT, '(15I8)') (I, I = 1, 3*NFACES)
         WRITE (LUNOUT, '(40I3)') ICOMP(1:NFACES)    ! Array may be larger

         CLOSE (LUNOUT)

      END IF

*     Save the (last) Cp distribution?

      IF (FILECPS(1:4) /= 'none' .AND. FILECPS(1:4) /= 'NONE' .AND.
     >    .NOT. SINGLE_TRIANGLE) THEN

         OPEN  (LUNCPS, FILE=FILECPS, STATUS='UNKNOWN')

         WRITE (LUNCPS, '(3I10)') 3*NFACES, 1, 1, 1
         WRITE (LUNCPS, '(8ES14.6)')
     >      (CP(J), CP(J), CP(J), J = 1, NFACES)

         CLOSE (LUNCPS)

         OPEN  (LUNCPS, FILE='Tecplotable_Cp.dat', STATUS='UNKNOWN')

         WRITE (LUNCPS, '(A)') 'TITLE = "Pressure distribution"',
     >      'VARIABLES = "x", "y", "z", "Cp"'
         WRITE (LUNCPS, '(A, I6, A, I6, A)')
     >      'ZONE T="ZONE 1", NODES=', 3*NFACES, ', ELEMENTS=', NFACES,
     >      ', DATAPACKING=POINT, ZONETYPE=FETRIANGLE'
         WRITE (LUNCPS, '(4ES15.7)')
     >      (XYZ(:,IV(1,J)), CP(J),
     >       XYZ(:,IV(2,J)), CP(J),
     >       XYZ(:,IV(3,J)), CP(J), J = 1, NFACES)
         WRITE (LUNCPS, '(3I8)') (I, I = 1, 3*NFACES)

         CLOSE (LUNCPS)

      END IF

      GO TO 999


*     Some of the error handling is reusable:

  900 CONTINUE
      WRITE (LUNERR, '(/, A, /, A, /, (A, I6))')
     >   ' Trouble reading surface normal flags at end of geometry.',
     >   ' Use NORMAL_FLAGS = F if none are present.',
     >   ' IOS:', IOS, ' IER:', IER

  999 CONTINUE

      END PROGRAM HYPER_AERO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE TRIANGULATE_PATCHES (NPATCH, NIJ, XYZ, IP, ICOMPONENT,
     >                                SWITCH_NORM, REF_LEN,
     >                                NPTS, NFACES, IV, ICOMP, SIGNS)
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
!
!     10/12/00  DAS  Initial adaptation of Mark Rimlinger's Struct2Unstruct
!                    routine from UV_MAP.
!
!     Author:  David Saunders, ELORET/NASA Ames Research Center, Mtn. View, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IMPLICIT NONE

!     Arguments:

      INTEGER, INTENT (IN) ::
     >   NPATCH                ! Number of patches

      INTEGER, INTENT (IN) ::
     >   NIJ(2,NPATCH)         ! Patch dimensions

      REAL, INTENT (IN) ::
     >   XYZ(3,*)              ! Surface point coordinates (packed)

      INTEGER, INTENT (IN) ::
     >   IP(NPATCH),           ! Indices in XYZ(1:3,*) of the start of
     >   ICOMPONENT(NPATCH)    ! each patch, and geometry component #s

      REAL, INTENT (IN) ::
     >   SWITCH_NORM(NPATCH),  ! +/-1. handles or right/left-handed patches
     >   REF_LEN               ! Reference length: cell edges less than
                               ! REF_LEN * 1.E-15 are considered to be
                               ! collapsed
      INTEGER, INTENT (OUT) ::
     >   NPTS                  ! Total mumber of points in XYZ(3,*):
                               ! all must be included if the triangulation
                               ! is written to disk in FAST format
      INTEGER, INTENT (OUT) ::
     >   NFACES                ! Number of triangles (tet. faces) found
                               ! (not counting collapsed triangles)
      INTEGER, INTENT (OUT) ::
     >   IV(3,*)               ! Pointers to the three vertices of
                               ! each triangle in XYZ(1:3,*)
      INTEGER, INTENT (OUT) ::
     >   ICOMP(*)              ! Flag for each triangle, set to the
                               ! appropriate ICOMPONENT element
      REAL, INTENT (OUT) ::
     >   SIGNS(*)              ! Flag for each triangle derived from
                               ! the SWITCH_NORM flag for its patch
!     Local constants:

      REAL, PARAMETER ::
     >   ZERO = 0.

!     Local variables:

      INTEGER
     >   I, I0, I1, I2, I3, I4, IC, J, J0, L, N, NI, NJ, NTRI

      REAL
     >   DIAG1, DIAG2, EDGE1, EDGE2, EDGE3, EDGE4, SN, TOL

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
      SUBROUTINE TRIANGLE_PROPERTIES (NPTS, NFACES, IV, XYZ, SIGNS,
     >                                FACE_AREA, XYZ_NORM, XYZ_CENTROID)
!
!     Calculate areas, components of unit normals, and centroids for all
!     faces of a triangulated surface, as needed for calculating Newtonian
!     aerodynamic coefficients.
!
!     08/25/00  J. Reuther   Original implementation.
!     10/18/00  D. Saunders  Argument-driven version.
!     12/13/00     "   "     Return FACE_AREA(*), not PROJECTED_AREA(1:3,*).
!
!     Sponsor:  NASA Ames Research Center, Mountain View, CA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IMPLICIT NONE

!     Arguments:

      INTEGER, INTENT (IN) ::
     >   NPTS,                   ! Total mumber of points in XYZ(3,*)
     >   NFACES,                 ! Number of triangles
     >   IV(3,*)                 ! Pointers to the three vertices of
                                 ! each triangle in XYZ(1:3,*)
      REAL, INTENT (IN) ::
     >   XYZ(3,NPTS),            ! Surface point coordinates
     >   SIGNS(NFACES)           ! +1. for outward normal, else -1.

      REAL, INTENT (OUT) ::
     >   FACE_AREA(NFACES)       ! Triangle areas >= 0.

      REAL, INTENT (OUT), DIMENSION (3,NFACES) ::
     >   XYZ_NORM,               ! Components of triangle outward unit normals
     >   XYZ_CENTROID            ! Triangle centroids

!     Local constants:

      REAL, PARAMETER ::
     >   HALF = 0.5, ONE = 1., SMALL = 1.E-30, THIRD = 1./ 3.

!     Local variables:

      INTEGER :: I, I1, I2, I3
      REAL    :: AREAX, AREAY, AREAZ, RAREA, SWITCH

!     Execution:

      DO I = 1, NFACES

         I1 = IV(1,I)
         I2 = IV(2,I)
         I3 = IV(3,I)
         SWITCH = HALF * SIGNS(I)

         AREAX = ((XYZ(2,I2) - XYZ(2,I1))  *
     >            (XYZ(3,I3) - XYZ(3,I1))  -
     >            (XYZ(3,I2) - XYZ(3,I1))  *
     >            (XYZ(2,I3) - XYZ(2,I1))) * SWITCH
         AREAY = ((XYZ(3,I2) - XYZ(3,I1))  *
     >            (XYZ(1,I3) - XYZ(1,I1))  -
     >            (XYZ(1,I2) - XYZ(1,I1))  *
     >            (XYZ(3,I3) - XYZ(3,I1))) * SWITCH
         AREAZ = ((XYZ(1,I2) - XYZ(1,I1))  *
     >            (XYZ(2,I3) - XYZ(2,I1))  -
     >            (XYZ(2,I2) - XYZ(2,I1))  *
     >            (XYZ(1,I3) - XYZ(1,I1))) * SWITCH

         FACE_AREA(I)  = SQRT (AREAX ** 2 + AREAY ** 2 + AREAZ ** 2)
         RAREA         = ONE / MAX (FACE_AREA(I), SMALL)

         XYZ_NORM(1,I) = AREAX * RAREA ! Outward unit normal components
         XYZ_NORM(2,I) = AREAY * RAREA
         XYZ_NORM(3,I) = AREAZ * RAREA

         XYZ_CENTROID(1,I) = (XYZ(1,I1) + XYZ(1,I2) + XYZ(1,I3)) * THIRD
         XYZ_CENTROID(2,I) = (XYZ(2,I1) + XYZ(2,I2) + XYZ(2,I3)) * THIRD
         XYZ_CENTROID(3,I) = (XYZ(3,I1) + XYZ(3,I2) + XYZ(3,I3)) * THIRD

      END DO

      END SUBROUTINE TRIANGLE_PROPERTIES
********************************************************************************
*
      SUBROUTINE HIGH_SPEED_AERO (NFACES, FACE_AREA, XYZ_NORM,
     >                            XYZ_CENTROID, XYZ_CG,
     >                            REF_LEN, REF_AREA, GAMMA,
     >                            IMPACT_METHOD, SHADOW_METHOD,
     >                            RMACH, ALPHA, BETA,
     >                            CL, CD, CS, CM_P, CM_Y, CM_R, CP)
*
*     This routine drives multiple Newtonian-type impact and shadow methods
*     of calculating aerodynamic forces and moments for the given triangular
*     elements at the given Mach, Alpha, and Beta.
*
*     The (x,y,z) convention is that of the Jameson flow solvers, namely
*     x > 0 downstream, y > 0 up, and z > 0 for left wing.
*
*     03/23/01  David Saunders  Switched from original approach (single method
*               ELORET/ARC      inside the loop over elements, with separate
*               Now ERC, Inc.   routines for each method) to all impact and
*                               shadow methods inside the loop over elements.
*                               One method of the appropriate type is applied
*                               to all elements; eventually, each element may
*                               be tagged with the desired method(s).
*
*     Sponsor:  TSA Branch, NASA Ames Research Center, Mountain View, CA.
*
********************************************************************************

      IMPLICIT NONE

*     Arguments:

      INTEGER, INTENT (IN) ::
     >   NFACES                  ! Number of triangles

      REAL, INTENT (IN) ::
     >   FACE_AREA(NFACES)       ! Triangle areas

      REAL, INTENT (IN), DIMENSION (3,NFACES) ::
     >   XYZ_NORM,               ! Components of triangle outward unit normals
     >   XYZ_CENTROID            ! Triangle centroids

      REAL, INTENT (IN) ::
     >   XYZ_CG(3),              ! Moment center
     >   REF_LEN,                ! Reference length
     >   REF_AREA,               ! Reference area
     >   GAMMA                   ! Ratio of specific heats

      INTEGER, INTENT (IN) ::
     >   IMPACT_METHOD,          ! See the Program Newtonian for the options
     >   SHADOW_METHOD           ! in both cases

      REAL, INTENT (IN) ::
     >   RMACH,                  ! Free-stream Mach number
     >   ALPHA, BETA             ! Angle of attack and yaw angle (degrees)

      REAL, INTENT (OUT) ::
     >   CL, CD, CS,             ! Aerodynamic coefficients
     >   CM_P, CM_Y, CM_R,
     >   CP(NFACES)              ! Pressure coefficients at each triangle

*     Local constants:

      REAL, PARAMETER ::
     >   HALF    = 0.5,  ONE = 1.,  R4TH = 0.25,  TWO = 2.,  ZERO = 0.,
     >   PI      = 3.14159265358979,
     >   DEG2RAD = 0.0174532925199433,
     >   RAD2DEG = 57.2957795130823,
     >   RTPI    = 1.7724538509055160,   ! sqrt(pi)
     >   RTPIR   = 0.5641895835477563    ! 1/sqrt(pi)

*     Local variables:

      INTEGER ::
     >   I, IFIRST, IPRANDTL, MER

      REAL ::
     >   ALFP, ALFW, CA, CB, CC, CI, CPSM1, CPSP1, CPSTAG, CP_BASE,
     >   CTHETA, CX, CY, CZ, DELTA, EMN, EMNS, EMP, ERFSC, EXPSC,
     >   GMM1, GMP1, GMQR, GMSQ, PC, PHI, PHIW, R, RAD, RMSQ,
     >   S, S1, S2, S3, SA, SB, SC, SCSQ, SS, STHETA,
     >   T1, T2, T3, TB, TMP1, TMP2, VC, XM, YM, ZM,
     >   ANGLE(3), AV(3), BS(8), FS(8)
      REAL ::
     >   PROJECTED_AREA, YMAX, YMIN, YROT

*     Execution:

      ANGLE  = ONE  ! Most free-stream properties are not used yet
      AV     = ONE
      FS     = ONE
      FS(6)  = RMACH
      IFIRST = 1    ! Activate the initial iteration in NEWTPM

      RAD = DEG2RAD * ALPHA
      SA  = SIN (RAD)
      CA  = COS (RAD)
      RAD = DEG2RAD * BETA
      SB  = SIN (RAD)
      CB  = COS (RAD)
      TB  = TAN (RAD)
      CC  = CB * CA
      CI  = CB * SA

      CX  = ZERO
      CY  = ZERO
      CZ  = ZERO
      XM  = ZERO
      YM  = ZERO
      ZM  = ZERO

      T1  = HALF * (GAMMA + ONE) * RMACH
      T2  = HALF * T1
      T3  = ONE / (T2 * RMACH)

      S    = RMACH * SQRT (HALF * GAMMA) ! Speed ratio for FMF method
      RMSQ = RMACH * RMACH
      GMSQ = GAMMA * RMSQ
      GMQR = TWO / GMSQ

      IF (GAMMA == ONE) THEN
         CPSTAG = GMQR * (RMSQ - ONE)
      ELSE
         GMP1   = GAMMA + ONE
         GMM1   = GAMMA - ONE
         TMP1   = GMSQ  * TWO - GMM1
         TMP2   = (GMP1 * RMACH) ** 2 / (TMP1 + TMP1)
         CPSTAG = GMQR * ((TMP2 ** (GAMMA/GMM1)) * (TMP1 / GMP1) - ONE)
      END IF

      CPSP1 = ONE + CPSTAG
      CPSM1 = ONE - CPSTAG

      IF (TB == ZERO) THEN
         IF (SA == ZERO) THEN
            PHIW = PI
         ELSE
            PHIW = ZERO
         END IF
      ELSE
         PHIW = ATAN2 (-TB, SA)
      END IF

      ALFW = ACOS (CC)

      RMSQ = ONE / RMSQ
      CP_BASE = RMSQ * (0.57 * RMSQ - ONE)

      PROJECTED_AREA = ZERO
      YMIN =  1.E+10
      YMAX = -YMIN

*     Loop over all triangles and integrate forces:

      DO I = 1, NFACES

*        Calculate the dot product of the inward pointing normal
*        of the triangle and the free-stream direction:

         CTHETA = -(CC * XYZ_NORM(1,I) +    ! Cosine (theta)
     >              CI * XYZ_NORM(2,I) +    ! Impact angle = pi/2 - theta
     >              SB * XYZ_NORM(3,I))

*        If the triangle points into the free-stream, apply an impact method:

         IF (CTHETA > ZERO) THEN

            PROJECTED_AREA = PROJECTED_AREA + FACE_AREA(I) * CTHETA

*           The working coodinate system has Y up:

            YROT = XYZ_CENTROID(1,I)*SA + XYZ_CENTROID(2,I)*CA
            YMAX = MAX (YROT, YMAX)  ! This is probably wrong if BETA /= 0
            YMIN = MIN (YROT, YMIN)

            SELECT CASE (IMPACT_METHOD)


               CASE (1) ! Modified Newtonian

                  CP(I) = CPSTAG * (CTHETA * CTHETA)


               CASE (2) ! Improved tangent-cone method from HABP

                  ANGLE(1) = 90. - RAD2DEG * ACOS (MIN (ONE, CTHETA)) ! Delta

                  CALL CONE (FS, AV, GAMMA, ANGLE, 0, CP(I), BS)


               CASE (3) ! Tangent-wedge empirical method from HABP

                  EMNS  = T1 * CTHETA + EXP (-T2 * CTHETA)
                  CP(I) = T3 * (EMNS * EMNS - ONE)


               CASE (4) ! Tangent-wedge oblique shock method from HABP

                  ANGLE(2) = 90. - RAD2DEG * ACOS (MIN (ONE, CTHETA)) ! Delta

                  CALL COMPR (FS, AV, GAMMA, ANGLE, 0, CPSTAG, IFIRST,
     >                        CP(I), BS)

                  IFIRST = 0 ! Suppress NEWTPM's initial iteration


               CASE (5) ! Newtonian + Prandtl-Meyer (NEWTPM in HABP)

                  ANGLE(2) = 90. - RAD2DEG * ACOS (MIN (ONE, CTHETA)) ! Delta

                  CALL NEWTPM (FS, AV, GAMMA, ANGLE, 0, CPSTAG, IFIRST,
     >                         ONE, CP(I), EMN, BS, IPRANDTL)

                  IFIRST = 0 ! Suppress NEWTPM's initial iteration


               CASE (6) ! Inclined circular cone (ACONE method in HABP)

                  ANGLE(1) = RAD2DEG * ASIN (-XYZ_NORM(1,I)) ! Cone angle

                  PHI = ZERO
                  IF (XYZ_NORM(3,I) /= ZERO) THEN
                     IF (XYZ_NORM(2,I) /= ZERO) THEN ! Meridian angle
                        PHI = ATAN2 (XYZ_NORM(3,I), -XYZ_NORM(2,I))
                     END IF
                  END IF

                  PHI  = PHI - PHIW                        ! Dev. fr. wind plane
                  VC   = CI * COS (PHI) - SB * SIN (PHI)
                  EMP  = RMACH * SQRT (CC**2 + VC**2)      ! Mach in PHI plane

                  ALFP = ZERO
                  IF (VC /= ZERO) THEN
                     IF (CC /= ZERO) ALFP = ATAN2 (VC, CC) ! Angle with X axis
                  END IF

                  CALL ACONE (FS, AV, GAMMA, ANGLE, ALFW, PHI, EMP,
     >                        ALFP, 0, CP(I), BS)


               CASE (7) ! Dahlem-Buck method from HABP

                  DELTA = 90. - RAD2DEG * ACOS (MIN (ONE, CTHETA))

                  IF (DELTA > 22.5) THEN
                     PC = TWO
                  ELSE
                     PC = MIN (5., ONE + ONE /
     >                    (SIN (ABS ((4. * DEG2RAD) * DELTA))) ** 0.75)
                  END IF


               CASE (8) ! Free molecular flow method from HABP

*                 Theta is in [0, pi], so sin (theta) cannot be negative.

                  STHETA = SQRT ((ONE + CTHETA) * (ONE - CTHETA))
                  SS     = S * STHETA             ! Above avoids trig. fns.
                  SC     = S * CTHETA
                  SCSQ   = SC * SC

                  EXPSC  = EXP (-SCSQ)

                  CALL CALERF (SC, 0, ERFSC)    ! Erf (SC)

                  R = (CPSP1 * SC * RTPIR   + CPSM1*R4TH) * EXPSC +
     >                (CPSP1 *(SCSQ + HALF) + CPSM1*R4TH  * RTPI * SC) *
     >                (ONE + ERFSC)

                  CP(I) = GMQR * (R - ONE)


               CASE DEFAULT

                  WRITE (*, *) 'Bad impact method: ', IMPACT_METHOD
                  STOP

            END SELECT


         ELSE ! Element does not see the free stream


            SELECT CASE (SHADOW_METHOD)


               CASE (1) ! Cp = 0.

                  CP(I) = ZERO


               CASE (2) ! Prandtl-Meyer expansion from free stream

                  ANGLE(2) = 90. - RAD2DEG * ACOS (MIN (ONE, CTHETA)) ! Delta
                  ANGLE(2) = ABS (ANGLE(2))

                  CALL EXPAND (FS, AV, GAMMA, ANGLE, 2, MER, CP(I), BS)


               CASE (3) ! Newtonian + Prandtl-Meyer (NEWTPM in HABP)

                  ANGLE(2) = 90. - RAD2DEG * ACOS (MIN (ONE, CTHETA)) ! Delta

                  CALL NEWTPM (FS, AV, GAMMA, ANGLE, 0, CPSTAG, IFIRST,
     >                         ONE, CP(I), EMN, BS, IPRANDTL)

                  IFIRST = 0 ! Suppress NEWTPM's initial iteration


               CASE (4) ! Base pressure relationship for Mach >= 1

                  CP(I) = CP_BASE


               CASE DEFAULT

                  WRITE (*, *) 'Bad shadow method: ', SHADOW_METHOD
                  STOP

            END SELECT

         END IF

*        Accumulate forces and moments:

         IF (ABS (CP(I)) > 1.E-7) THEN

            PC = -CP(I) * FACE_AREA(I)
            S1 = PC * XYZ_NORM(1,I)
            S2 = PC * XYZ_NORM(2,I)
            S3 = PC * XYZ_NORM(3,I)
            CX = CX + S1
            CY = CY + S2
            CZ = CZ + S3
            XM = XM - S2 * (XYZ_CENTROID(3,I) - XYZ_CG(3))
     >              + S3 * (XYZ_CENTROID(2,I) - XYZ_CG(2))
            YM = YM - S3 * (XYZ_CENTROID(1,I) - XYZ_CG(1))
     >              + S1 * (XYZ_CENTROID(3,I) - XYZ_CG(3))
            ZM = ZM - S1 * (XYZ_CENTROID(2,I) - XYZ_CG(2))
     >              + S2 * (XYZ_CENTROID(1,I) - XYZ_CG(1))

         END IF

      END DO


*     Transform to the wind axes:

      IF (REF_AREA > ZERO) THEN
         T1 = ONE / REF_AREA
      ELSE
         T1 = ONE / PROJECTED_AREA
      END IF

      CX = CX * T1
      CY = CY * T1
      CZ = CZ * T1

      CD =       CC * CX +      CI * CY + SB * CZ
      CL =      -SA * CX +      CA * CY
      CS = -CA * SB * CX - SA * SB * CY + CB * CZ

      IF (REF_LEN > ZERO) THEN
         T1 = T1 / REF_LEN
      ELSE
         T1 = T1 / (YMAX - YMIN)
      END IF

      CM_R = XM * T1
      CM_Y = YM * T1
      CM_P = ZM * T1

      END SUBROUTINE HIGH_SPEED_AERO
c-------------------------------------------------------------------------------
c
      subroutine acone (fs, av, gamma, angle, alf, phi, emp, alfp,
     >                  iflow, cp, bs)
c
c     This routine solves for the flow properties about a circular cone
c     at angle of attack using ARC CP # 792 and extended by D.N. Smyth.
c     See also the article by D.J.Jones in the AIAA Journal, Feb. 1972,
c     Vol. 10, No. 2, pp.234-236.
c
c     History:
c
c     c. 1968   HABP Mark III:                 A.E.Gentry, McDonnell Douglas
c     c. 1990   Improved tangent-cone method:  C.I.Cruz/G.J.Sova, LaRC/Rockwell;
c                                              D.N.Smyth (LaRC) co-author?
c     03/02/01  Partial restructuring:         D.A.Saunders, ELORET/NASA Ames
c
c-------------------------------------------------------------------------------

      implicit none

c     Arguments:
c     ----------

      real, intent (inout) ::      ! fs(6) is temporarily modified then restored
     >   fs(8)              ! freestream properties:
                            ! 1 = density, lbm/ft^3    2 = pressure, lbf/ft^2
                            ! 3 = temperature (units?) 4 = speed of sound?
                            ! 5 = viscosity?           6 = Mach number
                            ! 7 = M * speed of sound   8 = rho * M * a / mu

      real, intent (in) ::
     >   av(3)              ! freestream gas constants:
                            ! 1 = ??? temperature (units?)
                            ! 2 = ??? viscosity (units?)
                            ! 3 = ??? (units?)

      real, intent (in) ::
     >   gamma              ! ratio of specific heats

      real, intent (inout) ::
     >   angle(3)           ! angle(1) = impact angle (degrees, input)
                            ! angle(2) = |angle(1)|   (degrees, output)
                            ! angle(3) = shock angle? (degrees, output)

      real, intent (in) ::
     >   alf,               ! Angle of attack in the windward plane (radians)
     >   phi,               ! Deviation of meridian angle from wind plane (rad.)
     >   emp,               ! Mach number in plane of meridian angle
     >   alfp               ! Angle between emp and X axis (radians)

      integer, intent (in) ::
     >   iflow              ! iflow = 0 means return Cp only/no flow properties
                            ! iflow = 1 means return bs(1:8) & angle(3) as well

      real, intent (out) ::
     >   cp                 ! surface element pressure coefficient

      real, intent (out) ::
     >   bs(8)              ! local properties analogous to fs(*)

c     Local constants:
c     ----------------

      real, parameter ::
     >   deg2rad = 0.0174532925,
     >   rad2deg = 57.29577951,
     >   piby2   = 1.570796327

c     Local variables:
c     ----------------

      integer
     >   i, is, isdet, j, mer

      real
     >   ac, bk, ck, cosa, cosp, cos2p, cpcone, cpmin, cpstag, cpx, dcr,
     >   deli, dpsk, dsc, dsr, em, em2, em2i, ema, emc, emc2, emd, emns,
     >   ems, emsurf, fdp, fs6, g, gm1, gp1, p2, pcone, pmin, ps, pshk,
     >   pstag, psurf, scd, sina, sinp, t, tst1, x

      real
     >   a(18), aa(6) ! a(*) contains Fourier coefs. for parameters aa(*)
                      ! in the Jones formula (which assumes gamma = 1.4)

      data a / -0.076570,  1.477500,  0.064669,
     $          0.423390,  0.132410,  0.035871,
     1         -0.002083, -0.075797, -0.019230,
     2          0.298980, -0.100110,  0.295890,
     3         -0.997270, -0.417510,  0.068791,
     4         -0.039442,  0.104220,  0.063801/

c     Execution:
c     ----------

      g      = gamma
      gm1    = g - 1.
      gp1    = g + 1.
      cpcone = 0.0
      cpstag = 0.0
      cpx    = 0.0
      dpsk   = 0.0
      dsc    = 0.0
      em     = fs(6)
      ema    = em
      em2    = em**2
      cosa   = cos (alf)
      sina   = sin (alf)
      cosp   = cos (phi)
      sinp   = sin (phi)
      pmin   = 0.3
      cpmin  = (pmin - 1.) / (0.5 * g *em2)
      isdet  = 2
      dcr    = angle(1) * deg2rad

      if (dcr < 0.) go to 35

c     Test for zero cone angle:

      if (dcr < 1.e-4 ) go to 5
      if (dcr < 0.0524) go to 4

c     Jones' method (AIAA Journal, Feb. 1972, pp. 234-236):

      em2i   = 1./ em2
      cos2p  = cos (2. * phi)
      t      = sin (2.9 * dcr) * tan (dcr)
      ac     = alf / dcr
      ema    = em

      call cone (fs, av, g, angle, 1, cpcone, bs) ! Request flow conditions

      pcone  = 0.5 * g * em2 * cpcone + 1.
      ems    = bs(6)
      dsc    = angle(3) * deg2rad ! angle(3) is the shock angle from CONE

      do i = 1, 6
         j = 3*(i-1) + 1
         aa(i) = a(j) + a(j+1)*cosp + a(j+2)*cos2p
      end do

      cpx = ((aa(1)*t + aa(2)*t*em2i + aa(3)*em2i) +
     >       (aa(4)*t + aa(5)*t*em2i + aa(6)*em2i) * ac) * ac
      go to 40

c     Axial flow component:

    4 ema   = em*cosa
      fs6   = fs(6)
      fs(6) = ema

      call cone (fs, av, g, angle, 1, cpcone, bs)

      fs(6) = fs6
      pcone = 0.5 * g * ema**2 * cpcone + 1.
      ems   = bs(6)
      dsc   = angle(3) * deg2rad

    5 if (alf < 1.e-4) go to 40

c     Normal flow component:

      emc2  = em2 * sina**2
      if (emc2 > 1.) go to 20

c     Subsonic normal component:

      pstag = (1. + .5*gm1*emc2)**(g/gm1)
      go to 30

c     Supersonic normal component:

   20 pstag  = (0.5*gp1*emc2)**(g/gm1)*(gp1/(2.*g*emc2-gm1))**(1./gm1)

   30 cpstag = ((pstag - 1.)/(0.5*g*em2)) * cosp**2
      if (cpstag < 0.) cpstag = 0.

      if (dcr < 1.e-4) go to 40

c     Cross product term - empirical fit from cp #792:

      scd   = sin (dcr) * cos (dcr)
      x     = piby2 / (scd * sqrt (em2 - 1.))
      bk    = 1.95 + 0.07 * cos (x)
      cpx   = 2. * bk * scd * cosp * sina * cosa
      go to 40

c     Cone flow is detached or undefined (dcr < 0.).
c     Obtain impact angle (in meridian plane) and use impact methods:

   35 continue
      deli = alfp + dcr
      if (deli > 0.) go to 36

c     Expansion flow:

      angle(2) = rad2deg * abs (deli)
      fs(6)    = emp
      isdet    = 2 ! This means EXPAND returns just Cp and bs(6)

      call expand (fs, av, g, angle, isdet, mer, cp, bs)

      emsurf   =  bs(6)
      em2      =  emp**2
      psurf    =  0.5*g*em2*cp + 1.
      if (psurf < pmin) psurf = pmin
      go to 50

c     Compression flow - use tangent-cone empirical:

   36 ck    = (2. * gp1 / (g + 3.)) * emp * sin (deli)
      emns  = (ck + exp (-ck))**2
      cp    = (8. * gp1 * emns / ((3. * g + 5.) * emns + 2.)) *
     >        (sin (deli)) ** 2
      psurf = 0.5 * g * em2 * cp + 1.
      if (psurf < pmin) psurf = pmin
      pcone = psurf
      dpsk  = pcone - (2. * g * emns - gm1) / gp1
      go to 50

c     Combined surface pressure:

   40 psurf = 0.5 *g * em2 * (cpcone + cpstag + cpx) + 1.

      if (psurf < pmin) psurf = pmin

c     Shock angle dsr, relative to axis.  Use impact methods to obtain dpsk.
c     (Adjust to alpha = zero results.)

      deli  = dcr
      is    = 1
      fdp   = 1.
      emc   = ema
      if (emc <= 0.) emc = 0.025 * em
      go to 46

c     Calculate value from ncone:

   45 continue
      emns  = (ema * sin (dsc))**2
      ps    = (2. * g * emns - gm1) / gp1
      fdp   = ps / p2
      emc   = emp
      if (emc <= 0.) emc = 0.025 * em

c     Next, calculate value at actual impact angle:

      deli  = alfp + dcr
      dpsk  = 0.
      if (deli <= 0.) go to 50

      is    = 2
   46 emd   = emc * sin (deli)
      ck    = (2. * gp1 / (g + 3.)) * emd
      emns  = (ck + exp (-ck))**2
      p2    = (2. * g * emns - gm1) / gp1
      if (is == 1) go to 45

      dpsk  = psurf - p2 * fdp

   50 continue

      cp = (psurf - 1.) / (0.5 * g * em**2)

c     Calculate local flow properties?

      if (iflow /= 0) then

         if (dpsk < 0.) dpsk = 0.
         pshk  = psurf - dpsk
         if (pshk < 1.) pshk = 1.

         dsr    = (gp1 * pshk + gm1) / (2. * g * emp**2)
         if (dsr > 1.) then
            dsr = 1.
         else if (dsr < 0.) then
            dsr = 1./ emp**2
         end if
         dsr    = asin (sqrt (dsr))
         dsr    = dsr - alfp

         angle(3) = dsr * rad2deg ! Shock angle

c        Surface Mach:  The em2 used will be emp**2 in the EXPAND case.
c        Is this what was intended?

         emsurf = sqrt ((em2 * (gp1*pshk + gm1) - 2.*(pshk**2 - 1.)) /
     >                  (pshk* (gm1*pshk + gp1)))

         em2    = em**2
         if (emsurf <= 0.) emsurf = 0.025 * em
         bs(6)  = emsurf
         tst1   = (1. + 0.5*gm1*em2) / (1. + 0.5*gm1*emsurf**2)
         bs(1)  = fs(1) * psurf / tst1
         bs(2)  = fs(2) * psurf
         bs(3)  = fs(3) * tst1
         bs(4)  = fs(4) * sqrt (tst1)
         if (bs(3) >= av(1)) then
            bs(5) = 2.27e-8 * bs(3) * sqrt (bs(3)) / (bs(3) + 198.6)
         else
            bs(5) = av(2) * bs(3)**av(3)
         end if
         bs(7)  = bs(6) * bs(4)
         bs(8)  = bs(1) * bs(7) / bs(5)
      end if

      end subroutine acone
C-----------------------------------------------------------------------
C
      SUBROUTINE CALERF (ARG, JINT, RESULT)
C
C     This packet evaluates  erf(x),  erfc(x),  and  exp(x*x)*erfc(x)
C     for a real argument  x.  It contains three FUNCTION type
C     subprograms: ERF, ERFC, and ERFCX, and one SUBROUTINE type
C     subprogram, CALERF.  The calling statements for the primary
C     entries are:
C
C                   Y = ERF (X)
C
C                   Y = ERFC (X)
C     and
C                   Y = ERFCX (X)
C
C     The routine  CALERF  is intended for internal packet use only,
C     all computations within the packet being concentrated in this
C     routine.  The function subprograms invoke  CALERF  with the
C     statement
C
C          CALL CALERF (ARG, JINT, RESULT)
C
C     where the parameter usage is as follows:
C
C      Function                 Parameters for CALERF
C       call               ARG            JINT      Result
C
C     ERF (ARG)      ANY REAL ARGUMENT      0      ERF (ARG)
C     ERFC (ARG)     |ARG| < XBIG           1      ERFC (ARG)
C     ERFCX (ARG)    XNEG < ARG < XMAX      2      ERFCX (ARG)
C
C     The main computation evaluates near-minimax approximations
C     from "Rational Chebyshev approximations for the error function"
C     by W. J. Cody, Math. Comp., 1969, PP. 631-638.  This
C     transportable program uses rational functions that theoretically
C     approximate  erf(x)  and  erfc(x)  to at least 18 significant
C     decimal digits.  The accuracy achieved depends on the arithmetic
C     system, the compiler, the intrinsic functions, and proper
C     selection of the machine-dependent constants.
C
C-----------------------------------------------------------------------
C
C     Explanation of machine-dependent constants:
C
C     XMIN   = the smallest positive floating-point number.
C     XINF   = the largest positive finite floating-point number.
C     XNEG   = the largest negative argument acceptable to ERFCX;
C              the negative of the solution to the equation
C              2*exp(x*x) = XINF.
C     XSMALL = argument below which erf(x) may be represented by
C              2*x/sqrt(pi)  and above which  x*x  will not underflow.
C              A conservative value is the largest machine number X
C              such that   1.0 + X = 1.0   to machine precision.
C     XBIG   = largest argument acceptable to ERFC;  solution to
C              the equation:  W(x) * (1-0.5/x**2) = XMIN,  where
C              W(x) = exp(-x*x)/[x*sqrt(pi)].
C     XHUGE  = argument above which  1.0 - 1/(2*x*x) = 1.0  to
C              machine precision.  A conservative value is
C              1/[2*sqrt(XSMALL)]
C     XMAX   = largest acceptable argument to ERFCX; the minimum
C              of XINF and 1/[sqrt(pi)*XMIN].
C
C     Approximate values for some important machines are:
C
C                             XMIN       XINF        XNEG     XSMALL
C
C     CDC 7600      (S.P.)  3.13E-294   1.26E+322   -27.220  7.11E-15
C     CRAY-1        (S.P.)  4.58E-2467  5.45E+2465  -75.345  7.11E-15
C     IEEE (IBM/XT,
C       SUN, etc.)  (S.P.)  1.18E-38    3.40E+38     -9.382  5.96E-8
C     IEEE (IBM/XT,
C       SUN, etc.)  (D.P.)  2.23D-308   1.79D+308   -26.628  1.11D-16
C     IBM 195       (D.P.)  5.40D-79    7.23E+75    -13.190  1.39D-17
C     UNIVAC 1108   (D.P.)  2.78D-309   8.98D+307   -26.615  1.73D-18
C     VAX D-Format  (D.P.)  2.94D-39    1.70D+38     -9.345  1.39D-17
C     VAX G-Format  (D.P.)  5.56D-309   8.98D+307   -26.615  1.11D-16
C
C
C                             XBIG       XHUGE       XMAX
C
C     CDC 7600      (S.P.)  25.922      8.39E+6     1.80X+293
C     CRAY-1        (S.P.)  75.326      8.39E+6     5.45E+2465
C     IEEE (IBM/XT,
C       SUN, etc.)  (S.P.)   9.194      2.90E+3     4.79E+37
C     IEEE (IBM/XT,
C       SUN, etc.)  (D.P.)  26.543      6.71D+7     2.53D+307
C     IBM 195       (D.P.)  13.306      1.90D+8     7.23E+75
C     UNIVAC 1108   (D.P.)  26.582      5.37D+8     8.98D+307
C     VAX D-Format  (D.P.)   9.269      1.90D+8     1.70D+38
C     VAX G-Format  (D.P.)  26.569      6.71D+7     8.98D+307
C
C-----------------------------------------------------------------------
C
C     Error returns:
C
C     The program returns  ERFC = 0      for  ARG >= XBIG;
C
C                          ERFCX = XINF  for  ARG <  XNEG;
C        and
C                          ERFCX = 0     for  ARG >= XMAX.
C
C
C     Author: W. J. Cody
C             Mathematics and Computer Science Division
C             Argonne National Laboratory
C             Argonne, IL 60439
C
C     Latest modification by the author: March 19, 1990
C
C     Adjustments by D.A.Saunders, ELORET/NASA Ames, 02/02/01:
C
C        Source code was obtained from www.netlib.org/specfun/erf and
C        stripped of the CS/CD comments.  The REAL form here can be
C        compiled with the appropriate switch to give 32- or 64-bit
C        precision.  Machine-dependent constants are now supplied by
C        Fortran 90 intrinsics.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C-----------------------------------------------------------------------
C     Arguments:
C-----------------------------------------------------------------------

      REAL,    INTENT (IN)  ::
     >   ARG

      INTEGER, INTENT (IN)  ::
     >   JINT

      REAL,    INTENT (OUT) ::
     >   RESULT

C-----------------------------------------------------------------------
C     Local variables:
C-----------------------------------------------------------------------

      INTEGER
     >   I

      REAL
     >   DEL, FOUR, HALF, ONE, SIXTEN, SQRPI, TWO, THRESH,
     >   X, XDEN, XNUM, Y, YSQ, ZERO

      REAL
     >   A(5), B(4), C(9), D(8), P(6), Q(5)

      REAL, SAVE ::
     >   XBIG, XHUGE, XINF, XMAX, XMIN, XNEG, XSMALL

      LOGICAL, SAVE ::
     >   FIRST = .TRUE.

C     Extra variables needed by the Newton iteration for XBIG:

      INTEGER
     >   ITER, ITMAX
      REAL
     >   ALPHA, DV, F, FNORM, FNORM0, ROOTPI, T1, T2, TOL, V, VLAST, VSQ
      LOGICAL
     >   CONVGD, FAIL

C-----------------------------------------------------------------------
C     Mathematical constants:
C-----------------------------------------------------------------------

      DATA
     >   FOUR, ONE, HALF, TWO, ZERO /4.0E0, 1.0E0, 0.5E0, 2.0E0, 0.0E0/,
     >   SQRPI /5.6418958354775628695E-1/, THRESH /0.46875E0/,
     >   SIXTEN /16.0E0/

C-----------------------------------------------------------------------
C     Coefficients for approximation to  erf  in first interval:
C-----------------------------------------------------------------------

      DATA
     >   A /3.16112374387056560E00, 1.13864154151050156E02,
     >      3.77485237685302021E02, 3.20937758913846947E03,
     >      1.85777706184603153E-1/,
     >   B /2.36012909523441209E01, 2.44024637934444173E02,
     >      1.28261652607737228E03, 2.84423683343917062E03/

C-----------------------------------------------------------------------
C     Coefficients for approximation to  erfc  in second interval:
C-----------------------------------------------------------------------

      DATA
     >   C /5.64188496988670089E-1, 8.88314979438837594E0,
     >      6.61191906371416295E01, 2.98635138197400131E02,
     >      8.81952221241769090E02, 1.71204761263407058E03,
     >      2.05107837782607147E03, 1.23033935479799725E03,
     >      2.15311535474403846E-8/,
     >   D /1.57449261107098347E01, 1.17693950891312499E02,
     >      5.37181101862009858E02, 1.62138957456669019E03,
     >      3.29079923573345963E03, 4.36261909014324716E03,
     >      3.43936767414372164E03, 1.23033935480374942E03/

C-----------------------------------------------------------------------
C     Coefficients for approximation to  erfc  in third interval:
C-----------------------------------------------------------------------

      DATA
     >   P /3.05326634961232344E-1, 3.60344899949804439E-1,
     >      1.25781726111229246E-1, 1.60837851487422766E-2,
     >      6.58749161529837803E-4, 1.63153871373020978E-2/,
     >   Q /2.56852019228982242E00, 1.87295284992346047E00,
     >      5.27905102951428412E-1, 6.05183413124413191E-2,
     >      2.33520497626869185E-3/

C-----------------------------------------------------------------------
C     Execution:
C-----------------------------------------------------------------------

      IF (FIRST) THEN

C        Calculate precision-dependent quantities
C        (a minor price to pay for portability):

         FIRST  = .FALSE.
         XMIN   =  TINY (XMIN)
         XINF   =  HUGE (XINF)
         XNEG   = -SQRT (LOG (HALF * XINF))
         XSMALL =  EPSILON (XSMALL)
         XHUGE  =  HALF / SQRT (XSMALL)
C********XMAX   =  MIN (XINF, ONE / (XMIN * 1.772454E0))
         XDEN   =  0.8 * (XMIN * XINF) ! 0.8 < .5 sqrt(pi) avoids XMAX overflow
         XMAX   =  MIN (XINF, XDEN / XMIN)

C        Newton iteration for XBIG:  Solve the following for v:
C           f(v) = (exp(-v^2)) (v^2 - 0.5) / sqrt(pi) - XMIN v^3 = 0
C        where multiplying throughout by XINF helps avoid extremes, and we
C        use the existing variable SQRPI = 1/sqrt(pi).

         T1     = XMIN * XINF     ! O(1.0)
         ALPHA  = ZERO            ! For printout purposes only
         FNORM0 = 1.E+30          ! I.e., big to avoid iteration 0 test
         V      = -XNEG           ! Starting guess suggested by the above tables
         TOL    = 1.E-5 * ABS (V) ! No need for high precision, but must scale
         ITMAX  = 10              ! 32-bit arithmetic can hit the limit with
         CONVGD = .FALSE.         !    |f| stuck at ~1.E-3, but result is OK
         FAIL   = .FALSE.         ! For step-halving inner iteration

         DO ITER = 0, ITMAX

C           Evaluate f(v), the magnitude of which should converge to ~0.

            VSQ   = V * V
            T2    = SQRPI * XINF * EXP (-VSQ)
            F     = T2 * (VSQ - HALF) - V * VSQ * T1
            FNORM = ABS (F)

C           Halve the step until |f| is reduced (except first time through).

            DO WHILE (FNORM > FNORM0)
               IF (ALPHA > TWO * XSMALL) THEN ! Assume XSMALL = machine epsilon
                  ALPHA = HALF * ALPHA
                  V = VLAST - ALPHA * DV
                  VSQ   = V * V
                  T2    = SQRPI * XINF * EXP (-VSQ)
                  F     = T2 * (VSQ - HALF) - V * VSQ * T1
                  FNORM = ABS (F)
               ELSE
                  WRITE (*, '(A)') ' CALERF: Step halving failed.'
                  FAIL = .TRUE. ! To the cease outer iteration
                  EXIT
               END IF
            END DO

C****       WRITE (*, '(A, I3, A, ES9.2, A, ES9.2, A, ES14.6)')
C****>         ' CALERF:', ITER, '  |f|:', FNORM, '  step:', ALPHA,
C****>         '  v:', V

            IF (FNORM < TOL) CONVGD = .TRUE.

            IF (CONVGD .OR. FAIL) EXIT

C           Calculate step dv = f(v) / f'(v):

            FNORM0 = FNORM
            DV     = F / (T2 * V * (3. - TWO * VSQ) - 3. * VSQ * T1)
            VLAST  = V
            V      = V - DV
            ALPHA  = ONE

         END DO ! Next iteration

         IF (.NOT. CONVGD) THEN
            IF (DV > TOL .OR. FAIL) THEN
               WRITE (*, '(A)') ' CALERF: Iteration failed. XBIG <- 9.'
               V = 9.
            END IF
         END IF

         XBIG = V - V * TWO * XSMALL ! Play safe

C****    WRITE (*, '(A, ES25.15)')
C****>   ' XMIN:  ', XMIN,
C****>   ' XINF:  ', XINF,
C****>   ' XNEG:  ', XNEG,
C****>   ' XSMALL:', XSMALL,
C****>   ' XHUGE: ', XHUGE,
C****>   ' XMAX:  ', XMAX,
C****>   ' XBIG:  ', XBIG

      END IF

      X = ARG
      Y = ABS (X)

      IF (Y <= THRESH) THEN ! Evaluate  erf  for  |X| <= 0.46875:

         YSQ = ZERO
         IF (Y > XSMALL) YSQ = Y * Y
         XNUM = A(5)*YSQ
         XDEN = YSQ
         DO I = 1, 3
            XNUM = (XNUM + A(I)) * YSQ
            XDEN = (XDEN + B(I)) * YSQ
         END DO
         RESULT = X * (XNUM + A(4)) / (XDEN + B(4))
         IF (JINT /= 0) RESULT = ONE - RESULT
         IF (JINT == 2) RESULT = EXP (YSQ) * RESULT
         GO TO 800

      ELSE IF (Y <= FOUR) THEN ! Evaluate  erfc  for 0.46875 <= |X| <= 4.0:

         XNUM = C(9)*Y
         XDEN = Y
         DO I = 1, 7
            XNUM = (XNUM + C(I)) * Y
            XDEN = (XDEN + D(I)) * Y
         END DO
         RESULT = (XNUM + C(8)) / (XDEN + D(8))
         IF (JINT /= 2) THEN
            YSQ = AINT (Y*SIXTEN) / SIXTEN
            DEL = (Y - YSQ)*(Y + YSQ)
            RESULT = EXP (-YSQ*YSQ) * EXP (-DEL) * RESULT
         END IF

      ELSE ! Evaluate  erfc  for |X| > 4.0:

         RESULT = ZERO
         IF (Y >= XBIG) THEN
            IF ((JINT /= 2) .OR. (Y >= XMAX)) GO TO 300
            IF (Y >= XHUGE) THEN
               RESULT = SQRPI / Y
               GO TO 300
            END IF
         END IF
         YSQ = ONE / (Y * Y)
         XNUM = P(6)*YSQ
         XDEN = YSQ
         DO I = 1, 4
            XNUM = (XNUM + P(I)) * YSQ
            XDEN = (XDEN + Q(I)) * YSQ
         END DO
         RESULT = YSQ *(XNUM + P(5)) / (XDEN + Q(5))
         RESULT = (SQRPI -  RESULT) / Y
         IF (JINT /= 2) THEN
            YSQ = AINT (Y*SIXTEN) / SIXTEN
            DEL = (Y - YSQ)*(Y + YSQ)
            RESULT = EXP (-YSQ*YSQ) * EXP (-DEL) * RESULT
         END IF

      END IF

C     Fix up for negative argument, erf, etc.:

  300 IF (JINT == 0) THEN
         RESULT = (HALF - RESULT) + HALF
         IF (X < ZERO) RESULT = -RESULT
      ELSE IF (JINT == 1) THEN
         IF (X < ZERO) RESULT = TWO - RESULT
      ELSE
         IF (X < ZERO) THEN
            IF (X < XNEG) THEN
               RESULT = XINF
            ELSE
               YSQ = AINT (X*SIXTEN) / SIXTEN
               DEL = (X - YSQ)*(X + YSQ)
               Y = EXP (YSQ*YSQ) * EXP (DEL)
               RESULT = (Y+Y) - RESULT
            END IF
         END IF
      END IF

  800 RETURN

      END SUBROUTINE CALERF

C-----------------------------------------------------------------------
C
      REAL FUNCTION ERF (X)
C
C     This subprogram computes approximate values for erf(x).
C     (See comments heading CALERF.)
C
C     Author/date: W. J. Cody, January 8, 1985
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

      REAL X
      REAL RESULT

      CALL CALERF (X, 0, RESULT)

      ERF = RESULT

      END FUNCTION ERF

C-----------------------------------------------------------------------
C
      REAL FUNCTION ERFC (X)
C
C     This subprogram computes approximate values for erfc(x).
C     (See comments heading CALERF.)
C
C     Author/date: W. J. Cody, January 8, 1985
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

      REAL X
      REAL RESULT

      CALL CALERF (X, 1, RESULT)

      ERFC = RESULT

      END FUNCTION ERFC

C-----------------------------------------------------------------------
C
      REAL FUNCTION ERFCX (X)
C
C     This subprogram computes approximate values for exp(x*x) * erfc(x).
C     (See comments heading CALERF.)
C
C     Author/date: W. J. Cody, March 30, 1987
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

      REAL X
      REAL RESULT

      CALL CALERF (X, 2, RESULT)

      ERFCX = RESULT

      END FUNCTION ERFCX
c-------------------------------------------------------------------------------
c
      subroutine compr (fs, av, gamma, angle, iflow, cpstag, ifirst,
     >                  cp, bs)
c
c     Using the free stream Mach number and the equivalent wedge angle,
c     this routine computes the conditions behind the shock.
c
c     COMPR is called by subroutine FORCE one surface element at a time.
c     In HABP, COMPR is also called by FLOSEP, SHKEXP, and SKINFR, but
c     these usages have been eliminated here.
c
c     History:
c
c     c. 1953   Ames Research Center  Oblique shock relationships from
c                                     "Equations, Tables, and Charts for
c                                     Compressible Flow", NACA TR 1135.
c     Feb. 2001 David Saunders        Restructuring of (part of) program
c               ELORET/NASA Ames      HABP's FORTRAN 66 implementation.
c
c-------------------------------------------------------------------------------

      implicit none

c     Arguments:
c     ----------

      real, intent (in) ::
     >   fs(8)              ! freestream properties:
                            ! 1 = density, lbm/ft^3    2 = pressure, lbf/ft^2
                            ! 3 = temperature (units?) 4 = speed of sound?
                            ! 5 = viscosity?           6 = Mach number
                            ! 7 = M * speed of sound   8 = rho * M * a / mu

      real, intent (in) ::
     >   av(3)              ! freestream gas constants:
                            ! 1 = ??? temperature (units?)
                            ! 2 = ??? viscosity (units?)
                            ! 3 = ??? (units?)

      real, intent (in) ::
     >   gamma              ! ratio of specific heats


      real, intent (inout) ::
     >   angle(3)           ! angle(1) = is not used
                            ! angle(2) = impact angle (deg, input)
                            ! angle(3) = shock angle (deg, output if iflow /= 0)

      integer, intent (in) ::
     >   iflow              ! iflow = 0 means return Cp only/no flow properties
                            ! iflow = 1 means return bs(1:8) & angle(3) as well

      real, intent (in) ::
     >   cpstag             ! Peak Cp

      integer, intent (in) ::
     >   ifirst             ! ifirst = 1 means this is the first call for the
                            ! current gamma and free-stream Mach number;
                            ! ifirst = 0 means this is not the first such call,
                            ! and the iteration within NEWTPM is suppressed

      real, intent (out) ::
     >   cp                 ! Surface element pressure coefficient

      real, intent (out) ::
     >   bs(8)              ! local properties analogous to fs(*)

c     Local constants:
c     ----------------

      real, parameter ::
     >   deg2rad = 0.0174532925199433,
     >   rad2deg = 57.2957795130823

c     Local variables:
c     ----------------

      integer ::
     >   i, ier, iprandtl, iptr(3), itemp, ismall, j

      real ::
     >   b, bs6sq, c, cossq, d, deltar, emfssq, emn, etac, g, gm1, gp1,
     >   r(3), sinsq, small

c     Execution:
c     ----------

      if (abs (angle(2)) < 0.000001) then
         cp = 0.
         if (iflow /= 0) then
            angle(3) = 0.
            bs = fs ! Elements 1:8
         end if
         return
      end if

      g   = gamma
      gp1 = g + 1.
      gm1 = g - 1.
      ier = 1      ! This will be unset if shock is not detached

      emfssq = fs(6) ** 2
      deltar = angle(2) * deg2rad

      if (abs (angle(2)) <= 2.) then ! Cubic roots aren't reliable

c        Use weak oblique shock relationship (Liepman and Roshko, p.92).
c        Sometimes angle(2) < 0.   Don't know why it's trying to compress ...
c        gjs 6-13-88

         ier = 0
         angle(3) = 1./ emfssq + (0.5*gp1*deltar) / sqrt (emfssq - 1.)

      else if (angle(2) <= 55.) then

c        Set up a cubic to be solved for sin**2 theta (wedge shock angle).
c        The coefficients are defined in equation 150b of NACA TR 1135.

         sinsq  = (sin (deltar)) ** 2
         b = -(emfssq + 2.) / emfssq - g*sinsq
         c = (2.*emfssq + 1.) / emfssq**2 +
     >       (.25*gp1**2 + gm1/emfssq)*sinsq
         d = (sinsq - 1.) / emfssq**2

         call cubic (b, c, d, r, ier)

         if (ier == 0) then

c           The shock is not detached - three real roots were found.
c           The smallest root requires a decrease in entropy - not allowed.
c           The largest root is not attained in practice.  Therefore, pick
c           the middle root.  Bubble-sort the three roots into ascending order:

            iptr(1) = 1
            iptr(2) = 2
            iptr(3) = 3

            do i = 1, 2
               small = r(iptr(i))
               do j = i + 1, 3
                  if (r(iptr(j)) < small) then
                     small  = r(iptr(j))
                     ismall = j
                  end if
               end do
               if (small /= r(iptr(i))) then
                  itemp        = iptr(i)
                  iptr(i)      = iptr(ismall)
                  iptr(ismall) = itemp
               end if
            end do

            angle(3) = r(iptr(2)) ! Square of sin of wedge shock angle

         end if

      end if

c     Assign Cp using the selected (sin (theta))**2:

      if (ier == 0) then

         emn = emfssq * angle(3)
         if (emn < 1.01) emn = 1.01
         cp  = 4.*(emn - 1.) / (gp1*emfssq)
         iprandtl = 0

      else ! No cubic solution, or delta > 55 deg. - the shock is detached.

c        Calculate flow properties by the method of Kaufman (NEWTPM).
c        This method's initial iteration is a function of just gamma and the
c        free-stream Mach, so we avoid doing it more than once per condition.
c        If that iteration fails, NEWTPM stops, so convergence is not an issue.
c        However, NEWTPM either applies Prandtl-Meyer expansion techniques (via
c        EXPAND) or it cannot.  If not, it uses the detached shock method of
c        Smyth instead.
c
c        If local flow properties are requested, NEWTPM returns them only if
c        EXPAND is called.  We avoid reassigning them below in this case.
c        Otherwise, a modified Newtonian Cp is returned and local flow
c        properties are set below.

         etac = 1. ! Empirical correction factor

         call newtpm (fs, av, g, angle, iflow, cpstag, ifirst, etac,
     >                cp, emn, bs, iprandtl)

      end if

c     Remaining flow properties?  NEWTPM returns them if requested and
c     if it applied Prandtl-Meyer expansion techniques (iprandtl = 1),
c     so avoid reassigning them in this case:

      if (iflow /= 0 .and. iprandtl == 0) then

         if (angle(3) < 0.) then
            angle(3) = 0.
         else if (angle(3) > 1.) then
            angle(3) = 1.
         end if

         angle(3) = rad2deg * asin (sqrt (angle(3)))      ! shock angle

         if (emn < 1.01) emn = 1.01
         bs(1) = (fs(1) * gp1*emn) / (gm1*emn + 2.)       ! rho, eq. 129, TR1135
         bs(2) =  fs(2) * (2.*g*emn - gm1) / gp1          ! pressure, eq. 128
         r(3)  = (2.*g*emn - gm1)*(2. + gm1*emn) / (emn*gp1**2)
         bs(3) = fs(3) * r(3)                             ! temperature, eq. 130
         bs(4) = fs(4) * sqrt (r(3))                      ! speed of sound
         if (bs(3) < av(1)) then
            bs(5) = av(2) * bs(3)**av(3)                  ! viscosity
         else
            bs(5) = 2.27e-8 * bs(3) * sqrt (bs(3)) / (bs(3) + 198.6)
         end if
         bs6sq = (gp1*gp1*emfssq*emn - 4.*(emn-1.)*(g*emn+1.)) /
     >           ((2.*g*emn - gm1) * (2. + gm1*emn))
         if (bs6sq < 1.) then
            bs(6) = 1.01
         else
            bs(6) = sqrt (bs6sq)                          ! Mach number, eq. 132
         end if
         bs(7) = bs(4) * bs(6)                            ! velocity
         bs(8) = bs(1) * bs(7) / bs(5)                    ! Reynolds number / ft

      end if

      end subroutine compr
c-------------------------------------------------------------------------------

      subroutine cone (fs, av, gamma, angle, iflow, cp, bs)

c     This is a structured version of HABP's FORTRAN 66 subroutine cone.
c
c     The subroutine solves for the flow properties about a cone using
c     a combination of second-order slender body theory and the approximate
c     cone solution of Hammitt and Murthy.
c
c     It uses second-order slender-body theory for small values of the
c     similarity parameter, the approximate solution of Hammitt and Murthy
c     for large values, and a suitable transition function for the mid-range.
c     It is called by subroutine force one surface element at a time.
c
c     History:
c
c     c. 1968   HABP Mark III:                 A.E.Gentry, McDonnell Douglas
c     c. 1990   Improved tangent-cone method:  C.I.Cruz/G.J.Sova, LaRC/Rockwell;
c                                              D.N.Smyth (LaRC) co-author?
c     Jan 2001  Literal restructuring:         D.A.Saunders, ELORET/NASA Ames
c     02/15/01  Added "iflow" argument.          "    "
c
c-------------------------------------------------------------------------------

      implicit none

c     Arguments:
c     ----------

      real, intent (in) ::
     >   fs(8)              ! freestream properties:
                            ! 1 = density, lbm/ft^3    2 = pressure, lbf/ft^2
                            ! 3 = temperature (units?) 4 = speed of sound?
                            ! 5 = viscosity?           6 = Mach number
                            ! 7 = M * speed of sound   8 = rho * M * a / mu

      real, intent (in) ::
     >   av(3)              ! freestream gas constants:
                            ! 1 = ??? temperature (units?)
                            ! 2 = ??? viscosity (units?)
                            ! 3 = ??? (units?)

      real, intent (in) ::
     >   gamma              ! ratio of specific heats

      real, intent (inout) ::
     >   angle(3)           ! angle(1) = impact angle (degrees, input)
                            ! angle(2) = |angle(1)|   (degrees, output)
                            ! angle(3) = shock angle? (degrees, output)

      integer, intent (in) ::
     >   iflow              ! iflow = 0 means return Cp only/no flow properties
                            ! iflow = 1 means return bs(1:8) & angle(3) as well

      real, intent (out) ::
     >   cp                 ! surface element pressure coefficient

      real, intent (out) ::
     >   bs(8)              ! local properties analogous to fs(*)

c     Local constants:
c     ----------------

      real, parameter ::
     >   deg2rad = 0.0174532925,
     >   rad2deg = 57.29577951

c     Local variables:
c     ----------------

      integer
     >   ksol

      real
     >   a, am, beta, beta2, cphm, dcp2, dcr, dsr, dx, dxf, em, emc,
     >   emns, emnsf, emnssq, emsin, emsin0, emsinf, emsq, emssq,
     >   g, gm1, gmsq, gp1, pcp1, phi, pspc, sindc, tandc, tansq, tct1

      integer, save ::
     >   itc = 0, ihm = 0, isb = 0, icom = 0

c     Execution:
c     ----------

      itc = itc + 1

      em       = fs(6)
      angle(2) = abs (angle(1))
      dcr      = angle(2) * deg2rad
      sindc    = sin (dcr)

c     Test for a very small equivalent cone angle:
c     --------------------------------------------

      if (sindc <= 0.0001) then

c        Set freestream properties and return:

         cp = 0.0

         if (iflow /= 0) then
            bs = fs        ! Elements 1:8
            angle(3) = asin (1.0 / em) * rad2deg
            return
         end if

      end if

      emsin = em*sindc
      emsq  = em**2
      g     = gamma
      gm1   = g - 1.0
      gp1   = g + 1.0

c     Calculate upper transition point, emsinf:
c     -----------------------------------------

      if (em >= 10.0) then
         emsinf = 1.40
      else if (em > 1.5) then
         emsinf = 1.40 - 1.075 * exp (-0.8*(em - 1.5))
      else
         emsinf = 0.325
      end if

      if (emsin >= emsinf) then

c        Cone method of Hammitt and Murthy:
c        ----------------------------------

         call h_and_m () ! See internal procedures

         if (iflow /= 0) call flow_properties ()

         return

      end if

C     Calculate lower transition point, emsin0:
c     -----------------------------------------

      if (em > 3.0) then
         emsin0 = 0.3
      else
         emsin0 = 0.2
      end if

      if (emsin <= emsin0) then

c        2nd-order slender body theory:
c        ------------------------------

         call slender_body ()

         if (iflow /= 0) then
            dsr  = asin (1.0/em)
            pcp1 = 0.5*g*emsq*cp + 1.0
            tct1 = pcp1 ** (gm1/g)
            emc  = sqrt ((2.0/gm1) * ((1.0 + 0.5*gm1*emsq) / tct1-1.0))
            if (emc <= 0.0) emc = 0.025*em
         end if

      else

c        Combination of both 2nd-order slender body and Hammitt-Murthy:
c        --------------------------------------------------------------

         ksol  = 2
         sindc = emsinf / em
         dcr   = asin (sindc)

         call h_and_m ()

         if (ksol == 2) then

            icom   = icom + 1

            cphm   = cp
            dcr    = asin (emsin0/em)

            call slender_body ()

            dcp2   = tandc * (4.0*(phi-1.0) + tansq*(6.*phi*(2.*phi-1.0)
     >               * beta**2 - beta2 * (4.0*phi-1.0) + 4.0*am)) /
     >               (em * (cos (dcr))**3)
            dx     = emsin  - emsin0
            dxf    = emsinf - emsin0
            a      = (cphm - cp - dcp2*dxf) / dxf**2
            cp     = cp + dx*(dcp2 + a*dx)

            if (iflow /= 0) then
               pcp1   = 0.5*g*emsq*cp + 1.0
               emnsf  = em * sin (dsr)
               emns   = (emnsf - 1.0) * (emsin/emsinf)**2 + 1.0
               emnssq = emns**2
               gmsq   = g * emnssq
               emssq  = (gp1*em*emns)**2 - 4.0 * (emnssq - 1.0) *
     >                  (gmsq + 1.0)
               dsr    = asin (emns / em)
               emssq  = emssq / ((2.0*gmsq - gm1)*(gm1*emnssq + 2.0))
               pspc   = (2.0*gmsq - gm1) / (gp1 * pcp1)
               emc    = sqrt (((1.0 + 0.5*gm1*emssq) * pspc**(gm1/g)-1.)
     >                        * (2.0 / gm1))
               if (emc <= 0.0) emc = 0.025 * em
               tct1   = (2.0 + gm1*emsq) / (2.0 + gm1*emc**2)
            end if

         end if

      end if

      if (iflow /= 0) call flow_properties ()

      return


c     Internal procedures for subroutine cone:
c     ----------------------------------------

      contains

         subroutine h_and_m ()

!        Hammitt and Murthy approximate cone solution
!        --------------------------------------------

         real
     >      ck, elim, h1, h2, hs2, hx, pcps, psp1, rad, sin2, sinq, t1ts

         ihm  = ihm + 1

         sinq = sindc**2
         sin2 = sin (2.0*dcr)
         h2   = gm1*sinq + 2.0/emsq
         h1   = 2.0 - (g + 5.0)*sinq
         hx   = h2*h1/sin2**2

         if (hx <= -1.0) then ! Detached flow: use tangent-cone empirical
            ck   = 2.0 * emsin * gp1 / (g + 3.0)
            dsr  = asin (min (1.0, (ck + exp (-ck))/em))
            ksol = 1
         else if (abs (h1) > 0.001) then
            rad  = sqrt (1.0 + hx)
            dsr  = dcr - sin2 * (1.0 - rad) / h1
         else
            dsr  = dcr + 0.5*h2 * (1.0 - 0.25*hx*(1.0 - 0.5*hx)) / sin2
         end if

         emnssq = emsq * (sin (dsr))**2
         hs2  = 2.0 * (dsr - dcr)**2
         emc  = sqrt ((emsq - emnssq)*(1.0 + hs2) /
     >                (1.0 + 0.5*gm1 * (emnssq*(1.0 + hs2) - hs2*emsq)))
         if (emc <= 0.0) emc = 0.025*em

         elim = gm1 / gp1
         psp1 =  2.0 * (g / gp1) * emnssq - elim
         t1ts = (psp1 + elim) / ((1.0 + elim*psp1) * psp1)
         tct1 = (2.0 + gm1*emsq) / (2.0 + gm1*emc**2)
         pcps = (tct1*t1ts) ** (g/gm1)
         pcp1 = pcps * psp1
         cp   = 2.0 * (pcp1 - 1.0) / (g * emsq)

         end subroutine h_and_m

         subroutine slender_body ()

!        2nd-order slender body theory
!        -----------------------------

         isb   = isb + 1

         beta  = sqrt (emsq - 1.0)
         tandc = tan (dcr)
         tansq = tandc**2
         beta2 = 5.00*emsq - 1.0
         am    = 3.25*emsq + 0.5 + gp1*(emsq/beta)**2
         phi   = 0.69314718 - log (beta * tandc)
         cp    = tansq * (2.0*phi - 1.0 +
     >           tansq * (3.0*(beta*phi)**2 - beta2*phi + am))

         end subroutine slender_body

         subroutine flow_properties ()

!        Local flow properties output along with Cp by cone routine:
!        -----------------------------------------------------------

         angle(3) = dsr * rad2deg
         bs(1)    = fs(1) * pcp1 / tct1
         bs(2)    = fs(2) * pcp1
         bs(3)    = fs(3) * tct1
         bs(4)    = fs(4) * sqrt (tct1)
         if (bs(3) < av(1)) then
            bs(5) = av(2) * bs(3)**av(3)
         else
            bs(5) = 2.27e-8 * bs(3) * sqrt (bs(3)) / (bs(3) + 198.6)
         end if
         bs(6)    = emc
         bs(7)    = bs(6) * bs(4)
         bs(8)    = bs(1) * bs(7) / bs(5)

         end subroutine flow_properties

      end subroutine cone
************************************************************************
*
      SUBROUTINE CUBIC (A, B, C, ROOT, IER)
*
*     Calculate the roots of the cubic equation x^3 + ax^2 + bx + c = 0
*     for real coefficients a, b, c.
*
*     The method is that from Numerical Recipes:
*
*     Letting  Q = (a^2 - 3b) / 9  and  R = (2a^3 - 9ab + 27c) / 54,
*     then if  R^2 < Q^3 there are three real roots.
*
*     Note that this excludes Q = 0, so degenerate cases such as
*     (X - 1)^3 = 0 are evidently not handled.
*
*     Other cases are not treated here yet.
*
*     02/13/01  David Saunders  Initial implementation (real roots only).
*               ELORET/NASA Ames
*
************************************************************************

      IMPLICIT NONE

*     Arguments:

      REAL,    INTENT (IN) ::
     >   A, B, C             ! Cubic coefficients

      REAL,    INTENT (OUT) ::
     >   ROOT(3)             ! Three real roots (undefined if IER /= 0)

      INTEGER, INTENT (OUT) ::
     >   IER                 ! 0 means three real roots were found;
                             ! 1 means roots are complex (or possibly
                             ! not distinct) and are not returned

*     Local constants:

      REAL, PARAMETER ::
     >   THIRD = 1./3.,
     >   TWOPI = 6.28318530717958647692528

*     Local variables:

      REAL
     >   Q, R, ROOTQ, THETA

*     Execution:

      Q = (A * A  - 3.* B) / 9.
      R = (A * (2.* A * A - 9.* B) + 27.* C) / 54.

      IF (Q > 0. .AND. R ** 2 < Q ** 3) THEN
         ROOTQ = SQRT (Q)
         THETA = ACOS (R / (Q * ROOTQ))
         ROOTQ = -2.* ROOTQ

         ROOT(1) = ROOTQ * COS ( THETA          * THIRD) - (A * THIRD)
         ROOT(2) = ROOTQ * COS ((THETA + TWOPI) * THIRD) - (A * THIRD)
         ROOT(3) = ROOTQ * COS ((THETA - TWOPI) * THIRD) - (A * THIRD)

         IER = 0
      ELSE
         IER = 1
      END IF

      END SUBROUTINE CUBIC
c-------------------------------------------------------------------------------
c
      subroutine expand (fs, av, gamma, angle, isdet, mer, cp, bs)
c
c     Given the free stream conditions (fs) and the turning angle in degrees
c     (angle(2)), EXPAND performs an isentropic Prandtl-Meyer expansion
c     (angle(2) > 0) or compression (angle(2) < 0).
c
c     History:
c
c     Pre-1989   NASA Langley      Original code from program HABP (F66).
c     Feb. 2001  David Saunders    Restructured version (F90);
c                ELORET/NASA Ames  M = F(M) iteration replaced with Newton itn.
c
c-------------------------------------------------------------------------------

      implicit none

c     Arguments:
c     ----------

      real, intent (in) ::
     >   fs(8)              ! freestream properties:
                            ! 1 = density, lbm/ft^3    2 = pressure, lbf/ft^2
                            ! 3 = temperature (units?) 4 = speed of sound?
                            ! 5 = viscosity?           6 = Mach number
                            ! 7 = M * speed of sound   8 = rho * M * a / mu

      real, intent (in) ::
     >   av(3)              ! freestream gas constants:
                            ! 1 = ??? temperature (units?)
                            ! 2 = ??? viscosity (units?)
                            ! 3 = ??? (units?)

      real, intent (in) ::
     >   gamma              ! ratio of specific heats

      real, intent (inout) ::
     >   angle(3)           ! angle(1) = is not used
                            ! angle(2) = impact angle (deg, input)
                            ! angle(3) = shock angle (deg, output if iflow /= 0)

      integer, intent (in) ::
     >   isdet              ! isdet = 2 means just Cp and bs(6) are returned,
                            ! else bs(1:8) are returned but not Cp (!)

      integer, intent (out)::
     >   mer                ! mer = 2 if flow compressed to subsonic or
                            ! expanded to infinite Mach - Cp is not assigned;
                            ! otherwise, mer is not assigned

      real, intent (out) ::
     >   cp                 ! Surface element pressure coefficient; see isdet
                            ! and mer
      real, intent (out) ::
     >   bs(8)              ! Local properties analogous to fs(*); see isdet

c     Local constants:

      real, parameter ::
     >   deg2rad = 0.0174532925199433,
     >   rad2deg = 57.2957795130823

c     Local variables:

      integer
     >   iter
      real
     >   alpha, dm, emsq, f, fnorm, fnorm0, g, gm1, gp1, gr, m, mlast,
     >   nu1, nu1d, nu2, nu2d, rg, rootmsq, rtmsq, small, term, tol, z
      logical
     >   convgd, fail

c     Execution:

      g   = gamma
      gp1 = g + 1.
      gm1 = g - 1.
      gr  = sqrt (gp1 / gm1)
      rg  = 1./ gr

      if (fs(6) < 1.) then ! Input Mach is subsonic.
         emsq = 1.         ! For program continuity.
      else
         emsq = fs(6)**2
      end if

c     Calculate Prandtl-Meyer angle for free stream conditions using
c     equation 171c of TR 1135 (radians):

      rootmsq = sqrt (emsq - 1.)
      nu1     = gr * atan (rootmsq * rg) - atan (rootmsq)
      nu1d    = nu1  * rad2deg
      nu2d    = nu1d + angle(2)
      nu2     = nu2d * deg2rad ! Prandtl-Meyer angle after the expansion

c     Check if flow compressed to subsonic:

      if (nu2d <= 0.) then ! Return sonic conditions
         bs(6) = 1.
         mer   = 2
         return
      end if

c     Check if flow has expanded to an infinite Mach number,
c     taken as 100. for all practical purposes:

      if (nu2d >= 0.97815 * (gr - 1.) * 90.) then
         bs(6) = 100.
         mer   = 2
         return
      end if

c     Safeguarded Newton iteration to find Mach number downstream
c     corresp. to nu2, via solution of Eq. 171c (TR 1135) for M:

c     Starting guess (explain?):

      m = fs(6) * (1. + (nu2 - nu1)*(1. + .5*gm1*emsq) / rootmsq)
      if (m <= 1.0) m = 1.01

      alpha  = 0.              ! For printout purposes only
      fnorm0 = 1.e+30          ! I.e., big to avoid iteration 0 test
      tol    = 1.e-5 * abs (m) ! No need for high precision, but must scale
      small  = 2.* epsilon (small)
      convgd = .false.
      fail   = .false.         ! For step-halving inner iteration

      do iter = 0, 20

c        Evaluate f(m), the magnitude of which should converge to ~0.

         rtmsq = sqrt (m**2 - 1.)
         term  = rtmsq * rg
         f     = gr * atan (term) - atan (rtmsq) - nu2
         fnorm = abs (f)

c        Halve the step until |f| is reduced (except first time through):

         do while (fnorm > fnorm0)
            if (alpha > small) then
               alpha = .5 * alpha
               m     = mlast - alpha * dm
               rtmsq = sqrt (m**2 - 1.)
               term  = rtmsq * rg
               f     = gr * atan (term) - atan (rtmsq) - nu2
               fnorm = abs (f)
            else
               write (*, '(a)') ' EXPAND: Step halving failed.'
               fail  = .true.
               exit
            end if
         end do

*****    WRITE (*, '(A, I3, A, ES9.2, A, ES9.2, A, ES14.6)')
*****>      ' EXPAND:', ITER, '  |f|:', FNORM, '  step:', ALPHA,
*****>      '  M:', M

         if (fnorm < tol) convgd = .true.

         if (convgd .or. fail) exit

C        Calculate step dv = f(v) / f'(v):

         fnorm0 = fnorm
         dm     = f * rtmsq / (m / (1. + term**2) - 1./ m)
         mlast  = m
         m      = m - dm
         alpha  = 1.

      end do ! Next iteration

      if (.not. convgd) then
         if (dm > tol .or. fail) then
            write (*, *) ' EXPAND: Unconverged; using M = ', m
         end if
      end if

      bs(6) = m

c     Calculate final characteristics behind expansion fan:

      angle(3) = rad2deg * asin (1./ m)
      z = (2. + gm1 * emsq) / (2. + gm1 * m**2)

      if (isdet == 2) then
         cp = (z**(g / gm1) - 1.) / (.5*g*emsq)           ! Only
      else
         bs(1) = fs(1) * z**(1./ gm1)
         bs(2) = fs(2) * z**(g / gm1)
         bs(3) = fs(3) * z
         bs(4) = fs(4) * sqrt (z)
         if (bs(3) < av(1)) then
            bs(5) = av(2) * bs(3)**av(3)                  ! viscosity
         else
            bs(5) = 2.27e-8 * bs(3) * sqrt (bs(3)) / (bs(3) + 198.6)
         end if
         bs(7) = bs(4) * m                                ! velocity
         bs(8) = bs(1) * bs(7) / bs(5)                    ! Reynolds number / ft
      end if

      end subroutine expand
c-------------------------------------------------------------------------------
c
      subroutine newtpm (fs, av, gamma, angle, iflow, cpstag, ifirst,
     >                   etac, cp, emn, bs, iprandtl)
c
c     NEWTPM (Newtonian + Prandtl-Meyer) calculates surface conditions using
c     the blunt body shock-expansion technique of Kaufman, Journal of the
c     Astronautical Sciences, Vol X, No.2, Summer 1963.
c
c     Newtonian theory is used along the body until a point is reached where
c     the pressure and pressure gradient both match what would be calculated
c     by a continuing Prandtl-Meyer expansion.
c
c     History:
c
c     Pre-1989  NASA Langley          Original code from program HABP (F66).
c     Feb. 2001 David Saunders        Restructured version (F90).
c               ELORET/NASA Ames
c
c-------------------------------------------------------------------------------

      implicit none

c     Arguments:
c     ----------

      real, intent (inout) ::                      ! Mach may be temporarily set
     >   fs(8)              ! freestream properties:
                            ! 1 = density, lbm/ft^3    2 = pressure, lbf/ft^2
                            ! 3 = temperature (units?) 4 = speed of sound?
                            ! 5 = viscosity?           6 = Mach number
                            ! 7 = M * speed of sound   8 = rho * M * a / mu

      real, intent (in) ::
     >   av(3)              ! freestream gas constants:
                            ! 1 = ??? temperature (units?)
                            ! 2 = ??? viscosity (units?)
                            ! 3 = ??? (units?)

      real, intent (in) ::
     >   gamma              ! ratio of specific heats


      real, intent (inout) ::
     >   angle(3)           ! angle(1) = is not used
                            ! angle(2) = impact angle (deg, input)
                            ! angle(3) = shock angle (deg, output if iflow /= 0)

      integer, intent (in) ::
     >   iflow              ! iflow = 0 means return Cp only/no flow properties
                            ! iflow = 1 means return bs(1:8) & angle(3) as well
      real, intent (in) ::
     >   cpstag             ! Peak Cp

      integer, intent (in) ::
     >   ifirst             ! ifirst = 1 means this is the first call for the
                            ! current gamma and free-stream Mach number;
                            ! ifirst = 0 means this is not the first such call,
                            ! and the initial iteration here is suppressed

      real, intent (in) ::
     >   etac               ! Empirical correction factor applied to Cp

      real, intent (out) ::
     >   cp,                ! Surface element pressure coefficient
     >   emn                ! Normal Mach squared times square of shock angle

      real, intent (out) ::
     >   bs(8)              ! Local properties analogous to fs(*); see iprandtl

      integer, intent (out) ::
     >   iprandtl           ! Indicates which method was used here:
                            ! iprandtl = 1 means EXPAND was called; if flow
                            ! properties are requested, they are assigned here;
                            ! iprandtl = 0 means EXPAND was not called; a
                            ! Newtonian Cp is returned along with angle(3)
                            ! and emn (only, regardless of iflow)

c     Local constants:

      integer, parameter ::
     >   maxiter = 20

      real,    parameter ::
     >   gascp   = 6007.93, ! Sp. heat for air @ constant pr. (ft/sec)^2 deg R
     >   deg2rad = 0.0174532925199433,
     >   rad2deg = 57.2957795130823

c     Local variables:

      integer
     >   isdet, iter, mer

      real
     >   a, cp_unused, deltmu, emup, emlow, epsi, g, gm1, gp1,
     >   m1, m2, mach, machsq, mu, p, p1, p2, pc, ppfs, ppo, q,
     >   re, rho, t, tsubt, v

      real, save ::
     >   deltq, macho, mq, pcap

c     Execution:

      g   = gamma
      gm1 = g - 1.
      gp1 = g + 1.

      if (ifirst /= 0) then ! First time for this gamma & free-stream Mach:

         emlow = 0.91 + 0.3125*g
         emup  = emlow + 0.40
         macho = fs(6)

c        Free-stream static to stag. pressure ratio (1.0 / eq. 100, TR 1135):

         pcap  = (2./ (gp1*macho**2)) ** (g/gm1) *
     >           ((2.*g*macho**2 - gm1) / gp1) ** (1.0/gm1)

c        Set up a secant iteration to find the matching point Mach number.
c        The following revision avoids GO TOs at the expense of repeated code.

         m1 = emlow
         machsq = m1**2
         q  = (2./ (2. + gm1*machsq)) ** (g/gm1)  ! psubq / po; see Kaufman
         p1 = q*(1. - (q*(g*machsq)**2) / (4.*(machsq - 1.)*(1. - q)))

         m2 = emup
         machsq = m2**2
         q  = (2./ (2. + gm1*machsq)) ** (g/gm1)
         p2 = q*(1. - (q*(g*machsq)**2) / (4.*(machsq - 1.)*(1. - q)))

         do iter = 1, maxiter

            if (abs (p2 - p1) < 0.000001) then ! p vs. M is ~flat
               exit
            else ! Linear interpolation for Mach that gives pc = pcap
               mq = m1 + (pcap - p1)*(m2 - m1)/(p2 - p1)

               if (abs (mq - m2) < 0.0001) then ! Converged
                  exit
               else
                  if (mq > emup)  mq = emup
                  if (mq < emlow) mq = emlow

                  machsq = mq**2
                  q  = (2./ (2. + gm1*machsq)) ** (g/gm1)
                  pc = q*(1.- (q*(g*machsq)**2)/(4.*(machsq-1.)*(1.-q)))
                  p1 = p2
                  p2 = pc
                  m1 = m2
                  m2 = mq
               end if
            end if
         end do ! Next iteration

         if (iter >= maxiter) then
            write (*, 50) m1, m2, mq, emlow, emup, q, p1, p2, pc
   50       format (' *** NEWTPM: Matching pt. Mach iteration failed.',
     >              /, ' m1, m2, mq, emlow, emup: ', 5e12.3,
     >              /, ' q, p1, p2, pc:           ', 4e12.3)
            stop
         end if

         deltq = sqrt ((q - pcap) / (1. - pcap)) ! Matching pt. impact angle
         deltq = rad2deg * asin (deltq)          ! No need to force [0.,1.]?

      end if ! ifirst = 1 block


      deltmu = deltq - angle(2) ! Matching pt. expansion angle

c     Check if flow will expand at least to matching point:

      if (deltmu >= 0.) then

         fs(6) = mq       ! Temporarily set effective free-stream Mach number
         angle(2) = deltmu
         isdet = 0        ! Requests bs(*); Cp is not used.  (Explain?)

         call expand (fs, av, g, angle, isdet, mer, cp_unused, bs)

         fs(6) = macho    ! Restore it
         mach  = bs(6)    ! Surface mach number

c        Calculate surface pressure ratio (eq. 44 of tr 1135):

         ppo   = etac * (1. + .5*gm1*mach**2) ** (-g/gm1)
         ppfs  = (1./ pcap)*ppo
         cp    = (2./ (g*macho**2))*(ppfs - 1.)

         if (iflow /= 0) then

            iprandtl = 1  ! Suppresses setting of local properties in COMPR

            tsubt = fs(3) * (1. + .5*gm1*macho**2) ! Total temperature (eq. 43)
            t     = tsubt / (1. + .5*gm1*mach**2)  ! Temp. after expansion (d.R)
            p     = ppfs*fs(2)                     ! Surface pressure
            rho   = (g*p) / (gm1*gascp*t)          ! Density (eq. 26 of tr 1135)
            a     = sqrt (g*p/rho)                 ! Local speed of sound
            v     = mach*a                         ! Local velocity
            if (t < av(1)) then
               mu = av(2) * t**av(3)               ! Viscosity
            else
               mu = 2.27e-8 * t * sqrt (t) / (t + 198.6)
            end if
            re = rho*v/mu                          ! Reynolds # / ft

            bs(1) = rho
            bs(2) = p
            bs(3) = t
            bs(4) = a
            bs(5) = mu
            bs(6) = mach
            bs(7) = v
            bs(8) = re

         end if

      else

c        Flow has not reached the matching point.
c        Use Newtonian calculations and detached shock method of Smyth:

         cp = cpstag * (sin (deg2rad * angle(2))) ** 2

         if (iflow /= 0) then ! Effective flow values for detached conditions

            iprandtl = 0 ! Tells COMPR to fill in missing local properties

c           Square of Mach number normal to effective shock:

            machsq   = (0.5*cp * macho**2 - (2./gp1)) / (1. - gm1/gp1)

            epsi     = (gm1/gp1) * (1. + 2./(gm1*machsq)) ! Density ratio
            angle(3) = (0.5*cp)  / (1. - epsi)            ! Shock angle squared
            emn      = machsq * angle(3)                  ! Normal Mach sqrd. *
                                                          ! sqr. of shock angle
         end if

      end if

      end subroutine newtpm
C+----------------------------------------------------------------------
C
      SUBROUTINE SECOND (CPUSEC)
C
C PURPOSE:
C     Measures CPU time in same form as the CRAY utility.
C
C ARGUMENTS:
C    ARG      DIM  TYPE I/O/S DESCRIPTION
C    CPUSEC    -    R     O   CPU time used so far, in seconds.
C
C METHOD:
C     IRIS FORTRAN requires use of the C language's intrinsic clock
C     function, which returns CPU time used since the FIRST call to it.
C     However, because FORTRAN cannot access the intrinsic, we call a
C     C routine (iclock) that in turns calls clock.
C     Conversion to seconds is done here.
C
C KNOWN BUG:
C     The count of microseconds wraps around after about 36 minutes
C     for the original iclock form => the mclock form should be better.
C
C USAGE:
C        CALL SECOND (TIME1)
C        ::::::::::::::::::
C        CALL SECOND (TIME2)
C        TOTALT = TIME2 - TIME1
C
C ENVIRONMENT:  SGI IRIS, FORTRAN 77
C
C HISTORY:
C     08/31/82   Dan McKernan    VAX/VMS version.
C     05/30/90   David Saunders  IRIS version (iclock).
C     06/05/90   Dexter Hermstad IRIS version (continued).
C     03/21/97   D.Saunders      IRIS mclock form obtained from J.Reuther.
C-----------------------------------------------------------------------

C     Arguments:

      REAL CPUSEC

C     C routine to call a system utility:

C*****INTEGER iclock
      INTEGER mclock

C     Execution:

C*****CPUSEC = FLOAT (iclock()) * 1.E-6
      CPUSEC = FLOAT (mclock()) * 0.01

      RETURN
      END
