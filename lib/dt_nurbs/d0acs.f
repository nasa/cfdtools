      SUBROUTINE D0ACS (CMEM, IMEM, DMEM, NDEG, TOL, IPBF, ISBF, ICBF,
     +                  IER)

C     Approximate, by a spline space curve, a curve on a surface originally
C     given by the composition of a parameter space curve and the surface 
C     function.
C     
C     This version is a quick temporary version to be used until a proper
C     D2ACS subroutine can be built.   A proper D2ACS will find a proper set
C     of knots and refine the curve to meet the requested tolerance.  This
C     version just heuristically picks a set of points and interpolates them
C     without checking tolerance.
C
C     IER     Returned error code.  Zero implies no errors.
C             -1  = Dynamic memory is corupt or uninitialized
C             -2  = IPBF is a null pointer
C             -3  = IPBF is not a valid pointer
C             -4  = Ambiguous - IPBF is either not valid or deleted
C             -5  = IPBF points to a deleted entity
C             -6  = IPBF does not point to a B-spline Function entity
C             -7  = The spline array in IPBF is invalid
C             -8  = IPBF is not a curve
C             -9  = IPBF does not map into a 2D parameter space
C             -10 = ISBF is a null pointer
C             -11 = ISBF is not a valid pointer
C             -12 = Ambiguous - ISBF is either not valid or deleted
C             -13 = ISBF points to a deleted entity
C             -14 = ISBF does not point to a B-spline Function entity
C             -15 = The spline array in ISBF is invalid
C             -16 = ISBF is not a surface
C             -17 = ISBF does not map into 3D space
C             -18 = Unable to allocate space for ICBF entity
C             -19 = Unable to allocate sufficient workspace in DMEM
C             -20 = Unexpected error from DTSPVL
C             -21 = Unexpected error from DTNPVL
C             -22 = Unexpected error from DTNSI while doing x-coordinates
C             -23 = Unexpected error from DTNSI while doing y-coordinates
C             -24 = Unexpected error from DTNSI while doing z-coordinates 
C             -25 = Unexpected error from D2ERAS while erasing workspaces
C             
      CHARACTER CMEM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)
      INTEGER NDEG, IPBF, ISBF, ICBF, IER
      DOUBLE PRECISION TOL

      INTEGER ENTBF
      PARAMETER (ENTBF=246)
      INTEGER ENTDA
      PARAMETER (ENTDA=253)

      INTEGER IERX, NP, MRP, MP, KP, NCP, NS, MRS, MS, KS(2), NCS(2)
      INTEGER IPH, ISH, I, IPD, ISD, KC, NCC, ICD, IWDA, IWD
      INTEGER NWORK, JCC, LENC, LC, ITDA, ITD, NPT
      DOUBLE PRECISION PLO, PHI, SLO(2), SHI(2), UV(2), XYZ(3)
C****
      ITDA = 0
      IWDA = 0
      ICBF = 0
C     CALL DTERPT (0)

C     Check IPBF validity
      CALL D0PTR (CMEM, IMEM, DMEM, IPBF, ENTBF, 0, IPH, I, IER)
      IF (IER .NE. 0) GOTO 9900
      IPD = IMEM(IPH+3)
      CALL DTGET (DMEM(IPD), .TRUE., 1, NP, MRP, MP, KP, NCP, PLO, PHI,
     +            IERX)
      IF (IERX .NE. 0) THEN
         IER = -7
         GOTO 9900
      ELSEIF (NP .NE. 1) THEN
         IER = -8
         GOTO 9900
      ELSEIF (MP .NE. 2) THEN
         IER = -9
         GOTO 9900
      ENDIF

C     Check ISBF validity
      CALL D0PTR (CMEM, IMEM, DMEM, ISBF, ENTBF, -8, ISH, I, IER)
      IF (IER .NE. 0) GOTO 9900
      ISD = IMEM(ISH+3)
      CALL DTGET (DMEM(ISD), .TRUE., 1, NS, MRS, MS, KS, NCS, SLO, SHI,
     +            IERX)
      IF (IERX .NE. 0) THEN
         IER = -15
         GOTO 9900
      ELSEIF (NS .NE. 2) THEN
         IER = -16
         GOTO 9900
      ELSEIF (MS .NE. 3) THEN
         IER = -17
         GOTO 9900
      ENDIF

C     Compute parameters and allocate space for ICBF
      IF (NDEG .LT. 0) THEN
         KC = MAX (KP, MAX (KS(1), KS(2)))
      ELSE
         KC = NDEG + 1
      ENDIF
      NPT = (KC - 1) * (NCP - KP + 1) + 1
      NCC = NPT - 2 + KC
      LENC = 5 + KC + 4*NCC
      CALL D2BFDF (CMEM, IMEM, DMEM, LENC, ' ', ICBF, IERX)
      IF (IERX .NE. 0) THEN
         IER = -18
         GOTO 9900
      ENDIF
      ICD = IMEM(ICBF/256+3)


C     Compute and allocate workspaces ITDA and IWDA
      CALL D2DEFE (CMEM, IMEM, DMEM, ENTDA, 0, 0, 4*NPT, .FALSE., ITDA,
     +             IERX)
      IF (IERX .NE. 0) THEN
         IER = -19
         GOTO 9900
      ENDIF
      ITD = IMEM(ITDA/256+3)
      NWORK = MAX (KS(2)*(MS+2) + 3*MAX(KS(1),KS(2)) + NS + MS + 1,
     +             MAX (5*KP + MP - 1, NCC*3*KC + 2*KC**2 + 4*KC + 9))
      CALL D2DEFE (CMEM, IMEM, DMEM, ENTDA, 0, 0, NWORK, .FALSE., IWDA,
     +             IERX)
      IF (IERX .NE. 0) THEN
         IER = -19
         GOTO 9900
      ENDIF
      IWD = IMEM(IWDA/256+3)

C     Compute the NPT points in ITDA
      DO 100 I=0,NPT-1
         IF (NPT .GT. 1) THEN
            DMEM(ITD+I) = (PLO*DBLE(NPT-1-I) + PHI*DBLE(I))/DBLE (NPT-1)
         ELSE
            DMEM(ITD+I) = (PLO + PHI)/2.0D0
         ENDIF
         CALL DTSPVL (DMEM(ITD+I), DMEM(IPD), DMEM(IWD), NWORK, UV,
     +                IERX)
         IF (IERX .NE. 0) THEN
            IER = -20
            GOTO 9900
         ENDIF
         CALL DTNPVL (UV, 1, DMEM(ISD), DMEM(IWD), NWORK, XYZ, IERX)
         IF (IERX .NE. 0) THEN
            IER = -21
            GOTO 9900
         ENDIF
         DMEM(ITD+NPT+I) = XYZ(1)
         DMEM(ITD+2*NPT+I) = XYZ(2)
         DMEM(ITD+3*NPT+I) = XYZ(3)
  100 CONTINUE

C     Interpolate those points
      JCC = -2
      DO 200 I=1,3
         CALL DTNSI (NPT, DMEM(ITD), DMEM(ITD+I*NPT), KC-1, JCC,
     +               DMEM(IWD), NWORK, LENC, DMEM(ICD), LC, IERX)
         IF (IERX .NE. 0) THEN
            IER = -21 - I
            GOTO 9900
         ENDIF
  200 CONTINUE
      
C     Erase the workspace
      CALL D2ERAS (CMEM, IMEM, DMEM, IWDA, IERX)
      IF (IERX .NE. 0) THEN
         IER = -25
         GOTO 9900
      ENDIF
      CALL D2ERAS (CMEM, IMEM, DMEM, ITDA, IERX)
      IF (IERX .NE. 0) THEN
         IER = -25
         GOTO 9900
      ENDIF
      IER = 0

      RETURN

C     Error handling

 9900 CONTINUE
C     Erase workspace and result entities, suppressing error messages
      CALL DTERPT (0)
      IF (IWDA .NE. 0)  CALL D2ERAS (CMEM, IMEM, DMEM, IWDA, IERX)
      IF (ITDA .NE. 0)  CALL D2ERAS (CMEM, IMEM, DMEM, ITDA, IERX)
      IF (ICBF .NE. 0)  CALL D2ERAS (CMEM, IMEM, DMEM, ICBF, IERX)
      CALL DTERPT (1)
C     Report the error
      CALL DTERR (1, 'D0ACS ', IER, 0)
      ICBF = 0
      RETURN
      END
