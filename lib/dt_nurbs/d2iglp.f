      SUBROUTINE D2IGLP (CMEM, IMEM, DMEM, IGE, ITS, IXLP, ILP, IER)

C     Translate an IGES type 141 Boundary into a Loop entity.
C     Optionally, connect the Loop to a Trimmed Surface.
C
C   DYNAMIC MEMORY INPUT-OUTPUT
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     IGE     MEM pointer to an IGES Entity entity
C     ITS     MEM pointer to Trimmed Surface to which the Loop belongs, if 
C             any.  Otherwise, ITS must be zero.
C     IXLP    Index of Loop in the list of loop entities in ITS.
C             Ignored if ITS is zero.
C   OUTPUT:
C     ILP     MEM pointer to new Loop of a Trimmed Surface entity.
C     IER     Returned error code.  Zero implies no errors.  If negative,
C             the ILP pointer is set null (0).
C             +3  = Both problems +1 and +2 occurred.
C             +2  = The translation was completed, but the Edges are not
C                   connected consecutively and the Loop is technically
C                   invalid.  A perfectly valid IGES Boundary can yield a
C                   Loop not satisfying all the requirements.  See Manual.
C             +1  = Translation successful, but requested connection to
C                   trimmed surface was unsuccessful.
C             -1  = Dynamic memory is corrupt or uninitialized
C             -2  = IGE is a null pointer
C             -3  = Garbage value found in IGE
C             -4  = Ambiguous - garbage pointer or deleted entity
C             -5  = IGE points to a deleted entity
C             -6  = IGE does not point to an IGES Entity entity
C             -7  = IGE is not IGES Entity type 141 (Boundary)
C             -8  = IGE is not IGES Entity type 141 type 1 (i.e. the boundary
C                   is not given in the surface's parameter space)
C             -9  = Associated surface not given as IGES Entity type 128
C             -10 = Unsuccessful at allocating memory for Loop
C             -999= Unexpected error from lower-level routine before
C                   starting individual curve translations
C             -(1000+I) = Unexpected error from lower-level routine while
C                   translating Ith individual curve component
C
C   CALLS:
C     D0LABL
C     D0PTR
C     D2CNLT
C     D2LPDF
C     D2LPER
C     D2IGEG
C     D2IGFD
C     D2IGFI
C     D2IPFD
C     D2IPFI
C     DTERR
C
C   HISTORY:
C     26Aug93 P. Kraushar   Created
C*****
C
C   Long Name Alias:
C     ENTRY D2_IGES_TO_LOOP (CMEM, IMEM, DMEM, IGE, ITS, IXLP,
C    +                       ILP, IER)
C
      CHARACTER CMEM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)
      INTEGER IGE, ITS, IXLP, ILP, IER

      INTEGER ENTGE
      PARAMETER (ENTGE=238)

      CHARACTER LABEL*16
      INTEGER IHDR, ITYP, IERX, ISDE, ISGE, KNT(4), I, NCURV, NPCURV, JP
      INTEGER NTCURV, LENLBL, JER, J, IEDE, IEGE, IEG, K, M, IXEG
      DOUBLE PRECISION PBD(4), TOLU, TOLV, PIU, PIV, PBU, PBV, PEU, PEV
      CHARACTER SUBNAM*6

      DATA SUBNAM /'D2IGLP'/

C****
      JER = 0
      ILP = 0
C     Check pointer and type of IGE
      CALL D0PTR (CMEM, IMEM, DMEM, IGE, ENTGE, 0, IHDR, ITYP, IER)
      IF (IER .LT. 0) GOTO 9000
      IER = -999
      CALL D2IGFI (CMEM, IMEM, DMEM, IGE, 'D', 1, ITYP, IERX)
      IF (IERX .NE. 0) GOTO 9000
      IF (ITYP .NE. 141) THEN
         IER = -7
         GOTO 9000
      ENDIF
      CALL D2IGFI (CMEM, IMEM, DMEM, IGE, 'P', 1, ITYP, IERX)
      IF (IERX .NE. 0) GOTO 9000
      IF (ITYP .NE. 1) THEN
         IER = -8
         GOTO 9000
      ENDIF

C     Obtain, check type, and extract data from associated surface entity
      CALL D2IGFI (CMEM, IMEM, DMEM, IGE, 'P', 3, ISDE, IERX)
      IF (IERX .NE. 0) GOTO 9000
      CALL D2IGFI (CMEM, IMEM, DMEM, IGE, 'I', ISDE, ISGE, IERX)
      IF (IERX .NE. 0) GOTO 9000
      CALL D2IGFI (CMEM, IMEM, DMEM, ISGE, 'D', 1, ITYP, IERX)
      IF (IERX .NE. 0) GOTO 9000
      IF (ITYP .NE. 128) THEN
         IER = -9
         GOTO 9000
      ENDIF
C     Get 4 integer parameters: u index max, v index max, u degree, v degree
      CALL D2IPFI (CMEM, IMEM, DMEM, 4, ISGE, 1, 1, 1, KNT, IERX)
      IF (IERX .NE. 0) GOTO 9000
C     Change to DTNURBS style parameters: u ncoef, v ncoef, u order, v order
      DO 100 I=1,4
         KNT(I) = KNT(I) + 1
  100 CONTINUE
C     Get parameter range limits
      I = 10 + (KNT(1) + KNT(3)) + (KNT(2) + KNT(4)) + 4*KNT(1)*KNT(2)
      CALL D2IPFD (CMEM, IMEM, DMEM, 4, ISGE, I, 1, 1, PBD, IERX)
      IF (IERX .NE. 0) GOTO 9000
C     Compute absolute tolerances for u and v as relative to parameter ranges
      TOLU = (PBD(2) - PBD(1))*1.0E-6
      TOLV = (PBD(4) - PBD(3))*1.0E-6

C     Get number of edge components and allocate the Loop entity
      CALL D2IGFI (CMEM, IMEM, DMEM, IGE, 'P', 4, NCURV, IERX)
      IF (IERX .NE. 0) GOTO 9000
      JP = 7
      NTCURV = 0
      DO 200 I=1,NCURV
         CALL D2IGFI (CMEM, IMEM, DMEM, IGE, 'P', JP, NPCURV, IERX)
         IF (IERX .NE. 0) GOTO 9000
         NTCURV = NTCURV + NPCURV
         JP = JP + NPCURV + 3
  200 CONTINUE
      CALL D0LABL (CMEM, IMEM, DMEM, IHDR, LABEL, LENLBL)
      CALL D2LPDF (CMEM, IMEM, DMEM, NTCURV, 0, 0, LABEL, ILP, IERX)
      IF (IERX .NE. 0) THEN
         IER = -10
         GOTO 9000
      ENDIF
                   
C     Translate the boundary curve components into Edge entities
      JP = 7
      IXEG = 0
      DO 310 I=1,NCURV
         CALL D2IGFI (CMEM, IMEM, DMEM, IGE, 'P', JP, NPCURV, IERX)
         IF (IERX .NE. 0) GOTO 9000
         DO 300 J=1,NPCURV
            IXEG = IXEG + 1
            IEG = 0
            IER = -1000 - IXEG
            CALL D2IGFI (CMEM, IMEM, DMEM, IGE, 'P', JP+J, IEDE, IERX)
            IF (IERX .NE. 0) GOTO 9000
            CALL D2IGFI (CMEM, IMEM, DMEM, IGE, 'I', IEDE, IEGE, IERX)
            IF (IERX .NE. 0) GOTO 9000
            CALL D2IGEG (CMEM, IMEM, DMEM, IEGE, ILP, IXEG, IEG, IERX)
            IF (IERX .NE. 0) GOTO 9000
            IF (JER .EQ. 0) THEN
C              Obtain first control point as Edge beginning vertex
               CALL D2IGFI (CMEM, IMEM, DMEM, IEGE, 'P', 1, K, IERX)
               IF (IERX .NE. 0) GOTO 9000
               CALL D2IGFI (CMEM, IMEM, DMEM, IEGE, 'P', 2, M, IERX)
               IF (IERX .NE. 0) GOTO 9000
               M = 10 + (K + M) + K
               CALL D2IGFD (CMEM, IMEM, DMEM, IEGE, 'P', M, PBU, IERX)
               IF (IERX .NE. 0) GOTO 9000
               CALL D2IGFD (CMEM, IMEM, DMEM, IEGE, 'P', M+1, PBV, IERX)
               IF (IERX .NE. 0) GOTO 9000
C              Compare with previous end vertex
               IF (IXEG .EQ. 1) THEN
C                 Save initial vertex of first Edge for later
                  PIU = PBU
                  PIV = PBV
               ELSE IF (ABS(PBU-PEU) .GT. TOLU 
     +          .OR. ABS(PBV-PEV) .GT. TOLV) THEN
                  JER = 2
               ENDIF
               M = M + 3*K
C              Get last control point as Edge end vertex
               CALL D2IGFD (CMEM, IMEM, DMEM, IEGE, 'P', M, PEU, IERX)
               IF (IERX .NE. 0) GOTO 9000
               CALL D2IGFD (CMEM, IMEM, DMEM, IEGE, 'P', M+1, PEV, IERX)
               IF (IERX .NE. 0) GOTO 9000
               IF (IXEG .EQ. NTCURV) THEN
C                 Check if end of last Edge matches beginning of first Edge
                  IF (ABS(PEU-PIU) .GT. TOLU
     +             .OR. ABS(PEV-PIV) .GT. TOLV)  JER = 2
               ENDIF
            ENDIF
  300    CONTINUE
         JP = JP + NPCURV + 3
  310 CONTINUE

C     Finally, connect to Trimmed Surface, if selected.
      IF (ITS .NE. 0) THEN
         CALL D2CNLT (CMEM, IMEM, DMEM, IXLP, ILP, ITS, IERX)
         IF (IERX .NE. 0) THEN
            JER = JER + 1
         ENDIF
      ENDIF
      IER = JER
      IF (IER .GT. 0)  CALL DTERR (0, SUBNAM, IER, 0)
      RETURN

C     Error Handling

 9000 CONTINUE
C     Report error
      IF (IER .LE. -999) THEN
         CALL DTERR (4, SUBNAM, IER, 0)
      ELSE
         CALL DTERR (1, SUBNAM, IER, 0)
      ENDIF
C     Attempt to erase those entities that may have been generated
      IF (ILP .NE. 0) THEN 
         CALL DTERPT (0)
         CALL D2LPER (CMEM, IMEM, DMEM, ILP, IERX)
         CALL DTERPT (1)
         ILP = 0
      ENDIF
C     End of error handling
      RETURN
      END
