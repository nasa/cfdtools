      SUBROUTINE D2LPIG (CMEM, IMEM, DMEM, ILP, IGI, JDE, ISDE, IGE,
     +                   IER)

C     Translate a Loop entity into an IGES type 141 Boundary Entity
C     and insert that IGES Entity entity into the IGES Index.  Subordinate
C     components of the Loop are translated and placed in next unused
C     directory entry lines after the JDEth.
C
C   DYNAMIC MEMORY INPUT-OUTPUT
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     ILP     MEM pointer to a Loop entity
C     IGI     MEM pointer to an IGES Index entity
C     JDE     IGES directory entry line number at which to place the new
C             IGES Boundary (IGES type 141).  Referenced entities will
C             occupy the next available directory entry spots after JDE.
C     ISDE    IGES directory entry line number of surface for which this
C             Loop is part of the boundary.  (The surface must be translated
C             to IGES form before its boundary.  Note, however, that the
C             surface is not the Bounded Surface entity, but rather the
C             B-spline surface entity of type 128.)
C   OUTPUT:
C     IGE     MEM pointer to new IGES Entity entity (Boundary)
C     IER     Returned error code.  Zero implies no errors.  If negative,
C             the IGE pointer is set null (0).
C             -1  = Dynamic memory is corrupt or uninitialized
C             -2  = ILP is a null pointer
C             -3  = Garbage value found in ILP
C             -4  = Ambiguous - ILP is garbage pointer or deleted entity
C             -5  = ILP points to a deleted entity
C             -6  = ILP does not point to a Loop entity
C             -7  = IGI is a null pointer
C             -8  = Garbage value found in IGI
C             -9  = Ambiguous - IGI is garbage pointer or deleted entity
C             -10 = IGI points to a deleted entity
C             -11 = IGI does not point to an IGES Index entity
C             -12 = JDE is not a positive odd integer
C             -13 = The JDEth directory entry is not empty.
C             -14 = Insufficient room in IGES Index.
C             -998= Unexpected error while allocating and initializing
C                   IGES Entity for Boundary (IGES type 141).
C             -999= Unexpected error while obtaining associated B-spline
C                   surface from Trimmed Surface to which Loop belongs.
C             -(1000+I) = Unexpected error from lower-level routine while
C                   translating Ith non-null Edge entity of the Loop
C
C   CALLS:
C     D0LBIG
C     D0NXDE
C     D0PTR
C     D2BFIG
C     D2IGDE
C     D2IGED
C     D2IGFI
C     D2IGSI
C     DTERR
C
C   HISTORY:
C     31Aug93  P. Kraushar   Created
C*****
C
C   Long Name Alias:
C     ENTRY D2_LOOP_TO_IGES (CMEM, IMEM, DMEM, ILP, IGI, JDE, ISDE,
C    +        IGE, IER)
C
      CHARACTER CMEM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)
      INTEGER ILP, IGI, JDE, ISDE, IGE, IER

      INTEGER ENTGI
      PARAMETER (ENTGI=237)
      INTEGER ENTLP
      PARAMETER (ENTLP=242)
      INTEGER ENTTS
      PARAMETER (ENTTS=241)
      INTEGER ENTEG
      PARAMETER (ENTEG=243)
      INTEGER ENTBF
      PARAMETER (ENTBF=246)

      INTEGER IHDRLP, IHDRGI, IHDR, ITYP, IERX, ISBF, NEDG, NEG, IEG
      INTEGER IGDE, ITS, IMBF, IMDE, IMGE, IEBF, IEDE, IEGE, I, KB, JP
      CHARACTER SUBNAM*6

      DATA SUBNAM /'D2LPIG'/

C****
      IER = 0
      IGE = 0
C     Check JDE
      IF (JDE .LE. 0 .OR. MOD (JDE, 2) .NE. 1) THEN
         IER = -12
         GOTO 9000
      ENDIF
C     Check pointer and type of ILP
      CALL D0PTR (CMEM, IMEM, DMEM, ILP, ENTLP, 0, IHDRLP, ITYP, IER)
      IF (IER .LT. 0) GOTO 9000
C     Check pointer and type of IGI
      CALL D0PTR (CMEM, IMEM, DMEM, IGI, ENTGI, -5, IHDRGI, ITYP, IER)
      IF (IER .LT. 0) GOTO 9000
C     Check that JDE is an available directory position
      IGDE = JDE - 2
      CALL D0NXDE (CMEM, IMEM, DMEM, IHDRGI, IGDE)
      IF (IGDE .NE. JDE) THEN
         IER = -13
         GOTO 9000
      ENDIF
C     Count active edges in Loop
      NEDG = 0
      KB = IMEM(IHDRLP+2) + 2
      NEG = IMEM(KB)
      DO 100 I=KB+1,KB+NEG
         IF (IMEM(I) .NE. 0)  NEDG = NEDG + 1
  100 CONTINUE
C     Allocate the IGES Entity for the Boundary Entity (Type 141)
      IER = -998
      CALL D2IGDE (CMEM, IMEM, DMEM, 4+4*NEDG, 0, 0, 0, IGI, JDE,
     +             IGE, IERX)
      IF (IERX .NE. 0) GOTO 9000
      CALL D2IGSI (CMEM, IMEM, DMEM, 'I', 141, IGE, 'D', 1, IERX)
      IF (IERX .NE. 0) GOTO 9000
      CALL D0LBIG (CMEM, IMEM, DMEM, ILP, IGE)
      CALL D2IGSI (CMEM, IMEM, DMEM, 'I', 1, IGE, 'P', 1, IERX)
      IF (IERX .NE. 0) GOTO 9000
      CALL D2IGSI (CMEM, IMEM, DMEM, 'I', 2, IGE, 'P', 2, IERX)
      IF (IERX .NE. 0) GOTO 9000
      CALL D2IGSI (CMEM, IMEM, DMEM, 'I', ISDE, IGE, 'P', 3, IERX)
      IF (IERX .NE. 0) GOTO 9000
      CALL D2IGSI (CMEM, IMEM, DMEM, 'I', NEDG, IGE, 'P', 4, IERX)
      IF (IERX .NE. 0) GOTO 9000

C     Obtain associated surface B-spline function
      IER = -999
      ITS = IMEM(KB-2)
      CALL D0PTR (CMEM, IMEM, DMEM, ITS, ENTTS, 0, IHDR, ITYP, IERX)
      IF (IERX .NE. 0) GOTO 9000
      ISBF = IMEM(IMEM(IHDR+2))
      CALL D0PTR (CMEM, IMEM, DMEM, ISBF, ENTBF, 0, IHDR, ITYP, IERX)
      IF (IERX .NE. 0) GOTO 9000
      IF (DMEM(IMEM(IHDR+3)) .NE. 2.0D0) GOTO 9000
                   
C     Translate the Edge entities into IGES B-spline Curve Entities (Type 126)
      IEDE = JDE
      JP = 4
      DO 120 I=1,NEDG
         IER = -1000 - I
C        START LOOP (FIND NEXT EDGE ENTITY)
  110       CONTINUE
            KB = KB + 1
            IEG = IMEM(KB)
            IF (IEG .EQ. 0)  GOTO 110
C        END LOOP
         CALL D0PTR (CMEM, IMEM, DMEM, IEG, ENTEG, 0, IHDR, ITYP, IERX)
         IF (IERX .NE. 0) GOTO 9000
         IEBF = IMEM(IMEM(IHDR+2))
C        Construct a required Model space curve from the composition of the
C        edge and surface functions using default degrees and tolerances
         CALL D0ACS (CMEM, IMEM, DMEM, -1, -1.0D0, IEBF, ISBF, IMBF,
     +               IERX)
         IF (IERX .NE. 0) GOTO 9000
C        Find directory lines to put the model space curve and parameter curve
         CALL D0NXDE (CMEM, IMEM, DMEM, IHDRGI, IEDE)
         IMDE = IEDE
         CALL D0NXDE (CMEM, IMEM, DMEM, IHDRGI, IEDE)
         IF (IMDE .EQ. 0 .OR. IEDE .EQ. 0) THEN
            IER = -14
            GOTO 9000
         ENDIF
C        Store the data for this edge in the IGES Boundary entity
         CALL D2IGSI (CMEM, IMEM, DMEM, 'P', IMDE, IGE, 'P', JP+1, IERX)
         IF (IERX .NE. 0) GOTO 9000
         CALL D2IGSI (CMEM, IMEM, DMEM, 'I', 1, IGE, 'P', JP+2, IERX)
         IF (IERX .NE. 0) GOTO 9000
         CALL D2IGSI (CMEM, IMEM, DMEM, 'I', 1, IGE, 'P', JP+3, IERX)
         IF (IERX .NE. 0) GOTO 9000
         CALL D2IGSI (CMEM, IMEM, DMEM, 'P', IEDE, IGE, 'P', JP+4, IERX)
         IF (IERX .NE. 0) GOTO 9000
         JP = JP + 4
C        Translate the model space curve to IGES Entity type 126
         CALL D2BFIG (CMEM, IMEM, DMEM, IMBF, IGI, IMDE, IMGE, IERX)
         IF (IERX .NE. 0) GOTO 9000
C        Set IGES model space curve entity to physically dependent status
         CALL D2IGSI (CMEM, IMEM, DMEM, 'I', 1, IMGE, 'D', 22, IERX)
         IF (IERX .NE. 0) GOTO 9000
C        Translate the parameter space curve to IGES Entity type 126
         CALL D2BFIG (CMEM, IMEM, DMEM, IEBF, IGI, IEDE, IEGE, IERX)
         IF (IERX .NE. 0) GOTO 9000
C        Set IGES parameter space curve entity to physically dependent status
         CALL D2IGSI (CMEM, IMEM, DMEM, 'I', 1, IEGE, 'D', 22, IERX)
         IF (IERX .NE. 0) GOTO 9000
C        Set IGES parameter space curve entity to 2D-parametric status
         CALL D2IGSI (CMEM, IMEM, DMEM, 'I', 5, IEGE, 'D', 23, IERX)
         IF (IERX .NE. 0) GOTO 9000
C        Erase the model space curve B-spline Function entity
         CALL D2ERAS (CMEM, IMEM, DMEM, IMBF, IERX)
         IF (IERX .NE. 0) GOTO 9000
  120 CONTINUE
      IER = 0
      RETURN

C     Error Handling

 9000 CONTINUE
C     Report error
      IF (IER .LE. -998) THEN
         CALL DTERR (4, SUBNAM, IER, 0)
      ELSE
         CALL DTERR (1, SUBNAM, IER, 0)
      ENDIF
C     Attempt to erase those entities that may have been generated
      IF (IGE .NE. 0) THEN
         CALL DTERPT (0)
         CALL D2IGED (CMEM, IMEM, DMEM, IGE, IERX)
         CALL DTERPT (1)
         IGE = 0
      ENDIF
C     End of error handling
      RETURN
      END
