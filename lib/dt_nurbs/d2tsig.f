      SUBROUTINE D2TSIG (CMEM, IMEM, DMEM, ITS, IGI, JDE, IGE, IER)

C     Translate a Trimmed Surface entity into an IGES type 143 Bounded Surface
C     and insert that IGES Entity entity into the IGES Index.  Subordinate
C     components of the Trimmed Surface are translated and placed in next
C     unused directory entry lines after the JDEth.
C
C   DYNAMIC MEMORY INPUT-OUTPUT
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     ITS     MEM pointer to a Trimmed Surface entity
C     IGI     MEM pointer to an IGES Index entity
C     JDE     IGES directory entry line number at which to place the new
C             IGES Bounded Surface (IGES type 143).  Referenced entities
C             will occupy the next available directory entry spots after
C             JDE.
C   OUTPUT:
C     IGE     MEM pointer to new IGES Entity entity (Bounded Surface)
C     IER     Returned error code.  Zero implies no errors.  If negative,
C             the IGE pointer is set null (0).
C             -1  = Dynamic memory is corrupt or uninitialized
C             -2  = ITS is a null pointer
C             -3  = Garbage value found in ITS
C             -4  = Ambiguous - ITS is garbage pointer or deleted entity
C             -5  = ITS points to a deleted entity
C             -6  = ITS does not point to a Trimmed Surface entity
C             -7  = IGI is a null pointer
C             -8  = Garbage value found in IGI
C             -9  = Ambiguous - IGI is garbage pointer or deleted entity
C             -10 = IGI points to a deleted entity
C             -11 = IGI does not point to an IGES Index entity
C             -12 = JDE is not a positive odd integer
C             -13 = The JDEth directory entry is not empty.
C             -14 = Unable to find or make room in IGES Index.
C             -998= Unexpected error while allocating and initializing
C                   IGES Entity for Bounded Surface (IGES type 143)
C             -999= Unexpected error from lower-level routine while
C                   translating B-spline Surface into IGES type 128 Entity
C             -(1000+I) = Unexpected error from lower-level routine while
C                   translating Ith boundary Loop entity
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
C     D2LPIG
C     DTERR
C
C   HISTORY:
C     16Aug93 P. Kraushar   Created
C*****
C
C   Long Name Alias:
C     ENTRY D2_TRIMSURF_TO_IGES (CMEM, IMEM, DMEM, ITS, IGI, JDE, 
C    +        IGE, IER)
C
      CHARACTER CMEM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)
      INTEGER ITS, IGI, JDE, IGE, IER

      INTEGER ENTGI
      PARAMETER (ENTGI=237)
      INTEGER ENTTS
      PARAMETER (ENTTS=241)

      INTEGER IHDRTS, IHDRGI, ITYP, IERX, ISDE, ISGE, ISBF, NBNDRY
      INTEGER ILP, ILDE, ILGE, IGDE
      INTEGER I, KB, NLP
      CHARACTER SUBNAM*6

      DATA SUBNAM /'D2TSIG'/

C****
      IER = 0
      IGE = 0
      ISGE = 0
C     Check JDE
      IF (JDE .LE. 0 .OR. MOD (JDE, 2) .NE. 1) THEN
         IER = -12
         GOTO 9000
      ENDIF
C     Check pointer and type of ITS
      CALL D0PTR (CMEM, IMEM, DMEM, ITS, ENTTS, 0, IHDRTS, ITYP, IER)
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
C     Count active loops in trimmed surface
      NBNDRY = 0
      KB = IMEM(IHDRTS+2) + 4
      NLP = IMEM(KB)
      DO 100 I=KB+1,KB+NLP
         IF (IMEM(I) .NE. 0)  NBNDRY = NBNDRY + 1
  100 CONTINUE
C     Allocate the IGES Entity for the Bounded Surface Entity (Type 143)
      IER = -998
      CALL D2IGDE (CMEM, IMEM, DMEM, NBNDRY+3, 0, 0, 0, IGI, JDE,
     +             IGE, IERX)
      IF (IERX .NE. 0) GOTO 9000
      CALL D2IGSI (CMEM, IMEM, DMEM, 'I', 143, IGE, 'D', 1, IERX)
      IF (IERX .NE. 0) GOTO 9000
      CALL D0LBIG (CMEM, IMEM, DMEM, ITS, IGE)
      CALL D2IGSI (CMEM, IMEM, DMEM, 'I', 1, IGE, 'P', 1, IERX)
      IF (IERX .NE. 0) GOTO 9000
      CALL D2IGSI (CMEM, IMEM, DMEM, 'I', NBNDRY, IGE, 'P', 3, IERX)
      IF (IERX .NE. 0) GOTO 9000

C     Obtain and translate associated surface entity
      IER = -999
      ISBF = IMEM(KB-4)
      ISDE = IGDE
      CALL D0NXDE (CMEM, IMEM, DMEM, IHDRGI, ISDE)
      IF (ISDE .EQ.0) THEN
         IER = -14
         GOTO 9000
      ENDIF
      CALL D2BFIG (CMEM, IMEM, DMEM, ISBF, IGI, ISDE, ISGE, IERX)
      IF (IERX .NE. 0) GOTO 9000
      CALL D2IGSI (CMEM, IMEM, DMEM, 'P', ISDE, IGE, 'P', 2, IERX)
      IF (IERX .NE. 0) GOTO 9000
C     Mark IGES surface entity as physically dependent
      CALL D2IGSI (CMEM, IMEM, DMEM, 'I', 1, ISGE, 'D', 22, IERX)
      IF (IERX .NE. 0) GOTO 9000
                   
C     Translate the boundary Loop entities into IGES Entities (Type 141)
      ILDE = ISDE
      DO 120 I=1,NBNDRY
         IER = -1000 - I
C        START LOOP (FIND NEXT LOOP ENTITY)
  110       CONTINUE
            KB = KB + 1
            ILP = IMEM(KB)
            IF (ILP .EQ. 0)  GOTO 110
C        END LOOP
         CALL D0NXDE (CMEM, IMEM, DMEM, IHDRGI, ILDE)
         IF (ILDE .EQ. 0) THEN
            IER = -13
            GOTO 9000
         ENDIF
         CALL D2LPIG (CMEM, IMEM, DMEM, ILP, IGI, ILDE, ISDE, ILGE,
     +                IERX)
         IF (IERX .NE. 0)  GOTO 9000
         CALL D2IGSI (CMEM, IMEM, DMEM, 'P', ILDE, IGE, 'P', 3+I, IERX)
         IF (IERX .NE. 0)  GOTO 9000
C        Mark IGES Boundary entity as physically dependent
         CALL D2IGSI (CMEM, IMEM, DMEM, 'I', 1, ILGE, 'D', 22, IERX)
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
