      SUBROUTINE D2IGTS (CMEM, IMEM, DMEM, IGE, IJS, IXTS, ITS, IER)

C     Translate an IGES type 143 Bounded Surface into a Trimmed Surface entity.
C     Optionally, connect the Trimmed Surface to a Joined Surface.  Note that
C     any matching edges must be connected separately using D2CNEE.
C
C   DYNAMIC MEMORY INPUT-OUTPUT
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     IGE     MEM pointer to an IGES Entity entity
C     IJS     MEM pointer to Joined Surface to which ITS belongs, if any.
C             Otherwise, IJS must be zero.
C     IXTS    Index of ITS in the list of trimmed surface entities in IJS.
C             Ignored if IJS is zero.
C   OUTPUT:
C     ITS     MEM pointer to new Trimmed Surface entity
C     IER     Returned error code.  Zero implies no errors.  If negative,
C             the ITS pointer is set null (0).
C             +1  = Translation successful, but requested connection to
C                   joined surface was unsuccessful
C             -1  = Dynamic memory is corrupt or uninitialized
C             -2  = IGE is a null pointer
C             -3  = Garbage value found in IGE
C             -4  = Ambiguous - garbage pointer or deleted entity
C             -5  = IGE points to a deleted entity
C             -6  = IGE does not point to an IGES Entity entity
C             -7  = IGE is not IGES Entity type 143 (Bounded Surface)
C             -8  = IGE is not IGES Entity type 143 type 1 (i.e. its boundary
C                   is not given in the surface's parameter space)
C             -9  = Associated surface not given as IGES Entity type 128
C             -10 = Unsuccessful at allocating memory for Trimmed Surface
C             -999= Unexpected error from lower-level routine before
C                   starting boundary translations
C             -(1000+I) = Unexpected error from lower-level routine while
C                   translating Ith boundary component
C
C   CALLS:
C     DOLABL
C     D0PTR
C     D2CNTJ
C     D2IGBF
C     D2IGFI
C     D2IGLP
C     D2TSDF
C     D2TSER
C     DTERR
C
C   HISTORY:
C     15Aug93 P. Kraushar   Created
C*****
C
C   Long Name Alias:
C     ENTRY D2_IGES_TO_TRIMSURF (CMEM, IMEM, DMEM, IGE, IJS, IXTS,
C    +                           ITS, IER)
C
      CHARACTER CMEM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)
      INTEGER IGE, IJS, IXTS, ITS, IER

      INTEGER ENTGE
      PARAMETER (ENTGE=238)

      CHARACTER LABEL*16
      INTEGER IHDR, ITYP, IERX, ISDE, ISGE, ISBF, NBNDRY, IBDE, IBGE
      INTEGER ILP, I
      CHARACTER SUBNAM*6

      DATA SUBNAM /'D2IGTS'/

C****
      IER = 0
      ITS = 0
      ISBF = 0
C     Check pointer and type of IGE
      CALL D0PTR (CMEM, IMEM, DMEM, IGE, ENTGE, 0, IHDR, ITYP, IER)
      IF (IER .LT. 0) GOTO 9000
      IER = -999
      CALL D2IGFI (CMEM, IMEM, DMEM, IGE, 'D', 1, ITYP, IERX)
      IF (IERX .NE. 0) GOTO 9000
      IF (ITYP .NE. 143) THEN
         IER = -7
         GOTO 9000
      ENDIF
      CALL D2IGFI (CMEM, IMEM, DMEM, IGE, 'P', 1, ITYP, IERX)
      IF (IERX .NE. 0) GOTO 9000
      IF (ITYP .NE. 1) THEN
         IER = -8
         GOTO 9000
      ENDIF

C     Obtain, check type, and translate associated surface entity
      CALL D2IGFI (CMEM, IMEM, DMEM, IGE, 'P', 2, ISDE, IERX)
      IF (IERX .NE. 0) GOTO 9000
      CALL D2IGFI (CMEM, IMEM, DMEM, IGE, 'I', ISDE, ISGE, IERX)
      IF (IERX .NE. 0) GOTO 9000
      CALL D2IGFI (CMEM, IMEM, DMEM, ISGE, 'D', 1, ITYP, IERX)
      IF (IERX .NE. 0) GOTO 9000
      IF (ITYP .NE. 128) THEN
         IER = -9
         GOTO 9000
      ENDIF
      CALL D2IGBF (CMEM, IMEM, DMEM, ISGE, ISBF, IERX)
      IF (IERX .NE. 0) GOTO 9000

C     Get number of loop components and allocate the Trimmed Surface entity
      CALL D2IGFI (CMEM, IMEM, DMEM, IGE, 'P', 3, NBNDRY, IERX)
      IF (IERX .NE. 0) GOTO 9000
      CALL D0LABL (CMEM, IMEM, DMEM, IHDR, LABEL, I)
      CALL D2TSDF (CMEM, IMEM, DMEM, NBNDRY, ISBF, 0, 0, 1, LABEL,
     +             ITS, IERX)
      IF (IERX .NE. 0) THEN
         IER = -10
         GOTO 9000
      ENDIF
                   
C     Translate the boundary into Loop entities
      DO 100 I=1,NBNDRY
         IER = -1000 - I
         CALL D2IGFI (CMEM, IMEM, DMEM, IGE, 'P', I+3, IBDE, IERX)
         IF (IERX .NE. 0) GOTO 9000
         CALL D2IGFI (CMEM, IMEM, DMEM, IGE, 'I', IBDE, IBGE, IERX)
         IF (IERX .NE. 0) GOTO 9000
         CALL D2IGLP (CMEM, IMEM, DMEM, IBGE, ITS, I, ILP, IERX)
         IF (IERX .NE. 0) GOTO 9000
  100 CONTINUE
      IER = 0

C     Finally, connect to Joined Surface, if selected.
      IF (IJS .NE. 0) THEN
         CALL D2CNTJ (CMEM, IMEM, DMEM, IXTS, ITS, IJS, IERX)
         IF (IERX .NE. 0) THEN
            IER = 1
            CALL DTERR (0, SUBNAM, IER, 0)
         ENDIF
      ENDIF
      RETURN

C     Error Handling

 9000 CONTINUE
C     Report error
      IF (IER .LE. -999) THEN
         CALL DTERR (4, SUBNAM, IER, 0)
      ELSE
         CALL DTERR (1, SUBNAM, IER, 0)
      ENDIF
C     Attempt to erase any partially generated entity
      IF (ITS .NE. 0) THEN
         CALL DTERPT (0)
         CALL D2TSER (CMEM, IMEM, DMEM, ITS, IERX)
         CALL DTERPT (1)
         ITS = 0
      ENDIF
C     End of error handling
      RETURN
      END
