      SUBROUTINE D2IGED (CMEM, IMEM, DMEM, IGIE, IER)

C     Erase an IGES Entity or IGES Index entity.  Also locate and erase
C     certain IGES "physically dependent" entities depending on IGES type.  
C     If erasing an IGES Index, first erase all IGES Entities listed in 
C     that Index.  See D2IGER to restrict erasure to just the entity given.
C
C   DYNAMIC MEMORY INPUT-OUTPUT
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     IGIE    MEM pointer to IGES Entity or IGES Index entity
C   OUTPUT:
C     IER     Returned error code.  Zero implies no errors.
C             -1  = Dynamic memory is corrupt or uninitialized.
C             -3  = IGIE is a garbage value, not a pointer.
C             -6  = IGIE does not point to an IGES Entity or IGES Index 
C                   entity.
C
C   CALLS:
C     D0PTR
C     D1IGED
C     D2ERAS
C     D2IGER
C     DTERPT
C     DTERR
C
C   HISTORY:
C     08Sep93 P. Kraushar   Created
C*****
C
C   Long Name Alias:
C     ENTRY D2_IGES_ERASE_DEP (CMEM, IMEM, DMEM, IGIE, IER)
C
      CHARACTER CMEM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)
      INTEGER IGIE, IER

      INTEGER ENTGE
      PARAMETER (ENTGE=238)
      INTEGER ENTGI
      PARAMETER (ENTGI=237)

      INTEGER IHDR, ITYP, IERX, IGE, ICQ, IXI, I, NDEX, IXGI
      CHARACTER SUBNAM*6

      DATA SUBNAM /'D2IGED'/

C****
      IER = 0
      IF (IGIE .EQ. 0)  RETURN
C     Check pointer and type of IGIE
      CALL D0PTR (CMEM, IMEM, DMEM, IGIE, 0, 0, IHDR, ITYP, IERX)
      IF (IERX .LT. 0) THEN
         IF (IERX .EQ. -5 .OR. IERX .EQ. -4)  RETURN
         IER = IERX
         GOTO 9000
      ENDIF
      IF (ITYP .EQ. ENTGE) THEN
C        Identify IGES type
         CALL D1FEEI (CMEM, IMEM, DMEM, IHDR, 5, ITYP, IERX)
         IF (ITYP .EQ. 141) THEN
            CALL D0G141 (CMEM, IMEM, DMEM, IHDR, 0, IERX)
         ELSE IF (ITYP .EQ. 143) THEN
            CALL D0G143 (CMEM, IMEM, DMEM, IHDR, 0, IERX)
         ELSE
            CALL D0GEER (CMEM, IMEM, DMEM, IHDR, 0, IERX)
         ENDIF
      ELSE IF (ITYP .EQ. ENTGI) THEN
C        Extract critical data from IGES Index entity
         IXGI = IMEM(IHDR+2) + 1
         NDEX = IMEM(IXGI)/2
C        Delete any remaining IGES Entities
         DO 300 I=IXGI+1,IXGI+NDEX
            IF (IMEM(I) .NE. 0) THEN 
               CALL D0GEER (CMEM, IMEM, DMEM, IMEM(I), IGIE, IERX)
               IMEM(I) = 0
            ENDIF
  300    CONTINUE
C        Verify Global IGES Entity and delete it, ignoring any errors.
         IGE = IMEM(IXGI-1)
         CALL DTERPT (0)
         CALL D0PTR (CMEM, IMEM, DMEM, IGE, ENTGE, 0, IHDR, ITYP, IERX)
         IF (IERX .EQ. 0) THEN
            IXI = IMEM(IHDR+2)
            ICQ = IMEM(IXI+1)
C           Delete string values of global parameter
            IF (ICQ .NE. 0) CALL D2ERAS (CMEM, IMEM, DMEM, ICQ, IERX)
C           Delete Global IGES Entity entity
            CALL D2ERAS (CMEM, IMEM, DMEM, IGE, IERX)
         ENDIF
C        Delete IGES Index itself
         CALL D2ERAS (CMEM, IMEM, DMEM, IGIE, IERX)
C        Resume normal error processing
         CALL DTERPT (1)
      ELSE
C        Entity is neither IGES Entity nor IGES Index
         IER = -6
         GOTO 9000
      ENDIF
      RETURN

C     Error Handling

 9000 CONTINUE
C     Report error
      CALL DTERR (1, SUBNAM, IER, 0)
      RETURN
      END
