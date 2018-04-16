      SUBROUTINE D2IGER (CMEM, IMEM, DMEM, IGIE, IER)

C     Erase an IGES Entity or IGES Index entity.  Do not attempt to locate
C     and erase IGES "physically dependent" entities.  Do not erase an Index
C     that is not empty of IGES Entities.  See D2IGED to erase an entity
C     together with certain "physically dependent" entities.
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
C             -3  = IGIE is not a MEM pointer.
C             -6  = IGIE does not point to an IGES Entity or IGES Index 
C                   entity.
C             -7  = Attempted to erase IGES Index that is not empty.
C
C   CALLS:
C     D0GEER
C     D0PTR
C     D2ERAS
C     DTERPT
C     DTERR
C
C   HISTORY:
C     08Sep93 P. Kraushar   Created
C*****
C
C   Long Name Alias:
C     ENTRY D2_IGES_ERASE (CMEM, IMEM, DMEM, IGIE, IER)
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

      DATA SUBNAM /'D2IGER'/

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
         CALL D0GEER (CMEM, IMEM, DMEM, IHDR, 0, IER)
      ELSE IF (ITYP .EQ. ENTGI) THEN
C        Extract critical data from IGES Index entity
         IXGI = IMEM(IHDR+2)
         IGE = IMEM(IXGI)
         NDEX = IMEM(IXGI+1)/2
C        Verify that all IGES Entities have been deleted
         DO 300 I=IXGI+2,IXGI+1+NDEX
            IF (IMEM(I) .NE. 0) THEN
               IER = -7
               GOTO 9000
            ENDIF
  300    CONTINUE
C        Verify Global IGES Entity and delete it, ignoring any errors.
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
