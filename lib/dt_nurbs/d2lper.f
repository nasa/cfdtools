      SUBROUTINE D2LPER (CMEM, IMEM, DMEM, ILP, IER)

C     Erase a Loop entity and all component Edge entities, disconnecting 
C     from any Trimmed Surface of which the Loop might be a boundary.
C
C   DYNAMIC MEMORY INPUT-OUTPUT
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     ILP     MEM pointer to Loop entity
C   OUTPUT:
C     IER     Returned error code.  Zero implies no errors.
C             -1  = Dynamic memory is corrupt or uninitialized.
C             -3  = ILP is a garbage value, not a pointer.
C             -6  = ILP does not point to a Loop entity.
C
C   CALLS:
C     D0PTR
C     D2EGER
C     D2ERAS
C     DTERPT
C     DTERR
C
C   HISTORY:
C     08Sep93 P. Kraushar   Created
C*****
C
C   Long Name Alias:
C     ENTRY D2_LOOP_ERASE (CMEM, IMEM, DMEM, ILP, IER)
C
      CHARACTER CMEM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)
      INTEGER ILP, IER
                    
      INTEGER ENTLP
      PARAMETER (ENTLP=242)
      INTEGER ENTTS
      PARAMETER (ENTTS=241)

      INTEGER IHDR, ITYP, IERX, IXI, IEG, IXLP, ITS, NEG, IXL, I
      CHARACTER SUBNAM*6

      DATA SUBNAM /'D2LPER'/

C****
      IER = 0
      IF (ILP .EQ. 0)  RETURN
C     Check pointer and type of ILP
      CALL D0PTR (CMEM, IMEM, DMEM, ILP, ENTLP, 0, IHDR, ITYP, IERX)
      IF (IERX .LT. 0) THEN
         IF (IERX .EQ. -5 .OR. IERX .EQ. -4)  RETURN
         IER = IERX
         GOTO 9000
      ENDIF
      IXL = IMEM(IHDR+2)
      ITS  = IMEM(IXL)
      IXLP = IMEM(IXL+1)
      NEG  = IMEM(IXL+2)
C     If part of a trimmed surface, disconnect it.  Ignore bad data.
      IF (ITS .NE. 0) THEN
         CALL D0PTR (CMEM, IMEM, DMEM, ITS, ENTTS, 0, IHDR, ITYP, IERX)
         IF (IERX .EQ. 0) THEN
            IXI = IMEM(IHDR+2) + 4
            IF (IXLP .GT. 0 .AND. IXLP .LE. IMEM(IXI)) THEN
C              Don't reference IMEM(IXI+IXLP) if IXLP is not in range
               IF (IMEM(IXI+IXLP) .EQ. ILP) THEN
                  IMEM(IXI+IXLP) = 0
                  ITS = 0
               ENDIF
            ENDIF
            IF (ITS .NE. 0) THEN
C              Check all Loop locations in the Trimmed Surface entity
               DO 100 I=IXI+1,IXI+IMEM(IXI)
                  IF (IMEM(I) .EQ. ILP)  IMEM(I) = 0
  100          CONTINUE
            ENDIF
         ENDIF
      ENDIF
C     Now erase any and all component Edge entities, ignoring errors.
      CALL DTERPT (0)
      DO 200 I=IXL+3,IXL+2+NEG
         IEG = IMEM(I)
         IF (IEG .NE. 0)  CALL D2EGER (CMEM, IMEM, DMEM, IEG, IERX)
  200 CONTINUE
C     Finally, erase the Loop entity, ignoring errors.
      CALL D2ERAS (CMEM, IMEM, DMEM, ILP, IERX)
      CALL DTERPT (1)
      RETURN                    

C     Error Handling

 9000 CONTINUE
C     Report error
      CALL DTERR (1, SUBNAM, IER, 0)
      RETURN
      END
