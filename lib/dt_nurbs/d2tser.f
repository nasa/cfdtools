      SUBROUTINE D2TSER (CMEM, IMEM, DMEM, ITS, IER)

C     Erase a Trimmed Surface entity, its surface B-spline entity and all 
C     component Loop and Edge entities, and disconnect from the associated 
C     Joined Surface, if any.
C
C   DYNAMIC MEMORY INPUT-OUTPUT
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     ITS     MEM pointer to Trimmed Surface entity
C   OUTPUT:
C     IER     Returned error code.  Zero implies no errors.
C             -1  = Dynamic memory is corrupt or uninitialized.
C             -3  = ITS is a garbage value, not a pointer.
C             -6  = ITS does not point to a Trimmed Surface entity.
C
C   CALLS:
C     D0PTR
C     D2ERAS
C     D2LPER
C     DTERPT
C     DTERR
C
C   HISTORY:
C     08Sep93 P. Kraushar   Created
C*****
C
C   Long Name Alias:
C     ENTRY D2_TRIMSURF_ERASE (CMEM, IMEM, DMEM, ITS, IER)
C
      CHARACTER CMEM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)
      INTEGER ITS, IER
                    
      INTEGER ENTTS
      PARAMETER (ENTTS=241)
      INTEGER ENTJS
      PARAMETER (ENTJS=240)

      INTEGER IHDR, ITYP, IERX, IXI, ILP, IXTS, IJS, NLP, IXT, I, ISBF
      CHARACTER SUBNAM*6

      DATA SUBNAM /'D2TSER'/

C****
      IER = 0
      IF (ITS .EQ. 0)  RETURN
C     Check pointer and type of ITS
      CALL D0PTR (CMEM, IMEM, DMEM, ITS, ENTTS, 0, IHDR, ITYP, IERX)
      IF (IERX .LT. 0) THEN
         IF (IERX .EQ. -5 .OR. IERX .EQ. -4)  RETURN
         IER = IERX
         GOTO 9000
      ENDIF
      IXT = IMEM(IHDR+2)
      ISBF = IMEM(IXT)
      IJS  = IMEM(IXT+1)
      IXTS = IMEM(IXT+2)
C     Don't care about orientation in IXT+3
      NLP  = IMEM(IXT+4)
C     If part of a joined surface, disconnect it.  Ignore bad data.
      IF (IJS .NE. 0) THEN
         CALL D0PTR (CMEM, IMEM, DMEM, IJS, ENTJS, 0, IHDR, ITYP, IERX)
         IF (IERX .EQ. 0) THEN
            IXI = IMEM(IHDR+2) + 1
            IF (IXTS .GT. 0 .AND. IXTS .LE. IMEM(IXI)) THEN
C              Don't reference IMEM(IXI+IXTS) if IXTS is not in range
               IF (IMEM(IXI+IXTS) .EQ. ITS) THEN
                  IMEM(IXI+IXTS) = 0
                  IJS = 0
               ENDIF
            ENDIF
            IF (IJS .NE. 0) THEN
C              Check all Trimmed Surface locations in the Joined Surface
               DO 100 I=IXI+1,IXI+IMEM(IXI)
                  IF (IMEM(I) .EQ. ITS)  IMEM(I) = 0
  100          CONTINUE
            ENDIF
         ENDIF
      ENDIF
C     Now erase any and all component Loop entities, ignoring errors.
      CALL DTERPT (0)
      DO 200 I=IXT+5,IXT+4+NLP
         ILP = IMEM(I)
         IF (ILP .NE. 0)  CALL D2LPER (CMEM, IMEM, DMEM, ILP, IERX)
  200 CONTINUE
C     Now erase the B-spline surface, if any, ignoring errors.
      IF (ISBF .NE. 0)  CALL D2ERAS (CMEM, IMEM, DMEM, ISBF, IERX)
C     Finally, erase the Trimmed Surface entity itself, ignoring errors.
      CALL D2ERAS (CMEM, IMEM, DMEM, ITS, IERX)
      CALL DTERPT (1)
      RETURN                    

C     Error Handling

 9000 CONTINUE
C     Report error
      CALL DTERR (1, SUBNAM, IER, 0)
      RETURN
      END
