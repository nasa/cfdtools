      SUBROUTINE D2EGER (CMEM, IMEM, DMEM, IEG, IER)

C     Erase an Edge entity, disconnecting from any joined edge or Loop
C     and erasing the associated B-spline curve.
C
C   DYNAMIC MEMORY INPUT-OUTPUT
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     IEG     MEM pointer to Edge entity
C   OUTPUT:
C     IER     Returned error code.  Zero implies no errors.
C             -1  = Dynamic memory is corrupt or uninitialized.
C             -3  = IEG is a garbage value, not a pointer.
C             -6  = IEG does not point to an Edge entity.
C
C   CALLS:
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
C     ENTRY D2_EDGE_ERASE (CMEM, IMEM, DMEM, IEG, IER)
C
      CHARACTER CMEM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)
      INTEGER IEG, IER

      INTEGER ENTEG
      PARAMETER (ENTEG=243)
      INTEGER ENTLP
      PARAMETER (ENTLP=242)

      INTEGER IHDR, ITYP, IERX, ICBF, IXI, ILP, IXEG, JEG, I
      CHARACTER SUBNAM*6

      DATA SUBNAM /'D2EGER'/

C****
      IER = 0
      IF (IEG .EQ. 0)  RETURN
C     Check pointer and type of IEG
      CALL D0PTR (CMEM, IMEM, DMEM, IEG, ENTEG, 0, IHDR, ITYP, IERX)
      IF (IERX .LT. 0) THEN
         IF (IERX .EQ. -5 .OR. IERX .EQ. -4)  RETURN
         IER = IERX
         GOTO 9000
      ENDIF
      IXI = IMEM(IHDR+2)
      ICBF = IMEM(IXI)
      ILP  = IMEM(IXI+1)
      IXEG = IMEM(IXI+2)
      JEG  = IMEM(IXI+3)
C     If there is a joined edge, disconnect it.  Ignore bad data.
      IF (JEG .NE. 0) THEN
         CALL D0PTR (CMEM, IMEM, DMEM, JEG, ENTEG, 0, IHDR, ITYP, IERX)
         IF (IERX .EQ. 0) THEN
            IXI = IMEM(IHDR+1) + 3
            IF (IMEM(IXI) .EQ. IEG)  IMEM(IXI) = 0
         ENDIF
      ENDIF
C     If part of a loop, disconnect it.  Ignore bad data.
      IF (ILP .NE. 0) THEN
         CALL D0PTR (CMEM, IMEM, DMEM, ILP, ENTLP, 0, IHDR, ITYP, IERX)
         IF (IERX .EQ. 0) THEN
            IXI = IMEM(IHDR+2) + 2
            IF (IXEG .GT. 0 .AND. IXEG .LE. IMEM(IXI)) THEN
C              Don't reference IMEM(IXI+IXEG) if IXEG is not in range
               IF (IMEM(IXI+IXEG) .EQ. IEG) THEN
                  IMEM(IXI+IXEG) = 0
                  ILP = 0
               ENDIF
            ENDIF
            IF (ILP .NE. 0) THEN
C              Check all Edge locations in the Loop entity
               DO 100 I=IXI+1,IXI+IMEM(IXI)
                  IF (IMEM(I) .EQ. IEG)  IMEM(I) = 0
  100          CONTINUE
            ENDIF
         ENDIF
      ENDIF
C     Now erase the B-spline curve, and then the Edge entity, ignoring errors.
      CALL DTERPT (0)
      IF (ICBF .NE. 0)  CALL D2ERAS (CMEM, IMEM, DMEM, ICBF, IERX)
      CALL D2ERAS (CMEM, IMEM, DMEM, IEG, IERX)
      CALL DTERPT (1)
      RETURN                    

C     Error Handling

 9000 CONTINUE
C     Report error
      CALL DTERR (1, SUBNAM, IER, 0)
      RETURN
      END
