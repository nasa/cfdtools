      SUBROUTINE D0NXDE (CMEM, IMEM, DMEM, IHDR, JDE)

C     Locate the next empty Directory Entry slot in the IGES Index with
C     header index IHDR after the (presumed used) entry JDE.  Attempt to
C     expand the IGES Index if necessary.  Return JDE = 0 if failed.
C
C   DYNAMIC MEMORY INPUT-OUTPUT
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     IHDR    Header index to IGES Index
C   INPUT-OUTPUT:
C     JDE     On Input: Last IGES directory entry line number used.
C             On Output: Next unused directory entry line number, or
C             zero if there are no more.  JDE = -1 is valid input.
C
C   CALLS:    DTERPT, D2IGXI
C
C   HISTORY:
C     24Aug93  P Kraushar   Created
C*****
      CHARACTER CMEM*(*)
      INTEGER IMEM(*), IHDR, JDE
      DOUBLE PRECISION DMEM(*)

      INTEGER IXI0, IX
C****
      IXI0 = IMEM(IHDR+2) - 1
C     Start search at next entry (JDE+2), end with last entry in Index
      DO 100 IX = IXI0+(JDE+7)/2, IXI0+ABS(IMEM(IHDR+5))
         IF (IMEM(IX) .EQ. 0) THEN
            JDE = 2*(IX - IXI0) - 5
C           Update max entry in Index?
            IMEM(IXI0+2) = MAX (IMEM(IXI0+2), JDE+1)
            RETURN
         ENDIF
  100 CONTINUE
      JDE = 2*ABS(IMEM(IHDR+5)) - 3
      CALL DTERPT (0)
      CALL D2IGXI (CMEM, IMEM, DMEM, IMEM(IHDR), 0, 40, IX)
      CALL DTERPT (1)
      IF (IX .NE. 0)  JDE = 0
      RETURN
      END

