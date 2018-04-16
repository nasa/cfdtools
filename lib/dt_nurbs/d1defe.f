      SUBROUTINE D1DEFE (CMEM, IMEM, DMEM, ITYP, LENC, LENI, LEND,
     +    INITIZ, IDE, IER)

C   Do the work for D2DEFE and other allocation routines.  Check only to
C   ensure sufficient space in CMEM, IMEM, and DMEM.  If not, report error
C   under the name of the caller, which should be SUBNAM.
C
C   DYNAMIC MEMORY INPUT-OUTPUT:
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     ITYP    Data type id number for new entity
C     LENC    Character data length for new entity
C     LENI    Integer data length for new entity
C     LEND    Double precision data length for new entity
C     INITIZ  Logical value indicating whether to initialize the data
C             spaces of the entity to ' ', 0, and 0.0D0, respectively.
C     SUBNAM  Subroutine name string of the calling routine
C   OUTPUT:
C     IDE     Pointer to new data entity
C     IER     Returned error value.  Zero implies no errors.  Otherwise,
C             = -2  Insufficient space available for character data
C             = -3  Insufficient space available for integer data
C             = -4  Insufficient space available for double precision data
C
C   HISTORY:
C     03Feb92 P Kraushar    Created.
C     17Mar92 D. Parsons    Remove call to DTERR
C     22Mar92 D. Parsons    Adjusted IER values
C*****
C
      CHARACTER CMEM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)
      INTEGER ITYP, LENC, LENI, LEND, IDE, I, IER
      LOGICAL INITIZ

      INTEGER IHDR, NAVLBL
C****
C     Check for sufficient free space
      NAVLBL = IMEM(1) - IMEM(4) - LENC
      IF (NAVLBL .LT. 0) GO TO 9001
      NAVLBL = IMEM(16) - IMEM(5) - LENI
      IF (IMEM(18) .EQ. 0) NAVLBL = NAVLBL - 8
      IF (NAVLBL .LT. 0) GO TO 9002
      NAVLBL = IMEM(3) - IMEM(6) - LEND
      IF (NAVLBL .LT. 0) GO TO 9003

      IF (IMEM(18) .EQ. 0) THEN
C       Create new header block
        IMEM(16) = IMEM(16) - 8
        IHDR = IMEM(16)
        IMEM(19) = IMEM(19) + 1
        IF (IMEM(19) .GE. 256) IMEM(19) = 0
        IDE = IHDR*256 + IMEM(19)
      ELSE
C       Reuse existing free header block
        IHDR = IMEM(18)
        IMEM(18) = IMEM(IHDR+7)
        IDE = ABS(IMEM(IHDR)) + 1
        IF (IDE/256 .NE. IHDR) IDE = IHDR*256
      END IF
C     Initialize header block
      IMEM(IHDR) = IDE
      IMEM(IHDR+1) = IMEM(4)
      IMEM(IHDR+2) = IMEM(5)
      IMEM(IHDR+3) = IMEM(6)
      IMEM(IHDR+4) = LENC
      IMEM(IHDR+5) = LENI
      IMEM(IHDR+6) = LEND
      IMEM(IHDR+7) = IMEM(17)*256 + ITYP
C     Update master parameter block
      IMEM(17) = IHDR
      IMEM(4) = IMEM(4) + LENC
      IMEM(5) = IMEM(5) + LENI
      IMEM(6) = IMEM(6) + LEND
C     Initialize new entity, if requested
      IF (INITIZ) THEN
        IF (LENC .GT. 0) CMEM(IMEM(IHDR+1):IMEM(4)-1) = ' '
        DO 100 I = IMEM(IHDR+2),IMEM(5)-1
          IMEM(I) = 0
  100   CONTINUE
        DO 110 I = IMEM(IHDR+3),IMEM(6)-1
          DMEM(I) = 0.0D0
  110   CONTINUE
      END IF

      IER = 0
      RETURN

C   Error reporting section

 9001 IER = -2
      GO TO 9900

 9002 IER = -3
      GO TO 9900

 9003 IER = -4
      GO TO 9900

 9900 CONTINUE

      RETURN
      END
