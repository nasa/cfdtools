      SUBROUTINE D0ERAS (CMEM, IMEM, DMEM, IHDR)

C   Dynamic memory utility which carries out the actual erasure process
C   on the entity with header index IHDR.  Collects what it can of the
C   three data spaces.  Also collects the header itself, when possible,
C   or puts it in the free header list when all its data spaces are collected.
C
C   DYNAMIC MEMORY INPUT-OUTPUT:
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     IHDR    Index of header block for entity to be erased
C   OUTPUT:
C
C   CALLS:    D0AHDR
C   HISTORY:
C     04Feb92 P Kraushar   Created.
C     09Mar92 D. Parsons   Corrected error in substring ranges
C
C*****
C
      CHARACTER CMEM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)
      INTEGER IHDR

      INTEGER LENC, LENI, LEND, IBEG, IEND, LHDR, KHDR, I
      LOGICAL COLHDR
C****
      COLHDR = .TRUE.
      LENC = IMEM(IHDR+4)
      LENI = IMEM(IHDR+5)
      LEND = IMEM(IHDR+6)
C   Collect character data space, if present and unlocked (and used).
      IF (LENC .GT. 0) THEN
        IF (IMEM(IHDR+1) .GE. IMEM(7)) THEN
          IBEG = IMEM(IHDR+1) + LENC
          IEND = IMEM(4) - 1
          IF (IEND .GE. IBEG)
     +         CMEM(IBEG-LENC:IEND-LENC) = CMEM(IBEG:IEND)
        ELSE
          LENC = 0
          COLHDR = .FALSE.
        END IF
      END IF
C   Collect integer data space, if present and unlocked.
      IF (LENI .GT. 0) THEN
        IF (IMEM(IHDR+2) .GE. IMEM(8)) THEN
          IBEG = IMEM(IHDR+2) + LENI
          IEND = IMEM(5) - 1
          DO 200 I=IBEG,IEND
            IMEM(I-LENI) = IMEM(I)
  200     CONTINUE
        ELSE
          LENI = 0
          COLHDR = .FALSE.
        END IF
      END IF
C   Collect double precision data space, if present and unlocked.
      IF (LEND .GT. 0) THEN
        IF (IMEM(IHDR+3) .GE. IMEM(9)) THEN
          IBEG = IMEM(IHDR+3) + LEND
          IEND = IMEM(6) - 1
          DO 300 I=IBEG,IEND
            DMEM(I-LEND) = DMEM(I)
  300     CONTINUE
        ELSE
          LEND = 0
          COLHDR = .FALSE.
        END IF
      END IF
C   Now adjust all the preceding header blocks to reflect the shifts in
C   the data spaces.
      CALL D0AHDR (CMEM, IMEM, DMEM, IHDR, LENC, LENI, LEND, LHDR)
C   Collect the header block itself, if possible
      IF (COLHDR) THEN
        KHDR = IMEM(IHDR+7)/256
C       Remove it from the active list
        IF (LHDR .EQ.0) THEN
          IMEM(17) = KHDR
        ELSE
          IMEM(LHDR+7) = KHDR*256 + MOD(IMEM(LHDR+7),256)
        END IF
C       Collect it if it is the lowest header block, otherwise put it
C       on the free header list.
        IF (IHDR .EQ. IMEM(16)) THEN
          IMEM(16) = IMEM(16) + 8
          IMEM(IHDR) = 0
          IMEM(IHDR+7) = 0
        ELSE
          IMEM(IHDR+7) = IMEM(18)
          IMEM(18) = IHDR
        END IF
      END IF
C   Collect the lowest header block, if it is on the free header list.
      KHDR = IMEM(16)
      IF (IMEM(KHDR) .LT. 0) THEN
        IF (IMEM(KHDR+4) .EQ. 0 .AND. IMEM(KHDR+5) .EQ. 0 .AND.
     +      IMEM(KHDR+6) .EQ. 0) THEN
          IF (IMEM(18) .EQ. KHDR) THEN
            IMEM(18) = IMEM(KHDR+7)
          ELSE
            LHDR = IMEM(18)
            KHDR = IMEM(LHDR+7)
            DO 400 I=1,1+(IMEM(2)-IMEM(16))/8
              IF (KHDR .EQ. IMEM(16)) GO TO 410
              LHDR = KHDR
              KHDR = IMEM(LHDR+7)
  400       CONTINUE
C           The DO loop bound is chosen so that it CANNOT terminate here.
            STOP 'Dynamic memory fatally corrupted - D0ERAS'
  410       CONTINUE
            IMEM(LHDR+7) = IMEM(KHDR+7)
          END IF
          IMEM(16) = IMEM(16) + 8
          IMEM(KHDR) = 0
          IMEM(KHDR+7) = 0
        END IF
      END IF

      RETURN
      END
