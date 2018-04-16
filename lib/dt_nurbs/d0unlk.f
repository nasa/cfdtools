      SUBROUTINE D0UNLK (CMEM, IMEM, DMEM, IHDR, K, LOWER)

C   If locked, unlock the Kth data space of IHDR.  Find next lowest unlocked
C   location in that data space.  Also report lowest newly reclaimable
C   space found in LOWER.
C
C   DYNAMIC MEMORY INPUT-OUTPUT:
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     IHDR    Index of header block of entity to be unlocked
C     K       Data space to be unlocked
C             1 = Character data space
C             2 = Integer data space
C             3 = Double precision data space
C   OUTPUT:
C     LOWER   Lowest newly reclaimable data location in the data space.
C             This occurs when deleted entities were locked by being below
C             the formerly locked data space of IHDR.
C
C   CALLS:    none
C   HISTORY:
C     04Feb92 P Kraushar    Created.
C*****
C
      CHARACTER CMEM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)
      INTEGER IHDR, K, LOWER, I

      INTEGER KHDR, NHDR, MEMIN(3)

      DATA MEMIN / 2, 20, 2 /
C****
C     Default value for LOWER is current free space limit
      LOWER = IMEM(3+K)

      IF (IMEM(IHDR+3+K) .LT. 0) THEN
C       It was locked, so unlock it
        IMEM(IHDR+3+K) = -IMEM(IHDR+3+K)
        IF (IHDR .EQ. IMEM(9+K)) THEN
C         It was also the highest locked entity in the Kth data space,
C         so we need to search for the new highest locked entity.
          KHDR = IHDR
          DO 100 I=1,1+(IMEM(2)-IMEM(16))/8
            NHDR = IMEM(KHDR+7)/256
            IF (IMEM(KHDR) .LT. 0) THEN
C             This one is marked deleted, did it have space to reclaim?
              IF (IMEM(KHDR+3+K) .NE. 0) LOWER = IMEM(KHDR+K)
            END IF
            IF (IMEM(KHDR+K) .LE. IMEM(12+K)) THEN
C             End of search, this is the lowest one that might have been
C             locked or deleted.
              IMEM(6+K) = MEMIN(K)
              IMEM(9+K) = 0
              IMEM(12+K) = IMEM(K)
              GO TO 110
            ELSE IF (IMEM(NHDR+3+K) .LT. 0) THEN
C             End of search, the next one is locked.
              IMEM(6+K) = IMEM(KHDR+K)
              IMEM(9+K) = NHDR
              GO TO 110
            END IF
            KHDR = NHDR
  100     CONTINUE
C         The DO loop bound is chosen so that it CANNOT terminate here.
          STOP 'Dynamic memory fatally corrupted - D0UNLK'
  110     CONTINUE
        END IF
      END IF

      RETURN
      END
