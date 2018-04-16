      SUBROUTINE D2ERAS (CMEM, IMEM, DMEM, IDE, IER)

C   Erase the data entity pointed to by IDE.  Unlock it.  Reclaim its data
C   spaces for reuse, if unlocked.  Otherwise, mark it deleted.  The value
C   of IDE does not change, but it ceases to point to anything.
C
C   DYNAMIC MEMORY INPUT-OUTPUT:
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     IDE     Pointer to data entity to be erased
C   OUTPUT:
C     IER     Returned error value.  Zero implies no errors.  Otherwise,
C             -1 = Dynamic memory is corrupt or uninitialized.
C             -2 = Null pointer
C             -3 = IDE does not hold a pointer
C                = IDE is either not a pointer or points to a deleted
C             -4   entity (ambiguous case)
C             -5 = IDE points to a deleted entity
C             d2
C   CALLS:    D0ERAS, D0PTR, D0UNLK, DTERR
C   HISTORY:
C     28Jan92 P Kraushar    Created.
C*****
C
C   Long Name Alias:
C     ENTRY D2_ERASE_ENTITIY (CMEM, IMEM, DMEM, IDE, IER)

      CHARACTER CMEM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)

      INTEGER IDE, IER

      INTEGER IHDR, IDUM, KHDR, NHDR, LOWC, LOWI, LOWD, K, I

      CHARACTER*6 SUBNAM
      DATA SUBNAM /'D2ERAS'/

C****

      IER = 0

C     Decode pointer and test for lock state and space availability

      CALL D0PTR (CMEM, IMEM, DMEM, IDE, -1, 0, IHDR, IDUM, IER)
      IF (IER .NE. 0) GO TO 9900

      CALL D0UNLK (CMEM, IMEM, DMEM, IHDR, 1, LOWC)
      CALL D0UNLK (CMEM, IMEM, DMEM, IHDR, 2, LOWI)
      CALL D0UNLK (CMEM, IMEM, DMEM, IHDR, 3, LOWD)

C     If the data space is non-null, and cannot be collected right
C     now, ensure it is included in the area to be searched for future
C     collecting

      DO 100 K = 1, 3
         IF (IMEM(IHDR+3+K) .GT. 0 .AND. IMEM(IHDR+K) .LT. IMEM(6+K))
     +         IMEM(12+K) = MIN (IMEM(12+K), IMEM(IHDR+K))
  100 CONTINUE

C     Mark it erased, note the next header, then try to erase it.

      IMEM(IHDR) = -IMEM(IHDR)
      NHDR = IMEM(IHDR+7)/256
      CALL D0ERAS (CMEM, IMEM, DMEM, IHDR)

C     Did the unlocking expose other, previously deleted entities to
C     at least partial garbage collection?

      IF (LOWC .LT. IMEM(4) .OR. LOWI .LT. IMEM(5) .OR.
     +    LOWD .LT. IMEM(6)) THEN
        DO 110 I=1,1+(IMEM(2)-IMEM(16))/8
          KHDR = NHDR
          NHDR = IMEM(KHDR+7)/256
          IF (IMEM(KHDR) .LT. 0) CALL D0ERAS (CMEM, IMEM, DMEM, KHDR)
          IF (NHDR .EQ. 0) GO TO 120
          IF (IMEM(NHDR+1) .LT. LOWC .AND. IMEM(NHDR+2) .LT. LOWI
     +        .AND. IMEM(NHDR+3) .LT. LOWD) GO TO 120
  110   CONTINUE

C       The DO loop bound is chosen so that it CANNOT terminate here.
        STOP 'Dynamic memory fatally corrupted - D2ERAS'

  120   CONTINUE
      END IF

      IER = 0
      RETURN

C   Error reporting section

 9900 CALL DTERR (1, SUBNAM, IER, 0)
      RETURN
      END
