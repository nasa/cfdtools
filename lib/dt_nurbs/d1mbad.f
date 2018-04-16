      LOGICAL FUNCTION D1MBAD (CMEM, IMEM, DMEM, MEMAXC, MEMAXI,
     +    MEMAXD, IER)

C   Test whether dynamic memory structure is bad.  IER indicates the
C   nature of the failure.
C
C   DYNAMIC MEMORY INPUT-OUTPUT:
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     MEMAXC  If positive, the correct size of CMEM.  Otherwise, ignored.
C     MEMAXI  If positive, the correct size of IMEM.  Otherwise, ignored.
C     MEMAXD  If positive, the correct size of DMEM.  Otherwise, ignored.
C   OUTPUT:
C     D1MBAD   .TRUE. implies that the memory is bad, see IER for an
C                 explanation
C              .FALSE. implies that no errors were found
C
C     IER     Returned error value.  Zero implies no errors.  Triples refer
C             to same defect in CMEM, IMEM, or DMEM, respectively.
C
C     <NOTE THAT ERRORS -1 THROUGH -44 ARE RETURNED FROM D0MCHK>
C
C    -1,3   = Internal length doesn't match given length
C    -4,6   = Internal length not in valid range
C    -7     = Check value at end of IMEM is changed
C    -8     = Error check level flag at beginning of CMEM has invalid value
C    -9     = Check value at end of CMEM is changed
C    -10    = Check value at beginning of DMEM is changed
C    -11    = Check value at end of DMEM is changed
C    -12    = Lowest header (IMEM(16)) is out of range
C    -13    = Lowest header (IMEM(16)) is invalid (MOD 8 test)
C    -14,16 = Lowest free location is out of range
C    -17,19 = Header of highest lock is out of range
C    -20,22 = Header of highest lock is invalid (MOD 8 test)
C    -23,25 = Header of highest lock is marked deleted or zero
C    -26,28 = Lowest unlocked location does not match end of highest lock
C    -29,31 = Lowest search location is out of range
C    -32,34 = Lowest unlocked location is not at default when no locks
C    -35,37 = Lowest search location is not at default when no locks
C    -38    = Top of active list is out of range
C    -39    = Top of active list is invalid (MOD 8 test)
C    -40    = Pointer check value for top of active list is wrong
C    -41    = Top of free list is out of range
C    -42    = Top of free list is invalid (MOD 8 test)
C    -43    = Pointer check value for top of free list is wrong
C    -44    = Last new header sequence number is out of range (0-255)
C
C    -45    = Pointer check value in header block is wrong
C    -46    = List link value in header block is negative
C    -47,49 = Data space beginning location is out of range
C    -50,52 = Locked data space overlaps free space
C    -53,55 = Locked data space below lowest search location
C    -56,58 = Deleted data space below lowest search location
C    -59    = List link in header out of range
C    -60    = List link in header invalid (MOD 8 test)
C    -61    = Free header block found in active list
C    -62    = Active list loops back on itself
C    -63    = Character data spaces out of order or have gaps
C    -64    = Integer data spaces out of order or have gaps
C    -65    = Double precision data spaces out of order or have gaps
C    -66    = Active list ends before reaching all active header blocks
C    -67    = Active list continues too far
C    -68    = Active header block found in free list
C    -69    = Free list loops back on itself or enters active list
C    -70    = Free list ends before reaching all free header blocks
C    -71    = Free list continues too far
C
C   HISTORY:
C     11Feb92  P Kraushar     Created.
C      3Mar92  D. Parsons     Converted to a Lower-level routine
C                             from D2MBAD
C     11Aug92  D. Parsons     Split into two pieces (D0MCHK now checks
C                             the memory management section)
C*****
C
      CHARACTER CMEM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)

      INTEGER MEMAXC, MEMAXI, MEMAXD, IER

      INTEGER IHDR, NHDR, NACTV, NFREE, MEMIN(3), K, I

      DATA MEMIN / 2, 20, 2 /
C****
      D1MBAD = .FALSE.
      IER    = 0

C     Check the memory management section

      CALL D0MCHK (CMEM, IMEM, DMEM, MEMAXC, MEMAXI, MEMAXD, IER)
      IF (IER .NE. 0) GOTO 9900

C     Check individual headers

      NACTV = 0
      NFREE = 0

      DO 300 I=IMEM(16),IMEM(2)-8,8
        IF (ABS(IMEM(I))/256 .NE. I) GO TO 9045
        IF (IMEM(I+7) .LT. 0) GO TO 9046
        IF (IMEM(I).LT.0 .AND. IMEM(I+4).EQ.0 .AND. IMEM(I+5).EQ.0
     +      .AND. IMEM(I+6).EQ.0) THEN
          NFREE = NFREE + 1
          NHDR = IMEM(I+7)
        ELSE
          NACTV = NACTV + 1
          NHDR = IMEM(I+7)/256
          DO 250 K=1,3
            IF (IMEM(I+K).LT.MEMIN(K) .OR. IMEM(I+K).GT.IMEM(3+K))
     +          GO TO 9047
            IF (IMEM(I+3+K) .LT. 0) THEN
              IF (IMEM(I+K) + ABS(IMEM(I+3+K)) .GT. IMEM(6+K))
     +            GO TO 9050
              IF (IMEM(I+K) .LT. IMEM(12+K)) GO TO 9053
            END IF
            IF (IMEM(I) .LT. 0) THEN
              IF ((IMEM(I+K) .LT. IMEM(12+K))
     +            .AND. (IMEM(I+3+K) .NE. 0)) GO TO 9056
            END IF
  250     CONTINUE
        END IF
        IF (NHDR .NE. 0) THEN
          IF (NHDR.LT.IMEM(16) .OR. NHDR.GE.IMEM(2)) GO TO 9059
          IF (MOD (IMEM(2)-NHDR, 8) .NE. 0) GO TO 9060
        END IF
  300 CONTINUE

      NHDR = IMEM(17)
      DO 400 I=1,NACTV
        IHDR = NHDR
        IF (IMEM(IHDR).LT.0 .AND. IMEM(IHDR+4).EQ.0 .AND.
     +      IMEM(IHDR+5).EQ.0 .AND. IMEM(IHDR+6).EQ.0) GO TO 9061
        IF (IMEM(IHDR+7) .LT. 0) GO TO 9062
        NHDR = IMEM(IHDR+7)/256
        IMEM(IHDR+7) = -IMEM(IHDR+7)
        IF (NHDR .NE. 0) THEN
          IF (IMEM(NHDR+1) + ABS(IMEM(NHDR+4)) .NE. IMEM(IHDR+1))
     +        GO TO 9063
          IF (IMEM(NHDR+2) + ABS(IMEM(NHDR+5)) .NE. IMEM(IHDR+2))
     +        GO TO 9064
          IF (IMEM(NHDR+3) + ABS(IMEM(NHDR+6)) .NE. IMEM(IHDR+3))
     +        GO TO 9065
        ELSE
          IF (I .LT. NACTV) GO TO 9066
        END IF
  400 CONTINUE
      IF (NHDR .NE. 0) GO TO 9067

      NHDR = IMEM(18)
      DO 500 I=1,NFREE
        IHDR = NHDR
        IF (IMEM(IHDR).GE.0 .OR. IMEM(IHDR+4).NE.0 .OR.
     +      IMEM(IHDR+5).NE.0 .OR. IMEM(IHDR+6).NE.0) GO TO 9068
        IF (IMEM(IHDR+7) .LT. 0) GO TO 9069
        NHDR = IMEM(IHDR+7)
        IMEM(IHDR+7) = -IMEM(IHDR+7)
        IF (NHDR .EQ. 0) THEN
          IF (I .LT. NFREE) GO TO 9070
        END IF
  500 CONTINUE
      IF (NHDR .NE. 0) GO TO 9071

C     Restore list links to normal

      IER = 0

  600 CONTINUE

      DO 700 I=IMEM(16),IMEM(2)-8,8
        IMEM(I+7) = ABS(IMEM(I+7))
  700 CONTINUE

      IF (IER .NE. 0) GO TO 9900

      RETURN

C   Error reporting section
C     Triples distinguished by K refer to same defect in CMEM, IMEM, or DMEM

C     Pointer check value in header block is wrong
 9045 IER = -45
      GO TO 9900
C     List link value in header block is negative
 9046 IER = -46
      GO TO 9900
C     Data space beginning location is out of range
 9047 IER = -46 - K
      GO TO 9900
C     Locked data space overlaps free space
 9050 IER = -49 - K
      GO TO 9900
C     Locked data space below lowest search location
 9053 IER = -52 - K
      GO TO 9900
C     Deleted data space below lowest search location
 9056 IER = -55 - K
      GO TO 9900
C     List link in header out of range
 9059 IER = -59
      GO TO 9900
C     List link in header invalid (MOD 8 test)
 9060 IER = -60
      GO TO 9900

C   Next series of errors must restore header blocks before exiting

C     Free header block found in active list
 9061 IER = -61
      GO TO 600

C     Active list loops back on itself
 9062 IER = -62
      GO TO 600

C     Character data spaces out of order or have gaps
 9063 IER = -63
      GO TO 600

C     Integer data spaces out of order or have gaps
 9064 IER = -64
      GO TO 600

C     Double precision data spaces out of order or have gaps
 9065 IER = -65
      GO TO 600

C     Active list ends before reaching all active header blocks
 9066 IER = -66
      GO TO 600

C     Active list continues too far
 9067 IER = -67
      GO TO 600

C     Active header block found in free list
 9068 IER = -68
      GO TO 600

C     Free list loops back on itself or enters active list
 9069 IER = -69
      GO TO 600

C     Free list ends before reaching all free header blocks
 9070 IER = -70
      GO TO 600

C     Free list continues too far
 9071 IER = -71
      GO TO 600

 9900 CONTINUE

      IF (IER .NE. 0) THEN
         D1MBAD = .TRUE.
      ENDIF

      RETURN
      END
