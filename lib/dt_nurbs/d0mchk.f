      SUBROUTINE D0MCHK (CMEM, IMEM, DMEM, MEMAXC, MEMAXI, MEMAXD, IER)

C   Test whether dynamic memory structure is bad.
C
C   This routine only checks the memory management section.  It does
C   not check the individual entities.  See D1MBAD for a *COMPLETE*
C   check.
C
C   IER indicates the nature of the failure.
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
C     IER     Returned error value.  Zero implies no errors.  Triples refer
C             to same defect in CMEM, IMEM, or DMEM, respectively.
C     0     = No errors found
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
C   HISTORY:
C     11Aug92  D. Parsons     Extracted from D1MBAD
C*****

      CHARACTER CMEM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)

      INTEGER MEMAXC, MEMAXI, MEMAXD, IER

      INTEGER MX(3), MEMIN(3), K
      DATA MEMIN / 2, 20, 2 /

C****
      IER = 0

      MX(1) = MEMAXC
      MX(2) = MEMAXI
      MX(3) = MEMAXD
      DO 100 K=1,3
        IF (MX(K) .GE. MEMIN(K)) THEN
          IF (IMEM(K) .NE. MX(K)) GO TO 9001
        ELSE IF (IMEM(K) .LT. MEMIN(K)) THEN
          GO TO 9004
        END IF
        MX(K) = IMEM(K)
  100 CONTINUE

      IF (MX(1)+MX(2)+MX(3) .NE. IMEM(MX(2))) GO TO 9007
      IF (CMEM(1:1) .NE. 'L' .AND. CMEM(1:1) .NE. 'M' .AND.
     +    CMEM(1:1) .NE. 'H') GO TO 9008
      IF (CMEM(MX(1):MX(1)) .NE. '~') GO TO 9009
      IF (DMEM(1) .NE. 0.987654321D10) GO TO 9010
      IF (DMEM(MX(3)) .NE. 0.987654321D10) GO TO 9011

      IF (IMEM(16) .LT. MEMIN(2) .OR. IMEM(16) .GT. MX(2)) GO TO 9012
      IF (MOD (IMEM(2)-IMEM(16), 8) .NE. 0) GO TO 9013

      MX(2) = IMEM(16)
      DO 200 K=1,3
        IF (IMEM(3+K).LT.MEMIN(K) .OR. IMEM(3+K).GT.MX(K)) GO TO 9014
        IF (IMEM(9+K).NE.0) THEN
          IF (IMEM(9+K).LT.IMEM(16) .OR. IMEM(9+K).GT.IMEM(2)) GOTO 9017
          IF (MOD (IMEM(2)-IMEM(9+K), 8) .NE. 0) GO TO 9020
          IF (IMEM(IMEM(9+K)) .LE. 0) GO TO 9023
          IF (IMEM(6+K) .NE. IMEM(IMEM(9+K)+K) +
     +         ABS(IMEM(IMEM(9+K)+3+K)))
     +        GO TO 9026
          IF (IMEM(12+K).LT.MEMIN(K) .OR. IMEM(12+K).GE.IMEM(6+K))
     +        GO TO 9029
        ELSE
          IF (IMEM(6+K) .NE. MEMIN(K)) GO TO 9032
          IF (IMEM(12+K) .NE. IMEM(K)) GO TO 9035
        END IF
  200 CONTINUE

      IF (IMEM(17) .NE. 0) THEN
        IF (IMEM(17).LT.IMEM(16) .OR. IMEM(17).GE.IMEM(2)) GO TO 9038
        IF (MOD (IMEM(2)-IMEM(17), 8) .NE. 0) GO TO 9039
        IF (ABS (IMEM(IMEM(17)))/256 .NE. IMEM(17)) GO TO 9040
      END IF
      IF (IMEM(18) .NE. 0) THEN
        IF (IMEM(18).LT.IMEM(16) .OR. IMEM(18).GE.IMEM(2)) GO TO 9041
        IF (MOD (IMEM(2)-IMEM(18), 8) .NE. 0) GO TO 9042
        IF (ABS (IMEM(IMEM(18)))/256 .NE. IMEM(18)) GO TO 9043
      END IF
      IF (IMEM(19) .LT. 0 .OR. IMEM(19) .GT. 255) GO TO 9044

      RETURN

C   Error reporting section
C     Triples distinguished by K refer to same defect in CMEM, IMEM, 
C     or DMEM

C     Internal length doesn't match given length
 9001 IER = -K
      GO TO 9900

C     Internal length not in valid range
 9004 IER = -3 - K
      GO TO 9900

C     Check value at end of IMEM is changed
 9007 IER = -7
      GO TO 9900

C     Error check level flag at beginning of CMEM has invalid value
 9008 IER = -8
      GO TO 9900

C     Check value at end of CMEM is changed
 9009 IER = -9
      GO TO 9900

C     Check value at beginning of DMEM is changed
 9010 IER = -10
      GO TO 9900

C     Check value at end of DMEM is changed
 9011 IER = -11
      GO TO 9900

C     Lowest header (IMEM(16)) is out of range
 9012 IER = -12
      GO TO 9900

C     Lowest header (IMEM(16)) is invalid (MOD 8 test)
 9013 IER = -13
      GO TO 9900

C     Lowest free location is out of range
 9014 IER = -13 - K
      GO TO 9900

C     Header of highest lock is out of range
 9017 IER = -16 - K
      GO TO 9900

C     Header of highest lock is invalid (MOD 8 test)
 9020 IER = -19 - K
      GO TO 9900

C     Header of highest lock is marked deleted or zero
 9023 IER = -22 - K
      GO TO 9900

C     Lowest unlocked location does not match end of highest lock
 9026 IER = -25 - K
      GO TO 9900

C     Lowest search location is out of range
 9029 IER = -28 - K
      GO TO 9900

C     Lowest unlocked location is not at default when no locks
 9032 IER = -31 - K
      GO TO 9900

C     Lowest search location is not at default when no locks
 9035 IER = -34 - K
      GO TO 9900

C     Top of active list is out of range
 9038 IER = -38
      GO TO 9900

C     Top of active list is invalid (MOD 8 test)
 9039 IER = -39
      GO TO 9900

C     Pointer check value for top of active list is wrong
 9040 IER = -40
      GO TO 9900

C     Top of free list is out of range
 9041 IER = -41
      GO TO 9900

C     Top of free list is invalid (MOD 8 test)
 9042 IER = -42
      GO TO 9900

C     Pointer check value for top of free list is wrong
 9043 IER = -43
      GO TO 9900

C     Last new header sequence number is out of range (0-255)
 9044 IER = -44
      GO TO 9900

 9900 CONTINUE
      RETURN
      END
