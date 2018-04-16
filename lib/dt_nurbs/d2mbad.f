      LOGICAL FUNCTION D2MBAD (CMEM, IMEM, DMEM, MEMAXC, MEMAXI,
     +    MEMAXD, IER)

C   Test whether dynamic memory structure is bad.
C   If IER= 0, then there are no detected problems with the dynamic
C   memory structure, otherwise IER indicates the nature of the failure.
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
C   CALLS:    DTERR, D1MBAD
C   HISTORY:
C     11Feb92  P. Kraushar    Created.
C     3Mar92   D. Parsons     Changed so that most code is in a
C                             Lower-level routine (D1MBAD)
C*****
C
C   Long Name Alias:
C     ENTRY D2_MEMORY_BAD (CMEM, IMEM, DMEM, MEMAXC, MEMAXI,
C    +    MEMAXD, IER)

      EXTERNAL D1MBAD
      LOGICAL  D1MBAD

      CHARACTER SUBNAM*(6)
      CHARACTER CMEM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)
      INTEGER MEMAXC, MEMAXI, MEMAXD, IER

      DATA SUBNAM /'D2MBAD'/

      D2MBAD = D1MBAD (CMEM, IMEM, DMEM, MEMAXC, MEMAXI, MEMAXD, IER)
      IF (IER .NE. 0) THEN
         CALL DTERR (1, SUBNAM, IER, 0)
      ENDIF

      RETURN
      END
