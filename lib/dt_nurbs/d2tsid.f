      LOGICAL FUNCTION D2TSID (CMEM, IMEM, DMEM, IDE)

C   Test for implicit lock on the double precision space belonging
C     to the data entity pointed to by IDE.
C
C   DYNAMIC MEMORY INPUT-OUTPUT:
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     IDE     Pointer to data entity to be unlock
C   OUTPUT:
C     D2TSID   True if the double precision space is implicitly locked
C              False if the double precision space is not locked
C
C              Also, True if the IDE is invalid, or points to an
C              invalid entity.
C
C   CALLS:    D0PTR
C
C   HISTORY:
C     10March92   D. Parsons  created
C     29April92   D. Parsons  Set if IDE invalid
C*****
C
C   Long Name Alias:
C     ENTRY D2_TEST_IMPLICIT_LOCK_DMEM (CMEM, IMEM, DMEM, IDE)

      CHARACTER CMEM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)
      INTEGER IDE, IER

      INTEGER IHDR, IDUM

C****
C     Decode pointer

      CALL D0PTR (CMEM, IMEM, DMEM, IDE, -1, 0, IHDR, IDUM, IER)
      IF (IER .NE. 0) GO TO 9900

C     Test for locked double precision data space

      D2TSID = (IMEM(IHDR+3) .LT. IMEM(9))

      RETURN

C   Error section

 9900 D2TSID = .TRUE.

      RETURN
      END
