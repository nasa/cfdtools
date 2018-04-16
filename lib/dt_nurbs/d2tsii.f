      LOGICAL FUNCTION D2TSII (CMEM, IMEM, DMEM, IDE)

C   Test for implicit lock on the integer space belonging to the
C     data entity pointed to by IDE.
C
C   DYNAMIC MEMORY INPUT-OUTPUT:
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     IDE     Pointer to data entity to be unlock
C   OUTPUT:
C     D2TSII   True if the integer space is implicitly locked
C              False if the integer space is not locked
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
C     ENTRY D2_TEST_IMPLICIT_LOCK_IMEM (CMEM, IMEM, DMEM, IDE)

      CHARACTER CMEM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)
      INTEGER IDE, IER

      INTEGER IHDR, IDUM

C****
C     Decode pointer

      CALL D0PTR (CMEM, IMEM, DMEM, IDE, -1, 0, IHDR, IDUM, IER)
      IF (IER .NE. 0) GO TO 9900

C     Test for locked integer data space

      D2TSII = (IMEM(IHDR+2) .LT. IMEM(8))

      RETURN

C   Error section

 9900 D2TSII = .TRUE.

      RETURN
      END
