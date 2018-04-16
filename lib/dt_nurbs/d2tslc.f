      LOGICAL FUNCTION D2TSLC (CMEM, IMEM, DMEM, IDE)

C   Test for explicit lock on the character space belonging to the
C     data entity pointed to by IDE.
C
C   DYNAMIC MEMORY INPUT-OUTPUT:
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     IDE     Pointer to data entity to be unlock
C   OUTPUT:
C     D2TSLC   True if the character space is explicitly locked
C              False if the character space is not locked
C              (space may be implicitly locked by another entity)
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
C     ENTRY D2_TEST_EXPLICIT_LOCK_CMEM (CMEM, IMEM, DMEM, IDE)

      CHARACTER CMEM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)
      INTEGER IDE, IER

      INTEGER IHDR, IDUM

C****
C     Decode pointer

      CALL D0PTR (CMEM, IMEM, DMEM, IDE, -1, 0, IHDR, IDUM, IER)
      IF (IER .NE. 0) GO TO 9900

C     Test for locked character data space

      D2TSLC = (IMEM(IHDR+4) .LT. 0)

      RETURN

C   Error section

 9900 D2TSLC = .TRUE.

      RETURN
      END
