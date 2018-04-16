      SUBROUTINE D2DEFX (CMEM, IMEM, DMEM, IDE, MOREC, MOREI, MORED,
     +    INITIZ, IER)

C   Extend (or shrink) the data spaces allocated to an existing object.
C   This action fails if any of the data spaces being changed is within
C   the locked portion of its respective array.  This routine should be
C   avoided unless the data entity is very recently allocated.
C
C   DYNAMIC MEMORY INPUT-OUTPUT:
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     IDE     Pointer to entity to be extended (or shrunk)
C     MOREC   Additional characters (may be negative)
C     MOREI   Additional integers (may be negative)
C     MORED   Additional double precision values (may be negative)
C     INITIZ  Logical variable indicating whether to initialize any added
C             space to blanks or zeroes, as appropriate
C   OUTPUT:
C     IER     Returned error value.  Zero implies no errors.  Otherwise,
C             -1 = Dynamic memory is corrupt or uninitialized
C             -2 = Null pointer
C             -3 = IDE does not hold a pointer
C             -4 = IDE is either not a pointer or points to a deleted
C                  entity (ambiguous case)
C             -5 = IDE points to a deleted entity
C             -8 = Character data space is locked, cannot change it
C             -9 = Integer data space is locked, cannot change it
C            -10 = Double precision data space is locked, cannot
C                  change it
C            -11 = Insufficient character space available
C            -12 = Insufficient integer space available
C            -13 = Insufficient double precision space available
C            -14 = Attempt to reduce character space by more than LENC
C            -15 = Attempt to reduce integer space by more than LENI
C            -16 = Attempt to reduce double precision space by more
C                  than LEND
C
C   CALLS:
C     D0PTR    Pointer check
C     D1DEFX   Entity extend
C     DTERR    Error Handler
C
C   HISTORY:
C     04Feb92 P Kraushar    Created.
C     25Mar92 D. Parsons    Add error returns for over-shrinking
C     01Jun93 D. Parsons    Extracted guts to D1DEFX
C*****
C
C   Long Name Alias:
C     ENTRY D2_DEFINITION_EXTEND (CMEM, IMEM, DMEM, IDE, MOREC, MOREI,
C    +    MORED, INITIZ, IER)

      CHARACTER CMEM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)

      INTEGER IDE, MOREC, MOREI, MORED, IER
      LOGICAL INITIZ

      INTEGER IHDR, NEED, IDUM

      CHARACTER*6 SUBNAM
      DATA SUBNAM /'D2DEFX'/

C****

C     Decode pointer and test for lock state and space availability
      CALL D0PTR (CMEM, IMEM, DMEM, IDE, -1, 0, IHDR, IDUM, IER)
      IF (IER .NE. 0) GO TO 9900
C
      CALL D1DEFX (CMEM, IMEM, DMEM, IHDR, MOREC, MOREI, MORED,
     +      INITIZ, NEED, IER)
C
 9900 CONTINUE

      IF ( (IER .EQ. -11) .OR. (IER .EQ. -12) .OR. (IER .EQ. -13)) THEN
         CALL DTERR (2, SUBNAM, IER, NEED)
      ELSE IF (IER .NE. 0) THEN
         CALL DTERR (1, SUBNAM, IER, 0)
      ENDIF

      RETURN
      END
