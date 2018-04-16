      SUBROUTINE D0PTR (CMEM, IMEM, DMEM, IPTR, JTYP, IER0, IHDR,
     +    ITYP, IER)

C   Extract the header index from a pointer.  Check pointer validity.
C   Optionally extract and check the data type id number.  Utility routine
C   for dynamic memory subroutines.
C
C   DYNAMIC MEMORY INPUT-OUTPUT:
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     IPTR    Pointer to be checked and decoded
C     JTYP    Data type checking control
C             < 0  Ignore type checking.  Do not reference ITYP.
C             = 0  Report data type in ITYP.
C             > 0  Verify data type of entity is equal to JTYP.  Report
C                  as error if not.  Do not reference ITYP.
C     IER0    If errors occur, use error numbers negative from IER0 as
C             the zero point.  (For assigning different error numbers to
C             different erroneous pointer values in calling routine.)
C   OUTPUT:
C     IHDR    Index of header block corresponding to IPTR.
C     ITYP    Data type id when requested by JTYP = 0.
C     IER     Returned error value.  Zero implies no errors.  Otherwise,
C             = IER0 - (the natural error number for this routine)
C
C              -1    IMEM appears corrupted or uninitialized
C              -2    null pointer
C              -3    garbage value in pointer
C              -4    ambiguous--deleted entity or garbage pointer
C              -5    deleted entity
C              -6    entity is incorrect data type
C
C   CALLS:    D0MCHK
C
C   HISTORY:
C     28Jan92 P Kraushar    Created.
C*****
C
      CHARACTER CMEM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)
      INTEGER IPTR, JTYP, IER0, IHDR, ITYP, IER, IERX

      INTEGER KTYP
C****
C     Check the memory management section of the dynamic memory
      IF (CMEM(1:1) .NE. 'M' .AND. CMEM(1:1) .NE. 'L') THEN
         CALL D0MCHK (CMEM, IMEM, DMEM, 0, 0, 0, IERX)
         IF (IERX .NE. 0) GOTO 9001
      ENDIF
         
C     Null pointer is illegal in D0PTR
      IF (IPTR .EQ. 0) GO TO 9002

C     Extract header index from pointer value
      IHDR = IPTR/256

C     Do bounds checking first if error check level is high enough.
C     This prevents system aborts due to wild memory access violations.
      IF (CMEM(1:1) .NE. 'L') THEN
        IF (IHDR .LT. 20 .OR. IHDR .GE. IMEM(2)) GO TO 9003
        IF (CMEM(1:1) .NE. 'M') THEN
          IF ( MOD (IMEM(2)-IHDR, 8) .NE. 0) GO TO 9003
        END IF
        IF (IHDR .LT. IMEM(16)) GO TO 9004
      END IF

C     Always do basic pointer validity check
      IF (IMEM(IHDR) .NE. IPTR) THEN

C       Try to diagnose nature of failure, regardless of error check level.
        IF (IHDR .LT. 20 .OR. IHDR .GE. IMEM(2)) GO TO 9003
        IF (MOD (IMEM(2)-IHDR, 8) .NE. 0) GO TO 9003
        IF (IMEM(IHDR) .LT. 0) GO TO 9005
        IF (IMEM(IHDR)/256 .EQ. IHDR) GO TO 9004
        IF (IHDR .LT. IMEM(16)) GO TO 9004

        GO TO 9001

      END IF

C     Decode JTYP and do any requested data type id fetching or checking
      IF (JTYP .GE. 0) THEN
        KTYP = MOD( IMEM(IHDR+7), 256)
        IF (JTYP .EQ. 0) THEN
          ITYP = KTYP
        ELSE
          IF (JTYP .NE. KTYP) GO TO 9006
        END IF
      END IF

      IER = 0
      RETURN

C   Error reporting section

C     Dynamic memory array IMEM appears to be corrupted
 9001 IER = IER0 - 1
      RETURN

C     Null pointer
 9002 IER = IER0 - 2
      RETURN

C     Garbage value in pointer
 9003 IER = IER0 - 3
      RETURN

C     Ambiguous case: either a deleted entity or a garbage value
 9004 IER = IER0 - 4
      RETURN

C     Pointer refers to a deleted entity
 9005 IER = IER0 - 5
      RETURN

C     Entity is not of the correct data type
 9006 IER = IER0 - 6
      RETURN

      END
