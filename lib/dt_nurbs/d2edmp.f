      SUBROUTINE D2EDMP (CMEM, IMEM, DMEM, MEMAXC, MEMAXI, MEMAXD,
     +    LIMIT, LUNIT, IPTR)

C   Entity Dump for a dynamically allocated entity to a disk file.  A debugging
C   aid which attempts to make minimal assumptions about the consistency of
C   the dynamic memory arrays, while still dumping the entity contents in an
C   organized format.  Assumes the dump file is already open.
C
C   DYNAMIC MEMORY INPUT-OUTPUT:
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     MEMAXC  If positive, the actual size of CMEM.  Otherwise, D2EDMP will
C             assume the correct value is still in IMEM(1).
C     MEMAXI  If positive, the actual size of IMEM.  Otherwise, D2EDMP will
C             assume the correct value is still in IMEM(2).
C     MEMAXD  If positive, the actual size of DMEM.  Otherwise, D2EDMP will
C             assume the correct value is still in IMEM(3).
C     LIMIT   Maximum number of lines to print in the disk file for any one
C             data space of an entity.  If exceeded, prints the first LIMIT-2
C             and the last two lines, with a message in between.  If not
C             positive, LIMIT defaults to 20.
C     LUNIT   Logical unit number to use for the dump file.
C     IPTR    Pointer to the data entity to be dumped.
C
C   CALLS:    D1EDMP
C   HISTORY:
C     11Feb92 P Kraushar    Created.
C*****

C   Long Name Alias:
C     ENTRY D2_ENTITY_DUMP (CMEM, IMEM, DMEM, MEMAXC, MEMAXI, MEMAXD,
C    +    LIMIT, LUNIT, IPTR)

      CHARACTER CMEM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)
      INTEGER MEMAXC, MEMAXI, MEMAXD, LIMIT, LUNIT, IPTR

      INTEGER MXC, MXI, MXD, LIM, LOWHDR, I2, IHDR, MEMIN(3)

      DATA MEMIN / 2, 20, 2 /
C****
      IF (IPTR .EQ. 0) THEN
        WRITE (LUNIT, '(/A)') ' D2EDMP CALLED WITH NULL POINTER'
        RETURN
      END IF
      IHDR = IPTR/256
      I2 = MOD(ABS(IPTR),256)
      WRITE (LUNIT, 10) IPTR, IHDR, I2
   10 FORMAT(/' D2EDMP CALLED FOR POINTER',I13,' = ',I8,'#',I3.3)

C   Establish critical boundary parameters.
      IF (MEMAXC .GT. MEMIN(1)) THEN
        MXC = MEMAXC
        IF (MEMAXC .NE. IMEM(1)) WRITE (LUNIT, 20) MEMAXC
   20   FORMAT(' Given size for CMEM (',I12,') does not match IMEM(1)')
      ELSE IF (IMEM(1) .GT. MEMIN(1)) THEN
        MXC = IMEM(1)
      ELSE
        MXC = MEMIN(1)
        WRITE (LUNIT, '(A)') ' No valid bounds found for CMEM array'
      END IF

      IF (MEMAXI .GT. MEMIN(2)) THEN
        MXI = MEMAXI
        IF (MEMAXI .NE. IMEM(2)) WRITE (LUNIT, 30) MEMAXI
   30   FORMAT(' Given size for IMEM (',I12,') does not match IMEM(2)')
      ELSE IF (IMEM(2) .GT. MEMIN(2)) THEN
        MXI = IMEM(2)
      ELSE
        MXI = MEMIN(2)
        WRITE (LUNIT, '(A)') ' No valid bounds found for IMEM array'
      END IF

      IF (MEMAXD .GT. MEMIN(3)) THEN
        MXD = MEMAXD
        IF (MEMAXD .NE. IMEM(3)) WRITE (LUNIT, 40) MEMAXD
   40   FORMAT(' Given size for DMEM (',I12,') does not match IMEM(3)')
      ELSE IF (IMEM(3) .GT. MEMIN(3)) THEN
        MXD = IMEM(3)
      ELSE
        MXD = MEMIN(3)
        WRITE (LUNIT, '(A)') ' No valid bounds found for DMEM array'
      END IF

      LIM = LIMIT
      IF (LIMIT .LE. 0) LIM = 20

      IF (IMEM(16) .GT. MEMIN(2) .AND. IMEM(16) .LE. MXI) THEN
        LOWHDR = IMEM(16)
      ELSE IF (IMEM(5) .GT. MEMIN(2) .AND. IMEM(5) .LE. MXI) THEN
        LOWHDR = IMEM(5)
        WRITE (LUNIT, '(A,A)') ' Using IMEM(5) in place of IMEM(16) ',
     +      'as header low limit'
      ELSE
        LOWHDR = (MEMIN(2) + MXI)/2
        WRITE (LUNIT, '(A,A)') ' Cannot identify header low limit.  ',
     +      'Using midpoint of usable space.'
      END IF

C   Call D1EDMP to do the actual dumping.  (This part shared with D2MDMP.)
      CALL D1EDMP (CMEM, IMEM, DMEM, MXC, MXI, MXD, LIM, LOWHDR,
     +    LUNIT, IHDR)

      RETURN
      END
