      SUBROUTINE D2MDMP (CMEM, IMEM, DMEM, MEMAXC, MEMAXI, MEMAXD,
     +    LIMIT, LUNIT, IER)

C   Memory Dump of the dynamic memory arrays to a disk file.  A debugging
C   aid which attempts to make minimal assumptions about the consistency of
C   the dynamic memory arrays, while still dumping their contents in an
C   organized format.
C
C   DYNAMIC MEMORY INPUT-OUTPUT:
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     MEMAXC  If positive, the actual size of CMEM.  Otherwise, D2MDMP will
C             assume the correct value is still in IMEM(1).
C     MEMAXI  If positive, the actual size of IMEM.  Otherwise, D2MDMP will
C             assume the correct value is still in IMEM(2).
C     MEMAXD  If positive, the actual size of DMEM.  Otherwise, D2MDMP will
C             assume the correct value is still in IMEM(3).
C     LIMIT   Maximum number of lines to print in the disk file for any one
C             data space of an entity.  If exceeded, prints the first LIMIT-2
C             and the last two lines, with a message in between.  If not
C             positive, LIMIT defaults to 20.
C     LUNIT   Logical unit number to use for the dump file.
C   INPUT-OUTPUT:
C     IER     On input, if IER is less than zero, it is assumed to be the
C             output from a call to D2MBAD; if IER is zero, D2MBAD will be
C             called and its output recorded; and if IER is positive, D2MBAD
C             will not be called and no message is recorded in the dump file.
C             On output, IER is unchanged unless D2MBAD is called, in which
C             case it has the value returned by D2MBAD.
C
C   CALLS:    D1EDMP, D2MBAD
C   HISTORY:
C     10Feb92 P Kraushar   Created.
C     25Mar92 D. Parsons   Removed FILNAM parameter for portability.
C                          Calling program now has to open unit LUNIT.
C*****
C
C   Long Name Alias:
C     ENTRY D2_MEMORY_DUMP (CMEM, IMEM, DMEM, MEMAXC, MEMAXI, MEMAXD,
C    +    LIMIT, LUNIT, IER)

      CHARACTER CMEM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)
      INTEGER MEMAXC, MEMAXI, MEMAXD, LIMIT, LUNIT, IER

      INTEGER MXC, MXI, MXD, LIM, LOWHDR, I, IHDR, MEMIN(3)
      LOGICAL ISBAD, D2MBAD

      DATA MEMIN / 2, 20, 2 /
C****
      WRITE (LUNIT, '(A)') ' DYNAMIC MEMORY DUMP '

      IF (IER .EQ. 0) THEN
        ISBAD = D2MBAD( CMEM, IMEM, DMEM, MEMAXC, MEMAXI, MEMAXD, IER)
        IF (ISBAD) THEN
          WRITE (LUNIT, '(A,I4)') ' D2MBAD reports IER =',IER
        ELSE
          WRITE (LUNIT, '(A)') ' D2MBAD finds no structural flaws.'
        END IF
      ELSE IF (IER .LT. 0) THEN
        WRITE (LUNIT, '(A,I4)') ' Assumed call to D2MBAD reports IER =',
     +      IER
      END IF

C   Establish critical boundary parameters.
      IF (MEMAXC .GT. MEMIN(1)) THEN
        MXC = MEMAXC
        IF (MEMAXC .NE. IMEM(1)) WRITE (LUNIT, 10) MEMAXC
   10   FORMAT(' Given size for CMEM (',I12,') does not match IMEM(1)')
      ELSE IF (IMEM(1) .GT. MEMIN(1)) THEN
        MXC = IMEM(1)
      ELSE
        MXC = MEMIN(1)
        WRITE (LUNIT, '(A)') ' No valid bounds found for CMEM array'
      END IF

      IF (MEMAXI .GT. MEMIN(2)) THEN
        MXI = MEMAXI
        IF (MEMAXI .NE. IMEM(2)) WRITE (LUNIT, 20) MEMAXI
   20   FORMAT(' Given size for IMEM (',I12,') does not match IMEM(2)')
      ELSE IF (IMEM(2) .GT. MEMIN(2)) THEN
        MXI = IMEM(2)
      ELSE
        MXI = MEMIN(2)
        WRITE (LUNIT, '(A)') ' No valid bounds found for IMEM array'
      END IF

      IF (MEMAXD .GT. MEMIN(3)) THEN
        MXD = MEMAXD
        IF (MEMAXD .NE. IMEM(3)) WRITE (LUNIT, 30) MEMAXD
   30   FORMAT(' Given size for DMEM (',I12,') does not match IMEM(3)')
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

C   Record areas at beginnings and ends of CMEM, IMEM, DMEM
      WRITE (LUNIT, 40) (IMEM(I), I=1,MEMIN(2)-1)
   40 FORMAT(/' Master Parameter Block at beginning of IMEM:'/
     +    25X,'     CMEM         IMEM         DMEM'/
     +    ' Total length            ',3I13/
     +    ' Next free location      ',3I13/
     +    ' Lowest unlocked location',3I13/
     +    ' Header of highest lock  ',3I13/
     +    ' Lowest lock or delete   ',3I13/
     +    ' Last used header            ',I13/
     +    ' Top of active header list   ',I13/
     +    ' Top of free header list     ',I13/
     +    ' Last sequence number new hdr',I13/)
      ISBAD = (IMEM(1) + IMEM(2) + IMEM(3) .EQ. IMEM(MXI))
      WRITE (LUNIT, 50) MXI, IMEM(MXI), ISBAD, CMEM(1:1), MXC,
     +    CMEM(MXC:MXC), DMEM(1), MXD, DMEM(MXD)
   50 FORMAT(' Last IMEM(',I8,') =',I13, L1/
     +    ' First CMEM = ',A1,'  Last CMEM(',I8,') = ',A1/
     +    ' First DMEM = ',G20.12,'  Last DMEM(',I8,') =',G20.12/)

C   Next dump all possible headers (and therefore all entities)
      DO 100 IHDR = MXI-8,LOWHDR,-8
        CALL D1EDMP (CMEM, IMEM, DMEM, MXC, MXI, MXD, LIM, LOWHDR,
     +      LUNIT, IHDR)
  100 CONTINUE

      RETURN
      END
