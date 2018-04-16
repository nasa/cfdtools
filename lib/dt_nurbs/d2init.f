      SUBROUTINE D2INIT (CMEM, IMEM, DMEM, MEMAXC, MEMAXI, MEMAXD,
     +    TYPNAM, CHKLVL, IER)

C   Initialize dynamic memory arrays.  This routine must be called before
C   any other D2 subroutine.
C
C   DYNAMIC MEMORY INPUT-OUTPUT:
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     MEMAXC  Length of CMEM string
C     MEMAXI  Length of IMEM array
C     MEMAXD  Length of DMEM array
C     TYPNAM  String containing user-defined data type names for types
C             1, 2, 3, etc.  The first character is the delimiter which
C             is used to separate successive type names.  A null string
C             is indicated by two successive delimiters.
C     CHKLVL  Error check level to perform at.  Valid values are:
C             'H' = High: Perform all checks.
C             'M' = Medium: Perform user input checks, but skip internal
C                   dynamic memory consistency checks and computationally
C                   expensive checks.
C             'L' = Low: Skip all but essential checks.
C   OUTPUT:
C     IER     Error flag.  If negative, initialization was unsuccesful.
C             -1 = Inadequate CMEM space
C             -2 = Inadequate IMEM space
C             -3 = Inadequate DMEM space
C             -4 = Invalid value for CHKLVL
C             -5 = Error in TYPNAM or insufficient space to store it
C
C   CALLS:
C     D2CQDF
C     DTERR
C     DTERPT
C
C   HISTORY:
C     28Jan92 P Kraushar    Created.
C*****
C
C   Long Name Alias:
C     ENTRY D2_INITIALIZE (CMEM, IMEM, DMEM, MEMAXC, MEMAXI, MEMAXD,
C    +    TYPNAM, CHKLVL, IER)

      CHARACTER CMEM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)
      INTEGER MEMAXC, MEMAXI, MEMAXD, IER, I
      CHARACTER TYPNAM*(*), CHKLVL*1, SUBNAM*(6)

      INTEGER MEMIN(3)

      DATA MEMIN / 2, 20, 2/
      DATA SUBNAM /'D2INIT'/

C   Check absolute minimum requirements
      IF (MEMAXC .LT. MEMIN(1)) GO TO 9001
      IF (MEMAXI .LT. MEMIN(2)) GO TO 9002
      IF (MEMAXD .LT. MEMIN(3)) GO TO 9003
      IF (CHKLVL .NE. 'H' .AND. CHKLVL .NE. 'M' .AND. CHKLVL .NE. 'L')
     +    GO TO 9004

C   Initialize CMEM
      CMEM(1:1) = CHKLVL
      IF (CHKLVL .EQ. 'H') THEN
        DO 1 I=2,MEMAXC-1
    1     CMEM(I:I) = '*'
      END IF
      CMEM(MEMAXC:MEMAXC) = '~'

C   Initialize IMEM
C     Lengths of CMEM, IMEM, DMEM
      IMEM(1) = MEMAXC
      IMEM(2) = MEMAXI
      IMEM(3) = MEMAXD
C     Next free location in CMEM, IMEM, DMEM
      IMEM(4) = MEMIN(1)
      IMEM(5) = MEMIN(2)
      IMEM(6) = MEMIN(3)
C     Lowest unlocked location in CMEM, IMEM, DMEM
      IMEM(7) = IMEM(4)
      IMEM(8) = IMEM(5)
      IMEM(9) = IMEM(6)
C     Location of header holding highest lock in CMEM, IMEM, DMEM
      IMEM(10) = 0
      IMEM(11) = 0
      IMEM(12) = 0
C     Lowest location to check when compacting memory in CMEM, IMEM, DMEM
      IMEM(13) = IMEM(1)
      IMEM(14) = IMEM(2)
      IMEM(15) = IMEM(3)
C     Lowest used header location in IMEM
      IMEM(16) = IMEM(2)
C     Top of active header linked list in IMEM
      IMEM(17) = 0
C     Top of free header block linked list in IMEM
      IMEM(18) = 0
C     Last sequence number used in a new header block
      IMEM(19) = 0
      IF (CHKLVL .EQ. 'H') THEN
C       Set all unused to -1
        DO 2 I = 20, MEMAXI-1
    2     IMEM(I) = -1
      END IF
C     Checksum at end of IMEM
      IMEM(MEMAXI) = IMEM(1) + IMEM(2) + IMEM(3)

C   Initialize DMEM
      DMEM(1) = 0.987654321D10
      IF (CHKLVL .EQ. 'H') THEN
        DO 3 I = 2, MEMAXD-1
    3     DMEM(I) = -1.0D0
      END IF
      DMEM(MEMAXD) = DMEM(1)

C   Store user-defined type names for future reference
      CALL DTERPT(0)
      CALL D2CQDF (CMEM, IMEM, DMEM, 0, 0, TYPNAM, I, IER)
      CALL DTERPT(1)
      IF (IER .NE. 0) GO TO 9005

      IER = 0
      RETURN

C   Error reporting section

C     Inadequate CMEM space
 9001 IER = -1
      GO TO 9900

C     Inadequate IMEM space
 9002 IER = -2
      GO TO 9900

C     Inadequate DMEM space
 9003 IER = -3
      GO TO 9900

C     Bad error check level selection
 9004 IER = -4
      GO TO 9900

C     Problem with type names string or allocation
 9005 IER = -5
      GO TO 9900

 9900 CALL DTERR (1, SUBNAM, IER, 0)

      IF (MEMAXI .GE. 3) THEN
         IMEM(1) = 0
         IMEM(2) = 0
         IMEM(3) = 0
      ENDIF

      END
