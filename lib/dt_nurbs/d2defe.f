      SUBROUTINE D2DEFE (CMEM, IMEM, DMEM, ITYP, LENC, LENI, LEND,
     +    INITIZ, IDE, IER)

C   Define a new entity and return a pointer to it.  That is, allocate a
C   header block and any requested space in CMEM, IMEM and DMEM.
C
C   DYNAMIC MEMORY INPUT-OUTPUT:
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     ITYP    Data type id number for new entity
C     LENC    Character data length for new entity
C     LENI    Integer data length for new entity
C     LEND    Double precision length for new entity
C     INITIZ  Logical value indicating whether to initialize the data
C             spaces of the entity to ' ', 0, and 0.0D0, respectively.
C   OUTPUT:
C     IDE     Pointer to new data entity
C     IER     Returned error value.  Zero implies no errors.  Otherwise,
C             -1 = Dynamic memory is corrupt or uninitialized
C             -2 = Insufficient space available for character data
C             -3 = Insufficient space available for integer data
C             -4 = Insufficient space available for double precision data
C             -5 = Negative character data length requested
C             -6 = Negative integer data length requested
C             -7 = Negative double precision data length requested
C             -8 = Data type id value not between 1 and 255, inclusive.
C
C   CALLS:
C     D1DEFE
C     DTERR
C   HISTORY:
C     28Jan92 P Kraushar    Created.
C     17Mar92 D. Parsons    Lower level routine D2DEFE does not
C                           call DTERR
C     22Mar92 D. Parsons    Calls D1MBAD to check dynamic memory
C*****
C
C   Long Name Alias:
C     ENTRY D2_DEFINE_ENTITY (CMEM, IMEM, DMEM, ITYP, LENC, LENI,
C    +    LEND, INITIZ, IDE, IER)

      CHARACTER CMEM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)

      EXTERNAL D1MBAD
      LOGICAL  D1MBAD

      INTEGER ITYP, LENC, LENI, LEND, IDE, IER
      LOGICAL INITIZ

      CHARACTER*6 SUBNAM
      DATA SUBNAM /'D2DEFE'/

C****
      IER = 0

      IF (CMEM(1:1) .NE. 'M' .AND. CMEM(1:1) .NE. 'L') THEN
         IF (D1MBAD (CMEM, IMEM, DMEM, 0, 0, 0, IER)) 
     +         GOTO 9001
      ENDIF

      IF (LENC .LT. 0) GO TO 9005
      IF (LENI .LT. 0) GO TO 9006
      IF (LEND .LT. 0) GO TO 9007
      IF (ITYP .LT. 1 .OR. ITYP .GT. 255) GO TO 9008

      CALL D1DEFE (CMEM, IMEM, DMEM, ITYP, LENC, LENI, LEND, INITIZ,
     +    IDE, IER)
C         IER values of -2, -3, and -4 for "out of memory" errors are
C         detected by D1DEFE.
C
      IF (IER .NE. 0) GOTO 9900
C
      RETURN

C   Error reporting section

C     Corrupt dynamic memory
 9001 IER = -1
      GO TO 9900

C     Negative LENC input
 9005 IER = -5
      GO TO 9900

C     Negative LENI input
 9006 IER = -6
      GO TO 9900

C     Negative LEND input
 9007 IER = -7
      GO TO 9900

C     Bad data type id number
 9008 IER = -8
      GO TO 9900

 9900 CALL DTERR (1, SUBNAM, IER, 0)
      IDE = 0

      RETURN
      END
