      SUBROUTINE D2CQST (CMEM, IMEM, DMEM, ICQ, JSEQ, LENSTR, STR, IER)

C     Store a string into a Character Sequence (CQ) entity
C
C   DYNAMIC MEMORY INPUT-OUTPUT:
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     ICQ     Pointer to a Character String Sequence entity
C     JSEQ    Ordinal number of string to store to ICQ
C     LENSTR  Length of STR.  LENSTR = 0 to blank the field.
C     STR     String to store in JSEQth position in ICQ
C   OUTPUT:
C     IER     Returned error value.  Zero implies no errors.  Otherwise,
C             -1  = Dynamic memory is corrupt or uninitialized
C             -2  = ICQ is a null pointer
C             -3  = Garbage value found in ICQ
C             -4  = Ambiguous--garbage pointer or deleted entity
C             -5  = ICQ points to a deleted entity
C             -6  = ICQ does not point to a Character String Sequence
C                   entity
C             -8  = JSEQ is less than or equal to zero
C             -9  = LENSTR is less than zero
C             -10 = Entity ICQ is internally inconsistent
C             -11 = Character space is locked
C             -12 = Integer space is locked
C             -13 = Insufficient Character space for needed expansion
C             -14 = Insufficient Integer space for needed expansion
C
C   CALLS:
C     D0PTR    Pointer check
C     D1CQST   Character sequence store
C     DTERR    Error handler
C
C   HISTORY:
C     14Jun93 D. Parsons   Created
C*****
C
C   Long Name Alias:
C     ENTRY D2_CHAR_STRING_SEQ_STORE (CMEM, IMEM, DMEM, ICQ, JSEQ,
C    +      LENSTR, STR, IER)

      CHARACTER CMEM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)
      INTEGER ICQ, JSEQ, LENSTR, IER
      CHARACTER STR*(*)

      INTEGER ENTCQ
      PARAMETER (ENTCQ=249)

      INTEGER IHDR, NEED, ITYP

      CHARACTER*6 SUBNAM
      DATA SUBNAM /'D2CQST'/

C****
      IER = 0


C     Check pointer and type, then locate relevant data spaces
      CALL D0PTR (CMEM, IMEM, DMEM, ICQ, ENTCQ, 0, IHDR, ITYP, IER)
      IF (IER .LT. 0) GOTO 9900

      CALL D1CQST (CMEM, IMEM, DMEM, IHDR, JSEQ, LENSTR, STR,
     +                   NEED, IER)

C     Error reporting section

 9900 CONTINUE

      IF ((IER .EQ. -11) .OR. (IER .EQ. -12)) THEN
          CALL DTERR (2, SUBNAM, IER, NEED)
      ELSE IF (IER .NE. 0) THEN
          CALL DTERR (1, SUBNAM, IER, 0)
      ENDIF

      RETURN
      END

