      SUBROUTINE D2CQDF (CMEM, IMEM, DMEM, NSTR, NCHR, INSTR, ICQ, IER)

C   Define a character string sequence entity.  The space to allocate
C   may be given explicitly using NSTR and NCHR, or implicitly via the
C   initialization string INSTR.  INSTR contains a sequence of strings
C   marked off by delimiter characters.  Any character not present in
C   any of the strings may be used as a delimiter.  The subroutine
C   uses the first character in INSTR as the string delimiter.
C
C        Character String Sequence Entity:
C
C           CMEM  [1..len()]  Strings (with delimiters removed)
C
C           DMEM  (not used)
C
C           IMEM  [1]         Number of strings
C                 [2..1+No.]  Ending position of each string
C
C
C   DYNAMIC MEMORY INPUT-OUTPUT:
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     NSTR    Number of items in sequence.  The integer space allocated
C             is the maximum of this number and the number of strings
C             found in INSTR.
C     NCHR    Length of character space to allocate.  The character
C             space allocated is the maximum of this number and the
C             number of non-delimiter characters in INSTR.
C     INSTR   Initialization string.  First character is delimiter
C             used to separate successive strings.  Last character
C             should also be delimiter, but may be omitted.
C   OUTPUT:
C     ICQ     Pointer to character string sequence (CQ) entity
C     IER     Returned error value.  Zero implies no errors.  Otherwise,
C             -1  = Dynamic memory is corrupt or uninitialized
C             -2  = Insufficient space available for character data
C             -3  = Insufficient space available for integer data
C             -999= Unexpected error in called subroutine
C
C   CALLS:
C     D1CQDF   Define character sequence entity
C     DTERR    Error handler
C
C   HISTORY:
C     09Jun92 P Kraushar    Created.
C     28Jul92 D Parsons     Changed Length of string to trimmed length.
C     11Aug92 D Parsons     Change D2DEFE call to D1DEFE.  Add call
C                           to D1MBAD to test the dynamic memory.
C     25May93 D Parsons     Converted to call D1CQDF.
C
C*****
C
C   Long Name Alias:
C     ENTRY D2_CHAR_STRING_SEQ_DEFINE (CMEM, IMEM, DMEM, NSTR, NCHR,
C    +    INSTR, ICQ, IER)

      EXTERNAL D0TRMC, D1MBAD
      INTEGER  D0TRMC
      LOGICAL  D1MBAD

      CHARACTER CMEM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)

      INTEGER NSTR, NCHR, ICQ, IER
      CHARACTER INSTR*(*)

      INTEGER ENTCQ
      PARAMETER (ENTCQ=249)

      INTEGER IERX
C****
      IER = 0

C     Check Dynamic memory

      IF (CMEM(1:1) .NE. 'M' .AND. CMEM(1:1) .NE. 'L') THEN
         IF (D1MBAD (CMEM, IMEM, DMEM, 0, 0, 0, IERX)) THEN
            IER = -1
            GOTO 9900
         ENDIF
      ENDIF

C     Allocate the entity

      CALL D1CQDF (CMEM, IMEM, DMEM, NSTR, NCHR, INSTR, ICQ, IER)

C     Error reporting section

 9900 CONTINUE

      IF (IER .LT. 0) THEN
         IF (IER .GE. -3) THEN
            CALL DTERR (1, 'D2CQDF  ', IER, 0)
         ELSE
C           Unexpected error from called routine
            CALL DTERR (1, 'D1DEFE  ', IER, 0)
            IER = -999
            CALL DTERR (5, 'D2CQDF  ', IER, 0)
         ENDIF
      ENDIF

      RETURN
      END
