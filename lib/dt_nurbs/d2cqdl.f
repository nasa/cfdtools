      SUBROUTINE D2CQDL (CMEM, IMEM, DMEM, ICQ, JSEQ, IER)

C   Delete a string from a Character Sequence entity.
C
C   DYNAMIC MEMORY INPUT-OUTPUT:
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     ICQ     Pointer to a Character String Sequence entity
C     JSEQ    Ordinal number of string to delete from ICQ
C   OUTPUT:
C     IER     Returned error value.  Zero implies no errors.  
C             -1  = Dynamic memory is corrupt or uninitialized
C             -2  = ICQ is a null pointer
C             -3  = Garbage value found in ICQ
C             -4  = Ambiguous--garbage pointer or deleted entity
C             -5  = ICQ points to a deleted entity
C             -6  = ICQ does not point to a Character String Sequence
C                   entity
C             -8  = JSEQ is less or equal zero
C             -10 = Entity ICQ is internally inconsistent
C
C   CALLS:    D0PTR, DTERR, D1CQDL
C   HISTORY:
C     29Jun93 D. Parsons    Created.
C*****
C
C   Long Name Alias:
C     ENTRY D2_CHAR_STRING_SEQ_DELETE (CMEM, IMEM, DMEM, ICQ, JSEQ, IER)

      CHARACTER CMEM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)
      
      INTEGER ICQ, JSEQ, IER

      INTEGER ENTCQ
      PARAMETER (ENTCQ=249)

      INTEGER IHDR, ITYP
C****

      IER = 0
      
C     Check pointer and type, then locate relevant data spaces
      CALL D0PTR (CMEM, IMEM, DMEM, ICQ, ENTCQ, 0, IHDR, ITYP, IER)
      IF (IER .EQ. 0) THEN
         CALL D1CQDL (CMEM, IMEM, DMEM, IHDR, JSEQ, IER)
      ENDIF

C     Error reporting section

 9900 IF (IER .LT. 0) CALL DTERR (1, 'D2CQDL  ', IER, 0)
      
      RETURN
      END
