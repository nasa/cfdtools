      SUBROUTINE D1CQDL (CMEM, IMEM, DMEM, IHDR, JSEQ, IER)

C   Delete a string from a character sequence entity.
C
C   DYNAMIC MEMORY INPUT-OUTPUT:
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     IHDR    Header index of a Character String Sequence entity
C     JSEQ    Ordinal number of string to fetch from ICQ
C   OUTPUT:
C     IER     Returned error value.  Zero implies no errors.  
C             -8  = JSEQ is less or equal zero
C             -10 = Entity ICQ is internally inconsistent
C
C   CALLS: 
C     D1CQFE   Fetch string from Character Sequence
C     D1CQST   Store string to Character Sequence
C     D1FEEI   Fetch integer
C     D1STEI   Store integer
C
C   HISTORY:
C     29Jun93 D. Parsons    Created.
C*****

      CHARACTER      CMEM*(*)
      INTEGER        IMEM(*)
      DOUBLE PRECISION DMEM(*)
      
      INTEGER        JSEQ, LENSTR, IER, IHDR
      INTEGER        IERX, ENDPOS, NEED
      CHARACTER*(1)  STR

C****
      IER = 0
      
      IF (JSEQ .LE. 0) THEN
         IER = -8
         GOTO 9900
      ENDIF
      
      CALL D1CQFE (CMEM, IMEM, DMEM, IHDR, JSEQ, STR, LENSTR, IERX)
      IF (IERX .LT. 0) THEN
         IER = -10
         GOTO 9900
      ENDIF
      
C     If string is already deleted, or uninitialized.  Do nothing.
      
      IF (LENSTR .GE. 0) THEN
      
C        Store a "NULL" to free any used space.      
         IF (LENSTR .GT. 0) THEN
            CALL D1CQST (CMEM, IMEM, DMEM, IHDR, JSEQ, 0, ' ', NEED, 
     +         IERX)
            IF (IERX .LT. 0) THEN
               IER = -10
               GOTO 9900
            ENDIF
         ENDIF
         
C        Get the ending position for this string
         CALL D1FEEI (CMEM, IMEM, DMEM, IHDR, JSEQ+1, ENDPOS, IERX)
         IF (IERX .LT. 0) THEN
            IER = -10
            GOTO 9900
         ENDIF
         
C        Mark this field deleted.         
         IF (ENDPOS .GE. 0) THEN
            ENDPOS = -(ENDPOS+1)
            CALL D1STEI (CMEM, IMEM, DMEM, ENDPOS, JSEQ+1, IHDR, IERX)
            IF (IERX .LT. 0) THEN
               IER = -10
               GOTO 9900
            ENDIF
         ENDIF
      ENDIF
         
      
 9900 CONTINUE
 
      RETURN
      END
