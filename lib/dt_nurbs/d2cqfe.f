      SUBROUTINE D2CQFE (CMEM, IMEM, DMEM, ICQ, JSEQ, STR, LENSTR, IER)

C   Fetch a string from a Character Sequence Entity.
C
C   DYNAMIC MEMORY INPUT-OUTPUT:
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     ICQ     Pointer to a Character String Sequence entity
C     JSEQ    Ordinal number of string to fetch from ICQ
C   OUTPUT:
C     STR     The JSEQth string in ICQ
C     LENSTR  Used length of STR.  If the string is deleted or 
C             uninitialized LENSTR = -1.
C     IER     Returned error value.  Zero implies no errors.  
C             If IER < 0, LENSTR = DTJCON(1).
C             >0  = String was truncated to fit into STR.  True length
C                   is given by IER.
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
C   CALLS:    D0PTR, DTERR
C   HISTORY:
C     23Jun92 P Kraushar    Created.
C     28Jun93 D. Parsons    Add support for deleted & undefined
C                           strings.  These will be designated with
C                           -(endpos+1) in their associated IMEM
C                           positions.
C     29Jun93 D. Parsons    Extracted guts to D2CQFE.
C*****
C
C   Long Name Alias:
C     ENTRY D2_CHAR_STRING_SEQ_FETCH (CMEM, IMEM, DMEM, ICQ, JSEQ, STR,
C    +    LENSTR, IER)

      CHARACTER CMEM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)
      
      EXTERNAL DTJCON
      INTEGER  DTJCON
      
      INTEGER ICQ, JSEQ, LENSTR, IER
      CHARACTER STR*(*)

      INTEGER ENTCQ
      PARAMETER (ENTCQ=249)

      INTEGER IHDR, ITYP
C****

      IER = 0
      STR = ' '
      
C     Check pointer and type, then locate relevant data spaces
      CALL D0PTR (CMEM, IMEM, DMEM, ICQ, ENTCQ, 0, IHDR, ITYP, IER)
      IF (IER .EQ. 0) THEN
         CALL D1CQFE (CMEM, IMEM, DMEM, IHDR, JSEQ, STR, LENSTR, IER)
      ENDIF

C   Error reporting section

 9900 IF (IER .LT. 0) THEN
         CALL DTERR (1, 'D2CQFE  ', IER, 0)
         STR = ' '
         LENSTR = DTJCON(1)
      ENDIF
      
      RETURN
      END
