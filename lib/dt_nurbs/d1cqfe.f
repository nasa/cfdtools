      SUBROUTINE D1CQFE (CMEM, IMEM, DMEM, IHDR, JSEQ, STR, LENSTR, IER)

C   DYNAMIC MEMORY INPUT-OUTPUT:
C     CMEM    Dynamically managed character data
C     IMEM    Dynamically managed integer data
C     DMEM    Dynamically managed double precision data
C   INPUT:
C     IHDR    Header index of a Character String Sequence entity
C     JSEQ    Ordinal number of string to fetch from ICQ
C   OUTPUT:
C     STR     The JSEQth string in ICQ
C     LENSTR  Used length of STR.  If the string is deleted or 
C             uninitialized LENSTR = -1.
C     IER     Returned error value.  Zero implies no errors.  
C             If IER < 0, LENSTR = DTJCON(1).
C             >0  = String was truncated to fit into STR.  True length
C                   is given by IER.
C             -8  = JSEQ is less or equal zero
C             -10 = Entity ICQ is internally inconsistent
C
C   CALLS:
C     DTJCON   Integer machine constants
C     D0FITS   String copy
C
C   HISTORY:
C     23Jun92 P Kraushar    Created.
C     28Jun93 D. Parsons    Add support for deleted & undefined
C                           strings.  These will be designated with
C                           -(endpos+1) in their associated IMEM
C                           positions.
C     29Jun93 D. Parsons    Extracted guts from D2CQFE.
C*****

      CHARACTER CMEM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)
      
      EXTERNAL DTJCON
      INTEGER  DTJCON
      
      INTEGER JSEQ, LENSTR, IER
      CHARACTER STR*(*)

      INTEGER IHDR, NST, I, J, IXC0, IXI1
C****

      IER = 0
      STR = ' '
      
      IXC0 = IMEM(IHDR+1) - 1
      IXI1 = IMEM(IHDR+2)

C     Get current number of strings in sequence
      NST = IMEM(IXI1)

C     Consistency check
      IF (NST .GE. ABS (IMEM(IHDR+5))) THEN
        IER = -10
        GO TO 9900
      END IF

C     JSEQ range check
      IF (JSEQ .LE. 0) THEN
        IER = -8
        GO TO 9900
      ELSE IF (JSEQ .GT. NST) THEN
C       Uninitialized string
        IER = 0
        LENSTR = -1
        GO TO 9900
      END IF

C     Compute beginning and ending character offsets of selected string
      IF (JSEQ .EQ. 1) THEN
        I = 1
      ELSE
        I = IMEM(IXI1+JSEQ-1)
        IF (I .LT. 0) I = -(I+1)
        I = I+1
      END IF
      J = IMEM(IXI1+JSEQ)
      
C     Check for deleted string      
      IF (J .LT. 0) THEN
        IER = 0
        LENSTR = -1
        GO TO 9900
      END IF


C     Copy string into STR
      IF (J .GE. I) THEN
        CALL D0FITS (CMEM(IXC0+I:IXC0+J), STR, LENSTR, IER)
      ELSE
        STR = ' '
        LENSTR = 0
        IER = 0
      END IF

C   Error reporting section

 9900 CONTINUE
 
      IF (IER .LT. 0) THEN
         STR = ' '
         LENSTR = DTJCON(1)
      ENDIF
      
      RETURN
      END
