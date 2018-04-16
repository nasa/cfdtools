      SUBROUTINE D0IGRP (CLINE, LUNIT, IPCURL, IPCURP, EOP, EOR, 
     +                   EORFLG, TYPE, IVAL, RVAL, CVAL, ILEN)
C
C     SUBROUTINE TO READ A PARAMETER FROM THE IGES FILE
C
C     INPUT:
C        CLINE    The current line from the IGES File.
C        LUNIT    The unit number from which to fetch the next line,
C                 if necessary.
C        IPCURL   The pointer to the current line in the IGES file.
C        IPCURP   The pointer to the current position on the line.
C        EOP      The 'end-of-parameter' delimiter.
C        EOR      The 'end-of-record' delimiter.
C
C     OUTPUT:
C        CLINE    The current line from the IGES File.
C        IPCURL   The pointer to the current line in the IGES file.
C        IPCURP   The pointer to the current position on the line.
C        EORFLG   Whether or not the end of the IGES record was reached.
C        TYPE     Type of data found:
C                 'I' = Integer (default if no 'H','D','E', or '.' found)
C                 'R' = Real
C                 'C' = Character
C        IVAL     If TYPE = 'I', the desired integer, otherwise 0.
C        RVAL     If TYPE = 'R', the desired character string, otherwise 0.0.
C        CVAL     If TYPE = 'C', the desired character string, otherwise ' '.
C        ILEN     If TYPE = 'C', the length of the character string.
C
C     HISTORY
C        8/10/93  D. Parsons
C                 Extracted from HSRIGS and converted to read an actual
C                 file rather than an Array of IGES file lines.
C        9/5/93   D. Parsons
C                 Fixed problem reading real numbers without decimal
C                 points (eg. 1D0 as opposed to 1.0D0).
C
      INTEGER        IPCURL, IPCURP, MAXLEN, LUNIT
      CHARACTER*80   CLINE
      CHARACTER*(*)  CVAL
      CHARACTER*1    EOP, EOR, TYPE
      CHARACTER*9    FMTSTR
      INTEGER        IVAL, ILEN
      DOUBLE PRECISION RVAL
      LOGICAL        EORFLG
C
      INTEGER        I, IP, IPD, IDLEN, IOCHK
C
      EORFLG = .FALSE.
      CVAL = ' '
      IVAL = 0
      RVAL = 0.0
      ILEN = 0
      TYPE = 'I'
      MAXLEN = LEN(CVAL)
C
C     GET THE LENGTH OF THE DATA AREA
C
      IF (CLINE(73:73) .EQ. 'P') THEN
         IDLEN = 64
      ELSE
         IDLEN = 72
      ENDIF
C
C     LOOK FOR 'H', INDICATING THE END OF THE COUNT AND THE BEGINNING
C     OF THE STRING, OR 'D', 'E', OR '.' INDICATING A REAL VALUE.
C
      IP = IPCURP
      IPD = 0
   10 CONTINUE
      IF ((CLINE(IP:IP) .EQ. EOP) .OR. (CLINE(IP:IP) .EQ. EOR)) GOTO 20
      IF (CLINE(IP:IP) .EQ. 'H') THEN
         TYPE = 'C'
         GOTO 20
      ELSE IF (CLINE(IP:IP) .EQ. '.') THEN
         TYPE = 'R'
         IPD  = IP
      ELSE IF ((IPD .EQ. 0) .AND. 
     +         (CLINE(IP:IP) .EQ. 'D') .OR. 
     +         (CLINE(IP:IP) .EQ. 'E')) THEN
         TYPE = 'R'
      ENDIF
      IP = IP + 1
      IF (IP .GT. IDLEN) THEN
C
C        SINCE THE HOLERITH LENGTH OR THE NUMBER CAN NOT BE 'WRAPPED', 
C        THE STRING OR NUMBER HAS NOT YET BEGUN.  GET THE NEXT LINE.
C
         IPCURL = IPCURL + 1
         IPCURP = 1
         IP     = 1
         READ (LUNIT, '(A80)', REC=IPCURL, IOSTAT=IOCHK) CLINE
         IF (IOCHK .NE. 0) THEN
            EORFLG = .TRUE.
            GOTO 999
         ENDIF
      ENDIF
      GOTO 10
   20 CONTINUE
C
      IF (TYPE .EQ. 'C') THEN
      
C        READ IN THE CHARACTER STRING
C
C        IP POINTS AT THE 'H'
C
C        USE AN INTERNAL WRITE TO SET UP THE LENGTH OF THE NUMBER AND AN
C        INTERNAL READ TO GET THE LENGTH OF THE STRING
C
         WRITE (FMTSTR, '(2H(I,I2,1H))') IP-IPCURP
         READ (CLINE(IPCURP:IP-1), FMTSTR) ILEN
         DO 100 I = 1, MIN(ILEN, MAXLEN)
            IP = IP + 1
            IF (IP .GT. IDLEN) THEN
C        
C              GET THE NEXT LINE (STRING COVERS MORE THAN ONE LINE)
C        
               IPCURL = IPCURL + 1
               READ (LUNIT, '(A80)', REC=IPCURL, IOSTAT=IOCHK) CLINE
               IF (IOCHK .NE. 0) THEN
                  EORFLG = .TRUE.
                  GOTO 999
               ENDIF
               IP = 1
            ENDIF
            CVAL(I:I) = CLINE(IP:IP)
  100    CONTINUE
C
C        RESET THE CURRENT POSITION POINTER
C
         IF (CLINE(IP+1:IP+1) .EQ. EOR) EORFLG = .TRUE.
         IPCURP = IP+2
      ELSE
C
C        READ IN THE NUMBER STRING
C
C        IP POINTS AT THE EOP OR EOR
C
C        USE AN INTERNAL WRITE TO SET UP THE LENGTH OF THE NUMBER AND AN
C        INTERNAL READ TO PICK OFF THE NUMBER
C
         IF (IPD .GT. 0) THEN
C
C           DECIMAL FOUND, READ IN A REAL NUMBER WITH 'G' FORMAT SPECIFIER
C
            WRITE (FMTSTR, '(2H(G,I2,1H.,I2,1H))') IP-IPCURP, IP-IPD
            READ (CLINE(IPCURP:IP-1), FMTSTR) RVAL
            TYPE = 'R'
            
         ELSE IF ((TYPE .EQ. 'R') .AND. (IPD .EQ. 0)) THEN
C
C           NO DECIMAL FOUND (EG, 1D1), READ IN A REAL NUMBER WITH 
C           'G' FORMAT SPECIFIER
C
            WRITE (FMTSTR, '(2H(G,I2,1H.,I2,1H))') IP-IPCURP, 0
            READ (CLINE(IPCURP:IP-1), FMTSTR) RVAL
            
         ELSE
C
C           NO DECIMAL PLACE, USE 'I' FORMAT SPECIFIER
C
            IF (IP .EQ. IPCURP) THEN
               IVAL = 0
            ELSE        
               WRITE (FMTSTR, '(2H(I,I2,1H))') IP-IPCURP
               READ (CLINE(IPCURP:IP-1), FMTSTR) IVAL
               TYPE = 'I' 
            ENDIF
         ENDIF
C
C        RESET THE CURRENT POSITION POINTER, SKIPPING THE DELIMITER AT
C        THE END OF THE STRING
C
         IF (CLINE(IP:IP) .EQ. EOR) EORFLG = .TRUE.
         IPCURP = IP+1
      ENDIF
      
      IF (IPCURP .GT. IDLEN) THEN
         IPCURL = IPCURL + 1
         IPCURP = 1
         READ (LUNIT, '(A80)', REC=IPCURL, IOSTAT=IOCHK) CLINE
         IF (IOCHK .NE. 0) THEN
            EORFLG = .TRUE.
         ENDIF
      ENDIF 
         
  999 CONTINUE
      RETURN
      END
