      SUBROUTINE D0IGWP (LINE, LUNIT, IPCURL, IPCURP, EOP, EOR, 
     +                   EORFLG, TYPE, SECT, SEQ, DSEQ, 
     +                   IVAL, RVAL, CVAL)
C
C     SUBROUTINE TO WRITE A PARAMETER TO THE IGES FILE
C
C     INPUT:
C        LINE     The line being filled for the IGES File.
C        LUNIT    The unit number from which to fetch the next line,
C                 if necessary.
C        IPCURL   The pointer to the current line in the IGES file.
C        IPCURP   The pointer to the current position on the line.
C        EOP      The 'end-of-parameter' delimiter.
C        EOR      The 'end-of-record' delimiter.
C        EORFLG   Whether or not this is the last parameter in this 
C                 record.
C        TYPE     Type of data found:
C                 'I' = Integer (default if no 'H' or '.' found)
C                 'R' = Real
C                 'C' = Character
C                 'O' = Omitted Value, just put an EOP or EOR.
C        SECT     The Section for column 73
C        SEQ      The Sequence number for columns (74:80)
C        DSEQ     If SECT = 'P', the DE pointer.
C        IVAL     If TYPE = 'I', the integer to write.
C        RVAL     If TYPE = 'R', the real value to write.
C        CVAL     If TYPE = 'C', the character string to write.
C
C     OUTPUT:
C        LINE     The line being filled for the IGES File.
C        IPCURL   The pointer to the current line in the IGES file.
C        IPCURP   The pointer to the current position on the line.
C        SEQ      If the line was written, the SEQ is incremented.
C
C     HISTORY
C        8/18/93  D. Parsons  Created
C
      CHARACTER*(*)  LINE, CVAL
      INTEGER        LUNIT, IPCURL, IPCURP, SEQ, DSEQ, IVAL
      CHARACTER*1    EOP, EOR, TYPE, SECT
      DOUBLE PRECISION RVAL
      LOGICAL        EORFLG

      INTEGER        MAXLEN
C
      INTEGER        ILEN
      CHARACTER*8    FMTSTR
      
C     ***

      IF (SECT .EQ. 'P') THEN
         MAXLEN = 65
      ELSE
         MAXLEN = 72
      ENDIF
      
      IF (TYPE .EQ. 'C') THEN
C        (I2,'H',A)
         ILEN = LEN(CVAL)+3
      ELSE IF (TYPE .EQ. 'I') THEN
C        (In)
         IF (IVAL .EQ. 0) THEN
            ILEN = 1
         ELSE
            ILEN = INT(LOG10(FLOAT(IVAL)))+1
         ENDIF
      ELSE IF (TYPE .EQ. 'R') THEN
C        (E15.8)
         ILEN = 15
      ELSE
C        Omitted
         ILEN = 0    
      ENDIF                         
      
C     Is there room on the line?

      IF (IPCURP+ILEN+1 .GT. MAXLEN) THEN
         IF (SECT .EQ. 'P') THEN
            WRITE (LINE(66:80),'(I7.7,''P'',I7.7)') DSEQ, SEQ
         ELSE
            WRITE (LINE(73:80),'(A1,I7.7)') SECT, SEQ
         ENDIF
         WRITE (LUNIT, '(A80)', REC=IPCURL) LINE(1:80)
         LINE   = ' '
         SEQ    = SEQ+1
         IPCURP = 1
         IPCURL = IPCURL+1
      ENDIF
      
      IF (TYPE .EQ. 'C') THEN
         WRITE (LINE(IPCURP:IPCURP+ILEN-1),'(I2,''H'',A)') 
     +         LEN(CVAL), CVAL
      ELSE IF (TYPE .EQ. 'I') THEN
         WRITE (FMTSTR, '(''(I'',I2,'')'')') ILEN
         WRITE (LINE(IPCURP:IPCURP+ILEN-1),FMTSTR) IVAL
      ELSE IF (TYPE .EQ. 'R') THEN
         WRITE (LINE(IPCURP:IPCURP+ILEN-1),'(E15.8)') RVAL
      ENDIF
      
      IF (EORFLG) THEN
         WRITE (LINE(IPCURP+ILEN:IPCURP+ILEN), '(A1)') EOR
         IF (SECT .EQ. 'P') THEN
            WRITE (LINE(66:80),'(I7.7,''P'',I7.7)') DSEQ, SEQ
         ELSE
            WRITE (LINE(73:80),'(A1,I7.7)') SECT, SEQ
         ENDIF
         WRITE (LUNIT, '(A80)', REC=IPCURL) LINE(1:80)
         LINE   = ' '
         SEQ    = SEQ+1
         IPCURP = 1
         IPCURL = IPCURL+1
      ELSE
         WRITE (LINE(IPCURP+ILEN:IPCURP+ILEN), '(A1)') EOP
         IPCURP = IPCURP+ILEN+1
      ENDIF
      
  999 CONTINUE
      RETURN
      END
