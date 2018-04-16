      SUBROUTINE D0DMCQ (CMEM, IMEM, DMEM, IDE, LUNIT, IER)

C     Character Sequence Entity Check


      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*)
      DOUBLE PRECISION  DMEM(*)

      INTEGER IDE, LUNIT, IER, IERX, I1, I2
      INTEGER LENC, LENI, LEND
      INTEGER I, LPOS, EPOS, MXSEQ, IVAL

      IER   = 0
      IERX  = 0
      I1   = IDE/256
      I2   = MOD(ABS(IDE),256)

      CALL D2FELC (CMEM, IMEM, DMEM, IDE, LENC, IERX)
      IF (LENC .LT. 0) THEN
         IF (LUNIT .GE. 0) THEN
            WRITE (LUNIT, 100) I1, I2, LENC
            IER = -999
         ELSE
            IER = -1
            GOTO 9999
         ENDIF
      ENDIF

      CALL D2FELI (CMEM, IMEM, DMEM, IDE, LENI, IERX)
      IF (LENI .LT. 1) THEN
         IF (LUNIT .GE. 0) THEN
            WRITE (LUNIT, 110) I1, I2, LENI
            IER = -999
         ELSE
            IER = -2
            GOTO 9999
         ENDIF
      ENDIF

      CALL D2FELD (CMEM, IMEM, DMEM, IDE, LEND, IERX)
      IF (LEND .NE. 0) THEN
         IF (LUNIT .GE. 0) THEN
            WRITE (LUNIT, 120) I1, I2, LEND
            IER = -999
         ELSE
            IER = -3
            GOTO 9999
         ENDIF
      ENDIF

C     Check IMEM values

      LPOS = 0
      CALL D2FEEI (CMEM, IMEM, DMEM, IDE, 1, MXSEQ, IERX)
      DO 10 I = 1, MXSEQ
         CALL D2FEEI (CMEM, IMEM, DMEM, IDE, I+1, IVAL, IERX)
         IF (IVAL .LT. 0) THEN
            EPOS = -(IVAL+1)
         ELSE
            EPOS = IVAL
         ENDIF
         IF (EPOS .LT. LPOS) THEN
            IF (LUNIT .GE. 0) THEN
               WRITE (LUNIT, 130) I1, I2
               IER = -999
            ELSE
               IER = -4
               GOTO 9999
            ENDIF
         ENDIF
         LPOS = EPOS
   10 CONTINUE

      IF (LPOS .GT. LENC) THEN
         IF (LUNIT .GE. 0) THEN
            WRITE (LUNIT, 140) I1, I2, LPOS, LENC
            IER = -999
         ELSE
            IER = -5
            GOTO 9999
         ENDIF
      ENDIF



 9999 CONTINUE

      RETURN

  100 FORMAT (' IDE: ',I3,'#',I3.3,' CQ: LENC = ',I5,' < 0.')
  110 FORMAT (' IDE: ',I3,'#',I3.3,' CQ: LENI = ',I5,' < 1.')
  120 FORMAT (' IDE: ',I3,'#',I3.3,' CQ: LEND = ',I5,' != 0.')
  130 FORMAT (' IDE: ',I3,'#',I3.3,
     +        ' CQ: NON-INCREASING ENDING POSITIONS IN IMEM')
  140 FORMAT (' IDE: ',I3,'#',I3.3,
     +        ' CQ: LAST CMEM POSITION = ',I5,' < LENC = ',I5,'.')

      END
