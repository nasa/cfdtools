      SUBROUTINE D0DMGI (CMEM, IMEM, DMEM, IDE, LUNIT, IER)

C     IGES Index Entity Check


      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*)
      DOUBLE PRECISION  DMEM(*)

      INTEGER IDE, LUNIT, IER, IERX, I1, I2
      INTEGER LENC, LENI, LEND
      INTEGER NSL, NDE, IGE, IHDR, ITYP, I

C                          IGES Entity Data
      INTEGER     ENTGE
      PARAMETER  (ENTGE = 238)

      IER   = 0
      IERX  = 0
      I1   = IDE/256
      I2   = MOD(ABS(IDE),256)

      CALL D2FELC (CMEM, IMEM, DMEM, IDE, LENC, IERX)
      IF (LENC .LT. 3) THEN
         IF (LUNIT .GE. 0) THEN
            WRITE (LUNIT, 100) I1, I2, LENC
            IER = -999
         ELSE
            IER = -1
            GOTO 9999
         ENDIF
      ENDIF

      CALL D2FELI (CMEM, IMEM, DMEM, IDE, LENI, IERX)
      IF (LENI .LT. 2) THEN
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

      CALL D2IGLI (CMEM, IMEM, DMEM, IDE, NSL, NDE, IERX)

      IF (NSL .LT. 0) THEN
         IF (LUNIT .GE. 0) THEN
            WRITE (LUNIT, 130) I1, I2, NSL
            IER = -999
         ELSE
            IER = -4
            GOTO 9999
         ENDIF
      ENDIF

      IF (NDE .LT. 0) THEN
         IF (LUNIT .GE. 0) THEN
            WRITE (LUNIT, 140) I1, I2, NDE
            IER = -999
         ELSE
            IER = -5
            GOTO 9999
         ENDIF
      ENDIF

      IF (LENC .NE. (NSL*72+3)) THEN
         IF (LUNIT .GE. 0) THEN
            WRITE (LUNIT, 150) I1, I2, LENC, (NSL*72+3)
            IER = -999
         ELSE
            IER = -6
            GOTO 9999
         ENDIF
      ENDIF

      IF (LENI .NE. (NDE/2+2)) THEN
         IF (LUNIT .GE. 0) THEN
            WRITE (LUNIT, 160) I1, I2, LENI, (NDE/2+2)
            IER = -999
         ELSE
            IER = -7
            GOTO 9999
         ENDIF
      ENDIF

C     Check IGE is a valid IGES Entity pointer.
      CALL D2FEEI (CMEM, IMEM, DMEM, IDE, 1, IGE, IERX)
      CALL D0PTR (CMEM, IMEM, DMEM, IGE, ENTGE, 0, IHDR,
     +    ITYP, IERX)
      IF (IERX .NE. 0) THEN
         IF (LUNIT .GE. 0) THEN
            WRITE (LUNIT, 170) I1, I2, IGE, IERX
            IER = -999
         ELSE
            IER = -8
            GOTO 9999
         ENDIF
      ENDIF


C     Check IGES Entity pointers for validity.
      DO 10 I = 1, NDE/2
         CALL D2FEEI (CMEM, IMEM, DMEM, IDE, 2+I, IGE, IERX)
         IF (IGE .NE. 0) THEN
            CALL D0PTR (CMEM, IMEM, DMEM, IGE, ENTGE, 0, IHDR,
     +         ITYP, IERX)
            IF (IERX .NE. 0) THEN
               IF (LUNIT .GE. 0) THEN
                  WRITE (LUNIT, 180) I1, I2, I, IGE, IERX
                  IER = -999
               ELSE
                  IER = -9
                  GOTO 9999
               ENDIF
            ENDIF
         ENDIF
   10 CONTINUE

 9999 CONTINUE

      RETURN

  100 FORMAT (' IDE: ',I3,'#',I3.3,' GI: LENC = ',I5,' < 3.')
  110 FORMAT (' IDE: ',I3,'#',I3.3,' GI: LENI = ',I5,' < 2.')
  120 FORMAT (' IDE: ',I3,'#',I3.3,' GI: LEND = ',I5,' != 0.')
  130 FORMAT (' IDE: ',I3,'#',I3.3,' GI: NSL = ',I5,' < 0.')
  140 FORMAT (' IDE: ',I3,'#',I3.3,' GI: NDE = ',I5,' < 0.')
  150 FORMAT (' IDE: ',I3,'#',I3.3,' GI: LENC = ',I5,' != NSL*72+3=',
     +          I5,'.')
  160 FORMAT (' IDE: ',I3,'#',I3.3,' GI: LENI = ',I5,' != NDE/2+2=',
     +          I5,'.')
  170 FORMAT (' IDE: ',I3,'#',I3.3,' GI: GLOBAL IGE = ',I5,
     +        'NOT A VALID IGES ENTITY; D0PTR RETURNED',I5,'.')
  180 FORMAT (' IDE: ',I3,'#',I3.3,' GI: DE # ',I5,'; IGE = ',I5,
     +        'NOT A VALID IGES ENTITY; D0PTR RETURNED',I5,'.')

      END
