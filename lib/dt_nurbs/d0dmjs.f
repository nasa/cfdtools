      SUBROUTINE D0DMJS (CMEM, IMEM, DMEM, IDE, LUNIT, IER)

C     Joined Surface Entity check


      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*)
      DOUBLE PRECISION  DMEM(*)

      INTEGER IDE, LUNIT, IER, IERX, IDMY, IHDR
      INTEGER LENC, LENI, LEND
      INTEGER ICLS, MAXTS, I1, I2
      INTEGER ITS, IPTS, IJS, IXTS
      CHARACTER CDMY

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

      CALL D0JSFP (CMEM, IMEM, DMEM, IDE, CDMY, ICLS, MAXTS, IHDR,
     +      IERX)

      IF ((ICLS .LT. -1) .OR. (ICLS .GT. 1)) THEN
         IF (LUNIT .GE. 0) THEN
            WRITE (LUNIT, 125) I1, I2, ICLS
            IER = -999
         ELSE
            IER = -4
            GOTO 9999
         ENDIF
      ENDIF

      IF (MAXTS .LT. 1) THEN
         IF (LUNIT .GE. 0) THEN
            WRITE (LUNIT, 130) I1, I2, MAXTS
            IER = -999
         ELSE
            IER = -5
            GOTO 9999
         ENDIF
      ENDIF

      IF (LENI .NE. MAXTS+2) THEN
         IF (LUNIT .GE. 0) THEN
            WRITE (LUNIT, 140) I1, I2, LENI, MAXTS+2
            IER = -999
         ELSE
            IER = -2
            GOTO 9999
         ENDIF
      ENDIF

C     Loop through Trimmed Surface pointers, checking for validity

      DO 1000 ITS = 1, MAXTS
         CALL D2FEEI (CMEM, IMEM, DMEM, IDE, 2+ITS, IPTS, IERX)
         IF (IPTS .NE. 0) THEN
            CALL D0TSFP (CMEM, IMEM, DMEM, IPTS, CDMY, IDMY, IJS, IXTS,
     +            IDMY, IDMY, IHDR, IERX)
            IF (IERX .NE. 0) THEN
               IF (LUNIT .GE. 0) THEN
                  WRITE (LUNIT, 150) I1, I2, ITS
                  IER = -999
               ELSE
                  IER = -6
                  GOTO 9999
               ENDIF
            ELSE
               IF (IJS .NE. IDE) THEN
                  IF (LUNIT .GE. 0) THEN
                     WRITE (LUNIT, 160) I1, I2, ITS
                     IER = -999
                  ELSE
                     IER = -7
                     GOTO 9999
                  ENDIF
               ELSE IF (IXTS .NE. ITS) THEN
                  IF (LUNIT .GE. 0) THEN
                     WRITE (LUNIT, 170) I1, I2, ITS, IXTS
                     IER = -999
                  ELSE
                     IER = -8
                     GOTO 9999
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
 1000 CONTINUE

 9999 CONTINUE

      RETURN

  100 FORMAT (' IDE: ',I3,'#',I3.3,' JS: LENC = ',I5,' < 0.')
  110 FORMAT (' IDE: ',I3,'#',I3.3,' JS: LENI = ',I5,' < 1.')
  120 FORMAT (' IDE: ',I3,'#',I3.3,' JS: LEND = ',I5,' <> 0.')
  125 FORMAT (' IDE: ',I3,'#',I3.3,' JS: ICLS = ',I5,' <> (-1,0,+1).')
  130 FORMAT (' IDE: ',I3,'#',I3.3,' JS: MAXTS = ',I5,' < 0.')
  140 FORMAT (' IDE: ',I3,'#',I3.3,' JS: LENI = ',I5,' < ',I5,'.')
  150 FORMAT (' IDE: ',I3,'#',I3.3,
     +      ' JS: Trimmed Surf. ptr at index = ',I5,
     +      ' does not point to a valid Trimmed Surf.')
  160 FORMAT (' IDE: ',I3,'#',I3.3,
     +      ' JS: Trimmed Surface entity at index = ',I5,
     +      ' does not point back to this entity.')
  170 FORMAT (' IDE: ',I3,'#',I3.3,
     +      ' JS: Trimmed Surface entity at index = ',I5,
     +      ' has IXTS = ',I5,'.')

      END
