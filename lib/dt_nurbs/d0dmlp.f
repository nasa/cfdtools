      SUBROUTINE D0DMLP (CMEM, IMEM, DMEM, IDE, LUNIT, IER)

C     Loop Entity check


      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*)
      DOUBLE PRECISION  DMEM(*)

      INTEGER IDE, LUNIT, IER, IERX, IDMY, IHDR
      INTEGER LENC, LENI, LEND
      INTEGER ITS, IXLP, MAXEG
      INTEGER MAXLP, IELM, IEG, IPEG
      INTEGER ILP, IXEG, I1, I2
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
      IF (LENI .LT. 3) THEN
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

      CALL D0LPFP (CMEM, IMEM, DMEM, IDE, CDMY, ITS, IXLP, MAXEG, IHDR,
     +      IERX)

      IF (MAXEG .LT. 1) THEN
         IF (LUNIT .GE. 0) THEN
            WRITE (LUNIT, 130) I1, I2, MAXEG
            IER = -999
         ELSE
            IER = -9
            GOTO 9999
         ENDIF
      ENDIF

      IF (LENI .NE. MAXEG+3) THEN
         IF (LUNIT .GE. 0) THEN
            WRITE (LUNIT, 140) I1, I2, LENI, MAXEG+3
            IER = -999
         ELSE
            IER = -2
            GOTO 9999
         ENDIF
      ENDIF

      IF (ITS .NE. 0) THEN
         CALL D0TSFP (CMEM, IMEM, DMEM, ITS, CDMY, IDMY, IDMY, IDMY,
     +      IDMY, MAXLP, IHDR, IERX)
         IF (IERX .NE. 0) THEN
            IF (LUNIT .GE. 0) THEN
               WRITE (LUNIT, 150) I1, I2
               IER = -999
            ELSE
               IER = -4
               GOTO 9999
            ENDIF
         ENDIF
         IF (IXLP .EQ. 0) THEN
            IF (LUNIT .GE. 0) THEN
               WRITE (LUNIT, 160) I1, I2
               IER = -999
            ELSE
               IER = -5
               GOTO 9999
            ENDIF
         ELSE IF ((IXLP .GT. MAXLP) .OR. (IXLP .LT. 0)) THEN
            IF (LUNIT .GE. 0) THEN
               WRITE (LUNIT, 170) I1, I2, IXLP, MAXLP
               IER = -999
            ELSE
               IER = -7
               GOTO 9999
            ENDIF
         ELSE
            CALL D2FEEI (CMEM, IMEM, DMEM, ITS, 5+IXLP, IELM, IERX)
            IF (IELM .NE. IDE) THEN
               IF (LUNIT .GE. 0) THEN
                  WRITE (LUNIT, 180) I1, I2, IXLP
                  IER = -999
               ELSE
                  IER = -8
                  GOTO 9999
               ENDIF
            ENDIF
         ENDIF
      ELSE
         IF (IXLP .NE. 0) THEN
            IF (LUNIT .GE. 0) THEN
               WRITE (LUNIT, 190) I1, I2, IXLP
               IER = -999
            ELSE
               IER = -6
               GOTO 9999
            ENDIF
         ENDIF
      ENDIF

C     Loop through Edge pointers, checking for validity

      DO 1000 IEG = 1, MAXEG
         CALL D2FEEI (CMEM, IMEM, DMEM, IDE, 3+IEG, IPEG, IERX)
         IF (IPEG .NE. 0) THEN
            CALL D0EGFP (CMEM, IMEM, DMEM, IPEG, CDMY, IDMY, ILP,
     +            IXEG, IDMY, IHDR, IERX)
            IF (IERX .NE. 0) THEN
               IF (LUNIT .GE. 0) THEN
                  WRITE (LUNIT, 200) I1, I2, IEG
                  IER = -999
               ELSE
                  IER = -10
                  GOTO 9999
               ENDIF
            ELSE
               IF (ILP .NE. IDE) THEN
                  IF (LUNIT .GE. 0) THEN
                     WRITE (LUNIT, 210) I1, I2, IEG
                     IER = -999
                  ELSE
                     IER = -11
                     GOTO 9999
                  ENDIF
               ELSE IF (IXEG .NE. IEG) THEN
                  IF (LUNIT .GE. 0) THEN
                     WRITE (LUNIT, 220) I1, I2, IEG, IXEG
                     IER = -999
                  ELSE
                     IER = -12
                     GOTO 9999
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
 1000 CONTINUE

 9999 CONTINUE

      RETURN

  100 FORMAT (' IDE: ',I3,'#',I3.3,' LP: LENC = ',I5,' < 0.')
  110 FORMAT (' IDE: ',I3,'#',I3.3,' LP: LENI = ',I5,' < 3.')
  120 FORMAT (' IDE: ',I3,'#',I3.3,' LP: LEND = ',I5,' <> 0.')
  130 FORMAT (' IDE: ',I3,'#',I3.3,' LP: MAXEG = ',I5,' < 0.')
  140 FORMAT (' IDE: ',I3,'#',I3.3,' LP: LENI = ',I5,' < ',I5,'.')
  150 FORMAT (' IDE: ',I3,'#',I3.3,
     +    ' LP: ITS does not point to a valid Trimmed Surface entity.')
  160 FORMAT (' IDE: ',I3,'#',I3.3,
     +    ' LP: ITS points to a Trimmed Surface entity, ',
     +    'but IXLP = 0.')
  170 FORMAT (' IDE: ',I3,'#',I3.3,
     +    ' LP: IXLP = ',I5,' not between 0 and ITS(MAXLP) = ',I5,'.')
  180 FORMAT (' IDE: ',I3,'#',I3.3,
     +    ' LP: ITS(',I5,') does not point back to this Loop entity.')
  190 FORMAT (' IDE: ',I3,'#',I3.3,' LP: ITS = 0, but IXLP = ',I5,'.')
  200 FORMAT (' IDE: ',I3,'#',I3.3,' LP: Edge ptr at index = ',I5,
     +      ' does not point to a valid Edge entity.')
  210 FORMAT (' IDE: ',I3,'#',I3.3,
     +    ' LP: Edge entity at index = ',I5,
     +    ' does not point back to this Loop.')
  220 FORMAT (' IDE: ',I3,'#',I3.3,
     +    ' LP: Edge entity at index = ',I5,
     +    ' has IXEG = ',I5,'.')
      END
