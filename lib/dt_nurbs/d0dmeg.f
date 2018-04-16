      SUBROUTINE D0DMEG (CMEM, IMEM, DMEM, IDE, LUNIT, IER)

C     Edge Entity check


      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*)
      DOUBLE PRECISION  DMEM(*)

      INTEGER IDE, LUNIT, IER, IERX, IDMY, IHDR
      INTEGER LENC, LENI, LEND
      INTEGER IBFC, ILP, IXEG, JEG, JEGX, MAXEG, IELM, I1, I2
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
      IF (LENI .NE. 4) THEN
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

      CALL D0EGFP (CMEM, IMEM, DMEM, IDE, CDMY, IBFC, ILP, IXEG,
     +      JEG, IHDR, IERX)

      IF (IBFC .NE. 0) THEN
         CALL D0BFFP (CMEM, IMEM, DMEM, IBFC, CDMY, IDMY, IHDR, IERX)
         IF (IERX .NE. 0) THEN
            IF (LUNIT .GE. 0) THEN
               WRITE (LUNIT, 130) I1, I2
               IER = -999
            ELSE
               IER = -4
               GOTO 9999
            ENDIF
         ENDIF
      ENDIF

      IF (ILP .NE. 0) THEN
         CALL D0LPFP (CMEM, IMEM, DMEM, ILP, CDMY, IDMY, IDMY,
     +      MAXEG, IHDR, IERX)
         IF (IERX .NE. 0) THEN
            IF (LUNIT .GE. 0) THEN
               WRITE (LUNIT, 140) I1, I2
               IER = -999
            ELSE
               IER = -5
               GOTO 9999
            ENDIF
         ENDIF
         IF (IXEG .EQ. 0) THEN
            IF (LUNIT .GE. 0) THEN
               WRITE (LUNIT, 150) I1, I2
               IER = -999
            ELSE
               IER = -6
               GOTO 9999
            ENDIF
         ELSE IF ((IXEG .GT. MAXEG) .OR. (IXEG .LT. 0)) THEN
            IF (LUNIT .GE. 0) THEN
               WRITE (LUNIT, 160) I1, I2, IXEG, MAXEG
               IER = -999
            ELSE
               IER = -8
               GOTO 9999
            ENDIF
         ELSE
            CALL D2FEEI (CMEM, IMEM, DMEM, ILP, 3+IXEG, IELM, IERX)
            IF (IELM .NE. IDE) THEN
               IF (LUNIT .GE. 0) THEN
                  WRITE (LUNIT, 170) I1, I2, IXEG
                  IER = -999
               ELSE
                  IER = -9
                  GOTO 9999
               ENDIF
            ENDIF
         ENDIF
      ELSE
         IF (IXEG .NE. 0) THEN
            IF (LUNIT .GE. 0) THEN
               WRITE (LUNIT, 180) I1, I2, IXEG
               IER = -999
            ELSE
               IER = -7
               GOTO 9999
            ENDIF
         ENDIF
      ENDIF

      IF (JEG .NE. 0) THEN
         CALL D0EGFP (CMEM, IMEM, DMEM, JEG, CDMY, IDMY, IDMY, IDMY,
     +      JEGX, IHDR, IERX)
         IF (IERX .NE. 0) THEN
            IF (LUNIT .GE. 0) THEN
               WRITE (LUNIT, 190) I1, I2
               IER = -999
            ELSE
               IER = -10
               GOTO 9999
            ENDIF
         ENDIF
         IF (JEGX .NE. IDE) THEN
            IF (LUNIT .GE. 0) THEN
               WRITE (LUNIT, 200) I1, I2
               IER = -999
            ELSE
               IER = -11
               GOTO 9999
            ENDIF
         ENDIF
      ENDIF

 9999 CONTINUE

      RETURN

  100 FORMAT (' IDE: ',I3,'#',I3.3,' EG: LENC = ',I5,' < 0.')
  110 FORMAT (' IDE: ',I3,'#',I3.3,' EG: LENI = ',I5,' <> 4.')
  120 FORMAT (' IDE: ',I3,'#',I3.3,' EG: LEND = ',I5,' <> 0.')
  130 FORMAT (' IDE: ',I3,'#',I3.3,
     +      ' EG: IBFC does not point to a valid B-Spline ',
     +      'function entity.')
  140 FORMAT (' IDE: ',I3,'#',I3.3,
     +      ' EG: ILP does not point to a valid Loop entity.')
  150 FORMAT (' IDE: ',I3,'#',I3.3,
     +      ' EG: ILP points to a Loop entity, but IXEG = 0.')
  160 FORMAT (' IDE: ',I3,'#',I3.3,
     +      ' EG: IXEG = ',I5,' not between 0 and ',
     +      'ILP(MAXEG) = ',I5,'.')
  170 FORMAT (' IDE: ',I3,'#',I3.3,
     +      ' EG: ILP(',I5,') does not point back to this ',
     +      'Edge entity.')
  180 FORMAT (' IDE: ',I3,'#',I3.3,' EG: ILP = 0, but IXEG = ',I5,'.')
  190 FORMAT (' IDE: ',I3,'#',I3.3,
     +      ' EG: JEG does not point to a valid Edge entity.')
  200 FORMAT (' IDE: ',I3,'#',I3.3,
     +      ' EG: JEG entity does not point back to this Edge.')
      END
