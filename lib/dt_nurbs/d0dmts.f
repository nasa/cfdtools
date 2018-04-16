      SUBROUTINE D0DMTS (CMEM, IMEM, DMEM, IDE, LUNIT, IER)

C     Trimmed Surface Entity check


      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*)
      DOUBLE PRECISION  DMEM(*)

      INTEGER IDE, LUNIT, IER, IERX, IDMY, IHDR
      INTEGER LENC, LENI, LEND
      INTEGER IBFS, IJS, IXTS, IORIEN, MAXLP
      INTEGER MAXTS, IELM, ILP, IPLP, ITS, IXLP, I1, I2
      CHARACTER CDMY
      DOUBLE PRECISION DELM

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
      IF (LENI .LT. 4) THEN
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

      CALL D0TSFP (CMEM, IMEM, DMEM, IDE, CDMY, IBFS, IJS, IXTS,
     +      IORIEN, MAXLP, IHDR, IERX)

      IF (MAXLP .LT. 1) THEN
         IF (LUNIT .GE. 0) THEN
            WRITE (LUNIT, 130) I1, I2, MAXLP
            IER = -999
         ELSE
            IER = -11
            GOTO 9999
         ENDIF
      ENDIF

      IF (LENI .NE. MAXLP+5) THEN
         IF (LUNIT .GE. 0) THEN
            WRITE (LUNIT, 140) I1, I2, LENI, MAXLP+4
            IER = -999
         ELSE
            IER = -2
            GOTO 9999
         ENDIF
      ENDIF

      IF (IBFS .NE. 0) THEN
         CALL D0BFFP (CMEM, IMEM, DMEM, IBFS, CDMY, LEND, IHDR, IERX)
         IF ((IERX .EQ. 0) .AND. (LEND .GE. 1)) THEN
            CALL D2FEED (CMEM, IMEM, DMEM, IBFS, 1, DELM, IERX)
         ENDIF
         IF ((LEND .LT. 1) .OR.
     +         (DELM .NE. 2.0D0)
     +         .OR. (IERX .NE. 0)) THEN
            IF (LUNIT .GE. 0) THEN
               WRITE (LUNIT, 150) I1, I2
               IER = -999
            ELSE
               IER = -4
               GOTO 9999
            ENDIF
         ENDIF
      ENDIF

      IF (IJS .NE. 0) THEN
         CALL D0JSFP (CMEM, IMEM, DMEM, IJS, CDMY, IDMY, MAXTS, IHDR,
     +         IERX)
         IF (IERX .NE. 0) THEN
            IF (LUNIT .GE. 0) THEN
               WRITE (LUNIT, 160) I1, I2
               IER = -999
            ELSE
               IER = -5
               GOTO 9999
            ENDIF
         ENDIF
         IF (IXTS .EQ. 0) THEN
            IF (LUNIT .GE. 0) THEN
               WRITE (LUNIT, 170) I1, I2
               IER = -999
            ELSE
               IER = -6
               GOTO 9999
            ENDIF
         ELSE IF ((IXTS .GT. MAXTS) .OR. (IXTS .LT. 0)) THEN
            IF (LUNIT .GE. 0) THEN
               WRITE (LUNIT, 180) I1, I2, IXTS, MAXTS
               IER = -999
            ELSE
               IER = -8
               GOTO 9999
            ENDIF
         ELSE
            CALL D2FEEI (CMEM, IMEM, DMEM, IJS, 2+IXTS, IELM, IERX)
            IF (IELM .NE. IDE) THEN
               IF (LUNIT .GE. 0) THEN
                  WRITE (LUNIT, 190) I1, I2, IXTS
                  IER = -999
               ELSE
                  IER = -9
                  GOTO 9999
               ENDIF
            ENDIF
         ENDIF
      ELSE
         IF (IXTS .NE. 0) THEN
            IF (LUNIT .GE. 0) THEN
               WRITE (LUNIT, 200) I1, I2, IXTS
               IER = -999
            ELSE
               IER = -7
               GOTO 9999
            ENDIF
         ENDIF
      ENDIF

      IF (IABS(IORIEN) .NE. 1) THEN
         IF (LUNIT .GE. 0) THEN
            WRITE (LUNIT, 210) I1, I2, IORIEN
            IER = -999
         ELSE
            IER = -7
            GOTO 9999
         ENDIF
      ENDIF

C     Loop through Loop pointers, checking for validity

      DO 1000 ILP = 1, MAXLP
         CALL D2FEEI (CMEM, IMEM, DMEM, IDE, 5+ILP, IPLP, IERX)
         IF (IPLP .NE. 0) THEN
            CALL D0LPFP (CMEM, IMEM, DMEM, IPLP, CDMY, ITS, IXLP,
     +            IDMY, IHDR, IERX)
            IF (IERX .NE. 0) THEN
               IF (LUNIT .GE. 0) THEN
                  WRITE (LUNIT, 220) I1, I2, ILP
                  IER = -999
               ELSE
                  IER = -12
                  GOTO 9999
               ENDIF
            ELSE
               IF (ITS .NE. IDE) THEN
                  IF (LUNIT .GE. 0) THEN
                     WRITE (LUNIT, 230) I1, I2, ILP
                     IER = -999
                  ELSE
                     IER = -13
                     GOTO 9999
                  ENDIF
               ELSE IF (IXLP .NE. ILP) THEN
                  IF (LUNIT .GE. 0) THEN
                     WRITE (LUNIT, 240) I1, I2, ILP, IXLP
                     IER = -999
                  ELSE
                     IER = -14
                     GOTO 9999
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
 1000 CONTINUE

 9999 CONTINUE

      RETURN

  100 FORMAT (' IDE: ',I3,'#',I3.3,' TS: LENC = ',I5,' < 0.')
  110 FORMAT (' IDE: ',I3,'#',I3.3,' TS: LENI = ',I5,' < 4.')
  120 FORMAT (' IDE: ',I3,'#',I3.3,' TS: LEND = ',I5,' <> 0.')
  130 FORMAT (' IDE: ',I3,'#',I3.3,' TS: MAXLP = ',I5,' < 0.')
  140 FORMAT (' IDE: ',I3,'#',I3.3,' TS: LENI = ',I5,' < ',I5,'.')
  150 FORMAT (' IDE: ',I3,'#',I3.3,
     +      ' TS: IBFS does not point to a valid ',
     +      'Surface B-Spline entity.')
  160 FORMAT (' IDE: ',I3,'#',I3.3,
     +      ' TS: IJS does not point to a valid ',
     +      'Joined Surface entity.')
  170 FORMAT (' IDE: ',I3,'#',I3.3,
     +      ' TS: IJS points to a Joined Surface entity, ',
     +      'but IXTS = 0.')
  180 FORMAT (' IDE: ',I3,'#',I3.3,
     +      ' TS: IXTS = ',I5,' not between 0 and ',
     +      'IJS(MAXTS) = ',I5,'.')
  190 FORMAT (' IDE: ',I3,'#',I3.3,
     +      ' TS: IJS(',I5,') does not point back to this entity.')
  200 FORMAT (' IDE: ',I3,'#',I3.3,
     +      ' TS: IJS = 0, but IXTS = ',I5,'.')
  210 FORMAT (' IDE: ',I3,'#',I3.3,' TS: IORIEN = ',I5,' <> +-1.')
  220 FORMAT (' IDE: ',I3,'#',I3.3,
     +      ' TS: Loop pointer at index = ',I5,
     +      ' does not point to a valid Loop.')
  230 FORMAT (' IDE: ',I3,'#',I3.3,
     +      ' TS: Loop entity at index = ',
     +      I5,' does not point back to this entity.')
  240 FORMAT (' IDE: ',I3,'#',I3.3,
     +      ' TS: Loop entity at index = ',
     +      I5,' has IXLP = ',I5,'.')

      END
