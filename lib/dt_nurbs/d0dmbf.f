      SUBROUTINE D0DMBF (CMEM, IMEM, DMEM, IDE, LUNIT, IER)

C     B-Spline Function Entity check


      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*)
      DOUBLE PRECISION  DMEM(*)

      INTEGER IDE, LUNIT, IER, IERX, IDMY
      INTEGER LENC, LENI, LEND, ISTRT, I1, I2

      IER  = 0
      IERX = 0
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
      IF (LENI .NE. 0) THEN
         IF (LUNIT .GE. 0) THEN
            WRITE (LUNIT, 110) I1, I2, LENI
            IER = -999
         ELSE
            IER = -2
            GOTO 9999
         ENDIF
      ENDIF

      CALL D2FELD (CMEM, IMEM, DMEM, IDE, LEND, IERX)
      IF (LEND .LT. 0) THEN
         IF (LUNIT .GE. 0) THEN
            WRITE (LUNIT, 120) I1, I2, LEND
            IER = -999
         ELSE
            IER = -3
            GOTO 9999
         ENDIF
      ENDIF

      CALL D2FEBD (CMEM, IMEM, DMEM, IDE, ISTRT, IDMY, IERX)
      CALL D2UNLD (CMEM, IMEM, DMEM, IDE, IERX)

      CALL DTSCHK (DMEM(ISTRT), IERX)
      IF (IERX .NE. 0) THEN
         IF (LUNIT .GE. 0) THEN
            WRITE (LUNIT, 130) I1, I2
            IER = -999
         ELSE
            IER = -4
            GOTO 9999
         ENDIF
      ENDIF


 9999 CONTINUE

      RETURN

  100 FORMAT (' IDE: ',I3,'#',I3.3,' BF: LENC = ',I5,' < 0.')
  110 FORMAT (' IDE: ',I3,'#',I3.3,' BF: LENI = ',I5,' <> 0.')
  120 FORMAT (' IDE: ',I3,'#',I3.3,' BF: LEND = ',I5,' < 0.')
  130 FORMAT (' IDE: ',I3,'#',I3.3,
     +      ' BF: Spline contained in DMEM is invalid.')

      END
