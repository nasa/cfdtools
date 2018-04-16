      SUBROUTINE D0DMPJ (CMEM, IMEM, DMEM, IDE, LUNIT, IER)

C     Point on a Joined Surface Entity check


      EXTERNAL DTMCON
      DOUBLE PRECISION DTMCON

      CHARACTER         CMEM*(*)
      INTEGER           IMEM(*)
      DOUBLE PRECISION  DMEM(*)

      INTEGER IDE, LUNIT, IER, IERX, IDMY, IHDR
      INTEGER LENC, LENI, LEND
      INTEGER IEG, ITS, I1, I2, NWORK
      INTEGER IBFC, IBFS, ILP, ITSE
      CHARACTER CDMY
      DOUBLE PRECISION TOL, T, U, V, X, Y, Z, UE, VE, XE, YE, ZE
      DOUBLE PRECISION UERR, VERR, XERR, YERR, ZERR

      IER   = 0
      TOL   = 100*DTMCON(5)

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
      IF (LENI .NE. 2) THEN
         IF (LUNIT .GE. 0) THEN
            WRITE (LUNIT, 110) I1, I2, LENI
            IER = -999
         ELSE
            IER = -2
            GOTO 9999
         ENDIF
      ENDIF

      CALL D2FELD (CMEM, IMEM, DMEM, IDE, LEND, IERX)
      IF (LEND .NE. 6) THEN
         IF (LUNIT .GE. 0) THEN
            WRITE (LUNIT, 120) I1, I2, LEND
            IER = -999
         ELSE
            IER = -3
            GOTO 9999
         ENDIF
      ENDIF

      CALL D0PJFP (CMEM, IMEM, DMEM, IDE, CDMY, IEG, ITS, 
     +      T, U, V, X, Y, Z, IHDR, IERX)

      IF (IEG .NE. 0) THEN
         CALL D0EGFP (CMEM, IMEM, DMEM, IEG, CDMY, IBFC, ILP, IDMY,
     +         IDMY, IHDR, IERX)
         IF (IERX .NE. 0) THEN
            IF (LUNIT .GE. 0) THEN
               WRITE (LUNIT, 130) I1, I2, IEG
               IER = -999
            ELSE
               IER = -4
               GOTO 9999
            ENDIF
         ENDIF
         CALL D0BFFP (CMEM, IMEM, DMEM, IBFC, CDMY, IDMY, IHDR, IERX)
         IF (IERX .NE. 0) THEN
            IF (LUNIT .GE. 0) THEN
               WRITE (LUNIT, 150) I1, I2, IBFC
               IER = -999
            ELSE
               IER = -6
               GOTO 9999
            ENDIF
         ENDIF
         CALL D0LPFP (CMEM, IMEM, DMEM, ILP, CDMY, ITSE, IDMY,
     +         IDMY, IHDR, IERX)
         IF (IERX .NE. 0) THEN
            IF (LUNIT .GE. 0) THEN
               WRITE (LUNIT, 160) I1, I2, ILP
               IER = -999
            ELSE
               IER = -7
               GOTO 9999
            ENDIF
         ENDIF
         CALL D0TSFP (CMEM, IMEM, DMEM, ITSE, CDMY, IBFS, IDMY, 
     +         IDMY, IDMY, IDMY, IHDR, IERX)
         IF (IERX .NE. 0) THEN
            IF (LUNIT .GE. 0) THEN
               WRITE (LUNIT, 170) I1, I2, ITSE
               IER = -999
            ELSE
               IER = -8
               GOTO 9999
            ENDIF
         ENDIF
         CALL D0BFFP (CMEM, IMEM, DMEM, IBFS, CDMY, IDMY, IHDR, IERX)
         IF (IERX .NE. 0) THEN
            IF (LUNIT .GE. 0) THEN
               WRITE (LUNIT, 180) I1, I2, IBFS
               IER = -999
            ELSE
               IER = -9
               GOTO 9999
            ENDIF
         ENDIF
         IF (ITS .NE. 0) THEN
            IF (ITSE .NE. ITS) THEN
               IF (LUNIT .GE. 0) THEN
                  WRITE (LUNIT, 190) I1, I2, ITS, ITSE
                  IER = -999
               ELSE
                  IER = -10
                  GOTO 9999
               ENDIF
            ENDIF
         ELSE
            ITS = ITSE
         ENDIF

C        Find the (u,v) surface parameters

         CALL D0PJI1 (CMEM, IMEM, DMEM, IBFC, T, NWORK, UE, VE, IERX)
         IF (IERX .NE. 0) THEN
            IF       (IERX .EQ. -9) THEN
               IF (LUNIT .GE. 0) THEN
                  WRITE (LUNIT, 200) I1, I2
                  IER = -999
               ELSE
                  IER = -11
                  GOTO 9999
               ENDIF
            ELSE IF  (IERX .EQ. -10) THEN
               IF (LUNIT .GE. 0) THEN
                  WRITE (LUNIT, 210) I1, I2
                  IER = -999
               ELSE
                  IER = -12
                  GOTO 9999
               ENDIF
            ELSE IF  (IERX .EQ. -11) THEN
               IF (LUNIT .GE. 0) THEN
                  WRITE (LUNIT, 220) I1, I2
                  IER = -999
               ELSE
                  IER = -13
                  GOTO 9999
               ENDIF
            ELSE IF  (IERX .EQ. -12) THEN
               IF (LUNIT .GE. 0) THEN
                  WRITE (LUNIT, 230) I1, I2
                  IER = -999
               ELSE
                  IER = -14
                  GOTO 9999
               ENDIF
            ELSE IF  (IERX .EQ. -13) THEN
               IF (LUNIT .GE. 0) THEN
                  WRITE (LUNIT, 240) I1, I2
                  IER = -999
               ELSE
                  IER = -15
                  GOTO 9999
               ENDIF
            ELSE IF  (IERX .EQ. -14) THEN
               IF (LUNIT .GE. 0) THEN
                  WRITE (LUNIT, 250) I1, I2, T
                  IER = -999
               ELSE
                  IER = -16
                  GOTO 9999
               ENDIF
            ELSE IF  (IERX .EQ. -15) THEN
               IF (LUNIT .GE. 0) THEN
                  WRITE (LUNIT, 260) I1, I2, T
                  IER = -999
               ELSE
                  IER = -17
                  GOTO 9999
               ENDIF
            ELSE IF  (IERX .EQ. -16) THEN
               IF (LUNIT .GE. 0) THEN
                  WRITE (LUNIT, 270) I1, I2
                  IER = -999
               ELSE
                  IER = -18
                  GOTO 9999
               ENDIF
            ELSE IF  (IERX .EQ. -17) THEN
               IF (LUNIT .GE. 0) THEN
                  WRITE (LUNIT, 280) I1, I2
                  IER = -999
               ELSE
                  IER = -19
                  GOTO 9999
               ENDIF
            ELSE IF  (IERX .EQ. -18) THEN
               IF (LUNIT .GE. 0) THEN
                  WRITE (LUNIT, 290) I1, I2
                  IER = -999
               ELSE
                  IER = -20
                  GOTO 9999
               ENDIF
            ELSE 
               IF (LUNIT .GE. 0) THEN
                  WRITE (LUNIT, 150) I1, I2, IBFC
                  IER = -999
               ELSE
                  IER = -6
                  GOTO 9999
               ENDIF
            ENDIF
         ENDIF
         UERR = ABS((UE-U))
         VERR = ABS((VE-V))
         IF ( (UERR .GT. TOL) .OR. (VERR .GT. TOL) ) THEN
            IF (LUNIT .GE. 0) THEN
               WRITE (LUNIT, 390) I1, I2, T, UE, VE, U, V
               IER = -999
            ELSE
               IER = -30
               GOTO 9999
            ENDIF
         ENDIF
      ENDIF

      IF (ITS .NE. 0) THEN
         CALL D0TSFP (CMEM, IMEM, DMEM, ITS, CDMY, IBFS, IDMY, IDMY,
     +         IDMY, IDMY, IHDR, IERX)
         IF (IERX .NE. 0) THEN
            IF (LUNIT .GE. 0) THEN
               WRITE (LUNIT, 140) I1, I2, ITS
               IER = -999
            ELSE
               IER = -5
               GOTO 9999
            ENDIF
         ENDIF
         CALL D0BFFP (CMEM, IMEM, DMEM, IBFS, CDMY, IDMY, IHDR, IERX)
         IF (IERX .NE. 0) THEN
            IF (LUNIT .GE. 0) THEN
               WRITE (LUNIT, 180) I1, I2, IBFS
               IER = -999
            ELSE
               IER = -9
               GOTO 9999
            ENDIF
         ENDIF

C        Find the (x, y, z) 3-Space coordinates

         CALL D0PJI2 (CMEM, IMEM, DMEM, IBFS, U, V, NWORK, 
     +         XE, YE, ZE, IERX)
         IF (IERX .NE. 0) THEN
            IF       (IERX .EQ. -9) THEN
               IF (LUNIT .GE. 0) THEN
                  WRITE (LUNIT, 200) I1, I2
                  IER = -999
               ELSE
                  IER = -11
                  GOTO 9999
               ENDIF
            ELSE IF  (IERX .EQ. -19) THEN
               IF (LUNIT .GE. 0) THEN
                  WRITE (LUNIT, 300) I1, I2
                  IER = -999
               ELSE
                  IER = -21
                  GOTO 9999
               ENDIF
            ELSE IF  (IERX .EQ. -20) THEN
               IF (LUNIT .GE. 0) THEN
                  WRITE (LUNIT, 310) I1, I2
                  IER = -999
               ELSE
                  IER = -22
                  GOTO 9999
               ENDIF
            ELSE IF  (IERX .EQ. -21) THEN
               IF (LUNIT .GE. 0) THEN
                  WRITE (LUNIT, 320) I1, I2
                  IER = -999
               ELSE
                  IER = -23
                  GOTO 9999
               ENDIF
            ELSE IF  (IERX .EQ. -22) THEN
               IF (LUNIT .GE. 0) THEN
                  WRITE (LUNIT, 330) I1, I2
                  IER = -999
               ELSE
                  IER = -24
                  GOTO 9999
               ENDIF
            ELSE IF  (IERX .EQ. -23) THEN
               IF (LUNIT .GE. 0) THEN
                  WRITE (LUNIT, 340) I1, I2, U, V
                  IER = -999
               ELSE
                  IER = -25
                  GOTO 9999
               ENDIF
            ELSE IF  (IERX .EQ. -24) THEN
               IF (LUNIT .GE. 0) THEN
                  WRITE (LUNIT, 350) I1, I2, U, V
                  IER = -999
               ELSE
                  IER = -26
                  GOTO 9999
               ENDIF
            ELSE IF  (IERX .EQ. -25) THEN
               IF (LUNIT .GE. 0) THEN
                  WRITE (LUNIT, 360) I1, I2
                  IER = -999
               ELSE
                  IER = -27
                  GOTO 9999
               ENDIF
            ELSE IF  (IERX .EQ. -26) THEN
               IF (LUNIT .GE. 0) THEN
                  WRITE (LUNIT, 370) I1, I2
                  IER = -999
               ELSE
                  IER = -28
                  GOTO 9999
               ENDIF
            ELSE IF  (IERX .EQ. -27) THEN
               IF (LUNIT .GE. 0) THEN
                  WRITE (LUNIT, 380) I1, I2
                  IER = -999
               ELSE
                  IER = -29
                  GOTO 9999
               ENDIF
            ELSE 
               IF (LUNIT .GE. 0) THEN
                  WRITE (LUNIT, 180) I1, I2, IBFS
                  IER = -999
               ELSE
                  IER = -9
                  GOTO 9999
               ENDIF
            ENDIF
         ENDIF
         XERR = ABS((XE-X))
         YERR = ABS((YE-Y))
         ZERR = ABS((ZE-Z))
         IF ( (XERR .GT. TOL) 
     +         .OR. (YERR .GT. TOL) 
     +         .OR. (ZERR .GT. TOL) ) THEN
            IF (LUNIT .GE. 0) THEN
               WRITE (LUNIT, 400) I1, I2, U, V, XE, YE, ZE, X, Y, Z
               IER = -999
            ELSE
               IER = -31
               GOTO 9999
            ENDIF
         ENDIF
      ENDIF

 9999 CONTINUE

      RETURN

  100 FORMAT (' IDE: ',I3,'#',I3.3,' PJ: LENC = ',I5,' < 0.')
  110 FORMAT (' IDE: ',I3,'#',I3.3,' PJ: LENI = ',I5,' <> 2.')
  120 FORMAT (' IDE: ',I3,'#',I3.3,' PJ: LEND = ',I5,' <> 6.')
  130 FORMAT (' IDE: ',I3,'#',I3.3,
     +      ' PJ: IEG = ',I5,' Does not point to an Edge entity.')
  140 FORMAT (' IDE: ',I3,'#',I3.3,' PJ: ITS = ',I5,
     +      ' Does not point to a Trimmed Surface.')
  150 FORMAT (' IDE: ',I3,'#',I3.3,' PJ: IEG->IBFC = ',I5,' does not ',
     +      'point to a B-spline Function entity.')
  160 FORMAT (' IDE: ',I3,'#',I3.3,' PJ: IEG->ILP = ',I5,' does not ',
     +      'point to a Loop entity')
  170 FORMAT (' IDE: ',I3,'#',I3.3,' PJ: IEG->ILP->ITS = ',I5,' does ',
     +      'not point to a Trimmed Surface entity.')
  180 FORMAT (' IDE: ',I3,'#',I3.3,' PJ: IEG->ILP->ITS->IBFS = ',I5,
     +      ' does not point to a B-Spline function entity.')
  190 FORMAT (' IDE: ',I3,'#',I3.3,' PJ: ITS = ',I5,' is not zero and ',
     +      'IEG->ILP->ITS = ',I5,' are different.')
  200 FORMAT (' IDE: ',I3,'#',I3.3,' PJ: Insufficient DMEM for ',
     +      'working storage.')
  210 FORMAT (' IDE: ',I3,'#',I3.3,' PJ: The Curve B-spline function ',
     +      'entity has C(3) .LE. 0.')
  220 FORMAT (' IDE: ',I3,'#',I3.3,' PJ: The Curve B-spline function ',
     +      'entity has C(4) < C(3).')
  230 FORMAT (' IDE: ',I3,'#',I3.3,' PJ: The Curve B-spline function ',
     +      'entity has invalid knot set.')
  240 FORMAT (' IDE: ',I3,'#',I3.3,' PJ: The Curve B-spline function ',
     +      'entity has denominator = 0.')
  250 FORMAT (' IDE: ',I3,'#',I3.3,' PJ: T = ',F5.3,' is inside too ',
     +      'small an interval.')
  260 FORMAT (' IDE: ',I3,'#',I3.3,' PJ: T = ',F5.3,' out of range.')
  270 FORMAT (' IDE: ',I3,'#',I3.3,' PJ: The Curve B-spline function ',
     +      'entity has NDOM .NE. 1.')
  280 FORMAT (' IDE: ',I3,'#',I3.3,' PJ: The Curve B-spline function ',
     +      'entity has NRNG .LT. 1.')
  290 FORMAT (' IDE: ',I3,'#',I3.3,' PJ: The Curve B-spline function ',
     +      'entity does not have 2 dependent variables.')
  300 FORMAT (' IDE: ',I3,'#',I3.3,' PJ: The Surface B-spline ',
     +      'function entity has K(i) .LE. 0 for some i.')
  310 FORMAT (' IDE: ',I3,'#',I3.3,' PJ: The Surface B-spline ',
     +      'function entity has a number of B-spline coefficients ',
     +      'with respect to the ith independent variable is less ',
     +      'than Ki for some i.')
  320 FORMAT (' IDE: ',I3,'#',I3.3,' PJ: The Surface B-spline ',
     +      'function entity has invalid knot set.')
  330 FORMAT (' IDE: ',I3,'#',I3.3,' PJ: The Surface B-spline ',
     +      'function entity has denominator = 0.')
  340 FORMAT (' IDE: ',I3,'#',I3.3,' PJ: U,V = ',F5.3,',',F5.3,' is ',
     +      'inside too small an interval.')
  350 FORMAT (' IDE: ',I3,'#',I3.3,' PJ: U = ',F5.3,' or V = ',F5.3,
     +      ' is out of range.')
  360 FORMAT (' IDE: ',I3,'#',I3.3,' PJ: The Surface B-spline ',
     +      'function entity has NDOM .NE. 2.')
  370 FORMAT (' IDE: ',I3,'#',I3.3,' PJ: The Surface B-spline ',
     +      'function entity has NRNG .LT. 1.')
  380 FORMAT (' IDE: ',I3,'#',I3.3,' PJ: The Surface B-spline ',
     +      'function entity does not have 3 dependent variables.')
  390 FORMAT (' IDE: ',I3,'#',I3.3,' PJ: IBFC evaluated at T = (',
     +      F5.3,') = (',F5.3,',',F5.3,
     +      ') does not match given (U,V) = (',F5.3,',',F5.3,').')
  400 FORMAT (' IDE: ',I3,'#',I3.3,' PJ: IBFS evaluated at (U,V) = (',
     +      F5.3,',',F5.3,') = (',F5.3,',',F5.3,',',F5.3,
     +      ') does not match given (X,Y,Z) = (',
     +      F5.3,',',F5.3,',',F5.3,').')

      END
