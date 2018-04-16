      SUBROUTINE DTPLNE (PT, NORM, W, L, R, C, IER)
C
C     GENERATE A PLANE (RECTANGULAR), GIVEN A POINT, NORMAL, WIDTH,
C     LENGTH, AND ROTATION.
C
C
C     USAGE
C
C     DOUBLE PRECISION PT(3), NORM(3), W, L, R, C(MC)
C     CALL DTPLNE (PT, NORM, W, L, R, C, IER)
C
C     WHERE MC >= 28
C
C
C     INPUT
C
C     PT      THE POINT AT THE CENTER OF THE PLANE.
C
C     NORM    THE DESIRED NORMAL TO THE PLANE.
C
C     W       THE DESIRED WIDTH OF THE PLANE.
C
C     L       THE DESIRED LENGTH OF THE PLANE.
C
C     R       THE DESIRED X-Y ROTATION OF THE PLANE (COUNTERCLOCKWISE,
C             IN DEGREES).
C
C
C     OUTPUT
C
C     C       THE PLANE IN B-SPLINE FORM.
C
C     IER     SUCCESS/ERROR CODE.  IF IER < 0, THEN C(1) HAS BEEN
C             SET TO -1.
C
C             = 0     NO ERRORS.  RESULTS COMPUTED.
C
C             = -1    NORM(I) = 0 FOR ALL I.
C
C
C     WRITTEN BY
C         DEBORAH PARSONS
C         SEPTEMBER 12, 1989
C
      EXTERNAL DTERR, DTMCON
C
      DOUBLE PRECISION PT(*), NORM(*), W, L, R, C(*)
      DOUBLE PRECISION P(4,4), ROT(3), M(4,4), D2R, DTMCON
      DOUBLE PRECISION HALFL, HALFW, COSR, SINR, TEMP(4)
      INTEGER IER, I
C
      CHARACTER*8 SUBNAM
      DATA SUBNAM /'DTPLNE  '/
C
C     ERROR CHECKING
C
      IF (NORM(1)**2 + NORM(2)**2 + NORM(3)**2 .EQ. 0.0D0) THEN
          IER = -1
          GOTO 9900
      ENDIF
C
C     INITIALIZE DEGREES TO RADIANS CONVERSION
C
      D2R = DTMCON(15)
C
C     CREATE A RECTANGLE OF THE DESIRED DIMENSIONS AND ROTATION, 
C     CENTERED ABOUT THE ORIGIN OF THE X-Y PLANE
C
      HALFW = W * 0.5D0
      HALFL = L * 0.5D0
      COSR = COS(D2R*R)
      SINR = SIN(D2R*R)
C
      P(1,1) = -HALFW * COSR - HALFL * SINR
      P(2,1) = -HALFW * SINR + HALFL * COSR
      P(3,1) =  0.0D0
      P(4,1) =  1.0D0
C
      P(1,2) =  HALFW * COSR - HALFL * SINR
      P(2,2) =  HALFW * SINR + HALFL * COSR
      P(3,2) =  0.0D0
      P(4,2) =  1.0D0
C
      P(1,3) = -HALFW * COSR + HALFL * SINR
      P(2,3) = -HALFW * SINR - HALFL * COSR
      P(3,3) =  0.0D0
      P(4,3) =  1.0D0
C
      P(1,4) =  HALFW * COSR + HALFL * SINR
      P(2,4) =  HALFW * SINR - HALFL * COSR
      P(3,4) =  0.0D0
      P(4,4) =  1.0D0
C
C     ROTATE THE PLANE TO THE NORMAL AND TRANSLATE THE CENTER TO PT
C     STORING THE RESULTING COEFFICIENTS INTO THE C-ARRAY
C
      IF (NORM(1) .EQ. 0.0D0 .AND. NORM(3) .EQ. 0.0D0) THEN
          ROT(2) = 0.0D0
      ELSE
          ROT(2) = ATAN2(NORM(1),NORM(3)) / D2R
      ENDIF
      IF (ROT(2) .EQ. 0.0D0) THEN
          ROT(1) = ATAN2(-NORM(2),NORM(3)) / D2R
      ELSE
          IF (NORM(1) .EQ. 0.0D0 .AND. NORM(2) .EQ. 0.0D0) THEN
              ROT(1) = 0.0D0
          ELSE
              ROT(1) = ATAN2(-SIN(ROT(2))*NORM(2),NORM(1)) / D2R
          ENDIF
      ENDIF
      ROT(3) = 0.0D0
      CALL DTTROT (3, ROT, M, IER)
C
      DO 100 I = 1, 4
          TEMP(1) = PT(1)
          TEMP(2) = PT(2)
          TEMP(3) = PT(3)
          TEMP(4) = 0.0D0
          CALL DTMVPS (2, 4, 4, M, 4, P(1,I), TEMP, IER)
          C(16+I) = TEMP(1)
          C(20+I) = TEMP(2)
          C(24+I) = TEMP(3)
  100 CONTINUE
C
C     GENERATE THE REST OF THE C-ARRAY
C
      C(1)  = 2.0D0
      C(2)  = 3.0D0
      C(3)  = 2.0D0
      C(4)  = 2.0D0
      C(5)  = 2.0D0
      C(6)  = 2.0D0
      C(7)  = 0.0D0
      C(8)  = 0.0D0
      C(9)  = 0.0D0
      C(10) = 0.0D0
      C(11) = 1.0D0
      C(12) = 1.0D0
      C(13) = 0.0D0
      C(14) = 0.0D0
      C(15) = 1.0D0
      C(16) = 1.0D0
C
 9900 IF (IER .NE. 0) THEN
          CALL DTERR (1, SUBNAM, IER, 0)
          C(1) = -1.0D0
      ENDIF
C
 545  FORMAT (1X,4F10.5)
      RETURN
      END
