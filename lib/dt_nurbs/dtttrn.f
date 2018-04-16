      SUBROUTINE DTTTRN (NDIM, T, M, IER)
C
C     CREATE A TRANSFORMATION MATRIX TO TRANSLATE AN OBJECT.
C
C
C     USAGE
C
C         DOUBLE PRECISION T(NDIM), M(NDIM+1, NDIM+1)
C         CALL DTTTRN (NDIM, T, M, IER)
C
C
C     INPUT
C
C         NDIM    THE DIMENSION OF THE OBJECTS TO BE TRANSLATED.
C                 NDIM = 2 FOR PLANAR, OR NDIM = 3 FOR 3-D SPACE.
C
C         T       THE DISPLACEMENT AMOUNTS:  T(1) IS THE DISPLACEMENT
C                 ALONG THE X-AXIS, T(2) IS THE DISPLACEMENT ALONG THE
C                 Y-AXIS, AND (IF NDIM = 3) T(3) IS THE DISPLACEMENT
C                 ALONG THE Z-AXIS.
C
C
C     OUTPUT
C
C         M       THE TRANSLATION MATRIX IN HOMOGENEOUS COORDINATES.
C
C                         1   0   0   DX
C                 M =     0   1   0   DY
C                         0   0   1   DZ
C                         0   0   0   1
C
C                 TO TRANSLATE THE POINT (X,Y,Z), MULTIPLY
C
C                     1   0   0   DX          X
C                     0   1   0   DY    .     Y
C                     0   0   1   DZ          Z
C                     0   0   0   1           1
C
C                          X+DX
C                 TO GET   Y+DY
C                          Z+DZ
C                           1
C
C         IER     SUCCESS/ERROR CODE.
C                 FOR IER < 0, DTTTRN HAS SET M(1) TO DTMCON (1).
C
C                 IER =  0    NO ERRORS DETECTED
C
C                 IER = -1    NDIM < 2 OR NDIM > 3
C
C
C     WRITTEN BY
C         DEBORAH PARSONS
C         MAY 2, 1989
C
      EXTERNAL DTMCON
      INTEGER I, J, NDIM, IER
      DOUBLE PRECISION T(*), M(NDIM+1, *), DTMCON
C
      CHARACTER*8 SUBNAM
      DATA SUBNAM /'DTTTRN  '/
C
C     ERROR CHECKING
C
      IER = 0
      IF ((NDIM .LT. 2) .OR. (NDIM .GT. 3)) THEN
          IER = -1
          GOTO 9900
      ENDIF
C
C     INITIALIZE THE TRANSLATION MATRIX TO THE IDENTITY MATRIX
C
      DO 10 I = 1, NDIM + 1
          DO 10 J = 1, NDIM + 1
              IF (I .EQ. J) THEN
                  M(I,J) = 1.0
              ELSE
                  M(I,J) = 0.0
              ENDIF
   10 CONTINUE
C
C     FILL IN THE TRANSLATION PART
C
C             1   0   0   DX
C     M =     0   1   0   DY
C             0   0   1   DZ
C             0   0   0   1
C
      DO 20 I = 1, NDIM
          M(I,NDIM+1) = T(I)
   20 CONTINUE
      M(NDIM+1,NDIM+1) = 1.0
C
 9900 IF (IER .LT. 0) THEN
          M(1,1) = DTMCON(1)
          CALL DTERR (1, SUBNAM, IER, 0)
      ENDIF
      RETURN
      END
