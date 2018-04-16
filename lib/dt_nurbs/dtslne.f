      SUBROUTINE DTSLNE (NDIM, PS, PE, C, IER)
C
C     GENERATE A LINE IN SPLINE (C-ARRAY) FORM, GIVEN THE STARTING
C     AND ENDING POINTS.
C
C
C     USAGE
C
C         DOUBLE PRECISION PS(NDIM), PE(NDIM), C(NC)
C         CALL DTSLNE (NDIM, PS, PE, C, IER)
C
C         WHERE NC = 13 IF NDIM = 2
C            OR NC = 15 IF NDIM = 3
C
C
C     INPUT
C
C         NDIM    THE DIMENSION OF THE LINE.  NDIM = 2 (PLANAR) OR
C                 NDIM = 3 (3-SPACE).
C
C         PS      THE STARTING POINT OF THE LINE.
C
C         PE      THE ENDING POINT OF THE LINE.
C
C
C     OUTPUT
C
C         C       THE SPLINE ARRAY.
C
C         IER     SUCCESS/ERROR CODE.
C                 FOR IER < 0, DTSLNE HAS SET C(1) = -1.
C
C                 IER =  0    NO ERRORS DETECTED.
C
C                 IER = -1    NDIM < 2 OR NDIM > 3.
C
C
C     WRITTEN BY
C         DEBORAH PARSONS
C         AUGUST 14, 1989
C
      EXTERNAL DTERR
C
      INTEGER NDIM, IER
      DOUBLE PRECISION PS(*), PE(*), C(*)
C
      CHARACTER*8 SUBNAM
      DATA SUBNAM /'DTSLNE  '/
C
C     ERROR CHECKING
C
      IER = 0
      IF ((NDIM .LT. 2) .OR. (NDIM .GT. 3)) THEN
          IER = -1
          GOTO 9900
      ENDIF
C
C     LOAD THE C ARRAY
C
      C(1) = 1.0
      C(2) = NDIM
      C(3) = 2.0
      C(4) = 2.0
      C(5) = 2.0
      C(6) = 0.0
      C(7) = 0.0
      C(8) = 1.0
      C(9) = 1.0
      C(10) = PS(1)
      C(11) = PE(1)
      C(12) = PS(2)
      C(13) = PE(2)
      IF (NDIM .EQ. 3) THEN
          C(14) = PS(3)
          C(15) = PE(3)
      ENDIF
C
 9900 IF (IER .LT. 0) THEN
          C(1) = -1.0
          CALL DTERR (1, SUBNAM, IER, 0)
      ENDIF
      RETURN
      END
