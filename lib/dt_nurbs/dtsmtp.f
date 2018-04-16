      SUBROUTINE DTSMTP (X, NDIM, NX, IF, NIF, ICYC, XS, IER)
C
C     SMOOTH A SET OF POINTS.  IF IT IS NECESSARY THAT SOME POINTS
C     REMAIN FIXED, THE INDICES OF THOSE POINTS CAN BE INCLUDED IN
C     THE ARRAY IF.
C
C     SMOOTHED POINTS ARE COPIED TO XS AFTER A DELAY, SO THAT NO
C     SMOOTHED POINTS FROM THE CURRENT CYCLE ARE USED IN THAT CYCLE.
C
C
C     USAGE:
C
C         DOUBLE PRECISION X(NX,NDIM), XS(NX,NDIM)
C         INTEGER IF(NIF)
C         CALL DTSMTP (X, NDIM, NX, IF, NIF, ICYC, XS, IER)
C
C
C     INPUT:
C
C         X       ARRAY OF POINTS TO BE SMOOTHED
C         NDIM    DIMENSION OF POINTS TO BE SMOOTHED
C         NX      NUMBER OF POINTS TO BE SMOOTHED
C         IF      ARRAY OF INDICES OF FROZEN POINTS
C         NIF     NUMBER OF FROZEN POINTS
C         ICYC    NUMBER OF SMOOTHING CYCLES (ICYC .GE. 1)
C
C
C     OUTPUT:
C
C         XS      ARRAY OF SMOOTHED POINTS
C         IER     SUCCESS/ERROR CODE.
C                 FOR IER < 0, RESULTS HAVE NOT BEEN
C                 COMPUTED AND DTSMTP HAS SET XS(1,1) TO THE
C                 CLOBBER CONSTANT (DTMCON(1)).
C
C                 IER =  0    NO ERRORS DETECTED
C                 IER = -1    NDIM .LT. 1 OR NDIM .GT. 3
C                 IER = -2    NX .LT. 1
C                 IER = -3    NIF .LT. 0 OR NIF .GT. NX
C                 IER = -4    IF(I) .LE. 0 OR IF(I) .GT. NX
C                             FOR SOME I: 0 < I <= NX
C
C
C     WRITTEN BY
C                 DEBORAH PARSONS
C                 APRIL 13, 1989
C
C     ALGORITHM TAKEN FROM APPSMP, WRITTEN BY
C                 LLOYD R. TRACY
C                 SEPTEMBER 18, 1984
C
      EXTERNAL DTMCON
      INTEGER NDIM, NX, NIF, ICYC, IER, IF(*)
      DOUBLE PRECISION X(NX,*), XS(NX,*)
      DOUBLE PRECISION DTMCON, XT(3,3)
      INTEGER I, J, K
C
      CHARACTER*8 SUBNAM
      DATA SUBNAM /'DTSMTP  '/
C
C     ERROR CHECKING
C
      IF ((NDIM .LT. 1) .OR. (NDIM .GT. 3)) THEN
          IER = -1
          GOTO 9900
      ENDIF
C
      IF (NX .LT. 1) THEN
          IER = -2
          GOTO 9900
      ENDIF
C
      IF ((NIF .LT. 0) .OR. (NIF .GT. NX)) THEN
          IER = -3
          GOTO 9900
      ENDIF
C
C     COPY ALL OF THE INPUT POINTS INTO THE OUTPUT POINTS
C
      DO 100 I = 1, NX
          DO 100 J = 1, NDIM
              XS(I,J) = X(I,J)
  100 CONTINUE
C
C     IF THERE ARE FEWER THAN 3 POINTS, NO SMOOTHING CAN BE DONE
C
      IF (NX .LT. 3) GOTO 9000
C
C     IF ALL POINTS ARE FROZEN, SAVE TIME AND JUST EXIT
C
      IF (NIF .EQ. NX) GOTO 9000
C
C     START SMOOTHING
C
      DO 5000 I = 1, ICYC
C
C         XT(1) = XS(1)
C
          DO 200 K = 1, NDIM
              XT(1,K) = XS(1,K)
  200     CONTINUE
C
C         USE 3-POINT SMOOTHER ON 2ND POINT
C         XT(2) = ( XS(1) + 2*XS(2) + XS(3) ) / 4
C
          DO 300 K = 1, NDIM
              XT(2,K) = (XS(1,K) + 2.0*XS(2,K) + XS(3,K)) / 4.0
  300     CONTINUE
C
C         USE 5-POINT SMOOTHER ON INTERIOR POINTS
C         XT(3) = ( -XS(J-2) + 4*XS(J-1) + 10*XS(J) + 4*XS(J+1)
C                   -XS(J+2) ) / 16
C
          DO 700 J = 3, NX-2
              DO 400 K = 1, NDIM
                  XT(3,K) = (-XS(J-2,K) + 4.0*XS(J-1,K) + 10.0*XS(J,K)
     *                        + 4.0*XS(J+1,K) - XS(J+2,K))/16.0
  400         CONTINUE
C
C             SAVE THE SMOOTHED POINT XS(J-2)
C
              DO 500 K = 1, NDIM
                  XS(J-2,K) = XT(1,K)
  500         CONTINUE
C
C             ROTATE THE TEMPORARY POINTS
C
              DO 600 K = 1, NDIM
                  XT(1,K) = XT(2,K)
                  XT(2,K) = XT(3,K)
  600         CONTINUE
  700     CONTINUE
C
C         SAVE THE SMOOTHED POINT XS(NX-3)
C
          DO 800 K = 1, NDIM
              XS(NX-3,K) = XT(1,K)
  800     CONTINUE
C
C         DO THE SMOOTHED POINT XS(NX-1) IN PLACE
C
          DO 900 K = 1, NDIM
              XS(NX-1,K) = (XS(NX-2,K) + 2.0*XS(NX-1,K)
     *                      + XS(NX,K)) / 4.0
  900     CONTINUE
C
C         COPY XT2 TO XS(N-2)
C
          DO 1000 K = 1, NDIM
              XS(NX-2,K) = XT(2,K)
 1000     CONTINUE
C
C         COPY THE LAST POINT
C
          DO 1100 K = 1, NDIM
              XS(NX,K) = X(NX,K)
 1100     CONTINUE
C
C         PUT ALL OF THE FROZEN POINTS BACK
C
          DO 1300 J = 1, NIF
              IF ((IF(J) .LE. 0) .OR. (IF(J) .GT. NX)) THEN
                  IER = -4
                  GOTO 9900
              ENDIF
              DO 1200 K = 1, NDIM
                  XS(IF(J),K) = X(IF(J),K)
 1200         CONTINUE
 1300     CONTINUE
 5000 CONTINUE
C
 9000 CONTINUE
      RETURN
 9900 CONTINUE
      CALL DTERR (1, SUBNAM, IER, 0)
      XS(1,1) = DTMCON(1)
      RETURN
      END
