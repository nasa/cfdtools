      SUBROUTINE DTSEPP (C, NPTS, WORK, NWORK, V, IER)
C
C     EXTRACT POINTS FROM A SPLINE AT EQUALLY SPACED PARAMETER VALUES
C     IE., EVALUATE THE SPLINE AT T, T+DT, T+2DT, ...
C
C     WHERE DT = (TK - T0) / (NPTS-1)
C
C
C     USAGE
C
C         DOUBLE PRECISION C(NC), WORK(NWORK), V(NPTS, NDIM)
C         CALL DTSEPP (C, NPTS, WORK, NWORK, V, IER)
C
C         WHERE
C             NC IS THE LENGTH OF THE C ARRAY.
C             NDIM IS THE NUMBER OF DEPENDENT VARIABLES.
C                 NDIM = -1-C(2)  FOR A RATIONAL SPLINE (IE., C(2) < 0)
C                 NDIM = C(2)     OTHERWISE
C
C
C     INPUT
C
C         C       SPLINE DEFINITION ARRAY.
C
C         NPTS    THE NUMBER OF POINTS TO EVALUATE.  NPTS .GE. 1.
C
C
C     WORKING STORAGE
C
C         WORK    WORK ARRAY OF LENGTH NWORK.
C
C         NWORK   LENGTH OF ARRAY WORK; NWORK .GE. 5*C(3)-2.
C
C
C     OUTPUT
C
C         V       THE REQUESTED POINTS
C
C         IER     SUCCESS/ERROR CODE.
C                 FOR IER < 0, DTSEPP HAS SET V(1,1) TO DTMCON(1).
C
C                 IER =  0    NO ERRORS DETECTED.
C                 IER =  1    A POLE WAS FOUND AT THE I-TH POINT.
C                             DTSEPP HAS SET V(I) = DTMCON(3).
C                 IER = -1    C(3) .LE. 0.
C                 IER = -2    NPTS < 1
C                 IER = -3    NWORK TOO SMALL.
C                 IER = -6    C(4) .LT. C(3).
C                 IER = -8    INVALID KNOT SET.
C                 IER =-38    ATTEMPT TO EVALUATE AT POINT THAT IS
C                             INSIDE AN INTERVAL THAT IS TOO SMALL.
C                 IER =-51    C(1) .NE. 1.
C                 IER =-52    C(2) .EQ. 0. .OR. ABS(C(2)) .GT. MAXD
C
C
C     WRITTEN BY
C         DEBORAH PARSONS
C         APRIL 27, 1989
C
C     MODIFIED JAN 12 1991 BY P.KRAUSHAR TO ADD NDIMV ARGUMENT AND CORRECT 
C         ARRAY ERROR
C
      PARAMETER (MAXD=20)
      EXTERNAL DTSPVL, DTMCON
      INTEGER NPTS, NWORK, IER
      DOUBLE PRECISION C(*), WORK(*), V(NPTS,*), T, TSTART, TEND
      DOUBLE PRECISION VTMP(MAXD)
      DOUBLE PRECISION DT, DTMCON
      INTEGER IP, MODE, NEED, NDIM, I, J, ISTAT
C
      CHARACTER*8 SUBNAM
      DATA SUBNAM /'DTSEPV  '/
C
C     ERROR CHECKING
C
      IF (NPTS .LE. 0) THEN
          IER = -2
          GOTO 9900
      ENDIF
      NDIM = INT(C(2))
      IF (ABS(NDIM) .GT. MAXD) THEN
          IER = -52
          GOTO 9900
      ENDIF
      IF (NDIM .LT. 0) NDIM = -1 - NDIM
C
C     (ALL OTHER ERROR CHECKING HANDLED BY DTSPVL)
C
C     DON'T ALLOW ERROR MESSAGES FROM BELOW THIS LEVEL
C
      CALL DTERPT (0)
C
      IER    = 0
      IP     = 5 + C(3) + C(4)
      TSTART = C(6)
      TEND   = C(IP)
      DT = (TEND-TSTART) / (NPTS-1)
C
      DO 100 I = 1, NPTS
          T = TSTART + (I-1) * DT
          CALL DTSPVL (T, C, WORK, NWORK, VTMP, ISTAT)
C
          IF (ISTAT .EQ. 0) THEN
              DO 5 K = 1, NDIM
                  V(I,K) = VTMP(K)
   5          CONTINUE
          ELSE IF (ISTAT .EQ. -10) THEN
              IER = 1
              DO 10 J = 1, NDIM
                  V(I,J) = DTMCON(3)
   10         CONTINUE
          ELSE
              IER = ISTAT
              GOTO 9000
          ENDIF
  100 CONTINUE
C
 9000 CONTINUE
C
C     RE-ENABLE ERROR MESSAGES
C
      CALL DTERPT (1)
C
 9900 CONTINUE
      IF (IER .NE. 0) THEN
          IF (IER .GT. 0) THEN
              MODE = 0
          ELSE IF (IER .EQ. -3) THEN
              MODE = 2
              NEED = 5 * C(3) - 2
              V(1,1) = DTMCON(1)
          ELSE
              MODE = 1
              V(1,1) = DTMCON(1)
          ENDIF
          CALL DTERR(MODE, SUBNAM, IER, NEED)
      ENDIF
      RETURN
      END
