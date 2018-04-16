      SUBROUTINE DTSCHT (C, CHT, MINCD, MAXT, WORK, NWORK, T, NT, IER)
C
C     FIND POINTS ALONG A SPLINE, SATISFYING A CHORD-HEIGHT TOLERANCE.
C
C
C     USAGE
C
C         DOUBLE PRECISION C(NC), CHT, WORK(NWORK), T(MAXT)
C         CALL DTSCHT (C, CHT, MINCD, MAXT, WORK, NWORK, T, NT, IER)
C
C         WHERE NC IS THE LENGTH OF THE INPUT C ARRAY.
C
C
C     INPUT
C
C         C       A PARAMETRIC SPLINE IN B-SPLINE FORM.
C
C         CHT     THE DESIRED CHORD HEIGHT TOLERANCE.
C
C         MINCD   THE MINIMUM NUMBER OF CHORDS DESIRED.  MINCD .GE. 1.
C
C         MAXT    DIMENSION OF T.  THE MAXIMIMUM NUMBER OF PARAMETERS
C                 WHICH CAN BE RETURNED.
C
C
C     WORKING STORAGE
C
C         WORK    WORKING STORAGE OF DIMENSION NWORK.
C
C         NWORK   DIMENSION OF THE WORKING STORAGE ARRAY.
C                 NWORK = MAXT + 8
C
C
C     OUTPUT
C
C         T       THE SPLINE PARAMETERS.  TO GET THE COORDINATES OF
C                 THE POINTS, CALL DTSPVL.
C
C         NT      ABS(NT) = THE NUMBER OF PARAMETERS FOUND.  NT IS
C                 NEGATIVE IF AN ERROR OCCURRED DURING THE PROCESS.
C
C         IER     SUCCESS/ERROR CODE.
C                 FOR IER < 0, SOME PARAMETERS MAY HAVE BEEN COMPUTED.
C                 NT IS SET TO -NT
C
C                 IER = -1    C(1) .NE. 1.
C                 IER = -2    C(2) .GT. 6 OR C(2) .LT. -7.
C                 IER = -3    NWORK TOO SMALL.  NOT ENOUGH WORKING
C                             STORAGE
C                 IER = -4    MINCD .LE. 0.
C                 IER = -5    MAXT TOO SMALL.  NOT ENOUGH OUTPUT SPACE.
C                 IER = -6    C(3) .LE. 0.
C                 IER = -7    C(4) .LT. C(3).
C
C
C     HINT:   SINCE T ONLY RETURNS THE SPLINE PARAMETERS CORRESPONDING
C             TO THE REQUESTED POINTS, USE A LOOP SUCH AS THE FOLLOWING
C             TO GET THE COORDINATES OF THE POINTS:
C
C                   DO 10 I = 1, NT
C                       CALL DTSPVL(T(I), C, WORK, NWORK, V(I,1), IER)
C             10    CONTINUE
C
C
C     WRITTEN BY
C         DEBORAH PARSONS
C         APRIL 24, 1989
C
      EXTERNAL DTMCON
      INTEGER MINCD, NWORK, MAXT, IER, NDIM, MAXDIM
      PARAMETER (MAXDIM=6)
      DOUBLE PRECISION C(*), CHT, WORK(*), T(*), T0, T1, V
      DOUBLE PRECISION DTMCON, DMAX, TOL, DT
      INTEGER IP, ISTAT, NTWORK, ISTPTR, NT, MAXFUN, I
C
      CHARACTER*8 SUBNAM
      DATA SUBNAM /'DTSCHT  '/
C
C     ERROR CHECKING
C
      IER = 0
      NT  = 0
      IF (C(1) .NE. 1.0) THEN
          IER = -1
          GOTO 9000
      ENDIF
C
      NDIM = C(2)
      IF (NDIM .LT. 0) THEN
          NDIM = -1 - NDIM
      ENDIF
      IF ((NDIM .GT. MAXDIM) .OR. (NDIM .LT. 2)) THEN
          IER = -2
          GOTO 9000
      ENDIF
C
      IF (MINCD .LE. 0) THEN
          IER = -4
          GOTO 9000
      ENDIF
C
      IF (C(3) .LE. 0) THEN
          IER = -6
          GOTO 9000
      ENDIF
C
      IF (C(4) .LT. C(3)) THEN
          IER = -7
          GOTO 9000
      ENDIF
C
C     DIVIDE UP THE WORK SPACE
C     FIRST NTWORK ELEMENTS USED FOR WORK, REST USED FOR "STACK"
C
      NTWORK = 5 * C(3) - 2
      IF (NWORK .LT. NTWORK) THEN
          IER = -3
          GOTO 9000
      ENDIF
C
      ISTPTR = NTWORK
      MAXFUN = MINCD * 1000
      TOL    = MAX(SQRT(DTMCON(5))*CHT, DTMCON(5)**0.8)
C
C     SET T0 AND T1 TO THE STARTING AND ENDING PARAMETERS
C
      IP = 5 + C(3) + C(4)
      T0 = C(6)
      T1 = C(IP)
C
C     DIVIDE THE CURVE UP INTO MINCD SECTIONS
C
      IF ((ISTPTR+MINCD-1) .GT. NWORK) THEN
          IER = -3
          GOTO 9000
      ENDIF
C
      DT = (T1-T0) / MINCD
      DO 10 I = 1, MINCD-1
          WORK (ISTPTR+I) = T1
          T1 = T1-DT
   10 CONTINUE
      ISTPTR = ISTPTR+MINCD-1
C
C     ITERATE (PSEUDO-RECURSIVE) ACROSS THE SPLINE
C
 1000 CONTINUE
C
C         FIND THE MAXIMUM DISTANCE BETWEEN THE CURRENT SPLINE
C         SEGMENT AND THE CHORD
C
          CALL DTSCH1 (C, T0, T1, NDIM, MAXFUN, TOL, WORK, NTWORK,
     *                 DMAX, V, ISTAT)
          IF (ISTAT .EQ. 2) THEN
              IER = -4
              GOTO 9000
          ELSE IF (ISTAT .LT. 0) THEN
              IER = -99
              GOTO 9000
          ENDIF
          IF (DMAX .GT. CHT) THEN
C
C             PUSH PART OF THE SEGMENT ONTO A SAVE STACK
C
              ISTPTR = ISTPTR+1
              IF (ISTPTR .GT. NWORK) THEN
                  IER = -3
                  GOTO 9000
              ENDIF
              WORK(ISTPTR) = T1
              T1 = V
          ELSE
C
C             ONE SEGMENT SOLVED, POP THE NEXT SEGMENT OFF OF THE STACK
C
              NT = NT+1
              IF (NT .GT. MAXT-1) THEN
                  IER = -5
                  GOTO 9000
              ENDIF
              T(NT) = T0
              T0 = T1
              IF (ISTPTR .EQ. NTWORK) THEN
C
C                 NO MORE SEGMENTS LEFT--DONE.
C
                  NT = NT+1
                  T(NT) = T0
                  GOTO 9900
              ENDIF
              T1 = WORK(ISTPTR)
              ISTPTR = ISTPTR-1
          ENDIF
          GOTO 1000
 9000 CONTINUE
      NT = -NT
      CALL DTERR(1, SUBNAM, IER, 0)
 9900 CONTINUE
      RETURN
      END

      SUBROUTINE DTSCH1 (C, T0, T1, NDIM, MAXFUN, TOL, WORK, NTWORK,
     *                   DMAX, V, ISTAT)
C
C     FIND THE MAXIMUM DISTANCE BETWEEN THE SPLINE SEGMENT AND THE
C     CHORD.
C
C
C     LOCAL VARIABLES OF INTEREST:
C
C         EV      IS THE VALUE OF THE SPLINE AT PARAMETER V
C         ET0     IS THE VALUE OF THE SPLINE AT PARAMETER T0
C         ET1     IS THE VALUE OF THE SPLINE AT PARAMETER T1
C         C1MC0   IS THE VECTOR FROM C(T1) TO C(T0)
C         LENSQ   IS THE SQUARE OF THE LENGTH OF VECTOR C1MC0
C         C1MC    IS THE VECTOR FROM C(T1) TO C(V)
C         LNCMC1  IS THE LENGTH OF VECTOR C1MC
C         CMC0    IS THE VECTOR FROM C(V) TO C(T0)
C         LNCMC0  IS THE LENGTH OF VECTOR CMC0
C         P       IS THE POINT ON THE CHORD EVALUATED AT LAMBDA
C         PMC     IS THE VECTOR FROM P TO C(V)
C         LNPMC   IS THE LENGTH OF VECTOR PMC
C
C
C     WRITTEN BY
C         DEBORAH PARSONS
C         APRIL 19, 1989
C
      INTEGER NHOLD, MAXFUN, ISTAT, MAXDIM, NTWORK
      PARAMETER (NHOLD=16, MAXDIM=6)
      DOUBLE PRECISION HOLD(NHOLD), WORK(NTWORK), C(*), T0, T1, DMAX, V
      DOUBLE PRECISION TOLABS, TOLREL, STEP, LAMBDA, LENSQ, TOL
      DOUBLE PRECISION ET0(MAXDIM), ET1(MAXDIM), EV(MAXDIM)
      DOUBLE PRECISION C1MC0(MAXDIM), C1MC(MAXDIM), CMC0(MAXDIM)
      DOUBLE PRECISION P(MAXDIM), PMC(MAXDIM), LNCMC0, LNCMC1, LNPMC
      INTEGER NDIM, ICC, I
      EXTERNAL DTMCON, DDOT, DAXPY
      DOUBLE PRECISION DTMCON, DDOT, DAXPY
C
C     DON'T ALLOW ERROR MESSAGES FROM THIS LEVEL
C
      CALL DTERPT (0)
C
      TOLABS = TOL
      TOLREL = TOL
      STEP   = MAX ((T1-T0)/MAXFUN**2, SQRT(DTMCON(5)))
      ICC    = -2
C
C     EVALUATE C AT T0, AND T1
C
      CALL DTSPVL(T0, C, WORK, NTWORK, ET0, ISTAT)
      IF (ISTAT .NE. 0) THEN
          GOTO 9000
      ENDIF
      CALL DTSPVL(T1, C, WORK, NTWORK, ET1, ISTAT)
      IF (ISTAT .NE. 0) THEN
          GOTO 9000
      ENDIF
C
      DO 10 I = 1, NDIM
          C1MC0(I)  = ET1(I) - ET0(I)
   10 CONTINUE
C
      LENSQ  = DDOT (NDIM, C1MC0, 1, C1MC0, 1)
C
C     INITIALIZE V TO THE MIDPOINT OF THIS SEGMENT OF THE SPLINE C
C
      V      = (T0 + T1) / 2.0
C
C     ITERATION LOOP
C
 1000 CONTINUE
          CALL DTMP31 (MAXFUN, TOLABS, TOLREL, STEP, ICC, V, DMAX,
     *                 HOLD, NHOLD, ISTAT)
          IF (ISTAT .NE. 0) THEN
              GOTO 9000
          ENDIF
          IF (ICC .EQ. 1) THEN
C
C             EVALUATE THE DISTANCE TO THE CHORD AT THIS CURRENT
C             PARAMETER (V)
C
              IF ((V .GE. T0) .AND. (V .LE. T1)) THEN
                  CALL DTSPVL(V, C, WORK, NTWORK, EV, ISTAT)
                  IF (ISTAT .NE. 0) THEN
                      GOTO 9000
                  ENDIF
                  DO 1010 I = 1, NDIM
                      C1MC(I) = ET1(I) - EV(I)
                      CMC0(I) = EV(I) - ET0(I)
 1010             CONTINUE
                  LAMBDA = DDOT (NDIM, C1MC, 1, C1MC0, 1) / LENSQ
C
                  LNCMC0 = SQRT(DDOT(NDIM, CMC0, 1, CMC0, 1))
                  LNCMC1 = SQRT(DDOT(NDIM, C1MC, 1, C1MC, 1))
C
                  IF ((LAMBDA .GE. 0.0) .AND. (LAMBDA .LE. 1.0)) THEN
                      DO 1020 I = 1, NDIM
                          P(I) = LAMBDA * ET0(I) + (1.0-LAMBDA) * ET1(I)
                          PMC(I) = P(I) - EV(I)
 1020                 CONTINUE
                      LNPMC = SQRT(DDOT(NDIM, PMC, 1, PMC, 1))
                      DMAX  = - MIN (LNCMC0, LNPMC, LNCMC1)
                  ELSE
                      DMAX  = - MIN (LNCMC0, LNCMC1)
                  ENDIF
              ELSE
C
C                 V OUT OF THE INTERVAL, SET DMAX VERY HIGH
C
                  DMAX = SQRT(DTMCON(2))
              ENDIF
              GOTO 1000
          ENDIF
 9000 CONTINUE
C
C     CORRECT THE SIGN OF DMAX
C
      DMAX = -DMAX
      CALL DTERPT (1)
      RETURN
      END
