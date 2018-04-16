      SUBROUTINE DTMP31 (MAXFUN, TOLABS, TOLREL, STEP, ICC, X, F,
     +                                            HOLD, NHOLD, IER)
C
C     PURPOSE:
C
C        DTMP31 MINIMIZES A DOUBLE PRECISION FUNCTION F(X) OF A
C        SINGLE VARIABLE X WITHOUT USING DERIVATIVE VALUES.
C
C     USAGE:
C
C        CALL DTMP31 (MAXFUN, TOLABS, TOLREL, STEP, ICC, X, F,
C    +                                            HOLD, NHOLD, IER)
C
C     INPUT:
C
C        MAXFUN LIMIT ON THE NUMBER OF EVALUATIONS OF F(X).
C               AFTER MAXFUN EVALUATIONS ARE MADE AND THE SPECI-
C               FIED ACCURACY STILL NOT ACHIEVED, DTMP31 RETURNS WITH
C               F SET TO THE LOWEST VALUE OF F(X) OBTAINED.
C        TOLABS ABSOLUTE ACCURACY REQUIRED OF THE MINIMUM X.
C        TOLREL RELATIVE ACCURACY REQUIRED OF THE MINIMUM X.
C               IF THE CURRENT POSITION IS X(M) AND THE NEXT PRE-
C               DICTED POSITION IS X(M+1), THE MINIMUM X = X(M) IS
C               ACCEPTED IF
C               DABS(X(M) - X(M+1)) .LT. TOLABS            OR
C               DABS(X(M) - X(M+1)) .LT. DABS(TOLREL*X(M+1)).
C               TO SELECT ONLY ONE OF THESE ACCURACY REQUIREMENTS,
C               SET THE PARAMETER CORRESPONDING TO THE UNWANTED
C               ONE TO ZERO.
C        STEP   A REASONABLE CHANGE TO MAKE TO X TO START OFF
C               THE ITERATIONS. IT SHOULD BE AN ESTIMATE OF THE
C               ERROR IN THE INITIAL ESTIMATE X(0) TO THE TRUE
C               MINIMUM POSITION. A GROSS OVERESTIMATE MAY RESULT
C               IN A MINIMUM WHICH IS NOT THE ONE CLOSEST TO X(0).
C               ANY BAD ESTIMATE WILL AFFECT EFFICIENCY BUT NOT
C               FINAL CONVERGENCE.
C
C     INPUT/OUTPUT:
C
C        ICC    ON INPUT, ICC DEFINES WHICH COMPUTATION TO PERFORM
C               AS FOLLOWS:
C               INITIALIZE DTMP31:
C               ICC = -2 INITIALIZE ONLY. F IS NOT USED BUT
C                        THE NAME MUST APPEAR IN THE CALL.
C               ICC = -1 INITIALIZE AND PERFORM THE FIRST ITERATION
C                        STEP. ALL PARAMETERS IN THE CALL ARE USED.
C               PERFORM AN ITERATION STEP:
C               ICC = 1  PERFORM THE NEXT STEP.
C               ON OUTPUT, ICC HAS ONE OF THE FOLLOWING TWO VALUES:
C               ICC = 1  AN ITERATIVE STEP HAS BEEN PERFORMED BUT A
C                        MIMIMUM HAS NOT YET BEEN REACHED. EVALUATE F(X)
C                        FOR THE VALUE OF X RETURNED AND CALL
C                        DTMP31 AGAIN WITH ICC = 1.
C               ICC = 0  STOP. CHECK THE VALUE OF IER FOR RETURN MODE.
C        X      ON INPUT, AN ESTIMATE X(0) OF THE MINIMUM POSITION.
C               ON OUTPUT, THE VALUE OF THE REVISED ESTIMATE OF THE
C               MINIMUM COMPUTED WHEN AN ITERATIVE STEP IS PERFORMED.
C               ON A FINAL RETURN (ICC = 0, IER .GE. 0), X IS SET TO THE
C               POINT CORRESPONDING TO THE FUNCTION VALUE RETURNED IN F.
C               WHEN IER = 0, THIS IS THE REQUIRED MINIMUM. FOR A USER
C               ERROR RETURN (ICC = 0, IER .LT. 0), X IS UNCHANGED FROM
C               ITS INPUT VALUE.
C        F      ON INPUT, HOLDS THE VALUE OF
C               THE FUNCTION F(X).
C               IF INITIALIZING WITH ICC = -1, OR ON A RETURN TO THE
C               CALLING PROGRAM WITH ICC = 1, F MUST BE SET TO
C               THE VALUE OF F(X).
C               ON A FINAL RETURN (ICC = 0, IER .GE. 0), F WILL CONTAIN
C               THE LOWEST VALUE OF F(X) OBTAINED SO FAR.
C               FOR A USER ERROR RETURN
C               (ICC = 0, IER .LT. 0), F IS NOT SET.
C
C     WORKING STORAGE:
C
C        HOLD   AN ARRAY IN WHICH DTMP31 SAVES INTERNAL VARIABLE
C               VALUES REQUIRED FOR THE NEXT CALL.
C        NHOLD  THE LENGTH OF THE  HOLD  ARRAY (NHOLD .GE. 16).
C
C     OUTPUT:
C
C        IER    SUCCESS/ERROR CODE. RESULTS ARE RETURNED IN X AND F
C               FOR IER .GE. 0.
C               IER =  0 MINIMUM FOUND TO REQUIRED ACCURACY.
C               IER =  1 ACCURACY CANNOT BE ACHIEVED.
C               IER =  2 MORE THAN MAXFUN FUNCTION VALUES ARE NEEDED.
C               IER = -1 ICC .LT. -2 .OR. ICC .GT. 1 .OR. ICC = 0
C               IER = -2 MAXFUN .LE. 0
C               IER = -3 TOLABS .LT. ZERO .OR. TOLREL .LT. ZERO
C               IER = -4 TOLABS = ZERO = TOLREL
C               IER = -5 STEP = ZERO
C               IER = -6 NHOLD .LT. 16
C
C     SUBROUTINES CALLED:
C
C        DTERR
C
C
C     INTRINSIC FUNCTIONS USED:
C
C         DABS INT DBLE
C
C     AUTHOR:
C
C        DTMP31 - ADAPTED FROM HARWELL VD01A - MIKE HEALY
C
C     DATE:
C
C         4/08/86
C
      INTEGER MAXFUN, ICC, NHOLD, IER
      INTEGER IINC, IS, KK, MC, MODE, NEED
      DOUBLE PRECISION TOLABS, TOLREL, STEP, X, F
      DOUBLE PRECISION HOLD(*)
      DOUBLE PRECISION
     *      AXPXA,FA,FB,FC,D,DA,DB,DC,DADCD,
     *      XINC,XINCXC,XP,Z
      DOUBLE PRECISION P5, ZERO
C
      PARAMETER (P5 = 0.5D0, ZERO = 0.0D0)
C
      CHARACTER*8 SUBNAM
C
      DATA SUBNAM/'DTMP31'/
C
      KK = ICC
      ICC = 0
      IER = 0
      MODE = 1
      NEED = 0
C
C*********  ERROR CHECKING FOLLOWS:
C
      IF (KK .GE. -2 .AND. KK .LE. 1 .AND. KK .NE. 0) GO TO 100
            IER = -1
            GO TO 200
C
100   IF (MAXFUN .GT. 0) GO TO 110
            IER = -2
            GO TO 200
C
110   IF (TOLABS .GE. ZERO .AND. TOLREL .GE. ZERO)
     *      GO TO 120
            IER = -3
            GO TO 200
C
120   IF (TOLABS + TOLREL .NE. ZERO) GO TO 130
            IER = -4
            GO TO 200
C
130   IF (STEP .NE. ZERO) GO TO 140
            IER = -5
            GO TO 200
C
140   IF (NHOLD .GE. 16) GO TO 150
      IER = -6
      MODE = 2
      NEED = 16
      GO TO 200
C
  150 CONTINUE
C
C*********  END OF ERROR CHECKING.
C
      IF (KK .LT. 0) KK = KK + 4
      ICC = KK
      GO TO (1,2,2),ICC
    2 IS=6-ICC
      ICC=1
      IINC=1
      XINC=STEP+STEP
      MC=IS-3
      IF (MC .LE. 0) GO TO 4
      GO TO 300
    3 MC=MC+1
      IF (MAXFUN .GE. MC) GO TO 300
      ICC = 0
      IER = 2
      MODE = 0
   43 X=DB
      F=FB
      IF (FB .LE. FC) GO TO 200
      X=DC
      F=FC
      GO TO 200
C
    1 CONTINUE
      IS     = INT (HOLD( 1))
      IINC   = INT (HOLD( 2))
      MC     = INT (HOLD( 3))
      AXPXA  =       HOLD( 4)
      FA     =       HOLD( 5)
      FB     =       HOLD( 6)
      FC     =       HOLD( 7)
      D      =       HOLD( 8)
      DA     =       HOLD( 9)
      DB     =       HOLD(10)
      DC     =       HOLD(11)
      DADCD  =       HOLD(12)
      XINCXC =       HOLD(13)
      XINC   =       HOLD(14)
      XP     =       HOLD(15)
      Z      =       HOLD(16)
      GO TO (5,6,7,8),IS
    8 IS=3
    4 DC=X
      FC=F
      X=X+STEP
      GO TO 3
    7 IF (FC .LT. F)  GO TO 9
      IF (FC .GT. F) GO TO 11
   10 X=X+XINC
      XINC=XINC+XINC
      GO TO 3
    9 DB=X
      FB=F
      XINC=-XINC
      GO TO 13
   11 DB=DC
      FB=FC
      DC=X
      FC=F
   13 X=DC+DC-DB
      IS=2
      GO TO 3
    6 DA=DB
      DB=DC
      FA=FB
      FB=FC
   32 DC=X
      FC=F
      GO TO 14
    5 IF (FB .LT. FC) GO TO 16
      IF (F .GE. FB) GO TO 32
      FA=FB
      DA=DB
   19 FB=F
      DB=X
      GO TO 14
   16 IF (FA .LE. FC) GO TO 21
      XINC=FA
      FA=FC
      FC=XINC
      XINC=DA
      DA=DC
      DC=XINC
   21 XINC=DC
      IF ((D-DB)*(D-DC) .LT. ZERO) GO TO 32
      IF (F .GE. FA) GO TO 24
      FC=FB
      DC=DB
      GO TO 19
   24 FA=F
      DA=X
   14 IF (FB .GT. FC) GO TO 29
      IINC=2
      XINC=DC
      IF (FB .EQ. FC)  GO TO 45
   29 D=(FA-FB)/(DA-DB)-(FA-FC)/(DA-DC)
      IF(D*(DB-DC) .LE. ZERO) GO TO 33
      D=P5*(DB+DC-(FB-FC)/D)
      IF (DABS(D-X) .LE. TOLABS) GO TO 31
      IF (DABS(D-X) .GT. DABS(D*TOLREL)) GO TO 36
   31 ICC = 0
      IER = 0
      GO TO 43
   36 IS=1
      X=D
      DADCD = (DA - DC)*(DC - D)
      IF (DADCD .LT. ZERO) GO TO 3
      IF (DADCD .EQ. ZERO) GO TO 26
      IS=2
      GO TO (39,40),IINC
   39 IF (DABS(XINC) .LT. DABS(DC-D)) GO TO 41
      IF (DABS(XINC) .GE. DABS(DC-D)) GO TO 3
   33 IS=2
      GO TO (41,42),IINC
   41 X=DC
      GO TO 10
   40 IF (DABS(XINC-X) .GT. DABS(X-DC)) GO TO 3
   42 X=P5*(XINC+DC)
      XINCXC = (XINC - X)*(X - DC)
      IF (XINCXC .LE. ZERO) GO TO 26
      IF (XINCXC .GT. ZERO) GO TO 3
   45 X=P5*(DB+DC)
      IF ((DB-X)*(X-DC) .GT. ZERO) GO TO 3
   26 ICC = 0
      IER = 1
      MODE = 3
      GO TO 43
C
C*********  ERROR HANDLING CODE FOLLOWS.
C
200   CONTINUE
      IF (IER .EQ. 0) GO TO 300
      CALL DTERR (MODE, SUBNAM, IER, NEED)
      IF (IER.LT.0) RETURN
C
  300 CONTINUE
      HOLD( 1) = DBLE (IS    )
      HOLD( 2) = DBLE (IINC  )
      HOLD( 3) = DBLE (MC    )
      HOLD( 4) =       AXPXA
      HOLD( 5) =       FA
      HOLD( 6) =       FB
      HOLD( 7) =       FC
      HOLD( 8) =       D
      HOLD( 9) =       DA
      HOLD(10) =       DB
      HOLD(11) =       DC
      HOLD(12) =       DADCD
      HOLD(13) =       XINCXC
      HOLD(14) =       XINC
      HOLD(15) =       XP
      HOLD(16) =       Z
C
      RETURN
      END
