C<*>

      SUBROUTINE DTSPDR(X,NDER,C,WORK,NWORK,V,NDIMV,IER)
C
C  PURPOSE:  CALCULATE VALUES AND DERIVATIVES OF A UNIVARIATE SPLINE
C            GIVEN A PACKED C ARRAY.
C
C  INPUTS:
C           X       VALUE AT WHICH TO EVALUATE.
C
C           NDER    NUMBER OF DERIVATIVES TO EVALUATE,
C                   NDER .GE. 0.
C
C           C       SPLINE DEFINITION ARRAY.
C
C           NDIMV   FIRST DIMENSION OF V ARRAY,
C                   NDIMV .GE. ABS(C(2)).
C
C  WORKING STORAGE:
C
C           WORK    WORK ARRAY OF LENGTH NWORK.
C
C           NWORK   LENGTH OF WORK ARRAY,
C                   NWORK .GE. C(3)**2 + (3+NDER)*C(3)
C
C  OUTPUT:
C           V       ARRAY OF VALUES AND DERIVATIVES OF THE DEPENDENT
C                   VARIABLES.  J'TH DERIVATIVE OF I'TH VARIABLE =
C                             V(I,J+1).
C
C           IER     SUCCESS/ERROR CODE.
C                   IER=0     SUCCESS.
C                   IER=-1    C(3) .LE. 0.
C                   IER=-3    NWORK TOO SMALL.
C                   IER=-4    NDIMV .LT. C(2).
C                   IER=-6    C(4) .LT. C(3).
C                   IER=-8    KNOTS NOT ORDERED OR
C                             MULTIPLICITY TOO LARGE.
C                   IER=-10   DENOMINATOR .EQ. 0 ON A RATIONAL SPLINE
C                   IER=-38   ATTEMPT TO EVALUATE AT POINT THAT IS
C                             INSIDE AN INTERVAL THAT IS TOO SMALL.
C                   IER=-50   X OUT OF RANGE.
C                   IER=-51   C(1) .NE. 1.
C                   IER=-52   C(2) .EQ. 0. OR -1.
C                   IER=-53   NDER .LT. 0.
C
C  SUBROUTINES CALLED:  DTSPD1, DTMCON, DTERR
C
C  AUTHOR:  A.K. JONES
C
C  DATE:  1-NOV-84
C
C  MINOR REVISIONS BY FRITZ KLEIN, 21-JAN-85
C
C  MODIFIED BY DEBORAH PARSONS, MARCH 30, 1989
C     TO INCLUDE RATIONAL B-SPLINES
C
C  MODIFIED BY J. W. MANKE, 9/3/91, FOR RATIONAL B-SPLINE CERTIFICATION
C
      DOUBLE PRECISION X,C(*),DTMCON,V(NDIMV,*),WORK(*)
      DOUBLE PRECISION EPS, SUM, DIV
      CHARACTER*8 SUBNAM
      LOGICAL RATNL
      EXTERNAL DTSPD1,DTMCON,DTERR
C
      DATA SUBNAM /'DTSPDR  '/
C
      EPS = DTMCON(6)
C
C  BREAK UP INPUT C ARRAY.
C
      NDOM  = C(1)
      NRNG  = C(2)
      KORD  = C(3)
      NCOEF = C(4)
      ISPAN = C(5)
C
C  CHECK USER INPUT ERRORS
C
      IF (NRNG .LT. -1) THEN
          RATNL = .TRUE.
          NRAT  = - NRNG
          NRNG  = NRAT - 1
      ELSE IF (NRNG .GT. 0) THEN
          RATNL = .FALSE.
          NRAT  = 0
      ELSE
          IER = -52
          GOTO 9010
      ENDIF
C
      IF(KORD .GT. 0) GO TO 10
        IER = -1
        GO TO 9010
   10 CONTINUE
C
      IF(NDER .GE. 0) GO TO 20
        IER = -53
        GO TO 9010
   20 CONTINUE
C
      IF(NDIMV .GE. NRNG) GO TO 30
        IER = -4
        GO TO 9010
   30 CONTINUE
C
      IF(NCOEF .GE. KORD) GO TO 40
        IER = -6
        GO TO 9010
   40 CONTINUE
C
      IF(NDOM .EQ. 1) GO TO 50
        IER = -51
        GO TO 9010
   50 CONTINUE
C
C     ==================================================================
C     PROCESS B-SPLINE CASE
C     ==================================================================
C
      IF( .NOT.RATNL) THEN
C
C     ---------------------
C     ... CHECK VALID NWORK
C     ---------------------
C
      NEED = KORD**2 + (NDER+3) * KORD
      IF(NWORK .GE. NEED) GO TO 60
        IER = -3
        GO TO 9020
   60 CONTINUE
C
C     ---------------
C     ... CALL DTSPD1
C     ---------------
C
          IKTS  = 6
          ICOEF = IKTS + NCOEF + KORD
          INCK  = 1
          INCC  = 1
          NDIMC = NCOEF
C
        CALL DTSPD1(X,NDER,KORD,C(IKTS),INCK,NCOEF,C(ICOEF),INCC,
     1              NDIMC,NRNG,ISPAN,WORK,V,NDIMV,IER)
        IF (IER.EQ.-38) GO TO 9100
        IF (IER.NE.0) GO TO 9010
C
C     ==================================================================
C     PROCESS RATIONAL B-SPLINE CASE
C     ==================================================================
C
      ELSE
C
C     ---------------------
C     ... CHECK VALID NWORK
C     ---------------------
C
      NEED = KORD**2 + (NDER+3) * KORD + NRAT * (NDER+1)
      IF(NWORK .GE. NEED) GO TO 70
        IER = -3
        GO TO 9020
   70 CONTINUE
C
C     ---------------
C     ... CALL DTSPD1
C     ---------------
C
          IKTS  = 6
          ICOEF = IKTS + NCOEF + KORD
          INCK  = 1
          INCC  = 1
          NDIMC = NCOEF
          IW1   = 1
          IW2   = IW1 + KORD**2 + (NDER+3) * KORD
C
        CALL DTSPD1(X,NDER,KORD,C(IKTS),INCK,NCOEF,C(ICOEF),INCC,
     1              NDIMC,NRAT,ISPAN,WORK(IW1),WORK(IW2),NRAT,IER)
        IF (IER.EQ.-38) GO TO 9100
        IF (IER.NE.0) GO TO 9010
        DIV = WORK(IW2+NRAT-1)
        IF (ABS(DIV) .LE. EPS) THEN
            IER = -10
            GOTO 9010
        ENDIF
        DO 500 I = 1, NRNG
            V(I,1) = WORK(IW2+I-1) / DIV
            DO 100 J = 1, NDER
                WORK(IW1+J-1) = 1.0
  100       CONTINUE
            DO 400 J = 1, NDER
                SUM = WORK(IW2+J*NRAT+I-1)
                DO 200 K = 1, J
                    SUM = SUM - WORK(IW1+K-1)
     1                        * V(I,J-K+1) * WORK(IW2+(K+1)*NRAT-1)
  200           CONTINUE
                V(I,J+1) = SUM / DIV
                DO 300 K = 1, J
                    WORK(IW1+K-1) = DBLE(J+1)/DBLE(J+1-K)
     1                            * WORK(IW1+K-1)
  300           CONTINUE
  400       CONTINUE
  500   CONTINUE
C
C     ==================================================================
C
      ENDIF
C
C  NORMAL RETURN
C
 9000 CONTINUE
      IER = 0
      C(5) = ISPAN
      RETURN
C
C  ERROR RETURNS.
C
 9010 CONTINUE
      MODE = 1
      GO TO 9900
 9020 CONTINUE
      MODE = 2
      GO TO 9900
 9100 CONTINUE
      MODE = 3
 9900 CONTINUE
      V(1,1) = DTMCON(1)
      CALL DTERR(MODE,SUBNAM,IER,NEED)
      RETURN
C
      END
