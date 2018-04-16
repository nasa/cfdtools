C<*>

      SUBROUTINE DTSPVL(X,C,WORK,NWORK,V,IER)
C
C  PURPOSE:  CALCULATE VALUE ONLY OF A UNIVARIATE SPLINE
C            GIVEN A PACKED C ARRAY.
C
C  INPUTS:
C           X       VALUE AT WHICH TO EVALUATE.
C
C           C       SPLINE DEFINITION ARRAY.
C
C  WORKING STORAGE:
C
C           WORK    WORK ARRAY OF LENGTH NWORK.
C
C           NWORK   LENGTH OF WORK ARRAY,
C                   NWORK .GE. 5*C(3) - 2.
C
C  OUTPUT:
C           V       ARRAY OF VALUES OF THE DEPENDENT VARIABLES.
C                   V IS OF LENGTH ABS(C(2))
C
C           IER     SUCCESS/ERROR CODE.
C                   IER=0     SUCCESS.
C                   IER=-1    C(3) .LE. 0.
C                   IER=-3    NWORK TOO SMALL.
C                   IER=-6    C(4) .LT. C(3).
C                   IER=-8    INVALID KNOT SET.
C                   IER=-10   DENOMINATOR .EQ. 0 ON A RATIONAL SPLINE.
C                   IER=-38   ATTEMPT TO EVALUATE AT POINT THAT IS
C                             INSIDE AN INTERVAL THAT IS TOO SMALL.
C                   IER=-50   X OUT OF RANGE.
C                   IER=-51   C(1) .NE. 1.
C                   IER=-52   C(2) .EQ. 0. OR -1.
C
C  SUBROUTINES CALLED:  DTSPV1, DTMDON, DTERR
C
C  AUTHOR:  A.K. JONES
C
C  DATE:  1-NOV-84
C
C  MINOR REVISIONS BY FRITZ KLEIN, 21-JAN-85
C
C  MODIFIED BY DEBORAH PARSONS, MARCH 28, 1989
C     TO INCLUDE RATIONAL B-SPLINES
C
C  MODIFIED BY J. W. MANKE, 9/3/91, FOR RATIONAL B-SPLINE CERTIFICATION
C
      DOUBLE PRECISION X,C(*),DTMCON,V(*),WORK(*)
      DOUBLE PRECISION EPS, DIV
      EXTERNAL DTSPV1,DTMCON,DTERR
      CHARACTER*8 SUBNAM
      LOGICAL RATNL
C
      DATA SUBNAM /'DTSPVL  '/
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
C  ERROR-CHECK INPUTS.
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
      IF( KORD .GT. 0) GO TO 10
        IER = -1
      GO TO 9010
   10 CONTINUE
C
      IF(NCOEF .GE. KORD) GO TO 20
        IER = -6
      GO TO 9010
   20 CONTINUE
C
      IF(NDOM .EQ. 1) GO TO 30
        IER = -51
        GO TO 9010
   30 CONTINUE
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
          NEED = 5 * KORD - 2
          IF(NWORK .GE. NEED) GO TO 40
            IER = -3
          GO TO 9020
   40     CONTINUE
C
C     ---------------
C     ... CALL DTSPV1
C     ---------------
C
          IKTS  = 6
          ICOEF = IKTS + NCOEF + KORD
          INCK  = 1
          INCC  = 1
          NDIMC = NCOEF
C 
          CALL DTSPV1(X,KORD,C(IKTS),INCK,NCOEF,C(ICOEF),INCC,NDIMC,
     1              NRNG,ISPAN,WORK,V,IER)
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
          NEED = 5 * KORD - 2 + NRAT
          IF(NWORK .GE. NEED) GO TO 50
            IER = -3
          GO TO 9020
   50     CONTINUE
C
C     ---------------
C     ... CALL DTSPV1
C     ---------------
C
          IKTS  = 6
          ICOEF = IKTS + NCOEF + KORD
          INCK  = 1
          INCC  = 1
          NDIMC = NCOEF
          IW1   = 1
          IW2   = IW1 + 5 * KORD - 2
C
          CALL DTSPV1(X,KORD,C(IKTS),INCK,NCOEF,C(ICOEF),INCC,NDIMC,
     1                NRAT,ISPAN,WORK(IW1),WORK(IW2),IER)
          IF (IER.EQ.-38) GO TO 9100
          IF (IER.NE.0) GO TO 9010
          DIV = WORK(IW2+NRAT-1)
          IF (ABS(DIV) .LE. EPS) THEN
              IER = -10
              GOTO 9010
          ENDIF
          DO 100 I = 1, NRNG
              V(I) = WORK(IW2+I-1) / DIV
  100     CONTINUE
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
      CALL DTERR(MODE,SUBNAM,IER,NEED)
      V(1) = DTMCON(1)
      RETURN
C
      END
