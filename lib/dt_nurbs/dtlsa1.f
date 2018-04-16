      SUBROUTINE DTLSA1(NPTS,T,Y,IWT,WHT,NCOEF,KORD,XKNOTS,WORK,
     &                  ATA)
C
C=====================================================================
C
C  Parameters:
C
      DOUBLE PRECISION ZERO, ONE
      PARAMETER  (ZERO=0.D0,ONE=1.D0)
C
C  Arguments:
C
      INTEGER  IWT, KORD, NCOEF, NPTS
      DOUBLE PRECISION    ATA(KORD,*), WORK(*), T(NPTS), WHT(NPTS),
     &                    XKNOTS(*), Y(NPTS)
C
C  Internal:
C
      INTEGER  I, IK, IPOS, IR, I1, I2, J, JC, K, NDERV
      DOUBLE PRECISION    W, WW
C
C=====================================================================
C    Initialize variables and arrays.
C=====================================================================
C
      NDERV = 0
      DO 20 J = 1, NCOEF
        DO 10 I = 1, KORD
          ATA(I,J) = ZERO
   10   CONTINUE
   20 CONTINUE
C
C=======================================================================
C    Loop through points defining contribution to ATA.
C=======================================================================
C
      I1 = 1 + KORD
      I2 = I1 + KORD
      IPOS = KORD
      DO 70 K = 1, NPTS
C
C       Determine position parameter.
C
   30     IF( T(K) .LE. XKNOTS(IPOS+1) ) GO TO 40
          IPOS = IPOS + 1
          GO TO 30
C
C       Determine if point is to be used.
C
   40     W = ONE
          IF( IWT .NE. 0 ) W = WHT(K)
          IF( W .LE. ZERO ) GO TO 70
          WW   = W * W
C
C       Call DTBSP2 to evaluate b-splines.
C
          IK = -1
          CALL DTBSP2(XKNOTS,T(K),IPOS,IK,KORD,WORK(I1),WORK(I2),WORK)
C
C       Add contribution of this point.
C
          JC = IPOS - KORD
          DO 60 I = 1, KORD
            JC = JC + 1
            DO 50 J = 1, I
              IR = J - I + KORD
              ATA(IR,JC) = ATA(IR,JC) + WORK(J)*WORK(I)*WW
   50       CONTINUE
   60     CONTINUE
   70 CONTINUE
C
C=======================================================================
C    CHECK THAT ENOUGH POINTS ARE ACTIVE
C=======================================================================
C
      RETURN
      END
