      SUBROUTINE DTLSA2(NPTS,T,Y,IWT,WHT,NCOEF,KORD,XKNOTS,WORK,
     +                  ATB,IER)
C
C=====================================================================
C  Parameters:
C
      DOUBLE PRECISION ZERO, ONE
      PARAMETER  (ZERO=0.D0,ONE=1.D0)
C
C  Arguments:
C
      INTEGER IER, IWT, KORD, NCOEF, NPTS
      DOUBLE PRECISION    ATB(*), WORK(*), T(*), WHT(*)
      DOUBLE PRECISION    XKNOTS(*), Y(*)
C
C  Internal:
C
      INTEGER  I, IK, IPOS, I1, I2, J, JC, K, NEND
      DOUBLE PRECISION    W, WW
C
C=====================================================================
C    Initialize variables and arrays.
C=====================================================================
C
      NEND  = NCOEF + 1
      DO 10 J = 1, NCOEF
        ATB(J) = ZERO
10    CONTINUE
C
C=======================================================================
C    Loop through points defining contribution to ATB.
C=======================================================================
C
      I1 = 1 + KORD
      I2 = I1 + KORD
      IPOS = KORD
      DO 50 K = 1, NPTS
C
C             Verify T(K) is still in range.
C
        IF ( XKNOTS(KORD).GT.T(K) .OR. T(K).GT.XKNOTS(NEND) ) THEN
          IER = -37
          RETURN
        ENDIF
C
C             Determine position parameter.
C
   20   IF( T(K) .LE. XKNOTS(IPOS+1) ) GO TO 30
        IPOS = IPOS + 1
        GO TO 20
C
C             Determine if point is to be used.
C
   30   W = ONE
        IF( IWT .NE. 0 ) W = WHT(K)
        IF( W .LE. ZERO ) GO TO 50
        WW   = W * W
C
C             Call DTBSP2 to evaluate b-splines.
C
        IK = -1
        CALL DTBSP2(XKNOTS,T(K),IPOS,IK,KORD,WORK(I1),WORK(I2),WORK)
C
C             Add contribution of this point.
C
        JC = IPOS - KORD
        DO 40 I = 1, KORD
          JC = JC + 1
          ATB(JC) = ATB(JC) + WORK(I)*Y(K)*WW
   40   CONTINUE
   50 CONTINUE
C
C=======================================================================
C    Check that enough points are active.
C=======================================================================
C
      RETURN
      END
