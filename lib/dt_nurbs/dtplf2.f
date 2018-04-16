      SUBROUTINE DTPLF2(NPTS,T,Y,IW,WHT,NCOEF,KORD,XKNOTS,BSVAL,
     +                 WK1,WK2,ATB,IER)
C
C=====================================================================
C    PARAMETERS
C=====================================================================
C
      INTEGER IER, IW, KORD, NCOEF, NPTS
      DOUBLE PRECISION    ATB(*), BSVAL(*), T(*), WHT(*), WK1(*), WK2(*)
      DOUBLE PRECISION    XKNOTS(*), Y(*)
C
C=====================================================================
C    INTERNAL VARIABLES
C=====================================================================
C
      INTEGER I, IK, IPOS, J, JC, K, NACT, NEND, NDERV
      DOUBLE PRECISION    W, WW
C
C=====================================================================
C    INITIALIZE VARIBLES AND ARRAYS
C=====================================================================
C
      NEND  = NCOEF + 1
      NDERV = 0
      DO 10 J = 1, NCOEF
          ATB(J) = 0.0D0
10    CONTINUE
C
C=======================================================================
C    LOOP THROUGH POINTS DEFINING CONTRIBUTION TO ATA AND ATB
C=======================================================================
C
      IPOS = KORD
      DO 50 K = 1, NPTS
C
C . . . VERIFY T(K) IS STILL IN RANGE
C
          IF( XKNOTS(KORD) .LE. T(K) .AND. T(K) .LE. XKNOTS(NEND) )
     +    GO TO 20
          IER = -37
          GO TO 50
C
C . . . DETERMINE POSITION PARAMETER
C
20        IF( T(K) .LE. XKNOTS(IPOS+1) ) GO TO 30
          IPOS = IPOS + 1
          GO TO 20
C
C . . . DETERMINE IF POINT IS TO BE USED
C
30        W = 1.0D0
          IF( IW .NE. 0 ) W = WHT(K)
          IF( W .LE. 0.0D0 ) GO TO 50
          NACT = NACT + 1
          WW   = W * W
C
C . . . CALL DTBSP2 TO EVALUATE B-SPLINES
C
          IK = -1
          CALL DTBSP2(XKNOTS,T(K),IPOS,IK,KORD,WK1,WK2,BSVAL)
C
C . . . ADD CONTRIBUTION OF THIS POINT
C
          JC = IPOS - KORD
          DO 40 I = 1, KORD
              JC = JC + 1
              ATB(JC) = ATB(JC) + BSVAL(I) * Y(K) * WW
40        CONTINUE
50    CONTINUE
      RETURN
      END
