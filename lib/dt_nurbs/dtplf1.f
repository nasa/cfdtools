      SUBROUTINE DTPLF1(NPTS,T,Y,IW,WHT,NCOEF,KORD,XKNOTS,BSVAL,
     +                 WK1,WK2,ATA,ATB)
C
C=====================================================================
C    PARAMETERS
C=====================================================================
C
      INTEGER IW, KORD, NCOEF, NPTS
      DOUBLE PRECISION    ATA(KORD,*), ATB(*), BSVAL(*), T(*), WHT(*)
      DOUBLE PRECISION    WK1(*), WK2(*), XKNOTS(*), Y(*)
C
C=====================================================================
C    INTERNAL VARIABLES
C=====================================================================
C
      INTEGER I, IK, IPOS, IR, J, JC, K, NDERV
      DOUBLE PRECISION    W, WW
C
C=====================================================================
C    INITIALIZE VARIBLES AND ARRAYS
C=====================================================================
C
      NDERV = 0
      DO 20 J = 1, NCOEF
          ATB(J) = 0.0D0
          DO 10 I = 1, KORD
              ATA(I,J) = 0.0D0
10        CONTINUE
20    CONTINUE
C
C=======================================================================
C    LOOP THROUGH POINTS DEFINING CONTRIBUTION TO ATA AND ATB
C=======================================================================
C
      IPOS = KORD
      DO 70 K = 1, NPTS
C
C . . . DETERMINE POSITION PARAMETER
C
30        IF( T(K) .LE. XKNOTS(IPOS+1) ) GO TO 40
          IPOS = IPOS + 1
          GO TO 30
C
C . . . DETERMINE IF POINT IS TO BE USED
C
40        W = 1.0D0
          IF( IW .NE. 0 ) W = WHT(K)
          IF( W .LE. 0.0D0 ) GO TO 70
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
          DO 60 I = 1, KORD
              JC = JC + 1
              DO 50 J = 1, I
                  IR = J - I + KORD
                  ATA(IR,JC) = ATA(IR,JC) + BSVAL(J) * BSVAL(I) * WW
50            CONTINUE
              ATB(JC) = ATB(JC) + BSVAL(I) * Y(K) * WW
60        CONTINUE
70    CONTINUE
C
C=======================================================================
C    CHECK THAT ENOUGH POINTS ARE ACTIVE
C=======================================================================
C
      RETURN
      END
