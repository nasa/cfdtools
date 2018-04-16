      SUBROUTINE D0PIT1(N,X,Y,XKNTS,NCOEF,KORD,MLT,NLHS,INDLHS,
     +             NRHS,INDRHS,NDIM,A,WORK,BSPL,INDL,BCL,
     +             INDR,BCR)
C
C=======================================================================
C    PARAMETERS
C=======================================================================
C
      INTEGER INDL(*), INDLHS(*), INDR(*), INDRHS(*), KORD, MLT(*)
      INTEGER N, NCOEF, NDIM, NLHS, NRHS
      DOUBLE PRECISION A(NDIM,*), BCL(*), BCR(*), BSPL(KORD,*)
      DOUBLE PRECISION WORK(*), X(*)
      DOUBLE PRECISION XKNTS(*), Y(*)
C
C=======================================================================
C    INTERNAL PARAMETERS
C=======================================================================
C
      INTEGER I, ID, IR, IS, IPOS, JC, JR, J, KHALF, KP1
      INTEGER L, ML, M, NINT, NKNTS, NM2
      NKNTS = NCOEF + KORD
      KHALF = KORD / 2
      ML = KORD - 1
      M = 2 * ML + 1
      NM2 = N - 2
      DO 20 J = 1, NCOEF
          DO 10 I = 1, NDIM
             A(I,J) = 0.0D0
10        CONTINUE
20    CONTINUE
      KP1 = KORD + 1
C
C=======================================================================
C    CALL DTBSP1 TO EVALUATE B-SPLINES
C=======================================================================
C
      IPOS = KORD
      CALL DTBSP1(XKNTS,X(1),IPOS,KORD,ML,WORK,BSPL,KORD)
C
C=======================================================================
C    DETERMINE WHICH DERIVATIVE TO USE
C=======================================================================
C
      DO 30 I = 1, KHALF
          INDL(I) = 0
          BCL(I)  = 0.0D0
          IS = I + KHALF
          INDL(IS) = 1
          BCL(IS) = 0.0D0
30    CONTINUE
      INDL(1) = 1
      INDL(KORD) = 0
C
C=======================================================================
C    SOME LEFT-HAND BOUNDARY CONDITIONS PROVIDED RESET INDL AND BCL
C=======================================================================
C
      IF( NLHS .EQ. 0 ) GO TO 50
      DO 40 I = 1, NLHS
          ID = INDLHS(I) + 1
          INDL(ID) = 1
          IS = KP1 - ID
          INDL(IS) = 0
40    CONTINUE
C
C=======================================================================
C    DEFINE FIRST KHALF ROWS OF A
C=======================================================================
C
50    IR = 0
      DO 70 I = 1, KORD
          ID = KP1 - I
          IF( INDL(ID) .EQ. 0 ) GO TO 70
          IR = IR + 1
          DO 60 J = 1, KORD
              JR = IR - J + M
              A(JR,J) = BSPL(J,ID)
60        CONTINUE
70    CONTINUE
C
C=======================================================================
C    DEFINE THE INTERIOR ROWS OF A TO REPRESENT THE INTERPOLATION
C=======================================================================
C
      IF( NM2 .EQ. 0 ) GO TO 110
      DO 100 I = 1, NM2
          IF( MLT(I) .LT. 0 ) GO TO 100
          NINT = MLT(I) - 1
          IPOS = IPOS + MLT(I)
          CALL DTBSP1(XKNTS,X(I+1),IPOS,KORD,NINT,WORK,BSPL,KORD)
          NINT = MLT(I)
          DO 90 L = 1, NINT
              IR = IR + 1
              JC = IPOS - KORD
              DO 80 J = 1, KORD
                  JC = JC + 1
                  JR = IR - JC + M
                  A(JR,JC) = BSPL(J,L)
80            CONTINUE
90        CONTINUE
100   CONTINUE
C
C=======================================================================
C    DETERMINE WHICH DERIVATIVES ON RIGHT HAND SIDE
C=======================================================================
C
110   CALL DTBSP1(XKNTS,X(N),NCOEF,KORD,ML,WORK,BSPL,KORD)
      DO 120 I = 1, KHALF
          INDR(I) = 0
          BCR(I) = 0.0D0
          IS = KHALF + I
          INDR(IS) = 1
          BCR(IS) = 0.0D0
120   CONTINUE
      INDR(1) = 1
      INDR(KORD) = 0
C
C=======================================================================
C    SOME RIGHT HAND BOUNDARY CONDITIONS GIVEN RESET IND AND BC
C=======================================================================
C
      IF( NRHS .EQ. 0 ) GO TO 140
      DO 130 I = 1, NRHS
          ID = INDRHS(I) + 1
          INDR(ID) = 1
          IS = KP1 - ID
          INDR(IS) = 0
130   CONTINUE
C
C=======================================================================
C    NOW SET LAST KHALF ROWS OF A
C=======================================================================
C
140   DO 160 I = 1, KORD
          IF( INDR(I) .EQ. 0 ) GO TO 160
          JC = NCOEF - KORD
          IR = IR + 1
          DO 150 J = 1, KORD
              JC = JC + 1
              JR = IR - JC + M
              A(JR,JC) = BSPL(J,I)
150       CONTINUE
160   CONTINUE
C
C=======================================================================
C    RETURN
C=======================================================================
C
      RETURN
      END
