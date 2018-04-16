      SUBROUTINE D0PIT2(ALU,NDIM,Y,N,NCOEF,IPVT,KORD,MLT,NLHS,INDLHS,
     +                  BCLHS,INDL,BCL,NRHS,INDRHS,BCRHS,INDR,BCR,BSPC)
C
C=======================================================================
C     PARAMETERS
C=======================================================================
C
      INTEGER IPVT(*), INDL(*), INDLHS(*), INDR(*), INDRHS(*), KORD
      INTEGER N, NCOEF, NDIM, NLHS, NRHS, MLT(*)
      DOUBLE PRECISION ALU(NDIM,*), BCL(*), BCLHS(*), BCR(*), BCRHS(*)
      DOUBLE PRECISION BSPC(*), Y(*)
C
C=======================================================================
C     INTERNAL VARIABLES
C=======================================================================
C
      INTEGER I, ID, IR, JOB, KP1, L, ML, NINT, NM2, NR
C
C=======================================================================
C    DEFINE INTEGERS
C=======================================================================
C
      ML = KORD - 1
      KP1 = KORD + 1
      NM2 = N - 2
      JOB = 0
C
C=======================================================================
C    SET LEFT BOUNDARY CONDITIONS
C=======================================================================
C
      BCL(1) = Y(1)
      IF( NLHS .EQ. 0 ) GO TO 20
      DO 10 I = 1, NLHS
          ID = INDLHS(I) + 1
          BCL(ID) = BCLHS(I)
10    CONTINUE
20    IR = 0
      DO 30 I = 1, KORD
          ID = KP1 - I
          IF( INDL(ID) .EQ. 0 ) GO TO 30
          IR = IR + 1
          BSPC(IR) = BCL(ID)
30    CONTINUE
C
C=======================================================================
C    SET INTERPOLATION POINTS
C=======================================================================
C
      NR = 1
      IF( NM2 .EQ. 0 ) GO TO 70
      DO 60 I = 1, NM2
          IF( MLT(I) .LT. 0 ) GO TO 50
          NINT = MLT(I)
          DO 40 L = 1, NINT
              IR = IR + 1
              NR = NR + 1
              BSPC(IR) = Y(NR)
40        CONTINUE
          GO TO 60
50        NR = NR + IABS( MLT(I) )
60    CONTINUE
C
C=======================================================================
C    SET RIGHT BOUNDARY CONDIRIONS
C=======================================================================
C
70    BCR(1) = Y(NR+1)
      IF( NRHS .EQ. 0 ) GO TO 90
      DO 80 I = 1, NRHS
          ID = INDRHS(I) + 1
          BCR(ID) = BCRHS(I)
80    CONTINUE
90    DO 100 I = 1, KORD
          IF( INDR(I) .EQ. 0 ) GO TO 100
          IR = IR + 1
          BSPC(IR) = BCR(I)
100   CONTINUE
C
C=======================================================================
C     SOLVE FOR B-SPLINE COEFFICIENTS
C=======================================================================
C
      CALL DGBSL(ALU,NDIM,NCOEF,ML,ML,IPVT,BSPC,JOB)
      RETURN
      END
