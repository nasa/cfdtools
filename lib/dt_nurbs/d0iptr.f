      SUBROUTINE D0IPTR(A,M,N,WORK,NWORK,K,IER)
C
      INTEGER M,N,NWORK,K,IER,WORK(1)
      DOUBLE PRECISION A(1)
C
C     D0IPTR TRANSPOSES A DOUBLE PRECISION MATRIX A IN PLACE.  THE
C     MATRIX A(M,N) IS TRANSPOSED IN PLACES TO A(N,M).  ON INPUT A IS
C     INTERPRETED TO BE A RECTANGULAR MATRIX HAVING M ROW AND N
C     COLUMNS.  ON OUTPUT A IS INTERPRETED TO BE A RECTANGULAR
C     MATRIX HAVING N ROWS AND M COLUMNS.
C
C     WORK IS AN ONE-DIMENSIONAL ARRAY TO STORE INFORMATION TO
C     SPEED UP THE ALGORITHM.
C
C     D0IPTR IS THE RESULT OF CONVERTING ETRANS FROM EKS/LEVELTWO.
C     ETRANS IS DOCUMENTS IN NA-60.  ETRANS IS ALGORITHM 513 WHICH
C     WAS PUBLISHED IN ACM/TOMS.
C
C     D0IPTR CONVERTED BY ROGER G GRIMES  OCTOBER 19, 1979
C
C ... ERROR RETURNS
C
C      IER = 0  NORMAL RETURN, TRANSPOSITION COMPLETED
C          =-1  M .LE. 0
C          =-2  N .LE. 0
C          =-3  NWORK .LE. 0 WITH M.NE.N AND M.NE.1 AND N.NE.1
C          =-4  TRANSPOSITION FAILED.  SHOULD NEVER OCCUR
C
C ... INTERNAL VARIABLES
C
      CHARACTER*8 NAME
      DOUBLE PRECISION B,C,D
      INTEGER NCOUNT,I,IFIX,IR0,IR1,IR2,I1,I1C,I2,I2C,J,J1,KMI,MAX,MN,N1
C
      DATA NAME/'D0IPTR'/
C
C ... CHECK INPUT ARGUMENTS AND INITIALIZE
C
      IER = 0
      K = 0
      IF (M.LE.0) THEN
          IER = -1
          MODE = 1
 
      ELSE IF (N.LE.0) THEN
          IER = -2
          MODE = 1
 
      ELSE IF (M.NE.1 .AND. N.NE.1) THEN
          IF (NWORK.LE.0 .AND. M.NE.N) THEN
              IER = -3
              MODE = 2
 
          ELSE
 
C
              K = 1
              MN = M*N
              IF (M.EQ.N) THEN
C
C IF MATRIX IS SQUARE,EXCHANGE ELEMENTS A(I,J) AND A(J,I).
C
                  N1 = N - 1
                  DO 10 I = 1,N1
                      J1 = I + 1
                      DO 20 J = J1,N
                          I1 = I + (J-1)*N
                          I2 = J + (I-1)*M
                          B = A(I1)
                          A(I1) = A(I2)
                          A(I2) = B
   20                 CONTINUE
   10             CONTINUE
                  K = N1
 
              ELSE
 
                  NCOUNT = 2
                  K = MN - 1
                  IFIX = 0
                  DO 30 I = 1,NWORK
                      WORK(I) = 0
   30             CONTINUE
                  IF (M.GE.3 .AND. N.GE.3) THEN
C
C CALCULATE THE NUMBER OF FIXED POINTS, UCLIDS ALGORITHM
C FOR GCD(MN-1,N-1).
C
                      IR2 = K
                      IR1 = N - 1
   40                 CONTINUE
                      IR0 = MOD(IR2,IR1)
                      IR2 = IR1
                      IR1 = IR0
                      IF (IR0.NE.0) GO TO 40
                      NCOUNT = NCOUNT + IR2 - 1
                      IR1 = K/IR2
                      IFIX = IR1
                  END IF
C
C SET INITIAL VALUES FOR SEARCH
C
                  KMI = K - 1
                  I = 1
   50             CONTINUE
C
C AT LEAST ONE LOOP MUST BE RE-ARRANGED
C
C
C REARRANGE THE ELEMENTS OF A LOOP AND ITS COMPANION LOOP
C
                  I1 = I
                  B = A(I1+1)
                  I1C = KMI
                  C = A(I1C+1)
   60             CONTINUE
                  I2 = M*MOD(I1,N) + (I1/N)
                  I2C = K - I2
                  IF (I1.LE.NWORK) WORK(I1) = 2
                  IF (I1C.LE.NWORK) WORK(I1C) = 2
                  NCOUNT = NCOUNT + 2
                  IF (I2.EQ.I) THEN
                      GO TO 70
 
                  ELSE IF (I2.NE.KMI) THEN
                      A(I1+1) = A(I2+1)
                      A(I1C+1) = A(I2C+1)
                      I1 = I2
                      I1C = I2C
                      GO TO 60
 
                  END IF
C
C FINAL STORE AND TEST FOR FINISHED
C
                  D = B
                  B = C
                  C = D
   70             A(I1+1) = B
                  A(I1C+1) = C
                  IF (NCOUNT.LT.MN) THEN
   80                 CONTINUE
                      MAX = K - I
                      I = I + 1
                      KMI = K - I
                      IF (I.GT.MAX) THEN
                          GO TO 90
C
                      ELSE IF (I.EQ.IFIX) THEN
C
C SEARCH FOR LOOPS TO REARRANGE
C
                          IFIX = IFIX + IR1
                          GO TO 80
 
                      ELSE IF (I.GT.NWORK) THEN
                          I1 = I
  100                     CONTINUE
                          I2 = M*MOD(I1,N) + (I1/N)
                          IF (I2.GT.I .AND. I2.LT.MAX) THEN
                              I1 = I2
                              GO TO 100
 
                          END IF
 
                          IF (I2.NE.I) GO TO 80
 
                      ELSE IF (WORK(I).EQ.0) THEN
                          GO TO 50
 
                      ELSE
 
                          GO TO 80
 
                      END IF
 
                      GO TO 50
 
                  ELSE
                      GO TO 110
 
                  END IF
C
C ERROR RETURNS.
C
   90             IER = -4
                  K = I + 2
                  MODE = 3
 
C
C NORMAL RETURN
C
  110             IF (IER.EQ.0) K = I
 
              END IF
 
          END IF
 
      END IF
 
      IF (IER.NE.0) CALL DTERR(MODE,NAME,IER,1)
      RETURN
 
      END
