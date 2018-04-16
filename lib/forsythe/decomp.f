C+----------------------------------------------------------------------
C
      SUBROUTINE DECOMP(NDIM,N,A,COND,IPVT,WORK)
C
      INTEGER NDIM,N
      REAL A(NDIM,N),COND,WORK(N)
      INTEGER IPVT(N)
C
C     DECOMPOSES A SINGLE PRECISION MATRIX BY GAUSSIAN ELIMINATION
C     AND ESTIMATES THE CONDITION OF THE MATRIX.
C
C     USE SOLVE TO COMPUTE SOLUTIONS TO LINEAR SYSTEMS.
C
C     INPUT..
C
C        NDIM = DECLARED ROW DIMENSION OF THE ARRAY CONTAINING  A.
C
C        N = ORDER OF THE MATRIX.
C
C        A = MATRIX TO BE TRIANGULARIZED.
C
C     OUTPUT..
C
C        A  CONTAINS AN UPPER TRIANGULAR MATRIX  U  AND A PERMUTED
C          VERSION OF A LOWER TRIANGULAR MATRIX  I-L  SO THAT
C          (PERMUTATION MATRIX)*A = L*U .
C
C        COND = AN ESTIMATE OF THE CONDITION OF  A .
C           FOR THE LINEAR SYSTEM  A*X = B, CHANGES IN  A  AND  B
C           MAY CAUSE CHANGES  COND  TIMES AS LARGE IN  X .
C           IF COND+1 = COND, A IS SINGULAR TO WORKING PRECISION
C           COND = 1.0E+32  IF EXACT SINGULARITY IS DETECTED.
C
C        IPVT = THE PIVOT VECTOR.
C           IPVT(K) = THE INDEX OF THE K-TH PIVOT ROW
C           IPVT(N) = (-1)**(NUMBER OF INTERCHANGES)
C
C     WORK SPACE..  THE VECTOR  WORK  MUST BE DECLARED AND INCLUDED
C                   IN THE CALL.  ITS INPUT CONTENTS ARE IGNORED.
C                   ITS OUTPUT CONTENTS ARE USUALLY UNIMPORTANT.
C
C     THE DETERMINANT OF A CAN BE OBTAINED ON OUTPUT BY
C        DET(A) = IPVT(N) * A(1,1) * A(2,2) * ... * A(N,N).
C
C-----------------------------------------------------------------------
C
      REAL    EK, T, ANORM, YNORM, ZNORM
      INTEGER NM1, I, J, K, KP1, KB, KM1, M
      INTRINSIC ABS, SIGN
C
      IPVT(N) = 1
      IF (N .EQ. 1) GO TO 80
      NM1 = N - 1
C
C     COMPUTE 1-NORM OF A
C
      ANORM = 0.0E+0
      DO 10 J = 1, N
         T = 0.0E+0
         DO 5 I = 1, N
            T = T + ABS(A(I,J))
    5    CONTINUE
         IF (T .GT. ANORM) ANORM = T
   10 CONTINUE
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C
      DO 35 K = 1,NM1
         KP1= K+1
C
C        FIND PIVOT
C
         M = K
         DO 15 I = KP1,N
            IF (ABS(A(I,K)) .GT. ABS(A(M,K))) M = I
   15    CONTINUE
         IPVT(K) = M
         IF (M .NE. K) IPVT(N) = -IPVT(N)
         T = A(M,K)
         A(M,K) = A(K,K)
         A(K,K) = T
C
C        SKIP STEP IF PIVOT IS ZERO
C
         IF (T .EQ. 0.0E+0) GO TO 35
C
C        COMPUTE MULTIPLIERS
C
         DO 20 I = KP1,N
             A(I,K) = -A(I,K)/T
   20    CONTINUE
C
C        INTERCHANGE AND ELIMINATE BY COLUMNS
C
         DO 30 J = KP1,N
             T = A(M,J)
             A(M,J) = A(K,J)
             A(K,J) = T
             IF (T .EQ. 0.0E+0) GO TO 30
             DO 25 I = KP1,N
                A(I,J) = A(I,J) + A(I,K)*T
   25        CONTINUE
   30    CONTINUE
   35 CONTINUE
C
C     COND = (1-NORM OF A)*(AN ESTIMATE OF 1-NORM OF A-INVERSE)
C     ESTIMATE OBTAINED BY ONE STEP OF INVERSE ITERATION FOR THE
C     SMALL SINGULAR VECTOR.  THIS INVOLVES SOLVING TWO SYSTEMS
C     OF EQUATIONS, (A-TRANSPOSE)*Y = E  AND  A*Z = Y  WHERE  E
C     IS A VECTOR OF +1 OR -1 CHOSEN TO CAUSE GROWTH IN Y.
C     ESTIMATE = (1-NORM OF Z)/(1-NORM OF Y)
C
C     SOLVE (A-TRANSPOSE)*Y = E
C
      DO 50 K = 1, N
         T = 0.0E+0
         IF (K .EQ. 1) GO TO 45
         KM1 = K-1
         DO 40 I = 1, KM1
            T = T + A(I,K)*WORK(I)
   40    CONTINUE
   45    EK = 1.0E+0
         IF (T .LT. 0.0E+0) EK = -1.0E+0
         IF (A(K,K) .EQ. 0.0E+0) GO TO 90
         WORK(K) = -(EK + T)/A(K,K)
   50 CONTINUE
      DO 60 KB = 1, NM1
         K = N - KB
         T = 0.0E+0
         KP1 = K+1
         DO 55 I = KP1, N
            T = T + A(I,K)*WORK(K)
   55    CONTINUE
         WORK(K) = T
         M = IPVT(K)
         IF (M .EQ. K) GO TO 60
         T = WORK(M)
         WORK(M) = WORK(K)
         WORK(K) = T
   60 CONTINUE
C
      YNORM = 0.0E+0
      DO 65 I = 1, N
         YNORM = YNORM + ABS(WORK(I))
   65 CONTINUE
C
C     SOLVE A*Z = Y
C
      CALL SOLVE(NDIM, N, A, WORK, IPVT)
C
      ZNORM = 0.0E+0
      DO 70 I = 1, N
         ZNORM = ZNORM + ABS(WORK(I))
   70 CONTINUE
C
C     ESTIMATE CONDITION
C
      COND = ANORM*ZNORM/YNORM
      IF (COND .LT. 1.0E+0) COND = 1.0E+0
      RETURN
C
C     1-BY-1
C
   80 COND = 1.0E+0
      IF (A(1,1) .NE. 0.0E+0) RETURN
C
C     EXACT SINGULARITY
C
   90 COND = 1.0E+32
      RETURN
      END
