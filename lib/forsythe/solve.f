C+----------------------------------------------------------------------
C
      SUBROUTINE SOLVE(NDIM, N, A, B, IPVT)
C
      INTEGER NDIM, N, IPVT(N)
      REAL    A(NDIM,N),B(N)
C
C   SOLUTION OF LINEAR SYSTEM, A*X = B .
C   DO NOT USE IF DECOMP HAS DETECTED SINGULARITY.
C
C   INPUT..
C
C     NDIM = DECLARED ROW DIMENSION OF ARRAY CONTAINING A .
C
C     N = ORDER OF MATRIX.
C
C     A = TRIANGULARIZED MATRIX OBTAINED FROM DECOMP .
C
C     B = RIGHT HAND SIDE VECTOR.
C
C     IPVT = PIVOT VECTOR OBTAINED FROM DECOMP .
C
C   OUTPUT..
C
C     B = SOLUTION VECTOR, X .
C
C-----------------------------------------------------------------------
C
      INTEGER KB, KM1, NM1, KP1, I, K, M
      REAL    T
C
C     FORWARD ELIMINATION
C
      IF (N .EQ. 1) GO TO 50
      NM1 = N-1
      DO 20 K = 1, NM1
         KP1 = K+1
         M = IPVT(K)
         T = B(M)
         B(M) = B(K)
         B(K) = T
         DO 10 I = KP1, N
             B(I) = B(I) + A(I,K)*T
   10    CONTINUE
   20 CONTINUE
C
C     BACK SUBSTITUTION
C
      DO 40 KB = 1,NM1
         KM1 = N-KB
         K = KM1+1
         B(K) = B(K)/A(K,K)
         T = -B(K)
         DO 30 I = 1, KM1
             B(I) = B(I) + A(I,K)*T
   30    CONTINUE
   40 CONTINUE
   50 B(1) = B(1)/A(1,1)
      RETURN
      END
