C+----------------------------------------------------------------------
C
      SUBROUTINE HSOLVE ( MDIM, M, N, A, U, B, RESIDL )
C
C     LEAST SQUARES SOLUTION OF OVERDETERMINED SYSTEMS.
C     SOLVES TRIANGULAR SYSTEM FROM  HDECOM, GIVING  X  WHICH
C     MINIMIZES THE 2-NORM OF   AX - B.
C
C     USE  HSULVE  FOR UNDERDETERMINED SYSTEMS.
C
C     MDIM,M,N,A,U: RESULTS FROM  HDECOM.
C     B:            M-VECTOR, INPUT WITH RIGHT-HAND-SIDE;
C                   OUTPUT WITH SOLUTION  X  IN B(1) TO B(N),
C                   AND TRANSFORMED RESIDUAL IN REST OF B.
C     RESIDL:       OUTPUT WITH MINIMUM SUM OF SQUARES (=0 IF M=N).
C
C     NOTE:         DIVISION BY ZERO IMPLIES  A  IS RANK-DEFICIENT.
C
C     PROGRAMMER:   DAVID SAUNDERS, INFORMATICS INC, 1975.
C
C-----------------------------------------------------------------------
C
      IMPLICIT REAL ( A-H, O-Z )
C
      DIMENSION A(MDIM,N), U(M), B(M)
C
C  *  APPLY REFLECTIONS TO  B:
      DO 8 K=1,N
         T = A(K,K)
         A(K,K)= U(K)
         BETA  =-U(K)*T
         GAMMA = 0.E+0
         DO 6 I=K,M
            GAMMA = A(I,K)*B(I) + GAMMA
 6       CONTINUE
         GAMMA = GAMMA/BETA
         DO 7 I=K,M
            B(I) = B(I) - GAMMA*A(I,K)
 7       CONTINUE
         A(K,K) = T
 8    CONTINUE
C
C  *  PERFORM BACK-SUBSTITUTION:
      NP1 = N+1
      DO 10 KB=1,N
         K = NP1 - KB
         B(K) = B(K)/A(K,K)
         IF ( K.EQ.1 )  GOTO 10
            KM1 = K-1
            DO 9 I=1,KM1
               B(I) = B(I) - A(I,K)*B(K)
 9          CONTINUE
 10   CONTINUE
C
      RESIDL = 0.E+0
      IF ( M.EQ.N )  GOTO 99
C
      DO 11 I=NP1,M
         RESIDL = B(I)**2 + RESIDL
 11   CONTINUE
C
C  *  REQUIRED SOLUTION IN  B(I), I=1,...,N.
C  *  MINIMUM SUM OF SQUARES IN  RESIDL.
 99   RETURN
      END
