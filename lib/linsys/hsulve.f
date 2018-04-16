C+----------------------------------------------------------------------
C
      SUBROUTINE HSULVE ( NDIM, N, M, A, U, B )
C
C     LEAST SQUARES SOLUTION OF UNDERDETERMINED SYSTEMS, M.LE.N.
C     SOLVES TRIANGULAR SYSTEM FROM HDECOM GIVING MINIMAL-LENGTH
C     SOLUTION TO   AX = B,  WHERE THE CALL TO  HDECOM  INVOLVED
C     A-TRANSPOSE WITH  M  AND  N  REVERSED.  THUS, FOR THE CALL
C     TO  HSULVE,
C
C     PARAMETERS  NDIM,N,M,A,U  SHOULD MATCH EXACTLY THE LIST
C                 MDIM,M,N,A,U  IN THE CALL TO  HDECOM,  WITH
C                 B = N-VECTOR, INPUT  WITH RIGHT-HAND-SIDE
C                                      IN ELEMENTS B(1) TO B(M);
C                               OUTPUT WITH THE SOLUTION
C                                      IN ELEMENTS B(1) TO B(N);
C                 U = N-VECTOR, INPUT  WITH TRANSFORMATION INFO.
C                                      IN ELEMENTS U(1) TO U(M).
C
C     REFERENCE:    M. A. SAUNDERS, SYSTEMS OPTIMIZATION LAB, STANFORD.
C
C     PROGRAMMER:   D. A. SAUNDERS, INFORMATICS INC, 1975.
C
C-----------------------------------------------------------------------
C
      IMPLICIT REAL ( A-H, O-Z )
C
      DIMENSION A(NDIM,M), U(N), B(N)
C
C  *  FORWARD SUBSTITUTION FOR  LY = B.  OVERWRITE  B  WITH  Y
      DO 30 I=1,M
         T = 0.E+0
         IF ( I.EQ.1 )  GOTO 25
            IM1 = I-1
            DO 20 J=1,IM1
               T = A(J,I)*B(J) + T
 20         CONTINUE
 25      B(I) = ( B(I) - T )/A(I,I)
 30   CONTINUE
C
C  *  FOR SHORTEST LENGTH SOLUTION, FILL  B  OUT WITH ZEROES,
C  *  AND APPLY REFLECTIONS IN REVERSE ORDER:
      MP1 = M+1
      IF ( M.EQ.N )  GOTO 45
C
      DO 40 K=MP1,N
         B(K)=0.E+0
 40   CONTINUE
C
 45   DO 70 KB=1,M
         K = MP1 - KB
         T = A(K,K)
         A(K,K)= U(K)
         BETA  =-U(K)*T
         GAMMA = 0.E+0
         DO 50 I=K,N
            GAMMA = A(I,K)*B(I) + GAMMA
 50      CONTINUE
         GAMMA = GAMMA/BETA
         DO 60 I=K,N
            B(I) = B(I) - A(I,K)*GAMMA
 60      CONTINUE
         A(K,K)= T
 70   CONTINUE
C
C  *  SOLUTION IN B(I), I=1,...,N.
      RETURN
      END
