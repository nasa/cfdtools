C+---------------------------------------------------------------------
C
      SUBROUTINE DECSLV ( N, NDIM, AB, IERR )
C
C     DECSLV  SOLVES THE N*N SYSTEM   AX = B   BY  GAUSSIAN ELIMINATION
C     WITH PARTIAL PIVOTING.  THE RIGHT-HAND-SIDE  B  SHOULD BY APPEND-
C     ED AS THE (N+1)TH COLUMN OF THE ARRAY CONTAINING MATRIX  A.  THIS
C     COLUMN IS OVERWRITTEN WITH THE SOLUTION.   THE MATRIX  A  IS ALSO
C     OVERWRITTEN WITH THE  LU-FACTORIZATION.
C
C     IF THERE IS MORE THAN ONE RIGHT-HAND-SIDE FOR A GIVEN MATRIX, THE
C     USER SHOULD GO TO  DECOMP AND  SOLVE, FROM WHICH THIS ROUTINE WAS
C     DERIVED.
C
C     INPUT :
C           N   : ORDER OF MATRIX  A.    N>=2
C           NDIM: ROW DIM. OF ARRAY  AB  DECLARED IN CALLING PROGRAM
C           AB  : N*N+1 ARRAY CONTAINING  ( A B )  (I.E., MATRIX  A,
C                 AUGMENTED BY THE RIGHT-HAND-SIDE VECTOR,  B)
C
C     OUTPUT:
C           AB(1:N,N+1): SOLUTION VECTOR, X
C           IERR       : 0 ON RETURN MEANS  A  WAS FOUND TO BE SINGULAR
C
C     REFERENCE : FORSYTHE, MALCOLM, AND MOLER,  1972.
C     PROGRAMMER: D.A.SAUNDERS, INFORMATICS INC, 1979.
C
C----------------------------------------------------------------------
C
      IMPLICIT REAL ( A-H, O-Z )
C
      DIMENSION  AB(NDIM,1)
C
      NP1  = N+1
      NM1  = N-1
      IERR = 1
C
      DO 60 K=1,NM1
C
C  *     DETERMINE PIVOT ELEMENT:
C
         KP1 = K+1
         M = K
         DO 20 I=KP1,N
            IF ( ABS( AB(I,K) ).GT.ABS( AB(M,K) ) )  M = I
 20      CONTINUE
C
         T       = AB(M,K)
         AB(M,K) = AB(K,K)
         AB(K,K) = T
         IF ( T.EQ.0.E+0 )  GOTO 90
C
C  *     STORE THE MULTIPLIERS:
C
         T = -1.E+0/T
         DO 30 I=KP1,N
            AB(I,K) = AB(I,K)*T
 30      CONTINUE
C
C  *     APPLY THE MULTIPLIERS TO CURRENT SUBMATRIX, INCLUDING RHS:
C
         DO 50 J=KP1,NP1
            T       = AB(M,J)
            AB(M,J) = AB(K,J)
            AB(K,J) = T
            IF ( T.EQ.0.E+0 )  GOTO 50
            DO 40 I=KP1,N
               AB(I,J) = AB(I,J) + AB(I,K)*T
 40         CONTINUE
 50      CONTINUE
 60   CONTINUE
C
      IF ( AB(N,N).EQ.0.E+0 ) GO TO 90
C
C  *  BACK SUBSTITUTION:
C
      DO 80 KB=1,NM1
         KM1 = N-KB
         K = KM1+1
         T = AB(K,NP1)/AB(K,K)
         AB(K,NP1) = T
         DO 70 I=1,KM1
            AB(I,NP1) = AB(I,NP1) - AB(I,K)*T
 70      CONTINUE
 80   CONTINUE
C
      AB(1,NP1) = AB(1,NP1)/AB(1,1)
      RETURN
C
 90   CONTINUE
C  *  ERROR RETURN: MATRIX IS SINGULAR.
C
      IERR = 0
      RETURN
      END
