C+----------------------------------------------------------------------
C
      SUBROUTINE SOLBT ( M, MDIM, N, A, B, C, Y, IP )
C
C  SOLUTION OF BLOCK TRIDIAGONAL LINEAR SYSTEM.
C  COEFFICIENT MATRIX MUST HAVE BEEN PREVIOUSLY PROCESSED BY DECBT.
C  INPUT...
C     M    = ORDER OF EACH BLOCK.
C     MDIM = DECLARED ROW DIMENSION OF ALL MULTI-DIMENSIONAL ARRAYS.
C            (3-D ARRAYS MUST BE DECLARED WITH THE VALUE OF MDIM FOR
C            BOTH THE FIRST TWO DIMENSIONS.  MDIM PERMITS USE OF THE
C            SAME ARRAYS FOR SOLVING SYSTEMS WITH DIFFERENT SIZES OF
C            BLOCK, IN THE SAME PROGRAM.)
C     N    = NUMBER OF BLOCKS IN EACH DIRECTION OF MATRIX.
C     A,B,C= MDIM*MDIM*N ARRAYS CONTAINING BLOCK BIDIAGONAL
C            DECOMPOSITION OF COEFFICIENT MATRIX FROM DECBT.
C     IP   = MDIM*N INTEGER ARRAY OF PIVOT INFORMATION FROM DECBT.
C     Y    = MDIM*N ARRAY CONTAINING THE M*N R-H-SIDE BLOCK VECTOR
C
C  OUTPUT...
C     Y    = MDIM*N ARRAY CONTAINING SOLUTION VECTOR OF LENGTH M*N.
C
C  REFERENCE...
C     LAWRENCE LIVERMORE LABORATORY (ACQUIRED BY NASA AMES, 1980)
C
C  HISTORY...
C     03/25/81    PARAMETER MDIM INTRODUCED; ADAPTED TO USE SAME
C                 SOLVE MODULE ALREADY AVAILABLE.
C
C-----------------------------------------------------------------------
C
      IMPLICIT REAL ( A-H, O-Z )
C
      DIMENSION  A(MDIM,MDIM,N), B(MDIM,MDIM,N), C(MDIM,MDIM,N),
     +           IP(MDIM,N), Y(MDIM,N)
C
      NM1 = N-1
      NM2 = N-2
C.....FORWARD SOLUTION SWEEP.....
      CALL SOLVE ( M, MDIM, A, Y, IP )
      DO 30 K = 2,NM1
         KM1 = K-1
         DO 20 I = 1,M
            DP = 0.E+0
            DO 10 J = 1,M
               DP = DP + C(I,J,K)*Y(J,KM1)
 10         CONTINUE
            Y(I,K) = Y(I,K) - DP
 20      CONTINUE
         CALL SOLVE ( M, MDIM, A(1,1,K), Y(1,K), IP(1,K) )
 30   CONTINUE
      DO 50 I = 1,M
         DP = 0.E+0
         DO 40 J = 1,M
            DP = DP + C(I,J,N)*Y(J,NM1) + B(I,J,N)*Y(J,NM2)
 40      CONTINUE
         Y(I,N) = Y(I,N) - DP
 50   CONTINUE
      CALL SOLVE ( M, MDIM, A(1,1,N), Y(1,N), IP(1,N) )
C.....BACKWARD SOLUTION SWEEP.....
      DO 80 KB = 1,NM1
         K = N - KB
         KP1 = K+1
         DO 70 I = 1,M
            DP = 0.E+0
            DO 60 J = 1,M
               DP = DP + B(I,J,K)*Y(J,KP1)
 60         CONTINUE
            Y(I,K) = Y(I,K) - DP
 70      CONTINUE
 80   CONTINUE
      DO 100 I = 1,M
         DP = 0.E+0
         DO 90 J = 1,M
            DP = DP + C(I,J,1)*Y(J,3)
 90      CONTINUE
         Y(I,1) = Y(I,1) - DP
 100  CONTINUE
      RETURN
      END
