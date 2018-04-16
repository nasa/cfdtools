      SUBROUTINE DGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)
      INTEGER LDA,N,ML,MU,IPVT(*)
      DOUBLE PRECISION ABD(LDA,*),Z(*)
      DOUBLE PRECISION RCOND
C
C     DGBCO FACTORS A DOUBLE PRECISION BAND MATRIX BY GAUSSIAN
C     ELIMINATION AND ESTIMATES THE CONDITION OF THE MATRIX.
C
C     IF  RCOND  IS NOT NEEDED, DGBFA IS SLIGHTLY FASTER.
C     TO SOLVE  A*X = B , FOLLOW DGBCO BY DGBSL.
C     TO COMPUTE  INVERSE(A)*C , FOLLOW DGBCO BY DGBSL.
C     TO COMPUTE  DETERMINANT(A) , FOLLOW DGBCO BY DGBDI.
C
C     ON ENTRY
C
C        ABD     DOUBLE PRECISION(LDA, N)
C                CONTAINS THE MATRIX IN BAND STORAGE.  THE COLUMNS
C                OF THE MATRIX ARE STORED IN THE COLUMNS OF  ABD  AND
C                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS
C                ML+1 THROUGH 2*ML+MU+1 OF  ABD .
C                SEE THE COMMENTS BELOW FOR DETAILS.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  ABD .
C                LDA MUST BE .GE. 2*ML + MU + 1 .
C
C        N       INTEGER
C                THE ORDER OF THE ORIGINAL MATRIX.
C
C        ML      INTEGER
C                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.
C                0 .LE. ML .LT. N .
C
C        MU      INTEGER
C                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
C                0 .LE. MU .LT. N .
C                MORE EFFICIENT IF  ML .LE. MU .
C
C     ON RETURN
C
C        ABD     AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND
C                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        RCOND   DOUBLE PRECISION
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
C                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
C                           1.0 + RCOND .EQ. 1.0
C                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
C                UNDERFLOWS.
C
C        Z       DOUBLE PRECISION(N)
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
C                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C
C     BAND STORAGE
C
C           IF  A  IS A BAND MATRIX, THE FOLLOWING PROGRAM SEGMENT
C           WILL SET UP THE INPUT.
C
C                   ML = (BAND WIDTH BELOW THE DIAGONAL)
C                   MU = (BAND WIDTH ABOVE THE DIAGONAL)
C                   M = ML + MU + 1
C                   DO 20 J = 1, N
C                      I1 = MAX0(1, J-MU)
C                      I2 = MIN0(N, J+ML)
C                      DO 10 I = I1, I2
C                         K = I - J + M
C                         ABD(K,J) = A(I,J)
C                10    CONTINUE
C                20 CONTINUE
C
C           THIS USES ROWS  ML+1  THROUGH  2*ML+MU+1  OF  ABD .
C           IN ADDITION, THE FIRST  ML  ROWS IN  ABD  ARE USED FOR
C           ELEMENTS GENERATED DURING THE TRIANGULARIZATION.
C           THE TOTAL NUMBER OF ROWS NEEDED IN  ABD  IS  2*ML+MU+1 .
C           THE  ML+MU BY ML+MU  UPPER LEFT TRIANGLE AND THE
C           ML BY ML  LOWER RIGHT TRIANGLE ARE NOT REFERENCED.
C
C     EXAMPLE..  IF THE ORIGINAL MATRIX IS
C
C           11 12 13  0  0  0
C           21 22 23 24  0  0
C            0 32 33 34 35  0
C            0  0 43 44 45 46
C            0  0  0 54 55 56
C            0  0  0  0 65 66
C
C      THEN  N = 6, ML = 1, MU = 2, LDA .GE. 5  AND ABD SHOULD CONTAIN
C
C            *  *  *  +  +  +  , * = NOT USED
C            *  * 13 24 35 46  , + = USED FOR PIVOTING
C            * 12 23 34 45 56
C           11 22 33 44 55 66
C           21 32 43 54 65  *
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     LINPACK DGBFA
C     BLAS DAXPY,DDOT,DSCAL,DASUM,DRSCL
C     FORTRAN DABS,DMAX1,MAX0,MIN0,DSIGN
C
C     INTERNAL VARIABLES
C
      DOUBLE PRECISION DDOT,EK,T,WK,WKM
      DOUBLE PRECISION ANORM,S,DASUM,SM,YNORM
      INTEGER IS,INFO,J,JU,K,KB,KP1,L,LA,LM,LZ,M,MM
C
C
C     COMPUTE 1-NORM OF A
C
      ANORM = 0.0D0
      L = ML + 1
      IS = L + MU
      DO 10 J = 1, N
         ANORM = DMAX1(ANORM,DASUM(L,ABD(IS,J),1))
         IF (IS .GT. ML + 1) IS = IS - 1
         IF (J .LE. MU) L = L + 1
         IF (J .GE. N - ML) L = L - 1
   10 CONTINUE
C
C     FACTOR
C
      CALL DGBFA(ABD,LDA,N,ML,MU,IPVT,INFO)
C
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
C     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
C     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
C     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
C     OVERFLOW.
C
C     SOLVE TRANS(U)*W = E
C
      EK = 1.0D0
      DO 20 J = 1, N
         Z(J) = 0.0D0
   20 CONTINUE
      M = ML + MU + 1
      JU = 0
      DO 100 K = 1, N
         IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,-Z(K))
         IF (DABS(EK-Z(K)) .LE. DABS(ABD(M,K))) GO TO 30
            S = DABS(ABD(M,K))/DABS(EK-Z(K))
            CALL DSCAL(N,S,Z,1)
            EK = S*EK
   30    CONTINUE
         WK = EK - Z(K)
         WKM = -EK - Z(K)
         S = DABS(WK)
         SM = DABS(WKM)
         IF (ABD(M,K) .EQ. 0.0D0) GO TO 40
            WK = WK/ABD(M,K)
            WKM = WKM/ABD(M,K)
         GO TO 50
   40    CONTINUE
            WK = 1.0D0
            WKM = 1.0D0
   50    CONTINUE
         KP1 = K + 1
         JU = MIN0(MAX0(JU,MU+IPVT(K)),N)
         MM = M
         IF (KP1 .GT. JU) GO TO 90
            DO 60 J = KP1, JU
               MM = MM - 1
               SM = SM + DABS(Z(J)+WKM*ABD(MM,J))
               Z(J) = Z(J) + WK*ABD(MM,J)
               S = S + DABS(Z(J))
   60       CONTINUE
            IF (S .GE. SM) GO TO 80
               T = WKM - WK
               WK = WKM
               MM = M
               DO 70 J = KP1, JU
                  MM = MM - 1
                  Z(J) = Z(J) + T*ABD(MM,J)
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
         Z(K) = WK
  100 CONTINUE
      S = DASUM(N,Z,1)
      CALL DRSCL(N,S,Z,1)
C
C     SOLVE TRANS(L)*Y = W
C
      DO 120 KB = 1, N
         K = N + 1 - KB
         LM = MIN0(ML,N-K)
         IF (K .LT. N) Z(K) = Z(K) + DDOT(LM,ABD(M+1,K),1,Z(K+1),1)
         IF (DABS(Z(K)) .LE. 1.0D0) GO TO 110
            S = DABS(Z(K))
            CALL DRSCL(N,S,Z,1)
  110    CONTINUE
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
  120 CONTINUE
      S = DASUM(N,Z,1)
      CALL DRSCL(N,S,Z,1)
C
      YNORM = 1.0D0
C
C     SOLVE L*V = Y
C
      DO 140 K = 1, N
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
         LM = MIN0(ML,N-K)
         IF (K .LT. N) CALL DAXPY(LM,T,ABD(M+1,K),1,Z(K+1),1)
         IF (DABS(Z(K)) .LE. 1.0D0) GO TO 130
            S = DABS(Z(K))
            CALL DRSCL(N,S,Z,1)
            YNORM = YNORM/S
  130    CONTINUE
  140 CONTINUE
      S = DASUM(N,Z,1)
      CALL DRSCL(N,S,Z,1)
      YNORM = YNORM/S
C
C     SOLVE  U*Z = W
C
      DO 160 KB = 1, N
         K = N + 1 - KB
         IF (DABS(Z(K)) .LE. DABS(ABD(M,K))) GO TO 150
            S = DABS(ABD(M,K))/DABS(Z(K))
            CALL DSCAL(N,S,Z,1)
            YNORM = S*YNORM
  150    CONTINUE
         IF (ABD(M,K) .NE. 0.0D0) Z(K) = Z(K)/ABD(M,K)
         IF (ABD(M,K) .EQ. 0.0D0) Z(K) = 1.0D0
         LM = MIN0(K,M) - 1
         LA = M - LM
         LZ = K - LM
         T = -Z(K)
         CALL DAXPY(LM,T,ABD(LA,K),1,Z(LZ),1)
  160 CONTINUE
C     MAKE ZNORM = 1.0
      S = DASUM(N,Z,1)
      CALL DRSCL(N,S,Z,1)
      YNORM = YNORM/S
C
      IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0
      RETURN
      END

