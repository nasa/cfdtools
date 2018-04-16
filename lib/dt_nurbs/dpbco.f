      SUBROUTINE DPBCO(ABD,LDA,N,M,RCOND,Z,INFO)
      INTEGER LDA,N,M,INFO
      DOUBLE PRECISION ABD(LDA,*),Z(*)
      DOUBLE PRECISION RCOND
C
C     DPBCO FACTORS A DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE
C     MATRIX STORED IN BAND FORM AND ESTIMATES THE CONDITION OF THE
C     MATRIX.
C
C     IF  RCOND  IS NOT NEEDED, DPBFA IS SLIGHTLY FASTER.
C     TO SOLVE  A*X = B , FOLLOW DPBCO BY DPBSL.
C     TO COMPUTE  INVERSE(A)*C , FOLLOW DPBCO BY DPBSL.
C     TO COMPUTE  DETERMINANT(A) , FOLLOW DPBCO BY DPBDI.
C
C     ON ENTRY
C
C        ABD     DOUBLE PRECISION(LDA, N)
C                THE MATRIX TO BE FACTORED.  THE COLUMNS OF THE UPPER
C                TRIANGLE ARE STORED IN THE COLUMNS OF ABD AND THE
C                DIAGONALS OF THE UPPER TRIANGLE ARE STORED IN THE
C                ROWS OF ABD .  SEE THE COMMENTS BELOW FOR DETAILS.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  ABD .
C                LDA MUST BE .GE. M + 1 .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        M       INTEGER
C                THE NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
C                0 .LE. M .LT. N .
C
C     ON RETURN
C
C        ABD     AN UPPER TRIANGULAR MATRIX  R , STORED IN BAND
C                FORM, SO THAT  A = TRANS(R)*R .
C                IF  INFO .NE. 0 , THE FACTORIZATION IS NOT COMPLETE.
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
C                UNDERFLOWS.  IF INFO .NE. 0 , RCOND IS UNCHANGED.
C
C        Z       DOUBLE PRECISION(N)
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
C                IF  A  IS SINGULAR TO WORKING PRECISION, THEN  Z  IS
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C                IF  INFO .NE. 0 , Z  IS UNCHANGED.
C
C        INFO    INTEGER
C                = 0  FOR NORMAL RETURN.
C                = K  SIGNALS AN ERROR CONDITION.  THE LEADING MINOR
C                     OF ORDER  K  IS NOT POSITIVE DEFINITE.
C
C     BAND STORAGE
C
C           IF  A  IS A SYMMETRIC POSITIVE DEFINITE BAND MATRIX,
C           THE FOLLOWING PROGRAM SEGMENT WILL SET UP THE INPUT.
C
C                   M = (BAND WIDTH ABOVE DIAGONAL)
C                   DO 20 J = 1, N
C                      I1 = MAX0(1, J-M)
C                      DO 10 I = I1, J
C                         K = I-J+M+1
C                         ABD(K,J) = A(I,J)
C                10    CONTINUE
C                20 CONTINUE
C
C           THIS USES  M + 1  ROWS OF  A , EXCEPT FOR THE  M BY M
C           UPPER LEFT TRIANGLE, WHICH IS IGNORED.
C
C     EXAMPLE..  IF THE ORIGINAL MATRIX IS
C
C           11 12 13  0  0  0
C           12 22 23 24  0  0
C           13 23 33 34 35  0
C            0 24 34 44 45 46
C            0  0 35 45 55 56
C            0  0  0 46 56 66
C
C     THEN  N = 6 , M = 2  AND  ABD  SHOULD CONTAIN
C
C            *  * 13 24 35 46
C            * 12 23 34 45 56
C           11 22 33 44 55 66
C
C     LINPACK.  THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     LINPACK DPBFA
C     BLAS DAXPY,DDOT,DSCAL,DASUM,DRSCL
C     FORTRAN DABS,DMAX1,MAX0,MIN0,DREAL,DSIGN
C
C     INTERNAL VARIABLES
C
      DOUBLE PRECISION DDOT,EK,T,WK,WKM
      DOUBLE PRECISION ANORM,S,DASUM,SM,YNORM
      INTEGER I,J,J2,K,KB,KP1,L,LA,LB,LM,MU
C
C
C     FIND NORM OF A
C
      DO 30 J = 1, N
         L = MIN0(J,M+1)
         MU = MAX0(M+2-J,1)
         Z(J) = DASUM(L,ABD(MU,J),1)
         K = J - L
         IF (M .LT. MU) GO TO 20
         DO 10 I = MU, M
            K = K + 1
            Z(K) = Z(K) + DABS(ABD(I,J))
   10    CONTINUE
   20    CONTINUE
   30 CONTINUE
      ANORM = 0.0D0
      DO 40 J = 1, N
         ANORM = DMAX1(ANORM,Z(J))
   40 CONTINUE
C
C     FACTOR
C
      CALL DPBFA(ABD,LDA,N,M,INFO)
      IF (INFO .NE. 0) GO TO 180
C
C        RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C        ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E .
C        THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
C        GROWTH IN THE ELEMENTS OF W  WHERE  TRANS(R)*W = E .
C        THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
C
C        SOLVE TRANS(R)*W = E
C
         EK = 1.0D0
         DO 50 J = 1, N
            Z(J) = 0.0D0
   50    CONTINUE
         DO 110 K = 1, N
            IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,-Z(K))
            IF (DABS(EK-Z(K)) .LE. ABD(M+1,K)) GO TO 60
               S = ABD(M+1,K)/DABS(EK-Z(K))
               CALL DSCAL(N,S,Z,1)
               EK = S*EK
   60       CONTINUE
            WK = EK - Z(K)
            WKM = -EK - Z(K)
            S = DABS(WK)
            SM = DABS(WKM)
            WK = WK/ABD(M+1,K)
            WKM = WKM/ABD(M+1,K)
            KP1 = K + 1
            J2 = MIN0(K+M,N)
            I = M + 1
            IF (KP1 .GT. J2) GO TO 100
               DO 70 J = KP1, J2
                  I = I - 1
                  SM = SM + DABS(Z(J)+WKM*ABD(I,J))
                  Z(J) = Z(J) + WK*ABD(I,J)
                  S = S + DABS(Z(J))
   70          CONTINUE
               IF (S .GE. SM) GO TO 90
                  T = WKM - WK
                  WK = WKM
                  I = M + 1
                  DO 80 J = KP1, J2
                     I = I - 1
                     Z(J) = Z(J) + T*ABD(I,J)
   80             CONTINUE
   90          CONTINUE
  100       CONTINUE
            Z(K) = WK
  110    CONTINUE
         S = DASUM(N,Z,1)
         CALL DRSCL(N,S,Z,1)
C
C        SOLVE  R*Y = W
C
         DO 130 KB = 1, N
            K = N + 1 - KB
            IF (DABS(Z(K)) .LE. ABD(M+1,K)) GO TO 120
               S = ABD(M+1,K)/DABS(Z(K))
               CALL DSCAL(N,S,Z,1)
  120       CONTINUE
            Z(K) = Z(K)/ABD(M+1,K)
            LM = MIN0(K-1,M)
            LA = M + 1 - LM
            LB = K - LM
            T = -Z(K)
            CALL DAXPY(LM,T,ABD(LA,K),1,Z(LB),1)
  130    CONTINUE
         S = DASUM(N,Z,1)
         CALL DRSCL(N,S,Z,1)
C
         YNORM = 1.0D0
C
C        SOLVE TRANS(R)*V = Y
C
         DO 150 K = 1, N
            LM = MIN0(K-1,M)
            LA = M + 1 - LM
            LB = K - LM
            Z(K) = Z(K) - DDOT(LM,ABD(LA,K),1,Z(LB),1)
            IF (DABS(Z(K)) .LE. ABD(M+1,K)) GO TO 140
               S = ABD(M+1,K)/DABS(Z(K))
               CALL DSCAL(N,S,Z,1)
               YNORM = S*YNORM
  140       CONTINUE
            Z(K) = Z(K)/ABD(M+1,K)
  150    CONTINUE
         S = DASUM(N,Z,1)
         CALL DRSCL(N,S,Z,1)
         YNORM = YNORM/S
C
C        SOLVE  R*Z = W
C
         DO 170 KB = 1, N
            K = N + 1 - KB
            IF (DABS(Z(K)) .LE. ABD(M+1,K)) GO TO 160
               S = ABD(M+1,K)/DABS(Z(K))
               CALL DSCAL(N,S,Z,1)
               YNORM = S*YNORM
  160       CONTINUE
            Z(K) = Z(K)/ABD(M+1,K)
            LM = MIN0(K-1,M)
            LA = M + 1 - LM
            LB = K - LM
            T = -Z(K)
            CALL DAXPY(LM,T,ABD(LA,K),1,Z(LB),1)
  170    CONTINUE
C        MAKE ZNORM = 1.0
         S = DASUM(N,Z,1)
         CALL DRSCL(N,S,Z,1)
         YNORM = YNORM/S
C
         IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM
         IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0
  180 CONTINUE
      RETURN
      END
