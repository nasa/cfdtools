C+----------------------------------------------------------------------
C
      SUBROUTINE SVD(NM,M,N,A,W,MATU,U,MATV,V,IERR,RV1)
C
      INTEGER I,J,K,L,M,N,II,I1,KK,K1,LL,L1,MN,NM,ITS,IERR
      REAL A(NM,N),W(N),U(NM,N),V(NM,N),RV1(N)
      REAL C,F,G,H,S,X,Y,Z,SCALE,ANORM
      REAL SQRT,MAX,ABS,SIGN
      LOGICAL MATU,MATV
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE SVD,
C     NUM. MATH. 14, 403-420(1970) BY GOLUB AND REINSCH.
C     HANDBOOK FOR AUTO. COMP., VOL II-LINEAR ALGEBRA, 134-151(1971).
C
C     THIS SUBROUTINE DETERMINES THE SINGULAR VALUE DECOMPOSITION
C          T
C     A=USV  OF A REAL M BY N RECTANGULAR MATRIX.  HOUSEHOLDER
C     BIDIAGONALIZATION AND A VARIANT OF THE QR ALGORITHM ARE USED.
C
C     ON INPUT.
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.  NOTE THAT NM MUST BE AT LEAST
C          AS LARGE AS THE MAXIMUM OF M AND N.
C
C        M IS THE NUMBER OF ROWS OF A (AND U).
C
C        N IS THE NUMBER OF COLUMNS OF A (AND U) AND THE ORDER OF V.
C
C        A CONTAINS THE RECTANGULAR INPUT MATRIX TO BE DECOMPOSED.
C
C        MATU SHOULD BE SET TO .TRUE. IF THE U MATRIX IN THE
C          DECOMPOSITION IS DESIRED, AND TO .FALSE. OTHERWISE.
C
C        MATV SHOULD BE SET TO .TRUE. IF THE V MATRIX IN THE
C          DECOMPOSITION IS DESIRED, AND TO .FALSE. OTHERWISE.
C
C     ON OUTPUT.
C
C        A IS UNALTERED (UNLESS OVERWRITTEN BY U OR V).
C
C        W CONTAINS THE N (NON-NEGATIVE) SINGULAR VALUES OF A (THE
C          DIAGONAL ELEMENTS OF S).  THEY ARE UNORDERED.  IF AN
C          ERROR EXIT IS MADE, THE SINGULAR VALUES SHOULD BE CORRECT
C          FOR INDICES IERR+1,IERR+2,...,N.
C
C        U CONTAINS THE MATRIX U (ORTHOGONAL COLUMN VECTORS) OF THE
C          DECOMPOSITION IF MATU HAS BEEN SET TO .TRUE.  OTHERWISE
C          U IS USED AS A TEMPORARY ARRAY.  U MAY COINCIDE WITH A.
C          IF AN ERROR EXIT IS MADE, THE COLUMNS OF U CORRESPONDING
C          TO INDICES OF CORRECT SINGULAR VALUES SHOULD BE CORRECT.
C
C        V CONTAINS THE MATRIX V (ORTHOGONAL) OF THE DECOMPOSITION IF
C          MATV HAS BEEN SET TO .TRUE.  OTHERWISE V IS NOT REFERENCED.
C          V MAY ALSO COINCIDE WITH A IF U IS NOT NEEDED.  IF AN ERROR
C          EXIT IS MADE, THE COLUMNS OF V CORRESPONDING TO INDICES OF
C          CORRECT SINGULAR VALUES SHOULD BE CORRECT.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          K          IF THE K-TH SINGULAR VALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C        RV1 IS A TEMPORARY STORAGE ARRAY.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     MODIFIED TO ELIMINATE MACHEP
C
C-----------------------------------------------------------------------
C
      IERR = 0
C
      DO 100 I = 1, M
C
         DO 100 J = 1, N
            U(I,J) = A(I,J)
  100 CONTINUE
C     .......... HOUSEHOLDER REDUCTION TO BIDIAGONAL FORM ..........
      G = 0.0E+0
      SCALE = 0.0E+0
      ANORM = 0.0E+0
C
      DO 300 I = 1, N
         L = I + 1
         RV1(I) = SCALE * G
         G = 0.0E+0
         S = 0.0E+0
         SCALE = 0.0E+0
         IF (I .GT. M) GO TO 210
C
         DO 120 K = I, M
  120    SCALE = SCALE + ABS(U(K,I))
C
         IF (SCALE .EQ. 0.0E+0) GO TO 210
C
         DO 130 K = I, M
            U(K,I) = U(K,I) / SCALE
            S = S + U(K,I)**2
  130    CONTINUE
C
         F = U(I,I)
         G = -SIGN(SQRT(S),F)
         H = F * G - S
         U(I,I) = F - G
         IF (I .EQ. N) GO TO 190
C
         DO 150 J = L, N
            S = 0.0E+0
C
            DO 140 K = I, M
  140       S = S + U(K,I) * U(K,J)
C
            F = S / H
C
            DO 150 K = I, M
               U(K,J) = U(K,J) + F * U(K,I)
  150    CONTINUE
C
  190    DO 200 K = I, M
  200    U(K,I) = SCALE * U(K,I)
C
  210    W(I) = SCALE * G
         G = 0.0E+0
         S = 0.0E+0
         SCALE = 0.0E+0
         IF (I .GT. M .OR. I .EQ. N) GO TO 290
C
         DO 220 K = L, N
  220    SCALE = SCALE + ABS(U(I,K))
C
         IF (SCALE .EQ. 0.0E+0) GO TO 290
C
         DO 230 K = L, N
            U(I,K) = U(I,K) / SCALE
            S = S + U(I,K)**2
  230    CONTINUE
C
         F = U(I,L)
         G = -SIGN(SQRT(S),F)
         H = F * G - S
         U(I,L) = F - G
C
         DO 240 K = L, N
  240    RV1(K) = U(I,K) / H
C
         IF (I .EQ. M) GO TO 270
C
         DO 260 J = L, M
            S = 0.0E+0
C
            DO 250 K = L, N
  250       S = S + U(J,K) * U(I,K)
C
            DO 260 K = L, N
               U(J,K) = U(J,K) + S * RV1(K)
  260    CONTINUE
C
  270    DO 280 K = L, N
  280    U(I,K) = SCALE * U(I,K)
C
  290    ANORM = MAX(ANORM,ABS(W(I))+ABS(RV1(I)))
  300 CONTINUE
C     .......... ACCUMULATION OF RIGHT-HAND TRANSFORMATIONS ..........
      IF (.NOT. MATV) GO TO 410
C     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
      DO 400 II = 1, N
         I = N + 1 - II
         IF (I .EQ. N) GO TO 390
         IF (G .EQ. 0.0E+0) GO TO 360
C
         DO 320 J = L, N
C     .......... DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ..........
  320    V(J,I) = (U(I,J) / U(I,L)) / G
C
         DO 350 J = L, N
            S = 0.0E+0
C
            DO 340 K = L, N
  340       S = S + U(I,K) * V(K,J)
C
            DO 350 K = L, N
               V(K,J) = V(K,J) + S * V(K,I)
  350    CONTINUE
C
  360    DO 380 J = L, N
            V(I,J) = 0.0E+0
            V(J,I) = 0.0E+0
  380    CONTINUE
C
  390    V(I,I) = 1.0E+0
         G = RV1(I)
         L = I
  400 CONTINUE
C     .......... ACCUMULATION OF LEFT-HAND TRANSFORMATIONS ..........
  410 IF (.NOT. MATU) GO TO 510
C     ..........FOR I=MIN(M,N) STEP -1 UNTIL 1 DO -- ..........
      MN = N
      IF (M .LT. N) MN = M
C
      DO 500 II = 1, MN
         I = MN + 1 - II
         L = I + 1
         G = W(I)
         IF (I .EQ. N) GO TO 430
C
         DO 420 J = L, N
  420    U(I,J) = 0.0E+0
C
  430    IF (G .EQ. 0.0E+0) GO TO 475
         IF (I .EQ. MN) GO TO 460
C
         DO 450 J = L, N
            S = 0.0E+0
C
            DO 440 K = L, M
  440       S = S + U(K,I) * U(K,J)
C     .......... DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ..........
            F = (S / U(I,I)) / G
C
            DO 450 K = I, M
               U(K,J) = U(K,J) + F * U(K,I)
  450    CONTINUE
C
  460    DO 470 J = I, M
  470    U(J,I) = U(J,I) / G
C
         GO TO 490
C
  475    DO 480 J = I, M
  480    U(J,I) = 0.0E+0
C
  490    U(I,I) = U(I,I) + 1.0E+0
  500 CONTINUE
C     .......... DIAGONALIZATION OF THE BIDIAGONAL FORM ..........
C     .......... FOR K=N STEP -1 UNTIL 1 DO -- ..........
  510 DO 700 KK = 1, N
         K1 = N - KK
         K = K1 + 1
         ITS = 0
C     .......... TEST FOR SPLITTING.
C                   FOR L=K STEP -1 UNTIL 1 DO -- ..........
  520    DO 530 LL = 1, K
            L1 = K - LL
            L = L1 + 1
            IF (ABS(RV1(L)) + ANORM .EQ. ANORM) GO TO 565
C     .......... RV1(1) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                    THROUGH THE BOTTOM OF THE  LOOP ..........
            IF (ABS(W(L1)) + ANORM .EQ. ANORM) GO TO 540
  530    CONTINUE
C     .......... CANCELLATION OF RV1(L) IF L GREATER THAN 1 ..........
  540    C = 0.0E+0
         S = 1.0E+0
C
         DO 560 I = L, K
            F = S * RV1(I)
            RV1(I) = C * RV1(I)
            IF (ABS(F) + ANORM .EQ. ANORM) GO TO 565
            G = W(I)
            H = SQRT(F*F+G*G)
            W(I) = H
            C = G / H
            S = -F / H
            IF (.NOT. MATU) GO TO 560
C
            DO 550 J = 1, M
               Y = U(J,L1)
               Z = U(J,I)
               U(J,L1) = Y * C + Z * S
               U(J,I) = -Y * S + Z * C
  550       CONTINUE
C
  560    CONTINUE
C     .......... TEST FOR CONVERGENCE ..........
  565    Z = W(K)
         IF (L .EQ. K) GO TO 650
C     .......... SHIFT FROM BOTTOM 2 BY 2 MINOR ..........
         IF (ITS .EQ. 30) GO TO 1000
         ITS = ITS + 1
         X = W(L)
         Y = W(K1)
         G = RV1(K1)
         H = RV1(K)
         F=((Y - Z) * (Y + Z) + (G - H) * (G + H)) / (2.0E+0 * H * Y)
         G = SQRT(F*F+1.0E+0)
         F = ((X - Z) * (X + Z) + H * (Y / ( F + SIGN(G,F)) - H)) / X
C     .......... NEXT QR TRANSFORMATION ..........
         C = 1.0E+0
         S = 1.0E+0
C
         DO 600 I1 = L, K1
            I = I1 + 1
            G = RV1(I)
            Y = W(I)
            H = S * G
            G = C * G
            Z = SQRT(F*F+H*H)
            RV1(I1) = Z
            C = F / Z
            S = H / Z
            F = X * C + G * S
            G = -X * S + G * C
            H = Y * S
            Y = Y * C
            IF (.NOT. MATV) GO TO 575
C
            DO 570 J = 1, N
               X = V(J,I1)
               Z = V(J,I)
               V(J,I1) = X * C + Z * S
               V(J,I) = -X * S + Z * C
  570       CONTINUE
C
  575       Z = SQRT(F*F+H*H)
            W(I1) = Z
C     .......... ROTATION CAN BE ARBITRARY IF Z IS ZERO ..........
            IF (Z .EQ. 0.0E+0) GO TO 580
            C = F / Z
            S = H / Z
  580       F = C * G + S * Y
            X = -S * G + C * Y
            IF (.NOT. MATU) GO TO 600
C
            DO 590 J = 1, M
               Y = U(J,I1)
               Z = U(J,I)
               U(J,I1) = Y * C + Z * S
               U(J,I) = -Y * S + Z * C
  590       CONTINUE
C
  600    CONTINUE
C
         RV1(L) = 0.0E+0
         RV1(K) = F
         W(K) = X
         GO TO 520
C     .......... CONVERGENCE ..........
  650    IF (Z .GE. 0.0E+0) GO TO 700
C     .......... W(K) IS MADE NON-NEGATIVE ..........
         W(K) = -Z
         IF (.NOT. MATV) GO TO 700
C
         DO 690 J = 1, N
  690    V(J,K) = -V(J,K)
C
  700 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO A
C                   SINGULAR VALUE AFTER 30 ITERATIONS ..........
 1000 IERR = K
 1001 RETURN
      END
