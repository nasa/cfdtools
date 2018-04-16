      SUBROUTINE DTMMPS (IOP, M, L, N, A, NROWA, B, NROWB,
     1                   C, NROWC, IER)
C
C     *****************************************************************
C     *****************************************************************
C     ******                                                     ******
C     ******  DTMMPS  --  COMPUTE C = A * B  - STANDARD STORAGE  ******
C     ******                                                     ******
C     *****************************************************************
C     *****************************************************************
C     *                                                               *
C     * PURPOSE ...                                                   *
C     *                                                               *
C     *         COMPUTE THE MATRIX PRODUCT AB WHERE A IS A D.P.       *
C     *         M BY L MATRIX AND B IS A D.P. L BY N MATRIX. THE      *
C     *         PRODUCT MAY BE SET EQUAL TO C OR MAY OPTIONALLY       *
C     *         BE ADDED TO OR SUBTRACTED FROM C WHERE C IS A D.P.    *
C     *         M BY N MATRIX.                                        *
C     *                                                               *
C     *****************************************************************
C     *                                                               *
C     * PARAMETERS ...                                                *
C     *                                                               *
C     *         IOP      =  INTEGER                                   *
C     *                     OPTION CODE.                              *
C     *         M        =  INTEGER                                   *
C     *                     THE NUMBER OF ROWS IN THE A AND C         *
C     *                     MATRICES.                                 *
C     *         L        =  INTEGER                                   *
C     *                     THE NUMBER OF COLUMNS IN THE A MATRIX     *
C     *                     AND THE NUMBER OF ROWS IN THE B MATRIX.   *
C     *         N        =  INTEGER                                   *
C     *                     THE NUMBER OF COLUMNS IN THE B AND C      *
C     *                     MATRICES.                                 *
C     *         A        =  DOUBLE PRECISION                          *
C     *                     M BY L ARRAY IN WHICH THE A MATRIX IS     *
C     *                     STORED.                                   *
C     *         NROWA    =  INTEGER                                   *
C     *                     ROW DIMENSION OF ARRAY A.                 *
C     *         B        =  DOUBLE PRECISION                          *
C     *                     L BY N ARRAY IN WHICH THE B MATRIX IS     *
C     *                     STORED.                                   *
C     *         NROWB    =  INTEGER                                   *
C     *                     ROW DIMENSION OF ARRAY B.                 *
C     *         C        =  DOUBLE PRECISION                          *
C     *                     M BY N ARRAY IN WHICH THE C MATRIX IS     *
C     *                     STORED.                                   *
C     *         NROWC    =  INTEGER                                   *
C     *                     ROW DIMENSION OF ARRAY                  *
C     *         IER      =  INTEGER                                   *
C     *                     SUCCESS/ERROR CODE.                       *
C     *                                                               *
C     *****************************************************************
C
      INTEGER           IOP, M, L, N, NROWA, NROWB, NROWC, IER
C
      DOUBLE PRECISION  A(NROWA,*), B(NROWB,*), C(NROWC,*)
C
C     *****************************************************************
C     *                                                               *
C     * LOCAL VARIABLES ...                                           *
C
      DOUBLE PRECISION  SGN, BKJ, T
C
      CHARACTER*8              SUBNAM
C
      DATA SUBNAM / 'DTMMPS'  /
C     *                                                               *
C     *****************************************************************
C
      IER = 0
C
C     **********************
C     CHECK INPUT PARAMETERS
C     **********************
C
      ICAS = MOD ( IOP, 10 )
      ITRANS = IOP / 10
C
      IF ( ICAS   .LT. 1  .OR.  ICAS   .GT. 3 ) IER = -1
      IF ( ITRANS .NE. 0  .AND. ITRANS .NE. 1 ) IER = -1
C
      IF ( IER .NE. 0 ) GO TO 800
C
      IF ( NROWC .LT. M ) IER = -7
      IF ( NROWB .LT. L ) IER = -6
      IF ( ITRANS .EQ. 0  .AND.  NROWA .LT. M ) IER = -5
      IF ( ITRANS .EQ. 1  .AND.  NROWA .LT. L ) IER = -5
      IF ( N .LT. 0 ) IER = -4
      IF ( L .LT. 0 ) IER = -3
      IF ( M .LT. 0 ) IER = -2
C
      IF ( IER .NE. 0 ) GO TO 800
C
C     **********************
C     CHECK FOR ZERO PRODUCT
C     **********************
C
      IF (M .LE. 0 .OR. N .LE. 0) GO TO 900
C
      IF (L .GT. 0) GO TO 100
C
         IF ( ICAS .NE. 1 ) GO TO 900
C
            DO 50 J = 1, N
               DO 50 I = 1, M
                  C(I,J) = 0.0D0
   50       CONTINUE
C
            GO TO 900
C
C     ***************
C     COMPUTE PRODUCT
C     ***************
C
  100 IF ( ITRANS .EQ. 1 ) GO TO 300
C
C     ************
C     FORM  C = AB
C     ************
C
      SGN = 1.0D0
C
      IF (ICAS .EQ. 3) SGN = - SGN
C
      DO 250 J = 1, N
C
         IF (ICAS .GT. 1) GO TO 210
C
C        ... USE FIRST COLUMN OF A WHEN FORMING C = AB
C
            BKJ = SGN * B( 1, J)
            DO 200 I = 1, M
               C(I, J) =  A( I, 1) * BKJ
  200       CONTINUE
C
            IF (L .EQ. 1) GO TO 250
C
               K1  = 2
               GO TO 220
C
  210          K1 = 1
C
C ...... COMPUTE THE REMAINDER OF THE PRODUCT FOR THE J-TH COLUMN
C        OF C
C
  220    DO 240 K = K1, L
C
            BKJ = SGN * B( K, J)
            DO 230 I = 1, M
               C( I, J) = C( I, J) + A( I, K) * BKJ
  230       CONTINUE
C
  240    CONTINUE
C
  250 CONTINUE
C
      GO TO 900
C
C     *****************
C     FORM  C = (A**H)B
C     *****************
C
  300 SGN = 1.0D0
      IF ( ICAS .EQ. 3 ) SGN = - SGN
C
C ... FORM THE I-TH,J-TH ELEMENT OF
C
      DO 330 J = 1, N
         DO 330 I = 1, M
C
            T = 0.0D0
            DO 310 K = 1, L
               T = T +   ( A(K,I) ) * B(K,J)
  310       CONTINUE
C
            IF ( ICAS .NE. 1 ) GO TO 320
               C(I,J) = T
               GO TO 330
C
  320          C(I,J) = C(I,J) + SGN * T
C
  330 CONTINUE
C
      GO TO 900
C
C     ****************
C     ERROR PROCESSING
C     ****************
C
  800 CALL DTERR( 1, SUBNAM, IER, 0)
C
C     *****************
C     END OF PROCESSING
C     *****************
C
  900 RETURN
      END
