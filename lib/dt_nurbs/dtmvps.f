      SUBROUTINE DTMVPS (IOP, M, N, A, NROWA, X, Y, IER )
C
C     *****************************************************************
C     *****************************************************************
C     ******                                                     ******
C     ******  DTMVPS  --  COMPUTE Y = A * X  - STANDARD STORAGE  ******
C     ******                                                     ******
C     *****************************************************************
C     *****************************************************************
C     *                                                               *
C     * PURPOSE ...                                                   *
C     *                                                               *
C     *         COMPUTE THE MATRIX VECTOR PRODUCT AX WHERE A IS A     *
C     *         D.P. M BY N MATRIX AND X IS A D.P. N VECTOR. THE      *
C     *         PRODUCT MAY BE SET EQUAL TO Y OR MAY OPTIONALLY       *
C     *         BE ADDED TO OR SUBTRACTED FROM Y WHERE Y IS A D.P.    *
C     *         M VECTOR.                                             *
C     *                                                               *
C     *****************************************************************
C     *                                                               *
C     * PARAMETERS ...                                                *
C     *                                                               *
C     *         IOP      =  INTEGER                                   *
C     *                     OPTION CODE.                              *
C     *         M        =  INTEGER                                   *
C     *                     THE NUMBER OF ROWS IN THE A AND ELEMENTS  *
C     *                     IN Y.                                     *
C     *         N        =  INTEGER                                   *
C     *                     THE NUMBER OF COLUMNS IN THE A AND        *
C     *                     ELEMENTS IN X.                            *
C     *         A        =  DOUBLE PRECISION                          *
C     *                     M BY L ARRAY IN WHICH THE A MATRIX IS     *
C     *                     STORED.                                   *
C     *         NROWA    =  INTEGER                                   *
C     *                     ROW DIMENSION OF ARRAY A.                 *
C     *         X        =  DOUBLE PRECISION                          *
C     *                     N VECTOR TO BE MULTIPLIED BY A.           *
C     *         Y        =  DOUBLE PRECISION                          *
C     *                     M VECTOR WHERE RESULTS ARE STORED.        *
C     *         IER      =  INTEGER                                   *
C     *                     SUCCESS/ERROR CODE.                       *
C     *                                                               *
C     *****************************************************************
C
      INTEGER           IOP, M, N, NROWA, IER
C
      DOUBLE PRECISION  A(NROWA,*), X(*), Y(*)
C
C     *****************************************************************
C     *                                                               *
C     * LOCAL VARIABLES ...                                           *
C
      DOUBLE PRECISION  SGN, XJ, T
C
      CHARACTER*8              SUBNAM
C
      DATA SUBNAM / 'DTMVPS'  /
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
      IF ( ITRANS .EQ. 0  .AND.  NROWA .LT. M ) IER = -4
      IF ( ITRANS .EQ. 1  .AND.  NROWA .LT. N ) IER = -4
      IF ( N .LT. 0 ) IER = -3
      IF ( M .LT. 0 ) IER = -2
C
      IF ( IER .NE. 0 ) GO TO 800
C
C     **********************
C     CHECK FOR ZERO PRODUCT
C     **********************
C
      IF (M .LE. 0) GO TO 900
C
      IF (N .GT. 0) GO TO 100
C
         IF ( ICAS .NE. 1 ) GO TO 900
C
            DO 50 I = 1, M
               Y(I) = 0.0D0
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
C     FORM  Y = AX
C     ************
C
      SGN = 1.0D0
C
      IF (ICAS .EQ. 3) SGN = - SGN
C
      IF (ICAS .GT. 1) GO TO 210
C
C ...... IF FORMING Y = AX USE THE FIRST COLUMN OF A TO INITIALIZE Y.
C
          XJ = SGN * X(1)
          DO 200 I = 1, M
              Y(I) =  A( I, 1) * XJ
  200     CONTINUE
C
          IF (N .EQ. 1) GO TO 900
C
             K1  = 2
             GO TO 220
C
  210        K1 = 1
C
C ...... FORM THE REMAINDER OF THE PRODUCT.
C
  220  DO 240 J = K1, N
C
          XJ = SGN * X(J)
          DO 230 I = 1, M
             Y(I) = Y(I) + A(I,J) * XJ
  230     CONTINUE
C
  240 CONTINUE
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
      DO 330 I = 1, M
C
         T = 0.0D0
         DO 310 J = 1, N
            T = T +   ( A(J,I) ) * X(J)
  310    CONTINUE
C
         IF ( ICAS .NE. 1 ) GO TO 320
            Y(I) = T
            GO TO 330
C
  320       Y(I) = Y(I) + SGN * T
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
