      SUBROUTINE DTGESL ( A, NROW, N, B, WORK, NWORK, IER )
         INTEGER NROW, N, NWORK, IER, WORK(*)
         DOUBLE PRECISION    A(NROW,*), B(*)
C
C*********************************************************************
C
C  PURPOSE   DTGESL SOLVES AX=B FOR X WHERE A IS AN N BY N DOUBLE
C            GENERAL, DENSE MATRIX.  DTGESL IS A COMPANION ROUTINE TO
C            DTGELE FOR USE WHEN SEVERAL SYSTEMS OF EQUATIONS WITH
C            THE SAME COEFFICIENT MATRIX A, BUT DIFFERENT RIGHT HAND
C            SIDES, ARE TO BE SOLVED.  DTGELE IS USED TO SOLVE THE
C            FIRST SUCH SYSTEM- DTGESL CAN THEN SOLVE EACH OF THE
C            REMAINING SYSTEMS WITH MUCH LESS COST.
C
C  METHOD    IF THE FACTORS OF A AS COMPUTED BY DTGELE HAVE
C            NONZERO DIAGONAL ENTRIES THEN LINPACK SUBROUTINE
C            DGESL IS CALLED TO COMPUTE THE SOLUTION X.
C
C  USAGE     DOUBLE PRECISION  A(NROW,N), B(N), WORK(NWORK)
C            CALL DTGELE ( A, NROW, N, B, WORK, NWORK, RCOND, IER )
C                 "COMPUTE NEW RIGHT HAND SIDE IN B"
C            CALL DTGESL ( A, NROW, N, B, WORK, NWORK, IER )
C
C  INPUT     A         DOUBLY SUBSCRIPTED ARRAY WHICH CONTAINS THE
C                      FACTORS OF A AS COMPUTED BY DTGELE.  THIS
C                      ARRAY SHOULD NOT BE CHANGED BY THE USER
C                      BETWEEN THE CALL TO DTGELE AND THE USE OF
C                      DTGESL.
C
C            NROW      ROW DIMENSION OF THE ARRAY A WHICH MUST
C                      BE AT LEAST N.  (SEE INTRODUCTION TO LINEAR
C                      ALGEBRA SECTION FOR FURTHER DETAILS.)
C
C            N         ORDER OF THE MATRIX A AND THE LENGTH OF THE
C                      RIGHT HAND SIDE COLUMN VECTOR B.
C
C            B         VECTOR WHICH CONTAINS THE RIGHT HAND SIDE
C                      COLUMN VECTOR B.
C
C  WORKING   WORK      WORK VECTOR OF LENGTH NWORK.  THIS MUST
C  STORAGE             BE THE SAME WORK VECTOR AS USED BY DTGELE.
C                      THIS VECTOR SHOULD NOT BE CHANGED BY THE USER
C                      BETWEEN THE CALL TO DTGELE AND THE USE OF
C                      DTGESL.  IT CONTAINS INFORMATION REGARDING
C                      THE FACTORS OF A WHICH MUST BE PRESERVED TO
C                      SOLVE ADDITIONAL RIGHT HAND SIDES.
C
C            NWORK     THE LENGTH OF THE DOUBLE PRECISION ARRAY WORK
C                      WHICH MUST BE AT LEAST 2*N .
C
C  OUTPUT    B         IF IER = 0, THEN B HAS BEEN OVERWRITTEN WITH THE
C                      COMPUTED SOLUTION.  OTHERWISE IT IS UNCHANGED.
C
C            IER       SUCCESS/ERROR CODE WHICH COMMUNICATES TO THE
C                      USER SUCCESS, WARNINGS, OR ERRORS.  POSSIBLE
C                      RETURN VALUES ARE
C
C                      IER =  0, NORMAL RETURN
C                          =  1, SOLUTION HAS BEEN COMPUTED BUT A
C                                IS POORLY CONDITIONED
C                          =  2, SOLUTION HAS BEEN COMPUTED BUT A
C                                IS BADLY CONDITIONED
C                          = -1, N IS LESS THAN 1
C                          = -2, NROW IS LESS THAN N
C                          = -3, NWORK IS LESS THAN 2*N
C                          = -4, MATRIX IS EXACTLY SINGULAR, NO
C                                SOLUTION WAS COMPUTED
C                          = -5, INVALID ENTRY IN THE ARRAY W.
C
C
C  WRITTEN BY ROGER G. GRIMES ON SEPTEMBER 7, 1979
C ...        MODIFIED BY KAREN J. EKBLAD DECEMBER 1985
C
C*********************************************************************
C
C ... EXTERNAL SUBROUTINES CALLED ARE
C
C          DTERR, DGESL, DTJCON
C
C
C ... INTERNAL VARIABLES
C
      CHARACTER*8 NAME
      INTEGER I, IBGN, IEND, NEED, DTJCON
C
C ... INITIALIZE THE ERROR MESSAGE ARRAYS
C
      DATA NAME/'DTGESL'/
C
C
C ... TEST INITIAL PARAMETERS
C
         IER = 0
         IF (N .GE. 1) GO TO 10
            IER = -1
            CALL DTERR ( 1, NAME, IER, 1 )
            GO TO 70
C
  10     IF ( NROW .GE. N ) GO TO 20
            IER = -2
            CALL DTERR ( 1, NAME, IER, 1)
            GO TO 70
C
  20     IF ( NWORK .GE. 2*N  ) GO TO 30
            NEED = 2*N
            IER = -3
            CALL DTERR ( 2, NAME, IER, NEED )
            GO TO 70
C
C ... TEST FOR ZERO DIAGONAL ENTRIES
C
  30     CONTINUE
         DO 40 I = 1, N
            IF ( A(I,I) .NE. 0.0D0 ) GO TO 40
               IER = -4
               CALL DTERR ( 1, NAME, IER, 1 )
               GO TO 70
  40     CONTINUE
C
C ... TEST PIVOTING INFORMATION IN ARRAY WORK
C
         IBGN = 2 * DTJCON(14) * N + 1
         IEND = IBGN + N - 1
         DO 50 I = IBGN, IEND
            IPVT = IABS ( WORK(I) )
            IF ( IPVT .LT. 1 .OR. IPVT .GT. N ) THEN
               IER = -5
               CALL DTERR ( 1, NAME, IER, 1 )
               GO TO 70
            ENDIF
  50     CONTINUE
C
C ... TEST FOR ILL-CONDITIONING OF THE MATRIX
C
         IF ( WORK(1) .EQ. 2 ) THEN
            IER = 2
            CALL DTERR ( 0, NAME, IER, 1 )
C
         ELSE IF ( WORK(1) .EQ. 1) THEN
C
            IER = 1
            CALL DTERR ( 0, NAME, IER, 1 )
         ENDIF
C
C ... SOLVE THE MATRIX EQUATIONS
C
         CALL DGESL ( A, NROW, N, WORK(IBGN), B, 0 )
C
C ... END OF THE SUBROUTINE
C
  70     CONTINUE
         RETURN
      END
