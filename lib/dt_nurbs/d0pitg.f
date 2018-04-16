      SUBROUTINE D0PITG(N,X,Y,NDIM,NDEG,NRNG,ICC,MLT,MDIM,NLHS,INDLHS,
     +                  BCLHS,NRHS,INDRHS,BCRHS,HOLD,NHOLD,C,IER)
C
C=======================================================================
C
C  PURPOSE  D0PITG COMPUTES A SPLINE INTERPOLANT OF ODD DEGREE
C           TO USER SUPPLIED DOUBLE PRECISION DATA.
C
C  METHOD   THIS ROUTINE IS DESIGNED TO COMPUTE VERY GENERAL SPLINE
C           INTERPOLANTS OF ODD DEGREE TO USER SUPPLIED DATA INCLUD-
C           ING DERIVATIVES OF ORDER UP TO THE HALF ORDER MINUS ONE
C           OF THE SPLINE.  THIS ALSO INCLUDES GENERAL BOUNDARY CON-
C           DITIONS AT THE END POINTS. SPECIFICALLY, THE USER SUP-
C           PLIES ABSCISSAE T(I) FOR I = 1,...,N AND MULTIPLICITIES
C           M(I) FOR EACH INTERNAL ( 2 .LE. I .LE. N-1 ) ABSCISSA
C           WHERE M(I) .LE. K/2 AND K IS THE ORDER ( MUST BE EVEN )
C           OF THE SPLINE. THE MULTIPLICITIES INDICATE THE DERIVA-
C           TIVE INFORMATION THAT THE USER IS TO SUPPLY AT ANY
C           POINT.  THAT IS, THE USER IS TO SUPPLY POSITION AND
C           DERIVATIVE DATA UP TO ORDER M(I) - 1 AT EACH INTERNAL
C           ABSCISSA. IN ADDITION, THE USER SUPPLIES CERTAIN LEFT
C           AND RIGHT HAND BOUNDARY CONDITIONS WHICH INVLOVE DERI-
C           VATIVES OF ORDER UP TO K/2 - 1.  D0PITG THEN COMPUTES
C           THE COEFFICIENTS OF A SPLINE OF ORDER K WITH THE DIS-
C           TINCT KNOTS BEING THE USER SUPPLIED ABSCISSA AND WITH
C           THE INTERNAL KNOTS HAVING MULTIPLICITIES AS SUPPLIED
C           BY THE USER. THE SPLINE S(T) THEN SATISFIES THE FOL-
C           LOWWING CONDITIONS:
C
C           AT EACH INTERNAL KNOT T(I) ( 2 .LE. I .LE. N-1 )
C
C           S(J)(T(I)) = Y(J)(I) FOR 0 .LE. J .LE. M(I)-1 WHERE
C           THE NUMBERS Y(J)(I) ARE THE USER SUPPLIED POSITION AND
C           DERIVATIVE.
C
C           AT THE INITIAL POINT T(1), THE SPLINE S(T) SATISFIES
C
C           S(T(1)) = Y(1)
C
C           AND THE FOLLOWING DERIVATIVE BOUNDARY CONDITIONS. IF
C           THE NUMBER OF USER SUPPLIED LEFT HAND BOUNDARY CONDITIONS
C           (NLHS) EQUALS ZERO THEN S(T) SATISFIES THE SO CALLED
C           NATURAL BOUNDARY CONDITIONS
C
C           S(J)(T(1)) = Y(J)(1) FOR K/2 .LE. J .LE. K-2.
C
C           OTHERWISE, STARTING WITH S(K-2)(T(1)), NLHS OF THESE
C           NATURAL BOUNDARY CONDITIONS ARE REPLACED BY
C
C           S(INDLHS(I))(T(1)) = Y(INDLHS(I))(1)
C
C           WHERE INDLHS(I) IS A USER SUPPLIED INDEX AND
C           1 .LE. I .LE. NLHS.
C
C           SIMILARLY, AT THE FINAL POINT T(N), S(T) SATISFIES
C
C           S(T(N)) = Y(N)
C
C           AND DERIVATIVE AND BOUNDARY CONDITIONS ANALOGOUS TO THOSE
C           STATED ABOVE WITH NRHS AND INDRHS USED IN PLACE OF NLHS
C           AND INDLHS.
C
C           IN ORDER FOR THE INTERPOLATION PROBLEM TO BE WELL POSED
C           IT IS NECESSARY THAT NLHS + NRHS + MY .GE. K/2 WHERE
C
C            MY = SUM( ABS( MLT(I) ) + 2, I = 1, N-2.
C
C           S(J)(T(1)) = Y(J)(I) IF J IS A USER SUPPLIED INDEX
C                                LESS THAN K/2-1 OR
C           S(K-1-J)(T(1)) = 0   IF J IS NOT A USER SUPPLIED INDEX.
C
C           AT THE FINAL POINT T(N), S(T) SATISFIES
C
C
C  USAGE  DOUBLE PRECISION T(N),Y(NDIM,NRNG),BCLHS(MDIM,NRNG)
C         DOUBLE PRECISION BCRHS(MDIM,NRNG),C(MC), HOLD(NHOLD)
C         INTEGER MLT(N-2),INDLHS(NLHS),INDRHS(NRHS)
C         CALL D0PITG(N,T,Y,NDIM,NDEG,NRNG,ICC,MLT,MDIM,NLHS,INDLHS,
C                     BCLHS,NRHS,INDRHS,BCRHS,HOLD,NHOLD,C,IER)
C
C  INPUT  N       NUMBER OF DATA POINTS, N .GE. 2.
C
C         T       ARRAY OF N VALUES OF THE INDEPENDENT VARIABLE T IN
C                 ASCENDING ORDER.
C
C         Y       ARRAY OF VALUES CONTAINING THE POSITION AND DERIVATIVE
C                 VALUES FOR THE NRNG DEPENDENT VARIABLES. THE STRUC-
C                 OF THE J-TH COLUMN OF Y SHOULD BE AS FOLLOWS:
C
C                 Y(1,J) = Y(0)(1)
C
C                 Y(2,J) = Y(0)(2)
C
C                 Y(3,J) = Y(1)(2)
C                      .
C                      .
C                      .
C                 Y(MLT(1)+1,J) = Y(MLT(1)-1)(2)
C                      .
C                      .
C                      .
C                 Y(MY-MLT(M),J) = Y(0)(N-1)
C                      .
C                      .
C                      .
C                 Y(MY-MLT(M)+1,J) = Y(MLT(M)-1)(N-1)
C
C                 Y(MY,J) = Y(0)(N)
C
C
C         NDIM    DIMENSION CONSTANT FOR THE Y ARRAY, NDIM .GE. MY.
C
C         NDEG    DEGREE OF SPLINE, NDEG MUST BE ODD, NDEG GE. 1.
C
C         NRNG    NUMBER OF DEPENDENT VARIABLES, NRNG .GE. 1.
C
C         ICC     INITIAL/CONTINUE CALL CODE. ICC IS USED TO COM-
C                 MUNICATE FROM THE USER TO D0PITG WHETHER THIS
C                 IS AN INITIAL OR CONTINUE CALL.  SEE USAGE RE-
C                 MARKS FOR USAGE EXPLANATION.
C
C                 ICC .LT. 0  INITIAL CALL TO D0PITG.
C
C                 ICC .GT. 0  CONTINUE CALL TO D0PITG.
C
C         MLT     ARRAY OF N-2 VALUES CONTAINING THE MULTIPLICITY OF
C                 THE INTERIOR POINTS T(2),...,T(N-1),
C                 0 .LE. ABS(MLT(I)) .LE. K/2 WHERE K = NDEG + 1.
C
C                 MLT(I) .LT. 0  THE CORRESPONDING POINT AND ITS
C                                MLT(I)-1 DERIVATIVE VALUES WILL NOT
C                                BE INTERPOLATED.
C
C                 MLT(I) .GT. 0  THE CORRESPONDING POINT ANS ITS
C                                MLT(I)-1 DERIVATIVE VALUES WILL BE
C                                INTERPOLATED.
C
C         MDIM    DIMENSION CONSTANT FOR THE BCLHS AND BCRHS ARRAYS,
C                 MDIM .GE. MAX( NLHS, NRHS, 1 ).
C
C         NLHS    NUMBER OF DERIVATIVES AT T(1) SUPPLIED BY THE USER,
C                 0 .LE. NLHS .LE. K/2-1.
C
C         INDLHS  ARRAY OF INDICES OF THE DERIVATIVES SUPPLIED AT
C                 T(1) IN ASCENDING ORDER. 1 .LE. INDLHS(I) .LE.
C                 K/2-1 FOR 1 .LE. I .LE. NLHS.
C
C         BCLHS   ARRAY OF BOUNDARY CONDITIONS SUPPLIED AT T(1). THE
C                 I-TH, J-TH ELEMENT IS THE VALUE OF THE INDLHS(I)
C                 DERIVATIVE AT T(1) FOR THE J-TH DEPENDENT VARIABLE.
C
C         NRHS    NUMBER OF DERIVATIVES AT T(N) SUPPLIED BY THE USER,
C                 0 .LE. NRHS .LE. K/2-1.
C
C         INDRHS  ARRAY OF INDICES OF THE DERIVATIVE SUPPLIED AT
C                 T(N) IN ASCENDING ORDER. 1 .LE. INDRHS(I) .LE.
C                 K/2-1 FOR 1 .LE. I .LE. NRHS.
C
C         BCRHS   ARRAY OF BOUNDARY CONDITIONS SUPPLIED AT T(N). THE
C                 I-TH, J-TH ELEMENT IS THE VALUE OF THE INDRHS(I)
C                 DERIVATIVE AT T(N) FOR THE J-TH DEPENDENT VARIABLE.
C
C  WORKING  HOLD  WORK ARRAY OF LENGTH NHOLD.  THIS ARRAY MUST BE
C  STORAGE        UNCHANGED BY THE USER BETWEEN SUBSSEQUENT CALLS
C                 TO D0PITG FOR ICC .GT. 0.
C
C         NHOLD   THE LENGTH OF ARRAY HOLD,
C                 NHOLD .GE. NCOEF*(3*K-2) + MAX( 2*K*K, 2*NCOEF ) +
C                            + 6*K + 10
C                 WHERE NCOEF = MY + K.
C
C  OUTPUT ICC     INITIAL/CONTINUE CALL CODE:
C
C                 ICC .GT. 0  NORMALRETURN, ICC HAS BEEN SET TO THE
C                             ABSOLUTE VALUE OF ICC ON INPUT.  D0PITG
C                             MAY BE CALLED AGAIN FOR A DIFFERENT SET
C                             OF VALUES FOR THE BOUNDARY CONDITIONS
C                             AND POSITION AND DERIVATIVE DATA.
C
C                 ICC .EQ. 0  AN ERROR HAS BEEN DETECTED, IER .NE. 0,
C                             CHECK VALUE OF IER.
C
C         C       ARRAY OF MC ELEMENTS CONTAINING THE INFORMATION NEEDED
C                 TO EVALUATE THE SPLINE.  SEE HSSPVL ABSTRACT FOR A
C                 DESCRIPTION OF THE CONTENTS OF C,
C                 MC = 5 + ( NRNG + 1 ) * NCOEF + K.
C
C         IER     SUCCESS/ERROR CODE.  D0PITG CALLS THE STANDARD ERROR
C                 HANDLER DTERR TO PRINT OUT ERROR AND WARNING MESSAGES
C                 FOR IER .NE. 0. FOR IER .LT. 0, D0PITG SETS C(1) =
C                 -1.0D0.
C
C                 IER =  0   SUCCESS.
C
C                 IER =  1   RESULTS COMPUTED BUT THEY ARE SENSITIVE
C                            TO INACCURACIES IN THE ENTRIES OF T AND
C                            Y.
C
C                 IER =  2   RESULTS COMPUTED BUT THEY ARE STRONGLY
C                            SENSITIVE TO INACCURACIES IN THE ENTRIES
C                            OF T AND Y.
C
C                 IER = -1   NDEG .LT. 1 OF NDEG EVEN.
C
C                 IER = -2   N .LT. 2.
C
C                 IER = -3   NHOLD TO SMALL, THE NUMBER OF
C                             WORDS NEEDED IS GIVEN BY THE PRINT-
C                            ED ERROR MESSAGE.
C
C                 IER = -4   NDIM .LT. MY OR MDIM .LT. MAX( NLHS, NRHS,
C
C                 IER = -5   T(I) NOT IN ASCENDING ORDER.
C
C                 IER = -11  INVALID VALUE FOR ICC. THIS MAY BE CAUSED
C                            BY FAILING TO RESET ICC FOLLOWING AN
C                            UNSUCCESSFUL CALL TO D0PITG OR
C                            HOLD(1) .NE. 10.0D0 ON CONTINOUS CALL.
C
C                 IER = -12  NLHS + NRHS + MY .LT. K/2.
C
C                 IER = -20  THE MATRIX REPRESENTING THE INTERPOLATION
C                            PROBLEM IS SINGULAR.
C
C                 IER = -30   NLHS .LT. 0 OR NLHS .GT. K/2-1.
C
C                 IER = -31  INDLHS(J) .LT. 0 OR INDLHS(J) .GT. K/2-1
C                            OR INDLHS IS NOT IN ASCENDING ORDER.
C
C                 IER = -32  NRHS .LT. 0 OR NRHS .GT. K/2-1.
C
C                 IER = -33  INDRHS(J) .LT. 0 OR INDRHS(J) .GT. K/2-1
C                            OR INDRHS( IS NOT IN ASCENDING ORDER.
C
C                 IER = -34  ABS( MLT(I) ) .GT. K/2 OR .EQ. 0
C
C                 IER = -52  NRNG .LT. 1.
C
C  USAGE    THE INITIAL/CONTINUE CALL CODE ALLOWS THE USER TO EFFI-
C  REMARKS  CIENTLY COMPUTE SEVERAL SPLINE INTERPOLANTS BASED ON THE
C           SAME INDEPENDENT VARIABLE AND BOUNDARY CONDITIONS.  FOR
C           THE INITIAL CALL, ICC .LT. 0, THE USER SUPPLIES ALL INPUT
C           ARGUMENTS.  FOR SUBSEQUENT CALLS, ICC .GT. 0, THE USER CAN
C           REDEFINE BCRHS, BCRHS AND Y ARGUMENTS CORRESPONDING TO
C           THE INITIAL MLT, INDLHS, INDRHS AND X ARGUMENTS.  THE
C           MLT, INDRHS, INDLHS AND X ARGUMENTS MUST BE CHANGED
C           BETWEEN SUBSEQUENT CALLS TO D0PITG.
C
C  LOWER LEVEL ROUTINES
C
C       D0PIT1
C       D0PIT2
C
C  OTHER Library ROUTINES CALLED
C
C       DTBSP1
C       DTMCON
C       DTERR
C
C=======================================================================
C
C
C=======================================================================
C         PARAMETERS
C=======================================================================
C
      INTEGER ICC, IER, INDLHS(*), INDRHS(*), MLT(*), NDEG, NDIM, NLHS
      INTEGER NRHS, NHOLD, N, NRNG, MDIM
      DOUBLE PRECISION BCLHS(MDIM,*), BCRHS(MDIM,*), C(*), X(*)
      DOUBLE PRECISION Y(NDIM,*), HOLD(*)
C
C=======================================================================
C         INTERNAL VARIABLES
C=======================================================================
C
      DOUBLE PRECISION DIGCAL, DIGMAX, DTMCON, RCOND
      CHARACTER*8 SUBNAM
      INTEGER I1, I2, I3, I4, I5, I6, I7, I8, I, ICON, IL, J
      INTEGER KHALF, KHM1, KORD, LDA, ML, MODE, M, NCOEF, NEED, NLS1
      INTEGER NM1, NM2, NRS1
      DATA SUBNAM /'D0PITG  '/
C
C=======================================================================
C    CHECK INPUT USAGE ERRORS
C=======================================================================
C
      IER = 0
      MODE = 1
      NEED = 0
      NM2 = N - 2
C
C=======================================================================
C     CHECK INITIAL/CONTINUE CALL PARAMETER
C=======================================================================
      IF( ICC .LT. 0 ) GO TO 10
      IF( ICC .GT. 0  .AND.  HOLD(1) .EQ. 10.0D0 ) GO TO 310
      IER = -11
      GO TO 340
10    IF( NDEG .GE. 1 .AND. MOD(NDEG,2) .EQ. 1 ) GO TO 30
      IER = -1
      GO TO 340
30    IF( N .GE. 2 ) GO TO 40
      IER = -2
      GO TO 340
C
C=======================================================================
C    DEFINE ORDER AND SUM THE INTERIOR KNOTS
C=======================================================================
C
40    KORD = NDEG + 1
      KHALF = KORD / 2
      KHM1 = KHALF - 1
      NLS1 = NLHS - 1
      NRS1 = NRHS - 1
      NCOEF = KORD
      IF( N .EQ. 2 ) GO TO 60
      DO 50 I = 1, NM2
          IF( MLT(I) .LE. 0 ) GO TO 50
          NCOEF = NCOEF + MLT(I)
50    CONTINUE
C
C=======================================================================
C    CHECK THAT ENOUGH HOLDING STORAGE HAS BEEN ALLOCATED
C=======================================================================
C
60    NEED = NCOEF * ( 3 * KORD - 2 ) + MAX0( 2*KORD*KORD, 2*NCOEF )
     +       + 4 * KORD + 9
      IF( NHOLD .GE. NEED ) GO TO 70
      IER = -3
      MODE = 2
      GO TO 340
C
C=======================================================================
C    CHECK THE DIMENSION CONSTANTS
C=======================================================================
C
70    IF( NDIM .GE. ( NCOEF - KORD + 2 ) .AND.
     +    MDIM .GE. MAX0( NLHS, NRHS, 1 ) ) GO TO 80
      IER = -4
      GO TO 340
C
C=======================================================================
C    CHECK THAT X(I) IS STRICTLY INCREASING
C=======================================================================
C
80    NM1 = N - 1
      DO 90 I = 1, NM1
          IF( X(I) .LT. X(I+1) ) GO TO 90
          IER = -5
          GO TO 340
90    CONTINUE
C
C=======================================================================
C    CHECK THAT LEFT HAND BOUNDARIES ARE PROPERLY DEFINED
C=======================================================================
C
      IF( NLHS .GE. 0 .AND. NDEG .GE. ( 2 * NLHS + 1 ) ) GO TO 100
      IER = -30
      GO TO 340
100   IF( NLHS .EQ. 0 ) GO TO 130
      IF( NLHS .EQ. 1 ) GO TO 120
      DO 110 I = 1, NLS1
          IF( 1 .LE. INDLHS(I) .AND. INDLHS(I) .LE. KHM1
     +       .AND. INDLHS(I) .LT. INDLHS(I+1) ) GO TO 110
          IER = -31
          GO TO 340
110   CONTINUE
120   IF( 1 .LE. INDLHS(NLHS) .AND. INDLHS(NLHS) .LE. KHM1 ) GO TO 130
      IER = -31
      GO TO 340
C
C=======================================================================
C    CHECK THAT RIGHT HAND BOUNDARIES ARE PROPERLY DEFINED
C=======================================================================
C
130   IF( NRHS .GE. 0 .AND. NDEG .GE. ( 2 * NRHS + 1 ) ) GO TO 140
      IER = -32
      GO TO 340
140   IF( NRHS .EQ. 0 ) GO TO 170
      IF( NRHS .EQ. 1 ) GO TO 160
      DO 150 I = 1, NRS1
          IF( 1 .LE. INDRHS(I) .AND. INDRHS(I) .LE. KHM1
     +       .AND. INDRHS(I) .LT. INDRHS(I+1) ) GO TO 150
          IER = -33
          GO TO 340
150   CONTINUE
160   IF( 1 .LE. INDRHS(NRHS) .AND. INDRHS(NRHS) .LE. KHM1 ) GO TO 170
      IER = -33
      GO TO 340
C
C=======================================================================
C    CHECK THE MULTIPLICITIY OF INTERIOR KNOTS LESS THAN KHALF
C=======================================================================
C
170   IF( NM2 .EQ. 0 ) GO TO 190
      DO 180 I = 1, NM2
          IF( IABS( MLT(I) ) .LE. KHALF .AND. MLT(I) .NE. 0 ) GO TO 180
          IER = -34
          GO TO 340
180   CONTINUE
C
C=======================================================================
C    CHECK THAT SUFFICIENT INTERPOLATION CONDITIONS PROVIDED
C=======================================================================
C
190   ICON = NLHS + NRHS + NCOEF - KORD + 2
      IF( ICON .GE. ( KORD / 2 ) ) GO TO 200
      IER = -12
      GO TO 340
C
C=======================================================================
C    CHECK NUMBER OF DEPENDENT VARIABLES
C=======================================================================
C
200   IF( NRNG .GE. 1 ) GO TO 210
      IER = -52
      GO TO 340
C
C=======================================================================
C    DEFINE INTEGER PARAMETERS FOR LINPACK ROUTINES
C=======================================================================
C
210   ML = KORD - 1
      LDA = 3 * ML + 1
C
C=======================================================================
C    PARTITION THE HOLDING STORAGE IN THE FOLLOWING MANNER
C
C
C         HOLD                              LENGTH
C
C               I1  PROBLEM PARAMETERS       9
C               I2  LU DECOMP           NCOEF*LDA
C               I3  BDRY INDX LEFT      KORD
C               I4  BDRY COND LEFT      KORD
C               I5  BDRY INDX RIGHT     KORD
C               I6  BDRY COND RIGHT     KORD
C               I7  DTBSP1 WORK         KORD*KORD
C               I8  B-SPLINE VALUES     KORD*KORD
C
C
C=======================================================================
C
      I1 = 1
      I2 = 9  + I1
      I3 = I2 +  LDA * NCOEF
      I4 = I3 + KORD
      I5 = I4 +  KORD
      I6 = I5 + KORD
      I7 = I6 +  KORD
      I8 = I7 +  KORD * KORD
C
C=======================================================================
C    SET THE C-VECTOR
C=======================================================================
C
      C(1) = 1.0D0
      C(2) = DBLE( FLOAT( NRNG ) )
      C(3) = DBLE( FLOAT( KORD ) )
      C(4) = DBLE( FLOAT( NCOEF ) )
      C(5) = C(3)
      IL = 5
C
C=======================================================================
C    SET KNOTS
C=======================================================================
C
      DO 220 I = 1, KORD
          IL = IL + 1
          C(IL) = X(1)
220   CONTINUE
      IF( NM2 .EQ. 0 ) GO TO 250
      DO 240 I = 1, NM2
          IF( MLT(I) .LT. 0 ) GO TO 240
          M = MLT(I)
          DO 230 J = 1, M
              IL = IL + 1
              C(IL) = X(I+1)
230       CONTINUE
240   CONTINUE
250   DO 260 I = 1, KORD
          IL = IL + 1
          C(IL) = X(N)
260   CONTINUE
C
C=======================================================================
C    CALL D0PIT1 TO DEFINE THE INTERPOLATION MATRIX AND KNOWN VECTOR
C=======================================================================
C
      CALL D0PIT1(N,X,Y,C(6),NCOEF,KORD,MLT,NLHS,INDLHS,NRHS,
     +            INDRHS,LDA,HOLD(I2),HOLD(I7),HOLD(I8),
     +            HOLD(I3),HOLD(I4),HOLD(I5),HOLD(I6))
C
C=======================================================================
C    CALL LINPACK TO FACTOR AND SOLVE
C=======================================================================
C
      I8 = I7 +  NCOEF
      CALL DGBCO(HOLD(I2),LDA,NCOEF,ML,ML,HOLD(I8),RCOND,HOLD(I7))
C
C . . . TEST FOR EXACT SINGULARITY AND ILL-CONDITIONING OF THE MATRIX
C
      IF( RCOND .NE. 0.0E0 ) GO TO 280
      IER  = -20
      MODE = 3
      GO TO 340
C
280   DIGMAX = -DLOG10( DTMCON(5) )
      DIGCAL = -DLOG10( RCOND )
      MODE   = 0
      IF( 3.0E0 * DIGCAL .LT.         DIGMAX ) GO TO 300
      IF( 3.0E0 * DIGCAL .LT. 2.0E0 * DIGMAX ) GO TO 290
C
C . . . MATRIX IS BADLY CONDITION
C
          IER = 2
          CALL DTERR(MODE,SUBNAM,IER,NEED)
          GO TO 300
C
C . . . MATRIX IS POORLY CONDITIONED
C
290       IER = 1
          CALL DTERR(MODE,SUBNAM,IER,NEED)
C
C======================================================================
C    STORE PROBLEM PARAMETERS FOR SUBSEQUENT CALL
C======================================================================
C
300   HOLD(1) = FLOAT( I2 )
      HOLD(2) = FLOAT( I3 )
      HOLD(3) = FLOAT( I4 )
      HOLD(4) = FLOAT( I5 )
      HOLD(5) = FLOAT( I6 )
      HOLD(6) = FLOAT( I8 )
      HOLD(7) = FLOAT( NCOEF )
      HOLD(8) = FLOAT( LDA )
      HOLD(9) = FLOAT( IER )
      GO TO 320
C
C=======================================================================
C
C    SUBSEQUENT CALL
C
C=======================================================================
C
310   I2    = INT( HOLD(1) )
      I3    = INT( HOLD(2) )
      I4    = INT( HOLD(3) )
      I5    = INT( HOLD(4) )
      I6    = INT( HOLD(5) )
      I8    = INT( HOLD(6) )
      NCOEF = INT( HOLD(7) )
      LDA   = INT( HOLD(8) )
      IER   = INT( HOLD(9) )
      KORD  = IDINT( C(3) )
C
C . . . REPORT WARNING ERROR FROM PREVIOUS MATRIX FACTORIZATION
C
      IF( IER .EQ. 0 ) GO TO 320
      IF( IER .LT. 0 ) GO TO 340
      MODE = 0
      CALL DTERR(MODE,SUBNAM,IER,NEED)
C
C=======================================================================
C    LOOP THROUGH DEPENDENT VARIABLES DEFINING COEFFICIENTS
C=======================================================================
C
320   IL = 6 + KORD
      DO 330 I = 1, NRNG
          IL = IL + NCOEF
          CALL D0PIT2(HOLD(I2),LDA,Y(1,I),N,NCOEF,HOLD(I8),KORD,
     +                MLT,NLHS,INDLHS,BCLHS(1,I),HOLD(I3),HOLD(I4),NRHS,
     +                INDRHS,BCRHS(1,I),HOLD(I5),HOLD(I6),C(IL))
330   CONTINUE
      ICC = IABS( ICC )
      RETURN
340   C(1) = -1.0D0
      ICC = 0
      CALL DTERR(MODE,SUBNAM,IER,NEED)
      RETURN
      END
