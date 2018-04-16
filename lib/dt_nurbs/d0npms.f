      SUBROUTINE D0NPMS(N,X,Y,NDIM,NDEG,NRNG,ICC,HOLD,NHOLD,C,IER)
C
C=======================================================================
C
C  PURPOSE D0NPMS COMPUTES THE NATURAL SPLINE INTERPOLANT TO A
C          SET OF USER-SUPPLIED REAL DATA WITH SEVERAL DEPENDENT
C          VARIABLES.
C
C  METHOD  GIVEN A SET OF INDEPENDENT VARIABLE VALUES T(I),I=1,...,N
C          AND THE CORRESPONDING DEPENDENT VARIABLES
C          (Y(I,1),Y(I,2),...,Y(I,NRNG)),I=1,...,N D0NPMS SETS THE
C          INPUT PARAMETERS FOR HSPITG THAT DEFINES A NATURAL SPLINE
C          INTERPOLANT TO THE DATA.  SEE HSPITG ABSTRACT FOR A
C          MORE DETAILED DESCRIPTION OF THIS PROBLEM.
C
C  USAGE   DOUBLE PRECISION T(N),Y(NDIM,NRNG),C(MC), HOLD(NHOLD)
C          CALL D0NPMS(N,T,Y,NDIM,NDEG,NRNG,ICC,HOLD,NHOLD,C,IER)
C
C  INPUT   N       NUMBER OF DATA POINTS, N .GE. MAX( 2, ( NDEG + 1 )/2
C
C          T       ARRAY OF N VALUES OF THE INDEPENDENT VARIABLE T IN
C                  ASCENDING ORDER.
C
C          Y       TWO-DIMENSIONAL ARRAY CONTAINING THE VALUES OF THE
C                  NRNG DEPENDENT VARIABLES.  Y(I,J) IS THE I-TH POINT
C                  FOR THE J-TH DEPENDENT VARIABLE.
C
C          NDIM    DIMENSIONAL CONSTANT FOR THE Y ARRAY, NDIM .GE. N.
C
C          NDEG    DEGREE OF SPLINE, NDEG MUST BE ODD. NDEG GE. 1.
C
C          NRNG    NUMBER OF DEPENDENT VARIABLES, NRNG .GE. 1.
C
C          ICC     INITIAL/CONTINUE CALL CODE. ICC IS USED TO COM-
C                  MUNICATE FROM THE USER TO D0NPMS WHETHER THIS
C                  IS AN INITIAL OR CONTINUE CALL.  SEE USAGE RE-
C                  MARKS FOR USAGE EXPLANATION.
C
C                  ICC .LT. 0  INITIAL CALL TO D0NPMS.
C
C                  ICC .GT. 0  CONTINUE CALL TO D0NPMS.
C
C  WORKING  HOLD   WORK ARRAY OF LENGTH NHOLD.  THIS ARRAY MUST BE
C  STORAGE         UNCHANGED BY THE USER BETWEEN SUBSSEQUENT CALLS
C                  TO D0NPMS FOR ICC .GT. 0.
C
C          NHOLD   THE LENGTH OF ARRAY HOLD,
C                  NHOLD .GE. NCOEF*(3*K-2) + MAX( 2*K**2, 2*NCOEF )
C                             + 4*K + N + 7
C                  WHERE NCOEF = N + K - 2 AND K = NDEG + 1.
C
C  OUTPUT  ICC     INITIAL/CONTINUE CALL CODE:
C
C                  ICC .GT. 0  NORMAL RETURN, ICC HAS BEEN SET TO THE
C                              ABSOLUTE VALUE OF ICC ON INPUT.  D0NPMS
C                              MAY BE CALLED AGAIN FOR A DIFFERENT SET
C                              OF VALUES FOR THE DEPENDENT VARIABLE.
C
C                  ICC .EQ. 0  AN ERROR HAS BEEN DETECTED, IER .NE. 0,
C                              CHECK VALUE OF IER.
C
C          C       ARRAY OF MC ELEMENTS CONTAINING THE INFORMATION NEEDE
C                  TO EVALUATE THE SPLINE.  SEE HSSPVL ABSTRACT FOR A
C                  DESCRIPTION OF THE CONTENTS OF C,
C                  MC = 5 + (NRNG + 1)*NCOEF + K.
C
C          IER     SUCCESS/ERROR CODE.  D0NPMS CALLS THE STANDARD ERROR
C                  HANDLER DTERR TO PRINT OUT ERROR AND WARNING MESSAGES
C                  FOR IER .NE. 0. FOR IER .LT. 0, D0NPMS SETS C(1) =
C                  -1.0E0.
C
C                  IER =  0   SUCCESS.
C
C                  IER =  1   RESULTS COMPUTED BUT THEY ARE SENSITIVE
C                             TO INACCURACIES IN THE ENTRIES OF T AND
C                             Y.
C
C                  IER =  2   RESULTS COMPUTED BUT THEY ARE STRONGLY
C                             SENSITIVE TO INACCURACIES IN THE ENTRIES
C                             OF T AND Y.
C
C                  IER = -1   NDEG .LT. 1 OF NDEG EVEN.
C
C                  IER = -2   N .LT. MAX( 2, ( NDEG + 1 ) / 2 ).
C
C                  IER = -3   NHOLD TO SMALL, THE NUMBER OF
C                              WORDS NEEDED IS GIVEN BY THE PRINT-
C                             ED ERROR MESSAGE.
C
C                  IER = -4   NDIM .LT. N.
C
C                  IER = -5   X(I) NOT IN ASCENDING ORDER.
C
C                  IER = -11  INVALID VALUE FOR ICC. THIS MAY BE CAUSED
C                             BY FAILING TO RESET ICC FOLLOWING AN
C                             UNSUCCESSFUL CALL TO D0NPMS OR
C                             HOLD(I2) .NE. 10.0D0 ON CONTINOUS CALL
C                             WHERE I2 = N - 1 .
C
C                  IER = -52  NRNG .LT. 1.
C
C                  IER = -100 EXPECTED ERROR RETURN FORM HSPITG, SEE
C                             EXPLANATION OF MODE = 4 ERROR MESSAGE
C                             FROM DTERR.
C
C                  IER = -200 UNEXPECTED ERROR RETURN FROM HSPITG, SEE
C                             EXPLANATION OF MODE = 5 ERROR MESSAGE
C                             FROM DTERR.
C
C  USAGE    THE INITIAL/CONTINUE CALL CODE ALLOWS THE USER TO EFFI-
C  REMARKS  CIENTLY COMPUTE SEVERAL SPLINE INTERPOLANTS BASED ON THE
C           SAME INDEPENDENT VARIABLE VALUES.  FOR THE INITIAL CALL,
C           ICC .LT. 0, THE USER SUPPLIES ALL INPUT ARGUMENTS.  FOR
C           SUBSEQUENT CALLS, ICC .GT. 0, THE USER CAN REDEFINE THE
C           Y ARGUMENT CORRESPONDING TO THE INITIAL X ARGUMENT.
C
C ... MODIFIED BY: D. D. PARSONS 10/10/90
C         TO TEST FOR INSUFFICIENT NHOLD ON CONTINUOUS CASE
C
C  LOWER LEVEL ROUTINES
C
C       D0NPM1
C
C  OTHER Library ROUTINES CALLED
C
C       HSPITG
C       DTERR
C
C=======================================================================
C
C
C=======================================================================
C     PARAMETERS
C=======================================================================
C
      INTEGER ICC, IER, N, NDIM, NDEG, NHOLD, NRNG
      DOUBLE PRECISION C(*), X(*), Y(NDIM,*), HOLD(*)
C
C=======================================================================
C     INTERNAL VARIABLES
C=======================================================================
C
      CHARACTER*8 SUBNAM
      INTEGER I, I1, I2, K, MODE, NCOEF, NEED, NMIN, NM1
      DATA SUBNAM /'D0NPMS  '/
C
C=======================================================================
C     SET INTEGER PARAMETERS AND CHECK VALID INPUT
C=======================================================================
C
      IER = 0
      MODE = 1
      NEED = 0
      I2 = N - 1
      IF (I2 .GT. NHOLD) GO TO 10
      IF( ICC .LT. 0 ) GO TO 10
      IF( ICC .GT. 0  .AND.  HOLD(I2) .EQ. 10.0D0 ) GO TO 70
      IER = -11
      GO TO 90
10    K = NDEG + 1
      NM1 = N - 1
      NCOEF = N + K - 2
      NMIN = ( NDEG + 1 ) / 2
      IF( NDEG .GE. 1 .AND. ( MOD(NDEG,2) .EQ. 1 ) ) GO TO 20
      IER = -1
      GO TO 90
20    IF( N .GE. MAX0( NMIN, 2 ) ) GO TO 30
      IER = -2
      GO TO 90
C
C=======================================================================
C    CHECK NHOLD
C=======================================================================
C
30    NEED =  NCOEF * ( 3 * K - 2 ) + MAX0( 2*K*K, 2*NCOEF )
     +        + 4 * K + N + 7
      IF( NHOLD .GE. NEED ) GO TO 40
      IER = -3
      MODE = 2
      GO TO 90
C
C=======================================================================
C    CHECK VALID DIMENSIONAL CONSTANT
C=======================================================================
C
40    IF( NDIM .GE. N ) GO TO 50
      IER = -4
      GO TO 90
C
C=======================================================================
C    CHECK THAT X(I) ARE STRICTLY INCREASING
C=======================================================================
C
50    DO 60 I = 1, NM1
            IF( X(I) .LT. X(I+1) ) GO TO 60
            IER = -5
            GO TO 90
60    CONTINUE
C
C=======================================================================
C    CHECK VALID NUMBER OF DEPENDENT VARIABLES
C=======================================================================
C
      IF( NRNG .GE. 1 ) GO TO 70
      IER = -52
      GO TO 90
C
C=======================================================================
C     CALL D0NPM1 TO SET MLT AND THEN CALL D0PITG
C
C
C     HOLD                    LENGTH
C
C      I1  MULTIPLICITIES     N-2
C      I2  D0PITG HOLD        NCOEF*(3*K-2) + MAX(2K**2, 2*NCOEF)
C                             + 4*K + 9
C=======================================================================
C
70    I1 = 1
      I2 = I1 + N - 2
      NEED = NHOLD - N + 2
      CALL D0NPM1(N,X,Y,NDIM,NDEG,NRNG,ICC,HOLD(I1),HOLD(I2),NEED,
     +            C,IER)
      IF( IER .EQ. 0 ) GO TO 100
C
C=======================================================================
C    REPORT TRACEBACK
C=======================================================================
C
C
      IF( -20 .EQ. IER .OR. IER .GT. 0 ) GO TO 80
C
C . . . REPORT UNEXPECTED TRACEBACK
C
      MODE = 5
      IER = -200
      CALL DTERR(MODE,SUBNAM,IER,NEED)
      RETURN
C
C . . . REPORT EXPECTED TRACEBACK
C
80    MODE = 4
      IF( -20 .EQ. IER ) IER = -100
      CALL DTERR(MODE,SUBNAM,IER,NEED)
      RETURN
C
C=======================================================================
C    CLOBBER OUTPUT AND RETURN
C=======================================================================
C
90    C(1) = -1.0E0
      ICC  = 0
      CALL DTERR(MODE,SUBNAM,IER,NEED)
100   RETURN
      END
