      SUBROUTINE DTPLFT(NPTS,T,Y,NDEG,ICC,IW,WHT,NKNOTS,XKNOTS,HOLD,
     +                  NHOLD,C,IFAIL,IER)
C
C=======================================================================
C
C  PURPOSE  DTPLFT COMPUTES A WEIGHTED LEAST SQUARES SPLINE FIT TO
C           USER SPECIFIED DOUBLE PRECISION DATA.
C
C  METHOD   THE USER SUPPLIES THE ABSCISSAE, ORDINATES AND WEIGHTS
C           ( T(I), Y(I), W(I) )I = 1,...,N AND A SET OF INTERNAL KNOTS
C           X(1),...,X(M) AND A DEGREE NDEG.  THIS ROUTINE THEN COMPUTES
C           THE COEFFICIENTS OF A SPLINE S(T) OF DEGREE NDEG (ORDER
C           K = NDEG + 1) WITH INTERNAL KNOTS X(1),...,X(M) WHICH MINI-
C           MIZES THE ERROR DEFINED BY
C
C                 N
C           E = SIGMA W(I)*( Y(I) - S(T(I)) )**2
C                I=1
C
C           NOTE THAT THE INTERNAL KNOTS NEED NOT BE DISTINCT.  THUS,
C           ONE CAN FIT SPLINES WITH USER SPECIFIED MULTIPLE KNOTS.
C
C  USAGE    DOUBLE PRECISION T(N),Y(N),WHT(N),XKNOTS(NKNOTS),C(MC)
C           DOUBLE PRECISION  HOLD(NHOLD)
C           CALL DTPLFT(N,T,Y,NDEG,ICC,IWT,WHT,XKNOTS,NKNOTS,HOLD,
C                       NHOLD,C,IER)
C
C  INPUT    N       NUMBER OF DATAPOINTS, N .GE. NKNOTS + K.
C
C           T       ARRAY OF N VALUES OF THE INDEPENDENT VARIABLE T
C                   IN ASCENDING ORDER.
C
C           Y       ARRAY OF N VALUES OF THE DEPENDENT VARIABLE IN
C                   CORRESPONDENCE WITH THE T ARRAY.
C
C           NDEG    DEGREE OF THE SPLINE, NDEG .GE. 0.
C
C           ICC     INITIAL/CONTINUE CALL CODE.  ICC IS USED TO
C                   COMMUNICATE FROM THE USER TO DTPLFT WHETHER
C                   THIS IS AN INITIAL OR CONTINUE CALL.  SEE
C                   USAGE REMARKS FOR USAGE EXPLANATION.
C
C                   ICC .LT. 0    INITIAL CALL TO DTPLFT.
C
C                   ICC .GT. 0    CONTINUE CALL TO DTPLFT.
C
C           IWT     WEIGHT OPTION.
C
C                   IWT .EQ. 0  DTPLFT USES WHT(I) = 1.0 FOR ALL I AND
C                               THE WHT ARRAY IS NOT USED.
C
C                   IWT .NE. 0  DTPLFT USES WEIGHTS PROVIDED IN THE
C                               WHT ARRAY.
C
C           WHT     ARRAY OF N NONNEGATIVE WEIGHTS IN CORRESPONDENCE
C                   TO THE T ARRAY, WHT(I) .GE. 0. IF IWT .EQ. 0, WHT
C                   IS NOT USED AND IT MAY BE A DUMMY ARGUMENT.  IF
C                   IWT .NE. 0, AT LEAST NKNOTS + K OF THE WHT(I) MUST
C                   BE POSITIVE, WHT(I) .GT. 0.
C
C           NKNOTS  NUMBER OF VALUES IN THE XKNOTS ARRAY, NKNOTS .GE. 0.
C
C           XKNOTS  ARRAY OF NKNOTS VALUES CONTAINING THE INTERNAL
C                   KNOTS IN ASCENDING ORDER WITH MULTIPLICITY .LE. K.
C
C  WORKING  HOLD    WORK ARRAY OF LENGTH NHOLD.  THIS ARRAY MUST NOT BE
C                   CHANGED BY THE USER BETWEEN SUBSEQUENT CALLS TO
C                   DTPLFT FOR ICC .GT. 0.  IT CONTAINS INFORMATION
C                   WHICH MUST BE PRESERVED.
C
C           NHOLD   LENGTH OF ARRAY HOLD,
C                   NHOLD .GE. NCOEF*K + MAX( NCOEF, 3*K-2) + 7
C                   WHERE NCOEF = NKNOTS + K.
C
C  OUTPUT   ICC     INITIAL/CONTINUE CALL CODE:
C
C                   ICC .GT. 0   NORMAL RETURN, ICC HAS BEEN SET TO THE
C                                ABSOLUTE VALUE OF ICC ON INPUT.  DTPLFT
C                                MAY BE CALLED AGAIN FOR A DIFFERENT
C                                SET OF DEPENDENT VARIABLE DATA.
C
C                   ICC .EQ. 0   AN ERROR HAS BEEN DETECTED, IER .NE. 0,
C                                CHECK VALUE OF IER.
C
C           C       ARRAY OF MC ELEMENTS CONTAINING THE INFORMATION NEED
C                   TO EVALUATE THE SPLINE.  SEE HSSPVL FOR A DESCRIPTIO
C                   OF THE CONTENTS OF C, MC = 5 + NCOEF + K.
C
C           IFAIL   KNOT INTERVAL WHERE INTERLACING FAILS.  USED ONLY IF
C                   IER = -35.
C
C           IER     SUCCESS/ERROR CODE. DTPLFT CALLS THE STANDARD ERROR
C                   HANDLER DTERR TO PRINT OUT ERROR AND WARNING MESSAGE
C                   FOR IER .NE. 0.  FOR IER .LT. 0 DTPLFT SETS C(1) =
C                   -1.0.
C
C                   IER = 0     SUCCESS.
C
C                   IER =  1    RESULTS COMPUTED BUT THEY ARE SENSITIVE
C                               TO INACCURACIES IN THE ENTRIES OF T, Y,
C                               AND XKNOTS.
C
C                   IER =  2    RESULTS COMPUTED BUT THEY ARE STRONGLY
C                               SENSITIVE TO INACCURACIES IN THE ENTRIES
C                               OF  T, Y, AND XKNOTS.
C
C                   IER = -1    NDEG .LT. 0.
C
C                   IER = -2    N .LT. NKNOTS + N OR THE NUMBER OF POS-
C                               ITIVE WEIGHTS LESS THAN NKNOTS + K.
C
C                   IER = -3    NHOLD TOO SMALL, THE NUMBER OF ADDITION-
C                               AL REAL WORDS NEEDED IS GIVEN BY THE
C                               PRINTED ERROR MESSAGE.
C
C                   IER = -5    T(I) NOT IN ASCENDING ORDER.
C
C                   IER = -6    NKNOTS .LT. 0.
C
C                   IER = -8    XKNOTS(I) NOT IN ASCEDNDING ORDER OR THE
C                               MULTIPLICITY OF A KNOT EXCEEDS THE ORDER
C
C                   IER = -11   INVALID VALUE OF ICC.  THIS WAS PROBABLY
C                               CAUSED BY CALLING DTPLFT AFTER AN ERROR
C                               OCCURED WITHOUT RESETTING ICC OR
C                               HOLD(1) .NE. 8.0D0 ON CONTINOUS CALL.
C
C                   IER = -20   THE MATRIX REPRESENTING THE LEAST SQUARE
C                               PROBLEM IS SINGULAR.
C
C                   IER = -35   XKNOTS FAILED INTERLACING CONDITIONS.
C
C                   IER = -36   XKNOTS(I) .LT. T(1) OR
C                               XKNOTS(I) .GT. T(N)
C
C                   IER = -37   T(K) .LT. T(1) .OR. T(K) .GT. T(NPTS)
C                               FOR K = 2,...,NPTS-1.  IF THIS ERROR
C                               OCCURES ON A SUBSEQUENT CALL THAN
C                               THE C OR T VECTORS WERE PROBABLY
C                               ALTER BETWEEN CALLS.
C
C                   IER = -100  UNEXPECTED ERROR RETURN FORM DTILCK, SEE
C                               EXPLANATION OF MODE = 5 ERROR MESSAGE
C                               FROM DTERR.
C
C
C  USAGE  THE INITIAL/CONTINUE CALL CODE ALLOWS THE USER TO EFFICIENTLY
C  REMARK COMPUTE SEVERAL SPLINE APPROXIMATIONS BASED ON THE SAME INDEP-
C         ENDENT VARIABLE VALUES AND WEIGHT CONDITIONS.  FOR THE INITIAL
C         CALL, ICC .LT. 0, THE USER SUPPLIES ALL INPUT ARGUMENTS.  FOR
C         SUBSEQUENT CALLS, ICC .GT. 0, THE USER CAN REDEFINE THE Y
C         ARGUMENT CORRESPONDING TO THE INITIAL T, IWT AND WHT ARGUMENTS
C         THE T, IWT AND WHT ARGUMENTS MUST NOT BE CHANGED BETWEEN
C         SUBSEQUENT CALLS TO DTPLFT.
C
C  LOWER LEVEL ROUTINES
C
C     DTPLF1
C     DTPLF2
C
C  OTHER ROUTINES CALLED
C
C     DTILCK
C     DTBSPL
C     DTERR
C
C=======================================================================
C
C     PARAMETERS
C=======================================================================
C
      DOUBLE PRECISION C(*), T(*), WHT(*), XKNOTS(*), Y(*), HOLD(*)
      INTEGER ICC, IER, IW, NDEG, NKNOTS, NPTS, NHOLD
C
C=======================================================================
C     INTERNAL VARIABLES
C=======================================================================
C
      DOUBLE PRECISION DIGCAL, DIGMAX, DTMCON, RCOND
      CHARACTER*8 SUBNAM
      INTEGER IFAIL, INFO, I1, I2, I3, I4, I5, I, IL1, K, KM1
      INTEGER MODE, NACT, NCOEF, NEED, NOW, NTOTAL
      DATA SUBNAM /'DTPLFT  '/
C
C=======================================================================
C     DEFINE INTEGER PARAMETERS
C=======================================================================
C
      IER = 0
      NEED = 0
      MODE = 1
      IF( ICC .LT. 0 ) GO TO 10
      IF( ICC .GT. 0  .AND.  HOLD(1) .EQ. 8.0D0 ) GO TO 190
      IER = -11
      GO TO 210
10    K = NDEG + 1
      KM1 = K - 1
      NCOEF = NKNOTS + K
      NTOTAL = NCOEF + K
C
C=======================================================================
C     ERROR CHECKS ON INTEGER INPUT PARAMETERS
C=======================================================================
C
      IF( NDEG .GE. 0 ) GO TO 20
      IER = -1
      GO TO 210
C
C . . . CHECK NUMBER OF POINTS
C
20    IF( NPTS .GE. MAX0( NCOEF, 2 ) ) GO TO 30
      IER = -2
      GO TO 210
C
C . . . CHECK NUMBER OF POSITIVE WEIGHTS
C
30    IF( IW .EQ. 0 ) GO TO 50
      NACT = 0
      DO 40 I = 1, NPTS
          IF( WHT(I) .LE. 0.0E0 ) GO TO 40
          NACT = NACT + 1
40    CONTINUE
      IF( NACT .GE. MAX0( NCOEF, 2 ) ) GO TO 50
      IER = -2
      GO TO 210
C
C . . . CHECK NHOLD
C
50    NEED = NCOEF * K + MAX0( NCOEF, 3*K-2 ) + 7
      IF( NHOLD .GE. NEED ) GO TO 60
      IER = -3
      MODE = 2
      GO TO 210
C
C . . . CHECK NKNOTS
C
60    IF( NKNOTS .GE. 0 ) GO TO 70
      IER = -6
      GO TO 210
C
C=======================================================================
C     CHECK THAT T VECTOR ARE INCREASING
C=======================================================================
C
70    DO 80 I = 2, NPTS
          IF( T(I-1) .LE. T(I) ) GO TO 80
          IER = -5
          GO TO 210
80    CONTINUE
C
C=======================================================================
C     DEFINE COEFICIENT VECTOR
C=======================================================================
C
      C(1) = 1.0D0
      C(2) = 1.0D0
      C(3) = DBLE( FLOAT(K) )
      C(4) = DBLE( FLOAT(NCOEF) )
      C(5) = C(3)
      NOW = 5
      DO 90 I = 1, K
          NOW = NOW + 1
          C(NOW) = T(1)
90    CONTINUE
      IF( NKNOTS .EQ. 0 ) GO TO 120
      DO 110 I = 1, NKNOTS
          NOW = NOW + 1
          IF( T(1) .LT. XKNOTS(I) .AND. XKNOTS(I) .LT. T(NPTS) )
     +    GO TO 100
          IER = -36
          GO TO 210
100       C(NOW) = XKNOTS(I)
110   CONTINUE
120   DO 130 I = 1, K
          NOW = NOW + 1
          C(NOW) = T(NPTS)
130   CONTINUE
C
C=======================================================================
C     CALL DTILCK TO CHECK INTERLACING
C=======================================================================
C
      CALL DTILCK(NPTS,T,K,IW,WHT,NTOTAL,C(6),IFAIL,IER)
      IF( IER .EQ. 0 ) GO TO 140
      IER = -100
      MODE = 5
      GO TO 210
140   IF( IFAIL .EQ. 0 ) GO TO 150
      IF( IFAIL .LT. (-100) ) IER = -8
      IF( IFAIL .GT. 0 ) IER = -35
      GO TO 210
C
C=======================================================================
C     ALLOCATE STORAGE
C=======================================================================
C
C          HOLD                        LENGTH
C    I1  PROBLEM PARAMETERS           7
C    I2  BANDED MATRIX                K*NCOEF
C    I3  B-SPLINE VALUES              K
C    I4  DTBSP2 WK1                   K-1
C    I5  DTBSP2 WK2                   K-1
C
C=======================================================================
C
150   I1 = 1
      I2 = 7 + I1
      I3 = I2 +  K * NCOEF
      I4 = I3 +  K
      I5 = I4 +  K - 1
      IL1 = NCOEF + K + 6
C
C=======================================================================
C     CALL DTPLF1 TO SET UP EQUATIONS AND KNOWN VECTOR
C=======================================================================
C
      CALL DTPLF1(NPTS,T,Y,IW,WHT,NCOEF,K,C(6),HOLD(I3),HOLD(I4),
     +            HOLD(I5),HOLD(I2),C(IL1))
C
C=======================================================================
C     SOLVE THE NORMAL EQUATIONS VIA LINPACK
C=======================================================================
C
      RCOND = 0.0D0
      CALL DPBCO(HOLD(I2),K,NCOEF,KM1,RCOND,HOLD(I3),INFO)
C
C . . . TEST FOR EXACT SINGULARITY AND ILL-CONDITIONING OF THE MATRIX
C
      IF( RCOND .NE. 0.0E0 ) GO TO 160
      IER  = -20
      MODE = 3
      GO TO 210
C
160   DIGMAX = -DLOG10( DTMCON(5) )
      DIGCAL = -DLOG10( RCOND )
      MODE   = 0
      IF( 3.0E0 * DIGCAL .LT.         DIGMAX ) GO TO 180
      IF( 3.0E0 * DIGCAL .LT. 2.0E0 * DIGMAX ) GO TO 170
C
C . . . MATRIX IS BADLY CONDITION
C
          IER = 2
          CALL DTERR(MODE,SUBNAM,IER,NEED)
          GO TO 180
C
C . . . MATRIX IS POORLY CONDITIONED
C
170       IER = 1
          CALL DTERR(MODE,SUBNAM,IER,NEED)
C
C . . . SOLVE THE MATRIX EQUATIONS
C
180   CALL DPBSL(HOLD(I2),K,NCOEF,KM1,C(IL1))
C
C=======================================================================
C    SAVE PROBLEM PARAMETERS FOR POSSIBLE SUBSEQUENT CALL
C=======================================================================
C
      HOLD(1) = FLOAT( I2 )
      HOLD(2) = FLOAT( I3 )
      HOLD(3) = FLOAT( I4 )
      HOLD(4) = FLOAT( I5 )
      HOLD(5) = FLOAT( IER  )
      HOLD(6) = FLOAT( NCOEF )
      HOLD(7) = FLOAT( IL1 )
      GO TO 220
C
C======================================================================
C    SUBSEQUENT CALL
C=======================================================================
C
190   I2    = INT( HOLD(1) )
      I3    = INT( HOLD(2) )
      I4    = INT( HOLD(3) )
      I5    = INT( HOLD(4) )
      IER   = INT( HOLD(5) )
      NCOEF = INT( HOLD(6) )
      IL1   = INT( HOLD(7) )
      K     = IDINT( C(3) )
      KM1   = K - 1
C
C . . . REPORT WARNING ERROR FROM PREVIOUS MATRIX FACTORIZATION
C
      IF( IER .EQ. 0 ) GO TO 200
      IF( IER .LT. 0 ) GO TO 210
      MODE = 0
      CALL DTERR(MODE,SUBNAM,IER,NEED)
C
C=======================================================================
C    DEFINE RIGHT HAND SIDE AND SOLVE
C=======================================================================
C
200   CALL DTPLF2(NPTS,T,Y,IW,WHT,NCOEF,K,C(6),HOLD(I3),HOLD(I4),
     +            HOLD(I5),C(IL1),IER)
      IF( IER .LT. 0 ) GO TO 210
      CALL DPBSL(HOLD(I2),K,NCOEF,KM1,C(IL1))
      GO TO 220
C
C=======================================================================
C     REPORT ERRORS
C=======================================================================
C
210   C(1) = -1.0D0
      ICC  = 0
      CALL DTERR(MODE,SUBNAM,IER,NEED)
      GO TO 230
220   ICC = IABS( ICC )
230   RETURN
      END
