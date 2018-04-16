C<*>
      SUBROUTINE DTBSPL(XKNOTS,NKNOTS,X,IL,K,NDERV,NDIM,WORK,NWORK,
     +                  BSVAL,IER)
C***********************************************************************
C
C  PURPOSE  DTBSPL EVALUATES THE FUNCTION AND DERIVATIVE VALUES
C           OF THE NORMALIZED B-SPLINES ON A GIVEN KNOT SET.
C
C  METHOD   GIVEN AN ASCENDING KNOT SET AND A POINT AT WHICH
C           B-SPLINE FUNCTION AND DERIVATIVE VALUES ARE DESIRED
C           DTBSPL USES DE BOOR'S ALGORITHM TO EFFICIENTLY
C           EVALUATE THE B-SPLINES.
C           DTBSPL IS THE DOUBLE PRECISION VERSION OF HSBSPL.
C
C  USAGE    DOUBLE PRECISION XKNOTS(NKNOTS),BSVAL(NDIM,NDERV+1)
C           DIMENSION WORK(NWORK)
C           CALL DTBSPL(XKNOTS,NKNOTS,X,IL,K,NDERV,NDIM,WORK,NWORK,
C                       BSVAL,IER)
C
C  INPUT    XKNOTS  ARRAY OF NKNOTS VALUES CONTAINING THE KNOTS
C                   IN ASCENDING ORDER.
C
C           NKNOTS  NUMBER OF VALUES IN XKNOTS ARRAY, NKNOTS .GE. 2*K.
C
C           X       VALUE OF THE VARIABLE AT WHICH B-SPLINE EVALUATION
C                   IS DESIRED, XKNOTS(K) .LE. X .LE. XKNOTS(NKNOTS-K+1)
C
C           IL      KNOT AT WHICH SEARCH IS TO BEGIN.  IF DTBSPL IS
C                   TO BE CALLED REPEATEDLY TO EVALUATE B-SPLINES AT
C                   A MONOTONE SERIES OF POINTS, THEN IL SHOULD NOT
C                   BE CHANGED BETWEEN CALLS.
C
C           K       ORDER OF THE B-SPLINES, K .GE. 1.
C
C           NDERV   NUMBER OF DERIVATIVES DESIRED, NDERV .GE. 0.
C
C           NDIM    DIMENSION CONSTANT FOR ARRAY BSVAL, NDIM .GE. NBS =
C                   NKNOTS - K.
C
C  WORKING  WORK    WORK ARRAY OF LENGTH NWORK.
C  STORAGE
C           NWORK   THE LENGTH OF ARRAY WORK, NWORK .GE. K*K.
C
C  OUTPUT   BSVAL   ARRAY OF B-SPLINE VALUES WHERE BSVAL(I,J) IS THE
C                   (J-1)ST DERIVATIVE OF THE I-TH B-SPLINE,
C                   J = 1, ..., NDERV+1, I = 1, ..., NKNOTS-K.
C
C           IL      KNOT SUCH THAT XKNOTS(IL) .LE. X .LT. XKNOTS(IL+1)
C                   (IF X = XKNOTS(NKNOTS-K+1) THEN X = XKNOTS(IL+1))
C
C           IER     SUCCESS/ERROR CODE. DTBSPL CALLS THE STANDARD ERROR
C                   HANDLER DTERR TO PRINT OUT WARNING AND ERROR MESSAGE
C                   FOR IER .NE. 0.  FOR IER .LT. 0, DTBSPL SETS
C                   BSVAL(1,1) = DTMCON(1) AND IL = -1.
C
C                   IER =  0   SUCCESS.
C
C                   IER = -1   K .LT. 1.
C
C                   IER = -3   NWORK TOO SMALL, THE NUMBER OF
C                               WORDS NEEDED IS GIVEN BY THE PRINTED
C                              ERROR MESSAGE.
C
C                   IER = -4   NDIM .LT. NKNOTS - K.
C
C                   IER = -6  NKNOTS .LT. 2*K.
C
C                   IER = -7  NDERV .LT. 0
C
C                   IER = -8  XKNOTS ARRAY NOT IN ASCENDING ORDER OR
C                             THE MULTIPLICITY OF A KNOT EXCEEDS THE
C                             ORDER.
C
C                   IER = -9  X .LT. XKNOTS(K) OR X .GT.
C                             XKNOTS(NKNOTS-K+1).
C
C  LOWER LEVEL ROUTINES
C
C        DTBSP1
C        DTBSP2
C
C  OTHER LIBRARY ROUTINES CALLED
C
C        DTILC1
C        DTMCON
C        DTERR
C
C  AUTHORS:
C
C        DAVE FERGUSON, RICH MASTRO, AND FRITZ KLEIN
C        DECEMBER, 1984
C
C        BUG FIXED BY FRITZ KLEIN -- JUNE, 1992
C***********************************************************************
C
C=======================================================================
C    PARAMETERS
C=======================================================================
C
      INTEGER IER, K, NDERV, NDIM, NKNOTS, NWORK, IL
      DOUBLE PRECISION BSVAL(NDIM,*), X, XKNOTS(*), WORK(*)
C
C=======================================================================
C    INTERNAL VARIABLES
C=======================================================================
C
      CHARACTER*8 SUBNAM
      DOUBLE PRECISION DTMCON
      INTEGER IFAIL, ISTART, LOC, MODE, NBS, NDERV1, NEED, NKNL
      DATA SUBNAM /'DTBSPL  '/
C
C=======================================================================
C    INPUT ERROR CHECKS
C=======================================================================
C
      IER = 0
      NEED = 0
      MODE = 1
      IF( K .GE. 1 ) GO TO 10
      IER = -1
      GO TO 130
10    NEED = K*K
      IF( NWORK .GE. NEED ) GO TO 20
      IER = -3
      MODE = 2
      GO TO 130
C
C=======================================================================
C    CHECK INTEGER PARAMETERS THAT DEFINE B-SPLINES
C=======================================================================
C
20    IF( NKNOTS .GE. 2 * K ) GO TO 30
      IER = -6
      GO TO 130
30    NBS = NKNOTS - K
      IF( NDIM .GE. NBS ) GO TO 40
      IER = -4
      GO TO 130
40    IF( 0 .LE. NDERV ) GO TO 50
      IER = -7
      GO TO 130
C
C=======================================================================
C    CHECK KNOT VECTOR INCREASING AND VALID X VALUE
C=======================================================================
C
50    CALL DTILC1(XKNOTS,NKNOTS,K,IFAIL)
      IF( IFAIL .EQ. 0 ) GO TO 60
      IER = -8
      GO TO 130
60    NKNL = NKNOTS - K + 1
      IF( XKNOTS(K) .LE. X .AND. X .LE. XKNOTS(NKNL) ) GO TO 70
      IER = -9
      GO TO 130
C
70    CONTINUE
C
C=======================================================================
C    FIND IL SUCH THAT XKNOTS(IL) .LE. X .LT. XKNOTS(IL+1) UNLESS
C    X = XKNOTS(NBS+1).  IN THIS CASE, FIND IL SUCH THAT
C    XKNOTS(IL) .LT. X .LE. XKNOTS(IL+1).
C=======================================================================
C
      IL = MAX0(K, MIN0(IL, NBS))
      IF(X .LT. XKNOTS(IL)) GO TO 90
      ISTART = IL
      DO 80 I = ISTART, NBS
          IL = I
          IF(X .LT. XKNOTS(I+1)) GO TO 100
80    CONTINUE
      IL = NBS
      IF ( X .EQ. XKNOTS(NBS+1) ) THEN
85       CONTINUE
         IF ( XKNOTS(IL) .LT. XKNOTS(IL+1) ) GO TO 100
         IL = IL - 1
         GO TO 85
      ENDIF
      GO TO 100
90    IF(X .GE. XKNOTS(IL)) GO TO 100
      IL = IL - 1
      IF(IL .GT. K) GO TO 90
      IL = K
100   CONTINUE
C
C=======================================================================
C    INITIALIZE OUTPUT ARRAY TO ZERO AND CALL DTBSP1 TO COMPUTE
C    THE NONZERO B-SPLINE VALUES.
C=======================================================================
C
      NDERV1 = NDERV + 1
      DO 120 I = 1, NBS
        DO 110 J = 1, NDERV1
          BSVAL(I,J) = 0.D0
  110   CONTINUE
  120 CONTINUE
C
      NDV = MIN0 ( NDERV, K-1 )
      LOC    = 1
      ISTART = IL - K + 1
      CALL DTBSP1 (XKNOTS,X,IL,K,NDV,WORK(LOC),BSVAL(ISTART,1),NDIM)
      GO TO 140
C
C=======================================================================
C    REPORT USAGE ERROR AND CLOBBER OUTPUT
C=======================================================================
C
  130 CALL DTERR(MODE,SUBNAM,IER,NEED)
      BSVAL(1,1) = DTMCON(1)
      IL = -1
  140 RETURN
      END
