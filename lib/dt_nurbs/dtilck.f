      SUBROUTINE DTILCK(NPTS,T,KORD,IWT,WHT,NKNOTS,XKNOTS,IFAIL,IER)
C
C=======================================================================
C
C  PURPOSE  DTILCK WILL CHECK IF A SUBSEQUENCE OF A GIVEN INDEP-
C           PENDENT VARIABLE SATISFIES THE INTERLACING CONDITIONS
C           WITH RESPECT TO A GIVEN KNOT SET AND THAT THE KNOT
C           SET IS IN ASCENDING ORDER WITH ALL KNOTS HAVING A
C           MULTIPLICITY LESS THAN OR EQUAL TO THE THE ORDER.
C
C  USAGE    REAL XKNOTS(NKNOTS),T(N),WHT(N)
C           CALL DTILCK(N,T,KORD,IWT,WHT,NKNOTS,XKNOTS,IFAIL,IER)
C
C  INPUT    N       NUMBER OF VALUES IN THE T ARRAY,
C                   N .GE. NKNOTS - KORD.
C
C           T       ARRAY OF N VALUES OF THE INDEPENDENT VARIABLE
C                   IN ASCENDING ORDER.
C
C           KORD    ORDER OF THE SPLINE, KORD .GE. 1.
C
C           IWT     WEIGHT OPTION.
C
C                   IWT = 0         ALL WEIGHTS ASSUMED TO BE 1.
C                                   WHT ARRAY NOT USED.
C                   IWT .NE. 0      WEIGHTS SUPPLIED IN WHT ARRAY.
C
C           WHT     WEIGHT ARRAY.  USED ONLY IF IWT .NE. 0, AND THEN,
C                   ONLY TO DETERMINE WHICH WEIGHTS ARE POSITIVE.
C
C           NKNOTS  NUMBER OF VALUES IN THE XKNOTS ARRAY,
C                   NKNOTS .GE. 2*KORD.
C
C           XKNOTS  ARRAY OF NKNOTS VALUES CONTAINING THE KNOTS
C                   IN INCREASING ORDER WHERE THE MULTIPLICITY OF
C                   ANY ONE KNOT DOES NOT EXCEED KORD.
C
C  OUTPUT   IFAIL   INTEGER PARAMETER DESCRIBING THE RESULT OF THE
C                   INTERLACING AND KNOT VECTOR CHECK.
C
C                   IFAIL = 0       INTERLACING CONDITIONS SATISFIED.
C
C                   IFAIL = -1      T(1) .GT. XKNOTS(1)
C
C                   IFAIL = -2      T(N) .LT. XKNOTS(NKNOTS).
C
C                   IFAIL = -100-I  MULTIPLICITY OF XKNOTS(I) .GT. KORD
C                                   OR XKNOTS(I) .GT. XKNOTS(I+1).
C
C                   IFAIL = I       INTERLACING CONDITIONS FAILED AT
C                                   XKNOTS(I).
C
C           IER     SUCCESS/ ERROR CODE.  DTILCK CALLS THE STANDARD
C                   ERROR HANDLER DTERR TO PRINT OUT ERROR AND WARN-
C                   ING MESSAGES FOR IER .NE. 0.
C
C                   IER = 0   SUCCESS.
C
C                   IER = -1  KORD .LT. 1.
C
C                   IER = -2  N .GT. NKNOTS - KORD.
C
C                   IER = -5  T(I) NOT IN ASCENDING ORDER.
C
C                   IER = -6  NKNOTS .LT. 2*KORD.
C
C  OTHER ROUTINES CALLED
C
C     DTILC1
C     DTERR
C
C======================================================================
C
C=======================================================================
C     PARAMETERS
C=======================================================================
C
      EXTERNAL DTJCON
      INTEGER DTJCON
      DOUBLE PRECISION T(*), XKNOTS(*), WHT(*)
      INTEGER IER, IFAIL, KORD, NKNOTS, NPTS
C
C=======================================================================
C     INTERNAL VARIABLES
C=======================================================================
C
      CHARACTER*8 SUBNAM
      INTEGER I, IL1, IL2, IT, JNOW, MODE, NEED, NTOT
      DATA SUBNAM /'DTILCK  '/
C
C=======================================================================
C     DEFINE INTEGER PARAMETERS AND CHECK USAGE ERRORS
C=======================================================================
C
      IER = 0
      IFAIL = 0
      MODE = 1
      NEED = 0
      JNOW = 0
      JLAST = 0
      NTOT = NKNOTS - KORD
      IF( KORD .GE. 1 ) GO TO 10
      IER = -1
      GO TO 200
10    IF( NPTS .GE. MAX0( NTOT, 2 ) ) GO TO 20
      IER = -2
      GO TO 200
20    DO 30 I = 2, NPTS
          IF( T(I-1) .LE. T(I) ) GO TO 30
          IER = -5
          GO TO 200
30    CONTINUE
      IF( NKNOTS .GE. 2 * KORD ) GO TO 40
      IER = -6
      GO TO 200
C
C=======================================================================
C    CALL DTILC1 TO CHECK ASCENDING KNOT SEQUENCE WITH MULTIPLICITY
C    LESS THAN KORD
C=======================================================================
C
40    CALL DTILC1(XKNOTS,NKNOTS,KORD,IFAIL)
      IF( IFAIL .EQ. 0 ) GO TO 50
      IF( IFAIL .GT. 0 ) IFAIL = -100 - IFAIL
      GO TO 210
C
C=======================================================================
C     CHECK THAT KNOT INTERVAL CONTAINS T(1) AND T(NPTS)
C=======================================================================
C
50    IF( XKNOTS(1) .LE. T(1) ) GO TO 60
      IFAIL = -1
      GO TO 210
60    IF( T(NPTS) .LE. XKNOTS(NKNOTS) ) GO TO 70
      IFAIL = -2
      GO TO 210
C
C=======================================================================
C     CHECK FOR INTERLACING
C=======================================================================
C
70    DO 180 I = 1, NTOT
          IT = I
          IL1 = I + KORD
          IL2 = IL1 - 1
          IF( XKNOTS(I) .EQ. XKNOTS(IL2) ) GO TO 110
          IF( XKNOTS(I+1) .EQ. XKNOTS(IL1) ) GO TO 140
C
C=======================================================================
C      CASE WHEN XKNOTS(I) .LT. XKNOTS(I+K)
C=======================================================================
C
80        JNOW = JNOW + 1
          IF( JNOW .GT. NPTS ) GO TO 190
          IF( IWT .EQ. 0 ) GO TO 90
          IF( WHT(JNOW) .LE. 0.D0 ) GO TO 80
90        IF( JLAST .EQ. 0 ) GO TO 100
          IF( T(JNOW) .LE. T(JLAST) ) GO TO 80
100       IF( XKNOTS(I) .GE. T(JNOW) ) GO TO 80
          IF( XKNOTS(IL1) .LE. T(JNOW) ) GO TO 80
          GO TO 170
C
C=======================================================================
C     CASE WHEN XKNOTS(I) EQUALS XKNOTS(I+K-1)
C=======================================================================
C
110       JNOW = JNOW + 1
          IF( JNOW .GT. NPTS ) GO TO 190
          IF( IWT .EQ. 0 ) GO TO 120
          IF( WHT(JNOW) .LE. 0.D0 ) GO TO 110
120       IF( JLAST .EQ. 0 ) GO TO 130
          IF( T(JNOW) .LE. T(JLAST) ) GO TO 110
130       IF( XKNOTS(I) .GT. T(JNOW) ) GO TO 110
          IF( XKNOTS(IL1) .LE. T(JNOW) ) GO TO 110
          GO TO 170
C
C=======================================================================
C     CASE WHEN XKNOTS(I+1) EQUALS XKNOTS(I+KORD)
C=======================================================================
C
140       JNOW = JNOW + 1
          IF( JNOW .GT. NPTS ) GO TO 190
          IF( IWT .EQ. 0 ) GO TO 150
          IF( WHT(JNOW) .LE. 0.D0 ) GO TO 140
150       IF( JLAST .EQ. 0 ) GO TO 160
          IF( T(JNOW) .LE. T(JLAST) ) GO TO 140
160       IF( XKNOTS(I) .GE. T(JNOW) ) GO TO 140
          IF( XKNOTS(IL1) .LT. T(JNOW) ) GO TO 140
          GO TO 170
170       JLAST = JNOW
180   CONTINUE
C
C=======================================================================
C     REPORT ERROR AND/OR RETRUN
C=======================================================================
C
      GO TO 210
190   IFAIL = IT
      GO TO 210
200   CALL DTERR(MODE,SUBNAM,IER,NEED)
      IFAIL = DTJCON(1)
210   RETURN
      END
