      SUBROUTINE D0PR2X ( A, MDIM, M, N, INDX, WORK, NWORK, IER )
C
C----------------------------------------------------------------------
C
C ... PURPOSE  D0PRMX APPLIES THE PERMUTATION IN INDX TO THE COLUMNS
C              OF THE DOUBLE PRECISION TWO DIMENSIONAL ARRAY A.
C              COLUMN OF A(INDX(I)) IS MOVED TO COLUMN A(I).
C
C ... INPUT    MDIM     ROW DIMENSION OF THE STORAGE ARRAY A IN WHICH
C                       AN  M  BY  N  DATA ARRAY IS STORED.
C
C              M        NUMBER OF ROWS OF DATA IN THE TWO
C                       DIMENSIONAL ARRAY A.
C
C              N        NUMBER OF COLUMNS OF DATA IN THE TWO
C                       DIMENSIONAL ARRAY A.
C                       THE COLUMNS OF A ARE TO BE INTERCHANGED.
C
C              INDX     PERMUTATION ARRAY OF LENGTH N.
C
C ... WORK     WORK     A DOUBLE PRECISION WORK ARRAY OF LENGTH NWORK,
C                       TO BE USED TO TEMPORARILY HOLD A COLUMN OF TWO
C                       DIMENSIONAL ARRAY A WHILE THE COLUMNS ARE
C                       BEING INTERCHANGED.
C
C              NWORK    DIMENSION OF ARRAY WORK, WHERE
C                       NWORK IS AT LEAST AS LARGE AS M.
C
C ... INPUT/   A        DOUBLE PRECISION TWO DIMENSIONAL ARRAY,
C      OUTPUT           DIMENSIONED AT LEAST  M  BY  N , WHOSE
C                       COLUMNS ARE TO BE REORDERED.
C
C ... OUTPUT   IER      SUCCESS/ERROR FLAG.
C
C                       IER =  0    NORMAL RETURN
C                           = -1    M .LE. 0
C                           = -2    M .GT. MDIM
C                           = -3    N .LE. 0
C                           = -4    INDX(1) .LE. 0
C                           = -5    NWORK   .LT. M
C
C----------------------------------------------------------------------
C
C ... WRITTEN BY ALVIN C. MONG, BOEING COMPUTER SERVICES,
C                OCTOBER, 1987.
C
C----------------------------------------------------------------------
C
C ... GLOBAL VARIABLES
C
         INTEGER           M, N, INDX(*), MDIM, NWORK, IER
         DOUBLE PRECISION  A(MDIM,*), WORK(*)
C
C
C----------------------------------------------------------------------
C
C ... LOCAL VARIABLES
C
         CHARACTER*8       NAME
C
         INTEGER           I, II, NEXT, NOW
C
         DATA              NAME / 'D0PR2X' /
C
C
C----------------------------------------------------------------------
C
C ... CHECK INPUT
C
         IER = 0
            IF ( NWORK   .LT. M )                  IER = -5
            IF ( INDX(1) .LE. 0 )                  IER = -4
            IF ( N       .LE. 0 )                  IER = -3
            IF ( MDIM    .LT. M )                  IER = -2
            IF ( M       .LE. 0 )                  IER = -1
            IF ( IER     .EQ. 0 )                  GO TO 100
C
C ...... INPUT ERROR DETECTED.
C
               CALL DTERR ( 1, NAME, IER, 0 )
               GO TO 900
C
C ... APPLY THE PERMUTATION
C
  100    IF ( N .EQ. 1 ) GO TO 900
C
C ...... SEARCH FOR THE FIRST ENTRY NOT PERMUTED WHICH IS INDICATED BY
C        A NONNEGATIVE VALUE IN INDX.
C
            DO 400 I = 1, N
               IF ( INDX(I) .LE. 0 ) GO TO 400
C
C ...... INITIALIZE TO FOLLOW THE CURRENT CHAIN OF PERMUTATIONS
C
               NOW        =  I
               NEXT       =  INDX(NOW)
               INDX(NOW)  = -NEXT
               IF ( NOW .EQ. NEXT ) GO TO 400
               DO 150  II = 1, M
                  WORK(II)  =  A(II,NOW)
  150          CONTINUE
C
C ...... FOLLOW THE CHAIN - PERMUTE AS YOU GO UNTIL THE CHAIN ENDS
C
  200          IF ( INDX(NEXT) .LE. 0 ) GO TO 300
C
                  DO 250 II =  1, M
                     A(II,NOW) = A(II,NEXT)
  250             CONTINUE
                  NOW       =  NEXT
                  NEXT      =  INDX(NOW)
                  INDX(NOW) = -NEXT
                  GO TO 200
C
C ...... END OF THE CHAIN
C
  300          CONTINUE
               DO 350 II =  1, M
                  A(II,NOW) = WORK(II)
  350          CONTINUE
C
C ... END OF SEARCH LOOP FOR THE NEXT CHAIN
C
  400       CONTINUE
C
C----------------------------------------------------------------------
C
C ... PERMUTATION NOW FINISHED.  RESTORE THE INDX ARRAY.
C
            DO 500 I = 1, N
               INDX(I) = - INDX(I)
  500       CONTINUE
C
C----------------------------------------------------------------------
C
C ... END OF D0PR2X
C
  900    CONTINUE
         RETURN
      END
