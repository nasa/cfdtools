      SUBROUTINE D0PRMI ( A, N, INDX, IER )
C----------------------------------------------------------------------
C
C ... PURPOSE  D0PRMI APPLIES THE PERMUTATION IN INDX TO THE
C              INTEGER ARRAY A.  A(INDX(I)) IS MOVED
C              TO A(I).
C
C ... INPUT    N        NUMBER OF RECORDS IN A TO BE INTERCHANGED.
C
C              INDX      PERMUTATION ARRAY OF LENGTH N.
C
C ... INPUT/   A        INTEGER ARRAY TO BE REORDERED.
C      OUTPUT
C
C ... OUTPUT   IER      SUCCESS/ERROR FLAG.
C
C                       IER =  0     NORMAL RETURN
C                           = -1    N .LE. 0
C                           = -2    INDX(1) .LE. 0
C
C----------------------------------------------------------------------
C
C ... WRITTEN BY ROGER G GRIMES, BOEING COMPUTER SERVICES,
C                SEPTEMBER, 1983
C
C----------------------------------------------------------------------
C
C ... GLOBAL VARIABLES
C
         INTEGER  A(*)
C
         INTEGER           N, INDX(*), IER
C
C----------------------------------------------------------------------
C
C ... LOCAL VARIABLES
C
C
         CHARACTER*8 NAME
C
         INTEGER  ATEMP
C
         INTEGER           I, NEXT, NOW
C
C
         DATA              NAME / 'D0PRMI' /
C
C----------------------------------------------------------------------
C
C ... CHECK INPUT
C
         IER = 0
            IF ( INDX(1) .LE. 0 )                  IER = -2
            IF ( N .LE. 0 )                       IER = -1
            IF ( IER .EQ. 0 ) GO TO 100
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
               NOW         =  I
               NEXT        =  INDX(NOW)
               INDX(NOW)    = -NEXT
               IF ( NOW .EQ. NEXT ) GO TO 400
               ATEMP       =  A(NOW)
C
C ...... FOLLOW THE CHAIN - PERMUTE AS YOU GO UNTIL THE CHAIN ENDS
C
  200          IF ( INDX(NEXT) .LE. 0 ) GO TO 300
C
                  A(NOW)      =  A(NEXT)
                  NOW         =  NEXT
                  NEXT        =  INDX(NOW)
                  INDX(NOW)    = -NEXT
                  GO TO 200
C
C ...... END OF THE CHAIN
C
  300          A(NOW) = ATEMP
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
C ... END OF D0PRMI
C
  900    CONTINUE
         RETURN
      END
