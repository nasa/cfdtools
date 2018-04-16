C***** DMAX *****
      DOUBLE PRECISION FUNCTION DMAX(N,SX,INCX)
C
C     FINDS THE ELEMENT HAVING MAX. VALUE.
C
      DOUBLE PRECISION SX(1), DTMCON
      INTEGER I,INCX,INCXA,IX,N
C
      DMAX = - DTMCON(3)
      IF( N .LT. 1 ) RETURN
      DMAX = SX(1)
      IF(N.EQ.1)RETURN
      IF(INCX.EQ.1)GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      IX = 1
      INCXA = IABS ( INCX )
      IX = IX + INCXA
      DO 10 I = 2,N
         IF(SX(IX).LE.DMAX) GO TO 5
            DMAX = SX(IX)
    5    IX = IX + INCXA
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
   20 DO 30 I = 2,N
         IF(SX(I).LE.DMAX) GO TO 30
            DMAX = SX(I)
   30 CONTINUE
      RETURN
      END
C
C
C***** DMIN *****
      DOUBLE PRECISION FUNCTION DMIN(N,SX,INCX)
C
C     FINDS THE ELEMENT HAVING MIN. VALUE.
C
      DOUBLE PRECISION SX(1),DTMCON
      INTEGER I,INCX,INCXA,IX,N
C
      DMIN = DTMCON(3)
      IF( N .LT. 1 ) RETURN
      DMIN = SX(1)
      IF(N.EQ.1)RETURN
      IF(INCX.EQ.1)GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      IX = 1
      INCXA = IABS ( INCX )
      IX = IX + INCXA
      DO 10 I = 2,N
         IF(SX(IX).GE.DMIN) GO TO 5
            DMIN = SX(IX)
    5    IX = IX + INCXA
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
   20 DO 30 I = 2,N
         IF(SX(I).GE.DMIN) GO TO 30
            DMIN = SX(I)
   30 CONTINUE
      RETURN
      END
C
C
C***** DTPRMY *****
      SUBROUTINE DTPRMY ( A, N, INDX, IER )
C
C----------------------------------------------------------------------
C
C ... PURPOSE  DTPRMY APPLIES THE INVERSE OF THE PERMUTATION IN INDX
C              TO THE DOUBLE PRECISION ARRAY A.  A(I) IS MOVED TO
C              A(INDX(I)).
C
C ... INPUT    N        NUMBER OF RECORDS IN A TO BE INTERCHANGED.
C
C              INDX      PERMUTATION ARRAY OF LENGTH N.
C
C ... INPUT/   A        DOUBLE PRECISION TO BE REORDERED.
C      OUTPUT
C
C ... OUTPUT   IER      SUCCESS/ERROR FLAG.
C
C                       IER = 0     NORMAL RETURN
C                           = -1    N .LE. 0
C                           = -2    INDX(1) .LE. 0
C
C----------------------------------------------------------------------
C
C ... GLOBAL VARIABLES
C
         DOUBLE PRECISION  A(1)
C
         INTEGER           N, INDX(1), IER
C
C----------------------------------------------------------------------
C
C ... LOCAL VARIABLES
C
C
         CHARACTER*8       NAME
C
         DOUBLE PRECISION  ATEMP1, ATEMP2
C
         INTEGER           I, NEXT, NOW
C
         DATA              NAME / 'DTPRMY' /
C
C----------------------------------------------------------------------
C
C ... CHECK INPUT
C
         IER = 0
            IF ( INDX(1) .LE. 0 )                     IER = -2
            IF ( N .LE. 0 )                          IER = -1
            IF ( IER .EQ. 0 ) GO TO 100
C
C ...... INPUT ERROR DETECTED.
C
               CALL DTERR ( 1, NAME, IER, 0 )
               GO TO 900
C
C ... APPLY THE INVERSE OF THE PERMUTATION
C
  100    IF ( N .EQ. 1 ) GO TO 900
C
C ...... SEARCH FOR THE FIRST ENTRY NOT PERMUTED WHICH IS INDICATED BY
C        A NONNEGATIVE VALUE IN INDX.
C
            DO 300 I = 1, N
               IF ( INDX(I) .LE. 0 ) GO TO 300
C
C ...... INITIALIZE TO FOLLOW THE CURRENT CHAIN OF PERMUTATIONS
C
               NEXT        = INDX(I)
               ATEMP1      = A(I)
C
C ...... FOLLOW THE CHAIN - PERMUTE AS YOU GO UNTIL THE CHAIN ENDS
C
  200          IF ( INDX(NEXT) .LE. 0 ) GO TO 300
C
                  ATEMP2      =  ATEMP1
                  ATEMP1      =  A(NEXT)
                  A(NEXT)     =  ATEMP2
C
                  NOW         =  NEXT
                  NEXT        =  INDX(NOW)
                  INDX(NOW)    = -NEXT
                  GO TO 200
C
C ...... END OF THE CHAIN AND END OF SEARCH LOOP FOR THE NEXT CHAIN
C
  300       CONTINUE
C
C----------------------------------------------------------------------
C
C ... PERMUTATION NOW FINISHED.  RESTORE THE INDX ARRAY.
C
            DO 400 I = 1, N
               INDX(I) = - INDX(I)
  400       CONTINUE
C
C----------------------------------------------------------------------
C
C ... END OF DTPRMY
C
  900    CONTINUE
         RETURN
      END
C
C
C***** DTSORT *****
      SUBROUTINE DTSORT ( A, N, IER )
C
C----------------------------------------------------------------------
C
C  PURPOSE  DTSORT SORTS A DOUBLE PRECISION REAL ARRAY OF N ELEMENTS
C           INTO THE STANDARD COLLATING SEQUENCE.
C
C  METHOD   DTSORT USES A QUICKSORT ALGORITHM IN THE STYLE OF THE
C           CACM PAPER BY BOB SEDGEWICK, OCTOBER 1978
C
C  INPUT    N        NUMBER OF ELEMENTS TO BE SORTED.
C
C  INPUT/   A        REAL ARRAY THAT IS TO BE SORTED.
C    OUTPUT
C
C  OUTPUT   IER      SUCCESS/ERROR FLAG.
C
C              IER =  0    NORMAL RETURN
C              IER = -1    N .LE. 0
C              IER = -5    INTERNAL FAILURE
C----------------------------------------------------------------------
C
C ... GLOBAL VARIABLES
C
      INTEGER           N, IER
C
      DOUBLE PRECISION  A(1)
C
C----------------------------------------------------------------------
C
C ... LOCAL VARIABLES
C
      INTEGER           I, IP1, J, JM1, LEFT, LLEN, RIGHT, RLEN,
     1                  STACK(50), STKLEN, TINY, TOP
C
      DOUBLE PRECISION  ATEMP(2)
C
      CHARACTER*8       NAME
C
      DATA              STKLEN / 50 /,  TINY / 9 /,
     A                  NAME / 'DTSORT'/
C
C----------------------------------------------------------------------
C
C ... CHECK INPUT
C
      IER = 0
C
      IF ( N .GT. 0 ) GO TO 100
C
         IER = -1
         CALL DTERR ( 1, NAME, IER, 0 )
         GO TO 9000
C
C----------------------------------------------------------------------
C
C ... PROGRAM IS A DIRECT TRANSLATION INTO FORTRAN OF SEDGEWICK'S
C     PROGRAM 2, WHICH IS NON-RECURSIVE, IGNORES FILES OF LENGTH
C     LESS THAN 'TINY' DURING PARTITIONING, AND USES MEDIAN OF THREE
C     PARTITIONING.
C
  100 IF (N .EQ. 1)  GO TO 9000
C
      TOP = 1
      LEFT = 1
      RIGHT = N
      IF ( N .LE. TINY ) GO TO 2000
C
C     ===========================================================
C     QUICKSORT -- PARTITION THE FILE UNTIL NO SUBFILE REMAINS OF
C     LENGTH GREATER THAN 'TINY'
C     ===========================================================
C
C     ... WHILE NOT DONE DO ...
C
C
C         ... FIND MEDIAN OF LEFT, RIGHT AND MIDDLE ELEMENTS OF CURRENT
C             SUBFILE, WHICH IS  A(LEFT), ..., A(RIGHT)
C             (CORRECTION TO CACM ARTICLE INCORPORATED)
C
  300     I        = ( LEFT + RIGHT ) / 2
          ATEMP(1) = A (I)
          A (I)    = A (LEFT)
          A (LEFT) = ATEMP(1)
C
          IF (  A(LEFT+1) .LE. A(LEFT)  ) GO TO 400
              ATEMP(1)       = A (LEFT+1)
              A (LEFT+1)     = A (LEFT)
              A (LEFT)       = ATEMP(1)
C
  400     IF (  A(LEFT) .LE. A(RIGHT)  )  GO TO 600
              ATEMP(1)      = A (LEFT)
              A (LEFT)      = A (RIGHT)
              A (RIGHT)     = ATEMP(1)
C
              IF (  A (LEFT+1) .LE. A (LEFT)  )  GO TO 600
                  ATEMP(1)       = A (LEFT+1)
                  A (LEFT+1)     = A (LEFT)
                  A (LEFT)       = ATEMP(1)
C
  600     ATEMP(2)   = A (LEFT)
C
C         ... ATEMP(2) IS NOW THE MEDIAN VALUE OF THE THREE A-S.  NOW
C             MOVE FROM THE LEFT AND RIGHT ENDS SIMULTANEOUSLY,
C             EXCHANGING A-S UNTIL ALL A-S LESS THAN  ATEMP(2) ARE
C             PACKED TO THE LEFT, ALL A-S LARGER THAN  ATEMP(2) ARE
C             PACKED TO THE RIGHT.
C
          I = LEFT+1
          J = RIGHT
C
C         LOOP
C             REPEAT I = I+1 UNTIL A(I) >= ATEMP(2)  ;
C             REPEAT J = J-1 UNTIL A(J) <= ATEMP(2)  ;
C         EXIT IF J < I;
C             << EXCHANGE AS I AND J >>
C         END
C
  700     CONTINUE
  800         I  = I + 1
              IF ( A(I) .LT. ATEMP(2)   )  GO TO 800
C
  900         J = J - 1
              IF ( A(J) .GT. ATEMP(2)   )  GO TO 900
C
          IF (J .LT. I)  GO TO 1000
              ATEMP(1)  = A (I)
              A (I)     = A (J)
              A (J)     = ATEMP(1)
          GO TO 700
C
 1000     ATEMP(1)     = A (LEFT)
          A (LEFT)     = A (J)
          A (J)        = ATEMP(1)
C
C
C         ... WE HAVE NOW PARTITIONED THE FILE INTO TWO SUBFILES,
C             ONE IS (LEFT ... J-1)  AND THE OTHER IS (I...RIGHT).
C             PROCESS THE SMALLER NEXT.  STACK THE LARGER ONE.
C
          LLEN = J-LEFT
          RLEN = RIGHT - I + 1
          IF ( RLEN .GT. LLEN ) GO TO 1100
C
C         ... LEFT SUBFILE ( LEFT ... J-1 ) IS THE LARGE SUBFILE.
C             TEST IF IT IS SMALL ENOUGH NOT TO PROCESS FURTHER.
C
              IF ( LLEN .LE. TINY ) GO TO 1200
C
C             ... IF RIGHT SUBFILE IS SMALL THEN PROCESS THE LEFT.
C                 ELSE STACK THE LEFT AND PROCESS THE RIGHT.
C
                  IF ( RLEN .GT. TINY ) GO TO 1050
                       RIGHT = J - 1
                       GO TO 300
C
 1050                  IF ( TOP .GE. STKLEN ) GO TO 8000
                           STACK ( TOP    ) = LEFT
                           STACK ( TOP+1  ) = J - 1
                           TOP              = TOP + 2
                           LEFT             = I
                           GO TO 300
C
C        ... RIGHT SUBFILE ( I ... RIGHT ) IS THE LARGE SUBFILE.
C            TEST IF IT IS SMALL ENOUGH NOT TO PROCESS FURTHER.
C
 1100        IF ( RLEN .LE. TINY ) GO TO 1200
C
C            ... IF LEFT SUBFILE IS SMALL THEN PROCESS THE RIGHT.
C                ELSE STACK THE RIGHT AND PROCESS THE LEFT.
C
                 IF ( LLEN .GT. TINY ) GO TO 1150
                     LEFT = I
                     GO TO 300
C
 1150                IF ( TOP .GE. STKLEN ) GO TO 8000
                         STACK ( TOP     ) = I
                         STACK ( TOP + 1 ) = RIGHT
                         TOP               = TOP + 2
                         RIGHT             = J - 1
                         GO TO 300
C
C        ... BOTH LEFT AND RIGHT SUBFILE ARE SMALL.  POP THE STACK
C            TO GET THE NEXT SUBFILE TO PROCESS.
C
 1200        IF ( TOP .EQ. 1 ) GO TO 2000
                 TOP    = TOP - 2
                 LEFT   = STACK ( TOP )
                 RIGHT  = STACK ( TOP + 1 )
                 GO TO 300
C
C     ------------------------------------------------------------
C     INSERTION SORT THE ENTIRE FILE, WHICH CONSISTS OF A LIST
C     OF 'TINY' SUBFILES, LOCALLY OUT OF ORDER, GLOBALLY IN ORDER.
C     ------------------------------------------------------------
C
C     ... INSERTION SORT ... FOR I := N-1 STEP -1 TO 1 DO ...
C
 2000 I   = N - 1
      IP1 = N
C
 2100     IF (  A (I) .LE. A (IP1)  )  GO TO 2400
C
C             ... OUT OF ORDER ... MOVE UP TO CORRECT PLACE
C
              ATEMP(1)  = A (I)
              J         = IP1
              JM1       = I
C
C             ... REPEAT ... UNTIL 'CORRECT PLACE FOR ATEMP(1)   FOUND'
C
 2200             A (JM1)     = A (J)
                  JM1         = J
                  J           = J + 1
                  IF ( J .GT. N )  GO TO 2300
                  IF (  A (J) .LT. ATEMP(1)    )  GO TO 2200
C
 2300         A (JM1)     = ATEMP(1)
C
 2400     IP1 = I
          I   = I - 1
          IF ( I .GT. 0 )  GO TO 2100
C
      GO TO 9000
C
C---------------------------------------------------------------------
C
 8000 IER = -5
      CALL DTERR ( 3, NAME, IER, 0 )
C
C----------------------------------------------------------------------
C
C ... END OF DTSORT
C
 9000    CONTINUE
         RETURN
      END
C
C
C***** DTSRTN *****
      SUBROUTINE DTSRTN ( A, N, IOPT, ISTBLE, INDX, IER )
C
C----------------------------------------------------------------------
C
C  PURPOSE  DTSRTN SORTS A DOUBLE PRECISION ARRAY OF N ELEMENTS INTO
C           THE ASCENDING ORDER AND GENERATES THE PERMUTATION
C           ARRAY INDX.
C
C  METHOD   DTSRTN USES A QUICKSORT ALGORITHM IN THE STYLE OF THE
C           CACM PAPER BY BOB SEDGEWICK, OCTOBER 1978
C
C  INPUT    N        NUMBER OF ELEMENTS TO BE SORTED.
C
C           IOPT     SORTING OPTION.  IF IOPT = 0, THE ARRAY A
C                    IS RETURNED IN SORTED ORDER.  OTHERWISE A
C                    IS RETURNED IN ORIGINAL ORDER.
C
C           ISTBLE   STABILIZATION OPTION.  IF ISTBLE = 0, THE
C                    SORT IS STABILIZED.  OTHERWISE IT IS NOT.
C
C  INPUT/   A        ARRAY THAT IS TO BE SORTED.  IF IOPT = 0
C    OUTPUT          AND IER = 0 A IS RETURNED IN SORTED ORDER
C                    IF IER = -5, A IS PARTIALLY SORTED.  OTHERWISE
C                    A IS LEFT UNCHANGED.
C
C  OUTPUT   INDX      PERMUTATION ARRAY.
C
C           IER      SUCCESS/ERROR FLAG.
C
C              IER =  0    NORMAL RETURN
C              IER = -1    N .LE. 0
C              IER = -5    INTERNAL FAILURE
C
C
C ... GLOBAL VARIABLES
C
      INTEGER           N, INDX(1), IER, IOPT, ISTBLE
C
      DOUBLE PRECISION  A(1)
C
C----------------------------------------------------------------------
C
C ... LOCAL VARIABLES
C
      INTEGER     I, IP1, J, JM1, K, L, LEFT, LLEN, RIGHT,
     A            RLEN, STACK(50),TOP,STKLEN,TINY
C
      CHARACTER*8   NAME
      DOUBLE PRECISION   ATEMP(2)
C
C
      DATA        STKLEN, TINY /50,9 /,
     A            NAME / 'DTSRTN' /
C
C----------------------------------------------------------------------
C
C ... CHECK INPUT
C
      IER = 0
C
      IF ( N .LE. 0 ) IER = - 1
C
      IF ( IER .EQ. 0 ) GO TO 100
C
         INDX(1) = 0
         CALL DTERR ( 1, NAME, IER, 0 )
         GO TO 9000
C
C----------------------------------------------------------------------
C
C ... PROGRAM IS A DIRECT TRANSLATION INTO FORTRAN OF SEDGEWICK'S
C     PROGRAM 2, WHICH IS NON-RECURSIVE, IGNORES FILES OF LENGTH
C     LESS THAN 'TINY' DURING PARTITIONING, AND USES MEDIAN OF THREE
C     PARTITIONING.
C
C ... SET INDX TO THE INDENTITY PERMUTATION
C
  100 DO 200 I = 1, N
         INDX(I) = I
  200 CONTINUE
C
      IF (N .EQ. 1)  GO TO 9000
C
      TOP = 1
      LEFT = 1
      RIGHT = N
      IF ( N .LE. TINY) GO TO 2000
C
C     ===========================================================
C     QUICKSORT -- PARTITION THE FILE UNTIL NO SUBFILE REMAINS OF
C     LENGTH GREATER THAN 'TINY'
C     ===========================================================
C
C     ... WHILE NOT DONE DO ...
C
C
C         ... FIND MEDIAN OF LEFT, RIGHT AND MIDDLE ELEMENTS OF CURRENT
C             SUBFILE, WHICH IS  A(LEFT), ..., A(RIGHT)
C             (CORRECTION TO CACM ARTICLE INCORPORATED)
C
  300     I            = ( LEFT + RIGHT ) / 2
C
          ATEMP ( 1 )  = A(I)
          A(I)         = A ( LEFT )
          A ( LEFT )   = ATEMP ( 1 )
C
          K            = INDX (I)
          INDX (I)      = INDX (LEFT)
          INDX (LEFT)   = K
C
          IF ( A(LEFT+1) .LE. A(LEFT) ) GO TO 400
C
              ATEMP ( 1 )  = A ( LEFT+1 )
              A ( LEFT+1 ) = A ( LEFT )
              A ( LEFT )   = ATEMP ( 1 )
C
              K            = INDX (LEFT+1)
              INDX (LEFT+1) = INDX (LEFT)
              INDX (LEFT)  = K
C
  400     IF ( A(LEFT) .LE. A(RIGHT) ) GO TO 600
C
              ATEMP ( 1 ) = A ( LEFT )
              A ( LEFT )  = A ( RIGHT )
              A ( RIGHT ) = ATEMP ( 1 )
C
              K           = INDX (LEFT)
              INDX (LEFT)  = INDX (RIGHT)
              INDX (RIGHT) = K
C
              IF ( A(LEFT+1) .LE. A(LEFT) ) GO TO 600
C
                  ATEMP ( 1 )   = A ( LEFT+1 )
                  A ( LEFT+1 )  = A ( LEFT )
                  A ( LEFT )    = ATEMP ( 1 )
C
                  K            = INDX (LEFT+1)
                  INDX (LEFT+1) = INDX (LEFT)
                  INDX (LEFT)   = K
C
  600     ATEMP(2) = A (LEFT)
C
C       ... ATEMP(2) IS NOW THE MEDIAN VALUE OF THE THREE A-S.  NOW MOVE
C           FROM THE LEFT AND RIGHT ENDS SIMULTANEOUSLY, EXCHANGING
C           INDX UNTIL ALL A-S LESS THAN  ATEMP(2)  ARE SYMBOLLICALLY
C           PACKED TO THE LEFT, ALL A-S LARGER THAN  ATEMP(2)  ARE
C           PACKED TO THE RIGHT.
C
          I = LEFT+1
          J = RIGHT
C
C         LOOP
C             REPEAT I = I+1 UNTIL A(I) >= ATEMP;
C             REPEAT J = J-1 UNTIL A(J) <= ATEMP;
C         EXIT IF J < I;
C             << EXCHANGE AS I AND J >>
C         END
C
  700     CONTINUE
  800         I  = I + 1
              IF ( A(I) .LT. ATEMP(2) )  GO TO 800
C
  900         J  = J - 1
              IF ( A(J) .GT. ATEMP(2) )  GO TO 900
C
          IF (J .LT. I)  GO TO 1000
C
              ATEMP(1)  = A(I)
              A(I)      = A(J)
              A(J)      = ATEMP(1)
C
              K         = INDX (I)
              INDX (I)   = INDX (J)
              INDX (J)   = K
          GO TO 700
C
 1000     ATEMP(1)   = A (LEFT)
          A(LEFT)    = A (J)
          A(J)       = ATEMP(1)
C
          K          = INDX (LEFT)
          INDX (LEFT) = INDX (J)
          INDX (J)    = K
C
C
C         ... WE HAVE NOW PARTITIONED THE FILE INTO TWO SUBFILES,
C             ONE IS (LEFT ... J-1)  AND THE OTHER IS (I...RIGHT).
C             PROCESS THE SMALLER NEXT.  STACK THE LARGER ONE.
C
          LLEN = J-LEFT
          RLEN = RIGHT - I + 1
          IF ( RLEN .GT. LLEN ) GO TO 1100
C
C         ... LEFT SUBFILE ( LEFT ... J-1 ) IS THE LARGE SUBFILE.
C             TEST IF IT IS SMALL ENOUGH NOT TO PROCESS FURTHER.
C
              IF ( LLEN .LE. TINY ) GO TO 1200
C
C             ... IF RIGHT SUBFILE IS SMALL THEN PROCESS THE LEFT.
C                 ELSE STACK THE LEFT AND PROCESS THE RIGHT.
C
                  IF ( RLEN .GT. TINY ) GO TO 1050
                       RIGHT = J - 1
                       GO TO 300
C
 1050                  IF ( TOP .GE. STKLEN ) GO TO 8000
                           STACK ( TOP    ) = LEFT
                           STACK ( TOP+1  ) = J - 1
                           TOP              = TOP + 2
                           LEFT             = I
                           GO TO 300
C
C        ... RIGHT SUBFILE ( I ... RIGHT ) IS THE LARGE SUBFILE.
C            TEST IF IT IS SMALL ENOUGH NOT TO PROCESS FURTHER.
C
 1100        IF ( RLEN .LE. TINY ) GO TO 1200
C
C            ... IF LEFT SUBFILE IS SMALL THEN PROCESS THE RIGHT.
C                ELSE STACK THE RIGHT AND PROCESS THE LEFT.
C
                 IF ( LLEN .GT. TINY ) GO TO 1150
                     LEFT = I
                     GO TO 300
C
 1150                IF ( TOP .GE. STKLEN ) GO TO 8000
                         STACK ( TOP     ) = I
                         STACK ( TOP + 1 ) = RIGHT
                         TOP               = TOP + 2
                         RIGHT             = J - 1
                         GO TO 300
C
C        ... BOTH LEFT AND RIGHT SUBFILE ARE SMALL.  POP THE STACK
C            TO GET THE NEXT SUBFILE TO PROCESS.
C
 1200        IF ( TOP .EQ. 1 ) GO TO 2000
                 TOP    = TOP - 2
                 LEFT   = STACK ( TOP )
                 RIGHT  = STACK ( TOP + 1 )
                 GO TO 300
C
C     ------------------------------------------------------------
C     INSERTION SORT THE ENTIRE FILE, WHICH CONSISTS OF A LIST
C     OF 'TINY' SUBFILES, LOCALLY OUT OF ORDER, GLOBALLY IN ORDER.
C     ------------------------------------------------------------
C
C     ... INSERTION SORT ... FOR I := N-1 STEP -1 TO 1 DO ...
C
 2000 I   = N - 1
      IP1 = N
C
 2100     IF (  A (I) .LE. A (IP1) )  GO TO 2400
C
C             ... OUT OF ORDER ... MOVE UP TO CORRECT PLACE
C
              ATEMP(1)  = A   (I)
              K         = INDX (I)
              J         = IP1
              JM1       = I
C
C             ... REPEAT ... UNTIL 'CORRECT PLACE FOR ATEMP FOUND'
C
 2200             A   (JM1) = A   (J)
                  INDX (JM1) = INDX (J)
                  JM1       = J
                  J         = J + 1
                  IF ( J .GT. N )  GO TO 2300
                  IF (  A (J) .LT. ATEMP(1) )  GO TO 2200
C
 2300         A   (JM1) = ATEMP(1)
              INDX (JM1) = K
C
 2400     IP1 = I
          I   = I - 1
          IF ( I .GT. 0 )  GO TO 2100
C
C     -----------------------------------------------------------
C     FILE IS NOW SORTED AND INDX IS THE POINTER.
C     STABILIZE THE SORT IF REQUESTED.
C     -----------------------------------------------------------
C
      IF ( ISTBLE .NE. 0 ) GO TO 4000
C
C ... SEARCH FOR A SEQUENCE OF IDENTICAL ENTRIES
C
      I = 1
C
 3100 IP1 = I + 1
      IF ( A(I) .EQ. A(IP1) ) GO TO 3200
         I = IP1
         IF ( I .GE. N ) GO TO 4000
         GO TO 3100
C
C ... 2 ENTRIES ARE IDENTICAL.  FIND END OF SEQUENCE
C
 3200 J = IP1 + 1
C
 3300 IF ( J .GT. N ) GO TO 3400
         IF ( A(I) .NE. A(J) ) GO TO 3400
            J = J + 1
            GO TO 3300
C
C ... SORT INDX FOR IDENTICAL ENTRIES
C
 3400 L = J - I
      CALL DJSORT ( INDX(I), L, IER )
      IF ( IER .NE. 0 ) GO TO 8100
      I = J
      IF ( I .GE. N ) GO TO 4000
      GO TO 3100
C
C--------------------------------------------------------------------
C
C     ----------------------------------------
C     IF REQUESTED, RETURN A IN ORGINAL ORDER.
C     ----------------------------------------
C
 4000 IF ( IOPT .EQ. 0 ) GO TO 9000
C
      CALL DTPRMY ( A, N, INDX, IER )
      IF( IER .NE. 0 ) GO TO 8200
C
      GO TO 9000
C
C---------------------------------------------------------------------
C
C ... ERROR TRAPS
C
C ... INTERNAL STACK OVERFLOW
C
 8000 IER = -5
      INDX(1) = 0
      CALL DTERR ( 3, NAME, IER, 0 )
      GO TO 9000
C
C ... UNEXPECTED ERROR RETURN FROM DJSORT
C
 8100 IER = -6
      GO TO 8300
C
C ... UNEXPECTED ERROR RETURN FROM DTPRMX
C
 8200 IER = -7
 8300 INDX(1) = 0
      CALL DTERR ( 5, NAME, IER, 0 )
C
C----------------------------------------------------------------------
C
C ... END OF DTSRTN
C
 9000    CONTINUE
         RETURN
      END
C
C
C***** DJSORT *****
      SUBROUTINE DJSORT ( A, N, IER )
C
C----------------------------------------------------------------------
C
C  PURPOSE  DJSORT SORTS AN INTEGER ARRAY OF N ELEMENTS INTO THE
C           STANDARD COLLATING SEQUENCE.
C
C  METHOD   DJSORT USES A QUICKSORT ALGORITHM IN THE STYLE OF THE
C           CACM PAPER BY BOB SEDGEWICK, OCTOBER 1978
C
C  INPUT    N        NUMBER OF ELEMENTS TO BE SORTED.
C
C  INPUT/   A        ARRAY THAT IS TO BE SORTED.
C    OUTPUT
C
C  OUTPUT   IER      SUCCESS/ERROR FLAG.
C
C              IER =  0    NORMAL RETURN
C              IER = -1    N .LE. 0
C              IER = -5    INTERNAL FAILURE
C
C
C ... GLOBAL VARIABLES
C
      INTEGER        N, IER
C
      INTEGER        A(1)
C
C----------------------------------------------------------------------
C
C ... LOCAL VARIABLES
C
      INTEGER     I, IP1, J, JM1, LEFT, LLEN, RIGHT, RLEN,
     1            STACK(50), STKLEN, TINY, TOP
C
C
      INTEGER     ATEMP(2)
      CHARACTER*8  NAME
C
      DATA        STKLEN / 50 /,  TINY / 9 /,
     A            NAME / 'DJSORT'/
C
C----------------------------------------------------------------------
C
C ... CHECK INPUT
C
      IER = 0
C
      IF ( N .GT. 0 ) GO TO 100
C
         IER = -1
         CALL DTERR ( 1, NAME, IER, 0 )
         GO TO 9000
C
C----------------------------------------------------------------------
C
C ... PROGRAM IS A DIRECT TRANSLATION INTO FORTRAN OF SEDGEWICK'S
C     PROGRAM 2, WHICH IS NON-RECURSIVE, IGNORES FILES OF LENGTH
C     LESS THAN 'TINY' DURING PARTITIONING, AND USES MEDIAN OF THREE
C     PARTITIONING.
C
  100 IF (N .EQ. 1)  GO TO 9000
C
      TOP = 1
      LEFT = 1
      RIGHT = N
      IF ( N .LE. TINY ) GO TO 2000
C
C     ===========================================================
C     QUICKSORT -- PARTITION THE FILE UNTIL NO SUBFILE REMAINS OF
C     LENGTH GREATER THAN 'TINY'
C     ===========================================================
C
C     ... WHILE NOT DONE DO ...
C
C
C         ... FIND MEDIAN OF LEFT, RIGHT AND MIDDLE ELEMENTS OF CURRENT
C             SUBFILE, WHICH IS  A(LEFT), ..., A(RIGHT)
C             (CORRECTION TO CACM ARTICLE INCORPORATED)
C
  300     I        = ( LEFT + RIGHT ) / 2
          ATEMP(1) = A (I)
          A (I)    = A (LEFT)
          A (LEFT) = ATEMP(1)
C
          IF (  A(LEFT+1) .LE. A(LEFT)  ) GO TO 400
              ATEMP(1)       = A (LEFT+1)
              A (LEFT+1)     = A (LEFT)
              A (LEFT)       = ATEMP(1)
C
  400     IF (  A(LEFT) .LE. A(RIGHT)  )  GO TO 600
              ATEMP(1)      = A (LEFT)
              A (LEFT)      = A (RIGHT)
              A (RIGHT)     = ATEMP(1)
C
              IF (  A (LEFT+1) .LE. A (LEFT)  )  GO TO 600
                  ATEMP(1)       = A (LEFT+1)
                  A (LEFT+1)     = A (LEFT)
                  A (LEFT)       = ATEMP(1)
C
  600     ATEMP(2)   = A (LEFT)
C
C         ... ATEMP(2) IS NOW THE MEDIAN VALUE OF THE THREE A-S.  NOW
C             MOVE FROM THE LEFT AND RIGHT ENDS SIMULTANEOUSLY,
C             EXCHANGING A-S UNTIL ALL A-S LESS THAN  ATEMP(2) ARE
C             PACKED TO THE LEFT, ALL A-S LARGER THAN  ATEMP(2) ARE
C             PACKED TO THE RIGHT.
C
          I = LEFT+1
          J = RIGHT
C
C         LOOP
C             REPEAT I = I+1 UNTIL A(I) >= ATEMP(2)  ;
C             REPEAT J = J-1 UNTIL A(J) <= ATEMP(2)  ;
C         EXIT IF J < I;
C             << EXCHANGE AS I AND J >>
C         END
C
  700     CONTINUE
  800         I  = I + 1
              IF ( A(I) .LT. ATEMP(2)   )  GO TO 800
C
  900         J = J - 1
              IF ( A(J) .GT. ATEMP(2)   )  GO TO 900
C
          IF (J .LT. I)  GO TO 1000
              ATEMP(1)  = A (I)
              A (I)     = A (J)
              A (J)     = ATEMP(1)
          GO TO 700
C
 1000     ATEMP(1)     = A (LEFT)
          A (LEFT)     = A (J)
          A (J)        = ATEMP(1)
C
C
C         ... WE HAVE NOW PARTITIONED THE FILE INTO TWO SUBFILES,
C             ONE IS (LEFT ... J-1)  AND THE OTHER IS (I...RIGHT).
C             PROCESS THE SMALLER NEXT.  STACK THE LARGER ONE.
C
          LLEN = J-LEFT
          RLEN = RIGHT - I + 1
          IF ( RLEN .GT. LLEN ) GO TO 1100
C
C         ... LEFT SUBFILE ( LEFT ... J-1 ) IS THE LARGE SUBFILE.
C             TEST IF IT IS SMALL ENOUGH NOT TO PROCESS FURTHER.
C
              IF ( LLEN .LE. TINY ) GO TO 1200
C
C             ... IF RIGHT SUBFILE IS SMALL THEN PROCESS THE LEFT.
C                 ELSE STACK THE LEFT AND PROCESS THE RIGHT.
C
                  IF ( RLEN .GT. TINY ) GO TO 1050
                       RIGHT = J - 1
                       GO TO 300
C
 1050                  IF ( TOP .GE. STKLEN ) GO TO 8000
                           STACK ( TOP    ) = LEFT
                           STACK ( TOP+1  ) = J - 1
                           TOP              = TOP + 2
                           LEFT             = I
                           GO TO 300
C
C        ... RIGHT SUBFILE ( I ... RIGHT ) IS THE LARGE SUBFILE.
C            TEST IF IT IS SMALL ENOUGH NOT TO PROCESS FURTHER.
C
 1100        IF ( RLEN .LE. TINY ) GO TO 1200
C
C            ... IF LEFT SUBFILE IS SMALL THEN PROCESS THE RIGHT.
C                ELSE STACK THE RIGHT AND PROCESS THE LEFT.
C
                 IF ( LLEN .GT. TINY ) GO TO 1150
                     LEFT = I
                     GO TO 300
C
 1150                IF ( TOP .GE. STKLEN ) GO TO 8000
                         STACK ( TOP     ) = I
                         STACK ( TOP + 1 ) = RIGHT
                         TOP               = TOP + 2
                         RIGHT             = J - 1
                         GO TO 300
C
C        ... BOTH LEFT AND RIGHT SUBFILE ARE SMALL.  POP THE STACK
C            TO GET THE NEXT SUBFILE TO PROCESS.
C
 1200        IF ( TOP .EQ. 1 ) GO TO 2000
                 TOP    = TOP - 2
                 LEFT   = STACK ( TOP )
                 RIGHT  = STACK ( TOP + 1 )
                 GO TO 300
C
C     ------------------------------------------------------------
C     INSERTION SORT THE ENTIRE FILE, WHICH CONSISTS OF A LIST
C     OF 'TINY' SUBFILES, LOCALLY OUT OF ORDER, GLOBALLY IN ORDER.
C     ------------------------------------------------------------
C
C     ... INSERTION SORT ... FOR I := N-1 STEP -1 TO 1 DO ...
C
 2000 I   = N - 1
      IP1 = N
C
 2100     IF (  A (I) .LE. A (IP1)  )  GO TO 2400
C
C             ... OUT OF ORDER ... MOVE UP TO CORRECT PLACE
C
              ATEMP(1)  = A (I)
              J         = IP1
              JM1       = I
C
C             ... REPEAT ... UNTIL 'CORRECT PLACE FOR ATEMP(1)   FOUND'
C
 2200             A (JM1)     = A (J)
                  JM1         = J
                  J           = J + 1
                  IF ( J .GT. N )  GO TO 2300
                  IF (  A (J) .LT. ATEMP(1)    )  GO TO 2200
C
 2300         A (JM1)     = ATEMP(1)
C
 2400     IP1 = I
          I   = I - 1
          IF ( I .GT. 0 )  GO TO 2100
C
      GO TO 9000
C
C---------------------------------------------------------------------
C
 8000 IER = -5
      CALL DTERR ( 3, NAME, IER, 0 )
C
C----------------------------------------------------------------------
C
C ... END OF DJSORT
C
 9000    CONTINUE
         RETURN
      END
