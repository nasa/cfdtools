      FUNCTION D0TBU1 (XB,X,NX,NORD,NEXP,IS)
C ..
C     THIS FUNCTION IS A REWRITE OF THE OLD ETSRCD
C     THE DIFFERENCE (EXCEPT FOR NAME CHANGES) IS THAT THIS
C     RETURNS A NEGATIVE INDEX IF EXTRAPOLATION IS CALLED FOR.
C     THIS WAY THE INTERPOLATOR CAN ISSUE THE APPROPRIATE
C     WARNING.
C
C     INPUT
C
C        XB        - VALUE TO LOOK-UP
C        X         - ARRAY OF VALUE TO SEARCH.
C        NX        - NUMBER OF POINTS IN X.
C        NORD      - ORDER OF POLYNOMIAL THAT IS TO BE USED
C        NEXP      - =0 FOR EDGE VALUE, < 0 FOR EXTRAPOLATION
C
C     OUTPUT
C
C        D0TBU1    - INDICATOR 0 IF EXACT VALUE WAS FOUND,
C                              1 IF INTERPOLATION IS NEEDED.
C        IS        - THE TABLE INDEX FOR THE SMALLEST SUBSCRIPT
C                    TO BE USED IN THE POLYNOMIAL.
C
C     MODIFIED    BY TOM TOSCH - DECEMBER 1986.
C ...............................................................
      DOUBLE PRECISION XB     , X(*)        , DIF        ,
     +                 ASCEND
      INTEGER        D0TBU1   , NORD        , NX         ,
     +               IS       , IHALF       ,
     +               N        , IDLT        ,
     +               IMULT    , NEXP
C
      IHALF (N) = (N + 1) / 2
C
      N         = NX
      D0TBU1    = 1
      IMULT     = 1
      IS        = 1
C ...........................................N .EQ. 1 , DONE
      IF ( N .LE. 1 ) RETURN
C ...........................................EXTRAPOLATION ?
      IF ( X(1) .LE. XB  .AND.  XB .LE. X(N) ) GO TO 10
      IF ( X(N) .LE. XB  .AND.  XB .LE. X(1) ) GO TO 10
      IMULT     = -1
      IF ( NEXP .NE. 0 )                       GO TO 10
C ...........................................EXTRAP. NOT LEGAL
      IS        = -1
      IF ( X(N) .GT. X(1)  .AND. XB .GT. X(N) ) IS = -N
      IF ( X(N) .LT. X(1)  .AND. XB .LT. X(N) ) IS = -N
      RETURN
C ...........................................ASCEND. OR DESCEND.
10    ASCEND       = 1.0D0
      IF ( X(1) .GT. X(2) ) ASCEND = -ASCEND
C
C ...........................................BINARY SEARCH
C
      IS        = IHALF (N)
      IDLT      = IS
C
C ..........MAIN LOOP
C
20    IDLT      = IHALF (IDLT)
      DIF       = ASCEND*(X(IS) - XB)
      IF ( DIF .EQ. 0.D0  ) RETURN
      IF ( DIF .GT. 0.D0  ) THEN
C ......................................IS IS TOO LARGE
C                                       (UNLESS IS = 1)
         IF ( IS .LE. 1 ) GO TO 80
C ......................................IDLT TOO LARGE
C                                       (BECAUSE N NOT POWER OF 2)
         IF ( IS .LE. IDLT ) IDLT = IHALF (IDLT)
         IS     = IS - IDLT
         IS     = MAX0 (IS,1)
         GO TO 20
      END IF
C .....................................IS OK OR TOO SMALL
C                                      (UNLESS IS = N)
      IF ( IS .LT. N ) THEN
C .....................................XB NOT OUTSIDE RIGHT END
         DIF    = ASCEND*(X(IS+1) - XB)
         IF ( DIF .EQ. 0.D0 ) THEN
            IS  = IS + 1
            RETURN
         ELSE IF ( DIF .LT. 0.D0 ) THEN
C .....................................IS IS TOO SMALL
40          IS  = IS + IDLT
            IF ( IS .LE. N ) GO TO 20
            IS  = IS - IDLT
            IDLT= IHALF (IDLT)
            GO TO 40
         END IF
      END IF
C
C .....................................XB BRACKETED BY X(IS)
C                                      X(IS+1)
80    IF ( NORD .LE. 0 ) THEN
C .....................................NORD .LE. 0, GET NEAREST PT.
         IF ( IS .LT. N  .AND.  ABS(XB-X(IS  )) .GT.
     1                          ABS(XB-X(IS+1)) )     IS = IS + 1
      ELSE
C ....................................FIND NORD POINTS CENTERED (IF
C                                     POSSIBLE) AROUND XB
         IS     = MIN0( MAX0(1 , IS-(NORD-1)/2) , N-NORD)
         D0TBU1 = 0
      END IF
      IS        = IS * IMULT
      RETURN
      END
