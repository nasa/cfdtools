      SUBROUTINE DTRPUN(N,X)
C
C     DTRPUN IS A PORTABLE UNIFORM RANDOM NUMBER GENERATOR.
C     (PORTABILITY ASSUMES AT LEAST 16 BIT INTEGER ARITHMETIC.)
C
C PURPOSE
C
C     DTRPUN COMPUTES UNIFORMLY DISTRIBUTED RANDOM VARIABLES X WITH
C
C                        ( 0, T.LE.0      )
C         PROB(X.LT.T) = ( T, 0.LT.T.LT.1 ) .
C                        ( 1, T.GE.1      )
C
C USAGE
C
C     BEFORE CALLING THIS ROUTINE (OR ANY OTHER ROUTINE USING THIS
C     ROUTINE), THE INITIALIZATION ROUTINE DTRPST MUST BE CALLED ONCE.
C
C         REAL X(N)
C         INIT = ...
C         CALL DTRPST(INIT)
C         CALL DTRPUN(N,X)
C
C INPUT
C
C     N   THE NUMBER OF UNIFORM VARIABLES DESIRED
C
C OUTPUT
C
C     X   THE ARRAY OF UNIFORM VARIABLES WHICH SATISFIES 2**(-21).LE.X
C         .LE.1-2**(-21)  (IF N.LE.0, NO OUTPUT IS RETURNED.)
C
C COMMON
C
C     THE INTERNAL COMMUNICATION COMMON BLOCK DTRPCM WILL BE MODIFIED.
C
C METHOD
C
C     A MIXED CONGRUENTIAL GENERATOR IS USED.  THE MULTIPLIER IS
C     5**15.  THE CONSTANT IS 7261067085.  THE MODULUS IS 2**35.
C     THE ARITHMETIC IS DONE WITH MULTIPLE PRECISION BASE 2**7
C     INTEGERS.
C
C PROGRAMMED BY
C
C     STUART L. ANDERSON - 28 MAR. 1979
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      REAL    X(N)
      INTEGER ACCUM(5),A(5),C(5)
      REAL TWOPI, XNORM
      INTEGER INORM, ISEED
      COMMON/DTRPCM/TWOPI,XNORM,INORM,ISEED(5)
      SAVE /DTRPCM/
C
      INTEGER NLOCAL, I, J, MULT, ICARRY, K, ITEMP
      REAL XTEMP
C                                                 BASE 128 CONSTANTS
      DATA A(1),A(2),A(3),A(4),A(5)/113,87,117,19,13/
      DATA C(1),C(2),C(3),C(4),C(5)/ 27, 6, 44,46,77/
C
      NLOCAL = N
      IF(NLOCAL.LE.0)GO TO 200
      INDEX = 1
C
C****     FORM THE NEW SEED WITH 35 BIT ARITHMETIC
C
   10     DO 20 I=1,5
              ACCUM(I) = C(I)
   20     CONTINUE
          J = 5
   30         MULT   = ISEED(J)
              ICARRY = 0
              I      = 5
              K      = J
   40             ITEMP    = MULT*A(I)+ACCUM(K)+ICARRY
                  ICARRY   = ITEMP/128
                  ACCUM(K) = ITEMP-128*ICARRY
                  I        = I-1
                  K        = K-1
                  IF(K.GT.0)GO TO 40
              J = J-1
              IF(J.GT.0)GO TO 30
          DO 50 I=1,5
              ISEED(I) = ACCUM(I)
   50     CONTINUE
C
C****     FLOAT THE RESULT
C
          XTEMP = 0.
          I     = 3
   60         XTEMP = (XTEMP+FLOAT(ACCUM(I)))/128.
              I     = I-1
              IF(I.GT.0)GO TO 60
C
C****     IF XTEMP.NE.0, OUTPUT IT
C
          IF(XTEMP.EQ.0.)GO TO 70
          X(INDEX) = XTEMP
          INDEX    = INDEX+1
   70     IF(INDEX.LE.NLOCAL)GO TO 10
C
  200 RETURN
      END
