C+----------------------------------------------------------------------
C
      REAL FUNCTION URAND(IY)
C
      INTEGER  IY
C
C      URAND is a uniform random number generator based  on  theory  and
C  suggestions  given  in  D.E. Knuth (1969),  Vol  2.   The integer  IY
C  should be initialized to an arbitrary integer prior to the first call
C  to URAND.  The calling program should  not  alter  the  value  of  IY
C  between  subsequent calls to URAND.  Values of urand will be returned
C  in the interval (0,1).
C
C     This version matches the published version exactly:  its output is
C  in single precision but the internal floating point arithmetic is  in
C  double precision. (Actually a SAVE statement has been introduced here
C  to be sure of distinguishing first and subsequent calls.)
C
C     The companion book states that the least significant binary digits
C  of the (internal) IY are not particularly random,  but  when they are
C  converted to floating point numbers the least significant digits  are
C  usually not important. The user is nonetheless advised to pretest the
C  results by trying them on some problem similar to the  real  applica-
C  tion, for which the answer is known.
C
C     Further notes:
C
C  > Binary integer arithmetic is assumed.
C
C  > Use INT (K * URAND (IY)) to compute a random integer between  0 and
C    K - 1.
C
C  > The code appears to guard against integer overflow, yet VAX/VMS has
C    been shown to fail unless it is compiled with /CHECK=NOOVERFLOW.
C
C-----------------------------------------------------------------------
C
      INTEGER  IA,IC,ITWO,M2,M,MIC
      REAL  HALFM
      REAL  S
      REAL  ATAN,SQRT
      DATA  M2/0/,ITWO/2/
      SAVE  M2,ITWO

      IF (M2 .NE. 0) GO TO 20
C
C  IF FIRST ENTRY, COMPUTE MACHINE INTEGER WORD LENGTH
C
      M = 1
   10 M2 = M
      M = ITWO*M2
      IF (M .GT. M2) GO TO 10
      HALFM = M2
C
C  COMPUTE MULTIPLIER AND INCREMENT FOR LINEAR CONGRUENTIAL METHOD
C
      IA = 8*INT(HALFM*ATAN(1.E+0)/8.E+0) + 5
      IC = 2*INT(HALFM*(0.5E+0-SQRT(3.E+0)/6.E+0)) + 1
      MIC = (M2 - IC) + M2
C
C  S IS THE SCALE FACTOR FOR CONVERTING TO FLOATING POINT
C
      S = 0.5/HALFM
C
C  COMPUTE NEXT RANDOM NUMBER
C
   20 IY = IY*IA
C
C  THE FOLLOWING STATEMENT IS FOR COMPUTERS WHICH DO NOT ALLOW
C  INTEGER OVERFLOW ON ADDITION
C
      IF (IY .GT. MIC) IY = (IY - M2) - M2
C
      IY = IY + IC
C
C  THE FOLLOWING STATEMENT IS FOR COMPUTERS WHERE THE
C  WORD LENGTH FOR ADDITION IS GREATER THAN FOR MULTIPLICATION
C
      IF (IY/2 .GT. M2) IY = (IY - M2) - M2
C
C  THE FOLLOWING STATEMENT IS FOR COMPUTERS WHERE INTEGER
C  OVERFLOW AFFECTS THE SIGN BIT
C
      IF (IY .LT. 0) IY = (IY + M2) + M2
      URAND = FLOAT(IY)*S
      RETURN
      END
