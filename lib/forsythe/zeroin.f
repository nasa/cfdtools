C+----------------------------------------------------------------------
C
      REAL FUNCTION ZEROIN ( AX, BX, F, TOL, IER )
      INTEGER       IER
      REAL          AX,BX,F,TOL
      EXTERNAL      F
C
C  PURPOSE:
C     ZEROIN estimates a zero of the function F(X) in the interval AX,BX.
C
C  INPUT:
C
C  AX     Left endpoint of initial interval
C  BX     Right endpoint of initial interval
C  F      FUNCTION subprogram which evaluates F(X) for any X in
C         the interval AX,BX
C  TOL    Desired length of the interval of uncertainty of the
C         final result (TOL >= zero)
C
C  OUTPUT:
C
C  ZEROIN Abcissa approximating a zero of F in the interval AX,BX
C  IER    0 means no error detected;
C         1 means F(AX) and F(BX) have the same sign;
C         2 means AX >= BX.
C
C     ZEROIN originally assumed that F(AX) and F(BX) have opposite signs
C  without a check.  This version introduces the IER argument.
C
C     ZEROIN returns a zero X in the given interval (AX,BX) to within
C  a tolerance 4 * MACHEPS * ABS(X) + TOL, where MACHEPS is the relative
C  machine precision.
C
C     This FUNCTION subprogram is a slightly modified translation of
C  the ALGOL 60 procedure ZERO given in Richard Brent, Algorithms for
C  Minimization Without Derivatives, Prentice-Hall, Inc. (1973).
C
C-----------------------------------------------------------------------
C
      REAL  A,B,C,D,E,EPS,FA,FB,FC,TOL1,XM,P,Q,R,S
      REAL  ABS,SIGN
C
C  COMPUTE EPS, THE RELATIVE MACHINE PRECISION
C
      EPS = 1.0E+0
   10 EPS = EPS/2.0E+0
      TOL1 = 1.0E+0 + EPS
      IF (TOL1 .GT. 1.0E+0) GO TO 10
C
C INITIALIZATION
C
      IER = 2
      A = AX
      B = BX
      IF ( A .GE. B ) GO TO 99
C
      FA = F(A)
      FB = F(B)
C
      IER = 1
      IF ( FA .EQ. 0.E+0 .AND. FB .NE. 0.E+0 ) GO TO 15
      IF ( FB .EQ. 0.E+0 .AND. FA .NE. 0.E+0 ) GO TO 15
      IF ( FA * FB .GE. 0.E+0 ) GO TO 99
C
   15 IER = 0
C
C BEGIN STEP
C
   20 C = A
      FC = FA
      D = B - A
      E = D
   30 IF (ABS(FC) .GE. ABS(FB)) GO TO 40
      A = B
      B = C
      C = A
      FA = FB
      FB = FC
      FC = FA
C
C CONVERGENCE TEST
C
   40 TOL1 = 2.0E+0*EPS*ABS(B) + 0.5E+0*TOL
      XM = 0.5E+0*(C - B)
      IF (ABS(XM) .LE. TOL1) GO TO 90
      IF (FB .EQ. 0.0E+0) GO TO 90
C
C IS BISECTION NECESSARY
C
      IF (ABS(E) .LT. TOL1) GO TO 70
      IF (ABS(FA) .LE. ABS(FB)) GO TO 70
C
C IS QUADRATIC INTERPOLATION POSSIBLE
C
      IF (A .NE. C) GO TO 50
C
C LINEAR INTERPOLATION
C
      S = FB/FA
      P = 2.0E+0*XM*S
      Q = 1.0E+0 - S
      GO TO 60
C
C INVERSE QUADRATIC INTERPOLATION
C
   50 Q = FA/FC
      R = FB/FC
      S = FB/FA
      P = S*(2.0E+0*XM*Q*(Q - R) - (B - A)*(R - 1.0E+0))
      Q = (Q - 1.0E+0)*(R - 1.0E+0)*(S - 1.0E+0)
C
C ADJUST SIGNS
C
   60 IF (P .GT. 0.0E+0) Q = -Q
      P = ABS(P)
C
C IS INTERPOLATION ACCEPTABLE
C
      IF ((2.0E+0*P) .GE. (3.0E+0*XM*Q - ABS(TOL1*Q))) GO TO 70
      IF (P .GE. ABS(0.5E+0*E*Q)) GO TO 70
      E = D
      D = P/Q
      GO TO 80
C
C BISECTION
C
   70 D = XM
      E = D
C
C COMPLETE STEP
C
   80 A = B
      FA = FB
      IF (ABS(D) .GT. TOL1) B = B + D
      IF (ABS(D) .LE. TOL1) B = B + SIGN(TOL1, XM)
      FB = F(B)
      IF ((FB*(FC/ABS(FC))) .GT. 0.0E+0) GO TO 20
      GO TO 30
C
C DONE
C
   90 ZEROIN = B
   99 RETURN
      END
