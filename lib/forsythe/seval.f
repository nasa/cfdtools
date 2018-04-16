C+----------------------------------------------------------------------
C
      REAL FUNCTION SEVAL(N, U, X, Y, B, C, D)
      INTEGER N
      REAL  U, X(N), Y(N), B(N), C(N), D(N)
C
C  THIS SUBROUTINE EVALUATES THE CUBIC SPLINE FUNCTION
C
C    SEVAL = Y(I) + B(I)*(U-X(I)) + C(I)*(U-X(I))**2 + D(I)*(U-X(I))**3
C
C    WHERE  X(I) .LT. U .LT. X(I+1), USING HORNER'S RULE
C
C  IF  U .LT. X(1) THEN  I = 1  IS USED.
C  IF  U .GE. X(N) THEN  I = N  IS USED.
C
C  INPUT..
C
C    N = THE NUMBER OF DATA POINTS
C    U = THE ABSCISSA AT WHICH THE SPLINE IS TO BE EVALUATED
C    X,Y = THE ARRAYS OF DATA ABSCISSAS AND ORDINATES
C    B,C,D = ARRAYS OF SPLINE COEFFICIENTS COMPUTED BY SPLINE
C
C  IF  U  IS NOT IN THE SAME INTERVAL AS THE PREVIOUS CALL, THEN A
C  BINARY SEARCH IS PERFORMED TO DETERMINE THE PROPER INTERVAL.
C
C-----------------------------------------------------------------------
C
      INTEGER I, J, K
      REAL DX
      DATA I/1/
      IF ( I .GE. N ) I = 1
      IF ( U .LT. X(I) ) GO TO 10
      IF ( U .LE. X(I+1) ) GO TO 30
C
C  BINARY SEARCH
C
   10 I = 1
      J = N+1
   20 K = (I+J)/2
      IF ( U .LT. X(K) ) J = K
      IF ( U .GE. X(K) ) I = K
      IF ( J .GT. I+1 ) GO TO 20
C
C  EVALUATE SPLINE
C
   30 DX = U - X(I)
      SEVAL = Y(I) + DX*(B(I) + DX*(C(I) + DX*D(I)))
      RETURN
      END
