C+----------------------------------------------------------------------
C
      REAL FUNCTION AVINT (X, Y, N, XLO, XUP)
C
C     PURPOSE:
C        This function subprogram computes the approximate integral of
C     Y (represented by discrete points) with respect to X between the
C     limits XLO, XUP.  X and Y are N-dimensional arrays .  The X(I)
C     may be unequally spaced but must be given in ascending order.
C     XLO must be less than XUP, but otherwise there are no restrictions
C     on XLO and XUP except that for greatest accuracy the limits of
C     integration should not be too far from the first and last data
C     points.  The integral over each subinterval is computed by
C     averaging quadratic functions fit from left and right neighboring
C     intervals, with special treatment for the first and last
C     subintervals.
C 
C        This version of AVINT was taken almost verbatim from NUMERICAL
C     INTEGRATION by Davis and Rabinowitz, pg. 193.  They in turn
C     adapted it from Algorithm 77, CACM 5, 1962, pg. 96 by Hennion.
C
C        Robert Kennelly, Spring 1978 and Fall 1981.
C
C-----------------------------------------------------------------------


      DIMENSION X(N), Y(N)

      IB = 2
      DO 10 I = 1, N
         IF (X(I) .GE. XLO) GO TO 20
         IB = IB + 1
   10 CONTINUE
   20 CONTINUE

      IT = N
      DO 30 I = 1, N
         IF (X(IT) .LE. XUP) GO TO 40
         IT = IT - 1
   30 CONTINUE
   40 CONTINUE

      IT = IT - 1
      SUM = 0.E+0
      SYL = XLO
      DO 50 JM = IB, IT
         X1 = X(JM - 1)
         X2 = X(JM)
         X3 = X(JM + 1)
         TERM1 = Y(JM - 1)/((X1 - X2)*(X1 - X3))
         TERM2 = Y(JM)/((X2 - X1)*(X2 - X3))
         TERM3=Y(JM + 1)/((X3 - X1)*(X3 - X2))
         A = TERM1 + TERM2 + TERM3
         B = -(X2 + X3)*TERM1 - (X1 + X3)*TERM2 - (X1 + X2)*TERM3
         C = X2*X3*TERM1 + X1*X3*TERM2 + X1*X2*TERM3
         IF (JM .NE. IB) GO TO 60
            CA = A
            CB = B
            CC = C
            GO TO 70
   60    CONTINUE
            CA = .5E+0*(A + CA)
            CB = .5E+0*(B + CB)
            CC = .5E+0*(C + CC)
   70    CONTINUE
         SYU = X(JM)
         SUM = SUM + CA*(SYU**3 - SYL**3)/3.E+0 +
     >      CB*.5E+0*(SYU**2 - SYL**2) + CC*(SYU - SYL)
         CA = A
         CB = B
         CC = C
         SYL = SYU
   50 CONTINUE

      AVINT = SUM + CA*(XUP**3 - SYL**3)/3.E+0 +
     >   CB*.5E+0*(XUP**2 - SYL**2) + CC*(XUP - SYL)

      RETURN
      END
