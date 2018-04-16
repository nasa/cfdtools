C+----------------------------------------------------------------------
C
      SUBROUTINE LAPGRID ( IDIM, IMAX, JMAX, X, Y, NITERS, CONVRG )
C
C  PURPOSE:  LAPGRID uses the Laplace equation to compute an elliptic
C            solution to a 2-d grid.
C
C  ARGUMENTS:
C   NAME   DIM   TYPE   I/O/S   DESCRIPTION
C  IDIM     -     I       I     Declared row dimension of grid arrays.
C  IMAX,    -     I       I     Actual grid dimensions.
C  JMAX
C  X,Y  IMAX,JMAX R      I/O    Input :  Initial guesses for X and Y.
C                               Output:  The elliptic grid.
C  NITERS   -     I       I     Desired number of relaxation iterations.
C  CONVRG   -     R       O     Measure of convergence achieved: ratio
C                               of initial and final values of a quantity
C                               involving the square root of the sum of
C                               squares of dX and dY (between consecutive
C                               iterations) over all mesh points.
C  NOTES:
C   (1) IMPLICIT NONE is non-standard.
C
C  ENVIRONMENT:  VAX/VMS - FORTRAN 77
C
C  HISTORY:
C    Reese Sorenson, NASA Ames, 1980  Theory and practice (program GRAPE)
C    Jeff Cordova, Sterling, 1986     Initial adaptation as routine LAPGRID
C    John Melton/Ron Langhi, 1987     Argument list changes, etc.
C   
C-----------------------------------------------------------------------

C     Declarations
C     ------------

      IMPLICIT NONE

C ... Arguments:

      INTEGER
     >   IDIM, IMAX, JMAX, NITERS
      REAL
     >   CONVRG, X(IDIM,JMAX), Y(IDIM,JMAX)

C ... Local Constants:

      INTEGER
     >   MXWK
      REAL
     >   RELAX, ZERO, ONE, TWO, HALF
      PARAMETER
     >  (MXWK=500, ZERO=0.E+0, ONE=1.E+0, TWO=ONE+ONE, HALF=ONE/TWO,
     >   RELAX=0.9E+0)

C ... Local Variables:

      INTEGER
     >   I, J, N
      REAL
     >   A(MXWK), B(MXWK), C(MXWK), F1(MXWK), F2(MXWK),
     >   AL, BE, ELLRES, GA, RES1, RES2, RHSX, RHSY, XETA, YETA,
     >   XOLD, YOLD, XNRM, YNRM, XXI, YXI

C     Execution
C     ---------

      IF (JMAX .GT. MXWK) STOP 'LAPGRID: Not enough work-space.'

C ... Begin the iterative solution:

      DO 400, N = 1, NITERS

         XNRM = ZERO
         YNRM = ZERO

         DO 300, I = 2, IMAX-1

C ...       Boundary points:

            A(1)     = ZERO
            A(JMAX)  = ZERO
            B(1)     = ONE
            B(JMAX)  = ONE
            C(1)     = ZERO
            C(JMAX)  = ZERO      
            F1(1)    = X(I,1)
            F1(JMAX) = X(I,JMAX)
            F2(1)    = Y(I,1)
            F2(JMAX) = Y(I,JMAX)

C ...       Interior points:

            DO 100, J = 2, JMAX-1
               XXI  = X(I+1,J) - X(I-1,J)
               XETA = X(I,J+1) - X(I,J-1)
               YXI  = Y(I+1,J) - Y(I-1,J)
               YETA = Y(I,J+1) - Y(I,J-1)

               AL = XETA * XETA + YETA * YETA
               BE = (XXI * XETA  + YXI * YETA) * HALF
               GA = XXI * XXI   + YXI * YXI     
               A(J) = GA
               B(J) = -TWO * (AL + GA)
               C(J) = GA
               F1(J) = -AL * (X(I+1,J) + X(I-1,J)) + BE *
     >                 (X(I+1,J+1)-X(I-1,J+1)-X(I+1,J-1)+X(I-1,J-1))  
               F2(J) = -AL * (Y(I+1,J) + Y(I-1,J)) + BE *
     >                 (Y(I+1,J+1)-Y(I-1,J+1)-Y(I+1,J-1)+Y(I-1,J-1))  
 100        CONTINUE

C ...       Solve for both right hand sides together:

            CALL TRID2R ( JMAX, A, B, C, F1, F2 )

C ...       Under-relax coordinate values:

            IF (N .GT. 1 .AND. N .LT. NITERS) THEN

               DO 200, J = 2, JMAX-1 
                  X(I,J) = (ONE - RELAX) * X(I,J) + RELAX * F1(J)
                  Y(I,J) = (ONE - RELAX) * Y(I,J) + RELAX * F2(J)
  200          CONTINUE

            ELSE

C ...          Cheap check on convergence:

               DO 250, J = 2, JMAX-1 
                  XOLD = X(I,J)
                  YOLD = Y(I,J)
                  X(I,J) = (ONE - RELAX) * XOLD + RELAX * F1(J)
                  Y(I,J) = (ONE - RELAX) * YOLD + RELAX * F2(J)
                  XNRM = XNRM + (XOLD - X(I,J)) ** 2
                  YNRM = YNRM + (YOLD - Y(I,J)) ** 2
  250          CONTINUE
            END IF

  300    CONTINUE

         IF (N .EQ. 1) RES1 = SQRT (XNRM + YNRM)

 400  CONTINUE

      RES2 = SQRT (XNRM + YNRM)
      CONVRG = RES2 / RES1

      RETURN
      END
