      SUBROUTINE DTCMAT(X,N,XKNTS,NCOEF,KORD,MLT,NLHS,NRHS,
     +                  NDIM,A,WORK,BSPL)
C
C=======================================================================
C
C  PURPOSE:
C          DTCMAT constructs the coefficient matrix needed in solving
C          for the b-spline coefficients.  This is a low level 
C          routine which assumes that argument checking has occured
C          in the calling routine.
C
C  METHOD:
C          The b-spline functions are evaluated at each of the data
C          points, X.  An array is formed (in packed form) with each
C          row corresponding to a data interpolation constraint, and
C          each column the value of a b-spline function at that
C          point.
C
C  USAGE:
C          DOUBLE PRECISION X(N),XKNTS(*),A(NDIM,*),WORK(*),
C                           BSPL(KORD,KORD)
C          CALL DTCMAT(X,N,XKNTS,NCOEF,KORD,MLT,NLHS,NRHS,
C                      NDIM,A,WORK,BSPL)
C
C  INPUT:
C          N       Number of data points.
C
C          X       Array of N values of the independent variable in
C                  ascending order.
C
C          XKNTS   Array of knots.
C
C          NCOEF   Number of b-spline coefficients to solve for, 
C                  i.e. size of matrix.
C
C          KORD    Order of spline = DEGREE + 1.
C
C          MLT     Multiplicity of data at all points.
C
C          NLHS    Number of additional left endpoint conditions.
C
C          NRHS    Number of additional right endpoint conditions.
C
C          NDIM    First dimension of A array.
C
C  STORAGE:
C          WORK    Work array needed by DTBSP1.
C
C          BSPL    Work array for holding computation of b-splines
C                  at data points.
C
C
C  OTHER DT ROUTINES CALLED
C
C       DTBSP1
C
C     ******************************************************************
C
      DOUBLE PRECISION ZERO
      PARAMETER  (ZERO=0.D0)
C
C  Arguments:
C
      INTEGER  N, NCOEF, KORD, MLT, NLHS, NRHS, NDIM, NDER
      DOUBLE PRECISION  X(N), XKNTS(*), A(NDIM,*),
     &                  WORK(*), BSPL(KORD,*)
C
C Internal:
C
      INTEGER  I,IPOS,IR,J,JC,JR,KEVEN,KHALF,L,M,NM2,NLPT,NRPT
C
C     ******************************************************************
C
      KHALF = (KORD+1)/2
      KEVEN = 2*KHALF
      M = 2*(KORD-1) + 1
      NM2 = N - 2
      DO 20 J = 1, NCOEF
        DO 10 I = 1, NDIM
          A(I,J) = ZERO
   10   CONTINUE
   20 CONTINUE
C
C=======================================================================
C    Evaluate B-splines on left hand side.
C=======================================================================
C
      NLPT = MLT + NLHS
      IPOS = KORD
      NDER = KORD - NLPT
      CALL DTBSP1(XKNTS,X(1),IPOS,KORD,NDER,WORK,BSPL,KORD)
C
C=======================================================================
C    Define first KHALF rows of A.
C=======================================================================
C
      IR = 0
      DO 70 I=KEVEN-NLPT,KHALF+1,-1
        IR = IR + 1
        DO 60 J=1,KORD
          JR = IR - J + M
          A(JR,J) = BSPL(J,I)
   60   CONTINUE
   70 CONTINUE
      DO 90 I=1,NLPT
        IR = IR + 1
        DO 85 J =1,KORD
          JR = IR - J + M
          A(JR,J) = BSPL(J,I)
   85   CONTINUE
   90 CONTINUE
C
C=======================================================================
C    Define the interior rows of A to represent the interpolation.
C=======================================================================
C
      IF (MLT.EQ.0) GO TO 150
      DO 130 I = 1, NM2
        NDER = MLT - 1
        IPOS = IPOS + MLT
        CALL DTBSP1(XKNTS,X(I+1),IPOS,KORD,NDER,WORK,BSPL,KORD)
        DO 120 L=1,MLT
          IR = IR + 1
          JC = IPOS - KORD
          DO 110 J = 1, KORD
            JC = JC + 1
            JR = IR - JC + M
            A(JR,JC) = BSPL(J,L)
110       CONTINUE
120     CONTINUE
130   CONTINUE
C
C=======================================================================
C    Evaluate B-splines on right hand side.
C=======================================================================
C
150   CONTINUE
      NRPT = MLT + NRHS
      NDER = KORD - NRPT
      CALL DTBSP1(XKNTS,X(N),NCOEF,KORD,NDER,WORK,BSPL,KORD)
C
C=======================================================================
C    Now set last KHALF rows of A.
C=======================================================================
C
      DO 170 I=1,NRPT
        JC = NCOEF - KORD
        IR = IR + 1
        DO 160 J=1,KORD
          JC = JC + 1
          JR = IR - JC + M
          A(JR,JC) = BSPL(J,I)
 160    CONTINUE
 170  CONTINUE
      DO 190 I=KHALF+1,KEVEN-NRPT
        JC = NCOEF - KORD
        IR = IR + 1
        DO 185 J=1,KORD
          JC = JC + 1
          JR = IR - JC + M
          A(JR,JC) = BSPL(J,I)
 185    CONTINUE
 190  CONTINUE
C
      RETURN
      END
