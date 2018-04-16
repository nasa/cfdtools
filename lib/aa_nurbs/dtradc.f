C**********************************************************************
      SUBROUTINE DTRADC (C,T,NWORK,WORK,RADCV,IER)
C***********************************************************************
C
C FUNCTION-
C          COMPUTE THE RADIUS OF CURVATURE AT A SPECIFIED PARAMETER
C          VALUE FOR A (RATIONAL) SPLINE PLANAR CURVE.
C
C AUTHORS-
C          C.P.CHI             CREATION  OCTOBER   1989
C          D.A.Saunders 4/93   NDIMV=3 (not 10) suffices for planar curves.
C                              Return early from a bad DTSPDR return.
C
C INPUT-
C           C       = DTRC SPLINE DATA C-ARRAY
C           T       = PARAMETER VALUE FOR WHICH THE RADIUS OF CURVATURE
C                     IS TO BE CALCULATED
C           NWORK   = LENGTH OF WORK(*) >= ORDER * (ORDER + 5)
C
C OUTPUT-
C           RADCV   = RADIUS OF CURVATURE
C           IER     = SUCCESS/ERROR CODE.
C                     IER=0     SUCCESS.
C                     IER=-1    C(3) .LE. 0.
C                     IER=-3    NWORK TOO SMALL.
C                     IER=-6    C(4) .LT. C(3).
C                     IER=-8    INVALID KNOT SET.
C                     IER=-10   DENOMINATOR .EQ. 0 ON A RATIONAL SPLINE.
C                     IER=-38   ATTEMPT TO EVALUATE AT POINT THAT IS
C                               INSIDE AN INTERVAL THAT IS TOO SMALL.
C                     IER=-50   X OUT OF RANGE.
C                     IER=-51   C(1) .NE. 1.
C                     IER=-52   C(2) .EQ. 0.
C                     IER=-60   SPLINE DATA IS NOT A PLANAR CURVE
C                     IER=-61   RADIUS OF CURVATURE IS INFINITE
C
C TYPE/PARAMETER/DIMENSION/COMMON/EQUIVALENCE/DATA/FORMAT STATEMENTS-

        PARAMETER (NDIMV=3, NDER=2)

        DOUBLE PRECISION  C(*), WORK(*), T, V(NDIMV,NDER+1)
        DOUBLE PRECISION  DX, DY, DXX, DYY, DPHI, DARC, RADCV
        INTEGER NVAR, NWORK, IER

C ********************** START OF EXECUTABLE CODE **********************

C...CHECK GENERAL DATA IN C-ARRAY

      CALL DTSCHK(C,IER)
      IF (IER.NE.0) GOTO 999

C...CHECK WHETHER THE SPLINE DATA IS FOR A PLANAR CURVE

      NVAR = NINT (C(2))
      IF (NVAR.GT.2 .OR. NVAR.LT.-3) THEN
        IER=-60
        GOTO 999
      ENDIF

C...CALCULATE THE DERIVATIVES AT T

      CALL DTSPDR(T, NDER, C, WORK, NWORK, V, NDIMV, IER)
      IF (IER .NE. 0) GO TO 999

C...CALCULATE THE RADIUS OF CURVATURE AT T.
C   VALUES IN V(*) ARE AS FOLLOWS:
C     THE FIRST SUBSCRIPT OF V IDENTIFIES THE DEPENDENT VARIABLE
C     THE SECOND SUBSCRIPT OF V IDENTIFIES THE DERIVATIVE, SUBSCRIPT
C     VALUE N LOCATING THE (N-1)THE DERIVATIVE.  THE ZEROTH DERIVATIVE IS
C     THE FUNCTION VALUE.

      DX = V(1,2)
      DY = V(2,2)
      DXX= V(1,3)
      DYY= V(2,3)

      IF (NVAR.EQ.1) THEN
        DPHI = DXX 
        DARC = SQRT(1.D0 + DX ** 2)
      ELSE
        DPHI = DX*DYY - DY*DXX
        DARC = SQRT(DX ** 2 + DY ** 2)
      ENDIF

      IF (DPHI.NE.0.D0) THEN
        RADCV = DARC ** 3 / DPHI
      ELSE
        IER=-61
      ENDIF

  999 RETURN
      END
