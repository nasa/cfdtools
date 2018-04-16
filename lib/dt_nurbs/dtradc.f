C**********************************************************************
      SUBROUTINE DTRADC(C,AX,NWORK,WORK,RADCV,IER)
C***********************************************************************
C
C FUNCTION-
C          COMPUTE THE RADIUS OF CURVATURE AT A SPECIFIED PARAMETER
C          VALUE FOR A (RATIONAL) SPLINE PLANAR CURVE.
C
C AUTHORS-
C          C.P.CHI             CREATION  OCTOBER   1989
C
C INPUT-
C           C       = DTRC SPLINE DATA C-ARRAY
C           AX      = PARAMETER VALUE FOR WHICH THE RADIUS OF CURVATURE
C                     IS TO BE CALCULATED
C           NWORK   = DIMENSION OF THE WORKING ARRAY WORK
C
C OUTPUT-
C           RADCV   = RADIUS OF CURVATURE
C           IER     = SUCCESS/ERROR CODE.
C                     IER=0     SUCCESS.
C                     IER=-1    C(3) .LE. 0.
C                     IER=-3    NWORK TOO SMALL.
C                     IER=-4    NDIMV .LT. C(2)
C                     IER=-6    C(4) .LT. C(3).
C                     IER=-8    INVALID KNOT SET.
C                     IER=-10   DENOMINATOR .EQ. 0 ON A RATIONAL SPLINE.
C                     IER=-38   ATTEMPT TO EVALUATE AT POINT THAT IS
C                               INSIDE AN INTERVAL THAT IS TOO SMALL.
C                     IER=-50   X OUT OF RANGE.
C                     IER=-51   C(1) .NE. 1.
C                     IER=-52   C(2) .EQ. 0.
C                     IER=-53   NUMBER OF DERIVATIVES LESS THAN 0
C                     IER=-60   SPLINE DATA IS NOT A PLANAR CURVE
C                     IER=-61   RADIUS OF CURVATURE IS INFINITE
C
C TYPE/PARAMETER/DIMENSION/COMMON/EQUIVALENCE/DATA/FORMAT STATEMENTS-
C
        PARAMETER (NDIMV=10, NDER=2)
C
        DOUBLE PRECISION  C(*), WORK(*), AX, V(NDIMV,NDER+1)
        DOUBLE PRECISION  DX, DY, DXX, DYY, DPHI, DARC, RADCV
        INTEGER NWORK, IER
C
C ********************** START OF EXECUTABLE CODE **********************
C
C...CHECK GENERAL DATA IN C-ARRAY
C
      CALL DTSCHK(C,IER)
      IF(IER.NE.0) GOTO 999
C
C...CHECK WHETHER THE SPLINE DATA IS FOR A PLANAR CURVE
C
      IF((C(2).GT.2).OR.(C(2).LT.-3)) THEN
        IER=-60
        GOTO 999
      ENDIF
C
C...CALCULATE THE DERIVATIVES AT AX
C
      CALL DTSPDR(AX, NDER, C, WORK, NWORK, V, NDIMV, IER)
C
C...CALCULATE THE RADIUS OF CURVATURE AT AX.
C   VALUES IN THE V - ARRAY ARE AS FOLLOWS:
C     THE FIRST SUBSCRIPT OF V IDENTIFIES THE DEPENDENT VARIABLE
C     THE SECOND SUBSCRIPT OF V IDENTIFIES THE DERIVATIVE, SUBSCRIPT
C     VALUE N LOCATING THE (N-1)THE DERIVATIVE.  THE ZEROTH DERIVATIVE IS
C     THE FUNCTION VALUE
C
      DX = V(1,2)
      DY = V(2,2)
      DXX= V(1,3)
      DYY= V(2,3)
C
      IF(C(2).EQ.1) THEN
        DPHI=DXX 
        DARC=DSQRT(1.D0 + DX ** 2)
      ELSE
        DPHI = DX*DYY - DY*DXX
        DARC = DSQRT(DX ** 2 + DY ** 2)
      ENDIF
      IF(DPHI.EQ.0.D0) THEN
        IER=-61
        GOTO 999
      ENDIF
C
      RADCV = DARC ** 3 / DPHI
C
  999 CONTINUE
      RETURN
      END
