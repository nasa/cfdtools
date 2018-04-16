C**********************************************************************
      SUBROUTINE DTSREV(CRV, AXS, NDSF, CSF, ISIZE, IER)
C***********************************************************************
C
C FUNCTION-
C           GENERATE THE SPLINE ARRAY FOR A SURFACE OF REVOLUTION
C           BASED ON THE SPLINE OF A PLANAR CURVE (AND A UNIT CIRCLE).
C           THE FIRST (X) VARIABLE OF THE PLANAR CURVE GIVES POSITION
C           ALONG THE AXIS OF REVOLUTION AND THE SECOND (Y) VARIABLE 
C           OF THE PLANAR CURVE GIVES THE RADIUS OF THE SURFACE AT THAT
C           POSITION.
C
C AUTHORS-
C          C.P.CHI             CREATION  OCTOBER   1989
C
C INPUT-
C          CRV   = SPLINE ARRAY FOR THE PLANAR CURVE
C          AXS   = AXIS OF THE SURFACE OF REVOLUTION ('X', 'Y', OR 'Z')
C                  CONSIDER S AS THE INDEPENDENT VARIABLE FOR THE
C                  DIRECTION WHOSE CROSS SECTION IS ALWAYS A CIRCLE
C
C                  IF AXS = 'X'
C                     S = 0.0 IS IN THE XY - PLANE
C                  IF AXS = 'Y'
C                     S = 0.0 IS IN THE YZ - PLANE
C                  IF AXS = 'Z'
C                     S = 0.0 IS IN THE ZX - PLANE
C
C          NDSF  = LENGTH OF THE SPLINE ARRAY FOR THE SURFACE.
C                  MINIMUM REQUIRED VALUE FOR NDSF IS
C
C                    8+NKCIR+NKCRV+4*(NCP1*NCP2)
C
C                  NKCIR  = 12, LENGTH OF KNOT SEQUENCE FOR THE CIRCLE
C                  NKCRV  = LENGTH OF KNOT SEQUENCE FOR THE CURVE
C                  NCP1   = 9, NUMBER OF CONTROL POINTS FOR THE CIRCLE
C                  NCP2   = NUMBER OF CONTROL POINTS FOR THE CURVE
C
C OUTPUT-
C          CSF   = GENERATED SPLINE ARRAY FOR THE SURFACE OF REVOLUTION
C                  THE FIRST INDEPENDENT VARIABLE IS IN THE DIRECTION
C                  OF THE CIRCLE. THE SECOND INDEPENDENT VARIABLE IS
C                  IN THE DIRECTION OF THE CURVE. THE KNOT SEQUENCES
C                  AND COEFFICIENTS IN THE SPLINE ARRAY ARE ARRANGED
C                  ACCORDINGLY.
C          IER    = ERROR FLAG
C                 = 0   NO ERROR
C                 =-1   NUMBER OF INDEPENDENT VARIABLE NOT ACCEPTABLE
C                 =-2   NUMBER OF DEPENDENT VARIABLE NOT ACCEPTABLE
C                 =-3   ORDER OF THE SPLIN NOT ACCEPTABLE
C                 =-4   INCORRECT NUMBER OF KNOTS
C                 =-5   INCORRECT KNOT VECTOR
C                       MULTIPLICITY AT ENDS NOT CORRECT
C                       NOT IN ASCENDING ORDER
C                 =-66  SPLINE IS NOT FOR PLANAR CURVES
C                 =-67  LENGTH OF THE SPLINE ARRAY FOR THE SURFACE
C                       TOO SMALL
C
C REFERENCES-
C          1.
C
C NOTES-
C
C TYPE/PARAMETER/DIMENSION/COMMON/EQUIVALENCE/DATA/FORMAT STATEMENTS-
C
      CHARACTER AXS*1
C
      DOUBLE PRECISION CRV(*), CSF(NDSF)
      DOUBLE PRECISION CIRP(17),CIR1(9),CIR2(9),CIR3(9)
      DOUBLE PRECISION CCRV, CCIR, WCRV
C
      DATA (CIRP(I),I=1,17)
     $ /1.000D+00,          -4.000D+00,         3.000D+00,
     $   9.000D+00,         1.000D+00,          .000D+00,
     $  .000D+00,           .000D+00,           0.25000D+00,
     $  0.25000D+00,        0.50000D+00,        0.50000D+00,
     $  0.750000D+00,       0.750000D+00,       1.00000D+00,
     $  1.00000D+00,        1.00000D+00/
C
      DATA (CIR1(I),I=1,9)
     $ /1.000D+00,          .7071067D+00,        .000D+00,
     $  -.7071067D+00,      -1.000D+00,         -.7071067D+00,
     $   .000D+00,          .7071067D+00,       1.000D+00/
C
      DATA (CIR2(I),I=1,9)
     $ /.000D+00,          .7071067D+00,       1.000D+00,
     $ .7071067D+00,       .000D+00,           -.7071067D+00,
     $ -1.000D+00,         -.7071067D+00,      .000D+00/
C
      DATA (CIR3(I),I=1,9)
     $ /1.000D+00,          .7071067D+00,       1.000D+00,
     $  .7071067D+00,       1.000D+00,          .7071067D+00,
     $  1.000D+00,          .7071067D+00,       1.000D+00/
C
C ********************** START OF EXECUTABLE CODE **********************
C
      IER=0
C
C...INPUT SPLINES MUST BE PLANAR CURVES
C
      IF(CRV(1).NE.1.OR.(CRV(2).NE.2.AND.CRV(2).NE.-3)) GOTO 999
C
C...CHECK DATA OF THE INPUT SPLINES
C
C
      CALL DTSCHK(CRV, IER)
      IF(IER.NE.0) GOTO 999
C
C...CHECK REQUIRED ARRAY SIZE FOR THE SURFACE
C
      NKCIR=DINT(CIRP(3)+CIRP(4))
      NKCRV=DINT(CRV(3)+CRV(4))
      NCTRP=DINT(CIRP(4)*CRV(4))
C
C...GENERATE RATIONAL SPLINE SURFACE ONLY
C
      NCOEF=4
C
      ISIZE = 8+NKCIR+NKCRV+NCOEF*NCTRP
C
      IF(ISIZE.GT.NDSF) THEN
        IER=-67
        GOTO 999
      ENDIF
C
C...DEFINE THE GENERAL PARAMETER FOR THE SURFACE
C
      CSF(1) = 2.0D0
      CSF(2) = -4.0D0
      CSF(3) = CIRP(3)
      CSF(4) = CRV(3)
      CSF(5) = CIRP(4)
      CSF(6) = CRV(4)
      CSF(7) = CIRP(5)
      CSF(8) = CRV(5)
C
C...LOAD THE KNOT VECTORS
C
      IH1 = 5
      IKF1 = IH1+1
      IKL1 = IH1+NKCIR
      IK=8
      DO 100 I=IKF1, IKL1
      IK=IK+1
      CSF(IK) = CIRP(I)
  100 CONTINUE
C
      IH2 = 5
      IKF2= IH2+1
      IKL2= IH2+NKCRV
      DO 120 I=IKF2, IKL2
      IK=IK+1
      CSF(IK)=CRV(I)
  120 CONTINUE
C
C...NUMBER OF CONTROL POINTS
C
      NCP1 = DINT(CIRP(4))
      NCP2 = DINT(CRV(4))
C
C...SET UP POINTERS FOR THE COEF OF EACH INDEPENDENT VARIABLE
C   OF THE INPUT CURVE
C
      DO 130 N=1, 2
      IC = IKL2+(N-1)*NCP2
      IF(N.EQ.1) THEN
        IC1 = IC
      ELSE
        IC2 = IC
      ENDIF
  130 CONTINUE
C
C...GENERATE SPLINE COEF. FOR THE SURFACE
C
      ISF = IK
      NCXYZ = 3
C
C...LOOP FOR EACH COORDINATE
C
      DO 300 N=1, NCXYZ
C
C...THE COORDINATE OF DEPENDENT VARIABLES OF THE CURVE DEPEND ON
C   THE SPECIFIED AXIS OF THE SURFACE OF REVOLUTION.
C   THE COEF USED TO GENERATE THE SURFACE ARE TO BE SELECTED
C   ACCORDINGLY.
C
C              COORDINATE  COEF TO BE USED
C        AXIS          X          Y          Z
C      --------    --------------------------------
C         X         1ST COEF   2ND COEF   2ND COEF
C         Y         2ND COEF   1ST COEF   2ND COEF
C         Z         2ND COEF   2ND COEF   1ST COEF
C
      IF(AXS.EQ.'X') THEN
        IF(N.EQ.1) THEN
          ICA = IC1
        ELSEIF (N.EQ.2) THEN
          ICA = IC2
        ELSE
          ICA = IC2
        ENDIF
      ELSEIF(AXS.EQ.'Y') THEN
        IF(N.EQ.1) THEN
          ICA = IC2
        ELSEIF(N.EQ.2) THEN
          ICA = IC1
        ELSE
          ICA = IC2
        ENDIF
      ELSE
        IF(N.EQ.1) THEN
          ICA = IC2
        ELSEIF(N.EQ.2) THEN
          ICA = IC2
        ELSE
          ICA = IC1
        ENDIF
      ENDIF
C
      DO 200 I=1, NCP2
      ICA = ICA+1
      CCRV = CRV(ICA)
C
C...SELECT COEF OF THE CIRCLE TO BE USED TO GENERATE THE SURFACE
C
C                          COEF TO BE USED
C        AXIS          X          Y          Z
C      --------    --------------------------------
C         X           CIR3       CIR1       CIR2
C         Y           CIR2       CIR3       CIR1
C         Z           CIR1       CIR2       CIR3
C
      ICB = 0
      DO 150 J=1, NCP1
      ICB = ICB+1
C
      IF(AXS.EQ.'X') THEN
        IF(N.EQ.1) THEN
          CCIR = CIR3(ICB)
        ELSEIF(N.EQ.2) THEN
          CCIR = CIR1(ICB)
        ELSE
          CCIR = CIR2(ICB)
        ENDIF
      ELSEIF(AXS.EQ.'Y') THEN
        IF(N.EQ.1) THEN
          CCIR = CIR2(ICB)
        ELSEIF(N.EQ.2) THEN
          CCIR = CIR3(ICB)
        ELSE
          CCIR = CIR1(ICB)
        ENDIF
      ELSE
        IF(N.EQ.1) THEN
          CCIR = CIR1(ICB)
        ELSEIF(N.EQ.2) THEN
          CCIR = CIR2(ICB)
        ELSE
          CCIR = CIR3(ICB)
        ENDIF
      ENDIF
C
C...COEF FOR THE SURFACE
C
      ISF = ISF+1
      CSF(ISF) = CCRV*CCIR
C
  150 CONTINUE
  200 CONTINUE
  300 CONTINUE
C
C...PRODUCT OF THE WEIGHTS
C
      IB = IC2+NCP2
      DO 400 I=1,NCP2
      IB=IB+1
      IF(CRV(2).LT.0) THEN
        WCRV = CRV(IB)
      ELSE
        WCRV = 1.0D+00
      ENDIF
C
      IA = 0
      DO 350 J=1,NCP1
      IA=IA+1
      ISF=ISF+1
      CSF(ISF)=CIR3(IA)*WCRV
  350 CONTINUE
  400 CONTINUE
C
  999 CONTINUE
      RETURN
      END
