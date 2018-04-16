C+----------------------------------------------------------------------
C
C     NAME:  AKIMA (source module) or
C            IQHSCV (object module).
C            (See IQHSCV calling sequence below under USAGE.)
C
C     PURPOSE:  Two-dimensional data interpolation (fit plus evaluation)
C               for irregular data via triangularization and 2D spline
C               algorithms. (Adapted from IMSL's IQHSCV.)
C
C     REVISIONS:  
C
C     January 31, 1983 DLR
C        Now may input irregularly spaced data points at which to evaluate
C        the interpolant. Number of data points changed to be four or more.
C
C     March 17, 1983 DLR
C        Added parameter to turn extrapolation off or on.
C
C     March 11, 1986 DAS
C        All ZI(*) are set to 9998. at the outset because the package has
C        been observed to fail to set ZI(*) for an out-of-range point, even
C        if extrapolation is permitted.  This depended on the order of the
C        XI,YI points - evidently OK if the out-of-range point is first or
C        last, but not if it is somewhere else (???).
C
C        9999. is substituted where extrapolation would have been used if
C        it had not been turned off with IEX.NE.0.  (This is considered a
C        non-fatal error and is consistent with certain applications using
C        9999. for data points which are to be suppressed.)
C
C        9997. is substituted upon failure in computing the tension term
C        from a 3x3 system - yet to be seen.
C
C     ENVIRONMENT:  VAX/VMS FORTRAN
C
C     USAGE:   CALL IQHSCV (XD,YD,ZD,ND,XI,YI,ZI,NPI,IWK,WK,IER,IEX,TENS)
C
C     ARGUMENTS:
C        XD     - Vector of length ND. (INPUT)
C        YD     - Vector of length ND. (INPUT)
C        ZD     - Vector of length ND. (INPUT)
C                 ZD(I) is the function value at (XD(I),YD(I)).
C        ND     - Number of data points. ND must be >= 4. (INPUT)
C        XI     - Vector of length NPI. (INPUT)
C        YI     - Vector of length NPI. (INPUT)
C        ZI     - Vector containing the interpolated values. (OUTPUT)
C                 ZI(I) is the estimated value at the point (XI(I),YI(I)).
C                 See Revision history above for the meaning of the
C                 values 9999., 9998., 9997. in ZI(I).
C        NPI    - Number of points at which to interpolate. (INPUT)
C                 NPI must be at least 1.
C        IWK    - Integer work vector of length 31*ND+NPI.
C        WK     - Real work vector of length 6*ND.
C        IER    - Error parameter. (OUTPUT)
C                 IER = 0:   No error was detected.
C                 IER = 129: ND < 4, or NPI < 1.
C                 IER = 130: All data points are collinear.
C                 IER = 131: Some data points are identical.
C                 IER = 999: Bad value(s) 9997. or 9998. returned.
C        IEX    - Extrapolation parameter. IEX .NE. 0 signifies no 
C                 extrapolation. (Should have been other way round,
C                 but can't change it without affecting users.)
C        TENS   - Tension parameter.  Interpolated function value is
C                 TENS * linear term + (1-TENS) * nonlinear term.
C
C   REQD. ROUTINES:  IQHSD,IQHSE,IQHSF,IQHSG,IQHSH
C
C-----------------------------------------------------------------------
C
      SUBROUTINE IQHSCV (XD,YD,ZD,ND,XI,YI,ZI,NPI,IWK,WK,IER,IEX,TENS)
C
C                                  SPECIFICATIONS FOR ARGUMENTS
C
C
      INTEGER    ND,NPI,IER
      DIMENSION  IWK(1), WK(1)
      DIMENSION  XD(1),YD(1),ZD(1)
      DIMENSION  XI(1),YI(1)
      DIMENSION  ZI(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IL1,IL2,ITI,IZ,JIG0MN,JIG0MX,JIG1MN,
     1                   JIG1MX,JIGP,JNGP,JWIGP0,JWIGP,JWIPL,JWIPT,
     2                   JWIWL,JWIWP,JWNGP0,JWNGP,JWWPD,NDP0,NGP0,NGP1,
     3                   NL,NNGP,NT,NPI0
      INTEGER            ITPV
      REAL               DMMY(27)
      COMMON /IBCDPT/    DMMY,ITPV
C
C                                  SETTING OF SOME INPUT PARAMETERS TO
C                                    LOCAL VARIABLES.
C                                  FIRST EXECUTABLE STATEMENT
      LUNERR = 6
      IER = 0
C
C     The following is a way around an observed problem with out-of-range
C     points - see revision history above.
C
      DO 3 IZ = 1, NPI
         ZI(IZ) = 9998.
3     CONTINUE
C
      NDP0 = ND
      NPI0 = NPI
C
      IF (NDP0.LT.4 .OR. NPI0.LT.1) THEN
         IER = 129
         GO TO 9000
      END IF

      IWK(1) = NDP0
      IWK(3) = NPI0
C                                  ALLOCATION OF STORAGE AREAS IN THE
C                                    IWK ARRAY.
      JWIPT = 16
      JWIWL = 6*NDP0+1             
      JWNGP0 = JWIWL-1             
      JWIPL = 24*NDP0+1            
      JWIWP = 30*NDP0+1            
      JWIGP0 = 31*NDP0             
      JWWPD = 5*NDP0+1             
      CALL IQHSG (NDP0,XD,YD,NT,IWK(JWIPT),NL,IWK(JWIPL),IWK(JWIWL),
     1            IWK(JWIWP),WK,IER)
      IF (IER.GE.128) GO TO 9000
      IWK(5) = NT
      IWK(6) = NL
      IF (NT.EQ.0) GO TO 9005
C                                  SORTS OUTPUT GRID POINTS IN
C                                    ASCENDING ORDER OF THE TRIANGLE
C                                    NUMBER AND THE BORDER LINE SEGMENT
C                                    NUMBER.
      CALL IQHSH (XD,YD,NT,IWK(JWIPT),NL,IWK(JWIPL),NPI0,XI,YI,
     1 IWK(JWNGP0+1),IWK(JWIGP0+1))
C
C                                  ESTIMATES PARTIAL DERIVATIVES AT ALL
C                                    DATA POINTS.
C
C
      CALL IQHSE (NDP0,XD,YD,ZD,NT,IWK(JWIPT),WK,WK(JWWPD))
C
C                                  INTERPOLATES THE ZI VALUES.
C
      ITPV = 0
      JIG0MX = 0
      JIG1MN = NPI0+1
      NNGP = NT+2*NL
      DO 25 JNGP=1,NNGP
         ITI = JNGP
         IF (JNGP.LE.NT) GO TO 5
         IL1 = (JNGP-NT+1)/2
         IL2 = (JNGP-NT+2)/2
         IF (IL2.GT.NL) IL2 = 1
         ITI = IL1*(NT+NL)+IL2
    5    JWNGP = JWNGP0+JNGP
         NGP0 = IWK(JWNGP)
         IF (NGP0.EQ.0) GO TO 15
         JIG0MN = JIG0MX+1
         JIG0MX = JIG0MX+NGP0
         DO 10 JIGP=JIG0MN,JIG0MX
            JWIGP = JWIGP0+JIGP
            IZ = IWK(JWIGP)
            VXI=XI(IZ)
            VYI=YI(IZ)
            VIRZI=ZI(IZ)
            CALL IQHSF (XD,YD,ZD,NT,IWK(JWIPT),NL,IWK(JWIPL),WK,ITI,
     1      VXI,VYI,VIRZI,IEX,TENS)
            ZPRIME=VIRZI
            ZI(IZ) = ZPRIME
   10    CONTINUE
   15    JWNGP = JWNGP0+2*NNGP+1-JNGP
         NGP1 = IWK(JWNGP)
         IF (NGP1.EQ.0) GO TO 25
         JIG1MX = JIG1MN-1
         JIG1MN = JIG1MN-NGP1
         DO 20 JIGP=JIG1MN,JIG1MX
            JWIGP = JWIGP0+JIGP
            IZ = IWK(JWIGP)
            VIRZI=ZI(IZ)
            VIRXI=XI(IZ)
            VIRYI=YI(IZ)
            CALL IQHSF (XD,YD,ZD,NT,IWK(JWIPT),NL,IWK(JWIPL),WK,ITI,
     1      VIRXI,VIRYI,VIRZI,IEX,TENS)
            YI(IZ)=VIRYI
            XI(IZ)=VIRXI
            ZPRIME=VIRZI
            ZI(IZ) = ZPRIME
   20    CONTINUE
   25 CONTINUE
C
C****************************************************************
C
      GO TO 9005
C                                  ERROR EXIT
 9000 CONTINUE
      IF (IER.EQ.129) THEN
         WRITE (LUNERR,1001) 'less than 4 data points or 1 interp. pt.'
      ELSE IF (IER.EQ.130) THEN
         WRITE (LUNERR,1001) 'all data points are collinear.'
      ELSE IF (IER.EQ.131) THEN
         WRITE (LUNERR,1001) 'data points are not all distinct.'
      ELSE
         WRITE (LUNERR,1001) 'undefined error in IQHSCV. IER =', IER
      END IF
 9005 CONTINUE
C
C     Look for possible algorithm failure (as opposed to suppression
C     of extrapolation, where 9999. is used):
C
      NT = 0
      DO 9010 IZ = 1, NPI
         IF ( ZI(IZ) .EQ. 9998.  .OR.  ZI(IZ) .EQ. 9997. ) NT = NT+1
 9010 CONTINUE
C
      IF ( NT.GT.0 ) THEN
         WRITE (LUNERR,1001) 'number of pts. not interpolated =', NT
         IER = 999
      END IF
C
      RETURN
 1001 FORMAT ( '0IQHSCV failure: ', A, I6 )
      END
C-----------------------------------------------------------------------
C
C   COMPUTER            - DEC11/SINGLE
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY  SUBROUTINE IQHSCV
C
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - SINGLE/H36,H48,H60
C
C   REQD. ROUTINES - NONE                                          ISYD0
C
C-----------------------------------------------------------------------
C
      INTEGER FUNCTION IQHSD (X,Y,I1,I2,I3,I4)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            I1,I2,I3,I4
      DIMENSION            X(1),Y(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IDX
      REAL               A1SQ,A2SQ,A3SQ,A4SQ,B1SQ,B2SQ,B3SQ,B4SQ,C1SQ,
     1                   C2SQ,C3SQ,C4SQ,EPSLN,S1SQ,S2SQ,S3SQ,S4SQ,TOL,
     2                   U1,U2,U3,U4,X1,X2,X3,X4,Y1,Y2,Y3,Y4
      EQUIVALENCE        (C2SQ,C1SQ),(A3SQ,B2SQ),(B3SQ,A1SQ),
     1                   (A4SQ,B1SQ),(B4SQ,A2SQ),(C4SQ,C3SQ)
C                                  MACHINE PRECISION
      DATA               TOL/1.1921E-07/
C                                  FIRST EXECUTABLE STATEMENT
      EPSLN = TOL*100.0
C                                  PRELIMINARY PROCESSING
      X1 = X(I1)
      Y1 = Y(I1)
      X2 = X(I2)
      Y2 = Y(I2)
      X3 = X(I3)
      Y3 = Y(I3)
      X4 = X(I4)
      Y4 = Y(I4)
C                                  CALCULATION
      IDX = 0
      U3 = (Y2-Y3)*(X1-X3)-(X2-X3)*(Y1-Y3)
      U4 = (Y1-Y4)*(X2-X4)-(X1-X4)*(Y2-Y4)
      IF (U3*U4.LE.0.0) GO TO 5
      U1 = (Y3-Y1)*(X4-X1)-(X3-X1)*(Y4-Y1)
      U2 = (Y4-Y2)*(X3-X2)-(X4-X2)*(Y3-Y2)
      A1SQ = (X1-X3)**2+(Y1-Y3)**2
      B1SQ = (X4-X1)**2+(Y4-Y1)**2
      C1SQ = (X3-X4)**2+(Y3-Y4)**2
      A2SQ = (X2-X4)**2+(Y2-Y4)**2
      B2SQ = (X3-X2)**2+(Y3-Y2)**2
      C3SQ = (X2-X1)**2+(Y2-Y1)**2
      S1SQ = U1*U1/(C1SQ*AMAX1(A1SQ,B1SQ))
      S2SQ = U2*U2/(C2SQ*AMAX1(A2SQ,B2SQ))
      S3SQ = U3*U3/(C3SQ*AMAX1(A3SQ,B3SQ))
      S4SQ = U4*U4/(C4SQ*AMAX1(A4SQ,B4SQ))
      IF ((AMIN1(S3SQ,S4SQ)-AMIN1(S1SQ,S2SQ)).GT.EPSLN) IDX = 1
    5 IQHSD = IDX
      RETURN
      END
C
C     Old IMSL partial derivative estimator
C     (Old? What does this mean?)
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - DEC11/SINGLE
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY SUBROUTINE IQHSCV
C
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - SINGLE/H36,H48,H60
C
C   REQD. ROUTINES - NONE                                          ISYE0
C
C-----------------------------------------------------------------------
C
      SUBROUTINE IQHSE  (NDP,XD,YD,ZD,NT,IPT,PD,WK)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NDP,NT,IPT(1)
      DIMENSION    XD(1),YD(1),ZD(1)
      REAL         PD(1),WK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IDP,IPTI(3),IT,IV,JPD0,JPDMX,JPD,JPT0,JPT,NDP0,
     1                   NT0
      REAL               DX1,DX2,DY1,DY2,DZ1,DZ2,DZX1,DZX2,DZY1,DZY2,
     1                   EPSLN,TOL,XV(3),VPXX,VPXY,VPX,VPYX,VPYY,VPY,
     2                   VPZMN,VPZ,YV(3),ZV(3),ZXV(3),ZYV(3)
C                                  MACHINE PRECISION
      DATA               TOL/1.1921E-07/
C                                  FIRST EXECUTABLE STATEMENT
      EPSLN = TOL*100.0
C                                  PRELIMINARY PROCESSING
      NDP0 = NDP
      NT0 = NT
C                                  CLEARS THE PD ARRAY.
      JPDMX = 5*NDP0
      DO 5 JPD=1,JPDMX
         PD(JPD) = 0.0
    5 CONTINUE
      DO 10 IDP=1,NDP
         WK(IDP) = 0.0
   10 CONTINUE
C                                  ESTIMATES ZX AND ZY.
      DO 25 IT=1,NT0
         JPT0 = 3*(IT-1)
         DO 15 IV=1,3
            JPT = JPT0+IV
            IDP = IPT(JPT)
            IPTI(IV) = IDP
            XV(IV) = XD(IDP)
            YV(IV) = YD(IDP)
            ZV(IV) = ZD(IDP)
   15    CONTINUE
         DX1 = XV(2)-XV(1)
         DY1 = YV(2)-YV(1)
         DZ1 = ZV(2)-ZV(1)
         DX2 = XV(3)-XV(1)
         DY2 = YV(3)-YV(1)
         DZ2 = ZV(3)-ZV(1)
         VPX = DY1*DZ2-DZ1*DY2
         VPY = DZ1*DX2-DX1*DZ2
         VPZ = DX1*DY2-DY1*DX2
         VPZMN = ABS(DX1*DX2+DY1*DY2)*EPSLN
         IF (ABS(VPZ).LE.VPZMN) GO TO 25
         DO 20 IV=1,3
            IDP = IPTI(IV)
            JPD0 = 5*(IDP-1)+1
            PD(JPD0) = PD(JPD0)+VPX
            PD(JPD0+1) = PD(JPD0+1)+VPY
            WK(IDP) = WK(IDP)+VPZ
   20    CONTINUE
   25 CONTINUE
      DO 30 IDP=1,NDP0
         JPD0 = 5*(IDP-1)+1
         PD(JPD0) = -PD(JPD0)/WK(IDP)
         PD(JPD0+1) = -PD(JPD0+1)/WK(IDP)
   30 CONTINUE
C                                  ESTIMATES ZXX, ZXY, AND ZYY.
      DO 45 IT=1,NT0
         JPT0 = 3*(IT-1)
         DO 35 IV=1,3
            JPT = JPT0+IV
            IDP = IPT(JPT)
            IPTI(IV) = IDP
            XV(IV) = XD(IDP)
            YV(IV) = YD(IDP)
            JPD0 = 5*(IDP-1)+1
            ZXV(IV) = PD(JPD0)
            ZYV(IV) = PD(JPD0+1)
   35    CONTINUE
         DX1 = XV(2)-XV(1)
         DY1 = YV(2)-YV(1)
         DZX1 = ZXV(2)-ZXV(1)
         DZY1 = ZYV(2)-ZYV(1)
         DX2 = XV(3)-XV(1)
         DY2 = YV(3)-YV(1)
         DZX2 = ZXV(3)-ZXV(1)
         DZY2 = ZYV(3)-ZYV(1)
         VPXX = DY1*DZX2-DZX1*DY2
         VPXY = DZX1*DX2-DX1*DZX2
         VPYX = DY1*DZY2-DZY1*DY2
         VPYY = DZY1*DX2-DX1*DZY2
         VPZ = DX1*DY2-DY1*DX2
         VPZMN = ABS(DX1*DX2+DY1*DY2)*EPSLN
         IF (ABS(VPZ).LE.VPZMN) GO TO 45
         DO 40 IV=1,3
            IDP = IPTI(IV)
            JPD0 = 5*(IDP-1)+3
            PD(JPD0) = PD(JPD0)+VPXX
            PD(JPD0+1) = PD(JPD0+1)+VPXY+VPYX
            PD(JPD0+2) = PD(JPD0+2)+VPYY
   40    CONTINUE
   45 CONTINUE
      DO 50 IDP=1,NDP0
         JPD0 = 5*(IDP-1)+3
         PD(JPD0) = -PD(JPD0)/WK(IDP)
         PD(JPD0+1) = -PD(JPD0+1)/(2.0*WK(IDP))
         PD(JPD0+2) = -PD(JPD0+2)/WK(IDP)
   50 CONTINUE
      RETURN
      END
C
C   ROUTINE NAME   - IQHSF                                         ISYF0
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - DEC11/SINGLE
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY SUBROUTINE IQHSCV
C
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - SINGLE/H36,H48,H60
C
C   REQD. ROUTINES - NONE                                          ISYF0
C
C-----------------------------------------------------------------------
C
      SUBROUTINE IQHSF  (XD,YD,ZD,NT,IPT,NL,IPL,PDD,ITI,XII,YII,ZII,IEX,
     1                   TENS)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NT,NL,ITI,IPT(1),IPL(1)
      DIMENSION            XD(1),YD(1),ZD(1)
      REAL               CMAT(3,3),RHS(3),WKAREA(3),IPVT(3)
      REAL               PDD(1),XII,YII,ZII
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IDP,IL1,IL2,IT0,I,JIPL,JIPT,JPDD,JPD,KPD,NTL
      INTEGER            ITPV
      REAL               X0,Y0,AP,BP,CP,DP,P00,P10,P20,P30,P40,P50,P5,
     1                   P01,P11,P21,P31,P41,P02,P12,P22,P32,P03,P13,
     2                   P23,P04,P14,P05
      REAL               AA,AB,ACT2,AC,ADBC,AD,A,BB,BC,BDT2,B,CC,CD,
     1                   CSUV,C,DD,DLT,DX,DY,D,G1,G2,H1,H2,H3,LU,
     2                   LV,X(3),P0,P1,P2,P3,P4,PD(15),THSV,THUS,THUV,
     3                   THXU,U,V,Y(3),Z0,ZUU(3),ZUV(3),ZU(3),ZVV(3),
     4                   ZV(3),Z(3)
      COMMON /IBCDPT/    X0,Y0,AP,BP,CP,DP,P00,P10,P20,P30,P40,P50,
     1                   P01,P11,P21,P31,P41,P02,P12,P22,P32,P03,P13,
     2                   P23,P04,P14,P05,ITPV
      EQUIVALENCE        (P5,P50)
C                                  FIRST EXECUTABLE STATEMENT
C                                  PRELIMINARY PROCESSING
      IT0 = ITI
      NTL = NT+NL
      IF (IT0.LE.NTL) GO TO 5
      IL1 = IT0/NTL
      IL2 = IT0-IL1*NTL
      IF (IL1.EQ.IL2) GO TO 30
      GO TO 55
C                                  CALCULATION OF ZII BY INTERPOLATION.
C                                    CHECKS IF THE NECESSARY
C                                    COEFFICIENTS HAVE BEEN CALCULATED.
    5 IF (IT0.EQ.ITPV) GO TO 25
C                                  LOADS COORDINATE AND PARTIAL
C                                    DERIVATIVE VALUES AT THE VERTEXES.
      JIPT = 3*(IT0-1)
      JPD = 0
      DO 15 I=1,3
         JIPT = JIPT+1
         IDP = IPT(JIPT)
         X(I) = XD(IDP)
         Y(I) = YD(IDP)
         Z(I) = ZD(IDP)
         JPDD = 5*(IDP-1)
         DO 10 KPD=1,5
            JPD = JPD+1
            JPDD = JPDD+1
            PD(JPD) = PDD(JPDD)
   10    CONTINUE
   15 CONTINUE
C                                  DETERMINES THE COEFFICIENTS FOR THE
C                                    COORDINATE SYSTEM TRANSFORMATION
C                                    FROM THE X-Y SYSTEM TO THE U-V
C                                    SYSTEM AND VICE VERSA.
      X0 = X(1)
      Y0 = Y(1)
      A = X(2)-X0
      B = X(3)-X0
      C = Y(2)-Y0
      D = Y(3)-Y0
      AD = A*D
      BC = B*C
      DLT = AD-BC
      AP = D/DLT
      BP = -B/DLT
      CP = -C/DLT
      DP = A/DLT
C                                  CONVERTS THE PARTIAL DERIVATIVES AT
C                                    THE VERTEXES OF THE TRIANGLE FOR
C                                    THE U-V COORDINATE SYSTEM.
      AA = A*A
      ACT2 = 2.0*A*C
      CC = C*C
      AB = A*B
      ADBC = AD+BC
      CD = C*D
      BB = B*B
      BDT2 = 2.0*B*D
      DD = D*D
      DO 20 I=1,3
         JPD = 5*I
         ZU(I) = A*PD(JPD-4)+C*PD(JPD-3)
         ZV(I) = B*PD(JPD-4)+D*PD(JPD-3)
         ZUU(I) = AA*PD(JPD-2)+ACT2*PD(JPD-1)+CC*PD(JPD)
         ZUV(I) = AB*PD(JPD-2)+ADBC*PD(JPD-1)+CD*PD(JPD)
         ZVV(I) = BB*PD(JPD-2)+BDT2*PD(JPD-1)+DD*PD(JPD)
   20 CONTINUE
C                                  CALCULATES THE COEFFICIENTS OF THE
C                                    POLYNOMIAL.
      P00 = Z(1)
      P10 = ZU(1)
      P01 = ZV(1)
      P20 = 0.5*ZUU(1)
      P11 = ZUV(1)
      P02 = 0.5*ZVV(1)
      H1 = Z(2)-P00-P10-P20
      H2 = ZU(2)-P10-ZUU(1)
      H3 = ZUU(2)-ZUU(1)
      P30 = 10.0*H1-4.0*H2+0.5*H3
      P40 = -15.0*H1+7.0*H2-H3
      P50 = 6.0*H1-3.0*H2+0.5*H3
      H1 = Z(3)-P00-P01-P02
      H2 = ZV(3)-P01-ZVV(1)
      H3 = ZVV(3)-ZVV(1)
      P03 = 10.0*H1-4.0*H2+0.5*H3
      P04 = -15.0*H1+7.0*H2-H3
      P05 = 6.0*H1-3.0*H2+0.5*H3
      LU = SQRT(AA+CC)
      LV = SQRT(BB+DD)
      THXU = ATAN2(C,A)
      THUV = ATAN2(D,B)-THXU
      CSUV = COS(THUV)
      P41 = 5.0*LV*CSUV/LU*P50
      P14 = 5.0*LU*CSUV/LV*P05
      H1 = ZV(2)-P01-P11-P41
      H2 = ZUV(2)-P11-4.0*P41
      P21 = 3.0*H1-H2
      P31 = -2.0*H1+H2
      H1 = ZU(3)-P10-P11-P14
      H2 = ZUV(3)-P11-4.0*P14
      P12 = 3.0*H1-H2
      P13 = -2.0*H1+H2
      THUS = ATAN2(D-C,B-A)-THXU
      THSV = THUV-THUS
      AA = SIN(THSV)/LU
      BB = -COS(THSV)/LU
      CC = SIN(THUS)/LV
      DD = COS(THUS)/LV
      AC = AA*CC
      AD = AA*DD
      BC = BB*CC
      G1 = AA*AC*(3.0*BC+2.0*AD)
      G2 = CC*AC*(3.0*AD+2.0*BC)
      H1 = -AA*AA*AA*(5.0*AA*BB*P50+(4.0*BC+AD)*P41)-CC*CC*CC*(5.0*CC
     1*DD*P05+(4.0*AD+BC)*P14)
      H2 = 0.5*ZVV(2)-P02-P12
      H3 = 0.5*ZUU(3)-P20-P21
      P22 = (G1*H2+G2*H3-H1)/(G1+G2)
      P32 = H2-P22
      P23 = H3-P22
      ITPV = IT0
C                                  CONVERTS XII AND YII TO U-V SYSTEM.
   25 DX = XII-X0
      DY = YII-Y0
      U = AP*DX+BP*DY
      V = CP*DX+DP*DY
C                                  EVALUATES THE POLYNOMIAL.
      P0 = P00+V*(P01+V*(P02+V*(P03+V*(P04+V*P05))))
      P1 = P10+V*(P11+V*(P12+V*(P13+V*P14)))
      P2 = P20+V*(P21+V*(P22+V*P23))
      P3 = P30+V*(P31+V*P32)
      P4 = P40+V*P41
      ZII = P0+U*(P1+U*(P2+U*(P3+U*(P4+U*P5))))
      IF(TENS.GT.0. .AND. TENS.LE.1.)GOTO 70
      RETURN
C                                  CALCULATION OF ZII BY EXTRAPOLATION
C                                    IN THE RECTANGLE. CHECKS IF THE
C                                    NECESSARY COEFFICIENTS HAVE BEEN
C                                    CALCULATED.
   30 IF (IEX .NE. 0) THEN
         ZII = 9999.
         RETURN
      END IF
      IF (IT0.EQ.ITPV) GO TO 50
C                                  LOADS COORDINATE AND PARTIAL
C                                    DERIVATIVE VALUES AT THE END
C                                    POINTS OF THE BORDER LINE SEGMENT.
      JIPL = 3*(IL1-1)
      JPD = 0
      DO 40 I=1,2
         JIPL = JIPL+1
         IDP = IPL(JIPL)
         X(I) = XD(IDP)
         Y(I) = YD(IDP)
         Z(I) = ZD(IDP)
         JPDD = 5*(IDP-1)
         DO 35 KPD=1,5
            JPD = JPD+1
            JPDD = JPDD+1
            PD(JPD) = PDD(JPDD)
   35    CONTINUE
   40 CONTINUE
C                                  DETERMINES THE COEFFICIENTS FOR THE
C                                    COORDINATE SYSTEM TRANSFORMATION
C                                    FROM THE X-Y SYSTEM TO THE U-V
C                                    SYSTEM AND VICE VERSA.
      X0 = X(1)
      Y0 = Y(1)
      A = Y(2)-Y(1)
      B = X(2)-X(1)
      C = -B
      D = A
      AD = A*D
      BC = B*C
      DLT = AD-BC
      AP = D/DLT
      BP = -B/DLT
      CP = -BP
      DP = AP
C                                  CONVERTS THE PARTIAL DERIVATIVES AT
C                                    THE END POINTS OF THE BORDER LINE
C                                    SEGMENT FOR THE U-V COORDINATE
C                                    SYSTEM.
      AA = A*A
      ACT2 = 2.0*A*C
      CC = C*C
      AB = A*B
      ADBC = AD+BC
      CD = C*D
      BB = B*B
      BDT2 = 2.0*B*D
      DD = D*D
      DO 45 I=1,2
         JPD = 5*I
         ZU(I) = A*PD(JPD-4)+C*PD(JPD-3)
         ZV(I) = B*PD(JPD-4)+D*PD(JPD-3)
         ZUU(I) = AA*PD(JPD-2)+ACT2*PD(JPD-1)+CC*PD(JPD)
         ZUV(I) = AB*PD(JPD-2)+ADBC*PD(JPD-1)+CD*PD(JPD)
         ZVV(I) = BB*PD(JPD-2)+BDT2*PD(JPD-1)+DD*PD(JPD)
   45 CONTINUE
C                                  CALCULATES THE COEFFICIENTS OF THE
C                                    POLYNOMIAL.
      P00 = Z(1)
      P10 = ZU(1)
      P01 = ZV(1)
      P20 = 0.5*ZUU(1)
      P11 = ZUV(1)
      P02 = 0.5*ZVV(1)
      H1 = Z(2)-P00-P01-P02
      H2 = ZV(2)-P01-ZVV(1)
      H3 = ZVV(2)-ZVV(1)
      P03 = 10.0*H1-4.0*H2+0.5*H3
      P04 = -15.0*H1+7.0*H2-H3
      P05 = 6.0*H1-3.0*H2+0.5*H3
      H1 = ZU(2)-P10-P11
      H2 = ZUV(2)-P11
      P12 = 3.0*H1-H2
      P13 = -2.0*H1+H2
      P21 = 0.0
      P23 = -ZUU(2)+ZUU(1)
      P22 = -1.5*P23
      ITPV = IT0
C                                  CONVERTS XII AND YII TO U-V SYSTEM.
   50 DX = XII-X0
      DY = YII-Y0
      U = AP*DX+BP*DY
      V = CP*DX+DP*DY
C                                  EVALUATES THE POLYNOMIAL.
      P0 = P00+V*(P01+V*(P02+V*(P03+V*(P04+V*P05))))
      P1 = P10+V*(P11+V*(P12+V*P13))
      P2 = P20+V*(P21+V*(P22+V*P23))
      ZII = P0+U*(P1+U*P2)
      RETURN
C
C                                  CALCULATION OF ZII BY EXTRAPOLATION
C                                    IN THE TRIANGLE. CHECKS IF THE
C                                    NECESSARY COEFFICIENTS HAVE BEEN
C                                    CALCULATED.
   55 IF (IEX .NE. 0) THEN
         ZII = 9999.
         RETURN
      END IF
C
      IF (IT0.EQ.ITPV) GO TO 65
C                                  LOADS COORDINATE AND PARTIAL
C                                    DERIVATIVE VALUES AT THE VERTEX OF
C                                    THE TRIANGLE.
      JIPL = 3*IL2-2
      IDP = IPL(JIPL)
      X0 = XD(IDP)
      Y0 = YD(IDP)
      Z0 = ZD(IDP)
      JPDD = 5*(IDP-1)
      DO 60 KPD=1,5
         JPDD = JPDD+1
         PD(KPD) = PDD(JPDD)
   60 CONTINUE
C                                  CALCULATES THE COEFFICIENTS OF THE
C                                    POLYNOMIAL.
      P00 = Z0
      P10 = PD(1)
      P01 = PD(2)
      P20 = 0.5*PD(3)
      P11 = PD(4)
      P02 = 0.5*PD(5)
      ITPV = IT0
C                                  CONVERTS XII AND YII TO U-V SYSTEM.
   65 U = XII-X0
      V = YII-Y0
C                                  EVALUATES THE POLYNOMIAL.
      P0 = P00+V*(P01+V*P02)
      P1 = P10+V*P11
      ZII = P0+U*(P1+U*P20)
      RETURN
C
70    CONTINUE
C**************************************************************************
C
C     Add a tension term (interpolation only - not if extrapolating):
C
C**************************************************************************
C
C     Switch to barycentric coordinates
C
      DO 75 II=1,3
         CMAT(1,II)=X(II)
         CMAT(2,II)=Y(II)
         CMAT(3,II)=1.
75    CONTINUE
      RHS(1)=XII
      RHS(2)=YII
      RHS(3)=1.
      CALL DECOMP(3,3,CMAT,COND,IPVT,WKAREA)
      IF (COND.LT.1.E+8) THEN
         CALL SOLVE(3,3,CMAT,RHS,IPVT)
         BVL = RHS(1)*Z(1)+RHS(2)*Z(2)+RHS(3)*Z(3)
         ZII = TENS*BVL + (1.-TENS)*ZII
      ELSE
C        Distinguish this failure from no-extrapolation value (9999.):
         ZII = 9997.
      END IF
      RETURN
      END
C-----------------------------------------------------------------------
C
C   COMPUTER            - DEC11/SINGLE
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY SUBROUTINE IQHSCV
C
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - SINGLE/H36,H48,H60
C
C   REQD.  ROUTINES - IQHSD
C
C-----------------------------------------------------------------------
C
      SUBROUTINE IQHSG  (NDP,XD,YD,NT,IPT,NL,IPL,IWL,IWP,WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NDP,NT,NL,IER,IPT(1),IPL(1),IWL(1),IWP(1)
      DIMENSION         XD(1),YD(1)
      REAL                  WK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            ILF,ILIV,ILT3,ILVS,IL,IP1P1,IP1,IP2,IP3,IPL1,
     1                   IPL2,IPLJ1,IPLJ2,IPMN1,IPMN2,IPT1,IPT2,IPT3,
     2                   IPTI1,IPTI2,IPTI,IP,IREP,IT1T3,IT2T3,ITF(2),
     3                   ITS,ITT3R,ITT3,IT,IXVSPV,IXVS,JL1,JL2,JLT3,JP1,
     4                   JP2,JPC,JPMN,JPMX,JP,JWL1MN,JWL1,JWL,NDP0,
     5                   NDPM1,NL0,NLFC,NLFT2,NLF,NLNT3,NLN,NLSHT3,NLSH,
     6                   NLT3,NREP,NT0,NTF,IQHSD,NTT3P3,NTT3
      REAL               DSQI,DSQMN,EPSLN,SP,TOL,VP,X1,X2,X3,XDMP,Y1,Y2,
     1                   Y3,YDMP
      REAL               DSQF,SPDT,VPDT,U1,U2,U3,V1,V2,V3
      DATA               NREP/100/
C                                  MACHINE PRECISION
      DATA               TOL/1.1921E-07/
C                                  STATEMENT FUNCTIONS
      DSQF(U1,V1,U2,V2) = (U2-U1)**2+(V2-V1)**2
      SPDT(U1,V1,U2,V2,U3,V3) = (U2-U1)*(U3-U1)+(V2-V1)*(V3-V1)
      VPDT(U1,V1,U2,V2,U3,V3) = (V3-V1)*(U2-U1)-(U3-U1)*(V2-V1)
C                                  FIRST EXECUTABLE STATEMENT
      EPSLN = TOL*100.0
C                                  PRELIMINARY PROCESSING
      NDP0 = NDP
      NDPM1 = NDP0-1
C                                  DETERMINES THE CLOSEST PAIR OF DATA
C                                    POINTS AND THEIR MIDPOINT.
      DSQMN = DSQF(XD(1),YD(1),XD(2),YD(2))
      IPMN1 = 1
      IPMN2 = 2
      DO 10 IP1=1,NDPM1
         X1 = XD(IP1)
         Y1 = YD(IP1)
         IP1P1 = IP1+1
         DO 5 IP2=IP1P1,NDP0
            DSQI = DSQF(X1,Y1,XD(IP2),YD(IP2))
            IF (DSQI.EQ.0.0) GO TO 160
            IF (DSQI.GE.DSQMN) GO TO 5
            DSQMN = DSQI
            IPMN1 = IP1
            IPMN2 = IP2
    5    CONTINUE
   10 CONTINUE
      XDMP = (XD(IPMN1)+XD(IPMN2))/2.0
      YDMP = (YD(IPMN1)+YD(IPMN2))/2.0
C                                  SORTS THE OTHER (NDP-2) DATA POINTS
C                                    IN ASCENDING ORDER OF DISTANCE
C                                    FROM THE MIDPOINT AND STORES THE
C                                    SORTED DATA POINT NUMBERS IN THE
C                                    IWP ARRAY.
      JP1 = 2
      DO 15 IP1=1,NDP0
         IF (IP1.EQ.IPMN1.OR.IP1.EQ.IPMN2) GO TO 15
         JP1 = JP1+1
         IWP(JP1) = IP1
         WK(JP1) = DSQF(XDMP,YDMP,XD(IP1),YD(IP1))
   15 CONTINUE
      DO 25 JP1=3,NDPM1
         DSQMN = WK(JP1)
         JPMN = JP1
         DO 20 JP2=JP1,NDP0
            IF (WK(JP2).GE.DSQMN) GO TO 20
            DSQMN = WK(JP2)
            JPMN = JP2
   20    CONTINUE
         ITS = IWP(JP1)
         IWP(JP1) = IWP(JPMN)
         IWP(JPMN) = ITS
         WK(JPMN) = WK(JP1)
   25 CONTINUE
C                                  IF NECESSARY, MODIFIES THE ORDERING
C                                    IN SUCH A WAY THAT THE FIRST THREE
C                                    DATA POINTS ARE NOT COLLINEAR.
      X1 = XD(IPMN1)
      Y1 = YD(IPMN1)
      X2 = XD(IPMN2)
      Y2 = YD(IPMN2)
      DO 30 JP=3,NDP0
         IP = IWP(JP)
         SP = SPDT(XD(IP),YD(IP),X1,Y1,X2,Y2)
         VP = VPDT(XD(IP),YD(IP),X1,Y1,X2,Y2)
         IF (ABS(VP).GT.(ABS(SP)*EPSLN)) GO TO 35
   30 CONTINUE
      GO TO 165
   35 IF (JP.EQ.3) GO TO 45
      JPMX = JP
      DO 40 JPC=4,JPMX
         JP = JPMX+4-JPC
         IWP(JP) = IWP(JP-1)
   40 CONTINUE
      IWP(3) = IP
C                                  FORMS THE FIRST TRIANGLE. STORES
C                                    POINT NUMBERS OF THE VERTEXES OF
C                                    THE TRIANGLE IN THE IPT ARRAY, AND
C                                    STORES POINT NUMBERS OF THE
C                                    BORDER LINE SEGMENTS AND THE
C                                    TRIANGLE NUMBER IN THE IPL ARRAY.
   45 IP1 = IPMN1
      IP2 = IPMN2
      IP3 = IWP(3)
      IF (VPDT(XD(IP1),YD(IP1),XD(IP2),YD(IP2),XD(IP3),YD(IP3)).GE.0.0)
     1GO TO 50
      IP1 = IPMN2
      IP2 = IPMN1
   50 NT0 = 1
      NTT3 = 3
      IPT(1) = IP1
      IPT(2) = IP2
      IPT(3) = IP3
      NL0 = 3
      NLT3 = 9
      IPL(1) = IP1
      IPL(2) = IP2
      IPL(3) = 1
      IPL(4) = IP2
      IPL(5) = IP3
      IPL(6) = 1
      IPL(7) = IP3
      IPL(8) = IP1
      IPL(9) = 1
C                                  ADDS THE REMAINING (NDP-3) DATA
C                                    POINTS, ONE BY ONE.
      DO 150 JP1=4,NDP0
         IP1 = IWP(JP1)
         X1 = XD(IP1)
         Y1 = YD(IP1)
C                                  DETERMINES THE FIRST INVISIBLE AND
C                                    VISIBLE BORDER LINE SEGMENTS,
C                                    ILIV AND ILVS.
         DO 65 IL=1,NL0
            IP2 = IPL(3*IL-2)
            IP3 = IPL(3*IL-1)
            X2 = XD(IP2)
            Y2 = YD(IP2)
            X3 = XD(IP3)
            Y3 = YD(IP3)
            SP = SPDT(X1,Y1,X2,Y2,X3,Y3)
            VP = VPDT(X1,Y1,X2,Y2,X3,Y3)
            IF (IL.NE.1) GO TO 55
            IXVS = 0
            IF (VP.LE.(ABS(SP)*(-EPSLN))) IXVS = 1
            ILIV = 1
            ILVS = 1
            GO TO 65
   55       IXVSPV = IXVS
            IF (VP.GT.(ABS(SP)*(-EPSLN))) GO TO 60
            IXVS = 1
            IF (IXVSPV.EQ.1) GO TO 65
            ILVS = IL
            IF (ILIV.NE.1) GO TO 70
            GO TO 65
   60       IXVS = 0
            IF (IXVSPV.EQ.0) GO TO 65
            ILIV = IL
            IF (ILVS.NE.1) GO TO 70
   65    CONTINUE
         IF (ILIV.EQ.1.AND.ILVS.EQ.1) ILVS = NL0
   70    IF (ILVS.LT.ILIV) ILVS = ILVS+NL0
C                                  SHIFTS (ROTATES) THE IPL ARRAY TO
C                                    HAVE THE INVISIBLE BORDER LINE
C                                    SEGMENTS CONTAINED IN THE FIRST
C                                    PART OF THE IPL ARRAY.
         IF (ILIV.EQ.1) GO TO 85
         NLSH = ILIV-1
         NLSHT3 = NLSH*3
         DO 75 JL1=1,NLSHT3
            JL2 = JL1+NLT3
            IPL(JL2) = IPL(JL1)
   75    CONTINUE
         DO 80 JL1=1,NLT3
            JL2 = JL1+NLSHT3
            IPL(JL1) = IPL(JL2)
   80    CONTINUE
         ILVS = ILVS-NLSH
C                                  ADDS TRIANGLES TO THE IPT ARRAY,
C                                    UPDATES BORDER LINE SEGMENTS IN
C                                    THE IPL ARRAY, AND SETS FLAGS FOR
C                                    THE BORDER LINE SEGMENTS TO BE
C                                    REEXAMINED IN THE IWL ARRAY.
   85    JWL = 0
         DO 105 IL=ILVS,NL0
            ILT3 = IL*3
            IPL1 = IPL(ILT3-2)
            IPL2 = IPL(ILT3-1)
            IT = IPL(ILT3)
C                                  ADDS A TRIANGLE TO THE IPT
C                                    ARRAY.
            NT0 = NT0+1
            NTT3 = NTT3+3
            IPT(NTT3-2) = IPL2
            IPT(NTT3-1) = IPL1
            IPT(NTT3) = IP1
C                                  UPDATES BORDER LINE SEGMENTS IN
C                                    THE IPL ARRAY.
            IF (IL.NE.ILVS) GO TO 90
            IPL(ILT3-1) = IP1
            IPL(ILT3) = NT0
   90       IF (IL.NE.NL0) GO TO 95
            NLN = ILVS+1
            NLNT3 = NLN*3
            IPL(NLNT3-2) = IP1
            IPL(NLNT3-1) = IPL(1)
            IPL(NLNT3) = NT0
C                                  DETERMINES THE VERTEX THAT DOES
C                                    NOT LIE ON THE BORDER LINE
C                                    SEGMENTS.
   95       ITT3 = IT*3
            IPTI = IPT(ITT3-2)
            IF (IPTI.NE.IPL1.AND.IPTI.NE.IPL2) GO TO 100
            IPTI = IPT(ITT3-1)
            IF (IPTI.NE.IPL1.AND.IPTI.NE.IPL2) GO TO 100
            IPTI = IPT(ITT3)
C                                  CHECKS IF THE EXCHANGE IS
C                                    NECESSARY.
C
  100       IF (IQHSD(XD,YD,IP1,IPTI,IPL1,IPL2).EQ.0) GO TO 105
C
C                                  MODIFIES THE IPT ARRAY WHEN
C                                    NECESSARY.
            IPT(ITT3-2) = IPTI
            IPT(ITT3-1) = IPL1
            IPT(ITT3) = IP1
            IPT(NTT3-1) = IPTI
            IF (IL.EQ.ILVS) IPL(ILT3) = IT
            IF (IL.EQ.NL0.AND.IPL(3).EQ.IT) IPL(3) = NT0
C
C                                  SETS FLAGS IN THE IWL ARRAY.
            JWL = JWL+4
            IWL(JWL-3) = IPL1
            IWL(JWL-2) = IPTI
            IWL(JWL-1) = IPTI
            IWL(JWL) = IPL2
  105    CONTINUE
         NL0 = NLN
         NLT3 = NLNT3
         NLF = JWL/2
         IF (NLF.EQ.0) GO TO 150
C                                  IMPROVES TRIANGULATION.
         NTT3P3 = NTT3+3
         DO 145 IREP=1,NREP
            DO 135 ILF=1,NLF
               IPL1 = IWL(2*ILF-1)
               IPL2 = IWL(2*ILF)
C                                  LOCATES IN THE IPT ARRAY TWO
C                                    TRIANGLES ON BOTH SIDES OF THE
C                                    FLAGGED LINE SEGMENT.
               NTF = 0
               DO 110 ITT3R=3,NTT3,3
                  ITT3 = NTT3P3-ITT3R
                  IPT1 = IPT(ITT3-2)
                  IPT2 = IPT(ITT3-1)
                  IPT3 = IPT(ITT3)
                  IF (IPL1.NE.IPT1.AND.IPL1.NE.IPT2.AND.IPL1.NE.IPT3) GO
     1             TO 110
                  IF (IPL2.NE.IPT1.AND.IPL2.NE.IPT2.AND.IPL2.NE.IPT3) GO
     1             TO 110
                  NTF = NTF+1
                  ITF(NTF) = ITT3/3
                  IF (NTF.EQ.2) GO TO 115
  110          CONTINUE
               IF (NTF.LT.2) GO TO 135
C                                  DETERMINES THE VERTEXES OF THE
C                                    TRIANGLES THAT DO NOT LIE ON
C                                    THE LINE SEGMENT.
  115          IT1T3 = ITF(1)*3
               IPTI1 = IPT(IT1T3-2)
               IF (IPTI1.NE.IPL1.AND.IPTI1.NE.IPL2) GO TO 120
               IPTI1 = IPT(IT1T3-1)
               IF (IPTI1.NE.IPL1.AND.IPTI1.NE.IPL2) GO TO 120
               IPTI1 = IPT(IT1T3)
  120          IT2T3 = ITF(2)*3
               IPTI2 = IPT(IT2T3-2)
               IF (IPTI2.NE.IPL1.AND.IPTI2.NE.IPL2) GO TO 125
               IPTI2 = IPT(IT2T3-1)
               IF (IPTI2.NE.IPL1.AND.IPTI2.NE.IPL2) GO TO 125
               IPTI2 = IPT(IT2T3)
C                                  CHECKS IF THE EXCHANGE IS
C                                    NECESSARY.
C
  125          IF (IQHSD(XD,YD,IPTI1,IPTI2,IPL1,IPL2).EQ.0) GO TO 135
C                                  MODIFIES THE IPT ARRAY WHEN
C                                    NECESSARY.
               IPT(IT1T3-2) = IPTI1
               IPT(IT1T3-1) = IPTI2
               IPT(IT1T3) = IPL1
               IPT(IT2T3-2) = IPTI2
               IPT(IT2T3-1) = IPTI1
               IPT(IT2T3) = IPL2
C                                  SETS NEW FLAGS.
               JWL = JWL+8
               IWL(JWL-7) = IPL1
               IWL(JWL-6) = IPTI1
               IWL(JWL-5) = IPTI1
               IWL(JWL-4) = IPL2
               IWL(JWL-3) = IPL2
               IWL(JWL-2) = IPTI2
               IWL(JWL-1) = IPTI2
               IWL(JWL) = IPL1
               DO 130 JLT3=3,NLT3,3
                  IPLJ1 = IPL(JLT3-2)
                  IPLJ2 = IPL(JLT3-1)
                  IF ((IPLJ1.EQ.IPL1.AND.IPLJ2.EQ.IPTI2).OR.(IPLJ2.EQ.IP
     1            L1.AND.IPLJ1.EQ.IPTI2)) IPL(JLT3)
     2            = ITF(1)
                  IF ((IPLJ1.EQ.IPL2.AND.IPLJ2.EQ.IPTI1).OR.(IPLJ2.EQ.IP
     1            L2.AND.IPLJ1.EQ.IPTI1)) IPL(JLT3)
     2            = ITF(2)
  130          CONTINUE
  135       CONTINUE
            NLFC = NLF
            NLF = JWL/2
            IF (NLF.EQ.NLFC) GO TO 150
C                                  RESETS THE IWL ARRAY FOR THE
C                                    NEXT ROUND.
            JWL1MN = 2*NLFC+1
            NLFT2 = NLF*2
            DO 140 JWL1=JWL1MN,NLFT2
               JWL = JWL1+1-JWL1MN
               IWL(JWL) = IWL(JWL1)
  140       CONTINUE
            NLF = JWL/2
  145    CONTINUE
  150 CONTINUE
C                                  REARRANGES THE IPT ARRAY SO THAT THE
C                                    VERTEXES OF EACH TRIANGLE ARE
C                                    LISTED COUNTER-CLOCKWISE.
      DO 155 ITT3=3,NTT3,3
         IP1 = IPT(ITT3-2)
         IP2 = IPT(ITT3-1)
         IP3 = IPT(ITT3)
         IF (VPDT(XD(IP1),YD(IP1),XD(IP2),YD(IP2),XD(IP3),YD(IP3)).GE.0.
     1   0) GO TO 155
         IPT(ITT3-2) = IP2
         IPT(ITT3-1) = IP1
  155 CONTINUE
      NT = NT0
      NL = NL0
      RETURN
C                                  ERROR EXIT
  160 IER = 131
      RETURN
  165 IER = 130
      RETURN
      END
C-----------------------------------------------------------------------
C
C   COMPUTER            - DEC11/SINGLE
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY SUBROUTINE IQHSCV
C
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - SINGLE/H36,H48,H60
C
C   REQD. ROUTINES - NONE
C
C-----------------------------------------------------------------------
C
      SUBROUTINE IQHSH  (XD,YD,NT,IPT,NL,IPL,NPI,XI,YI,NGP,IGP)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NT,NL,NXI,NYI,IPT(1),IPL(1),NGP(1),IGP(1)
      DIMENSION            XD(1),YD(1)
      DIMENSION            XI(1),YI(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IL0T3,IL0,ILP1T3,ILP1,INSD,IP1,IP2,IP3,IT0T3,
     1                   IT0,IZI,JIGP0,JIGP1I,JIGP1,
     2                   JNGP0,JNGP1,L,NGP0,NGP1,NL0,NT0,NXINYI,
     3                   NPI0
      REAL               X1,X2,X3,XII,XIMN,XIMX,XMN,XMX,Y1,Y2,Y3,YII,
     1                   YIMN,YIMX,YMN,YMX
      REAL               SPDT,VPDT,U1,U2,U3,V1,V2,V3
C                                  STATEMENT FUNCTIONS
      SPDT(U1,V1,U2,V2,U3,V3) = (U1-U2)*(U3-U2)+(V1-V2)*(V3-V2)
      VPDT(U1,V1,U2,V2,U3,V3) = (U1-U3)*(V2-V3)-(V1-V3)*(U2-U3)
C                                  FIRST EXECUTABLE STATEMENT
C
      NT0 = NT
      NL0 = NL
      NPI0 = NPI
      NXINYI = NPI
      XIMN = AMIN1(XI(1),XI(NPI0))
      XIMX = AMAX1(XI(1),XI(NPI0))
      YIMN = AMIN1(YI(1),YI(NPI0))
      YIMX = AMAX1(YI(1),YI(NPI0))
C                                  DETERMINES GRID POINTS INSIDE THE
C                                    DATA AREA.
      JNGP0 = 0
      JNGP1 = 2*(NT0+2*NL0)+1        
      JIGP0 = 0
      JIGP1 = NXINYI+1
C
C     For each triangle determine the set of data points inside
C
      DO 80 IT0=1,NT0
         NGP0 = 0
         NGP1 = 0
         IT0T3 = IT0*3
         IP1 = IPT(IT0T3-2)
         IP2 = IPT(IT0T3-1)
         IP3 = IPT(IT0T3)
C
C        The vertices of the triangle
C
         X1 = XD(IP1)
         Y1 = YD(IP1)
         X2 = XD(IP2)
         Y2 = YD(IP2)
         X3 = XD(IP3)
         Y3 = YD(IP3)
C
C        Min and max values in both variables of the triangle
C
         XMN = AMIN1(X1,X2,X3)
         XMX = AMAX1(X1,X2,X3)
         YMN = AMIN1(Y1,Y2,Y3)
         YMX = AMAX1(Y1,Y2,Y3)
   15    DO 70 IPI=1,NPI0
            YII = YI(IPI)
            XII = XI(IPI)
C
C           Gross check to see if there is any possibility of the point
C           lying inside the triangle
C
            IF (YII.LT.YMN.OR.YII.GT.YMX) GO TO 70
            IF (XII.LT.XMN.OR.XII.GT.XMX) GO TO 70
C
C           Check if the point lies in the triangle
C
            L = 0
            IF (VPDT(X1,Y1,X2,Y2,XII,YII)) 70,20,25
   20       L = 1
   25       IF (VPDT(X2,Y2,X3,Y3,XII,YII)) 70,30,35
   30       L = 1
   35       IF (VPDT(X3,Y3,X1,Y1,XII,YII)) 70,40,45
   40       L = 1
C
C           Point lies in the triangle
C
   45       CONTINUE
C
C           If L=1 point lies on a segment of the triangle
C
            IF (L.EQ.1) GO TO 50
            NGP0 = NGP0+1
            JIGP0 = JIGP0+1
            IGP(JIGP0) = IPI
            GO TO 70
   50       IF (JIGP1.GT.NXINYI) GO TO 60
            DO 55 JIGP1I=JIGP1,NXINYI
               IF (IPI.EQ.IGP(JIGP1I)) GO TO 70
   55       CONTINUE
   60       NGP1 = NGP1+1
            JIGP1 = JIGP1-1
            IGP(JIGP1) = IPI
   70    CONTINUE
C
C        End of YI checks
C
   75    JNGP0 = JNGP0+1
         NGP(JNGP0) = NGP0
         JNGP1 = JNGP1-1
         NGP(JNGP1) = NGP1
   80 CONTINUE
C                                  DETERMINES GRID POINTS OUTSIDE THE
C                                    DATA AREA. - IN SEMI-INFINITE
C                                    RECTANGULAR AREA.
      DO 225 IL0=1,NL0
         NGP0 = 0
         NGP1 = 0
         IL0T3 = IL0*3
         IP1 = IPL(IL0T3-2)
         IP2 = IPL(IL0T3-1)
         X1 = XD(IP1)
         Y1 = YD(IP1)
         X2 = XD(IP2)
         Y2 = YD(IP2)
         XMN = XIMN
         XMX = XIMX
         YMN = YIMN
         YMX = YIMX
         IF (Y2.GE.Y1) XMN = AMIN1(X1,X2)
         IF (Y2.LE.Y1) XMX = AMAX1(X1,X2)
         IF (X2.LE.X1) YMN = AMIN1(Y1,Y2)
         IF (X2.GE.X1) YMX = AMAX1(Y1,Y2)
C
   95    DO 150 IPI=1,NPI0
            YII = YI(IPI)
            IF (YII.LT.YMN.OR.YII.GT.YMX) GO TO 150
            XII = XI(IPI)
            L = 0
            IF (VPDT(X1,Y1,X2,Y2,XII,YII)) 105,100,150
  100       L = 1
  105       IF (SPDT(X2,Y2,X1,Y1,XII,YII)) 150,110,115
  110       L = 1
  115       IF (SPDT(X1,Y1,X2,Y2,XII,YII)) 150,120,125
  120       L = 1
  125       IZI = IPI
            IF (L.EQ.1) GO TO 130
            NGP0 = NGP0+1
            JIGP0 = JIGP0+1
            IGP(JIGP0) = IZI
            GO TO 150
  130       IF (JIGP1.GT.NXINYI) GO TO 140
            DO 135 JIGP1I=JIGP1,NXINYI
               IF (IZI.EQ.IGP(JIGP1I)) GO TO 150
  135       CONTINUE
  140       NGP1 = NGP1+1
            JIGP1 = JIGP1-1
            IGP(JIGP1) = IZI
  150    CONTINUE
  155    JNGP0 = JNGP0+1
         NGP(JNGP0) = NGP0
         JNGP1 = JNGP1-1
         NGP(JNGP1) = NGP1
C                                  - IN SEMI-INFINITE TRIANGULAR AREA.
         NGP0 = 0
         NGP1 = 0
         ILP1 = MOD(IL0,NL0)+1
         ILP1T3 = ILP1*3
         IP3 = IPL(ILP1T3-1)
         X3 = XD(IP3)
         Y3 = YD(IP3)
         XMN = XIMN
         XMX = XIMX
         YMN = YIMN
         YMX = YIMX
         IF (Y3.GE.Y2.AND.Y2.GE.Y1) XMN = X2
         IF (Y3.LE.Y2.AND.Y2.LE.Y1) XMX = X2
         IF (X3.LE.X2.AND.X2.LE.X1) YMN = Y2
         IF (X3.GE.X2.AND.X2.GE.X1) YMX = Y2
  170    DO 215 IPI=1,NPI0
            YII = YI(IPI)
            IF (YII.LT.YMN.OR.YII.GT.YMX) GO TO 215
            XII = XI(IPI)
            L = 0
            IF (SPDT(X1,Y1,X2,Y2,XII,YII)) 180,175,215
  175       L = 1
  180       IF (SPDT(X3,Y3,X2,Y2,XII,YII)) 190,185,215
  185       L = 1
  190       IZI = IPI
            IF (L.EQ.1) GO TO 195
            NGP0 = NGP0+1
            JIGP0 = JIGP0+1
            IGP(JIGP0) = IZI
            GO TO 215
  195       IF (JIGP1.GT.NXINYI) GO TO 205
            DO 200 JIGP1I=JIGP1,NXINYI
               IF (IZI.EQ.IGP(JIGP1I)) GO TO 215
  200       CONTINUE
  205       NGP1 = NGP1+1
            JIGP1 = JIGP1-1
            IGP(JIGP1) = IZI
  215    CONTINUE
  220    JNGP0 = JNGP0+1
         NGP(JNGP0) = NGP0
         JNGP1 = JNGP1-1
         NGP(JNGP1) = NGP1
  225 CONTINUE
      RETURN
      END
