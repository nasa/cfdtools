C+-------------------------------------------------------------------
C
      PROGRAM FOURDIGIT
C
C     (formerly PROGRAM ANALIN)
C
C  PURPOSE:  To calculate the ordinates and surface slope for both
C            symmetrical and cambered airfoils in the NACA 4-digit, 
C            5-digit, and 16-series airfoil families.
C
C  NAMELIST "INPUTS":
C
C   VAR     DIM   TYPE    DESCRIPTION
C  NAME      -    C*80    Title desired on printed and plotted output 
C  PROFILE   -    C*10    NACA airfoil family (4-digit)
C  CAMBER    -    C*10    Camber line (2-digit or 3-digit; upper case;
C                         left justified.  E.g. '4-DIGIT   '
C  TOC       -     R      Thickness-chord ratio of airfoil.
C  RLE       -     R      Leading-edge radius. Not used with 4-digit but
C                         may be used with 4-digit modified.
C  DX        -     R      Basic chordwise increment in x/c for computing 
C                         ordinates.  Usually set to 0.01.
C  CHD       -     R      Model chord used for listing ordinates in
C                         dimensional units.
C  XM        -     R      Nondimensional chordwise location of maximum thick-
C                         ness.  Used for 4-digit modified airfoils only.
C  D1        -     R      Trailing-edge slope of 4-digit modified airfoils.
C  CMB       -     R      Maximum camber ordinate-to-chord ratio (e.g., .04).
C  CM        -     R      Location of maximum camber position (e.g., 0.4).
C                         Note: CMB=K1 and CM=r for 5-digit airfoils, where
C                               K1 and r are defined in NASA TM X-3284.
C  K2OK1     -     R      Defined in NASA TM X-3284.
C
C  ENVIRONMENT:  DEC VAX/VMS FORTRAN (this version)
C
C  REFERENCE:  NASA Technical Memorandum TM X-3284 (November, 1975)
C              Charles Ladson and Cuyler Brooks
C              Langley Research Center
C
C------------------------------------------------------------------------

      COMMON XU,XL,YU,YL,YUPR(200),YLPR(200),I
      COMMON /MAIN/ YSTART(3),CHD,KON
      DIMENSION XU(200), XL(200), YU(200), YL(200)
      DIMENSION COEFFS(4)
      DIMENSION XAU(200),YAU(200),XAL(200),YAL(200)
      DIMENSION CLI(10), A(10), TANTH0(10), YCMB(10), TANTH(10),
     +          YCP2(10), IF6XA(10) 
      COMPLEX ROOTS(3)
      EQUIVALENCE (CLI(1),CMB)
      REAL K2OK1
      CHARACTER PROFILE*10, CAMBER*10, NAME*80 
      NAMELIST /INPUTS/ NAME, PROFILE, CAMBER, TOC, RLE, DX, CHD, XM,
     >                  D1, CMB, CM, K2OK1, CLI, A
      PARAMETER ( LUNSAV=13, LUNRD = 5, LUNOUT=6 )

      PI=3.141592654
      YSTART(1)=1.0
      YSTART(2)=4.0
      YSTART(3)=7.0
      KON=0
      D1=1.0
      XM=0.5

   20 READ (LUNRD,INPUTS,END=999) 
      KON=KON+1
      ICKY=0
      FRAC=1.0
CCCCC D1=1.0
      WRITE(LUNOUT, 660) PROFILE,CAMBER
CCCCC XM=0.5
      IF (CAMBER.EQ.'6-SERIES  ' .OR. CAMBER.EQ.'6A-SERIES ') CMB=CLI(1)
      WRITE(LUNOUT,690) NAME
      IF (ICKY.LE.1) ICKY=1
      E=1.E-6
      DO 50 I=1,10
      IF6XA(I)=0
   50 CONTINUE
      IF (PROFILE.EQ.'4-DIGIT   ') WRITE(LUNOUT,700) TOC,RLE,DX,CHD
      IF (PROFILE.EQ.'4-DIGITMOD') WRITE(LUNOUT,710) TOC,RLE,DX,CHD,
     >                                               XM,D1
      IF (CAMBER .EQ.'2-DIGIT   ') WRITE(LUNOUT,720) CMB,CM
      IF (CAMBER .EQ.'3-DIGIT   ') WRITE(LUNOUT,730) CMB,CM
      IF (CAMBER .EQ.'3-DIGITREF') WRITE(LUNOUT,740) CMB,CM,K2OK1
      IF (CAMBER .EQ.'6-SERIES  ' .OR. CAMBER.EQ.'6A-SERIES')
     >   WRITE(LUNOUT,750)
      IF (CAMBER .EQ.'6-SERIES  ' .OR. CAMBER.EQ.'6A-SERIES')
     >   WRITE(LUNOUT,760) (CLI(I),A(I),I=1,ICKY)
      IF (PROFILE.EQ.'4-DIGIT   ') GO TO 110

C     COMPUTED CONSTANTS

      A0=SQRT(2.0*RLE)*0.2/TOC
      D0=0.002
      IF (D1.GT.0.0) GO TO 100
      D1=.1*(2.24-5.42*XM+12.3*XM**2)/(1.-0.878*XM)
      COEFFS(1)=4.5*(1.-XM)**2
      COEFFS(2)=-1.274*(1.-XM)
      COEFFS(3)=.5*(1.-XM)**2+.086436
      COEFFS(4)=-.098*(1.-XM)
      J=1
      D=D1
      WRITE(LUNOUT,770) J,IERR,D1
      DO 80 I=1,3
      IF (AIMAG(ROOTS(I)).GT.1.E-20) GO TO 80
      GO TO (60,70,70), J
   60 TEST=ABS(REAL(ROOTS(I))-D1)
      D=REAL(ROOTS(I))
   70 IF (ABS(REAL(ROOTS(I))-D1).LT.TEST) D=REAL(ROOTS(I))
      IF (I.EQ.3) GO TO 90
      IF (ABS(REAL(ROOTS(I))-D1).LT.TEST) TEST=ABS(REAL(ROOTS(I))-D1)
   80 CONTINUE  
      J=J+1
      IF (J.EQ.4) WRITE(LUNOUT,780)
      GO TO 100
   90 D1=D
      WRITE(LUNOUT,770) J,IERR,D1
  100 CONTINUE
      D3=(3.*D1-0.588/(1.-XM))/(3.*(1.-XM)**2)
      D2=-1.5*(1.-XM)*D3-.5*D1/(1.-XM)
      A3=0.1/XM**3+(2.*D1*(1.-XM)-0.588)/(2.*XM*(1.-XM)**2)-3.*A0/(8.*XM
     +**2.5)
      A2=-0.10/XM**2+.5*A0/XM**1.5-2.*XM*A3
      A1=-.5*A0/XM**.5-2.*XM*A2-3.*XM**2*A3
C     RC IS RADIUS OF CURVATURE AT X=XM
      RC=((1.-XM)**2/(2.*D1*(1.-XM)-0.588))*.2/TOC
      WRITE(LUNOUT,790) A0,A1,A2,A3,D0,D1,D2,D3,RC
C
C     PROFILE, X LE XM
C
  110 CONTINUE
      IF (CMB.LE.1.E-6) WRITE(LUNOUT,800)
      IF (CMB.GT.1.E-6) WRITE(LUNOUT,810)
      X=0.0
      Y=0.0
      XC=0.0
      YC=0.0
      XU(1)=0.0
      YU(1)=0.0
      XL(1)=0.0
      YL(1)=0.0
      XUC=0.0
      YUC=0.0
      XLC=0.0
      YLC=0.0
      XAU(1)=0.0
      YAU(1)=0.0
      XAL(1)=0.0
      YAL(1)=0.0
      K=2
      IF (CAMBER.EQ.'2-DIGIT   ') GO TO 120
      IF (CAMBER.EQ.'3-DIGIT   ') GO TO 130
      IF (CAMBER.EQ.'3-DIGITREF') GO TO 140
      IF (CAMBER.EQ.'6-SERIES  ') GO TO 150
      IF (CAMBER.EQ.'6A-SERIES ') GO TO 150
      WRITE(LUNOUT,820) CAMBER
      GO TO 20
  120 TANTH0(1)=2.*CMB/CM
      IF (CMB.LT.E) TANTH0(1)=E
      YP=1.E10
      YPP=1.E10
      YUP=-1./TANTH0(1)
      YLP=-1./TANTH0(1)
      GO TO 230
  130 TANTH0(1)=CMB*CM**2*(3.0-CM)/6.0
      IF (CMB.LT.E) TANTH0(1)=E
      YP=1.E10
      YPP=1.E10
      YUP=-1./TANTH0(1)
      YLP=-1./TANTH0(1)
      GO TO 230
  140 TANTH0(1)=CMB*(3.*CM**2-K2OK1*(1-CM)**3-CM**3)/6.0
      IF (CMB.LT.E) TANTH0(1)=E
      YP=1.E10
      YPP=1.E10
      YUP=-1./TANTH0(1)
      YLP=-1./TANTH0(1)
      GO TO 230
  150 L=0
      CLIS=CLI(1)
      AS=A(1)
  160 L=L+1
      A(1)=A(L)
      CLI(1)=CLI(L)
      K=2
      U=0.005
      V=-(A(1)-U)/ABS(A(1)-U)
      OMXL=(1.-U)*ALOG(1.-U)
      AMXL=(A(1)-U)*ALOG(ABS(A(1)-U))
      OMXL1=-ALOG(1.-U)-1.
      AMXL1=-ALOG(ABS(A(1)-U))+V
      OMXL2=1./(1.-U)
      AMXL2=-V/ABS(A(1)-U)
      IF (A(1).LT.E.OR.ABS(1.-A(1)).LT.E) GO TO 170
      G=-(A(1)**2*(.5*ALOG(A(1))-0.25)+0.25)/(1.-A(1))
      Q=1.0
      H=(0.5*(1.-A(1))**2*ALOG(1.-A(1))-0.25*(1.-A(1))**2)/(1.-A(1))+G
      Z=.5*(A(1)-U)*AMXL-.5*(1.-U)*OMXL-.25*(A(1)-U)**2+.25*(1.-U)**2
      Z1=.5*((A(1)-U)*AMXL1-AMXL-(1.-U)*OMXL1+OMXL+(A(1)-U)-(1.-U))
      Z2=.5*(A(1)-U)*AMXL2-AMXL1-.5*(1.-U)*OMXL2+OMXL1
  170 CONTINUE
      IF (A(1).LT.E) GO TO 180
      IF (ABS(A(1)-1.).LT.E) GO TO 190
  180 H=-.5
      Q=1.0
      Z1=U*ALOG(U)-.5*U-.5*(1.-U)*OMXL1+.5*OMXL-.5
      GO TO 200
  190 H=0.0
      Q=H
      Z1=-OMXL1
      GO TO 200
  200 TANTH0(L)=CLI(1)*(Z1/(1.-Q*A(1))-1.-ALOG(U)-H)/PI/(A(1)+1.)/2.0
      IF (ICKY.GT.1.AND.L.LT.ICKY) GO TO 160
      IF (ICKY.EQ.1) GO TO 220
      DO 210 J=2,ICKY
  210 TANTH0(1)=TANTH0(1)+TANTH0(J)
  220 CONTINUE
      IF (CMB.LT.E) TANTH0(1)=E
      YP=1.E10
      YPP=1.E10
      YUP=-1./TANTH0(1)
      YLP=-1./TANTH0(1)
  230 CONTINUE
      I=1
      IF (CMB.GT.1.E-6) WRITE(LUNOUT,850) X,XU(I),YU(I),XUC,YUC,YUP,
     >                                    XL(I),YL(I),XLC,YLC,YLP
      IF (CMB.LE.1.E-6) WRITE(LUNOUT,830) X,Y,YP,YPP,XC,YC
      X=.00025
  240 CONTINUE
      IF (PROFILE.EQ.'4-DIGIT   ') GO TO 250
      IF (PROFILE.EQ.'4-DIGITMOD') GO TO 260
      WRITE(LUNOUT,860) PROFILE
      GO TO 20
  250 Y=0.29690*SQRT(X)-0.12600*X-0.35160*X**2+0.28430*X**3-0.1015*X**4
      YP=.5*.2969/SQRT(X)-.126-2.*.3516*X+3.*.2843*X*X-4.*.1015*X**3
      YPP=-.5*.5*.2969/SQRT(X**3)-2.*.3516+2.*3.*.2843*X-3.*4.*.1015*X*X
      GO TO 270
  260 Y=A0*X**.5+A1*X+A2*X**2+A3*X**3
      YP=.5*A0/X**.5+A1+2.*A2*X+3.*A3*X**2
      YPP=-.25*A0/X**1.5+2.*A2+6.*A3*X
  270 CONTINUE
      Y=Y*TOC/.2
      YP=YP*TOC/.2
      YPP=YPP*TOC/.2
      XC=X*CHD
      YC=Y*CHD
      IF(CMB.LT.E) CM=0.5
      IF (CAMBER.EQ.'2-DIGIT   ') GO TO 280
      IF (CAMBER.EQ.'3-DIGIT   ') GO TO 290
      IF (CAMBER.EQ.'3-DIGITREF') GO TO 300
      IF (CAMBER.EQ.'6-SERIES  ') GO TO 310
      IF (CAMBER.EQ.'6A-SERIES ') GO TO 310
      WRITE(LUNOUT,820) CAMBER
      GO TO 20
  280 YCMB(1)=CMB*(2.0*CM*X-X*X)/CM**2
      TANTH(1)=2.*CMB*(1.-X/CM)/CM
      IF (X.GT.CM) YCMB(1)=CMB*(1.-2.*CM+2.*CM*X-X*X)/(1.-CM)**2
      IF (X.GT.CM) TANTH(1)=(2.*CM-2.*X)*CMB/(1.-CM)**2
      F=SQRT(1.+TANTH(1)**2)
      THP=-2.*CMB/CM**2/F**2
      IF (X.GT.CM) THP=-2.*CMB/(1.-CM)**2/F**2
      GO TO 480
  290 YCMB(1)=CMB*(X**3-3.*CM*X**2+CM**2*(3.-CM)*X)/6.
      TANTH(1)=CMB*(3.*X**2-6.*CM*X+CM**2*(3.-CM))/6.
      IF (X.GT.CM) YCMB(1)=CMB*CM**3*(1.-X)/6.
      IF (X.GT.CM) TANTH(1)=-CMB*CM**3/6.
      F=SQRT(1.+TANTH(1)**2)
      THP=CMB*(X-CM)/F**2
      IF (X.GT.CM) THP=0.0
      GO TO 480
  300 YCMB(1)=CMB*((X-CM)**3-K2OK1*(1-CM)**3*X-CM**3*X+CM**3)/6
      TANTH(1)=CMB*(3.*(X-CM)**2-K2OK1*(1-CM)**3-CM**3)/6.
      IF (X.GT.CM) YCMB(1)=CMB*(K2OK1*(X-CM)**3-K2OK1*(1-CM)**3*X-CM**3*
     +X+CM**3)/6.
      IF (X.GT.CM) 
     +     TANTH(1)=CMB*(3*K2OK1*(X-CM)**2-K2OK1*(1-CM)**3-CM**3)/6.
      F=SQRT(1.+TANTH(1)**2)
      THP=CMB*(X-CM)/F**2
      IF (X.GT.CM) THP=K2OK1*CMB*(X-CM)/F**2
      GO TO 480
  310 L=0
      A(1)=AS
      CLI(1)=CLIS
  320 L=L+1
      CLI(1)=CLI(L)
      XC=X*CHD
      YC=Y*CHD
      XLL=X*ALOG(X)
      Q=1.0
      IF (ABS(1.-A(1)).LT.E.AND.ABS(1.-X).LT.E) GO TO 370
      IF (A(1).LT.E.AND.(1.-X).LT.E) GO TO 380
      IF (ABS(A(1)-X).LT.E) GO TO 330
      IF (ABS(1.-X).LT.E) GO TO 350
      IF (ABS(A(1)-1.).LT.E) GO TO 360
      V=-(A(1)-X)/ABS(A(1)-X)
      OMXL=(1.-X)*ALOG(1.-X)
      AMXL=(A(1)-X)*ALOG(ABS(A(1)-X))
      OMXL1=-ALOG(1.-X)-1.
      AMXL1=-ALOG(ABS(A(1)-X))+V
      OMXL2=1./(1.-X)
      AMXL2=1./(A(1)-X)
      Z=.5*(A(1)-X)*AMXL-.5*(1.-X)*OMXL-.25*(A(1)-X)**2+.25*(1.-X)**2
      Z1=.5*((A(1)-X)*AMXL1-AMXL-(1.-X)*OMXL1+OMXL+(A(1)-X)-(1.-X))
      Z2=.5*(A(1)-X)*AMXL2-AMXL1-.5*(1.-X)*OMXL2+OMXL1
      G=-(A(1)*A(1)*(.5*ALOG(A(1))-0.25)+0.25)/(1.-A(1))
      H=(0.5*(1.-A(1))**2*ALOG(1.-A(1))-0.25*(1.-A(1))**2)/(1.-A(1))+G
      GO TO 390
  330 Z=-.5*(1.-X)**2*ALOG(1.-X)+0.25*(1.-X)**2
      Z1=-.5*(1.-X)*(-ALOG(1.-X)-1.)+.5*(1.-X)*ALOG(1.-X)-.5*(1.-X)
      Z2=-ALOG(1.-X)-0.5
      G=-(A(1)**2*(.5*ALOG(A(1))-0.25)+0.25)/(1.-A(1))
      H=(0.5*(1.-A(1))**2*ALOG(1.-A(1))-0.25*(1.-A(1))**2)/(1.-A(1))+G
      GO TO 390
  340 G=-.25
      GO TO 390
  350 CONTINUE
      G=-(A(1)**2*(.5*ALOG(A(1))-0.25)+0.25)/(1.-A(1))
      H=(0.5*(1.-A(1))**2*ALOG(1.-A(1))-0.25*(1.-A(1))**2)/(1.-A(1))+G
      Z=.5*(A(1)-1.)**2*ALOG(ABS(A(1)-1.))-0.25*(A(1)-1.)**2
      Z1=-(A(1)-1.)*ALOG(ABS(A(1)-1.))
      Z2=-1.E10
      GO TO 390
  360 G=0.0
      H=G
      Q=G
      Z=-(1.-X)*ALOG(1.-X)
      Z1=ALOG(1.-X)+1.
      Z2=-1./(1.-X)
      GO TO 390
  370 Z=0.0
      G=Z
      H=Z
      Q=Z
      Z1=-1.E10       
      Z2=-1.E10
      GO TO 390
  380 G=-.25
      H=-.5
      Q=1.0
      Z=-.25
      Z1=0.0
      Z2=-10.E10
      GO TO 390
  390 YCMB(L)=CLI(1)*(Z/(1.-Q*A(1))-XLL+G-H*X)/PI/(A(1)+1.)/2.
      XSV=X
      IF (X.LT.0.005) X=0.005
      TANTH(L)=CLI(1)*(Z1/(1.-Q*A(1))-1.-ALOG(X)-H)/PI/(A(1)+1.)/2.0
      X=XSV
      IF (IF6XA(L).EQ.1) TANTH(L)=-5.
      IF (X.GT.0.005) GO TO 400
      YCP2(L)=0.0
      GO TO 420
  400 CONTINUE
      IF (ABS(1.-X).GT.E) GO TO 410
      YCP2(L)=1./E
      GO TO 420
  410 PIA=PI*(A(1)+1.)*2.
      YCP2(L)=CLI(1)*(Z2/(1.-Q*A(1))-1./X)/PIA
  420 CONTINUE
C
C       MODIFIED CAMBERLINE OPTION
C
      IF (CAMBER.EQ.'6A-SERIES') GO TO 430
      GO TO 450
  430 YCMB(L)=YCMB(L)*0.97948
      TANTH(L)=TANTH(L)*0.97948
      YCP2(L)=YCP2(L)*0.97948
      IF (ABS(A(1)-.8).LT.E) GO TO 440
      WRITE(LUNOUT,840)
      IF (KON.EQ.3) KON=0
      GO TO 20
  440 CONTINUE
      IF (TANTH(L).LE.-.24521*CLI(1)) YCMB(L)=0.24521*CLI(1)*(1.-X)
      IF (TANTH(L).LE.-.24521*CLI(1)) YCP2(L)=0.0
      IF (TANTH(L).LE.-.24521*CLI(1)) TANTH(L)=-0.24521*CLI(1)
      IF (TANTH(L).LE.-.24521*CLI(1)) IF6XA(L)=1
  450 CONTINUE
      IF (ICKY.GT.1.AND.L.LT.ICKY) GO TO 320
      IF (ICKY.EQ.1) GO TO 470
      DO 460 J=2,ICKY
      YCMB(1)=YCMB(1)+YCMB(J)
      TANTH(1)=TANTH(1)+TANTH(J)
      YCP2(1)=YCP2(1)+YCP2(J)
  460 CONTINUE
  470 CONTINUE
      F=SQRT(1.+TANTH(1)**2)
      THP=YCP2(1)/F**2
  480 CONTINUE
      IF(X.GT.XM) GO TO 600
      IF(ABS(X-XM).LT.E) GO TO 600
  490 CONTINUE
      SINTH=TANTH(1)/F
      COSTH=1./F
      I=I+1
      XU(I)=X-Y*SINTH
      YU(I)=YCMB(1)+Y*COSTH
      XL(I)=X+Y*SINTH
      YL(I)=YCMB(1)-Y*COSTH
      XAU(I) = XU(I)
      YAU(I) = YU(I)
      XAL(I) = XL(I)
      YAL(I) = YL(I)
  500 CONTINUE
      XUC=XU(I)*CHD
      YUC=YU(I)*CHD
      XLC=XL(I)*CHD
      YLC=YL(I)*CHD
      IF (CMB.LE.1.E-6) GO TO 510
      YUP=0.0
      YLP=YUP
      IF (ABS(TANTH(1)).LT.1.E-6) GO TO 510
      YUP=(TANTH(1)*F+YP-TANTH(1)*Y*THP)/(F-YP*TANTH(1)-Y*THP)
      YLP=(TANTH(1)*F-YP+TANTH(1)*Y*THP)/(F+YP*TANTH(1)+Y*THP)
      YUPR(I) = YUP
      YLPR(I) = YLP
  510 CONTINUE
      IF (X .LE. 0.0975) FRAC = 0.25
      IF (X .LE. 0.0124) FRAC = 0.025
C     IF (X .LE.0.00225) FRAC = 0.025
      IF (CMB.GT.1.E-6) WRITE(LUNOUT,850) X,XU(I),YU(I),XUC,YUC,YUP,     
     +                             XL(I),YL(I),XLC,YLC,YLP 
      IF (CMB.LE.1.E-6) WRITE(LUNOUT,830) X,Y,YP,YPP,XC,YC
      X=X+FRAC*DX
      FRAC=1.0
      IF(ABS(X-XM).LT.E) GO TO 520
      IF(X.LT.XM) GO TO 240
C
C     PROFILE - X GE XM
C
      X=XM
  520 CONTINUE
      IF (PROFILE.EQ.'4-DIGIT   ') GO TO 530
      IF (PROFILE.EQ.'4-DIGITMOD') GO TO 540
      WRITE(LUNOUT,860) PROFILE
      GO TO 20
  530 Y=0.29690*SQRT(X)-0.12600*X-0.35160*X**2+0.28430*X**3-0.1015*X**4
      YP=.5*.2969/SQRT(X)-.126-2.*.3516*X+3.*.2843*X*X-4.*.1015*X**3
      YPP=-.5*.5*.2969/SQRT(X**3)-2.*.3516+2.*3.*.2843*X-3.*4.*.1015*X*X
      GO TO 550
  540 Y=D0+D1*(1.-X)+D2*(1.-X)**2+D3*(1.-X)**3
      YP=-D1-2.*D2*(1.-X)-3.*D3*(1.-X)**2
      YPP=2.*D2+6.*D3*(1.-X)
  550 CONTINUE
      Y=Y*TOC/.2
      YP=YP*TOC/.2
      YPP=YPP*TOC/.2
      XC=X*CHD
      YC=Y*CHD
      IF (CAMBER.EQ.'2-DIGIT   ') GO TO 560
      IF (CAMBER.EQ.'3-DIGIT   ') GO TO 570
      IF (CAMBER.EQ.'3-DIGITREF') GO TO 580
      IF (CAMBER.EQ.'6-SERIES  ') GO TO 590
      IF (CAMBER.EQ.'6A-SERIES ') GO TO 590
      WRITE(LUNOUT,820) CAMBER
      GO TO 20
  560 YCMB(1)=CMB*(2.0*CM*X-X*X)/CM**2
      TANTH(1)=2.*CMB*(1.-X/CM)/CM
      IF (X.GT.CM) YCMB(1)=CMB*(1.-2.*CM+2.*CM*X-X*X)/(1.-CM)**2
      IF (X.GT.CM) TANTH(1)=(2.*CM-2.*X)*CMB/(1.-CM)**2
      F=SQRT(1.+TANTH(1)**2)
      THP=-2.*CMB/CM**2/F**2
      IF (X.GT.CM) THP=-2.*CMB/(1.-CM)**2/F**2
      GO TO 610
  570 YCMB(1)=CMB*(X**3-3.*CM*X**2+CM**2*(3.-CM)*X)/6.
      TANTH(1)=CMB*(3.*X**2-6.*CM*X+CM**2*(3.-CM))/6.
      IF (X.GT.CM) YCMB(1)=CMB*CM**3*(1.-X)/6.
      IF (X.GT.CM) TANTH(1)=-CMB*CM**3/6.
      F=SQRT(1.+TANTH(1)**2)
      THP=CMB*(X-CM)/F**2
      IF (X.GT.CM) THP=0.0
      GO TO 610
  580 YCMB(1)=CMB*((X-CM)**3-K2OK1*(1-CM)**3*X-CM**3*X+CM**3)/6.
      TANTH(1)=CMB*(3.*(X-CM)**2-K2OK1*(1-CM)**3-CM**3)/6.
      IF (X.GT.CM) YCMB(1)=CMB*(K2OK1*(X-CM)**3-K2OK1*(1-CM)**3*X-CM**3*
     +X+CM**3)/6.
      IF (X.GT.CM) TANTH(1)=
     +    CMB*(3*K2OK1*(X-CM)**2-K2OK1*(1-CM)**3-CM**3)/6.
      F=SQRT(1.+TANTH(1)**2)
      THP=CMB*(X-CM)/F**2
      IF (X.GT.CM) THP=K2OK1*CMB*(X-CM)/F**2
      GO TO 610
  590 GO TO 310
  600 CONTINUE
  610 CONTINUE
      SINTH=TANTH(1)/F
      COSTH=1./F
      I=I+1
      XU(I)=X-Y*SINTH
      YU(I)=YCMB(1)+Y*COSTH
      XL(I)=X+Y*SINTH
      YL(I)=YCMB(1)-Y*COSTH
  620 CONTINUE
      XUC=XU(I)*CHD
      YUC=YU(I)*CHD
      XLC=XL(I)*CHD
      YLC=YL(I)*CHD
      IF (CMB.LE.1.E-6) GO TO 630
      YUP=0.0
      YLP=YUP
      IF (ABS(TANTH(1)).LT.1.E-10) GO TO 630
      YUP=TANTH(1)*(F+YP/TANTH(1)-Y*THP)/(F-YP*TANTH(1)-Y*THP)
      YLP=TANTH(1)*(F-YP/TANTH(1)+Y*THP)/(F+YP*TANTH(1)+Y*THP)
      YUPR(I) = YUP
      YLPR(I) = YLP
  630 CONTINUE
      IF (CMB.GT.1.E-6) WRITE(LUNOUT,850) X,XU(I),YU(I),XUC,YUC,YUP,
     >                                    XL(I),YL(I),XLC,YLC,YLP
      IF (CMB.LE.1.E-6) WRITE(LUNOUT,830) X,Y,YP,YPP,XC,YC
      X=X+DX
      XAU(I) = XU(I)
      YAU(I) = YU(I)
      XAL(I) = XL(I)
      YAL(I) = YL(I)
CCCCC IF (X.LE.1.0) GO TO 520
      IF (X.LE.1.000001) GO TO 520  ! Avoid losing the trailing edge point
      WRITE(LUNOUT,690) NAME
  632 CONTINUE
C      DO 633 II = 1, I, 4
C         I1     = II + 3
C         WRITE (LUN,870) (XAU(K),YAU(K),K=II,I1)
  633 CONTINUE
C      DO 634 II = 1, I, 4
C         I1     = II + 3
C         WRITE (LUN,870) (XAL(K),YAL(K),K=II,I1)
  634 CONTINUE
  635 GO TO 20
C
  640 FORMAT (I3/(8F10.0))
  645 FORMAT (I5)
  660 FORMAT (' PROFILE  ',A,'  CAMBER  ',A)
  680 FORMAT (8F10.0)
  690 FORMAT (/10X,A,10X,A/)
  700 FORMAT (' PROFILE PARAMETERS'/' T/C=',F10.5/' L.E.RADIUS=',F10.5,/
     1        ' BASIC X INTERVAL=',F10.5/' CHORD=',F10.5,/)
  710 FORMAT (' PROFILE PARAMETERS'/' T/C=',F10.5/' L.E.RADIUS=',F10.5,/
     1        ' BASIC X INTERVAL=',F10.5/' CHORD=',F10.5,/
     2        ' POSITION OF MAXIMUM THICKNESS, XM=',F10.5,/
     3        ' CONSTANT D1=',F10.5,/)
  720 FORMAT (' CAMBER LINE PARAMETERS',/' CAMBER(YCMAX) =',F10.5,/
     1        ' POSITION OF MAXIMUM CAMBER=',F10.5,/)
  730 FORMAT (' CAMBER LINE PARAMETERS',/' CAMBER PARAMETER K1=',F10.5,/
     1        ' POSITION OF ZERO CAMBER LINE CURVATURE=',F10.5,/)
  740 FORMAT (' CAMBER LINE PARAMETERS',/' CAMBER PARAMETER K1=',F10.5,/
     1        ' POSITION OF ZERO CAMBER LINE CURVATURE=',F10.5,/
     2        ' RATIO OF AFT TO FORWARD CAMBER LINE CURVATURE FACTOR,',
     3        ' K2OK1=',F10.5,/)
  750 FORMAT (' CAMBER LINE PARAMETERS',/7X,'CLI',9X,'A')
  760 FORMAT (2F10.3)
  770 FORMAT (' J=',I2,/' IERR=',I2,/' D1=',E13.6)
  780 FORMAT (' NO REAL ROOTS')
  790 FORMAT (' A0,1,2,3=',4F13.6,/' D0,1,2,3=',4F13.6,/' RC=',F13.3,//)
  800 FORMAT (9X,'X/C',10X,'Y/C',8X,'DY/DX',6X,'D2Y/DX2',22X,'X',12X,'Y'
     1        /)
  810 FORMAT (/,' UNCAMBERED',25X,'UPPER SURFACE VALUES',39X,
     1        'LOWER SURFACE VALUES',/5X,'X/C',17X,'XU/C',6X,'YU/C',5X,
     2        '  XU        YU      DYU/DXU',13X,'XL/C',4X,'YL/C',5X,
     3        '  XL     YL   DYL/DXL')
  820 FORMAT (' BAD HOLLERITH CAMBER SPECIFICATION',/A)
  830 FORMAT (4F13.6,10X,2F13.6)
  840 FORMAT ('  MODIFIED CAMBER LINE OPTION ONLY ALLOWED IF A=0.8')
  850 FORMAT (F10.5,10X,4F10.5,E11.2,6X,4F10.5,E11.2)
  860 FORMAT (' BAD HOLLERITH PROFILE SPECIFICATION',/A)
  870 FORMAT (8F10.6)

  999 CONTINUE
      OPEN (LUNSAV, FILE='airfoil.dat', STATUS='unknown')
      CALL PRWRIT ( LUNSAV, 200, NAME, I, I, XU, XL, YU, YL, 1, 1 ) 
      END
C+---------------------------------------------------------------------
C
      SUBROUTINE PRWRIT (LUN, MAXPTS, TITLE, NU, NL, XU, XL, YU, YL,
     >                   FORMAT, PRECISION)
C
C  ACRONYM: Program PRofile: WRITe one airfoil profile
C                   --       ----
C  PURPOSE: PRWRIT writes one airfoil profile to the indicated file, in
C           one of four formats - "standard",  wraparound (either way),
C           or "three-column" - described in PRREAD.
C
C           PRWRIT  may also be used to save camber/thickness distribu-
C           tions,  or to save second derivative information for reuse.
C           Both of these extended uses employ  "standard"  format (not
C           wrap-around), with suitable labeling.
C
C           This version allows for three types of precision,  prompted
C           by the need to retain more digits for manipulating the all-
C           important leading edge region effectively, and to deal with
C           large magnitudes such as when the coordinates are in milli-
C           meters, or when the "ordinates" are really 2nd derivatives.
C
C  METHOD:  Values are passed to PRWRIT as separate surfaces,  and  are
C           returned untouched.  Note that the wrap-around formats omit
C           one of the two leading edge points, which are assumed to be
C           the same point.
C           
C  PARAMETERS:
C   ARG     DIM   TYPE  I/O/S   DESCRIPTION
C   LUN      -      I     I     Logical unit for file being written to.
C   MAXPTS   -      I     I     Max. no. of points on any one surface.
C   TITLE    -     C*(*)  I     Variable-length title for profile.
C   NU,NL    -      I     I     Number of data points for upper and lower
C                               surfaces (NU > 0; NL=0 is OK if FORMAT=1;
C                               NU=NL=no. of pts. in 3-column format and
C                               the camber/thickness distributions - FORMATs
C                               4 & 5).
C   XU     MAXPTS   R     I     Upper surface abscissas.
C   XL     MAXPTS   R     I     Lower surface abscissas (if any).
C   YU,YL  MAXPTS   R     I     Corresponding ordinates, y" values, or
C                               camber/thickness values, accdg. to FORMAT.
C   FORMAT   -      I     I     Requested format for output where:
C                               = 1 means standard PROFILE format (coords.)
C                               = 2 means clockwise wrap-around format
C                               = 3 means counterclockwise wrap-around
C                               = 4 means 3-column format
C                               = 5 means 2nd derivatives (standard format)
C                               = 6 means camber/thickness (standard format)
C   PRECISION -     I     I     Controls number of digits in saved values:
C                               = 1 means "full" (see PURPOSE above);
C                               = 2 means "engineering" or "flow code" (F10.6);
C                               = 3 means "low" - appropriate for y" values.
C  FILES USED:
C   UNIT     I/O/S    DESCRIPTION
C   LUN        O      File (assumed open) to contain one or more datasets
C
C  ENVIRONMENT:    VAX/VMS FORTRAN
C
C  AUTHOR:         Leslie Collins, Informatics, Palo Alto, CA.
C
C  HISTORY:
C  12/23/82   LJC   Coding adapted from PRREAD.
C  10/30/83   DAS   Introduced alternative E-format so that PRWRIT
C                   can be used to save 2nd derivatives similarly.
C  07/09/84   LJC   Added writing of wraparound formats.
C  02/11/84   DAS   Handled thickness/camber - hard to avoid dupli-
C                   cation of code now.  But PRWRIT is still handy
C                   for getting this stuff out of the main program.
C  02/28/85   LJC   Added 3-column format. (FORMAT=4)
C  04/09/87   DAS   Values except X are written with F10.7 format
C                   instead of F10.6.   (May help difficulties in
C                   critical leading edge region; will not affect
C                   formatted 2-column reads, but may be a little
C                   misleading.)
C  04/23/87   DAS   Erred further in direction of more precision:
C                   use E format except on basically normalized
C                   data; go to F12.8 for normalized airfoils.
C  04/24/87   DAS   Retained F10.6 option after all for old flow
C                   codes (FORMAT=7; standard PROFILE format only.
C  04/27/87   DAS   Reluctantly introduced PRECISION argument after
C                   the above failed to handle FORMAT="flowcode" and
C                   COUNTERCLOCKWISE both.
C  11/03/87   DAS   Switched to G formats for "full" and "low" precision.
C  12/12/91   DAS   Had to "comment out" the text following NU, NL
C                   because of how RDXYZ works now.  (PROFILE uses
C                   PRWRIT to save 2nd derivatives, which are read back
C                   via RDXYZ.)
C
C----------------------------------------------------------------------

      IMPLICIT NONE

C  *  Arguments:

      CHARACTER
     >   TITLE * (*)
      INTEGER
     >   FORMAT, LUN, MAXPTS, NL, NU, PRECISION
      REAL
     >   XU (MAXPTS), XL (MAXPTS), YU (MAXPTS), YL (MAXPTS)

C  *  Local variables:

      INTEGER
     >   I
      CHARACTER
     >   FMT * 13

C  *  Execution:

      IF (PRECISION .EQ. 1) THEN        ! "Full" precision.
         FMT = '(2(1X,G15.7))'          ! Gw.d needs w >= d + 8 (basically).

      ELSE IF (PRECISION .EQ. 2) THEN   ! "Engineering" precision.
         FMT = '(2F10.6)     '          ! Traditional for many flow codes.

      ELSE                              ! "Low" precision.
         FMT = '(2(1X,G13.5))'          ! Appropriate for y" values (FORMAT=5).
      END IF

      IF (FORMAT .EQ. 4) FMT (2:2) = '3'

      WRITE (LUN, '(A)') TITLE

      IF (FORMAT .EQ. 1) THEN           ! Standard PROFILE format.

         WRITE (LUN, 1001) NU, 'Upper Coordinates'
         WRITE (LUN, FMT) (XU (I), YU (I), I = 1, NU)

         IF (NL .GT. 0) THEN
            WRITE (LUN, 1001) NL, 'Lower Coordinates'
            WRITE (LUN, FMT) (XL (I), YL (I), I = 1, NL)
         END IF            

      ELSE IF (FORMAT .EQ. 2) THEN      ! Wrap-around (lower surface first).

         WRITE (LUN, 1001) NU + NL - 1, 'Coordinates Clockwise'
         WRITE (LUN, FMT) (XL (I), YL (I), I = NL, 1, -1)
         WRITE (LUN, FMT) (XU (I), YU (I), I = 2, NU)

      ELSE IF (FORMAT .EQ. 3) THEN      ! Wrap-around (upper surface first).

         WRITE (LUN, 1001) NU + NL - 1, 'Coordinates Counter-clockwise'
         WRITE (LUN, FMT) (XU (I), YU (I), I = NU, 1, -1)
         WRITE (LUN, FMT) (XL (I), YL (I), I = 2, NL)

      ELSE IF (FORMAT .EQ. 4) THEN      ! Three-column format.

         WRITE (LUN, 1001) NU, 'Coordinates per surface'
         WRITE (LUN, FMT) (XU (I), YU (I), YL (I), I = 1, NU)

      ELSE IF (FORMAT .EQ. 5) THEN      ! 2nd derivatives (both surfaces).
                                        ! Use standard format.
         WRITE (LUN, 1001) NU, 'Upper 2nd Derivatives'
         WRITE (LUN, FMT) (XU (I), YU (I), I = 1, NU)

         WRITE (LUN, 1001) NL, 'Lower 2nd Derivatives'
         WRITE (LUN, FMT) (XL (I), YL (I), I = 1, NL)

      ELSE  ! FORMAT = 6: Camber and Thickness distributions in Standard format.

         WRITE (LUN, 1001) NU, 'Camber'
         WRITE (LUN, FMT) (XU (I), YU (I), I = 1, NU)

         WRITE (LUN, 1001) NL, 'Thickness'
         WRITE (LUN, FMT) (XL (I), YL (I), I = 1, NL)

      END IF        

      RETURN

C  *  Formats:

 1001 FORMAT (I4, ' ! ', A)

      END
