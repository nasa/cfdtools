C+-----------------------------------------------------------------------
C
      PROGRAM SIXSERIES
C
C    (formerly PROGRAM LADSON)
C
C  PURPOSE: To produce the ordinates for airfoils of any thickness,
C           thickness distribution, or camber in the NACA 6- and 6A-series.
C
C  NAMELIST "INPUTS":
C
C   VAR     DIM   TYPE    DESCRIPTION
C  NAME      -    C*80    Title desired on printed and plotted output 
C  SERIET    -     I      NACA airfoil family (series 63 thru 67 and
C                         63A thru 65A)
C  SERIEC    -     I      Camber line (series 63 thru 66 and 63A thru 65A)
C  TOC       -     R      Thickness-chord ratio of airfoil
C  RLE       -     R      Leading-edge radius may be entered if desired
C                         (not used in program)
C  CMBNMR    -     R      Number of mean lines to be summed (>=1.)
C  CHD       -     R      Model chord used for listing ordinates in dimens-
C                         ional units
C  CLI       10    R      Design lift coefficient; set to 0.0 for a symmetrical
C                         airfoil. (Additional coefficients for up to nine
C                         mean lines may be added.)
C   A        10    R      Mean line chordwise loading (use 0.8 for 6A-series)
C
C  ENVIRONMENT:  DEC VAX/VMS FORTRAN (this version)
C
C  REFERENCE:    NASA Technical Memorandum TM X-3069 (September, 1974) 
C                Charles Ladson and Cuyler Brooks
C                Langley Research Center
C  HISTORY:
C    November 1989   Liam Hardy/NASA Ames   Interpolation routine FTLUP replaced
C                                           by LCSFIT to eliminate non-standard
C                                           %LOC usage; finite difference
C                                           derivatives (FUNCTION DIF) replaced
C                                           with FDI2K.
C-------------------------------------------------------------------------------

      CHARACTER
     >   METHOD*1, NAME*80, SERIEC*3, SERIET*3

      INTEGER
     >   I, KON, MXCOMB, NEVAL, NX

      LOGICAL
     >   NEW

      PARAMETER
     >   (LUNRD=5, LUNOUT=6, LUNSAV=10, MXCOMB=10, NX=200)

      REAL
     >   CHD, YUPR(NX), YLPR(NX),
     >   XU(NX), XL(NX), YU(NX), YL(NX),
     >   XAU(NX), YAU(NX), XAL(NX), YAL(NX),
     >   XT(NX+1), YT(NX+1), YTP(NX+1), YTPP(NX+1),
     >   PHI(NX+1), EPS(NX+1), PSI(NX+1),
     >   CLI(MXCOMB), A(MXCOMB), TANTH0(MXCOMB), YCMB(MXCOMB),
     >   TANTH(MXCOMB), YCP2(MXCOMB),
     >   AYP(NX+1), AYPP(NX+1), AX(NX+1), AXUC(NX+1), AYUC(NX+1),
     >   AYUP(NX+1), AXLC(NX+1), AYLC(NX+1), AYLP(NX+1)

      NAMELIST /INPUTS/ A, CHD, CLI, CMBNMR, NAME, SERIEC, SERIET, TOC
                      

      E=0.1D-10
      PI=ATAN(1.) * 4.
      KON=0
      DX=0.01

C     Input parameters normalized by the chord (CHD)
C     TOC - T/C, Thickness, RLE - Leading edge radius, XM - X(YMAX)/CHORD
C     DX - Interval/Chord, CHD - Chord in desired units
C     CMBNMR - Number of mean lines to be summed (>=1.)

  100 READ ( LUNRD, NML=INPUTS, END=999 )
      ICKY=CMBNMR
      IF(ICKY.LT.1) ICKY=1
      ICKYP=ICKY+1
      DO 110 J=ICKYP,10
         CLI(J)=0.0
         A(J)=0.0
  110 CONTINUE
      IF6XA=0
      KON=KON+1
      FRAC=1.0
C     Print Inputs:
      WRITE(LUNOUT,NML=INPUTS)

C     Slope of camberline at origin, TANTH0:
      CLIS=CLI(1)
      AS=A(1)
      L=0
  200 CONTINUE
         L=L+1
         A(1)=A(L)
         CLI(1)=CLI(L)
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
         U=0.005
         V=-(A(1)-U)/ABS(A(1)-U)
         OMXL=(1.-U)*ALOG(1.-U)
         AMXL=(A(1)-U)*ALOG(ABS(A(1)-U))
         OMXL1=-ALOG(1.-U)-1.
         AMXL1=-ALOG(ABS(A(1)-U))+V
         OMXL2=1./(1.-U)
         AMXL2=-V/ABS(A(1)-U)
         IF(A(1).GE.E.AND.ABS(1.-A(1)).GE.E) THEN
            G=-(A(1)*A(1)*(.5*ALOG(A(1))-0.25)+0.25)/(1.-A(1))
            Q=1.0
            H=(0.5*(1.-A(1))**2*ALOG(1.-A(1))-0.25*(1.-A(1))**2)/
     >        (1.-A(1))+G
            Z=.5*(A(1)-U)*AMXL-.5*(1.-U)*OMXL-.25*(A(1)-U)**2+.25*
     >        (1.-U)**2
            Z1=.5*((A(1)-U)*AMXL1-AMXL-(1.-U)*OMXL1+OMXL+(A(1)-U)-
     >         (1.-U))
            Z2=.5*(A(1)-U)*AMXL2-AMXL1-.5*(1.-U)*OMXL2+OMXL1
         END IF
         IF(A(1).LT.E) THEN
            H=-.5
            Q=1.0
            Z1=U*ALOG(U)-.5*U-.5*(1.-U)*OMXL1+.5*OMXL-.5
         ELSE IF(ABS(A(1)-1.).LT.E) THEN
            H=0.0
            Q=0.0
            Z1=-OMXL1
         END IF
         TANTH0(L)=CLI(1)*(Z1/(1.-Q*A(1))-1.-ALOG(U)-H)/PI/
     >             (A(1)+1.)/2.0
         IF(L.LT.ICKY.AND.ICKY.GT.1)
     >GO TO 200

      IF(ICKY.NE.1) THEN
         DO 210 J=2,ICKY
            TANTH0(1)=TANTH0(1)+TANTH0(J)
  210    CONTINUE
      END IF

C     Slope of profile at origin, Upper and Lower:
      YP=10.E10
      YPP=10.E10
      IF(TANTH0(1).NE.0.) THEN
         YUP=-1./TANTH0(1)
         YLP=-1./TANTH0(1)
      ELSE
         YUP=0.
         YLP=0.
      END IF

C     First station aft of origin on uncambered profile:
      I=1
      X=.00025
C     Start loop for X increment:
  300 CONTINUE
C       Skip thickness computation after first pass:
      IF(I.GT.1) GO TO 550
C       Select series:
      IF (SERIET.EQ.'63')  CALL PHEP63  (PHI,EPS)
      IF (SERIET.EQ.'64')  CALL PHEP64  (PHI,EPS)
      IF (SERIET.EQ.'65')  CALL PHEP65  (PHI,EPS)
      IF (SERIET.EQ.'66')  CALL PHEP66  (PHI,EPS)
      IF (SERIET.EQ.'67')  CALL PHEP67  (PHI,EPS)
      IF (SERIET.EQ.'63')  CALL PHPS63  (PHI,PSI)
      IF (SERIET.EQ.'64')  CALL PHPS64  (PHI,PSI)
      IF (SERIET.EQ.'65')  CALL PHPS65  (PHI,PSI)
      IF (SERIET.EQ.'66')  CALL PHPS66  (PHI,PSI)
      IF (SERIET.EQ.'67')  CALL PHPS67  (PHI,PSI)
      IF (SERIET.EQ.'63A') CALL PHEP63A (PHI,EPS)
      IF (SERIET.EQ.'64A') CALL PHEP64A (PHI,EPS)
      IF (SERIET.EQ.'65A') CALL PHEP65A (PHI,EPS)
      IF (SERIET.EQ.'63A') CALL PHPS63A (PHI,PSI)
      IF (SERIET.EQ.'64A') CALL PHPS64A (PHI,PSI)
      IF (SERIET.EQ.'65A') CALL PHPS65A (PHI,PSI)
      RAT=1.0
      IT=0
      ACRAT=1.0

C       Loop start for thickness iteration:
  400 CONTINUE
         IT=IT+1
         WRITE(LUNOUT,1003)IT,RAT
         ACRAT=ACRAT*RAT
         YMAX=0.0
         DO 410 J=1,NX
            XT(J)=-2.0*COSH(PSI(J)*ACRAT)*COS(PHI(J)-EPS(J)*ACRAT)
            YT(J)= 2.0*SINH(PSI(J)*ACRAT)*SIN(PHI(J)-EPS(J)*ACRAT)
            IF(YT(J).GT.YMAX) XYM=XT(J)
            IF(YT(J).GT.YMAX) YMAX=YT(J)
  410    CONTINUE

C        Estimate first and second derivatives by finite differencing:

         CALL FD12K (NX, XT, YT, YTP, YTPP, YTPP)

C        Estimate location of maximum thickness:

         XTP=1.0
         DO 420 J=3, NX
            IF (YTP(J).LT.0.0.AND.YTP(J-1).GE.0.0)
     >         XTP=XT(J-1)+YTP(J-1)*(XT(J)-XT(J-1))/(YTP(J-1)-YTP(J))
  420    CONTINUE

         CALL LCSFIT (NX+1, XT, YT, .TRUE., 'B', 1, XTP, YM, YM)
         XO=XT(1)
         XL(1)=XT(NX)
         TR=2.*YM/(XL(1)-XO)
         RAT=TOC/TR
         SF=RAT
         IF(TOC.GT.E .AND. ABS(RAT-1.0).GT.0.0001 .AND. IT.LE.10)
     >GO TO 400

      IF(I.EQ.1) THEN
         DO 500 J=1,NX+1
            XT(J)=(XT(J)-XO)/(XL(1)-XO)
C           Scale linearly to exact thickness:
            YT  (J)=SF*YT  (J)/(XL(1)-XO)
            YTP (J)=SF*YTP(J)
            YTPP(J)=SF*YTPP(J)*(XL(1)-XO)
  500    CONTINUE
      END IF
      XTP=(XTP-XO)/(XL(1)-XO)
      YMAX=YMAX*SF/(XL(1)-XO)
      YM=YM*SF/(XL(1)-XO)
      XYM=(XYM-XO)/(XL(1)-XO)
      XL(1)=0.0
      IF(TOC.GT.E) THEN
C        Fit tilted ellipse at eleventh profile point:
         CN=2.*YTP(11)-YT(11)/XT(11)+0.1
         AN=XT(11)*(YTP(11)*XT(11)-YT(11))/(XT(11)*
     >      (2.*YTP(11)-CN)-YT(11))
         BN=SQRT((YT(11)-CN*XT(11))**2/(1.-(XT(11)-AN)**2/AN**2))
         DO 510 J=1,10
            YT(J)=BN*SQRT(1.-(XT(J)-AN)**2/AN**2)+CN*XT(J)
            IF(XT(J).LE.E) GO TO 510
            YTP(J)=BN**2*(AN-XT(J))/AN**2/(YT(J)-CN*XT(J))+CN
            YTPP(J)=-BN**4/AN**2/(YT(J)-CN*XT(J))**3
  510    CONTINUE
         RNP=BN**2/AN
         IF (I.EQ.1) WRITE(LUNOUT,1004)XYM,YMAX,XTP,YM,XT(11),YT(11),
     >               YTP(11),RNP,RAT,ACRAT,IT
      END IF
      X=0.0
      ALI=ABS(CLI(1))
      IF(ALI.LE.E.AND.ICKY.EQ.1) THEN
C       Print uncambered column headings:
        WRITE(LUNOUT,1001)
      ELSE
C       Print cambered column headings:
        WRITE(LUNOUT,1002)
      END IF

      X=0.00025
      XL(1)=0.0
  550 CONTINUE
      YUPR(I) = YUP
      YLPR(I) = YLP

C     Interpolate for thickness and derivatives at desired values of X:
      NEW = .TRUE.
      METHOD = 'B'
      NEVAL = 1.
      CALL LCSFIT (NX, XT, YT, NEW, METHOD, NEVAL, X, Y, Y) 
      CALL LCSFIT (NX, XT, YTP, NEW, METHOD, NEVAL, X, YP, YP)
      CALL LCSFIT (NX, XT, YTPP, NEW, METHOD, NEVAL, X, YPP, YPP)

C     Compute camberline:
      A(1)=AS
      CLI(1)=CLIS
      L=0
  600 CONTINUE
         L=L+1
         A(1)=A(L)
         CLI(1)=CLI(L)
         XC=X*CHD
         YC=Y*CHD
         XLL=X*ALOG(X)
         Q=1.0
         IF(ABS(1.-A(1)).LT.E.AND.ABS(1.-X).LT.E) THEN
            G=0.0
            H=0.0
            Q=0.0
            Z=0.0
            Z1=-10.E10
            Z2=-10.E10
         ELSE IF(A(1).LT.E.AND.(1.-X).LT.E) THEN
            G=-.25
            H=-.5
            Q=1.0
            Z=-.25
            Z1=0.0
            Z2=-10.E10
         ELSE IF(ABS(A(1)-X).LT.E) THEN
            Z=-.5*(1.-X)**2*ALOG(1.-X)+0.25*(1.-X)**2
            Z1=-.5*(1.-X)*(-ALOG(1.-X)-1.)+.5*(1.-X)*ALOG(1.-X)-.5*
     >         (1.-X)
            Z2=-ALOG(1.-X)-0.5
            G=-(A(1)**2*(.5*ALOG(A(1))-0.25)+0.25)/(1.-A(1))
            H=(0.5*(1.-A(1))**2*ALOG(1.-A(1))-0.25*(1.-A(1))**2)/
     >        (1.-A(1))+G
         ELSE IF(ABS(1.-X).LT.E) THEN
            G=-(A(1)**2*(.5*ALOG(A(1))-0.25)+0.25)/(1.-A(1))
            H=(0.5*(1.-A(1))**2*ALOG(1.-A(1))-0.25*(1.-A(1))**2)/
     >        (1.-A(1))+G
            Z=.5*(A(1)-1.)**2*ALOG(ABS(A(1)-1.))-0.25*(A(1)-1.)**2
            Z1=-(A(1)-1.)* ALOG(ABS(A(1)-1.))
            Z2=-10.E10
         ELSE IF(ABS(A(1)-1.).LT.E) THEN
            G=0.0
            H=0.0
            Q=0.0
            Z=-(1.-X)*ALOG(1.-X)
            Z1=ALOG(1.-X)+1.
            Z2=-1./(1.-X)
         ELSE
            V=-(A(1)-X)/ABS(A(1)-X)
            OMXL=(1.-X)*ALOG(1.-X)
            AMXL=(A(1)-X)*ALOG(ABS(A(1)-X))
            OMXL1=-ALOG(1.-X)-1.
            AMXL1=-ALOG(ABS(A(1)-X))-1.
            OMXL2=1./(1.-X)
            AMXL2=1./(A(1)-X)
            Z=.5*(A(1)-X)*AMXL-.5*(1.-X)*OMXL-.25*(A(1)-X)**2+.25*(1.-X)
     >        **2
            Z1=.5*((A(1)-X)*AMXL1-AMXL-(1.-X)*OMXL1+OMXL+(A(1)-X)-
     >         (1.-X))
            Z2=.5*(A(1)-X)*AMXL2-AMXL1-.5*(1.-X)*OMXL2+OMXL1
            IF(A(1).LE.E) THEN
               G=-.25
               H=-.5
            ELSE
               G=-(A(1)*A(1)*(.5*ALOG(A(1))-0.25)+0.25)/(1.-A(1))
               H=(0.5*(1.-A(1))**2*ALOG(1.-A(1))-0.25*(1.-A(1))**2)/
     >           (1.-A(1))+G
            END IF
         END IF

         YCMB(L)=CLI(1)*(Z/(1.-Q*A(1))-XLL+G-H*X)/PI/(A(1)+1.)/2.
         XSV=X
         IF(X.LT.0.005) X=0.005
         TANTH(L)=CLI(1)*(Z1/
     >   (1.-Q*A(1))-1.-ALOG(X)-H)/PI/(A(1)+1.)/2.0
         X=XSV
         IF(IF6XA.EQ.1) TANTH(L)=-2.0
         IF(X.LE.0.005) THEN
               YCP2(L)=0.0
            ELSE IF (ABS(1.-X).LE.E) THEN
                  YCP2(L)=1./E
               ELSE
                  PIA=PI*(A(1)+1.)*2
                  YCP2(L)=CLI(1)*(Z2/(1.-Q*A(1))-1./X)/PIA
         END IF  
C        Modified camberline option:
         IF(SERIEC.EQ.'63A' .OR. SERIEC.EQ.'64A' .OR. SERIEC.EQ.'65A')
     >      THEN
            YCMB(L)=YCMB(L)*0.97948
            TANTH(L)=TANTH(L)*0.97948
            IF(TANTH(L).LE.-.24521*CLI(1)) YCMB(L)=0.24521*CLI(1)*
     >      (1.-X)
            IF(TANTH(L).LE.-.24521*CLI(1)) YCP2(L)=0.0
            IF(TANTH(L).LE.-.24521*CLI(1)) TANTH(L)=-0.24521*CLI(1)
            IF(TANTH(L).LE.-.24521*CLI(1)) IF6XA=1
         END IF
         IF(ICKY.GT.1.AND.L.LT.ICKY)
     > GO TO 600

      IF(ICKY.EQ.1) GO TO 620
      DO 610 J=2,ICKY
         YCMB (1)=YCMB (1)+YCMB (J)
         TANTH(1)=TANTH(1)+TANTH(J)
         YCP2 (1)=YCP2 (1)+YCP2 (J)
  610 CONTINUE
  620 CONTINUE
      F=SQRT(1.+TANTH(1)**2)
      THP=YCP2(1)/F**2
      SINTH=TANTH(1)/F
      COSTH=1./F
C     Camberline and derivatives computed:
      I=I+1
C     Combine thickness distributuion and camberline:
      XU(I)=X-Y*SINTH
      YU(I)=YCMB(1)+Y*COSTH
      XL(I)=X+Y*SINTH
      YL(I)=YCMB(1)-Y*COSTH
      IF (X.GE. .815) THEN
         IF(SERIET.EQ.'63A' .OR. SERIET.EQ.'64A' .OR. SERIET.EQ.'65A') 
     >   THEN
            IF (X.LE. .825) THEN
               X2 = 1.0
               X1 = XU(I)
               Y2 = 0.
               Y1 = YU(I)
               S1 = (Y2 - Y1) / (X2 - X1)
               S2 = (Y2 - YL(I)) / (X2 - XL(I))
               B1 = Y2 - S1 * X2
               B2 = Y2 - S2 * X2
            END IF
         YU(I) = S1 * XU(I) + B1
         YL(I) = S2 * XL(I) + B2
         END IF
      END IF
C       Multiply by chord:
      XUC=XU(I)*CHD
      YUC=YU(I)*CHD
      XLC=XL(I)*CHD
      YLC=YL(I) *CHD
      IF(ALI.GT.E.OR.ICKY.NE.1) THEN
C        Find local slope of cambered profile:
         YUP=(TANTH(1)*F+YP-TANTH(1)*Y*THP)/(F-YP*TANTH(1)-Y*THP)
         YLP=(TANTH(1)*F-YP+TANTH(1)*Y*THP)/(F+YP*TANTH(1)+Y*THP)
      END IF

C     Find X increment:
      IF(X.LE.0.0975) FRAC=0.25
      IF(X.LE.0.01225) FRAC=0.025
C     Store profile in appropriate arrays:
      AXUC(I) = XUC
      AYUC(I) = YUC
      IF (ALI.GE.E .OR. ICKY.NE.1) THEN
         AX(I) = X
         AXLC(I) = XLC
         AYLC(I) = YLC
      END IF

C     Increment X and return to start of X loop:
      X=X+FRAC*DX
      FRAC=1.0
      XAU(I) = XUC
      YAU(I) = YUC
      XAL(I) = XLC
      YAL(I) = YLC
      IF(X.LE.1.0) GO TO 300

C     Calculate derivatives of final coordinates; tabulate and save results.

      OPEN (LUNSAV, FILE='airfoil.dat', STATUS='UNKNOWN')

      IF (ALI.GE.E .OR. ICKY.NE.1) THEN
         CALL FD12K (I, XU, YU, AYUP, AYUP, AYUP)
         CALL FD12K (I, XL, YL, AYLP, AYLP, AYLP)
         WRITE (LUNOUT,1005) (AX(M), XU(M), YU(M), AXUC(M), AYUC(M),
     >                       AYUP(M), XL(M), YL(M), AXLC(M), AYLC(M),
     >                       AYLP(M), M = 1, I)
C        Save coordinates in standard "PROFILE" form:
         CALL PRWRIT (LUNSAV, NX, NAME, I, I, XU, XL, YU, YL, 1, 1)
      ELSE
         CALL FD12K (I, XU, YU, AYUP, AYPP, AYPP)
         WRITE (LUNOUT,1006) (XU(M),YU(M),AYP(M),AYPP(M),AXUC(M),
     >                       AYUC(M), M = 1, I)
         CALL PRWRIT (LUNSAV, NX, NAME, I, I, XU, XL, YU, YL, 1, 1)
      END IF

C     Return to read for next case:
      GO TO 100

  999 STOP ' '

C     Formats:

 1001 FORMAT(9X, 'X/C', 10X, 'Y/C', 9X, 'DY/DX', 9X, 'D2Y/DX2',
     >       22X, 'X', 14X,'Y'/)
 1002 FORMAT(/, 'UNCAMBERED                         UPPER SURFACE VALUES
     >                                     LOWER SURFACE VALUES '/
     >       6X, 'X/C',11X, 'XU/C', 7X, 'YU/C', 8X, 'XU', 9X, 'YU',
     >       6X, 'DYU/DXU', 12X, 'XL/C', 7X, 'YL/C', 8X, 'XL', 9X, 'YL',
     >       6X, 'DYL/DXL')
 1003 FORMAT(' RAT('I2') = ', F10.5)
 1004 FORMAT(//
     >       ' PEAK IS AT X/C = ', F10.6/
     >       ' MAXIMUM Y/C = ', F10.6/
     >       ' SLOPE CHANGES SIGN AT X/C, Y/C = ', 2F10.6/
     >       ' X/C FIT OF ELLIPSE = ', F10.6/
     >       ' Y/C FIT OF ELLIPSE = ', F10.6/
     >       ' SLOPE FIT OF ELLIPSE = ', F10.6/
     >       ' RADIUS AT ORIGIN OF ELLIPSE THRU XT(11)/C, YT(11)/C = ',
     >         F10.6/
     >       ' RATIO OF T/C INPUT TO T/C COMPUTED = ', F10.6/
     >       ' CUMULATIVE SCALING OF EPS, PSI = ', F10.6/
     >       ' NUMBER OF ITERATIONS = ', I10)
 1005 FORMAT (F10.6, 4X, 4F11.6, E13.4, 4X, 4F11.6, E13.4)
 1006 FORMAT (F13.6, F15.6, 2E15.6, 10X, 2F15.6)

      END

C+----------------------------------------------------------------------
C
      SUBROUTINE FD12K (N, X, F, FP, FPP, FK)
C
C  PURPOSE: FD12K returns estimates of 1st and 2nd derivatives and of
C           curvature, by finite differencing, for each of the points
C           (X(I), F(I)), I = 1 : N.  The abscissas are assumed to be
C           nonuniform, and they must be monotonic.
C
C           This routine combines calls to FDCNTR, FD1SID, FDCURV for
C           the common case of needing results for N >= 2 points.
C
C           If (say) curvature is wanted at a single point only, call
C           either FDCNTR or FD1SID and FDCURV directly.
C
C  INPUTS:  X(*) & F(*) are N-vectors defining some curve in 2-space.
C           For N > 2, the 3-pt formulas are used for all I (with the
C                      one-sided forms used at the end-points).
C           For N = 2, the 2-pt formulas are used.
C
C  OUTPUTS: FP, FPP, FK are N-vectors representing 1st and 2nd deriv-
C           atives and curvature respectively.  These are assigned in
C           reverse order (FK, FPP, FP) so that a call such as
C
C                     CALL FD12K (N, X, Y, YP, YP, YP)
C
C           can be used if just 1st derivatives are desired, to avoid
C           declaring storage for FPP and FK. (Similarly for the case
C           when 1st and 2nd derivatives are desired but curvature is
C           not. The unnecessary arithmetic in these cases is consid-
C           ered preferable to another argument and extra logic.)
C
C  METHOD:  Central differencing is used at all interior points, with
C           one-sided 3-point formulas used at each end-point.
C
C           The curvature formula is safeguarded against overflow  in
C           the presence of large 1st derivatives.  The basic formula
C           used here is:
C
C               FK (I) = FPP (I) / (1 + FP(I) ** 2) ** 3/2
C
C           Note that if X is not necessarily monotonic, curvature is
C           defined as
C
C               CURVATURE = (X" ** 2  +  Y" ** 2) ** 1/2
C
C           where " means 2nd derivative with respect to  arc-length.
C           See modules CURV2D and CURV3D for these parametric cases.
C
C  NOTES:   1. Finite differencing errors can be large if the delta-X
C              values are too small,  especially if the precision  in
C              the function values is less than full.
C           2. Nevertheless, finite differences have been observed to
C              behave better than the analytic derivatives of splines
C              in airfoil geometry applications.
C
C  EXTERNALS:
C           FDCNTR modularizes the central 3-point formulas for first
C                  and second derivatives by finite differences.
C           FDCURV modularizes the curvature formula (safe-guarded).
C           FD1SID modularizes the 1-sided forward and backward 3-pt.
C                  formulas for first and second derivatives.
C
C  HISTORY:
C           09/15/83   DAS   Initial implementation (interior pts only).
C           12/27/85   DAS   End points are handled by FD1SID now.
C           09/18/87   DAS   The N=2 case is handled now.
C           08/21/89   DAS   Formulation revised to use separate dF/dX
C                            terms instead of a common denominator.
C           08/17/91   DAS   Introduced FDCNTR and FDCURV when it was
C                            found that FD12K did not lend itself to
C                            application to one point at a time.
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT   NONE

C     Arguments:

      INTEGER    N
      REAL       X (N), F (N), FP (N), FPP (N), FK (N)

C     Local constants:

      REAL       ZERO
      PARAMETER (ZERO = 0.E+0)

C     Local variables:

      INTEGER    I
      REAL       FPI, FPPI

C     Procedures:

      EXTERNAL   FDCNTR, FDCURV, FD1SID

C     Execution:

C     Assign values in the order  curvature, f", f'  so that the
C     user can pass the same array for, say, curvature and f" if
C     curvature is not wanted:

      IF (N .EQ. 2) THEN

         FK (1)  = ZERO
         FPP (1) = ZERO
         FP (1)  = (F (2) - F (1)) / (X (2) - X (1))
         FK (2)  = ZERO
         FPP (2) = ZERO
         FP (2)  = FP (1)

      ELSE

C        Forward 3-pt. differencing for the first point:

         CALL FD1SID (1, 1, X, F, FPI, FPPI)
         CALL FDCURV (FPI, FPPI, FK (1))
         FPP (1) = FPPI
         FP (1)  = FPI
      
C        Central 3-pt. differencing for the bulk of the points:

         DO 20 I = 2, N-1
            CALL FDCNTR (I, X, F, FPI, FPPI)
            CALL FDCURV (FPI, FPPI, FK (I))
            FPP (I) = FPPI
            FP (I)  = FPI
   20    CONTINUE

C        Backward differencing for the last point:

         CALL FD1SID (N, -1, X, F, FPI, FPPI)
         CALL FDCURV (FPI, FPPI, FK (N))
         FPP (N) = FPPI
         FP (N)  = FPI
      END IF

      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE FDCNTR (I, X, F, FP, FPP)
C
C  PURPOSE: FDCNTR returns central 3-point finite-difference estimates
C           of the 1st and 2nd derivatives at the point  (X(I), F(I)).
C           Use FD1SID for the end-point cases.
C
C  ARGS:    Obvious from PURPOSE.
C
C  METHOD:  FPP is computed first, in case only FP is desired, so that
C           the same item may be passed for both arguments.
C
C  HISTORY: 08/17/91  DAS  FDCNTR adapted from FD12K's in-line code,
C                          as for FD1SID, for the case of a single I
C                          at a time (which FD12K can't do).
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT   NONE

C     Arguments:

      INTEGER    I
      REAL       X (*), F (*), FP, FPP

C     Local constants:

      REAL       ONE
      PARAMETER (ONE = 1.E+0)

C     Local variables:

      REAL     DEL1, DEL2, DIV, DX1, DX2, W

C     Execution:

      DX1  =  X (I) - X (I-1)
      DEL1 = (F (I) - F (I-1)) / DX1
      DX2  =  X (I+1) - X (I)
      DEL2 = (F (I+1) - F (I)) / DX2
      DIV  = ONE / (DX1 + DX2)
      W    = DX2 * DIV
      FPP  = (DEL2 - DEL1) * (DIV + DIV)
      FP   = W * DEL1 + (ONE - W) * DEL2

      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE FDCURV (DYDX, D2YDX2, KAPPA)
C
C  PURPOSE: FDCURV (Finite-Difference CURVature estimate) returns a
C           safe-guarded value for curvature at one point on the
C           curve Y = Y (X) given the 1st and 2nd derivatives at
C           that point, using the formula
C
C               KAPPA = D2YDX2 / (1 + DYDX ** 2) ** 3/2
C
C           The sign of KAPPA is clearly that of the 2nd derivative.
C           The derivatives could well be obtained from a spline, but
C           experience shows finite differencing can be preferable.
C
C           See modules CURV2D and CURV3D for the parametric cases.
C
C  ARGUMENTS:  Obvious from the description.  KAPPA is REAL.
C
C  HISTORY: 08/17/91  Derived FDCURV from FD12K, along with FDCNTR
C                     when it was found that FD12K did not lend
C                     itself to application to one point at a time.
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT   NONE

C     Arguments:

      REAL       DYDX, D2YDX2, KAPPA

C     Local constants:

      REAL       DYDXMAX, ONE
      PARAMETER (DYDXMAX = 1.E+10, ONE = 1.E+0)

C     Local variables:

      REAL       TERM

C     Execution:

      TERM = ONE + (MIN (ABS (DYDX), DYDXMAX)) ** 2
      KAPPA = D2YDX2 / (TERM * SQRT (TERM))

      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE FD1SID (I, INC, X, F, FP, FPP)
C
C  PURPOSE: FD1SID returns one-sided 3-point finite-difference estimates
C           of the 1st and 2nd derivatives at the point  ( X(I), F(I) ).
C           If INC = 1, points I, I+1, I+2 are used,  while if INC = -1,
C           points  I-2, I-1, I are used. The abscissas need not be uni-
C           formly spaced.  
C
C  ARGS:    Obvious from PURPOSE.
C
C  METHOD:  FPP is computed first,  in case only FP is desired,  so that
C           the same item may be passed for both arguments. The formula-
C           tion is similar to that of central differencing - see FD12K.
C
C  HISTORY: 12/27/85  DAS  Initial implementation  (prompted by the need
C                          to approximate an airfoil leading edge with a
C                          cubic having specified slope at one end).
C           08/21/89  DAS  Formulation revised as for centrals (FD12K).
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT   NONE

C     Arguments:

      INTEGER    INC, I
      REAL       X (*), F (*), FP, FPP

C     Local constants:

      REAL       ONE
      PARAMETER (ONE = 1.E+0)

C     Local variables:

      INTEGER  I1, I2
      REAL     DEL1, DEL2, DIV, DX1, DX2, W

C     Execution:

C     The minus signs take care of themselves for the backward case.

      I1   = I  + INC
      I2   = I1 + INC
      DX1  =  X (I1) - X (I)
      DEL1 = (F (I1) - F (I)) / DX1
      DX2  =  X (I2) - X (I1)
      DEL2 = (F (I2) - F (I1)) / DX2
      DIV  = ONE / (DX1 + DX2)
      W    = -DX1 * DIV
      FPP  = (DEL2 - DEL1) * (DIV + DIV)
      FP   = W * DEL2 + (ONE - W) * DEL1

      RETURN
      END
C+----------------------------------------------------------------------
C
      SUBROUTINE LCSFIT (NDATA, X, Y, NEW, METHOD, NEVAL, XEVAL, YEVAL,
     &   YPEVAL)
C
C     Two-liner:  Storage-efficient local cubic spline fit (2-space)
C     ----------  (monotonic and piecewise linear options too)
C
C     Description and usage:
C     ----------------------
C
C        LCSFIT is the non-parametric analog of PLSFIT (parametric).
C     It is intended for spline applications which do not require the
C     spline coefficients as output.  It is efficient for repeated
C     calls with the same data, so repeated use with NEVAL = 1 may be
C     preferable to storing vectors of results.
C
C        LCSFIT offers monotonic spline and piecewise linear options
C     also.  And it returns an interpolated first derivative along
C     with the function value.  (The second derivative is omitted
C     because Y" is not guaranteed to be continuous by local methods.)
C
C        See PLSFIT for more details on local methods.  As with most
C     numerical methods, scaling of the data to the unit interval (and
C     unscaling of the result) is recommended to avoid unnecessary
C     effects attributable to the data units.  Utilities GETSCALE and
C     USESCALE from the present authors are appropriate.  The data
C     abscissas should be distinct and either ascending or descending.
C     PROTECT is available to check this.  Extrapolation is permitted
C     (mainly in case of round-off; it is normally inadvisable).
C
C        The CSFIT/CSEVAL or CSDVAL pair are probably preferable if
C     efficiency is not an issue, since CSFIT gives Y" continuity.
C
C     Arguments:
C     ----------
C
C     Name    Type/Dimension  I/O/S  Description
C     NDATA   I               I      Length of X, Y input data arrays.
C
C     X,      R (NDATA)       I      Input data coordinates.  The Xs
C     Y                              must be distinct and monotonic,
C                                    either ascending or descending.
C                                    (No check here.) 
C
C     NEW     L               I      If control flag NEW is .TRUE., the
C                                    search for a bracket starts from
C                                    scratch, otherwise locally-saved
C                                    search and fit information will be
C                                    assumed to be correct. If calling
C                                    LCSFIT from within a loop, set
C                                    NEW = .FALSE. after the first call.
C
C     METHOD   C*1            I      (Uppercase) Type of fit to be used:
C                                    'M' means Monotonic piecewise cubics;
C                                    'B' means non-monotonic "Bessel"-type
C                                        piecewise cubics (looser fit);
C                                    'L' means piecewise Linear fit;
C                                    'C' means Cyclic (periodic) end
C                                        conditions: loose fit assumed.
C
C     NEVAL   I               I      Number of interpolations requested.
C                                    NEVAL >= 1.  One call per result
C                                    (NEVAL = 1) may save storage, and is
C                                    not too inefficient as long as NEW
C                                    is set to .FALSE. after the first.
C
C     XEVAL   R (NEVAL)       I      Abscissa(s) to interpolate to.  These
C                                    are normally in the data range, but
C                                    extrapolation - probably due to
C                                    round-off - is not prevented.
C
C     YEVAL   R (NEVAL)       O      Interpolated function value(s).
C
C     YPEVAL  R (NEVAL)       O      Interpolated 1st derivative value(s).
C                                    Pass the same storage as for YEVAL
C                                    if no derivatives are required.
C
C     Significant local variables:
C     ----------------------------
C
C     MEMORY         Indicates that coefficients are correct for the
C                    current point.
C
C     H, DEL         Delta X and forward difference derivative arrays.
C
C     B, C, D        Coefficients of cubic on the bracketing interval.
C
C     Procedures:
C     -----------
C
C     INTERVAL  1-D "interpolation" search.
C     BESSEL    First derivative (central 3-point formula).
C     BRODLIE   First derivative (central), adjusted for monotonicity.
C     BUTLAND   First derivative (non-central), adjusted for monotonicity.
C     THREEPT   First derivative (non-central 3-point formula).
C
C     Environment:  FORTRAN 90
C     ------------
C
C     Error handling:  None
C     ---------------
C
C     Notes:
C     ------
C
C     (1)  Since many of the calculations must be repeated at both ends
C          of an interval, the various finite difference quantities used
C          are stored as arrays. The following "map" of a typical interior
C          interval and its neighbors should help in understanding the
C          notation.  The local array indices are all numbered relative
C          to the left-hand end of the interval which brackets the point
C          to be evaluated.
C
C                                  LEFT       RIGHT
C
C          Point         -1          0         +1          +2
C
C          Data           x -------- x -------- x --------- x
C
C          Interval      -1          0         +1
C
C
C     Author: Robert Kennelly, Sterling Software/NASA Ames  (PLSFIT)
C     -------
C
C     History:
C     --------
C
C     27 Feb. 1987  R.A.Kennelly  Initial implementation of PLSFIT.
C     23 Aug. 1989  D.A.Saunders  LCSFIT adapted as non-parametric form,
C                                 for embedding in other utilities where
C                                 minimizing work-space is desirable.
C     20 June 1991    "    "      THREEPT (monotonic) renamed BUTLAND;
C                                 THREEPT (pure 3-pt. formula) now used
C                                 for nonmonotonic end-point handling;
C                                 METHOD='C' case belatedly added, as
C                                 needed by PLSINTRP for closed curves.
C     23 July 1991    "    "      The tests for being in the same interval
C                                 as before were not allowing for the
C                                 descending-Xs case.
C     06 May  1998    "    "      Minor Fortran 90 updates.
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER,   INTENT (IN)  :: NDATA, NEVAL
      REAL,      INTENT (IN)  :: X (NDATA), Y (NDATA), XEVAL (NEVAL)
      REAL,      INTENT (OUT) :: YEVAL (NEVAL), YPEVAL (NEVAL)
      LOGICAL,   INTENT (IN)  :: NEW
      CHARACTER, INTENT (IN)  :: METHOD * 1

C     Local constants:

      REAL, PARAMETER :: ZERO = 0., ONE = 1., TWO = 2., THREE = 3.

C     Local variables:

      LOGICAL
     &   CYCLIC, LINEAR, MEMORY, MONO
      INTEGER
     &   IEVAL, J, K, LEFT, RIGHT
      REAL
     &   ARROW, B (0:1), C, DELY (-1:1), D, DX, H (-1:1), XBYARROW, XE

C     Procedures:

      REAL, EXTERNAL ::
     &   BESSEL, BRODLIE, BUTLAND, THREEPT

C     Storage:

      SAVE
     &   ARROW, B, C, D, LEFT, RIGHT

C     Execution:
C     ----------

      MONO   = METHOD == 'M'
      CYCLIC = METHOD == 'C'
      LINEAR = METHOD == 'L'

      IF (CYCLIC) THEN
         IF (Y (NDATA) /= Y (1)) STOP 'LCSFIT: End points must match.'
      END IF

C     Initialize search or avoid it if possible:

      IF (NEW) THEN
         MEMORY = .FALSE.
         ARROW  = SIGN (ONE, X (2) - X (1))
         LEFT   = 1
      END IF

      IEVAL = 1
      XE = XEVAL (1)
      XBYARROW = XE * ARROW

      IF (.NOT. NEW) THEN
      
C        We can save a lot of time when LCSFIT is being called from within
C        a loop by setting MEMORY if possible. The out-of-range checking
C        relies on the fact that RIGHT = LEFT + 1 after the last return.
C        Cater to the more likely case of XE in the previous, interior
C        interval.

         MEMORY = XBYARROW >= X (LEFT)  * ARROW .AND.
     &            XBYARROW <  X (RIGHT) * ARROW

         IF (.NOT. MEMORY) THEN
            MEMORY =
     &         LEFT  == 1     .AND. XBYARROW <  X (RIGHT) * ARROW
     &         .OR.
     &         RIGHT == NDATA .AND. XBYARROW >= X (LEFT)  * ARROW
         END IF

      END IF

      IF (MEMORY) GO TO 70 ! Skip the bulk of the computation

C     Loop over evaluation points requiring a new search:
C     ---------------------------------------------------

   10 CONTINUE

C        Interpolation search for bracketing interval:
C        ---------------------------------------------

         CALL INTERVAL (NDATA, X, XE, ARROW, LEFT)

         RIGHT = LEFT + 1

C         -------------------------------------------
C        |                                           |
C        |   1 <= LEFT < RIGHT = LEFT + 1 <= NDATA   |
C        |                                           |
C         -------------------------------------------

C        Compute derivatives by finite-differences:
C        ------------------------------------------

         IF (NDATA > 2 .AND. .NOT. LINEAR) THEN

C           Interval and derivative approximations:
C           ---------------------------------------

C           The following duplicates more code than PLSFIT's approach,
C           but eliminates some indirection - no need to wrap-around here.
C           Handle the end conditions first to minimize testing LEFT, RIGHT.

            IF (LEFT == 1) THEN

               H (0) = X (2) - X (1)
               DELY (0) = (Y (2) - Y (1)) / H (0)
               H (1) = X (3) - X (2)
               DELY (1) = (Y (3) - Y (2)) / H (1)

               IF (CYCLIC) THEN ! Loose fit assumed
                  H (-1) = X (NDATA) - X (NDATA - 1)
                  DELY (-1) = (Y (NDATA) - Y (NDATA - 1)) / H (-1)
                  B (0) = BESSEL (0, H, DELY)
                  B (1) = BESSEL (1, H, DELY)
               ELSE
                  IF (MONO) THEN
                     B (0) = BUTLAND (0, H, DELY)
                     B (1) = BRODLIE (1, H, DELY)
                  ELSE
                     B (0) = THREEPT (0, H, DELY)
                     B (1) = BESSEL  (1, H, DELY)
                  END IF
               END IF

            ELSE IF (RIGHT == NDATA) THEN

               H(-1) = X (LEFT) - X (LEFT-1)
               DELY(-1) = (Y (LEFT) - Y (LEFT-1)) / H (-1)
               H (0) = X (RIGHT) - X (LEFT)
               DELY (0) = (Y (RIGHT) - Y (LEFT))  / H (0)

               IF (CYCLIC) THEN
                  H (1) = X (2) - X (1)
                  DELY (1) = (Y (2) - Y (1)) / H (1)
                  B (0) = BESSEL (0, H, DELY)
                  B (1) = BESSEL (1, H, DELY)
               ELSE

                  IF (MONO) THEN
                     B (0) = BRODLIE (0, H, DELY)
                     B (1) = BUTLAND (1, H, DELY)
                  ELSE
                     B (0) = BESSEL  (0, H, DELY)
                     B (1) = THREEPT (1, H, DELY)
                  END IF
               END IF

            ELSE

               K = LEFT
               DO J = -1, +1
                  H (J)    =  X (K) - X (K-1)
                  DELY (J) = (Y (K) - Y (K-1)) / H (J)
                  K = K + 1
               END DO

C              Select interpolation scheme:
C              ----------------------------

C              Compute (possibly adjusted) first derivatives at both
C              left- and right-hand endpoints of the interval.

               IF (MONO) THEN

C                 Monotone - use Brodlie modification of Butland's
C                 formula to adjust the derivatives at the knots.

                  B (0) = BRODLIE (0, H, DELY)
                  B (1) = BRODLIE (1, H, DELY)

               ELSE ! METHOD = 'B'

C                 Bessel - use central difference formula at the knots.

                  B (0) = BESSEL (0, H, DELY)
                  B (1) = BESSEL (1, H, DELY)

               END IF

            END IF

C           Compute the remaining cubic coefficients.

            C = (THREE * DELY (0) - TWO * B (0) - B (1)) / H (0)
            D = ( -TWO * DELY (0) + B (0) + B (1)) / H (0) ** 2

         ELSE ! NDATA = 2 .OR. METHOD = 'L'

C           Degenerate case (linear).
C           -------------------------

            B (0) = (Y (RIGHT) - Y (LEFT)) / (X (RIGHT) - X (LEFT))
            C     = ZERO
            D     = ZERO

         END IF

C        Evaluate the cubic (derivative first in case only YEVAL is reqd.):
C        ------------------------------------------------------------------

   70    CONTINUE ! Start of same-interval loop inside new-interval loop

            DX = XE - X (LEFT)
            YPEVAL (IEVAL) = B (0) + DX * (TWO * C + DX * THREE * D)
            YEVAL  (IEVAL) = Y (LEFT) + DX * (B (0) + DX * (C + DX * D))

C           The next evaluation point may be in the same interval:
C           ------------------------------------------------------

            IF (IEVAL < NEVAL) THEN ! Skips this if NEVAL = 1

               IEVAL = IEVAL + 1
               XE = XEVAL (IEVAL)
               XBYARROW = XE * ARROW
               IF (XBYARROW >= X (LEFT)  * ARROW  .AND.
     &             XBYARROW <  X (RIGHT) * ARROW) GO TO 70
            
               GO TO 10 ! Else much more work required.

            END IF

C     Termination:
C     ------------

      RETURN

      END SUBROUTINE LCSFIT
C+----------------------------------------------------------------------
C
      SUBROUTINE INTERVAL (NX, X, XFIND, ARROW, LEFT)
C
C     One-liner: Interpolation search for interval containing a point.
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        Written primarily for interval-based interpolations such as
C     piecewise linear or cubic spline, INTERVAL performs a search to
C     locate the best interval for evaluating the interpolant at a
C     given point. The normal case returns the "left-hand" endpoint of
C     the interval bracketing the point, but for the out-of-range cases
C     below or above the range of the knots, the interval to be used is
C     the first or last. The array of knots must be monotonic, either
C     increasing or decreasing. Diagrammatically, LEFT is returned as
C     shown below for the normal case (no extrapolation):
C
C          X (1)  ...   X (LEFT)   X (LEFT+1)   ...      X (NX)
C                                ^
C                              XFIND
C
C     And for extrapolation:
C
C                     X (LEFT = 1)  ...   X (NX)
C             ^
C           XFIND
C
C     or,
C                X (1)  ...   X (LEFT = NX-1)    X (NX)
C                                                           ^
C                                                         XFIND
C
C     If the point to be bracketed (XFIND) matches one of the knots, the
C     index of that knot is returned as LEFT, i.e., the condition for a
C     bracket of an interior point is:
C
C        X (LEFT) <= XFIND < X (LEFT+1)  if  ARROW = +1.0,  or
C        X (LEFT) >= XFIND > X (LEFT+1)  if  ARROW = -1.0.
C
C        This is a low-level routine with minimal error checking. The
C     calling program is assumed to have verified the following:
C
C     (1)  NX >= 2
C     (2)  X strictly monotonic
C     (3)  ARROW = +1.0 or -1.0
C
C     Subroutine PROTECT is available from the author for easily checking
C     conditions (2) and (3). LEFT is verified on input, but efficiency in
C     loops will benefit from passing the best estimate available, usually
C     just the result of the last call.
C
C        INTERVAL was originally written for use with CSEVAL and TABLE1.
C     The interpolation search was adapted from ideas in Sedgewick's book
C     referenced below.
C
C     Arguments:
C     ----------
C
C     Name  Dimension  Type  I/O/S  Description
C     NX                I    I      Number of points in array X; must
C                                   be >= 2 (no check performed).
C
C     X        NX       R    I      Array of points defining the set
C                                   of intervals to be examined. Only
C                                   the first NX-1 points are required.
C
C     XFIND             R    I      The point for which a bracketing
C                                   interval is sought.
C
C     ARROW             R    I      Monotonicity indicator for input
C                                   array X:
C                                     -1.0  strictly decreasing
C                                      0.0  NOT ALLOWED!
C                                     +1.0  strictly increasing
C                                   Supplied by the calling routine for
C                                   reasons of speed (not checked).
C
C     LEFT              I    I/O    Input: guessed index of left-hand
C                                   endpoint of the interval containing
C                                   the specified point.
C
C                                   Output: index of the largest array
C                                   value <= specified point (if ARROW=+1.0).
C                                   Special case for data out of range:
C                                   return left endpoint of closest interval.
C                                   Thus, LEFT = 1 for XFIND < X (2), and
C                                   LEFT = NX-1 for XFIND >= X (NX-1).
C                                   (If ARROW=-1.0, reverse the inequalities.)
C
C     Environment:  Digital VAX-11/780, VMS FORTRAN
C     ------------  Apple Macintosh, Absoft MacFORTRAN/020 v2.3
C
C     Notes:
C     ------
C
C     (1)  IMPLICIT NONE and eight character symbols are not (yet) standard.
C
C     (2)  In speed-critical applications, it might be a good idea to build
C          this algorithm in-line since it is typically called many times
C          from within a loop. Another potential speed-up is removal of the
C          ARROW multiplies, which restricts the method to increasing data.
C          So far, the simplicity of separating out the messy search details
C          and the generality of bi-directional searching have outweighed
C          the modest speed penalty incurred.
C
C     Bibliography:
C     -------------
C
C     (1) Sedgewick, R.  Algorithms.  Reading: Addison-Wesley, 1983.
C            (Chap. 14)
C
C     Author:  Robert Kennelly and David Saunders, Sterling Federal Systems
C     -------
C
C     Development history:
C     --------------------
C
C     20 Oct. 1987    RAK    Interpolation search adapted (with mods.
C                            for bidirectional search and some minor
C                            repair) from CSEVAL (RAK) and TABLE1 (DAS).
C     08 Aug. 1988    DAS    Clarified descriptions of bracketing, where
C                            the inequalities depend upon ARROW.
C
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     &   ONE
      PARAMETER
     &  (ONE = 1.0E+0)

C     Arguments.

      INTEGER
     &   LEFT, NX
      REAL
     &   ARROW, X (NX), XFIND

C     Local variables.

      INTEGER
     &   LENGTH, NXLESS1, RIGHT, TRIAL
      REAL
     &   XBYARROW

C     Execution.
C     ----------

      XBYARROW = XFIND * ARROW

C     Simplify things by disposing of two important special cases so that
C     X (LEFT) and X (RIGHT) can really bracket XFIND. As a by-product,
C     this also takes care of the NX = 2, 3 cases.

      NXLESS1 = NX - 1

      IF (XBYARROW .GE. X (NXLESS1) * ARROW) THEN
         LEFT = NXLESS1
         GO TO 990
      ELSE IF (XBYARROW .LT. X (2) * ARROW) THEN
         LEFT = 1
         GO TO 990
      END IF

C      ---------------------------------
C     |                                 |
C     |   X (2) <= XFIND < X (NX - 1)   |
C     |            - or -               |
C     |   X (2) > XFIND >= X (NX - 1)   |
C     |                                 |
C     |   NX > 3                        |
C     |                                 |
C      ---------------------------------

C     Adjust the pointers. We hope that the calling routine has provided
C     a reasonable guess (since it's probably working on an ordered array
C     of points to evaluate), but check anyway.

      LEFT = MIN (MAX (2, LEFT), NX - 2)

      IF (XBYARROW .GE. X (LEFT) * ARROW) THEN
         IF (XBYARROW .LT. X (LEFT + 1) * ARROW) THEN

C           XFIND is in the original guessed-at interval.

            GO TO 990
         ELSE

C           We'll look farther to the right. Evidently LEFT was < NX - 2.

            RIGHT = NXLESS1
            LEFT  = LEFT + 1
         END IF
      ELSE

C        Look to the left of the guess. Evidently LEFT was > 2.

         RIGHT = LEFT
         LEFT  = 2
      END IF

C      ----------------------------------
C     |                                  |
C     |   2 <= LEFT < RIGHT <= NX - 1    |
C     |                                  |
C      ----------------------------------

C     The interval length must decrease each time through - terminate
C     when the correct interval is found or when the interval length
C     cannot be decreased.

   10 CONTINUE
         LENGTH = RIGHT - LEFT
         IF (LENGTH .GT. 1) THEN

C           The trial value is a "linear" estimate of the left-hand endpoint
C           of the interval bracketing the target XFIND, with protection
C           against round-off (which can affect convergence).

            TRIAL = MIN (RIGHT - 1, LEFT + MAX (0, INT (REAL (LENGTH) *
     &         (XFIND - X (LEFT)) / (X (RIGHT) - X (LEFT)))))

C            ------------------------------------------
C           |                                          |
C           |   2 <= LEFT <= TRIAL < RIGHT <= NX - 1   |
C           |                                          |
C            ------------------------------------------

C           Adjust pointers. Increase LEFT or decrease RIGHT until done.

            IF (XBYARROW .GE. X (TRIAL + 1) * ARROW) THEN
               LEFT  = TRIAL + 1
            ELSE IF (XBYARROW .LT. X (TRIAL) * ARROW) THEN
               RIGHT = TRIAL
            ELSE

C              We're done: XFIND is in the interval [X (TRIAL), X (TRIAL+1)).

               LEFT  = TRIAL
               GO TO 990
            END IF
            GO TO 10

         END IF

C     Termination.
C     ------------

  990 CONTINUE
      RETURN
      END
C+----------------------------------------------------------------------
C
      FUNCTION BESSEL (J, H, DEL)
C
C     One-liner: First derivative using central 3-point formula
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        Computes a first derivative approximation using the central
C     3-point formula.  The data must be in the form of arrays containing
C     finite difference interval lengths and 2-point forward difference
C     derivatives.  BESSEL is intended to be used by PLSFIT for determin-
C     ing end conditions on an interval for (non-monotonic) interpolation
C     by piecewise cubics.  See the PLSFIT header for more details.
C
C     Arguments:
C     ----------
C
C     Name    Type/Dimension  I/O/S  Description
C     J       I               I      Indicates at which end of the
C                                    interval the derivative is to be
C                                    estimated. J = 0 means left-hand
C                                    side, J = 1 means right.
C
C     H       R (-1:1)        I      Array of interval lengths. The 0th
C                                    element is the length of the interval
C                                    on which the cubic is to be deter-
C                                    mined.
C
C     DEL     R (-1:1)        I      Array of derivative estimates. The
C                                    0th element is the forward difference
C                                    derivative over the interval on which
C                                    the cubic is to be determined.
C                                     
C     BESSEL  R                 O    The function value is the adjusted
C                                    derivative.
C
C     Notes:
C     ------
C
C     (1)  IMPLICIT NONE is non-standard.
C
C     Author:  Robert Kennelly, Sterling Federal Systems/NASA-Ames
C     -------
C
C     Development history:
C     --------------------
C
C     18 Feb. 1987    RAK    Initial design and coding.
C
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     &   ONE
      PARAMETER
     &  (ONE = 1.0E+0)

C     Arguments.

      INTEGER
     &   J
      REAL
     &   H (-1:1), DEL (-1:1), BESSEL

C     Local variables.

      REAL
     &   WEIGHT

C     Execution.
C     ----------

C     Estimate first derivative on left (J = 0) or right side (J = 1) of
C     an interval.

      WEIGHT = H (J) / (H (J) + H (J - 1))
      BESSEL = WEIGHT * DEL (J - 1) + (ONE - WEIGHT) * DEL (J)

C     Termination.
C     ------------

      RETURN
      END
C+----------------------------------------------------------------------
C
      FUNCTION BRODLIE (J, H, DEL)
C
C     One-liner: First derivative, adjusted for monotonicity
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        BRODLIE is intended to be used by PLSFIT for determining end
C     conditions on an interval for monotonic interpolation by piecewise
C     cubics. The data must be in the form of arrays containing finite
C     difference interval lengths and 2-point forward difference deriva-
C     tives. See the PLSFIT header for more details.
C
C        The method is due to Brodlie, Butland, Carlson, and Fritsch,
C     as referenced in the PLSFIT header.
C
C     Arguments:
C     ----------
C
C     Name    Type/Dimension  I/O/S  Description
C     J       I               I      Indicates at which end of the
C                                    interval the derivative is to be
C                                    estimated. J = 0 means left-hand
C                                    side, J = 1 means right.
C
C     H       R (-1:1)        I      Array of interval lengths. The 0th
C                                    element is the length of the interval
C                                    on which the cubic is to be deter-
C                                    mined.
C
C     DEL     R (-1:1)        I      Array of derivative estimates. The
C                                    0th element is the forward difference
C                                    derivative over the interval on which
C                                    the cubic is to be determined.
C
C     BRODLIE R                 O    The function value is the adjusted
C                                    derivative.
C
C     Notes:
C     ------
C
C     (1)  IMPLICIT NONE and 8-character symbolic names are non-standard.
C
C     Author:  Robert Kennelly, Sterling Federal Systems/NASA-Ames
C     -------
C
C     Development history:
C     --------------------
C
C     18 Feb. 1987    RAK    Initial design and coding.
C     04 Dec. 2002    DAS    SIGN work-around for Intel compiler bug.
C
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     &   ZERO, ONE, THIRD
      PARAMETER
     &  (ZERO   = 0.0E+0,
     &   ONE    = 1.0E+0,
     &   THIRD  = ONE / 3.0E+0)

C     Arguments.

      INTEGER
     &   J
      REAL
     &   BRODLIE, H (-1:1), DEL (-1:1)

C     Local variables.

      REAL
     &   ALPHA, PRODUCT

C     Execution.
C     ----------

C     Compare the algebraic signs of the two DEL's.  Have to test that
C     at least one is positive to avoid a zero denominator (this fancy
C     test permits one term to be zero, but the answer below is zero
C     anyway in these cases).  The trick is to work around the SIGN
C     function, which returns positive even if its 2nd argument is zero.

C**** NO:  SIGN misbehaves on Intel systems when the 2nd argument is zero.

      PRODUCT = DEL (J - 1) * DEL (J)

      IF (PRODUCT == ZERO) THEN

         BRODLIE = ZERO

      ELSE IF (SIGN (ONE, -DEL (J - 1)) .NE. SIGN (ONE, DEL (J))) THEN

C        Form "weighted harmonic mean" of the two finite-difference
C        derivative approximations.  Note that we try to avoid overflow
C        by not multiplying them together directly.

         ALPHA   = THIRD * (ONE + H (J) / (H (J - 1) + H (J)))
         BRODLIE = PRODUCT / (ALPHA * DEL (J) +
     &      (ONE - ALPHA) * DEL (J - 1))
      ELSE

C        The signs differ, so make this point a local extremum.

         BRODLIE = ZERO
      END IF

C     Termination.
C     ------------

      RETURN
      END
C+----------------------------------------------------------------------
C
      FUNCTION BUTLAND (J, H, DEL)
C
C     One-liner: First derivative, non-central 3-point formula, adjusted
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        Computes a first derivative approximation for PLSFIT over an
C     interval at a data boundary, using a modified forward or backward
C     3-point formula.  The data must be in the form of arrays containing
C     finite difference interval lengths and 2-point forward difference
C     derivatives, and the differencing direction is controlled by a flag.
C     See PLSFIT for more details, or THREEPT for the pure 3-pt. formula.
C
C        The "shape preserving adjustments" are from PCHIP, a monotone
C     piecewise cubic interpolation package by F. N. Fritsch.
C
C     Arguments:
C     ----------
C
C     Name    Type/Dimension  I/O/S  Description
C     J       I               I      Indicates at which end of the
C                                    interval the derivative is to be
C                                    estimated. J = 0 means left-hand
C                                    side, J = 1 means right.
C
C     H       R (-1:1)        I      Array of interval lengths. The 0th
C                                    element is the length of the interval
C                                    on which the cubic is to be deter-
C                                    mined.
C
C     DEL     R (-1:1)        I      Array of derivative estimates. The
C                                    0th element is the forward difference
C                                    derivative over the interval on which
C                                    the cubic is to be determined.
C
C     BUTLAND R                 O    The function value is the adjusted
C                                    derivative.
C
C     Environment:  Fortran 90
C     ------------
C
C     Author:  Robert Kennelly, Sterling Federal Systems/NASA-Ames
C     -------
C
C     History:
C     --------
C
C     18 Feb. 1987    RAK    Initial design and coding, as THREEPT.
C     20 June 1991    DAS    Monotonic form renamed BUTLAND; THREEPT
C                            is now the pure 3-point formula.
C     04 Dec. 2002     "     SIGN work-arounds for Intel compiler bug.
C
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     &   J
      REAL
     &   H (-1:1), DEL (-1:1), BUTLAND

C     Local constants.

      REAL, PARAMETER ::
     &   ZERO  = 0.0E+0,
     &   ONE   = 1.0E+0,
     &   THREE = 3.0E+0

C     Local variables.

      INTEGER
     &   STEP
      REAL
     &   DMAX, WEIGHT
      LOGICAL
     &   CONSTRAIN

C     Execution.
C     ----------

C     Estimate first derivative on a left-hand boundary using a 3-point
C     forward difference (STEP = +1), or with a backward difference for
C     the right-hand boundary (STEP = -1).

      STEP = 1 - J - J   ! J here is consistent with related modules.

C     In {H, DEL} form, the derivative looks like a weighted average.

C***  Avoid zero as the second argument of SIGN: Intel compiler misbehaves.

      IF (DEL (0) == ZERO) THEN

         BUTLAND = ZERO

      ELSE ! BUTLAND below cannot be zero

         WEIGHT  = -H (0) / (H (0) + H (STEP))
         BUTLAND = WEIGHT * DEL (STEP) + (ONE - WEIGHT) * DEL (0)

C        Shape-preserving adjustments.  Note that we try to avoid overflow
C        by not multiplying quantities directly.

         IF (SIGN (ONE, BUTLAND) /= SIGN (ONE, DEL (0))) THEN

C           Defer to the estimate closest to the boundary.

            BUTLAND = ZERO

C******  ELSE IF (SIGN (ONE, DEL (0)) .NE. SIGN (ONE, DEL (STEP))) THEN
         ELSE

            IF (DEL (STEP) == ZERO) THEN
               CONSTRAIN = DEL (0) < ZERO
            ELSE
               CONSTRAIN = SIGN (ONE, DEL (0)) /= SIGN (ONE, DEL (STEP))
            END IF

            IF (CONSTRAIN) THEN

C              If the monotonicity switches, may need to bound the estimate.

               DMAX = THREE * DEL (0)
               IF (ABS (BUTLAND) > ABS (DMAX)) BUTLAND = DMAX
            END IF

         END IF

      END IF

C     Termination.
C     ------------

      END
C+----------------------------------------------------------------------
C
      FUNCTION THREEPT (J, H, DEL)
C
C     One-liner: First derivative, non-central 3-point formula
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        Computes a first derivative approximation for PLSFIT over an
C     interval at a data boundary, using a forward or backward 3-point
C     formula.  The data must be in the form of arrays containing finite
C     difference interval lengths and 2-point forward difference deriva-
C     tives, and the differencing direction is controlled by a flag. See
C     PLSFIT for more details.
C
C        See module BUTLAND for a version with "shape-preserving"
C     adjustments.
C
C     Arguments:
C     ----------
C
C     Name    Type/Dimension  I/O/S  Description
C     J       I               I      Indicates at which end of the
C                                    interval the derivative is to be
C                                    estimated. J = 0 means left-hand
C                                    side, J = 1 means right. 
C
C     H       R (-1:1)        I      Array of interval lengths. The 0th
C                                    element is the length of the interval
C                                    on which the cubic is to be deter-
C                                    mined.
C
C     DEL     R (-1:1)        I      Array of derivative estimates. The
C                                    0th element is the forward difference
C                                    derivative over the interval on which
C                                    the cubic is to be determined.
C                                     
C     THREEPT R                 O    The function value is the derivative.
C
C     Environment:  VAX/VMS; FORTRAN 77
C     ------------
C
C     IMPLICIT NONE and 8-character symbolic names are non-standard.
C
C     Author:  Robert Kennelly, Sterling Federal Systems/NASA-Ames
C     -------
C
C     History:
C     --------
C
C     18 Feb. 1987    RAK    Initial design and coding.
C     06 June 1991    DAS    Original THREEPT renamed BUTLAND; THREEPT
C                            now gives unmodified 1-sided 3-pt. results.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     &   J
      REAL
     &   H (-1:1), DEL (-1:1), THREEPT

C     Local constants.

      REAL
     &   ONE
      PARAMETER
     &  (ONE = 1.0E+0)

C     Local variables.

      INTEGER
     &   STEP
      REAL
     &   WEIGHT

C     Execution.
C     ----------

C     Estimate first derivative on a left-hand boundary using a 3-point
C     forward difference (STEP = +1), or with a backward difference for
C     the right-hand boundary (STEP = -1).

      STEP = 1 - J - J   ! J here is consistent with related modules.

C     In {H, DEL} form, the derivative looks like a weighted average.

      WEIGHT  = -H (0) / (H (0) + H (STEP))
      THREEPT = WEIGHT * DEL (STEP) + (ONE - WEIGHT) * DEL (0)

C     Termination.
C     ------------

      RETURN
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
