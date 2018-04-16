      SUBROUTINE DTGTP1(NPT,X,NDIM,Y,NYDIM,NDEG,NDOM,NDEP,IWT,WHT,
     &                  NDIMK,WORK,NWORK,C,IER)
C
C=======================================================================
C
C  PURPOSE:
C
C----- This is the least squares fitting routine for tensor product
C      splines on a grid of data.
C
C     ******************************************************************
C
C  Arguments:
C
      INTEGER  NPT(*), NDIM, NYDIM(*), NDEG(*), NDOM, NDEP,
     &         IWT, NDIMK, NWORK, IER
      DOUBLE PRECISION  X(NDIM,*), Y(*), WHT(NDIM,*),
     &                  WORK(*), C(*)
C
C  Internal:
C
      INTEGER IS, IV, IX, IZ, NSYS, ICC, IPREAM
      INTEGER IW1SIZ, IW2, IW2SIZ, IW3, IW3SIZ, KNOTST,
     &        IFAIL, NKNOTS, IERMAX
C
C     ******************************************************************
C
C----- Compute pointers and work space sizes.
C
      NSYS = NDEP
      IPREAM = 2 + 3 * NDOM
      IW1SIZ = 0
      IW2SIZ = 0
      IW3SIZ = 0
      IERMAX = 0
      DO 110 IV = 1,NDOM
        NSYS = NSYS * NPT(IV)
        KORD = C(2+IV)
        NCOEF = C(2+NDOM+IV)
        IPREAM = IPREAM + KORD + NCOEF
        IW1SIZ = MAX(IW1SIZ, NPT(IV))
        IW2SIZ = MAX(IW2SIZ,5+2*NCOEF+KORD)
        IW3SIZ = MAX(IW3SIZ,KORD*(NCOEF+3))
  110 CONTINUE
      IW2 = IW1SIZ + 1
      IW3 = IW2 + IW2SIZ
      IW4 = IW3 + IW3SIZ
      IWPSIZ = NWORK - IW1SIZ - IW2SIZ - IW3SIZ
C
C----- Copy and compress the data into the work vector.
C
      IWPT = IW4 - 1
      IYPT = 0
      NLOOP = NSYS/NPT(1)
      DO 140 J = 1,NLOOP
C
        DO 120 I=1,NPT(1)
          WORK(IWPT+I) = Y(IYPT+I)
  120   CONTINUE
C
        IWPT = IWPT + NPT(1)
        IYPT = IYPT + NYDIM(1)
C
        JREDUC = J
        IPRODY = NYDIM(1)
        DO 130 IX = 2,NDOM
          IF (MOD(JREDUC,NPT(IX)) .NE. 0) GO TO 140
          JREDUC = JREDUC/NPT(IX)
          IYPT = IYPT + IPRODY*(NYDIM(IX) - NPT(IX))
          IPRODY = IPRODY*NYDIM(IX)
  130   CONTINUE
C
  140 CONTINUE
C
C----- Loop through each variable.
C
      ISKPOL = 1
      KNOTST = 3 + 3 * NDOM
      DO 220 IV = 1,NDOM
        ICC = -1
        NSYS = NSYS/NPT(IV)
        NCOEF = NINT( C(2+NDOM+IV) )
        KORD = NINT( C(2+IV) )
        NKNOTS = NCOEF + KORD
        IKNOTS = NKNOTS - 2*KORD
        IKNST = KNOTST + KORD
        IZ = IW2 + 5 + NKNOTS
        IC = IW4
        IW = IC
        ISKP1 = NPT(IV)*ISKPOL
        ISKP2 = NCOEF*ISKPOL
C
C----- Solve each of the univariate systems for this pass.
C
        NLOOP = NSYS/ISKPOL
        DO 210 IS = 1,NLOOP
          DO 205 JS = 0,ISKPOL-1
            CALL DCOPY(NPT(IV),WORK(IC+JS),ISKPOL,WORK,1)
            CALL DTLSAA(NPT(IV),X(1,IV),WORK,NDEG(IV),ICC,IWT,
     &                  WHT(1,IV),C(IKNST),IKNOTS,WORK(IW3),IW3SIZ,
     &                  IW2SIZ,WORK(IW2),NCD,IFAIL,IER)
            IF (IER .LT. 0) THEN
              IER = -200
              RETURN
            ENDIF
            IERMAX = MAX (IER, IERMAX)
            CALL DCOPY(NCOEF,WORK(IZ),1,WORK(IW+JS),ISKPOL)
  205     CONTINUE
          IC = IC + ISKP1
          IW = IW + ISKP2
  210   CONTINUE
        NSYS = NSYS * NCOEF
        ISKPOL = ISKP2
        KNOTST = KNOTST + NKNOTS
  220 CONTINUE
      IER = IERMAX
C
C             Place coefficients into C vector.
C
      IC = IPREAM
      IW = IW4 - 1
      DO 340 I = 1,NSYS
        C(IC+I) = WORK(IW+I)
  340 CONTINUE
C
      RETURN
      END
