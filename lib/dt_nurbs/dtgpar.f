      SUBROUTINE DTGPAR(NPT,Y,NDOM,NDEP,WORK,NWORK,NDIM,X,IER)
C
C=======================================================================
C
C  PURPOSE:
C          Construct a sequence of parameter values using
C          normalized chord length.
C
C  USAGE:
C          DOUBLE PRECISION  Y(NPT(1),..,NPT(NDOM),NDEP),
C                            X(NDIM,NDOM),WORK(NWORK)
C          INTEGER  NPT(NDOM)
C          CALL DTGPAR(NPT,Y,NDOM,NDEP,WORK,NWORK,NDIM,X,IER)
C
C  INPUTS:
C          NPT     Array containing the number of data points
C                  in each dimension.
C
C          Y       Array of data points with each column
C                  corresponding to coordinate axis.
C
C          NDOM    Number of independent variables.
C
C          NDEP    Number of dependent variables.
C
C          NDIM    Row dimension of X-array.  Must be at least as
C                  big as the largest number of data points in
C                  any dimension.
C
C  STORAGE:
C          WORK    Work array of length NWORK.
C
C          NWORK   Length of work array,
C                  NWORK .GE. MAX(NPT(I))
C
C  OUTPUT:
C          X       Array of parameter values in each direction.
C
C          IER     Success/error code.
C
C                  IER =  0  Success.
C
C                  IER = -2  NPT(.) .LT. 1.
C
C                  IER = -3  NWORK too small.
C
C                  IER = -4  NDIM .LT. MAX( NPT(.) )
C
C                  IER = -8  NDEP.LT.1
C
C                  IER = -9  NDOM.LT.1
C
C
C     ******************************************************************
C
      CHARACTER*8  SUBNAM
      DOUBLE PRECISION ZERO, ONE
      PARAMETER  (ZERO=0.D0,ONE=1.D0,SUBNAM='DTGPAR  ')
C
C  Arguments:
      INTEGER  NDOM,NPT(NDOM),NDEP,NDIM,IER
      DOUBLE PRECISION  Y(*),X(NDIM,NDOM),WORK(NWORK)
C
C  Internal:
      DOUBLE PRECISION  SUM,WN,RAVG
      INTEGER  I,IS,ISKPOL,ISKP2,IV,J0,J1,J2,JS,MODE,NAVG
      INTEGER  NEED, MAXPTS, NSYS, IY0, NLOOP, IY, K
C
C     ******************************************************************
C
C  Initialize.
C
      IER = 0
      MODE = 1
C
C=======================================================================
C    Check NDOM.
C=======================================================================
C
      IF (NDOM.LT.1) THEN
        IER = -9
        CALL DTERR(MODE,SUBNAM,IER,NEED)
        RETURN
      ENDIF
C
C=======================================================================
C    Check NDEP.
C=======================================================================
C
      IF (NDEP.LT.1) THEN
        IER = -8
        CALL DTERR(MODE,SUBNAM,IER,NEED)
        RETURN
      ENDIF
C
C=======================================================================
C    Check NPT.
C=======================================================================
C
      MAXPTS = 0
      NSYS = 1
      DO 10 I=1,NDOM
        IF (NPT(I).LT.1) THEN
          MODE = 1
          IER = -2
          CALL DTERR(MODE,SUBNAM,IER,NEED)
          RETURN
        ENDIF
        NSYS = NSYS*NPT(I)
        MAXPTS = MAX(MAXPTS,NPT(I))
   10 CONTINUE
C
C=======================================================================
C    Check NDIM.
C=======================================================================
C
      IF (NDIM.LT.MAXPTS) THEN
        MODE = 1
        IER = -4
        CALL DTERR(MODE,SUBNAM,IER,NEED)
        RETURN
      ENDIF
C
C=======================================================================
C    Check NWORK.
C=======================================================================
C
      NEED = MAXPTS
      IF (NWORK.LT.NEED) THEN
        MODE = 2
        IER = -3
        CALL DTERR(MODE,SUBNAM,IER,NEED)
        RETURN
      ENDIF
C
C=======================================================================
C    Compute chord lengths in each dimension.
C=======================================================================
C
      ISKPOL = 1
      DO 400 IV=1,NDOM
C
        DO 200 I=1,NPT(IV)
          X(I,IV) = ZERO
  200   CONTINUE
        IY0 = 1
        ISKP2 = ISKPOL*NPT(IV)
        NLOOP = NSYS/ISKP2
C
        DO 350 IS=1,NLOOP
          DO 340 JS=0,ISKPOL-1
            IY = IY0 + JS
C
C             Compute distances along current line.
C
            WORK(1) = ZERO
            DO 300 I=2,NPT(IV)
              SUM = ZERO
              J0 = IY + (I-2)*ISKPOL
              DO 250 K=1,NDEP
                J1 = J0 + (K-1)*NSYS
                J2 = J1 + ISKPOL
                SUM = SUM + (Y(J2)-Y(J1))**2
  250         CONTINUE
              WORK(I) = WORK(I-1) + SQRT(SUM)
  300       CONTINUE
C
C             Normalize distances by total to reduce interval to [0,1].
C
            WN = WORK(NPT(IV))
            DO 320 I=2,NPT(IV)-1
              X(I,IV) = X(I,IV) + WORK(I)/WN
  320       CONTINUE
C
  340     CONTINUE
          IY0 = IY0 + ISKP2
  350   CONTINUE
C
C             Average parameter lengths in current direction.
C
        X(1,IV) = ZERO
        NAVG = NLOOP*ISKPOL
        RAVG = ONE/FLOAT(NAVG)
        DO 370 I=2,NPT(IV)-1
          X(I,IV) = RAVG*X(I,IV)
  370   CONTINUE
        X(NPT(IV),IV) = ONE
C
        ISKPOL = ISKP2
  400 CONTINUE
C
      RETURN
      END
