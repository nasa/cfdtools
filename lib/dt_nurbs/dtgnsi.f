      SUBROUTINE DTGNSI(X,NDIM,NDOM,NPT,Y,NYDIM,NDEP,NDEG,
     &                  WORK,NWORK,MAXC,C,NC,IER)
C
C=======================================================================
C
C  PURPOSE:
C          DTGNSI is used to produce multivariate natural spline 
C          interpolants to data on a regular grid.
C
C  USAGE:
C          DOUBLE PRECISION X(NDIM,NDOM),Y(NPT(1),...,NPT(NDOM),NDEP),
C                           WORK(NWORK),C(MAXC)
C          CALL DTGNSI(X,NDIM,NDOM,NPT,Y,NYDIM,NDEP,NDEG,
C                      WORK,NWORK,MAXC,C,NC,IER)
C
C  INPUT:
C          X       A real array of dimension X(NDIM,*) which contains
C                  the values of the independent variables at the grid
C                  points.
C
C          NDIM    The leading declared dimension of X.
C
C          NDOM    The number of independent variables.
C
C          NPT     An integer array containing the number of data points
C                  in each of the NDOM independent variables.
C
C          Y       The array of the dependent function evaluated on the
C                  grid.
C
C          NYDIM   Array of dimensional constants for Y.
C
C          NDEP    The number of dependent functions.
C
C          NDEG    Array containing the degree of the interpolant to
C                  be used in each independent variable.
C
C          MAXC    Maximum storage allocated to C array.
C
C  STORAGE:
C          WORK    The work array.
C
C          NWORK   Size of the work array.
C                  NWORK .GE.  5 + 3*MO + 3*MP + 3*MO*MP + 5*MO**2
C                    where  MO = 1 + MAX(NDEG(I))
C                    and    MP = MAX(NPT(I))
C
C  OUTPUT:
C          C       The spline definition array for the result.
C
C          NC      Length of spline array.
C
C          IER     Success/error code.
C
C                  IER = 0    Success: results computed.
C
C                  IER = 1    Results are computed, but they are
C                             sensitive to inaccuracies in the data.
C
C                  IER = 2    Results are computed, but they are
C                             strongly sensitive to inaccuracies in
C                             the data.
C
C                  IER = -1   NDEG(I).LT.0.
C
C                  IER = -2   NPT(I) .LT. MAX (2,(NDEG(I)+1)/2)
C
C                  IER = -3   NWORK too small. The number of elements
C                             needed is given by the printed message.
C
C                  IER = -4   NYDIM(I) .LT. NPT(I) for I = 1,...,NDOM.
C
C                  IER = -5   The X values are not strictly increasing.
C
C                  IER = -9   NDOM .LT. 1.
C
C                  IER = -14  The user has not supplied enough storage
C                             for the spline vector, C.
C
C                  IER = -20  The coefficient matrix for the interpolation
C                             problem is singular.
C
C
C  OTHER DT ROUTINES CALLED
C
C        DTERR
C        DTGNS1
C
C     ******************************************************************
C
C  Parameters:
C
      CHARACTER*8 SUBNAM
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.D0,SUBNAM = 'DTGNSI  ')
C
C  Arguments:
C
      INTEGER  NPT(*), NYDIM(*), NDEG(*), NDIM, NDEP, NDOM,
     &         NWORK, IER
      DOUBLE PRECISION X(NDIM,*), Y(*), WORK(NWORK), C(MAXC)
C
C  Internal:
C
      INTEGER  IX,IV,MODE,NEED,MAXNPT,MAXORD
      INTEGER  NCOEF, I, KORD, KHALF, NCI
C
C     ******************************************************************
C
      C(1) = -ONE
      MODE = 1
      NEED = 0
C
C=======================================================================
C    Check NDOM.
C=======================================================================
C
      IF (NDOM .LT. 1) THEN
        IER = -9
        CALL DTERR(MODE,SUBNAM,IER,NEED)
        RETURN
      ENDIF
C
      MAXNPT = 0
      MAXORD = 0
      DO 20 IV = 1,NDOM
        MAXNPT = MAX(MAXNPT,NPT(IV))
C
C=======================================================================
C    Check NDEG.
C=======================================================================
C
        IF (NDEG(IV).LT.0) THEN
          IER = -1
          CALL DTERR(MODE,SUBNAM,IER,NEED)
          RETURN
        ENDIF
        MAXORD = MAX(MAXORD,1+NDEG(IV))
C
C=======================================================================
C    Check NPT().
C=======================================================================
C
        IF ( NPT(IV) .LT. MAX(2,(NDEG(IV)+1)/2) ) THEN
          IER = -2
          CALL DTERR(MODE,SUBNAM,IER,NEED)
          RETURN
        ENDIF
C
C=======================================================================
C    Check NYDIM.
C=======================================================================
C
        IF (NYDIM(IV) .LT. NPT(IV)) THEN
          IER = -4
          CALL DTERR(MODE,SUBNAM,IER,NEED)
          RETURN
        ENDIF
C
        DO 10 IX = 2,NPT(IV)
          IF (X(IX,IV) .LE. X(IX-1,IV)) THEN
            IER = -5
            CALL DTERR(MODE,SUBNAM,IER,NEED)
            RETURN
          ENDIF
   10   CONTINUE
   20 CONTINUE
C
C=======================================================================
C    Check NWORK.
C=======================================================================
C
      NEED = MAXNPT + (5 + 2*MAXNPT + 3*MAXORD)
     &       + MAXORD*(3*MAXNPT + 5*MAXORD)
      IF (NWORK.LT.NEED) THEN
        IER = -3
        MODE = 2
        CALL DTERR(MODE,SUBNAM,IER,NEED)
        RETURN
      ENDIF
C
C=======================================================================
C    Check MAXC.
C=======================================================================
C
      NC = 2 + 3*NDOM
      NCOEF = NDEP
      DO 50 I=1,NDOM
        KORD = 1 + NDEG(I)
        KHALF = (1 + KORD)/2
        NCI = NPT(I) - 2 + 2*KHALF
        NC = NC + (NCI+KORD)
        NCOEF = NCOEF*NCI
   50 CONTINUE
      NC = NC + NCOEF
      IF (MAXC.LT.NC) THEN
        IER = -14
        MODE = 2
        CALL DTERR(MODE,SUBNAM,IER,NEED)
        RETURN
      ENDIF
C
C=======================================================================
C    Call the appropriate fitting routine.
C=======================================================================
C
       CALL DTGNS1(NPT,X,NDIM,Y,NYDIM,NDEG,NDOM,NDEP,
     &                WORK,NWORK,C,IER)
      IF (IER.EQ.0)
     &  THEN
          C(1) = NDOM
        ELSEIF (IER.GT.0) THEN
          CALL DTERR(0,SUBNAM,IER,0)
        ELSE
          CALL DTERR(5,SUBNAM,IER,0)
      ENDIF
C
      RETURN
      END
