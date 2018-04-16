      SUBROUTINE DTGCSI(X,NDIM,NDOM,NPT,Y,NYDIM,YLHS,NLHS,YRHS,NRHS,
     &                  NBDIM,NDEP,NDEG,WORK,NWORK,MAXC,C,NC,IER)
C
C=======================================================================
C
C  PURPOSE:
C          DTGCSI is used to produce multivariate complete spline 
C          interpolants to data on a regular grid.
C
C  USAGE:
C          DOUBLE PRECISION X(NDIM,NDOM),WORK(NWORK),C(MAXC),
C                           Y(NYDIM(1),...,NYDIM(NDOM),NDEP),
C                           YLHS(NBDIM(1),...,NBDIM(3),NDOM),
C                           YRHS(NBDIM(1),...,NBDIM(3),NDOM),
C          INTEGER  NPT(NDOM),NYDIM(NDOM),NLHS(NDOM),NRHS(NDOM),
C                   NBDIM(3),NDEG(NDOM)
C          CALL DTGCSI(X,NDIM,NDOM,NPT,Y,NYDIM,YLHS,NLHS,YRHS,NRHS,
C                      NBDIM,NDEP,NDEG,WORK,NWORK,MAXC,C,NC,IER)
C
C  INPUT:
C          X       A real array which contains the values of the
C                  independent variables at the grid points.
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
C          YLHS    Array containing boundary conditions on left hand
C                  hyperplanes.  YLHS(i,j,k,m) is the ith partial
C                  derivative of dependent function k with respect 
C                  to variable m at the jth point of the hyperplane
C                  perpendicular to the direction m.
C
C          NLHS    Number of partial derivative boundary conditions 
C                  on left hyperplane in each dimension.
C
C          YRHS    Array of boundary conditions on right hand 
C                  hyperplanes in same format as YLHS.
C
C          NRHS    Number of partial derivative boundary conditions 
C                  on right hyperplane in each dimension.
C
C          NBDIM   Array of dimensional constants for YLHS and YRHS.
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
C                  NWORK .GE. 3*(MO+1)*NC + 2MO*MO + MO + 5
C                    where  MO = 1 + MAX(NDEG(I))
C                           NC = MAX(1+NDEG(I)+NPT(I))
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
C                  IER = -1   NDEG.LT.0.
C
C                  IER = -2   NPT(I) .LT. MAX (2,(NDEG+1)/2)
C
C                  IER = -3   NWORK too small. The number of elements
C                             needed is given by the printed message.
C
C                  IER = -4   NYDIM(I) .LT. NPT(I) for I = 1,...,NDOM.
C
C                  IER = -5   The X values are not strictly increasing.
C
C                  IER = -6   NBDIM(1) .LT. MAX(NLHS(I),NRHS(I)).
C
C                  IER = -7   NBDIM(2) is too small to contain all
C                             of the boundary points in some
C                             direction.
C
C                  IER = -8   NBDIM(3) .LT. NDEP.
C
C                  IER = -9   NDOM .LT. 1.
C
C                  IER = -10  MAX (NLHS(I), NRHS(I)) .GT. NDEG(I)/2 
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
C        DTGCS1
C
C     ******************************************************************
C
C  Parameters:
C
      CHARACTER*8 SUBNAM
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.D0,SUBNAM = 'DTGCSI  ')
C
C  Arguments:
C
      INTEGER  NPT(*), NYDIM(*), NDEG(*), NLHS(*), NRHS(*), NDIM,
     &         NBDIM(*), NDEP, NDOM, NWORK, IER
      DOUBLE PRECISION  X(NDIM,*), Y(*), YLHS(*), YRHS(*), 
     &                  WORK(NWORK), C(MAXC)
C
C  Internal:
C
      INTEGER  IX,IV,MODE,NEED,MAXNPT,MAXDEG,MAXLHS,MAXRHS,KORD,KHALF,
     &         NC,NCI,NCOEF,MINNPT,NPTOT,MAXORD,MAXCOF
      INTEGER  I
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
      MINNPT = NPT(1)
      NPTOT  = 1
      MAXDEG = 0
      MAXLHS = 0
      MAXRHS = 0
      DO 20 IV = 1,NDOM
        MAXNPT = MAX(MAXNPT,NPT(IV))
        MINNPT = MIN(MINNPT,NPT(IV))
        MAXLHS = MAX(MAXLHS,NLHS(IV))
        MAXRHS = MAX(MAXRHS,NRHS(IV))
        NPTOT = NPTOT*NPT(IV)
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
        MAXDEG = MAX(MAXDEG,NDEG(IV))
C
C=======================================================================
C    Check too many boundary conditions
C=======================================================================
        IF (MAX(NLHS(IV),NRHS(IV)).GT.NDEG(IV)/2) THEN
          IER = -10
          CALL DTERR(MODE,SUBNAM,IER,NEED)
          RETURN
        ENDIF
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
10      CONTINUE
20    CONTINUE
C
C=======================================================================
C    Check NBDIM.
C=======================================================================
C
      IF (NBDIM(1).LT.MAXLHS .OR. NBDIM(1).LT.MAXRHS) THEN
        IER = -6
        CALL DTERR(MODE,SUBNAM,IER,NEED)
        RETURN
      ENDIF
      IF (NBDIM(2) .LT. NPTOT/MINNPT) THEN
        IER = -7
        CALL DTERR(MODE,SUBNAM,IER,NEED)
        RETURN
      ENDIF
      IF (NBDIM(3) .LT. NDEP) THEN
        IER = -8
        CALL DTERR(MODE,SUBNAM,IER,NEED)
        RETURN
      ENDIF
C
C=======================================================================
C    Check NWORK.
C=======================================================================
C
      MAXORD = 1 + MAXDEG
      MAXCOF = MAXNPT + MAXDEG + 1
      NEED = 3*(MAXORD+1)*MAXCOF + 5 + MAXORD*(2*MAXORD+1)
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
        CALL DTERR(MODE,SUBNAM,IER,NC)
        RETURN
      ENDIF
C
C=======================================================================
C    Call the appropriate fitting routine.
C=======================================================================
C
       CALL DTGCS1(NPT,X,NDIM,Y,NYDIM,YLHS,NLHS,YRHS,NRHS,NBDIM,NDEG,
     &             NDOM,NDEP,WORK,NWORK,C,IER)
      IF (IER.EQ.0) 
     &  THEN
          C(1) = NDOM
        ELSEIF (IER .GT. 0) THEN
          CALL DTERR(0,SUBNAM,IER,0)
        ELSE
          CALL DTERR(5,SUBNAM,IER,0)
      ENDIF
C
      RETURN
      END
