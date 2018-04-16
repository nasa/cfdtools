      SUBROUTINE DTGTPA(X,NDIM,NDOM,NPT,Y,NYDIM,NDEP,NDEG,
     &                  IWT,WHT,XKNOTS,NKNOTS,NDIMK,
     &                  WORK,NWORK,MAXC,C,NC,IFAIL,IER)
C
C=======================================================================
C
C  PURPOSE:
C          DTGTPA is used to produce multivariate least square spline
C          approximations to data on a regular grid.
C
C  USAGE:
C          DOUBLE PRECISION X(NDIM,NDOM),Y(NPT(1),...,NPT(NDOM),NDEP),
C                           WHT(NDIM,NDOM),XKNOTS(NDIMK,NDOM),
C                           WORK(NWORK),C(MAXC)
C          INTEGER  NPT(NDOM),NYDIM(NDOM),NDEG(NDOM),NKNOTS(NDOM)
C          CALL DTGTPA(X,NDIM,NDOM,NPT,Y,NYDIM,NDEP,NDEG,IWT,WHT,
C                      XKNOTS,NKNOTS,NDIMK,WORK,NWORK,MAXC,C,NC,
C                      IFAIL,IER)
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
C          NDEG    Array containing the degree of the spline to
C                  be used for each independent variable.
C
C          IWT     Weight option.
C
C                  IWT .EQ. 0  DTGTPA uses WHT(I,J) = 1.0 for all I
C                              and J=1,NDOM.  The WHT array is not used.
C
C                  IWT .NE. 0  DTGTPA uses weights provided in the
C                              WHT array.
C
C          WHT     Array of nonnegative weights in correspondence
C                  to the X array, WHT(I,J) .GE. 0.  If IWT .EQ. 0, 
C                  WHT is not used and it may be a dummy argument.  If
C                  IWT .NE. 0, at least NKNOTS + K of the WHT(I,J) must
C                  be positive for each J.
C
C          XKNOTS  Array containing the internal knots in each 
C                  direction in ascending order with multiplicity.
C
C          NKNOTS  Array containing the number of knots in each
C                  independent variable direction.
C
C          NDIMK   Leading dimension of the XKNOTS array.
C
C  STORAGE:
C          WORK    The work array.
C
C          NWORK   Size of the work array.
C                  NWORK .GE.  MP + 5 + 2*NC + MO + (NC+3)*MO + ND
C                    where  MO = 1 + MAX(NDEG(I))
C                           MP = MAX(NPT(I))
C                           NC = MAX(1+NDEG(I)+NKNOTS(I))
C                           ND = number of data points (the Y array
C                                is copied and thus not destroyed).
C
C  OUTPUT:
C          C       The spline definition array for the result.
C
C          NC      Length of spline array.
C
C          IFAIL   Knot interval where interlacing fails.  Used only if
C                  IER = -35.
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
C                  IER = -6   NKNOTS(I) .LT. 0.
C
C                  IER = -8   XKNOTS(.,I) not in ascending order or the
C                             multiplicity of a knot exceeds the order.
C
C                  IER = -9   NDOM .LT. 1.
C
C                  IER = -14  The user has not supplied enough storage
C                             for the spline vector, C.
C
C                  IER = -20  The matrix representing the least square
C                             problem is singular.
C
C                  IER = -35  XKNOTS failed interlacing conditions.
C
C                  IER = -36  An element of the knot array XKNOTS is
C                             outside the interval defined by the
C                             1st and last data point.
C
C
C  OTHER DT ROUTINES CALLED
C
C        DTERR
C        DTGTP1
C
C     ******************************************************************
C
C  Parameters:
C
      CHARACTER*8 SUBNAM
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.D0,SUBNAM = 'DTGTPA  ')
C
C  Arguments:
C
      INTEGER  NPT(*), NYDIM(*), NDEG(*), NDIM, NDEP, NDOM,
     &         IWT, NKNOTS(*), NDIMK, NWORK, IER
      DOUBLE PRECISION  X(NDIM,*), Y(*), WHT(NDIM,*), XKNOTS(NDIMK,*), 
     &                  WORK(NWORK), C(MAXC)
C
C  Internal:
C
      INTEGER  IX,IV,MODE,NEED,MAXNPT,MAXORD,MAXCOF,NTOT,NDATA
      INTEGER  KORD, NCOEF, I, NCI, IK, IK0, NPTS, K
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
      MAXCOF = 0
      NDATA = 1
      DO 20 IV = 1,NDOM
        MAXNPT = MAX(MAXNPT,NPT(IV))
        NDATA = NDATA*NPT(IV)
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
        KORD = 1 + NDEG(IV)
        MAXORD = MAX(MAXORD,KORD)
        MAXCOF = MAX(MAXCOF,NKNOTS(IV)+KORD)
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
C
   20 CONTINUE
C
C=======================================================================
C    Check NWORK.
C=======================================================================
C
      NEED = MAXNPT + (5 + 2*MAXCOF + MAXORD)
     &       + MAXORD*(MAXCOF + 3) + NDATA
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
        NCI = NKNOTS(I) + KORD
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
C    Set the C-vector.
C=======================================================================
C
      C(2) = DBLE( FLOAT(NDEP) )
      IK = 3*NDOM + 2
C
      DO 200 IV=1,NDOM
        K = NDEG(IV) + 1
        C(2+IV) = DBLE( FLOAT(K) )
        NCOEF = NKNOTS(IV) + K
        C(2+NDOM+IV) = DBLE( FLOAT(NCOEF) )
        C(2+2*NDOM+IV) = C(2+NDOM+IV)
        IK0 = IK + 1
C
C=======================================================================
C    Set knots.
C=======================================================================
C
C  Left endpoint:
C
        DO 90 I = 1, K
          IK = IK + 1
          C(IK) = X(1,IV)
   90   CONTINUE
C
C  Interior:
C
        NPTS = NPT(IV)
        DO 110 I = 1, NKNOTS(IV)
          IF (X(1,IV).LT.XKNOTS(I,IV) .AND. XKNOTS(I,IV).LT.X(NPTS,IV))
     &      THEN
              IK = IK + 1
              C(IK) = XKNOTS(I,IV)
            ELSE
              IER = -36
              CALL DTERR(MODE,SUBNAM,IER,NEED)
              RETURN
          ENDIF
  110   CONTINUE
C
C  Right endpoint:
C
        DO 130 I = 1, K
          IK = IK + 1
          C(IK) = X(NPTS,IV)
  130   CONTINUE
C
C=======================================================================
C     Call DTILCK to check interlacing.
C=======================================================================
C
        NTOT = NCOEF + K
        CALL DTILCK(NPTS,X(1,IV),K,IWT,WHT(1,IV),NTOT,C(IK0),IFAIL,IER)
        IF (IER.NE.0) THEN
          IER = -100
          MODE = 5
          CALL DTERR(MODE,SUBNAM,IER,NEED)
          RETURN
        ENDIF
        IF (IFAIL.NE.0) THEN
          IF( IFAIL .LT. -100 ) IER = -8
          IF( IFAIL .GT. 0 ) IER = -35
          CALL DTERR(MODE,SUBNAM,IER,NEED)
          RETURN
        ENDIF
C
  200 CONTINUE
C
C=======================================================================
C    Call the fitting routine.
C=======================================================================
C
       CALL DTGTP1(NPT,X,NDIM,Y,NYDIM,NDEG,NDOM,NDEP,IWT,WHT,
     &             NDIMK,WORK,NWORK,C,IER)
C
      IF (IER .EQ. 0)
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
