      SUBROUTINE D0NPM1(N,X,Y,NDIM,NDEG,NRNG,ICC,MLT,HOLD,NHOLD,C,IER)
C
C=======================================================================
C
C  PURPOSE  D0NPM1 WILL DEFINE THE MULTIPLICITY VECTOR FOR A SIMPLE
C          SPLINE INTERPOLANT.
C
C  USAGE   DOUBLE PRECISION X(N),Y(NDIM,N),HOLD(NHOLD),C(MC)
C          INTEGER MLT(N-2)
C          CALL D0NPM1(N,X,Y,NDIM,NDEG,NRNG,ICC,MLT,HOLD,NHOLD,C,IER)
C
C  OTHER Library ROUTINES CALLED
C
C       D0PITG
C       HSBSPL
C       DTERR
C
C=======================================================================
C
C
C=======================================================================
C     PARAMETERS
C=======================================================================
C
      INTEGER ICC, IER, MLT(*), N, NDIM, NDEG, NHOLD, NRNG
      DOUBLE PRECISION C(*), X(*), Y(*), HOLD(*)
C
C=======================================================================
C      INTERNAL VARIABLES
C=======================================================================
C
      DOUBLE PRECISION BCLHS(1), BCRHS(1)
      INTEGER I, INDLHS(1), INDRHS(1), NLHS, NM2, NRHS, MDIM
C
C=======================================================================
C    SET MULTIPLICITY VECTOR AND CALL D0PITG
C=======================================================================
C
      NM2 = N - 2
      IF( NM2 .EQ. 0 ) GO TO 20
      DO 10 I = 1, NM2
          MLT(I) = 1
10    CONTINUE
20    NRHS = 0
      NLHS = 0
      MDIM = 1
      CALL D0PITG(N,X,Y,NDIM,NDEG,NRNG,ICC,MLT,MDIM,NLHS,INDLHS,BCLHS,
     +            NRHS,INDRHS,BCRHS,HOLD,NHOLD,C,IER)
      RETURN
      END
