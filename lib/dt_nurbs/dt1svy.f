      SUBROUTINE DT1SVY(C,MC,NORD,NBREAK,NKNOTS,NCOEFS,IER)
      DOUBLE PRECISION C(*)
      INTEGER MC,NORD(*),NBREAK(*),NKNOTS(*),NCOEFS,IER
C
C     YANK KNOT SETS AND MULTIPLICITIES FROM C ARRAY
C
C  **
C  **               IER      SUCCESS/ERROR CODE
C  **                        =    0 SUCCESS
C  **                        =   -1 NORD TOO SMALL
C  **                        =   -8 INVALID XKNT ARRAY
C  **                        =  -15 TOO FEW KNOTS
C  **                        =  -51 INVALID C(1)
C  **                        =  -52 INVALID C(2)
C  **
C  **    SUBPROGRAMS REFERENCED:  DTERR
C  **
      INTEGER I,J,MULTI,N,NBRK,NKNTS,NCFS,NRD
      DOUBLE PRECISION CREAL,XKNOT
      EXTERNAL DTERR
C
      CHARACTER*8 SUBNAM
      PARAMETER ( SUBNAM = 'DT1SVY' )
C
      IER = 0
      N = INT(C(1))
      CREAL = DBLE(N)
      IF ((N .LT. 1) .OR. (CREAL .NE. C(1))) IER = -51
      M = INT(C(2))
      CREAL = DBLE(M)
      IF ((M .LT. 1) .OR. (CREAL .NE. C(2))) IER = -52
      IF (IER .LT. 0) GO TO 900
C
      MC = 3*N + 2
      NCOEFS = 1
      DO 120 I = 1,N
          NRD = INT( C(2+I) )
          NCFS = INT( C(2+N+I) )
          NKNTS = NRD + NCFS
          NCFS = MAX(1,NCFS)
          NCOEFS = NCOEFS*NCFS
          NBRK = 1
          XKNOT = C(MC+1)
          MULTI = 1
          DO 140 J = 2,NKNTS
              IF (C(MC+J) .NE. XKNOT) THEN
                  NBRK = NBRK + 1
                  XKNOT = C(MC+J)
                  MULTI = 1
              ELSE
                  MULTI = MULTI + 1
                  IF (MULTI .GT. NRD) IER = -8
              END IF
              IF (C(MC+J) .LT. C(MC+J-1)) IER = -8
  140     CONTINUE
          CREAL = DBLE( NRD )
          IF ((NRD .LT. 1) .OR. (CREAL .NE. C(2+I))) IER = -1
          IF (NKNTS .LT. 2*NRD) IER = -15
          MC = MC + NKNTS
          NBREAK(I) = NBRK
          NKNOTS(I) = NKNTS
          NORD(I) = NRD
          IF (IER .LT. 0) GO TO 900
 120  CONTINUE
      MC = MC + M*NCOEFS
 900  IF (IER .LT. 0)  THEN
          CALL DTERR(1,SUBNAM,IER,0)
          MC = -1
          NORD(1) = -1
          NBREAK(1) = -1
          NKNOTS(1) = -1
          NCOEFS = -1
      END IF
      RETURN
      END
