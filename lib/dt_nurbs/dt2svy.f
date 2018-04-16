      SUBROUTINE DT2SVY(C,IDOM,BREAK,MULTI,NBREAK,IKNOTS,NKNOTS,IER)
      DOUBLE PRECISION C(*),BREAK(*)
      INTEGER IDOM,MULTI(*),NBREAK,IKNOTS,NKNOTS,IER
C
C     YANK KNOT SETS AND MULTIPLICITIES FROM C ARRAY
C
C  **
C  **               IER      SUCCESS/ERROR CODE
C  **                        =    0 SUCCESS
C  **                        =  -14 IDOM OUT OF RANGE
C  **                        =  -51 INVALID C ARRAY
C  **
C  **   SUBPROGRAMS REFERENCED:  DTMCON   DTERR
C  **
      INTEGER I,N,NORD
      DOUBLE PRECISION DTMCON,XTEMP
      EXTERNAL DTMCON
      EXTERNAL DTERR
C
      CHARACTER*8 SUBNAM
      PARAMETER ( SUBNAM = 'DT2SVY' )
C
      IER = 0
      N = INT(C(1))
      IF ((N .LT. IDOM) .OR. (IDOM .LE. 0)) IER = -14
      IF (N .LT. 1) IER = -51
      IF (IER .LT. 0) THEN
          CALL DTERR(1,SUBNAM,IER,0)
          MULTI(1) = -1
          BREAK(1) = DTMCON(1)
          NBREAK = -1
          IKNOTS = -1
          NKNOTS = -1
          RETURN
      END IF
C
      NORD = C(2+IDOM)
      NKNOTS = C(2+IDOM) + C(2+N+IDOM)
      IKNOTS = 3*N + 2
      DO 20 I = 1,IDOM-1
          IKNOTS = IKNOTS + C(2+I) + C(2+N+I)
  20  CONTINUE
      XTEMP = C(IKNOTS + 1)
      BREAK(1) = XTEMP
      MULTI(1) = 1
      NBREAK = 1
      DO 40 I = 2,NKNOTS
          IF (C(IKNOTS+I) .NE. XTEMP) THEN
              NBREAK = NBREAK + 1
              XTEMP = C(IKNOTS+I)
              BREAK(NBREAK) = XTEMP
              MULTI(NBREAK) = 1
          ELSE
              MULTI(NBREAK) = MULTI(NBREAK) + 1
          END IF
  40  CONTINUE
      IKNOTS = IKNOTS + 1
      RETURN
      END
