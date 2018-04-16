      SUBROUTINE DTSUT3 ( C, RANGE )
C
C     =================================================================
C
C     --------------
C     ... PARAMETERS
C     --------------
C
      DOUBLE PRECISION  C(*)         , RANGE(2,*)
C
C     ----------------------
C     ... INTERNAL VARIABLES
C     ----------------------
C
      DOUBLE PRECISION  YMIN         , YMAX         , DENOM        ,
     +                  VALUE        , DTMCON       , DMIN         ,
     +                  DMAX
      INTEGER           ICOEF        , IDOM         , KORD         , 
     1                  NBS          , NCOEF        , NDOM         ,
     1                  NRNG         , IRNG
      LOGICAL           RATNL
C
C    ==================================================================
C
      NDOM  = C(1)
      NRNG  = C(2)
      RATNL = (NRNG .LT. 0)
      IF (RATNL) NRNG = -(NRNG + 1)
      ICOEF = 2 + 3 * NDOM
      NCOEF = 1
C
      DO 10 IDOM = 1, NDOM
C
          KORD  = C(2+IDOM)
          NBS   = C(2+IDOM+NDOM)
          ICOEF = ICOEF + KORD + NBS
          NCOEF = NCOEF * NBS

 10   CONTINUE
C
      IF (RATNL) THEN
        DO 30 IRNG = 1, NRNG
          YMIN = DTMCON(2)
          YMAX = -YMIN
          DO 20 IDOM = 1, NCOEF
            DENOM = C(ICOEF+NRNG*NCOEF+IDOM)
            IF (DENOM .NE. 0.D0) THEN
              VALUE = C(ICOEF+(IRNG-1)*NCOEF+IDOM) / DENOM
            ELSE
              VALUE = C(ICOEF+(IRNG-1)*NCOEF+IDOM)
              ENDIF
            YMIN = MIN (YMIN, VALUE)
            YMAX = MAX (YMAX, VALUE)
 20         CONTINUE
          RANGE(1,IRNG) = YMIN
          RANGE(2,IRNG) = YMAX
 30       CONTINUE
      ELSE
        DO 40 IRNG = 1, NRNG
          RANGE(1,IRNG) = DMIN (NCOEF, C(ICOEF+1), 1)
          RANGE(2,IRNG) = DMAX (NCOEF, C(ICOEF+1), 1)
          ICOEF = ICOEF + NCOEF
 40       CONTINUE
        ENDIF
C
C     =================================================================
C
      RETURN
      END
