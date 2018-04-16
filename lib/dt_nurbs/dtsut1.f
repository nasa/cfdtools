      SUBROUTINE DTSUT1 ( C, CSIZE )
C
C     ================================================================
C
C     --------------
C     ... PARAMETERS
C     --------------
C
      INTEGER           CSIZE
C
      DOUBLE PRECISION  C(CSIZE)
C
C     ----------------------
C     ... INTERNAL VARIABLES
C     ----------------------
C
      INTEGER           IDOM          , NCOEF        , NKNTS       ,
     1                  NDOM          , NRNG
C
C     ================================================================
C
      NDOM = INT ( C(1) )
      NRNG = INT ( ABS( C(2) ) )
      NCOEF = 1
      NKNTS = 0
C
      DO 10 IDOM = 1, NDOM
C
          NCOEF = NCOEF * INT ( C(IDOM+NDOM+2) )
          NKNTS = NKNTS + INT ( C(IDOM+NDOM+2) )
     1                  + INT ( C(IDOM+2) )
C
 10   CONTINUE                                       

C
      CSIZE = 2 + 3 * NDOM + NKNTS + NRNG * NCOEF
C
C     ================================================================
C
      RETURN
      END
