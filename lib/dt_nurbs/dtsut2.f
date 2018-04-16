      SUBROUTINE DTSUT2 ( C, IDOM, SXMIN, SXMAX )
C
C     ==================================================================
C
C     --------------
C     ... PARAMETERS
C     --------------
C
      INTEGER           IDOM
C
      DOUBLE PRECISION  C(*)          , SXMIN        , SXMAX
C
C     ----------------------
C     ... INTERNAL VARIABLES
C     ----------------------
C
      INTEGER           I             , ILK          , KORD          , 
     1                  NBS          , NDOM
C
C     ==================================================================
C
      NDOM = INT ( C(1) )
C
      ILK  = 3 + 3 * NDOM
C
C     ------------------------------------
C     ... INCREMENT TO POINT TO FIRST KNOT
C     ------------------------------------
C
      DO 10 I = 1, IDOM - 1
C
          KORD = INT ( C(2+I) )
          NBS  = INT ( C(2+NDOM+I) )
          ILK  = ILK + NBS + KORD
C
 10   CONTINUE
C
      KORD = INT ( C(2+IDOM) )
      NBS   = INT ( C(2+NDOM+IDOM) )
      SXMIN = C(ILK+KORD-1)
      SXMAX = C(ILK+NBS)
C
C     ==================================================================
C
      RETURN
      END
