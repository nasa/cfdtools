      SUBROUTINE DTSPRM ( CIN, INDX, IWORK, NWORK, COUT, IER )
C
C     =================================================================
C
C     --------------
C     ... PARAMETERS
C     --------------
C
      INTEGER           INDX(*)       , IWORK(*)       , NWORK        ,
     1                  IER
C
      DOUBLE PRECISION  CIN(*)        , COUT(*)
C
C     ----------------------
C     ... INTERNAL VARIABLES
C     ----------------------
C
      INTEGER           I             , IDOM           , IFAIL        ,
     1                  IKIN          , IKOT           , ILAST        ,
     2                  INCX          , ISHFT          , ISTRT        , 
     3                  JDOM          , MODE           , NCOEF        ,
     4                  NDOM          , NEED           , NKNTS        ,
     5                  NPRM          , NRNG           ,  NSHFT
C                                                                 
      CHARACTER*8       SUBNAM
C
C     ------------
C     ... EXTERNAL
C     ------------
C
      DOUBLE PRECISION  DTMCON
C 
C     =================================================================
C
      SUBNAM = 'DTSPRM  '
      IER    = 0
      MODE   = 3
      NEED   = 0
C
C     ------------------------
C     ... CHECK VALID C VECTOR
C     ------------------------
C
      CALL DTSCHK ( CIN, IFAIL )
C
      IF ( IFAIL .NE. 0 ) THEN
C
          IER     = -1
          COUT(1) = DTMCON(2)
          CALL DTERR ( MODE, SUBNAM, IER, NEED )
          RETURN
C
      END IF
C
      NDOM = INT ( CIN(1) )
      NRNG = INT ( CIN(2) )
C
C     --------------------
C     ... CHECK VALID INDX
C     --------------------
C
      DO 10 IDOM = 1, NDOM
C
          IF ( INDX(IDOM) .LT. 1 .OR. INDX(IDOM) .GT. NDOM ) THEN
C
              IER     = -2
              COUT(1) = DTMCON(2)
              CALL DTERR ( MODE, SUBNAM, IER, NEED )
              RETURN
C
          END IF
C 
 10   CONTINUE
C
C     ------------------------------
C     ... DETERMINE WORKSPACE NEEDED
C     ------------------------------
C
      NCOEF = 1
      NPRM  = 0
      NKNTS = 0
      INCX  = 1
C
      DO 20 IDOM = 1, NDOM
C
          NCOEF = NCOEF * CIN(2+IDOM+NDOM)
          NKNTS = NKNTS + CIN(2+IDOM) + CIN(2+IDOM+NDOM)
          IF ( INDX(IDOM) .NE. IDOM) NPRM = NPRM + 1
C
 20   CONTINUE
C
C     --------------------------------------------
C     ... IF NO PERMUTATION SPECIFIED PERFORM COPY
C     --------------------------------------------
C
      IF ( NPRM .EQ. 0 ) THEN
C
         NPRM  = 2 + 3 * NDOM + NKNTS + NRNG * NCOEF
         NKNTS = 1
         CALL DCOPY ( NPRM, CIN, INCX, COUT, INCX )
         RETURN
C
      END IF
C
      NEED =  NCOEF
C
      IF ( NWORK .LT. NEED ) THEN
C
          IER     = -5
          COUT(1) = DTMCON(2)
          MODE    = 4
          CALL DTERR ( MODE, SUBNAM, IER, NEED )
C
      END IF 
C
C     ----------------------------
C     ... PERMUTE KNOT INFORMATION
C     ----------------------------
C
      COUT(1) = CIN(1)
      COUT(2) = CIN(2)
C
      DO 30 IDOM = 1, NDOM
C
          COUT(2+IDOM)       = CIN(2+INDX(IDOM))
          COUT(2+IDOM+NDOM)  = CIN(2+INDX(IDOM)+NDOM)
          COUT(2+IDOM+2*NDOM) = CIN(2+INDX(IDOM)+2*NDOM)
C
 30   CONTINUE
C
      IKOT = 3 + 3 * NDOM
C
      DO 50 IDOM = 1, NDOM
C
          NKNTS = COUT(2+IDOM) + COUT(2+IDOM+NDOM)
          IKIN  = 3 + 3 * NDOM
C
          DO 40 JDOM = 1, INDX(IDOM) - 1
C
              IKIN = IKIN + INT ( CIN(2+JDOM) + CIN(2+JDOM+NDOM) )
C
 40       CONTINUE
C
          CALL DCOPY ( NKNTS, CIN(IKIN), INCX, COUT(IKOT),
     1                 INCX )
          IKOT = IKOT + NKNTS
C
 50   CONTINUE
C
C     ---------------------------------------------
C     ... DEFINE PERMUTATION INDEX FOR COEFFICIENTS
C     ---------------------------------------------
C
      ISTRT    = 1
      ILAST    = 1
      ISHFT    = 1
      IWORK(1) = 1
      IKIN     = 2 + 3 * NDOM
C
      DO 80 IDOM = 1, NDOM
C
          ISTRT = ILAST + 1
          ILAST = ILAST * INT ( COUT(2+IDOM+NDOM) )
          NSHFT = 1
          IKIN  = IKIN + INT (CIN(2+IDOM) + CIN(2+IDOM+NDOM))
C
          DO 60 JDOM = 1, INDX(IDOM) - 1
C
              NSHFT = NSHFT * INT ( CIN(2+JDOM+NDOM) )
C
 60       CONTINUE
C
          DO 70 I = ISTRT, ILAST
C        
              IWORK(I) = IWORK(I-ISHFT) + NSHFT
C
 70       CONTINUE
C
          ISHFT = ILAST
C
 80   CONTINUE
C
C     ----------------------------
C     ... PERMUTE THE COEFFICIENTS
C     ----------------------------
C
      DO 100 IDOM = 1, NRNG
C
          DO 90 I = 1, NCOEF
C
              COUT(IKIN+I) = CIN(IKIN+IWORK(I))
C          
 90       CONTINUE
C
          IKIN = IKIN + NCOEF
C
 100  CONTINUE
C
C     ========================================================
C
      RETURN
      END
