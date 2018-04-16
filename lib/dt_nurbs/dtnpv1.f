C<*>
      SUBROUTINE DTNPV1 ( X,      INCX,   NRNG,   C,
     1                    NCOEF,  COEF,   NDIM,   INDX,
     2                    AT,     BSVAL,  WORK,   V,
     3                    IER       )
C
C     ==================================================================
C
C     --------------
C     ... PARAMETERS
C     --------------
C
      INTEGER      INCX,     NRNG,     NCOEF,    NDIM,
     1             INDX   ,  IER
C
      DOUBLE PRECISION     X(*),     C(*),     COEF(NCOEF,*),
     1                     AT(NDIM,*),  BSVAL(*), WORK(*),  V(*)
C
C     ----------------------
C     ... INTERNAL VARIABLES
C     ----------------------
C
      INTEGER      I,        IK,       IKPT,     INCZ,     IPOSX,
     1             IPOSY,    IX,       J,        KX,       KY,
     2             NBSX,     NBSY
C
C     -------------
C     ... EXTERNALS
C     -------------
C
      DOUBLE PRECISION       DDOT
      EXTERNAL     DTBSP2,   DTNPV3,   DDOT
C
C     ==================================================================
C
      INCZ = 1
C
C     ---------------
C     ... LOCATE X(1)
C     ---------------
C
      IKPT  = 9
      KX    = INT ( C(3) )
      NBSX  = INT ( C(5) )
      IPOSX = INT ( C(7) )
C
      CALL DTNPV3 ( X(1), C(IKPT), NBSX, KX, IPOSX, IER )
C
      C(7) = DBLE ( IPOSX )
C
      IF ( IER .NE. 0) THEN
C
          RETURN
C
      END IF
C
C     ------------------------------
C     ... CHECK FOR INVALID KNOT SET
C     ------------------------------
C
      CALL DTILC1( C(IKPT + IPOSX - KX), 2 * KX, KX, IFAIL)
C
      IF ( IFAIL .NE. 0) THEN
C
          IER = -8
C
          RETURN
C
      END IF
C
C     ----------------
C     ... LOCATE X(IX)
C     ----------------
C
      IKPT  = IKPT + NBSX + KX
      KY    = INT ( C(4) )
      NBSY  = INT ( C(6) )
      IPOSY = INT ( C(8) )
      IX    = 1 + INCX
C
      CALL DTNPV3 ( X(IX), C(IKPT), NBSY, KY, IPOSY, IER )
C
      C(8) = DBLE ( IPOSY )
C
      IF ( IER .NE. 0) THEN
C
          RETURN
C
      END IF
C
C     ------------------------------
C     ... CHECK FOR INVALID KNOT SET
C     ------------------------------
C
      CALL DTILC1( C(IKPT + IPOSY - KY), 2 * KY, KY, IFAIL)
C
      IF ( IFAIL .NE. 0) THEN
C
          IER = -8
C
          RETURN
C
      END IF
C
C     ---------------------------
C     ... EVALUATE IN X DIRECTION
C     ---------------------------
C
      IK = -1
C
      CALL DTBSP2 ( C(9),       X(1), IPOSX, IK, KX, WORK(1),
     1              WORK(KX), BSVAL                         )
C
C     ------------------------
C     ... COLLAPSE X DIRECTION
C     ------------------------
C
      INDX = ( IPOSY - KY ) * NBSX + IPOSX - KX + 1
C
      DO 20 I = 1, KY
C
          DO 10 J = 1, NRNG

              AT(I,J) = DDOT (KX, COEF(INDX,J), INCZ, BSVAL, INCZ)
C
 10       CONTINUE
C
          INDX = INDX + NBSX
C
 20   CONTINUE
C
C     ---------------------------
C     ... EVALUATE IN Y DIRECTION
C     ---------------------------
C
      IK = -1
C
      CALL DTBSP2( C(IKPT),    X(IX), IPOSY, IK, KY, WORK(1),
     1             WORK(KY), BSVAL                             )
C
      DO 30 J = 1, NRNG
C
          V(J) = DDOT ( KY, AT(1,J), INCZ, BSVAL, INCZ )
C
 30   CONTINUE
C
C     ==================================================================
C
C     ----------
C     ... RETURN
C     ----------
C
      RETURN
      END
