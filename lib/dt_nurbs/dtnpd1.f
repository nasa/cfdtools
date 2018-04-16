C<*>
      SUBROUTINE DTNPD1 ( X, IDER1, IDER2, C, AT, KMAX, BSVAL,
     *                  WORK, V, IER)
C
C  **           THIS SUBROUTINE CALCULATES THE PARTIAL DERIVATIVES
C  **           OF A BIVARIATE M-COMPONENT SPLINE FUNCTION.
C  **
C  **           INPUT:          X       = (X(1), X(2)) THE POINT
C  **                                      OF EVALUATION
C  **
C  **                           IDER1   = THE INDEX OF THE FIRST
C  **                                     DESIRED DERIVATIVE
C  **
C  **                           IDER2   = THE INDEX OF THE SECOND
C  **                                     DESIRED DERIVATIVE
C  **
C  **                           C       = ARRAY OF SPLINE PARAMETERS
C  **                                     C(1) = 2 AND C(2) > 0
C  **
C  **                           AT      = A SCRATCH ARRAY OF DIMENSION
C  **                                     KMAX X KMAX X C(2)
C  **
C  **                           KMAX    = MAX(C(3), C(4))
C  **
C  **                           BSVAL   = A SCRATCH ARRAY OF DIMENSION
C  **                                     KMAX X MAXDR WHERE MAXDR =
C  **                                     MAX(IDER(1)+1, IDER(2)+1))
C  **
C  **                           WORK    = A SCRATCH ARRAY OF DIMENSION
C  **                                     KMAX X KMAX
C  **
C  **           OUTPUT:         V       = AN ARRAY OF DIMENSION
C  **                                     IDER1+1 X IDER2+1 X C(2)
C  **
C  **                           IER     = ERROR FLAG
C  **                                     0  SUCCESS
C  **                                     ELSE COMPUTATIONAL FAILURE
C  **
C  ** WRITTEN BY
C  **     DAVE FERGUSON       JULY, 1989
C     ==================================================================
C
C     --------------
C     ... PARAMETERS
C     --------------
C
      INTEGER      IDER1, IDER2,  KMAX,  IER
C
      DOUBLE PRECISION       X(*),     C(*),
     1             AT(KMAX, KMAX, *),  BSVAL(KMAX, *)
     2             ,WORK(*),  V(IDER1+1, IDER2+1, *)
C
C     ----------------------
C     ... INTERNAL VARIABLES
C     ----------------------
C
      INTEGER      I,        IKPT,     INCZ,     IPOSX,
     1             IPOSY,    IX,       KX,       KY,
     2             NBSX,     NBSY,     NRNG,     INCX,     INDX,
     3             INDR,     IOVER,    NCOEF
      INTEGER      NDERX, NDERY, JD, JR, JDX, JDY, J, K
C
C     -------------
C     ... EXTERNALS
C     -------------
C
      DOUBLE PRECISION       DDOT
      EXTERNAL     DTNPD3,   DTNPD4,   DDOT
C
C     ==================================================================
C
      INCZ = 1
      INCX = 1
C
      NRNG = C(2)
C
C     INITIALIZE V TO ZERO
C
      DO 5 I = 1, IDER1+1
          DO 5 J = 1, IDER2+1
              DO 5 K = 1, NRNG
                  V(I,J,K) = 0.0
    5 CONTINUE
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
      CALL DTNPD3 ( X(1), C(IKPT), NBSX, KX, IPOSX, IER )
C
      C(7) = DBLE( IPOSX )
C
      IF ( IER .NE. 0 ) THEN
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
      CALL DTNPD3 ( X(IX), C(IKPT), NBSY, KY, IPOSY, IER )
C
      C(8) = DBLE( IPOSY )
C
      IF ( IER .NE. 0 ) THEN
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
C     ----------------------
C     ... EVALUATE B-SPLINES
C     ----------------------
C
      CALL DTNPD4 ( C(9), X(1), IPOSX, KX, IDER1, WORK, KMAX, BSVAL)
C
      INDX = ( IPOSY - KY ) * NBSX + IPOSX - KX + 1
C
      NDERX = MIN0(IDER1, KX-1)
      IOVER = 8 + NBSX + NBSY + KX + KY
      NCOEF = NBSX * NBSY
        DO 100 I = 1, KY
                DO 110 JD = 0, NDERX
                        DO 120 JR = 1, NRNG
                        INDR = NCOEF * (JR - 1) + INDX + IOVER
                        AT(I, JD+1, JR) = DDOT(KX, C(INDR), INCZ
     *                                     ,BSVAL(1, JD+1), INCZ)
120                     CONTINUE
110             CONTINUE
        INDX = INDX + NBSX
100     CONTINUE
C     ---------------------------
C     ... EVALUATE IN Y DIRECTION
C     ---------------------------
C
C     ----------------------
C     ... EVALUATE B-SPLINES
C     ----------------------
C
C
      CALL DTNPD4 ( C(IKPT), X(IX), IPOSY, KY, IDER2, WORK, KMAX,
     1              BSVAL)
C
        NDERY = MIN0(IDER2, KY-1)
        DO 200 JDX = 0, NDERX
                DO 210 JDY = 0, NDERY
                    DO 220 JR = 1, NRNG
                    V(JDX+1, JDY+1, JR) = DDOT(KY, AT(1, JDX+1, JR)
     *                               ,INCZ, BSVAL(1, JDY+1), INCZ)
220                 CONTINUE
210             CONTINUE
200     CONTINUE
C
C     ==================================================================
C
C     ----------
C     ... RETURN
C     ----------
C
      RETURN
      END
