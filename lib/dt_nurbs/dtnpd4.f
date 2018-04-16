C<*>
      SUBROUTINE DTNPD4( KNOTS, X, IPOS, KORD, NDER, WORK, NDIM, BSVAL )
C
C     ==================================================================
C
C     --------------
C     ... PARAMETERS
C     --------------
C
      INTEGER           IPOS    , KORD    , NDER    , NDIM
C
      DOUBLE PRECISION  KNOTS(*), X       , WORK(*) , BSVAL(NDIM,*)
C
C     -----------------------
C     ... INTERNAL PARAMETERS
C     -----------------------
C
      INTEGER           INCX    , INCY    , ND
C
      DOUBLE PRECISION  ZERO
C
C     -------------
C     ... EXTERNALS
C     -------------
C
      EXTERNAL          DTBSP1
C
C     ==================================================================
C
C     BUG FIX 05/31/90: REPAIRED HANDLING OF REQUESTS FOR DERIVATIVES
C     HIGHER THAN THE DEGREE OF THE SPLINE. - ERIC BRECHNER
C
C     ==================================================================
C
      INCX = 0
      INCY = 1
      ZERO = 0.0D0
C
C     -------------------------------
C     ... DEFINE B-SPLINE DERIVATIVES
C     -------------------------------
C
      ND = MIN (KORD - 1, NDER)
C
      CALL DTBSP1 ( KNOTS, X, IPOS, KORD, ND, WORK, BSVAL, NDIM )
C
C     ==================================================================
C
C     ----------
C     ... RETURN
C     ----------
C
      RETURN
      END
