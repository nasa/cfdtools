      SUBROUTINE  d0xcpy(N,D,istx,INCX,isty,INCY)
C
C     COPIES part of A VECTOR to another part.
C
c     a call to this routine should be equivalent to
c     'call dcopy( n, d(istx), incx, d(isty), incy )'
c     if dcopy is not vectorized

      DOUBLE PRECISION D(*)
      INTEGER I,INCX,INCY,istx,isty,IX,IY,N
C
      IF(N.LE.0)RETURN
      IX = istx
      IY = isty
      IF(INCX.LT.0)IX = (-N+1)*INCX + istx
      IF(INCY.LT.0)IY = (-N+1)*INCY + isty
      DO 10 I = 1,N
        D(IY) = D(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
      END
