*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     File  blas2.f
*     NAG versions of the Level 2  BLAS Matrix-vector routines
*
*     dgemv           dger            dsymv           dsyr     
*     dtrmv           dtrsv
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine dgemv ( trans, m, n, alpha, a, lda, x, incx,
     $                   beta, y, incy )

C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     Modified by PEG 8/31/94.
C     .. Entry Points ..
C     ENTRY      F06PAF( TRANS, M, N, ALPHA, A, LDA, X, INCX,
C    $                   BETA, Y, INCY )

C     .. Scalar Arguments ..
      double precision   alpha, beta
      integer            incx, incy, lda, m, n
      character*1        trans
C     .. Array Arguments ..
      double precision   a( lda, * ), x( * ), y( * )
C     ..
C
C  Purpose
C  =======
C
C  DGEMV  performs one of the matrix-vector operations
C
C     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
C
C  where alpha and beta are scalars, x and y are vectors and A is an
C  m by n matrix.
C
C  Parameters
C  ==========
C
C  TRANS  - CHARACTER*1.
C           On entry, TRANS specifies the operation to be performed as
C           follows:
C
C              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
C
C              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
C
C              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
C
C           Unchanged on exit.
C
C  M      - INTEGER.
C           On entry, M specifies the number of rows of the matrix A.
C           M must be at least zero.
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the number of columns of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  ALPHA  - DOUBLE PRECISION.
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
C           Before entry, the leading m by n part of the array A must
C           contain the matrix of coefficients.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, m ).
C           Unchanged on exit.
C
C  X      - DOUBLE PRECISION array of DIMENSION at least
C           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
C           and at least
C           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
C           Before entry, the incremented array X must contain the
C           vector x.
C           Unchanged on exit.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C  BETA   - DOUBLE PRECISION.
C           On entry, BETA specifies the scalar beta. When BETA is
C           supplied as zero then Y need not be set on input.
C           Unchanged on exit.
C
C  Y      - DOUBLE PRECISION array of DIMENSION at least
C           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
C           and at least
C           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
C           Before entry with BETA non-zero, the incremented array Y
C           must contain the vector y. On exit, Y is overwritten by the
C           updated vector y.
C
C  INCY   - INTEGER.
C           On entry, INCY specifies the increment for the elements of
C           Y. INCY must not be zero.
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C
C
C     .. Parameters ..
      double precision   one         , zero
      parameter        ( one = 1.0d+0, zero = 0.0d+0 )
C     .. Local Scalars ..
      double precision   temp
      integer            i, info, ix, iy, j, jx, jy, kx, ky, lenx, leny
C     .. External Subroutines ..
      external           f06aaz
C     .. Intrinsic Functions ..
      intrinsic          max
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      info = 0
      if     ( .not.(trans.eq.'N' .or. trans.eq.'n').and.
     $         .not.(trans.eq.'T' .or. trans.eq.'t').and.
     $         .not.(trans.eq.'C' .or. trans.eq.'c')      )then
         info = 1
      else if( m.lt.0 )then
         info = 2
      else if( n.lt.0 )then
         info = 3
      else if( lda.lt.max( 1, m ) )then
         info = 6
      else if( incx.eq.0 )then
         info = 8
      else if( incy.eq.0 )then
         info = 11
      end if
      if( info.ne.0 )then
         call f06aaz( 'dgemv/f06paf ', info )
         return
      end if
C
C     Quick return if possible.
C
      if( ( m.eq.0 ).or.( n.eq.0 ).or.
     $    ( ( alpha.eq.zero ).and.( beta.eq.one ) ) )
     $   return
C
C     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
C     up the start points in  X  and  Y.
C
      if( (trans.eq.'N' .or. trans.eq.'n') )then
         lenx = n
         leny = m
      else
         lenx = m
         leny = n
      end if
      if( incx.gt.0 )then
         kx = 1
      else
         kx = 1 - ( lenx - 1 )*incx
      end if
      if( incy.gt.0 )then
         ky = 1
      else
         ky = 1 - ( leny - 1 )*incy
      end if
C
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through A.
C
C     First form  y := beta*y.
C
      if( beta.ne.one )then
         if( incy.eq.1 )then
            if( beta.eq.zero )then
               do 10, i = 1, leny
                  y( i ) = zero
   10          continue
            else
               do 20, i = 1, leny
                  y( i ) = beta*y( i )
   20          continue
            end if
         else
            iy = ky
            if( beta.eq.zero )then
               do 30, i = 1, leny
                  y( iy ) = zero
                  iy      = iy   + incy
   30          continue
            else
               do 40, i = 1, leny
                  y( iy ) = beta*y( iy )
                  iy      = iy           + incy
   40          continue
            end if
         end if
      end if
      if( alpha.eq.zero )
     $   return
      if( (trans.eq.'N' .or. trans.eq.'n') )then
C
C        Form  y := alpha*A*x + y.
C
         jx = kx
         if( incy.eq.1 )then
            do 60, j = 1, n
               if( x( jx ).ne.zero )then
                  temp = alpha*x( jx )
                  do 50, i = 1, m
                     y( i ) = y( i ) + temp*a( i, j )
   50             continue
               end if
               jx = jx + incx
   60       continue
         else
            do 80, j = 1, n
               if( x( jx ).ne.zero )then
                  temp = alpha*x( jx )
                  iy   = ky
                  do 70, i = 1, m
                     y( iy ) = y( iy ) + temp*a( i, j )
                     iy      = iy      + incy
   70             continue
               end if
               jx = jx + incx
   80       continue
         end if
      else
C
C        Form  y := alpha*A'*x + y.
C
         jy = ky
         if( incx.eq.1 )then
            do 100, j = 1, n
               temp = zero
               do 90, i = 1, m
                  temp = temp + a( i, j )*x( i )
   90          continue
               y( jy ) = y( jy ) + alpha*temp
               jy      = jy      + incy
  100       continue
         else
            do 120, j = 1, n
               temp = zero
               ix   = kx
               do 110, i = 1, m
                  temp = temp + a( i, j )*x( ix )
                  ix   = ix   + incx
  110          continue
               y( jy ) = y( jy ) + alpha*temp
               jy      = jy      + incy
  120       continue
         end if
      end if
C
C     end of dgemv (f06paf)
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine dger  ( m, n, alpha, x, incx, y, incy, a, lda )

C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     Modified by PEG 8/31/94.
C     .. Entry Points ..
C     ENTRY      F06PMF( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )

C     .. Scalar Arguments ..
      double precision   alpha
      integer            incx, incy, lda, m, n
C     .. Array Arguments ..
      double precision   a( lda, * ), x( * ), y( * )
C     ..
C
C  Purpose
C  =======
C
C  DGER   performs the rank 1 operation
C
C     A := alpha*x*y' + A,
C
C  where alpha is a scalar, x is an m element vector, y is an n element
C  vector and A is an m by n matrix.
C
C  Parameters
C  ==========
C
C  M      - INTEGER.
C           On entry, M specifies the number of rows of the matrix A.
C           M must be at least zero.
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the number of columns of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  ALPHA  - DOUBLE PRECISION.
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  X      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( m - 1 )*abs( INCX ) ).
C           Before entry, the incremented array X must contain the m
C           element vector x.
C           Unchanged on exit.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C  Y      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCY ) ).
C           Before entry, the incremented array Y must contain the n
C           element vector y.
C           Unchanged on exit.
C
C  INCY   - INTEGER.
C           On entry, INCY specifies the increment for the elements of
C           Y. INCY must not be zero.
C           Unchanged on exit.
C
C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
C           Before entry, the leading m by n part of the array A must
C           contain the matrix of coefficients. On exit, A is
C           overwritten by the updated matrix.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, m ).
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C
C
C     .. Parameters ..
      double precision   zero
      parameter        ( zero = 0.0d+0 )
C     .. Local Scalars ..
      double precision   temp
      integer            i, info, ix, j, jy, kx
C     .. External Subroutines ..
      external           f06aaz
C     .. Intrinsic Functions ..
      intrinsic          max
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      info = 0
      if     ( m.lt.0 )then
         info = 1
      else if( n.lt.0 )then
         info = 2
      else if( incx.eq.0 )then
         info = 5
      else if( incy.eq.0 )then
         info = 7
      else if( lda.lt.max( 1, m ) )then
         info = 9
      end if
      if( info.ne.0 )then
         call f06aaz( 'dger/f06pmf  ', info )
         return
      end if
C
C     Quick return if possible.
C
      if( ( m.eq.0 ).or.( n.eq.0 ).or.( alpha.eq.zero ) )
     $   return
C
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through A.
C
      if( incy.gt.0 )then
         jy = 1
      else
         jy = 1 - ( n - 1 )*incy
      end if
      if( incx.eq.1 )then
         do 20, j = 1, n
            if( y( jy ).ne.zero )then
               temp = alpha*y( jy )
               do 10, i = 1, m
                  a( i, j ) = a( i, j ) + x( i )*temp
   10          continue
            end if
            jy = jy + incy
   20    continue
      else
         if( incx.gt.0 )then
            kx = 1
         else
            kx = 1 - ( m - 1 )*incx
         end if
         do 40, j = 1, n
            if( y( jy ).ne.zero )then
               temp = alpha*y( jy )
               ix   = kx
               do 30, i = 1, m
                  a( i, j ) = a( i, j ) + x( ix )*temp
                  ix        = ix        + incx
   30          continue
            end if
            jy = jy + incy
   40    continue
      end if
C
C     end of dger (f06pmf)
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine dsymv ( uplo, n, alpha, a, lda, x, incx,
     $                   beta, y, incy )

C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     Modified by PEG 8/31/94.
C     .. Entry Points ..
C     ENTRY      F06PCF( UPLO, N, ALPHA, A, LDA, X, INCX,
C    $                   BETA, Y, INCY )

C     .. Scalar Arguments ..
      double precision   alpha, beta
      integer            incx, incy, lda, n
      character*1        uplo
C     .. Array Arguments ..
      double precision   a( lda, * ), x( * ), y( * )
C     ..
C
C  Purpose
C  =======
C
C  DSYMV  performs the matrix-vector  operation
C
C     y := alpha*A*x + beta*y,
C
C  where alpha and beta are scalars, x and y are n element vectors and
C  A is an n by n symmetric matrix.
C
C  Parameters
C  ==========
C
C  UPLO   - CHARACTER*1.
C           On entry, UPLO specifies whether the upper or lower
C           triangular part of the array A is to be referenced as
C           follows:
C
C              UPLO = 'U' or 'u'   Only the upper triangular part of A
C                                  is to be referenced.
C
C              UPLO = 'L' or 'l'   Only the lower triangular part of A
C                                  is to be referenced.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the order of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  ALPHA  - DOUBLE PRECISION.
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
C           Before entry with  UPLO = 'U' or 'u', the leading n by n
C           upper triangular part of the array A must contain the upper
C           triangular part of the symmetric matrix and the strictly
C           lower triangular part of A is not referenced.
C           Before entry with UPLO = 'L' or 'l', the leading n by n
C           lower triangular part of the array A must contain the lower
C           triangular part of the symmetric matrix and the strictly
C           upper triangular part of A is not referenced.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, n ).
C           Unchanged on exit.
C
C  X      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCX ) ).
C           Before entry, the incremented array X must contain the n
C           element vector x.
C           Unchanged on exit.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C  BETA   - DOUBLE PRECISION.
C           On entry, BETA specifies the scalar beta. When BETA is
C           supplied as zero then Y need not be set on input.
C           Unchanged on exit.
C
C  Y      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCY ) ).
C           Before entry, the incremented array Y must contain the n
C           element vector y. On exit, Y is overwritten by the updated
C           vector y.
C
C  INCY   - INTEGER.
C           On entry, INCY specifies the increment for the elements of
C           Y. INCY must not be zero.
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C
C
C     .. Parameters ..
      double precision   one         , zero
      parameter        ( one = 1.0d+0, zero = 0.0d+0 )
C     .. Local Scalars ..
      double precision   temp1, temp2
      integer            i, info, ix, iy, j, jx, jy, kx, ky
C     .. External Subroutines ..
      external           f06aaz
C     .. Intrinsic Functions ..
      intrinsic          max
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      info = 0
      if     ( .not.(uplo.eq.'U' .or. uplo.eq.'u').and.
     $         .not.(uplo.eq.'L' .or. uplo.eq.'l')      )then
         info = 1
      else if( n.lt.0 )then
         info = 2
      else if( lda.lt.max( 1, n ) )then
         info = 5
      else if( incx.eq.0 )then
         info = 7
      else if( incy.eq.0 )then
         info = 10
      end if
      if( info.ne.0 )then
         call f06aaz( 'dsymv/f06pcf ', info )
         return
      end if
C
C     Quick return if possible.
C
      if( ( n.eq.0 ).or.( ( alpha.eq.zero ).and.( beta.eq.one ) ) )
     $   return
C
C     Set up the start points in  X  and  Y.
C
      if( incx.gt.0 )then
         kx = 1
      else
         kx = 1 - ( n - 1 )*incx
      end if
      if( incy.gt.0 )then
         ky = 1
      else
         ky = 1 - ( n - 1 )*incy
      end if
C
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through the triangular part
C     of A.
C
C     First form  y := beta*y.
C
      if( beta.ne.one )then
         if( incy.eq.1 )then
            if( beta.eq.zero )then
               do 10, i = 1, n
                  y( i ) = zero
   10          continue
            else
               do 20, i = 1, n
                  y( i ) = beta*y( i )
   20          continue
            end if
         else
            iy = ky
            if( beta.eq.zero )then
               do 30, i = 1, n
                  y( iy ) = zero
                  iy      = iy   + incy
   30          continue
            else
               do 40, i = 1, n
                  y( iy ) = beta*y( iy )
                  iy      = iy           + incy
   40          continue
            end if
         end if
      end if
      if( alpha.eq.zero )
     $   return
      if( (uplo.eq.'U' .or. uplo.eq.'u') )then
C
C        Form  y  when A is stored in upper triangle.
C
         if( ( incx.eq.1 ).and.( incy.eq.1 ) )then
            do 60, j = 1, n
               temp1 = alpha*x( j )
               temp2 = zero
               do 50, i = 1, j - 1
                  y( i ) = y( i ) + temp1*a( i, j )
                  temp2  = temp2  + a( i, j )*x( i )
   50          continue
               y( j ) = y( j ) + temp1*a( j, j ) + alpha*temp2
   60       continue
         else
            jx = kx
            jy = ky
            do 80, j = 1, n
               temp1 = alpha*x( jx )
               temp2 = zero
               ix    = kx
               iy    = ky
               do 70, i = 1, j - 1
                  y( iy ) = y( iy ) + temp1*a( i, j )
                  temp2   = temp2   + a( i, j )*x( ix )
                  ix      = ix      + incx
                  iy      = iy      + incy
   70          continue
               y( jy ) = y( jy ) + temp1*a( j, j ) + alpha*temp2
               jx      = jx      + incx
               jy      = jy      + incy
   80       continue
         end if
      else
C
C        Form  y  when A is stored in lower triangle.
C
         if( ( incx.eq.1 ).and.( incy.eq.1 ) )then
            do 100, j = 1, n
               temp1  = alpha*x( j )
               temp2  = zero
               y( j ) = y( j )       + temp1*a( j, j )
               do 90, i = j + 1, n
                  y( i ) = y( i ) + temp1*a( i, j )
                  temp2  = temp2  + a( i, j )*x( i )
   90          continue
               y( j ) = y( j ) + alpha*temp2
  100       continue
         else
            jx = kx
            jy = ky
            do 120, j = 1, n
               temp1   = alpha*x( jx )
               temp2   = zero
               y( jy ) = y( jy )       + temp1*a( j, j )
               ix      = jx
               iy      = jy
               do 110, i = j + 1, n
                  ix      = ix      + incx
                  iy      = iy      + incy
                  y( iy ) = y( iy ) + temp1*a( i, j )
                  temp2   = temp2   + a( i, j )*x( ix )
  110          continue
               y( jy ) = y( jy ) + alpha*temp2
               jx      = jx      + incx
               jy      = jy      + incy
  120       continue
         end if
      end if
C
C     end of dsymv (f06pcf)
C
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine dsyr  ( uplo, n, alpha, x, incx, a, lda )

C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     Modified by PEG 8/31/94.
C     .. Entry Points ..
C     ENTRY      F06PPF( UPLO, N, ALPHA, X, INCX, A, LDA )

C     .. Scalar Arguments ..
      double precision   alpha
      integer            incx, lda, n
      character*1        uplo
C     .. Array Arguments ..
      double precision   a( lda, * ), x( * )
C     ..
C
C  Purpose
C  =======
C
C  DSYR   performs the symmetric rank 1 operation
C
C     A := alpha*x*x' + A,
C
C  where alpha is a real scalar, x is an n element vector and A is an
C  n by n symmetric matrix.
C
C  Parameters
C  ==========
C
C  UPLO   - CHARACTER*1.
C           On entry, UPLO specifies whether the upper or lower
C           triangular part of the array A is to be referenced as
C           follows:
C
C              UPLO = 'U' or 'u'   Only the upper triangular part of A
C                                  is to be referenced.
C
C              UPLO = 'L' or 'l'   Only the lower triangular part of A
C                                  is to be referenced.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the order of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  ALPHA  - DOUBLE PRECISION.
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  X      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCX ) ).
C           Before entry, the incremented array X must contain the n
C           element vector x.
C           Unchanged on exit.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
C           Before entry with  UPLO = 'U' or 'u', the leading n by n
C           upper triangular part of the array A must contain the upper
C           triangular part of the symmetric matrix and the strictly
C           lower triangular part of A is not referenced. On exit, the
C           upper triangular part of the array A is overwritten by the
C           upper triangular part of the updated matrix.
C           Before entry with UPLO = 'L' or 'l', the leading n by n
C           lower triangular part of the array A must contain the lower
C           triangular part of the symmetric matrix and the strictly
C           upper triangular part of A is not referenced. On exit, the
C           lower triangular part of the array A is overwritten by the
C           lower triangular part of the updated matrix.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, n ).
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C
C
C     .. Parameters ..
      double precision   zero
      parameter        ( zero = 0.0d+0 )
C     .. Local Scalars ..
      double precision   temp
      integer            i, info, ix, j, jx, kx
C     .. External Subroutines ..
      external           f06aaz
C     .. Intrinsic Functions ..
      intrinsic          max
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      info = 0
      if     ( .not.(uplo.eq.'U' .or. uplo.eq.'u').and.
     $         .not.(uplo.eq.'L' .or. uplo.eq.'l')      )then
         info = 1
      else if( n.lt.0 )then
         info = 2
      else if( incx.eq.0 )then
         info = 5
      else if( lda.lt.max( 1, n ) )then
         info = 7
      end if
      if( info.ne.0 )then
         call f06aaz( 'dsyr/f06ppf  ', info )
         return
      end if
C
C     Quick return if possible.
C
      if( ( n.eq.0 ).or.( alpha.eq.zero ) )
     $   return
C
C     Set the start point in X if the increment is not unity.
C
      if( incx.le.0 )then
         kx = 1 - ( n - 1 )*incx
      else if( incx.ne.1 )then
         kx = 1
      end if
C
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through the triangular part
C     of A.
C
      if( (uplo.eq.'U' .or. uplo.eq.'u') )then
C
C        Form  A  when A is stored in upper triangle.
C
         if( incx.eq.1 )then
            do 20, j = 1, n
               if( x( j ).ne.zero )then
                  temp = alpha*x( j )
                  do 10, i = 1, j
                     a( i, j ) = a( i, j ) + x( i )*temp
   10             continue
               end if
   20       continue
         else
            jx = kx
            do 40, j = 1, n
               if( x( jx ).ne.zero )then
                  temp = alpha*x( jx )
                  ix   = kx
                  do 30, i = 1, j
                     a( i, j ) = a( i, j ) + x( ix )*temp
                     ix        = ix        + incx
   30             continue
               end if
               jx = jx + incx
   40       continue
         end if
      else
C
C        Form  A  when A is stored in lower triangle.
C
         if( incx.eq.1 )then
            do 60, j = 1, n
               if( x( j ).ne.zero )then
                  temp = alpha*x( j )
                  do 50, i = j, n
                     a( i, j ) = a( i, j ) + x( i )*temp
   50             continue
               end if
   60       continue
         else
            jx = kx
            do 80, j = 1, n
               if( x( jx ).ne.zero )then
                  temp = alpha*x( jx )
                  ix   = jx
                  do 70, i = j, n
                     a( i, j ) = a( i, j ) + x( ix )*temp
                     ix        = ix        + incx
   70             continue
               end if
               jx = jx + incx
   80       continue
         end if
      end if
C
C     end of dsyr (f06ppf)
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine dtrmv ( uplo, trans, diag, n, a, lda, x, incx )

C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     Modified by PEG 8/31/94.
C     .. Entry Points ..
C     ENTRY      F06PFF( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )

C     .. Scalar Arguments ..
      integer            incx, lda, n
      character*1        diag, trans, uplo
C     .. Array Arguments ..
      double precision   a( lda, * ), x( * )
C     ..
C
C  Purpose
C  =======
C
C  DTRMV  performs one of the matrix-vector operations
C
C     x := A*x,   or   x := A'*x,
C
C  where x is n element vector and A is an n by n unit, or non-unit,
C  upper or lower triangular matrix.
C
C  Parameters
C  ==========
C
C  UPLO   - CHARACTER*1.
C           On entry, UPLO specifies whether the matrix is an upper or
C           lower triangular matrix as follows:
C
C              UPLO = 'U' or 'u'   A is an upper triangular matrix.
C
C              UPLO = 'L' or 'l'   A is a lower triangular matrix.
C
C           Unchanged on exit.
C
C  TRANS  - CHARACTER*1.
C           On entry, TRANS specifies the operation to be performed as
C           follows:
C
C              TRANS = 'N' or 'n'   x := A*x.
C
C              TRANS = 'T' or 't'   x := A'*x.
C
C              TRANS = 'C' or 'c'   x := A'*x.
C
C           Unchanged on exit.
C
C  DIAG   - CHARACTER*1.
C           On entry, DIAG specifies whether or not A is unit
C           triangular as follows:
C
C              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
C
C              DIAG = 'N' or 'n'   A is not assumed to be unit
C                                  triangular.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the order of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
C           Before entry with  UPLO = 'U' or 'u', the leading n by n
C           upper triangular part of the array A must contain the upper
C           triangular matrix and the strictly lower triangular part of
C           A is not referenced.
C           Before entry with UPLO = 'L' or 'l', the leading n by n
C           lower triangular part of the array A must contain the lower
C           triangular matrix and the strictly upper triangular part of
C           A is not referenced.
C           Note that when  DIAG = 'U' or 'u', the diagonal elements of
C           A are not referenced either, but are assumed to be unity.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, n ).
C           Unchanged on exit.
C
C  X      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCX ) ).
C           Before entry, the incremented array X must contain the n
C           element vector x. On exit, X is overwritten with the
C           tranformed vector x.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C
C
C     .. Parameters ..
      double precision   zero
      parameter        ( zero = 0.0d+0 )
C     .. Local Scalars ..
      double precision   temp
      integer            i, info, ix, j, jx, kx
      logical            nounit
C     .. External Subroutines ..
      external           f06aaz
C     .. Intrinsic Functions ..
      intrinsic          max
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      info = 0
      if     ( .not.(uplo .eq.'U' .or. uplo .eq.'u').and.
     $         .not.(uplo .eq.'L' .or. uplo .eq.'l')      )then
         info = 1
      else if( .not.(trans.eq.'N' .or. trans.eq.'n').and.
     $         .not.(trans.eq.'T' .or. trans.eq.'t').and.
     $         .not.(trans.eq.'C' .or. trans.eq.'c')      )then
         info = 2
      else if( .not.(diag .eq.'U' .or. diag .eq.'u').and.
     $         .not.(diag .eq.'N' .or. diag .eq.'n')      )then
         info = 3
      else if( n.lt.0 )then
         info = 4
      else if( lda.lt.max( 1, n ) )then
         info = 6
      else if( incx.eq.0 )then
         info = 8
      end if
      if( info.ne.0 )then
         call f06aaz( 'dtrmv/f06pff ', info )
         return
      end if
C
C     Quick return if possible.
C
      if( n.eq.0 )
     $   return
C
      nounit = (diag.eq.'N' .or. diag.eq.'n')
C
C     Set up the start point in X if the increment is not unity. This
C     will be  ( N - 1 )*INCX  too small for descending loops.
C
      if( incx.le.0 )then
         kx = 1 - ( n - 1 )*incx
      else if( incx.ne.1 )then
         kx = 1
      end if
C
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through A.
C
      if( (trans.eq.'N' .or. trans.eq.'n') )then
C
C        Form  x := A*x.
C
         if( (uplo.eq.'U' .or. uplo.eq.'u') )then
            if( incx.eq.1 )then
               do 20, j = 1, n
                  if( x( j ).ne.zero )then
                     temp = x( j )
                     do 10, i = 1, j - 1
                        x( i ) = x( i ) + temp*a( i, j )
   10                continue
                     if( nounit )
     $                  x( j ) = x( j )*a( j, j )
                  end if
   20          continue
            else
               jx = kx
               do 40, j = 1, n
                  if( x( jx ).ne.zero )then
                     temp = x( jx )
                     ix   = kx
                     do 30, i = 1, j - 1
                        x( ix ) = x( ix ) + temp*a( i, j )
                        ix      = ix      + incx
   30                continue
                     if( nounit )
     $                  x( jx ) = x( jx )*a( j, j )
                  end if
                  jx = jx + incx
   40          continue
            end if
         else
            if( incx.eq.1 )then
               do 60, j = n, 1, -1
                  if( x( j ).ne.zero )then
                     temp = x( j )
                     do 50, i = n, j + 1, -1
                        x( i ) = x( i ) + temp*a( i, j )
   50                continue
                     if( nounit )
     $                  x( j ) = x( j )*a( j, j )
                  end if
   60          continue
            else
               kx = kx + ( n - 1 )*incx
               jx = kx
               do 80, j = n, 1, -1
                  if( x( jx ).ne.zero )then
                     temp = x( jx )
                     ix   = kx
                     do 70, i = n, j + 1, -1
                        x( ix ) = x( ix ) + temp*a( i, j )
                        ix      = ix      - incx
   70                continue
                     if( nounit )
     $                  x( jx ) = x( jx )*a( j, j )
                  end if
                  jx = jx - incx
   80          continue
            end if
         end if
      else
C
C        Form  x := A'*x.
C
         if( (uplo.eq.'U' .or. uplo.eq.'u') )then
            if( incx.eq.1 )then
               do 100, j = n, 1, -1
                  temp = x( j )
                  if( nounit )
     $               temp = temp*a( j, j )
                  do 90, i = j - 1, 1, -1
                     temp = temp + a( i, j )*x( i )
   90             continue
                  x( j ) = temp
  100          continue
            else
               jx = kx + ( n - 1 )*incx
               do 120, j = n, 1, -1
                  temp = x( jx )
                  ix   = jx
                  if( nounit )
     $               temp = temp*a( j, j )
                  do 110, i = j - 1, 1, -1
                     ix   = ix   - incx
                     temp = temp + a( i, j )*x( ix )
  110             continue
                  x( jx ) = temp
                  jx      = jx   - incx
  120          continue
            end if
         else
            if( incx.eq.1 )then
               do 140, j = 1, n
                  temp = x( j )
                  if( nounit )
     $               temp = temp*a( j, j )
                  do 130, i = j + 1, n
                     temp = temp + a( i, j )*x( i )
  130             continue
                  x( j ) = temp
  140          continue
            else
               jx = kx
               do 160, j = 1, n
                  temp = x( jx )
                  ix   = jx
                  if( nounit )
     $               temp = temp*a( j, j )
                  do 150, i = j + 1, n
                     ix   = ix   + incx
                     temp = temp + a( i, j )*x( ix )
  150             continue
                  x( jx ) = temp
                  jx      = jx   + incx
  160          continue
            end if
         end if
      end if
C
C     end of dtrmv (f06pff)
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine dtrsv ( uplo, trans, diag, n, a, lda, x, incx )

C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     Modified by PEG 8/31/94.
C     .. Entry Points ..
C     ENTRY      F06PJF( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )

C     .. Scalar Arguments ..
      integer            incx, lda, n
      character*1        diag, trans, uplo
C     .. Array Arguments ..
      double precision   a( lda, * ), x( * )
C     ..
C
C  Purpose
C  =======
C
C  DTRSV  solves one of the systems of equations
C
C     A*x = b,   or   A'*x = b,
C
C  where b and x are n element vectors and A is an n by n unit, or
C  non-unit, upper or lower triangular matrix.
C
C  No test for singularity or near-singularity is included in this
C  routine. Such tests must be performed before calling this routine.
C
C  Parameters
C  ==========
C
C  UPLO   - CHARACTER*1.
C           On entry, UPLO specifies whether the matrix is an upper or
C           lower triangular matrix as follows:
C
C              UPLO = 'U' or 'u'   A is an upper triangular matrix.
C
C              UPLO = 'L' or 'l'   A is a lower triangular matrix.
C
C           Unchanged on exit.
C
C  TRANS  - CHARACTER*1.
C           On entry, TRANS specifies the equations to be solved as
C           follows:
C
C              TRANS = 'N' or 'n'   A*x = b.
C
C              TRANS = 'T' or 't'   A'*x = b.
C
C              TRANS = 'C' or 'c'   A'*x = b.
C
C           Unchanged on exit.
C
C  DIAG   - CHARACTER*1.
C           On entry, DIAG specifies whether or not A is unit
C           triangular as follows:
C
C              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
C
C              DIAG = 'N' or 'n'   A is not assumed to be unit
C                                  triangular.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the order of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
C           Before entry with  UPLO = 'U' or 'u', the leading n by n
C           upper triangular part of the array A must contain the upper
C           triangular matrix and the strictly lower triangular part of
C           A is not referenced.
C           Before entry with UPLO = 'L' or 'l', the leading n by n
C           lower triangular part of the array A must contain the lower
C           triangular matrix and the strictly upper triangular part of
C           A is not referenced.
C           Note that when  DIAG = 'U' or 'u', the diagonal elements of
C           A are not referenced either, but are assumed to be unity.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, n ).
C           Unchanged on exit.
C
C  X      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCX ) ).
C           Before entry, the incremented array X must contain the n
C           element right-hand side vector b. On exit, X is overwritten
C           with the solution vector x.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C
C
C     .. Parameters ..
      double precision   zero
      parameter        ( zero = 0.0d+0 )
C     .. Local Scalars ..
      double precision   temp
      integer            i, info, ix, j, jx, kx
      logical            nounit
C     .. External Subroutines ..
      external           f06aaz
C     .. Intrinsic Functions ..
      intrinsic          max
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      info = 0
      if     ( .not.(uplo .eq.'U' .or. uplo .eq.'u').and.
     $         .not.(uplo .eq.'L' .or. uplo .eq.'l')      )then
         info = 1
      else if( .not.(trans.eq.'N' .or. trans.eq.'n').and.
     $         .not.(trans.eq.'T' .or. trans.eq.'t').and.
     $         .not.(trans.eq.'C' .or. trans.eq.'c')      )then
         info = 2
      else if( .not.(diag .eq.'U' .or. diag .eq.'u').and.
     $         .not.(diag .eq.'N' .or. diag .eq.'n')      )then
         info = 3
      else if( n.lt.0 )then
         info = 4
      else if( lda.lt.max( 1, n ) )then
         info = 6
      else if( incx.eq.0 )then
         info = 8
      end if
      if( info.ne.0 )then
         call f06aaz( 'dtrsv/f06pjf ', info )
         return
      end if
C
C     Quick return if possible.
C
      if( n.eq.0 )
     $   return
C
      nounit = (diag.eq.'N' .or. diag.eq.'n')
C
C     Set up the start point in X if the increment is not unity. This
C     will be  ( N - 1 )*INCX  too small for descending loops.
C
      if( incx.le.0 )then
         kx = 1 - ( n - 1 )*incx
      else if( incx.ne.1 )then
         kx = 1
      end if
C
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through A.
C
      if( (trans.eq.'N' .or. trans.eq.'n') )then
C
C        Form  x := inv( A )*x.
C
         if( (uplo.eq.'U' .or. uplo.eq.'u') )then
            if( incx.eq.1 )then
               do 20, j = n, 1, -1
                  if( x( j ).ne.zero )then
                     if( nounit )
     $                  x( j ) = x( j )/a( j, j )
                     temp = x( j )
                     do 10, i = j - 1, 1, -1
                        x( i ) = x( i ) - temp*a( i, j )
   10                continue
                  end if
   20          continue
            else
               jx = kx + ( n - 1 )*incx
               do 40, j = n, 1, -1
                  if( x( jx ).ne.zero )then
                     if( nounit )
     $                  x( jx ) = x( jx )/a( j, j )
                     temp = x( jx )
                     ix   = jx
                     do 30, i = j - 1, 1, -1
                        ix      = ix      - incx
                        x( ix ) = x( ix ) - temp*a( i, j )
   30                continue
                  end if
                  jx = jx - incx
   40          continue
            end if
         else
            if( incx.eq.1 )then
               do 60, j = 1, n
                  if( x( j ).ne.zero )then
                     if( nounit )
     $                  x( j ) = x( j )/a( j, j )
                     temp = x( j )
                     do 50, i = j + 1, n
                        x( i ) = x( i ) - temp*a( i, j )
   50                continue
                  end if
   60          continue
            else
               jx = kx
               do 80, j = 1, n
                  if( x( jx ).ne.zero )then
                     if( nounit )
     $                  x( jx ) = x( jx )/a( j, j )
                     temp = x( jx )
                     ix   = jx
                     do 70, i = j + 1, n
                        ix      = ix      + incx
                        x( ix ) = x( ix ) - temp*a( i, j )
   70                continue
                  end if
                  jx = jx + incx
   80          continue
            end if
         end if
      else
C
C        Form  x := inv( A' )*x.
C
         if( (uplo.eq.'U' .or. uplo.eq.'u') )then
            if( incx.eq.1 )then
               do 100, j = 1, n
                  temp = x( j )
                  do 90, i = 1, j - 1
                     temp = temp - a( i, j )*x( i )
   90             continue
                  if( nounit )
     $               temp = temp/a( j, j )
                  x( j ) = temp
  100          continue
            else
               jx = kx
               do 120, j = 1, n
                  temp = x( jx )
                  ix   = kx
                  do 110, i = 1, j - 1
                     temp = temp - a( i, j )*x( ix )
                     ix   = ix   + incx
  110             continue
                  if( nounit )
     $               temp = temp/a( j, j )
                  x( jx ) = temp
                  jx      = jx   + incx
  120          continue
            end if
         else
            if( incx.eq.1 )then
               do 140, j = n, 1, -1
                  temp = x( j )
                  do 130, i = n, j + 1, -1
                     temp = temp - a( i, j )*x( i )
  130             continue
                  if( nounit )
     $               temp = temp/a( j, j )
                  x( j ) = temp
  140          continue
            else
               kx = kx + ( n - 1 )*incx
               jx = kx
               do 160, j = n, 1, -1
                  temp = x( jx )
                  ix   = kx
                  do 150, i = n, j + 1, -1
                     temp = temp - a( i, j )*x( ix )
                     ix   = ix   - incx
  150             continue
                  if( nounit )
     $               temp = temp/a( j, j )
                  x( jx ) = temp
                  jx      = jx   - incx
  160          continue
            end if
         end if
      end if
C
C     end of dtrsv (f06pjf)
      end
