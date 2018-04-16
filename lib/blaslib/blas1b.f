*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     File  blas1.f
*     NAG versions of the Level 1  BLAS Vector routines.
*
*     daxpy           dcopy           ddot            dnrm2    
*     dscal           dswap           idamax          drot
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine daxpy ( n, alpha, x, incx, y, incy )

C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     Modified by PEG 9/25/88. 8/31/94.
C     .. Entry Points ..
C     ENTRY      F06ECF( N, ALPHA, X, INCX, Y, INCY )
C     .. Scalar Arguments ..
      double precision   alpha
      integer            incx, incy, n
C     .. Array Arguments ..
      double precision   x( * ), y( * )
C     ..
C
C  daxpy/F06ECF performs the operation
C
C     y := alpha*x + y
C
C
C  Nag Fortran 77 version of the Blas routine DAXPY.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 3-September-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      double precision   zero
      parameter        ( zero = 0.0d+0 )
C     .. Local Scalars ..
      integer            i, ix, iy
C     ..
C     .. Executable Statements ..
      if( n.gt.0 )then
         if( alpha.ne.zero )then
            if( ( incx.eq.incy ).and.( incx.gt.0 ) )then
               do 10, ix = 1, 1 + ( n - 1 )*incx, incx
                  y( ix ) = alpha*x( ix ) + y( ix )
   10          continue
            else
               if( incy.ge.0 )then
                  iy = 1
               else
                  iy = 1 - ( n - 1 )*incy
               end if
               if( incx.gt.0 )then
                  do 20, ix = 1, 1 + ( n - 1 )*incx, incx
                     y( iy ) = alpha*x( ix ) + y( iy )
                     iy      = iy            + incy
   20             continue
               else
                  ix = 1 - ( n - 1 )*incx
                  do 30, i = 1, n
                     y( iy ) = alpha*x( ix ) + y( iy )
                     ix      = ix            + incx
                     iy      = iy            + incy
   30             continue
               end if
            end if
         end if
      end if
C
C     end of daxpy (f06ecf)
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine dcopy ( n, x, incx, y, incy )

C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     Modified by PEG 8/31/94.
C     .. Entry Points ..
C     ENTRY      F06EFF( N, X, INCX, Y, INCY )
C     .. Scalar Arguments ..
      integer            incx, incy, n
C     .. Array Arguments ..
      double precision   x( * ), y( * )
C     ..
C
C  dcopy/F06EFF performs the operation
C
C     y := x
C
C
C  Nag Fortran 77 version of the Blas routine DCOPY.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 26-November-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      integer            i, ix, iy
C     ..
C     .. Executable Statements ..
      if( n.gt.0 )then
         if( ( incx.eq.incy ).and.( incy.gt.0 ) )then
            do 10, iy = 1, 1 + ( n - 1 )*incy, incy
               y( iy ) = x( iy )
   10       continue
         else
            if( incx.ge.0 )then
               ix = 1
            else
               ix = 1 - ( n - 1 )*incx
            end if
            if( incy.gt.0 )then
               do 20, iy = 1, 1 + ( n - 1 )*incy, incy
                  y( iy ) = x( ix )
                  ix      = ix      + incx
   20          continue
            else
               iy = 1 - ( n - 1 )*incy
               do 30, i = 1, n
                  y( iy ) = x( ix )
                  iy      = iy      + incy
                  ix      = ix      + incx
   30          continue
            end if
         end if
      end if
C
C     end of dcopy (f06eff)
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      double precision function ddot  ( n, x, incx, y, incy )

C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     Modified by PEG 8/31/94.
C     .. Entry Points ..
C     DOUBLE PRECISION          F06EAF
C     ENTRY                     F06EAF( N, X, INCX, Y, INCY )
C     .. Scalar Arguments ..
      integer                           incx, incy, n
C     .. Array Arguments ..
      double precision                  x( * ), y( * )
C     ..
C
C  ddot/F06EAF returns the value
C
C     F06EAF = x'y
C
C
C  Nag Fortran 77 version of the Blas routine DDOT.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 21-September-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      double precision      zero
      parameter           ( zero = 0.0d+0 )
C     .. Local Scalars ..
      double precision      sum
      integer               i, ix, iy
C     ..
C     .. Executable Statements ..
      sum = zero
      if( n.gt.0 )then
         if( ( incx.eq.incy ).and.( incx.gt.0 ) )then
            do 10, ix = 1, 1 + ( n - 1 )*incx, incx
               sum = sum + x( ix )*y( ix )
   10       continue
         else
            if( incy.ge.0 )then
               iy = 1
            else
               iy = 1 - ( n - 1 )*incy
            end if
            if( incx.gt.0 )then
               do 20, ix = 1, 1 + ( n - 1 )*incx, incx
                  sum = sum + x( ix )*y( iy )
                  iy  = iy  + incy
   20          continue
            else
               ix = 1 - ( n - 1 )*incx
               do 30, i = 1, n
                  sum = sum + x( ix )*y( iy )
                  ix  = ix  + incx
                  iy  = iy  + incy
   30          continue
            end if
         end if
      end if
C
      ddot   = sum
C     F06EAF = SUM
C
C     end of ddot (f06eaf)
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      double precision function dnrm2 ( n, x, incx )

C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     Modified by PEG 8/31/94.
C     .. Entry Points ..
C     DOUBLE PRECISION          F06EJF
C     ENTRY                     F06EJF( N, X, INCX )
C     .. Scalar Arguments ..
      integer                           incx, n
C     .. Array Arguments ..
      double precision                  x( * )
C     ..
C
C  dnrm2/F06EJF returns the euclidean norm of a vector via the function
C  name, so that
C
C     dnrm2/F06EJF := sqrt( x'*x )
C
C
C  Nag Fortran 77 version of the Blas routine DNRM2.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 25-October-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      double precision      one         , zero
      parameter           ( one = 1.0d+0, zero = 0.0d+0 )
C     .. Local Scalars ..
      double precision      norm, scale, ssq
C     .. External Functions ..
      double precision      f06bmf
      external              f06bmf
C     .. External Subroutines ..
      external              f06fjf
C     .. Intrinsic Functions ..
      intrinsic             abs
C     ..
C     .. Executable Statements ..
      if( n.lt.1 )then
         norm  = zero
      else if( n.eq.1 )then
         norm  = abs( x( 1 ) )
      else
         scale = zero
         ssq   = one
         call f06fjf( n, x, incx, scale, ssq )
         norm  = f06bmf( scale, ssq )
      end if
C
      dnrm2  = norm
C     F06EJF = NORM
C
C     end of dnrm2 (f06ejf)
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine dscal ( n, alpha, x, incx )

C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     Modified by PEG 8/31/94.
C     .. Entry Points ..
C     ENTRY      f06EDF( N, ALPHA, X, INCX )
C     .. Scalar Arguments ..
      double precision   alpha
      integer            incx, n
C     .. Array Arguments ..
      double precision   x( * )
C     ..
C
C  dscal/F06EDF performs the operation
C
C     x := alpha*x
C
C
C  Nag Fortran 77 version of the Blas routine DSCAL.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 26-November-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      double precision   one         , zero
      parameter        ( one = 1.0d+0, zero = 0.0d+0 )
C     .. Local Scalars ..
      integer            ix
C     ..
C     .. Executable Statements ..
      if( n.gt.0 )then
         if( alpha.eq.zero )then
            do 10, ix = 1, 1 + ( n - 1 )*incx, incx
               x( ix ) = zero
   10       continue
         else if( alpha.eq.( -one ) )then
            do 20, ix = 1, 1 + ( n - 1 )*incx, incx
               x( ix ) = -x( ix )
   20       continue
         else if( alpha.ne.one )then
            do 30, ix = 1, 1 + ( n - 1 )*incx, incx
               x( ix ) = alpha*x( ix )
   30       continue
         end if
      end if
C
C     end of dscal (f06edf)
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine dswap ( n, x, incx, y, incy )

C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     Modified by PEG 8/31/94.
C     .. Entry Points ..
C     ENTRY      F06EGF( N, X, INCX, Y, INCY )
C     .. Scalar Arguments ..
      integer            incx, incy, n
C     .. Array Arguments ..
      double precision   x( * ), y( * )
C     ..
C
C  dswap/F06EGF performs the operations
C
C     temp := x,   x := y,   y := temp.
C
C
C  Nag Fortran 77 version of the Blas routine DSWAP.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 26-November-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      double precision   temp
      integer            i, ix, iy
C     ..
C     .. Executable Statements ..
      if( n.gt.0 )then
         if( ( incx.eq.incy ).and.( incy.gt.0 ) )then
            do 10, iy = 1, 1 + ( n - 1 )*incy, incy
               temp    = x( iy )
               x( iy ) = y( iy )
               y( iy ) = temp
   10       continue
         else
            if( incx.ge.0 )then
               ix = 1
            else
               ix = 1 - ( n - 1 )*incx
            end if
            if( incy.gt.0 )then
               do 20, iy = 1, 1 + ( n - 1 )*incy, incy
                  temp    = x( ix )
                  x( ix ) = y( iy )
                  y( iy ) = temp
                  ix      = ix      + incx
   20          continue
            else
               iy = 1 - ( n - 1 )*incy
               do 30, i = 1, n
                  temp    = x( ix )
                  x( ix ) = y( iy )
                  y( iy ) = temp
                  iy      = iy      + incy
                  ix      = ix      + incx
   30          continue
            end if
         end if
      end if
C
C     end of dswap (f06egf)
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine drot  ( n, x, incx, y, incy, c, s )

C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     Modified by PEG 8/31/94.
C     .. Entry Points ..
C     ENTRY      F06EPF( N, X, INCX, Y, INCY, C, S )
C     .. Scalar Arguments ..
      double precision   c, s
      integer            incx, incy, n
C     .. Array Arguments ..
      double precision   x( * ), y( * )
C     ..
C
C  drot/F06EPF performs the plane rotation
C
C     ( x  y ) = ( x  y )*( c  -s ).
C                         ( s   c )
C
C
C  Nag Fortran 77 version of the Blas routine DROT.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 23-January-1984.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      double precision   one         , zero
      parameter        ( one = 1.0d+0, zero = 0.0d+0 )
C     .. Local Scalars ..
      double precision   temp1
      integer            i, ix, iy
C     ..
C     .. Executable Statements ..
      if( n.gt.0 )then
         if( ( s.ne.zero ).or.( c.ne.one ) )then
            if( ( c.eq.zero ).and.( s.eq.one ) )then
               if( ( incx.eq.incy ).and.( incx.gt.0 ) )then
                  do 10, ix = 1, 1 + ( n - 1 )*incx, incx
                     temp1   = -x( ix )
                     x( ix ) =  y( ix )
                     y( ix ) =  temp1
   10             continue
               else
                  if( incy.ge.0 )then
                     iy = 1
                  else
                     iy = 1 - ( n - 1 )*incy
                  end if
                  if( incx.gt.0 )then
                     do 20, ix = 1, 1 + ( n - 1 )*incx, incx
                        temp1   = -x( ix )
                        x( ix ) =  y( iy )
                        y( iy ) =  temp1
                        iy      =  iy       + incy
   20                continue
                  else
                     ix = 1 - ( n - 1 )*incx
                     do 30, i = 1, n
                        temp1   = -x( ix )
                        x( ix ) =  y( iy )
                        y( iy ) =  temp1
                        ix      =  ix      + incx
                        iy      =  iy      + incy
   30                continue
                  end if
               end if
            else if( ( c.eq.zero ).and.( s.eq.( -one ) ) )then
               if( ( incx.eq.incy ).and.( incx.gt.0 ) )then
                  do 40, ix = 1, 1 + ( n - 1 )*incx, incx
                     temp1   =  x( ix )
                     x( ix ) = -y( ix )
                     y( ix ) =  temp1
   40             continue
               else
                  if( incy.ge.0 )then
                     iy = 1
                  else
                     iy = 1 - ( n - 1 )*incy
                  end if
                  if( incx.gt.0 )then
                     do 50, ix = 1, 1 + ( n - 1 )*incx, incx
                        temp1   =  x( ix )
                        x( ix ) = -y( iy )
                        y( iy ) =  temp1
                        iy      =  iy       + incy
   50                continue
                  else
                     ix = 1 - ( n - 1 )*incx
                     do 60, i = 1, n
                        temp1   =  x( ix )
                        x( ix ) = -y( iy )
                        y( iy ) =  temp1
                        ix      =  ix      + incx
                        iy      =  iy      + incy
   60                continue
                  end if
               end if
            else
               if( ( incx.eq.incy ).and.( incx.gt.0 ) )then
                  do 70, ix = 1, 1 + ( n - 1 )*incx, incx
                     temp1   = x( ix )
                     x( ix ) = s*y( ix ) + c*temp1
                     y( ix ) = c*y( ix ) - s*temp1
   70             continue
               else
                  if( incy.ge.0 )then
                     iy = 1
                  else
                     iy = 1 - ( n - 1 )*incy
                  end if
                  if( incx.gt.0 )then
                     do 80, ix = 1, 1 + ( n - 1 )*incx, incx
                        temp1   = x( ix )
                        x( ix ) = s*y( iy ) + c*temp1
                        y( iy ) = c*y( iy ) - s*temp1
                        iy      = iy        + incy
   80                continue
                  else
                     ix = 1 - ( n - 1 )*incx
                     do 90, i = 1, n
                        temp1   = x( ix )
                        x( ix ) = s*y( iy ) + c*temp1
                        y( iy ) = c*y( iy ) - s*temp1
                        ix      = ix        + incx
                        iy      = iy        + incy
   90                continue
                  end if
               end if
            end if
         end if
      end if
C
C     end of drot (f06epf)
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      integer function idamax( n, x, incx )

C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     Modified by PEG 8/31/94.
C     .. Entry Points ..
C     INTEGER          F06JLF
C     ENTRY            F06JLF( N, X, INCX )
C     .. Scalar Arguments ..
      integer                  incx, n
C     .. Array Arguments ..
      double precision         x( * )
C     ..
C
C  F06JLF returns the smallest value of i such that
C
C     abs( x( i ) ) = max( abs( x( j ) ) )
C                      j
C
C
C  Nag Fortran 77 version of the Blas routine IDAMAX.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 31-May-1983.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      double precision         xmax
      integer                  i, imax, ix
C     .. Intrinsic Functions ..
      intrinsic                abs
C     ..
C     .. Executable Statements ..
      if( n.gt.0 )then
         imax = 1
         if( n.gt.1 )then
            xmax = abs( x( 1 ) )
            ix   = 1
            do 10, i = 2, n
               ix = ix + incx
               if( xmax.lt.abs( x( ix ) ) )then
                  xmax = abs( x( ix ) )
                  imax = i
               end if
   10       continue
         end if
      else
         imax = 0
      end if
C
      idamax = imax
C     F06JLF = IMAX
C
C     end of idamax (f06jlf)
      end

