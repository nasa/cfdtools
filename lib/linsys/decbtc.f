C+----------------------------------------------------------------------
C
      SUBROUTINE DECBTC ( NDIM,NSIZE,NBLOKS, L,D,U, IP, NWRK, WRK, IER )
C
C  PURPOSE:  DECBTC may be used with SOLBTC to solve cyclic block-tri-
C            diagonal systems  Ax = b  of the form
C
C                              XX    X  |     |
C                              XXX      |     |
C                               XXX     |     |
C                                XXX    |  =  |
C                                 XXX   |     |
C                                  XXX  |     |
C                              X    XX  |     |
C
C            for one or more right hand sides. DECBTC should be called
C            once for a given left hand side, to do the factorization,
C            then SOLBTC should be called as many times as  there  are
C            right hand sides.
C
C            All blocks are assumed to be of the same size, NSIZE, for
C            any given system, but provision is made for handling sys-
C            tems with differing block sizes using the same workspace.
C
C  METHOD:   Partition matrix A to work with true block-tridiagonal A1.
C
C                                  |
C                             A1   | V  y     r
C                                  |       =
C                            ------+--  -     -
C                             U    |A2  z     s
C
C            Each of U and V shown here has just 2 non-zero blocks.
C
C            (1) Factorize  A1  in the usual way
C            (2) Solve    A1 W  =  V
C            (3) Form        C  =  A2 - U W  and factorize it
C            (4) Solve    A1 t  =  r
C            (5) Solve     C z  =  s  - U t
C            (6) Solve    A1 y  =  r  - V z
C
C            Steps (1):(3) are performed by DECBTC.
C            Steps (4):(6) are performed by SOLBTC for each RHS.
C
C  PARAMETERS:
C  VARIABLE  TYPE  I/O/S  DIMENSIONS  DESCRIPTION
C
C    NDIM     INT    I       -        Declared row dimension of all
C                                     multi-dimensional arrays being
C                                     passed.  (For 3-D arrays, the
C                                     first two dimensions must be
C                                     declared the same in the calling
C                                     program.)  NDIM >= NSIZE.  This
C                                     permits use of the same arrays for
C                                     systems with differing block size.
C    NSIZE    INT    I       -        Size of each block in the system
C                                     being solved.
C    NBLOKS   INT    I       -        Number of blocks along the main
C                                     diagonal of matrix A.
C    L        REAL  I/O   NDIM,NDIM,  Input with the NSIZE*NSIZE lower
C                         NBLOKS      diagonal blocks in L(,,K) for
C                                     K=2:NBLOKS, and the block in the
C                                     (1,NBLOKS) position in L(,,1).
C                                     Output along with arrays D and U
C                                     with the necessary factorizations.
C    D        REAL  I/O   NDIM,NDIM,  Input with the diagonal blocks.
C                         NBLOKS      D(I,J,K) contains element (I,J)
C                                     of the Kth main diagonal block.
C    U        REAL  I/O   NDIM,NDIM,  Input with the upper diagonal
C                         NBLOKS      blocks in U(,,K) for K=1:NBLOKS-1
C                                     and the block in the (NBLOKS,1)
C                                     position in U(,,NBLOKS).
C    IP       INT    O    NDIM,NBLOKS Block vector for pivot informa-
C                                     tion.  Kept separate from WRK(*)
C                                     because IP(NSIZE,NBLOKS) has to
C                                     be checked here...
C    NWRK     INT    I       -        Amount of workspace provided.
C                                     NWRK >= NDIM*( 2*NDIM+NBLOKS-1 )
C    WRK      REAL  S/O   NWRK        Work space for DECBT/SOLBT, etc.
C    IER      INT    O       -        Error return code.
C                                     IER = 0 means no problems;
C                                         =-1 means NSIZE or NBLOKS bad;
C                                         =-2 means NWRK is too small;
C                                         =-3 means matrix C singular;
C                                         = K means the (modified) Kth
C                                             diagonal block singular.
C
C  ENVIRONMENT:  VAX/VMS, VAX-11 FORTRAN.  (Generic SQRT and ABS
C                are the only FORTRAN-77 features taken advantage of.)
C
C  EXTERNAL REFERENCES:
C    NAME    DESCRIPTION
C   DECBT    Factorizes block-tridiagonal matrix (A1 here).
C   SOLBT    Uses factorization of DECBT for one or more rh sides.
C   DECOMP   Factorizes dense matrix (C here; also used by DECBT).
C
C  NOTES:
C       (1)  Direct factorization of matrix A is another possible way
C            to go, but it involves fill-in in the triangular factors
C            of the submatrices U, V indicated above....more storage.
C
C       (2)  The bulk of the work is still done in factorizing matrix
C            A1; the multiple solves involving A1, and the lesser one
C            with matrix C are not as expensive as they may seem.
C
C       (3)  Do not be misled by the dimensioning of WRK() below.  It
C            is treated as WRK(NDIM,NDIM,*) for convenience here  and
C            in SOLBTC.  This is due to the fact that DECBT/SOLBT ex-
C            pect 2 extra blocks which have to be supplied (as zeros)
C            in the space originally containing the nonzeros of V...
C
C  AUTHOR:   David Saunders, Informatics Inc, October 1981.
C
C-----------------------------------------------------------------------
C
      IMPLICIT REAL ( A-H, O-Z )
C
      REAL     L(NDIM,NDIM,NBLOKS), D(NDIM,NDIM,NBLOKS),
     +         U(NDIM,NDIM,NBLOKS), WRK(NDIM,NDIM,1)
      INTEGER  IP(NDIM,NBLOKS)
C
C
      IER = -2
      IF ( NWRK .LT. NDIM*( 2*NDIM + NBLOKS - 1 ) ) GO TO 999
C
      NBM1 = NBLOKS-1
      NBM2 = NBLOKS-2
C
C  *  Prepare for factoring the reduced, true block-tridiagonal matrix,
C  *  A1, by saving the non-zero blocks of submatrix V in the first and
C  *  second NDIM*NDIM segments of WRK(*).
C  *  Zero out the extra blocks expected by DECBT/SOLBT also.
C
      DO 220 J = 1, NSIZE
         DO 210 I = 1, NSIZE
            WRK(I,J,1) = L(I,J,1)
            L(I,J,1)   = 0.E+0
            WRK(I,J,2) = U(I,J,NBM1)
            U(I,J,NBM1)= 0.E+0
 210     CONTINUE
 220  CONTINUE
C
C  *  Factorize A1.
C
      CALL DECBT ( NSIZE, NDIM, NBM1, D, U, L, IP, IER )
      IF ( IER.NE.0 ) GO TO 999
C
C  *  Solve  A1 W = V  for 1 column of W at a time,
C  *  and form 1 column of  C = A2 - U W  as we go...
C
      IW1 = 2*NDIM*NDIM
      IW2 = IW1 + NBM2*NDIM
C
      DO 300 JCOL = 1, NSIZE
C
C  *     Set up the RHS block vector (only first, last blocks nonzero).
C
         DO 230 I = 1, NSIZE
            WRK(I+IW1,1,1) = WRK(I,JCOL,1)
            WRK(I+IW2,1,1) = WRK(I,JCOL,2)
 230     CONTINUE
C
         IW = IW1
         DO 250 J = 2, NBM2
            IW = IW + NDIM
            DO 240 I = 1, NSIZE
               WRK(I+IW,1,1) = 0.E+0
 240        CONTINUE
 250     CONTINUE
C
C  *     Solve for this column of W.
C
         CALL SOLBT ( NSIZE, NDIM, NBM1, D, U, L, WRK(1,1,3), IP )
C
C  *     Set up corresponding column of  C = A2 - U W  in place of A2.
C  *     (Only 2 blocks of U are nonzero.)
C
         DO 270 I = 1, NSIZE
            SUM1  = 0.E+0
            DO 260 J = 1, NSIZE
               SUM1  = U(I,J,NBLOKS) * WRK(J+IW1,1,1)  +
     +                 L(I,J,NBLOKS) * WRK(J+IW2,1,1)  +  SUM1
 260        CONTINUE
            D(I,JCOL,NBLOKS) = D(I,JCOL,NBLOKS) - SUM1
 270     CONTINUE
 300  CONTINUE
C
C  *  Factorize the Schur complement, A2, of A1 (in place).
C
      CALL DECOMP ( NSIZE, NDIM, D(1,1,NBLOKS), IP(1,NBLOKS) )
      IF ( IP(NSIZE,NBLOKS).NE.0 ) GO TO 999
C
      IER = -3
C
C
 999  RETURN
      END
