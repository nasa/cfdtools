C+----------------------------------------------------------------------
C
      SUBROUTINE COLSLV ( NDIM,NDIM2,NROWS,NBLOKS,A,RHS,TOL,PMIN,ISING )
C
C     PURPOSE:
C       General solver for the linear systems arising from collocation
C       techniques. These can be thought of as bidiagonal systems, and
C       will always have an even number of equations in each block.
C       They can be illustrated by the following block array  A(2,4,4)
C       with right hand side array  RHS(2,4):
C
C         0000                    0
C         00XX                    X
C           XXXX                  X
C           XXXX                  X
C             XXXX                X
C             XXXX                X
C               XX00              X
C               0000              0
C
C     PARAMETERS:
C     VAR       TYPE  I/O/S  DIMENSION  DESCRIPTION
C
C     NDIM      INT     I      -        Declared row dimension of block
C                                       array A. NDIM>=NROWS (Permits
C                                       use of same arrays for systems
C                                       with different block sizes)
C     NDIM2     INT     I      -        Declared column dimension of
C                                       block array A. NDIM2>=2*NROWS
C     NROWS     INT     I      -        Number of rows in each block
C                                       of the system being solved.
C                                       NROWS<=NDIM<=NDIM2/2
C     NBLOKS    INT     I      -        Number of blocks in the system
C                                       of the type illustrated above.
C     A         REAL   I/S   NDIM,      Input with coefficients of
C                            NDIM2,     the block-bidiagonal matrix
C                            NBLOKS     illustrated. The zeroes shown do
C                                       not need to be entered.
C     RHS       REAL   I/O   NDIM,      Input with the initial right-
C                            NBLOKS     hand side values. Output with 
C                                       the final solution. Note the
C                                       unused portions of the first
C                                       and last blocks.
C     TOL       REAL    I      -        Minimum acceptable pivot size.
C                                       Enter 0. if no better tolerance
C                                       is known.
C     PMIN      REAL    O      -        Minimum pivot found during  
C                                       row interchanges.
C     ISING     INT     O      -        Positive K indicates block K is
C                                       singular. Negative K shows that
C                                       block K contains minimum pivot 
C                                       element.
C     EXTERNAL REFERENCES: None
C
C  ENVIRONMENT:  VAX/VMS, VAX-11 FORTRAN.  (Generic SQRT and ABS
C                are the only FORTRAN-77 features taken advantage of.)
C
C     ALGORITHM:  P.M. Prenter, Colorado State University
C     CODING:     Rosalie Lefkowitz, Informatics Inc., May 1981.
C
C     NOTES:
C     1: The user must initialize NDIM,NDIM2,NROWS and NBLOKS as
C        well as set up block arrays A and RHS. COLSLV uses only
C        the lower right quadrant of A(*,*,1) and the upper left
C        quadrant of A(*,*,NBLOKS).   The upper half of RHS(*,1)
C        and the lower half of RHS(*,NBLOKS) are unused.
C     2: COLSLV overwrites the right hand side with the solution
C        and also overwrites all entries of A() during the modi-
C        fied Gauss elimination process, that permits row inter-
C        changes not leading to fill-in of the original  matrix.
C     3: The calling program must check the error flag ISING  on
C        return from COLSLV.  A positive integer value points to
C        a singular block or one with pivot element less than or
C        equal to TOL.  ISING < 0 points to the block containing 
C        the smallest pivot element, in case this is of interest
C        to the user.
C     4: Matrix singularity in the first block will often be due
C        to the boundary conditions in collocation applications.
C        This may be rectifiable by interchanging columns in the
C        calling program,  with appropriate variable reordering.
C     5: Note that this version operates on the  right-hand-side
C        during the triangularization.   It is not intended  for
C        efficient handling of multiple right-hand-sides  for  a
C        given left-hand-side.
C
C----------------------------------------------------------------------
C
      IMPLICIT REAL ( A-H, O-Z )
C
      DIMENSION   A(NDIM,NDIM2,NBLOKS), RHS(NDIM,NBLOKS)
C
      IF ( TOL.LT.0.E+0 ) TOL = 0.E+0
      PMIN  = 1.0E+35
C
C ... Set row and column identifying indices
      M0 = NROWS/2
      M1 = M0 + 1
      M2 = NROWS + 1
      M4 = NROWS + M0
      M5 = NROWS + M1
      NCOLS  = NROWS + NROWS
      NBLKM1 = NBLOKS - 1
C
C ... Forward Sweep: triangularization with partial pivoting
      DO 220 K=1,NBLOKS
C ...    Operate on sub-block defined by columns M1-M4
         JEND = M4
         IF ( K.EQ.NBLOKS ) JEND = NROWS
         JENDM1 = JEND - 1
         IF ( NROWS.GT.2 ) GO TO 30
C ...    NROWS=2 is a special case
         IF ( K.EQ.1 .OR. K.EQ.NBLOKS ) GO TO 145
 30      JSTART = M1
         IF ( K.EQ.1 ) JSTART = M2
         DO 140 J=JSTART,JENDM1
C ...       Find row index equivalent to column index
            JROW   = J - M0
            JROW1  = JROW + 1
            M      = JROW
C
C ...       Find pivot element
            IEND = NROWS
            IF ( K.EQ.NBLOKS ) IEND = M0
            DO 40 I=JROW1,IEND
               IF ( ABS(A(I,J,K)).GT.ABS(A(M,J,K)) ) M=I
 40         CONTINUE
C
            JXEND  = NCOLS
            IF ( K.EQ.NBLOKS ) JXEND = NROWS
            IF ( M.EQ.JROW ) GO TO 70
C ...       Exchange rows
            DO 60 JX = J,JXEND
               T            = A(M,JX,K)
               A(M,JX,K)    = A(JROW,JX,K)
               A(JROW,JX,K) = T
 60         CONTINUE
            T           = RHS(M,K)
            RHS(M,K)    = RHS(JROW,K)
            RHS(JROW,K) = T
C
C ...       Check for singular matrix and minimum pivot
 70         T = ABS( A(JROW,J,K) )
            IF ( T.GT.PMIN ) GO TO 80
            PMIN  = T
            ISING = K
            IF ( PMIN.LT.TOL ) GO TO 990
C
C ...       The following applies to triangularizing the sub-block
C           defined by K and columns M1-M4.
 80         T    = 1.E+0/A(JROW,J,K)
            J1   = J + 1
            DO 120 I = JROW1,IEND
               AMULT = A(I,J,K)*T
               DO 100 JX = J1,JXEND
                  A(I,JX,K) = A(I,JX,K) - AMULT*A(JROW,JX,K)
 100           CONTINUE
               RHS(I,K) = RHS(I,K) - AMULT*RHS(JROW,K)
 120        CONTINUE
 140      CONTINUE
C
C ...    Test last diagonal element of sub-block in block K
 145     T  = ABS( A(JEND-M0,JEND,K) )
         IF ( T.GT.PMIN ) GO TO 150
         PMIN  = T
         ISING = K
         IF ( PMIN.LE.TOL ) GO TO 990
C
C ...    The following applies to zeroing out the sub-block defined
C        by K+1 and columns 1-M0.
 150     IF ( K.EQ.NBLOKS ) GO TO 220
         DO 200 J = M2,M4
            JROW   = J - M0
            T      = 1.E+0/A(JROW,J,K)
            JXSTRT = J - NROWS + 1
            NEND = NROWS
            IF ( K.EQ.NBLKM1 ) NEND = M0
            DO 180 I=1,NEND
               AMULT = A(I,J-NROWS,K+1)*T
               DO 160 JX = JXSTRT,NROWS
                  A(I,JX,K+1) = A(I,JX,K+1) - AMULT*A(JROW,JX+NROWS,K)
 160           CONTINUE
               RHS(I,K+1) = RHS(I,K+1) - AMULT*RHS(JROW,K)
 180        CONTINUE
 200     CONTINUE
 220  CONTINUE
C
C ... Back Substitution
C
C ... NBLOKS is a special case
      K = NBLOKS
      IF ( NROWS.EQ.2 ) GO TO 340
      NM1 = M0 - 1
      DO 320 JB = 1,NM1
         JM1      = M0 - JB
         J        = JM1 + 1
         JROW     = J + M0
         RHS(J,K) = RHS(J,K)/A(J,JROW,K)
         DO 300 I = 1,JM1
            RHS(I,K) = RHS(I,K) - A(I,JROW,K)*RHS(J,K)
 300     CONTINUE
 320  CONTINUE
 340  RHS(1,K) = RHS(1,K)/A(1,M1,K)
C
C ... Do back substitution for general K
      DO 420 KB = 1,NBLKM1
         K      = NBLKM1 - KB + 1
         ISTART = 1
         IF ( K.EQ.1 ) ISTART = M1
         DO 400 IB = ISTART,NROWS
            I    = NROWS - IB + ISTART
            SUM  = 0.E+0
            IF ( I.EQ.NROWS ) GO TO 370
            JROW = I + M1
            DO 360 J = JROW,M4
               SUM = SUM + A(I,J,K)*RHS(J-M0,K)
 360        CONTINUE
 370        DO 380 J = M5,NCOLS
               SUM = SUM + A(I,J,K)*RHS(J-M4,K+1)
 380        CONTINUE
            RHS(I,K) = (RHS(I,K) - SUM )/A(I,I+M0,K)
 400     CONTINUE
 420  CONTINUE
C ... System not singular. Indicate block with minimum pivot.
      ISING = -ISING
C
C ... Error exit
 990  RETURN
      END
