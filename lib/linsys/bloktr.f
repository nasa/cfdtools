C+----------------------------------------------------------------------
C
      SUBROUTINE BLOKTR (NDIM, NSIZE, NBLOKS, L,D,U, B, A, R, IP, IER)
C
C  PURPOSE:
C     BLOKTR solves a block-tridiagonal system of equations (one right
C     hand side, with all blocks the same size in any one system).
C     It overwrites the right hand side vector with the solution.
C     Upper diagonal blocks are also overwritten.
C
C  ARGUMENTS:
C  VARIABLE  TYPE  I/O/S  DIMENSIONS  DESCRIPTION
C
C    NDIM     INT    I       -        Declared row dimension of all
C                                     multi-dimensional arrays being
C                                     passed.  (For 3-D arrays, the
C                                     first two dimensions must be
C                                     declared the same in the calling
C                                     program.) NDIM >= NSIZE. This
C                                     permits use of the same arrays for
C                                     systems with differing block-size.
C    NSIZE    INT    I       -        Size of each block in the system
C                                     being solved.
C    NBLOKS   INT    I       -        Number of blocks along the main
C                                     diagonal.
C    L        REAL   I   NDIM,NDIM,   Input with the lower diagonal
C                        NBLOKS       blocks.  Not overwritten.
C                                     ( L(*,*,1) unused here. )
C    D        REAL   I   NDIM,NDIM,   Input with the blocks
C                        NBLOKS       along the main diagonal.
C                                     Not overwritten.
C    U        REAL  I/S  NDIM,NDIM,   Input with the upper diagonal
C                        NBLOKS       blocks; overwritten.
C                                     ( U(*,*,NBLOKS) unused here. )
C    B        REAL  I/O  NDIM,NBLOKS  Input with the right hand side
C                                     vector.  Output with the solution.
C    A        REAL   S   NDIM,NDIM    Work space for dense system
C                                     solver  DECOMP/SOLVE.
C    R        REAL   S   NDIM,NDIM+1  Workspace for  SOLVE.
C    IP       INT    S   NDIM         Pivot vector for dense systems.
C    IER      INT    O       -        0 on return if solution OK;
C                                     N on return if (modified) diagonal
C                                     diagonal block N is singular.
C  EXTERNAL REFERENCES:
C     DECOMP     Factorizes dense systems (LINSYS library, not FORSYTHE)
C     SOLVE      Completes solution of dense systems
C
C  ENVIRONMENT:  VAX/VMS, FORTRAN 77.
C
C  NOTES:
C     (1)  This implementation avoids explicit inverses.
C
C     (2)  For multiple right-hand-sides,  with a single left-hand-
C          side, use  DECBT and SOLBT  instead. These routines also
C          permit a third block top and bottom,  not allowed for in
C          this version.
C
C     (3)  Note that the third dimension of  L  and  U  is the same
C          as for  D, even though blocks  L(1)  and  U(NBLOKS)  are
C          unused here.  This is intended for convenience when  the
C          block-tridiagonal system is being set up by the user....
C
C     (4)  THe calling program should always test  error  flag  IER
C          upon return from BLOKTR. Any non-zero value indicates an
C          early termination of the solution due to apparent matrix
C          singularity. (This usually indicates bad input.)
C
C  REFERENCE:  Isaacson and Keller, 1966    Pages 58-60
C
C  HISTORY:  Sept. '78  DS/MZ  Original implementation.
C            Sept. '91  DAS    FORTRAN 77-ized; avoid some tests in
C                              inner loops by treating first block
C                              separately.
C
C  AUTHORS:  David Saunders, Mark Ziring, NASA Ames/Informatics Inc.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   NDIM, NSIZE, NBLOKS, IP(NSIZE), IER
      REAL
     >   L(NDIM,NDIM,NBLOKS), D(NDIM,NDIM,NBLOKS), U(NDIM,NDIM,NBLOKS),
     >   B(NDIM,NBLOKS), A(NDIM,NSIZE), R(NDIM,NSIZE+1)

C     Local constants:

      REAL
     >   ZERO
      PARAMETER
     >  (ZERO = 0.E+0)

C     Local variables:

      INTEGER
     >   I, J, K, M, NRHS
      REAL
     >   SUM

C     Execution:

      IER  = 0
      NRHS = NSIZE + 1

C  *  Forward pass:

      DO 500, I = 1, NBLOKS

C  *     Modify Ith diagonal block:

         IF (I .EQ. 1) THEN

            DO 130, K = 1, NSIZE
               DO 120, J = 1, NSIZE
                  A(J,K) = D(J,K,1)
  120          CONTINUE
  130       CONTINUE

         ELSE

            DO 160, K = 1, NSIZE
               DO 150, J = 1, NSIZE
                  SUM = ZERO
                  DO 140, M = 1, NSIZE
                     SUM = L(J,M,I)*U(M,K,I-1) + SUM
  140             CONTINUE
                  A(J,K) = D(J,K,I) - SUM
  150          CONTINUE
  160       CONTINUE

         END IF

C  *     Solve dense system for intermediate results.
C        First, LU decomposition by Gaussian elimination:

         CALL DECOMP (NSIZE, NDIM, A, IP)

C  *     Fail and return if  DECOMP  detects singularity in
C        modified diagonal block. (Could mean bad input.)

         IF (IP(NSIZE) .EQ. 0) THEN
            IER = I
            GO TO 999
         END IF

C  *     Set up composite right-hand-side:

         IF (I .EQ. 1) THEN

            DO 320, K = 1, NSIZE
               DO 310, J = 1, NSIZE
                  R(J,K+1) = U(J,K,I)
  310          CONTINUE
               R(K,1) = B(K,I)
  320       CONTINUE

         ELSE

            DO 340, K = 1, NSIZE
               SUM = ZERO
               DO 330, J = 1, NSIZE
                  R(J,K+1) = U(J,K,I)     ! Redundant if I=NBLOKS)
                  SUM = L(K,J,I)*B(J,I-1) + SUM
  330          CONTINUE
               R(K,1) = B(K,I) - SUM
  340       CONTINUE

         END IF

         IF (I .EQ. NBLOKS) NRHS = 1

         DO 400, J = 1, NRHS
            CALL SOLVE (NSIZE, NDIM, A, R(1,J), IP)
  400    CONTINUE

C  *     Transfer intermediate results, overwriting original
C        upper block and right-hand-side vector:

         DO 430, K = 1, NSIZE
            B(K,I) = R(K,1)
            IF (I .LT. NBLOKS ) THEN
               DO 420, J = 1, NSIZE
                  U(J,K,I) = R(J,K+1)
  420          CONTINUE
            END IF
  430    CONTINUE

  500 CONTINUE

C  *  End of forward pass.
C     Now perform back substitution:

      DO 600, I = 2, NBLOKS
         J = NBLOKS + 1 - I
         DO 530, K = 1, NSIZE
            SUM = ZERO
            DO 520, M = 1, NSIZE
               SUM = SUM + U(K,M,J)*B(M,J+1)
  520       CONTINUE
            B(K,J) = B(K,J) - SUM
  530    CONTINUE
  600 CONTINUE

  999 RETURN
      END
