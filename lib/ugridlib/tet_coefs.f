C+------------------------------------------------------------------------------
C
      SUBROUTINE TET_COEFS (X1, X2, X3, X4, XT, PQRS)
C
C     One-liner:  Interpolation coefficients for tetrahedral grid cells
C
C     Description:
C
C        TET_COEFS calculates coefficients p, q, r, s such that, for the given
C     vertices X1, X2, X3, X4 of a tetrahedron and a target point XT,
C
C        p X1  +  q X2  +  r X3  +  s X4  =  XT
C          p   +    q   +    r   +    s   =  1
C
C     Eliminating s leads to a 3 x 3 linear system solved by LU decomposition.
C     A degenerate cell (two or more common points) would lead to singularity.
C     This should be guarded against at a higher level.
C
C     The coefficients may be used for linear interpolation of other quanitities
C     within the tetrahedron.  If XT is outside the cell, some coefficients may
C     be negative.  This could be employed as a means of searching for the cell
C     containing a target point.  However, for large meshes, nearest-cell
C     searches based on Alternating Digital Tree techniques are recommended.
C
C     08/11/04  DAS  Adaptation of PROJECT3 (triangle utility).
C
C     Author:   David Saunders, ELORET/NASA Ames Research Center, CA
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      REAL, INTENT (IN), DIMENSION (3) ::
     >   X1, X2, X3, X4,     ! (x,y,z) coordinates of the tetrahedron vertices
     >   XT                  ! Target (x,y,z) point

      REAL, INTENT (OUT), DIMENSION (4) ::
     >   PQRS                ! Interpolation coefficients p, q, r, s; sum = 1

C     Local variables:

      INTEGER
     >   I, IER

      REAL
     >   A(3,3)

C     Execution:

      DO I = 1, 3
         A(I,1)  = X1(I) - X4(I)
         A(I,2)  = X2(I) - X4(I)
         A(I,3)  = X3(I) - X4(I)
         PQRS(I) = XT(I) - X4(I)  ! RHS b
      END DO

C     Solution of Ax = b by LU decomposition.  RHS is overwritten with solution.

      CALL LUSOLVE (3, 3, A, PQRS, IER)

      IF (IER /= 0) THEN
         WRITE (*, '(/, A)')
     >      ' TET_COEFS:  Singular matrix.  Collapsed cell?  Aborting.'
         WRITE (*, '(A, 1P, 3E19.11)')
     >      ' X1:', X1, ' X2:', X2, ' X3:', X3, ' X4:', X4, 'XT:', XT
         STOP
      END IF

      PQRS(4) = 1.E+0 - (PQRS(1) + PQRS(2) + PQRS(3))

      END SUBROUTINE TET_COEFS
