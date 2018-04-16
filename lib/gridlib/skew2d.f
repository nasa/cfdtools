C+------------------------------------------------------------------------------
C
      SUBROUTINE SKEW2D (IDIM, JDIM, I1, I2, J1, J2, K, X, Y, NBAD,
     >                   MAXPRINT, LUNOUT)
C
C  ONE-LINER: Checks a 2-D mesh for skew (nonconvex) cells
C
C  DESCRIPTION:
C
C        SKEW2D performs one type of grid quality check for the indicated
C     portion of a 2-D mesh.  Nonconvex cells are detected.  Other tests
C     for acceptable orthogonality are not considered.
C
C        If any bad cells are detected during a first pass, the search is
C     repeated with print-out turned on unless LUNOUT < 0.  (This approach
C     should be more efficient for the common case of good grids, at the
C     expense of duplicate code.)  No more than MAXPRINT bad cells will be
C     listed for a given call, however.
C
C        The method is that of the PLOT3D tools utility in program PROPER,
C     extended to cover both cases of nonconvexity instead of just one.
C                                                                c
C        For a typical cell as shown, the strategy is        4--------2
C     to check whether vertices 3 and 4 are on opposite   d /        / b
C     sides of diagonal e = 1-2, and similarly for the     1--------3
C     vertices 1 and 2 and diagonal 3-4.                  i,j   a  i+1,j
C
C        For a "good" cell, cross-products e x a and e x d should have
C     opposite directions.  Thus if their third components have the
C     same sign, the cell is bad.  Repeat for the other case.
C     
C        Note that, to be fail-safe, a negative-area check should be
C     performed first.  This would prevent the possibility of two
C     opposite vertices, both on their respective wrong sides of a
C     diagonal, from passing the present test.
C
C  ENVIRONMENT:  FORTRAN 77 with minor extensions
C
C  HISTORY:
C
C     21 May 1994   DAS   Initial implementation, intended to prevent
C                         subsequent radial-line redistribution from failing
C                         ungracefully.  Collapsed cells are permitted.
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER IDIM, JDIM      ! (I) Grid array dimensions in calling program
      INTEGER I1, I2,         ! (I) Indication of the subset of cells to check
     >        J1, J2
      INTEGER K               ! (I) For the case of multiple planar meshes
      REAL    X (IDIM, JDIM), ! (I) Grid coordinates
     >        Y (IDIM, JDIM)
      INTEGER NBAD            ! (O) Count of nonconvex cells found
      INTEGER MAXPRINT        ! (O) Limit on # bad cells listed if LUNOUT > 0
      INTEGER LUNOUT          ! (I) LUNOUT < 0 suppresses printout of bad cells

C     Local constants.

      REAL
     >   ZERO
      PARAMETER
     >  (ZERO = 0.0E+0)

C     Local variables.

      INTEGER
     >   I, J, NPRINT
      REAL
     >   U1, U2, V1, V2, W1, W2
      LOGICAL
     >   BAD

C     Statement functions.

      REAL
     >   A1, A2, B1, B2, CROSS3  ! 3rd component of A x B

      CROSS3 (A1, A2, B1, B2) = A1 * B2 - A2 * B1

C     Execution.

      NBAD = 0 

      DO J = J1, J2 - 1
         DO I = I1, I2 - 1

C           First, check pts. 3, 4 against diagonal 1-2 (origin at pt. 1):

            W1 = X (I, J + 1)     - X (I, J)  ! edge d = vector (W1, W2)
            V1 = X (I + 1, J)     - X (I, J)  ! edge a = vector (V1, V2)
            U1 = X (I + 1, J + 1) - X (I, J)  ! diag e = vector (U1, U2)
            W2 = Y (I, J + 1)     - Y (I, J)
            V2 = Y (I + 1, J)     - Y (I, J)
            U2 = Y (I + 1, J + 1) - Y (I, J)

C           Third components of e x a and e x d:

            IF (CROSS3 (U1, U2, V1, V2) *
     >          CROSS3 (U1, U2, W1, W2) .GT. ZERO) THEN

               NBAD = NBAD + 1

            ELSE               ! Repeat for the other diagonal

C              Check pts. 1, 2 against diagonal 3-4 (origin at pt. 3).
C              Edge a = (-V1, -V2) (negative of vector a above).

               U1 = X (I, J + 1)     - X (I + 1, J)  ! diag f = (U1, U2)
               W1 = X (I + 1, J + 1) - X (I + 1, J)  ! edge b = (W1, W2)
               U2 = Y (I, J + 1)     - Y (I + 1, J)
               W2 = Y (I + 1, J + 1) - Y (I + 1, J)

C              Third components of f x a and f x b:

               IF (CROSS3 (U1, U2, -V1, -V2) *
     >             CROSS3 (U1, U2,  W1,  W2) .GT. ZERO) NBAD = NBAD + 1

            END IF

         END DO
      END DO

C     Keeping the I/O out of the loops makes the most common NBAD = 0 case
C     more efficient.  Relocate bad cells if printout is requested.

      IF (NBAD   .EQ. 0) GO TO 999   ! Saves extra indenting.
      IF (LUNOUT .LT. 0) GO TO 999   !  " "
      NPRINT = 0

      DO J = J1, J2 - 1
         DO I = I1, I2 - 1
            W1 = X (I, J + 1)     - X (I, J)
            V1 = X (I + 1, J)     - X (I, J)
            U1 = X (I + 1, J + 1) - X (I, J)
            W2 = Y (I, J + 1)     - Y (I, J)
            V2 = Y (I + 1, J)     - Y (I, J)
            U2 = Y (I + 1, J + 1) - Y (I, J)

            BAD = CROSS3 (U1, U2, V1, V2) *
     >            CROSS3 (U1, U2, W1, W2) .GT. ZERO

            IF (BAD) THEN
               NPRINT = NPRINT + 1
            ELSE
               U1 = X (I, J + 1)     - X (I + 1, J)
               W1 = X (I + 1, J + 1) - X (I + 1, J)
               U2 = Y (I, J + 1)     - Y (I + 1, J)
               W2 = Y (I + 1, J + 1) - Y (I + 1, J)

               BAD = CROSS3 (U1, U2, -V1, -V2) *
     >               CROSS3 (U1, U2,  W1,  W2) .GT. ZERO

               IF (BAD) NPRINT = NPRINT + 1
            END IF

            IF (BAD) THEN
               WRITE (LUNOUT,'(A,3I5)') ' Nonconvex cell (I,J,K):',I,J,K
               IF (NPRINT .GE. MAXPRINT) GO TO 999
            END IF

         END DO
      END DO

  999 RETURN
      END
