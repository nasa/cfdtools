C+------------------------------------------------------------------------------
C
      SUBROUTINE TRILINT (IDIM, JDIM, KDIM, X, Y, Z, XTARGET, YTARGET,
     >                    ZTARGET, ITARGET, JTARGET, KTARGET, EPS,
     >                    P, Q, R, STATUS)
C
C ACRONYM: TRILinear INTerpolation (quasi-rectangular data in 3 dimensions)
C          ----      ---
C
C DESCRIPTION:
C
C         Interpolation within a quasi-rectangular (sub)mesh involves:
C     (1) locating the cell (if any) containing the current target point, and
C     (2) deriving a function value from the values of the function at that
C     cell's vertices (and there may be more than one function).
C
C         This routine serves both purposes: if the cell indicated by (i,j,k)
C     contains the target point (x,y,z), an interpolated function value is given
C     by a generalization of the 8-point linear formula for rectangular data:
C
C         F(x,y,z) = (1-p)(1-q)(1-r)Fi,j,k     +  p(1-q)(1-r)Fi+1,j,k     +
C                    (1-p)  q  (1-r)Fi,j+1,k   +  p  q  (1-r)Fi+1,j+1,k   +
C                    (1-p)(1-q)  r  Fi,j,k+1   +  p(1-q)  r  Fi+1,j,k+1   +
C                    (1-p)  q    r  Fi,j+1,k+1 +  p  q    r  Fi+1,j+1,k+1
C
C     for some p, q & r in [0, 1].  If not, a direction to move in is indicated.
C
C         Three nonlinear equations in p, q & r are obtained by applying the
C     above to X, Y & Z.  A Newton iteration is initiated with starting guesses
C     of 0.5, which should lead to rapid convergence if the cell is the correct
C     one.  Iterates outside [-4, +4] terminate the calculation immediately to
C     save time and guard against overflow: as with convergence to a solution
C     outside [0, 1], a direction for moving i, j & k is implied.
C
C         Singularity in the Jacobian matrix is probably an indication of
C     singularity in the mesh.  Proceeding to the next cell in the higher
C     level's search strategy is probably appropriate.  Failure with all
C     cells means either the target point is outside the mesh or the mesh
C     is corrupted.  Handling this is application-dependent.
C
C         For a comprehensive treatment of this 3D case, see the INT3D program
C     in Pieter Buning's PLOT3D Tools collection at NASA Ames.  It supports
C     PLOT3D's IBLANK feature.  This routine is an attempt to provide a more
C     general-purpose lower level module without too much loss in efficiency.
C
C         Since the functions associated with the mesh may be stored either as
C     F(i,j,k,1:n) PLOT3D-style or as F(1:n,i,j,k) and similar uncertainty
C     applies to storing the interpolated function values, the evaluations of
C     F must be left to the higher level.
C
C ARGUMENTS:
C
C     ARG    DIM   TYPE I/O/S DESCRIPTION
C     IDIM,         I     I   Dimensions of X, Y & Z (*, *) in the calling
C     JDIM,                   program (all cells not necessarily active)
C     KDIM
C     X, IDIM,JDIM,KDIM R I   Coordinates of the quasi-rectangular mesh;
C     Y,                      those for the specified target cell are assumed
C     Z                       to be meaningful
C     XTARGET,      R     I   Coordinates of the target point
C     YTARGET,
C     ZTARGET
C     ITARGET,      I     I   These indicate the cell to be processed on this
C     JTARGET,                call.  See STATUS.
C     KTARGET
C     EPS           R     I   Tolerance used to determine whether the computed
C                             P, Q & R are effectively within [0, 1], and also
C                             to help end the Newton iteration.  Suggestion:
C                             1.E-6 or 1.D-12 for 32- or 64-bit arithmetic.
C     P,            R     O   Interpolation coefficients as indicated above.
C     Q,                      See STATUS.
C     R
C     STATUS        I     O   0:  I/J/KTARGET define the enclosing cell and
C                                 P, Q & R are inside [0-EPS, 1+EPS] and are
C                                 usable as trilinear interplation coefficients.
C                             1:  P, Q & R are outside [0, 1] and their values
C                                 suggest how to adjust I/J/KTARGET on the next
C                                 call.  E.g.:
C                                 IF (P .GT. 1.) ITARGET = MIN (ITARGET + 1, I2)
C                                 Avoiding an infinite loop may be tricky.
C                            -1:  No P, Q, & R could be calculated because of
C                                 matrix singularity.  See discussion above.
C
C PROCEDURES:
C
C     LUSOLVE    Square system solution by LU decomposition with pivoting
C
C ERROR HANDLING:
C
C     None, for efficiency.  I/J/KTARGET are assumed to be valid.
C
C ENVIRONMENT:   FORTRAN 77 with minor extensions
C
C HISTORY:
C
C   07/28/93   DAS    Initial implementation adapted from 2D version (BILINT)
C                     following discussion with James Reuther and Scott
C                     Thomas, and perusal of Pieter Buning's INT3D.
C   11/24/93    "     Replaced NI, NJ, NK with KDIM.  (All are actually
C                     redundant, but usage with submeshes is clearer.)
C   03/25/04    "     Boundary layer grid searches suggested bumping up
C                     TOOBIG from 4. to 10000., and MAXITER from 4 to 8.
C
C AUTHOR: David Saunders, Sterling Software/NASA Ames, Moffett Field, CA.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   IDIM, JDIM, KDIM, ITARGET, JTARGET, KTARGET, STATUS
      REAL
     >   X (IDIM, JDIM, KDIM), Y (IDIM, JDIM, KDIM),
     >   Z (IDIM, JDIM, KDIM), XTARGET, YTARGET, ZTARGET, EPS, P, Q, R

C     Local constants:

      INTEGER
     >   MAXITER
      REAL
     >   HALF, TOOBIG
      PARAMETER
     >  (MAXITER = 8, HALF = 0.5E+0, TOOBIG = 1.E+4)

C     Local variables:

      INTEGER
     >   I, J, K, L, IER
      REAL
     >   AJAC (3, 3), DP (3), PL, QL, RL, XP, YP, ZP, XQ, YQ, ZQ,
     >   XR, YR, ZR, XPQ, YPQ, ZPQ, XQR, YQR, ZQR, XRP, YRP, ZRP,
     >   XPQR, YPQR, ZPQR, XT, YT, ZT

C     Procedures:

      EXTERNAL
     >   LUSOLVE

C     Execution:

      STATUS = 1      ! Assume most calls are searching for the right cell

      I = ITARGET
      J = JTARGET
      K = KTARGET

      XT = X (I, J, K) - XTARGET
      YT = Y (I, J, K) - YTARGET
      ZT = Z (I, J, K) - ZTARGET

      XP = X (I + 1, J, K) - X (I, J, K)
      YP = Y (I + 1, J, K) - Y (I, J, K)
      ZP = Z (I + 1, J, K) - Z (I, J, K)

      XQ = X (I, J + 1, K) - X (I, J, K)
      YQ = Y (I, J + 1, K) - Y (I, J, K)
      ZQ = Z (I, J + 1, K) - Z (I, J, K)

      XR = X (I, J, K + 1) - X (I, J, K)
      YR = Y (I, J, K + 1) - Y (I, J, K)
      ZR = Z (I, J, K + 1) - Z (I, J, K)

      XPQ = X (I + 1, J + 1, K) - X (I, J + 1, K) - XP
      YPQ = Y (I + 1, J + 1, K) - Y (I, J + 1, K) - YP
      ZPQ = Z (I + 1, J + 1, K) - Z (I, J + 1, K) - ZP

      XRP = X (I + 1, J, K + 1) - X (I, J, K + 1) - XP
      YRP = Y (I + 1, J, K + 1) - Y (I, J, K + 1) - YP
      ZRP = Z (I + 1, J, K + 1) - Z (I, J, K + 1) - ZP

      XQR = X (I, J + 1, K + 1) - X (I, J, K + 1) - XQ
      YQR = Y (I, J + 1, K + 1) - Y (I, J, K + 1) - YQ
      ZQR = Z (I, J + 1, K + 1) - Z (I, J, K + 1) - ZQ

      XPQR = X (I + 1, J + 1, K + 1) - X (I, J + 1, K + 1) -
     >       X (I + 1, J,     K + 1) + X (I, J, K + 1)     - XPQ
      YPQR = Y (I + 1, J + 1, K + 1) - Y (I, J + 1, K + 1) -
     >       Y (I + 1, J,     K + 1) + Y (I, J, K + 1)     - YPQ
      ZPQR = Z (I + 1, J + 1, K + 1) - Z (I, J + 1, K + 1) -
     >       Z (I + 1, J,     K + 1) + Z (I, J, K + 1)     - ZPQ

C     For any usable solution to P, Q & R, the Newton iteration should not
C     need safeguarding - hence no step-halving inner iteration is included
C     to ensure that ||f|| decreases the way one would in general.

      L = 0
      PL = HALF
      QL = HALF
      RL = HALF

  300 CONTINUE         ! Jacobian (i, j) = partial df(i)/dp(j)

         AJAC (1, 1) = QL * (RL * XPQR + XPQ) + RL * XRP + XP
         AJAC (2, 1) = QL * (RL * YPQR + YPQ) + RL * YRP + YP
         AJAC (3, 1) = QL * (RL * ZPQR + ZPQ) + RL * ZRP + ZP

         AJAC (1, 2) = RL * (PL * XPQR + XQR) + PL * XPQ + XQ
         AJAC (2, 2) = RL * (PL * YPQR + YQR) + PL * YPQ + YQ
         AJAC (3, 2) = RL * (PL * ZPQR + ZQR) + PL * ZPQ + ZQ

         AJAC (1, 3) = PL * (QL * XPQR + XRP) + QL * XQR + XR
         AJAC (2, 3) = PL * (QL * YPQR + YRP) + QL * YQR + YR
         AJAC (3, 3) = PL * (QL * ZPQR + ZRP) + QL * ZQR + ZR

C        RHS elements are f (1), f (2), f (3):

         DP (1) = PL * AJAC (1, 1) + QL * (RL * XQR + XQ) + RL * XR + XT
         DP (2) = PL * AJAC (2, 1) + QL * (RL * YQR + YQ) + RL * YR + YT
         DP (3) = PL * AJAC (3, 1) + QL * (RL * ZQR + ZQ) + RL * ZR + ZT

         CALL LUSOLVE (3, 3, AJAC, DP, IER)       ! Solve J dp = f

         IF (IER .NE. 0) THEN
            STATUS = -1       ! Singular matrix
            GO TO 999
         END IF

         PL = PL - DP (1)     ! New solution
         QL = QL - DP (2)
         RL = RL - DP (3)

C        Assume that more often than not, we're searching for the right I, J, K
C        rather than being in the eventual enclosing cell. Exit early if we can.

         IF (MAX (ABS (PL), ABS (QL), ABS (RL)) .GT. TOOBIG) GO TO 900

         IF (MAX (ABS (DP (1)), ABS (DP (2)), ABS (DP (3))).GT.EPS) THEN
            L = L + 1
            IF (L .LE. MAXITER) GO TO 300
         END IF

C     Drop through - either converged for the correct cell, converged for
C     the wrong cell (a nearby one?), or unconverged.  As with the diverged
C     case, the current solution in the latter two cases should still point
C     to a better I, J & K for the next call.

      IF (MAX (ABS (PL - HALF), ABS (QL - HALF), ABS (RL - HALF)) .LT.
     >    HALF + EPS) STATUS = 0

  900 P = PL
      Q = QL
      R = RL

  999 RETURN
      END
