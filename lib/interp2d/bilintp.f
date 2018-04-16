C+------------------------------------------------------------------------------
C
      SUBROUTINE BILINTP (IDIM, JDIM, XY, XTARGET, YTARGET,
     >                    ITARGET, JTARGET, EPS, P, Q, STATUS)
C
C ACRONYM: BILinear INTerpolation (2-D quasi-rectangular data, Pairs form)
C          ---      ---                                        -
C
C DESCRIPTION:
C
C         BILINTP is a variant of BILINT for an (X,Y) grid stored as pairs.
C
C         Interpolation within a quasi-rectangular (sub)mesh involves:
C     (1) locating the cell (if any) containing the current target point, and
C     (2) deriving a function value from the values of the function at that
C     cell's vertices (and there may be more than one function).
C
C         This routine serves both requirements: if the cell indicated by (i,j)
C     contains the target point (x,y), an interpolated function value is given
C     by a generalization of the 4-point bilinear formula for rectangular data:
C
C         F(x,y) = (1-p)(1-q)Fi,j + p(1-q)Fi+1,j + (1-p)qFi,j+1 + pqFi+1,j+1
C
C     for some p and q in [0, 1].  If not, a direction to move in is indicated.
C
C         Two nonlinear equations in p and q are obtained by applying the above
C     to X and Y:
C
C     f1 = (1-p)(1-q)Xi,j + p(1-q)Xi+1,j + (1-p)qXi,j+1 + pqXi+1,j+1 - Xtarg = 0
C     f2 = (1-p)(1-q)Yi,j + p(1-q)Yi+1,j + (1-p)qYi,j+1 + pqYi+1,j+1 - Ytarg = 0
C
C         Solution by a Newton iteration generalizes nicely to the trilinear
C     case.  Starting guesses of 0.5 should lead to rapid convergence if the
C     cell is the correct one.  Iterates outside [-4, +4] terminate the
C     calculation immediately to save time and guard against overflow: as with
C     convergence to a solution outside [0, 1], a direction for moving i and j
C     is implied.
C
C         Singularity in the Jacobian matrix is probably an indication of
C     singularity in the mesh.  Proceeding to the next cell in the higher
C     level's search strategy is probably appropriate.  Failure with all
C     cells means either the target point is outside the mesh or the mesh
C     is corrupted.  Handling this is application-dependent.
C
C         For a comprehensive treatment of the 3D case, see the INT3D program
C     in Pieter Buning's PLOT3D Tools collection at NASA Ames.  It supports
C     PLOT3D's IBLANK feature.  This routine is an attempt to provide a more
C     general-purpose lower level module without too much loss in efficiency,
C     and with the 3D analogue in mind.
C
C         Since the functions associated with the mesh may be stored either as
C     F(i,j,1:n) PLOT3D-style or as F(1:n,i,j) and similar uncertainty applies
C     to storing the interpolated function values, the evaluations of F must
C     be left to the higher level.
C
C ARGUMENTS:
C
C     ARG    DIM   TYPE I/O/S DESCRIPTION
C     IDIM,         I     I   Dimensions of X and Y in the calling program
C     JDIM                    (all cells not necessarily active)
C     XY 2,IDIM,JDIM R    I   Coordinates of the quasi-rectangular (sub)mesh;
C                             those for the target cell are assumed meaningful
C     XTARGET,      R     I   Coordinates of the target point
C     YTARGET
C     ITARGET,      I     I   These indicate the cell to be processed on this
C     JTARGET                 call.  See STATUS.
C     EPS           R     I   Tolerance used to determine whether the computed
C                             P & Q are effectively within [0, 1], and also to
C                             help terminate the Newton iteration.  Suggestion:
C                             1.E-6 or 1.E-12 for 32- or 64-bit arithmetic.
C     P,            R     O   Interpolation coefficients as indicated above.
C     Q                       See STATUS.
C     STATUS        I     O   0:  ITARGET & JTARGET define the enclosing cell
C                                 and P & Q are inside [0-EPS, 1+EPS] and are
C                                 usable as bilinear interplation coefficients.
C                             1:  P & Q are outside [0, 1] and their values
C                                 suggest how to adjust ITARGET & JTARGET on
C                                 the next call.  E.g.:
C                                 IF (P .GT. 1.) ITARGET = MIN (ITARGET + 1, I2)
C                                 Avoiding an infinite loop may be tricky.
C                            -1:  No P & Q could be calculated because of
C                                 matrix singularity.  See discussion above.
C
C PROCEDURES:
C
C     LUSOLVE    Square system solution by LU decomposition with pivoting
C                (now in-line for efficiency)
C
C ERROR HANDLING:
C
C     None, for efficiency.  ITARGET & JTARGET are assumed to be valid.
C
C HISTORY:
C
C   07/16/93   DAS    Initial implementation, with 3D version in mind,
C                     following discussion with James Reuther and Scott
C                     Thomas, and perusal of Pieter Buning's INT3D.
C   11/24/93    "     NI, NJ arguments replaced by JDIM argument.
C                     All are actually redundant, but submesh usage is clearer.
C   11/01/96    "     In-lined LUSOLVE.
C   12/08/97    "     FNORM wasn't being tested upon hitting MXITER, yet
C                     p and q could both be in [0, 1].  Non-convex cells can
C                     still meet both convergence tests with non-unique
C                     solutions - no work-around is provided here. MXITER
C                     has been raised from 4 (really 5) to 8 in view of
C                     the missing step-halving and possibly tiny EPS.
C   01/08/98    "     FNORM was unintentionally being updated as ||correction||.
C   01/10/98    "     Application to unnormalized data requires a reference
C                     length, but we can't change the arguments at this stage,
C                     so a combination of current cell size and target (X,Y)
C                     units is used to determine a residual scale factor.
C   01/12/01    "     BILINTP adapted from BILINT for use with RIPPLE2P.
C
C AUTHOR: David Saunders, Sterling Software/NASA Ames, Moffett Field, CA.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   IDIM, JDIM, ITARGET, JTARGET, STATUS
      REAL
     >   XY (2, IDIM, JDIM), XTARGET, YTARGET, EPS, P, Q

C     Local constants:

      INTEGER
     >   MAXITER
      REAL
     >   HALF, TOOBIG, TOOLOW, ZERO
      PARAMETER
     >  (MAXITER = 8, HALF = 0.5E+0, TOOBIG = 4.E+0, TOOLOW = -3.E+0,
     >   ZERO = 0.E+0)

C     Local variables:

      INTEGER
     >   I, J, K, M
      REAL
     >   AJAC (2, 2), DP (2), FNORM, PK, QK, RTOL, T,
     >   XP, YP, XQ, YQ, XPQ, YPQ, XT, YT

C     Procedures:

C**** EXTERNAL
C****>   LUSOLVE

C     Execution:

      STATUS = 1      ! Assume most calls are searching for the right cell
      I = ITARGET
      J = JTARGET
      XT = XY (1, I, J) - XTARGET
      YT = XY (2, I, J) - YTARGET
      XQ = XY (1, I, J + 1) - XY (1, I, J)
      YQ = XY (2, I, J + 1) - XY (2, I, J)
      XP = XY (1, I + 1, J) - XY (1, I, J)
      YP = XY (2, I + 1, J) - XY (2, I, J)
      XPQ = XY (1, I + 1, J + 1) - XY (1, I, J + 1) - XP
      YPQ = XY (2, I + 1, J + 1) - XY (2, I, J + 1) - YP
      RTOL = MAX (ABS (XP), ABS (YP), ABS (XTARGET), ABS (YTARGET)) *EPS

C     For any usable solution to P & Q, the Newton iteration should not
C     need safeguarding - hence no step-halving inner iteration is included
C     to ensure that ||f|| decreases the way one would in general.

      K = 0
      PK = HALF
      QK = HALF

  300 CONTINUE

         AJAC (1, 1) = XP + QK * XPQ      ! J (i, j) = partial df(i)/dp(j)
         AJAC (2, 1) = YP + QK * YPQ
         AJAC (1, 2) = XQ + PK * XPQ
         AJAC (2, 2) = YQ + PK * YPQ

         DP (1) = PK * AJAC (1, 1) + QK * XQ + XT     ! f(1) and f (2)
         DP (2) = PK * AJAC (2, 1) + QK * YQ + YT

         FNORM = MAX (ABS (DP (1)), ABS (DP (2)))

CCC      WRITE (9, '(A, I2, A, 1P, E21.14, A, 2E9.2, A, 2E22.14, A, 2I4)')
CCC     >   ' BILINTP:', K, '  ||f||:', FNORM, '  X,YT:', XTARGET, YTARGET,
CCC     >   '  P, Q:', PK, QK, '  I,JT:', ITARGET, JTARGET

C*****   CALL LUSOLVE (2, 2, AJAC, DP, IER)       ! Solve J dp = f

C        Cray inlining doesn't seem to help, so do the 2x2 case explicitly
C        to avoid an external reference and to avoid degenerate loops:

C        Perform the LU factorization, solving L y = b as we go:

C        Determine pivot element:

         M = 1
         T = AJAC (1, 1)
         IF (ABS (AJAC (2, 1)) .GT. ABS (T)) THEN
            AJAC (1, 1) = AJAC (2, 1)
            AJAC (2, 1) = T
            M = 2
         END IF

         IF (AJAC (1, 1) .EQ. ZERO) THEN  ! Singular matrix
            STATUS = -1
            GO TO 999
         END IF

C        Store the multiplier:

         AJAC (2, 1) = -AJAC (2, 1) / AJAC (1, 1)

C        Apply the multiplier to current submatrix and RHS:

         T = AJAC (M, 2)
         AJAC (M, 2) = AJAC (1, 2)
         AJAC (1, 2) = T
         AJAC (2, 2) = AJAC (2, 2) + AJAC (2, 1) * T

         T = DP (M)
         DP (M) = DP (1)
         DP (1) = T
         DP (2) = DP (2) + AJAC (2, 1) * T

         IF (AJAC (2, 2) .EQ. ZERO) THEN
            STATUS = -1       ! Singular matrix
            GO TO 999
         END IF

C        Back substitution (solution of U x = y):

         DP (2) = DP (2) / AJAC (2, 2)
         DP (1) = (DP (1) - AJAC (1, 2) * DP (2)) / AJAC (1, 1)


         PK = PK - DP (1)     ! New solution
         QK = QK - DP (2)

C        Assume that more often than not, we're searching for the right I, J
C        rather than being in the eventual enclosing cell. Exit early if we can.

         IF (MAX (PK, QK) .GT. TOOBIG) GO TO 900
         IF (MIN (PK, QK) .LT. TOOLOW) GO TO 900

         IF (FNORM .GT. RTOL) THEN
            K = K + 1
            IF (K .LT. MAXITER) GO TO 300
         ELSE
C           Drop through - either converged for the correct cell, converged for
C           the wrong cell (a nearby one?), or unconverged.  As for a diverged
C           case, the current solution in the latter two cases should still
C           point to a better I and J for the next call.

CCC      K = K + 1
CCC      WRITE (9, '(A, I2, A, 1P, E21.14, A, 2E9.2, A, 2E22.14, A, 2I4)')
CCC     >   ' BILINTP:', K, '  ||f||:', FNORM, '  X,YT:', XTARGET, YTARGET,
CCC     >   '  P, Q:', PK, QK, '  I,JT:', ITARGET, JTARGET

            IF (MAX (ABS (PK - HALF), ABS (QK - HALF)) .LT. HALF + EPS)
     >         STATUS = 0

         END IF


  900 P = PK
      Q = QK

  999 RETURN
      END
