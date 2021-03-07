C+------------------------------------------------------------------------------
C
      SUBROUTINE BILINT (IDIM, JDIM, X, Y, XTARGET, YTARGET,
     >                   ITARGET, JTARGET, EPS, P, Q, STATUS)
C
C ACRONYM: BILinear INTerpolation (quasi-rectangular data in 2 dimensions)
C          ---      ---
C
C DESCRIPTION:
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
C     X,  IDIM,JDIM R     I   Coordinates of the quasi-rectangular (sub)mesh;
C     Y                       those for the target cell are assumed meaningful
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
C
C AUTHOR: David Saunders, Sterling Software/NASA Ames, Moffett Field, CA.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   IDIM, JDIM, ITARGET, JTARGET, STATUS
      REAL
     >   X (IDIM, JDIM), Y (IDIM, JDIM), XTARGET, YTARGET, EPS, P, Q

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
      XT = X (I, J) - XTARGET
      YT = Y (I, J) - YTARGET
      XQ = X (I, J + 1) - X (I, J)
      YQ = Y (I, J + 1) - Y (I, J)
      XP = X (I + 1, J) - X (I, J)
      YP = Y (I + 1, J) - Y (I, J)
      XPQ = X (I + 1, J + 1) - X (I, J + 1) - XP
      YPQ = Y (I + 1, J + 1) - Y (I, J + 1) - YP
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
CCC     >   ' BILINT:', K, '  ||f||:', FNORM, '  X,YT:', XTARGET, YTARGET,
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
CCC     >   ' BILINT:', K, '  ||f||:', FNORM, '  X,YT:', XTARGET, YTARGET,
CCC     >   '  P, Q:', PK, QK, '  I,JT:', ITARGET, JTARGET

            IF (MAX (ABS (PK - HALF), ABS (QK - HALF)) .LT. HALF + EPS)
     >         STATUS = 0

         END IF


  900 P = PK
      Q = QK

  999 RETURN
      END
C+------------------------------------------------------------------------------
C
      SUBROUTINE CHORDS3D (N, X, Y, Z, NORMALIZ, TOTAL, CHORD)
C
C  ONE-LINER: Cumulative [relative] chord-lengths for 3-space geometric curve
C
C  DESCRIPTION:
C
C        CHORDS3D computes the cumulative Euclidean distances for two or
C     more points on a 3-space curve represented by three arrays, with an
C     option to normalize all values by the TOTAL distance.  Thus CHORD(1)
C     is returned as 0. and CHORD(N) may be either 1. or TOTAL depending
C     on whether NORMALIZ is true or false.
C
C        CHORDS3D was introduced in spite of the existing CHORD3D function
C     because use of CHORD3D with the double precision DTNURBS library in a
C     single precision application such as SMOOTH would clash with prior
C     use of the single precision version of CHORD3D.  The functionality
C     is a little different (multiple values per call, and an option to
C     normalize), and since NURBS are intended for geometric data (i.e.,
C     X and Y expected to have similar units), the careful safeguarding
C     of CHORD3D is eschewed.
C
C  ARGUMENTS:
C
C     Name    Type/Dimension  I/O/S  Description
C
C     N         I             I      Number of points on the curve. N >= 2.
C
C     X,        R (N)         I      Coordinates of data points.
C     Y,
C     Z
C
C     NORMALIZ  L             I      .TRUE. means normalize results
C                                    to the interval [0, 1].
C
C     TOTAL     R               O    Total chord length, returned
C                                    because it is otherwise lost
C                                    if results are normalized.
C
C     CHORD     R (N)           O    Cumulative chord lengths as
C                                    described above.
C
C  ENVIRONMENT:  FORTRAN 77 with minor extensions
C
C  HISTORY:
C
C     13 Mar. 1992   DAS   Analog of CHORDS2D for use with DTNURBS.
C     08 Dec. 1993    "    Added TOTAL argument.
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     &   N
      REAL
     &   CHORD (N), TOTAL, X (N), Y (N), Z (N)
      LOGICAL
     &   NORMALIZ

C     Local constants.

      REAL
     &   ZERO, ONE
      PARAMETER
     &  (ZERO = 0.0E+0,
     &   ONE  = 1.0E+0)

C     Local variables.

      INTEGER
     &   I
      REAL
     &   DINV

C     Execution.

      CHORD (1) = ZERO

      DO 10, I = 2, N
         CHORD (I) = CHORD (I - 1) +
     &      SQRT ((X (I) - X (I - 1)) ** 2 + (Y (I) - Y (I - 1)) ** 2 +
     &            (Z (I) - Z (I - 1)) ** 2)
   10 CONTINUE

      TOTAL = CHORD (N)

      IF (NORMALIZ) THEN
         DINV = ONE / TOTAL
         DO 20, I = 2, N
            CHORD (I) = CHORD (I) * DINV
   20    CONTINUE
         CHORD (N) = ONE
      END IF

      RETURN
      END
C+------------------------------------------------------------------------------
C
      SUBROUTINE DELQ3D (IMIN, IMAX, JMIN, JMAX, KMIN, KMAX,
     >                   I1, I2, J1, J2, K1, K2, X0, Y0, Z0, S0,
     >                   DFACEI, DFACEJ, DFACEK, X, Y, Z)
C
C  ONE-LINER: 3-space surface perturbations (not XYZs) given new corners
C
C  DESCRIPTION:
C
C        DELQ3D performs stage 1 of the WARPQ3D 3-space surface grid
C     perturbation in a form which is reusable by WARP3D.  It returns
C     face perturbations rather than perturbed face coordinates, to
C     avoid work-space problems in WARP3D's first stage.  See WARP2D
C     and WARPQ3D for further details of these two-stage algorithms.
C
C        The three cases of a block face are handled here by three similar
C     code sections.  Eliminating WARP2D's special handling of the case of
C     fixed corners helps keep the bulk down.
C
C  HISTORY:
C
C     12/20/94  DAS/JJR  Adaptation of WARPQ3D for WARP3D.
C     12/23/94    DAS    DELQ3D is now used by WARPQ3D.
C     02/08/96     "     DELQ3D does only stage 1 now. Use on a subgrid with
C                        full-grid S0 requires transforming S0 here.
C
C  AUTHOR:  David Saunders/James Reuther, NASA Ames, Mt. View, CA.
C
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER IMIN, IMAX, JMIN, JMAX, KMIN, KMAX      ! I Grid array dimensions.

      INTEGER I1, I2, J1, J2, K1, K2                  ! I Define active face,
                                                      !   one pair being equal.

      REAL    X0 (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX),   ! I Original face coords.
     >        Y0 (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX),   !   in appropriate places.
     >        Z0 (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX)

      REAL    S0 (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX, 3) ! I Relative arc-lengths
                                                      !   for the I,J,K lines -
                                                      !   see PARAMXYZ.  If from
                                                      !   a grid larger than the
                                                      !   active subgrid, S0 is
                                                      !   transformed here.

      REAL    DFACEI (3, JMIN:JMAX, KMIN:KMAX),    ! O Reqd. face perturbations:
     >        DFACEJ (3, IMIN:IMAX, KMIN:KMAX),    !   DFACEI (1:3,J1:J2,K1:K2)
     >        DFACEK (3, IMIN:IMAX, JMIN:JMAX)     !   = dX,dY,dZ on an I face
                                                   !   of a 3-space grid block.

      REAL    X (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX), ! I Grid coordinates:
     >        Y (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX), !   new edges of a face in;
     >        Z (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX)  !   unchanged on output.

C-------------------------------------------------------------------------------

C     Local constants.

      REAL    EPS, ONE

      PARAMETER (EPS = 1.E-8, ONE = 1.E+0) ! EPS safeguards a divide by zero -
                                           ! presumably only if result is zero.
C     Local variables.

      INTEGER I, J, K
      REAL    DELI, DELJ, DELK, S0I1, S0J1, S0K1, SRANGEI, SRANGEJ,
     >        SRANGEK, WTI1, WTI2, WTJ1, WTJ2, WTK1, WTK2

C     Execution.
C     ----------


      IF (I1 .EQ. I2) THEN

C        I plane case:
C        -------------

         I = I1

C        Set up the corner perturbations:

         DO K = K1, K2, K2 - K1
            DO J = J1, J2, J2 - J1
               DFACEI (1, J, K) = X (I, J, K) - X0 (I, J, K)
               DFACEI (2, J, K) = Y (I, J, K) - Y0 (I, J, K)
               DFACEI (3, J, K) = Z (I, J, K) - Z0 (I, J, K)
            END DO
         END DO

C        Set up intermediate edge perturbations corresponding to the
C        final corners but otherwise derived from the original edges.

         DO J = J1, J2, J2 - J1
            S0K1 = S0 (I, J, K1, 3)
            SRANGEK = ONE / (S0 (I, J, K2, 3) - S0K1)
            DO K = K1 + 1, K2 - 1
               WTK2 = (S0 (I, J, K, 3) - S0K1) * SRANGEK
               WTK1 = ONE - WTK2
               DFACEI (1, J, K) = WTK1 * DFACEI (1, J, K1) +
     >                            WTK2 * DFACEI (1, J, K2)
               DFACEI (2, J, K) = WTK1 * DFACEI (2, J, K1) +
     >                            WTK2 * DFACEI (2, J, K2)
               DFACEI (3, J, K) = WTK1 * DFACEI (3, J, K1) +
     >                            WTK2 * DFACEI (3, J, K2)
            END DO
         END DO

         DO K = K1, K2, K2 - K1
            S0J1 = S0 (I, J1, K, 2)
            SRANGEJ = ONE / (S0 (I, J2, K, 2) - S0J1)
            DO J = J1 + 1, J2 - 1
               WTJ2 = (S0 (I, J, K, 2) - S0J1) * SRANGEJ
               WTJ1 = ONE - WTJ2
               DFACEI (1, J, K) = WTJ1 * DFACEI (1, J1, K) +
     >                            WTJ2 * DFACEI (1, J2, K)
               DFACEI (2, J, K) = WTJ1 * DFACEI (2, J1, K) +
     >                            WTJ2 * DFACEI (2, J2, K)
               DFACEI (3, J, K) = WTJ1 * DFACEI (3, J1, K) +
     >                            WTJ2 * DFACEI (3, J2, K)
            END DO
         END DO

C        Interpolate the intermediate perturbations of interior points.
C        The contributions from each pair of edges are not independent.

         DO K = K1 + 1, K2 - 1
            DO J = J1 + 1, J2 - 1
               WTJ2 = (S0 (I, J,  K, 2) - S0 (I, J1, K, 2)) /
     >                (S0 (I, J2, K, 2) - S0 (I, J1, K, 2))
               WTJ1 = ONE - WTJ2
               WTK2 = (S0 (I, J, K,  3) - S0 (I, J, K1, 3)) /
     >                (S0 (I, J, K2, 3) - S0 (I, J, K1, 3))
               WTK1 = ONE - WTK2

               DELJ = WTJ1 * DFACEI (1, J1, K) + WTJ2 * DFACEI (1, J2,K)
               DELK = WTK1 * DFACEI (1, J, K1) + WTK2 * DFACEI (1, J,K2)

               DFACEI (1, J, K) = (ABS (DELJ) *DELJ + ABS (DELK) *DELK)/
     >                             MAX (ABS (DELJ) + ABS (DELK), EPS)

               DELJ = WTJ1 * DFACEI (2, J1, K) + WTJ2 * DFACEI (2, J2,K)
               DELK = WTK1 * DFACEI (2, J, K1) + WTK2 * DFACEI (2, J,K2)

               DFACEI (2, J, K) = (ABS (DELJ) *DELJ + ABS (DELK) *DELK)/
     >                             MAX (ABS (DELJ) + ABS (DELK), EPS)

               DELJ = WTJ1 * DFACEI (3, J1, K) + WTJ2 * DFACEI (3, J2,K)
               DELK = WTK1 * DFACEI (3, J, K1) + WTK2 * DFACEI (3, J,K2)

               DFACEI (3, J, K) = (ABS (DELJ) *DELJ + ABS (DELK) *DELK)/
     >                             MAX (ABS (DELJ) + ABS (DELK), EPS)
            END DO
         END DO

      ELSE IF (J1 .EQ. J2) THEN

C        J plane case:
C        -------------

         J = J1

C        Corner perturbations:

         DO K = K1, K2, K2 - K1
            DO I = I1, I2, I2 - I1
               DFACEJ (1, I, K) = X (I, J, K) - X0 (I, J, K)
               DFACEJ (2, I, K) = Y (I, J, K) - Y0 (I, J, K)
               DFACEJ (3, I, K) = Z (I, J, K) - Z0 (I, J, K)
            END DO
         END DO

C        Intermediate edge perturbations:

         DO I = I1, I2, I2 - I1
            S0K1 = S0 (I, J, K1, 3)
            SRANGEK = ONE / (S0 (I, J, K2, 3) - S0K1)
            DO K = K1 + 1, K2 - 1
               WTK2 = (S0 (I, J, K, 3) - S0K1) * SRANGEK
               WTK1 = ONE - WTK2
               DFACEJ (1, I, K) = WTK1 * DFACEJ (1, I, K1) +
     >                            WTK2 * DFACEJ (1, I, K2)
               DFACEJ (2, I, K) = WTK1 * DFACEJ (2, I, K1) +
     >                            WTK2 * DFACEJ (2, I, K2)
               DFACEJ (3, I, K) = WTK1 * DFACEJ (3, I, K1) +
     >                            WTK2 * DFACEJ (3, I, K2)
            END DO
         END DO

         DO K = K1, K2, K2 - K1
            S0I1 = S0 (I1, J, K, 1)
            SRANGEI = ONE / (S0 (I2, J, K, 1) - S0I1)
            DO I = I1 + 1, I2 - 1
               WTI2 = (S0 (I, J, K, 1) - S0I1) * SRANGEI
               WTI1 = ONE - WTI2
               DFACEJ (1, I, K) = WTI1 * DFACEJ (1, I1, K) +
     >                            WTI2 * DFACEJ (1, I2, K)
               DFACEJ (2, I, K) = WTI1 * DFACEJ (2, I1, K) +
     >                            WTI2 * DFACEJ (2, I2, K)
               DFACEJ (3, I, K) = WTI1 * DFACEJ (3, I1, K) +
     >                            WTI2 * DFACEJ (3, I2, K)
            END DO
         END DO

C        Intermediate perturbations of interior points:

         DO K = K1 + 1, K2 - 1
            DO I = I1 + 1, I2 - 1
               WTI2 = (S0 (I,  J, K, 1) - S0 (I1, J, K, 1)) /
     >                (S0 (I2, J, K, 1) - S0 (I1, J, K, 1))
               WTI1 = ONE - WTI2
               WTK2 = (S0 (I, J, K,  3) - S0 (I, J, K1, 3)) /
     >                (S0 (I, J, K2, 3) - S0 (I, J, K1, 3))
               WTK1 = ONE - WTK2

               DELI = WTI1 * DFACEJ (1, I1, K) + WTI2 * DFACEJ (1, I2,K)
               DELK = WTK1 * DFACEJ (1, I, K1) + WTK2 * DFACEJ (1, I,K2)

               DFACEJ (1, I, K) = (ABS (DELI) *DELI + ABS (DELK) *DELK)/
     >                             MAX (ABS (DELI) + ABS (DELK), EPS)

               DELI = WTI1 * DFACEJ (2, I1, K) + WTI2 * DFACEJ (2, I2,K)
               DELK = WTK1 * DFACEJ (2, I, K1) + WTK2 * DFACEJ (2, I,K2)

               DFACEJ (2, I, K) = (ABS (DELI) *DELI + ABS (DELK) *DELK)/
     >                             MAX (ABS (DELI) + ABS (DELK), EPS)

               DELI = WTI1 * DFACEJ (3, I1, K) + WTI2 * DFACEJ (3, I2,K)
               DELK = WTK1 * DFACEJ (3, I, K1) + WTK2 * DFACEJ (3, I,K2)

               DFACEJ (3, I, K) = (ABS (DELI) *DELI + ABS (DELK) *DELK)/
     >                             MAX (ABS (DELI) + ABS (DELK), EPS)
            END DO
         END DO

      ELSE IF (K1 .EQ. K2) THEN

C        K plane case:
C        -------------

         K = K1

C        Corner perturbations:

         DO J = J1, J2, J2 - J1
            DO I = I1, I2, I2 - I1
               DFACEK (1, I, J) = X (I, J, K) - X0 (I, J, K)
               DFACEK (2, I, J) = Y (I, J, K) - Y0 (I, J, K)
               DFACEK (3, I, J) = Z (I, J, K) - Z0 (I, J, K)
            END DO
         END DO

C        Intermediate edge perturbations:

         DO I = I1, I2, I2 - I1
            S0J1 = S0 (I, J1, K, 2)
            SRANGEJ = ONE / (S0 (I, J2, K, 2) - S0J1)
            DO J = J1 + 1, J2 - 1
               WTJ2 = (S0 (I, J, K, 2) - S0J1) * SRANGEJ
               WTJ1 = ONE - WTJ2
               DFACEK (1, I, J) = WTJ1 * DFACEK (1, I, J1) +
     >                            WTJ2 * DFACEK (1, I, J2)
               DFACEK (2, I, J) = WTJ1 * DFACEK (2, I, J1) +
     >                            WTJ2 * DFACEK (2, I, J2)
               DFACEK (3, I, J) = WTJ1 * DFACEK (3, I, J1) +
     >                            WTJ2 * DFACEK (3, I, J2)
            END DO
         END DO

         DO J = J1, J2, J2 - J1
            S0I1 = S0 (I1, J, K, 1)
            SRANGEI = ONE / (S0 (I2, J, K, 1) - S0I1)
            DO I = I1 + 1, I2 - 1
               WTI2 = (S0 (I, J, K, 1) - S0I1) * SRANGEI
               WTI1 = ONE - WTI2
               DFACEK (1, I, J) = WTI1 * DFACEK (1, I1, J) +
     >                            WTI2 * DFACEK (1, I2, J)
               DFACEK (2, I, J) = WTI1 * DFACEK (2, I1, J) +
     >                            WTI2 * DFACEK (2, I2, J)
               DFACEK (3, I, J) = WTI1 * DFACEK (3, I1, J) +
     >                            WTI2 * DFACEK (3, I2, J)
            END DO
         END DO

C        Intermediate perturbations of interior points:

         DO J = J1 + 1, J2 - 1
            DO I = I1 + 1, I2 - 1
               WTI2 = (S0 (I,  J, K, 1) - S0 (I1, J, K, 1)) /
     >                (S0 (I2, J, K, 1) - S0 (I1, J, K, 1))
               WTI1 = ONE - WTI2
               WTJ2 = (S0 (I, J,  K, 2) - S0 (I, J1, K, 2)) /
     >                (S0 (I, J2, K, 2) - S0 (I, J1, K, 2))
               WTJ1 = ONE - WTJ2

               DELI = WTI1 * DFACEK (1, I1, J) + WTI2 * DFACEK (1, I2,J)
               DELJ = WTJ1 * DFACEK (1, I, J1) + WTJ2 * DFACEK (1, I,J2)

               DFACEK (1, I, J) = (ABS (DELI) *DELI + ABS (DELJ) *DELJ)/
     >                             MAX (ABS (DELI) + ABS (DELJ), EPS)

               DELI = WTI1 * DFACEK (2, I1, J) + WTI2 * DFACEK (2, I2,J)
               DELJ = WTJ1 * DFACEK (2, I, J1) + WTJ2 * DFACEK (2, I,J2)

               DFACEK (2, I, J) = (ABS (DELI) *DELI + ABS (DELJ) *DELJ)/
     >                             MAX (ABS (DELI) + ABS (DELJ), EPS)

               DELI = WTI1 * DFACEK (3, I1, J) + WTI2 * DFACEK (3, I2,J)
               DELJ = WTJ1 * DFACEK (3, I, J1) + WTJ2 * DFACEK (3, I,J2)

               DFACEK (3, I, J) = (ABS (DELI) *DELI + ABS (DELJ) *DELJ)/
     >                             MAX (ABS (DELI) + ABS (DELJ), EPS)
            END DO
         END DO

      END IF

      RETURN
      END
C+------------------------------------------------------------------------------
C
      SUBROUTINE FMINRC (AX, BX, XMIN, FMIN, TOL, NUMFUN, CALLER,
     >                   LUNOUT, ISTAT)
C
C     Oneliner:
C
C     Function MINimization (1-D; no derivatives; Reverse Communication)
C     -        ---                                -       -
C
C     Description:
C
C        FMINRC estimates the point where a function of 1 variable has a
C     minimum value on the interval (AX, BX).  The corresponding function
C     value is also returned upon convergence.  Function derivatives are
C     not required.
C
C        This is an adaptation of FMIN77, which is itself a restructured,
C     subroutine version of Brent's FMIN function routine, referenced below.
C     FMINRC avoids the problems of passing as arguments the names of modules
C     with assumed calling sequences (such as FMIN77's FUN (X)) at the expense
C     of having to return to the calling program for every function evaluation. 
C     Such reverse communication is usually preferable to the use of common
C     blocks forced on nontrivial function routines by the likes of FMIN.
C
C
C     Arguments:
C
C     Name   Type   I/O/S   Description
C
C     AX      R     I       Left endpoint of initial interval.
C
C     BX      R     I       Right endpoint of initial interval.
C
C     XMIN    R       O     XMIN is output with the point at which the
C                           next function evaluation is required, until
C                           termination, when XMIN is the best estimate
C                           of the location of the minimum.
C
C     FMIN    R     I/O     After the initial call, FMIN should be input
C                           with the function value corresponding to XMIN.
C                           Upon termination, FMIN is output with the lowest
C                           function value found, corresponding to the
C                           final value of XMIN.
C
C     TOL     R     I       Desired length of the interval of uncertainty
C                           of the final result.  TOL >= 0.
C
C     NUMFUN  I     I/O     At the start of a minimization, NUMFUN should
C                           be input with the maximum number of function
C                           evaluations to be permitted during this
C                           minimization.  On output, NUMFUN contains
C                           the number of function evaluations so far.
C
C     CALLER  C*(*) I       Name of the calling program, printed if LUNOUT > 0
C                           and also used to warn of the (rare) possibility of
C                           nested minimizations, which are not supported
C                           by this version of the FMIN algorithm.  No more
C                           than 8 characters are significant.
C
C     LUNOUT  I     I       Logical unit number for output of the iteration
C                           history.  Use LUNOUT < 0 if you don't want the
C                           evaluations printed.  Any diagnostics go to
C                           the unit defined by ABS (LUNOUT).
C
C     ISTAT   I     I/O     Status flag used as follows:
C                           ISTAT = +2 on input means this call initiates
C                                      a new minimization.
C                           ISTAT = +1 on return means the calling program
C                                      should evaluate the function at the
C                                      current value of XMIN, and call FMINRC
C                                      again with this value in argument FMIN
C                                      and ISTAT = +1 still.
C                           ISTAT =  0 on output means a minimum has been
C                                      found to the specified tolerance.
C                           ISTAT = -1 on output means the limit on the
C                                      number of function evaluations
C                                      has been reached but the convergence
C                                      criteria have not been satisfied.
C                                      XMIN, FMIN are the best values found.
C                           ISTAT = -2 means AX > BX.  This is checked only
C                                      at the start of a minimization.
C                           ISTAT = -3 means the input value of NUMFUN on
C                                      the first call is either non-positive
C                                      or absurdly large.
C                           ISTAT = -4 means the calling program's name
C                                      has changed in the middle of a
C                                      minimization.  Such nested use
C                                      of FMINRC is not supported because
C                                      local variables are reused until a
C                                      minimization is complete.
C
C     Usage:
C
C        ISTAT = 2          ! Initialize the minimization
C        AX = ...
C        BX = ...
C        NUMFUN = ...       ! Max. no. of fn. evals. allowed
C                           <Define TOL and SUBNAME as well.>
C     10 CONTINUE
C
C           CALL FMINRC (AX, BX, XMIN, FMIN, TOL, NUMFUN, SUBNAME,
C       >                LUNOUT, ISTAT)
C
C           IF (ISTAT .LT. -1) THEN      ! Fatal error
C
C              <Handle it>
C
C           ELSE IF (ISTAT .LT. 0) THEN  ! Iteration limit reached
C
C              <Unconverged, but XMIN, FMIN may still be usable>
C
C           ELSE IF (ISTAT .GT. 0) THEN  ! Evaluate the function
C
C              CALL fun (XMIN, FMIN)     ! Or whatever
C              GO TO 10
C
C           ELSE ! ISTAT = 0 (success).  Ensure that FMIN is at XMIN.
C
C              CALL fun (XMIN, FMIN)     ! Or whatever
C
C           END IF
C
C 
C     Notes on the FMIN algorithm:
C
C        The method used is a combination of Golden Section search and
C     successive parabolic interpolation.  Convergence is never much slower
C     than that for a Fibonacci search.  If the function has a continuous
C     second derivative which is positive at the minimum (which is not at
C     AX or BX), then convergence is superlinear, and usually of the order
C     of about 1.324.
C
C        The function is never evaluated at two points closer together
C     than RTEPS * ABS (FMIN) + (TOL/3), where RTEPS is approximately the
C     square root of the relative machine precision.  If the function is
C     unimodal and the computed values are always unimodal when separated
C     by at least RTEPS * ABS (X) + (TOL/3), then the FMIN algorithm
C     approximates the abcissa of the global minimum of the function on the
C     interval (AX, BX) with an error less than 3 * RTEPS * ABS (FMIN) + TOL.
C     If the function is not unimodal, then FMIN may approximate a local, but
C     perhaps non-global, minimum to the same accuracy.
C
C        FMIN77 is a modified version of the ALGOL 60 procedure LOCALMIN
C     in Richard Brent, Algorithms for Minimization without Derivatives,
C     Prentice-Hall, Inc. (1973).  The FORTRAN 77 translation started from
C     the FORTRAN 66 version, FMIN, in Forsythe, Malcolm, and Moler but
C     also made use of the original ALGOL.
C
C        FMINRC is a minimal modification which takes advantage of the SAVE
C     feature of FORTRAN 77 to ensure that local variables are preserved
C     between related calls.
C
C
C     Further implementation notes:
C
C        (1)  IMPLICIT NONE is non-standard.
C
C        (2)  A more definitive adaptation permitting nested minimizations
C             would require the local variables to be transferred to and
C             from an additional work-space array argument.
C
C     History:
C
C     12/31/83  R.A.Kennelly  FMIN restructured and (partially) cleaned up
C               NASA Ames     as FMIN77.  Notes on unanswered questions are
C                             sprinkled throughout the code.
C
C     04/10/92  D.A.Saunders  FMINRC adapted from FMIN77 to avoid the FUN
C         Sterling Software/  argument problems.  Added the argument for
C               NASA Ames     counting/limiting the no. of function evals.
C                             RTEPS is now calculated during initialization.
C
C     02/28/93  DAS           Ensuring that FMIN is at XMIN was not shown
C                             in the code example above.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      REAL
     >   AX, BX, XMIN, FMIN, TOL
      INTEGER
     >   NUMFUN, LUNOUT, ISTAT
      CHARACTER
     >   CALLER * (*)

C     Local constants:

C     RATIO is ~(3 - SQRT (5)) / 2, the squared inverse of the Golden Ratio.

      REAL
     >   ZERO, ONE, TWO, THREE, HALF, THIRD, RATIO, SMALL
      PARAMETER
     >  (ZERO = 0.0E+0, ONE = 1.0E+0, TWO = 2.0E+0, THREE = 3.0E+0,
     >   HALF = ONE / TWO, THIRD = ONE / THREE,
     >   RATIO = .381966011250105E+0, SMALL = 1.E+0 / 2048.E+0)

C     Local variables:

      REAL
     >   A, B, DELTA, E, F, FU, FV, FW, FX, MIDDLE, P, Q, R, RTEPS,
     >   TOL1, TOL2, U, V, W, X
      INTEGER
     >   LENNAM, LUNERR, MAXFUN
      LOGICAL
     >   GOLDEN, TEST
      CHARACTER
     >   NAME * 8

C     System functions:

      INTRINSIC
     >   ABS, SIGN, SQRT

      SAVE         ! Vital for repeated calls during one minimization


C     Execution:
C     ----------

CRAK
CRAK  How about checking all inputs, e.g. TOL, etc.?
CRAK

C     Notation:  At the start of a cycle, we have
C
C        A, B   Left and right endpoints of current interval of uncertainty
C        X      Point at which the function is smallest, or the most recent
C               point if there was a tie
C        U      Last point at which the function was evaluated
C        W      Second best point thus far
C        V      Third best point thus far
C        DELTA  Proposed offset from latest point (U = X + DELTA)
C        E      A measure of the previous step size (precisely old DELTA when
C               a parabolic step was taken)


      IF (ISTAT .GT. 1) THEN

C        Initialize the iteration.
C        First, compute the square root of the relative machine precision:

         RTEPS = SMALL
    5    CONTINUE
            RTEPS = RTEPS / TWO
            TOL1 = ONE + RTEPS
            IF (TOL1 .GT. ONE)
     >   GO TO 5

         RTEPS = SQRT (RTEPS)

C        Initialize a new minimization.
C        ------------------------------
               
         LUNERR = ABS (LUNOUT)
         LENNAM = MIN (LEN (CALLER), LEN (NAME))
         NAME = CALLER
         IF (LUNOUT .GT. 0) WRITE (LUNOUT, 1010) NAME

         IF (AX .GT. BX) THEN
            ISTAT = -2
            WRITE (LUNERR, 1030) AX, BX

         ELSE IF (NUMFUN .LE. 1 .OR. NUMFUN .GT. 200) THEN
            ISTAT = -3
            WRITE (LUNERR, 1040) NUMFUN

         ELSE
            ISTAT = 1
            MAXFUN = NUMFUN
            NUMFUN = 0
            A  = AX
            B  = BX
            E  = ZERO
            X  = A + RATIO * (B - A)
            W  = X
            V  = X
            XMIN = X

C           Return early for the first function evaluation.
         END IF

         GO TO 99


      ELSE IF (NUMFUN .EQ. 0) THEN

C        Second call to FMINRC after the first function evaluation.

         NUMFUN = 1
         FX = FMIN
         FW = FX
         FV = FX
         IF (LUNOUT .GT. 0) WRITE (LUNOUT, 1020) NUMFUN, X, FX

         TEST = .TRUE.
      END IF


      IF (CALLER (1 : LENNAM) .NE. NAME (1 : LENNAM)) THEN  ! Fatal error
         ISTAT = -4
         WRITE (LUNERR, 1050) CALLER
         GO TO 99
      END IF


C     Start of a regular iteration.
C     -----------------------------

   10 CONTINUE

      IF (TEST) THEN   ! "TEST" was introduced to allow jumping out for a
C                      ! function value then continuing on the next entry.

C        Test the two stopping criteria.
C        -------------------------------

C        The distance from X to the midpoint of the interval of uncertainty
C        plus half the width of the interval must be less than twice the
C        tolerance.

         MIDDLE = HALF * (A + B)
         TOL1 = RTEPS * ABS (X) + THIRD * TOL
         TOL2 = TWO * TOL1

CRAK     Note Brent doesn't use THIRD above (added in FM&M) ?
CRAK     Equivalent to:  X-A .LE. TOL2  when X > MIDDLE   ...and...
CRAK                     B-X .LE. TOL2  when X < MIDDLE
CRAK     This isn't the same as the error estimate in the header. (It's
CRAK     actually a little tighter.)

         IF (ABS (X - MIDDLE) .LE. (TOL2 - HALF * (B - A))) THEN ! Success
            ISTAT = 0
            GO TO 20
         ELSE IF (NUMFUN .EQ. MAXFUN) THEN ! Too many function evaluations
            ISTAT = -1
            GO TO 20
         END IF


C        Compute a trial step.
C        ---------------------

C        Use parabolic interpolation if possible.

         GOLDEN = (NUMFUN .LT. 3 .OR. ABS (E) .LE. TOL1)
         IF (.NOT.GOLDEN) THEN


C           Fit a parabola.
C           ---------------

            R = (X - W) * (FX - FV)
            Q = (X - V) * (FX - FW)
            P = (X - V) * Q - (X - W) * R
            Q = TWO * (Q - R)

            IF (Q .GT. ZERO) THEN
               P = -P
            ELSE
               Q = -Q
            END IF

C           Is the parabola acceptable ?  Require that the proposed step be:
C
C           (a) Less than half of the second-to-last one (measured
C               roughly by E in case the last step was by Golden Section),
C               *** I don't like this one !!!  No good reason for it that
C               I can see, and not discussed in Brent except a mention that
C               parabolic steps should be small near a minimum, p.74  *** 
C
C           (b) Smaller than the steps to the left or right endpoints of
C               the current interval of uncertainty.  *** Note that (b)
C               also tosses out the parabolic choice on the second step,
C               as must happen since the fit is ill-defined. ***  NUMFUN
C               makes this redundant, though. ***

            IF (ABS (P) .LT. ABS (HALF * Q * E) .AND.
     >         P .GT. Q * (A - X) .AND.
     >         P .LT. Q * (B - X)) THEN


C              Take a parabolic interpolation step.
C              ------------------------------------

C              Here, E is saved as the step taken last iteration.

CRAK           This is suspect, since when parab. follows golden, we end up
CRAK           using E in the first test, then saving corresponding DELTA as
CRAK           next E, thus lagging SAME THING in for TWO iterations !?

               E = DELTA
               DELTA = P / Q
               U = X + DELTA

C              The function must not be evaluated too close to the endpoints.

               IF ((U - A) .LT. TOL2 .OR. (B - U) .LT. TOL2)
     >            DELTA = SIGN (TOL1, MIDDLE - X)

            ELSE

C              Proposed parabolic point failed.

               GOLDEN = .TRUE.
            END IF
         END IF

         IF (GOLDEN) THEN


C           Use a Golden-Section step.
C           --------------------------

C           Here, E is the step from X to the more distant endpoint.

            IF (X .GE. MIDDLE) THEN
               E = A - X
            ELSE
               E = B - X
            END IF
            DELTA = RATIO * E
         END IF

C        The function must not be evaluated too close to X.

         IF (ABS (DELTA) .LT. TOL1) THEN
            U = X + SIGN (TOL1, DELTA)
         ELSE
            U = X + DELTA
         END IF


         XMIN = U
         TEST = .FALSE.

C        Return to caller to evaluate the function at the chosen point.
C        --------------------------------------------------------------


      ELSE

C        Pick up where we left off, with the new function value in hand.
C        ---------------------------------------------------------------

         FU = FMIN
         NUMFUN = NUMFUN + 1
         IF (LUNOUT .GT. 0) WRITE (LUNOUT, 1020) NUMFUN, U, FU

         IF (FU .LE. FX) THEN

C           The new value is an improvement.

            IF (U .LT. X) THEN
               B = X
            ELSE
               A = X
            END IF

            V  = W
            FV = FW
            W  = X
            FW = FX
            X  = U
            FX = FU

         ELSE

C           The new value is not as good as previous best.

            IF (U .LT. X) THEN
               A = U
            ELSE
               B = U
            END IF

            IF (FU .LE. FW .OR. W .EQ. X) THEN
               V  = W
               FV = FW
               W  = U
               FW = FU
            ELSE IF (FU .LE .FV .OR. V .EQ. X .OR. V .EQ. W) THEN
               V  = U
               FV = FU
            ELSE

C              What is left over here ???  Looks like:
C
C                 FU > FX, FW, FV  &  W >< X  &  V >< X, W
C
C              This DOES happen !  Why is it not discussed ?

            END IF

         END IF

         TEST = .TRUE.          ! Go test for convergence.
         GO TO 10

      END IF

      GO TO 99


C     Terminate.
C     ----------

   20 XMIN = X
      FMIN = FX


   99 RETURN

C     Formats:

 1010 FORMAT (/, ' FMINRC:  Iteration history for calls by ', A, //
     >        ' Iter', 23X, 'X', 21X, 'F(X)')
 1020 FORMAT (1X, I3, 1P, 2E25.15)
 1030 FORMAT (/, ' FMINRC:  AX > BX is illegal: ', 1P, 2E24.16)
 1040 FORMAT (/, ' FMINRC:  Suspicious input for (max.) NUMFUN: ', I10)
 1050 FORMAT (/, ' FMINRC:  Cannot start a new minimization before the',
     >       ' previous one is done.', /, ' Offending module: ', A)

      END
C+----------------------------------------------------------------------
C
      SUBROUTINE INTERVAL (NX, X, XFIND, ARROW, LEFT)
C
C     One-liner: Interpolation search for interval containing a point.
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        Written primarily for interval-based interpolations such as
C     piecewise linear or cubic spline, INTERVAL performs a search to
C     locate the best interval for evaluating the interpolant at a
C     given point. The normal case returns the "left-hand" endpoint of
C     the interval bracketing the point, but for the out-of-range cases
C     below or above the range of the knots, the interval to be used is
C     the first or last. The array of knots must be monotonic, either
C     increasing or decreasing. Diagrammatically, LEFT is returned as
C     shown below for the normal case (no extrapolation):
C
C          X (1)  ...   X (LEFT)   X (LEFT+1)   ...      X (NX)
C                                ^
C                              XFIND
C
C     And for extrapolation:
C
C                     X (LEFT = 1)  ...   X (NX)
C             ^
C           XFIND
C
C     or,
C                X (1)  ...   X (LEFT = NX-1)    X (NX)
C                                                           ^
C                                                         XFIND
C
C     If the point to be bracketed (XFIND) matches one of the knots, the
C     index of that knot is returned as LEFT, i.e., the condition for a
C     bracket of an interior point is:
C
C        X (LEFT) <= XFIND < X (LEFT+1)  if  ARROW = +1.0,  or
C        X (LEFT) >= XFIND > X (LEFT+1)  if  ARROW = -1.0.
C
C        This is a low-level routine with minimal error checking. The
C     calling program is assumed to have verified the following:
C
C     (1)  NX >= 2
C     (2)  X strictly monotonic
C     (3)  ARROW = +1.0 or -1.0
C
C     Subroutine PROTECT is available from the author for easily checking
C     conditions (2) and (3). LEFT is verified on input, but efficiency in
C     loops will benefit from passing the best estimate available, usually
C     just the result of the last call.
C
C        INTERVAL was originally written for use with CSEVAL and TABLE1.
C     The interpolation search was adapted from ideas in Sedgewick's book
C     referenced below.
C
C     Arguments:
C     ----------
C
C     Name  Dimension  Type  I/O/S  Description
C     NX                I    I      Number of points in array X; must
C                                   be >= 2 (no check performed).
C
C     X        NX       R    I      Array of points defining the set
C                                   of intervals to be examined. Only
C                                   the first NX-1 points are required.
C
C     XFIND             R    I      The point for which a bracketing
C                                   interval is sought.
C
C     ARROW             R    I      Monotonicity indicator for input
C                                   array X:
C                                     -1.0  strictly decreasing
C                                      0.0  NOT ALLOWED!
C                                     +1.0  strictly increasing
C                                   Supplied by the calling routine for
C                                   reasons of speed (not checked).
C
C     LEFT              I    I/O    Input: guessed index of left-hand
C                                   endpoint of the interval containing
C                                   the specified point.
C
C                                   Output: index of the largest array
C                                   value <= specified point (if ARROW=+1.0).
C                                   Special case for data out of range:
C                                   return left endpoint of closest interval.
C                                   Thus, LEFT = 1 for XFIND < X (2), and
C                                   LEFT = NX-1 for XFIND >= X (NX-1).
C                                   (If ARROW=-1.0, reverse the inequalities.)
C
C     Environment:  Digital VAX-11/780, VMS FORTRAN
C     ------------  Apple Macintosh, Absoft MacFORTRAN/020 v2.3
C
C     Notes:
C     ------
C
C     (1)  IMPLICIT NONE and eight character symbols are not (yet) standard.
C
C     (2)  In speed-critical applications, it might be a good idea to build
C          this algorithm in-line since it is typically called many times
C          from within a loop. Another potential speed-up is removal of the
C          ARROW multiplies, which restricts the method to increasing data.
C          So far, the simplicity of separating out the messy search details
C          and the generality of bi-directional searching have outweighed
C          the modest speed penalty incurred.
C
C     Bibliography:
C     -------------
C
C     (1) Sedgewick, R.  Algorithms.  Reading: Addison-Wesley, 1983.
C            (Chap. 14)
C
C     Author:  Robert Kennelly and David Saunders, Sterling Federal Systems
C     -------
C
C     Development history:
C     --------------------
C
C     20 Oct. 1987    RAK    Interpolation search adapted (with mods.
C                            for bidirectional search and some minor
C                            repair) from CSEVAL (RAK) and TABLE1 (DAS).
C     08 Aug. 1988    DAS    Clarified descriptions of bracketing, where
C                            the inequalities depend upon ARROW.
C
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     &   ONE
      PARAMETER
     &  (ONE = 1.0E+0)

C     Arguments.

      INTEGER
     &   LEFT, NX
      REAL
     &   ARROW, X (NX), XFIND

C     Local variables.

      INTEGER
     &   LENGTH, NXLESS1, RIGHT, TRIAL
      REAL
     &   XBYARROW

C     Execution.
C     ----------

      XBYARROW = XFIND * ARROW

C     Simplify things by disposing of two important special cases so that
C     X (LEFT) and X (RIGHT) can really bracket XFIND. As a by-product,
C     this also takes care of the NX = 2, 3 cases.

      NXLESS1 = NX - 1

      IF (XBYARROW .GE. X (NXLESS1) * ARROW) THEN
         LEFT = NXLESS1
         GO TO 990
      ELSE IF (XBYARROW .LT. X (2) * ARROW) THEN
         LEFT = 1
         GO TO 990
      END IF

C      ---------------------------------
C     |                                 |
C     |   X (2) <= XFIND < X (NX - 1)   |
C     |            - or -               |
C     |   X (2) > XFIND >= X (NX - 1)   |
C     |                                 |
C     |   NX > 3                        |
C     |                                 |
C      ---------------------------------

C     Adjust the pointers. We hope that the calling routine has provided
C     a reasonable guess (since it's probably working on an ordered array
C     of points to evaluate), but check anyway.

      LEFT = MIN (MAX (2, LEFT), NX - 2)

      IF (XBYARROW .GE. X (LEFT) * ARROW) THEN
         IF (XBYARROW .LT. X (LEFT + 1) * ARROW) THEN

C           XFIND is in the original guessed-at interval.

            GO TO 990
         ELSE

C           We'll look farther to the right. Evidently LEFT was < NX - 2.

            RIGHT = NXLESS1
            LEFT  = LEFT + 1
         END IF
      ELSE

C        Look to the left of the guess. Evidently LEFT was > 2.

         RIGHT = LEFT
         LEFT  = 2
      END IF

C      ----------------------------------
C     |                                  |
C     |   2 <= LEFT < RIGHT <= NX - 1    |
C     |                                  |
C      ----------------------------------

C     The interval length must decrease each time through - terminate
C     when the correct interval is found or when the interval length
C     cannot be decreased.

   10 CONTINUE
         LENGTH = RIGHT - LEFT
         IF (LENGTH .GT. 1) THEN

C           The trial value is a "linear" estimate of the left-hand endpoint
C           of the interval bracketing the target XFIND, with protection
C           against round-off (which can affect convergence).

            TRIAL = MIN (RIGHT - 1, LEFT + MAX (0, INT (REAL (LENGTH) *
     &         (XFIND - X (LEFT)) / (X (RIGHT) - X (LEFT)))))

C            ------------------------------------------
C           |                                          |
C           |   2 <= LEFT <= TRIAL < RIGHT <= NX - 1   |
C           |                                          |
C            ------------------------------------------

C           Adjust pointers. Increase LEFT or decrease RIGHT until done.

            IF (XBYARROW .GE. X (TRIAL + 1) * ARROW) THEN
               LEFT  = TRIAL + 1
            ELSE IF (XBYARROW .LT. X (TRIAL) * ARROW) THEN
               RIGHT = TRIAL
            ELSE

C              We're done: XFIND is in the interval [X (TRIAL), X (TRIAL+1)).

               LEFT  = TRIAL
               GO TO 990
            END IF
            GO TO 10

         END IF

C     Termination.
C     ------------

  990 CONTINUE
      RETURN
      END
C+------------------------------------------------------------------------------
C
      SUBROUTINE INTSEC5 (IDIM, JDIM, I1, I2, J1, J2, XS, YS, ZS,
     >                    US, VS, IS, JS, EPS, K1, K2, XC, YC, ZC,
     >                    TC, KC, TTOTAL, METHOD, CLOSED,
     >                    TINT, UINT, VINT, XINT, YINT, ZINT,
     >                    LUNOUT, IER)
C
C ONE-LINER: INTerSECtion of a curve and a surface in 3-space (bilinear version)
C
C PURPOSE:
C
C        INTSEC5 determines where a curve in 3-space meets a surface.  The
C     curve and surface are defined as discrete points, and extrapolation of
C     the curve (but not the surface) is permitted.
C
C        This is an adaptation of INTSEC4, using parametric bilinear
C     interpolation in place of bicubic, with analytic derivatives from
C     the formulation of BILINT for variation with respect to p and q.
C
C        This version also treats the curve as piecewise linear.  See INTSEC4
C     for further details.
C
C HISTORY:
C     09/23/96  DAS  Adaptation of INTSEC4.
C     09/24/97   "   Replaced PLSCURVE with PLSCRV3D, specifying 'L'
C                    (appropriate for wing surface grids) rather than
C                    adding a METHOD argument.
C     11/28/97   "   Generalized WBINTR forced a METHOD argument.
C     02/06/98   "   Convergence test now uses relative errors in each of x,y,z.
C     08/06/99   "   PBILINT now has analytic derivatives, allowing removal of
C                    the finite differencing originally implemented here.
C
C AUTHOR: David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   IDIM, JDIM,           ! Max. # pts. provided for in the surface arrays
     >   I1, I2, J1, J2        ! Grid index range eligible for searching;
                               ! 1 <= I1 < I1 + 1 < I2 <= IDIM, etc.

      REAL, INTENT (IN), DIMENSION (IDIM, JDIM) ::
     >   XS, YS, ZS            ! Surface grid coordinates

      REAL, INTENT (INOUT), DIMENSION (IDIM, JDIM) ::
     >   US, VS                ! Parametric variables at the surface grid pts.:
                               ! calculated internally if requested via
                               ! US (I1, J1) < 0 on input, else input from a
                               ! prior call to INTSEC5 or PARAM2D (or PLBICUBE)

      INTEGER, INTENT (INOUT) ::
     >   IS, JS                ! Indices of the surface cell containing the
                               ! intersection point:
                               ! on input:  the values from a previous
                               !            call, or I1, J1, if no better
                               !            estimate is known;
                               ! on output: the "lower left" indices of
                               !            the cell (unless IER > 0);
                               ! I1 <= IS < I2 input & output; likewise for JS

      REAL, INTENT (IN) ::
     >   EPS                   ! Tolerance for convergence test here and in
                               ! PBILINT/RIPPLE2D, relative to 1.;
                               ! try 5.*machine eps for 32-bit arithmetic,
                               ! or 1.E-7 should be OK for 64-bit "

      INTEGER, INTENT (IN) ::
     >   K1, K2                ! Index range for the curve arrays: 1 <= K1 < K2

      REAL, INTENT (IN), DIMENSION (K2) ::
     >   XC, YC, ZC            ! Curve coordinates

      REAL, INTENT (INOUT) ::
     >   TC (K2)               ! Parameterization of the curve (normalized
                               ! cumulative chord lengths) - see TTOTAL.
                               ! TC (K1) = 0 and TC (K2) = 1.  TC is an
                               ! argument to avoid likely recalculation.

      INTEGER, INTENT (INOUT) ::
     >   KC                    ! Index of curve interval containing the
                               ! intersection point.  Use K1 if no better
                               ! estimate is known on input.  On output,
                               ! K1 <= KC < K2 by INTERVAL's convention.

      REAL, INTENT (INOUT) ::
     >   TTOTAL                ! Total cumulative chord length of the
                               ! curve (K1 : K2); accompanies TC (*)
                               ! because the normalization loses it.
                               ! Input TTOTAL = 0. if the calling program
                               ! is not supplying TC (*) and TTOTAL, in
                               ! which case INTSEC5 calculates them.

      CHARACTER, INTENT (IN) ::
     >   METHOD * 1            ! Curve fit method; use 'L', 'M', or 'B'
                               ! for piecewise linear, monotonic (tight)
                               ! cubic, or "Bessel" (smooth) cubic resp.

      LOGICAL, INTENT  (IN) ::
     >   CLOSED                ! .TRUE. means the curve is closed, with
                               ! input end pts. matching - see PLSCRV3D.

      REAL, INTENT (INOUT) ::
     >   TINT, UINT, VINT      ! Estimated point of intersection in the
                               ! parameter space.  On input, enter 0. if
                               ! unknown, in which case TC (KC) and
                               ! US/VS (IS, JS) are used as starting
                               ! guesses, else enter the values from a
                               ! previous call.  All are in [0, 1].

      REAL, INTENT (OUT) ::
     >   XINT, YINT, ZINT      ! Estimated point of intersection in real space

      INTEGER, INTENT (IN) ::
     >   LUNOUT                ! Logical unit for displaying iterations,
                               ! which are suppressed if LUNOUT < 0.
                               ! |LUNOUT| is used for error messages.

      INTEGER, INTENT (OUT) ::
     >   IER                   ! 0 means no problem was encountered;
                               ! 1 means trouble in the surface interpolation
                               !   but results may still be usable;
                               ! 2 means the safeguarded iteration failed to
                               !   satisfy EPS but did satisfy 10*EPS;
                               ! 3 means the iteration limit was reached without
                               !   (close) convergence;
                               ! 4 means the step-halving failed - somehow the
                               !   iteration diverged;
                               ! 5 means the linearized system was singular:
                               !   the surface and curve could be parallel.

C     Procedures:

      EXTERNAL
     >   CHORDS3D,             ! Cumulative chord lengths for a line
     >   PARAM2D,              ! Cumulative chord lengths for a surface mesh
     >   PBILINT,              ! Parametric bilinear interpolation
     >   PLSCRV3D,             ! Local cubic spline utility for a 3-space curve
     >   LUSOLVE               ! Square system solver (LU decomposition)

C-------------------------------------------------------------------------------

C     Local constants:

      INTEGER, PARAMETER ::
     >   ITMAX   = 10      ! Outer iteration limit.  (Halving the step
                           ! is the inner iteration.)
      REAL, PARAMETER ::
     >   HALF    = 0.5,
     >   ONE     = 1.0,
     >   STEPMIN = 0.001,  ! Limit on repeated step halvings
     >   ZERO    = 0.0

C     Local variables:

      INTEGER
     >   I, ITER, LUNERR, NPTS

      REAL
     >   A (3, 3), DT (3), F (3), T (3), TLAST (3),
     >   ABSF, ALPHA, FNORM, FNORM0, P, Q, TOLER, XINTC, YINTC, ZINTC,
     >   XINTS, YINTS, ZINTS, XSCALE, YSCALE, ZSCALE

      LOGICAL
     >   NEW


C     Initialization:
C     ---------------

      LUNERR = ABS (LUNOUT)
      NPTS = K2 - K1 + 1
      NEW = .TRUE.     ! For the curve
      A (1, 3) = ZERO  ! -999. would suppress curve derivatives

C     Parameterize the surface?
C     -------------------------

      IF (US (I1, J1) < ZERO) THEN ! Must be a first call

         CALL PARAM2D (IDIM, JDIM, I1, I2, J1, J2, XS, YS, ZS, US, VS)

      END IF


C     Parameterize the curve?
C     -----------------------

      IF (TTOTAL == ZERO) THEN ! Specify normalized chord lengths

         CALL CHORDS3D (NPTS, XC (K1), YC (K1), ZC (K1), .TRUE., TTOTAL,
     >                  TC (K1))
      END IF


C     Safeguarded Newton iteration to solve the 3x3 nonlinear system:
C     ---------------------------------------------------------------

      ITER   = 0
      ALPHA  = ZERO      ! Step length initialized for print purposes only
      FNORM0 = 1.E+20    ! I.e., big to avoid iteration 0 test
      TOLER  = 3.* EPS   ! Used for || f || test; avoids divide by 3

      T (1) = UINT       ! Starting guesses for u, v, t
      T (2) = VINT
      T (3) = TINT
      IF (T (1) == ZERO) T (1) = US (IS, JS)
      IF (T (2) == ZERO) T (2) = VS (IS, JS)      
      IF (T (3) == ZERO) T (3) = TC (KC)

  300 CONTINUE

C        Evaluate the surface and its u/v derivatives at the current (u,v):

         CALL PBILINT (1, IDIM, JDIM, I1, I2, J1, J2, XS, YS, ZS,
     >                 US, VS, T (1), T (2), IS, JS, EPS, P, Q,
     >                 XINTS, YINTS, ZINTS, A (1, 1), A (1, 2), IER)

         IF (IER /= 0) WRITE (LUNERR, 1000) T

C        Evaluate the curve coordinates and derivatives at the current t:

         CALL PLSCRV3D (NPTS, XC (K1), YC (K1), ZC (K1), TC (K1),
     >                  METHOD, NEW, CLOSED, T (3), KC,
     >                  XINTC, YINTC, ZINTC, A (1, 3))

         IF (ITER == 0) THEN  ! Get the scaling from the initial curve point:
            NEW = .FALSE.
            XSCALE = MIN (ONE, ONE / MAX (EPS, ABS (XINTC)))
            YSCALE = MIN (ONE, ONE / MAX (EPS, ABS (YINTC)))
            ZSCALE = MIN (ONE, ONE / MAX (EPS, ABS (ZINTC)))
         END IF

C        Find the norm of the residual vector f, which should converge to ~0.
C        Set up the RHS vector for the system J dT = f in the process.

         F (1) = XINTS - XINTC
         F (2) = YINTS - YINTC
         F (3) = ZINTS - ZINTC
         FNORM = ABS (F (1)) * XSCALE + ABS (F (2)) * YSCALE +
     >           ABS (F (3)) * ZSCALE

C        Halve the step until || f || is reduced (except first time through).
C        See comment below about constraining u and v but not t.

         IF (FNORM >= FNORM0) THEN
            IF (ALPHA > STEPMIN) THEN
               ALPHA = HALF * ALPHA
               T (1) = MAX (ZERO, MIN (ONE, TLAST (1) - ALPHA * DT (1)))
               T (2) = MAX (ZERO, MIN (ONE, TLAST (2) - ALPHA * DT (2)))
               T (3) = TLAST (3) - ALPHA * DT (3)
               GO TO 300
            END IF
            GO TO 830
         END IF

         IF (LUNOUT > 0) THEN
            WRITE (LUNOUT, '(A, I3, A, 1P, E14.6, A, E9.2, A, 3E14.6)')
     >         ' INTSEC5:', ITER, '  ||f||:', FNORM, '  step:', ALPHA,
     >         '  t:', T
         END IF

         IF (FNORM < TOLER) GO TO 900 ! Success


         ITER = ITER + 1
         IF (ITER == ITMAX) GO TO 810

         FNORM0 = FNORM

C        Fix up the LHS matrix for the next iteration.  The first
C        two columns are already the X/Y/Z vectors of surface derivatives:

         A (1, 3) = -A (1, 3)
         A (2, 3) = -A (2, 3)
         A (3, 3) = -A (3, 3)

C        Solve  J dt = f; dt overwrites f.

         CALL LUSOLVE (3, 3, A, F, IER)

         IF (IER /= 0) GO TO 840

         DO I = 1, 3
            DT (I) = F (I)
            TLAST (I) = T (I)
            T (I) = T (I) - DT (I)
         END DO

C        Extrapolation by the curve is permitted:

         T (1) = MAX (ZERO, MIN (ONE, T (1)))
         T (2) = MAX (ZERO, MIN (ONE, T (2)))
         ALPHA = ONE

      GO TO 300                           ! Do another iteration


C     Error handling:
C     ---------------

  810 IER = 2
      WRITE (LUNERR, 1010) 'Iteration limit reached.'
      IF (FNORM < 10.* TOLER) GO TO 900

      IER = 3
      GO TO 999

  830 IER = 4
      WRITE (LUNERR, 1010) 'Step-halving iteration failed.'
      GO TO 999

  840 IER = 5
      WRITE (LUNERR, 1010) 'Singular system. Parallel?'
      GO TO 999
      
  900 CONTINUE
C     Wrap up, having converged or almost converged.
C     Trying to average the surface & curve results somehow is too awkward.

      UINT = T (1)
      VINT = T (2)
      TINT = T (3)
      XINT = XINTS
      YINT = YINTS
      ZINT = ZINTS

  999 RETURN

C     Formats:

 1000 FORMAT (/, ' INTSEC5:  PBILINT target out of range.  Proceeding',
     >        ' with u, v, t =', 3F10.6)
 1010 FORMAT (/, ' INTSEC5: ', A)

      END SUBROUTINE INTSEC5
C+------------------------------------------------------------------------------
C
      SUBROUTINE LUSOLVE (N, NDIM, A, B, IER)
C
C     ONE-LINER:  Square-system solution by LU decomposition
C
C     DESCRIPTION:
C
C        LUSOLVE solves the N*N system  A x = b  by Gaussian elimination
C     with partial pivoting.  The right-hand side  b  is overwritten by the
C     solution  x.  (Matrix  A  is also overwritten.)
C
C        For more than one right-hand side, use instead DECOMP (once) and
C     SOLVE once per RHS  b.  See also DECSLV for the single-right-hand-side
C     case:  it is slightly more efficient than LUSOLVE in its processing of
C     b  as column N+1 of the matrix, but may be less convenient.
C
C        LUSOLVE was adapted from DECSLV, itself a merger of the original
C     DECOMP and SOLVE (reference below).  Note that one version of DECOMP
C     and SOLVE provides an estimate of the matrix condition number, while
C     another version avoids this extra work for applications where cond (A)
C     of no interest.  All of these variations vary the first index in the
C     inner loops.
C
C     INPUT:
C        N   : Order of matrix  A;  N >= 2
C        NDIM: Row dim. of array  A  declared in the calling program
C        A   : N*N array containing matrix  A
C        B   : Right-hand-side N-vector  b
C
C     OUTPUT:
C        B   : Solution vector  x
C        IER : 0 means no problem was detected;
C              1 means  A  was found to be singular
C
C     REFERENCE:  Forsythe, Malcolm & Moler, 1972.
C
C     HISTORY:    07/09/93  DAS  Revision of DECSLV to avoid having to
C                                store b in an extra column of A.
C                                IER = 0 now means no error (as it should
C                                have from the start).
C
C     PROGRAMMING: David Saunders, Sterling Software/NASA Ames, Mt. View, CA
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER
     >   N, NDIM, IER
      REAL
     >   A (NDIM, N), B (N)

C     Local constants:

      REAL
     >   ONE, ZERO
      PARAMETER
     >  (ONE = 1.E+0, ZERO = 0.E+0)

C     Local variables:

      INTEGER
     >   I, J, K, M
      REAL
     >   T

C     Execution:

C     Perform the LU factorization, solving L y = b as we go:

      DO 60 K = 1, N - 1

C        Determine pivot element:

         M = K
         DO 20 I = K + 1, N
            IF (ABS (A (I, K)) .GT. ABS (A (M, K))) M = I
 20      CONTINUE

         T = A (M, K)
         A (M, K) = A (K, K)
         A (K, K) = T
         IF (T .EQ. ZERO) GO TO 90

C        Store the multipliers:

         T = -ONE / T
         DO 30 I = K + 1, N
            A (I, K) = A (I, K) * T
 30      CONTINUE

C        Apply the multipliers to current submatrix, including RHS:

         DO 50 J = K + 1, N
            T = A (M, J)
            A (M, J) = A (K, J)
            A (K, J) = T

            DO 40 I = K + 1, N
               A (I, J) = A (I, J) + A (I, K) * T
 40         CONTINUE

 50      CONTINUE

         T = B (M)
         B (M) = B (K)
         B (K) = T

         DO 55 I = K + 1, N
            B (I) = B (I) + A (I, K) * T
 55      CONTINUE

 60   CONTINUE

      IF (A (N, N) .EQ. ZERO) GO TO 90

C     Back substitution (solution of U x = y):

      DO 80 K = N, 2, -1
         T = B (K) / A (K, K)
         B (K) = T
         DO 70 I = 1, K - 1
            B (I) = B (I) - A (I, K) * T
 70      CONTINUE
 80   CONTINUE

      B (1) = B (1) / A (1, 1)
      IER = 0

      GO TO 99


 90   IER = 1  ! Matrix was singular

 99   RETURN
      END
C+------------------------------------------------------------------------------
C
      SUBROUTINE MORPH_LINE_3D (I1, I2, S1, S2, X0, Y0, Z0, X, Y, Z)
C
C  ONE-LINER: Variant of NULINE3D to control end points AND end-point slopes
C
C  DESCRIPTION:
C
C        MORPH_LINE_3D perturbs the interior points of a 3-space curve given
C     perturbed end points and desired end-point slopes.  More precisely, any
C     pair of points I1, I2 may be the controlling points; each point between
C     these two is moved according to an arc-length-based combination of the
C     distances between the original and new points corresponding to I1 and I2.
C
C        This first stage (that of NULINE3D) is further adjusted to obtain the
C     specified end-point slopes in a similar arc-length-based way which turns
C     out to be overdetermined:  N - 1 intervals in which 2-point derivative
C     estimates are specified but only N - 2 interior points to play with.
C
C        The new coordinates of the two end points should be input in the
C     desired output curve.
C
C  ENVIRONMENT:  Fortran 90
C
C  HISTORY:
C
C     03/16/98  DAS  Original NULINE3D, used here as the first stage.
C     02/18/04   "   MORPH_LINE_3D applies a linear least squares technique
C                    to achieve some degree of slope control at the cost of
C                    losing the original relative spacing.  Therefore ...
C     02/24/04   "   Interpolate back to the original relative spacings, and
C                    iterate to optimize some measure of goodness.
C     02/26/04   "   Replaced dense HDECOM/HSOLVE pair with specialized sparse
C                    pair BIDECOM/BISOLVE, and weighted the end-point slopes.
C                    Excessive weighting can lose slopes in the middle.
C                    Suppress the redistribution of the stage 1 result as
C                    hardly worth the extra arithmetic.
C
C  AUTHOR:  David Saunders, ELORET/NASA Ames Research Center, Moffet Field, CA.
C
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   I1, I2                ! Indices of given perturbed "end" points
      REAL, INTENT (IN) ::
     >   S1(3), S2(3)          ! Unit tangents desired at points I1, I2
      REAL, INTENT (IN) ::
     >   X0(*), Y0(*), Z0(*)   ! Original curve coordinates
      REAL, INTENT (INOUT) ::
     >   X(*), Y(*), Z(*)      ! Desired curve, with points I1, I2 input

C-------------------------------------------------------------------------------

C     Local constants:

      INTEGER, PARAMETER ::
     >   LUNOUT = -6           ! Suppress FMINRC iterations; diagnostics to +6
      REAL, PARAMETER ::
     >   ONE = 1., ZERO = 0.
      CHARACTER, PARAMETER ::
     >   METHOD * 1 = 'B',     ! Plain Hermite cubics of X, Y, Z vs. arc
     >   CALLER * 8 = 'MORPH_3D'
      LOGICAL, PARAMETER ::
     >   CLOSED = .FALSE.

C     Local variables:

      INTEGER
     >   I, IEVAL, ISTAT, J, LUNERR, M, N, NDATA, NUMFUN
      REAL, DIMENSION (I1:I2) ::
     >   ARC, ARC0, B1, B2, B3, C, D, S, U, W, XI, YI, ZI
      REAL
     >   AV, BV, DX1, DX2, DY1, DY2, DZ1, DZ2, OBJ, T, TOL, TLENGTH1,
     >   V, W1, W2, DERIVS(3)
      LOGICAL
     >   FIRST, NEW

C     Procedures:

      EXTERNAL
     >   BIDECOM,  ! Specialized linear least squares pair
     >   BISOLVE,
     >   FMINRC,   ! Reverse-communication 1-D minimizer
     >   PLSCRV3D  ! Parametric spline interpolation utility

C     Execution:

CCC   write (6, '(/, (a, 3f12.7))') ' s1: ', s1, ' s2: ', s2

      M = I2 - I1;  N = M - 1;  NDATA = M + 1
      DERIVS(1) = -999.  ! Suppresses derivative outputs from interpolations

C     First, apply the simple end-point-location algorithm of NULINE3D:

      ARC0(I1) = ZERO
      DO I = I1 + 1, I2
         ARC0(I) = ARC0(I - 1) + SQRT (
     >     (X0(I) - X0(I - 1)) ** 2 +(Y0(I) - Y0(I - 1)) ** 2 +
     >     (Z0(I) - Z0(I - 1)) ** 2)
      END DO

      T = ONE / ARC0(I2)
      ARC0(I1:I2) = ARC0(I1:I2) * T ! Normalize original arc lengths

      XI(I1) = X(I1)    ! End points of coordinates to be iterated on
      YI(I1) = Y(I1)
      ZI(I1) = Z(I1)
      XI(I2) = X(I2)
      YI(I2) = Y(I2)
      ZI(I2) = Z(I2)

      DX1 = X(I1) - X0(I1)
      DX2 = X(I2) - X0(I2)
      DY1 = Y(I1) - Y0(I1)
      DY2 = Y(I2) - Y0(I2)
      DZ1 = Z(I1) - Z0(I1)
      DZ2 = Z(I2) - Z0(I2)

      DO I = I1 + 1, I2 - 1
         W2 = ARC0(I)
         W1 = ONE - W2
         XI(I) = X0(I) + W1 * DX1 + W2 * DX2
         YI(I) = Y0(I) + W1 * DY1 + W2 * DY2
         ZI(I) = Z0(I) + W1 * DZ1 + W2 * DZ2
      END DO

C     This stage 1 result now needs to be adjusted to obtain the specified
C     end-point slopes.
C     The interim relative spacing is NOT necessarily the same as for the
C     original curve - only if it's a straight line, we believe.
C     Therefore, update the current relative spacing and total length:

      ARC(I1) = ZERO
      DO I = I1 + 1, I2
         ARC(I) = ARC(I - 1) + SQRT (
     >      (XI(I) - XI(I - 1)) ** 2 + (YI(I) - YI(I - 1)) ** 2 +
     >      (ZI(I) - ZI(I - 1)) ** 2)
      END DO

      TLENGTH1 = ARC(I2)
      T = ONE / TLENGTH1

      ARC(I1:I2) = ARC(I1:I2) * T ! Normalized interim arc lengths

CCC   write (6, '(a, f12.7)') ' Stage 1 total length:', tlength1
CCC   write (10, '(a)') ' Stage 1 arc lengths:'
CCC   write (10, '(i3, 2f12.8)') (i, arc0(i), arc(i), i = i1, i2)

C     Impose the original relative spacing precisely on the stage 1 result:
C     NO - this doesn't seem to improve things significantly.

CCC   NEW = .TRUE.
CCC   IEVAL = 1

CCC   DO I = I1 + 1, I2 - 1
CCC      CALL PLSCRV3D (NDATA, XI(I1), YI(I1), ZI(I1), ARC(I1),
CCC  >                  METHOD, NEW, CLOSED, ARC0(I), IEVAL,
CCC  >                  X(I), Y(I), Z(I), DERIVS)
CCC      NEW = .FALSE.
CCC   END DO

CCC   DO I = I1 + 1, I2 - 1 ! Can't do this till PLSCRV3D is done
CCC      XI(I) = X(I)
CCC      YI(I) = Y(I)
CCC      ZI(I) = Z(I)
CCC   END DO
 
C     Copying ARC0(*) here is very close to recalculating arcs:

CCC   DO I = I1 + 1, I2
CCC      ARC(I) = ARC0(I)
CCC   END DO

C     Begin stage 2, where we attempt to impose slopes in ALL intervals,
C     including the first and last interval.  This tends to conflict with
C     imposing the original relative spacing, but we recover that at the
C     end at the expense of giving up a little on the end-point slopes.
C     Two-point derivatives keep the overdetermined system very simple.

      V = ONE ! Reasonable estimate for the unknown multiplier of the
              ! stage 1 total length that gives the final total length
      T = V * TLENGTH1

      FIRST = .TRUE.
      ISTAT = 2        ! Initialize the minimization
      TOL = MAX (EPSILON (TOL), 1.E-10)
      AV = V * 0.5
      BV = V + V
      NUMFUN = 30      ! Limit; FMINRC typically takes about 6 iterations
      LUNERR = ABS (LUNOUT)

   10 CONTINUE

         CALL FMINRC (AV, BV, V, OBJ, TOL, NUMFUN, CALLER, LUNOUT,
     >                ISTAT)

         IF (ISTAT < -1) THEN

            WRITE (LUNERR, '(/, 2A)') CALLER, ': FMINRC fatal error'
            STOP

         ELSE IF (ISTAT < 0) THEN ! Iteration limit; may be usable

            WRITE (LUNERR, '(/, 2A)') CALLER, ': Iteration limit.'

         ELSE IF (ISTAT > 0) THEN ! Evaluate the objective function

            CALL OBJECTIVE
            GO TO 10

         ELSE ! ISTAT = 0 (success).

         END IF

C     Ensure that everything matches V(best), not V(last):

      CALL OBJECTIVE

C     Calculate the optimized arc lengths:

      DO I = I1 + 1, I2
         ARC(I) = ARC(I - 1) + SQRT (
     >      (XI(I) - XI(I - 1)) ** 2 + (YI(I) - YI(I - 1)) ** 2 +
     >      (ZI(I) - ZI(I - 1)) ** 2)
      END DO

CCC   write (6, '(a, f12.7)') ' Current total length:', arc(i2)

      T = ONE / ARC(I2)

      ARC(I1:I2) = T * ARC(I1:I2) ! Renormalize

CCC   write (10, '(a)') ' Penultimate arc lengths:'
CCC   write (10, '(i3, 2f12.8)') (i, arc0(i), arc(i), i = i1, i2)

C     Recover the original relative spacing precisely:

      NEW = .TRUE.
      IEVAL = 1

      DO I = I1 + 1, I2 - 1
         CALL PLSCRV3D (NDATA, XI(I1), YI(I1), ZI(I1), ARC(I1),
     >                  METHOD, NEW, CLOSED, ARC0(I), IEVAL,
     >                  X(I), Y(I), Z(I), DERIVS)
         NEW = .FALSE.
      END DO

CCC   write (6, '(a, 2f12.7)') ' Target slopes: ',
CCC  >   s1(2)/s1(1), s2(2)/s2(1)
CCC   write (6, '(a, 2f12.7)') ' Final  slopes: ',
CCC  >   (y(i1+1) - y(i1)) / (x(i1+1) - x(i1)),
CCC  >   (y(i2) - y(i2-1)) / (x(i2) - x(i2-1))

C     Done.

C     Internal procedure for 2nd stage of 2-stage MORPH_LINE_3D algorithm:

      CONTAINS

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         SUBROUTINE OBJECTIVE

!        Evaluate a function of variable V being minimized to achieve specified
!        [end-point] slopes as well as possible.  V is a multiplier on the stage
!        1 total arc length that produces the final unknown arc length.
!        Actually, the optimal V turns out to differ, and it's not clear why.
!        Including a term aimed at preserving the original relative spacing as
!        well appears to be counter-productive.

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!        Local constants:

         REAL, PARAMETER ::
     >      HALF   = 0.5,
     >      WRANGE = 1.5 ! Highest : lowest weight = 1 : 1 / WRANGE
                         ! with parabolic variation from end to end; 10 too big!
!        Local variables:

         REAL
     >      A, B, DS, RESIDUAL1, RESIDUAL2, RESIDUAL3

!        Execution:

         IF (FIRST) THEN ! Factorize the matrix once only

            FIRST = .FALSE.

            B = ONE / WRANGE   ! Coeffs. of parabola on [0, 1] touching B at 0.5
            A = 4. * (ONE - B)

            DO I = I1, I2 - 1  ! Avoid starting W at I1+1; only the RH sides do
               W2 = (ARC(I) + ARC(I+1)) * HALF  ! Preserve any data symmetry
               W(I) = A * (W2 - HALF) ** 2 + B
            END DO

            CALL BIDECOM (NDATA, W, D, U, C, S) ! Upper bidiagonal + transforms

         END IF

C        The three coordinates present three right-hand-sides:

         T = V * TLENGTH1 ! Using V rather than T as the variable eases
                          ! bracketing the solution for any curve
         DO I = I1 + 1, I2
            W2 = (ARC(I-1) + ARC(I)) * HALF ! Normalized
            W1 = ONE - W2
            DS = (ARC(I) - ARC(I-1)) * T    ! Unnormalized
            B1(I) = DS * (W1 * S1(1) + W2 * S2(1))
            B2(I) = DS * (W1 * S1(2) + W2 * S2(2))
            B3(I) = DS * (W1 * S1(3) + W2 * S2(3))
         END DO

         I = I2
         B1(I) = B1(I) - XI(I)
         B2(I) = B2(I) - YI(I)
         B3(I) = B3(I) - ZI(I)

         I = I1 + 1
         B1(I) = B1(I) + XI(I1)
         B2(I) = B2(I) + YI(I1)
         B3(I) = B3(I) + ZI(I1)
 
         CALL BISOLVE (NDATA, W, D, U, C, S, B1(I), RESIDUAL1)
         CALL BISOLVE (NDATA, W, D, U, C, S, B2(I), RESIDUAL2)
         CALL BISOLVE (NDATA, W, D, U, C, S, B3(I), RESIDUAL3)

         OBJ = RESIDUAL1 + RESIDUAL2 + RESIDUAL3

         DO I = I1 + 1, I2 - 1
            XI(I) = B1(I)
            YI(I) = B2(I)
            ZI(I) = B3(I)
         END DO

CCC      write (6, '(a, 2f12.7)') ' Target   slopes: ',
CCC  >      s1(2)/s1(1), s2(2)/s2(1)
CCC      write (6, '(a, 2f12.7)') ' Iterated slopes: ',
CCC  >      (yi(i1+1) - yi(i1)) / (xi(i1+1) - xi(i1)),
CCC  >      (yi(i2) - yi(i2-1)) / (xi(i2) - xi(i2-1))

CCC      write (6, '(a, f16.13, 1p, 4e19.11)')
CCC  >      ' V/R123/OBJ: ',
CCC  >      v, residual1, residual2, residual3, obj

CCC      write (11, '(a, f16.13, 1p, 4e19.11)')
CCC  >      ' V/R123/OBJ: ',
CCC  >      v, residual1, residual2, residual3, obj

         END SUBROUTINE OBJECTIVE

      END SUBROUTINE MORPH_LINE_3D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine bidecom (n, w, d, u, c, s)
!
!        BIDECOM and BISOLVE treat a specialized weighted linear least squares
!     problem arising from perturbing a grid line in 3-space while controlling
!     the slopes at its new end points (and everywhere in between).  For n grid
!     points, n - 1 equations express the desired 2-point gradients but there
!     are only n - 2 interior grid points to adjust once the end point locations
!     have been imposed.
!
!        If the slopes in the n - 1 interior intervals are weighted by weights
!     w(1), w(2), ..., w(n-1), the relevant (n-1) x (n-2) matrix is:
!
!                       |  w1                              |
!                       | -w2  w2                          |
!                       |     -w3  w3                      |
!                       !         -w4  w4                  |
!                       :              :    :              |
!                       :                   :    :         |
!                       |                        :  w(n-2) |
!                       |                          -w(n-1) |
!
!        Applying a sequence of Givens matrices (symmetric plane rotations)
!     converts this to an upper bidiagonal system defined by vectors d(*) and
!     u(*).  BISOLVE applies the sequence (saved here in vectors c(*) and s(*))
!     to a given right-hand-side then solves the upper triangular system (or
!     rather the upper (n-2) x (n-2) portion of it) to produce the desired
!     linear least squares solution.  The minimized residual (2-norm squared)
!     is the square of the last element of the transformed right-hand side.
!
!     02/25/04  DAS  Initial implementation for MORPH_LINE_3D.
!
!     Author:  David Saunders, ELORET/NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

!     Arguments:

      integer, intent (in) :: n   ! Underlying number of grid points defining
                                  ! n - 1 intervals and an (n-1) x (n-2) matrix.
      real, intent (in) :: w(n-1) ! Weights scaling the matrix as shown.

      real, intent (out), dimension (n-1) :: d, u, c, s  ! See outline above.

!     Local variables:

      integer i

      real    gi, ci, si, wi

!     Execution:

      d(1) = w(1)   ! Initialize the forward pass

      do i = 2, n - 1

!        Construct the orthogonal matrix that zeros out -w(i), update the
!        diagonal above it and alongside it, and assign the corresponding
!        triangle factor element above the diagonal:

         wi = w(i)
         gi = sqrt (d(i-1)**2 + wi**2);  ! "Gamma"
         ci = d(i-1) / gi;               c(i) = ci
         si =   -wi  / gi;               s(i) = si
         d(i-1) =  gi
         u(i-1) =  si * wi
         d(i)   = -ci * wi

      end do

      end subroutine bidecom

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine bisolve (n, w, d, u, c, s, b, rsq)
!
!        BISOLVE completes the solution of a specialized weighted linear least
!     squares problem for a given right-hand-side vector b.  See BIDECOM for
!     further details.
!
!     02/26/04  DAS  Initial implementation for MORPH_LINE_3D.
!
!     Author:  David Saunders, ELORET/NASA Ames Research Center, CA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

!     Arguments:

      integer, intent (in) :: n   ! Underlying number of grid points defining
                                  ! n - 1 intervals and an (n-1) x (n-2) matrix
                                  ! which has been triangularized by BIDECOM
      real, intent (in) :: w(n-1) ! Weights for scaling the RHS vector b

      real, intent (in), dimension (n-1) :: d, u, c, s  ! Factorization info.
                                                        ! from BIDECOM
      real, intent (inout) :: b(n-1) ! Input with unweighted RHS;
                                     ! output with least squares solution
      real, intent (out) :: rsq      ! Minimized sum of squares 

!     Local variables:

      integer i

      real    bim1, ci, si

!     Execution:

!     Apply the weights to the RHS:

      do i = 1, n - 1
         b(i) = w(i) * b(i)
      end do

!     Transform the RHS proper:

      do i = 2, n - 1
         bim1   = b(i-1)
         b(i-1) = c(i) * bim1 + s(i) * b(i)
         b(i)   = s(i) * bim1 - c(i) * b(i)
      end do

      rsq = b(n - 1)**2

!     Solve the bidiagonal upper triangular system:

      i = n - 2
      b(i) = b(i) / d(i)

      do i = n - 3, 1, -1
         b(i) = (b(i) - u(i) * b(i+1)) / d(i)
      end do

      end subroutine bisolve
C+------------------------------------------------------------------------------
C
      SUBROUTINE NULINE3D (I1, I2, X0, Y0, Z0, X, Y, Z)
C
C  ONE-LINER: Adjust interior points of a 3-space line given new end points
C
C  DESCRIPTION:
C
C        NULINE3D perturbs the interior points of a 3-space curve given
C     perturbed end points.  More precisely, any pair of points I1, I2
C     may be the controlling points; each point between these two is moved
C     according to an arc-length-based combination of the distances between
C     the original and new points corresponding to I1 and I2.  The new
C     coordinates of these two points should be input in the desired output
C     curve.  In-place adjustment is permitted.
C
C  ENVIRONMENT:  Fortran 90
C
C  HISTORY:
C
C     03/04/94  DAS  Implemented for James Reuther's grid perturbation scheme.
C     01/09/95   "   Renamed NULINE3D from WARP3D since WARP3D now applies
C                    to all points of a 3-space grid block; I1, I2, not 1, N.
C     03/16/98   "   Eliminated ARC argument - make it an automatic array.
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER I1, I2                  !  I  Points given as perturbed
      REAL    X0 (*), Y0 (*), Z0 (*)  !  I  Original curve coordinates
      REAL    X (*), Y (*), Z (*)     ! I/O Desired curve, with points
                                      !     I1, I2 input

C-------------------------------------------------------------------------------

C     Local constants:

      REAL, PARAMETER :: ONE = 1., ZERO = 0.

C     Local variables:

      INTEGER I
      REAL    ARC(I1:I2), DX1, DX2, DY1, DY2, DZ1, DZ2, SINV, W1, W2

C     Execution:

      ARC (I1) = ZERO
      DO I = I1 + 1, I2
         ARC (I) = ARC (I - 1) + SQRT (
     >      (X0 (I) - X0 (I - 1)) ** 2 + (Y0 (I) - Y0 (I - 1)) ** 2 +
     >      (Z0 (I) - Z0 (I - 1)) ** 2)
      END DO

      DX1 = X (I1) - X0 (I1)
      DX2 = X (I2) - X0 (I2)
      DY1 = Y (I1) - Y0 (I1)
      DY2 = Y (I2) - Y0 (I2)
      DZ1 = Z (I1) - Z0 (I1)
      DZ2 = Z (I2) - Z0 (I2)
      SINV = ONE / ARC (I2)

      DO I = I1 + 1, I2 - 1
         W2 = ARC (I) * SINV
         W1 = ONE - W2
         X (I) = X0 (I) + W1 * DX1 + W2 * DX2
         Y (I) = Y0 (I) + W1 * DY1 + W2 * DY2
         Z (I) = Z0 (I) + W1 * DZ1 + W2 * DZ2
      END DO

      END SUBROUTINE NULINE3D
C+------------------------------------------------------------------------------
C
      SUBROUTINE PARAM2D (IDIM, JDIM, I1, I2, J1, J2, X, Y, Z, U, V)
C
C ONE-LINER: PARAMeterization in 2-space of an XYZ surface
C            -----               -
C PURPOSE:
C
C        PARAM2D parameterizes the given (sub)mesh on a surface using
C     the usual chord-length approximation to arc-length between grid
C     points.  The arc length of each row and column is generally
C     normalized to 1, but this version can suppress normalization.
C     U(I1,J1) should be input as 0. or -999. respectively.
C     This version also handles degenerate lines by inserting uniform
C     u or v for the normalized case.
C
C        It is hoped that tying the dimensions of U and V to those of
C     X and Y (to keep the argument list modest) is not unwise, though
C     it could be if small submeshes are commonly involved.
C
C ENVIRONMENT:
C
C     FORTRAN 77 + IMPLICIT NONE, trailing ! comments, and 8-char. names
C
C HISTORY:
C
C     11/19/93  DAS  Initial implementation.
C     12/03/93   "   Normalized each row and column.
C     03/20/95   "   Suppress normalization kludge: U(I1,J1) = -999.
C     02/06/98   "   Handled degenerate lines by inserting uniform u or v.
C
C AUTHOR: David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER IDIM, JDIM       ! (I) Grid dimensions in the calling program.

      INTEGER I1, I2, J1, J2   ! (I) Define the submesh to be parameterized.

      REAL    X (IDIM, JDIM),  ! (I) The surface grid coordinates.
     >        Y (IDIM, JDIM),
     >        Z (IDIM, JDIM)

      REAL    U (IDIM, JDIM),  ! (O) The chord-length-based parameterization:
     >        V (IDIM, JDIM)   !     U (I, J) = SUM (K=2:I) dS (K) / S (J) where
                               !     dS (K) ** 2 = (X (K, J) - X (K-1, J) ** 2 +
                               !                   (Y (K, J) - Y (K-1, J) ** 2 +
                               !                   (Z (K, J) - Z (K-1, J) ** 2
                               !     and S (J) is the total length of row J;
                               !     similarly for V (I, J).
                               ! (I) KLUDGE: If U(I1,J1) = -999. on input,
                               !     the normalization is suppressed.
C-------------------------------------------------------------------------------

C     Local constants.

      REAL      EPS, ONE, FLAG, ZERO
      PARAMETER
     >  (EPS = 1.E-6, FLAG = -999., ONE = 1., ZERO = 0.)

C     Local variables.

      INTEGER   I, J
      REAL      RLENGTH
      LOGICAL   NORM

C     Execution.

      NORM = U (I1, J1) .NE. FLAG

      DO J = J1, J2

         U (I1, J) = ZERO

         DO I = I1 + 1, I2
            U (I, J) = U (I - 1, J)  +  SQRT (
     >         (X (I, J) - X (I - 1, J)) ** 2 +
     >         (Y (I, J) - Y (I - 1, J)) ** 2 +
     >         (Z (I, J) - Z (I - 1, J)) ** 2)
         END DO

         IF (NORM) THEN

            IF (U (I2, J) .GT. EPS) THEN
               RLENGTH = ONE / U (I2, J)
               DO I = I1 + 1, I2 - 1
                  U (I, J) = U (I, J) * RLENGTH
               END DO
            ELSE
               RLENGTH = ONE / REAL (I2 - I1)
               DO I = I1 + 1, I2 - 1
                  U (I, J) = REAL (I - I1) * RLENGTH
               END DO
            END IF
            U (I2, J) = ONE

         END IF

      END DO


      DO I = I1, I2

         V (I, J1) = ZERO

         DO J = J1 + 1, J2
            V (I, J) = V (I, J - 1)  +  SQRT (
     >         (X (I, J) - X (I, J - 1)) ** 2 +
     >         (Y (I, J) - Y (I, J - 1)) ** 2 +
     >         (Z (I, J) - Z (I, J - 1)) ** 2)
         END DO

         IF (NORM) THEN

            IF (V (I, J2) .GT. EPS) THEN
               RLENGTH = ONE / V (I, J2)
               DO J = J1 + 1, J2 - 1
                  V (I, J) = V (I, J) * RLENGTH
               END DO
            ELSE
               RLENGTH = ONE / REAL (J2 - J1)
               DO J = J1 + 1, J2 - 1
                  V (I, J) = REAL (J - J1) * RLENGTH
               END DO
            END IF
            V (I, J2) = ONE

         END IF

      END DO

      RETURN
      END
C+------------------------------------------------------------------------------
C
      SUBROUTINE PARAMXYZ (IMIN, IMAX, JMIN, JMAX, KMIN, KMAX,
     >                     I1, I2, J1, J2, K1, K2, X, Y, Z, S)
C
C  ONE-LINER: Relative arc-lengths for all lines of a 3-space subgrid
C
C  DESCRIPTION:
C
C        PARAMXYZ parameterizes the indicated portion of a regular
C     3-space grid by setting up the normalized arc-length increments
C     in all three index directions.  If a single plane is specified
C     (e.g. I1 = I2), the degenerate direction is assumed to be unused.
C     (All 1s are returned since no attempt is made to suppress them.)
C
C  ENVIRONMENT:  FORTRAN 77 with minor extensions
C
C     12/06/94  DAS/JJR  Separated set-up of original arc-lengths from
C                        subsequent grid perturbations; normalized them.
C                        Including all the face lines proved tedious!
C     12/16/94   "   "   Allowed for sub-volumes (or planes).
C     08/27/96    DAS    Degenerate directions are now handled properly.
C
C  AUTHOR:  David Saunders/James Reuther, NASA Ames, Mt. View, CA.
C
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER IMIN, IMAX, JMIN, JMAX, KMIN, KMAX     ! I Grid array dimensions.
      INTEGER I1, I2, J1, J2, K1, K2                 ! I Indicate active volume.
      REAL    X (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX),   ! I Grid coordinates.
     >        Y (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX),
     >        Z (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX)
      REAL    S (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX, 3) ! O Relative arc-lengths
                                                     !   for the I,J,K lines:
                                                     !   S(I1,J,K,1) = 0.,
                                                     !   S(I,J1,K,2) = 0.,
                                                     !   S(I,J,K1,3) = 0.,
                                                     !   S(I2,J,K,1) = 1., etc.
C-------------------------------------------------------------------------------

C     Local constants.

      REAL    EPS, ONE, ZERO
      PARAMETER (EPS = 1.E-30, ONE = 1.E+0, ZERO = 0.E+0)

C     Local variables.

      INTEGER I, J, K
      REAL    RTOTAL

C     Local functions.

      REAL    DELI, DELJ, DELK

      DELI (I, J, K) = SQRT ((X (I, J, K) - X (I - 1, J, K)) ** 2 +
     >                       (Y (I, J, K) - Y (I - 1, J, K)) ** 2 +
     >                       (Z (I, J, K) - Z (I - 1, J, K)) ** 2)

      DELJ (I, J, K) = SQRT ((X (I, J, K) - X (I, J - 1, K)) ** 2 +
     >                       (Y (I, J, K) - Y (I, J - 1, K)) ** 2 +
     >                       (Z (I, J, K) - Z (I, J - 1, K)) ** 2)

      DELK (I, J, K) = SQRT ((X (I, J, K) - X (I, J, K - 1)) ** 2 +
     >                       (Y (I, J, K) - Y (I, J, K - 1)) ** 2 +
     >                       (Z (I, J, K) - Z (I, J, K - 1)) ** 2)

C     Execution.
C     ----------

C     Zero the three low-end faces (or edges if one plane is specified).

      DO K = K1, K2
         DO J = J1, J2
            S (I1, J, K, 1) = ZERO
         END DO

         DO I = I1, I2
            S (I, J1, K, 2) = ZERO
         END DO
      END DO

      DO J = J1, J2
         DO I = I1, I2
            S (I, J, K1, 3) = ZERO
         END DO
      END DO

C     Set up the low-end edge lines because they are missed by the
C     following loops over most of the low-end faces:

      DO I = I1 + 1, I2
         S (I, J1, K1, 1) =  S (I - 1, J1, K1, 1) + DELI (I, J1, K1)
      END DO

      DO J = J1 + 1, J2
         S (I1, J, K1, 2) =  S (I1, J - 1, K1, 2) + DELJ (I1, J, K1)
      END DO

      DO K = K1 + 1, K2
         S (I1, J1, K, 3) =  S (I1, J1, K - 1, 3) + DELK (I1, J1, K)
      END DO

C     Set up the rest of the low-end face lines because they are
C     missed by the the main loop over most of the volume.

      DO K = K1 + 1, K2
         DO J = J1 + 1, J2
            S (I1, J, K, 2) =  S (I1, J - 1, K, 2) + DELJ (I1, J, K)
            S (I1, J, K, 3) =  S (I1, J, K - 1, 3) + DELK (I1, J, K)
         END DO

         DO I = I1 + 1, I2
            S (I, J1, K, 1) =  S (I - 1, J1, K, 1) + DELI (I, J1, K)
            S (I, J1, K, 3) =  S (I, J1, K - 1, 3) + DELK (I, J1, K)
         END DO
      END DO

      DO J = J1 + 1, J2
         DO I = I1 + 1, I2
            S (I, J, K1, 1) =  S (I - 1, J, K1, 1) + DELI (I, J, K1)
            S (I, J, K1, 2) =  S (I, J - 1, K1, 2) + DELJ (I, J, K1)
         END DO
      END DO

C     Traverse the block just once for all lines except those within
C     the low-end faces.

      DO K = K1 + 1, K2
         DO J = J1 + 1, J2
            DO I = I1 + 1, I2
               S (I, J, K, 1) =  S (I - 1, J, K, 1) + DELI (I, J, K)
               S (I, J, K, 2) =  S (I, J - 1, K, 2) + DELJ (I, J, K)
               S (I, J, K, 3) =  S (I, J, K - 1, 3) + DELK (I, J, K)
            END DO
         END DO
      END DO

C     Normalizing requires another pass through the volume.
C     Degenerate cases must be avoided, so separate the 3 directions:

      DO K = K1, K2
         DO J = J1, J2
            RTOTAL = ONE / MAX (S (I2, J, K, 1), EPS)
            DO I = I1 + 1, I2 - 1  ! No trips if I1 = I2
               S (I, J, K, 1) = S (I, J, K, 1) * RTOTAL
            END DO
            S (I2, J, K, 1) = ONE
         END DO
      END DO

      DO K = K1, K2
         DO I = I1, I2
            RTOTAL = ONE / MAX (S (I, J2, K, 2), EPS)
            DO J = J1 + 1, J2 - 1
               S (I, J, K, 2) = S (I, J, K, 2) * RTOTAL
            END DO
            S (I, J2, K, 2) = ONE
         END DO
      END DO

      DO J = J1, J2
         DO I = I1, I2
            RTOTAL = ONE / MAX (S (I, J, K2, 3), EPS)
            DO K = K1 + 1, K2 - 1
               S (I, J, K, 3) = S (I, J, K, 3) * RTOTAL
            END DO
            S (I, J, K2, 3) = ONE
         END DO
      END DO

      RETURN
      END
C+------------------------------------------------------------------------------
C
      SUBROUTINE PBILINT (MODE, IDIM, JDIM, I1, I2, J1, J2, XS, YS, ZS,
     >                    US, VS, UTARG, VTARG, IS, JS, EPS, P, Q,
     >                    XINT, YINT, ZINT, XYZU, XYZV, IER)
C
C     Parametric bilinear interpolation within a regular 3-space surface mesh.
C     Arguments are much as for PLBICUBE, q.v.  New arguments P, Q are outputs.
C
C     Jan. 95   DAS  Initial implementation (no derivatives).
C     08/06/99   "   Analytic derivatives can come from the BILINT formulation.
C                    Also, return the relevant (p,q) that is sometimes needed.
C
C     Author:  David Saunders  Sterling Software/NASA Ames, Mt. View, CA.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   MODE, IDIM, JDIM, I1, I2, J1, J2

      REAL, INTENT (IN), DIMENSION (IDIM, JDIM) ::
     >   XS, YS, ZS, US, VS 

      REAL, INTENT (IN) ::
     >   UTARG, VTARG, EPS

      INTEGER, INTENT (INOUT) ::
     >   IS, JS

      REAL, INTENT (OUT) ::
     >   P, Q, XINT, YINT, ZINT, XYZU (3), XYZV (3)

      INTEGER, INTENT (OUT) ::
     >   IER

C     Local constants:

      REAL, PARAMETER ::
     >   ONE = 1.

C     Local variables:

      REAL
     >   A (2, 2), AJ (2, 2), DP (2), FP, FQ, FPQ, PM1, QM1

C     Procedures:

      EXTERNAL
     >   RIPPLE2D

C     Execution:
      
C     Locate the enclosing cell:

      CALL RIPPLE2D (IDIM, JDIM, I1, I2, J1, J2, US, VS,
     >               UTARG, VTARG, IS, JS, EPS, P, Q, IER)

C     Ignore IER - the nearest cell is always returned.

      PM1  = ONE - P
      QM1  = ONE - Q

      XINT = QM1 * (PM1 * XS (IS, JS)   + P * XS (IS+1, JS)) +
     >         Q * (PM1 * XS (IS, JS+1) + P * XS (IS+1, JS+1))
      YINT = QM1 * (PM1 * YS (IS, JS)   + P * YS (IS+1, JS)) +
     >         Q * (PM1 * YS (IS, JS+1) + P * YS (IS+1, JS+1))
      ZINT = QM1 * (PM1 * ZS (IS, JS)   + P * ZS (IS+1, JS)) +
     >         Q * (PM1 * ZS (IS, JS+1) + P * ZS (IS+1, JS+1))

C     Interpolated derivatives as well?

      IF (MODE /= 0) THEN

C        The BILINT formulation provides analytic partial dx/dp, etc.
C        The chain rule converts p/q derivatives to u/v derivatives.

C        Set up the 2 x 2 LHS matrix common to derivatives of x, y, and z.
C        The algebra saves a few operations over the more obvious form.

         FP  = US (IS + 1, JS) - US (IS, JS)
         FQ  = US (IS, JS + 1) - US (IS, JS)
         FPQ = US (IS + 1, JS + 1) - US (IS, JS + 1) - FP

         AJ (1, 1) = FP + Q * FPQ ! du/dp
         AJ (2, 1) = FQ + P * FPQ ! du/dq

         FP  = VS (IS + 1, JS) - VS (IS, JS)
         FQ  = VS (IS, JS + 1) - VS (IS, JS)
         FPQ = VS (IS + 1, JS + 1) - VS (IS, JS + 1) - FP

         AJ (1, 2) = FP + Q * FPQ ! dv/dp
         AJ (2, 2) = FQ + P * FPQ ! dv/dq

C        For X ...

         FP  = XS (IS + 1, JS) - XS (IS, JS)
         FQ  = XS (IS, JS + 1) - XS (IS, JS)
         FPQ = XS (IS + 1, JS + 1) - XS (IS, JS + 1) - FP

         DP (1) = FP + Q * FPQ ! dx/dp
         DP (2) = FQ + P * FPQ ! dx/dq

         A = AJ ! Could factorize once & solve 3 times, but LUSOLVE
                ! is in use elsewhere, while DECOMP & SOLVE are not.

         CALL LUSOLVE (2, 2, A, DP, IER) ! Solve J dfuv = dfpq

         XYZU (1) = DP (1) ! dx/du
         XYZV (1) = DP (2) ! dx/dv

C        For Y ...

         FP  = YS (IS + 1, JS) - YS (IS, JS)
         FQ  = YS (IS, JS + 1) - YS (IS, JS)
         FPQ = YS (IS + 1, JS + 1) - YS (IS, JS + 1) - FP

         DP (1) = FP + Q * FPQ ! dy/dp
         DP (2) = FQ + P * FPQ ! dy/dq

         A = AJ

         CALL LUSOLVE (2, 2, A, DP, IER)

         XYZU (2) = DP (1) ! dy/du
         XYZV (2) = DP (2) ! dy/dv

C        For Z ...

         FP  = ZS (IS + 1, JS) - ZS (IS, JS)
         FQ  = ZS (IS, JS + 1) - ZS (IS, JS)
         FPQ = ZS (IS + 1, JS + 1) - ZS (IS, JS + 1) - FP

         DP (1) = FP + Q * FPQ ! dz/dp
         DP (2) = FQ + P * FPQ ! dz/dq

         A = AJ

         CALL LUSOLVE (2, 2, A, DP, IER)

         XYZU (3) = DP (1) ! dz/du
         XYZV (3) = DP (2) ! dz/dv

      END IF

      END SUBROUTINE PBILINT
C+------------------------------------------------------------------------------
C
      SUBROUTINE PLSCRV3D (NDATA, X, Y, Z, T, METHOD, NEW, CLOSED,
     >                     TEVAL, IEVAL, XEVAL, YEVAL, ZEVAL, DERIVS)
C
C     Description:
C
C        PLSCRV3D interpolates X, Y, Z along a curve at a specified TEVAL using
C     piecewise linear or local cubic spline techniques.  It is an adaptation of
C     the earlier PLSFIT3D and PLSCURVE intended to avoid internal arc length
C     calculations by assuming that the cumulative chord lengths associated with
C     the data points are supplied by the application. (PLSFIT3D avoided storing
C     all the Ts as part of its space-efficient graphics-oriented approach, but
C     in practical grid generation applications the Ts are needed anyway.)
C
C        PLSCRV3D restores the METHOD options that were dropped by PLSCURVE,
C     in anticipation of a requirement for piecewise linear results.  The option
C     to return derivatives for intersection calculations influenced the choice
C     between one and many target Ts: a single T wins.  Use of normalized Ts is
C     recommended for intersection calculations, so T is known to be in [0, 1],
C     but other applications may prefer not to normalize, and this is permitted.
C     Arithmetic for unwanted derivatives is avoided by calls with input
C     DERIVS(1) = -999.
C
C        See Robert Kennelly's original 2-space PLSFIT for more on compact local
C     spline methods.
C
C     Environment:  Fortran 90
C
C     History:
C
C     09/22/97  DAS  Initial adaptation of PLSFIT3D/PLSCURVE with more efficient
C                    curve/surface intersections and wing surface grids in mind,
C                    especially NDATA = 2 and piecewise linear cases.
C     11/17/97   "   Allowed for non-normalized input Ts.
C     05/06/98   "   Minor Fortran 90 updates.
C
C     Author:  David Saunders, Sterling Software/NASA Ames, Mtn. View, CA
C
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   NDATA               ! Length of X, Y, Z input data arrays.

      REAL, INTENT (IN) ::
     >   X(NDATA),           ! Input curve coordinates. Successive
     >   Y(NDATA),           ! data points must be distinct.
     >   Z(NDATA)            ! If CLOSED, then first and last points
                             ! must agree (not checked here).
      REAL, INTENT (IN) ::
     >   T(NDATA)            ! Parametric variable values associated
                             ! with the data points, usually normalized
                             ! cumulative chord lengths, but they need
                             ! not be normalized except for intersection
                             ! calculations via INTSEC4 or -5.
      CHARACTER, INTENT (IN) ::
     >   METHOD*1            ! The type of parametric fit to be used:
                             !    'M' means monotonic piecewise cubics;
                             !    'B' means non-monotonic "Bessel"-type
                             !        piecewise cubics (looser fit);
                             !    'L' means piecewise linear (only way
                             !        to guarantee no overshoots).
      LOGICAL, INTENT (IN) ::
     >   NEW                 ! If NEW = T, the search for a bracket starts
                             ! from scratch, otherwise locally-saved search
                             ! and fit information will be assumed to be
                             ! correct. If calling PLSCRV3D from within a
                             ! loop, set NEW = F after the first call.
      LOGICAL, INTENT (IN) ::
     >   CLOSED              ! CLOSED = T means periodic boundary conditions
                             ! are to be used. The curve should wrap around
                             ! smoothly on itself. The calling routine must
                             ! ensure that the end points agree.
      REAL, INTENT (IN) ::
     >   TEVAL               ! Target value of T at which to evaluate X,Y,Z.

      INTEGER, INTENT (INOUT) ::
     >   IEVAL               ! Input with the data index at which to start
                             ! searching for TEVAL; output with actual index.

      REAL, INTENT (OUT) ::
     >   XEVAL, YEVAL, ZEVAL ! Interpolated coordinates.

      REAL, INTENT (INOUT) ::
     >   DERIVS(3)           ! Partial derivatives of X, Y, Z with respect to T.
                             ! Enter DERIVS(1) = -999. to suppress these,
                             ! else  DERIVS(1) = 0. (say) on the first call
                             ! and something other than -999. thereafter.
C     Procedures:

      REAL, EXTERNAL ::
     >   BESSEL,             ! 1st deriv. (central 3-point formula)
     >   BRODLIE,            ! 1st deriv. (central) adjusted for monotonicity
     >   BUTLAND,            ! 1st deriv. (non-central)  "    "    "
     >   THREEPT             ! 1st deriv. (non-central 3-point formula)
      EXTERNAL
     >   INTERVAL            ! Interpolatory 1-D search utility
C-------------------------------------------------------------------------------

C     Local constants:

      REAL, PARAMETER :: ONE = 1., TWO = 2., THREE = 3.

C     Local variables:

      LOGICAL    DERIV, LINEAR, LOOSE, MEMORY, TWOPTS
      INTEGER    IND (-1:2), J, LEFT, RIGHT
      REAL       BX (0:1), BY (0:1), BZ (0:1), CX, CY, CZ, DT,
     >           DELX (-1:1), DELY (-1:1), DELZ (-1:1), DX, DY, DZ,
     >           H (-1:1), RH, TLEFT, TRIGHT

C     Storage:

      SAVE       BX, BY, BZ, CX, CY, CZ, DX, DY, DZ, LEFT, RIGHT,
     >           TLEFT, TRIGHT, DERIV, LINEAR, LOOSE, TWOPTS

C     Execution:

C     The 2-data-point case is common enough for intersections that it
C     is treated specially.  Minor code duplication is OK.

      IF (NEW) THEN

         DERIV  = DERIVS (1) /= -999.
         TWOPTS = NDATA == 2

         IF (TWOPTS) THEN

            RH     = ONE / (T (2) - T (1))
            BX (0) = (X (2) - X (1)) * RH
            BY (0) = (Y (2) - Y (1)) * RH
            BZ (0) = (Z (2) - Z (1)) * RH

            DT    = TEVAL - T (1)
            XEVAL = X (1) + DT * BX (0)
            YEVAL = Y (1) + DT * BY (0)
            ZEVAL = Z (1) + DT * BZ (0)

            IF (DERIV) THEN
               DERIVS (1) = BX (0)
               DERIVS (2) = BY (0)
               DERIVS (3) = BZ (0)
            END IF

            GO TO 999

         ELSE
            MEMORY = .FALSE.
            LINEAR = METHOD == 'L'
            LOOSE  = METHOD == 'B' ! Else monotonic
         END IF

      ELSE IF (TWOPTS) THEN

         DT    = TEVAL - T (1)
         XEVAL = X (1) + DT * BX (0)
         YEVAL = Y (1) + DT * BY (0)
         ZEVAL = Z (1) + DT * BZ (0)

         IF (DERIV) THEN  ! INTSEC5 overwrites previous derivatives
            DERIVS (1) = BX (0)
            DERIVS (2) = BY (0)
            DERIVS (3) = BZ (0)
         END IF

         GO TO 999

      ELSE ! NEW = F and NDATA > 2

C        We can save time when PLSCRV3D is being called from within a loop
C        by setting MEMORY if possible.  The out-of-range checking relies on
C        the fact that RIGHT = LEFT + 1 after the last return.  Cater to the
C        more likely case of TEVAL in the previous (same) interior interval.

         MEMORY = TEVAL >= TLEFT .AND. TEVAL < TRIGHT

         IF (.NOT. MEMORY) THEN
            MEMORY = (LEFT  == 1)     .AND. (TEVAL <  TRIGHT) .OR.
     >               (RIGHT == NDATA) .AND. (TEVAL >= TLEFT)
         END IF

      END IF

      IF (.NOT. MEMORY) THEN

C        Interpolation search for bracketing interval:

         LEFT = IEVAL

         CALL INTERVAL (NDATA, T, TEVAL, ONE, LEFT)

         IEVAL  = LEFT
         TLEFT  = T (LEFT)
         RIGHT  = LEFT + 1
         TRIGHT = T (RIGHT)

C         -------------------------------------------
C        |   1 <= LEFT < RIGHT = LEFT + 1 <= NDATA   |
C        |   TLEFT <= TEVAL <= TRIGHT                |
C         -------------------------------------------

         IF (LINEAR) THEN ! Piecewise linear, NDATA > 2

            RH     = ONE / (TRIGHT - TLEFT)
            BX (0) = (X (RIGHT) - X (LEFT)) * RH
            BY (0) = (Y (RIGHT) - Y (LEFT)) * RH
            BZ (0) = (Z (RIGHT) - Z (LEFT)) * RH

         ELSE ! Piecewise cubic

C           Three cases are handled together: (1) both endpoints are
C           interior to the data range, or the interval is at the
C           beginning/end of the range with either (2) periodic,
C           or (3) free end conditions.  For special cases (2) and (3),
C           the initial calculations will be overridden below.

            IND (-1) = LEFT - 1
            IND ( 0) = LEFT
            IND (+1) = RIGHT
            IND (+2) = RIGHT + 1

C           Patch index array for periodic case (2). (Later, the free end
C           case will be overridden again, following derivative calculation.)

            IF (LEFT == 1) THEN
               IND (-1) = NDATA - 1 ! Left side wrap-around boundary condition
            ELSE IF (RIGHT == NDATA) THEN
               IND (2) = 2          ! Right side
            END IF

C           Interval and derivative approximations:

            H (-1) = T (IND (-1) + 1) - T (IND (-1))
            H ( 0) = T (RIGHT) - T (LEFT)
            H (+1) = T (IND (2)) - T (IND (2) - 1)

            DO J = -1, +1
               RH       = ONE / H (J)
               DELX (J) = (X (IND (J + 1)) - X (IND (J))) * RH
               DELY (J) = (Y (IND (J + 1)) - Y (IND (J))) * RH
               DELZ (J) = (Z (IND (J + 1)) - Z (IND (J))) * RH
            END DO

C           Select the interpolation scheme, and compute [adjusted] first
C           derivatives at both left- and right-hand endpoints of the interval:

            IF (LOOSE) THEN ! Bessel - use the central formula at 0, +1

               BX (0) = BESSEL (0, H, DELX)
               BX (1) = BESSEL (1, H, DELX)
               BY (0) = BESSEL (0, H, DELY)
               BY (1) = BESSEL (1, H, DELY)
               BZ (0) = BESSEL (0, H, DELZ)
               BZ (1) = BESSEL (1, H, DELZ)

            ELSE ! Monotone - use the Brodlie modification of Butland's formula

               BX (0) = BRODLIE (0, H, DELX)
               BX (1) = BRODLIE (1, H, DELX)
               BY (0) = BRODLIE (0, H, DELY)
               BY (1) = BRODLIE (1, H, DELY)
               BZ (0) = BRODLIE (0, H, DELZ)
               BZ (1) = BRODLIE (1, H, DELZ)

            END IF

C           Patch initial/final derivatives if not periodic, case (3):

            IF (.NOT. CLOSED) THEN
               IF (LEFT == 1) THEN
                  IF (LOOSE) THEN
                     BX (0) = THREEPT (0, H, DELX)
                     BY (0) = THREEPT (0, H, DELY)
                     BZ (0) = THREEPT (0, H, DELZ)
                  ELSE
                     BX (0) = BUTLAND (0, H, DELX)
                     BY (0) = BUTLAND (0, H, DELY)
                     BZ (0) = BUTLAND (0, H, DELZ)
                  END IF
               ELSE IF (RIGHT == NDATA) THEN
                  IF (LOOSE) THEN
                     BX (1) = THREEPT (1, H, DELX)
                     BY (1) = THREEPT (1, H, DELY)
                     BZ (1) = THREEPT (1, H, DELZ)
                  ELSE
                     BX (1) = BUTLAND (1, H, DELX)
                     BY (1) = BUTLAND (1, H, DELY)
                     BZ (1) = BUTLAND (1, H, DELZ)
                  END IF
               END IF
            END IF

C           Compute the remaining cubic coefficients relative
C           to the left-hand endpoint.

            RH = ONE / H (0)
            CX = (THREE * DELX (0) - TWO * BX (0) - BX (1)) * RH
            CY = (THREE * DELY (0) - TWO * BY (0) - BY (1)) * RH
            CZ = (THREE * DELZ (0) - TWO * BZ (0) - BZ (1)) * RH
            DX = (-TWO * DELX (0) + BX (0) + BX (1)) * RH ** 2
            DY = (-TWO * DELY (0) + BY (0) + BY (1)) * RH ** 2
            DZ = (-TWO * DELZ (0) + BZ (0) + BZ (1)) * RH ** 2

         END IF

      END IF

C     Evaluate the identified polynomials:

      DT = TEVAL - TLEFT

      IF (LINEAR) THEN

         XEVAL = X (LEFT) + DT * BX (0)
         YEVAL = Y (LEFT) + DT * BY (0)
         ZEVAL = Z (LEFT) + DT * BZ (0)

         IF (DERIV) THEN
            DERIVS (1) = BX (0)
            DERIVS (2) = BY (0)
            DERIVS (3) = BZ (0)
         END IF

      ELSE ! Parametric cubics

         XEVAL = X (LEFT) + DT * (BX (0) + DT * (CX + DT * DX))
         YEVAL = Y (LEFT) + DT * (BY (0) + DT * (CY + DT * DY))
         ZEVAL = Z (LEFT) + DT * (BZ (0) + DT * (CZ + DT * DZ))

         IF (DERIV) THEN
            DERIVS (1) = BX (0) + DT * (TWO * CX + DT * THREE * DX)
            DERIVS (2) = BY (0) + DT * (TWO * CY + DT * THREE * DY)
            DERIVS (3) = BZ (0) + DT * (TWO * CZ + DT * THREE * DZ)
         END IF

      END IF

  999 RETURN

      END SUBROUTINE PLSCRV3D
C+------------------------------------------------------------------------------
C
      SUBROUTINE PROJECT4 (X1, X2, X3, X4, XTARG, X0, DSQ, P, Q, IER)
C
C     One-liner: Orthogonal projection of a point onto a 4-sided cell in 3-space
C
C     Purpose:
C
C        PROJECT4 projects point XTARG onto the 4-sided cell defined by
C     X1, X2, X3, and X4, with a "counterclockwise" convention employed.
C     More precisely, since the quadrilateral may not be planar, the method
C     determines the (p,q) which minimizes the distance of XTARG to the
C     quadratic surface defined by
C
C        X0(p,q) = (1-p)(1-q) X1 + p(1-q) X2 + pq X3 + (1-p)q X4
C
C     where point X1 is taken as the origin for p and q.  Also returned are
C     the projected point, X0, and the (squared) distance to X0.
C
C        This version solves the nonlinear least squares problem via a
C     Newton iteration.  Starting guesses of 0.5 are used for p and q.
C
C     Usage:
C
C        Treating the quadrilateral as two triangles may be preferable (using
C     PROJECT3 twice) if the nearest cell to the target point is being sought,
C     because projection to a triangle leads to a noniterative linear method.
C     PROJECT4 allows the (p,q) from PROJECT3 to be converted to equivalent
C     fractional coordinates referred to the associated quadrilateral.  If
C     the quad. is a planar parallelogram, the two (p,q)s are the same, but
C     in general the quad. may be nonplanar, yet a best solution is provided
C     for that cell which can be adjusted for bidirectional parametric
C     interpolation within a structured surface grid.
C
C        PROJECT4 can be thought of as solving the inverse of the problem
C     solved by BILINT:  for a given target (u,v), BILINT calculates (p,q)
C     within the specified cell in 2-space (u/v space), from which (x,y,z)
C     can be interpolated, whereas here we need to determine the (p,q) or
C     (u,v) corresponding to a given (x,y,z) point near a cell in 3-space.
C     
C     History:
C
C     Nov. 1998  J. Reuther  Initial implementation of Newton method.
C      "    "    D. Saunders Gauss-Newton iteration analogous to the Newton
C                            iteration of BILINT seems OK for ~planar cells.
C     Sep. 1999    "     "   Clean-up and revised arguments after belatedly
C                            comparing the two approaches and finding that
C                            they match for planar quads., but full Newton
C                            is much more robust for nonplanar cells.
C                            Added step-halving to cover extreme cases,
C                            including reversal of uphill directions.
C     10/22/99     "     "   Spurious IER = 3 resulted if initial DP < EPS.
C                            Therefore, distinguish this from the step-halving
C                            convergence test, which also checks for DP < EPS.
C     07/25/03     "     "   Instead of reversing search direction, ensure that
C                            the Hessian is positive definite to begin with.
C                            Three-way convergence test is more thorough too.
C
C     AUTHOR: James Reuther/David Saunders, NASA Ames Research Center, CA.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      REAL, INTENT (IN), DIMENSION (3) ::
     >   X1, X2, X3, X4, ! (x,y,z) coordinates of the quadrilateral vertices
     >   XTARG           ! Target (x,y,z) point being projected

      REAL, INTENT (OUT) ::
     >   X0(3),          ! Coordinates of the projected point
     >   DSQ,            ! Corresponding squared distance from XTARG
     >   P, Q            ! Coefficients usable for bidirectional interpolation
                         ! as described above, referred to vertex X1

      INTEGER, INTENT (OUT) ::
     >   IER             ! 0: The iteration converged with p & q in [0,1];
                         ! 1: The iteration converged; p & q not in [0,1];
                         ! 2: The iteration limit was reached;
                         ! 3: The step halving failed - no convergence

C     Local constants:

      INTEGER, PARAMETER ::
     >   MAXITER = 40    ! Normally much less; positive definiteness needs more

      REAL, PARAMETER ::
     >   BIG = 1.E+30, HALF = 0.5, ONE = 1., ZERO = 0.

C     Local variables:

      INTEGER
     >   N

      REAL
     >   A(2,2), DV(2), G(2), ALPHA, DELTA, DET_G, DET_MIN,
     >   DV_NORM, G_NORM, V_NORM, F0, P0, Q0, ROOT_EPS,
     >   TOL_F, TOL_F_HALF, TOL_F_3RD,
     >   AX1, AX2, AX3, AX4, DX, DX_DP, DX_DQ,
     >   AY1, AY2, AY3, AY4, DY, DY_DP, DY_DQ,
     >   AZ1, AZ2, AZ3, AZ4, DZ, DZ_DP, DZ_DQ

      LOGICAL
     >   CONVGD, REVERSED

C     Execution:

      AX1 = X1(1)
      AX2 = X2(1) - AX1
      AX3 = X4(1) - AX1
      AX4 = X3(1) - X2(1) - AX3

      AY1 = X1(2)
      AY2 = X2(2) - AY1
      AY3 = X4(2) - AY1
      AY4 = X3(2) - X2(2) - AY3

      AZ1 = X1(3)
      AZ2 = X2(3) - AZ1
      AZ3 = X4(3) - AZ1
      AZ4 = X3(3) - X2(3) - AZ3

      IF (PRECISION (TOL_F) > 10) THEN ! # decimal digits
         ROOT_EPS   = 1.E-7    ! True sqrt (eps) allows too many halvings
         TOL_F      = 1.E-12   ! Seek 12 digits in function DSQ
         TOL_F_HALF = 1.E-6
         TOL_F_3RD  = 1.E-4
      ELSE
         ROOT_EPS   = 1.E-4
         TOL_F      = 1.E-6
         TOL_F_HALF = 1.E-3
         TOL_F_3RD  = 1.E-2
      END IF

      DET_MIN = TOL_F ! Unclear what is best for det (Hessian)
      DV_NORM = BIG   ! Initialize measure of convergence and ...
      F0      = BIG   ! ... measure of reduced function value
      P       = HALF
      Q       = HALF
      IER     = 0
      ALPHA   = ZERO

C     Minimize the projected distance to the quadrilateral:

      DO N = 1, MAXITER
 
  100    CONTINUE ! Step-halving subiteration starts here

         X0(1)  = AX1 + AX2 * P + (AX3 + AX4 * P) * Q 
         X0(2)  = AY1 + AY2 * P + (AY3 + AY4 * P) * Q
         X0(3)  = AZ1 + AZ2 * P + (AZ3 + AZ4 * P) * Q

         DX     = X0(1) - XTARG(1)
         DY     = X0(2) - XTARG(2)
         DZ     = X0(3) - XTARG(3)

         DSQ    = DX ** 2 + DY ** 2 + DZ ** 2

         DX_DP  = AX2 + AX4 * Q
         DY_DP  = AY2 + AY4 * Q
         DZ_DP  = AZ2 + AZ4 * Q

         DX_DQ  = AX3 + AX4 * P
         DY_DQ  = AY3 + AY4 * P
         DZ_DQ  = AZ3 + AZ4 * P

C        Gradient vector:

         G(1)   = DX * DX_DP + DY * DY_DP + DZ * DZ_DP
         G(2)   = DX * DX_DQ + DY * DY_DQ + DZ * DZ_DQ

         G_NORM = MAX (ABS (G(1)), ABS (G(2)))
         V_NORM = MAX (ABS (P),    ABS (Q))

CCCCC    write (6, '(i3, a, 1p, e14.7, a, 0p, f9.6, a, 1p, 3e22.14)')
CCCCC>      n, '  ||g||:', g_norm, '  Step:', alpha,
CCCCC>      '  dsq, p, q:', dsq, p, q

C        Convergence tests from Practical Optimization, 8.2.3.2:

         IF (F0 - DSQ < TOL_F      * (ONE + DSQ)     .AND.
     >       G_NORM   < TOL_F_3RD  * (ONE + DSQ)     .AND.
     >       DV_NORM  < TOL_F_HALF * (ONE + V_NORM)) THEN
            EXIT
         END IF

C        Halve the step if the projected distance has not reduced:

         IF (DSQ > F0) THEN ! Not when N = 1

            IF (DV_NORM > ROOT_EPS * (ONE + V_NORM)) THEN

               ALPHA =  ALPHA * HALF
               P = P0 - ALPHA * DV(1)
               Q = Q0 - ALPHA * DV(2)
               DV_NORM = HALF * DV_NORM
               GO TO 100

            END IF

            IER = 3
            EXIT ! Step-halving failed - give up

         END IF

C        Calculate a search direction:

         A(1,1) = DX_DP ** 2 + DY_DP ** 2 + DZ_DP ** 2
         A(2,2) = DX_DQ ** 2 + DY_DQ ** 2 + DZ_DQ ** 2
         A(1,2) = DX * AX4 + DX_DP * DX_DQ +
     >            DY * AY4 + DY_DP * DY_DQ +
     >            DZ * AZ4 + DZ_DP * DZ_DQ
         A(2,1) = A(1,2) 

C        Ensure positive definiteness for a descent direction:

         DET_G  = A(1,1) * A(2,2) - A(1,2) ** 2

         IF (DET_G < DET_MIN) THEN
            DELTA = (DET_MIN - DET_G) / (A(1,1) + A(2,2))
CCCCC       write (6, '(A, 1P, 5E20.10)')
CCCCC>         '  DET_G, DELTA, A11, A22, A21: ',
CCCCC>            DET_G, DELTA, A(1,1), A(2,2), A(1,2)
            A(1,1) = A(1,1) + DELTA
            A(2,2) = A(2,2) + DELTA
         END IF

         CALL LUSOLVE (2, 2, A, G, IER) ! Solve G dv = g

         DV      = G
         DV_NORM = MAX (ABS (DV(1)), ABS (DV(2)))
         P0      = P
         Q0      = Q
         P       = P0 - DV(1)
         Q       = Q0 - DV(2)
         F0      = DSQ
         ALPHA   = ONE

      END DO


      IF (N >= MAXITER) THEN
         IER = 2            ! Too many iterations
      ELSE IF (IER /= 3) THEN
         IF (MAX (ABS (P - HALF), ABS (Q - HALF)) > HALF + ROOT_EPS)
     >      IER = 1         ! Converged to a point outside the cell
      END IF                ! Else IER = 0 from LUSOLVE or initialization

      END SUBROUTINE PROJECT4
C+------------------------------------------------------------------------------
C
      SUBROUTINE RIPPLE2D (MDIM, NDIM, M1, M2, N1, N2, X, Y, XP, YP,
     >                     M, N, EPS, P, Q, STATUS)
C
C     ONE-LINER:  "RIPPLE" search of a 2D quasi-rectangular grid
C                  ------              --
C     PURPOSE:
C
C        RIPPLE2D performs a search within a quasi-rectangular 2D grid for
C     the cell containing a target point by using a "ripple" pattern.  The
C     search should be efficient for typical applications where the specified
C     initial cell is a good starting guess (probably the result of a previous
C     call).  It also returns the corresponding bilinear interpolation
C     coefficients (though they may not apply to, say, bicubic applications).
C
C     METHOD:
C
C        Points are tested around the starting guess in rectangular "layers"
C     of increasing size.  If the cell defined by (M, N) does not contain
C     the target point, then the next search domain would be the surrounding
C     grid layer defined by the corner cells (M + 1, N + 1), (M + 1, N - 1),
C     (M - 1, N - 1), and (M - 1, N + 1).  If the cell is not found then the
C     next domain would be the layer defined by (M + 2, N + 2), (M + 2, N - 2),
C     (M - 2, N - 2), and (M - 2, N + 2), and so on.  The indices for each of
C     the corners are incremented or decremented until they reach the indicated
C     grid limits (M1 <= I <= M2 - 1 and N1 <= J <= N2 - 1) after which they
C     are held constant.  The search stops when a match is found or after all
C     corners of the layer have reached the corners of the (sub)grid, by which
C     time all specified cells will have been tested.
C
C        Within the rectangular layer, the corner points and the points interior
C     to them, here called edges, are treated separately for the sake of
C     symmetry.  An edge is no longer tested after the boundary of a grid has
C     been encountered by it.  Both corners and edges are cycled through in
C     clockwise directions as follows (where L is the current layer or level):
C
C
C          (M - L, N + L)                           (M + L, N + L)
C
C                          4 +. . . . . . . . .+ 1
C                            .        4        .  \
C                            .                 .   \
C                            .                 .    Corner
C                            .                 .
C                            . 3      +      1<------Edge
C                            .         (M,N)   .
C                            .                 .
C                            .                 .
C                            .        2        .
C                          3 +. . . . . . . . .+ 2
C
C          (M - L, N - L)                           (M + L, N - L)
C
C
C        The subroutine BILINT is used to perform the test.  BILINT returns
C     bilinear interpolation coefficients which are not used here but are
C     returned to the higher level.
C
C     ARGUMENTS:
C
C     NAME    DIM    TYPE  I/O/S     DESCRIPTION
C     ----    ---    ----  -----     -----------
C     MDIM,    -       I   I         Dimensions of the grid arrays in the
C     NDIM                           calling program
C     M1,      -       I   I         Indices specifying the (sub)grid to be
C     M2,                            searched as points (M1:M2, N1:N2)
C     N1,
C     N2
C     X,   (MDIM,NDIM) R   I         Grid coordinates
C     Y
C     XP,      -       R   I         Coordinates of the target point
C     YP
C     M,       -       I   I O       Input indices of starting guess and
C     N                              output indices of "lower left corner"
C                                    of the enclosing cell (or nearest to it)
C     EPS      -       R   I         Search tolerance used by BILINT, q.v.
C     P,       -       R   O         Bilinear interpolation coefficients;
C     Q                              see BILINT for their usage
C     STATUS   -       I   O         0 = match found;
C                                    1 = failed to find match; M, N, P, Q
C                                        represent the smallest mismatch
C                                    (Note that for the sake of efficiency, the
C                                    matrix singularity condition checked for
C                                    by BILINT is not tested for here.)
C
C     PROCEDURES:
C
C        BILINT     Performs test for enclosing cell and returns P, Q
C
C     ENVIRONMENT:
C
C        FORTRAN 77 with minor extensions
C
C     HISTORY:
C
C        11/10/93  M.Wong      Original design and coding.
C        11/24/93  D.Saunders  Arguments changed to allow search of sub-
C                              grids, etc.
C        11/28/93      "       Ensured "lower left corner" indices (not
C                              "upper" or "right" corner corresponding to
C                              P or Q = 1.), as for 1-D utility INTERVAL.
C        08/15/94      "       Egregious omission rectified: P, Q are now
C                              returned for bilinear applications (not re-
C                              quired for the original bicubic applications).
C                              Also: success on the initial BILINT call was
C                              skipping the possible adjustment of 11/28/93.
C        01/23/96      "       Returned the nearest M, N, P, Q to treat possible
C                              difficulties near the edges of (sub)regions.
C                              This can help establish an appropriate EPS too.
C        01/10/98      "       The 11/28/93 adjustment is dubious - should at
C                              least update P or Q via BILINT.
C     AUTHOR:
C
C        Michael D. Wong, Sterling Software, NASA/Ames, Moffett Field, CA
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER    MDIM, NDIM, M1, M2, N1, N2, M, N, STATUS
      REAL       X (MDIM, NDIM), Y (MDIM, NDIM), XP, YP, EPS, P, Q

C     Local constants.

      REAL       ONE, ZERO
      PARAMETER (ONE = 1.E+0, ZERO = 0.E+0)

C     Local variables.

      INTEGER    I, J, L, ICORNER, M1P1, M2M1, M2M2, N1P1, N2M1, N2M2,
     >           NEDGE, ISTART, IEND, JEND, JSTART, MBEST, NBEST
      REAL       EBEST, ERR, PBEST, PERR, QBEST, QERR
      LOGICAL    ITEST, JTEST, TEST

C     Execution.


C     Begin with the initial guess.

      CALL BILINT (MDIM, NDIM, X, Y, XP, YP, M, N, EPS, P, Q, STATUS)

      IF (STATUS .EQ. 0) GO TO 90    ! M, N, P, Q may need adjusting

      IF (P .LT. ZERO) THEN
         PERR = -P
      ELSE IF (P .GT. ONE) THEN
         PERR = P - ONE
      ELSE
         PERR = ZERO
      END IF

      IF (Q .LT. ZERO) THEN
         QERR = -Q
      ELSE IF (Q .GT. ONE) THEN
         QERR = Q - ONE
      ELSE
         QERR = ZERO
      END IF

      EBEST = PERR + QERR
      PBEST = P
      QBEST = Q
      MBEST = M
      NBEST = N

C     Enter the rectangular layer search pattern.

      M1P1 = M1 + 1
      M2M1 = M2 - 1
      M2M2 = M2M1 - 1
      N1P1 = N1 + 1
      N2M1 = N2 - 1
      N2M2 = N2M1 - 1

      L = 0
   10 CONTINUE

         L = L + 1  ! Increment "ripple" level.

C        Define and test the corner points.

         DO ICORNER = 1, 4

            ITEST = .TRUE.
            JTEST = .TRUE.

            IF (ICORNER .EQ. 1) THEN

               I = M + L
               J = N + L

               IF (I .GT. M2M1) THEN
                  I     = M2M1
                  ITEST = .FALSE.
               END IF

               IF (J .GT. N2M1) THEN
                  J     = N2M1
                  JTEST = .FALSE.
               END IF

            ELSE IF (ICORNER .EQ. 2) THEN

               I = M + L
               J = N - L

               IF (I .GT. M2M1) THEN
                  I     = M2M1
                  ITEST = .FALSE.
               END IF

               IF (J .LT. N1) THEN
                  J     = N1
                  JTEST = .FALSE.
               END IF

            ELSE IF (ICORNER .EQ. 3) THEN

               I = M - L
               J = N - L

               IF (I .LT. M1) THEN
                  I     = M1
                  ITEST = .FALSE.
               END IF

               IF (J .LT. N1) THEN
                  J     = N1
                  JTEST = .FALSE.
               END IF

            ELSE  ! ICORNER = 4

               I = M - L
               J = N + L

               IF (I .LT. M1) THEN
                  I     = M1
                  ITEST = .FALSE.
               END IF

               IF (J .GT. N2M1) THEN
                  J     = N2M1
                  JTEST = .FALSE.
               END IF

            END IF

C           Avoid repeated tests.

            IF (ITEST .OR. JTEST) THEN

               CALL BILINT (MDIM, NDIM, X, Y, XP, YP, I, J, EPS,
     >                      P, Q, STATUS)

               IF (STATUS .EQ. 0) THEN
                  M = I
                  N = J
                  GO TO 90
               END IF

               IF (P .LT. ZERO) THEN
                  PERR = -P
               ELSE IF (P .GT. ONE) THEN
                  PERR = P - ONE
               ELSE
                  PERR = ZERO
               END IF

               IF (Q .LT. ZERO) THEN
                  QERR = -Q
               ELSE IF (Q .GT. ONE) THEN
                  QERR = Q - ONE
               ELSE
                  QERR = ZERO
               END IF

               ERR = PERR + QERR
               IF (ERR .LT. EBEST) THEN
                  EBEST = ERR
                  PBEST = P
                  QBEST = Q
                  MBEST = I
                  NBEST = J
               END IF
            END IF
         END DO

C        Cycle through interior points of the connector lines one time only.

         DO NEDGE = 1, 4

            IF (NEDGE .EQ. 1) THEN

               ISTART = M + L
               IEND   = ISTART
               JSTART = N - L + 1
               JEND   = N + L - 1
               TEST   = IEND .LE. M2M1

            ELSE IF (NEDGE .EQ. 2) THEN

               ISTART = M - L + 1
               IEND   = M + L - 1
               JSTART = N - L
               JEND   = JSTART
               TEST   = JSTART .GE. N1

            ELSE IF (NEDGE .EQ. 3) THEN

               ISTART = M - L
               IEND   = ISTART
               JSTART = N - L + 1
               JEND   = N + L - 1
               TEST   = ISTART .GE. M1

            ELSE ! Edge 4.

               ISTART = M - L + 1
               IEND   = M + L - 1
               JSTART = N + L
               JEND   = JSTART
               TEST   = JEND .LE. N2M1

            END IF

            IF (TEST) THEN

C              Limit the search to the interior points of the edge.

               IF (ISTART .NE. IEND) THEN

                  IF (ISTART .LT. M1P1) ISTART = M1P1
                  IF (IEND   .GT. M2M2) IEND   = M2M2

               ELSE IF (JSTART .NE. JEND) THEN

                  IF (JSTART .LT. N1P1) JSTART = N1P1
                  IF (JEND   .GT. N2M2) JEND   = N2M2

               END IF

               DO I = ISTART, IEND
                  DO J = JSTART, JEND

                     CALL BILINT (MDIM, NDIM, X, Y, XP, YP, I, J,
     >                            EPS, P, Q, STATUS)

                     IF (STATUS .EQ. 0) THEN
                        M = I
                        N = J
                        GO TO 90
                     END IF

                     IF (P .LT. ZERO) THEN
                        PERR = -P
                     ELSE IF (P .GT. ONE) THEN
                        PERR = P - ONE
                     ELSE
                        PERR = ZERO
                     END IF

                     IF (Q .LT. ZERO) THEN
                        QERR = -Q
                     ELSE IF (Q .GT. ONE) THEN
                        QERR = Q - ONE
                     ELSE
                        QERR = ZERO
                     END IF

                     ERR = PERR + QERR
                     IF (ERR .LT. EBEST) THEN
                        EBEST = ERR
                        PBEST = P
                        QBEST = Q
                        MBEST = I
                        NBEST = J
                     END IF
                  END DO
               END DO
            END IF
         END DO

C        Check to see if the entire grid has been tested.

         IF (M - L .GT. M1 .OR. M + L .LT. M2M1 .OR.
     >       N - L .GT. N1 .OR. N + L .LT. N2M1)
     >GO TO 10  ! Next search level.


C     Dropped through - entire (sub)grid was tested without success.
C     Return the nearest result - may be adequate near a boundary:

      M = MBEST
      N = NBEST
      P = PBEST
      Q = QBEST


   90 CONTINUE

C     Ensure "lower left" as for INTERVAL's 1-D search:

      IF (P .GE. ONE) THEN
         IF (M + 1 .LT. M2) THEN
            M = M + 1
            CALL BILINT (MDIM, NDIM, X, Y, XP, YP, M, N, EPS, P, Q,
     >                   STATUS)
         END IF
      END IF

      IF (Q .GE. ONE) THEN
         IF (N + 1 .LT. N2) THEN
            N = N + 1
            CALL BILINT (MDIM, NDIM, X, Y, XP, YP, M, N, EPS, P, Q,
     >                   STATUS)
         END IF
      END IF

      RETURN
      END
C+------------------------------------------------------------------------------
C
      SUBROUTINE WARP3D (IMIN, IMAX, JMIN, JMAX, KMIN, KMAX,
     >                   I1, I2, J1, J2, K1, K2, X0, Y0, Z0, S0,
     >                   DFACEI, DFACEJ, DFACEK, X, Y, Z)
C
C  ONE-LINER: Perturb interior of a 3-space grid given new faces
C
C  DESCRIPTION:
C
C        WARP3D perturbs the interior points of a 3-space grid block
C     given the original grid and perturbed faces. All six faces are
C     assumed to be perturbed (though some may of course be fixed).
C     If all edges are unperturbed, considerably less work is needed,
C     but here some corners and/or edges are assumed to have been
C     perturbed, requiring a three-stage algorithm. Stage 1 accounts
C     for corner motion, stage 2 for edge motion, and stage 3 matches
C     the desired faces.
C
C        The perturbed faces should be input as faces of the desired
C     output grid.  The original relative arc-length increments in
C     each index direction should also be input.  See PARAMXYZ for
C     setting them up in preparation for multiple perturbations.
C
C  ENVIRONMENT:  FORTRAN 77 with minor extensions
C
C  HISTORY:
C
C     03/04/94  DAS/JJR  WARP3D written to perturb a single 3-space line,
C                        for James's wing/body design code.
C     12/07/94   "   "   Rewritten to apply the idea at a higher level,
C                        for efficient perturbation of a whole block.
C     02/01/96    DAS    Use on a subgrid with full-grid S0 must transform S0.
C     02/10/96     "     It's now a 3-stage algorithm (corners, edges, faces)
C                        where originally the first 2 stages were blurred.
C     02/12/96     "     Storing all the intermediate face perturbations allows
C                        a single pass through the volume points.
C     03/29/96  JJR/DAS  Stage 2 is done more consistently now by distinguishing
C                        the two types of perturbation affecting each face.
C                        This required adding DFACE*(*,*,*,*,4).
C
C  AUTHOR:  David Saunders/James Reuther, NASA Ames, Mt. View, CA.
C
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER IMIN, IMAX, JMIN, JMAX, KMIN, KMAX     ! I Grid array dimensions.

      INTEGER I1, I2, J1, J2, K1, K2                 ! I Define active volume.

      REAL    X0 (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX),  ! I Original grid coords.
     >        Y0 (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX),
     >        Z0 (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX)

      REAL    S0 (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX, 3)! I Relative arc-lengths
                                                     !   for the I,J,K lines
                                                     !   see PARAMXYZ.  If from
                                                     !   a grid larger than the
                                                     !   active subgrid, S0 is
                                                     !   transformed here.

      REAL    DFACEI (3, JMIN:JMAX, KMIN:KMAX, 2, 4),! S For face perturbations:
     >        DFACEJ (3, IMIN:IMAX, KMIN:KMAX, 2, 4),!   DFACEI(1:3,*,*,1,M) =
     >        DFACEK (3, IMIN:IMAX, JMIN:JMAX, 2, 4) !   dX,dY,dZ on the I = I1
                                                     !   face for stage M, etc.;
                                                     !   M = 1, 2a, 2b, 3.

      REAL    X (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX),   !I/O Grid coords.: input
     >        Y (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX),   !   with faces perturbed;
     >        Z (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX)    !   output fully perturbed.

C-------------------------------------------------------------------------------

C     Local constants.

      REAL    EPS, ONE

      PARAMETER (EPS = 1.E-8, ONE = 1.E+0) ! EPS safeguards a divide by zero -
                                           ! presumably only if result is zero.
C     Local variables.

      INTEGER I, J, K, L
      REAL    DEL1, DEL2, DELI, DELJ, DELK,
     >        DELIJ, DELIK, DELJI, DELJK, DELKI, DELKJ,
     >        WTI1, WTI2, WTJ1, WTJ2, WTK1, WTK2


C     Execution.
C     ----------

C     Set up the face perturbations for all three stages:

C     Stage 1:  Corner motion (only). This stage is reusable by WARPQ3D.
C     --------

C     I = I1 and I2 intermediate faces.

      CALL DELQ3D (IMIN, IMAX, JMIN, JMAX, KMIN, KMAX,
     >             I1, I1, J1, J2, K1, K2, X0, Y0, Z0, S0,
     >             DFACEI (1,JMIN,KMIN,1,1), DFACEJ, DFACEK, X, Y, Z)

      CALL DELQ3D (IMIN, IMAX, JMIN, JMAX, KMIN, KMAX,
     >             I2, I2, J1, J2, K1, K2, X0, Y0, Z0, S0,
     >             DFACEI (1,JMIN,KMIN,2,1), DFACEJ, DFACEK, X, Y, Z)

C     J = J1 and J2 intermediate faces.

      CALL DELQ3D (IMIN, IMAX, JMIN, JMAX, KMIN, KMAX,
     >             I1, I2, J1, J1, K1, K2, X0, Y0, Z0, S0,
     >             DFACEI, DFACEJ (1,IMIN,KMIN,1,1), DFACEK, X, Y, Z)

      CALL DELQ3D (IMIN, IMAX, JMIN, JMAX, KMIN, KMAX,
     >             I1, I2, J2, J2, K1, K2, X0, Y0, Z0, S0,
     >             DFACEI, DFACEJ (1,IMIN,KMIN,2,1), DFACEK, X, Y, Z)

C     K = K1 and K2 intermediate faces.

      CALL DELQ3D (IMIN, IMAX, JMIN, JMAX, KMIN, KMAX,
     >             I1, I2, J1, J2, K1, K1, X0, Y0, Z0, S0,
     >             DFACEI, DFACEJ, DFACEK (1,IMIN,JMIN,1,1), X, Y, Z)

      CALL DELQ3D (IMIN, IMAX, JMIN, JMAX, KMIN, KMAX,
     >             I1, I2, J1, J2, K2, K2, X0, Y0, Z0, S0,
     >             DFACEI, DFACEJ, DFACEK (1,IMIN,JMIN,2,1), X, Y, Z)


C     Stage 2:  Handle edge motion from above interim edges to final edges.
C     --------

C     James's insight here: consider the 4 faces affected by the motion of
C     4 edges in a given (index) direction; furthermore, keep the two
C     directions of the resulting face perturbations separated for proper
C     weighted combination during the final interpolation into the interior.
C     The 3 index directions are then independent of each other (added).

C     I1 and I2 faces:

      L = 1
      DO I = I1, I2, I2 - I1

C        K = K1 and K2 edge perturbations (J direction edges):

         DO J = J1 + 1, J2 - 1
            DFACEI(1,J,K1,L,2) = X(I,J,K1)-X0(I,J,K1)-DFACEI(1,J,K1,L,1)
            DFACEI(2,J,K1,L,2) = Y(I,J,K1)-Y0(I,J,K1)-DFACEI(2,J,K1,L,1)
            DFACEI(3,J,K1,L,2) = Z(I,J,K1)-Z0(I,J,K1)-DFACEI(3,J,K1,L,1)
            DFACEI(1,J,K2,L,2) = X(I,J,K2)-X0(I,J,K2)-DFACEI(1,J,K2,L,1)
            DFACEI(2,J,K2,L,2) = Y(I,J,K2)-Y0(I,J,K2)-DFACEI(2,J,K2,L,1)
            DFACEI(3,J,K2,L,2) = Z(I,J,K2)-Z0(I,J,K2)-DFACEI(3,J,K2,L,1)
         END DO

C        J = J1 and J2 edge perturbations (K direction edges):

         DO K = K1 + 1, K2 - 1
            DFACEI(1,J1,K,L,3) = X(I,J1,K)-X0(I,J1,K)-DFACEI(1,J1,K,L,1)
            DFACEI(2,J1,K,L,3) = Y(I,J1,K)-Y0(I,J1,K)-DFACEI(2,J1,K,L,1)
            DFACEI(3,J1,K,L,3) = Z(I,J1,K)-Z0(I,J1,K)-DFACEI(3,J1,K,L,1)
            DFACEI(1,J2,K,L,3) = X(I,J2,K)-X0(I,J2,K)-DFACEI(1,J2,K,L,1)
            DFACEI(2,J2,K,L,3) = Y(I,J2,K)-Y0(I,J2,K)-DFACEI(2,J2,K,L,1)
            DFACEI(3,J2,K,L,3) = Z(I,J2,K)-Z0(I,J2,K)-DFACEI(3,J2,K,L,1)
         END DO

C        Interpolate stage 2 interior points for this I face, keeping the
C        two index directions separated.

         DO K = K1 + 1, K2 - 1
            DO J = J1 + 1, J2 - 1
               WTJ2 = (S0 (I, J,  K, 2) - S0 (I, J1, K, 2)) /
     >                (S0 (I, J2, K, 2) - S0 (I, J1, K, 2))
               WTJ1 = ONE - WTJ2
               WTK2 = (S0 (I, J, K,  3) - S0 (I, J, K1, 3)) /
     >                (S0 (I, J, K2, 3) - S0 (I, J, K1, 3))
               WTK1 = ONE - WTK2

               DFACEI(1,J,K,L,3) = WTJ1*DFACEI(1,J1,K,L,3) +
     >                             WTJ2*DFACEI(1,J2,K,L,3)
               DFACEI(2,J,K,L,3) = WTJ1*DFACEI(2,J1,K,L,3) +
     >                             WTJ2*DFACEI(2,J2,K,L,3)
               DFACEI(3,J,K,L,3) = WTJ1*DFACEI(3,J1,K,L,3) +
     >                             WTJ2*DFACEI(3,J2,K,L,3)

               DFACEI(1,J,K,L,2) = WTK1*DFACEI(1,J,K1,L,2) +
     >                             WTK2*DFACEI(1,J,K2,L,2)
               DFACEI(2,J,K,L,2) = WTK1*DFACEI(2,J,K1,L,2) +
     >                             WTK2*DFACEI(2,J,K2,L,2)
               DFACEI(3,J,K,L,2) = WTK1*DFACEI(3,J,K1,L,2) +
     >                             WTK2*DFACEI(3,J,K2,L,2)
            END DO
         END DO
         L = 2
      END DO

C     J1 and J2 faces, stage 2:

      L = 1
      DO J = J1, J2, J2 - J1

C        K = K1 and K2 edge perturbations (I direction):

         DO I = I1 + 1, I2-1
            DFACEJ(1,I,K1,L,2) = X(I,J,K1)-X0(I,J,K1)-DFACEJ(1,I,K1,L,1)
            DFACEJ(2,I,K1,L,2) = Y(I,J,K1)-Y0(I,J,K1)-DFACEJ(2,I,K1,L,1)
            DFACEJ(3,I,K1,L,2) = Z(I,J,K1)-Z0(I,J,K1)-DFACEJ(3,I,K1,L,1)
            DFACEJ(1,I,K2,L,2) = X(I,J,K2)-X0(I,J,K2)-DFACEJ(1,I,K2,L,1)
            DFACEJ(2,I,K2,L,2) = Y(I,J,K2)-Y0(I,J,K2)-DFACEJ(2,I,K2,L,1)
            DFACEJ(3,I,K2,L,2) = Z(I,J,K2)-Z0(I,J,K2)-DFACEJ(3,I,K2,L,1)
         END DO

C        I = I1 and I2 edge perturbations (K direction):

         DO K = K1 + 1, K2 - 1
            DFACEJ(1,I1,K,L,3) = X(I1,J,K)-X0(I1,J,K)-DFACEJ(1,I1,K,L,1)
            DFACEJ(2,I1,K,L,3) = Y(I1,J,K)-Y0(I1,J,K)-DFACEJ(2,I1,K,L,1)
            DFACEJ(3,I1,K,L,3) = Z(I1,J,K)-Z0(I1,J,K)-DFACEJ(3,I1,K,L,1)
            DFACEJ(1,I2,K,L,3) = X(I2,J,K)-X0(I2,J,K)-DFACEJ(1,I2,K,L,1)
            DFACEJ(2,I2,K,L,3) = Y(I2,J,K)-Y0(I2,J,K)-DFACEJ(2,I2,K,L,1)
            DFACEJ(3,I2,K,L,3) = Z(I2,J,K)-Z0(I2,J,K)-DFACEJ(3,I2,K,L,1)
         END DO

C        Interpolate stage 2 interior points for this J face.

         DO K = K1 + 1, K2 - 1
            DO I = I1 + 1, I2 - 1
               WTI2 = (S0 (I,  J, K, 1) - S0 (I1, J, K, 1)) /
     >                (S0 (I2, J, K, 1) - S0 (I1, J, K, 1))
               WTI1 = ONE - WTI2
               WTK2 = (S0 (I, J, K,  3) - S0 (I, J, K1, 3)) /
     >                (S0 (I, J, K2, 3) - S0 (I, J, K1, 3))
               WTK1 = ONE - WTK2

               DFACEJ(1,I,K,L,3) = WTI1*DFACEJ(1,I1,K,L,3) +
     >                             WTI2*DFACEJ(1,I2,K,L,3)
               DFACEJ(2,I,K,L,3) = WTI1*DFACEJ(2,I1,K,L,3) +
     >                             WTI2*DFACEJ(2,I2,K,L,3)
               DFACEJ(3,I,K,L,3) = WTI1*DFACEJ(3,I1,K,L,3) +
     >                             WTI2*DFACEJ(3,I2,K,L,3)

               DFACEJ(1,I,K,L,2) = WTK1*DFACEJ(1,I,K1,L,2) +
     >                             WTK2*DFACEJ(1,I,K2,L,2)
               DFACEJ(2,I,K,L,2) = WTK1*DFACEJ(2,I,K1,L,2) +
     >                             WTK2*DFACEJ(2,I,K2,L,2)
               DFACEJ(3,I,K,L,2) = WTK1*DFACEJ(3,I,K1,L,2) +
     >                             WTK2*DFACEJ(3,I,K2,L,2)
            END DO
         END DO
         L = 2
      END DO

C     K1 and K2 faces, stage 2:

      L = 1
      DO K = K1, K2, K2 - K1

C        J = J1 and J2 edge perturbations (I direction):

         DO I = I1 + 1, I2 - 1
            DFACEK(1,I,J1,L,2) = X(I,J1,K)-X0(I,J1,K)-DFACEK(1,I,J1,L,1)
            DFACEK(2,I,J1,L,2) = Y(I,J1,K)-Y0(I,J1,K)-DFACEK(2,I,J1,L,1)
            DFACEK(3,I,J1,L,2) = Z(I,J1,K)-Z0(I,J1,K)-DFACEK(3,I,J1,L,1)
            DFACEK(1,I,J2,L,2) = X(I,J2,K)-X0(I,J2,K)-DFACEK(1,I,J2,L,1)
            DFACEK(2,I,J2,L,2) = Y(I,J2,K)-Y0(I,J2,K)-DFACEK(2,I,J2,L,1)
            DFACEK(3,I,J2,L,2) = Z(I,J2,K)-Z0(I,J2,K)-DFACEK(3,I,J2,L,1)
         END DO

C        I = I1 and I2 edge perturbations (J direction):

         DO J = J1 + 1, J2 - 1
            DFACEK(1,I1,J,L,3) = X(I1,J,K)-X0(I1,J,K)-DFACEK(1,I1,J,L,1)
            DFACEK(2,I1,J,L,3) = Y(I1,J,K)-Y0(I1,J,K)-DFACEK(2,I1,J,L,1)
            DFACEK(3,I1,J,L,3) = Z(I1,J,K)-Z0(I1,J,K)-DFACEK(3,I1,J,L,1)
            DFACEK(1,I2,J,L,3) = X(I2,J,K)-X0(I2,J,K)-DFACEK(1,I2,J,L,1)
            DFACEK(2,I2,J,L,3) = Y(I2,J,K)-Y0(I2,J,K)-DFACEK(2,I2,J,L,1)
            DFACEK(3,I2,J,L,3) = Z(I2,J,K)-Z0(I2,J,K)-DFACEK(3,I2,J,L,1)
         END DO

C        Interpolate stage 2 interior points for this K face.

         DO J = J1 + 1, J2 - 1
            DO I = I1 + 1, I2 - 1
               WTI2 = (S0 (I,  J, K, 1) - S0 (I1, J, K, 1)) /
     >                (S0 (I2, J, K, 1) - S0 (I1, J, K, 1))
               WTI1 = ONE - WTI2
               WTJ2 = (S0 (I, J,  K, 2) - S0 (I, J1, K, 2)) /
     >                (S0 (I, J2, K, 2) - S0 (I, J1, K, 2))
               WTJ1 = ONE - WTJ2

               DFACEK(1,I,J,L,3) = WTI1*DFACEK(1,I1,J,L,3) +
     >                             WTI2*DFACEK(1,I2,J,L,3)
               DFACEK(2,I,J,L,3) = WTI1*DFACEK(2,I1,J,L,3) +
     >                             WTI2*DFACEK(2,I2,J,L,3)
               DFACEK(3,I,J,L,3) = WTI1*DFACEK(3,I1,J,L,3) +
     >                             WTI2*DFACEK(3,I2,J,L,3)

               DFACEK(1,I,J,L,2) = WTJ1*DFACEK(1,I,J1,L,2) +
     >                             WTJ2*DFACEK(1,I,J2,L,2)
               DFACEK(2,I,J,L,2) = WTJ1*DFACEK(2,I,J1,L,2) +
     >                             WTJ2*DFACEK(2,I,J2,L,2)
               DFACEK(3,I,J,L,2) = WTJ1*DFACEK(3,I,J1,L,2) +
     >                             WTJ2*DFACEK(3,I,J2,L,2)
            END DO
         END DO
         L = 2
      END DO


C     Stage 3:  Handle face motion from above interim faces to final faces.
C     --------

      DO K = K1 + 1, K2 - 1
         DO J = J1 + 1, J2 - 1
            DFACEI(1,J,K,1,4) = X(I1,J,K) - X0(I1,J,K) -
     >         DFACEI(1,J,K,1,1) - DFACEI(1,J,K,1,2) - DFACEI(1,J,K,1,3)
            DFACEI(2,J,K,1,4) = Y(I1,J,K) - Y0(I1,J,K) -
     >         DFACEI(2,J,K,1,1) - DFACEI(2,J,K,1,2) - DFACEI(2,J,K,1,3)
            DFACEI(3,J,K,1,4) = Z(I1,J,K) - Z0(I1,J,K) -
     >         DFACEI(3,J,K,1,1) - DFACEI(3,J,K,1,2) - DFACEI(3,J,K,1,3)
            DFACEI(1,J,K,2,4) = X(I2,J,K) - X0(I2,J,K) -
     >         DFACEI(1,J,K,2,1) - DFACEI(1,J,K,2,2) - DFACEI(1,J,K,2,3)
            DFACEI(2,J,K,2,4) = Y(I2,J,K) - Y0(I2,J,K) -
     >         DFACEI(2,J,K,2,1) - DFACEI(2,J,K,2,2) - DFACEI(2,J,K,2,3)
            DFACEI(3,J,K,2,4) = Z(I2,J,K) - Z0(I2,J,K) -
     >         DFACEI(3,J,K,2,1) - DFACEI(3,J,K,2,2) - DFACEI(3,J,K,2,3)
         END DO

         DO I = I1 + 1, I2 - 1
            DFACEJ(1,I,K,1,4) = X(I,J1,K) - X0(I,J1,K) -
     >         DFACEJ(1,I,K,1,1) - DFACEJ(1,I,K,1,2) - DFACEJ(1,I,K,1,3)
            DFACEJ(2,I,K,1,4) = Y(I,J1,K) - Y0(I,J1,K) -
     >         DFACEJ(2,I,K,1,1) - DFACEJ(2,I,K,1,2) - DFACEJ(2,I,K,1,3)
            DFACEJ(3,I,K,1,4) = Z(I,J1,K) - Z0(I,J1,K) -
     >         DFACEJ(3,I,K,1,1) - DFACEJ(3,I,K,1,2) - DFACEJ(3,I,K,1,3)
            DFACEJ(1,I,K,2,4) = X(I,J2,K) - X0(I,J2,K) -
     >         DFACEJ(1,I,K,2,1) - DFACEJ(1,I,K,2,2) - DFACEJ(1,I,K,2,3)
            DFACEJ(2,I,K,2,4) = Y(I,J2,K) - Y0(I,J2,K) -
     >         DFACEJ(2,I,K,2,1) - DFACEJ(2,I,K,2,2) - DFACEJ(2,I,K,2,3)
            DFACEJ(3,I,K,2,4) = Z(I,J2,K) - Z0(I,J2,K) -
     >         DFACEJ(3,I,K,2,1) - DFACEJ(3,I,K,2,2) - DFACEJ(3,I,K,2,3)
         END DO
      END DO

      DO J = J1 + 1, J2 - 1
         DO I = I1 + 1, I2 - 1
            DFACEK(1,I,J,1,4) = X(I,J,K1) - X0(I,J,K1) -
     >         DFACEK(1,I,J,1,1) - DFACEK(1,I,J,1,2) - DFACEK(1,I,J,1,3)
            DFACEK(2,I,J,1,4) = Y(I,J,K1) - Y0(I,J,K1) -
     >         DFACEK(2,I,J,1,1) - DFACEK(2,I,J,1,2) - DFACEK(2,I,J,1,3)
            DFACEK(3,I,J,1,4) = Z(I,J,K1) - Z0(I,J,K1) -
     >         DFACEK(3,I,J,1,1) - DFACEK(3,I,J,1,2) - DFACEK(3,I,J,1,3)
            DFACEK(1,I,J,2,4) = X(I,J,K2) - X0(I,J,K2) -
     >         DFACEK(1,I,J,2,1) - DFACEK(1,I,J,2,2) - DFACEK(1,I,J,2,3)
            DFACEK(2,I,J,2,4) = Y(I,J,K2) - Y0(I,J,K2) -
     >         DFACEK(2,I,J,2,1) - DFACEK(2,I,J,2,2) - DFACEK(2,I,J,2,3)
            DFACEK(3,I,J,2,4) = Z(I,J,K2) - Z0(I,J,K2) -
     >         DFACEK(3,I,J,2,1) - DFACEK(3,I,J,2,2) - DFACEK(3,I,J,2,3)
         END DO
      END DO


C     Perturb the interior volume points.
C     All stages are performed at once via interpolation from the face
C     perturbations stored for each stage.  Note that the three stages
C     accumulate the contributions from the three subscript directions
C     with varying degrees of independence.

      DO K = K1 + 1, K2 - 1
         DO J = J1 + 1, J2 - 1
            DO I = I1 + 1, I2 - 1
               WTI2 = (S0 (I,  J, K, 1) - S0 (I1, J, K, 1)) /
     >                (S0 (I2, J, K, 1) - S0 (I1, J, K, 1))
               WTI1 = ONE - WTI2
               WTJ2 = (S0 (I, J,  K, 2) - S0 (I, J1, K, 2)) /
     >                (S0 (I, J2, K, 2) - S0 (I, J1, K, 2))
               WTJ1 = ONE - WTJ2
               WTK2 = (S0 (I, J, K,  3) - S0 (I, J, K1, 3)) /
     >                (S0 (I, J, K2, 3) - S0 (I, J, K1, 3))
               WTK1 = ONE - WTK2

C              Stage 1:

               DELI = WTI1*DFACEI (1,J,K,1,1) + WTI2*DFACEI (1,J,K,2,1)
               DELJ = WTJ1*DFACEJ (1,I,K,1,1) + WTJ2*DFACEJ (1,I,K,2,1)
               DELK = WTK1*DFACEK (1,I,J,1,1) + WTK2*DFACEK (1,I,J,2,1)

               DEL1 = (ABS(DELI)*DELI + ABS(DELJ)*DELJ + ABS(DELK)*DELK)
     >                / MAX (ABS (DELI) + ABS (DELJ) + ABS (DELK), EPS)

C              Stage 2:

               DELIJ= WTJ1*DFACEJ (1,I,K,1,2) + WTJ2*DFACEJ (1,I,K,2,2)
               DELIK= WTK1*DFACEK (1,I,J,1,2) + WTK2*DFACEK (1,I,J,2,2)
               DELI = (ABS (DELIJ) * DELIJ + ABS (DELIK) * DELIK) /
     >                 MAX (ABS (DELIJ) + ABS (DELIK), EPS)

               DELJI= WTI1*DFACEI (1,J,K,1,2) + WTI2*DFACEI (1,J,K,2,2)
               DELJK= WTK1*DFACEK (1,I,J,1,3) + WTK2*DFACEK (1,I,J,2,3)
               DELJ = (ABS (DELJI) * DELJI + ABS (DELJK) * DELJK) /
     >                 MAX (ABS (DELJI) + ABS (DELJK), EPS)

               DELKI= WTI1*DFACEI (1,J,K,1,3) + WTI2*DFACEI (1,J,K,2,3)
               DELKJ= WTJ1*DFACEJ (1,I,K,1,3) + WTJ2*DFACEJ (1,I,K,2,3)
               DELK = (ABS (DELKI) * DELKI + ABS (DELKJ) * DELKJ) /
     >                 MAX (ABS (DELKI) + ABS (DELKJ), EPS)

               DEL2 = DELI + DELJ + DELK

C              Stage 3:

               DELI = WTI1*DFACEI (1,J,K,1,4) + WTI2*DFACEI (1,J,K,2,4)
               DELJ = WTJ1*DFACEJ (1,I,K,1,4) + WTJ2*DFACEJ (1,I,K,2,4)
               DELK = WTK1*DFACEK (1,I,J,1,4) + WTK2*DFACEK (1,I,J,2,4)

               X (I, J, K) = X0 (I, J, K) + DEL1 + DEL2 + DELI+DELJ+DELK

C              Repeat all stages for Y perturbations:

               DELI = WTI1*DFACEI (2,J,K,1,1) + WTI2*DFACEI (2,J,K,2,1)
               DELJ = WTJ1*DFACEJ (2,I,K,1,1) + WTJ2*DFACEJ (2,I,K,2,1)
               DELK = WTK1*DFACEK (2,I,J,1,1) + WTK2*DFACEK (2,I,J,2,1)

               DEL1 = (ABS(DELI)*DELI + ABS(DELJ)*DELJ + ABS(DELK)*DELK)
     >                / MAX (ABS (DELI) + ABS (DELJ) + ABS (DELK), EPS)

               DELIJ= WTJ1*DFACEJ (2,I,K,1,2) + WTJ2*DFACEJ (2,I,K,2,2)
               DELIK= WTK1*DFACEK (2,I,J,1,2) + WTK2*DFACEK (2,I,J,2,2)
               DELI = (ABS (DELIJ) * DELIJ + ABS (DELIK) * DELIK) /
     >                 MAX (ABS (DELIJ) + ABS (DELIK), EPS)

               DELJI= WTI1*DFACEI (2,J,K,1,2) + WTI2*DFACEI (2,J,K,2,2)
               DELJK= WTK1*DFACEK (2,I,J,1,3) + WTK2*DFACEK (2,I,J,2,3)
               DELJ = (ABS (DELJI) * DELJI + ABS (DELJK) * DELJK) /
     >                 MAX (ABS (DELJI) + ABS (DELJK), EPS)

               DELKI= WTI1*DFACEI (2,J,K,1,3) + WTI2*DFACEI (2,J,K,2,3)
               DELKJ= WTJ1*DFACEJ (2,I,K,1,3) + WTJ2*DFACEJ (2,I,K,2,3)
               DELK = (ABS (DELKI) * DELKI + ABS (DELKJ) * DELKJ) /
     >                 MAX (ABS (DELKI) + ABS (DELKJ), EPS)

               DEL2 = DELI + DELJ + DELK

               DELI = WTI1*DFACEI (2,J,K,1,4) + WTI2*DFACEI (2,J,K,2,4)
               DELJ = WTJ1*DFACEJ (2,I,K,1,4) + WTJ2*DFACEJ (2,I,K,2,4)
               DELK = WTK1*DFACEK (2,I,J,1,4) + WTK2*DFACEK (2,I,J,2,4)

               Y (I, J, K) = Y0 (I, J, K) + DEL1 + DEL2 + DELI+DELJ+DELK

C              ... and for Z perturbations:

               DELI = WTI1*DFACEI (3,J,K,1,1) + WTI2*DFACEI (3,J,K,2,1)
               DELJ = WTJ1*DFACEJ (3,I,K,1,1) + WTJ2*DFACEJ (3,I,K,2,1)
               DELK = WTK1*DFACEK (3,I,J,1,1) + WTK2*DFACEK (3,I,J,2,1)

               DEL1 = (ABS(DELI)*DELI + ABS(DELJ)*DELJ + ABS(DELK)*DELK)
     >                / MAX (ABS (DELI) + ABS (DELJ) + ABS (DELK), EPS)

               DELIJ= WTJ1*DFACEJ (3,I,K,1,2) + WTJ2*DFACEJ (3,I,K,2,2)
               DELIK= WTK1*DFACEK (3,I,J,1,2) + WTK2*DFACEK (3,I,J,2,2)
               DELI = (ABS (DELIJ) * DELIJ + ABS (DELIK) * DELIK) /
     >                 MAX (ABS (DELIJ) + ABS (DELIK), EPS)

               DELJI= WTI1*DFACEI (3,J,K,1,2) + WTI2*DFACEI (3,J,K,2,2)
               DELJK= WTK1*DFACEK (3,I,J,1,3) + WTK2*DFACEK (3,I,J,2,3)
               DELJ = (ABS (DELJI) * DELJI + ABS (DELJK) * DELJK) /
     >                 MAX (ABS (DELJI) + ABS (DELJK), EPS)

               DELKI= WTI1*DFACEI (3,J,K,1,3) + WTI2*DFACEI (3,J,K,2,3)
               DELKJ= WTJ1*DFACEJ (3,I,K,1,3) + WTJ2*DFACEJ (3,I,K,2,3)
               DELK = (ABS (DELKI) * DELKI + ABS (DELKJ) * DELKJ) /
     >                 MAX (ABS (DELKI) + ABS (DELKJ), EPS)

               DEL2 = DELI + DELJ + DELK

               DELI = WTI1*DFACEI (3,J,K,1,4) + WTI2*DFACEI (3,J,K,2,4)
               DELJ = WTJ1*DFACEJ (3,I,K,1,4) + WTJ2*DFACEJ (3,I,K,2,4)
               DELK = WTK1*DFACEK (3,I,J,1,4) + WTK2*DFACEK (3,I,J,2,4)

               Z (I, J, K) = Z0 (I, J, K) + DEL1 + DEL2 + DELI+DELJ+DELK
            END DO
         END DO
      END DO

      RETURN
      END
C+------------------------------------------------------------------------------
C
      SUBROUTINE WARPQ3D (IMIN, IMAX, JMIN, JMAX, KMIN, KMAX,
     >                    I1, I2, J1, J2, K1, K2, X0, Y0, Z0, S0,
     >                    DFACEI, DFACEJ, DFACEK, X, Y, Z)
C
C  ONE-LINER: Perturb interior of a 3-space grid face given new edges
C
C  DESCRIPTION:
C
C        WARPQ3D is the quasi-3D form of WARP2D.  It perturbs the
C     interior of one face of a 3-space grid block given perturbed
C     edges of that face, which is indicated by one pair of equal index
C     arguments.  (E.g.: I2 = I1 means a face in the J/K subspace.)
C
C        The two-stage algorithm uses an intermediate perturbation to
C     account for any corner motion, involving edges derived from the
C     original edges, then a second perturbation to account for any
C     differences between the intermediate edges and the specified
C     new edges.
C
C        The perturbed edges should be input as edges of the desired
C     output face.  The original relative arc-length increments in
C     each index direction should also be input.  See PARAMXYZ for
C     setting them up in preparation for multiple perturbations.
C
C  ENVIRONMENT:  FORTRAN 77 with minor extensions
C
C  HISTORY:
C
C     12/16/94  DAS/JDR  Quasi-3D analogue of WARP2D's two-stage algorithm,
C                        with stage 1 of WARP3D's algorithm in mind.
C     12/23/94    DAS    Lots of code replaced by DELQ3D, which is the form
C                        of WARPQ3D needed by WARP3D.
C     02/08/96     "     DELQ3D does stage 1 only now (all that WARP3D needs).
C                        S0(*) values are renormalized before use, in case
C                        the routine is being applied to a subgrid.
C
C  AUTHOR:  David Saunders/James Reuther, NASA Ames, Mt. View, CA.
C
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER IMIN, IMAX, JMIN, JMAX, KMIN, KMAX      ! I Grid array dimensions.

      INTEGER I1, I2, J1, J2, K1, K2                  ! I Define active face,
                                                      !   one pair being equal.

      REAL    X0 (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX),   ! I Original face coords.
     >        Y0 (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX),   !   in appropriate places.
     >        Z0 (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX)

      REAL    S0 (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX, 3) ! I Relative arc-lengths
                                                      !   for the I,J,K lines
                                                      !   (see PARAMXYZ), now
                                                      !   renormalized here.

      REAL    DFACEI (3, JMIN:JMAX, KMIN:KMAX),    ! S For face perturbations:
     >        DFACEJ (3, IMIN:IMAX, KMIN:KMAX),    !   DFACEI (1:3, J1:J2, K:K2)
     >        DFACEK (3, IMIN:IMAX, JMIN:JMAX)     !   = dX,dY,dZ on an I face,
                                                   !   etc.
      REAL    X (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX), ! I/O Grid coordinates:
     >        Y (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX), !     new edges of a face in;
     >        Z (IMIN:IMAX, JMIN:JMAX, KMIN:KMAX)  !     full new face out.

C-------------------------------------------------------------------------------

      REAL    ONE
      PARAMETER (ONE = 1.E+0)

C     Local variables.

      INTEGER I, J, K
      REAL    DELI, DELJ, DELK, WTI1, WTI2, WTJ1, WTJ2, WTK1, WTK2

C     Execution.
C     ----------

C     Stage 1:
C     Handle any corner motion by generating an intermediate face with
C     the final corners but otherwise derived from the original edges.
C     Actually, just set up the appropriate face perturbations.

      CALL DELQ3D (IMIN, IMAX, JMIN, JMAX, KMIN, KMAX,
     >             I1, I2, J1, J2, K1, K2, X0, Y0, Z0, S0,
     >             DFACEI, DFACEJ, DFACEK, X, Y, Z)


C     Stage 2:
C     Set up the perturbations from the intermediate edges to the final
C     edges, then interpolate them into the interior points.

      IF (I1 .EQ. I2) THEN          ! I plane case:

         I = I1

C        J = J1 and J2 edge perturbations:

         DO K = K1 + 1, K2 - 1
            DFACEI (1,J1,K) = X (I,J1,K) - X0 (I,J1,K) - DFACEI (1,J1,K)
            DFACEI (2,J1,K) = Y (I,J1,K) - Y0 (I,J1,K) - DFACEI (2,J1,K)
            DFACEI (3,J1,K) = Z (I,J1,K) - Z0 (I,J1,K) - DFACEI (3,J1,K)
            DFACEI (1,J2,K) = X (I,J2,K) - X0 (I,J2,K) - DFACEI (1,J2,K)
            DFACEI (2,J2,K) = Y (I,J2,K) - Y0 (I,J2,K) - DFACEI (2,J2,K)
            DFACEI (3,J2,K) = Z (I,J2,K) - Z0 (I,J2,K) - DFACEI (3,J2,K)
         END DO

C        K = K1 and K2 edge perturbations:

         DO J = J1 + 1, J2 - 1
            DFACEI (1,J,K1) = X (I,J,K1) - X0 (I,J,K1) - DFACEI (1,J,K1)
            DFACEI (2,J,K1) = Y (I,J,K1) - Y0 (I,J,K1) - DFACEI (2,J,K1)
            DFACEI (3,J,K1) = Z (I,J,K1) - Z0 (I,J,K1) - DFACEI (3,J,K1)
            DFACEI (1,J,K2) = X (I,J,K2) - X0 (I,J,K2) - DFACEI (1,J,K2)
            DFACEI (2,J,K2) = Y (I,J,K2) - Y0 (I,J,K2) - DFACEI (2,J,K2)
            DFACEI (3,J,K2) = Z (I,J,K2) - Z0 (I,J,K2) - DFACEI (3,J,K2)
         END DO

C        Interior points: accumulate the (independent) contributions.

         DO K = K1 + 1, K2 - 1
            DO J = J1 + 1, J2 - 1
               WTJ2 = (S0 (I, J,  K, 2) - S0 (I, J1, K, 2)) /
     >                (S0 (I, J2, K, 2) - S0 (I, J1, K, 2))
               WTJ1 = ONE - WTJ2
               WTK2 = (S0 (I, J, K,  3) - S0 (I, J, K1, 3)) /
     >                (S0 (I, J, K2, 3) - S0 (I, J, K1, 3))
               WTK1 = ONE - WTK2

               DELJ = WTJ1 * DFACEI (1, J1, K) + WTJ2 * DFACEI (1, J2,K)
               DELK = WTK1 * DFACEI (1, J, K1) + WTK2 * DFACEI (1, J,K2)

               X (I, J, K) = X0 (I, J, K) + DFACEI (1,J,K) + DELJ + DELK

               DELJ = WTJ1 * DFACEI (2, J1, K) + WTJ2 * DFACEI (2, J2,K)
               DELK = WTK1 * DFACEI (2, J, K1) + WTK2 * DFACEI (2, J,K2)

               Y (I, J, K) = Y0 (I, J, K) + DFACEI (2,J,K) + DELJ + DELK

               DELJ = WTJ1 * DFACEI (3, J1, K) + WTJ2 * DFACEI (3, J2,K)
               DELK = WTK1 * DFACEI (3, J, K1) + WTK2 * DFACEI (3, J,K2)

               Z (I, J, K) = Z0 (I, J, K) + DFACEI (3,J,K) + DELJ + DELK
            END DO
         END DO

      ELSE IF (J1 .EQ. J2) THEN     ! J plane case:

         J = J1

C        I = I1 and I2 edge perturbations:

         DO K = K1 + 1, K2 - 1
            DFACEJ (1,I1,K) = X (I1,J,K) - X0 (I1,J,K) - DFACEJ (1,I1,K)
            DFACEJ (2,I1,K) = Y (I1,J,K) - Y0 (I1,J,K) - DFACEJ (2,I1,K)
            DFACEJ (3,I1,K) = Z (I1,J,K) - Z0 (I1,J,K) - DFACEJ (3,I1,K)
            DFACEJ (1,I2,K) = X (I2,J,K) - X0 (I2,J,K) - DFACEJ (1,I2,K)
            DFACEJ (2,I2,K) = Y (I2,J,K) - Y0 (I2,J,K) - DFACEJ (2,I2,K)
            DFACEJ (3,I2,K) = Z (I2,J,K) - Z0 (I2,J,K) - DFACEJ (3,I2,K)
         END DO

C        K = K1 and K2 edge perturbations:

         DO I = I1 + 1, I2 - 1
            DFACEJ (1,I,K1) = X (I,J,K1) - X0 (I,J,K1) - DFACEJ (1,I,K1)
            DFACEJ (2,I,K1) = Y (I,J,K1) - Y0 (I,J,K1) - DFACEJ (2,I,K1)
            DFACEJ (3,I,K1) = Z (I,J,K1) - Z0 (I,J,K1) - DFACEJ (3,I,K1)
            DFACEJ (1,I,K2) = X (I,J,K2) - X0 (I,J,K2) - DFACEJ (1,I,K2)
            DFACEJ (2,I,K2) = Y (I,J,K2) - Y0 (I,J,K2) - DFACEJ (2,I,K2)
            DFACEJ (3,I,K2) = Z (I,J,K2) - Z0 (I,J,K2) - DFACEJ (3,I,K2)
         END DO

C        Interior points: accumulate the (independent) contributions.

         DO K = K1 + 1, K2 - 1
            DO I = I1 + 1, I2 - 1
               WTI2 = (S0 (I,  J, K, 1) - S0 (I1, J, K, 1)) /
     >                (S0 (I2, J, K, 1) - S0 (I1, J, K, 1))
               WTI1 = ONE - WTI2
               WTK2 = (S0 (I, J, K,  3) - S0 (I, J, K1, 3)) /
     >                (S0 (I, J, K2, 3) - S0 (I, J, K1, 3))
               WTK1 = ONE - WTK2

               DELI = WTI1 * DFACEJ (1, I1, K) + WTI2 * DFACEJ (1, I2,K)
               DELK = WTK1 * DFACEJ (1, I, K1) + WTK2 * DFACEJ (1, I,K2)

               X (I, J, K) = X0 (I, J, K) + DFACEJ (1,I,K) + DELI + DELK

               DELI = WTI1 * DFACEJ (2, I1, K) + WTI2 * DFACEJ (2, I2,K)
               DELK = WTK1 * DFACEJ (2, I, K1) + WTK2 * DFACEJ (2, I,K2)

               Y (I, J, K) = Y0 (I, J, K) + DFACEJ (2,I,K) + DELI + DELK

               DELI = WTI1 * DFACEJ (3, I1, K) + WTI2 * DFACEJ (3, I2,K)
               DELK = WTK1 * DFACEJ (3, I, K1) + WTK2 * DFACEJ (3, I,K2)

               Z (I, J, K) = Z0 (I, J, K) + DFACEJ (3,I,K) + DELI + DELK
            END DO
         END DO

      ELSE IF (K1 .EQ. K2) THEN     ! K plane case:

         K = K1

C        I = I1 and I2 edge perturbations:

         DO J = J1 + 1, J2 - 1
            DFACEK (1,I1,J) = X (I1,J,K) - X0 (I1,J,K) - DFACEK (1,I1,J)
            DFACEK (2,I1,J) = Y (I1,J,K) - Y0 (I1,J,K) - DFACEK (2,I1,J)
            DFACEK (3,I1,J) = Z (I1,J,K) - Z0 (I1,J,K) - DFACEK (3,I1,J)
            DFACEK (1,I2,J) = X (I2,J,K) - X0 (I2,J,K) - DFACEK (1,I2,J)
            DFACEK (2,I2,J) = Y (I2,J,K) - Y0 (I2,J,K) - DFACEK (2,I2,J)
            DFACEK (3,I2,J) = Z (I2,J,K) - Z0 (I2,J,K) - DFACEK (3,I2,J)
         END DO

C        J = J1 and J2 edge perturbations:

         DO I = I1 + 1, I2 - 1
            DFACEK (1,I,J1) = X (I,J1,K) - X0 (I,J1,K) - DFACEK (1,I,J1)
            DFACEK (2,I,J1) = Y (I,J1,K) - Y0 (I,J1,K) - DFACEK (2,I,J1)
            DFACEK (3,I,J1) = Z (I,J1,K) - Z0 (I,J1,K) - DFACEK (3,I,J1)
            DFACEK (1,I,J2) = X (I,J2,K) - X0 (I,J2,K) - DFACEK (1,I,J2)
            DFACEK (2,I,J2) = Y (I,J2,K) - Y0 (I,J2,K) - DFACEK (2,I,J2)
            DFACEK (3,I,J2) = Z (I,J2,K) - Z0 (I,J2,K) - DFACEK (3,I,J2)
         END DO

C        Interior points: accumulate the (independent) contributions.

         DO J = J1 + 1, J2 - 1
            DO I = I1 + 1, I2 - 1
               WTI2 = (S0 (I,  J, K, 1) - S0 (I1, J, K, 1)) /
     >                (S0 (I2, J, K, 1) - S0 (I1, J, K, 1))
               WTI1 = ONE - WTI2
               WTJ2 = (S0 (I, J,  K, 2) - S0 (I, J1, K, 2)) /
     >                (S0 (I, J2, K, 2) - S0 (I, J1, K, 2))
               WTJ1 = ONE - WTJ2

               DELI = WTI1 * DFACEK (1, I1, J) + WTI2 * DFACEK (1, I2,J)
               DELJ = WTJ1 * DFACEK (1, I, J1) + WTJ2 * DFACEK (1, I,J2)

               X (I, J, K) = X0 (I, J, K) + DFACEK (1,I,J) + DELI + DELJ

               DELI = WTI1 * DFACEK (2, I1, J) + WTI2 * DFACEK (2, I2,J)
               DELJ = WTJ1 * DFACEK (2, I, J1) + WTJ2 * DFACEK (2, I,J2)

               Y (I, J, K) = Y0 (I, J, K) + DFACEK (2,I,J) + DELI + DELJ

               DELI = WTI1 * DFACEK (3, I1, J) + WTI2 * DFACEK (3, I2,J)
               DELJ = WTJ1 * DFACEK (3, I, J1) + WTJ2 * DFACEK (3, I,J2)

               Z (I, J, K) = Z0 (I, J, K) + DFACEK (3,I,J) + DELI + DELJ
            END DO
         END DO

      END IF

      RETURN
      END
C+----------------------------------------------------------------------
C
      FUNCTION BESSEL (J, H, DEL)
C
C     One-liner: First derivative using central 3-point formula
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        Computes a first derivative approximation using the central
C     3-point formula.  The data must be in the form of arrays containing
C     finite difference interval lengths and 2-point forward difference
C     derivatives.  BESSEL is intended to be used by PLSFIT for determin-
C     ing end conditions on an interval for (non-monotonic) interpolation
C     by piecewise cubics.  See the PLSFIT header for more details.
C
C     Arguments:
C     ----------
C
C     Name    Type/Dimension  I/O/S  Description
C     J       I               I      Indicates at which end of the
C                                    interval the derivative is to be
C                                    estimated. J = 0 means left-hand
C                                    side, J = 1 means right.
C
C     H       R (-1:1)        I      Array of interval lengths. The 0th
C                                    element is the length of the interval
C                                    on which the cubic is to be deter-
C                                    mined.
C
C     DEL     R (-1:1)        I      Array of derivative estimates. The
C                                    0th element is the forward difference
C                                    derivative over the interval on which
C                                    the cubic is to be determined.
C                                     
C     BESSEL  R                 O    The function value is the adjusted
C                                    derivative.
C
C     Notes:
C     ------
C
C     (1)  IMPLICIT NONE is non-standard.
C
C     Author:  Robert Kennelly, Sterling Federal Systems/NASA-Ames
C     -------
C
C     Development history:
C     --------------------
C
C     18 Feb. 1987    RAK    Initial design and coding.
C
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     &   ONE
      PARAMETER
     &  (ONE = 1.0E+0)

C     Arguments.

      INTEGER
     &   J
      REAL
     &   H (-1:1), DEL (-1:1), BESSEL

C     Local variables.

      REAL
     &   WEIGHT

C     Execution.
C     ----------

C     Estimate first derivative on left (J = 0) or right side (J = 1) of
C     an interval.

      WEIGHT = H (J) / (H (J) + H (J - 1))
      BESSEL = WEIGHT * DEL (J - 1) + (ONE - WEIGHT) * DEL (J)

C     Termination.
C     ------------

      RETURN
      END
C+----------------------------------------------------------------------
C
      FUNCTION BRODLIE (J, H, DEL)
C
C     One-liner: First derivative, adjusted for monotonicity
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        BRODLIE is intended to be used by PLSFIT for determining end
C     conditions on an interval for monotonic interpolation by piecewise
C     cubics. The data must be in the form of arrays containing finite
C     difference interval lengths and 2-point forward difference deriva-
C     tives. See the PLSFIT header for more details.
C
C        The method is due to Brodlie, Butland, Carlson, and Fritsch,
C     as referenced in the PLSFIT header.
C
C     Arguments:
C     ----------
C
C     Name    Type/Dimension  I/O/S  Description
C     J       I               I      Indicates at which end of the
C                                    interval the derivative is to be
C                                    estimated. J = 0 means left-hand
C                                    side, J = 1 means right.
C
C     H       R (-1:1)        I      Array of interval lengths. The 0th
C                                    element is the length of the interval
C                                    on which the cubic is to be deter-
C                                    mined.
C
C     DEL     R (-1:1)        I      Array of derivative estimates. The
C                                    0th element is the forward difference
C                                    derivative over the interval on which
C                                    the cubic is to be determined.
C
C     BRODLIE R                 O    The function value is the adjusted
C                                    derivative.
C
C     Notes:
C     ------
C
C     (1)  IMPLICIT NONE and 8-character symbolic names are non-standard.
C
C     Author:  Robert Kennelly, Sterling Federal Systems/NASA-Ames
C     -------
C
C     Development history:
C     --------------------
C
C     18 Feb. 1987    RAK    Initial design and coding.
C     04 Dec. 2002    DAS    SIGN work-around for Intel compiler bug.
C
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     &   ZERO, ONE, THIRD
      PARAMETER
     &  (ZERO   = 0.0E+0,
     &   ONE    = 1.0E+0,
     &   THIRD  = ONE / 3.0E+0)

C     Arguments.

      INTEGER
     &   J
      REAL
     &   BRODLIE, H (-1:1), DEL (-1:1)

C     Local variables.

      REAL
     &   ALPHA, PRODUCT

C     Execution.
C     ----------

C     Compare the algebraic signs of the two DEL's.  Have to test that
C     at least one is positive to avoid a zero denominator (this fancy
C     test permits one term to be zero, but the answer below is zero
C     anyway in these cases).  The trick is to work around the SIGN
C     function, which returns positive even if its 2nd argument is zero.

C**** NO:  SIGN misbehaves on Intel systems when the 2nd argument is zero.

      PRODUCT = DEL (J - 1) * DEL (J)

      IF (PRODUCT == ZERO) THEN

         BRODLIE = ZERO

      ELSE IF (SIGN (ONE, -DEL (J - 1)) .NE. SIGN (ONE, DEL (J))) THEN

C        Form "weighted harmonic mean" of the two finite-difference
C        derivative approximations.  Note that we try to avoid overflow
C        by not multiplying them together directly.

         ALPHA   = THIRD * (ONE + H (J) / (H (J - 1) + H (J)))
         BRODLIE = PRODUCT / (ALPHA * DEL (J) +
     &      (ONE - ALPHA) * DEL (J - 1))
      ELSE

C        The signs differ, so make this point a local extremum.

         BRODLIE = ZERO
      END IF

C     Termination.
C     ------------

      RETURN
      END
C+----------------------------------------------------------------------
C
      FUNCTION BUTLAND (J, H, DEL)
C
C     One-liner: First derivative, non-central 3-point formula, adjusted
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        Computes a first derivative approximation for PLSFIT over an
C     interval at a data boundary, using a modified forward or backward
C     3-point formula.  The data must be in the form of arrays containing
C     finite difference interval lengths and 2-point forward difference
C     derivatives, and the differencing direction is controlled by a flag.
C     See PLSFIT for more details, or THREEPT for the pure 3-pt. formula.
C
C        The "shape preserving adjustments" are from PCHIP, a monotone
C     piecewise cubic interpolation package by F. N. Fritsch.
C
C     Arguments:
C     ----------
C
C     Name    Type/Dimension  I/O/S  Description
C     J       I               I      Indicates at which end of the
C                                    interval the derivative is to be
C                                    estimated. J = 0 means left-hand
C                                    side, J = 1 means right.
C
C     H       R (-1:1)        I      Array of interval lengths. The 0th
C                                    element is the length of the interval
C                                    on which the cubic is to be deter-
C                                    mined.
C
C     DEL     R (-1:1)        I      Array of derivative estimates. The
C                                    0th element is the forward difference
C                                    derivative over the interval on which
C                                    the cubic is to be determined.
C
C     BUTLAND R                 O    The function value is the adjusted
C                                    derivative.
C
C     Environment:  Fortran 90
C     ------------
C
C     Author:  Robert Kennelly, Sterling Federal Systems/NASA-Ames
C     -------
C
C     History:
C     --------
C
C     18 Feb. 1987    RAK    Initial design and coding, as THREEPT.
C     20 June 1991    DAS    Monotonic form renamed BUTLAND; THREEPT
C                            is now the pure 3-point formula.
C     04 Dec. 2002     "     SIGN work-arounds for Intel compiler bug.
C
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     &   J
      REAL
     &   H (-1:1), DEL (-1:1), BUTLAND

C     Local constants.

      REAL, PARAMETER ::
     &   ZERO  = 0.0E+0,
     &   ONE   = 1.0E+0,
     &   THREE = 3.0E+0

C     Local variables.

      INTEGER
     &   STEP
      REAL
     &   DMAX, WEIGHT
      LOGICAL
     &   CONSTRAIN

C     Execution.
C     ----------

C     Estimate first derivative on a left-hand boundary using a 3-point
C     forward difference (STEP = +1), or with a backward difference for
C     the right-hand boundary (STEP = -1).

      STEP = 1 - J - J   ! J here is consistent with related modules.

C     In {H, DEL} form, the derivative looks like a weighted average.

C***  Avoid zero as the second argument of SIGN: Intel compiler misbehaves.

      IF (DEL (0) == ZERO) THEN

         BUTLAND = ZERO

      ELSE ! BUTLAND below cannot be zero

         WEIGHT  = -H (0) / (H (0) + H (STEP))
         BUTLAND = WEIGHT * DEL (STEP) + (ONE - WEIGHT) * DEL (0)

C        Shape-preserving adjustments.  Note that we try to avoid overflow
C        by not multiplying quantities directly.

         IF (SIGN (ONE, BUTLAND) /= SIGN (ONE, DEL (0))) THEN

C           Defer to the estimate closest to the boundary.

            BUTLAND = ZERO

C******  ELSE IF (SIGN (ONE, DEL (0)) .NE. SIGN (ONE, DEL (STEP))) THEN
         ELSE

            IF (DEL (STEP) == ZERO) THEN
               CONSTRAIN = DEL (0) < ZERO
            ELSE
               CONSTRAIN = SIGN (ONE, DEL (0)) /= SIGN (ONE, DEL (STEP))
            END IF

            IF (CONSTRAIN) THEN

C              If the monotonicity switches, may need to bound the estimate.

               DMAX = THREE * DEL (0)
               IF (ABS (BUTLAND) > ABS (DMAX)) BUTLAND = DMAX
            END IF

         END IF

      END IF

C     Termination.
C     ------------

      END
C+----------------------------------------------------------------------
C
      FUNCTION THREEPT (J, H, DEL)
C
C     One-liner: First derivative, non-central 3-point formula
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        Computes a first derivative approximation for PLSFIT over an
C     interval at a data boundary, using a forward or backward 3-point
C     formula.  The data must be in the form of arrays containing finite
C     difference interval lengths and 2-point forward difference deriva-
C     tives, and the differencing direction is controlled by a flag. See
C     PLSFIT for more details.
C
C        See module BUTLAND for a version with "shape-preserving"
C     adjustments.
C
C     Arguments:
C     ----------
C
C     Name    Type/Dimension  I/O/S  Description
C     J       I               I      Indicates at which end of the
C                                    interval the derivative is to be
C                                    estimated. J = 0 means left-hand
C                                    side, J = 1 means right. 
C
C     H       R (-1:1)        I      Array of interval lengths. The 0th
C                                    element is the length of the interval
C                                    on which the cubic is to be deter-
C                                    mined.
C
C     DEL     R (-1:1)        I      Array of derivative estimates. The
C                                    0th element is the forward difference
C                                    derivative over the interval on which
C                                    the cubic is to be determined.
C                                     
C     THREEPT R                 O    The function value is the derivative.
C
C     Environment:  VAX/VMS; FORTRAN 77
C     ------------
C
C     IMPLICIT NONE and 8-character symbolic names are non-standard.
C
C     Author:  Robert Kennelly, Sterling Federal Systems/NASA-Ames
C     -------
C
C     History:
C     --------
C
C     18 Feb. 1987    RAK    Initial design and coding.
C     06 June 1991    DAS    Original THREEPT renamed BUTLAND; THREEPT
C                            now gives unmodified 1-sided 3-pt. results.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     &   J
      REAL
     &   H (-1:1), DEL (-1:1), THREEPT

C     Local constants.

      REAL
     &   ONE
      PARAMETER
     &  (ONE = 1.0E+0)

C     Local variables.

      INTEGER
     &   STEP
      REAL
     &   WEIGHT

C     Execution.
C     ----------

C     Estimate first derivative on a left-hand boundary using a 3-point
C     forward difference (STEP = +1), or with a backward difference for
C     the right-hand boundary (STEP = -1).

      STEP = 1 - J - J   ! J here is consistent with related modules.

C     In {H, DEL} form, the derivative looks like a weighted average.

      WEIGHT  = -H (0) / (H (0) + H (STEP))
      THREEPT = WEIGHT * DEL (STEP) + (ONE - WEIGHT) * DEL (0)

C     Termination.
C     ------------

      RETURN
      END
