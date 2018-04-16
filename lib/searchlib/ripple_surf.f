C+------------------------------------------------------------------------------
C
      SUBROUTINE RIPPLE_SURF (MDIM, NDIM, M1, M2, N1, N2, X, Y, Z,
     >                        XP, YP, ZP, M, N, DTOL, DWEIGHT, P, Q,
     >                        DIST, ERROR, STATUS)
C
C     ONE-LINER:  "RIPPLE" search of a structured SURFace grid
C                  ------                         ----
C     PURPOSE:
C
C        This is an adaptation of RIPPLE2D.  RIPPLE2D is typically used to find
C     a given (u,v) in a parameterization of a 3-space surface grid (normalized
C     arc lengths in the two index directions).  Then the corresponding (x,y,z)
C     on the surface can be evaluated via parametric bilinear interpolation
C     within the cell that the search identified, using its coefficients (p,q).
C
C        RIPPLE_SURF effectively solves the inverse problem:  for a given target
C     (x,y,z) it determines the surface cell and coefficients (p,q) that most
C     nearly match (x,y,z) via bilinear interpolation within that cell.  Within
C     a given cell, this is the nonlinear least squares problem solved by the
C     utility PROJECT4.  RIPPLE_SURF serves to drive PROJECT4 in a reasonably
C     efficient search pattern, just as RIPPLE2D drives the underlying BILINT.
C
C        Since PROJECT4 determines the distance to the "plane" of the cell,
C     which could be zero even when the target point is well outside the cell,
C     a composite measure of goodness is required.  This is now a weighted sum
C     of the combined deviation of p and q from the unit interval and the
C     projected distance.  Actually, to balance these properly, the p and q
C     contributions are scaled by the averaged edge lengths of the cell.
C
C        Thus, if a (p,q) is found inside the unit square (to within the small
C     tolerance of PROJECT4), and the distance tolerance is satisfied, then the
C     desired cell has been located.  Otherwise, the cell with the lowest value
C     of
C               (dp * edgei + dq * edgej) * wpq + dist * (1 - wpq)
C
C     is returned as the best estimate.  This could still be large, meaning some
C     other surface patch really contains the target point, or the point is
C     outside the surface altogether.
C
C        Functions associated with the target (xp,yp,zp) may be interpolated as
C     follows, where (i,j) = the output (M,N):
C
C                f(p,q) = (1-q) ((1-p)f(i,j)   + pf(i+1,j)) +
C                            q  ((1-p)f(i,j+1) + pf(i+1,j+1))
C
C     METHOD (as for RIPPLE2D):
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
C     X,   (MDIM,NDIM) R   I         Surface grid coordinates
C     Y,
C     Z
C     XP,      -       R   I         Coordinates of the target point
C     YP,
C     ZP
C     M,       -       I   I O       Input indices of starting guess and
C     N                              output indices of "lower left corner"
C                                    of the enclosing cell (or nearest to it)
C     DTOL     -       R   I         Distance tolerance needed to prevent
C                                    accepting a (p,q) in some far away cell.
C     DWEIGHT  -       R   I         Weight to be applied to the DIST portion
C                                    of the closeness measure;  1 - DWEIGHT
C                                    is applied to the combined deviation of
C                                    (denormalized) p & q from being inside a
C                                    cell; try DWEIGHT = 0.5
C     P,       -       R     O       Bilinear interpolation coefficients;
C     Q                              see above for their usage
C     DIST     -       R     O       The distance of the target point to the
C                                    point in the ["plane" of] cell (M,N)
C                                    represented by the output (M,N) and (P,Q);
C                                    DIST >= 0.
C     ERROR    -       R     O       The best composite measure of goodness as
C                                    explained above;  ERROR >= 0.
C     STATUS   -       I     O       0 = match found; P & Q are both in [0,1];
C                                    1 = failed to find an enclosing cell
C                                        with (P,Q) in the unit square;
C                                        M, N, P, Q, and DSQ represent the
C                                        nearest results found; this could
C                                        mean some other surface patch may be
C                                        preferable, or the target is off the
C                                        edge of the data; since the search
C                                        will have checked every cell before
C                                        quitting, this case can become
C                                        expensive.  Adding an artificial outer
C                                        layer to the surface being searched may
C                                        help (suggestion by James Reuther).
C     PROCEDURES:
C
C        PROJECT4   Performs test for enclosing cell and returns best P, Q.
C
C     HISTORY:
C
C        11/10/93  M.Wong      Original design and coding of RIPPLE2D.
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
C        04/08/04  D.Saunders  Initial adaptation of RIPPLE2D as RIPPLE_SURF.
C        04/12/04      "       Added DTOL argument and included (DSQ - DTOL**2)
C                              in the measure of EBEST.
C        04/15/04   DAS/JJR    James Reuther suggesting denormalizing the p & q
C                              deviations to balance them with the distance
C                              deviation.  Go with unsquared distances to do
C                              this properly.  Also:  weight the two types of
C                              contribution, and return the best composite
C                              measure of closeness to save recovering it in
C                              the calling program.  Also: skip RIPPLE2D's
C                              possible adjustment of p & q to force p = 1
C                              to be p = 0 in the next cell, because this
C                              ignores the projected distance issue.
C     AUTHORS:
C
C        Michael D. Wong, Sterling Software, NASA/Ames, Moffett Field, CA
C        David Saunders,  ELORET/NASA Ames.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN)                      :: MDIM, NDIM,
     >                                             M1, M2, N1, N2

      INTEGER, INTENT (INOUT)                   :: M, N

      INTEGER, INTENT (OUT)                     :: STATUS

      REAL, INTENT (IN), DIMENSION (MDIM, NDIM) :: X, Y, Z

      REAL, INTENT (IN)                         :: XP, YP, ZP

      REAL, INTENT (IN)                         :: DTOL, DWEIGHT

      REAL, INTENT (OUT)                        :: P, Q, DIST, ERROR

C     Local constants:

      REAL, PARAMETER    :: HALF = 0.5E+0, ONE = 1.E+0, ZERO = 0.E+0
      LOGICAL, PARAMETER :: FALSE = .FALSE., TRUE = .TRUE.

C     Local variables:

      INTEGER :: I, J, L, IP, JP,
     >           ICORNER, M1P1, M2M1, M2M2, N1P1, N2M1, N2M2,
     >           NEDGE, ISTART, IEND, JEND, JSTART, MBEST, NBEST
      REAL    :: DSQ, DBEST, EBEST, EDGE1, EDGE2, ERR,
     >           PBEST, PERR, QBEST, QERR, PQWEIGHT
      REAL    :: X0(3), X1(3), X2(3), X3(3), X4(3), XTARG(3)
      LOGICAL :: ITEST, JTEST, TEST

C     Execution:

C     Begin with the initial guess.

      XTARG(1) = XP;  XTARG(2) = YP;  XTARG(3) = ZP;  I = M;  J = N

      CALL SETUP  ! Internal procedure sets up for projecting to current cell

      CALL PROJECT4 (X1, X2, X3, X4, XTARG, X0, DSQ, P, Q, STATUS)

      DIST = SQRT (DSQ)

      IF (STATUS == 0 .AND. DIST < DTOL) THEN ! Avoid denormalizing (p,q)
         ERROR = DWEIGHT * DIST
         GO TO 99  ! Success
      END IF 

      PQWEIGHT = ONE - DWEIGHT
      EBEST    = 1.E+32

      CALL MEASURE  ! Compute a composite measure of closeness; save the best

      IF (STATUS == 1) THEN ! Is PROJECT4's tolerance tighter than reqd. here?
         IF (EBEST < DTOL) THEN
            STATUS = 0
            GO TO 90
         END IF
      END IF

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

            ITEST = TRUE;  JTEST = TRUE

            SELECT CASE (ICORNER)

            CASE (1)

               I = M + L;  J = N + L

               IF (I > M2M1) THEN
                   I = M2M1;  ITEST = FALSE
               END IF

               IF (J > N2M1) THEN
                   J = N2M1;  JTEST = FALSE
               END IF

            CASE (2)

               I = M + L;  J = N - L

               IF (I > M2M1) THEN
                   I = M2M1;  ITEST = FALSE
               END IF

               IF (J < N1) THEN
                   J = N1;  JTEST = FALSE
               END IF

            CASE (3)

               I = M - L;  J = N - L

               IF (I < M1) THEN
                   I = M1;  ITEST = FALSE
               END IF

               IF (J < N1) THEN
                   J = N1;  JTEST = FALSE
               END IF

            CASE (4)

               I = M - L;  J = N + L

               IF (I < M1) THEN
                   I = M1;  ITEST = FALSE
               END IF

               IF (J > N2M1) THEN
                   J = N2M1;  JTEST = FALSE
               END IF

            END SELECT

C           Avoid repeated tests.

            IF (ITEST .OR. JTEST) THEN

               CALL SETUP

               CALL PROJECT4 (X1, X2, X3, X4, XTARG, X0, DSQ, P, Q,
     >                        STATUS)

               DIST = SQRT (DSQ)

               IF (STATUS == 0 .AND. DIST < DTOL) THEN
                  M = I
                  N = J
                  ERROR = DWEIGHT * DIST
                  GO TO 99  ! Success
               END IF

               CALL MEASURE

               IF (STATUS == 1) THEN
                  IF (EBEST < DTOL) THEN
                     STATUS = 0
                     GO TO 90
                  END IF
               END IF

            END IF

         END DO ! Next corner

C        Cycle through interior points of the connector lines one time only.

         DO NEDGE = 1, 4

            SELECT CASE (NEDGE)

            CASE (1)
               ISTART = M + L
               IEND   = ISTART
               JSTART = N - L + 1
               JEND   = N + L - 1
               TEST   = IEND <= M2M1
            CASE (2)
               ISTART = M - L + 1
               IEND   = M + L - 1
               JSTART = N - L
               JEND   = JSTART
               TEST   = JSTART >= N1
            CASE (3)
               ISTART = M - L
               IEND   = ISTART
               JSTART = N - L + 1
               JEND   = N + L - 1
               TEST   = ISTART >= M1
            CASE (4)
               ISTART = M - L + 1
               IEND   = M + L - 1
               JSTART = N + L
               JEND   = JSTART
               TEST   = JEND <= N2M1
            END SELECT

            IF (TEST) THEN

C              Limit the search to the interior points of the edge.

               IF (ISTART /= IEND) THEN
                  IF (ISTART < M1P1) ISTART = M1P1
                  IF (IEND   > M2M2) IEND   = M2M2
               ELSE IF (JSTART /= JEND) THEN
                  IF (JSTART < N1P1) JSTART = N1P1
                  IF (JEND   > N2M2) JEND   = N2M2
               END IF

               DO I = ISTART, IEND

                  DO J = JSTART, JEND

                     CALL SETUP

                     CALL PROJECT4 (X1, X2, X3, X4, XTARG, X0, DSQ,
     >                              P, Q, STATUS)

                     DIST = SQRT (DSQ)

                     IF (STATUS == 0 .AND. DIST < DTOL) THEN
                        M = I
                        N = J
                        ERROR = DWEIGHT * DIST
                        GO TO 99  ! Success
                     END IF

                     CALL MEASURE

                     IF (STATUS == 1) THEN
                        IF (EBEST < DTOL) THEN
                           STATUS = 0
                           GO TO 90
                        END IF
                     END IF

                  END DO

               END DO

            END IF

         END DO ! Next edge

C        Check to see if the entire grid has been tested.

         IF (M - L > M1 .OR. M + L < M2M1 .OR.
     >       N - L > N1 .OR. N + L < N2M1)
     >GO TO 10  ! Next search level


   90 CONTINUE

C     Return the best result - it may be adequate near a boundary:

      M     = MBEST
      N     = NBEST
      P     = PBEST
      Q     = QBEST
      DIST  = DBEST
      ERROR = EBEST

   99 CONTINUE


C     Internal procedures for RIPPLE_SURF:
C     ------------------------------------

      CONTAINS

C        ---------------------------------------------------
         SUBROUTINE SETUP ! Set up a quad. cell for PROJECT4
C        ---------------------------------------------------

C        PROJECT4 expects the cell to be specified in counterclockwise
C        order from the lower-left vertex.

         IP = I + 1;  JP = J + 1
         X1(1) = X(I,J);   X1(2) = Y(I,J);   X1(3) = Z(I,J)
         X2(1) = X(IP,J);  X2(2) = Y(IP,J);  X2(3) = Z(IP,J)
         X3(1) = X(IP,JP); X3(2) = Y(IP,JP); X3(3) = Z(IP,JP)
         X4(1) = X(I,JP);  X4(2) = Y(I,JP);  X4(3) = Z(I,JP)

         END SUBROUTINE SETUP

C        --------------------------------------------------------------
         SUBROUTINE MEASURE ! Evaluate a composite measure of closeness
C        --------------------------------------------------------------

         IF (P < ZERO) THEN
            PERR = -P
         ELSE IF  (P > ONE) THEN
            PERR = P - ONE
         ELSE
            PERR = ZERO
         END IF

         IF (PERR > ZERO) THEN
            EDGE1 = (X1(1)-X2(1))**2 +(X1(2)-X2(2))**2 +(X1(3)-X2(3))**2
            EDGE2 = (X3(1)-X4(1))**2 +(X3(2)-X4(2))**2 +(X3(3)-X4(3))**2
            PERR  = PERR * HALF * (SQRT (EDGE1) + SQRT (EDGE2))
         END IF

         IF (Q < ZERO) THEN
            QERR = -Q
         ELSE IF  (Q > ONE) THEN
            QERR = Q - ONE
         ELSE
            QERR = ZERO
         END IF

         IF (QERR > ZERO) THEN
            EDGE1 = (X1(1)-X4(1))**2 +(X1(2)-X4(2))**2 +(X1(3)-X4(3))**2
            EDGE2 = (X2(1)-X3(1))**2 +(X2(2)-X3(2))**2 +(X2(3)-X3(3))**2
            QERR  = QERR * HALF * (SQRT (EDGE1) + SQRT (EDGE2))
         END IF

         ERR = PQWEIGHT * (PERR + QERR) + DWEIGHT * DIST

         IF (ERR < EBEST) THEN
            EBEST = ERR
            DBEST = DIST
            PBEST = P
            QBEST = Q
            MBEST = I
            NBEST = J
         END IF

         END SUBROUTINE MEASURE

      END SUBROUTINE RIPPLE_SURF
