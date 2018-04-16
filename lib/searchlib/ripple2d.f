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
