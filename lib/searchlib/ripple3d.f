C+------------------------------------------------------------------------------
C
      SUBROUTINE RIPPLE3D (IDIM, JDIM, KDIM, I1, I2, J1, J2, K1, K2,
     >                     X, Y, Z, XP, YP, ZP, ICELL, JCELL, KCELL,
     >                     EPS, P, Q, R, STATUS)
C
C     ONE-LINER:  "RIPPLE" search of a 3D quasi-rectangular grid
C                  ------              --
C     PURPOSE:
C
C        RIPPLE3D performs a search within a quasi-rectangular 3D grid for
C     the cell containing a target point by using a "ripple" pattern.  It also
C     returns the corresponding coefficients needed for trilinear interpolation.
C     (See TRILINT.)  The search should be efficient for typical applications
C     where the specified initial cell is a good starting guess (probably the
C     result of a previous call).
C
C     METHOD:
C
C        Points are tested around the starting guess in rectangular "layers"
C     of increasing size.  If the cell defined by (ICELL, JCELL, KCELL) does
C     not contain the target point, then the next search domain would be the
C     surrounding grid layer, then the layer surrounding that, and so on.
C
C        Subroutine TRILINT is used to perform the test for containment. The
C     resulting coefficients P, Q, R, may be used at the higher level if
C     trilinear interpolation is adequate.  This implementation is simpler
C     than that of RIPPLE2D, which became more "involved" than expected
C     in order to avoid some repetition once boundaries are encountered.
C     Algorithms involving more clever searches which step towards the
C     target cell in 2- or 3-space are certainly more efficient yet.
C     The present approach has the virtue of simplicity.
C
C     ARGUMENTS:
C
C     NAME    DIM    TYPE  I/O/S     DESCRIPTION
C     ----    ---    ----  -----     -----------
C     IDIM,    -       I   I         Dimensions of the grid arrays in the
C     JDIM,                          calling program
C     KDIM 
C     I1,I2,   -       I   I         The (sub)grid to be searched is
C     J1,J2,                         indicated by (I1:I2, J1:J2, K1:K2)
C     K1,K2
C     X,    (IDIM,     R   I         Grid coordinates
C     Y,     JDIM,
C     Z      KDIM)
C     XP,      -       R   I         Coordinates of the target point
C     YP,
C     ZP
C     ICELL,   -       I   I O       Input indices of starting guess and
C     JCELL,                         output indices of "lower left corner"
C     KCELL                          of the enclosing cell (or nearest to it)
C     EPS      -       R   I         Search tolerance used by TRILINT, q.v.
C     P,       -       R     O       Trilinear interpolation coefficients
C     Q,                             to be used as described in TRILINT
C     R
C     STATUS   -       I     O       0 = match found;
C                                    1 = failed to find match; I/J/KCELL and
C                                        P, Q, R represent the smallest mismatch
C                                    (Note that for the sake of efficiency, the
C                                    matrix singularity condition checked for
C                                    by TRILINT is not tested for here.)
C
C     PROCEDURES:
C
C        TRILINT     Performs test for enclosing cell
C
C     ENVIRONMENT:
C
C        FORTRAN 77 with minor extensions
C
C     HISTORY:
C
C        08/08/94   DAS   Initial implementation derived from J. Reuther's
C                         programming in the OPT67 wing/body design code, but
C                         in reusable form analogous to M. Wong's RIPPLE2D.
C        01/23/96   DAS   Returned the nearest results to handle possible
C                         difficulties near the boundaries of (sub)regions.
C                         This may help establish an appropriate EPS too.
C        01/10/98   DAS   Call TRILINT again after adjusting for P/Q/R >= 1.
C
C     AUTHOR:
C
C        David Saunders, Sterling Software/NASA Ames, Moffett Field, CA
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   IDIM, JDIM, KDIM, I1, I2, J1, J2, K1, K2, ICELL, JCELL, KCELL,
     >   STATUS
      REAL
     >   X (IDIM, JDIM, KDIM), Y (IDIM, JDIM, KDIM),
     >   Z (IDIM, JDIM, KDIM), XP, YP, ZP, EPS, P, Q, R

C     Local constants.

      REAL
     >   ONE, ZERO
      PARAMETER
     >  (ONE  = 1.E+0,
     >   ZERO = 0.E+0)

C     Local variables.

      INTEGER
     >   I, J, K, L, LMAX, INDEX1, INDEX2, IBEST, JBEST, KBEST
      REAL
     >   EBEST, ERR, PBEST, PERR, QBEST, QERR, RBEST, RERR

C     Execution.


C     Begin with the initial guess.

      CALL TRILINT (IDIM, JDIM, KDIM, X, Y, Z, XP, YP, ZP, ICELL, JCELL,
     >              KCELL, EPS, P, Q, R, STATUS)

      IF (STATUS .EQ. 0) GO TO 90  ! For possible special-case adjustment

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

      IF (R .LT. ZERO) THEN
         RERR = -R
      ELSE IF (R .GT. ONE) THEN
         RERR = R - ONE
      ELSE
         RERR = ZERO
      END IF

      EBEST = PERR + QERR + RERR
      PBEST = P
      QBEST = Q
      RBEST = R
      IBEST = ICELL
      JBEST = JCELL
      KBEST = KCELL

C     Enter the rectangular layer search pattern.
C     The order K, J, I for the sublayers searched is arbitrary.

      LMAX = MAX (ICELL - I1, I2 - ICELL - 1,
     >            JCELL - J1, J2 - JCELL - 1,
     >            KCELL - K1, K2 - KCELL - 1)

      L = 1       ! Layer #
   10 CONTINUE

C        Two full K sub-planes of size 2L + 1 by 2L + 1 (except at an
C        I or J boundary).  We really need a case statement.
C        Note: Once a K boundary is encountered, the outer increment is
C        no longer 2L, and the next layer repeats all of this part of
C        the boundary.  For simplicity, we do not attempt to avoid that.

         INDEX1 = MAX (KCELL - L, K1)
         INDEX2 = MIN (KCELL + L, K2 - 1)

         DO K = INDEX1, INDEX2, INDEX2 - INDEX1
            DO J = MAX (J1, JCELL - L), MIN (J2 - 1, JCELL + L)
               DO I = MAX (I1, ICELL - L), MIN (I2 - 1, ICELL + L)

                  CALL TRILINT (IDIM, JDIM, KDIM, X, Y, Z, XP, YP,
     >                          ZP, I, J, K, EPS, P, Q, R, STATUS)

                  IF (STATUS .EQ. 0) THEN
                     ICELL = I
                     JCELL = J
                     KCELL = K
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

                  IF (R .LT. ZERO) THEN
                     RERR = -R
                  ELSE IF (R .GT. ONE) THEN
                     RERR = R - ONE
                  ELSE
                     RERR = ZERO
                  END IF

                  ERR = PERR + QERR + RERR
                  IF (ERR .LT. EBEST) THEN
                     EBEST = ERR
                     PBEST = P
                     QBEST = Q
                     RBEST = R
                     IBEST = I
                     JBEST = J
                     KBEST = K
                  END IF
               END DO
            END DO
         END DO

C        Two J sub-planes of size 2L - 1 by 2L + 1:

         INDEX1 = MAX (JCELL - L, J1)
         INDEX2 = MIN (JCELL + L, J2 - 1)

         DO J = INDEX1, INDEX2, INDEX2 - INDEX1
            DO K = MAX (K1, KCELL - L + 1), MIN (K2 - 1, KCELL + L - 1)
               DO I = MAX (I1, ICELL - L), MIN (I2 - 1, ICELL + L)

                  CALL TRILINT (IDIM, JDIM, KDIM, X, Y, Z, XP, YP,
     >                          ZP, I, J, K, EPS, P, Q, R, STATUS)

                  IF (STATUS .EQ. 0) THEN
                     ICELL = I
                     JCELL = J
                     KCELL = K
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

                  IF (R .LT. ZERO) THEN
                     RERR = -R
                  ELSE IF (R .GT. ONE) THEN
                     RERR = R - ONE
                  ELSE
                     RERR = ZERO
                  END IF

                  ERR = PERR + QERR + RERR
                  IF (ERR .LT. EBEST) THEN
                     EBEST = ERR
                     PBEST = P
                     QBEST = Q
                     RBEST = R
                     IBEST = I
                     JBEST = J
                     KBEST = K
                  END IF
               END DO
            END DO
         END DO

C        Two I sub-planes of size 2L - 1 by 2L - 1:

         INDEX1 = MAX (ICELL - L, I1)
         INDEX2 = MIN (ICELL + L, I2 - 1)

         DO I = INDEX1, INDEX2, INDEX2 - INDEX1
            DO J = MAX (J1, JCELL - L + 1), MIN (J2 - 1, JCELL + L - 1)
               DO K = MAX (K1, KCELL - L + 1), MIN (K2 - 1, KCELL + L-1)

                  CALL TRILINT (IDIM, JDIM, KDIM, X, Y, Z, XP, YP,
     >                          ZP, I, J, K, EPS, P, Q, R, STATUS)

                  IF (STATUS .EQ. 0) THEN
                     ICELL = I
                     JCELL = J
                     KCELL = K
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

                  IF (R .LT. ZERO) THEN
                     RERR = -R
                  ELSE IF (R .GT. ONE) THEN
                     RERR = R - ONE
                  ELSE
                     RERR = ZERO
                  END IF

                  ERR = PERR + QERR + RERR
                  IF (ERR .LT. EBEST) THEN
                     EBEST = ERR
                     PBEST = P
                     QBEST = Q
                     RBEST = R
                     IBEST = I
                     JBEST = J
                     KBEST = K
                  END IF
               END DO
            END DO
         END DO

         L = L + 1  ! Increment "ripple" level.

         IF (L .LE. LMAX)
     >GO TO 10


C     Dropped through - entire (sub)grid was tested without success.
C     Return the nearest result - may be adequate near a boundary:

      ICELL = IBEST
      JCELL = JBEST
      KCELL = KBEST
      P     = PBEST
      Q     = QBEST
      R     = RBEST


   90 CONTINUE

C     Ensure "lower left" indices as for INTERVAL's 1-D search:

      IF (P .GE. ONE) THEN
         IF (ICELL + 1 .LT. I2) THEN
            ICELL = ICELL + 1
            CALL TRILINT (IDIM, JDIM, KDIM, X, Y, Z, XP, YP, ZP,
     >                    ICELL, JCELL, KCELL, EPS, P, Q, R, STATUS)
         END IF
      END IF

      IF (Q .GE. ONE) THEN
         IF (JCELL + 1 .LT. J2) THEN
            JCELL = JCELL + 1
            CALL TRILINT (IDIM, JDIM, KDIM, X, Y, Z, XP, YP, ZP,
     >                    ICELL, JCELL, KCELL, EPS, P, Q, R, STATUS)
         END IF
      END IF

      IF (R .GE. ONE) THEN
         IF (KCELL + 1 .LT. K2) THEN
            KCELL = KCELL + 1
            CALL TRILINT (IDIM, JDIM, KDIM, X, Y, Z, XP, YP, ZP,
     >                    ICELL, JCELL, KCELL, EPS, P, Q, R, STATUS)
         END IF
      END IF

      RETURN
      END
