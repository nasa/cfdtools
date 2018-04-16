C+----------------------------------------------------------------------
C
      SUBROUTINE DISTPTOS (P0, P1, P2, DIST, T, IER)
C
C     One-liner:  Computes the distance from a point to a line segment.
C     ----------
C
C     Description and usage:
C     ----------------------
C
C     DISTPTOS carefully computes the shortest distance from a point in
C     3-space to a (finite) line segment. There are three cases:
C
C        (1) The perpendicular from P0 to the segment P1, P2 intersects
C            the segment, i.e., the parametric distance from P1 along the
C            segment to the intersection point is >= 0 and <= dist (P1, P2).
C            DIST is just the length of the perpendicular, and T is the
C            corresponding parametric distance.
C
C        (2) The intersection's parametric distance along the segment is
C            < 0 (we have run off one end of the segment). DIST is the
C            Euclidean distance from P0 to P1, the nearest point on the
C            segment, and T is returned as calculated.
C
C        (3) We have run off the other end. DIST = dist (P0, P2), T left
C            as calculated (> dist (P1, P2) in this case).
C
C     An error flag is provided to indicate the degenerate case of P1 = P2,
C     which may or may not actually be a problem for the calling routine.
C     DIST is set to dist (P0, P1), and T is zero. Special cases (2)
C     and (3), above, are also indicated by IER. (Case (3) in particular
C     is easier to check at this level, if it is of interest.)
C
C     DISTPTOS and companion DISTPTOP were originally written for use by
C     MODGRID, which is used to prepare input for a finite element grid
C     generator where local cell size depends on the distance from certain
C     special line segments.
C
C     Arguments:
C     ----------
C
C     Name     Dimension  Type  I/O/S  Description
C     P0        3          R    I      Field point.
C
C     P1        3          R    I      First of two points which are the
C                                      endpoints of the reference segment.
C
C     P2        3          R    I      Second reference point.
C
C     DIST                 R      O    Distance from point P0 to the line
C                                      segment P1, P2.
C
C     T                    R      O    Unnormalized parametric distance
C                                      from P1 along the (extended) line
C                                      containing the segment to the point
C                                      closest to P0. Positive T is in the
C                                      direction from P1 to P2.
C
C     IER                  I      O    Error flag.
C                                         0: No problems, case (1) above;
C                                        +1: P1 is identical to P2; DIST
C                                            is distance from P0 to P1,
C                                            and T is set to zero.
C                                        +2: Intersection point is off
C                                            the P1 end (T <= zero);
C                                        +3: Intersection point is off
C                                            P2 end (T >= dist (P1, P2)).
C                                      It is not necessary to check IER
C                                      if the calling routine only needs
C                                      the distance to the segment and is
C                                      not affected by the special cases.
C
C     External modules:
C     -----------------
C
C     DISTPTOP  Computes the distance between two points in space.
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN 77 v4.7
C     ------------
C
C     Notes:
C     ------
C
C     (1)  IMPLICIT NONE and eight character symbols are not (yet) standard.
C
C     (2)  The formula for parametric distance from P1 to the perpendicular
C          intersection point was derived by differentiating the expression
C          for the distance from point P0 to an arbitrary point on the
C          line containing segment P1, P2 (in parametric form). Setting the
C          derivative equal to zero and solving for T, using the fact that
C          L ** 2 + M ** 2 + N ** 2 = 1.0, gives the result below.
C
C     Author:  Robert Kennelly, Sterling Federal Systems/NASA-Ames
C     -------
C
C     Development history:
C     --------------------
C
C     18 Jan. 1989    RAK    Initial design and coding.
C
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Constants.

      REAL
     &   ZERO, ONE
      PARAMETER
     &  (ZERO = 0.0E+0,
     &   ONE  = 1.0E+0)

C     Arguments.

      INTEGER
     &   IER
      REAL
     &   DIST, P0 (3), P1 (3), P2 (3), T

C     Local variables.

      REAL
     &   INVPP, L, M, N, PINT (3), P1TOP2

C     Execution.
C     ----------

      IER = 0
      CALL DISTPTOP (P1, P2, P1TOP2)

      IF (P1TOP2 .EQ. ZERO) THEN

C        The line segment is degenerate. Compute what we can and quit.

         IER = +1
         T   = ZERO
         CALL DISTPTOP (P0, P1, DIST)
         GO TO 990
      END IF

C     Form the (normalized) direction cosines for the line segment.

      INVPP = ONE / P1TOP2

      L = (P2 (1) - P1 (1)) * INVPP
      M = (P2 (2) - P1 (2)) * INVPP
      N = (P2 (3) - P1 (3)) * INVPP

C     Compute the parametric distance from P1 to the intersection point
C     of P0 with the line through P1, P2. See Note (2) above.

      T = L * (P0 (1) - P1 (1)) +
     &    M * (P0 (2) - P1 (2)) +
     &    N * (P0 (3) - P1 (3))

C     Off the end? If so, compute distance to nearest endpoint and quit.
C     Note that in these cases T doesn't "match" the point on the segment
C     from which the distance was calculated.

      IF (T .LT. ZERO) THEN
         IER = +2
         CALL DISTPTOP (P0, P1, DIST)
      ELSE IF (T .GT. P1TOP2) THEN
         IER = +3
         CALL DISTPTOP (P0, P2, DIST)
      ELSE

C        Plug in T to get the intersection point, which must actually lie
C        on the finite-length segment in this case, then compute distance
C        to P0 and we're done.

         PINT (1) = P1 (1) + T * L
         PINT (2) = P1 (2) + T * M
         PINT (3) = P1 (3) + T * N
         CALL DISTPTOP (P0, PINT, DIST)
      END IF

C     Termination.
C     ------------

  990 CONTINUE
      RETURN
      END

