C+----------------------------------------------------------------------
C
      SUBROUTINE FIXWINDO (N, X, Y, XMIN, XMAX, YMIN, YMAX, FLAG, EDGES)
C
C ONE-LINER:  FIX up a partially-defined plotting WINDOw
C
C PURPOSE:
C
C        FIXWINDO deals with the problem of determining a well-chosen
C     plotting window which is only partially defined to begin with.
C     Use of the global data limits is not smart enough: the window may
C     be smaller than this because only the data satisfying the partial
C     specification should be eligible, not all of it.
C
C METHOD:
C        Very little is done if the entire window is already specified.
C     Otherwise, the global data limits are determined.  If neither Y
C     limit is specified, then the X window cannot be modified, and
C     this case is handled for efficiency reasons: the X limits are set
C     to the global data extrema for X unless they are specified.
C     The case of neither X limit's being specified is handled analogously.
C
C        If either or both Y limits are specified, the X window is reduced
C     to the extremes of the data points within the Y window, and vice
C     versa.  These cases require more work.
C
C        The strategy is for LAYOUT (or equivalent) to use EDGES (1:2)
C     and EDGES (3:4) in its calls to LINAX instead of the full X (*)
C     and Y (*), thus leaving LINAX intact.  Modularizing the determin-
C     ation of EDGES (*) here leaves LAYOUT (already non-trivial)
C     largely unchanged.
C
C ERROR HANDLING:
C
C        If a specified window edge (in X, say) excludes all the data,
C     the edge(s) being determined (Ymin and/or Ymax) are assigned the
C     the appropriate global data value rather than values suggesting
C     a very narrow window.  This situation is assumed to be handled at
C     a higher level.
C
C ENVIRONMENT:
C     VAX/VMS; FORTRAN 77, with...
C     > IMPLICIT NONE
C     > Trailing ! comments
C     > Names up to 8 characters
C
C HISTORY:
C     10/06/91  D.A.Saunders   Initial implementation to solve a LAYOUT
C                              problem with partially-windowed data.
C
C AUTHOR: David Saunders, Sterling Software/NASA Ames, Mountain View, CA.
C
C ----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER   N         ! (I) Number of data points. N >= 1.
      REAL      X (N)     ! (I) Given abscissas ...
      REAL      Y (N)     ! (I) ... and ordinates.
      REAL      XMIN      ! (I) Partial (?) definition of plotting window
      REAL      XMAX      ! (I)    as input to (say) LAYOUT, q.v.
      REAL      YMIN      ! (I)
      REAL      YMAX      ! (I)
      REAL      FLAG      ! (I) Flag for indicating that XMIN, XMAX, etc.
                          !        is not yet specified.
      REAL      EDGES (4) ! (O) Xmin, Xmax, Ymin, Ymax determined as
                          !        defining a well-chosen plot window.

C-----------------------------------------------------------------------

C     Local variables.

      INTEGER   I, I1
      REAL      BOUNDHI, BOUNDLO, WMAX, WMIN, XHI, XLO, YHI, YLO
      LOGICAL   HARDER

C     Execution.

      EDGES (1) = XMIN
      EDGES (2) = XMAX
      EDGES (3) = YMIN
      EDGES (4) = YMAX

      HARDER = .FALSE.
      DO I = 1, 4
         IF (EDGES (I) .EQ. FLAG) HARDER = .TRUE.
      END DO

      IF (HARDER) THEN

C        Find the global data range, which may or may not be narrowed.

         XLO = X (1)
         XHI = XLO
         YLO = Y (1)
         YHI = YLO
         DO 200, I = 2, N
            IF (X (I) .GT. XHI) THEN   ! It can't also be a new minimum...
               XHI = X (I)
            ELSE IF (X (I) .LT. XLO) THEN
               XLO = X (I)
            END IF
            IF (Y (I) .GT. YHI) THEN
               YHI = Y (I)
            ELSE IF (Y (I) .LT. YLO) THEN
               YLO = Y (I)
            END IF
  200    CONTINUE

C        Consider the case of specifying YMIN:
C        XMIN and XMAX should allow for only those points I with
C        Y (I) >= YMIN, not for the full X range.
C        But if neither Y limit were specified, then global data limits
C        for X will do where an X limit was not specified.

         BOUNDLO = EDGES (3)    ! YMIN
         IF (BOUNDLO .EQ. FLAG) BOUNDLO = YLO
         BOUNDHI = EDGES (4)    ! YMAX
         IF (BOUNDHI .EQ. FLAG) BOUNDHI = YHI

         IF (BOUNDLO .GT. YLO .OR. BOUNDHI .LT. YHI) THEN

            ! The Y range was restricted.  The X range may be affected too.

            WMIN = FLAG            ! No valid values can be guaranteed
            WMAX = FLAG
            DO 300, I = 1, N
               I1 = I
               IF (BOUNDLO .LE. Y (I) .AND. Y (I) .LE. BOUNDHI) THEN
                  WMIN = X (I)
                  WMAX = WMIN
                  GO TO 310
               END IF
  300       CONTINUE
  310       DO 320, I = I1 + 1, N
               IF (BOUNDLO .LE. Y (I) .AND. Y (I) .LE. BOUNDHI) THEN
                  WMIN = MIN (WMIN, X (I))
                  WMAX = MAX (WMAX, X (I))
               END IF
  320       CONTINUE

C           Only assign reduced limits if they weren't specified.

            IF (EDGES (1) .EQ. FLAG) EDGES (1) = WMIN  ! May still be undefined.
            IF (EDGES (2) .EQ. FLAG) EDGES (2) = WMAX

         END IF


C        Analogous treatment of Y where X range may or may not be limited.

         BOUNDLO = EDGES (1)
         IF (BOUNDLO .EQ. FLAG) BOUNDLO = XLO
         BOUNDHI = EDGES (2)
         IF (BOUNDHI .EQ. FLAG) BOUNDHI = XHI

         IF (BOUNDLO .GT. XLO .OR. BOUNDHI .LT. XHI) THEN

            ! The X range was restricted.  The Y range may be affected too.

            WMIN = FLAG
            WMAX = FLAG

            DO 500, I = 1, N
               I1 = I
               IF (BOUNDLO .LE. X (I) .AND. X (I) .LE. BOUNDHI) THEN
                  WMIN = Y (I)
                  WMAX = WMIN
                  GO TO 510
               END IF
  500       CONTINUE
  510       DO 520, I = I1 + 1, N
               IF (BOUNDLO .LE. X (I) .AND. X (I) .LE. BOUNDHI) THEN
                  WMIN = MIN (WMIN, Y (I))
                  WMAX = MAX (WMAX, Y (I))
               END IF
  520       CONTINUE

            IF (EDGES (3) .EQ. FLAG) EDGES (3) = WMIN
            IF (EDGES (4) .EQ. FLAG) EDGES (4) = WMAX

         END IF

C        Finally, ensure that all limits are defined.

         IF (EDGES (1) .EQ. FLAG) EDGES (1) = XLO
         IF (EDGES (2) .EQ. FLAG) EDGES (2) = XHI
         IF (EDGES (3) .EQ. FLAG) EDGES (3) = YLO
         IF (EDGES (4) .EQ. FLAG) EDGES (4) = YHI

      END IF

      RETURN
      END
