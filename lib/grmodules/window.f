C+----------------------------------------------------------------------
C
      SUBROUTINE WINDOW (NDATA, X, Y, X1, X2, Y1, Y2, NBEGIN, NFINAL,
     >   IER)
C
C     One-liner:  Find a segment of curve (X, Y) enclosed by a rectangle
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        WINDOW is a simple method for finding the portions of an (X, Y)
C     curve which lie inside a rectangular window, for the purpose of
C     plotting parametric curves.  We wish to determine which points must
C     be connected in order to save time by not fitting the portions of
C     the curve which will not be visible.  The beginning index returned
C     is thus either the first one (if the initial point lies inside the
C     window) or indicates the point immediately preceding the first point
C     which does lie inside.  The final index is set analogously.  For a
C     segment to be identified as interior, at least one point (X, Y) must
C     lie inside the window.
C
C        Indices NBEGIN and NFINAL must be initialized by the calling
C     routine to delimit the range of interest.  WINDOW resets these on
C     return, and the user must update them before calling WINDOW again.
C     An error flag is set if the number of data points supplied is
C     less than 2, or if no portion of the curve is visible.
C
C     Arguments:
C     ----------
C
C     Name     Type/Dimension  I/O/S  Description
C     NDATA    I               I      Length of the X, Y data arrays.
C
C     X        R (NDATA)       I      Array of input points' first
C                                     components (abscissas).
C
C     Y        R (NDATA)       I      Second components (ordinates).
C
C     X1,X2    R               I      "Side" coordinates of the rectangular
C                                     window of interest.  The values need
C                                     not be ordered - each point's "X"
C                                     coordinate will be checked to see if
C                                     it is between X1 and X2.
C
C     Y1,Y2    R               I      "Top" and "bottom" coordinates of the
C                                     window.
C
C     NBEGIN   I               I/O    Set by caller to index of first
C                                     point to be examined (usually
C                                     just 1).  On output, NBEGIN is
C                                     the index of the first point before
C                                     the window, i.e.,  point NBEGIN + 1
C                                     is inside.  NBEGIN = 1 if the
C                                     curve originates inside.
C
C     NFINAL   I               I/O    Set by caller to index of last
C                                     point to be examined (usually
C                                     just NDATA).  On output, NFINAL is
C                                     the index of the first point beyond
C                                     the window, i.e.,  point NFINAL - 1
C                                     is inside.  NFINAL = NDATA if the
C                                     curve terminates inside.
C
C     IER      I                 O    Error flag:
C                                       0 = normal return (a non-empty
C                                           segment was found);
C                                       1 = NDATA < 2, or NFINAL - NBEGIN
C                                           < 1 (either a blunder, or
C                                           we've run out of data);
C                                       2 = no segment lies within the
C                                           window.
C
C     Significant local variables:
C     ----------------------------
C
C     CURROUT  Status of current point, .TRUE. means outside the window
C     PREVOUT  Status of previous point
C     EDGES    Counts edge-crossings
C     OUTSIDE  Statement function which returns status of a point
C
C     Environment:  Digital VAX-11/780 VMS FORTRAN 77
C     ------------  Apple Macintosh Absoft MacFORTRAN/020
C
C     Notes:
C     ------
C
C     (1)  IMPLICIT NONE, 8-character symbolic names, and use of "!" for
C          comments are non-standard (until FORTRAN-8X).
C
C     (2)  For the time being, we do NOT attempt to handle the case
C          where two adjacent points in the data arrays are both outside
C          the window, but where the fitted curve loops inside.  We also
C          ignore the possiblity that a portion of the fitted curve
C          loops outside and then back - such a segment is considered
C          to be entirely enclosed by the window.  The rationale is that
C          the cost of examining such cases outweighs their benefit,
C          since the first case is presumably very rare, and the second
C          merely means that the interior parts of the reentrant curve
C          segment are fitted with somewhat fewer points than they would
C          be if each portion were fitted separately.  Time will tell if
C          this is acceptable in practice - a next level of sophisti-
C          cation might be to implement a LINEAR clipping algorithm and
C          check for intersection of the successive chords with the
C          window.  (See, for example CLIPPER in Foley & van Dam.)
C
C     Author:  Robert Kennelly, NASA-Ames/Sterling Federal Systems
C     -------
C
C     Development history:
C     --------------------
C
C     12 Mar. 1987    RAK    Initial design and coding.
C     20 Apr. 1987    RAK    Generalized - X1,X2 and Y1,Y2 pairs need not
C                            be ordered.
C
C-----------------------------------------------------------------------

C     Declarations.
C     -------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   IER, NBEGIN, NDATA, NFINAL
      REAL
     >   X1, X2, X (NDATA), Y1, Y2, Y (NDATA)

C     Local variables.

      LOGICAL
     >   CURROUT, PREVOUT
      INTEGER
     >   EDGES, N
      REAL
     >   XLEFT, XRIGHT, YBOTTOM, YTOP

C     Procedures.

      LOGICAL
     >   OUTSIDE
      REAL
     >   XDUMMY, YDUMMY

      OUTSIDE (XDUMMY, YDUMMY) = (XDUMMY .LT. XLEFT .OR.
     >                            XDUMMY .GT. XRIGHT .OR.
     >                            YDUMMY .LT. YBOTTOM .OR.
     >                            YDUMMY .GT. YTOP)

C     Execution.
C     ----------

      IER    = 0
      NBEGIN = MAX (1, NBEGIN)
      NFINAL = MIN (NDATA, NFINAL)

      IF (NDATA .LT. 2 .OR. (NFINAL - NBEGIN) .LT. 0) THEN

C        There is not enought data to draw a line through, or there are
C        no points to be checked.

         IER = +1
         GO TO 990
      END IF

C      ----------------------------------------------------------------
C     |                                                                |
C     |  1 <= NBEGIN <= NFINAL <= NDATA                                |
C     |  2 <= NDATA                                                    |
C     |                                                                |
C      ----------------------------------------------------------------

      XLEFT   = MIN (X1, X2)
      XRIGHT  = MAX (X1, X2)
      YBOTTOM = MIN (Y1, Y2)
      YTOP    = MAX (Y1, Y2)

C     Loop over the data (if there is more than one point), looking for
C     edge crossings.

      CURROUT = OUTSIDE (X (NBEGIN), Y (NBEGIN))
      EDGES   = 0

      DO 10, N = NBEGIN + 1, NFINAL

C        Save the status of the previous point.

         PREVOUT = CURROUT
         CURROUT = OUTSIDE (X (N), Y (N))
         IF (CURROUT .NEQV. PREVOUT) THEN

C           An edge has been crossed, do the bookkeeping.

            EDGES = EDGES + 1
            IF (CURROUT) THEN

C              It's the end of the line as soon as a transition is made
C              from inside to outside.  Set "right-hand" index and quit.

               NFINAL = MIN (NDATA, N + 1)
               GO TO 990
            ELSE

C              We just entered the visible region from outside.  Set
C              "left-hand" index and continue looping to look for an
C              exit transition.

               NBEGIN = N - 1
            END IF
         END IF
   10 CONTINUE

C     Clean up some special cases.

      IF (EDGES .EQ. 0 .AND. CURROUT) THEN

C        Evidently the entire curve lies outside the window.

         IER = +2
         GO TO 990
      ELSE IF (EDGES .EQ. 1 .AND. .NOT.CURROUT) THEN

C        The final point is inside.

         NFINAL = NDATA
      END IF

C     Termination.
C     ------------

  990 CONTINUE
      RETURN
      END
