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
