C+------------------------------------------------------------------------------
C
      SUBROUTINE CHORDSRF (IDIM, JDIM, IROW, J1, J2, X, Y, Z, NORMALIZ,
     >                     TOTAL, CHORD)
C
C  ONE-LINER: [Relative] chord-lengths for the Ith [sub]row of a 3-space surface
C
C  DESCRIPTION:
C
C        CHORDSRF is a variant of CHORDS3D for treating a row of a surface grid.
C     (CHORDS3D handles either a line or a column of a grid.)
C
C  ARGUMENTS:
C
C     Name    Type/Dimension  I/O/S  Description
C
C     IDIM,     I             I      Surface grid dimensions in calling program
C     JDIM
C     IROW      I             I      Index of grid row to be processed
C     J1,       I             I      Portion of the row to be processed; J2 > J1
C     J2
C     X,        R (*)         I      Coordinates of data points
C     Y,
C     Z
C     NORMALIZ  L             I      .TRUE. means normalize results to [0,1]
C
C     TOTAL     R               O    Total chord length, returned because it is
C                                    otherwise lost if results are normalized
C     CHORD     R (*)           O    CHORD(J1:J2) contains cumulative
C                                    chord lengths, possibly normalized
C
C  HISTORY:
C
C     21 Sep. 1996   DAS   Adaptation of CHORDS3D.
C
C  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     &   IDIM, JDIM, IROW, J1, J2
      REAL
     &   X (IDIM, JDIM), Y (IDIM, JDIM), Z (IDIM, JDIM), CHORD (*),
     &   TOTAL
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
     &   I, J
      REAL
     &   DINV

C     Execution.

      I = IROW
      CHORD (J1) = ZERO

      DO J = J1 + 1, J2
         CHORD (J) = CHORD (J - 1) +
     &      SQRT ((X (I, J) - X (I, J - 1)) ** 2 +
     &            (Y (I, J) - Y (I, J - 1)) ** 2 +
     &            (Z (I, J) - Z (I, J - 1)) ** 2)
      END DO

      TOTAL = CHORD (J2)

      IF (NORMALIZ) THEN
         DINV = ONE / TOTAL
         DO J = J1 + 1, J2 - 1
            CHORD (J) = CHORD (J) * DINV
         END DO
         CHORD (J2) = ONE
      END IF

      RETURN
      END
