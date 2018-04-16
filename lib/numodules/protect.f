C+----------------------------------------------------------------------
C
      SUBROUTINE PROTECT (NX, X, Y, ARROW, DISTINCT)
C
C     One-liner: Test data for monotonicity and distinctness
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        PROTECT scans a curve represented by arrays X and Y to check
C     strict monotonicity (increasing or decreasing) of the X data and
C     to verify that no two successive points match. One or the other
C     of these tests is required by conventional or parametric spline
C     fitting methods. Since it is often inefficient to perform these
C     tests in the fitting routine, which may be called more than once
C     with the same data or with data known to be good, the present
C     modular approach was chosen.
C
C        The initial application is program SMOOTH, where a single set
C     of data points may be fit by several different techniques - the
C     idea was to check each dataset just once as it was read in, and
C     then test flags ARROW and DISTINCT as needed.
C
C     Arguments:
C     ----------
C
C     Name     Type/Dimension  I/O/S  Description
C     NX       I               I      Dimension of X and Y arrays. If
C                                     NX <= 1, then ARROW will be
C                                     0.0 and DISTINCT will be .FALSE.
C
C     X        R (NX)          I      Array of abscissas.
C
C     Y        R (NX)          I      Array of ordinates.
C
C     ARROW    R                 O    Monotonicity indicator:
C                                       -1.0  strictly decreasing
C                                        0.0  neither
C                                       +1.0  strictly increasing
C
C     DISTINCT L                 O    Indicates whether successive points
C                                     are distinct in both X and Y.
C
C     Notes:
C     ------
C
C     (1)  IMPLICIT NONE and 8-character symbolic names are non-standard.
C
C     (2)  There is no provision for roundoff error in the monotonicity
C          test, i.e., all interval lengths are compared to zero.
C
C     Author:  Robert Kennelly, Sterling Federal Systems/NASA-Ames
C     -------
C
C     Development history:
C     --------------------
C
C     22 Feb. 1987    RAK    Initial design and coding.
C     14 Apr. 1988    RAK    Corrected value returned by DISTINCT when
C                            NX <= 1. Revised comments.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Constants.

      REAL
     &   ZERO, ONE
      PARAMETER
     &  (ZERO   = 0.0E+0,
     &   ONE    = 1.0E+0)

C     Arguments.

      LOGICAL
     &   DISTINCT
      INTEGER
     &   NX
      REAL
     &   ARROW, X (NX), Y (NX)

C     Local variables.

      LOGICAL
     &   FIRST
      INTEGER
     &   J
      REAL
     &   DX

C     Execution.
C     ----------

C     Set the output flags and bail out if the input data is trivial
C     or improper.

      ARROW    = ZERO
      DISTINCT = .FALSE.
      IF (NX .LE. 1) GO TO 990

C     Reset the direction flag according to the first interval, and reset
C     distinctness flag prior to testing.

      DX = X (2) - X (1)
      IF (DX .NE. ZERO) THEN
         FIRST = (DX .GT. ZERO)
         ARROW = SIGN (ONE, DX)
      END IF

      DISTINCT = .TRUE.

C     The approach is to try to set ARROW = ZERO and DISTINCT = .FALSE.
C     as we go, and quit as soon as possible. No harm is done if ARROW is
C     set repeatedly - it's not worth testing for.

      DO 10, J = 1, NX - 1
         DX = X (J + 1) - X (J)
         IF (DX .NE. ZERO) THEN

C           Compare the sign of the increment in this interval to that of
C           the first interval.

            IF ((DX .GT. ZERO) .NEQV. FIRST) ARROW = ZERO
         ELSE

C           If a pair of X's match, the data is not strictly monotonic. We
C           must still check whether the Y's are distinct.

            ARROW = ZERO
            IF (Y (J + 1) - Y (J) .EQ. ZERO) THEN

C              The data is neither monotonic nor distinct - time to die.

               DISTINCT = .FALSE.
               GO TO 990

            END IF
         END  IF
   10 CONTINUE

C     Termination.
C     ------------

  990 CONTINUE
      RETURN
      END
