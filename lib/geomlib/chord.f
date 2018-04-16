C+----------------------------------------------------------------------
C
      FUNCTION CHORD (X, Y, I1, I2)
C
C     One-liner: Summed chord-lengths for X-Y curve over range of indices
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        Computes the sum of the Euclidean distance between adjacent
C     points in a curve represented by two arrays. The calculation is
C     rearranged so as to avoid (almost) all chance of overflow, and
C     any unnecessary loss of precision when one component of the
C     distance from one point to the next is small relative to the other.
C     The calling routine must supply beginning and ending indices for
C     the summation (this is intended to facilitate operations with
C     packed data arrays). The result does not depend on the order of
C     I1 and I2.
C
C        CHORD was originally written for use with PLSFIT, which performs
C     parametric cubic interpolation with cumulative chord length as the
C     curve parameter. In use, it is a good idea to try to use all
C     available information to avoid (expensively) re-calculating the
C     lengths of the same intervals over and over; CHORD should be
C     thought of as providing length increments.
C
C     Arguments:
C     ----------
C
C     Name    Type/Dimension  I/O/S  Description
C     X         R (*)         I      Array of abscissas.
C
C     Y         R (*)         I      Array of ordinates.
C
C     I1,I2                   I      Indices for summation. The loop
C                                    runs from MIN(I1,I2) to MAX(I1,I2)
C                                    so the result is independent of
C                                    order.
C
C     CHORD   R                 O    Function value is the sum of the
C                                    chord lengths along the curve
C                                    defined by the X and Y arrays
C                                    between indices I1 and I2.
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
C      8 Apr. 1988  RAK/DAS  Reformulated to reduce chance of overflow
C                            or unnecessary underflow. Result is not
C                            dependent on order of I1, I2.
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
     &   I1, I2
      REAL
     &   CHORD, X (*), Y (*)

C     Local variables.

      INTEGER
     &   I
      REAL
     &   DENOM, DX, DY

C     Execution.
C     ----------

      CHORD = ZERO

      DO 10, I = MIN (I1, I2), MAX (I1, I2) - 1
         DX    = ABS (X (I + 1) - X (I))
         DY    = ABS (Y (I + 1) - Y (I))
         DENOM = MAX (DX, DY)
         IF (DENOM .GT. ZERO) CHORD = CHORD +
     &      DENOM * SQRT (ONE + (MIN (DX, DY) / DENOM) ** 2)
   10 CONTINUE

C     Termination.
C     ------------

      RETURN
      END

