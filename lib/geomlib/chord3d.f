C+----------------------------------------------------------------------
C
      FUNCTION CHORD3D (X, Y, Z, I1, I2)
C
C     One-liner: Summed chord-lengths for X-Y-Z curve over range of indices
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        CHORD3D is the XYZ analog of CHORD (XY), needed by PLSFIT3D.
C     It computes the sum of the Euclidean distance between adjacent
C     points in a curve represented by three arrays. The calculation is
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
C     thought of as providing length increments.  The same applies to
C     CHORD3D.
C
C     Arguments:
C     ----------
C
C     Name  Type/Dimension  I/O/S  Description
C
C     X,Y,Z    R (*)        I      Coordinates representing the curve.
C
C     I1,I2    I            I      Indices for summation. The loop runs
C                                  from MIN(I1,I2) to MAX(I1,I2) so the
C                                  result is independent of order.
C
C     CHORD3D  R              O    Function value is the sum of the
C                                  chord lengths along the curve between
C                                  indices I1 and I2.
C
C     Environment:
C     ------------
C
C     VAX/VMS; FORTRAN 77; IMPLICIT NONE is non-standard.
C
C     Author:  Robert Kennelly, Sterling Federal Systems/NASA-Ames
C     -------
C
C     History:
C     --------
C
C     18 Feb. 1987    RAK    Initial design and coding of CHORD.
C     15 Aug. 1989  DAS/RAK  CHORD adapted as CHORD3D for XYZ case.
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
     &   CHORD3D, X (*), Y (*), Z (*)

C     Local variables.

      INTEGER
     &   I
      REAL
     &   DMAX, DX, DY, DZ, R

C     Execution.
C     ----------

      CHORD3D = ZERO

      DO 10, I = MIN (I1, I2), MAX (I1, I2) - 1
         DX   = ABS (X (I + 1) - X (I))
         DY   = ABS (Y (I + 1) - Y (I))
         DZ   = ABS (Z (I + 1) - Z (I))
         DMAX = MAX (DX, DY, DZ)
         IF (DMAX .GT. ZERO) THEN
            R = ONE / DMAX
            CHORD3D = CHORD3D + DMAX *
     &         SQRT ((DX * R) ** 2 + (DY * R) ** 2 + (DZ * R) ** 2)
         END IF
   10 CONTINUE

C     Termination.
C     ------------

      RETURN
      END

