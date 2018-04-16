C+----------------------------------------------------------------------
C
      FUNCTION THREEPT (J, H, DEL)
C
C     One-liner: First derivative, non-central 3-point formula
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        Computes a first derivative approximation for PLSFIT over an
C     interval at a data boundary, using a forward or backward 3-point
C     formula.  The data must be in the form of arrays containing finite
C     difference interval lengths and 2-point forward difference deriva-
C     tives, and the differencing direction is controlled by a flag. See
C     PLSFIT for more details.
C
C        See module BUTLAND for a version with "shape-preserving"
C     adjustments.
C
C     Arguments:
C     ----------
C
C     Name    Type/Dimension  I/O/S  Description
C     J       I               I      Indicates at which end of the
C                                    interval the derivative is to be
C                                    estimated. J = 0 means left-hand
C                                    side, J = 1 means right. 
C
C     H       R (-1:1)        I      Array of interval lengths. The 0th
C                                    element is the length of the interval
C                                    on which the cubic is to be deter-
C                                    mined.
C
C     DEL     R (-1:1)        I      Array of derivative estimates. The
C                                    0th element is the forward difference
C                                    derivative over the interval on which
C                                    the cubic is to be determined.
C                                     
C     THREEPT R                 O    The function value is the derivative.
C
C     Environment:  VAX/VMS; FORTRAN 77
C     ------------
C
C     IMPLICIT NONE and 8-character symbolic names are non-standard.
C
C     Author:  Robert Kennelly, Sterling Federal Systems/NASA-Ames
C     -------
C
C     History:
C     --------
C
C     18 Feb. 1987    RAK    Initial design and coding.
C     06 June 1991    DAS    Original THREEPT renamed BUTLAND; THREEPT
C                            now gives unmodified 1-sided 3-pt. results.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     &   J
      REAL
     &   H (-1:1), DEL (-1:1), THREEPT

C     Local constants.

      REAL
     &   ONE
      PARAMETER
     &  (ONE = 1.0E+0)

C     Local variables.

      INTEGER
     &   STEP
      REAL
     &   WEIGHT

C     Execution.
C     ----------

C     Estimate first derivative on a left-hand boundary using a 3-point
C     forward difference (STEP = +1), or with a backward difference for
C     the right-hand boundary (STEP = -1).

      STEP = 1 - J - J   ! J here is consistent with related modules.

C     In {H, DEL} form, the derivative looks like a weighted average.

      WEIGHT  = -H (0) / (H (0) + H (STEP))
      THREEPT = WEIGHT * DEL (STEP) + (ONE - WEIGHT) * DEL (0)

C     Termination.
C     ------------

      RETURN
      END
