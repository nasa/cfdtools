C+------------------------------------------------------------------------------
C
      FUNCTION DERIV3 (J, JL, JR, H, DEL)
C
C     One-liner: First derivative using 3-point formulae (central or not)
C     ----------
C
C     Description and usage:
C     ----------------------
C
C        DERIV3 modularizes the choice of central or left- or right-handed
C     3-point first derivative formulae, including the degenerate 2-point
C     case.  It works with the same 3-interval/4-point data structures as
C     the earlier BESSEL and THREEPT utilities, q.v., simplifying their use
C     for bicubic interpolation.  The data organization is as follows:
C
C        Point:       -1    0     1     2
C                  ---+-----+-----+-----+---
C        Interval:      -1     0     1
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
C     JL,     I               I      Pointers to the left and right
C     JR                             boundaries of the 4-point structure.
C                                    JL = -1 and JR = 2 except when a
C                                    boundary (or two) is involved.
C
C     H       R (-1:1)        I      Array of interval lengths.
C
C     DEL     R (-1:1)        I      Array of 2-point slopes for the intervals.
C                                     
C     DERIV3  R                 O    The function value is the desired
C                                    derivative.
C
C     Procedures:
C     -----------
C
C     BESSEL   Central 3-point first derivative
C     THREEPT  One-sided   "    "    "    "    (left or right)
C
C     Environment:  FORTRAN 77 with minor extensions
C     ------------
C
C     History:
C     --------
C
C     11/26/93   DAS   Adaptation of BESSEL for use in bicubic interpolation.
C
C     Author:  David Saunders, Sterling Software/NASA Ames, Moffett Field, CA
C     -------
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER
     >   J, JL, JR
      REAL
     >   H (-1:1), DEL (-1:1), DERIV3

C     Procedures:

      REAL
     >   BESSEL, THREEPT
      EXTERNAL
     >   BESSEL, THREEPT

C     Execution.
C     ----------

      IF (JL .EQ. 0 .AND. JR .EQ. 1) THEN  ! Degenerate case

         DERIV3 = DEL (0)

      ELSE IF (J .EQ. 0) THEN  ! Left end of middle interval

         IF (JL .EQ. 0) THEN
            DERIV3 = THREEPT (0, H, DEL)
         ELSE
            DERIV3 = BESSEL (0, H, DEL)
         END IF

      ELSE ! Right end of middle interval (J = 1)

         IF (JR .EQ. 1) THEN
            DERIV3 = THREEPT (1, H, DEL)
         ELSE
            DERIV3 = BESSEL (1, H, DEL)
         END IF

      END IF

      RETURN
      END
