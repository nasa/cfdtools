C+------------------------------------------------------------------------------
C
      SUBROUTINE BSAFNORM (C, XLE, YLE, CHORD, THICK,
     >                     XLENEW, YLENEW, CHORDNEW, THICKNEW)

C     ONE-LINER: B-Spline AirFoil NORMalizing utility
C                - -      -  -    ----
C     PURPOSE:
C
C        BSAFNORM adjusts the chord, leading edge, and maximum thickness
C     of an airfoil in DT_NURBS B-spline curve form.
C
C     METHOD:
C
C        The curve is assumed to represent the airfoil in clockwise
C     wrap-around form, from lower trailing edge to upper trailing edge.
C     The control points are shifted and/or scaled in-place, using the
C     general purpose utility BSCRVMOD.  BSAFGEOM can provide the input
C     curve's geometric properties.  The relevant adjustment formulas are:
C
C        X  <--  XLENEW  +  (X - XLE) * (CHORDNEW / CHORD)
C        Y  <--  YLENEW  +  (Y - YLE) * (THICKNEW / THICK)
C
C     where X and Y refer to the control point coordinates.
C
C     ENVIRONMENT:
C
C     VAX/VMS; FORTRAN 77, with...
C        > IMPLICIT NONE
C        > Trailing ! comments
C        > Names up to 8 characters
C
C     HISTORY:
C
C     02/26/93  D.A.Saunders  Initial implementation.
C     02/01/94       "        Dispensed with DOUBLE PRECISION.  Compile with
C                             /REAL_LENGTH=64 or -r8 switches as necessary.
C
C     AUTHOR: David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments
C     ---------

      REAL   C (*)      ! (I/O) DT_NURBS B-spline curve vector
                        !       representing an airfoil.  See BSAFGEOM.

      REAL   XLE, YLE   ! (I)   Present leading edge X and Y coordinates

      REAL   CHORD      ! (I)   Present chord = max (XTEL, XTEU) - XLE

      REAL   THICK      ! (I)   Present thickness = max (Yupper - Ylower)

      REAL   XLENEW, YLENEW     ! (I)   Desired leading edge coordinates

      REAL   CHORDNEW, THICKNEW ! (I)   Desired chord and thickness

C     Procedures
C     ----------

      EXTERNAL BSCRVMOD ! Modifies a 2-space B-spline curve.

C-------------------------------------------------------------------------------

C     Local variables
C     ---------------

      REAL   P (2)      ! Shift/scale parameters

C     Execution
C     ---------

C     Avoid possible roundoff associated with unnecessary scaling.

C     First, the chord and leading edge X.
C     BSCRVMOD's "LINEARX" tranformation is  X <-- P1 + P2 * X.

      IF (CHORDNEW .NE. CHORD) THEN
         P (2) = CHORDNEW / CHORD
         P (1) = XLENEW - XLE * P (2)

         CALL BSCRVMOD ('LINEARX', P, C)

      ELSE IF (XLENEW .NE. XLE) THEN
         P (1) = XLENEW - XLE

         CALL BSCRVMOD ('SHIFTX', P, C)

      END IF

C     Likewise for maximum thickness and leading edge Y:

      IF (THICKNEW .NE. THICK) THEN
         P (2) = THICKNEW / THICK
         P (1) = YLENEW - YLE * P (2)

         CALL BSCRVMOD ('LINEARY', P, C)

      ELSE IF (YLENEW .NE. YLE) THEN
         P (1) = YLENEW - YLE

         CALL BSCRVMOD ('SHIFTY', P, C)

      END IF

      RETURN
      END
