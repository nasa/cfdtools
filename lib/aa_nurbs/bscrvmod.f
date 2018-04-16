C+------------------------------------------------------------------------------
C
      SUBROUTINE BSCRVMOD (ACTION, PARAMS, C)
C
C     One-liner:
C
C        Modifies (the control points of) a 2-space B-spline curve (in-place)
C
C     Description:
C
C           BSCRVMOD modifies a B-spline XY curve by applying the indicated
C        transformation to the control points of its DT_NURBS representation.
C        Scaling, shifting, rotating, etc., are the intended transformations.
C        Some transformations apply to the X or the Y control points; others
C        apply to both.  Numerical values for the transformation are expected
C        in the PARAMS argument as detailed below.
C
C     Method:
C
C           ACTION takes the form SCALE|SHIFT|ROTATE|...[X|Y|XY] for
C        ease of use if not efficiency.  An invalid string here causes
C        the subroutine to STOP since such a programming error should
C        be corrected.  An error return argument is thereby avoided.
C        No error checking of C (*) is done at this level.
C
C           ACTION and PARAMS are to be used as follows, where X means
C        all X coordinates of the control points, etc.
C
C           SCALEX     X <-- PARAMS (1) * X
C           SCALEY     Y <-- PARAMS (1) * Y
C           SCALEXY    X <-- PARAMS (1) * X  &  Y <-- PARAMS (1) * Y
C
C           SHIFTX     X <-- PARAMS (1) + X
C           SHIFTY     Y <-- PARAMS (1) + Y
C
C           LINEARX    X <-- PARAMS (1) + PARAMS (2) * X
C           LINEARY    Y <-- PARAMS (1) + PARAMS (2) * Y
C
C           ROTATE     Each (X, Y) control point is rotated PARAM (1)
C                      degrees about the point (PARAMS (2), PARAMS (3)).
C                      Positive is anticlockwise.
C
C           <To be extended as the needs arise>
C
C     Environment:
C
C        VAX/VMS; FORTRAN 77 with minor extensions.
C
C     History:
C
C        20 Feb. 1993   DAS   Initial implementation, with lofting of
C                             wing sections in mind.
C        25 Feb. 1993    "    Allowing for distinct input and output
C                             C (*) arrays was inefficient - do the
C                             copying at a higher level only if needed.
C        01 Feb. 1994    "    Dispensed with DOUBLE PRECISION.  Compile with
C                             /REAL_LENGTH=64 or -r8 switches as necessary.
C
C     Author:  David Saunders, Sterling Software/NASA Ames, Mt. View, CA
C
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      CHARACTER ACTION * (*)   ! Self-descriptive upper case string
      REAL      PARAMS (*)     ! Input parameters of the transformation
      REAL      C (*)          ! Input/output B-spline curve

C-------------------------------------------------------------------------------

C     Local variables:

      INTEGER   I, J, IX1, IX2, IY1, IY2, KORDER, NCPTS
      REAL      ANGLE, CI, SI, XP, YP

C     Execution:

      KORDER = INT (C (3))          ! Degree + 1
      NCPTS  = INT (C (4))          ! # control pts.
      IX1 = 6 + NCPTS + KORDER      ! Index of first control pt. X
      IY1 = IX1 + NCPTS             ! .......................... Y
      IX2 = IY1 - 1                 ! .........last............. X
      IY2 = IX2 + NCPTS             ! .......................... Y

C     Duplicating some code here makes it easier to catch an invalid ACTION:

      IF (ACTION .EQ. 'SCALEX') THEN

         DO I = IX1, IX2
            C (I) = PARAMS (1) * C (I)
         END DO

      ELSE IF (ACTION .EQ. 'SCALEY') THEN

         DO I = IY1, IY2
            C (I) = PARAMS (1) * C (I)
         END DO

      ELSE IF (ACTION .EQ. 'SCALEXY') THEN

         DO I = IX1, IY2
            C (I) = PARAMS (1) * C (I)
         END DO

      ELSE IF (ACTION .EQ. 'SHIFTX') THEN

         DO I = IX1, IX2
            C (I) = PARAMS (1) + C (I)
         END DO

      ELSE IF (ACTION .EQ. 'SHIFTY') THEN

         DO I = IY1, IY2
            C (I) = PARAMS (1) + C (I)
         END DO

      ELSE IF (ACTION .EQ. 'LINEARX') THEN

         DO I = IX1, IX2
            C (I) = PARAMS (1) + PARAMS (2) * C (I)
         END DO

      ELSE IF (ACTION .EQ. 'LINEARY') THEN

         DO I = IY1, IY2
            C (I) = PARAMS (1) + PARAMS (2) * C (I)
         END DO

      ELSE IF (ACTION .EQ. 'ROTATE') THEN

         ANGLE = PARAMS (1) * ASIN (1.D+0) / 90.D+0
         CI = COS (ANGLE)
         SI = SIN (ANGLE)

         DO I = IX1, IX2
            J = I + NCPTS
            XP = C (I) - PARAMS (2)
            YP = C (J) - PARAMS (3)
            C (I) = XP * CI - YP * SI + PARAMS (2)
            C (J) = XP * SI + YP * CI + PARAMS (3)
         END DO

      ELSE

         STOP 'BSCRVMOD: Invalid ACTION'

      END IF

C     Termination:

      RETURN
      END
