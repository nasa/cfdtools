C+------------------------------------------------------------------------------
C
      SUBROUTINE MODGRID3 (N, X0, Y0, T0, INDEX, C, TINT, X, Y, LUNOUT,
     >                     IER)
C
C ONE-LINER: Modify a 2-space grid line after a B-spline curve bndry. perturbn.
C
C PURPOSE:
C
C        MODGRID3 adjusts one radial line of a grid in 2-space given the
C     original grid line and the perturbed boundary at the "low" end, which
C     is represented by a B-spline curve in DT_NURBS form.  It was developed
C     for airfoil design-by-optimization applications.
C
C METHOD:
C
C        The intersection of the new boundary and the original grid line
C     (possibly extrapolated) is determined as the solution of two nonlinear
C     equations.  The original point distribution is then applied to the
C     modified arc-length of the radial line using parametric local spline
C     techniques which avoid storage of spline coefficients.
C
C        For efficient use in a loop over many grid lines, the results of
C     one intersection calculation may be used as starting guesses for the
C     next - see the INDEX and TINT arguments.
C
C ENVIRONMENT:
C
C     VAX/VMS; FORTRAN 77, with...
C     > IMPLICIT NONE
C     > Trailing ! comments
C     > Names up to 8 characters
C
C HISTORY:
C
C     06/29/93  D.A.Saunders  Initial implementation as MODGRD2D.
C     07/31/93       "        Renamed MODGRID3 when the discretized
C                             boundary (not B-spline curve) case was
C                             needed - see MODGRID2.
CC
C AUTHOR: David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER N                ! (I) Number of radial line grid points.

      REAL    X0 (N), Y0 (N)   ! (I) Original radial grid line, with the
                               !     boundary assumed to be near the low end.

      REAL    T0 (N)           ! (S) Workspace for the corresponding
                               !     cumulative chord lengths.

      INTEGER INDEX            ! (I) Input with 0 on the first call for a
                               !     given grid; normally input with the
                               !     output from the previous call otherwise.
                               !     (Points to the interval on the original
                               !     radial line at which the search for the
                               !     current point of intersection will start.)
                               ! (O) Points to the actual interval, for possible
                               !     reuse as a starting guess on the next call.

      REAL    C (*)            ! (I) B-spline curve (polynomial or rational) in
                               !     DT_NURBS form representing the boundary.

      REAL    TINT             ! (I) Estimate of the B-spline curve's dependent
                               !     variable at the point of intersection,
                               !     typically 0. on the first call and the
                               !     output from the previous call thereafter.
                               ! (O) The actual value for the present call.

      REAL    X (N), Y (N)     ! (O) The revised radial grid points.  MUST NOT
                               !     be the same locations as X0/Y0 (*).

      INTEGER LUNOUT           ! (I) Logical unit for displaying iterations,
                               !     which are suppressed if LUNOUT < 0.
                               !     |LUNOUT| is used for error messages.

      INTEGER IER              ! (O) 0 means no problem was encountered;
                               !     other values as for INTSEC3, q.v.

C     Procedures.

      REAL      CHORD
      EXTERNAL  CHORD          ! Chord-length utility.
      EXTERNAL  LCSFIT         ! Local spline utility.  "Loose" fits are used.

C-------------------------------------------------------------------------------

C     Local constants.

      REAL      ONE, ZERO
      LOGICAL   CALCT2, NEW
      CHARACTER METHOD * 1
      PARAMETER
     >  (CALCT2 = .FALSE.,     ! Avoids duplication by INTSEC3.
     >   METHOD = 'B',         ! "Bessel" = "loose" spline fits.
     >   ONE    = 1.E+0,
     >   NEW    = .TRUE.,      ! Always alternating between X, Y, etc.
     >   ZERO   = 0.E+0)

C     Local variables.

      INTEGER   I, I2, LUNERR
      REAL      SCALE, TANEW, TI, XINT, YINT

C     Execution.

      IER = 0
      LUNERR = ABS (LUNOUT)

      IF (INDEX .EQ. 0) THEN
         INDEX = 1
         TINT = ZERO
      END IF

C     Set up the cumulative chord lengths for the original radial curve.

      T0 (1) = ZERO
      DO 200, I = 2, N
         T0 (I) = T0 (I - 1) + CHORD (X0, Y0, I - 1, I)
  200 CONTINUE

C     Determine where it intersects the boundary curve.

      CALL INTSEC3 (C, CALCT2, N, X0, Y0, T0, I2, XINT, YINT, TINT,
     >              LUNOUT, IER)
      IF (IER .NE. 0) GO TO 800

C     Recover the parametric variable for the radial line at its new end.

      TANEW = T0 (I2) +
     >        SQRT ((XINT - X0 (I2)) ** 2 + (YINT - Y0 (I2)) ** 2)
      INDEX = I2

C     Redistribute the grid points in the same relative way.

      X (1) = XINT   ! Make sure of the boundary point
      Y (1) = YINT
      SCALE = (ONE - TANEW / T0 (N))

      DO 300, I = 2, N
         TI = TANEW + SCALE * T0 (I)

         CALL LCSFIT (N, T0, X0, NEW, METHOD, 1, TI, X (I), X (I))
         CALL LCSFIT (N, T0, Y0, NEW, METHOD, 1, TI, Y (I), Y (I))

  300 CONTINUE

      GO TO 999


C     Error handling.

  800 WRITE (LUNERR, 1010) ' MODGRID3:  Bad return from INTSEC3'

  999 RETURN

C     Formats.

 1010 FORMAT (/, A)

      END
