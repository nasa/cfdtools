C+------------------------------------------------------------------------------
C
      SUBROUTINE MODGRID2 (NGRID, XGRID0, YGRID0, TGRID0, IGRID,
     >                     NBOUND, XBOUND, YBOUND, TBOUND, IBOUND,
     >                     CALCTB, XGRID, YGRID, LUNOUT, IER)
C
C ONE-LINER: Modify a 2-space grid line after a discrete boundary perturbation
C
C PURPOSE:
C
C        MODGRID2 adjusts one radial line of a grid in 2-space given the
C     original grid line and the perturbed boundary at the "low" end, which
C     is represented by discrete points.  It was developed for airfoil
C     design-by-optimization applications.
C
C METHOD:
C
C        The intersection of the new boundary and the original grid line
C     (possibly extrapolated) is determined as the solution of two nonlinear
C     equations, using parametric local spline techniques which avoid storage
C     of spline coefficients.  The original point distribution is then applied
C     to the modified arc-length of the radial line.
C
C        For efficient use in a loop over many grid lines, the results of
C     one intersection calculation may be used as starting guesses for the
C     next - see the IGRID and IBOUND arguments.
C
C        Also, the cumulative chord lengths TBOUND (*) may best be provided
C     by the higher level - see the CALCTB argument.
C
C ENVIRONMENT:  FORTRAN 77 with minor extensions.
C
C HISTORY:
C
C     07/31/93  D.A.Saunders  MODGRID2 adapted from the earlier MODGRID3.
C     12/08/93    "    "      Replaced CHORD usage with CHORDS2D, etc.
C
C AUTHOR: David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER NGRID            ! (I) Number of radial line grid points.

      REAL    XGRID0 (NGRID),  ! (I) Original radial grid line, with the
     >        YGRID0 (NGRID)   !     boundary assumed to be near the low end.

      REAL    TGRID0 (NGRID)   ! (S) Workspace for the corresponding
                               !     cumulative chord lengths.

      INTEGER IGRID            ! (I) Points to the interval on the original
                               !     radial line at which the search for the
                               !     current point of intersection will start.
                               !     For the first call, use IGRID = 1.
                               ! (O) Points to the actual interval, for possible
                               !     reuse as a starting guess on the next call.

      INTEGER NBOUND           ! (I) Number of boundary curve points.

      REAL    XBOUND (NBOUND), ! (I) The boundary represented as discrete pts.
     >        YBOUND (NBOUND)

      REAL    TBOUND (NBOUND)  ! (S) Workspace for the corresponding
                               !     cumulative chord lengths.  See CALCTB.

      INTEGER IBOUND           ! (I/O) Boundary pointer analogous to IGRID, q.v.

      LOGICAL CALCTB           ! (I) .FALSE. means TBOUND (*) is provided by
                               !     the calling program (since the boundary
                               !     curve may be fixed for many calls to
                               !     MODGRID2).  .TRUE. on a first call and
                               !     .FALSE. thereafter may be appropriate,
                               !     but be sure that TBOUND is SAVEd at the
                               !     higher level (or in COMMON).

      REAL    XGRID (NGRID),   ! (O) The revised radial grid points.  MUST NOT
     >        YGRID (NGRID)    !     be the same locations as XGRID0/YGRID0 (*).

      INTEGER LUNOUT           ! (I) Logical unit for displaying iterations,
                               !     which are suppressed if LUNOUT < 0.
                               !     |LUNOUT| is used for error messages.

      INTEGER IER              ! (O) 0 means no problem was encountered;
                               !     other values as for INTSEC2, q.v.

C     Procedures.

      EXTERNAL  CHORDS2D       ! Chord-length utility.
      EXTERNAL  INTSEC2        ! Curve/curve intersection utility.
      EXTERNAL  LCSFIT         ! Local spline utility.  "Loose" fits are used.

C-------------------------------------------------------------------------------

C     Local constants.

      REAL      ONE, ZERO
      LOGICAL   CALCT2, NEW
      CHARACTER METHOD * 1
      PARAMETER
     >  (CALCT2 = .FALSE.,     ! Avoids duplication by INTSEC2.
     >   METHOD = 'B',         ! "Bessel" = "loose" spline fits.
     >   ONE    = 1.E+0,
     >   NEW    = .TRUE.,      ! Always alternating between X, Y, etc.
     >   ZERO   = 0.E+0)

C     Local variables.

      INTEGER   I, I2, LUNERR
      REAL      SCALE, TANEW, TI, TTOTAL, XINT, YINT

C     Execution.

      IER = 0
      LUNERR = ABS (LUNOUT)

C     Set up the cumulative chord lengths for the original radial curve.

      CALL CHORDS2D (NGRID, XGRID0, YGRID0, .FALSE., TTOTAL, TGRID0)

C     Determine where it intersects the boundary curve.

      CALL INTSEC2 (NBOUND, XBOUND, YBOUND, TBOUND, IBOUND, CALCTB,
     >              NGRID,  XGRID0, YGRID0, TGRID0, IGRID,  CALCT2,
     >              XINT, YINT, LUNOUT, IER)
      IF (IER .NE. 0) GO TO 800

C     Recover the parametric variable for the radial line at its new end.

      TANEW = TGRID0 (IGRID) +
     >        SQRT ((XINT - XGRID0 (IGRID)) ** 2 +
     >              (YINT - YGRID0 (IGRID)) ** 2)

C     Redistribute the grid points in the same relative way.

      XGRID (1) = XINT   ! Make sure of the boundary point
      YGRID (1) = YINT
      SCALE = (ONE - TANEW / TTOTAL)

      DO 300, I = 2, NGRID
         TI = TANEW + SCALE * TGRID0 (I)

         CALL LCSFIT (NGRID, TGRID0, XGRID0, NEW, METHOD, 1, TI,
     >                XGRID (I), XGRID (I))
         CALL LCSFIT (NGRID, TGRID0, YGRID0, NEW, METHOD, 1, TI,
     >                YGRID (I), YGRID (I))
  300 CONTINUE

      GO TO 999


C     Error handling.

  800 WRITE (LUNERR, 1010) ' MODGRID2:  Bad return from INTSEC2'

  999 RETURN

C     Formats.

 1010 FORMAT (/, A)

      END
