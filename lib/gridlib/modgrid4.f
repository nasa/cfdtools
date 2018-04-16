C+------------------------------------------------------------------------------
C
      SUBROUTINE MODGRID4 (NGRID, XGRID0, YGRID0, ZGRID0, TGRID0, IGRID,
     >                     IDIM, JDIM, I1, I2, J1, J2, XSURF, YSURF,
     >                     ZSURF, USURF, VSURF, ISURF, JSURF, EPS,
     >                     XGRID, YGRID, ZGRID, LUNOUT, IER)
C
C ONE-LINER: Modify a 3-space grid line after a discrete boundary perturbation
C
C PURPOSE:
C
C        MODGRID4 adjusts one radial line of a grid in 3-space given the
C     original grid line and (at its "low" end) a perturbed boundary surface
C     represented by a regular mesh of geometry data.  The intersection of
C     the original radial line with the new surface becomes a computational
C     grid point on the new surface.  Wing/body design-by-optimization is
C     the initial application.
C
C METHOD:
C
C        The intersection of the new surface and the original grid line
C     (possibly extrapolated) is determined as the solution of 3 nonlinear
C     equations, using parametric local spline techniques which avoid storage
C     of spline coefficients.  The original point distribution is then applied
C     to the modified arc-length of the radial line (and the modified surface
C     grid is defined by the revised end points of all radial grid lines).
C
C        For efficient use in a loop over many grid lines, the results of
C     one intersection calculation may be used as starting guesses for the
C     next - see the IGRID and ISURF/JSURF arguments.
C
C ENVIRONMENT:  FORTRAN 77, with minor enhancements
C
C HISTORY:
C
C     12/08/93  D.A.Saunders  MODGRID4 adapted from the earlier MODGRID2.
C
C AUTHOR: David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER NGRID            ! (I) Number of radial line grid points.

      REAL    XGRID0 (NGRID),  ! (I) Original radial grid line, with the
     >        YGRID0 (NGRID),  !     surface assumed to be near the low end.
     >        ZGRID0 (NGRID)

      REAL    TGRID0 (NGRID)   ! (S) Workspace for the corresponding
                               !     parameterization (normalized chords).

      INTEGER IGRID            ! (I) Points to the interval on the original
                               !     radial line at which the search for the
                               !     current point of intersection will start.
                               !     For the first call, use IGRID=1.
                               ! (O) Points to the actual interval, for possible
                               !     reuse as a starting guess on the next call.

      INTEGER IDIM, JDIM       ! (I)   Max. # points provided for in the
                               !       surface data arrays.

      INTEGER I1, I2, J1, J2   ! (I)   First & last indices in each grid
                               !       direction eligible for searching;
                               !       1 <= I1 < I1 + 1 < I2 <= IDIM, etc.

      REAL    XSURF (IDIM, JDIM)!(I)   Surface grid coordinates.
      REAL    YSURF (IDIM, JDIM)
      REAL    ZSURF (IDIM, JDIM)

      REAL    USURF (IDIM, JDIM)!(I/O) Parametric variables at the surface
      REAL    VSURF (IDIM, JDIM)!      grid points: calculated internally if
                               !       reqd. via USURF (I1, J1) < 0 on input,
                               !       else input from a prior call to MODGRID4
                               !       or PARAM2D.

      INTEGER ISURF, JSURF     ! (I/O) Indices of the surface cell containing
                               !       the point of intersection:
                               !       on input:  the values from a previous
                               !                  call, or I1, J1, if no better
                               !                  estimate is known;
                               !       on output: the "lower left" indices of
                               !                  the cell (unless IER > 0);
                               !       I1 <= ISURF < I2 on input & output;
                               !       likewise for JSURF.

      REAL    EPS              ! (I)   Tolerance for several convergence tests;
                               !       try 1.E-6 or 1.E-14 for 32- or 64-bit
                               !       arithmetic.

      REAL    XGRID (NGRID),   ! (O) The revised radial grid points.  MUST NOT
     >        YGRID (NGRID),   !     be the same locations as XGRID0/YGRID0 (*).
     >        ZGRID (NGRID)    !     X/Y/ZGRID (1) becomes part of the revised
                               !     computational mesh for the surface.

      INTEGER LUNOUT           ! (I) Logical unit for displaying iterations,
                               !     which are suppressed if LUNOUT < 0.
                               !     |LUNOUT| is used for error messages.

      INTEGER IER              ! (O) 0 means no problem was encountered;
                               !     other values as for INTSEC4, q.v.

C     Procedures.

      EXTERNAL  CHORDS3D       ! Chord-length utility.
      EXTERNAL  INTSEC4        ! Surface/curve intersection utility.
      EXTERNAL  PARAM2D        ! Parameterizes a surface (chord lengths).
      EXTERNAL  PLSCURVE       ! Local spline utility.  "Loose" fits are used.

C-------------------------------------------------------------------------------

C     Local constants.

      REAL      ONE, ZERO
      LOGICAL   CLOSED
      PARAMETER
     >  (ONE    = 1.E+0,
     >   ZERO   = 0.E+0,
     >   CLOSED = .FALSE.)

C     Local variables.

      INTEGER   I, IPREVS, J, JPREVS, LUNERR
      REAL      DERIVS (3), SCALE, SMATRIX (4, 4, 3), TI, TTOTAL,
     >          TINT, UINT, VINT, XINT, YINT, ZINT
      LOGICAL   NEW

C     Storage.

      SAVE      IPREVS, JPREVS, ! Viable since we can assume only one surface
     >          SMATRIX, TINT,  ! at a time is being processed
     >          UINT, VINT

C     Execution.

      IER = 0
      LUNERR = ABS (LUNOUT)

C     Parameterize the surface?

      IF (USURF (I1, J1) .LT. ZERO) THEN

         CALL PARAM2D (IDIM, JDIM, I1, I2, J1, J2, XSURF, YSURF, ZSURF,
     >                 USURF, VSURF)
         IPREVS = I1 - 1
         JPREVS = J1 - 1
         TINT = ZERO
         UINT = ZERO
         VINT = ZERO
      ELSE                  ! To be sure they're defined, this can't hurt
         IF (IGRID .EQ. 1)  TINT = ZERO
         IF (ISURF .EQ. I1) UINT = ZERO
         IF (JSURF .EQ. J1) VINT = ZERO
      END IF

C     Set up normalized cumulative chord lengths for the original radial curve.

      CALL CHORDS3D (NGRID, XGRID0, YGRID0, ZGRID0, .TRUE., TTOTAL,
     >               TGRID0)

C     Determine where it intersects the surface.

      CALL INTSEC4 (IDIM, JDIM, I1, I2, J1, J2, XSURF, YSURF, ZSURF,
     >              USURF, VSURF, ISURF, JSURF, IPREVS, JPREVS,
     >              EPS, SMATRIX, 1, NGRID, XGRID0, YGRID0, ZGRID0,
     >              TGRID0, IGRID, TTOTAL, CLOSED, TINT, UINT, VINT,
     >              XINT, YINT, ZINT, LUNOUT, IER)
      IF (IER .NE. 0) GO TO 800

C     Redistribute the grid points in the same relative way.

      XGRID (1) = XINT   ! Make sure of the new surface grid point
      YGRID (1) = YINT
      ZGRID (1) = ZINT
      SCALE = (ONE - TINT / TGRID0 (NGRID))
      NEW = .TRUE.

      DO I = 2, NGRID
         TI = TINT + SCALE * TGRID0 (I)

         CALL PLSCURVE (NGRID, XGRID0, YGRID0, ZGRID0, NEW, CLOSED,
     >                  TI, TTOTAL, J, XGRID (I), YGRID (I), ZGRID (I),
     >                  DERIVS)
      END DO

      GO TO 999


C     Error handling.

  800 WRITE (LUNERR, 1010) ' MODGRID4:  Bad return from INTSEC4'

  999 RETURN

C     Formats.

 1010 FORMAT (/, A)

      END
