C+------------------------------------------------------------------------------
C
      SUBROUTINE BSAFGEOM (C, XLE, YLE, TLE, XTEL, YTEL, XTEU, YTEU,
     >                     CHORD, THICK, XTHICK, CAMBER, XCAMBER, AREA,
     >                     RADIUS, LUNOUT, IER)

C ONE-LINER: B-Spline AirFoil GEOMetry calculations
C            - -      -  -    ----
C PURPOSE:
C
C        BSAFGEOM determines some standard geometric properties of an
C     airfoil in DT_NURBS B-spline curve form.
C
C METHOD:
C
C        The curve is assumed to represent the airfoil in clockwise
C     wrap-around form, from lower trailing edge to upper trailing edge.
C     General purpose utility BSZEROMN is used to locate the leading
C     edge by minimizing X with respect to T.  For the thickness and
C     camber, FMINRC is used to maximize | Yupper(X) +/- Ylower(X) |.
C     These calculations use BSEVAL on each surface to determine Y
C     as a function of X (itself a zero-finding problem).
C
C        CAMBER is returned as the average of the upper and lower surface
C     Ys; THICK is not halved; and neither is divided by CHORD.
C
C        NOTE: These are the usual approximations to maximum airfoil
C     thickness and camber.  More precise determinations involving the
C     true mean-line and appropriate normals may well be tractable with
C     the B-spline curve representation, but not this time around...
C
C ENVIRONMENT:
C     VAX/VMS; FORTRAN 77, with...
C     > IMPLICIT NONE
C     > Trailing ! comments
C     > Names up to 8 characters
C
C HISTORY:
C     05/16/92  D.A.Saunders  Initial implementation.
C     05/28/92       "        Added area calculation.
C     04/28/93       "        Added leading edge radius, and printing option.
C     02/01/94       "        Dispensed with DOUBLE PRECISION.  Compile with
C                             /REAL_LENGTH=64 or -r8 switches as necessary.
C
C AUTHOR: David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments
C     ---------

      REAL    C (*)     ! (I)   DT_NURBS B-spline curve vector:
                        !       C(1) = 1 (# parametric variables);
                        !       C(2) = 2 (for X and Y);
                        !       C(3) = polynomial order k (4 for cubics);
                        !       C(4) = # control points, N;
                        !       C(5) = pointer for efficient evaluations;
                        !       C(6 : 5 + N + k) = knot vector,
                        !       followed by X coords. of control pts.,
                        !       followed by Y coords. of control pts.

      REAL    XLE, YLE, TLE     ! (O)   Leading edge X, Y, and T coordinates.

      REAL    XTEL, YTEL !(O)   Lower trailing edge X and Y coordinates.

      REAL    XTEU, YTEU !(O)   Upper trailing edge X and Y coordinates.

      REAL    CHORD     ! (O)   Chord = max (XTEL, XTEU) - XLE.

      REAL    THICK     ! (O)   Max. thickness = max (Yupper - Ylower).

      REAL    XTHICK    ! (O)   Corresponding X.

      REAL    CAMBER    ! (O)   Max. camber = max | Yupper + Ylower| / 2
                        !       (relative to the X-axis, with appropriate sign).

      REAL    XCAMBER   ! (O)   Corresponding X.

      REAL    AREA      ! (O)   Sectional area.

      REAL    RADIUS    ! (O)   Leading edge radius of curvature.

      INTEGER LUNOUT    ! (I)   Logical unit for display of results.
                        !       LUNOUT < 0 suppresses the writes.
                        !       Error messages go to |LUNOUT|.

      INTEGER IER       ! (O)   Success/error code:
                        !        0 means no problem was encountered;
                        !        ? <to be completed>

C     Procedures
C     ----------

      EXTERNAL BSAREA   ! Integrate a B-spline curve (Y w.r.t. X).
      EXTERNAL BSEVAL   ! Evaluates a B-spline curve (Y for given X).
      EXTERNAL BSZEROMN ! Finds a zero, min, or max of a B-spline curve.
      EXTERNAL DTRADC   ! Radius of curvature utility.
      EXTERNAL FMINRC   ! 1-D minimizer.

C-------------------------------------------------------------------------------

C     Local constants
C     ---------------

      INTEGER
     >   MAXFUN, MAXORDER, NWORK
      REAL
     >   HALF, TOL
      CHARACTER
     >   NAME * 8
      PARAMETER
     >  (MAXFUN = 50,          ! Limit on minimizer's # evaluations
     >   MAXORDER = 7,        ! Max. degree here is 6 for work-space reasons
     >   NAME   = 'BSAFGEOM',
     >   NWORK  = MAXORDER * (MAXORDER + 5),  ! Enough for 2nd derivs. of X & Y
     >   HALF   = 0.5E+0,
     >   TOL    = 0.E+0)       ! Minimizer tolerance

C     Local variables
C     ---------------

      INTEGER
     >   I, ISTAT, IX, IY, LUNERR, NCPTS, NUMFUN
      REAL
     >   F, QUAD, T, T1, T2, WORK (NWORK), X, X1, X2, XY (2), YL, YU
      LOGICAL
     >   NEW

C     Execution
C     ---------

      NCPTS = INT (C (4))              ! # control pts.
      IX = 5 + INT (C (3)) + NCPTS     ! Offsets of control pt. Xs and Ys
      IY = IX + NCPTS
      T1 = C (6)                       ! First knot
      T2 = C (IX)                      ! Last knot
      LUNERR = ABS (LUNOUT)

C     Find the leading edge value of T by minimizing X w.r.t. T:

      CALL BSZEROMN ('MIN', 1, C, T1, T2, TLE, XY, LUNERR, IER)

      IF (IER .NE. 0) THEN
         WRITE (LUNERR, 1000) 'BSZEROMN'
         ISTAT = -3     ! Not meangingful, but error msg. needs a value.
         GO TO 900
      END IF

      XLE  = XY (1)
      YLE  = XY (2)
      XTEL = C (IX + 1)
      XTEU = C (IX + NCPTS)
      CHORD = MAX (XTEL, XTEU) - XY (1)
      YTEL = C (IY + 1)
      YTEU = C (IY + NCPTS)

C     The maximum-finding iterations for thickness and camber are similar.

      X1 = XLE
      X2 = MIN (XTEL, XTEU)
      NEW = .TRUE.             ! Spline, that is (for BSEVAL)

      DO 300, I = 1, 2

         ISTAT = 2             ! 2 initializes the iteration
         NUMFUN = MAXFUN

  200    CONTINUE

            CALL FMINRC (X1, X2, X, F, TOL, NUMFUN, NAME,-LUNERR, ISTAT)

            IF (ISTAT .EQ. -1) THEN  ! Probable application error

               IER = 1
               WRITE (LUNERR, 1000) 'FMINRC iteration limit reached'
               GO TO 900

            ELSE IF (ISTAT .LT. 0) THEN  ! Some other fatal error

               IER = 2
               WRITE (LUNERR, 1000) 'FMINRC'
               GO TO 900

            ELSE IF (ISTAT .GT. 0) THEN  ! Evaluate the function

               CALL BSEVAL ('Y', C, NEW, 1, X, T1, TLE, T, YL, LUNERR,
     >                      IER)
               IF (IER .NE. 0) THEN
                  WRITE (LUNERR, 1000) 'BSEVAL (lower)'
                  GO TO 900
               END IF

               NEW = .FALSE.

               CALL BSEVAL ('Y', C, NEW, 1, X, TLE, T2, T, YU, LUNERR,
     >                      IER)
               IF (IER .NE. 0) THEN
                  WRITE (LUNERR, 1000) 'BSEVAL (upper)'
                  GO TO 900
               END IF

               IF (I .EQ. 1) THEN    ! F is -Thickness
                  F = YL - YU
               ELSE                  ! F is -|Camber|
                  F = -ABS (YU + YL) 
               END IF

               GO TO 200

            ELSE   ! ISTAT = 0 (success).  Ensure everything is at the best X.

               CALL BSEVAL ('Y', C, NEW, 1, X, T1, TLE, T, YL, LUNERR,
     >                      IER)

               CALL BSEVAL ('Y', C, NEW, 1, X, TLE, T2, T, YU, LUNERR,
     >                      IER)

               IF (I .EQ. 1) THEN
                  THICK  = YU - YL
                  XTHICK = X
               ELSE
                  CAMBER  = (YU + YL) * HALF
                  XCAMBER = X
               END IF

            END IF

  300 CONTINUE


C     Estimate the area using  Integral y(x) dx  =  Integral y(t) x'(t) dt

      CALL BSAREA (C, TLE, T2, AREA, NWORK, WORK, LUNERR, IER)

      IF (IER .NE. 0) THEN
         WRITE (LUNERR, 1000) 'BSAREA (upper)'
         GO TO 900
      END IF

      CALL BSAREA (C, TLE, T1, QUAD, NWORK, WORK, LUNERR, IER)

      IF (IER .NE. 0) THEN
         WRITE (LUNERR, 1000) 'BSAREA (lower)'
         GO TO 900
      END IF

      AREA = AREA - QUAD    ! Various signs work out this way


C     Radius of curvature at the leading edge:

      CALL DTRADC (C, TLE, NWORK, WORK, RADIUS, IER)

      IF (IER .NE. 0) THEN
         WRITE (LUNERR, 1000) 'DTRADC'
         GO TO 900
      ELSE
         RADIUS = ABS (RADIUS)
      END IF

      IF (LUNOUT .GT. 0) WRITE (LUNOUT, 1020)
     >   'Leading edge X:       ', XLE,
     >   'Leading edge Y:       ', YLE,
     >   'Leading edge radius:  ', RADIUS,
     >   'Lower trailing edge X:', XTEL,
     >   'Lower trailing edge Y:', YTEL,
     >   'Upper trailing edge X:', XTEU,
     >   'Upper trailing edge Y:', YTEU,
     >   'Chord:                ', CHORD,
     >   'Maximum thickness:    ', THICK,
     >   'Corresponding X:      ', XTHICK,
     >   'Maximum camber:       ', CAMBER,
     >   'Corresponding X:      ', XCAMBER,
     >   'Area:                 ', AREA,
     >   ' '

      GO TO 999


C     Error handling.

  900 WRITE (LUNERR, 1010) IER, ISTAT

  999 RETURN


C     Formats.

 1000 FORMAT (/, ' *** BSAFGEOM: Bad return from ', A)
 1010 FORMAT (' IER: ', I4, '   ISTAT: ', I3)
 1020 FORMAT (/, (1X, A, 1X, F24.16))

      END
