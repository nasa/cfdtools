C+----------------------------------------------------------------------
C
      SUBROUTINE CURVDIS3D (NDAT, XDAT, YDAT, ZDAT, N, POWER, ISMOOTH,
     >                      LUNOUT, ARCS_ONLY, X, Y, Z, IER)
C
C     Acronym:  CURVature-based DIStribution: 3-Dimensional version
C               ----            ---           - -
C
C     Purpose:
C
C        CURVDIS3D redistributes points along a curve in 3-space which
C     is represented by a discrete dataset, not necessarily monotonic,
C     such that the local spacing is inversely proportional to the local
C     curvature (with adjustments to handle the zero curvature regions).
C     CURVDIS3D is intended for geometric applications where the coord-
C     inates have comparable units, although it may be used for non-
C     geometric purposes as well (e.g. plotting applications).  In both
C     cases, if the units of the curve coordinates are large (producing
C     small values of curvatures), the coordinates should be normalized
C     to the range [0,1].  For geometric curves, only one of the coord-
C     inates should be scaled to [0,1] and the others should be scaled
C     such that the true geometric shape is retained.
C
C
C     Method:
C
C        Given the availability of subroutine ARBDIS for translating
C     arbitrary "shape" functions into 1-D distributions whose spacing
C     at any abscissa is proportional to the shape function at that
C     abscissa, the problem then is to set up a shape function related
C     to curvature; ARBDIS does the rest.  (Abscissas in this case are
C     measures of arc length along the curve, and the shape function is
C     defined at each of the original data points.)
C        The subroutine is passed NDAT coordinates (XDAT, YDAT, ZDAT),
C     of the curve, for which it computes the cumulative arc lengths
C     (actually the cumulative chord lengths), s(n).  Next, FD12K is
C     called thrice to compute the 2nd derivatives (d2x/ds2, d2y/ds2 and
C     d2z/ds2) in order to approximate the local curvatures' magnitudes,
C     k(s), which are given by:
C
C     (1)  k(s) = sqrt ((d2x/ds2) ** 2 + (d2y/ds2) ** 2 + (d2z/ds2) ** 2
C
C        A direct inverse proportionality would yield a shape function
C     with smaller spacings in regions of larger curvature.  However,
C     for straight segments of the curve where k(s) = 0., a simple
C     inverse proportionality such as shape(s) = 1./k(s) would produce
C     an infinite local spacing.  Therefore, some arbitrary constant
C     (say 1.) should be added to k(s) in the denominator.
C        Clustering may be controlled by raising this fraction to some
C     exponent in the range [0,1].  Zero would produce shape(s) = 1., or
C     uniform spacing; 1 would yield "direct" inverse proportionality
C     spacing; and an exponent in between (0,1) would produce spacings
C     between these two extremes.  Thus, our spacing function becomes:
C
C        (2)   shape(s) = [1./(k(s) + 1.)] ** exponent  
C  
C        After the relative spacings (defining the shape function) are
C     calculated from eqn. (2), an option is provided to smooth the
C     shape function to reduce abrupt changes observed in sphere/cone
C     capsule defining sections.  The ARBDIS utility is then called to
C     compute redistributed arc lengths for the specified number N of
C     output points.  A further option is provided to smooth those
C     arc-lengths (with index-based abscissas) as may also be desirable.
C     LCSFIT then serves to evaluate the new X-Y-Z coordinates corres-
C     ponding to the new arc lengths.
C
C        Note:  Workspace is used/reused in the following sequence:
C
C     (1 : NDAT)                       Input curve chord lengths
C     (  NDAT + 1 : 2*NDAT)            FP values for FD12K calls
C     (2*NDAT + 1 : 3*NDAT)            d2x/ds2 values
C     (3*NDAT + 1 : 4*NDAT)            d2y/ds2 values
C     (4*NDAT + 1 : 5*NDAT)            d2z/ds2 values
C     (  NDAT + 1 : 2*NDAT)            xdat-ydat-zdat relative spacings
C                                      defining the shape function
C     (2*NDAT + 1 : 2*NDAT + N)        New curve chord lengths
C     (2*NDAT + N + 1 : 5*NDAT + 6*N)  Workspace for ARBDIS
C     (2*NDAT + N + 1 : 2*NDAT + 2*N)  Reused for index-based abscissas
C                                      if ARBDIS results are smoothed
C
C
C     Procedures:
C
C        ARBDIS           Distributes points in an interval according to
C                         an arbitrary shape function
C        FD12K            Derivatives by finite differencing
C        CHORDS3D         Cumulative chord lengths along a curve
C        DETECT_VERTICES  Identifies curvature spikes likely due to sharp
C                         vertices where some clustering is desirable
C        LCSFIT           Local cubic spline for the interpolations
C        SMOOTH1D_NITERS  Explicit smoothing utility with iteration control
C        VERTEX_CURVATURE Broadens and moderates the height of curvature spikes
C                         to enable some clustering towards sharp vertices
C
C
C     Error Handling:
C
C        See IER description below.  It is the user's responsibility
C     to check IER on return from CURVDIS3D.  Error messages are sent to
C     |LUNOUT|.
C
C
C     Environment:
C        FORTRAN 77 + minor extensions (originally)
C        Fortran 90 (now)
C
C     History:
C        08/12/88    D.A.Saunders   Initial design of CURVDIS (X-Y case)
C        09/09/88    B.A.Nishida    Initial implementation
C        12/20/88    M.D.Wong       Generalized to 3 dimensions (X-Y-Z)
C        02/15/96    DAS            Parametric cubic interpolation now,
C                                   instead of linear.
C        11/29/13    DAS            ISMOOTH argument introduced as done
C                                   earlier for CURVDIS.  Sharp vertices
C                                   are also handled the same way.
C
C     Authors:  David Saunders, Brian Nishida, Michael Wong
C               Sterling Software/NASA Ames Research Center, Mt. View, CA
C
C ----------------------------------------------------------------------


      IMPLICIT NONE

C     Arguments:
C     ----------

      INTEGER, INTENT (IN)  :: NDAT       ! Number of input curve points
      REAL,    INTENT (IN)  :: XDAT(NDAT) ! X coordinates of input curve
      REAL,    INTENT (IN)  :: YDAT(NDAT) ! Y coordinates of input curve
      REAL,    INTENT (IN)  :: ZDAT(NDAT) ! Z coordinates of input curve
      INTEGER, INTENT (IN)  :: N          ! Number of output curve points
      REAL,    INTENT (IN)  :: POWER      ! Exponent for spacing function;
                                          ! controls the point clustering.
                                          ! POWER must be in [0., 1.];
                                          ! 0.0 produces uniform spacing;
                                          ! 1.0 maximizes the curvature effect;
                                          ! 0.5 is suggested, but start with 1.0
                                          ! if CURVDIS3D2 is being employed, as
                                          ! is now recommended
      INTEGER, INTENT (IN)  :: ISMOOTH    ! 0 => no smoothing;
                                          ! 1 => smooth shape function only;
                                          ! 2 => smooth redistributed arcs only;
                                          ! 3 => perform both smoothings;
                                          ! use ISMOOTH = -1, -2, or -3 to allow
                                          ! automated treatment of apparently
                                          ! sharp vertices which otherwise are
                                          ! missed (1-pt. curvature spikes are
                                          ! essentially invisible to the scheme)
      INTEGER, INTENT (IN)  :: LUNOUT     ! Logical unit for showing convergence
                                          ! history; LUNOUT < 0 suppresses it
      LOGICAL, INTENT (IN)  :: ARCS_ONLY  ! T => skip calculations of X & Y(1:N)
                                          ! but return revised arcs as X(1:N)
      REAL,    INTENT (OUT) :: X(N),      ! X-Y-Z coordinates of output curve,
     >                         Y(N),      ! unless ARCS_ONLY = T, in which case
     >                         Z(N)       ! just the arc-lengths are returned as
                                          ! X(1:N)
      INTEGER, INTENT (OUT) :: IER        ! 0 => no errors;
                                          ! 1 => POWER is out of range;
                                          ! 2 => failure in ARBDIS


C----------------------------------------------------------------------


C     Local Constants:
C     ----------------

      INTEGER,   PARAMETER :: NVERTMAX = 12,     ! Limit on # vertices handled
     >                        NVMIN    = 5,      ! Limits on # vertex neighbors
     >                        NVMAX    = 15      ! each side of curvature spikes
      REAL,      PARAMETER :: ALPHA    = 1.E+0,  ! In case curvature = 0.
     >                        EXPONENT = 3.,     ! Controls shape fun. widths
     >                        FRACTION = 0.10,   ! Applied to spike heights
     >                        FRACTN   = 0.035,  ! Applied to NDAT for NV
     >                        ONE      = 1.E+0,
     >                        HALF     = 0.5+0,
     >                        ZERO     = 0.E+0
      LOGICAL,   PARAMETER :: FALSE    = .FALSE.,
     >                        TRUE     = .TRUE.
      CHARACTER, PARAMETER :: LINEAR*1 = 'L',
     >                        TIGHT*1  = 'M',
     >                        NAME*11  = 'CURVDIS3D: '

C     Local Variables:
C     ----------------

      INTEGER :: I, ISMTH, IV, LUNERR, NITERS, NV, NVERTEX
      INTEGER :: IVERTEX(NVERTMAX)
      REAL    :: DSQ, DSQMIN, ONEOVERN, STOTAL, XD, YD, ZD
      REAL    :: WORK(5*NDAT + 6*N)
      LOGICAL :: TREAT_VERTICES


C     Procedures:
C     -----------

      EXTERNAL :: ARBDIS, CHORDS3D, FD12K, LCSFIT, SMOOTH1D_NITERS
      EXTERNAL :: DETECT_VERTICES, VERTEX_CURVATURE


C     Execution:
C     ----------

      TREAT_VERTICES = ISMOOTH < 0             ! Retrofited option
      NVERTEX = 0                              ! See LCSFIT test near the end

      ISMTH  = ABS (ISMOOTH)
      LUNERR = ABS (LUNOUT)

      IF (POWER < ZERO .OR. POWER > ONE) THEN
         WRITE (LUNERR, '(/, 2A, ES12.4)')
     >      NAME, 'Bad POWER input outside [0, 1]:', POWER
         IER = 1
         GO TO 999
      END IF


C     Compute the cumulative chord lengths of the input data:

      CALL CHORDS3D (NDAT, XDAT, YDAT, ZDAT, FALSE, STOTAL, WORK(1))


C     Calculate the 2nd derivatives w.r.t. chord length:

      CALL FD12K (NDAT, WORK(1), XDAT, WORK(NDAT + 1), 
     >            WORK(2*NDAT + 1), WORK(2*NDAT + 1))

      CALL FD12K (NDAT, WORK(1), YDAT, WORK(NDAT + 1), 
     >            WORK(3*NDAT + 1), WORK(3*NDAT + 1))

      CALL FD12K (NDAT, WORK(1), ZDAT, WORK(NDAT + 1), 
     >            WORK(4*NDAT + 1), WORK(4*NDAT + 1))

 
C     Local curvature magnitudes:

      DO I = 1, NDAT
         WORK(NDAT + I) = SQRT (WORK (2*NDAT + I) ** 2 +
     >                          WORK (3*NDAT + I) ** 2 +
     >                          WORK (4*NDAT + I) ** 2)
      END DO


C     Attempt to cluster towards any vertices that will otherwise be missed?

      IF (TREAT_VERTICES) THEN

         CALL DETECT_VERTICES (NDAT, WORK(NDAT+1), NVERTMAX, IVERTEX,
     >                         NVERTEX)

         IF (NVERTEX > 0) THEN  ! Adjust curvature in-place; x/y data untouched
            NV = MIN (NVMAX, MAX (NVMIN, NINT (FRACTN * REAL (NDAT))))
            WRITE (LUNERR, '(A, I3, A, I3)')
     >         '# vertices found:', NVERTEX,
     >         '   # broadening neighbors either side:', NV

            DO IV = 1, NVERTEX
               CALL VERTEX_CURVATURE (NDAT, WORK(1), WORK(NDAT+1),
     >                              IVERTEX(IV), NV, EXPONENT, FRACTION)
            END DO
         END IF
      END IF

C     Moderate the local curvatures:

      IF (POWER == HALF) THEN  ! Avoid logarithms
         DO I = 1, NDAT
            WORK(NDAT+I) = ONE / SQRT (WORK(NDAT+I) + ALPHA)
         END DO
      ELSE
         DO I = 1, NDAT
            WORK(NDAT+I) = (WORK(NDAT+I) + ALPHA)**(-POWER)
         END DO
      END IF

CCCC  write (50, '(a)') '# Unsmoothed shape function'
CCCC  write (50, '(2es14.6)') (work(i), work(ndat+i), i = 1, ndat)


      IF (ISMTH == 1 .OR. ISMTH == 3) THEN

C        Smooth the shape function further with an explicit method:

         NITERS = MAX (4, MIN (NDAT/40, 15))
CCCC     write (*, '(a, i5)') 'NITERS for shape function:', NITERS

         CALL SMOOTH1D_NITERS (NITERS, 1, NDAT, WORK(1), WORK(NDAT+1))

CCCC     write (51, '(a)') '# smoothed shape function'
CCCC     write (51, '(2es14.6)') (work(i), work(ndat+i), i = 1, ndat)

      END IF


C     Redistribute the arc lengths based on the shape function:

      CALL ARBDIS (N, ZERO, STOTAL, NDAT, WORK(1), WORK(NDAT + 1), 'A',
     >             LUNOUT, WORK(2*NDAT + N + 1), WORK(2*NDAT + 1), IER)

      IF (IER /= 0) THEN
         WRITE (LUNERR, '(/, 2A, I4)') NAME, 'ARBDIS IER:', IER
         IER = 2
         GO TO 999
      END IF


      IF (ISMTH >= 2) THEN  ! Smooth the redistributed arc-lengths

         ONEOVERN = ONE / REAL (N)  ! Index-based abscissas in ARBDIS workspace
         DO I = 1, N
            WORK(2*NDAT+N+I) = REAL (I) * ONEOVERN
         END DO

         NITERS = MAX (5, MIN (N/20, 20))
CCCC     write (*, '(a, i5)') 'NITERS for redistributed arc lengths:',
CCCC >      NITERS

CCCC     write (52, '(a)') '# Unsmoothed redistributed arc lengths'
CCCC     write (52, '(2es14.6)')
CCCC >      (work(2*ndat+n+i), work(2*ndat+i), i = 1, n)

         CALL SMOOTH1D_NITERS (NITERS, 1, N, WORK(2*NDAT+N+1),
     >                         WORK(2*NDAT+1))
CCCC     write (53, '(a)') '# smoothed redistributed arc lengths'
CCCC     write (53, '(2es14.6)')
CCCC >      (work(2*ndat+n+i), work(2*ndat+i), i = 1, n)

      END IF


C     Redistribute the curve coordinates, unless only revised arcs are wanted:

      IF (ARCS_ONLY) THEN
         X(1:N) = WORK(2*NDAT+1:2*NDAT+N)
      ELSE

         CALL LCSFIT (NDAT, WORK (1), XDAT, TRUE, TIGHT, N,
     >                WORK(2*NDAT+1), X, WORK(2*NDAT+1+N))  ! Unused derivatives

         CALL LCSFIT (NDAT, WORK (1), YDAT, TRUE, TIGHT, N,
     >                WORK(2*NDAT+1), Y, WORK(2*NDAT+1+N))  !    "     "

         CALL LCSFIT (NDAT, WORK (1), ZDAT, TRUE, TIGHT, N,
     >                WORK(2*NDAT+1), Z, WORK(2*NDAT+1+N))  !    "     "

         IF (NVERTEX > 0) THEN  ! Even monotonic cubics are dubious at a vertex

            DO IV = 1, NVERTEX  ! Messy!  Find grid pt. nearest to data vertex
               I = IVERTEX(IV)  ! Data index
               XD = XDAT(I)
               YD = YDAT(I)
               ZD = ZDAT(I)
               DSQMIN = 1.E+6

               DO I = 1, N
                  DSQ = (X(I) - XD)**2 + (Y(I) - YD)**2 + (Z(I) - ZD)**2
                  IF (DSQ < DSQMIN) THEN
                      DSQMIN = DSQ
                      IVERTEX(IV) = I  ! Now a redistributed index
                  END IF
               END DO

               I = IVERTEX(IV) - 1  ! Reinterpolate XYZ linearly at I-1, I, I+1

               CALL LCSFIT (NDAT, WORK(1), XDAT, TRUE, LINEAR, 3,
     >                      WORK(2*NDAT+I), X(I), WORK(2*NDAT+1+N))
               CALL LCSFIT (NDAT, WORK(1), YDAT, TRUE, LINEAR, 3,
     >                      WORK(2*NDAT+I), Y(I), WORK(2*NDAT+1+N))
               CALL LCSFIT (NDAT, WORK(1), ZDAT, TRUE, LINEAR, 3,
     >                      WORK(2*NDAT+I), Z(I), WORK(2*NDAT+1+N))
            END DO

         END IF

      END IF


  999 RETURN

      END SUBROUTINE CURVDIS3D
