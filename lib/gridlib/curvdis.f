C+------------------------------------------------------------------------------
C
      SUBROUTINE CURVDIS (NDAT, XDAT, YDAT, N, POWER, ISMOOTH,
     >                    LUNOUT, ARCS_ONLY, X, Y, IER)
C
C     Acronym:  CURVature-based DIStribution (2-space)
C               ----            ---
C     Purpose:
C
C        CURVDIS redistributes points along a curve in 2-space which is
C     represented by a discrete dataset, not necessarily monotonic, such
C     that the local spacing is inversely proportional to the local
C     curvature (with adjustments to handle the zero curvature regions).
C     CURVDIS is intended for geometric applications where the coord-
C     inates have comparable units, although it may be used for non-
C     geometric purposes as well (e.g. plotting applications).  In both
C     cases, if the units of the curve coordinates are large (producing
C     small values of curvatures), the coordinates should be normalized
C     to the range [0,1].  For geometric curves, only one of the coord-
C     inate axes should be scaled to [0,1] and the other axis should be
C     scaled such that the true geometric shape is retained.
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
C        The subroutine is passed NDAT coordinates (XDAT, YDAT), of the
C     curve, for which it computes the cumulative arc lengths (or more
C     precisely, the cumulative chord lengths), s(n).  Next, FD12K is
C     called twice to compute the 2nd derivatives (d2x/ds2 and d2y/ds2)
C     in order to approximate the local curvatures' magnitudes, k(s),
C     which are given by:
C
C        (1)   k(s) = sqrt((d2x/ds2)**2 + (d2y/ds2)**2)
C
C        A direct inverse proportionality would yield a shape function
C     with smaller spacings in regions of larger curvature.  However,
C     for straight segments of the curve where k(s) = 0., a simple
C     inverse proportionality such as shape(s) = 1./k(s) would produce
C     an infinite local spacing.  Therefore, some arbitrary constant
C     (say 1.) should be added to k(s) in the denominator.
C        Clustering may be controlled by raising this fraction to some
C     exponent in the range [0,1].  0 would produce shape(s) = 1., or
C     uniform spacing; 1 would yield "direct" inverse proportionality
C     spacing; and an exponent in between (0,1) would produce spacings
C     between these two extremes.  Thus, our spacing function becomes:
C
C         (2)   shape(s) = [1./(k(s) + 1.)]**exponent
C
C        After the relative spacings (defining the shape function) are
C     calculated from eqn. (2), an option is provided to smooth the
C     shape function to reduce abrupt changes observed in sphere/cone
C     capsule defining sections.  The ARBDIS utility is then called to
C     compute redistributed arc lengths for the specified number N of
C     output points.  A further option is provided to smooth those
C     arc-lengths (with index-based abscissas) as also found desirable
C     for sphere/cone applications.  LCSFIT then serves to evaluate the
C     new X-Y coordinates corresponding to the new arc lengths.
C
C        This version also has the option to constrain cell growth rates
C     by smoothing them until they lie within [0.8, 1.2].  In order not
C     to change the calling sequence, ISMOOTH < 0 invokes that option,
C     as it also does for the vertex-capturing option.
C
C        NO: Constraining growth rates is a bad idea, defeating the
C     purpose of curvature-based spacing.  ISMOOTH < 0 is now used
C     instead to specify special treatment of sharp corners.
C
C        LATER:  ISMOOTH < 0 also controls broadening of high-curvature
C     regions onto flat/low curvature segments as found desirable in
C     flow calculations for entry capsules at the shoulder.
C
C     Work-space Usage:
C
C        (1)  WORK(1:NDAT)                 Input curve chord lengths
C        (2)  WORK(NDAT+1:2*NDAT)          FP values for FD12K calls
C        (3)  WORK(2*NDAT+1:3*NDAT)        d2x/ds2 values
C        (4)  WORK(3*NDAT+1:4*NDAT)        d2y/ds2 values
C        (5)  WORK(NDAT+1:2*NDAT)          xdat-ydat relative spacings
C                                          defining the shape function
C        (6)  WORK(2*NDAT+1:2*NDAT+N)      New curve chord lengths
C        (7)  WORK(2*NDAT+N+1:5*NDAT+6*N)  Work area used by ARBDIS
C        (8)  WORK(2*NDAT+N+1:2*NDAT+2*N)  Reused for index-based abscissas
C                                           if ARBDIS results are smoothed
C     Procedures:
C
C        ARBDIS           Distributes points in an interval according to an
C                         arbitrary shape function
C        FD12K            Calculates derivatives by finite differencing
C        CHORDS2D         Computes the cumulative chord lengths of a curve
C        DETECT_VERTICES  Identifies curvature spikes likely due to sharp
C                         vertices where some clustering is desirable
C        LCSFIT           Interpolates x and y to the redistributed arc
C                         lengths along the parameterized curve
C        SMOOTH1D_NITERS  Explicit smoothing utility with iteration control
C        VERTEX_CURVATURE Broadens and moderates the height of curvature spikes
C                         to enable some clustering towards sharp vertices
C     Error Handling:
C
C        See IER description below.  It is the user's responsibility
C     to check IER on return from CURVDIS.  Error messages are sent to |LUNOUT|.
C
C     History:
C        08/12/88    D.A.Saunders   Initial design.
C        09/09/88    B.A.Nishida    Initial implementation.
C        06/05/89    DAS            Refined the description a little.
C        02/15/96     "             Cubic interpolation rather than linear.
C        12/01/10    DAS, ERC, Inc. Inserted explicit smoothing of the shape
C                                   function, as suggested by the defining
C                                   section of a sphere/cone capsule forebody.
C        12/03/10     "    "        Index-based smoothing of the redistributed
C                                   arc-lengths helps for cases with extreme
C                                   jumps in curvature, such as a sphere/cone.
C                                   New argument ISMOOTH allows suppression of
C                                   either or both smoothings.
C                                   The work-space is no longer an argument, as
C                                   an automatic array is now preferable.
C        12/06/10     "    "        Option to return only the redistributed arc-
C                                   length distribution (as X(1:N)), as is more
C                                   convenient for the HEAT_SHIELD program.
C        08/01/12     "    "        Very large grids can mean the original limit
C                                   of 30 smoothing iterations may not be enough
C                                   yet n/8 may still be too many.  Compromise
C                                   by using n/10 with a limit of 100.
C        10/23/13     "    "        Now that CURVDIS2 has been introduced on top
C                                   of CURVDIS to normalize the data before any
C                                   curvature calculations, less smoothing seems
C                                   to be desirable.  (A loop in CURVDIS2 that
C                                   lowers the exponent by 0.1 if CURVDIS does
C                                   not converge also helps adhering to the
C                                   curvature-based shape function as much as
C                                   possible without failing.)
C        11/01/13     "    "        Option to try and resolve apparent sharp
C                                   vertices in the curve (ISMOOTH < 0).
C        11/04/13     "    "        NV = 6 needs to scale with NDAT.
C        11/07/13     "    "        Option to smooth cell growth rates until
C                                   they're within [0.8, 1.2] (ISMOOTH < 0 too).
C        11/08/13     "    "        NO!  Such smoothing can lose the connection
C                                   between the final arc length distribution
C                                   and the original data curvature distribu-
C                                   tion!  Suppress this misguided addition.
C        04/18/16    DAS, AMA, Inc. Poor resolution of a small shoulder radius
C                                   by CAPSULE_GRID prompted reverting to more
C                                   smoothing iterations that should allow a
C                                   higher exponent to converge.
C        06/27/16     "    "        The 04/18/16 change turned out to be overly
C                                   influenced by a generatrix from a CAD system
C                                   that had unusually dense resolution on the
C                                   shoulder.  There is no way of safeguarding
C                                   any smoothing algorithm from wildly varying
C                                   numbers of data points!
C        06/29/16     "    "        The latest attempt to improve capsule grid
C                                   point distributions in the presence of flat
C                                   segments (sphere/cone) seeks to spread the
C                                   redistributed points further on to the low
C                                   curvature regions as found desirable for
C                                   good resolution of spikes in surface heat
C                                   flux at the shoulder.  New utilities
C                                   DETECT_FLATNESS and HANDLE_FLATNESS are
C                                   activated as for the earlier attempt at
C                                   treating sharp corners (ISMOOTH < 0).
C
C     Authors:  David Saunders, Brian Nishida, Sterling Software
C               NASA Ames Research Center, Moffett Field, CA
C
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:
C     ----------

      INTEGER, INTENT (IN)  :: NDAT       ! Number of input curve points
      REAL,    INTENT (IN)  :: XDAT(NDAT) ! X coordinates of input curve
      REAL,    INTENT (IN)  :: YDAT(NDAT) ! Y coordinates of input curve
      INTEGER, INTENT (IN)  :: N          ! Number of output curve points
      REAL,    INTENT (IN)  :: POWER      ! Exponent for spacing function;
                                          ! controls the point clustering.
                                          ! POWER must be in [0., 1.];
                                          ! 0.0 produces uniform spacing;
                                          ! 1.0 maximizes the curvature effect;
                                          ! 0.5 is suggested, but start with 1.0
                                          ! if CURVDIS2 is being employed, as is
                                          ! now recommended
      INTEGER, INTENT (IN)  :: ISMOOTH    ! 0 => no smoothing;
                                          ! 1 => smooth shape function only;
                                          ! 2 => smooth redistributed arcs only;
                                          ! 3 => perform both smoothings;
                                          ! use ISMOOTH = -1, -2, or -3 to allow
                                          ! automated treatment of apparently
                                          ! sharp vertices which otherwise are
                                          ! missed (1-pt. curvature spikes are
                                          ! essentially invisible to the scheme)
                                          ! and also to constrain cell growth
                                          ! rates via further smoothing of them.
                                          ! NO: the latter is a BAD idea; see
                                          ! the History above.
      INTEGER, INTENT (IN)  :: LUNOUT     ! Logical unit for showing convergence
                                          ! history; LUNOUT < 0 suppresses it
      LOGICAL, INTENT (IN)  :: ARCS_ONLY  ! T => skip calculations of X & Y(1:N)
                                          ! but return revised arcs as X(1:N)
      REAL,    INTENT (OUT) :: X(N), Y(N) ! X & Y coordinates of output curve,
                                          ! unless ARCS_ONLY = T, in which case
                                          ! just the arc-lengths are returned as
                                          ! X(1:N)
      INTEGER, INTENT (OUT) :: IER        ! 0 => no errors;
                                          ! 1 => POWER is out of range;
                                          ! 2 => failure in ARBDIS

C------------------------------------------------------------------------------

C     Local Constants:
C     ----------------

      INTEGER,   PARAMETER :: NFLMAX   = 12,     ! Limit on # flat regions
     >                        NVERTMAX = 12,     ! Limit on # vertices handled
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
     >                        NAME*9   = 'CURVDIS: '

C     Local Variables:
C     ----------------

      INTEGER :: I, ISMTH, IV, LUNERR, NFL, NITERS, NV, NVERTEX
      INTEGER :: IFL(NFLMAX), IFR(NFLMAX), IVERTEX(NVERTMAX)
      REAL    :: DSQ, DSQMIN, ONEOVERN, STOTAL, XD, YD
      REAL    :: WORK(5*NDAT + 6*N)
      LOGICAL :: CONSTRAIN_GROWTH_RATES, TREAT_VERTICES

C     Procedures:
C     -----------

      EXTERNAL :: ARBDIS, CHORDS2D, FD12K, LCSFIT, SMOOTH1D_NITERS
      EXTERNAL :: DETECT_VERTICES, VERTEX_CURVATURE
      EXTERNAL :: DETECT_FLATNESS, HANDLE_FLATNESS

C     Execution:
C     ----------

      TREAT_VERTICES = ISMOOTH < 0             ! Retrofited option, now
                                               ! invoking better blending of
                                               ! zero-curvature regions likewise
CCC   CONSTRAIN_GROWTH_RATES = TREAT_VERTICES  ! Kludge for now
      CONSTRAIN_GROWTH_RATES = .FALSE.         ! See History above
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

      CALL CHORDS2D (NDAT, XDAT, YDAT, FALSE, STOTAL, WORK(1))

C     Calculate the 2nd derivatives, d2x/ds2 and d2y/ds2:

      CALL FD12K (NDAT, WORK(1), XDAT, WORK(NDAT+1), WORK(2*NDAT+1),
     >            WORK(2*NDAT+1))
      CALL FD12K (NDAT, WORK(1), YDAT, WORK(NDAT+1), WORK(3*NDAT+1),
     >            WORK(3*NDAT+1))

C     Local curvature magnitudes:

      DO I = 1, NDAT
         WORK(NDAT+I) = SQRT (WORK(2*NDAT+I)**2 + WORK(3*NDAT+I)**2)
      END DO

C     Attempt to cluster towards any vertices that will otherwise be missed?
C     Broadening of high curvature onto flat/low-curvature regions is also
C     automated here now, as needed for atmospheric entry capsule grids.

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

         CALL DETECT_FLATNESS (NDAT, WORK(NDAT+1), NFLMAX, IFL, IFR,
     >                         NFL)
         IF (NFL > 0) THEN
            WRITE (LUNERR, '(A, I3)')
     >         '# flat regions detected:', NFL, 'IFL, IFR:'
            WRITE (LUNERR, '(2I6)') (IFL(I), IFR(I), I = 1, NFL)

            CALL HANDLE_FLATNESS (NDAT, WORK(1), WORK(NDAT+1),
     >                            NFL, IFL, IFR)
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

CCCC     NITERS = MAX (5, MIN (NDAT/10, 30))
CCCC     NITERS = MAX (4, MIN (NDAT/40, 15))
         NITERS = MAX (9, MIN (NDAT/20, 20))
CCCC     write (*, '(a, i5)') 'NITERS for shape function:', NITERS

         CALL SMOOTH1D_NITERS (NITERS, 1, NDAT, WORK(1), WORK(NDAT+1))

CCCC     write (51, '(a)') '# smoothed shape function'
CCCC     write (51, '(2es14.6)') (work(i), work(ndat+i), i = 1, ndat)

      END IF

C     Redistribute the arc lengths based on the shape function:

      STOTAL = WORK(NDAT)

      CALL ARBDIS (N, ZERO, STOTAL, NDAT, WORK(1), WORK(NDAT+1), 'A',
     >             LUNOUT, WORK(2*NDAT+N+1), WORK(2*NDAT+1), IER)

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

CCCC     NITERS = MAX (5, MIN (N/10, 30))
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

         IF (CONSTRAIN_GROWTH_RATES) THEN   ! NOT A GOOD IDEA -- suppressed now
            CALL ADJUST_GROWTHS (N, WORK(2*NDAT+1))         ! Internal procedure
         END IF

         CALL LCSFIT (NDAT, WORK(1), XDAT, TRUE, TIGHT, N,
     >                WORK(2*NDAT+1), X, WORK(2*NDAT+1+N))  ! Unused derivatives

         CALL LCSFIT (NDAT, WORK(1), YDAT, TRUE, TIGHT, N,
     >                WORK(2*NDAT+1), Y, WORK(2*NDAT+1+N))  !    "     "

         IF (NVERTEX > 0) THEN  ! Even monotonic cubics are dubious at a vertex

            DO IV = 1, NVERTEX  ! Messy!  Find grid pt. nearest to data vertex
               I = IVERTEX(IV)  ! Data index
               XD = XDAT(I)
               YD = YDAT(I)
               DSQMIN = 1.E+6

               DO I = 1, N
                  DSQ = (X(I) - XD)**2 + (Y(I) - YD)**2
                  IF (DSQ < DSQMIN) THEN
                      DSQMIN = DSQ
                      IVERTEX(IV) = I  ! Now a redistributed index
                  END IF
               END DO

               I = IVERTEX(IV) - 1  ! Reinterpolate X, Y linearly at I-1, I, I+1

               CALL LCSFIT (NDAT, WORK(1), XDAT, TRUE, LINEAR, 3,
     >                      WORK(2*NDAT+I), X(I), WORK(2*NDAT+1+N))
               CALL LCSFIT (NDAT, WORK(1), YDAT, TRUE, LINEAR, 3,
     >                      WORK(2*NDAT+I), Y(I), WORK(2*NDAT+1+N))
            END DO
         END IF

      END IF

C     Error handling:
C     ---------------

  999 RETURN

C     Internal procedure for CURVDIS if ISMOOTH < 0:

      CONTAINS  ! No longer recommended; retained for posterity

C        -----------------------------------------------------------------------
C
         SUBROUTINE ADJUST_GROWTHS (N, S)  ! Arguments simplify nomenclature
C
C        Evaluate cell growth rates 2:N-1.  If any lie outside [GMIN, GMAX],
C        apply one smoothing iteration to them and calculate the corresponding
C        adjusted arc lengths.  Repeat as often as needed, up to a limit.
C
C        Later:  This is a bad idea because it can lose the connection between
C        the redistributed arcs and the original data curvature distribution.
C        It is suppressed now, but left in-line just in case it's of interest.
C        -----------------------------------------------------------------------

C        Arguments:

         INTEGER, INTENT (IN)    :: N      ! # points in arc-length distribution
         REAL,    INTENT (INOUT) :: S(N)   ! Arc lengths adjusted in-place

C        Local constants:

         INTEGER, PARAMETER :: ITMAX = 10  ! Limit on adjustment iterations
         REAL, PARAMETER    :: GMAX = 1.2, ! Cell growth rate limits
     >                         GMIN = 0.8

C        Local variables:

         INTEGER :: I, ITER, NBAD
         REAL    :: B(N), C(N), G(N), R(N)

C        Execution:

         DO ITER = 1, ITMAX  ! Until all growth rates are acceptable

            NBAD = 0
            G(1) = ONE
            DO I = 2, N - 1
               G(I) = (S(I+1) - S(I)) / (S(I) - S(I-1))
               IF (G(I) < GMIN .OR. G(I) > GMAX) NBAD = NBAD + 1
            END DO
            G(N) = ONE
CCCC        write (*, '(a, 2i4)') '  Itn. & # bad growth rates:', iter, nbad

C           Note that if growth rates are recovered later from the
C           redistributed (X, Y)s, they won't be identical to what we see here.

            IF (NBAD == 0) EXIT

            CALL SMOOTH1D_NITERS (1, 1, N, S, G)

C           Back out the arc lengths that correspond to the adjusted growths:

            DO I = 2, N - 1
CCC            A(I) =  G(I)          ! Subdiagonal is just G(2:N-1)
               B(I) = -(ONE + G(I))  ! Diagonal
               C(I) =  ONE           ! Superdiagonal
               R(I) =  ZERO          ! RHS, except for first & last elements
            END DO
            R(2)   = -S(1)*G(2)
            R(N-1) = -S(N)

            CALL TRDIAG (G(2), B(2), C(2), R(2), S(2), N-2)

         END DO

         END SUBROUTINE ADJUST_GROWTHS

      END SUBROUTINE CURVDIS
