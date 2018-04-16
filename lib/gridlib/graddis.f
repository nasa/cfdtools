C+------------------------------------------------------------------------------
C
      SUBROUTINE GRADDIS (NDAT, XDAT, YDAT, SDAT, GDAT, N, POWER,
     >                    ISMOOTH, LUNOUT, X, Y, S, IER)
C
C     Acronym:  GRADient-based DIStribution (2-space)
C               ----           ---
C     Purpose:
C
C        GRADDIS is a gradient-based analogue of curvature-based CURVDIS
C     that is expected to be used in a loop over POWER by GRADDIS2, q.v.
C
C        GRADDIS redistributes points along a curve in 2-space which is
C     represented by a discrete dataset, not necessarily monotonic in x,
C     such that the local spacing is inversely proportional to the local
C     gradient magnitude (with appropriate handling of zero gradients).
C     GRADDIS is intended for geometric applications where the coord-
C     inates have comparable units, although it may be used for non-
C     geometric purposes as well (e.g. plotting applications).  Calling
C     it via GRADDIS2 is strongly recommended, as it deals with data
C     scaling, gradient calculation, and retries if the iterative method
C     of redistribution does not converge.
C
C     Method:
C
C        Given the availability of subroutine ARBDIS for translating
C     arbitrary "shape" functions into 1-D distributions whose spacing
C     at any abscissa is proportional to the shape function at that
C     abscissa, the problem then is to set up a shape function related
C     to gradient; ARBDIS does the rest.  (Abscissas in this case are
C     measures of arc length along the curve, and the shape function is
C     defined at each of the original data points.)
C        The same shape function as used by CURVDIS is appropriate,
C     except curvature is replaced by the input gradient data.
C
C                shape(s) = [1./(gdat(s) + 1.)]**exponent
C
C     where gdat(:) = |dF/ds| for some some function f.
C
C        After the relative spacings (defining the shape function) are
C     calculated as shown above, an option is provided to smooth the
C     shape function, which may contain spikes due to shock waves, say.
C     The ARBDIS utility is then called to compute redistributed arc
C     lengths for the specified number N of output points.  A further
C     option is provided to smooth those arc-lengths (with index-based
C     abscissas), although this seems less desirable for flow field
C     temperature profiles across shocks than it does for curvature-
C     based redistributions.  LCSFIT then serves to evaluate the new
C     X-Y coordinates corresponding to the new arc lengths.
C
C        The work space is used/reused in the following sequence:
C
C     (1)  WORK(1:NDAT)               xdat-ydat relative spacings
C                                     defining the shape function
C     (2)  WORK(NDAT+1:4*NDAT+5*N)    Work area used by ARBDIS
C     (3)  WORK(NDAT+1:NDAT+N)        Reused for index-based abscissas
C                                     if ARBDIS results are smoothed
C
C     Procedures:
C
C        ARBDIS           Distributes points in an interval according to an
C                         arbitrary shape function
C        LCSFIT           Interpolates x and y to the redistributed arc
C                         lengths along the parameterized curve
C        SMOOTH1D_NITERS  Explicit smoothing utility with iteration control
C
C     Error Handling:
C
C        See IER description below.  It is the user's responsibility
C     to check IER on return from GRADDIS.  Error messages are sent to |LUNOUT|.
C
C     History:
C        08-09/88    D.A.Saunders,  Initial implementation of CURVDIS.
C                    B.A.Nishida
C        10/23/13    D.A.Saunders   Version of CURVDIS used as GRADDIS start.
C        11/20/13      "     "      Initial adaptation of CURVDIS as GRADDIS.
C        12/02/13      "     "      Not returning the updated arc lengths
C                                   was poor design: recovering them from the
C                                   new (x,y)s is not identical to using the
C                                   precise arc lengths that produced those
C                                   new (x,y)s (for function interpolations).
C
C     Author:  David Saunders, ERC, Inc. at NASA Ames Research Center, CA
C
C ------------------------------------------------------------------------------


      IMPLICIT NONE

C     Arguments:
C     ----------

      INTEGER, INTENT (IN)  :: NDAT       ! Number of input curve points
      REAL,    INTENT (IN)  :: XDAT(NDAT) ! X coordinates of input curve
      REAL,    INTENT (IN)  :: YDAT(NDAT) ! Y coordinates of input curve
      REAL,    INTENT (IN)  :: SDAT(NDAT) ! Corresponding arc lengths
      REAL,    INTENT (IN)  :: GDAT(NDAT) ! |dF/ds| as explained above
      INTEGER, INTENT (IN)  :: N          ! Number of output curve points
      REAL,    INTENT (IN)  :: POWER      ! Exponent for spacing function;
                                          ! controls the point clustering.
                                          ! POWER is in [0., 1.] (or higher);
                                          ! 0.0 produces uniform spacing;
                                          ! 1.0 maximizes the gradient effect;
                                          ! starting with 1.0 via use of
                                          ! GRADDIS2 is recommended
      INTEGER, INTENT (IN)  :: ISMOOTH    ! 0 => no smoothing;
                                          ! 1 => smooth shape function only;
                                          ! 2 => smooth redistributed arcs only;
                                          ! 3 => perform both smoothings;
                                          ! 1 is recommended for profiles
                                          ! involving shocks and boundary layers
      INTEGER, INTENT (IN)  :: LUNOUT     ! Logical unit for showing convergence
                                          ! history; LUNOUT < 0 suppresses it
      REAL,    INTENT (OUT) :: X(N), Y(N) ! X & Y coordinates of output curve
      REAL,    INTENT (OUT) :: S(N)       ! Corresponding new arc lengths
      INTEGER, INTENT (OUT) :: IER        ! 0 => no errors;
                                          ! 1 => POWER is out of range;
                                          ! 2 => failure in ARBDIS


C------------------------------------------------------------------------------


C     Local Constants:
C     ----------------

      REAL,      PARAMETER :: ALPHA    = 1.E+0,  ! In case gradient = 0.
     >                        ONE      = 1.E+0,
     >                        TWO      = 2.E+0,
     >                        HALF     = 0.5+0,
     >                        ZERO     = 0.E+0
      LOGICAL,   PARAMETER :: FALSE    = .FALSE.,
     >                        TRUE     = .TRUE.
      CHARACTER, PARAMETER :: TIGHT*1  = 'M',
     >                        NAME*9   = 'GRADDIS: '


C     Local Variables:
C     ----------------

      INTEGER :: I, LUNERR, NITERS
      REAL    :: ONEOVERN, STOTAL
      REAL    :: WORK(4*NDAT + 5*N)

C     Procedures:
C     -----------

      EXTERNAL :: ARBDIS, LCSFIT, SMOOTH1D_NITERS


C     Execution:
C     ----------

      LUNERR = ABS (LUNOUT)

      IF (POWER < ZERO .OR. POWER > TWO) THEN
         WRITE (LUNERR, '(/, 2A, ES12.4)')
     >      NAME, 'Bad POWER input outside [0., 2.]:', POWER
         IER = 1
         GO TO 999
      END IF

C     Moderate the local gradient magnitudes:

      IF (POWER == HALF) THEN  ! Avoid logarithms
         DO I = 1, NDAT
            WORK(I) = ONE / SQRT (GDAT(I) + ALPHA)
         END DO
      ELSE
         DO I = 1, NDAT
            WORK(I) = (GDAT(I) + ALPHA)**(-POWER)
         END DO
      END IF

CCCC  write (50, '(a)') '# Unsmoothed shape function'
CCCC  write (50, '(2es14.6)') (sdat(i), work(i), i = 1, ndat)

      IF (ISMOOTH == 1 .OR. ISMOOTH == 3) THEN

C        Smooth the shape function further with an explicit method:

         NITERS = MAX (4, MIN (NDAT/40, 15))
CCCC     write (*, '(a, i5)') 'NITERS for shape function:', NITERS

         CALL SMOOTH1D_NITERS (NITERS, 1, NDAT, SDAT, WORK(1))

CCCC     write (51, '(a)') '# smoothed shape function'
CCCC     write (51, '(2es14.6)') (sdat(i), work(i), i = 1, ndat)

      END IF

C     Redistribute the arc lengths based on the shape function:

      STOTAL = SDAT(NDAT)

      CALL ARBDIS (N, ZERO, STOTAL, NDAT, SDAT, WORK(1), 'A',
     >             LUNOUT, WORK(NDAT+1), S, IER)

      IF (IER /= 0) THEN
         WRITE (LUNERR, '(/, 2A, I4)') NAME, 'ARBDIS IER:', IER
         IER = 2
         GO TO 999
      END IF

      IF (ISMOOTH >= 2) THEN  ! Smooth the redistributed arc-lengths

         ONEOVERN = ONE / REAL (N)  ! Index-based abscissas in ARBDIS workspace
         DO I = 1, N
            WORK(NDAT+I) = REAL (I) * ONEOVERN
         END DO

         NITERS = MAX (5, MIN (N/20, 20))
CCCC     write (*, '(a, i5)') 'NITERS for redistributed arc lengths:',
CCCC >      NITERS

CCCC     write (52, '(a)') '# Unsmoothed redistributed arc lengths'
CCCC     write (52, '(2es14.6)') (work(ndat+i), s(i), i = 1, n)

         CALL SMOOTH1D_NITERS (NITERS, 1, N, WORK(NDAT+1), S)

CCCC     write (53, '(a)') '# smoothed redistributed arc lengths'
CCCC     write (53, '(2es14.6)') (work(ndat+i), s(i), i = 1, n)

      END IF

C     Redistribute the curve coordinates:
                                                         ! Unused derivs.
      CALL LCSFIT (NDAT, SDAT, XDAT, TRUE, TIGHT, N, S, X, WORK(NDAT+1))
      CALL LCSFIT (NDAT, SDAT, YDAT, TRUE, TIGHT, N, S, Y, WORK(NDAT+1))

  999 RETURN

      END SUBROUTINE GRADDIS
