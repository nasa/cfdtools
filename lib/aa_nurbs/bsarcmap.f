C+------------------------------------------------------------------------------
C
      SUBROUTINE BSARCMAP (C, T1, T2, NPTS, TNORM, ARC, T, LUNERR, IER)
C
C ONE-LINER: B-Spline curve ARC-length-based MAPping of T
C            - -            ---              ---
C PURPOSE:
C
C        BSARCMAP provides an approximate solution to the problem of
C     distributing points on a B-spline curve such that their distances
C     along the curve vary in some desired way.  (Distributing the
C     parametric variable T directly may not be good enough, since T
C     does not necessarily vary uniformly with arc length.)
C
C        This version is limited to 2- or 3-space curves for work-space
C     reasons.
C
C METHOD:
C
C        The mapping is numerical, using cumulative chord length to
C     approximate arc length, and using local spline interpolation
C     techniques to determine the output distribution of T from the
C     correspondence between the (denormalized) input distribution and
C     the associated arc-length-like distances calculated here.  As long
C     as the curve is not grossly nonuniform and TNORM is not too coarse,
C     a good approximation to the desired arc-length variation should be
C     obtained.
C
C ENVIRONMENT:
C     VAX/VMS; FORTRAN 77, with...
C     > IMPLICIT NONE
C     > Trailing ! comments
C     > Names up to 8 characters
C
C HISTORY:
C     04/26/93  D.A.Saunders  Initial implementation.
C     02/01/94       "        Dispensed with DOUBLE PRECISION.  Compile with
C                             /REAL_LENGTH=64 or -r8 switches as necessary.
C     07/16/97       "        Normalizing ARC before splining is better than
C                             denormalizing TNORM (I) for each interpolation.
C
C AUTHOR: David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
C
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments
C     ---------

      REAL    C (*)     ! (I)   DT_NURBS B-spline curve vector:
                        !       C(1) = 1 (# parametric variables);
                        !       C(2) = 2 or 3 (for X and Y [and Z]);
                        !       C(3) = polynomial order k = degree + 1;
                        !       C(4) = # control points, N;
                        !       C(5) = pointer for efficient evaluations;
                        !       C(6 : 5 + N + k) = knot vector,
                        !       followed by X coords. of control pts.,
                        !       ........... Y .......
                        !       [.......... Z .......]

      REAL    T1, T2    ! (I)   T values corresponding to the ends
                        !       of the interval being distributed.  T1 < T2

      INTEGER NPTS      ! (I)   Desired number of points.

      REAL    TNORM (NPTS)      ! (I)   Normalized form of desired distribution,
                                !       on [0, 1]:
                                !       TNORM (1) <-> T1 and TNORM (NPTS) <-> T2

      REAL    ARC (NPTS)! (S)   Work-space for chord-length values.

      REAL    T (NPTS)  ! (O)   Distribution of T which provides an
                        !       approximation to the desired distribution
                        !       of arc length.  T (1 & NPTS) = T1 & T2.
                           
      INTEGER LUNERR    ! (I)   Logical unit for error messages.

      INTEGER IER       ! (O)   Success/error code:
                        !        0 means no problem was encountered;
                        !       12 means <to be completed>

C     Procedures
C     ----------

      EXTERNAL DTSPVL   ! Evaluates a B-spline curve at the given T.
      EXTERNAL LCSFIT   ! Local cubic spline interpolation in 2-space.

C-------------------------------------------------------------------------------

C     Local constants
C     ---------------

      INTEGER
     >   MAXDIM, MAXK, NWORK
      REAL
     >   ONE, ZERO
      CHARACTER
     >   METHOD * 1

      PARAMETER
     >  (MAXDIM = 3,              ! Arc lengths don't mean much in 4-space ...
     >   MAXK   = 11,             ! Local storage handles polynomial deg. <= 10
     >   METHOD = 'B',            ! "Bessel" fit (no need for monotonic)
     >   NWORK  = 5 * MAXK - 2,   ! This is why the degree is limited
     >   ONE    = 1.,
     >   ZERO   = 0.)

C     Local variables
C     ---------------

      INTEGER
     >   I, J, NDIM
      REAL
     >   RTOTAL, DT, SUM, TI, WORK (NWORK), XY (MAXDIM), XY0 (MAXDIM)

C     Execution
C     ---------

C     Most of the error checking is left to the lower level utilities.

      NDIM = NINT (C (2))

      IF (NDIM .GT. MAXDIM) THEN
         IER = 2
         WRITE (LUNERR, 1000) 'Too many dimensions.'
         GO TO 900
      END IF

      IF (NINT (C (3)) .GT. MAXK) THEN
         IER = 3
         WRITE (LUNERR, 1000) 'Too high a degree.'
         GO TO 900
      END IF

C     Evaluate the initial end point:

      CALL DTSPVL (T1, C, WORK, NWORK, XY0, IER)

      IF (IER .NE. 0) THEN
         WRITE (LUNERR, 1000) 'DTSPVL failure. I and T:', 1, T1
         GO TO 900
      END IF

      DT = T2 - T1
      ARC (1) = ZERO

C     For each further point in the normalized distribution ...

      DO I = 2, NPTS

C        Evaluate the curve at the denormalized T:

         TI = T1 + DT * TNORM (I)
         IF (I .EQ. NPTS) TI = T2  ! Make sure of it

         CALL DTSPVL (TI, C, WORK, NWORK, XY, IER)

         IF (IER .NE. 0) THEN
            WRITE (LUNERR, 1000) 'DTSPVL failure. I and T:', I, TI
            GO TO 900
         END IF

C        Find the approximate arc-length increment:

         SUM = ZERO
         DO J = 1, NDIM
            SUM = SUM + (XY (J) - XY0 (J)) ** 2
            XY0 (J) = XY (J)       ! Save the previous coordinates
         END DO

         ARC (I) = ARC (I - 1) + SQRT (SUM)
      END DO

C     Normalize the evaluated arc lengths:

      RTOTAL = ONE / ARC (NPTS)

      DO I = 2, NPTS
         ARC (I) = ARC (I) * RTOTAL
      END DO
      ARC (NPTS) = ONE

C     For each point of the desired normalized arc-length distribution,
C     interpolate a T within TNORM (*) then denormalize it:

      CALL LCSFIT (NPTS, ARC, TNORM, .TRUE., METHOD, NPTS - 2,
     >             TNORM (2), T (2), T (2))

      T (1) = T1
      DO I = 2, NPTS - 1
         T (I) = T1 + DT * T (I)
      END DO
      T (NPTS) = T2

      GO TO 999


C     Error handling:

  900 WRITE (LUNERR, 1010) IER


C     Termination:

  999 RETURN

 
C     Formats:

 1000 FORMAT (/, ' *** BSARCMAP trouble: ', A, I6)
 1010 FORMAT ('  IER: ', I4)

      END
