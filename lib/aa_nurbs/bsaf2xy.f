C+------------------------------------------------------------------------------
C
      SUBROUTINE BSAF2XY (C, NU, NL, T, X, Y, LUNERR, IER)
C
C ONE-LINER: B-Spline AirFoil conversion to (X, Y) points
C            - -      -  -               -   -  -
C PURPOSE:
C
C        BSAF2XY discretizes an airfoil given in DT_NURBS B-spline form.
C     It is specialized to airfoils because it assumes 2-space and needs
C     to identify the leading edge (a minimization-in-T problem) then
C     use appropriate distributions of T on each surface.
C
C        BSAF2XY is intended for compatibility with applications which
C     expect discrete (X,Y) airfoil definitions, or for plotting.  It
C     discretizes the parametric variable rather than arc length.  See
C     BSARCMAP as used by program BSPROFILE if precise distributions of
C     arc length are essential.  BSAF2XY is cheaper than using BSARCMAP,
C     and it should easily suffice for the intended purposes.
C
C        See also BSAFY4X for discretization at specified Xs.
C
C METHOD:
C
C        The curve is assumed to represent the airfoil in clockwise
C     wrap-around form, from lower trailing edge to upper trailing edge.
C     General purpose utility BSZEROMN is used to locate the value T*
C     of the parametric variable T corresponding to the leading edge.
C     Then the intervals [T1, T*] and [T*, TN] are distributed via the
C     polynomial/sinusoidal scheme of FOILGRD using what is normally a
C     good choice of inputs for airfoils, where T1 and TN are the
C     first and last knots in the knot vector contained in C (*).
C
C ENVIRONMENT:
C     DEC/OpenVMS; SGI/IRIX; FORTRAN 77, with...
C     > IMPLICIT NONE
C     > Trailing ! comments
C     > Names up to 8 characters
C
C HISTORY:
C     05/05/92  D.A.Saunders  Initial design and implementation, with
C                             application to existing flow solvers in
C                             mind:  we don't want to rethink their
C                             grid generation at this stage, so we need
C                             to discretize the B-spline.  (Vinokur
C                             distributions used.)
C     02/01/94       "        Dispensed with DOUBLE PRECISION.  Compile with
C                             /REAL_LENGTH=64 or -r8 switches as necessary.
C     07/21/97       "        Switched from Vinokur distributions to the
C                             more convenient distributions of FOILGRD.
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

      INTEGER NU, NL    ! (I)   Number of upper and lower surface
                        !       points to generate.

      REAL    T (*)     ! (S)   Work-space for the T values corresp. to
                        !       X (*) and Y (*) (NL + NU - 1 of them).

      REAL    X (*),    ! (O)   The desired coordinates in wrap-around
     >        Y (*)     !       order as for the B-spline curve.
                        !       Elements 1 : NL + NU - 1 are returned.

      INTEGER LUNERR    ! (I)   Logical unit for error messages.

      INTEGER IER       ! (O)   Success/error code:
                        !        0 means no problem was encountered;
                        !       12 means <to be completed>

C     Procedures
C     ----------

      EXTERNAL BSZEROMN ! Finds a zero, min, or max of a B-spline curve.
      EXTERNAL DTSPVL   ! Evaluates a B-spline curve at the given T.
      EXTERNAL FOILGRD  ! Linear/quadratic/sine/cosine hybrid distribution

C-------------------------------------------------------------------------------

C     Local constants
C     ---------------

      INTEGER
     >   MAXK, NWORK
      REAL
     >   C1, C2, C3, C4
      PARAMETER
     >  (MAXK  = 11,            ! Local storage handles deg. <= 10
     >   NWORK = 5 * MAXK - 2,  ! DTSPVL warns if this is exceeded
     >   C1    = 0.04,          ! FOILGRD coefs. recommended for airfoils
     >   C2    = 0.0,           !    when applied to arc lengths, and
     >   C3    = 0.3,           !    normally adequate when applied to
     >   C4    = 0.66)          !    the T ranges of an airfoil curve

C     Local variables
C     ---------------

      INTEGER
     >   I, IX, IY, NCPTS, NDEG
      REAL
     >   T1, T2, TLE, WORK (NWORK), XY (2)

C     Execution
C     ---------

      NDEG  = NINT (C (3)) - 1
      NCPTS = NINT (C (4))
      IX = 6 + NDEG + NCPTS  ! Offset of X coordinates of control pts.
      IY = IX + NCPTS        !   "   "   Y
      T2 = C (IX)            ! Last knot
      T1 = C (6)             ! First knot

C     Find the leading edge value of T by minimizing X w.r.t. T:

      CALL BSZEROMN ('MIN', 1, C, T1, T2, TLE, XY, LUNERR, IER)
      IF (IER .NE. 0) GO TO 800

C     Generate a distribution of Ts for the lower surface.
C     -NL tells FOILGRD to switch the LE to TE order:

      CALL FOILGRD (-NL, T1, TLE, C1, C2, C3, C4, T (1))

C     ... and for the upper surface:

      CALL FOILGRD ( NU, TLE, T2, C1, C2, C3, C4, T (NL))

C     Evaluate X and Y at each T:

      X (1) = C (1 + IX)  ! Avoid possible round-off at the end points
      Y (1) = C (1 + IY)

      DO I = 2, NL + NU - 2

         CALL DTSPVL (T (I), C, WORK, NWORK, XY, IER)
         IF (IER. NE. 0) GO TO 810

         X (I) = XY (1)
         Y (I) = XY (2)
      END DO

      I = NL + NU - 1
      X (I) = C (NCPTS + IX)
      Y (I) = C (NCPTS + IY)

      GO TO 999


C     Error handling.

  800 WRITE (LUNERR, 1010) 'BSZEROMN', IER
      GO TO 999
      
  810 WRITE (LUNERR, 1010) 'DTSPVL', IER


  999 RETURN

C     Formats.

 1010 FORMAT ('0BSAF2XY: Bad return from ', A, '.  IER: ', I4)

      END
