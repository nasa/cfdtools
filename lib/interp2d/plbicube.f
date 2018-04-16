C+------------------------------------------------------------------------------
C
      SUBROUTINE PLBICUBE (MODE, IDIM, JDIM, I1, I2, J1, J2, XD, YD, ZD,
     >                     UD, VD, UE, VE, IE, JE, IMAT, JMAT, XYZMAT,
     >                     EPS, P, Q, XE, YE, ZE, XYZU, XYZV, IER)
C
C ACRONYM: Parametric Local BICUBic interpolation
C          -          -     -----
C PURPOSE:
C
C        PLBICUBE performs parametric bicubic spline interpolation of an XYZ
C     surface for the given point (UE,VE) in parameter space.  In the spirit of
C     1-D local spline techniques, it avoids storing the bulk of the spline
C     coefficients - they are generated on the fly once the cell containing the
C     target point is identified.
C
C        The surface data points must form a regular mesh, not necessarily
C     rectangular (although results are more reliable for parallelogram cells).
C     The mesh can be parameterized using chord-length approximations to arc
C     lengths along the grid lines, which are usually normalized.  (Some appli-
C     cations benefit from retaining the arc-length-like nature of unnormalized
C     (u,v)s.  For instance, the singular nose point of an aircraft fuselage
C     presents problems for normalized data.)
C
C        First partial derivatives with respect to u and v are also returned
C     if requested (as would be needed for surface intersection calculations).
C     Efficient application to more than one surface at a time is provided for.
C
C        For the nonparametric bicubic form, see LCSFIT2D.
C
C METHOD:
C
C        The Hermite-type interpolation of Coon's bicubic scheme is used, giving
C     first derivative continuity at cell boundaries.  Since the bicubic form is
C     a Cartesian product of 1-D forms in parametric variables U and V, the
C     method employs the 4-point formulas of the 1-D techniques.  The added
C     requirement of cross derivatives is accommodated by the same 4x4 stencil
C     (possibly 3x* or 2x* at a boundary).  As few as 2 points in each direction
C     are handled, but periodic end conditions are not.
C
C        Providing for interleaved application to more than surface is awkward:
C     local variables cannot be SAVEd for reuse on the next call because they
C     may apply to another surface.  Arguments must be used instead as explained
C     now.
C
C        (1) The application may signal the first call for a new surface
C     (sub)grid by setting UD (I1, J1) = -1., causing PLBICUBE to parameterize
C     it with a call to PARAM2D.  Alternatively, the application may choose to
C     call PARAM2D explicitly.
C
C     NOTE: = -1. was originally < 0., when u, v were assumed in [0,1].
C             -1. was chosen for consistency with existing applications.
C             If -1. is a legitimate value for UD (I1, J1), this later
C             generalization for arbitrary u, v has caused a problem.
C
C        (2) Whether the present target cell for a given surface is the same
C     cell as on the previous call for that surface (i.e., whether the input
C     FMATRIX is reusable) is determined via the arguments IMAT and JMAT.
C     These need not be initialized unless PARAM2D is invoked outside PLBICUBE,
C     when I1 - 1, J1 - 1 are suggested.  They are checked for a match once the
C     current target cell is identified, and updated to the latter's indices
C     prior to returning, ready for the next call for the same surface.
C
C        (3) Note that applications involving more than one surface must be
C     sure to pass surface-specific variables for many of the arguments,
C     including IMAT and JMAT.
C
C        An outline of the procedure for a given surface follows:
C
C     >  IF <new surface, signalled by ud (i1,j1) .eq. -1.> THEN
C
C        >  parameterize it with a call to PARAM2D
C
C     >  Identify cell (i,j) = (ie,je) containing (ue,ve), and corresponding
C        p and q fractional coordinates.
C
C     >  IF <different from cell on previous call for this surface> THEN
C
C        >  Set up 3 matrices (F = X, Y, & Z) of grid-pt.-related information:
C
C                      | F | Fq |                 | F (i,j)   F (i,j+1) |
C           XYZMAT  =  |---+----|  where | F | =  |                     |
C                      | Fp| Fpq|                 |F(i+1,j)   F(i+1,j+1)|
C
C           and similarly for the first & cross derivatives w.r.t. p & q.
C
C           Note that these derivatives are not Fu, Fv, Fuv - a source of
C           confusion during the initial implementation.
C
C     >  Set up products U N and N^ V where U = [p^3 p^2 p 1] and N is
C
C           | 2 -2  1  1|
C           |-3  3 -2 -1|  (from a rearrangement of the Hermite formula
C           | 0  0  1  0|   for the cubic with specified f and f' at two
C           | 1  0  0  0|   end points)
C
C        and V is a similar function of q.  (N^ is N-transpose here.)
C
C     >  Then  F (u,v) = U N XYZMAT N^ V  and (if requested)
C              Fu      = U'N XYZMAT N^ V  (derivative w.r.t. u)
C              Fv      = U N XYZMAT N^ V' (derivative w.r.t. v)
C        for each F = X, Y, & Z.
C
C HISTORY:
C
C   11/22/93    DAS   Initial development in conjunction with LCSFIT2D.
C   08/15/94     "    RIPPLE2D now has arguments P, Q (not used here).
C   03/20/95     "    Call PARAM2D if UD (I1, J1) = -1., not < 0.;
C                     appended SDOT4 here in place of BLAS SDOT utility.
C   01/23/96     "    RIPPLE2D now returns usable results if it fails.
C   11/17/97     "    The P and Q from RIPPLE2D are what should be used,
C                     not (UE - UD(IT,JT)) / (UD(IT+1,JT) - UD(IT,JT)), etc.
C                     Test data with ~rectangular (u,v)s obscured the error.
C   06/23/99     "    BICUBMAT was making a simplifying assumption which
C                     was exact only for uniform (u,v) spacing.  Somehow,
C                     resulting function values still tended to be on the
C                     surface, but cases with highly non-uniform spacing
C                     showed that evaluations at uniform intervals within
C                     cells gave values that were shifted along the surface.
C                     The only changes to PLBICUBE are in the header, except
C                     that SDOT4 is now a Fortran 90 internal procedure.
C   07/22/99     "    The derivatives had been correct only for rectangular
C                     cells.  Now we use the chain rule appropriately for
C                     arbitrary cell shapes, although parallelograms are best.
C   08/06/99     "    PLBICUBE now returns (p,q), needed for PLBICUT to return.
C
C AUTHOR: David Saunders, Sterling Software/NASA Ames, Moffett Field, CA.
C
C ------------------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT  (IN) ::
     >   MODE                  ! 0 means return a function value only;
                               ! 1 means return first partial derivatives too
      INTEGER, INTENT (IN) ::
     >   IDIM, JDIM,           ! Max. # pts. provided for in the surface arrays
     >   I1, I2, J1, J2        ! Grid index range eligible for searching;
                               ! 1 <= I1 < I1 + 1 < I2 <= IDIM, etc.

      REAL, INTENT (IN), DIMENSION (IDIM, JDIM) ::
     >   XD, YD, ZD            ! Surface grid coordinates

      REAL, INTENT (INOUT), DIMENSION (IDIM, JDIM) ::
     >   UD, VD                ! Parametric variables at the surface grid pts.:
                               ! calculated internally if requested via
                               ! UD (I1, J1) < 0 on input, else input from a
                               ! prior call to PLBICUBE or PARAM2D
      REAL, INTENT (IN) ::
     >   UE, VE                ! Target (u,v) point (E = evaluation)

      INTEGER, INTENT (INOUT) ::
     >   IE, JE                ! Indices of the surface cell containing the
                               ! target point:
                               ! on input:  the values from a previous
                               !            call, or I1, J1, if no better
                               !            estimate is known;
                               ! on output: the "lower left" indices of
                               !            the cell (unless IER > 0);
                               ! I1 <= IE < I2 input & output; likewise for JE
      INTEGER, INTENT (INOUT) ::
     >   IMAT, JMAT            ! Surface-specific target cell indices on a
                               ! previous ! call, employed for efficiency.
                               ! Initialize as explained above if PARAM2D is
                               ! called EXPLICITLY.  DO NOT CHANGE THEIR VALUES
                               ! THEREAFTER either way.
      REAL, INTENT (INOUT) ::
     >   XYZMAT (4, 4, 3)      ! Cell-specific intermediate values: a 4x4 matrix
                               ! for each of X, Y, & Z, passed in and out for
                               ! efficiency reasons as explained above
      REAL, INTENT (IN) ::
     >   EPS                   ! Tolerance for the BILINT search utility;
                               ! try 5.*machine eps for 32-bit arithmetic,
                               ! or 1.E-7 should be OK for 64-bit "
      REAL, INTENT (OUT) ::
     >   P, Q                  ! Fractional cell coords. (p,q) <-> (UE,VE)

      REAL, INTENT (OUT) ::
     >   XE, YE, ZE,           ! Interpolated surface coordinates <-> (UE,VE)
     >   XYZU (3),             ! Partial derivs. of X, Y, & Z w.r.t. U and ...
     >   XYZV (3)              ! ... V at (UE,VE) (if MODE = 1); likewise for

      INTEGER, INTENT (OUT) ::
     >   IER                   ! 0 if the interpolation was successful;
                               ! 1 if the target point was out of range;
                               ! results may still be usable near a boundary
C     Procedures:

      EXTERNAL
     >   BICUBMAT,             ! Ancillary, common to LCSFIT2D and PLBICUBE
     >   PARAM2D,              ! Parameterizes a surface grid (chord lengths)
     >   RIPPLE2D              ! 2D grid search utility

C-------------------------------------------------------------------------------

C     Local constants:

      REAL, PARAMETER ::
     >   ONE = 1.0, TWO = 2.0, THREE = 3.0, FOUR = 4.0, SIX = 6.0

C     Local variables:

      INTEGER
     >   I, IT, J, JT

      REAL
     >   A (2, 2), AJ (2, 2), DP (2), FP, FQ, FPQ, SQ,
     >   UN (4), UNF (4, 3), UPRIMEN (4), VN (4), VPRIMEN (4), XYZE (3)

C     Execution:

C     Parameterize the (sub)grid?
C     ---------------------------

      IF (UD (I1, J1) == -ONE) THEN ! Originally < 0, but now u, v
                                    ! are not necessarily in [0,1].

         CALL PARAM2D (IDIM, JDIM, I1, I2, J1, J2, XD, YD, ZD, UD, VD)

         IMAT = I1 - 1 ! So they are defined when checked below
         JMAT = J1 - 1
      END IF


C     Identify the cell containing the target point:
C     ----------------------------------------------

      CALL RIPPLE2D (IDIM, JDIM, I1, I2, J1, J2, UD, VD, UE, VE, IE, JE,
     >               EPS, P, Q, IER)

C****IF (IER /= 0) GO TO 99 ! Search failed, but RIPPLE2D now gives
                            ! the best results that it found

      IT = IE ! To avoid excessive argument references
      JT = JE


C     New cell-specific matrices for X, Y, Z and their derivatives?
C     -------------------------------------------------------------

      IF (IT /= IMAT .OR. JT /= JMAT) THEN

         IMAT = IT ! N.B.: These are surface-specific
         JMAT = JT

         CALL BICUBMAT (IDIM, JDIM, I1, I2, J1, J2, UD, VD, XD, IT, JT,
     >                  XYZMAT (1, 1, 1))
         CALL BICUBMAT (IDIM, JDIM, I1, I2, J1, J2, UD, VD, YD, IT, JT,
     >                  XYZMAT (1, 1, 2))
         CALL BICUBMAT (IDIM, JDIM, I1, I2, J1, J2, UD, VD, ZD, IT, JT,
     >                  XYZMAT (1, 1, 3))
      END IF


C     The rest depends on the target values of the parametric variables:
C     ------------------------------------------------------------------

C     Set up product [U N] where U = [p^3 p^2 p 1] and N is
C
C           | 2 -2  1  1|
C           |-3  3 -2 -1|  (from a rearrangement of the Hermite formula
C           | 0  0  1  0|   for the cubic with specified f and f' at two
C           | 1  0  0  0|   end points)

      SQ = P * P
      UN (2) = SQ * (THREE - TWO * P)
      UN (1) = ONE - UN (2)
      UN (3) = SQ * (P - TWO) + P
      UN (4) = SQ * (P - ONE)

C     Similarly for [N^ V]:

      SQ = Q * Q
      VN (2) = SQ * (THREE - TWO * Q)
      VN (1) = ONE - VN (2)
      VN (3) = SQ * (Q - TWO) + Q
      VN (4) = SQ * (Q - ONE)

C     For each of X, Y, & Z:
C     ----------------------

      DO J = 1, 3

C        [U N] [XYZMAT] is a row vector:

         DO I = 1, 4
            UNF (I, J) = SDOT4 (UN, XYZMAT (1, I, J))
         END DO

C        [U N] [XYZMAT] [N^ V] is a scalar (the desired interpolation):

         XYZE (J) = SDOT4 (UNF (1, J), VN)
C        --------

      END DO

      XE = XYZE (1)
      YE = XYZE (2)
      ZE = XYZE (3)


C     Interpolated derivatives as well?
C     ---------------------------------

      IF (MODE /= 0) THEN

C        The chain rule converts p/q derivatives to u/v derivatives.
C        Set up the 2 x 2 LHS matrix common to derivatives of x, y, and z:

         FP  = UD (IT + 1, JT) - UD (IT, JT)
         FQ  = UD (IT, JT + 1) - UD (IT, JT)
         FPQ = UD (IT + 1, JT + 1) - UD (IT, JT + 1) - FP

         AJ (1, 1) = FP + Q * FPQ ! du/dp
         AJ (2, 1) = FQ + P * FPQ ! du/dq

         FP  = VD (IT + 1, JT) - VD (IT, JT)
         FQ  = VD (IT, JT + 1) - VD (IT, JT)
         FPQ = VD (IT + 1, JT + 1) - VD (IT, JT + 1) - FP

         AJ (1, 2) = FP + Q * FPQ ! dv/dp
         AJ (2, 2) = FQ + P * FPQ ! dv/dq

C        U' is [3p^2, 2p, 1, 0] (derivative w.r.t. p, not u).  Form [U'N]:

         UPRIMEN (1) = SIX * P * (P - ONE)
         UPRIMEN (2) = -UPRIMEN (1)
         UPRIMEN (3) = P * (P * THREE - FOUR) + ONE
         UPRIMEN (4) = P * (P * THREE - TWO)

C        Similarly for [N^ V']:

         VPRIMEN (1) = SIX * Q * (Q - ONE)
         VPRIMEN (2) = -VPRIMEN (1)
         VPRIMEN (3) = Q * (Q * THREE - FOUR) + ONE
         VPRIMEN (4) = Q * (Q * THREE - TWO)

C        For each of X, Y, & Z:
C        ----------------------

         DO J = 1, 3

C           [U N XYZMAT] [N^ V'] is the derivative w.r.t. q.

            DP (2) = SDOT4 (UNF (1, J), VPRIMEN) ! d/dq

C           [U'N XYZMAT] [N^ V]  (derivative w.r.t. p):

            DO I = 1, 4
               UNF (I, J) = SDOT4 (UPRIMEN, XYZMAT (1, I, J))
            END DO

            DP (1) = SDOT4 (UNF (1, J), VN) ! d/dp

            A = AJ ! Could factorize once & solve 3 times, but LUSOLVE
                   ! is in use elsewhere, while DECOMP & SOLVE are not.

            CALL LUSOLVE (2, 2, A, DP, IER) ! Solve J dfuv = dfpq

            XYZU (J) = DP (1) ! d/du
            XYZV (J) = DP (2) ! d/dv

         END DO

      END IF

C     Internal procedure for PLBICUBE:

      CONTAINS

!        ---------------------------------------------------------------
         REAL FUNCTION SDOT4 (SX, SY)
!
!        Dot product of two 4-vectors, to avoid single/double problems
!        with the BLAS utility SDOT originally used by PLBICUBE.
!        ---------------------------------------------------------------

         REAL SX (4), SY (4)

         SDOT4 = SX (1) * SY (1) + SX (2) * SY (2) +
     >           SX (3) * SY (3) + SX (4) * SY (4)

         END FUNCTION SDOT4

      END SUBROUTINE PLBICUBE
