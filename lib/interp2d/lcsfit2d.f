C+----------------------------------------------------------------------
C
      SUBROUTINE LCSFIT2D (MODE, IDIM, JDIM, I1, I2, J1, J2, XD, YD, FD,
     >                     XEVAL, YEVAL, IEVAL, JEVAL, IMAT, JMAT,
     >                     FMATRIX, EPS, P, Q, FEVAL, FXEVAL, FYEVAL,
     >                     IER)
C
C ACRONYM: Local Cubic Spline FIT in 2 Dimensions
C          -     -     -      ---    - -
C PURPOSE:
C
C        LCSFIT2D performs bicubic local spline interpolation of a
C     discretized function of X and Y at the given point (XEVAL, YEVAL).
C     In the spirit of 1-D local spline techniques, it avoids storing
C     the bulk of the spline coefficients - they are generated on the
C     fly once the cell containing the target point is identified.
C
C        The data points are expected to form a regular grid, not
C     necessarily uniform or rectangular, although results are most
C     reliable for parallelogram cells.  F should be a single-valued
C     function of X and Y.  Otherwise, use the parametric bicubic form,
C     PLBICUBE.
C
C        First partial derivatives are also returned if requested.
C     Efficient application to more than one function at a time is
C     provided for.
C
C METHOD:
C
C        The Hermite-type interpolation of Coon's bicubic scheme is
C     used, giving first derivative continuity at cell boundaries.
C     Since the bicubic form is basically a Cartesian product of 1-D
C     forms in variables X and Y, the method employs the 4-point
C     formulas of the 1-D techniques.  The added requirement of cross
C     derivatives is accommodated by the same 4x4 stencil (possibly
C     3x* or 2x* at a boundary).  As few as 2 points in each direction
C     are handled.
C
C        Providing for interleaved application to more than one dataset
C     is awkward: local variables cannot be SAVEd for reuse on the next
C     call because they may apply to another dataset.  Arguments must
C     be used instead, as explained now.
C
C        (1) Whether the present target cell for a given dataset is the
C     same cell as on the previous call for that dataset (i.e., whether
C     the input FMATRIX is reusable) is determined via the arguments
C     IMAT and JMAT.  On the first call for the surface, they should
C     be input with out-of-range values.  They are checked for a match
C     when the current target cell is identified, and updated to the
C     latter's indices prior to returning, ready for the next call.
C     DO NOT CHANGE THEIR VALUES BETWEEN CALLS.
C
C        (2) Note that applications involving more than one dataset must
C     be careful to pass dataset-specific variables for numerous arguments,
C     including IMAT and JMAT.
C
C        An outline of the procedure for a given dataset follows:
C
C     >  Identify the cell (i,j) = (ieval,jeval) containing (xeval,yeval)
C        with corresponding (p,q) fractional coordinates.
C
C     >  IF <different from cell on previous call for this dataset> THEN
C
C        >  Set up the matrix of grid-point-related information:
C
C                       | F | Fq |               | F (i,j)   F (i,j+1) |
C           FMATRIX  =  |---+----|  where  F  =  |                     |
C                       | Fp| Fpq|               |F(i+1,j)   F(i+1,j+1)|
C
C           and similarly for the first & cross derivatives w.r.t. p & q.
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
C     >  Then  F (x,y) = U N FMATRIX N^ V  and (if requested)
C              Fp      = U'N FMATRIX N^ V  (derivative w.r.t. p)
C              Fq      = U N FMATRIX N^ V' (derivative w.r.t. q)
C
C     >  Convert Fp, Fq to Fx, Fy via the chain rule (solve a 2 x 2 system).
C
C HISTORY:
C
C     11/22/93  DAS  Initial development in conjunction with PLBICUBE.
C     08/15/94   "   RIPPLE2D now has arguments P, Q (not used here).
C     10/12/95   "   Clarified applicability, which is more limited than
C                    originally realized. PLBICUBE is usually preferable.
C     01/23/96   "   RIPPLE2D failure still gives possibly usable results.
C     11/17/97   "   The P and Q from RIPPLE2D are what should be used,
C                    not P = (XEVAL - XD(IT,JT))/(XD(IT+1,JT) - XD(IT,JT)).
C     07/21/99   "   Derivatives were correct only for rectangular cells.
C                    Solve a 2 x 2 system to handle cell skewing.
C     08/06/99   "   Return (p,q) for consistency with analogous routines.
C
C AUTHOR: David Saunders, Sterling Software/NASA Ames, Moffett Field, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments:

      INTEGER, INTENT (IN) ::
     >   MODE,            ! MODE = 0 means return a function value only;
                          !      = 1 means return first partial derivatives too
     >   IDIM,            ! Max. # points provided for in the data arrays.
     >   JDIM,            ! (JDIM is superfluous but retained for symmetry.)
     >   I1,              ! First & last points in each grid direction
     >   I2,              ! eligible for the search and interpolation:
     >   J1,              !    1 <= I1 < I1 + 1 < I2 <= IDIM, etc.
     >   J2               ! Be careful not to include singular grid cells.

      REAL, INTENT (IN), DIMENSION (IDIM, JDIM) ::
     >   XD, YD,          ! (I1:I2, J1:J2) are the meaningful data coordinates
     >   FD               ! FD may be the third coordinate of a surface or it
                          ! it may be some other function of X and Y (D = Data)
      REAL, INTENT (IN) ::
     >   XEVAL, YEVAL     ! Coordinates of the target point

      INTEGER, INTENT (INOUT) ::
     >   IEVAL, JEVAL     ! Indices of the cell containing the target point;
                          ! on input:  the values from a previous call, or
                          !            I1, J1 if no better estimate is known;
                          ! on output: the "lower left" indices of the cell
                          !            containing the target;
                          ! I1 <= IEVAL < I2 on input & output; sim. for JEVAL

      INTEGER, INTENT (INOUT) ::
     >   IMAT, JMAT       ! Dataset-specific target cell indices on a previous
                          ! call, employed for efficiency.  Initialize as
                          ! explained above; DO NOT CHANGE THEREAFTER.

      REAL, INTENT (INOUT) ::
     >   FMATRIX (4, 4)   ! Cell-specific intermediate values isolated for
                          ! efficiency reasons as explained above

      REAL, INTENT (IN) ::
     >   EPS              ! Tolerance used by the search utility - see BILINT

      REAL, INTENT (OUT) ::
     >   P, Q,            ! Fractional cell coords. (p,q) <-> (XEVAL, YEVAL)
     >   FEVAL,           ! Interpolated function value    "      "     "
     >   FXEVAL,          ! Partial df/dx and df/dy (if MODE = 1) "     "
     >   FYEVAL

      INTEGER, INTENT (OUT) ::
     >   IER              ! IER = 0 if the interpolation was successful;
                          !     = 1 if the target point was out of range,
                          !         but results may still be adequate
                     
C     Procedures:

      EXTERNAL
     >   BICUBMAT,        ! Ancillary routine common to LCSFIT2D and PLBICUBE
     >   LUSOLVE,         ! For 2 x 2 solution if derivatives are requested
     >   RIPPLE2D         ! 2D grid search utility

C     Local constants:

      REAL, PARAMETER ::
     >   ONE = 1., TWO = 2., THREE = 3., FOUR = 4., SIX = 6.

C     Local variables:

      INTEGER
     >   I, IT, JT

      REAL
     >   AJ (2, 2), DP (2), FP, FQ, FPQ, SQ

      REAL, DIMENSION (4) ::
     >   UN, UNF, UPRIMEN, VN, VPRIMEN


C     Execution:

C     Identify the cell containing the target point:
C     ----------------------------------------------

      CALL RIPPLE2D (IDIM, JDIM, I1, I2, J1, J2, XD, YD, XEVAL, YEVAL,
     >               IEVAL, JEVAL, EPS, P, Q, IER)

C*****IF (IER /= 0) GO TO 99   ! Search failed, but results may be usable

      IT = IEVAL    ! To avoid excessive argument references
      JT = JEVAL


C     New cell-specific matrix of F, Fx, Fy, and Fxy elements?
C     --------------------------------------------------------

      IF (IT /= IMAT .OR. JT /= JMAT) THEN

         IMAT = IT   ! N.B.: These are dataset-specific
         JMAT = JT

         CALL BICUBMAT (IDIM, JDIM, I1, I2, J1, J2, XD, YD, FD, IT, JT,
     >                  FMATRIX)
      END IF


C     The rest depends on the target (X, Y).
C     --------------------------------------

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

C     [U N] [FMATRIX] is a row vector:

      DO I = 1, 4
         UNF (I) = SDOT4 (UN, FMATRIX (1, I))
      END DO

C     [U N] [FMATRIX] [N^ V] is a scalar (the desired interpolation):

      FEVAL = SDOT4 (UNF, VN)
C     -----

C     Interpolated derivatives as well?
C     ---------------------------------

      IF (MODE /= 0) THEN

C        U' is [3p^2, 2p, 1, 0] (derivative w.r.t. p, not x).  Form [U'N]:

         UPRIMEN (1) = SIX * P * (P - ONE)
         UPRIMEN (2) = -UPRIMEN (1)
         UPRIMEN (3) = P * (P * THREE - FOUR) + ONE
         UPRIMEN (4) = P * (P * THREE - TWO)

C        Similarly for [N^ V']:

         VPRIMEN (1) = SIX * Q * (Q - ONE)
         VPRIMEN (2) = -VPRIMEN (1)
         VPRIMEN (3) = Q * (Q * THREE - FOUR) + ONE
         VPRIMEN (4) = Q * (Q * THREE - TWO)

C        [U N FMATRIX] [N^ V'] is the derivative w.r.t. q.

         DP (2) = SDOT4 (UNF, VPRIMEN) ! dF/dq

C        [U'N FMATRIX] [N^ V] is the derivative w.r.t. p:

         DO I = 1, 4
            UNF (I) = SDOT4 (UPRIMEN, FMATRIX (1, I))
         END DO

         DP (1) = SDOT4 (UNF, VN) ! dF/dp

C        Use the chain rule to convert to derivatives w.r.t. x and y:

         FP  = XD (IT + 1, JT) - XD (IT, JT)
         FQ  = XD (IT, JT + 1) - XD (IT, JT)
         FPQ = XD (IT + 1, JT + 1) - XD (IT, JT + 1) - FP

         AJ (1, 1) = FP + Q * FPQ ! dx/dp
         AJ (2, 1) = FQ + P * FPQ ! dx/dq

         FP  = YD (IT + 1, JT) - YD (IT, JT)
         FQ  = YD (IT, JT + 1) - YD (IT, JT)
         FPQ = YD (IT + 1, JT + 1) - YD (IT, JT + 1) - FP

         AJ (1, 2) = FP + Q * FPQ ! dy/dp
         AJ (2, 2) = FQ + P * FPQ ! dy/dq

         CALL LUSOLVE (2, 2, AJ, DP, IER) ! Solve J dfxy = dfpq

         FXEVAL = DP (1)          ! df/dx
         FYEVAL = DP (2)          ! df/dy
 
      END IF

C     Internal procedure for LCSFIT2D:

      CONTAINS

!        ---------------------------------------------------------------
         REAL FUNCTION SDOT4 (SX, SY)
!
!        Dot product of two 4-vectors, to avoid single/double problems
!        with the BLAS utility SDOT originally used.
!        ---------------------------------------------------------------

         REAL SX (4), SY (4)

         SDOT4 = SX (1) * SY (1) + SX (2) * SY (2) +
     >           SX (3) * SY (3) + SX (4) * SY (4)

         END FUNCTION SDOT4

      END SUBROUTINE LCSFIT2D
