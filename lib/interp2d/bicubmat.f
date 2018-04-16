C+----------------------------------------------------------------------
C
      SUBROUTINE BICUBMAT (IDIM, JDIM, I1, I2, J1, J2, X, Y, F,
     >                     IEVAL, JEVAL, FMATRIX)
C
C ACRONYM:  BICUBic interpolation MATrix setup
C           -----                 ---
C PURPOSE:
C
C        BICUBMAT sets up the 4x4 matrix of values of a function of two
C     variables and its first and cross derivatives at the four vertices
C     of the grid cell defined by IEVAL and JEVAL, as needed for Hermite-
C     type cubic interpolation in two dimensions on a regular grid.  For
C     the parametric case (PLBICUBE), it is called three times for the
C     target cell, namely for each of X, Y, and Z as functions of U and V.
C     For the nonparametric case (LCSFIT2D), it is called once for F as a
C     function of X and Y.
C
C        Some of the code was not worth trying to avoid repeating here for
C     the parametric case.  Note that the 3-point derivatives are exact for
C     (bi)quadratic data but not for bicubic data, even though we're
C     constructing a bicubic spline.
C
C ARGUMENTS:  These are obvious from a look at PLBICUBE or LCSFIT2D.
C
C HISTORY:
C     11/28/93  DAS  Modularized to suit parametric and nonparametric cases.
C     02/21/93   "   Mix-up found in all cross-derivative calls to DERIV3.
C                    Baffling how it didn't affect the function values
C                    during original testing.
C     06/23/99   "   Highly nonuniform data revealed erroneous simplification:
C                    dX(i,j) and dY(i,j) should be applied at ALL vertices
C                    to convert dF/dX to dF/dp, etc. (not dX(i+1,j), etc.).
C     07/31/99   "   Non-parallelogram cells revealed the need for proper
C                    use of 2 dX and 2 dY values at the target cell edges.
C
C AUTHOR: David Saunders, Sterling Software/NASA Ames, Moffett Field, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments (obvious from a look at LCSFIT2D & PLBICUBE):

      INTEGER, INTENT (IN) ::
     >   IDIM, JDIM, I1, I2, J1, J2, IEVAL, JEVAL

      REAL, INTENT (IN), DIMENSION (IDIM, JDIM) ::
     >   X, Y, F

      REAL, INTENT (OUT) ::
     >   FMATRIX (4, 4)

C     Local constants:

      REAL, PARAMETER ::
     >  HALF = 0.5

C     Local variables:

      INTEGER
     >   I, ILEFT, IRIGHT, IT, J, JLEFT, JRIGHT, JT

      REAL
     >   DELX (-1:1, -1:2), DELY (-1:1, -1:2), ! dX, dY for 3 intervals/4 lines
     >   SLPX (-1:1, -1:2), SLPY (-1:1, -1:2), ! 2-pt. dF/dX, dF/dy values
     >   GRDX ( 0:1, -1:2), GRDY (-1:2,  0:1), ! 3-pt.   "      "
     >   DXY, DYX, GIJ, HIJ

C     Procedures:

      REAL, EXTERNAL ::
     >   DERIV3

C     Execution:

      IT = IEVAL ! Avoid excessive argument references
      JT = JEVAL ! The searching ensures that IEVAL <= I2 -1 and JEVAL <= J2 -1

C     Cell neighborhood information is stored as follows in each direction:
C
C        Point:      -1     0     1     2
C                 - - +-----+=====+-----+ - -
C        Interval:      -1     0     1

      ILEFT  =  -1
      IF (IT == I1)     ILEFT = 0
      IRIGHT =  2
      IF (IT == I2 - 1) IRIGHT = 1

      JLEFT  =  -1
      IF (JT == J1)     JLEFT = 0
      JRIGHT =  2
      IF (JT == J2 - 1) JRIGHT = 1

C     Normally there are 3 intervals in I at 4 J lines:

      DO J = JLEFT, JRIGHT
         DO I = ILEFT, IRIGHT - 1
            HIJ = X (IT+I+1, JT+J) - X (IT+I, JT+J)
            DELX (I, J) = HIJ
            SLPX (I, J) = (F (IT+I+1, JT+J) - F (IT+I, JT+J)) / HIJ
         END DO
      END DO

C     Similarly, there are normally 3 J intervals at 4 I lines.
C     The derivative utilities expect VECTORS of data, so the
C     Y-direction subscripts must be interchanged.

      DO I = ILEFT, IRIGHT
         DO J = JLEFT, JRIGHT - 1
            GIJ = Y (IT+I, JT+J+1) - Y (IT+I, JT+J)
            DELY (J, I) = GIJ
            SLPY (J, I) = (F (IT+I, JT+J+1) - F (IT+I, JT+J)) / GIJ
         END DO
      END DO

C     The cross derivatives require 4 x 2 or 2 x 4 stencils of 1st partials.
C     First, the X derivatives at the 0 & 1 I points for J lines -1 to 2.

      DO J = JLEFT, JRIGHT
         DO I = 0, 1
            GRDX (I, J) = DERIV3 (I, ILEFT, IRIGHT, DELX (-1, J),
     >                            SLPX (-1, J))
         END DO
      END DO

C     Now, the Y derivatives at the 0 & 1 J points for I lines -1 to 2:

      DO J = 0, 1
         DO I = ILEFT, IRIGHT
            GRDY (I, J) = DERIV3 (J, JLEFT, JRIGHT, DELY (-1, I),
     >                            SLPY (-1, I))
         END DO
      END DO

C     Set up the cell-specific/function-specific 4 x 4 matrix:

C     First, the upper left 2 x 2 (function values):

      FMATRIX (1, 1) = F (IT, JT)
      FMATRIX (2, 1) = F (IT + 1, JT)
      FMATRIX (1, 2) = F (IT, JT + 1)
      FMATRIX (2, 2) = F (IT + 1, JT + 1)

C     Next, the lower left 2 x 2 (dF/dp values):

      FMATRIX (3, 1) = DELX (0, 0) * GRDX (0, 0)
      FMATRIX (4, 1) = DELX (0, 0) * GRDX (1, 0)
      FMATRIX (3, 2) = DELX (0, 1) * GRDX (0, 1)
      FMATRIX (4, 2) = DELX (0, 1) * GRDX (1, 1)

C     Next, the upper right 2 x 2 (dF/dq values):

      FMATRIX (1, 3) = DELY (0, 0) * GRDY (0, 0)
      FMATRIX (1, 4) = DELY (0, 0) * GRDY (0, 1)
      FMATRIX (2, 3) = DELY (0, 1) * GRDY (1, 0)
      FMATRIX (2, 4) = DELY (0, 1) * GRDY (1, 1)

C     Now for the cross derivatives (lower right 2 x 2).
C     Replace the 2-pt. F slopes with 2-pt. "slopes" of the 1st derivatives
C     (normally 3 intervals along each of the 0 and 1 lines):

      DO J = 0, 1
         DO I = ILEFT, IRIGHT - 1
            SLPX (I, J) = (GRDY (I+1, J) - GRDY (I, J)) / DELX (I, J)
         END DO
      END DO

      DO I = 0, 1
         DO J = JLEFT, JRIGHT - 1
            SLPY (J, I) = (GRDX (I, J+1) - GRDX (I, J)) / DELY (J, I)
         END DO
      END DO

C     Average each pair of possible estimates.
C     DXY gives d/dp (dF/dq); DYX gives d/dq (dF/dp):

C     0,0:
      DXY = DERIV3 (0, ILEFT, IRIGHT, DELX (-1, 0), SLPX (-1, 0))
      DYX = DERIV3 (0, JLEFT, JRIGHT, DELY (-1, 0), SLPY (-1, 0))
      FMATRIX (3, 3) = (DXY + DYX) * HALF * DELX (0, 0) * DELY (0, 0)

C     1,0:
      DXY = DERIV3 (1, ILEFT, IRIGHT, DELX (-1, 0), SLPX (-1, 0))
      DYX = DERIV3 (1, JLEFT, JRIGHT, DELY (-1, 0), SLPY (-1, 0))
      FMATRIX (4, 3) = (DXY + DYX) * HALF * DELX (0, 0) * DELY (0, 1)

C     0,1:
      DXY = DERIV3 (0, ILEFT, IRIGHT, DELX (-1, 1), SLPX (-1, 1))
      DYX = DERIV3 (0, JLEFT, JRIGHT, DELY (-1, 1), SLPY (-1, 1))
      FMATRIX (3, 4) = (DXY + DYX) * HALF * DELX (0, 1) * DELY (0, 0)

C     1,1:
      DXY = DERIV3 (1, ILEFT, IRIGHT, DELX (-1, 1), SLPX (-1, 1))
      DYX = DERIV3 (1, JLEFT, JRIGHT, DELY (-1, 1), SLPY (-1, 1))
      FMATRIX (4, 4) = (DXY + DYX) * HALF * DELX (0, 1) * DELY (0, 1)

      END SUBROUTINE BICUBMAT
