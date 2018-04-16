C+----------------------------------------------------------------------
C
      SUBROUTINE PLGRID (NDIM, JDIM, XGRID, YGRID, JMIN, JINC, JMAX,
     >                   KMIN, KINC, KMAX, LINTYP, LINTHK, LINCOL,
     >                   MXWORK, RWORK)
C
C ACRONYM: Plot GRID (2D)
C               ----
C PURPOSE:
C    PLGRID is a general-purpose routine for plotting a 2D computational
C    mesh.  The coordinate arrays may be 1- or 2-dimensional, according
C    to whether the region is truly or only pseudo-rectangular.
C    PLGRID does not handle plotting of the boundary of the mesh or of
C    axes, etc.
C
C METHOD:
C    Assume the following prior set-up for this (sub)plot:
C      *  physical origin is set;
C      *  axes are scaled and drawn.
C    One grid line at a time is set up in work-space and plotted.
C
C ARGUMENTS:
C     ARG   TYPE I/O/S  DIM     DESCRIPTION
C    NDIM     I    I     -      No. of dimensions of XGRID and YGRID.
C    JDIM     I    I     -      Effective row dimension of XGRID(*,*)
C                               and YGRID(*,*) in the routine setting
C                               up this array.
C    XGRID,   R    I  JDIM,*    The coordinates are:
C    YGRID                      (  XGRID(J), YGRID(K) ) if NDIM=1;
C                               (XGRID(J,K),YGRID(J,K)) if NDIM=2.
C                               (* = 1 or KMAX.)
C    JMIN,JINC,JMAX,    -       Define  a  window  within  the  grid
C    KMIN,KINC,KMAX             arrays for the portion to be plotted
C                               and allow thinning of  dense  meshes.
C                               Set JMIN,KMIN,JINC,KINC = 1 for most
C                               cases.
C    LINTYP   I    I    2       LINTYP(1) and LINTYP(2) define the type
C                               of the J and K lines respectively as
C                               supplied by POLYLINE:
C                                   3  ...  Solid
C                                   4  ...  Dotted
C                                   5  ...  Dashed
C                                   6  ...  Chain-dotted
C                                   7  ...  Chain-dashed
C    LINTHK   I    I     2      LINTHK(1) and LINTHK(2) contain the line
C                               thickness multiples of the J and K lines
C                               respectively.  The actual line thickness
C                               is the product of LINTHK and .01 inches.
C    LINCOL   R    I     2      LINCOL(1) and LINCOL(2) specify the
C                               colors of the J and K lines respectively.
C                                   1. ...  Black
C                                   2. ...  Magenta
C                                   3. ...  Red
C                                   4. ...  Yellow
C                                   5. ...  Green
C                                   6. ...  Cyan
C                                   7. ...  Blue
C                                   8. ...  White
C                                with colors between the six primaries
C                                specified using real numbers in [2.,7.]
C    MXWORK   I    I     -       Calling program is responsible for
C                                checking that MXWORK is at least
C                                MAX(JMAX,KMAX)*2.
C    RWORK    R    S   MXWORK    Workspace for actual plotting arrays.
C                      
C
C EXTERNAL REFERENCES:
C    POLYLINE    Plots curve of indicated type/thickness/color.
C
C ENVIRONMENT:
C    VAX/VMS; FORTRAN 77; CA-DISSPLA graphics.
C
C HISTORY:
C    June 85    RCL      Original coding of PLGRID.
C    Mar. 87    RCL      Identical to the original PLGRID, but calls
C                        FMCURV to take advantage of custom line types
C                        and extended color selection.
C    DEc. 89    M.Wong   Calls to FMCURV replaced with calls to the
C                        generic curve plotting routine POLYLINE.
C   01/13/91  D.Saunders FMGRID renamed back to PLGRID for general use.
C   04/09/91   "   "     NSYM was undefined; K-line subscripts in RWORK
C                        args. to POLYLINE used KMIN/KMAX instead of
C                        JMIN/JMAX.
C
C AUTHOR: Rosalie Lefkowitz, NASA Ames/Sterling Software, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C     Arguments.

      INTEGER    LINTHK (2), LINTYP (2), JDIM, JINC, JMIN, JMAX, KINC, 
     >           KMIN, KMAX, MXWORK, NDIM 
      REAL       LINCOL (2), XGRID (JDIM, *), YGRID (JDIM, *),
     >           RWORK (MXWORK)

C     Local constants.

      INTEGER    NSYM
      PARAMETER (NSYM = -1)

C     Local variables.

      INTEGER    IER, J, K, LINE, NPTS, THICK
      REAL       COLOR

C  *  Execution:

      COLOR = LINCOL(1)
      IF (COLOR .EQ. 8.) COLOR = 0.
      NPTS = KMAX - KMIN + 1

C  *  Draw J lines first:

      DO 50, J = JMIN, JMAX, JINC
         IF (NDIM .EQ. 1) THEN
            DO 20, K = KMIN, KMAX
               RWORK (K) = XGRID (J, 1)
   20       CONTINUE
            IF (J .EQ. JMIN) THEN
               DO 30, K = KMIN, KMAX
                  RWORK (K + KMAX) = YGRID (K, 1)
   30          CONTINUE
            END IF
         ELSE
            DO 40, K = KMIN, KMAX
               RWORK (K)        = XGRID (J, K)
               RWORK (K + KMAX) = YGRID (J, K)
   40       CONTINUE
         END IF

         CALL POLYLINE (NPTS, RWORK (KMIN), RWORK (KMIN + KMAX),
     >                  LINTYP (1), LINTHK (1), NSYM, COLOR, IER)
   50 CONTINUE

C  *  Now the K lines: 

      COLOR = LINCOL (2)
      IF (COLOR .EQ. 8.) COLOR = 0.
      NPTS = JMAX - JMIN + 1

      DO 150, K = KMIN, KMAX, KINC
         IF (NDIM .EQ. 1) THEN
            DO 120, J = JMIN, JMAX
               RWORK (J + JMAX) = YGRID (K, 1)
  120       CONTINUE
            IF (K .EQ. KMIN) THEN
               DO 130, J = JMIN, JMAX
                  RWORK (J)= XGRID (J, 1)
  130          CONTINUE
            END IF
         ELSE
            DO 140, J = JMIN, JMAX
               RWORK (J)        = XGRID (J, K)
               RWORK (J + JMAX) = YGRID (J, K)
  140       CONTINUE
         END IF

         CALL POLYLINE (NPTS, RWORK (JMIN), RWORK (JMIN + JMAX),
     >                  LINTYP (2), LINTHK (2), NSYM, COLOR, IER)
  150 CONTINUE

      RETURN
      END
