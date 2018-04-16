C+----------------------------------------------------------------------
C
      SUBROUTINE LINTRP ( NY, YM, YP, XM, XP, X, Y )
C
C PURPOSE:  LINTRP  interpolates linearly for all elements of the given
C           parallel arrays YM(*), YP(*), into the third parallel array
C           Y(*),  according to the abscissas XM, XP, X that are common
C           to all elements of the corresponding array...   (Originally
C           developed for interpolating wing pressure data from comput-
C           ational span stations to arbitrary span stations.)
C
C PARAMETERS:
C   ARG   DIM   TYPE I/O/S DESCRIPTION
C   NY     -      I    I   No. of interpolated elements required (>=1).
C   YM,YP  NY     R    I   Parallel arrays between which Y(*) is reqd.
C   XM,XP  -      R    I   Abscissas applying to YM(*), YP(*) resp.
C   X      -      R    I   Abscissa applying to desired array Y(*).
C                          Normally XM < X < XP, but extrapolation is
C                          permissible with this code.
C   Y      NY     R    O   Interpolated value(s).
C    
C ENVIRONMENT: FORTRAN IV
C
C AUTHOR: David Saunders, Informatics, Palo Alto, CA.   (07/30/82)
C
C-----------------------------------------------------------------------
C
      DIMENSION  YM(NY), YP(NY), Y(NY)
C
      RATIO = (X-XM) / (XP-XM)
C
      DO 20 I = 1, NY
         Y(I) = YM(I) + (YP(I)-YM(I))*RATIO
 20   CONTINUE
C
      RETURN
      END
