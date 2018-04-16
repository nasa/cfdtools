C+----------------------------------------------------------------------
C
      SUBROUTINE PLPOLY ( NPTS, XINCH, YINCH, THKNSS, ICOLOR )
C
C  PURPOSE:  PLPOLY draws the outline of a polygon and blanks out the
C            interior. Color and line thickness are user options.
C            
C  ARGUMENTS:
C    ARG    TYPE  I/O/S   DIM     DESCRIPTION
C   NPTS     I      I      -      No. of points in curve to be drawn.
C                                 A connection between the first and
C                                 last point is assumed. No repetition
C                                 point is needed.
C   XINCH,   R      I    NPTS     Abscissas and ordinates of curve to
C    YINCH                        be drawn. (In inches NOT data values.)
C   THKNSS   R      I      -      Thickness of curve in inches.
C                                 THKNSS<0. means take system default
C                                 pen width. 
C   ICOLOR   I      I      -      0 means use current color ( may be
C                                 default black ); color routines will
C                                 not be called here. Else
C                                    1   ...  Black
C                                    2   ...  Magenta
C                                    3   ...  Red
C                                    4   ...  Yellow
C                                    5   ...  Green
C                                    6   ...  Cyan
C                                    7   ...  Blue
C                                    8   ...  White
C
C  EXTERNAL REFERENCES:
C    DISSPLA Routines:
C    BLPOLY - Draws outline of closed curve and blanks interior.
C    RESET  - Resets specified option to default value.
C    SETCLR - Sets hardware color.
C
C  ENVIRONMENT:  VAX/11 FORTRAN V3.0
C
C  DEVELOPMENT HISTORY:
C    DATE   INITIALS  DESCRIPTION
C  09/01/89  MDW      Updated to run on DISSPLA version 11.0.
C  02/25/85  RCL      Original Design and Coding
C
C  AUTHOR:  Rosalie Lefkowitz  Informatics General Corporation
C
C-----------------------------------------------------------------------
C
      LOGICAL    LCOLOR

      DIMENSION  XINCH(NPTS), YINCH(NPTS)

      CHARACTER  COLARY(8)*7

C  *  These are the colors available from DISSPLA.

      DATA       COLARY /'BLACK','MAGENTA','RED','YELLOW','GREEN',
     +                   'CYAN','BLUE','WHITE'/

C
C  *  Check for color request.
C
      LCOLOR = ( ICOLOR.GT.0 .AND. ICOLOR.LE.8 )

      IF ( LCOLOR ) CALL SETCLR ( COLARY(ICOLOR) )

      CALL BLPOLY ( XINCH, YINCH, NPTS, THKNSS )
C
      IF ( LCOLOR ) CALL RESET ('SETCLR')
C
      RETURN
      END
