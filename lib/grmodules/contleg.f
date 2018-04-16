C+**********************************************************************
C
      SUBROUTINE CONTLEG ( NLEVS, RLEVS, NPRLEV, NDGLBL, XLEGND, YLEGND, 
     >                     TXLGND, HTLGND, LINCOL, LINTHK, LINTYP )
C
C  PURPOSE:  
C     CONTLEG produces a contour plot legend which shows a list of contour level
C     numbers and values along with the associated line type/colors.  The legend
C     format used here is contour number, line type and value.  For example,
C
C                               1 _____ 1.000
C                               2 ..... 2.000
C                               etc.
C
C  ARGUMENTS:
C    NAME   TYPE  I/O/S   DIM    DESCRIPTION
C    NLEVS   I      I      -     Number of contour levels.
C    RLEVS   R      I    NLEVS   Actual contour levels.
C    NPRLEV  I      I    NLEVS   Number of contour lines per level.
C    NDGLBL  I      I      -     Number of digits for contour labels.
C    XLEGND  R      I      -     Horizontal and vertical position of 
C    YLEGND  R      I      -     upper left corner of legend. (inches)
C    TXLGND  C*10   I      -     Legend text (10 characters).
C    HTLGND  R      I      -     Height of legend and labels.
C                                0.0 means take default of 0.07.
C    LINCOL  R      I    NLEVS   LINCOL(I) specifies the color of the
C                                Ith contour line  plotted.
C                                Any number outside the range 1.-8.
C                                will be ignored by POLYLINE, which will
C                                use the output device default color.
C                                  1.   ...  Black
C                                  2.   ...  Magenta
C                                  3.   ...  Red
C                                  4.   ...  Yellow
C                                  5.   ...  Green
C                                  6.   ...  Cyan
C                                  7.   ...  Blue
C                                  8.   ...  White
C                               with colors between the six primaries
C                               specified using real numbers in [2.,7.]
C    LINTHK   I    I     NLEVS  Defines the thickness of the line drawn
C                               in multiples of .01 inch. 1:9 is valid.
C    LINTYP   I    I     NLEVS  The line type corresponding to the current
C                               contour level.  (Note:  Use of line types 
C                               supplied by POLYLINE requires line type
C                               3 to correspond to contour level 1, etc.)
C                                Contour Level    Line Type      Pattern
C                                -------------    ---------      -------
C                                     1              3           Solid line
C                                     2              4           Dots
C                                     3              5           Dashes
C                                     4              6           Chain-dots
C                                     5              7           Chain-dashes
C                                     6              8           Long-dashes
C
C  EXTERNAL REFERENCES:
C    ANGLE  - For ensuring horizontal legend lines.           DISSPLA
C    GRACE  - Sets plot area grace margin.                    DISSPLA
C    HEIGHT - Sets height of text.                            DISSPLA
C    INTNO  - Writes integer number (from physical origin).   DISSPLA
C    MESSAG - Writes message (from physical origin).          DISSPLA
C    REALNO - Writes real number (from physical origin).      DISSPLA
C    XMESS  - Outputs length of message.                      DISSPLA
C    XINVRS, YINVRS - Convert plot position in inches to      DISSPLA
C             plot position in data units.
C    POLYLINE - Plots curve of indicated type/thickness/color.
C
C  STANDARDS VIOLATIONS: None.
C
C  ENVIRONMENT:  VAX/VMS -- FORTRAN 77
C                CA Graphics DISSPLA version 11.0
C
C  HISTORY:
C  06/12/89    MDW      Call to FMCURV replaced by call to POLYLINE
C  06/11/89    MDW      Updated to run on DISSPLA version 11.0.
C                       Took out %REF's and Hollerith 
C  03/23/87    RCL      LINCOL is now type REAL.
C  10/07/86    RCL      Added output of line types/colors. Clarified
C                       use of XLEGND, YLEGND.
C  01/09/81    CJG      Original coding.
C
C  AUTHOR: Carol Green, Ron Langhi, Rosalie Lefkowitz, Informatics Inc.
C
C-**********************************************************************
 
      IMPLICIT NONE

C  *  Arguments

      CHARACTER  
     >   TXLGND * 10
      INTEGER
     >   NLEVS, LINTHK (NLEVS), LINTYP (NLEVS), NPRLEV (NLEVS), 
     >   NDGLBL
      LOGICAL    
     >   NOLINE
      REAL       
     >   LINCOL (NLEVS), RLEVS (NLEVS), XLEGND, YLEGND, HTLGND

C  *  Local variables

      INTEGER
     >   LEVEL, IER
      REAL
     >   ALENG, COLOR, HLEG, X, XDAT (2), XINVRS, XLINE1, XLINE2, 
     >   XMESS, XTEXT, YDAT (2), YINC, YINVRS, YLINE, YTEXT
 
C  *  Extend grace margin beyond edge of page to permit the showing 
C  *  of lines in the legend.

      CALL GRACE ( 20.)
      CALL ANGLE ( 0.0 )
      HLEG   = HTLGND
      IF ( HLEG.LE.0. ) HLEG = .07
      CALL HEIGHT ( HLEG )
      YINC   = HLEG + HLEG
      ALENG  = XMESS( 'CONTOUR', 7 )
      XLINE1 = XLEGND + .4*ALENG
      XLINE2 = XLEGND + 1.3*ALENG
      XTEXT  = XLEGND + 1.5 * ALENG
      YTEXT  = YLEGND - HLEG
      YLINE  = YLEGND - 0.3*HLEG

C  *  Print legend headings:

      CALL MESSAG ( 'CONTOUR', 7, XLEGND, YTEXT )
      CALL MESSAG ( TXLGND, 10, XTEXT, YTEXT )

C  *  Fill in legend values:

      NOLINE = .FALSE.

      DO 10 LEVEL = 1, NLEVS

C  *     Level number:

         YTEXT = YTEXT - YINC
         CALL INTNO ( LEVEL, XLEGND, YTEXT )

C  *     Show line type:

         YLINE = YLINE - YINC
         YDAT(1) = YINVRS(XLINE1,YLINE)
         YDAT(2) = YINVRS(XLINE2,YLINE)
         XDAT(1) = XINVRS(XLINE1,YLINE)
         XDAT(2) = XINVRS(XLINE2,YLINE)

         COLOR = LINCOL (LEVEL)
         IF (LINCOL(LEVEL) .EQ. 8)  COLOR = 0    ! White
         CALL POLYLINE (2, XDAT, YDAT, LINTYP (LEVEL), 
     >      LINTHK (LEVEL), -1, COLOR, IER)

C  *     Contour value:

         CALL REALNO ( RLEVS( LEVEL ), NDGLBL, XTEXT, YTEXT )

C  *     Mark with asterisk if level not found:

         IF ( NPRLEV( LEVEL ).EQ.0 ) THEN
            X = XTEXT + ALENG / 7. * ( NDGLBL + 3 )
            CALL MESSAG( ' *', 2, X, YTEXT )
            NOLINE = .TRUE.
         END IF
   10 CONTINUE

      IF ( NOLINE ) THEN
         YTEXT = YTEXT - 1.5 * YINC
         CALL MESSAG ( '*  NO CONTOUR LINES FOUND', 25, XLEGND, YTEXT )
      END IF

C  *  Restore appropriate grace margin.

      CALL GRACE ( 0.0 )

      RETURN
      END
