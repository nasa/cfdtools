C+----------------------------------------------------------------------
C
      SUBROUTINE SCLEAN (I)
C
C  PURPOSE:  SCLEAN writes escape sequences which erase the screen of a
C            DEI VT640  "green screen"  graphics terminal. The graphics
C            display  can  be left  in  place or  erased along with the
C            VT100  display.  In both cases  the  terminal is placed in
C            VT100  mode.   A  separate call can  be used to return the
C            terminal to alpha mode.
C
C            <This is obscure.  Did Rosalie mean the I=2 case, which
C            returns the terminal to vector graphics mode?  DAS>
C
C            Logical unit 6 is hard-coded.
C
C  ARGUMENTS:
C   VAR  TYPE  DIM  I/O/S  DESCRIPTION
C    I    I     -     I    1 means  go from  alpha mode  to transparent
C                            mode,  erase all of the VT100 display, but
C                            leave the graphics display with the cursor
C                            set to home.
C                          2 means return to TEK 4014 vector mode.
C                          3 means wipe screen and return to transparent
C                            (i.e. VT100) mode.
C
C  ENVIRONMENT:  DEC VAX/VMS; FORTRAN; VT640 and Selanar termiinals.
C         Also:  SGI IRIS/IRIX; GraphOn GO-200 terminals.
C
C  HISTORY:
C    Oct. 1983    RCL     Original design and coding.
C    Dec. 1985    RCL     Changed SCLEAN(3) to be compatible with 
C                         Selanar terminals as well as VT640's.
C    May 1 '91 D.Saunders Writing integer ASCII codes didn't work on the
C                         IRIS - use CHAR(n) instead of n in the I/O list.
C
C  AUTHOR: Rosalie Lefkowitz, Informatics General Corp.
C
C-----------------------------------------------------------------------

C  *  Arguments:

      INTEGER I

C  *  Execution:

      IF (I .EQ. 1) THEN

C  *     Go from alpha mode to transparent (VT100) mode, erase all of the
C        VT100 display (leaving the graphics display), and home the cursor.
C        The sequence is:       ESC, ESC, FF, CAN, ESC, [, 2, J, ESC, [, H

         WRITE (6, 1001) CHAR (27), CHAR (27), CHAR (12), CHAR (24),
     >      CHAR (27), CHAR (91), CHAR (50), CHAR (74), CHAR (27),
     >      CHAR (91), CHAR (72)

      ELSE IF (I .EQ. 2) THEN

C  *     Return to TEK 4014 vector mode.      
C                        ESC,  ESC, GS

         WRITE (6, 1001) CHAR (27), CHAR (27), CHAR (29)

      ELSE IF (I .EQ. 3) THEN

C  *     Clear graphics screen and return to transparent (VT100) mode.
C                        GS, ESC, FF, CAN

         WRITE (6, 1001) CHAR (29), CHAR (27), CHAR (12), CHAR (24)
      END IF

 999  RETURN

1001  FORMAT ('+', 11A1)
      END
