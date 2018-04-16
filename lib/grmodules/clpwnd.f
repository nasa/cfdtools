C+----------------------------------------------------------------------
      LOGICAL FUNCTION CLPWND (CLIPXL, CLIPXR, CLIPYB, CLIPYT,
     1                         X1, X2, Y1, Y2)
C
C  PURPOSE -
C    Performs line clipping on the rectangular window bounded by the four
C    variables CLIP**, sets CLPWND to true if at least part of the line
C    is in the window, and sets (X1, Y1) and (X2, Y2) to the endpoints of
C    the clipped line.
C
C  PARAMETERS -
C    ARG  TYPE I/O/S  DIM         DESCRIPTION
C   CLPWND L   RETURN  -       TRUE = at least part of line lies in the window.
C			       FALSE= line is completely outside window.
C   CLIPXL R     I     -       Minimum X-value of window
C   CLIPXR R     I     -       Maximum X-value of window
C   CLIPYB R     I     -       Minimum Y-value of window
C   CLIPYT R     I     -       Maximum Y-value of window
C   X1     R    I/O    -       Original first endpoint X of line on input.
C   X2     R    I/O    -       Original second endpoint X of line on input.
C   Y1     R    I/O    -       Original first endpoint Y of line on input.
C   Y2     R    I/O    -       Original second endpoint Y of line on input.
C				  IF CLPWND is true then on output
C				  (X1, Y1) is first endpoint of the
C				  clipped line, ELSE (X1, Y1) have
C				  unpredictable values.
C				  IF CLPWND is true then on output
C				  (X2, Y2) is second endpoint of the
C				  clipped line, ELSE (X2, Y2) have
C				  unpredictable values.
C
C  NOTES:
C	This is a FORTRAN implementation of an algorithm invented by Dan Cohen
C	and Ivan Sutherland.  The original Pascal code was copied, with minor
C	revision for implementation on VAX, from Principles of Interactive
C	Computer Graphics (Second Edition) New York, 1979, pp 66-67.
C
C	Comments beginning with "CF" are Pascal source statements left in
C	this version for clarity.
C
C  NON-STANDARD CODE:
C
C	The Pascal set type is implemented by using a subset of the powers
C	of two and checking for their presence with a potentially machine-
C	specific bitwise operation: the intrinsic function IAND.
C
C  EXTERNAL REFERENCES -
C    CODE	(in this source file) sets a non-zero code in C
C		 if point (X,Y) is out of window.
C
C  ENVIRONMENT:  VAX 11/780, VMS V4.4, FORTRAN V4.5 (and earlier)
C
C  AUTHOR - Dexter L. Hermstad		Informatics General Corporation
C		04/15/83		Palo Alto, California
C
C  REVISIONS -
C    DATE      PERSON  STATEMENT OF CHANGES
C    05/22/87	DLH	FORTRAN version of original Pascal source.
C    05/22/84	DLH	Corrected subroutine header for argument order
C    04/15/83	DLH	Original coding
C
C-----------------------------------------------------------------------
C
      REAL CLIPXL, CLIPXR, CLIPYB, CLIPYT, X1, X2, Y1, Y2
CF   type edge=(left,right,bottom,top); outcode = set of edge;
      INTEGER LEFT, RIGHT, BOTTOM, TOP
      PARAMETER (LEFT=1, RIGHT=2, BOTTOM=4, TOP=8)
CF   var c,c1,c2:outcode; x,y:real;
      INTEGER C, C1, C2
      REAL X, Y
C
CF   procedure Code (x,y:real; var c:outcode);
CF   begin
CF      c:= [];
CF      if x<Clipxl then c:= [left] else if x>Clipxr then c:= [right];
CF      if y<Clipyb then c:= c+ [bottom] else if y>Clipyt then c:= c+ [top];
CF   end;
      CLPWND= .FALSE.
      CALL CODE (CLIPXL, CLIPXR, CLIPYB, CLIPYT, X1, Y1, C1)
      CALL CODE (CLIPXL, CLIPXR, CLIPYB, CLIPYT, X2, Y2, C2)
C
CF   while (c1<>[]) or (c2<>[]) do begin
  100 CONTINUE
      IF (C1.NE.0 .OR. C2.NE.0) THEN
CF         if (c1*c2)<>[] then goto 800;
         IF (IAND (C1, C2).NE.0) GO TO 800
         C= C1
         IF (C.EQ.0) C= C2
CF      if left in c then begin {Crosses left edge}
         IF (IAND (LEFT, C).NE.0) THEN
            Y= Y1+ (Y2-Y1)* (CLIPXL-X1)/(X2-X1)
            X= CLIPXL
CF      if right in c then begin {Crosses right edge}
         ELSE IF (IAND (RIGHT, C).NE.0) THEN
            Y= Y1+ (Y2-Y1)* (CLIPXR-X1)/(X2-X1)
            X= CLIPXR
CF      if bottom in c then begin {Crosses bottom edge}
         ELSE IF (IAND (BOTTOM, C).NE.0) THEN
            X= X1+ (X2-X1)* (CLIPYB-Y1)/(Y2-Y1)
            Y= CLIPYB
CF      if top in c then begin {Crosses top edge}
         ELSE IF (IAND (TOP, C).NE.0) THEN
            X= X1+ (X2-X1)* (CLIPYT-Y1)/(Y2-Y1)
            Y= CLIPYT
         ENDIF
         IF (C.EQ.C1) THEN
            X1= X
            Y1= Y
            CALL CODE (CLIPXL, CLIPXR, CLIPYB, CLIPYT, X, Y, C1)
         ELSE
            X2= X
            Y2= Y
            CALL CODE (CLIPXL, CLIPXR, CLIPYB, CLIPYT, X, Y, C2)
         ENDIF
CF   end;
      GO TO 100
      ENDIF
CF   { If we reach here, the line from (x1,y1) to (x2,y2) is visible}
      CLPWND= .TRUE.
800   CONTINUE
      RETURN
      END
C+
CF   procedure Code (x,y:real; var c:outcode);
      SUBROUTINE CODE (CLIPXL, CLIPXR, CLIPYB, CLIPYT, X, Y, C)
C
C   PURPOSE: CODE sets a non-zero code in C if point (X,Y) is out of window.
C
C   NOTES:
C     CLIPXL, CLIPXR, CLIPYB, CLIPYT added because in FORTRAN this is not an
C     internal subroutine that shares use of dummy arguments passed to CLPWND.
C-
      REAL CLIPXL, CLIPXR, CLIPYB, CLIPYT
      REAL X, Y
      INTEGER C
C-
      INTEGER LEFT, RIGHT, BOTTOM, TOP
      PARAMETER (LEFT=1, RIGHT=2, BOTTOM=4, TOP=8)
C
      C= 0
CF      if x<Clipxl then c:= [left] else if x>Clipxr then c:= [right];
      IF (X.LT.CLIPXL) C= LEFT
      IF (X.GT.CLIPXR) C= RIGHT
CF      if y<Clipyb then c:= c+ [bottom] else if y>Clipyt then c:= c+ [top];
      IF (Y.LT.CLIPYB) C= C+ BOTTOM
      IF (Y.GT.CLIPYT) C= C+ TOP
      RETURN
      END
