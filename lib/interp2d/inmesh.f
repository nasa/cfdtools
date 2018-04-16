C+----------------------------------------------------------------------
C
      LOGICAL FUNCTION INMESH ( X0, Y0, MAXI, NI, NJ, X, Y, IC, JC )
C
C PURPOSE: INMESH  is returned TRUE if the point (X0,Y0) is found to lie
C          in one of the cells of the given pseudo-rectangular mesh.  It
C          indicates which cell was identified, but leaves anything fur-
C          ther up to the application.  If no such cell is found, INMESH
C          is returned as FALSE.
C
C          Applications for which this capability has been found  useful
C          include computation of velocity profiles along normals to  an
C          airfoil surface, where the mesh wraps around the airfoil; and
C          tracking of the wake behind a rotor blade in forward  flight,
C          where the wake is represented as a curvilinear mesh,  and the
C          problem becomes one of determining whether or not given comp-
C          utational mesh points lie in the wake.
C
C ARGUMENTS:
C ARG    DIM     TYPE  I/O/S  DESCRIPTION
C X0,Y0   -        R     I    Coordinates of point being tested
C MAXI    -        I     I    Allows for using less than full mesh array
C NI,NJ   -        I     I    Indicate portions of mesh array in use
C X,Y   NI,NJ      R     I    X,Y coordinates of the mesh
C IC,JC   -        I     O    If (X0,Y0) is found to lie in a cell, then
C                             that cell has defined by mesh points IC,JC
C                             and IC-1,JC and IC,JC-1 and IC-1,JC-1
C
C METHOD:  For each cell of the mesh, treated as a quadrilateral,  the 4
C          triangles formed by connecting the 4 vertices to (X0,Y0)  are
C          considered. Determinants related to the formulas for the con-
C          ditions that  (X0,Y0)  lies on each side of the quadrilateral
C          will all have the same (positive) sign if  (X0,Y0) is in fact
C          inside the quadrilateral.
C
C SOURCE:  Len Wallitt/Lyn King, FAR Branch, NASA Ames, c. 1979.
C          This version: David Saunders, Informatics, 11/09/83.
C
C-----------------------------------------------------------------------
C
      DIMENSION X(MAXI,NJ), Y(MAXI,NJ)
C
      INMESH = .FALSE.
C
      DO 300 J = 2, NJ
         DO 200 I = 2, NI
            IF ( ( X0-X(I,J) )*( Y(I,J-1)-Y(I,J) ) -
     +           ( Y0-Y(I,J) )*( X(I,J-1)-X(I,J) ) .LT. 0.E0 )
     +         GO TO 200
            IF ( ( X0-X(I-1,J) )*( Y(I,J)-Y(I-1,J) ) -
     +           ( Y0-Y(I-1,J) )*( X(I,J)-X(I-1,J) ) .LT. 0.E0 )
     +         GO TO 200
            IF ( ( X0-X(I,J-1) )*( Y(I-1,J-1)-Y(I,J-1) ) -
     +           ( Y0-Y(I,J-1) )*( X(I-1,J-1)-X(I,J-1) ) .LT. 0.E0 )
     +         GO TO 200
            IF ( ( X0-X(I-1,J-1) )*( Y(I-1,J)-Y(I-1,J-1) ) -
     +           ( Y0-Y(I-1,J-1) )*( X(I-1,J)-X(I-1,J-1) ) .LT. 0.E0 )
     +         GO TO 200
            INMESH = .TRUE.
            IC = I
            JC = J
            GO TO 999
  200    CONTINUE
  300 CONTINUE
C
  999 RETURN
      END
