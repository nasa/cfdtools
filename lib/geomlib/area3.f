C+-----------------------------------------------------------------------
C
      REAL FUNCTION AREA3 (X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3)
C
C ONE-LINER:
C     This function routine finds the area of a triangle in 3-space.
C
C PURPOSE:
C     AREA3 calculates the area of a triangle P1, P2 and P3 in space. 
C     Vertices may be entered in any order.
C 
C ARGUMENTS:
C     ARG       DIM     TYPE    I/O      DESCRIPTION
C     X1         -        R      I       X coordinate of P1
C     Y1         -        R      I       Y     "         "
C     Z1         -        R      I       Z     "         "
C     .                   .      .       .
C     .                   .      .       .
C     .                   .      .       .
C     Z3         -        R      I       Z coordinate of P3
C 
C METHOD:
C     The function calculates the area of a triangle using the cross product 
C     method. (Reference: "Standard Mathematical Tables", Beyer, P. 204).
C
C ERROR HANDLING:  None
C
C EXTERNAL REFERENCES:  None
C
C ENVIRONMENT:  VAX/VMS FORTRAN 77
C
C HISTORY:
C     4/3/89     M.Wong     Initial design and implementation
C
C AUTHOR:   Michael D. Wong, Sterling Software, Palo Alto, CA
C
C------------------------------------------------------------------------

      IMPLICIT NONE

C Arguments
C ---------
      REAL       X1,X2,X3,
     >           Y1,Y2,Y3,
     >           Z1,Z2,Z3


C Execution
C ---------
      
      AREA3 = .5 * SQRT ((Y1 * (Z2 - Z3)  -  Z1 * (Y2 - Y3)  +
     >                   (Y2 * Z3 - Y3 * Z2)) ** 2  
     >                            +
     >                   (Z1 * (X2 - X3) - X1 * (Z2 - Z3) +
     >                   (Z2 * X3 - Z3 * X2)) ** 2          
     >                            +
     >                   (X1 * (Y2 - Y3) - Y1 * (X2 - X3) + 
     >                   (X2 * Y3 - X3 * Y2)) ** 2 )


      RETURN
      END
