C+**********************************************************************
C
      SUBROUTINE COPY ( NX, X1, X2 )
C
C PURPOSE:  COPY copies the elements of X1(*) into X2(*).
C
C PARAMETERS:
C    ARG   DIM   TYPE I/O/S DESCRIPTION 
C     NX    -      I    I   No. of elements in input and output arrays
C     X1   NX      R    I   Array to be copied
C     X2   NX      R    O   Copy of X1(*)
C    
C ENVIRONMENT: FORTRAN IV
C
C AUTHOR: David Saunders, Informatics, Palo Alto, CA.  (07/30/82)
C
C-**********************************************************************
C
      DIMENSION  X1(NX), X2(NX)
C
      DO 20 I=1,NX
         X2(I) = X1(I)
 20   CONTINUE
C
      RETURN
      END
