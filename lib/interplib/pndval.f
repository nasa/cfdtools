C+----------------------------------------------------------------------
C
      SUBROUTINE PNDVAL (NDEG, COEFS, NPTS, X, Y, YP, YPP)
C
C ACRONYM:  PolyNomial and Derivatives eVALuation
C           -   -          -            ---
C PURPOSE:  PNDVAL evaluates the polynomial curve defined by NDEG and
C           COEFS(*), at the abscissas supplied in X(*).  It also re-
C           turns values of the 1st and 2nd derivatives at each X(I).
C           (Use PNEVAL if you don't need derivatives.)
C
C METHOD:   Horner's rule (valid for NDEG>=2 for YPP, YP, and Y, with
C           special treatment for NDEG=1 or 0), for each abscissa.
C     
C ARGUMENTS:
C    ARG    DIM      TYPE I/O/S DESCRIPTION 
C   NDEG     -        I     I   Degree of polynomial (>=0)
C   COEFS  0:NDEG     R     I   Coefficients in the order
C                               C(0), C(1), ..., C(NDEG)
C   NPTS     -        I     I   No. of points to evaluate curve at (>=1)
C   X       NPTS      R     I   Abscissas
C   Y       NPTS      R     O   Ordinates on the curve
C   YP      NPTS      R     O   Corresponding 1st derivatives
C   YPP     NPTS      R     O   Corresponding 2nd derivatives.  (Pass YP
C                               here if YPP is not needed.)
C    
C ENVIRONMENT:  VAX/VMS FORTRAN; IMPLICIT NONE is nonstandard.
C
C HISTORY:
C   June '83  D.A.Saunders  Adapted from PNEVAL.
C   09/30/88    "     "     Revised to assign YPP, YP, Y in that order.
C   07/18/89    "     "     Cosmetics in preparation for Mac. translation.
C
C AUTHOR: David Saunders, Informatics/NASA Ames, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C  *  Arguments:

      INTEGER
     >   NDEG, NPTS
      REAL
     >   COEFS (0:NDEG), X (NPTS), Y (NPTS), YP (NPTS), YPP (NPTS)

C  *  Local constants:

      REAL
     >   ZERO
      PARAMETER
     >  (ZERO = 0.E+0)

C  *  Local variables:

      INTEGER
     >   I, J
      REAL
     >   XI

C  *  Execution:

C  *  Treat NDEG=0 or 1 separately so that Horner's rule can safely
C     be applied to each of Y, YP, YPP for the general case:

      IF (NDEG .EQ. 0) THEN

         DO 20 I = 1, NPTS
            YPP (I) = ZERO
            YP (I)  = ZERO
            Y (I)   = COEFS (0)
   20    CONTINUE

      ELSE IF (NDEG .EQ. 1) THEN

         DO 30 I = 1, NPTS
            YPP (I) = ZERO
            YP (I)  = COEFS (1)
            Y (I)   = COEFS (1) * X (I) + COEFS (0)
   30    CONTINUE

      ELSE

         DO 50 I = 1, NPTS
            XI      = X (I)
            YPP (I) = ZERO
            YP (I)  = COEFS (NDEG) * REAL (NDEG)
            Y (I)   = COEFS (NDEG) * XI + COEFS (NDEG-1)
            DO 40 J = NDEG, 2, -1
               YPP (I) = YPP (I) * XI + COEFS (J) * REAL (((J-1)*J))
               YP (I)  = YP (I)  * XI + COEFS (J-1) * REAL ((J-1))
               Y (I)   = Y (I)   * XI + COEFS (J-2)
   40       CONTINUE
   50    CONTINUE

      END IF

      RETURN
      END
