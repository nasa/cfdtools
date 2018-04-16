C+---------------------------------------------------------------------
C
      SUBROUTINE HARDY ( XD,YD,ZD,ND,XI,YI,ZI,NI,IWK,WK,R,IER,COND )
C
C  PURPOSE: Does bivariate interpolation via Hardy's multiquadratic 
C           method:  basically fitting as many linear parameters as
C           there are data points, with basis functions of the form
C
C           PHI(K,X,Y) = SQRT( (X-XD(K))**2 + (Y-YD(K))**2 + R**2 )
C
C           The data can then be represented by the function
C
C                    G(X,Y) = SUM (K=1:N) A(K)*PHI(K,X,Y)
C
C           where  A(K)  is the Kth linear parameter computed here.
C
C           This function is evaluated at the given point(s).
C
C  NOTES:   The original data points need not be uniformly  spaced,
C           but since the function used interpolates the data, this
C           method is not likely to be appropriate for noisy data.
C
C  PARAMETERS:
C       ARG    TYPE  I/O/S   DIM     DESCRIPTION
C
C        XD    R*4     I      ND     Input X locations of data points.
C
C        YD    R*4     I      ND     Input Y locations of data points.
C
C        ZD    R*4     I      ND     Input function values at data points.
C
C        ND     I      I      -      Number of data points.
C
C        XI    R*4     I      NI     X locations of points at which to 
C                                    approximate.
C
C        YI    R*4     I      NI     Y locations of points at which to 
C                                    approximate.
C
C        ZI    R*4     O      NI     Approximated value at points (XI(I),YI(I))
C
C        NI     I      I      -      Number of interpolation points.
C
C        IWK   I*4     S      ND     Pivot vector for solving square system.
C
C        WK    R*4     S   ND*(ND+1) For setting up and solving square system.
C
C         R    R*4     I      -      Shape parameter. Try R in the range
C                                    0 to 2.
C
C        IER    I      O      -      Error return code. 0 means no error.
C                                    129 means matrix is singular. The 
C                                    probable cause is that two of the
C                                    data points are the same.
C
C        COND   R*4    O      -      Condition number of matrix involved.
C
C  COMMONS USED: None
C
C  FILES USED: None
C
C  EXTERNAL REFERENCES: DECOMP,SOLVE
C
C  ENVIRONMENT: VAX 11/780; FORTRAN
C
C  KNOWN MACHINE DEPENDENCIES: END DO is a VAX extension.
C
C  AUTHOR:   Douglas L. Richter      Informatics General Corporation
C
C  DEVELOPMENT HISTORY:
C       DATE   INITIALS  DESCRIPTION
C     04/26/83    DLR    Original design and coding
C     09/28/83    DLR    Used DECOMP and SOLVE instead of IMSL routine
C     12/28/83    DAS    Expanded method description; made HGK in-line;
C                        avoided overwriting the original function values.
C
C-----------------------------------------------------------------------
C
      DIMENSION XD(ND), YD(ND), ZD(ND), XI(NI), YI(NI), ZI(NI), IWK(ND),
     1          WK(ND,ND+1)
C
C     Statement function for Kth basis function:
C
      HGK ( X, Y, K ) = SQRT ( (X-XD(K))**2 + (Y-YD(K))**2 + RSQ )
C
C     Note that HGK ( XD(I), YD(I), J ) = HGK ( XD(J), YD(J), I )
C     so the left hand side of the system being solved is symmetric:
C
      RSQ = R*R
C
      DO I = 1, ND
         DO J = 1, I
            WK(I,J) = HGK ( XD(J), YD(J), I )
            WK(J,I) = WK(I,J)
         END DO
      END DO
C
C     Factorize the matrix:
C
      CALL DECOMP ( ND, ND, WK, COND, IWK, WK(1,ND+1) )
C
      IER = 0
      IF ( COND.GE.1.E+32 ) THEN
         IER = 129
         RETURN
      END IF
C
C     Else set up the right hand side and solve:
C
      DO I = 1, ND
         WK(I,ND+1) = ZD(I)
      END DO
C
      CALL SOLVE ( ND, ND, WK, WK(1,ND+1), IWK )
C
C     Desired coefficients overwrite the RHS.
C     Evaluate function at approximation points:
C
      DO I = 1, NI
         VAL = 0.
         DO J = 1, ND
            VAL = VAL + WK(J,ND+1) * HGK ( XI(I), YI(I), J )
         END DO
         ZI(I) = VAL
      END DO
C
      RETURN
      END
