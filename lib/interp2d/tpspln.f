C+---------------------------------------------------------------------
C
      SUBROUTINE TPSPLN ( XD,YD,ZD,ND,XI,YI,ZI,NI,IWK,WK,IER,COND )
C
C  PURPOSE: TPSPLN ("thin-plate spline") does bivariate interpolation
C           for irregular data by a method analogous to Hardy's  (see
C           subroutine HARDY) except that the basis functions are of
C           the form
C
C                PHI(K,X,Y) = D(K)**2 * LOG ( D(K) )   where
C
C              D(K)**2 = ( X - XD(K) )**2  +  ( Y - YD(K) )**2.
C
C           This function is derived from the differential equation
C           for small deflections in a thin plate.   (Actually, the
C           function used is double the one above,  since advantage
C           is taken of the fact that LOG (F**0.5) = 0.5 LOG (F).)
C
C           The data can then be represented by the function
C
C                    G(X,Y) = SUM (K=1:N) A(K)*PHI(K,X,Y)
C
C           where A(K) is the Kth linear parameter computed here.
C
C           This function is evaluated at the given point(s).
C
C  ARGUMENTS:
C       ARG    TYPE  I/O/S   DIM     DESCRIPTION
C        XD    R*4     I      ND     Input X locations of data points.
C        YD    R*4     I      ND     Input Y locations of data points.
C        ZD    R*4     I      ND     Input function values at data points.
C        ND     I      I      -      Number of data points.
C        XI    R*4     I      NI     X locations of points at which to 
C                                    approximate.
C        YI    R*4     I      NI     Y locations of points at which to 
C                                    approximate.
C        ZI    R*4     O      NI     Approximated value at points (XI(I),YI(I))
C        NI     I      I      -      Number of interpolation points.
C        IWK   I*4     S      ND     Pivot vector for solving square system.
C        WK    R*4     S   ND*(ND+1) For setting up and solving square system.
C        IER    I      O      -      Error return code. 0 means no error.
C                                    129 means matrix is singular. The 
C                                    probable cause is that two of the
C                                    data points are the same.
C        COND   R*4    O      -      Condition number of matrix involved.
C
C  COMMONS USED: None
C
C  FILES USED: None
C
C  EXTERNAL REFERENCES: DECOMP, SOLVE
C
C  ENVIRONMENT: VAX 11/780; FORTRAN
C
C  KNOWN SYSTEM DEPENDENCIES:
C     END DO is a VAX extension.
C
C  AUTHOR:   Douglas L. Richter   Informatics General Corporation
C
C  DEVELOPMENT HISTORY:
C       DATE   INITIALS  DESCRIPTION
C     04/26/83    DLR    Original design and coding.
C     12/28/83    DAS    Expanded as for subroutine HARDY.
C
C-----------------------------------------------------------------------
C
      DIMENSION XD(ND), YD(ND), ZD(ND), XI(NI), YI(NI), ZI(NI), IWK(ND),
     1          WK(ND,ND+1)
C
C     Statement function for square of distance...
C
      DSQ ( X, Y, K ) = ( X - XD(K) )**2 + ( Y - YD(K) )**2
C
C
C     Set up the left hand side matrix (symmetric):
C
      DO I = 1, ND
         DO J = 1, I
            DIJSQ   = DSQ ( XD(J), YD(J), I )
	    IF ( DIJSQ .GT. 0.E+0 ) THEN
               WK(I,J) = DIJSQ * ALOG ( DIJSQ )
            ELSE
               WK(I,J) = 0.E+0
            END IF
            WK(J,I) = WK(I,J)
         END DO
      END DO
C
C     Factorize the matrix:
C
      CALL DECOMP ( ND, ND, WK, COND, IWK, WK(1,ND+1) )
C
      IER = 0
      IF ( COND .GE. 1.E+32 ) THEN
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
         VAL = 0.E+0
         DO J = 1, ND
            DIJSQ = DSQ ( XI(I), YI(I), J )
            IF ( DIJSQ .GT. 0.E+0 ) THEN
               VAL = VAL + WK(J,ND+1) * DIJSQ * ALOG ( DIJSQ )
            END IF
         END DO
         ZI(I) = VAL
      END DO
C
      RETURN
      END
