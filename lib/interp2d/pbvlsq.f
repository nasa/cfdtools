C+---------------------------------------------------------------------
C
      SUBROUTINE PBVLSQ (M,ND,NI,XD,YD,ZD,XI,YI,ZI,C,WK,IER)
C
C  ACRONYM: Polynomial BiVariate Least SQuares
C
C  PURPOSE: Fits a bivariate least squares polynomial of degree M
C           to the points (XD(I),YD(I),ZD(I)) and evaluates it at
C           the points (XI(I),YI(I)).
C
C  ARGUMENTS:
C       ARG    TYPE  I/O/S   DIM     DESCRIPTION
C        M     I*4     I      -      Degree of polynomial fit desired.
C        ND    I*4     I      -      Number of data points.
C        NI    I*4     I      -      Number of points at which to evaluate.
C        XD    R*4     I      ND     Input X locations of data points.
C        YD    R*4     I      ND     Input Y locations of data points.
C        ZD    R*4     I      ND     Input function values at data points.
C        XI    R*4     I      NI     X location of points at which to
C                                    evaluate.
C        YI    R*4     I      NI     Y location of points at which to
C                                    evaluate.
C        ZI    R*4     O      NI     Approximated value at points (XI(I),YI(I))
C
C        C     R*4     O      MT     The coefficients of each monomial in
C                                    the order of all terms involved with
C                                    increasing powers of X. For example, for
C                                    a polynomial of degree 2, the terms
C                                    1,Y,Y**2,X,XY,X**2 would be associated
C                                    with C(1) through C(6) in that order.
C
C        WK    R*4     S             A work array for the routine. The
C                                    dimension of WK is as follows:
C                                    Let M be the degree of the polynomial.
C                                    Then MT, the number of terms in the
C                                    polynomial is (M+1)*(M+2)/2. Then WK
C                                    should be dimensioned at least
C                                    (MT+2)*ND.
C
C        IER   I*4      O     -      Error return code.
C                                    129 = Incorrect arguments.
C                                          Either ND or MT is <= 0.
C
C                                    130 = Unable to compute a solution
C                                          because the rank of the matrix
C                                          is less than the number of equality
C                                          contraints.
C
C                                    131 = Iterative refinement failed. A
C                                          solution has been obtained, but
C                                          it may be inaccurate.  Not used.
C
C  LOCAL VARIABLES:
C       VAR    TYPE    DIM     DESCRIPTION
C       MT     I*4      -      Column dimension of linear system.
C       IND    I*4      -      Keeps track of work array space.
C
C  ENVIRONMENT: VAX 11/780; FORTRAN
C
C  EXTERNAL REFERENCES:
C     PBASF  Generates bivariate polynomial terms
C     HDESOL Linear least squares (LINSYS)
C
C  AUTHOR:   Douglas L. Richter    Informatics General Corporation
C
C  DEVELOPMENT HISTORY:
C       DATE   INITIALS  DESCRIPTION
C
C     03/29/83    DLR    Original design and coding
C     May, 1987   DAS    Tidied header up somewhat
C     10/06/97    DLH    Changed IMSL call LLBQF to non-IMSL HDESOL.  Removed
C			 CV.  Added error checking previously done by IMSL.
C
C----------------------------------------------------------------------
C
      DIMENSION XD(1),YD(1),ZD(1),XI(1),YI(1),ZI(1),C(1),WK(1)
C
C     Check for bad degree or no data
C
      IF (M.LE.0 .OR. ND.LE.0) THEN
         IER=129
         RETURN
      ENDIF
C
C     Initialize column dimension and work array counter
C
      MT = M*(M+3)/2+1
      IND=1
C
C     Check for not enough data (rank too low)
C
      IF (ND.LT.MT) THEN
         IER=130
         RETURN
      ENDIF
C
C     Form linear ND x MT system
C
      DO I=0,M
         DO J=0,M-I
            DO K=1,ND
               WK(IND)=PBASF(I,J,XD(K),YD(K))
               IND=IND+1
            ENDDO
         ENDDO
      ENDDO
C
C     Add the right-hand-side as column MT+1
C
      DO K=1,ND
         WK(IND)=ZD(K)
         IND=IND+1
      ENDDO
C
C     Solve
C
      CALL HDESOL(ND,ND,MT+1,WK(1),WK(IND),SSQMIN)
C
C     Copy the solution
C
      DO I=1,MT
         C(I)=WK(IND+I-1)
      ENDDO
C
C     Evaluate the polynomial at XI,YI points
C
      IF(NI.EQ.0)RETURN
      DO I=1,NI
         ZI(I)=0.
      ENDDO
      IND=1
      DO I=0,M
         DO J=0,M-I
            DO K=1,NI
               ZI(K)=ZI(K)+C(IND)*PBASF(I,J,XI(K),YI(K))
            ENDDO
            IND=IND+1
         ENDDO
      ENDDO
      RETURN
      END
C
C     PBASF returns the value of the monomial X**I * Y**J
C
      FUNCTION PBASF(I,J,X,Y)
C
C     Calculate X**I
C
      XV=1.
      DO K=1,I
         XV=X*XV
      ENDDO
C
C     Calculate Y**J
C
      YV=1.
      DO K=1,J
         YV=Y*YV
      ENDDO
C
      PBASF=XV*YV
C
      RETURN
      END
