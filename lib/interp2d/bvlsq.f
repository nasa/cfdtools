C+---------------------------------------------------------------------
C
      SUBROUTINE BVLSQ(F,M,XD,YD,ZD,ND,XI,YI,ZI,NI,C,WK,IER)
C
C  ACRONYM: BiVariate Least SQuares fit
C
C  PURPOSE: Fits data in least squares sense with user-supplied basis
C           functions, and evaluates result at specified points.
C
C  ARGUMENTS:
C       ARG    TYPE  I/O/S   DIM     DESCRIPTION
C
C       F      Real    I      -      Input bivariate basis functions.
C             Function               Function must be of the form F(I,X,Y)
C                                    where I specifies the basis function to
C                                    use, for I = 1:M.  For example:
C
C                                        FUNCTION FXY1 (I,X,Y)
C                                        GO TO (10,20,30) I
C                                    10  FXY1 = 1.
C                                        GO TO 40
C                                    20  FXY1 = X
C                                        GO TO 40
C                                    30  FXY1 = Y
C                                    40  RETURN
C                                        END
C
C                                    This is a bivariate polynomial of degree
C                                    1 (a plane).
C
C        M     I*4     I     -       Number of basis functions.
C       XD     R*4     I     ND      Input X coordinates of data.
C       YD     R*4     I     ND      Input Y coordinates of data.
C       ZD     R*4     I     ND      Input function values at point (X,Y).
C       ND     I*4     I     -       Number of data points.
C       XI     R*4     I     NI      Input X coordinates for evaluation.
C       YI     R*4     I     NI      Input Y coordinates for evaluation.
C       ZI     R*4     O     NI      Output evaluations of interpolant at
C                                    point (XI(I),YI(I)) for I=1:NI.
C       NI     I*4     I      -      Number of evaluation points.
C        C     R*4     O      M      Output coefficients of basis functions.
C                                    C(I) corresponds to F(I,X,Y).
C       WK     R*4     S      *      Work array should be dimensioned
C                                    (M+2)*ND where M is
C                                    the number of basis functions and
C                                    ND is the number of data points.
C      IER     I*4     O      -      Return error code.
C                                    129 = Incorrect arguments.
C                                          Either ND or M is <= 0.
C                                    130 = Unable to compute a solution
C                                          because the rank of the matrix
C                                          is less than the number of equality
C                                          contraints.
C                                    131 = Iterative refinement failed. A
C                                          solution has been obtained, but
C                                          it may be inaccurate.  Not used.
C
C  EXTERNAL REFERENCES: HDESOL (Linear least squares solution, LINSYS)
C
C  ENVIRONMENT: VAX 11/780; FORTRAN
C
C  KNOWN SYSTEM DEPENDENCIES: END DO is non-standard.
C
C  AUTHOR:   Douglas L. Richter    Informatics General Corporation
C
C  DEVELOPMENT HISTORY:
C       DATE   INITIALS  DESCRIPTION
C     03/21/83    DLR    Original design and coding
C     May, 1987   DAS    Tidied header up somewhat
C     10/06/97    DLH    Changed IMSL call LLBQF to non-IMSL HDESOL.  Removed
C			 CV.  Added error checking previously done by IMSL.
C
C----------------------------------------------------------------------
C
      EXTERNAL  F
      DIMENSION XD(1),YD(1),ZD(1),XI(1),YI(1),ZI(1),C(1),WK(1)
C
C     Check for bad degree or no data
C
      IF (M.LE.0 .OR. ND.LE.0) THEN
         IER=129
         RETURN
      ENDIF
C
C     Check for not enough data (rank too low)
C
      IF (ND.LT.M) THEN
         IER=130
         RETURN
      ENDIF
C
C     Set up ND x M system of equations
C
      IND=1
      DO I=1,M
         DO K=1,ND
            WK(IND)=F(I,XD(K),YD(K))
            IND=IND+1
         ENDDO
      ENDDO
C
C     Add the right-hand-side as column M+1
C
      DO K=1,ND
         WK(IND)=ZD(K)
         IND=IND+1
      ENDDO
C
C     Solve
C
      CALL HDESOL(ND,ND,M+1,WK(1),WK(IND),SSQMIN)
C
C     Copy the solution
C
      DO I=1,M
         C(I)=WK(IND+I-1)
      ENDDO
C
C     Evaluate the interpolant at XI,YI points
C
      IF(NI.EQ.0)RETURN
      DO I=1,NI
         ZI(I)=0.
      ENDDO
      IND=1
      DO I=1,M
         DO K=1,NI
            ZI(K)=ZI(K)+C(I)*F(I,XI(K),YI(K))
         ENDDO
      ENDDO
      RETURN
      END
