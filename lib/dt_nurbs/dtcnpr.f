        SUBROUTINE DTCNPR( S, T, IOPT, WORK, NWORK, C, IER)
C===============================================================================
C  PURPOSE:  EXTRACT CONSTANT PARAMETER CURVES FROM A SURFACE
C
C  USAGE:
C
C  INPUTS:  
C  
C             S       SPLINE ARRAY FOR THE SURFACE
C  
C             T       CONSTANT PARAMETER VALUE
C  
C             IOPT    OPTION FLAG
C                     = 1     FIRST PARAMETER IS FIXED
C                     = 2     SECOND PARAMETER IS FIXED
C  
C
C
C
C  WORKING STORAGE:
C            WORK    WORK ARRAY OF LENGTH NWORK
C  
C             NWORK   LENGTH OF WORK ARRAY, NWORK .GE. K ** 2 + NCOEF
C                     WHERE K = MAX( S(3), S(4)) AND 
C                     NCOEF = MAX( S(5), S(6))
C  
C
C
C  OUTPUT:  
C  
C             C       SPLINE ARRAY FOR CONSTANT PARAMETER CURVE
C  
C             IER     SUCCESS/ERROR FLAG
C                     =  0   COMPUTATION WAS SUCCESSFUL
C                            ERRORS BETWEEN -1 AND -9 ARE RETURNED
C                            FROM DTBSPL
C                     = -10  IOPT .NE. 1 OR 2
C                     = -11  NUMBER OF INDEPENDENT VARIABLES .NE. 2
C                     = -12  NUMBER OF DEPENDENT VARIABLES .EQ. 0
C                     = -13  T IS OUTSIDE THE PARAMETER RANGE
C                     = -14  WORK IS TOO SMALL
C
C  ROUTINES CALLED:  DTBSPL
C
C  PROGRAMMER:
C
C  CREATION DATE: 
C
C  UPDATE DATE:  JUNE 10, 1993               UPDATE VERSION #: 2.02
C  UPDATE NOTICE: BUG FIX: - D. R. FERGUSON
C                Fixed and indexing problem.
C
C===============================================================================
C  **           I / O VARIABLE TYPES
C  **
        DOUBLE PRECISION S(*), WORK(*), C(*)
        DOUBLE PRECISION T
        INTEGER IOPT, NWORK, IER
C  **
C  **           INTERNAL VARIABLE TYPES
C  **
        INTEGER NKU, NKV, N1, N2, NOWS, NOWC, MABS
        INTEGER I, IV, J, KU, KV, IL, ND, NEED
C  **
C  **           ERROR CHECKING
C  **
        IER = -10
        IF( IOPT .LT. 1 .OR. IOPT .GT. 2)RETURN
        IER = -11
        IF (S(1) .NE. 2.D0) RETURN
        IER = -12
        IF (S(2) .EQ. 0.D0) RETURN
        KU  = S(3)
        KV  = S(4)
        N1  = S(5)
        N2  = S(6)
        NKU = N1 + KU
        NKV = N2 + KV
        IER = -13
        IF( IOPT .EQ. 1 .AND. (S(9) .GT. T .OR. S(8 + NKU) .LT. T))
     *      RETURN
        IF( IOPT .EQ. 2 .AND. (S(9 + NKU) .GT. T 
     *      .OR. S(8 + NKU + NKV) .LT. T))
     *      RETURN
        IER = -14
        NEED = MAX(N1, N2) + MAX(KU, KV)
        IF( NEED .GT. NWORK) RETURN
        IER = 0
C  **
C  **           EVALUATE B-SPLINES IN U (OR V) DEPENDING ON IOPT
C  **
        IL = 1
        ND = 0
        IF( IOPT .EQ. 1 ) 
     *       CALL DTBSPL(S(9), NKU, T, IL, KU, ND, N1, WORK, NWORK
     *       , WORK(1 + KU ** 2), IER)
        IF( IOPT .EQ. 2 ) 
     *       CALL DTBSPL(S(9 + NKU), NKV, T, IL, KV, ND, N2, WORK, NWORK
     *       , WORK(1 + KV ** 2), IER)
        IF( IER .NE. 0 ) RETURN
        IF( IOPT .EQ. 2 )GO TO 10
C  **
C  **           BEGIN CONSTRUCTION OF C FOR CONSTANT 
C  **           FIRST PARAMETER VALUE
C  **
        C(1) = 1.D0
        C(2) = S(2)
        C(3) = S(4)
        C(4) = S(6)
        C(5) = S(8)
C  **
C  **           SET KNOTS
C  **
        NOWS = 8 + NKU
        DO 20 I = 1, NKV
                NOWS = NOWS + 1
                C(5 + I) = S(NOWS)
20      CONTINUE
C  **
C  **           COMPUTE THE COEFFICIENTS
C  **
        NOWC = 5 + NKV
        NOWS = 8 + NKV + NKU
        MABS = ABS( S(2) )
        DO 30 IV = 1, MABS
                DO 40 J = 1, N2
                        NOWC = NOWC + 1
                        C(NOWC) = 0.D0
                        DO 50 I = 1, N1
                                NOWS = NOWS + 1
                                C(NOWC) = C(NOWC) + 
     *                                    S(NOWS) * WORK(I + KU ** 2)
     
50                      CONTINUE
40              CONTINUE
30      CONTINUE
        RETURN
10      CONTINUE
C  **
C  **           BEGIN CONSTRUCTION OF C FOR CONSTANT 
C  **           FIRST PARAMETER VALUE
C  **
        C(1) = 1.D0
        C(2) = S(2)
        C(3) = S(3)
        C(4) = S(5)
        C(5) = S(7)
C  **
C  **           SET KNOTS
C  **
        NOWS = 8
        DO 60 I = 1, NKU
                NOWS = NOWS + 1
                C(5 + I) = S(NOWS)
60      CONTINUE
C  **
C  **           COMPUTE THE COEFFICIENTS
C  **
        NOWC = 5 + NKU
        NOWS = 8 + NKV + NKU
        MABS = ABS( S(2) )
        DO 70 IV = 1, MABS
                DO 80 I = 1, N1
                        NOWC = NOWC + 1
                        C(NOWC) = 0.D0
                        DO 90 J = 1, N2
                           C(NOWC) = C(NOWC) 
     *                              + s(nows+(iv-1)*n1*n2
     *                              + (j-1)*n1+i) * WORK(J + KV ** 2)
90                      CONTINUE
80              CONTINUE
c        NOWS = NOWS + N1 * N2
70      CONTINUE
        RETURN
        END
