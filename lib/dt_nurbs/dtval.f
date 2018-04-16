        FUNCTION DTVAL(X,CC,WORK,NWORK)
        DOUBLE PRECISION X, DTVAL
        DOUBLE PRECISION CC(*), WORK(*), V(10, 2)
C  **
C  **           ARC-LENGTH FUNCTION BOX
C  **
        NDIMV = 10
        NC = 5 + CC(4) + CC(3) + CC(2)*CC(4)
        M = CC(2)
        IF(M .LT. 0)M = -M - 1
        NDER = 1
        CALL DTSPDR( X, NDER, CC, WORK, NWORK, V, NDIMV, IER)
        DTVAL = 0.
        DO 10 I = 1, M
        DTVAL = DTVAL + V(I, 2)**2
10      CONTINUE
        DTVAL = SQRT( DTVAL )
        RETURN
        END
 
