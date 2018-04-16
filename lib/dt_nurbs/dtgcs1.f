      SUBROUTINE DTGCS1(NPT,X,NDIM,Y,NYDIM,YLHS,NLHS,YRHS,NRHS,NBDIM, 
     &                  NDEG,NDOM,NDEP,WORK,NWORK,C,IER) 
C 
C======================================================================= 
C 
C  PURPOSE: 
C 
C        This is the basic fitting routine for tensor product, complete 
C        spline interpolants. 
C 
C        Note:  The interpolating c-vector computed below is not 
C        ----   unique for the input data.  Additional conditions can 
C               be imposed on NDOM-2,...,0 dimensional hyperplanes 
C               by giving the desired values of the mixed higher order 
C               partials on these hyperplanes.  Currently all such 
C               values are assumed zero. 
C 
C     ****************************************************************** 
C 
      INTEGER MAXDOM 
      DOUBLE PRECISION  ZERO 
      PARAMETER  (MAXDOM=10, ZERO=0.D0) 
C 
C  Arguments: 
C 
      INTEGER  NPT(*), NDIM, NYDIM(*), NDEG(*), NDOM, NDEP, 
     &         NWORK, IER, NLHS(*), NRHS(*), NBDIM(*) 
      DOUBLE PRECISION  X(NDIM,*), Y(*), WORK(*), C(*), YLHS(*), 
     &                  YRHS(*) 
C 
C  Internal: 
C 
      INTEGER NCDIM(MAXDOM+1),IDX(MAXDOM+1),INCY(MAXDOM+1), 
     +    INCC(MAXDOM+1),NDAT(MAXDOM+1),II,IBC,JBC,IPROD 
      INTEGER IS, IV, IX, ICPT, IZ, NSYS, ICC, IPREAM, JREDUC 
      INTEGER IW1SIZ, IW2, IW2SIZ, IW3, IW3SIZ, KNOTST, IPRODC, IPRODY, 
     &        IERMAX, KORD, KHALF, NCIX, NKNOTS, NCOEF, NSYSBC, NCTOT, 
     &        ISKP1, ISKP2, ISKPOL, NBC, NDATA, NLOOP, ILHS, IRHS 
C 
C     ****************************************************************** 
C 
C             Compute pointers and work space sizes. 
C 
      C(2) = NDEP 
      NSYS = NDEP 
      NSYSBC = NDEP 
      IPREAM = 2 + 3 * NDOM 
      NCTOT = NDEP 
      IW1SIZ = 0 
      IW2SIZ = 0 
      IERMAX = 0 
      DO 110 IV = 1,NDOM 
        KORD = NDEG(IV) + 1 
        KHALF = (KORD+1)/2 
        NCOEF = (NPT(IV) - 2) + 2*KHALF 
        NCDIM(IV) = NCOEF 
        C(2+IV) = KORD 
        C(2+NDOM+IV) = NCOEF 
        C(2+2*NDOM+IV) = 0 
        NSYS = NSYS * NPT(IV) 
        NDATA = NPT(IV) + NLHS(IV) + NRHS(IV) 
        NDAT(IV) = NDATA 
        NSYSBC = NSYSBC*NDATA 
        NCTOT = NCTOT*NCOEF 
        IPREAM = IPREAM + KORD + NCOEF 
        IW1SIZ = MAX (IW1SIZ, NDATA) 
        IW2SIZ = MAX(IW2SIZ,5+2*NCOEF+KORD) 
  110 CONTINUE 
      NDAT(NDOM+1) = NDEP 
      NCDIM(NDOM+1) = NDEP 
C 
      INCC(1) = 1 
      INCY(1) = 1 
      IDX(1) = 0 
      DO 111 IV=2,NDOM+1 
        INCC(IV) = INCC(IV-1)*NCDIM(IV-1) 
        INCY(IV) = INCY(IV-1)*NYDIM(IV-1) 
        IDX(IV) = 0 
  111 CONTINUE 
      IW2 = IW1SIZ + 1 
      IW3 = IW2 + IW2SIZ 
      IW3SIZ = NWORK - IW1SIZ - IW2SIZ 
C 
C     ------------------------------------------------------------------ 
C             Copy the data into the C vector. 
C     ------------------------------------------------------------------ 
C 
      ICPT = IPREAM + 1 
      IYPT = 1 
      IBC = 0 
C     START LOOP (PASS THROUGH ALL INDICES OF NDOM+1 DIMENSIONAL ARRAYS 
C       Y AND C) 
C       NOTE THAT THE INDICES IN IDX() START AT ZERO 
  200   CONTINUE 
        IF (ICPT .GT. IPREAM+NCTOT) STOP 'INTERNAL CHECK 1 FAILED' 
        IF (IBC .EQ. 0) THEN 
C         AT INDICES BELONGING TO GRID POINT 
          C(ICPT) = Y(IYPT) 
        ELSE IF (IBC .GT. 0) THEN 
C         AT INDICES BELONGING TO A BOUNDARY CONDITION 
C         COMPUTE BOUNDARY PT INDEX ON NDOM-1 DIMENSIONAL HYPERPLANE 
C           EXCLUDING IBC DIRECTION 
          JBC = 0 
          IPROD = 1 
          DO 205 IV=1,NDOM 
            IF (IV.EQ.IBC) GO TO 205 
            JBC = JBC + IDX(IV)*IPROD 
            IPROD = IPROD * NPT(IV) 
  205     CONTINUE 
          IF (IDX(IBC) .LT. NPT(IBC) + NLHS(IBC)) THEN 
C           AT LEFT HAND BOUNDARY CONDITION 
            II = IDX(IBC) + 1 - NPT(IBC) 
            IS = II + NBDIM(1)*(JBC + NBDIM(2)*(IDX(NDOM+1) +  
     +          NBDIM(3)*(IBC-1))) 
            C(ICPT) = YLHS(IS) 
          ELSE 
C           AT RIGHT HAND BOUNDARY CONDITION 
            II = IDX(IBC) + 1 - NPT(IBC) - NLHS(IBC) 
            IS = II + NBDIM(1)*(JBC + NBDIM(2)*(IDX(NDOM+1) +  
     +          NBDIM(3)*(IBC-1))) 
            C(ICPT) = YRHS(IS) 
          END IF 
        ELSE 
C         AT INDICES BELONGING TO DOUBLE (OR HIGHER) BOUNDARY 
          C(ICPT) = ZERO 
        END IF 
C       MOVE TO NEXT INDICES, UPDATING ICPT AND IYPT TO MATCH 
        II = 1 
C       START LOOP (FIND INDEX WHICH INCREMENTS) 
  220     CONTINUE 
          IDX(II) = IDX(II) + 1 
          IF (IDX(II) .LT. NDAT(II)) GO TO 230 
          IYPT = IYPT - (IDX(II)-1)*INCY(II) 
          ICPT = ICPT - (IDX(II)-1)*INCC(II) 
          IDX(II) = 0 
          IF (II .LE. NDOM) THEN 
            IF (NPT(II) .LT. NDAT(II)) THEN 
C             LEAVING BOUNDARY INDEX ON II  
              IF (IBC .LT. 0) THEN 
                IBC = IBC + NDOM + 1 
              ELSE 
                IBC = 0 
              END IF 
            END IF 
          END IF 
          II = II + 1 
          IF (II .GT. NDOM+1) GO TO 290 
          GO TO 220 
C       END LOOP (FIND INDEX WHICH INCREMENTS) 
  230   CONTINUE 
        IYPT = IYPT + INCY(II) 
C         NOTE THAT IYPT IS INVALID WHILE AT BOUNDARY CONDITIONS 
        ICPT = ICPT + INCC(II) 
        IF (II .LE. NDOM) THEN 
          IF (IDX(II) .EQ. NPT(II)) THEN 
C           ENTERING BOUNDARY INDEX ON II 
            IF (IBC .EQ. 0) THEN 
              IBC = II 
            ELSE 
              IBC = IBC - NDOM - 1 
            END IF 
          END IF 
        END IF 
        GO TO 200 
C     END LOOP (PASS THROUGH ALL VALID INDICES) 
  290 CONTINUE 
C 
C     ------------------------------------------------------------------ 
C             Loop through each variable. 
C     ------------------------------------------------------------------ 
C 
      KNOTST = 3 + 3 * NDOM 
      NSYS = NSYSBC 
      DO 390 IV = 1,NDOM 
        ICC = -1 
        NDATA = NDAT(IV) 
        NSYS = NSYS/NDATA 
        NCOEF = NINT( C(2+NDOM+IV) ) 
        NKNOTS = NCOEF + NINT( C(2+IV) ) 
        IZ = IW2 + 5 + NKNOTS 
        ICPT = IPREAM + 1 
        ILHS = NPT(IV) + 1 
        IRHS = ILHS + NLHS(IV) 
C 
C----- Solve each of the univariate systems for this pass. 
C 
        NLOOP = NSYS/INCC(IV) 
C       START LOOP (PASS THROUGH ALL OUTER INDEX COMBINATIONS) 
  300     CONTINUE 
          DO 310 JS = 0,INCC(IV)-1 
            CALL DCOPY(NDATA,C(ICPT+JS),INCC(IV),WORK,1) 
            CALL DTCGEN(NDEG(IV),X(1,IV),NPT(IV),WORK,1,WORK(ILHS), 
     &                 NLHS(IV),WORK(IRHS),NRHS(IV),ICC,WORK(IW3), 
     &                 IW3SIZ,IW2SIZ,WORK(IW2),NCD,IER) 
            IF (IER .LT. 0) THEN 
              IER = -200 
              RETURN 
            ENDIF 
            IERMAX = MAX(IER,IERMAX) 
            CALL DCOPY(NCOEF,WORK(IZ),1,C(ICPT+JS),INCC(IV)) 
  310     CONTINUE 
          II = IV + 1 
C         START LOOP (FIND INDEX TO INCREMENT) 
  320       CONTINUE 
            IDX(II) = IDX(II) + 1 
            IF (IDX(II) .LT. NDAT(II)) GO TO 330 
            ICPT = ICPT - (IDX(II)-1)*INCC(II) 
            IDX(II) = 0 
            IF (II.EQ.NDOM+1) GO TO 380 
            II = II + 1 
            GO TO 320 
C         END LOOP (FIND INDEX TO INCREMENT) 
  330     CONTINUE 
          ICPT = ICPT + INCC(II) 
          GO TO 300 
C       END LOOP (PASS THROUGH ALL OUTER INDICES) 
  380   CONTINUE 
        NSYS = NSYS * NCOEF 
        CALL DCOPY(NKNOTS,WORK(IW2+5),1,C(KNOTST),1) 
        KNOTST = KNOTST + NKNOTS 
  390 CONTINUE 
      IF (NSYS .NE. NCTOT) STOP 'INTERNAL CHECK 2 FAILED' 
      IER = IERMAX 
C 
      RETURN 
      END
