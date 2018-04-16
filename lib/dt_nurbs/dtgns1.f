       SUBROUTINE DTGNS1(NPT,X,NDIM,Y,NYDIM,NDEG,NDOM,NDEP, 
     &                  WORK,NWORK,C,IER) 
C 
C======================================================================= 
C 
C  PURPOSE: 
C 
C----- This is the basic fitting routine for tensor product, natural  
C      spline interpolants on a grid. 
C 
C 
C     ****************************************************************** 
      INTEGER MAXDOM 
      PARAMETER (MAXDOM=10) 
C 
C  Arguments: 
C 
      INTEGER  NPT(*), NDIM, NYDIM(*), NDEG(*), NDOM, NDEP, 
     &         NWORK, IER 
      DOUBLE PRECISION  X(NDIM,*), Y(*), WORK(*), C(*) 
C 
C  Internal: 
C 
      INTEGER NCDIM(MAXDOM+1), IDX(MAXDOM+1), INCC(MAXDOM+1),  
     +    INCY(MAXDOM+1), NDAT(MAXDOM+1), II 
      INTEGER IS, IV, IX, ICPT, IZ, NSYS, ICC, IPREAM 
      INTEGER IW1SIZ, IW2, IW2SIZ, IW3, IW3SIZ, KNOTST, IPRODC, 
     &        IERMAX, KORD, KHALF, NCIX, NKNOTS 
      DOUBLE PRECISION  BCLHS,BCRHS 
C 
C     ****************************************************************** 
C 
C             Compute pointers and work space sizes. 
C 
      IF (NDOM .GT. MAXDOM) THEN 
C       THE CHANCE OF THIS LIMITATION BEING EXCEEDED BY A LEGITIMATE 
C       DATA SET IS NEGLIGIBLE.  SHOULD THAT EVER HAPPEN, THE REMEDY IS 
C       TO INCREASE MAXDOM AND RECOMPILE. 
        IER = -201 
        RETURN 
      END IF 
      C(2) = NDEP 
      NSYS = NDEP 
      IPREAM = 2 + 3 * NDOM 
      IW1SIZ = 0 
      IW2SIZ = 0 
      IERMAX = 0 
      DO 110 IV = 1,NDOM 
        KORD = NDEG(IV) + 1 
        KHALF = (KORD+1)/2 
        NCOEF = (NPT(IV) - 2) + 2*KHALF 
        C(2+IV) = KORD 
        C(2+NDOM+IV) = NCOEF 
        C(2+2*NDOM+IV) = 0 
        NSYS = NSYS * NPT(IV) 
        NCDIM(IV) = NCOEF 
        NDAT(IV) = NPT(IV) 
        IPREAM = IPREAM + KORD + NCOEF 
        IW1SIZ = MAX (IW1SIZ, NPT(IV)) 
        IW2SIZ = MAX(IW2SIZ,5+2*NCOEF+KORD) 
110   CONTINUE 
C 
      NCDIM(NDOM+1) = NDEP 
      NDAT(NDOM+1) = NDEP 
      IDX(1) = 0 
      INCC(1) = 1 
      INCY(1) = 1 
      DO 111 IV=2,NDOM+1 
        IDX(IV) = 0 
        INCC(IV) = INCC(IV-1)*NCDIM(IV-1) 
        INCY(IV) = INCY(IV-1)*NYDIM(IV-1) 
111   CONTINUE 
C 
      IW2 = IW1SIZ + 1 
      IW3 = IW2 + IW2SIZ 
      IW3SIZ = NWORK - IW1SIZ - IW2SIZ 
C 
C     ------------------------------------------------------------------ 
C     Interpolate raw data in direction of first parameter. 
C     ------------------------------------------------------------------ 
C 
      NKNOTS = NINT(C(3)) + NCDIM(1) 
      IZ = IW2 + 5 + NKNOTS 
      ICPT = IPREAM+1 
      IYPT = 1 
      ICC = -1 
      NLOOP = NSYS/NDAT(1) 
C     START LOOP (Take each sequence of grid data in the first parameter 
C     direction, compute B-spline coefficients, save in C array) 
200     CONTINUE 
        IF (NLOOP .LT. 0) STOP 'DTGNS1 INTERNAL CHECK 1 FAILED' 
        NLOOP = NLOOP - 1 
        CALL DCOPY (NDAT(1), Y(IYPT), 1, WORK, 1) 
        CALL DTCGEN (NDEG(1), X(1,1), NDAT(1), WORK, 1, BCLHS, 0,  
     +      BCRHS, 0, ICC, WORK(IW3), IW3SIZ, IW2SIZ, WORK(IW2), 
     +      NCD, IER) 
        IF (IER .LT. 0) THEN 
          IER = -200 
          RETURN 
        END IF 
        IERMAX = MAX(IER, IERMAX) 
        CALL DCOPY (NCDIM(1), WORK(IZ), 1, C(ICPT), 1) 
C       Increment to next sequence of grid data 
        II = 2 
C       START LOOP (Find indices and location in simulated  
C       NDOM+1 dimensional arrays Y and C) 
220       CONTINUE 
          IDX(II) = IDX(II) + 1 
          IF (IDX(II) .LT. NDAT(II)) GO TO 230 
C         Finished II-th dimension, go to next dimension 
          IYPT = IYPT - (IDX(II)-1)*INCY(II) 
          ICPT = ICPT - (IDX(II)-1)*INCC(II) 
          IDX(II) = 0 
          IF (II .EQ. NDOM+1) GO TO 290 
          II = II + 1 
          GO TO 220 
C       END LOOP (find indices and location ...) 
230     CONTINUE 
        IYPT = IYPT + INCY(II) 
        ICPT = ICPT + INCC(II) 
        GO TO 200 
C     END LOOP (each sequence in first parameter direction) 
290   CONTINUE 
C     Save knots of first parameter 
      KNOTPT = 3 + 3*NDOM 
      CALL DCOPY( NKNOTS, WORK(IW2+5), 1, C(KNOTPT), 1) 
      KNOTPT = KNOTPT + NKNOTS 
C     Adjust total count of coefficients at this stage 
      NSYS = NSYS/NDAT(1)*NCDIM(1) 
C 
C     ------------------------------------------------------------------ 
C     Repeat this process in each remaining parameter direction. 
C     ------------------------------------------------------------------ 
C 
      DO 390 IV = 2,NDOM 
        ICC = -1 
        NSYS = NSYS/NDAT(IV) 
        NKNOTS = NCDIM(IV) + NINT( C(2+IV) ) 
        IZ = IW2 + 5 + NKNOTS 
        ICPT = IPREAM + 1 
C 
C----- Solve each of the univariate systems for this pass. 
C 
        NLOOP = NSYS/INCC(IV) 
C       START LOOP (Take each full block of B-spline coefficients in the 
C       IV-th parameter direction, compute the next level B-spline coef- 
C       ficients, and save them in the C array) 
300       CONTINUE 
          IF (NLOOP .LT. 0) STOP 'DTGNS1 INTERNAL CHECK 2 FAILED' 
          NLOOP = NLOOP - 1 
          DO 310 JS = 0,INCC(IV)-1 
            CALL DCOPY(NDAT(IV),C(ICPT+JS),INCC(IV),WORK,1) 
            CALL DTCGEN(NDEG(IV),X(1,IV),NDAT(IV),WORK,1,BCLHS, 
     &                  0,BCRHS,0,ICC,WORK(IW3),IW3SIZ, 
     &                  IW2SIZ,WORK(IW2),NCD,IER) 
            IF (IER .LT. 0) THEN 
              IER = -200 
              RETURN 
            ENDIF 
            IERMAX = MAX (IER, IERMAX) 
            CALL DCOPY(NCDIM(IV),WORK(IZ),1,C(ICPT+JS),INCC(IV)) 
310       CONTINUE 
C         Increment to next block's location 
          II = IV + 1 
C         START LOOP (Increment indices and locate next block) 
320         CONTINUE 
            IDX(II) = IDX(II) + 1 
            IF (IDX(II) .LT. NDAT(II)) GO TO 330 
C           Finished II-th dimension, go to next dimension 
            ICPT = ICPT - (IDX(II)-1)*INCC(II) 
            IDX(II) = 0 
            IF (II .EQ. NDOM+1) GO TO 380 
            II = II + 1 
            GO TO 320 
C         END LOOP (find indices and location ...) 
330       CONTINUE 
          ICPT = ICPT + INCC(II) 
          GO TO 300 
C       END LOOP (each block in IV-th parameter direction) 
380     CONTINUE 
C       Adjust total count of coefficients at this stage 
        NSYS = NSYS * NCDIM(IV) 
C       Save knots of IV-th parameter 
        CALL DCOPY(NKNOTS,WORK(IW2+5),1,C(KNOTPT),1) 
        KNOTPT = KNOTPT + NKNOTS 
390   CONTINUE 
      IER = IERMAX 
C 
      RETURN 
      END
