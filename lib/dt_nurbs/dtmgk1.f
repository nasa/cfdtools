      SUBROUTINE DTMGK1(C,ICPTR,NC,T,MLTT,NT,XKNOTS,MULTI,INDX,
     1           TOUT,NTOUT,NDIM,NDNEED)
      DOUBLE PRECISION C(*),T(*),XKNOTS(*),TOUT(*)
      INTEGER ICPTR(*),NC,MLTT(*),NT,MULTI(*),INDX(*),
     1           NTOUT,NDIM,NDNEED
C
C     TO SET UP AND INTERLACE SEVERAL KNOT VECTORS
C
C     INPUT    C         ARRAY OF SPLINE VECTORS
C              ICPTR     POINTER TO START OF EACH SPLINE VECTOR
C              NC        NUMBER OF SPLINE VECTORS
C              T         KNOTS TO BE ADDED
C              MLTT      MULTIPLICITIES OF T
C              NT        NUMBER OF ELEMENTS IN T AND MLT
C     WORKING
C     STORAGE  XKNOTS    WORK ARRAY
C              MULTI    WORK ARRAY
C              INDX     WORK ARRAY
C     OUTPUT   TOUT      MERGED KNOT VECTOR
C              NTOUT     NUMBER OF ELEMENTS IN TOUT
C
C     SUBROUTINES USED:
C              DT2SVY
C
      DOUBLE PRECISION A,B,TEMP
      INTEGER I,IKNOTS,KKNOTS,NKNOTS
C
C
      NDIM2 = (NDIM*NC + 1)/2
      A = C(ICPTR(1) + 5)
      I = INT(C(ICPTR(1)+2) + C(ICPTR(1)+3))
      B = C(ICPTR(1) + I + 4)
      NDNEED = INT(C(ICPTR(1)+3) + C(ICPTR(1)+2))
      IF (NDNEED .GT. NDIM2) GO TO 900
      INDX(1) = 1
      KORD = INT(C(ICPTR(1) + 2))
      DO 100 I = 1,NC
          CALL DT2SVY(C(ICPTR(I)),1,XKNOTS(INDX(I)),
     1            MULTI(INDX(I)),NKNOTS,IKNOTS,KKNOTS,IER)
          IF (I .NE. NC) THEN
              NDNEED = NDNEED+INT(C(ICPTR(I+1)+3)+C(ICPTR(I+1)+2))
              IF (NDNEED .GT. NDIM2) GO TO 900
              INDX(I+1) = INDX(I) + NKNOTS
          END IF
          INDX(I) = INDX(I) + 1
  100 CONTINUE
      IF (KORD .GT. NDIM) GO TO 900
      DO 200 I = 1,KORD
          TOUT(I) = A
  200 CONTINUE
      NTOUT = KORD
      IT = 1
C
      IF (NT .GE. IT) THEN
          IF (T(IT) .EQ. TOUT(NTOUT)) IT = 2
      END IF
C
  350 MLT = 0
      TEMP = B
      DO 300 I = 1,NC
          TEMP = MIN(TEMP,XKNOTS(INDX(I)))
  300 CONTINUE
      IF (NT .GE. IT) THEN
          IF (T(IT) .LE. TEMP) THEN
              TEMP = T(IT)
              MLT = MLTT(IT)
              IT = IT + 1
          END IF
      END IF
      IF (TEMP .LT. B) THEN
          DO 320 I = 1,NC
             IF (XKNOTS(INDX(I)) .EQ. TEMP) THEN
                 MLT = MAX(MLT,MULTI(INDX(I)))
                 INDX(I) = INDX(I) + 1
             END IF
  320    CONTINUE
         DO 340 I = 1,MLT
             NTOUT = NTOUT+1
             IF (NTOUT .GT. NDIM) GO TO 900
             TOUT(NTOUT) = TEMP
  340    CONTINUE
      GO TO 350
      END IF
      DO 400 I = 1,KORD
          NTOUT = NTOUT + 1
          IF (NTOUT .GT. NDIM) GO TO 900
          TOUT(NTOUT) = B
  400 CONTINUE
  900 RETURN
      END

