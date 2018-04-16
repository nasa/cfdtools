      SUBROUTINE D2IGCY (CMEM, IMEM, DMEM, IGE, LENC, C, IER)
      
C   Convert an IGES Entity type 126 (B-Spline Curve) or 128 (B-Spline
C   Surface) to a standard Library Spline Array (C).
C
C   DYNAMIC MEMORY INPUT-OUTPUT:
C     CMEM     Dynamically managed character data
C     IMEM     Dynamically managed integer data
C     DMEM     Dynamically managed double precision data
C
C   INPUT:
C     IGE      Pointer to the IGES Spline Entity (type 126 or 128)
C     LENC     Maximum defined length of the C array.
C
C   OUTPUT:
C     C        Converted spline.
C     IER      Returned error value.  Zero implies no errors.  
C              -1  = Dynamic memory is corrupt or uninitialized
C              -2  = IGE is a null pointer
C              -3  = Garbage value found in IGE
C              -4  = Ambiguous--garbage pointer or deleted entity
C              -5  = IGE points to a deleted entity
C              -6  = IGE does not point to an IGES Entity
C              -7  = IGE is not a B-Spline (type 126 or 128)
C              -8  = Data in Entity is inconsistent with IGES Standard
C              -9  = LENC insufficient
C
C
C   CALLS:
C     D2IGFI
C     D2IPFI
C     D2IGFD
C     D2IPFD
C     D0PTR
C     DTERPT
C     DTERR
C
C   HISTORY:
C     9/15/93  D. Parsons  Created from HSRIGS.
C
C*****
C
C   Long Name Alias:
C     ENTRY D2_IGES_TO_CARRAY (CMEM, IMEM, DMEM, IGE, LENC, C, IER)

      CHARACTER CMEM*(*)
      INTEGER IMEM(*)
      DOUBLE PRECISION DMEM(*)
      
      INTEGER     ENTGE
      PARAMETER  (ENTGE = 238)

      INTEGER  IGE, LENC, IER
      DOUBLE PRECISION C(*), DVAL
      INTEGER  ITYPE, IHDR, ITYP, IERX
      INTEGER  IDATA(9)
       
      INTEGER  IP(3), IQ(4), I, J, LN, M
      INTEGER  INDEX1, INDEX2, IDEG1, IDEG2, K1, K2, N1, N2
      INTEGER  IA, IB, IC, IN1, IN2
      INTEGER  IPROP3, R, NEED, NUM, INCDA, JPAR, INCPAR 
      
      CHARACTER*6    SUBNAM
      DATA SUBNAM   /'D2IGCY'/
      
C****

      IER = 0
      IERX = 0
      
C     Check pointer and type

      CALL D0PTR (CMEM, IMEM, DMEM, IGE, ENTGE, 0, IHDR, ITYP, IER)
      
      IF (IER .NE. 0) GOTO 9900
      
      CALL D2IGFI (CMEM, IMEM, DMEM, IGE, 'D', 1, ITYPE, IERX)
      
      IF (ITYPE .EQ. 126) THEN
      
C        Convert B-Spline Curve Entity

         CALL DTERPT(0)
         CALL D2IPFI (CMEM, IMEM, DMEM, 6, IGE, 1, 1, 1, IDATA, IERX)
         CALL DTERPT(1)
         IF (IERX .NE. 0) THEN
            IER = -8
            GOTO 9900
         ENDIF
         INDEX1 = IDATA(1)
         IDEG1  = IDATA(2)
         K1     = IDEG1 + 1
         N1     = INDEX1 + 1
         LN     = 1
         IPROP3 = IDATA(5)
         
C        Load number of Independent/Dependent variables,
C             order of the splines,
C             number of B-Splines coefficients

         C(1) = LN
         
         IF (IPROP3 .EQ. 0) THEN
            M = 4
            C(2) = -M
         ELSE
            M = 3
            C(2) = M
         ENDIF

         C(3) = K1
         C(4) = N1
         C(5) = 0

C        Calculate the necessary length of the C array

         NEED = 3*LN+2
         NEED = NEED + C(3) + C(4)
         R = C(4)
         NEED = NEED + M * R
            
         IF (NEED .GT. LENC) THEN
            IER = -9
            GOTO 9900
         ENDIF
         
C        Set up C-array pointers
C
         IP(1) = 5
         IP(2) = IP(1) + C(4) + C(3)
         
         IQ(1) = IP(2)
         
         DO 110 J = 2, M
            IQ(J) = IQ(J-1) + C(4)
  110    CONTINUE
C
C        Set up IGES Array pointer
C
         IA = INDEX1 + IDEG1 + 1
C
C        Load C array
C
C        --  knots  --
C
         NUM    = IA+1
         JPAR   = 7
         INCPAR = 1
         INCDA  = 1
         CALL DTERPT(0)
         CALL D2IPFD (CMEM, IMEM, DMEM, NUM, IGE, JPAR, INCPAR, INCDA, 
     +      C(IP(1)+1), IERX)
         CALL DTERPT(1)
         IF (IERX .NE. 0) THEN
            IER = -8
            GOTO 9900
         ENDIF
         
C        --  control points --

         IF (IPROP3 .EQ. 0) THEN
C           Rational
            DO 130 I = 1, 3
               NUM    = INDEX1+1
               JPAR   = 8+IA+INDEX1+I
               INCPAR = 3
               INCDA  = 1
               CALL DTERPT(0)
               CALL D2IPFD (CMEM, IMEM, DMEM, NUM, IGE, JPAR, INCPAR, 
     +            INCDA, C(IQ(I)+1), IERX)
               CALL DTERPT(1)
               IF (IERX .NE. 0) THEN
                  IER = -8
                  GOTO 9900
               ENDIF
         
               DO 130 J = 0, INDEX1
                  CALL DTERPT(0)
                  CALL D2IGFD (CMEM, IMEM, DMEM, IGE, 'P', 8+IA+J, 
     +                  DVAL, IERX)
                  CALL DTERPT(1) 
                  IF (IERX .NE. 0) THEN
                     IER = -8
                     GOTO 9900
                  ENDIF
                  C(IQ(I)+1+J) = C(IQ(I)+1+J)*DVAL
  130       CONTINUE
   
            NUM    = INDEX1+1
            JPAR   = 8+IA
            INCPAR = 1
            INCDA  = 1
            CALL DTERPT(0)
            CALL D2IPFD (CMEM, IMEM, DMEM, NUM, IGE, JPAR, INCPAR, 
     +         INCDA, C(IQ(4)+1), IERX)
            CALL DTERPT(1)
            IF (IERX .NE. 0) THEN
               IER = -8
               GOTO 9900
            ENDIF
         ELSE
C           Non-Rational
            DO 150 I = 1, 3
               NUM    = INDEX1+1
               JPAR   = 8+IA+INDEX1+I
               INCPAR = 3
               INCDA  = 1
               CALL DTERPT(0)
               CALL D2IPFD (CMEM, IMEM, DMEM, NUM, IGE, JPAR, INCPAR, 
     +            INCDA, C(IQ(I)+1), IERX)
               CALL DTERPT(1)
               IF (IERX .NE. 0) THEN
                  IER = -8
                  GOTO 9900
               ENDIF
  150       CONTINUE
         ENDIF

      ELSE IF (ITYPE .EQ. 128) THEN
      
C        Convert B-Spline Surface Entity

         CALL DTERPT(0)
         CALL D2IPFI (CMEM, IMEM, DMEM, 9, IGE, 1, 1, 1, IDATA, IERX)
         CALL DTERPT(1)
         IF (IERX .NE. 0) THEN
            IER = -8
            GOTO 9900
         ENDIF

         INDEX1 = IDATA(1)
         INDEX2 = IDATA(2)
         IDEG1  = IDATA(3)
         IDEG2  = IDATA(4)
         K1     = IDEG1 + 1
         K2     = IDEG2 + 1
         N1     = INDEX1 + 1
         N2     = INDEX2 + 1
         LN     = 2
         IPROP3 = IDATA(7)

C        Load number of Independent/Dependent variables,
C             order of the splines,
C             number of B-Splines coefficients

         C(1) = LN
            
         IF (IPROP3 .EQ. 0) THEN
            M = 4
            C(2) = -M
         ELSE
            M = 3
            C(2) = M
         ENDIF
            
         C(3) = K1
         C(4) = K2
         C(5) = N1
         C(6) = N2
         C(7) = 0
         C(8) = 0

C        Calculate the necessary length of the C array

         NEED = 3*LN+2
      
         DO 210 I = 2, LN+1
            NEED = NEED + C(1+I) + C(LN+1+I)
  210    CONTINUE
   
         R = 1
         DO 220 J = LN+3, 2*LN+2
            R = R * C(J)
  220    CONTINUE 
         
         NEED = NEED + M * R
            
         IF (NEED .GT. LENC) THEN
            IER = -9
            GOTO 9900
         ENDIF
         
C        Set up C-Array Pointers
C
         IP(1) = 8
         IP(2) = IP(1) + C(3) + C(5)
         IP(3) = IP(2) + C(4) + C(6) 
         
         IQ(1) = IP(3)
         
         DO 230 J = 2, M
            IQ(J) = IQ(J-1) + C(5)*C(6)
  230    CONTINUE

C        Set up IGES Index Pointers

         IN1 = INDEX1 + 1 - IDEG1
         IN2 = INDEX2 + 1 - IDEG2
         IA  = IN1 + 2*IDEG1
         IB  = IN2 + 2*IDEG2
         IC  = (INDEX1+1)*(INDEX2+1)

C        Load C-Array
C
C        --  knots  --

         NUM    = IA+1
         JPAR   = 10
         INCPAR = 1
         INCDA  = 1
         CALL DTERPT(0)
         CALL D2IPFD (CMEM, IMEM, DMEM, NUM, IGE, JPAR, INCPAR, INCDA, 
     +      C(IP(1)+1), IERX)
         CALL DTERPT(1)
         IF (IERX .NE. 0) THEN
            IER = -8
            GOTO 9900
         ENDIF

         NUM    = IB+1
         JPAR   = 11+IA
         INCPAR = 1
         INCDA  = 1
         CALL DTERPT(0)
         CALL D2IPFD (CMEM, IMEM, DMEM, NUM, IGE, JPAR, INCPAR, INCDA, 
     +      C(IP(2)+1), IERX)
         CALL DTERPT(1)
         IF (IERX .NE. 0) THEN
            IER = -8
            GOTO 9900
         ENDIF

C        --  control points --
C
         IF (IPROP3 .EQ. 0) THEN
            DO 240 I = 1, 3
               NUM    = (INDEX1+1)*(INDEX2+1)
               JPAR   = 11+IA+IB+IC+I
               INCPAR = 3
               INCDA  = 1
               CALL DTERPT(0)
               CALL D2IPFD (CMEM, IMEM, DMEM, NUM, IGE, JPAR, INCPAR, 
     +            INCDA, C(IQ(I)+1), IERX)
               CALL DTERPT(1)
               IF (IERX .NE. 0) THEN
                  IER = -8
                  GOTO 9900
               ENDIF
         
               DO 240 J = 0, (INDEX1+1)*(INDEX2+1)-1
                  CALL DTERPT(0)
                  CALL D2IGFD (CMEM, IMEM, DMEM, IGE, 'P', 12+IA+IB+J, 
     +                  DVAL, IERX)
                  CALL DTERPT(1) 
                  IF (IERX .NE. 0) THEN
                     IER = -8
                     GOTO 9900
                  ENDIF
                  C(IQ(I)+1+J) = C(IQ(I)+1+J)*DVAL
  240       CONTINUE
  
            NUM    = (INDEX1+1)*(INDEX2+1)
            JPAR   = 12+IA+IB
            INCPAR = 1
            INCDA  = 1
            CALL DTERPT(0)
            CALL D2IPFD (CMEM, IMEM, DMEM, NUM, IGE, JPAR, INCPAR, 
     +         INCDA, C(IQ(4)+1), IERX)
            CALL DTERPT(1)
            IF (IERX .NE. 0) THEN
               IER = -8
               GOTO 9900
            ENDIF
         ELSE
            DO 250 I = 1, 3
               NUM    = (INDEX1+1)*(INDEX2+1)
               JPAR   = 11+IA+IB+IC+I
               INCPAR = 3
               INCDA  = 1
               CALL DTERPT(0)
               CALL D2IPFD (CMEM, IMEM, DMEM, NUM, IGE, JPAR, INCPAR, 
     +            INCDA, C(IQ(I)+1), IERX)
               CALL DTERPT(1)
               IF (IERX .NE. 0) THEN
                  IER = -8
                  GOTO 9900
               ENDIF
  250       CONTINUE
         ENDIF
      ELSE
         IER = -7
      ENDIF 
      
C     Error reporting section

 9900 CONTINUE

      IF (IER .NE. 0) THEN
          CALL DTERR (1, SUBNAM, IER, 0) 
          C(1) = 0
      ENDIF

      RETURN
      END
