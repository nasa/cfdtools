      SUBROUTINE DTPRNT (C)
C
C     TRIVIAL SUBROUTINE TO PRINT OUT THE CONTENTS OF A C ARRAY
C
      EXTERNAL DTSCHK, DTMCON, DTDCLB
C
      INTEGER I, J, N, M, IP1, IP2, IER, DTDCLB
      DOUBLE PRECISION C(*), R, DTMCON
      CHARACTER*12 CFMT
C
      IF (DTDCLB(C(1)) .EQ. 1) THEN
          WRITE (6,*) ' << (DTPRNT): Invalid C array >>'
          RETURN
      ENDIF
C
      CALL DTSCHK (C,IER)
      IF (IER .NE. 0) THEN
          WRITE (6,*) ' << (DTPRNT): Invalid C array >>'
          RETURN
      ENDIF
      N = INT(C(1))
      M = ABS(INT(C(2)))
C
      WRITE (6,*) ' Independent variables (n):'
      WRITE (6,'(5X,F8.3)') C(1)
      WRITE (6,*) ' '
C
      WRITE (6,*) ' Dependent variables (m):'
      WRITE (6,'(5X,F8.3)') C(2)
      WRITE (6,*) ' '
C
      WRITE (6,*) ' Order of the splines (k) for each indep. var.:'
      WRITE (6,'(5X,9F8.3)') (C(I),I=3,N+2)
      WRITE (6,*) ' '
C
      WRITE (6,*) ' Number of B-spline coeff. for each indep. var.:'
      WRITE (6,'(5X,9F8.3)') (C(I),I=N+3,2*N+2)
      WRITE (6,*) ' '
C
      WRITE (6,*) ' Index of last span (jspan) for each indep. var.:'
      WRITE (6,'(5X,9F8.3)') (C(I),I=2*N+3,3*N+2)
      WRITE (6,*) ' '
C
      WRITE (6,*) ' Knots:'
      IP1 = 3*N+2
      IP2 = IP1+C(3)+C(N+3)
      DO 10 J = 1, N
          WRITE (6,'(5X,9F8.3)') (C(I),I=IP1+1,IP2)
          WRITE (6,*) ' '
          IP1 = IP2
          IP2 = IP2+C(3+J)+C(N+3+J)
   10 CONTINUE
C
      WRITE (6,*) ' Coefficients:'
      R = 1.0
      DO 15 I = N+3, 2*N+2
          R = R*C(I)
   15 CONTINUE
      IP2 = IP1+R
      DO 20 J = 1, M
          WRITE (CFMT,'(''(5X,'',I3,''F8.3)'')') INT(C(N+3))
          WRITE (6,CFMT) (C(I),I=IP1+1,IP2)
          WRITE (6,*) ' '
          IP1 = IP2
          IP2 = IP2+R
   20 CONTINUE
C
      RETURN
      END

