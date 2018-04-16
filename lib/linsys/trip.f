C+----------------------------------------------------------------------
C
      SUBROUTINE TRIP ( A, B, C, F, Q, S, J1, J2 )
C
C  PURPOSE:  TRIP solves one scalar periodic tridiagonal system.
C            This module was extracted from program ARC2D (CFD Branch,
C            NASA Ames) - author unknown.
C
C  PARAMS:   Presumably the usual restrictions regarding diagonal
C            dominance apply. There was no description for the argument
C            list; the following is an educated guess:
C
C            J1, J2 are the first and last indices defining the system.
C                   All of the following are vectors of length J2.
C            A(*)   is the subdiagonal;
C            C(*)   is the superdiagonal;  A(J1) and C(J2) are the
C                   elements which make the system periodic.
C            B(*)   is the main diagonal.
C            F(*)   is input with the right-hand-side vector;
C              "    " output with the solution.
C            Q(*)   is work-space.
C            S(*)   is work-space.
C
C-----------------------------------------------------------------------

      DIMENSION  A(J2), B(J2), C(J2), F(J2), Q(J2), S(J2)

      JA = J1 + 1
      FN = F(J2)
C
C         FORWARD ELIMINATION SWEEP
C
      Q(J1) = -C(J1)/B(J1)
      F(J1) = F(J1)/B(J1)
      S(J1) = - A(J1)/B(J1)
      DO 10 J=JA,J2
         P =1./( B(J) + A(J)*Q(J-1))
         Q(J) = - C(J)*P
         F(J) = ( F(J) - A(J)*F(J-1))*P
         S(J) = - A(J)*S(J-1)*P
   10 CONTINUE
C
C         BACKWARD PASS
C
      JJ = J1 + J2
      Q(J2) = 0.
      S(J2) = 1.
      DO 11 I=JA,J2
         J = JJ - I
         S(J) = S(J) + Q(J)*S(J+1)
         Q(J) = F(J) + Q(J)*Q(J+1)
11    CONTINUE

      F(J2) = ( FN - C(J2)*Q(J1) - A(J2)*Q(J2-1))/(C(J2)*S(J1) +
     1        A(J2)*S(J2-1)  +B(J2))
C
C         BACKWARD ELIMINATION PASS
C
      DO 12 I=JA,J2
         J = JJ -I
         F(J) = F(J2)*S(J) + Q(J)
12    CONTINUE

      RETURN
      END
