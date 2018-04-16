C+----------------------------------------------------------------------
C
      SUBROUTINE DTDLSQ ( N, A, B, C, D, R, S, SSQ )
C
C ACRONYM: Diagonal + TriDiagonal system; Least SQuares solution
C          -          -  -                -     --
C
C PURPOSE: DTDLSQ solves one overdetermined linear system of the form
C
C                     |       |  |   |       |   |
C                     |   D   |  | x |       | r |
C                     |       |  |   |       |   |
C                     | - - - |          =   | - |
C                     |       |              |   |
C                     |   T   |              | s |
C                     |       |              |   |
C
C          where  D  is diagonal, and  T  is tridiagonal.
C
C METHOD:  The QR decomposition of the left-hand-side matrix is used,
C          where  R  is readily shown to be upper tridiagonal.  R  is
C          produced by a sequence of symmetric plane rotations, which
C          are also applied to the right-hand-side vector.  The least
C          squares solution follows in the usual way by  solving  the
C          (square upper tridiagonal) portion of the transformed sys-
C          tem.
C
C          Annihilation of the tridiagonal lower portion is  done  in
C          3-part steps involving a(i), c(i+1), b(i) in that order as
C          indicated for i=1:
C
C                           | d         |     | r r r     |
C                           |   d       |     |   x       |
C                           |     d     |     |     d     |
C                           |       d   |     |       d   |
C                           |         d |     |         d |
C                  Q3 Q2 Q1 | - - - - - |  =  | - - - - - |
C                           | a b       |     | 0 0       |
C                           | c a b     |     | 0 x x     |
C                           |   c a b   |     |   c a b   |
C                           |     c a b |     |     c a b |
C                           |       c a |     |       c a |
C
C          Note that all elements of a, b, c, d  are numbered by row.
C          Vectors a, b, c are overwritten with the  upper triangular
C          factor R.  The right-hand-side vectors are transformed in-
C          place, and the least squares solution overwrites the upper
C          half, r.   The minimum sum of squares is derived from  the
C          rest of the transformed right-hand-side vector s.
C
C ARGUMENTS:
C   ARG  DIM TYPE I/O/S DESCRIPTION
C    N    -    I    I   Order of the system (N >= 3).
C    A    N    R   I/S  Input with main diagonal of  tridiagonal part
C                       of the matrix; destroyed upon return.
C    B    N    R   I/S  Input with upper diagonal of tridiagonal part
C                       of the matrix in B(1):B(N-1);  B(N) is needed
C                       to avoid special handling of i=N-1 iteration.
C    C    N    R   I/S  Input with lower diagonal of tridiagonal part
C                       of the matrix in  C(2):C(N);  destroyed  upon
C                       return.
C    D    N    R   I/S  Input with upper diagonal of the matrix;  de-
C                       stroyed upon return.
C    R    N    R  I/S/O Input with upper half of right-hand side vec-
C                       tor; returned with the least squares solution
C                       of the overdetermined system.
C    S    N    R   I/S  Input with lower half of right-hand-side vec-
C                       tor; destroyed upon return.
C    SSQ  -    R    O   Minimum sum of squares corresponding  to  the
C                       calculated solution.
C
C ERROR HANDLING:  None.  Divide by zero can mean the matrix is sing-
C                  ular but is more likely to mean  unwitting  errors
C                  in the call to DTDLSQ.
C
C ENVIRONMENT:  VAX/VMS FORTRAN 77
C
C DEVELOPMENT HISTORY:
C     DATE  INITIALS  DESCRIPTION
C   11/01/84   MAS    Algorithm description.
C   Nov. '84   DAS    Initial implementation,  for efficient solution
C                     of the thickness/curvature  refinement  problem
C                     in airfoil design.
C
C AUTHORS: Michael Saunders, Stanford University, Palo Alto, CA.
C          David Saunders,   Informatics General, Palo Alto, CA.
C
C-----------------------------------------------------------------------

      IMPLICIT  NONE

C ... Arguments:

      INTEGER   N
      REAL      A(N), B(N), C(N), D(N), R(N), S(N), SSQ

C ... Local variables:

      INTEGER   I
      REAL      A2, B1, B2, CS, RI, R1, R11, R12, R13, R2, R22,
     +          SI, SN, S1, T1, T2

C ... Intrinsics:

      INTRINSIC SQRT

C ... Execution:

      DO 20 I = 1, N-1

C ...    Eliminate a(i):

         T1 = SQRT ( D(I) ** 2 + A(I) ** 2 )
         CS = D(I) / T1
         SN = A(I) / T1
         T2 = SN * B(I)
         B1 =-CS * B(I)

C ...    Modify affected elements of RHS:

         R1 = R(I)
         S1 = S(I)
         RI = CS * R1 + SN * S1
         SI = SN * R1 - CS * S1

C ...    Eliminate c(i+1):

         R11 = SQRT ( T1 ** 2 + C(I+1) ** 2 )
         CS  = T1 / R11
         SN  = C(I+1) / R11
         R12 = CS * T2 + SN * A(I+1)
         A2  = SN * T2 - CS * A(I+1)
         R13 =           SN * B(I+1)
         B2  =          -CS * B(I+1)
         R(I)   = CS * RI + SN * S(I+1)
         S(I+1) = SN * RI - CS * S(I+1)

C ...    Eliminate b(i):

         R22 = SQRT ( D(I+1) ** 2 + B1 ** 2 )
         CS  = D(I+1) / R22
         SN  = B1 / R22
         R2  = R(I+1)
         R(I+1) = CS * R2 + SN * SI
         S(I)   = SN * R2 - CS * SI

C ...    Here is the ith row of upper triangular factor R ...

         A(I) = R11
         B(I) = R12
         C(I) = R13

C ...    ... and the modified (i+1)th row of T and D:

         A(I+1) = A2
         B(I+1) = B2
         D(I+1) = R22

   20 CONTINUE

C ... Nth iteration eliminates just one element of T:

      T1 = SQRT ( D(N) ** 2 + A(N) ** 2 )
      CS = D(N) / T1
      SN = A(N) / T1
      A(N) = T1
      RI   = R(N)
      R(N) = CS * RI + SN * S(N)
      S(N) = SN * RI - CS * S(N)

C ... Solve the upper portion of the triangularized system:

      R(N) = R(N) / A(N)
      R(N-1) = ( R(N-1) - B(N-1) * R(N) ) / A(N-1)

      DO 30 I = N-2, 1, -1
         R(I) = ( R(I) - B(I) * R(I+1) - C(I) * R(I+2) ) / A(I)
   30 CONTINUE

C ... Minimum residual is contained in lower part of modified RHS:

      SSQ = S(1)
      DO 40 I = 2, N
         SSQ = S(I) ** 2 + SSQ
   40 CONTINUE

      RETURN
      END
