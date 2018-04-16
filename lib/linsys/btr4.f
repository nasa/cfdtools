C+----------------------------------------------------------------------
C
      SUBROUTINE BTR4 ( MXBLOK, IL, IU, A, B, C, F )
C
C     BTR4 solves one block tridiagonal system with 4x4 blocks.
C     Block inversions use nonpivoted LU decomposition.
C
C     Arguments:
C     MXBLOK     Leading dimension of block arrays in calling program
C     A,B,C      4*4 sub-, main-, and super-diagonal blocks, numbered
C                by rows; B is overloaded
C     F          Forcing function; solution is output in F
C     IL, IU     Starting and finishing block indices
C
C     Author:    Joe Steger et al., NASA Ames, 1970s
C                COMMON block driven except for IL, IU
C     Mods:      David Saunders, Sterling Software, 6/9/86
C                Tidied up somewhat; argument driven except for /LUD44 /
C
C-----------------------------------------------------------------------

      REAL          A(MXBLOK,4,4), B(MXBLOK,4,4), C(MXBLOK,4,4),
     1              F(MXBLOK,4)

      COMMON/LUD44 /L11,L21,L22,L31,L32,L33,L41,L42,L43,L44,V1,V2,V3,V4,
     1              U12,U13,U14,U23,U24,U34
      REAL          L11,L21,L22,L31,L32,L33,L41,L42,L43,L44

      REAL          G(4,4)

      IS = IL +1
      I = IL
      DO 11 N=1,4
         G(N,1)=B(I,N,1)
         G(N,2)=B(I,N,2)
         G(N,3)=B(I,N,3)
         G(N,4)=B(I,N,4)
   11 CONTINUE

      CALL LUDEC(G)

      D1 = V1*F(I,1)
      D2 = V2*( F(I,2) - L21*D1)
      D3 = V3*( F(I,3) - L31*D1 - L32*D2)
      D4 = V4*( F(I,4) - L41*D1 - L42*D2 - L43*D3)
      F(I,4) = D4
      F(I,3) = D3 - U34*D4
      F(I,2) = D2 - U24*D4 - U23*F(I,3)
      F(I,1) = D1 - U14*D4 - U13*F(I,3) - U12*F(I,2)

      DO 12 M=1,4
         D1 = V1*C(I,1,M)
         D2 = V2*( C(I,2,M) - L21*D1)
         D3 = V3*( C(I,3,M) - L31*D1 - L32*D2)
         D4 = V4*( C(I,4,M) -L41*D1 - L42*D2 - L43*D3)
         B(I,4,M) = D4
         B(I,3,M) = D3 - U34*D4
         B(I,2,M) = D2 - U24*D4 - U23*B(I,3,M)
         B(I,1,M) = D1 - U14*D4 - U13*B(I,3,M) - U12*B(I,2,M)
   12 CONTINUE

      DO 13 I=IS,IU
         IR = I -1
         DO 14 N=1,4
            F(I,N) = F(I,N) - A(I,N,1)*F(IR,1) - A(I,N,2)*F(IR,2)
     1                      - A(I,N,3)*F(IR,3) - A(I,N,4)*F(IR,4)
C****       DO 14 M=1,4
C****          G(N,M) = B(I,N,M)
C****1                  - A(I,N,1)*B(IR,1,M) - A(I,N,2)*B(IR,2,M)
C****2                  - A(I,N,3)*B(IR,3,M) - A(I,N,4)*B(IR,4,M)

            G(N,1) = B(I,N,1)
     1               - A(I,N,1)*B(IR,1,1) - A(I,N,2)*B(IR,2,1)
     2               - A(I,N,3)*B(IR,3,1) - A(I,N,4)*B(IR,4,1)
            G(N,2) = B(I,N,2)
     1               - A(I,N,1)*B(IR,1,2) - A(I,N,2)*B(IR,2,2)
     2               - A(I,N,3)*B(IR,3,2) - A(I,N,4)*B(IR,4,2)
            G(N,3) = B(I,N,3)
     1               - A(I,N,1)*B(IR,1,3) - A(I,N,2)*B(IR,2,3)
     2               - A(I,N,3)*B(IR,3,3) - A(I,N,4)*B(IR,4,3)
            G(N,4) = B(I,N,4)
     1               - A(I,N,1)*B(IR,1,4) - A(I,N,2)*B(IR,2,4)
     2               - A(I,N,3)*B(IR,3,4) - A(I,N,4)*B(IR,4,4)
   14    CONTINUE

         CALL LUDEC(G)

         D1 = V1*F(I,1)
         D2 = V2*( F(I,2) - L21*D1)
         D3 = V3*( F(I,3) - L31*D1 - L32*D2)
         D4 = V4*( F(I,4) - L41*D1 - L42*D2 - L43*D3)
         F(I,4) = D4
         F(I,3) = D3 - U34*D4
         F(I,2) = D2 - U24*D4 - U23*F(I,3)
         F(I,1) = D1 - U14*D4 - U13*F(I,3) - U12*F(I,2)

         IF (I .EQ. IU) GO TO 13

         DO 15 M=1,4
            D1 = V1*C(I,1,M)
            D2 = V2*( C(I,2,M) - L21*D1)
            D3 = V3*( C(I,3,M) - L31*D1 - L32*D2)
            D4 = V4*( C(I,4,M) - L41*D1 - L42*D2 - L43*D3)
            B(I,4,M) = D4
            B(I,3,M) = D3 - U34*D4
            B(I,2,M) = D2 - U24*D4 - U23*B(I,3,M)
            B(I,1,M) = D1 - U14*D4 - U13*B(I,3,M) - U12*B(I,2,M)
   15    CONTINUE
   13 CONTINUE

      IT = IL + IU
      DO 21 II = IS,IU
         I = IT - II
         IP = I+1
         DO 22 N=1,4
            F(I,N) = F(I,N) - B(I,N,1)*F(IP,1) - B(I,N,2)*F(IP,2)
     1                      - B(I,N,3)*F(IP,3) - B(I,N,4)*F(IP,4)
   22    CONTINUE
   21 CONTINUE
      RETURN
      END
C+----------------------------------------------------------------------
      SUBROUTINE LUDEC(A)
C     LU decomposition for 4x4 matrix (no pivoting).   Steger, et al.
C-----------------------------------------------------------------------
      REAL          A(4,4)
      COMMON/LUD44 /L11,L21,L22,L31,L32,L33,L41,L42,L43,L44,V1,V2,V3,V4,
     1              U12,U13,U14,U23,U24,U34
      REAL          L11,L21,L22,L31,L32,L33,L41,L42,L43,L44
      L11 = A(1,1)
      V1  = 1./L11
      U12 = V1*A(1,2)
      U13 = V1*A(1,3)
      U14 = V1*A(1,4)
      L21 = A(2,1)
      L22 = A(2,2) - L21*U12
      V2  = 1./L22
      U23 = ( A(2,3) - L21*U13)* V2
      U24 = ( A(2,4) - L21*U14)* V2
      L31 = A(3,1)
      L32 = A(3,2) - L31*U12
      L33 = A(3,3) - L31*U13 - L32*U23
      V3  = 1./L33
      U34 = ( A(3,4) - L31*U14 - L32*U24)* V3
      L41 = A(4,1)
      L42 = A(4,2) - L41*U12
      L43 = A(4,3) - L41*U13 - L42*U23
      L44 = A(4,4) - L41*U14 - L42*U24 - L43*U34
      V4  = 1./L44
      RETURN
      END
C+---------------------------------------------------------------------
C
      SUBROUTINE BTR4P ( MXBLOK, IL, IU, A, B, C, F, H, U )
C
C     BTR4P solves one periodic block tridiagonal system (4x4 blocks).
C     This version does not overload the block matrices.
C     Block inversions use nonpivoted LU decomposition.
C
C     Arguments:
C     MXBLOK     Leading dimension of block arrays in calling program
C     A,B,C      4*4 sub-, main-, and super-diagonal blocks, numbered
C                by rows
C     F          Forcing function; solution is output in F
C     H,U        Work space
C     IL, IU     Starting and finishing block indices
C
C     Author:    Joe Steger et al., NASA Ames, 1970s
C                COMMON block driven except for IL, IU
C     Mods:      David Saunders, Sterling Software, 6/9/86
C                Tidied up somewhat; argument driven except for /LUD44 /
C
C-----------------------------------------------------------------------

      REAL          A(MXBLOK,4,4), B(MXBLOK,4,4), C(MXBLOK,4,4),
     1              H(MXBLOK,4,4), U(MXBLOK,4,4), F(MXBLOK,4)

      COMMON/LUD44 /L11,L21,L22,L31,L32,L33,L41,L42,L43,L44,V1,V2,V3,V4,
     1              U12,U13,U14,U23,U24,U34
      REAL          L11,L21,L22,L31,L32,L33,L41,L42,L43,L44

      REAL          G(4,4), E(4,4), Q(4,4), FN(4), S(4,4)
      PARAMETER     ( MN=4 )

      IS = IL +1
      IM = IU -1
      IN = IU -2

C  FORWARD SWEEP
      I = IL
      DO 10 N=1,MN
         DO 10 M=1,MN
            E(N,M) = C(IU,N,M)
   10       G(N,M) = B(I,N,M)

      CALL LUDEC(G)

      D1 = V1*F(I,1)
      D2 = V2*( F(I,2) - L21*D1)
      D3 = V3*( F(I,3) - L31*D1 - L32*D2)
      D4 = V4*( F(I,4) - L41*D1 - L42*D2 - L43*D3)
      F(I,4) = D4
      F(I,3) = D3 - U34*D4
      F(I,2) = D2 - U24*D4 - U23*F(I,3)
      F(I,1) = D1 - U14*D4 - U13*F(I,3) - U12*F(I,2)

      DO 11 M=1,MN
         D1 = V1*C(I,1,M)
         D2 = V2*( C(I,2,M) - L21*D1)
         D3 = V3*( C(I,3,M) - L31*D1 - L32*D2)
         D4 = V4*( C(I,4,M) - L41*D1 - L42*D2 - L43*D3)
         U(I,4,M) = D4
         U(I,3,M) = D3 - U34*D4
         U(I,2,M) = D2 - U24*D4 - U23*U(I,3,M)
         U(I,1,M) = D1 - U14*D4 - U13*U(I,3,M) - U12*U(I,2,M)

         D1 = V1*A(I,1,M)
         D2 = V2*( A(I,2,M) - L21*D1)
         D3 = V3*( A(I,3,M) - L31*D1 - L32*D2)
         D4 = V4*( A(I,4,M) -L41*D1 - L42*D2 - L43*D3)
         H(I,4,M) = D4
         H(I,3,M) = D3 - U34*D4
         H(I,2,M) = D2 - U24*D4 - U23*H(I,3,M)
         H(I,1,M) = D1 - U14*D4 - U13*H(I,3,M) - U12*H(I,2,M)
   11 CONTINUE

      DO 13 N=1,MN
         FN(N) = F(IU,N) -E(N,1)*F(I,1) -E(N,2)*F(I,2) -E(N,3)*F(I,3)
     1                   -E(N,4)*F(I,4)
         DO 13 M=1,MN
            Q(N,M) = B(IU,N,M) -E(N,1)*H(I,1,M) -E(N,2)*H(I,2,M)
     1                         -E(N,3)*H(I,3,M) -E(N,4)*H(I,4,M)
   13 CONTINUE

      DO 14 I=IS,IN
         L = I-1
         DO 15 N=1,MN
            F(I,N) = F(I,N) - A(I,N,1)*F(L,1) -A(I,N,2)*F(L,2)
     1                      - A(I,N,3)*F(L,3) -A(I,N,4)*F(L,4)
            DO 15 M=1,MN
               G(N,M) = B(I,N,M) -A(I,N,1)*U(L,1,M) -A(I,N,2)*U(L,2,M)
     1                           -A(I,N,3)*U(L,3,M) -A(I,N,4)*U(L,4,M)
               S(N,M) = -E(N,1)*U(L,1,M) -E(N,2)*U(L,2,M)
     1                  -E(N,3)*U(L,3,M) -E(N,4)*U(L,4,M)
               H(I,N,M) = -A(I,N,1)*H(L,1,M) -A(I,N,2)*H(L,2,M)
     1                    -A(I,N,3)*H(L,3,M) -A(I,N,4)*H(L,4,M)
   15    CONTINUE

         DO 7 N=1,MN
            DO 7 M=1,MN
    7          E(N,M) = S(N,M)

         CALL LUDEC(G)

         D1 = V1*F(I,1)
         D2 = V2*( F(I,2) - L21*D1)
         D3 = V3*( F(I,3) - L31*D1 - L32*D2)
         D4 = V4*( F(I,4) - L41*D1 - L42*D2 - L43*D3)
         F(I,4) = D4
         F(I,3) = D3 - U34*D4
         F(I,2) = D2 - U24*D4 - U23*F(I,3)
         F(I,1) = D1 - U14*D4 - U13*F(I,3) - U12*F(I,2)

         DO 16 M=1,MN
            D1 = V1*C(I,1,M)
            D2 = V2*( C(I,2,M) - L21*D1)
            D3 = V3*( C(I,3,M) - L31*D1 - L32*D2)
            D4 = V4*( C(I,4,M) -L41*D1 - L42*D2 - L43*D3)
            U(I,4,M) = D4
            U(I,3,M) = D3 - U34*D4
            U(I,2,M) = D2 - U24*D4 - U23*U(I,3,M)
            U(I,1,M) = D1 - U14*D4 - U13*U(I,3,M) - U12*U(I,2,M)

            D1 = V1*H(I,1,M)
            D2 = V2*( H(I,2,M) - L21*D1)
            D3 = V3*( H(I,3,M) - L31*D1 - L32*D2)
            D4 = V4*( H(I,4,M) -L41*D1 - L42*D2 - L43*D3)
            H(I,4,M) = D4
            H(I,3,M) = D3 - U34*D4
            H(I,2,M) = D2 - U24*D4 - U23*H(I,3,M)
            H(I,1,M) = D1 - U14*D4 - U13*H(I,3,M) - U12*H(I,2,M)
   16    CONTINUE

         DO 12 N=1,MN
            FN(N) = FN(N) -E(N,1)*F(I,1) -E(N,2)*F(I,2) -E(N,3)*F(I,3)
     1                    -E(N,4)*F(I,4)
            DO 12 M=1,MN
               Q(N,M) = Q(N,M) -E(N,1)*H(I,1,M) -E(N,2)*H(I,2,M)
     1                         -E(N,3)*H(I,3,M) -E(N,4)*H(I,4,M)
   12    CONTINUE
   14 CONTINUE
C
      I = IM
      L = I-1
      DO 17 N=1,MN
         F(I,N) = F(I,N) - A(I,N,1)*F(L,1) - A(I,N,2)*F(L,2)
     1                   - A(I,N,3)*F(L,3) - A(I,N,4)*F(L,4)
         DO 17 M=1,MN
            G(N,M) = B(I,N,M) -A(I,N,1)*U(L,1,M) -A(I,N,2)*U(L,2,M)
     1                        -A(I,N,3)*U(L,3,M) -A(I,N,4)*U(L,4,M)
            S(N,M) = -E(N,1)*U(L,1,M) -E(N,2)*U(L,2,M) -E(N,3)*U(L,3,M)
     1               -E(N,4)*U(L,4,M)  + A(IU,N,M)
            H(I,N,M) = -A(I,N,1)*H(L,1,M) -A(I,N,2)*H(L,2,M)
     1                 -A(I,N,3)*H(L,3,M) -A(I,N,4)*H(L,4,M) + C(I,N,M)
   17 CONTINUE

      DO 8 N=1,MN
         DO 8 M=1,MN
    8       E(N,M) = S(N,M)

      CALL LUDEC(G)

      D1 = V1*F(I,1)
      D2 = V2*( F(I,2) - L21*D1)
      D3 = V3*( F(I,3) - L31*D1 - L32*D2)
      D4 = V4*( F(I,4) - L41*D1 - L42*D2 - L43*D3)
      F(I,4) = D4
      F(I,3) = D3 - U34*D4
      F(I,2) = D2 - U24*D4 - U23*F(I,3)
      F(I,1) = D1 - U14*D4 - U13*F(I,3) - U12*F(I,2)

      DO 18 M=1,MN
         D1 = V1*H(I,1,M)
         D2 = V2*( H(I,2,M) - L21*D1)
         D3 = V3*( H(I,3,M) - L31*D1 - L32*D2)
         D4 = V4*( H(I,4,M) -L41*D1 - L42*D2 - L43*D3)
         H(I,4,M) = D4
         H(I,3,M) = D3 - U34*D4
         H(I,2,M) = D2 - U24*D4 - U23*H(I,3,M)
         H(I,1,M) = D1 - U14*D4 - U13*H(I,3,M) - U12*H(I,2,M)
   18 CONTINUE

      DO 19 N=1,MN
         FN(N) = FN(N) -E(N,1)*F(I,1) -E(N,2)*F(I,2) -E(N,3)*F(I,3)
     1                 -E(N,4)*F(I,4)
         DO 19 M=1,MN
            Q(N,M) = Q(N,M) -E(N,1)*H(I,1,M) -E(N,2)*H(I,2,M)
     1                      -E(N,3)*H(I,3,M) -E(N,4)*H(I,4,M)
   19 CONTINUE

      I = IU
      CALL LUDEC(Q)

      D1 = V1*FN( 1)
      D2 = V2*( FN( 2) - L21*D1)
      D3 = V3*( FN( 3) - L31*D1 - L32*D2)
      D4 = V4*( FN( 4) - L41*D1 - L42*D2 - L43*D3)
      F(I,4) = D4
      F(I,3) = D3 - U34*D4
      F(I,2) = D2 - U24*D4 - U23*F(I,3)
      F(I,1) = D1 - U14*D4 - U13*F(I,3) - U12*F(I,2)

C  BACKWARD SWEEP
      K = IU
      I = IM
      DO 20 N=1,MN
         F(I,N) = F(I,N) -H(I,N,1)*F(K,1) -H(I,N,2)*F(K,2)
     1                   -H(I,N,3)*F(K,3) -H(I,N,4)*F(K,4)
   20 CONTINUE

      IT = IL + IN
      DO 22 J=IL,IN
         I = IT - J
         L = I+1
         DO 21 N=1,MN
            F(I,N) = F(I,N) -U(I,N,1)*F(L,1) -U(I,N,2)*F(L,2)
     1               -U(I,N,3)*F(L,3) -U(I,N,4)*F(L,4) -H(I,N,1)*F(K,1)
     2               -H(I,N,2)*F(K,2) -H(I,N,3)*F(K,3) -H(I,N,4)*F(K,4)
   21    CONTINUE
   22 CONTINUE
      RETURN
      END
