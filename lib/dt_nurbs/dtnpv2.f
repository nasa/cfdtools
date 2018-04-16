C<*>
      SUBROUTINE DTNPV2 ( X,      INCX,   NDOM,   NRNG,   C,
     1                    NCOEF,  COEF,   NDIM,   NZERO,  INDX,
     2                    IPOS,   AT,     BSVAL,  WORK,   V,
     3                    IER       )
C
C     ==================================================================
C
C     --------------
C     ... PARAMETERS
C     --------------
C
      INTEGER      INCX,     NDOM,     NRNG,     NCOEF,    NDIM,
     1             NZERO,    IPOS(*),  INDX(*),  IER
C
      DOUBLE PRECISION       X(*),     C(*),     COEF(NCOEF,*),
     1             AT(NDIM,*),  BSVAL(*), WORK(*),  V(*)
C
C     ----------------------
C     ... INTERNAL VARIABLES
C     ----------------------
C
      INTEGER      I,        IK,       IKPT,     IL,       ILAST,
     2             ILPT,     INCZ,     INPT,     IOPT,     ISHFT,
     3             ISTRT,    IX,       J,        KORD,     N,
     4             NBS,      NSHFT
C
C     -------------
C     ... EXTERNALS
C     -------------
C
      DOUBLE PRECISION    DDOT
      EXTERNAL            DTBSP2,   DTNPV3,   DCOPY,    DDOT
C
C     ==================================================================
C
C     ---------------------------
C     ...DEFINE LOCATION POINTERS
C     ---------------------------
C
      IOPT = 2
      INPT = IOPT + NDOM
      ILPT = INPT + NDOM
      IKPT = 3 + 3 * NDOM
      IX   = 1 - INCX
      INCZ = 1
C
C     -------------------------------------------
C     ... DEFINE INDEX VECTOR AND SET IPOS VECTOR
C     -------------------------------------------
C
      INDX(1) = 1
      NSHFT   = 1
C
      DO 10 N = 1, NDOM
C
          KORD    = INT ( C(IOPT+N) )
          NBS     = INT ( C(INPT+N) )
          IPOS(N) = INT ( C(ILPT+N) )
          IX      = IX + INCX
C
C     ----------------
C     ... LOCATE X(IX)
C     ----------------
C
          CALL DTNPV3 ( X(IX), C(IKPT), NBS, KORD, IPOS(N), IER )
C
          C(ILPT+N) = DBLE ( IPOS(N) )
C
          IF ( IER .NE. 0 ) THEN
C
              RETURN
C
          END IF
C
C      ------------------------------
C      ... CHECK FOR INVALID KNOT SET
C      ------------------------------
C
          CALL DTILC1( C(IKPT + IPOS(N) - KORD), 2 * KORD, KORD, IFAIL)
C
          IF ( IFAIL .NE. 0) THEN
C
              IER = -8
C
              RETURN
C
          END IF
C
          INDX(1) = INDX(1) + ( IPOS(N) - KORD ) * NSHFT
          NSHFT   = NSHFT * NBS
          IKPT    = IKPT + NBS + KORD
C
 10   CONTINUE
C
C     --------------------------------
C     ... BUILD REMAINING INDEX VECTOR
C     --------------------------------
C
      NSHFT = INT ( C(INPT+1) )
      ILAST = 1
      ISHFT = 1
C
      DO 30 N = 2, NDOM
C
          ISTRT = ILAST + 1
          ILAST = ILAST * INT ( C(IOPT+N) )
C
          DO 20 I = ISTRT, ILAST
C
              INDX(I) = INDX(I-ISHFT) + NSHFT
C
 20       CONTINUE
C
          ISHFT = ILAST
          NSHFT = NSHFT * INT ( C(INPT+N) )
C
 30   CONTINUE
C
C     ==================================================================
C
C     ---------------------------
C     ... EVALUATE IN X DIRECTION
C     ---------------------------
C
C     -------------------------
C     ... SET SPLINE PARAMETERS
C     -------------------------
C
      KORD = INT( C(IOPT+1) )
      NBS  = INT( C(INPT+1) )
      IKPT = 3 + 3 * NDOM
      IX   = 1
C
C     ----------------------
C     ... EVALUATE B-SPLINES
C     ----------------------
C
      IK = -1
      CALL DTBSP2 ( C(IKPT),    X(1), IPOS(1), IK, KORD, WORK(1),
     1              WORK(KORD), BSVAL                         )
C
C     ------------------------
C     ... COLLAPSE X DIRECTION
C     ------------------------
C
      DO 50 I = 1, NZERO
C
          DO 40 J = 1, NRNG

              AT(I,J) = DDOT (KORD, COEF(INDX(I),J), INCZ, BSVAL, INCZ)
C
 40       CONTINUE
C
 50   CONTINUE
C
      IKPT      = IKPT + KORD + NBS
C
C     ==================================================================
C
C     ----------------------------------
C     ... LOOP THROUGH HIGHER DIMENSIONS
C     ----------------------------------
C
      DO 80 N = 2, NDOM
C
C     ------------------
C     ... SET PARAMETERS
C     ------------------
C
          KORD  = INT ( C(IOPT+N) )
          NBS   = INT ( C(INPT+N) )
          NZERO = NZERO / KORD
          IX    = IX + INCX
C
C     ----------------------
C     ... EVALUATE B-SPLINES
C     ----------------------
C
          IK = -1
          CALL DTBSP2( C(IKPT),    X(IX), IPOS(N), IK, KORD, WORK(1),
     1                 WORK(KORD), BSVAL                             )
C
C     ---------------------------
C     ... COLLAPSE N-TH DIRECTION
C     ---------------------------
C
          IL = 1 - KORD
C
          DO 70 I = 1, NZERO
C
              IL = IL + KORD
C
              DO 60 J = 1, NRNG
C
                  AT(I,J) = DDOT ( KORD, AT(IL,J), INCZ, BSVAL, INCZ )
C
 60           CONTINUE
C
 70       CONTINUE
C
          IKPT      = IKPT + KORD + NBS
C
 80   CONTINUE
C
C     ----------------
C     ... COPY RESULTS
C     ----------------
C
      CALL DCOPY ( NRNG, AT(1,1), NDIM, V, INCZ )
C
C     ==================================================================
C
C     ----------
C     ... RETURN
C     ----------
C
      RETURN
      END
